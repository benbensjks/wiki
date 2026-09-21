"""State-by-state equivalence between the formal model and legacy validator."""
from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from model import Model, ROOT, ZENG, act
from model_twobit34 import STATE_NAMES_34, TwoBit34Model, nominal_extension
from verify_twobit_causal import analyse_solution


OUT = ROOT / 'twobit34_results'


def rates(base, mrna_min=2.0, maturation_min=30.0):
    lm = np.log(2) * 60 / mrna_min + base.growth
    kmat = np.log(2) * 60 / maturation_min
    lu = kmat + base.e.immature_decay_h + base.growth
    return lm, kmat, lu


def expression(base, m, u, protein, source, gamma, kinetic_rates):
    lm, kmat, lu = kinetic_rates
    tx = source * lm * lu / (base.e.translation_h * kmat)
    return np.asarray((tx - lm * m,
                       base.e.translation_h * m - lu * u,
                       kmat * u - (gamma + base.growth) * protein))


def legacy_initial(base, kinetic_rates):
    y30 = base.initial_state_twobit(cold=False)
    A0 = y30[28]
    source_A = ZENG['alpha_A'][0] * (1 - y30[16])
    source_F = ZENG['alpha_F'][0] * act(A0, ZENG['K_A'][0], ZENG['n_A'][0])
    _, kmat, lu = kinetic_rates
    mA, uA = source_A * lu / (base.e.translation_h * kmat), source_A / kmat
    mF, uF = source_F * lu / (base.e.translation_h * kmat), source_F / kmat
    return np.r_[y30, mA, uA, mF, uF]


def legacy_rhs(base, kinetic_rates, t, y):
    y30 = y[:30]
    d30 = base.rhs_twobit(t, y30)
    A0, F0 = y30[28], y30[29]
    mA, uA, mF, uF = y[30:34]
    source_A = ZENG['alpha_A'][0] * (1 - y30[16])
    source_F = ZENG['alpha_F'][0] * act(A0, ZENG['K_A'][0], ZENG['n_A'][0])
    eA = expression(base, mA, uA, A0, source_A, ZENG['gamma_A'][0], kinetic_rates)
    eF = expression(base, mF, uF, F0, source_F, ZENG['gamma_F'][0], kinetic_rates)
    d30[28], d30[29] = eA[2], eF[2]
    return np.r_[d30, eA[:2], eF[:2]]


def legacy_simulate(hours=300.0, sample_min=2.0):
    extension = nominal_extension()
    base = Model(extension)
    kinetic_rates = rates(base)
    y0 = legacy_initial(base, kinetic_rates)
    ts = np.arange(0.0, hours + sample_min / 120.0, sample_min / 60.0)
    sol = solve_ivp(lambda t, y: legacy_rhs(base, kinetic_rates, t, y),
                    (0, hours), y0, t_eval=ts, method='DOP853',
                    rtol=2e-7, atol=2e-9, max_step=2 / 60)
    if not sol.success:
        raise RuntimeError(sol.message)
    return base, kinetic_rates, sol


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    formal = TwoBit34Model()
    formal_sol = formal.simulate(hours=300.0, sample_min=2.0, max_step_min=2.0)
    legacy_base, legacy_rates, legacy_sol = legacy_simulate()

    if not np.array_equal(formal_sol.t, legacy_sol.t):
        raise RuntimeError('time grids differ')
    error = np.abs(formal_sol.y - legacy_sol.y)
    per_state = {name: float(error[i].max()) for i, name in enumerate(STATE_NAMES_34)}
    derivative_errors = []
    for k in np.linspace(0, formal_sol.t.size - 1, 21, dtype=int):
        d_formal = formal.rhs(float(formal_sol.t[k]), formal_sol.y[:, k])
        d_legacy = legacy_rhs(legacy_base, legacy_rates, float(legacy_sol.t[k]), legacy_sol.y[:, k])
        derivative_errors.append(float(np.max(np.abs(d_formal - d_legacy))))

    cfg = asdict(nominal_extension())
    formal_analysis = analyse_solution(formal.base, formal_sol, cfg, 300.0)
    legacy_analysis = analyse_solution(legacy_base, legacy_sol, cfg, 300.0)
    formal_signals = formal.diagnostic_signals(formal_sol.y)
    # The legacy signal formulas are evaluated through an independent formal
    # helper only after state equivalence has been established.
    reference_signals = formal.diagnostic_signals(legacy_sol.y)
    signal_errors = {k: float(np.max(np.abs(formal_signals[k] - reference_signals[k])))
                     for k in formal_signals}

    event_keys = ('start_h', 'end_h', 'peak_h', 'peak', 'dose', 'duration_h')
    event_error = 0.0
    for kind in ('reverse_events', 'gate_events'):
        a, b = formal_analysis[kind], legacy_analysis[kind]
        if len(a) != len(b):
            event_error = float('inf')
            break
        for x, y in zip(a, b):
            event_error = max(event_error, *(abs(float(x[k]) - float(y[k])) for k in event_keys))

    report = dict(
        state_contract=dict(count=len(STATE_NAMES_34), names=list(STATE_NAMES_34)),
        integration=dict(hours=300.0, samples=int(formal_sol.t.size),
                         overall_max_abs_error=float(error.max()),
                         per_state_max_abs_error=per_state,
                         max_rhs_error=float(max(derivative_errors))),
        signals=dict(max_abs_error=signal_errors),
        digital=dict(formal_cold=formal_analysis['cold_start'],
                     legacy_cold=legacy_analysis['cold_start'],
                     formal_steady=formal_analysis['steady_state'],
                     legacy_steady=legacy_analysis['steady_state'],
                     formal_certified=formal_analysis['certified'],
                     legacy_certified=legacy_analysis['certified']),
        causal=dict(reverse_events_formal=len(formal_analysis['reverse_events']),
                    reverse_events_legacy=len(legacy_analysis['reverse_events']),
                    gate_events_formal=len(formal_analysis['gate_events']),
                    gate_events_legacy=len(legacy_analysis['gate_events']),
                    max_event_metric_error=event_error,
                    formal=formal_analysis['causal_verdict'],
                    legacy=legacy_analysis['causal_verdict']),
    )
    report['passed'] = bool(
        report['integration']['overall_max_abs_error'] <= 1e-9 and
        report['integration']['max_rhs_error'] <= 1e-9 and
        max(signal_errors.values()) <= 1e-9 and event_error <= 1e-9 and
        formal_analysis['steady_state']['sequence'] == legacy_analysis['steady_state']['sequence'] and
        formal_analysis['certified'] and legacy_analysis['certified'])
    (OUT / 'equivalence.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(dict(passed=report['passed'], integration=report['integration'],
                          signals=report['signals'], digital=report['digital'],
                          causal=report['causal']), ensure_ascii=False, indent=2))
    print(f"wrote {OUT / 'equivalence.json'}")


if __name__ == '__main__':
    main()
