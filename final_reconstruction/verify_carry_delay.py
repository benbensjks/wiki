"""Test a biochemical A0/F0 expression delay without changing ZENG values.

The baseline model treats A0 and F0 as directly produced mature proteins.  This
variant replaces only those two production terms by explicit mRNA -> immature
protein -> mature protein chains.  Four states are appended to the 30-state
two-bit system; A0/F0 themselves remain at indices 28/29.  The new half-lives
are explicitly uncalibrated Extension assumptions and are scanned here.
"""
from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from model import Extension, Model, ZENG, act
from verify_bit0_part2 import CAND, OUT
from verify_twobit_causal import analyse_solution


def _rates(model, mrna_half_life_min, maturation_half_life_min):
    lm = np.log(2) * 60 / mrna_half_life_min + model.growth
    kmat = np.log(2) * 60 / maturation_half_life_min
    lu = kmat + model.e.immature_decay_h + model.growth
    return lm, kmat, lu


def _expression(model, m, u, protein, target_mature_source, gamma, rates):
    lm, kmat, lu = rates
    transcript_source = (target_mature_source * lm * lu /
                         (model.e.translation_h * kmat))
    return np.asarray((transcript_source - lm * m,
                       model.e.translation_h * m - lu * u,
                       kmat * u - (gamma + model.growth) * protein))


def initial_state(model, rates):
    y30 = model.initial_state_twobit(cold=False)
    A0, F0 = y30[28], y30[29]
    source_A = ZENG['alpha_A'][0] * (1 - y30[16])
    source_F = ZENG['alpha_F'][0] * act(A0, ZENG['K_A'][0], ZENG['n_A'][0])
    lm, kmat, lu = rates

    def precursor_ss(source):
        m = source * lu / (model.e.translation_h * kmat)
        u = source / kmat
        return m, u

    mA, uA = precursor_ss(source_A)
    mF, uF = precursor_ss(source_F)
    return np.r_[y30, mA, uA, mF, uF]


def simulate(cfg, carry_mrna_half_life_min, carry_maturation_half_life_min,
             hours=300.0, sample_min=2.0, upstream_initial7=None,
             carry_repressor_maturation_half_life_min=None):
    model = Model(Extension(**cfg))
    repressor_mat = (carry_maturation_half_life_min
                     if carry_repressor_maturation_half_life_min is None
                     else carry_repressor_maturation_half_life_min)
    rates_A = _rates(model, carry_mrna_half_life_min, carry_maturation_half_life_min)
    rates_F = _rates(model, carry_mrna_half_life_min, repressor_mat)

    def rhs(t, y):
        y30 = y[:30]
        d30 = model.rhs_twobit(t, y30)
        A0, F0 = y30[28], y30[29]
        mA, uA, mF, uF = y[30:34]
        source_A = ZENG['alpha_A'][0] * (1 - y30[16])
        source_F = ZENG['alpha_F'][0] * act(A0, ZENG['K_A'][0], ZENG['n_A'][0])
        eA = _expression(model, mA, uA, A0, source_A, ZENG['gamma_A'][0], rates_A)
        eF = _expression(model, mF, uF, F0, source_F, ZENG['gamma_F'][0], rates_F)
        d30[28], d30[29] = eA[2], eF[2]
        return np.r_[d30, eA[:2], eF[:2]]

    ts = np.arange(0.0, hours + sample_min / 120.0, sample_min / 60.0)
    y0 = model.initial_state_twobit(cold=False)
    A0, F0 = y0[28], y0[29]
    source_A = ZENG['alpha_A'][0] * (1 - y0[16])
    source_F = ZENG['alpha_F'][0] * act(A0, ZENG['K_A'][0], ZENG['n_A'][0])

    def precursor_ss(source, rates):
        _, kmat, lu = rates
        return source * lu / (model.e.translation_h * kmat), source / kmat

    mA, uA = precursor_ss(source_A, rates_A)
    mF, uF = precursor_ss(source_F, rates_F)
    y0 = np.r_[y0, mA, uA, mF, uF]
    if upstream_initial7 is not None:
        upstream_initial7 = np.asarray(upstream_initial7, dtype=float)
        if upstream_initial7.shape != (7,):
            raise ValueError('upstream_initial7 must contain six oscillator states and C31 mRNA')
        y0[:7] = upstream_initial7
    sol = solve_ivp(rhs, (0, hours), y0, t_eval=ts,
                    method='DOP853', rtol=2e-7, atol=2e-9, max_step=2 / 60)
    if not sol.success:
        raise RuntimeError(sol.message)
    result = analyse_solution(model, sol, cfg, hours)
    result['carry_expression_extension'] = dict(
        added_states=['M_A0', 'A0_u', 'M_F0', 'F0_u'],
        mrna_half_life_min=float(carry_mrna_half_life_min),
        maturation_half_life_min=float(carry_maturation_half_life_min),
        activator_maturation_half_life_min=float(carry_maturation_half_life_min),
        repressor_maturation_half_life_min=float(repressor_mat),
        zeng_parameters_changed=False)
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--uM', type=float, default=6.0)
    ap.add_argument('--maturation-min', default='10,20,30,45,60,90,120')
    ap.add_argument('--mrna-min', type=float, default=2.0)
    ap.add_argument('--hours', type=float, default=300.0)
    ap.add_argument('--output', type=Path, default=OUT / 'carry_delay_scan.json')
    args = ap.parse_args()

    base = asdict(Extension())
    base.update(CAND)
    base['uM_per_au'] = args.uM
    rows = []
    args.output.parent.mkdir(parents=True, exist_ok=True)
    for maturation in (float(x) for x in args.maturation_min.split(',')):
        row = simulate(base, args.mrna_min, maturation, hours=args.hours)
        rows.append(row)
        args.output.write_text(json.dumps(rows, ensure_ascii=False, indent=2), encoding='utf-8')
        c = row['causal_verdict']
        print(f"uM={args.uM:g} tmat={maturation:g}min cold={row['cold_start']['passed']} "
              f"steady={row['steady_state']['passed']} one_to_one={c['exactly_one_gate_and_flip_per_late_reverse']} "
              f"order={c['causal_order_after_reverse_start']} certified={row['certified']} "
              f"seq={row['steady_state']['sequence']}", flush=True)
    print(f'wrote {args.output}')


if __name__ == '__main__':
    main()
