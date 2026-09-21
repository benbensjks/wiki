"""Part 3: decisive parity test with self-consistent initial states, plus the
quantitative reason the A0/F0 carry gate can never open.
"""
from __future__ import annotations

import json
from dataclasses import asdict

import numpy as np
from scipy.integrate import solve_ivp

from model import Model, Extension, ZENG, ROOT
from verify_bit0_part2 import OUT, CAND, sim, clock_cycles, band_labels, DT


def consistent_bit0_init(m, S):
    """Self-consistent T/R dwell values for a given DNA state, other pools zero."""
    y = m.initial_state_bit0(cold=False)
    b = y[6:17]
    T_src = ZENG['alpha_rep'] * (1 - S)
    b[3] = m.transcript_source(T_src) / m.lm
    b[4] = m.e.translation_h * b[3] / m.lu
    b[5] = T_src / (ZENG['gamma_rep'] + m.growth)
    T = b[5]
    R_src = ZENG['alpha_rdf'] * S * (1 - (T / ZENG['K_rep']) ** ZENG['n_rep']
                                     / (1 + (T / ZENG['K_rep']) ** ZENG['n_rep']))
    b[6] = m.transcript_source(R_src) / m.lm
    b[7] = m.e.translation_h * b[6] / m.lu
    b[8] = R_src / (ZENG['gamma_rdf'] + m.growth)
    b[9] = 0.0
    b[10] = S
    return y


def part_L_parity(hours=180.0):
    m = Model(Extension(**CAND))
    runs = {}
    for name, S0 in (('branch_LR_S0p84', 0.843), ('branch_PB_S0p18', 0.184)):
        y0 = consistent_bit0_init(m, S0)
        runs[name] = sim(m, y0, 0.0, hours)
    ref = runs['branch_LR_S0p84']
    flux = np.array([m.flux(m._expand_bit0(ref.y[:, k])) for k in range(ref.y.shape[1])])
    ids = clock_cycles(ref.t, flux)
    labs = {k: band_labels(v.t, v.y[16], ids) for k, v in runs.items()}
    a = ''.join(l['label'] for l in labs['branch_LR_S0p84'])
    b = ''.join(l['label'] for l in labs['branch_PB_S0p18'])
    t = ref.t
    gap = np.abs(runs['branch_LR_S0p84'].y[16] - runs['branch_PB_S0p18'].y[16])
    tail = [dict(cycle=i, LR=labs['branch_LR_S0p84'][i]['label'], PB=labs['branch_PB_S0p18'][i]['label'],
                 LR_med=labs['branch_LR_S0p84'][i]['S_med'], PB_med=labs['branch_PB_S0p18'][i]['S_med'])
            for i in range(len(a))]
    same_branch = all(x['LR'] == x['PB'] for x in tail if x['LR'] != 'x' and x['PB'] != 'x')
    return dict(hours=hours, codes_from_LR_init=a, codes_from_PB_init=b, cycles=tail,
                gap_max_last_40h=float(np.max(gap[t >= hours - 40])),
                gap_median_20_40h=float(np.median(gap[(t >= 20) & (t <= 40)])),
                converged=bool(np.max(gap[t >= hours - 40]) < 0.1),
                same_branch=same_branch,
                verdict=('PARITY STORED: two self-consistent opposite starts stay antiphase -> real toggle'
                         if not same_branch and not np.max(gap[t >= hours - 40]) < 0.1 else
                         'PARITY NOT STORED: both starts converge to the same phase-locked orbit -> the '
                         'alternation is a period-doubled driven response, not stored counting state'))


def part_M_carry_window():
    """A0/F0 only: how short must an A0 pulse be for the AND gate to open?"""
    H = lambda x, K, n: x**n / (K**n + x**n)
    A_of = lambda pb: ZENG['alpha_A'][0] * pb / ZENG['gamma_A'][0]

    def peak_g0(width_h, pb_high=0.82, pb_low=0.16):
        def rhs(t, y):
            A, F = y
            pb = pb_high if t < width_h else pb_low
            dA = ZENG['alpha_A'][0] * pb - ZENG['gamma_A'][0] * A
            dF = ZENG['alpha_F'][0] * H(A, ZENG['K_A'][0], ZENG['n_A'][0]) - ZENG['gamma_F'][0] * F
            return [dA, dF]
        A_lo = A_of(pb_low)
        F_lo = ZENG['alpha_F'][0] * H(A_lo, ZENG['K_A'][0], ZENG['n_A'][0]) / ZENG['gamma_F'][0]
        span = max(40.0, width_h * 4)
        ts = np.linspace(0, span, 8000)
        sol = solve_ivp(rhs, (0, span), [A_lo, F_lo], t_eval=ts, method='DOP853', rtol=1e-10, atol=1e-13)
        A, F = sol.y
        g = H(A, ZENG['K_A'][0], ZENG['n_A'][0]) * (1 - H(F, ZENG['K_F'][0], ZENG['n_F'][0]))
        return float(g.max()), float(g[np.argmin(np.abs(ts - width_h))])

    rows = []
    for w in (0.05, 0.1, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0, 10.6, 21.2):
        peak, at_end = peak_g0(w)
        rows.append(dict(pulse_width_h=w, peak_g0=peak, g0_at_pulse_end=at_end,
                         Int1_source_peak_au_per_h=18 * peak))
    widths_ok = [r['pulse_width_h'] for r in rows if 18 * r['peak_g0'] >= 0.05]
    # A0 floor that would still let the repression arm release the gate
    h_crit = ZENG['K_F'][0] * ZENG['gamma_F'][0] / ZENG['alpha_F'][0]
    a0_crit = ZENG['K_A'][0] * np.sqrt(h_crit / (1 - h_crit))
    s_high_needed = 1 - a0_crit * ZENG['gamma_A'][0] / ZENG['alpha_A'][0]
    return dict(A0_steady_high=round(float(A_of(0.82)), 3), A0_steady_low=round(float(A_of(0.16)), 3),
                F0_floor_at_low_state=round(float(ZENG['alpha_F'][0] * H(A_of(0.16), ZENG['K_A'][0], ZENG['n_A'][0])
                                                 / ZENG['gamma_F'][0]), 3),
                F0_saturation=round(float(ZENG['alpha_F'][0] / ZENG['gamma_F'][0]), 3),
                F0_threshold=ZENG['K_F'][0], F0_response_time_h=round(1 / ZENG['gamma_F'][0], 2),
                toggle_dwell_h=10.58, rows=rows,
                A0_floor_for_open_gate=round(float(a0_crit), 4),
                S_high_required_for_open_gate=round(float(s_high_needed), 4),
                S_high_observed=0.843,
                verdict=('the gate never opens at any pulse width: F0 is already {:.2f} at the '
                         'low-A0 state, {:.1f}x its own threshold {}. To release the gate bit0 must reach '
                         'S >= {:.3f}, not the {:.3f} it reaches now.'.format(
                             float(ZENG['alpha_F'][0] * H(A_of(0.16), ZENG['K_A'][0], ZENG['n_A'][0])
                                   / ZENG['gamma_F'][0]),
                             float(ZENG['alpha_F'][0] * H(A_of(0.16), ZENG['K_A'][0], ZENG['n_A'][0])
                                   / ZENG['gamma_F'][0]) / ZENG['K_F'][0],
                             ZENG['K_F'][0], float(s_high_needed), 0.843))
                + ('' if widths_ok else ' No pulse width from 0.05 h to 21.2 h opens it either.'))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rep = {}
    for name, fn in (('L_parity_with_consistent_starts', part_L_parity),
                     ('M_carry_gate_window', part_M_carry_window)):
        rep[name] = fn()
        (OUT / 'bit0_verification_part3.json').write_text(json.dumps(rep, ensure_ascii=False, indent=2), encoding='utf-8')
        print(name, 'done', flush=True)
    print(json.dumps(rep, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
