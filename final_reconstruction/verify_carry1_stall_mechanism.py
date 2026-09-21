"""Step-0 mechanistic check of the bit2 stall at the best stage-1 point.

Point: carry1 mRNA 4 min, A1/F1 maturation 60 min, 34-state prefix and ZENG
frozen.  Question: is the stall a genuine forward/reverse balance, and is it
pinned by the disappearance of the RDF2 pool (the Kinh/(Kinh+R2) suppression)
rather than by the pulse dose?

Two tests:

1.  At the deepest point of each carry excursion, evaluate the forward and
    reverse rates, the two Hill terms and Kinh/(Kinh+R2).  A balance shows up as
    forward ~ reverse with |dS2/dt| ~ 0.
2.  Thought experiment: hold R2 fixed at a value R and solve
        k_fwd (1-S) H(I2;K_D_int2,2) Kinh/(Kinh+R) = k_rev S H(q I2 R; K_C, 2)
    for S.  The resulting S*(R) curve shows how deep bit2 can be written for a
    given RDF2 background; S*(R->0) = k_fwd/(k_fwd+k_rev) = 7/12 = 0.583 is the
    floor claim under test.

    python verify_carry1_stall_mechanism.py
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

from model import ROOT, ZENG
from model_threebit51 import ThreeBit51Model, ThreeBitCarryParameters
from model_twobit34 import CarryExpressionParameters
from verify_twobit_causal import GATE_ON, MIN_GATE_DOSE_H, _segments

MRNA, MAT = 4.0, 60.0
HOURS = 600.0
SKIP_H = 100.0
K_FWD, K_REV = ZENG['k_fwd'], ZENG['k_rev']
KD_INT2, N_INT2 = ZENG['K_D_int'][2], 2.0
KINH = ZENG['K_inh']


def hill(x, K, n):
    x = np.maximum(np.asarray(x, dtype=float), 0.0)
    return x ** n / (K ** n + x ** n)


def build():
    c1 = CarryExpressionParameters(mrna_half_life_min=MRNA,
                                   activator_maturation_half_life_min=MAT,
                                   repressor_maturation_half_life_min=MAT)
    return ThreeBit51Model(carry=ThreeBitCarryParameters(carry1=c1))


def balance_residual(S, I2, R2, Kc, q):
    fwd = K_FWD * (1 - S) * hill(I2, KD_INT2, N_INT2) * KINH / (KINH + max(R2, 0.0))
    rev = K_REV * S * hill(q * I2 * R2, Kc, 2)
    return rev - fwd


def main():
    global MRNA, MAT
    ap = argparse.ArgumentParser()
    ap.add_argument('--mrna', type=float, default=MRNA)
    ap.add_argument('--mat', type=float, default=MAT)
    ap.add_argument('--tag', default=None)
    ap.add_argument('--hours', type=float, default=HOURS)
    args = ap.parse_args()
    MRNA, MAT = args.mrna, args.mat
    tag = args.tag or f'mrna{MRNA:g}_mat{MAT:g}'.replace('.', 'p')
    out = ROOT / 'threebit51_results' / f'carry1_stall_{tag}'
    out.mkdir(parents=True, exist_ok=True)
    model = build()
    sol = model.simulate(hours=args.hours, sample_min=2.0, max_step_min=2.0)
    t = sol.t
    sig = model.diagnostic_signals(sol.y)
    S2, I2, R2 = sig['S2'], sol.y[36], sol.y[42]
    C2 = sol.y[43]
    fwd, rev = sig['J_fwd2'], sig['J_rev2']
    g1 = sig['g1']
    q = model.base.q
    Kc = model.base.K_complex

    gates = [g for g in _segments(t, g1, GATE_ON, min_area=MIN_GATE_DOSE_H)
             if g['start_h'] > SKIP_H]
    cycles = []
    for g in gates:
        a = int(np.searchsorted(t, g['start_h'] - 1.0))
        b = int(np.searchsorted(t, g['start_h'] + 10.0))
        sl = slice(a, min(b, len(t) - 1))
        k = a + int(np.argmin(S2[sl]))
        cycles.append(dict(
            gate_start_h=g['start_h'], gate_duration_h=g['duration_h'],
            gate_peak=g['peak'], s2_min=float(S2[k]), s2_min_h=float(t[k]),
            s2_before=float(S2[a]),
            forward_at_min=float(fwd[k]), reverse_at_min=float(rev[k]),
            balance_ratio=float(rev[k] / fwd[k]) if fwd[k] > 0 else None,
            ds2dt_at_min=float((S2[k + 1] - S2[k - 1]) / (t[k + 1] - t[k - 1])),
            I2_at_min=float(I2[k]), R2_at_min=float(R2[k]), C2_at_min=float(C2[k]),
            H_I2=float(hill(I2[k], KD_INT2, N_INT2)),
            H_C2=float(hill(C2[k], Kc, 2)),
            kinh_factor=float(KINH / (KINH + max(R2[k], 0.0))),
            R2_before=float(R2[a]), I2_before=float(I2[a])))

    # S*(R2) thought experiment at the measured pulse peak
    I2_pulse = float(np.percentile(I2[t > SKIP_H], 99))
    grid = np.logspace(-3, 1.6, 160)
    s_star, s_fail = [], []
    for R in grid:
        try:
            s_star.append(brentq(balance_residual, 1e-9, 1 - 1e-9, args=(I2_pulse, R, Kc, q),
                                 xtol=1e-10))
        except ValueError:
            s_star.append(np.nan)
    s_star = np.array(s_star)
    floor = K_FWD / (K_FWD + K_REV)

    observed = np.array([c['s2_min'] for c in cycles])
    r_at_min = np.array([c['R2_at_min'] for c in cycles])
    summary = dict(point=dict(carry1_mrna_min=MRNA, carry1_maturation_min=MAT),
                   hours=HOURS, gates_after_skip=len(gates),
                   pulse_peak_I2=I2_pulse, K_complex=float(Kc), q=float(q),
                   floor_if_R2_to_zero=float(floor),
                   observed_s2_min_median=float(np.median(observed)),
                   observed_s2_min_min=float(observed.min()),
                   R2_at_min_median=float(np.median(r_at_min)),
                   median_balance_ratio=float(np.nanmedian(
                       [c['balance_ratio'] for c in cycles if c['balance_ratio']])),
                   median_abs_ds2dt=float(np.median([abs(c['ds2dt_at_min']) for c in cycles])),
                   median_H_I2=float(np.median([c['H_I2'] for c in cycles])),
                   median_H_C2=float(np.median([c['H_C2'] for c in cycles])),
                   median_kinh_factor=float(np.median([c['kinh_factor'] for c in cycles])),
                   S_star_at_observed_R2=float(np.interp(np.median(r_at_min), grid, s_star)),
                   cycles=cycles)

    # ------------------------------------------------ figure: one late cycle
    g = gates[-3] if len(gates) >= 3 else gates[-1]
    w = slice(int(np.searchsorted(t, g['start_h'] - 6.0)),
              int(np.searchsorted(t, g['start_h'] + 16.0)))
    fig, ax = plt.subplots(4, 1, figsize=(11, 13), sharex=True, constrained_layout=True)
    ax[0].plot(t[w], S2[w], c='#0072B2')
    ax[0].axhspan(0.0, 0.30, color='#0072B2', alpha=.12)
    ax[0].axhspan(0.70, 1.0, color='#D55E00', alpha=.12)
    ax[0].axhline(floor, color='k', ls='--', lw=1,
                  label=f'floor k_fwd/(k_fwd+k_rev) = {floor:.3f}')
    ax[0].set_ylabel('S2'); ax[0].legend(fontsize=8); ax[0].set_title(
        f'bit2 at carry1 mRNA {MRNA:g} min / maturation {MAT:g} min - one carry window')
    ax[1].plot(t[w], I2[w], label='I2 (Int2)'); ax[1].plot(t[w], R2[w], label='R2 (RDF2)')
    ax[1].plot(t[w], C2[w] * 20, label='C2 x20 (Int2 x RDF2 complex)')
    ax[1].set_ylabel('a.u.'); ax[1].legend(fontsize=8)
    ax[2].plot(t[w], fwd[w], label='forward J_fwd2'); ax[2].plot(t[w], rev[w], label='reverse J_rev2')
    ax[2].set_ylabel('h$^{-1}$'); ax[2].legend(fontsize=8)
    ax[3].plot(t[w], g1[w], c='#009E73', label='g1 carry gate')
    ax[3].plot(t[w], hill(I2[w], KD_INT2, N_INT2), '--', label='H(I2;0.6,2)')
    ax[3].plot(t[w], hill(C2[w], Kc, 2), ':', label='H(C2;K_C,2)')
    ax[3].plot(t[w], KINH / (KINH + R2[w]), '-.', label='Kinh/(Kinh+R2)')
    ax[3].set_ylabel('dimensionless'); ax[3].set_xlabel('time (h)'); ax[3].legend(fontsize=8)
    for a in ax:
        a.axvline(g['start_h'], color='grey', lw=.8)
        for c in cycles:
            a.axvline(c['s2_min_h'], color='red', lw=.5, alpha=.5)
        a.grid(alpha=.25)
    for ext in ('png', 'pdf'):
        fig.savefig(out /f'stall_cycle.{ext}', dpi=200, bbox_inches='tight', facecolor='white')
    plt.close(fig)

    # ------------------------------------------------ figure: S*(R2) curve
    fig2, ax2 = plt.subplots(1, 2, figsize=(12, 4.6), constrained_layout=True)
    ax2[0].semilogx(grid, s_star, c='#0072B2')
    ax2[0].axhline(0.30, color='k', ls=':', label='lower read band 0.30')
    ax2[0].axhline(floor, color='k', ls='--', label=f'R2->0 floor {floor:.3f}')
    ax2[0].scatter(r_at_min, observed, s=14, c='red', zorder=5,
                   label='observed per carry cycle')
    ax2[0].set_xlabel('RDF2 background R2 (a.u.)'); ax2[0].set_ylabel('balanced S2*')
    ax2[0].set_title('how deep can bit2 be written vs RDF2 background')
    ax2[0].legend(fontsize=8); ax2[0].grid(alpha=.25)
    ax2[1].hist(observed, bins=12, color='#D55E00', alpha=.85)
    ax2[1].axvline(0.30, color='k', ls=':', label='lower read band 0.30')
    ax2[1].axvline(floor, color='k', ls='--', label=f'floor {floor:.3f}')
    ax2[1].set_xlabel('per-cycle S2 minimum'); ax2[1].set_ylabel('carry cycles')
    ax2[1].set_title('observed stall depth'); ax2[1].legend(fontsize=8); ax2[1].grid(alpha=.25)
    for ext in ('png', 'pdf'):
        fig2.savefig(out /f'stall_floor.{ext}', dpi=200, bbox_inches='tight', facecolor='white')
    plt.close(fig2)

    (out /'stall_mechanism.json').write_text(
        json.dumps(dict(summary=summary, s_star_curve=dict(R2=grid.tolist(),
                                                           S_star=s_star.tolist())),
                   ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(summary, ensure_ascii=False, indent=2, default=str))


if __name__ == '__main__':
    main()

