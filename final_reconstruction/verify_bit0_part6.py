"""Part 6: the A0/F0 -> bit1 carry interface, tested at the conversions where
bit0 dwells high enough for the AND gate to release at all.

Two-bit system = bit0 + bit1 + A0 + F0 (30 states). Nothing in ZENG is changed.
"""
from __future__ import annotations

import json
from dataclasses import asdict

import numpy as np
import pandas as pd
from scipy.signal import find_peaks

from model import Model, Extension, ZENG, ROOT
from verify_bit0_part2 import OUT, CAND, band_labels, clock_cycles


def H(x, K, n):
    x = np.maximum(x, 0.0)
    return x**n / (K**n + x**n)


def analyse(cfg, hours=180.0, sample_min=2.0):
    m = Model(Extension(**cfg))
    sol = m.simulate_twobit(hours=hours, sample_min=sample_min, max_step_min=2.0)
    t = sol.t
    S0, S1 = sol.y[6 + 10], sol.y[17 + 10]      # bit0_S, bit1_S
    A0, F0 = sol.y[28], sol.y[29]
    flux = np.array([m.flux(m._expand_twobit(sol.y[:, k])) for k in range(sol.y.shape[1])])
    ids = clock_cycles(t, flux)
    g0 = H(A0, ZENG['K_A'][0], ZENG['n_A'][0]) * (1 - H(F0, ZENG['K_F'][0], ZENG['n_F'][0]))
    lab0 = band_labels(t, S0, ids)
    lab1 = band_labels(t, S1, ids)
    seq = ''.join(l0['label'] + l1['label'] for l0, l1 in zip(lab0, lab1))
    # falling edges of bit0 (LR -> PB), measured on the state, not on a threshold crossing
    fall = [k for k in range(1, len(lab0)) if lab0[k - 1]['label'] == '1' and lab0[k]['label'] == '0']
    rises = [k for k in range(1, len(lab0)) if lab0[k - 1]['label'] == '0' and lab0[k]['label'] == '1']
    thr = 0.1 * (g0.max() - g0.min()) + g0.min()
    open_windows = []
    inside = False
    for k, v in enumerate(g0):
        if v > thr and not inside:
            start = k; inside = True
        elif v <= thr and inside:
            open_windows.append((start, k - 1)); inside = False
    if inside:
        open_windows.append((start, len(g0) - 1))
    big = [w for w in open_windows if (w[1] - w[0]) * (t[1] - t[0]) > 0.2]
    per_fall = []
    for k in fall:
        mid = (t[ids[k - 1]] + t[ids[k]]) / 2
        per_fall.append(dict(to_cycle=k, transition_h=float(mid),
                             gate_openings_in_previous_half_period=int(sum(
                                 1 for a, b in big if mid - 5.3 <= t[a] <= mid + 5.3)),
                             gate_openings_in_next_half_period=int(sum(
                                 1 for a, b in big if mid < t[a] <= mid + 10.6)),
                             bit1_label_before=lab1[k - 1]['label'], bit1_label_after=lab1[k]['label']))
    ok1 = lab1[-6:]
    bit1_alternates = ('x' not in ''.join(l['label'] for l in ok1)
                       and all(ok1[i + 1]['label'] != ok1[i]['label'] for i in range(len(ok1) - 1)))
    return dict(uM_per_au=cfg['uM_per_au'], add_growth=cfg['add_growth'],
                bit0_codes=''.join(l['label'] for l in lab0),
                bit0_commit=float(np.median([l['commitment'] for l in lab0[2:]])),
                bit0_dwell=[float(np.percentile(S0, 2)), float(np.percentile(S0, 98))],
                bit1_codes=''.join(l['label'] for l in lab1),
                bit1_commit=float(np.median([l['commitment'] for l in lab1[2:]])),
                bit1_dwell=[float(np.percentile(S1, 2)), float(np.percentile(S1, 98))],
                bit1_alternates=bit1_alternates,
                two_bit_sequence=seq,
                A0_range=[float(A0.min()), float(A0.max())], F0_range=[float(F0.min()), float(F0.max())],
                g0_range=[float(g0.min()), float(g0.max())], Int1_source_peak=float(18 * g0.max()),
                F0_below_threshold_fraction=float(np.mean(F0 < ZENG['K_F'][0])),
                gate_open_windows=[(round(float(t[a]), 2), round(float(t[b]), 2)) for a, b in big],
                n_fallings=len(fall), n_risings=len(rises), per_falling_edge=per_fall,
                bit1_S_range=[float(S1.min()), float(S1.max())])


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    base = asdict(Extension()); base.update(CAND)
    rows = []
    for u in (6.0, 8.0, 9.0, 10.0):
        cfg = base.copy(); cfg['uM_per_au'] = u
        r = analyse(cfg)
        rows.append(r)
        print(f"uM={u}: bit0={r['bit0_codes']} ({r['bit0_commit']:.2f}) bit1={r['bit1_codes']} "
              f"({r['bit1_commit']:.2f}) g0max={r['g0_range'][1]:.4f} Int1src={r['Int1_source_peak']:.3f} "
              f"F0<0.6 frac={r['F0_below_threshold_fraction']:.2f}", flush=True)
        (OUT / 'twobit_by_conversion.json').write_text(json.dumps(rows, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(rows, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
