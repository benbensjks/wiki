"""Part 7: long verification of the two-bit readout at the conversion where the
carry gate actually opens. Read point = flux trough (dwell midpoint), the same
read point the archived decoder uses.
"""
from __future__ import annotations

import json
from dataclasses import asdict

import numpy as np
import pandas as pd
from scipy.signal import find_peaks

from model import Model, Extension, ZENG, ROOT
from verify_bit0_part2 import OUT, CAND

H = lambda x, K, n: np.maximum(x, 0.0)**n / (K**n + np.maximum(x, 0.0)**n)


def read_two_bit(cfg, hours=300.0):
    m = Model(Extension(**cfg))
    sol = m.simulate_twobit(hours=hours, sample_min=2.0, max_step_min=2.0)
    t = sol.t
    S0, S1 = sol.y[16], sol.y[27]
    A0, F0 = sol.y[28], sol.y[29]
    flux = np.array([m.flux(m._expand_twobit(sol.y[:, k])) for k in range(sol.y.shape[1])])
    pk, _ = find_peaks(flux, prominence=max(1.0, .1 * np.ptp(flux)), distance=150)
    troughs = [a + int(np.argmin(flux[a:b])) for a, b in zip(pk[:-1], pk[1:])]
    samples = []
    for k in troughs:
        b0 = 1 if S0[k] >= 0.5 else 0
        b1 = 1 if S1[k] >= 0.5 else 0
        samples.append(dict(t=float(t[k]), S0=float(S0[k]), S1=float(S1[k]), b0=b0, b1=b1,
                            value=2 * b1 + b0))
    vals = [s['value'] for s in samples]
    counts = [v for v in vals if 0 <= v <= 3]
    increments_ok = all((counts[i + 1] - counts[i]) % 4 == 1 for i in range(len(counts) - 1))
    # one gate opening per bit0 LR -> PB interval, measured on the state
    g0 = H(A0, ZENG['K_A'][0], ZENG['n_A'][0]) * (1 - H(F0, ZENG['K_F'][0], ZENG['n_F'][0]))
    thr = 0.1 * g0.max()
    opens, inside = [], False
    for k, v in enumerate(g0):
        if v > thr and not inside:
            start = k; inside = True
        elif v <= thr and inside:
            opens.append((t[start], t[k])); inside = False
    fall = np.flatnonzero((S0[:-1] >= .5) & (S0[1:] < .5))
    per_fall = [int(sum(1 for a, b in opens if 0 <= t[i] - a <= 12)) for i in fall]
    return dict(uM_per_au=cfg['uM_per_au'], hours=hours, reads=len(samples),
                sequence=''.join(str(v) for v in vals),
                samples=samples,
                increments_correct=bool(increments_ok),
                distinct_values=sorted(set(vals)),
                bit0_troughs=[round(s['S0'], 3) for s in samples],
                bit1_troughs=[round(s['S1'], 3) for s in samples],
                gate_openings=len(opens),
                openings_per_bit0_falling_edge=per_fall,
                exactly_one_per_falling_edge=bool(per_fall and all(n == 1 for n in per_fall)),
                g0_peak=float(g0.max()), Int1_source_peak=float(18 * g0.max()),
                bit1_flip_lead_h=[float(t[i]) for i in fall])


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    base = asdict(Extension()); base.update(CAND)
    rows = []
    for u in (5.0, 6.0, 7.0, 8.0):
        cfg = base.copy(); cfg['uM_per_au'] = u
        r = read_two_bit(cfg)
        rows.append(r)
        print(f"uM={u}: reads={r['reads']} seq={r['sequence']} increments_ok={r['increments_correct']} "
              f"gate_openings={r['gate_openings']} one_per_fall={r['exactly_one_per_falling_edge']} "
              f"g0max={r['g0_peak']:.3f}", flush=True)
        (OUT / 'twobit_readout.json').write_text(json.dumps(rows, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(rows, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
