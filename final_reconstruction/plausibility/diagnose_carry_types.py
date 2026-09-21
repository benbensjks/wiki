"""Re-aggregate the three-segment leak metric by CARRY TYPE.

Why
---
`leak_ratio = median(far_off_rev) / median(gate_on_rev)` is not a usable
statistic for this circuit.  The bit2 DNA flips on only every other carry, so
consecutive carries alternate between two distinct types:

    R-type : the carry performs the bit2 reverse action   -> gate_on_rev ~ 1.0
             and leaves far_off_rev ~ 0.0215, far_off_fwd ~ 2.4e-06
    F-type : the carry performs the forward flip instead  -> gate_on_rev ~ 0
             and leaves far_off_rev ~ 0.0023, far_off_fwd ~ 6.5e-03

Measured on the eight initial states (plausibility/eight_initial_diagnostic/):
the per-window pattern is exactly FRFRFR... , and over a 600 h horizon an
initial state whose phase makes the train end on an F gives 15 windows
(8 F + 7 R) while the balanced case gives 14 (7 F + 7 R).  A median over an ODD
number of windows therefore lands ON a carry type instead of between the two,
and the ratio moves by four orders of magnitude - 3.5e+04 for an 8-F set against
0.0121 for an 8-R set - with no change whatsoever in the circuit.

Consequence: the leak ratio must not be compared across runs with different
window counts, and the earlier n=5 vs n=6 leak interval was computed from
medians whose carry-type composition was never checked.  This script measures
that composition and replaces the statistic with a parity-robust one:

    pair consecutive windows (one R + one F = one full 8-read super-period) and
    integrate gate_on_rev, tail_rev, far_off_rev and far_off_fwd over the pair.

The paired ratio sum(far_off_rev)/sum(gate_on_rev) no longer depends on where the
600 h window happens to cut the carry train, because every pair contains exactly
one of each type.

    python diagnose_carry_types.py --n 5 --hours 600
    python diagnose_carry_types.py --n 6 --hours 600
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plausibility_common import (OUT, frozen_carry0, frozen_extension,  # noqa: E402
                                 gate_windows, leak_report)
from model_threebit51 import ThreeBit51Model, ThreeBitCarryParameters  # noqa: E402
from model_twobit34 import CarryExpressionParameters  # noqa: E402

R_THRESHOLD = 1e-3      # gate_on_rev above this = the carry did the bit2 reverse action


def build(n, mrna=2.0, mat=32.5):
    carry1 = CarryExpressionParameters(mrna_half_life_min=mrna,
                                       activator_maturation_half_life_min=mat,
                                       repressor_maturation_half_life_min=mat)
    return ThreeBit51Model(extension=frozen_extension(),
                           carry=ThreeBitCarryParameters(carry0=frozen_carry0(),
                                                         carry1=carry1),
                           n_A1_gate=n)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=float, required=True)
    ap.add_argument('--hours', type=float, default=600.0)
    ap.add_argument('--sample-min', type=float, default=2.0)
    ap.add_argument('--max-step-min', type=float, default=2.0)
    args = ap.parse_args()

    model = build(args.n)
    sol = model.simulate(hours=args.hours, sample_min=args.sample_min,
                         max_step_min=args.max_step_min)
    t, y = sol.t, sol.y
    sig = model.diagnostic_signals(y)
    win = gate_windows(t, sig['g1'])
    rep = leak_report(t, sig['g1'], sig['J_rev2'], sig['J_fwd2'], y[36], win)

    tag = f'n{args.n:g}'
    outdir = OUT / 'carry_type_diagnostic'
    outdir.mkdir(parents=True, exist_ok=True)

    rows = []
    for k, r in enumerate(rep[0.05]['per_window']):
        rows.append(dict(k=k, gate_on_rev=r['gate_on_rev'], tail_rev=r['tail_rev'],
                         far_off_rev=r['far_off_rev'], far_off_fwd=r['far_off_fwd'],
                         gate_on_h=r['gate_on_h'], tail_h=r['tail_h'],
                         far_off_h=r['far_off_h'], evaluable=r['far_off_evaluable'],
                         type='R' if r['gate_on_rev'] > R_THRESHOLD else 'F'))
    tbl = pd.DataFrame(rows)
    tbl.to_csv(outdir / f'windows_{tag}.csv', index=False, encoding='utf-8')

    nR = int((tbl.type == 'R').sum())
    nF = int((tbl.type == 'F').sum())
    print(f'=== n_A1_gate = {args.n:g}  600 h ===')
    print(f'windows {len(tbl)}   R-type {nR}   F-type {nF}   '
          f'pattern {"".join(tbl.type)}')
    print(f'window count is {"EVEN (balanced)" if len(tbl) % 2 == 0 else "ODD (unbalanced -> raw ratio unreliable)"}')

    print('\n--- per-window table ---')
    print(tbl.to_string(index=False))

    per_type = {}
    for ty in ('R', 'F'):
        sub = tbl[tbl.type == ty]
        if sub.empty:
            per_type[ty] = None
            continue
        per_type[ty] = dict(
            windows=int(len(sub)),
            gate_on_rev_median=float(sub.gate_on_rev.median()),
            gate_on_rev_min=float(sub.gate_on_rev.min()),
            gate_on_rev_max=float(sub.gate_on_rev.max()),
            far_off_rev_median=float(sub.far_off_rev.median()),
            far_off_fwd_median=float(sub.far_off_fwd.median()),
            tail_rev_median=float(sub.tail_rev.median()))
    print('\n--- per carry type (medians) ---')
    print(json.dumps(per_type, ensure_ascii=False, indent=2))

    # parity-robust statistic: consecutive pairs = one full 8-read super-period
    pairs = []
    for i in range(0, len(tbl) - 1, 2):
        a, b = tbl.iloc[i], tbl.iloc[i + 1]
        if not (a.evaluable and b.evaluable):
            continue
        g = a.gate_on_rev + b.gate_on_rev
        pairs.append(dict(pair=i // 2, types=f'{a.type}{b.type}',
                          gate_on_rev=g, tail_rev=a.tail_rev + b.tail_rev,
                          far_off_rev=a.far_off_rev + b.far_off_rev,
                          far_off_fwd=a.far_off_fwd + b.far_off_fwd,
                          ratio=(a.far_off_rev + b.far_off_rev) / g if g else None))
    pdf = pd.DataFrame(pairs)
    pdf.to_csv(outdir / f'paired_{tag}.csv', index=False, encoding='utf-8')
    print('\n--- paired (one R + one F = one 8-read super-period) ---')
    print(pdf.to_string(index=False))

    paired = dict(
        n_pairs=int(len(pdf)),
        all_pairs_are_one_of_each_type=bool(len(pdf) and (pdf.types == 'RF').all()
                                           or len(pdf) and (pdf.types == 'FR').all()),
        gate_on_rev_median=float(pdf.gate_on_rev.median()) if len(pdf) else None,
        far_off_rev_median=float(pdf.far_off_rev.median()) if len(pdf) else None,
        far_off_fwd_median=float(pdf.far_off_fwd.median()) if len(pdf) else None,
        tail_rev_median=float(pdf.tail_rev.median()) if len(pdf) else None,
        ratio_median=float(pdf.ratio.median()) if len(pdf) else None,
        ratio_min=float(pdf.ratio.min()) if len(pdf) else None,
        ratio_max=float(pdf.ratio.max()) if len(pdf) else None)
    print('\n--- parity-robust paired aggregate ---')
    print(json.dumps(paired, ensure_ascii=False, indent=2))

    result = dict(n_A1_gate=args.n, hours=args.hours, windows=len(tbl),
                  n_R=nR, n_F=nF, pattern=''.join(tbl.type),
                  parity_balanced=bool(len(tbl) % 2 == 0),
                  raw_leak_ratio_5pct=rep[0.05]['leak_ratio'],
                  raw_gate_on_rev_5pct=rep[0.05]['gate_on_rev'],
                  raw_far_off_rev_5pct=rep[0.05]['far_off_rev'],
                  raw_far_off_fwd_5pct=rep[0.05]['far_off_fwd'],
                  per_type=per_type, paired=paired,
                  note=('the raw median-of-medians ratio is parity-unstable because consecutive '
                        'carries alternate between two types; the paired statistic integrates one '
                        'full 8-read super-period and is the comparable quantity'))
    (outdir / f'summary_{tag}.json').write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding='utf-8')
    print(f'\nwrote {outdir}')


if __name__ == '__main__':
    main()
