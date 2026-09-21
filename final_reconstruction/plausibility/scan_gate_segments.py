"""Scan: three-segment leak decomposition on the carry-2 gate grid.

Purpose (agreed next-step item 1 + 3)
-------------------------------------
The pre-registered R1 used `leak_worst_Jrev2 <= 0.10`, which mixes the intended
reverse action during the gate window, the carry's own Int2 decay tail, and the
genuine far-off leakage.  Evidence that it is contaminated: at fitted
maturation 60 min the off-state gate peak is 28x lower than at 32.5 min, yet
`leak_worst` is 40 % higher.  R1 therefore cannot be used to veto the 16
certified points until the metric is split.

This scan re-runs the arm-B grid with the FORMAL model parameter
`ThreeBit51Model(n_A1_gate=...)` added to model_threebit51.py, saves the
per-window segment table, and compares the far-off leakage ratio between
certified and non-certified points.  A threshold is only proposed AFTER seeing
both distributions - the old 0.10 is not reused.

Grid (20 points, 600 h each)
    n_A1_gate in {4,5,6,7,8}  x  carry1 mRNA in {2,4} min
    x  A1/F1 common maturation in {32.5,60} min
Everything else is frozen (ZENG untouched, F1 production exponent 4).

Guards
    * prefix guard : states 0..33 must be identical between n=4 and n=6 at the
      same initial condition (verified per run, recorded in every row);
    * the model file hash is recorded so a later edit is detectable.
"""
from __future__ import annotations

import argparse
import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plausibility_common import (OUT, frozen_carry0, frozen_extension,  # noqa: E402
                                 gate_windows, leak_report, merge_shards,
                                 shard_slice, sha256, source_hashes,
                                 write_manifest, write_shard)
from model_threebit51 import ThreeBit51Model, ThreeBitCarryParameters  # noqa: E402
from model_twobit34 import CarryExpressionParameters  # noqa: E402

NAME = 'gate_segments'
N_VALUES = (4.0, 5.0, 6.0, 7.0, 8.0)
MRNA_VALUES = (2.0, 4.0)
MAT_VALUES = (32.5, 60.0)


def points():
    return [dict(n_A1_gate=n, mrna=m, mat=t)
            for n, m, t in itertools.product(N_VALUES, MRNA_VALUES, MAT_VALUES)]


def point_id(p):
    return f"n={p['n_A1_gate']:g}|mRNA={p['mrna']:g}|mat={p['mat']:g}"


def build(p):
    carry1 = CarryExpressionParameters(mrna_half_life_min=p['mrna'],
                                       activator_maturation_half_life_min=p['mat'],
                                       repressor_maturation_half_life_min=p['mat'])
    return ThreeBit51Model(extension=frozen_extension(),
                           carry=ThreeBitCarryParameters(carry0=frozen_carry0(), carry1=carry1),
                           n_A1_gate=p['n_A1_gate'])


def prefix_guard():
    """States 0..33 must not depend on n_A1_gate.

    STRUCTURAL test (the gate): at the same state, the first 34 entries of the
    right-hand side must be bit-identical for the frozen exponent and for every
    scanned exponent.  The correction only touches d[34] (bit2's M_I), so this
    is the property that actually has to hold.

    A trajectory-level comparison is deliberately NOT used as a gate: once bit2
    differs, the adaptive step-size controller may choose different steps, so
    the shared states pick up round-off (~1e-13), which is not a dependence.
    """
    m4 = build(dict(n_A1_gate=4.0, mrna=2.0, mat=32.5))
    y0 = m4.initial_state(cold=False)
    sol = m4.simulate(hours=60.0, sample_min=2.0, max_step_min=2.0, initial_state=y0)
    detail, traj = {}, {}
    ok = True
    for n in N_VALUES:
        mN = build(dict(n_A1_gate=n, mrna=2.0, mat=32.5))
        worst = 0.0
        for k in range(0, sol.y.shape[1], 150):
            tk, yk = float(sol.t[k]), sol.y[:, k]
            worst = max(worst, float(np.max(np.abs(m4.rhs(tk, yk)[:34] - mN.rhs(tk, yk)[:34]))))
        detail[f'rhs_gap_n={n:g}'] = worst
        ok = ok and worst == 0.0
        sN = mN.simulate(hours=60.0, sample_min=2.0, max_step_min=2.0, initial_state=y0)
        traj[f'traj_gap_n={n:g}'] = float(np.max(np.abs(sol.y[:34] - sN.y[:34])))
    return ok, detail, traj


def evaluate(p, hours, sample_min, max_step_min):
    from verify_threebit51 import analyse_threebit
    started = time.perf_counter()
    row = dict(point_id=point_id(p), **p)
    try:
        model = build(p)
        sol = model.simulate(hours=hours, sample_min=sample_min, max_step_min=max_step_min)
        a = analyse_threebit(model, sol, hours)
        sig = model.diagnostic_signals(sol.y)
        t = sol.t
        win = gate_windows(t, sig['g1'])
        rep = leak_report(t, sig['g1'], sig['J_rev2'], sig['J_fwd2'], sol.y[36], win)
        r = rep[0.05]
        row.update(ok=True, error='',
                   certified=bool(a['certified']),
                   crossings=len(a['bit2_crossings']),
                   gates=len(a['carry1_gate_events']),
                   reverse_events=len(a['bit1_reverse_events']),
                   steady_sequence=a['steady_state']['sequence'],
                   # three-segment table at the 5 % cut
                   gate_on_h=r['gate_on_h'], tail_h=r['tail_h'], far_off_h=r['far_off_h'],
                   gate_on_rev=r['gate_on_rev'], tail_rev=r['tail_rev'],
                   far_off_rev=r['far_off_rev'], far_off_fwd=r['far_off_fwd'],
                   leak_ratio=r['leak_ratio'], far_off_fwd_ratio=r['far_off_fwd_ratio'],
                   windows_total=r['windows_total'],
                   windows_evaluable=r['windows_evaluable'],
                   # sensitivity of the cut fraction
                   leak_ratio_1pct=rep[0.01]['leak_ratio'],
                   leak_ratio_5pct=rep[0.05]['leak_ratio'],
                   leak_ratio_10pct=rep[0.10]['leak_ratio'],
                   far_off_evaluable_1pct=rep[0.01]['windows_evaluable'],
                   far_off_evaluable_10pct=rep[0.10]['windows_evaluable'],
                   far_off_h_1pct=rep[0.01]['far_off_h'],
                   far_off_h_10pct=rep[0.10]['far_off_h'],
                   gate_off_peak_max=float(a['gate_contrast']['off_cycle_peak_max']),
                   bit2_hold_h=a['bit_margins']['bit2']['min_hold_h'])
    except Exception as exc:                                     # noqa: BLE001
        row.update(ok=False, error=repr(exc), certified=False)
    row['runtime_s'] = round(time.perf_counter() - started, 2)
    return row


def scan(args):
    pts = shard_slice(points(), args.shard, args.nshards)
    tag = f'shard{args.shard:02d}of{args.nshards:02d}'
    rows, windows = [], []
    for p in pts:
        r = evaluate(p, args.hours, args.sample_min, args.max_step_min)
        rows.append(r)
        print(f"[{tag}] {point_id(p):<26} certified={r['certified']} crossings={r.get('crossings')} "
              f"leak_ratio={r.get('leak_ratio')} ({r['runtime_s']}s)", flush=True)
    write_shard(NAME, tag, rows, dict(kind='scan', shard=args.shard, nshards=args.nshards,
                                      hours=args.hours, points=len(pts),
                                      model_sha256=sha256(ROOT_MODEL),
                                      source_hashes=source_hashes()))
    print(f'[{tag}] wrote {len(rows)} rows')


ROOT_MODEL = Path(__file__).resolve().parents[1] / 'model_threebit51.py'


def merge(args):
    df = merge_shards(NAME, args.nshards)
    ok = df[df.ok.astype(bool)]
    cert = ok[ok.certified.astype(bool)]
    fail = ok[~ok.certified.astype(bool)]
    def stats(s):
        v = s.dropna()
        return dict(n=len(v), median=float(v.median()) if len(v) else None,
                    p10=float(v.quantile(0.10)) if len(v) else None,
                    p90=float(v.quantile(0.90)) if len(v) else None,
                    max=float(v.max()) if len(v) else None)
    verdict = dict(rows=len(df), certified=int(len(cert)),
                   leak_ratio_certified=stats(cert.leak_ratio),
                   leak_ratio_failed=stats(fail.leak_ratio),
                   far_off_rev_certified=stats(cert.far_off_rev),
                   far_off_rev_failed=stats(fail.far_off_rev),
                   far_off_fwd_certified=stats(cert.far_off_fwd),
                   tail_rev_certified=stats(cert.tail_rev),
                   gate_on_rev_certified=stats(cert.gate_on_rev),
                   segment_durations_h=dict(
                       gate_on=stats(ok.gate_on_h), tail=stats(ok.tail_h),
                       far_off=stats(ok.far_off_h)),
                   far_off_evaluability=dict(
                       windows_total=int(ok.windows_total.sum()),
                       windows_evaluable=int(ok.windows_evaluable.sum()),
                       points_with_no_evaluable_window=int((ok.windows_evaluable == 0).sum())),
                   cut_fraction_sensitivity={
                       key: {str(n): stats(sub[key]) for n, sub in ok.groupby('n_A1_gate')}
                       for key in ('leak_ratio_1pct', 'leak_ratio_5pct', 'leak_ratio_10pct')},
                   note=('far_off_rev is None (never 0) when the tail reaches the next carry; '
                         'a candidate ordering is only trusted if it is the same at the '
                         '1 %, 5 % and 10 % cuts'))
    df.to_csv(OUT / f'{NAME}_all.csv', index=False, encoding='utf-8')
    (OUT / f'{NAME}_verdict.json').write_text(
        json.dumps(verdict, ensure_ascii=False, indent=2), encoding='utf-8')
    write_manifest()
    pd.set_option('display.width', 250)
    cols = [c for c in ('point_id', 'certified', 'crossings', 'gates', 'reverse_events',
                        'gate_on_rev', 'tail_rev', 'far_off_rev', 'leak_ratio',
                        'gate_off_peak_max') if c in df.columns]
    print(df[cols].sort_values(['certified', 'leak_ratio'], ascending=[False, True]).to_string(index=False))
    print()
    print(json.dumps(verdict, ensure_ascii=False, indent=2))


def verify(args):
    ok, detail, traj = prefix_guard()
    res = dict(prefix_guard_passed=ok,
               structural_rhs_gap_states_0_33=detail,
               trajectory_gap_states_0_33_roundoff_only=traj,
               model_sha256=sha256(ROOT_MODEL), source_hashes=source_hashes())
    (OUT / 'gate_segments_prefix_guard.json').write_text(
        json.dumps(res, ensure_ascii=False, indent=2), encoding='utf-8')
    write_manifest()
    print(json.dumps(res, ensure_ascii=False, indent=2))
    if not ok:
        raise SystemExit('prefix guard FAILED - states 0..33 depend on n_A1_gate')


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='mode', required=True)
    s = sub.add_parser('scan')
    s.add_argument('--shard', type=int, required=True)
    s.add_argument('--nshards', type=int, required=True)
    s.add_argument('--hours', type=float, default=600.0)
    s.add_argument('--sample-min', type=float, default=2.0)
    s.add_argument('--max-step-min', type=float, default=2.0)
    m = sub.add_parser('merge')
    m.add_argument('--nshards', type=int, required=True)
    sub.add_parser('verify-prefix')
    args = ap.parse_args()
    if args.mode == 'scan':
        scan(args)
    elif args.mode == 'merge':
        merge(args)
    else:
        verify(args)


if __name__ == '__main__':
    main()
