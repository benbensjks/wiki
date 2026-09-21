"""Stage-2 probe: does the bit2 crossing count keep rising past 60 min maturation?

Seven points, 600 h each, sharded over seven single-point shards:

    carry1 mRNA   = 3, 4, 6 min   x   A1/F1 maturation = 70, 85 min   (6 points)
    control       = the stage-1 best point, 4 min / 60 min            (1 point)

Only the new A1/F1 expression parameters move.  The scanner itself is imported
from scan_threebit51_carry1_dsh (unchanged), so the criterion, the model wrapper
and every recorded field are identical to stage 1.

Decision rule fixed before running:
    crossings keep rising and a point certifies  -> fill the maturation grid
    crossings saturate (~11-12) regardless       -> the R2-floor binds, stop tuning
    gate contrast degrades further               -> the fix must sharpen, not widen

    python probe_carry1_stage2.py scan  --shard 0 --nshards 7 --hours 600
    python probe_carry1_stage2.py merge --nshards 7 --hours 600
"""
from __future__ import annotations

import argparse
import csv
import itertools
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from model import ROOT
from scan_threebit51_carry1_dsh import evaluate, sha256

OUT = ROOT / 'threebit51_results' / 'carry1_probe2'
MRNA = (3.0, 4.0, 6.0)
MAT = (70.0, 85.0)
CONTROL = dict(carry1_mrna_min=4.0, carry1_maturation_min=60.0)
STAGE1_BEST = dict(mrna=4.0, mat=60.0, crossings=12, gate_width_h=1.333,
                   off_on=0.0163, s2_cycle_min_median=0.629)
FROZEN_FILES = ('model.py', 'model_twobit34.py', 'model_threebit51.py',
                'verify_threebit51.py', 'verify_twobit_causal.py',
                'scan_threebit51_carry1_dsh.py')


def jobs():
    out = [dict(carry1_mrna_min=a, carry1_maturation_min=b)
           for a, b in itertools.product(MRNA, MAT)]
    out.append(dict(CONTROL))
    return out


def scan(args):
    OUT.mkdir(parents=True, exist_ok=True)
    allj = jobs()
    mine = allj[args.shard::args.nshards]
    tag = f'shard{args.shard:02d}of{args.nshards:02d}'
    csv_path = OUT / f'probe_{tag}.csv'
    rows = []
    for i, job in enumerate(mine):
        r = evaluate(job, args.hours)
        r.pop('cycle_metrics', None)
        rows.append(r)
        fields = sorted({k for x in rows for k in x})
        with csv_path.open('w', newline='', encoding='utf-8') as fh:
            w = csv.DictWriter(fh, fieldnames=fields, extrasaction='ignore')
            w.writeheader()
            w.writerows(rows)
        print(f"[{tag}] {i + 1}/{len(mine)} mRNA={job['carry1_mrna_min']:<4g} "
              f"mat={job['carry1_maturation_min']:<5g} cert={r['certified']} "
              f"cross={r.get('bit2_crossings')} err={r.get('crossing_error')} "
              f"contrast={r.get('gate_contrast_passed')} "
              f"gateW={r.get('gate_duration_median_h'):.2f} "
              f"S2min_med={r.get('s2_cycle_min_median'):.3f}", flush=True)
    meta = dict(shard=args.shard, nshards=args.nshards, hours=args.hours,
                points=len(mine), expected_total=len(allj),
                grid=dict(carry1_mrna_min=MRNA, carry1_maturation_min=MAT,
                          control=CONTROL),
                criterion='verify_threebit51.analyse_threebit via '
                          'scan_threebit51_carry1_dsh.evaluate (unchanged)',
                sha256={f.replace('.py', ''): sha256(ROOT / f) for f in FROZEN_FILES},
                csv=csv_path.name, csv_sha256=sha256(csv_path))
    (OUT / f'meta_{tag}.json').write_text(json.dumps(meta, ensure_ascii=False, indent=2),
                                          encoding='utf-8')
    print(f'[{tag}] wrote {csv_path}', flush=True)


def _bool(s):
    return s.map(lambda x: str(x).strip().lower() in ('true', '1'))


def merge(args):
    frames = []
    for i in range(args.nshards):
        p = OUT / f'probe_shard{i:02d}of{args.nshards:02d}.csv'
        if not p.exists():
            raise SystemExit(f'missing shard artefact {p}')
        frames.append(pd.read_csv(p))
    d = pd.concat(frames, ignore_index=True)
    for c in ('solver_success', 'finite', 'cold_passed', 'steady_passed', 'certified',
              'one_to_one', 'causal_order', 'directions_alternate', 'gate_contrast_passed'):
        if c in d.columns:
            d[c] = _bool(d[c])
    d = d.sort_values(['carry1_mrna_min', 'carry1_maturation_min']).reset_index(drop=True)

    expected = {(j['carry1_mrna_min'], j['carry1_maturation_min']) for j in jobs()}
    got = {(float(r.carry1_mrna_min), float(r.carry1_maturation_min)) for r in d.itertuples()}
    completeness = dict(rows=len(d), unique=len(got), expected=len(expected),
                        duplicates=len(d) - len(got), missing=sorted(expected - got),
                        unexpected=sorted(got - expected))
    if completeness['rows'] != len(expected) or completeness['duplicates'] \
            or completeness['missing'] or completeness['unexpected']:
        raise SystemExit(f'incomplete probe: {completeness}')
    d.to_csv(OUT / 'probe_all.csv', index=False, encoding='utf-8')

    # compare the control against its stage-1 value
    ctl = d[(d.carry1_mrna_min == CONTROL['carry1_mrna_min'])
            & (d.carry1_maturation_min == CONTROL['carry1_maturation_min'])]
    control_check = None
    if len(ctl):
        row = ctl.iloc[0]
        control_check = dict(
            crossings=float(row.bit2_crossings),
            stage1_crossings=STAGE1_BEST['crossings'],
            reproduces=bool(float(row.bit2_crossings) == STAGE1_BEST['crossings']),
            gate_width_h=float(row.gate_duration_median_h),
            stage1_gate_width_h=STAGE1_BEST['gate_width_h'],
            off_on=float(row.off_on_peak_ratio),
            s2_cycle_min_median=float(row.s2_cycle_min_median))

    m70 = d[d.carry1_maturation_min == 70.0]
    m85 = d[d.carry1_maturation_min == 85.0]
    lift = dict(
        crossings_70=[float(x) for x in m70.bit2_crossings],
        crossings_85=[float(x) for x in m85.bit2_crossings],
        mrna_70=[float(x) for x in m70.carry1_mrna_min],
        mrna_85=[float(x) for x in m85.carry1_mrna_min],
    )
    stage1_at_60 = {3.0: None, 4.0: 12.0, 6.0: None}   # 3 and 6 were not in stage 1
    trend = dict(
        best_crossings=float(d.bit2_crossings.max()),
        best_rows=d.sort_values('bit2_crossings', ascending=False)
        .head(3)[['carry1_mrna_min', 'carry1_maturation_min', 'bit2_crossings',
                  'crossing_error', 'certified', 'gate_contrast_passed',
                  'gate_duration_median_h', 'off_on_peak_ratio',
                  's2_cycle_min_median', 'steady_sequence']].to_dict(orient='records'),
        certified=int(d.certified.sum()),
        gate_contrast_passed=int(d.gate_contrast_passed.sum()),
        stage1_reference=STAGE1_BEST)

    if trend['certified']:
        verdict = 'trend confirmed: at least one point certifies -> fill the maturation grid'
    elif trend['best_crossings'] >= 14:
        verdict = 'crossings reach 14 but certification still fails -> inspect causal/contrast'
    elif trend['best_crossings'] > STAGE1_BEST['crossings']:
        verdict = 'trend continues upward but has not reached 14 -> extend maturation further'
    else:
        verdict = ('saturation: crossings do not exceed the stage-1 best of %d -> the R2 '
                   'floor binds, stop tuning maturation' % STAGE1_BEST['crossings'])

    summary = dict(completeness=completeness, control_check=control_check,
                   lift=lift, stage1_at_60=stage1_at_60, trend=trend, verdict=verdict,
                   stage2_mrna_vs_stage1={3.0: 'new', 4.0: 'same as stage-1 best', 6.0: 'new'})
    (OUT / 'probe_summary.json').write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, default=str), encoding='utf-8')

    cols = [c for c in ('carry1_mrna_min', 'carry1_maturation_min', 'certified', 'crossing_error',
                        'bit2_crossings', 'bit1_reverse_events', 'gate_events',
                        'gate_contrast_passed', 'off_on_peak_ratio', 'gate_duration_median_h',
                        's2_cycle_min_median', 'abs_net_ds2_median', 'rev_over_fwd_median',
                        'bit2_low_cycles', 'steady_sequence') if c in d.columns]
    print(json.dumps({k: summary[k] for k in ('completeness', 'control_check', 'verdict')},
                     ensure_ascii=False, indent=2, default=str))
    print(d[cols].to_string(index=False))
    print(f'wrote {OUT / "probe_all.csv"}')


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='mode', required=True)
    s = sub.add_parser('scan')
    s.add_argument('--shard', type=int, required=True)
    s.add_argument('--nshards', type=int, default=7)
    s.add_argument('--hours', type=float, default=600.0)
    m = sub.add_parser('merge')
    m.add_argument('--nshards', type=int, default=7)
    m.add_argument('--hours', type=float, default=600.0)
    a = ap.parse_args()
    scan(a) if a.mode == 'scan' else merge(a)


if __name__ == '__main__':
    main()
