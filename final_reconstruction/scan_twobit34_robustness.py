"""Sharded 75-point robustness scan for the formal 34-state two-bit model.

One script owns both calculation and aggregation so the acceptance criterion
cannot drift between the shard runner and the summary.

Examples
--------
python scan_twobit34_robustness.py scan --shard 0 --nshards 15
python scan_twobit34_robustness.py merge --nshards 15
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import platform
import sys
import time
import traceback
from dataclasses import asdict, replace
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy

from model import ROOT
from model_twobit34 import (CarryExpressionParameters, TwoBit34Model,
                            nominal_extension)
from verify_twobit_causal import analyse_solution


UM_VALUES = (5.5, 5.75, 6.0, 6.25, 6.5)
CARRY_MAT_VALUES = (25.0, 27.5, 30.0, 32.5, 35.0)
CARRY_MRNA_VALUES = (1.0, 2.0, 4.0)
GRID_KEYS = ('uM_per_au', 'carry_maturation_half_life_min',
             'carry_mrna_half_life_min')
OUT = ROOT / 'twobit34_results' / 'robustness'


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest().upper()


def grid_points():
    return [dict(zip(GRID_KEYS, values)) for values in
            itertools.product(UM_VALUES, CARRY_MAT_VALUES, CARRY_MRNA_VALUES)]


def _failure_reason(analysis):
    steady = analysis['steady_state']
    causal = analysis['causal_verdict']
    if not steady['valid_windows']:
        return 'readout_unlabelled_window'
    if not steady['increments_mod4']:
        return 'readout_wrong_increment'
    if not causal['exactly_one_gate_and_flip_per_late_reverse']:
        return 'carry_not_one_to_one'
    if not causal['causal_order_after_reverse_start']:
        return 'causal_order_failed'
    if not causal['bit1_directions_alternate']:
        return 'bit1_direction_failed'
    return ''


def evaluate(point, hours=300.0):
    started = time.perf_counter()
    row = dict(point)
    try:
        extension = replace(nominal_extension(), uM_per_au=point['uM_per_au'])
        carry = CarryExpressionParameters(
            mrna_half_life_min=point['carry_mrna_half_life_min'],
            activator_maturation_half_life_min=point['carry_maturation_half_life_min'],
            repressor_maturation_half_life_min=point['carry_maturation_half_life_min'])
        model = TwoBit34Model(extension, carry)
        sol = model.simulate(hours=hours, sample_min=2.0, max_step_min=2.0)
        analysis = analyse_solution(model.base, sol, asdict(extension), hours)
        c = analysis['causal_verdict']; cold = analysis['cold_start']; steady = analysis['steady_state']
        setup = analysis['margins']['bit1']['min_setup_h']
        hold = analysis['margins']['bit1']['min_hold_h']
        row.update(
            solver_success=True, error='', finite=bool(np.isfinite(sol.y).all()),
            minimum_state=float(sol.y.min()), maximum_state=float(sol.y.max()),
            samples=int(sol.t.size), reads=int(cold['reads']),
            cold_passed=bool(cold['passed']), cold_sequence=cold['sequence'],
            steady_passed=bool(steady['passed']), steady_sequence=steady['sequence'],
            minimum_commitment=float(steady['minimum_commitment']),
            one_to_one=bool(c['exactly_one_gate_and_flip_per_late_reverse']),
            causal_order=bool(c['causal_order_after_reverse_start']),
            directions_alternate=bool(c['bit1_directions_alternate']),
            unassigned_gate_events=len(c['unassigned_gate_events']),
            unassigned_bit1_crossings=len(c['unassigned_bit1_crossings']),
            reverse_events=len(analysis['reverse_events']),
            gate_events=len(analysis['gate_events']),
            certified=bool(analysis['certified']), failure_reason=_failure_reason(analysis),
            bit0_min=float(analysis['signal_ranges']['S0'][0]),
            bit0_max=float(analysis['signal_ranges']['S0'][1]),
            bit1_min=float(analysis['signal_ranges']['S1'][0]),
            bit1_max=float(analysis['signal_ranges']['S1'][1]),
            Jrev0_peak=float(analysis['signal_ranges']['Jrev0'][1]),
            g0_peak=float(analysis['signal_ranges']['g0'][1]),
            Int1_source_peak=float(analysis['signal_ranges']['Int1_source_peak']),
            bit1_setup_h=float(setup) if setup is not None else np.nan,
            bit1_hold_h=float(hold) if hold is not None else np.nan,
            timing_margin_h=float(min(setup, hold)) if setup is not None and hold is not None else np.nan)
    except Exception as exc:  # noqa: BLE001
        row.update(solver_success=False, error=repr(exc), finite=False,
                   cold_passed=False, steady_passed=False, certified=False,
                   failure_reason='integration_failed', error_trace=traceback.format_exc(limit=3))
    row['runtime_s'] = round(time.perf_counter() - started, 3)
    return row


def scan(args):
    OUT.mkdir(parents=True, exist_ok=True)
    points = grid_points()
    mine = points[args.shard::args.nshards]
    tag = f'shard{args.shard:02d}of{args.nshards:02d}'
    csv_path = OUT / f'robustness_{tag}.csv'
    rows = []
    for i, point in enumerate(mine):
        row = evaluate(point, args.hours)
        rows.append(row)
        fields = sorted({k for r in rows for k in r})
        with csv_path.open('w', newline='', encoding='utf-8') as fh:
            writer = csv.DictWriter(fh, fieldnames=fields, extrasaction='ignore')
            writer.writeheader(); writer.writerows(rows)
        print(f"[{tag}] {i + 1}/{len(mine)} uM={point['uM_per_au']:.2f} "
              f"mat={point['carry_maturation_half_life_min']:.1f} "
              f"mRNA={point['carry_mrna_half_life_min']:.1f} "
              f"certified={row['certified']} reason={row['failure_reason'] or 'pass'} "
              f"margin={row.get('timing_margin_h', float('nan')):.3f}", flush=True)

    meta = dict(shard=args.shard, nshards=args.nshards, hours=args.hours,
                points=len(mine), expected_total=len(points), grid=dict(
                    uM_per_au=UM_VALUES,
                    carry_maturation_half_life_min=CARRY_MAT_VALUES,
                    carry_mrna_half_life_min=CARRY_MRNA_VALUES),
                criterion='verify_twobit_causal.analyse_solution',
                python=sys.version, platform=platform.platform(),
                versions=dict(numpy=np.__version__, pandas=pd.__version__, scipy=scipy.__version__),
                sha256=dict(model_py=sha256(ROOT / 'model.py'),
                            model_twobit34_py=sha256(ROOT / 'model_twobit34.py'),
                            verifier_py=sha256(ROOT / 'verify_twobit_causal.py'),
                            scanner_py=sha256(Path(__file__)), csv=sha256(csv_path)))
    (OUT / f'meta_{tag}.json').write_text(
        json.dumps(meta, ensure_ascii=False, indent=2), encoding='utf-8')
    print(f'[{tag}] wrote {csv_path}')


def _as_bool(series):
    return series.map(lambda x: str(x).strip().lower() in ('true', '1'))


def _heatmaps(data):
    fig, axes = plt.subplots(2, 3, figsize=(14, 8), sharex=True, sharey=True,
                             constrained_layout=True)
    for col, mrna in enumerate(CARRY_MRNA_VALUES):
        sub = data[data.carry_mrna_half_life_min == mrna]
        passed = sub.pivot(index='carry_maturation_half_life_min',
                           columns='uM_per_au', values='certified').astype(float)
        margin = sub.pivot(index='carry_maturation_half_life_min',
                           columns='uM_per_au', values='timing_margin_h').astype(float)
        im0 = axes[0, col].imshow(passed.values, origin='lower', aspect='auto',
                                  vmin=0, vmax=1, cmap='RdYlGn')
        finite = margin.values[np.isfinite(margin.values)]
        vmax = max(float(finite.max()) if finite.size else 1.0, .1)
        im1 = axes[1, col].imshow(margin.values, origin='lower', aspect='auto',
                                  vmin=0, vmax=vmax, cmap='viridis')
        for r in range(passed.shape[0]):
            for c in range(passed.shape[1]):
                axes[0, col].text(c, r, 'PASS' if passed.iloc[r, c] else 'FAIL',
                                  ha='center', va='center', fontsize=7,
                                  color='black' if passed.iloc[r, c] else 'white')
                value = margin.iloc[r, c]
                axes[1, col].text(c, r, '—' if not np.isfinite(value) else f'{value:.2f}',
                                  ha='center', va='center', fontsize=8,
                                  color='white' if np.isfinite(value) and value < .55*vmax else 'black')
        for ax in axes[:, col]:
            ax.set_xticks(range(len(passed.columns)), [f'{x:g}' for x in passed.columns])
            ax.set_yticks(range(len(passed.index)), [f'{x:g}' for x in passed.index])
            ax.set_xlabel('uM per a.u.')
        axes[0, col].set_title(f'carry mRNA half-life = {mrna:g} min')
    axes[0, 0].set_ylabel('A0/F0 maturation half-life (min)')
    axes[1, 0].set_ylabel('A0/F0 maturation half-life (min)')
    fig.suptitle('Formal 34-state robustness: certification and timing margin', fontsize=15)
    fig.colorbar(im0, ax=axes[0, :], shrink=.7, label='certified')
    fig.colorbar(im1, ax=axes[1, :], shrink=.7, label='min(bit1 setup, hold), h')
    fig.savefig(OUT / 'robustness_heatmaps.png', dpi=220, bbox_inches='tight', facecolor='white')
    fig.savefig(OUT / 'robustness_heatmaps.pdf', bbox_inches='tight', facecolor='white')
    plt.close(fig)


def merge(args):
    frames, metas = [], []
    for shard in range(args.nshards):
        tag = f'shard{shard:02d}of{args.nshards:02d}'
        csv_path, meta_path = OUT / f'robustness_{tag}.csv', OUT / f'meta_{tag}.json'
        if not csv_path.exists() or not meta_path.exists():
            raise SystemExit(f'missing shard artefact for {tag}')
        frames.append(pd.read_csv(csv_path)); metas.append(json.loads(meta_path.read_text(encoding='utf-8')))
    data = pd.concat(frames, ignore_index=True)
    for key in GRID_KEYS:
        data[key] = data[key].astype(float)
    for key in ('solver_success', 'finite', 'cold_passed', 'steady_passed', 'one_to_one',
                'causal_order', 'directions_alternate', 'certified'):
        data[key] = _as_bool(data[key])
    data = data.sort_values(list(GRID_KEYS)).reset_index(drop=True)

    expected = {tuple(p[k] for k in GRID_KEYS) for p in grid_points()}
    got = {tuple(float(getattr(row, k)) for k in GRID_KEYS) for row in data.itertuples()}
    completeness = dict(rows=len(data), unique_points=len(got), expected_points=len(expected),
                        duplicates=len(data) - len(got), missing=sorted(expected - got),
                        unexpected=sorted(got - expected))
    if completeness['rows'] != 75 or completeness['duplicates'] or completeness['missing'] or completeness['unexpected']:
        raise SystemExit(f'incomplete scan: {completeness}')

    all_path = OUT / 'robustness_all.csv'
    data.to_csv(all_path, index=False, encoding='utf-8')
    passed = data[data.certified].copy()
    best = passed.sort_values(['timing_margin_h', 'minimum_commitment', 'g0_peak'],
                              ascending=False).head(10)
    columns = list(GRID_KEYS) + ['timing_margin_h', 'bit1_setup_h', 'bit1_hold_h',
                                'minimum_commitment', 'g0_peak', 'Int1_source_peak',
                                'cold_sequence']
    best_rows = [{k: (bool(v) if isinstance(v, (np.bool_, bool)) else
                      float(v) if isinstance(v, (np.floating, float)) else v)
                  for k, v in row.items()} for row in best[columns].to_dict(orient='records')]
    hash_keys = ('model_py', 'model_twobit34_py', 'verifier_py', 'scanner_py')
    environment_uniform = all(len({m['sha256'][key] for m in metas}) == 1 for key in hash_keys)
    summary = dict(
        completeness=completeness, solver_failures=int((~data.solver_success).sum()),
        nonfinite_runs=int((~data.finite).sum()),
        cold_passed=int(data.cold_passed.sum()), steady_passed=int(data.steady_passed.sum()),
        certified=int(data.certified.sum()), rate=float(data.certified.mean()),
        failure_reasons={str(k): int(v) for k, v in data.failure_reason.fillna('').value_counts().items()},
        by_carry_mrna={str(v): dict(points=int((data.carry_mrna_half_life_min == v).sum()),
                                    certified=int(data.loc[data.carry_mrna_half_life_min == v, 'certified'].sum()),
                                    rate=float(data.loc[data.carry_mrna_half_life_min == v, 'certified'].mean()))
                       for v in CARRY_MRNA_VALUES},
        best_points=best_rows,
        environment=dict(shards=args.nshards, uniform=environment_uniform,
                         python=sorted({m['python'] for m in metas}),
                         platform=sorted({m['platform'] for m in metas}),
                         versions={k: sorted({m['versions'][k] for m in metas})
                                   for k in ('numpy', 'pandas', 'scipy')}))
    (OUT / 'robustness_summary.json').write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
    _heatmaps(data)

    lines = ['# 正式34状态模型联合鲁棒性扫描', '',
             f"- 完整性：{completeness['rows']}/{completeness['expected_points']}，重复{completeness['duplicates']}，缺失{len(completeness['missing'])}",
             f"- 积分失败：{summary['solver_failures']}；非有限轨迹：{summary['nonfinite_runs']}",
             f"- 冷启动通过：{summary['cold_passed']}/75；稳态通过：{summary['steady_passed']}/75",
             f"- 完整因果认证：**{summary['certified']}/75 = {summary['rate']:.1%}**", '',
             '## carry mRNA半衰期分层', '',
             '| mRNA半衰期/min | 通过/总数 | 通过率 |', '|---:|---:|---:|']
    for value, row in summary['by_carry_mrna'].items():
        lines.append(f"| {value} | {row['certified']}/{row['points']} | {row['rate']:.1%} |")
    lines += ['', '## 最佳点（按最小setup/hold裕量排序）', '',
              '| uM/a.u. | A/F成熟/min | carry mRNA/min | 最小裕量/h | setup/h | hold/h |',
              '|---:|---:|---:|---:|---:|---:|']
    for row in best_rows:
        lines.append(f"| {row['uM_per_au']:g} | {row['carry_maturation_half_life_min']:g} | "
                     f"{row['carry_mrna_half_life_min']:g} | {row['timing_margin_h']:.3f} | "
                     f"{row['bit1_setup_h']:.3f} | {row['bit1_hold_h']:.3f} |")
    lines += ['', '![联合鲁棒性热图](robustness_heatmaps.png)', '',
              '判据与正式模型保持不变；本扫描只改变三个新增、未标定的接口参数，不修改ZENG。', '']
    (OUT / 'robustness_summary.md').write_text('\n'.join(lines), encoding='utf-8')

    manifest = {}
    for path in sorted(OUT.glob('*')):
        if path.is_file() and not path.name.startswith('SHA256SUMS'):
            manifest[str(path.relative_to(ROOT))] = sha256(path)
    for name in ('model.py', 'model_twobit34.py', 'verify_twobit_causal.py',
                 'scan_twobit34_robustness.py'):
        path = ROOT / name; manifest[str(path.relative_to(ROOT))] = sha256(path)
    (OUT / 'SHA256SUMS.json').write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True), encoding='utf-8')
    (OUT / 'SHA256SUMS.txt').write_text(
        ''.join(f'{value}  {name}\n' for name, value in sorted(manifest.items())), encoding='utf-8')
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print(f'wrote {all_path}')


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='mode', required=True)
    scan_ap = sub.add_parser('scan')
    scan_ap.add_argument('--shard', type=int, required=True)
    scan_ap.add_argument('--nshards', type=int, default=15)
    scan_ap.add_argument('--hours', type=float, default=300.0)
    merge_ap = sub.add_parser('merge')
    merge_ap.add_argument('--nshards', type=int, default=15)
    args = ap.parse_args()
    scan(args) if args.mode == 'scan' else merge(args)


if __name__ == '__main__':
    main()
