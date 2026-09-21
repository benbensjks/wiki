"""Merge the 16 shard CSVs into one 800-point table and summarise.

Reports cold-start and steady-state pass rates separately, cross-tabulates
against the archived rule on the same trajectories and against the archived
100 h scan, and emits a SHA256 manifest for every artefact.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from model import ROOT
from rescan800_readwindow import GRID_KEYS, grid_configs, sha256

SCAN_DIR = ROOT / 'bit0_results' / 'rescan800'
KEYS = ['uM_per_au', 'maturation_half_life_min', 'complex_on_au_inv_h', 'complex_off_h', 'add_growth']
ARCHIVE = ROOT / 'bit0_results' / 'extension_parameter_scan.csv'


def verdict_from_labels(labels, drop, minimum_reads=8):
    chosen = labels[drop:]
    valid = len(chosen) >= minimum_reads and all(l != 'x' for l in chosen)
    increments = bool(valid and all((int(chosen[i + 1]) - int(chosen[i])) % 2 == 1
                                    for i in range(len(chosen) - 1)))
    return valid and increments


def failure_reason(labels, drop, solver_ok, minimum_reads=8):
    """Why a window set fails: one label per point, in priority order."""
    if not solver_ok:
        return 'integration_failed'
    chosen = labels[drop:]
    if len(chosen) < minimum_reads:
        return 'too_few_reads'
    if all(l == '1' for l in chosen):
        return 'stuck_high'
    if all(l == '0' for l in chosen):
        return 'stuck_low'
    if any(l == 'x' for l in chosen):
        return 'unlabelled_window_x'
    if not all((int(chosen[i + 1]) - int(chosen[i])) % 2 == 1 for i in range(len(chosen) - 1)):
        return 'no_alternation'
    return ''


def audit_versions(nshards):
    metas = []
    for s in range(nshards):
        p = SCAN_DIR / f'meta_shard{s:02d}of{nshards:02d}.json'
        metas.append(json.loads(p.read_text(encoding='utf-8')))
    python_strings = sorted({m['python'] for m in metas})
    platforms = sorted({m['platform'] for m in metas})
    hashes = {}
    for key in ('model_py', 'bit0_diagnostic_py', 'verify_twobit_causal_py', 'scanner_py', 'upstream_py'):
        vals = sorted({m['sha256'][key] for m in metas})
        hashes[key] = dict(unique_values=len(vals), value=vals[0] if len(vals) == 1 else None,
                           all_values=vals if len(vals) > 1 else None)
    return dict(shards=len(metas), distinct_python_strings=python_strings,
                distinct_platforms=platforms, shard_hashes_uniform=all(
                    v['unique_values'] == 1 for v in hashes.values()), hashes=hashes,
                hours=sorted({m['hours'] for m in metas}),
                sample_min=sorted({m['sample_min'] for m in metas}),
                max_step_min=sorted({m['max_step_min'] for m in metas}),
                grid_identical_across_shards=len({json.dumps(m['grid'], sort_keys=True) for m in metas}) == 1)


def load_all(nshards=16):
    frames = []
    for s in range(nshards):
        p = SCAN_DIR / f'rescan800_shard{s:02d}of{nshards:02d}.csv'
        if not p.exists():
            raise SystemExit(f'missing shard file {p}')
        frames.append(pd.read_csv(p))
    d = pd.concat(frames, ignore_index=True)
    return d


def check_completeness(d):
    expected = {(round(c['uM_per_au'], 10), round(c['maturation_half_life_min'], 10),
                 round(c['complex_on_au_inv_h'], 10), round(c['complex_off_h'], 10),
                 bool(c['add_growth'])) for c in grid_configs()}
    got = {(round(r.uM_per_au, 10), round(r.maturation_half_life_min, 10),
            round(r.complex_on_au_inv_h, 10), round(r.complex_off_h, 10),
            bool(r.add_growth)) for r in d.itertuples()}
    return dict(rows=len(d), unique_points=len(got), expected_points=len(expected),
                duplicates=len(d) - len(got), missing=sorted(expected - got),
                unexpected=sorted(got - expected))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--nshards', type=int, default=16)
    args = ap.parse_args()

    d = load_all(args.nshards)
    completeness = check_completeness(d)
    d['cold_start_sequence'] = d['cold_start_sequence'].fillna('')
    d['steady_state_sequence'] = d['steady_state_sequence'].fillna('')

    # burn-in sweep straight off the full label string
    for drop in range(0, 13):
        d[f'pass_drop{drop}'] = [verdict_from_labels(list(s), drop) for s in d.cold_start_sequence]
    passmat = d[[f'pass_drop{k}' for k in range(13)]].to_numpy(dtype=bool)
    d['pass_at_drop2'] = passmat[:, 2]
    d['pass_at_drop4'] = passmat[:, 4]
    d['pass_at_drop6'] = passmat[:, 6]
    d['first_pass_drop'] = np.where(passmat.any(axis=1), passmat.argmax(axis=1), -1)
    d['burn_in_only'] = d.first_pass_drop > 0
    solver_ok = d.solver.astype(bool).to_numpy()
    d['failure_reason_cold'] = [failure_reason(list(s), 0, bool(ok))
                                for s, ok in zip(d.cold_start_sequence, solver_ok)]
    d['failure_reason_steady'] = [failure_reason(list(s), 4, bool(ok))
                                  for s, ok in zip(d.cold_start_sequence, solver_ok)]

    arch = pd.read_csv(ARCHIVE)
    for frame in (arch, d):
        for k in KEYS:
            if k == 'add_growth':
                frame[k] = frame[k].map(lambda v: str(v).strip().lower() in ('true', '1'))
            else:
                frame[k] = frame[k].astype(float)
    merged = d.merge(arch[KEYS + ['stable_alternation', 'codes', 'dynamic_range',
                                  'alternation_accuracy', 'confident_fraction']]
                     .rename(columns={'stable_alternation': 'archived100h_stable',
                                      'codes': 'archived100h_codes',
                                      'dynamic_range': 'archived100h_dynamic_range',
                                      'alternation_accuracy': 'archived100h_alternation_accuracy',
                                      'confident_fraction': 'archived100h_confident_fraction'}),
                     on=KEYS, how='left', validate='one_to_one')
    merged = merged.sort_values(KEYS).reset_index(drop=True)

    out_csv = SCAN_DIR / 'rescan800_all.csv'
    merged.to_csv(out_csv, index=False, encoding='utf-8')

    solved = merged[merged.solver == True]                       # noqa: E712
    n = len(merged)
    cold = int(merged.cold_start_passed.sum())
    steady = int(merged.steady_state_passed.sum())
    meta0 = json.loads((SCAN_DIR / f'meta_shard00of{args.nshards:02d}.json').read_text(encoding='utf-8'))

    new = merged.steady_state_passed.to_numpy(dtype=bool)
    old300 = merged.archived_rule_stable.fillna(False).to_numpy(dtype=bool)
    old100 = merged.archived100h_stable.fillna(False).to_numpy(dtype=bool)

    def contingency(new_pass, old_pass):
        return dict(both_pass=int((new_pass & old_pass).sum()),
                    new_only=int((new_pass & ~old_pass).sum()),
                    old_only=int((~new_pass & old_pass).sum()),
                    both_fail=int((~new_pass & ~old_pass).sum()))

    numeric = merged.select_dtypes(include=[np.number])
    nonfinite, all_nan = {}, {}
    for c in numeric.columns:
        v = numeric[c].to_numpy(dtype=float)
        if np.all(np.isnan(v)):
            all_nan[c] = int(len(v))       # e.g. the empty `error` column of solved runs
            continue
        k = int((~np.isfinite(v)).sum())
        if k:
            nonfinite[c] = k

    summary = dict(
        hours=float(meta0['hours']),
        sample_min=float(meta0['sample_min']),
        max_step_min=float(meta0['max_step_min']),
        points=n, solver_failures=int((~merged.solver.astype(bool)).sum()),
        nonfinite_columns=nonfinite,
        all_nan_columns_excluded=all_nan,
        all_runs_finite=bool(not nonfinite),
        completeness=completeness,
        criterion=dict(name='finite read window (verify_twobit_causal geometry, mod-2 reduction)',
                       window=('read window centred on the flux trough; total width = '
                               'READ_WINDOW_FRACTION x clock period = 20% of the period, '
                               'i.e. trough +/- 10% (measured: 2.133 h on a 10.600 h cycle)'),
                       bands=[0.30, 0.70], min_occupancy=0.80, min_reads=8,
                       cold_drop=0, steady_drop=4,
                       source_module='verify_twobit_causal._read_windows (imported verbatim)',
                       reduction='mod 2 instead of mod 4 for a single bit',
                       burn_in_definition='drop the first k reads, then require >=8 labelled alternating reads'),
        pass_rates=dict(cold_start=dict(passed=cold, rate=cold / n),
                        steady_state=dict(passed=steady, rate=steady / n),
                        drop2=int(merged.pass_at_drop2.sum()),
                        drop4=int(merged.pass_at_drop4.sum()),
                        drop6=int(merged.pass_at_drop6.sum()),
                        burn_in_only=int(merged.burn_in_only.sum()),
                        never_passed=int((merged.first_pass_drop < 0).sum())),
        criterion_contingency=dict(
            new_steady_vs_archived_rule_same_300h_trajectory=contingency(new, old300),
            new_steady_vs_archived_100h_scan=contingency(new, old100)),
        failure_reasons=dict(
            steady_state={k: int(v) for k, v in
                          merged.failure_reason_steady.value_counts().items()},
            cold_start={k: int(v) for k, v in
                        merged.failure_reason_cold.value_counts().items()}),
        first_pass_drop_histogram={str(k): int(v) for k, v in
                                   merged.first_pass_drop.value_counts().sort_index().items()},
        environment_audit=audit_versions(args.nshards),
        burn_in_sweep={f'drop{k}': dict(passed=int(merged[f'pass_drop{k}'].sum()),
                                        rate=float(merged[f'pass_drop{k}'].mean()))
                       for k in range(13)},
        by_axis={k: {str(v): dict(points=int((merged[k] == v).sum()),
                                  cold_rate=float(merged.loc[merged[k] == v, 'cold_start_passed'].mean()),
                                  steady_rate=float(merged.loc[merged[k] == v, 'steady_state_passed'].mean()))
                     for v in sorted(merged[k].unique())} for k in KEYS},
        criterion_change=dict(
            archived_100h_stable=int(merged.archived100h_stable.fillna(False).sum()),
            archived_rule_on_300h=int(merged.archived_rule_stable.fillna(False).sum()),
            archived_100h_stable_but_steady_pass=int(
                ((merged.archived100h_stable == True) & merged.steady_state_passed).sum()),  # noqa: E712
            archived_100h_not_stable_but_steady_pass=int(
                ((merged.archived100h_stable != True) & merged.steady_state_passed).sum()),  # noqa: E712
            steady_pass_but_archived_rule_300h_fail=int(
                (merged.steady_state_passed & (merged.archived_rule_stable != True)).sum()),  # noqa: E712
            archived_rule_300h_pass_but_steady_fail=int(
                ((merged.archived_rule_stable == True) & ~merged.steady_state_passed).sum())),  # noqa: E712
        dwell_levels=dict(S0_dwell_low_median=float(merged.S0_dwell_low.median()),
                          S0_dwell_high_median=float(merged.S0_dwell_high.median())),
        interpretation_caveat=(
            'This grid only contains uM_per_au in {0.1, 0.3, 1, 3, 10} and maturation half-life in '
            '{1, 2.5, 5, 10, 20} min. The independently verified working band (uM_per_au ~ 6-10, '
            'maturation 12-32 min) is mostly NOT inside this grid, so the pass rate below measures '
            'the original grid, not the true width of the feasible window.'),
        passing_points=[dict(uM_per_au=float(r.uM_per_au),
                             maturation_half_life_min=float(r.maturation_half_life_min),
                             complex_on_au_inv_h=float(r.complex_on_au_inv_h),
                             complex_off_h=float(r.complex_off_h),
                             add_growth=bool(r.add_growth),
                             first_pass_drop=int(r.first_pass_drop),
                             archived_rule_stable_300h=bool(r.archived_rule_stable),
                             min_commitment=float(r.steady_state_min_commitment))
                        for r in merged[merged.steady_state_passed].itertuples()],
    )
    out_json = SCAN_DIR / 'rescan800_summary.json'
    out_json.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')

    # short markdown companion
    bi = summary['burn_in_sweep']
    cc = summary['criterion_contingency']
    c1 = cc['new_steady_vs_archived_rule_same_300h_trajectory']
    c2 = cc['new_steady_vs_archived_100h_scan']
    ea = summary['environment_audit']
    md = [f"# 原 800 点重扫（有限读窗判据，{int(summary['hours'])} h）", '',
          '## 判据定义（逐字复用验证器，只做两位→一位的 mod-2 归约）', '',
          f"- 读窗：以通量谷为中心，**总宽 = READ_WINDOW_FRACTION × 时钟周期 = 周期的 20%**，"
          f"即谷值前后各约 10%（实测 2.133 h / 10.600 h = 0.202）。",
          f"- 标签：窗口内 ≥80% 采样 ≥0.70 记 1，≤0.30 记 0，否则记 `x`。",
          f"- 冷启动：drop 0；稳态：丢弃前 4 个读窗。两者都要求剩余读窗 ≥8、全部有标签、标签严格交替（mod 2）。",
          f"- 来源：`{summary['criterion']['source_module']}`；参数 `sample_min={summary['sample_min']}`、"
          f"`max_step_min={summary['max_step_min']}`（与归档 800 点扫描相同）。", '',
          '## 数据完整性', '',
          f"- 总行数 **{n}**；唯一参数组合 **{completeness['unique_points']}**（期望 "
          f"{completeness['expected_points']}）；重复 {completeness['duplicates']}；"
          f"缺失 {len(completeness['missing'])}；多余 {len(completeness['unexpected'])}。",
          f"- 积分失败：**{summary['solver_failures']}**。",
          f"- 非有限值列：{summary['nonfinite_columns'] if summary['nonfinite_columns'] else '无'}"
          f"（全部 {n} 个点的状态极值/驻留/commitment 均有限：{summary['all_runs_finite']}）；"
          f"全空列已排除：{list(summary['all_nan_columns_excluded'])}。",
          f"- 16 个分片环境一致：python {ea['distinct_python_strings'][0].splitlines()[0]}；"
          f"模型/脚本哈希跨分片一致 = {ea['shard_hashes_uniform']}；网格跨分片一致 = {ea['grid_identical_across_shards']}。", '',
          '## 通过率', '',
          f"- **冷启动通过 {cold}/{n} = {cold / n:.1%}**",
          f"- **稳态通过 {steady}/{n} = {steady / n:.1%}**",
          f"- 丢弃 2 / 4 / 6 个读窗后通过：{summary['pass_rates']['drop2']} / "
          f"{summary['pass_rates']['drop4']} / {summary['pass_rates']['drop6']}",
          f"- 只在更长烧入后通过（burn-in only）：{summary['pass_rates']['burn_in_only']}",
          f"- 任何烧入长度都不通过：{summary['pass_rates']['never_passed']}", '',
          '## 烧入期敏感性（丢弃前 k 个读窗后的通过数）', '',
          '| drop | ' + ' | '.join(str(k) for k in range(13)) + ' |',
          '|---|' + '---|' * 13,
          '| passed | ' + ' | '.join(str(bi[f'drop{k}']['passed']) for k in range(13)) + ' |', '',
          '## 新旧判据 2×2 列联表', '',
          '| 对照 | 都通过 | 仅新判据 | 仅旧判据 | 都不通过 |', '|---|---|---|---|---|',
          f"| 新稳态 vs 归档判据（同一 300 h 轨迹） | {c1['both_pass']} | {c1['new_only']} | "
          f"{c1['old_only']} | {c1['both_fail']} |",
          f"| 新稳态 vs 归档 100 h 扫描结果 | {c2['both_pass']} | {c2['new_only']} | "
          f"{c2['old_only']} | {c2['both_fail']} |", '',
          '## 失败原因分布（稳态）', '',
          '| 原因 | 点数 |', '|---|---|']
    for k, v in sorted(summary['failure_reasons']['steady_state'].items(), key=lambda kv: -kv[1]):
        md.append(f'| {k or "(通过)"} | {v} |')
    md += ['', '## 首次通过所需烧入长度（-1 = 从未通过）', '',
           '| first_pass_drop | 点数 |', '|---|---|']
    for k, v in summary['first_pass_drop_histogram'].items():
        md.append(f'| {k} | {v} |')
    md += ['', '## 各轴通过率（稳态）', '']
    for k in KEYS:
        md += [f'### {k}', '', '| 取值 | 点数 | 冷启动通过率 | 稳态通过率 |', '|---|---|---|---|']
        for v, r in summary['by_axis'][k].items():
            md.append(f"| {v} | {r['points']} | {r['cold_rate']:.1%} | {r['steady_rate']:.1%} |")
        md.append('')
    md += ['## 本次扫描能回答 / 不能回答', '',
           '能回答：bit0 在统一有限读窗判据下，**原网格内部**的可行窗口有多宽、各处通过率与失败原因。', '',
           '不能回答：γ 是否包含生长稀释；A0/F0 延迟是否生物学合理；两比特接口是否通过；'
           '是否可以接入 bit2。这四项必须独立保留，不能因为通过率上升就合并成“计数器已经成功”。', '',
           '**重要限定**：本网格只含 `uM_per_au ∈ {0.1, 0.3, 1, 3, 10}` 与成熟半衰期 '
           '`{1, 2.5, 5, 10, 20} min`。此前独立验证过的可用带（`uM_per_au ≈ 6–10`、成熟 12–32 min）'
           '大部分不在此网格内，因此上面的通过率描述的是**原网格**，不是可行窗口的真实宽度。', '']
    out_md = SCAN_DIR / 'rescan800_summary.md'
    out_md.write_text('\n'.join(md), encoding='utf-8')

    manifest = {}
    for p in sorted(SCAN_DIR.glob('*')):
        if p.is_file() and not p.name.startswith('SHA256SUMS'):
            # the manifest cannot contain a stable hash of itself
            manifest[str(p.relative_to(ROOT))] = sha256(p)
    for name in ('model.py', 'bit0_diagnostic.py', 'verify_twobit_causal.py',
                 'rescan800_readwindow.py', 'run_rescan800.ps1', 'merge_rescan800.py',
                 'verify_rescan_sample.py'):
        p = ROOT / name
        if p.exists():
            manifest[str(p.relative_to(ROOT))] = sha256(p)
    man_path = SCAN_DIR / 'SHA256SUMS.json'
    man_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True), encoding='utf-8')
    lines = [f'{v}  {k}' for k, v in sorted(manifest.items(), key=lambda kv: kv[1])]
    (SCAN_DIR / 'SHA256SUMS.txt').write_text('\n'.join(lines) + '\n', encoding='utf-8')

    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print(f'wrote {out_csv}')
    print(f'wrote {out_json}')
    print(f'wrote {out_md}')
    print(f'wrote {man_path}')


if __name__ == '__main__':
    main()
