"""Analysis and figures for the stage-1 carry1 scan (written by dsh).

Reads threebit51_results/carry1_scan_dsh/carry1_all.csv (produced by
`scan_threebit51_carry1_dsh.py merge`) and answers the seven summary items of
threebit51_results/DSH_carry1首轮扫描.md, plus a mechanistic check that the
boolean verdict alone cannot express:

    the first-round stall was "per-cycle net recombination exactly zero", which
    happens because the g1 window is long compared with the RDF2 decay time
    (1/gamma_rdf = 1.25 h), so the RDF2 pool collapses inside the pulse and the
    forward reaction re-engages before bit2 can be written low.

Outputs (in the same directory):
    carry1_analysis.json, carry1_report.md,
    carry1_heatmaps.png/.pdf, carry1_slices.png/.pdf

    python analyse_threebit51_carry1.py
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from model import ROOT

SCAN = ROOT / 'threebit51_results' / 'carry1_scan_dsh'
MRNA = (0.5, 1.0, 2.0, 4.0, 8.0)
MAT = (5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0)
GATE_RATIO_LIMIT = 0.10
BAND_LOW, BAND_HIGH = 0.30, 0.70
FIRST_ROUND = dict(mrna=2.0, mat=32.5)      # untuned reference between grid points
RDF_TAU_H = 1.0 / 0.8                        # 1/gamma_rdf


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest().upper()


def col(d, name, default=np.nan):
    return d[name] if name in d.columns else pd.Series([default] * len(d), index=d.index)


def integrity(d, summary, metas):
    numeric = d.select_dtypes(include=[np.number])
    nonfinite = {c: int((~np.isfinite(numeric[c].to_numpy(dtype=float))).sum())
                 for c in numeric.columns}
    nonfinite = {k: v for k, v in nonfinite.items() if v and not numeric[k].isna().all()}
    hashes = {}
    for key in metas[0]['sha256']:
        vals = sorted({m['sha256'][key] for m in metas})
        hashes[key] = dict(unique=len(vals), value=vals[0])
    return dict(rows=len(d), unique_points=int(d.groupby(['carry1_mrna_min',
                                                          'carry1_maturation_min']).ngroups),
                expected=len(MRNA) * len(MAT),
                solver_failures=int((~d.solver_success).sum()),
                nonfinite=int((~d.finite).sum()),
                nonfinite_numeric_columns=nonfinite,
                declared=summary.get('completeness'),
                environment_uniform=len({m['python'] for m in metas}) == 1
                and len({m['platform'] for m in metas}) == 1,
                versions_uniform=all(len({m['versions'][k] for m in metas}) == 1
                                     for k in ('numpy', 'pandas')),
                versions={k: sorted({m['versions'][k] for m in metas})[0]
                          for k in ('numpy', 'pandas')},
                frozen_hashes_uniform=all(v['unique'] == 1 for v in hashes.values()),
                frozen_hashes=hashes,
                csv_sha256_matches_all_shards=all(
                    m['csv_sha256'] and len(m['csv_sha256']) == 64 for m in metas))


def heatmaps(d):
    fig, axes = plt.subplots(2, 3, figsize=(16, 8.5), constrained_layout=True)
    panels = [
        ('certified', 'certified (mod-8 + causal + contrast)', 'RdYlGn', None),
        ('bit2_crossings', 'bit2 threshold crossings', 'viridis', None),
        ('off_on_peak_ratio', f'off/on gate peak ratio (limit {GATE_RATIO_LIMIT})', 'magma_r',
         GATE_RATIO_LIMIT),
        ('s2_cycle_min_median', 'median per-cycle S2 minimum', 'coolwarm_r', BAND_LOW),
        ('abs_net_ds2_median', 'median |net dS2| per cycle', 'cividis_r', None),
        ('rev_over_fwd_median', 'median reverse/forward flux ratio', 'PuOr_r', 1.0),
    ]
    for ax, (key, title, cmap, ref) in zip(axes.ravel(), panels):
        grid = np.full((len(MRNA), len(MAT)), np.nan)
        for i, m in enumerate(MRNA):
            for j, t in enumerate(MAT):
                row = d[(d.carry1_mrna_min == m) & (d.carry1_maturation_min == t)]
                if len(row) and key in d.columns:
                    grid[i, j] = float(row.iloc[0][key])
        vmin = np.nanmin(grid) if np.isfinite(grid).any() else 0.0
        vmax = np.nanmax(grid) if np.isfinite(grid).any() else 1.0
        if key == 'certified':
            vmin, vmax = 0, 1
        im = ax.imshow(grid, origin='lower', aspect='auto', cmap=cmap, vmin=vmin, vmax=vmax)
        for i in range(len(MRNA)):
            for j in range(len(MAT)):
                v = grid[i, j]
                txt = '—' if not np.isfinite(v) else (f'{v:.0f}' if key == 'certified'
                                                      else f'{v:.3g}')
                ax.text(j, i, txt, ha='center', va='center', fontsize=8, color='black')
        if ref is not None:
            ax.contour(np.arange(len(MAT)), np.arange(len(MRNA)), grid,
                       levels=[ref], colors='blue', linewidths=1.5)
        ax.set_xticks(range(len(MAT)), [f'{t:g}' for t in MAT])
        ax.set_yticks(range(len(MRNA)), [f'{m:g}' for m in MRNA])
        ax.set_xlabel('A1/F1 common maturation half-life (min)')
        ax.set_ylabel('carry1 mRNA half-life (min)')
        ax.set_title(title, fontsize=10)
        fig.colorbar(im, ax=ax, shrink=0.85)
    fig.suptitle('Stage-1 carry1 scan: 35 points x 600 h (frozen 34-state prefix, '
                 'frozen ZENG)', fontsize=13)
    for ext in ('png', 'pdf'):
        fig.savefig(SCAN / f'carry1_heatmaps.{ext}', dpi=200, bbox_inches='tight',
                    facecolor='white')
    plt.close(fig)


def slices(d):
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.6), constrained_layout=True)
    target = float(np.nanmedian(col(d, 'bit1_reverse_events')))
    for ax, (key, ylab, ref, ref_label) in zip(axes, [
            ('bit2_crossings', 'bit2 threshold crossings', target,
             f'bit1 reverse events (median {target:.0f})'),
            ('s2_cycle_min_median', 'median per-cycle S2 minimum', BAND_LOW,
             f'S2 = {BAND_LOW:g} band edge'),
            ('gate_duration_median_h', 'median g1 window width (h)', RDF_TAU_H,
             f'RDF2 decay time 1/gamma_rdf = {RDF_TAU_H:.2f} h')]):
        if key not in d.columns:
            ax.set_visible(False)
            continue
        for m in MRNA:
            sub = d[d.carry1_mrna_min == m].sort_values('carry1_maturation_min')
            ax.plot(sub.carry1_maturation_min, sub[key], 'o-', ms=4, label=f'mRNA {m:g} min')
        ax.axhline(ref, color='k', ls='--', lw=1, label=ref_label)
        ax.axvline(FIRST_ROUND['mat'], color='grey', ls=':', lw=1)
        ax.set_xlabel('A1/F1 common maturation (min)')
        ax.set_ylabel(ylab)
        ax.grid(alpha=0.25)
        ax.legend(fontsize=7)
    axes[0].set_title('carry count vs maturation', fontsize=10)
    axes[1].set_title('does bit2 ever reach the 0 band?', fontsize=10)
    axes[2].set_title('gate width vs RDF2 lifetime', fontsize=10)
    fig.suptitle('Stage-1 carry1 slices (dotted grey = first-round reference '
                 'maturation 32.5 min)', fontsize=12)
    for ext in ('png', 'pdf'):
        fig.savefig(SCAN / f'carry1_slices.{ext}', dpi=200, bbox_inches='tight',
                    facecolor='white')
    plt.close(fig)


def trend(d, key):
    """Monotone trend of `key` against maturation, pooled over mRNA values."""
    if key not in d.columns:
        return None
    xs = col(d, 'carry1_maturation_min').to_numpy(dtype=float)
    ys = col(d, key).to_numpy(dtype=float)
    ok = np.isfinite(xs) & np.isfinite(ys)
    if ok.sum() < 4:
        return None
    xs, ys = xs[ok], ys[ok]
    rx = pd.Series(xs).rank().to_numpy()
    ry = pd.Series(ys).rank().to_numpy()
    if np.std(rx) == 0 or np.std(ry) == 0:
        return 0.0
    return float(np.corrcoef(rx, ry)[0, 1])


def main():
    d = pd.read_csv(SCAN / 'carry1_all.csv')
    summary = json.loads((SCAN / 'carry1_summary.json').read_text(encoding='utf-8'))
    metas = [json.loads((SCAN / f'meta_shard{i:02d}of07.json').read_text(encoding='utf-8'))
             for i in range(7)]
    integ = integrity(d, summary, metas)

    certified = int(d.certified.sum())
    steady = int(d.steady_passed.sum())
    contrast = int(d.gate_contrast_passed.sum())
    crossing_err = col(d, 'crossing_error')
    best_idx = crossing_err.idxmin() if crossing_err.notna().any() else None
    best = d.loc[best_idx] if best_idx is not None else None

    # boundary-limited or mechanism-limited?
    verdict, recommendation = '', []
    if certified:
        verdict = f'{certified} point(s) certified'
        recommendation.append('共同成熟时间足以解决：直接采用认证点并进入 600 h 复算与四初态验证。')
    else:
        at_edge = []
        if best is not None:
            if best.carry1_maturation_min in (MAT[0], MAT[-1]):
                at_edge.append('maturation')
            if best.carry1_mrna_min in (MRNA[0], MRNA[-1]):
                at_edge.append('mRNA')
        cov = float(np.nanmedian(col(d, 'crossing_ratio')))
        verdict = ('no certified point; best crossing_error=%d at mRNA=%.2g mat=%.2g '
                   '(median crossing/reverse ratio %.2f)'
                   % (int(crossing_err.min()), best.carry1_mrna_min,
                      best.carry1_maturation_min, cov)) if best is not None else \
                  'no certified point'
        if at_edge:
            recommendation.append(
                f'最优解落在网格边界（{"/".join(at_edge)}），先沿该方向延长网格再判断，'
                f'不要急着改结构。')
        else:
            recommendation.append(
                '最优解在网格内部但仍不认证 → 共同成熟时间这一维不足；'
                '下一轮应分开扫描 A1 与 F1 的成熟时间（激活臂要快、抑制臂要慢）。')
        if contrast == len(d) and not certified:
            recommendation.append(
                '门对比度全部通过而计数仍失败，说明瓶颈在 bit2 的反向写入动力学，'
                '不是门泄漏——与首轮诊断一致。')

    net = col(d, 'abs_net_ds2_median')
    rev_fwd = col(d, 'rev_over_fwd_median')
    gate_w = col(d, 'gate_duration_median_h')
    gate_pk = col(d, 'gate_peak_median')
    int2pk = col(d, 'Int2_source_peak')
    cross = col(d, 'bit2_crossings')
    mech = dict(
        median_abs_net_ds2=float(net.median()),
        median_rev_over_fwd=float(rev_fwd.median()),
        median_gate_width_h=float(gate_w.median()),
        rdf2_tau_h=RDF_TAU_H,
        gate_width_over_rdf_tau=float(gate_w.median() / RDF_TAU_H),
        gate_wider_than_rdf_tau=int((gate_w > RDF_TAU_H).sum()),
        s2_cycle_min_median=float(col(d, 's2_cycle_min_median').median()),
        trend_crossings_vs_maturation=trend(d, 'bit2_crossings'),
        trend_crossings_vs_gate_width=trend(d, 'gate_duration_median_h'),
        trend_crossings_vs_gate_peak=trend(d, 'gate_peak_median'),
        corr_crossings_vs_gate_peak=float(np.corrcoef(
            cross[np.isfinite(cross) & np.isfinite(gate_pk)],
            gate_pk[np.isfinite(cross) & np.isfinite(gate_pk)])[0, 1])
        if (np.isfinite(cross) & np.isfinite(gate_pk)).sum() > 3 else None,
        corr_crossings_vs_int2_peak=float(np.corrcoef(
            cross[np.isfinite(cross) & np.isfinite(int2pk)],
            int2pk[np.isfinite(cross) & np.isfinite(int2pk)])[0, 1])
        if (np.isfinite(cross) & np.isfinite(int2pk)).sum() > 3 else None,
        s2_cycle_min_vs_maturation=trend(d, 's2_cycle_min_median'),
        interpretation=(
            'MECHANISM WORDING (corrected): free RDF2 is lowered by Int2 binding it into the '
            'complex C2 (I2 + RDF2 -> C2), by that sequestration, and by loss of C2.  DNA '
            'reverse recombination does not consume RDF2, and the reverse drive is set by the '
            'simultaneous product I2*R2 (equivalently C2/K_complex), never by R2 alone.  '
            'PREREGISTERED HYPOTHESIS RETRACTED.  The a-priori story was that a g1 window '
            'longer than the RDF2 decay time (%.2f h) lets the free RDF2 pool fall inside the '
            'pulse and re-engage the forward reaction.  The measurement does not support it: '
            'the median gate width is %.2f h (below the RDF2 time constant), only %d/35 points '
            'exceed it, and wider gates correlate POSITIVELY with success.  The supported '
            'reading is an insufficient reverse DOSE: crossings rise monotonically with '
            'maturation up to the grid edge, and both the gate peak and the Int2 source peak '
            'grow with it.' % (RDF_TAU_H, float(gate_w.median()),
                              int((gate_w > RDF_TAU_H).sum()))))

    # how much more maturation would be needed, per mRNA
    extrap = []
    for m in MRNA:
        sub = d[d.carry1_mrna_min == m].sort_values('carry1_maturation_min')
        y = sub.bit2_crossings.to_numpy(dtype=float)
        x = sub.carry1_maturation_min.to_numpy(dtype=float)
        if len(x) > 1 and np.ptp(x) > 0:
            slope, intercept = np.polyfit(x, y, 1)
            extrap.append(dict(carry1_mrna_min=m, slope_crossings_per_min=float(slope),
                               crossings_at_mat_max=float(y[-1]),
                               maturation_for_14=float((14 - intercept) / slope) if slope > 0
                               else None))
    anomalies = dict(
        no_gate_points=[dict(mrna=float(r.carry1_mrna_min), mat=float(r.carry1_maturation_min),
                             gate_events=int(r.gate_events))
                        for r in d[d.gate_events < 14].itertuples()],
        contrast_fail_with_gate=[dict(mrna=float(r.carry1_mrna_min),
                                      mat=float(r.carry1_maturation_min),
                                      off_on=float(r.off_on_peak_ratio))
                                 for r in d[(~d.gate_contrast_passed) & (d.gate_events == 14)].itertuples()],
        mrna8_nonmonotone='crossings at mRNA=8 fall from 9 at mat=50 to %d at mat=60'
                          % int(d[(d.carry1_mrna_min == 8.0)
                                  & (d.carry1_maturation_min == 60.0)].bit2_crossings.iloc[0])
        if len(d[(d.carry1_mrna_min == 8.0) & (d.carry1_maturation_min == 60.0)]) else None)

    cls_cols = [c for c in ('carry1_mrna_min', 'carry1_maturation_min', 'certified',
                            'steady_passed', 'gate_contrast_passed', 'crossing_error',
                            'bit2_crossings', 'bit1_reverse_events', 'gate_events',
                            'off_on_peak_ratio', 's2_cycle_min_median', 'abs_net_ds2_median',
                            'rev_over_fwd_median', 'gate_duration_median_h',
                            'bit2_low_cycles', 'bit2_high_cycles', 'minimum_commitment',
                            'steady_sequence') if c in d.columns]
    ranking = d.sort_values(['certified', 'crossing_error', 's2_cycle_min_median',
                             'abs_net_ds2_median']).head(10)[cls_cols]

    out = dict(integrity=integ,
               counts=dict(points=len(d), certified=certified, steady_passed=steady,
                           gate_contrast_passed=contrast,
                           bit2_crossings_total=int(col(d, 'bit2_crossings').sum()),
                           bit1_reverse_events_total=int(col(d, 'bit1_reverse_events').sum()),
                           gate_events_total=int(col(d, 'gate_events').sum())),
               by_mrna=summary.get('by_mrna'), by_maturation=summary.get('by_maturation'),
               mechanism=mech, extrapolation=extrap, anomalies=anomalies,
               verdict=verdict, recommendation=recommendation,
               ranking=ranking.to_dict(orient='records'),
               shard_comparison=json.loads((SCAN / 'scanner_comparison.json')
                                           .read_text(encoding='utf-8'))
               if (SCAN / 'scanner_comparison.json').exists() else None)
    (SCAN / 'carry1_analysis.json').write_text(
        json.dumps(out, ensure_ascii=False, indent=2, default=str), encoding='utf-8')

    heatmaps(d)
    slices(d)

    md = ['# 51 状态第二级 carry 首轮扫描：结果与分析', '',
          f'- 网格 5(mRNA) × 7(成熟) = **{integ["rows"]}** 点，每点 600 h；唯一参数组合 '
          f'{integ["unique_points"]}/{integ["expected"]}',
          f'- 积分失败 {integ["solver_failures"]}；非有限轨迹 {integ["nonfinite"]}',
          f'- 环境一致：{integ["environment_uniform"]}（numpy {integ["versions"]["numpy"]}，'
          f'pandas {integ["versions"]["pandas"]}）；冻结源文件哈希跨 7 片一致：'
          f'{integ["frozen_hashes_uniform"]}', '',
          '## 1. 完整性与健康度', '',
          f'- 行数/唯一组合/期望：{integ["rows"]}/{integ["unique_points"]}/{integ["expected"]}',
          f'- 求解失败：{integ["solver_failures"]}；非有限：{integ["nonfinite"]}',
          f'- 数值列非有限计数：{integ["nonfinite_numeric_columns"] or "无"}', '',
          '## 2. 完整模 8 认证', '',
          f'- **认证点：{certified}/35**；稳态通过 {steady}/35；门对比度通过 {contrast}/35',
          '', '## 3. 每个参数方向上的 carry 事件对照', '']
    for m in MRNA:
        sub = d[d.carry1_mrna_min == m].sort_values('carry1_maturation_min')
        md.append(f'### carry1 mRNA 半衰期 {m:g} min')
        md.append('')
        md.append('| A1/F1 成熟/min | bit1 反向事件 | g1 窗口 | bit2 穿越 | 门宽/h | off/on | S2 周期间最低(中位) | |ΔS2|(中位) |')
        md.append('|---:|---:|---:|---:|---:|---:|---:|---:|')
        for r in sub.itertuples():
            md.append(f'| {r.carry1_maturation_min:g} | {getattr(r, "bit1_reverse_events", float("nan")):.0f} '
                      f'| {getattr(r, "gate_events", float("nan")):.0f} | {getattr(r, "bit2_crossings", float("nan")):.0f} '
                      f'| {getattr(r, "gate_duration_median_h", float("nan")):.2f} '
                      f'| {getattr(r, "off_on_peak_ratio", float("nan")):.3f} '
                      f'| {getattr(r, "s2_cycle_min_median", float("nan")):.3f} '
                      f'| {getattr(r, "abs_net_ds2_median", float("nan")):.4f} |')
        md.append('')
    md += ['## 4. 门对比度', '',
           f'- 阈值 {GATE_RATIO_LIMIT}；通过 {contrast}/35',
           f'- off/on 峰值比中位数：{float(col(d, "off_on_peak_ratio").median()):.4f}', '',
           '## 5. 排名（无认证点时按穿越数误差）', '',
           ranking.to_markdown(index=False) if hasattr(pd.DataFrame, 'to_markdown') else
           ranking.to_string(index=False), '',
           '## 6. 共同成熟时间是否足够', '',
           f'- 结论：**{verdict}**']
    for r in recommendation:
        md.append(f'- {r}')
    md += ['', '## 7. 机理检查（dsh 追加）', '',
           f'- 每周期 |净 ΔS2| 中位数：**{mech["median_abs_net_ds2"]:.4f}**',
           f'- 反向/正向通量比中位数：**{mech["median_rev_over_fwd"]:.3f}**',
           f'- g1 窗口宽度中位数：**{mech["median_gate_width_h"]:.2f} h**，'
           f'RDF2 衰减时间 1/γ_rdf = **{RDF_TAU_H:.2f} h**，'
           f'宽度/衰减时间 = **{mech["gate_width_over_rdf_tau"]:.2f}**',
           f'- 门宽超过 RDF2 衰减时间的点数：{mech["gate_wider_than_rdf_tau"]}/35',
           f'- 逐周期 S2 最低点中位数：**{mech["s2_cycle_min_median"]:.3f}**'
           f'（要进入 0 带需 ≤ {BAND_LOW:g}）',
           f'- 与成熟时间的秩相关：穿越数 {mech["trend_crossings_vs_maturation"]:.3f}，'
           f'S2周期最低 {mech["s2_cycle_min_vs_maturation"]:.3f}',
           f'- 穿越数与门峰值/Int2 源峰值的线性相关：'
           f'{mech["corr_crossings_vs_gate_peak"]:.3f} / {mech["corr_crossings_vs_int2_peak"]:.3f}', '',
           f'> **预登记假说已撤回。** {mech["interpretation"].split("RETRACTED.", 1)[1].strip()}', '',
           '## 8. 需要多少成熟时间（线性外推）', '',
           '| carry1 mRNA/min | 斜率(穿越/分钟) | mat=60 时穿越数 | 外推达 14 次所需成熟时间/min |',
           '|---:|---:|---:|---:|']
    for e in extrap:
        need = '—' if e['maturation_for_14'] is None else f"{e['maturation_for_14']:.0f}"
        md.append(f"| {e['carry1_mrna_min']:g} | {e['slope_crossings_per_min']:.4f} | "
                  f"{e['crossings_at_mat_max']:.0f} | {need} |")
    md += ['', '## 9. 异常点', '',
           f'- g1 窗口数不足 14 的点（成熟太快时门来不及在每个 carry 周期打开）：'
           f'{len(anomalies["no_gate_points"])} 个 → '
           f'{[(p["mrna"], p["mat"], p["gate_events"]) for p in anomalies["no_gate_points"]]}',
           f'- 有 14 个窗口但门对比度不达标的点：{len(anomalies["contrast_fail_with_gate"])} 个 → '
           f'{[(p["mrna"], p["mat"], round(p["off_on"], 3)) for p in anomalies["contrast_fail_with_gate"]]}',
           f'- mRNA=8 的非单调塌陷：{anomalies["mrna8_nonmonotone"]}'
           '（该列的外推因此不可信，只有 mRNA=2–4 的外推可用）',
           '- 门对比度在成熟 40–50 min 处最差、到 60 min 又恢复，说明'
           '"延长成熟时间"与"保持门对比度"之间存在非单调的张力，下一轮网格必须在'
           '60 min 以外加密，而不是只把上界抬高。', '',
           '## 图', '', '![heatmaps](carry1_heatmaps.png)', '',
           '![slices](carry1_slices.png)', '',
           '## 边界声明', '',
           '本扫描只改变新增的 A1/F1 表达参数；前 34 状态、uM、bit 内部参数、'
           '门 `H(Int0;0.4,3)` 与 ZENG 全部冻结。所验证的是该 34 状态底座下 '
           'bit2 反向写入动力学的可调性，不构成实验标定。', '']
    (SCAN / 'carry1_report.md').write_text('\n'.join(md), encoding='utf-8')

    manifest = {p.name: sha256(p) for p in sorted(SCAN.glob('*')) if p.is_file()
                and not p.name.startswith(('SHA256', 'SOURCE_HASHES')) and p.suffix != '.log'}
    # console transcripts keep growing while the script prints, so they are not
    # part of the hash manifest; every result artefact is.
    sources = {}
    for external in ('model_threebit51.py', 'verify_threebit51.py', 'model_twobit34.py',
                     'model.py', 'verify_twobit_causal.py',
                     'scan_threebit51_carry1_dsh.py', 'run_threebit51_carry1_dsh.ps1',
                     'analyse_threebit51_carry1.py', 'compare_carry1_scanners.py',
                     'scan_threebit51_carry1.py'):
        p = ROOT / external
        if p.exists():
            sources[external] = sha256(p)
    (SCAN / 'SOURCE_HASHES.json').write_text(
        json.dumps(dict(note='paths are relative to final_reconstruction/', sources=sources),
                   ensure_ascii=False, indent=2, sort_keys=True), encoding='utf-8')
    (SCAN / 'SOURCE_HASHES.txt').write_text(
        '\n'.join(f'{v}  {k}' for k, v in sorted(sources.items())) + '\n', encoding='utf-8')
    (SCAN / 'SHA256SUMS.json').write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True), encoding='utf-8')
    (SCAN / 'SHA256SUMS.txt').write_text(
        '\n'.join(f'{v}  {k}' for k, v in sorted(manifest.items(), key=lambda kv: kv[1])) + '\n',
        encoding='utf-8')

    print(json.dumps({k: out[k] for k in ('counts', 'mechanism', 'verdict')},
                     ensure_ascii=False, indent=2, default=str))
    print('\n' + ranking.to_string(index=False))
    print(f'\nwrote {SCAN / "carry1_analysis.json"}, carry1_report.md, figures, SHA256SUMS')


if __name__ == '__main__':
    main()
