"""Analysis and figures for the stage-3 split-maturation scan (63 points).

Answers the pre-registered question: does separating the A1 (activator) and F1
(repressor) maturation times let the carry pulse deliver reverse dose while the
free RDF2 pool is still present -- or does the same substrate/timing wall bind?

Mechanism wording: free RDF2 falls because Int2 binds it into C2, sequesters it,
and C2 is lost; DNA reverse recombination does not consume RDF2.  The reverse
drive therefore depends on the simultaneous product I2*R2 (equivalently
C2/K_complex), not on R2 alone.

    python analyse_carry1_split.py
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

SCAN = ROOT / 'threebit51_results' / 'carry1_split_scan'
A1 = (10.0, 20.0, 30.0, 40.0)
F1 = (20.0, 30.0, 45.0, 60.0, 80.0)
MRNA = (3.0, 4.0, 6.0)
CONTROLS = [(60.0, 60.0, 4.0), (65.0, 65.0, 4.0), (70.0, 70.0, 4.0)]
CLASS_ORDER = ['success', 'partial_reverse_write', 'flux_tug_of_war', 'substrate_timing',
               'gate_not_formed', 'unclassified_failure', 'integration_failed', 'gate_leak']
LOW_BAND = 0.30


def sha256(p):
    return hashlib.sha256(p.read_bytes()).hexdigest().upper()


def reclassify(d, expected_gates):
    """Corrected classifier, a pure function of the columns recorded in the CSV.

    Two thresholds in the scanner were mis-specified and are corrected here:

    1. gate_not_formed compared the episode count against 14, the number of gate
       windows in a whole 600 h run, while episodes are only collected after the
       100 h burn-in (the true count is 12 for every point, because bit1 is
       frozen).  The corrected rule uses three quarters of the observed count.
    2. gate_leak used "the dwell moves S2 by more than 0.10", but a working
       toggle moves S2 across its whole range during every dwell, so that
       criterion fires on everything (observed 0.11..0.43).  It is dropped; the
       classes below are outcome based instead.

    Classes:
        success                 certified
        gate_not_formed         fewer than 3/4 of the expected carry windows
        substrate_timing        max(I2*R2)/K_D_comp never reaches 1
        partial_reverse_write   substrate reached and bit2 does enter the 0.3
                                band at least once, but the count is not certified
        flux_tug_of_war         substrate reached, bit2 never enters the 0.3 band
    """
    out = []
    for r in d.itertuples():
        if bool(r.certified):
            cls = 'success'
        elif r.gate_events < 0.75 * expected_gates:
            cls = 'gate_not_formed'
        elif not np.isfinite(r.max_I2R2_over_KD_best) or r.max_I2R2_over_KD_best < 1.0:
            cls = 'substrate_timing'
        elif r.episodes_entering_low_band >= 1:
            cls = 'partial_reverse_write'
        else:
            cls = 'flux_tug_of_war'
        out.append(cls)
    return out


def grid_of(d, mrna, column):
    g = np.full((len(A1), len(F1)), np.nan)
    sub = d[d.carry1_mrna_min == mrna]
    for i, a in enumerate(A1):
        for j, f in enumerate(F1):
            row = sub[(sub.carry1_a1_mat_min == a) & (sub.carry1_f1_mat_min == f)]
            if len(row) and column in d.columns:
                g[i, j] = float(row.iloc[0][column])
    return g


def main():
    d = pd.read_csv(SCAN / 'split_all.csv')
    summary = json.loads((SCAN / 'split_summary.json').read_text(encoding='utf-8'))
    metas = [json.loads((SCAN / f'meta_shard{i:02d}of09.json').read_text(encoding='utf-8'))
             for i in range(9)]
    expected_gates = int(d.gate_events.max())
    d['classification_raw'] = d['classification']
    d['classification'] = reclassify(d, expected_gates)
    frozen_ok = all(
        all(hashlib.sha256((ROOT / f'{k}.py').read_bytes()).hexdigest().upper() == m['sha256'][k]
            for k in m['sha256']) for m in metas)

    ctrl = d[(d.carry1_mrna_min == 4.0) &
             (d.carry1_a1_mat_min == d.carry1_f1_mat_min) &
             (d.carry1_a1_mat_min.isin([60.0, 65.0, 70.0]))].sort_values('carry1_a1_mat_min')

    # ---------------------------------------------------------------- figures
    panels = [('episodes_entering_low_band', 'carry episodes entering S2 <= 0.3', 'viridis'),
              ('max_I2R2_over_KD_best', 'max I2*R2 / K_D_comp', 'magma'),
              ('S2_min_in_episode_best', 'deepest episode S2', 'coolwarm_r'),
              ('bit2_crossings', 'bit2 threshold crossings', 'cividis')]
    fig, axes = plt.subplots(len(MRNA), len(panels), figsize=(18, 10.5),
                             constrained_layout=True, squeeze=False)
    for i, m in enumerate(MRNA):
        for j, (key, title, cmap) in enumerate(panels):
            ax = axes[i][j]
            g = grid_of(d, m, key)
            vmin = np.nanmin(g) if np.isfinite(g).any() else 0
            vmax = np.nanmax(g) if np.isfinite(g).any() else 1
            if key == 'S2_min_in_episode_best':
                vmin, vmax = 0.0, 1.0
            im = ax.imshow(g, origin='lower', aspect='auto', cmap=cmap, vmin=vmin, vmax=vmax)
            for a in range(len(A1)):
                for b in range(len(F1)):
                    v = g[a, b]
                    ax.text(b, a, '—' if not np.isfinite(v) else
                            (f'{v:.0f}' if key in ('bit2_crossings', 'episodes_entering_low_band')
                             else f'{v:.2f}'), ha='center', va='center', fontsize=8)
            if key == 'max_I2R2_over_KD_best':
                ax.contour(np.arange(len(F1)), np.arange(len(A1)), g, levels=[1.0],
                           colors='cyan', linewidths=1.6)
            if key == 'S2_min_in_episode_best':
                ax.contour(np.arange(len(F1)), np.arange(len(A1)), g, levels=[LOW_BAND],
                           colors='black', linewidths=1.6)
            ax.set_xticks(range(len(F1)), [f'{x:g}' for x in F1])
            ax.set_yticks(range(len(A1)), [f'{x:g}' for x in A1])
            ax.set_xlabel('F1 maturation (min)')
            ax.set_ylabel('A1 maturation (min)')
            if i == 0:
                ax.set_title(title, fontsize=10)
            if j == 0:
                ax.text(-0.32, 0.5, f'carry1 mRNA {m:g} min', transform=ax.transAxes,
                        rotation=90, va='center', ha='center', fontsize=11)
            fig.colorbar(im, ax=ax, shrink=0.8)
    fig.suptitle('Stage-3 split-maturation scan (63 points x 600 h): '
                 'does the pulse beat the RDF2 substrate clock?', fontsize=13)
    for ext in ('png', 'pdf'):
        fig.savefig(SCAN / f'split_heatmaps.{ext}', dpi=190, bbox_inches='tight',
                    facecolor='white')
    plt.close(fig)

    fig2, axes2 = plt.subplots(1, 3, figsize=(16.5, 4.8), constrained_layout=True)
    ax = axes2[0]
    ax.scatter(d.max_I2R2_over_KD_best, d.episodes_entering_low_band, s=26,
               c=d.bit2_crossings, cmap='viridis')
    ax.axvline(1.0, color='k', ls='--', label='substrate threshold 1.0')
    for a, f, m in CONTROLS:
        r = d[(d.carry1_mrna_min == m) & (d.carry1_a1_mat_min == a) & (d.carry1_f1_mat_min == f)]
        if len(r):
            ax.scatter(r.max_I2R2_over_KD_best, r.episodes_entering_low_band, s=70,
                       marker='*', c='red', zorder=5,
                       label=f'control {a:g}/{f:g}/{m:g}')
    ax.set_xlabel('max I2*R2 / K_D_comp reached in a carry episode')
    ax.set_ylabel('episodes entering S2 <= 0.3')
    ax.set_title('reaching the substrate threshold vs actually writing bit2 low')
    ax.legend(fontsize=7)
    ax.grid(alpha=.25)

    ax = axes2[1]
    if len(ctrl):
        x = np.arange(len(ctrl))
        ax.bar(x - 0.2, ctrl.bit2_crossings, width=0.4, label='bit2 crossings')
        ax.bar(x + 0.2, ctrl.episodes_entering_low_band, width=0.4,
               label='episodes entering S2 <= 0.3')
        ax.set_xticks(x, [f'{a:g}/{a:g}' for a in ctrl.carry1_a1_mat_min])
        ax.axhline(14, color='k', ls='--', lw=1, label='needed: 14')
        ax.set_xlabel('common A1/F1 maturation (min), mRNA 4')
        ax.set_title('common-maturation controls (60-70 min)')
        ax.legend(fontsize=7)
        ax.grid(alpha=.25, axis='y')

    ax = axes2[2]
    for cls in CLASS_ORDER:
        n = int((d.classification == cls).sum())
        if n:
            ax.barh(cls, n, color='#D55E00' if cls != 'success' else '#009E73')
            ax.text(n, cls, f' {n}', va='center', fontsize=9)
    ax.set_xlabel('points')
    ax.set_title('pre-registered classification (63 points)')
    ax.grid(alpha=.25, axis='x')
    for ext in ('png', 'pdf'):
        fig2.savefig(SCAN / f'split_synthesis.{ext}', dpi=200, bbox_inches='tight',
                     facecolor='white')
    plt.close(fig2)

    # --------------------------------------------------------- decision logic
    any_low = int((d.episodes_entering_low_band > 0).sum())
    best_low = int(d.episodes_entering_low_band.max())
    best_sub = float(d.max_I2R2_over_KD_best.max())
    certified = int(d.certified.sum())
    ctrl_best_cross = int(ctrl.bit2_crossings.max()) if len(ctrl) else 0
    split_best_cross = int(d[d.carry1_a1_mat_min != d.carry1_f1_mat_min].bit2_crossings.max())
    deepest = float(d.S2_min_in_episode_best.min())
    substrate_pts = int((d.max_I2R2_over_KD_best >= 1.0).sum())
    concord = float(((d.max_I2R2_over_KD_best >= 1.0) ==
                     (d.episodes_entering_low_band > 0)).mean())
    if certified:
        verdict = (f'{certified} point(s) certified: split maturation works. Freeze the best '
                   'point and proceed to the 600 h re-check and the eight initial states.')
    elif any_low:
        verdict = (
            f'PARTIAL: {any_low}/63 points do write bit2 into the S2 <= 0.3 band (best '
            f'{best_low} of 12 carry episodes), so the pre-registered "no point enters the low '
            f'band" trigger does NOT fire. But none of them counts: the best split point reaches '
            f'{split_best_cross} crossings against {ctrl_best_cross} for the un-split 60/60 '
            f'control, and the deepest write anywhere is only {deepest:.2f} against the 0.30 '
            f'band edge. Expression-time tuning is therefore saturated at partial reverse '
            f'writes, not exhausted at zero: one bounded local refinement is defensible, but the '
            f'evidence already points at structure (RDF2/T2 feedback or an independent gate '
            f'exponent) as the lever that can move 4/12 writes to 12/12.')
    else:
        verdict = ('No point enters the S2 <= 0.3 band: expression-time tuning is exhausted '
                   'under the frozen 34-state base. Per the pre-registered rule, stop and move '
                   'to structure (RDF2/T2 feedback) or an independent gate exponent n_A1_gate.')
    cls_counts = {str(k): int(v) for k, v in d.classification.value_counts().items()}
    cls_raw = {str(k): int(v) for k, v in d.classification_raw.value_counts().items()}

    out = dict(points=len(d), certified=certified, classification=cls_counts,
               classification_as_recorded_by_scanner=cls_raw,
               expected_gates_after_burn_in=expected_gates,
               points_reaching_substrate=substrate_pts,
               substrate_lowband_concordance=concord,
               best_split_crossings=split_best_cross,
               best_common_control_crossings=ctrl_best_cross,
               deepest_write=deepest,
               classifier_note=('Two scanner thresholds were mis-specified and are corrected '
                                'here: (1) gate_not_formed compared against the whole-run gate '
                                'count 14 while episodes start after the 100 h burn-in (true '
                                'count %d); (2) gate_leak used "the dwell moves S2 by more '
                                'than 0.10", which is true for a working toggle as well '
                                '(observed 0.11-0.43), so it is dropped and the classes are '
                                'outcome based. Both label sets are reported.'
                                % expected_gates),
               points_entering_low_band=any_low, best_episodes_entering_low_band=best_low,
               max_substrate_ratio_reached=best_sub,
               controls=ctrl[['carry1_a1_mat_min', 'bit2_crossings',
                              'episodes_entering_low_band', 'max_I2R2_over_KD_best',
                              'S2_min_in_episode_best', 'classification']]
               .to_dict(orient='records'),
               frozen_files_unchanged=frozen_ok,
               environment=dict(python=sorted({m['python'][:22] for m in metas}),
                                numpy=sorted({m['versions']['numpy'] for m in metas})),
               preregistered=metas[0]['preregistered'],
               verdict=verdict,
               files={p.name: sha256(p) for p in sorted(SCAN.glob('*')) if p.is_file()
                      and p.suffix not in ('.log',) and not p.name.startswith(('SHA256',
                                                                              'split_shard',
                                                                              'meta_'))})
    (SCAN / 'split_analysis.json').write_text(json.dumps(out, ensure_ascii=False, indent=2,
                                                         default=str), encoding='utf-8')

    md = ['# 第三轮扫描：A1/F1 分开成熟（63 点 × 600 h）', '',
          f'- 完整性 {summary["completeness"]["rows"]}/{summary["completeness"]["expected"]}，'
          f'重复 {summary["completeness"]["duplicates"]}，缺失 {len(summary["completeness"]["missing"])}',
          f'- 积分失败 {summary["solver_failures"]}；非有限 {summary["nonfinite"]}',
          f'- 冻结源文件跨 9 片一致且未被改动：**{frozen_ok}**', '',
          '## 预登记分类结果', '',
          f'修正判据下（`gate_not_formed` 阈值改为实测的烧入后窗口数 {expected_gates}）：']
    for cls in CLASS_ORDER:
        n = cls_counts.get(cls, 0)
        if n:
            md.append(f'- `{cls}`：**{n}** 点')
    md += ['', '扫描器内预登记判据（阈值误用整段 14）记录的标签：']
    for cls, n in cls_raw.items():
        md.append(f'- `{cls}`：{n} 点')
    md += ['', '> 判据修正说明（两条）：',
           f'> 1. `gate_not_formed` 原本拿整段 600 h 的 14 个窗口当阈值，但事件只从 100 h '
           f'烧入之后收集，真实值对所有点都是 **{expected_gates}**（bit1 冻结）；'
           '已修正为"少于其 75%"。',
           '> 2. `gate_leak` 原本用"驻留期 S2 移动 >0.10"，但**正常翻转在每个驻留期都会'
           '跨越整个范围**（实测 0.11–0.43），该判据对一切点都成立，已废弃；'
           '改用基于结果的分类。',
           '> 两套标签（扫描器内预登记 / 修正后）都已列出。', '']
    md += ['', '## 关键问句', '',
           f'- 有 **{any_low}** 个点至少一个周期把 S2 写进 ≤0.3（单点最多 {best_low}/12 个周期）',
           f'- 达到底物阈值 `max(I2·R2)/K_D_comp ≥ 1` 的点：**{substrate_pts}/63**',
           f'- 底物阈值与"进入低带"的一致率：**{concord:.1%}**（本轮最强的机制证据）',
           f'- 最深写入：**{deepest:.3f}**（带边 0.30）',
           f'- 分开成熟的最好穿越数：**{split_best_cross}**；未分开的 60/60 对照：'
           f'**{ctrl_best_cross}**（需要 14）', '',
           '## 共同成熟对照（mRNA=4）', '',
           '| A1/F1 | bit2 穿越 | 进入低带周期数 | max(I2·R2)/KD | 最深 S2 | 分类 |',
           '|---:|---:|---:|---:|---:|---|']
    for r in out['controls']:
        md.append(f"| {r['carry1_a1_mat_min']:g}/{r['carry1_a1_mat_min']:g} | "
                  f"{r['bit2_crossings']:.0f} | {r['episodes_entering_low_band']:.0f} | "
                  f"{r['max_I2R2_over_KD_best']:.3f} | {r['S2_min_in_episode_best']:.3f} | "
                  f"`{r['classification']}` |")
    md += ['', '## 判定', '', verdict, '',
           '## 机制表述（已修正）', '',
           '自由 RDF2 的下降来自 `I2 + RDF2 → C2` 的结合隔离与复合物降解，'
           '**DNA 反向重组不消耗 RDF2**；反向速率取决于同一时刻的乘积 `I2·R2`'
           '（等价 `C2/K_complex`），因此**不能单用 R2 阈值判断"底物不足"**。'
           '本轮所有底物指标都按事件对齐记录：g1 开启时的 R2、I2 峰值时的 R2、'
           '`max(I2·R2/1.2)`、`max(C2/K_complex)`、`∫J_rev2`、`∫J_fwd2` 与净变化。', '',
           '## 图', '', '![heatmaps](split_heatmaps.png)', '',
           '![synthesis](split_synthesis.png)', '',
           '## 边界声明', '',
           '本轮只改变新增的 A1/F1 mRNA 半衰期与成熟时间；前 34 状态、uM、bit 内部参数、'
           '门 `H(Int0;0.4,3)` 与 ZENG 全部冻结，判据仍为冻结的 '
           '`verify_threebit51.analyse_threebit`。', '']
    (SCAN / 'split_report.md').write_text('\n'.join(md), encoding='utf-8')

    manifest = {p.name: sha256(p) for p in sorted(SCAN.glob('*')) if p.is_file()
                and p.suffix != '.log' and not p.name.startswith(('SHA256', 'SOURCE_HASHES'))}
    sources = {}
    for name in ('scan_carry1_split_dsh.py', 'run_carry1_split_dsh.ps1',
                 'analyse_carry1_split.py', 'model_threebit51.py', 'verify_threebit51.py'):
        if (ROOT / name).exists():
            sources[name] = sha256(ROOT / name)
    (SCAN / 'SOURCE_HASHES.json').write_text(
        json.dumps(dict(note='relative to final_reconstruction/', sources=sources),
                   ensure_ascii=False, indent=2, sort_keys=True), encoding='utf-8')
    (SCAN / 'SHA256SUMS.json').write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True), encoding='utf-8')
    (SCAN / 'SHA256SUMS.txt').write_text(
        '\n'.join(f'{v}  {k}' for k, v in sorted(manifest.items())) + '\n', encoding='utf-8')

    print(json.dumps({k: out[k] for k in ('points', 'certified', 'classification',
                                          'points_entering_low_band',
                                          'best_episodes_entering_low_band',
                                          'max_substrate_ratio_reached',
                                          'frozen_files_unchanged', 'verdict')},
                     ensure_ascii=False, indent=2, default=str))
    print('controls:')
    print(ctrl[['carry1_a1_mat_min', 'bit2_crossings', 'episodes_entering_low_band',
                'max_I2R2_over_KD_best', 'S2_min_in_episode_best', 'classification']]
          .to_string(index=False))
    print('wrote', SCAN / 'split_report.md')


if __name__ == '__main__':
    main()
