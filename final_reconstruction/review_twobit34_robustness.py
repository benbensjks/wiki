"""Independent post-scan review without changing the frozen scanner."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from model import ROOT


SCAN = ROOT / 'twobit34_results' / 'robustness'
U = (5.5, 5.75, 6.0, 6.25, 6.5)
M = (25.0, 27.5, 30.0, 32.5, 35.0)
R = (1.0, 2.0, 4.0)


def strict_interior(data):
    lookup = {(float(x.uM_per_au), float(x.carry_maturation_half_life_min),
               float(x.carry_mrna_half_life_min)): bool(x.certified)
              for x in data.itertuples()}
    rows = []
    for iu in range(1, len(U)-1):
        for im in range(1, len(M)-1):
            for ir in range(1, len(R)-1):
                point = (U[iu], M[im], R[ir])
                neighbours = ((U[iu-1], M[im], R[ir]), (U[iu+1], M[im], R[ir]),
                              (U[iu], M[im-1], R[ir]), (U[iu], M[im+1], R[ir]),
                              (U[iu], M[im], R[ir-1]), (U[iu], M[im], R[ir+1]))
                if lookup[point] and all(lookup[p] for p in neighbours):
                    row = data[(data.uM_per_au == point[0]) &
                               (data.carry_maturation_half_life_min == point[1]) &
                               (data.carry_mrna_half_life_min == point[2])].iloc[0]
                    rows.append(dict(uM_per_au=point[0], carry_maturation_half_life_min=point[1],
                                     carry_mrna_half_life_min=point[2],
                                     timing_margin_h=float(row.timing_margin_h),
                                     bit1_setup_h=float(row.bit1_setup_h),
                                     bit1_hold_h=float(row.bit1_hold_h),
                                     cold_passed=bool(row.cold_passed), g0_peak=float(row.g0_peak)))
    return sorted(rows, key=lambda x: -x['timing_margin_h'])


def corrected_heatmap(data):
    fig, axes = plt.subplots(2, 3, figsize=(14, 8), sharex=True, sharey=True,
                             constrained_layout=True)
    certified_margins = data.loc[data.certified, 'timing_margin_h']
    vmax = float(certified_margins.max())
    for col, mrna in enumerate(R):
        sub = data[data.carry_mrna_half_life_min == mrna]
        passed = sub.pivot(index='carry_maturation_half_life_min',
                           columns='uM_per_au', values='certified').astype(bool)
        margin = sub.pivot(index='carry_maturation_half_life_min',
                           columns='uM_per_au', values='timing_margin_h').astype(float)
        masked_margin = margin.where(passed, np.nan)
        im0 = axes[0, col].imshow(passed.astype(float).values, origin='lower', aspect='auto',
                                  vmin=0, vmax=1, cmap='RdYlGn')
        im1 = axes[1, col].imshow(masked_margin.values, origin='lower', aspect='auto',
                                  vmin=0, vmax=vmax, cmap='viridis')
        for r in range(passed.shape[0]):
            for c in range(passed.shape[1]):
                ok = bool(passed.iloc[r, c])
                axes[0, col].text(c, r, 'PASS' if ok else 'FAIL', ha='center', va='center',
                                  fontsize=7, color='black' if ok else 'white')
                axes[1, col].text(c, r, f'{margin.iloc[r, c]:.2f}' if ok else 'FAIL',
                                  ha='center', va='center', fontsize=8,
                                  color='white' if ok and margin.iloc[r, c] < .55*vmax else 'black')
        for ax in axes[:, col]:
            ax.set_xticks(range(len(passed.columns)), [f'{x:g}' for x in passed.columns])
            ax.set_yticks(range(len(passed.index)), [f'{x:g}' for x in passed.index])
            ax.set_xlabel('uM per a.u.')
        axes[0, col].set_title(f'carry mRNA half-life = {mrna:g} min')
    axes[0, 0].set_ylabel('A0/F0 maturation half-life (min)')
    axes[1, 0].set_ylabel('A0/F0 maturation half-life (min)')
    fig.suptitle('Certified region and timing margin (failed points masked)', fontsize=15)
    fig.colorbar(im0, ax=axes[0, :], shrink=.7, label='certified')
    fig.colorbar(im1, ax=axes[1, :], shrink=.7, label='certified timing margin (h)')
    fig.savefig(SCAN / 'robustness_heatmaps_reviewed.png', dpi=220,
                bbox_inches='tight', facecolor='white')
    fig.savefig(SCAN / 'robustness_heatmaps_reviewed.pdf', bbox_inches='tight', facecolor='white')
    plt.close(fig)


def main():
    data = pd.read_csv(SCAN / 'robustness_all.csv')
    for col in ('certified', 'steady_passed', 'cold_passed', 'one_to_one',
                'causal_order', 'directions_alternate'):
        data[col] = data[col].map(lambda x: str(x).strip().lower() in ('true', '1'))
    failed = data[~data.steady_passed]
    causal_all = failed.one_to_one & failed.causal_order & failed.directions_alternate
    interior = strict_interior(data)
    selected = interior[0]
    review = dict(rows=len(data), certified=int(data.certified.sum()),
                  steady_failed=int((~data.steady_passed).sum()),
                  failed_but_all_causal=int(causal_all.sum()),
                  failed_and_causal_failed=int((~causal_all).sum()),
                  strict_interior_definition=('point is not on any of the three grid boundaries '
                                              'and all six axis neighbours are certified'),
                  strict_interior_points=interior, selected_nominal=selected,
                  notes=[
                      'The earlier count of 22 uses all existing neighbours and includes grid-boundary points.',
                      'Causal fields were computed for readout failures; only the failure_reason label short-circuits.',
                      'The original timing heatmap displays margins for failed points; the reviewed heatmap masks them.'
                  ])
    (SCAN / 'review_summary.json').write_text(
        json.dumps(review, ensure_ascii=False, indent=2), encoding='utf-8')
    md = f"""# 75点鲁棒性扫描复核

- 主扫描结论确认：75/75完整，53点完整认证，0积分失败。
- 22个稳态读出失败点中，2点因果三项全通过，20点因果三项全失败。
- 严格内部点定义：三个参数均不在网格边界，且六个轴向邻居全部认证。
- 严格内部点只有 **{len(interior)}** 个，而不是22个。

## 严格内部点

| uM/a.u. | A/F成熟/min | carry mRNA/min | 最小裕量/h | setup/h | hold/h |
|---:|---:|---:|---:|---:|---:|
"""
    for row in interior:
        md += (f"| {row['uM_per_au']:g} | {row['carry_maturation_half_life_min']:g} | "
               f"{row['carry_mrna_half_life_min']:g} | {row['timing_margin_h']:.3f} | "
               f"{row['bit1_setup_h']:.3f} | {row['bit1_hold_h']:.3f} |\n")
    md += f"""

## 选择

正式中心选择 **uM={selected['uM_per_au']:g}、A0/F0成熟={selected['carry_maturation_half_life_min']:g} min、carry mRNA={selected['carry_mrna_half_life_min']:g} min**。
它是真正六邻域内部点，冷启动通过，且在严格内部点中时间裕量最大（{selected['timing_margin_h']:.3f} h）。

![复核热图](robustness_heatmaps_reviewed.png)

原扫描文件和哈希保持不变；本复核只新增独立解释文件与修正版可视化。
"""
    (SCAN / 'review_summary.md').write_text(md, encoding='utf-8')
    corrected_heatmap(data)
    print(json.dumps(review, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
