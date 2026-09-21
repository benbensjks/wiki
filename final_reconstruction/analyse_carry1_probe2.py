"""Analysis and figures for the stage-2 probe (maturation 70/85 min).

Combines the stage-1 grid (maturation 5-60 min) with the stage-2 probe
(70/85 min) so the response of bit2 to the common A1/F1 maturation time can be
read as a whole, and tests the pre-registered decision rule.

    python analyse_carry1_probe2.py
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

S1 = ROOT / 'threebit51_results' / 'carry1_scan_dsh'
S2 = ROOT / 'threebit51_results' / 'carry1_probe2'
RDF_TAU_H = 1.0 / 0.8          # 1/gamma_rdf
STAGE1_BEST_CROSSINGS = 12


def sha256(p):
    return hashlib.sha256(p.read_bytes()).hexdigest().upper()


def series(frame, column):
    out = {}
    for m in sorted(set(frame.carry1_mrna_min)):
        sub = frame[frame.carry1_mrna_min == m].sort_values('carry1_maturation_min')
        out[m] = list(zip(sub.carry1_maturation_min, sub[column]))
    return out


def main():
    s1 = pd.read_csv(S1 / 'carry1_all.csv')
    s2 = pd.read_csv(S2 / 'probe_all.csv')
    probe = json.loads((S2 / 'probe_summary.json').read_text(encoding='utf-8'))
    control = probe['control_check']
    new_points = s2[s2.carry1_maturation_min >= 70.0]
    best_probe = float(new_points.bit2_crossings.max())
    control_crossings = float(control['crossings'])

    mrna_colors = {0.5: '#999999', 1.0: '#56B4E9', 2.0: '#009E73',
                   3.0: '#0072B2', 4.0: '#D55E00', 6.0: '#CC79A7', 8.0: '#F0E442'}
    all_mrna = sorted(set(s1.carry1_mrna_min) | set(s2.carry1_mrna_min))

    def merged(column):
        out = {}
        for m in all_mrna:
            pts = series(s1, column).get(m, []) + series(s2, column).get(m, [])
            if pts:
                out[m] = sorted(pts)
        return out

    fig, axes = plt.subplots(1, 3, figsize=(16.5, 4.8), constrained_layout=True)

    ax = axes[0]
    for m, pts in merged('bit2_crossings').items():
        ax.plot([p[0] for p in pts], [p[1] for p in pts], 'o-', ms=4,
                color=mrna_colors.get(m, 'k'), label=f'mRNA {m:g} min')
    ax.axhline(14, color='k', ls='--', lw=1, label='needed: 14 = bit1 reverse events')
    ax.axvline(60, color='grey', ls=':', lw=1)
    ax.set_xlabel('A1/F1 common maturation (min)')
    ax.set_ylabel('bit2 threshold crossings')
    ax.set_title('the rise is followed by a cliff, not a plateau')
    ax.legend(fontsize=7)
    ax.grid(alpha=.25)

    ax = axes[1]
    for m, pts in merged('gate_duration_median_h').items():
        pts = [(x, y) for x, y in pts if np.isfinite(y)]
        ax.plot([p[0] for p in pts], [p[1] for p in pts], 'o-', ms=4,
                color=mrna_colors.get(m, 'k'), label=f'mRNA {m:g} min')
    ax.axhline(RDF_TAU_H, color='k', ls='--', lw=1, label=f'RDF2 decay time {RDF_TAU_H:.2f} h')
    ax.axvspan(60, 70, color='red', alpha=.10)
    ax.set_xlabel('A1/F1 common maturation (min)')
    ax.set_ylabel('median g1 window width (h)')
    ax.set_title('the collapse sits just past the RDF2 decay time')
    ax.legend(fontsize=7)
    ax.grid(alpha=.25)

    ax = axes[2]
    for m, pts in merged('RDF2_median').items():
        ax.semilogy([p[0] for p in pts], np.maximum([p[1] for p in pts], 1e-3), 'o-', ms=4,
                    color=mrna_colors.get(m, 'k'), label=f'mRNA {m:g} min')
    ax.set_xlabel('A1/F1 common maturation (min)')
    ax.set_ylabel('median RDF2 (a.u.)')
    ax.set_title('the RDF2 pool that reverse recombination needs')
    ax.legend(fontsize=7)
    ax.grid(alpha=.25, which='both')

    fig.suptitle('Stage-2 probe: common A1/F1 maturation is a narrow window, not a monotone knob',
                 fontsize=12)
    for ext in ('png', 'pdf'):
        fig.savefig(S2 / f'probe_response.{ext}', dpi=200, bbox_inches='tight', facecolor='white')
    plt.close(fig)

    rule = [
        ('crossings keep rising and a point certifies (new points only)',
         'NOT FIRED' if not probe['trend']['certified'] else 'FIRED'),
        (f'crossings saturate without exceeding the stage-1 best of {STAGE1_BEST_CROSSINGS}',
         (f'FIRED (worse: new points collapse to {best_probe:.0f}, control stays at '
          f'{control_crossings:.0f})') if best_probe <= STAGE1_BEST_CROSSINGS else 'not fired'),
        ('gate contrast degrades further',
         'not fired (contrast improves: %.4f -> %.4f)'
         % (control['off_on'], float(new_points.off_on_peak_ratio.median()))),
    ]

    out = dict(control_reproduces=control['reproduces'], control=control,
               probe_best_crossings=best_probe,
               stage1_best_crossings=STAGE1_BEST_CROSSINGS,
               probe_grid=sorted({(float(r.carry1_mrna_min), float(r.carry1_maturation_min))
                                  for r in new_points.itertuples()}),
               stage2_rows=s2[['carry1_mrna_min', 'carry1_maturation_min', 'bit2_crossings',
                               's2_cycle_min_median', 's2_cycle_min_min', 'I2_peak',
                               'Int2_source_peak', 'RDF2_median', 'RDF2_max',
                               'gate_duration_median_h', 'gate_delay_median_h',
                               'rev_over_fwd_median', 'off_on_peak_ratio',
                               'gate_contrast_passed', 'steady_sequence']]
               .to_dict(orient='records'),
               decision_rule=rule,
               verdict=('Common maturation is exhausted: 60 min is the optimum of the scanned '
                        'range and 70/85 min collapse (crossings 12 -> 1, reverse/forward flux '
                        'ratio 1.13 -> 0.003-0.027, free RDF2 median 0.42 -> 0.05-0.12) while '
                        'the Int2 pulse grows. The next lever must sharpen the carry pulse '
                        'rather than widen it, i.e. split the A1 and F1 maturation times. '
                        'Mechanism wording: free RDF2 is lowered by Int2 binding it into C2, by '
                        'sequestration and by complex loss, not by DNA reverse recombination; '
                        'the reverse drive depends on the simultaneous product I2*R2, so a '
                        'falling R2 alone is not proof of substrate shortage.'),
               files={p.name: sha256(p) for p in sorted(S2.glob('*')) if p.is_file()
                      and p.suffix != '.log'})
    (S2 / 'probe_analysis.json').write_text(json.dumps(out, ensure_ascii=False, indent=2,
                                                       default=str), encoding='utf-8')

    md = ['# 第二轮探针（A1/F1 共同成熟 70 / 85 min）：结果与判定', '',
          f'- 对照点 4/60 复现：**{control["reproduces"]}**（穿越数 {control["crossings"]:.0f}，'
          f'门宽 {control["gate_width_h"]:.3f} h，与第一轮逐值一致）',
          f'- 探针 6 个新点全部 **穿越数 = 1**（第一轮最优为 12）', '',
          '## 与第一轮合并后的响应（bit2 穿越数）', '',
          '| mRNA/min | 成熟 5→60 段（第一轮） | 70 | 85 |', '|---|---|---|---|']
    for m in all_mrna:
        a = s1[s1.carry1_mrna_min == m].sort_values('carry1_maturation_min')
        b = s2[(s2.carry1_mrna_min == m) & (s2.carry1_maturation_min == 70)]
        c = s2[(s2.carry1_mrna_min == m) & (s2.carry1_maturation_min == 85)]
        if not len(a):
            continue
        md.append(f"| {m:g} | {' '.join(f'{int(x)}' for x in a.bit2_crossings)} | "
                  f"{int(b.bit2_crossings.iloc[0]) if len(b) else '—'} | "
                  f"{int(c.bit2_crossings.iloc[0]) if len(c) else '—'} |")
    md += ['', '## 预登记判据的触发情况', '', '| 判据 | 结果 |', '|---|---|']
    for cond, res in rule:
        md.append(f'| {cond} | {res} |')
    md += ['', '## 结论', '', out['verdict'], '',
           '塌陷的证据（同一批 600 h 运行）：', '',
           '1. 反向/正向通量比中位数：**1.13 → 0.003–0.027**（这是由同一时刻的'
           '`I2·R2` / `C2` 直接算出的通量比，不是靠 R2 单值推断）',
           '2. 自由 RDF2 中位数：**0.42 → 0.05–0.12**。注意机制表述：自由 RDF2 是被 '
           '`I2 + RDF2 → C2` 的结合所隔离、以及复合物降解而下降的，'
           '**DNA 反向重组本身不消耗 RDF2**；而且反向速率取决于同一时刻的乘积 `I2·R2`'
           '（等价地 `C2/K_complex`），不能只看 R2',
           f'3. 门宽：**1.33 h → 1.47–1.73 h**，而 RDF2 衰减时间为 **{RDF_TAU_H:.2f} h**——'
           '脉冲一旦明显长于底物寿命，脉冲越强反而越无效',
           '',
           '注意 I2 峰值与 Int2 源峰值在塌陷点是**升高**的（2.37→2.5–3.2 / 9.25→9.5–11.2），'
           '所以这不是"剂量不足"，而是脉冲与底物存量的时序错配。'
           '另外：停滞时刻观测到的 `I2 = 0` 只是脉冲结束后的静止状态，'
           '不能用来解释脉冲期间为什么失败——这正是第三轮扫描改用'
           '事件对齐指标（门开启时的 R2、I2 峰值时的 R2、`max(I2·R2/1.2)`、'
           '`max(C2/K_complex)`）的原因。',
           '', '![response](probe_response.png)', '',
           '## 边界声明', '',
           '本探针只改变新增的 A1/F1 表达参数；前 34 状态、uM、bit 内部参数、'
           '门 `H(Int0;0.4,3)` 与 ZENG 全部冻结；第一轮扫描器与其判据未做任何修改。', '']
    (S2 / 'probe_report.md').write_text('\n'.join(md), encoding='utf-8')
    print(json.dumps({'control_reproduces': control['reproduces'],
                      'probe_best_crossings': best_probe,
                      'stage1_best_crossings': STAGE1_BEST_CROSSINGS,
                      'decision_rule': rule}, ensure_ascii=False, indent=2))
    print('wrote', S2 / 'probe_report.md')


if __name__ == '__main__':
    main()
