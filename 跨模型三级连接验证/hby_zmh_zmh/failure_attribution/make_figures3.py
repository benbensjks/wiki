"""Corrected figures after the review. Does not overwrite the withdrawn ones.

  fig3_mechanism_v2        same-time flux statistic, and why it does NOT
                           separate pass from fail (replaces fig3_mechanism)
  fig4_eight_states_v2     the corrected eight DIGITAL-STATE continuations
                           (replaces fig4_eight_phases)

The earlier fig3_mechanism / fig4_eight_phases / fig6_phase_events are left on
disk unchanged and are marked withdrawn in figures/FIGURE_STATUS.md, so the
record of what was wrong stays visible.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np

sys.dont_write_bytecode = True

ROOT = Path(__file__).resolve().parent
FIG = ROOT / 'figures'
RES = ROOT / 'results'
OK, BAD, MID = '#2e7d32', '#c62828', '#ef9a00'


def newest(pat):
    d = sorted(RES.glob(pat), key=lambda p: p.name)
    if not d:
        raise SystemExit(f'no results for {pat}')
    return d[-1]


def main():
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8')
    FIG.mkdir(exist_ok=True)
    r5 = newest('round5_*')
    r5b = newest('round5b_*')

    with (r5 / 'partC_flux.csv').open(encoding='utf-8') as fh:
        rows = list(csv.DictReader(fh))
    for r in rows:
        r['ok'] = (r['counting'] == 'True' and r['events'] == 'True')
        for k in ('max_I2_times_max_RDF2', 'max_same_time_IR',
                  'max_same_time_IR_over_threshold', 'int_ratio',
                  'v_rev_over_v_fwd_at_max', 'clock_K'):
            r[k] = float(r[k])

    # =====================================================================
    # fig 3 v2
    # =====================================================================
    fig, axes = plt.subplots(1, 2, figsize=(16, 7), layout='constrained')
    ax = axes[0]
    for ok, mk, col, lab in ((True, 'o', OK, 'counts (counting AND events)'),
                             (False, 'X', BAD, 'fails')):
        sel = [r for r in rows if r['ok'] is ok]
        ax.scatter([r['max_same_time_IR_over_threshold'] for r in sel],
                   [r['max_I2_times_max_RDF2'] / max(r['max_same_time_IR'], 1e-9)
                    for r in sel],
                   s=170, marker=mk, color=col, edgecolor='k', linewidth=1.1,
                   zorder=5, label=lab)
    for r in rows:
        if r['clock_K'] in (0.40, 0.25, 0.20) and r['bypass'] == 'False' \
                or r['bypass'] == 'True':
            tag = 'bypass' if r['bypass'] == 'True' else f"K={r['clock_K']:g}"
            ax.annotate(f"{tag} {r['arm']}",
                        (r['max_same_time_IR_over_threshold'],
                         r['max_I2_times_max_RDF2'] / max(r['max_same_time_IR'], 1e-9)),
                        textcoords='offset points', xytext=(7, 6), fontsize=10)
    ax.axhline(1.0, color='gray', ls=':', lw=1.6)
    ax.text(0.11, 1.15, 'no inflation', fontsize=11, color='gray')
    ax.set_xscale('log')
    ax.set_xlabel(r'same-time drive / threshold:  $\max_t(I_2 RDF_2)\,/\,K_{D,comp}$',
                  fontsize=15)
    ax.set_ylabel('peak x peak  /  same-time', fontsize=15)
    ax.set_title('(a) how much the old indicator overstated the drive', fontsize=14)
    ax.grid(alpha=.22); ax.tick_params(labelsize=13)
    ax.legend(fontsize=12, loc='upper right')
    ax.text(.03, .06,
            'the old rule "drive must reach 5-10x threshold" is refuted:\n'
            'bypass/on counts at 2.18 while K=0.25 fails at 4.44',
            transform=ax.transAxes, fontsize=12, color=BAD, fontweight='bold')

    ax = axes[1]
    for ok, mk, col, lab in ((True, 'o', OK, 'counts'), (False, 'X', BAD, 'fails')):
        sel = [r for r in rows if r['ok'] is ok]
        ax.scatter([r['int_ratio'] for r in sel],
                   [r['max_same_time_IR_over_threshold'] for r in sel],
                   s=170, marker=mk, color=col, edgecolor='k', linewidth=1.1,
                   zorder=5, label=lab)
    ax.axvline(0.7, color=MID, ls='--', lw=2.2)
    ax.text(0.705, 0.6, 'gap 0.63 - 0.76', fontsize=11, color=MID, rotation=90)
    ax.set_yscale('log')
    ax.set_xlabel(r'integrated flux ratio  $\int v_{rev}\,dt \,/\, \int v_{fwd}\,dt$',
                  fontsize=15)
    ax.set_ylabel(r'same-time drive / threshold', fontsize=15)
    ax.set_title('(b) only the integrated ratio separates - and it restates '
                 'the outcome', fontsize=14)
    ax.grid(alpha=.22); ax.tick_params(labelsize=13)
    ax.legend(fontsize=12, loc='lower right')
    ax.text(.03, .06,
            r'$\int v_{rev}/\int v_{fwd}$ is essentially the reverse/forward flip-count '
            'ratio,\nso it must not be sold as an independent mechanism criterion',
            transform=ax.transAxes, fontsize=11)
    fig.suptitle('Corrected mechanism figure: no same-time scalar in this family '
                 'discriminates pass from fail\n'
                 'replaces fig3_mechanism, which multiplied maxima attained at '
                 'different times', fontsize=15)
    for ext in ('png', 'svg'):
        fig.savefig(FIG / f'fig3_mechanism_v2.{ext}', dpi=150)
    plt.close(fig)
    print('fig3_v2 done')

    # =====================================================================
    # fig 4 v2
    # =====================================================================
    with (r5b / 'eight_states.csv').open(encoding='utf-8') as fh:
        st = list(csv.DictReader(fh))
    st = sorted(st, key=lambda r: int(r['value']))
    seqs = [r['steady_sequence'] for r in st]
    n = max(len(s) for s in seqs)
    M = np.full((len(seqs), n), np.nan)
    for i, s in enumerate(seqs):
        for j, ch in enumerate(s):
            M[i, j] = np.nan if ch == 'x' else int(ch)

    fig, ax = plt.subplots(figsize=(13.5, 7), layout='constrained')
    cmap = ListedColormap(plt.colormaps['viridis'](np.linspace(0, 1, 8)))
    cmap.set_bad('#f2f2f2')
    ax.imshow(M, cmap=cmap, vmin=-.5, vmax=7.5, aspect='auto')
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if np.isnan(M[i, j]):
                ax.text(j, i, 'x', ha='center', va='center', fontsize=11,
                        color=BAD, fontweight='bold')
    ax.set_yticks(range(len(st)))
    ax.set_yticklabels(
        [f"start {r['value']}   t0 = {float(r['t_start_h']):.1f} h   "
         f"snapshot S0/S1/S2 = {float(r['snapshot_S0']):.2f}/"
         f"{float(r['snapshot_S1']):.2f}/{float(r['snapshot_S2']):.2f}"
         for r in st], fontsize=11)
    ax.set_xlabel('steady read window index (after drop-8)', fontsize=15)
    ax.tick_params(labelsize=12)
    cb = fig.colorbar(ax.images[0], ax=ax, ticks=range(8), pad=.02)
    cb.set_label('decoded value', fontsize=14)
    cb.ax.tick_params(labelsize=12)
    allpass = all(r['counting'] == 'True' and r['events'] == 'True' for r in st)
    ax.set_title(
        ('Eight DIGITAL-STATE initial conditions, each continued with its own '
         'upstream time  dy/dt = f(t + t_start, y)\n'
         + ('all eight read a mod-8 sequence from their own start: 8/8 counting, '
            '8/8 event causality' if allpass else 'not all pass'))
        + '\nreplaces fig4_eight_phases, which sampled 9.20 h inside one 10.52 h '
          'clock cycle (all one digital value) and restarted the upstream at t=0',
        fontsize=14)
    for ext in ('png', 'svg'):
        fig.savefig(FIG / f'fig4_eight_states_v2.{ext}', dpi=150)
    plt.close(fig)
    print('fig4_v2 done')

    # =====================================================================
    # status note
    # =====================================================================
    lines = [
        '# 图件状态',
        '',
        '## 当前可用',
        '',
        '| 文件 | 内容 | 依据 |',
        '|---|---|---|',
        '| `fig1_before_after.*` | `clock_K=0.40` 失败 vs `0.10` 修复 | round 2/3，300 h，单初态 |',
        '| `fig2_clock_K_window.*` | `clock_K` 工作区与两种失效模式 | round 2，含锚点复现校验 |',
        '| **`fig3_mechanism_v2.*`** | 同刻通量统计及其**不判别** | round 5 Part C |',
        '| **`fig4_eight_states_v2.*`** | **八数字初态**（1–7、0）各自带上游时间续算，8/8 | round 5b |',
        '| `fig5_bypass.*` | `pulse_gate` 存在/移除对比，`clock_K` 均为 0.40 | round 4b，**300 h、单初态** |',
        '',
        '## 已撤回（保留在盘上，不得引用）',
        '',
        '| 文件 | 撤回原因 |',
        '|---|---|',
        '| `fig3_mechanism.*` | 指标用 `max(I2)·max(RDF2)` —— 两个不同时刻的极值相乘，不是模型到过的反应状态；最多夸大 6.55 倍。由 `fig3_mechanism_v2` 取代 |',
        '| `fig4_eight_phases.*` | 八个"相位"采样跨度仅 9.20 h = 0.875 个时钟周期，**全部落在同一个数字值 2**（phase 7 未定），且续算把上游时间重置为 0。由 `fig4_eight_states_v2` 取代 |',
        '| `fig6_phase_events.*` | 同上，"八相位"构造不成立 |',
        '',
        '## 仍需注明的适用条件',
        '',
        '- `fig5_bypass` 的通过是 **300 h、单初态**；600 h 验证属于 **保留时钟门、`clock_K=0.10`、负自馈 off** 的另一个配置。两者验证范围不同，不得混写。',
        '- 旁路后不同 `clock_K` 结果完全相同，是因为该参数**已不再进入动力学**；这只验证旁路实现正确，**不能**说明时钟门在所有工况下"除了伤害什么都没做"。',
    ]
    (FIG / 'FIGURE_STATUS.md').write_text('\n'.join(lines) + '\n',
                                          encoding='utf-8')
    print('wrote', FIG / 'FIGURE_STATUS.md')


if __name__ == '__main__':
    main()
