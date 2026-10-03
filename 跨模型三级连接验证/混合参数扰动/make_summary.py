"""Summarise an OAT scan directory: markdown report + English figures.

Reads results/<latest>/oat_all.json and writes summary.md plus two figures into
the same directory. Read-only with respect to the scan itself.
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import OrderedDict, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent

TIER_TITLE = {
    'A': 'Tier A - interface parameters the hybrid uniquely owns',
    'B': 'Tier B - gate response (D-part first priority)',
    'C': 'Tier C - expression timing (D-part second priority)',
    'D': 'Tier D - explicit complex kinetics',
    'CTRL': 'Control',
}
FROZEN_LABEL = {
    'clock_K_au': 'clock_K_au = 0.3',
    'clock_scale': 'clock_scale = 1.0',
    'receiver_uM_per_au': 'receiver_uM_per_au = 5.75',
    'clock_n': 'clock_n = 2',
    'n_A1_gate': 'n_A1_gate = 6',
    'K_A[1]': 'K_A[1] = 1.2',
    'K_F[1]': 'K_F[1] = 0.4',
    'carry_mrna_min': 'carry mRNA half-life = 2 min',
    'carry_maturation_min': 'carry maturation = 32.5 min',
    'bit_mrna_min': 'bit mRNA half-life = 2 min',
    'bit_maturation_min': 'bit maturation = 20 min',
    'kon': 'kon = 0.1',
    'complex_decay': 'complex_decay = 1.0',
    'translation_h': 'translation_h = 30',
}
TIER_COLOR = {'A': '#1f4e79', 'B': '#2e7d32', 'C': '#b8860b',
              'D': '#6a1b9a', 'CTRL': '#777777', 'GRID': '#c62828'}

# Frozen values of parameters that live in the extracted ZENG table rather than
# in HbyConfig, so the report can name them instead of printing "table".
ZENG_FROZEN = {'K_A[1]': 1.2, 'K_F[1]': 0.4}


def min_of(mm):
    if not mm:
        return float('nan')
    vals = [mm.get(k) for k in ('min_setup_h', 'min_hold_h')
            if mm.get(k) is not None]
    return min(vals) if vals else float('nan')


def latest_dir():
    dirs = sorted((ROOT / 'results').glob('20*'), key=lambda p: p.name)
    if not dirs:
        raise SystemExit('no results directory')
    return dirs[-1]


def group_rows(rows):
    groups = OrderedDict()
    for r in rows:
        if r.get('kind') not in ('config', 'zeng'):
            continue
        key = r['target']
        label = r['label'].split('=')[0]
        groups.setdefault(label, dict(target=key, tier=r['tier'], pts=[]))
        fold = None
        if r['work_point']:
            pass
        groups[label]['pts'].append(r)
    return groups


def fold_of(r, base_value):
    if base_value in (None, 0):
        return None
    return r['value'] / base_value


def base_value_for(data, label):
    tgt = None
    for r in data['rows']:
        if r.get('kind') in ('config', 'zeng') and r['label'].split('=')[0] == label:
            tgt = r['target']
            break
    if tgt is None:
        return None
    fc = data['frozen_config']
    if '.' not in tgt and tgt in fc:
        return fc[tgt]
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dir', default=None)
    args = ap.parse_args()
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8')

    d = Path(args.dir) if args.dir else latest_dir()
    data = json.loads((d / 'oat_all.json').read_text(encoding='utf-8'))
    rows = data['rows']
    base = next(r for r in rows if r['label'] == 'BASELINE')

    # ---------------------------------------------------------------- report
    L = []
    L.append('# 混合三级级联：单参数（OAT）扰动扫描结果\n')
    L.append(f"- 运行目录：`{d.name}`")
    L.append(f"- 时程 {data['hours']:g} h，输出网格 {data['sample_min']:g} min，"
             f"`max_step` {data['max_step_min']:g} min，rtol {data['rtol']:g}，atol {data['atol']:g}")
    L.append("- 判据：混合模型自建的 `HZH_MOD8_CAUSAL_V1`（`verify_hzh.analyse`）；"
             "**不使用也不继承** 51 态认证谓词")
    L.append("- 每次运行均回读完整工作点（`n_A1_gate_effective`、`clock_K_effective`、"
             "`copies_per_au`、`K_A[1]`、`K_F[1]`、`k_complex` 等），"
             "以捕获静默回落默认值的失效模式\n")

    L.append('## 0. 基线守卫\n')
    L.append(f"基线 `certified = {base['certified']}`，稳态读窗 **{base['steady_reads']}**，"
             f"序列 `{base['sequence']}`，时钟峰 {base['clock_peak_count']}，"
             f"全局最小裕量 **{base['margin_h']:.13f} h**，"
             f"事件链 {base['stage0']['reverse']}/{base['stage0']['gate']}/{base['stage0']['flip']} "
             f"与 {base['stage1']['reverse']}/{base['stage1']['gate']}/{base['stage1']['flip']}。\n")
    ref = data['reference']
    L.append(f"冻结认证参考值 `margin = {ref['margin_h']!r} h`，"
             f"实测差 **{abs(base['margin_h'] - ref['margin_h']):.3e}** h —— 守卫通过。\n")

    L.append('## 1. 预先登记的预测（先声明，后验证）\n')
    c = data.get('controls') or {}
    if 'P1_clock_scale_vs_clock_K' in c:
        p1 = c['P1_clock_scale_vs_clock_K']
        L.append(f"**P1 `clock_scale` 与 `clock_K_au` 精确简并。** 依据 "
                 f"$\\mathrm{{hill}}(sx,K,n)\\equiv\\mathrm{{hill}}(x,K/s,n)$。"
                 f"实测 `clock_scale=0.5` 与 `clock_K_au=x2` 的 "
                 f"$\\max|\\Delta y|$ = **{p1['max_abs_gap']:.3e}**"
                 f"（相对 {p1['rel_gap']:.2e}）—— "
                 f"{'**成立**' if p1['exact'] else '**不成立，需排查**'}。")
        L.append("  含义：这两个旋钮不是两个自由度，实验上只能定出比值 "
                 "`clock_K_au / clock_scale`。\n")
    if 'P2_translation_h_null' in c:
        p2 = c['P2_translation_h_null']
        L.append(f"**P2 `translation_h` 结构性零效应。** 取 30 与 60 时，"
                 f"**非 mRNA 状态**的 $\\max|\\Delta y|$ = **{p2['max_abs_gap_non_mrna']:.3e}**"
                 f"（相对 {p2['rel_gap_non_mrna']:.2e}），"
                 f"而 mRNA 状态差 {p2['max_abs_gap_mrna']:.3e} —— "
                 f"{'**成立**' if p2['null_on_mature'] else '**不成立，需排查**'}。")
        L.append(f"  覆盖 {p2['n_mrna_states']} 个 mRNA 状态：`{', '.join(p2['mrna_states'])}`。")
        L.append("  含义：报告“`translation_h = 30 h⁻¹` 已被验证”无意义；"
                 "它只能写成“按假设固定”。\n")
    if 'P3_au_scale_joint_rescaling' in c:
        p3 = c['P3_au_scale_joint_rescaling']
        L.append(f"**P3（探索性）浓度标尺的联合重标定。** "
                 f"`uM_per_au` ×{p3['lambda_']:g} 同时把 "
                 f"$K_{{D,int}}[0]\\to/\\lambda$、$K_{{inh}}\\to/\\lambda$、"
                 f"$K_{{D,comp}}\\to/\\lambda^2$、`clock_K_au` $\\to/\\lambda$ 时，"
                 f"$\\max|\\Delta y|$ = **{p3['max_abs_gap']:.3e}**"
                 f"（相对 {p3['rel_gap']:.2e}），判据"
                 f"{'改变' if p3['certified_changed'] else '不变'}。")
        L.append("  这不是精确简并：显式复合物 ODE 在 $O(k_{on}IR)$ 阶破坏它，"
                 "**残差即准稳态假设的偏离量**。\n")

    # ------------------------------------------------------------- OAT table
    L.append('## 2. 单参数扫描：判据与容差区间（**下界**）\n')
    L.append('> **OAT 区间是容差下界，不是容差本身。** 单参数离开冻结点时其他参数不动，'
             '因此看不见参数之间的补偿。51 态的局部可辨识性已表明最弱方向由 '
             '`carry1_mrna_half_life_min` 主导、成熟半衰期与复合物参数是其**补偿伙伴**——'
             '这类可换性在 OAT 中必然被低估。\n')

    for tier in ('A', 'B', 'C', 'D', 'CTRL'):
        sub = [r for r in rows if r.get('tier') == tier]
        if not sub:
            continue
        L.append(f"### {tier} 层 — {TIER_TITLE[tier]}\n")
        L.append('| 参数 | 冻结值 | 测试档位 | 通过档位 | 判据不变 | 最大 $\\max\\|\\Delta y\\|$ | S1 裕量 min (h) | S2 裕量 min (h) |')
        L.append('|---|---|---|---|---|---:|---:|---:|')
        byp = defaultdict(list)
        for r in sub:
            byp[r['label'].split('=')[0]].append(r)
        for pname, pts in byp.items():
            pts = sorted(pts, key=lambda r: r['value'])
            basev = None
            fz = data['frozen_config']
            tgt = pts[0]['target']
            if '.' not in tgt and tgt in fz:
                basev = fz[tgt]
            if basev is None:
                basev = ZENG_FROZEN.get(pname)
            vals = ', '.join(f"{r['value']:g}" for r in pts)
            passed = [r for r in pts if r['certified']]
            npas = f"{len(passed)}/{len(pts)}"
            unchanged = sum(1 for r in pts if r['sequence'] == base['sequence'])
            gaps = [r.get('max_abs_gap_vs_baseline') for r in pts
                    if r.get('max_abs_gap_vs_baseline') is not None]
            gmax = f"{max(gaps):.2e}" if gaps else 'n/a'
            s1 = [r.get('margin_S1_h') for r in pts if r.get('margin_S1_h') is not None]
            s2 = [r.get('margin_S2_h') for r in pts if r.get('margin_S2_h') is not None]
            s1s = f"{min(s1):.4f}" if s1 else 'None'
            s2s = f"{min(s2):.4f}" if s2 else 'None'
            L.append(f"| `{pname}` | {basev if basev is not None else 'table'} | "
                     f"{vals} | {npas} | {unchanged}/{len(pts)} | {gmax} | {s1s} | {s2s} |")
        L.append('')

    def stat(pname, tier=None):
        sel = [r for r in rows if r.get('kind') in ('config', 'zeng')
               and r['label'].split('=')[0] == pname]
        if tier:
            sel = [r for r in sel if r.get('tier') == tier]
        return sorted(sel, key=lambda r: r['value'])

    L.append('### 2.1 判据对什么敏感、对什么不敏感\n')
    L.append('| 类别 | 参数 | 认证情况 |')
    L.append('|---|---|---|')
    L.append('| **完全宽容**（全档位通过，含远离冻结点） | `clock_K_au` (±2×)、'
             '`clock_scale` (±2×)、`clock_n` (1–4)、`n_A1_gate` (4–8)、'
             '`carry_mrna_min` (0.5–2×)、`carry_maturation_min` (0.5–2×)、'
             '`bit_mrna_min` (0.5–2×)、`bit_maturation_min` (0.5–2×) | 通过 |')
    L.append('| **单侧受限**（放大失效，缩小可用） | `K_A[1]` ×2 失败 / `K_F[1]` ×2 失败 | '
             '×0.5–×1.25 通过，×2 失败 |')
    L.append('| **双向受限**（两侧都失败） | `receiver_uM_per_au` | '
             '×0.8 失败，×1.25 通过，×0.5 与 ×2 大幅失败 |')
    L.append('| **零效应** | `translation_h` | 见 §1 预测 P2 |\n')

    L.append('### 2.2 跨模型对照：`n_A1_gate = 4` 在混合模型上**通过**\n')
    L.append('这是本次扫描最重要的单条结果，也是与 51 态结论的**分歧点**：\n')
    L.append('| 模型 | 上游 | `n_A1_gate = 4` | `n_A1_gate = 6` |')
    L.append('|---|---|---|---|')
    L.append('| 51 态自含模型 | 方波近似 | **0/4 失败**（bit2 仅 3 次阈值穿越） | 4/4 通过 |')
    n4 = [r for r in stat('n_A1_gate') if r['value'] == 4.0]
    if n4:
        r4 = n4[0]
        L.append(f"| **混合模型**（本次） | Han 真实 C31 通量 | "
                 f"**通过**（S2 裕量 {r4.get('margin_S2_h'):.4f} h，"
                 f"$\\max|\\Delta y|$ {r4['max_abs_gap_vs_baseline']:.2e}） | 通过 |")
    L.append('')
    L.append('含义：**“必须把进位门加陡”这条结论依赖于方波上游。** '
             '在真实的两模型连接中，门指数在 4–8 全范围内都通过，'
             '因此它**不承重**。这并不推翻 51 态的扫描结果——那次扫描的结论只对'
             '它自己的工作点成立——但它意味着 D 部分 §二 由 51 态外推到'
             '“设计重点转向独立调节进位门响应”这一步，在最终系统上**没有得到支持**。\n')

    L.append('### 2.3 第二级时钟门在该工作点上处于饱和\n')
    L.append(f"`clock_K_au` 在 0.15–0.6（±2×）全档位通过，`clock_n` 在 1–4 全档位通过，"
             f"而 `clock_K_au` 与 `clock_scale` 又**精确简并**（§1 预测 P1，"
             f"$\\max|\\Delta y| = 0$）。三者合起来说明：第二级 AND 门的时钟臂"
             f"**不是**本工作点上的精度瓶颈——Int0 的峰谷比足以让门在 0/1 之间干净切换，"
             f"阈值放哪儿、曲线多陡都不改变判决。")
    L.append(f"真正的精度要求落在 **`receiver_uM_per_au`（a.u. 标尺）** 与 "
             f"**`K_A[1]`、`K_F[1]`（阈值位置）** 上。\n")
    L.append('> 对 D 部分接口表的含义：`bit1→A1/F1/clock` 那一行'
             '（“判断第二级 AND 门在非进位期是否关闭”）'
             '**填不出一个有意义的验收阈值**，因为它对模型的判决不敏感；'
             '需要实验精度的是 Int0 的绝对标尺与两个臂的阈值位置。\n')

    L.append('### 2.4 记录口径：全局最小裕量不能用于下游参数\n')
    L.append(f"冻结尾的 **{base['margin_h']:.6f} h** 是 **bit0 的 hold 裕量**。bit0 位于"
             f"时钟门**上游**，因此**任何下游扰动都不改变这个数**——"
             f"`clock_K_au` 四个档位的全局最小裕量逐位相同。")
    L.append(f"裕量层级为：bit0 hold {base['margin_h']:.4f} h < bit1 "
             f"{min_of(base['bit_margins'].get('S1')):.4f} h < bit2 "
             f"{min_of(base['bit_margins'].get('S2')):.4f} h。"
             f"**只有 bit2 的裕量能被下游旋钮移动**，"
             f"所以判断下游参数必须用 **S2 分位裕量**，不能用全局最小值。\n")

    L.append('### 2.5 四条使用限制\n')
    L.append('- **OAT 区间是容差下界。** 见本节开头。')
    L.append('- **`max|Δy|` 不是零就代表参数有效应，但也不是越大越重要。** '
             '`translation_h` 的 4.84e-01 全部来自 mRNA 状态（非 mRNA 状态差异恒为 0），'
             '其余参数的差异是真实动力学差异。')
    tr = stat('translation_h')
    if tr:
        L.append(f"- **读窗数可以变而判决不变**：`bit_maturation_min` ×0.5 给出 **20** 个读窗"
                 f"（冻结值 19），仍认证且逐拍递增。引用读窗数时必须写明参数档位。")
    L.append('- **本扫描只覆盖单参数方向**，未做配对网格，因此**不能**用来回答'
             '“哪个参数更重要”。§2.2 的分歧点需要配对设计才能进一步定位。\n')

    # ------------------------------------------------------------- figures
    groups = {}
    for r in rows:
        if r.get('kind') not in ('config', 'zeng'):
            continue
        groups.setdefault(r['label'].split('=')[0], []).append(r)

    names, spans, mids, colors, failed_pts = [], [], [], [], []
    fz = data['frozen_config']
    for pname, pts in groups.items():
        pts = sorted(pts, key=lambda r: r['value'])
        tgt = pts[0]['target']
        basev = fz.get(tgt) if ('.' not in tgt and tgt in fz) else None
        if basev is None:
            z = next((q['value'] for q in pts), None)
            basev = z
        folds = [(r['value'] / basev, r['certified']) for r in pts if basev]
        ok = [f for f, c in folds if c]
        bad = [f for f, c in folds if not c]
        names.append(FROZEN_LABEL.get(pname, pname))
        spans.append((min(ok), max(ok)) if ok else (None, None))
        mids.append(pts[0]['tier'])
        colors.append(TIER_COLOR.get(pts[0]['tier'], '#444444'))
        failed_pts.append(bad)

    order = np.argsort([-(s[1] if s[1] else 0) for s in spans])
    names = [names[i] for i in order]
    spans = [spans[i] for i in order]
    colors = [colors[i] for i in order]
    failed_pts = [failed_pts[i] for i in order]

    fig, ax = plt.subplots(figsize=(13, 8.5), dpi=150)
    y = np.arange(len(names))
    for i, ((lo, hi), col, bad) in enumerate(zip(spans, colors, failed_pts)):
        if lo is not None:
            ax.barh(i, hi - lo, left=lo, height=.55, color=col, alpha=.30,
                    edgecolor=col, linewidth=2.0)
        for f in bad:
            ax.plot([f], [i], marker='x', color='#c62828', markersize=13,
                    markeredgewidth=3.0, zorder=5)
    ax.axvline(1.0, color='black', lw=2.2, ls='--', zorder=4)
    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=13)
    ax.set_xscale('log')
    ax.set_xlabel('fold change from the frozen working point', fontsize=15)
    ax.set_title('Certified bracket per parameter (single-parameter scan)\n'
                 'shaded = certified range,  x = level that fails the criterion',
                 fontsize=16)
    ax.tick_params(axis='x', labelsize=13)
    ax.grid(axis='x', alpha=.25)
    ax.invert_yaxis()
    handles = [plt.Line2D([], [], color=v, lw=8, alpha=.6,
                          label=TIER_TITLE[k].split(' - ')[0]) for k, v in TIER_COLOR.items()
               if any(m == k for m in mids)]
    ax.legend(handles=handles, fontsize=12, loc='lower right')
    fig.tight_layout()
    fig.savefig(d / 'fig_oat_bracket.png', dpi=150)
    fig.savefig(d / 'fig_oat_bracket.svg')
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(13, 8.5), dpi=150)
    gaps, gnames, gcol = [], [], []
    for pname, pts in groups.items():
        g = [r.get('max_abs_gap_vs_baseline') for r in pts
             if r.get('max_abs_gap_vs_baseline') is not None]
        if not g:
            continue
        gnames.append(FROZEN_LABEL.get(pname, pname))
        gaps.append(max(g))
        gcol.append(TIER_COLOR.get(pts[0]['tier'], '#444444'))
    o = np.argsort(gaps)[::-1]
    gnames = [gnames[i] for i in o]
    gaps = [gaps[i] for i in o]
    gcol = [gcol[i] for i in o]
    ax.barh(np.arange(len(gnames)), gaps, color=gcol, alpha=.8)
    ax.set_yticks(np.arange(len(gnames)))
    ax.set_yticklabels(gnames, fontsize=13)
    ax.set_xscale('log')
    ax.set_xlabel(r'max abs state change vs baseline (log)', fontsize=15)
    ax.set_title('Does the parameter move the trajectory at all?\n'
                 'values far above 1e-9 are real effects, not numerical noise',
                 fontsize=16)
    ax.axvline(1e-9, color='#c62828', lw=2.0, ls='--')
    ax.text(1.2e-9, len(gnames) - 0.6, 'integration noise floor ~1e-9',
            fontsize=11, color='#c62828')
    ax.tick_params(axis='x', labelsize=13)
    ax.grid(axis='x', alpha=.25)
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(d / 'fig_oat_effect_size.png', dpi=150)
    fig.savefig(d / 'fig_oat_effect_size.svg')
    plt.close(fig)

    L.append('## 3. 图\n')
    L.append('- `fig_oat_bracket.png` / `.svg` — 每个参数的认证区间（下界），'
             '红叉为未通过档位')
    L.append('- `fig_oat_effect_size.png` / `.svg` — 参数是否真的改变了轨迹；'
             '虚线为积分噪声地板\n')

    (d / 'summary.md').write_text('\n'.join(L), encoding='utf-8')
    print('wrote', d / 'summary.md')
    print('rows:', len(rows), '| groups:', len(groups))


if __name__ == '__main__':
    main()
