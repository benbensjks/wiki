"""Combine the three stages into one report plus the paired-grid figure.

v2: two attributions from v1 are retracted (see the 修正记录 section), the
criterion is reported as three separate tiers, and every "range" statement is
stated as a set of tested discrete levels.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
RES = HERE / 'results'

ZENG_FROZEN = {'K_A[1]': 1.2, 'K_F[1]': 0.4}
# Parameters that act on expression timing only (Tier C).
EXPRESSION_PARAMS = ['carry_mrna_min', 'carry_maturation_min',
                     'bit_mrna_min', 'bit_maturation_min']


def newest(pat):
    d = sorted(RES.glob(pat), key=lambda p: p.name)
    if not d:
        raise SystemExit(f'no results for {pat}')
    return d[-1]


def load(p):
    return json.loads(p.read_text(encoding='utf-8'))


def tiers(r):
    """The criterion decomposes exactly as certified = counting AND events."""
    s0 = r.get('stage0') or {}
    s1 = r.get('stage1') or {}
    counting = bool(r.get('increments_mod8'))
    events = bool(s0 and s1
                  and s0.get('one_to_one') and s0.get('order') and s0.get('alt')
                  and s1.get('one_to_one') and s1.get('order') and s1.get('alt'))
    return counting, events, bool(r.get('certified'))


def main():
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8')
    s1d, s2d, s3d = newest('2026*'), newest('refine_*'), newest('threshold_*')
    s1, s2, s3 = load(s1d / 'oat_all.json'), load(s2d / 'refine_all.json'), load(s3d / 'threshold_sweep.json')
    base = next(r for r in s1['rows'] if r['label'] == 'BASELINE')

    L = []
    A = L.append
    A('# 混合三级级联参数扰动分析：结论汇总')
    A('')
    A('> 本目录为一个**新增的独立实验目录**，未修改 `hby_zmh_hby/`、`hybrid_model.py`、')
    A('> `model_hzh.py`、`verify_hzh.py`、`run_han_comparison.py` 或任何既有主文件。')
    A('')
    A('## 0. 口径与可复现性')
    A('')
    A('- 判据：混合模型自建的 `HZH_MOD8_CAUSAL_V1`（`verify_hzh.analyse`）；'
      '**不使用也不继承** 51 态认证谓词')
    A(f'- 时程 {s1["hours"]:g} h，输出网格 {s1["sample_min"]:g} min，'
      f'`max_step` {s1["max_step_min"]:g} min，rtol {s1["rtol"]:g}，atol {s1["atol"]:g}')
    A(f'- 运行目录：`{s1d.name}`（OAT）、`{s2d.name}`（边界夹逼 + 配对网格）、'
      f'`{s3d.name}`（阈值扫描，零次重积分）')
    A('- 每次运行均回读完整工作点（`n_A1_gate_effective`、`clock_K_effective`、'
      '`copies_per_au`、`K_A[1]`、`K_F[1]`、`k_complex`），以捕获静默回落默认值')
    A('')
    A(f'**基线守卫通过**：稳态读窗 {base["steady_reads"]}，序列 `{base["sequence"]}`，'
      f'全局最小裕量 {base["margin_h"]:.13f} h，事件链 14/14/14 与 7/7/7；'
      f'与冻结认证值（{s1["reference"]["margin_h"]!r} h）相差 '
      f'{abs(base["margin_h"] - s1["reference"]["margin_h"]):.2e} h。')
    A('')
    A('### 0.1 判据按三层分别报告')
    A('')
    A('`HZH_MOD8_CAUSAL_V1` 的判决在代码里就是两部分的与：'
      '`certified_v1 = steady.passed and all(event.passed)`。'
      '因此本报告一律分开给出三层，**不得把事件层失败写成“电路没有计数”**：')
    A('')
    A('| 层 | 定义 | 字段 |')
    A('|---|---|---|')
    A('| `counting_passed` | 稳态读窗的码值逐拍模 8 递增 | `steady.increments_mod8` |')
    A('| `event_causality_passed` | 在**预先声明的事件阈值**下，两条事件链一一对应、'
      '因果有序、方向交替 | `events.*.{one_to_one,causal_order,alternating_directions}` |')
    A('| `full_certified` | 两者同时通过 | `certified_v1` |')
    A('')
    A('### 0.2 采样口径')
    A('')
    A('**本报告所有“范围”一律是离散采样结果，不是连续区间。** '
      '档位之间为线性/对数插值假设，**未验证**。凡出现“通过/失败”均指被扫过的那个具体档位。')
    A('')

    # ---------------------------------------------------------- predictions
    c = s1['controls']
    p1, p2 = c['P1_clock_scale_vs_clock_K'], c['P2_translation_h_null']
    A('## 1. 预先登记的预测（先声明，后验证）')
    A('')
    A('只有下面两条是**代数恒等式**级别的预登记预测，写在扫描之前，且可被证伪。')
    A('')
    A(f"**P1 `clock_scale` 与 `clock_K_au` 精确简并 —— 成立（逐位精确）。** "
      f"依据 $\\mathrm{{hill}}(sx,K,n)\\equiv\\mathrm{{hill}}(x,K/s,n)$。"
      f"`clock_scale=0.5` 与 `clock_K_au=×2` 的 $\\max|\\Delta y|$ = "
      f"**{p1['max_abs_gap']:.3e}**（相对 {p1['rel_gap']:.2e}）。"
      f"OAT 表独立重现了整条简并：两者在每一个对应档位上的 $\\max|\\Delta y|$、"
      f"S2 裕量、事件链**全部逐位相同**。")
    A('')
    A('  → 这两个旋钮在模型中**不是两个自由度**；实验上只能定出比值 '
      '`clock_K_au / clock_scale`。把二者分别标定并相加误差的写法不成立。')
    A('')
    A(f"**P2 `translation_h` 结构性零效应 —— 成立（逐位精确）。** "
      f"取 30 与 60 时，**非 mRNA 状态**的 $\\max|\\Delta y|$ = "
      f"**{p2['max_abs_gap_non_mrna']:.3e}**，而 mRNA 状态差 "
      f"{p2['max_abs_gap_mrna']:.3e}（覆盖 {p2['n_mrna_states']} 个状态："
      f"`{', '.join(p2['mrna_states'])}`）。")
    A('')
    A('  → 报告“`translation_h = 30` 已被验证”无意义；只能写成“按假设固定”。'
      '注意该参数的 $\\max|\\Delta y|$ 读数（4.84e-01）**不是零**，'
      '它**全部来自 mRNA 状态**。')
    A('')

    A('### 1.1 探索性检查（**非预登记预测**，结论已收窄）')
    A('')
    A('扫描脚本里这条检查带的是 `exploratory=True`，**不是**与 P1/P2 并列的预登记预测。')
    A('')
    A(f"**尝试的一组局部联合重标定未能保持轨迹不变。** 做法是把 `uM_per_au` ×2，"
      f"同时把 $K_{{D,int}}[0]\\to/\\lambda$、$K_{{inh}}\\to/\\lambda$、"
      f"$K_{{D,comp}}\\to/\\lambda^2$、`clock_K_au` $\\to/\\lambda$，"
      f"测得 $\\max|\\Delta y|$ = **{c['P3_au_scale_joint_rescaling']['max_abs_gap']:.3e}**。")
    A('')
    A('> **该变换不是当前方程的一组严格候选对称性**，因此上面的残差只能说明'
      '“这套具体重标定规则不等价”。**不能**由此推出：')
    A('>')
    A('> - 准稳态假设失效；')
    A('> - a.u. 标尺原则上无法被其他参数吸收；')
    A('> - 根因是复合物弛豫率与 Int 周转率同量级。')
    A('')
    A('原因是 `uM_per_au` 只经 `copies_per_au` 改变 bit0 的 **C31 mRNA 源**'
      '（`model_hzh.py:49,81`），而 bit0 的 RDF 由 '
      '`alpha_rdf*s*(1-hill(tr,K_rep,n_rep))` 驱动 —— 取决于 $S$ 与 $T$，'
      '**与 C31 通量无关**。所以 $R$ 并不随该标尺缩放，'
      '脚本却按 $1/\\lambda$ 缩放了 $K_{inh}$、按 $1/\\lambda^2$ 缩放了 $K_{D,comp}$，'
      '这两步都缺少依据；`bind = kon·i·r` 也应按 $1/\\lambda$ 而非 $1/\\lambda^2$ 变化。')
    A('')
    A('（$k_{off}+\\delta_C = 2.0\\ \\mathrm{h^{-1}}$ 与 '
      '$\\gamma_{\\mathrm{int}}[0] = 2.0\\ \\mathrm{h^{-1}}$ 同量级是一句'
      '**独立的算术事实**，但它**不是**上面残差的解释，不得连写。）')
    A('')

    # ------------------------------------------------------------ paired grid
    grid = [r for r in s2['rows'] if r.get('kind') == 'grid']
    gates = sorted({r['n_A1_gate'] for r in grid})
    folds = sorted({r['uM_fold'] for r in grid})

    def cell(g, f):
        for r in grid:
            if r['n_A1_gate'] == g and r['uM_fold'] == f:
                return r
        return None

    A('## 2. 核心结果：配对网格')
    A('')
    A('这是唯一能回答“两个旋钮能否互相补偿”的设计'
      '（单参数扫描在构造上看不见可换性）。网格为 '
      f'`n_A1_gate` ∈ {{{", ".join(f"{g:g}" for g in gates)}}} × '
      f'`uM_per_au` ∈ {{{", ".join("×"+format(f, "g") for f in folds)}}}。')
    A('')
    A('| `n_A1_gate` \\ `uM_per_au` | ' + ' | '.join(f'×{f:g}' for f in folds) + ' |')
    A('|---|' + '---|' * len(folds))
    for g in gates:
        cells = []
        for f in folds:
            r = cell(g, f)
            if r is None:
                cells.append('**通过**（冻结工作点）')
            else:
                cnt, ev, full = tiers(r)
                mark = '通过' if full else ('事件层失败' if cnt else '计数失败')
                cells.append(mark)
        A(f'| **n = {g:g}** | ' + ' | '.join(cells) + ' |')
    A('')
    A('**三行完全相同。** 在 $n_{A1,gate} = 4, 6, 8$ 三个被扫点上，'
      '`uM_per_au` 的通过/失败格局是**同一个**：'
      '×1.0 与 ×1.25 通过，×0.8 与 ×2 失败。')
    A('')
    A('> **结论（限定在已扫描范围内）**：在当前混合架构和已扫描范围内，'
      '将第二级门指数从 4 提高到 6 或 8，**没有扩大浓度标尺的可计数区间**；'
      '门指数与浓度标尺在这组网格上**不是可以相互补偿的两个旋钮**。')
    A('>')
    A('> 这**不等于**“任何情况下加陡都无效”——它只是这组离散网格的结果。'
      '也**不推翻**独立 51 状态模型中选择 `n_A1_gate = 6` 的结果：')
    A('>')
    A('> - **独立 51 状态模型**：`n=6` 是该模型泄漏、裕量与协同性之间的工程折中；')
    A('> - **混合模型**：已测的 `n=3–8` 都能计数，`n=6` 可以保留为兼容工作点，'
      '但**不再是混合模型成功的承重条件**。')
    A('')
    A('两端失效的机制不同：')
    A('')
    A(f"- **×2 端**：`uM_per_au` 变大 ⇒ `copies_per_au` 变大 ⇒ bit0 的 C31 mRNA 源"
      f"按 a.u. 变小，时钟门峰值从 0.770 掉到 **0.455**，第二级门 $g_1$ 未达事件阈值"
      f"（stage1 事件 0/0/0），**计数失败**。")
    A(f"- **×0.8 端**：时钟峰值升到 **0.839**，同一个反向通量波形被固定阈值"
      f'切成 **28 段**，而 bit1 实际翻转仍是 **14** 次、码串仍严格模 8 —— '
      f'**这是事件分割假象，不是过度触发，也不是计数失败**（见 §3）。')
    A('')

    # --------------------------------------------------------------- tiers
    A('## 3. 计数层与事件层必须分开：哪些“失败”不是计数失败')
    A('')
    A('| 案例 | 码串 | `counting` | `event` | `full` | 性质 |')
    A('|---|---|---|---|---|---|')
    cases = [
        ('uM ×0.8', 'receiver_uM_per_au=x0.8'),
        ('clock_K_au ×4', 'clock_K_au=x4'),
        ('uM ×1.5', 'receiver_uM_per_au=x1.5'),
        ('uM ×1.75', 'receiver_uM_per_au=x1.75'),
        ('K_A[1] ×1.75', 'K_A[1]=x1.75'),
        ('n_A1_gate = 2', 'n_A1_gate=2'),
    ]
    pool = s1['rows'] + [r for r in s2['rows'] if r.get('kind') != 'grid']
    for lab, key in cases:
        r = next((x for x in pool if x['label'] == key), None)
        if r is None:
            raise SystemExit(f'summary case not found in results: {key!r}')
        cnt, ev, full = tiers(r)
        note = {
            'receiver_uM_per_au=x0.8': '计数正确，反向通量被切成 28 段（分割假象）',
            'clock_K_au=x4': '计数正确，门峰值未过默认 0.05 事件阈值（分割假象）',
            'receiver_uM_per_au=x1.5': '**真失效**：第二级门从不点火',
            'receiver_uM_per_au=x1.75': '**真失效**',
            'K_A[1]=x1.75': '**真失效：退化成模 4**',
            'n_A1_gate=2': '事件链全对，**位状态读不出**（信号层）',
        }[key]
        A(f"| {lab} | `{r['sequence']}` | {'✓' if cnt else '✗'} | "
          f"{'✓' if ev else '✗'} | {'✓' if full else '✗'} | {note} |")
    A('')
    A(f"**验证（`{s3d.name}`）**：对每条轨迹**只积分一次**，然后固定轨迹、"
      f"只改两个分段阈值（`JREV_THRESHOLD_H`、`GATE_THRESHOLD`），"
      f"每案例扫 {len(s3['jrev_sweep'])}×{len(s3['gate_sweep'])} 个设定"
      f"（去重后 {len(s3['cases']['uM_per_au=x0.8']['sweeps'])} 个）：")
    A('')
    A('| 案例 | 仍计数 | 变为 `full_certified` | 判定 |')
    A('|---|---|---|---|')
    for k, lab in (('uM_per_au=x0.8', '`uM ×0.8`'),
                   ('clock_K_au=x4', '`clock_K_au ×4`'),
                   ('uM_per_au=x1.5', '`uM ×1.5`（对照）')):
        rec = s3['cases'][k]
        n = len(rec['sweeps'])
        ncnt = sum(1 for s in rec['sweeps'] if s['counts'])
        ncert = sum(1 for s in rec['sweeps'] if s['certified'])
        verd = ('**真失效 —— 改阈值也不能恢复计数** ✓' if ncert == 0 and ncnt == 0
                else '**事件分割假象 —— 阈值一改即恢复**')
        A(f'| {lab} | {ncnt}/{n} | {ncert}/{n} | {verd} |')
    A('')
    ex = next(s for s in s3['cases']['uM_per_au=x0.8']['sweeps'] if s['certified'])
    A(f"例：`uM ×0.8` 在 `jrev={ex['jrev']}`、`gate={ex['gate']}` 下，"
      f"stage0 恢复 14/14/14、stage1 7/7/7、码串 `{ex['sequence']}`、"
      f"`full_certified` 通过。")
    A('')
    A('> **引用要求**：默认事件阈值没有识别出事件，**不能**写成“电路没有计数”。'
      '报告任何容差结论时必须同时给出上面三层。§3 的阈值扫描只说明**判据对阈值敏感**，'
      '**不说明哪个阈值“正确”**；冻结认证使用的是 `JREV=0.1`、`GATE=0.05`，'
      '本分析没有也不应改写该认证，只是标注其敏感度。')
    A('')

    # ----------------------------------------------------------- discrete pts
    A('## 4. 各参数的已测离散点结果')
    A('')
    A('下表为**离散采样结果**，不是连续区间。每格列出被测试的具体档位。')
    A('')
    A('| 参数 | 已测通过的档位 | 已测失败的档位 | 正确解释 |')
    A('|---|---|---|---|')
    def pts_of(name):
        sel = [r for r in pool if r.get('kind') in ('config', 'zeng')
               and r['label'].split('=')[0] == name]
        return sorted(sel, key=lambda r: r['value'])
    def fmt(vals):
        return ', '.join(f'{v:g}' for v in vals) if vals else '—'
    rows_spec = [
        ('receiver_uM_per_au', '×0.8、0.9、1、1.1、1.25', '×0.5、1.5、1.75、2',
         '唯一明显呈双侧边界的接口参数'),
        ('K_A[1]', '×0.5、0.8、1、1.25', '×1.5、1.75、2',
         '上侧边界位于 ×1.25 与 ×1.5 之间；单侧上限'),
        ('K_F[1]', '×0.5、0.8、1、1.25', '×1.5、1.75、2', '同样为单侧上限'),
        ('n_A1_gate', '已测 3–8', '2', '低协同性会失败；3 以上本轮没有区分度'),
        ('clock_K_au', '所测档位计数均正确', '×4 仅默认事件判据失败',
         '计数本身不敏感；事件识别对阈值敏感'),
    ]
    for name, ok, bad, note in rows_spec:
        A(f'| `{name}` | {ok} | {bad} | {note} |')
    A(f"| 表达时间参数（{len(EXPRESSION_PARAMS)} 个："
      + '、'.join(f'`{p}`' for p in EXPRESSION_PARAMS) + '） | '
      '单参数 ×0.5–2 各档位均通过 | 未观察到失败 | '
      '仅说明**单参数扫描在这些点上不紧**，不代表联合扰动或实验鲁棒区间 |')
    A('| `translation_h` | 全档位通过 | — | 结构性零效应（预测 P2） |')
    A('| `kon`、`complex_decay` | 各 1 档通过 | — | 采样太少，不作推断 |')
    A('')
    A('### 4.1 三层判据下的完整清单')
    A('')
    A('| 参数 | 档位 | `counting` | `event` | `full` |')
    A('|---|---|---|---|---|')
    groups = {}
    for r in pool:
        if r.get('kind') in ('config', 'zeng'):
            groups.setdefault(r['label'].split('=')[0], []).append(r)
    for name, pts in groups.items():
        pts = sorted(pts, key=lambda r: r['value'])
        vals = ', '.join(f'{r["value"]:g}' for r in pts)
        ncnt = sum(1 for r in pts if tiers(r)[0])
        nev = sum(1 for r in pts if tiers(r)[1])
        nfull = sum(1 for r in pts if tiers(r)[2])
        A(f'| `{name}` | {vals} | {ncnt}/{len(pts)} | {nev}/{len(pts)} | '
          f'{nfull}/{len(pts)} |')
    A('')
    A('### 4.2 敏感性分级（限已测离散点）')
    A('')
    A('| 分级 | 参数 | 已测依据 |')
    A('|---|---|---|')
    A('| **已测点全部通过，未见边界** | `clock_scale`（与 `clock_K_au` 精确简并）、'
      '`clock_n`、`carry_mrna_min`、`carry_maturation_min`、`bit_mrna_min`、'
      '`bit_maturation_min` | 见 §4.1 |')
    A('| **仅低端失败** | `n_A1_gate`（2 失败，3–8 通过） | 见 §4.1 |')
    A('| **仅高端失败** | `K_A[1]`、`K_F[1]` | ×1.25 通过，×1.5 起失败 |')
    A('| **双侧失败** | `receiver_uM_per_au` | ×0.5 与 ×1.5 起失败 |')
    A('| **零效应** | `translation_h` | 逐位精确（P2） |')
    A('')

    # ---------------------------------------------------------- engineering
    A('## 5. 对 Engineering D 的修改建议')
    A('')
    A('### 5.1 设计优先级（替换 D 部分 §二的优先级列表）')
    A('')
    A('1. 先标定接收端浓度尺度 `uM_per_au`；')
    A('2. 再控制 `K_A[1]`、`K_F[1]` 所代表的激活/抑制阈值；')
    A('3. 确保 `n_A1_gate` 不落入低协同性失败区（已测为 `n < 3`）；')
    A('4. 表达时间参数用于动态匹配，而不是当前最主要的修复旋钮。')
    A('')
    A('### 5.2 替换原句（“与其寻找苛刻的时间窗口，不如提升反馈灵敏性”）')
    A('')
    A('> 独立 51 状态模型表明，提高第二级门的有效协同性可以修复该模型中的持续泄漏；'
      '但在最终混合架构中，`n_A1_gate = 3–8` 的已测点均能连续计数，'
      '提高到 6 或 8 没有扩大浓度标尺的成功区间。'
      '因此，门协同性**不应被视为适用于所有架构的通用修复手段**。'
      '对混合模型而言，更直接的工程任务是标定接收端浓度尺度，'
      '并使 A1/F1 调控阈值落入已观察到的工作区；'
      '表达和成熟时间则用于保证写入、关断与恢复能够在一个时钟周期内完成。')
    A('')
    A('### 5.3 需要改写的原有表述')
    A('')
    A('- **“把 `n_A1_gate` 从 5 提到 6 使泄漏降低约 68%”** —— 该 20 点扫描的结论'
      '只对独立 51 状态模型成立。混合模型上必须另行表述为 §2 的结论。')
    A('- **§三 接口表 `bit1→A1/F1/clock` 一行**填不出有意义的验收阈值：'
      '`clock_K_au` 在所测档位上计数均正确，且与 `clock_scale` 精确简并，'
      '该行只能要求“记录比值”，不能要求“达到某个阈值”。')
    A('- **§一 把 `n=6` 作为参考陡度**可以保留，但应补一句它**不承重**，'
      '且 D 部分当前把 `n` 列为首要杠杆的写法需要按 §5.1 调整顺序。')
    A('')
    A('### 5.4 新增实验要求')
    A('')
    A('**须分开表征“进位脉冲的峰值”与“时钟门的开度”，并联合判断。** '
      '本次扫描的失效来自两者不匹配：时钟门开度不足（×2 端）导致第二级完全不动作；'
      '而开度足够时又会让同一个反向通量波形被固定阈值切成多段（×0.8 端）。'
      '这对应 D 部分 §三 接口表中 `上游→bit0` 与 `Int1→bit1` 两行的实测项，'
      '但**不能各自独立验收**。')
    A('')

    # -------------------------------------------------------------- boundary
    A('## 6. 本分析的边界（不得越界引用）')
    A('')
    A('1. **所有“范围”都是离散采样结果。** 档位之间的插值假设未验证。')
    A('2. **OAT 是容差下界，不是容差本身。** 单参数离开冻结点时其余不动，'
      '看不见补偿；§2 的配对网格只覆盖了一对参数。')
    A('3. **未做随机性、参数分布、质粒丢失、资源竞争或实体串扰。** '
      '全部为确定性单轨迹结果；`full_certified` 是布尔值，**不是“成功率”**。')
    A('4. **表达时间参数“已测点全过”不等于鲁棒区间。** 需联合扰动、'
      '参数分布与多初态边界扫描后才能谈鲁棒性。')
    A('5. **§3 只说明判据对事件阈值敏感，不说明哪个阈值正确。** '
      '冻结认证使用 `JREV=0.1`、`GATE=0.05`，本分析未改写该认证。')
    A('6. **`n_A1_gate ≥ 3 即通过`不表示陡度在实验中不重要。** '
      '模型里 $n$ 是**单个 Hill 项的指数**；实体实现能否达到该有效陡度仍是开放项。')
    A('7. **`n=4` 的跨模型差异，机理未定位。** 独立 51 状态模型与混合模型'
      '**使用同一个韩亚轩上游**（`final_reconstruction/model.py:21` 与 '
      '`run_han_comparison.py:26` 指向同一文件、同一 v53d；前者把上游状态内联积分，'
      '后者预解后插值）。两者真正的主要差别在**中间一级**：'
      '前者用自己的展开式 A0/F0、bit1 与旧参数表（`ZENG`），'
      '后者用曾同学现行 Python 参数的约化 A0/F0 与 bit1，两边再连接同一个 bit2 接收端。'
      '因此目前只能说：**对门协同性的要求依赖于中间 bit1 模块产生的输入波形与接口动力学。** '
      '要定位原因，需在保持同一上游与同一 bit2 的条件下只替换中间 bit1 模块。')
    A('8. **本文档 v1 的两条归因已撤回**，见 §7。')
    A('9. **未修改任何既有文件。** 全部新增内容位于 `混合参数扰动/`。')
    A('')

    A('## 7. 修正记录（本版相对 v1）')
    A('')
    A('| # | v1 的写法 | 状态 | 依据 |')
    A('|---|---|---|---|')
    A('| 1 | “51 态上游是方波近似，故 n=4 失败依赖方波上游” | **撤回** | '
      '`final_reconstruction/model.py:21` 与 `run_han_comparison.py:26` '
      '指向同一个 `Flux_Driven_Translation_Burden_Model.py`（v53d）；'
      '`model.py` 中 `square`/`amplitude` 检索零命中。差别在中间一级，不在上游 |')
    A('| 2 | “P3 残差 5.698 证明准稳态假设失效 / a.u. 标尺不可被吸收” | **撤回** | '
      '该联合重标定不是方程的严格候选对称性：`uM_per_au` 只改 bit0 的 C31 mRNA 源，'
      '而 bit0 的 RDF 由 $S$、$T$ 驱动，$R$ 不随该标尺缩放，'
      '故缩放 $K_{inh}$、$K_{D,comp}$ 缺少依据 |')
    A('| 3 | 把 P3 列为“第三个预登记预测” | **改口径** | '
      '脚本产物中该条为 `exploratory=True`；已移至 §1.1 并标注非预登记 |')
    A('| 4 | “`uM ×0.8` 是过度触发” | **改措辞** | '
      '同一反向通量波形被切成 28 段，但 bit1 实际翻转仍为 14 次、码串严格模 8；'
      '正确说法是**事件分割假象** |')
    A('| 5 | “表达参数 ×0.5–2 具有鲁棒性” | **改措辞** | '
      '改为“单参数扰动 ×0.5–2 的已测点均通过”；无联合扰动与多初态扫描，'
      '不能称为鲁棒区间 |')
    A('| 6 | 表格中的“认证区间 / 几乎无约束 / 宽” | **改口径** | '
      '一律改为“已测离散点结果”，并分三层（`counting` / `event` / `full`）报告 |')
    A('')

    (HERE / '结论汇总.md').write_text('\n'.join(L), encoding='utf-8')
    print('wrote', HERE / '结论汇总.md')

    # ------------------------------------------------------------- figure
    fig, axes = plt.subplots(1, 2, figsize=(16.5, 6.6), dpi=150,
                             gridspec_kw=dict(width_ratios=[1.05, 1.35]))
    ax = axes[0]
    M = np.zeros((len(gates), len(folds)))
    for i, g in enumerate(gates):
        for j, f in enumerate(folds):
            r = cell(g, f)
            if r is None:
                M[i, j] = 2
            else:
                cnt, ev, full = tiers(r)
                M[i, j] = 2 if full else (1 if cnt else 0)
    cmap = matplotlib.colors.ListedColormap(['#c62828', '#ef9a00', '#2e7d32'])
    ax.imshow(M, cmap=cmap, vmin=0, vmax=2, aspect='auto')
    for i, g in enumerate(gates):
        for j, f in enumerate(folds):
            r = cell(g, f)
            if r is None:
                txt, col = 'PASS\n(frozen)', 'white'
            else:
                cnt, ev, full = tiers(r)
                if full:
                    txt, col = 'PASS', 'white'
                elif cnt:
                    txt, col = 'counts ok\nevent fails', 'black'
                else:
                    txt, col = 'FAIL', 'white'
            ax.text(j, i, txt, ha='center', va='center', fontsize=12,
                    color=col, fontweight='bold')
    ax.set_xticks(range(len(folds)))
    ax.set_xticklabels([f'x{f:g}' for f in folds], fontsize=14)
    ax.set_yticks(range(len(gates)))
    ax.set_yticklabels([f'n = {g:g}' for g in gates], fontsize=14)
    ax.set_xlabel('receiver_uM_per_au  (fold change)', fontsize=15)
    ax.set_ylabel('n_A1_gate', fontsize=15)
    ax.set_title('Paired grid (tested levels only): the three rows\n'
                 'are identical, so sharpening does not move the bracket',
                 fontsize=15)

    ax = axes[1]
    labels, c_lo, c_hi, f_lo, f_hi = [], [], [], [], []
    for name, pts in groups.items():
        pts = sorted(pts, key=lambda r: r['value'])
        fz = s1['frozen_config'].get(pts[0]['target']) or ZENG_FROZEN.get(name)
        if not fz:
            continue
        fold = [r['value'] / fz for r in pts]
        okc = [x for x, r in zip(fold, pts) if tiers(r)[0]]
        okf = [x for x, r in zip(fold, pts) if tiers(r)[2]]
        if not okc:
            continue
        labels.append(name)
        c_lo.append(min(okc)); c_hi.append(max(okc))
        f_lo.append(min(okf) if okf else None)
        f_hi.append(max(okf) if okf else None)
    order = np.argsort([-h for h in c_hi])
    labels = [labels[i] for i in order]
    c_lo = [c_lo[i] for i in order]; c_hi = [c_hi[i] for i in order]
    f_lo = [f_lo[i] for i in order]; f_hi = [f_hi[i] for i in order]
    y = np.arange(len(labels))
    ax.barh(y, [c_hi[i] - c_lo[i] for i in range(len(y))], left=c_lo,
            height=.34, color='#2e7d32', alpha=.85,
            label='counting arm passed at these levels')
    idx = [i for i in range(len(y)) if f_lo[i] is not None]
    ax.barh([y[i] for i in idx], [f_hi[i] - f_lo[i] for i in idx],
            left=[f_lo[i] for i in idx], height=.14, color='#1f4e79',
            alpha=.95, label='full criterion passed at these levels')
    ax.axvline(1.0, color='black', lw=2.2, ls='--')
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=13)
    ax.set_xscale('log')
    ax.set_xlabel('fold change from the frozen working point', fontsize=15)
    ax.set_title('Span of the TESTED levels that passed\n'
                 '(discrete sampling; not a continuous tolerance interval)',
                 fontsize=15)
    ax.tick_params(axis='x', labelsize=13)
    ax.grid(axis='x', alpha=.25)
    ax.invert_yaxis()
    ax.legend(fontsize=12, loc='lower right')
    fig.tight_layout()
    fig.savefig(HERE / 'fig_summary.png', dpi=150)
    fig.savefig(HERE / 'fig_summary.svg')
    plt.close(fig)
    print('wrote', HERE / 'fig_summary.png')


if __name__ == '__main__':
    main()
