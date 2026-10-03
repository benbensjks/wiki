"""Hand-drawing drafts for the hybrid counter: bit1 (middle) and carry1 -> bit2.

Two figures
-----------
fig_bit1_internal    the 6-state ZMH middle block: the carry-in gate
                     (A0 / F0 / x / I1) and the switch pb1 with its two output
                     chains (T1, RDF1) and the ALGEBRAIC complex (I1*RDF1)^2.
fig_carry1_wired     HBY carry1 (A1 / F1, two 3-state expression chains) wired
                     from bit1's pb1 into the g1 gate, plus the clock tap that
                     reads bit0's mature integrase, driving bit2.

Every connector corresponds to one term in the RHS of
`跨模型三级连接验证/hby_zmh_hby/model_hzh.py`:
  middle_rhs  :65-77   (A0, F0, pb1, I1, T1, RDF1)   -- uses donor self.z
  rhs         :87-90   (wiring: y[10]=b0_S -> middle; y[13]=pb1 and
                        y[2]=b0_I -> tail)
and of `跨模型三级连接验证/hybrid_model.py`:
  HbyReceiver.rhs     :200-218  (carry1 = A1/F1 chains + bit2)
  HbyReceiver.signals :189-198  (g1 gate and the clock)

TWO PARAMETER TABLES ARE IN PLAY and must not be mixed up:
  * A0/F0/pb1/I1/T1/RDF1 use the DONOR's own `p` (self.z, from
    `wiki任务/完整二级级联.py`)
  * A1/F1/bit2 use the ZENG table (self.tail.p, from
    `final_reconstruction/model.py`)
The symbols k_fwd / k_rev / K_D_comp / K_inh exist in BOTH with different
values (donor 0.118 / 0.08 / 3.2 / 0.1026 vs ZENG 7.0 / 5.0 / 1.2 / 0.1).

This is a WIRING DRAFT, not a data figure: no simulation output is read.

Glyph policy: only characters present in Microsoft YaHei (same as the bit0
drafts).  The bare "left tack" U+22A3 is NOT in YaHei and is written as
"抑" (or "⊣" replaced by the word) throughout.

Usage
-----
    & 'D:\\aconade\\python.exe' -B .\\make_bit1_carry1_sketches.py
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

HERE = Path(__file__).resolve().parent
BIT0 = HERE.parent / 'figures_bit0'
sys.path.insert(0, str(BIT0))

# reuse the verified Draft / geometry-check machinery from the bit0 drafts
# rather than forking a second copy that could drift out of sync.
import make_bit0_sketches as B
from tempo_style import (BLUE_GRAY, DEEP_BLUE, INK, LIGHT_BLUE, LIGHT_GRAY,
                         LIGHT_ORANGE, LIGHT_TEAL, ORANGE, RISK_ORANGE, TEAL,
                         WHITE, new_canvas)

matplotlib.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'DejaVu Sans']
matplotlib.rcParams['axes.unicode_minus'] = False
matplotlib.rcParams['svg.fonttype'] = 'none'

OUT = HERE / 'out'

Draft = B.Draft
_arc_points = B._arc_points


def save_pair(fig, stem):
    """Write both PNG and SVG; the SVG keeps editable text (svg.fonttype=none)."""
    path = OUT / f'{stem}.png'
    fig.savefig(path, dpi=220, bbox_inches='tight',
                facecolor=fig.get_facecolor())
    fig.savefig(OUT / f'{stem}.svg', bbox_inches='tight',
                facecolor=fig.get_facecolor())
    plt.close(fig)
    return path


def land_on(d, name, k):
    """Exact landing point + local tangent at sample `k` of connector `name`.

    A repression bar must sit ON a connector when it blocks a flux rather than
    a box; placing it by eye leaves it floating (that was the r9 defect in the
    bit0 drafts).  The tangent is used to lay the bar ACROSS the target.
    """
    for nm, p0, p1, rad in d.lines:
        if nm == name:
            pts = _arc_points(p0, p1, rad)
            k = max(2, min(k, len(pts) - 3))
            q0, q1 = pts[k - 2], pts[k + 2]
            tx, ty = q1[0] - q0[0], q1[1] - q0[1]
            nn = math.hypot(tx, ty) or 1.0
            return pts[k], (tx / nn, ty / nn)
    raise KeyError(name)


# --------------------------------------------------------------------------
# extra check: crossing connectors
# --------------------------------------------------------------------------
def _seg_cross(a, b, c, dd):
    """True if segment ab properly crosses segment cdd."""
    def cr(o, p, q):
        return ((p[0] - o[0]) * (q[1] - o[1]) -
                (p[1] - o[1]) * (q[0] - o[0]))
    d1, d2 = cr(c, dd, a), cr(c, dd, b)
    d3, d4 = cr(a, b, c), cr(a, b, dd)
    return ((d1 > 0) != (d2 > 0)) and ((d3 > 0) != (d4 > 0))


def check_crossings(d, tol=0.09):
    """Report connector pairs that cross each other.

    A deliberate T-junction -- one connector ENDS on another, as the repressor
    taps do -- is not a crossing, so a pair whose endpoint lies within `tol` of
    the other polyline is skipped.  This check does not exist for the bit0
    drafts: crossings are invisible to the overlap and dangling checks.
    """
    paths = {nm: _arc_points(p0, p1, rad) for nm, p0, p1, rad in d.lines}
    problems = []
    names = list(paths)
    for i in range(len(names)):
        for k in range(i + 1, len(names)):
            na, nb = names[i], names[k]
            pa, pb = paths[na], paths[nb]
            join = False
            for e in (pa[0], pa[-1]):
                if min(math.hypot(e[0] - q[0], e[1] - q[1]) for q in pb) < tol:
                    join = True
            for e in (pb[0], pb[-1]):
                if min(math.hypot(e[0] - q[0], e[1] - q[1]) for q in pa) < tol:
                    join = True
            if join:
                continue
            hit = False
            for s in range(len(pa) - 1):
                for t in range(len(pb) - 1):
                    if _seg_cross(pa[s], pa[s + 1], pb[t], pb[t + 1]):
                        hit = True
                        break
                if hit:
                    break
            if hit:
                problems.append(f'CONNECTORS CROSS: {na} x {nb}')
    print(f'--- crossings: {len(problems)} problem(s)')
    for p in problems:
        print('    ', p)
    return problems


def check_all(d, title):
    problems = B.check(d, title) + check_crossings(d)
    print(f'=== {title}: TOTAL {len(problems)} problem(s)')
    return problems


# --------------------------------------------------------------------------
# Figure 1: bit1 (ZMH middle, 6 states) internal flow
# --------------------------------------------------------------------------
def fig_bit1_internal():
    """Two panels.

    A single panel is over-constrained: RDF1 must both TAP the forward flux
    (a backward edge) and FEED the algebraic complex, and every placement that
    satisfies both forces a connector to cross another.  Splitting into
    "gate" and "switch + outputs" makes each subgraph planar, and the panel
    boundary is bridged with the net-label convention the bit0 drafts already
    use for long links (`I1 >` / `< I1`).
    """
    B.LAST_DRAFTS.clear()
    fig, ax = new_canvas(20.4, 13.6)
    d = Draft(ax, s=1.0)

    # ================= PANEL 1: the carry-in gate ==========================
    d.box('A0', 2.70, 11.70, 2.10, 0.74, 'A0\n进位激活臂', fc=LIGHT_TEAL,
          ec=TEAL, fs=9.4)
    d.box('F0', 2.70, 10.00, 2.10, 0.74, 'F0\n延迟抑制臂', fc=LIGHT_TEAL,
          ec=TEAL, fs=9.4)
    d.circle('g0', 5.30, 10.85, 0.40, '×', fc=WHITE, ec=TEAL, fs=14)
    d.box('I1', 7.90, 10.85, 2.10, 0.74, 'I1\n整合酶', fc=LIGHT_TEAL,
          ec=TEAL, fs=9.4)

    d.arrow('cin', (0.55, 11.70), (1.64, 11.70), color=DEEP_BLUE)
    d.label(0.45, 12.14, 'b0_S  (= PB0)', fs=9.4, color=DEEP_BLUE, ha='left')
    d.label(0.45, 11.26, 'A0 源 = 0.1413·(1 − b0_S)', fs=8.6,
            color=BLUE_GRAY, ha='left')
    d.arrow('a0g', (3.75, 11.70), (5.20, 11.14), color=TEAL)
    d.arrow('a0f', (2.70, 11.33), (2.70, 10.42), color=TEAL)
    d.repress('f0g', (3.75, 10.00), (5.22, 10.56), color=RISK_ORANGE)
    d.arrow('gI', (5.70, 10.85), (6.81, 10.85), color=TEAL)
    d.label(5.30, 9.34,
            'I1 源 = 0.3788·H(A0;1.0,2)·(1 − H(F0;2.6454,5.4332))\n'
            'F0 源 = 0.0979·H(A0;1.0,2)',
            fs=8.6, color=BLUE_GRAY, ha='center')
    d.label(9.55, 11.45, 'I1 >', fs=10.0, color=DEEP_BLUE, weight='bold',
            ha='center')
    d.label(9.55, 10.25, 'I1 >', fs=10.0, color=DEEP_BLUE, weight='bold',
            ha='center')

    # ================= PANEL 2: the switch and its outputs =================
    d.circle('pb1', 14.20, 9.60, 0.98, 'pb1', fc=WHITE, ec=DEEP_BLUE, fs=13)
    d.box('RDF1', 11.60, 7.60, 2.10, 0.74, 'RDF1', fc=LIGHT_ORANGE,
          ec=ORANGE, fs=9.4)
    d.box('T1', 11.60, 5.20, 2.10, 0.74, 'T1\nRep', fc=LIGHT_ORANGE,
          ec=ORANGE, fs=9.4)
    d.circle('cx', 14.20, 5.60, 0.46, '×', fc=WHITE, ec=BLUE_GRAY, fs=13)

    # the two inputs that cross the panel boundary, by net label
    d.arrow('vf', (11.00, 10.20), (13.26, 9.88), color=DEEP_BLUE)
    d.label(10.90, 10.54, '< I1', fs=10.0, color=DEEP_BLUE, weight='bold',
            ha='right')
    d.label(11.30, 10.78, 'vf（正向）', fs=8.6, color=DEEP_BLUE, ha='left')
    d.arrow('i1cx', (11.00, 5.60), (13.72, 5.60), color=DEEP_BLUE)
    d.label(10.90, 5.94, '< I1', fs=10.0, color=DEEP_BLUE, weight='bold',
            ha='right')
    d.label(15.30, 5.60, '代数复体 (I1·RDF1)²（无 C 状态）', fs=8.6,
            color=BLUE_GRAY, ha='left')
    # pb1 fans out to the two output species by net label
    d.arrow('pb1r', (9.90, 7.60), (10.53, 7.60), color=ORANGE)
    d.label(9.80, 7.94, '< PB1', fs=10.0, color=ORANGE, weight='bold',
            ha='right')
    d.arrow('pb1t', (9.90, 5.20), (10.53, 5.20), color=ORANGE)
    d.label(9.80, 4.95, '< PB1', fs=10.0, color=ORANGE, weight='bold',
            ha='right')
    d.label(15.60, 9.60, 'PB1 >', fs=10.0, color=ORANGE, weight='bold',
            ha='left')
    # RDF1 blocks the forward flux: the bar lands ON vf, laid across it
    _p, _tan = land_on(d, 'vf', 46)
    d.repress('rdfv', (10.90, 7.97), _p, color=RISK_ORANGE, rad=-0.04,
              bar_dir=(-_tan[1], _tan[0]), bar_at_tip=True)
    d.label(11.60, 11.28, 'RDF1 抑制 vf（同 bit0 的 ⑨）', fs=8.6,
            color=RISK_ORANGE, ha='center')
    # RDF1 feeds the complex, the complex drives the reverse flux
    d.arrow('Rcx', (12.20, 7.23), (13.76, 5.78), color=ORANGE)
    d.arrow('vr', (14.20, 6.06), (14.20, 8.62), color=ORANGE)
    d.label(14.72, 7.34, 'vr（反向）', fs=8.6, color=ORANGE, ha='left')
    # T1 represses RDF1
    d.repress('Trdf', (11.60, 5.57), (11.60, 7.23), color=RISK_ORANGE)
    d.label(8.60, 6.40, 'T1 抑制 RDF1', fs=8.6, color=RISK_ORANGE, ha='right')

    ax.text(10.2, 13.28,
            'bit1 —— ZMH 约化位（6 状态）的内部流程机制',
            ha='center', va='center', fontsize=16.5, color=DEEP_BLUE,
            fontweight='bold')
    ax.text(10.2, 12.88,
            '每个框 = model_hzh.py middle_rhs 的一个状态（:73-77）'
            '　·　参数取自 donor 自己的 p（完整二级级联.py）',
            ha='center', va='center', fontsize=10.0, color=BLUE_GRAY)
    ax.text(5.30, 12.40, '① 进位门：A0 激活、F0 延迟抑制', ha='center',
            va='center', fontsize=11.5, color=TEAL, fontweight='bold')
    ax.text(13.20, 12.40, '② 开关 pb1 与两条输出链', ha='center',
            va='center', fontsize=11.5, color=DEEP_BLUE, fontweight='bold')
    ax.plot([10.30, 10.30], [4.80, 11.90], color=BLUE_GRAY, lw=1.1,
            ls=(0, (5, 4)))

    y0 = 4.10
    ax.plot([0.6, 19.4], [4.60, 4.60], color=BLUE_GRAY, lw=1.2,
            ls=(0, (5, 4)))
    ax.text(0.85, 4.48, '图例：本图画的连接', ha='left', va='top',
            fontsize=12.5, color=DEEP_BLUE, fontweight='bold')
    col1 = [
        '【进位输入（来自 bit0）】',
        'b0_S → A0：源项 0.1413·(1 − b0_S)，A0 唯一天然驱动',
        '   · 与 carry0 的 E 是同一条 PB0 线（PB0 = 1 − b0_S）',
        '',
        '【非相干前馈环 I1-FFL（面板 ①）】',
        'A0 同时激活 I1 与 F0；F0 延迟地抑制 I1',
        '   · I1 源 = 0.3788·H(A0;1.0,2)·(1 − H(F0;2.6454,5.4332))',
        '   · "延迟臂"：A0 先开出 I1，F0 随后把 I1 关掉',
    ]
    col2 = [
        '【开关 pb1（= PB1，面板 ②）】',
        'dpb1/dt = −vf + vr（正向减 pb1、反向加 pb1）',
        'vf = 0.118·pb1·H(I1;1.8,4)·0.1026/(0.1026+RDF1)',
        'vr = 0.08·(1−pb1)·(I1·RDF1)²/(3.2² + (I1·RDF1)²)',
        '',
        '【为什么分成两个面板】',
        'RDF1 既要抑制 vf（反向边）、又要喂复体，单面板必然交叉；',
        '拆开后各自可平面化，跨面板的线用 net label 衔接。',
        '【net label  > / < = 同名标签表示同一条连线】',
    ]
    col3 = [
        '【两条输出链】',
        'PB1 驱动 T1：源项 0.2·pb1，衰减 0.005·T1',
        'LR1 = 1 − pb1 驱动 RDF1，且 T1 抑制 RDF1',
        'RDF1 反过来抑制 vf（与 bit0 的 ⑨ 同义）',
        '',
        '【与 bit0 的结构差别】',
        'bit1 没有 C 状态：复体是代数量 (I1·RDF1)²',
        'bit0/bit2 有显式 C 状态 + k_complex = q·K_D_comp',
    ]
    for col, x in ((col1, 0.95), (col2, 7.10), (col3, 13.30)):
        for i, line in enumerate(col):
            ax.text(x, y0 - 0.05 - i * 0.315, line, ha='left', va='center',
                    fontsize=9.3, color=INK)

    check_all(d, 'fig_bit1_internal')
    return save_pair(fig, 'fig_bit1_internal')


# --------------------------------------------------------------------------
# Figure 2: carry1 (HBY) -> bit2
# --------------------------------------------------------------------------
def fig_carry1_wired():
    B.LAST_DRAFTS.clear()
    fig, ax = new_canvas(26.0, 15.0)
    d = Draft(ax, s=1.0)

    def container(x, y, w, h, ec, title, tx=None, ty=None, fs=10.2):
        ax.add_patch(FancyBboxPatch(
            (x, y), w, h, boxstyle='round,pad=0,rounding_size=0.25',
            fc='none', ec=ec, lw=1.6, zorder=1))
        ax.text(tx if tx is not None else x + 0.28,
                ty if ty is not None else y + h - 0.30, title, fontsize=fs,
                color=ec, ha='left', va='center', zorder=2)

    # ---- blocks -----------------------------------------------------------
    d.box('b1.blk', 3.40, 10.20, 3.60, 3.60,
          'bit1（ZMH 约化位）\n\npb1 / I1\nT1 / RDF1\n\n'
          '输出\nPB1 = pb1', fc=LIGHT_TEAL, ec=TEAL, fs=9.6)
    d.box('b2.blk', 23.00, 10.20, 3.80, 3.60,
          'bit2（HBY 11 状态）\n\nM_I / I_u / I\nM_T / T_u / T\n'
          'M_R / R_u / R / C / S\n\n输入\nα_Int[1]·g1', fc=LIGHT_BLUE,
          ec=DEEP_BLUE, fs=9.6)

    container(7.00, 5.60, 13.60, 8.20, TEAL, 'carry1（HBY · 6 状态）')

    # A1 chain (vertical, top = mRNA)
    d.box('M_A1', 10.00, 12.60, 2.30, 0.72, 'M_A1\nA1 的 mRNA', fc=LIGHT_GRAY,
          ec=TEAL, fs=9.2)
    d.box('A1_u', 10.00, 11.40, 2.30, 0.72, 'A1_u\n未成熟 A1', fc=LIGHT_GRAY,
          ec=TEAL, fs=9.2)
    d.box('A1', 10.00, 10.20, 2.30, 0.72, 'A1\n（激活臂）', fc=LIGHT_TEAL,
          ec=TEAL, fs=9.4)
    # F1 chain, placed ENTIRELY BELOW the A1 row so the A1 -> gate line runs
    # clear above it (a chain in the same row would be crossed by that line).
    d.box('M_F1', 13.30, 9.20, 2.30, 0.72, 'M_F1\nF1 的 mRNA', fc=LIGHT_GRAY,
          ec=TEAL, fs=9.2)
    d.box('F1_u', 13.30, 8.00, 2.30, 0.72, 'F1_u\n未成熟 F1', fc=LIGHT_GRAY,
          ec=TEAL, fs=9.2)
    d.box('F1', 13.30, 6.80, 2.30, 0.72, 'F1\n（抑制臂）', fc=LIGHT_TEAL,
          ec=TEAL, fs=9.4)
    d.circle('g1', 17.10, 10.20, 0.42, '×', fc=WHITE, ec=TEAL, fs=14)
    d.box('ck', 17.80, 7.00, 4.60, 1.00,
          'clock 门  H(b0_I ; 0.3, 2.0)', fc=LIGHT_TEAL, ec=TEAL, fs=9.6)

    # ---- connectors -------------------------------------------------------
    # pb1 -> A1 chain source, with A1's negative autoregulation as a T-bar
    d.arrow('pb_A1', (5.20, 11.30), (8.82, 12.60), color=DEEP_BLUE, rad=-0.10)
    d.label(5.60, 11.02, 'PB1 = pb1', fs=9.0, color=DEEP_BLUE, ha='left')
    _p, _tan = land_on(d, 'pb_A1', 30)
    # exits A1's LEFT edge so the tap stays clear of the A1_u box above
    d.repress('auto', (8.87, 10.34), _p, color=RISK_ORANGE, rad=0.05,
              bar_dir=(-_tan[1], _tan[0]), bar_at_tip=True)
    d.label(5.70, 10.32, 'A1 负自调  1 − H(A1;0.6,2)', fs=8.6,
            color=RISK_ORANGE, ha='left')
    # the two chains
    d.arrow('A1_u_A1', (10.00, 12.24), (10.00, 11.76), color=TEAL)
    d.arrow('M_A1_A1_u', (10.00, 11.04), (10.00, 10.56), color=TEAL)
    d.arrow('M_F1_F1_u', (13.30, 8.84), (13.30, 8.36), color=TEAL)
    d.arrow('F1_u_F1', (13.30, 7.64), (13.30, 7.16), color=TEAL)
    d.label(10.00, 13.32,
            'A1 源 = 16.0·pb1·(1 − H(A1;0.6,2))', fs=8.8, color=BLUE_GRAY,
            ha='center')
    d.label(13.30, 12.60,
            'F1 源 = 5.0·H(A1;1.2,4)', fs=8.8, color=BLUE_GRAY, ha='center')
    # A1 -> F1 chain
    d.arrow('A1_MF1', (10.90, 9.84), (12.13, 9.20), color=TEAL, rad=-0.10)
    # A1 -> gate, along the clear lane between the A1 row and the F1 chain
    d.arrow('A1_g1', (11.15, 10.20), (16.66, 10.20), color=TEAL)
    d.label(9.60, 8.55,
            'A1 进门用 n = 6.0；驱动 F1 时用 n = 4.0\n'
            '同一个 K_A[1] = 1.2，两个不同的 Hill 系数',
            fs=8.8, color=BLUE_GRAY, ha='center')
    d.arrow('F1_g1', (14.45, 6.80), (16.80, 9.90), color=TEAL, rad=-0.10)
    d.arrow('g1_b2', (17.52, 10.20), (21.08, 10.20), color=DEEP_BLUE)
    # clock
    d.arrow('ck_g1', (17.40, 7.50), (17.12, 9.80), color=TEAL)
    d.arrow('bi_ck', (14.80, 7.00), (15.46, 7.00), color=DEEP_BLUE,
            ls=(0, (5, 4)), lw=1.6, ms=12)
    d.label(14.70, 7.36, '< b0_I', fs=9.4, color=DEEP_BLUE, weight='bold',
            ha='right')

    # ---- formula call-outs ------------------------------------------------
    d.label(17.30, 5.95,
            'g1 = H(A1 ; 1.2, 6) · (1 − H(F1 ; 0.4, 4)) · clock\n'
            'Int2 源 = α_Int[1]·g1 = 38.0·g1',
            fs=9.0, color=TEAL, ha='center')

    ax.text(13.0, 14.50,
            'bit1 接上 carry1 与 bit2 之后的完整流程机制',
            ha='center', va='center', fontsize=17, color=DEEP_BLUE,
            fontweight='bold')
    ax.text(13.0, 14.02,
            'carry1 = HBY 的 6 状态块（hybrid_model.py:200-218）'
            '　·　门 g1 与时钟见 :189-198',
            ha='center', va='center', fontsize=10.2, color=BLUE_GRAY)

    ax.plot([0.6, 25.4], [5.00, 5.00], color=BLUE_GRAY, lw=1.2,
            ls=(0, (5, 4)))
    ax.text(0.85, 4.88, '图例：本图画的连接与参数出处', ha='left', va='top',
            fontsize=12.5, color=DEEP_BLUE, fontweight='bold')
    col1 = [
        '【本图画的连接】',
        'PB1 = pb1 → A1 链的源（A1 负自调压住这条源）',
        'A1 → F1 链的源：F1 由 A1 激活（H(A1;1.2,4)）',
        'A1 → g1、F1 抑制 g1、clock → g1、g1 → bit2',
        '',
        '【两条链的三级】',
        'M → 未成熟 → 成熟，与 bit0 的 a→b→c 同构',
        'carry 链：mRNA 半衰期 2.0 min、成熟 32.5 min',
        'bit2 链：mRNA 半衰期 2.0 min、成熟 20.0 min',
    ]
    col2 = [
        '【三因子 AND 门 g1】',
        'g1 = H(A1;1.2,6) · (1 − H(F1;0.4,4)) · clock',
        '· A1 进门用 n = 6.0（n_A1_gate），驱动 F1 时用 n = 4.0',
        '· K_A[1] = 1.2 两处相同；K_F[1] = 0.4、n_F[1] = 4.0',
        '· F1 只抑制门，不抑制 A1',
        '· 无时钟时 g1 = 0 —— 这就是第三级的时序闸门',
        '',
        '【Int2 去 bit2】',
        'α_Int[1]·g1 = 38.0·g1 → bit2 的 I 链（K_D_int[2] = 0.6）',
    ]
    col3 = [
        '【clock 回绕（本图只画接口）】',
        'clock = H(b0_I ; 0.3, 2.0)',
        '   · 读的是 bit0 的成熟游离整合酶 b0_I',
        '   · 跨过 carry0 直接门控 carry1，是本模型的关键设计',
        '   · 命名沿用 HBY：K = 0.3, n = 2.0（frozen extension）',
        '',
        '【本图未画】',
        'bit0、carry0、bit2 的内部流程；各位的清除/衰减项',
    ]
    for col, x in ((col1, 0.95), (col2, 9.10), (col3, 18.10)):
        for i, line in enumerate(col):
            ax.text(x, 4.50 - i * 0.315, line, ha='left', va='center',
                    fontsize=9.3, color=INK)

    check_all(d, 'fig_carry1_wired')
    return save_pair(fig, 'fig_carry1_wired')


def main():
    for o in (fig_bit1_internal(), fig_carry1_wired()):
        print('wrote', o)


if __name__ == '__main__':
    main()
