"""Hand-drawing drafts for the hybrid three-level counter: bit0.

Two figures
-----------
fig_bit0_internal          bit0 only: the 11 species and every reaction, with a
                           legend mapping each abstract label to its real
                           biochemical component.
fig_bit0_wired             the same bit0 block wired to the upstream oscillator
                           and to the downstream carry0 -> bit1, plus the clock
                           tap that later feeds carry1.

Every connector corresponds to one term in the RHS of
`跨模型三级连接验证/hby_zmh_hby/model_hzh.py` (bit0) or
`跨模型三级连接验证/hybrid_model.py` (carry0 / bit1); the source line ranges are
named in the legend so the sketch can be audited.

This is a WIRING DRAFT, not a data figure: no simulation output is read.

Glyph policy: only characters present in Microsoft YaHei are used.  In
particular `⇄ ⊣ ⊗ ▸ ◂ ⑪ ⑫` are NOT available and have been replaced by
`/`, `抑制`, `×`, `>`/`<` and `(11)`/`(12)`.

Usage
-----
    & 'D:\\aconade\\python.exe' -B .\\make_bit0_sketches.py
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from tempo_style import (BLUE_GRAY, DEEP_BLUE, INK, LIGHT_BLUE, LIGHT_GRAY,
                         LIGHT_ORANGE, LIGHT_TEAL, ORANGE, RISK_ORANGE, TEAL,
                         WHITE, new_canvas, save)

matplotlib.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'DejaVu Sans']
matplotlib.rcParams['axes.unicode_minus'] = False
matplotlib.rcParams['svg.fonttype'] = 'none'   # keep real <text> in the SVG

OUT = HERE / 'out'


def save_pair(fig, stem):
    """Write both PNG and SVG; the SVG keeps editable text (svg.fonttype=none)."""
    path = OUT / f'{stem}.png'
    fig.savefig(path, dpi=220, bbox_inches='tight',
                facecolor=fig.get_facecolor())
    fig.savefig(OUT / f'{stem}.svg', bbox_inches='tight',
                facecolor=fig.get_facecolor())
    plt.close(fig)
    return path

# --------------------------------------------------------------------------
# bit0 geometry in module-local units (1 unit = 1 inch at scale 1)
# --------------------------------------------------------------------------
SPECIES = {
    'a': (2.60, 10.90, 2.00, 0.70, 'a', 'M_I'),
    'b': (2.60,  9.80, 2.00, 0.70, 'b', 'I_u'),
    'c': (2.60,  8.70, 2.00, 0.70, 'c', 'I'),
    'd': (11.00, 10.70, 2.00, 0.70, 'd', 'M_T'),
    'e': (11.00,  9.75, 2.00, 0.70, 'e', 'T_u'),
    'f': (11.00,  8.80, 2.00, 0.70, 'f', 'T'),
    'g': (11.00,  7.15, 2.00, 0.70, 'g', 'M_R'),
    'h': (11.00,  6.20, 2.00, 0.70, 'h', 'R_u'),
    'i': (11.00,  5.25, 2.00, 0.70, 'i', 'R'),
    'j': ( 5.40,  5.40, 2.40, 0.90, 'j', 'C'),
}
K_CENTRE = (7.20, 7.90)
K_RADIUS = 0.95
R8_RAD = -0.18          # curvature of the forward-flux arrow (shared with r9)

FC = {'a': LIGHT_GRAY, 'b': LIGHT_GRAY, 'c': LIGHT_TEAL, 'd': LIGHT_GRAY,
      'e': LIGHT_GRAY, 'f': LIGHT_ORANGE, 'g': LIGHT_GRAY, 'h': LIGHT_GRAY,
      'i': LIGHT_ORANGE, 'j': LIGHT_BLUE}


LAST_DRAFTS: list = []   # geometry of the most recently drawn figure(s)


class Draft:
    """Draws and simultaneously records geometry for the self-check."""

    def __init__(self, ax, s=1.0):
        self.ax, self.s = ax, s
        self.boxes = []          # (name, cx, cy, w, h) in DATA coordinates
        self.lines = []          # (name, p0, p1, rad)
        self.fs = max(6.0, 9.4 * s)
        LAST_DRAFTS.append(self)

    # ---- primitives ------------------------------------------------------
    def box(self, name, cx, cy, w, h, text, fc=LIGHT_GRAY, ec=DEEP_BLUE,
            fs=None, ls='-', zorder=3):
        self.ax.add_patch(FancyBboxPatch(
            (cx - w / 2, cy - h / 2), w, h,
            boxstyle='round,pad=0,rounding_size=0.12', fc=fc, ec=ec,
            lw=1.6, ls=ls, zorder=zorder))
        self.ax.text(cx, cy, text, ha='center', va='center',
                     fontsize=fs or self.fs, color=INK, zorder=zorder + 1,
                     linespacing=1.3)
        self.boxes.append((name, cx, cy, w, h))

    def circle(self, name, cx, cy, r, text, fc=LIGHT_TEAL, ec=DEEP_BLUE,
               fs=None, zorder=5):
        self.ax.add_patch(Circle((cx, cy), r, fc=fc, ec=ec, lw=1.8,
                                 zorder=zorder))
        self.ax.text(cx, cy, text, ha='center', va='center',
                     fontsize=fs or self.fs, color=ec, zorder=zorder + 1)
        self.boxes.append((name, cx, cy, 2 * r, 2 * r))

    def arrow(self, name, p0, p1, color=DEEP_BLUE, rad=0.0, lw=1.8, ms=13,
              ls='-', zorder=4):
        self.ax.add_patch(FancyArrowPatch(
            p0, p1, arrowstyle='-|>', mutation_scale=ms, lw=lw, color=color,
            linestyle=ls, zorder=zorder, shrinkA=1.5, shrinkB=1.5,
            connectionstyle=f'arc3,rad={rad}'))
        self.lines.append((name, p0, p1, rad))

    def repress(self, name, p0, p1, color=DEEP_BLUE, rad=0.0, lw=1.8,
                zorder=4, bar=0.14, bar_dir=None, bar_at_tip=False):
        dx, dy = p1[0] - p0[0], p1[1] - p0[1]
        n = math.hypot(dx, dy)
        ux, uy = dx / n, dy / n
        if bar_at_tip:
            # shaft end == bar centre == p1.  Needed when the tip lands exactly
            # ON another connector: backing the shaft off by `bar` would push
            # the bar off that connector and make it read as a parallel line.
            ex, ey = p1[0], p1[1]
        else:
            ex, ey = p1[0] - ux * bar, p1[1] - uy * bar
        self.ax.add_patch(FancyArrowPatch(
            p0, (ex, ey), arrowstyle='-', lw=lw, color=color, zorder=zorder,
            shrinkA=1.5, shrinkB=0.0, connectionstyle=f'arc3,rad={rad}'))
        if bar_dir is not None:
            # explicit bar direction: used when the tip lands ON another
            # connector, so the bar sits squarely across the flux it blocks
            # whatever angle this line happens to arrive from.  Supersedes the
            # old bar_tangent flag (removed): this covers that case and every
            # other one, and unlike a flag it is actually exercised by a figure.
            bn = math.hypot(*bar_dir) or 1.0
            px, py = bar_dir[0] / bn, bar_dir[1] / bn
        else:
            px, py = -uy, ux
        self.ax.plot([ex - px * bar, ex + px * bar],
                     [ey - py * bar, ey + py * bar], color=color, lw=lw + 0.6,
                     solid_capstyle='round', zorder=zorder + 1)
        self.lines.append((name, p0, p1, rad))

    def bind(self, name, p0, p1, color=INK, rad=0.0, lw=1.5, zorder=4):
        self.ax.add_patch(FancyArrowPatch(
            p0, p1, arrowstyle='<|-|>', mutation_scale=11, lw=lw, color=color,
            zorder=zorder, shrinkA=1.5, shrinkB=1.5,
            connectionstyle=f'arc3,rad={rad}'))
        self.lines.append((name, p0, p1, rad))

    def label(self, x, y, text, fs=None, color=INK, ha='center', va='center',
              style='normal', weight='normal', zorder=8):
        self.ax.text(x, y, text, ha=ha, va=va, fontsize=fs or self.fs * 0.88,
                     color=color, style=style, fontweight=weight, zorder=zorder,
                     linespacing=1.3)


def _pt(key, s, ox, oy, dx=0.0, dy=0.0):
    cx, cy = SPECIES[key][0], SPECIES[key][1]
    return (ox + s * (cx + dx), oy + s * (cy + dy))


def _rect(key, s, ox, oy):
    cx, cy, w, h = SPECIES[key][:4]
    return (ox + s * cx, oy + s * cy, s * w, s * h)


# --------------------------------------------------------------------------
# the bit0 module (shared by both figures)
# --------------------------------------------------------------------------
def draw_bit0(d, ox=0.0, oy=0.0, compact=False):
    s, fs = d.s, d.fs
    for key, (cx, cy, w, h, letter, code) in SPECIES.items():
        X, Y, W, H = _rect(key, s, ox, oy)
        d.box(f'b0.{key}', X, Y, W, H, f'{letter}  {code}', fc=FC[key], fs=fs)

    kx, ky = ox + s * K_CENTRE[0], oy + s * K_CENTRE[1]
    r = K_RADIUS * s
    d.circle('b0.k', kx, ky, r, 'k', fc=LIGHT_TEAL)
    d.label(kx, ky - r - 0.30 * s, 'k :  S (LR)  /  PB = 1 − S',
            fs=fs * 0.92, color=DEEP_BLUE, weight='bold')

    # (1) external drive into M_I
    d.arrow('r1', _pt('a', s, ox, oy, 0, 0.92), _pt('a', s, ox, oy, 0, 0.40),
            color=ORANGE)
    d.label(*_pt('a', s, ox, oy, -1.60, 0.70), '① 转录\n(外源)', fs=fs * 0.82,
            color=ORANGE)
    # (2) (3) integrase expression chain
    d.arrow('r2', _pt('a', s, ox, oy, 0, -0.36), _pt('b', s, ox, oy, 0, 0.36),
            color=TEAL)
    d.label(*_pt('a', s, ox, oy, 1.32, -0.55), '② 翻译', fs=fs * 0.84,
            color=TEAL, ha='left')
    d.arrow('r3', _pt('b', s, ox, oy, 0, -0.36), _pt('c', s, ox, oy, 0, 0.36),
            color=TEAL)
    d.label(*_pt('b', s, ox, oy, 1.32, -0.55), '③ 成熟', fs=fs * 0.84,
            color=TEAL, ha='left')

    # (4) Rep chain, (5) RDF chain
    for a_, b_ in (('d', 'e'), ('e', 'f')):
        d.arrow(f'rep{a_}{b_}', _pt(a_, s, ox, oy, 0, -0.36),
                _pt(b_, s, ox, oy, 0, 0.36), color=TEAL)
    for a_, b_ in (('g', 'h'), ('h', 'i')):
        d.arrow(f'rdf{a_}{b_}', _pt(a_, s, ox, oy, 0, -0.36),
                _pt(b_, s, ox, oy, 0, 0.36), color=TEAL)
    d.label(*_pt('d', s, ox, oy, -1.68, -0.10), '④ Rep 链\n(PB 驱动)',
            fs=fs * 0.80, color=TEAL, ha='right')
    d.label(*_pt('g', s, ox, oy, -1.68, -0.10), '⑤ RDF 链\n(LR 驱动)',
            fs=fs * 0.80, color=TEAL, ha='right')

    # (6) T represses RDF expression
    d.repress('r6', _pt('f', s, ox, oy, 0.30, -0.40),
              _pt('g', s, ox, oy, 0.30, 0.40), color=RISK_ORANGE)
    if compact:
        d.label(*_pt('f', s, ox, oy, -1.10, -0.92), '⑥ T 抑制 RDF',
                fs=fs * 0.80, color=RISK_ORANGE, ha='right')
    else:
        d.label(*_pt('f', s, ox, oy, 1.42, -1.00),
                '⑥ T 抑制 RDF 表达\n   K_rep = 0.85, n_rep = 3.9',
                fs=fs * 0.80, color=RISK_ORANGE, ha='left')

    # (7) complex formation
    d.bind('r7a', _pt('c', s, ox, oy, 0, -0.38),
           _pt('j', s, ox, oy, -0.62, 0.48), color=INK, rad=-0.14)
    d.bind('r7b', _pt('i', s, ox, oy, -1.03, -0.04),
           _pt('j', s, ox, oy, 1.23, 0.04), color=INK, rad=0.12)
    if not compact:
        d.label(*_pt('j', s, ox, oy, 0.0, -0.80),
                '⑦ 可逆结合   k_on = 0.1, k_off = 1, δ_C = 1',
                fs=fs * 0.80, color=INK)

    # (8) forward flux into the switch, (9) RDF represses it.
    # ⑨ must LAND on the ⑧ flux: a floating stub next to box i reads as
    # "RDF represses something unnamed" and was invisible to the overlap
    # checks, so the T-bar is now placed on the r8 arc itself.
    _r8p0 = _pt('c', s, ox, oy, 1.05, 0.06)
    _r8p1 = (kx - r - 0.05 * s, ky + 0.24 * s)
    d.arrow('r8', _r8p0, _r8p1, color=DEEP_BLUE, rad=R8_RAD)
    d.label(*_pt('c', s, ox, oy, 2.45, 1.00), '⑧ 正向重组\n(I 激活)',
            fs=fs * 0.80, color=DEEP_BLUE, ha='left')
    _p8 = _arc_points(_r8p0, _r8p1, R8_RAD)
    _land8 = _p8[18]                                    # t = 0.30 along r8
    _tan8 = (_p8[20][0] - _p8[16][0], _p8[20][1] - _p8[16][1])
    # ⑨ arrives at 130.7 deg to r8's tangent, so a bar drawn perpendicular to
    # ⑨ itself crosses the flux at only 40.7 deg and reads as a skewed slash.
    # Laying the bar along r8's LOCAL NORMAL instead makes it straddle the flux
    # squarely at any arrival angle, and bar_at_tip keeps it centred on r8
    # rather than pushed off it by the usual shaft-end offset.
    d.repress('r9', _pt('i', s, ox, oy, 1.00, 0.06), _land8,
              color=RISK_ORANGE, rad=-0.30,
              bar_dir=(-_tan8[1], _tan8[0]), bar_at_tip=True)
    if compact:
        d.label(_land8[0], _land8[1] - 0.78, '⑨', fs=fs * 0.80,
                color=RISK_ORANGE)
    else:
        d.label(_land8[0], _land8[1] - 0.78, '⑨ RDF 抑制正向',
                fs=fs * 0.80, color=RISK_ORANGE)

    # (10) reverse flux from the complex
    d.arrow('r10', _pt('j', s, ox, oy, 0.62, 0.50),
            (kx - 0.34 * s, ky - r - 0.05 * s), color=ORANGE, rad=-0.26)
    d.label(*_pt('j', s, ox, oy, 2.55, 0.72), '⑩ 反向重组\n(C 驱动)',
            fs=fs * 0.80, color=ORANGE, ha='left')

    # (11) (12) the switch drives both expression chains
    d.arrow('r11', (kx + 0.55 * s, ky + r * 0.82),
            _pt('d', s, ox, oy, -1.12, 0.0), color=DEEP_BLUE, rad=-0.20)
    d.label(*_pt('d', s, ox, oy, 0.05, 1.10), '(11) PB', fs=fs * 0.82,
            color=DEEP_BLUE)
    d.arrow('r12', (kx + r * 0.94, ky - 0.18 * s),
            _pt('g', s, ox, oy, -1.12, 0.0), color=DEEP_BLUE, rad=0.18)
    d.label(*_pt('g', s, ox, oy, 0.05, 1.10), '(12) LR', fs=fs * 0.82,
            color=DEEP_BLUE)


# --------------------------------------------------------------------------
# self-check
# --------------------------------------------------------------------------
def _arc_points(p0, p1, rad, n=60):
    x1, y1 = p0
    x2, y2 = p1
    mx, my = (x1 + x2) / 2.0, (y1 + y2) / 2.0
    dx, dy = x2 - x1, y2 - y1
    cx, cy = mx + rad * dy, my - rad * dx
    pts = []
    for k in range(n + 1):
        t = k / n
        u = 1 - t
        pts.append((u * u * x1 + 2 * u * t * cx + t * t * x2,
                    u * u * y1 + 2 * u * t * cy + t * t * y2))
    return pts


def _inside(px, py, box, pad=0.0):
    _, bx, by, bw, bh = box
    return abs(px - bx) < bw / 2 + pad and abs(py - by) < bh / 2 + pad


def _rect_distance(px, py, box):
    """Distance from a point to a box rectangle (0.0 when inside)."""
    _, bx, by, bw, bh = box
    dx = max(abs(px - bx) - bw / 2, 0.0)
    dy = max(abs(py - by) - bh / 2, 0.0)
    return math.hypot(dx, dy)


# A connector tip within TIP_TOL_BOX of a box, or within TIP_TOL_LINE of
# another connector, counts as terminated.  A tip that satisfies neither is
# DANGLING: it points at empty space, which the overlap checks below cannot
# see.  Names in OPEN_ENDED_TIPS are deliberately open by design (a source
# stub whose tail leaves from nowhere) and must carry a label saying so.
TIP_TOL_BOX = 0.30
TIP_TOL_LINE = 0.12
OPEN_ENDED_TIPS: set = set()


def check(d, title):
    problems = []
    for i in range(len(d.boxes)):
        n1, x1, y1, w1, h1 = d.boxes[i]
        for k in range(i + 1, len(d.boxes)):
            n2, x2, y2, w2, h2 = d.boxes[k]
            if (abs(x1 - x2) < (w1 + w2) / 2 and
                    abs(y1 - y2) < (h1 + h2) / 2):
                problems.append(f'BOX OVERLAP: {n1} <-> {n2}')
    for name, p0, p1, rad in d.lines:
        # a box is the source/target if one of the endpoints sits on it
        attached = {b[0] for b in d.boxes
                    if _inside(p0[0], p0[1], b, 0.20)
                    or _inside(p1[0], p1[1], b, 0.20)}
        pts = _arc_points(p0, p1, rad)[4:-4]
        for b in d.boxes:
            if b[0] in attached:
                continue
            hit = sum(1 for (px, py) in pts
                      if _inside(px, py, b, -0.06))
            if hit:
                problems.append(f'ARROW THROUGH BOX: {name} -> {b[0]} '
                                f'({hit} samples)')
    # every tip must land on a box or on another connector
    paths = {nm: _arc_points(o0, o1, orad) for nm, o0, o1, orad in d.lines}
    for name, p0, p1, rad in d.lines:
        if name in OPEN_ENDED_TIPS:
            continue
        dbox = min((_rect_distance(p1[0], p1[1], b) for b in d.boxes),
                   default=float('inf'))
        dline = min((min(math.hypot(p1[0] - qx, p1[1] - qy)
                         for qx, qy in pts)
                     for nm, pts in paths.items() if nm != name),
                    default=float('inf'))
        if dbox > TIP_TOL_BOX and dline > TIP_TOL_LINE:
            problems.append(
                f'DANGLING TIP: {name} ends at ({p1[0]:.2f}, {p1[1]:.2f}); '
                f'nearest box {dbox:.2f}, nearest connector {dline:.2f}')
    print(f'--- {title}: {len(d.boxes)} boxes, {len(d.lines)} connectors, '
          f'{len(problems)} problem(s)')
    for p in problems:
        print('    ', p)
    return problems


# --------------------------------------------------------------------------
# Figure 1
# --------------------------------------------------------------------------
def fig_internal():
    LAST_DRAFTS.clear()
    fig, ax = new_canvas(20.0, 13.6)
    d = Draft(ax, s=1.0)
    draw_bit0(d)

    ax.text(10.0, 13.05, 'bit0 —— 单个 HBY 位（11 状态）的内部流程机制',
            ha='center', va='center', fontsize=17, color=DEEP_BLUE,
            fontweight='bold')
    ax.text(10.0, 12.56,
            'HBY 早期基础表参数　·　每条箭头 = model_hzh.py:39-63 中的一项 RHS',
            ha='center', va='center', fontsize=10.5, color=BLUE_GRAY)

    y0 = 4.15
    ax.plot([0.6, 19.4], [y0 + 0.44, y0 + 0.44], color=BLUE_GRAY, lw=1.2,
            ls=(0, (5, 4)))
    ax.text(0.85, y0 + 0.20, '图例：抽象序号 → 真实生化元件', ha='left',
            va='top', fontsize=12.5, color=DEEP_BLUE, fontweight='bold')

    col1 = [
        'a — M_I：整合酶 mRNA（ϕC31 整合酶转录本）',
        'b — I_u：未成熟整合酶（翻译产物）',
        'c — I：成熟游离整合酶（重组酶本体）',
        'd — M_T：BM3R1 / Rep 的 mRNA',
        'e — T_u：未成熟 BM3R1',
        'f — T：成熟 BM3R1 / Rep（阻遏蛋白）',
        'g — M_R：RDF 的 mRNA',
        'h — R_u：未成熟 RDF',
        'i — R：成熟游离 RDF（决定重组方向）',
        'j — C：Int·RDF 复合物（显式结合态）',
        'k — S：LR 构象比例；PB = 1 − S',
    ]
    col2 = [
        '① 转录：由上游 C31 启动子活性驱动（本图中是外源）',
        '② 翻译：M_I → I_u',
        '③ 成熟：I_u → I（k_mat）',
        '④ Rep 表达链：由 PB 驱动，源 α_rep·(1−S)',
        '⑤ RDF 表达链：由 LR 驱动，源 α_rdf·S·(1−H(T))',
        '⑥ T 抑制 RDF 表达：K_rep = 0.85，n_rep = 3.9',
        '⑦ 可逆结合：I + R 生成 C；k_on = 0.1，k_off = 1，δ_C = 1',
        '⑧ 正向重组：v_f = k_fwd·H(I;K_D_int,2)·K_inh/(K_inh+R)',
        '⑨ 游离 RDF 抑制正向 —— 就是 ⑧ 分母里的 K_inh/(K_inh+R)',
        '⑩ 反向重组：v_r = k_rev·H(C;K_complex,2)',
        '(11)(12) k 的两条输出：PB 与 LR 分别驱动两条表达链',
    ]
    col3 = [
        '开关（本图的核心）：',
        'dS/dt = v_f·(1−S) − v_r·S',
        '· 正向 → LR（S 上升）；反向 → PB（S 下降）',
        '· S + PB = 1；S 是 DNA 构象比例，不是浓度',
        '',
        'K_complex = q·K_D_comp，q = k_on/(k_off+δ_C) = 0.05',
        '· 因此 K_complex = 0.05 × 1.2 = 0.06 a.u.',
        '',
        '图上未逐条画的项：每个分子各自按 λ（mRNA 半衰期 2 min）',
        '或 γ（蛋白清除）衰减；C 另有 δ_C 的复合物降解。',
        '→ 本图只画 bit0；上游与下游 carry 见 fig_bit0_wired。',
    ]
    for col, x in ((col1, 0.95), (col2, 7.15), (col3, 14.05)):
        for i, line in enumerate(col):
            ax.text(x, y0 - 0.22 - i * 0.335, line, ha='left', va='center',
                    fontsize=9.7, color=INK)
    ax.add_patch(plt.Rectangle((13.90, y0 - 1.32), 5.5, 1.36, fill=False,
                               ec=ORANGE, lw=1.4, ls=(0, (4, 3)), zorder=6))

    check(d, 'fig_bit0_internal')
    return save_pair(fig, 'fig_bit0_internal')


# --------------------------------------------------------------------------
# Figure 2
# --------------------------------------------------------------------------
def fig_wired():
    LAST_DRAFTS.clear()
    fig, ax = new_canvas(32.0, 18.6)
    d = Draft(ax, s=1.0)
    OX, OY = 8.60, 4.60

    def container(x, y, w, h, ec, title, tx=None, ty=None, fs=10.2):
        ax.add_patch(FancyBboxPatch((x, y), w, h,
                                    boxstyle='round,pad=0,rounding_size=0.25',
                                    fc='none', ec=ec, lw=1.6, zorder=1))
        ax.text(tx if tx is not None else x + 0.28,
                ty if ty is not None else y + h - 0.30, title, fontsize=fs,
                color=ec, ha='left', va='center', zorder=2)

    # --- model scope ------------------------------------------------------
    ax.add_patch(FancyBboxPatch((7.90, 7.60), 24.0, 10.20,
                                boxstyle='round,pad=0,rounding_size=0.25',
                                fc='none', ec=BLUE_GRAY, lw=1.8,
                                ls=(0, (7, 5)), zorder=1))
    ax.text(8.25, 17.50, '34-state hybrid model（本次工作的模型范围）',
            fontsize=11.5, color=BLUE_GRAY, ha='left', va='center', zorder=2)

    # --- upstream ---------------------------------------------------------
    container(0.45, 12.05, 7.05, 4.85, BLUE_GRAY,
              '上游：韩亚轩 v53d（外生输入，不在 34 态内）')
    for key, y, t in (('A', 15.80, 'A  TetR'), ('B', 14.70, 'B  CI'),
                      ('C', 13.60, 'C  LacI')):
        d.box(f'up.{key}', 3.85, y, 2.80, 0.86, t, fc=LIGHT_ORANGE, ec=ORANGE,
              fs=10.0)
    d.repress('up.ab', (2.22, 15.80), (2.22, 14.70), color=ORANGE, rad=0.45)
    d.repress('up.bc', (2.22, 14.70), (2.22, 13.60), color=ORANGE, rad=0.45)
    d.repress('up.ca', (5.48, 13.60), (5.48, 15.80), color=ORANGE, rad=0.55)
    d.label(1.15, 14.70, '抑制环', fs=9.2, color=ORANGE)
    d.box('up.D', 3.85, 12.62, 3.60, 0.86, 'D  C31 启动子',
          fc=LIGHT_ORANGE, ec=ORANGE, fs=10.0)

    # --- bit0 -------------------------------------------------------------
    draw_bit0(d, ox=OX, oy=OY, compact=True)
    ax.text(OX + 6.70, 16.90, 'bit0（HBY，11 状态）', fontsize=12.5,
            color=DEEP_BLUE, ha='center', va='center', fontweight='bold')

    d.arrow('w1', (5.65, 12.62), (OX + 1.60, OY + 10.90), color=ORANGE,
            rad=-0.16)
    d.label(8.30, 12.55, '① 转录\n(C31 启动子活性)', fs=9.6, color=ORANGE,
            ha='left')

    # --- carry0 -----------------------------------------------------------
    container(21.00, 11.20, 6.40, 6.40, TEAL,
              'carry0（ZMH 现行 Python 代码）')
    d.box('c0.E', 23.20, 16.30, 3.60, 0.88, 'E  A0（激活臂）', fc=LIGHT_TEAL,
          ec=TEAL, fs=10.0)
    d.box('c0.F', 23.20, 14.72, 3.60, 0.88, 'F  F0（抑制臂）', fc=LIGHT_TEAL,
          ec=TEAL, fs=10.0)
    d.circle('c0.g', 26.10, 13.85, 0.42, '×', fc=WHITE, ec=TEAL, fs=14)
    # E is a FORK: the activator drives the gate directly AND, with a delay,
    # builds its own repressor F.  Without this edge F had NO incoming arrow at
    # all, so the incoherent feed-forward loop was invisible and the panel read
    # as a plain AND of two independent inputs.  The E-F gap is deliberately
    # wide enough for the arrowhead to read at page scale.
    d.arrow('c0.ef', (23.20, 15.86), (23.20, 15.16), color=TEAL)
    d.label(23.45, 15.51, '延迟臂', fs=9.0, color=TEAL, ha='left')
    d.arrow('c0.eg', (25.00, 16.20), (26.10, 14.27), color=TEAL)
    d.repress('c0.fg', (24.80, 14.30), (25.68, 13.85), color=RISK_ORANGE)
    d.label(22.80, 13.28,
            'g0 = H(A0;1.0,2)·G(F0;2.6454,5.4332)\n无时钟', fs=9.0, color=TEAL)
    d.box('c0.G', 22.80, 12.00, 3.40, 0.86, 'G  Int1', fc=LIGHT_TEAL, ec=TEAL,
          fs=10.0)
    d.arrow('c0.gG', (25.90, 13.48), (24.10, 12.43), color=TEAL)

    # The PB0 net label needs a LEAD-IN ARROW into E as well as a name.  Without
    # one, box E had in-degree 0 by drawing: "< PB0" was only a floating text
    # and "(11) PB" (from draw_bit0) sat right next to it, so the reader saw two
    # PB labels and no way in.  Solid, not dashed like w3 -> ck: the legend
    # reserves dashes for links NOT expanded in this figure, and this one is
    # listed under 【本图画的连接】.
    d.arrow('w4', (20.62, 16.30), (21.42, 16.30), color=DEEP_BLUE, lw=1.8,
            ms=13)

    # --- bit1 -------------------------------------------------------------
    container(27.60, 11.80, 4.10, 5.80, DEEP_BLUE,
              'bit1（ZMH 约化位）', tx=27.88)
    d.box('b1.blk', 29.65, 14.60, 3.55, 3.60,
          'H\npb1 / I1\nT1 / RDF1\n\n4 状态\n无 C 复合物',
          fc=LIGHT_BLUE, ec=DEEP_BLUE, fs=9.4)
    d.arrow('c0.GH', (24.50, 12.00), (27.88, 13.10), color=TEAL)

    # --- clock tap (net label, not a drawn long line) ---------------------
    d.box('ck', 24.30, 8.60, 7.00, 1.05,
          'I  —  clock 门   H(b0_I ; K = 0.3, n = 2)', fc=LIGHT_TEAL, ec=TEAL,
          fs=10.2)
    d.arrow('w3', (19.70, 8.60), (20.80, 8.60), color=TEAL, ls=(0, (5, 4)),
            lw=1.6, ms=12)

    # --- net labels -------------------------------------------------------
    d.label(OX + 8.90, OY + 8.70, 'PB0 >', fs=10.6, color=DEEP_BLUE,
            weight='bold', ha='left')
    d.label(20.90, 16.62, '< PB0', fs=10.6, color=DEEP_BLUE, weight='bold',
            ha='right')
    d.label(OX + 1.60, OY + 6.55, 'b0_I >', fs=10.6, color=TEAL,
            weight='bold', ha='right')
    d.label(19.60, 8.95, '< b0_I', fs=10.6, color=TEAL, weight='bold',
            ha='right')

    # --- title + legend ---------------------------------------------------
    ax.text(16.0, 18.20, 'bit0 接上上游与下游 carry 之后的完整流程机制',
            ha='center', va='center', fontsize=17, color=DEEP_BLUE,
            fontweight='bold')
    ax.plot([0.6, 31.4], [4.55, 4.55], color=BLUE_GRAY, lw=1.2,
            ls=(0, (5, 4)))
    ax.text(0.85, 4.44, '图例：抽象序号 → 真实生化元件', ha='left', va='top',
            fontsize=12.5, color=DEEP_BLUE, fontweight='bold')

    col1 = [
        '【上游 · 韩亚轩 v53d，外生】',
        'A — TetR：振荡器阻遏蛋白（总蛋白）',
        'B — CI',
        'C — LacI',
        'D — C31 启动子（PLtetO1 型，输出活性 h31(t)）',
        '',
        '【bit0 · HBY 早期基础表】',
        'a - k 同 fig_bit0_internal 的 11 个状态',
        '',
        '【carry0 + bit1 · ZMH 现行 Python 代码】',
        'E — A0：进位激活臂（由 bit0 的 PB 驱动）',
        'F — F0：进位抑制臂（曾代码里的变量名是 R0）',
        'G — Int1：bit1 的整合酶（由 g0 门产生）',
        'H — bit1 约化模块：pb1 / I1 / T1 / RDF1',
    ]
    col2 = [
        '【本图画的连接】',
        '① 上游 C31 启动子活性 → bit0 的 a（M_I 转录）',
        '   注：翻译系数 0.5/min 也取自上游，但反应 ② 在 bit0 内部',
        '',
        '(11) bit0 的 k → carry0 的 E（PB0 驱动 A0）',
        '   · dA0 = α_A0·(1 − b0_S) − γ_A0·A0',
        '',
        'carry0 内部：A0 激活、F0 抑制，二者在 × 节点处相乘',
        '   · g0 = H(A0;1.0,2)·G(F0;2.6454,5.4332)',
        '   · 与 carry1 不同，g0 只有两项、没有时钟',
        '',
        'G → H：Int1 进入 bit1，驱动 bit1 的 pb1 翻转',
        '   · v_f = k_fwd·pb1·H(I1;1.8,4)·Kinh/(Kinh+RDF1)',
    ]
    col3 = [
        '【clock 回绕：本图只画接口】',
        'I — clock 门 H(b0_I; K = 0.3, n = 2)',
        '   · 它读的是 bit0 的成熟游离整合酶 c（= b0_I）',
        '   · 这条线一路绕到第三级，参与 carry1 的三因子 AND：',
        '     g1 = H(A1;1.2,6) · G(F1;0.4,4) · clock',
        '   · 本图不展开 carry1 与 bit2',
        '',
        '【虚线 = 尚未在本图展开的连接】',
        '【net label  > / <  = 同名标签表示同一条连线】',
        '   用于避免跨越整张图的长线',
        '',
        '【本图未画的项】各位的清除/衰减项，以及 carry1、bit2',
    ]
    for col, x in ((col1, 0.95), (col2, 10.60), (col3, 21.30)):
        for i, line in enumerate(col):
            ax.text(x, 4.05 - i * 0.295, line, ha='left', va='center',
                    fontsize=9.5, color=INK)

    check(d, 'fig_bit0_wired')
    return save_pair(fig, 'fig_bit0_wired')


def main():
    for o in (fig_internal(), fig_wired()):
        print('wrote', o)


if __name__ == '__main__':
    main()
