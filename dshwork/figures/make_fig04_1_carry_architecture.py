r"""Figure 4-1 -- Carry architecture   (dshwork, per the 精简交付版 brief)

Brief: `dshwork/精简交付版/02_作图要求精简版.md` section "Figure 4-1 -- Carry
architecture".  Requirements:

  * top strip: the whole chain `Clock -> Bit 0 -> Carry 0 -> Bit 1 -> Carry 1 -> Bit 2`;
  * below, the two carry stages SIDE BY SIDE:
        Carry 0: PB0 -> A0 -> Int1 expression
                         \-> F0 -| Int1 expression
        Carry 1: PB1 -> A1 -> Int2 expression
                         \-> F1 -| Int2 expression
                      A1 -| its own production
                      Int0 clock gate -> AND -> Int2 expression
  * ONE inset explaining `mRNA -> immature -> mature`, instead of repeating the full
    expression chain next to every arrow;
  * MUST be right: (a) only stage 2 has the clock AND gate; (b) F represses the
    OUTPUT, and the gate multiplies the activation and the DE-REPRESSED RESPONSE --
    it must NOT read as "A concentration x F concentration"; (c) `n_A1_gate=6` is
    labelled only on the A1 -> Int2 gate, while A1 -> F1 production keeps 4;
  * labels verbatim: `Carry 0`, `Carry 1`, `Delayed repression`,
    `Negative autoregulation`, `Clock gate`, `Next-bit integrase`,
    `Modelled interface; molecular implementation pending`.

Everything drawn is taken from the model source and re-asserted at run time against
the frozen artefacts, so the picture cannot silently drift from the model:

    model.py:118  g0 = act(A0, K_A[0], n_A[0]) * rep(F0, K_F[0], n_F[0])
    model.py:119  clock = act(y[8], clock_K_au, clock_n)      # mature free Int0
    model.py:120  g1 = act(A1, K_A[1], n_A[1]) * rep(F1, ...) * clock
    model.py:129  sources = (None, alpha_Int[0]*g0, alpha_Int[1]*g1)
    model.py:157  auto = rep(A, K_auto1, n_auto1) if j == 1 else 1.0
    model.py:158  d[A_j] = alpha_A[j] * PB_j * auto - loss*A
    model.py:159  d[F_j] = alpha_F[j] * act(A_j, K_A[j], n_A[j]) - loss*F

  => A is driven by PB; F is driven by A and represses the OUTPUT (the integrase
     expression); the gate is a product of HILL RESPONSES; the clock enters only g1;
     negative autoregulation exists only for A1.

Self-checks (the agent has no image input): labels verbatim; the two stage
asymmetries; each F T-bar lands on the OUTPUT box; PB->A and A->F directions; the
chain order left to right; the two exponents cross-checked against the frozen
working point; no label over another label, over a box, or over any straight arrow /
connector; everything inside the canvas; the SVG carries every label as <text>.

Run:  & 'D:\aconade\python.exe' .\make_fig04_1_carry_architecture.py
Out:  ./out/fig04_1_carry_architecture.png / .svg
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                                          # noqa: E402
import numpy as np                                                       # noqa: E402

from tempo_style import (BLUE_GRAY, DEEP_BLUE, INK, LIGHT_BLUE, LIGHT_GRAY,
                         LIGHT_ORANGE, LIGHT_TEAL, ORANGE, RISK_ORANGE, TEAL,
                         and_node, arrow, box, new_canvas, note, repression,
                         scope)                                          # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE / 'out'
STEM = 'fig04_1_carry_architecture'
FR = HERE.parents[1] / 'final_reconstruction'
PROFILE = FR / 'plausibility' / 'threebit51_selected_v1.json'

SVG_FONTTYPE = 'none'
REQUIRED_LABELS = ('Carry 0', 'Carry 1', 'Delayed repression',
                   'Negative autoregulation', 'Clock gate', 'Next-bit integrase',
                   'Modelled interface; molecular implementation pending')
GATE_LABEL = 'n_A1_gate=6'          # the brief writes it without spaces
F1_LABEL = 'n=4'
CLOCK_LABELS = ('Clock', 'Int0 (clock)', 'Clock gate')
CHAIN = ('Clock', 'Bit 0', 'Carry 0', 'Bit 1', 'Carry 1', 'Bit 2')
CHAIN_Y = 9.60
BLK0 = (0.30, 3.90, 7.90, 8.60)     # x0, y0, x1, y1
BLK1 = (8.30, 3.90, 16.60, 8.60)


def facts() -> dict:
    """The two exponents, read from the frozen artefacts (never hard-coded)."""
    prof = json.loads(PROFILE.read_text(encoding='utf-8'))
    n_gate = prof['frozen_parameters']['n_A1_gate']
    sys.path.insert(0, str(FR))
    import model as M                      # the ZENG table value actually in force
    return dict(n_A1_gate=float(n_gate), n_A_table=float(M.ZENG['n_A'][1]))


def draw_repression(ax, p0, p1, **kw):
    """Repression T-bar coordinates.

    `tempo_style.repression()` returns None, so capture the Line2D it draws as the
    bar -- that way the check reads the artist that was really added.
    """
    k = len(ax.lines)
    repression(ax, p0, p1, **kw)
    bar = ax.lines[k]
    return (np.asarray(bar.get_xdata(), dtype=float),
            np.asarray(bar.get_ydata(), dtype=float))


def draw():
    fig, ax = new_canvas(16.8, 10.6)
    S = {}

    # ------------------------------------------------------- top: whole chain
    xs = [2.10 + 2.60 * i for i in range(6)]
    for i, (name, cx) in enumerate(zip(CHAIN, xs)):
        fc, ec = ((LIGHT_ORANGE, ORANGE) if name.startswith('Carry')
                  else (LIGHT_TEAL, TEAL) if name == 'Clock'
                  else (LIGHT_BLUE, DEEP_BLUE))
        box(ax, cx, CHAIN_Y, 2.00, 0.72, name, fc=fc, ec=ec, fontsize=11.5,
            bold=True)
        if i:
            arrow(ax, (xs[i - 1] + 1.00, CHAIN_Y), (cx - 1.00, CHAIN_Y),
                  color=DEEP_BLUE, lw=2.0)

    # ======================================= stage 0 (left): no clock, no autoreg.
    scope(ax, BLK0[0], BLK0[1], BLK0[2] - BLK0[0], BLK0[3] - BLK0[1], 'Carry 0',
          color=BLUE_GRAY, fontsize=11.5)
    S['pb0'] = box(ax, 2.60, 7.70, 1.50, 0.72, 'PB0', fc=LIGHT_GRAY, ec=DEEP_BLUE,
                   fontsize=11)
    S['a0'] = box(ax, 4.50, 7.70, 1.40, 0.72, 'A0', fc=LIGHT_TEAL, ec=TEAL,
                  fontsize=11.5, bold=True)
    S['f0'] = box(ax, 4.50, 5.90, 1.40, 0.72, 'F0', fc=LIGHT_TEAL, ec=TEAL,
                  fontsize=11.5, bold=True)
    S['out0'] = box(ax, 6.20, 5.30, 2.60, 0.80, 'Int1 expression', fc=LIGHT_ORANGE,
                    ec=ORANGE, fontsize=11)
    S['pb0_a0'] = arrow(ax, (3.35, 7.70), (3.80, 7.70), color=TEAL, lw=2.0)
    S['a0_f0'] = arrow(ax, (4.50, 7.34), (4.50, 6.26), color=TEAL, lw=2.0)
    ax.plot([5.20, 6.20], [7.70, 7.70], color=TEAL, lw=2.0)
    S['a0_out'] = arrow(ax, (6.20, 7.70), (6.20, 5.72), color=TEAL, lw=2.0)
    note(ax, 2.95, 6.80, 'Delayed repression', fontsize=9, color=INK, ha='center',
         va='center', style='normal')
    ax.plot([4.50, 4.50], [5.54, 4.50], color=DEEP_BLUE, lw=1.9)
    ax.plot([4.50, 6.20], [4.50, 4.50], color=DEEP_BLUE, lw=1.9)
    S['f0_bar'] = draw_repression(ax, (6.20, 4.50), (6.20, 4.90), color=DEEP_BLUE,
                                  lw=1.9)
    note(ax, 6.20, 4.16, 'Next-bit integrase', fontsize=9, color=ORANGE,
         ha='center', va='center', style='normal')

    # =============================== stage 1 (right): clock AND + auto-repression
    scope(ax, BLK1[0], BLK1[1], BLK1[2] - BLK1[0], BLK1[3] - BLK1[1], 'Carry 1',
          color=BLUE_GRAY, fontsize=11.5)
    S['pb1'] = box(ax, 9.80, 7.70, 1.50, 0.72, 'PB1', fc=LIGHT_GRAY, ec=DEEP_BLUE,
                   fontsize=11)
    S['a1'] = box(ax, 11.60, 7.70, 1.40, 0.72, 'A1', fc=LIGHT_TEAL, ec=TEAL,
                  fontsize=11.5, bold=True)
    S['f1'] = box(ax, 11.60, 5.90, 1.40, 0.72, 'F1', fc=LIGHT_TEAL, ec=TEAL,
                  fontsize=11.5, bold=True)
    S['out1'] = box(ax, 13.20, 5.30, 2.60, 0.80, 'Int2 expression', fc=LIGHT_ORANGE,
                    ec=ORANGE, fontsize=11)
    S['clk'] = box(ax, 15.70, 6.90, 1.80, 0.72, 'Int0 (clock)', fc=LIGHT_BLUE,
                   ec=DEEP_BLUE, fontsize=10)
    S['pb1_a1'] = arrow(ax, (10.55, 7.70), (10.90, 7.70), color=TEAL, lw=2.0)
    S['a1_f1'] = arrow(ax, (11.60, 7.34), (11.60, 6.26), color=TEAL, lw=2.0)
    # the clock AND gate: A1 (from above) and Int0 (from the right) meet at the node
    and_node(ax, 13.20, 6.90, r=0.30)
    ax.plot([12.30, 13.20], [7.70, 7.70], color=TEAL, lw=2.0)
    S['a1_and'] = arrow(ax, (13.20, 7.70), (13.20, 7.22), color=TEAL, lw=2.0)
    S['clk_and'] = arrow(ax, (14.80, 6.90), (13.52, 6.90), color=TEAL, lw=2.0)
    S['and_out'] = arrow(ax, (13.20, 6.60), (13.20, 5.72), color=TEAL, lw=2.0)
    note(ax, 14.15, 6.62, 'Clock gate', fontsize=9, color=TEAL, ha='center',
         va='center', style='normal')
    note(ax, 13.36, 7.46, GATE_LABEL, fontsize=8.8, color=RISK_ORANGE, ha='left',
         va='center', style='normal')
    note(ax, 12.40, 6.62, F1_LABEL, fontsize=8.8, color=RISK_ORANGE, ha='left',
         va='center', style='normal')
    # A1 represses its OWN production -- stage 1 only (model.py:157)
    draw_repression(ax, (11.95, 8.10), (10.72, 7.86), color=DEEP_BLUE, lw=1.8,
                    rad=0.45)
    note(ax, 11.05, 8.40, 'Negative autoregulation', fontsize=9, color=DEEP_BLUE,
         ha='center', va='center', style='normal')
    note(ax, 10.55, 6.85, 'Delayed repression', fontsize=9, color=INK, ha='center',
         va='center', style='normal')
    ax.plot([11.60, 11.60], [5.54, 4.50], color=DEEP_BLUE, lw=1.9)
    ax.plot([11.60, 13.20], [4.50, 4.50], color=DEEP_BLUE, lw=1.9)
    S['f1_bar'] = draw_repression(ax, (13.20, 4.50), (13.20, 4.90), color=DEEP_BLUE,
                                  lw=1.9)
    note(ax, 15.00, 5.30, 'Next-bit integrase', fontsize=9, color=ORANGE,
         ha='left', va='center', style='normal')

    # --------------------------------------------------------- inset: expression
    scope(ax, 0.30, 0.55, 7.30, 2.75, 'Expression chain', color=BLUE_GRAY,
          fontsize=10.5)
    for cx, name in zip((1.55, 3.55, 5.70), ('mRNA', 'immature', 'mature')):
        box(ax, cx, 2.30, 1.50, 0.70, name, fc=LIGHT_BLUE, ec=DEEP_BLUE,
            fontsize=10.5)
    arrow(ax, (2.30, 2.30), (2.80, 2.30), color=TEAL, lw=1.8)
    arrow(ax, (4.30, 2.30), (4.95, 2.30), color=TEAL, lw=1.8)
    note(ax, 0.60, 1.58,
         'Every "expression" node above is this one chain\n'
         '(transcription \u2192 maturation); it is not repeated per arrow.\n'
         'The delay of the repression arm comes from it.',
         fontsize=9, color=INK, ha='left', va='top', style='normal')

    # ------------------------------------------- gate definitions and caveats
    note(ax, 8.90, 2.95, 'Carry 0 gate  =  act(A0) \u00d7 (1 \u2212 act(F0))',
         fontsize=9.5, color=INK, ha='left', va='center', style='normal')
    note(ax, 8.90, 2.45,
         'Carry 1 gate  =  act(A1) \u00d7 (1 \u2212 act(F1)) \u00d7 clock',
         fontsize=9.5, color=INK, ha='left', va='center', style='normal')
    note(ax, 8.90, 2.00,
         'act() is a Hill response: the gate multiplies RESPONSES,\n'
         'not A and F concentrations.',
         fontsize=9, color=RISK_ORANGE, ha='left', va='top', style='normal')
    note(ax, 8.90, 1.12,
         'Carry 1 is the only stage with the clock factor, and the only one\n'
         'with negative autoregulation (model.py:120 and :157).',
         fontsize=9, color=INK, ha='left', va='top', style='normal')
    note(ax, 16.55, 1.35,
         'Modelled interface; molecular implementation pending',
         fontsize=9, color='#5a6b70', ha='right', va='center', style='normal')
    return fig, ax, S


def collect_texts(fig):
    return [t for t in fig.findobj(matplotlib.text.Text)
            if isinstance(t.get_text(), str) and t.get_text().strip()]


def straight_inks(ax, fig):
    """Straight drawn segments in display coords, from the artists themselves.

    Line2D pieces are always straight.  For arrows, the endpoints come from the
    patch; diagonal ones (the autoregulation arc) are skipped and counted, because a
    sampled straight segment would misrepresent a curve.
    """
    segs, skipped = [], 0
    for ln in ax.lines:
        pts = ax.transData.transform(np.column_stack([
            np.asarray(ln.get_xdata(), dtype=float),
            np.asarray(ln.get_ydata(), dtype=float)]))
        segs += [(pts[i], pts[i + 1]) for i in range(len(pts) - 1)]
    for p in ax.patches:
        if p.__class__.__name__ != 'FancyArrowPatch':
            continue
        pts = getattr(p, '_posA_posB', None)
        if not pts:
            skipped += 1
            continue
        a, b = np.asarray(pts[0], float), np.asarray(pts[1], float)
        if abs(a[0] - b[0]) > 1e-9 and abs(a[1] - b[1]) > 1e-9:
            skipped += 1                 # diagonal / curved connector
            continue
        segs.append(tuple(ax.transData.transform(np.vstack([a, b]))))
    return segs, skipped


def verify(fig, ax, S, svg_path) -> int:
    fails = []

    def check(cond, msg):
        print(('  PASS  ' if cond else '  FAIL  ') + msg)
        if not cond:
            fails.append(msg)

    texts = [t.get_text() for t in collect_texts(fig)]

    print('-- 1. required labels, verbatim')
    for lab in REQUIRED_LABELS:
        check(lab in texts, f'label present: {lab!r}')

    print('-- 2. the two stage asymmetries the brief insists on')
    # exact label matches: a substring test would also catch the gate-formula line,
    # which mentions the word clock
    clock_texts = [t for t in collect_texts(fig) if t.get_text() in CLOCK_LABELS]
    check(len(clock_texts) == 3,
          f'the three clock elements are present (got '
          f'{[t.get_text() for t in clock_texts]})')
    check(len([t for t in clock_texts if t.get_text() == 'Clock gate']) == 1,
          "'Clock gate' is drawn exactly once")
    for t in clock_texts:
        if t.get_text() == 'Clock':
            continue                     # the chain box, deliberately at the top
        x, y = t.get_position()
        check(not (BLK0[0] <= x <= BLK0[2] and BLK0[1] <= y <= BLK0[3]),
              f'no clock element inside Carry 0 (found {t.get_text()!r})')
    x, y = [t for t in clock_texts if t.get_text() == 'Clock gate'][0].get_position()
    check(BLK1[0] <= x <= BLK1[2] and BLK1[1] <= y <= BLK1[3],
          'the clock gate is inside Carry 1')
    autos = [t for t in collect_texts(fig)
             if t.get_text() == 'Negative autoregulation']
    check(len(autos) == 1,
          f"'Negative autoregulation' appears exactly once (got {len(autos)})")
    if autos:
        ax_, _ay = autos[0].get_position()
        check(BLK1[0] <= ax_ <= BLK1[2],
              'negative autoregulation is drawn in Carry 1 only (model.py:157)')
    check(len([t for t in collect_texts(fig)
               if t.get_text() == 'Delayed repression']) == 2,
          'both stages show the delayed repression arm')
    ands = [t for t in collect_texts(fig) if t.get_text() == '\u00d7']
    check(len(ands) == 1, f'exactly one AND node is drawn (got {len(ands)})')
    if ands:
        check(BLK1[0] <= ands[0].get_position()[0] <= BLK1[2],
              'the AND node is inside Carry 1')

    print("-- 3. each F represses the OUTPUT box, not A")
    for tag, (bx, by), out in (('F0', S['f0_bar'], S['out0']),
                               ('F1', S['f1_bar'], S['out1'])):
        b = out.get_bbox()
        ex, ey = float(np.mean(bx)), float(np.mean(by))
        check(b.x0 - 0.3 <= ex <= b.x1 + 0.3 and b.y0 - 0.3 <= ey <= b.y1 + 0.3,
              f"the {tag} T-bar sits on the output box edge "
              f"({tag} bar at ({ex:.2f}, {ey:.2f}))")
        fb = S[tag.lower()].get_bbox()
        check(not (fb.x0 <= ex <= fb.x1 and fb.y0 <= ey <= fb.y1),
              f'the {tag} T-bar is NOT on the {tag} box itself')

    print('-- 4. directions: PB drives A; A drives F and the output')
    for tag, ar, src, dst in (('PB0->A0', S['pb0_a0'], S['pb0'], S['a0']),
                              ('PB1->A1', S['pb1_a1'], S['pb1'], S['a1'])):
        a, b = ar._posA_posB
        sb, db = src.get_bbox(), dst.get_bbox()
        check(b[0] > a[0], f'{tag} points left to right')
        check(abs(a[0] - sb.x1) < 0.25 and abs(b[0] - db.x0) < 0.25,
              f'{tag} runs from the PB box to the A box')
    for tag, ar, src, dst in (('A0->F0', S['a0_f0'], S['a0'], S['f0']),
                              ('A1->F1', S['a1_f1'], S['a1'], S['f1'])):
        a, b = ar._posA_posB
        sb, db = src.get_bbox(), dst.get_bbox()
        check(b[1] < a[1], f'{tag} points downwards (A drives F production)')
        check(abs(a[1] - sb.y0) < 0.25 and abs(b[1] - db.y1) < 0.25,
              f'{tag} runs from the A box down to the F box')

    print("-- 5. the two exponents are the model's, not decoration")
    check(S['n_gate'] == 6.0,
          f"frozen working point n_A1_gate = {S['n_gate']:g} "
          f"(profile frozen_parameters)")
    check(S['n_table'] == 4.0,
          f"ZENG table exponent used by F production = {S['n_table']:g}")
    check(GATE_LABEL in texts, f'the gate arm carries {GATE_LABEL!r}')
    check(F1_LABEL in texts, f'the F1 production arm carries {F1_LABEL!r}')
    gl = [t for t in collect_texts(fig) if t.get_text() == GATE_LABEL][0]
    fl = [t for t in collect_texts(fig) if t.get_text() == F1_LABEL][0]
    check(gl.get_position()[0] > S['f1'].get_bbox().x1,
          f'{GATE_LABEL} is placed on the A1 -> gate arm')
    check(fl.get_position()[1] < S['a1'].get_bbox().y0,
          f'{F1_LABEL} is placed on the A1 -> F1 arm (below A1)')

    print('-- 6. the top chain is in the required left-to-right order')
    order = []
    for name in CHAIN:
        hits = [t.get_position()[0] for t in collect_texts(fig)
                if t.get_text() == name and abs(t.get_position()[1] - CHAIN_Y) < 0.5]
        if hits:
            order.append((min(hits), name))
    order.sort()
    check([n for _x, n in order] == list(CHAIN),
          f'the chain reads {[n for _x, n in order]}')

    print('-- 7. layout: no label over a label, a box, or a connector')
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()

    def clash(a, b, tol=0.5):
        inter = matplotlib.transforms.Bbox.intersection(a, b)
        return inter is not None and inter.width > tol and inter.height > tol

    patch_centres, box_bboxes = [], []
    for p in ax.patches:
        if p.__class__.__name__ != 'FancyBboxPatch' or p.get_width() <= 0.5:
            continue
        # A `scope()` block border is a CONTAINER: transparent face, dashed edge.
        # Annotations are supposed to live inside it, so it is not a "box" here.
        container = (float(p.get_facecolor()[3]) == 0.0
                     and p.get_linestyle() not in ('-', 'solid'))
        if container:
            continue
        bb = p.get_bbox()
        patch_centres.append(((bb.x0 + bb.x1) / 2, (bb.y0 + bb.y1) / 2))
        box_bboxes.append(p.get_window_extent(renderer=renderer))
    items = [(t, t.get_window_extent(renderer=renderer))
             for t in collect_texts(fig) if t.get_visible()]
    pairs = [(items[i][0].get_text()[:22], items[j][0].get_text()[:22])
             for i in range(len(items)) for j in range(i + 1, len(items))
             if clash(items[i][1], items[j][1])]
    check(not pairs, f'no two labels overlap (found {len(pairs)}: {pairs[:4]})')

    over_box = []
    for t, bb in items:
        cx, cy = t.get_position()
        if any(abs(cx - px) < 0.05 and abs(cy - py) < 0.05
               for px, py in patch_centres):
            continue                        # this text IS the label of that box
        for pbb in box_bboxes:
            if clash(pbb, bb, tol=2.0):
                over_box.append(t.get_text()[:26])
                break
    check(not over_box,
          f'no annotation sits on a box (found {len(over_box)}: {over_box[:4]})')

    segs, skipped = straight_inks(ax, fig)
    on_ink = []
    for t, bb in items:
        for a, b in segs:
            for f in np.linspace(0.0, 1.0, 40):
                pt = a + f * (b - a)
                if bb.x0 < pt[0] < bb.x1 and bb.y0 < pt[1] < bb.y1:
                    on_ink.append((t.get_text()[:26],))
                    break
            else:
                continue
            break
    check(not on_ink,
          f'no label sits on a connector (found {len(on_ink)}: {on_ink[:4]}) '
          f'[{skipped} diagonal connector(s) not sampled]')
    figbb = fig.get_window_extent(renderer=renderer)
    outside = [t.get_text() for t, bb in items
               if bb.x0 < figbb.x0 - 1 or bb.x1 > figbb.x1 + 1
               or bb.y0 < figbb.y0 - 1 or bb.y1 > figbb.y1 + 1]
    check(not outside, f'every label is inside the canvas (outside: {outside})')

    print('-- 8. the exported SVG carries the labels as real text')
    svg = svg_path.read_text(encoding='utf-8')
    check('<text' in svg, "SVG uses <text> (svg.fonttype = 'none')")
    for lab in REQUIRED_LABELS + (GATE_LABEL, F1_LABEL):
        esc = lab.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')
        check(esc in svg, f'SVG carries the label: {lab!r}')

    print()
    print(f'FAILURES: {len(fails)}')
    for f in fails:
        print('  -', f)
    return len(fails)


def main() -> int:
    # The gate formulas contain U+2212 (minus sign) and U+00D7; the Windows console
    # is GBK and would raise UnicodeEncodeError while printing a check message.
    try:
        sys.stdout.reconfigure(encoding='utf-8', errors='replace')
    except Exception:                                            # noqa: BLE001
        pass
    plt.rcParams['svg.fonttype'] = SVG_FONTTYPE
    fct = facts()
    print(f"from the frozen artefacts: n_A1_gate = {fct['n_A1_gate']:g} "
          f"(profile), ZENG n_A table = {fct['n_A_table']:g} (model.py)")
    fig, ax, S = draw()
    S['n_gate'] = fct['n_A1_gate']
    S['n_table'] = fct['n_A_table']
    OUT.mkdir(parents=True, exist_ok=True)
    png, svg = OUT / f'{STEM}.png', OUT / f'{STEM}.svg'
    for p in (png, svg):
        fig.savefig(p, dpi=200, bbox_inches='tight',
                    facecolor=fig.get_facecolor())
    fails = verify(fig, ax, S, svg)
    plt.close(fig)
    print(f'wrote {png}')
    print(f'wrote {svg}')
    return 1 if fails else 0


if __name__ == '__main__':
    sys.exit(main())
