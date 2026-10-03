r"""Figure 3-1 -- Single-bit memory   (dshwork, per the 精简交付版 brief)

Brief: `dshwork/精简交付版/02_作图要求精简版.md` section "Figure 3-1 -- Single-bit
memory".  Composition required there:

  * PB (0) and LR (1) DNA configurations side by side, LEFT and RIGHT;
  * the round trip between them labelled `Int` and `Int-RDF complex`;
  * ONE regulatory block below, containing exactly
        PB -> BM3R1 -| RDF expression
        LR -> RDF expression
        Int + RDF <=> Int-RDF complex
  * beside each state one sentence:
        PB side:  Low RDF favours forward switching
        LR side:  Accumulated RDF supports reverse switching
  * the complex is an EFFECTIVE STATE of the model, so it is a plain box, never a
    resolved molecular structure;
  * the forbidden statement -- "DNA has no memory / memory lives in the RDF pool" --
    must NOT appear anywhere in the figure.

Why this script verifies itself
-------------------------------
The agent drawing this has no image input, so "it looks right" cannot be checked by
looking.  Two checks replace that:

  1. TEXT  -- the figure is exported as SVG with `svg.fonttype = 'none'`, which
     writes real <text> elements instead of glyph outlines.  That is also what the
     art team wants (labels stay editable), and it lets us assert that every
     required label is literally present in the exported file, and that no
     forbidden wording is.  The PNG is written as the visual reference.
  2. LAYOUT -- every text artist's rendered bounding box is computed and checked for
     pairwise overlap and for falling outside the canvas -- the closest available
     substitute for noticing overlapping labels by eye.

State geometry is read back from the patches that were actually created, so
"PB is on the left" is checked on the drawing rather than re-stated.

Run:  & 'D:\aconade\python.exe' .\make_fig03_1_single_bit_memory.py
Out:  ./out/fig03_1_single_bit_memory.png   (visual reference)
      ./out/fig03_1_single_bit_memory.svg   (editable, for the art team)
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                                          # noqa: E402

from tempo_style import (BLUE_GRAY, DEEP_BLUE, INK, LIGHT_BLUE, LIGHT_ORANGE,
                         LIGHT_TEAL, ORANGE, TEAL, arrow, binding, box,
                         new_canvas, repression, scope)                  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE / 'out'
STEM = 'fig03_1_single_bit_memory'

# Labels required by the brief, verbatim (note the EN DASH in the complex name).
REQUIRED_LABELS = (
    'PB / 0',
    'LR / 1',
    'BM3R1',
    'RDF',
    'Int',
    'Int\u2013RDF complex',
    'Forward switching',
    'Reverse switching',
    'Low RDF favours forward switching',
    'Accumulated RDF supports reverse switching',
)
# The brief forbids "the DNA has no memory / memory lives in the RDF pool".  The
# whole word family is forbidden here, which is the strict reading of that rule.
FORBIDDEN = ('memor', 'epigenetic')

SVG_FONTTYPE = 'none'   # keep <text> in the SVG: editable AND checkable


def draw():
    """Draw the figure; return (fig, ax, shapes)."""
    fig, ax = new_canvas(13.0, 9.0)
    S = {}

    # ---------------------------------------------------------------- states
    S['pb'] = box(ax, 3.00, 7.40, 3.4, 1.40, 'PB / 0', fc=LIGHT_BLUE,
                  ec=DEEP_BLUE, fontsize=14, bold=True)
    S['lr'] = box(ax, 10.00, 7.40, 3.4, 1.40, 'LR / 1', fc=LIGHT_BLUE,
                  ec=DEEP_BLUE, fontsize=14, bold=True)

    # ------------------------------------------- the round trip between them
    S['fwd_arrow'] = arrow(ax, (4.75, 7.78), (8.25, 7.78), color=TEAL, lw=2.4,
                           ms=18)
    ax.text(6.50, 8.58, 'Int', ha='center', va='center', fontsize=11.5,
            color=TEAL, fontweight='bold', zorder=8)
    ax.text(6.50, 8.16, 'Forward switching', ha='center', va='center',
            fontsize=10, color=INK, zorder=8)

    S['rev_arrow'] = arrow(ax, (8.25, 7.02), (4.75, 7.02), color=DEEP_BLUE,
                           lw=2.4, ms=18)
    ax.text(6.50, 6.52, 'Int\u2013RDF complex', ha='center', va='center',
            fontsize=11, color=DEEP_BLUE, fontweight='bold', zorder=8)
    ax.text(6.50, 6.12, 'Reverse switching', ha='center', va='center',
            fontsize=10, color=INK, zorder=8)

    # ------------------------------------------------ one sentence per state
    # The two sentences are wide and there is only ~6.5 data units of free space
    # between the two vertical connectors, so they are staggered vertically rather
    # than set on one line: on one line they overlapped by ~1 inch (found by the
    # text-bbox check below, which is why that check exists).
    ax.text(3.45, 5.80, 'Low RDF favours forward switching', ha='left',
            va='center', fontsize=10, color=INK, zorder=8)
    ax.text(9.35, 5.28, 'Accumulated RDF supports reverse switching', ha='right',
            va='center', fontsize=10, color=INK, zorder=8)

    # ------------------------------------------- ONE regulatory block below
    scope(ax, 0.70, 0.50, 11.60, 4.40, '', color=BLUE_GRAY)

    # PB -> BM3R1
    S['bm3r1'] = box(ax, 3.20, 4.30, 2.60, 1.00, 'BM3R1', fc=LIGHT_ORANGE,
                     ec=ORANGE, fontsize=11.5, bold=True)
    S['pb_to_bm3r1'] = arrow(ax, (3.10, 6.68), (3.15, 4.82), color=TEAL, lw=2.0)

    # LR -> RDF expression -> RDF
    S['lr_to_rdf'] = arrow(ax, (9.60, 6.68), (9.60, 2.97), color=TEAL, lw=2.0)
    ax.text(9.85, 3.90, 'RDF expression', ha='left', va='center', fontsize=10,
            color=INK, zorder=8)
    S['rdf'] = box(ax, 9.60, 2.50, 2.60, 1.00, 'RDF', fc=LIGHT_ORANGE,
                   ec=ORANGE, fontsize=11.5, bold=True)

    # BM3R1 -| RDF expression: the T-bar crosses the vertical production arrow
    S['repression'] = repression(ax, (4.50, 4.30), (9.60, 4.30),
                                 color=DEEP_BLUE, lw=2.0)

    # Int + RDF  <=>  Int-RDF complex  (an effective state: a plain box, no structure)
    S['int_rdf'] = box(ax, 3.10, 1.30, 2.60, 0.90, 'Int + RDF', fc=LIGHT_TEAL,
                       ec=TEAL, fontsize=11)
    S['equilibrium'] = binding(ax, (4.50, 1.30), (6.62, 1.30), color=INK, lw=1.8,
                               ms=14)
    S['complex'] = box(ax, 8.40, 1.30, 3.40, 0.90, 'Int\u2013RDF complex',
                       fc=LIGHT_TEAL, ec=INK, fontsize=11)
    return fig, ax, S


def patch_centre(patch):
    """Centre of a rounded box, read back from the artist (not from the input)."""
    b = patch.get_bbox()
    return ((b.x0 + b.x1) / 2.0, (b.y0 + b.y1) / 2.0)


def arrow_endpoints(art):
    """Explicit endpoints of a FancyArrowPatch, if this matplotlib exposes them."""
    pts = getattr(art, '_posA_posB', None)
    return pts


def verify(fig, ax, S, svg_path: Path) -> int:
    fails = []

    def check(cond, msg):
        print(('  PASS  ' if cond else '  FAIL  ') + msg)
        if not cond:
            fails.append(msg)

    texts = [t.get_text() for t in ax.texts if t.get_text().strip()]

    print('-- 1. required labels, verbatim, as standalone text elements')
    for lab in REQUIRED_LABELS:
        check(lab in texts, f'label present: {lab!r}')

    print('-- 2. forbidden wording absent')
    for bad in FORBIDDEN:
        hits = sorted({t for t in texts if bad.lower() in t.lower()})
        check(not hits, f'no forbidden wording {bad!r} (hits: {hits})')

    print('-- 3. PB is left of LR; the round trip points the right way')
    pbx, pby = patch_centre(S['pb'])
    lrx, lry = patch_centre(S['lr'])
    check(pbx < lrx, f'PB card LEFT of LR card (x {pbx:.2f} < {lrx:.2f})')
    check(abs(pby - lry) < 1e-9, 'the two states are drawn at the same height')
    fwd, rev = arrow_endpoints(S['fwd_arrow']), arrow_endpoints(S['rev_arrow'])
    check(fwd is not None and fwd[1][0] > fwd[0][0],
          'forward arrow runs PB -> LR (left to right)')
    check(rev is not None and rev[1][0] < rev[0][0],
          'reverse arrow runs LR -> PB (right to left)')
    check(fwd is not None and fwd[0][0] > pbx and fwd[1][0] < lrx,
          'the forward arrow starts at PB and ends at LR')

    print('-- 4. BM3R1 represses the RDF-expression step, not something else')
    rdf_x = arrow_endpoints(S['lr_to_rdf'])[0][0]
    check(abs(rdf_x - 9.60) < 1e-9,
          f'RDF-expression arrow is the vertical line at x = {rdf_x:.2f}')
    check(patch_centre(S['bm3r1'])[0] < rdf_x,
          'BM3R1 sits to the LEFT of the arrow it represses')
    check(abs(patch_centre(S['rdf'])[0] - rdf_x) < 1e-9,
          'the RDF box is directly below the RDF-expression arrow')
    check(patch_centre(S['rdf'])[1] < patch_centre(S['bm3r1'])[1],
          'RDF sits BELOW BM3R1 (the repression acts on the production step)')

    print('-- 5. layout: no overlapping labels, nothing outside the canvas')
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    spans = [(t.get_text(), t.get_window_extent(renderer=renderer))
             for t in ax.texts if t.get_text().strip()]
    overlaps = []
    for i in range(len(spans)):
        for j in range(i + 1, len(spans)):
            (ta, ba), (tb, bb) = spans[i], spans[j]
            # Bbox.intersection is a staticmethod in this matplotlib version.
            inter = matplotlib.transforms.Bbox.intersection(ba, bb)
            if inter is not None and inter.width > 0.5 and inter.height > 0.5:
                overlaps.append((ta[:26], tb[:26], round(inter.width, 1),
                                 round(inter.height, 1)))
    check(not overlaps,
          f'no two labels overlap (found {len(overlaps)}: {overlaps[:3]})')
    axbb = ax.get_window_extent(renderer=renderer)
    outside = [t for t, bb in spans
               if bb.x0 < axbb.x0 - 1 or bb.x1 > axbb.x1 + 1
               or bb.y0 < axbb.y0 - 1 or bb.y1 > axbb.y1 + 1]
    check(not outside, f'every label is inside the canvas (outside: {outside})')

    print('-- 6. the exported SVG carries the labels as real text')
    svg = svg_path.read_text(encoding='utf-8')
    check('<text' in svg, "SVG uses <text> (svg.fonttype = 'none')")
    for lab in REQUIRED_LABELS:
        esc = (lab.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;'))
        check(esc in svg, f'SVG carries the label: {lab!r}')
    for bad in FORBIDDEN:
        check(bad.lower() not in svg.lower(),
              f'SVG has no forbidden wording: {bad!r}')

    print()
    print(f'FAILURES: {len(fails)}')
    for f in fails:
        print('  -', f)
    return len(fails)


def main() -> int:
    plt.rcParams['svg.fonttype'] = SVG_FONTTYPE
    fig, ax, S = draw()
    OUT.mkdir(parents=True, exist_ok=True)
    png, svg = OUT / f'{STEM}.png', OUT / f'{STEM}.svg'
    for p in (png, svg):
        fig.savefig(p, dpi=220, bbox_inches='tight',
                    facecolor=fig.get_facecolor())
    fails = verify(fig, ax, S, svg)
    plt.close(fig)
    print(f'wrote {png}')
    print(f'wrote {svg}')
    return 1 if fails else 0


if __name__ == '__main__':
    sys.exit(main())
