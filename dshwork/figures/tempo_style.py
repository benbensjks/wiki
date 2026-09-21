"""TEMPO schematic style: palette, line conventions and drawing primitives.

Colours and line styles follow the project-wide visual specification already
used by the Oscillator chapter, so every mechanism figure on the Wiki belongs
to one visual family:

  * repression        -> solid line ending in a T-bar
  * activation        -> solid arrow
  * reversible binding-> double-headed arrow
  * AND / multiplication -> the factors meet at a circle marked with a cross
  * threshold         -> dashed TEMPO-orange horizontal line
  * model scope       -> blue-gray dashed rounded box
  * not yet modelled  -> blue-gray dashed hollow arrow

ALL LABELS ARE ENGLISH ONLY.  The figure briefs require copy-ready English
labels, and English avoids any CJK font dependency in headless rendering.
"""
from __future__ import annotations

import matplotlib
matplotlib.use('Agg')            # headless; must precede pyplot import
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch
from matplotlib.lines import Line2D

# ----------------------------------------------------------------- palette
DEEP_BLUE = '#304B53'      # structure, DNA state, repression
ORANGE = '#D29144'         # timing, thresholds, risk and emphasis
TEAL = '#4F9194'           # recombinase, results, already-passed paths
INK = '#33383A'            # body and axis text
WHITE = '#FFFFFF'
LIGHT_GRAY = '#F3F5F6'
LIGHT_BLUE = '#DCE5E7'
LIGHT_ORANGE = '#F2E2CF'
LIGHT_TEAL = '#DCEBEC'
BLUE_GRAY = '#8CA0A5'
RISK_ORANGE = '#B86F35'

FONT = 'DejaVu Sans'       # ships with matplotlib; no CJK needed


def new_canvas(w=12.0, h=8.0, xlim=None, ylim=None, facecolor=WHITE):
    """Create a blank figure with hidden axes in data coordinates."""
    fig, ax = plt.subplots(figsize=(w, h), dpi=200)
    fig.patch.set_facecolor(facecolor)
    ax.set_xlim(*(xlim if xlim else (0, w)))
    ax.set_ylim(*(ylim if ylim else (0, h)))
    ax.set_aspect('equal', adjustable='box')
    ax.axis('off')
    return fig, ax


# -------------------------------------------------------------- primitives
def box(ax, cx, cy, w, h, text, fc=LIGHT_GRAY, ec=DEEP_BLUE, tc=INK,
        fontsize=10, bold=False, lw=1.6, ls='-', rounding=0.14, zorder=3,
        linespacing=1.35, alpha=1.0):
    """Rounded rectangle with centred text."""
    patch = FancyBboxPatch((cx - w / 2, cy - h / 2), w, h,
                           boxstyle=f'round,pad=0,rounding_size={rounding}',
                           fc=fc, ec=ec, lw=lw, ls=ls, alpha=alpha, zorder=zorder)
    ax.add_patch(patch)
    ax.text(cx, cy, text, ha='center', va='center', fontsize=fontsize,
            color=tc, fontweight='bold' if bold else 'normal',
            zorder=zorder + 1, linespacing=linespacing)
    return patch


def arrow(ax, p0, p1, color=DEEP_BLUE, lw=1.8, ms=15, ls='-', rad=0.0,
          zorder=2, alpha=1.0):
    """Solid arrow from p0 to p1 (activation / production)."""
    a = FancyArrowPatch(p0, p1, arrowstyle='-|>', mutation_scale=ms, lw=lw,
                        color=color, linestyle=ls, zorder=zorder, alpha=alpha,
                        shrinkA=1.5, shrinkB=1.5,
                        connectionstyle=f'arc3,rad={rad}')
    ax.add_patch(a)
    return a


def repression(ax, p0, p1, color=DEEP_BLUE, lw=1.8, bar=0.16, rad=0.0, zorder=2):
    """Repression: a line that ends in a T-bar instead of an arrow head."""
    dx, dy = p1[0] - p0[0], p1[1] - p0[1]
    length = (dx * dx + dy * dy) ** 0.5
    ux, uy = dx / length, dy / length
    ex, ey = p1[0] - ux * 0.16, p1[1] - uy * 0.16
    ax.add_patch(FancyArrowPatch(p0, (ex, ey), arrowstyle='-', lw=lw,
                                 color=color, zorder=zorder, shrinkA=1.5,
                                 shrinkB=0.0,
                                 connectionstyle=f'arc3,rad={rad}'))
    px, py = -uy, ux
    ax.plot([ex - px * bar, ex + px * bar], [ey - py * bar, ey + py * bar],
            color=color, lw=lw + 0.6, solid_capstyle='round', zorder=zorder + 1)


def binding(ax, p0, p1, color=INK, lw=1.5, ms=11, rad=0.0, ls='-', zorder=2):
    """Reversible binding: double-headed arrow."""
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle='<|-|>', mutation_scale=ms,
                                 lw=lw, color=color, linestyle=ls, zorder=zorder,
                                 shrinkA=1.5, shrinkB=1.5,
                                 connectionstyle=f'arc3,rad={rad}'))


def hollow_arrow(ax, p0, p1, color=BLUE_GRAY, lw=1.5, rad=0.0, zorder=2):
    """Dashed hollow arrow for 'not yet modelled' couplings."""
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle='-|>', mutation_scale=13,
                                 lw=lw, color=color, linestyle=(0, (5, 4)),
                                 zorder=zorder, shrinkA=2, shrinkB=2,
                                 connectionstyle=f'arc3,rad={rad}'))


def and_node(ax, cx, cy, r=0.34, fc=LIGHT_TEAL, ec=DEEP_BLUE, glyph='×'):
    """AND node: the factors meet here and are multiplied."""
    ax.add_patch(Circle((cx, cy), r, fc=fc, ec=ec, lw=1.7, zorder=5))
    ax.text(cx, cy, glyph, ha='center', va='center', fontsize=13,
            color=ec, zorder=6)


def hline(ax, y, x0, x1, color=ORANGE, lw=1.6, ls='--', label=None,
          label_x=None, label_va='bottom', fontsize=9, zorder=1):
    """Threshold reference line."""
    ax.plot([x0, x1], [y, y], color=color, lw=lw, ls=ls, zorder=zorder)
    if label:
        ax.text(label_x if label_x is not None else (x0 + x1) / 2, y,
                label, ha='center', va=label_va, fontsize=fontsize,
                color=color, zorder=zorder + 1)


def scope(ax, x0, y0, w, h, label, color=BLUE_GRAY, lw=1.6, fontsize=9):
    """Model-scope boundary box."""
    ax.add_patch(FancyBboxPatch((x0, y0), w, h,
                                boxstyle='round,pad=0,rounding_size=0.22',
                                fc='none', ec=color, lw=lw,
                                ls=(0, (7, 5)), zorder=1))
    ax.text(x0 + 0.18, y0 + h - 0.26, label, fontsize=fontsize, color=color,
            ha='left', va='top', zorder=2)


def note(ax, x, y, text, color=INK, fontsize=9.5, ha='left', va='top',
         style='italic', wrap_width=None):
    """Small annotation text."""
    ax.text(x, y, text, ha=ha, va=va, fontsize=fontsize, color=color,
            style=style, zorder=7, linespacing=1.4)


def legend_lines(ax, x, y, items, fontsize=9, dy=0.42):
    """Draw a small line-convention legend.  items: (kind, label)."""
    for i, (kind, label) in enumerate(items):
        yy = y - i * dy
        if kind == 'repression':
            ax.plot([x, x + 0.55], [yy, yy], color=DEEP_BLUE, lw=1.8, zorder=7)
            ax.plot([x + 0.55, x + 0.55], [yy - 0.09, yy + 0.09],
                    color=DEEP_BLUE, lw=2.4, zorder=7)
        elif kind == 'activation':
            ax.add_patch(FancyArrowPatch((x, yy), (x + 0.55, yy),
                                         arrowstyle='-|>', mutation_scale=12,
                                         lw=1.8, color=TEAL, zorder=7))
        elif kind == 'binding':
            ax.add_patch(FancyArrowPatch((x, yy), (x + 0.55, yy),
                                         arrowstyle='<|-|>', mutation_scale=10,
                                         lw=1.5, color=INK, zorder=7))
        elif kind == 'threshold':
            ax.plot([x, x + 0.55], [yy, yy], color=ORANGE, lw=1.6,
                    ls='--', zorder=7)
        elif kind == 'scope':
            ax.add_patch(FancyBboxPatch((x, yy - 0.12), 0.55, 0.24,
                                        boxstyle='round,pad=0,rounding_size=0.06',
                                        fc='none', ec=BLUE_GRAY, lw=1.4,
                                        ls=(0, (5, 4)), zorder=7))
        elif kind == 'hollow':
            ax.add_patch(FancyArrowPatch((x, yy), (x + 0.55, yy),
                                         arrowstyle='-|>', mutation_scale=11,
                                         lw=1.4, color=BLUE_GRAY,
                                         linestyle=(0, (4, 3)), zorder=7))
        ax.text(x + 0.72, yy, label, ha='left', va='center',
                fontsize=fontsize, color=INK, zorder=7)


def save(fig, path):
    import pathlib
    path = pathlib.Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=220, bbox_inches='tight', facecolor=fig.get_facecolor())
    plt.close(fig)
    return path
