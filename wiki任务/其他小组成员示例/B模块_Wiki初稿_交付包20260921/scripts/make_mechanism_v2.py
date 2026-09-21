# -*- coding: utf-8 -*-
"""TEMPO B module - redesigned single-bit switch mechanism figure (v2).

Concept: a two-state cycle (PB <-> LR) that makes the state-dependent logic explicit:
  - which protein is produced in each state (PB -> BM3R1; LR -> RDF);
  - how full the RDF pool is at pulse arrival in each state;
  - which reaction the arriving clock pulse therefore drives.

All labels English (group convention). TEMPO palette.

Run:  python make_mechanism_v2.py
Out:  ../fig03_mechanism_v2_sketch.png
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle

DEEP = "#304B53"
ORANGE = "#D29144"
TEAL = "#4F9194"
LIGHT_BLUE = "#DCE5E7"
LIGHT_ORANGE = "#F2E2CF"
LIGHT_TEAL = "#DCEBEC"
GRAY = "#8CA0A5"
CHARCOAL = "#33383A"

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "fig03_mechanism_v2_sketch.png")


def rbox(ax, x, y, w, h, fc, ec, lw=1.8, rad=0.16, z=2):
    p = FancyBboxPatch((x, y), w, h,
                       boxstyle=f"round,pad=0.06,rounding_size={rad}",
                       facecolor=fc, edgecolor=ec, linewidth=lw, zorder=z)
    ax.add_patch(p)


def arrow(ax, p1, p2, color, lw=1.9, rad=0.0, style="-|>", mutation=15, z=3):
    a = FancyArrowPatch(p1, p2, arrowstyle=style, mutation_scale=mutation,
                        color=color, lw=lw,
                        connectionstyle=f"arc3,rad={rad}", zorder=z)
    ax.add_patch(a)


def tbar(ax, p1, p2, color, lw=1.9, z=3):
    a = FancyArrowPatch(p1, p2, arrowstyle="-", color=color, lw=lw, zorder=z)
    ax.add_patch(a)
    dx, dy = p2[0] - p1[0], p2[1] - p1[1]
    n = np.hypot(dx, dy)
    px, py = -dy / n, dx / n
    half = 0.16
    ax.plot([p2[0] - px * half, p2[0] + px * half],
            [p2[1] - py * half, p2[1] + py * half],
            color=color, lw=lw, zorder=z, solid_capstyle="butt")


def txt(ax, x, y, s, size=10, color=CHARCOAL, ha="center", va="center",
        weight="normal", style="normal", z=4):
    ax.text(x, y, s, fontsize=size, color=color, ha=ha, va=va,
            weight=weight, style=style, zorder=z)


def pulse_glyph(ax, x, y, w=0.55, h=0.22, color=TEAL, lw=1.8):
    xs = [x, x + 0.18 * w, x + 0.18 * w, x + 0.62 * w, x + 0.62 * w, x + w]
    ys = [y, y, y + h, y + h, y, y]
    ax.plot(xs, ys, color=color, lw=lw, zorder=4)


def dna_icon(ax, cx, cy, left_label, right_label, caption):
    w, h = 1.55, 0.62
    gap = 0.10
    rbox(ax, cx - w - gap / 2, cy - h / 2, w, h, LIGHT_BLUE, DEEP, lw=1.4, rad=0.10)
    rbox(ax, cx + gap / 2, cy - h / 2, w, h, LIGHT_TEAL, TEAL, lw=1.4, rad=0.10)
    txt(ax, cx - w / 2 - gap / 2, cy, left_label, size=10, weight="bold")
    txt(ax, cx + w / 2 + gap / 2, cy, right_label, size=10, weight="bold")
    # site arrowheads suggesting the recombination sites
    arrow(ax, (cx - w - gap / 2 - 0.55, cy + 0.02), (cx - w - gap / 2 - 0.12, cy + 0.02),
          DEEP, lw=1.2, mutation=9)
    arrow(ax, (cx + w + gap / 2 + 0.55, cy + 0.02), (cx + w + gap / 2 + 0.12, cy + 0.02),
          TEAL, lw=1.2, mutation=9)
    txt(ax, cx, cy - h / 2 - 0.34, caption, size=9.5, color=GRAY)


def pool_gauge(ax, x, y, w, h, frac, fill_color):
    rbox(ax, x, y, w, h, "white", GRAY, lw=1.3, rad=0.08)
    fw = max(0.05, w * frac)
    rbox(ax, x + 0.03, y + 0.05, fw - 0.06, h - 0.10, fill_color, fill_color,
         lw=0.5, rad=0.06)


fig = plt.figure(figsize=(14.6, 8.8), dpi=200)
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 16)
ax.set_ylim(0, 9.6)
ax.axis("off")

# title + legend
txt(ax, 8, 9.35, "Single-bit switch: state-dependent memory",
    size=15.5, weight="bold")
lx = 0.55
for label, color in [("Integrase (Int)", TEAL), ("RDF", ORANGE),
                     ("BM3R1", DEEP), ("Int\u2013RDF complex", GRAY)]:
    ax.plot([lx, lx + 0.42], [9.00, 9.00], color=color, lw=3.0)
    txt(ax, lx + 0.50, 9.00, label, size=9.5, ha="left")
    lx += 0.62 + 0.115 * len(label)

# ---------------- PB card ----------------
rbox(ax, 0.35, 1.55, 7.10, 6.40, "white", DEEP, lw=2.2, rad=0.22)
txt(ax, 3.90, 7.52, "PB state  \u2014  bit = 0", size=13.5, weight="bold", color=DEEP)
dna_icon(ax, 3.90, 6.55, "attP", "attB", "DNA: attP \u00d7 attB")

rbox(ax, 1.00, 4.62, 2.30, 0.92, LIGHT_BLUE, DEEP, lw=1.5)
txt(ax, 2.15, 5.08, "BM3R1", size=10.5, weight="bold")
arrow(ax, (3.20, 6.18), (2.15, 5.56), DEEP, lw=1.5, rad=0.15)
txt(ax, 2.10, 6.22, "expressed", size=9, color=DEEP, ha="right")

rbox(ax, 4.45, 4.62, 2.30, 0.92, LIGHT_ORANGE, ORANGE, lw=1.5)
txt(ax, 5.60, 5.08, "RDF gene", size=10.5, weight="bold")
txt(ax, 5.60, 4.42, "pBM3R1 repressed \u2192 RDF OFF", size=9, color=ORANGE)
tbar(ax, (3.30, 5.08), (4.45, 5.08), ORANGE, lw=2.0)
txt(ax, 3.90, 5.30, "repressed", size=9, color=ORANGE)

txt(ax, 3.90, 3.88, "RDF pool (memory)", size=10.5, weight="bold")
pool_gauge(ax, 1.20, 3.22, 5.40, 0.55, 0.03, ORANGE)
txt(ax, 3.90, 2.92, "\u2248 0.002 \u00b5M at pulse arrival", size=9.5, color=CHARCOAL)
txt(ax, 3.90, 2.30, "Pool emptied between pulses\n(BM3R1 keeps RDF production off)",
    size=9.5, color=GRAY)

# ---------------- LR card ----------------
rbox(ax, 8.55, 1.55, 7.10, 6.40, "white", TEAL, lw=2.2, rad=0.22)
txt(ax, 12.10, 7.52, "LR state  \u2014  bit = 1", size=13.5, weight="bold", color=TEAL)
dna_icon(ax, 12.10, 6.55, "attL", "attR", "DNA: attL \u00d7 attR")

rbox(ax, 9.20, 4.62, 2.30, 0.92, LIGHT_BLUE, GRAY, lw=1.4)
txt(ax, 10.35, 5.08, "BM3R1", size=10.5, weight="bold", color=GRAY)
txt(ax, 10.35, 4.40, "decays after the flip", size=9, color=GRAY)

rbox(ax, 12.65, 4.62, 2.30, 0.92, LIGHT_ORANGE, ORANGE, lw=1.5)
txt(ax, 13.80, 5.08, "RDF gene", size=10.5, weight="bold")
txt(ax, 13.80, 4.42, "de-repressed \u2192 RDF ON", size=9, color=ORANGE)
arrow(ax, (12.45, 6.18), (13.80, 5.56), ORANGE, lw=1.5, rad=-0.15)
txt(ax, 13.45, 6.08, "produced", size=9, color=ORANGE, ha="left")

txt(ax, 12.10, 3.88, "RDF pool (memory)", size=10.5, weight="bold")
pool_gauge(ax, 9.40, 3.22, 5.40, 0.55, 0.92, ORANGE)
txt(ax, 12.10, 2.92, "\u2248 4.5 \u00b5M at pulse arrival", size=9.5, color=CHARCOAL)
txt(ax, 12.10, 2.30, "Pool fills during the LR dwell\n(no BM3R1 left to repress it)",
    size=9.5, color=GRAY)

# ---------------- transitions ----------------
arrow(ax, (4.60, 7.95), (11.40, 7.95), TEAL, lw=2.8, rad=0.16, mutation=20)
pulse_glyph(ax, 5.00, 8.68)
txt(ax, 5.65, 8.74,
    "clock pulse 1: free Int drives forward recombination (PB \u2192 LR)",
    size=10.5, color=TEAL, weight="bold", ha="left")

arrow(ax, (11.40, 1.55), (4.60, 1.55), DEEP, lw=2.8, rad=0.16, mutation=20)
pulse_glyph(ax, 5.00, 0.74, color=DEEP)
txt(ax, 5.65, 0.80,
    "clock pulse 2: Int\u2013RDF complex drives reverse recombination (LR \u2192 PB)",
    size=10.5, color=DEEP, weight="bold", ha="left")

txt(ax, 8.00, 0.46,
    "Memory = the RDF pool: \u2248 4.5 \u00b5M (LR) vs \u2248 0.002 \u00b5M (PB) at pulse arrival, "
    "\u2248 2500\u00d7 contrast.",
    size=9.5, color=CHARCOAL)
txt(ax, 8.00, 0.20,
    "Constraint: RDF delay > pulse width \u2014 otherwise RDF unlocks during the pulse "
    "and the DNA flips back within the same cycle.",
    size=9.5, color=CHARCOAL)

fig.savefig(OUT, dpi=200, facecolor="white", bbox_inches="tight")
print("saved:", os.path.abspath(OUT))
