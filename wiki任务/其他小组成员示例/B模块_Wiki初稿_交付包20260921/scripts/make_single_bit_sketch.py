# -*- coding: utf-8 -*-
"""TEMPO B module - sketch for the single-bit counter mechanism figure.

Panel A: molecular mechanism of the Int / RDF / BM3R1 switch.
Panel B: one-pulse-one-flip timing logic and the RDF delay constraint.

This is a *sketch* for the art team, not the final wiki figure.
All labels are English (group convention: figures ship in English).

Run:  python make_single_bit_sketch.py
Out:  ../figures/fig03_1_single_bit_mechanism_sketch.png
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

# TEMPO palette (kept identical to the group visual spec)
DEEP = "#304B53"
ORANGE = "#D29144"
TEAL = "#4F9194"
LIGHT_BLUE = "#DCE5E7"
LIGHT_ORANGE = "#F2E2CF"
LIGHT_TEAL = "#DCEBEC"
GRAY = "#8CA0A5"
CHARCOAL = "#33383A"
PURPLE = "#7E6BB5"

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "figures", "fig03_1_single_bit_mechanism_sketch.png")


def box(ax, x, y, w, h, text, fc, ec, fontsize=10.5):
    p = FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.10,rounding_size=0.18",
        linewidth=1.6, edgecolor=ec, facecolor=fc, zorder=2,
    )
    ax.add_patch(p)
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
            fontsize=fontsize, color=CHARCOAL, weight="bold", zorder=3)


def arrow(ax, p1, p2, color, lw=1.8, rad=0.0, ls="-", style="-|>",
          mutation=14, zorder=1):
    a = FancyArrowPatch(
        p1, p2, arrowstyle=style, mutation_scale=mutation,
        color=color, lw=lw, linestyle=ls,
        connectionstyle=f"arc3,rad={rad}", zorder=zorder,
    )
    ax.add_patch(a)


def repression(ax, p1, p2, color, lw=2.0):
    """Repression line: solid line with a T-bar end at p2."""
    arrow(ax, p1, p2, color, lw=lw, style="-", zorder=1)
    dx, dy = p2[0] - p1[0], p2[1] - p1[1]
    n = np.hypot(dx, dy)
    px, py = -dy / n, dx / n
    half = 0.18
    ax.plot([p2[0] - px * half, p2[0] + px * half],
            [p2[1] - py * half, p2[1] + py * half],
            color=color, lw=lw, zorder=2, solid_capstyle="butt")


def label(ax, x, y, text, fontsize=9, color=CHARCOAL, ha="center", va="center"):
    ax.text(x, y, text, fontsize=fontsize, color=color, ha=ha, va=va, zorder=4)


fig = plt.figure(figsize=(14.0, 7.0), dpi=200)
axA = fig.add_axes([0.015, 0.04, 0.52, 0.85])
axB = fig.add_axes([0.585, 0.14, 0.40, 0.72])

# ----------------------------------------------------------------------------
# Panel A - mechanism
# ----------------------------------------------------------------------------
axA.set_xlim(0, 11.5)
axA.set_ylim(0, 10)
axA.axis("off")

box(axA, 0.25, 7.00, 2.70, 1.35, "$\\phi$C31 Int pulse\n(from the clock)", LIGHT_BLUE, DEEP)
box(axA, 3.65, 7.00, 2.70, 1.35, "Int\u2013RDF\ncomplex", LIGHT_TEAL, GRAY)
box(axA, 7.05, 7.00, 2.90, 1.35, "RDF pool\n(memory)", LIGHT_ORANGE, ORANGE)
box(axA, 2.60, 1.15, 3.40, 1.55, "DNA state\nPB (1 \u2212 S)  \u21c4  LR (S)", LIGHT_BLUE, DEEP)
box(axA, 7.60, 1.10, 2.70, 1.20, "BM3R1 (T)", LIGHT_BLUE, DEEP)

# binding equilibria
arrow(axA, (2.95, 7.68), (3.65, 7.68), CHARCOAL, style="<|-|>", mutation=11)
arrow(axA, (6.35, 7.68), (7.05, 7.68), CHARCOAL, style="<|-|>", mutation=11)
label(axA, 3.30, 8.72, "reversible\nbinding", fontsize=8.4)
label(axA, 6.70, 8.72, "reversible\nbinding", fontsize=8.4)

# forward / reverse recombination
arrow(axA, (2.20, 7.00), (3.90, 2.70), TEAL, lw=2.2, rad=-0.12)
label(axA, 2.00, 4.62, "Forward recombination\n(PB \u2192 LR)\nfree Int drives",
      fontsize=8.9, color=TEAL, ha="right")
arrow(axA, (5.00, 7.00), (5.00, 2.70), DEEP, lw=2.2)
label(axA, 5.15, 4.62, "Reverse recombination\n(LR \u2192 PB)\ncomplex drives",
      fontsize=8.9, color=DEEP, ha="left")

# state-dependent production split
arrow(axA, (6.00, 1.90), (7.60, 1.75), DEEP, lw=1.6)
label(axA, 6.80, 2.48, "PB drives BM3R1", fontsize=8.5, color=DEEP)
arrow(axA, (5.90, 2.70), (7.30, 7.00), ORANGE, lw=1.6, rad=-0.12)
label(axA, 6.85, 5.58, "LR state:\nRDF production\n(de-repressed)",
      fontsize=8.5, color=ORANGE, ha="left")

# BM3R1 represses RDF production
repression(axA, (8.95, 2.30), (8.95, 7.00), ORANGE)
label(axA, 9.15, 4.62, "BM3R1 represses\nRDF production",
      fontsize=8.8, color=ORANGE, ha="left")

# footnote - the three points that must not be drawn wrong
label(axA, 0.15, 0.52,
      "1. Memory is stored in the RDF pool, not in Int.\n"
      "2. Reverse recombination is driven by the Int\u2013RDF complex.\n"
      "3. PB drives BM3R1; low BM3R1 in LR lets RDF accumulate.",
      fontsize=8.8, ha="left", color=CHARCOAL)

axA.text(5.75, 9.72, "A  Single-bit switch: Int / RDF / BM3R1",
         fontsize=12.5, weight="bold", color=CHARCOAL, ha="center")

# ----------------------------------------------------------------------------
# Panel B - one-pulse-one-flip timing
# ----------------------------------------------------------------------------
t = np.linspace(0, 11.6, 1200)

pulse = (np.exp(-((t - 1.0) / 0.42) ** 2)
         + np.exp(-((t - 11.0) / 0.42) ** 2))
s_lr = (0.5 * (1 + np.tanh((t - 1.75) / 0.30))
        * 0.5 * (1 + np.tanh((11.45 - t) / 0.30)))
bm3 = (1.0 / (1.0 + np.exp((t - 2.05) / 0.75))
       * 1.0 / (1.0 + np.exp(-(t - 11.35) / 0.50)))
rdf = (0.88 / (1.0 + np.exp(-(t - 4.05) / 0.85))
       * (1 - 0.55 * np.exp(-((t - 11.05) / 0.45) ** 2)))

axB.plot(t, pulse, color=TEAL, lw=2.4, label="Int input pulse (from clock)")
axB.plot(t, s_lr, color=DEEP, lw=2.4, label="LR fraction $S$ (bit state)")
axB.plot(t, bm3, color=ORANGE, lw=2.2, label="BM3R1 ($T$, delay)")
axB.plot(t, rdf, color=PURPLE, lw=2.2, label="RDF pool ($R$, memory)")

axB.axvline(1.0, color=GRAY, ls="--", lw=1.1)
axB.axvline(11.0, color=GRAY, ls="--", lw=1.1)
axB.axvspan(3.0, 10.6, color=LIGHT_ORANGE, alpha=0.45, zorder=0)
axB.text(6.2, 0.09, "RDF accumulation window (memory)", fontsize=8.6,
         color="#8a5a20", ha="center")

axB.annotate("pulse 1: RDF pool empty\n\u2192 PB \u2192 LR", xy=(1.02, 0.72),
             xytext=(2.15, 0.17), fontsize=8.8, color=CHARCOAL, ha="center",
             arrowprops=dict(arrowstyle="->", color=GRAY, lw=1.0))
axB.annotate("RDF rise starts only after BM3R1 decays:\ndelay > pulse width",
             xy=(3.45, 0.22), xytext=(5.35, 0.38), fontsize=8.8, color=CHARCOAL,
             ha="center",
             arrowprops=dict(arrowstyle="->", color=GRAY, lw=1.0))
axB.annotate("pulse 2: complex forms\n\u2192 LR \u2192 PB", xy=(11.02, 0.52),
             xytext=(9.35, 0.16), fontsize=8.8, color=CHARCOAL, ha="center",
             arrowprops=dict(arrowstyle="->", color=GRAY, lw=1.0))

axB.set_xlim(0, 11.6)
axB.set_ylim(-0.05, 1.42)
axB.set_xlabel("time (schematic, one clock cycle)", fontsize=9.5)
axB.set_ylabel("normalized level (schematic)", fontsize=9.5)
axB.set_xticks([])
axB.set_yticks([])
for spine in ("top", "right"):
    axB.spines[spine].set_visible(False)
axB.legend(fontsize=8.2, loc="upper right", frameon=False)
axB.set_title("B  One-pulse-one-flip timing logic", fontsize=12.5,
              weight="bold", color=CHARCOAL, pad=8)

fig.suptitle("Single-bit counter mechanism \u2014 sketch for the TEMPO wiki "
             "(B module, 2026-09-21)", fontsize=14, weight="bold",
             color=CHARCOAL, y=0.99)
fig.savefig(OUT, dpi=200, facecolor="white", bbox_inches="tight")
print("saved:", os.path.abspath(OUT))
