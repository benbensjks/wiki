# -*- coding: utf-8 -*-
"""TEMPO B module - design-window concept figure (Chapter 3).

Panel A: why the Int-degradation-tag window has edges
         (three regimes with schematic S traces: multi-flip / pass / stuck).
Panel B: the joint expression window in the (beta_R, beta_B) plane (schematic),
         with the failure mode labelled on each side and the leak shift arrow.

Concept figure only; the measured maps are the data figures.
All labels English. TEMPO palette.

Run:  python make_design_window_v1.py
Out:  ../fig03_design_window_v1.png
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
GREEN = "#8FBF9F"
GREEN_DARK = "#5E9E77"
RED = "#C97A5B"

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "fig03_design_window_v1.png")


def txt(ax, x, y, s, size=10, color=CHARCOAL, ha="center", va="center",
        weight="normal", style="normal"):
    ax.text(x, y, s, fontsize=size, color=color, ha=ha, va=va,
            weight=weight, style=style, zorder=5)


def pulse_glyph(ax, x, y, w, h, color=TEAL, lw=1.6):
    xs = [x, x + 0.22 * w, x + 0.22 * w, x + 0.70 * w, x + 0.70 * w, x + w]
    ys = [y, y, y + h, y + h, y, y]
    ax.plot(xs, ys, color=color, lw=lw, zorder=4)


fig, (axA, axB) = plt.subplots(1, 2, figsize=(14.6, 6.4), dpi=200,
                               gridspec_kw={"width_ratios": [1.15, 1.0], "wspace": 0.10})

# ============================================================ Panel A
axA.set_xlim(0, 10)
axA.set_ylim(0, 10)
axA.axis("off")
txt(axA, 5.0, 9.72, "A  Why the tag window has edges", size=13, weight="bold")

# window bar (0-30 h-1 mapped to x 0.5-9.5)
x0, x1 = 0.5, 9.5


def kx(k):
    return x0 + (k / 30.0) * (x1 - x0)


bar_y, bar_h = 7.55, 0.55
axA.add_patch(Rectangle((x0, bar_y), x1 - x0, bar_h, facecolor="white",
                        edgecolor=GRAY, lw=1.4, zorder=1))
axA.add_patch(Rectangle((kx(0), bar_y), kx(3) - kx(0), bar_h, facecolor=LIGHT_ORANGE,
                        edgecolor="none", zorder=2))
axA.add_patch(Rectangle((kx(24), bar_y), kx(30) - kx(24), bar_h, facecolor=LIGHT_ORANGE,
                        edgecolor="none", zorder=2))
axA.add_patch(Rectangle((kx(3), bar_y), kx(24) - kx(3), bar_h, facecolor=GREEN,
                        edgecolor="none", zorder=2))
axA.add_patch(Rectangle((kx(6), bar_y), kx(18) - kx(6), bar_h, facecolor=GREEN_DARK,
                        edgecolor="none", zorder=2))
axA.add_patch(Rectangle((x0, bar_y), x1 - x0, bar_h, facecolor="none",
                        edgecolor=GRAY, lw=1.4, zorder=3))

for k, lab in [(0, "0"), (3, "3"), (6, "6"), (18, "18"), (24, "24"), (30, "30")]:
    axA.plot([kx(k), kx(k)], [bar_y - 0.12, bar_y], color=GRAY, lw=1.0, zorder=3)
    txt(axA, kx(k), bar_y - 0.36, lab, size=8.5, color=GRAY)

txt(axA, kx(1.5), bar_y + 0.95, "too weak", size=9.5, color=ORANGE, weight="bold")
txt(axA, kx(27), bar_y + 0.95, "too strong", size=9.5, color=ORANGE, weight="bold")
txt(axA, kx(13.5), bar_y + 0.98, "deterministic pass  [3, 24] h$^{-1}$",
    size=9.5, color=DEEP, weight="bold")
txt(axA, kx(12), bar_y + bar_h / 2, "stochastic [6, 18]", size=8.5, color="white",
    weight="bold")
txt(axA, kx(13.5), bar_y - 0.85,
    "Int degradation-tag rate $k_{tag,int}$", size=10, color=CHARCOAL)

# mini S traces -------------------------------------------------------------
def s_clean(t):
    out = np.zeros_like(t)
    for i, ti in enumerate(t):
        ph = ti % 10.0
        out[i] = 1.0 if 2.5 < ph < 8.0 else 0.0
    return out


def s_multiflip(t):
    out = np.zeros_like(t)
    for i, ti in enumerate(t):
        ph = ti % 10.0
        if ph < 2.0:
            out[i] = 0.0
        elif ph < 3.2:
            out[i] = 1.0
        elif ph < 4.4:
            out[i] = 0.0
        elif ph < 5.6:
            out[i] = 1.0
        elif ph < 7.0:
            out[i] = 0.35
        else:
            out[i] = 0.75
    return out


def s_stuck(t):
    out = np.zeros_like(t)
    for i, ti in enumerate(t):
        ph = ti % 10.0
        out[i] = 1.0 if ph > 2.5 else 0.0
    return out


def mini(ax, xc, trace, caption, cap_color):
    w, h = 2.5, 1.05
    x_left, y_base = xc - w / 2, 3.35
    t = np.linspace(0, 20, 800)
    s = trace(t)
    ax.plot(x_left + (t / 20.0) * w, y_base + s * h, color=DEEP, lw=1.8, zorder=4)
    # pulse glyphs above
    for k in range(2):
        pulse_glyph(ax, x_left + k * w / 2 + 0.18, y_base + h + 0.45, w / 2 - 0.36, 0.30)
    ax.plot([x_left, x_left + w], [y_base, y_base], color=GRAY, lw=0.9, zorder=2)
    txt(ax, xc, y_base - 0.42, caption, size=8.8, color=cap_color, weight="bold")


mini(axA, 1.95, s_multiflip, "too weak:\ninter-pulse Int \u2192 extra flips", ORANGE)
mini(axA, 5.00, s_clean, "pass:\none pulse = one flip", GREEN_DARK)
mini(axA, 8.05, s_stuck, "too strong:\npeak clipped \u2192 stuck in LR", ORANGE)

txt(axA, 5.0, 1.55,
    "Schematic S(t) traces (two pulse periods each; flat line = S = 0).",
    size=8.6, color=GRAY)

# ============================================================ Panel B
axB.set_xlim(0, 10)
axB.set_ylim(0, 10)
axB.axis("off")
txt(axB, 5.0, 9.72, "B  The joint expression window (schematic)", size=13, weight="bold")

# axes arrows
axB.add_patch(FancyArrowPatch((1.0, 1.2), (9.4, 1.2), arrowstyle="-|>",
                              mutation_scale=16, color=GRAY, lw=1.6))
axB.add_patch(FancyArrowPatch((1.0, 1.2), (1.0, 8.9), arrowstyle="-|>",
                              mutation_scale=16, color=GRAY, lw=1.6))
txt(axB, 9.4, 0.82, "RDF expression strength  ($\\beta_R$)", size=10, ha="right")
txt(axB, 0.72, 8.9, "BM3R1 expression\nstrength  ($\\beta_B$)", size=10, ha="right")

# pass plateau
plateau = FancyBboxPatch((3.4, 3.4), 3.4, 3.2,
                         boxstyle="round,pad=0.10,rounding_size=0.35",
                         facecolor=GREEN, edgecolor=GREEN_DARK, lw=2.0, zorder=3)
axB.add_patch(plateau)
txt(axB, 5.1, 5.35, "pass plateau", size=11, weight="bold", color="white")
txt(axB, 5.1, 4.72, "$\\beta_R$ 4\u20138  \u00d7  $\\beta_B$ 0.75\u20136",
    size=9, color="white")
txt(axB, 5.1, 4.18, "(measured)", size=8.2, color="white")

# untested range above
axB.add_patch(Rectangle((3.6, 6.95), 3.0, 1.05, facecolor="white",
                        edgecolor=GRAY, lw=1.2, linestyle="--", zorder=2))
txt(axB, 5.1, 7.47, "beyond tested range", size=8.6, color=GRAY)

# failure labels
txt(axB, 2.15, 5.3, "RDF too weak:\nreverse flip\nincomplete", size=9, color=ORANGE,
    weight="bold")
txt(axB, 8.05, 5.3, "RDF too strong:\nleak floor amplified\n\u2192 spontaneous flips",
    size=9, color=ORANGE, weight="bold")
txt(axB, 5.1, 2.55, "BM3R1 too weak: pool not emptied \u2192 multi-flip",
    size=9, color=ORANGE, weight="bold")

# leak shift arrow
axB.add_patch(FancyArrowPatch((8.8, 7.05), (6.95, 6.35), arrowstyle="-|>",
                              mutation_scale=15, color=DEEP, lw=2.0,
                              linestyle="--", connectionstyle="arc3,rad=0.20"))
txt(axB, 9.5, 8.05,
    "higher promoter leak\nshrinks the window\ntoward lower RDF strength",
    size=8.8, color=DEEP, ha="right")

txt(axB, 5.0, 0.55,
    "Schematic; the measured pass/fail map is the feasible-region data figure.",
    size=8.6, color=GRAY)

fig.savefig(OUT, dpi=200, facecolor="white", bbox_inches="tight")
print("saved:", os.path.abspath(OUT))
