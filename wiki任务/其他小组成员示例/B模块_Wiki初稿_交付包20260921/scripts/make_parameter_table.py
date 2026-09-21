# -*- coding: utf-8 -*-
"""TEMPO B module - parameter table figure for the Chapter 3 wiki page.

Follows the group rule "parameter tables should be delivered as images".
All labels English (group convention).

Run:  python make_parameter_table.py
Out:  ../figures/fig03_7_parameter_table.png
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DEEP = "#304B53"
ORANGE = "#D29144"
TEAL = "#4F9194"
GRAY = "#8CA0A5"
CHARCOAL = "#33383A"
LIGHT = "#F3F5F6"
LIGHT_ORANGE = "#F2E2CF"

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "figures", "fig03_7_parameter_table.png")

rows = [
    # (Quantity, Value / range, Nature, Evidence)
    ("Model conditions", "", "", ""),
    ("Transcription coefficients $\\alpha_B=\\alpha_R$", "120 h⁻¹", "model condition", "literature-informed"),
    ("mRNA effective loss $\\gamma_B=\\gamma_R$", "4 h⁻¹", "model condition", "literature-informed"),
    ("Growth dilution $\\mu$ (50 min doubling)", "0.8318 h⁻¹", "model condition", "team condition"),
    ("Total DNA $D_{tot}$", "0.017 µM", "model condition", "literature value"),
    ("Cell volume / conversion", "1 fL; 1 µM ≈ 602 molecules", "assumption", "pending calibration"),
    ("Primary scheme (non-zero RDF promoter leak)", "", "", ""),
    ("BM3R1 translation $\\beta_B$", "2 h⁻¹", "model-effective parameter", "pending calibration"),
    ("RDF translation $\\beta_R$", "8 h⁻¹", "model-effective parameter", "pending calibration"),
    ("Int input flux scale $s_I$", "0.30", "design variable", "paired with $k_I$"),
    ("Int extra clearance $k_I$ (tag)", "12 h⁻¹", "model-effective parameter", "pending calibration"),
    ("RDF promoter leak floor $\\ell_R$", "0.008 (0.8%)", "literature prior (Cello B1)", "pending calibration"),
    ("RDF return fraction $\\eta$", "1 (working hypothesis)", "mechanistic assumption", "to be tested"),
    ("Design targets (measurable form)", "", "", ""),
    ("BM3R1 max synthesis flux / PB steady", "1.02 µM/h; 1.226 µM ≈ 740 copies/cell", "model target", "pending calibration"),
    ("RDF max synthesis flux / de-repressed steady", "4.08 µM/h; 4.905 µM ≈ 2950 copies/cell", "model target", "pending calibration"),
    ("BM3R1 : RDF maximum capacity", "≈ 1 : 4", "model target", "mechanism-inferred"),
    ("Operating windows", "", "", ""),
    ("Int degradation tag $k_{tag,int}$", "deterministic [3, 24] h⁻¹; stochastic-stable 6–18 h⁻¹", "operating window", "mechanism-inferred"),
    ("BM3R1 half-repression constant $K_B$ / Hill $n$", "18.6 nM (12–50 nM) / 3.4 (2.9–3.4)", "operating window", "mechanism-inferred"),
    ("RDF relative leak $\\ell_R$", "≤0.8% design target; ≤1% tolerance", "operating window", "mechanism-inferred"),
    ("Start-up condition", "known initial state + trough entry into first pulse", "operating condition", "mechanism-inferred"),
    ("ϕC31 RBS scale (input interface)", "0.45 (zero-leak point); 0.30 (non-zero-leak scheme)", "design variable", "pending calibration"),
]

fig = plt.figure(figsize=(13.6, 9.2), dpi=200)
ax = fig.add_axes([0.02, 0.02, 0.96, 0.92])
ax.axis("off")

col_widths = [0.42, 0.34, 0.14, 0.10]
x0 = 0.0
xs = [x0]
for w in col_widths[:-1]:
    xs.append(xs[-1] + w)

headers = ["Quantity", "Value / range", "Nature", "Evidence"]
n = len(rows)
y_top = 0.96
row_h = (y_top - 0.02) / (n + 1)

# header
for x, w, htxt in zip(xs, col_widths, headers):
    ax.add_patch(plt.Rectangle((x, y_top - row_h), w, row_h,
                               facecolor=DEEP, edgecolor="white", lw=1.2,
                               transform=ax.transAxes, zorder=2))
    ax.text(x + 0.012, y_top - row_h / 2, htxt, transform=ax.transAxes,
            fontsize=12, color="white", weight="bold", va="center", zorder=3)

# rows
for i, (q, v, nat, ev) in enumerate(rows):
    y = y_top - row_h * (i + 2)
    section = (v == "" and nat == "")
    bg = LIGHT_ORANGE if section else (LIGHT if i % 2 == 0 else "white")
    for x, w in zip(xs, col_widths):
        ax.add_patch(plt.Rectangle((x, y), w, row_h, facecolor=bg,
                                   edgecolor="#DDE2E4", lw=0.8,
                                   transform=ax.transAxes, zorder=1))
    if section:
        ax.text(xs[0] + 0.012, y + row_h / 2, q, transform=ax.transAxes,
                fontsize=12, color=DEEP, weight="bold", va="center", zorder=3)
        continue
    ax.text(xs[0] + 0.012, y + row_h / 2, q, transform=ax.transAxes,
            fontsize=10.5, color=CHARCOAL, va="center", zorder=3)
    ax.text(xs[1] + 0.012, y + row_h / 2, v, transform=ax.transAxes,
            fontsize=10.5, color=CHARCOAL, va="center", zorder=3)
    ax.text(xs[2] + 0.012, y + row_h / 2, nat, transform=ax.transAxes,
            fontsize=9.5, color=TEAL, va="center", zorder=3)
    ev_color = GRAY if "pending" in ev else (DEEP if "literature" in ev else ORANGE)
    ax.text(xs[3] + 0.012, y + row_h / 2, ev, transform=ax.transAxes,
            fontsize=9.5, color=ev_color, va="center", zorder=3)

ax.text(0.0, 0.995,
        "Chapter 3 - single-bit counter: parameters, targets and operating windows",
        transform=ax.transAxes, fontsize=13.5, weight="bold", color=CHARCOAL,
        va="bottom")
ax.text(0.0, 0.0,
        "Model conditions: 50 min doubling, unloaded clock input (period 10.59 h). "
        "Values are model quantities, not experimental specifications; "
        "all ranges carry the definitions given in the Chapter 3 text. "
        "Evidence labels: literature value / model-effective parameter / design variable / pending calibration.",
        transform=ax.transAxes, fontsize=8.6, color=GRAY, va="top")

fig.savefig(OUT, dpi=200, facecolor="white", bbox_inches="tight")
print("saved:", os.path.abspath(OUT))
