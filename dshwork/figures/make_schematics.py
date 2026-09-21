"""Generate the code-drawable schematic figures of the two figure briefs.

Covers, from `04_第3章_作图说明.md` and `05_第4章_作图说明.md` plus the
mechanism-figure list in `01_9.20交付_...md`:

   fig01_modeling_workflow      Design -> Build -> Test -> Learn -> Decision loop
   fig03_single_bit_mechanism   Int / RDF / BM3R1 switch (chapter 3, Figure 3-1)
   fig04_ffl_carry              I1-FFL carry module (chapter 4, Figure 4-1)
   fig04_6_bm3r1_shutdown_interface   BM3R1 shared interface (chapter 4, Figure 4-6)
   fig07_full_system_coupling   whole-system coupling (chapter 1 / overview)
   fig08_design_panel           modeling-guided design panel

Data figures (readout strips, heatmaps, dose-response) are NOT produced here:
they need the simulation outputs and are handled separately.

Usage
-----
    python make_schematics.py --out ./out
    python make_schematics.py --out ./out --only fig04_ffl_carry
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from tempo_style import (BLUE_GRAY, DEEP_BLUE, INK, LIGHT_BLUE, LIGHT_GRAY,
                         LIGHT_ORANGE, LIGHT_TEAL, ORANGE, RISK_ORANGE, TEAL,
                         and_node, arrow, binding, box, hline, hollow_arrow,
                         legend_lines, new_canvas, note, repression, save,
                         scope)


# --------------------------------------------------------------------------
# 1. Modeling workflow
# --------------------------------------------------------------------------
def fig01_modeling_workflow():
    fig, ax = new_canvas(14.0, 6.6)
    stages = [
        ('Design', 'engineering question\n\u2192 computable criterion', LIGHT_BLUE, DEEP_BLUE),
        ('Build', 'explicit continuous\nODE model', LIGHT_GRAY, DEEP_BLUE),
        ('Test', 'criteria, scans,\nperturbations', LIGHT_TEAL, TEAL),
        ('Learn', 'mechanism, limits,\nuncalibrated quantities', LIGHT_ORANGE, ORANGE),
        ('Decision', 'design rules,\nor redesign', LIGHT_GRAY, RISK_ORANGE),
    ]
    xs = [1.7, 4.4, 7.1, 9.8, 12.5]
    y = 4.3
    for (title, caption, fc, ec), cx in zip(stages, xs):
        box(ax, cx, y, 2.3, 1.05, title, fc=fc, ec=ec, fontsize=13, bold=True)
        note(ax, cx, y - 0.78, caption, color=INK, fontsize=8.8,
             ha='center', va='top', style='normal')
    for a, b in zip(xs[:-1], xs[1:]):
        arrow(ax, (a + 1.15, y), (b - 1.15, y), color=DEEP_BLUE, lw=2.0)

    # feedback loop back to Design
    yb = 2.05
    ax.plot([xs[-1], xs[-1], xs[0], xs[0]], [y - 0.53, yb, yb, y - 0.53],
            color=DEEP_BLUE, lw=2.0, zorder=2, solid_capstyle='round')
    arrow(ax, (xs[0], yb + 0.28), (xs[0], y - 0.55), color=DEEP_BLUE, lw=2.0)
    note(ax, (xs[0] + xs[-1]) / 2, yb + 0.12, 'Learned design rule feeds back',
         color=DEEP_BLUE, fontsize=9.5, ha='center', va='bottom')

    note(ax, 7.0, 6.15, 'TEMPO modeling workflow', color=DEEP_BLUE,
         fontsize=15, ha='center', va='center', style='normal')
    note(ax, 0.3, 1.25,
         '"Test" here means in-model validation, not wet-lab experiments.',
         color=BLUE_GRAY, fontsize=9)
    return save(fig, OUT / 'fig01_modeling_workflow.png')


# --------------------------------------------------------------------------
# 2. Single-bit switch mechanism  (chapter 3, Figure 3-1)
# --------------------------------------------------------------------------
def fig03_single_bit_mechanism():
    fig, ax = new_canvas(14.0, 9.6)

    box(ax, 2.5, 7.2, 2.5, 0.95, 'Integrase (Int)', fc=LIGHT_TEAL, ec=TEAL,
        fontsize=11, bold=True)
    box(ax, 7.0, 7.2, 2.9, 1.05, 'Int\u2013RDF complex', fc=LIGHT_TEAL,
        ec=INK, fontsize=11, bold=True)
    box(ax, 11.5, 7.2, 2.9, 1.05, 'RDF pool\n(memory)', fc=LIGHT_ORANGE,
        ec=ORANGE, fontsize=11, bold=True)
    box(ax, 11.5, 5.15, 2.9, 0.9, 'BM3R1', fc=LIGHT_BLUE, ec=DEEP_BLUE,
        fontsize=11, bold=True)
    box(ax, 7.0, 2.75, 5.8, 1.35,
        'DNA state\nPB  (1 \u2212 S)          LR  (S)', fc=LIGHT_GRAY,
        ec=DEEP_BLUE, fontsize=11.5, bold=True)

    binding(ax, (3.8, 7.2), (5.5, 7.2), color=INK)
    binding(ax, (8.5, 7.2), (10.0, 7.2), color=INK)
    note(ax, 4.65, 7.42, 'reversible\nbinding', color=INK, fontsize=8,
         ha='center', va='bottom', style='normal')

    arrow(ax, (2.5, 6.72), (5.35, 3.45), color=TEAL, lw=2.2, rad=-0.08)
    note(ax, 3.05, 4.95, 'Forward recombination\n(PB \u2192 LR)', color=TEAL,
         fontsize=9.5, ha='left', va='center', style='normal')

    arrow(ax, (7.0, 6.67), (7.0, 3.45), color=DEEP_BLUE, lw=2.2)
    note(ax, 7.28, 5.05, 'Reverse recombination\n(LR \u2192 PB)', color=DEEP_BLUE,
         fontsize=9.5, ha='left', va='center', style='normal')

    arrow(ax, (9.85, 3.25), (11.5, 6.65), color=ORANGE, lw=2.0, rad=-0.18)
    note(ax, 11.85, 4.85, 'S drives\nRDF production', color=ORANGE, fontsize=9.5,
         ha='left', va='center', style='normal')

    arrow(ax, (8.8, 2.10), (10.6, 4.72), color=DEEP_BLUE, lw=1.8, rad=-0.22)
    note(ax, 9.15, 1.62, 'PB fraction drives BM3R1', color=DEEP_BLUE,
         fontsize=9.5, ha='left', va='center', style='normal')

    repression(ax, (11.5, 5.62), (11.5, 6.65), color=DEEP_BLUE, lw=2.0)
    note(ax, 11.78, 6.14, 'BM3R1 represses\nRDF production', color=DEEP_BLUE,
         fontsize=9.5, ha='left', va='center', style='normal')

    note(ax, 7.0, 9.15, 'Single-bit switch', color=DEEP_BLUE, fontsize=15,
         ha='center', va='center', style='normal')
    legend_lines(ax, 0.35, 8.75,
                 [('repression', 'Repression'), ('activation', 'Activation'),
                  ('binding', 'Reversible binding')], fontsize=8.8, dy=0.45)
    note(ax, 0.35, 1.05,
         '1.  Memory is stored in the RDF pool, not in Int\n'
         '2.  Reverse recombination is driven by the complex\n'
         '3.  S high \u2192 BM3R1 low \u2192 RDF high',
         color=INK, fontsize=9.5, va='top')
    return save(fig, OUT / 'fig03_single_bit_mechanism.png')


# --------------------------------------------------------------------------
# 3. I1-FFL carry module  (chapter 4, Figure 4-1)
# --------------------------------------------------------------------------
def _mini_panel(ax, x0, y0, w, h, title, curves, note_text, fill=None):
    """Small schematic time-course drawn directly in data coordinates."""
    ax.add_patch(plt_rect(x0, y0, w, h))
    note(ax, x0 + w / 2, y0 + h + 0.12, title, color=DEEP_BLUE, fontsize=9.5,
         ha='center', va='bottom', style='normal')
    if fill is not None:
        xs, ys = fill
        ax.fill_between(x0 + xs * w, y0, y0 + ys * h, color=LIGHT_TEAL,
                        zorder=4.5, linewidth=0)
    for xs, ys, color, lw, ls in curves:
        ax.plot(x0 + xs * w, y0 + ys * h, color=color, lw=lw, ls=ls, zorder=5)
    note(ax, x0 + w - 0.06, y0 + 0.06, 'time', color=BLUE_GRAY, fontsize=8.0,
         ha='right', va='bottom', style='normal')
    note(ax, x0 + w / 2, y0 - 0.16, note_text, color=INK, fontsize=8.5,
         ha='center', va='top', style='normal')


def plt_rect(x, y, w, h):
    from matplotlib.patches import Rectangle
    return Rectangle((x, y), w, h, fc='white', ec=BLUE_GRAY, lw=1.0, zorder=4)


def fig04_ffl_carry():
    fig, ax = new_canvas(15.6, 9.4)

    box(ax, 4.0, 8.55, 4.6, 0.95, 'Previous bit DNA state\nPB fraction = 1 \u2212 S',
        fc=LIGHT_GRAY, ec=DEEP_BLUE, fontsize=10.5, bold=True)
    box(ax, 4.0, 7.05, 3.4, 0.9, 'Activator A\n(coherent arm)', fc=LIGHT_ORANGE,
        ec=ORANGE, fontsize=10.5, bold=True)
    box(ax, 4.0, 4.95, 3.4, 0.95, 'Repressor F\n(incoherent, delayed arm)',
        fc=LIGHT_BLUE, ec=BLUE_GRAY, fontsize=10.5, bold=True)
    and_node(ax, 7.7, 6.0, r=0.36)
    box(ax, 11.6, 6.0, 3.4, 0.9, 'Clock gate\nH(Int0; 0.4, 3)', fc=LIGHT_TEAL,
        ec=DEEP_BLUE, fontsize=10.5, bold=True)
    box(ax, 7.7, 4.45, 4.6, 0.9, 'Carry promoter output  g', fc=LIGHT_TEAL,
        ec=DEEP_BLUE, fontsize=11, bold=True)
    box(ax, 7.7, 3.10, 5.2, 0.95, 'Next-bit integrase pulse\n(through M_I / I_u / I)',
        fc=LIGHT_GRAY, ec=TEAL, fontsize=10.5)
    box(ax, 7.7, 1.75, 3.8, 0.9, 'Next bit:  S flip', fc=LIGHT_GRAY,
        ec=DEEP_BLUE, fontsize=10.5, bold=True)

    arrow(ax, (4.0, 8.06), (4.0, 7.52), color=DEEP_BLUE, lw=2.0)
    arrow(ax, (4.0, 6.58), (4.0, 5.45), color=ORANGE, lw=2.0, rad=-0.30)
    note(ax, 4.28, 6.02, 'Expression delay', color=ORANGE, fontsize=9.5,
         ha='left', va='center', style='normal')

    arrow(ax, (5.75, 7.10), (7.45, 6.22), color=ORANGE, lw=2.0, rad=0.12)
    arrow(ax, (5.75, 4.98), (7.45, 5.80), color=BLUE_GRAY, lw=2.0, rad=-0.12)
    arrow(ax, (9.85, 6.0), (8.10, 6.0), color=DEEP_BLUE, lw=2.0)
    note(ax, 6.50, 6.95, 'A and F are\nmultiplied (AND)', color=DEEP_BLUE,
         fontsize=9.0, ha='center', va='center', style='normal')

    arrow(ax, (7.7, 5.60), (7.7, 4.92), color=DEEP_BLUE, lw=2.2)
    arrow(ax, (7.7, 3.98), (7.7, 3.60), color=DEEP_BLUE, lw=2.2)
    arrow(ax, (7.7, 2.60), (7.7, 2.22), color=DEEP_BLUE, lw=2.2)

    # ---- two small time-course insets on the right -----------------------
    # Panel 1: no delayed arm -> g simply tracks A, i.e. a plateau.
    # Panel 2: delayed arm    -> A rises fast, F catches up slowly, g is a pulse.
    t = np.linspace(0, 1, 500)
    g_plateau = np.where(t < 0.25, 0.0, 1.0)
    a_fast = 1 - np.exp(-np.clip(t - 0.15, 0, None) / 0.06)
    f_slow = 1 - np.exp(-np.clip(t - 0.15, 0, None) / 0.22)
    g_pulse = a_fast * (1 - f_slow)

    _mini_panel(ax, 11.0, 7.55, 3.9, 1.25, 'Without the delayed arm',
                [(t, g_plateau, TEAL, 2.4, '-')],
                'g tracks A: a plateau', fill=(t, g_plateau))
    _mini_panel(ax, 11.0, 4.30, 3.9, 1.25, 'With the delayed arm',
                [(t, a_fast, ORANGE, 1.9, '-'), (t, f_slow, BLUE_GRAY, 1.9, '-'),
                 (t, g_pulse, TEAL, 2.4, '-')],
                'F catches up: g is a narrow pulse')

    note(ax, 14.62, 8.58, 'g', color=TEAL, fontsize=9.5, ha='right',
         va='center', style='normal')
    note(ax, 13.30, 5.34, 'A', color=ORANGE, fontsize=9.5, ha='left',
         va='center', style='normal')
    note(ax, 14.60, 5.20, 'F', color=BLUE_GRAY, fontsize=9.5, ha='right',
         va='center', style='normal')
    note(ax, 12.12, 5.04, 'g', color=TEAL, fontsize=9.5, ha='left',
         va='center', style='normal')

    note(ax, 7.8, 9.05, 'I1-FFL carry module', color=DEEP_BLUE, fontsize=15,
         ha='center', va='center', style='normal')
    note(ax, 0.3, 1.15,
         '\u2022  A and F are multiplied (AND), not added\n'
         '\u2022  F lags A \u2014 this is what makes a pulse instead of a plateau\n'
         '\u2022  All three factors are required: A, F, and the clock gate',
         color=INK, fontsize=9.8, va='top')
    return save(fig, OUT / 'fig04_ffl_carry.png')


# --------------------------------------------------------------------------
# 4. BM3R1 shared interface  (chapter 4, Figure 4-6)
# --------------------------------------------------------------------------
def fig04_6_bm3r1_shutdown_interface():
    fig, ax = new_canvas(15.2, 8.4)

    box(ax, 7.0, 6.55, 5.4, 1.05, 'Shared pool:  BM3R1', fc=LIGHT_ORANGE,
        ec=ORANGE, lw=2.0, ls=(0, (6, 3)), fontsize=13, bold=True)
    box(ax, 2.9, 4.05, 4.8, 1.35,
        'Counter side\nBM3R1 represses RDF production', fc=LIGHT_TEAL,
        ec=TEAL, fontsize=10.5)
    box(ax, 11.3, 4.05, 5.2, 1.55,
        'Shutdown side\nBM3R1 represses shutdown sRNA\n(molecular timer)',
        fc=LIGHT_BLUE, ec=DEEP_BLUE, fontsize=10.5)

    repression(ax, (5.35, 6.05), (4.35, 4.80), color=DEEP_BLUE, lw=2.0, rad=-0.15)
    repression(ax, (8.65, 6.05), (9.75, 4.90), color=DEEP_BLUE, lw=2.0, rad=0.15)
    note(ax, 3.05, 5.62, 'represses RDF production', color=DEEP_BLUE,
         fontsize=9.2, ha='left', va='center', style='normal')
    note(ax, 9.05, 5.62, 'represses shutdown sRNA', color=DEEP_BLUE,
         fontsize=9.2, ha='left', va='center', style='normal')

    hollow_arrow(ax, (13.9, 4.60), (9.75, 6.55), color=BLUE_GRAY, lw=1.6, rad=0.28)
    note(ax, 14.05, 5.70, 'Titration of\nfree BM3R1', color=BLUE_GRAY,
         fontsize=9.2, ha='left', va='center', style='normal')
    note(ax, 14.05, 5.15, 'Not yet modelled', color=BLUE_GRAY, fontsize=9.2,
         ha='left', va='center', style='normal')

    note(ax, 11.3, 3.05, 'Timer input: residual BM3R1\nafter the PB \u2192 LR transition',
         color=DEEP_BLUE, fontsize=9.2, ha='center', va='top', style='normal')

    note(ax, 7.6, 7.95, 'BM3R1 shared interface', color=DEEP_BLUE, fontsize=15,
         ha='center', va='center', style='normal')
    legend_lines(ax, 0.35, 7.60,
                 [('repression', 'Repression'), ('hollow', 'Prospective coupling')],
                 fontsize=9.0, dy=0.48)
    note(ax, 0.35, 1.85,
         '\u2022  The shared molecule is BM3R1, not the integrase or the RDF pool\n'
         '\u2022  Titration of free BM3R1 by the shutdown promoter is not yet modelled\n'
         '\u2022  Whether both chapters use the same BM3R1 sub-model is not yet verified\n'
         '\u2022  Prospective coupling requiring joint model validation',
         color=INK, fontsize=9.8, va='top')
    return save(fig, OUT / 'fig04_6_bm3r1_shutdown_interface.png')


# --------------------------------------------------------------------------
# 5. Full-system coupling
# --------------------------------------------------------------------------
def fig07_full_system_coupling():
    fig, ax = new_canvas(15.0, 8.4)

    box(ax, 6.5, 7.30, 7.6, 0.95,
        'Clock module  \u2014  C31 translation waveform  (period 10.59 h)',
        fc=LIGHT_BLUE, ec=DEEP_BLUE, fontsize=11, bold=True)

    bits_x = [1.9, 6.5, 11.1]
    carry_x = [4.2, 8.8]
    for i, cx in enumerate(bits_x):
        label = f'Bit {i}\nInt / RDF / BM3R1' if i < 2 else 'Bit 2\nInt / RDF / BM3R1'
        box(ax, cx, 4.60, 2.1, 1.15, label, fc=LIGHT_GRAY, ec=DEEP_BLUE,
            fontsize=9.6, bold=False)
    for i, cx in enumerate(carry_x):
        box(ax, cx, 4.60, 2.0, 0.95, f'Carry {i}\nI1-FFL', fc=LIGHT_TEAL,
            ec=TEAL, fontsize=9.4)

    # shared clock bus
    ybus = 6.35
    ax.plot([bits_x[0], bits_x[-1]], [ybus, ybus], color=DEEP_BLUE, lw=1.8, zorder=2)
    arrow(ax, (6.5, 6.80), (6.5, ybus), color=DEEP_BLUE, lw=1.8)
    for cx in bits_x:
        arrow(ax, (cx, ybus), (cx, 5.20), color=DEEP_BLUE, lw=1.6)
    note(ax, 3.2, 6.50, 'shared input', color=DEEP_BLUE, fontsize=9.0,
         ha='center', va='bottom', style='normal')

    for a, b in [(2.95, 3.20), (5.20, 5.50), (7.50, 7.80), (9.80, 10.05)]:
        arrow(ax, (a, 4.60), (b, 4.60), color=DEEP_BLUE, lw=1.8)

    box(ax, 11.1, 2.30, 5.6, 1.10,
        'Output duration module\nDelayed sRNA shutdown triggered by BM3R1',
        fc=LIGHT_ORANGE, ec=ORANGE, fontsize=10)
    arrow(ax, (11.1, 4.00), (11.1, 2.90), color=ORANGE, lw=2.0)
    note(ax, 11.35, 3.45, 'residual BM3R1 (timer)', color=ORANGE, fontsize=9.2,
         ha='left', va='center', style='normal')

    scope(ax, 0.55, 3.80, 13.0, 1.95, 'Counter (this work)')

    note(ax, 7.5, 8.05, 'Full-system coupling', color=DEEP_BLUE, fontsize=15,
         ha='center', va='center', style='normal')
    note(ax, 0.35, 1.85,
         '\u2022  Modules are unidirectional: there is no connection from a later bit back to an earlier one\n'
         '\u2022  The clock is a shared input, not a relayed one\n'
         '\u2022  The output module reads the shared BM3R1 pool; it does not modify the integrase or RDF pools',
         color=INK, fontsize=9.8, va='top')
    return save(fig, OUT / 'fig07_full_system_coupling.png')


# --------------------------------------------------------------------------
# 6. Modeling-guided design panel
# --------------------------------------------------------------------------
def fig08_design_panel():
    fig, ax = new_canvas(14.4, 8.6)
    rows = [
        ('Clock period', '10.59 h (given by the clock module)', 'given', LIGHT_BLUE, DEEP_BLUE),
        ('Recombinase dose per cycle', '\u2248 3.9 a.u. (peak-to-peak integral)',
         'measured interface', LIGHT_BLUE, DEEP_BLUE),
        ('Concentration-scale factor\n(a.u. \u2192 \u00b5M)',
         'finite, non-convex window around the centre;\nno single \u00b1% figure before a 2D scan',
         'calibration required', LIGHT_ORANGE, ORANGE),
        ('Carry maturation half-life', '30\u201335 min (centre 32.5 min)',
         'certified', LIGHT_TEAL, TEAL),
        ('Carry mRNA half-life', '1\u20132 min', 'certified', LIGHT_TEAL, TEAL),
        ('Gate quality',
         'peak ratio L_peak plus the integrated\nleakage dose outside the main window',
         'acceptance metric', LIGHT_GRAY, BLUE_GRAY),
        ('Level-three counting',
         'not achieved; expression-time tuning is saturated\n(best: 4 of 12 carry writes)',
         'limitation', LIGHT_ORANGE, RISK_ORANGE),
    ]
    x0, x1 = 0.35, 14.05
    c1, c2 = 5.35, 10.45
    ytop = 7.35
    dy = 0.92
    note(ax, x0 + 0.05, ytop + 0.42,
         'Modeling-guided design rules', color=DEEP_BLUE, fontsize=15,
         ha='left', va='center', style='normal')
    for cx, head in ((x0, 'Design quantity'), (c1, 'Recommended range'),
                     (c2, 'Evidence strength')):
        note(ax, cx + 0.1, ytop + 0.02, head, color=DEEP_BLUE, fontsize=10.5,
             ha='left', va='center', style='normal')
        ax.plot([x0, x1], [ytop - 0.22, ytop - 0.22], color=DEEP_BLUE, lw=1.2)
    for i, (q, r, tag, fc, ec) in enumerate(rows):
        yc = ytop - 0.62 - i * dy
        if i % 2 == 1:
            ax.add_patch(plt_rect(x0, yc - dy / 2 + 0.06, x1 - x0, dy - 0.12))
            ax.patches[-1].set_facecolor(LIGHT_GRAY)
            ax.patches[-1].set_edgecolor('none')
            ax.patches[-1].set_zorder(0)
        note(ax, x0 + 0.1, yc, q, color=INK, fontsize=10, ha='left',
             va='center', style='normal')
        note(ax, c1 + 0.1, yc, r, color=INK, fontsize=9.4, ha='left',
             va='center', style='normal')
        box(ax, c2 + 1.5, yc, 3.1, 0.52, tag, fc=fc, ec=ec, fontsize=9.2,
            rounding=0.1)
    ax.plot([x1, x1], [ytop - 0.22, ytop - 0.62 - (len(rows) - 1) * dy - 0.45],
            color=BLUE_GRAY, lw=0.8, ls=(0, (4, 4)))
    note(ax, x0 + 0.05, 0.78,
         'Evidence strength: certified = verified by the frozen acceptance criteria; '
         'calibration required = uncalibrated factor that the design depends on;\n'
         'limitation = reported as a boundary, not as a success.',
         color=BLUE_GRAY, fontsize=9.2, va='top')
    return save(fig, OUT / 'fig08_design_panel.png')


FIGURES = {
    'fig01_modeling_workflow': fig01_modeling_workflow,
    'fig03_single_bit_mechanism': fig03_single_bit_mechanism,
    'fig04_ffl_carry': fig04_ffl_carry,
    'fig04_6_bm3r1_shutdown_interface': fig04_6_bm3r1_shutdown_interface,
    'fig07_full_system_coupling': fig07_full_system_coupling,
    'fig08_design_panel': fig08_design_panel,
}

OUT = Path(__file__).resolve().parent / 'out'


def main():
    global OUT
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=Path, default=OUT)
    ap.add_argument('--only', default=None,
                    help='comma-separated figure names, default: all')
    args = ap.parse_args()
    OUT = args.out
    names = list(FIGURES) if not args.only else [s.strip() for s in args.only.split(',')]
    for name in names:
        if name not in FIGURES:
            raise SystemExit(f'unknown figure: {name}')
        path = FIGURES[name]()
        print(f'wrote {path}')


if __name__ == '__main__':
    main()
