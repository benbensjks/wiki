"""Figure B2: where the carry-1 gate pulse actually dies.

NEW FILE (added 2026-09-29).  Reads only `results/probe_B2_gate_decomposition.json`
and writes `figures/probe_B2_gate_decomposition.{png,svg,pdf}`.

Style matches `make_probe_figures.py`: tempo_style palette, English in-figure
labels only, SVG with `svg.fonttype = 'none'`.

WHAT IT SHOWS (and why the obvious reading is wrong)
  The gate is g0 = act(A0) * rep(F0).  The first-round report said the pulse
  collapsed; the natural guess is "the activator arm weakened".  Panel A shows
  the opposite: act(A0) at the gate peak barely moves and can even RISE, while
  rep(F0) collapses by 13x and then 2477x, and the product reproduces the
  measured g0 peak exactly.  Panel B shows why: F0 is never cleared, because A0
  no longer returns to ~0.  Panel C shows the consequence: in the baseline the
  gate peak fires 0.43 h after the bit0 reset; in the perturbed cases it fires
  20.6-30.8 h after it, i.e. it is no longer reset-triggered at all.

Usage:  python make_figure_B2_gate.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')

HERE = Path(__file__).resolve().parent
RES = HERE / 'results'
FIG = HERE / 'figures'
sys.path.insert(0, str(Path(r'C:\Users\18633\Desktop\wiki\dshwork\figures')))

from tempo_style import BLUE_GRAY, DEEP_BLUE, INK, ORANGE, RISK_ORANGE, TEAL  # noqa: E402

SVG_FONTTYPE = 'none'
DPI = 200
CHECKS = []
LABELS = {'baseline_frozen_zeng': 'baseline\n(frozen ZENG)',
          'zmh_bit1_block': 'zmh bit1 block',
          'zmh_bit0_block': 'zmh bit0 block'}


def check(ok, msg):
    CHECKS.append((bool(ok), msg))
    print(('  OK   ' if ok else '  FAIL ') + msg, flush=True)


def style(ax):
    ax.grid(axis='y', ls=':', alpha=0.5)
    ax.set_axisbelow(True)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    ax.tick_params(labelsize=8)


def main():
    plt.rcParams['svg.fonttype'] = SVG_FONTTYPE
    plt.rcParams['font.family'] = 'DejaVu Sans'
    d = json.loads((RES / 'probe_B2_gate_decomposition.json').read_text(encoding='utf-8'))
    rows = d['rows']
    check(len(rows) == 3, f'three cases (found {len(rows)})')
    tags = [r['tag'] for r in rows]
    xs = np.arange(len(rows))
    names = [LABELS.get(t, t) for t in tags]

    # the product identity is the load-bearing check of this whole figure
    for r in rows:
        gp = r['gate_peak']
        prod = gp['act_of_A0'] * gp['rep_of_F0']
        check(abs(prod - gp['g0']) <= 1e-12 * max(1.0, abs(gp['g0'])),
              f"{r['tag']}: act(A0) x rep(F0) reproduces the measured g0 peak "
              f"({prod:.8g} vs {gp['g0']:.8g})")
    check(rows[0]['per_cycle_gate_peaks'][0]['delay_after_reset_h'] < 1.0,
          'baseline gate peak fires within 1 h of the reset (reset-triggered pulse)')
    check(all(max(p['delay_after_reset_h'] for p in r['per_cycle_gate_peaks']) > 10.0
              for r in rows[1:]),
          'both perturbed cases fire the gate peak >10 h after the reset (no longer triggered)')

    fig, axes = plt.subplots(1, 3, figsize=(16.2, 5.0), constrained_layout=True)
    fig.patch.set_facecolor('white')

    w = 0.26
    series = [('act(A0) at gate peak', 'act_of_A0', TEAL),
              ('rep(F0) at gate peak', 'rep_of_F0', ORANGE),
              ('product = g0 peak', 'g0', DEEP_BLUE)]
    for k, (label, key, col) in enumerate(series):
        vals = [r['gate_peak'][key] for r in rows]
        bars = axes[0].bar(xs + (k - 1) * w, vals, w, label=label, color=col,
                           edgecolor=DEEP_BLUE, lw=0.8)
        for rect, v in zip(bars, vals):
            axes[0].text(rect.get_x() + rect.get_width() / 2, v, f'{v:.3g}', ha='center',
                         va='bottom', fontsize=6.8, color=INK)
    axes[0].set_yscale('log')
    axes[0].set_ylim(1e-4, 8)
    axes[0].set_xticks(xs)
    axes[0].set_xticklabels(names, fontsize=8)
    axes[0].set_ylabel('value at the gate peak', fontsize=9)
    axes[0].legend(fontsize=7.2, frameon=False, loc='lower left')
    axes[0].set_title('A  The collapse is in the repressor arm, not the activator arm',
                      fontsize=9.5, color=DEEP_BLUE, pad=8)
    style(axes[0])

    for i, r in enumerate(rows):
        lo, hi = r['factor_extremes']['F0']
        axes[1].plot([i, i], [lo, hi], color=BLUE_GRAY, lw=2.4, zorder=2)
        axes[1].plot([i], [lo], 'o', color=RISK_ORANGE, ms=8, zorder=3)
        axes[1].text(i + 0.12, lo, f'min {lo:.3f}', fontsize=7.5, color=RISK_ORANGE,
                     va='center')
        axes[1].plot([i], [r['gate_peak']['F0']], 's', color=DEEP_BLUE, ms=6, zorder=3)
        axes[1].text(i + 0.12, r['gate_peak']['F0'], f'at peak {r["gate_peak"]["F0"]:.3f}',
                     fontsize=7.5, color=DEEP_BLUE, va='center')
    axes[1].axhline(0.6, color=ORANGE, ls='--', lw=1.3)
    axes[1].text(len(rows) - 0.55, 0.6, ' K_F[0] = 0.6', fontsize=7.5, color=ORANGE,
                 va='bottom')
    axes[1].set_xticks(xs)
    axes[1].set_xticklabels(names, fontsize=8)
    axes[1].set_xlim(-0.5, len(rows) - 0.35)
    axes[1].set_ylabel('F0 (a.u.)', fontsize=9)
    axes[1].set_title('B  F0 is never cleared: its minimum rises 0.43 -> 4.64 a.u.',
                      fontsize=9.5, color=DEEP_BLUE, pad=8)
    style(axes[1])

    for i, r in enumerate(rows):
        dl = [p['delay_after_reset_h'] for p in r['per_cycle_gate_peaks']]
        axes[2].plot(np.full(len(dl), i) + np.linspace(-0.13, 0.13, len(dl)), dl, 'o',
                     color=TEAL if i == 0 else ORANGE, ms=4.5, alpha=0.85, zorder=3)
        axes[2].plot([i - 0.2, i + 0.2], [np.median(dl)] * 2, color=DEEP_BLUE, lw=2.2,
                     zorder=4)
        axes[2].text(i + 0.22, np.median(dl), f'median {np.median(dl):.2f} h', fontsize=7.5,
                     color=INK, va='center')
    axes[2].axhline(1.0, color=BLUE_GRAY, ls='--', lw=1.2)
    axes[2].text(len(rows) - 0.55, 1.0, ' 1 h: a reset-triggered pulse would sit here',
                 fontsize=7.5, color=BLUE_GRAY, va='bottom')
    axes[2].set_xticks(xs)
    axes[2].set_xticklabels(names, fontsize=8)
    axes[2].set_xlim(-0.5, len(rows) - 0.3)
    axes[2].set_ylabel('gate-peak delay after the bit0 reset (h)', fontsize=9)
    axes[2].set_title('C  The pulse stops being reset-triggered at all',
                      fontsize=9.5, color=DEEP_BLUE, pad=8)
    style(axes[2])

    fig.suptitle('Probe B2: decomposing the carry-1 gate collapse '
                 '(600 h, frozen 51-state working point)', fontsize=10.5, color=DEEP_BLUE)
    FIG.mkdir(parents=True, exist_ok=True)
    outs = []
    for ext in ('png', 'svg', 'pdf'):
        p = FIG / f'probe_B2_gate_decomposition.{ext}'
        fig.savefig(p, dpi=DPI, bbox_inches='tight', facecolor=fig.get_facecolor())
        outs.append(p)
    plt.close(fig)
    check(all(p.exists() and p.stat().st_size > 8000 for p in outs),
          'figure B2 written as PNG + SVG + PDF with non-trivial size')
    svg = [p for p in outs if p.suffix == '.svg'][0].read_text(encoding='utf-8')
    check('<text' in svg, "figure B2 SVG uses <text> (svg.fonttype = 'none')")

    bad = [m for ok, m in CHECKS if not ok]
    print(f'\n{len(CHECKS) - len(bad)}/{len(CHECKS)} checks passed')
    if bad:
        raise SystemExit('FAILED CHECKS:\n' + '\n'.join(bad))


if __name__ == '__main__':
    main()
