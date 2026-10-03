"""Figures for the version probe (B) and the hybrid boundary probes (P1-P4).

NEW FILE (added 2026-09-28).  Reads ONLY the two JSON files this package
produced; it does not re-run a model and does not touch any existing artefact.

Style follows the project-wide visual specification: tempo_style palette,
English in-figure labels only (no CJK font dependency in headless rendering),
SVG written with `svg.fonttype = 'none'` so the text stays real text.

TWO TRAPS THIS SCRIPT DELIBERATELY AVOIDS
  1. In the failed version-probe cases bit1 and bit2 never cross 0.5, so no
     setup/hold margin EXISTS for them.  A "binding margin" aggregated over bits
     would then silently report bit0's margin and look healthy.  Panel B plots
     the per-bit margins and marks the missing ones "no crossing" instead.
  2. The 51-state frozen gate-duty targets (0.305 / 0.475) belong to the 51-state
     gate signal.  This hybrid's g1 is a different quantity on a different scale
     (duty ~0.018), so those lines are NOT drawn as if they were a target here.

Usage:  python make_probe_figures.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')

HERE = Path(__file__).resolve().parent
RES = HERE / 'results'
FIG = HERE / 'figures'
sys.path.insert(0, str(Path(r'C:\Users\18633\Desktop\wiki\dshwork\figures')))

from tempo_style import BLUE_GRAY, DEEP_BLUE, INK, LIGHT_ORANGE, LIGHT_TEAL  # noqa: E402
from tempo_style import ORANGE, RISK_ORANGE, TEAL                              # noqa: E402

SVG_FONTTYPE = 'none'
DPI = 200
CHECKS = []


def check(ok, msg):
    CHECKS.append((bool(ok), msg))
    print(('  OK   ' if ok else '  FAIL ') + msg, flush=True)


def load(name):
    p = RES / name
    if not p.exists():
        raise SystemExit(f'missing input: {p}')
    return json.loads(p.read_text(encoding='utf-8'))


def save(fig, stem):
    FIG.mkdir(parents=True, exist_ok=True)
    outs = []
    for ext in ('png', 'svg', 'pdf'):
        p = FIG / f'{stem}.{ext}'
        fig.savefig(p, dpi=DPI, bbox_inches='tight', facecolor=fig.get_facecolor())
        outs.append(p)
    plt.close(fig)
    return outs


def style(ax):
    ax.grid(axis='y', ls=':', alpha=0.5)
    ax.set_axisbelow(True)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    ax.tick_params(labelsize=8)


def arms_panel(ax, cols, arms, values, title):
    # values is one list per CASE; imshow wants one row per ARM.
    m = np.array([[1.0 if v else 0.0 for v in row] for row in values]).T
    if m.shape != (len(arms), len(cols)):
        raise ValueError(f'arm matrix is {m.shape}, expected {(len(arms), len(cols))}')
    ax.imshow(m, cmap=mcolors.ListedColormap([LIGHT_ORANGE, LIGHT_TEAL]),
              vmin=0, vmax=1, aspect='auto')
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(cols, rotation=25, ha='right', fontsize=8)
    ax.set_yticks(range(len(arms)))
    ax.set_yticklabels(arms, fontsize=8)
    for i in range(len(arms)):
        for j in range(len(cols)):
            ax.text(j, i, 'pass' if m[i, j] else 'FAIL', ha='center', va='center',
                    fontsize=7.5, color=INK if m[i, j] else RISK_ORANGE,
                    fontweight='bold' if not m[i, j] else 'normal')
    ax.set_xticks(np.arange(-.5, len(cols), 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(arms), 1), minor=True)
    ax.grid(which='minor', color='white', lw=1.4)
    ax.tick_params(which='minor', length=0)
    ax.set_title(title, fontsize=9.5, color=DEEP_BLUE, pad=8)


# --------------------------------------------------------------- figure A
def figure_version_probe():
    d = load('probe_B_zeng_version_51state.json')
    rows, verdicts = d['rows'], d['verdicts']
    check(len(rows) == 3, f'version probe has 3 cases (found {len(rows)})')
    check(rows[0]['certified'], 'version probe baseline is certified')
    check(abs(rows[0]['margins']['binding_min_hold_h'] - 2.913722775395758) < 1e-9,
          'baseline bit1 hold margin reproduces 2.913722775395758 h')
    check(rows[0]['carry1_gate_events'] == 14 and rows[0]['bit2_crossings'] == 14,
          'baseline has 14 carry-1 gate events and 14 bit2 crossings')
    check(all(r['carry1_gate_events'] == 0 for r in rows[1:]),
          'both current-code variants produce ZERO carry-1 gate events')
    gaps = d['restore_gap_after_each_run']
    check(max(gaps.values()) == 0.0,
          f'ZENG restored bit-identically after every run (max gap {max(gaps.values()):g})')

    tags = [r['tag'].replace('_', ' ') for r in rows]
    fig, axes = plt.subplots(1, 3, figsize=(15.6, 4.8), constrained_layout=True)
    fig.patch.set_facecolor('white')

    arms = ['steady (mod 8)', 'one-to-one', 'causal order', 'alternating', 'gate contrast']
    vals = [[r['steady_passed'], r['exactly_one'], r['causal_order'], r['alternating'],
             r['contrast_passed']] for r in rows]
    arms_panel(axes[0], tags, arms, vals, 'A  Frozen 51-state predicate arms')

    w, xs = 0.26, np.arange(len(rows))
    for k, (bit, off, col) in enumerate((('bit0', -w, BLUE_GRAY), ('bit1', 0.0, TEAL),
                                         ('bit2', w, DEEP_BLUE))):
        vals_b = []
        for r in rows:
            v = r['margins'].get(bit, {}).get('min_hold_h')
            vals_b.append(np.nan if v is None else v)
        bars = axes[1].bar(xs + off, [0 if np.isnan(v) else v for v in vals_b], w,
                           label=f'{bit} hold margin', color=col, edgecolor=DEEP_BLUE, lw=0.7)
        for rect, v in zip(bars, vals_b):
            if np.isnan(v):
                axes[1].text(rect.get_x() + rect.get_width() / 2, 0.06, 'no\ncrossing',
                             ha='center', va='bottom', fontsize=6.5, color=RISK_ORANGE)
            else:
                axes[1].text(rect.get_x() + rect.get_width() / 2, v, f'{v:.2f}', ha='center',
                             va='bottom', fontsize=7, color=INK)
    axes[1].axhline(2.913722775395758, color=ORANGE, ls='--', lw=1.2)
    axes[1].text(len(rows) - 0.5, 2.913722775395758, ' baseline binding 2.914 h',
                 ha='right', va='bottom', fontsize=7.5, color=ORANGE)
    axes[1].set_xticks(xs)
    axes[1].set_xticklabels(tags, rotation=25, ha='right', fontsize=8)
    axes[1].set_ylabel('min hold margin (h)', fontsize=9)
    axes[1].set_ylim(0, 5.4)
    axes[1].legend(fontsize=7, frameon=False, ncol=3, loc='upper left')
    axes[1].set_title('B  Per-bit hold margin (missing = no crossing at all)',
                      fontsize=9.5, color=DEEP_BLUE, pad=8)
    style(axes[1])

    g1 = [verdicts[r['tag']]['signal_ranges']['g1'][1] for r in rows]
    bars = axes[2].bar(xs, g1, 0.55, color=[TEAL] + [ORANGE] * (len(rows) - 1),
                       edgecolor=DEEP_BLUE, lw=0.9)
    axes[2].set_yscale('log')
    axes[2].axhline(0.05, color=RISK_ORANGE, ls='--', lw=1.4)
    axes[2].text(len(rows) - 0.5, 0.05, ' GATE_ON = 0.05', ha='right', va='bottom',
                 fontsize=7.5, color=RISK_ORANGE)
    for i, (r, v) in enumerate(zip(rows, g1)):
        axes[2].text(i, v, f'{v:.2e}\n{r["carry1_gate_events"]} gate\nevents', ha='center',
                     va='bottom', fontsize=7, color=INK)
    axes[2].set_xticks(xs)
    axes[2].set_xticklabels(tags, rotation=25, ha='right', fontsize=8)
    axes[2].set_ylabel('carry-1 gate peak (a.u.)', fontsize=9)
    axes[2].set_ylim(1e-4, 1.2)
    axes[2].set_title('C  The gate never opens under the current code values',
                      fontsize=9.5, color=DEEP_BLUE, pad=8)
    style(axes[2])

    fig.suptitle('Version probe B: frozen 51-state working point under the current code '
                 'parameters (600 h, clock 0.3/2.0, n_A1_gate = 6)',
                 fontsize=10.5, color=DEEP_BLUE)
    outs = save(fig, 'probe_B_version_51state')
    check(all(p.exists() and p.stat().st_size > 8000 for p in outs),
          'figure A written as PNG + SVG + PDF with non-trivial size')
    svg = [p for p in outs if p.suffix == '.svg'][0].read_text(encoding='utf-8')
    check('<text' in svg, "figure A SVG uses <text> (svg.fonttype = 'none')")


# --------------------------------------------------------------- figure B
def figure_hybrid_probes():
    d = load('probe_P1_P4_hybrid.json')
    rows, inv = d['rows'], d['inventories']
    check(len(rows) == 7, f'hybrid probe has 7 cases (found {len(rows)})')
    check(rows[0]['certified_v1'], 'hybrid probe baseline is certified_v1')
    check(rows[0]['steady_reads'] == 19,
          f'baseline has 19 steady read windows (found {rows[0]["steady_reads"]})')
    check(abs(rows[0]['global_min_timing_margin_h'] - 1.178460901856866) < 1e-12,
          'baseline global min timing margin reproduces the frozen 1.178460901856866 h')
    check(all(r['certified_v1'] for r in rows),
          'every pre-registered boundary probe still certifies')
    check(all(r['boundary_clips'] == 0 for r in rows),
          'no probe clips a read window (failure mode is not a sampling artefact)')

    tags = [r['tag'].replace('_', ' ') for r in rows]
    fig, axes = plt.subplots(1, 4, figsize=(21.0, 5.0), constrained_layout=True)
    fig.patch.set_facecolor('white')

    arms = ['steady (mod 8)', 'chain bit0 -> bit1', 'chain bit1 -> bit2']
    vals = [[r['steady_passed'], r['events']['bit0_to_bit1']['passed'],
             r['events']['bit1_to_bit2']['passed']] for r in rows]
    arms_panel(axes[0], tags, arms, vals, 'A  HZH_MOD8_CAUSAL_V1 arms')

    read = [r['steady_reads'] for r in rows]
    clips = [r['boundary_clips'] for r in rows]
    axes[1].bar(range(len(rows)), read, 0.6, color=[TEAL] + [ORANGE] * (len(rows) - 1),
                edgecolor=DEEP_BLUE, lw=0.9)
    axes[1].axhline(read[0], color=BLUE_GRAY, ls='--', lw=1.2)
    for i, (v, c) in enumerate(zip(read, clips)):
        axes[1].text(i, v, f'{v}\n{int(c)} clips', ha='center', va='bottom', fontsize=7.5,
                     color=INK)
    axes[1].set_xticks(range(len(rows)))
    axes[1].set_xticklabels(tags, rotation=25, ha='right', fontsize=8)
    axes[1].set_ylabel('steady read windows / 300 h', fontsize=9)
    axes[1].set_ylim(0, max(read) * 1.28)
    axes[1].set_title('B  Read windows and boundary clips', fontsize=9.5, color=DEEP_BLUE, pad=8)
    style(axes[1])

    duty = [r['duty_g1_gt_0p1'] for r in rows]
    axes[2].bar(range(len(rows)), duty, 0.6, color=[TEAL] + [ORANGE] * (len(rows) - 1),
                edgecolor=DEEP_BLUE, lw=0.9)
    for i, v in enumerate(duty):
        axes[2].text(i, v, f'{v:.4f}', ha='center', va='bottom', fontsize=7.5, color=INK)
    axes[2].set_xticks(range(len(rows)))
    axes[2].set_xticklabels(tags, rotation=25, ha='right', fontsize=8)
    axes[2].set_ylabel('duty(g1 > 0.1)', fontsize=9)
    axes[2].set_ylim(0, max(duty) * 1.45)
    axes[2].set_title('C  Gate duty (per-case comparison only)', fontsize=9.5,
                      color=DEEP_BLUE, pad=8)
    axes[2].text(0.02, 0.97, 'no inherited target line:\nthe 51-state duty targets\n(0.305 / 0.475)'
                             '\nbelong to a different gate signal',
                 transform=axes[2].transAxes, fontsize=6.8, color=BLUE_GRAY, va='top')
    style(axes[2])

    pools = ('b2_I', 'F1')
    w, xs = 0.36, np.arange(len(rows))
    for k, (pool, col) in enumerate(zip(pools, (TEAL, DEEP_BLUE))):
        vals_p = [inv[r['tag']]['per_pool'][pool]['gate_phase_minimum']['copies'] for r in rows]
        axes[3].bar(xs + (k - 0.5) * w, vals_p, w, label=f'{pool} gate-phase minimum',
                    color=col, edgecolor=DEEP_BLUE, lw=0.8)
        for i, v in enumerate(vals_p):
            axes[3].text(i + (k - 0.5) * w, v, f'{v:.2f}', ha='center', va='bottom',
                         fontsize=6.8, color=INK)
    axes[3].set_yscale('log')
    axes[3].axhline(1.0, color=BLUE_GRAY, ls='--', lw=1.2)
    axes[3].axhline(100.0, color=RISK_ORANGE, ls='--', lw=1.4)
    axes[3].axhline(20.0, color=ORANGE, ls='-.', lw=1.2)
    axes[3].text(len(rows) - 0.5, 1.0, ' 1 copy', ha='right', va='bottom', fontsize=7,
                 color=BLUE_GRAY)
    axes[3].text(len(rows) - 0.5, 100.0, ' 100 copies: >= 10 % intrinsic noise', ha='right',
                 va='bottom', fontsize=7, color=RISK_ORANGE)
    axes[3].text(len(rows) - 0.5, 20.0, ' 20 copies: the frozen 51-state b2_I minimum',
                 ha='right', va='bottom', fontsize=7, color=ORANGE)
    axes[3].set_xticks(xs)
    axes[3].set_xticklabels(tags, rotation=25, ha='right', fontsize=8)
    axes[3].set_ylabel('gate-phase minimum (copies)', fontsize=9)
    axes[3].set_ylim(0.5, 300)
    axes[3].legend(fontsize=7, frameon=False, loc='upper left')
    axes[3].set_title('D  Smallest working pool (P4): b2_R and b2_C are below 1 copy, '
                      'so they carry no CV', fontsize=9.5, color=DEEP_BLUE, pad=8)
    style(axes[3])

    fig.suptitle('Boundary probes P1-P4 on the frozen 34-state hybrid (Han v53d upstream, '
                 '300 h, clock from this model\'s own Int0)', fontsize=10.5, color=DEEP_BLUE)
    outs = save(fig, 'probe_P1_P4_hybrid')
    check(all(p.exists() and p.stat().st_size > 8000 for p in outs),
          'figure B written as PNG + SVG + PDF with non-trivial size')
    svg = [p for p in outs if p.suffix == '.svg'][0].read_text(encoding='utf-8')
    check('<text' in svg, "figure B SVG uses <text> (svg.fonttype = 'none')")


def main():
    plt.rcParams['svg.fonttype'] = SVG_FONTTYPE
    plt.rcParams['font.family'] = 'DejaVu Sans'
    print('Figure A: version probe B')
    figure_version_probe()
    print('Figure B: hybrid boundary probes')
    figure_hybrid_probes()
    bad = [m for ok, m in CHECKS if not ok]
    print(f'\n{len(CHECKS) - len(bad)}/{len(CHECKS)} checks passed')
    if bad:
        raise SystemExit('FAILED CHECKS:\n' + '\n'.join(bad))


if __name__ == '__main__':
    main()
