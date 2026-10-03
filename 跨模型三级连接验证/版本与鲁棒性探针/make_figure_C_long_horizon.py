"""Figures for probe C: the frozen hybrid over 1200 h, plus a mid-run zoom.

NEW FILE (added 2026-09-29).  Reads only what probe C produced
(`results/probe_C_long_horizon.json` and `results/probe_C_long_horizon_series.npz`)
and writes two figures:

  figures/probe_C_long_horizon_overview.{png,svg,pdf}   -- the whole 1200 h
  figures/probe_C_midwindow_zoom.{png,svg,pdf}          -- 42 h around t = 600 h

Conventions follow the established data-figure scripts in `dshwork/figures/`:
tempo_style palette, English in-figure labels only, an annotation band above the
plot rows and a footer band below them so no annotation lands on a curve,
PNG + SVG + PDF at 200 dpi, and SVG written with `svg.fonttype = 'none'`.

Usage:  python make_figure_C_long_horizon.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.transforms import Bbox

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')

HERE = Path(__file__).resolve().parent
RES = HERE / 'results'
FIG = HERE / 'figures'
sys.path.insert(0, str(Path(r'C:\Users\18633\Desktop\wiki\dshwork\figures')))

from tempo_style import (BLUE_GRAY, DEEP_BLUE, INK, LIGHT_BLUE,  # noqa: E402
                         LIGHT_GRAY, ORANGE, RISK_ORANGE, TEAL)

SVG_FONTTYPE = 'none'
DPI = 200
BAND_LOW, BAND_HIGH = 0.30, 0.70
GATE_ON = 0.05
JREV_THRESHOLD = 0.10
HORIZON_MARKS = (300.0, 600.0)

OVERVIEW_LABELS = (
    'C31 mRNA (shared, upstream)',
    'C31 translation flux',
    'DNA state', '(LR fraction)',
    'Low band', 'High band',
    'Decoded read', 'window value',
    'Integrase pools',
    'Time (h)',
    '300 h', '600 h',
)
ZOOM_LABELS = (
    'DNA state', '(LR fraction)',
    'Integrase pools',
    'Carry gates',
    'Recombination', 'flux (per h)',
    'read window', 'reset event', 'gate open', 'bit-2 flip',
    'Time (h)',
)
FOOTER_PREFIXES = ('Model:', 'Input:', 'Horizon:')
CHECKS = []


def check(ok, msg):
    CHECKS.append((bool(ok), msg))
    print(('  OK   ' if ok else '  FAIL ') + msg, flush=True)


def style(ax):
    ax.grid(alpha=0.18)
    for side in ('top',):
        ax.spines[side].set_visible(False)
    ax.tick_params(labelsize=8.5)


def normalised(v):
    v = np.asarray(v, dtype=float)
    peak = float(np.max(np.abs(v)))
    return v / peak if peak > 0 else v


def bands(ax):
    ax.axhspan(0.0, BAND_LOW, color=TEAL, alpha=0.07, lw=0)
    ax.axhspan(BAND_HIGH, 1.0, color=DEEP_BLUE, alpha=0.07, lw=0)


def load():
    js = json.loads((RES / 'probe_C_long_horizon.json').read_text(encoding='utf-8'))
    npz = np.load(RES / 'probe_C_long_horizon_series.npz', allow_pickle=False)
    return js, {k: npz[k] for k in npz.files}


def band_annotation(ax, js):
    meta = js['meta']
    row = js['horizons'][-1]
    ax.text(0.0, 1.0,
            'Frozen hybrid HBY bit0 / ZMH bit1 / HBY bit2, one continuous run',
            transform=ax.transAxes, ha='left', va='top', fontsize=13,
            color=INK, fontweight='bold')
    ax.text(0.0, 0.56,
            f"{meta['hours']:.0f} h, {meta['sample_min']:g} min grid, DOP853 "
            f"rtol {meta['solver']['rtol']:g} / atol {meta['solver']['atol']:g}   |   "
            f"clock {meta['receiver_config']['clock_K_au']:g}/"
            f"{meta['receiver_config']['clock_n']:g}, "
            f"n_A1_gate {meta['receiver_config']['n_A1_gate']:g}, "
            f"uM_per_au {meta['receiver_config']['receiver_uM_per_au']:g}",
            transform=ax.transAxes, ha='left', va='top', fontsize=9.5, color=INK)
    ax.text(0.0, 0.06,
            f"certified_v1 = {row['certified_v1']} at every probed horizon   |   "
            f"steady read windows {row['steady']['reads']} (mod-8 increments "
            f"{row['steady']['increments_mod8']}, boundary clips "
            f"{row['steady']['boundary_clips']})   |   min timing margin "
            f"{row['global_min_timing_margin_h']:.4f} h",
            transform=ax.transAxes, ha='left', va='bottom', fontsize=9.5, color=TEAL)


def footer(ax, js, title):
    meta = js['meta']
    row = js['horizons'][-1]
    ax.text(0.0, 1.0, f'Model: {title}', transform=ax.transAxes, ha='left', va='top',
            fontsize=8.5, color=INK)
    ax.text(0.0, 0.60,
            'Input: Han v53d upstream, true integration to the full horizon '
            '(no end-point interpolation, normalisation window unchanged)',
            transform=ax.transAxes, ha='left', va='top', fontsize=8.5, color=INK)
    ax.text(0.0, 0.20,
            f"Horizon: {meta['sub_horizons'][0]:.0f} / {meta['sub_horizons'][1]:.0f} / "
            f"{meta['sub_horizons'][2]:.0f} h from a single run; the 300 h slice "
            f"reproduces the frozen certification exactly",
            transform=ax.transAxes, ha='left', va='top', fontsize=8.5, color=INK)


# ------------------------------------------------------------------ figure 1
def overview(js, S):
    t = S['ov_time']
    fig = plt.figure(figsize=(13.5, 12.2), constrained_layout=True)
    gs = fig.add_gridspec(6, 1, height_ratios=[0.42, 1.0, 1.0, 0.85, 0.85, 0.26],
                          hspace=0.34)
    axN = fig.add_subplot(gs[0]); axN.axis('off')
    axA = fig.add_subplot(gs[1])
    axB = fig.add_subplot(gs[2], sharex=axA)
    axC = fig.add_subplot(gs[3], sharex=axA)
    axD = fig.add_subplot(gs[4], sharex=axA)
    axE = fig.add_subplot(gs[5]); axE.axis('off')
    band_annotation(axN, js)

    axA.plot(t, S['ov_m31'], lw=1.3, color=DEEP_BLUE, label='C31 mRNA (shared, upstream)')
    axA.set_ylabel('C31 mRNA\n(copies / cell)', fontsize=9.5, color=DEEP_BLUE)
    axA.tick_params(axis='y', labelsize=8.5, colors=DEEP_BLUE)
    axA2 = axA.twinx()
    axA2.plot(t, S['ov_flux'], lw=1.2, color=ORANGE, alpha=0.9,
              label='C31 translation flux')
    axA2.set_ylabel('translation flux\n(copies / cell / min)', fontsize=9.5, color=ORANGE)
    axA2.tick_params(axis='y', labelsize=8.5, colors=ORANGE)
    axA2.spines['top'].set_visible(False)
    lines = axA.get_lines() + axA2.get_lines()
    axA.legend(lines, [l.get_label() for l in lines], ncol=2, loc='upper left',
               fontsize=8.5, framealpha=0.9)
    axA.set_title('A  Upstream over the whole horizon (a real integration, not '
                  'interpolated end points)', fontsize=10, color=DEEP_BLUE, pad=6)
    style(axA)

    for key, label, col in (('ov_S0', 'S0 (HBY bit0)', DEEP_BLUE),
                            ('ov_S1', 'S1 (ZMH bit1)', ORANGE),
                            ('ov_S2', 'S2 (HBY bit2)', TEAL)):
        axB.plot(t, S[key], lw=1.2, color=col, label=label)
    bands(axB)
    axB.set_ylim(-0.05, 1.05)
    axB.set_ylabel('DNA state\n(LR fraction)', fontsize=9.5)
    axB.legend(ncol=3, loc='upper left', fontsize=8.5, framealpha=0.9)
    axB.text(0.995, 0.06, 'Low band', transform=axB.transAxes, ha='right', va='bottom',
             fontsize=8.5, color=TEAL)
    axB.text(0.995, 0.94, 'High band', transform=axB.transAxes, ha='right', va='top',
             fontsize=8.5, color=DEEP_BLUE)
    axB.set_title('B  All three bits keep flipping for the whole run',
                  fontsize=10, color=DEEP_BLUE, pad=6)
    style(axB)

    trough, value = S['read_trough'], S['read_value']
    ok = np.isfinite(value)
    axC.step(trough, np.where(ok, value, np.nan), where='mid', lw=1.3, color=INK)
    axC.plot(trough[ok], value[ok], 'o', ms=3.0, color=INK)
    axC.plot(trough[~ok], np.full((~ok).sum(), 3.5), 'x', ms=5, color=RISK_ORANGE)
    axC.set_yticks(range(8))
    axC.set_ylim(-0.6, 7.6)
    axC.set_ylabel('Decoded read\nwindow value', fontsize=9.5)
    axC.set_title(f'C  {int(ok.sum())} of {len(value)} read windows decode; the '
                  f'mod-8 sweep never breaks', fontsize=10, color=DEEP_BLUE, pad=6)
    style(axC)

    for key, label, col in (('ov_Int0', 'Int0 (HBY bit0)', DEEP_BLUE),
                            ('ov_Int1', 'Int1 (ZMH bit1)', ORANGE),
                            ('ov_Int2', 'Int2 (HBY bit2)', TEAL)):
        axD.plot(t, S[key], lw=1.1, color=col, label=label)
    axD.set_ylabel('Integrase pools\n(a.u.)', fontsize=9.5)
    axD.legend(ncol=3, loc='upper left', fontsize=8.5, framealpha=0.9)
    axD.set_title('D  The three integrase pulses stay separated and repeating',
                  fontsize=10, color=DEEP_BLUE, pad=6)
    style(axD)
    axD.set_xlabel('Time (h)', fontsize=10)

    for ax in (axA, axB, axC, axD):
        for x in HORIZON_MARKS:
            ax.axvline(x, color=BLUE_GRAY, ls=':', lw=1.1)
        ax.set_xlim(0, js['meta']['hours'])
    for x in HORIZON_MARKS:
        axB.text(x, 1.02, f'{x:.0f} h', ha='center', va='bottom', fontsize=8.5,
                 color=BLUE_GRAY)
    footer(axE, js, 'HBY bit0 + ZMH current-code A0/F0 and bit1 + HBY A1/F1 and bit2, '
                    '34 states')
    return fig


# ------------------------------------------------------------------ figure 2
def zoom(js, S):
    t = S['zw_time']
    lo, hi = float(S['zoom'][0]), float(S['zoom'][1])
    fig = plt.figure(figsize=(13.5, 11.0), constrained_layout=True)
    gs = fig.add_gridspec(6, 1, height_ratios=[0.55, 1.0, 0.95, 1.0, 1.0, 0.24],
                          hspace=0.30)
    axN = fig.add_subplot(gs[0]); axN.axis('off')
    axA = fig.add_subplot(gs[1])
    axB = fig.add_subplot(gs[2], sharex=axA)
    axC = fig.add_subplot(gs[3], sharex=axA)
    axD = fig.add_subplot(gs[4], sharex=axA)
    axE = fig.add_subplot(gs[5]); axE.axis('off')

    row = js['horizons'][-1]
    axN.text(0.0, 1.0, f'Mid-run zoom: t = {lo:.0f} - {hi:.0f} h of the 1200 h run',
             transform=axN.transAxes, ha='left', va='top', fontsize=13, color=INK,
             fontweight='bold')
    axN.text(0.0, 0.60,
             'sampling 1 min; the same frozen candidate, no re-initialisation, '
             'no parameter change at the window boundary',
             transform=axN.transAxes, ha='left', va='top', fontsize=9.5, color=INK)
    axN.text(0.0, 0.03,
             f"read windows here: "
             f"{int(((S['read_trough'] >= lo) & (S['read_trough'] <= hi)).sum())}"
             f"   |   bit0->bit1 reverses (whole run): "
             f"{row['events']['bit0_to_bit1']['reverse_events']}"
             f"   |   bit1->bit2 reverses (whole run): "
             f"{row['events']['bit1_to_bit2']['reverse_events']}",
             transform=axN.transAxes, ha='left', va='bottom', fontsize=9.5, color=TEAL)

    idx = {name: i for i, name in enumerate(S['state_names'])}
    y = S['zw_y']
    S0, S1, S2 = y[idx['b0_S']], 1.0 - y[idx['pb1_zmh']], y[idx['b2_S']]

    inside = (S['read_start'] <= hi) & (S['read_end'] >= lo)
    for a, b in zip(S['read_start'][inside], S['read_end'][inside]):
        axA.axvspan(a, b, color=LIGHT_BLUE, alpha=0.55, lw=0, zorder=0)
    axA.text(0.995, 0.62, 'read window', transform=axA.transAxes, ha='right',
             va='top', fontsize=8.5, color=BLUE_GRAY)
    for series, label, col in ((S0, 'S0 (HBY bit0)', DEEP_BLUE),
                               (S1, 'S1 (ZMH bit1)', ORANGE),
                               (S2, 'S2 (HBY bit2)', TEAL)):
        axA.plot(t, series, lw=1.7, color=col, label=label)
    for series, col in ((S0, DEEP_BLUE), (S1, ORANGE), (S2, TEAL)):
        axA.axhline(BAND_LOW, color=col, ls=':', lw=0.8, alpha=0.5)
        axA.axhline(BAND_HIGH, color=col, ls=':', lw=0.8, alpha=0.5)
    for key, col in (('S0', DEEP_BLUE), ('S1', ORANGE), ('S2', TEAL)):
        times = [c['time_h'] for c in row['verdict']['crossings'][key]
                 if lo <= c['time_h'] <= hi]
        axA.plot(times, np.full(len(times), -0.045), '|', ms=7, color=col)
    axA.set_ylim(-0.08, 1.12)
    axA.set_ylabel('DNA state\n(LR fraction)', fontsize=9.5)
    axA.legend(ncol=3, loc='upper left', fontsize=8.5, framealpha=0.9)
    axA.set_title('A  One read window per clock period: the arithmetic that the '
                  'counter is supposed to perform', fontsize=10, color=DEEP_BLUE, pad=6)
    style(axA)

    for series, label, col in ((y[idx['b0_I']], 'Int0 (HBY bit0)', DEEP_BLUE),
                               (y[idx['I1_zmh']], 'Int1 (ZMH bit1)', ORANGE),
                               (y[idx['b2_I']], 'Int2 (HBY bit2)', TEAL)):
        axB.plot(t, series, lw=1.7, color=col, label=label)
    axB.set_ylabel('Integrase pools\n(a.u.)', fontsize=9.5)
    axB.legend(ncol=3, loc='upper left', fontsize=8.5, framealpha=0.9)
    axB.set_title('B  Every carrier pulse is one narrow event, not a level',
                  fontsize=10, color=DEEP_BLUE, pad=6)
    style(axB)

    axC.plot(t, S['zw_g0'], lw=1.8, color=DEEP_BLUE, label='g0 (bit0 -> bit1 gate)')
    axC.plot(t, S['zw_g1'], lw=1.8, color=TEAL, label='g1 (bit1 -> bit2 gate)')
    axC.plot(t, S['zw_clock'], lw=1.4, color=BLUE_GRAY, label='clock gate')
    axC.axhline(GATE_ON, color=RISK_ORANGE, ls='--', lw=1.2)
    axC.text(hi, GATE_ON, f' gate threshold {GATE_ON:g}', ha='right', va='bottom',
             fontsize=8.5, color=RISK_ORANGE)
    for a in zoom_events(js, 'reverse_start_h', lo, hi):
        axC.axvline(a, color=RISK_ORANGE, ls='-', lw=1.0, alpha=0.55)
    for a in zoom_events(js, 'gate_start_h', lo, hi):
        axC.axvline(a, color=TEAL, ls='-', lw=1.0, alpha=0.55)
    for a in zoom_events(js, 'flip_h', lo, hi):
        axC.axvline(a, color=DEEP_BLUE, ls='--', lw=1.0, alpha=0.55)
    axC.plot([], [], color=RISK_ORANGE, lw=1.4, label='reset event')
    axC.plot([], [], color=TEAL, lw=1.4, label='gate open')
    axC.plot([], [], color=DEEP_BLUE, ls='--', lw=1.4, label='bit-2 flip')
    axC.set_ylabel('Carry gates', fontsize=9.5)
    axC.legend(ncol=3, loc='upper left', fontsize=8.5, framealpha=0.9)
    axC.set_title('C  Each reset is followed by exactly one gate opening and one '
                  'bit-2 flip', fontsize=10, color=DEEP_BLUE, pad=6)
    style(axC)

    axD.plot(t, S['zw_Jrev0'], lw=1.6, color=DEEP_BLUE, label='J_rev0 (bit0 reset flux)')
    axD.plot(t, S['zw_Jrev1'], lw=1.6, color=ORANGE, label='J_rev1 (bit1 reset flux)')
    axD.plot(t, S['zw_Jfwd2'], lw=1.3, color=TEAL, label='J_fwd2 (bit2 write)')
    axD.plot(t, S['zw_Jrev2'], lw=1.3, color=BLUE_GRAY, label='J_rev2 (bit2 reset)')
    axD.axhline(JREV_THRESHOLD, color=RISK_ORANGE, ls='--', lw=1.2)
    axD.text(hi, JREV_THRESHOLD, f' J_rev threshold {JREV_THRESHOLD:g}', ha='right',
             va='bottom', fontsize=8.5, color=RISK_ORANGE)
    axD.set_ylabel('Recombination\nflux (per h)', fontsize=9.5)
    axD.legend(ncol=2, loc='upper left', fontsize=8.5, framealpha=0.9)
    axD.set_title('D  The fluxes are the events: they return to zero between carries',
                  fontsize=10, color=DEEP_BLUE, pad=6)
    style(axD)
    axD.set_xlabel('Time (h)', fontsize=10)

    for ax in (axA, axB, axC, axD):
        ax.set_xlim(lo, hi)
    footer(axE, js, 'HBY bit0 + ZMH current-code A0/F0 and bit1 + HBY A1/F1 and bit2, '
                    '34 states')
    return fig


def zoom_events(js, field, lo, hi):
    row = js['horizons'][-1]
    out = []
    for chain in row['verdict']['events'].values():
        for a in chain['associations']:
            v = a.get(field)
            if v is not None and lo <= v <= hi:
                out.append(float(v))
    return sorted(out)


def s0_crossings(js, key):
    return js['horizons'][-1]['verdict']['crossings'][key]


# ------------------------------------------------------------------ checks
def collect_texts(fig):
    out = []
    for ax in fig.axes:
        out.extend(ax.texts)
        if ax.get_legend() is not None:
            out.extend(ax.get_legend().get_texts())
    return out


def verify_layout(fig):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    tick_ids = set()
    for ax in fig.axes:
        for t in list(ax.get_xticklabels()) + list(ax.get_yticklabels()):
            tick_ids.add(id(t))
    placed = [(t, t.get_window_extent(renderer=renderer)) for t in collect_texts(fig)
              if t.get_visible() and id(t) not in tick_ids]
    clashes = []
    for i in range(len(placed)):
        for j in range(i + 1, len(placed)):
            inter = Bbox.intersection(placed[i][1], placed[j][1])
            if inter is not None and inter.width > 1.0 and inter.height > 1.0:
                clashes.append((placed[i][0].get_text()[:24],
                                placed[j][0].get_text()[:24]))
    check(not clashes, f'no two figure labels overlap (found {len(clashes)}: {clashes[:3]})')
    figbb = fig.get_window_extent(renderer=renderer)
    outside = [t.get_text()[:28] for t, bb in placed
               if bb.x0 < figbb.x0 - 1 or bb.x1 > figbb.x1 + 1
               or bb.y0 < figbb.y0 - 1 or bb.y1 > figbb.y1 + 1]
    check(not outside, f'every figure label is inside the canvas (outside: {outside})')


def save(fig, stem, labels):
    FIG.mkdir(parents=True, exist_ok=True)
    outs = []
    for ext in ('png', 'svg', 'pdf'):
        p = FIG / f'{stem}.{ext}'
        fig.savefig(p, dpi=DPI, bbox_inches='tight', facecolor=fig.get_facecolor())
        outs.append(p)
    check(all(p.exists() and p.stat().st_size > 8000 for p in outs),
          f'{stem}: written as PNG + SVG + PDF with non-trivial size')
    svg = [p for p in outs if p.suffix == '.svg'][0].read_text(encoding='utf-8')
    check('<text' in svg, f'{stem}: SVG uses <text> (svg.fonttype = none)')
    for lab in labels:
        esc = lab.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')
        check(esc in svg, f'{stem}: SVG carries {lab!r}')
    for prefix in FOOTER_PREFIXES:
        check(prefix in svg, f'{stem}: SVG carries the footer element {prefix!r}')
    return outs


def main():
    plt.rcParams['svg.fonttype'] = SVG_FONTTYPE
    plt.rcParams['font.family'] = 'DejaVu Sans'
    js, S = load()
    row = js['horizons'][-1]
    check(js['meta']['hours'] >= 600.0,
          f"run horizon is {js['meta']['hours']:.0f} h (the point of this probe)")
    check(row['certified_v1'], 'the 1200 h horizon is certified_v1')
    check(all(h['certified_v1'] for h in js['horizons']),
          'certified_v1 holds at 300 / 600 / 1200 h')
    check(row['steady']['reads'] > 100,
          f"{row['steady']['reads']} steady read windows (vs 19 at 300 h, 48 in the "
          '51-state certification)')
    check(S['zw_time'][0] >= 0.35 * js['meta']['hours'] and
          S['zw_time'][-1] <= 0.65 * js['meta']['hours'],
          'the zoom window sits in the middle of the run')
    inside = ((S['read_trough'] >= S['zoom'][0]) & (S['read_trough'] <= S['zoom'][1])).sum()
    check(inside >= 3, f'the zoom window contains {int(inside)} complete read windows')

    print('\n-- overview figure')
    fig = overview(js, S)
    verify_layout(fig)
    save(fig, 'probe_C_long_horizon_overview', OVERVIEW_LABELS)
    plt.close(fig)

    print('\n-- mid-run zoom figure')
    fig = zoom(js, S)
    verify_layout(fig)
    save(fig, 'probe_C_midwindow_zoom', ZOOM_LABELS)
    plt.close(fig)

    bad = [m for ok, m in CHECKS if not ok]
    print(f'\n{len(CHECKS) - len(bad)}/{len(CHECKS)} checks passed')
    if bad:
        raise SystemExit('FAILED CHECKS:\n' + '\n'.join(bad))


if __name__ == '__main__':
    main()
