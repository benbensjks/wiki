"""Flip-behaviour diagnostics at the frozen three-bit working point.

Parameter set used (and why)
    threebit51_selected_v1.json -> n_A1_gate = 6, carry1 mRNA half-life 2 min,
    A1/F1 maturation 32.5 min, F1 production exponent 4, uM_per_au 5.75.
    This is the only parameter set with a complete evidence chain: 14 bit1
    reverse events = 14 g1 gate windows = 14 bit2 flips, 56 read windows strictly
    mod-8 (48 after the drop-8 cut), and a perfect R/F carry alternation with
    zero same-type adjacency. n=7 and n=8 have slightly lower far-off leakage but
    the improvement is saturated (+1.22x / +1.06x) at a higher cooperativity, so
    n=6 is the frozen choice.

What is plotted
    fig1_counting_600h      the full run: S0/S1/S2, the decoded 3-bit value, the
                            carry gate, and the bit2 pools
    fig2_R_vs_F_carry       one R-type and one F-type carry side by side, i.e.
                            the two flip directions and their causal chains
    fig3_read_windows       every read window against the 0.30/0.70 bands, with
                            commitment and the value sequence
    fig4_carry_alternation  the R/F super-period structure and the far-off leak

Output goes to threebit51_results/flip_diagnostics/ together with the trajectory
CSV, the read-window table and a working-point record.

    python plot_flip_behaviour.py [--hours 600]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]          # final_reconstruction/
OUT = Path(__file__).resolve().parent               # threebit51_results/flip_diagnostics/
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'plausibility'))

from plausibility_common import frozen_carry0, frozen_extension, restore  # noqa: E402
from plausibility_common import build_threebit, gate_windows, sha256  # noqa: E402
from scan_carry_pairing import segment_rows  # noqa: E402
from working_point import working_point_block  # noqa: E402
from verify_threebit51 import analyse_threebit, read_windows_three  # noqa: E402

BAND_LOW, BAND_HIGH = 0.30, 0.70
C_F = '#1F78B4'      # F-type carry
C_R = '#E31A1C'      # R-type carry
C_DNA = ('#6A3D9A', '#33A02C', '#FF7F00')


def simulate(hours, sample_min, max_step_min):
    model, patch = build_threebit(carry1=frozen_carry0())
    try:
        sol = model.simulate(hours=hours, sample_min=sample_min,
                             max_step_min=max_step_min)
        sig = model.diagnostic_signals(sol.y)
        analysis = analyse_threebit(model, sol, hours)
        wp = working_point_block(model)
    finally:
        restore(patch)
    return model, sol, sig, analysis, wp


def table(sol, sig, analysis, flux):
    t, y = sol.t, sol.y
    reads = read_windows_three(t, flux, (y[16], y[27], y[44]))
    rows = []
    for r in reads:
        rows.append(dict(cycle=r['cycle'], trough_h=r['trough_h'], value=r['value'],
                         S0=r['bits'][0]['label'], S1=r['bits'][1]['label'],
                         S2=r['bits'][2]['label'],
                         S0_val=float(y[16][np.argmin(np.abs(t - r['trough_h']))]),
                         S1_val=float(y[27][np.argmin(np.abs(t - r['trough_h']))]),
                         S2_val=float(y[44][np.argmin(np.abs(t - r['trough_h']))]),
                         commitment_min=float(min(b['commitment'] for b in r['bits'])),
                         window_start_h=r['window_start_h'],
                         window_end_h=r['window_end_h']))
    return pd.DataFrame(rows)


model_flux_cache = [None]


def fig1(t, y, sig, reads, wp):
    fig, ax = plt.subplots(3, 1, figsize=(13, 9), sharex=True, constrained_layout=True)
    for i, (idx, name, col) in enumerate(zip((16, 27, 44), ('bit0  S0', 'bit1  S1', 'bit2  S2'),
                                             C_DNA)):
        ax[0].plot(t, y[idx], color=col, lw=1.4, label=name)
    ax[0].axhspan(-0.05, BAND_LOW, color='0.85', zorder=0)
    ax[0].axhspan(BAND_HIGH, 1.05, color='0.85', zorder=0)
    ax[0].axhline(BAND_LOW, color='k', ls=':', lw=0.9)
    ax[0].axhline(BAND_HIGH, color='k', ls=':', lw=0.9)
    ax[0].set_ylim(-0.05, 1.05)
    ax[0].set_ylabel('DNA state')
    ax[0].set_title('Three-bit counter at the frozen working point '
                    '(n_A1_gate=%g, carry1 mRNA %g min, maturation %g min)\n'
                    'grey bands are the non-read zones below 0.30 and above 0.70'
                    % (wp['n_A1_gate_effective'], wp['carry1_mrna_half_life_min'],
                       wp['A1_maturation_half_life_min']), fontsize=10)
    ax[0].legend(loc='center right', fontsize=9, ncol=3)

    rt = reads.trough_h.to_numpy()
    rv = reads.value.to_numpy(dtype=float)
    rv = np.where(np.isnan(rv), np.nan, rv)
    ax[1].step(rt, rv, where='mid', color='k', lw=1.6)
    ax[1].scatter(rt, rv, s=14, color='k', zorder=3)
    ax[1].set_yticks(range(8))
    ax[1].set_ylim(-0.6, 7.6)
    ax[1].set_ylabel('read value (mod 8)')
    ax[1].set_title('Finite-window readout: %d windows, all committed, strictly mod-8'
                    % len(reads), fontsize=10)
    valid = reads.commitment_min.to_numpy()
    ax[1].text(0.01, 0.06, 'min commitment = %.2f   unlabelled = %d'
               % (valid.min(), int((reads.value.isna()).sum())),
               transform=ax[1].transAxes, fontsize=9)

    g1 = np.asarray(sig['g1'])
    for k, w in enumerate(gate_windows(t, g1)):
        s2 = float(y[44][w['i0']])
        col = C_F if s2 < 0.5 else C_R
        ax[2].axvspan(w['start_h'], w['end_h'], color=col, alpha=0.45, lw=0)
    ax[2].plot(t, g1 / max(g1.max(), 1e-12), color='k', lw=1.2, label='g1 (scaled)')
    for idx, name, col in ((36, 'Int2', '#0072B2'), (42, 'RDF2', '#D55E00'),
                           (43, 'C2', '#009E73')):
        v = y[idx]
        ax[2].plot(t, v / max(v.max(), 1e-12), color=col, lw=1.1, label=name)
    ax[2].set_ylabel('normalised')
    ax[2].set_xlabel('time (h)')
    ax[2].legend(loc='upper right', fontsize=8, ncol=4)
    ax[2].set_title('carry gate windows: blue = F-type (forward write), '
                    'red = R-type (reverse write)', fontsize=10)
    fig.savefig(OUT / 'fig1_counting_600h.png', dpi=190, bbox_inches='tight',
                facecolor='white')
    plt.close(fig)


def _carry_zoom(ax, t, y, sig, w, kind, reads):
    lo = w['start_h'] - 5.0
    hi = w['end_h'] + 10.0
    m = (t >= lo) & (t <= hi)
    col = C_F if kind == 'F' else C_R
    ax[0].axvspan(w['start_h'], w['end_h'], color=col, alpha=0.18, lw=0)
    ax[0].plot(t[m], y[16][m], color=C_DNA[0], lw=1.3, label='S0')
    ax[0].plot(t[m], y[27][m], color=C_DNA[1], lw=1.3, label='S1')
    ax[0].plot(t[m], y[44][m], color=C_DNA[2], lw=1.8, label='S2')
    ax[0].axhline(BAND_LOW, color='k', ls=':', lw=0.8)
    ax[0].axhline(BAND_HIGH, color='k', ls=':', lw=0.8)
    ax[0].set_ylim(-0.06, 1.06)
    ax[0].set_ylabel('DNA state')
    ax[0].set_title('%s-type carry  (S2 at gate start = %.3f -> %s write)'
                    % (kind, w['S2_start'], 'FORWARD' if kind == 'F' else 'REVERSE'),
                    fontsize=10)
    ax[0].legend(loc='center left', fontsize=8, ncol=3)

    ax[1].axvspan(w['start_h'], w['end_h'], color=col, alpha=0.18, lw=0)
    ax[1].plot(t[m], sig['g1'][m], color='k', lw=1.3, label='g1 (gate)')
    j1 = np.asarray(sig['J_rev1'])
    ax[1].plot(t[m], j1[m] / max(j1.max(), 1e-12), color='#8C564B', lw=1.1,
               label='J_rev1 (bit1 reverse, scaled)')
    ax[1].set_ylabel('gate / flux')
    ax[1].legend(loc='upper right', fontsize=8)
    ax[1].set_title('the carry is triggered by the bit1 reverse pulse through the gate',
                    fontsize=9)

    ax[2].axvspan(w['start_h'], w['end_h'], color=col, alpha=0.18, lw=0)
    for idx, name, c in ((36, 'Int2', '#0072B2'), (42, 'RDF2', '#D55E00'),
                         (43, 'C2', '#009E73')):
        ax[2].plot(t[m], y[idx][m], color=c, lw=1.2, label=name)
    ax[2].set_ylabel('bit2 pools (a.u.)')
    ax[2].set_xlabel('time (h)')
    ax[2].legend(loc='upper right', fontsize=8)
    ax[2].set_title('Int2 / RDF2 / complex drive the flip', fontsize=9)

    for a in ax:
        for _, r in reads.iterrows():
            if lo <= r.window_start_h <= hi:
                a.axvspan(r.window_start_h, r.window_end_h, color='0.55', alpha=0.18, lw=0)
    yv = reads[(reads.window_start_h >= lo) & (reads.window_end_h <= hi)]
    for _, r in yv.iterrows():
        ax[0].text(r.trough_h, 1.045, str(int(r.value)), ha='center', fontsize=8)
    ax[0].set_xlim(lo, hi)


def fig2(t, y, sig, rows, reads):
    Fs = [r for r in rows if r['type'] == 'F']
    Rs = [r for r in rows if r['type'] == 'R']
    if not Fs or not Rs:
        print('fig2 skipped: need one F-type and one R-type carry, found %d/%d '
              '(run more hours)' % (len(Fs), len(Rs)))
        return
    F, R = Fs[-1], Rs[-1]
    fig, axes = plt.subplots(3, 2, figsize=(14, 9), constrained_layout=True)
    _carry_zoom(axes[:, 0], t, y, sig, F, 'F', reads)
    _carry_zoom(axes[:, 1], t, y, sig, R, 'R', reads)
    fig.suptitle('The two flip directions at the frozen working point  '
                 '(grey bands = read windows, numbers = the read value)',
                 fontsize=11)
    fig.savefig(OUT / 'fig2_R_vs_F_carry.png', dpi=190, bbox_inches='tight',
                facecolor='white')
    plt.close(fig)


def fig3(reads):
    fig, ax = plt.subplots(2, 1, figsize=(13, 6.5), sharex=True, constrained_layout=True)
    for _, r in reads.iterrows():
        for i, (col, key) in enumerate(zip(C_DNA, ('S0_val', 'S1_val', 'S2_val'))):
            ax[0].plot([r.window_start_h, r.window_end_h], [r[key]] * 2,
                       color=col, lw=2.6, solid_capstyle='butt')
    ax[0].axhspan(-0.05, BAND_LOW, color='0.88', zorder=0)
    ax[0].axhspan(BAND_HIGH, 1.05, color='0.88', zorder=0)
    ax[0].axhline(BAND_LOW, color='k', ls=':', lw=0.9)
    ax[0].axhline(BAND_HIGH, color='k', ls=':', lw=0.9)
    ax[0].set_ylim(-0.05, 1.05)
    ax[0].set_ylabel('state inside each\nread window')
    ax[0].set_title('Every read window sits inside a committed band '
                    '(purple bit0, green bit1, orange bit2)', fontsize=10)
    ax[1].plot(reads.trough_h, reads.commitment_min, 'o-', color='k', ms=4, lw=1)
    ax[1].axhline(0.80, color=C_R, ls='--', lw=1.2)
    ax[1].text(reads.trough_h.iloc[1], 0.83, 'acceptance floor 0.80', color=C_R, fontsize=9)
    ax[1].set_ylim(0.75, 1.02)
    ax[1].set_ylabel('worst-bit\ncommitment')
    ax[1].set_xlabel('time (h)')
    ax[1].set_title('commitment of the least-committed bit in each window '
                    '(1.00 = every sample in band)', fontsize=10)
    fig.savefig(OUT / 'fig3_read_windows.png', dpi=190, bbox_inches='tight',
                facecolor='white')
    plt.close(fig)


def fig4(t, y, rows, reads):
    """Two panels, because one cannot serve both scales.

    Panel 1 zooms on FOUR consecutive carries: at the full-run scale the gate
    window is 1.6 h out of a 42.9 h period, so the coloured bars degenerate into
    slivers.  Panel 2 then shows the carry-type sequence for the whole run.
    """
    ev = list(rows)
    k0 = max(0, len(ev) // 2 - 2)
    zoom = ev[k0:k0 + 4]
    lo = zoom[0]['start_h'] - 2.0
    hi = (zoom[-1]['i1'] + 2) if len(zoom) < len(ev) else len(t) - 1
    hi_h = float(t[min(hi, len(t) - 1)]) + 3.0

    fig, ax = plt.subplots(2, 1, figsize=(13, 7.5),
                           gridspec_kw=dict(height_ratios=[2.2, 1]), constrained_layout=True)

    for i, r in enumerate(zoom):
        col = C_F if r['type'] == 'F' else C_R
        ax[0].barh(i, r['width_h'], left=r['start_h'], color=col, height=0.55)
        j0 = r['i1']
        ax[0].barh(i, r['tail_h'], left=float(t[j0]), color='0.45', height=0.55)
        nxt = float(t[ev[k0 + i + 1]['i0']]) if (k0 + i + 1) < len(ev) else hi_h
        ax[0].barh(i, max(nxt - (float(t[j0]) + r['tail_h']), 0.0),
                   left=float(t[j0]) + r['tail_h'], color='0.86', height=0.55)
        ax[0].text(r['start_h'] + r['width_h'] / 2, i + 0.42,
                   'gate %.2f h' % r['width_h'], ha='center', fontsize=7.5, color=col)
        ax[0].text(float(t[j0]) + r['tail_h'] / 2, i - 0.42, 'tail %.2f h' % r['tail_h'],
                   ha='center', fontsize=7.5, color='0.3')
        ax[0].text((float(t[j0]) + r['tail_h'] + nxt) / 2, i,
                   'far-off %.1f h' % r['far_off_h'], ha='center', fontsize=8, color='0.45')
    ax[0].set_xlim(lo, hi_h)
    ax[0].set_ylim(-0.7, len(zoom) - 0.3)
    ax[0].set_yticks(range(len(zoom)))
    ax[0].set_yticklabels(['carry %d\n(%s)' % (k0 + i, r['type'])
                           for i, r in enumerate(zoom)], fontsize=8)
    ax[0].set_xlabel('time (h)')
    ax[0].set_title('One carry period in three segments: gate window (coloured) + '
                    'Int2 tail (dark grey) + far-off dwell (light grey)\n'
                    'blue = F-type (forward write), red = R-type (reverse write); '
                    'the types alternate strictly across the whole run', fontsize=10)

    types = ''.join(r['type'] for r in ev)
    ax[1].step(range(len(types)), [1 if x == 'F' else 0 for x in types], where='mid',
               color='k', lw=1.2)
    ax[1].scatter(range(len(types)), [1 if x == 'F' else 0 for x in types],
                  c=[C_F if x == 'F' else C_R for x in types], s=45, zorder=3)
    ax[1].set_yticks([0, 1])
    ax[1].set_yticklabels(['R-type\n(reverse)', 'F-type\n(forward)'])
    ax[1].set_ylim(-0.6, 1.6)
    ax[1].set_xlim(-0.6, len(types) - 0.4)
    ax[1].set_xlabel('carry index (0 .. %d over %g h)' % (len(types) - 1, t[-1]))
    ax[1].set_title('carry type sequence  %s   ->  0 same-type adjacencies'
                    % types, fontsize=10)
    fig.savefig(OUT / 'fig4_carry_alternation.png', dpi=190, bbox_inches='tight',
                facecolor='white')
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--hours', type=float, default=600.0)
    ap.add_argument('--sample-min', type=float, default=2.0)
    ap.add_argument('--max-step-min', type=float, default=2.0)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    model, sol, sig, analysis, wp = simulate(args.hours, args.sample_min,
                                             args.max_step_min)
    t, y = sol.t, sol.y
    from plausibility_common import flux_array
    flux = flux_array(model, sol)
    reads = table(sol, sig, analysis, flux)

    g1 = np.asarray(sig['g1'])
    wins = gate_windows(t, g1)
    rows = segment_rows(t, sig['J_rev2'], sig['J_fwd2'], sig['S2'], y[36], wins, 0.05)

    traj = pd.DataFrame(dict(
        time_h=t, S0=y[16], S1=y[27], S2=y[44],
        I0=y[8], R0=y[14], C0=y[15], I1=y[19], R1=y[25], C1=y[26],
        I2=y[36], R2=y[42], C2=y[43], A1=y[45], F1=y[46],
        g1=g1, J_rev1=sig['J_rev1'], J_rev2=sig['J_rev2'], J_fwd2=sig['J_fwd2']))
    traj.to_csv(OUT / 'flip_trajectory.csv', index=False, encoding='utf-8')
    reads.to_csv(OUT / 'flip_read_windows.csv', index=False, encoding='utf-8')

    fig1(t, y, sig, reads, wp)
    fig2(t, y, sig, rows, reads)
    fig3(reads)
    fig4(t, y, rows, reads)

    types = ''.join(r['type'] for r in rows)
    summary = dict(
        working_point=wp,
        hours=args.hours, sample_min=args.sample_min,
        read_windows=len(reads),
        read_windows_unlabelled=int(reads.value.isna().sum()),
        min_commitment=float(reads.commitment_min.min()),
        steady_sequence=analysis['steady_state']['sequence'],
        certified=bool(analysis['certified']),
        crossings=len(analysis['bit2_crossings']),
        gate_windows=len(analysis['carry1_gate_events']),
        reverse_events=len(analysis['bit1_reverse_events']),
        carry_type_sequence=types,
        same_type_adjacency=sum(1 for a, b in zip(types, types[1:]) if a == b),
        gate_on_h_median=float(np.median([r['width_h'] for r in rows])),
        tail_h_median=float(np.median([r['tail_h'] for r in rows])),
        far_off_h_median=float(np.median([r['far_off_h'] for r in rows])),
        figures=['fig1_counting_600h.png', 'fig2_R_vs_F_carry.png',
                 'fig3_read_windows.png', 'fig4_carry_alternation.png'],
        model_sha256=wp['model_sha256'],
        note=('the parameter set is the frozen working point from '
              'plausibility/threebit51_selected_v1.json; nothing was tuned here'),
    )
    (OUT / 'flip_diagnostics_summary.json').write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')

    print(json.dumps({k: v for k, v in summary.items() if k != 'working_point'},
                     ensure_ascii=False, indent=2))
    print()
    print('working point: n_A1_gate=%g  mRNA=%g  mat=%g  uM_per_au=%g'
          % (wp['n_A1_gate_effective'], wp['carry1_mrna_half_life_min'],
             wp['A1_maturation_half_life_min'], wp['uM_per_au']))
    print('wrote figures + tables to', OUT)


if __name__ == '__main__':
    main()
