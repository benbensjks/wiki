"""Figures for the HBY-ZMH-ZMH failure attribution.

Style follows the wiki rule: English labels, large fonts, thick lines.

  fig1_before_after        the failure and the single-parameter fix, same window
  fig2_clock_K_window      the operating window along clock_K, with the two
                           distinct failure modes marked
  fig3_mechanism           cross-experiment discriminator: reverse drive product
                           over its Hill threshold, vs the reverse Hill value
  fig4_eight_phases        eight orbit-phase initial states all lock to mod 8

Reads only saved artifacts (round 2 and round 3 trajectories and the phase CSV);
writes figures into failure_attribution/figures/. No existing file is modified.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
import numpy as np

sys.dont_write_bytecode = True

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT.parent))
sys.path.insert(0, str(ROOT.parent.parent))

import diagnose_round1 as D                                       # noqa: E402
import verify_hzh as V                                            # noqa: E402
from diagnose_round2 import band_margins                          # noqa: E402

FIG = ROOT / 'figures'
RES = ROOT / 'results'
LOW, HIGH = V.LOW, V.HIGH

CLR = dict(s0='#2685ad', s1='#d97726', s2='#43875a', ok='#2e7d32',
           bad='#c62828', mid='#ef9a00', gate='#d97726', clock='#6a1b9a',
           int0='#7550a1', int2='#2685ad', rdf2='#b8860b',
           fwd='#2685ad', rev='#c62828')


def newest(pat):
    d = sorted(RES.glob(pat), key=lambda p: p.name)
    if not d:
        raise SystemExit(f'no results for {pat}')
    return d[-1]


def style(ax, xlabel=None, ylabel=None, title=None, legend=False):
    ax.grid(alpha=.22)
    ax.tick_params(labelsize=13)
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=15)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=15)
    if title:
        ax.set_title(title, fontsize=15)
    if legend:
        ax.legend(fontsize=12, loc='best', framealpha=.9)


def load_traj(path):
    z = np.load(path)
    return z['time_h'], z['states']


def main():
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8')
    FIG.mkdir(exist_ok=True)
    r2 = newest('round2_*')
    r3 = newest('round3_*')
    print('round2 =', r2.name, '| round3 =', r3.name)

    # signals need a model instance; the Han upstream is ~5 s
    han = D.HanInput(D.HOURS)
    prefix = D.HZHModel('han', han)

    def sig_of(states, autoreg, clock_K):
        m = D.HZZModel(prefix, autoreg)
        m.clock_K = clock_K
        return m.signals(states), m

    # =====================================================================
    # fig 1: before / after
    # =====================================================================
    tB, yB = load_traj(r2 / 'traj_K0.40_off.npz')
    tA, yA = load_traj(r3 / 'traj_work300.npz')
    sB, _ = sig_of(yB, False, 0.40)
    sA, _ = sig_of(yA, False, 0.10)

    win = (150.0, 215.0)
    mB = (tB >= win[0]) & (tB <= win[1])
    mA = (tA >= win[0]) & (tA <= win[1])

    panels = [
        ('Int0 and clock gate', 'Int0 (a.u.)', 'clock gate'),
        ('DNA LR fraction', 'LR fraction', None),
        ('Decoded value', 'value', None),
        ('Gate g1 and bit2 flux', 'g1', 'J (/h)'),
        ('Bit2 species', 'Int2 (a.u.)', 'RDF2 (a.u.)'),
    ]
    fig, axes = plt.subplots(5, 2, figsize=(15, 13.5), sharex=True,
                             layout='constrained')
    for col, (tt, yy, ss, tag, k) in enumerate([
            (tB, yB, sB, 'BEFORE  clock_K = 0.40  (as published)', 0.40),
            (tA, yA, sA, 'AFTER  clock_K = 0.10  (single-lever fix)', 0.10)]):
        m = mB if col == 0 else mA
        xt = tt[m]
        axes[0, col].set_title(tag, fontsize=15,
                               color=CLR['bad'] if col == 0 else CLR['ok'])
        ax = axes[0, col]
        ax.plot(xt, ss['int0'][m], color=CLR['int0'], lw=2.4, label='Int0')
        ax.set_ylabel(panels[0][1], fontsize=15)
        a2 = ax.twinx()
        a2.plot(xt, ss['clock'][m], color=CLR['clock'], lw=2.0, ls='--',
                label='clock')
        a2.set_ylabel(panels[0][2], fontsize=15)
        a2.tick_params(labelsize=13)
        ln = ax.get_lines() + a2.get_lines()
        ax.legend(ln, [l.get_label() for l in ln], fontsize=12, loc='upper right')
        ax.grid(alpha=.22)
        ax.tick_params(labelsize=13)

        ax = axes[1, col]
        for b, c in (('S0', CLR['s0']), ('S1', CLR['s1']), ('S2', CLR['s2'])):
            ax.plot(xt, ss[b][m], color=c, lw=2.2, label=b)
        ax.axhspan(0, LOW, color=CLR['s2'], alpha=.09)
        ax.axhspan(HIGH, 1, color=CLR['s0'], alpha=.09)
        ax.axhline(LOW, ls='--', color='gray', lw=1.2)
        ax.axhline(HIGH, ls='--', color='gray', lw=1.2)
        ax.set_ylim(-.05, 1.05)
        style(ax, ylabel=panels[1][1], legend=True)

        ax = axes[2, col]
        rw = V.read_windows(tt, ss)[1]
        rt = np.array([r['trough_h'] for r in rw])
        vv = np.array([np.nan if r['value'] is None else r['value'] for r in rw])
        k2 = (rt >= win[0]) & (rt <= win[1])
        ax.step(rt[k2], vv[k2], where='mid', color='#33383a', lw=2.4)
        ax.scatter(rt[k2], vv[k2], c=vv[k2], cmap='viridis', vmin=0, vmax=7,
                   s=60, zorder=5, edgecolor='k', linewidth=.6)
        bad = np.isnan(vv[k2])
        if bad.any():
            ax.scatter(rt[k2][bad], np.full(bad.sum(), -0.9), marker='x',
                       color=CLR['bad'], s=90, linewidth=3.0,
                       label='unlabelled window')
            ax.legend(fontsize=12, loc='lower right')
        ax.set_yticks(range(8))
        ax.set_ylim(-1.4, 7.6)
        style(ax, ylabel=panels[2][1])

        ax = axes[3, col]
        ax.plot(xt, ss['g1'][m], color=CLR['gate'], lw=2.4, label='g1')
        ax.axhline(V.GATE_THRESHOLD, color='gray', ls=':', lw=1.6)
        ax.text(win[0] + 1, V.GATE_THRESHOLD + .004, 'event threshold',
                fontsize=11, color='gray')
        ax.set_ylabel(panels[3][1], fontsize=15)
        a3 = ax.twinx()
        a3.plot(xt, ss['J_fwd2'][m], color=CLR['fwd'], lw=1.8, alpha=.85,
                label='forward')
        a3.plot(xt, ss['J_rev2'][m], color=CLR['rev'], lw=1.8, alpha=.85,
                label='reverse')
        a3.set_ylabel(panels[3][2], fontsize=15)
        a3.tick_params(labelsize=13)
        ln = ax.get_lines()[:1] + a3.get_lines()
        ax.legend(ln, [l.get_label() for l in ln], fontsize=12, loc='upper right')
        ax.grid(alpha=.22)
        ax.tick_params(labelsize=13)

        ax = axes[4, col]
        ax.plot(xt, yy[20][m], color=CLR['int2'], lw=2.4, label='mature Int2')
        ax.set_ylabel(panels[4][1], fontsize=15)
        a4 = ax.twinx()
        a4.plot(xt, yy[22][m], color=CLR['rdf2'], lw=2.4, label='RDF2')
        a4.set_ylabel(panels[4][2], fontsize=15)
        a4.tick_params(labelsize=13)
        ln = ax.get_lines() + a4.get_lines()
        ax.legend(ln, [l.get_label() for l in ln], fontsize=12, loc='upper right')
        ax.grid(alpha=.22)
        ax.tick_params(labelsize=13)

        seq = V.read_verdict(rw, V.DROP)['sequence']
        axes[4, col].text(.01, .05, 'steady: ' + seq, transform=axes[4, col].transAxes,
                          fontsize=13, color=CLR['bad'] if col == 0 else CLR['ok'],
                          fontweight='bold')
        for a in axes[:, col]:
            a.set_xlim(*win)

    axes[-1, 0].set_xlabel('Time (h)', fontsize=15)
    axes[-1, 1].set_xlabel('Time (h)', fontsize=15)
    fig.suptitle('HBY-ZMH-ZMH three-bit cascade: the failure is one clock '
                 'threshold, not the architecture\n'
                 'same Han upstream, same HBY bit0, same ZMH middle, same '
                 'exploratory ZMH tail, A1 autoregulation OFF', fontsize=16)
    for ext in ('png', 'svg'):
        fig.savefig(FIG / f'fig1_before_after.{ext}', dpi=150)
    plt.close(fig)
    print('fig1 done')

    # =====================================================================
    # fig 2: the clock_K operating window
    # =====================================================================
    with (r2 / 'round2_all.csv').open(encoding='utf-8') as fh:
        rows = list(csv.DictReader(fh))

    def g(r, k):
        try:
            return float(r[k])
        except (TypeError, ValueError):
            return np.nan

    fig, axes = plt.subplots(3, 1, figsize=(11.5, 12), sharex=True,
                             layout='constrained')
    for ax, key, lab, colr in (
            (axes[0], 'gate_dose_max', 'worst gate-pulse dose (h)', CLR['gate']),
            (axes[1], 'RDF2_peak', 'RDF2 peak (a.u.)', CLR['rdf2']),
            (axes[2], 'S2_band_margin', 'S2 band margin', CLR['s2'])):
        for arm, mk, ls in (('off', 'o', '-'), ('on', 's', '--')):
            pts = sorted([r for r in rows if r['arm'] == arm],
                         key=lambda r: g(r, 'clock_K'))
            x = [g(r, 'clock_K') for r in pts]
            y = [g(r, key) for r in pts]
            ax.plot(x, y, marker=mk, ls=ls, lw=2.4, ms=9, color=colr,
                    alpha=.95 if arm == 'off' else .6,
                    label=f'autoregulation {arm}')
            for r, xx, yy in zip(pts, x, y):
                ok = r['counting'] == 'True' and r['event'] == 'True'
                ax.scatter([xx], [yy], s=210, facecolor='none',
                           edgecolor=CLR['ok'] if ok else CLR['bad'],
                           linewidth=2.6, zorder=6)
        style(ax, ylabel=lab, legend=True)
    axes[0].axhline(V.GATE_MIN_DOSE_H, color='gray', ls=':', lw=1.8)
    axes[0].text(0.021, V.GATE_MIN_DOSE_H + .003,
                 'event-rule minimum dose 0.02 h', fontsize=11, color='gray')
    axes[2].axhline(0.0, color='black', lw=1.6)
    axes[2].text(0.021, .01, 'band margin 0', fontsize=11)
    axes[2].set_xscale('log')
    axes[2].set_xticks([0.02, 0.05, 0.07, 0.1, 0.2, 0.25, 0.3, 0.35, 0.4])
    axes[2].set_xticklabels(['0.02', '0.05', '0.07', '0.10', '0.20', '0.25',
                             '0.30', '0.35', '0.40'], fontsize=13)
    axes[2].set_xlabel('clock_K  (log scale)', fontsize=15)
    axes[2].scatter([], [], s=210, facecolor='none', edgecolor=CLR['ok'],
                    linewidth=2.6, label='counting AND event pass')
    axes[2].scatter([], [], s=210, facecolor='none', edgecolor=CLR['bad'],
                    linewidth=2.6, label='fails')
    axes[2].legend(fontsize=12, loc='lower right')
    fig.suptitle('clock_K operating window: one-sided, and two different '
                 'failure modes\n'
                 'upper boundary in (0.20, 0.25) for OFF and (0.25, 0.30) for ON;'
                 ' no lower boundary down to 0.02', fontsize=15)
    for ext in ('png', 'svg'):
        fig.savefig(FIG / f'fig2_clock_K_window.{ext}', dpi=150)
    plt.close(fig)
    print('fig2 done')

    # =====================================================================
    # fig 3: cross-experiment discriminator
    # =====================================================================
    fig, ax = plt.subplots(figsize=(11.5, 8), layout='constrained')
    xs, ys, ok = [], [], []
    for r in rows:
        ir = g(r, 'i_times_r')
        kd = 0.8
        xs.append(ir / kd)
        ys.append(g(r, 'reverse_hill'))
        ok.append(r['counting'] == 'True' and r['event'] == 'True')
    xs, ys, ok = np.array(xs), np.array(ys), np.array(ok)
    ax.scatter(xs[ok], ys[ok], s=200, color=CLR['ok'], edgecolor='k',
               linewidth=1.2, zorder=5, label='HZZ tail: counts (this work)')
    ax.scatter(xs[~ok], ys[~ok], s=200, color=CLR['bad'], edgecolor='k',
               linewidth=1.2, zorder=5, marker='X',
               label='HZZ tail: fails (this work)')
    for x, y, o, r in zip(xs, ys, ok, rows):
        if o or g(r, 'clock_K') in (0.4, 0.25, 0.35):
            ax.annotate(f"K={g(r,'clock_K'):g} {r['arm']}", (x, y),
                        textcoords='offset points', xytext=(8, 7), fontsize=11)

    # ZMH bit0 driven by Han upstream with its normalised interface: the same
    # discriminator, from 韩亚轩上游替换_结果分析.md (2026-09-27).
    zmh = [(1.5732, 0.7122, 'ZMH bit0, Zeng upstream (worked)'),
           (0.4213, 0.1508, 'ZMH bit0, Han upstream (bit0 never returned)')]
    for x, y, lab in zmh:
        ax.scatter([x], [y], s=260, marker='D', facecolor='none',
                   edgecolor='#1f4e79', linewidth=2.8, zorder=6)
        ax.annotate(lab, (x, y), textcoords='offset points', xytext=(10, -18),
                    fontsize=11, color='#1f4e79')

    xx = np.logspace(-2, 2, 400)
    ax.plot(xx, xx ** 2 / (1 + xx ** 2), color='gray', lw=2.0, ls='--',
            label=r'$x^2/(1+x^2)$,  x = drive / threshold')
    ax.axhline(0.5, color='gray', ls=':', lw=1.6)
    ax.text(1.05, 0.52, 'Hill half-point', fontsize=11, color='gray')
    ax.axvline(1.0, color='black', lw=1.4, ls=':')
    ax.set_xscale('log')
    ax.set_xlabel('reverse drive product / its Hill threshold  '
                  r'($i\cdot r / K_{D,comp}$,  or  $Int\cdot RDF/3.2$)', fontsize=15)
    ax.set_ylabel('reverse Hill value (fraction of maximum)', fontsize=15)
    ax.set_ylim(-.03, 1.06)
    ax.set_title('One discriminator across two independent experiments\n'
                 'every passing point sits at drive/threshold >~ 6; every '
                 'failing point below ~5', fontsize=15)
    style(ax, legend=True)
    ax.legend(fontsize=12, loc='lower right')
    for ext in ('png', 'svg'):
        fig.savefig(FIG / f'fig3_mechanism.{ext}', dpi=150)
    plt.close(fig)
    print('fig3 done')

    # =====================================================================
    # fig 4: eight orbit-phase initial states
    # =====================================================================
    with (r3 / 'phases.csv').open(encoding='utf-8') as fh:
        ph = list(csv.DictReader(fh))
    seqs = [p['sequence'] for p in ph]
    n = max(len(s) for s in seqs)
    M = np.full((len(seqs), n), np.nan)
    for i, s in enumerate(seqs):
        for j, ch in enumerate(s):
            M[i, j] = np.nan if ch == 'x' else int(ch)

    fig, ax = plt.subplots(figsize=(13, 6.5), layout='constrained')
    cmap = ListedColormap(plt.cm.viridis(np.linspace(0, 1, 8)))
    cmap.set_bad('#f2f2f2')
    ax.imshow(M, cmap=cmap, vmin=-.5, vmax=7.5, aspect='auto')
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if np.isnan(M[i, j]):
                ax.text(j, i, 'x', ha='center', va='center', fontsize=11,
                        color=CLR['bad'], fontweight='bold')
    ax.set_yticks(range(len(ph)))
    ax.set_yticklabels([f"phase {p['phase_index']}  (t0 = "
                        f"{float(p['source_time_h']):.1f} h)" for p in ph],
                       fontsize=12)
    ax.set_xlabel('steady read window index (after drop-8)', fontsize=15)
    ax.set_title('Eight orbit-phase initial states, all from ONE converged '
                 'period at clock_K = 0.10\n'
                 'every phase reproduces the identical mod-8 sequence'
                 if all(p['full'] == 'True' for p in ph) else
                 'Eight orbit-phase initial states at clock_K = 0.10',
                 fontsize=15)
    ax.tick_params(labelsize=12)
    cb = fig.colorbar(ax.images[0], ax=ax, ticks=range(8), pad=.02)
    cb.set_label('decoded value', fontsize=14)
    cb.ax.tick_params(labelsize=12)
    for ext in ('png', 'svg'):
        fig.savefig(FIG / f'fig4_eight_phases.{ext}', dpi=150)
    plt.close(fig)
    print('fig4 done')

    print('FIGURES', FIG)


if __name__ == '__main__':
    main()
