"""Round-4 figures: the pulse_gate bypass, and the eight-phase event accounting.

  fig5_bypass        clock_K = 0.40 with pulse_gate present vs removed
  fig6_phase_events  per-phase reverse / gate / flip counts for both stages,
                     with the single failing event located

Reads saved trajectories and result JSON only. English labels, large fonts,
thick lines. Writes into failure_attribution/figures/.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

sys.dont_write_bytecode = True

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT.parent))
sys.path.insert(0, str(ROOT.parent.parent))

import diagnose_round1 as D                                       # noqa: E402
import verify_hzh as V                                            # noqa: E402

FIG = ROOT / 'figures'
RES = ROOT / 'results'
CLR = dict(s0='#2685ad', s1='#d97726', s2='#43875a', ok='#2e7d32',
           bad='#c62828', gate='#d97726', clock='#6a1b9a',
           int0='#7550a1', int2='#2685ad', rdf2='#b8860b',
           fwd='#2685ad', rev='#c62828')


def newest(pat):
    d = sorted(RES.glob(pat), key=lambda p: p.name)
    if not d:
        raise SystemExit(f'no results for {pat}')
    return d[-1]


def main():
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8')
    FIG.mkdir(exist_ok=True)
    r2 = newest('round2_*')
    r4 = newest('round4_*')
    r4b = newest('round4b_*')
    han = D.HanInput(D.HOURS)
    prefix = D.HZHModel('han', han)

    # =====================================================================
    # fig 5: pulse_gate present vs removed, both at clock_K = 0.40
    # =====================================================================
    zB = np.load(r2 / 'traj_K0.40_off.npz')
    zA = np.load(r4b / 'traj_bypass_K0.40_off.npz')
    tB, yB = zB['time_h'], zB['states']
    tA, yA = zA['time_h'], zA['states']
    mB = D.HZZModel(prefix, False); mB.clock_K = 0.40
    mA = D.HZZModel(prefix, False); mA.clock_K = 0.40
    sB, sA = mB.signals(yB), mA.signals(yA)

    win = (150.0, 215.0)
    kk = lambda t: (t >= win[0]) & (t <= win[1])

    fig, axes = plt.subplots(3, 2, figsize=(15, 11), sharex=True,
                             layout='constrained')
    for col, (t, y, s, tag, col_ok) in enumerate([
            (tB, yB, sB, 'pulse_gate PRESENT  (as published)  -> fails', False),
            (tA, yA, sA, 'pulse_gate REMOVED  (clock term = 1)  -> counts', True)]):
        m = kk(t)
        xt = t[m]
        rr = V.read_windows(t, s)[1]
        seq = V.read_verdict(rr, V.DROP)['sequence']

        ax = axes[0, col]
        ax.set_title(tag, fontsize=15,
                     color=CLR['ok'] if col_ok else CLR['bad'])
        for b, c in (('S0', CLR['s0']), ('S1', CLR['s1']), ('S2', CLR['s2'])):
            ax.plot(xt, s[b][m], color=c, lw=2.4, label=b)
        ax.axhspan(0, V.LOW, color=CLR['s2'], alpha=.09)
        ax.axhspan(V.HIGH, 1, color=CLR['s0'], alpha=.09)
        ax.axhline(V.LOW, ls='--', color='gray', lw=1.2)
        ax.axhline(V.HIGH, ls='--', color='gray', lw=1.2)
        ax.set_ylim(-.05, 1.05)
        ax.grid(alpha=.22); ax.tick_params(labelsize=13)
        ax.set_ylabel('DNA LR fraction', fontsize=15)
        ax.legend(fontsize=12, ncol=3, loc='lower right')
        ax.text(.01, .04, 'steady: ' + seq, transform=ax.transAxes, fontsize=13,
                color=CLR['ok'] if col_ok else CLR['bad'], fontweight='bold')

        ax = axes[1, col]
        ax.plot(xt, s['g1'][m], color=CLR['gate'], lw=2.4, label='g1 (gate)')
        ax.axhline(V.GATE_THRESHOLD, color='gray', ls=':', lw=1.6)
        ax.set_ylabel('gate g1', fontsize=15)
        a2 = ax.twinx()
        a2.plot(xt, s['clock'][m], color=CLR['clock'], lw=2.2, ls='--',
                label='clock factor')
        a2.set_ylim(-.05, 1.08)
        a2.set_ylabel('clock factor', fontsize=15)
        a2.tick_params(labelsize=13)
        ln = ax.get_lines()[:1] + a2.get_lines()
        ax.legend(ln, [l.get_label() for l in ln], fontsize=12, loc='upper right')
        ax.grid(alpha=.22); ax.tick_params(labelsize=13)

        ax = axes[2, col]
        ax.plot(xt, y[20][m], color=CLR['int2'], lw=2.4, label='mature Int2')
        ax.set_ylabel('Int2 (a.u.)', fontsize=15)
        a3 = ax.twinx()
        a3.plot(xt, y[22][m], color=CLR['rdf2'], lw=2.4, label='RDF2')
        a3.set_ylabel('RDF2 (a.u.)', fontsize=15)
        a3.tick_params(labelsize=13)
        ln = ax.get_lines() + a3.get_lines()
        ax.legend(ln, [l.get_label() for l in ln], fontsize=12, loc='upper right')
        ax.grid(alpha=.22); ax.tick_params(labelsize=13)
        ax.set_xlabel('Time (h)', fontsize=15)
        for a in axes[:, col]:
            a.set_xlim(*win)

    fig.suptitle('pulse_gate must not act as a discriminator against the real '
                 'Int0\n'
                 'clock_K left at its source value 0.40; only the third factor '
                 'of the gate is changed', fontsize=16)
    for ext in ('png', 'svg'):
        fig.savefig(FIG / f'fig5_bypass.{ext}', dpi=150)
    plt.close(fig)
    print('fig5 done')

    # =====================================================================
    # fig 6: per-phase event accounting
    # =====================================================================
    rev = json.loads((RES / f'phase_reeval_{r4.name}.json').read_text(encoding='utf-8'))
    unm = json.loads((RES / f'phase_unmatched_both_stages_{r4.name}.json')
                     .read_text(encoding='utf-8'))
    rows = rev['rows']
    n = len(rows)
    x = np.arange(n)
    fig, axes = plt.subplots(2, 1, figsize=(13, 9), sharex=True,
                             layout='constrained')
    for ax, stage, key, lab in (
            (axes[0], 'stage0', 'stage0', 'bit0 -> bit1'),
            (axes[1], 'stage1', 'stage1', 'bit1 -> bit2')):
        w = 0.26
        r_ = [r[key][0] for r in rows]
        g_ = [r[key][1] for r in rows]
        f_ = [r[key][2] for r in rows]
        ax.bar(x - w, r_, w, color='#6a1b9a', label='reverse events')
        ax.bar(x, g_, w, color=CLR['gate'], label='gate events')
        ax.bar(x + w, f_, w, color=CLR['s2'], label='bit flips')
        for i in range(n):
            if not rows[i]['events']:
                ax.axvspan(i - .42, i + .42, color=CLR['bad'], alpha=.13, zorder=0)
        ax.set_ylabel('event count', fontsize=15)
        ax.set_title(f'{lab}', fontsize=14)
        ax.grid(alpha=.22, axis='y')
        ax.tick_params(labelsize=13)
        if stage == 'stage0':
            ax.legend(fontsize=12, ncol=3, loc='upper right')
    axes[1].set_xticks(x)
    axes[1].set_xticklabels([f"phase {r['phase']}\n{r['events'] and 'pass' or 'FAIL'}"
                             for r in rows], fontsize=12)
    axes[1].set_xlabel('orbit-phase initial state', fontsize=15)
    axes[1].set_ylim(0, max(max(r['stage1']) for r in rows) * 1.45)

    bad = unm['phase7']['stage0']['unmatched']
    offend = [u for u in bad if not u['in_dropped_burnin']]
    if offend:
        u = offend[0]
        axes[0].annotate(
            f"reverse #{u['reverse_id']} at {u['reverse_peak_h']:.1f} h:\n"
            f"gate opened, bit1 did NOT flip\n"
            f"({u['gap_to_end_h']:.0f} h before the window end, so this is "
            f"NOT a window-edge effect)",
            xy=(7 - 0.26, rows[7]['stage0'][0]),
            xytext=(3.1, rows[7]['stage0'][0] * 0.62),
            fontsize=12, color=CLR['bad'],
            arrowprops=dict(arrowstyle='->', color=CLR['bad'], lw=2.2))
    axes[1].text(.5, .06,
                 'counting passes in 8/8 phases; the event arm fails in 1/8, '
                 'and the failure is at stage0, not at bit2',
                 transform=axes[1].transAxes, ha='center', fontsize=13,
                 fontweight='bold')
    fig.suptitle('Eight orbit-phase initial states at clock_K = 0.10, 400 h\n'
                 'shaded = the phase whose event arm fails', fontsize=15)
    for ext in ('png', 'svg'):
        fig.savefig(FIG / f'fig6_phase_events.{ext}', dpi=150)
    plt.close(fig)
    print('fig6 done')
    print('FIGURES', FIG)


if __name__ == '__main__':
    main()
