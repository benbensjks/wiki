"""M-26: two-stage carry waveform and causal timing figure for HBY-ZMH-HBY.

Draws one complete second-stage carry period (97.0-139.4 h) of the CERTIFIED
HZH 34-state run, first stage above and second stage below on a shared time axis.

Variants:
  *_annotated  curves + event verticals + causal arrows + measured offsets
  *_clean      curves only

Direction convention (important, and the first version got it wrong):
  S0 = b0_S and S1 = 1 - pb1_zmh are LR fractions.  A carry is the bit going
  LR -> PB, i.e. the S signal crossing 0.5 DOWNWARD.  `associations[].flip_h`
  records every 0.5 crossing of the next bit, alternating down and up, so only
  the direction=='down' rows are carries.  Cross-checked against
  bit1_to_bit2.associations[].reverse_peak_h, which lists the same events.

Everything is read from the saved certified trajectory and the certified event
table; no re-integration.  The script asserts that re-analysing the saved
trajectory reproduces certified_v1 / 19 reads / the certified sequence, and that
the window holds exactly two low-order carries and one high-order carry.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

sys.dont_write_bytecode = True

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent
for _d in (PARENT, PARENT / 'hby_zmh_hby',
           PARENT / 'hby_zmh_hby' / 'certification'):
    sys.path.insert(0, str(_d))

import verify_hzh as V                                 # noqa: E402
from model_hzh import HZHModel, NAMES, IDX             # noqa: E402
from run_han_comparison import HanInput                # noqa: E402

CERT = (PARENT / 'hby_zmh_hby' / 'certification' / 'results'
        / '20260928_141723_738780' / 'baseline_saved' / 'verdict.json')
TRAJ = (PARENT / 'hby_zmh_hby' / 'results'
        / '20260927_233646_924908' / 'han' / 'trajectory.npz')

WIN = (97.0, 139.4)
HOURS = 300.0

C = dict(pb='#1f4e79', A='#d97726', F='#2e7d32', g='#6a1b9a',
         Int='#8d6e14', clk='#6a1b9a', low='#c2185b', gate='#ef6c00',
         hinge='#111111', ink='#3a3f42')

R = dict(pb0=0, af0=1, g0=2, int1=3, pb1=4, af1=5, g1=6, int2=7)
ROWS = ['PB0  (DNA)', 'A0 and F0  (normalised)', 'g0  (gate)', 'Int1  (a.u.)',
        'PB1  (DNA)', 'A1 and F1  (normalised)', 'g1  (gate)', 'Int2  (a.u.)']


def half_cross(t, x, lo, hi):
    k = (t >= lo) & (t <= hi)
    if k.sum() < 5:
        return None, None
    tt, xx = t[k], np.asarray(x, dtype=float)[k]
    amp = float(xx.max())
    if amp <= 0:
        return None, amp
    hit = tt[xx >= 0.5 * amp]
    return (float(hit.min()) if hit.size else None), amp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--variant', default='both',
                    choices=['annotated', 'clean', 'both'])
    args = ap.parse_args()
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8')
    out = HERE / 'out'
    out.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------ data
    verdict = json.loads(CERT.read_text(encoding='utf-8'))
    m = HZHModel('han', HanInput(HOURS))
    z = np.load(TRAJ)
    t, y = z['time_h'], z['states']
    assert tuple(z['state_names']) == tuple(NAMES)

    vv, _ = V.analyse(t, y, m.z, m.tail)
    guard = dict(certified=bool(vv['certified_v1']),
                 reads=int(vv['steady']['reads']) == 19,
                 sequence=vv['steady']['sequence'] == '1234567012345670123',
                 margin=abs(vv['global_min_timing_margin_h']
                            - 1.1784609018563117) < 1e-6)
    print('certified reproduction guard:', json.dumps(guard), flush=True)
    if not all(guard.values()):
        raise SystemExit(f'saved trajectory no longer reproduces the certified '
                         f'verdict: {guard}')

    sig = V.signals(y, m.z, m.tail)
    sig['int0'] = y[IDX['b0_I']]
    A0t = verdict['events']['bit0_to_bit1']['associations']
    A1t = verdict['events']['bit1_to_bit2']['associations']

    lo, hi = WIN
    ins = lambda x: lo <= x <= hi

    # carries: LR -> PB is the DOWNWARD 0.5 crossing of an LR-fraction signal
    s0_down = [c['time_h'] for c in V.crossings(t, sig['S0'])
               if c['direction'] == 'down' and ins(c['time_h'])]
    s1_down = [c['time_h'] for c in V.crossings(t, sig['S1'])
               if c['direction'] == 'down' and ins(c['time_h'])]
    hinge_tab = [a['flip_h'] for a in A0t
                 if a['flip_h'] is not None
                 and a.get('flip_direction') == 'down' and ins(a['flip_h'])]
    hinge_jrev = [a['reverse_peak_h'] for a in A1t
                  if a['reverse_peak_h'] is not None and ins(a['reverse_peak_h'])]
    g0_start = [a['gate_start_h'] for a in A0t
                if a['gate_start_h'] is not None and ins(a['gate_start_h'])]
    g1_start = [a['gate_start_h'] for a in A1t
                if a['gate_start_h'] is not None and ins(a['gate_start_h'])]
    x2 = [c['time_h'] for c in V.crossings(t, sig['S2']) if ins(c['time_h'])]

    counts = dict(low_order=len(s0_down), g0_open=len(g0_start),
                  high_order=len(hinge_tab), g1_open=len(g1_start),
                  bit2_flip=len(x2))
    print('window event counts:', json.dumps(counts), flush=True)
    print('  S0 down :', [round(v, 3) for v in s0_down], flush=True)
    print('  S1 down :', [round(v, 3) for v in s1_down], flush=True)
    print('  A0t down-flip :', [round(v, 3) for v in hinge_tab], flush=True)
    print('  A1t J_rev1 peak:', [round(v, 3) for v in hinge_jrev], flush=True)
    print('  g0 open :', [round(v, 3) for v in g0_start], flush=True)
    print('  g1 open :', [round(v, 3) for v in g1_start], flush=True)
    print('  S2 cross:', [round(v, 3) for v in x2], flush=True)
    if len(s0_down) != 2 or len(x2) != 1 or len(hinge_tab) != 1:
        raise SystemExit('window must hold exactly 2 low-order carries, '
                         f'1 high-order carry and 1 bit2 flip: {counts}')
    if len(hinge_jrev) != 1:
        raise SystemExit('bit1 reverse-flux list disagrees with the S1 '
                         f'down-crossing list: {hinge_jrev}')

    # --------------------------------------------------- measured quantities
    def segs(sig_key):
        s = V.segments(t, sig[sig_key], V.GATE_THRESHOLD,
                       min_dose=V.GATE_MIN_DOSE_H)
        return [q for q in s if ins(q['peak_h'])]

    seg0, seg1 = segs('g0'), segs('g1')

    def stage_timing(Acol, Fcol, ss):
        A = np.asarray(y[Acol], dtype=float)
        F = np.asarray(y[Fcol], dtype=float)
        rows = []
        for s in ss:
            ta, ampa = half_cross(t, A, s['start_h'] - 8, s['peak_h'] + 14)
            tf, ampf = half_cross(t, F, s['start_h'] - 8, s['peak_h'] + 14)
            rows.append(dict(gate_start_h=s['start_h'], gate_peak_h=s['peak_h'],
                             gate_peak=s['peak'], gate_dose=s['dose'],
                             gate_duration_h=s['end_h'] - s['start_h'],
                             A_half_h=ta, F_half_h=tf,
                             A_to_F_h=None if None in (ta, tf) else tf - ta,
                             A_max=ampa, F_max=ampf))
        return rows

    tim0 = stage_timing(IDX['A0_zmh'], IDX['F0_zmh'], seg0)
    tim1 = stage_timing(IDX['A1'], IDX['F1'], seg1)
    med = lambda v: float(np.median([x for x in v if x is not None]))

    # why the second low-order carry does not propagate: bit1 is already PB
    gate_ctx = []
    for x in g0_start:
        j = int(np.argmin(np.abs(t - x)))
        s1 = float(sig['S1'][j])
        gate_ctx.append(dict(gate_open_h=x, S1=s1,
                             bit1_state=('PB' if s1 <= V.LOW else
                                         'LR' if s1 >= V.HIGH else 'undecided')))
    nxt = [c['time_h'] for c in V.crossings(t, sig['S1'])
           if c['direction'] == 'down' and c['time_h'] > g0_start[-1]]
    print('  bit1 state at each g0 opening:',
          [(round(q['gate_open_h'], 2), q['bit1_state']) for q in gate_ctx],
          '| next bit1 LR→PB at', round(nxt[0], 3) if nxt else None, flush=True)
    if (len(gate_ctx) > 1
            and gate_ctx[0]['bit1_state'] != 'LR'
            or (len(gate_ctx) > 1 and gate_ctx[-1]['bit1_state'] != 'PB')):
        raise SystemExit('the "gate opens but bit1 is already PB" annotation is '
                         f'not supported by the trajectory: {gate_ctx}')

    summary = dict(
        window=list(WIN), counts=counts,
        gate_context=gate_ctx,
        next_bit1_carry_after_last_gate=(nxt[0] if nxt else None),
        convention=('S0 and S1 are LR fractions; a carry is the DOWNWARD 0.5 '
                    'crossing. associations[].flip_h lists every crossing, so '
                    'only direction=="down" rows are carries.'),
        stage0=dict(A_to_F_h=med([r['A_to_F_h'] for r in tim0]),
                    gate_peak=med([r['gate_peak'] for r in tim0]),
                    gate_duration_h=med([r['gate_duration_h'] for r in tim0]),
                    gate_dose=med([r['gate_dose'] for r in tim0]), detail=tim0),
        stage1=dict(A_to_F_h=med([r['A_to_F_h'] for r in tim1]),
                    gate_peak=med([r['gate_peak'] for r in tim1]),
                    gate_duration_h=med([r['gate_duration_h'] for r in tim1]),
                    gate_dose=med([r['gate_dose'] for r in tim1]), detail=tim1),
        event_lines=dict(bit0_LR_to_PB=s0_down, g0_open=g0_start,
                         bit1_LR_to_PB=hinge_tab,
                         bit1_LR_to_PB_from_Jrev1=hinge_jrev,
                         g1_open=g1_start, bit2_flip=x2),
        source=dict(cert=str(CERT.relative_to(PARENT)),
                    traj=str(TRAJ.relative_to(PARENT)), guard=guard))
    (out / 'm26_measurements.json').write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
    print('Δt(A→F): stage0 %.3f h, stage1 %.3f h'
          % (summary['stage0']['A_to_F_h'], summary['stage1']['A_to_F_h']),
          flush=True)

    # -------------------------------------------------------------- drawing
    for variant in (['annotated', 'clean'] if args.variant == 'both'
                    else [args.variant]):
        ann = variant == 'annotated'
        fig, axes = plt.subplots(8, 1, figsize=(14, 16), sharex=True,
                                 gridspec_kw=dict(height_ratios=[1, 1.15, 1, 1,
                                                                 1, 1.15, 1, 1]))
        fig.subplots_adjust(left=0.085, right=0.915, top=0.918, bottom=0.150,
                            hspace=0.22)
        k = (t >= lo) & (t <= hi)
        tt = t[k]

        axes[R['pb0']].plot(tt, 1 - np.asarray(y[IDX['b0_S']], dtype=float)[k],
                            color=C['pb'], lw=2.6, label='PB0 = 1 - b0_S')
        A0, F0 = (np.asarray(y[IDX['A0_zmh']], dtype=float),
                  np.asarray(y[IDX['F0_zmh']], dtype=float))
        nA0, nF0 = float(A0[k].max()), float(F0[k].max())
        axes[R['af0']].plot(tt, A0[k] / nA0, color=C['A'], lw=2.6,
                            label=f'A0 / {nA0:.3f}')
        axes[R['af0']].plot(tt, F0[k] / nF0, color=C['F'], lw=2.6,
                            label=f'F0 / {nF0:.3f}')
        axes[R['g0']].plot(tt, np.asarray(sig['g0'], dtype=float)[k],
                           color=C['g'], lw=2.8, label='g0')
        axes[R['int1']].plot(tt, np.asarray(y[IDX['I1_zmh']], dtype=float)[k],
                             color=C['Int'], lw=2.6, label='Int1 = I1_zmh')

        axes[R['pb1']].plot(tt, np.asarray(y[IDX['pb1_zmh']], dtype=float)[k],
                            color=C['pb'], lw=2.6, label='PB1 = pb1_zmh')
        A1, F1 = (np.asarray(y[IDX['A1']], dtype=float),
                  np.asarray(y[IDX['F1']], dtype=float))
        nA1, nF1 = float(A1[k].max()), float(F1[k].max())
        axes[R['af1']].plot(tt, A1[k] / nA1, color=C['A'], lw=2.6,
                            label=f'A1 / {nA1:.3f}')
        axes[R['af1']].plot(tt, F1[k] / nF1, color=C['F'], lw=2.6,
                            label=f'F1 / {nF1:.3f}')
        axes[R['g1']].plot(tt, np.asarray(sig['g1'], dtype=float)[k],
                           color=C['g'], lw=2.8, label='g1')
        aclk = axes[R['g1']].twinx()
        aclk.plot(tt, np.asarray(sig['clock'], dtype=float)[k], color=C['clk'],
                  lw=1.6, ls='--', alpha=.75, label='clock factor')
        aclk.set_ylim(-0.02, 1.05)
        aclk.set_ylabel('clock', fontsize=11, color=C['clk'])
        aclk.tick_params(labelsize=10, colors=C['clk'])
        axes[R['int2']].plot(tt, np.asarray(y[IDX['b2_I']], dtype=float)[k],
                             color=C['Int'], lw=2.6, label='Int2 = b2_I')

        for r, lab in enumerate(ROWS):
            axes[r].set_ylabel(lab, fontsize=12)
            axes[r].grid(alpha=.18)
            axes[r].tick_params(labelsize=11)
            axes[r].legend(fontsize=10.5, loc='upper right', framealpha=.92)
        axes[R['af0']].legend(fontsize=10.5, loc='upper left', framealpha=.92)
        axes[R['af1']].legend(fontsize=10.5, loc='upper left', framealpha=.92)
        aclk.legend(fontsize=10, loc='lower right', framealpha=.92)
        for r in (R['pb0'], R['pb1']):
            axes[r].set_ylim(-.03, 1.03)
        for r in (R['af0'], R['af1']):
            axes[r].set_ylim(-.03, 1.08)
        for r in (R['g0'], R['g1']):
            axes[r].set_ylim(-.015, 0.68)
            axes[r].axhline(V.GATE_THRESHOLD, color='gray', ls=':', lw=1.4)
            axes[r].text(WIN[0] + 0.25, V.GATE_THRESHOLD + .012,
                         'event threshold 0.05', fontsize=9.5, color='gray')
        for r in (R['int1'], R['int2']):
            axes[r].set_ylim(-.06, 3.05)
        axes[-1].set_xlabel('Time (h)', fontsize=13)
        axes[-1].set_xlim(*WIN)
        axes[-1].set_xticks(np.arange(100, 140, 5))
        axes[R['pb0']].text(0.006, 0.86, 'FIRST CARRY   bit0 → bit1',
                            transform=axes[R['pb0']].transAxes, fontsize=12,
                            fontweight='bold', color=C['ink'])
        axes[R['pb1']].text(0.006, 0.86, 'SECOND CARRY   bit1 → bit2',
                            transform=axes[R['pb1']].transAxes, fontsize=12,
                            fontweight='bold', color=C['ink'])

        if ann:
            for r in range(0, 4):
                for x in s0_down:
                    axes[r].axvline(x, color=C['low'], lw=1.8, alpha=.85, zorder=1)
                for x in g0_start:
                    axes[r].axvline(x, color=C['gate'], lw=1.8, ls='--',
                                    alpha=.9, zorder=1)
            for r in range(4, 8):
                for x in g1_start:
                    axes[r].axvline(x, color=C['gate'], lw=1.8, ls='--',
                                    alpha=.9, zorder=1)
                for x in x2:
                    axes[r].axvline(x, color=C['hinge'], lw=1.8, ls=':',
                                    alpha=.85, zorder=1)
            for r in range(8):
                axes[r].axvline(hinge_tab[0], color=C['hinge'], lw=1.2,
                                alpha=.32, zorder=0)

            # continuous hinge line through all eight rows
            fig.canvas.draw()
            xf = axes[R['pb0']].transData.transform((hinge_tab[0], 0))[0]
            xf /= fig.bbox.width
            fig.add_artist(Line2D([xf, xf],
                                  [axes[-1].get_position().y0,
                                   axes[R['pb0']].get_position().y1],
                                  transform=fig.transFigure, color=C['hinge'],
                                  lw=2.4, alpha=.92))

            handles = [Line2D([], [], color=C['low'], lw=1.8,
                              label='bit0 LR→PB'),
                       Line2D([], [], color=C['gate'], lw=1.8, ls='--',
                              label='carry gate opens'),
                       Line2D([], [], color=C['hinge'], lw=2.4,
                              label='bit1 LR→PB  (stage-0 out = stage-1 in)'),
                       Line2D([], [], color=C['hinge'], lw=1.8, ls=':',
                              label='bit2 flips')]
            axes[R['pb0']].legend(handles=handles, fontsize=10.5,
                                  loc='lower left', framealpha=.95)

            # arrows
            axes[R['pb0']].annotate(
                '', xy=(g0_start[0] - 0.35, 0.05), xytext=(g0_start[0] - 1.1, .80),
                arrowprops=dict(arrowstyle='->', color=C['ink'], lw=1.6,
                                ls='--', shrinkA=1, shrinkB=1))
            axes[R['af0']].annotate(
                'PB0 drives A0', xy=(g0_start[0] - 1.3, 0.32),
                xytext=(g0_start[0] - 11.5, 0.60), fontsize=10.5, color=C['ink'],
                arrowprops=dict(arrowstyle='->', color=C['ink'], lw=1.5, ls='--'))
            p0 = seg0[0]['peak_h']
            axes[R['g0']].annotate(
                '', xy=(p0, 0.52), xytext=(p0, 0.30),
                arrowprops=dict(arrowstyle='->', color=C['ink'], lw=1.7))
            axes[R['g0']].annotate(
                'one gate pulse → one Int1 pulse', xy=(p0 + 0.5, 0.33),
                xytext=(p0 + 2.6, 0.50), fontsize=10.5, color=C['ink'],
                arrowprops=dict(arrowstyle='->', color=C['ink'], lw=1.4, ls='--'))
            axes[R['int2']].annotate(
                '', xy=(x2[0] - 0.2, 2.35), xytext=(hinge_tab[0] + 0.2, 2.35),
                arrowprops=dict(arrowstyle='->', color=C['ink'], lw=1.7))
            axes[R['int2']].annotate(
                'bit1 carry → bit2 flips %.2f h later' % (x2[0] - hinge_tab[0]),
                xy=(hinge_tab[0] + 0.35, 2.45),
                xytext=(hinge_tab[0] + 2.8, 2.72), fontsize=10.5, color=C['ink'],
                arrowprops=dict(arrowstyle='->', color=C['ink'], lw=1.3, ls='--'))

            for r, tim, tag in ((R['af0'], tim0, 'stage 0: ZMH reduced A0/F0'),
                                (R['af1'], tim1, 'stage 1: HBY A1/F1 chain')):
                q = tim[0]
                ta, tf = q['A_half_h'], q['F_half_h']
                if ta is None or tf is None:
                    continue
                yb = 0.88
                axes[r].plot([ta], [0.5], marker='o', ms=7, mfc='none',
                             mec=C['A'], mew=2.0, zorder=6)
                axes[r].plot([tf], [0.5], marker='o', ms=7, mfc='none',
                             mec=C['F'], mew=2.0, zorder=6)
                for xv in (ta, tf):
                    axes[r].plot([xv, xv], [0.5, yb], color='gray', lw=1.1,
                                 alpha=.7, zorder=5)
                axes[r].annotate('', xy=(ta, yb), xytext=(tf, yb),
                                 arrowprops=dict(arrowstyle='<->', color='gray',
                                                 lw=1.7))
                axes[r].text((ta + tf) / 2 + 0.6, yb + 0.02,
                             'Δt(A→F) = %.3f h' % (tf - ta), fontsize=10.5,
                             color='#444444')
                axes[r].text(0.997, 0.06, tag, transform=axes[r].transAxes,
                             ha='right', fontsize=10, color='#707070')

            if len(g0_start) > 1:
                axes[R['g0']].annotate(
                    'gate opens; no bit1 flip follows\n(bit1 is already PB)',
                    xy=(g0_start[1] - 0.4, 0.52),
                    xytext=(g0_start[1] - 18.0, 0.46), fontsize=10.5,
                    color='#8a5a00',
                    arrowprops=dict(arrowstyle='->', color='#8a5a00', lw=1.5))

            s0, s1 = summary['stage0'], summary['stage1']
            box = ('window %.1f–%.1f h = one complete second-stage carry period'
                   '   |   low-order carries %d, high-order carry %d\n'
                   'Δt(A→F)  stage 0  %.3f h   (ZMH reduced A0/F0)    '
                   'stage 1  %.3f h   (HBY A1/F1, +%.0f%%)\n'
                   'gate peak  g0 %.3f   g1 %.3f   |   g0 opens %.2f h before '
                   'the J_rev0 peak;  g1 opens %.2f h after the J_rev1 peak'
                   % (WIN[0], WIN[1], counts['low_order'], counts['bit2_flip'],
                      s0['A_to_F_h'], s1['A_to_F_h'],
                      100 * (s1['A_to_F_h'] / s0['A_to_F_h'] - 1),
                      s0['gate_peak'], s1['gate_peak'],
                      g0_start[0] - [a['reverse_peak_h'] for a in A0t
                                     if ins(a['reverse_peak_h'])][0],
                      g1_start[0] - hinge_jrev[0]))
            fig.text(0.5, 0.022, box, ha='center', va='bottom', fontsize=10,
                     color='#444444', linespacing=1.6)

        fig.suptitle('Two-stage carry timing, HBY–ZMH–HBY (HZH, 34-state)',
                     fontsize=15, y=0.982)
        fig.text(0.5, 0.955,
                 'Han v53d upstream (peak_load_fraction = 0); frozen '
                 'selected_v1 + HbyConfig; 300 h, 1 min grid, DOP853 '
                 'rtol 2e-7/atol 2e-9; criterion HZH_MOD8_CAUSAL_V1.',
                 ha='center', fontsize=9.5, color='#6a6a6a')
        fig.text(0.5, 0.940,
                 'A/F curves are each normalised to their own maximum inside '
                 'the window; the number after the slash in each legend entry '
                 'is that maximum.',
                 ha='center', fontsize=9.5, color='#6a6a6a')

        fig.canvas.draw()
        bb = fig.bbox
        bad = []
        for txt in fig.findobj(matplotlib.text.Text):
            if not txt.get_text().strip() or not txt.get_visible():
                continue
            try:
                e = txt.get_window_extent(fig.canvas.get_renderer())
            except Exception:
                continue
            if (e.x0 < bb.x0 - 1 or e.x1 > bb.x1 + 1
                    or e.y0 < bb.y0 - 1 or e.y1 > bb.y1 + 1):
                bad.append((txt.get_text()[:44], round(e.x0), round(e.x1),
                            round(e.y0), round(e.y1)))
        print(f'[{variant}] overflow check: {len(bad)} text item(s) outside the '
              f'canvas', flush=True)
        for b in bad[:10]:
            print('    ', b, flush=True)

        stem = out / f'M26_two_stage_carry_timing_{variant}'
        for ext in ('png', 'svg', 'pdf'):
            fig.savefig(f'{stem}.{ext}', dpi=160)
        plt.close(fig)
        print(f'[{variant}] wrote {stem.name}.*', flush=True)

    print('OUTPUT', out, flush=True)


if __name__ == '__main__':
    main()
