"""M-21 three-arm figure: kept / removed / removed-with-total-input-matched.

Panels (4; the 3-panel guideline is deliberately exceeded because the reviewer
asked for all four quantities on one page):
  (a) the Int2 source waveform u_2(t) = alpha_Int[1]*g1(t), one carry window
  (b) the resulting mature Int2, same window
  (c) bit2's distance to the 0.30/0.70 band edge, per steady read window
  (d) far-dwell gate activity per carry cycle, with the in-pulse dose beside it

Reads saved trajectories only -- no integration, no smoothing, no resampling.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'DejaVu Sans']
matplotlib.rcParams['axes.unicode_minus'] = False

sys.dont_write_bytecode = True

HERE = Path(__file__).resolve().parent
SCAN = HERE.parent / '时钟工作区扫描'
for _d in (str(SCAN), str(HERE)):
    if _d not in sys.path:
        sys.path.insert(0, _d)

import common as C                                                # noqa: E402
from run_m21_clock_bypass import BypassReceiver, gate_no_clock_array  # noqa: E402

DEEP_BLUE, ORANGE, TEAL, PURPLE = '#304B53', '#D29144', '#4F9194', '#6a1b9a'
GRAY = '#8ca0a5'
STYLE = dict(kept=dict(c=TEAL, ls='-', lab='KEPT  (clock present)'),
             removed=dict(c=ORANGE, ls='--', lab='REMOVED  (clock = 1, unmatched)'),
             matched=dict(c=PURPLE, ls='-.', lab='MATCHED  (clock = 1, total input matched)'))
ALPHA0 = 38.0
ALPHA_M = 23.435159
SRC2 = 'm21_20261001_225706'


def newest(pat):
    d = sorted([p for p in (HERE / 'results').glob(pat)
                if p.is_dir() and 'SMOKE' not in p.name and 'SUPERSEDED' not in p.name])
    if not d:
        raise SystemExit('no results for ' + pat)
    return d[-1]


def main():
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8')
    d2 = HERE / 'results' / SRC2
    dm = newest('m21matched_*')
    print('2-arm:', d2.name, ' matched:', dm.name, flush=True)

    han = C.HanInput(300.0)
    mk = C.HZHModel('han', han)
    mr = C.HZHModel('han', han); mr.tail = BypassReceiver(mr.tail.c)
    mm = C.HZHModel('han', han)
    recv = BypassReceiver(mm.tail.c); recv.p = dict(recv.p)
    s = list(recv.p['alpha_Int']); s[1] = ALPHA_M; recv.p['alpha_Int'] = tuple(s)
    mm.tail = recv

    def load(path, model, bp, alpha):
        z = np.load(path)
        t, y = z['time_h'], z['states']
        sig = C.sig_of(model, y, 'hzh')
        if bp:
            g = gate_no_clock_array(y, model)
            sig['clock'] = np.ones_like(np.asarray(sig['clock'], dtype=float))
            sig['g1'] = g
            sig['_u2'] = g * float(alpha)
        return t, y, sig

    tk, yk, sk = load(d2 / 'traj_kept_K0.30.npz', mk, False, ALPHA0)
    tr, yr, sr = load(d2 / 'traj_removed_K0.30.npz', mr, True, ALPHA0)
    tm, ym, sm = load(dm / 'traj_matched_K0.30.npz', mm, True, ALPHA_M)

    DROP = int(C.VH.DROP)
    _, reads = C.VH.read_windows(tk, sk)
    sel = reads[DROP:]
    win8 = reads[DROP:DROP + 8]                       # the matching window
    tw = (240.0, 250.0)                               # one steady carry window
    m = (tk >= tw[0]) & (tk <= tw[1])

    def band_per_window(t, y):
        S2 = y[C.VH.IDX['b2_S']]
        out = []
        for r in sel:
            mm_ = (t >= r['start_h']) & (t <= r['end_h'])
            v = S2[mm_]
            lab = r['labels'][2]
            out.append(C.VH.LOW - float(v.max()) if lab == 0
                       else float(v.min()) - C.VH.HIGH if lab == 1
                       else float('nan'))
        return np.array(out)

    bk, br, bm = band_per_window(tk, yk), band_per_window(tr, yr), band_per_window(tm, ym)
    mech = dict(kept=C.mechanism(tk, sk, sel), removed=C.mechanism(tr, sr, sel),
                matched=C.mechanism(tm, sm, sel))

    fig, axes = plt.subplots(2, 2, figsize=(17.5, 11.5), layout='constrained')

    # ---------------- (a) Int2 source waveform ----------------
    ax = axes[0][0]
    for k in ('kept', 'removed', 'matched'):
        t_, sig_ = {'kept': (tk, sk), 'removed': (tr, sr), 'matched': (tm, sm)}[k]
        a = {'kept': ALPHA0, 'removed': ALPHA0, 'matched': ALPHA_M}[k]
        ax.plot(t_[m], a * np.asarray(sig_['g1'], float)[m], color=STYLE[k]['c'],
                ls=STYLE[k]['ls'], lw=2.4, label=STYLE[k]['lab'])
    ax.set_ylabel(r'Int2 source  $u_2=\alpha_{Int,1}\,g_1$  (a.u./h)', fontsize=13)
    ax.set_xlabel('Time (h)', fontsize=12)
    ax.set_xlim(*tw)
    ax.set_title('(a) the Int2 SOURCE waveform:  the matched arm is scaled down so its\n'
                 'total integral over one mod-8 cycle equals the kept arm', fontsize=13)
    ax.legend(fontsize=10.5, loc='upper right')
    ax.grid(alpha=.22); ax.tick_params(labelsize=11)

    # ---------------- (b) mature Int2 ----------------
    ax = axes[0][1]
    for k, (t_, y_) in (('kept', (tk, yk)), ('removed', (tr, yr)), ('matched', (tm, ym))):
        ax.plot(t_[m], y_[C.VH.IDX['b2_I']][m], color=STYLE[k]['c'],
                ls=STYLE[k]['ls'], lw=2.4, label=STYLE[k]['lab'])
    ax.set_ylabel('mature Int2  (a.u.)', fontsize=13)
    ax.set_xlabel('Time (h)', fontsize=12)
    ax.set_xlim(*tw)
    ax.set_title('(b) the resulting MATURE Int2:  matching the source integral does not\n'
                 'make the three Int2 trajectories coincide', fontsize=13)
    ax.legend(fontsize=10.5, loc='upper right')
    ax.grid(alpha=.22); ax.tick_params(labelsize=11)

    # ---------------- (c) band-edge distance per steady window ----------------
    ax = axes[1][0]
    x = np.arange(len(bk))
    for k, b in (('kept', bk), ('removed', br), ('matched', bm)):
        ax.plot(x, b, color=STYLE[k]['c'], ls=STYLE[k]['ls'], lw=2.4, marker='o',
                ms=5, label='%s   min %.4f' % (STYLE[k]['lab'].split('  ')[0],
                                               np.nanmin(b)))
    ax.set_xlabel('steady read window index (after drop-%d)' % DROP, fontsize=12)
    ax.set_ylabel('bit2 distance to band edge', fontsize=13)
    ax.set_xlim(-0.5, len(bk) - 0.5)
    ax.set_xticks(list(range(0, len(bk), 2)))
    ax.set_title('(c) the discriminating quantity:  matching the total input recovers most\n'
                 'of the lost band margin  (%.4f -> %.4f unmatched -> %.4f matched)'
                 % (np.nanmin(bk), np.nanmin(br), np.nanmin(bm)), fontsize=13)
    ax.legend(fontsize=10.5, loc='lower left')
    ax.grid(alpha=.22); ax.tick_params(labelsize=11)

    # ---------------- (d) far-dwell activity + in-cycle source integral -------
    # NOTE alpha_Int[1] does NOT enter g1 (g1 is built from A1/F1, and the
    # coupling is one-way), so the removed and matched arms have IDENTICAL raw
    # gate far-dwell/dose.  The quantity that actually reaches bit2 is the
    # SOURCE alpha*g1, so panel (d) is drawn in the source domain -- otherwise
    # two of the three arms would coincide and the panel would mislead.
    ax = axes[1][1]
    w = 0.26
    n = min(len(mech[k]['per_cycle']) for k in mech)
    xx = np.arange(n)
    a2 = ax.twinx()
    # set the scale BEFORE the ticks: set_yscale installs a fresh locator and
    # silently discards any set_yticks applied earlier (that is why a 10^-14
    # label kept landing below the canvas)
    ax.set_yscale('log')
    vals = []
    for i, k in enumerate(('kept', 'removed', 'matched')):
        a = {'kept': ALPHA0, 'removed': ALPHA0, 'matched': ALPHA_M}[k]
        far = np.array([c['g1_far_dwell'] for c in mech[k]['per_cycle']][:n]) * a
        dose = np.array([c['g1_dose'] for c in mech[k]['per_cycle']][:n]) * a
        vals += [float(v) for v in far if v > 0]
        ax.bar(xx + (i - 1) * w, far, w, color=STYLE[k]['c'], alpha=.88,
               label='%s  source at far dwell' % STYLE[k]['lab'].split('  ')[0])
        a2.plot(xx, dose, color=STYLE[k]['c'], ls=STYLE[k]['ls'], lw=2.0,
                marker='.', ms=8)
    lo, hi = min(vals), max(vals)
    ax.set_ylim(lo * 0.3, hi * 3.0)
    # explicit decade ticks strictly inside the view
    e0 = math.floor(math.log10(lo * 0.3))
    e1 = math.ceil(math.log10(hi * 3.0))
    ax.set_yticks([10.0 ** e for e in range(e0, e1 + 1)])
    ax.minorticks_off()
    a2.minorticks_off()
    a2.set_ylabel(r'source per cycle  $\alpha_{Int,1}\!\int g_1\,dt$  (a.u.)',
                  fontsize=11.5, color=GRAY)
    a2.tick_params(labelsize=10, colors=GRAY)
    ax.set_xlabel('steady carry cycle index', fontsize=12)
    ax.set_ylabel(r'source at far dwell  $\alpha_{Int,1}\,g_1$', fontsize=12)
    ax.set_title('(d) at MATCHED total input the unmatched-in-time arm still sits orders of\n'
                 'magnitude above KEPT in the far dwell -- the gate is doing temporal screening',
                 fontsize=13)
    ax.legend(fontsize=10, loc='lower right')
    ax.grid(alpha=.22, axis='y'); ax.tick_params(labelsize=11)

    s = json.loads((dm / 'summary.json').read_text(encoding='utf-8'))
    fig.suptitle('Does the Int0 clock gate act by TEMPORAL SCREENING or just by changing the '
                 'total Int2 input?\n300 h, 34-state hybrid, clock_K = 0.30/2;  matching window '
                 '= read windows 8-15 (%.1f-%.1f h), alpha_Int[1] %.6f -> %.6f'
                 % (s['match_window']['t_start_h'], s['match_window']['t_end_h'],
                    ALPHA0, ALPHA_M), fontsize=13.5)
    fig.text(0.5, 0.004,
             'all three arms: certified = True, 19 steady windows, sequence '
             '1234567012345670123;  matching verified to %.0e relative'
             % s['total_source']['rel_gap_matched_vs_kept'],
             ha='center', va='bottom', fontsize=10.5, color=GRAY)

    out = HERE / 'out'
    out.mkdir(exist_ok=True)
    for ext in ('png', 'svg'):
        fig.savefig(out / f'fig_M21_three_arm.{ext}', dpi=150, facecolor='white')
    plt.close(fig)
    print('wrote', out / 'fig_M21_three_arm.png', flush=True)

    (out / 'fig_M21_three_arm_provenance.json').write_text(json.dumps(dict(
        source_2arm=d2.name, source_matched=dm.name,
        arms={k: dict(band_edge_min=float(np.nanmin(b)),
                      band_edge_per_window=[None if np.isnan(v) else float(v) for v in b])
              for k, b in (('kept', bk), ('removed', br), ('matched', bm))},
        alpha_Int1=dict(kept=ALPHA0, removed=ALPHA0, matched=ALPHA_M),
        match_window=s['match_window'], total_source=s['total_source'],
        bit_margins=s['bit_margins'],
        note=('no integration in this script; source/mature traces rebuilt from saved '
              'states.  far-dwell = median of the lowest 20 % of g1 within a cycle '
              '(common.mechanism), dose = in-cycle integral of g1.'),
    ), indent=1, ensure_ascii=False), encoding='utf-8')


if __name__ == '__main__':
    main()
