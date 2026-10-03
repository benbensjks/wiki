"""M-21 bit-state comparison: clock gate KEPT / REMOVED / REMOVED-and-matched.

Style follows `wiki_submission_assets/plot_hybrid_submission.py` (the certified
hybrid bit-state figure): same rcParams, same band shading (0-0.30 / 0.70-1) and
0.30/0.70 guide lines, same dpi, same state definitions --

    Bit 0 = b0_S ,  Bit 1 = 1 - pb1_zmh ,  Bit 2 = b2_S
    Int0 translation = 30 * b0_M_I

Panel layout is per-bit (rather than the reference's single combined bit panel)
because that is what makes the result legible: bit0 and bit1 are BIT-IDENTICAL
across the three arms (the shared 23-state prefix), while bit2 is the only one
that differs.  The script ASSERTS that identity rather than merely drawing it.

Three zoom variants are written from the same trajectories:

    fig_M21_bitstates_3arm    full 300 h (context / record)
    fig_M21_bitstates_1p5cyc  1.5 carry periods
    fig_M21_bitstates_1cyc    1   carry period   <- most legible for display

One "carry period" = one bit1->bit2 carry cycle = 2 clock read windows (~21 h),
i.e. the period of the phenomenon under study (bit2 flips once per carry).  The
zoom is anchored on the steady window where the REMOVED arm's band-edge distance
is smallest, so the contrast is at its maximum inside the frame.

Reads saved trajectories only -- no integration, no smoothing, no resampling.
The other two M-21 figures (fig_M21_clock_gate_bypass.*, fig_M21_three_arm.*)
are left untouched.
"""
from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.dont_write_bytecode = True

HERE = Path(__file__).resolve().parent
SCAN = HERE.parent / '时钟工作区扫描'
ROOT = HERE.parent
for _d in (str(SCAN), str(HERE)):
    if _d not in sys.path:
        sys.path.insert(0, _d)

import common as C                                                # noqa: E402
from verify_hzh import IDX                                        # noqa: E402

SRC2 = 'm21_20261001_225706'
INT0_COLOR = '#6853a3'
ARM = dict(kept=dict(c='#4F9194', ls='-', lab='KEPT  (clock present)'),
           removed=dict(c='#D29144', ls='--',
                        lab='REMOVED  (clock = 1, total input +62 %)'),
           matched=dict(c='#6a1b9a', ls='-.',
                        lab='MATCHED  (clock = 1, total Int2 input matched)'))
VARIANTS = [('fig_M21_bitstates_3arm', 0.0),
            ('fig_M21_bitstates_1p5cyc', 1.5),
            ('fig_M21_bitstates_1cyc', 1.0)]
BAND_EDGE = dict(kept=0.2821, removed=0.1857, matched=0.2524)


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest().upper()


def newest(pat: str) -> Path:
    d = sorted([p for p in (HERE / 'results').glob(pat)
                if p.is_dir() and 'SMOKE' not in p.name and 'SUPERSEDED' not in p.name])
    if not d:
        raise SystemExit('no results for ' + pat)
    return d[-1]


def render(name, n_periods, t, st, traj, reads, checks, t0, period):
    """Draw one variant.  n_periods = 0 means the full 300 h."""
    kept_arm = 'kept'
    plt.rcParams.update({'svg.fonttype': 'none', 'font.size': 11})
    fig, ax = plt.subplots(5, 1, figsize=(12, 13.5), sharex=True,
                           constrained_layout=True,
                           gridspec_kw={'height_ratios': [.70, 1.0, 1.0, 1.5, .95]})

    a = ax[0]
    for k in ('kept', 'removed', 'matched'):
        a.plot(t, 30.0 * traj[k][1][IDX['b0_M_I']], color=INT0_COLOR,
               ls=ARM[k]['ls'], lw=1.3, alpha=.9)
    a.set_ylabel('Int0 translation\n(a.u./h)', fontsize=10)
    a.set_title('(a) shared driver — the three arms are bit-identical here '
                '(same 23-state prefix, one-way coupling)', fontsize=11)

    for i, (b, nm) in enumerate((('S0', 'Bit 0  (biochemical, b0_S)'),
                                 ('S1', 'Bit 1  (reduced, 1 − pb1_zmh)')), start=1):
        a = ax[i]
        a.axhspan(0, 0.3, color='#329b78', alpha=.06)
        a.axhspan(0.7, 1, color='#2878b5', alpha=.06)
        a.axhline(.3, color='gray', ls='--', lw=.7)
        a.axhline(.7, color='gray', ls='--', lw=.7)
        for k in ('kept', 'removed', 'matched'):
            a.plot(t, st[k][b], color=ARM[k]['c'], ls=ARM[k]['ls'], lw=1.5,
                   label=ARM[k]['lab'] if i == 1 else None)
        a.set_ylabel(b + '  (LR fraction)', fontsize=10)
        a.set_ylim(-.04, 1.13)
        a.set_title('(%s) %s — identical in all three arms (max|Δ| = %.0e)'
                    % ('bc'[i - 1], nm, checks[b + '_kept_vs_removed_maxabs']),
                    fontsize=11)
        if i == 1:
            a.legend(loc='upper center', ncol=3, fontsize=8.5)

    a = ax[3]
    a.axhspan(0, 0.3, color='#329b78', alpha=.06)
    a.axhspan(0.7, 1, color='#2878b5', alpha=.06)
    a.axhline(.3, color='gray', ls='--', lw=.7)
    a.axhline(.7, color='gray', ls='--', lw=.7)
    for k in ('kept', 'removed', 'matched'):
        a.plot(t, st[k]['S2'], color=ARM[k]['c'], ls=ARM[k]['ls'], lw=1.5)
    for r in reads[int(C.VH.DROP):]:
        a.axvspan(r['start_h'], r['end_h'], color='#303030', alpha=.05, lw=0)
    a.set_ylabel('Bit 2  S2  (LR fraction)', fontsize=10)
    a.set_ylim(-.04, 1.13)
    a.set_title('(d) Bit 2 — the ONLY bit that differs:  max|removed − kept| = %.3f, '
                'max|matched − kept| = %.3f\nshaded = steady read windows;  '
                'band-edge distance %.4f (kept) / %.4f (removed) / %.4f (matched)'
                % (checks['S2_kept_vs_removed_maxabs'],
                   checks['S2_kept_vs_matched_maxabs'],
                   BAND_EDGE['kept'], BAND_EDGE['removed'], BAND_EDGE['matched']),
                fontsize=11)

    a = ax[4]
    a.step(reads_t, reads_v, where='post', color='#303030', lw=1.4)
    a.scatter(reads_t, reads_v, c=reads_v, cmap='viridis', s=27, zorder=3)
    a.set_yticks(range(8))
    a.set_ylabel('Decoded value', fontsize=10)
    a.set_xlabel('Time (h)', fontsize=11)
    a.set_title('(e) decoded value — all three arms read the same mod-8 sequence '
                '1234567012345670123 (19 steady windows after drop-%d)'
                % int(C.VH.DROP), fontsize=11)

    for a in ax:
        a.grid(alpha=.12)
    if n_periods == 0:
        for a in ax:
            a.set_xlim(0, 300)
        span_l1 = 'full 300 h run'
        span_l2 = ''
    else:
        x1 = t0 + n_periods * period
        for a in ax:
            a.set_xlim(t0, x1)
        # explicit ticks strictly inside the view: the auto locator puts a tick
        # just outside the limit whose label is then half off the canvas
        start = math.ceil(t0 / 5.0) * 5.0
        ticks = [v for v in np.arange(start, x1 - 1e-9, 5.0)]
        ax[-1].set_xticks(ticks)
        span_l1 = ('%g carry period%s = %.1f h   (%.1f – %.1f h)'
                   % (n_periods, '' if n_periods == 1 else 's', x1 - t0, t0, x1))
        span_l2 = 'one carry period = one bit1→bit2 carry cycle = 2 clock read windows'
    fig.suptitle('Int0 clock gate kept vs removed vs input-matched:  a BIT-STATE comparison\n'
                 '34-state hybrid HZH_MOD8_CAUSAL_V1, clock_K = 0.30/2 — '
                 + span_l1 + ('\n' + span_l2 if span_l2 else ''),
                 fontsize=12.5)
    out = HERE / 'out'
    out.mkdir(exist_ok=True)
    for ext in ('png', 'svg'):
        fig.savefig(out / f'{name}.{ext}', dpi=180, facecolor='white')
    plt.close(fig)
    print('wrote', out / f'{name}.png', flush=True)
    return span_l1


def main():
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8')
    d2 = HERE / 'results' / SRC2
    dm = newest('m21matched_*')
    print('2-arm:', d2.name, ' matched:', dm.name, flush=True)

    src = dict(kept=d2 / 'traj_kept_K0.30.npz',
               removed=d2 / 'traj_removed_K0.30.npz',
               matched=dm / 'traj_matched_K0.30.npz')
    traj, meta = {}, {}
    for k, p in src.items():
        z = np.load(p)
        traj[k] = (z['time_h'], z['states'])
        meta[k] = dict(file=str(p.relative_to(ROOT)), sha256=sha(p))
    t = traj['kept'][0]
    st = {k: dict(S0=v[1][IDX['b0_S']], S1=1.0 - v[1][IDX['pb1_zmh']],
                  S2=v[1][IDX['b2_S']]) for k, v in traj.items()}

    checks = {}
    for b in ('S0', 'S1', 'S2'):
        checks[b + '_kept_vs_removed_maxabs'] = float(
            np.max(np.abs(st['removed'][b] - st['kept'][b])))
        checks[b + '_kept_vs_matched_maxabs'] = float(
            np.max(np.abs(st['matched'][b] - st['kept'][b])))
    identity_ok = (checks['S0_kept_vs_removed_maxabs'] == 0.0
                   and checks['S1_kept_vs_removed_maxabs'] == 0.0)
    print('    bit0/bit1 identical across arms:', identity_ok, flush=True)
    print('    bit2 max|removed-kept| = %.6f ; max|matched-kept| = %.6f'
          % (checks['S2_kept_vs_removed_maxabs'],
             checks['S2_kept_vs_matched_maxabs']), flush=True)
    if not identity_ok:
        raise SystemExit('bit0/bit1 are NOT identical across arms -- this figure '
                         'would be making a false claim; aborting')

    # ---- read windows (shared: they come from the shared Int0) ----
    han = C.HanInput(300.0)
    mk = C.HZHModel('han', han)
    sk = C.sig_of(mk, traj['kept'][1], 'hzh')
    _, reads = C.VH.read_windows(t, sk)
    global reads_t, reads_v
    reads_t = np.array([r['trough_h'] for r in reads])
    reads_v = np.array([np.nan if r['value'] is None else r['value'] for r in reads])

    # ---- anchor: the steady window with the SMALLEST removed-arm band margin ----
    DROP = int(C.VH.DROP)
    sel = reads[DROP:]
    per_win = []
    S2r = st['removed']['S2']
    for r in sel:
        m = (t >= r['start_h']) & (t <= r['end_h'])
        v = S2r[m]
        lab = r['labels'][2]
        per_win.append(C.VH.LOW - float(v.max()) if lab == 0
                       else float(v.min()) - C.VH.HIGH if lab == 1 else np.nan)
    i_anchor = int(np.nanargmin(per_win))
    clock_period = float(np.median(np.diff([r['cycle_start_h'] for r in reads])))
    period = 2.0 * clock_period          # one carry period = 2 clock windows
    t0 = float(sel[i_anchor]['cycle_start_h'])
    print('    anchor: steady[%d] cycle %d  t0 = %.3f h   clock period %.4f h  '
          'carry period %.4f h  (removed margin there %.4f)'
          % (i_anchor, sel[i_anchor]['cycle'], t0, clock_period, period,
             per_win[i_anchor]), flush=True)

    spans = {}
    for name, nper in VARIANTS:
        spans[name] = render(name, nper, t, st, traj, reads, checks, t0, period)

    out = HERE / 'out'
    (out / 'fig_M21_bitstates_provenance.json').write_text(json.dumps(dict(
        style_reference='wiki_submission_assets/plot_hybrid_submission.py',
        inputs=meta, hours=float(t[-1]), n_samples=int(len(t)),
        state_definitions=dict(bit0='b0_S', bit1='1 - pb1_zmh', bit2='b2_S',
                               int0_translation='30 * b0_M_I'),
        identity_checks=checks,
        identity_asserted=bool(identity_ok),
        read_count=int(np.sum(~np.isnan(reads_v))),
        sequence=''.join('x' if np.isnan(v) else str(int(v)) for v in reads_v),
        band_edge_distance=BAND_EDGE,
        zoom=dict(definition='one carry period = one bit1->bit2 carry cycle = '
                             '2 clock read windows',
                  clock_period_h=clock_period, carry_period_h=period,
                  anchor_steady_index=i_anchor,
                  anchor_cycle=int(sel[i_anchor]['cycle']),
                  anchor_t0_h=t0,
                  anchor_removed_margin=float(per_win[i_anchor]),
                  variants=spans),
        note=('saved trajectories only; no integration, smoothing or resampling.  '
              'bit0/bit1 identity across arms is ASSERTED in this script, not drawn.'),
        outputs={p.name: sha(p) for p in sorted(out.glob('fig_M21_bitstates_*'))
                 if p.suffix in ('.png', '.svg')},
    ), ensure_ascii=False, indent=2), encoding='utf-8')
    print('wrote provenance', flush=True)


if __name__ == '__main__':
    main()
