"""M-21 pair figure: HZH clock gate KEPT vs REMOVED, one complete mod-8 cycle.

Requested layout
----------------
* bit0 / bit1 / bit2 go into ONE panel (colour = bit), not one panel per bit;
* each arm gets exactly ONE figure (3 rows: Int0 translation / three bit states /
  decoded value);
* the two arms are then rendered side by side as a single spliced figure;
* the window is one COMPLETE mod-8 cycle, so the decoded staircase visits all
  eight states 0..7.

Window
------
steady read windows 7..14 = carry cycles 15..22 = 161.933 - 246.667 h (84.733 h,
8 read windows), decoded values 0,1,2,3,4,5,6,7.  This is the first complete
steady mod-8 cycle that starts on value 0, and it also contains the stretch
where the REMOVED arm's band margin is smallest (cycles 21-22), so the contrast
is at its maximum inside the frame.

Style follows `wiki_submission_assets/plot_hybrid_submission.py`: same rcParams,
band shading (0-0.30 / 0.70-1), 0.30/0.70 guide lines, dpi=180, and the same
state definitions (Bit 0 = b0_S, Bit 1 = 1 - pb1_zmh, Bit 2 = b2_S,
Int0 translation = 30 * b0_M_I).

Reads saved trajectories only -- no integration, no smoothing, no resampling.
All previously written figures are left untouched.
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
BIT_C = ['#2878b5', '#dc8b22', '#329b78']
BIT_L = ['Bit 0  (b0_S)', 'Bit 1  (1 − pb1_zmh)', 'Bit 2  (b2_S)']
INT0_C = '#6853a3'
ARM_C = dict(kept='#4F9194', removed='#D29144')
ARM_T = dict(kept='KEPT      clock = H(Int0; 0.30, 2.0)',
             removed='REMOVED      clock ≡ 1')
BAND_EDGE = dict(kept=0.2821, removed=0.1857)
WIN = (161.933, 246.667)


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest().upper()


def load():
    d2 = HERE / 'results' / SRC2
    src = dict(kept=d2 / 'traj_kept_K0.30.npz', removed=d2 / 'traj_removed_K0.30.npz')
    traj, meta = {}, {}
    for k, p in src.items():
        z = np.load(p)
        traj[k] = (z['time_h'], z['states'])
        meta[k] = dict(file=str(p.relative_to(ROOT)), sha256=sha(p))
    han = C.HanInput(300.0)
    mk = C.HZHModel('han', han)
    sk = C.sig_of(mk, traj['kept'][1], 'hzh')
    _, reads = C.VH.read_windows(traj['kept'][0], sk)
    return traj, meta, reads


def arm_panels(fig, axes, key, traj, reads, show_ylabels, t, st):
    """Draw the three rows of ONE arm into the given axes column."""
    x = (t >= WIN[0]) & (t <= WIN[1])
    a0, a1, a2 = axes

    # ---- row 1: Int0 translation ----
    a0.plot(t[x], 30.0 * traj[key][1][IDX['b0_M_I']][x], color=INT0_C, lw=1.4)
    a0.set_ylim(0, None)
    if show_ylabels:
        a0.set_ylabel('Int0 translation\n(a.u./h)', fontsize=10)
    a0.set_title(ARM_T[key], fontsize=12.5, color=ARM_C[key], fontweight='bold')

    # ---- row 2: all three bit states in ONE panel ----
    a1.axhspan(0, 0.3, color='#329b78', alpha=.06)
    a1.axhspan(0.7, 1, color='#2878b5', alpha=.06)
    a1.axhline(.3, color='gray', ls='--', lw=.7)
    a1.axhline(.7, color='gray', ls='--', lw=.7)
    for i, b in enumerate(('S0', 'S1', 'S2')):
        a1.plot(t[x], st[key][b][x], color=BIT_C[i], lw=1.5, label=BIT_L[i])
    for r in reads[int(C.VH.DROP):]:
        a1.axvspan(r['start_h'], r['end_h'], color='#303030', alpha=.055, lw=0)
    a1.set_ylim(-.04, 1.13)
    if show_ylabels:
        a1.set_ylabel('DNA LR fraction', fontsize=10)
    a1.legend(loc='lower right', ncol=3, fontsize=8.5, framealpha=.9)

    # ---- row 3: decoded value ----
    rt = np.array([r['trough_h'] for r in reads])
    rv = np.array([np.nan if r['value'] is None else r['value'] for r in reads])
    inwin = (rt >= WIN[0]) & (rt <= WIN[1])
    a2.step(rt[inwin], rv[inwin], where='post', color='#303030', lw=1.5)
    a2.scatter(rt[inwin], rv[inwin], c=rv[inwin], cmap='viridis', s=45, zorder=3,
               edgecolor='k', linewidth=.4)
    a2.set_yticks(range(8))
    a2.set_ylim(-.4, 7.4)
    if show_ylabels:
        a2.set_ylabel('Decoded value', fontsize=10)
    a2.set_xlabel('Time (h)', fontsize=11)


def finish(fig, axes_cols, twocol):
    """Common cosmetics: grid, x ticks inside range, one shared x label."""
    for col in axes_cols:
        for a in col:
            a.grid(alpha=.12)
            a.set_xlim(*WIN)
    start = math.ceil(WIN[0] / 10.0) * 10.0
    ticks = [v for v in np.arange(start, WIN[1] - 1e-9, 10.0)]
    for col in axes_cols:
        col[-1].set_xticks(ticks)
    fig.suptitle('HZH clock gate KEPT vs REMOVED — one complete mod-8 cycle '
                 '(8 read windows = 8 states 0…7)\n'
                 '34-state hybrid HZH_MOD8_CAUSAL_V1, clock_K = 0.30/2;  '
                 't = %.1f – %.1f h' % (WIN[0], WIN[1]), fontsize=12.5)


def main():
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8')
    traj, meta, reads = load()
    t = traj['kept'][0]
    st = {k: dict(S0=v[1][IDX['b0_S']], S1=1.0 - v[1][IDX['pb1_zmh']],
                  S2=v[1][IDX['b2_S']]) for k, v in traj.items()}
    # the claim this pair rests on: bit0/bit1 identical, only bit2 differs
    d0 = float(np.max(np.abs(st['removed']['S0'] - st['kept']['S0'])))
    d1 = float(np.max(np.abs(st['removed']['S1'] - st['kept']['S1'])))
    d2 = float(np.max(np.abs(st['removed']['S2'] - st['kept']['S2'])))
    print('    identity: max|dS0|=%.0e  max|dS1|=%.0e  max|dS2|=%.6f' % (d0, d1, d2),
          flush=True)
    if d0 != 0.0 or d1 != 0.0:
        raise SystemExit('bit0/bit1 differ across arms -- layout assumes they do not')

    plt.rcParams.update({'svg.fonttype': 'none', 'font.size': 11})
    out = HERE / 'out'
    out.mkdir(exist_ok=True)
    written = []

    # ---------- per-arm figures (each arm gets exactly ONE figure) ----------
    for key in ('kept', 'removed'):
        fig, axes = plt.subplots(3, 1, figsize=(9.5, 11), sharex=True,
                                 constrained_layout=True,
                                 gridspec_kw={'height_ratios': [.7, 1.5, 1.0]})
        arm_panels(fig, axes, key, traj, reads, True, t, st)
        finish(fig, [axes], False)
        name = f'fig_M21_pair_{key}'
        for ext in ('png', 'svg'):
            fig.savefig(out / f'{name}.{ext}', dpi=180, facecolor='white')
        plt.close(fig)
        written.append(name)
        print('wrote', out / f'{name}.png', flush=True)

    # ---------- spliced figure: the two arms side by side ----------
    fig, grid = plt.subplots(3, 2, figsize=(18.5, 11), sharex=True,
                             constrained_layout=True,
                             gridspec_kw={'height_ratios': [.7, 1.5, 1.0]})
    arm_panels(fig, grid[:, 0], 'kept', traj, reads, True, t, st)
    arm_panels(fig, grid[:, 1], 'removed', traj, reads, True, t, st)
    # make the two columns directly comparable: identical y limits
    for r in range(3):
        lo = min(grid[r][0].get_ylim()[0], grid[r][1].get_ylim()[0])
        hi = max(grid[r][0].get_ylim()[1], grid[r][1].get_ylim()[1])
        grid[r][0].set_ylim(lo, hi)
        grid[r][1].set_ylim(lo, hi)
    finish(fig, [grid[:, 0], grid[:, 1]], True)

    # short call-outs in the REMOVED column (the only place numbers differ)
    g1 = grid[1][1]
    g1.annotate('bit2 settles shallower:\nmin band-edge distance %.4f -> %.4f'
                % (BAND_EDGE['kept'], BAND_EDGE['removed']),
                xy=(243.0, 0.90), xytext=(196.0, 0.60), fontsize=10.5,
                color='#c62828', fontweight='bold',
                arrowprops=dict(arrowstyle='->', color='#c62828', lw=1.6))
    grid[2][0].text(.02, .96, 'bit0 / bit1 identical in both columns\n(max|Δ| = 0)',
                    transform=grid[2][0].transAxes, fontsize=10, va='top',
                    color='#545454')
    name = 'fig_M21_pair_combined'
    for ext in ('png', 'svg'):
        fig.savefig(out / f'{name}.{ext}', dpi=180, facecolor='white')
    plt.close(fig)
    written.append(name)
    print('wrote', out / f'{name}.png', flush=True)

    (out / 'fig_M21_pair_provenance.json').write_text(json.dumps(dict(
        style_reference='wiki_submission_assets/plot_hybrid_submission.py',
        inputs=meta,
        window=dict(read_windows='steady 7..14', carry_cycles='15..22',
                    t_start_h=WIN[0], t_end_h=WIN[1], span_h=WIN[1] - WIN[0],
                    n_read_windows=8,
                    decoded_values=[0, 1, 2, 3, 4, 5, 6, 7],
                    why=('first complete steady mod-8 cycle that starts on value 0; '
                         'also contains the stretch where the REMOVED arm band '
                         'margin is smallest (cycles 21-22)')),
        state_definitions=dict(bit0='b0_S', bit1='1 - pb1_zmh', bit2='b2_S',
                               int0_translation='30 * b0_M_I'),
        identity_checks=dict(bit0_maxabs=d0, bit1_maxabs=d1, bit2_maxabs=d2),
        band_edge_distance=BAND_EDGE,
        note=('saved trajectories only; no integration, smoothing or resampling.  '
              'The two columns share identical y limits so they are directly '
              'comparable.  bit0/bit1 identity is ASSERTED before drawing.'),
        outputs={p.name: sha(p) for p in sorted(out.glob('fig_M21_pair_*'))
                 if p.suffix in ('.png', '.svg')},
    ), ensure_ascii=False, indent=2), encoding='utf-8')
    print('wrote provenance', flush=True)


if __name__ == '__main__':
    main()
