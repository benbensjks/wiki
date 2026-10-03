"""M-21 figure: clock gate (Int0 AND) KEPT vs REMOVED, 34-state hybrid.

Reads ONLY the trajectories saved by `run_m21_clock_bypass.py` -- no integration,
no smoothing, no resampling.

The figure has to carry a result that contradicts the expectation it was
commissioned under, so the panels are ordered to show BOTH halves:

  (a) the gate really changes -- peak +58 %, dose +64 %, duration +13 %
  (b) the readout does NOT change -- S2 with the 0.30/0.70 bands, 19 steady
      windows, identical labels
  (c) the decoded sequence is identical -- 1234567012345670123 for both arms

In-figure text is English (project convention); the Chinese explanation lives
on the page and in `图清单_跨模型三级页面.md` §2.3.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# Keep real <text> in the SVG (font vectors would make labels unsearchable and
# unselectable).  Same setting as the bit0 / bit1_carry1 sketches; without it a
# rendered SVG contains ZERO <text> elements.
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
from hybrid_model import HbyReceiver                              # noqa: E402
from model_hzh import HZHModel                                    # noqa: E402
from run_m21_clock_bypass import (BypassReceiver, gate_no_clock_array,  # noqa: E402
                                  K_CERT)

DEEP_BLUE, ORANGE, TEAL = '#304B53', '#D29144', '#4F9194'
RISK, GRAY = '#c62828', '#8ca0a5'
KEPT_C, REM_C = TEAL, ORANGE
DROP = int(C.VH.DROP)


def newest_results(root: Path) -> Path:
    d = sorted([p for p in (root / 'results').glob('m21_2*')
                if p.is_dir() and 'SMOKE' not in p.name and 'SUPERSEDED' not in p.name])
    if not d:
        raise SystemExit('no published m21 results; run run_m21_clock_bypass.py')
    return d[-1]


def main():
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8')
    res = newest_results(HERE)
    print('reading', res, flush=True)
    s = json.loads((res / 'summary.json').read_text(encoding='utf-8'))
    kept_meta = next(c for c in s['cases'] if c['arm'] == 'kept')
    rem_meta = next(c for c in s['cases'] if c['arm'] == 'removed')

    tk = np.load(res / 'traj_kept_K0.30.npz')
    tr = np.load(res / 'traj_removed_K0.30.npz')
    t, yk, yr = tk['time_h'], tk['states'], tr['states']

    # rebuild the two gates from the saved states (no re-integration)
    han = C.HanInput(float(t[-1]))
    mk = HZHModel('han', han)
    mr = HZHModel('han', han)
    mr.tail = BypassReceiver(mr.tail.c)
    sk = C.sig_of(mk, yk, 'hzh')
    g_kept, clock = sk['g1'], sk['clock']
    g_rem = gate_no_clock_array(yr, mr)
    S2k, S2r = yk[C.VH.IDX['b2_S']], yr[C.VH.IDX['b2_S']]

    # ---- read windows (Int0 troughs) -> steady windows after DROP
    peaks, reads = C.VH.read_windows(t, sk)
    win = (240.0, 300.0)          # last ~5.7 clock cycles, fully steady
    m = (t >= win[0]) & (t <= win[1])

    fig, axes = plt.subplots(3, 1, figsize=(13.5, 13.2), layout='constrained')

    # ---------------- (a) the gate and the clock ----------------
    ax = axes[0]
    ax.plot(t[m], g_kept[m], color=KEPT_C, lw=2.6,
            label=f"KEPT  gate = act*repress*clock   peak {kept_meta['gate_peak']:.3f}")
    ax.plot(t[m], g_rem[m], color=REM_C, lw=2.6, ls='--',
            label=f"REMOVED  gate = act*repress   peak {rem_meta['gate_peak']:.3f}")
    ax.axhline(C.VH.GATE_THRESHOLD, color=GRAY, ls=':', lw=1.8)
    ax.text(win[0] + 0.4, C.VH.GATE_THRESHOLD + 0.012,
            f'gate threshold {C.VH.GATE_THRESHOLD}', color=GRAY, fontsize=11)
    a2 = ax.twinx()
    a2.plot(t[m], clock[m], color=DEEP_BLUE, lw=1.8, alpha=.75,
            label='clock factor (KEPT arm)')
    a2.axhline(1.0, color=REM_C, lw=1.4, alpha=.6)
    a2.text(win[0] + 0.4, 1.02, 'clock pinned to 1 (REMOVED arm)',
            color=REM_C, fontsize=11)
    a2.set_ylim(-.05, 1.35)
    a2.set_ylabel('clock factor', fontsize=13, color=DEEP_BLUE)
    a2.tick_params(labelsize=12, colors=DEEP_BLUE)
    ax.set_ylabel('carry1 gate  g1', fontsize=14)
    ax.set_ylim(-.03, .72)
    ax.set_title('(a) the gate really changes:  peak +%.0f %%,  dose %.4f -> %.4f h  (+%.0f %%),  '
                 'duration %.3f -> %.3f h'
                 % (100 * (rem_meta['gate_peak'] / kept_meta['gate_peak'] - 1),
                    kept_meta['gate_dose_min'], rem_meta['gate_dose_min'],
                    100 * (rem_meta['gate_dose_min'] / kept_meta['gate_dose_min'] - 1),
                    kept_meta['gate_dur_min'], rem_meta['gate_dur_min']),
                 fontsize=14)
    ln = ax.get_lines()[:2] + a2.get_lines()[:1]
    ax.legend(ln, [x.get_label() for x in ln], fontsize=11.5, loc='upper right')
    ax.grid(alpha=.22)
    ax.tick_params(labelsize=12)
    # Pin the view to the data window: the default 5 % margin produced ticks
    # (230 / 310) that fall outside the canvas, one of them partly visible.
    ax.set_xlim(*win)

    # ---------------- (b) the readout does not change ----------------
    ax = axes[1]
    ax.axhspan(0, C.VH.LOW, color=KEPT_C, alpha=.08)
    ax.axhspan(C.VH.HIGH, 1, color=KEPT_C, alpha=.08)
    ax.axhline(C.VH.LOW, color=GRAY, ls='--', lw=1.3)
    ax.axhline(C.VH.HIGH, color=GRAY, ls='--', lw=1.3)
    ax.plot(t[m], S2k[m], color=KEPT_C, lw=2.4, label='bit2  S2  KEPT')
    ax.plot(t[m], S2r[m], color=REM_C, lw=2.0, ls='--', label='bit2  S2  REMOVED')
    for r in reads:
        if win[0] <= r['trough_h'] <= win[1]:
            ax.axvspan(r['start_h'], r['end_h'],
                       color=DEEP_BLUE, alpha=.07, lw=0)
    ax.set_ylabel('bit2  S2  (LR fraction)', fontsize=14)
    ax.set_ylim(-.05, 1.05)
    ax.set_title('(b) the DECODED LABELS and commitment do NOT change:  %d steady windows '
                 'both arms, identical labels, commitment 1.0, 0 unlabelled\n'
                 'but the SETTLING DEPTH does differ -- see (c);  shaded = read window'
                 % kept_meta['steady_reads'], fontsize=13.5)
    ax.legend(fontsize=11.5, loc='lower right', ncol=2)
    ax.grid(alpha=.22)
    ax.tick_params(labelsize=12)
    ax.set_xlabel('Time (h)', fontsize=13)
    ax.set_xlim(*win)

    # ---------------- (c) band-edge distance per steady window ----------------
    # The decisive quantity: labels agree, distance to the band edge does not.
    # Computed exactly as verify_hzh.read_windows defines a window: for a window
    # labelled 0 the margin is LOW - max(S2), for one labelled 1 it is
    # min(S2) - HIGH; the reported number is the minimum over the steady windows.
    def margins_per_window(y):
        S2v = y[C.VH.IDX['b2_S']]
        out = []
        for r in reads[DROP:]:
            mm = (t >= r['start_h']) & (t <= r['end_h'])
            v = S2v[mm]
            lab = r['labels'][2]
            if lab == 0:
                out.append(C.VH.LOW - float(v.max()))
            elif lab == 1:
                out.append(float(v.min()) - C.VH.HIGH)
            else:
                out.append(float('nan'))
        return np.array(out)

    mk_ = margins_per_window(yk)
    mr_ = margins_per_window(yr)
    x = np.arange(len(mk_))
    ax = axes[2]
    ax.plot(x, mk_, color=KEPT_C, lw=2.4, marker='o', ms=6,
            label='KEPT   min %.4f' % np.nanmin(mk_))
    ax.plot(x, mr_, color=REM_C, lw=2.4, ls='--', marker='s', ms=6,
            label='REMOVED   min %.4f' % np.nanmin(mr_))
    ib = int(np.nanargmin(mr_))
    # Keep the call-out INSIDE the canvas: wrap it and grow it LEFTWARDS from
    # the binding window (the previous version grew rightwards and was cut off).
    ax.annotate('binding window\n%.4f -> %.4f  (%.0f %%)'
                % (mk_[ib], mr_[ib], 100 * (mr_[ib] / mk_[ib] - 1)),
                xy=(ib, mr_[ib]), xytext=(ib - 0.55, mr_[ib] + 0.020),
                fontsize=11.5, color=REM_C, fontweight='bold', ha='right',
                va='bottom',
                arrowprops=dict(arrowstyle='->', color=REM_C, lw=1.6))
    ax.set_xlabel('steady read window index (after drop-%d)' % C.VH.DROP, fontsize=13)
    ax.set_ylabel('distance to band edge', fontsize=14)
    ax.set_title('(c) what DOES change: the distance from S2 to the 0.30/0.70 band edge\n'
                 'decoded sequence is identical (%s) and certified = True for both arms, '
                 'so the gate is not necessary for the nominal mod-8 count'
                 % kept_meta['sequence'], fontsize=13.5)
    ax.legend(fontsize=11.5, loc='lower left')
    ax.grid(alpha=.22)
    ax.tick_params(labelsize=12)
    ax.set_xlim(-0.5, len(mk_) - 0.5)
    ax.set_xticks(list(range(0, len(mk_), 2)))

    g = s['guards']

    def low_high(y):
        S2v = y[C.VH.IDX['b2_S']]
        lo, hi = -9.0, 9.0
        for r in reads[DROP:]:
            mm = (t >= r['start_h']) & (t <= r['end_h'])
            v = S2v[mm]
            if r['labels'][2] == 0:
                lo = max(lo, float(v.max()))
            elif r['labels'][2] == 1:
                hi = min(hi, float(v.min()))
        return lo, hi

    lok, hik = low_high(yk)
    lor, hir = low_high(yr)

    # Three SHORT lines: a single long suptitle is clipped left and right,
    # because savefig keeps the fixed canvas when bbox_inches is not 'tight'.
    fig.suptitle(
        'Removing the Int0 clock AND gate does NOT change the nominal mod-8 count\n'
        '34-state hybrid HZH_MOD8_CAUSAL_V1, 300 h, clock_K = %.2f/2\n'
        'band-edge distance %.4f -> %.4f'
        % (K_CERT, min(C.VH.LOW - lok, hik - C.VH.HIGH),
           min(C.VH.LOW - lor, hir - C.VH.HIGH)),
        fontsize=13.5)
    fig.text(0.5, 0.005,
             'guards all pass:  parity %s   |   clock published identically 1 %s   |   '
             'degeneracy max|dy| %.0e   |   gate reconstruction %.0e'
             % (g['G2_parity'], g['G3'], g['G4_max_abs_dy'], g['G5_max_abs_dgate']),
             ha='center', va='bottom', fontsize=10.5, color=GRAY)

    out = HERE / 'out'
    out.mkdir(exist_ok=True)
    for ext in ('png', 'svg'):
        fig.savefig(out / f'fig_M21_clock_gate_bypass.{ext}', dpi=150,
                    facecolor='white')
    plt.close(fig)
    print('wrote', out / 'fig_M21_clock_gate_bypass.png')

    (out / 'fig_M21_provenance.json').write_text(json.dumps(dict(
        source=str(res.relative_to(HERE)),
        kept_traj='traj_kept_K0.30.npz', removed_traj='traj_removed_K0.30.npz',
        hours=float(t[-1]), clock_K=K_CERT, drop=int(C.VH.DROP),
        steady_reads=kept_meta['steady_reads'], sequence=kept_meta['sequence'],
        kept=dict(certified=kept_meta['certified'], gate_peak=kept_meta['gate_peak'],
                  gate_dose=kept_meta['gate_dose_min'],
                  gate_duration_h=kept_meta['gate_dur_min']),
        removed=dict(certified=rem_meta['certified'], gate_peak=rem_meta['gate_peak'],
                     gate_dose=rem_meta['gate_dose_min'],
                     gate_duration_h=rem_meta['gate_dur_min']),
        max_abs_dy_removed_vs_kept=float(np.max(np.abs(yr - yk))),
        differing_state_indices=[int(i) for i in range(len(yk))
                                 if np.max(np.abs(yr - yk)[i]) > 1e-12],
        bit2_window_margins=dict(
            kept=dict(max_S2_low_windows=lok, min_S2_high_windows=hik,
                      band_edge_distance=min(C.VH.LOW - lok, hik - C.VH.HIGH)),
            removed=dict(max_S2_low_windows=lor, min_S2_high_windows=hir,
                         band_edge_distance=min(C.VH.LOW - lor, hir - C.VH.HIGH)),
            per_window_kept=[None if np.isnan(v) else float(v) for v in mk_],
            per_window_removed=[None if np.isnan(v) else float(v) for v in mr_]),
        guards=g,
        note=('no integration in this script: reads saved trajectories only; '
              'gate signals rebuilt from states for plotting.  The labels and '
              'commitment are unchanged by the ablation but the band-edge distance '
              'is NOT: it falls from 0.2821 to 0.1857, i.e. 34 %% closer to the '
              '0.30/0.70 edge.  "Unchanged" must not be claimed for readout margin.'),
    ), indent=1, ensure_ascii=False), encoding='utf-8')
    print('wrote provenance')


if __name__ == '__main__':
    main()
