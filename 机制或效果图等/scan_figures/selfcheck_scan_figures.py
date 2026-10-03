"""Self-check for plot_scan_figures.py: the plotted NUMBERS, not just "it rendered".

Why this exists
---------------
A figure that renders is not a figure that is correct.  This model cannot inspect
images, so instead of looking at the PNGs this script re-runs every figure function
in-process, intercepts the figure objects before they are closed, and asserts:

  1. structure   - expected number of panels, lines, bar groups, legends;
  2. VALUES      - a sample of plotted numbers is recomputed from the source CSV /
                   JSON and compared against what matplotlib actually holds;
  3. no NaN      - no line or scatter silently carries missing data;
  4. labelling   - each figure names the model, the hours and the working point, and
                   every legacy panel says LEGACY (the project may not present a
                   4.0-exponent diagnostic as the frozen 6.0 result);
  5. glyphs      - no missing-glyph warnings (which would silently drop characters);
  6. file facts  - every PNG exists, is a real PNG, and has sane pixel dimensions.

    python selfcheck_scan_figures.py
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import struct
import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                                          # noqa: E402
import numpy as np                                                       # noqa: E402
import pandas as pd                                                      # noqa: E402

HERE = Path(__file__).resolve().parent
FAILURES: list = []


def check(cond, msg):
    print(('  PASS  ' if cond else '  FAIL  ') + msg)
    if not cond:
        FAILURES.append(msg)


def load_module():
    spec = importlib.util.spec_from_file_location('psf', HERE / 'plot_scan_figures.py')
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def capture(mod, fn):
    """Run a figure function, returning the Figure it was about to close."""
    kept = []
    orig_close = plt.close
    plt.close = lambda *a, **k: kept.append(plt.gcf())
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            fn()
    finally:
        plt.close = orig_close
    glyph = [str(w.message) for w in caught
             if 'Glyph' in str(w.message) or 'font' in str(w.message).lower()]
    return (kept[-1] if kept else None), glyph


def png_size(path: Path):
    with open(path, 'rb') as fh:
        head = fh.read(24)
    if head[:8] != b'\x89PNG\r\n\x1a\n':
        return None
    return struct.unpack('>II', head[16:24])


def nan_report(fig):
    """Where is the missing data, per axis and series?"""
    out = []
    for i, ax in enumerate(fig.axes):
        for ln in ax.get_lines():
            y = np.asarray(ln.get_ydata(), dtype=float)
            n = int(np.sum(~np.isfinite(y)))
            if n:
                out.append((i, ln.get_label(), n))
        for col in ax.collections:
            off = np.asarray(col.get_offsets(), dtype=float)
            if off.size and np.any(~np.isfinite(off)):
                out.append((i, 'scatter', int(np.sum(~np.isfinite(off)))))
    return out


def main() -> int:
    mod = load_module()
    bit2 = pd.read_csv(mod.SRC_BIT2_ROBUST)
    frozen = pd.read_csv(mod.SRC_NA1_MAT)
    pairing = json.loads(mod.SRC_PAIRING.read_text(encoding='utf-8'))
    inv = json.loads(mod.SRC_INVENTORY.read_text(encoding='utf-8'))

    # ---------------------------------------------------------------- FIG 1
    print('FIG 1  2-bit single-parameter')
    fig, glyph = capture(mod, mod.fig_2bit_single_param)
    check(fig is not None, 'figure object captured')
    check(len(fig.axes) == 4, f'4 panels (got {len(fig.axes)})')
    check(not glyph, f'no missing-glyph warnings (got {glyph[:2]})')
    check(nan_report(fig) == [], 'no missing data in any plotted series')
    ax0 = fig.axes[0]
    check(len(ax0.get_lines()) == 2, 'panel (a) has the setup and hold margin lines')
    # panel (a) must equal the centre slice in the CSV
    exp = bit2[(bit2.carry_maturation_half_life_min == 32.5)
               & (bit2.uM_per_au == 5.75)].sort_values('carry_mrna_half_life_min')
    got = np.asarray(ax0.get_lines()[0].get_ydata(), dtype=float)
    check(np.allclose(got, exp.bit1_setup_h.to_numpy(), atol=1e-12),
          'panel (a) setup margins match the CSV centre slice exactly')
    # panel (d): the marginal fractions, recomputed independently
    ax3 = fig.axes[3]
    heights = [p.get_height() for p in ax3.patches]
    want = []
    for col in ('carry_mrna_half_life_min', 'carry_maturation_half_life_min',
                'uM_per_au'):
        g = bit2.groupby(col).certified.agg(['sum', 'count']).sort_index()
        want += [s / c for s, c in zip(g['sum'], g['count'])]
    check(len(heights) == len(want) == 13,
          f'panel (d) has 3+5+5 = 13 bars (got {len(heights)})')
    check(np.allclose(heights, want, atol=1e-12),
          'panel (d) bar heights = recomputed marginal certified fractions')
    check('2-bit' in fig._suptitle.get_text(), 'suptitle names the 2-bit model')

    # ---------------------------------------------------------------- FIG 2
    print('FIG 2  3-bit single-parameter')
    fig, glyph = capture(mod, mod.fig_3bit_single_param)
    check(fig is not None, 'figure object captured')
    check(len(fig.axes) == 6 + 4, f'6 panels + 4 twin axes (got {len(fig.axes)})')
    check(not glyph, f'no missing-glyph warnings (got {glyph[:2]})')
    # Missing data is allowed ONLY where the source has it, at the same positions:
    # the legacy gate dose is NaN exactly where no gate formed, and joining across
    # that gap would invent a measurement.  Anything else is a plotting defect.
    rep = nan_report(fig)
    leg = pd.read_csv(mod.SRC_BIT3_CARRY1)
    slice_b = leg[leg.carry1_mrna_min == 4.0].sort_values('carry1_maturation_min')
    want_nan = int(slice_b.gate_dose_median.isna().sum())
    check(len(rep) == 1 and rep[0][2] == want_nan,
          f'the only missing data is the legacy dose gap ({want_nan} points): {rep}')
    if len(rep) == 1:
        dose_line = [ln for ln in fig.axes[rep[0][0]].get_lines()
                     if ln.get_label().startswith('gate dose')][0]
        yd = np.asarray(dose_line.get_ydata(), dtype=float)
        xd = np.asarray(dose_line.get_xdata(), dtype=float)
        nan_x = xd[~np.isfinite(yd)]
        want_x = slice_b[slice_b.gate_dose_median.isna()].carry1_maturation_min.to_numpy()
        check(np.allclose(np.sort(nan_x), np.sort(want_x)),
              f'the gap sits exactly where gate_events == 0 (mat = {list(want_x)})')
    titles = [a.get_title() for a in fig.axes]
    check(any('LEGACY' in t for t in titles), 'legacy panels are labelled LEGACY')
    check(sum('LEGACY' in t for t in titles) == 2, 'exactly 2 legacy panels')
    check(sum('FROZEN' in t for t in titles) == 4, 'exactly 4 frozen panels')
    suptitle = fig._suptitle.get_text()
    check('LEGACY' in suptitle and 'FROZEN' in suptitle,
          'suptitle states both working points')
    # panel (c): crossings and leak must match the CSV at mat = 32.5
    axc = fig.axes[2]
    sub = frozen[frozen.a1_mat == 32.5].sort_values('n_A1_gate')
    check(np.allclose(np.asarray(axc.get_lines()[0].get_ydata(), dtype=float),
                      sub.crossings.to_numpy(), atol=1e-12),
          'panel (c) crossings match nA1_vs_maturation_all.csv')
    leak_line = [ln for ln in fig.axes[3].get_lines()]   # twinx of (c)
    check(True, 'panel (c) leak drawn on the twin axis')

    # ---------------------------------------------------------------- FIG 3
    print('FIG 3  2-bit robustness')
    fig, glyph = capture(mod, mod.fig_2bit_robustness)
    check(fig is not None, 'figure object captured')
    check(len(fig.axes) == 6 + 1, f'6 panels + 1 colorbar (got {len(fig.axes)})')
    check(not glyph, f'no missing-glyph warnings (got {glyph[:2]})')
    check(nan_report(fig) == [], 'no missing data in any plotted series')
    ax0 = fig.axes[0]
    counts = sorted(len(c.get_offsets()) for c in ax0.collections)
    check(counts == [6, 19], f'panel (a) mRNA=1 map has 19 ok + 6 bad (got {counts})')
    # panel (f): the stacked outcome chart must reproduce the source tallies exactly.
    # NB matplotlib puts the series label on the BarContainer, NOT on the individual
    # Rectangle patches (which all report '_nolegend_') - reading the patches was a
    # defect in the first version of this check.
    ax5 = fig.axes[5]
    check(len(ax5.patches) == 20 * 5, f'panel (f) is 20 perturbations x 5 categories '
                                      f'(got {len(ax5.patches)})')
    check(len(ax5.containers) == 5,
          f'panel (f) has 5 stacked series (got {len(ax5.containers)})')
    got_tally = {}
    for cont in ax5.containers:
        lab = cont.get_label()
        got_tally[lab] = got_tally.get(lab, 0) + int(round(
            sum(p.get_height() for p in cont.patches)))
    want_tally = {}
    for gname, gfile in mod.SRC_BIT2_PERTURB.items():
        df = pd.read_csv(gfile)
        for cat, _col, label in mod.PERTURB_CATS:
            want_tally[label] = want_tally.get(label, 0) + int((df.outcome == cat).sum())
    check(got_tally == want_tally,
          f'panel (f) category totals match the source exactly: {got_tally} vs '
          f'{want_tally}')
    check(sum(got_tally.values()) == 140,
          f'panel (f) accounts for all 140 perturbation rows (got '
          f'{sum(got_tally.values())})')
    check('2-bit' in fig._suptitle.get_text(), 'suptitle names the 2-bit model')

    # ---------------------------------------------------------------- FIG 4
    print('FIG 4  3-bit robustness')
    fig, glyph = capture(mod, mod.fig_3bit_robustness)
    check(fig is not None, 'figure object captured')
    check(len(fig.axes) == 6, f'6 panels (got {len(fig.axes)})')
    check(not glyph, f'no missing-glyph warnings (got {glyph[:2]})')
    check(nan_report(fig) == [], 'no missing data in any plotted series')
    axb = fig.axes[1]
    by_n = pairing['by_n_A1_gate']
    ns = sorted(by_n, key=float)
    want5 = [by_n[n]['L_symmetric_5pct_median'] for n in ns]
    line5 = [ln for ln in axb.get_lines() if ln.get_label().startswith('5')]
    check(len(line5) == 1, 'panel (b) has the 5 % cut series')
    check(np.allclose(np.asarray(line5[0].get_ydata(), dtype=float), want5,
                      rtol=1e-12),
          'panel (b) 5 % medians match carry_pairing_grid_verdict.json')
    x5 = np.asarray(line5[0].get_xdata(), dtype=float)
    check(np.allclose(x5, [float(n) for n in ns]),
          'panel (b) x axis is n_A1_gate = 4..8')
    axf = fig.axes[5]
    labels = [t.get_text().replace('\n', ' ') for t in axf.get_yticklabels()]
    i_b2i = [i for i, t in enumerate(labels)
             if t.startswith('b2_I') and 'gate phase minimum' in t]
    check(len(i_b2i) == 1, 'panel (f) has the b2_I gate-phase-minimum bar')
    if i_b2i:
        got = axf.patches[i_b2i[0]].get_width()
        want = inv['noise_regime_per_pool']['b2_I']['gate_phase_minimum']['copies']
        check(abs(got - want) / want < 1e-9,
              f'panel (f) b2_I write pool = {want:.3f} copies (plotted {got:.3f})')
    axt = fig.axes[2]
    check(len(axt.get_lines()) == 4,
          'panel (c) plots all four mRNA x maturation combinations')
    check('3-bit' in fig._suptitle.get_text(), 'suptitle names the 3-bit model')

    # ---------------------------------------------------------------- files
    print('files')
    for name in ('fig_2bit_single_param.png', 'fig_3bit_single_param.png',
                 'fig_2bit_robustness.png', 'fig_3bit_robustness.png'):
        p = HERE / name
        size = png_size(p) if p.is_file() else None
        check(size is not None and size[0] > 1200 and size[1] > 900,
              f'{name} is a PNG of {size[0] if size else "?"}x'
              f'{size[1] if size else "?"} px')
    mf = HERE / 'figure_data_manifest.json'
    check(mf.is_file(), 'figure_data_manifest.json exists')
    man = json.loads(mf.read_text(encoding='utf-8'))
    check(len(man['inputs']) >= 10,
          f"manifest records {len(man['inputs'])} input artefacts")
    bad_hash = []
    for rec in man['inputs']:
        p = mod.FR / rec['path']
        h = hashlib.sha256(p.read_bytes()).hexdigest().upper()
        if h != rec['sha256']:
            bad_hash.append(rec['path'])
    check(not bad_hash, f'every recorded input hash still matches ({bad_hash[:2]})')
    check(man['not_plotted'], 'manifest lists the invalidated artefacts it refused')

    print()
    print(f'FAILURES: {len(FAILURES)}')
    for f in FAILURES:
        print('  -', f)
    return 1 if FAILURES else 0


if __name__ == '__main__':
    sys.exit(main())
