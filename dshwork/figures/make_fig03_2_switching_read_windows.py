r"""Figure 3-2 -- Switching and read windows   (dshwork, per the 精简交付版 brief)

Brief: `dshwork/精简交付版/02_作图要求精简版.md` section "Figure 3-2 -- Switching and
read windows".  Requirements:

  * two CONSECUTIVE clock cycles from the LATER part of the selected model's run,
    on ONE shared time axis, three stacked panels:
        top    : Int input and the S0 trajectory, with the 0-0.3 and 0.7-1 read
                 bands shaded, and the read window on each trough;
        middle : RDF0 and the complex C0 -- different magnitudes, so either separate
                 axes or an explicit statement that they are normalised;
        bottom : the ACTUAL forward and reverse recombination fluxes J_fwd0, J_rev0;
  * the read window is 20 % of the clock period in TOTAL (trough +-10 %);
  * REAL trajectory only -- never an idealised square wave;
  * a trailing window with insufficient data must NOT be padded or extended;
  * labels, verbatim: `Time (h)`, `LR fraction S0`, `Read window`, `Low band`,
    `High band`, `RDF0`, `Complex C0`, `Forward flux`, `Reverse flux`.

Data
----
The frozen two-bit (34-state) run: `twobit34_results/trajectories.csv` (300 h, 2 min
grid) and `twobit34_results/read_windows.csv`, which holds the windows the project's
own verifier produced.  The windows are USED AS GIVEN (not recomputed) and then
checked to be trough +-10 % of the period; nothing is invented here.

Chosen cycles: the last two complete ones, cycles 25 and 26 (troughs 271.833 and
282.433 h).  Both of their read windows and BOTH flips fall inside
[266.533, 287.733] h = exactly two periods, and both edges lie inside the data, so
no window needs padding.  Cycle 25 reads 0 and cycle 26 reads 1, so one dwell in
each band is visible.

Why this script verifies itself
-------------------------------
The agent drawing this has no image input.  The checks below are the substitute:

  1. labels  -- export as SVG with `svg.fonttype='none'` (real <text>, editable for
     the art team) and assert every required label is literally in the file;
  2. REAL data -- assert the plotted arrays are element-wise equal to the CSV values
     (no smoothing, no idealisation);
  3. not a square wave -- assert the trajectory actually spends samples inside the
     forbidden middle band, and report the measured 10-90 % transition width;
  4. 20 % total -- assert every drawn window's width equals 0.20 x period to within
     one sample step, and that the window is centred on its trough;
  5. no padding -- assert the plotted range and every window lie inside the data, and
     that the windows are the verifier's own numbers unchanged;
  6. bands -- assert the shaded spans really are [0, 0.30] and [0.70, 1.0];
  7. layout -- text bounding boxes must not overlap each other, must stay inside the
     canvas, and must not sit on top of any plotted curve.

Run:  & 'D:\aconade\python.exe' .\make_fig03_2_switching_read_windows.py
Out:  ./out/fig03_2_switching_read_windows.png / .svg
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                                          # noqa: E402
import numpy as np                                                       # noqa: E402
import pandas as pd                                                      # noqa: E402

from tempo_style import (DEEP_BLUE, INK, LIGHT_BLUE, ORANGE, RISK_ORANGE,  # noqa: E402
                         TEAL)

HERE = Path(__file__).resolve().parent
OUT = HERE / 'out'
STEM = 'fig03_2_switching_read_windows'
FR = HERE.parents[1] / 'final_reconstruction'
TRAJ = FR / 'twobit34_results' / 'trajectories.csv'
WINDOWS = FR / 'twobit34_results' / 'read_windows.csv'

BAND_LOW, BAND_HIGH = 0.30, 0.70        # verify_twobit_causal.py
READ_WINDOW_FRACTION = 0.20             # TOTAL width, trough +-10 %
SVG_FONTTYPE = 'none'
C_LOW_BAND, C_HIGH_BAND = '#FBEDED', '#EAF6EC'
C_WINDOW = LIGHT_BLUE
REQUIRED_LABELS = ('Time (h)', 'LR fraction S0', 'Read window', 'Low band',
                   'High band', 'RDF0', 'Complex C0', 'Forward flux',
                   'Reverse flux')


def load():
    """The frozen trajectory + the verifier's own read windows, and the 2 cycles."""
    traj = pd.read_csv(TRAJ)
    win = pd.read_csv(WINDOWS)
    period = float(np.median(np.diff(win.trough_h.to_numpy())))
    # the last two consecutive COMPLETE cycles
    c1, c2 = win.iloc[-3], win.iloc[-2]
    x0 = float(c1.trough_h - period / 2.0)
    x1 = float(c2.trough_h + period / 2.0)
    used = win[(win.trough_h >= x0) & (win.trough_h <= x1)]
    return traj, win, used, period, x0, x1


def emptiest_x(sel, y, halfwidth, lo, hi, n=60):
    """Pick the x at which a label sitting at level `y` clears BOTH curves.

    Label placement must be data-driven, not guessed: at a fixed x the S0 trace may
    pass through the band label, or the (twin-axis) Int input pulse may reach that
    height.  Score each candidate by the smaller of the two clearances and take the
    best.  Deterministic, so the same data gives the same placement.
    """
    imax = float(sel.C31_flux_au_h.max())
    best, best_score = lo + halfwidth, -1.0
    for xc in np.linspace(lo + halfwidth, hi - halfwidth, n):
        w = sel[(sel.time_h >= xc - halfwidth) & (sel.time_h <= xc + halfwidth)]
        if len(w) < 2:
            continue
        s0_clear = float(np.min(np.abs(w.S0 - y)))          # distance to the S0 curve
        in_clear = (y * imax - float(w.C31_flux_au_h.max())) / max(imax, 1e-12)
        score = min(s0_clear, in_clear)
        if score > best_score:
            best, best_score = float(xc), score
    return best


def draw(traj, used, period, x0, x1):
    sel = traj[(traj.time_h >= x0) & (traj.time_h <= x1)].reset_index(drop=True)
    t = sel.time_h.to_numpy()
    series = {}

    fig, (axA, axB, axC) = plt.subplots(
        3, 1, figsize=(11.5, 10.0), sharex=True,
        gridspec_kw=dict(height_ratios=[1.25, 1.0, 1.0], hspace=0.14))

    # ---------------------------------------------------------------- panel A
    BANDS = ((0.0, BAND_LOW, C_LOW_BAND, 'Low band', 0.15),
             (BAND_HIGH, 1.0, C_HIGH_BAND, 'High band', 0.90))
    for lo, hi, col, lab, ty in BANDS:
        axA.axhspan(lo, hi, color=col, zorder=0)
    for lo, hi, col, lab, ty in BANDS:
        xc = emptiest_x(sel, ty, 0.75, x0, x1)
        axA.text(xc, ty, lab, fontsize=9.5,
                 color='#7a2f2f' if lab == 'Low band' else '#2f5d3a',
                 ha='center', va='center', zorder=6)
    for i, w in enumerate(used.itertuples()):
        axA.axvspan(w.window_start_h, w.window_end_h, color=C_WINDOW, zorder=1)
        if i == 0:
            axA.text((w.window_start_h + w.window_end_h) / 2, 0.60, 'Read window',
                     ha='center', va='center', fontsize=9.5, color='#3d5a63',
                     zorder=6)
        axA.text(w.trough_h, 0.50, str(int(w.bit0)), ha='center', va='center',
                 fontsize=15, color=DEEP_BLUE, fontweight='bold', zorder=6)
        axA.plot([w.trough_h, w.trough_h], [0.0, 1.0], color='#3d5a63', lw=0.8,
                 ls=':', zorder=2)
    axA.plot(t, sel.S0.to_numpy(), color=DEEP_BLUE, lw=2.0, label='LR fraction S0',
             zorder=5)
    series[axA] = [('LR fraction S0', t, sel.S0.to_numpy())]
    axA.set_ylim(0, 1)
    axA.set_ylabel('LR fraction S0', fontsize=11)
    axAin = axA.twinx()
    axAin.plot(t, sel.C31_flux_au_h.to_numpy(), color=TEAL, lw=1.7,
               label='Int input', zorder=4)
    series[axAin] = [('Int input', t, sel.C31_flux_au_h.to_numpy())]
    axAin.set_ylabel('Int input (C31 flux, a.u./h)', fontsize=10, color=TEAL)
    axAin.tick_params(axis='y', colors=TEAL)
    # Legends sit OUTSIDE the axes, above each panel: inside they would land on the
    # curves (in the upper right S0 is at ~0.98, exactly where a legend would go).
    axA.legend(handles=[axA.get_lines()[0], axAin.get_lines()[0]],
               labels=['LR fraction S0', 'Int input'],
               loc='lower left', bbox_to_anchor=(0.0, 1.0), ncol=2, fontsize=9,
               frameon=False)

    # ---------------------------------------------------------------- panel B
    rdf = axB.plot(t, sel.b0_R.to_numpy(), color=ORANGE, lw=1.9, label='RDF0',
                   zorder=4)[0]
    series[axB] = [('RDF0', t, sel.b0_R.to_numpy())]
    axB.set_ylabel('RDF0 (a.u.)', fontsize=11, color=ORANGE)
    axB.tick_params(axis='y', colors=ORANGE)
    axBin = axB.twinx()
    cpx = axBin.plot(t, sel.b0_C.to_numpy(), color=DEEP_BLUE, lw=1.9,
                     label='Complex C0', zorder=4)[0]
    series[axBin] = [('Complex C0', t, sel.b0_C.to_numpy())]
    axBin.set_ylabel('Complex C0 (a.u.)', fontsize=10, color=DEEP_BLUE)
    axBin.tick_params(axis='y', colors=DEEP_BLUE)
    axB.legend(handles=[rdf, cpx], labels=['RDF0', 'Complex C0'], loc='lower left',
               bbox_to_anchor=(0.0, 1.0), ncol=2, fontsize=9, frameon=False)
    axB.text(1.0, 1.045, 'RDF0 and Complex C0 on separate axes - not normalised',
             transform=axB.transAxes, ha='right', va='bottom', fontsize=8.5,
             color='#5a6b70', zorder=6)

    # ---------------------------------------------------------------- panel C
    fwd = axC.plot(t, sel.J_fwd0.to_numpy(), color=TEAL, lw=1.9,
                   label='Forward flux', zorder=4)[0]
    rev = axC.plot(t, sel.J_rev0.to_numpy(), color=RISK_ORANGE, lw=1.9,
                   label='Reverse flux', zorder=4)[0]
    series[axC] = [('Forward flux', t, sel.J_fwd0.to_numpy()),
                   ('Reverse flux', t, sel.J_rev0.to_numpy())]
    axC.set_ylabel('recombination flux (a.u./h)', fontsize=10)
    axC.legend(handles=[fwd, rev], labels=['Forward flux', 'Reverse flux'],
               loc='lower left', bbox_to_anchor=(0.0, 1.0), ncol=2, fontsize=9,
               frameon=False)
    axC.set_xlabel('Time (h)', fontsize=11)

    for ax in (axA, axB, axC):
        ax.grid(alpha=0.18)
        ax.set_xlim(x0, x1)
    return fig, (axA, axB, axC), series, sel, BANDS


def collect_texts(fig):
    """Every Text artist in the figure, including legend entries and axis labels."""
    return [t for t in fig.findobj(matplotlib.text.Text)
            if isinstance(t.get_text(), str) and t.get_text().strip()]


def annotation_texts(fig):
    """Texts a human placed: everything except the axis-managed tick labels.

    Tick labels of a twin axis are still present as Text artists even when they are
    not drawn, so including them produced thousands of phantom 'overlaps'.
    """
    tick_ids = set()
    for ax in fig.axes:
        for t in list(ax.get_xticklabels()) + list(ax.get_yticklabels()):
            tick_ids.add(id(t))
    return [t for t in collect_texts(fig)
            if id(t) not in tick_ids and t.get_visible()]


def verify(fig, axes, series, traj, win, used, period, x0, x1, bands, svg_path) -> int:
    fails = []

    def check(cond, msg):
        print(('  PASS  ' if cond else '  FAIL  ') + msg)
        if not cond:
            fails.append(msg)

    texts = [t.get_text() for t in collect_texts(fig)]

    print('-- 1. required labels, verbatim, in the rendered text')
    for lab in REQUIRED_LABELS:
        check(lab in texts, f'label present: {lab!r}')

    print('-- 2. the plotted curves ARE the CSV values (no idealisation)')
    sel = traj[(traj.time_h >= x0) & (traj.time_h <= x1)].reset_index(drop=True)
    pairs = [('LR fraction S0', 'S0'), ('Int input', 'C31_flux_au_h'),
             ('RDF0', 'b0_R'), ('Complex C0', 'b0_C'),
             ('Forward flux', 'J_fwd0'), ('Reverse flux', 'J_rev0')]
    for label, col in pairs:
        line = None
        for ax in fig.axes:
            for ln in ax.get_lines():
                if ln.get_label() == label:
                    line = ln
        if line is None:
            check(False, f'curve for {label!r} exists')
            continue
        want = sel[col].to_numpy()
        got = np.asarray(line.get_ydata(), dtype=float)
        check(got.shape == want.shape and np.array_equal(got, want),
              f'{label} is the raw {col} column, element-wise identical '
              f'({len(got)} samples)')

    print('-- 3. real switching, not a square wave')
    s = sel.S0.to_numpy()
    tt = sel.time_h.to_numpy()
    mid = int(np.sum((s > BAND_LOW) & (s < BAND_HIGH)))
    check(mid > 0, f'{mid} samples lie strictly between {BAND_LOW} and {BAND_HIGH} '
                   f'(an ideal step would have none)')
    check(bool(np.all((s > 0) & (s < 1))),
          'S0 is never exactly 0 or 1 (continuous dynamics, no idealised levels)')
    cross = np.flatnonzero((s[:-1] - 0.5) * (s[1:] - 0.5) < 0)
    check(len(cross) == 2, f'exactly 2 flips (S0 = 0.5 crossings) in the frame '
                           f'(got {len(cross)})')
    widths = []
    for i in cross:
        j = i
        while j > 0 and 0.1 < s[j] < 0.9:
            j -= 1
        k = i + 1
        while k < len(s) - 1 and 0.1 < s[k] < 0.9:
            k += 1
        widths.append(float(tt[k] - tt[j]))
    check(all(0.0 < w < period for w in widths),
          f'each flip takes a finite fraction of the period (10-90 % widths '
          f'{ [round(w, 2) for w in widths] } h of a {period:.2f} h period)')
    step = float(np.max(np.abs(np.diff(s))))
    check(step < 0.5, f'the flip is resolved on the 2 min grid, not drawn as a jump '
                      f'(max |dS0| between samples = {step:.4f})')

    print('-- 4. the read window is 20 % of the period in TOTAL, centred on the trough')
    for w in used.itertuples():
        width = float(w.window_end_h - w.window_start_h)
        centre = (w.window_start_h + w.window_end_h) / 2.0
        check(abs(width - READ_WINDOW_FRACTION * period) <= 1.5 * 0.0333,
              f'cycle {w.cycle}: width {width:.4f} h = {width / period:.4f} of the '
              f'{period:.3f} h period (target {READ_WINDOW_FRACTION})')
        check(abs(centre - w.trough_h) <= 1.5 * 0.0333,
              f'cycle {w.cycle}: window is centred on its trough '
              f'(off by {abs(centre - w.trough_h):.4f} h)')

    print('-- 5. no trailing window is padded or extended')
    check(len(used) == 2, f'exactly 2 read windows are drawn (got {len(used)})')
    check(x0 >= traj.time_h.min() and x1 <= traj.time_h.max(),
          f'the plotted range [{x0:.3f}, {x1:.3f}] h lies inside the data '
          f'[{traj.time_h.min():.1f}, {traj.time_h.max():.1f}] h')
    check(abs((x1 - x0) - 2 * period) <= 1.5 * 0.0333,
          f'the range spans exactly 2 periods ({x1 - x0:.4f} vs {2 * period:.4f} h)')
    for w in used.itertuples():
        check(w.window_start_h >= x0 and w.window_end_h <= x1,
              f'cycle {w.cycle}: the whole window is inside the plotted range')
    orig = win[(win.trough_h >= x0) & (win.trough_h <= x1)]
    check(np.array_equal(used.window_start_h.to_numpy(), orig.window_start_h.to_numpy())
          and np.array_equal(used.window_end_h.to_numpy(), orig.window_end_h.to_numpy()),
          'the windows are the verifier\'s own numbers, unchanged (nothing recomputed)')
    tail = win[win.window_end_h > traj.time_h.max()]
    check(len(tail) == 0,
          f'no window in the whole run falls outside the data (would need padding): '
          f'{len(tail)}')

    print('-- 6. the shaded bands are 0-0.30 and 0.70-1.00')
    want_bands = [(0.0, BAND_LOW), (BAND_HIGH, 1.0)]
    got_bands = sorted((round(float(lo), 4), round(float(hi), 4))
                       for lo, hi, _c, _l, _y in bands)
    check(got_bands == want_bands,
          f'the drawn band spans are [0, 0.30] and [0.70, 1.0] (got {got_bands})')
    # ... and the spans really exist as patches on the axes, at those levels
    patch_bands = sorted((round(float(p.get_y()), 4),
                          round(float(p.get_y() + p.get_height()), 4))
                         for p in axes[0].patches
                         if p.__class__.__name__ == 'Rectangle'
                         and round(float(p.get_height()), 4) in (BAND_LOW, 0.30))
    check(patch_bands == want_bands,
          f'shaded patches on the axes match (got {patch_bands})')

    print('-- 7. layout: no overlapping labels, inside the canvas, off the curves')
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()

    def clash(a, b, tol=0.5):
        """True when two display bboxes share more than `tol` pixels each way."""
        inter = matplotlib.transforms.Bbox.intersection(a, b)
        return inter is not None and inter.width > tol and inter.height > tol

    placed = annotation_texts(fig)
    items = [(t, t.get_window_extent(renderer=renderer)) for t in placed]
    overlaps = []
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            (ta, ba), (tb, bb) = items[i], items[j]
            if clash(ba, bb):
                overlaps.append((ta.get_text()[:22], tb.get_text()[:22]))
    check(not overlaps, f'no two placed labels overlap (found {len(overlaps)}: '
                        f'{overlaps[:4]})')
    # x tick labels of the shared axis must be readable too
    xt = [t for t in axes[2].get_xticklabels() if t.get_text() and t.get_visible()]
    xb = [t.get_window_extent(renderer=renderer) for t in xt]
    tick_clash = [i for i in range(len(xb) - 1) if clash(xb[i], xb[i + 1])]
    check(not tick_clash,
          f'the {len(xt)} x tick labels do not collide (clashes: {tick_clash})')
    figbb = fig.get_window_extent(renderer=renderer)
    outside = [t.get_text() for t, bb in items
               if bb.x0 < figbb.x0 - 1 or bb.x1 > figbb.x1 + 1
               or bb.y0 < figbb.y0 - 1 or bb.y1 > figbb.y1 + 1]
    check(not outside, f'every label is inside the canvas (outside: {outside})')

    # a label must not sit on top of a drawn curve in its own axes
    intruders = []
    for ax, curves in series.items():
        axbb = ax.get_window_extent(renderer=renderer)
        for t in placed:
            bb = t.get_window_extent(renderer=renderer)
            if not (bb.x0 >= axbb.x0 - 1 and bb.x1 <= axbb.x1 + 1
                    and bb.y0 >= axbb.y0 - 1 and bb.y1 <= axbb.y1 + 1):
                continue                      # not inside this axes
            inv = ax.transData.inverted()
            (dx0, dy0) = inv.transform((bb.x0, bb.y0))
            (dx1, dy1) = inv.transform((bb.x1, bb.y1))
            for name, xs, ys in curves:
                inside = (xs >= min(dx0, dx1)) & (xs <= max(dx0, dx1)) & \
                         (ys >= min(dy0, dy1)) & (ys <= max(dy0, dy1))
                if inside.any():
                    intruders.append((t.get_text()[:22], name, int(inside.sum())))
    check(not intruders,
          f'no label sits on a curve (found {len(intruders)}: {intruders[:4]})')

    print('-- 8. the exported SVG carries the labels as real text')
    svg = svg_path.read_text(encoding='utf-8')
    check('<text' in svg, "SVG uses <text> (svg.fonttype = 'none')")
    for lab in REQUIRED_LABELS:
        esc = lab.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')
        check(esc in svg, f'SVG carries the label: {lab!r}')

    print()
    print(f'FAILURES: {len(fails)}')
    for f in fails:
        print('  -', f)
    return len(fails)


def main() -> int:
    plt.rcParams['svg.fonttype'] = SVG_FONTTYPE
    traj, win, used, period, x0, x1 = load()
    print(f'cycles used: {list(used.cycle)}  troughs {list(used.trough_h)}  '
          f'bit0 {list(used.bit0)}')
    print(f'period {period:.4f} h | range [{x0:.3f}, {x1:.3f}] h | '
          f'{x1 - x0:.4f} h = {(x1 - x0) / period:.4f} periods')
    fig, axes, series, sel, bands = draw(traj, used, period, x0, x1)
    OUT.mkdir(parents=True, exist_ok=True)
    png, svg = OUT / f'{STEM}.png', OUT / f'{STEM}.svg'
    for p in (png, svg):
        fig.savefig(p, dpi=200, bbox_inches='tight',
                    facecolor=fig.get_facecolor())
    fails = verify(fig, axes, series, traj, win, used, period, x0, x1, bands, svg)
    plt.close(fig)
    print(f'wrote {png}')
    print(f'wrote {svg}')
    return 1 if fails else 0


if __name__ == '__main__':
    sys.exit(main())
