r"""Figure 4-3 -- A sharper gate enables three-bit counting   (dshwork, per the brief)

Brief: `dshwork/精简交付版/02_作图要求精简版.md` section "Figure 4-3 -- A sharper gate
enables three-bit counting".  Requirements:

  * Panel A: under the SAME expression settings, the S2 trajectory and read windows at
    n=4 and at n=6; mark the former as failing and the latter as passing; show only
    ONE common time span, keeping the real transitions and the low state;
  * Panel B: the paired `L_symmetric` for n = 4, 5, 6, 7, 8 -- four ACTUAL points per
    exponent plus its median, with a star on the selected n=6;
  * the caption must carry `5% tail cut; 2 min output grid; 4 model configurations per
    exponent`;
  * MUST be right: those four points are four PARAMETER COMBINATIONS, not replicate
    experiments; n=4 fails through real dynamics, not merely through a scoring error;
    n=5 is only the first passing TESTED integer, not a located exact continuous
    boundary;
  * labels verbatim: `Original gate (n=4)`, `Independent gate exponent (n=6)`,
    `Not certified`, `Certified`, `Paired wrong-direction flux ratio`, `Selected`.

Data (read-only):
  n=6 S2 + windows  threebit51_results/wiki_n6_20260924/{trajectories,read_windows}.csv
  n=4 S2 + windows  data/n4_matched_600h.{csv,json}   <- produced by
                    make_n4_matched_trajectory.py (a NEW matched run; its provenance
                    JSON is checked here, including the realised gate exponent)
  Panel B           plausibility/carry_pairing_grid_verdict.json (per_point + by_n)

The bit-2 curve is read from the column named `b2_S` in BOTH files -- the model's own
name for that state -- never from a bare `S2`.

CORRECTION (this version)
-------------------------
An earlier version of `make_n4_matched_trajectory.py` wrote its columns as
`S0 = y[44], S1 = y[45], S2 = y[46]`, i.e. it put **b2_S, A1 and F1** under the names
S0, S1, S2 (the real indices are S0 = 16, S1 = 27, S2 = 44).  Panel A therefore drew
bit 2's **F1 protein** as the n=4 "S2" curve next to bit 2's real S2 for n=6 -- two
different physical quantities under one axis label -- and the note "S2 peaks at 2.25;
the axis is clipped at 1.05" was describing F1, a concentration in a.u. rather than a
bounded fraction.  Both files are fixed: the generator now resolves state indices by
name and re-verifies its own output, and this figure additionally refuses to draw a
curve that leaves [0, 1].

What did NOT change: the verdict numbers.  `crossings = 5`, `gates = 14`, the 17
never-committing windows and `certified = False` all come from
`verify_threebit51.analyse_threebit`, which reads `states = (y[16], y[27], y[44])` and so
was always correct.  They are unaffected by the plotting-column bug.

Common span: [30, 150] h.  At n=4 bit 2 makes its one real 0 -> 1 transition at about
38 h and then never returns to the low state; n=6 alternates cleanly over the same span.

Self-checks: labels verbatim; the S2 curves are the raw `b2_S` columns over the span;
neither curve leaves [0, 1] (a wrong column such as F1 fails this immediately); the n=4
file passes `check_state_column_semantics.check_file`, which replays the recorded
read-window bit-2 labels from the plotted column; the n=4 run is the MATCHED one (same
extension/carry, only the exponent differs) and it is genuinely uncertified with real
dynamics rather than scored as a failure; Panel B has 4 points per n plus the median,
medians equal the frozen values 0.4489/0.0521/0.0169/0.0139/0.0131, the star is on n=6,
and the "parameter combinations, not replicates" and "first certified tested integer"
notes are present; layout and SVG checks.

Run:  & 'D:\aconade\python.exe' .\make_fig04_3_sharper_gate.py
Out:  ./out/fig04_3_sharper_gate.png / .svg
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                                          # noqa: E402
import numpy as np                                                       # noqa: E402
import pandas as pd                                                      # noqa: E402

from tempo_style import DEEP_BLUE, INK, LIGHT_BLUE, ORANGE, RISK_ORANGE, TEAL  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE / 'out'
DATA = HERE / 'data'
STEM = 'fig04_3_sharper_gate'
FR = HERE.parents[1] / 'final_reconstruction'
N6DIR = FR / 'threebit51_results' / 'wiki_n6_20260924'
N4CSV = DATA / 'n4_matched_600h.csv'
N4JSON = DATA / 'n4_matched_600h.json'
PAIRING = FR / 'plausibility' / 'carry_pairing_grid_verdict.json'

SVG_FONTTYPE = 'none'
REQUIRED_LABELS = ('Original gate (n=4)', 'Independent gate exponent (n=6)',
                   'Not certified', 'Certified',
                   'Paired wrong-direction flux ratio', 'Selected')
CAPTION_LINE = ('5% tail cut; 2 min output grid; 4 model configurations per exponent')
MEDIANS = {4: 0.4489, 5: 0.0521, 6: 0.0169, 7: 0.0139, 8: 0.0131}
# The model's own name for bit 2's DNA-conformation state.  Read the curve from this
# column in BOTH files; a bare `S2` header is exactly what went wrong before.
S2_COLUMN = 'b2_S'
T0, T1 = 30.0, 150.0         # the one common time span
HOURS_FULL = 600.0           # the run length every whole-run claim is measured over
SPAN_NOTE = ('span [30, 150] h: a startup-phase example; the whole-run (600 h) '
             'evidence is in the caption and the provenance JSON')
BAND_LOW, BAND_HIGH, MIN_OCCUPANCY = 0.30, 0.70, 0.80


def _bitstr(v) -> str:
    """Normalise a recorded bit label to '0', '1' or 'x'.

    The labels come from two different files (a JSON list and a CSV whose `x` entries
    make the column object-typed), so a value can arrive as an int, a float, a numpy
    scalar or a string.
    """
    s = str(v)
    if s in ('0', '1'):
        return s
    try:
        f = float(s)
    except (TypeError, ValueError):
        return 'x'
    return '0' if f == 0.0 else ('1' if f == 1.0 else 'x')


def load():
    if not N4CSV.is_file():
        raise SystemExit(f'FATAL: {N4CSV} is missing; run '
                         f'make_n4_matched_trajectory.py first (it is the matched '
                         f'n=4 run this figure compares against)')
    n4 = pd.read_csv(N4CSV)
    n4meta = json.loads(N4JSON.read_text(encoding='utf-8'))
    n4w = pd.DataFrame(n4meta['read_windows'])
    n6 = pd.read_csv(N6DIR / 'trajectories.csv')
    n6w = pd.read_csv(N6DIR / 'read_windows.csv')
    pairing = json.loads(PAIRING.read_text(encoding='utf-8'))
    return n4, n4meta, n4w, n6, n6w, pairing


def figure(D):
    n4, n4meta, n4w, n6, n6w, pairing = D
    fig = plt.figure(figsize=(13.4, 9.8), constrained_layout=True)
    # Row 0 is a dedicated annotation band: the panel titles, the status words and the
    # qualifiers live there, so no annotation can land on a curve or on a title.  Its
    # height is not cosmetic: y = 0.86 / 0.56 / 0.26 are axes fractions, so a band that
    # is too short makes the title and the status word collide (the layout check below
    # fires on exactly that).
    gs = fig.add_gridspec(3, 2, height_ratios=[0.26, 1.0, 1.05], hspace=0.24,
                          wspace=0.20)
    axN = fig.add_subplot(gs[0, :])
    axN.axis('off')
    ax1 = fig.add_subplot(gs[1, 0])
    ax2 = fig.add_subplot(gs[1, 1], sharey=ax1)
    ax3 = fig.add_subplot(gs[2, :])

    counts6 = pairing['per_point']['n=6|mRNA=2|mat=32.5']
    exc = n4meta['excursion']
    for ax, traj, win, n, title, cert, xa, cross, gates, undef in (
            (ax1, n4, n4w, 4, 'Original gate (n=4)', bool(n4meta['certified']), 0.02,
             int(n4meta['crossings']), int(n4meta['gates']),
             int(n4meta['undefined_bit_labels'])),
            (ax2, n6, n6w, 6, 'Independent gate exponent (n=6)',
             bool(counts6['certified']), 0.52,
             int(counts6['crossings']), int(counts6['gates']), 0)):
        m = (traj.time_h >= T0) & (traj.time_h <= T1)
        t = traj.time_h[m]
        s2 = traj[S2_COLUMN][m]
        ax.plot(t, s2, color=DEEP_BLUE if n == 6 else RISK_ORANGE, lw=1.6)
        for w in win.itertuples():
            if w.window_end_h < T0 or w.window_start_h > T1:
                continue
            ax.axvspan(max(w.window_start_h, T0), min(w.window_end_h, T1),
                       color=LIGHT_BLUE, zorder=0)
        ax.axhspan(0.0, 0.30, color='#FBEDED', zorder=0)
        ax.axhspan(0.70, 1.0, color='#EAF6EC', zorder=0)
        ax.set_ylim(-0.05, 1.05)
        ax.set_xlim(T0, T1)
        ax.set_xlabel('Time (h)', fontsize=10)
        ax.grid(alpha=0.18)
        ax.set_title('', fontsize=1)
        # annotations go in the band above the panels (see the gridspec comment);
        # each panel owns one column of the band, so the texts cannot run into the
        # other panel's annotations
        axN.text(xa, 0.86, title, transform=axN.transAxes, ha='left', va='center',
                 fontsize=11.5, color=INK, fontweight='bold')
        axN.text(xa, 0.56, 'Certified' if cert else 'Not certified',
                 transform=axN.transAxes, ha='left', va='center', fontsize=11,
                 color=TEAL if cert else RISK_ORANGE, fontweight='bold')
        # the readout summary INSIDE the plotted span, counted from the same window
        # labels the certification used -- this is the panel's visual claim
        in_span = [_bitstr(getattr(w, 'bit2')) for w in win.itertuples()
                   if w.window_end_h >= T0 and w.window_start_h <= T1]
        n1 = sum(1 for v in in_span if v == '1')
        n0 = sum(1 for v in in_span if v == '0')
        nx = len(in_span) - n1 - n0
        if cert:
            qual = (f'every one of {gates} carries flips bit 2\n'
                    f'{undef} read windows never commit\n'
                    f'in this span: {n1} read 1, {n0} read 0, {nx} never commit')
        else:
            # The whole-run evidence, not just this span: from the provenance JSON's
            # `excursion` block, which the data generator measures over all 600 h and
            # refuses to write if the "never returns to the low band" clause stops
            # holding.  Three lines only -- the annotation band is narrow.
            qual = (f'S2 flips on only {cross} of {gates} carries; '
                    f'{undef} read windows never commit\n'
                    f'enters the high band once (t = {exc["t_first_high_h"]:.1f} h) and '
                    f'never returns to the low band in {HOURS_FULL:.0f} h '
                    f'(min {exc["s2_min_after_first_high"]:.2f})\n'
                    f'in this span: {n1} read 1, {n0} read 0, {nx} never commit')
        axN.text(xa, 0.26, qual, transform=axN.transAxes, ha='left', va='top',
                 fontsize=8.4, color='#3d5a63')
    # The span note goes at the FOOT of the figure, not in the annotation band: the band
    # is only tall enough for the panel titles, the status words and the three-line
    # qualifiers, and putting it there made it collide with both qualifiers (caught by
    # the layout check below).
    fig.supxlabel(SPAN_NOTE, fontsize=8.8, color='#3d5a63')
    ax1.set_ylabel('S2 (LR fraction of bit 2)', fontsize=10)

    # ---------------------------------------------------------------- Panel B
    pp = pairing['per_point']
    for n in range(4, 9):
        xs = [n] * 4
        ys = [pp[k]['5pct_L_symmetric_median'] for k in pp
              if k.startswith(f'n={n}|')]
        ax3.scatter(xs, ys, s=46, facecolor='white', edgecolor=DEEP_BLUE,
                    linewidth=1.3, zorder=4)
        med = pairing['by_n_A1_gate'][str(n)]['L_symmetric_5pct_median']
        ax3.plot([n - 0.22, n + 0.22], [med, med], color=INK, lw=2.4, zorder=5)
        ax3.annotate(f'{med:.4f}', xy=(n, med), xytext=(6, -3),
                     textcoords='offset points', fontsize=8.5, color=INK)
    ax3.scatter([6], [pairing['by_n_A1_gate']['6']['L_symmetric_5pct_median']],
                s=340, marker='*', facecolor='#f59f00', edgecolor='k',
                linewidth=0.9, zorder=7, label='Selected')
    ax3.set_yscale('log')
    ax3.set_xlabel('n_A1_gate  (independent gate exponent)', fontsize=10)
    ax3.set_ylabel('Paired wrong-direction flux ratio', fontsize=10)
    ax3.set_title('Paired R/F super-period L_symmetric at the frozen settings: four '
                  'actual configurations per exponent (open circles) and their median '
                  '(bar)', fontsize=11)
    ax3.legend(fontsize=9, loc='upper right')
    ax3.grid(alpha=0.22, which='both')
    ax3.text(0.015, 0.06,
             CAPTION_LINE + '\n'
             'The four points per exponent are four PARAMETER COMBINATIONS, not '
             'replicate experiments.\n'
             'n=4 fails through real dynamics, not a scoring rule; n=5 is the first '
             'certified TESTED integer,\nnot a located exact continuous boundary.',
             transform=ax3.transAxes, ha='left', va='bottom', fontsize=8.6,
             color='#3d5a63')
    return fig, dict(ax1=ax1, ax2=ax2, ax3=ax3, n4=n4, n6=n6, n4meta=n4meta,
                     pairing=pairing)


def collect_texts(fig):
    return [t for t in fig.findobj(matplotlib.text.Text)
            if isinstance(t.get_text(), str) and t.get_text().strip()]


def verify(fig, S, svg_path) -> int:
    fails = []

    def check(cond, msg):
        print(('  PASS  ' if cond else '  FAIL  ') + msg)
        if not cond:
            fails.append(msg)

    texts = [t.get_text() for t in collect_texts(fig)]
    print('-- 1. required labels, verbatim')
    for lab in REQUIRED_LABELS:
        check(lab in texts, f'label present: {lab!r}')
    check(any(CAPTION_LINE in t for t in texts),
          f'the figure carries the caption line {CAPTION_LINE!r}')

    print('-- 2. Panel A: a MATCHED n=4 run, and a real dynamical failure')
    meta = S['n4meta']
    check(meta['gate_exponent']['requested'] == 4.0
          and meta['gate_exponent']['realised_by_model'] == 4.0,
          'the n=4 run realised gate exponent 4.0 (asserted in its provenance)')
    check('frozen' in meta['expression_settings']['source'],
          'the n=4 run used the same frozen expression settings as the n=6 run')
    check(meta['certified'] is False,
          'the matched n=4 run is NOT certified (its own verdict)')
    # The matched run DOES form gates and does flip bit2 a few times; what fails is
    # that the flips do not track the carries and many read windows never commit.
    # (An earlier version of this check asserted 0 crossings, which came from a 30 h
    # probe that was simply too short - the 600 h run has crossings=5, gates=14.)
    check(meta['crossings'] < meta['gates'],
          f"n=4 completes {meta['crossings']} bit-2 crossings against "
          f"{meta['gates']} carry gates - the flips do not track the carries")
    check(meta['undefined_bit_labels'] > 0,
          f"n=4 leaves {meta['undefined_bit_labels']} bit labels undefined (the level "
          f"never commits inside those windows) - a dynamical failure, not a scoring "
          f"artefact")
    # The caption's mechanistic clause must rest on the whole-run measurement, not on
    # this 120 h window.  Both halves are checked against the provenance JSON, and the
    # data generator refuses to write that block unless the first half holds.
    exc = meta['excursion']
    check(exc['returns_to_low_band_after_first_high'] is False,
          f"S2 never returns to the low band after first entering the high band at "
          f"t = {exc['t_first_high_h']:.2f} h (whole-run min "
          f"{exc['s2_min_after_first_high']:.4f} > {exc['band_low']})")
    check(exc['int_J_rev2_after_first_high'] is not None
          and exc['int_J_rev2_whole'] is not None,
          'the whole-run reverse-flux integrals are recorded - the evidence behind the '
          'claim that the reverse write is insufficient')
    check(meta['status'].startswith('NEW DETERMINISTIC RUN'),
          'the n=4 data file labels itself as a new run, not a pre-existing artefact')
    n6cert = S['pairing']['per_point']['n=6|mRNA=2|mat=32.5']['certified']
    check(bool(n6cert) is True, 'the n=6 configuration is certified')

    print('-- 3. the two curves are the raw b2_S columns, and they are FRACTIONS')
    n4, n6 = S['n4'], S['n6']
    for name, traj, ax in (('n=4', n4, S['ax1']), ('n=6', n6, S['ax2'])):
        m = (traj.time_h >= T0) & (traj.time_h <= T1)
        want = traj[S2_COLUMN][m].to_numpy(dtype=float)
        ln = ax.get_lines()[0]
        got = np.asarray(ln.get_ydata(), dtype=float)
        check(got.shape == want.shape and np.array_equal(got, want),
              f'{name} curve is the raw {S2_COLUMN} column over [{T0}, {T1}] h '
              f'({len(got)} samples)')
        check(np.allclose(ax.get_xlim(), (T0, T1)),
              f'{name} shows the same time span [{T0}, {T1}] h')
        # A DNA-conformation fraction cannot leave [0, 1].  The previous version plotted
        # F1 here (max 2.39 a.u.); this one line catches that entire class of error.
        check(float(got.min()) >= -1e-9 and float(got.max()) <= 1.0 + 1e-9,
              f'{name} curve stays inside [0, 1] as a fraction must '
              f'(min {got.min():.4f}, max {got.max():.4f})')
        check(want.min() < 0.30 and want.max() > 0.70,
              f'{name} keeps both the low state and the high state '
              f'(min {want.min():.3f}, max {want.max():.3f})')

    print('-- 3b. the n=6 artefact agrees with itself: S2 == b2_S')
    if 'S2' in n6.columns and S2_COLUMN in n6.columns:
        d = float(np.abs(n6['S2'].to_numpy(dtype=float)
                         - n6[S2_COLUMN].to_numpy(dtype=float)).max())
        check(d == 0.0,
              f'the derived n=6 S2 column equals b2_S entry-wise (max|d| = {d:g})')
    else:
        check(False, f'the n=6 artefact carries both S2 and {S2_COLUMN}')

    print('-- 3c. COLUMN SEMANTICS: each file means what its header claims')
    sys.path.insert(0, str(HERE))
    from check_state_column_semantics import (check_file, state_index_map,  # noqa: E402
                                              windows_from_json)
    n4_gate = float(meta['gate_exponent']['realised_by_model'])
    idx4, model4 = state_index_map(n4_gate)
    sem_fails = check_file(N4CSV, windows_from_json(N4JSON), model4, idx4,
                           n_gate=n4_gate)
    check(sem_fails == 0,
          f'the n=4 matched file passes the column-semantics checker '
          f'({sem_fails} failure(s); see its output above)')

    print('-- 4. Panel B: four real points per exponent, medians, and the star')
    ax3 = S['ax3']
    pp = S['pairing']['per_point']
    for n in range(4, 9):
        want = sorted(pp[k]['5pct_L_symmetric_median'] for k in pp
                      if k.startswith(f'n={n}|'))
        check(len(want) == 4, f'n={n} has 4 configurations in the artefact')
        med = S['pairing']['by_n_A1_gate'][str(n)]['L_symmetric_5pct_median']
        check(abs(med - MEDIANS[n]) < 5e-5,
              f'n={n} median = {med:.4f} matches the frozen value {MEDIANS[n]:.4f}')
    # one scatter per exponent (4 configurations each), plus the single-point star
    sized = [(len(c.get_offsets()), c) for c in ax3.collections
             if c.get_offsets() is not None]
    n_points = sum(k for k, _c in sized if k == 4)
    check(n_points == 20,
          f'20 configuration points are drawn, four per exponent (got {n_points}; '
          f'sizes {[k for k, _c in sized]})')
    scatters = [c for k, c in sized if k == 4]
    if scatters:
        xs = np.concatenate([np.asarray(c.get_offsets())[:, 0] for c in scatters])
        check(sorted(set(np.round(xs, 3))) == [4, 5, 6, 7, 8],
              f'points are placed at n = 4..8 (got {sorted(set(xs))})')
    stars = [c for c in ax3.collections
             if c.get_offsets() is not None and len(c.get_offsets()) == 1]
    check(len(stars) == 1, 'exactly one star is drawn')
    if stars:
        sx, sy = np.asarray(stars[0].get_offsets())[0]
        want = S['pairing']['by_n_A1_gate']['6']['L_symmetric_5pct_median']
        check(abs(sx - 6.0) < 1e-9 and abs(sy - want) < 1e-12,
              f'the star marks n=6 at its median ({want:.4f})')
    leg = ax3.get_legend()
    check(leg is not None
          and 'Selected' in [t.get_text() for t in leg.get_texts()],
          'the legend names the star as Selected')
    joined = ' '.join(texts)
    for must in ('four PARAMETER COMBINATIONS, not',
                 'n=5 is the first certified TESTED integer'):
        check(must in joined, f'the required caveat is present: {must!r}')

    print('-- 5. layout')
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()

    def clash(a, b, tol=0.5):
        inter = matplotlib.transforms.Bbox.intersection(a, b)
        return inter is not None and inter.width > tol and inter.height > tol

    tick_ids = set()
    for ax in fig.axes:
        for t in list(ax.get_xticklabels()) + list(ax.get_yticklabels()):
            tick_ids.add(id(t))
    placed = [(t, t.get_window_extent(renderer=renderer))
              for t in collect_texts(fig)
              if t.get_visible() and id(t) not in tick_ids]
    pairs = [(placed[i][0].get_text()[:24], placed[j][0].get_text()[:24])
             for i in range(len(placed)) for j in range(i + 1, len(placed))
             if clash(placed[i][1], placed[j][1])]
    check(not pairs, f'no two labels overlap (found {len(pairs)}: {pairs[:4]})')
    figbb = fig.get_window_extent(renderer=renderer)
    outside = [t.get_text() for t, bb in placed
               if bb.x0 < figbb.x0 - 1 or bb.x1 > figbb.x1 + 1
               or bb.y0 < figbb.y0 - 1 or bb.y1 > figbb.y1 + 1]
    check(not outside, f'every figure label is inside the canvas (outside: {outside})')

    print('-- 6. the SVG carries the labels as real text')
    svg = svg_path.read_text(encoding='utf-8')
    check('<text' in svg, "SVG uses <text> (svg.fonttype = 'none')")
    for lab in REQUIRED_LABELS:
        esc = lab.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')
        check(esc in svg, f'SVG carries: {lab!r}')
    check('5% tail cut' in svg, 'the SVG carries the caption line')
    check('never returns to the low band' in svg,
          'the SVG states the WHOLE-RUN claim (never returns to the low band), not just '
          'what the plotted 120 h window shows')
    check('startup-phase example' in svg,
          'the SVG labels the plotted span as a startup-phase example')

    print()
    print(f'FAILURES: {len(fails)}')
    for f in fails:
        print('  -', f)
    return len(fails)


def main() -> int:
    try:
        sys.stdout.reconfigure(encoding='utf-8', errors='replace')
    except Exception:                                                # noqa: BLE001
        pass
    plt.rcParams['svg.fonttype'] = SVG_FONTTYPE
    D = load()
    print(f"n=4 matched run: certified={D[1]['certified']} "
          f"crossings={D[1]['crossings']} gates={D[1]['gates']}")
    fig, S = figure(D)
    OUT.mkdir(parents=True, exist_ok=True)
    png, svg = OUT / f'{STEM}.png', OUT / f'{STEM}.svg'
    for p in (png, svg):
        fig.savefig(p, dpi=200, bbox_inches='tight',
                    facecolor=fig.get_facecolor())
    fails = verify(fig, S, svg)
    plt.close(fig)
    print(f'wrote {png}')
    print(f'wrote {svg}')
    return 1 if fails else 0


if __name__ == '__main__':
    sys.exit(main())
