r"""Figure 4-2a -- Two-bit counting under the real Han oscillator  (dshwork)

Why this figure
---------------
The deliverable set had no panel showing Han's oscillator itself running together with the
two-bit counter completing whole mod-4 cycles.  Figure 4-2 ("From modulo 4 to modulo 8")
plots only `S0`/`S1`/`S2` and the read strip -- it contains no upstream trace at all -- and
`twobit34_results/figures/01_counter_overview` shows the upstream only as a C31 flux line,
with the oscillator's own state variables absent.  This figure closes that gap.

Three panels (the brief allows at most three), one shared time axis:

  A  the upstream oscillator: `TetR total`, `CI`, `LacI` (each normalized to its own
     maximum) with the C31 translation flux that drives bit 0 on a second y axis;
  B  the two DNA-conformation states `b0_S` / `b1_S` with the 0-0.30 / 0.70-1 read bands
     and the 28 read windows;
  C  the decoded value, i.e. the mod-4 count itself, as a staircase over those windows.

Data (read-only):  data/twobit_selected_300h.{csv,json}, produced by
`make_twobit_under_upstream_run.py` at the SHARED selected profile (uM_per_au = 5.75,
carry A/F maturation = 32.5 min, carry mRNA = 2 min).  That script asserts every field of
the model against `twobit34_results/selected_profile.json` and asserts the run is not on
the bare `nominal_extension()` default (6.0 / 30.0), which is what `run_twobit34.py` uses.

Figure discipline (§0.3 of the outline): every result figure carries model / input /
parameter version.  Here that is a footer block printed from the provenance JSON, not
typed by hand.

Self-checks: curves identical to the raw CSV columns; S0/S1 inside [0, 1]; the staircase
equals `read_windows.value` element-wise; the sequence is 7 whole mod-4 cycles; the
provenance JSON says the working point is the selected profile and that the run
reproduces the 75-point grid row; the three-element footer is present verbatim; layout and
SVG-text checks.  No leakage number appears anywhere: the two-bit certification has no
leak arm, so there is nothing to quote.

Run:  & 'D:\aconade\python.exe' .\make_fig04_2a_twobit_under_upstream.py
Out:  ./out/fig04_2a_twobit_under_upstream.png / .svg / .pdf
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

from tempo_style import (DEEP_BLUE, INK, LIGHT_BLUE, ORANGE,  # noqa: E402
                         RISK_ORANGE, TEAL)

HERE = Path(__file__).resolve().parent
OUT = HERE / 'out'
DATA = HERE / 'data'
STEM = 'fig04_2a_twobit_under_upstream'
CSV = DATA / 'twobit_selected_300h.csv'
META = DATA / 'twobit_selected_300h.json'

SVG_FONTTYPE = 'none'
BAND_LOW, BAND_HIGH = 0.30, 0.70
OSCILLATOR = (('TetR_total', 'TetR total'), ('CI', 'CI'), ('LacI', 'LacI'))
EXPECTED_SEQUENCE = '1230' * 7

REQUIRED_LABELS = (
    'Han oscillator (upstream)',
    'C31 translation flux',
    'Two-bit DNA state',
    'Low band', 'High band',
    'Decoded value',
    'Time (h)',
)
FOOTER_PREFIXES = ('Model:', 'Input:', 'Parameters:')


def load():
    if not CSV.is_file():
        raise SystemExit(f'FATAL: {CSV} is missing; run '
                         f'make_twobit_under_upstream_run.py first')
    return (pd.read_csv(CSV), json.loads(META.read_text(encoding='utf-8')))


def normalised(v):
    v = np.asarray(v, dtype=float)
    peak = float(np.max(np.abs(v)))
    return v / peak if peak > 0 else v


def figure(df, meta):
    t = df.time_h.to_numpy(dtype=float)
    fig = plt.figure(figsize=(13.0, 10.6), constrained_layout=True)
    # Row 0 is a dedicated annotation band (figure title, working point, status); the
    # bottom band carries the three-element footer.  Keeping both out of the plot rows
    # means no annotation can land on a curve.
    gs = fig.add_gridspec(5, 1, height_ratios=[0.20, 1.0, 1.0, 0.85, 0.22],
                          hspace=0.30)
    axN = fig.add_subplot(gs[0]); axN.axis('off')
    axA = fig.add_subplot(gs[1])
    axB = fig.add_subplot(gs[2], sharex=axA)
    axC = fig.add_subplot(gs[3], sharex=axA)
    axF = fig.add_subplot(gs[4]); axF.axis('off')

    # ------------------------------------------------------------------ Panel A
    for (col, label), colour in zip(OSCILLATOR, (DEEP_BLUE, ORANGE, TEAL)):
        axA.plot(t, normalised(df[col]), lw=1.5, color=colour, label=label)
    axA.set_ylabel('upstream, each trace\nnormalized to its own max', fontsize=10)
    axA.set_ylim(-0.05, 1.12)
    axA.legend(ncol=3, loc='upper left', fontsize=9, framealpha=0.9)
    axA.text(0.995, 0.94, 'Han oscillator (upstream)', transform=axA.transAxes,
             ha='right', va='top', fontsize=11.5, color=INK, fontweight='bold')
    axA.grid(alpha=0.18)

    axA2 = axA.twinx()
    axA2.plot(t, df.C31_flux_au_h, lw=1.2, color=RISK_ORANGE, alpha=0.85)
    axA2.set_ylabel('C31 translation flux\n(a.u./h)', fontsize=10, color=RISK_ORANGE)
    axA2.tick_params(axis='y', labelsize=9, colors=RISK_ORANGE)
    axA2.text(0.995, 0.80, 'C31 translation flux', transform=axA2.transAxes,
              ha='right', va='top', fontsize=9.5, color=RISK_ORANGE)

    # ------------------------------------------------------------------ Panel B
    for col, label, colour in (('b0_S', 'bit 0 (S0)', DEEP_BLUE),
                               ('b1_S', 'bit 1 (S1)', ORANGE)):
        axB.plot(t, df[col], lw=1.6, color=colour, label=label)
    for w in meta['read_windows']:
        axB.axvspan(w['window_start_h'], w['window_end_h'], color=LIGHT_BLUE, zorder=0)
    axB.axhspan(0.0, BAND_LOW, color='#FBEDED', zorder=0)
    axB.axhspan(BAND_HIGH, 1.0, color='#EAF6EC', zorder=0)
    axB.set_ylim(-0.05, 1.12)
    axB.set_ylabel('LR fraction S', fontsize=10)
    axB.text(0.995, 0.94, 'Two-bit DNA state  (read windows shaded)', transform=axB.transAxes,
             ha='right', va='top', fontsize=11.5, color=INK, fontweight='bold')
    # The band labels sit at the RIGHT edge: the legend owns the upper left of this panel,
    # and putting them there made them collide with the 'bit 0 (S0)' legend entry.
    axB.text(0.995, 0.10, 'Low band', transform=axB.transAxes, ha='right', va='bottom',
             fontsize=9, color='#8a6b6b')
    axB.text(0.995, 0.78, 'High band', transform=axB.transAxes, ha='right', va='bottom',
             fontsize=9, color='#5d7d63')
    axB.legend(ncol=2, loc='upper left', fontsize=9, framealpha=0.9)
    axB.grid(alpha=0.18)

    # ------------------------------------------------------------------ Panel C
    troughs = np.asarray([w['trough_h'] for w in meta['read_windows']], dtype=float)
    values = np.asarray([np.nan if w['value'] is None else w['value']
                         for w in meta['read_windows']], dtype=float)
    edges = np.r_[troughs[0] - (troughs[1] - troughs[0]) / 2,
                  (troughs[:-1] + troughs[1:]) / 2,
                  troughs[-1] + (troughs[-1] - troughs[-2]) / 2]
    axC.step(edges, np.r_[values, values[-1]], where='post', lw=1.9, color=INK)
    axC.scatter(troughs, values, s=44, facecolor='white', edgecolor=TEAL, linewidth=1.6,
                zorder=5)
    for x, v in zip(troughs, values):
        if np.isfinite(v):
            axC.annotate(f'{int(v)}', xy=(x, v), xytext=(0, 8), textcoords='offset points',
                         ha='center', fontsize=8.4, color=INK)
    axC.set_yticks([0, 1, 2, 3])
    # Headroom above 3: the per-window value labels are offset upward from y = 3, and with
    # a tight top they ran into this panel's title.
    axC.set_ylim(-0.55, 4.40)
    axC.set_ylabel('Decoded value', fontsize=10)
    axC.set_xlabel('Time (h)', fontsize=10)
    axC.text(0.995, 0.94, 'Decoded value: seven whole mod-4 cycles', transform=axC.transAxes,
             ha='right', va='top', fontsize=11.5, color=INK, fontweight='bold')
    seq = meta['cold_start']['sequence']
    # Top-left of this panel: the value labels stick up from y = 0 along the bottom, and
    # the title is right-aligned, so the top-left corner is the free one.
    axC.text(0.012, 0.97, f'read sequence ({len(seq)} windows): {seq}',
             transform=axC.transAxes, ha='left', va='top', fontsize=8.6,
             color='#3d5a63', family='monospace')
    axC.grid(alpha=0.18)

    # ------------------------------------------------------------- bands and footer
    wp = meta['working_point']
    axN.text(0.0, 0.98, 'The two-bit counter counts mod 4 under the real Han oscillator',
             transform=axN.transAxes, ha='left', va='top', fontsize=15.5, color=INK,
             fontweight='bold')
    status = ('Certified' if meta['certified'] else 'Not certified')
    axN.text(0.0, 0.40,
             f"{status}  |  {meta['cold_start']['reads']} read windows over "
             f"{meta['hours']:g} h  |  sequence {seq[:12]}...  |  "
             f"bit 1 setup {meta['margins']['bit1']['min_setup_h']:.2f} h / "
             f"hold {meta['margins']['bit1']['min_hold_h']:.2f} h  |  "
             f"minimum commitment {meta['cold_start']['minimum_commitment']:.2f}",
             transform=axN.transAxes, ha='left', va='top', fontsize=10.5,
             color=TEAL if meta['certified'] else RISK_ORANGE, fontweight='bold')

    footer = (
        'Model: 34-state two-bit counter (model_twobit34.py), deterministic ODE, '
        'no noise term\n'
        'Input: Han oscillator C31 translation flux -- the real upstream, coupled '
        'one-way at the unloaded limit (rho = 1)\n'
        f"Parameters: uM_per_au = {wp['uM_per_au']:g}, carry A/F maturation = "
        f"{wp['carry_maturation_half_life_min']:g} min, carry mRNA = "
        f"{wp['carry_mrna_half_life_min']:g} min, clock = "
        f"{wp['clock_K_au']:g} / {wp['clock_n']:.1f}, extra dilution term = 0 "
        f"(selected_profile.json)")
    axF.text(0.0, 1.0, footer, transform=axF.transAxes, ha='left', va='top',
             fontsize=9.2, color='#3d5a63', linespacing=1.6)

    return fig, dict(axA=axA, axA2=axA2, axB=axB, axC=axC, axN=axN, axF=axF,
                     df=df, meta=meta, values=values, troughs=troughs, edges=edges)


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
    # Panel titles and the footer are single multi-line text objects, so labels are
    # matched as substrings, not as whole strings.
    joined_texts = '\n'.join(texts)
    footer_lines = [ln.strip() for t in texts for ln in t.splitlines()]
    df, meta = S['df'], S['meta']

    print('-- 1. required labels, verbatim')
    for lab in REQUIRED_LABELS:
        check(lab in joined_texts, f'label present: {lab!r}')

    print('-- 2. curves are the raw CSV columns')
    for col, label, ax in (('b0_S', 'bit 0 (S0)', S['axB']),
                           ('b1_S', 'bit 1 (S1)', S['axB'])):
        want = df[col].to_numpy(dtype=float)
        hits = [ln for ln in ax.get_lines()
                if np.array_equal(np.asarray(ln.get_ydata(), dtype=float), want)]
        check(len(hits) == 1, f'{label} is the raw {col} column ({len(hits)} match)')
    for col, label in OSCILLATOR:
        want = normalised(df[col])
        hits = [ln for ln in S['axA'].get_lines()
                if np.array_equal(np.asarray(ln.get_ydata(), dtype=float), want)]
        check(len(hits) == 1,
              f'upstream trace {label} is the normalized raw {col} column '
              f'({len(hits)} match)')
    flux = df.C31_flux_au_h.to_numpy(dtype=float)
    hits = [ln for ln in S['axA2'].get_lines()
            if np.array_equal(np.asarray(ln.get_ydata(), dtype=float), flux)]
    check(len(hits) == 1, 'the flux trace is the raw C31_flux_au_h column')

    print('-- 3. the states are fractions and the panel keeps both bands')
    for col in ('b0_S', 'b1_S'):
        v = df[col].to_numpy(dtype=float)
        check(v.min() >= -1e-12 and v.max() <= 1.0 + 1e-9,
              f'{col} stays inside [0, 1] (min {v.min():.4f}, max {v.max():.4f})')
        check(v.min() < BAND_LOW and v.max() > BAND_HIGH,
              f'{col} visits both the low and the high band '
              f'(min {v.min():.3f}, max {v.max():.3f})')

    print('-- 4. the decoded staircase is the artefact, unchanged')
    want_vals = np.asarray([np.nan if w['value'] is None else w['value']
                            for w in meta['read_windows']], dtype=float)
    got = np.asarray(S['values'], dtype=float)
    check(got.shape == want_vals.shape and np.array_equal(got, want_vals),
          f'the staircase carries all {len(want_vals)} read_windows.value entries '
          f'unchanged')
    check(np.all(np.isfinite(got)), 'no unlabelled window in this run')
    seq = meta['cold_start']['sequence']
    check(seq == EXPECTED_SEQUENCE,
          f'the sequence is seven whole mod-4 cycles ({seq!r})')
    check(np.all(np.diff(got.astype(int)) % 4 == 1),
          'every consecutive decoded value is +1 mod 4')

    print('-- 5. the working point and the provenance, not typed by hand')
    wp = meta['working_point']
    check(abs(wp['uM_per_au'] - 5.75) < 1e-12,
          f"uM_per_au is the selected profile's 5.75 (got {wp['uM_per_au']})")
    check(abs(wp['carry_maturation_half_life_min'] - 32.5) < 1e-12,
          f"carry A/F maturation is the selected profile's 32.5 min "
          f"(got {wp['carry_maturation_half_life_min']})")
    check(wp['source'].endswith('selected_profile.json'),
          'the provenance names selected_profile.json as the source of the point')
    g = meta.get('grid_row_cross_check') or {}
    check(g.get('matches') is True,
          'this run reproduces the 75-point robustness grid row for that point '
          f"(setup {g.get('bit1_setup_h')} vs {g.get('this_run_setup_h')})")

    print('-- 6. the three-element footer (model / input / parameters)')
    for prefix in FOOTER_PREFIXES:
        check(any(ln.startswith(prefix) for ln in footer_lines),
              f'the footer carries a {prefix!r} element')

    print('-- 7. no leakage claim is made')
    joined = joined_texts.lower()
    check('leak' not in joined and 'l_symmetric' not in joined,
          'the figure makes NO leakage claim (the two-bit certification has no leak arm)')

    print('-- 8. layout')
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
    pairs = [(placed[i][0].get_text()[:28], placed[j][0].get_text()[:28])
             for i in range(len(placed)) for j in range(i + 1, len(placed))
             if clash(placed[i][1], placed[j][1])]
    check(not pairs, f'no two labels overlap (found {len(pairs)}: {pairs[:4]})')
    figbb = fig.get_window_extent(renderer=renderer)
    outside = [t.get_text()[:30] for t, bb in placed
               if bb.x0 < figbb.x0 - 1 or bb.x1 > figbb.x1 + 1
               or bb.y0 < figbb.y0 - 1 or bb.y1 > figbb.y1 + 1]
    check(not outside, f'every figure label is inside the canvas (outside: {outside})')

    print('-- 9. the SVG carries the labels as real text')
    svg = svg_path.read_text(encoding='utf-8')
    check('<text' in svg, "SVG uses <text> (svg.fonttype = 'none')")
    for lab in REQUIRED_LABELS:
        esc = lab.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')
        check(esc in svg, f'SVG carries: {lab!r}')
    for prefix in FOOTER_PREFIXES:
        check(prefix in svg, f'SVG carries the footer element {prefix!r}')

    print()
    print(f'FAILURES: {len(fails)}')
    for f in fails:
        print('  -', f)
    return len(fails)


def main() -> int:
    try:
        sys.stdout.reconfigure(encoding='utf-8', errors='replace')
    except Exception:                                                   # noqa: BLE001
        pass
    plt.rcParams['svg.fonttype'] = SVG_FONTTYPE
    df, meta = load()
    print(f"run: {meta['status']}")
    print(f"working point: uM_per_au={meta['working_point']['uM_per_au']} "
          f"carry_mat={meta['working_point']['carry_maturation_half_life_min']} "
          f"certified={meta['certified']} reads={meta['cold_start']['reads']}")
    fig, S = figure(df, meta)
    OUT.mkdir(parents=True, exist_ok=True)
    png, svg, pdf = (OUT / f'{STEM}.png', OUT / f'{STEM}.svg', OUT / f'{STEM}.pdf')
    for p in (png, svg, pdf):
        fig.savefig(p, dpi=200, bbox_inches='tight', facecolor=fig.get_facecolor())
    fails = verify(fig, S, svg)
    plt.close(fig)
    for p in (png, svg, pdf):
        print(f'wrote {p}')
    return 1 if fails else 0


if __name__ == '__main__':
    sys.exit(main())
