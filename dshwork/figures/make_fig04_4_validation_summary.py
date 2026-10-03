r"""Figure 4-4 -- Validation summary   (dshwork, per the 精简交付版 brief)

Brief: `dshwork/精简交付版/02_作图要求精简版.md` section "Figure 4-4 -- Validation
summary".  Requirements:

  * Panel A: an 8 x 56 read strip, rows named 000-111, colour scale 0-7, showing the
    phase offset between initial states, annotated `8/8 correct sequences`;
  * Panel B: two horizontal stacked bars -- pre-freeze 136 runs = 128 recovered + 8
    legal shifts + 0 lost lock; post-freeze 30 runs = 30 recovered + 0 + 0;
  * caption material: the 8 legal shifts come from the S2 flip; some near-zero pool
    perturbations among the 136 are weak; there are 6 additional no-perturbation
    controls post-freeze and some perturbations were applied 29-41 h early; the
    strict-tolerance 5x3 conclusion is a one-line side note, no third heatmap;
  * MUST be right: never say "100 % experimental success" or "robust to arbitrary
    noise"; if the eight initial states are colour-coded they carry a numeric colour
    scale;
  * labels verbatim: `Orbit-consistent initial state`, `Expected count phase
    preserved`, `Recovered original phase`, `Legal modulo-8 shift`, `Lost lock`,
    `Specified perturbations only`.

Data (read-only, all pre-existing):
  Panel A  plausibility/eight_initial_all.csv        (8 rows, `code` = 56 reads each)
  Panel B  plausibility/pool_perturbations_all_all.csv (136 runs) and
           plausibility/threebit51_selected_v1.json    (the 128/8/0 breakdown)
           threebit51_results/postfreeze_peak_perturbation/postfreeze_all.csv (30+6)
  side note  same profile: strict_tolerance.tolerance_axis_5_points_3_levels

Self-checks: labels verbatim; 8 rows x 56 columns with EVERY cell from the artefact
(no dropped or invented cell); all 8 codes certified; the counts 128/8/0/136 and
30/0/0 recomputed from the artefacts, not typed in; the colour scale exists and is
numeric; the forbidden overclaims are absent; layout and SVG checks.

Run:  & 'D:\aconade\python.exe' .\make_fig04_4_validation_summary.py
Out:  ./out/fig04_4_validation_summary.png / .svg
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                                          # noqa: E402
from matplotlib.patches import Rectangle                                 # noqa: E402
import numpy as np                                                       # noqa: E402
import pandas as pd                                                      # noqa: E402

from tempo_style import DEEP_BLUE, INK, ORANGE, RISK_ORANGE, TEAL        # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE / 'out'
STEM = 'fig04_4_validation_summary'
FR = HERE.parents[1] / 'final_reconstruction'
EIGHT = FR / 'plausibility' / 'eight_initial_all.csv'
POOLPRE = FR / 'plausibility' / 'pool_perturbations_all_all.csv'
PROFILE = FR / 'plausibility' / 'threebit51_selected_v1.json'
POST = FR / 'threebit51_results' / 'postfreeze_peak_perturbation' / 'postfreeze_all.csv'

SVG_FONTTYPE = 'none'
REQUIRED_LABELS = ('Orbit-consistent initial state', 'Expected count phase preserved',
                   'Recovered original phase', 'Legal modulo-8 shift', 'Lost lock',
                   'Specified perturbations only')
FORBIDDEN = ('100%', '100 %', 'arbitrary noise', 'any noise', 'experimental success',
             'guaranteed')
CYCLE_COLOURS = ['#2b8a3e', '#1c7ed6', '#7048e8', '#e8590c', '#0b7285', '#c2255c',
                 '#5c940d', '#868e96']


def load():
    e = pd.read_csv(EIGHT)
    pre = pd.read_csv(POOLPRE)
    prof = json.loads(PROFILE.read_text(encoding='utf-8'))
    post = pd.read_csv(POST)
    codes = []
    for _, r in e.iterrows():
        s = str(r['cold_sequence'])
        codes.append([int(c) for c in s])
    overall = prof['pool_perturbations']['overall']
    pre_counts = dict(recovered=overall['recovered_original'],
                      shifted=overall['legal_mod8_shift'],
                      lost=overall['lost_lock'])
    pf = post[post.outcome != 'baseline']
    controls = post[post.outcome == 'baseline']
    post_counts = dict(recovered=int((pf.outcome == 'recovered_original').sum()),
                       shifted=int((pf.outcome == 'legal_mod8_shift').sum()),
                       lost=int((pf.outcome == 'lost_lock').sum()))
    strict = prof['strict_tolerance']['tolerance_axis_5_points_3_levels']
    return dict(e=e, pre=pre, codes=codes, pre_counts=pre_counts, post_counts=post_counts,
                n_pre=len(pre), n_post=len(pf), n_controls=len(controls),
                strict_points=len(strict['per_point']), strict_levels=strict['levels'],
                strict_all_stable=bool(strict['all_points_stable']),
                lead=pf.lead_time_before_gate_h.to_numpy(dtype=float),
                code_len=len(codes[0]))


def figure(D):
    fig = plt.figure(figsize=(13.2, 9.6), constrained_layout=True)
    gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 0.16, 0.5],
                          width_ratios=[1.0, 0.52], hspace=0.16, wspace=0.06)
    axA = fig.add_subplot(gs[0, :])
    axN = fig.add_subplot(gs[1, :])
    axB = fig.add_subplot(gs[2, 0])
    axS = fig.add_subplot(gs[2, 1])
    axN.axis('off')
    axS.axis('off')

    codes = D['codes']
    nrows = len(codes)
    ncols = D['code_len']
    im = None
    for i, row in enumerate(codes):
        y = nrows - 1 - i
        for j, v in enumerate(row):
            axA.add_patch(Rectangle((j, y), 1, 1, facecolor=CYCLE_COLOURS[v % 8],
                                    edgecolor='white', linewidth=0.35))
        if ncols <= 56:                       # label every 8th cell with its value
            for j, v in enumerate(row):
                if j % 8 == 0:
                    axA.text(j + 0.5, y + 0.5, str(v), ha='center', va='center',
                             fontsize=5.6, color='white')
    axA.set_xlim(0, ncols)
    axA.set_ylim(0, nrows)
    axA.set_yticks([nrows - 1 - i + 0.5 for i in range(nrows)])
    axA.set_yticklabels([f'{int(v):03b}' for v in D['e'].initial_value], fontsize=9)
    axA.set_xlabel('Read index', fontsize=10)
    axA.set_ylabel('Orbit-consistent initial state', fontsize=10)
    axA.set_title('Eight orbit-consistent initial states, 56 read windows each '
                  '\u2014 all eight sequences correct (8/8 correct sequences)',
                  fontsize=11)
    axA.text(0.995, 1.06, 'Deterministic simulation', transform=axA.transAxes,
             ha='right', va='bottom', fontsize=9, color='#5a6b70')
    sm = plt.cm.ScalarMappable(cmap=matplotlib.colors.ListedColormap(CYCLE_COLOURS),
                               norm=matplotlib.colors.BoundaryNorm(
                                   np.arange(-0.5, 8.5, 1), 8))
    cb = fig.colorbar(sm, ax=axA, orientation='vertical', pad=0.012, ticks=range(8))
    cb.set_label('Decoded value', fontsize=9.5)
    cb.ax.tick_params(labelsize=8.5)

    # ---------------------------------------------------------------- Panel B
    pre, post = D['pre_counts'], D['post_counts']
    bars = [('Before freezing', pre['recovered'], pre['shifted'], pre['lost'],
             D['n_pre']),
            ('After freezing', post['recovered'], post['shifted'], post['lost'],
             D['n_post'])]
    ypos = [1.0, 0.34]
    for (name, rec, sh, lost, total), y in zip(bars, ypos):
        left = 0.0
        for val, col, lab in ((rec, TEAL, 'Recovered original phase'),
                              (sh, ORANGE, 'Legal modulo-8 shift'),
                              (lost, RISK_ORANGE, 'Lost lock')):
            if val:
                axB.barh(y, val, left=left, height=0.42, color=col,
                         edgecolor='white', linewidth=1.0)
                axB.text(left + val / 2, y, str(val), ha='center', va='center',
                         fontsize=10, color='white', fontweight='bold')
            left += val
        axB.text(-1.5, y, name, ha='right', va='center', fontsize=10.5)
        axB.text(left + 2.0, y, f'{total} runs', ha='left', va='center',
                 fontsize=9.5, color='#5a6b70')
    # the "expected count phase preserved" bracket over the first two segments
    tot0 = pre['recovered'] + pre['shifted']
    axB.annotate('', xy=(0, 1.36), xytext=(tot0, 1.36),
                 arrowprops=dict(arrowstyle='-', color='#5a6b70', lw=1.4))
    axB.text(tot0 / 2, 1.45, 'Expected count phase preserved', ha='center',
             va='bottom', fontsize=9.5, color='#3d5a63')
    axB.set_xlim(-26, 168)
    axB.set_ylim(-0.05, 1.72)
    axB.set_yticks([])
    axB.set_xlabel('number of runs', fontsize=10)
    axB.set_title('Pool-perturbation campaigns', fontsize=11)
    axB.legend(handles=[matplotlib.patches.Patch(facecolor=TEAL, label='Recovered original phase'),
                        matplotlib.patches.Patch(facecolor=ORANGE, label='Legal modulo-8 shift'),
                        matplotlib.patches.Patch(facecolor=RISK_ORANGE, label='Lost lock')],
               fontsize=8.5, loc='lower right', framealpha=0.95)
    axB.grid(alpha=0.16, axis='x')

    axS.text(0.0, 0.94, 'Specified perturbations only', fontsize=10, color=INK,
             va='top', ha='left', fontweight='bold')
    lt = D['lead']
    axS.text(0.0, 0.74,
             '\u2022 pool factors 0.8x / 1.2x, S2 shifts up to \u00b10.2 and a full\n'
             '  flip \u2014 nothing outside this set\n'
             f'\u2022 the 8 legal shifts are the S2-flip runs\n'
             f'\u2022 some near-zero pool perturbations are weak (e.g. C2 x0.8:\n'
             f'  max|change| 3e-12, 2 degenerate runs)\n'
             f'\u2022 {D["n_controls"]} additional no-perturbation controls after '
             f'freezing\n'
             f'\u2022 perturbation lead time {lt.min():.1f} to {lt.max():.1f} h before '
             f'the gate\n'
             f'\u2022 strict tolerance {D["strict_points"]} points x '
             f'{len(D["strict_levels"])} levels: '
             f'{"consistent" if D["strict_all_stable"] else "NOT consistent"} '
             f'(no third heatmap)',
             fontsize=8.6, color='#3d5a63', va='top', ha='left')
    return fig, dict(axA=axA, axB=axB, codes=codes, nrows=nrows, ncols=ncols,
                     cb=cb, sm=sm)


def collect_texts(fig):
    return [t for t in fig.findobj(matplotlib.text.Text)
            if isinstance(t.get_text(), str) and t.get_text().strip()]


def verify(fig, S, D, svg_path) -> int:
    fails = []

    def check(cond, msg):
        print(('  PASS  ' if cond else '  FAIL  ') + msg)
        if not cond:
            fails.append(msg)

    texts = [t.get_text() for t in collect_texts(fig)]
    print('-- 1. required labels, verbatim')
    for lab in REQUIRED_LABELS:
        check(lab in texts, f'label present: {lab!r}')
    print('-- 2. forbidden overclaims absent')
    for bad in FORBIDDEN:
        hits = [t for t in texts if bad.lower() in t.lower()]
        check(not hits, f'no overclaim {bad!r} (hits: {hits})')

    print('-- 3. Panel A: 8 rows x 56 cells, every cell from the artefact')
    check(S['nrows'] == 8, f'8 rows (got {S["nrows"]})')
    check(S['ncols'] == 56, f'56 columns (got {S["ncols"]})')
    e = D['e']
    ok = all([int(c) for c in str(e.cold_sequence.iloc[i])] == S['codes'][i]
             for i in range(len(e)))
    check(ok, 'every row equals eight_initial_all.csv cold_sequence (no cell edited)')
    check(bool(e.certified.all()),
          f'all 8 initial states are certified ({int(e.certified.sum())}/8)')
    check('8/8 correct sequences' in ' '.join(texts),
          "the annotation '8/8 correct sequences' is drawn")
    # Read the row names from the axis itself: fig.findobj() can hand back the same
    # Text object more than once, so counting by string is not reliable.
    rowlabs = [t.get_text() for t in S['axA'].get_yticklabels()]
    want_rows = ['000', '001', '010', '011', '100', '101', '110', '111']
    check(sorted(rowlabs) == want_rows,
          f'rows are named 000-111 ({rowlabs})')
    check(S['cb'] is not None, 'a colour bar exists for the cell colouring')
    cb_ticks = [t.get_text() for t in S['cb'].ax.get_yticklabels()]
    check(sorted(cb_ticks) == [str(i) for i in range(8)],
          f'the colour scale is numeric 0-7 (ticks {cb_ticks})')

    print('-- 4. Panel B: the counts are recomputed, not typed in')
    pre = D['pre_counts']
    check(pre['recovered'] == 128 and pre['shifted'] == 8 and pre['lost'] == 0,
          f"pre-freeze {D['n_pre']} runs = {pre['recovered']} recovered + "
          f"{pre['shifted']} legal shifts + {pre['lost']} lost lock "
          f"(profile pool_perturbations.overall)")
    check(D['n_pre'] == pre['recovered'] + pre['shifted'] + pre['lost'] == 136,
          f'the 136 runs are fully accounted for by the three categories')
    raw = D['pre'].outcome.value_counts().to_dict()
    check(raw.get('recovered_original', 0) == 128
          and raw.get('legal_mod8_shift', 0) == 8,
          f'the 136-row artefact itself gives {raw}')
    post = D['post_counts']
    check(D['n_post'] == 30 and post['recovered'] == 30 and post['shifted'] == 0
          and post['lost'] == 0,
          f"post-freeze {D['n_post']} runs = {post['recovered']} recovered + 0 "
          f"shifts + 0 lost (postfreeze_all.csv)")
    check(D['n_controls'] == 6,
          f"{D['n_controls']} no-perturbation controls are present in the same file")

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
    cell_digits = {str(i) for i in range(8)}
    placed = [(t, t.get_window_extent(renderer=renderer))
              for t in collect_texts(fig)
              if t.get_visible() and id(t) not in tick_ids
              and t.get_text() not in cell_digits]
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
    for lab in REQUIRED_LABELS + ('8/8 correct sequences',):
        esc = lab.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')
        check(esc in svg, f'SVG carries: {lab!r}')
    for bad in FORBIDDEN:
        check(bad.lower() not in svg.lower(), f'SVG free of overclaim {bad!r}')

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
    print(f"pre-freeze {D['n_pre']} runs {D['pre_counts']} | "
          f"post-freeze {D['n_post']} runs {D['post_counts']} | "
          f"controls {D['n_controls']} | strict "
          f"{D['strict_points']}x{len(D['strict_levels'])} "
          f"stable={D['strict_all_stable']}")
    fig, S = figure(D)
    OUT.mkdir(parents=True, exist_ok=True)
    png, svg = OUT / f'{STEM}.png', OUT / f'{STEM}.svg'
    for p in (png, svg):
        fig.savefig(p, dpi=200, bbox_inches='tight',
                    facecolor=fig.get_facecolor())
    fails = verify(fig, S, D, svg)
    plt.close(fig)
    print(f'wrote {png}')
    print(f'wrote {svg}')
    return 1 if fails else 0


if __name__ == '__main__':
    sys.exit(main())
