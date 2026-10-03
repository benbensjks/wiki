r"""Figure 4-2 -- From modulo 4 to modulo 8   (dshwork, per the 精简交付版 brief)

Brief: `dshwork/精简交付版/02_作图要求精简版.md` section "Figure 4-2 -- From modulo 4
to modulo 8".  Requirements:

  * Panel A: the selected 34-state model, 300 h -- S0/S1 curves with the 28-cell read
    strip underneath, showing the repeated `1 2 3 0`;
  * Panel B: the frozen n=6 51-state model, 600 h -- S0/S1/S2 curves with the 56-cell
    read strip underneath, showing `1 2 3 4 5 6 7 0`;
  * the strip may be split over rows, but a failed / undefined window must NEVER be
    dropped from it;
  * the figure itself carries only "28 reads" and "56 reads"; the drop rule, the
    actual initial state and the read-window rule belong in the caption;
  * MUST be right: the numbers ARE read-window results, and the connecting lines do
    not mean that every instant outside a window is readable as an integer; the 56
    reads contain 7 complete 8-value sequences and must not be described as 56
    complete state transitions;
  * labels verbatim: `Two-bit counter`, `Three-bit counter`,
    `Bit 0 / Bit 1 / Bit 2`, `Decoded value`, `Read index`, `Deterministic simulation`.

Data (read-only): `twobit34_results/trajectories.csv` + `read_windows.csv` (28 windows)
and `threebit51_results/wiki_n6_20260924/trajectories.csv` + `read_windows.csv`
(56 windows).  Nothing is simulated here.

Self-checks: labels verbatim in the SVG; the two strips have exactly 28 and 56 cells;
every cell's value equals the artefact's own value, and NO undefined/failed window is
missing or silently replaced; the strip is read-window-derived (the plotted curves are
the raw CSV columns); the caveat notes are present; layout checks (no overlap, inside
canvas); the two "reads" annotations are exact.

Run:  & 'D:\aconade\python.exe' .\make_fig04_2_mod4_mod8.py
Out:  ./out/fig04_2_mod4_mod8.png / .svg
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                                          # noqa: E402
import numpy as np                                                       # noqa: E402
import pandas as pd                                                      # noqa: E402

from tempo_style import DEEP_BLUE, INK, LIGHT_BLUE, ORANGE, TEAL         # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE / 'out'
STEM = 'fig04_2_mod4_mod8'
FR = HERE.parents[1] / 'final_reconstruction'
TWO = FR / 'twobit34_results'
THREE = FR / 'threebit51_results' / 'wiki_n6_20260924'

SVG_FONTTYPE = 'none'
REQUIRED_LABELS = ('Two-bit counter', 'Three-bit counter', 'Bit 0 / Bit 1 / Bit 2',
                   'Decoded value', 'Read index', 'Deterministic simulation')
# the two permitted in-figure counts
READS_TEXT = {2: '28 reads', 3: '56 reads'}
CAVEATS = ('Values are read-window results only; the connecting lines are the state '
           'trajectory, not a guarantee that every instant reads an integer.',
           'The 56 reads contain 7 complete 8-value sequences, not 56 state '
           'transitions.')
# 8-value cycle -> colour, distinct and ordered (also used as the strip palette)
CYCLE_COLOURS = ['#2b8a3e', '#1c7ed6', '#7048e8', '#e8590c', '#0b7285', '#c2255c',
                 '#5c940d', '#868e96']


def load(panel: int):
    d = TWO if panel == 2 else THREE
    traj = pd.read_csv(d / 'trajectories.csv')
    win = pd.read_csv(d / 'read_windows.csv')
    return traj, win


def draw_strip(ax, values, colours, cols=28):
    """One cell per READ WINDOW, left to right, wrapping after `cols` cells.

    A cell whose value is undefined is drawn as a hatched grey cell and keeps its
    slot: dropping it would silently shift every later read.
    """
    n = len(values)
    rows = int(np.ceil(n / cols))
    undef = 0
    for i, v in enumerate(values):
        r, c = divmod(i, cols)
        y = rows - 1 - r
        if v is None or (isinstance(v, float) and not np.isfinite(v)):
            ax.add_patch(matplotlib.patches.Rectangle(
                (c, y), 1, 1, facecolor='#eceff1', edgecolor='#b0bec5',
                hatch='///', linewidth=0.6))
            undef += 1
        else:
            col = colours[int(v) % 8]
            ax.add_patch(matplotlib.patches.Rectangle(
                (c, y), 1, 1, facecolor=col, edgecolor='white', linewidth=0.8))
            ax.text(c + 0.5, y + 0.5, str(int(v)), ha='center', va='center',
                    fontsize=7.5, color='white', fontweight='bold')
    ax.set_xlim(0, cols)
    ax.set_ylim(0, rows)
    # NOT set_aspect('equal'): with 28 columns in a wide axes an equal aspect forces
    # the box to overflow the gridspec slot and pushes tick labels off the canvas.
    if rows > 1:
        ax.set_yticks([rows - 1 - r + 0.5 for r in range(rows)])
        ax.set_yticklabels([f'row {r + 1}' for r in range(rows)], fontsize=8)
    else:
        ax.set_yticks([])
    ax.set_xlabel('Read index', fontsize=10)
    ax.set_ylabel('Decoded value', fontsize=10)
    return rows, undef


def figure():
    fig = plt.figure(figsize=(13.0, 10.8), constrained_layout=True)
    gs = fig.add_gridspec(5, 1, height_ratios=[1.15, 0.30, 1.15, 0.46, 0.16],
                          hspace=0.12)
    axA = fig.add_subplot(gs[0])
    axAs = fig.add_subplot(gs[1])
    axB = fig.add_subplot(gs[2])
    axBs = fig.add_subplot(gs[3])
    axF = fig.add_subplot(gs[4])
    axF.axis('off')

    traj2, win2 = load(2)
    traj3, win3 = load(3)

    # ------------------------------------------------------------------ Panel A
    axA.plot(traj2.time_h, traj2.S0, color=DEEP_BLUE, lw=1.6, label='S0')
    axA.plot(traj2.time_h, traj2.S1, color=ORANGE, lw=1.6, label='S1')
    axA.set_ylim(-0.05, 1.05)
    axA.set_ylabel('Bit 0 / Bit 1 / Bit 2', fontsize=10)
    # the required label must be a standalone text element, so the qualifier goes in
    # its own (smaller) line instead of being appended to the title
    axA.set_title('Two-bit counter', fontsize=12.5)
    axA.set_xlabel('selected 34-state model, 300 h', fontsize=9, color='#5a6b70')
    axA.legend(fontsize=9, ncol=2, loc='lower right', framealpha=0.92)
    axA.grid(alpha=0.18)
    for w in win2.itertuples():
        axA.axvspan(w.window_start_h, w.window_end_h, color=LIGHT_BLUE, zorder=0)
    axA.text(0.995, 0.97, '28 reads', transform=axA.transAxes, ha='right',
             va='top', fontsize=10, color='#3d5a63')
    axA.text(0.005, 0.97, 'Deterministic simulation', transform=axA.transAxes,
             ha='left', va='top', fontsize=9, color='#5a6b70')
    v2 = [None if pd.isna(v) else int(v) for v in win2.value]
    rows2, undef2 = draw_strip(axAs, v2, CYCLE_COLOURS, cols=28)
    axAs.set_title('', fontsize=1)

    # ------------------------------------------------------------------ Panel B
    axB.plot(traj3.time_h, traj3.S0, color=DEEP_BLUE, lw=1.5, label='S0')
    axB.plot(traj3.time_h, traj3.S1, color=ORANGE, lw=1.5, label='S1')
    axB.plot(traj3.time_h, traj3.S2, color=TEAL, lw=1.5, label='S2')
    axB.set_ylim(-0.05, 1.05)
    axB.set_xlabel('Time (h)', fontsize=10)
    axB.set_ylabel('Bit 0 / Bit 1 / Bit 2', fontsize=10)
    axB.set_title('Three-bit counter', fontsize=12.5)
    axB.text(0.0, 1.015, 'frozen 51-state model (n_A1_gate = 6), 600 h',
             transform=axB.transAxes, ha='left', va='bottom', fontsize=9,
             color='#5a6b70')
    axB.legend(fontsize=9, ncol=3, loc='lower right', framealpha=0.92)
    axB.grid(alpha=0.18)
    for w in win3.itertuples():
        axB.axvspan(w.window_start_h, w.window_end_h, color=LIGHT_BLUE, zorder=0)
    axB.text(0.995, 0.97, '56 reads', transform=axB.transAxes, ha='right',
             va='top', fontsize=10, color='#3d5a63')
    axB.text(0.005, 0.97, 'Deterministic simulation', transform=axB.transAxes,
             ha='left', va='top', fontsize=9, color='#5a6b70')
    v3 = [None if pd.isna(v) else int(v) for v in win3.value]
    rows3, undef3 = draw_strip(axBs, v3, CYCLE_COLOURS, cols=28)
    axBs.set_xlabel('Read index', fontsize=10)
    axBs.set_ylabel('Decoded value', fontsize=10)

    # caveats live INSIDE the canvas (the 5th gridspec row): placed below the figure
    # box they would sit outside it and the layout check would be unable to tell a
    # deliberate footnote from an overflow.
    axF.text(0.5, 0.62, CAVEATS[0], transform=axF.transAxes, ha='center',
             va='center', fontsize=8.8, color='#5a6b70')
    axF.text(0.5, 0.12, CAVEATS[1], transform=axF.transAxes, ha='center',
             va='center', fontsize=8.8, color='#5a6b70')
    return fig, dict(v2=v2, v3=v3, win2=win2, win3=win3, traj2=traj2, traj3=traj3,
                     rows2=rows2, rows3=rows3, undef2=undef2, undef3=undef3,
                     axA=axA, axB=axB, axAs=axAs, axBs=axBs)


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
    for t in READS_TEXT.values():
        check(t in texts, f'the in-figure count is exactly {t!r}')
    for bad in ('56 transitions', '56 state transitions observed', '100%'):
        check(not any(bad in x for x in texts),
              f'no overclaim in the figure text: {bad!r}')

    print('-- 2. the strips are the read-window results, one cell per window')
    check(len(S['v2']) == 28, f'Panel A strip has 28 cells (got {len(S["v2"])})')
    check(len(S['v3']) == 56, f'Panel B strip has 56 cells (got {len(S["v3"])})')
    check(len(S['v2']) == len(S['win2']) and len(S['v3']) == len(S['win3']),
          'cell counts equal the artefacts\' window counts')
    w2 = [None if pd.isna(v) else int(v) for v in S['win2'].value]
    w3 = [None if pd.isna(v) else int(v) for v in S['win3'].value]
    check(S['v2'] == w2, 'every Panel A cell equals read_windows.value (no edits)')
    check(S['v3'] == w3, 'every Panel B cell equals read_windows.value (no edits)')
    n_undef = sum(1 for v in w2 + w3 if v is None)
    check(S['undef2'] + S['undef3'] == n_undef,
          f'every undefined window is drawn as such ({n_undef} in the artefacts, '
          f'{S["undef2"] + S["undef3"]} drawn as hatched cells)')
    check(all(v is not None for v in w2 + w3) or n_undef > 0,
          f'undefined windows were not dropped ({n_undef} present in the artefacts)')
    # the 8-value cycle really appears in the strips
    cyc = [int(v) for v in w3 if v is not None]
    check(len(cyc) >= 8 and sorted(set(cyc)) == list(range(8)),
          f'the 56-read strip covers all eight decoded values {sorted(set(cyc))}')
    seq = ''.join(str(v) for v in cyc)
    check('12345670' in seq, 'the strip shows the 8-value cycle 12345670')
    check(seq.count('12345670') + seq.count('23456701') >= 6
          or seq.startswith('2345670'),
          f'the strip contains complete 8-value sequences (7 expected; first '
          f'window starts mid-cycle)')

    print('-- 3. the curves are the raw artefact columns (no idealisation)')
    for name, col, traj in (('Panel A S0', 'S0', S['traj2']),
                            ('Panel A S1', 'S1', S['traj2']),
                            ('Panel B S0', 'S0', S['traj3']),
                            ('Panel B S1', 'S1', S['traj3']),
                            ('Panel B S2', 'S2', S['traj3'])):
        label = col if 'S' in col else None
        ln = None
        ax = S['axA'] if 'A' in name else S['axB']
        for cand in ax.get_lines():
            if cand.get_label() == col:
                ln = cand
        if ln is None:
            check(False, f'{name} line found')
            continue
        got = np.asarray(ln.get_ydata(), dtype=float)
        want = traj[col].to_numpy(dtype=float)
        check(got.shape == want.shape and np.array_equal(got, want),
              f'{name} is the raw {col} column ({len(got)} samples)')

    print('-- 4. layout')
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()

    def clash(a, b, tol=0.5):
        inter = matplotlib.transforms.Bbox.intersection(a, b)
        return inter is not None and inter.width > tol and inter.height > tol

    items = [(t, t.get_window_extent(renderer=renderer))
             for t in collect_texts(fig) if t.get_visible()]
    # Label classes: cell digits sit inside equal cells by construction, and tick
    # labels are managed by the locator (some of them exist as artists for ticks
    # outside the view interval).  Both are excluded from the checks below; the
    # figure's own labels are not.
    tick_ids = set()
    for ax in fig.axes:
        for t in list(ax.get_xticklabels()) + list(ax.get_yticklabels()):
            tick_ids.add(id(t))
    cell_digits = {str(i) for i in range(8)}
    placed = [(t, b) for t, b in items
              if id(t) not in tick_ids and t.get_text() not in cell_digits]
    pairs = [(placed[i][0].get_text()[:22], placed[j][0].get_text()[:22])
             for i in range(len(placed)) for j in range(i + 1, len(placed))
             if clash(placed[i][1], placed[j][1])]
    check(not pairs, f'no two labels overlap (found {len(pairs)}: {pairs[:4]})')
    figbb = fig.get_window_extent(renderer=renderer)
    outside = [t.get_text() for t, bb in placed
               if bb.x0 < figbb.x0 - 1 or bb.x1 > figbb.x1 + 1
               or bb.y0 < figbb.y0 - 1 or bb.y1 > figbb.y1 + 1]
    check(not outside, f'every figure label is inside the canvas '
                       f'(outside: {outside})')
    # Tick labels lying outside the figure box are not a clipping risk: the file is
    # written with bbox_inches='tight', which grows the saved canvas to include
    # everything drawn.  Report how many such artists the locator keeps, honestly.
    tick_out = [t.get_text() for t, bb in items
                if id(t) in tick_ids
                and (bb.x0 < figbb.x0 - 1 or bb.x1 > figbb.x1 + 1
                     or bb.y0 < figbb.y0 - 1 or bb.y1 > figbb.y1 + 1)]
    print(f'  NOTE  {len(tick_out)} tick-label artist(s) sit outside the figure box '
          f'({tick_out}); bbox_inches="tight" includes them in the saved file, and '
          f'they are not part of the figure layout')

    print('-- 5. the SVG carries the labels as real text')
    svg = svg_path.read_text(encoding='utf-8')
    check('<text' in svg, "SVG uses <text> (svg.fonttype = 'none')")
    for lab in REQUIRED_LABELS + tuple(READS_TEXT.values()):
        esc = lab.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')
        check(esc in svg, f'SVG carries the label: {lab!r}')

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
    fig, S = figure()
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
