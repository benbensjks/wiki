"""M-18: two-terminal clock working-region map, drawn as a DISCRETE parameter map.

Data: `results/clock_axis_20261001_180035_369646/all.csv` (60 rows =
10 clock_K x 2 Hill n x 3 arms), 300 h per point, Han upstream.

Three variants, so the two design choices can be compared directly:

    fig_M18_map_3cell           3 criteria cells per point (counting / events /
                                certified), single panel  -- "strict three cells"
    fig_M18_map_2cell           2 criteria cells (counting / events), single
                                panel  -- isolates the cell count
    fig_M18_map_2cell_failmode  2 cells + a failure-mode heatmap panel
                                -- the default proposal

Design rules (all deliberate, see the manifest entry):
  * x is an ORDINAL axis: one equally spaced column per MEASURED clock_K.  The
    measured levels are far from uniform (0.075 -> 0.15 spans 0.075 while
    0.20 -> 0.30 spans 0.02), so a linear axis would crush the boundary points
    together.  The cost is that column spacing is NOT physical distance -- the
    caption must say so.
  * a top secondary axis carries r = clock_K / Int0_peak (0.5484815380117275),
    the only cross-terminal comparable number.
  * every point is a stacked glyph; NOTHING is drawn between measured points, so
    no colour band can ever imply an unscanned interval.
  * pass/fail is carried by BOTH colour and a glyph (sqrt / cross), per the
    project's "do not rely on colour alone" rule.
  * the terminal boundary is drawn as a BRACKET over the unscanned gap between
    the last passing and the first failing measured level, labelled as an
    interval -- the positive way to say "this range was not scanned".

`audit()` measures every text bbox against every other text bbox and against
every map cell in display space after a draw, because I cannot look at the
figure; it is the numerical substitute for visual inspection.

Read-only with respect to every existing artefact.
"""
from __future__ import annotations

import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Rectangle
from matplotlib.transforms import Bbox

# Without this the SVG backend turns every glyph into a path, the figure becomes
# unsearchable and the layout check sees 0 <text> nodes.
matplotlib.rcParams['svg.fonttype'] = 'none'

sys.dont_write_bytecode = True

HERE = Path(__file__).resolve().parent
RES = HERE / 'results' / 'clock_axis_20261001_180035_369646'
OUT = HERE / 'out'
SRC = RES / 'all.csv'
INT0_PEAK = 0.5484815380117275
EXPECTED = '1234567012345670123'
OKC, BADC = '#2e7d32', '#c62828'
GRAY, INKC = '#8ca0a5', '#304B53'
ARM_SHORT = {('hzh', 'base'): 'HZH  (HBY terminal)',
             ('hzz', 'off'): 'HZZ  self-fb OFF',
             ('hzz', 'on'): 'HZZ  self-fb ON'}
MODE_TXT = {'s': 'every labelled window reads correctly; only unlabelled '
                 'windows remain',
            'r': 'the x-free core is a fixed 4-state cycle',
            'u': 'otherwise -- largely unlabelled / non-discriminating'}
# geometry of the map, in cell units
CELL_H, COL_W, LANE_GAP, ARM_GAP = 1.0, 0.84, 0.34, 1.90


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest().upper()


def failure_mode(seq: str) -> str:
    """Classify a failing steady sequence; the rule is printed in the caption so
    a reader can re-derive it from the raw strings."""
    if seq == EXPECTED:
        return ''                      # not a failure
    nx = seq.count('x')
    if nx and all(seq[i] == EXPECTED[i] for i in range(len(seq)) if seq[i] != 'x'):
        return 's'
    if nx > len(seq) // 2:
        return 'u'
    core = seq.replace('x', '')
    if len(core) >= 8 and all(core[i] == core[i % 4] for i in range(len(core))):
        return 'r'
    return 'u'


def load():
    with SRC.open(encoding='utf-8') as fh:
        rows = list(csv.DictReader(fh))
    Ks = sorted({float(r['clock_K']) for r in rows})
    pts = {}
    for r in rows:
        key = (r['kind'], r['arm'], round(float(r['n']), 1), float(r['clock_K']))
        cnt = r['counting_passed'] == 'True'
        ev = r['event_causality_passed'] == 'True'
        cer = r['certified_v1'] == 'True'
        assert cer == (cnt and ev), key      # asserted, not assumed
        pts[key] = dict(counting=cnt, events=ev, certified=cer,
                        mode=failure_mode(r['steady_sequence']),
                        sequence=r['steady_sequence'],
                        reads=int(r['steady_reads']),
                        margin=float(r['global_min_timing_margin_h']))
    return rows, Ks, pts


def draw_map(ax, Ks, pts, ncells):
    """ncells = 3 -> (counting, events, certified); 2 -> (counting, events)."""
    fields = [('counting', 'counting'), ('events', 'events'),
              ('certified', 'certified')][:ncells]
    glyphs = []                       # (data x, data y, char) for the SVG check
    y, bands, cells = 0.0, [], []
    arms = [('hzh', 'base'), ('hzz', 'off'), ('hzz', 'on')]
    for ai, (kind, arm) in enumerate(arms):
        y0 = y
        for n in (2.0, 3.0):
            for i, (f, _) in enumerate(fields):
                for j, K in enumerate(Ks):
                    p = pts[(kind, arm, n, K)]
                    ok = p[f]
                    r = Rectangle((j - COL_W / 2, y - (i + 1) * CELL_H), COL_W,
                                  CELL_H * 0.90, facecolor=OKC if ok else BADC,
                                  edgecolor='white', lw=1.2, zorder=3)
                    ax.add_patch(r)
                    cells.append(r)
                    glyph = '√' if ok else ('×' if f != 'counting'
                                            else (p['mode'] or '×'))
                    ax.text(j, y - (i + 0.5) * CELL_H, glyph, ha='center',
                            va='center', fontsize=12, color='white',
                            fontweight='bold', zorder=4)
            ax.text(-0.55, y - ncells * CELL_H / 2, 'n=%d' % int(n), ha='right',
                    va='center', fontsize=10, color=INKC)
            y -= ncells * CELL_H + LANE_GAP
        bands.append((y0, y + LANE_GAP, kind, arm))
        if ai < len(arms) - 1:
            ax.plot([-1.45, len(Ks) - 0.38], [y - LANE_GAP - 0.06] * 2,
                    color=GRAY, ls=(0, (4, 3)), lw=1.0, zorder=1)
            y -= ARM_GAP
    return bands, y


def mark_boundary(ax, Ks, pts, kind, arm, band_end):
    """Bracket the UNSCANNED gap between last-pass and first-fail measured level."""
    per_col = [all(pts[(kind, arm, n, Ks[j])]['counting'] for n in (2.0, 3.0))
               for j in range(len(Ks))]
    last_ok = max((j for j, v in enumerate(per_col) if v), default=None)
    first_bad = next((j for j, v in enumerate(per_col) if not v), None)
    if last_ok is None or first_bad is None:
        return None
    yb = band_end - 0.10
    ax.plot([last_ok, first_bad], [yb, yb], color=INKC, lw=2.0, zorder=5)
    for x in (last_ok, first_bad):
        ax.plot([x, x], [yb - 0.16, yb + 0.16], color=INKC, lw=2.0, zorder=5)
    ax.text((last_ok + first_bad) / 2, yb - 0.24,
            'terminal boundary in (%.2f, %.2f]  —  NOT scanned'
            % (Ks[last_ok], Ks[first_bad]), ha='center', va='top', fontsize=9.5,
            color=INKC, fontweight='bold', zorder=5)
    return dict(last_passing_K=Ks[last_ok], first_failing_K=Ks[first_bad],
                interval=[Ks[last_ok], Ks[first_bad]],
                r_last=Ks[last_ok] / INT0_PEAK, r_first=Ks[first_bad] / INT0_PEAK)


def finish_map(ax, Ks, ybot, title, notes):
    ax.set_xlim(-1.45, len(Ks) - 0.38)
    ax.set_ylim(ybot - 1.70, 1.35)
    ax.set_yticks([])
    for s in ('top', 'right', 'left'):
        ax.spines[s].set_visible(False)
    ax.spines['bottom'].set_color(GRAY)
    ax.set_xticks(range(len(Ks)))
    ax.set_xticklabels(['%.3f' % K for K in Ks], fontsize=10)
    ax.set_xlabel('clock_K  (a.u.)   —  one equally spaced column per MEASURED '
                  'level; column spacing is NOT physical distance', fontsize=10.5)
    top = ax.secondary_xaxis('top')
    top.set_xticks(range(len(Ks)))
    top.set_xticklabels(['%.3f' % (K / INT0_PEAK) for K in Ks], fontsize=9)
    top.set_xlabel('r = clock_K / Int0_peak      (Int0_peak = %.6f)'
                   % INT0_PEAK, fontsize=10.5)
    ax.set_title(title, fontsize=13, pad=34)
    ax.tick_params(axis='both', length=0)
    for t in notes:
        t.set_clip_on(False)


def legend(fig, ncells, y=0.105, x=0.135):
    """Kept to short lines: a long single line runs off the canvas (the bbox
    audit measures this, my eye cannot)."""
    lines = ['√  pass   (green)          ×  fail   (red; the counting cell '
             'carries the failure-mode letter)',
             'certified ≡ counting ∧ events, verified 60/60 — the %s cell carries '
             'no independent information'
             % ('third' if ncells == 3 else 'omitted'),
             'failure modes:   s = ' + MODE_TXT['s'],
             '                          r = ' + MODE_TXT['r'],
             '                          u = ' + MODE_TXT['u']]
    for i, s in enumerate(lines):
        fig.text(x, y - i * 0.019, s, ha='left', fontsize=9.5, color=INKC)


def audit(fig, tag):
    """Display-space bbox audit: text vs text, and text vs map cells.

    A text whose centre falls inside a cell is that cell's own glyph (or a mode
    letter) and is expected there.  Everything else must be disjoint.
    """
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    texts = []
    for t in fig.findobj(matplotlib.text.Text):
        if not t.get_visible() or not t.get_text().strip():
            continue
        try:
            texts.append((t.get_text(), t.get_window_extent(rend)))
        except Exception:                                    # pragma: no cover
            pass
    rects = [p.get_window_extent(rend) for p in fig.findobj(Rectangle)
             if p.get_zorder() == 3]
    tt = []
    for i in range(len(texts)):
        for j in range(i + 1, len(texts)):
            ov = Bbox.intersection(texts[i][1], texts[j][1])
            if ov is not None and ov.width > 1.0 and ov.height > 1.0:
                tt.append((texts[i][0][:22], texts[j][0][:22],
                           round(ov.width, 1), round(ov.height, 1)))
    tc = []
    for txt, bb in texts:
        cx, cy = (bb.x0 + bb.x1) / 2, (bb.y0 + bb.y1) / 2
        for rb in rects:
            if rb.x0 <= cx <= rb.x1 and rb.y0 <= cy <= rb.y1:
                break
        else:
            for rb in rects:
                ov = Bbox.intersection(bb, rb)
                if ov is not None and ov.width > 1.0 and ov.height > 1.0:
                    tc.append((txt[:22], round(ov.width, 1), round(ov.height, 1)))
                    break
    fb = fig.get_window_extent(rend)
    off = [(t[:22], round(max(fb.x0 - b.x0, b.x1 - fb.x1, fb.y0 - b.y0,
                              b.y1 - fb.y1), 1))
           for t, b in texts
           if b.x0 < fb.x0 - 1 or b.x1 > fb.x1 + 1 or b.y0 < fb.y0 - 1
           or b.y1 > fb.y1 + 1]
    print('  audit %-28s texts=%3d  TEXT-TEXT=%d  TEXT-CELL=%d  OUT=%d'
          % (tag, len(texts), len(tt), len(tc), len(off)), flush=True)
    for row in tt[:8]:
        print('      TT', row)
    for row in tc[:8]:
        print('      TC', row)
    for row in off[:8]:
        print('      OUT', row)
    return dict(texts=len(texts), text_text=tt, text_cell=tc, outside=off)


def draw_heat(axh, Ks, pts):
    fails = []
    for kind, arm in (('hzz', 'off'), ('hzz', 'on')):
        for n in (2.0, 3.0):
            for K in Ks:
                p = pts[(kind, arm, n, K)]
                if not p['counting']:
                    fails.append((kind, arm, n, K, p))
    fails.sort(key=lambda r: (r[0], r[1], r[3], r[2]))
    n = len(fails[0][4]['sequence'])
    M = np.full((len(fails), n), 8.0)                 # 8 = unlabelled
    for i, (_, _, _, _, p) in enumerate(fails):
        for j, ch in enumerate(p['sequence']):
            M[i, j] = 8.0 if ch == 'x' else float(ch)
    cmap = ListedColormap(list(plt.colormaps['viridis'](np.linspace(0, 1, 8)))
                          + ['#e9e9e9'])
    axh.imshow(M, cmap=cmap, vmin=-0.5, vmax=8.5, aspect='auto')
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if M[i, j] == 8.0:
                axh.text(j, i, 'x', ha='center', va='center', fontsize=7.5,
                         color=BADC, fontweight='bold')
    axh.set_yticks(range(len(fails)))
    axh.set_yticklabels(['%s   n=%d   K=%.3f   %s'
                         % (('OFF' if a == 'off' else 'ON '), int(nn), K,
                            pts[(k, a, nn, K)]['mode'])
                         for k, a, nn, K, _ in fails], fontsize=9)
    axh.set_xticks(range(0, n, 2))
    axh.set_xlabel('steady read-window index   —   grey "x" = window never '
                   'labelled, colour = decoded value', fontsize=10.5)
    axh.set_title('Failure modes of the %d failing points:  "s" leaves a clean '
                  'short tail, "r" repeats a 4-state cycle, "u" is largely '
                  'unlabelled' % len(fails), fontsize=11.5, pad=8)
    return fails


def build(name, ncells, with_heat, title, height, meta):
    OUT.mkdir(exist_ok=True)
    bounds, note = {}, []
    if with_heat:
        fig = plt.figure(figsize=(13.6, height))
        ax = fig.add_axes([0.135, 0.420, 0.815, 0.520])
        axh = fig.add_axes([0.135, 0.165, 0.800, 0.185])
        cax = fig.add_axes([0.955, 0.165, 0.013, 0.185])
    else:
        fig = plt.figure(figsize=(13.6, height))
        ax = fig.add_axes([0.135, 0.235, 0.815, 0.680])
        axh = cax = None
    bands, ybot = draw_map(ax, Ks, pts, ncells)
    fields = ['counting', 'events', 'certified'][:ncells]
    for kind, arm in (('hzz', 'off'), ('hzz', 'on')):
        b = mark_boundary(ax, Ks, pts, kind, arm,
                          next(x[1] for x in bands if x[2:] == (kind, arm)))
        if b:
            bounds['%s/%s' % (kind, arm)] = b
    y0, y1 = next((x[0], x[1]) for x in bands if x[2:] == ('hzh', 'base'))
    note.append(ax.text(4.5, y1 - 0.45,
                        'HZH passes at ALL %d measured levels (%.3f – %.3f):\n'
                        'neither boundary was reached' % (len(Ks), Ks[0], Ks[-1]),
                        ha='center', va='top', fontsize=10.5, color=OKC,
                        fontweight='bold', zorder=5))
    finish_map(ax, Ks, ybot, title, note)
    note.append(ax.text(-1.42, 1.02, 'each lane, top to bottom:   '
                        + '  /  '.join(fields),
                        ha='left', va='center', fontsize=10.5, color=INKC,
                        fontweight='bold', zorder=5))
    # arm labels live in the left figure margin: inside the axes they would be
    # pushed off the canvas (found by the bbox audit, not by eye)
    fh = fig.get_figheight() * fig.dpi
    for (kind, arm) in (('hzh', 'base'), ('hzz', 'off'), ('hzz', 'on')):
        y0, y1 = next((x[0], x[1]) for x in bands if x[2:] == (kind, arm))
        fig.text(0.128, ax.transData.transform((0, (y0 + y1) / 2))[1] / fh,
                 ARM_SHORT[(kind, arm)], ha='right', va='center', fontsize=11,
                 color=INKC, fontweight='bold')
    if axh is not None:
        fails = draw_heat(axh, Ks, pts)
        cb = fig.colorbar(axh.images[0], cax=cax, ticks=list(range(8)) + [8])
        cb.ax.set_yticklabels([str(i) for i in range(8)] + ['x'], fontsize=9)
        cb.outline.set_visible(False)
    legend(fig, ncells, y=0.098 if axh is not None else 0.135)
    a = audit(fig, name)
    meta['audit_' + name] = a
    for ext in ('png', 'svg'):
        fig.savefig(OUT / ('%s.%s' % (name, ext)), dpi=170, facecolor='white')
    plt.close(fig)
    print('  wrote %s.png' % name, flush=True)
    return bounds


def main():
    global Ks, pts
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8')
    rows, Ks, pts = load()
    print('loaded %d points; %d K levels' % (len(rows), len(Ks)), flush=True)
    meta = dict(source=str(SRC.relative_to(HERE)), source_sha256=sha(SRC),
                points=len(rows), K_levels=[float(k) for k in Ks],
                r_levels=[float(k / INT0_PEAK) for k in Ks],
                Int0_peak=INT0_PEAK, n_levels=[2.0, 3.0],
                arms=['hzh/base', 'hzz/off', 'hzz/on'], hours=300.0,
                certified_equals_counting_and_events=True,
                failure_mode_rule=('s: every labelled window equals the expected '
                                   'sequence (only unlabelled remain); r: the '
                                   'x-free core is a fixed 4-state cycle; '
                                   'u: otherwise'),
                note=('ordinal x axis -- column spacing is NOT physical distance; '
                      'only measured levels are drawn; boundaries are intervals, '
                      'not points'))
    meta['boundaries'] = build(
        'fig_M18_map_3cell', 3, False,
        'M-18  clock working region — discrete map, THREE criteria per point',
        9.6, meta)
    build('fig_M18_map_2cell', 2, False,
          'M-18  clock working region — discrete map, TWO criteria per point',
          9.6, meta)
    meta['boundaries_hzz_on'] = build(
        'fig_M18_map_2cell_failmode', 2, True,
        'M-18  clock working region — discrete map + failure-mode heatmap',
        14.6, meta)
    meta['outputs'] = {p.name: sha(p) for p in sorted(OUT.glob('fig_M18_*'))
                       if p.suffix in ('.png', '.svg')}
    (OUT / 'fig_M18_map_provenance.json').write_text(
        json.dumps(meta, ensure_ascii=False, indent=2), encoding='utf-8')
    print('wrote provenance', flush=True)


if __name__ == '__main__':
    main()
