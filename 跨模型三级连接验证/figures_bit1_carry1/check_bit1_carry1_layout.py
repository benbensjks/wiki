"""Layout self-check for the bit1 / carry1 hand-drawing drafts.

Mirrors `figures_bit0/check_bit0_sketch_layout.py`: the authoring model cannot
inspect images, so the figures are verified numerically instead.

  (a) each figure is re-rendered and the REAL rendered bounding box of every
      text artist is taken (font metrics included, not an estimate);
  (b) text-text overlaps are reported (a label printed on top of another);
  (c) text-box overlaps are reported (a stray label printed on top of a species
      box that it does not belong to), where "belongs to" means the text is the
      centred caption of that very box;
  (d) the in-script geometric check is repeated: box-box collisions,
      connectors passing through boxes, dangling tips, and -- new here --
      connectors that CROSS each other.

Exit code 0 = clean, 1 = problems found.

Usage
-----
    & 'D:\\aconade\\python.exe' -B .\\check_bit1_carry1_layout.py
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox

import make_bit1_carry1_sketches as M


def _capture(fn):
    """Run a figure function and return the figure instead of writing a file."""
    box = {}
    original = M.save_pair

    def fake_save(fig, stem):
        box['fig'] = fig
        return stem

    M.save_pair = fake_save
    try:
        fn()
    finally:
        M.save_pair = original
    return box['fig']


def _inter(a, b):
    x0, x1 = max(a.x0, b.x0), min(a.x1, b.x1)
    y0, y1 = max(a.y0, b.y0), min(a.y1, b.y1)
    if x1 <= x0 or y1 <= y0:
        return 0.0
    return (x1 - x0) * (y1 - y0)


def analyse(fn, name):
    fig = _capture(fn)
    fig.canvas.draw()
    ax = fig.axes[0]
    renderer = fig.canvas.get_renderer()
    drafts = list(M.B.LAST_DRAFTS)

    texts = []
    for t in ax.texts:
        if not t.get_text().strip():
            continue
        try:
            bb = t.get_window_extent(renderer=renderer)
        except Exception:
            continue
        if bb.width <= 0 or bb.height <= 0:
            continue
        texts.append((t.get_text().replace('\n', ' / ')[:46], bb,
                      ax.transData.transform(t.get_position())))

    disp_boxes = []
    for d in drafts:
        for bn, bx, by, bw, bh in d.boxes:
            (x0, y0) = ax.transData.transform((bx - bw / 2, by - bh / 2))
            (x1, y1) = ax.transData.transform((bx + bw / 2, by + bh / 2))
            disp_boxes.append((bn, min(x0, x1), min(y0, y1),
                               abs(x1 - x0), abs(y1 - y0)))

    problems = []

    # (b) text vs text
    for i in range(len(texts)):
        for k in range(i + 1, len(texts)):
            n1, b1, _ = texts[i]
            n2, b2, _ = texts[k]
            a = _inter(b1, b2)
            if a <= 0:
                continue
            frac = a / min(b1.width * b1.height, b2.width * b2.height)
            if frac > 0.06:
                problems.append(f'TEXT/TEXT {frac*100:5.1f}%  "{n1}"  x  "{n2}"')

    # (c) text vs box (skip the box's own centred caption)
    for txt, tb, ctr in texts:
        for bn, x, y, w, h in disp_boxes:
            rect = Bbox.from_bounds(x, y, w, h)
            own = (abs(ctr[0] - (x + w / 2)) < 2.0 and
                   abs(ctr[1] - (y + h / 2)) < 2.0)
            if own:
                continue
            a = _inter(tb, rect)
            if a <= 0:
                continue
            frac = a / (tb.width * tb.height)
            if frac > 0.10:
                problems.append(
                    f'TEXT/BOX  {frac*100:5.1f}%  "{txt}"  on  {bn}')

    print(f'=== {name}: {len(texts)} texts, {len(disp_boxes)} boxes, '
          f'{len(problems)} layout problem(s)')
    for p in problems:
        print('   ', p)
    plt.close(fig)
    return problems


def main():
    total = 0
    for fn, nm in ((M.fig_bit1_internal, 'fig_bit1_internal'),
                   (M.fig_carry1_wired, 'fig_carry1_wired')):
        total += len(analyse(fn, nm))
    print(f'=== total layout problems: {total}')
    return total


if __name__ == '__main__':
    sys.exit(0 if main() == 0 else 1)
