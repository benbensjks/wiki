"""Layout / containment check for the M-18 map figures.

Verifies, numerically (I cannot view images):
  * every <text> element lies inside the SVG canvas  (OUT must be 0)
  * the drawn map reproduces all.csv exactly          (no silent data drift)
  * the pass/fail counts and the boundary intervals agree with the scan summary
"""
from __future__ import annotations

import csv
import hashlib
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RES = HERE / 'results' / 'clock_axis_20261001_180035_369646'
OUT = HERE / 'out'
EXPECTED = '1234567012345670123'


def main() -> int:
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8')
    fails = []

    # ---- 1. SVG text containment -------------------------------------------
    print('== SVG text containment ==')
    for p in sorted(OUT.glob('fig_M18_map_*.svg')):
        s = p.read_text(encoding='utf-8')
        mt = re.search(r'viewBox="([-0-9. ]+)"', s)
        if not mt:
            fails.append('%s: no viewBox' % p.name)
            continue
        x0, y0, w, h = [float(v) for v in mt.group(1).split()]
        tot = bad = 0
        worst = 0.0
        pat = re.compile(r'<text[^>]*?x="([-0-9.eE]+)"[^>]*?y="([-0-9.eE]+)"')
        for m in pat.finditer(s):
            tot += 1
            x, y = float(m.group(1)), float(m.group(2))
            dx = max(x0 - x, x - (x0 + w), 0.0)
            dy = max(y0 - y, y - (y0 + h), 0.0)
            if dx > 0 or dy > 0:
                bad += 1
                worst = max(worst, max(dx, dy))
        print('  %-32s texts=%3d  OUT=%d  worst=%.2f  canvas=%.0fx%.0f'
              % (p.name, tot, bad, worst, w, h))
        if bad:
            fails.append('%s: %d text elements outside the canvas' % (p.name, bad))
        if tot < 20:
            fails.append('%s: only %d <text> elements (svg.fonttype unset?)'
                         % (p.name, tot))

    # ---- 2. provenance vs source -------------------------------------------
    print('\n== provenance vs all.csv ==')
    meta = json.loads((OUT / 'fig_M18_map_provenance.json').read_text(encoding='utf-8'))
    src = RES / 'all.csv'
    if meta['source_sha256'] != hashlib.sha256(src.read_bytes()).hexdigest().upper():
        fails.append('provenance sha does not match all.csv')
    rows = list(csv.DictReader(src.open(encoding='utf-8')))
    Ks = sorted({float(r['clock_K']) for r in rows})
    if [round(k, 6) for k in meta['K_levels']] != [round(k, 6) for k in Ks]:
        fails.append('K levels differ from all.csv')
    n_ok = sum(1 for r in rows if r['counting_passed'] == 'True')
    n_ev = sum(1 for r in rows if r['event_causality_passed'] == 'True')
    print('  rows=%d  K=%d  counting_pass=%d  events_pass=%d'
          % (len(rows), len(Ks), n_ok, n_ev))

    # certified == counting and events, re-checked here
    viol = [r['label'] for r in rows
            if (r['certified_v1'] == 'True')
            != (r['counting_passed'] == 'True' and r['event_causality_passed'] == 'True')]
    print('  certified != counting and events : %d violations' % len(viol))
    if viol:
        fails.append('certified_v1 is not counting and events for %d rows' % len(viol))

    # every all-pass sequence equals the expected staircase
    seqs = {r['steady_sequence'] for r in rows if r['counting_passed'] == 'True'}
    print('  distinct passing sequences: %s' % sorted(seqs))
    if seqs != {EXPECTED}:
        fails.append('passing sequences are not all the expected staircase')

    # ---- 3. boundaries ------------------------------------------------------
    print('\n== boundaries ==')
    for arm, b in meta['boundaries'].items():
        kind, a = arm.split('/')
        per_col = []
        for K in Ks:
            vals = [r['counting_passed'] == 'True' for r in rows
                    if r['kind'] == kind and r['arm'] == a
                    and abs(float(r['clock_K']) - K) < 1e-9]
            per_col.append(all(vals) and len(vals) == 2)
        last_ok = max(j for j, v in enumerate(per_col) if v)
        first_bad = min(j for j, v in enumerate(per_col) if not v)
        ok = (abs(b['last_passing_K'] - Ks[last_ok]) < 1e-9
              and abs(b['first_failing_K'] - Ks[first_bad]) < 1e-9)
        print('  %-9s last pass %.3f  first fail %.3f  r=(%.4f, %.4f]  %s'
              % (arm, b['last_passing_K'], b['first_failing_K'],
                 b['r_last'], b['r_first'], 'OK' if ok else 'MISMATCH'))
        if not ok:
            fails.append('%s: boundary interval disagrees with all.csv' % arm)

    # HZH must have no boundary, and must pass at both ends
    hzh = [r for r in rows if r['kind'] == 'hzh']
    if 'hzh/base' in meta['boundaries']:
        fails.append('HZH was given a boundary although it passes everywhere')
    print('  hzh/base  passes %d/%d, K range %.3f-%.3f, no boundary (both sides '
          'unscanned)' % (sum(1 for r in hzh if r['counting_passed'] == 'True'),
                          len(hzh), min(float(r['clock_K']) for r in hzh),
                          max(float(r['clock_K']) for r in hzh)))

    # ---- 4. failure-mode classifier ----------------------------------------
    print('\n== failure modes ==')
    sys.path.insert(0, str(HERE))
    from make_fig_M18_clock_map import failure_mode
    counts = {}
    bad = []
    for r in rows:
        if r['counting_passed'] == 'True':
            if failure_mode(r['steady_sequence']) != '':
                bad.append(r['label'])
            continue
        m = failure_mode(r['steady_sequence'])
        counts[m] = counts.get(m, 0) + 1
        if m not in ('s', 'r', 'u'):
            bad.append(r['label'])
    print('  failing points: %d  modes %s' % (20, counts))
    if bad:
        fails.append('unclassifiable rows: %s' % bad)

    # ---- 5. glyph content, re-derived and read back out of the SVG ---------
    # Independent of the figure script: rebuild the expected glyph for every cell
    # straight from all.csv, then count the white bold glyphs actually emitted in
    # the SVG, grouped by column.  Catches a column/row mapping or transposition
    # error that a "figure ran fine" check cannot see.
    print('\n== glyph content, SVG vs all.csv ==')
    COL = {'counting': 'counting_passed', 'events': 'event_causality_passed',
           'certified': 'certified_v1'}
    gly = re.compile(r'<text[^>]*font-size: 12px[^>]*fill: #ffffff[^>]*'
                     r'x="([-0-9.]+)"[^>]*>([^<]+)</text>')
    for ncells in (2, 3):
        name = ('fig_M18_map_2cell' if ncells == 2 else 'fig_M18_map_3cell')
        p = OUT / (name + '.svg')
        if not p.exists():
            fails.append('%s missing' % p.name)
            continue
        fields = ['counting', 'events', 'certified'][:ncells]
        exp = [{} for _ in Ks]
        for r in rows:
            j = Ks.index(float(r['clock_K']))
            for f in fields:
                if r[COL[f]] == 'True':
                    g = '√'
                elif f == 'counting':
                    g = failure_mode(r['steady_sequence']) or '×'
                else:
                    g = '×'
                exp[j][g] = exp[j].get(g, 0) + 1
        got = {}
        for m in gly.finditer(p.read_text(encoding='utf-8')):
            got.setdefault(round(float(m.group(1)), 1), {})
            d = got[round(float(m.group(1)), 1)]
            d[m.group(2)] = d.get(m.group(2), 0) + 1
        ok = len(got) == len(Ks)
        cols = [got[k] for k in sorted(got)]
        for j in range(min(len(cols), len(exp))):
            if cols[j] != exp[j]:
                ok = False
                print('    col %d (K=%.3f) SVG=%s  expected=%s'
                      % (j, Ks[j], cols[j], exp[j]))
        tot_g = sum(sum(c.values()) for c in cols)
        print('  %-26s columns_in_svg=%d/%d  glyphs=%d  expect %d cells  %s'
              % (name, len(got), len(Ks), tot_g, 6 * ncells * len(Ks),
                 'OK' if ok else 'MISMATCH'))
        if not ok:
            fails.append('%s: glyph map does not match all.csv' % name)

    print('\n== result ==')
    if fails:
        for f in fails:
            print('  FAIL', f)
        return 1
    print('  ALL CHECKS PASS')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
