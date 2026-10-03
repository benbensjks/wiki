"""Report the M-18 figure hashes and check the heatmap panel content."""
from __future__ import annotations

import csv
import hashlib
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / 'out'
RES = HERE / 'results' / 'clock_axis_20261001_180035_369646'


def main() -> int:
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8')
    rows = list(csv.DictReader((RES / 'all.csv').open(encoding='utf-8')))
    fails = [r for r in rows if r['counting_passed'] != 'True'
             and r['kind'] == 'hzz']
    n_x = sum(r['steady_sequence'].count('x') for r in fails)
    n_lab = sum(len(r['steady_sequence']) - r['steady_sequence'].count('x')
                for r in fails)
    print('failing HZZ points: %d   unlabelled windows: %d   labelled: %d'
          % (len(fails), n_x, n_lab))

    s = (OUT / 'fig_M18_map_2cell_failmode.svg').read_text(encoding='utf-8')
    xs = re.findall(r'<text[^>]*font-size: 7.5px[^>]*>x</text>', s)
    print('"x" markers drawn in the heatmap: %d  ->  %s'
          % (len(xs), 'OK' if len(xs) == n_x else 'MISMATCH'))
    ylab = re.findall(r'<text[^>]*font-size: 9px[^>]*>([^<]*(?:OFF|ON )[^<]*)</text>', s)
    print('heatmap row labels: %d  ->  %s'
          % (len(ylab), 'OK' if len(ylab) == len(fails) else 'MISMATCH'))
    for t in ylab[:4] + ylab[-2:]:
        print('   ', t)

    meta = json.loads((OUT / 'fig_M18_map_provenance.json').read_text(encoding='utf-8'))
    print('\n== outputs ==')
    for k, v in sorted(meta['outputs'].items()):
        p = OUT / k
        print('  %-34s %9d B  %s' % (k, p.stat().st_size, v[:16]))
    print('\n== audit summary recorded in provenance ==')
    for k in sorted(meta):
        if k.startswith('audit_'):
            a = meta[k]
            print('  %-30s texts=%3d TT=%d TC=%d OUT=%d'
                  % (k[6:], a['texts'], len(a['text_text']), len(a['text_cell']),
                     len(a['outside'])))
    print('\nboundary intervals:')
    for arm, b in meta['boundaries'].items():
        print('  %-9s (%.3f, %.3f]   r in (%.4f, %.4f]'
              % (arm, b['interval'][0], b['interval'][1], b['r_last'], b['r_first']))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
