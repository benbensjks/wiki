"""Local 150-point map around the add_growth=True bit0 regime.

Stage A scans uM_per_au x maturation time for each of the six kon/koff pairs
that passed in the archived grid.  Each point uses the same 300 h finite read
window criterion as rescan800_readwindow.py.  ZENG is unchanged.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import sys
from dataclasses import asdict
from pathlib import Path

from model import Extension, ROOT
from rescan800_readwindow import evaluate
from verify_bit0_part2 import CAND


UM_VALUES = (0.6, 0.8, 1.0, 1.2, 1.4)
MAT_VALUES = (6.0, 8.0, 10.0, 12.0, 14.0)
KINETIC_PAIRS = ((0.1, 3.0), (0.3, 3.0), (1.0, 1.0),
                 (1.0, 3.0), (3.0, 3.0), (3.0, 10.0))


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest().upper()


def configs():
    base = asdict(Extension())
    base.update(CAND)
    rows = []
    for u, mat, (kon, koff) in itertools.product(UM_VALUES, MAT_VALUES, KINETIC_PAIRS):
        cfg = base.copy()
        cfg.update(uM_per_au=u, maturation_half_life_min=mat,
                   complex_on_au_inv_h=kon, complex_off_h=koff,
                   add_growth=True)
        rows.append(cfg)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--shard', type=int, default=0)
    ap.add_argument('--nshards', type=int, default=1)
    ap.add_argument('--hours', type=float, default=300.0)
    ap.add_argument('--out-dir', type=Path,
                    default=ROOT / 'bit0_results' / 'growth_true_local')
    args = ap.parse_args()

    all_cfg = configs()
    mine = all_cfg[args.shard::args.nshards]
    args.out_dir.mkdir(parents=True, exist_ok=True)
    tag = f'shard{args.shard:02d}of{args.nshards:02d}'
    out = args.out_dir / f'growth_local_{tag}.csv'
    rows = []
    for i, cfg in enumerate(mine):
        row = evaluate(cfg, args.hours)
        rows.append(row)
        fields = sorted({k for r in rows for k in r})
        with out.open('w', newline='', encoding='utf-8') as fh:
            writer = csv.DictWriter(fh, fieldnames=fields, extrasaction='ignore')
            writer.writeheader()
            writer.writerows(rows)
        print(f"[{tag}] {i + 1}/{len(mine)} uM={cfg['uM_per_au']:.1f} "
              f"mat={cfg['maturation_half_life_min']:.0f} "
              f"kon/koff={cfg['complex_on_au_inv_h']:g}/{cfg['complex_off_h']:g} "
              f"steady={row['steady_state_passed']} seq={row['steady_state_sequence']}", flush=True)

    meta = dict(shard=args.shard, nshards=args.nshards, hours=args.hours,
                points=len(mine), total_points=len(all_cfg),
                grid=dict(uM_per_au=UM_VALUES,
                          maturation_half_life_min=MAT_VALUES,
                          kinetic_pairs=KINETIC_PAIRS, add_growth=True),
                python=sys.version,
                sha256=dict(model_py=sha256(ROOT / 'model.py'),
                            verifier_py=sha256(ROOT / 'verify_twobit_causal.py'),
                            evaluator_py=sha256(ROOT / 'rescan800_readwindow.py'),
                            scanner_py=sha256(Path(__file__)),
                            csv=sha256(out)))
    (args.out_dir / f'meta_{tag}.json').write_text(
        json.dumps(meta, ensure_ascii=False, indent=2), encoding='utf-8')
    print(f'[{tag}] wrote {out}')


if __name__ == '__main__':
    main()
