"""Full 4x4 kon/koff audit at two centres of the growth-true ridge."""
from __future__ import annotations

import argparse
import csv
import itertools
from dataclasses import asdict
from pathlib import Path

from model import Extension, ROOT
from rescan800_readwindow import evaluate
from verify_bit0_part2 import CAND


CENTRES = ((0.8, 14.0), (1.0, 10.0))
KON = (0.1, 0.3, 1.0, 3.0)
KOFF = (1.0, 3.0, 10.0, 30.0)


def configs():
    base = asdict(Extension())
    base.update(CAND)
    out = []
    for (u, mat), kon, koff in itertools.product(CENTRES, KON, KOFF):
        cfg = base.copy()
        cfg.update(uM_per_au=u, maturation_half_life_min=mat,
                   complex_on_au_inv_h=kon, complex_off_h=koff,
                   add_growth=True)
        out.append(cfg)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--shard', type=int, default=0)
    ap.add_argument('--nshards', type=int, default=1)
    ap.add_argument('--hours', type=float, default=300.0)
    ap.add_argument('--out-dir', type=Path,
                    default=ROOT / 'bit0_results' / 'growth_true_kinetics')
    args = ap.parse_args()
    mine = configs()[args.shard::args.nshards]
    args.out_dir.mkdir(parents=True, exist_ok=True)
    tag = f'shard{args.shard:02d}of{args.nshards:02d}'
    path = args.out_dir / f'growth_kinetics_{tag}.csv'
    rows = []
    for i, cfg in enumerate(mine):
        row = evaluate(cfg, args.hours)
        rows.append(row)
        fields = sorted({k for r in rows for k in r})
        with path.open('w', newline='', encoding='utf-8') as fh:
            w = csv.DictWriter(fh, fieldnames=fields, extrasaction='ignore')
            w.writeheader(); w.writerows(rows)
        print(f"[{tag}] {i + 1}/{len(mine)} uM/mat={cfg['uM_per_au']}/{cfg['maturation_half_life_min']} "
              f"kon/koff={cfg['complex_on_au_inv_h']}/{cfg['complex_off_h']} "
              f"pass={row['steady_state_passed']}", flush=True)
    print(f'[{tag}] wrote {path}')


if __name__ == '__main__':
    main()
