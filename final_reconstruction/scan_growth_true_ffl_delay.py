"""Scan separate A0/F0 maturation times in the growth-true regime."""
from __future__ import annotations

import argparse
import csv
import itertools
from dataclasses import asdict
from pathlib import Path

from model import Extension, ROOT
from verify_bit0_part2 import CAND
from verify_carry_delay import simulate


CANDIDATES = (
    dict(candidate='kinetic_centre', uM_per_au=0.8,
         maturation_half_life_min=14.0,
         complex_on_au_inv_h=1.0, complex_off_h=3.0),
    dict(candidate='highest_S_high', uM_per_au=0.8,
         maturation_half_life_min=14.0,
         complex_on_au_inv_h=1.0, complex_off_h=1.0),
)
A_MAT = (2.0, 5.0, 10.0, 20.0)
F_MAT = (10.0, 20.0, 30.0, 45.0, 60.0)


def jobs():
    return [(c, a, f) for c, a, f in itertools.product(CANDIDATES, A_MAT, F_MAT)
            if f > a]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--shard', type=int, default=0)
    ap.add_argument('--nshards', type=int, default=1)
    ap.add_argument('--hours', type=float, default=300.0)
    ap.add_argument('--out-dir', type=Path,
                    default=ROOT / 'bit0_results' / 'growth_true_ffl_delay')
    args = ap.parse_args()
    mine = jobs()[args.shard::args.nshards]
    args.out_dir.mkdir(parents=True, exist_ok=True)
    tag = f'shard{args.shard:02d}of{args.nshards:02d}'
    path = args.out_dir / f'ffl_delay_{tag}.csv'
    base = asdict(Extension()); base.update(CAND); base['add_growth'] = True
    rows = []
    for i, (candidate, amat, fmat) in enumerate(mine):
        cfg = base.copy(); cfg.update({k: v for k, v in candidate.items() if k != 'candidate'})
        result = simulate(cfg, carry_mrna_half_life_min=2.0,
                          carry_maturation_half_life_min=amat,
                          carry_repressor_maturation_half_life_min=fmat,
                          hours=args.hours, sample_min=2.0)
        c = result['causal_verdict']; sig = result['signal_ranges']
        row = dict(candidate=candidate['candidate'], A0_maturation_min=amat,
                   F0_maturation_min=fmat, uM_per_au=cfg['uM_per_au'],
                   bit_maturation_min=cfg['maturation_half_life_min'],
                   kon=cfg['complex_on_au_inv_h'], koff=cfg['complex_off_h'],
                   g0_peak=sig['g0'][1], Int1_source_peak=sig['Int1_source_peak'],
                   S1_min=sig['S1'][0], S1_max=sig['S1'][1],
                   cold_passed=result['cold_start']['passed'],
                   steady_passed=result['steady_state']['passed'],
                   steady_sequence=result['steady_state']['sequence'],
                   one_to_one=c['exactly_one_gate_and_flip_per_late_reverse'],
                   causal_order=c['causal_order_after_reverse_start'],
                   directions_alternate=c['bit1_directions_alternate'],
                   certified=result['certified'],
                   bit1_setup_h=result['margins']['bit1']['min_setup_h'],
                   bit1_hold_h=result['margins']['bit1']['min_hold_h'])
        rows.append(row)
        with path.open('w', newline='', encoding='utf-8') as fh:
            w = csv.DictWriter(fh, fieldnames=list(row)); w.writeheader(); w.writerows(rows)
        print(f"[{tag}] {i + 1}/{len(mine)} {candidate['candidate']} A/F={amat:g}/{fmat:g} "
              f"g0={row['g0_peak']:.3g} steady={row['steady_passed']} "
              f"order={row['causal_order']} certified={row['certified']}", flush=True)
    print(f'[{tag}] wrote {path}')


if __name__ == '__main__':
    main()
