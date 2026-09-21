"""Upstream-phase robustness check for the selected carry-delay candidate."""
from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path

import numpy as np

from model import Extension, Model
from verify_bit0_part2 import CAND, OUT
from verify_carry_delay import simulate


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--phases', default='0,2,4,6,8,10')
    ap.add_argument('--uM', type=float, default=6.0)
    ap.add_argument('--maturation-min', type=float, default=30.0)
    ap.add_argument('--hours', type=float, default=300.0)
    ap.add_argument('--output', type=Path, default=OUT / 'carry_delay_phase_scan.json')
    args = ap.parse_args()

    phases = [float(x) for x in args.phases.split(',')]
    base = asdict(Extension())
    base.update(CAND)
    base['uM_per_au'] = args.uM

    # The upstream is autonomous and downstream-independent.  Take six
    # different points from one settled upstream cycle, while resetting every
    # downstream/carry state to the identical PB initial condition.
    model = Model(Extension(**base))
    warm_end = 120.0 + max(phases)
    warm = model.simulate_twobit(hours=warm_end, sample_min=2.0, max_step_min=2.0)

    rows = []
    args.output.parent.mkdir(parents=True, exist_ok=True)
    for phase in phases:
        k = int(np.argmin(np.abs(warm.t - (120.0 + phase))))
        row = simulate(base, 2.0, args.maturation_min, hours=args.hours,
                       upstream_initial7=warm.y[:7, k])
        row['upstream_phase_offset_h'] = phase
        rows.append(row)
        args.output.write_text(json.dumps(rows, ensure_ascii=False, indent=2), encoding='utf-8')
        print(f"phase={phase:g}h cold={row['cold_start']['passed']} "
              f"steady={row['steady_state']['passed']} order={row['causal_verdict']['causal_order_after_reverse_start']} "
              f"certified={row['certified']} seq={row['steady_state']['sequence']}", flush=True)
    print(f'wrote {args.output}')


if __name__ == '__main__':
    main()
