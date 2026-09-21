"""Serial shard re-scan of the original 800-point Extension grid.

Read-window criterion (single-bit reduction of verify_twobit_causal):
  * window  = flux trough +/- READ_WINDOW_FRACTION/2 of the clock cycle (-_-)
  * label   = '1' if >=80% of the window is >=0.70, '0' if >=80% is <=0.30, else 'x'
  * cold    = every window labelled and labels alternate (mod 2) from read 0
  * steady  = same after dropping STEADY_DROP_READS leading reads
The window geometry, bands and occupancy rule are imported verbatim from
verify_twobit_causal; only the "two bits -> one bit" reduction is new here.

The archived rule (single trough sample, 0.8/0.2 bands) is also evaluated on the
SAME trajectory so the two criteria can be compared point by point.

Model, ZENG and the grid definition are never touched.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import platform
import sys
import time
import traceback
from dataclasses import asdict
from pathlib import Path

import numpy as np

from model import Extension, Model, ROOT
from bit0_diagnostic import GRID, cycle_metrics, score_cycles
from verify_twobit_causal import (MIN_OCCUPANCY, READ_WINDOW_FRACTION,
                                  STEADY_DROP_READS, _read_windows)
from verify_bit0_part2 import CAND, clock_cycles

VERIFIER = ROOT / 'verify_twobit_causal.py'
SCANNER = Path(__file__).resolve()
GRID_KEYS = list(GRID)
HOURS_DEFAULT = 300.0
SAMPLE_MIN = 2.0
MAX_STEP_MIN = 4.0            # identical to the archived 800-point scan numerics


def grid_configs():
    base = asdict(Extension())
    base.update(CAND)
    out = []
    for vals in itertools.product(*(GRID[k] for k in GRID_KEYS)):
        cfg = base.copy()
        cfg.update(dict(zip(GRID_KEYS, vals)))
        out.append(cfg)
    return out


def single_bit_verdict(reads, drop):
    """Faithful mod-2 reduction of verify_twobit_causal._counter_verdict."""
    chosen = reads[drop:]
    labels = [r['bit0']['label'] for r in chosen]
    valid = bool(len(labels) >= 8 and all(l != 'x' for l in labels))
    increments = bool(valid and all((int(labels[i + 1]) - int(labels[i])) % 2 == 1
                                    for i in range(len(labels) - 1)))
    return dict(drop_reads=drop, reads=len(labels), valid_windows=valid,
                increments_mod2=increments, passed=bool(valid and increments),
                sequence=''.join(labels),
                minimum_commitment=(float(min(r['bit0']['commitment'] for r in chosen))
                                    if chosen else 0.0),
                median_commitment=(float(np.median([r['bit0']['commitment'] for r in chosen]))
                                   if chosen else 0.0))


def evaluate(cfg, hours):
    t0 = time.perf_counter()
    m = Model(Extension(**cfg))
    row = {k: cfg[k] for k in ('uM_per_au', 'maturation_half_life_min', 'complex_on_au_inv_h',
                               'complex_off_h', 'add_growth')}
    try:
        sol = m.simulate_bit0(hours=hours, sample_min=SAMPLE_MIN, max_step_min=MAX_STEP_MIN)
        t = sol.t
        S0 = sol.y[16]
        flux = np.array([m.flux(m._expand_bit0(sol.y[:, k])) for k in range(sol.y.shape[1])])
        reads = _read_windows(t, flux, S0, S0)      # bit1 slot reuses S0 and is ignored
        cold = single_bit_verdict(reads, 0)
        steady = single_bit_verdict(reads, STEADY_DROP_READS)
        cycles = cycle_metrics(m, sol)[2]
        arch = score_cycles(cycles)
        row.update(dict(solver=True, error='',
                        windows=len(reads),
                        clock_period_h=float(np.median(np.diff(t[clock_cycles(t, flux)]))),
                        cold_start_passed=cold['passed'],
                        cold_start_sequence=cold['sequence'],
                        cold_start_min_commitment=cold['minimum_commitment'],
                        steady_state_passed=steady['passed'],
                        steady_state_sequence=steady['sequence'],
                        steady_state_min_commitment=steady['minimum_commitment'],
                        steady_state_median_commitment=steady['median_commitment'],
                        steady_reads=steady['reads'],
                        unlabelled_windows=int(sum(1 for r in reads if 'x' in r['bit0']['label'])),
                        S0_min=float(S0.min()), S0_max=float(S0.max()),
                        S0_dwell_low=float(np.percentile(S0, 2)),
                        S0_dwell_high=float(np.percentile(S0, 98)),
                        archived_rule_codes=arch['codes'],
                        archived_rule_stable=bool(arch['stable_alternation']),
                        archived_rule_confident_fraction=float(arch['confident_fraction']),
                        archived_rule_dynamic_range=float(arch['dynamic_range'])))
    except Exception as exc:                                  # noqa: BLE001
        row.update(dict(solver=False, error=repr(exc), cold_start_passed=False,
                        steady_state_passed=False, cold_start_sequence='',
                        steady_state_sequence='', archived_rule_stable=False,
                        error_trace=traceback.format_exc(limit=2)))
    row['runtime_s'] = round(time.perf_counter() - t0, 2)
    return row


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest().upper()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--shard', type=int, default=0)
    ap.add_argument('--nshards', type=int, default=1)
    ap.add_argument('--hours', type=float, default=HOURS_DEFAULT)
    ap.add_argument('--out-dir', type=Path, default=ROOT / 'bit0_results' / 'rescan800')
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    configs = grid_configs()
    mine = configs[args.shard::args.nshards]
    tag = f'shard{args.shard:02d}of{args.nshards:02d}'
    csv_path = args.out_dir / f'rescan800_{tag}.csv'

    rows = []
    for i, cfg in enumerate(mine):
        rows.append(evaluate(cfg, args.hours))
        print(f'[{tag}] {i + 1}/{len(mine)} uM={cfg["uM_per_au"]:<5g} '
              f'tmat={cfg["maturation_half_life_min"]:<5g} kon={cfg["complex_on_au_inv_h"]:<4g} '
              f'koff={cfg["complex_off_h"]:<5g} growth={str(cfg["add_growth"]):<5} '
              f'cold={rows[-1]["cold_start_passed"]} steady={rows[-1]["steady_state_passed"]} '
              f'seq={rows[-1].get("steady_state_sequence", "")} ({rows[-1]["runtime_s"]}s)', flush=True)
        all_fields = sorted({k for r in rows for k in r})
        with csv_path.open('w', newline='', encoding='utf-8') as fh:
            w = csv.DictWriter(fh, fieldnames=all_fields, extrasaction='ignore')
            w.writeheader()
            w.writerows(rows)

    meta = dict(shard=args.shard, nshards=args.nshards, hours=args.hours,
                sample_min=SAMPLE_MIN, max_step_min=MAX_STEP_MIN,
                points=len(mine), grid_keys=GRID_KEYS,
                grid=GRID,
                criterion=dict(source='verify_twobit_causal',
                               band_low=0.30, band_high=0.70,
                               minimum_occupancy=MIN_OCCUPANCY,
                               read_window_fraction=READ_WINDOW_FRACTION,
                               steady_drop_reads=STEADY_DROP_READS,
                               reduction='mod 2 instead of mod 4'),
                python=sys.version, platform=platform.platform(),
                sha256=dict(model_py=sha256(ROOT / 'model.py'),
                            bit0_diagnostic_py=sha256(ROOT / 'bit0_diagnostic.py'),
                            verify_twobit_causal_py=sha256(VERIFIER),
                            scanner_py=sha256(SCANNER),
                            upstream_py=sha256(ROOT.parent / '其他小组成员任务/Week4_振荡器-C31建模/code'
                                               / 'Flux_Driven_Translation_Burden_Model.py')),
                csv=csv_path.name, csv_sha256=sha256(csv_path))
    (args.out_dir / f'meta_{tag}.json').write_text(
        json.dumps(meta, ensure_ascii=False, indent=2), encoding='utf-8')
    print(f'[{tag}] wrote {csv_path} ({len(rows)} rows)', flush=True)


if __name__ == '__main__':
    main()
