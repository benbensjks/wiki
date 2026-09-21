"""Post-scan spot check: recompute selected points in one fresh process and
compare the full label string, commitments, state extremes and archived-rule
codes against the merged 800-point CSV.

Selection is deterministic (fixed seed) and always includes
  * boundary points on the extension axes (both add_growth values),
  * points where the new criterion and the archived rule disagree (both ways),
  * five random interior points.
"""
from __future__ import annotations

import argparse
import json
import random
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from model import ROOT
from rescan800_readwindow import evaluate, sha256

SCAN_DIR = ROOT / 'bit0_results' / 'rescan800'
KEYS = ['uM_per_au', 'maturation_half_life_min', 'complex_on_au_inv_h', 'complex_off_h', 'add_growth']
NUMERIC = KEYS[:-1]
FLOAT_TOL = 1e-9
COMPARE = ['cold_start_sequence', 'steady_state_sequence', 'archived_rule_codes',
           'cold_start_min_commitment', 'steady_state_min_commitment',
           'steady_state_median_commitment', 'S0_min', 'S0_max',
           'S0_dwell_low', 'S0_dwell_high', 'clock_period_h', 'windows']
BOOLS = ['cold_start_passed', 'steady_state_passed', 'archived_rule_stable']


def pick_points(merged, n_random=5, seed=20260919):
    picks = []

    def add(row, why):
        key = (round(row.uM_per_au, 10), round(row.maturation_half_life_min, 10),
               round(row.complex_on_au_inv_h, 10), round(row.complex_off_h, 10), bool(row.add_growth))
        if key not in {p['key'] for p in picks}:
            picks.append(dict(key=key, why=why, row=row))

    for axis in NUMERIC:
        for growth in (False, True):
            sub = merged[merged.add_growth == growth]
            if not len(sub):
                continue
            add(sub.loc[sub[axis].idxmin()], f'boundary {axis}=min growth={growth}')
            add(sub.loc[sub[axis].idxmax()], f'boundary {axis}=max growth={growth}')
    new = merged.steady_state_passed.to_numpy(dtype=bool)
    old = merged.archived100h_stable.fillna(False).to_numpy(dtype=bool)
    for mask, why in (((new & ~old), 'new pass / archived 100 h fail'),
                      ((~new & old), 'archived 100 h pass / new fail')):
        idx = np.flatnonzero(mask)
        for i in idx[:2]:
            add(merged.iloc[i], why)
    rng = random.Random(seed)
    for i in rng.sample(range(len(merged)), min(n_random, len(merged))):
        add(merged.iloc[i], 'random interior')
    return picks


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--hours', type=float, default=300.0)
    ap.add_argument('--n-random', type=int, default=5)
    ap.add_argument('--seed', type=int, default=20260919)
    args = ap.parse_args()

    merged = pd.read_csv(SCAN_DIR / 'rescan800_all.csv')
    picks = pick_points(merged, args.n_random, args.seed)
    print(f'checking {len(picks)} points '
          f'({", ".join(sorted({p["why"].split()[0] for p in picks}))})', flush=True)

    results = []
    for p in picks:
        row = p['row']
        cfg = {k: (bool(row[k]) if k == 'add_growth' else float(row[k])) for k in KEYS}
        rec = evaluate(cfg, args.hours)
        diffs = {}
        for col in COMPARE:
            a, b = row.get(col), rec.get(col)
            if isinstance(b, str) or isinstance(a, str):
                same = str(a) == str(b)
            else:
                same = bool(np.isfinite(a) and np.isfinite(b) and abs(float(a) - float(b)) <= FLOAT_TOL)
            if not same:
                diffs[col] = dict(csv=a, recomputed=b)
        for col in BOOLS:
            if bool(row.get(col)) != bool(rec.get(col)):
                diffs[col] = dict(csv=bool(row.get(col)), recomputed=bool(rec.get(col)))
        results.append(dict(point=cfg, why=p['why'], matches=not diffs, differences=diffs,
                            csv_cold=bool(row.cold_start_passed),
                            csv_steady=bool(row.steady_state_passed),
                            recomputed_cold=bool(rec['cold_start_passed']),
                            recomputed_steady=bool(rec['steady_state_passed']),
                            sequence=rec.get('cold_start_sequence', '')))
        status = 'OK ' if not diffs else 'DIFF'
        print(f"[{status}] {p['why']:<34} uM={cfg['uM_per_au']:<5g} tmat={cfg['maturation_half_life_min']:<5g} "
              f"kon={cfg['complex_on_au_inv_h']:<4g} koff={cfg['complex_off_h']:<5g} "
              f"growth={str(cfg['add_growth']):<5} seq={rec.get('cold_start_sequence', '')}"
              + (f' diffs={list(diffs)}' if diffs else ''), flush=True)

    report = dict(hours=args.hours, seed=args.seed, checked=len(results),
                  matched=sum(1 for r in results if r['matches']),
                  mismatched=[r for r in results if not r['matches']],
                  environment=dict(python=sys.version,
                                   numpy=np.__version__,
                                   pandas=pd.__version__),
                  sha256=dict(model_py=sha256(ROOT / 'model.py'),
                              verifier_py=sha256(ROOT / 'verify_twobit_causal.py'),
                              scanner_py=sha256(ROOT / 'rescan800_readwindow.py'),
                              merged_csv=sha256(SCAN_DIR / 'rescan800_all.csv'),
                              this_script=sha256(Path(__file__))),
                  points=results)
    out = SCAN_DIR / 'rescan800_spotcheck.json'
    out.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
    print(f'\nmatched {report["matched"]}/{report["checked"]} -> wrote {out}')


if __name__ == '__main__':
    main()
