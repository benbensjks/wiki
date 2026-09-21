"""Strict numerical convergence audit for the two bit0 operating regimes.

The selected points include four add_growth=True passes, two adjacent failures,
and two add_growth=False passes.  ZENG and model.py are not modified.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import asdict
from pathlib import Path

import numpy as np

from model import Extension, Model, ROOT
from rescan800_readwindow import single_bit_verdict
from verify_bit0_part2 import CAND, clock_cycles
from verify_twobit_causal import _read_windows


POINTS = [
    dict(name='growth_pass_low_kon', uM_per_au=1.0, maturation_half_life_min=10.0,
         complex_on_au_inv_h=0.1, complex_off_h=3.0, add_growth=True),
    dict(name='growth_pass_fast_binding', uM_per_au=1.0, maturation_half_life_min=10.0,
         complex_on_au_inv_h=1.0, complex_off_h=1.0, add_growth=True),
    dict(name='growth_pass_high_kon', uM_per_au=1.0, maturation_half_life_min=10.0,
         complex_on_au_inv_h=3.0, complex_off_h=10.0, add_growth=True),
    dict(name='growth_pass_mid', uM_per_au=1.0, maturation_half_life_min=10.0,
         complex_on_au_inv_h=0.3, complex_off_h=3.0, add_growth=True),
    dict(name='growth_fail_koff_neighbor', uM_per_au=1.0, maturation_half_life_min=10.0,
         complex_on_au_inv_h=0.1, complex_off_h=1.0, add_growth=True),
    dict(name='growth_fail_maturation_neighbor', uM_per_au=1.0, maturation_half_life_min=5.0,
         complex_on_au_inv_h=0.1, complex_off_h=3.0, add_growth=True),
    dict(name='nogrowth_pass_reference', uM_per_au=10.0, maturation_half_life_min=10.0,
         complex_on_au_inv_h=0.1, complex_off_h=1.0, add_growth=False),
    dict(name='nogrowth_pass_slow_maturation', uM_per_au=10.0, maturation_half_life_min=20.0,
         complex_on_au_inv_h=0.3, complex_off_h=10.0, add_growth=False),
]

SETTINGS = [
    dict(name='archived_4min', rtol=2e-7, atol=2e-9, max_step_min=4.0),
    dict(name='default_2min', rtol=2e-7, atol=2e-9, max_step_min=2.0),
    dict(name='default_1min', rtol=2e-7, atol=2e-9, max_step_min=1.0),
    dict(name='strict_1min', rtol=1e-9, atol=1e-11, max_step_min=1.0),
]


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest().upper()


def run_setting(cfg, setting, hours=300.0):
    model = Model(Extension(**cfg))
    sol = model.simulate_bit0(hours=hours, sample_min=2.0,
                              rtol=setting['rtol'], atol=setting['atol'],
                              max_step_min=setting['max_step_min'])
    t, S = sol.t, sol.y[16]
    flux = np.asarray([model.flux(model._expand_bit0(sol.y[:, k]))
                       for k in range(sol.y.shape[1])])
    reads = _read_windows(t, flux, S, S)
    cold = single_bit_verdict(reads, 0)
    steady = single_bit_verdict(reads, 4)
    peaks = clock_cycles(t, flux)
    return dict(setting=setting['name'], rtol=setting['rtol'], atol=setting['atol'],
                max_step_min=setting['max_step_min'], solver_success=bool(sol.success),
                windows=len(reads), clock_period_h=float(np.median(np.diff(t[peaks]))),
                cold_passed=cold['passed'], cold_sequence=cold['sequence'],
                steady_passed=steady['passed'], steady_sequence=steady['sequence'],
                steady_min_commitment=steady['minimum_commitment'],
                steady_median_commitment=steady['median_commitment'],
                S_min=float(S.min()), S_max=float(S.max()),
                S_dwell_low=float(np.percentile(S, 2)),
                S_dwell_high=float(np.percentile(S, 98)))


def compare_to_strict(runs):
    ref = runs[-1]
    numeric_fields = ('clock_period_h', 'steady_min_commitment',
                      'steady_median_commitment', 'S_min', 'S_max',
                      'S_dwell_low', 'S_dwell_high')
    rows = []
    for run in runs[:-1]:
        diffs = {k: abs(float(run[k]) - float(ref[k])) for k in numeric_fields}
        rows.append(dict(setting=run['setting'],
                         same_cold_sequence=run['cold_sequence'] == ref['cold_sequence'],
                         same_steady_sequence=run['steady_sequence'] == ref['steady_sequence'],
                         same_verdict=(run['cold_passed'] == ref['cold_passed'] and
                                       run['steady_passed'] == ref['steady_passed']),
                         absolute_differences=diffs,
                         max_state_metric_difference=max(diffs[k] for k in
                                                         ('S_min', 'S_max', 'S_dwell_low',
                                                          'S_dwell_high'))))
    converged = bool(all(r['same_cold_sequence'] and r['same_steady_sequence'] and
                         r['same_verdict'] and r['max_state_metric_difference'] <= 0.005
                         for r in rows))
    return dict(reference=ref['setting'], comparisons=rows, converged=converged)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--points', default='all', help='all or comma-separated point indices')
    ap.add_argument('--hours', type=float, default=300.0)
    ap.add_argument('--output', type=Path,
                    default=ROOT / 'bit0_results' / 'verification' /
                    'growth_regime_convergence.json')
    args = ap.parse_args()
    ids = range(len(POINTS)) if args.points == 'all' else [int(x) for x in args.points.split(',')]

    base = asdict(Extension())
    base.update(CAND)
    results = []
    args.output.parent.mkdir(parents=True, exist_ok=True)
    for i in ids:
        spec = POINTS[i]
        cfg = base.copy()
        cfg.update({k: v for k, v in spec.items() if k != 'name'})
        runs = []
        for setting in SETTINGS:
            run = run_setting(cfg, setting, args.hours)
            runs.append(run)
            print(f"point={i} {spec['name']} setting={setting['name']} "
                  f"steady={run['steady_passed']} seq={run['steady_sequence']}", flush=True)
        result = dict(index=i, point=spec, runs=runs, convergence=compare_to_strict(runs))
        results.append(result)
        args.output.write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding='utf-8')
        print(f"point={i} converged={result['convergence']['converged']}", flush=True)

    meta = dict(hours=args.hours, points=results,
                sha256=dict(model_py=sha256(ROOT / 'model.py'),
                            verifier_py=sha256(ROOT / 'verify_twobit_causal.py'),
                            this_script=sha256(Path(__file__))))
    args.output.write_text(json.dumps(meta, ensure_ascii=False, indent=2), encoding='utf-8')
    print(f'wrote {args.output}')


if __name__ == '__main__':
    main()
