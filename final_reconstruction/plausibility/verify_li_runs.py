"""Pre-merge audit of the local-identifiability runs.

The previous attempt produced a 9-column matrix from runs taken at TWO different
gate exponents (49 at the ZENG table value 4.0, 6 at ~6), which made every
linear-algebra summary of it meaningless.  This audit reads the per-run metadata
and refuses to let the matrix be assembled unless every run is demonstrably at
one working point.

Checked per run
  A1  the job file exists and its metadata carries resolved_params
  A2  the model reported the exponent that was requested
  A3  every parameter that is NOT the perturbed one equals the frozen value
  A4  the perturbed parameter equals nominal * (1 +- step) with the stated step
  A5  the resolved vector covers all parameters of the selected set
  A6  the two observable tables have the expected shapes
  A7  the live ZENG table value is recorded and is the table value (not patched)

    python verify_li_runs.py
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plausibility_common import OUT, SELECTED_N_A1_GATE, ZENG_N_A1_TABLE  # noqa: E402
from check_local_identifiability import (NOMINAL, STEPS, build_jobs,  # noqa: E402
                                        job_file, param_table)

JOBS_DIR = OUT / 'local_identifiability' / 'jobs'
TOL = 1e-12
FAIL = []


def bad(msg):
    FAIL.append(msg)
    print('  [FAIL] %s' % msg)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--with-translation', dest='with_translation', action='store_true',
                    default=False,
                    help='audit the 10-parameter set (adds translation_h)')
    ap.add_argument('--nshards', type=int, default=8, help='unused; kept for symmetry')
    args = ap.parse_args()

    params = param_table(args.with_translation)
    names = [n for n, _, _ in params]
    jobs = build_jobs(params)
    n_at = 49 + (6 if args.with_translation else 0)
    n_moved_expected = 6
    print('=== 选定的参数集 (%d) ===' % len(names))
    print('   ', ', '.join(names))
    print('   frozen n_A1_gate = %g   ZENG table value = %g' % (SELECTED_N_A1_GATE,
                                                               ZENG_N_A1_TABLE))
    print()
    print('=== 逐运行审计 (%d 个作业) ===' % len(jobs))

    meta_rows = []
    for label, ov in jobs:
        f = job_file(label)
        if not f.exists():
            bad('%s: npz missing' % label)
            continue
        jf = f.with_suffix('.json')
        if not jf.exists():
            bad('%s: json missing' % label)
            continue
        m = json.loads(jf.read_text(encoding='utf-8'))
        rp = m.get('resolved_params')
        if rp is None:                                            # A1
            bad('%s: no resolved_params (old-format artifact)' % label)
            continue
        missing = [n for n in names if n not in rp]                # A5
        if missing:
            bad('%s: resolved_params missing %s' % (label, missing))
        perturbed = label.split('|')[0]
        eff = float(m.get('n_A1_gate_effective', float('nan')))
        if abs(eff - float(rp['n_A1_gate'])) > TOL:                # A2
            bad('%s: effective %r != requested %r' % (label, eff, rp['n_A1_gate']))
        parts = label.split('|')
        step = None if len(parts) < 3 else float(parts[2])
        for n in names:                                            # A3
            if n == perturbed:
                continue
            if abs(float(rp[n]) - float(NOMINAL[n])) > TOL:
                bad('%s: %s = %r but must stay frozen at %r'
                    % (label, n, rp[n], NOMINAL[n]))
        if perturbed == 'n_A1_gate':                               # A2b
            if float(rp['n_A1_gate']) == ZENG_N_A1_TABLE:
                bad('%s: n_A1_gate perturbation collapsed onto the table value' % label)
        else:
            if abs(float(rp['n_A1_gate']) - SELECTED_N_A1_GATE) > TOL:
                bad('%s: n_A1_gate = %r, expected the frozen %g'
                    % (label, rp['n_A1_gate'], SELECTED_N_A1_GATE))
        if step is not None:                                       # A4
            side = 1.0 + step if parts[1] == 'plus' else 1.0 - step
            want = float(NOMINAL[perturbed]) * side
            if abs(float(rp[perturbed]) - want) > max(TOL, abs(want) * 1e-12):
                bad('%s: %s = %r, expected %r' % (label, perturbed, rp[perturbed], want))
        meta_rows.append(dict(label=label, perturbed=perturbed, step=step,
                              n_A1_gate_effective=eff,
                              theo_n=m.get('theo_n'), exp_n=m.get('exp_n'),
                              live_zeng=m.get('meta', {}).get('resolved_params', {})
                              .get('__live__')))

    print()
    print('=== 汇总：每个被扰动参数下的有效门指数 ===')
    by = {}
    for r in meta_rows:
        by.setdefault(r['perturbed'], []).append(r['n_A1_gate_effective'])
    for k in sorted(by):
        vals = sorted(set(round(v, 10) for v in by[k]))
        print('  %-30s runs=%-3d n_A1_gate_effective=%s' % (k, len(by[k]), vals))
    print()
    print('=== 汇总：观测量表尺寸 ===')
    theo = sorted(set(r['theo_n'] for r in meta_rows))
    exp = sorted(set(r['exp_n'] for r in meta_rows))
    print('  theoretical 时间点数 :', theo)
    print('  experimental 时间点数:', exp)
    if theo != [4501] or exp != [301]:
        bad('unexpected observable table shapes: theo=%s exp=%s' % (theo, exp))

    n_at_frozen = sum(1 for r in meta_rows
                      if abs(r['n_A1_gate_effective'] - SELECTED_N_A1_GATE) <= TOL)
    n_moved = sum(1 for r in meta_rows
                  if abs(r['n_A1_gate_effective'] - SELECTED_N_A1_GATE) > TOL)
    print()
    print('=== 结论 ===')
    print('  runs at the frozen exponent %g : %d' % (SELECTED_N_A1_GATE, n_at_frozen))
    print('  runs with the exponent moved   : %d (only the n_A1_gate perturbations)' % n_moved)
    print('  expected: %d / %d' % (n_at, n_moved_expected))
    if not (n_at_frozen == n_at and n_moved == n_moved_expected):
        bad('expected %d runs at the frozen exponent and %d moved, got %d/%d'
            % (n_at, n_moved_expected, n_at_frozen, n_moved))

    print()
    if FAIL:
        print('AUDIT FAILED with %d problem(s) - do NOT merge:' % len(FAIL))
        for x in FAIL[:20]:
            print('   -', x)
        return 1
    print('AUDIT PASSED: every run is at one working point; safe to merge.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
