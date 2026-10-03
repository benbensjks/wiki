"""Symmetric version probe (B) for the frozen 51-state counter.

NEW FILE (added 2026-09-28).  It does NOT modify any existing script, model
file, frozen artefact or manifest.  Every output goes under
`跨模型三级连接验证/版本与鲁棒性探针/results/`.

WHY THIS PROBE EXISTS
---------------------
The 34-state hybrid model already carries ONE half of this axis: the arm in
`hby_zmh_hby/pdf_parameters/` puts the EARLY table values back into the middle
block, and the count collapses (29/5/1 instead of 29/14/7).  The 51-state model
has no symmetric arm -- it runs entirely on the early `ZENG` table -- so on its
own it cannot say whether its verdict is robust or conditional on that version.

This script runs the missing half: at the FROZEN 51-state working point, swap in
the collaborator's CURRENT `wiki任务/完整二级级联.py` values, for the entries
that are genuinely 1:1 counterparts, and report whether the verdict and the
binding timing margin survive.

PRE-REGISTERED MAPPING (fixed before the runs; see MAPPED / NOT_MAPPED below)
---------------------------------------------------------------------------
Her file is in min^-1, the 51-state table is in h^-1, so every rate is x60.

MAPPED (1:1 counterparts)
  * k_fwd 0.118 -> 7.08 /h        (table 7.0)
  * k_rev 0.08  -> 4.8  /h        (table 5.0)
  * K_rep 1.8, n_rep 3.0          (table 0.85 / 3.9)
  * alpha_rep, gamma_rep          per-response: ZENG holds ONE shared slot for
                                  all three bits, her file gives DIFFERENT values
                                  for bit0 and bit1, so this script runs BOTH
                                  variants instead of silently picking one
  * alpha_rdf, gamma_rdf          same shared-slot ambiguity, same two variants
  * K_D_int[i] 1.8                her K_D_int0 == K_D_int1 == 1.8; index 2 has
                                  no counterpart and is LEFT AT THE TABLE VALUE
  * gamma_int[1] = 0.08 -> 4.8 /h her `gamma_int1`; index 0 and 2 are NOT mapped
                                  (see below)

NOT MAPPED, and why (each is declared, not silently dropped)
  * K_D_comp 3.2      Her threshold acts on the PRODUCT int*rdf; in this model
                      the same chemistry is converted through
                      q = kon/(koff + complex_decay + growth) into an explicit
                      complex threshold K_complex = q*K_D_comp = 0.06 acting on
                      state C.  Different physical quantities.
  * gamma_int[0]      Her Int0 is a square-wave-driven species cleared at
                      gamma_int = 2.0 /min (=120 /h); this model's bit0 `I` is
                      the mature free integrase produced from the real C31
                      translation flux.  Not the same species.
  * gamma_int[2]      Her file has no bit2.
  * K_D_int[2]        Her file has no bit2.
  * Kinh 0.1026       Form is 1:1 identical, but it differs from the table value
                      0.1 by 2.6 %; left alone so the probe changes one thing at
                      a time.
  * alpha_A0/gamma_A0/alpha_R0/gamma_R0/n_R0/K_R0
                      Her A0/R0 block is the I1-FFL; this model's A0/F0 arm is
                      built from alpha_A/gamma_A/K_A/n_A with F driven through
                      act(A, K_A, n_A).  Not a 1:1 map.
  * alpha_Int1 0.3788 Her Int1 is an a.u.-scale species with its own
                      gamma_int1; this model's alpha_Int[1] feeds a translation
                      cascade.  Units and role differ.
  * k_int 6.0         Not referenced by this model's rhs at all.

PRE-REGISTERED ACCEPTANCE (fixed before the runs)
  primary   : `certified` (the frozen predicate, unchanged) stays True
  secondary : steady-state read window count, minimum commitment, and the
              per-bit setup/hold margins, all compared against the baseline run
              executed in the SAME process (never against a re-typed number)
  reported  : gate contrast, event counts, code sequence, elapsed time
No threshold is invented here: the script reports margins as numbers.

INTEGRITY GUARDS
  1. the frozen working point is asserted (clock 0.3/2.0, n_A1_gate 6, uM 5.75)
  2. the baseline run must reproduce `certified = True`, else the script aborts
     rather than comparing against a mis-specified working point
  3. after restoring ZENG, rhs(0, y0) must be bit-identical to the pre-patch
     snapshot (H2: an un-restored global patch would contaminate later runs)
  4. one model instance is alive at a time
"""
from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')

WIKI = Path(r'C:\Users\18633\Desktop\wiki')
FR = WIKI / 'final_reconstruction'
HERE = Path(__file__).resolve().parent
OUT = HERE / 'results'
ZMH_PY = WIKI / 'wiki任务' / '完整二级级联.py'

sys.path.insert(0, str(FR))
sys.path.insert(0, str(FR / 'plausibility'))

import model as M                                   # noqa: E402  (frozen, imported only)
import plausibility_common as PC                    # noqa: E402
from verify_threebit51 import analyse_threebit      # noqa: E402

HOURS = 600.0
SAMPLE_MIN = 2.0
MIN_PER_H = 60.0
EXPECTED_CLOCK = (0.3, 2.0)
EXPECTED_N_GATE = 6.0
EXPECTED_UM_PER_AU = 5.75

# ---- the two variants that resolve the shared-slot ambiguity (values are min^-1)
VARIANT_BLOCKS = {
    'zmh_bit1_block': dict(alpha_rep=0.2, gamma_rep=0.005,
                           alpha_rdf=0.0265, gamma_rdf=0.005),
    'zmh_bit0_block': dict(alpha_rep=0.1821, gamma_rep=0.0119,
                           alpha_rdf=0.0325, gamma_rdf=0.0055),
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest().upper()


def jsonable(o):
    if isinstance(o, dict):
        return {str(k): jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [jsonable(v) for v in o]
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (bool, int, float, str)) or o is None:
        return o
    return str(o)


def build_patch(block_key):
    """The 1:1 map, expressed as absolute ZENG keys (h^-1)."""
    b = VARIANT_BLOCKS[block_key]
    return {
        'k_fwd': 0.118 * MIN_PER_H,
        'k_rev': 0.08 * MIN_PER_H,
        'K_rep': 1.8,
        'n_rep': 3.0,
        'alpha_rep': b['alpha_rep'] * MIN_PER_H,
        'gamma_rep': b['gamma_rep'] * MIN_PER_H,
        'alpha_rdf': b['alpha_rdf'] * MIN_PER_H,
        'gamma_rdf': b['gamma_rdf'] * MIN_PER_H,
        'K_D_int': (1.8, 1.8, float(M.ZENG['K_D_int'][2])),
        'gamma_int': (float(M.ZENG['gamma_int'][0]), 0.08 * MIN_PER_H,
                      float(M.ZENG['gamma_int'][2])),
    }


def assert_working_point(model):
    problems = []
    if (model.e.clock_K_au, model.e.clock_n) != EXPECTED_CLOCK:
        problems.append(f'clock gate is {(model.e.clock_K_au, model.e.clock_n)}, '
                        f'expected {EXPECTED_CLOCK}')
    if float(model.n_A1_gate_effective) != EXPECTED_N_GATE:
        problems.append(f'n_A1_gate_effective = {model.n_A1_gate_effective}')
    if float(model.e.uM_per_au) != EXPECTED_UM_PER_AU:
        problems.append(f'uM_per_au = {model.e.uM_per_au}')
    if model.e.add_growth:
        problems.append('add_growth is True')
    if problems:
        raise RuntimeError('working-point assertion failed: ' + '; '.join(problems))


def margins_of(verdict):
    out = {}
    for bit, m in verdict['bit_margins'].items():
        out[bit] = {'min_setup_h': m.get('min_setup_h'),
                    'min_hold_h': m.get('min_hold_h')}
    holds = [v['min_hold_h'] for v in out.values() if v['min_hold_h'] is not None]
    setups = [v['min_setup_h'] for v in out.values() if v['min_setup_h'] is not None]
    out['binding_min_hold_h'] = min(holds) if holds else None
    out['binding_min_setup_h'] = min(setups) if setups else None
    return out


def run_case(tag, patch):
    """Build one 51-state model with `patch` in force, integrate, verify."""
    original = {k: M.ZENG[k] for k in patch}
    row = {'tag': tag, 'patch': jsonable(patch), 'requested_keys': sorted(patch)}
    started = time.perf_counter()
    try:
        for key, value in patch.items():
            M.ZENG[key] = value
        model, residual = PC.build_threebit(n_A1_gate=EXPECTED_N_GATE)
        if residual:
            raise RuntimeError(f'unexpected residual patch {residual}')
        assert_working_point(model)
        sol = model.simulate(hours=HOURS, sample_min=SAMPLE_MIN,
                             rtol=2e-7, atol=2e-9, max_step_min=2.0)
        verdict = analyse_threebit(model, sol, HOURS)
    finally:
        PC.restore(original)
    row['runtime_s'] = time.perf_counter() - started
    cold, steady = verdict['cold_start'], verdict['steady_state']
    causal = verdict['causal_verdict']
    contrast = verdict['gate_contrast']
    row.update(
        certified=bool(verdict['certified']),
        steady_passed=bool(steady['passed']),
        steady_reads=int(steady['reads']),
        steady_sequence=steady['sequence'],
        minimum_commitment=float(steady['minimum_commitment']),
        cold_reads=int(cold['reads']),
        cold_sequence=cold['sequence'],
        carry1_gate_events=len(verdict['carry1_gate_events']),
        bit1_reverse_events=len(verdict['bit1_reverse_events']),
        bit2_crossings=len(verdict['bit2_crossings']),
        exactly_one=bool(causal['exactly_one_gate_and_flip_per_late_reverse']),
        causal_order=bool(causal['causal_order_after_reverse_start']),
        alternating=bool(causal['bit2_directions_alternate']),
        off_to_on_peak_ratio=contrast['off_peak_to_on_peak'],
        contrast_passed=bool(contrast['passed']),
        on_cycle_peak_median=contrast['on_cycle_peak_median'],
        margins=margins_of(verdict),
        signal_ranges=verdict['signal_ranges'],
    )
    return row, verdict


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    meta = {
        'script': Path(__file__).name,
        'script_sha256': sha256(Path(__file__).resolve()),
        'hours': HOURS,
        'sample_min': SAMPLE_MIN,
        'solver': dict(method='DOP853', rtol=2e-7, atol=2e-9, max_step_min=2.0),
        'working_point': dict(clock_K_au=EXPECTED_CLOCK[0], clock_n=EXPECTED_CLOCK[1],
                              n_A1_gate=EXPECTED_N_GATE, uM_per_au=EXPECTED_UM_PER_AU,
                              add_growth=False,
                              source='plausibility_common.frozen_extension() + '
                                     'threebit51_results/wiki_n6_20260924/parameters.json'),
        'zeng_baseline': jsonable({k: v for k, v in M.ZENG.items()}),
        'variant_blocks_min_inv': VARIANT_BLOCKS,
        'source_sha256': PC.source_hashes(),
    }
    meta['source_sha256']['zmh_complete_two_bit_py'] = sha256(ZMH_PY)

    # guard 3: a pre-patch rhs snapshot, used again after every restore
    ref_model, ref_patch = PC.build_threebit(n_A1_gate=EXPECTED_N_GATE)
    assert not ref_patch
    assert_working_point(ref_model)
    y_ref = ref_model.initial_state()
    rhs_ref = np.asarray(ref_model.rhs(0.0, y_ref)).copy()
    del ref_model

    cases = [('baseline_frozen_zeng', {})]
    for key in VARIANT_BLOCKS:
        cases.append((key, build_patch(key)))

    rows, verdicts, restore_gaps = [], {}, {}
    for tag, patch in cases:
        print(f'[run] {tag} ...', flush=True)
        row, verdict = run_case(tag, patch)
        rows.append(row)
        verdicts[tag] = jsonable(verdict)
        print(f'      certified={row["certified"]} reads={row["steady_reads"]} '
              f'binding_hold={row["margins"]["binding_min_hold_h"]} '
              f'({row["runtime_s"]:.1f} s)', flush=True)
        check, check_patch = PC.build_threebit(n_A1_gate=EXPECTED_N_GATE)
        assert not check_patch
        gap = float(np.max(np.abs(np.asarray(check.rhs(0.0, y_ref)) - rhs_ref)))
        restore_gaps[tag] = gap
        del check

    # guard 2: the baseline must reproduce the frozen verdict
    baseline = rows[0]
    if not baseline['certified']:
        raise RuntimeError('baseline is not certified at the asserted working point; '
                           'refusing to compare a version probe against it')

    result = dict(meta=meta, rows=rows, verdicts=verdicts,
                  restore_gap_after_each_run=restore_gaps,
                  h2_note='rhs(0, y0) recomputed after every restore; a non-zero gap '
                          'means a ZENG patch leaked into a later model instance')
    (OUT / 'probe_B_zeng_version_51state.json').write_text(
        json.dumps(jsonable(result), ensure_ascii=False, indent=1), encoding='utf-8')

    hdr = ['tag', 'certified', 'steady_reads', 'minimum_commitment', 'binding_min_hold_h',
           'binding_min_setup_h', 'carry1_gate_events', 'bit2_crossings',
           'exactly_one', 'causal_order', 'alternating', 'off_to_on_peak_ratio',
           'runtime_s']
    lines = [','.join(hdr)]
    for r in rows:
        lines.append(','.join(str(x) for x in [
            r['tag'], r['certified'], r['steady_reads'], r['minimum_commitment'],
            r['margins']['binding_min_hold_h'], r['margins']['binding_min_setup_h'],
            r['carry1_gate_events'], r['bit2_crossings'], r['exactly_one'],
            r['causal_order'], r['alternating'], r['off_to_on_peak_ratio'],
            r['runtime_s']]))
    (OUT / 'probe_B_zeng_version_51state.csv').write_text('\n'.join(lines) + '\n',
                                                          encoding='utf-8')
    print('\n'.join(lines))
    print('restore gaps:', restore_gaps)
    print('written to', OUT)


if __name__ == '__main__':
    main()
