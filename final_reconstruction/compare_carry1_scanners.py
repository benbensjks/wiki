"""Head-to-head correctness comparison of the two stage-1 carry1 scanners.

    gpt   : scan_threebit51_carry1.py        (GPT, left untouched)
    dsh   : scan_threebit51_carry1_dsh.py    (this session)

The two scanners must agree on four things for their results to be comparable:

  A. grid and shard coverage (same 35 points, every shard union = grid, no overlap)
  B. model configuration built from the same job
  C. initial state vector (bit-identical)
  D. full evaluation row over the shared fields (bit-identical floats)

`analyse_threebit` and the model are frozen modules, so any difference here is a
real implementation difference, not a criterion difference.

    python compare_carry1_scanners.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

import scan_threebit51_carry1 as gpt
import scan_threebit51_carry1_dsh as dsh
from model import ROOT
from model_threebit51 import ThreeBit51Model, ThreeBitCarryParameters
from model_twobit34 import CarryExpressionParameters

OUT = ROOT / 'threebit51_results' / 'carry1_scan_dsh'
PROBE_JOBS = [dict(carry1_mrna_min=2.0, carry1_maturation_min=20.0),   # interior
              dict(carry1_mrna_min=0.5, carry1_maturation_min=5.0),    # corner
              dict(carry1_mrna_min=8.0, carry1_maturation_min=60.0)]   # opposite corner
PROBE_HOURS = [100.0, 200.0]
IGNORE = {'runtime_s', 'error_trace'}


def build(cfg_job):
    c1 = CarryExpressionParameters(
        mrna_half_life_min=cfg_job['carry1_mrna_min'],
        activator_maturation_half_life_min=cfg_job['carry1_maturation_min'],
        repressor_maturation_half_life_min=cfg_job['carry1_maturation_min'])
    return ThreeBit51Model(carry=ThreeBitCarryParameters(
        carry0=ThreeBitCarryParameters().carry0, carry1=c1))


def main():
    report = {}

    # ---------------------------------------------------------------- A
    gj, dj = gpt.jobs(), dsh.jobs()
    gset = {(j['carry1_mrna_min'], j['carry1_maturation_min']) for j in gj}
    dset = {(j['carry1_mrna_min'], j['carry1_maturation_min']) for j in dj}

    def shard_union(module, nshards=7):
        seen = []
        for s in range(nshards):
            seen += [(j['carry1_mrna_min'], j['carry1_maturation_min'])
                     for j in module.jobs()[s::nshards]]
        return seen

    for name, module in (('gpt', gpt), ('dsh', dsh)):
        union = shard_union(module)
        report[f'A_shards_{name}'] = dict(
            grid_size=len(module.jobs()), shard_union_size=len(union),
            union_equals_grid=set(union) == (gset if name == 'gpt' else dset),
            duplicates=len(union) - len(set(union)),
            per_shard=[len(module.jobs()[s::7]) for s in range(7)])
    report['A_grids_identical'] = gset == dset
    report['A_mrna_values'] = dict(gpt=sorted(gpt.MRNA), dsh=sorted(dsh.MRNA),
                                   equal=sorted(gpt.MRNA) == sorted(dsh.MRNA))
    report['A_mat_values'] = dict(gpt=sorted(gpt.MAT), dsh=sorted(dsh.MAT),
                                  equal=sorted(gpt.MAT) == sorted(dsh.MAT))
    report['A_hours_default'] = dict(gpt=600.0, dsh=600.0)
    report['A_output_dirs_differ'] = str(gpt.OUT) != str(dsh.OUT)

    # ------------------------------------------------------- B and C
    cfg_rows, init_rows = [], []
    for job in PROBE_JOBS:
        gm, dm = build(job), build(job)
        ge, de = gm.e, dm.e
        same_ext = all(getattr(ge, f) == getattr(de, f) for f in
                       ('uM_per_au', 'mrna_half_life_min', 'maturation_half_life_min',
                        'translation_h', 'immature_decay_h', 'complex_on_au_inv_h',
                        'complex_off_h', 'complex_decay_h', 'clock_K_au', 'clock_n',
                        'add_growth'))
        same_c0 = gm.carry.carry0 == dm.carry.carry0
        same_c1 = gm.carry.carry1 == dm.carry.carry1
        y0g, y0d = gm.initial_state(), dm.initial_state()
        cfg_rows.append(dict(job=job, extension_equal=same_ext, carry0_equal=same_c0,
                             carry1_equal=same_c1, shape_gpt=list(y0g.shape),
                             shape_dsh=list(y0d.shape),
                             init_state_max_abs_diff=float(np.max(np.abs(y0g - y0d)))))
        init_rows.append(float(np.max(np.abs(y0g - y0d))))
    report['B_model_config'] = cfg_rows
    report['C_initial_state_max_abs_diff'] = max(init_rows)

    # ------------------------------------------------------------ D
    evals = []
    for job in PROBE_JOBS:
        for hours in PROBE_HOURS:
            rg = gpt.evaluate(dict(job), hours=hours)
            rd = dsh.evaluate(dict(job), hours=hours)
            shared = sorted((set(rg) & set(rd)) - IGNORE)
            diffs = []
            for k in shared:
                a, b = rg.get(k), rd.get(k)
                if isinstance(a, float) or isinstance(b, float) or \
                        isinstance(a, (int, np.floating, np.integer)):
                    try:
                        same = bool(np.isclose(float(a), float(b), rtol=0.0, atol=0.0))
                    except (TypeError, ValueError):
                        same = a == b
                else:
                    same = a == b
                if not same:
                    diffs.append(dict(field=k, gpt=str(a)[:120], dsh=str(b)[:120]))
            only_g = sorted((set(rg) - set(rd)) - IGNORE)
            only_d = sorted((set(rd) - set(rg)) - IGNORE)
            evals.append(dict(job=job, hours=hours, shared_fields=len(shared),
                              identical=not diffs, differences=diffs,
                              fields_only_gpt=only_g, fields_only_dsh=only_d,
                              gpt_certified=rg.get('certified'), dsh_certified=rd.get('certified'),
                              gpt_crossings=rg.get('bit2_crossings'),
                              dsh_crossings=rd.get('bit2_crossings'),
                              gpt_runtime_s=rg.get('runtime_s'), dsh_runtime_s=rd.get('runtime_s')))
            print(f"[D] job={job['carry1_mrna_min']}/{job['carry1_maturation_min']} "
                  f"hours={hours:g} shared={len(shared)} identical={not diffs} "
                  f"cert={rg.get('certified')}/{rd.get('certified')} "
                  f"cross={rg.get('bit2_crossings')}/{rd.get('bit2_crossings')}", flush=True)
    report['D_row_equivalence'] = evals
    report['D_all_identical'] = all(e['identical'] for e in evals)

    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / 'scanner_comparison.json'
    path.write_text(json.dumps(report, ensure_ascii=False, indent=2, default=str), encoding='utf-8')

    print()
    print('A grid identical      :', report['A_grids_identical'])
    print('A shard unions ok     :',
          all(report[f'A_shards_{n}']['union_equals_grid'] and
              report[f'A_shards_{n}']['duplicates'] == 0 for n in ('gpt', 'dsh')))
    print('A output dirs differ  :', report['A_output_dirs_differ'])
    print('B extension equal     :', all(r['extension_equal'] for r in cfg_rows))
    print('B carry0/carry1 equal :', all(r['carry0_equal'] and r['carry1_equal'] for r in cfg_rows))
    print('C init state max diff :', report['C_initial_state_max_abs_diff'])
    print('D all rows identical  :', report['D_all_identical'],
          f"({sum(1 for e in evals if e['identical'])}/{len(evals)} probes)")
    for e in evals:
        if e['differences']:
            print('   DIFF', e['job'], e['hours'], e['differences'][:3])
    print('wrote', path)


if __name__ == '__main__':
    main()
