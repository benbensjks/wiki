"""Check 3 - absolute scale and expression burden.

Question
--------
The model is written in "a.u." with a.u. Hill thresholds.  The only bridge to
real units is Extension.uM_per_au, and by construction
`copies_per_au = 602.214 * uM_per_au`, i.e. 1 a.u. = uM_per_au micromolar.
At the frozen working point uM_per_au = 5.75, so every a.u. number implies a
physical concentration.  Nobody has ever printed that table.

This script prints it and asks two plausibility questions:

 1. Are the implied absolute levels (and the implied Hill thresholds, which are
    also in a.u.) inside a defensible range for a bacterial cell?
 2. Does the circuit's total protein synthesis rate fit inside the cell's
    translation budget, or does the design demand an implausible share of it?

Pre-registered flagging rules
-----------------------------
FLAG a free pool whose late-run maximum exceeds REF_FREE_TF_UM (30 uM).
FLAG the design if the summed synthesis demand exceeds REF_BURDEN_FRACTION
(10 %) of PROTEIN_SYNTHESIS_BUDGET_PER_H.
Both reference numbers are order-of-magnitude assumptions and are printed with
the result so a reader can re-scale them.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plausibility_common import (COPIES_PER_UM_PER_FL, OUT,  # noqa: E402
                                 PROTEIN_SYNTHESIS_BUDGET_PER_H, REF_BURDEN_FRACTION,
                                 REF_FREE_TF_UM, build_threebit, flux_array,
                                 frozen_carry0, restore, write_manifest)

BIT_OFFSETS = {'bit0': 6, 'bit1': 17, 'bit2': 34}
POOL_NAMES = ['M_I', 'I_u', 'I', 'M_T', 'T_u', 'T', 'M_R', 'R_u', 'R', 'C', 'S']


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--hours', type=float, default=300.0)
    ap.add_argument('--sample-min', type=float, default=2.0)
    ap.add_argument('--max-step-min', type=float, default=2.0)
    args = ap.parse_args()

    model, patch = build_threebit(carry1=frozen_carry0())
    try:
        sol = model.simulate(hours=args.hours, sample_min=args.sample_min,
                             max_step_min=args.max_step_min)
    finally:
        restore(patch)

    from model import ZENG
    t, y = sol.t, sol.y
    cut = t > 0.5 * args.hours
    um = model.e.uM_per_au

    # ---------------------------------------------------------- pool table
    pools = []
    for name, off in BIT_OFFSETS.items():
        for k, sub in enumerate(POOL_NAMES):
            arr = y[off + k][cut]
            pools.append(dict(pool=f'{name}.{sub}', kind='DNA fraction' if sub == 'S' else 'concentration',
                              au_max=float(arr.max()), au_median=float(np.median(arr)),
                              uM_max=float(arr.max()) * um,
                              copies_per_cell_max=float(arr.max()) * um * COPIES_PER_UM_PER_FL))
    for name, idx in (('A0', 28), ('F0', 29), ('A1', 45), ('F1', 46)):
        arr = y[idx][cut]
        pools.append(dict(pool=name, kind='concentration', au_max=float(arr.max()),
                          au_median=float(np.median(arr)), uM_max=float(arr.max()) * um,
                          copies_per_cell_max=float(arr.max()) * um * COPIES_PER_UM_PER_FL))
    pool_df = pd.DataFrame(pools)
    conc = pool_df[pool_df.kind == 'concentration'].copy()
    conc['flagged_above_30uM'] = conc.uM_max > REF_FREE_TF_UM

    # ------------------------------------------------------ threshold table
    thr = [('K_D_int bit0/1/2', ZENG['K_D_int']),
           ('K_inh', (ZENG['K_inh'],)),
           ('K_complex (derived q*1.2)', (model.base.K_complex,)),
           ('K_rep', (ZENG['K_rep'],)),
           ('K_A0 / K_A1', ZENG['K_A']),
           ('K_F0 / K_F1', ZENG['K_F']),
           ('K_auto1', (ZENG['K_auto1'],)),
           ('clock_K', (model.e.clock_K_au,))]
    thr_rows = [dict(threshold=name, au=float(v), uM=float(v) * um,
                     copies_per_cell=float(v) * um * COPIES_PER_UM_PER_FL)
                for name, values in thr for v in values]

    # ----------------------------------------------------- synthesis burden
    src = []
    for i, (name, off) in enumerate(BIT_OFFSETS.items()):
        S, T = y[off + 10], y[off + 5]
        src.append((f'{name}.T', ZENG['alpha_rep'] * (1.0 - S)))
        src.append((f'{name}.RDF', ZENG['alpha_rdf'] * S *
                    (1.0 - (np.maximum(T, 0.0) / ZENG['K_rep']) ** ZENG['n_rep'] /
                     (1.0 + (np.maximum(T, 0.0) / ZENG['K_rep']) ** ZENG['n_rep']))))
    flux = flux_array(model, sol) / um                     # a.u./h
    sig = model.diagnostic_signals(y)
    from plausibility_common import hill as vhill, repression as vrep
    src.append(('bit0.Int (real C31 flux)', flux))
    src.append(('bit1.Int (alpha_Int0 * g0)', ZENG['alpha_Int'][0] * sig['g0']))
    src.append(('bit2.Int (alpha_Int1 * g1)', ZENG['alpha_Int'][1] * sig['g1']))
    src.append(('A0', ZENG['alpha_A'][0] * (1.0 - y[16])))
    src.append(('F0', ZENG['alpha_F'][0] * vhill(y[28], ZENG['K_A'][0], ZENG['n_A'][0])))
    src.append(('A1', ZENG['alpha_A'][1] * (1.0 - y[27]) *
                vrep(y[45], ZENG['K_auto1'], ZENG['n_auto1'])))
    src.append(('F1', ZENG['alpha_F'][1] * vhill(y[45], ZENG['K_A'][1], ZENG['n_A'][1])))
    burden_rows = [dict(species=name, au_per_h_max=float(np.max(v[cut])),
                        uM_per_h_max=float(np.max(v[cut])) * um,
                        copies_per_cell_h_max=float(np.max(v[cut])) * um * COPIES_PER_UM_PER_FL)
                   for name, v in src]
    burden_df = pd.DataFrame(burden_rows)
    total = float(burden_df.copies_per_cell_h_max.sum())
    burden_fraction = total / PROTEIN_SYNTHESIS_BUDGET_PER_H

    report = dict(hours=args.hours, uM_per_au=float(um),
                  reference=dict(free_tf_uM=REF_FREE_TF_UM,
                                 burden_fraction=REF_BURDEN_FRACTION,
                                 protein_synthesis_budget_per_h=PROTEIN_SYNTHESIS_BUDGET_PER_H),
                  flagged_pools=conc[conc.flagged_above_30uM].pool.tolist(),
                  max_free_pool_uM=float(conc.uM_max.max()),
                  total_synthesis_copies_per_cell_h=total,
                  burden_fraction=float(burden_fraction),
                  burden_flagged=bool(burden_fraction > REF_BURDEN_FRACTION))
    OUT.mkdir(parents=True, exist_ok=True)
    pool_df.to_csv(OUT / 'absolute_scale_pools.csv', index=False, encoding='utf-8')
    pd.DataFrame(thr_rows).to_csv(OUT / 'absolute_scale_thresholds.csv', index=False, encoding='utf-8')
    burden_df.to_csv(OUT / 'absolute_scale_burden.csv', index=False, encoding='utf-8')
    (OUT / 'absolute_scale.json').write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
    write_manifest()

    pd.set_option('display.width', 220)
    print('--- free pools (late-half maximum) ---')
    print(conc[['pool', 'au_max', 'uM_max', 'copies_per_cell_max', 'flagged_above_30uM']]
          .sort_values('uM_max', ascending=False).to_string(index=False))
    print('\n--- Hill thresholds in physical units ---')
    print(pd.DataFrame(thr_rows).to_string(index=False))
    print('\n--- synthesis demand ---')
    print(burden_df.sort_values('copies_per_cell_h_max', ascending=False).to_string(index=False))
    print(f"\ntotal {total:.3e} copies/cell/h = {burden_fraction:.1%} of the assumed budget "
          f"({PROTEIN_SYNTHESIS_BUDGET_PER_H:.1e}); flagged = {report['burden_flagged']}")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
