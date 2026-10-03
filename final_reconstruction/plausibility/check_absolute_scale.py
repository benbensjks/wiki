"""Check 3 - absolute scale / burden, with an explicit interface-scale sensitivity.

Two things this script now reports that the first version did not:

  1. FLAGS CARRY MULTIPLES.  `max_free_pool_uM = 42.9` against a reference of
     30 uM is 1.43x, and that reference is itself declared to be an
     order-of-magnitude assumption.  A bare list of "flagged pools" reads like
     seven violations when the honest statement is "1.43x a reference whose own
     uncertainty is a factor of a few".  Every flag now carries its multiple.

  2. SENSITIVITY TO uM_per_au.  The Hill thresholds (K_D_int, K_inh, K_complex,
     K_rep, K_A, K_F) are held in model a.u. and are NOT rescaled when
     uM_per_au changes, so the uM figures are not obviously invariant.  Since
     uM_per_au is the one parameter the whole project marks as an UNCALIBRATED
     interface scale, a single uM number is not enough.  The whole analysis is
     repeated at several uM values and the a.u.-level maxima are compared, which
     separates "the physical number scales with the conversion" from "the a.u.
     dynamics themselves move".

    python check_absolute_scale.py --hours 300
    python check_absolute_scale.py --hours 300 --uM-sweep 5.5,5.75,6.5
"""
from __future__ import annotations

import argparse
import json
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plausibility_common import (COPIES_PER_UM_PER_FL, OUT,  # noqa: E402
                                 PROTEIN_SYNTHESIS_BUDGET_PER_H, REF_BURDEN_FRACTION,
                                 REF_FREE_TF_UM, build_threebit, flux_array,
                                 frozen_carry0, frozen_extension, restore,
                                 write_manifest)
from working_point import working_point_block  # noqa: E402

BIT_OFFSETS = {'bit0': 6, 'bit1': 17, 'bit2': 34}
POOL_NAMES = ['M_I', 'I_u', 'I', 'M_T', 'T_u', 'T', 'M_R', 'R_u', 'R', 'C', 'S']


def analyse(um, hours, sample_min, max_step_min):
    """Run the whole check at one uM_per_au and return the tables + scalars.

    Returns the a.u.-level maxima as well as the converted ones, because the
    a.u. maxima are what show whether the dynamics themselves depend on uM.
    """
    ext = replace(frozen_extension(), uM_per_au=float(um))
    model, patch = build_threebit(carry1=frozen_carry0(), extension=ext)
    try:
        sol = model.simulate(hours=hours, sample_min=sample_min,
                             max_step_min=max_step_min)
        from model import ZENG
        t, y = sol.t, sol.y
        cut = t > 0.5 * hours

        pools = []
        for name, off in BIT_OFFSETS.items():
            for k, sub in enumerate(POOL_NAMES):
                arr = y[off + k][cut]
                pools.append(dict(pool=f'{name}.{sub}',
                                  kind='DNA fraction' if sub == 'S' else 'concentration',
                                  au_max=float(arr.max()), au_median=float(np.median(arr)),
                                  uM_max=float(arr.max()) * um,
                                  copies_per_cell_max=float(arr.max()) * um
                                  * COPIES_PER_UM_PER_FL))
        for name, idx in (('A0', 28), ('F0', 29), ('A1', 45), ('F1', 46)):
            arr = y[idx][cut]
            pools.append(dict(pool=name, kind='concentration', au_max=float(arr.max()),
                              au_median=float(np.median(arr)), uM_max=float(arr.max()) * um,
                              copies_per_cell_max=float(arr.max()) * um
                              * COPIES_PER_UM_PER_FL))
        pool_df = pd.DataFrame(pools)
        conc = pool_df[pool_df.kind == 'concentration'].copy()
        conc['multiple_of_reference'] = conc.uM_max / REF_FREE_TF_UM
        conc['flagged_above_reference'] = conc.uM_max > REF_FREE_TF_UM

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

        src = []
        for i, (name, off) in enumerate(BIT_OFFSETS.items()):
            S, T = y[off + 10], y[off + 5]
            src.append((f'{name}.T', ZENG['alpha_rep'] * (1.0 - S)))
            src.append((f'{name}.RDF', ZENG['alpha_rdf'] * S *
                        (1.0 - (np.maximum(T, 0.0) / ZENG['K_rep']) ** ZENG['n_rep'] /
                         (1.0 + (np.maximum(T, 0.0) / ZENG['K_rep']) ** ZENG['n_rep']))))
        flux = flux_array(model, sol) / um
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
                            copies_per_cell_h_max=float(np.max(v[cut])) * um
                            * COPIES_PER_UM_PER_FL)
                       for name, v in src]
        burden_df = pd.DataFrame(burden_rows)
        total = float(burden_df.copies_per_cell_h_max.sum())
        burden_fraction = total / PROTEIN_SYNTHESIS_BUDGET_PER_H
        wp = working_point_block(model)
        flagged = conc[conc.flagged_above_reference]
        scalars = dict(
            uM_per_au=float(um),
            max_free_pool_au=float(conc.au_max.max()),
            max_free_pool_uM=float(conc.uM_max.max()),
            max_free_pool_multiple_of_reference=float(conc.multiple_of_reference.max()),
            flagged_pool_count=int(len(flagged)),
            flagged_pools=[dict(pool=r.pool, uM_max=float(r.uM_max),
                                multiple_of_reference=float(r.multiple_of_reference))
                           for r in flagged.sort_values('uM_max', ascending=False).itertuples()],
            burden_fraction=float(burden_fraction),
            total_synthesis_copies_per_cell_h=total,
            # a.u.-level fingerprint of the trajectory itself
            au_fingerprint={r.pool: round(float(r.au_max), 12)
                            for r in conc.sort_values('pool').itertuples()})
    finally:
        restore(patch)
    return dict(pool_df=pool_df, conc=conc, thr_rows=thr_rows, burden_df=burden_df,
                scalars=scalars, working_point=wp)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--hours', type=float, default=300.0)
    ap.add_argument('--sample-min', type=float, default=2.0)
    ap.add_argument('--max-step-min', type=float, default=2.0)
    ap.add_argument('--uM', type=float, default=None,
                    help='interface scale for the PRIMARY report (default: the frozen profile value)')
    ap.add_argument('--uM-sweep', type=str, default='5.5,5.75,6.5',
                    help='comma-separated uM values for the sensitivity table')
    args = ap.parse_args()

    frozen_um = float(frozen_extension().uM_per_au)
    report_um = frozen_um if args.uM is None else float(args.uM)
    sweep = sorted({float(x) for x in args.uM_sweep.split(',') if x.strip()} | {report_um})

    primary = analyse(report_um, args.hours, args.sample_min, args.max_step_min)
    pool_df, thr_rows, burden_df = primary['pool_df'], primary['thr_rows'], primary['burden_df']
    conc = primary['conc']

    sens_rows = [primary['scalars']]
    for um in sweep:
        if abs(um - report_um) < 1e-12:
            continue
        sens_rows.append(analyse(um, args.hours, args.sample_min,
                                 args.max_step_min)['scalars'])
    sens_rows.sort(key=lambda r: r['uM_per_au'])
    au_fps = [r['au_fingerprint'] for r in sens_rows]
    au_invariant = all(fp == au_fps[0] for fp in au_fps)

    report = dict(
        hours=args.hours,
        uM_per_au=float(report_um),
        uM_per_au_is_frozen_profile_value=bool(abs(report_um - frozen_um) < 1e-12),
        reference=dict(free_tf_uM=REF_FREE_TF_UM, burden_fraction=REF_BURDEN_FRACTION,
                       protein_synthesis_budget_per_h=PROTEIN_SYNTHESIS_BUDGET_PER_H,
                       reference_caveat=('REF_FREE_TF_UM = 30 uM is declared an '
                                         'order-of-magnitude reference, so a multiple of a '
                                         'few is NOT a violation')),
        max_free_pool_uM=primary['scalars']['max_free_pool_uM'],
        max_free_pool_multiple_of_reference=primary['scalars']['max_free_pool_multiple_of_reference'],
        flagged_pool_count=primary['scalars']['flagged_pool_count'],
        flagged_pools=primary['scalars']['flagged_pools'],
        order_of_magnitude_violation=bool(
            primary['scalars']['max_free_pool_multiple_of_reference'] >= 10.0),
        total_synthesis_copies_per_cell_h=primary['scalars']['total_synthesis_copies_per_cell_h'],
        burden_fraction=primary['scalars']['burden_fraction'],
        burden_flagged=bool(primary['scalars']['burden_fraction'] > REF_BURDEN_FRACTION),
        uM_sensitivity=dict(
            values=sorted(r['uM_per_au'] for r in sens_rows),
            per_value={f'{r["uM_per_au"]:g}': dict(
                max_free_pool_au=r['max_free_pool_au'],
                max_free_pool_uM=r['max_free_pool_uM'],
                max_free_pool_multiple_of_reference=r['max_free_pool_multiple_of_reference'],
                flagged_pool_count=r['flagged_pool_count'],
                burden_fraction=r['burden_fraction']) for r in sens_rows},
            au_level_maxima_invariant=bool(au_invariant),
            reading=('the a.u.-level pool maxima are IDENTICAL across the swept uM values, so '
                     'uM_per_au acts purely as a reporting conversion here and the uM figures '
                     'scale exactly linearly with it'
                     if au_invariant else
                     'the a.u.-level pool maxima MOVE with uM_per_au, so the Hill thresholds '
                     'held in a.u. do shift the dynamics; the uM figures are NOT a pure '
                     'rescaling and the sweep must be quoted alongside them')),
        working_point=primary['working_point'],
    )

    OUT.mkdir(parents=True, exist_ok=True)
    pool_df.to_csv(OUT / 'absolute_scale_pools.csv', index=False, encoding='utf-8')
    pd.DataFrame(thr_rows).to_csv(OUT / 'absolute_scale_thresholds.csv', index=False,
                                  encoding='utf-8')
    burden_df.to_csv(OUT / 'absolute_scale_burden.csv', index=False, encoding='utf-8')
    pd.DataFrame([{k: v for k, v in r.items() if k != 'au_fingerprint'} for r in sens_rows]
                 ).to_csv(OUT / 'absolute_scale_uM_sensitivity.csv', index=False,
                          encoding='utf-8')
    (OUT / 'absolute_scale.json').write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
    write_manifest()

    pd.set_option('display.width', 220)
    print('--- free pools (late-half maximum), uM_per_au = %g ---' % report_um)
    print(conc[['pool', 'au_max', 'uM_max', 'multiple_of_reference', 'flagged_above_reference']]
          .sort_values('uM_max', ascending=False).to_string(index=False))
    print('\n--- uM_per_au sensitivity ---')
    print(pd.DataFrame([{k: v for k, v in r.items() if k != 'au_fingerprint'}
                        for r in sens_rows]).to_string(index=False))
    print('\n  a.u.-level maxima invariant across uM:', au_invariant)
    print(f"\ntotal {report['total_synthesis_copies_per_cell_h']:.3e} copies/cell/h = "
          f"{report['burden_fraction']:.1%} of the assumed budget; "
          f"flagged = {report['burden_flagged']}")
    print('\n  max free pool = %.2f uM = %.2fx the %g uM reference; order-of-magnitude '
          'violation = %s' % (report['max_free_pool_uM'],
                              report['max_free_pool_multiple_of_reference'],
                              REF_FREE_TF_UM, report['order_of_magnitude_violation']))
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
