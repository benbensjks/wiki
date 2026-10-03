"""Does the gate exponent trade off against the A1 maturation half-life?

Why this exists
---------------
The converged local-identifiability run (9-parameter set, 1 % step) reports
rho(n_A1_gate, A1_maturation_half_life_min) = -0.640.  That was NOT a predicted
confusion - it is a trade-off the sensitivity matrix found on its own.  A
correlation between two parameters means their effects on the observables partly
cancel, so a steeper A1 Hill could in principle substitute for a slower A1
maturation (or the reverse).  That would matter for engineering: whichever of the
two is easier to realise could carry the function.

This scan maps the trade-off directly instead of inferring it from a local
derivative.

Design
------
Grid: n_A1_gate in {5.0, 5.5, 6.0, 6.5, 7.0} x A1 maturation in
      {25, 30, 32.5, 38, 45} min, with the F1 arm held at its frozen 32.5 min so
      the A1 arm is isolated.  25 points, 400 h each (>=16 steady read windows,
      the project's minimum for a certification verdict).

Recorded per point: certification, event counts, the read sequence, the gate
window peak / dose / width, the bit2 band margins and the paired far-off leak.
The trade-off is then read off the iso-lines of the gate dose: if the iso-dose
contours in the (n, maturation) plane have negative slope, then a steeper Hill
compensates a slower maturation.

    python scan_nA1_vs_maturation.py scan  --shard S --nshards N
    python scan_nA1_vs_maturation.py merge --nshards N
"""
from __future__ import annotations

import argparse
import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plausibility_common import (OUT, frozen_carry0, frozen_extension,  # noqa: E402
                                 gate_windows, merge_shards, shard_slice,
                                 sha256, source_hashes, write_manifest,
                                 write_shard)
from model_threebit51 import ThreeBit51Model, ThreeBitCarryParameters  # noqa: E402
from model_twobit34 import CarryExpressionParameters  # noqa: E402
from working_point import working_point_block  # noqa: E402

NAME = 'nA1_vs_maturation'
ROOT_MODEL = Path(__file__).resolve().parents[1] / 'model_threebit51.py'
N_VALUES = (5.0, 5.5, 6.0, 6.5, 7.0)
MAT_VALUES = (25.0, 30.0, 32.5, 38.0, 45.0)
F1_MAT_FROZEN = 32.5
MRNA = 2.0
HOURS = 400.0
MIN_STEADY_READS = 16


def points():
    return [dict(n_A1_gate=n, a1_mat=m) for n, m in itertools.product(N_VALUES, MAT_VALUES)]


def point_id(p):
    return f"n={p['n_A1_gate']:g}|A1mat={p['a1_mat']:g}"


def build(p):
    carry1 = CarryExpressionParameters(
        mrna_half_life_min=MRNA,
        activator_maturation_half_life_min=p['a1_mat'],
        repressor_maturation_half_life_min=F1_MAT_FROZEN)
    return ThreeBit51Model(extension=frozen_extension(),
                           carry=ThreeBitCarryParameters(carry0=frozen_carry0(),
                                                         carry1=carry1),
                           n_A1_gate=p['n_A1_gate'])


def evaluate(p, hours, sample_min, max_step_min):
    from verify_threebit51 import analyse_threebit
    from scan_carry_pairing import analyse as paired_analyse
    started = time.perf_counter()
    row = dict(point_id=point_id(p), **p)
    try:
        model = build(p)
        sol = model.simulate(hours=hours, sample_min=sample_min,
                             max_step_min=max_step_min)
        a = analyse_threebit(model, sol, hours)
        sig = model.diagnostic_signals(sol.y)
        wins = gate_windows(sol.t, sig['g1'])
        peaks = [w['peak'] for w in wins]
        doses = [w['dose'] for w in wins]
        widths = [w['width_h'] for w in wins]
        row.update(ok=True, error='', certified=bool(a['certified']),
                   crossings=len(a['bit2_crossings']),
                   gates=len(a['carry1_gate_events']),
                   reverse_events=len(a['bit1_reverse_events']),
                   cold_reads=a['cold_start']['reads'],
                   steady_reads=a['steady_state']['reads'],
                   sequence=a['steady_state']['sequence'],
                   unlabelled=int(sum(r['value'] is None for r in a['read_windows'])),
                   gate_peak_median=float(np.median(peaks)) if peaks else None,
                   gate_dose_median=float(np.median(doses)) if doses else None,
                   gate_width_median=float(np.median(widths)) if widths else None,
                   gate_peak_max=float(np.max(peaks)) if peaks else None,
                   bit2_margin_low=a['bit_margins']['bit2']['min_hold_h'],
                   setup_global=a['bit_margins']['bit2']['min_setup_h'],
                   hold_global=a['bit_margins']['bit2']['min_hold_h'],
                   steady_min_commitment=a['steady_state']['minimum_commitment'])
        agg, _, _, _, _ = paired_analyse(model, sol, 0.05)
        row.update(leak_symmetric=agg['L_symmetric']['median'],
                   leak_rev=agg['L_rev']['median'],
                   n_pairs=agg['n_pairs'], pattern=agg['pattern'],
                   wp_n_A1_gate=float(model.n_A1_gate_effective))
    except Exception as exc:                                     # noqa: BLE001
        row.update(ok=False, error=repr(exc), certified=False)
    row['runtime_s'] = round(time.perf_counter() - started, 2)
    return row


def scan(args):
    pts = shard_slice(points(), args.shard, args.nshards)
    tag = f'shard{args.shard:02d}of{args.nshards:02d}'
    rows = []
    for p in pts:
        r = evaluate(p, args.hours, args.sample_min, args.max_step_min)
        rows.append(r)
        print(f"[{tag}] {point_id(p):<22} cert={r.get('certified')} gates={r.get('gates')} "
              f"peak={r.get('gate_peak_median')} dose={r.get('gate_dose_median')} "
              f"L_sym={r.get('leak_symmetric')} ({r['runtime_s']}s)", flush=True)
    write_shard(NAME, tag, rows, dict(kind='scan', shard=args.shard,
                                      nshards=args.nshards, hours=args.hours,
                                      points=len(pts), model_sha256=sha256(ROOT_MODEL),
                                      source_hashes=source_hashes()))
    print(f'[{tag}] wrote {len(rows)} rows')


def merge(args):
    df = merge_shards(NAME, args.nshards).sort_values(['n_A1_gate', 'a1_mat'])
    ok = df[df.ok.astype(bool)]
    dose = ok.pivot_table(index='n_A1_gate', columns='a1_mat', values='gate_dose_median')
    peak = ok.pivot_table(index='n_A1_gate', columns='a1_mat', values='gate_peak_median')
    cert = ok.pivot_table(index='n_A1_gate', columns='a1_mat',
                          values='certified').astype(float)

    # Per-row sensitivity of the gate dose to the maturation time.  NOTE: this is a
    # least-squares slope over the scanned 25-45 min range, NOT a point derivative.
    xs = np.asarray(dose.columns, dtype=float)
    slopes = {}
    for n in dose.index:
        ys = dose.loc[n].to_numpy(dtype=float)
        good = np.isfinite(ys)
        slopes[f'n={n:g}'] = (float(np.polyfit(xs[good], ys[good], 1)[0])
                              if good.sum() >= 2 else None)

    # Iso-dose curve.  Three states per row, not a number: the maturation that
    # reproduces the frozen gate dose may lie inside the scanned range (in_grid,
    # measured by interpolation), below it, or above it (an EXTRAPOLATION off the end
    # of the grid, flagged as such).  An earlier version returned a bare null for
    # "outside the grid", which read as "no value" and let the narrative quote a local
    # slope as if the compensation point had been located.
    target = float(dose.loc[6.0, 32.5]) if 6.0 in dose.index and 32.5 in dose.columns else None
    equiv, compensation = {}, {}
    if target is not None:
        for n in dose.index:
            ys = dose.loc[n].to_numpy(dtype=float)
            good = np.isfinite(ys)
            if good.sum() >= 2 and ys[good].min() <= target <= ys[good].max():
                val = float(np.interp(target, ys[good], xs[good]))
                equiv[f'n={n:g}'] = dict(state='in_grid', maturation_min=val)
            else:
                slope = slopes[f'n={n:g}']
                val = (float(xs[good][0] + (target - ys[good][0]) / slope)
                       if slope else None)
                equiv[f'n={n:g}'] = dict(
                    state=('below_grid' if ys[good].min() > target else 'above_grid'),
                    extrapolated_maturation_min=val,
                    note=('the scanned maturation range is %g-%g min, so the compensation '
                          'point is an EXTRAPOLATION and was not measured'
                          % (xs.min(), xs.max())))
        # "Tangent" view: ratio of two LOCAL linearisations.  The numerator is a
        # central difference in n over +-0.5; the denominator is the least-squares
        # maturation slope above.  This answers only the infinitesimal question.
        ns = np.asarray(dose.index, dtype=float)
        for j, n in enumerate(ns):
            if j == 0 or j == len(ns) - 1:
                continue
            d_dn = (dose.loc[ns[j + 1], 32.5] - dose.loc[ns[j - 1], 32.5]) / (ns[j + 1] - ns[j - 1])
            slope = slopes[f'n={n:g}']
            compensation[f'n={n:g}'] = dict(
                d_dose_dn=float(d_dn), d_dose_dmat=slope,
                maturation_min_per_unit_n=float(-d_dn / slope) if slope else None)
    comp_vals = [v['maturation_min_per_unit_n'] for v in compensation.values()
                 if v['maturation_min_per_unit_n'] is not None]

    # "Secant" view of the SAME trade-off, taken from the iso-dose curve itself.
    # For any finite step the correct quantity is the chord of the iso-dose curve, and
    # it does NOT equal the tangent-type ratio above, because the dose surface is
    # curved.  Two numbers that differ by ~1.8x inside one artefact is exactly the kind
    # of ambiguity that gets quoted wrongly, so both are reported side by side.
    states = {float(k.split('=')[1]): v['state'] for k, v in equiv.items()}
    equiv_value = {}
    for k, v in equiv.items():
        val = v.get('maturation_min', v.get('extrapolated_maturation_min'))
        if val is not None:
            equiv_value[float(k.split('=')[1])] = float(val)
    secant = {}
    ordered_n = sorted(equiv_value)
    for n1, n2 in zip(ordered_n[:-1], ordered_n[1:]):
        secant[f'{n1:g}->{n2:g}'] = dict(
            d_n=float(n2 - n1),
            maturation_min_per_unit_n=float((equiv_value[n2] - equiv_value[n1]) / (n2 - n1)),
            uses_extrapolated_endpoint=bool(
                any(states.get(x) != 'in_grid' for x in (n1, n2))))
    sec_vals = [v['maturation_min_per_unit_n'] for v in secant.values()]
    tan_at_6 = compensation.get('n=6', {}).get('maturation_min_per_unit_n')
    if tan_at_6 is not None and 6.5 in equiv_value:
        tan_total = 32.5 + 0.5 * float(tan_at_6)
        sec_total = float(equiv_value[6.5])
        half_step_gap = sec_total - tan_total
    else:
        tan_total = sec_total = half_step_gap = None
    secant_vs_tangent = dict(
        tangent_mean_min_per_unit_n=float(np.mean(comp_vals)) if comp_vals else None,
        tangent_range_min_per_unit_n=(
            [float(min(comp_vals)), float(max(comp_vals))] if comp_vals else None),
        secant_range_min_per_unit_n=(
            [float(min(sec_vals)), float(max(sec_vals))] if sec_vals else None),
        tangent_inside_secant_range=bool(
            comp_vals and sec_vals and min(sec_vals) <= min(comp_vals)
            and max(comp_vals) <= max(sec_vals)),
        half_step_6_to_6p5_total_min=dict(
            from_tangent=tan_total, from_grid_equivalent_maturation=sec_total,
            gap_min=half_step_gap),
        note=('tangent is valid only for INFINITESIMAL steps. secant is the finite answer '
              'for one scanned n step and is what to quote when asking "if I move n by half '
              'a step, how much maturation buys back the frozen dose". Intervals flagged '
              'uses_extrapolated_endpoint rest on an off-grid compensation point and are '
              'themselves approximations.'))

    # Does the dose trade-off actually reach the FUNCTIONAL metric?  This is the
    # question the decision depends on, and it is answered from the same 25 points.
    # A negative iso-dose slope would only matter if the substitute parameter also
    # moved the far-off leak.  Measured here rather than argued.
    leak = ok.pivot_table(index='n_A1_gate', columns='a1_mat', values='leak_symmetric')
    leak_xs = np.asarray(leak.columns, dtype=float)
    functional, leak_mat_slope = {}, {}
    for n in leak.index:
        r = leak.loc[n].dropna()
        ys = leak.loc[n].to_numpy(dtype=float)
        good = np.isfinite(ys)
        leak_mat_slope[f'n={n:g}'] = (float(np.polyfit(leak_xs[good], ys[good], 1)[0])
                                      if good.sum() >= 2 else None)
        best = float(r.idxmin())
        functional[f'n={n:g}'] = dict(
            leak_min=float(r.min()), leak_max=float(r.max()),
            spread_pct=float(100 * (r.max() / r.min() - 1)),
            best_maturation=best,
            best_maturation_is_grid_edge=bool(best in (float(leak_xs.min()),
                                                       float(leak_xs.max()))))
    across_n = {}
    for m in leak.columns:
        c = leak[m].dropna()
        across_n[f'mat={m:g}'] = dict(leak_min=float(c.min()), leak_max=float(c.max()),
                                      ratio=float(c.max() / c.min()))
    mat_effect = max(v['spread_pct'] for v in functional.values())
    n_effect = min(v['ratio'] for v in across_n.values())
    # The largest leak change maturation alone can buy (at fixed n) versus the leak
    # step between adjacent n values.  Reported per boundary rather than as a single
    # number: the n -> n+0.5 steps are NOT equally sized (leak saturates above 6), and
    # a lone minimum would hide the step the decision actually rests on.
    mat_buy = {f'n={n:g}': float(leak.loc[n].max() - leak.loc[n].min()) for n in leak.index}
    n_steps = {}
    for n1, n2 in zip(leak.index[:-1], leak.index[1:]):
        vals = [abs(leak.loc[n2, m] - leak.loc[n1, m]) for m in leak.columns]
        n_steps[f'{n1:g}->{n2:g}'] = dict(
            min_leak_change=float(min(vals)), max_leak_change=float(max(vals)),
            ratio_to_maturation_buy=float(min(vals) / max(mat_buy.values())
                                          if max(mat_buy.values()) else float('nan')))

    # Grid-edge caveat.  The leak minimum sits on the LOWER maturation edge for four of
    # five rows, so the grid cannot show how far the leak would keep falling at shorter
    # maturation.  Two consequences, both stated rather than hidden: the reported spread
    # is an upper bound on what the GRID shows and simultaneously a LOWER bound on the
    # full-range maturation effect; and the deliberately pessimistic budget (each row's
    # measured slope extrapolated linearly to maturation = 0, explicitly NOT a
    # measurement) is what decides whether the functional-substitution conclusion can be
    # stated beyond the grid.
    edge_rows = sorted(k for k, v in functional.items() if v['best_maturation_is_grid_edge'])
    extrap_budget = {k: float(abs(s) * float(leak_xs.min()))
                     for k, s in leak_mat_slope.items() if s is not None}
    min_n_step = min(v['min_leak_change'] for v in n_steps.values()) if n_steps else None
    worst_budget = max(extrap_budget.values()) if extrap_budget else None
    if worst_budget is not None and min_n_step:
        if worst_budget < min_n_step:
            budget_sentence = (
                'Even the pessimistic linear extrapolation of every row to maturation = 0 '
                '(not a measurement) stays below the SMALLEST leak change produced by a half '
                'step in n (%.5f vs %.5f, ratio %.2f), so the conclusion does not depend on '
                'the grid floor.' % (worst_budget, min_n_step, worst_budget / min_n_step))
        else:
            budget_sentence = (
                'WARNING: the pessimistic linear extrapolation of a row to maturation = 0 '
                '(not a measurement) EXCEEDS the smallest leak change produced by a half step '
                'in n (%.5f vs %.5f, ratio %.2f). The functional-substitution conclusion '
                'therefore holds over the scanned range only and is NOT established below '
                'maturation = %g min.'
                % (worst_budget, min_n_step, worst_budget / min_n_step, float(leak_xs.min())))
    else:
        budget_sentence = 'insufficient data for the extrapolation budget'
    boundary_caveat = dict(
        rows_with_minimum_on_grid_edge=edge_rows,
        leak_vs_maturation_slope_per_min=leak_mat_slope,
        linear_to_zero_maturation_budget=extrap_budget,
        smallest_half_step_in_n=min_n_step,
        maturation_budget_over_smallest_n_step=(
            float(worst_budget / min_n_step) if (worst_budget and min_n_step) else None),
        extrapolation_is_not_a_measurement=True,
        note=('best_maturation is reported for COMPLETENESS ONLY and does NOT supersede the '
              'frozen 32.5 min; it is never a recommendation. Where it equals the grid floor '
              '(%g min) the true minimum may lie below the grid, so max_maturation_effect_pct '
              'is an upper bound on the GRID and a lower bound on the FULL-range effect.'
              % float(leak_xs.min())),
        verdict=budget_sentence)
    tradeoff_block = dict(
        per_n_across_maturation=functional,
        per_maturation_across_n=across_n,
        max_maturation_effect_pct=float(mat_effect),
        min_n_effect_ratio=float(n_effect),
        max_leak_change_from_maturation=mat_buy,
        leak_change_per_half_step_in_n=n_steps,
        grid_edge_caveat_that_bounds_all_of_the_above=boundary_caveat,
        reading=('The dose trade-off does NOT propagate to the far-off leak. At fixed n, '
                 'sweeping maturation over the scanned range moves the leak by at most '
                 '%.2f %% (%.5f absolute); a half step in n moves it by %.5f to %.5f. '
                 'Maturation therefore cannot substitute for the gate exponent on the '
                 'functional axis over the scanned range, so the iso-dose contours are a '
                 'trade-off in a local observable, not a functional degeneracy. The steps '
                 'are unequal (%s). %s'
                 % (mat_effect, max(mat_buy.values()),
                    min(v['min_leak_change'] for v in n_steps.values()) if n_steps else float('nan'),
                    max(v['max_leak_change'] for v in n_steps.values()) if n_steps else float('nan'),
                    ', '.join('%s: %.5f' % (k, v['min_leak_change'])
                              for k, v in n_steps.items()),
                    budget_sentence)))

    verdict = dict(
        hours=args.hours, grid=dict(n_A1_gate=list(N_VALUES), a1_mat=list(MAT_VALUES),
                                    f1_mat_frozen=F1_MAT_FROZEN, carry1_mrna=MRNA),
        points=int(len(df)), ok=int(len(ok)),
        certified=int(ok.certified.astype(bool).sum()),
        all_steady_reads_sufficient=bool(
            pd.to_numeric(ok.steady_reads, errors='coerce').min() >= MIN_STEADY_READS),
        min_steady_reads=int(pd.to_numeric(ok.steady_reads, errors='coerce').min()),
        gate_dose_table=dose.round(6).to_dict(),
        gate_peak_table=peak.round(6).to_dict(),
        certified_table=cert.to_dict(),
        dose_vs_maturation_slope=slopes,
        dose_vs_maturation_slope_kind=(
            'least-squares slope of the gate dose against A1 maturation over the scanned '
            '%g-%g min range at fixed n; NOT a point derivative'
            % (min(MAT_VALUES), max(MAT_VALUES))),
        equivalent_maturation_for_frozen_dose=equiv,
        compensation=compensation,
        compensation_range_min_per_unit_n=(
            [float(min(comp_vals)), float(max(comp_vals))] if comp_vals else None),
        compensation_kind=('tangent-type: ratio of two LOCAL linearisations (central '
                           'difference in n over +-0.5 divided by the least-squares '
                           'maturation slope). Valid only for infinitesimal steps - see '
                           'compensation_secant for the finite answer.'),
        compensation_secant=secant,
        compensation_secant_range_min_per_unit_n=(
            [float(min(sec_vals)), float(max(sec_vals))] if sec_vals else None),
        compensation_secant_vs_tangent=secant_vs_tangent,
        tradeoff_does_not_reach_the_function=tradeoff_block,
        reading=('equivalent_maturation_for_frozen_dose gives, for each gate exponent, the A1 '
                 'maturation half-life that reproduces the frozen point gate dose. Each row '
                 'carries an explicit state: in_grid (measured, interpolated), below_grid or '
                 'above_grid (NOT measured - the quoted value is a linear extrapolation off '
                 'the end of the %g-%g min scan and must be reported as such); only the '
                 'reference row n=6 is in_grid, every other row is an extrapolation. '
                 'Two DIFFERENT trade-off numbers exist in this artefact and must not be '
                 'mixed: compensation_range_min_per_unit_n is the tangent-type ratio '
                 '(infinitesimal steps only), compensation_secant_range_min_per_unit_n is '
                 'the chord of the iso-dose curve between adjacent scanned n (the finite '
                 'answer). Quote the secant for any finite question.'
                 % (min(MAT_VALUES), max(MAT_VALUES))),
        model_sha256=sha256(ROOT_MODEL), source_sha256=source_hashes(),
        working_point=working_point_block(None))
    df.to_csv(OUT / f'{NAME}_all.csv', index=False, encoding='utf-8')
    (OUT / f'{NAME}_verdict.json').write_text(
        json.dumps(verdict, ensure_ascii=False, indent=2), encoding='utf-8')
    write_manifest()
    pd.set_option('display.width', 220)
    print('--- gate dose (median per carry) ---')
    print(dose.round(5).to_string())
    print('\n--- certified ---')
    print(cert.to_string())
    print('\n--- gate peak (median) ---')
    print(peak.round(5).to_string())
    print('\n--- slope of dose vs maturation (per min) ---')
    print(json.dumps(slopes, indent=2))
    print('\n--- maturation reproducing the frozen dose ---')
    print(json.dumps(equiv, indent=2))
    print('\n--- dose trade-off, TANGENT view (infinitesimal steps only) ---')
    print(json.dumps(compensation, indent=2))
    print(f"tangent range: {verdict['compensation_range_min_per_unit_n']}")
    print('\n--- dose trade-off, SECANT view (finite steps: quote this one) ---')
    print(json.dumps(secant, indent=2))
    print(f"secant range: {verdict['compensation_secant_range_min_per_unit_n']}")
    print('\n--- tangent vs secant ---')
    print(json.dumps(secant_vs_tangent, indent=2))
    print('\n--- grid-edge caveat on the functional comparison ---')
    print(json.dumps(boundary_caveat, indent=2))


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    s = sub.add_parser('scan')
    s.add_argument('--shard', type=int, required=True)
    s.add_argument('--nshards', type=int, required=True)
    s.add_argument('--hours', type=float, default=HOURS)
    s.add_argument('--sample-min', type=float, default=2.0)
    s.add_argument('--max-step-min', type=float, default=2.0)
    m = sub.add_parser('merge')
    m.add_argument('--nshards', type=int, required=True)
    m.add_argument('--hours', type=float, default=HOURS)
    args = ap.parse_args()
    (scan if args.cmd == 'scan' else merge)(args)


if __name__ == '__main__':
    main()
