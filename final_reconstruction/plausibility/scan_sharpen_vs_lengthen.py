"""Check 5 - the main scan: "lengthen the carry pulse" versus "sharpen the gate".

Both arms are evaluated with ONE criterion set, so the comparison cannot drift.

  Arm A (lengthen): carry1 mRNA half-life x A1/F1 common maturation half-life.
      This is the route proposed after the first 35-point scan: the optimum sat
      on the maturation grid edge and crossings rose monotonically with it.
      Only already-existing slots are moved.

  Arm B (sharpen): the exponent of the A1 arm inside the carry-2 gate,
      n_A1_gate in {4,5,6,7,8}, at the frozen maturation (32.5 min) and at
      60 min.  This adds ONE new, explicitly labelled interface parameter and
      does not touch Zeng's table in the production model.
      Mechanistically it is a different lever: lengthening adds dose, sharpening
      removes the off-state floor.  Because the complex Hill is quadratic at low
      drive, removing the floor is far more effective than adding dose.

IMPORTANT - what the run-time patch does
    The diagnostic implementation patches model.ZENG['n_A'][1], which is used
    both by the gate's A1 arm and by the Hill that drives F1.  At the operating
    point the second use is inert (F1's steady state stays far below K_F1 either
    way), so the observed effect is attributable to the gate arm, but a
    production implementation must expose the gate exponent as its own new
    parameter and re-run this comparison before the claim is final.

Pre-registered decision rules (fixed before any run)
    R1 SUCCESS      : a point has certified=True, crossings == gates == 14 and
                      leak_worst_Jrev2 <= LEAK_BUDGET.
    R2 SATURATION   : within arm A at fixed mRNA, crossings(mat=120) <=
                      crossings(mat=60)  ->  more maturation no longer helps.
    R3 SHARPEN WINS : arm B certifies at the frozen 32.5 min while arm A has no
                      certified point at any maturation  ->  the lever is
                      sharpness, not length.
    R4 BOTH FAIL    : no certified point in either arm  ->  neither lever works
                      and the interface needs a structural change.
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
from plausibility_common import (FROZEN_N_A1, OUT, ZENG_N_A1_TABLE,  # noqa: E402
                                 build_decoupled,
                                 build_threebit, frozen_carry0, gate_windows,
                                 leak_budget, load_valid_decoupling, merge_shards,
                                 restore, shard_slice, source_hashes,
                                 verify_decoupling, write_manifest, write_shard)
from model_twobit34 import CarryExpressionParameters  # noqa: E402

NAME = 'sharpen_vs_lengthen'
LEAK_BUDGET = 0.10          # pre-registered: worst inter-carry reverse dose
CROSSINGS_REQUIRED = 14

GRID_A = [(m, t) for m in (2.0, 3.0, 4.0, 5.0) for t in (60.0, 70.0, 80.0, 90.0, 100.0, 120.0)]
GRID_B = [(n, m, t) for n in (4.0, 5.0, 6.0, 7.0, 8.0)
          for m in (2.0, 4.0) for t in (32.5, 60.0)]


def points(arm):
    if arm == 'A':
        return [dict(arm='A', mrna=m, mat=t, n_A1=None) for m, t in GRID_A]
    if arm == 'B':
        return [dict(arm='B', mrna=m, mat=t, n_A1=n) for n, m, t in GRID_B]
    return points('A') + points('B')


def point_id(p):
    if p['arm'] == 'A':
        return f"A|mRNA={p['mrna']:g}|mat={p['mat']:g}"
    return f"B|nA1={p['n_A1']:g}|mRNA={p['mrna']:g}|mat={p['mat']:g}"


def evaluate(p, hours, sample_min, max_step_min):
    from verify_threebit51 import analyse_threebit
    started = time.perf_counter()
    row = dict(point_id=point_id(p), arm=p['arm'], mrna=p['mrna'], mat=p['mat'],
               n_A1_gate=p['n_A1'], decoupled=bool(p['arm'] == 'B'),
               n_f1_drive=(FROZEN_N_A1 if p['arm'] == 'B' else None))
    carry1 = CarryExpressionParameters(mrna_half_life_min=p['mrna'],
                                       activator_maturation_half_life_min=p['mat'],
                                       repressor_maturation_half_life_min=p['mat'])
    if p['arm'] == 'B':
        # gate exponent moved on its own; the F1 production Hill stays frozen
        model, patch = build_decoupled(n_A1_gate=p['n_A1'], carry1=carry1)
    else:
        # ARCHIVE FIDELITY: when this scan was RUN, build_threebit left the gate
        # exponent unset, so arm A executed at the ZENG TABLE value (4.0) - the
        # configuration the 20-point grid shows always fails.  The helper now
        # defaults to the frozen working point (6.0), which would silently turn
        # arm A into a DIFFERENT experiment.  Pin the historical value so the
        # archived run stays reproducible.
        # This scan is SUPERSEDED (see INVALIDATED_ARTIFACTS.json); the lever
        # question is answered by scan_gate_segments.py / scan_carry_pairing.py.
        model, patch = build_threebit(carry1=carry1, n_A1_gate=ZENG_N_A1_TABLE)
    try:
        sol = model.simulate(hours=hours, sample_min=sample_min,
                             max_step_min=max_step_min)
        analysis = analyse_threebit(model, sol, hours)
        sig = model.diagnostic_signals(sol.y)
        t = sol.t
        windows = gate_windows(t, sig['g1'])
        budget = leak_budget(t, sig['g1'], sig['J_rev2'], windows)
        S2 = sig['S2']
        late = t > 0.5 * hours
        k = int(np.flatnonzero(late)[int(np.argmin(S2[late]))])
        row.update(
            ok=True, error='',
            certified=bool(analysis['certified']),
            cold_passed=bool(analysis['cold_start']['passed']),
            steady_passed=bool(analysis['steady_state']['passed']),
            cold_sequence=analysis['cold_start']['sequence'],
            steady_sequence=analysis['steady_state']['sequence'],
            reads=len(analysis['read_windows']),
            crossings=len(analysis['bit2_crossings']),
            gates=len(analysis['carry1_gate_events']),
            reverse_events=len(analysis['bit1_reverse_events']),
            one_to_one=bool(analysis['causal_verdict']['exactly_one_gate_and_flip_per_late_reverse']),
            causal_order=bool(analysis['causal_verdict']['causal_order_after_reverse_start']),
            directions_alternate=bool(analysis['causal_verdict']['bit2_directions_alternate']),
            gate_on_peak=float(analysis['gate_contrast']['on_cycle_peak_median']),
            gate_off_peak_max=float(analysis['gate_contrast']['off_cycle_peak_max']),
            gate_off_on_ratio=analysis['gate_contrast']['off_peak_to_on_peak'],
            leak_worst_Jrev2=budget['worst_Jrev2'],
            leak_median_Jrev2=budget['median_Jrev2'],
            gate_width_median_h=float(np.median([w['width_h'] for w in windows])) if windows else np.nan,
            bit2_setup_h=analysis['bit_margins']['bit2']['min_setup_h'],
            bit2_hold_h=analysis['bit_margins']['bit2']['min_hold_h'],
            S2_min_late=float(S2[k]), R2_at_S2min=float(sol.y[42][k]),
            I2_at_S2min=float(sol.y[36][k]), C2_at_S2min=float(sol.y[43][k]),
            G1_below_band=bool(S2[k] <= 0.30),
            patched_previous_n_A=repr(patch.get('n_A', None)),
        )
    except Exception as exc:                                     # noqa: BLE001
        row.update(ok=False, error=repr(exc), certified=False, steady_passed=False,
                   crossings=0, gates=0, reverse_events=0)
    finally:
        restore(patch)
    row['runtime_s'] = round(time.perf_counter() - started, 2)
    return row


def scan(args):
    if args.arm in ('B', 'AB'):
        art = load_valid_decoupling()
        if art is None:
            raise SystemExit(
                'REFUSING to scan the sharpen arm.\n'
                'The gate exponent must be decoupled from the F1 production exponent first, '
                'and a fresh, passing decoupling artefact must exist.\n'
                'Run:  python scan_sharpen_vs_lengthen.py verify-decoupling\n'
                'Then re-run this scan.')
        print(f"decoupling verified: initial-state gap={art['initial_state_max_abs_gap']:.1e}, "
              f"trajectory gap={art['trajectory_max_abs_gap']:.1e}, "
              f"frozen n_A1={art['frozen_n_A1']}", flush=True)
    pts = shard_slice(points(args.arm), args.shard, args.nshards)
    tag = f'shard{args.shard:02d}of{args.nshards:02d}'
    rows = []
    for p in pts:
        rows.append(evaluate(p, args.hours, args.sample_min, args.max_step_min))
        print(f"[{tag}] {point_id(p):<28} certified={rows[-1]['certified']} "
              f"crossings={rows[-1]['crossings']} leak={rows[-1].get('leak_worst_Jrev2')} "
              f"({rows[-1]['runtime_s']}s)", flush=True)
    write_shard(NAME, tag, rows, dict(kind='scan', arm=args.arm, shard=args.shard,
                                      nshards=args.nshards, hours=args.hours,
                                      points=len(pts),
                                      criterion=dict(leak_budget=LEAK_BUDGET,
                                                     crossings_required=CROSSINGS_REQUIRED,
                                                     source='verify_threebit51.analyse_threebit'),
                                      criticality='slow',
                                      source_hashes=source_hashes()))
    print(f'[{tag}] wrote {len(rows)} rows')


def decide(df, baseline_reverse_events=None):
    ok = df[df['ok'].astype(bool)] if 'ok' in df.columns else df
    winners = ok[(ok.certified.astype(bool)) &
                 (ok.crossings >= CROSSINGS_REQUIRED) &
                 (ok.leak_worst_Jrev2.fillna(9e9) <= LEAK_BUDGET)] if len(ok) else ok
    verdict = dict(R1_success=bool(len(winners)),
                   R1_points=winners.point_id.tolist() if len(winners) else [])
    if len(ok):
        a = ok[ok.arm == 'A']
        sat = []
        for m, sub in a.groupby('mrna'):
            sub = sub.sort_values('mat')
            if len(sub) >= 2:
                first = int(sub.crossings.iloc[0]); last = int(sub.crossings.iloc[-1])
                sat.append(dict(mrna=float(m), crossings_at_first=first, crossings_at_last=last,
                                saturated=bool(last <= first),
                                overshoots_required=int(last) > CROSSINGS_REQUIRED))
        verdict['R2_arm_A_saturation'] = sat
        verdict['R2_arm_A_overshoot'] = dict(
            points_above_required=int((a.crossings > CROSSINGS_REQUIRED).sum()),
            points_at_required=int((a.crossings == CROSSINGS_REQUIRED).sum()),
            note=('arm A does not saturate below the required count; past mat~90 it '
                  'overshoots to 20-27 crossings, i.e. ringing rather than counting'))
        b = ok[(ok.arm == 'B') & (ok.mat == 32.5)]
        verdict['R3_sharpen_certifies_at_frozen_maturation'] = bool(
            len(b) and b.certified.astype(bool).any())
        verdict['R3_arm_A_has_any_certified'] = bool(
            len(a) and a.certified.astype(bool).any())
        verdict['R4_both_arms_fail'] = bool(not len(ok) or not ok.certified.astype(bool).any())
    if baseline_reverse_events is not None:
        bad = ok[ok.reverse_events != baseline_reverse_events]
        verdict['projection_guard'] = dict(
            expected_reverse_events=int(baseline_reverse_events),
            violations=bad.point_id.tolist())
    return verdict


def merge(args):
    df = merge_shards(NAME, args.nshards)
    # Projection guard: the frozen configuration at the SAME duration is already
    # one of the scanned points (arm B, gate exponent at the frozen value).
    # Using it avoids comparing a short reference run against longer scan rows.
    frozen_id = 'B|nA1=4|mRNA=2|mat=32.5'
    expected = None
    ids = set(df.point_id) if 'point_id' in df.columns else set()
    if frozen_id in ids:
        expected = int(df.loc[df.point_id == frozen_id, 'reverse_events'].iloc[0])
        print(f'projection guard reference = in-scan frozen point {frozen_id} '
              f'(reverse_events={expected}, same {args.nshards}-shard duration as the scan)')
    elif args.projection_guard:
        try:
            ref = evaluate(dict(arm='A', mrna=2.0, mat=32.5, n_A1=None),
                           args.guard_hours, args.sample_min, args.max_step_min)
            expected = int(ref['reverse_events'])
            print(f'projection guard reference = live run at {args.guard_hours} h '
                  f'(reverse_events={expected}) - NOTE this must match the scan duration')
        except Exception as exc:                                 # noqa: BLE001
            print(f'projection guard reference failed: {exc!r}')
    verdict = decide(df, expected)
    df.to_csv(OUT / f'{NAME}_all.csv', index=False, encoding='utf-8')
    (OUT / f'{NAME}_verdict.json').write_text(
        json.dumps(dict(rows=len(df), verdict=verdict), ensure_ascii=False, indent=2),
        encoding='utf-8')
    write_manifest()
    pd.set_option('display.width', 260)
    cols = [c for c in ('point_id', 'certified', 'crossings', 'gates', 'reverse_events',
                        'leak_worst_Jrev2', 'gate_off_peak_max', 'S2_min_late',
                        'bit2_hold_h', 'steady_sequence') if c in df.columns]
    print(df[cols].sort_values(['certified', 'crossings'], ascending=False).to_string(index=False))
    print()
    print(json.dumps(verdict, ensure_ascii=False, indent=2))


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='mode', required=True)
    s = sub.add_parser('scan')
    s.add_argument('--arm', choices=['A', 'B', 'AB'], default='AB')
    s.add_argument('--shard', type=int, required=True)
    s.add_argument('--nshards', type=int, required=True)
    s.add_argument('--hours', type=float, default=600.0)
    s.add_argument('--sample-min', type=float, default=2.0)
    s.add_argument('--max-step-min', type=float, default=2.0)
    m = sub.add_parser('merge')
    m.add_argument('--nshards', type=int, required=True)
    m.add_argument('--projection-guard', action='store_true',
                   help='run one frozen-configuration reference and compare bit1 event counts')
    m.add_argument('--guard-hours', type=float, default=200.0)
    m.add_argument('--sample-min', type=float, default=2.0)
    m.add_argument('--max-step-min', type=float, default=2.0)
    v = sub.add_parser('verify-decoupling',
                       help='prove the gate exponent is independent of the F1 drive exponent')
    v.add_argument('--hours', type=float, default=100.0)
    args = ap.parse_args()
    if args.mode == 'scan':
        scan(args)
    elif args.mode == 'merge':
        merge(args)
    else:
        res = verify_decoupling(hours=args.hours)
        write_manifest()
        print(json.dumps({k: v for k, v in res.items() if k != 'source_sha256'},
                         ensure_ascii=False, indent=2))
        if not res['passed']:
            raise SystemExit('decoupling verification FAILED - do not scan arm B')


if __name__ == '__main__':
    main()
