"""Strict-tolerance stability of the selected three-bit candidate.

Two axes are separated on purpose, because they answer different questions:

  tolerance axis (5 representative points)
      rtol / atol / max_step are tightened together while the OUTPUT grid stays
      at 2 min.  Keeping the grid fixed means the read windows are cut from the
      same sample instants, so the codes, event counts and margins are
      comparable run to run and only the integrator precision moves.

  output-grid axis (the candidate point only)
      sample_min and max_step are halved at fixed tight tolerance.  A finer
      output grid changes where the flux trough is detected, so this is NOT a
      precision check - it tests whether the readout depends on the 2 min grid
      the whole project has used so far.  It is reported separately.

Representative points (as agreed)
    n=6 mRNA=2 mat=32.5   the candidate centre
    n=5 mRNA=2 mat=32.5   the discrete pass boundary
    n=4 mRNA=2 mat=32.5   a certified failure
    n=7 mRNA=2 mat=32.5   the plateau
    n=6 mRNA=2 mat=60     a maturation boundary

Reported per level, and compared against that point's baseline
    complete 56-read code string, the bit0/bit1 prefix string and the bit2 string
    event counts (crossings / gates / reverse events) and the R/F pattern
    commitment and 0.30 / 0.70 band-edge margins
    setup / hold for all three bits
    paired L_rev / L_fwd / L_symmetric at the 1 %, 5 % and 10 % tail cuts
    the frozen 34-state projection (readout-level identity, plus raw trajectory gap)

The retired single-window leak_ratio is deliberately NOT carried into these
tables: `threebit51_leak_correction.json` retires it, and mixing it with the
paired metric would repeat the mistake it corrects.

    python scan_strict_tolerance.py scan --mode tol  --shard P --nshards 5
    python scan_strict_tolerance.py scan --mode grid --shard 0 --nshards 1
    python scan_strict_tolerance.py merge --mode tol  --nshards 5
    python scan_strict_tolerance.py merge --mode grid --nshards 1
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plausibility_common import (OUT, flux_array, frozen_carry0,  # noqa: E402
                                 frozen_extension, merge_shards, shard_slice,
                                 sha256, source_hashes, write_manifest,
                                 write_shard)
from model_threebit51 import ThreeBit51Model, ThreeBitCarryParameters  # noqa: E402
from model_twobit34 import CarryExpressionParameters  # noqa: E402
from check_read_commitment import characterise, window_metrics  # noqa: E402
from scan_carry_pairing import analyse as paired_analyse  # noqa: E402

ROOT_MODEL = Path(__file__).resolve().parents[1] / 'model_threebit51.py'

POINTS = (
    dict(label='candidate_centre',  n_A1_gate=6.0, mrna=2.0, mat=32.5),
    dict(label='pass_boundary',     n_A1_gate=5.0, mrna=2.0, mat=32.5),
    dict(label='certified_failure', n_A1_gate=4.0, mrna=2.0, mat=32.5),
    dict(label='plateau',           n_A1_gate=7.0, mrna=2.0, mat=32.5),
    dict(label='maturation_edge',   n_A1_gate=6.0, mrna=2.0, mat=60.0),
)

LEVELS = (
    dict(level='baseline', rtol=2e-7, atol=2e-9, max_step_min=2.0, sample_min=2.0),
    dict(level='tight', rtol=1e-9, atol=1e-11, max_step_min=1.0, sample_min=2.0),
    dict(level='ultra', rtol=1e-11, atol=1e-13, max_step_min=0.5, sample_min=2.0),
)

GRID_LEVELS = (
    dict(level='grid_1min', rtol=1e-9, atol=1e-11, max_step_min=0.5, sample_min=1.0),
    dict(level='grid_0p5min', rtol=1e-9, atol=1e-11, max_step_min=0.25, sample_min=0.5),
)

# fields taken from characterise(); its leak columns are intentionally NOT taken
CHAR_FIELDS = (
    'certified', 'crossings', 'gates', 'reverse_events',
    'cold_reads', 'steady_reads', 'cold_valid', 'steady_valid',
    'cold_sequence', 'steady_sequence', 'steady_min_commitment_all_bits',
    'unlabelled_reads_all', 'clock_period_h', 'read_window_h_median',
    'read_windows_clipped',
    'causal_exactly_one_gate_and_flip', 'causal_order', 'causal_dirs_alternate',
    'causal_unassigned_gates', 'causal_unassigned_crossings',
    'gate_contrast_passed', 'gate_contrast_off_on_peak_ratio',
    'bit0_setup_min_h', 'bit0_hold_min_h', 'bit1_setup_min_h', 'bit1_hold_min_h',
    'bit2_setup_min_h', 'bit2_hold_min_h', 'setup_min_h_global', 'hold_min_h_global',
    'drop8_bit0_commitment_min', 'drop8_bit1_commitment_min', 'drop8_bit2_commitment_min',
    'drop8_bit0_margin_to_band_low', 'drop8_bit1_margin_to_band_low',
    'drop8_bit2_margin_to_band_low',
    'drop8_bit0_margin_to_band_high', 'drop8_bit1_margin_to_band_high',
    'drop8_bit2_margin_to_band_high',
)


def build(p):
    carry1 = CarryExpressionParameters(mrna_half_life_min=p['mrna'],
                                       activator_maturation_half_life_min=p['mat'],
                                       repressor_maturation_half_life_min=p['mat'])
    return ThreeBit51Model(extension=frozen_extension(),
                           carry=ThreeBitCarryParameters(carry0=frozen_carry0(),
                                                         carry1=carry1),
                           n_A1_gate=p['n_A1_gate'])


def readout_strings(model, sol):
    from verify_threebit51 import read_windows_three
    flux = flux_array(model, sol)
    reads = read_windows_three(sol.t, flux, (sol.y[16], sol.y[27], sol.y[44]))
    return dict(
        code=''.join('x' if r['value'] is None else str(r['value']) for r in reads),
        prefix_code=''.join(r['bits'][0]['label'] + r['bits'][1]['label'] for r in reads),
        bit2_code=''.join(r['bits'][2]['label'] for r in reads))


def run_level(p, lvl, hours):
    model = build(p)
    sol = model.simulate(hours=hours, sample_min=lvl['sample_min'],
                         rtol=lvl['rtol'], atol=lvl['atol'],
                         max_step_min=lvl['max_step_min'])
    row = dict(label=p['label'], n_A1_gate=p['n_A1_gate'], mrna=p['mrna'], mat=p['mat'],
               hours=hours, **lvl)
    ch = characterise(model, sol, hours, p=None, crosscheck=False)
    for k in CHAR_FIELDS:
        row[k] = ch.get(k)
    row.update(readout_strings(model, sol))
    agg, _, pairs, dropped, anomalies = paired_analyse(model, sol, 0.05)
    row.update(n_pairs=agg['n_pairs'], n_R=agg['n_R'], n_F=agg['n_F'],
               n_dropped=agg['n_dropped'], pattern=agg['pattern'],
               same_type_adjacent=sum(1 for a in anomalies
                                      if a['reason'] == 'same type adjacent'))
    for f in (0.01, 0.05, 0.10):
        a, _, _, _, _ = paired_analyse(model, sol, f)
        pc = int(round(f * 100))
        for key in ('L_rev', 'L_fwd', 'L_symmetric'):
            row[f'{pc}pct_{key}'] = a[key]['median']
            row[f'{pc}pct_{key}_min'] = a[key]['min']
            row[f'{pc}pct_{key}_max'] = a[key]['max']
    return row, sol


def _aligned_gap(sol_a, sol_b):
    """Trajectory gap on a common time base.

    Different levels use different OUTPUT grids (2 / 1 / 0.5 min), so the state
    arrays have different lengths and cannot be subtracted directly.  The grids
    are nested, but NOT by an integer ratio of array lengths (18001 * 4 = 72004
    != 72001 because of the inclusive endpoint), so the alignment is done by TIME
    VALUE: sample the finer solution at the coarser solution's instants and
    require an exact match.  Returns None when the instants do not coincide.
    """
    coarse, fine = (sol_a, sol_b) if sol_a.t.size <= sol_b.t.size else (sol_b, sol_a)
    idx = np.clip(np.searchsorted(fine.t, coarse.t), 0, fine.t.size - 1)
    if not np.allclose(fine.t[idx], coarse.t, rtol=0.0, atol=1e-9):
        return None, None, False
    fa, ca = fine.y[:, idx], coarse.y
    if fa.shape != ca.shape:
        return None, None, False
    return (float(np.max(np.abs(fa[:34] - ca[:34]))),
            float(np.max(np.abs(fa - ca))), True)


def compare(row, base_row, base_sol, sol):
    """Identity-of-answer checks against the baseline run of the same point.

    Both `sol` and `base_sol` are solve_ivp result objects; the state arrays are
    reached through `.y`.  (An OdeResult is a dict subclass, so slicing the
    object itself raises KeyError(slice) rather than an obvious TypeError.)
    """
    prefix_gap, full_gap, comparable = _aligned_gap(sol, base_sol)
    out = dict(
        code_identical=bool(row['code'] == base_row['code']),
        steady_code_identical=bool(row['steady_sequence'] == base_row['steady_sequence']),
        prefix_code_identical=bool(row['prefix_code'] == base_row['prefix_code']),
        bit2_code_identical=bool(row['bit2_code'] == base_row['bit2_code']),
        events_identical=bool(row['crossings'] == base_row['crossings'] and
                              row['gates'] == base_row['gates'] and
                              row['reverse_events'] == base_row['reverse_events']),
        pattern_identical=bool(row['pattern'] == base_row['pattern']),
        certified_identical=bool(row['certified'] == base_row['certified']),
        traj_gap_comparable=bool(comparable),
        prefix_traj_gap=prefix_gap,
        full_traj_gap=full_gap)
    for key in ('drop8_bit2_commitment_min', 'drop8_bit2_margin_to_band_low',
                'drop8_bit2_margin_to_band_high', 'setup_min_h_global',
                'hold_min_h_global', '5pct_L_symmetric', '5pct_L_rev', '5pct_L_fwd'):
        try:
            out[f'{key}_delta'] = float(row[key] - base_row[key])
        except (TypeError, KeyError):
            out[f'{key}_delta'] = None
    return out


def scan(args):
    if args.mode == 'tol':
        pts = shard_slice(list(POINTS), args.shard, args.nshards)
        levels = list(LEVELS)
    else:
        pts = [POINTS[0]]
        levels = list(GRID_LEVELS)
    if getattr(args, 'levels', ''):
        want = {x.strip() for x in args.levels.split(',') if x.strip()}
        levels = [l for l in levels if l['level'] in want]
        if not levels:
            raise SystemExit(f'no level matched {sorted(want)}')
        print(f'restricted to levels: {[l["level"] for l in levels]}', flush=True)
    tag = f'shard{args.shard:02d}of{args.nshards:02d}'
    name = f'strict_tolerance_{args.mode}'
    rows = []
    for p in pts:
        base_row, base_sol = None, None
        for lvl in levels:
            started = time.perf_counter()
            try:
                row, sol = run_level(p, lvl, args.hours)
                if base_row is None:
                    row.update(compare(row, row, sol, sol))
                    base_row, base_sol = row, sol
                else:
                    # A comparison failure must NOT discard the run's own metrics.
                    try:
                        row.update(compare(row, base_row, base_sol, sol))
                        row['compare_error'] = ''
                    except Exception as cerr:                    # noqa: BLE001
                        row['compare_error'] = repr(cerr)
                    row['same_as_baseline'] = bool(
                        row['code_identical'] and row['events_identical'] and
                        row['pattern_identical'] and row['certified_identical'])
                    del sol
                row.update(ok=True, error='')
            except Exception as exc:                             # noqa: BLE001
                row = dict(label=p['label'], level=lvl['level'], ok=False,
                           error=repr(exc), certified=None)
            row['runtime_s'] = round(time.perf_counter() - started, 2)
            rows.append(row)
            print(f"[{tag}] {row['label']:<18} {lvl['level']:<12} "
                  f"cert={row.get('certified')} code_eq={row.get('code_identical')} "
                  f"pairs={row.get('n_pairs')} pattern={row.get('pattern')} "
                  f"same_adj={row.get('same_type_adjacent')} "
                  f"L_sym5={row.get('5pct_L_symmetric')} "
                  f"prefix_gap={row.get('prefix_traj_gap')} ({row['runtime_s']}s)", flush=True)
    write_shard(name, tag, rows, dict(mode=args.mode, shard=args.shard,
                                      nshards=args.nshards, hours=args.hours,
                                      points=len(pts), model_sha256=sha256(ROOT_MODEL),
                                      source_hashes=source_hashes()))
    print(f'[{tag}] wrote {len(rows)} rows')


def merge(args):
    name = f'strict_tolerance_{args.mode}'
    df = merge_shards(name, args.nshards)
    ok = df[df.ok.astype(bool)]
    verdict = dict(mode=args.mode, rows=len(df), ok=int(len(ok)),
                   hours=args.hours if 'hours' in df.columns else None,
                   levels=[l['level'] for l in (LEVELS if args.mode == 'tol' else GRID_LEVELS)],
                   model_sha256=sha256(ROOT_MODEL), source_sha256=source_hashes(),
                   per_point={})
    for label, sub in ok.groupby('label'):
        sub = sub.sort_values('runtime_s')
        base = sub[sub.level == 'baseline']
        base = base.iloc[0] if len(base) else sub.iloc[0]
        entry = dict(levels={}, selection_verdict={})
        for _, r in sub.iterrows():
            entry['levels'][r.level] = dict(
                rtol=r.rtol, atol=r.atol, max_step_min=r.max_step_min,
                sample_min=r.sample_min, runtime_s=r.runtime_s,
                certified=_b(r.certified), crossings=_i(r.crossings),
                gates=_i(r.gates), reverse_events=_i(r.reverse_events),
                n_pairs=_i(r.n_pairs), n_R=_i(r.n_R), n_F=_i(r.n_F),
                same_type_adjacent=_i(r.same_type_adjacent), pattern=r.get('pattern'),
                code=r.get('code'), steady_sequence=r.get('steady_sequence'),
                code_identical=_b(r.get('code_identical')),
                prefix_code_identical=_b(r.get('prefix_code_identical')),
                bit2_code_identical=_b(r.get('bit2_code_identical')),
                events_identical=_b(r.get('events_identical')),
                pattern_identical=_b(r.get('pattern_identical')),
                same_as_baseline=_b(r.get('same_as_baseline')),
                prefix_traj_gap=_f(r.get('prefix_traj_gap')),
                full_traj_gap=_f(r.get('full_traj_gap')),
                commitment_bit2=_f(r.get('drop8_bit2_commitment_min')),
                margin_low=_f(r.get('drop8_bit2_margin_to_band_low')),
                margin_high=_f(r.get('drop8_bit2_margin_to_band_high')),
                setup_global=_f(r.get('setup_min_h_global')),
                hold_global=_f(r.get('hold_min_h_global')),
                L_rev_1pct=_f(r.get('1pct_L_rev')), L_rev_5pct=_f(r.get('5pct_L_rev')),
                L_rev_10pct=_f(r.get('10pct_L_rev')),
                L_fwd_5pct=_f(r.get('5pct_L_fwd')),
                L_symmetric_1pct=_f(r.get('1pct_L_symmetric')),
                L_symmetric_5pct=_f(r.get('5pct_L_symmetric')),
                L_symmetric_10pct=_f(r.get('10pct_L_symmetric')))
        non_base = sub[sub.level != 'baseline']
        entry['selection_verdict'] = dict(
            all_code_identical=bool(non_base.code_identical.all()) if len(non_base) else None,
            all_prefix_identical=bool(non_base.prefix_code_identical.all())
            if len(non_base) else None,
            all_bit2_identical=bool(non_base.bit2_code_identical.all())
            if len(non_base) else None,
            all_events_identical=bool(non_base.events_identical.all())
            if len(non_base) else None,
            all_pattern_identical=bool(non_base.pattern_identical.all())
            if len(non_base) else None,
            certified_baseline=_b(base.certified),
            certified_under_tightening=bool(non_base.certified.apply(_b).all())
            if len(non_base) else None,
            max_prefix_traj_gap=_f(non_base.prefix_traj_gap.max()) if len(non_base) else None,
            max_full_traj_gap=_f(non_base.full_traj_gap.max()) if len(non_base) else None,
            L_symmetric_5pct_spread=(_f(non_base['5pct_L_symmetric'].max() -
                                        non_base['5pct_L_symmetric'].min())
                                     if len(non_base) else None))
        verdict['per_point'][label] = entry
    verdict['all_points_stable'] = bool(all(
        (v['selection_verdict']['all_code_identical'] is not False) and
        (v['selection_verdict']['all_events_identical'] is not False) and
        (v['selection_verdict']['all_prefix_identical'] is not False)
        for v in verdict['per_point'].values()))
    df.to_csv(OUT / f'{name}_all.csv', index=False, encoding='utf-8')
    (OUT / f'{name}_verdict.json').write_text(
        json.dumps(verdict, ensure_ascii=False, indent=2), encoding='utf-8')
    write_manifest()
    print(json.dumps(verdict, ensure_ascii=False, indent=2))


def _f(v):
    return None if v is None or (not isinstance(v, str) and pd.isna(v)) else float(v)


def _b(v):
    if v is None:
        return None
    if isinstance(v, str):
        return v
    try:
        if pd.isna(v):
            return None
    except (TypeError, ValueError):
        pass
    if isinstance(v, (bool, np.bool_)):
        return bool(v)
    return bool(v)


def _i(v):
    if v is None or (not isinstance(v, str) and pd.isna(v)):
        return None
    try:
        return int(v)
    except (TypeError, ValueError):
        return None


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='mode_cmd', required=True)
    s = sub.add_parser('scan')
    s.add_argument('--mode', choices=('tol', 'grid'), required=True)
    s.add_argument('--shard', type=int, required=True)
    s.add_argument('--nshards', type=int, required=True)
    s.add_argument('--hours', type=float, default=600.0)
    s.add_argument('--levels', type=str, default='',
                   help='comma-separated level names to run (default: all)')
    m = sub.add_parser('merge')
    m.add_argument('--mode', choices=('tol', 'grid'), required=True)
    m.add_argument('--nshards', type=int, required=True)
    m.add_argument('--hours', type=float, default=600.0)
    args = ap.parse_args()
    (scan if args.mode_cmd == 'scan' else merge)(args)


if __name__ == '__main__':
    main()
