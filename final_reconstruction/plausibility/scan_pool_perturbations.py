"""Molecular-pool perturbations of the selected three-bit candidate.

Scope
-----
The 34-state base already has its own robustness scan, so only the pools that are
NEW in the three-bit extension, or that directly drive it, are perturbed:

    S2          [44]                 bit2 DNA
    C2          [43]                 bit2 explicit complex
    RDF2_full   [40, 41, 42]        b2_M_R, b2_R_u, b2_R
    Int2_full   [34, 35, 36]        b2_M_I, b2_I_u, b2_I
    A1_full     [45, 47, 48]        A1 plus its two expression precursors
    F1_full     [46, 49, 50]        F1 plus its two expression precursors

Every perturbation is applied to all EIGHT orbit-derived initial states, so the
question answered is not "does one state survive" but "is the attractor basin
wide enough that all eight phase-consistent states survive".

Perturbation set
    medium  0.8x / 1.2x on every pool x 8 initial states       (96 runs)
    s2add   S2 +-0.1, +-0.2 and a full flip x 8 initial states (40 runs)

S2 is a DNA occupancy fraction, so an additive shift is clipped into [0, 1] and
every clip is recorded; an unphysical S > 1 is not a state of this model.

Outcome classification per run
    lost_lock            not certified (mod-8 broken, <16 valid reads, or the
                         causal / gate-contrast conditions fail)
    recovered_original   certified AND the steady code string equals the
                         unperturbed one for that initial state
    legal_mod8_shift     certified but the code is offset - the perturbation
                         advanced or delayed the counter, which is still a legal
                         mod-8 count
Separately flagged, because they are different failures:
    readout_ok_no_causal strict mod-8 readout but the causal verdict fails
    alternation_broken   same-type carry adjacency (the failure fingerprint)

The paired L_rev / L_fwd / L_symmetric at the 1 %, 5 % and 10 % cuts are recorded
so that phase dependence of the leak can be checked across the eight states.
The retired single-window leak_ratio is not computed at all.

    python scan_pool_perturbations.py scan  --phase medium|s2add|all --shard S --nshards N
    python scan_pool_perturbations.py merge --phase medium|s2add|all --nshards N
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
from plausibility_common import (OUT, frozen_carry0, frozen_extension,  # noqa: E402
                                 merge_shards, shard_slice, sha256,
                                 source_hashes, write_manifest, write_shard)
from scan_strict_tolerance import build, readout_strings  # noqa: E402
from check_read_commitment import characterise  # noqa: E402
from scan_carry_pairing import analyse as paired_analyse  # noqa: E402

ROOT_MODEL = Path(__file__).resolve().parents[1] / 'model_threebit51.py'
STATES_JSON = OUT / 'eight_initial_states.json'
BASELINE_CSV = OUT / 'eight_initial_all.csv'
SELECTED = dict(label='candidate_centre', n_A1_gate=6.0, mrna=2.0, mat=32.5)
CANONICAL_PATTERN = 'FRFRFRFRFRFRFR'

POOLS = (
    dict(pool='S2', idx=(44,)),
    dict(pool='C2', idx=(43,)),
    dict(pool='RDF2_full', idx=(40, 41, 42)),
    dict(pool='Int2_full', idx=(34, 35, 36)),
    dict(pool='A1_full', idx=(45, 47, 48)),
    dict(pool='F1_full', idx=(46, 49, 50)),
)
MULT = (0.8, 1.2)
S2_ADD = (0.1, -0.1, 0.2, -0.2, None)      # None means flip


def jobs(phase):
    out = []
    if phase in ('medium', 'all'):
        for p in POOLS:
            for f in MULT:
                for v in range(8):
                    out.append(dict(pool=p['pool'], idx=list(p['idx']), mode='mult',
                                    factor=f, init=v))
    if phase in ('s2add', 'all'):
        for spec in S2_ADD:
            for v in range(8):
                out.append(dict(pool='S2', idx=[44],
                                mode='flip' if spec is None else 'add',
                                factor=None if spec is None else spec, init=v))
    return out


def load_states():
    payload = json.loads(STATES_JSON.read_text(encoding='utf-8'))
    return {int(k, 2): np.asarray(v, dtype=float)
            for k, v in payload['initial_states'].items()}


def load_baselines():
    df = pd.read_csv(BASELINE_CSV)
    return {int(r.initial_value): dict(steady_sequence=r.steady_sequence,
                                       cold_sequence=r.cold_sequence,
                                       certified=bool(r.certified))
            for _, r in df.iterrows()}


def apply_perturbation(y, job):
    y = np.asarray(y, dtype=float).copy()
    idx = list(job['idx'])
    before = y[idx].copy()
    if job['mode'] == 'mult':
        y[idx] = y[idx] * float(job['factor'])
    elif job['mode'] == 'add':
        y[idx] = y[idx] + float(job['factor'])
    else:
        y[idx] = 1.0 - y[idx]
    clipped = 0
    if job['pool'] == 'S2':
        raw = y[idx].copy()
        y[idx] = np.clip(raw, 0.0, 1.0)
        clipped = int(np.sum(raw != y[idx]))
    return y, before, y[idx].copy(), clipped


def run_one(job, states, baselines, hours, sample_min, max_step_min, rtol, atol):
    started = time.perf_counter()
    row = dict(kind=job['pool'], mode=job['mode'], factor=job['factor'],
               init=job['init'], init_bits=f"{job['init']:03b}", idx=json.dumps(job['idx']))
    try:
        y0, before, after, clipped = apply_perturbation(states[job['init']], job)
        model = build(SELECTED)
        sol = model.simulate(hours=hours, sample_min=sample_min,
                             max_step_min=max_step_min, rtol=rtol, atol=atol,
                             initial_state=y0)
        ch = characterise(model, sol, hours, p=None, crosscheck=False)
        row.update(certified=ch['certified'], crossings=ch['crossings'],
                   gates=ch['gates'], reverse_events=ch['reverse_events'],
                   cold_reads=ch['cold_reads'], steady_reads=ch['steady_reads'],
                   cold_sequence=ch['cold_sequence'],
                   steady_sequence=ch['steady_sequence'],
                   causal_exactly_one_gate_and_flip=ch['causal_exactly_one_gate_and_flip'],
                   causal_order=ch['causal_order'],
                   causal_dirs_alternate=ch['causal_dirs_alternate'],
                   gate_contrast_passed=ch['gate_contrast_passed'],
                   steady_min_commitment_all_bits=ch['steady_min_commitment_all_bits'],
                   bit2_commitment_min=ch['drop8_bit2_commitment_min'],
                   bit2_margin_to_band_low=ch['drop8_bit2_margin_to_band_low'],
                   bit2_margin_to_band_high=ch['drop8_bit2_margin_to_band_high'],
                   setup_min_h_global=ch['setup_min_h_global'],
                   hold_min_h_global=ch['hold_min_h_global'])
        row.update(readout_strings(model, sol))   # code / prefix_code / bit2_code
        row.update(values_before=json.dumps([float(x) for x in before]),
                   values_after=json.dumps([float(x) for x in after]),
                   clipped=clipped)

        agg, _, _, dropped, anomalies = paired_analyse(model, sol, 0.05)
        row.update(n_pairs=agg['n_pairs'], n_R=agg['n_R'], n_F=agg['n_F'],
                   n_dropped=agg['n_dropped'], pattern=agg['pattern'],
                   same_type_adjacent=sum(1 for a in anomalies
                                          if a['reason'] == 'same type adjacent'),
                   windows_total=agg['n_windows'])
        for f in (0.01, 0.05, 0.10):
            a, _, _, _, _ = paired_analyse(model, sol, f)
            pc = int(round(f * 100))
            for key in ('L_rev', 'L_fwd', 'L_symmetric'):
                row[f'{pc}pct_{key}'] = a[key]['median']
                row[f'{pc}pct_{key}_min'] = a[key]['min']
                row[f'{pc}pct_{key}_max'] = a[key]['max']

        base = baselines[job['init']]
        readout_ok = bool(ch['cold_valid'] and ch['steady_valid'])
        causal_ok = bool(ch['causal_exactly_one_gate_and_flip'] and ch['causal_order'] and
                         ch['causal_dirs_alternate'])
        if not ch['certified']:
            outcome = 'lost_lock'
        elif ch['steady_sequence'] == base['steady_sequence']:
            outcome = 'recovered_original'
        else:
            outcome = 'legal_mod8_shift'
        row.update(readout_ok=readout_ok, causal_ok=causal_ok,
                   readout_ok_no_causal=bool(readout_ok and not causal_ok),
                   alternation_broken=bool(row['same_type_adjacent'] > 0),
                   alternation_restored=bool(row['pattern'] == CANONICAL_PATTERN),
                   outcome=outcome, ok=True, error='')
    except Exception as exc:                                     # noqa: BLE001
        row.update(ok=False, error=repr(exc), outcome='error')
    row['runtime_s'] = round(time.perf_counter() - started, 2)
    return row


def scan(args):
    all_jobs = jobs(args.phase)
    mine = shard_slice(all_jobs, args.shard, args.nshards)
    tag = f'shard{args.shard:02d}of{args.nshards:02d}'
    name = f'pool_perturbations_{args.phase}'
    states = load_states()
    baselines = load_baselines()
    rows = []
    for job in mine:
        r = run_one(job, states, baselines, args.hours, args.sample_min,
                    args.max_step_min, args.rtol, args.atol)
        rows.append(r)
        print(f"[{tag}] {r['kind']:<10} {r['mode']:<5} f={str(r['factor']):<5} "
              f"init={r['init_bits']} outcome={r['outcome']:<19} "
              f"cert={r.get('certified')} pairs={r.get('n_pairs')} "
              f"pattern={r.get('pattern')} L_sym5={r.get('5pct_L_symmetric')} "
              f"({r['runtime_s']}s)", flush=True)
    write_shard(name, tag, rows, dict(phase=args.phase, shard=args.shard,
                                      nshards=args.nshards, hours=args.hours,
                                      points=len(mine), model_sha256=sha256(ROOT_MODEL),
                                      states_json_sha256=sha256(STATES_JSON),
                                      baseline_csv_sha256=sha256(BASELINE_CSV),
                                      source_hashes=source_hashes()))
    print(f'[{tag}] wrote {len(rows)} rows')


def _change(vb, va):
    """How far the perturbation actually moved the targeted states.

    A pool that happens to sit at ~0 at the sampled instant is unchanged by a
    multiplicative factor, so such a run is NOT a test of that pool.  Recording
    the realised change keeps those runs from being read as evidence.
    """
    if not isinstance(vb, str) or not isinstance(va, str):
        return None
    try:
        a = np.asarray(json.loads(vb), dtype=float)
        b = np.asarray(json.loads(va), dtype=float)
    except (ValueError, TypeError):
        return None
    return float(np.max(np.abs(b - a))) if a.size else None


def merge(args):
    name = f'pool_perturbations_{args.phase}'
    df = merge_shards(name, args.nshards)
    df['max_abs_change'] = [_change(b, a) for b, a in zip(df.values_before, df.values_after)]
    df['changed_state'] = df.max_abs_change.apply(
        lambda v: None if v is None else bool(v > 1e-9))
    ok = df[df.ok.astype(bool)]
    verdict = dict(phase=args.phase, rows=len(df), ok=int(len(ok)),
                   hours=float(df.hours.dropna().iloc[0]) if len(df) and 'hours' in df else None,
                   model_sha256=sha256(ROOT_MODEL), source_sha256=source_hashes(),
                   baseline_csv_sha256=sha256(BASELINE_CSV),
                   canonical_pattern=CANONICAL_PATTERN,
                   overall=dict(
                       lost_lock=int((ok.outcome == 'lost_lock').sum()),
                       recovered_original=int((ok.outcome == 'recovered_original').sum()),
                       legal_mod8_shift=int((ok.outcome == 'legal_mod8_shift').sum()),
                       readout_ok_no_causal=int(ok.readout_ok_no_causal.astype(bool).sum()),
                       alternation_broken=int(ok.alternation_broken.astype(bool).sum()),
                       errors=int((~df.ok.astype(bool)).sum())),
                   by_perturbation={})

    for (pool, mode, factor), sub in ok.groupby(['kind', 'mode', 'factor'], dropna=False):
        key = f'{pool}|{mode}|{factor}'
        ls = pd.to_numeric(sub['5pct_L_symmetric'], errors='coerce').dropna()
        entry = dict(
            pool=pool, mode=mode,
            factor=(None if pd.isna(factor) else float(factor)),
            runs=int(len(sub)),
            recovered_original=int((sub.outcome == 'recovered_original').sum()),
            legal_mod8_shift=int((sub.outcome == 'legal_mod8_shift').sum()),
            lost_lock=int((sub.outcome == 'lost_lock').sum()),
            alternation_restored=int(sub.alternation_restored.astype(bool).sum()),
            alternation_broken=int(sub.alternation_broken.astype(bool).sum()),
            readout_ok_no_causal=int(sub.readout_ok_no_causal.astype(bool).sum()),
            clips_recorded=int(pd.to_numeric(sub.clipped, errors='coerce').fillna(0).sum()),
            L_symmetric_5pct_median=float(ls.median()) if len(ls) else None,
            L_symmetric_5pct_min=float(ls.min()) if len(ls) else None,
            L_symmetric_5pct_max=float(ls.max()) if len(ls) else None,
            L_symmetric_5pct_relative_spread=(float((ls.max() - ls.min()) / ls.median())
                                              if len(ls) and ls.median() else None),
            L_rev_5pct_median=float(pd.to_numeric(sub['5pct_L_rev'], errors='coerce').median()),
            L_fwd_5pct_median=float(pd.to_numeric(sub['5pct_L_fwd'], errors='coerce').median()),
            certified=int(sub.certified.astype(bool).sum()),
            degenerate_runs=int((sub.changed_state == False).sum()),  # noqa: E712
            max_abs_change_min=float(pd.to_numeric(sub.max_abs_change, errors='coerce').min())
            if len(sub) else None,
            max_abs_change_max=float(pd.to_numeric(sub.max_abs_change, errors='coerce').max())
            if len(sub) else None)
        entry['phase_dependent_leak'] = bool(
            entry['L_symmetric_5pct_relative_spread'] is not None and
            entry['L_symmetric_5pct_relative_spread'] > 0.05)
        verdict['by_perturbation'][key] = entry

    sens = {k: v for k, v in verdict['by_perturbation'].items()
            if v['lost_lock'] or v['readout_ok_no_causal'] or v['alternation_broken']}
    verdict['sensitive_perturbations'] = sorted(sens)
    verdict['any_failure'] = bool(len(sens) > 0)
    verdict['all_eight_intact_under_medium'] = bool(
        all(v['recovered_original'] + v['legal_mod8_shift'] == v['runs']
            for k, v in verdict['by_perturbation'].items()
            if v['mode'] == 'mult'))
    df.to_csv(OUT / f'{name}_all.csv', index=False, encoding='utf-8')
    (OUT / f'{name}_verdict.json').write_text(
        json.dumps(verdict, ensure_ascii=False, indent=2), encoding='utf-8')
    write_manifest()
    print(json.dumps(verdict, ensure_ascii=False, indent=2))


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    s = sub.add_parser('scan')
    s.add_argument('--phase', choices=('medium', 's2add', 'all'), required=True)
    s.add_argument('--shard', type=int, required=True)
    s.add_argument('--nshards', type=int, required=True)
    s.add_argument('--hours', type=float, default=600.0)
    s.add_argument('--sample-min', type=float, default=2.0)
    s.add_argument('--max-step-min', type=float, default=2.0)
    s.add_argument('--rtol', type=float, default=2e-7)
    s.add_argument('--atol', type=float, default=2e-9)
    m = sub.add_parser('merge')
    m.add_argument('--phase', choices=('medium', 's2add', 'all'), required=True)
    m.add_argument('--nshards', type=int, required=True)
    args = ap.parse_args()
    (scan if args.cmd == 'scan' else merge)(args)


if __name__ == '__main__':
    main()
