"""Eight phase-consistent digital initial states of the 51-state three-bit counter.

Construction (mirrors `verify_twobit34_initial_states.py`, the ratified two-bit
precedent, so the two initial-state sets are built the same way)
----------------------------------------------------------------------------
1. Warm up the n_A1_gate = 6 model and take the late stable period-8 orbit.
2. Sample one COMPLETE 51-state biochemical vector at the flux trough of each
   digital value 0..7 - not just S0/S1/S2.  Because every flux trough sits at the
   same phase of the upstream oscillator, the eight samples are already
   phase-consistent modulo one clock cycle; the measured spread is recorded.
3. Reset the six upstream oscillator states AND the shared C31 mRNA to ONE common
   clock phase: indices 0..6, taken from the value-0 sample.
   Index 6 is `b0_M_I`, and `model.py` states "Bit0 M_I IS the upstream C31 mRNA,
   not a second transcript pool", so it is a shared quantity and must not carry
   eight different phases.  S0/S1/S2 are NOT edited - they come from the orbit and
   already encode the digital value.

Read-window convention: the state is sampled at a trough, so the first newly
observable read after the restart is one clock period later, i.e. the sequence
must be (value + i + 1) mod 8 for i = 0, 1, 2, ...

Stages
------
    guard                              structural no-feedback / prefix identity
    extract --warmup-h 300 --cutoff-h 180
    verify --shard V --nshards 8 --hours 600
    merge --nshards 8

`verify` reuses `characterise()` from `check_read_commitment.py`, so the
commitment, band-edge, setup/hold, causal 1:1 and three-segment leak definitions
are byte-for-byte the same instrument that produced the two-candidate margins.
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
                                 merge_shards, sha256, source_hashes,
                                 write_manifest, write_shard)
from model_threebit51 import (STATE_NAMES_51, ThreeBit51Model,  # noqa: E402
                              ThreeBitCarryParameters)
from model_twobit34 import CarryExpressionParameters  # noqa: E402
from check_read_commitment import characterise  # noqa: E402

NAME = 'eight_initial'
N_VALUES = 8
ROOT_MODEL = Path(__file__).resolve().parents[1] / 'model_threebit51.py'
STATES_JSON = OUT / 'eight_initial_states.json'
STATES_CSV = OUT / 'eight_initial_states.csv'
PREFIX_GUARD = OUT / 'eight_initial_prefix_guard.json'

# the selected carry-2 candidate, fixed before this stage
SELECTED = dict(n_A1_gate=6.0, mrna=2.0, mat=32.5)
PHASE_RESET_INDICES = tuple(range(7))   # 6 oscillator states + shared C31 mRNA


def build_selected():
    carry1 = CarryExpressionParameters(
        mrna_half_life_min=SELECTED['mrna'],
        activator_maturation_half_life_min=SELECTED['mat'],
        repressor_maturation_half_life_min=SELECTED['mat'])
    return ThreeBit51Model(extension=frozen_extension(),
                           carry=ThreeBitCarryParameters(carry0=frozen_carry0(),
                                                         carry1=carry1),
                           n_A1_gate=SELECTED['n_A1_gate'])


def decode_initial(y):
    """3-bit digital value read from S0, S1, S2 at the same instant."""
    return (int(y[16] >= 0.5) + 2 * int(y[27] >= 0.5) + 4 * int(y[44] >= 0.5))


def expected_sequence(initial_value, length):
    return ''.join(str((initial_value + i + 1) % 8) for i in range(length))


# ------------------------------------------------------- structural prefix guard
def structural_guard(sample_min=2.0, max_step_min=2.0, hours=60.0):
    """States 0..33 must have no dependence on states 34..50, at all.

    Two exact identities are tested, both of which must hold with gap == 0.0:
      (a) the prefix of the 51-state derivative IS the frozen 34-state
          derivative evaluated on the same prefix;
      (b) changing states 34..50 by scaling AND by additive offsets leaves the
          prefix derivative bit-identical.
    """
    model = build_selected()
    m34 = model.twobit                       # the very object rhs() delegates to
    sol = model.simulate(hours=hours, sample_min=sample_min, max_step_min=max_step_min)
    probes = [model.initial_state(cold=False), sol.y[:, 0],
              sol.y[:, len(sol.t) // 2], sol.y[:, -1]]

    worst_identity = 0.0
    worst_feedback = 0.0
    rows = []
    for i, p in enumerate(probes):
        p = np.asarray(p, dtype=float)
        d_ref = model.rhs(0.0, p)[:34]
        d_frozen = m34.rhs(0.0, p[:34])
        identity = float(np.max(np.abs(d_ref - d_frozen)))
        worst_identity = max(worst_identity, identity)
        feedback = 0.0
        for scale in (0.0, 0.5, 1.5, 3.0):
            y = p.copy()
            y[34:] *= scale
            feedback = max(feedback, float(np.max(np.abs(model.rhs(0.0, y)[:34] - d_ref))))
        for offset in (0.0, 1.0, 10.0, 1e3):
            y = p.copy()
            y[34:] += offset
            feedback = max(feedback, float(np.max(np.abs(model.rhs(0.0, y)[:34] - d_ref))))
        worst_feedback = max(worst_feedback, feedback)
        rows.append(dict(probe=i, prefix_identity_gap=identity,
                         feedback_scaling_and_offset_gap=feedback))
    result = dict(passed=bool(worst_identity == 0.0 and worst_feedback == 0.0),
                  prefix_derivative_identity_gap=worst_identity,
                  feedback_gap_scaling_and_offset=worst_feedback,
                  probes=rows, n_state=len(STATE_NAMES_51),
                  prefix_states=34,
                  model_sha256=sha256(ROOT_MODEL), source_sha256=source_hashes(),
                  note=('rhs() computes states 0..33 by calling TwoBit34Model.rhs on y[:34] '
                        'itself, so the no-feedback property is structural, not merely '
                        'observed. Both gaps are exact zeros.'))
    OUT.mkdir(parents=True, exist_ok=True)
    PREFIX_GUARD.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding='utf-8')
    write_manifest()
    return result


# ------------------------------------------------------------------ extraction
def extract(args):
    from verify_threebit51 import analyse_threebit

    model = build_selected()
    sol = model.simulate(hours=args.warmup_h, sample_min=args.sample_min,
                         max_step_min=args.max_step_min)
    analysis = analyse_threebit(model, sol, args.warmup_h)

    selected = {}
    for read in analysis['read_windows']:
        value = read['value']
        if value is None or read['trough_h'] < args.cutoff_h:
            continue
        if int(value) in selected:
            continue
        k = int(np.argmin(np.abs(sol.t - read['trough_h'])))
        selected[int(value)] = dict(value=int(value), time_h=float(sol.t[k]), index=k,
                                    cycle=int(read['cycle']),
                                    trough_h=float(read['trough_h']),
                                    window_h=[float(read['window_start_h']),
                                              float(read['window_end_h'])],
                                    state=sol.y[:, k].copy())
    if set(selected) != set(range(N_VALUES)):
        raise RuntimeError(f'could not extract all eight digital values: {sorted(selected)}')

    before = np.vstack([selected[v]['state'] for v in range(N_VALUES)])
    decode_before = [decode_initial(selected[v]['state']) for v in range(N_VALUES)]
    raw7 = before[:, :7]
    scale = np.maximum(np.max(np.abs(raw7), axis=0), 1e-12)
    spread_abs = np.ptp(raw7, axis=0)
    spread_rel = spread_abs / scale

    # one actually simulated upstream state, not the numerical mean
    common7 = selected[0]['state'][:7].copy()
    edits = []
    for v in range(N_VALUES):
        selected[v]['state_before_phase_reset'] = selected[v]['state'].copy()
        delta = common7 - selected[v]['state'][:7]
        edits.append(float(np.max(np.abs(delta))))
        selected[v]['state'][:7] = common7

    audit = dict(
        common_phase_source_value=0,
        phase_reset_indices=list(PHASE_RESET_INDICES),
        phase_reset_state_names=[STATE_NAMES_51[i] for i in PHASE_RESET_INDICES],
        sampled_times_h=[selected[v]['time_h'] for v in range(N_VALUES)],
        sampled_cycles=[selected[v]['cycle'] for v in range(N_VALUES)],
        decode_of_extracted_state=decode_before,
        decode_matches_value=bool(all(decode_before[v] == v for v in range(N_VALUES))),
        max_upstream_absolute_spread=float(spread_abs.max()),
        max_upstream_relative_spread=float(spread_rel.max()),
        upstream_absolute_spread=dict(zip([STATE_NAMES_51[i] for i in PHASE_RESET_INDICES],
                                          map(float, spread_abs))),
        upstream_relative_spread=dict(zip([STATE_NAMES_51[i] for i in PHASE_RESET_INDICES],
                                          map(float, spread_rel))),
        phase_edit_max_abs=float(max(edits)),
        note=('the eight troughs are one clock period apart, so the upstream block is '
              'already at a common phase; the edit magnitude above shows how much was '
              'actually changed. S0/S1/S2 are never edited.'))

    payload = dict(selected_profile=SELECTED, phase_audit=audit,
                   state_names=list(STATE_NAMES_51),
                   warmup_h=args.warmup_h, cutoff_h=args.cutoff_h,
                   initial_states={f'{v:03b}': [float(x) for x in selected[v]['state']]
                                   for v in range(N_VALUES)},
                   model_sha256=sha256(ROOT_MODEL), source_sha256=source_hashes(),
                   construction=('late orbit-derived complete 51-state biochemical states; '
                                 'six oscillator states plus shared C31 mRNA (index 6, '
                                 'b0_M_I) reset to one common clock phase'))
    STATES_JSON.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding='utf-8')

    table = {'state': list(STATE_NAMES_51)}
    for v in range(N_VALUES):
        table[f'initial_{v:03b}'] = selected[v]['state']
    pd.DataFrame(table).to_csv(STATES_CSV, index=False, encoding='utf-8')
    write_manifest()
    print(json.dumps(audit, ensure_ascii=False, indent=2))
    print(f'wrote {STATES_JSON.name} and {STATES_CSV.name}')


def load_states():
    if not STATES_JSON.exists():
        raise SystemExit(f'missing {STATES_JSON} - run `extract` first')
    payload = json.loads(STATES_JSON.read_text(encoding='utf-8'))
    return payload, {int(k, 2): np.asarray(v, dtype=float)
                     for k, v in payload['initial_states'].items()}


# ------------------------------------------------------------------ verification
def prefix_check(model, y0, sol51, hours, sample_min, max_step_min):
    """Does the frozen 34-state prefix evolve as the frozen model says it should?"""
    sol34 = model.twobit.simulate(hours=hours, sample_min=sample_min,
                                  max_step_min=max_step_min,
                                  initial_state=np.asarray(y0)[:34])
    gap = float(np.max(np.abs(sol51.y[:34] - sol34.y[:34])))
    return dict(prefix_states=34, prefix_traj_gap=gap,
                prefix_note=('the 51-state rhs evaluates states 0..33 with the frozen '
                             'TwoBit34Model derivative, so any residual here is the '
                             'adaptive step-size controller, not a coupling'))


def run_one(value, state, hours, sample_min, max_step_min, prefix=True):
    started = time.perf_counter()
    row = dict(initial_value=value, initial_bits=f'{value:03b}')
    try:
        model = build_selected()
        y0 = np.asarray(state, dtype=float)
        sol = model.simulate(hours=hours, sample_min=sample_min,
                             max_step_min=max_step_min, initial_state=y0)
        row.update(characterise(model, sol, hours, p=None, crosscheck=False))
        decoded = decode_initial(y0)
        cold = row['cold_sequence']
        expected = expected_sequence(value, len(cold))
        row.update(decoded_initial=decoded,
                   first_new_read=cold[:1],
                   first_new_read_correct=bool(cold[:1] == str((value + 1) % 8)),
                   expected_sequence=expected,
                   sequence_relation_correct=bool(decoded == value and cold == expected))
        if prefix:
            row.update(prefix_check(model, y0, sol, hours, sample_min, max_step_min))
        row.update(ok=True, error='')
        row['passed'] = bool(row['sequence_relation_correct'] and row['cold_valid'] and
                             row['steady_valid'] and row['certified'] and
                             row.get('causal_exactly_one_gate_and_flip') and
                             row.get('causal_order') and row.get('causal_dirs_alternate') and
                             row.get('causal_unassigned_gates') == 0 and
                             row.get('causal_unassigned_crossings') == 0)
    except Exception as exc:                                     # noqa: BLE001
        row.update(ok=False, error=repr(exc), passed=False)
    row['runtime_s'] = round(time.perf_counter() - started, 2)
    return row


def verify(args):
    _, states = load_states()
    values = list(range(args.shard, N_VALUES, args.nshards))
    tag = f'shard{args.shard:02d}of{args.nshards:02d}'
    rows = []
    for v in values:
        r = run_one(v, states[v], args.hours, args.sample_min, args.max_step_min,
                    prefix=not args.no_prefix_check)
        rows.append(r)
        print(f"[{tag}] {v:03b} decoded={r.get('decoded_initial')} "
              f"expected={r.get('expected_sequence')} observed={r.get('cold_sequence')} "
              f"passed={r.get('passed')} cert={r.get('certified')} "
              f"setup={r.get('setup_min_h_global')} hold={r.get('hold_min_h_global')} "
              f"leak5%={r.get('leak5pct_ratio')} ({r['runtime_s']}s)", flush=True)
    write_shard(NAME, tag, rows, dict(kind='scan', shard=args.shard, nshards=args.nshards,
                                      hours=args.hours, points=len(values),
                                      selected_profile=SELECTED,
                                      model_sha256=sha256(ROOT_MODEL),
                                      states_json_sha256=sha256(STATES_JSON),
                                      source_hashes=source_hashes()))
    print(f'[{tag}] wrote {len(rows)} rows')


def merge(args):
    df = merge_shards(NAME, args.nshards).sort_values('initial_value')
    ok = df[df.ok.astype(bool)]
    verdict = dict(rows=len(df), ok=int(len(ok)),
                   passed=int(ok.passed.astype(bool).sum()),
                   all_eight_pass=bool(len(ok) == N_VALUES and ok.passed.astype(bool).all()),
                   initial_states_sha256=sha256(STATES_JSON),
                   prefix_guard_sha256=(sha256(PREFIX_GUARD) if PREFIX_GUARD.exists() else None),
                   model_sha256=sha256(ROOT_MODEL), source_sha256=source_hashes(),
                   selected_profile=SELECTED,
                   per_value={}, )
    for _, r in df.iterrows():
        v = int(r.initial_value)
        verdict['per_value'][f'{v:03b}'] = dict(
            decoded_initial=None if pd.isna(r.get('decoded_initial')) else int(r.decoded_initial),
            expected=r.get('expected_sequence'), observed=r.get('cold_sequence'),
            first_new_read_correct=_b(r.get('first_new_read_correct')),
            sequence_relation_correct=_b(r.get('sequence_relation_correct')),
            certified=_b(r.get('certified')), passed=_b(r.get('passed')),
            cold_reads=None if pd.isna(r.get('cold_reads')) else int(r.cold_reads),
            steady_reads=None if pd.isna(r.get('steady_reads')) else int(r.steady_reads),
            steady_sequence=r.get('steady_sequence'),
            crossings=None if pd.isna(r.get('crossings')) else int(r.crossings),
            gates=None if pd.isna(r.get('gates')) else int(r.gates),
            reverse_events=None if pd.isna(r.get('reverse_events')) else int(r.reverse_events),
            setup_min_h_global=_f(r.get('setup_min_h_global')),
            hold_min_h_global=_f(r.get('hold_min_h_global')),
            bit0_setup_min_h=_f(r.get('bit0_setup_min_h')), bit0_hold_min_h=_f(r.get('bit0_hold_min_h')),
            bit1_setup_min_h=_f(r.get('bit1_setup_min_h')), bit1_hold_min_h=_f(r.get('bit1_hold_min_h')),
            bit2_setup_min_h=_f(r.get('bit2_setup_min_h')), bit2_hold_min_h=_f(r.get('bit2_hold_min_h')),
            bit0_commitment_min=_f(r.get('drop8_bit0_commitment_min')),
            bit1_commitment_min=_f(r.get('drop8_bit1_commitment_min')),
            bit2_commitment_min=_f(r.get('drop8_bit2_commitment_min')),
            bit0_margin_to_band_low=_f(r.get('drop8_bit0_margin_to_band_low')),
            bit1_margin_to_band_low=_f(r.get('drop8_bit1_margin_to_band_low')),
            bit2_margin_to_band_low=_f(r.get('drop8_bit2_margin_to_band_low')),
            bit0_margin_to_band_high=_f(r.get('drop8_bit0_margin_to_band_high')),
            bit1_margin_to_band_high=_f(r.get('drop8_bit1_margin_to_band_high')),
            bit2_margin_to_band_high=_f(r.get('drop8_bit2_margin_to_band_high')),
            leak1pct_ratio=_f(r.get('leak1pct_ratio')),
            leak5pct_ratio=_f(r.get('leak5pct_ratio')),
            leak10pct_ratio=_f(r.get('leak10pct_ratio')),
            leak5pct_far_off_rev=_f(r.get('leak5pct_far_off_rev')),
            leak5pct_far_off_fwd=_f(r.get('leak5pct_far_off_fwd')),
            prefix_traj_gap=_f(r.get('prefix_traj_gap')),
            error=r.get('error'))
    verdict['spread'] = dict(
        setup_min_h_global=_spread(ok.get('setup_min_h_global')),
        hold_min_h_global=_spread(ok.get('hold_min_h_global')),
        leak5pct_ratio=_spread(ok.get('leak5pct_ratio')),
        prefix_traj_gap=_spread(ok.get('prefix_traj_gap')))
    df.to_csv(OUT / f'{NAME}_all.csv', index=False, encoding='utf-8')
    (OUT / f'{NAME}_verdict.json').write_text(
        json.dumps(verdict, ensure_ascii=False, indent=2), encoding='utf-8')
    write_manifest()
    print(json.dumps(verdict, ensure_ascii=False, indent=2))


def _f(v):
    return None if v is None or pd.isna(v) else float(v)


def _b(v):
    """Boolean that survives a missing column or a NaN instead of becoming True."""
    if v is None:
        return None
    if isinstance(v, str):
        return v
    try:
        if pd.isna(v):
            return None
    except (TypeError, ValueError):
        pass
    return bool(v)


def _spread(s):
    if s is None:
        return None
    v = pd.to_numeric(s, errors='coerce').dropna()
    if v.empty:
        return None
    return dict(n=int(v.size), min=float(v.min()), median=float(v.median()),
                max=float(v.max()))


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='mode', required=True)
    g = sub.add_parser('guard')
    g.add_argument('--hours', type=float, default=60.0)
    e = sub.add_parser('extract')
    e.add_argument('--warmup-h', type=float, default=300.0)
    e.add_argument('--cutoff-h', type=float, default=180.0)
    e.add_argument('--sample-min', type=float, default=2.0)
    e.add_argument('--max-step-min', type=float, default=2.0)
    v = sub.add_parser('verify')
    v.add_argument('--shard', type=int, required=True)
    v.add_argument('--nshards', type=int, required=True)
    v.add_argument('--hours', type=float, default=600.0)
    v.add_argument('--sample-min', type=float, default=2.0)
    v.add_argument('--max-step-min', type=float, default=2.0)
    v.add_argument('--no-prefix-check', action='store_true')
    m = sub.add_parser('merge')
    m.add_argument('--nshards', type=int, required=True)
    args = ap.parse_args()
    if args.mode == 'guard':
        result = structural_guard(hours=args.hours)
        print(json.dumps(result, ensure_ascii=False, indent=2))
        if not result['passed']:
            raise SystemExit('structural prefix guard FAILED')
    elif args.mode == 'extract':
        extract(args)
    elif args.mode == 'verify':
        verify(args)
    else:
        merge(args)


if __name__ == '__main__':
    main()
