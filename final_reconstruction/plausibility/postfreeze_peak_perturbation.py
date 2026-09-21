"""POST-FREEZE validation: bit2 pulse-phase molecular-pool perturbations.

This is an ATTACHMENT to the frozen working point, not a modification of it.
Nothing in plausibility/threebit51_selected_v1.json, threebit51_provisional_n6.json
or threebit51_leak_correction.json is touched or rewritten; all output goes to
threebit51_results/postfreeze_peak_perturbation/.

Why it exists
-------------
The freeze-time sweep (plausibility/pool_perturbations_all_all.csv) sampled its
initial conditions at the FLUX TROUGH, where the Int2, RDF2 and complex pools sit
near their own minima.  The realised perturbation magnitudes there were as low as
1.3e-07 (Int2), 1.4e-05 (RDF2) and 3.0e-12 (C2) - those pools were never actually
tested.  This supplement perturbs each pool at a phase where it is LARGE.

Design
------
Events: from the LATE stable period-8 orbit of the frozen working point
(n_A1_gate = 6, carry1 mRNA 2 min, A1/F1 maturation 32.5 min), one F-type carry
(S2 low at the carry start -> the carry performs the forward write) and one R-type
carry (S2 high -> the reverse write).  The startup segment is never used.

Perturbation instant per object, recorded explicitly:
    Int2_full  (M_I2, I2_u, I2)   at the Int2 mature-protein peak inside the carry
    RDF2_full  (M_R2, R2_u, R2)   at the g1 main-window START
    C2         (explicit complex) at the C2 peak inside the carry

Factors 0.5 / 0.8 / 1.2 / 1.5 / 2.0 x both carry directions
    -> 3 objects x 2 directions x 5 factors = 30 perturbed runs
plus one unperturbed reference from each of the 6 (object, direction) instants
    -> 36 runs total, sharing ONE 600 h warm-up run done in `prepare`.

Horizon 450 h, not the suggested 400 h: the binding requirement is ">= 32 read
windows after the drop-8 cut", the clock period is 10.60 h, so 32 steady reads
need >= 40 windows ~= 429 h.  400 h would give only ~29.

Degeneracy rule (as specified): if the pool value at the perturbation instant is
below 1 % of that pool's own cycle peak, the run is marked degenerate and is
excluded from the robustness conclusion.

Success criterion: `recovered_original` is primary; `legal_mod8_shift` is reported
SEPARATELY and never merged into it.

Stages
    prepare                      one warm-up run, locate the events, freeze y0
    scan --shard S --nshards N   integrate one shard of the 36 jobs
    merge --nshards N            classify and answer the seven questions
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
from plausibility_common import (frozen_carry0, frozen_extension,  # noqa: E402
                                 gate_windows, sha256, source_hashes)
from scan_strict_tolerance import build, readout_strings  # noqa: E402
from check_read_commitment import characterise  # noqa: E402
from scan_carry_pairing import analyse as paired_analyse, segment_rows  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
PLAUS = ROOT / 'plausibility'
OUT = ROOT / 'threebit51_results' / 'postfreeze_peak_perturbation'
EVENTS_FILE = OUT / 'carry_events_selected.json'
SELECTED = dict(n_A1_gate=6.0, mrna=2.0, mat=32.5)
FROZEN_INPUTS = ('threebit51_selected_v1.json', 'threebit51_provisional_n6.json',
                 'threebit51_leak_correction.json')
MODEL_FILE = ROOT / 'model_threebit51.py'

HOURS = 450.0
SAMPLE_MIN = 2.0
MAX_STEP_MIN = 2.0
RTOL, ATOL = 2e-7, 2e-9
WARMUP_H = 250.0
MIN_STEADY_READS = 32
DEGENERACY_FRACTION = 1e-3

OBJECTS = (
    dict(obj='Int2_full', idx=(34, 35, 36), instant='int2_peak'),
    dict(obj='RDF2_full', idx=(40, 41, 42), instant='g1_window_start'),
    dict(obj='C2', idx=(43,), instant='c2_peak'),
)
DIRECTIONS = ('F', 'R')
FACTORS = (0.5, 0.8, 1.2, 1.5, 2.0)


def jobs():
    out = []
    for o in OBJECTS:
        for d in DIRECTIONS:
            out.append(dict(obj=o['obj'], idx=list(o['idx']), instant=o['instant'],
                            direction=d, factor=None))
    for o in OBJECTS:
        for d in DIRECTIONS:
            for f in FACTORS:
                out.append(dict(obj=o['obj'], idx=list(o['idx']), instant=o['instant'],
                                direction=d, factor=f))
    return out


def job_id(j):
    f = 'baseline' if j['factor'] is None else f"x{j['factor']:g}"
    return f"{j['obj']}|{j['direction']}|{f}"


def strict_mod8(seq):
    if not isinstance(seq, str) or not seq or 'x' in seq:
        return False
    v = [int(c) for c in seq]
    return all((v[i + 1] - v[i]) % 8 == 1 for i in range(len(v) - 1))


# ------------------------------------------------------------------- prepare
def prepare(args):
    OUT.mkdir(parents=True, exist_ok=True)
    model = build(SELECTED)
    sol = model.simulate(hours=args.warmup_h_total, sample_min=SAMPLE_MIN,
                         max_step_min=MAX_STEP_MIN)
    t, y = sol.t, sol.y
    sig = model.diagnostic_signals(y)
    wins = gate_windows(t, sig['g1'], threshold=0.01, min_duration_h=0.05)
    rows = segment_rows(t, sig['J_rev2'], sig['J_fwd2'], sig['S2'], y[36], wins, 0.05)
    late = [r for r in rows if r['start_h'] >= WARMUP_H and r['far_off_evaluable']]
    if not late:
        raise SystemExit('no late evaluable carry window found')
    chosen = {}
    for want in DIRECTIONS:
        cand = [r for r in late if r['type'] == want]
        if not cand:
            raise SystemExit(f'no late {want}-type carry found')
        chosen[want] = cand[-1]

    tail = t > 0.3 * t[-1]
    cycle_peak = {str(i): float(np.max(y[i][tail])) for i in range(y.shape[0])}

    events = {}
    for d, r in chosen.items():
        i0, i1 = r['i0'], r['i1']
        g1pk = int(i0 + np.argmax(sig['g1'][i0:i1 + 1]))
        rec = dict(direction=d, carry_type=r['type'], window_start_h=r['start_h'],
                   window_end_h=r['end_h'], S2_at_carry_start=r['S2_start'],
                   i0=int(i0), i1=int(i1), starts_after_h=WARMUP_H,
                   from_startup_segment=False, g1_peak_h=float(t[g1pk]))
        j = r['k'] - 1
        prev_i1 = int(wins[j]['i1']) if j >= 0 else 0
        rec['carry_period_start_h'] = float(t[prev_i1])
        for o in OBJECTS:
            # Documented rule: perturb each pool at its OWN peak inside the carry
            # period, defined as [end of the previous gate window, end of this one].
            # The instants named in the task spec (g1 window start, in-gate peak) are
            # reported as alternatives; for RDF2 and C2 during an F-type carry they
            # are degenerate because those pools are built up during the PRIOR dwell
            # and consumed by the F carry itself.
            best_k, best_v = None, -np.inf
            for i in o['idx']:
                kk = prev_i1 + int(np.argmax(y[i][prev_i1:i1 + 1]))
                if float(y[i][kk]) > best_v:
                    best_k, best_v = int(kk), float(y[i][kk])
            k = best_k
            peaks = [float(cycle_peak[str(x)]) for x in o['idx']]
            vals = [float(y[x][k]) for x in o['idx']]
            pmax = max(peaks)

            def ratio_at(a, b, _o=o, _pmax=pmax):
                v = max(float(np.max(y[x][a:b + 1])) for x in _o['idx'])
                return (v / _pmax) if _pmax > 0 else None

            rec[o['obj']] = dict(
                definition='self_peak_in_carry_period',
                definition_note=('peak of this pool over the carry period '
                                 '[previous gate window end, this gate window end]; the '
                                 'per-object instants named in the task specification are '
                                 'reported as ratio_at_* alternatives'),
                index=int(k), time_h=float(t[k]),
                carry_period_start_h=float(t[prev_i1]), carry_period_end_h=float(t[i1]),
                state=[float(v) for v in y[:, k]],
                values=vals, cycle_peaks=peaks,
                max_value=float(max(vals)), max_cycle_peak=pmax,
                ratio_to_cycle_peak=(best_v / pmax) if pmax > 0 else None,
                ratio_at_g1_window_start=ratio_at(i0, i0),
                ratio_at_gate_in_window_peak=ratio_at(i0, i1),
                ratio_at_g1_peak=ratio_at(g1pk, g1pk),
                at_window_boundary=bool(k in (i0, i1, prev_i1)))
        rec['rdf2_at_g1_peak'] = float(y[42][g1pk])
        events[d] = rec

    EVENTS_FILE.write_text(json.dumps(dict(
        selected=SELECTED, warmup_h=WARMUP_H, warmup_total_h=args.warmup_h_total,
        hours=HOURS, sample_min=SAMPLE_MIN, min_steady_reads=MIN_STEADY_READS,
        degeneracy_fraction=DEGENERACY_FRACTION,
        horizon_note=('450 h rather than the suggested 400 h: >= 32 steady reads after the '
                      'drop-8 cut need >= 40 windows at the 10.60 h clock period'),
        model_sha256=sha256(MODEL_FILE),
        events=events, cycle_peaks=cycle_peak), ensure_ascii=False, indent=2),
        encoding='utf-8')
    print(json.dumps({k: {kk: vv for kk, vv in v.items() if kk != 'state'}
                      for k, v in events.items()}, ensure_ascii=False, indent=2))
    print(f'wrote {EVENTS_FILE}')


def load_events():
    if not EVENTS_FILE.exists():
        raise SystemExit(f'missing {EVENTS_FILE} - run `prepare` first')
    return json.loads(EVENTS_FILE.read_text(encoding='utf-8'))


# ---------------------------------------------------------------------- scan
def evaluate(job, events, hours):
    started = time.perf_counter()
    ev = events['events'][job['direction']][job['obj']]
    y0 = np.asarray(ev['state'], dtype=float).copy()
    idx = list(job['idx'])
    before = y0[idx].copy()
    if job['factor'] is not None:
        y0[idx] = y0[idx] * float(job['factor'])
    after = y0[idx].copy()
    peaks = np.asarray(ev['cycle_peaks'], dtype=float)
    ratio = float(np.max(after / np.where(peaks > 0, peaks, np.inf)))

    row = dict(job_id=job_id(job), obj=job['obj'], direction=job['direction'],
               factor=job['factor'], definition=ev['definition'],
               source_index=int(ev['index']), source_time_h=ev['time_h'],
               lead_time_before_gate_h=float(
                   events['events'][job['direction']]['window_start_h'] - ev['time_h']),
               S2_at_carry_start=ev.get('S2_at_carry_start'),
               values_before=json.dumps([float(x) for x in before]),
               values_after=json.dumps([float(x) for x in after]),
               cycle_peaks=json.dumps([float(x) for x in peaks]),
               ratio_to_cycle_peak=ratio,
               degenerate=bool(job['factor'] is not None
                               and ratio < DEGENERACY_FRACTION),
               hours=hours)
    try:
        m = build(SELECTED)
        s = m.simulate(hours=hours, sample_min=SAMPLE_MIN, max_step_min=MAX_STEP_MIN,
                       rtol=RTOL, atol=ATOL, initial_state=y0)
        ch = characterise(m, s, hours, p=None, crosscheck=False)
        row.update(certified=ch['certified'], crossings=ch['crossings'],
                   gates=ch['gates'], reverse_events=ch['reverse_events'],
                   cold_reads=ch['cold_reads'], steady_reads=ch['steady_reads'],
                   cold_sequence=ch['cold_sequence'],
                   steady_sequence=ch['steady_sequence'],
                   steady_min_commitment_all_bits=ch['steady_min_commitment_all_bits'],
                   unlabelled_reads_all=ch['unlabelled_reads_all'],
                   bit2_margin_to_band_low=ch['drop8_bit2_margin_to_band_low'],
                   bit2_margin_to_band_high=ch['drop8_bit2_margin_to_band_high'],
                   setup_min_h_global=ch['setup_min_h_global'],
                   hold_min_h_global=ch['hold_min_h_global'],
                   causal_exactly_one_gate_and_flip=ch['causal_exactly_one_gate_and_flip'],
                   causal_order=ch['causal_order'],
                   causal_dirs_alternate=ch['causal_dirs_alternate'],
                   causal_unassigned_gates=ch['causal_unassigned_gates'],
                   causal_unassigned_crossings=ch['causal_unassigned_crossings'],
                   gate_contrast_passed=ch['gate_contrast_passed'])
        row.update(readout_strings(m, s))

        # frozen 34-state projection: integrate the frozen prefix model itself
        s34 = m.twobit.simulate(hours=hours, sample_min=SAMPLE_MIN,
                                max_step_min=MAX_STEP_MIN, rtol=RTOL, atol=ATOL,
                                initial_state=y0[:34])
        row['prefix_traj_gap'] = float(np.max(np.abs(s.y[:34] - s34.y)))
        row['prefix_states'] = 34

        for f in (0.01, 0.05, 0.10):
            a, _, _, _, anom = paired_analyse(m, s, f)
            pc = int(round(f * 100))
            for key in ('L_rev', 'L_fwd', 'L_symmetric'):
                row[f'{pc}pct_{key}'] = a[key]['median']
            row[f'{pc}pct_n_pairs'] = a['n_pairs']
            row[f'{pc}pct_same_type_adjacent'] = sum(
                1 for x in anom if x['reason'] == 'same type adjacent')
            row[f'{pc}pct_pattern'] = a['pattern']

        readout_ok = bool(strict_mod8(row['steady_sequence']))
        causal_ok = bool(ch['causal_exactly_one_gate_and_flip'] and ch['causal_order'] and
                         ch['causal_dirs_alternate'] and
                         ch['causal_unassigned_gates'] == 0 and
                         ch['causal_unassigned_crossings'] == 0)
        row.update(readout_ok=readout_ok, causal_ok=causal_ok,
                   readout_ok_but_causal_fail=bool(readout_ok and not causal_ok),
                   enough_steady_reads=bool(ch['steady_reads'] >= MIN_STEADY_READS),
                   ok=True, error='')
    except Exception as exc:                                     # noqa: BLE001
        row.update(ok=False, error=repr(exc))
    row['runtime_s'] = round(time.perf_counter() - started, 2)
    return row


def scan(args):
    events = load_events()
    allj = jobs()
    mine = allj[args.shard::args.nshards]
    tag = f'shard{args.shard:02d}of{args.nshards:02d}'
    rows = [evaluate(j, events, args.hours) for j in mine]
    for r in rows:
        print(f"[{tag}] {r['job_id']:<28} cert={r.get('certified')} "
              f"reads={r.get('steady_reads')} degen={r.get('degenerate')} "
              f"ratio={r.get('ratio_to_cycle_peak')} "
              f"L_sym5={r.get('5pct_L_symmetric')} ({r['runtime_s']}s)", flush=True)
    OUT.mkdir(parents=True, exist_ok=True)
    csv = OUT / f'postfreeze_{tag}.csv'
    pd.DataFrame(rows).to_csv(csv, index=False, encoding='utf-8')
    meta = dict(kind='postfreeze_peak_perturbation', shard=args.shard,
                nshards=args.nshards, jobs=len(rows), hours=args.hours,
                model_sha256=sha256(MODEL_FILE),
                selected_v1_sha256=sha256(PLAUS / 'threebit51_selected_v1.json'),
                provisional_sha256=sha256(PLAUS / 'threebit51_provisional_n6.json'),
                correction_sha256=sha256(PLAUS / 'threebit51_leak_correction.json'),
                events_file_sha256=sha256(EVENTS_FILE),
                source_hashes=source_hashes(),
                csv=csv.name, csv_sha256=sha256(csv))
    (OUT / f'meta_postfreeze_{tag}.json').write_text(
        json.dumps(meta, ensure_ascii=False, indent=2), encoding='utf-8')
    print(f'[{tag}] wrote {len(rows)} rows')


# --------------------------------------------------------------------- merge
def classify(rows):
    def is_base(r):
        """Baseline rows carry factor = NaN after the CSV round-trip, not None."""
        f = r.get('factor')
        if f is None:
            return True
        try:
            return bool(pd.isna(f))
        except (TypeError, ValueError):
            return False

    base = {(r['direction'], r['obj']): r for r in rows if is_base(r) and r.get('ok')}
    for r in rows:
        if is_base(r):
            r['outcome'] = 'baseline'
            r['baseline_code'] = r.get('steady_sequence')
            continue
        if not r.get('ok'):
            r['outcome'] = 'error'
            continue
        b = base.get((r['direction'], r['obj']))
        r['baseline_code'] = b.get('steady_sequence') if b else None
        r['baseline_certified'] = bool(b.get('certified')) if b else None
        r['code_matches_baseline'] = (None if b is None
                                      else bool(str(r.get('code')) == str(b.get('code'))))
        if b is None or not b.get('certified'):
            r['outcome'] = 'no_baseline'
        elif not r.get('certified'):
            r['outcome'] = 'lost_lock'
        elif r.get('steady_sequence') == b.get('steady_sequence'):
            r['outcome'] = 'recovered_original'
        else:
            r['outcome'] = 'legal_mod8_shift'
    return rows


def phase_structure():
    """Where each pool actually sits during each carry type.

    This is why the perturbation instant is the pool's own peak inside the carry
    period rather than the in-gate instants named in the task spec: during an
    F-type carry, RDF2 and C2 are already depleted, so an in-gate perturbation
    would be degenerate by the stated 1 % rule.  They are built up during the
    PRIOR dwell and consumed by the F carry itself.
    """
    if not EVENTS_FILE.exists():
        return None
    ev = json.loads(EVENTS_FILE.read_text(encoding='utf-8'))['events']
    out = {}
    for d, rec in ev.items():
        gate = rec['window_start_h']
        for obj in ('Int2_full', 'RDF2_full', 'C2'):
            v = rec[obj]
            out[f'{obj}|{d}'] = dict(
                chosen_instant_h=v['time_h'],
                chosen_ratio_to_cycle_peak=v['ratio_to_cycle_peak'],
                lead_time_before_gate_h=float(gate - v['time_h']),
                ratio_at_g1_window_start=v['ratio_at_g1_window_start'],
                ratio_at_gate_in_window_peak=v['ratio_at_gate_in_window_peak'],
                degenerate_at_spec_named_instant=bool(
                    (v['ratio_at_g1_window_start'] or 0) < DEGENERACY_FRACTION
                    and (v['ratio_at_gate_in_window_peak'] or 0) < DEGENERACY_FRACTION))
    return out


def merge(args):
    frames = []
    for s in range(args.nshards):
        p = OUT / f'postfreeze_shard{s:02d}of{args.nshards:02d}.csv'
        if not p.exists():
            raise SystemExit(f'missing shard {p}')
        frames.append(pd.read_csv(p))
    df = pd.DataFrame(classify(pd.concat(frames, ignore_index=True).to_dict('records')))
    df.to_csv(OUT / 'postfreeze_all.csv', index=False, encoding='utf-8')

    pert = df[(df.factor.notna()) & (df.ok.astype(bool))]
    eff = pert[~pert.degenerate.astype(bool)]
    answers = {}
    for o in OBJECTS:
        for d in DIRECTIONS:
            sub = eff[(eff.obj == o['obj']) & (eff.direction == d)]
            rec = sub[sub.outcome == 'recovered_original']
            answers[f"{o['obj']}|{d}"] = dict(
                effective_runs=int(len(sub)),
                recovered_original=int(len(rec)),
                legal_mod8_shift=int((sub.outcome == 'legal_mod8_shift').sum()),
                lost_lock=int((sub.outcome == 'lost_lock').sum()),
                max_recovered_factor=(float(rec.factor.max()) if len(rec) else None),
                min_recovered_factor=(float(rec.factor.min()) if len(rec) else None),
                all_tested_factors_recovered=bool(len(rec) == len(sub) and len(sub) > 0))
    losses = int((eff.outcome == 'lost_lock').sum())
    shifts = int((eff.outcome == 'legal_mod8_shift').sum())
    q = {}
    for o in OBJECTS:
        for d in DIRECTIONS:
            q[f"{o['obj']} ({d}-type carry) max recoverable factor"] = \
                answers[f"{o['obj']}|{d}"]['max_recovered_factor']
    q['2_most_sensitive_pool'] = ('none - no pool lost lock at any tested effective factor'
                                  if losses == 0 else 'see by_object_and_direction')
    q['3_direction_asymmetry'] = {
        o['obj']: dict(F=answers[f"{o['obj']}|F"]['max_recovered_factor'],
                       R=answers[f"{o['obj']}|R"]['max_recovered_factor'])
        for o in OBJECTS}
    q['4_legal_phase_shifts'] = shifts
    q['5_readout_ok_but_causal_fail'] = int(
        eff.readout_ok_but_causal_fail.astype(bool).sum())
    q['6_selected_v1_needs_v2'] = bool(losses > 0)
    q['7_hash_traceability'] = dict(
        frozen_files_modified=False,
        selected_v1=sha256(PLAUS / 'threebit51_selected_v1.json'),
        provisional=sha256(PLAUS / 'threebit51_provisional_n6.json'),
        correction=sha256(PLAUS / 'threebit51_leak_correction.json'),
        model_threebit51=sha256(MODEL_FILE),
        attachment_manifest='SHA256SUMS.json in this directory')

    verdict = dict(
        status='post_freeze_validation',
        task='51-state three-bit model: bit2 pulse-phase molecular-pool perturbation supplement',
        attachment_to='plausibility/threebit51_selected_v1.json',
        frozen_inputs=q['7_hash_traceability'],
        design=dict(hours=HOURS, sample_min=SAMPLE_MIN, warmup_h=WARMUP_H,
                    min_steady_reads=MIN_STEADY_READS,
                    degenerate_threshold=DEGENERACY_FRACTION,
                    horizon_note=('450 h rather than the suggested 400 h: >= 32 steady reads '
                                  'after the drop-8 cut need >= 40 windows at the 10.60 h '
                                  'clock period'),
                    outcomes=dict(
                        total_runs=int(len(df)), ok=int(df.ok.astype(bool).sum()),
                        baselines=int(df.factor.isna().sum()), perturbed=int(len(pert)),
                        degenerate=int(pert.degenerate.astype(bool).sum()),
                        effective=int(len(eff)),
                        recovered_original=int((eff.outcome == 'recovered_original').sum()),
                        legal_mod8_shift=shifts, lost_lock=losses,
                        readout_ok_but_causal_fail=q['5_readout_ok_but_causal_fail'],
                        all_enough_steady_reads=bool(eff.enough_steady_reads.astype(bool).all()),
                        min_steady_reads_observed=int(pd.to_numeric(
                            eff.steady_reads, errors='coerce').min()))),
        by_object_and_direction=answers,
        pool_phase_structure=phase_structure(),
        answers_to_the_seven_questions=q,
        per_run={r['job_id']: dict(
            obj=r['obj'], direction=r['direction'], factor=r['factor'],
            definition=r['definition'], ratio_to_cycle_peak=_f(r.get('ratio_to_cycle_peak')),
            degenerate=_b(r.get('degenerate')), outcome=r.get('outcome'),
            certified=_b(r.get('certified')), steady_reads=_i(r.get('steady_reads')),
            steady_sequence=r.get('steady_sequence'), code=r.get('code'),
            readout_ok=_b(r.get('readout_ok')), causal_ok=_b(r.get('causal_ok')),
            readout_ok_but_causal_fail=_b(r.get('readout_ok_but_causal_fail')),
            commitment=_f(r.get('steady_min_commitment_all_bits')),
            margin_low=_f(r.get('bit2_margin_to_band_low')),
            margin_high=_f(r.get('bit2_margin_to_band_high')),
            setup=_f(r.get('setup_min_h_global')), hold=_f(r.get('hold_min_h_global')),
            same_type_adjacent_5pct=_i(r.get('5pct_same_type_adjacent')),
            L_rev_5pct=_f(r.get('5pct_L_rev')), L_fwd_5pct=_f(r.get('5pct_L_fwd')),
            L_symmetric_1pct=_f(r.get('1pct_L_symmetric')),
            L_symmetric_5pct=_f(r.get('5pct_L_symmetric')),
            L_symmetric_10pct=_f(r.get('10pct_L_symmetric')),
            prefix_traj_gap=_f(r.get('prefix_traj_gap')),
            error=r.get('error')) for _, r in df.iterrows()},
    )
    (OUT / 'postfreeze_verdict.json').write_text(
        json.dumps(verdict, ensure_ascii=False, indent=2), encoding='utf-8')
    write_local_manifest()
    print(json.dumps({k: verdict[k] for k in
                      ('design', 'by_object_and_direction',
                       'answers_to_the_seven_questions')}, ensure_ascii=False, indent=2))


def _f(v):
    return None if v is None or (not isinstance(v, str) and pd.isna(v)) else float(v)


def _i(v):
    if v is None or (not isinstance(v, str) and pd.isna(v)):
        return None
    try:
        return int(v)
    except (TypeError, ValueError):
        return None


def _b(v):
    if v is None or isinstance(v, str):
        return v
    try:
        if pd.isna(v):
            return None
    except (TypeError, ValueError):
        pass
    return bool(v)


def write_local_manifest():
    OUT.mkdir(parents=True, exist_ok=True)
    man = {}
    for p in sorted(OUT.rglob('*')):
        if p.is_file() and not p.name.startswith('SHA256SUMS'):
            man[p.name] = sha256(p)
    for name, h in source_hashes().items():
        man[f'frozen_source::{name}'] = h
    for name in FROZEN_INPUTS:
        q = PLAUS / name
        if q.exists():
            man[name] = sha256(q)
    man['model_threebit51.py'] = sha256(MODEL_FILE)
    (OUT / 'SHA256SUMS.json').write_text(
        json.dumps(man, ensure_ascii=False, indent=2, sort_keys=True), encoding='utf-8')
    (OUT / 'SHA256SUMS.txt').write_text(
        ''.join(f'{v}  {k}\n' for k, v in sorted(man.items())), encoding='utf-8')
    return man


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('prepare')
    p.add_argument('--warmup-h-total', type=float, default=600.0)
    s = sub.add_parser('scan')
    s.add_argument('--shard', type=int, required=True)
    s.add_argument('--nshards', type=int, required=True)
    s.add_argument('--hours', type=float, default=HOURS)
    m = sub.add_parser('merge')
    m.add_argument('--nshards', type=int, required=True)
    args = ap.parse_args()
    if args.cmd == 'prepare':
        prepare(args)
    elif args.cmd == 'scan':
        scan(args)
    else:
        merge(args)


if __name__ == '__main__':
    main()
