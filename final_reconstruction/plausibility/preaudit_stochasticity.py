"""Stochasticity PRE-AUDIT: what can be decided WITHOUT a stochastic simulator.

STATUS: written, NOT YET RUN.  Running it is a separate, explicit decision; the
"decision rule" below says in advance what a result would license.

Why a pre-audit at all
----------------------
Every claim in this project is a deterministic ODE prediction.  Before paying for
a Gillespie/tau-leaping implementation of a 51-state model with Hill kinetics
(which needs propensity choices the model does not currently make), it is worth
asking what a stochastic study could still change, and whether part of the answer
is already visible from deterministic information.

  CAN be answered here (cheap, deterministic):
    (1) Is the design in a SMALL-NUMBER regime at the decision instant?  Copy
        numbers follow from the model's OWN conversion `copies_per_au`
        (model.py:90, `602.214076 * uM_per_au` for a 1 fL cell), so this is a
        measurement, not an assumption.  A pool at 10 copies cannot be described
        by a continuous concentration.
    (2) How large a STATIC perturbation of the drive does the frozen design
        tolerate before certification is lost?  A fluctuation slow compared with
        one carry period (42.4 h) is quasi-static, so this bounds such
        fluctuations - from the optimistic side only.

  CANNOT be answered here (must be stated, never implied):
    (3) FAST fluctuation (period <~ one carry period): noise acting inside a
        single write/erase decision.  Needs a time-dependent or stochastic model.
    (4) PHASE DIFFUSION of the upstream clock.  An ODE has one phase; noise gives
        a phase distribution growing in time, and the read criterion samples a
        narrow window (trough +-10 %, 2.13 h).  Nothing deterministic bounds this.
    (5) PLASMID SEGREGATION / low-copy extinction, and bursty transcription.
    (6) Single-cell success RATE.  A deterministic run is one idealised cell; the
        engineering claim "a population counts" is a distribution statement.

So this script CANNOT clear the stochasticity question.  Its purpose is to decide
whether a full stochastic study is REQUIRED (small numbers, or a narrow static
tolerance) or merely DESIRABLE (large numbers, a wide tolerance).

Design
------
Three subcommands, all deterministic:

  plan        print the registered grid, cost and rules.  Runs nothing.

  inventory   one 600 h run at the frozen working point; convert every
              concentration-like state to copies via `copies_per_au` and report
              the smallest decision-relevant pool over the post-burn-in phase and
              during the gate windows.

  tolerance   a graded STATIC perturbation of ONE knob at a time, factor f in
              F_GRID around 1.0:
                conc_scale : f * uM_per_au       global gain: rescales every a.u.
                                                 concentration relative to the
                                                 production terms (model.py:90,136)
                clock_K    : f * clock_K_au      clock threshold (duty cycle)
                K_A1       : f * ZENG['K_A'][1]  gate threshold (ZENG patch)
              recording certification, crossings, gate count, read sequence, bit2
              margins and the paired leak at every f.

              `K_A1` patches the GLOBAL ZENG slot that BOTH the carry gate and
              the F1 production term read (see plausibility_common's decoupling
              note), so its effect is NOT attributable to one arm.  That is
              acceptable here because the question is "how much perturbation can
              the design absorb", not "which arm moved it" - but it is recorded,
              and the patch is restored immediately after each run.

  h2check     the three f = 1.0 construction paths run back-to-back IN ONE
              PROCESS and compared as full arrays: they must agree bit-for-bit
              (`max|d state| == 0.0`).  If they do not, a ZENG patch leaked
              between models (H2) and the tolerance scan is void.  This is a
              separate subcommand because it cannot be checked across shards.

Pre-registered decision rule (fixed before any run, not after)
--------------------------------------------------------------
Let W be the largest SYMMETRIC, CONNECTED certified factor interval around 1.0,
taken over the three knobs (the worst knob wins).

  * W <  0.05   -> SMALL: nothing but f = 1.0 certifies.  Static tolerance is
                   below one grid step; a stochastic study is REQUIRED before any
                   counting claim.
  * 0.05 <= W < 0.10 -> MODERATE: report the interval; do NOT claim noise
                   robustness.
  * W >= 0.10   -> WIDE: the deterministic conclusion tolerates slow
                   perturbations of at least this size; a stochastic study stays
                   DESIRABLE for (3)-(6), and `inventory` decides how urgent.

  CONSEQUENCE OF THE GRID (stated so the rule is not oversold): on
  F_GRID = (0.80, 0.90, 0.95, 1.00, 1.05, 1.10, 1.20) with steps of 0.05 around
  1.0, W can only take the values 0.00, 0.05, 0.10 or 0.20.  The rule is a
  three-way classification on a discrete grid, NOT a continuous measurement of a
  tolerance edge.  A finer grid would be needed to locate the edge itself.

And, from `inventory`, on the smallest decision-relevant pool N_min (copies):

  * N_min <  100   -> intrinsic relative noise >~ 10 %; a continuous description
                      of that pool is questionable.
  * N_min >= 1000  -> intrinsic relative noise <~ 3 %; the deterministic
                      description is defensible for that pool.

Both rules are diagnostic.  Neither licenses the sentence "the counter works in
real cells" - that needs (3)-(6).

Usage (NOT run yet)
-------------------
    python preaudit_stochasticity.py plan
    python preaudit_stochasticity.py inventory
    python preaudit_stochasticity.py h2check --hours 100
    python preaudit_stochasticity.py tolerance --shard S --nshards N
    python preaudit_stochasticity.py merge --nshards N

Cost: 3 x 7 = 21 tolerance runs + 1 inventory + 3 h2check runs.  A 600 h run is
~250 CPU-s, so the tolerance scan is ~90 CPU-min total and shards well; the
h2check is registered at 100 h to keep it cheap (it compares implementations, not
behaviour).
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve().parent
# Both directories are needed and the ORDER does not matter as the names are
# disjoint: `plausibility_common` / `working_point` / `scan_carry_pairing` live in
# _HERE, while `model*.py` / `verify_*.py` live one level up.  The rest of the
# package relies on sys.path[0] (the script's own directory) for the first group,
# which only holds when the script is launched from its own folder; inserting both
# explicitly makes this script runnable from anywhere.
sys.path.insert(0, str(_HERE))
sys.path.insert(0, str(_HERE.parent))
from plausibility_common import (OUT, build_threebit, frozen_extension,  # noqa: E402
                                 gate_windows, merge_shards, restore,
                                 shard_slice, sha256, source_hashes,
                                 write_manifest, write_shard)
from model_threebit51 import STATE_NAMES_51  # noqa: E402
from working_point import working_point_block  # noqa: E402

NAME = 'preaudit_stochasticity'
INV_NAME = 'preaudit_stochasticity_inventory'
H2_NAME = 'preaudit_stochasticity_h2check'
ROOT_MODEL = Path(__file__).resolve().parents[1] / 'model_threebit51.py'

# ------------------------------------------------------------ pre-registered
HOURS = 600.0                 # 200 h gives only ~10 steady windows: not certifiable
H2_HOURS = 100.0              # the h2check compares code paths, not behaviour
SAMPLE_MIN = 2.0
MAX_STEP_MIN = 2.0
MIN_STEADY_READS = 16         # the project's minimum for a certification verdict
BURN_IN_H = 100.0             # inventory statistics start after the cold phase

F_GRID = (0.80, 0.90, 0.95, 1.00, 1.05, 1.10, 1.20)
KNOBS = ('conc_scale', 'clock_K', 'K_A1')

# Inventory: concentration-like pools only.  Deliberately EXCLUDES the storage
# variable b*_S (a latched state, not a pool - converting it to copies would be a
# category error) and the immature forms b*_I_u / b*_R_u / A0_u / F0_u / A1_u /
# F1_u (transient intermediates).  The exclusions are reported in the artefact.
INVENTORY_STATES = (
    'm_TetR', 'TetR_total', 'm_CI', 'CI', 'm_LacI', 'LacI',
    'b0_I', 'b0_R', 'b0_C', 'b1_I', 'b1_R', 'b1_C', 'b2_I', 'b2_R', 'b2_C',
    'b0_M_I', 'b0_M_R', 'b1_M_I', 'b1_M_R', 'b2_M_I', 'b2_M_R',
    'A0', 'F0', 'A1', 'F1', 'M_A1', 'M_F1',
)
# Pools whose copy number decides a WRITE (the gate's integrase and its RDF).
DECISION_STATES = ('b2_I', 'b2_R', 'b2_C', 'A1', 'F1')

SMALL_NUMBER_THRESHOLD = 100.0     # copies; below this, >=10 % intrinsic noise
LARGE_NUMBER_THRESHOLD = 1000.0    # copies; above this, <=3 % intrinsic noise
WIDTH_SMALL = 0.05
WIDTH_MODERATE = 0.10


def points():
    return [dict(knob=k, factor=f) for k in KNOBS for f in F_GRID]


def point_id(p):
    return f"{p['knob']}|f={p['factor']:g}"


def copies_per_au(model):
    """The model's OWN a.u. -> copies factor, with the access path recorded.

    model.py:90 defines `copies_per_au = 602.214076 * uM_per_au` for a 1 fL cell.
    ThreeBit51Model wraps TwoBit34Model wraps Model, so the attribute may sit on
    any of the three depending on version: each path is tried, the value is
    type-checked to be a number (not a bound method), and the path used is
    returned so the artefact records provenance instead of assuming it.
    """
    for path in ('copies_per_au', 'base.copies_per_au', 'inner.copies_per_au',
                 'base.base.copies_per_au'):
        obj, found = model, True
        for part in path.split('.'):
            if hasattr(obj, part):
                obj = getattr(obj, part)
            else:
                found = False
                break
        if found and isinstance(obj, (int, float, np.floating)) and obj > 0:
            return float(obj), path
    raise RuntimeError('cannot determine copies_per_au: no known attribute path matched')


def copies_scale_from(copies_per_au_value):
    """Invert model.py:90 (`copies_per_au = 602.214076 * uM_per_au`)."""
    return float(copies_per_au_value) / 602.214076


def build_for(p):
    """Build the frozen model with exactly ONE knob perturbed.

    Returns (model, patched, realised_factor).  `realised_factor` is computed from
    the local parameter objects rather than read back off the model, so it does
    not depend on which wrapper exposes the extension.
    """
    import model as M
    frozen = frozen_extension()
    knob, f = p['knob'], float(p['factor'])
    if knob == 'conc_scale':
        ext = replace(frozen, uM_per_au=frozen.uM_per_au * f)
        model, patched = build_threebit(extension=ext)
        return model, patched, float(ext.uM_per_au) / float(frozen.uM_per_au)
    if knob == 'clock_K':
        ext = replace(frozen, clock_K_au=frozen.clock_K_au * f)
        model, patched = build_threebit(extension=ext)
        return model, patched, float(ext.clock_K_au) / float(frozen.clock_K_au)
    if knob == 'K_A1':
        baseline = float(M.ZENG['K_A'][1])
        target = baseline * f
        model, patched = build_threebit(k_A1=target)
        return model, patched, target / baseline
    raise ValueError(f"unknown knob {knob!r}")


def analyse_run(model, sol, hours):
    """The certification/leak readout shared by every subcommand.

    The frozen verifier's predicate is

        certified = steady and one and order and alt and contrast

    so recording only `certified` is not enough to interpret a failure: four of the
    first run's non-certified points had a PERFECT digital readout (14/14 crossings,
    the exact 12345670 sequence, no unlabelled window) and still failed.  The five
    arms and the contrast ratio are therefore recorded separately, so a failure can
    be attributed instead of guessed.
    """
    from verify_threebit51 import analyse_threebit
    from scan_carry_pairing import analyse as paired_analyse
    a = analyse_threebit(model, sol, hours)
    sig = model.diagnostic_signals(sol.y)
    wins = gate_windows(sol.t, sig['g1'])
    agg, _, _, _, _ = paired_analyse(model, sol, 0.05)
    steady = float(a['steady_state']['reads'])
    cv = a['causal_verdict']
    contrast = a['gate_contrast']
    arms = dict(steady=bool(a['steady_state']['passed']),
                one=bool(cv['exactly_one_gate_and_flip_per_late_reverse']),
                order=bool(cv['causal_order_after_reverse_start']),
                alt=bool(cv['bit2_directions_alternate']),
                contrast=bool(contrast['passed']))
    return dict(
        certified=bool(a['certified']),
        arm_steady=arms['steady'], arm_one=arms['one'], arm_order=arms['order'],
        arm_alt=arms['alt'], arm_contrast=arms['contrast'],
        failed_arms=','.join(k for k, ok in arms.items() if not ok) or '',
        off_on_gate_peak_ratio=contrast['off_peak_to_on_peak'],
        off_on_gate_peak_threshold=contrast['threshold'],
        unassigned_gates=len(cv['unassigned_gate_events']),
        unassigned_crossings=len(cv['unassigned_bit2_crossings']),
        crossings=len(a['bit2_crossings']),
        gates=len(a['carry1_gate_events']),
        reverse_events=len(a['bit1_reverse_events']),
        cold_reads=a['cold_start']['reads'],
        steady_reads=a['steady_state']['reads'],
        sequence=a['steady_state']['sequence'],
        unlabelled=int(sum(r['value'] is None for r in a['read_windows'])),
        gate_peak_median=float(np.median([w['peak'] for w in wins])) if wins else None,
        gate_dose_median=float(np.median([w['dose'] for w in wins])) if wins else None,
        gate_width_median=float(np.median([w['width_h'] for w in wins])) if wins else None,
        n_gate_windows=len(wins),
        setup_min_h=a['bit_margins']['bit2']['min_setup_h'],
        hold_min_h=a['bit_margins']['bit2']['min_hold_h'],
        min_commitment=a['steady_state']['minimum_commitment'],
        leak_symmetric=agg['L_symmetric']['median'],
        leak_pairs=agg['n_pairs'],
        steady_reads_sufficient=bool(steady >= MIN_STEADY_READS),
    )


def evaluate(p, hours, sample_min, max_step_min):
    """One run of one perturbed point.  Strictly serial: one live model at a time."""
    import model as M
    started = time.perf_counter()
    row = dict(point_id=point_id(p), **p)
    baseline_kA1 = float(M.ZENG['K_A'][1])
    patched = None
    try:
        model, patched, realised = build_for(p)
        row['realised_factor'] = realised
        row['n_A1_gate_effective'] = float(model.n_A1_gate_effective)
        sol = model.simulate(hours=hours, sample_min=sample_min,
                             max_step_min=max_step_min)
        row.update(ok=True, error='', **analyse_run(model, sol, hours))
        row['final_state_sum_abs'] = float(np.sum(np.abs(sol.y[:, -1])))
    except Exception as exc:                                     # noqa: BLE001
        row.update(ok=False, error=repr(exc), certified=False)
    finally:
        if patched:
            restore(patched)
        row['zeng_K_A1_restored'] = bool(
            abs(float(M.ZENG['K_A'][1]) - baseline_kA1) == 0.0)
    row['runtime_s'] = round(time.perf_counter() - started, 2)
    return row


def scan(args):
    pts = shard_slice(points(), args.shard, args.nshards)
    tag = f'shard{args.shard:02d}of{args.nshards:02d}'
    rows = []
    for p in pts:
        r = evaluate(p, args.hours, args.sample_min, args.max_step_min)
        print(f"[{tag}] {point_id(p):<22} cert={r.get('certified')} "
              f"cross={r.get('crossings')} gates={r.get('gates')} "
              f"L={r.get('leak_symmetric')} ({r['runtime_s']}s)", flush=True)
        rows.append(r)
    write_shard(NAME, tag, rows, dict(kind='tolerance', shard=args.shard,
                                      nshards=args.nshards, hours=args.hours,
                                      points=len(pts), model_sha256=sha256(ROOT_MODEL),
                                      source_hashes=source_hashes()))
    print(f'[{tag}] wrote {len(rows)} rows')


def _symmetric_certified_width(df, knob):
    """Largest symmetric connected certified factor interval around 1.0.

    Both edges are found by walking OUTWARD from f = 1.0 and stopping at the first
    uncertified grid factor.  That is the definition of "the contiguous certified
    run containing 1.0"; an inward walk from the grid edges is wrong whenever the
    certified set does not touch a grid edge (it would then report lower = 1.0 for
    a knob that certifies at, say, 0.95 as well).

    Returns (width, detail).  width is None when the f = 1.0 row is missing or not
    certified, because then there is no certified point to widen from.
    """
    sub = df[(df.knob == knob) & (df.ok.astype(bool))].sort_values('factor')
    if sub.empty:
        return None, 'no usable row for this knob'
    at_one = sub[sub.factor == 1.0]
    if at_one.empty or not bool(at_one.certified.astype(bool).any()):
        return None, 'f=1.0 is missing or not certified: nothing to widen from'
    fac = {float(r.factor): bool(r.certified) for r in sub.itertuples()}
    lo = 1.0
    for f in sorted((x for x in fac if x < 1.0), reverse=True):
        if fac[f]:
            lo = f
        else:
            break
    hi = 1.0
    for f in sorted(x for x in fac if x > 1.0):
        if fac[f]:
            hi = f
        else:
            break
    width = float(min(1.0 - lo, hi - 1.0))
    return width, dict(lower=lo, upper=hi,
                       certified_factors=[f for f in sorted(fac) if fac[f]],
                       uncertified_factors=[f for f in sorted(fac) if not fac[f]],
                       grid_floor=min(fac) if fac else None,
                       grid_ceiling=max(fac) if fac else None,
                       touches_grid_edge=bool(lo == min(fac) or hi == max(fac)),
                       width_is_grid_limited=bool(lo == min(fac) or hi == max(fac)))


def merge(args):
    df = merge_shards(NAME, args.nshards).sort_values(['knob', 'factor'])
    ok = df[df.ok.astype(bool)]
    widths, details = {}, {}
    for knob in KNOBS:
        w, d = _symmetric_certified_width(df, knob)
        widths[knob] = w
        details[knob] = d
    worst = None
    if all(w is not None for w in widths.values()):
        worst = float(min(widths.values()))
    if worst is None:
        size_verdict = 'UNDECIDED: f=1.0 is not certified for at least one knob'
        worst_knob = None
    elif worst < WIDTH_SMALL:
        worst_knob = min((k for k in widths if widths[k] == worst))
        size_verdict = (f'SMALL: symmetric static tolerance = {worst:g} (< '
                        f'{WIDTH_SMALL:g}). A stochastic study is REQUIRED before '
                        f'any counting claim.')
    elif worst < WIDTH_MODERATE:
        worst_knob = min((k for k in widths if widths[k] == worst))
        size_verdict = (f'MODERATE: symmetric static tolerance = {worst:g} '
                        f'({WIDTH_SMALL:g}..{WIDTH_MODERATE:g}). Report the '
                        f'interval; do not claim noise robustness.')
    else:
        worst_knob = min((k for k in widths if widths[k] == worst))
        size_verdict = (f'WIDE: symmetric static tolerance = {worst:g} (>= '
                        f'{WIDTH_MODERATE:g}). The deterministic conclusion '
                        f'tolerates slow perturbations of at least this size; a '
                        f'stochastic study stays DESIRABLE for fast noise, phase '
                        f'diffusion, plasmid loss and single-cell rates.')
    # The f=1.0 rows must at least agree with each other; bit-identity is settled
    # by the h2check subcommand, which runs all three paths in ONE process.
    f1 = df[(df.factor == 1.0) & df.ok.astype(bool)].sort_values('knob')
    f1_agreement = dict(n=int(len(f1)))
    if len(f1):
        f1_agreement.update(
            sequences_all_equal=bool(f1.sequence.nunique() == 1),
            crossings_all_equal=bool(f1.crossings.nunique() == 1),
            gates_all_equal=bool(f1.gates.nunique() == 1),
            leak_spread=float(f1.leak_symmetric.max() - f1.leak_symmetric.min()))
    noncert = ok[~ok.certified.astype(bool)]
    # Attribute every failure to the arm(s) of
    #   certified = steady and one and order and alt and contrast
    # that failed.  Without this, "certified=False" is uninterpretable: the first
    # run had four points with a perfect digital readout that failed anyway.
    arm_cols = ('arm_steady', 'arm_one', 'arm_order', 'arm_alt', 'arm_contrast')
    failure_arms = {c.replace('arm_', ''): int((~noncert[c].astype(bool)).sum())
                    for c in arm_cols if c in noncert.columns}
    contrast_fail_only = int(sum(
        1 for r in noncert.itertuples()
        if not bool(r.arm_contrast)
        and all(bool(getattr(r, c)) for c in arm_cols if c != 'arm_contrast')))
    digital_readout_perfect_but_uncertified = int(sum(
        1 for r in noncert.itertuples()
        if int(r.crossings) == int(r.gates) == 14 and int(r.unlabelled) == 0
        and bool(r.arm_steady) and bool(r.arm_one)))
    verdict = dict(
        hours=args.hours, grid=dict(factors=list(F_GRID), knobs=list(KNOBS),
                                    note=('W can only take 0.00/0.05/0.10/0.20 on '
                                          'this grid: a classification, not a '
                                          'measured edge')),
        points=int(len(df)), ok=int(len(ok)),
        certified=int(ok.certified.astype(bool).sum()),
        min_steady_reads=(int(pd.to_numeric(ok.steady_reads, errors='coerce').min())
                          if len(ok) else None),
        steady_reads_sufficient_for_certification=bool(
            len(ok) and pd.to_numeric(ok.steady_reads, errors='coerce').min()
            >= MIN_STEADY_READS),
        symmetric_certified_width=widths,
        symmetric_certified_interval=details,
        worst_knob=worst_knob, worst_width=worst, size_verdict=size_verdict,
        certification_arms=dict(
            predicate='certified = steady and one and order and alt and contrast',
            n_non_certified=int(len(noncert)),
            failures_by_arm=failure_arms,
            contrast_only_failures=contrast_fail_only,
            digital_readout_perfect_but_uncertified=digital_readout_perfect_but_uncertified,
            note=('The contrast arm is the frozen verifier\'s off/on GATE PEAK ratio '
                  'criterion (verify_threebit51.MAX_OFF_ON_GATE_RATIO = 0.10). H4 '
                  'demoted exactly this off/on peak-ratio metric for LEAKAGE '
                  'quantification, but it is still part of the certification '
                  'predicate, so a point can count perfectly and still fail '
                  'certification through that arm.')),
        off_on_gate_peak_ratio_table=ok.pivot_table(
            index='knob', columns='factor', values='off_on_gate_peak_ratio').to_dict(),
        f1_cross_knob_agreement=f1_agreement,
        f1_bit_identity_checked_by=H2_NAME,
        non_certified_points=[
            dict(point_id=r.point_id, factor=r.factor, crossings=r.crossings,
                 gates=r.gates, setup_min_h=r.setup_min_h, hold_min_h=r.hold_min_h,
                 unlabelled=r.unlabelled, steady_reads=r.steady_reads,
                 min_commitment=r.min_commitment, sequence=r.sequence,
                 leak_symmetric=r.leak_symmetric,
                 failed_arms=r.failed_arms,
                 off_on_gate_peak_ratio=r.off_on_gate_peak_ratio,
                 arm_steady=bool(r.arm_steady), arm_one=bool(r.arm_one),
                 arm_order=bool(r.arm_order), arm_alt=bool(r.arm_alt),
                 arm_contrast=bool(r.arm_contrast))
            for r in noncert.itertuples()],
        certification_table=ok.pivot_table(index='knob', columns='factor',
                                           values='certified').to_dict(),
        leak_table=ok.pivot_table(index='knob', columns='factor',
                                  values='leak_symmetric').round(6).to_dict(),
        reading=('A STATIC perturbation tolerance is NOT a stochastic result. It '
                 'bounds only fluctuations slow compared with one carry period '
                 '(42.4 h). Fast noise, clock phase diffusion, plasmid segregation '
                 'and single-cell success rates are not measured here and must not '
                 'be inferred - see (3)-(6) in the module docstring.'),
        model_sha256=sha256(ROOT_MODEL), source_sha256=source_hashes(),
        working_point=working_point_block(None))
    df.to_csv(OUT / f'{NAME}_all.csv', index=False, encoding='utf-8')
    (OUT / f'{NAME}_verdict.json').write_text(
        json.dumps(verdict, ensure_ascii=False, indent=2), encoding='utf-8')
    write_manifest()
    pd.set_option('display.width', 200)
    print('--- certified (rows=knob, cols=factor) ---')
    print(ok.pivot_table(index='knob', columns='factor',
                         values='certified').to_string())
    print('\n--- paired leak L_symmetric ---')
    print(ok.pivot_table(index='knob', columns='factor',
                         values='leak_symmetric').round(5).to_string())
    print('\n--- symmetric certified width per knob ---')
    print(json.dumps(widths, indent=2))
    print('\n--- f=1.0 cross-knob agreement ---')
    print(json.dumps(f1_agreement, indent=2))
    print('\n--- size verdict ---')
    print(size_verdict)
    if len(noncert):
        print(f'\n--- {len(noncert)} non-certified point(s) ---')
        for r in noncert.itertuples():
            print(f'  {r.point_id:<22} cross={r.crossings:>2} gates={r.gates:>2} '
                  f'setup={r.setup_min_h:.4g} hold={r.hold_min_h:.4g} '
                  f'off/on={r.off_on_gate_peak_ratio} failed=[{r.failed_arms}]')
        print('\n--- failures by arm ---')
        print(json.dumps(failure_arms, indent=2))
        print(f'contrast-only failures: {contrast_fail_only}')
        print(f'perfect digital readout but uncertified: '
              f'{digital_readout_perfect_but_uncertified}')


def h2check(args):
    """All three f = 1.0 construction paths, compared as full arrays, one process."""
    import model as M
    frozen = frozen_extension()
    baseline_kA1 = float(M.ZENG['K_A'][1])
    runs, problems = {}, []
    order = ('conc_scale', 'clock_K', 'K_A1')
    for knob in order:
        patched = None
        try:
            model, patched, realised = build_for(dict(knob=knob, factor=1.0))
            sol = model.simulate(hours=args.hours, sample_min=args.sample_min,
                                 max_step_min=args.max_step_min)
            runs[knob] = dict(t=np.array(sol.t, copy=True),
                              y=np.array(sol.y, copy=True),
                              realised_factor=realised)
            del model
        except Exception as exc:                                 # noqa: BLE001
            problems.append(f'{knob}: {exc!r}')
        finally:
            if patched:
                restore(patched)
            if abs(float(M.ZENG['K_A'][1]) - baseline_kA1) != 0.0:
                problems.append(f'{knob}: ZENG["K_A"][1] was not restored')
    comparison = {}
    ref = order[0]
    if ref in runs:
        for knob in order[1:]:
            if knob not in runs:
                continue
            same_t = bool(np.array_equal(runs[ref]['t'], runs[knob]['t']))
            if same_t:
                gap = float(np.max(np.abs(runs[ref]['y'] - runs[knob]['y'])))
            else:
                gap = None
                problems.append(f'{knob}: time grid differs from {ref}')
            comparison[f'{ref}_vs_{knob}'] = dict(
                time_grid_identical=same_t, max_abs_state_gap=gap,
                bit_identical=bool(gap == 0.0))
    if not comparison:
        problems.append('fewer than two construction paths produced a trajectory')
    elif not all(v['bit_identical'] for v in comparison.values()):
        problems.append('the three f=1.0 construction paths are NOT bit-identical: '
                        'a ZENG patch leaked between models (H2) - the tolerance '
                        'scan would be void')
    artefact = dict(
        kind='h2_cross_path_bit_identity', hours=args.hours,
        paths=list(order), realised_factors={k: v['realised_factor'] for k, v in runs.items()},
        comparison=comparison, problems=problems, passed=bool(not problems),
        reading=('The three ways of stating the frozen working point (extension with '
                 'the frozen uM_per_au, extension with the frozen clock_K_au, and a '
                 'ZENG K_A patch at its own frozen value) must produce the SAME model, '
                 'hence bit-identical trajectories. Any non-zero gap means a global '
                 'ZENG patch contaminated another model instance (H2).'),
        model_sha256=sha256(ROOT_MODEL), source_sha256=source_hashes(),
        working_point=working_point_block(None))
    (OUT / f'{H2_NAME}.json').write_text(
        json.dumps(artefact, ensure_ascii=False, indent=2), encoding='utf-8')
    write_manifest()
    print(json.dumps(artefact, ensure_ascii=False, indent=2))
    return artefact


def inventory(args):
    """One frozen run: copy numbers of every concentration-like pool."""
    row, patched = None, None
    started = time.perf_counter()
    try:
        model, patched, _ = build_for(dict(knob='conc_scale', factor=1.0))
        cpa, cpa_path = copies_per_au(model)
        sol = model.simulate(hours=args.hours, sample_min=args.sample_min,
                             max_step_min=args.max_step_min)
        a = analyse_run(model, sol, args.hours)
        sig = model.diagnostic_signals(sol.y)
        wins = gate_windows(sol.t, sig['g1'])
        post = np.asarray(sol.t) >= BURN_IN_H
        in_gate = np.zeros(len(sol.t), dtype=bool)
        for w in wins:
            in_gate[w['i0']:w['i1'] + 1] = True
        decision = in_gate & post
        per_state = {}
        for name in INVENTORY_STATES:
            if name not in STATE_NAMES_51:
                raise RuntimeError(f'{name} is not a state of the 51-state model')
            y = np.asarray(sol.y[STATE_NAMES_51.index(name)], dtype=float) * cpa
            per_state[name] = dict(
                copies_min_post_burnin=float(y[post].min()),
                copies_p05_post_burnin=float(np.percentile(y[post], 5)),
                copies_median_post_burnin=float(np.median(y[post])),
                copies_max_post_burnin=float(y[post].max()),
                copies_min_in_gate=(float(y[decision].min())
                                    if decision.any() else None),
                copies_median_in_gate=(float(np.median(y[decision]))
                                       if decision.any() else None),
                copies_p05_in_gate=(float(np.percentile(y[decision], 5))
                                    if decision.any() else None),
                is_decision_pool=bool(name in DECISION_STATES))

        # ------------------------------------------------------------------
        # THREE rules, all kept, because the first two both turned out to be
        # mis-specified and the revisions must stay visible.
        #
        # Rule A (pre-registered in the module docstring): thresholds applied to the
        # smallest MINIMUM over the decision pools.  DEGENERATE.  It fired on b2_C at
        # 5.1e-9 copies - a pool that is structurally ABSENT at that instant.  "Absent"
        # is not a small-number regime, so Rule A named the wrong pool.
        #
        # Rule B (added after run 1): thresholds applied to the smallest pool whose
        # whole-cycle MEDIAN is >= 1 copy.  ALSO DEGENERATE, in a second way: b2_I is
        # zero for most of the cycle, so its median (1.97) describes the OFF phase, not
        # the level at which it does its job; and the per-pool uM_per_au sensitivity
        # still folded the whole-cycle minimum (0.0011 copies) back in, yielding a
        # nonsense requirement of 5.4e5 x the frozen conversion.
        #
        # Two failed revisions in one session is evidence that the AGGREGATE STATISTIC
        # is the wrong instrument, not that the threshold needs tuning.  Every pool
        # here is pulsatile, so no single number ("the smallest pool") is defensible.
        # Rule C therefore REFUSES TO NAME ONE POOL and reports per-pool, per-phase
        # numbers, with the noise implication (CV ~ 1/sqrt(N)) at each.  It is the
        # reported rule; A and B are retained only as the record of the revisions.
        # ------------------------------------------------------------------
        dec = {k: v for k, v in per_state.items() if v['is_decision_pool']}

        def _verdict(n_min, pool, tag):
            if n_min < SMALL_NUMBER_THRESHOLD:
                return (f'{tag}: the smallest decision-relevant pool ({pool}) reaches '
                        f'{n_min:.3g} copies, below {SMALL_NUMBER_THRESHOLD:g}. '
                        f'Intrinsic relative noise is >~ 10 % and a continuous '
                        f'description of that pool is questionable.')
            if n_min >= LARGE_NUMBER_THRESHOLD:
                return (f'{tag}: the smallest decision-relevant pool ({pool}) stays at '
                        f'{n_min:.3g} copies, above {LARGE_NUMBER_THRESHOLD:g}. '
                        f'Intrinsic relative noise is <~ 3 % and the deterministic '
                        f'description is defensible for that pool.')
            return (f'{tag}: the smallest decision-relevant pool ({pool}) reaches '
                    f'{n_min:.3g} copies, i.e. 10-30 % intrinsic noise. Report the '
                    f'number; do not claim either way.')

        smallest_a = min(dec, key=lambda k: dec[k]['copies_min_post_burnin'])
        n_min_a = float(dec[smallest_a]['copies_min_post_burnin'])
        present_by_median = {k: v for k, v in dec.items()
                             if v['copies_median_post_burnin'] >= 1.0}

        # Rule C: per-pool, per-phase, no aggregate label.
        uM = float(getattr(args, 'uM_per_au', None)
                   or copies_scale_from(cpa))
        noise_regime = {}
        for k, v in dec.items():
            entries = {}
            phases = {
                'gate_phase_minimum': v['copies_min_in_gate'],
                'gate_phase_p05': v['copies_p05_in_gate'],
                'gate_phase_median': v['copies_median_in_gate'],
                'whole_cycle_trough': v['copies_min_post_burnin'],
                'whole_cycle_median': v['copies_median_post_burnin'],
            }
            for label, val in phases.items():
                if val is None:
                    continue
                if val < 1.0:
                    entries[label] = dict(
                        copies=float(val), regime='STRUCTURALLY ABSENT',
                        implied_cv_pct=None,
                        note=('below one copy: this phase is a zero, not a '
                              'noise-limited level, so it carries no CV'))
                else:
                    entries[label] = dict(
                        copies=float(val),
                        implied_cv_pct=float(100.0 / np.sqrt(val)),
                        uM_per_au_for_100_copies=float(uM * 100.0 / val),
                        factor_vs_frozen=float(100.0 / val),
                        regime=('small numbers' if val < SMALL_NUMBER_THRESHOLD
                                else ('intermediate' if val < LARGE_NUMBER_THRESHOLD
                                      else 'large numbers')))
            noise_regime[k] = entries
        # the binding numbers, named explicitly rather than min()-ed across phases
        binding = []
        for k, entries in noise_regime.items():
            for label, e in entries.items():
                if e.get('copies') and 1.0 <= e['copies'] < SMALL_NUMBER_THRESHOLD:
                    binding.append(dict(pool=k, statistic=label,
                                        copies=e['copies'],
                                        implied_cv_pct=e['implied_cv_pct']))
        binding.sort(key=lambda d: d['copies'])
        # Build the three verdict strings in plain statements.  The nested-call
        # version of this was 60 characters of parentheses deep and shipped with an
        # unbalanced bracket; plain statements cannot hide that.
        rule_a_text = _verdict(
            n_min_a, smallest_a, 'RULE A (degenerate, kept for the record)')
        if present_by_median:
            b_pool = min(present_by_median,
                         key=lambda k: present_by_median[k]['copies_median_post_burnin'])
            b_val = float(present_by_median[b_pool]['copies_median_post_burnin'])
            rule_b_text = _verdict(b_val, b_pool, 'RULE B (superseded)')
        else:
            rule_b_text = ('RULE B (superseded): no decision pool has a whole-cycle '
                           'median >= 1 copy')
        if binding:
            binding_text = '; '.join(
                f"{b['pool']} {b['statistic']} = {b['copies']:.3g} copies "
                f"(CV ~ {b['implied_cv_pct']:.0f} %)" for b in binding[:4])
        else:
            binding_text = 'none below 100 copies'
        row = dict(
            ok=True, error='', hours=args.hours, burn_in_h=BURN_IN_H,
            copies_per_au=cpa, copies_per_au_source=cpa_path,
            cell_volume_fl=1.0, per_state=per_state,
            decision_pools=list(DECISION_STATES),
            number_verdict_rule_a_pre_registered=rule_a_text,
            number_verdict_rule_b_superseded=rule_b_text,
            number_verdict=(
                'NO SINGLE-POOL VERDICT. Every decision pool is pulsatile, so "the '
                'smallest pool" is not a well-defined instrument: two revisions of '
                'that statistic both misfired (see rule_revision_reason). The '
                'binding numbers are the per-pool, per-phase entries in '
                'noise_regime_per_pool; the smallest ones are: '
                + binding_text + '.'),
            noise_regime_per_pool=noise_regime,
            binding_small_number_entries=binding,
            rule_b_is_a_post_hoc_revision=True,
            rule_revision_reason=(
                'Rule A (as pre-registered) fired on b2_C at 5.1e-9 copies - a pool '
                'structurally ABSENT at that instant, not noise-limited. Rule B (added '
                'after run 1) used the whole-cycle median and then folded the '
                'whole-cycle minimum back into its uM_per_au sensitivity, producing a '
                'requirement of 5.4e5x the frozen conversion for b2_I, which is '
                'nonsense because that minimum is again a zero. Rule C reports '
                'per-pool, per-phase numbers and refuses to name a single pool. Only '
                'Rule C is offered as a result; A and B are kept as the record.'),
            uM_per_au_used_for_the_conversion=uM,
            uM_per_au_sensitivity_note=(
                'copy numbers scale LINEARLY with uM_per_au, which is an uncalibrated '
                'interface assumption, so every copy number here inherits it. '
                'noise_regime_per_pool lists, per statistic, the uM_per_au that would '
                'be needed to lift that statistic to 100 copies.'),
            excluded_states=dict(
                storage=['b0_S', 'b1_S', 'b2_S'],
                storage_reason=('the storage variable is a latched state, not a '
                                'pool; converting it to copies would be a category '
                                'error'),
                immature=['b0_I_u', 'b0_R_u', 'b1_I_u', 'b1_R_u', 'b2_I_u',
                          'b2_R_u', 'A0_u', 'F0_u', 'A1_u', 'F1_u'],
                immature_reason='transient intermediates, not the functional pool'),
            certification=dict(certified=a['certified'], crossings=a['crossings'],
                               gates=a['gates'], steady_reads=a['steady_reads'],
                               sequence=a['sequence'],
                               steady_reads_sufficient=a['steady_reads_sufficient']),
            runtime_s=round(time.perf_counter() - started, 2),
            model_sha256=sha256(ROOT_MODEL), source_sha256=source_hashes(),
            working_point=None)
        row['working_point'] = working_point_block(model)
        del model
    except Exception as exc:                                     # noqa: BLE001
        row = dict(ok=False, error=repr(exc), hours=args.hours,
                   runtime_s=round(time.perf_counter() - started, 2))
    finally:
        if patched:
            restore(patched)
    (OUT / f'{INV_NAME}.json').write_text(
        json.dumps(row, ensure_ascii=False, indent=2), encoding='utf-8')
    write_manifest()
    print(json.dumps(row, ensure_ascii=False, indent=2))
    return row


def plan(_args):
    """Print the registered grid, cost and rules without simulating."""
    n = len(points())
    print('pre-registered plan (this subcommand runs no simulation)')
    print(f'  tolerance hours/run : {HOURS}')
    print(f'  h2check hours/run   : {H2_HOURS}')
    print(f'  factor grid         : {list(F_GRID)}')
    print(f'  knobs               : {list(KNOBS)}')
    print(f'  points              : {n} = {len(KNOBS)} knobs x {len(F_GRID)} factors')
    print(f'  estimated cost      : {n} x ~250 CPU-s ~= {n * 250 / 60:.0f} CPU-min, '
          f'+ 3 x {H2_HOURS:g} h for the h2check, + 1 run for inventory')
    print(f'  decision rule       : W < {WIDTH_SMALL} SMALL / < {WIDTH_MODERATE} '
          f'MODERATE / >= {WIDTH_MODERATE} WIDE')
    print(f'  reachable W on grid : 0.00 / 0.05 / 0.10 / 0.20 only')
    print(f'  inventory rule      : < {SMALL_NUMBER_THRESHOLD:g} copies small, '
          f'>= {LARGE_NUMBER_THRESHOLD:g} copies large')
    print('  NOT answered here   : fast noise, clock phase diffusion, plasmid loss, '
          'single-cell rates')
    for fn in (f'{NAME}_all.csv', f'{NAME}_verdict.json', f'{INV_NAME}.json',
               f'{H2_NAME}.json'):
        print(f'  output              : {OUT / fn}')


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest='cmd', required=True)
    sub.add_parser('plan', help='print the registered plan; runs nothing')
    for name, helptext in (('inventory', 'one frozen run: copy-number inventory'),
                           ('h2check', 'three f=1.0 paths compared bit-for-bit')):
        p = sub.add_parser(name, help=helptext)
        p.add_argument('--hours', type=float,
                       default=HOURS if name == 'inventory' else H2_HOURS)
        p.add_argument('--sample-min', type=float, default=SAMPLE_MIN)
        p.add_argument('--max-step-min', type=float, default=MAX_STEP_MIN)
    t = sub.add_parser('tolerance', help='sharded static-tolerance scan')
    t.add_argument('--shard', type=int, required=True)
    t.add_argument('--nshards', type=int, required=True)
    t.add_argument('--hours', type=float, default=HOURS)
    t.add_argument('--sample-min', type=float, default=SAMPLE_MIN)
    t.add_argument('--max-step-min', type=float, default=MAX_STEP_MIN)
    m = sub.add_parser('merge', help='merge the tolerance shards and judge')
    m.add_argument('--nshards', type=int, required=True)
    m.add_argument('--hours', type=float, default=HOURS)
    args = ap.parse_args()
    if args.cmd == 'plan':
        plan(args)
    elif args.cmd == 'inventory':
        inventory(args)
    elif args.cmd == 'h2check':
        h2check(args)
    elif args.cmd == 'tolerance':
        scan(args)
    else:
        merge(args)


if __name__ == '__main__':
    main()
