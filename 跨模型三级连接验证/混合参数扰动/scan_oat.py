"""OAT parameter perturbation scan on the frozen HZH three-bit hybrid cascade.

Scope and discipline
--------------------
* Read-only with respect to every existing model/verifier/certification file.
  The hybrid is built by swapping the receiver's HbyConfig (and, for threshold
  parameters, the extracted ZENG table dict) on a shallow copy of a prototype
  model. Nothing under hby_zmh_hby/ is written.
* Criterion is the hybrid's own HZH_MOD8_CAUSAL_V1 (verify_hzh.analyse); the
  51-state predicate is neither used nor inherited.
* Output per run: the criterion boolean AND the continuous minimum timing
  margin, because a boolean alone cannot distinguish "just passed" from
  "comfortably passed".
* Every run records the fully resolved config read back from the receiver, so
  a silent fallback to a default (the failure mode that invalidated the earlier
  n=4 runs) is caught rather than trusted.

Two pre-registered predictions are asserted, not merely observed:
  P1  clock_scale and clock_K_au are EXACTLY degenerate, because
      hill(s*x, K, n) == hill(x, K/s, n) algebraically.
      => clock_scale=0.5 and clock_K_au=0.6 must give identical trajectories.
  P2  translation_h is structurally unobservable (the translation constant
      cancels from the mature-protein response), so scaling it must leave
      every state except mRNA unchanged.
A failure of either assertion means the harness is wrong, not the model.
"""
from __future__ import annotations

import argparse
import copy
import json
import platform
import sys
import time
from dataclasses import asdict
from datetime import datetime
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parent
PARENT = ROOT.parent
for _p in (PARENT, PARENT / 'hby_zmh_hby', PARENT / 'hby_zmh_hby' / 'certification'):
    sys.path.insert(0, str(_p))

from run_han_comparison import HanInput, HAN_PATH, HAN_SHA          # noqa: E402
from hybrid_model import (HbyConfig, HbyReceiver, sha256,           # noqa: E402
                          source_hashes, require_hash)
from model_hzh import HZHModel, NAMES, IDX                          # noqa: E402
from verify_hzh import analyse                                      # noqa: E402

HOURS = 300.0
SAMPLE_MIN = 1.0
MAX_STEP_MIN = 1.0
RTOL = 2e-7
ATOL = 2e-9

FROZEN = HbyConfig()
ZENG = HbyReceiver().p

# Frozen certified reference (certification/results/20260928_141723_738780).
REF = dict(certified=True, steady_reads=19, clock_peak_count=28,
           sequence='1234567012345670123',
           margin_h=1.1784609018563117,
           stage0=dict(reverse=14, gate=14, flip=14),
           stage1=dict(reverse=7, gate=7, flip=7))


def dump(path, obj):
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False),
                    encoding='utf-8')


# --------------------------------------------------------------------------
# harness
# --------------------------------------------------------------------------
def build_receiver(cfg, zeng_patch=None):
    """Receiver with the given config, plus optional ZENG-table overrides.

    k_complex is re-derived after any patch because HbyReceiver.__init__
    computes it from p['K_D_comp'].
    """
    rec = HbyReceiver(cfg)
    if zeng_patch:
        rec.p = dict(rec.p)
        for (key, slot), value in zeng_patch.items():
            current = rec.p[key]
            if slot is None:
                rec.p[key] = float(value)
            else:
                seq = list(current)
                seq[slot] = float(value)
                rec.p[key] = tuple(seq)
        rec.k_complex = rec.q * rec.p['K_D_comp']
    return rec


def variant(prototype, cfg, zeng_patch=None):
    """Shallow copy sharing h31/han/z, with a fresh receiver.

    Equivalent to constructing HZHModel('han', han) with that config: h31 and
    the ZMH donor table do not depend on HbyConfig, and copies_per_au is
    re-derived exactly as __init__ does.
    """
    m = copy.copy(prototype)
    m.tail = build_receiver(cfg, zeng_patch)
    m.p = m.tail.p
    m.copies_per_au = 602.214076 * m.tail.c.receiver_uM_per_au
    return m


def integrate(model, hours=HOURS, sample_min=SAMPLE_MIN,
              max_step_min=MAX_STEP_MIN, rtol=RTOL, atol=ATOL):
    count = round(hours * 60 / sample_min)
    t = np.linspace(0, hours, count + 1)
    if len(t) < 2 or not np.all(np.diff(t) > 0):
        raise ValueError('Bad sampling grid')
    sol = solve_ivp(model.rhs, (0, hours), model.initial_state(), t_eval=t,
                    method='DOP853', rtol=rtol, atol=atol,
                    max_step=max_step_min / 60)
    if not sol.success or not np.all(np.isfinite(sol.y)):
        raise RuntimeError(sol.message)
    return sol.t, sol.y


def work_point_readback(model):
    c = model.tail.c
    return dict(asdict(c), n_A1_gate_effective=float(c.n_A1_gate),
                clock_K_effective=float(c.clock_K_au / c.clock_scale),
                k_complex=float(model.tail.k_complex),
                K_A1=float(model.p['K_A'][1]), K_F1=float(model.p['K_F'][1]),
                K_D_int0=float(model.p['K_D_int'][0]),
                K_inh=float(model.p['K_inh']),
                K_D_comp=float(model.p['K_D_comp']),
                copies_per_au=float(model.copies_per_au))


def evaluate(prototype, label, cfg, zeng_patch=None, ref_traj=None):
    """One run.

    ref_traj is the baseline trajectory. The max|dy| gap is reported alongside
    the verdict because the global minimum timing margin is pinned by bit0
    (1.178460901856866 h, a hold margin upstream of the clock gate) and is
    therefore blind to every downstream parameter. A downstream perturbation
    must be judged on the per-bit margins, not on the global minimum.
    """
    row = dict(label=label)
    tic = time.perf_counter()
    try:
        m = variant(prototype, cfg, zeng_patch)
        row['work_point'] = work_point_readback(m)
        t, y = integrate(m)
        v, sig = analyse(t, y, m.z, m.tail)
        bm = v['bit_margins']
        row.update(
            certified=bool(v['certified_v1']),
            steady_reads=int(v['steady']['reads']),
            cold_reads=int(v['cold']['reads']),
            increments_mod8=bool(v['steady']['increments_mod8']),
            sequence=str(v['steady']['sequence']),
            cold_sequence=str(v['cold']['sequence']),
            minimum_commitment=(None if v['steady']['minimum_commitment'] is None
                                else float(v['steady']['minimum_commitment'])),
            boundary_clips=int(v['steady']['boundary_clips']),
            margin_h=(None if v['global_min_timing_margin_h'] is None
                      else float(v['global_min_timing_margin_h'])),
            clock_peak_count=int(v['clock_peak_count']),
            bit_margins={b: {k: (None if x is None else float(x))
                             for k, x in mm.items()}
                         for b, mm in bm.items()},
            stage0=dict(reverse=int(v['events']['bit0_to_bit1']['reverse_events']),
                        gate=int(v['events']['bit0_to_bit1']['gate_events']),
                        flip=int(v['events']['bit0_to_bit1']['flip_events']),
                        one_to_one=bool(v['events']['bit0_to_bit1']['one_to_one']),
                        order=bool(v['events']['bit0_to_bit1']['causal_order']),
                        alt=bool(v['events']['bit0_to_bit1']['alternating_directions'])),
            stage1=dict(reverse=int(v['events']['bit1_to_bit2']['reverse_events']),
                        gate=int(v['events']['bit1_to_bit2']['gate_events']),
                        flip=int(v['events']['bit1_to_bit2']['flip_events']),
                        one_to_one=bool(v['events']['bit1_to_bit2']['one_to_one']),
                        order=bool(v['events']['bit1_to_bit2']['causal_order']),
                        alt=bool(v['events']['bit1_to_bit2']['alternating_directions'])),
            error=None)
        if ref_traj is not None:
            row['max_abs_gap_vs_baseline'] = float(
                np.max(np.abs(y - ref_traj)))
            row['max_abs_gap_int0'] = float(np.max(np.abs(y[2] - ref_traj[2])))
        else:
            row['max_abs_gap_vs_baseline'] = 0.0
            row['max_abs_gap_int0'] = 0.0
        row['_traj'] = y
    except Exception as exc:                                   # keep failures
        row.update(certified=False, error=repr(exc), steady_reads=0,
                   sequence='', cold_sequence='', increments_mod8=False,
                   minimum_commitment=None, boundary_clips=0, margin_h=None,
                   clock_peak_count=0, stage0=None, stage1=None,
                   bit_margins=None, max_abs_gap_vs_baseline=None,
                   max_abs_gap_int0=None)
    row['runtime_s'] = time.perf_counter() - tic
    return row


def min_of(mm, keys=('min_setup_h', 'min_hold_h')):
    if not mm:
        return None
    vals = [mm[k] for k in keys if mm.get(k) is not None]
    return min(vals) if vals else None


# --------------------------------------------------------------------------
# perturbation specification
# --------------------------------------------------------------------------
def mult(x):
    return [('x%g' % f, FROZEN_LOOKUP(x) * f) for f in (0.5, 0.8, 1.25, 2.0)]


def FROZEN_LOOKUP(name):
    return getattr(FROZEN, name)


def spec_rows():
    """(label, tier, kind, target, levels)  kind in {config, zeng}."""
    out = []

    def cfg(label, tier, name, levels):
        for tag, val in levels:
            out.append((f'{label}={tag}', tier, 'config', name, val))

    def zeng(label, tier, key, slot, base, factors):
        for f in factors:
            out.append((f'{label}=x{f:g}', tier, 'zeng', (key, slot), base * f))

    # --- Tier A: interface parameters the hybrid uniquely owns -------------
    cfg('clock_K_au', 'A', 'clock_K_au', mult('clock_K_au'))
    cfg('clock_scale', 'A', 'clock_scale', mult('clock_scale'))
    cfg('receiver_uM_per_au', 'A', 'receiver_uM_per_au', mult('receiver_uM_per_au'))
    cfg('clock_n', 'A', 'clock_n', [('1', 1.0), ('3', 3.0), ('4', 4.0)])

    # --- Tier B: gate response (the D-part's stated first priority) -------
    cfg('n_A1_gate', 'B', 'n_A1_gate', [('4', 4.0), ('5', 5.0), ('7', 7.0), ('8', 8.0)])
    zeng('K_A[1]', 'B', 'K_A', 1, ZENG['K_A'][1], (0.5, 0.8, 1.25, 2.0))
    zeng('K_F[1]', 'B', 'K_F', 1, ZENG['K_F'][1], (0.5, 0.8, 1.25, 2.0))

    # --- Tier C: expression timing (the D-part's second priority) ---------
    cfg('carry_mrna_min', 'C', 'carry_mrna_min', mult('carry_mrna_min'))
    cfg('carry_maturation_min', 'C', 'carry_maturation_min', mult('carry_maturation_min'))
    cfg('bit_mrna_min', 'C', 'bit_mrna_min', mult('bit_mrna_min'))
    cfg('bit_maturation_min', 'C', 'bit_maturation_min', mult('bit_maturation_min'))

    # --- explicit complex kinetics (NOT reducible to q in an explicit ODE) -
    cfg('kon', 'D', 'kon', [('x2', 0.2)])
    cfg('complex_decay', 'D', 'complex_decay', [('x2', 2.0)])

    # --- control: predicted exactly null ---------------------------------
    cfg('translation_h', 'CTRL', 'translation_h', [('x0.5', 15.0), ('x2', 60.0)])

    return out


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--runs', type=int, default=0, help='limit number of OAT rows')
    ap.add_argument('--no-controls', action='store_true')
    args = ap.parse_args()
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8')

    out = ROOT / 'results' / datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    out.mkdir(parents=True, exist_ok=False)
    dump(out / 'status.json', dict(status='RUNNING'))

    sources = {str(p): sha256(p) for p in
               (Path(__file__), PARENT / 'hybrid_model.py',
                PARENT / 'hby_zmh_hby' / 'model_hzh.py',
                PARENT / 'hby_zmh_hby' / 'certification' / 'verify_hzh.py',
                PARENT / 'run_han_comparison.py', HAN_PATH)}
    sources.update(source_hashes())

    print('building Han upstream once (config-independent, reused for all runs)',
          flush=True)
    tic = time.perf_counter()
    han = HanInput(HOURS)
    han_s = time.perf_counter() - tic
    assert han.p.peak_load_fraction == 0.0, 'upstream load feedback appeared'
    tic = time.perf_counter()
    prototype = HZHModel('han', han)
    proto_s = time.perf_counter() - tic
    print(f'  HanInput {han_s:.2f} s ; prototype (incl. h31) {proto_s:.2f} s',
          flush=True)

    rows = []

    # ---- baseline guard: must reproduce the frozen certification ----------
    base = evaluate(prototype, 'BASELINE', FROZEN)
    b = base
    checks = {
        'certified': b['certified'] == REF['certified'],
        'steady_reads': b['steady_reads'] == REF['steady_reads'],
        'clock_peak_count': b['clock_peak_count'] == REF['clock_peak_count'],
        'sequence': b['sequence'] == REF['sequence'],
        'increments_mod8': b['increments_mod8'] is True,
        'stage0_events': (b['stage0']['reverse'], b['stage0']['gate'],
                          b['stage0']['flip']) == (14, 14, 14),
        'stage1_events': (b['stage1']['reverse'], b['stage1']['gate'],
                          b['stage1']['flip']) == (7, 7, 7),
        'margin_matches_frozen': abs(b['margin_h'] - REF['margin_h']) < 1e-6,
        'n_A1_gate_readback': b['work_point']['n_A1_gate_effective'] == 6.0,
        'clock_K_readback': abs(b['work_point']['clock_K_effective'] - 0.3) < 1e-12,
    }
    print('baseline guard:', json.dumps(checks), flush=True)
    if not all(checks.values()):
        dump(out / 'status.json', dict(status='FAILED',
             error='baseline does not reproduce frozen certification',
             checks=checks, baseline={k: v for k, v in b.items()
                                      if not k.startswith('_')}))
        raise RuntimeError(f'Baseline guard failed: {checks}')
    print(f'  baseline margin {b["margin_h"]!r} ; runtime {b["runtime_s"]:.1f} s',
          flush=True)
    rows.append(base)

    # ---- pre-registered degeneracy / null controls -----------------------
    controls = {}
    if not args.no_controls:
        print('degeneracy + null controls', flush=True)
        cs = evaluate(prototype, 'DEGEN clock_scale=x0.5',
                      HbyConfig(**{**asdict(FROZEN), 'clock_scale': 0.5}),
                      ref_traj=base['_traj'])
        ck = evaluate(prototype, 'DEGEN clock_K_au=x2',
                      HbyConfig(**{**asdict(FROZEN), 'clock_K_au': 0.6}),
                      ref_traj=base['_traj'])
        gap_scale = float(np.max(np.abs(cs['_traj'] - ck['_traj'])))
        rel_scale = gap_scale / float(np.max(np.abs(base['_traj'])))
        controls['P1_clock_scale_vs_clock_K'] = dict(
            max_abs_gap=gap_scale, rel_gap=rel_scale,
            exact=bool(rel_scale < 1e-8),
            a=cs['label'], b=ck['label'],
            note='hill(s*x,K,n) == hill(x,K/s,n) algebraically')
        print(f'  P1 clock_scale=0.5 vs clock_K_au=0.6 : max|dy| = {gap_scale:.3e} '
              f'(rel {rel_scale:.2e})', flush=True)

        th = evaluate(prototype, 'CTRL translation_h=x2',
                      HbyConfig(**{**asdict(FROZEN), 'translation_h': 60.0}),
                      ref_traj=base['_traj'])
        # mRNA states: b0_M_*, b2_M_* (contain '_M_') and M_A1 / M_F1.
        is_mrna = [i for i, n in enumerate(NAMES)
                   if n.startswith('M_') or '_M_' in n]
        is_mature = [i for i in range(len(NAMES)) if i not in is_mrna]
        y_scale = float(np.max(np.abs(base['_traj'])))
        gap_mat = float(np.max(np.abs(th['_traj'][is_mature]
                                      - base['_traj'][is_mature])))
        gap_mrna = float(np.max(np.abs(th['_traj'][is_mrna]
                                       - base['_traj'][is_mrna])))
        controls['P2_translation_h_null'] = dict(
            max_abs_gap_non_mrna=gap_mat, max_abs_gap_mrna=gap_mrna,
            rel_gap_non_mrna=gap_mat / y_scale, state_scale=y_scale,
            null_on_mature=bool(gap_mat / y_scale < 1e-6),
            differs_in_mrna=bool(gap_mrna / y_scale > 1e-3),
            n_mrna_states=len(is_mrna), mrna_states=[NAMES[i] for i in is_mrna])
        print(f'  P2 translation_h=x2 : max|dy| non-mRNA = {gap_mat:.3e} '
              f'(rel {gap_mat / y_scale:.2e}) ; mRNA = {gap_mrna:.3e}', flush=True)

        # Exploratory: the a.u. scale of the Int0 pool trades against the
        # thresholds it feeds. Under uM_per_au -> lam*uM_per_au the bit0 chain
        # scales as I -> I/lam and (at quasi-steady state) C -> C/lam**2, so
        # invariance requires K_D_int[0] -> /lam, K_inh -> /lam,
        # K_D_comp -> /lam**2 (k_complex -> /lam**2) and clock_K_au -> /lam.
        # The explicit complex ODE breaks this at O(kon*I*R), so the residual
        # gap measures the departure from the quasi-steady-state assumption.
        lam = 2.0
        joint_cfg = HbyConfig(**{**asdict(FROZEN),
                                 'receiver_uM_per_au': FROZEN.receiver_uM_per_au * lam,
                                 'clock_K_au': FROZEN.clock_K_au / lam})
        joint_patch = {('K_D_int', 0): ZENG['K_D_int'][0] / lam,
                       ('K_inh', None): ZENG['K_inh'] / lam,
                       ('K_D_comp', None): ZENG['K_D_comp'] / lam ** 2}
        jr = evaluate(prototype, 'DEGEN uM_per_au joint rescale', joint_cfg,
                      joint_patch, ref_traj=base['_traj'])
        gap_j = float(np.max(np.abs(jr['_traj'] - base['_traj'])))
        controls['P3_au_scale_joint_rescaling'] = dict(
            lambda_=lam, max_abs_gap=gap_j, rel_gap=gap_j / y_scale,
            exact=False, exploratory=True,
            traded_against={'K_D_int[0]': '/lam', 'K_inh': '/lam',
                            'K_D_comp': '/lam**2', 'clock_K_au': '/lam'},
            note=('algebraically exact only under complex quasi-steady state; '
                  'the residual is the explicit-complex correction'),
            certified_changed=bool(jr['certified'] != base['certified']))
        print(f'  P3 a.u. joint rescaling : max|dy| = {gap_j:.3e}', flush=True)
        rows += [cs, ck, th, jr]

    # ---- OAT scan --------------------------------------------------------
    spec = spec_rows()
    if args.runs:
        spec = spec[:args.runs]
    print(f'OAT rows: {len(spec)}', flush=True)
    for i, (label, tier, kind, target, value) in enumerate(spec, 1):
        if kind == 'config':
            cfg = HbyConfig(**{**asdict(FROZEN), target: float(value)})
            patch = None
        else:
            cfg = FROZEN
            patch = {target: float(value)}
        r = evaluate(prototype, label, cfg, patch, ref_traj=base['_traj'])
        r.update(tier=tier, kind=kind, target=str(target), value=float(value))
        bm = r.get('bit_margins') or {}
        s1m = min_of(bm.get('S1'))
        s2m = min_of(bm.get('S2'))
        r['margin_S1_h'] = s1m
        r['margin_S2_h'] = s2m
        rows.append(r)
        print(f'[{i}/{len(spec)}] {label:34s} tier={tier:4s} '
              f'cert={str(r["certified"]):5s} reads={r["steady_reads"]:3d} '
              f'S1m={("None" if s1m is None else format(s1m, ".4f")):>8s} '
              f'S2m={("None" if s2m is None else format(s2m, ".4f")):>8s} '
              f'gap={r["max_abs_gap_vs_baseline"]:.3e} {r["runtime_s"]:.1f}s',
              flush=True)

    # ---- persist ---------------------------------------------------------
    slim = [{k: v for k, v in r.items() if not k.startswith('_')} for r in rows]
    dump(out / 'oat_all.json', dict(hours=HOURS, sample_min=SAMPLE_MIN,
         max_step_min=MAX_STEP_MIN, rtol=RTOL, atol=ATOL,
         frozen_config=asdict(FROZEN), reference=REF,
         controls=controls, harness_timing=dict(han_s=han_s, prototype_s=proto_s),
         rows=slim))

    import csv
    def bmval(r, bit, key):
        mm = (r.get('bit_margins') or {}).get(bit) or {}
        return mm.get(key)

    cols = ['label', 'tier', 'kind', 'target', 'value', 'certified',
            'steady_reads', 'cold_reads', 'increments_mod8', 'sequence',
            'minimum_commitment', 'boundary_clips', 'margin_h',
            'margin_S1_h', 'margin_S2_h',
            'S0_setup_h', 'S0_hold_h', 'S1_setup_h', 'S1_hold_h',
            'S2_setup_h', 'S2_hold_h',
            'clock_peak_count', 's0_rev', 's0_gate', 's0_flip',
            's1_rev', 's1_gate', 's1_flip',
            'max_abs_gap_vs_baseline', 'max_abs_gap_int0',
            'runtime_s', 'error']
    with (out / 'oat_all.csv').open('w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for r in rows:
            s0 = r.get('stage0') or {}
            s1 = r.get('stage1') or {}
            bm = r.get('bit_margins') or {}
            s1m = min_of(bm.get('S1'))
            s2m = min_of(bm.get('S2'))
            w.writerow([r['label'], r.get('tier', 'CTRL'), r.get('kind', ''),
                        r.get('target', ''), r.get('value', ''), r['certified'],
                        r['steady_reads'], r['cold_reads'], r['increments_mod8'],
                        r['sequence'], r['minimum_commitment'], r['boundary_clips'],
                        r['margin_h'], s1m, s2m,
                        bmval(r, 'S0', 'min_setup_h'), bmval(r, 'S0', 'min_hold_h'),
                        bmval(r, 'S1', 'min_setup_h'), bmval(r, 'S1', 'min_hold_h'),
                        bmval(r, 'S2', 'min_setup_h'), bmval(r, 'S2', 'min_hold_h'),
                        r['clock_peak_count'],
                        s0.get('reverse'), s0.get('gate'), s0.get('flip'),
                        s1.get('reverse'), s1.get('gate'), s1.get('flip'),
                        r.get('max_abs_gap_vs_baseline'),
                        r.get('max_abs_gap_int0'),
                        round(r['runtime_s'], 3), r['error']])

    for r in rows:
        r.pop('_traj', None)
        r.pop('_signals', None)
    dump(out / 'status.json', dict(status='COMPLETED', rows=len(rows),
         controls_ok=all(v.get('exact', v.get('null_on_mature', True))
                         for v in controls.values()) if controls else None))
    print('OUTPUT', out, flush=True)


if __name__ == '__main__':
    main()
