"""Boundary probes P1-P4 for the frozen 34-state hybrid counter.

NEW FILE (added 2026-09-28).  It does NOT modify any existing script, model,
parameter file or result folder.  Everything it produces goes under
`跨模型三级连接验证/版本与鲁棒性探针/results/`.

WHY PROBES AND NOT A SCAN
-------------------------
The unknown is not in the parameter space, it is in the INTERFACES.  A 300 h run
of this model costs about 16 s, so a small pre-registered set of single-lever
perturbations is cheap; a grid would still answer the wrong question.

PRE-REGISTERED CASES (fixed before the runs)
  baseline            frozen candidate, Han v53d, no change
  p1_rdf1_x2          non-self-consistent initial state: RDF1_zmh x2 only
  p1_rdf1_x0p5        non-self-consistent initial state: RDF1_zmh x0.5 only
  p1_a1_x0p5          non-self-consistent initial state: A1 x0.5 only
  p2_load_0p10        Han upstream `peak_load_fraction` 0 -> 0.10
  p3_um_5p00          receiver_uM_per_au 5.75 -> 5.00
  p3_um_6p50          receiver_uM_per_au 5.75 -> 6.50

DELIBERATELY NOT PROBED
  * clock_scale -- the clock is driven by this model's own mature free Int0
    (b0_I), so the scale is an identity, not a free lever.  Sweeping it would be
    a fake question; that is exactly why earlier configurations failed.
  * n_A1_gate -- two completed experiments already localise the failure to the
    bit1/interface side, not the bit2 receiver.  Turning this knob to obtain a
    success would be buying a result with a parameter.

WHAT EVERY CASE MUST REPORT (the hard requirement)
--------------------------------------------------
300 h yields only ~19 steady read windows, so a bare pass/fail cannot separate
"the dynamics failed" from "there were too few windows to judge".  Each case
therefore reports, separately: the sweep verdict, the number of read windows,
the boundary-clip count, the code sequences, both event-chain counts with their
one-to-one / order / alternation arms, and the timing margins in hours.

P4 (molecule numbers, deterministic, no SSA)
-------------------------------------------
Mirrors the 51-state Rule C method (see plausibility/preaudit_stochasticity.py):
copies = a.u. * copies_per_au with copies_per_au = 602.214076 * uM_per_au;
statistics split into gate-phase and whole-cycle; burn-in 100 h; < 1 copy is
STRUCTURALLY ABSENT and carries no CV; < 100 copies is "small numbers"
(>= 10 % intrinsic noise), >= 1000 is large.  Like the 51-state package this
REFUSES to name a single binding pool, because every pool here is pulsatile.
It is a pre-check for whether a stochastic study is needed at all, not a
substitute for one.
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from dataclasses import asdict, replace
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')

WIKI = Path(r'C:\Users\18633\Desktop\wiki')
CROSS = WIKI / '跨模型三级连接验证'
HZH = CROSS / 'hby_zmh_hby'
HERE = Path(__file__).resolve().parent
OUT = HERE / 'results'

for _p in (CROSS, HZH, HZH / 'certification', WIKI / 'final_reconstruction',
           WIKI / 'final_reconstruction' / 'plausibility'):
    sys.path.insert(0, str(_p))

from hybrid_model import HbyConfig, HbyReceiver, require_hash, sha256, source_hashes  # noqa: E402
from model_hzh import HZHModel, IDX                                                   # noqa: E402
from run_han_comparison import HAN_PATH, HAN_SHA, HanInput                            # noqa: E402
from verify_hzh import analyse, signals as raw_signals                                # noqa: E402
from plausibility_common import gate_windows                                          # noqa: E402

HOURS = 300.0
SAMPLE_MIN = 1.0
MAX_STEP_MIN = 1.0
RTOL, ATOL = 2e-7, 2e-9
BURN_IN_H = 100.0
COPIES_PER_UM_PER_FL = 602.214076
COLUMNS = 1.0  # 1 fL, same convention as the frozen packages
SMALL_NUMBER_THRESHOLD = 100.0
LARGE_NUMBER_THRESHOLD = 1000.0
INVENTORY_STATES = ('b2_I', 'b2_R', 'b2_C', 'A1', 'F1')
LOAD_FRACTION = 0.10
PROBED_IC = (('p1_rdf1_x2', 'RDF1_zmh', 2.0),
             ('p1_rdf1_x0p5', 'RDF1_zmh', 0.5),
             ('p1_a1_x0p5', 'A1', 0.5))
PROBED_UM = (('p3_um_5p00', 5.0), ('p3_um_6p50', 6.5))
CASE_TAGS = (['baseline'] + [t for t, _, _ in PROBED_IC] + ['p2_load_0p10']
             + [t for t, _ in PROBED_UM])


def jsonable(o):
    if isinstance(o, dict):
        return {str(k): jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [jsonable(v) for v in o]
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (bool, int, float, str)) or o is None:
        return o
    return str(o)


def integrate(model, hours, sample_min, max_step_min, rtol, atol):
    count = round(hours * 60 / sample_min)
    t = np.linspace(0, hours, count + 1)
    if len(t) < 2 or not np.all(np.diff(t) > 0):
        raise ValueError('Bad sampling grid')
    sol = solve_ivp(model.rhs, (0, hours), model.initial_state(), t_eval=t,
                    method='DOP853', rtol=rtol, atol=atol, max_step=max_step_min / 60)
    if not sol.success or not np.all(np.isfinite(sol.y)):
        raise RuntimeError(sol.message)
    return sol.t, sol.y


def han_with_load(hours, load):
    """Build a Han upstream with a non-zero resource load.

    The reference demand peak is computed exactly as the upstream module's own
    driver does it (its `main()`), not re-invented here.
    """
    require_hash(HAN_PATH, HAN_SHA)
    spec = importlib.util.spec_from_file_location('probe_han_reference', HAN_PATH)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    duration = max(3000.0, 60 * hours)
    t = np.r_[np.arange(0, duration, 1.0), duration]
    base = module.make_reference_parameters()
    if load == 0.0:
        reference_peak = 1.0          # rho is exactly 1 when load == 0
    else:
        _, _, _, _, reference_demand, _ = module.simulate(base, t, 1.0)
        keep = t >= module.SETTINGS['transient_min']
        reference_peak = float(np.max(reference_demand[keep]))
    p = replace(base, peak_load_fraction=float(load))
    sol, *_ = module.simulate(p, t, reference_peak)
    han = type('LoadedHan', (), {})()
    han.module, han.p, han.t, han.sol = module, p, t, sol
    han.normalization_mask = (t >= 1000) & (t <= 3000)
    han.c31_max = float(sol.y[7, han.normalization_mask].max())
    han.flux_copies_min = p.c31_translation_per_mrna_per_min * sol.y[6]
    han.copies_per_uM = COPIES_PER_UM_PER_FL
    han.reference_peak = reference_peak
    han.interp = lambda tt, sig, _h=han: float(np.interp(tt, _h.t, sig))
    return han


def molecule_inventory(t, y, g1, copies_per_au):
    """The 51-state Rule C statistics, applied to this model's pools."""
    wins = gate_windows(t, g1, threshold=0.01, min_duration_h=0.05)
    post = np.asarray(t) >= BURN_IN_H
    in_gate = np.zeros(len(t), dtype=bool)
    for w in wins:
        in_gate[w['i0']:w['i1'] + 1] = True
    decision = in_gate & post
    pools, binding = {}, []
    for name in INVENTORY_STATES:
        if name not in IDX:
            raise RuntimeError(f'{name} is not a state of the HZH layout')
        v = np.asarray(y[IDX[name]], dtype=float) * copies_per_au
        entry = dict(copies_min_gate=(float(v[decision].min()) if decision.any() else None),
                     copies_p05_gate=(float(np.percentile(v[decision], 5)) if decision.any() else None),
                     copies_median_gate=(float(np.median(v[decision])) if decision.any() else None),
                     copies_min_post=float(v[post].min()),
                     copies_median_post=float(np.median(v[post])))
        for label, val in (('gate_phase_minimum', entry['copies_min_gate']),
                           ('gate_phase_p05', entry['copies_p05_gate']),
                           ('gate_phase_median', entry['copies_median_gate']),
                           ('whole_cycle_minimum', entry['copies_min_post']),
                           ('whole_cycle_median', entry['copies_median_post'])):
            if val is None:
                continue
            if val < 1.0:
                entry[label] = dict(copies=val, regime='STRUCTURALLY ABSENT',
                                    implied_cv_pct=None)
            else:
                entry[label] = dict(copies=val,
                                    implied_cv_pct=float(100.0 / np.sqrt(val)),
                                    regime=('small numbers' if val < SMALL_NUMBER_THRESHOLD
                                            else ('intermediate' if val < LARGE_NUMBER_THRESHOLD
                                                  else 'large numbers')))
                if val < SMALL_NUMBER_THRESHOLD:
                    binding.append(dict(pool=name, statistic=label, copies=val,
                                        implied_cv_pct=float(100.0 / np.sqrt(val))))
        pools[name] = entry
    binding.sort(key=lambda d: d['copies'])
    return dict(gate_windows=len(wins), gate_samples=int(decision.sum()),
                burn_in_h=BURN_IN_H, copies_per_au=copies_per_au,
                per_pool=pools, binding_small_number_entries=binding,
                rule=('per-pool, per-phase; NO single-pool verdict -- every pool here '
                      'is pulsatile, exactly as in the 51-state Rule C'))


def perturbed_initial_state(model, name, factor):
    if name not in IDX:
        raise RuntimeError(f'{name} is not a state of the HZH layout')
    base_initial = model.initial_state

    def scaled():
        y = base_initial()
        y[IDX[name]] = y[IDX[name]] * factor
        return y
    model.initial_state = scaled
    return model


def with_receiver_scale(model, um):
    model.tail = HbyReceiver(replace(HbyConfig(), receiver_uM_per_au=float(um)))
    model.copies_per_au = COPIES_PER_UM_PER_FL * float(um) * COLUMNS
    return model


def build_cases(han0, han_loaded):
    """(tag, factory -> (model, han_used), extra -> dict)."""
    cases = [('baseline', lambda: (HZHModel('han', han0), han0),
              lambda m, h: dict(interface=dict(mode='han', peak_load_fraction=0.0,
                                               receiver_uM_per_au=float(m.tail.c.receiver_uM_per_au))))]
    for tag, name, factor in PROBED_IC:
        cases.append((tag,
                      (lambda n=name, f=factor: (perturbed_initial_state(HZHModel('han', han0), n, f), han0)),
                      (lambda m, h, n=name, f=factor: dict(
                          interface=dict(mode='han', peak_load_fraction=0.0,
                                         receiver_uM_per_au=float(m.tail.c.receiver_uM_per_au)),
                          initial_state_change=dict(state=n, factor=f)))))
    cases.append(('p2_load_0p10', lambda: (HZHModel('han', han_loaded), han_loaded),
                  lambda m, h: dict(interface=dict(
                      mode='han', peak_load_fraction=float(h.p.peak_load_fraction),
                      reference_C31_demand_peak=float(h.reference_peak),
                      receiver_uM_per_au=float(m.tail.c.receiver_uM_per_au)))))
    for tag, um in PROBED_UM:
        cases.append((tag,
                      (lambda u=um: (with_receiver_scale(HZHModel('han', han0), u), han0)),
                      (lambda m, h, u=um: dict(interface=dict(
                          mode='han', peak_load_fraction=0.0, receiver_uM_per_au=u)))))
    return cases


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    meta = dict(script=Path(__file__).name,
                script_sha256=sha256(Path(__file__).resolve()),
                hours=HOURS, sample_min=SAMPLE_MIN, max_step_min=MAX_STEP_MIN,
                solver=dict(method='DOP853', rtol=RTOL, atol=ATOL),
                han_sha256=sha256(HAN_PATH), han_expected=HAN_SHA,
                receiver_config=jsonable(asdict(HbyConfig())),
                load_fraction=LOAD_FRACTION, probed_initial_states=PROBED_IC,
                probed_uM_per_au=PROBED_UM,
                source_hashes=source_hashes(),
                pre_registered_cases=CASE_TAGS,
                not_probed=['clock_scale (same-source clock: identity, not a lever)',
                            'n_A1_gate (failure already localised away from bit2)'])
    print('Preparing Han upstream (load = 0)...', flush=True)
    han0 = HanInput(HOURS)
    print('Preparing Han upstream (load = %.2f)...' % LOAD_FRACTION, flush=True)
    han_loaded = han_with_load(HOURS, LOAD_FRACTION)

    rows, verdicts, inventories = [], {}, {}
    for tag, factory, extra_factory in build_cases(han0, han_loaded):
        print(f'[run] {tag} ...', flush=True)
        tic = time.perf_counter()
        model, han_used = factory()
        t, y = integrate(model, HOURS, SAMPLE_MIN, MAX_STEP_MIN, RTOL, ATOL)
        runtime = time.perf_counter() - tic
        extra = extra_factory(model, han_used)
        verdict, sig = analyse(t, y, model.z, model.tail)
        steady, cold = verdict['steady'], verdict['cold']
        row = dict(tag=tag,
                   certified_v1=bool(verdict['certified_v1']),
                   steady_reads=int(steady['reads']),
                   steady_sequence=steady['sequence'],
                   steady_passed=bool(steady['passed']),
                   boundary_clips=int(steady['boundary_clips']),
                   minimum_commitment=steady['minimum_commitment'],
                   cold_reads=int(cold['reads']),
                   cold_sequence=cold['sequence'],
                   clock_peak_count=int(verdict['clock_peak_count']),
                   global_min_timing_margin_h=verdict['global_min_timing_margin_h'],
                   bit_margins=verdict['bit_margins'],
                   events={k: {j: x[j] for j in ('reverse_events', 'gate_events', 'flip_events',
                                                 'one_to_one', 'causal_order',
                                                 'alternating_directions', 'passed')}
                           for k, x in verdict['events'].items()},
                   signal_ranges=verdict['signal_ranges'],
                   g1_peak=float(np.max(sig['g1'])),
                   g1_median=float(np.median(sig['g1'])),
                   duty_g1_gt_0p5=float(np.mean(sig['g1'] > 0.5)),
                   duty_g1_gt_0p1=float(np.mean(sig['g1'] > 0.1)),
                   runtime_s=runtime)
        row.update(extra)
        rows.append(row)
        verdicts[tag] = jsonable(verdict)
        cpa = COPIES_PER_UM_PER_FL * float(model.tail.c.receiver_uM_per_au) * COLUMNS
        inventories[tag] = molecule_inventory(t, y, sig['g1'], cpa)
        print(f'      certified={row["certified_v1"]} reads={row["steady_reads"]} '
              f'clips={row["boundary_clips"]} margin={row["global_min_timing_margin_h"]} '
              f'({runtime:.1f} s)', flush=True)
        del model

    if not rows[0]['certified_v1']:
        raise RuntimeError('baseline is not certified_v1; refusing to report probes '
                           'against a mis-specified working point')

    (OUT / 'probe_P1_P4_hybrid.json').write_text(
        json.dumps(jsonable(dict(meta=meta, rows=rows, verdicts=verdicts,
                                 inventories=inventories)), ensure_ascii=False, indent=1),
        encoding='utf-8')

    hdr = ['tag', 'certified_v1', 'steady_reads', 'boundary_clips', 'minimum_commitment',
           'global_min_timing_margin_h', 'clock_peak_count', 'chain01_passed',
           'chain12_passed', 'duty_g1_gt_0p1', 'g1_peak', 'runtime_s']
    lines = [','.join(hdr)]
    for r in rows:
        lines.append(','.join(str(x) for x in [
            r['tag'], r['certified_v1'], r['steady_reads'], r['boundary_clips'],
            r['minimum_commitment'], r['global_min_timing_margin_h'],
            r['clock_peak_count'], r['events']['bit0_to_bit1']['passed'],
            r['events']['bit1_to_bit2']['passed'], r['duty_g1_gt_0p1'], r['g1_peak'],
            r['runtime_s']]))
    (OUT / 'probe_P1_P4_hybrid.csv').write_text('\n'.join(lines) + '\n', encoding='utf-8')
    print('\n'.join(lines))
    print('written to', OUT)


if __name__ == '__main__':
    main()
