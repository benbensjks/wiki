"""Round 1 failure attribution for the HBY-ZMH-ZMH (HZZ) 23-state cascade.

READ-ONLY WITH RESPECT TO EXISTING FILES
----------------------------------------
This script defines its own probes; it does not modify model_hzz.py,
run_comparison.py, the published results directory, or any other existing
file. Everything it writes goes under failure_attribution/results/<stamp>/.

WHAT IT TESTS
-------------
The published round is a whole-module tail replacement (expression layer, DNA
recombination response, parameter table and the clock input source all change
at once), so its failure cannot be attributed to one parameter AS RUN. This
script performs the missing single-lever decomposition.

Pre-registered fault tree (from reading the published artifacts):

  F1  RDF2 pool too small.  K_rep = 0.45 (ZENG uses 0.85) with n = 3.9, while
      tr settles at alpha_rep*pb2/gamma_rep = 3.5*0.182/0.6 = 1.063, gives a
      suppression 1/(1+(1.063/0.45)^3.9) = 1/29.6, hence
      r_ss = 6*0.818*0.0338/0.8 = 0.207 -- which is where the published RDF2
      peak (0.2071) actually sits. The reverse Hill
      (i*r)^2/(K_D_comp^2+(i*r)^2) then sits at ~1% of saturation, so the
      reverse drive is only ~1.06x the forward drive and bit2 cannot flip back.

  F2  Gate dose too small.  Per-pulse dose 0.0137 (off) / 0.0156 (on) against
      the 0.02 event threshold and against 0.237 per cycle in the certified
      HZH tail. Two sub-causes:
      (a) pulse_gate = int0^3/(0.4^3+int0^3) is faithful to the source file,
          but the source drives it with a SQUARE WAVE whose plateau is
          k_int/gamma_int = 6/2 = 3.0 (pulse_gate = 0.998), whereas the HBY
          bit0 Int0 peaks at 0.5485 (pulse_gate = 0.7205);
      (b) F1 has no expression chain in the 6-state tail, so it crosses
          K_R1 = 0.4 within ~0.17 h and shuts the gate almost immediately.

Probes (single lever unless combined; every other value frozen):
  base           as published
  krep85         p['K_rep']      0.45 -> 0.85      (tests F1)
  aint38         p['alpha_Int2'] 28   -> 38        (tests F2 dose)
  clk020         clock_K         0.40 -> 0.20      (tests F2a)
  clk007         clock_K         0.40 -> 0.07      (reproduces the source
                 file's near-unity pulse_gate at the HBY Int0 peak)
  krep85_aint38  F1 + dose
  krep85_clk007  F1 + F2a

PRE-REGISTERED READING RULES
----------------------------
  * If krep85 alone restores reverse flips and mod-8, F1 is the dominant cause.
  * If only clk007 (or aint38) restores it, F2 dominates.
  * If nothing in this round restores it, the tail's STRUCTURE (no expression
    chain, no explicit complex) is implicated beyond any single parameter, and
    the next lever is structural, not parametric.
  * In every case the three-tier split is reported: counting / event / full.
    A run whose code string is right but whose event arm fails is NOT a
    counting failure.
"""
from __future__ import annotations

import csv
import json
import platform
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import solve_ivp

sys.dont_write_bytecode = True

ROOT = Path(__file__).resolve().parent          # failure_attribution
HZZ = ROOT.parent                               # hby_zmh_zmh
PARENT = HZZ.parent                             # 跨模型三级连接验证
for _d in (HZZ, PARENT, PARENT / 'hby_zmh_hby',
           PARENT / 'hby_zmh_hby' / 'certification'):
    sys.path.insert(0, str(_d))

from model_hzz import (HZZModel, HZHModel, NAMES, IDX,          # noqa: E402
                       TAIL_SOURCE, structural_checks)
from run_han_comparison import HanInput, HAN_PATH               # noqa: E402
from hybrid_model import sha256, source_hashes                  # noqa: E402
import verify_hzh as V                                          # noqa: E402

HOURS = 300.0
SOLVER = dict(method='DOP853', rtol=2e-7, atol=2e-9, max_step=1 / 60)

# The published round, used as the reproduction guard.
PUBLISHED = {
    'off': dict(steady='xxx4567456745674567', events=[7, 0, 1], s2_flips=1),
    'on': dict(steady='5674567456745674567', events=[7, 0, 1], s2_flips=1),
}

# (probe id, overrides, arms)  overrides: p-dict keys or the 'clock_K' pseudo-key
PROBES = [
    ('base',          {},                                   ('off', 'on')),
    ('krep85',        {'K_rep': 0.85},                      ('off', 'on')),
    ('aint38',        {'alpha_Int2': 38.0},                 ('off', 'on')),
    ('clk020',        {'clock_K': 0.20},                    ('on',)),
    ('clk007',        {'clock_K': 0.07},                    ('on',)),
    ('krep85_aint38', {'K_rep': 0.85, 'alpha_Int2': 38.0},  ('on',)),
    ('krep85_clk007', {'K_rep': 0.85, 'clock_K': 0.07},     ('on',)),
]


def dump(path, obj):
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2,
                               allow_nan=False), encoding='utf-8')


def integrate(model):
    t = np.linspace(0, HOURS, int(HOURS * 60) + 1)
    sol = solve_ivp(model.rhs, (0, HOURS), model.initial_state(), t_eval=t,
                    **SOLVER)
    if not sol.success or not np.isfinite(sol.y).all():
        raise RuntimeError(sol.message)
    return sol.t, sol.y


# --------------------------------------------------------------------------
# criterion (identical to the published harness, so results are comparable)
# --------------------------------------------------------------------------
def evaluate(t, sig):
    peaks, reads = V.read_windows(t, sig)
    crosses = {b: V.crossings(t, sig[b]) for b in ('S0', 'S1', 'S2')}
    events = {}
    for j in (0, 1):
        rev = V.segments(t, sig[f'J_rev{j}'], V.JREV_THRESHOLD_H)
        gate = V.segments(t, sig[f'g{j}'], V.GATE_THRESHOLD,
                          min_dose=V.GATE_MIN_DOSE_H)
        events[f'bit{j}_to_bit{j+1}'] = V.associate(rev, gate,
                                                    crosses[f'S{j+1}'], HOURS)
    steady = V.read_verdict(reads, V.DROP)
    event_pass = all(x['passed'] for x in events.values())
    perbit = {}
    for i, b in enumerate(('S0', 'S1', 'S2')):
        sel = reads[V.DROP:]
        perbit[b] = dict(minimum_commitment=min(r['commitment'][i] for r in sel),
                         unlabelled=sum(r['labels'][i] is None for r in sel),
                         minimum=float(sig[b].min()),
                         maximum=float(sig[b].max()),
                         timing=V.margins(sel, crosses[b]))
    return dict(version='HZZ_ROUND1_ATTRIBUTION_V1', certified=None,
                counting_passed=bool(steady['passed']),
                event_causality_passed=bool(event_pass),
                diagnostic_passed=bool(steady['passed'] and event_pass),
                cold=V.read_verdict(reads, 0), steady=steady, per_bit=perbit,
                crossings=crosses, events=events, read_windows=reads,
                clock_peaks=len(peaks))


# --------------------------------------------------------------------------
# diagnostics that discriminate the fault tree
# --------------------------------------------------------------------------
def raw_gate_segments(t, g1):
    """Every g1 segment above the event threshold, with its own duration/dose.

    V.segments is called with zero minimums so NOTHING is filtered out; the
    two acceptance thresholds are then reported per segment rather than
    silently dropping the segments that fail them.
    """
    segs = V.segments(t, g1, V.GATE_THRESHOLD, min_duration_h=0.0, min_dose=0.0)
    out = []
    for s in segs:
        dur = s['end_h'] - s['start_h']
        out.append(dict(start_h=s['start_h'], end_h=s['end_h'],
                        peak_h=s['peak_h'], peak=s['peak'], dose=s['dose'],
                        duration_h=dur,
                        meets_duration=bool(dur >= V.PULSE_DURATION_H),
                        meets_dose=bool(s['dose'] >= V.GATE_MIN_DOSE_H)))
    return out


def mechanism(y, p, clock_K):
    """Closed-form check of the F1 hypothesis against the measured trajectory."""
    i_pk = float(y[IDX['I2_zmh']].max())
    r_pk = float(y[IDX['RDF2_zmh']].max())
    tr_min = float(y[IDX['T2_zmh']].min())
    pb_min = float(y[IDX['pb2_zmh']].min())
    supp = 1.0 / (1.0 + (tr_min / p['K_rep']) ** p['n'])
    r_ss_pred = p['alpha_rdf'] * (1 - pb_min) * supp / p['gamma_rdf']
    ir = i_pk * r_pk
    hill = ir ** 2 / (p['K_D_comp'] ** 2 + ir ** 2)
    vr = p['k_rev'] * (1 - pb_min) * hill
    vf = (p['k_fwd'] * pb_min
          * (i_pk ** 2 / (p['K_D_int2'] ** 2 + i_pk ** 2))
          * (p['Kinh'] / (p['Kinh'] + r_pk)))
    return dict(
        I2_peak=i_pk, RDF2_peak=r_pk, T2_min=tr_min, pb2_min=pb_min,
        K_rep=p['K_rep'], n=p['n'], alpha_Int2=p['alpha_Int2'],
        K_A1=p['K_A1'], K_R1=p['K_R1'], clock_K=clock_K,
        rdf_suppression=supp, RDF2_steady_predicted=r_ss_pred,
        RDF2_pred_vs_measured=(r_ss_pred / r_pk if r_pk > 0 else None),
        i_times_r_peak=ir, K_D_comp=p['K_D_comp'], reverse_hill_value=hill,
        v_reverse_at_peak_per_h=vr, v_forward_at_peak_per_h=vf,
        reverse_over_forward=(vr / vf if vf > 0 else None),
        reverse_sufficient=bool(vr > vf))


def summarise(t, y, sig, res, clock_K, overrides):
    segs = raw_gate_segments(t, sig['g1'])
    mech = mechanism(y, sig['_p'], clock_K)
    ev = res['events']['bit1_to_bit2']
    return dict(
        clock_K=clock_K, overrides=overrides,
        counting_passed=res['counting_passed'],
        event_causality_passed=res['event_causality_passed'],
        diagnostic_passed=res['diagnostic_passed'],
        cold=res['cold']['sequence'], steady=res['steady']['sequence'],
        stage0=[res['events']['bit0_to_bit1']['reverse_events'],
                res['events']['bit0_to_bit1']['gate_events'],
                res['events']['bit0_to_bit1']['flip_events']],
        stage1=[ev['reverse_events'], ev['gate_events'], ev['flip_events']],
        s2_flips=len(res['crossings']['S2']),
        gate_segment_count=len(segs),
        gate_duration_range=[min(s['duration_h'] for s in segs),
                             max(s['duration_h'] for s in segs)] if segs else None,
        gate_dose_range=[min(s['dose'] for s in segs),
                         max(s['dose'] for s in segs)] if segs else None,
        gate_segments_meeting_duration=sum(s['meets_duration'] for s in segs),
        gate_segments_meeting_dose=sum(s['meets_dose'] for s in segs),
        gate_segments_meeting_both=sum(s['meets_duration'] and s['meets_dose']
                                       for s in segs),
        g1_peak=float(sig['g1'].max()),
        clock_peak=float(sig['clock'].max()),
        clock_duty_gt_05=float(np.mean(sig['clock'] > 0.5)),
        int2_source_peak=float(sig['int2_source'].max()),
        J_fwd2_total=float(np.trapezoid(np.maximum(sig['J_fwd2'], 0), t)),
        J_rev2_total=float(np.trapezoid(np.maximum(sig['J_rev2'], 0), t)),
        mechanism=mech)


def main():
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8')
    out = ROOT / 'results' / datetime.now().strftime('%Y%m%d_%H%M%S')
    out.mkdir(parents=True, exist_ok=False)
    dump(out / 'status.json', dict(status='running'))

    paths = [HZZ / 'model_hzz.py', HZZ / 'run_comparison.py', TAIL_SOURCE,
             PARENT / 'hybrid_model.py', PARENT / 'run_han_comparison.py',
             PARENT / 'hby_zmh_hby' / 'model_hzh.py',
             PARENT / 'hby_zmh_hby' / 'certification' / 'verify_hzh.py',
             HAN_PATH]
    before = {str(p): sha256(p) for p in paths}
    before.update(source_hashes())

    print('build shared %g h Han input' % HOURS, flush=True)
    han = HanInput(HOURS)
    prefix = HZHModel('han', han)

    checks = structural_checks(prefix)
    dump(out / 'structural_checks.json', checks)
    print('structural checks:', json.dumps(checks), flush=True)

    # ---- certified HZH baseline guard (the prefix must still be the good one)
    tb, yb = integrate(prefix)
    bv, _ = V.analyse(tb, yb, prefix.z, prefix.tail)
    guard = dict(certified=bv['certified_v1'],
                 reads=bv['steady']['reads'] == 19,
                 sequence=bv['steady']['sequence'] == '1234567012345670123',
                 margin=abs(bv['global_min_timing_margin_h']
                            - 1.1784609018563117) < 1e-6)
    print('HZH prefix guard:', json.dumps(guard), flush=True)
    if not all(guard.values()):
        dump(out / 'status.json', dict(status='FAILED', guard=guard))
        raise RuntimeError(f'HZH prefix guard failed: {guard}')

    rows = []
    reproduces = {}
    for probe, overrides, arms in PROBES:
        for arm in arms:
            label = f'{probe}|{arm}'
            model = HZZModel(prefix, arm == 'on')
            # single-lever overrides; each HZZModel re-parses the source dict,
            # so mutating model.p here cannot leak into another run
            clock_K = float(model.clock_K)
            for k, v in overrides.items():
                if k == 'clock_K':
                    model.clock_K = float(v)
                    clock_K = float(v)
                else:
                    model.p[k] = v
            t, y = integrate(model)
            sig = model.signals(y)
            sig['_p'] = model.p
            if y.min() < -1e-7:
                raise RuntimeError(f'{label}: nonphysical state')
            res = evaluate(t, sig)
            res['prefix_trajectory_gap_vs_HZH'] = float(np.max(np.abs(y[:17] - yb[:17])))
            res['state_extrema'] = {n: [float(y[i].min()), float(y[i].max())]
                                    for i, n in enumerate(NAMES)}
            row = dict(label=label, probe=probe, arm=arm, **overrides)
            row.update(summarise(t, y, sig, res, clock_K, overrides))
            row['prefix_gap'] = res['prefix_trajectory_gap_vs_HZH']
            row['per_bit'] = res['per_bit']
            rows.append(row)

            if probe == 'base':
                want = PUBLISHED[arm]
                ok = (row['steady'] == want['steady']
                      and row['stage1'] == want['events']
                      and row['s2_flips'] == want['s2_flips'])
                reproduces[arm] = dict(
                    ok=bool(ok), want=want,
                    got=dict(steady=row['steady'], stage1=row['stage1'],
                             s2_flips=row['s2_flips']))

            np.savez_compressed(out / f'traj_{probe}_{arm}.npz',
                                time_h=t, states=y, state_names=np.array(NAMES))
            dump(out / f'verdict_{probe}_{arm}.json', res)

            m = row['mechanism']
            print(f"[{label:16s}] K={clock_K:.2f} Krep={m['K_rep']:.2f} "
                  f"cnt={int(row['counting_passed'])} "
                  f"ev={int(row['event_causality_passed'])} "
                  f"S2flips={row['s2_flips']} "
                  f"seq={row['steady'][:20]:20s} "
                  f"g1pk={row['g1_peak']:.4f} "
                  f"dur={row['gate_duration_range']} "
                  f"dose={row['gate_dose_range']} "
                  f"RDF2pk={m['RDF2_peak']:.4f} "
                  f"pred={m['RDF2_steady_predicted']:.4f} "
                  f"vr/vf={m['reverse_over_forward']:.3f}", flush=True)

    dump(out / 'reproduction_guard.json', reproduces)
    if not all(v['ok'] for v in reproduces.values()):
        dump(out / 'status.json', dict(status='FAILED',
                                       error='published round not reproduced',
                                       reproduction=reproduces))
        raise RuntimeError(f'Reproduction guard failed: {reproduces}')
    print('published round reproduced for both arms', flush=True)

    # ------------------------------------------------------------- verdict
    base_off = next(r for r in rows if r['label'] == 'base|off')
    verdicts = {}
    for r in rows:
        if r['probe'] == 'base':
            continue
        verdicts[r['label']] = dict(
            flips_restored=bool(r['s2_flips'] > PUBLISHED[r['arm']]['s2_flips']),
            counting_restored=bool(r['counting_passed']),
            event_restored=bool(r['event_causality_passed']),
            full_restored=bool(r['diagnostic_passed']),
            reverse_over_forward=r['mechanism']['reverse_over_forward'],
            RDF2_peak=r['mechanism']['RDF2_peak'],
            gate_segments_meeting_both=r['gate_segments_meeting_both'])

    dump(out / 'round1_summary.json', dict(
        hours=HOURS, solver=SOLVER, n_states=len(NAMES), state_names=list(NAMES),
        published=PUBLISHED, reproduction_guard=reproduces,
        hzh_prefix_guard=guard, structural_checks=checks,
        probes=[dict(probe=p, overrides=o, arms=list(a)) for p, o, a in PROBES],
        reading_rules=(
            'krep85 alone restoring reverse flips => F1 (RDF2 pool) dominates; '
            'only clk007/aint38 restoring => F2 (gate dose) dominates; nothing '
            'restoring => the tail STRUCTURE is implicated beyond any single '
            'parameter. counting / event / full are always reported separately.'),
        baseline_off=base_off, verdicts=verdicts, rows=rows,
        source_hashes=before, python=sys.version, numpy=np.__version__,
        scipy=scipy.__version__, platform=platform.platform()))

    cols = ['label', 'probe', 'arm', 'clock_K', 'K_rep', 'alpha_Int2',
            'counting_passed', 'event_causality_passed', 'diagnostic_passed',
            'cold', 'steady', 's2_flips', 'stage0_r/g/f', 'stage1_r/g/f',
            'gate_segment_count', 'gates_meet_duration', 'gates_meet_dose',
            'gates_meet_both', 'gate_dur_min', 'gate_dur_max',
            'gate_dose_min', 'gate_dose_max', 'g1_peak', 'clock_peak',
            'clock_duty_gt_05', 'int2_source_peak', 'J_fwd2_total',
            'J_rev2_total', 'I2_peak', 'RDF2_peak', 'T2_min', 'pb2_min',
            'rdf_suppression', 'RDF2_predicted', 'pred_over_measured',
            'i_times_r_peak', 'reverse_hill', 'v_rev', 'v_fwd', 'vr_over_vf',
            'prefix_gap']
    with (out / 'round1_all.csv').open('w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for r in rows:
            m = r['mechanism']
            dur = r['gate_duration_range'] or [None, None]
            dose = r['gate_dose_range'] or [None, None]
            w.writerow([
                r['label'], r['probe'], r['arm'], r['clock_K'], m['K_rep'],
                m['alpha_Int2'], r['counting_passed'],
                r['event_causality_passed'], r['diagnostic_passed'], r['cold'],
                r['steady'], r['s2_flips'],
                '%d/%d/%d' % tuple(r['stage0']), '%d/%d/%d' % tuple(r['stage1']),
                r['gate_segment_count'],
                r['gate_segments_meeting_duration'],
                r['gate_segments_meeting_dose'],
                r['gate_segments_meeting_both'], dur[0], dur[1], dose[0],
                dose[1], r['g1_peak'], r['clock_peak'], r['clock_duty_gt_05'],
                r['int2_source_peak'], r['J_fwd2_total'], r['J_rev2_total'],
                m['I2_peak'], m['RDF2_peak'], m['T2_min'], m['pb2_min'],
                m['rdf_suppression'], m['RDF2_steady_predicted'],
                m['RDF2_pred_vs_measured'], m['i_times_r_peak'],
                m['reverse_hill_value'], m['v_reverse_at_peak_per_h'],
                m['v_forward_at_peak_per_h'], m['reverse_over_forward'],
                r['prefix_gap']])

    after = {str(p): sha256(p) for p in paths}
    after.update(source_hashes())
    if before != after:
        raise RuntimeError('an existing source file changed during the run')
    dump(out / 'source_hashes.json', before)
    dump(out / 'status.json', dict(status='COMPLETED', rows=len(rows)))
    dump(out / 'SHA256SUMS.json',
         {str(p.relative_to(out)): sha256(p)
          for p in sorted(out.rglob('*')) if p.is_file()})
    print('OUTPUT', out, flush=True)


if __name__ == '__main__':
    main()
