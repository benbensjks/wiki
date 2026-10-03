"""Round 3: acceptance-grade re-check of the clock_K working point.

Round 2 bracketed the window: clock_K <= 0.20 passes on the A1-autoregulation
OFF arm, the upper boundary lies in (0.20, 0.25), and no lower boundary was
found down to 0.02. This round takes K = 0.10 -- inside the window with a 2.5x
margin, and with a better S2 band margin than 0.20 -- and re-checks it the way
the HZH work was checked:

  1. 600 h long run (300 h trajectory is also kept, for the phase source)
  2. eight phase-shifted initial states taken from ONE converged orbit period
  3. event-threshold sensitivity on the SAME trajectories (no re-integration),
     applied both to the passing point and to a known failing point

No existing file is modified. Output goes to
failure_attribution/results/round3_<stamp>/.
"""
from __future__ import annotations

import csv
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

sys.dont_write_bytecode = True

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT.parent))
sys.path.insert(0, str(ROOT.parent.parent))

import diagnose_round1 as D                                       # noqa: E402
import verify_hzh as V                                            # noqa: E402
from diagnose_round2 import band_margins                          # noqa: E402

K_WORK = 0.10
ARM = 'off'
HOURS_LONG = 600.0
N_PHASE = 8
# failing point, for the threshold-sensitivity contrast
K_FAIL = 0.25


def integrate_from(model, y0, hours):
    t = np.linspace(0, hours, int(hours * 60) + 1)
    sol = solve_ivp(model.rhs, (0, hours), np.asarray(y0, dtype=float),
                    t_eval=t, **D.SOLVER)
    if not sol.success or not np.isfinite(sol.y).all():
        raise RuntimeError(sol.message)
    return sol.t, sol.y


def evaluate(t, sig, hours):
    """Same criterion as the published harness, with an explicit horizon."""
    peaks, reads = V.read_windows(t, sig)
    crosses = {b: V.crossings(t, sig[b]) for b in ('S0', 'S1', 'S2')}
    events = {}
    for j in (0, 1):
        rev = V.segments(t, sig[f'J_rev{j}'], V.JREV_THRESHOLD_H)
        gate = V.segments(t, sig[f'g{j}'], V.GATE_THRESHOLD,
                          min_dose=V.GATE_MIN_DOSE_H)
        events[f'bit{j}_to_bit{j+1}'] = V.associate(rev, gate,
                                                    crosses[f'S{j+1}'], hours)
    steady = V.read_verdict(reads, V.DROP)
    perbit = {}
    for i, b in enumerate(('S0', 'S1', 'S2')):
        sel = reads[V.DROP:]
        perbit[b] = dict(minimum_commitment=min(r['commitment'][i] for r in sel),
                         unlabelled=sum(r['labels'][i] is None for r in sel),
                         minimum=float(sig[b].min()),
                         maximum=float(sig[b].max()),
                         timing=V.margins(sel, crosses[b]))
    return dict(counting_passed=bool(steady['passed']),
                event_causality_passed=bool(all(x['passed'] for x in events.values())),
                full_passed=bool(steady['passed']
                                 and all(x['passed'] for x in events.values())),
                cold=V.read_verdict(reads, 0), steady=steady, per_bit=perbit,
                crossings=crosses, events=events, read_windows=reads,
                clock_peaks=len(peaks))


def sens(t, y, model, hours):
    """Event-threshold sweep on a fixed trajectory (zero re-integration).

    Only JREV_THRESHOLD_H and GATE_THRESHOLD are varied: evaluate() passes both
    explicitly to V.segments, so patching the module globals takes effect.
    PULSE_DURATION_H is a default argument of V.segments and is therefore NOT
    reachable this way; it is left untouched and reported as such.
    """
    sig = model.signals(y)
    jrev0, gate0 = V.JREV_THRESHOLD_H, V.GATE_THRESHOLD
    try:
        out = []
        for jr in (0.02, 0.05, 0.1, 0.2, 0.4):
            for gt in (0.02, 0.05, 0.1, 0.2):
                V.JREV_THRESHOLD_H, V.GATE_THRESHOLD = jr, gt
                r = evaluate(t, sig, hours)
                out.append(dict(
                    jrev=jr, gate=gt,
                    counting=r['counting_passed'],
                    events=r['event_causality_passed'],
                    full=r['full_passed'],
                    sequence=r['steady']['sequence'],
                    stage1=[r['events']['bit1_to_bit2'][k] for k in
                            ('reverse_events', 'gate_events', 'flip_events')]))
    finally:
        V.JREV_THRESHOLD_H, V.GATE_THRESHOLD = jrev0, gate0
    return out


def main():
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8')
    out = ROOT / 'results' / ('round3_' + datetime.now().strftime('%Y%m%d_%H%M%S'))
    out.mkdir(parents=True, exist_ok=False)
    D.dump(out / 'status.json', dict(status='running'))

    print('build shared %g h Han input' % HOURS_LONG, flush=True)
    han = D.HanInput(HOURS_LONG)
    prefix = D.HZHModel('han', han)

    tb, yb = D.integrate(prefix)
    bv, _ = V.analyse(tb, yb, prefix.z, prefix.tail)
    guard = dict(certified=bv['certified_v1'],
                 reads=bv['steady']['reads'] == 19,
                 sequence=bv['steady']['sequence'] == '1234567012345670123',
                 margin=abs(bv['global_min_timing_margin_h']
                            - 1.1784609018563117) < 1e-6)
    print('HZH prefix guard:', json.dumps(guard), flush=True)
    if not all(guard.values()):
        D.dump(out / 'status.json', dict(status='FAILED', guard=guard))
        raise RuntimeError(guard)

    # ---- 1. working-point 300 h, used as the phase source and reference
    m300 = D.HZZModel(prefix, ARM == 'on')
    m300.clock_K = K_WORK
    t300, y300 = D.integrate(m300)
    s300 = m300.signals(y300)
    r300 = evaluate(t300, s300, D.HOURS)
    peaks, _ = V.read_windows(t300, s300)
    pt = t300[peaks]
    period = float(np.median(np.diff(pt)))
    print(f'working point K={K_WORK} {ARM}: {D.HOURS:g} h period={period:.4f} h '
          f'peaks={len(peaks)} seq={r300["steady"]["sequence"]}', flush=True)

    # ---- 2. 600 h long run
    print('600 h run', flush=True)
    m600 = D.HZZModel(prefix, ARM == 'on')
    m600.clock_K = K_WORK
    t600, y600 = integrate_from(m600, m600.initial_state(), HOURS_LONG)
    s600 = m600.signals(y600)
    r600 = evaluate(t600, s600, HOURS_LONG)
    print(f'  600 h: counting={r600["counting_passed"]} '
          f'events={r600["event_causality_passed"]} '
          f'reads={r600["steady"]["reads"]} seq={r600["steady"]["sequence"]}',
          flush=True)

    # ---- 3. eight phase-shifted initial states from one orbit period
    #      taken from the tail of the 300 h run, one full period apart
    t_hi = pt[-2]
    t_lo = t_hi - period
    targets = t_lo + np.arange(N_PHASE) * period / N_PHASE
    idx = [int(np.argmin(np.abs(t300 - x))) for x in targets]
    phase_align_gap = float(np.max(np.abs(t300[idx] - targets)))
    assert phase_align_gap <= 1 / 60 + 1e-9, (phase_align_gap, targets, t300[idx])
    phase_rows = []
    for i, j in enumerate(idx):
        mp = D.HZZModel(prefix, ARM == 'on')
        mp.clock_K = K_WORK
        tp, yp = integrate_from(mp, y300[:, j], D.HOURS)
        sp = mp.signals(yp)
        rp = evaluate(tp, sp, D.HOURS)
        bm = band_margins(rp)
        row = dict(phase_index=i, source_time_h=float(t300[j]),
                   counting=rp['counting_passed'],
                   events=rp['event_causality_passed'],
                   full=rp['full_passed'], sequence=rp['steady']['sequence'],
                   s2_flips=len(rp['crossings']['S2']),
                   S2_band_margin=bm['S2']['band_margin'],
                   S2_unlabelled=bm['S2']['n_unlabelled'],
                   S2_commitment=rp['per_bit']['S2']['minimum_commitment'],
                   stage1=[rp['events']['bit1_to_bit2'][k] for k in
                           ('reverse_events', 'gate_events', 'flip_events')])
        phase_rows.append(row)
        print(f'  phase {i}: t0={row["source_time_h"]:7.2f} '
              f'cnt={int(row["counting"])} ev={int(row["events"])} '
              f'S2fl={row["s2_flips"]} bandM={row["S2_band_margin"]:+.4f} '
              f'seq={row["sequence"]}', flush=True)
        np.savez_compressed(out / f'traj_phase{i}.npz', time_h=tp, states=yp,
                            state_names=np.array(D.NAMES))

    # ---- 4. event-threshold sensitivity, both on the passing and failing points
    print('threshold sensitivity', flush=True)
    mf = D.HZZModel(prefix, ARM == 'on')
    mf.clock_K = K_FAIL
    tf, yf = D.integrate(mf)
    sf = mf.signals(yf)
    sens_pass = sens(t600, y600, m600, HOURS_LONG)
    sens_fail = sens(tf, yf, mf, D.HOURS)
    for name, s in (('pass', sens_pass), ('fail', sens_fail)):
        n_full = sum(x['full'] for x in s)
        print(f'  {name}: {len(s)} settings, {n_full} fully certified', flush=True)

    np.savez_compressed(out / 'traj_work300.npz', time_h=t300, states=y300,
                        state_names=np.array(D.NAMES))
    np.savez_compressed(out / 'traj_work600.npz', time_h=t600, states=y600,
                        state_names=np.array(D.NAMES))
    np.savez_compressed(out / 'traj_fail_K025.npz', time_h=tf, states=yf,
                        state_names=np.array(D.NAMES))
    for nm, (tt, yy, ss, mm) in dict(
            work300=(t300, y300, s300, m300), work600=(t600, y600, s600, m600),
            fail=(tf, yf, sf, mf)).items():
        D.dump(out / f'signals_{nm}.json',
               {k: [float(np.min(v)), float(np.max(v))]
                for k, v in ss.items() if np.ndim(v) == 1})

    D.dump(out / 'round3_summary.json', dict(
        working_point=dict(arm=ARM, clock_K=K_WORK),
        hours_long=HOURS_LONG, hours_short=D.HOURS, n_phase=N_PHASE,
        period_h=period, phase_align_gap_h=phase_align_gap, prefix_guard=guard,
        work300=dict(counting=r300['counting_passed'],
                     events=r300['event_causality_passed'],
                     sequence=r300['steady']['sequence'],
                     reads=r300['steady']['reads']),
        work600=dict(counting=r600['counting_passed'],
                     events=r600['event_causality_passed'],
                     sequence=r600['steady']['sequence'],
                     reads=r600['steady']['reads'],
                     s2_flips=len(r600['crossings']['S2']),
                     stage1=[r600['events']['bit1_to_bit2'][k] for k in
                             ('reverse_events', 'gate_events', 'flip_events')]),
        phases=phase_rows,
        phases_all_pass=bool(all(p['full'] for p in phase_rows)),
        phase_sequences=sorted({p['sequence'] for p in phase_rows}),
        threshold_sensitivity_pass=sens_pass,
        threshold_sensitivity_fail=sens_fail,
        sensitivity_pass_full_count=sum(x['full'] for x in sens_pass),
        sensitivity_fail_full_count=sum(x['full'] for x in sens_fail),
        note=('phase initial states are taken from ONE converged orbit period '
              'at the tail of the 300 h run, so they are genuine orbit phases '
              'rather than random perturbations')))

    with (out / 'phases.csv').open('w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh)
        w.writerow(['phase_index', 'source_time_h', 'counting', 'events', 'full',
                    'sequence', 's2_flips', 'S2_band_margin', 'S2_unlabelled',
                    'S2_commitment', 'stage1_r_g_f'])
        for p in phase_rows:
            w.writerow([p['phase_index'], p['source_time_h'], p['counting'],
                        p['events'], p['full'], p['sequence'], p['s2_flips'],
                        p['S2_band_margin'], p['S2_unlabelled'],
                        p['S2_commitment'], '%d/%d/%d' % tuple(p['stage1'])])

    D.dump(out / 'status.json', dict(
        status='COMPLETED', phases=len(phase_rows),
        phases_all_pass=bool(all(p['full'] for p in phase_rows))))
    D.dump(out / 'SHA256SUMS.json',
           {str(p.relative_to(out)): D.sha256(p)
            for p in sorted(out.rglob('*')) if p.is_file()})
    print('OUTPUT', out, flush=True)


if __name__ == '__main__':
    main()
