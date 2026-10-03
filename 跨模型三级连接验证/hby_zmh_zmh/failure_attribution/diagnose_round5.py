"""Round 5: the two corrections the review demanded.

PART C -- same-time flux statistics (NO re-integration).
  Rounds 1-4 computed the reverse-drive indicator as max(I2) * max(RDF2), the
  product of two maxima attained at DIFFERENT times, and combined it with
  min(pb2), a third time. None of those is a reaction state the model ever
  visited. This part recomputes, from every saved trajectory, the same-time
  quantities:
      max_t ( I2(t) * RDF2(t) )
      max_t v_rev(t),  max_t v_fwd(t),  and integral(v_rev)/integral(v_fwd)
  and checks whether any of them actually separates the passing from the
  failing configurations. If none does, no scalar of this family may be quoted
  as a criterion.

PART D -- corrected eight digital-state continuation (WITH integration).
  The round-3/4 "eight phases" were eight samples inside ONE clock cycle
  (span 9.20 h against a 10.52 h period), all holding the same digital value,
  and each was re-integrated from t = 0 so the upstream oscillator phase
  restarted. This part instead:
    * locates eight consecutive read windows covering values 0-7 in a converged
      orbit, and takes the state at each window's trough;
    * continues each with its OWN upstream time,  dy/dt = f(t + t_start, y),
      by shifting the time argument handed to the Han interpolation;
    * reports the continued sequence, the three-tier verdict and the fluxes.
  The Han upstream is integrated long enough to cover t_start + horizon.
"""
from __future__ import annotations

import csv
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np

sys.dont_write_bytecode = True

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT.parent))
sys.path.insert(0, str(ROOT.parent.parent))

import diagnose_round1 as D                                       # noqa: E402
import verify_hzh as V                                            # noqa: E402
from diagnose_round3 import evaluate, integrate_from              # noqa: E402
from diagnose_round2 import band_margins                          # noqa: E402

RES = ROOT / 'results'
I2, RDF2, PB2 = 20, 22, 19          # NAMES indices in the 23-state layout
KD_COMP = 0.8
KD_INT2 = 1.0
KINH = 0.1

# ---- PART C: every saved trajectory, with its known configuration ----------
def trajectories():
    r1 = sorted(RES.glob('2026*'))
    r1b = sorted(RES.glob('round1b_*'))
    r2 = sorted(RES.glob('round2_*'))
    r4b = sorted(RES.glob('round4b_*'))
    out = []
    if r1:
        base = r1[-1]
        for probe, k in (('base', 0.40), ('krep85', 0.40), ('aint38', 0.40),
                         ('clk020', 0.20), ('clk007', 0.07),
                         ('krep85_aint38', 0.40), ('krep85_clk007', 0.07)):
            for arm in ('off', 'on'):
                p = base / f'traj_{probe}_{arm}.npz'
                if p.exists():
                    out.append((p, arm, k, False))
    if r1b:
        for arm, k in (('off', 0.20), ('off', 0.07)):
            p = r1b[-1] / f'traj_clk{int(round(k*100)):03d}_{arm}.npz'
            if p.exists():
                out.append((p, arm, k, False))
    if r2:
        for k in (0.40, 0.20, 0.25, 0.30, 0.35, 0.05, 0.02):
            p = r2[-1] / f'traj_K{k:.2f}_off.npz'
            if p.exists():
                out.append((p, 'off', k, False))
        for k in (0.25, 0.30, 0.35):
            p = r2[-1] / f'traj_K{k:.2f}_on.npz'
            if p.exists():
                out.append((p, 'on', k, False))
    if r4b:
        for k in (0.40, 0.20, 0.10):
            p = r4b[-1] / f'traj_bypass_K{k:.2f}_off.npz'
            if p.exists():
                out.append((p, 'off', k, True))
        for k in (0.40, 0.10):
            p = r4b[-1] / f'traj_bypass_K{k:.2f}_on.npz'
            if p.exists():
                out.append((p, 'on', k, True))
    return out


def fluxes(y, t):
    i = np.asarray(y[I2], dtype=float)
    r = np.asarray(y[RDF2], dtype=float)
    pb = np.asarray(y[PB2], dtype=float)
    ir = i * r
    j = int(np.argmax(ir))
    vr = 7.0 * (1 - pb) * ir ** 2 / (KD_COMP ** 2 + ir ** 2)
    vf = (7.0 * pb * i ** 2 / (KD_INT2 ** 2 + i ** 2)
          * KINH / (KINH + r))
    return dict(
        max_I2=float(i.max()), max_RDF2=float(r.max()),
        max_I2_times_max_RDF2=float(i.max() * r.max()),
        max_same_time_IR=float(ir[j]),
        max_same_time_IR_at_h=float(t[j]),
        max_same_time_IR_over_threshold=float(ir[j] / KD_COMP),
        max_pb2_used_before=float(pb.min()),
        max_v_rev=float(vr.max()), max_v_fwd=float(vf.max()),
        v_rev_over_v_fwd_at_max=float(vr.max() / vf.max()),
        int_v_rev=float(np.trapezoid(np.maximum(vr, 0), t)),
        int_v_fwd=float(np.trapezoid(np.maximum(vf, 0), t)),
        int_ratio=float(np.trapezoid(np.maximum(vr, 0), t)
                        / max(np.trapezoid(np.maximum(vf, 0), t), 1e-30)))


def main():
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8')
    out = ROOT / 'results' / ('round5_' + datetime.now().strftime('%Y%m%d_%H%M%S'))
    out.mkdir(parents=True, exist_ok=False)
    D.dump(out / 'status.json', dict(status='running'))

    print('build Han upstream', flush=True)
    han = D.HanInput(D.HOURS)
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

    # =====================================================================
    # PART C
    # =====================================================================
    print('\n--- PART C: same-time flux statistics (no re-integration) ---',
          flush=True)
    rows = []
    for path, arm, k, bypass in trajectories():
        z = np.load(path)
        t, y = z['time_h'], z['states']
        m = D.HZZModel(prefix, arm == 'on')
        m.clock_K = k
        if bypass:
            from diagnose_round4 import BypassClock
            m = BypassClock(prefix, arm == 'on')
            m.clock_K = k
        sig = m.signals(y)
        res = evaluate(t, sig, float(t[-1]))
        r = dict(source=f'{path.parent.name}/{path.name}', arm=arm, clock_K=k,
                 bypass=bypass, counting=res['counting_passed'],
                 events=res['event_causality_passed'])
        r.update(fluxes(y, t))
        rows.append(r)
        flag = 'PASS' if (r['counting'] and r['events']) else 'fail'
        print(f"  {r['source'][:44]:44s} K={k:.2f} {'byp' if bypass else '   '} "
              f"{flag}  peakxpeak={r['max_I2_times_max_RDF2']:8.4f} "
              f"sameT={r['max_same_time_IR']:8.4f} "
              f"({r['max_same_time_IR_over_threshold']:6.2f}) "
              f"vr/vf(max)={r['v_rev_over_v_fwd_at_max']:9.2f} "
              f"int_ratio={r['int_ratio']:.4f}", flush=True)

    passing = [r for r in rows if r['counting'] and r['events']]
    failing = [r for r in rows if not (r['counting'] and r['events'])]
    sep = {}
    for key in ('max_same_time_IR_over_threshold', 'v_rev_over_v_fwd_at_max',
                'int_ratio'):
        pv = [r[key] for r in passing]
        fv = [r[key] for r in failing]
        overlap = bool(pv and fv and min(pv) <= max(fv))
        sep[key] = dict(passing_min=min(pv) if pv else None,
                        passing_max=max(pv) if pv else None,
                        failing_min=min(fv) if fv else None,
                        failing_max=max(fv) if fv else None,
                        separates=not overlap,
                        overlap_value=(max(fv) if overlap and fv else None))
        print(f"  SEPARATION {key}: pass[{sep[key]['passing_min']:.4f},"
              f"{sep[key]['passing_max']:.4f}] fail[{sep[key]['failing_min']:.4f},"
              f"{sep[key]['failing_max']:.4f}] separates={sep[key]['separates']}",
              flush=True)
    if passing:
        raw = [r['max_I2_times_max_RDF2'] for r in passing]
        raw_ex = [r['max_same_time_IR'] for r in passing]
        print(f"  inflation of peak-x-peak over same-time on passing runs: "
              f"max {max(a/b for a, b in zip(raw, raw_ex)):.2f}x", flush=True)

    D.dump(out / 'partC_flux.json', dict(
        note=('peak_x_peak multiplies maxima attained at different times and is '
              'NOT a reaction state; same_time is max_t I2(t)*RDF2(t)'),
        n_passing=len(passing), n_failing=len(failing),
        separation=sep, rows=rows))

    # =====================================================================
    # PART D
    # =====================================================================
    print('\n--- PART D: eight digital-state continuations with upstream time ---',
          flush=True)
    K_WORK, HOURS_SRC, HOURS_CONT = 0.10, 300.0, 150.0
    han_long = D.HanInput(HOURS_SRC + HOURS_CONT + 20.0)
    prefix_long = D.HZHModel('han', han_long)

    ms = D.HZZModel(prefix_long, False)
    ms.clock_K = K_WORK
    ts, ys = D.integrate(ms)
    ss = ms.signals(ys)
    reads = V.read_windows(ts, ss)[1]
    steady = reads[V.DROP:]
    # find a consecutive run of eight windows covering all values 0-7
    start = None
    for a in range(len(steady) - 8 + 1):
        vals = [steady[a + i]['value'] for i in range(8)]
        if None not in vals and sorted(vals) == list(range(8)):
            start = a
            break
    if start is None:
        raise RuntimeError('no consecutive eight-window run covering 0-7')
    chosen = steady[start:start + 8]
    print('  source run: ' + ''.join(
        'x' if r['value'] is None else str(r['value']) for r in steady))
    print('  chosen windows: ' +
          ', '.join(f"v={r['value']}@t={r['trough_h']:.2f}h" for r in chosen),
          flush=True)

    class TimeShifted(D.HZZModel):
        """Continue a snapshot with its own upstream time: dy/dt = f(t+t0, y).

        HZZModel.rhs forwards t to prefix.bit0_rhs, which interpolates the Han
        upstream at 60*t. Shifting the argument keeps the oscillator phase that
        the snapshot was taken at, instead of restarting it at zero.
        """
        def __init__(self, prefix, autoregulation, t_shift, clock_K):
            super().__init__(prefix, autoregulation)
            self.t_shift = float(t_shift)
            self.clock_K = clock_K

        def rhs(self, t, y):
            return super().rhs(t + self.t_shift, y)

    cont = []
    for r in chosen:
        t0 = float(r['trough_h'])
        j = int(np.argmin(np.abs(ts - t0)))
        snap = ys[:, j]
        mc = TimeShifted(prefix_long, False, t0, K_WORK)
        tc, yc = integrate_from(mc, snap, HOURS_CONT)
        sc = mc.signals(yc)
        rc = evaluate(tc, sc, HOURS_CONT)
        bm = band_margins(rc)
        # value actually read from the snapshot itself
        s0, s1, s2 = float(sc['S0'][0]), float(sc['S1'][0]), float(sc['S2'][0])
        row = dict(target_value=r['value'], t_start_h=t0, horizon_h=HOURS_CONT,
                   snapshot_S0=s0, snapshot_S1=s1, snapshot_S2=s2,
                   counting=rc['counting_passed'],
                   events=rc['event_causality_passed'],
                   full=rc['full_passed'],
                   steady_sequence=rc['steady']['sequence'],
                   increments_mod8=rc['steady']['increments_mod8'],
                   reads=rc['steady']['reads'],
                   s2_flips=len(rc['crossings']['S2']),
                   stage1=[rc['events']['bit1_to_bit2'][k] for k in
                           ('reverse_events', 'gate_events', 'flip_events')],
                   S2_band_margin=bm['S2']['band_margin'],
                   S2_unlabelled=bm['S2']['n_unlabelled'])
        cont.append(row)
        print(f"  start value {r['value']}: t0={t0:7.2f} h  "
              f"snap S0/S1/S2={s0:.3f}/{s1:.3f}/{s2:.3f}  "
              f"cnt={int(row['counting'])} ev={int(row['events'])} "
              f"inc8={int(row['increments_mod8'])} reads={row['reads']} "
              f"stage1={row['stage1']} seq={row['steady_sequence']}", flush=True)
        np.savez_compressed(out / f'traj_digit{r["value"]}.npz',
                            time_h=tc, states=yc, state_names=np.array(D.NAMES))

    n_cnt = sum(r['counting'] for r in cont)
    n_ev = sum(r['events'] for r in cont)
    n_inc = sum(r['increments_mod8'] for r in cont)
    print(f'  -> counting {n_cnt}/8, events {n_ev}/8, increments_mod8 {n_inc}/8',
          flush=True)

    D.dump(out / 'partD_eight_states.json', dict(
        clock_K=K_WORK, source_hours=HOURS_SRC, continuation_hours=HOURS_CONT,
        upstream_hours=HOURS_SRC + HOURS_CONT + 20.0,
        construction=('eight consecutive read windows covering values 0-7 in a '
                      'converged orbit; snapshot taken at each window trough; '
                      'continued with its own upstream time via dy/dt = '
                      'f(t + t_start, y)'),
        chosen=[dict(value=r['value'], trough_h=r['trough_h']) for r in chosen],
        continuations=cont,
        counting_passed=n_cnt, event_passed=n_ev, increments_mod8=n_inc,
        note=('supersedes the round-3/4 "eight phases", which sampled 9.20 h '
              'inside a single 10.52 h clock cycle (all one digital value) and '
              'restarted the upstream at t=0'),
        prefix_guard=guard))

    with (out / 'partC_flux.csv').open('w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh)
        w.writerow(['source', 'arm', 'clock_K', 'bypass', 'counting', 'events',
                    'max_I2', 'max_RDF2', 'max_I2_times_max_RDF2',
                    'max_same_time_IR', 'max_same_time_IR_over_threshold',
                    'max_v_rev', 'max_v_fwd', 'v_rev_over_v_fwd_at_max',
                    'int_v_rev', 'int_v_fwd', 'int_ratio'])
        for r in rows:
            w.writerow([r['source'], r['arm'], r['clock_K'], r['bypass'],
                        r['counting'], r['events'], r['max_I2'],
                        r['max_RDF2'], r['max_I2_times_max_RDF2'],
                        r['max_same_time_IR'],
                        r['max_same_time_IR_over_threshold'], r['max_v_rev'],
                        r['max_v_fwd'], r['v_rev_over_v_fwd_at_max'],
                        r['int_v_rev'], r['int_v_fwd'], r['int_ratio']])
    with (out / 'partD_eight_states.csv').open('w', newline='',
                                               encoding='utf-8') as fh:
        w = csv.writer(fh)
        w.writerow(['value', 't_start_h', 'snapshot_S0', 'snapshot_S1',
                    'snapshot_S2', 'counting', 'events', 'full',
                    'increments_mod8', 'reads', 's2_flips', 'stage1_r_g_f',
                    'S2_band_margin', 'S2_unlabelled', 'steady_sequence'])
        for r in cont:
            w.writerow([r['target_value'], r['t_start_h'], r['snapshot_S0'],
                        r['snapshot_S1'], r['snapshot_S2'], r['counting'],
                        r['events'], r['full'], r['increments_mod8'],
                        r['reads'], r['s2_flips'],
                        '%d/%d/%d' % tuple(r['stage1']), r['S2_band_margin'],
                        r['S2_unlabelled'], r['steady_sequence']])

    D.dump(out / 'status.json', dict(status='COMPLETED',
                                     n_flux_rows=len(rows),
                                     eight_states_counting=n_cnt,
                                     eight_states_events=n_ev,
                                     eight_states_increments=n_inc))
    D.dump(out / 'SHA256SUMS.json',
           {str(p.relative_to(out)): D.sha256(p)
            for p in sorted(out.rglob('*')) if p.is_file()})
    print('OUTPUT', out, flush=True)


if __name__ == '__main__':
    main()
