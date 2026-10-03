"""Round 5b: corrected eight digital-state continuation, with a usable horizon.

Round 5 PART D built the eight starting states correctly -- the snapshot values
covered 1..7 and 0 with clean band separation, and every continued sequence
incremented mod 8 from its own start -- but the continuation was only 150 h.
That yields about 5 read windows, while the criterion needs DROP = 8 plus
MIN_STEADY = 16, i.e. at least 24 windows, i.e. about 252 h at a 10.52 h period.
So the run reported 0/8 purely for lack of windows.

This script repeats PART D with a 300 h continuation and an upstream integrated
far enough to cover t_start + 300 h. Nothing else changes; the construction is
identical to round 5 PART D.
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

K_WORK = 0.10
HOURS_SRC = 300.0
HOURS_CONT = 300.0
UPSTREAM_HOURS = HOURS_SRC + HOURS_CONT + 20.0


class TimeShifted(D.HZZModel):
    """dy/dt = f(t + t_start, y): keeps the upstream phase of the snapshot."""

    def __init__(self, prefix, autoregulation, t_shift, clock_K):
        super().__init__(prefix, autoregulation)
        self.t_shift = float(t_shift)
        self.clock_K = clock_K

    def rhs(self, t, y):
        return super().rhs(t + self.t_shift, y)


def main():
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8')
    out = ROOT / 'results' / ('round5b_' + datetime.now().strftime('%Y%m%d_%H%M%S'))
    out.mkdir(parents=True, exist_ok=False)
    D.dump(out / 'status.json', dict(status='running'))

    print(f'build Han upstream to {UPSTREAM_HOURS:g} h', flush=True)
    han = D.HanInput(UPSTREAM_HOURS)
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

    ms = D.HZZModel(prefix, False)
    ms.clock_K = K_WORK
    ts, ys = D.integrate(ms)
    ss = ms.signals(ys)
    reads = V.read_windows(ts, ss)[1]
    steady = reads[V.DROP:]
    start = None
    for a in range(len(steady) - 8 + 1):
        vals = [steady[a + i]['value'] for i in range(8)]
        if None not in vals and sorted(vals) == list(range(8)):
            start = a
            break
    if start is None:
        raise RuntimeError('no consecutive eight-window run covering 0-7')
    chosen = steady[start:start + 8]
    print('  source steady sequence: ' + ''.join(
        'x' if r['value'] is None else str(r['value']) for r in steady))
    print('  chosen: ' + ', '.join(
        f"v={r['value']}@{r['trough_h']:.2f}h" for r in chosen), flush=True)

    cont = []
    for r in chosen:
        t0 = float(r['trough_h'])
        j = int(np.argmin(np.abs(ts - t0)))
        mc = TimeShifted(prefix, False, t0, K_WORK)
        tc, yc = integrate_from(mc, ys[:, j], HOURS_CONT)
        sc = mc.signals(yc)
        rc = evaluate(tc, sc, HOURS_CONT)
        bm = band_margins(rc)
        s0, s1, s2 = float(sc['S0'][0]), float(sc['S1'][0]), float(sc['S2'][0])
        row = dict(target_value=r['value'], t_start_h=t0,
                   horizon_h=HOURS_CONT,
                   snapshot_S0=s0, snapshot_S1=s1, snapshot_S2=s2,
                   counting=rc['counting_passed'],
                   events=rc['event_causality_passed'],
                   full=rc['full_passed'],
                   all_reads=len(rc['read_windows']),
                   steady_reads=rc['steady']['reads'],
                   steady_sequence=rc['steady']['sequence'],
                   increments_mod8=rc['steady']['increments_mod8'],
                   s2_flips=len(rc['crossings']['S2']),
                   stage1=[rc['events']['bit1_to_bit2'][k] for k in
                           ('reverse_events', 'gate_events', 'flip_events')],
                   S2_band_margin=bm['S2']['band_margin'],
                   S2_unlabelled=bm['S2']['n_unlabelled'],
                   S2_commitment=rc['per_bit']['S2']['minimum_commitment'])
        cont.append(row)
        np.savez_compressed(out / f'traj_digit{r["value"]}.npz',
                            time_h=tc, states=yc, state_names=np.array(D.NAMES))
        print(f"  start {r['value']}: t0={t0:7.2f}  snap={s0:.3f}/{s1:.3f}/{s2:.3f}  "
              f"cnt={int(row['counting'])} ev={int(row['events'])} "
              f"reads={row['steady_reads']}/{row['all_reads']} "
              f"inc8={int(row['increments_mod8'])} "
              f"stage1={row['stage1']} seq={row['steady_sequence']}", flush=True)

    n_cnt = sum(r['counting'] for r in cont)
    n_ev = sum(r['events'] for r in cont)
    n_inc = sum(r['increments_mod8'] for r in cont)
    n_full = sum(r['full'] for r in cont)
    starts = [(r['target_value'], r['steady_sequence'][:1]) for r in cont]
    print(f'  -> counting {n_cnt}/8, events {n_ev}/8, increments_mod8 {n_inc}/8, '
          f'full {n_full}/8', flush=True)
    print(f'  first steady value vs requested start: {starts}', flush=True)

    D.dump(out / 'eight_states.json', dict(
        clock_K=K_WORK, source_hours=HOURS_SRC, continuation_hours=HOURS_CONT,
        upstream_hours=UPSTREAM_HOURS, prefix_guard=guard,
        construction=('eight consecutive read windows covering values 0-7 in a '
                      'converged 300 h orbit; snapshot at each window trough; '
                      'continued with dy/dt = f(t + t_start, y) so the upstream '
                      'oscillator keeps its phase'),
        chosen=[dict(value=r['value'], trough_h=r['trough_h']) for r in chosen],
        continuations=cont, counting_passed=n_cnt, event_passed=n_ev,
        increments_mod8=n_inc, full_passed=n_full,
        supersedes=('round 5 PART D (150 h: only ~5 windows, so DROP=8 + '
                    'MIN_STEADY=16 could not be satisfied) and the round-3/4 '
                    '"eight phases" (9.20 h inside one 10.52 h clock cycle, all '
                    'one digital value, upstream restarted at t=0)')))
    with (out / 'eight_states.csv').open('w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh)
        w.writerow(['value', 't_start_h', 'snapshot_S0', 'snapshot_S1',
                    'snapshot_S2', 'counting', 'events', 'full', 'all_reads',
                    'steady_reads', 'increments_mod8', 's2_flips',
                    'stage1_r_g_f', 'S2_band_margin', 'S2_unlabelled',
                    'S2_commitment', 'steady_sequence'])
        for r in cont:
            w.writerow([r['target_value'], r['t_start_h'], r['snapshot_S0'],
                        r['snapshot_S1'], r['snapshot_S2'], r['counting'],
                        r['events'], r['full'], r['all_reads'],
                        r['steady_reads'], r['increments_mod8'], r['s2_flips'],
                        '%d/%d/%d' % tuple(r['stage1']), r['S2_band_margin'],
                        r['S2_unlabelled'], r['S2_commitment'],
                        r['steady_sequence']])
    D.dump(out / 'status.json', dict(status='COMPLETED', counting=n_cnt,
                                     events=n_ev, increments=n_inc,
                                     full=n_full))
    D.dump(out / 'SHA256SUMS.json',
           {str(p.relative_to(out)): D.sha256(p)
            for p in sorted(out.rglob('*')) if p.is_file()})
    print('OUTPUT', out, flush=True)


if __name__ == '__main__':
    main()
