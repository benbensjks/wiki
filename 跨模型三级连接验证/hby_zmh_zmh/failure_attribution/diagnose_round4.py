"""Round 4: two open questions from round 3.

PART A -- is the 7/8 phase result a horizon effect?
  Round 3 took 8 orbit-phase initial states at clock_K = 0.10 and ran each for
  300 h. Counting passed in 8/8 but the event-causality arm passed in only 7/8:
  phases 0 and 7 both had 7 reverse / 6 gate / 6 flip instead of 7/7/7, and only
  phase 7 failed the causal test. The sequences were still valid mod-8 strings,
  so the suspicion is that one fewer complete carry fits inside a 300 h window
  depending on where the phase starts. This part re-runs ALL eight phases for
  400 h and additionally records how many clock peaks (cycles) fall inside the
  window and how far the last carry event sits from the window end.

PART B -- is pulse_gate a cycle marker that must not act as a discriminator?
  The source file 前馈三级级联.py drives int0 with a square wave whose plateau is
  k_int/gamma_int = 3.0, so its pulse_gate = int0^3/(0.4^3+int0^3) is ~0.998
  during a pulse and ~0 between pulses: a cycle marker, not a threshold test.
  Round 2 showed the count only works while clock_K stays low enough that this
  term is NOT the limiting factor. This part forces the clock term to 1 and
  re-runs the original clock_K = 0.40. If the count then works, the term is
  indeed doing harm as a discriminator rather than being needed as a gate.

Nothing outside failure_attribution/ is written.

INVALIDATION NOTE
-----------------
The FIRST run of this script (results/round4_20261001_165636/) has a usable
PART A but an INVALID PART B. BypassClock.tail_signals wrote q['gate'] while
HZZModel.tail_signals publishes the gate under the key 'g1', so the parent's
clock-multiplied g1 survived into the reported signals. The dynamics were
correct (int2_source WAS overridden, which is why the count was restored and
RDF2 jumped 35x), but the event arm was judged on the un-bypassed g1 and
therefore reported a spurious gate=0 failure at clock_K = 0.40.
PART B was re-run by diagnose_round4b.py with the key fixed. Only that re-run
may be cited for the event arm.
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
from diagnose_round3 import evaluate, integrate_from              # noqa: E402

K_WORK = 0.10
HOURS_PHASE = 400.0
N_PHASE = 8
BYPASS_CASES = [(0.40, 'off'), (0.20, 'off'), (0.10, 'off'),
                (0.40, 'on'), (0.10, 'on')]


class BypassClock(D.HZZModel):
    """clock term forced to 1, i.e. pulse_gate removed.

    Everything else -- act, repress, the Int2 source, the whole recombination
    block -- is left exactly as in the parent class. Only the third factor of
    the gate is replaced.

    NOTE ON THE DICT KEY: HZZModel.tail_signals returns the gate under the key
    'g1', NOT 'gate'. An earlier version of this class wrote q['gate'], which
    silently created a new key and left q['g1'] carrying the clock-multiplied
    value from the parent. The dynamics still used the bypassed int2_source, so
    the COUNT was right, but the REPORTED g1 was the un-bypassed one and the
    event arm was judged on the wrong signal. Both keys are now set.
    """

    def tail_signals(self, y, int0):
        q = super().tail_signals(y, int0)
        arr = np.asarray(int0, dtype=float)
        clk = np.ones_like(arr) if arr.ndim else 1.0
        gate = q['activation'] * q['repression']
        q['clock'] = clk
        q['g1'] = gate
        q['gate'] = gate
        q['int2_source'] = self.p['alpha_Int2'] * gate
        return q


def main():
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8')
    out = ROOT / 'results' / ('round4_' + datetime.now().strftime('%Y%m%d_%H%M%S'))
    out.mkdir(parents=True, exist_ok=False)
    D.dump(out / 'status.json', dict(status='running'))

    print('build shared %g h Han input' % HOURS_PHASE, flush=True)
    han = D.HanInput(HOURS_PHASE)
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

    # ---- phase source: identical construction to round 3
    m300 = D.HZZModel(prefix, False)
    m300.clock_K = K_WORK
    t300, y300 = D.integrate(m300)
    s300 = m300.signals(y300)
    peaks, _ = V.read_windows(t300, s300)
    pt = t300[peaks]
    period = float(np.median(np.diff(pt)))
    targets = (pt[-2] - period) + np.arange(N_PHASE) * period / N_PHASE
    idx = [int(np.argmin(np.abs(t300 - x))) for x in targets]
    align = float(np.max(np.abs(t300[idx] - targets)))
    assert align <= 1 / 60 + 1e-9, align
    print(f'period={period:.4f} h  phases aligned within {align*60:.3f} min',
          flush=True)

    # =====================================================================
    # PART A
    # =====================================================================
    print('\n--- PART A: eight phases at %g h ---' % HOURS_PHASE, flush=True)
    ph = []
    for i, j in enumerate(idx):
        mp = D.HZZModel(prefix, False)
        mp.clock_K = K_WORK
        tp, yp = integrate_from(mp, y300[:, j], HOURS_PHASE)
        sp = mp.signals(yp)
        rp = evaluate(tp, sp, HOURS_PHASE)
        bm = band_margins(rp)
        pk, _ = V.read_windows(tp, sp)
        s1 = rp['events']['bit1_to_bit2']
        last_rev = max([e['peak_h'] for e in
                        V.segments(tp, sp['J_rev1'], V.JREV_THRESHOLD_H)],
                       default=None)
        row = dict(
            phase=i, source_time_h=float(t300[j]),
            counting=rp['counting_passed'], events=rp['event_causality_passed'],
            full=rp['full_passed'], sequence=rp['steady']['sequence'],
            reads=rp['steady']['reads'],
            s2_flips=len(rp['crossings']['S2']),
            stage1=[s1['reverse_events'], s1['gate_events'], s1['flip_events']],
            stage0=[rp['events']['bit0_to_bit1'][k] for k in
                    ('reverse_events', 'gate_events', 'flip_events')],
            clock_peaks=len(pk),
            last_rev_h=last_rev,
            gap_to_end_h=(None if last_rev is None else
                          float(HOURS_PHASE - last_rev)),
            S2_band_margin=bm['S2']['band_margin'],
            S2_unlabelled=bm['S2']['n_unlabelled'],
            S2_commitment=rp['per_bit']['S2']['minimum_commitment'])
        ph.append(row)
        np.savez_compressed(out / f'traj_p{i}_{int(HOURS_PHASE)}h.npz',
                            time_h=tp, states=yp,
                            state_names=np.array(D.NAMES))
        gap = row['gap_to_end_h']
        gap_str = 'None' if gap is None else f'{gap:.2f}'
        print(f"  phase {i}: cnt={int(row['counting'])} "
              f"ev={int(row['events'])} reads={row['reads']} "
              f"peaks={row['clock_peaks']} "
              f"stage1={row['stage1']} s2fl={row['s2_flips']} "
              f"gapToEnd={gap_str} h "
              f"seq={row['sequence']}", flush=True)

    all_cnt = all(p['counting'] for p in ph)
    all_ev = all(p['events'] for p in ph)
    print(f'  -> counting {sum(p["counting"] for p in ph)}/{len(ph)}, '
          f'event {sum(p["events"] for p in ph)}/{len(ph)}', flush=True)

    # =====================================================================
    # PART B
    # =====================================================================
    print('\n--- PART B: pulse_gate bypassed (clock term = 1) ---', flush=True)
    byp = []
    for clock_K, arm in BYPASS_CASES:
        mb = BypassClock(prefix, arm == 'on')
        mb.clock_K = clock_K
        t, y = D.integrate(mb)
        sig = mb.signals(y)
        sig['_p'] = mb.p
        res = evaluate(t, sig, D.HOURS)
        bm = band_margins(res)
        s1 = res['events']['bit1_to_bit2']
        segs = D.raw_gate_segments(t, sig['g1'])
        mech = D.mechanism(y, mb.p, clock_K)
        dose = [s['dose'] for s in segs] or [None]
        dur = [s['duration_h'] for s in segs] or [None]
        row = dict(arm=arm, clock_K=clock_K, bypass=True,
                   counting=res['counting_passed'],
                   events=res['event_causality_passed'],
                   full=res['full_passed'],
                   sequence=res['steady']['sequence'],
                   reads=res['steady']['reads'],
                   s2_flips=len(res['crossings']['S2']),
                   stage1=[s1['reverse_events'], s1['gate_events'],
                           s1['flip_events']],
                   gate_segments=len(segs),
                   gate_dose=[min(dose), max(dose)],
                   gate_duration=[min(dur), max(dur)],
                   meets_both=sum(s['meets_duration'] and s['meets_dose']
                                  for s in segs),
                   g1_peak=float(np.max(sig['g1'])),
                   S2_band_margin=bm['S2']['band_margin'],
                   S2_unlabelled=bm['S2']['n_unlabelled'],
                   mechanism=mech)
        byp.append(row)
        np.savez_compressed(out / f'traj_bypass_K{clock_K:.2f}_{arm}.npz',
                            time_h=t, states=y, state_names=np.array(D.NAMES))
        print(f"  bypass K={clock_K:.2f} {arm:3s}: cnt={int(row['counting'])} "
              f"ev={int(row['events'])} stage1={row['stage1']} "
              f"dose=[{dose[0]:.4f},{dose[-1]:.4f}] "
              f"meetBoth={row['meets_both']} "
              f"RDF2pk={mech['RDF2_peak']:.4f} "
              f"vr/vf={mech['reverse_over_forward']:.1f} "
              f"seq={row['sequence']}", flush=True)

    D.dump(out / 'round4_summary.json', dict(
        partA=dict(hours=HOURS_PHASE, n_phase=N_PHASE, period_h=period,
                   align_min=align * 60, clock_K=K_WORK, phases=ph,
                   all_counting=all_cnt, all_events=all_ev,
                   reading=('if every phase now passes, the round-3 7/8 was a '
                            'window-boundary effect; if phase 7 still fails, the '
                            'result is phase-dependent')),
        partB=dict(question=('does forcing the clock term to 1 (removing '
                             'pulse_gate) restore the count at the original '
                             'clock_K = 0.40?'),
                   cases=byp,
                   restored_at_original=[r for r in byp
                                         if r['clock_K'] == 0.40
                                         and r['counting'] and r['events']]),
        prefix_guard=guard))

    for name, rows_, cols in (
            ('phases_400h.csv', ph,
             ['phase', 'source_time_h', 'counting', 'events', 'full', 'reads',
              'clock_peaks', 'stage0_r_g_f', 'stage1_r_g_f', 's2_flips',
              'last_rev_h', 'gap_to_end_h', 'S2_band_margin', 'S2_unlabelled',
              'S2_commitment', 'sequence']),
            ('bypass.csv', byp,
             ['arm', 'clock_K', 'counting', 'events', 'full', 'reads',
              's2_flips', 'stage1_r_g_f', 'gate_segments', 'dose_min',
              'dose_max', 'dur_min', 'dur_max', 'meets_both', 'g1_peak',
              'S2_band_margin', 'S2_unlabelled', 'I2_peak', 'RDF2_peak',
              'reverse_hill', 'vr_over_vf', 'sequence'])):
        with (out / name).open('w', newline='', encoding='utf-8') as fh:
            w = csv.writer(fh)
            w.writerow(cols)
            for r in rows_:
                if name.startswith('phases'):
                    w.writerow([r['phase'], r['source_time_h'], r['counting'],
                                r['events'], r['full'], r['reads'],
                                r['clock_peaks'], '%d/%d/%d' % tuple(r['stage0']),
                                '%d/%d/%d' % tuple(r['stage1']), r['s2_flips'],
                                r['last_rev_h'], r['gap_to_end_h'],
                                r['S2_band_margin'], r['S2_unlabelled'],
                                r['S2_commitment'], r['sequence']])
                else:
                    m = r['mechanism']
                    w.writerow([r['arm'], r['clock_K'], r['counting'],
                                r['events'], r['full'], r['reads'],
                                r['s2_flips'], '%d/%d/%d' % tuple(r['stage1']),
                                r['gate_segments'], r['gate_dose'][0],
                                r['gate_dose'][1], r['gate_duration'][0],
                                r['gate_duration'][1], r['meets_both'],
                                r['g1_peak'], r['S2_band_margin'],
                                r['S2_unlabelled'], m['I2_peak'],
                                m['RDF2_peak'], m['reverse_hill_value'],
                                m['reverse_over_forward'], r['sequence']])

    D.dump(out / 'status.json', dict(status='COMPLETED',
                                     partA_all_counting=all_cnt,
                                     partA_all_events=all_ev,
                                     partB_restored=len(
                                         [r for r in byp if r['clock_K'] == 0.40
                                          and r['counting'] and r['events']])))
    D.dump(out / 'SHA256SUMS.json',
           {str(p.relative_to(out)): D.sha256(p)
            for p in sorted(out.rglob('*')) if p.is_file()})
    print('OUTPUT', out, flush=True)


if __name__ == '__main__':
    main()
