"""Round 4b: re-run PART B of round 4 with the signal-key bug fixed.

Round 4's BypassClock wrote q['gate'], but HZZModel.tail_signals publishes the
gate under the key 'g1'. The parent's clock-multiplied g1 therefore survived
into the reported signals, so round 4's PART B judged the event arm on the
un-bypassed gate and reported a spurious gate = 0. The dynamics were correct;
only the reported signal was wrong.

This script imports the FIXED BypassClock from diagnose_round4 and re-runs the
same five cases, so the two runs can be compared directly. Round 4's PART A
(eight phases at 400 h) is unaffected and is not repeated.
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
from diagnose_round2 import band_margins                          # noqa: E402
from diagnose_round3 import evaluate                              # noqa: E402
import diagnose_round4 as R4                                      # noqa: E402


def main():
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8')
    out = ROOT / 'results' / ('round4b_' + datetime.now().strftime('%Y%m%d_%H%M%S'))
    out.mkdir(parents=True, exist_ok=False)
    D.dump(out / 'status.json', dict(status='running'))

    print('build shared %g h Han input' % D.HOURS, flush=True)
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

    rows = []
    for clock_K, arm in R4.BYPASS_CASES:
        mb = R4.BypassClock(prefix, arm == 'on')
        mb.clock_K = clock_K
        t, y = D.integrate(mb)
        sig = mb.signals(y)
        sig['_p'] = mb.p
        res = evaluate(t, sig, D.HOURS)
        bm = band_margins(res)
        s1 = res['events']['bit1_to_bit2']
        s0 = res['events']['bit0_to_bit1']
        segs = D.raw_gate_segments(t, sig['g1'])
        mech = D.mechanism(y, mb.p, clock_K)
        dose = [s['dose'] for s in segs] or [None]
        dur = [s['duration_h'] for s in segs] or [None]
        row = dict(
            arm=arm, clock_K=clock_K, bypass=True,
            counting=res['counting_passed'],
            events=res['event_causality_passed'],
            full=res['full_passed'],
            sequence=res['steady']['sequence'], reads=res['steady']['reads'],
            s2_flips=len(res['crossings']['S2']),
            stage0=[s0['reverse_events'], s0['gate_events'], s0['flip_events']],
            stage1=[s1['reverse_events'], s1['gate_events'], s1['flip_events']],
            gate_segments=len(segs),
            dose_min=dose[0], dose_max=dose[-1],
            dur_min=dur[0], dur_max=dur[-1],
            meets_both=sum(s['meets_duration'] and s['meets_dose'] for s in segs),
            g1_peak=float(np.max(sig['g1'])),
            clock_signal_max=float(np.max(sig['clock'])),
            S2_band_margin=bm['S2']['band_margin'],
            S2_unlabelled=bm['S2']['n_unlabelled'],
            S2_commitment=res['per_bit']['S2']['minimum_commitment'],
            mechanism=mech)
        rows.append(row)
        np.savez_compressed(out / f'traj_bypass_K{clock_K:.2f}_{arm}.npz',
                            time_h=t, states=y, state_names=np.array(D.NAMES))
        print(f"  bypass K={clock_K:.2f} {arm:3s}: "
              f"cnt={int(row['counting'])} ev={int(row['events'])} "
              f"stage1={row['stage1']} dose=[{dose[0]:.4f},{dose[-1]:.4f}] "
              f"meetBoth={row['meets_both']} g1pk={row['g1_peak']:.4f} "
              f"clockmax={row['clock_signal_max']:.3f} "
              f"RDF2pk={mech['RDF2_peak']:.4f} "
              f"vr/vf={mech['reverse_over_forward']:.1f} "
              f"S2bandM={bm['S2']['band_margin']:+.4f} "
              f"seq={row['sequence']}", flush=True)

    orig = [r for r in rows if r['clock_K'] == 0.40]
    D.dump(out / 'round4b_summary.json', dict(
        question=('with pulse_gate removed (clock term forced to 1), does the '
                  'source file original clock_K = 0.40 work?'),
        reference_round4_invalid_partB=('round4_20261001_165636 PART B reported '
                                        'gate=0 at K=0.40; that was the signal-key '
                                        'bug, not a real failure'),
        cases=rows,
        restored_at_original=orig,
        all_restored=bool(all(r['counting'] and r['events'] for r in orig)),
        prefix_guard=guard))

    with (out / 'bypass.csv').open('w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh)
        w.writerow(['arm', 'clock_K', 'counting', 'events', 'full', 'reads',
                    's2_flips', 'stage0_r_g_f', 'stage1_r_g_f', 'gate_segments',
                    'dose_min', 'dose_max', 'dur_min', 'dur_max', 'meets_both',
                    'g1_peak', 'clock_signal_max', 'S2_band_margin',
                    'S2_unlabelled', 'S2_commitment', 'I2_peak', 'RDF2_peak',
                    'reverse_hill', 'vr_over_vf', 'sequence'])
        for r in rows:
            m = r['mechanism']
            w.writerow([r['arm'], r['clock_K'], r['counting'], r['events'],
                        r['full'], r['reads'], r['s2_flips'],
                        '%d/%d/%d' % tuple(r['stage0']),
                        '%d/%d/%d' % tuple(r['stage1']), r['gate_segments'],
                        r['dose_min'], r['dose_max'], r['dur_min'], r['dur_max'],
                        r['meets_both'], r['g1_peak'], r['clock_signal_max'],
                        r['S2_band_margin'], r['S2_unlabelled'],
                        r['S2_commitment'], m['I2_peak'], m['RDF2_peak'],
                        m['reverse_hill_value'], m['reverse_over_forward'],
                        r['sequence']])

    D.dump(out / 'status.json', dict(status='COMPLETED', rows=len(rows),
                                     all_restored_at_original=bool(
                                         all(r['counting'] and r['events']
                                             for r in orig))))
    D.dump(out / 'SHA256SUMS.json',
           {str(p.relative_to(out)): D.sha256(p)
            for p in sorted(out.rglob('*')) if p.is_file()})
    print('OUTPUT', out, flush=True)


if __name__ == '__main__':
    main()
