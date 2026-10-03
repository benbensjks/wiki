"""Ground the M-26 design in the actual certified HZH numbers. Read-only.

Prints:
  * the certified event association table (low-order and high-order carries)
  * candidate time windows that contain two low-order carries and exactly one
    high-order carry
  * the measured activation-to-repression offset for BOTH stages, so the
    figure caption can state how far apart they actually are
  * the quantities available per stage, with their ranges in the window
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent            # .../hby_zmh_zmh/failure_attribution
HZZ = ROOT.parent                                 # .../hby_zmh_zmh
PARENT = HZZ.parent                               # .../跨模型三级连接验证
for _d in (PARENT, PARENT / 'hby_zmh_hby',
           PARENT / 'hby_zmh_hby' / 'certification'):
    sys.path.insert(0, str(_d))

import verify_hzh as V                                            # noqa: E402
from hybrid_model import HbyReceiver                              # noqa: E402
from model_hzh import HZHModel, NAMES, IDX                        # noqa: E402
from run_han_comparison import HanInput                           # noqa: E402

CERT = PARENT / 'hby_zmh_hby' / 'certification' / 'results' / '20260928_141723_738780'
TRAJ = PARENT / 'hby_zmh_hby' / 'results' / '20260927_233646_924908' / 'han' / 'trajectory.npz'


def main():
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8')
    print('certified verdict :', (CERT / 'baseline_saved' / 'verdict.json').exists())
    print('baseline traj     :', TRAJ.exists())

    v = json.loads((CERT / 'baseline_saved' / 'verdict.json').read_text(encoding='utf-8'))
    print('certified_v1 =', v['certified_v1'],
          '| steady reads =', v['steady']['reads'],
          '| seq =', v['steady']['sequence'])
    print('global min timing margin =', v['global_min_timing_margin_h'])
    print()

    han = HanInput(300.0)
    m = HZHModel('han', han)
    z = np.load(TRAJ)
    t, y = z['time_h'], z['states']
    assert tuple(z['state_names']) == tuple(NAMES)
    sig = V.signals(y, m.z, m.tail)
    sig['int0'] = y[IDX['b0_I']]

    # sanity: the saved trajectory must reproduce the certified verdict
    vv, _ = V.analyse(t, y, m.z, m.tail)
    print('re-analysed certified_v1 =', vv['certified_v1'],
          '| reads =', vv['steady']['reads'],
          '| seq =', vv['steady']['sequence'])
    print()

    def table(stage):
        rows = v['events'][stage]['associations']
        print(f'--- {stage}: reverse / gate / flip per association ---')
        print('  id  rev_start  rev_peak   gate_start  flip_h     dir     '
              'g_cnt f_cnt')
        for r in rows:
            print('  %2d  %9.3f %9.3f  %10s  %9s  %-6s  %5d %5d' % (
                r['reverse_id'], r['reverse_start_h'], r['reverse_peak_h'],
                'None' if r['gate_start_h'] is None else f"{r['gate_start_h']:.3f}",
                'None' if r['flip_h'] is None else f"{r['flip_h']:.3f}",
                str(r['flip_direction']), r['gate_count'], r['flip_count']))
        return rows

    r0 = table('bit0_to_bit1')
    r1 = table('bit1_to_bit2')
    print()

    # ---- candidate windows: two low-order carries and exactly one high-order
    print('--- candidate windows (need 2 low-order, 1 high-order carry) ---')
    gates0 = V.segments(t, sig['g0'], V.GATE_THRESHOLD, min_dose=V.GATE_MIN_DOSE_H)
    gates1 = V.segments(t, sig['g1'], V.GATE_THRESHOLD, min_dose=V.GATE_MIN_DOSE_H)
    x0 = V.crossings(t, sig['S0'])
    x1 = V.crossings(t, sig['S1'])
    x2 = V.crossings(t, sig['S2'])
    for a in r1:
        if a['gate_count'] != 1 or a['flip_count'] != 1 or a['flip_h'] is None:
            continue
        anchor = a['flip_h']
        for width in (48.0, 50.0, 53.0):
            lo, hi = anchor - width * 0.55, anchor + width * 0.45
            n0 = sum(1 for c in x0 if lo <= c['time_h'] <= hi and c['direction'] == 'down')
            n1 = sum(1 for c in x1 if lo <= c['time_h'] <= hi and c['direction'] == 'down')
            n2 = sum(1 for c in x2 if lo <= c['time_h'] <= hi)
            ng0 = sum(1 for g in gates0 if lo <= g['peak_h'] <= hi)
            ng1 = sum(1 for g in gates1 if lo <= g['peak_h'] <= hi)
            print(f'  anchor flip@{anchor:8.3f} width={width:4.1f} '
                  f'[{lo:7.2f},{hi:7.2f}]  bit0LR->PB={n0} gates0={ng0} '
                  f'bit1LR->PB={n1} gates1={ng1} bit2flips={n2}')
    print()

    # ---- activation vs repression timing, both stages
    print('--- activation vs repression timing (per carry pulse) ---')
    print('  stage  pulse_anchor   A_cross_half  F_cross_half  A->F offset (h)')
    for tag, Acol, Fcol, gates in (
            ('stage0 (ZMH reduced A0/F0)', IDX['A0_zmh'], IDX['F0_zmh'], gates0),
            ('stage1 (HBY A1/F1 chain)  ', IDX['A1'], IDX['F1'], gates1)):
        A = np.asarray(y[Acol], dtype=float)
        F = np.asarray(y[Fcol], dtype=float)
        for g in gates:
            lo, hi = g['start_h'] - 6.0, g['start_h'] + 12.0
            k = (t >= lo) & (t <= hi)
            if k.sum() < 10:
                continue
            tt = t[k]
            a_amp = float(A[k].max())
            f_amp = float(F[k].max())
            a_half = tt[A[k] >= 0.5 * a_amp]
            f_half = tt[F[k] >= 0.5 * f_amp]
            a0 = float(a_half.min()) if a_half.size else float('nan')
            f0 = float(f_half.min()) if f_half.size else float('nan')
            print('  %s  %9.3f   %10.3f   %10.3f   %+8.3f' %
                  (tag, g['start_h'], a0, f0, f0 - a0))
    print()

    # ---- what is available per stage, with ranges in a sample window
    lo, hi = 161.9, 246.7
    k = (t >= lo) & (t <= hi)
    print(f'--- signal ranges over the M-28 window [{lo}, {hi}] ---')
    for name, arr in (('PB0 = 1-b0_S', 1 - y[IDX['b0_S']]),
                      ('A0_zmh', y[IDX['A0_zmh']]),
                      ('F0_zmh', y[IDX['F0_zmh']]),
                      ('g0 (signal)', sig['g0']),
                      ('Int1 = I1_zmh', y[IDX['I1_zmh']]),
                      ('RDF1_zmh', y[IDX['RDF1_zmh']]),
                      ('PB1 = pb1_zmh', y[IDX['pb1_zmh']]),
                      ('A1', y[IDX['A1']]),
                      ('F1', y[IDX['F1']]),
                      ('g1 (signal)', sig['g1']),
                      ('Int2 = b2_I', y[IDX['b2_I']]),
                      ('clock factor', sig['clock']),
                      ('J_rev0', sig['J_rev0']),
                      ('J_rev1', sig['J_rev1']),
                      ('J_rev2', sig['J_rev2'])):
        a = np.asarray(arr, dtype=float)[k]
        print('  %-18s min=%12.5f  max=%12.5f' % (name, a.min(), a.max()))


if __name__ == '__main__':
    main()
