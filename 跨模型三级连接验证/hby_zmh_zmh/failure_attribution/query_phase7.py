"""Locate the unmatched event for phase 7 at 400 h, BOTH stages.

An earlier version of this query looked only at stage1 (bit1 -> bit2) and
therefore wrongly reported phase 7 as passing. The fresh-process re-evaluation
shows phase 7 has stage0 = [19, 18, 18] and stage1 = [9, 8, 8], so the failure
is in stage0 -- bit0's reverse event without a qualifying bit1 gate.

Uses saved trajectories only; no re-integration.
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT.parent))
sys.path.insert(0, str(ROOT.parent.parent))

import diagnose_round1 as D                                       # noqa: E402
import verify_hzh as V                                            # noqa: E402

HOURS = 400.0
d = sorted((ROOT / 'results').glob('round4_*'))[-1]
han = D.HanInput(HOURS)
prefix = D.HZHModel('han', han)

out = {}
for i in (0, 1, 7):
    z = np.load(d / f'traj_p{i}_{int(HOURS)}h.npz')
    t, y = z['time_h'], z['states']
    m = D.HZZModel(prefix, False)
    m.clock_K = 0.10
    sig = m.signals(y)
    crosses = {b: V.crossings(t, sig[b]) for b in ('S0', 'S1', 'S2')}
    rec = {}
    for j, bit in ((0, 'S1'), (1, 'S2')):
        rev = V.segments(t, sig[f'J_rev{j}'], V.JREV_THRESHOLD_H)
        gate = V.segments(t, sig[f'g{j}'], V.GATE_THRESHOLD,
                          min_dose=V.GATE_MIN_DOSE_H)
        a = V.associate(rev, gate, crosses[bit], HOURS)
        rows = a['associations']
        bad = [r for r in rows
               if not (r['gate_count'] == 1 and r['flip_count'] == 1)]
        rec[f'stage{j}'] = dict(
            bit=bit, reverse=len(rev), gate=len(gate), flips=len(crosses[bit]),
            passed=a['passed'], one_to_one=a['one_to_one'],
            order=a['causal_order'], alt=a['alternating_directions'],
            late_rows=a['late_associations'],
            dropped_first_2=[r['reverse_id'] for r in rows[:2]],
            unmatched=[dict(reverse_id=r['reverse_id'],
                            reverse_peak_h=r['reverse_peak_h'],
                            gap_to_end_h=float(HOURS - r['reverse_peak_h']),
                            in_dropped_burnin=bool(r['reverse_id'] < 2),
                            gate_count=r['gate_count'],
                            flip_count=r['flip_count']) for r in bad],
            unassigned_gate=a['unassigned_gate_events'],
            unassigned_flip=a['unassigned_flip_events'],
            last_reverse_h=float(rev[-1]['peak_h']) if rev else None,
            last_gate_start_h=float(gate[-1]['start_h']) if gate else None)
    out[f'phase{i}'] = rec

for k, rec in out.items():
    print('=' * 72)
    print(k)
    for st, v in rec.items():
        print(f"  {st} ({v['bit']}): rev={v['reverse']} gate={v['gate']} "
              f"flips={v['flips']} passed={v['passed']} "
              f"o2o={v['one_to_one']} ord={v['order']} alt={v['alt']} "
              f"late={v['late_rows']} dropped={v['dropped_first_2']} "
              f"unassigned_g={v['unassigned_gate']} unassigned_f={v['unassigned_flip']}")
        for r in v['unmatched']:
            print(f"     unmatched reverse #{r['reverse_id']} at "
                  f"{r['reverse_peak_h']:.2f} h "
                  f"({r['gap_to_end_h']:.2f} h before the end) "
                  f"in_dropped_burnin={r['in_dropped_burnin']} "
                  f"gate={r['gate_count']} flip={r['flip_count']}")
        print(f"     last reverse {v['last_reverse_h']} h, "
              f"last gate start {v['last_gate_start_h']} h")

(ROOT / 'results' / f'phase_unmatched_both_stages_{d.name}.json').write_text(
    json.dumps(out, ensure_ascii=False, indent=2), encoding='utf-8')
