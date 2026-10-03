"""Re-evaluate all eight saved 400 h phase trajectories in a FRESH process.

Round 4's in-process Part A recorded phase 7 as event_causality_passed = False,
but re-evaluating the same saved trajectory in a fresh process gives True. This
script settles which is right and, if they differ, finds the cause.
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
from diagnose_round3 import evaluate                              # noqa: E402

HOURS = 400.0
d = sorted((ROOT / 'results').glob('round4_*'))[-1]
print('dir', d.name)
print('module thresholds at import:',
      'JREV =', V.JREV_THRESHOLD_H, ' GATE =', V.GATE_THRESHOLD,
      ' PULSE =', V.PULSE_DURATION_H, ' DOSEMIN =', V.GATE_MIN_DOSE_H)

han = D.HanInput(HOURS)
prefix = D.HZHModel('han', han)

rows = []
for i in range(8):
    z = np.load(d / f'traj_p{i}_{int(HOURS)}h.npz')
    t, y = z['time_h'], z['states']
    m = D.HZZModel(prefix, False)
    m.clock_K = 0.10
    sig = m.signals(y)
    r = evaluate(t, sig, HOURS)
    s1 = r['events']['bit1_to_bit2']
    s0 = r['events']['bit0_to_bit1']
    rows.append(dict(
        phase=i, counting=r['counting_passed'],
        events=r['event_causality_passed'],
        stage0=[s0['reverse_events'], s0['gate_events'], s0['flip_events']],
        stage1=[s1['reverse_events'], s1['gate_events'], s1['flip_events']],
        late1=s1['late_associations'],
        unassigned_gate1=s1['unassigned_gate_events'],
        unassigned_flip1=s1['unassigned_flip_events'],
        one_to_one=s1['one_to_one'], order=s1['causal_order'],
        alt=s1['alternating_directions'],
        sequence=r['steady']['sequence']))
    print(f"  phase {i}: cnt={int(r['counting_passed'])} "
          f"ev={int(r['event_causality_passed'])} "
          f"stage0={rows[-1]['stage0']} stage1={rows[-1]['stage1']} "
          f"o2o={s1['one_to_one']} ord={s1['causal_order']} "
          f"alt={s1['alternating_directions']} "
          f"unassigned_g={s1['unassigned_gate_events']} "
          f"unassigned_f={s1['unassigned_flip_events']}")

print()
print('counting %d/8, events %d/8' % (sum(r['counting'] for r in rows),
                                      sum(r['events'] for r in rows)))
(ROOT / 'results' / f'phase_reeval_{d.name}.json').write_text(
    json.dumps(dict(dir=d.name, rows=rows), ensure_ascii=False, indent=2),
    encoding='utf-8')
