"""Stage 3: is a criterion failure a counting failure, or a threshold artifact?

HZH_MOD8_CAUSAL_V1 has two independent arms:
  (a) steady mod-8 counting  -> verdict['steady']['passed'] / increments_mod8
  (b) both event chains one-to-one, causally ordered, direction alternating
Arm (b) segments the reverse-flux and gate signals with two hard thresholds
(JREV_THRESHOLD_H = 0.1 /h, GATE_THRESHOLD = 0.05). A run can therefore count
perfectly and still fail, if a flux pulse is split in two or a gate pulse never
crosses the threshold.

This script integrates each case ONCE and then re-analyses the identical
trajectory under a sweep of those two thresholds. No re-integration, so the
trajectory is held fixed and only the criterion changes.

Cases:
  uM x0.8      - counts perfectly (19/19 windows, mod-8 increments, commitment
                 1.0) but reports 28 reverse events instead of 14
  clock_K x4   - counts perfectly but reports 0 second-stage gate events
  uM x1.5      - genuine failure: bit1->bit2 never fires (0 gate, 0 reverse)
                 included as the control that a threshold change must NOT rescue
"""
from __future__ import annotations

import json
import sys
from dataclasses import asdict
from datetime import datetime
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

import scan_oat as S                                              # noqa: E402
import verify_hzh as VH                                           # noqa: E402
from hybrid_model import HbyConfig                                # noqa: E402

FROZEN = S.FROZEN

CASES = [
    ('uM_per_au=x0.8', dict(receiver_uM_per_au=FROZEN.receiver_uM_per_au * 0.8)),
    ('clock_K_au=x4', dict(clock_K_au=FROZEN.clock_K_au * 4.0)),
    ('uM_per_au=x1.5', dict(receiver_uM_per_au=FROZEN.receiver_uM_per_au * 1.5)),
]
JREV_SWEEP = [0.02, 0.05, 0.1, 0.2, 0.4]
GATE_SWEEP = [0.02, 0.05, 0.1, 0.2]


def summarise(v):
    return dict(
        certified=bool(v['certified_v1']),
        counts=bool(v['steady']['increments_mod8']),
        reads=int(v['steady']['reads']),
        sequence=str(v['steady']['sequence']),
        min_commitment=(None if v['steady']['minimum_commitment'] is None
                        else float(v['steady']['minimum_commitment'])),
        s0=dict(reverse=int(v['events']['bit0_to_bit1']['reverse_events']),
                gate=int(v['events']['bit0_to_bit1']['gate_events']),
                flip=int(v['events']['bit0_to_bit1']['flip_events']),
                one_to_one=bool(v['events']['bit0_to_bit1']['one_to_one']),
                order=bool(v['events']['bit0_to_bit1']['causal_order']),
                alt=bool(v['events']['bit0_to_bit1']['alternating_directions'])),
        s1=dict(reverse=int(v['events']['bit1_to_bit2']['reverse_events']),
                gate=int(v['events']['bit1_to_bit2']['gate_events']),
                flip=int(v['events']['bit1_to_bit2']['flip_events']),
                one_to_one=bool(v['events']['bit1_to_bit2']['one_to_one']),
                order=bool(v['events']['bit1_to_bit2']['causal_order']),
                alt=bool(v['events']['bit1_to_bit2']['alternating_directions'])),
    )


def main():
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8')
    out = ROOT / 'results' / ('threshold_' + datetime.now().strftime('%Y%m%d_%H%M%S_%f'))
    out.mkdir(parents=True, exist_ok=False)
    S.dump(out / 'status.json', dict(status='RUNNING'))

    han = S.HanInput(S.HOURS)
    assert han.p.peak_load_fraction == 0.0
    prototype = S.HZHModel('han', han)

    # record the defaults we are about to override
    defaults = dict(JREV_THRESHOLD_H=float(VH.JREV_THRESHOLD_H),
                    GATE_THRESHOLD=float(VH.GATE_THRESHOLD),
                    GATE_MIN_DOSE_H=float(VH.GATE_MIN_DOSE_H))
    print('verifier thresholds (defaults):', defaults, flush=True)

    results = {}
    for label, kw in CASES:
        cfg = HbyConfig(**{**asdict(FROZEN), **kw})
        m = S.variant(prototype, cfg)
        t, y = S.integrate(m)
        rec = dict(config={k: float(v) for k, v in kw.items()},
                   clock_range=None, sweeps=[])

        def analyse_at(jrev, gate):
            VH.JREV_THRESHOLD_H = jrev
            VH.GATE_THRESHOLD = gate
            v, sig = VH.analyse(t, y, m.z, m.tail)
            return v

        base_v = analyse_at(defaults['JREV_THRESHOLD_H'], defaults['GATE_THRESHOLD'])
        rec['clock_range'] = [float(x) for x in base_v['signal_ranges']['clock']]
        rec['at_default'] = summarise(base_v)
        print(f'\n=== {label} ===', flush=True)
        print('  default:', json.dumps(rec['at_default']['s0']),
              json.dumps(rec['at_default']['s1']),
              'cert=', rec['at_default']['certified'],
              'counts=', rec['at_default']['counts'], flush=True)

        for jr in JREV_SWEEP:
            for gt in GATE_SWEEP:
                if jr == defaults['JREV_THRESHOLD_H'] and gt == defaults['GATE_THRESHOLD']:
                    continue
                v = analyse_at(jr, gt)
                sm = summarise(v)
                rec['sweeps'].append(dict(jrev=jr, gate=gt, **sm))
        # restore
        VH.JREV_THRESHOLD_H = defaults['JREV_THRESHOLD_H']
        VH.GATE_THRESHOLD = defaults['GATE_THRESHOLD']

        ok = [s for s in rec['sweeps'] if s['certified']]
        cnt = [s for s in rec['sweeps'] if s['counts']]
        print(f'  sweep: {len(rec["sweeps"])} settings, '
              f'{len(cnt)} still count, {len(ok)} become certified', flush=True)
        for s in ok[:6]:
            print(f'    CERTIFIED at jrev={s["jrev"]} gate={s["gate"]} '
                  f'seq={s["sequence"]} s0={s["s0"]["reverse"]}/{s["s0"]["gate"]}/'
                  f'{s["s0"]["flip"]} s1={s["s1"]["reverse"]}/{s["s1"]["gate"]}/'
                  f'{s["s1"]["flip"]}', flush=True)
        results[label] = rec

    S.dump(out / 'threshold_sweep.json', dict(
        defaults=defaults, jrev_sweep=JREV_SWEEP, gate_sweep=GATE_SWEEP,
        cases=results, hours=S.HOURS, note=(
            'The trajectory is integrated once per case and held fixed; only '
            'the two segmentation thresholds of the event arm change.')))
    S.dump(out / 'status.json', dict(status='COMPLETED', cases=len(results)))
    print('\nOUTPUT', out, flush=True)


if __name__ == '__main__':
    main()
