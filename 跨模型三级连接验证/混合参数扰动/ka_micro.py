"""K_A[1] micro-scan: bracket the n=4 boundary so the shift can be quantified.

dense_all.csv gives n=4 pass at fold x1.00 and fail at x1.25, and n=6 pass at
x1.30 / fail at x1.35 -- so the only defensible statement is "the upper bound
moves up", with no percentage. This fills x1.05-1.20 at both n to measure it.

Also records per-bit commitment and band margin, because dense showed that at
the n=6 uM x1.35 point the counting failure is a state-commitment failure
(commitment 0.0) with both event chains still complete.
"""
from __future__ import annotations

import csv, json, sys, time
from dataclasses import asdict
from datetime import datetime
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

import scan_oat as S                                              # noqa: E402
import verify_hzh as VH                                           # noqa: E402
import perbit_margin as PB                                        # noqa: E402
from hybrid_model import HbyConfig                                # noqa: E402

FROZEN, ZENG = S.FROZEN, S.ZENG


def main():
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8')
    out = ROOT / 'results' / ('kamicro_' + datetime.now().strftime('%Y%m%d_%H%M%S_%f'))
    out.mkdir(parents=True, exist_ok=True)
    S.dump(out / 'status.json', dict(status='RUNNING'))

    print('building Han upstream once', flush=True)
    han = S.HanInput(S.HOURS)
    assert han.p.peak_load_fraction == 0.0
    prototype = S.HZHModel('han', han)

    spec = [(n, f) for n in (4.0, 6.0)
            for f in (1.00, 1.05, 1.10, 1.15, 1.20, 1.25)]
    rows = []
    for i, (n, f) in enumerate(spec, 1):
        label = f'n={n:g} K_A[1]=x{f:.2f}'
        tic = time.perf_counter()
        cfg = HbyConfig(**{**asdict(FROZEN), 'n_A1_gate': n})
        value = ZENG['K_A'][1] * f
        m = S.variant(prototype, cfg, {('K_A', 1): value})
        t, y = S.integrate(m)
        v, sig = VH.analyse(t, y, m.z, m.tail)
        sig['_t'] = t
        reads = v['read_windows']
        pb = PB.per_bit(sig, reads, t)
        bm = v['bit_margins']
        s0, s1 = v['events']['bit0_to_bit1'], v['events']['bit1_to_bit2']
        row = dict(
            label=label, n_A1_gate=n, fold=f, K_A1_value=value,
            gate_readback=float(m.tail.c.n_A1_gate),
            K_A1_readback=float(m.p['K_A'][1]),
            counting_passed=bool(v['steady']['increments_mod8']),
            event_causality_passed=bool(all(e['passed'] for e in v['events'].values())),
            full_certified=bool(v['certified_v1']),
            sequence=str(v['steady']['sequence']),
            commitment_global=(None if v['steady']['minimum_commitment'] is None
                               else float(v['steady']['minimum_commitment'])),
            per_bit={b: dict(pb[b], setup_h=bm[b]['min_setup_h'],
                             hold_h=bm[b]['min_hold_h']) for b in PB.BITS},
            s0=dict(reverse=int(s0['reverse_events']), gate=int(s0['gate_events']),
                    flip=int(s0['flip_events']), one_to_one=bool(s0['one_to_one'])),
            s1=dict(reverse=int(s1['reverse_events']), gate=int(s1['gate_events']),
                    flip=int(s1['flip_events']), one_to_one=bool(s1['one_to_one'])),
            g1_peak=float(np.max(sig['g1'])),
            clock_peak=float(np.max(sig['clock'])),
            runtime_s=time.perf_counter() - tic)
        rows.append(row)
        print(f'[{i}/{len(spec)}] {label:20s} val={value:.3f} '
              f'cnt={int(row["counting_passed"])} ev={int(row["event_causality_passed"])} '
              f'full={int(row["full_certified"])} cmt={row["commitment_global"]} '
              f'| S2 cmt={pb["S2"]["commitment"]} bandM={pb["S2"]["band_margin"]} '
              f'g1pk={row["g1_peak"]:.4f}', flush=True)

    S.dump(out / 'kamicro.json', dict(
        hours=S.HOURS, sample_min=S.SAMPLE_MIN, rtol=S.RTOL, atol=S.ATOL,
        frozen_config=asdict(FROZEN), rows=rows))
    with (out / 'kamicro.csv').open('w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh)
        cols = ['label', 'n_A1_gate', 'fold', 'K_A1_value', 'K_A1_readback',
                'counting_passed', 'event_causality_passed', 'full_certified',
                'sequence', 'commitment_global']
        for b in PB.BITS:
            cols += [f'{b}_commitment', f'{b}_band_margin', f'{b}_hold_h']
        cols += ['s0_rev', 's0_gate', 's0_flip', 's1_rev', 's1_gate', 's1_flip',
                 'g1_peak', 'clock_peak', 'runtime_s']
        w.writerow(cols)
        for r in rows:
            vals = [r['label'], r['n_A1_gate'], r['fold'], r['K_A1_value'],
                    r['K_A1_readback'], r['counting_passed'],
                    r['event_causality_passed'], r['full_certified'],
                    r['sequence'], r['commitment_global']]
            for b in PB.BITS:
                p = r['per_bit'][b]
                vals += [p['commitment'], p['band_margin'], p['hold_h']]
            vals += [r['s0']['reverse'], r['s0']['gate'], r['s0']['flip'],
                     r['s1']['reverse'], r['s1']['gate'], r['s1']['flip'],
                     r['g1_peak'], r['clock_peak'], round(r['runtime_s'], 3)]
            w.writerow(vals)
    S.dump(out / 'status.json', dict(status='COMPLETED', rows=len(rows)))
    print('OUTPUT', out, flush=True)


if __name__ == '__main__':
    main()
