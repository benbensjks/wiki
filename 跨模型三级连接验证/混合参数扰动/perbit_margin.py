"""Item 1 follow-up: per-bit band margins and per-bit timing margins.

scan_dense records only the GLOBAL minimum band margin, which is pinned by an
n-independent bit (the same trap as the global timing margin, which is pinned
by bit0's hold). n=3/4/6/8 gave bit-identical global band margins at the same
uM, so that number cannot show whether a steeper gate improves bit2's margin.

This script records, separately for S0, S1 and S2:
  * band margin  - how far the closest sample stays inside its 0.3/0.7 band
  * commitment   - worst in-band occupancy over the read windows
  * setup / hold - timing margins
so the third possible outcome ("the counting boundary does not move, but the
operating margin improves with n") becomes measurable rather than hidden behind
a pass/fail boolean.
"""
from __future__ import annotations

import csv
import json
import sys
import time
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
LOW, HIGH = VH.LOW, VH.HIGH
BITS = ('S0', 'S1', 'S2')


def per_bit(sig, reads, t):
    """Per-bit occupancy and band margin over the read windows.

    Unlabelled windows are NOT skipped. `verify_hzh.read_windows` labels a
    window None exactly when neither band reaches OCC, so skipping them would
    drop precisely the windows that answer "which bit loses commitment first".
    Occupancy is therefore always max(frac<=LOW, frac>=HIGH), and the count of
    unlabelled windows is reported separately.

    band_margin is measured under the majority side (low if frac<=LOW is at
    least frac>=HIGH, else high). For an unlabelled window that margin is
    necessarily negative or below the band, which is the intended signal.
    """
    out = {}
    for b in BITS:
        worst, worst_commit, n, n_unlab = None, None, 0, 0
        for r in reads:
            m = (t >= r['start_h']) & (t <= r['end_h'])
            v = sig[b][m]
            if v.size == 0:
                continue
            lab = r['labels'][BITS.index(b)]
            frac_lo = float(np.mean(v <= LOW))
            frac_hi = float(np.mean(v >= HIGH))
            occ = max(frac_lo, frac_hi)
            side = 0 if frac_lo >= frac_hi else 1
            if lab is None:
                n_unlab += 1
            n += int(v.size)
            d = float(LOW - v.max()) if side == 0 else float(v.min() - HIGH)
            worst = d if worst is None else min(worst, d)
            worst_commit = occ if worst_commit is None else min(worst_commit, occ)
        out[b] = dict(band_margin=worst, commitment=worst_commit,
                      n_samples=n, n_unlabelled=n_unlab, n_windows=len(reads))
    return out


def main():
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8')
    out = ROOT / 'results' / ('perbit_' + datetime.now().strftime('%Y%m%d_%H%M%S_%f'))
    out.mkdir(parents=True, exist_ok=True)
    S.dump(out / 'status.json', dict(status='RUNNING'))

    print('building Han upstream once', flush=True)
    han = S.HanInput(S.HOURS)
    assert han.p.peak_load_fraction == 0.0
    prototype = S.HZHModel('han', han)

    spec = [(n, f) for n in (3.0, 4.0, 6.0, 8.0) for f in (0.80, 1.00, 1.25)]
    rows = []
    for i, (n, f) in enumerate(spec, 1):
        label = f'n={n:g} uM=x{f:.2f}'
        tic = time.perf_counter()
        cfg = HbyConfig(**{**asdict(FROZEN), 'n_A1_gate': n,
                           'receiver_uM_per_au': FROZEN.receiver_uM_per_au * f})
        m = S.variant(prototype, cfg)
        t, y = S.integrate(m)
        v, sig = VH.analyse(t, y, m.z, m.tail)
        sig['_t'] = t
        reads = v['read_windows']
        pb = per_bit(sig, reads, t)
        bm = v['bit_margins']
        row = dict(
            label=label, n_A1_gate=n, uM_fold=f,
            gate_readback=float(m.tail.c.n_A1_gate),
            counting_passed=bool(v['steady']['increments_mod8']),
            event_causality_passed=bool(all(e['passed'] for e in v['events'].values())),
            full_certified=bool(v['certified_v1']),
            sequence=str(v['steady']['sequence']),
            commitment_global=(None if v['steady']['minimum_commitment'] is None
                               else float(v['steady']['minimum_commitment'])),
            band_margin_global=min(
                (pb[b]['band_margin'] for b in BITS
                 if pb[b]['band_margin'] is not None), default=None),
            per_bit={b: dict(pb[b],
                             setup_h=bm[b]['min_setup_h'],
                             hold_h=bm[b]['min_hold_h']) for b in BITS},
            clock_peak_count=int(v['clock_peak_count']),
            runtime_s=time.perf_counter() - tic)
        rows.append(row)
        print(f'[{i}/{len(spec)}] {label:16s} cnt={int(row["counting_passed"])} '
              f'full={int(row["full_certified"])} | '
              + ' '.join(
                  f'{b}: bandM={("None" if pb[b]["band_margin"] is None else format(pb[b]["band_margin"], "+.4f"))}'
                  f' cmt={("None" if pb[b]["commitment"] is None else format(pb[b]["commitment"], ".3f"))}'
                  f' hold={("None" if bm[b]["min_hold_h"] is None else format(bm[b]["min_hold_h"], ".3f"))}'
                  for b in BITS), flush=True)

    S.dump(out / 'perbit.json', dict(
        hours=S.HOURS, sample_min=S.SAMPLE_MIN, rtol=S.RTOL, atol=S.ATOL,
        frozen_config=asdict(FROZEN),
        note=('global band margin is pinned by an n-independent bit; only the '
              'per-bit values can show whether a steeper gate improves bit2'),
        rows=rows))

    with (out / 'perbit.csv').open('w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh)
        cols = ['label', 'n_A1_gate', 'uM_fold', 'gate_readback',
                'counting_passed', 'event_causality_passed', 'full_certified',
                'sequence']
        for b in BITS:
            cols += [f'{b}_band_margin', f'{b}_commitment', f'{b}_n_unlabelled',
                     f'{b}_setup_h', f'{b}_hold_h']
        cols += ['band_margin_global', 'commitment_global', 'clock_peak_count',
                 'runtime_s']
        w.writerow(cols)
        for r in rows:
            vals = [r['label'], r['n_A1_gate'], r['uM_fold'], r['gate_readback'],
                    r['counting_passed'], r['event_causality_passed'],
                    r['full_certified'], r['sequence']]
            for b in BITS:
                p = r['per_bit'][b]
                vals += [p['band_margin'], p['commitment'], p['n_unlabelled'],
                         p['setup_h'], p['hold_h']]
            vals += [r['band_margin_global'], r['commitment_global'],
                     r['clock_peak_count'], round(r['runtime_s'], 3)]
            w.writerow(vals)

    S.dump(out / 'status.json', dict(status='COMPLETED', rows=len(rows)))
    print('OUTPUT', out, flush=True)


if __name__ == '__main__':
    main()
