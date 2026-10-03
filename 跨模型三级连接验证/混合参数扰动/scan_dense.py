"""Stage 4: densify the two boundary families and record continuous diagnostics.

Answers verification items 1 and 2:
  (1) for each n_A1_gate in {3,4,6,8}, locate the lower (uM fold 0.50-0.80) and
      upper (1.25-1.50) counting boundary of receiver_uM_per_au at 0.05 steps;
  (2) locate the K_A[1] and K_F[1] boundaries over 1.25-1.50 at 0.05 steps,
      at n_A1_gate = 4 and 6.

A boolean pass/fail hides the third possible outcome, so every run also records
continuous diagnostics:
  * band_margin_h  - how far the closest sample stays inside its read band
  * commitment     - worst in-band occupancy over the read windows
  * g1 peak / width-at-half-max / per-cycle integral / far-dwell residual
  * bit2 forward and reverse flux integrals
  * per-bit setup and hold margins
"""
from __future__ import annotations

import argparse
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
from hybrid_model import HbyConfig                                # noqa: E402
import verify_hzh as VH                                           # noqa: E402

FROZEN = S.FROZEN
ZENG = S.ZENG
LOW, HIGH = VH.LOW, VH.HIGH


def band_margin(sig, reads):
    """Smallest distance of any sample from the band edge it must lie inside.

    For a bit read as LOW the samples must stay below LOW; for HIGH, above
    HIGH. Positive = inside the band; negative = a sample crossed the edge.
    Returns None when no bit was decisively labelled.
    """
    worst = None
    for r in reads:
        mask = (sig['_t'] >= r['start_h']) & (sig['_t'] <= r['end_h'])
        for name, label in zip(('S0', 'S1', 'S2'), r['labels']):
            if label is None:
                continue
            v = sig[name][mask]
            if v.size == 0:
                continue
            d = float(LOW - v.max()) if label == 0 else float(v.min() - HIGH)
            worst = d if worst is None else min(worst, d)
    return worst


def cycle_stats(sig, reads, name):
    """Per-cycle peak, width at half maximum and integral of a signal."""
    t, y = sig['_t'], sig[name]
    peaks, widths, integrals, residuals = [], [], [], []
    for r in reads:
        m = (t >= r['cycle_start_h']) & (t <= r['cycle_end_h'])
        if m.sum() < 3:
            continue
        tt, yy = t[m], y[m]
        pk = float(yy.max())
        peaks.append(pk)
        if pk > 0:
            half = pk / 2.0
            above = tt[yy >= half]
            widths.append(float(above.max() - above.min()) if above.size > 1 else 0.0)
        integrals.append(float(np.trapezoid(np.maximum(yy, 0.0), tt)))
        # far-dwell residual: the 20 % of the cycle furthest from the peak
        k = max(1, int(0.2 * tt.size))
        idx = np.argsort(yy)[:k]
        residuals.append(float(np.median(yy[idx])))
    return dict(peak_max=(max(peaks) if peaks else None),
                peak_min=(min(peaks) if peaks else None),
                width_max=(max(widths) if widths else None),
                integral_min=(min(integrals) if integrals else None),
                integral_max=(max(integrals) if integrals else None),
                far_dwell_median=(float(np.median(residuals)) if residuals else None))


def evaluate(prototype, label, cfg, zeng_patch=None, extra=None):
    row = dict(label=label)
    if extra:
        row.update(extra)
    tic = time.perf_counter()
    try:
        m = S.variant(prototype, cfg, zeng_patch)
        row['work_point'] = S.work_point_readback(m)
        t, y = S.integrate(m)
        v, sig = VH.analyse(t, y, m.z, m.tail)
        sig['_t'] = t
        reads = v['read_windows']
        bm = v['bit_margins']
        s0, s1 = (v['events']['bit0_to_bit1'], v['events']['bit1_to_bit2'])

        def ev(e):
            return dict(reverse=int(e['reverse_events']), gate=int(e['gate_events']),
                        flip=int(e['flip_events']),
                        one_to_one=bool(e['one_to_one']),
                        order=bool(e['causal_order']),
                        alt=bool(e['alternating_directions']))
        events_ok = all(x['passed'] for x in (s0, s1))
        counting = bool(v['steady']['increments_mod8'])
        row.update(
            counting_passed=counting,
            event_causality_passed=bool(events_ok),
            full_certified=bool(v['certified_v1']),
            steady_reads=int(v['steady']['reads']),
            sequence=str(v['steady']['sequence']),
            commitment=(None if v['steady']['minimum_commitment'] is None
                        else float(v['steady']['minimum_commitment'])),
            boundary_clips=int(v['steady']['boundary_clips']),
            band_margin_h=band_margin(sig, reads),
            margin_global_h=(None if v['global_min_timing_margin_h'] is None
                             else float(v['global_min_timing_margin_h'])),
            margin_S0_h=S.min_of(bm.get('S0')),
            margin_S1_h=S.min_of(bm.get('S1')),
            margin_S2_h=S.min_of(bm.get('S2')),
            clock_peak_count=int(v['clock_peak_count']),
            g1=cycle_stats(sig, reads, 'g1'),
            clock=cycle_stats(sig, reads, 'clock'),
            j_fwd2_integral=float(np.trapezoid(np.maximum(sig['J_fwd2'], 0.0), t)),
            j_rev2_integral=float(np.trapezoid(np.maximum(sig['J_rev2'], 0.0), t)),
            stage0=ev(s0), stage1=ev(s1),
            signal_ranges={k: [float(x) for x in v['signal_ranges'][k]]
                           for k in ('g0', 'g1', 'clock', 'int0')},
            error=None)
    except Exception as exc:
        row.update(counting_passed=False, event_causality_passed=False,
                   full_certified=False, steady_reads=0, sequence='',
                   commitment=None, boundary_clips=0, band_margin_h=None,
                   margin_global_h=None, margin_S0_h=None, margin_S1_h=None,
                   margin_S2_h=None, clock_peak_count=0, g1=None, clock=None,
                   j_fwd2_integral=None, j_rev2_integral=None, stage0=None,
                   stage1=None, signal_ranges=None, error=repr(exc))
    row['runtime_s'] = time.perf_counter() - tic
    return row


def spec():
    out = []
    for n in (3.0, 4.0, 6.0, 8.0):
        for f in (0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80):
            out.append((f'n={n:g} uM=x{f:.2f}', 'DENSE_LOW', 'config',
                        'receiver_uM_per_au', FROZEN.receiver_uM_per_au * f, n))
        for f in (1.25, 1.30, 1.35, 1.40, 1.45, 1.50):
            out.append((f'n={n:g} uM=x{f:.2f}', 'DENSE_HIGH', 'config',
                        'receiver_uM_per_au', FROZEN.receiver_uM_per_au * f, n))
    for n in (4.0, 6.0):
        for f in (1.25, 1.30, 1.35, 1.40, 1.45, 1.50):
            out.append((f'n={n:g} K_A[1]=x{f:.2f}', 'DENSE_KA', 'zeng',
                        ('K_A', 1), ZENG['K_A'][1] * f, n))
            out.append((f'n={n:g} K_F[1]=x{f:.2f}', 'DENSE_KF', 'zeng',
                        ('K_F', 1), ZENG['K_F'][1] * f, n))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--families', default='ALL')
    ap.add_argument('--runs', type=int, default=0)
    args = ap.parse_args()
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8')

    out = ROOT / 'results' / ('dense_' + datetime.now().strftime('%Y%m%d_%H%M%S_%f'))
    out.mkdir(parents=True, exist_ok=False)
    S.dump(out / 'status.json', dict(status='RUNNING'))

    print('building Han upstream once', flush=True)
    han = S.HanInput(S.HOURS)
    assert han.p.peak_load_fraction == 0.0
    prototype = S.HZHModel('han', han)

    rows = []
    base = evaluate(prototype, 'BASELINE', FROZEN)
    checks = dict(counting=base['counting_passed'] is True,
                  events=base['event_causality_passed'] is True,
                  full=base['full_certified'] is True,
                  reads=base['steady_reads'] == S.REF['steady_reads'],
                  sequence=base['sequence'] == S.REF['sequence'],
                  margin=abs(base['margin_global_h'] - S.REF['margin_h']) < 1e-6,
                  gate_readback=base['work_point']['n_A1_gate_effective'] == 6.0)
    print('baseline guard:', json.dumps(checks), flush=True)
    if not all(checks.values()):
        S.dump(out / 'status.json', dict(status='FAILED', checks=checks))
        raise RuntimeError(f'Baseline guard failed: {checks}')
    rows.append(base)
    print(f"  baseline: band_margin {base['band_margin_h']:.4f} | "
          f"commitment {base['commitment']:.4f} | "
          f"g1 peak_max {base['g1']['peak_max']:.4f} | "
          f"g1 integral_min {base['g1']['integral_min']:.3e}", flush=True)

    sp = spec()
    if args.families != 'ALL':
        sp = [r for r in sp if r[1] in args.families.split(',')]
    if args.runs:
        sp = sp[:args.runs]
    print(f'rows: {len(sp)}', flush=True)
    for i, (label, fam, kind, target, value, n) in enumerate(sp, 1):
        cfg = HbyConfig(**{**asdict(FROZEN), 'n_A1_gate': n,
                           **({'receiver_uM_per_au': float(value)}
                              if kind == 'config' else {})})
        patch = None if kind == 'config' else {target: float(value)}
        r = evaluate(prototype, label, cfg, patch,
                     extra=dict(family=fam, kind=kind, target=str(target),
                                value=float(value), n_A1_gate=n))
        rows.append(r)
        print(f'[{i}/{len(sp)}] {label:22s} {fam:10s} '
              f'cnt={int(r["counting_passed"])} ev={int(r["event_causality_passed"])} '
              f'full={int(r["full_certified"])} '
              f'bandM={("None" if r["band_margin_h"] is None else format(r["band_margin_h"], "+.4f")):>7s} '
              f'S2m={("None" if r["margin_S2_h"] is None else format(r["margin_S2_h"], ".3f")):>5s} '
              f'{r["runtime_s"]:.0f}s', flush=True)

    S.dump(out / 'dense_all.json', dict(
        hours=S.HOURS, sample_min=S.SAMPLE_MIN, max_step_min=S.MAX_STEP_MIN,
        rtol=S.RTOL, atol=S.ATOL, frozen_config=asdict(FROZEN),
        reference=S.REF, baseline_guard=checks, rows=rows))

    cols = ['label', 'family', 'n_A1_gate', 'kind', 'target', 'value',
            'counting_passed', 'event_causality_passed', 'full_certified',
            'steady_reads', 'sequence', 'commitment', 'boundary_clips',
            'band_margin_h', 'margin_global_h', 'margin_S0_h', 'margin_S1_h',
            'margin_S2_h', 'clock_peak_count',
            'g1_peak_max', 'g1_peak_min', 'g1_width_max', 'g1_integral_min',
            'g1_integral_max', 'g1_far_dwell_median',
            'clock_peak_max', 'clock_width_max', 'clock_integral_min',
            'j_fwd2_integral', 'j_rev2_integral',
            's0_rev', 's0_gate', 's0_flip', 's0_o2o', 's0_ord', 's0_alt',
            's1_rev', 's1_gate', 's1_flip', 's1_o2o', 's1_ord', 's1_alt',
            'g1_range_lo', 'g1_range_hi', 'clock_range_lo', 'clock_range_hi',
            'runtime_s', 'error']
    with (out / 'dense_all.csv').open('w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for r in rows:
            g1 = r.get('g1') or {}
            ck = r.get('clock') or {}
            s0 = r.get('stage0') or {}
            s1 = r.get('stage1') or {}
            sr = r.get('signal_ranges') or {}
            gr = sr.get('g1') or [None, None]
            cr = sr.get('clock') or [None, None]
            w.writerow([
                r['label'], r.get('family', 'BASE'), r.get('n_A1_gate'),
                r.get('kind', ''), r.get('target', ''), r.get('value', ''),
                r['counting_passed'], r['event_causality_passed'],
                r['full_certified'], r['steady_reads'], r['sequence'],
                r['commitment'], r['boundary_clips'], r['band_margin_h'],
                r['margin_global_h'], r['margin_S0_h'], r['margin_S1_h'],
                r['margin_S2_h'], r['clock_peak_count'],
                g1.get('peak_max'), g1.get('peak_min'), g1.get('width_max'),
                g1.get('integral_min'), g1.get('integral_max'),
                g1.get('far_dwell_median'),
                ck.get('peak_max'), ck.get('width_max'), ck.get('integral_min'),
                r['j_fwd2_integral'], r['j_rev2_integral'],
                s0.get('reverse'), s0.get('gate'), s0.get('flip'),
                s0.get('one_to_one'), s0.get('order'), s0.get('alt'),
                s1.get('reverse'), s1.get('gate'), s1.get('flip'),
                s1.get('one_to_one'), s1.get('order'), s1.get('alt'),
                gr[0], gr[1], cr[0], cr[1],
                round(r['runtime_s'], 3), r['error']])

    S.dump(out / 'status.json', dict(status='COMPLETED', rows=len(rows)))
    print('OUTPUT', out, flush=True)


if __name__ == '__main__':
    main()
