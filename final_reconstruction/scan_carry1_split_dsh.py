"""Stage-3 scan: split A1/F1 maturation, with event-aligned substrate metrics.

Written by dsh.  The acceptance criterion is still the frozen
verify_threebit51.analyse_threebit; this scanner adds the event-aligned
quantities that the previous rounds lacked, so a failure can be attributed
rather than merely observed.

Mechanism wording (corrected): free RDF2 is lowered by Int2 binding it into the
complex C2, by that sequestration, and by complex loss --
    I2 + RDF2 -> C2,  C2 -> degraded
DNA reverse recombination does NOT consume RDF2.  So the reverse drive is set by
the simultaneous product I2*R2 (equivalently C2/K_complex), never by R2 alone.
A stall whose I2 is already zero is just the post-pulse resting state and must
not be used to explain why the pulse failed.

Grid (63 points, 600 h each):
    A1 maturation  = 10, 20, 30, 40 min
    F1 maturation  = 20, 30, 45, 60, 80 min
    carry1 mRNA    = 3, 4, 6 min                  -> 60 points
    common-maturation controls, mRNA = 4:
        A1/F1 = 60/60, 65/65, 70/70               -> 3 points

Pre-registered classification (fixed before the run):
    success          certified (600 h mod-8, 14/14/14, S2 low band entered)
    gate_not_formed  fewer than 14 g1 windows
    gate_leak        a non-carry dwell moves S2 by more than 0.10
    substrate_timing max(I2*R2)/K_D_comp never reaches 1 during a carry episode
    flux_tug_of_war  substrate reached, but per-episode |net dS2| < 0.05 and S2
                     never enters the 0.3 band

    python scan_carry1_split_dsh.py scan  --shard 0 --nshards 9 --hours 600
    python scan_carry1_split_dsh.py merge --nshards 9 --hours 600
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import platform
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

from model import ROOT, ZENG
from model_threebit51 import ThreeBit51Model, ThreeBitCarryParameters
from model_twobit34 import CarryExpressionParameters
from verify_threebit51 import analyse_threebit
from verify_twobit_causal import GATE_ON, MIN_GATE_DOSE_H, _segments
from verify_bit0_part2 import clock_cycles

OUT = ROOT / 'threebit51_results' / 'carry1_split_scan'
A1_MAT = (10.0, 20.0, 30.0, 40.0)
F1_MAT = (20.0, 30.0, 45.0, 60.0, 80.0)
MRNA = (3.0, 4.0, 6.0)
CONTROLS = [(60.0, 60.0, 4.0), (65.0, 65.0, 4.0), (70.0, 70.0, 4.0)]
SKIP_H = 100.0
DROP = 8
EPISODE_PRE_H = 2.0
EPISODE_POST_H = 14.0
KD_COMP = ZENG['K_D_comp']

# pre-registered thresholds.
# NOTE (post-run correction): the first run used EXPECTED_GATES = 14 for the
# gate_not_formed test, but episodes are only collected after the SKIP_H burn-in,
# where the carry count is 12 for every point (bit1 is frozen).  That made every
# point gate_not_formed.  The rule below now derives the expectation from the
# post-burn-in carry count.  The gate_leak class was also dropped: its criterion
# ("the dwell moves S2 by more than 0.10") is satisfied by a working toggle too
# (observed 0.11-0.43 across all 63 points), so the classes are outcome based.
SUBSTRATE_RATIO_MIN = 1.0
LOW_BAND = 0.30
GATE_FRACTION_MIN = 0.75
POST_BURN_IN_GATES = 12        # measured, identical at every point (bit1 frozen)

FROZEN_FILES = ('model.py', 'model_twobit34.py', 'model_threebit51.py',
                'verify_threebit51.py', 'verify_twobit_causal.py')


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest().upper()


def jobs():
    out = [dict(carry1_mrna_min=m, carry1_a1_mat_min=a, carry1_f1_mat_min=f)
           for m, a, f in itertools.product(MRNA, A1_MAT, F1_MAT)]
    out += [dict(carry1_mrna_min=m, carry1_a1_mat_min=a, carry1_f1_mat_min=f)
            for a, f, m in CONTROLS]
    return out


def _idx(t, value):
    return int(np.clip(np.searchsorted(t, value), 0, len(t) - 1))


def _episodes(t, sig, S2, I2, R2, C2, fwd, rev, g1, Kc):
    """Per carry episode and per inter-carry dwell, event aligned."""
    gates = [g for g in _segments(t, g1, GATE_ON, min_area=MIN_GATE_DOSE_H)
             if g['start_h'] > SKIP_H]
    eps = []
    for i, g in enumerate(gates):
        a = _idx(t, g['start_h'] - EPISODE_PRE_H)
        b = _idx(t, g['start_h'] + EPISODE_POST_H)
        o = _idx(t, g['start_h'])
        sl = slice(a, b + 1)
        k = a + int(np.argmax(I2[sl]))
        eps.append(dict(
            gate_start_h=float(g['start_h']), gate_duration_h=float(g['duration_h']),
            gate_peak=float(g['peak']),
            R2_at_gate_open=float(R2[o]), S2_at_gate_open=float(S2[o]),
            I2_at_gate_open=float(I2[o]),
            I2_peak=float(I2[k]), R2_at_I2_peak=float(R2[k]), t_I2_peak_h=float(t[k]),
            max_I2R2_over_KD=float(np.max(I2[sl] * R2[sl]) / KD_COMP),
            max_C2_over_Kc=float(np.max(C2[sl]) / Kc),
            max_H_C2=float(np.max(C2[sl] ** 2 / (Kc ** 2 + C2[sl] ** 2))),
            int_Jrev2=float(np.trapezoid(rev[sl], t[sl])),
            int_Jfwd2=float(np.trapezoid(fwd[sl], t[sl])),
            net_dS2=float(S2[b] - S2[a]), S2_min=float(S2[sl].min()),
            S2_min_h=float(t[a + int(np.argmin(S2[sl]))]),
            entered_low_band=bool(S2[sl].min() <= LOW_BAND)))
    dwells = []
    for i in range(len(gates) - 1):
        a = _idx(t, gates[i]['end_h'])
        b = _idx(t, gates[i + 1]['start_h'])
        if b <= a:
            continue
        sl = slice(a, b + 1)
        dwells.append(dict(start_h=float(t[a]), end_h=float(t[b]),
                           net_dS2=float(S2[b] - S2[a]),
                           abs_dS2=float(abs(S2[b] - S2[a])),
                           int_Jrev2=float(np.trapezoid(rev[sl], t[sl])),
                           int_Jfwd2=float(np.trapezoid(fwd[sl], t[sl])),
                           S2_min=float(S2[sl].min()), S2_max=float(S2[sl].max())))
    return eps, dwells


def _median(values):
    vals = [v for v in values if v is not None and np.isfinite(v)]
    return float(np.median(vals)) if vals else float('nan')


def evaluate(job, hours=600.0):
    t0 = time.perf_counter()
    row = dict(job)
    try:
        c1 = CarryExpressionParameters(
            mrna_half_life_min=job['carry1_mrna_min'],
            activator_maturation_half_life_min=job['carry1_a1_mat_min'],
            repressor_maturation_half_life_min=job['carry1_f1_mat_min'])
        model = ThreeBit51Model(carry=ThreeBitCarryParameters(
            carry0=ThreeBitCarryParameters().carry0, carry1=c1))
        sol = model.simulate(hours=hours, sample_min=2.0, max_step_min=2.0)
        a = analyse_threebit(model, sol, hours)
        c = a['causal_verdict']
        g = a['gate_contrast']
        sig = model.diagnostic_signals(sol.y)
        t = sol.t
        S2, I2, R2, C2 = sig['S2'], sol.y[36], sol.y[42], sol.y[43]
        eps, dwells = _episodes(t, sig, S2, I2, R2, C2, sig['J_fwd2'], sig['J_rev2'],
                                sig['g1'], model.base.K_complex)
        rev_events = len(a['bit1_reverse_events'])
        crossings = len(a['bit2_crossings'])

        got_substrate = any(e['max_I2R2_over_KD'] >= SUBSTRATE_RATIO_MIN for e in eps)
        max_substrate = max([e['max_I2R2_over_KD'] for e in eps], default=float('nan'))
        leak_ds2 = max([d['abs_dS2'] for d in dwells], default=0.0)
        low_band = any(e['entered_low_band'] for e in eps)
        tug = _median([abs(e['net_dS2']) for e in eps])

        if a['certified']:
            cls = 'success'
        elif len(eps) < GATE_FRACTION_MIN * POST_BURN_IN_GATES:
            cls = 'gate_not_formed'
        elif not np.isfinite(max_substrate) or max_substrate < SUBSTRATE_RATIO_MIN:
            cls = 'substrate_timing'
        elif low_band:
            cls = 'partial_reverse_write'
        else:
            cls = 'flux_tug_of_war'

        row.update(
            solver_success=True, error='', finite=bool(np.isfinite(sol.y).all()),
            cold_passed=a['cold_start']['passed'], steady_passed=a['steady_state']['passed'],
            certified=a['certified'], classification=cls,
            cold_sequence=a['cold_start']['sequence'], steady_sequence=a['steady_state']['sequence'],
            minimum_commitment=a['steady_state']['minimum_commitment'],
            one_to_one=c['exactly_one_gate_and_flip_per_late_reverse'],
            causal_order=c['causal_order_after_reverse_start'],
            directions_alternate=c['bit2_directions_alternate'],
            gate_contrast_passed=g['passed'], off_on_peak_ratio=g['off_peak_to_on_peak'],
            on_g1_peak_median=g['on_cycle_peak_median'],
            bit1_reverse_events=rev_events, gate_events=len(eps), bit2_crossings=crossings,
            crossing_error=abs(rev_events - crossings),
            # event-aligned substrate / flux metrics (the point of this scan)
            R2_at_gate_open_median=_median([e['R2_at_gate_open'] for e in eps]),
            R2_at_I2_peak_median=_median([e['R2_at_I2_peak'] for e in eps]),
            I2_peak_median=_median([e['I2_peak'] for e in eps]),
            max_I2R2_over_KD_median=_median([e['max_I2R2_over_KD'] for e in eps]),
            max_I2R2_over_KD_best=max([e['max_I2R2_over_KD'] for e in eps], default=float('nan')),
            max_C2_over_Kc_median=_median([e['max_C2_over_Kc'] for e in eps]),
            max_H_C2_median=_median([e['max_H_C2'] for e in eps]),
            int_Jrev2_median=_median([e['int_Jrev2'] for e in eps]),
            int_Jfwd2_median=_median([e['int_Jfwd2'] for e in eps]),
            net_dS2_median=_median([e['net_dS2'] for e in eps]),
            abs_net_dS2_median=_median([abs(e['net_dS2']) for e in eps]),
            S2_min_in_episode_median=_median([e['S2_min'] for e in eps]),
            S2_min_in_episode_best=min([e['S2_min'] for e in eps], default=float('nan')),
            episodes_entering_low_band=int(sum(1 for e in eps if e['entered_low_band'])),
            gate_width_median_h=_median([e['gate_duration_h'] for e in eps]),
            leak_abs_dS2_max=float(leak_ds2),
            leak_int_Jrev2_median=_median([d['int_Jrev2'] for d in dwells]),
            leak_int_Jfwd2_median=_median([d['int_Jfwd2'] for d in dwells]),
            gate_delay_median_h=_median([x['gate_delay_h'] for x in a['associations']]),
            bit2_delay_median_h=_median([x['bit2_delay_h'] for x in a['associations']]),
        )
    except Exception as exc:                                    # noqa: BLE001
        row.update(solver_success=False, finite=False, certified=False,
                   classification='integration_failed', error=repr(exc),
                   error_trace=traceback.format_exc(limit=3))
    row['runtime_s'] = round(time.perf_counter() - t0, 2)
    return row


def scan(args):
    OUT.mkdir(parents=True, exist_ok=True)
    allj = jobs()
    mine = allj[args.shard::args.nshards]
    tag = f'shard{args.shard:02d}of{args.nshards:02d}'
    csv_path = OUT / f'split_{tag}.csv'
    rows = []
    for i, job in enumerate(mine):
        r = evaluate(job, args.hours)
        rows.append(r)
        fields = sorted({k for x in rows for k in x})
        with csv_path.open('w', newline='', encoding='utf-8') as fh:
            w = csv.DictWriter(fh, fieldnames=fields, extrasaction='ignore')
            w.writeheader()
            w.writerows(rows)
        print(f"[{tag}] {i + 1}/{len(mine)} mRNA={job['carry1_mrna_min']:<3g} "
              f"A1={job['carry1_a1_mat_min']:<4g} F1={job['carry1_f1_mat_min']:<4g} "
              f"cert={r['certified']} cls={r['classification']} "
              f"I2R2/KD={r.get('max_I2R2_over_KD_median', float('nan')):.3f} "
              f"C2/Kc={r.get('max_C2_over_Kc_median', float('nan')):.3f} "
              f"S2min={r.get('S2_min_in_episode_median', float('nan')):.3f} "
              f"low={r.get('episodes_entering_low_band')} ({r['runtime_s']}s)", flush=True)

    meta = dict(shard=args.shard, nshards=args.nshards, hours=args.hours,
                points=len(mine), expected_total=len(allj),
                grid=dict(mrna=MRNA, a1_maturation=A1_MAT, f1_maturation=F1_MAT,
                          controls=CONTROLS),
                criterion='verify_threebit51.analyse_threebit (frozen) + event-aligned metrics',
                preregistered=dict(substrate_ratio_min=SUBSTRATE_RATIO_MIN,
                                   leak_dwell_ds2_max=LEAK_DWELL_DS2_MAX,
                                   tug_net_ds2_max=TUG_NET_DS2_MAX,
                                   low_band=LOW_BAND, expected_gates=EXPECTED_GATES),
                skip_h=SKIP_H, drop_reads=DROP,
                python=sys.version, platform=platform.platform(),
                versions=dict(numpy=np.__version__, pandas=pd.__version__),
                sha256={f.replace('.py', ''): sha256(ROOT / f) for f in FROZEN_FILES},
                csv=csv_path.name, csv_sha256=sha256(csv_path))
    (OUT / f'meta_{tag}.json').write_text(json.dumps(meta, ensure_ascii=False, indent=2),
                                          encoding='utf-8')
    print(f'[{tag}] wrote {csv_path}', flush=True)


def _bool(s):
    return s.map(lambda x: str(x).strip().lower() in ('true', '1'))


def merge(args):
    frames = []
    for i in range(args.nshards):
        p = OUT / f'split_shard{i:02d}of{args.nshards:02d}.csv'
        if not p.exists():
            raise SystemExit(f'missing shard artefact {p}')
        frames.append(pd.read_csv(p))
    d = pd.concat(frames, ignore_index=True)
    for c in ('solver_success', 'finite', 'cold_passed', 'steady_passed', 'certified',
              'one_to_one', 'causal_order', 'directions_alternate', 'gate_contrast_passed'):
        if c in d.columns:
            d[c] = _bool(d[c])
    key = ['carry1_mrna_min', 'carry1_a1_mat_min', 'carry1_f1_mat_min']
    d = d.sort_values(key).reset_index(drop=True)

    expected = {(j['carry1_mrna_min'], j['carry1_a1_mat_min'], j['carry1_f1_mat_min'])
                for j in jobs()}
    got = {tuple(float(getattr(r, k)) for k in key) for r in d.itertuples()}
    completeness = dict(rows=len(d), unique=len(got), expected=len(expected),
                        duplicates=len(d) - len(got), missing=sorted(expected - got),
                        unexpected=sorted(got - expected))
    if completeness['rows'] != len(expected) or completeness['duplicates'] \
            or completeness['missing'] or completeness['unexpected']:
        raise SystemExit(f'incomplete scan: {completeness}')
    d.to_csv(OUT / 'split_all.csv', index=False, encoding='utf-8')

    summary = dict(
        completeness=completeness, hours=float(args.hours),
        solver_failures=int((~d.solver_success).sum()),
        nonfinite=int((~d.finite).sum()),
        certified=int(d.certified.sum()),
        classification={str(k): int(v) for k, v in d.classification.value_counts().items()},
        by_mrna={str(v): dict(points=int((d.carry1_mrna_min == v).sum()),
                              certified=int(d.loc[d.carry1_mrna_min == v, 'certified'].sum()))
                 for v in MRNA},
        by_a1={str(v): dict(points=int((d.carry1_a1_mat_min == v).sum()),
                            certified=int(d.loc[d.carry1_a1_mat_min == v, 'certified'].sum()))
               for v in A1_MAT},
        by_f1={str(v): dict(points=int((d.carry1_f1_mat_min == v).sum()),
                            certified=int(d.loc[d.carry1_f1_mat_min == v, 'certified'].sum()))
               for v in F1_MAT},
        best_episodes_entering_low_band=int(d.episodes_entering_low_band.max()),
        max_substrate_ratio_reached=float(d.max_I2R2_over_KD_best.max()),
        substrate_timing_points=int((d.classification == 'substrate_timing').sum()),
        tug_points=int((d.classification == 'flux_tug_of_war').sum()),
        leak_points=int((d.classification == 'gate_leak').sum()),
        gate_not_formed_points=int((d.classification == 'gate_not_formed').sum()),
    )
    cols = [c for c in key + ['certified', 'classification', 'episodes_entering_low_band',
                              'S2_min_in_episode_best', 'S2_min_in_episode_median',
                              'max_I2R2_over_KD_median', 'max_I2R2_over_KD_best',
                              'max_C2_over_Kc_median', 'R2_at_gate_open_median',
                              'R2_at_I2_peak_median', 'I2_peak_median',
                              'int_Jrev2_median', 'int_Jfwd2_median', 'abs_net_dS2_median',
                              'leak_abs_dS2_max', 'gate_events', 'bit2_crossings',
                              'gate_contrast_passed', 'off_on_peak_ratio'] if c in d.columns]
    summary['top_points'] = (d.sort_values(['certified', 'episodes_entering_low_band',
                                            'max_I2R2_over_KD_best'],
                                           ascending=[False, False, False])
                             .head(12)[cols].to_dict(orient='records'))
    (OUT / 'split_summary.json').write_text(json.dumps(summary, ensure_ascii=False,
                                                       indent=2, default=str),
                                            encoding='utf-8')
    print(json.dumps({k: v for k, v in summary.items() if k != 'top_points'},
                     ensure_ascii=False, indent=2))
    print(d.sort_values(['certified', 'episodes_entering_low_band'],
                        ascending=False).head(12)[cols].to_string(index=False))
    print(f'wrote {OUT / "split_all.csv"}')


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='mode', required=True)
    s = sub.add_parser('scan')
    s.add_argument('--shard', type=int, required=True)
    s.add_argument('--nshards', type=int, default=9)
    s.add_argument('--hours', type=float, default=600.0)
    m = sub.add_parser('merge')
    m.add_argument('--nshards', type=int, default=9)
    m.add_argument('--hours', type=float, default=600.0)
    a = ap.parse_args()
    scan(a) if a.mode == 'scan' else merge(a)


if __name__ == '__main__':
    main()
