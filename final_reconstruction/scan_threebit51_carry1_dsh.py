"""Stage-1 scan of the new A1/F1 expression parameters in the 51-state model.

Written by dsh.  It does NOT overwrite or import scan_threebit51_carry1.py; the
acceptance criterion is the frozen verify_threebit51.analyse_threebit, used
verbatim, so the two scanners are directly comparable.

Grid (35 points, 600 h each):
    carry1 mRNA half-life      0.5, 1, 2, 4, 8 min
    A1/F1 common maturation    5, 10, 20, 30, 40, 50, 60 min

Only the new A1/F1 expression parameters move.  The frozen 34-state prefix, the
frozen extension (uM 5.75, bit maturation 20 min, carry0 32.5 min, clock gate
0.4/3, add_growth False) and ZENG are all untouched.

Beyond the first-round fields this scanner also records, per late read cycle,
the bit2 forward/reverse flux integrals and the net S2 change, because the
first-round failure was "per-cycle net recombination exactly zero", which the
boolean verdict alone cannot rank.

    python scan_threebit51_carry1_dsh.py scan  --shard 0 --nshards 7 --hours 600
    python scan_threebit51_carry1_dsh.py merge --nshards 7
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

from model import ROOT
from model_threebit51 import (ThreeBit51Model, ThreeBitCarryParameters,
                              selected_threebit_extension)
from model_twobit34 import CarryExpressionParameters
from verify_threebit51 import analyse_threebit

OUT = ROOT / 'threebit51_results' / 'carry1_scan_dsh'
MRNA = (0.5, 1.0, 2.0, 4.0, 8.0)
MAT = (5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0)
DROP = 8                       # one whole mod-8 period, per the acceptance rule

FROZEN_FILES = ('model.py', 'model_twobit34.py', 'model_threebit51.py',
                'verify_threebit51.py', 'verify_twobit_causal.py',
                'test_threebit51.py', 'scan_threebit51_carry1.py')


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest().upper()


def jobs():
    return [dict(carry1_mrna_min=a, carry1_maturation_min=b)
            for a, b in itertools.product(MRNA, MAT)]


def _cycle_flux_metrics(t, sig, reads, drop=DROP):
    """Per late cycle: bit2 flux integrals, net S2 change, band occupancy."""
    rows = []
    for r in reads[drop:]:
        a = int(np.searchsorted(t, r['cycle_start_h']))
        b = int(np.searchsorted(t, r['cycle_end_h']))
        if b <= a:
            continue
        sl = slice(a, b + 1)
        S2 = sig['S2'][sl]
        fwd = float(np.trapezoid(sig['J_fwd2'][sl], t[sl]))
        rev = float(np.trapezoid(sig['J_rev2'][sl], t[sl]))
        rows.append(dict(cycle=r['cycle'], s2_min=float(S2.min()), s2_max=float(S2.max()),
                         s2_start=float(S2[0]), s2_end=float(S2[-1]),
                         net_ds2=float(S2[-1] - S2[0]),
                         fwd2_integral=fwd, rev2_integral=rev,
                         rev_over_fwd=(rev / fwd if fwd > 0 else np.nan),
                         bit2_label=r['bits'][2]['label'],
                         bit2_commitment=float(r['bits'][2]['commitment'])))
    return rows


def _summarise_cycles(cycles):
    if not cycles:
        return dict(cycles=0)
    s2min = np.array([c['s2_min'] for c in cycles])
    s2max = np.array([c['s2_max'] for c in cycles])
    net = np.array([c['net_ds2'] for c in cycles])
    rr = np.array([c['rev_over_fwd'] for c in cycles], dtype=float)
    fwd = np.array([c['fwd2_integral'] for c in cycles])
    rev = np.array([c['rev2_integral'] for c in cycles])
    labels = [c['bit2_label'] for c in cycles]
    return dict(cycles=len(cycles),
                s2_cycle_min_median=float(np.median(s2min)),
                s2_cycle_min_min=float(s2min.min()),
                s2_cycle_max_median=float(np.median(s2max)),
                abs_net_ds2_median=float(np.median(np.abs(net))),
                abs_net_ds2_max=float(np.abs(net).max()),
                fwd2_integral_median=float(np.median(fwd)),
                rev2_integral_median=float(np.median(rev)),
                rev_over_fwd_median=float(np.nanmedian(rr)),
                bit2_labelled_cycles=int(sum(1 for l in labels if l != 'x')),
                bit2_high_cycles=int(sum(1 for l in labels if l == '1')),
                bit2_low_cycles=int(sum(1 for l in labels if l == '0')))


def evaluate(job, hours=600.0):
    t0 = time.perf_counter()
    row = dict(job)
    try:
        c1 = CarryExpressionParameters(
            mrna_half_life_min=job['carry1_mrna_min'],
            activator_maturation_half_life_min=job['carry1_maturation_min'],
            repressor_maturation_half_life_min=job['carry1_maturation_min'])
        carry = ThreeBitCarryParameters(carry0=ThreeBitCarryParameters().carry0, carry1=c1)
        model = ThreeBit51Model(carry=carry)
        sol = model.simulate(hours=hours, sample_min=2.0, max_step_min=2.0)
        a = analyse_threebit(model, sol, hours)
        c = a['causal_verdict']
        g = a['gate_contrast']
        sig = model.diagnostic_signals(sol.y)
        cycles = _cycle_flux_metrics(sol.t, sig, a['read_windows'])
        summ = _summarise_cycles(cycles)
        rev_events = len(a['bit1_reverse_events'])
        crossings = len(a['bit2_crossings'])

        def _med(values):
            vals = [v for v in values if v is not None and np.isfinite(v)]
            return float(np.median(vals)) if vals else float('nan')

        gates_ev = a['carry1_gate_events']
        rev_ev = a['bit1_reverse_events']
        assoc = a['associations']
        # The first-round stall is a balance point between the reverse push and the
        # collapse of the RDF2 pool (gamma_rdf = 0.8/h, tau = 1.25 h), so the gate
        # window width relative to that time constant is recorded explicitly.
        row.update(
            gate_duration_median_h=_med([g['duration_h'] for g in gates_ev]),
            gate_dose_median=_med([g['dose'] for g in gates_ev]),
            gate_peak_median=_med([g['peak'] for g in gates_ev]),
            rev_event_duration_median_h=_med([e['duration_h'] for e in rev_ev]),
            rev_event_peak_median=_med([e['peak'] for e in rev_ev]),
            gate_delay_median_h=_med([x['gate_delay_h'] for x in assoc]),
            bit2_delay_median_h=_med([x['bit2_delay_h'] for x in assoc]),
            rdf2_tau_over_gate_width=(
                (1.0 / 0.8) / _med([g['duration_h'] for g in gates_ev])
                if gates_ev and _med([g['duration_h'] for g in gates_ev]) > 0 else float('nan')),
        )
        row.update(
            solver_success=True, error='', finite=bool(np.isfinite(sol.y).all()),
            cold_passed=a['cold_start']['passed'], steady_passed=a['steady_state']['passed'],
            cold_sequence=a['cold_start']['sequence'], steady_sequence=a['steady_state']['sequence'],
            minimum_commitment=a['steady_state']['minimum_commitment'],
            certified=a['certified'],
            one_to_one=c['exactly_one_gate_and_flip_per_late_reverse'],
            causal_order=c['causal_order_after_reverse_start'],
            directions_alternate=c['bit2_directions_alternate'],
            unassigned_gates=len(c['unassigned_gate_events']),
            unassigned_crossings=len(c['unassigned_bit2_crossings']),
            gate_contrast_passed=g['passed'],
            off_on_peak_ratio=g['off_peak_to_on_peak'],
            off_median_to_on_peak=g['off_median_to_on_peak'],
            on_g1_peak_median=g['on_cycle_peak_median'],
            g1_peak=sig['g1'].max(), Int2_source_peak=float(sig['Int2_source'].max()),
            bit1_reverse_events=rev_events, gate_events=len(a['carry1_gate_events']),
            bit2_crossings=crossings,
            crossing_error=abs(rev_events - crossings),
            crossing_ratio=(crossings / rev_events if rev_events else np.nan),
            bit2_setup_h=a['bit_margins']['bit2']['min_setup_h'],
            bit2_hold_h=a['bit_margins']['bit2']['min_hold_h'],
            S2_min=float(sig['S2'].min()), S2_max=float(sig['S2'].max()),
            RDF2_median=float(np.median(sol.y[42])), RDF2_max=float(sol.y[42].max()),
            RDF2_min=float(sol.y[42].min()),
            I2_peak=float(sol.y[36].max()),
            windows=len(a['read_windows']), **summ)
        row['cycle_metrics'] = cycles
    except Exception as exc:                                    # noqa: BLE001
        row.update(solver_success=False, finite=False, certified=False,
                   error=repr(exc), error_trace=traceback.format_exc(limit=3))
    row['runtime_s'] = round(time.perf_counter() - t0, 2)
    return row


def scan(args):
    OUT.mkdir(parents=True, exist_ok=True)
    allj = jobs()
    mine = allj[args.shard::args.nshards]
    tag = f'shard{args.shard:02d}of{args.nshards:02d}'
    csv_path = OUT / f'carry1_{tag}.csv'
    detail_path = OUT / f'cycles_{tag}.json'
    rows, details = [], {}
    for i, job in enumerate(mine):
        r = evaluate(job, args.hours)
        details[f"{job['carry1_mrna_min']}_{job['carry1_maturation_min']}"] = r.pop('cycle_metrics', [])
        rows.append(r)
        fields = sorted({k for x in rows for k in x})
        with csv_path.open('w', newline='', encoding='utf-8') as fh:
            w = csv.DictWriter(fh, fieldnames=fields, extrasaction='ignore')
            w.writeheader()
            w.writerows(rows)
        detail_path.write_text(json.dumps(details, ensure_ascii=False), encoding='utf-8')
        print(f"[{tag}] {i + 1}/{len(mine)} mRNA={job['carry1_mrna_min']:<4g} "
              f"mat={job['carry1_maturation_min']:<5g} cert={r['certified']} "
              f"steady={r.get('steady_passed')} contrast={r.get('gate_contrast_passed')} "
              f"rev={r.get('bit1_reverse_events')} gates={r.get('gate_events')} "
              f"cross={r.get('bit2_crossings')} S2min_med={r.get('s2_cycle_min_median')} "
              f"net|dS2|={r.get('abs_net_ds2_median')} ({r['runtime_s']}s)", flush=True)

    meta = dict(shard=args.shard, nshards=args.nshards, hours=args.hours,
                points=len(mine), expected_total=len(allj),
                grid=dict(carry1_mrna_min=MRNA, carry1_maturation_min=MAT),
                criterion='verify_threebit51.analyse_threebit (frozen)',
                drop_reads=DROP, sample_min=2.0, max_step_min=2.0,
                frozen_prefix_Layout='0:34 TwoBit34Model prefix',
                extension=dict(uM_per_au=5.75, bit_maturation_half_life_min=20.0,
                               carry0_maturation_half_life_min=32.5,
                               clock_K_au=0.4, clock_n=3.0, add_growth=False),
                python=sys.version, platform=platform.platform(),
                versions=dict(numpy=np.__version__, pandas=pd.__version__),
                sha256={f.replace('.py', ''): sha256(ROOT / f) for f in FROZEN_FILES},
                csv=csv_path.name, csv_sha256=sha256(csv_path))
    (OUT / f'meta_{tag}.json').write_text(json.dumps(meta, ensure_ascii=False, indent=2),
                                          encoding='utf-8')
    print(f'[{tag}] wrote {csv_path}', flush=True)


def _as_bool(s):
    return s.map(lambda x: str(x).strip().lower() in ('true', '1'))


def merge(args):
    frames = []
    for i in range(args.nshards):
        p = OUT / f'carry1_shard{i:02d}of{args.nshards:02d}.csv'
        if not p.exists():
            raise SystemExit(f'missing shard artefact {p}')
        frames.append(pd.read_csv(p))
    d = pd.concat(frames, ignore_index=True)
    for c in ('solver_success', 'finite', 'cold_passed', 'steady_passed', 'certified',
              'one_to_one', 'causal_order', 'directions_alternate', 'gate_contrast_passed'):
        if c in d.columns:
            d[c] = _as_bool(d[c])
    d = d.sort_values(['carry1_mrna_min', 'carry1_maturation_min']).reset_index(drop=True)

    expected = {(a, b) for a, b in itertools.product(MRNA, MAT)}
    got = {(float(r.carry1_mrna_min), float(r.carry1_maturation_min)) for r in d.itertuples()}
    completeness = dict(rows=len(d), unique=len(got), expected=len(expected),
                        duplicates=len(d) - len(got), missing=sorted(expected - got),
                        unexpected=sorted(got - expected))
    if completeness['rows'] != len(expected) or completeness['duplicates'] \
            or completeness['missing'] or completeness['unexpected']:
        raise SystemExit(f'incomplete scan: {completeness}')

    d.to_csv(OUT / 'carry1_all.csv', index=False, encoding='utf-8')
    summary = dict(
        completeness=completeness,
        hours=float(args.hours),
        solver_failures=int((~d.solver_success).sum()),
        nonfinite=int((~d.finite).sum()),
        certified=int(d.certified.sum()),
        steady_passed=int(d.steady_passed.sum()),
        gate_contrast_passed=int(d.gate_contrast_passed.sum()),
        one_to_one=int(d.one_to_one.sum()),
        by_mrna={str(v): dict(points=int((d.carry1_mrna_min == v).sum()),
                              certified=int(d.loc[d.carry1_mrna_min == v, 'certified'].sum()))
                 for v in MRNA},
        by_maturation={str(v): dict(points=int((d.carry1_maturation_min == v).sum()),
                                    certified=int(d.loc[d.carry1_maturation_min == v, 'certified'].sum()))
                       for v in MAT},
    )
    ranking_cols = ['carry1_mrna_min', 'carry1_maturation_min', 'crossing_error', 'bit2_crossings',
                    'bit1_reverse_events', 'gate_events', 'off_on_peak_ratio', 's2_cycle_min_median',
                    'abs_net_ds2_median', 'rev_over_fwd_median', 'bit2_low_cycles',
                    'bit2_high_cycles', 'minimum_commitment', 'steady_sequence']
    cols = [c for c in ranking_cols if c in d.columns]
    best = d.sort_values(['certified', 'crossing_error', 's2_cycle_min_median',
                          'abs_net_ds2_median', 'off_on_peak_ratio'],
                         ascending=[False, True, True, True, True]).head(10)
    summary['best_points'] = best[cols].to_dict(orient='records')
    (OUT / 'carry1_summary.json').write_text(json.dumps(summary, ensure_ascii=False,
                                                        indent=2, default=str), encoding='utf-8')
    print(json.dumps({k: v for k, v in summary.items() if k != 'best_points'},
                     ensure_ascii=False, indent=2))
    print(best[cols].to_string(index=False))
    print(f'wrote {OUT / "carry1_all.csv"}')


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='mode', required=True)
    s = sub.add_parser('scan')
    s.add_argument('--shard', type=int, required=True)
    s.add_argument('--nshards', type=int, default=7)
    s.add_argument('--hours', type=float, default=600.0)
    m = sub.add_parser('merge')
    m.add_argument('--nshards', type=int, default=7)
    m.add_argument('--hours', type=float, default=600.0)
    a = ap.parse_args()
    scan(a) if a.mode == 'scan' else merge(a)


if __name__ == '__main__':
    main()
