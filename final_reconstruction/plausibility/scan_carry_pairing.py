"""Parity-robust, direction-symmetric carry leak metric.

The defect being fixed
----------------------
`leak_ratio = median(far_off_rev) / median(gate_on_rev)` is statistically
aliased.  bit2 is written on only every other carry, so consecutive carries
alternate between two physically distinct types:

    R type : S2 is HIGH when the carry opens  -> the carry writes the bit to 0
             (the intended REVERSE action)
    F type : S2 is LOW  when the carry opens  -> the carry writes the bit to 1
             (the intended FORWARD action)

Each carry leaves the DNA in the state it just wrote, so the pool that stays
active during the following quiet dwell is the OPPOSITE direction.  The genuine
errors are therefore:

    wrong REVERSE leak  =  far_off_rev during an F-type period
    wrong FORWARD leak  =  far_off_fwd during an R-type period

A median taken over an odd number of windows lands ON one of the two types
instead of between them.  Measured on three initial states of the SAME circuit
(plausibility/eight_initial_diagnostic/paired_three_values.json): the raw ratio
varied by 2.87e+06x while the paired statistic below varied by 0.08 %.

The metric implemented here
---------------------------
Windows are classified by S2 AT THE CARRY START - the physical definition - not
by the size of the reverse dose.  One R window plus one F window form one full
8-read super-period.  Incomplete single windows at either end are DROPPED and
recorded; they are never padded with zeros.  Alle sums are formed INSIDE each
pair first, and only then are statistics taken over pairs:

    L_rev       = far_off_rev(F) / gate_on_rev(R)
    L_fwd       = far_off_fwd(R) / gate_on_fwd(F)
    L_symmetric = [far_off_rev(F) + far_off_fwd(R)]
                  / [gate_on_rev(R) + gate_on_fwd(F)]

tail is kept independent (`tail_rev_R`, `tail_fwd_F`) and is never folded back
into far-off.  The whole decomposition is repeated at the 1 %, 5 % and 10 % tail
cuts, and the old aliased ratio is recorded alongside for the erratum.

Per-window and per-pair tables are WRITTEN TO DISK for every run, so a future
re-aggregation never has to re-integrate.

    python scan_carry_pairing.py scan  --kind grid|eight --shard S --nshards N --hours 600
    python scan_carry_pairing.py merge --kind grid|eight --nshards N
"""
from __future__ import annotations

import argparse
import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plausibility_common import (OUT, frozen_carry0, frozen_extension,  # noqa: E402
                                 gate_windows, merge_shards, shard_slice,
                                 sha256, source_hashes, write_manifest,
                                 write_shard)
from model_threebit51 import ThreeBit51Model, ThreeBitCarryParameters  # noqa: E402
from model_twobit34 import CarryExpressionParameters  # noqa: E402

ROOT_MODEL = Path(__file__).resolve().parents[1] / 'model_threebit51.py'
PAIR_DIR = OUT / 'carry_pairing'
FRACTIONS = (0.01, 0.05, 0.10)
MIN_FAR_OFF_H = 0.10
S2_SPLIT = 0.5
R_DOSE_FLOOR = 1e-3      # only for the classification cross-check, never for pairing

N_VALUES = (4.0, 5.0, 6.0, 7.0, 8.0)
MRNA_VALUES = (2.0, 4.0)
MAT_VALUES = (32.5, 60.0)
SELECTED_EIGHT = dict(n_A1_gate=6.0, mrna=2.0, mat=32.5)


def grid_points():
    return [dict(n_A1_gate=n, mrna=m, mat=t)
            for n, m, t in itertools.product(N_VALUES, MRNA_VALUES, MAT_VALUES)]


def point_id(p):
    return f"n={p['n_A1_gate']:g}|mRNA={p['mrna']:g}|mat={p['mat']:g}"


def slug(p):
    return f"n{p['n_A1_gate']:g}_mRNA{p['mrna']:g}_mat{p['mat']:g}"


def build(p):
    carry1 = CarryExpressionParameters(mrna_half_life_min=p['mrna'],
                                       activator_maturation_half_life_min=p['mat'],
                                       repressor_maturation_half_life_min=p['mat'])
    return ThreeBit51Model(extension=frozen_extension(),
                           carry=ThreeBitCarryParameters(carry0=frozen_carry0(),
                                                         carry1=carry1),
                           n_A1_gate=p['n_A1_gate'])


# ------------------------------------------------------------- window segments
def segment_rows(t, J_rev2, J_fwd2, S2, I2, wins, tail_fraction,
                 min_far_off_h=MIN_FAR_OFF_H):
    """Three segments per carry window, with BOTH directions integrated."""
    rows = []
    n = len(t)
    for k, w in enumerate(wins):
        i0, i1 = w['i0'], w['i1']
        I2pk = float(np.max(I2[i0:i1 + 1])) if i1 >= i0 else 0.0
        thr = tail_fraction * I2pk
        j = i1
        while j + 1 < n and I2[j + 1] > thr:
            j += 1
        nxt = wins[k + 1]['i0'] if k + 1 < len(wins) else n - 1
        far_h = float(t[nxt] - t[j]) if nxt > j else 0.0
        ev = bool(k + 1 < len(wins) and far_h >= min_far_off_h)
        rows.append(dict(
            k=k, start_h=w['start_h'], end_h=w['end_h'], width_h=w['width_h'],
            i0=int(i0), i1=int(i1), S2_start=float(S2[i0]),
            type=('R' if S2[i0] >= S2_SPLIT else 'F'),
            gate_on_rev=float(np.trapezoid(J_rev2[i0:i1 + 1], t[i0:i1 + 1])),
            gate_on_fwd=float(np.trapezoid(J_fwd2[i0:i1 + 1], t[i0:i1 + 1])),
            tail_h=float(t[j] - t[i1]),
            tail_rev=float(np.trapezoid(J_rev2[i1:j + 1], t[i1:j + 1])),
            tail_fwd=float(np.trapezoid(J_fwd2[i1:j + 1], t[i1:j + 1])),
            far_off_h=far_h,
            far_off_rev=(float(np.trapezoid(J_rev2[j:nxt], t[j:nxt])) if ev else None),
            far_off_fwd=(float(np.trapezoid(J_fwd2[j:nxt], t[j:nxt])) if ev else None),
            far_off_evaluable=ev, I2_peak=I2pk))
    return rows


def pair_rows(rows):
    """One R + one F per super-period.  Singles at the ends are dropped, not padded."""
    pairs, dropped, anomalies = [], [], []
    i, n = 0, len(rows)
    while i < n - 1:
        a, b = rows[i], rows[i + 1]
        if a['type'] != b['type'] and a['far_off_evaluable'] and b['far_off_evaluable']:
            R = a if a['type'] == 'R' else b
            F = b if a['type'] == 'R' else a
            denom_rev = R['gate_on_rev']
            denom_fwd = F['gate_on_fwd']
            denom_sym = denom_rev + denom_fwd
            pairs.append(dict(
                pair=len(pairs), k_R=R['k'], k_F=F['k'],
                order=f"{a['type']}{b['type']}",
                gate_on_rev_R=R['gate_on_rev'], gate_on_fwd_F=F['gate_on_fwd'],
                far_off_rev_F=F['far_off_rev'], far_off_fwd_R=R['far_off_fwd'],
                tail_rev_R=R['tail_rev'], tail_fwd_F=F['tail_fwd'],
                L_rev=(F['far_off_rev'] / denom_rev if denom_rev else None),
                L_fwd=(R['far_off_fwd'] / denom_fwd if denom_fwd else None),
                L_symmetric=((F['far_off_rev'] + R['far_off_fwd']) / denom_sym
                             if denom_sym else None)))
            i += 2
        else:
            anomalies.append(dict(
                k_a=a['k'], type_a=a['type'], k_b=b['k'], type_b=b['type'],
                evaluable_a=a['far_off_evaluable'], evaluable_b=b['far_off_evaluable'],
                reason=('same type adjacent' if a['type'] == b['type']
                        else 'not evaluable')))
            i += 1
    if i == n - 1:
        dropped.append(dict(k=rows[i]['k'], type=rows[i]['type'],
                            reason='incomplete single window at the end - dropped, not padded'))
    return pairs, dropped, anomalies


def _stat(pairs, key):
    v = [p[key] for p in pairs if p[key] is not None]
    if not v:
        return dict(n=0, median=None, min=None, max=None)
    a = np.asarray(v, dtype=float)
    return dict(n=int(a.size), median=float(np.median(a)), min=float(a.min()),
                max=float(a.max()))


def analyse(model, sol, tail_fraction):
    t, y = sol.t, sol.y
    sig = model.diagnostic_signals(y)
    J_rev2, J_fwd2, S2 = sig['J_rev2'], sig['J_fwd2'], sig['S2']
    I2 = y[36]
    wins = gate_windows(t, sig['g1'], threshold=0.01, min_duration_h=0.05)
    rows = segment_rows(t, J_rev2, J_fwd2, S2, I2, wins, tail_fraction)
    pairs, dropped, anomalies = pair_rows(rows)

    nR = int(sum(r['type'] == 'R' for r in rows))
    nF = int(sum(r['type'] == 'F' for r in rows))
    # cross-check: does the S2-at-carry-start rule agree with "large reverse dose"?
    agree = int(sum((r['type'] == 'R') == (r['gate_on_rev'] > R_DOSE_FLOOR) for r in rows))

    # the OLD aliased statistic, reproduced for the erratum
    ev = [r for r in rows if r['far_off_evaluable']]
    go = [r['gate_on_rev'] for r in rows]
    old_far = [r['far_off_rev'] for r in ev]
    old_ratio = (float(np.median(old_far)) / float(np.median(go))
                 if ev and np.median(go) else None)

    agg = dict(
        tail_fraction=tail_fraction, n_windows=len(rows), n_R=nR, n_F=nF,
        n_pairs=len(pairs), n_dropped=len(dropped), n_anomalies=len(anomalies),
        windows_evaluable=len(ev), parity_balanced=bool(len(rows) % 2 == 0),
        classification_agreement=(agree / len(rows) if rows else None),
        old_aliased_leak_ratio=old_ratio,
        L_rev=_stat(pairs, 'L_rev'), L_fwd=_stat(pairs, 'L_fwd'),
        L_symmetric=_stat(pairs, 'L_symmetric'),
        gate_on_rev_R=_stat(pairs, 'gate_on_rev_R'),
        gate_on_fwd_F=_stat(pairs, 'gate_on_fwd_F'),
        far_off_rev_F=_stat(pairs, 'far_off_rev_F'),
        far_off_fwd_R=_stat(pairs, 'far_off_fwd_R'),
        tail_rev_R=_stat(pairs, 'tail_rev_R'),
        tail_fwd_F=_stat(pairs, 'tail_fwd_F'),
        pattern=''.join(r['type'] for r in rows))
    return agg, rows, pairs, dropped, anomalies


# ------------------------------------------------------------------- scanning
def evaluate(kind, p, value, hours, sample_min, max_step_min):
    from verify_threebit51 import analyse_threebit

    started = time.perf_counter()
    ident = point_id(p) if kind == 'grid' else f'init={value:03b}'
    row = dict(kind=kind, identity=ident, **p)
    try:
        model = build(p)
        if value is None:
            sol = model.simulate(hours=hours, sample_min=sample_min,
                                 max_step_min=max_step_min)
        else:
            states = json.loads((OUT / 'eight_initial_states.json').read_text(encoding='utf-8'))
            y0 = np.asarray(states['initial_states'][f'{value:03b}'], dtype=float)
            sol = model.simulate(hours=hours, sample_min=sample_min,
                                 max_step_min=max_step_min, initial_state=y0)
        a = analyse_threebit(model, sol, hours)
        row.update(ok=True, error='', certified=bool(a['certified']),
                   crossings=len(a['bit2_crossings']),
                   gates=len(a['carry1_gate_events']),
                   reverse_events=len(a['bit1_reverse_events']),
                   steady_sequence=a['steady_state']['sequence'])
        PAIR_DIR.mkdir(parents=True, exist_ok=True)
        base = slug(p) if kind == 'grid' else f'init{value:03b}'
        for f in FRACTIONS:
            agg, rows, pairs, dropped, anomalies = analyse(model, sol, f)
            pc = int(round(f * 100))
            pd.DataFrame(rows).to_csv(PAIR_DIR / f'{kind}_windows_{base}_{pc}pct.csv',
                                      index=False, encoding='utf-8')
            pd.DataFrame(pairs).to_csv(PAIR_DIR / f'{kind}_pairs_{base}_{pc}pct.csv',
                                       index=False, encoding='utf-8')
            if f == 0.05:
                row.update(dropped=json.dumps(dropped, ensure_ascii=False),
                           anomalies=json.dumps(anomalies, ensure_ascii=False))
            for key in ('n_windows', 'n_R', 'n_F', 'n_pairs', 'n_dropped', 'n_anomalies',
                        'windows_evaluable', 'parity_balanced', 'pattern',
                        'classification_agreement', 'old_aliased_leak_ratio'):
                row[f'{pc}pct_{key}'] = agg[key]
            for key in ('L_rev', 'L_fwd', 'L_symmetric', 'gate_on_rev_R',
                        'gate_on_fwd_F', 'far_off_rev_F', 'far_off_fwd_R',
                        'tail_rev_R', 'tail_fwd_F'):
                s = agg[key]
                row[f'{pc}pct_{key}_median'] = s['median']
                row[f'{pc}pct_{key}_min'] = s['min']
                row[f'{pc}pct_{key}_max'] = s['max']
                row[f'{pc}pct_{key}_n'] = s['n']
    except Exception as exc:                                     # noqa: BLE001
        row.update(ok=False, error=repr(exc))
    row['runtime_s'] = round(time.perf_counter() - started, 2)
    return row


def scan(args):
    if args.kind == 'grid':
        pts = [(p, None) for p in grid_points()]
    else:
        pts = [(dict(SELECTED_EIGHT), v) for v in range(8)]
    mine = shard_slice(pts, args.shard, args.nshards)
    tag = f'shard{args.shard:02d}of{args.nshards:02d}'
    name = f'carry_pairing_{args.kind}'
    rows = []
    for p, v in mine:
        r = evaluate(args.kind, p, v, args.hours, args.sample_min, args.max_step_min)
        rows.append(r)
        print(f"[{tag}] {r['identity']:<24} cert={r.get('certified')} "
              f"pairs={r.get('5pct_n_pairs')} R={r.get('5pct_n_R')} F={r.get('5pct_n_F')} "
              f"L_sym={r.get('5pct_L_symmetric_median')} "
              f"old={r.get('5pct_old_aliased_leak_ratio')} ({r['runtime_s']}s)", flush=True)
    write_shard(name, tag, rows, dict(kind=args.kind, shard=args.shard,
                                      nshards=args.nshards, hours=args.hours,
                                      points=len(mine), model_sha256=sha256(ROOT_MODEL),
                                      source_hashes=source_hashes()))
    print(f'[{tag}] wrote {len(rows)} rows')


def merge(args):
    name = f'carry_pairing_{args.kind}'
    df = merge_shards(name, args.nshards)
    ok = df[df.ok.astype(bool)]
    keys = ('n_pairs', 'n_R', 'n_F', 'n_dropped', 'n_anomalies', 'parity_balanced',
            'pattern', 'classification_agreement', 'old_aliased_leak_ratio',
            'L_rev_median', 'L_fwd_median', 'L_symmetric_median',
            'L_symmetric_min', 'L_symmetric_max', 'gate_on_rev_R_median',
            'gate_on_fwd_F_median', 'far_off_rev_F_median', 'far_off_fwd_R_median',
            'tail_rev_R_median', 'tail_fwd_F_median')
    verdict = dict(kind=args.kind, rows=len(df), ok=int(len(ok)),
                   model_sha256=sha256(ROOT_MODEL), source_sha256=source_hashes(),
                   metric_definition=dict(
                       classification='R if S2 at carry start >= 0.5 else F',
                       pairing='one R + one F per 8-read super-period; end singles dropped',
                       L_rev='far_off_rev(F) / gate_on_rev(R)',
                       L_fwd='far_off_fwd(R) / gate_on_fwd(F)',
                       L_symmetric='[far_off_rev(F)+far_off_fwd(R)] / '
                                   '[gate_on_rev(R)+gate_on_fwd(F)]',
                       note='all sums inside a pair first, statistics over pairs second'),
                   per_point={})

    for _, r in ok.iterrows():
        d = {'certified': bool(r.certified),
             'crossings': None if pd.isna(r.crossings) else int(r.crossings),
             'gates': None if pd.isna(r.gates) else int(r.gates),
             'reverse_events': None if pd.isna(r.reverse_events) else int(r.reverse_events),
             'steady_sequence': r.get('steady_sequence')}
        for pc in (1, 5, 10):
            for k in keys:
                v = r.get(f'{pc}pct_{k}')
                if v is None or (not isinstance(v, str) and pd.isna(v)):
                    d[f'{pc}pct_{k}'] = None
                elif isinstance(v, str):
                    d[f'{pc}pct_{k}'] = v
                elif isinstance(v, (bool, np.bool_)):
                    d[f'{pc}pct_{k}'] = bool(v)
                else:
                    d[f'{pc}pct_{k}'] = float(v)
        verdict['per_point'][str(r.identity)] = d

    if args.kind == 'grid':
        by_n = {}
        for n in N_VALUES:
            sub = ok[ok.n_A1_gate == n]
            entry = {
                'points': int(len(sub)),
                'certified': int(sub.certified.astype(bool).sum()),
                'L_rev_5pct_median': float(sub['5pct_L_rev_median'].median())
                if len(sub) else None,
                'L_fwd_5pct_median': float(sub['5pct_L_fwd_median'].median())
                if len(sub) else None,
            }
            for pc in (1, 5, 10):
                col = f'{pc}pct_L_symmetric_median'
                s = pd.to_numeric(sub[col], errors='coerce').dropna() if col in sub else None
                entry[f'L_symmetric_{pc}pct_median'] = float(s.median()) if s is not None and len(s) else None
                entry[f'L_symmetric_{pc}pct_min'] = float(s.min()) if s is not None and len(s) else None
                entry[f'L_symmetric_{pc}pct_max'] = float(s.max()) if s is not None and len(s) else None
            by_n[f'{n:g}'] = entry
        verdict['by_n_A1_gate'] = by_n
        for pc in (1, 5, 10):
            order = sorted(N_VALUES, key=lambda n: by_n[f'{n:g}'][f'L_symmetric_{pc}pct_median']
                           if by_n[f'{n:g}'][f'L_symmetric_{pc}pct_median'] is not None else 1e9)
            verdict[f'ordering_{pc}pct'] = [f'{n:g}' for n in order]
        verdict['ordering_stable_across_cuts'] = bool(
            verdict['ordering_1pct'] == verdict['ordering_5pct'] == verdict['ordering_10pct'])
        verdict['n6_beats_n5'] = bool(
            by_n['6']['L_symmetric_5pct_median'] is not None and
            by_n['5']['L_symmetric_5pct_median'] is not None and
            by_n['6']['L_symmetric_5pct_median'] < by_n['5']['L_symmetric_5pct_median'])
    else:
        vals = [d['5pct_L_symmetric_median'] for d in verdict['per_point'].values()
                if d['5pct_L_symmetric_median'] is not None]
        a = np.asarray(vals, dtype=float)
        verdict['eight_state_consistency'] = dict(
            n=len(vals), min=float(a.min()) if a.size else None,
            max=float(a.max()) if a.size else None,
            median=float(np.median(a)) if a.size else None,
            spread=(float(a.max() - a.min()) if a.size else None),
            relative_spread=(float((a.max() - a.min()) / np.median(a)) if a.size else None),
            all_certified=bool(ok.certified.astype(bool).all()),
            any_dropped_pairs=bool((ok['5pct_n_dropped'] > 0).any()),
            any_anomalies=bool((ok['5pct_n_anomalies'] > 0).any()))
    df.to_csv(OUT / f'{name}_all.csv', index=False, encoding='utf-8')
    (OUT / f'{name}_verdict.json').write_text(
        json.dumps(verdict, ensure_ascii=False, indent=2), encoding='utf-8')
    write_manifest()
    print(json.dumps(verdict, ensure_ascii=False, indent=2))


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='mode', required=True)
    s = sub.add_parser('scan')
    s.add_argument('--kind', choices=('grid', 'eight'), required=True)
    s.add_argument('--shard', type=int, required=True)
    s.add_argument('--nshards', type=int, required=True)
    s.add_argument('--hours', type=float, default=600.0)
    s.add_argument('--sample-min', type=float, default=2.0)
    s.add_argument('--max-step-min', type=float, default=2.0)
    m = sub.add_parser('merge')
    m.add_argument('--kind', choices=('grid', 'eight'), required=True)
    m.add_argument('--nshards', type=int, required=True)
    args = ap.parse_args()
    (scan if args.mode == 'scan' else merge)(args)


if __name__ == '__main__':
    main()
