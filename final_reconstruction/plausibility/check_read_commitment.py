"""Read-window commitment, band-edge and setup/hold margins for the carry-2 candidates.

Why this exists
---------------
The three-segment scan certified 16/20 points by the steady mod-8 rule.  That
rule already implies ">= 80 % of every read window sits inside one 0.30/0.70
band", but certification is a BINARY statement: it says nothing about how far
from the 0.80 occupancy edge the point sits, how close the DNA pools come to the
0.30/0.70 band edges, or how much setup/hold time the reads have.

`characterise()` below is the single-point instrument for all of that.  It is
shared with `scan_eight_initial_states.py` so the eight-initial-state stage uses
exactly the same definitions as this two-point margin check.

Version history (provenance - keep this block accurate)
------------------------------------------------------
v1  B508EB62CF878F285C8FB4DA52FDCAEF3DA8FD54B7BBC8CF983ABC9FD4903421
    Produced `read_commitment_shard*of02.csv`.  Its merge path had a bug: the
    `drop8_`/`drop4_` column prefix was missing, so `merge` raised KeyError.
    The scan/evaluate path was unaffected; the shard rows were correct.  A
    byte-identical copy of this source and of both v1 shard CSVs is kept under
    `plausibility/read_commitment_v1/`.
v2  769BD9DF26B8D1F234375C5187A4D1FB6D2CF6BD8504A9660DEAB72DDFCAE146
    Fixes only the merge column name.  No integration, no data change.
v3  (this file) adds: setup/hold margins for all three bits, per-bit band-edge
    margins, the authoritative `analyse_threebit` verdict fields and the
    three-segment leak metrics, so the 8-initial-state stage can reuse it.

Reproduce
---------
    python check_read_commitment.py scan  --shard 0 --nshards 2 --hours 600
    python check_read_commitment.py scan  --shard 1 --nshards 2 --hours 600
    python check_read_commitment.py merge --nshards 2

Two drops are reported: 8 (the value `verify_threebit51.counter_verdict` uses,
i.e. the one that produced the certification) and 4
(`verify_twobit_causal.STEADY_DROP_READS`, reported for contrast).

Setup/hold come from `verify_threebit51.analyse_threebit` and are computed over
ALL read windows, not only the steady set - that is the upstream definition and
it is not re-derived here.

No new acceptance threshold is invented by this script.  It reports margins; the
freeze decision stays with the reviewer.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plausibility_common import (OUT, flux_array, frozen_carry0,  # noqa: E402
                                 frozen_extension, gate_windows, leak_report,
                                 merge_shards, shard_slice, sha256,
                                 source_hashes, write_manifest, write_shard)
from model_threebit51 import ThreeBit51Model, ThreeBitCarryParameters  # noqa: E402
from model_twobit34 import CarryExpressionParameters  # noqa: E402

NAME = 'read_commitment'
ROOT_MODEL = Path(__file__).resolve().parents[1] / 'model_threebit51.py'

# The two points the reviewer asked for, verbatim.  Nothing else is scanned here.
POINTS = [dict(n_A1_gate=5.0, mrna=2.0, mat=32.5),
          dict(n_A1_gate=6.0, mrna=2.0, mat=32.5)]

PRIMARY_DROP = 8      # verify_threebit51.counter_verdict  -> produced the certification
ALT_DROP = 4          # verify_twobit_causal.STEADY_DROP_READS -> contrast only
MIN_READS = 16        # the >=16-reads rule that forces >=254 h
MIN_OCCUPANCY = 0.80
BAND_LOW = 0.30
BAND_HIGH = 0.70
REFERENCE = OUT / 'gate_segments_all.csv'


def point_id(p):
    return f"n={p['n_A1_gate']:g}|mRNA={p['mrna']:g}|mat={p['mat']:g}"


def build(p):
    """Identical construction to scan_gate_segments.build."""
    carry1 = CarryExpressionParameters(mrna_half_life_min=p['mrna'],
                                       activator_maturation_half_life_min=p['mat'],
                                       repressor_maturation_half_life_min=p['mat'])
    return ThreeBit51Model(extension=frozen_extension(),
                           carry=ThreeBitCarryParameters(carry0=frozen_carry0(),
                                                         carry1=carry1),
                           n_A1_gate=p['n_A1_gate'])


# ------------------------------------------------------------------- helpers
def _stats(values):
    v = np.asarray([x for x in values if x is not None], dtype=float)
    if v.size == 0:
        return dict(n=0, min=None, median=None, max=None)
    return dict(n=int(v.size), min=float(v.min()), median=float(np.median(v)),
                max=float(v.max()))


def _flat(prefix, s):
    return {f'{prefix}_n': s['n'], f'{prefix}_min': s['min'],
            f'{prefix}_median': s['median'], f'{prefix}_max': s['max']}


def window_metrics(reads, drop):
    """Band/occupancy margins over the steady set, for one drop value.

    Every bit gets the same treatment so the frozen bit0/bit1 columns can be
    read as extra evidence that the 34-state prefix is unaffected.  The bit2
    keys keep the exact names used by v1 so the re-run can be compared
    column-by-column against the archived shard CSVs.
    """
    chosen = reads[drop:]
    out = dict(reads_total=len(reads), reads_steady=len(chosen),
               reads_rule_satisfied=bool(len(chosen) >= MIN_READS))
    for i in range(3):
        w = [r['bits'][i] for r in chosen]
        low = [x for x in w if x['label'] == '0']
        high = [x for x in w if x['label'] == '1']
        out[f'bit{i}_low_reads'] = len(low)
        out[f'bit{i}_high_reads'] = len(high)
        out[f'bit{i}_unlabelled_reads'] = int(sum(x['label'] == 'x' for x in w))
        out.update(_flat(f'bit{i}_commitment', _stats([x['commitment'] for x in w])))
        out.update(_flat(f'bit{i}_low_commitment', _stats([x['commitment'] for x in low])))
        out.update(_flat(f'bit{i}_high_commitment', _stats([x['commitment'] for x in high])))
        out.update(_flat(f'bit{i}_low_occupancy', _stats([x['low_occupancy'] for x in low])))
        out.update(_flat(f'bit{i}_high_occupancy', _stats([x['high_occupancy'] for x in high])))
        # worst excursion of the pool inside a committed read of that state
        out.update(_flat(f'bit{i}_low_window_max', _stats([x['maximum'] for x in low])))
        out.update(_flat(f'bit{i}_high_window_min', _stats([x['minimum'] for x in high])))
        # margins: positive means strictly inside the band
        out[f'bit{i}_margin_commitment_low'] = (out[f'bit{i}_low_commitment_min'] - MIN_OCCUPANCY
                                                if low else None)
        out[f'bit{i}_margin_commitment_high'] = (out[f'bit{i}_high_commitment_min'] - MIN_OCCUPANCY
                                                 if high else None)
        out[f'bit{i}_margin_to_band_low'] = (BAND_LOW - out[f'bit{i}_low_window_max_max']
                                             if low else None)
        out[f'bit{i}_margin_to_band_high'] = (out[f'bit{i}_high_window_min_min'] - BAND_HIGH
                                              if high else None)

    # v1-compatible names: these four are the BIT2 margins
    out['margin_commitment_low'] = out['bit2_margin_commitment_low']
    out['margin_commitment_high'] = out['bit2_margin_commitment_high']
    out['margin_to_band_low'] = out['bit2_margin_to_band_low']
    out['margin_to_band_high'] = out['bit2_margin_to_band_high']
    return out


def geometry(reads, peaks, t):
    if not reads:
        return {}
    widths = np.asarray([r['window_end_h'] - r['window_start_h'] for r in reads])
    starts = np.asarray([r['window_start_h'] for r in reads])
    ends = np.asarray([r['window_end_h'] for r in reads])
    cyc_a = np.asarray([r['cycle_start_h'] for r in reads])
    cyc_b = np.asarray([r['cycle_end_h'] for r in reads])
    clipped = np.isclose(starts, cyc_a) | np.isclose(ends, cyc_b)
    period = np.median(np.diff(t[peaks])) if len(peaks) > 1 else np.nan
    return dict(clock_period_h=float(period),
                read_window_h_median=float(np.median(widths)),
                read_window_h_min=float(widths.min()),
                read_windows_clipped=int(clipped.sum()),
                read_windows_clipped_fraction=float(clipped.mean()),
                read_windows_total=int(len(reads)),
                catch_up_h=float(t[peaks[0]] - t[0]) if len(peaks) else None,
                guard_note=('a window clipped at a clock-cycle edge is evaluated on '
                            'fewer samples than the nominal 20 % width'))


def reference_row(p, hours):
    """The matching row of gate_segments_all.csv, for a cross-check only."""
    if not REFERENCE.exists() or hours != 600.0:
        return {}
    ref = pd.read_csv(REFERENCE)
    hit = ref[ref.point_id == point_id(p)]
    if hit.empty:
        return {}
    r = hit.iloc[0]
    return dict(ref_crossings=int(r.crossings), ref_gates=int(r.gates),
                ref_reverse_events=int(r.reverse_events),
                ref_certified=bool(r.certified),
                ref_leak_ratio=float(r.leak_ratio),
                ref_gate_on_rev=float(r.gate_on_rev),
                ref_bit2_hold_h=float(r.bit2_hold_h))


# ------------------------------------------------------- single-point instrument
def characterise(model, sol, hours, p=None, crosscheck=True):
    """Full single-point characterisation.  Shared with the 8-initial-state scan.

    Everything authoritative comes from `verify_threebit51.analyse_threebit`, the
    same function the certification scan used, so the definitions cannot drift.
    """
    from verify_threebit51 import analyse_threebit
    from verify_bit0_part2 import clock_cycles

    t, y = sol.t, sol.y
    a = analyse_threebit(model, sol, hours)
    flux = flux_array(model, sol)
    peaks = clock_cycles(t, flux)
    sig = model.diagnostic_signals(y)
    reads = a['read_windows']

    row = dict(hours=hours, **geometry(reads, peaks, t))

    # ---- certification / readout, identical definitions to the segment scan
    row.update(
        certified=bool(a['certified']),
        crossings=len(a['bit2_crossings']),
        gates=len(a['carry1_gate_events']),
        reverse_events=len(a['bit1_reverse_events']),
        cold_reads=a['cold_start']['reads'],
        steady_reads=a['steady_state']['reads'],
        cold_valid=a['cold_start']['valid_windows'],
        steady_valid=a['steady_state']['valid_windows'],
        cold_sequence=a['cold_start']['sequence'],
        steady_sequence=a['steady_state']['sequence'],
        steady_min_commitment_all_bits=a['steady_state']['minimum_commitment'],
        unlabelled_reads_all=int(sum(r['value'] is None for r in reads)),
        # 1:1 correspondence between bit1 reverse events, g1 windows and bit2 flips
        causal_exactly_one_gate_and_flip=bool(
            a['causal_verdict']['exactly_one_gate_and_flip_per_late_reverse']),
        causal_order=bool(a['causal_verdict']['causal_order_after_reverse_start']),
        causal_dirs_alternate=bool(a['causal_verdict']['bit2_directions_alternate']),
        causal_unassigned_gates=len(a['causal_verdict']['unassigned_gate_events']),
        causal_unassigned_crossings=len(a['causal_verdict']['unassigned_bit2_crossings']),
        gate_contrast_off_on_peak_ratio=a['gate_contrast']['off_peak_to_on_peak'],
        gate_contrast_passed=bool(a['gate_contrast']['passed']),
        gate_contrast_off_cycle_peak_max=a['gate_contrast']['off_cycle_peak_max'],
    )

    # ---- setup/hold for all three bits (ALL windows, upstream definition)
    for i in range(3):
        m = a['bit_margins'][f'bit{i}']
        row[f'bit{i}_setup_min_h'] = m['min_setup_h']
        row[f'bit{i}_hold_min_h'] = m['min_hold_h']
    setups = [row[f'bit{i}_setup_min_h'] for i in range(3)]
    holds = [row[f'bit{i}_hold_min_h'] for i in range(3)]
    row['setup_min_h_global'] = float(min(v for v in setups if v is not None)) if any(
        v is not None for v in setups) else None
    row['hold_min_h_global'] = float(min(v for v in holds if v is not None)) if any(
        v is not None for v in holds) else None

    # ---- band/occupancy margins at both drops
    for drop, tag in ((PRIMARY_DROP, 'drop8'), (ALT_DROP, 'drop4')):
        for k, v in window_metrics(reads, drop).items():
            row[f'{tag}_{k}'] = v

    # ---- three-segment leak decomposition
    win = gate_windows(t, sig['g1'])
    rep = leak_report(t, sig['g1'], sig['J_rev2'], sig['J_fwd2'], y[36], win)
    for f, tag in ((0.01, 'leak1pct'), (0.05, 'leak5pct'), (0.10, 'leak10pct')):
        rr = rep[f]
        row[f'{tag}_gate_on_rev'] = rr['gate_on_rev']
        row[f'{tag}_tail_rev'] = rr['tail_rev']
        row[f'{tag}_far_off_rev'] = rr['far_off_rev']
        row[f'{tag}_far_off_fwd'] = rr['far_off_fwd']
        row[f'{tag}_ratio'] = rr['leak_ratio']
        row[f'{tag}_gate_on_h'] = rr['gate_on_h']
        row[f'{tag}_tail_h'] = rr['tail_h']
        row[f'{tag}_far_off_h'] = rr['far_off_h']
        row[f'{tag}_windows_total'] = rr['windows_total']
        row[f'{tag}_windows_evaluable'] = rr['windows_evaluable']

    # ---- v1-compatible aliases (recomputed here, compared against the CSV)
    row.update(ref_crossings=row['crossings'], ref_gates=row['gates'],
               ref_leak_ratio=row['leak5pct_ratio'],
               ref_gate_on_rev=row['leak5pct_gate_on_rev'],
               ref_bit2_hold_h=row['bit2_hold_min_h'])
    if crosscheck and p is not None:
        ref = reference_row(p, hours)
        for k, v in ref.items():
            row[f'crosscheck_{k}'] = v
        if 'ref_leak_ratio' in ref:
            row['crosscheck_leak_ratio_gap'] = abs(row['ref_leak_ratio'] - ref['ref_leak_ratio'])
            row['crosscheck_hold_gap'] = abs(row['ref_bit2_hold_h'] - ref['ref_bit2_hold_h'])
        if 'ref_crossings' in ref:
            row['crosscheck_crossings_gap'] = abs(row['ref_crossings'] - ref['ref_crossings'])
    return row


def evaluate(p, hours, sample_min, max_step_min):
    started = time.perf_counter()
    row = dict(point_id=point_id(p), **p)
    try:
        model = build(p)
        sol = model.simulate(hours=hours, sample_min=sample_min,
                             max_step_min=max_step_min)
        row.update(ok=True, error='')
        row.update(characterise(model, sol, hours, p))
    except Exception as exc:                                     # noqa: BLE001
        row.update(ok=False, error=repr(exc))
    row['runtime_s'] = round(time.perf_counter() - started, 2)
    return row


def scan(args):
    pts = shard_slice(POINTS, args.shard, args.nshards)
    tag = f'shard{args.shard:02d}of{args.nshards:02d}'
    rows = []
    for p in pts:
        r = evaluate(p, args.hours, args.sample_min, args.max_step_min)
        rows.append(r)
        print(f"[{tag}] {point_id(p):<24} ok={r['ok']} cert={r.get('certified')} "
              f"reads={r.get('drop8_reads_steady')} "
              f"bit2_low_commit_min={r.get('drop8_bit2_low_commitment_min')} "
              f"bit2_high_commit_min={r.get('drop8_bit2_high_commitment_min')} "
              f"setup={r.get('setup_min_h_global')} hold={r.get('hold_min_h_global')} "
              f"({r['runtime_s']}s)", flush=True)
    write_shard(NAME, tag, rows, dict(kind='check', shard=args.shard, nshards=args.nshards,
                                      hours=args.hours, points=len(pts),
                                      model_sha256=sha256(ROOT_MODEL),
                                      steady_drop=PRIMARY_DROP, alt_drop=ALT_DROP,
                                      source_hashes=source_hashes()))
    print(f'[{tag}] wrote {len(rows)} rows')


MARGIN_KEYS = (('margin_commitment_low', 'margin_commitment_high',
                'margin_to_band_low', 'margin_to_band_high',
                'bit2_low_commitment_min', 'bit2_high_commitment_min',
                'bit2_low_window_max_max', 'bit2_high_window_min_min') +
               tuple(f'{b}_{k}' for b in ('bit0', 'bit1')
                     for k in ('commitment_min', 'low_commitment_min', 'high_commitment_min',
                               'margin_to_band_low', 'margin_to_band_high',
                               'low_window_max_max', 'high_window_min_min')))
D4_KEYS = ('bit2_low_commitment_min', 'bit2_high_commitment_min',
           'margin_to_band_low', 'margin_to_band_high')


def merge(args):
    df = merge_shards(NAME, args.nshards)
    ok = df[df.ok.astype(bool)].sort_values('n_A1_gate')
    verdict = dict(rows=len(df), ok=int(len(ok)), primary_drop=PRIMARY_DROP,
                   alt_drop=ALT_DROP, min_reads_rule=MIN_READS,
                   min_occupancy_rule=MIN_OCCUPANCY,
                   bands=dict(low=BAND_LOW, high=BAND_HIGH),
                   source_sha256=source_hashes(),
                   model_sha256=sha256(ROOT_MODEL),
                   note=('Margins are reported, not thresholded. A negative margin means the '
                         'point violates that specific edge; no new acceptance threshold is '
                         'introduced by this script.'))

    scalar = ('certified', 'crossings', 'gates', 'reverse_events', 'clock_period_h',
              'read_window_h_median', 'read_windows_clipped', 'read_windows_total',
              'cold_reads', 'steady_reads', 'steady_min_commitment_all_bits',
              'unlabelled_reads_all', 'causal_exactly_one_gate_and_flip', 'causal_order',
              'causal_dirs_alternate', 'causal_unassigned_gates',
              'causal_unassigned_crossings', 'gate_contrast_passed',
              'gate_contrast_off_on_peak_ratio', 'setup_min_h_global', 'hold_min_h_global',
              'bit0_setup_min_h', 'bit0_hold_min_h', 'bit1_setup_min_h', 'bit1_hold_min_h',
              'bit2_setup_min_h', 'bit2_hold_min_h', 'steady_sequence',
              'crosscheck_leak_ratio_gap', 'crosscheck_crossings_gap', 'crosscheck_hold_gap',
              'leak5pct_ratio', 'leak5pct_far_off_rev', 'leak5pct_far_off_fwd',
              'leak1pct_ratio', 'leak10pct_ratio')

    per_point = {}
    for _, r in ok.iterrows():
        n = float(r.n_A1_gate)
        d = {'reads_steady_drop8': int(r['drop8_reads_steady']),
             'reads_rule_satisfied': bool(r['drop8_reads_rule_satisfied']),
             'reads_steady_drop4': int(r['drop4_reads_steady']),
             'bit2_low_reads': int(r['drop8_bit2_low_reads']),
             'bit2_high_reads': int(r['drop8_bit2_high_reads']),
             'bit2_unlabelled_reads': int(r['drop8_bit2_unlabelled_reads'])}
        for k in scalar:
            v = r.get(k)
            if v is None or (not isinstance(v, str) and pd.isna(v)):
                d[k] = None
            elif isinstance(v, str):
                d[k] = v
            elif isinstance(v, (bool, np.bool_)):
                d[k] = bool(v)
            else:
                d[k] = float(v)
        for k in MARGIN_KEYS:
            v = r.get(f'drop{PRIMARY_DROP}_{k}')
            d[k] = None if v is None or pd.isna(v) else float(v)
        for k in D4_KEYS:
            v = r.get(f'drop{ALT_DROP}_{k}')
            d[f'drop4_{k}'] = None if v is None or pd.isna(v) else float(v)
        per_point[f'n={n:g}'] = d
    verdict['per_point'] = per_point

    if len(ok) == 2:
        a = per_point[[k for k in per_point if k.startswith('n=5')][0]]
        b = per_point[[k for k in per_point if k.startswith('n=6')][0]]
        cmp_keys = ('bit2_low_commitment_min', 'bit2_high_commitment_min',
                    'margin_commitment_low', 'margin_commitment_high',
                    'margin_to_band_low', 'margin_to_band_high',
                    'bit2_low_window_max_max', 'bit2_high_window_min_min',
                    'bit0_setup_min_h', 'bit0_hold_min_h',
                    'bit1_setup_min_h', 'bit1_hold_min_h',
                    'bit2_setup_min_h', 'bit2_hold_min_h',
                    'setup_min_h_global', 'hold_min_h_global')
        cmp_rows, wins = {}, []
        for k in cmp_keys:
            va, vb = a.get(k), b.get(k)
            if va is None or vb is None:
                cmp_rows[k] = dict(n5=va, n6=vb, larger=None, difference=None)
                continue
            bigger = 'n=6' if vb > va else ('n=5' if va > vb else 'tie')
            cmp_rows[k] = dict(n5=float(va), n6=float(vb), larger=bigger,
                               difference=float(vb - va))
            wins.append(bigger)
        verdict['margin_comparison'] = cmp_rows
        verdict['larger_count'] = dict(n5=wins.count('n=5'), n6=wins.count('n=6'),
                                       tie=wins.count('tie'))
        verdict['margin_note'] = ('For band/occupancy margins a LARGER number is farther from '
                                  'the edge; for bit2_low_window_max_max a SMALLER number is '
                                  'farther from the 0.30 edge. Read the row, not the tally, '
                                  'for that one.')

    df.to_csv(OUT / f'{NAME}_all.csv', index=False, encoding='utf-8')
    (OUT / f'{NAME}_verdict.json').write_text(
        json.dumps(verdict, ensure_ascii=False, indent=2), encoding='utf-8')
    write_manifest()
    print(json.dumps(verdict, ensure_ascii=False, indent=2))


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='mode', required=True)
    s = sub.add_parser('scan')
    s.add_argument('--shard', type=int, required=True)
    s.add_argument('--nshards', type=int, required=True)
    s.add_argument('--hours', type=float, default=600.0)
    s.add_argument('--sample-min', type=float, default=2.0)
    s.add_argument('--max-step-min', type=float, default=2.0)
    m = sub.add_parser('merge')
    m.add_argument('--nshards', type=int, required=True)
    args = ap.parse_args()
    if args.mode == 'scan':
        scan(args)
    else:
        merge(args)


if __name__ == '__main__':
    main()
