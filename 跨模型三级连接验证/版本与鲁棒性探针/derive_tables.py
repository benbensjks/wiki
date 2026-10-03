"""Derive the report CSVs from the stored JSON verdicts, with no model runs.

NEW FILE (added 2026-09-29).  Written because the first-round CSVs flattened a
missing quantity into a number: in the two failed version-probe cases bit1 and
bit2 NEVER CROSS 0.5, so their setup/hold margins are undefined (the JSON
correctly stores null), and the aggregate `binding_min_hold_h` column then
silently reported BIT0's margin (4.2388 h / 3.3089 h).  Read on its own, that
column looks like a healthy margin.  Same failure class the project has already
been bitten by: "which bit is the binding one".

This script makes the derived tables self-describing:
  * `margins_defined_for`  -- which bits actually have a margin, e.g. "bit0"
  * `all_bits_crossed`     -- True only when all three bits have one
  * `binding_is_the_binding_bit` -- True only when the reported binding value is
                              not a fallback onto a non-binding bit
  * `bit*_hold_h` / `bit*_setup_h` -- explicit per-bit columns, empty when null

It also RE-DERIVES not a single new number: every value comes from the JSON
already on disk.  If the JSON and the CSV ever disagree, the JSON is the
authority; this script is a view, not a measurement.

Usage:  python derive_tables.py
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')

HERE = Path(__file__).resolve().parent
RES = HERE / 'results'

CHECKS = []


def check(ok, msg):
    CHECKS.append((bool(ok), msg))
    print(('  OK   ' if ok else '  FAIL ') + msg, flush=True)


def load(name):
    return json.loads((RES / name).read_text(encoding='utf-8'))


def cell(v):
    return '' if v is None else v


def version_probe_table():
    d = load('probe_B_zeng_version_51state.json')
    rows = d['rows']
    cols = ['tag', 'certified', 'steady_passed', 'steady_reads', 'minimum_commitment',
            'bit0_hold_h', 'bit1_hold_h', 'bit2_hold_h',
            'bit0_setup_h', 'bit1_setup_h', 'bit2_setup_h',
            'margins_defined_for', 'all_bits_crossed', 'binding_min_hold_h',
            'binding_is_the_binding_bit', 'carry1_gate_events', 'bit1_reverse_events',
            'bit2_crossings', 'exactly_one', 'causal_order', 'alternating',
            'off_to_on_peak_ratio', 'runtime_s']
    out = []
    for r in rows:
        m = r['margins']
        defined = [b for b in ('bit0', 'bit1', 'bit2')
                   if m.get(b, {}).get('min_hold_h') is not None]
        out.append({
            'tag': r['tag'], 'certified': r['certified'],
            'steady_passed': r['steady_passed'], 'steady_reads': r['steady_reads'],
            'minimum_commitment': r['minimum_commitment'],
            'bit0_hold_h': cell(m.get('bit0', {}).get('min_hold_h')),
            'bit1_hold_h': cell(m.get('bit1', {}).get('min_hold_h')),
            'bit2_hold_h': cell(m.get('bit2', {}).get('min_hold_h')),
            'bit0_setup_h': cell(m.get('bit0', {}).get('min_setup_h')),
            'bit1_setup_h': cell(m.get('bit1', {}).get('min_setup_h')),
            'bit2_setup_h': cell(m.get('bit2', {}).get('min_setup_h')),
            'margins_defined_for': '+'.join(defined) if defined else 'none',
            'all_bits_crossed': len(defined) == 3,
            'binding_min_hold_h': cell(m.get('binding_min_hold_h')),
            # The binding bit is the one whose hold margin was the frozen constraint:
            # bit1.  If bit1 has no margin at all, the aggregate column is a fallback.
            'binding_is_the_binding_bit': m.get('bit1', {}).get('min_hold_h') is not None,
            'carry1_gate_events': r['carry1_gate_events'],
            'bit1_reverse_events': r['bit1_reverse_events'],
            'bit2_crossings': r['bit2_crossings'],
            'exactly_one': r['exactly_one'], 'causal_order': r['causal_order'],
            'alternating': r['alternating'],
            'off_to_on_peak_ratio': cell(r['off_to_on_peak_ratio']),
            'runtime_s': r['runtime_s']})
    write_csv('probe_B_zeng_version_51state.csv', cols, out)
    check(out[0]['all_bits_crossed'] is True, 'version probe baseline: all three bits crossed')
    check(out[1]['all_bits_crossed'] is False and out[1]['margins_defined_for'] == 'bit0',
          'zmh_bit1_block: only bit0 has a margin, and the table now says so')
    check(out[1]['binding_is_the_binding_bit'] is False,
          'zmh_bit1_block: aggregate margin flagged as a fallback onto bit0')
    check(out[1]['bit1_hold_h'] == '' and out[1]['bit2_hold_h'] == '',
          'zmh_bit1_block: undefined margins serialise as empty, not as a number')
    return out


def hybrid_table():
    d = load('probe_P1_P4_hybrid.json')
    rows = d['rows']
    cols = ['tag', 'certified_v1', 'steady_reads', 'boundary_clips',
            'minimum_commitment', 'global_min_timing_margin_h',
            'bitS0_hold_h', 'bitS1_hold_h', 'bitS2_hold_h',
            'margins_defined_for', 'clock_peak_count', 'chain01_passed', 'chain12_passed',
            'chain01_reverse_events', 'chain12_reverse_events', 'duty_g1_gt_0p1', 'g1_peak',
            'runtime_s']
    out = []
    for r in rows:
        bm = r['bit_margins']
        defined = [b for b in ('S0', 'S1', 'S2') if bm.get(b, {}).get('min_hold_h') is not None]
        ev = r['events']
        out.append({
            'tag': r['tag'], 'certified_v1': r['certified_v1'],
            'steady_reads': r['steady_reads'], 'boundary_clips': r['boundary_clips'],
            'minimum_commitment': r['minimum_commitment'],
            'global_min_timing_margin_h': cell(r['global_min_timing_margin_h']),
            'bitS0_hold_h': cell(bm.get('S0', {}).get('min_hold_h')),
            'bitS1_hold_h': cell(bm.get('S1', {}).get('min_hold_h')),
            'bitS2_hold_h': cell(bm.get('S2', {}).get('min_hold_h')),
            'margins_defined_for': '+'.join(defined) if defined else 'none',
            'clock_peak_count': r['clock_peak_count'],
            'chain01_passed': ev['bit0_to_bit1']['passed'],
            'chain12_passed': ev['bit1_to_bit2']['passed'],
            'chain01_reverse_events': ev['bit0_to_bit1']['reverse_events'],
            'chain12_reverse_events': ev['bit1_to_bit2']['reverse_events'],
            'duty_g1_gt_0p1': r['duty_g1_gt_0p1'], 'g1_peak': r['g1_peak'],
            'runtime_s': r['runtime_s']})
    write_csv('probe_P1_P4_hybrid.csv', cols, out)
    check(all(r['margins_defined_for'] == 'S0+S1+S2' for r in out),
          'hybrid: every case has all three bit margins defined')
    check(all(r['global_min_timing_margin_h'] not in ('', None) for r in out),
          'hybrid: global minimum timing margin is defined in every case')
    return out


def write_csv(name, cols, rows):
    path = RES / name
    with path.open('w', encoding='utf-8', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f'  wrote {path.name} ({len(rows)} rows, {len(cols)} columns)')


def main():
    print('Deriving report tables from the stored JSON verdicts')
    version_probe_table()
    hybrid_table()
    bad = [m for ok, m in CHECKS if not ok]
    print(f'\n{len(CHECKS) - len(bad)}/{len(CHECKS)} checks passed')
    if bad:
        raise SystemExit('FAILED CHECKS:\n' + '\n'.join(bad))


if __name__ == '__main__':
    main()
