"""Causal and read-window verification for the bit0 -> bit1 interface.

This verifier deliberately does not modify Model or any ZENG parameter.  It
uses three independent objects:

1. the physical bit0 reverse-recombination flux J_rev0;
2. an absolute carry-promoter opening (g0 >= 0.05), including its dose;
3. the resulting bit1 DNA-state crossing.

Digital readout is evaluated over a finite low-clock window around each flux
trough.  The same 0.3/0.7 bands and 80% occupancy rule are used for both bits.
Cold-start and post-transient (steady) verdicts are reported separately.
"""
from __future__ import annotations

import argparse
import csv
import json
from dataclasses import asdict
from pathlib import Path

import numpy as np

from model import Extension, Model, ROOT, ZENG, act, rep
from verify_bit0_part2 import CAND, OUT, clock_cycles


BAND_LOW = 0.30
BAND_HIGH = 0.70
MIN_OCCUPANCY = 0.80
READ_WINDOW_FRACTION = 0.20
STEADY_DROP_READS = 4

# Engineering acceptance limits, not fitted biochemical constants.  Keeping
# these absolute prevents a tiny numerical ripple from becoming a "carry".
JREV0_ON_PER_H = 0.10
GATE_ON = 0.05
MIN_EVENT_H = 0.20
MIN_GATE_DOSE_H = 0.02
EVENT_TOL_H = 2.0 / 60.0


def _act_array(x, K, n):
    x = np.maximum(np.asarray(x, dtype=float), 0.0)
    return x**n / (K**n + x**n)


def _segments(t, signal, threshold, min_duration_h=MIN_EVENT_H, min_area=0.0):
    """Return contiguous above-threshold events with peaks and total dose."""
    mask = np.asarray(signal) >= threshold
    changes = np.diff(np.r_[False, mask, False].astype(np.int8))
    starts = np.flatnonzero(changes == 1)
    stops = np.flatnonzero(changes == -1) - 1
    events = []
    for a, b in zip(starts, stops):
        duration = float(t[b] - t[a])
        if duration < min_duration_h:
            continue
        dose = float(np.trapezoid(signal[a:b + 1], t[a:b + 1]))
        if dose < min_area:
            continue
        p = a + int(np.argmax(signal[a:b + 1]))
        events.append(dict(i0=int(a), i1=int(b), start_h=float(t[a]),
                           end_h=float(t[b]), peak_h=float(t[p]),
                           peak=float(signal[p]), dose=dose,
                           duration_h=duration))
    return events


def _crossings(t, signal, threshold=0.5):
    """Linearly interpolated threshold crossings, including direction."""
    z = np.asarray(signal) - threshold
    ids = np.flatnonzero(z[:-1] * z[1:] < 0)
    out = []
    for i in ids:
        frac = -z[i] / (z[i + 1] - z[i])
        tc = float(t[i] + frac * (t[i + 1] - t[i]))
        out.append(dict(time_h=tc, direction='up' if z[i + 1] > z[i] else 'down'))
    return out


def _window_label(signal, a, b):
    s = np.asarray(signal[a:b + 1])
    hi = float(np.mean(s >= BAND_HIGH))
    lo = float(np.mean(s <= BAND_LOW))
    label = '1' if hi >= MIN_OCCUPANCY else ('0' if lo >= MIN_OCCUPANCY else 'x')
    return dict(label=label, high_occupancy=hi, low_occupancy=lo,
                commitment=max(hi, lo), minimum=float(s.min()),
                maximum=float(s.max()), median=float(np.median(s)))


def _read_windows(t, flux, S0, S1):
    peaks = clock_cycles(t, flux)
    rows = []
    for cycle, (a, b) in enumerate(zip(peaks[:-1], peaks[1:])):
        trough = a + int(np.argmin(flux[a:b + 1]))
        half = max(1, int(round(READ_WINDOW_FRACTION * (b - a) / 2)))
        wa, wb = max(a, trough - half), min(b, trough + half)
        l0, l1 = _window_label(S0, wa, wb), _window_label(S1, wa, wb)
        value = None if 'x' in (l0['label'], l1['label']) else 2 * int(l1['label']) + int(l0['label'])
        rows.append(dict(cycle=cycle, cycle_start_h=float(t[a]), cycle_end_h=float(t[b]),
                         trough_h=float(t[trough]), window_start_h=float(t[wa]),
                         window_end_h=float(t[wb]), bit0=l0, bit1=l1, value=value))
    return rows


def _counter_verdict(reads, drop=0):
    chosen = reads[drop:]
    vals = [r['value'] for r in chosen]
    valid = bool(len(vals) >= 8 and all(v is not None for v in vals))
    increments = bool(valid and all((vals[i + 1] - vals[i]) % 4 == 1
                                    for i in range(len(vals) - 1)))
    return dict(drop_reads=drop, reads=len(vals), valid_windows=valid,
                increments_mod4=increments, passed=bool(valid and increments),
                sequence=''.join('x' if v is None else str(v) for v in vals),
                minimum_commitment=(float(min(min(r['bit0']['commitment'],
                                                   r['bit1']['commitment']) for r in chosen))
                                    if chosen else 0.0))


def _setup_hold_margins(reads, crossings):
    times = np.asarray([c['time_h'] for c in crossings], dtype=float)
    setup, hold = [], []
    for r in reads:
        before = times[times <= r['window_start_h']]
        after = times[times >= r['window_end_h']]
        if before.size:
            setup.append(float(r['window_start_h'] - before[-1]))
        if after.size:
            hold.append(float(after[0] - r['window_end_h']))
    return dict(min_setup_h=min(setup) if setup else None,
                min_hold_h=min(hold) if hold else None)


def analyse(cfg, hours=300.0, sample_min=2.0):
    m = Model(Extension(**cfg))
    sol = m.simulate_twobit(hours=hours, sample_min=sample_min, max_step_min=2.0)
    return analyse_solution(m, sol, cfg, hours)


def analyse_solution(m, sol, cfg, hours):
    """Analyse a compatible trajectory whose first 30 states are two-bit states."""
    t, y = sol.t, sol.y
    S0, S1 = y[16], y[27]
    A0, F0 = y[28], y[29]

    expanded = [m._expand_twobit(y[:, k]) for k in range(y.shape[1])]
    flux = np.asarray([m.flux(v) for v in expanded])
    I0, R0, C0 = y[8], y[14], y[15]
    vr0 = ZENG['k_rev'] * _act_array(C0, m.K_complex, 2)
    jrev0 = vr0 * S0
    g0 = (_act_array(A0, ZENG['K_A'][0], ZENG['n_A'][0]) *
          (1.0 - _act_array(F0, ZENG['K_F'][0], ZENG['n_F'][0])))

    reverse_events = _segments(t, jrev0, JREV0_ON_PER_H)
    gate_events = _segments(t, g0, GATE_ON, min_area=MIN_GATE_DOSE_H)
    s0_cross = _crossings(t, S0)
    s1_cross = _crossings(t, S1)
    reads = _read_windows(t, flux, S0, S1)

    # Pair by the nearest major reverse episode (midpoints between successive
    # reverse peaks).  Ordering is evaluated separately so an early carry/flip
    # remains visible instead of disappearing from the association table.
    associations = []
    used_gates, used_s1 = set(), set()
    for ri, rev in enumerate(reverse_events):
        left = (0.0 if ri == 0 else
                (reverse_events[ri - 1]['peak_h'] + rev['peak_h']) / 2)
        right = (hours if ri + 1 == len(reverse_events) else
                 (rev['peak_h'] + reverse_events[ri + 1]['peak_h']) / 2)
        gates = [gi for gi, gate in enumerate(gate_events)
                 if left <= gate['peak_h'] < right]
        gi = gates[0] if len(gates) == 1 else None
        if gi is not None:
            used_gates.add(gi)
            gate = gate_events[gi]
        else:
            gate = None
        flips = [si for si, cross in enumerate(s1_cross)
                 if left <= cross['time_h'] < right]
        si = flips[0] if len(flips) == 1 else None
        if si is not None:
            used_s1.add(si)
        associations.append(dict(reverse_event=ri, gate_event=gi, bit1_crossing=si,
                                 gate_count=len(gates), bit1_crossing_count=len(flips),
                                 gate_after_reverse_start=(bool(gate['start_h'] >= rev['start_h'] - EVENT_TOL_H)
                                                           if gate else None),
                                 bit1_after_reverse_start=(bool(s1_cross[si]['time_h'] >= rev['start_h'] - EVENT_TOL_H)
                                                           if si is not None else None),
                                 gate_delay_from_reverse_start_h=(gate['start_h'] - rev['start_h']
                                                                  if gate else None),
                                 gate_delay_from_reverse_peak_h=(gate['start_h'] - rev['peak_h']
                                                                 if gate else None),
                                 bit1_delay_from_gate_h=(s1_cross[si]['time_h'] - gate['start_h']
                                                         if si is not None and gate else None),
                                 bit1_delay_from_reverse_start_h=(s1_cross[si]['time_h'] - rev['start_h']
                                                                 if si is not None else None),
                                 bit1_direction=(s1_cross[si]['direction'] if si is not None else None)))

    late = associations[2:] if len(associations) > 4 else associations
    one_to_one = bool(late and all(a['gate_count'] == 1 and a['bit1_crossing_count'] == 1 for a in late)
                      and not (set(range(len(gate_events))) - used_gates)
                      and not (set(range(len(s1_cross))) - used_s1))
    causal_order = bool(late and all(a['gate_after_reverse_start'] and a['bit1_after_reverse_start']
                                    for a in late))
    directions = [a['bit1_direction'] for a in late if a['bit1_direction']]
    alternating_directions = bool(len(directions) == len(late) and
                                  all(directions[i] != directions[i + 1]
                                      for i in range(len(directions) - 1)))

    cold = _counter_verdict(reads, 0)
    steady = _counter_verdict(reads, STEADY_DROP_READS)
    margins0 = _setup_hold_margins(reads, s0_cross)
    margins1 = _setup_hold_margins(reads, s1_cross)
    return dict(
        uM_per_au=float(cfg['uM_per_au']), hours=float(hours), add_growth=bool(cfg['add_growth']),
        thresholds=dict(band_low=BAND_LOW, band_high=BAND_HIGH,
                        minimum_occupancy=MIN_OCCUPANCY,
                        read_window_fraction=READ_WINDOW_FRACTION,
                        jrev0_on_per_h=JREV0_ON_PER_H, gate_on=GATE_ON,
                        minimum_event_h=MIN_EVENT_H, minimum_gate_dose_h=MIN_GATE_DOSE_H),
        signal_ranges=dict(S0=[float(S0.min()), float(S0.max())],
                           S1=[float(S1.min()), float(S1.max())],
                           Jrev0=[float(jrev0.min()), float(jrev0.max())],
                           g0=[float(g0.min()), float(g0.max())],
                           Int1_source_peak=float(ZENG['alpha_Int'][0] * g0.max())),
        reverse_events=reverse_events, gate_events=gate_events,
        s0_crossings=s0_cross, s1_crossings=s1_cross,
        associations=associations,
        causal_verdict=dict(exactly_one_gate_and_flip_per_late_reverse=one_to_one,
                            causal_order_after_reverse_start=causal_order,
                            bit1_directions_alternate=alternating_directions,
                            unassigned_gate_events=sorted(set(range(len(gate_events))) - used_gates),
                            unassigned_bit1_crossings=sorted(set(range(len(s1_cross))) - used_s1)),
        read_windows=reads, cold_start=cold, steady_state=steady,
        margins=dict(bit0=margins0, bit1=margins1),
        certified=bool(one_to_one and causal_order and alternating_directions and steady['passed']))


def _parse_values(text):
    if ':' in text:
        start, stop, step = map(float, text.split(':'))
        return [round(float(x), 10) for x in np.arange(start, stop + step / 2, step)]
    return [float(x) for x in text.split(',')]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--uM', default='5.5:7.5:0.25', help='comma list or start:stop:step')
    ap.add_argument('--hours', type=float, default=300.0)
    ap.add_argument('--sample-min', type=float, default=2.0)
    ap.add_argument('--output', type=Path, default=OUT / 'twobit_causal_scan.json')
    args = ap.parse_args()

    base = asdict(Extension())
    base.update(CAND)
    rows = []
    args.output.parent.mkdir(parents=True, exist_ok=True)
    for u in _parse_values(args.uM):
        cfg = base.copy()
        cfg['uM_per_au'] = u
        row = analyse(cfg, hours=args.hours, sample_min=args.sample_min)
        rows.append(row)
        args.output.write_text(json.dumps(rows, ensure_ascii=False, indent=2), encoding='utf-8')
        print(f"uM={u:.2f} cold={row['cold_start']['passed']} "
              f"steady={row['steady_state']['passed']} causal={row['causal_verdict']['exactly_one_gate_and_flip_per_late_reverse']} "
              f"order={row['causal_verdict']['causal_order_after_reverse_start']} "
              f"certified={row['certified']} seq={row['steady_state']['sequence']}", flush=True)

    summary_path = args.output.with_name(args.output.stem + '_summary.csv')
    with summary_path.open('w', newline='', encoding='utf-8-sig') as fh:
        fields = ['uM_per_au', 'cold_start_passed', 'steady_state_passed', 'causal_passed', 'causal_order',
                  'directions_alternate', 'certified', 'steady_sequence', 'min_commitment',
                  'g0_peak', 'Jrev0_peak', 'bit1_min_setup_h', 'bit1_min_hold_h']
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for r in rows:
            writer.writerow(dict(uM_per_au=r['uM_per_au'],
                                 cold_start_passed=r['cold_start']['passed'],
                                 steady_state_passed=r['steady_state']['passed'],
                                 causal_passed=r['causal_verdict']['exactly_one_gate_and_flip_per_late_reverse'],
                                 causal_order=r['causal_verdict']['causal_order_after_reverse_start'],
                                 directions_alternate=r['causal_verdict']['bit1_directions_alternate'],
                                 certified=r['certified'], steady_sequence=r['steady_state']['sequence'],
                                 min_commitment=r['steady_state']['minimum_commitment'],
                                 g0_peak=r['signal_ranges']['g0'][1],
                                 Jrev0_peak=r['signal_ranges']['Jrev0'][1],
                                 bit1_min_setup_h=r['margins']['bit1']['min_setup_h'],
                                 bit1_min_hold_h=r['margins']['bit1']['min_hold_h']))
    print(f'wrote {args.output}')
    print(f'wrote {summary_path}')


if __name__ == '__main__':
    main()
