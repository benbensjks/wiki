"""Per-window diagnostic for the eight-initial-state leak anomaly.

Two things in `eight_initial_all.csv` need explaining before any of it is quoted:

  * values 3 (011) and 7 (111) produce 15 gate windows / 15 reverse events /
    15 bit2 crossings over the same 600 h, where the other six produce 14;
  * value 3 has `leak5pct_gate_on_rev` = 3.46e-07 against ~0.4975 for the other
    five 14-window runs, i.e. the reverse dose inside its gate windows is
    essentially zero, which is what turns the ratio into 34705.

This script re-runs one initial state and prints the raw tables the summary row
is derived from: the read windows, the two different gate-window definitions
(threshold 0.01 used by the leak decomposition, threshold 0.05 used by
`analyse_threebit`), the reverse events, and the per-window three-segment table.

    python diagnose_eight_leak.py --value 3 --hours 600
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plausibility_common import OUT, flux_array, gate_windows, leak_report  # noqa: E402
from scan_eight_initial_states import (STATES_JSON, build_selected,  # noqa: E402
                                       decode_initial)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--value', type=int, required=True)
    ap.add_argument('--hours', type=float, default=600.0)
    ap.add_argument('--sample-min', type=float, default=2.0)
    ap.add_argument('--max-step-min', type=float, default=2.0)
    ap.add_argument('--tag', type=str, default='')
    args = ap.parse_args()

    from verify_threebit51 import analyse_threebit
    from verify_bit0_part2 import clock_cycles

    payload = json.loads(STATES_JSON.read_text(encoding='utf-8'))
    y0 = np.asarray(payload['initial_states'][f'{args.value:03b}'], dtype=float)
    model = build_selected()
    sol = model.simulate(hours=args.hours, sample_min=args.sample_min,
                         max_step_min=args.max_step_min, initial_state=y0)
    t, y = sol.t, sol.y
    a = analyse_threebit(model, sol, args.hours)
    sig = model.diagnostic_signals(y)
    flux = flux_array(model, sol)
    peaks = clock_cycles(t, flux)

    outdir = OUT / 'eight_initial_diagnostic'
    outdir.mkdir(parents=True, exist_ok=True)
    tag = args.tag or f'{args.value:03b}'

    print(f'=== value {args.value:03b}  decoded={decode_initial(y0)} '
          f'certified={a["certified"]} ===')
    print(f'reverse events {len(a["bit1_reverse_events"])}  '
          f'gate events(0.05) {len(a["carry1_gate_events"])}  '
          f'bit2 crossings {len(a["bit2_crossings"])}  '
          f'read windows {len(a["read_windows"])}')
    print(f'causal verdict: {json.dumps(a["causal_verdict"], ensure_ascii=False)}')
    print('gate_contrast:', json.dumps({k: v for k, v in a['gate_contrast'].items()
                                        if k != 'passed'}, ensure_ascii=False))

    reads = pd.DataFrame([dict(cycle=r['cycle'], trough_h=r['trough_h'],
                               value=r['value'],
                               S0=r['bits'][0]['label'], S1=r['bits'][1]['label'],
                               S2=r['bits'][2]['label'],
                               win_h=r['window_end_h'] - r['window_start_h'])
                          for r in a['read_windows']])
    reads.to_csv(outdir / f'reads_{tag}.csv', index=False, encoding='utf-8')
    print('\n--- read windows ---')
    print(reads.to_string(index=False))

    rev = a['bit1_reverse_events']
    ge = a['carry1_gate_events']
    win01 = gate_windows(t, sig['g1'], threshold=0.01, min_duration_h=0.05)
    win05 = gate_windows(t, sig['g1'], threshold=0.05, min_duration_h=0.05)
    print(f'\n--- event counts ---')
    print(f'reverse events            {len(rev)}')
    print(f'gate episodes thr=0.05    {len(ge)}')
    print(f'gate windows  thr=0.05    {len(win05)}')
    print(f'gate windows  thr=0.01    {len(win01)}   <- used by the leak decomposition')
    print(f'bit2 crossings            {len(a["bit2_crossings"])}')

    rows = []
    for k, w in enumerate(win01):
        rows.append(dict(k=k, start_h=w['start_h'], end_h=w['end_h'],
                         width_h=w['width_h'], peak=w['peak'], dose=w['dose'],
                         i0=w['i0'], i1=w['i1'], samples=w['i1'] - w['i0'] + 1))
    gw = pd.DataFrame(rows)
    gw.to_csv(outdir / f'gate_windows_{tag}.csv', index=False, encoding='utf-8')
    print('\n--- gate windows (threshold 0.01, the leak decomposition input) ---')
    print(gw.to_string(index=False))

    rep = leak_report(t, sig['g1'], sig['J_rev2'], sig['J_fwd2'], y[36], win01)
    for f in (0.01, 0.05, 0.10):
        tbl = pd.DataFrame(rep[f]['per_window'])
        tbl.to_csv(outdir / f'segments_{int(f*100)}pct_{tag}.csv', index=False,
                   encoding='utf-8')
    print('\n--- three-segment table at the 5 % cut ---')
    tbl = pd.DataFrame(rep[0.05]['per_window'])
    pd.set_option('display.width', 260)
    print(tbl.to_string(index=False))
    print('\n--- aggregate at each cut ---')
    for f in (0.01, 0.05, 0.10):
        r = rep[f]
        print(f'  {int(f*100):>2}%  windows={r["windows_total"]:>3} '
              f'evaluable={r["windows_evaluable"]:>3} '
              f'gate_on_rev={r["gate_on_rev"]!r} tail_rev={r["tail_rev"]!r} '
              f'far_off_rev={r["far_off_rev"]!r} far_off_fwd={r["far_off_fwd"]!r} '
              f'ratio={r["leak_ratio"]!r}')

    # does each gate window actually contain the reverse event it is paired with?
    rev_peak = np.asarray([e['peak_h'] for e in rev])
    print('\n--- pairing check: reverse-event peak inside which gate window ---')
    pair = []
    for e in rev:
        inside = [k for k, w in enumerate(win01) if w['start_h'] <= e['peak_h'] <= w['end_h']]
        pair.append(dict(reverse_peak_h=e['peak_h'], reverse_dose=e['dose'],
                         gate_windows_containing=inside))
    print(pd.DataFrame(pair).to_string(index=False))
    print(f'\nwrote diagnostic tables to {outdir}')


if __name__ == '__main__':
    main()
