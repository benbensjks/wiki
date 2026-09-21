"""Verify all four phase-consistent digital initial states of TwoBit34Model.

The initial states are not made by changing S alone.  They are sampled from
the late stable period-4 orbit, including every biochemical pool, then their
six oscillator states plus shared C31 mRNA are reset to one common clock phase.
"""
from __future__ import annotations

import argparse
import json
from dataclasses import asdict, replace
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from model import ROOT
from model_twobit34 import (STATE_NAMES_34, CarryExpressionParameters,
                            TwoBit34Model, nominal_extension)
from verify_twobit_causal import analyse_solution


OUT = ROOT / 'twobit34_results' / 'initial_states'


def extract_orbit_states(model, warmup_h=300.0, cutoff_h=180.0):
    sol = model.simulate(hours=warmup_h, sample_min=2.0, max_step_min=2.0)
    analysis = analyse_solution(model.base, sol, asdict(model.e), warmup_h)
    selected = {}
    for read in analysis['read_windows']:
        value = read['value']
        if value is not None and read['trough_h'] >= cutoff_h and value not in selected:
            k = int(np.argmin(np.abs(sol.t - read['trough_h'])))
            selected[int(value)] = dict(time_h=float(sol.t[k]), index=k,
                                        read=read, state=sol.y[:, k].copy())
    if set(selected) != {0, 1, 2, 3}:
        raise RuntimeError(f'could not extract all four values: {sorted(selected)}')

    # All flux troughs represent the same oscillator phase modulo one cycle.
    # Use one actually simulated upstream state rather than their numerical mean.
    common_upstream7 = selected[0]['state'][:7].copy()
    raw = np.vstack([selected[v]['state'][:7] for v in range(4)])
    scale = np.maximum(np.max(np.abs(raw), axis=0), 1e-12)
    phase_spread_abs = np.ptp(raw, axis=0)
    phase_spread_rel = phase_spread_abs / scale
    for value in range(4):
        selected[value]['state_before_phase_reset'] = selected[value]['state'].copy()
        selected[value]['state'][:7] = common_upstream7
    audit = dict(common_phase_source_value=0,
                 sampled_times_h=[selected[v]['time_h'] for v in range(4)],
                 max_upstream_absolute_spread=float(phase_spread_abs.max()),
                 max_upstream_relative_spread=float(phase_spread_rel.max()),
                 upstream_absolute_spread=dict(zip(STATE_NAMES_34[:7], map(float, phase_spread_abs))),
                 upstream_relative_spread=dict(zip(STATE_NAMES_34[:7], map(float, phase_spread_rel))))
    return selected, audit


def decode_initial(y):
    b0 = int(y[16] >= .5); b1 = int(y[27] >= .5)
    return 2*b1 + b0


def expected_sequence(initial_value, length):
    # Initial states are sampled at a read trough.  The first newly observable
    # read after integration is one clock period later.
    return ''.join(str((initial_value + i + 1) % 4) for i in range(length))


def plot_sequences(results, out_dir):
    fig, axes = plt.subplots(4, 1, figsize=(12, 8), sharex=True, constrained_layout=True)
    colors = ('#6A3D9A', '#1F78B4', '#33A02C', '#E31A1C')
    for value, ax in enumerate(axes):
        reads = results[value]['read_windows']
        t = np.asarray([r['trough_h'] for r in reads])
        observed = np.asarray([r['value'] for r in reads], dtype=float)
        expected = np.asarray([(value + i + 1) % 4 for i in range(len(reads))])
        ax.step(t, expected, where='mid', color='0.65', ls='--', lw=2,
                label='expected')
        ax.step(t, observed, where='mid', color=colors[value], lw=1.5,
                label='observed')
        ax.scatter(t, observed, color=colors[value], s=18, zorder=3)
        ax.set_yticks((0, 1, 2, 3)); ax.set_ylim(-.35, 3.35)
        ax.set_ylabel(f'init {value:02b}')
        ax.legend(loc='upper right', ncol=2, fontsize=8)
    axes[0].set_title('Four phase-consistent initial states: expected vs observed mod-4 count')
    axes[-1].set_xlabel('time after restart (h)')
    fig.savefig(out_dir / 'four_initial_states.png', dpi=220, bbox_inches='tight', facecolor='white')
    fig.savefig(out_dir / 'four_initial_states.pdf', bbox_inches='tight', facecolor='white')
    plt.close(fig)


def plot_initial_composition(selected, out_dir):
    matrix = np.vstack([selected[v]['state'] for v in range(4)]).T
    scale = np.maximum(np.max(np.abs(matrix), axis=1, keepdims=True), 1e-12)
    normalized = matrix / scale
    fig, ax = plt.subplots(figsize=(8, 12), constrained_layout=True)
    im = ax.imshow(normalized, aspect='auto', cmap='viridis', vmin=0, vmax=1)
    ax.set_xticks(range(4), ('00', '01', '10', '11'))
    ax.set_yticks(range(34), STATE_NAMES_34, fontsize=7)
    ax.set_xlabel('initial digital state'); ax.set_title(
        'Self-consistent 34-state initial conditions\n(each molecular state normalized across four branches)')
    fig.colorbar(im, ax=ax, shrink=.6, label='normalized state value')
    fig.savefig(out_dir / 'initial_state_composition.png', dpi=220, bbox_inches='tight', facecolor='white')
    fig.savefig(out_dir / 'initial_state_composition.pdf', bbox_inches='tight', facecolor='white')
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--uM', type=float, default=6.0)
    ap.add_argument('--carry-mat-min', type=float, default=30.0)
    ap.add_argument('--carry-mrna-min', type=float, default=2.0)
    ap.add_argument('--warmup-h', type=float, default=300.0)
    ap.add_argument('--hours', type=float, default=300.0)
    ap.add_argument('--out-dir', type=Path, default=OUT)
    args = ap.parse_args()

    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    extension = replace(nominal_extension(), uM_per_au=args.uM)
    carry = CarryExpressionParameters(
        mrna_half_life_min=args.carry_mrna_min,
        activator_maturation_half_life_min=args.carry_mat_min,
        repressor_maturation_half_life_min=args.carry_mat_min)
    model = TwoBit34Model(extension, carry)
    selected, phase_audit = extract_orbit_states(model, args.warmup_h)

    state_table = {'state': STATE_NAMES_34}
    for value in range(4):
        state_table[f'initial_{value:02b}'] = selected[value]['state']
    pd.DataFrame(state_table).to_csv(out_dir / 'four_initial_states.csv', index=False)

    results = {}
    for value in range(4):
        y0 = selected[value]['state']
        decoded = decode_initial(y0)
        sol = model.simulate(hours=args.hours, sample_min=2.0, max_step_min=2.0,
                             initial_state=y0)
        analysis = analyse_solution(model.base, sol, asdict(extension), args.hours)
        observed = analysis['cold_start']['sequence']
        expected = expected_sequence(value, len(observed))
        relation_ok = decoded == value and observed == expected
        result = dict(initial_value=value, initial_bits=f'{value:02b}',
                      source_time_h=selected[value]['time_h'], decoded_initial=decoded,
                      expected_sequence=expected, observed_sequence=observed,
                      sequence_relation_correct=relation_ok,
                      cold_start=analysis['cold_start'], steady_state=analysis['steady_state'],
                      certified=analysis['certified'], causal_verdict=analysis['causal_verdict'],
                      margins=analysis['margins'], signal_ranges=analysis['signal_ranges'],
                      reverse_events=len(analysis['reverse_events']),
                      gate_events=len(analysis['gate_events']),
                      read_windows=analysis['read_windows'])
        result['passed'] = bool(relation_ok and analysis['cold_start']['passed'] and
                                analysis['certified'])
        results[value] = result
        print(f"initial={value:02b} decoded={decoded:02b} passed={result['passed']} "
              f"observed={observed}", flush=True)

    report = dict(parameters=dict(extension=asdict(extension), carry=asdict(carry),
                                  warmup_h=args.warmup_h, verification_h=args.hours),
                  construction=('late orbit-derived complete biochemical states; '
                                'oscillator plus C31 mRNA reset to one common phase'),
                  phase_audit=phase_audit,
                  results={str(k): v for k, v in results.items()},
                  all_four_pass=bool(all(r['passed'] for r in results.values())))
    (out_dir / 'four_initial_states.json').write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
    plot_sequences(results, out_dir)
    plot_initial_composition(selected, out_dir)

    rows = []
    for value, r in results.items():
        rows.append(dict(initial=f'{value:02b}', passed=r['passed'],
                         expected=r['expected_sequence'], observed=r['observed_sequence'],
                         bit0_setup_h=r['margins']['bit0']['min_setup_h'],
                         bit0_hold_h=r['margins']['bit0']['min_hold_h'],
                         bit1_setup_h=r['margins']['bit1']['min_setup_h'],
                         bit1_hold_h=r['margins']['bit1']['min_hold_h']))
    pd.DataFrame(rows).to_csv(out_dir / 'four_initial_states_summary.csv', index=False)
    if not report['all_four_pass']:
        raise SystemExit('one or more initial states failed; inspect four_initial_states.json')
    print(f'all four passed -> wrote {out_dir}')


if __name__ == '__main__':
    main()
