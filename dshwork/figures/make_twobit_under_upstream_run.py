r"""One 300 h run of the two-bit counter under Han's oscillator, at the SELECTED profile.

Why this file exists
--------------------
No figure in the deliverable set shows Han's oscillator actually oscillating together with
the two-bit counter completing whole mod-4 cycles.  Figure 4-2 (the two-bit readout) plots
only `S0` / `S1` and the read strip -- it contains no upstream trace at all -- and
`twobit34_results/figures/01_counter_overview` shows the upstream only as a C31 flux line,
with the oscillator's own state variables absent.

It also fixes a working-point problem while it is at it.  `run_twobit34.py` builds the
model as a bare `TwoBit34Model()`, which silently lands on `nominal_extension()`:

    nominal_extension()          uM_per_au = 6.0   carry A/F maturation = 30.0 min
    selected_profile.json        uM_per_au = 5.75  carry A/F maturation = 32.5 min

so the main 34-state deliverable (its `parameters.json` records 6.0 / 30.0) runs at a
DIFFERENT point from the one the outline labels it with ("5.75 / 32.5 / 2") and from the
one the 51-state chain, the 140-point perturbation scan and `initial_states_selected` all
use.  This script builds the model EXPLICITLY, ASSERTS every field against
`selected_profile.json`, and ASSERTS that it is not on the bare constructor default --
the guard whose absence let that split happen.

Run:  & 'D:\aconade\python.exe' .\make_twobit_under_upstream_run.py
Out:  ./data/twobit_selected_300h.csv   time_h + 34 named states + fluxes/gate signals
      ./data/twobit_selected_300h.json  parameters, hashes, verdict, read windows,
                                        and a cross-check against the 75-point grid row
"""
from __future__ import annotations

import hashlib
import json
import sys
import time
from dataclasses import asdict, replace
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
FR = HERE.parents[1] / 'final_reconstruction'
DATA = HERE / 'data'
STEM = 'twobit_selected_300h'

HOURS = 300.0
SAMPLE_MIN = 2.0
MAX_STEP_MIN = 2.0

# The selected robust centre, verbatim from selected_profile.json.  Written out here so
# that the assertion below has something independent to check the file against.
SELECTED = dict(uM_per_au=5.75, carry_mrna_half_life_min=2.0,
                carry_maturation_half_life_min=32.5)

DOCUMENTED_INDEX = {'b0_S': 16, 'b1_S': 27, 'A0': 28, 'F0': 29}
OSCILLATOR_NAMES = ('m_TetR', 'TetR_total', 'm_CI', 'CI', 'm_LacI', 'LacI')
PROFILE = FR / 'twobit34_results' / 'selected_profile.json'
GRID = FR / 'twobit34_results' / 'robustness' / 'robustness_all.csv'


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for blk in iter(lambda: fh.read(1 << 20), b''):
            h.update(blk)
    return h.hexdigest().upper()


def build_model():
    """Explicitly construct the selected profile and prove it is the selected profile."""
    sys.path.insert(0, str(FR))
    from model_twobit34 import (CarryExpressionParameters, TwoBit34Model,
                                nominal_extension)

    profile = json.loads(PROFILE.read_text(encoding='utf-8'))
    bare_e = nominal_extension()
    bare_c = CarryExpressionParameters()

    ext = replace(bare_e, uM_per_au=float(profile['extension']['uM_per_au']))
    carry = CarryExpressionParameters(
        mrna_half_life_min=float(profile['carry_expression']['mrna_half_life_min']),
        activator_maturation_half_life_min=float(
            profile['carry_expression']['activator_maturation_half_life_min']),
        repressor_maturation_half_life_min=float(
            profile['carry_expression']['repressor_maturation_half_life_min']))
    model = TwoBit34Model(ext, carry)

    # (a) every field of the constructed model equals selected_profile.json, field by field
    problems = []
    for key, want in profile['extension'].items():
        got = getattr(model.e, key)
        if got != want:
            problems.append(f'extension.{key}: model has {got!r}, profile says {want!r}')
    for key, want in profile['carry_expression'].items():
        got = getattr(model.carry, key)
        if got != want:
            problems.append(f'carry_expression.{key}: model has {got!r}, profile says '
                            f'{want!r}')
    if problems:
        raise SystemExit('FATAL: the constructed model is not the selected profile:\n  '
                         + '\n  '.join(problems))

    # (b) the point is NOT the bare constructor default -- this is the guard whose absence
    #     let the main deliverable be labelled 5.75/32.5 while running at 6.0/30.0
    if model.e.uM_per_au == bare_e.uM_per_au:
        raise SystemExit(f'FATAL: uM_per_au is still the bare nominal_extension default '
                         f'{bare_e.uM_per_au}; this run is not on the selected profile')
    if (model.carry.activator_maturation_half_life_min
            == bare_c.activator_maturation_half_life_min):
        raise SystemExit('FATAL: carry maturation is still the bare '
                         'CarryExpressionParameters default; not the selected profile')

    # (c) the states we are about to slice are where the project says they are
    index = {n: i for i, n in enumerate(model.state_names)}
    for name, want in DOCUMENTED_INDEX.items():
        if index.get(name) != want:
            raise SystemExit(f'FATAL: state_names[{name!r}] = {index.get(name)}, '
                             f'documented layout says {want}')

    print(f'shared selected profile asserted: uM_per_au = {model.e.uM_per_au} '
          f'(bare default {bare_e.uM_per_au}), carry mRNA = '
          f'{model.carry.mrna_half_life_min} min, carry A/F maturation = '
          f'{model.carry.activator_maturation_half_life_min} min '
          f'(bare default {bare_c.activator_maturation_half_life_min} min), '
          f'clock = {model.e.clock_K_au}/{model.e.clock_n}, '
          f'add_growth = {model.e.add_growth}')
    return model, index, profile


def main() -> int:
    try:
        sys.stdout.reconfigure(encoding='utf-8', errors='replace')
    except Exception:                                                   # noqa: BLE001
        pass
    started = time.perf_counter()
    model, index, profile = build_model()
    from verify_twobit_causal import analyse_solution
    from model_twobit34 import STATE_NAMES_34

    sol = model.simulate(hours=HOURS, sample_min=SAMPLE_MIN,
                         max_step_min=MAX_STEP_MIN)
    t, y = sol.t, sol.y
    flux = np.asarray([model.flux(y[:, k]) for k in range(y.shape[1])])
    sig = model.diagnostic_signals(y)
    analysis = analyse_solution(model.base, sol, asdict(model.e), HOURS)

    # The signal dict repeats two states under unprefixed names; prove that identity
    # rather than shipping a second, ambiguously named copy of the same trace.
    state_of_signal = {'S0': 'b0_S', 'S1': 'b1_S'}
    dup = {s: float(np.abs(np.asarray(sig[s]) - y[index[st]]).max())
           for s, st in state_of_signal.items()}
    if max(dup.values()) != 0.0:
        raise SystemExit(f'FATAL: diagnostic_signals S0/S1 are not identical to '
                         f'b0_S/b1_S: {dup}')
    extra_signals = ('J_fwd0', 'J_rev0', 'g0', 'Int1_source')

    import pandas as pd
    rows = dict(time_h=t, C31_flux_au_h=flux)
    for name in STATE_NAMES_34:
        rows[name] = y[index[name]]
    for name in extra_signals:
        rows[name] = sig[name]
    DATA.mkdir(parents=True, exist_ok=True)
    csv = DATA / f'{STEM}.csv'
    pd.DataFrame(rows).to_csv(csv, index=False, encoding='utf-8')

    windows = [dict(cycle=int(w['cycle']), trough_h=float(w['trough_h']),
                    window_start_h=float(w['window_start_h']),
                    window_end_h=float(w['window_end_h']), value=w['value'],
                    bit0=w['bit0']['label'], bit1=w['bit1']['label'],
                    bit0_commitment=float(w['bit0']['commitment']),
                    bit1_commitment=float(w['bit1']['commitment']))
               for w in analysis['read_windows']]

    # ------------------------------------------------ cross-check the 75-point grid row
    grid_row = None
    if GRID.is_file():
        g = pd.read_csv(GRID)
        m = ((g.uM_per_au.round(6) == SELECTED['uM_per_au'])
             & (g.carry_maturation_half_life_min == SELECTED['carry_maturation_half_life_min'])
             & (g.carry_mrna_half_life_min == SELECTED['carry_mrna_half_life_min']))
        if int(m.sum()) == 1:
            r = g[m].iloc[0]
            grid_row = dict(
                row_found=True,
                uM_per_au=float(r.uM_per_au),
                carry_maturation_half_life_min=float(r.carry_maturation_half_life_min),
                carry_mrna_half_life_min=float(r.carry_mrna_half_life_min),
                certified=bool(r.certified), reads=int(r.reads),
                cold_sequence=str(r.cold_sequence),
                bit1_setup_h=float(r.bit1_setup_h), bit1_hold_h=float(r.bit1_hold_h),
                this_run_setup_h=float(analysis['margins']['bit1']['min_setup_h']),
                this_run_hold_h=float(analysis['margins']['bit1']['min_hold_h']),
                this_run_reads=int(analysis['cold_start']['reads']),
                this_run_sequence=str(analysis['cold_start']['sequence']),
            )
            grid_row['matches'] = bool(
                grid_row['certified'] == bool(analysis['certified'])
                and grid_row['reads'] == grid_row['this_run_reads']
                and grid_row['cold_sequence'] == grid_row['this_run_sequence']
                and abs(grid_row['bit1_setup_h'] - grid_row['this_run_setup_h']) < 1e-9
                and abs(grid_row['bit1_hold_h'] - grid_row['this_run_hold_h']) < 1e-9)

    doc = dict(
        purpose=('one 300 h run of the 34-state two-bit counter driven by the real Han '
                 'oscillator, at the selected profile, for the figure that shows the '
                 'oscillator and the counter together'),
        status='NEW DETERMINISTIC RUN made for this figure - not a pre-existing artefact',
        working_point=dict(
            uM_per_au=float(model.e.uM_per_au),
            carry_mrna_half_life_min=float(model.carry.mrna_half_life_min),
            carry_maturation_half_life_min=float(
                model.carry.activator_maturation_half_life_min),
            clock_K_au=float(model.e.clock_K_au), clock_n=float(model.e.clock_n),
            add_growth=bool(model.e.add_growth),
            source='twobit34_results/selected_profile.json',
            asserted=('every extension and carry_expression field is compared to that '
                      'profile, and the run is asserted NOT to sit on the bare '
                      'nominal_extension() / CarryExpressionParameters() defaults'),
            bare_defaults=dict(uM_per_au=6.0, carry_maturation_half_life_min=30.0,
                               note=('run_twobit34.py uses the bare default; the split '
                                     'between this profile and that default is why this '
                                     'script asserts which one it is on'))),
        extension={k: v for k, v in vars(model.e).items()},
        carry_expression={k: v for k, v in vars(model.carry).items()},
        hours=HOURS, sample_min=SAMPLE_MIN, max_step_min=MAX_STEP_MIN,
        n_samples=int(t.size),
        state_index={n: index[n] for n in STATE_NAMES_34},
        signal_columns=dict(
            kept=list(extra_signals),
            dropped_reason=('model_twobit34.diagnostic_signals also returns S0 and S1, '
                            'which are bit-identical to b0_S and b1_S; the identity is '
                            'asserted above and the duplicate columns are not shipped, '
                            'so no bare S0/S1 label can be mistaken for a state')),
        certified=bool(analysis['certified']),
        causal_verdict=analysis['causal_verdict'],
        cold_start=dict(reads=analysis['cold_start']['reads'],
                        sequence=analysis['cold_start']['sequence'],
                        passed=bool(analysis['cold_start']['passed']),
                        minimum_commitment=float(
                            analysis['cold_start']['minimum_commitment'])),
        steady_state=dict(reads=analysis['steady_state']['reads'],
                          sequence=analysis['steady_state']['sequence'],
                          passed=bool(analysis['steady_state']['passed']),
                          drop_reads=analysis['steady_state']['drop_reads']),
        margins=analysis['margins'],
        reverse_events=len(analysis['reverse_events']),
        gate_events=len(analysis['gate_events']),
        signal_ranges=analysis['signal_ranges'],
        read_windows=windows,
        grid_row_cross_check=grid_row,
        oscillator_states=list(OSCILLATOR_NAMES),
        model_sha256=sha256(FR / 'model_twobit34.py'),
        source_sha256={n: sha256(FR / n) for n in
                       ('model.py', 'model_twobit34.py', 'verify_twobit_causal.py')},
        profile_sha256=sha256(PROFILE),
        csv_sha256=sha256(csv),
        runtime_s=round(time.perf_counter() - started, 1),
    )
    (DATA / f'{STEM}.json').write_text(
        json.dumps(doc, ensure_ascii=False, indent=2), encoding='utf-8')

    print(f"certified={doc['certified']} reads={doc['cold_start']['reads']} "
          f"seq={doc['cold_start']['sequence']} "
          f"setup={analysis['margins']['bit1']['min_setup_h']:.6f} "
          f"hold={analysis['margins']['bit1']['min_hold_h']:.6f}")
    if grid_row is not None:
        print(f"grid row cross-check: matches={grid_row['matches']} "
              f"(grid setup {grid_row['bit1_setup_h']:.6f}, hold "
              f"{grid_row['bit1_hold_h']:.6f})")
        if not grid_row['matches']:
            raise SystemExit('FATAL: this run does not reproduce the 75-point grid row '
                             'for the selected profile; do not plot it')
    print(f'wrote {csv}')
    print(f'wrote {DATA / (STEM + ".json")}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
