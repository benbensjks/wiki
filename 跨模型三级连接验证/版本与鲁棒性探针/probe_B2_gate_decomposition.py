"""Probe B2: decompose the carry-1 gate collapse into act(A0) and rep(F0).

NEW FILE (added 2026-09-29).  Second-round probe, written after the first-round
report was reviewed.  It does NOT modify any existing script, model file, frozen
artefact or result: it only reads, re-runs the SAME three working points as
probe B, and writes new outputs under this folder.

WHY IT EXISTS
-------------
Probe B showed the version swap kills the counter, and located the break at the
carry-1 gate: g0 peak 0.442834 -> 0.0229591 -> 0.000222469 while bit0's own
reverse flux J_rev0 only falls 0.931613 -> 0.268935 -> 0.199769.  Those two
amplitudes do not match, so "weaker reset -> weaker pulse" is NOT a sufficient
mechanism.  The frozen diagnostic set does not expose A0 or F0, so the first
report could not say WHICH factor of

    g0 = act(A0 ; K_A[0], n_A[0]) * rep(F0 ; K_F[0], n_F[0])

collapsed.  This probe measures it directly.

The mapping is imported from `probe_B_zeng_version_51state.build_patch`, so the
two probes cannot drift apart; the working point, solver settings and the
assertions are the same as probe B.

PRE-REGISTERED QUESTIONS (fixed before the runs)
  Q1 does the gate fail through the activator arm act(A0) or the repressor arm
     rep(F0)?  Reported as the value of each factor AT the g0 peak, for every
     case, plus each factor's own extremes over the run.
  Q2 is the gate peak still phase-locked to bit0's PB->LR transition (i.e. to the
     reset event), or has the phase relationship moved?
  Q3 does the drive itself change?  clock_gate and the upstream period are the
     same in all three cases by construction, so this is a control, not a result.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')

WIKI = Path(r'C:\Users\18633\Desktop\wiki')
FR = WIKI / 'final_reconstruction'
HERE = Path(__file__).resolve().parent
OUT = HERE / 'results'
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(FR))
sys.path.insert(0, str(FR / 'plausibility'))

import model as M                                     # noqa: E402
import plausibility_common as PC                      # noqa: E402
import probe_B_zeng_version_51state as PB             # noqa: E402

from model import act, rep                            # noqa: E402

HOURS = 600.0
SAMPLE_MIN = 2.0
EXPECTED_CLOCK = PB.EXPECTED_CLOCK
EXPECTED_N_GATE = PB.EXPECTED_N_GATE
EXPECTED_UM_PER_AU = PB.EXPECTED_UM_PER_AU
CROSS_FRACTION = 0.5


def jsonable(o):
    if isinstance(o, dict):
        return {str(k): jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [jsonable(v) for v in o]
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (bool, int, float, str)) or o is None:
        return o
    return str(o)


def run_case(tag, patch):
    original = {k: M.ZENG[k] for k in patch}
    row = {'tag': tag, 'patch': jsonable(patch)}
    started = time.perf_counter()
    try:
        for key, value in patch.items():
            M.ZENG[key] = value
        model, residual = PC.build_threebit(n_A1_gate=EXPECTED_N_GATE)
        if residual:
            raise RuntimeError(f'unexpected residual patch {residual}')
        if (model.e.clock_K_au, model.e.clock_n) != EXPECTED_CLOCK:
            raise RuntimeError('clock gate is not the frozen working point')
        if float(model.n_A1_gate_effective) != EXPECTED_N_GATE:
            raise RuntimeError('gate exponent is not the frozen working point')
        if float(model.e.uM_per_au) != EXPECTED_UM_PER_AU:
            raise RuntimeError('uM_per_au is not the frozen working point')
        sol = model.simulate(hours=HOURS, sample_min=SAMPLE_MIN,
                             rtol=2e-7, atol=2e-9, max_step_min=2.0)
        t, y = sol.t, sol.y
    finally:
        PC.restore(original)
    row['runtime_s'] = time.perf_counter() - started

    S0 = y[16]
    PB0 = 1.0 - S0
    A0, F0 = y[28], y[29]
    R0, T0, C0 = y[14], y[11], y[15]
    I0 = y[8]
    act0 = np.array([act(a, M.ZENG['K_A'][0], M.ZENG['n_A'][0]) for a in A0])
    rep0 = np.array([rep(f, M.ZENG['K_F'][0], M.ZENG['n_F'][0]) for f in F0])
    g0 = act0 * rep0
    inhibit = M.ZENG['K_inh'] / (M.ZENG['K_inh'] + np.maximum(R0, 0.0))

    peak = int(np.argmax(g0))
    # PB0 rising edges: the reset events that are supposed to launch a carry pulse
    z = PB0 - CROSS_FRACTION
    edges = np.flatnonzero(z[:-1] * z[1:] < 0)
    rising = [int(i) for i in edges if z[i + 1] > 0]
    gaps = np.diff(t[rising]) if len(rising) > 1 else np.array([np.nan])

    # per-cycle gate peak and the factors at that peak
    peaks = []
    for a, b in zip(rising[:-1], rising[1:]):
        j = a + int(np.argmax(g0[a:b + 1]))
        peaks.append(dict(cycle_start_h=float(t[a]), peak_h=float(t[j]),
                          g0=float(g0[j]), act_of_A0=float(act0[j]), rep_of_F0=float(rep0[j]),
                          A0=float(A0[j]), F0=float(F0[j]), S0=float(S0[j]),
                          delay_after_reset_h=float(t[j] - t[a])))

    row.update(
        gate_peak=dict(g0=float(g0.max()), act_of_A0=float(act0[peak]),
                       rep_of_F0=float(rep0[peak]), A0=float(A0[peak]), F0=float(F0[peak]),
                       S0=float(S0[peak]), PB0=float(PB0[peak]), t_h=float(t[peak])),
        factor_extremes=dict(
            act_of_A0=[float(act0.min()), float(act0.max())],
            rep_of_F0=[float(rep0.min()), float(rep0.max())],
            A0=[float(A0.min()), float(A0.max())],
            F0=[float(F0.min()), float(F0.max())],
            PB0=[float(PB0.min()), float(PB0.max())],
            S0=[float(S0.min()), float(S0.max())],
            R0=[float(R0.min()), float(R0.max())],
            T0=[float(T0.min()), float(T0.max())],
            C0=[float(C0.min()), float(C0.max())],
            forward_inhibition_Kinh_over_Kinh_plus_R0=[float(inhibit.min()), float(inhibit.max())]),
        reset_events=dict(count=len(rising),
                          median_period_h=float(np.nanmedian(gaps)) if len(rising) > 1 else None,
                          pb0_duty_above_half=float(np.mean(PB0 >= CROSS_FRACTION))),
        per_cycle_gate_peaks=peaks,
    )
    return row


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cases = [('baseline_frozen_zeng', {})] + [(k, PB.build_patch(k)) for k in PB.VARIANT_BLOCKS]
    rows = []
    for tag, patch in cases:
        print(f'[run] {tag} ...', flush=True)
        row = run_case(tag, patch)
        rows.append(row)
        gp = row['gate_peak']
        print(f"      g0 peak {gp['g0']:.6g} at t={gp['t_h']:.2f} h | "
              f"act(A0) {gp['act_of_A0']:.6g} | rep(F0) {gp['rep_of_F0']:.6g} | "
              f"reset events {row['reset_events']['count']} "
              f"({row['runtime_s']:.1f} s)", flush=True)

    base = rows[0]
    result = dict(
        meta=dict(script=Path(__file__).name,
                  script_sha256=PB.sha256(Path(__file__).resolve()),
                  mapping_source='probe_B_zeng_version_51state.build_patch (imported)',
                  hours=HOURS, sample_min=SAMPLE_MIN,
                  working_point=dict(clock_K_au=EXPECTED_CLOCK[0], clock_n=EXPECTED_CLOCK[1],
                                     n_A1_gate=EXPECTED_N_GATE, uM_per_au=EXPECTED_UM_PER_AU),
                  source_sha256=PC.source_hashes()),
        rows=rows,
        reading=dict(
            q1=('Which arm of g0 = act(A0)*rep(F0) collapses: compare gate_peak.act_of_A0 '
                'and gate_peak.rep_of_F0 against the baseline row.'),
            q2=('Every per_cycle_gate_peaks entry carries delay_after_reset_h, so the '
                'phase-lock to the PB0 rising edge is visible per cycle.'),
            q3=('clock_gate is not recomputed here; probe B already reported it as '
                'unchanged (0.769725 / 0.769549 / 0.769397).')))
    (OUT / 'probe_B2_gate_decomposition.json').write_text(
        json.dumps(jsonable(result), ensure_ascii=False, indent=1), encoding='utf-8')

    hdr = ('tag,g0_peak,act_of_A0_at_peak,rep_of_F0_at_peak,A0_at_peak,F0_at_peak,'
           'A0_max,F0_max,act_min,rep_min,reset_events,pb0_duty_above_half,'
           'gate_peaks_per_cycle_median,runtime_s')
    lines = [hdr]
    for r in rows:
        gp, fx, rs = r['gate_peak'], r['factor_extremes'], r['reset_events']
        gpv = [p['g0'] for p in r['per_cycle_gate_peaks']]
        lines.append(','.join(str(x) for x in [
            r['tag'], gp['g0'], gp['act_of_A0'], gp['rep_of_F0'], gp['A0'], gp['F0'],
            fx['A0'][1], fx['F0'][1], fx['act_of_A0'][0], fx['rep_of_F0'][0],
            rs['count'], rs['pb0_duty_above_half'],
            float(np.median(gpv)) if gpv else '', r['runtime_s']]))
    (OUT / 'probe_B2_gate_decomposition.csv').write_text('\n'.join(lines) + '\n',
                                                        encoding='utf-8')
    print('\n'.join(lines))
    print('written to', OUT)


if __name__ == '__main__':
    main()
