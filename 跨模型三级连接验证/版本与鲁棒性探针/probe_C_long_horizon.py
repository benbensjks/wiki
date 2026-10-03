"""Probe C: extend the frozen 34-state hybrid candidate from 300 h to 1200 h.

NEW FILE (added 2026-09-29).  It does NOT modify any existing script, model
file, frozen artefact or result folder.  Everything goes under this folder's
`results/` and (via the companion plotting script) `figures/`.

WHY THIS IS RUN, AND WHY IT IS CHEAP
------------------------------------
The 300 h candidate leaves a stated boundary: "only 300 h / 19 steady read
windows".  The Han upstream here is a TRUE integration, not end-point
interpolation, so extending the horizon only integrates more of the same
system.  One 1200 h run therefore removes that boundary and lines the model up
with the 51-state certification's 48 steady windows.

Because the output grid is equidistant (1 min), a longer run keeps the shorter
horizons inside it exactly: analysing `t[:n]` of the 1200 h run is the same
experiment as a standalone run of that length.  So ONE run yields three
horizons (300 / 600 / 1200 h), and the 300 h slice doubles as a bit-level
reproduction check against the frozen certification.

PRE-REGISTERED QUESTIONS
  Q1 does `certified_v1` still hold at 600 h and 1200 h?
  Q2 does the mod-8 sweep continue without drift (read-window spacing, and the
     decoded sequence over every steady window)?
  Q3 do the two event chains stay one-to-one / ordered / alternating as the run
     gets longer, or does a rare failure eventually appear?
  Q4 what does the extension actually cost (Han integration vs downstream
     integration, measured separately)?

INTEGRITY GUARDS
  G1 frozen working point asserted (HbyConfig defaults are the frozen ones).
  G2 the Han upstream covers the whole horizon and is NOT extrapolated.
  G3 the normalisation window is still the upstream module's own 1000-3000 min
     reference horizon: extending the run must not move `c31_max`.
  G4 the 0-300 h slice must reproduce the frozen certification exactly
     (19 steady reads, sequence, S0 hold margin, clock peak count) and must
     agree with the frozen saved trajectory.
  G5 every citation of "300 h / 19 windows" in the report is superseded by this
     run, not deleted: both numbers stay in the output.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')

WIKI = Path(r'C:\Users\18633\Desktop\wiki')
CROSS = WIKI / '跨模型三级连接验证'
HZH = CROSS / 'hby_zmh_hby'
HERE = Path(__file__).resolve().parent
OUT = HERE / 'results'

for _p in (CROSS, HZH, HZH / 'certification'):
    sys.path.insert(0, str(_p))

from hybrid_model import HbyConfig, sha256, source_hashes            # noqa: E402
from model_hzh import HZHModel, IDX, NAMES                           # noqa: E402
from run_han_comparison import HAN_PATH, HAN_SHA, HanInput           # noqa: E402
from verify_hzh import analyse, signals as raw_signals               # noqa: E402

HOURS = 1200.0
SAMPLE_MIN = 1.0
MAX_STEP_MIN = 1.0
RTOL, ATOL = 2e-7, 2e-9
SUB_HORIZONS = (300.0, 600.0, 1200.0)
ZOOM = (578.0, 620.0)
OV_STEP_MIN = 5.0
FROZEN_SUMMARY = HZH / 'certification' / 'results' / '20260928_141723_738780' / 'summary.json'
FROZEN_TRAJ = HZH / 'results' / '20260927_233646_924908' / 'han' / 'trajectory.npz'

CHECKS = []


def check(ok, msg):
    CHECKS.append((bool(ok), msg))
    print(('  OK   ' if ok else '  FAIL ') + msg, flush=True)
    return bool(ok)


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


def integrate(model, hours, sample_min, max_step_min, rtol, atol):
    count = round(hours * 60 / sample_min)
    t = np.linspace(0, hours, count + 1)
    if len(t) < 2 or not np.all(np.diff(t) > 0):
        raise ValueError('Bad sampling grid')
    sol = solve_ivp(model.rhs, (0, hours), model.initial_state(), t_eval=t,
                    method='DOP853', rtol=rtol, atol=atol, max_step=max_step_min / 60)
    if not sol.success or not np.all(np.isfinite(sol.y)):
        raise RuntimeError(sol.message)
    return sol.t, sol.y


def horizon_report(tt, yy, model, hours):
    verdict, sig = analyse(tt, yy, model.z, model.tail)
    steady, cold = verdict['steady'], verdict['cold']
    troughs = np.array([r['trough_h'] for r in verdict['read_windows']], dtype=float)
    gaps = np.diff(troughs)
    q = max(1, len(gaps) // 5)
    return dict(
        hours=hours,
        certified_v1=bool(verdict['certified_v1']),
        cold=dict(reads=cold['reads'], sequence=cold['sequence'], passed=cold['passed'],
                  minimum_commitment=cold['minimum_commitment'],
                  boundary_clips=cold['boundary_clips']),
        steady=dict(reads=steady['reads'], sequence=steady['sequence'],
                    passed=steady['passed'], complete=steady['complete'],
                    increments_mod8=steady['increments_mod8'],
                    minimum_commitment=steady['minimum_commitment'],
                    boundary_clips=steady['boundary_clips']),
        events={k: {j: x[j] for j in ('reverse_events', 'gate_events', 'flip_events',
                                      'one_to_one', 'causal_order',
                                      'alternating_directions', 'passed')}
                for k, x in verdict['events'].items()},
        bit_margins=verdict['bit_margins'],
        global_min_timing_margin_h=verdict['global_min_timing_margin_h'],
        clock_peak_count=int(verdict['clock_peak_count']),
        read_window_spacing=dict(n=len(gaps), median_h=float(np.median(gaps)),
                                 min_h=float(gaps.min()), max_h=float(gaps.max()),
                                 first_quintile_median_h=float(np.median(gaps[:q])),
                                 last_quintile_median_h=float(np.median(gaps[-q:])),
                                 drift_h=float(np.median(gaps[-q:]) - np.median(gaps[:q]))),
        signal_ranges={k: [float(np.min(v)), float(np.max(v))] for k, v in sig.items()},
        duty_g0_gt_0p05=float(np.mean(sig['g0'] > 0.05)),
        duty_g1_gt_0p05=float(np.mean(sig['g1'] > 0.05)),
        duty_clock_gt_0p05=float(np.mean(sig['clock'] > 0.05)),
        clock_gate_peak=float(np.max(sig['clock'])),
        verdict=jsonable(verdict))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    frozen = json.loads(FROZEN_SUMMARY.read_text(encoding='utf-8'))['run_results'][0]
    config = HbyConfig()
    meta = dict(
        script=Path(__file__).name, script_sha256=sha256(Path(__file__).resolve()),
        hours=HOURS, sample_min=SAMPLE_MIN, max_step_min=MAX_STEP_MIN,
        solver=dict(method='DOP853', rtol=RTOL, atol=ATOL),
        sub_horizons=list(SUB_HORIZONS), zoom_h=list(ZOOM),
        han_sha256=sha256(HAN_PATH), han_expected=HAN_SHA,
        receiver_config=jsonable(config.__dict__),
        source_hashes=source_hashes(),
        frozen_reference=str(FROZEN_SUMMARY.relative_to(WIKI)),
        purpose=('Remove the stated "only 300 h / 19 steady windows" boundary by '
                 'integrating the same frozen candidate for 1200 h; the Han upstream '
                 'is a true integration, not end-point interpolation.'))

    print('Integrating Han upstream to %.0f h...' % HOURS, flush=True)
    tic = time.perf_counter()
    han = HanInput(HOURS)
    han_seconds = time.perf_counter() - tic

    check(han.t[-1] >= 60 * HOURS - 1e-9,
          f'Han upstream covers the full horizon (t_end = {han.t[-1]:.1f} min '
          f'>= {60 * HOURS:.1f} min)')
    check(han.sol.y.shape[1] == han.t.size, 'Han solution grid is complete')
    nidx = np.flatnonzero(han.normalization_mask)
    check(nidx.size == 2001 and han.t[nidx[0]] == 1000.0 and han.t[nidx[-1]] == 3000.0,
          'normalisation mask is still exactly minutes 1000-3000 (2001 samples)')
    check(abs(han.c31_max - float(han.sol.y[7][nidx].max())) < 1e-12,
          f'c31_max {han.c31_max:.6f} is the max over that same reference window, '
          'so extending the run did not move the normalisation')
    last_mrna = han.interp(60 * HOURS, han.sol.y[6])
    check(np.isfinite(last_mrna), 'Han interpolation is defined at the last hour '
                                  '(no extrapolation)')

    print('Building model and integrating downstream to %.0f h...' % HOURS, flush=True)
    tic = time.perf_counter()
    model = HZHModel('han', han)
    build_seconds = time.perf_counter() - tic
    check((model.tail.c.clock_K_au, model.tail.c.clock_n) == (0.3, 2.0) and
          model.tail.c.n_A1_gate == 6.0 and model.tail.c.receiver_uM_per_au == 5.75,
          'frozen working point: clock 0.3/2.0, n_A1_gate 6, uM_per_au 5.75')
    tic = time.perf_counter()
    t, y = integrate(model, HOURS, SAMPLE_MIN, MAX_STEP_MIN, RTOL, ATOL)
    run_seconds = time.perf_counter() - tic
    check(len(t) == round(HOURS * 60 / SAMPLE_MIN) + 1,
          f'downstream trajectory has {len(t)} samples on a 1 min grid')

    rows = []
    for hours in SUB_HORIZONS:
        n = round(hours * 60 / SAMPLE_MIN) + 1
        print(f'Analysing horizon {hours:.0f} h ...', flush=True)
        rows.append(horizon_report(t[:n], y[:, :n], model, hours))
    by_hours = {r['hours']: r for r in rows}

    # ---- G4: the 300 h slice must reproduce the frozen certification exactly
    sub300 = by_hours[300.0]
    check(sub300['certified_v1'] == frozen['certified_v1'],
          f"300 h slice certified_v1 = {sub300['certified_v1']} "
          f"(frozen {frozen['certified_v1']})")
    check(sub300['steady']['reads'] == frozen['steady']['reads'],
          f"300 h slice steady reads = {sub300['steady']['reads']} "
          f"(frozen {frozen['steady']['reads']})")
    check(sub300['steady']['sequence'] == frozen['steady']['sequence'],
          '300 h slice steady sequence is byte-identical to the frozen one')
    check(sub300['cold']['sequence'] == frozen['cold']['sequence'],
          '300 h slice cold sequence is byte-identical to the frozen one')
    for bit in ('S0', 'S1', 'S2'):
        got = sub300['bit_margins'][bit]['min_hold_h']
        want = frozen['bit_margins'][bit]['min_hold_h']
        check(got is not None and abs(got - want) < 1e-9,
              f'300 h slice {bit} min hold margin {got!r} matches frozen {want!r}')
    check(sub300['clock_peak_count'] == frozen['clock_peak_count'],
          f"300 h slice clock peak count = {sub300['clock_peak_count']} "
          f"(frozen {frozen['clock_peak_count']})")
    for chain in ('bit0_to_bit1', 'bit1_to_bit2'):
        got = sub300['events'][chain]
        want = frozen['event_chains'][chain]
        check(all(got[k] == want[k] for k in ('reverse_events', 'gate_events',
                                              'flip_events', 'one_to_one',
                                              'causal_order', 'alternating_directions')),
              f'300 h slice chain {chain} matches frozen '
              f'({want["reverse_events"]}/{want["gate_events"]}/{want["flip_events"]})')

    # ---- G4b: and with the frozen saved trajectory itself
    traj_gap = None
    if FROZEN_TRAJ.exists():
        z = np.load(FROZEN_TRAJ)
        tf, yf = z['time_h'], z['states']
        if tuple(z['state_names']) == tuple(NAMES) and np.allclose(tf, t[:len(tf)], rtol=0, atol=1e-10):
            traj_gap = float(np.max(np.abs(y[:, :len(tf)] - yf)))
            check(traj_gap < 1e-4,
                  f'0-300 h states agree with the frozen saved trajectory '
                  f'(max abs gap {traj_gap:.3e})')
        else:
            check(False, 'frozen saved trajectory has a different layout or grid')
    else:
        print('  note  frozen saved trajectory not found; state-level cross-check skipped')

    # ---- Q2: drift over the longest horizon
    long = by_hours[1200.0]
    sp = long['read_window_spacing']
    check(abs(sp['drift_h']) < 0.05,
          f"read-window spacing drift (last vs first quintile) = {sp['drift_h']:+.4f} h")
    check(long['steady']['increments_mod8'] and long['steady']['boundary_clips'] == 0,
          f"1200 h: mod-8 increments hold over all {long['steady']['reads']} steady "
          f"windows with 0 boundary clips")
    check(min(r['certified_v1'] for r in rows),
          'certified_v1 holds at every probed horizon (300 / 600 / 1200 h)')

    # ---- series for the figures (downsampled overview + full-resolution zoom)
    sig = raw_signals(y, model.z, model.tail)
    ov = slice(None, None, int(round(OV_STEP_MIN / SAMPLE_MIN)))
    zm = (t >= ZOOM[0]) & (t <= ZOOM[1])
    series = dict(
        hours=HOURS, zoom=np.array(ZOOM), ov_step_min=OV_STEP_MIN,
        ov_time=t[ov], ov_S0=sig['S0'][ov], ov_S1=sig['S1'][ov], ov_S2=sig['S2'][ov],
        ov_g0=sig['g0'][ov], ov_g1=sig['g1'][ov], ov_clock=sig['clock'][ov],
        ov_Jrev0=sig['J_rev0'][ov], ov_Jrev1=sig['J_rev1'][ov],
        ov_Jfwd2=sig['J_fwd2'][ov], ov_Jrev2=sig['J_rev2'][ov],
        ov_Int0=y[IDX['b0_I']][ov], ov_Int1=y[IDX['I1_zmh']][ov],
        ov_Int2=y[IDX['b2_I']][ov],
        ov_m31=np.interp(60 * t[ov], han.t, han.sol.y[6]),
        ov_flux=np.interp(60 * t[ov], han.t, han.flux_copies_min),
        zw_time=t[zm], zw_y=y[:, zm].astype(np.float32),
        zw_g0=sig['g0'][zm], zw_g1=sig['g1'][zm], zw_clock=sig['clock'][zm],
        zw_Jrev0=sig['J_rev0'][zm], zw_Jrev1=sig['J_rev1'][zm],
        zw_Jfwd2=sig['J_fwd2'][zm], zw_Jrev2=sig['J_rev2'][zm],
        state_names=np.array(NAMES),
        read_trough=np.array([r['trough_h'] for r in long['verdict']['read_windows']]),
        read_start=np.array([r['start_h'] for r in long['verdict']['read_windows']]),
        read_end=np.array([r['end_h'] for r in long['verdict']['read_windows']]),
        read_value=np.array([np.nan if r['value'] is None else r['value']
                             for r in long['verdict']['read_windows']]),
    )
    np.savez_compressed(OUT / 'probe_C_long_horizon_series.npz', **series)

    result = dict(
        meta=meta,
        timing=dict(han_integration_s=han_seconds, model_build_s=build_seconds,
                    downstream_integration_s=run_seconds,
                    total_s=han_seconds + build_seconds + run_seconds,
                    frozen_300h_downstream_s=15.788040000014007,
                    note=('the frozen 300 h figure is the downstream integration only, '
                          'measured in results/20260927_233646_924908/han/summary.json')),
        horizons=rows,
        frozen_reference=jsonable(frozen),
        cross_check_against_frozen_trajectory_max_abs_gap=traj_gap,
        closing_statement=(
            'The "only 300 h / 19 steady windows" boundary is removed: the same frozen '
            'candidate is reported at 300 / 600 / 1200 h from a single run, and the '
            '300 h slice reproduces the frozen certification exactly.'))
    (OUT / 'probe_C_long_horizon.json').write_text(
        json.dumps(jsonable(result), ensure_ascii=False, indent=1), encoding='utf-8')

    hdr = ['hours', 'certified_v1', 'cold_reads', 'steady_reads', 'steady_sequence',
           'min_commitment', 'boundary_clips', 'chain01_rev', 'chain01_pass',
           'chain12_rev', 'chain12_pass', 'min_timing_margin_h', 'clock_peaks',
           'median_window_h', 'spacing_drift_h']
    lines = [','.join(hdr)]
    for r in rows:
        lines.append(','.join(str(x) for x in [
            r['hours'], r['certified_v1'], r['cold']['reads'], r['steady']['reads'],
            r['steady']['sequence'], r['steady']['minimum_commitment'],
            r['steady']['boundary_clips'],
            r['events']['bit0_to_bit1']['reverse_events'],
            r['events']['bit0_to_bit1']['passed'],
            r['events']['bit1_to_bit2']['reverse_events'],
            r['events']['bit1_to_bit2']['passed'],
            r['global_min_timing_margin_h'], r['clock_peak_count'],
            r['read_window_spacing']['median_h'], r['read_window_spacing']['drift_h']]))
    (OUT / 'probe_C_long_horizon.csv').write_text('\n'.join(lines) + '\n', encoding='utf-8')
    print('\n'.join(lines))

    bad = [m for ok, m in CHECKS if not ok]
    print(f'\n{len(CHECKS) - len(bad)}/{len(CHECKS)} checks passed')
    print('written to', OUT)
    if bad:
        raise SystemExit('FAILED CHECKS:\n' + '\n'.join(bad))


if __name__ == '__main__':
    main()
