r"""Matched n=4 trajectory for Figure 4-3 Panel A   (data generator, not a figure)

Why this file exists
--------------------
Figure 4-3 Panel A compares the S2 trajectory of the SAME expression settings at
two gate exponents.  The frozen run (n_A1_gate = 6) already exists as
`threebit51_results/wiki_n6_20260924/`.  No n=4 run with the identical expression
settings existed -- the legacy n=4 scans predate `n_A1_gate` and differ in their
other settings -- so this script produces exactly that one missing run.

It is a NEW deterministic run made for the figure.  It is not a pre-existing
artefact and must not be cited as one.  Everything that identifies it (model hash,
the realised gate exponent, the expression settings, hours, grid) is written next to
the data in `data/n4_matched_600h.json`.

The build goes through `plausibility_common.build_threebit(n_A1_gate=4.0)`, which
ASSERTS the exponent the model actually reports (H1): leaving the constructor
argument out would silently fall back to the ZENG table value and there would be no
way to tell the two apart from the output.

Column -> state mapping
-----------------------
State indices are looked up BY NAME in `model.state_names`; they are never written as
bare integers here.  The lookup is then asserted against the documented 51-state
layout, because `state_names` itself is what is being trusted:

    b0_S = 16   b1_S = 27   b2_S = 44   A1 = 45   F1 = 46

An earlier version of this script hardcoded ``S0=sol.y[44], S1=sol.y[45],
S2=sol.y[46]``, which wrote **b2_S, A1 and F1** under the names S0, S1, S2.  Figure 4-3
then plotted bit 2's F1 protein as the n=4 "S2" curve (F1 reaches 2.39 a.u., so it is
not even a fraction), and the caption's "S2 peaks at 2.25, outside the physical range"
was describing F1.  The script now also re-checks its own output before returning:
fractions must lie in [0, 1], the saved b2_S column must REPLAY the recorded read-window
bit-2 labels, and the saved A1/F1 must reproduce the gate identity

    g1 = act(A1, K_A[1], n_gate) * rep(F1, K_F[1], n_F[1]) * clock_gate

which ties the columns back to the physics rather than to their own header.

Run:  & 'D:\aconade\python.exe' .\make_n4_matched_trajectory.py
Out:  ./data/n4_matched_600h.csv    time_h, b0_S, b1_S, b2_S, A1, F1, g1,
                                    clock_gate, C31_flux
      ./data/n4_matched_600h.json   provenance + verdict + read windows + index map

Cross-checker: `check_state_column_semantics.py` re-verifies the same semantics from
disk alone, for this file and for the wiki n=6 artefact.
"""
from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np                                                       # noqa: E402

HERE = Path(__file__).resolve().parent
DATA = HERE / 'data'
FR = HERE.parents[1] / 'final_reconstruction'
STEM = 'n4_matched_600h'
HOURS = 600.0
SAMPLE_MIN = 2.0
MAX_STEP_MIN = 2.0
GATE = 4.0

# The documented 51-state layout, asserted AGAINST model.state_names below.  These are
# not the values used to slice the state vector -- the name lookup is -- they are the
# independent statement of where the states are supposed to live.
DOCUMENTED_INDEX = {'b0_S': 16, 'b1_S': 27, 'b2_S': 44, 'A1': 45, 'F1': 46}
SAVED_STATES = ('b0_S', 'b1_S', 'b2_S', 'A1', 'F1', 'b2_I', 'b2_R', 'b2_C')
# bit 2's pools and its two recombination fluxes: the whole-run evidence behind the
# figure's claim that the reverse write is insufficient.  The wiki n=6 artefact already
# carries all of these, so saving them here also makes the n=4 run usable for the
# low-copy-phase write-flux analysis without a second integration.
SAVED_FLUXES = ('J_fwd2', 'J_rev2')
BAND_LOW, BAND_HIGH, MIN_OCCUPANCY = 0.30, 0.70, 0.80


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest().upper()


def main() -> int:
    try:
        sys.stdout.reconfigure(encoding='utf-8', errors='replace')
    except Exception:                                                # noqa: BLE001
        pass
    sys.path.insert(0, str(FR))
    sys.path.insert(0, str(FR / 'plausibility'))
    from plausibility_common import build_threebit, frozen_extension
    from verify_threebit51 import analyse_threebit

    started = time.perf_counter()
    model, patched = build_threebit(n_A1_gate=GATE)
    realised = float(model.n_A1_gate_effective)
    if realised != GATE:
        raise SystemExit(f'FATAL: requested n_A1_gate={GATE} but the model reports '
                         f'{realised}')
    print(f'n_A1_gate realised = {realised} (asserted by the shared builder)')
    sol = model.simulate(hours=HOURS, sample_min=SAMPLE_MIN,
                         max_step_min=MAX_STEP_MIN)
    a = analyse_threebit(model, sol, HOURS)
    sig = model.diagnostic_signals(sol.y)

    # ------------------------------------------------ column -> state, BY NAME
    # Never slice with a bare integer: that is how this file previously wrote b2_S, A1
    # and F1 under the names S0, S1 and S2.
    index = {name: i for i, name in enumerate(model.state_names)}
    for name, want in DOCUMENTED_INDEX.items():
        got = index.get(name)
        if got != want:
            raise SystemExit(f'FATAL: model.state_names[{name!r}] = {got}, but the '
                             f'documented 51-state layout says {want}. Refusing to '
                             f'write a trajectory whose columns cannot be trusted.')
    print(f'state index map (from model.state_names, asserted): '
          f'{ {n: index[n] for n in SAVED_STATES} }')

    import pandas as pd
    # `diagnostic_signals` has no C31-flux key; the upstream flux is model.flux(y).
    # Record it because the read-window geometry follows its peaks.
    c31 = np.array([model.flux(sol.y[:, k]) for k in range(sol.y.shape[1])])
    rows = dict(time_h=sol.t, g1=sig['g1'], clock_gate=sig['clock_gate'], C31_flux=c31)
    for name in SAVED_STATES:
        rows[name] = sol.y[index[name]]
    for name in SAVED_FLUXES:
        rows[name] = sig[name]
    DATA.mkdir(parents=True, exist_ok=True)
    csv = DATA / f'{STEM}.csv'
    pd.DataFrame(rows).to_csv(csv, index=False, encoding='utf-8')

    # ------------------------------------------- whole-run excursion evidence
    # The figure's caption claims that bit 2 enters the high band and never returns.
    # Measure that on the saved column, over the WHOLE run, so the claim is machine-
    # checked rather than asserted.  `J_fwd2` / `J_rev2` are the only two channels that
    # can move S2, so they are what "the reverse write is insufficient" has to rest on.
    import numpy as _np
    tt = _np.asarray(sol.t, dtype=float)
    S2 = _np.asarray(rows['b2_S'], dtype=float)
    hi = _np.flatnonzero(S2 >= BAND_HIGH)
    t_first_high = float(tt[hi[0]]) if hi.size else None
    after = S2[hi[0]:] if hi.size else S2
    after_t = tt[hi[0]:] if hi.size else tt
    returns_low = bool((after <= BAND_LOW).any())
    excursion = dict(
        band_low=BAND_LOW, band_high=BAND_HIGH,
        t_first_high_h=t_first_high,
        s2_min_after_first_high=float(after.min()) if hi.size else None,
        s2_max_after_first_high=float(after.max()) if hi.size else None,
        returns_to_low_band_after_first_high=returns_low,
        n_samples_after_first_high=int(after.size),
        # trapezoidal integrals over the whole run and over the post-first-high part
        int_J_fwd2_whole=float(_np.trapezoid(rows['J_fwd2'], tt)),
        int_J_rev2_whole=float(_np.trapezoid(rows['J_rev2'], tt)),
        int_J_fwd2_after_first_high=(float(_np.trapezoid(rows['J_fwd2'][hi[0]:], after_t))
                                     if hi.size else None),
        int_J_rev2_after_first_high=(float(_np.trapezoid(rows['J_rev2'][hi[0]:], after_t))
                                     if hi.size else None),
        max_J_rev2_whole=float(_np.max(rows['J_rev2'])),
        max_J_rev2_after_first_high=(float(_np.max(rows['J_rev2'][hi[0]:]))
                                     if hi.size else None),
        bit2_crossings=[dict(time_h=float(c['time_h']), direction=c['direction'])
                        for c in a['bit2_crossings']],
        carry1_gate_events=[dict(start_h=float(g['start_h']), peak_h=float(g['peak_h']),
                                 end_h=float(g['end_h']))
                            for g in a['carry1_gate_events']],
    )
    print(f"excursion: first high at t={t_first_high} h; returns to low band afterwards="
          f"{returns_low}; min after = {excursion['s2_min_after_first_high']}")

    def bitnum(label):
        """Read a window's bit label: 'x' means the window never committed.

        At n=4 some windows are UNDEFINED, so the label is the string 'x'.  That is a
        result, not an error, and must be recorded as null rather than dropped.
        """
        s = str(label)
        return int(s) if s.isdigit() else None

    windows = [dict(trough_h=w['trough_h'], window_start_h=w['window_start_h'],
                    window_end_h=w['window_end_h'], value=w['value'],
                    bit0=bitnum(w['bits'][0]['label']),
                    bit1=bitnum(w['bits'][1]['label']),
                    bit2=bitnum(w['bits'][2]['label']),
                    raw_labels=[str(b['label']) for b in w['bits']])
               for w in a['read_windows']]
    n_undef_bits = sum(1 for w in windows for b in (w['bit0'], w['bit1'], w['bit2'])
                       if b is None)
    doc = dict(
        purpose=('matched n=4 comparison run for Figure 4-3 Panel A: the frozen '
                 'expression settings with n_A1_gate = 4 instead of 6'),
        status='NEW DETERMINISTIC RUN made for this figure - not a pre-existing artefact',
        gate_exponent=dict(requested=GATE, realised_by_model=realised,
                           checked_by='plausibility_common.build_threebit assert'),
        excursion=excursion,
        excursion_reading=('the figure caption states that bit 2 enters the high band and '
                           'never returns to the low band, and that the reverse write is '
                           'insufficient. Both clauses are measured here on the saved '
                           'column over the WHOLE 600 h run: see '
                           'returns_to_low_band_after_first_high and the J_rev2 '
                           'integrals. This script refuses to finish if the first clause '
                           'does not hold.'),
        state_index={n: index[n] for n in SAVED_STATES},
        state_index_source=('model.state_names, looked up by name and asserted against '
                            'the documented 51-state layout (b0_S=16, b1_S=27, b2_S=44, '
                            'A1=45, F1=46); no bare integer slicing anywhere in this '
                            'script'),
        expression_settings=dict(
            source='plausibility_common.frozen_extension() + frozen_carry0() '
                   '(identical to the frozen n=6 run)',
            extension={k: v for k, v in vars(frozen_extension()).items()}),
        hours=HOURS, sample_min=SAMPLE_MIN, max_step_min=MAX_STEP_MIN,
        n_samples=int(len(sol.t)),
        certified=bool(a['certified']), crossings=len(a['bit2_crossings']),
        gates=len(a['carry1_gate_events']),
        reverse_events=len(a['bit1_reverse_events']),
        steady_reads=a['steady_state']['reads'],
        unlabelled_reads=int(sum(1 for r in a['read_windows']
                                 if r['value'] is None)),
        undefined_bit_labels=n_undef_bits,
        read_windows=windows,
        model_sha256=sha256(FR / 'model_threebit51.py'),
        source_sha256={n: sha256(FR / n) for n in
                       ('model.py', 'model_twobit34.py', 'model_threebit51.py',
                        'verify_threebit51.py')},
        csv_sha256=sha256(csv),
        runtime_s=round(time.perf_counter() - started, 1))
    (DATA / f'{STEM}.json').write_text(
        json.dumps(doc, ensure_ascii=False, indent=2), encoding='utf-8')
    print(f"certified={doc['certified']} crossings={doc['crossings']} "
          f"gates={doc['gates']} reads={doc['steady_reads']} "
          f"unlabelled={doc['unlabelled_reads']}")
    print(f'wrote {csv}')
    print(f'wrote {DATA / (STEM + ".json")}')

    # ---------------------------------------------------------------- self-check
    # Do not return 0 just because two files were written: re-read them off disk and
    # verify the COLUMN SEMANTICS, not merely that the columns exist.  This is the layer
    # whose absence let the previous version ship b2_S/A1/F1 under the names S0/S1/S2.
    print()
    print('--- semantic self-check on the files just written ---')
    from check_state_column_semantics import check_file, windows_from_json
    fails = check_file(csv, windows_from_json(DATA / f'{STEM}.json'), model, index)
    if fails:
        raise SystemExit(f'FATAL: {fails} column-semantics check(s) failed; the data '
                         f'files are NOT trustworthy and must not be plotted.')
    print('column semantics OK')

    # The figure's caption asserts this; if it ever stops holding the figure is wrong.
    if excursion['returns_to_low_band_after_first_high']:
        raise SystemExit(
            'FATAL: bit 2 DOES return to the low band after first entering the high band '
            f"(min {excursion['s2_min_after_first_high']:.4f} <= {BAND_LOW}). The Figure "
            '4-3 caption claims it does not, so the caption would be false. Fix the '
            'caption before publishing this run.')
    print(f"excursion claim holds: after first entering the high band at t="
          f"{excursion['t_first_high_h']:.2f} h, S2 never returns to <= {BAND_LOW} "
          f"(min {excursion['s2_min_after_first_high']:.4f} over "
          f"{excursion['n_samples_after_first_high']} samples)")
    return 0


if __name__ == '__main__':
    sys.exit(main())
