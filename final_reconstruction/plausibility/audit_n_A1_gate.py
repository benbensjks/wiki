"""Exhaustive audit of what `n_A1_gate` does to the 51-state model.

Purpose
-------
`n_A1_gate` is an added interface parameter that changes how a PUBLISHED table
entry (ZENG['n_A'][1]) is used: it splits that single Hill exponent into a
gate-arm exponent and an F1-production exponent.  Before anything is built on
top of it, four things must be settled by measurement rather than by reading:

  Q1  Where is ZENG['n_A'][1] actually read?
  Q2  Which state derivatives does the gate reach?  Is index 34 the ONLY channel?
  Q3  Is the `d[34] += (g1_gate - g1_frozen) * scale` correction EXACT, or is it
      a linearisation that only happens to be close?
  Q4  Does the INITIAL STATE depend on the gate exponent?  If it does, then the
      class-parameter model starts from a state consistent with one exponent and
      integrates with another - an internal inconsistency that must be recorded.
  Q5  Is the frozen prefix (states 0..33) really untouched?

Every question is answered by running the model, not by inspecting the source.
Read-only with respect to the frozen artefacts: this script only constructs
models in memory and evaluates derivatives.

    python audit_n_A1_gate.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plausibility_common import (OUT, frozen_carry0, frozen_extension,  # noqa: E402
                                 restore, sha256, source_hashes)
from model_threebit51 import (STATE_NAMES_51, ThreeBit51Model,  # noqa: E402
                              ThreeBitCarryParameters)
from model_twobit34 import CarryExpressionParameters  # noqa: E402
import model as M  # noqa: E402

MODEL_FILE = Path(__file__).resolve().parents[1] / 'model_threebit51.py'
CARGO = dict(mrna_half_life_min=2.0, activator_maturation_half_life_min=32.5,
             repressor_maturation_half_life_min=32.5)


def build(n_gate):
    return ThreeBit51Model(extension=frozen_extension(),
                           carry=ThreeBitCarryParameters(carry0=frozen_carry0(),
                                                         carry1=CarryExpressionParameters(**CARGO)),
                           n_A1_gate=n_gate)


def probes(n_probes=5):
    """Realistic states taken from a real n=6 trajectory."""
    m = build(6.0)
    sol = m.simulate(hours=200.0, sample_min=2.0, max_step_min=2.0)
    idx = np.linspace(1000, len(sol.t) - 1, n_probes).astype(int)
    return m, sol.t[idx], [sol.y[:, k].copy() for k in idx]


def map43_to_51(i):
    """43-layout index -> 51-state index, per ThreeBit51Model._embed_base43.

    out[:28] = y[:28];  out[28:39] = y[34:45];  out[39:41] = y[28:30];  out[41:43] = y[45:47]
    so the inverse is piecewise, NOT a constant offset:
        0..27   -> 0..27
        28..38  -> 34..44      (the bit2 eleven-state module)
        39,40   -> 28,29       (A0, F0)
        41,42   -> 45,46       (A1, F1)
    """
    if i <= 27:
        return i
    if 28 <= i <= 38:
        return i + 6
    if i in (39, 40):
        return i - 11
    if i in (41, 42):
        return i + 4
    raise ValueError(i)


# Which 43-layout entries are actually CARRIED into the 51-state derivative.
# ThreeBit51Model.rhs takes d[34:45] = d43[28:39] and then overwrites A1/F1 with
# its own _expression()/_carry1_sources() calls, so d43[41] and d43[42] are
# DISCARDED.  That is exactly why the F1 production exponent stays pinned.
USED_IN_51_RHS = set(range(28, 39))


def q2_channel(y_states):
    """Which entries of the 43-layout right-hand side does the gate reach?"""
    base = M.Model(frozen_extension())
    rows = []
    for y in y_states:
        y43 = ThreeBit51Model._embed_base43(y)
        saved = M.ZENG['n_A']
        M.ZENG['n_A'] = (saved[0], 4.0)
        try:
            d4 = base.rhs(0.0, y43).copy()
        finally:
            M.ZENG['n_A'] = saved
        M.ZENG['n_A'] = (saved[0], 6.0)
        try:
            d6 = base.rhs(0.0, y43).copy()
        finally:
            M.ZENG['n_A'] = saved
        diff = np.abs(d4 - d6)
        hit = [int(i) for i in np.flatnonzero(diff > 0)]
        rows.append(dict(
            max_gap=float(diff.max()), indices_touched_43=hit,
            indices_touched_51=[map43_to_51(i) for i in hit],
            names_51=[STATE_NAMES_51[map43_to_51(i)] for i in hit],
            carried_into_51_rhs=[i in USED_IN_51_RHS for i in hit],
            discarded=[map43_to_51(i) for i in hit if i not in USED_IN_51_RHS]))
    return rows


def q3_exactness(y_states):
    """Is the d[34] correction exact, or a linearisation?

    Independently measures the LOCAL SLOPE of the bit2 M_I source with respect to
    the gate value (by nudging A1 slightly), and compares it with the `scale`
    factor the correction uses.
    """
    m4, m6 = build(4.0), build(6.0)
    base = m6.base
    scale = (M.ZENG['alpha_Int'][1] * base.lm * base.lu /
             (m6.e.translation_h * base.kmat))
    rows = []
    for y in y_states:
        y43 = ThreeBit51Model._embed_base43(y)
        # g1 and the M_I source at two nearby A1 values, ZENG at the TABLE value
        saved = M.ZENG['n_A']
        M.ZENG['n_A'] = (saved[0], 4.0)
        try:
            vals = []
            for f in (1.0, 1.0005):
                y43b = y43.copy()
                y43b[41] = y43b[41] * f
                _, g1, _ = base.carry_promoters(y43b)
                d = base.rhs(0.0, y43b)
                vals.append((float(np.asarray(g1)), float(d[28])))
            slope_num = (vals[1][1] - vals[0][1]) / (vals[1][0] - vals[0][0])
        finally:
            M.ZENG['n_A'] = saved
        # what the correction actually applies
        d4 = m4.rhs(0.0, y)
        d6 = m6.rhs(0.0, y)
        measured = float(d6[34] - d4[34])
        _, g1_4, _ = base.carry_promoters(y43)          # ZENG restored -> n=4
        g1_6 = float(m6.diagnostic_signals(y[:, None])['g1'][0])
        predicted_scale = float(g1_6 - g1_4) * scale
        rows.append(dict(
            numerical_source_slope=slope_num, correction_scale=float(scale),
            relative_slope_error=abs(slope_num - scale) / max(abs(scale), 1e-300),
            measured_d34_change=measured, predicted_d34_change=predicted_scale,
            relative_d34_error=(abs(measured - predicted_scale) /
                                max(abs(predicted_scale), 1e-300)),
            other_indices_changed=[int(i) for i in np.flatnonzero(
                np.abs(d6 - d4) > 1e-14)]))
    return rows, float(scale)


def q4_initial_state():
    m4, m6 = build(4.0), build(6.0)
    rows = []
    for cold in (False, True):
        y4, y6 = m4.initial_state(cold=cold), m6.initial_state(cold=cold)
        d = np.abs(y4 - y6)
        rows.append(dict(cold=cold, max_gap=float(d.max()),
                         indices_changed=[int(i) for i in np.flatnonzero(d > 0)],
                         names_changed=[STATE_NAMES_51[i] for i in np.flatnonzero(d > 0)]))
    return rows


def q5_prefix(y_states):
    m4, m6 = build(4.0), build(6.0)
    rows = []
    for y in y_states:
        d4, d6 = m4.rhs(0.0, y), m6.rhs(0.0, y)
        rows.append(dict(prefix_gap=float(np.max(np.abs(d4[:34] - d6[:34]))),
                         tail_gap=float(np.max(np.abs(d4[34:] - d6[34:])))))
    return rows


def q6_effective_semantics():
    """Confirm the fallback and the equality shortcut."""
    rows = []
    for ng in (None, 4.0, 6.0):
        m = build(ng)
        rows.append(dict(requested=ng, effective=m.n_A1_gate_effective,
                         correction_branch_active=bool(
                             m.n_A1_gate_effective != M.ZENG['n_A'][1])))
    return rows


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    _, times, ys = probes()
    report = dict(
        model_sha256=sha256(MODEL_FILE), source_sha256=source_hashes(),
        zeng_n_A=tuple(float(x) for x in M.ZENG['n_A']),
        probe_times_h=[float(t) for t in times],
        Q1_read_sites={
            'note': 'static: ZENG["n_A"][1] is read in exactly two functional places',
            'gate_arm': 'model.py:120  g1 = act(A1, K_A[1], ZENG["n_A"][1]) * rep(F1, ...) * clock',
            'f1_production': 'model_threebit51.py:96  source_F1 = alpha_F[1] * act(A1, K_A[1], ZENG["n_A"][1])',
            'plus': ['model_threebit51.py:90 (effective property)',
                     'model_threebit51.py:143 (correction branch condition)'],
        },
        Q2_gate_channel=q2_channel(ys),
        Q4_initial_state=q4_initial_state(),
        Q5_prefix=q5_prefix(ys),
        Q6_effective_semantics=q6_effective_semantics(),
    )
    rows, scale = q3_exactness(ys)
    report['Q3_correction_exactness'] = dict(
        correction_scale_used=scale, rows=rows,
        max_relative_slope_error=max(r['relative_slope_error'] for r in rows),
        max_relative_d34_error=max(r['relative_d34_error'] for r in rows),
        verdict=('EXACT' if all(r['relative_d34_error'] < 1e-9 for r in rows) else 'APPROXIMATE'))
    (OUT / 'n_A1_gate_audit.json').write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')

    print('ZENG n_A =', report['zeng_n_A'])
    print()
    print('Q2  gate reaches which 43-layout derivative entries:')
    for r in report['Q2_gate_channel']:
        print('    gap=%.3e  43-idx=%s -> 51-idx %s %s' % (
            r['max_gap'], r['indices_touched_43'], r['indices_touched_51'], r['names_51']))
        print('        carried into the 51-state rhs: %s   discarded: %s' % (
            r['carried_into_51_rhs'], r['discarded'] or 'none'))
    print()
    print('Q3  correction exactness (numerical source slope vs the scale used):')
    print('    scale used by the correction = %.10g' % scale)
    for r in rows:
        print('    slope_num=%.10g  rel_err=%.3e | d34 measured=%.6e predicted=%.6e rel_err=%.3e | other idx changed=%s'
              % (r['numerical_source_slope'], r['relative_slope_error'],
                 r['measured_d34_change'], r['predicted_d34_change'],
                 r['relative_d34_error'], r['other_indices_changed']))
    print('    ->', report['Q3_correction_exactness']['verdict'])
    print()
    print('Q4  initial state depends on the gate exponent?')
    for r in report['Q4_initial_state']:
        print('    cold=%-5s max_gap=%.3e  changed=%s' % (
            r['cold'], r['max_gap'], r['names_changed'] or 'none'))
    print()
    print('Q5  frozen prefix under a gate change:')
    for r in report['Q5_prefix']:
        print('    prefix(0:34) gap=%.3e   tail(34:51) gap=%.3e' % (r['prefix_gap'], r['tail_gap']))
    print()
    print('Q6  effective semantics:')
    for r in report['Q6_effective_semantics']:
        print('    requested=%-5s effective=%-6s correction_branch_active=%s' % (
            r['requested'], r['effective'], r['correction_branch_active']))
    print()
    print('wrote', OUT / 'n_A1_gate_audit.json')


if __name__ == '__main__':
    main()
