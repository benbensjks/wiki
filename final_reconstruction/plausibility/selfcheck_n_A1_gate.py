"""Self-check for the n_A1_gate plumbing in plausibility_common.

Run this after ANY change to the model builders.  Every check is an assertion on
measured behaviour, not on source text.  Exits non-zero on the first failure.

Checks
  C1  the named constants have the expected values and the legacy alias is intact
  C2  build_threebit() with no gate argument lands on the FROZEN point (6.0),
      NOT on the ZENG table fallback (4.0) - this is the defect being fixed
  C3  the three resolution rules for build_threebit's gate exponent
  C4  build_decoupled() lands on the frozen point by default
  C5  build_decoupled refuses None instead of silently using the table value
  C6  ZENG is restored after every build (no global leak between calls)
  C7  the constructor path and the global-patch path agree at the selected exponent
  C8  states 0..33 are identical for every gate exponent (the frozen prefix)
  C9  the initial state does not depend on the gate exponent
  C10 the direct-construction subclass (check_rdf_collapse) lands on the frozen point
  C11 no script in this directory constructs the 51-state model without naming the
      gate exponent somewhere in the file (lint against this defect class)

    python selfcheck_n_A1_gate.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plausibility_common import (FROZEN_N_A1, SELECTED_N_A1_GATE,  # noqa: E402
                                 ZENG_N_A1_TABLE, build_decoupled,
                                 build_threebit, frozen_extension, restore)
import model as M  # noqa: E402

RESULTS = []


def check(name, ok, detail=''):
    RESULTS.append((name, bool(ok), detail))
    print('  [%s] %-6s %s' % ('PASS' if ok else 'FAIL', name, detail))
    return bool(ok)


def main():
    print('=== C1 constants ===')
    check('C1a', ZENG_N_A1_TABLE == 4.0, 'ZENG_N_A1_TABLE=%g' % ZENG_N_A1_TABLE)
    check('C1b', SELECTED_N_A1_GATE == 6.0, 'SELECTED_N_A1_GATE=%g' % SELECTED_N_A1_GATE)
    check('C1c', FROZEN_N_A1 == ZENG_N_A1_TABLE,
          'FROZEN_N_A1=%g (legacy alias of the TABLE value)' % FROZEN_N_A1)
    check('C1d', M.ZENG['n_A'][1] == 4.0, 'live ZENG n_A[1]=%g' % M.ZENG['n_A'][1])

    print('=== C2 the defect being fixed ===')
    m, p = build_threebit()
    try:
        got = float(m.n_A1_gate_effective)
        check('C2', got == SELECTED_N_A1_GATE,
              'build_threebit() -> n_A1_gate_effective=%g (was 4.0 before the fix)' % got)
    finally:
        restore(p)

    print('=== C3 resolution rules ===')
    cases = [
        ('no n_A1, no n_A1_gate', dict(), SELECTED_N_A1_GATE),
        ('n_A1=5.0 only (shared slot)', dict(n_A1=5.0), 5.0),
        ('n_A1_gate=8.0 only', dict(n_A1_gate=8.0), 8.0),
        ('n_A1=5.0 + n_A1_gate=6.0', dict(n_A1=5.0, n_A1_gate=6.0), 6.0),
    ]
    for label, kw, want in cases:
        m, p = build_threebit(**kw)
        try:
            got = float(m.n_A1_gate_effective)
            check('C3:' + label, got == want, 'got %g want %g' % (got, want))
        finally:
            restore(p)

    print('=== C4/C5 build_decoupled ===')
    m, p = build_decoupled()
    try:
        got = float(m.n_A1_gate_effective)
        check('C4', got == SELECTED_N_A1_GATE, 'default -> %g' % got)
    finally:
        restore(p)
    try:
        build_decoupled(n_A1_gate=None)
        check('C5', False, 'None was accepted - it must raise')
    except ValueError as exc:
        check('C5', True, 'None raises ValueError: %s' % str(exc)[:60])
    except Exception as exc:                                     # noqa: BLE001
        check('C5', False, 'raised the wrong type: %r' % exc)

    print('=== C6 no global leak ===')
    before = M.ZENG['n_A']
    for kw in (dict(n_A1=5.0), dict(n_A1_gate=8.0)):
        m, p = build_threebit(**kw)
        restore(p)
    m, p = build_decoupled(n_A1_gate=7.0)
    restore(p)
    check('C6', M.ZENG['n_A'] == before, 'ZENG n_A=%s after all builds' % (M.ZENG['n_A'],))

    print('=== C7 constructor path vs global-patch path at the selected exponent ===')
    a, pa = build_threebit(n_A1_gate=SELECTED_N_A1_GATE)
    try:
        sa = a.simulate(hours=60.0, sample_min=2.0, max_step_min=2.0)
    finally:
        restore(pa)
    b, pb = build_decoupled(n_A1_gate=SELECTED_N_A1_GATE)
    try:
        sb = b.simulate(hours=60.0, sample_min=2.0, max_step_min=2.0)
    finally:
        restore(pb)
    gap = float(np.max(np.abs(sa.y - sb.y)))
    check('C7', gap < 1e-9, 'trajectory max|A-B| = %.3e over 60 h' % gap)

    print('=== C8 frozen prefix is exponent-independent ===')
    ref = None
    ok = True
    detail = []
    for ng in (4.0, 5.0, 6.0, 7.0, 8.0):
        m, p = build_threebit(n_A1_gate=ng)
        try:
            d = m.rhs(0.0, sa.y[:, 500])
            if ref is None:
                ref = d[:34].copy()
            else:
                g = float(np.max(np.abs(d[:34] - ref)))
                ok = ok and g == 0.0
                detail.append('n=%g:%s' % (ng, '0' if g == 0.0 else '%.1e' % g))
        finally:
            restore(p)
    check('C8', ok, 'rhs[0:34] gap vs n=4: ' + ' '.join(detail))

    print('=== C9 initial state is exponent-independent ===')
    gaps = []
    for cold in (False, True):
        y4 = None
        worst = 0.0
        for ng in (4.0, 6.0, 8.0):
            m, p = build_threebit(n_A1_gate=ng)
            try:
                y = m.initial_state(cold=cold)
            finally:
                restore(p)
            if y4 is None:
                y4 = y
            worst = max(worst, float(np.max(np.abs(y - y4))))
        gaps.append('cold=%s:%.1e' % (cold, worst))
    check('C9', all(g.endswith('0.0e+00') for g in gaps), ' '.join(gaps))

    print('=== C10 the direct-construction subclass lands on the frozen point ===')
    try:
        from check_rdf_collapse import Bit2RdfProbe
        probe = Bit2RdfProbe(extension=frozen_extension(), carry=None,
                             k_rep=None, n_rep=None)
        got = float(probe.n_A1_gate_effective)
        check('C10', got == SELECTED_N_A1_GATE,
              'Bit2RdfProbe (direct construction) -> %g' % got)
    except Exception as exc:                                     # noqa: BLE001
        check('C10', False, 'could not build the probe: %r' % exc)

    print('=== C11 no script constructs the 51-state model without naming the exponent ===')
    offenders = []
    here = Path(__file__).resolve().parent
    for f in sorted(here.glob('*.py')):
        s = f.read_text(encoding='utf-8')
        constructs = ('ThreeBit51Model(' in s or 'super().__init__(extension=' in s)
        if not constructs:
            continue
        if 'n_A1_gate' not in s:
            offenders.append(f.name)
    check('C11', not offenders,
          'all clear' if not offenders else 'OFFENDERS: ' + ', '.join(offenders))

    failed = [n for n, ok_, _ in RESULTS if not ok_]
    print()
    print('checked %d, failed %d %s' % (len(RESULTS), len(failed),
                                        ('-> ' + ', '.join(failed)) if failed else ''))
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
