"""Check 2 - identifiability / degeneracy audit.

Two questions the project has never tested explicitly.

Part q  (complex parameters)
    The bit's explicit complex enters through binding = kon*I*R, unbinding =
    koff*C and loss = (koff + delta_C + growth)*C, so in quasi-steady state the
    reverse Hill only sees q = kon/(koff + delta_C + growth).  If trajectories
    with different (kon, koff, delta_C) but equal q agree, then only ONE of the
    three is identifiable and the report must say "q is calibrated", not
    "three complex parameters are calibrated".

Part delay  (carry expression chain)
    The chain is QSS-matched, so the mature protein's steady state equals
    source/gamma regardless of the chain parameters; the chain only contributes
    a delay.  If (mRNA half-life, maturation half-life) pairs with equal
    step-response half-time give equal dynamics, the chain contributes roughly
    one effective degree of freedom, not two.

Pre-registered criteria (fixed before any run)
    Degenerate if max |Delta S| over the late half of the run < 0.02.
    The q part also runs a deliberate q-MISMATCHED control; if that control
    does not move by more than 0.10 the whole test is declared insensitive and
    no conclusion is drawn.
"""
from __future__ import annotations

import argparse
import json
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plausibility_common import (OUT, build_threebit, frozen_carry0,  # noqa: E402
                                 frozen_extension, write_manifest)
from model_twobit34 import CarryExpressionParameters  # noqa: E402

DEGENERATE_TOL = 0.02
CONTROL_MIN_SHIFT = 0.10

# (kon, koff, complex_decay) - first entry is the frozen profile
Q_SET = [
    (0.10, 1.00, 1.00),   # frozen / reference
    (0.05, 0.10, 0.90),
    (0.15, 2.00, 1.00),
    (0.20, 3.00, 1.00),
    (0.25, 4.00, 1.00),
    (0.30, 5.00, 1.00),
]
Q_MISMATCHED = (0.10, 1.00, 0.50)      # q = 0.0667 instead of 0.05

DELAY_GRID = [(m, t) for m in (0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0)
              for t in (10.0, 20.0, 32.5, 45.0, 60.0, 80.0, 100.0, 120.0)]


def chain_t50(mrna_min, mat_min, translation=30.0, gamma=1.0):
    """Half-time of the mature protein's step response for one expression chain."""
    lm = np.log(2) * 60 / mrna_min
    kmat = np.log(2) * 60 / mat_min
    lu = kmat
    source = 1.0
    tx = source * lm * lu / (translation * kmat)

    def f(_t, y):
        m, u, p = y
        return [tx - lm * m, translation * m - lu * u, kmat * u - gamma * p]

    sol = solve_ivp(f, (0.0, 120.0), [0.0, 0.0, 0.0], method='DOP853',
                    rtol=1e-10, atol=1e-14, dense_output=True)
    ts = np.linspace(0.0, 120.0, 24001)
    p = sol.sol(ts)[2]
    target = 0.5 * source / gamma
    return float(ts[int(np.argmax(p >= target))])


def run_threebit(hours, extension=None, carry1=None, sample_min=2.0, max_step_min=2.0):
    model, patch = build_threebit(carry1=carry1 or frozen_carry0(), extension=extension)
    try:
        sol = model.simulate(hours=hours, sample_min=sample_min, max_step_min=max_step_min)
    finally:
        from plausibility_common import restore
        restore(patch)
    return model, sol


def late_gap(a, b, hours, which=(16, 27, 44)):
    """max |Delta state| over the late half, for the three DNA states."""
    n = min(a.y.shape[1], b.y.shape[1])
    cut = n // 2
    return float(max(np.max(np.abs(a.y[i, cut:n] - b.y[i, cut:n])) for i in which))


def part_q(hours):
    base_ext = frozen_extension()
    _, ref = run_threebit(hours, extension=base_ext)
    rows = []
    for kon, koff, dC in Q_SET + [Q_MISMATCHED]:
        ext = replace(base_ext, complex_on_au_inv_h=kon, complex_off_h=koff,
                      complex_decay_h=dC)
        _, sol = run_threebit(hours, extension=ext)
        q = kon / (koff + dC)
        rows.append(dict(part='q', kon=kon, koff=koff, complex_decay=dC, q=q,
                         max_dS_vs_frozen=late_gap(sol, ref, hours),
                         is_control=bool((kon, koff, dC) == Q_MISMATCHED)))
    return rows


def parse_conv(text):
    return [float(x) for x in text.split(',')]


def part_delay(hours, samples='2.0,60;3.0,45;1.0,120;4.0,32.5'):
    """Compare pairs of chain settings chosen to have (nearly) equal t50."""
    pairs = []
    for chunk in samples.split(';'):
        m1, t1, m2, t2 = parse_conv(chunk.replace(';', ','))
        pairs.append(((m1, t1), (m2, t2)))
    rows = []
    ref_model, ref_sol = None, None
    for (a, b) in pairs:
        ta, tb = chain_t50(*a), chain_t50(*b)
        rel = abs(ta - tb) / max(ta, tb)
        c1 = CarryExpressionParameters(mrna_half_life_min=a[0],
                                       activator_maturation_half_life_min=a[1],
                                       repressor_maturation_half_life_min=a[1])
        c2 = CarryExpressionParameters(mrna_half_life_min=b[0],
                                       activator_maturation_half_life_min=b[1],
                                       repressor_maturation_half_life_min=b[1])
        ma, sa = run_threebit(hours, carry1=c1)
        mb, sb = run_threebit(hours, carry1=c2)
        rows.append(dict(part='delay',
                         set_A=f'mRNA{a[0]:g}/mat{a[1]:g}', t50_A=ta,
                         set_B=f'mRNA{b[0]:g}/mat{b[1]:g}', t50_B=tb,
                         rel_t50_gap=rel,
                         max_dS_between_pair=late_gap(sa, sb, hours),
                         degenerate=bool(rel < 0.05 and late_gap(sa, sb, hours) < DEGENERATE_TOL)))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--part', choices=['q', 'delay', 'both'], default='both')
    ap.add_argument('--hours', type=float, default=100.0)
    args = ap.parse_args()

    rows = []
    if args.part in ('q', 'both'):
        rows += part_q(args.hours)
    if args.part in ('delay', 'both'):
        rows += part_delay(args.hours)

    OUT.mkdir(parents=True, exist_ok=True)
    import pandas as pd
    df = pd.DataFrame(rows)
    df.to_csv(OUT / 'identifiability.csv', index=False, encoding='utf-8')

    verdict = {}
    q = df[df.part == 'q']
    if len(q):
        qm = q[~q.is_control].max_dS_vs_frozen.max()
        ctrl = q[q.is_control].max_dS_vs_frozen.max() if q.is_control.any() else float('nan')
        verdict['q_part'] = dict(
            max_shift_within_equal_q=float(qm), control_shift=float(ctrl),
            degenerate=bool(qm < DEGENERATE_TOL),
            test_sensitive=bool(ctrl > CONTROL_MIN_SHIFT))
    d = df[df.part == 'delay']
    if len(d):
        close = d[d.rel_t50_gap < 0.05]
        verdict['delay_part'] = dict(
            pairs_tested=int(len(d)), pairs_within_5pct_t50=int(len(close)),
            max_shift_among_close_pairs=(float(close.max_dS_between_pair.max())
                                         if len(close) else None),
            degenerate=bool(len(close) and close.max_dS_between_pair.max() < DEGENERATE_TOL))
    (OUT / 'identifiability.json').write_text(
        json.dumps(dict(criteria=dict(degenerate_tol=DEGENERATE_TOL,
                                      control_min_shift=CONTROL_MIN_SHIFT),
                        verdict=verdict, rows=rows), ensure_ascii=False, indent=2),
        encoding='utf-8')
    write_manifest()
    pd.set_option('display.width', 220)
    print(df.to_string(index=False))
    print()
    print(json.dumps(verdict, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
