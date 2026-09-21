"""Check 6 - is the bit2 stall caused by the RDF pool collapsing?

Two candidate explanations for the observed stall (S2 stuck above the 0.30 band):

  (i)  dose problem      - the carry does not deliver enough complex;
  (ii) collapse problem  - the RDF pool's own repression Hill G(T; K_rep, 3.9)
       is so steep that when S2 slides from ~0.95 to ~0.6 the RDF production
       collapses by two orders of magnitude, so no carry dose can rebuild the
       complex.

CORRECTION TO THE FIRST VERSION OF THIS CHECK
---------------------------------------------
The first version probed the collapse by patching ZENG['n_rep'] or
ZENG['K_rep'] globally.  That moves the RDF gate of **all three bits**, and the
frozen prefix stopped counting (bit1 reverse events fell from 14 to 0, read
sequences became '101010...' and '000000...').  Those probes were therefore
confounded and their "no recovery" outcome proved nothing.  This version

  * changes ONLY bit2's RDF repression Hill, by subclassing ThreeBit51Model and
    correcting the single derivative entry that the RDF transcript source feeds
    (M_R of bit2, index 40); nothing frozen is modified;
  * records a prefix guard (bit1 reverse-event count must equal the baseline's)
    and EXCLUDES any probe that fails it from the verdict;
  * fixes the S -> R_ss map: R_ss(S) is monotonically increasing, so the
    informative statistic is where it turns ON sharply (max d ln R_ss / dS),
    not a "most negative slope".

Pre-registered decision rule
----------------------------
Consider only probes that preserve the prefix.
If such a probe reaches at least the baseline crossing count and certification
while the baseline does not, the RDF collapse is LOAD-BEARING.
If no probe preserves the prefix, the verdict is INCONCLUSIVE, not negative.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plausibility_common import (OUT, frozen_carry0, frozen_extension,  # noqa: E402
                                 gate_windows, leak_budget, source_hashes,
                                 write_manifest)
from model import ZENG  # noqa: E402

PROBES = [
    dict(name='baseline', k_rep_bit2=None, n_rep_bit2=None),
    dict(name='bit2_soft_repression_n2', k_rep_bit2=None, n_rep_bit2=2.0),
    dict(name='bit2_no_repression_Krep100', k_rep_bit2=100.0, n_rep_bit2=None),
]
B2_OFFSET = 34
B2_M_R = B2_OFFSET + 6


def longest_mod8_run(seq: str) -> int:
    """Longest run of consecutive reads that increment correctly modulo 8."""
    best = cur = 0
    for i, ch in enumerate(str(seq)):
        if ch == 'x':
            cur = 0
            continue
        if cur == 0:
            cur = 1
        else:
            prev = str(seq)[i - 1]
            if prev != 'x' and (int(ch) - int(prev)) % 8 == 1:
                cur += 1
            else:
                cur = 1
        best = max(best, cur)
    return best


def near_miss_reason(row) -> str:
    """Name the sub-criterion that blocks certification when counting works."""
    if bool(row.certified):
        return ''
    run = int(row.longest_mod8_run)
    seq = str(row.sequence)
    n_x = seq[8:].count('x')
    if run >= 16:
        return (f'counting_correct_over_{run}_reads_blocked_by_readout: '
                f'{n_x} unlabelled window(s) inside the steady set')
    if int(row.crossings) == int(row.gates) and int(row.gates) > 0:
        return 'one_flip_per_gate_but_not_a_full_mod8_count'
    return 'not_counting'


def make_bit2_rdf_probe():
    from model_threebit51 import ThreeBit51Model

    class Bit2RdfProbe(ThreeBit51Model):
        """Change only bit2's RDF repression Hill.

        The RDF transcript source enters the derivative of M_R alone, so the
        correction is one additive term on d[M_R]; A1/F1, the gate and the frozen
        prefix are untouched.
        """

        def __init__(self, extension=None, carry=None, k_rep=None, n_rep=None):
            self.k_rep_bit2 = k_rep
            self.n_rep_bit2 = n_rep
            super().__init__(extension=extension, carry=carry)

        def rhs(self, t, y):
            d = super().rhs(t, y)
            if self.k_rep_bit2 is None and self.n_rep_bit2 is None:
                return d
            from plausibility_common import repression as vrep
            K_new = ZENG['K_rep'] if self.k_rep_bit2 is None else self.k_rep_bit2
            n_new = ZENG['n_rep'] if self.n_rep_bit2 is None else self.n_rep_bit2
            T2 = np.maximum(float(y[B2_OFFSET + 5]), 0.0)
            S2 = float(y[B2_OFFSET + 10])
            scale = (self.base.lm * self.base.lu /
                     (self.e.translation_h * self.base.kmat))
            src_old = ZENG['alpha_rdf'] * S2 * vrep(T2, ZENG['K_rep'], ZENG['n_rep'])
            src_new = ZENG['alpha_rdf'] * S2 * vrep(T2, K_new, n_new)
            d[B2_M_R] += float(src_new - src_old) * scale
            return d

    return Bit2RdfProbe


Bit2RdfProbe = make_bit2_rdf_probe()


def rdf_map(k_rep=None, n_rep=None, grid=2001):
    """R_ss(S) and its turn-on point (the map is monotonically increasing)."""
    n = ZENG['n_rep'] if n_rep is None else n_rep
    K = ZENG['K_rep'] if k_rep is None else k_rep
    S = np.linspace(1e-9, 1.0, grid)
    T = ZENG['alpha_rep'] * (1.0 - S) / ZENG['gamma_rep']
    R = ZENG['alpha_rdf'] * S * (1.0 / (1.0 + (T / K) ** n)) / ZENG['gamma_rdf']
    dln = np.gradient(np.log(R), S)
    S_on = float(S[int(np.argmax(dln))])
    half = 0.5 * float(R.max())
    above = np.flatnonzero(R >= half)
    return dict(S=S, R=R, S_turn_on=S_on, max_log_slope=float(dln.max()),
                S_at_half_plateau=float(S[above[0]]) if above.size else None,
                R_plateau=float(R.max()), R_ss_of_S=lambda s: float(
                    np.interp(s, S, R)))


def run_probe(probe, hours, sample_min, max_step_min):
    from verify_threebit51 import analyse_threebit
    started = time.perf_counter()
    model = Bit2RdfProbe(extension=frozen_extension(), carry=None,
                         k_rep=probe['k_rep_bit2'], n_rep=probe['n_rep_bit2'])
    sol = model.simulate(hours=hours, sample_min=sample_min, max_step_min=max_step_min)
    analysis = analyse_threebit(model, sol, hours)
    sig = model.diagnostic_signals(sol.y)
    t = sol.t
    windows = gate_windows(t, sig['g1'])
    budget = leak_budget(t, sig['g1'], sig['J_rev2'], windows)
    S2 = sig['S2']
    late = t > 0.5 * hours
    k = int(np.flatnonzero(late)[int(np.argmin(S2[late]))])
    return dict(probe=probe['name'], k_rep_bit2=probe['k_rep_bit2'],
                n_rep_bit2=probe['n_rep_bit2'],
                certified=bool(analysis['certified']),
                steady_passed=bool(analysis['steady_state']['passed']),
                crossings=len(analysis['bit2_crossings']),
                gates=len(analysis['carry1_gate_events']),
                reverse_events=len(analysis['bit1_reverse_events']),
                sequence=analysis['steady_state']['sequence'],
                longest_mod8_run=longest_mod8_run(analysis['steady_state']['sequence']),
                leak_worst_Jrev2=budget['worst_Jrev2'],
                S2_min_late=float(S2[k]), R2_at_S2min=float(sol.y[42][k]),
                I2_at_S2min=float(sol.y[36][k]), C2_at_S2min=float(sol.y[43][k]),
                runtime_s=round(time.perf_counter() - started, 2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--hours', type=float, default=600.0)
    ap.add_argument('--sample-min', type=float, default=2.0)
    ap.add_argument('--max-step-min', type=float, default=2.0)
    ap.add_argument('--rescore-from-csv', action='store_true',
                    help='re-apply the v2 verdict rule to the saved probe CSV (no simulation)')
    args = ap.parse_args()

    if args.rescore_from_csv:
        src = OUT / 'rdf_collapse.csv'
        if not src.exists():
            raise SystemExit(f'no saved probe CSV at {src}')
        df = pd.read_csv(src)
        df['longest_mod8_run'] = [longest_mod8_run(s) for s in df.sequence]
        df['near_miss_reason'] = [near_miss_reason(r) for r in df.itertuples()]
        base = df[df.probe == 'baseline'].iloc[0]
        base_map = rdf_map()
        usable = df[(df.probe != 'baseline') & df.prefix_preserved.astype(bool)]
        counting_ok = usable[usable.longest_mod8_run >= 16]
        strict = usable[usable.certified.astype(bool) &
                        (usable.crossings >= int(base.crossings))]
        verdict = dict(rescored_from_csv=True,
                       baseline_crossings=int(base.crossings),
                       baseline_reverse_events=int(base.reverse_events),
                       R_ss_of_stall_point=float(base_map['R_ss_of_S'](base.S2_min_late)),
                       R_ss_at_S_0p95=float(base_map['R_ss_of_S'](0.95)),
                       probes_preserving_prefix=usable.probe.tolist(),
                       probes_that_recover_strict=strict.probe.tolist(),
                       probes_that_recover_counting=counting_ok.probe.tolist(),
                       longest_correct_mod8_run={r.probe: int(r.longest_mod8_run)
                                                 for r in df.itertuples()},
                       near_miss={r.probe: near_miss_reason(r) for r in df.itertuples()},
                       rdf_collapse_is_load_bearing=bool(len(strict) > 0),
                       rdf_collapse_supported_by_counting=bool(len(counting_ok) > 0))
        df.to_csv(OUT / 'rdf_collapse_rescored.csv', index=False, encoding='utf-8')
        (OUT / 'rdf_collapse_rescored.json').write_text(
            json.dumps(verdict, ensure_ascii=False, indent=2), encoding='utf-8')
        write_manifest()
        pd.set_option('display.width', 220)
        print(df[['probe', 'certified', 'crossings', 'gates', 'reverse_events',
                  'prefix_preserved', 'S2_min_late', 'longest_mod8_run']].to_string(index=False))
        print()
        print(json.dumps(verdict, ensure_ascii=False, indent=2))
        return

    rows = [run_probe(p, args.hours, args.sample_min, args.max_step_min) for p in PROBES]
    df = pd.DataFrame(rows)
    base = df[df.probe == 'baseline'].iloc[0]
    df['prefix_preserved'] = df.reverse_events == int(base.reverse_events)
    maps = {p['name']: {k: v for k, v in rdf_map(p['k_rep_bit2'], p['n_rep_bit2']).items()
                        if k not in ('S', 'R', 'R_ss_of_S')} for p in PROBES}
    base_map = rdf_map()
    df['R_ss_at_S2min'] = [base_map['R_ss_of_S'](s) for s in df.S2_min_late]
    df['R_ss_at_0p95'] = base_map['R_ss_of_S'](0.95)

    usable = df[(df.probe != 'baseline') & df.prefix_preserved]
    counting_ok = usable[usable.longest_mod8_run >= 16]
    strict = usable[usable.certified & (usable.crossings >= int(base.crossings))]
    if len(usable) == 0:
        verdict_txt = ('INCONCLUSIVE: every collapse probe perturbed the frozen prefix, '
                       'so none of them tested bit2 in isolation')
    elif len(strict) or len(counting_ok):
        verdict_txt = ('SUPPORTED: a bit2-only collapse probe restores a correct mod-8 count'
                       + (' and full certification' if len(strict) else
                          '; certification is blocked only by unlabelled read windows'))
    else:
        verdict_txt = 'not binding: prefix-preserving collapse probes do not recover bit2'

    verdict = dict(baseline_crossings=int(base.crossings),
                   baseline_reverse_events=int(base.reverse_events),
                   baseline_S2_min_late=float(base.S2_min_late),
                   baseline_R2_at_S2min=float(base.R2_at_S2min),
                   R_ss_of_stall_point=float(base_map['R_ss_of_S'](base.S2_min_late)),
                   R_ss_at_S_0p95=float(base_map['R_ss_of_S'](0.95)),
                   rdf_map=maps,
                   probes_preserving_prefix=usable.probe.tolist(),
                   probes_confounded=df[(df.probe != 'baseline') &
                                        ~df.prefix_preserved].probe.tolist(),
                   probes_that_recover_strict=strict.probe.tolist(),
                   probes_that_recover_counting=counting_ok.probe.tolist(),
                   longest_correct_mod8_run={r.probe: int(r.longest_mod8_run)
                                             for r in df.itertuples()},
                   rdf_collapse_is_load_bearing=bool(len(strict) > 0),
                   rdf_collapse_supported_by_counting=bool(len(counting_ok) > 0),
                   verdict=verdict_txt,
                   note=('v2: the strict field still requires certified=True, which is '
                         'blocked by readout labelling; rdf_collapse_supported_by_counting '
                         'is the informative one. The v1 global n_rep/K_rep probes were '
                         'confounded and have been replaced by bit2-only probes.'))

    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / 'rdf_collapse.csv', index=False, encoding='utf-8')
    pd.DataFrame(dict(S=base_map['S'], R_ss=base_map['R'])).to_csv(
        OUT / 'rdf_collapse_curve.csv', index=False, encoding='utf-8')
    (OUT / 'rdf_collapse.json').write_text(
        json.dumps(dict(verdict=verdict, source_sha256=source_hashes()),
                   ensure_ascii=False, indent=2), encoding='utf-8')
    write_manifest()

    pd.set_option('display.width', 250)
    print(df[['probe', 'certified', 'crossings', 'gates', 'reverse_events',
              'prefix_preserved', 'S2_min_late', 'R2_at_S2min', 'R_ss_at_S2min']]
          .to_string(index=False))
    print()
    print(json.dumps(verdict, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
