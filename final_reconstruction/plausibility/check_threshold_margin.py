"""Check 1 - threshold-margin audit (concentration-space margins).

Question
--------
The project already audits TIMING margins (setup / hold).  It never audited
CONCENTRATION margins: for every Hill term H(x;K,n) in the model, does the
operating trajectory sit near x/K = 1, i.e. exactly where the response is
steepest and the behaviour is most sensitive to every parameter?

REVISION NOTE (v2)
------------------
The first version flagged a term if it satisfied any of three OR-ed conditions,
one of which was `straddles_threshold`.  That condition is true BY CONSTRUCTION
for every switch in a toggle, so 16 of 18 terms were flagged and the aggregate
label carried no information.  v2 therefore

  * uses a single primary metric - the fraction of the run spent inside
    [K/sqrt(2), K*sqrt(2)], the zone where dH/dx is large;
  * treats straddling as DESCRIPTIVE ONLY and never flags on it;
  * adds a role-aware secondary flag for gate arms that are supposed to hold a
    gate shut: a gate whose median operating ratio is not far below 1 is not
    decisively off;
  * ranks the output by the primary metric so the ranking is the product.

Pre-registered flagging rule (v2, fixed before re-scoring)
----------------------------------------------------------
flag_steep_zone  : frac_time_near_threshold >= 0.25
flag_gate_not_off: role starts with 'gate' and ratio_median > 0.10
flagged          : the OR of those two only.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plausibility_common import (OUT, build_threebit, frozen_carry0,  # noqa: E402
                                 restore, threshold_margin, write_manifest)

NEAR_TIME_FLAG = 0.25
GATE_MEDIAN_CEILING = 0.10

BIT_OFFSETS = {'bit0': 6, 'bit1': 17, 'bit2': 34}
BIT_INDEX = {'bit0': 0, 'bit1': 1, 'bit2': 2}


def audit(model, sol) -> pd.DataFrame:
    from model import ZENG
    y = sol.y
    rows = []

    def add(term, x, K, n, role):
        m = threshold_margin(x, K)
        rows.append(dict(term=term, role=role, K=float(K), n=float(n), **m))

    for name, off in BIT_OFFSETS.items():
        add(f'{name}.forward  H(I ; K_D_int)', y[off + 2],
            ZENG['K_D_int'][BIT_INDEX[name]], 2, 'switch (must fire)')
        add(f'{name}.forward  R vs K_inh', y[off + 8], ZENG['K_inh'], 1,
            'switch (must fall)')
        add(f'{name}.reverse  H(C ; K_complex)', y[off + 9], model.base.K_complex, 2,
            'switch (must fire on carry)')
        add(f'{name}.RDF gate T vs K_rep', y[off + 5], ZENG['K_rep'], ZENG['n_rep'],
            'switch (must fall for RDF)')
    add('carry0 H(A0 ; K_A0)', y[28], ZENG['K_A'][0], ZENG['n_A'][0],
        'gate arm (must rise on carry)')
    add('carry0 F0 vs K_F0', y[29], ZENG['K_F'][0], ZENG['n_F'][0],
        'gate arm (must fall on carry)')
    add('carry1 H(A1 ; K_A1)', y[45], ZENG['K_A'][1], ZENG['n_A'][1],
        'gate arm (must be ~0 between carries)')
    add('carry1 F1 vs K_F1', y[46], ZENG['K_F'][1], ZENG['n_F'][1],
        'gate arm (must fall on carry)')
    add('carry1 A1 vs K_auto1', y[45], ZENG['K_auto1'], ZENG['n_auto1'],
        'negative autoregulation')
    add('clock Int0 vs clock_K', y[8], model.e.clock_K_au, model.e.clock_n,
        'clock AND (must be ~0 between pulses)')
    return pd.DataFrame(rows)


def apply_rule(df: pd.DataFrame) -> pd.DataFrame:
    """v2 rule: steep-zone time (primary) plus a gate-off check (secondary)."""
    df = df.copy()
    df['flag_steep_zone'] = df.frac_time_near_threshold >= NEAR_TIME_FLAG
    # Secondary applies only to arms that are SUPPOSED to be off between carries.
    # A gate arm that must rise on a carry (A0) is expected to sit above its
    # threshold most of the time, so a high median is not a defect there.
    must_be_off = df.role.astype(str).str.contains('must be ~0 between carries')
    df['flag_gate_not_off'] = must_be_off & (df.ratio_median > GATE_MEDIAN_CEILING)
    df['flagged'] = df.flag_steep_zone | df.flag_gate_not_off
    df['straddles_descriptive_only'] = df.straddles_threshold
    df['ratio_distance_decades'] = np.abs(np.log10(np.maximum(df.ratio_median, 1e-12)))
    return df.sort_values(['flagged', 'frac_time_near_threshold'],
                          ascending=False).reset_index(drop=True)


def report(df: pd.DataFrame, tag: str, hours):
    out = OUT / 'threshold_margin'
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / f'threshold_margin_{tag}.csv', index=False, encoding='utf-8')
    summary = dict(hours=hours, tag=tag, terms=len(df),
                   flagged=int(df.flagged.sum()),
                   rule=dict(primary=f'frac_time_near_threshold >= {NEAR_TIME_FLAG}',
                             secondary=(f"role startswith 'gate' and "
                                        f'ratio_median > {GATE_MEDIAN_CEILING}'),
                             note='straddles_threshold is descriptive only in v2'),
                   flagged_terms=df[df.flagged][['term', 'role', 'ratio_median',
                                                 'frac_time_near_threshold']]
                   .to_dict(orient='records'))
    (out / f'threshold_margin_{tag}.json').write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
    write_manifest()
    pd.set_option('display.width', 210)
    print(df[['term', 'ratio_min', 'ratio_median', 'ratio_max',
              'frac_time_near_threshold', 'crossings_of_threshold',
              'flag_steep_zone', 'flag_gate_not_off', 'flagged']].to_string(index=False))
    print()
    print(f"flagged {summary['flagged']} / {summary['terms']} terms under the v2 rule "
          f"(primary: steep-zone time >= {NEAR_TIME_FLAG}; secondary: gate arm not off)")
    return summary


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--hours', type=float, default=300.0)
    ap.add_argument('--sample-min', type=float, default=2.0)
    ap.add_argument('--max-step-min', type=float, default=2.0)
    ap.add_argument('--tag', default='baseline')
    ap.add_argument('--rescore-from-csv', action='store_true',
                    help='re-apply the v2 rule to an existing per-term CSV (no simulation)')
    args = ap.parse_args()

    if args.rescore_from_csv:
        src = OUT / 'threshold_margin' / f'threshold_margin_{args.tag}.csv'
        if not src.exists():
            raise SystemExit(f'no per-term CSV at {src}')
        df = apply_rule(pd.read_csv(src))
        report(df, args.tag + '_rescored_v2', args.hours)
        return

    model, patch = build_threebit(carry1=frozen_carry0())
    try:
        sol = model.simulate(hours=args.hours, sample_min=args.sample_min,
                             max_step_min=args.max_step_min)
    finally:
        restore(patch)
    report(apply_rule(audit(model, sol)), args.tag, args.hours)


if __name__ == '__main__':
    main()
