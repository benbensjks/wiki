"""Test whether robust add_growth=True bit0 centres can drive bit1."""
from __future__ import annotations

import json
from dataclasses import asdict

from model import Extension, ROOT
from verify_bit0_part2 import CAND
from verify_carry_delay import simulate as simulate_with_delay
from verify_twobit_causal import analyse


CENTRES = [
    dict(name='ridge_low_conversion', uM_per_au=0.8,
         maturation_half_life_min=14.0),
    dict(name='archived_grid_centre', uM_per_au=1.0,
         maturation_half_life_min=10.0),
]


def compact(result):
    c = result['causal_verdict']
    return dict(cold_start=result['cold_start'], steady_state=result['steady_state'],
                certified=result['certified'], signal_ranges=result['signal_ranges'],
                causal_verdict=c, associations=result['associations'],
                margins=result['margins'])


def main():
    base = asdict(Extension())
    base.update(CAND)
    base.update(complex_on_au_inv_h=1.0, complex_off_h=3.0,
                add_growth=True)
    rows = []
    out = ROOT / 'bit0_results' / 'verification' / 'growth_true_cascade.json'
    for centre in CENTRES:
        cfg = base.copy()
        cfg.update({k: v for k, v in centre.items() if k != 'name'})
        direct = analyse(cfg, hours=300.0, sample_min=2.0)
        delayed = simulate_with_delay(cfg, carry_mrna_half_life_min=2.0,
                                      carry_maturation_half_life_min=30.0,
                                      hours=300.0, sample_min=2.0)
        row = dict(name=centre['name'], parameters=cfg,
                   direct_A0_F0=compact(direct),
                   explicit_A0_F0_expression=compact(delayed))
        rows.append(row)
        out.write_text(json.dumps(rows, ensure_ascii=False, indent=2), encoding='utf-8')
        for label, result in (('direct', direct), ('delay30', delayed)):
            print(f"{centre['name']} {label}: g0={result['signal_ranges']['g0'][1]:.6g} "
                  f"Int1src={result['signal_ranges']['Int1_source_peak']:.6g} "
                  f"steady={result['steady_state']['passed']} certified={result['certified']} "
                  f"seq={result['steady_state']['sequence']}", flush=True)
    print(f'wrote {out}')


if __name__ == '__main__':
    main()
