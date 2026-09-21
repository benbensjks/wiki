"""Part 4: how wide is the toggle region, and does it survive numerics?

Reports both criteria for every point:
  commitment  - per cycle, >=80% of the cycle inside one 0.7/0.3 band, bands alternate
  instant     - the archived rule: S sampled at the flux trough, >=0.8 / <=0.2
"""
from __future__ import annotations

import itertools
import json
from dataclasses import asdict

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks

from model import Model, Extension, ROOT
from verify_bit0_part2 import OUT, CAND, clock_cycles, band_labels


def dual_code(m, sol):
    flux = np.array([m.flux(m._expand_bit0(sol.y[:, k])) for k in range(sol.y.shape[1])])
    ids = clock_cycles(sol.t, flux)
    lab = band_labels(sol.t, sol.y[16], ids)
    late = lab[2:] if len(lab) > 6 else lab
    codes = ''.join(l['label'] for l in late)
    commit = float(np.median([l['commitment'] for l in late])) if late else 0.0
    ok_commit = bool(len(codes) >= 6 and 'x' not in codes
                     and all(codes[i + 1] != codes[i] for i in range(len(codes) - 1)) and commit >= 0.8)
    troughs = np.array([a + int(np.argmin(flux[a:b])) for a, b in zip(ids[:-1], ids[1:])])
    vals = sol.y[16][troughs]
    ic = ''.join('1' if v >= 0.8 else ('0' if v <= 0.2 else 'x') for v in vals[2:])
    ok_instant = bool(len(ic) >= 5 and 'x' not in ic[-5:]
                      and all(ic[-5:][i + 1] != ic[-5:][i] for i in range(4)))
    return dict(codes_commitment=codes, median_commitment=commit, stored_commitment=ok_commit,
                codes_instant=ic, stored_instant=ok_instant,
                dwell_low=float(np.percentile(sol.y[16], 2)), dwell_high=float(np.percentile(sol.y[16], 98)),
                S_lo=float(sol.y[16].min()), S_hi=float(sol.y[16].max()))


def run_one(cfg, hours=120.0, tol=None):
    m = Model(Extension(**cfg))
    rtol, atol, ms = tol or (2e-7, 2e-9, 2.0)
    ts = np.arange(0.0, hours + 2.0 / 120.0, 2.0 / 60.0)
    sol = solve_ivp(m.rhs_bit0, (0, hours), m.initial_state_bit0(False), t_eval=ts,
                    method='DOP853', rtol=rtol, atol=atol, max_step=ms / 60.0)
    if not sol.success:
        raise RuntimeError(sol.message)
    return {**cfg, **dual_code(m, sol)}


def part_N_kon_koff():
    base = asdict(Extension()); base.update(CAND)
    rows = []
    for kon, koff in itertools.product((0.03, 0.1, 0.3, 1.0), (0.3, 1.0, 3.0, 10.0)):
        c = base.copy(); c['complex_on_au_inv_h'] = kon; c['complex_off_h'] = koff
        rows.append(run_one(c))
        print(f"kon={kon:5.2f} koff={koff:5.1f} commit={rows[-1]['codes_commitment']} "
              f"({rows[-1]['median_commitment']:.2f}) instant={rows[-1]['codes_instant']}", flush=True)
        pd.DataFrame(rows).to_csv(OUT / 'bit0_kon_koff_map.csv', index=False)
    return dict(points=len(rows),
                stored_by_commitment=int(sum(r['stored_commitment'] for r in rows)),
                stored_by_instant=int(sum(r['stored_instant'] for r in rows)))


def part_O_conversion_sweep():
    base = asdict(Extension()); base.update(CAND)
    vals = (10.0, 12.0, 16.0, 20.0, 30.0, 40.0, 60.0, 100.0)
    rows = []
    for u in vals:
        c = base.copy(); c['uM_per_au'] = u
        rows.append(run_one(c))
        print(f"uM_per_au={u:6.1f} commit={rows[-1]['codes_commitment']} ({rows[-1]['median_commitment']:.2f}) "
              f"dwell={rows[-1]['dwell_low']:.3f}/{rows[-1]['dwell_high']:.3f} instant={rows[-1]['codes_instant']}",
              flush=True)
        pd.DataFrame(rows).to_csv(OUT / 'bit0_conversion_sweep.csv', index=False)
    ok = [r['uM_per_au'] for r in rows if r['stored_commitment']]
    return dict(points=len(rows), stored_by_commitment=int(len(ok)),
                stored_uM_per_au=ok, rows=rows)


def part_P_growth():
    base = asdict(Extension()); base.update(CAND)
    rows = []
    for g in (False, True):
        c = base.copy(); c['add_growth'] = g
        rows.append(run_one(c))
        print(f"add_growth={g} commit={rows[-1]['codes_commitment']} ({rows[-1]['median_commitment']:.3f}) "
              f"dwell={rows[-1]['dwell_low']:.3f}/{rows[-1]['dwell_high']:.3f}", flush=True)
    return dict(rows=rows,
                verdict=('toggle survives only when Zeng gamma is read as the total clearance '
                         'including growth dilution' if rows[0]['stored_commitment'] and not rows[1]['stored_commitment']
                         else 'inspect rows'))


def part_Q_numerics():
    base = asdict(Extension()); base.update(CAND)
    tols = {'default_2e-7_2min': (2e-7, 2e-9, 2.0),
            'tight_1e-9_1min': (1e-9, 1e-11, 1.0),
            'tighter_1e-11_0p25min': (1e-11, 1e-13, 0.25)}
    rows = {}
    for name, tol in tols.items():
        rows[name] = run_one(base, tol=tol)
        print(name, rows[name]['codes_commitment'], round(rows[name]['median_commitment'], 4),
              round(rows[name]['dwell_low'], 4), round(rows[name]['dwell_high'], 4), flush=True)
    codes = {r['codes_commitment'] for r in rows.values()}
    dwell = max(abs(rows[a]['dwell_high'] - rows[b]['dwell_high']) for a in rows for b in rows)
    return dict(rows=rows, codes_agree=len(codes) == 1, max_dwell_difference=float(dwell),
                verdict=('NUMERICALLY CONVERGED' if len(codes) == 1 and dwell < 0.02 else 'TOLERANCE-SENSITIVE'))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rep = {}
    for name, fn in (('N_kon_koff_region', part_N_kon_koff),
                     ('O_conversion_sweep', part_O_conversion_sweep),
                     ('P_growth_interpretation', part_P_growth),
                     ('Q_numerical_convergence', part_Q_numerics)):
        rep[name] = fn()
        (OUT / 'bit0_verification_part4.json').write_text(json.dumps(rep, ensure_ascii=False, indent=2), encoding='utf-8')
        print(name, 'done', flush=True)
    print(json.dumps(rep, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
