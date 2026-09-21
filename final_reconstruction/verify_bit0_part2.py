"""Part 2: is the bit0 alternation a stored bit or a drive-locked period-2 response?

Also diagnoses why the bit0 -> bit1 carry module produces nothing, and maps how
wide the surviving parameter window is. Zeng parameters are never touched.
"""
from __future__ import annotations

import itertools
import json
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks

from model import Model, Extension, ZENG, ROOT

OUT = ROOT / 'bit0_results' / 'verification'
CONF = ROOT / 'bit0_results' / 'confirmed_candidate'
CAND = dict(uM_per_au=10.0, maturation_half_life_min=20.0,
            complex_on_au_inv_h=0.1, complex_off_h=1.0, add_growth=False)
DT = 2.0 / 60.0


def sim(m, y0, t0, t1, sample_min=2.0, rtol=2e-7, atol=2e-9, max_step_min=2.0):
    ts = np.arange(t0, t1 + sample_min / 120.0, sample_min / 60.0)
    sol = solve_ivp(m.rhs_bit0, (t0, t1), y0, t_eval=ts, method='DOP853',
                    rtol=rtol, atol=atol, max_step=max_step_min / 60.0)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol


def clock_cycles(t, flux):
    ids, _ = find_peaks(flux, prominence=max(1.0, .1 * np.ptp(flux)), distance=int(5 / (t[1] - t[0])))
    return ids


def band_labels(t, S, ids):
    lab = []
    for a, b in zip(ids[:-1], ids[1:]):
        s = S[a:b + 1]
        hi = float(np.mean(s >= 0.7)); lo = float(np.mean(s <= 0.3))
        lab.append(dict(hi=hi, lo=lo, commitment=max(hi, lo),
                        label='1' if hi >= 0.7 else ('0' if lo >= 0.7 else 'x'),
                        S_start=float(s[0]), S_end=float(s[-1]), S_med=float(np.median(s)),
                        S_lo=float(s.min()), S_hi=float(s.max())))
    return lab


def part_H_memory_from_scratch(hours=150.0):
    """Identical worlds differing only in the DNA state at t=0."""
    m = Model(Extension(**CAND))
    out = {}
    traj = {}
    for label, S0 in (('S1', 1.0), ('S0', 0.0), ('Smid', 0.5)):
        y0 = m.initial_state_bit0(cold=False); y0[16] = S0
        sol = sim(m, y0, 0.0, hours)
        traj[label] = sol
        out[label] = dict(S_initial=S0, S_end=float(sol.y[16][-1]),
                          S_median_last_20h=float(np.median(sol.y[16][sol.t >= hours - 20])))
    t = traj['S1'].t
    pairs = {}
    for a, b in itertools.combinations(('S1', 'S0', 'Smid'), 2):
        gap = np.abs(traj[a].y[16] - traj[b].y[16])
        pairs[f'{a}_vs_{b}'] = dict(gap_at_20h=float(np.interp(20, t, gap)),
                                    gap_median_20_40h=float(np.median(gap[(t >= 20) & (t <= 40)])),
                                    gap_max_last_30h=float(np.max(gap[t >= hours - 30])))
    max_res = max(v['gap_max_last_30h'] for v in pairs.values())
    return dict(hours=hours, per_run=out, pairwise=pairs,
                verdict=('PARITY IS STORED: different initial DNA states stay on different branches'
                         if max_res > 0.1 else
                         'PARITY IS NOT STORED: the initial DNA state is forgotten; the alternation is locked to the clock'))


def part_I_kick_test(ref_hours=160.0, kicks=(80.0, 85.0, 90.0, 95.0)):
    """Force the DNA state to the opposite value mid-dwell; does the parity shift?"""
    m = Model(Extension(**CAND))
    y0 = m.initial_state_bit0(cold=False)
    ref = sim(m, y0, 0.0, ref_hours)
    ref_flux = np.array([m.flux(m._expand_bit0(ref.y[:, k])) for k in range(ref.y.shape[1])])
    ids = clock_cycles(ref.t, ref_flux)
    ref_lab = band_labels(ref.t, ref.y[16], ids)
    rows = []
    for tk in kicks:
        pre = sim(m, y0, 0.0, tk)
        y1 = pre.y[:, -1].copy(); y1[16] = 1.0 - y1[16]
        post = sim(m, y1, tk, ref_hours)
        t = np.r_[pre.t, post.t[1:]]; S = np.r_[pre.y[16], post.y[16][1:]]
        lab = band_labels(t, S, ids)
        # parity of every cycle after the kick, compared with the reference
        tail = [(i, lab[i]['label'], ref_lab[i]['label']) for i in range(len(lab) - 4, len(lab))
                if 0 <= i < len(ref_lab)]
        shifted = all(l == ('1' if r == '0' else '0' if r == '1' else 'x') for _, l, r in tail if l != 'x' and r != 'x')
        same = all(l == r for _, l, r in tail)
        rows.append(dict(kick_h=tk, S_before=float(1.0 - y1[16]),
                         S_after_kick=float(y1[16]), tail=tail,
                         parity_flipped=bool(shifted), unchanged=bool(same),
                         commitment_after=float(np.median([lab[i]['commitment'] for i in range(len(lab) - 4, len(lab))]))))
    n_flip = sum(r['parity_flipped'] for r in rows)
    return dict(reference_cycles=[dict(cycle=i, **l) for i, l in enumerate(ref_lab)],
                kicks=rows, kicks_that_flip_parity=n_flip, kicks_that_change_nothing=sum(r['unchanged'] for r in rows),
                verdict=('COUNTING-CAPABLE: a forced flip persists as a parity shift' if n_flip else
                         'NOT COUNTING: a forced flip is erased, the parity is set by the clock'))


def part_J_carry_module():
    """Why does bit1 never move? Read the twobit trajectory and decompose g0."""
    df = pd.read_csv(CONF / 'twobit_candidate_trajectory.csv')
    t = df.time_h.to_numpy(); S0 = df.b0_S.to_numpy()
    A0 = df.A0.to_numpy(); F0 = df.F0.to_numpy()
    H = lambda x, K, n: x**n / (K**n + x**n)
    g0 = H(A0, ZENG['K_A'][0], ZENG['n_A'][0]) * (1 - H(F0, ZENG['K_F'][0], ZENG['n_F'][0]))
    drive = ZENG['alpha_A'][0] * (1 - S0)
    rows = dict(
        A0_range=[float(A0.min()), float(A0.max())], A0_median=float(np.median(A0)),
        F0_range=[float(F0.min()), float(F0.max())], F0_median=float(np.median(F0)),
        a0_hill_term_range=[float(H(A0, ZENG['K_A'][0], ZENG['n_A'][0]).min()),
                            float(H(A0, ZENG['K_A'][0], ZENG['n_A'][0]).max())],
        f0_repression_term_range=[float((1 - H(F0, ZENG['K_F'][0], ZENG['n_F'][0])).min()),
                                  float((1 - H(F0, ZENG['K_F'][0], ZENG['n_F'][0])).max())],
        g0_recomputed_range=[float(g0.min()), float(g0.max())],
        g0_times_alpha_Int1=[float(18 * g0.min()), float(18 * g0.max())],
        F0_response_time_h=float(1 / ZENG['gamma_F'][0]),
        A0_drive_period_h=float(np.median(np.diff(t[clock_cycles(t, df.upstream_flux_uM_h.to_numpy())]))),
        S0_dwell_bands=[float(np.percentile(S0, 5)), float(np.percentile(S0, 95))],
        bit1_mature_Int_max=float(df.b1_I.max()),
    )
    rows['achilles_heel'] = (
        'F0 responds in {:.2f} h while the A0 drive itself has a {:.2f} h period, so the '
        'incoherent arm always tracks the coherent arm: the AND gate never opens '
        '(g0 peaks at {:.2e}, i.e. {:.2e} a.u./h into Int1).'.format(
            rows['F0_response_time_h'], rows['A0_drive_period_h'], rows['g0_recomputed_range'][1],
            rows['g0_times_alpha_Int1'][1]))
    rows['verdict'] = ('CARRY MODULE IS STRUCTURALLY SHUT: g0 stays at noise level, so bit1 cannot move'
                       if rows['g0_recomputed_range'][1] < 0.01 else
                       'carry promoter does open; re-examine bit1')
    return rows


def _region_one(cfg, hours):
    m = Model(Extension(**cfg))
    try:
        sol = m.simulate_bit0(hours=hours, sample_min=2.0, max_step_min=2.0)
        flux = np.array([m.flux(m._expand_bit0(sol.y[:, k])) for k in range(sol.y.shape[1])])
        ids = clock_cycles(sol.t, flux)
        lab = band_labels(sol.t, sol.y[16], ids)
        late = lab[2:] if len(lab) > 6 else lab
        codes = ''.join(l['label'] for l in late)
        commit = float(np.median([l['commitment'] for l in late])) if late else 0.0
        ok = (len(codes) >= 6 and 'x' not in codes
              and all(codes[i + 1] != codes[i] for i in range(len(codes) - 1))
              and commit >= 0.8)
        return dict(**cfg, codes=codes, median_commitment=commit, stored_bit=bool(ok),
                    S_lo=float(sol.y[16].min()), S_hi=float(sol.y[16].max()), error='')
    except Exception as exc:
        return dict(**cfg, codes='', median_commitment=0.0, stored_bit=False, error=repr(exc))


def part_K_region(workers=6, hours=120.0):
    grid = dict(uM_per_au=(6.0, 8.0, 9.0, 10.0, 11.0, 12.0, 14.0, 18.0),
                maturation_half_life_min=(12.0, 16.0, 18.0, 20.0, 22.0, 26.0, 32.0))
    base = asdict(Extension()); base.update(CAND)
    cfgs = []
    for u, mt in itertools.product(grid['uM_per_au'], grid['maturation_half_life_min']):
        c = base.copy(); c['uM_per_au'] = u; c['maturation_half_life_min'] = mt; cfgs.append(c)
    with ProcessPoolExecutor(max_workers=workers) as pool:
        rows = list(pool.map(_region_one, cfgs, itertools.repeat(hours), chunksize=2))
    d = pd.DataFrame(rows)
    d.to_csv(OUT / 'bit0_region_map.csv', index=False)
    ok = d[d.stored_bit]
    pivot = d.pivot_table(index='uM_per_au', columns='maturation_half_life_min', values='stored_bit')
    return dict(hours=hours, points=len(d), stored_bit_points=int(len(ok)),
                stored_region=([] if not len(ok) else dict(
                    uM_per_au=[float(ok.uM_per_au.min()), float(ok.uM_per_au.max())],
                    maturation_half_life_min=[float(ok.maturation_half_life_min.min()),
                                               float(ok.maturation_half_life_min.max())])),
                table=pivot.astype(int).to_dict(),
                frame=rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rep = {}
    for name, fn in (('H_initial_state_memory', part_H_memory_from_scratch),
                     ('I_kick_test', part_I_kick_test),
                     ('J_carry_module', part_J_carry_module)):
        rep[name] = fn()
        (OUT / 'bit0_verification_part2.json').write_text(
            json.dumps(rep, ensure_ascii=False, indent=2), encoding='utf-8')
        print(name, 'done', flush=True)
    print(json.dumps(rep, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
