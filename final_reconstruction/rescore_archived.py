"""Re-score archived trajectories with a dwell-based bit criterion.

The archived runs decoded S at one instant per clock cycle with fixed 0.2/0.8
bands. That rule calls a clean two-level toggle "ambiguous" whenever its dwell
levels sit slightly inside the bands. This script scores the same trajectories
by how much of each clock cycle the DNA actually spends inside one band.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy.signal import find_peaks

from model import ROOT

OUT = ROOT / 'bit0_results' / 'verification'
BAND = 0.7          # a state is read as 0/1 when it sits outside this ball around 0.5
MIN_DWELL = 0.8     # ... for at least this fraction of the clock cycle


def cycles_from_flux(t, flux):
    ids, _ = find_peaks(flux, prominence=max(1.0, .1 * np.ptp(flux)),
                        distance=int(5 / (t[1] - t[0])))
    return ids


def dwell_codes(t, S, ids):
    rows = []
    for n, (a, b) in enumerate(zip(ids[:-1], ids[1:]), 1):
        s = S[a:b + 1]
        hi = float(np.mean(s >= BAND)); lo = float(np.mean(s <= 1 - BAND))
        lab = '1' if hi >= MIN_DWELL else ('0' if lo >= MIN_DWELL else 'x')
        rows.append(dict(cycle=n, label=lab, frac_high=hi, frac_low=lo,
                         commitment=max(hi, lo), S_start=float(s[0]), S_end=float(s[-1]),
                         S_lo=float(s.min()), S_hi=float(s.max())))
    return rows


def score(t, S, ids, drop=2):
    rows = dwell_codes(t, S, ids)
    late = rows[drop:]
    codes = ''.join(r['label'] for r in late)
    commit = float(np.median([r['commitment'] for r in late])) if late else 0.0
    alternates = (len(codes) >= 6 and 'x' not in codes
                  and all(codes[i + 1] != codes[i] for i in range(len(codes) - 1)))
    return dict(codes=codes, median_commitment=commit, alternates=alternates,
                cycles=len(rows), frame=rows)


def rescore(path, fluence_col, s_cols, sample_min=1.0):
    d = pd.read_csv(path)
    t = d.time_h.to_numpy(); flux = d[fluence_col].to_numpy()
    ids = cycles_from_flux(t, flux)
    return {name: score(t, d[col].to_numpy(), ids) for name, col in s_cols.items()}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rep = {}
    rep['archived_full_3bit'] = rescore(
        ROOT / 'results' / 'trajectories.csv', 'C31_translation_flux_uM_h',
        {'bit0': 'b0_S', 'bit1': 'b1_S', 'bit2': 'b2_S'})
    rep['confirmed_bit0_candidate'] = rescore(
        ROOT / 'bit0_results' / 'confirmed_candidate' / 'bit0_180h_trajectory.csv',
        'flux_uM_h', {'bit0': 'b0_S'}, sample_min=2.0)
    rep['twobit_candidate'] = rescore(
        ROOT / 'bit0_results' / 'confirmed_candidate' / 'twobit_candidate_trajectory.csv',
        'upstream_flux_uM_h', {'bit0': 'b0_S', 'bit1': 'b1_S'}, sample_min=2.0)
    print(json.dumps({k: {kk: {kkk: vvv for kkk, vvv in vv.items() if kkk != 'frame'}
                          for kk, vv in v.items()} for k, v in rep.items()},
                     ensure_ascii=False, indent=2))
    (OUT / 'archived_rescore.json').write_text(json.dumps(rep, ensure_ascii=False, indent=2), encoding='utf-8')


if __name__ == '__main__':
    main()
