"""verify_audit.py — independent re-verification of the 核查报告 (2026-08-01) disputed claims.

Checks (all at operating point krep=15/krdf=200/tag=12 unless noted):
  1. residual baseline 68 / 60 / 48 copies -> score (claim: 68 & 60 fail, 48 works)
  2. krep=6 and krep=10 (krdf=200/tag=12/bsub=0.26) -> score (claim: both work, contra "krep<10 all fail")
  3. K=12 nM (n=3.4) -> score (claim: collapse but NOT zero)
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from multiprocessing import Pool
import zhao_core as zc
from analysis_lib import simulate, make_source_baseline, toggle_score
from couple_oscillator import pulse_times

_, TROUGHS = pulse_times()

BASE = 0.28  # waveform baseline µM used by all scan scripts

JOBS = [
    # (label, krep, krdf, tag, bsub, K_rep, n_rep)
    ("resid 68cp (bsub=0.167)", 15, 200, 12, BASE - 68/602, 0.0186, 3.4),
    ("resid 60cp (bsub=0.180)", 15, 200, 12, BASE - 60/602, 0.0186, 3.4),
    ("resid 48cp (bsub=0.200)", 15, 200, 12, 0.20,         0.0186, 3.4),
    ("krep=6  (bsub=0.26)",      6, 200, 12, 0.26,         0.0186, 3.4),
    ("krep=10 (bsub=0.26)",     10, 200, 12, 0.26,         0.0186, 3.4),
    ("K=12nM n=3.4 (bsub=0.26)", 15, 200, 12, 0.26,        0.012,  3.4),
    ("K=8nM  n=3.4 (bsub=0.26)", 15, 200, 12, 0.26,        0.008,  3.4),
]


def one(job):
    label, krep, krdf, tag, b, K, n = job
    P = zc.default_params()
    P.update(K_rep=K, n_rep=n, krep_tsl=krep, krdf_tsl=krdf, k_tag_int=tag)
    src = make_source_baseline(P['k_dil'], baseline_sub=b, scale=1.0)
    y0 = zc.y0_PB(P, rep_mrna=P['k_tscr'] * P['Dtot'] / P['k_rna'],
                  rep=krep * P['k_tscr'] * P['Dtot'] / P['k_rna'] / P['k_dil'])
    sol = simulate(P, src, y0)
    LRf = zc.LR_total(sol.y.T) / P['Dtot']
    score, H, L, fid, _ = toggle_score(sol.t, LRf, TROUGHS)
    return label, score, H, L, fid


if __name__ == '__main__':
    with Pool(7) as p:
        out = p.map(one, JOBS)
    for label, score, H, L, fid in out:
        Hs = f"{H:.3f}" if H == H else "  -  "
        Ls = f"{L:.3f}" if L == L else "  -  "
        print(f"{label:28s} score={score:.3f}  H={Hs} L={Ls} fid={fid:.2f}")
