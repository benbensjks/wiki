"""scan_joint.py — joint scan (krep_tsl × krdf_tsl × k_tag) at low residual baseline."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from multiprocessing import Pool
import zhao_core as zc
from analysis_lib import simulate, make_source_baseline, toggle_score
from couple_oscillator import pulse_times

_, TROUGHS = pulse_times()
BSUB = float(os.environ.get('BASELINE_SUB', '0.26'))
K_REP = float(os.environ.get('K_REP', '0.0186'))
KREPS = [3.0, 10.0, 30.0]
KRDFS = [4.0, 8.0, 40.0, 200.0]
TAGS = [0.0, 4.0, 8.0]


def one(job):
    krep, krdf, tag = job
    P = zc.default_params()
    P.update(K_rep=K_REP, n_rep=3.4, krep_tsl=krep, krdf_tsl=krdf,
             k_tag_int=tag)
    src = make_source_baseline(P['k_dil'], baseline_sub=BSUB, scale=1.0)
    y0 = zc.y0_PB(P, rep_mrna=P['k_tscr'] * P['Dtot'] / P['k_rna'],
                  rep=krep * P['k_tscr'] * P['Dtot'] / P['k_rna'] / P['k_dil'])
    sol = simulate(P, src, y0)
    LRf = zc.LR_total(sol.y.T) / P['Dtot']
    score, H, L, fid, _ = toggle_score(sol.t, LRf, TROUGHS)
    return krep, krdf, tag, score, H, L, fid


if __name__ == '__main__':
    jobs = [(a, b, c) for a in KREPS for b in KRDFS for c in TAGS]
    with Pool(12) as p:
        out = p.map(one, jobs)
    out.sort(key=lambda r: -r[3])
    print(f"{'krep':>5} {'krdf':>6} {'tag':>4} {'score':>7} {'H':>6} {'L':>6} {'fid':>5}")
    for krep, krdf, tag, score, H, L, fid in out:
        Hs = f"{H:.2f}" if H == H else "  - "
        Ls = f"{L:.2f}" if L == L else "  - "
        print(f"{krep:5.0f} {krdf:6.0f} {tag:4.0f} {score:7.3f} {Hs:>6} {Ls:>6} {fid:5.2f}")
