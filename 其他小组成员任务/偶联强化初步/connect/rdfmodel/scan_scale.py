"""scan_scale.py — add Int-promoter strength (scale) to lower the RDF requirement."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from multiprocessing import Pool
import zhao_core as zc
from analysis_lib import simulate, make_source_baseline, toggle_score
from couple_oscillator import pulse_times

_, TROUGHS = pulse_times()
BSUB = 0.26
KREPS = [15.0, 30.0]
KRDFS = [50.0, 100.0, 200.0]
TAGS = [8.0, 12.0]
SCALES = [1.0, 2.0, 4.0]


def one(job):
    krep, krdf, tag, sc = job
    P = zc.default_params()
    P.update(K_rep=0.0186, n_rep=3.4, krep_tsl=krep, krdf_tsl=krdf,
             k_tag_int=tag)
    src = make_source_baseline(P['k_dil'], baseline_sub=BSUB, scale=sc)
    y0 = zc.y0_PB(P, rep_mrna=P['k_tscr'] * P['Dtot'] / P['k_rna'],
                  rep=krep * P['k_tscr'] * P['Dtot'] / P['k_rna'] / P['k_dil'])
    sol = simulate(P, src, y0)
    Y = sol.y.T
    LRf = zc.LR_total(Y) / P['Dtot']
    RDF = zc.rdf_total(Y)
    score, H, L, fid, _ = toggle_score(sol.t, LRf, TROUGHS)
    return krep, krdf, tag, sc, score, H, L, fid, RDF.max()


if __name__ == '__main__':
    jobs = [(a, bb, c, d) for a in KREPS for bb in KRDFS for c in TAGS for d in SCALES]
    with Pool(12) as p:
        out = p.map(one, jobs)
    out.sort(key=lambda r: -r[4])
    print(f"{'krep':>5} {'krdf':>6} {'tag':>4} {'scale':>6} {'score':>7} {'H':>6} {'L':>6} {'fid':>5} {'RDFmax':>8}")
    for krep, krdf, tag, sc, score, H, L, fid, rmax in out:
        Hs = f"{H:.2f}" if H == H else "  - "
        Ls = f"{L:.2f}" if L == L else "  - "
        print(f"{krep:5.0f} {krdf:6.0f} {tag:4.0f} {sc:6.1f} {score:7.3f} {Hs:>6} {Ls:>6} {fid:5.2f} {rmax:8.1f}")
