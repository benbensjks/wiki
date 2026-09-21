"""scan_refine.py — refine around the discovered working regime + robustness to leak."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from multiprocessing import Pool
import zhao_core as zc
from analysis_lib import simulate, make_source_baseline, toggle_score
from couple_oscillator import pulse_times

_, TROUGHS = pulse_times()

KREPS = [15.0, 30.0, 60.0]
KRDFS = [100.0, 200.0, 400.0]
TAGS = [6.0, 8.0, 12.0]
BSUBS = [0.20, 0.24, 0.26]  # residual 0.08 / 0.04 / 0.02 µM


def one(job):
    krep, krdf, tag, b = job
    P = zc.default_params()
    P.update(K_rep=0.0186, n_rep=3.4, krep_tsl=krep, krdf_tsl=krdf,
             k_tag_int=tag)
    src = make_source_baseline(P['k_dil'], baseline_sub=b, scale=1.0)
    y0 = zc.y0_PB(P, rep_mrna=P['k_tscr'] * P['Dtot'] / P['k_rna'],
                  rep=krep * P['k_tscr'] * P['Dtot'] / P['k_rna'] / P['k_dil'])
    sol = simulate(P, src, y0)
    LRf = zc.LR_total(sol.y.T) / P['Dtot']
    score, H, L, fid, _ = toggle_score(sol.t, LRf, TROUGHS)
    return krep, krdf, tag, b, score, H, L, fid


if __name__ == '__main__':
    jobs = [(a, bb, c, d) for a in KREPS for bb in KRDFS for c in TAGS for d in BSUBS]
    with Pool(12) as p:
        out = p.map(one, jobs)
    out.sort(key=lambda r: -r[4])
    print(f"{'krep':>5} {'krdf':>6} {'tag':>4} {'b_sub':>6} {'score':>7} {'H':>6} {'L':>6} {'fid':>5}")
    for krep, krdf, tag, b, score, H, L, fid in out[:30]:
        Hs = f"{H:.2f}" if H == H else "  - "
        Ls = f"{L:.2f}" if L == L else "  - "
        print(f"{krep:5.0f} {krdf:6.0f} {tag:4.0f} {b:6.2f} {score:7.3f} {Hs:>6} {Ls:>6} {fid:5.2f}")
    import numpy as np
    np.savez('scan_refine.npz', out=np.array([[a, bb, c, d, s] for a, bb, c, d, s, H, L, f in out]))
