"""scan_K.py — sensitivity to BM3R1 repression threshold K and Hill n.

At a chosen operating point (krep_tsl, krdf_tsl, baseline_sub), scan
K_rep in [5, 50] nM and n_rep in {2.9, 3.4} (Cello B2/B1 BM3R1 fits).
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from multiprocessing import Pool
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import zhao_core as zc
from analysis_lib import simulate, make_source_baseline, toggle_score
from couple_oscillator import pulse_times

_, TROUGHS = pulse_times()

KREP = float(os.environ.get('KREP_TSL', '15'))
KRDF = float(os.environ.get('KRDF_TSL', '200'))
BSUB = float(os.environ.get('BASELINE_SUB', '0.26'))
TAG = float(os.environ.get('K_TAG', '12'))
Ks = [0.005, 0.008, 0.012, 0.0186, 0.025, 0.035, 0.05]
NS = [2.9, 3.4]


def one(job):
    K, n = job
    P = zc.default_params()
    P.update(K_rep=K, n_rep=n, krep_tsl=KREP, krdf_tsl=KRDF, k_tag_int=TAG)
    src = make_source_baseline(P['k_dil'], baseline_sub=BSUB, scale=1.0)
    y0 = zc.y0_PB(P, rep_mrna=P['k_tscr'] * P['Dtot'] / P['k_rna'],
                  rep=KREP * P['k_tscr'] * P['Dtot'] / P['k_rna'] / P['k_dil'])
    sol = simulate(P, src, y0)
    LRf = zc.LR_total(sol.y.T) / P['Dtot']
    score, H, L, fid, _ = toggle_score(sol.t, LRf, TROUGHS)
    return K, n, score, H, L


if __name__ == '__main__':
    jobs = [(K, n) for K in Ks for n in NS]
    with Pool(12) as p:
        out = p.map(one, jobs)
    for K, n, score, H, L in out:
        print(f"K={K*1000:5.1f} nM  n={n}  score={score:.3f}  H={H:.2f} L={L:.2f}")
    fig, ax = plt.subplots(figsize=(7, 5))
    for n in NS:
        xs = [K * 1000 for K, nn, s, H, L in out if nn == n]
        ys = [s for K, nn, s, H, L in out if nn == n]
        ax.plot(xs, ys, 'o-', label=f'n = {n}')
    ax.set_xlabel('K_BM3R1 (nM)'); ax.set_ylabel('toggle score')
    ax.set_title(f'K sensitivity at krep_tsl={KREP}, krdf_tsl={KRDF}, baseline -{BSUB} µM')
    ax.legend(); ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig('../figures/scan_K.png', dpi=150)
    print('saved figures/scan_K.png')
