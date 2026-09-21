"""fig_baseline_req.py — toggle score vs residual Int baseline (operating point)."""
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
BS = [0.0, 0.08, 0.12, 0.16, 0.20, 0.22, 0.24, 0.25, 0.26, 0.27]


def one(b):
    P = zc.default_params()
    P.update(K_rep=0.0186, n_rep=3.4, krep_tsl=15.0, krdf_tsl=200.0, k_tag_int=12.0)
    src = make_source_baseline(P['k_dil'], baseline_sub=b, scale=1.0)
    y0 = zc.y0_PB(P, rep_mrna=P['k_tscr'] * P['Dtot'] / P['k_rna'],
                  rep=15.0 * P['k_tscr'] * P['Dtot'] / P['k_rna'] / P['k_dil'])
    sol = simulate(P, src, y0)
    LRf = zc.LR_total(sol.y.T) / P['Dtot']
    s, H, L, f, _ = toggle_score(sol.t, LRf, TROUGHS)
    return b, s


if __name__ == '__main__':
    with Pool(10) as p:
        out = p.map(one, BS)
    out.sort()
    resid = [max(0, 0.28 - b) * 602 for b, s in out]
    scores = [s for b, s in out]
    for b, s in out:
        print(f"b_sub={b:.2f} residual={max(0,0.28-b)*602:.0f} copies score={s:.3f}")
    fig, ax = plt.subplots(figsize=(7.5, 5))
    ax.plot(resid, scores, 'o-', color='tab:blue')
    ax.axvline(170, color='r', ls='--', label="A-module v36 baseline (170 copies)")
    ax.axvspan(0, 50, alpha=0.15, color='g', label='working region')
    ax.set_xlabel('residual Int baseline (copies/cell)')
    ax.set_ylabel('toggle score')
    ax.set_title('Promoter-tightness requirement at the operating point\n'
                 '(krep_tsl=15, krdf_tsl=200, k_tag=12)')
    ax.legend(); ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig('../figures/fig_baseline_requirement.png', dpi=150)
    print('saved figures/fig_baseline_requirement.png')
