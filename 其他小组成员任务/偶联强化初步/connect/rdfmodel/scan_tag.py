"""scan_tag.py — ssrA-tagged integrase rescue scan.

k_tag_int: extra Int degradation (h^-1). Source amplitude is compensated by
factor (k_dil+k_tag)/k_dil so the Int peak stays ~4 µM (stronger promoter
compensates the tag — an experimental design lever).
baseline_sub: residual Int leak removal (promoter tightness).
Fixed: krep_tsl=10, krdf_tsl=8, K=18.6 nM, n=3.4.
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
TAGS = [0.0, 2.0, 4.0, 8.0, 12.0, 16.0]
BSUBS = [0.0, 0.12, 0.20, 0.26]
KREP = float(os.environ.get('KREP_TSL', '15'))
KRDF = float(os.environ.get('KRDF_TSL', '200'))


def one(job):
    tag, b = job
    P = zc.default_params()
    P.update(K_rep=0.0186, n_rep=3.4, krep_tsl=KREP, krdf_tsl=KRDF,
             k_tag_int=tag)
    # NO peak compensation: the tag itself suppresses the baseline by
    # (k_dil+k_tag)/k_dil; the peak stays >= ~0.8 µM, well above threshold.
    src = make_source_baseline(P['k_dil'], baseline_sub=b, scale=1.0)
    y0 = zc.y0_PB(P, rep_mrna=P['k_tscr'] * P['Dtot'] / P['k_rna'],
                  rep=KREP * P['k_tscr'] * P['Dtot'] / P['k_rna'] / P['k_dil'])
    sol = simulate(P, src, y0)
    LRf = zc.LR_total(sol.y.T) / P['Dtot']
    score, H, L, fid, _ = toggle_score(sol.t, LRf, TROUGHS)
    return tag, b, score, H, L, fid


if __name__ == '__main__':
    jobs = [(t, b) for t in TAGS for b in BSUBS]
    with Pool(12) as p:
        out = p.map(one, jobs)
    Z = np.zeros((len(TAGS), len(BSUBS)))
    for tag, b, score, H, L, fid in out:
        i, j = TAGS.index(tag), BSUBS.index(b)
        Z[i, j] = score
        print(f"k_tag={tag:5.1f}  b_sub={b:.2f}  score={score:.3f}  "
              f"H={H if H==H else -1:.2f} L={L if L==L else -1:.2f} fid={fid:.2f}")
    np.savez('scan_tag.npz', tags=TAGS, bsubs=BSUBS, Z=Z)
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    im = ax.imshow(Z, origin='lower', aspect='auto', cmap='viridis',
                   extent=[0, len(BSUBS) - 1, 0, len(TAGS) - 1])
    ax.set_xticks(range(len(BSUBS)))
    ax.set_xticklabels([f'-{b}\n(resid {max(0,0.28-b)*602:.0f} cp)' for b in BSUBS],
                       fontsize=8)
    ax.set_yticks(range(len(TAGS)))
    ax.set_yticklabels([f'{t:.0f} (t½ {np.log(2)/(2+t)*60:.0f} min)' for t in TAGS],
                       fontsize=8)
    ax.set_xlabel('Int baseline removed (µM) / residual copies')
    ax.set_ylabel('k_tag_int (h⁻¹) / Int half-life')
    ax.set_title(f'Tagged-integrase rescue: toggle score\n(krep_tsl={KREP}, krdf_tsl={KRDF})')
    for i in range(len(TAGS)):
        for j in range(len(BSUBS)):
            ax.text(j, i, f'{Z[i, j]:.2f}', ha='center', va='center',
                    color='w' if Z[i, j] < 0.55 else 'k', fontsize=9)
    fig.colorbar(im, label='toggle score')
    fig.tight_layout()
    fig.savefig('../figures/scan_tag.png', dpi=150)
    print('saved figures/scan_tag.png')
