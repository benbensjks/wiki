"""scan_waveform.py — how tight/strong must the Int-driving promoter be?

Scan A-module waveform modifications:
  baseline_sub: subtract a constant from C31(t) [µM] — proxy for reducing the
                PLtetO1 leak / Int baseline (promoter tightness)
  scale       : multiply whole waveform — proxy for promoter/RBS strength
Metric: toggle_score at C31 troughs (clean alternation H/L after transient).
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

BASELINES = [0.0, 0.02, 0.05, 0.08, 0.12, 0.16, 0.20, 0.24]
SCALES = [0.5, 0.75, 1.0, 1.5, 2.0]


def one(job):
    b, s = job
    P = zc.default_params()
    P.update(K_rep=0.0186, n_rep=3.4, krep_tsl=0.3)
    src = make_source_baseline(P['k_dil'], baseline_sub=b, scale=s)
    sol = simulate(P, src, zc.y0_PB_ss(P))
    LRf = zc.LR_total(sol.y.T) / P['Dtot']
    score, H, L, fid, _ = toggle_score(sol.t, LRf, TROUGHS)
    return b, s, score, H, L, fid


if __name__ == '__main__':
    jobs = [(b, s) for b in BASELINES for s in SCALES]
    with Pool(12) as p:
        out = p.map(one, jobs)
    Z = np.zeros((len(BASELINES), len(SCALES)))
    Hz = np.full_like(Z, np.nan); Lz = np.full_like(Z, np.nan)
    for b, s, score, H, L, fid in out:
        i, j = BASELINES.index(b), SCALES.index(s)
        Z[i, j] = score
        Hz[i, j] = H; Lz[i, j] = L
        print(f"baseline_sub={b:.2f}  scale={s:.2f}  score={score:.3f}  H={H:.2f} L={L:.2f}")
    np.savez('scan_waveform.npz', baselines=BASELINES, scales=SCALES, Z=Z, H=Hz, L=Lz)

    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    im = ax.imshow(Z, origin='lower', aspect='auto', cmap='viridis',
                   extent=[min(SCALES), max(SCALES),
                           min(BASELINES) * 602, max(BASELINES) * 602])
    ax.set_xlabel('waveform amplitude scale (×)')
    ax.set_ylabel('Int baseline removed (copies/cell)')
    ax.set_title('Toggle score vs Int waveform properties\n(BM3R1 delay, K=18.6 nM, n=3.4)')
    for i, b in enumerate(BASELINES):
        for j, s in enumerate(SCALES):
            ax.text(s, b * 602, f'{Z[i, j]:.2f}', ha='center', va='center',
                    color='w' if Z[i, j] < 0.5 else 'k', fontsize=8)
    fig.colorbar(im, label='toggle score')
    fig.tight_layout()
    fig.savefig('../figures/scan_waveform.png', dpi=150)
    print('saved figures/scan_waveform.png')
