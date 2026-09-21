"""scan_delay.py — BM3R1/RDF expression-rate scan (the delay-circuit design space).

krep_tsl: BM3R1 translation rate (h^-1)  — proxy for Pc promoter x RBS driving bm3r1
krdf_tsl: RDF translation rate (h^-1)    — proxy for RBS driving rdf
Waveform: A-module v36 with Int baseline removed (baseline_sub) — required
promoter tightness identified in scan_waveform.
Metric: toggle_score (clean H/L alternation at troughs).
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

KREP = [0.3, 1.0, 3.0, 6.0, 10.0, 15.0, 30.0]
KRDF = [1.0, 2.0, 4.0, 8.0, 16.0, 40.0]
BASELINE_SUB = float(os.environ.get('BASELINE_SUB', '0.12'))
SCALE = float(os.environ.get('SCALE', '1.0'))


def one(job):
    krep, krdf = job
    P = zc.default_params()
    P.update(K_rep=0.0186, n_rep=3.4, krep_tsl=krep, krdf_tsl=krdf)
    src = make_source_baseline(P['k_dil'], baseline_sub=BASELINE_SUB, scale=SCALE)
    sol = simulate(P, src, zc.y0_PB_ss(P))
    LRf = zc.LR_total(sol.y.T) / P['Dtot']
    score, H, L, fid, _ = toggle_score(sol.t, LRf, TROUGHS)
    return krep, krdf, score, H, L, fid


if __name__ == '__main__':
    jobs = [(a, b) for a in KREP for b in KRDF]
    with Pool(12) as p:
        out = p.map(one, jobs)
    Z = np.zeros((len(KREP), len(KRDF)))
    for krep, krdf, score, H, L, fid in out:
        i, j = KREP.index(krep), KRDF.index(krdf)
        Z[i, j] = score
        print(f"krep={krep:5.1f}  krdf={krdf:5.1f}  score={score:.3f}  "
              f"H={H if H==H else float('nan'):.2f} L={L if L==L else float('nan'):.2f} fid={fid:.2f}")
    np.savez(f'scan_delay_b{BASELINE_SUB}.npz', krep=KREP, krdf=KRDF, Z=Z,
             baseline_sub=BASELINE_SUB, scale=SCALE)

    fig, ax = plt.subplots(figsize=(8, 5.5))
    im = ax.imshow(Z, origin='lower', aspect='auto', cmap='viridis',
                   extent=[np.log10(min(KRDF)), np.log10(max(KRDF)),
                           np.log10(min(KREP)), np.log10(max(KREP))])
    ax.set_xlabel('log10 krdf_tsl (h$^{-1}$)')
    ax.set_ylabel('log10 k_BM3R1_tsl (h$^{-1}$)')
    ax.set_title(f'Toggle score: BM3R1 vs RDF expression rates\n'
                 f'(Int baseline -{BASELINE_SUB} µM, scale {SCALE}, K=18.6 nM, n=3.4)')
    for i, a in enumerate(KREP):
        for j, b in enumerate(KRDF):
            ax.text(np.log10(b), np.log10(a), f'{Z[i, j]:.2f}', ha='center',
                    va='center', color='w' if Z[i, j] < 0.55 else 'k', fontsize=8)
    fig.colorbar(im, label='toggle score')
    fig.tight_layout()
    fn = f'../figures/scan_delay_b{BASELINE_SUB}.png'
    fig.savefig(fn, dpi=150)
    print('saved', fn)
