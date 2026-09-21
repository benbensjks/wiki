"""final_run.py — operating-point simulation + deliverable numbers.

Chosen operating point from scans (see report): BM3R1 delay circuit,
A-module v36 waveform with Int baseline reduced to ≤0.02 µM (12 copies).
Outputs: full timecourse figure, toggle metrics, max expression ratio,
promoter-strength conversions for the experimental group.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import zhao_core as zc
from analysis_lib import make_source_baseline, toggle_score
from couple_oscillator import load_c31_spline, pulse_times

KREP = float(os.environ.get('KREP_TSL', '15'))
KRDF = float(os.environ.get('KRDF_TSL', '200'))
BSUB = float(os.environ.get('BASELINE_SUB', '0.26'))
TAG = float(os.environ.get('K_TAG', '12'))
K_REP = float(os.environ.get('K_REP', '0.0186'))
N_REP = float(os.environ.get('N_REP', '3.4'))

_, TROUGHS = pulse_times()
PEAKS, _ = pulse_times()

P = zc.default_params()
P.update(K_rep=K_REP, n_rep=N_REP, krep_tsl=KREP, krdf_tsl=KRDF, k_tag_int=TAG)
C = zc._rate_constants(P)
src = make_source_baseline(P['k_dil'], baseline_sub=BSUB, scale=1.0)
spl, _ = load_c31_spline(1.0)

# initial condition: PB with BM3R1 at its PB steady state
y0 = zc.y0_PB(P, rep_mrna=P['k_tscr'] * P['Dtot'] / P['k_rna'],
              rep=KREP * P['k_tscr'] * P['Dtot'] / P['k_rna'] / P['k_dil'])

t_eval = np.arange(0, 100, 0.01)
sol = solve_ivp(lambda t, y: zc.rhs(t, y, P, C, src), (0, 100), y0,
                method='LSODA', rtol=1e-7, atol=1e-12, max_step=0.05,
                t_eval=t_eval)
Y = sol.y.T
LRf = zc.LR_total(Y) / P['Dtot']
PBf = zc.PB_total(Y) / P['Dtot']
INT = zc.int_total(Y)
RDF = zc.rdf_total(Y)
REP = Y[:, 37]

score, H, L, fid, samples = toggle_score(sol.t, LRf, TROUGHS)
print(f"operating point: krep_tsl={KREP}, krdf_tsl={KRDF}, k_tag_int={TAG}, "
      f"baseline_sub={BSUB}, K={K_REP*1000:.1f} nM, n={N_REP}")
print(f"toggle score = {score:.3f}  (H={H:.3f}, L={L:.3f}, alternation fidelity={fid:.2f})")
print("trough samples:", ' '.join(f'{v:.2f}' for _, v in samples))
print(f"Int peak = {INT.max():.2f} µM, Int trough = {INT[-1]:.3f} µM")
print(f"RDF max = {RDF.max():.2f} µM ({RDF.max()*602:.0f} copies), "
      f"BM3R1 max = {REP.max():.2f} µM")

# ---- expression ratio of the first-stage switch
hi_samples = [v for _, v in samples if v > 0.6]
lo_samples = [v for _, v in samples if v < 0.4]
if hi_samples and lo_samples:
    on = np.mean(hi_samples); off = np.mean(lo_samples)
    print(f"first-stage switch expression ratio: ON={on:.3f}, OFF={off:.3f}, "
          f"ON/OFF={on/off:.1f}x, max LR fraction={max(hi_samples):.3f}")

# ---- promoter strength conversions
mrna_ss = P['k_tscr'] * P['Dtot'] / P['k_rna']
bm3_prod = KREP * mrna_ss            # µM/h
print(f"\nBM3R1 mRNA steady state (PB) = {mrna_ss:.3f} µM = {mrna_ss*602:.0f} copies")
print(f"BM3R1 production = {bm3_prod:.2f} µM/h = {bm3_prod*602:.0f} copies/h")
print(f"BM3R1 steady state = {KREP*mrna_ss/P['k_dil']:.2f} µM")
print(f"delay after flip ~ ln(SS/K)/k_dil = "
      f"{np.log(KREP*mrna_ss/P['k_dil']/K_REP)/P['k_dil']:.2f} h")
int1_ss = KREP * mrna_ss * H / P['k_dil'] if H == H else np.nan
print(f"same promoter flipped to LR drives Int1: production "
      f"{KREP*mrna_ss*H:.2f} µM/h -> Int1 SS ~ {int1_ss:.2f} µM")

# ---- figure
fig, ax = plt.subplots(4, 1, figsize=(13, 11), sharex=True)
ax[0].plot(sol.t, np.maximum(0, spl(sol.t) - BSUB), 'm', lw=1)
ax[0].set_ylabel('Int input (µM)'); ax[0].set_title('φC31 waveform (baseline-corrected)')
ax[1].plot(sol.t, INT, 'g', lw=1); ax[1].set_ylabel('Int total (µM)')
ax[2].plot(sol.t, LRf, 'r', lw=1.2, label='LR')
ax[2].plot(sol.t, PBf, 'b', lw=1.2, label='PB')
tt, ss = zip(*samples)
ax[2].plot(tt, ss, 'ko', ms=4)
ax[2].set_ylabel('DNA fraction'); ax[2].set_ylim(-0.05, 1.05); ax[2].legend()
ax[3].plot(sol.t, RDF, 'k', lw=1, label='RDF')
ax[3].plot(sol.t, REP, 'c', lw=1, label='BM3R1')
ax[3].set_ylabel('µM'); ax[3].set_xlabel('time (h)'); ax[3].legend()
for a in ax:
    a.grid(alpha=0.3)
fig.suptitle(f'BM3R1 counter @ krep_tsl={KREP}, krdf_tsl={KRDF}, baseline -{BSUB} µM, '
             f'score={score:.2f}')
fig.tight_layout()
fig.savefig('../figures/final_operating_point.png', dpi=150)
print('\nsaved figures/final_operating_point.png')
