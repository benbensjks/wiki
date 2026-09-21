"""diagnose3.py — timecourse with tagged Int, watch LR→PB failure mechanism."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import zhao_core as zc
from analysis_lib import make_source_baseline
from couple_oscillator import load_c31_spline

TAG = float(os.environ.get('K_TAG', '8'))
BSUB = float(os.environ.get('BASELINE_SUB', '0.12'))
KREP = float(os.environ.get('KREP_TSL', '10'))
KRDF = float(os.environ.get('KRDF_TSL', '8'))

P = zc.default_params()
P.update(K_rep=0.0186, n_rep=3.4, krep_tsl=KREP, krdf_tsl=KRDF, k_tag_int=TAG)
C = zc._rate_constants(P)
scale = (P['k_dil'] + TAG) / P['k_dil']
src = make_source_baseline(P['k_dil'], baseline_sub=BSUB, scale=scale)
spl, _ = load_c31_spline(1.0)
y0 = zc.y0_PB(P, rep_mrna=P['k_tscr'] * P['Dtot'] / P['k_rna'],
              rep=KREP * P['k_tscr'] * P['Dtot'] / P['k_rna'] / P['k_dil'])
t_eval = np.arange(0, 30, 0.01)
sol = solve_ivp(lambda t, y: zc.rhs(t, y, P, C, src), (0, 30), y0,
                method='LSODA', rtol=1e-7, atol=1e-12, max_step=0.05,
                t_eval=t_eval)
Y = sol.y.T
LRf = zc.LR_total(Y) / P['Dtot']
INT = zc.int_total(Y)
RDF = zc.rdf_total(Y)
REP = Y[:, 37]

fig, ax = plt.subplots(4, 1, figsize=(13, 11), sharex=True)
ax[0].plot(sol.t, INT, 'g'); ax[0].set_ylabel('Int total µM')
ax[1].plot(sol.t, LRf, 'r'); ax[1].set_ylabel('LR frac'); ax[1].set_ylim(-0.05, 1.05)
ax[2].plot(sol.t, RDF, 'k'); ax[2].set_ylabel('RDF µM')
ax[3].plot(sol.t, REP, 'c'); ax[3].set_ylabel('BM3R1 µM'); ax[3].set_xlabel('h')
for a in ax:
    a.grid(alpha=0.3)
fig.suptitle(f'k_tag={TAG}, b_sub={BSUB}, krep={KREP}, krdf={KRDF}')
fig.tight_layout()
fig.savefig('../figures/diagnose3.png', dpi=140)

for h in np.arange(0, 30, 0.5):
    i = np.argmin(np.abs(sol.t - h))
    print(f"t={h:5.1f}  Int={INT[i]:6.3f}  LR={LRf[i]:.3f}  RDF={RDF[i]:6.3f}  BM3={REP[i]:7.4f}")
