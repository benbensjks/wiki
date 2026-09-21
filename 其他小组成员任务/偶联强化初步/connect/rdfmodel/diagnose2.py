"""diagnose2.py — detailed dynamics at krep_tsl=10, krdf_tsl=8, baseline_sub=0.12."""
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

KREP = float(os.environ.get('KREP_TSL', '10'))
KRDF = float(os.environ.get('KRDF_TSL', '8'))
BSUB = float(os.environ.get('BASELINE_SUB', '0.12'))

P = zc.default_params()
P.update(K_rep=0.0186, n_rep=3.4, krep_tsl=KREP, krdf_tsl=KRDF)
C = zc._rate_constants(P)
src = make_source_baseline(P['k_dil'], baseline_sub=BSUB, scale=1.0)
spl, _ = load_c31_spline(1.0)

t_eval = np.arange(0, 40, 0.01)
sol = solve_ivp(lambda t, y: zc.rhs(t, y, P, C, src), (0, 40),
                zc.y0_PB_ss(P), method='LSODA', rtol=1e-7, atol=1e-12,
                max_step=0.05, t_eval=t_eval)
Y = sol.y.T
LRf = zc.LR_total(Y) / P['Dtot']
INT = zc.int_total(Y)
RDF = zc.rdf_total(Y)
REP = Y[:, 37]

fig, ax = plt.subplots(4, 1, figsize=(13, 11), sharex=True)
ax[0].plot(sol.t, np.maximum(0, spl(sol.t) - BSUB), 'm'); ax[0].set_ylabel('Int input µM')
ax[1].plot(sol.t, INT, 'g'); ax[1].set_ylabel('Int total µM')
ax[2].plot(sol.t, LRf, 'r'); ax[2].set_ylabel('LR frac'); ax[2].set_ylim(-0.05, 1.05)
ax[3].plot(sol.t, RDF, 'k', label='RDF')
ax[3].plot(sol.t, REP, 'c', label='BM3R1')
ax[3].legend(); ax[3].set_ylabel('µM'); ax[3].set_xlabel('h')
for a in ax:
    a.grid(alpha=0.3)
fig.suptitle(f'krep_tsl={KREP}, krdf_tsl={KRDF}, baseline_sub={BSUB}')
fig.tight_layout()
fig.savefig('../figures/diagnose2.png', dpi=140)

for h in np.arange(0, 40, 1.0):
    i = np.argmin(np.abs(sol.t - h))
    print(f"t={h:4.1f}  Int={INT[i]:6.3f}  LR={LRf[i]:.3f}  RDF={RDF[i]:6.3f}  BM3={REP[i]:7.4f}")
