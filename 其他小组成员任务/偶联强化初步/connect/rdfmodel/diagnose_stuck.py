"""diagnose_stuck.py — why does the counter get stuck at LR with the oscillator input?"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import zhao_core as zc
from couple_oscillator import make_int_source

P = zc.default_params()
P.update(K_rep=0.0186, n_rep=3.4, krep_tsl=0.3)
source, spl = make_int_source(P['k_dil'])
C = zc._rate_constants(P)
t_eval = np.arange(0, 30, 0.005)
sol = solve_ivp(lambda t, y: zc.rhs(t, y, P, C, source), (0, 30),
                zc.y0_PB_ss(P), method='LSODA', rtol=1e-7, atol=1e-12,
                max_step=0.02, t_eval=t_eval)
Y = sol.y.T
LRf = zc.LR_total(Y) / P['Dtot']
INT = zc.int_total(Y)
RDF = zc.rdf_total(Y)
REP = Y[:, 37]
RDFm = Y[:, 32]

fig, ax = plt.subplots(4, 1, figsize=(12, 11), sharex=True)
ax[0].plot(sol.t, spl(sol.t), 'm'); ax[0].set_ylabel('C31 input µM')
ax[1].plot(sol.t, INT, 'g'); ax[1].set_ylabel('Int total µM')
ax[2].plot(sol.t, LRf, 'r'); ax[2].set_ylabel('LR frac'); ax[2].set_ylim(-0.05, 1.05)
ax[3].plot(sol.t, RDF, 'k', label='RDF protein')
ax[3].plot(sol.t, RDFm, 'k--', label='RDF mRNA')
ax[3].plot(sol.t, REP, 'c', label='BM3R1')
ax[3].legend(); ax[3].set_ylabel('µM'); ax[3].set_xlabel('h')
for a in ax:
    a.grid(alpha=0.3)
fig.tight_layout()
fig.savefig('../figures/diagnose_stuck.png', dpi=140)

# print samples each hour
for h in range(0, 30):
    i = np.argmin(np.abs(sol.t - h))
    print(f"t={h:2d}h  C31={spl(h):.2f}  Int={INT[i]:.3f}  LR={LRf[i]:.3f}  "
          f"RDF={RDF[i]:.3f}  RDFm={RDFm[i]:.3f}  BM3={REP[i]:.4f}")
