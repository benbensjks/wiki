"""fig_verification.py — Fig 1: port verification & BM3R1 swap with paper's pulses.
(a) Python port of the TetR model, 12-min pulses/24 h (cf. paper Fig. S11A)
(b) BM3R1 delay circuit, same pulses
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import zhao_core as zc

fig, axes = plt.subplots(2, 2, figsize=(13, 7), sharex=True)
for col, (name, K, n) in enumerate([('TetR (Zhao et al. 2019)', 0.01, 2.0),
                                    ('BM3R1 (TEMPO)', 0.0186, 3.4)]):
    P = zc.default_params()
    P.update(K_rep=K, n_rep=n)
    C = zc._rate_constants(P)
    src = lambda t: P['k_int'] * zc.square_pulse(t, P)
    t_eval = np.arange(0, 24 * 4, 0.02)
    sol = solve_ivp(lambda t, y: zc.rhs(t, y, P, C, src), (0, 24 * 4),
                    zc.y0_PB_ss(P), method='LSODA', rtol=1e-7, atol=1e-12,
                    max_step=0.05, t_eval=t_eval)
    Y = sol.y.T
    LRf = zc.LR_total(Y) / P['Dtot']
    PBf = zc.PB_total(Y) / P['Dtot']
    INT = zc.int_total(Y)
    RDF = zc.rdf_total(Y)
    REP = Y[:, 37]
    puls = np.array([zc.square_pulse(t, P) for t in sol.t])
    ax = axes[0, col]
    ax.plot(sol.t, puls, 'm:', lw=1, label='arabinose pulse')
    ax.plot(sol.t, LRf, 'r', lw=1.2, label='LR/Dtot')
    ax.plot(sol.t, PBf, 'b', lw=1.2, label='PB/Dtot')
    ax.set_ylim(-0.05, 1.1); ax.legend(fontsize=8, loc='right')
    ax.set_title(f'{name}: 12-min pulses every 24 h')
    ax = axes[1, col]
    ax.plot(sol.t, INT, 'g', lw=1, label='Int')
    ax.plot(sol.t, RDF, 'k', lw=1, label='RDF')
    ax.plot(sol.t, 5 * REP, 'c', lw=1, label=f'5×{name.split()[0]}')
    ax.legend(fontsize=8); ax.set_xlabel('time (h)'); ax.set_ylabel('µM')
axes[0, 0].set_ylabel('fraction')
fig.tight_layout()
fig.savefig('../figures/fig1_verification.png', dpi=150)
print('saved figures/fig1_verification.png')
