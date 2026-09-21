"""leak_rate.py — how low must the Int baseline be?

Constant Int input at baseline levels; measure PB->LR flip time (no RDF) and
LR->PB flip time (RDF at steady state) as functions of the baseline Int
concentration. The inter-pulse gap of the v36 oscillator is ~4.8 h, so the
memory requires flip time >> 5 h in BOTH directions.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import zhao_core as zc

LEVELS = [0.005, 0.01, 0.02, 0.03, 0.05, 0.08, 0.12, 0.17, 0.28]
T_MAX = 12.0


def flip_time(P, y0, level, direction):
    C = zc._rate_constants(P)
    src = lambda t: level * P['k_dil'] if t >= 0 else 0.0  # sustains [int]=level
    # production = k_dil*level keeps free int at 'level' absent consumption
    t_eval = np.arange(0, T_MAX, 0.02)
    sol = solve_ivp(lambda t, y: zc.rhs(t, y, P, C, src), (0, T_MAX), y0,
                    method='LSODA', rtol=1e-7, atol=1e-12, max_step=0.05,
                    t_eval=t_eval)
    LRf = zc.LR_total(sol.y.T) / P['Dtot']
    target = 0.5 if direction == 'up' else 0.5
    if direction == 'up':
        idx = np.where(LRf >= target)[0]
    else:
        idx = np.where(LRf <= 1 - target)[0]
    return sol.t[idx[0]] if len(idx) else np.inf


P = zc.default_params()
P.update(K_rep=0.0186, n_rep=3.4, krep_tsl=10.0, krdf_tsl=8.0)

print(f"{'Int baseline (µM)':>18} {'copies':>7} {'PB->LR t1/2 (h)':>16} {'LR->PB t1/2 (h)':>16}")
rows = []
for lv in LEVELS:
    # PB start, repressor at PB steady state, no RDF
    y0_up = zc.y0_PB(P, rep_mrna=0.51, rep=P['krep_tsl'] * 0.51 / P['k_dil'])
    t_up = flip_time(P, y0_up, lv, 'up')
    # LR start, RDF at LR steady state (repressor ~0)
    y0_dn = zc.y0_LR_ss(P)
    t_dn = flip_time(P, y0_dn, lv, 'down')
    rows.append((lv, lv * 602, t_up, t_dn))
    fmt = lambda x: f"{x:15.2f}" if np.isfinite(x) else f"{' >12':>15}"
    print(f"{lv:18.3f} {lv*602:7.0f} {fmt(t_up)} {fmt(t_dn)}")

rows = np.array(rows)
fig, ax = plt.subplots(figsize=(7.5, 5))
ax.plot(rows[:, 1], np.minimum(rows[:, 2], T_MAX), 'o-', label='PB→LR (RDF=0)')
ax.plot(rows[:, 1], np.minimum(rows[:, 3], T_MAX), 's-', label='LR→PB (RDF=SS)')
ax.axhline(4.8, color='r', ls='--', label='inter-pulse gap ≈ 4.8 h')
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlabel('Int baseline (copies/cell)'); ax.set_ylabel('flip half-time (h)')
ax.set_title('State-holding constraint on the Int baseline\n(krep_tsl=10, krdf_tsl=8, K=18.6 nM)')
ax.legend(); ax.grid(alpha=0.3, which='both')
fig.tight_layout()
fig.savefig('../figures/leak_rate.png', dpi=150)
print('saved figures/leak_rate.png')
