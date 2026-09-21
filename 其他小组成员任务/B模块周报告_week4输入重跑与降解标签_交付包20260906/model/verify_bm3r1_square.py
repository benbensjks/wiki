"""verify_bm3r1_square.py — BM3R1 delay circuit with the paper's square pulses.

Control experiment: replace TetR by BM3R1 (K=18.6 nM, n=3.4), keep the paper's
12-min pulses every 24 h. Expect the same clean toggling as the TetR version
(LR alternating ~0.3 / ~0.9).
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.integrate import solve_ivp
import zhao_core as zc

for name, K, n in [('TetR (paper)', 0.01, 2.0), ('BM3R1', 0.0186, 3.4)]:
    P = zc.default_params()
    P.update(K_rep=K, n_rep=n)
    C = zc._rate_constants(P)
    src = lambda t: P['k_int'] * zc.square_pulse(t, P)
    t_eval = np.arange(0, 24 * 5, 0.02)
    sol = solve_ivp(lambda t, y: zc.rhs(t, y, P, C, src), (0, 24 * 5),
                    zc.y0_PB_ss(P), method='LSODA', rtol=1e-6, atol=1e-11,
                    max_step=0.05, t_eval=t_eval)
    LRf = zc.LR_total(sol.y.T) / P['Dtot']
    samples = []
    for k in range(1, 6):
        i = np.argmin(np.abs(sol.t - (k * 24 + 12)))  # mid inter-pulse
        samples.append(LRf[i])
    print(f"{name}: LR at mid-gap of cycles 1-5: "
          + ' '.join(f'{v:.2f}' for v in samples))
