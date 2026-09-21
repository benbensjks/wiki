"""verify_tetr.py — verify the Python port against Zhao 2019 published behaviour.

Test 1 (Fig. 5E/F): single 12-min arabinose pulse at 0.5 h, TetR delay circuit.
  - PB start: efficient PB->LR within ~1 h
  - LR start: efficient LR->PB
Test 2 (Fig. S11A): repeated 12-min pulses every 24 h over 5 days.
  - counter toggles: LR fraction alternates between ~30 % (low) and ~88 % (high)
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from scipy.integrate import solve_ivp
import zhao_core as zc

P = zc.default_params()
C = zc._rate_constants(P)

def int_src(t):
    return P['k_int'] * zc.square_pulse(t, P)

def run(y0, t_span, label):
    sol = solve_ivp(lambda t, y: zc.rhs(t, y, P, C, int_src),
                    t_span, y0, method='LSODA',
                    rtol=1e-6, atol=1e-11, max_step=0.02,
                    dense_output=False, t_eval=np.arange(t_span[0], t_span[1], 0.01))
    Y = sol.y.T
    LR = zc.LR_total(Y) / P['Dtot']
    PB = zc.PB_total(Y) / P['Dtot']
    print(f"--- {label} ---")
    print(f"  initial LR frac = {LR[0]:.3f}, final LR frac = {LR[-1]:.3f}")
    return sol, LR, PB

# Test 1a: PB start, single pulse (Fig 5E)
sol1, LR1, _ = run(zc.y0_PB(P), (0, 2), "Fig5E: PB start, 12-min pulse at 0.5 h")
i_pulse_end = np.argmin(np.abs(sol1.t - 1.0))
print(f"  LR frac at t=1.0 h (30 min after pulse end): {LR1[i_pulse_end]:.3f}  (paper: near-complete)")

# Test 1b: LR start, single pulse (Fig 5F)
y0lr = np.zeros(38); y0lr[0] = P['Dtot']; y0lr[22] = 0.5
sol2, LR2, _ = run(y0lr, (0, 2), "Fig5F: LR start, 12-min pulse at 0.5 h")

# Test 2: repeated pulsing 5 days (Fig S11A)
sol3, LR3, _ = run(zc.y0_PB_ss(P), (0, 24 * 5), "FigS11A: 12-min pulses every 24 h, 5 days")
# sample LR just before each pulse (t = k*24 + 0.4) and after switching
for k in range(1, 6):
    t_b, t_a = k * 24 - 0.1, k * 24 + 2
    i_b, i_a = np.argmin(np.abs(sol3.t - t_b)), np.argmin(np.abs(sol3.t - t_a))
    print(f"  cycle {k}: LR before pulse {LR3[i_b]:.2f}, after pulse {LR3[i_a]:.2f}")

# mass conservation check
dna = zc.LR_total(sol3.y.T) + zc.PB_total(sol3.y.T)
print(f"  DNA conservation: min={dna.min():.5f}, max={dna.max():.5f} (Dtot={P['Dtot']})")
