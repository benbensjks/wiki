"""analysis_lib.py — shared simulation & metric utilities for the single-stage scans."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.integrate import solve_ivp
import zhao_core as zc
from couple_oscillator import load_c31_spline, pulse_times

T_END = 100.0


def simulate(P, source, y0, t_end=T_END, dt=0.05):
    C = zc._rate_constants(P)
    t_eval = np.arange(0, t_end, dt)
    sol = solve_ivp(lambda t, y: zc.rhs(t, y, P, C, source),
                    (0, t_end), y0, method='LSODA',
                    rtol=1e-6, atol=1e-11, max_step=0.1, t_eval=t_eval)
    return sol


def make_source_baseline(k_dil, baseline_sub=0.0, scale=1.0):
    """Int source from A's waveform: w(t) = scale · max(0, C31(t) - baseline_sub)
    (baseline subtraction = tighter promoter; scale = promoter strength).
    Production S(t) = dw/dt + k_dil·w, floored at 0.  Note scale amplifies the
    residual baseline too — it cannot improve the peak:leak ratio.
    """
    spl, dspl = load_c31_spline(1.0)

    def source(t):
        c = float(spl(t))
        w = max(0.0, c - baseline_sub)
        dwdt = float(dspl(t)) if c > baseline_sub else 0.0
        return scale * max(0.0, dwdt + k_dil * w)

    return source


def toggle_score(t, LRf, troughs_h, skip=3):
    """Score clean one-toggle-per-period behaviour.

    Sample LR fraction at each C31 trough (mid inter-pulse gap). A perfect
    counter alternates H/L/H/L after transient.
    Returns (score, H_mean, L_mean, alternation_fidelity, samples).
      score = min(H_mean, 1-L_mean) * alternation_fidelity
    """
    tr = np.array(troughs_h)
    tr = tr[tr <= t[-1]]
    if len(tr) < skip + 4:
        return 0.0, np.nan, np.nan, 0.0, []
    s = np.interp(tr, t, LRf)[skip:]
    hi = s > 0.6
    lo = s < 0.4
    mid = ~(hi | lo)
    # alternation fidelity: fraction of consecutive samples that strictly alternate
    alt = 0
    n = 0
    for k in range(len(s) - 1):
        if mid[k] or mid[k + 1]:
            continue
        n += 1
        if hi[k] != hi[k + 1]:
            alt += 1
    fid = alt / n if n else 0.0
    H = s[hi].mean() if hi.any() else np.nan
    L = s[lo].mean() if lo.any() else np.nan
    core = min(np.nan_to_num(H, nan=0.0), 1 - np.nan_to_num(L, nan=1.0))
    return core * fid, H, L, fid, list(zip(tr[skip:], s))
