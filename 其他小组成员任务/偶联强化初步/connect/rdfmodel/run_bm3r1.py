"""run_bm3r1.py — BM3R1 delay-circuit counter driven by the oscillator waveform.

Compares:
  (a) TetR delay circuit (paper's repressor, control)  K=0.01 µM, n=2
  (b) BM3R1 delay circuit (TEMPO design)               K=0.0186 µM, n=3.4
      K_BM3 derived so that steady-state repression leak of the RDF promoter
      equals the Cello B1_BM3R1 measured leak (ymin/ymax = 0.004/0.5 = 0.8 %):
      (S/K)^3.4 = 1/0.008 - 1 with S = 0.0765 µM (PB-state steady state).
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import zhao_core as zc
from couple_oscillator import make_int_source, pulse_times, COPIES_PER_UM

OUT = os.path.join(_HERE := os.path.dirname(os.path.abspath(__file__)), '..', 'figures')
os.makedirs(OUT, exist_ok=True)

T_END = 100.0  # h (full CSV range)


def simulate(P, source, y0, t_end=T_END, t_eval=None):
    C = zc._rate_constants(P)
    if t_eval is None:
        t_eval = np.arange(0, t_end, 0.02)
    sol = solve_ivp(lambda t, y: zc.rhs(t, y, P, C, source),
                    (0, t_end), y0, method='LSODA',
                    rtol=1e-6, atol=1e-11, max_step=0.05, t_eval=t_eval)
    return sol


def toggle_metrics(t, LRf, peaks_h, troughs_h, skip=2):
    """Per-pulse directional efficiencies, using troughs as sampling points."""
    tr = troughs_h
    if len(tr) < 4:
        return np.nan, np.nan, []
    LR_at = np.interp(tr, t, LRf)
    fwd, rev = [], []
    for k in range(skip, len(tr) - 1):
        before, after = LR_at[k], LR_at[k + 1]
        if before < 0.5:
            fwd.append(after)          # PB->LR: want after high
        else:
            rev.append(1 - after)      # LR->PB: want after low
    f = np.mean(fwd) if fwd else np.nan
    r = np.mean(rev) if rev else np.nan
    return f, r, list(zip(tr[skip:-1], LR_at[skip:-1]))


def run_case(name, P, y0):
    source, spl = make_int_source(P['k_dil'])
    sol = simulate(P, source, y0)
    Y = sol.y.T
    LRf = zc.LR_total(Y) / P['Dtot']
    PBf = zc.PB_total(Y) / P['Dtot']
    INT = zc.int_total(Y)
    RDF = zc.rdf_total(Y)
    REP = Y[:, 37]
    pk, tr = pulse_times()
    f, r, samples = toggle_metrics(sol.t, LRf, pk, tr)
    print(f"[{name}]  fwd(PB->LR) eff = {f:.3f}, rev(LR->PB) eff = {r:.3f}, "
          f"min = {np.nanmin([f, r]):.3f}")
    print(f"    LR frac at troughs (after transient): "
          + ' '.join(f'{v:.2f}' for _, v in samples[:12]))
    return dict(name=name, t=sol.t, LRf=LRf, PBf=PBf, INT=INT, RDF=RDF,
                REP=REP, spl=spl, eff=(f, r))


P_tet = zc.default_params()                                  # TetR control
P_bm3 = zc.default_params()
P_bm3.update(K_rep=0.0186, n_rep=3.4, krep_tsl=0.3)          # BM3R1 delay

res_tet = run_case('TetR control', P_tet, zc.y0_PB_ss(P_tet))
res_bm3 = run_case('BM3R1', P_bm3, zc.y0_PB_ss(P_bm3))

# ---- figure
fig, axes = plt.subplots(3, 2, figsize=(13, 10), sharex=True)
for col, res in enumerate([res_tet, res_bm3]):
    t = res['t']
    ax = axes[0, col]
    ax.plot(t, res['spl'](t), 'm:', lw=1, label='C31 input (µM)')
    ax.plot(t, res['INT'], 'g', lw=1, label='Int total')
    ax.set_title(f"{res['name']}: Int input")
    ax.legend(fontsize=8); ax.set_ylabel('µM')
    ax = axes[1, col]
    ax.plot(t, res['LRf'], 'r', lw=1.2, label='LR/Dtot')
    ax.plot(t, res['PBf'], 'b', lw=1.2, label='PB/Dtot')
    ax.set_title(f"{res['name']}: DNA state   (fwd={res['eff'][0]:.2f}, rev={res['eff'][1]:.2f})")
    ax.legend(fontsize=8); ax.set_ylabel('fraction')
    ax.set_ylim(-0.05, 1.05)
    ax = axes[2, col]
    ax.plot(t, res['RDF'], 'k', lw=1, label='RDF total')
    ax.plot(t, res['REP'], 'c', lw=1, label=f"{res['name']} protein")
    ax.set_title(f"{res['name']}: RDF & repressor")
    ax.legend(fontsize=8); ax.set_ylabel('µM'); ax.set_xlabel('time (h)')
fig.suptitle('Single-input counter driven by repressilator φC31 waveform (v36)')
fig.tight_layout()
fig.savefig(os.path.join(OUT, 'fig_bm3r1_vs_tetr_oscillator.png'), dpi=150)
print('saved', os.path.join(OUT, 'fig_bm3r1_vs_tetr_oscillator.png'))
