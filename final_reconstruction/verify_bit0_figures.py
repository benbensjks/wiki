"""Figures for the bit0 report: operating point comparison and the carry gate."""
from __future__ import annotations

import json
from dataclasses import asdict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from model import Model, Extension, ZENG, ROOT
from verify_bit0_part2 import OUT, CAND, clock_cycles

H = lambda x, K, n: np.maximum(x, 0.0)**n / (K**n + np.maximum(x, 0.0)**n)


def load(cfg, hours=80.0):
    m = Model(Extension(**cfg))
    sol = m.simulate_twobit(hours=hours, sample_min=2.0, max_step_min=2.0)
    flux = np.array([m.flux(m._expand_twobit(sol.y[:, k])) for k in range(sol.y.shape[1])])
    return m, sol, flux


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    base = asdict(Extension()); base.update(CAND)
    runs = {}
    for u in (6.0, 10.0):
        cfg = base.copy(); cfg['uM_per_au'] = u
        runs[u] = load(cfg)
    m, sol, flux = runs[6.0]
    t = sol.t
    S0, S1 = sol.y[16], sol.y[27]
    A0, F0 = sol.y[28], sol.y[29]
    g0 = H(A0, ZENG['K_A'][0], ZENG['n_A'][0]) * (1 - H(F0, ZENG['K_F'][0], ZENG['n_F'][0]))
    ids = clock_cycles(t, flux)

    fig, ax = plt.subplots(5, 1, figsize=(13, 14), sharex=True, constrained_layout=True)
    for u, ls in ((6.0, '-'), (10.0, '--')):
        s = runs[u][1]
        ax[0].plot(s.t, s.y[16], ls, label=f'bit0 S, uM_per_au={u:g}')
    ax[0].axhspan(0.3, 0.7, color='grey', alpha=.15, label='ambiguous')
    ax[0].set_ylabel('S0'); ax[0].legend(fontsize=8)
    ax[0].set_title('bit0 dwells in committed bands only at the lower a.u.->uM conversion')
    ax[1].plot(t, S0, label='bit0 S'); ax[1].plot(t, S1, label='bit1 S')
    ax[1].set_ylabel('DNA state, uM_per_au=6'); ax[1].legend(fontsize=8)
    ax[2].plot(t, A0, label='A0 (carry activator)'); ax[2].plot(t, F0, label='F0 (carry repressor)')
    ax[2].axhline(ZENG['K_F'][0], color='k', ls=':', label=f"F0 threshold {ZENG['K_F'][0]}")
    ax[2].set_ylabel('a.u.'); ax[2].legend(fontsize=8)
    ax[3].plot(t, g0, c='#D55E00', label='carry promoter g0 = H(A0)*G(F0)')
    u10 = runs[10.0]
    g0b = H(u10[1].y[28], ZENG['K_A'][0], ZENG['n_A'][0]) * (1 - H(u10[1].y[29], ZENG['K_F'][0], ZENG['n_F'][0]))
    ax[3].plot(u10[1].t, g0b, c='#0072B2', ls='--', label='same gate at uM_per_au=10')
    ax[3].set_ylabel('promoter activity'); ax[3].legend(fontsize=8)
    for a in ax[:4]:
        for i in ids:
            a.axvline(t[i], color='grey', lw=.4, alpha=.4)
        a.grid(alpha=.2)
    ax[4].plot(t, np.gradient(S0, t), label='dS0/dt'); ax[4].plot(t, np.gradient(S1, t), label='dS1/dt')
    ax[4].set_xlabel('time (h)'); ax[4].set_ylabel('h$^{-1}$'); ax[4].legend(fontsize=8); ax[4].grid(alpha=.2)
    fig.savefig(OUT / 'bit0_operating_point.png', dpi=160); plt.close(fig)

    pd.DataFrame({'time_h': t, 'bit0_S': S0, 'bit1_S': S1, 'A0': A0, 'F0': F0, 'carry_g0': g0,
                  'flux_uM_h': flux}).to_csv(OUT / 'twobit_uM6_trajectory.csv', index=False)
    print('wrote', OUT / 'bit0_operating_point.png')


if __name__ == '__main__':
    main()
