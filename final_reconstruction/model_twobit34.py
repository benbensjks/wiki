"""Formal 34-state two-bit biochemical counter.

State contract
--------------
6 Han-Yaxuan repressilator/C31 states
+ 11 original biochemical states for bit0
+ 11 original biochemical states for bit1
+ mature A0/F0 carry regulators
+ M_A0/A0_u/M_F0/F0_u expression states
= 34 continuous ODE states.

ZENG is imported read-only from model.py.  The nominal profile interprets the
published gamma values as total effective clearance (add_growth=False), so
growth dilution is not added a second time downstream.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from model import (BIT_NAMES, OSC_NAMES, TWOBIT_STATE_NAMES, Extension, Model,
                   ZENG, act)


STATE_NAMES_34 = (TWOBIT_STATE_NAMES +
                  ('M_A0', 'A0_u', 'M_F0', 'F0_u'))
N_STATE_34 = len(STATE_NAMES_34)


@dataclass(frozen=True)
class CarryExpressionParameters:
    """New, uncalibrated expression parameters outside Zeng's table."""

    mrna_half_life_min: float = 2.0
    activator_maturation_half_life_min: float = 30.0
    repressor_maturation_half_life_min: float = 30.0

    def __post_init__(self):
        for name in ('mrna_half_life_min', 'activator_maturation_half_life_min',
                     'repressor_maturation_half_life_min'):
            value = getattr(self, name)
            if not np.isfinite(value) or value <= 0:
                raise ValueError(f'{name} must be positive and finite')


def nominal_extension():
    """Selected two-bit profile; ZENG values themselves remain untouched."""
    return Extension(uM_per_au=6.0,
                     mrna_half_life_min=2.0,
                     maturation_half_life_min=20.0,
                     translation_h=30.0,
                     immature_decay_h=0.0,
                     complex_on_au_inv_h=0.1,
                     complex_off_h=1.0,
                     complex_decay_h=1.0,
                     clock_K_au=0.3,
                     clock_n=2.0,
                     add_growth=False)


class TwoBit34Model:
    """Production two-bit model matching the validated delayed-carry variant."""

    state_names = STATE_NAMES_34
    n_state = N_STATE_34

    def __init__(self, extension=None, carry=None):
        self.e = extension or nominal_extension()
        self.carry = carry or CarryExpressionParameters()
        self.base = Model(self.e)
        self.rates_A = self._rates(self.carry.activator_maturation_half_life_min)
        self.rates_F = self._rates(self.carry.repressor_maturation_half_life_min)

    def _rates(self, maturation_half_life_min):
        lm = np.log(2) * 60 / self.carry.mrna_half_life_min + self.base.growth
        kmat = np.log(2) * 60 / maturation_half_life_min
        lu = kmat + self.e.immature_decay_h + self.base.growth
        return lm, kmat, lu

    def _expression(self, m, u, protein, target_mature_source, gamma, rates):
        lm, kmat, lu = rates
        transcript_source = (target_mature_source * lm * lu /
                             (self.e.translation_h * kmat))
        return np.asarray((transcript_source - lm * m,
                           self.e.translation_h * m - lu * u,
                           kmat * u - (gamma + self.base.growth) * protein))

    def _precursor_ss(self, source, rates):
        _, kmat, lu = rates
        m = source * lu / (self.e.translation_h * kmat)
        u = source / kmat
        return m, u

    def initial_state(self, cold=False):
        y30 = self.base.initial_state_twobit(cold=cold)
        if cold:
            return np.r_[y30, np.zeros(4)]
        A0 = y30[28]
        source_A = ZENG['alpha_A'][0] * (1 - y30[16])
        source_F = ZENG['alpha_F'][0] * act(A0, ZENG['K_A'][0], ZENG['n_A'][0])
        mA, uA = self._precursor_ss(source_A, self.rates_A)
        mF, uF = self._precursor_ss(source_F, self.rates_F)
        return np.r_[y30, mA, uA, mF, uF]

    def rhs(self, t, y):
        y = np.asarray(y)
        if y.shape != (N_STATE_34,):
            raise ValueError(f'y must have shape ({N_STATE_34},)')
        y30 = y[:30]
        d30 = self.base.rhs_twobit(t, y30)
        A0, F0 = y30[28], y30[29]
        mA, uA, mF, uF = y[30:34]
        source_A = ZENG['alpha_A'][0] * (1 - y30[16])
        source_F = ZENG['alpha_F'][0] * act(A0, ZENG['K_A'][0], ZENG['n_A'][0])
        eA = self._expression(mA, uA, A0, source_A, ZENG['gamma_A'][0], self.rates_A)
        eF = self._expression(mF, uF, F0, source_F, ZENG['gamma_F'][0], self.rates_F)
        d30[28], d30[29] = eA[2], eF[2]
        return np.r_[d30, eA[:2], eF[:2]]

    def simulate(self, hours=300.0, cold=False, sample_min=2.0,
                 rtol=2e-7, atol=2e-9, max_step_min=2.0, initial_state=None):
        if hours <= 0 or sample_min <= 0 or max_step_min <= 0:
            raise ValueError('positive duration and time steps required')
        y0 = self.initial_state(cold=cold) if initial_state is None else np.asarray(initial_state, dtype=float)
        if y0.shape != (N_STATE_34,):
            raise ValueError(f'initial_state must have shape ({N_STATE_34},)')
        ts = np.arange(0.0, hours + sample_min / 120.0, sample_min / 60.0)
        sol = solve_ivp(self.rhs, (0, hours), y0, t_eval=ts, method='DOP853',
                        rtol=rtol, atol=atol, max_step=max_step_min / 60.0)
        if not sol.success:
            raise RuntimeError(sol.message)
        return sol

    def flux(self, y):
        return self.base.flux(self.base._expand_twobit(np.asarray(y)[:30]))

    def diagnostic_signals(self, y):
        """Return physical signals used by the causal acceptance test."""
        y = np.asarray(y)
        if y.ndim == 1:
            y = y[:, None]
        S0, S1 = y[16], y[27]
        I0, R0, C0 = y[8], y[14], y[15]
        A0, F0 = y[28], y[29]

        def hill(x, K, n):
            x = np.maximum(np.asarray(x), 0.0)
            return x**n / (K**n + x**n)

        vf0 = (ZENG['k_fwd'] * hill(I0, ZENG['K_D_int'][0], 2) *
               ZENG['K_inh'] / (ZENG['K_inh'] + np.maximum(R0, 0.0)))
        vr0 = ZENG['k_rev'] * hill(C0, self.base.K_complex, 2)
        j_fwd0 = vf0 * (1 - S0)
        j_rev0 = vr0 * S0
        g0 = hill(A0, ZENG['K_A'][0], ZENG['n_A'][0]) * \
            (1 - hill(F0, ZENG['K_F'][0], ZENG['n_F'][0]))
        return dict(S0=S0, S1=S1, J_fwd0=j_fwd0, J_rev0=j_rev0,
                    g0=g0, Int1_source=ZENG['alpha_Int'][0] * g0)
