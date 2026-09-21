"""Formal 51-state three-bit extension with the frozen 34-state prefix.

Layout:
  0:34   exact TwoBit34Model state order
  34:45  bit2 eleven-state biochemical module
  45:51  A1/F1 mature proteins and their expression precursors

The lower subsystem is evaluated by TwoBit34Model.rhs itself, making the
no-feedback prefix property structural rather than merely observational.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.integrate import solve_ivp

from model import BIT_NAMES, Extension, Model, ZENG, act, rep
from model_twobit34 import (STATE_NAMES_34, CarryExpressionParameters,
                            TwoBit34Model)


STATE_NAMES_51 = (STATE_NAMES_34 + tuple(f'b2_{name}' for name in BIT_NAMES) +
                  ('A1', 'F1', 'M_A1', 'A1_u', 'M_F1', 'F1_u'))
N_STATE_51 = len(STATE_NAMES_51)


def selected_threebit_extension():
    """Frozen two-bit profile plus the code-confirmed bit2 clock gate."""
    return Extension(uM_per_au=5.75, mrna_half_life_min=2.0,
                     maturation_half_life_min=20.0, translation_h=30.0,
                     immature_decay_h=0.0, complex_on_au_inv_h=0.1,
                     complex_off_h=1.0, complex_decay_h=1.0,
                     clock_K_au=0.4, clock_n=3.0, add_growth=False)


@dataclass(frozen=True)
class ThreeBitCarryParameters:
    carry0: CarryExpressionParameters = field(default_factory=lambda:
        CarryExpressionParameters(mrna_half_life_min=2.0,
                                  activator_maturation_half_life_min=32.5,
                                  repressor_maturation_half_life_min=32.5))
    carry1: CarryExpressionParameters = field(default_factory=lambda:
        CarryExpressionParameters(mrna_half_life_min=2.0,
                                  activator_maturation_half_life_min=32.5,
                                  repressor_maturation_half_life_min=32.5))


class ThreeBit51Model:
    state_names = STATE_NAMES_51
    n_state = N_STATE_51

    def __init__(self, extension=None, carry=None, n_A1_gate=None):
        self.e = extension or selected_threebit_extension()
        self.carry = carry or ThreeBitCarryParameters()
        # NEW, explicitly labelled interface parameter: the Hill exponent of the
        # A1 arm INSIDE the carry-2 gate only.
        #   None  -> use the frozen table value ZENG['n_A'][1] (previous behaviour)
        #   value -> use it for g1, while the F1 production Hill keeps the frozen
        #            exponent, so the two uses are decoupled by construction.
        # Neither ZENG nor model.py is modified; only bit2's M_I derivative is
        # corrected, so states 0..33 are untouched.
        self.n_A1_gate = n_A1_gate
        self.twobit = TwoBit34Model(self.e, self.carry.carry0)
        self.base = Model(self.e)
        self.rates_A1 = self._rates(self.carry.carry1.mrna_half_life_min,
                                    self.carry.carry1.activator_maturation_half_life_min)
        self.rates_F1 = self._rates(self.carry.carry1.mrna_half_life_min,
                                    self.carry.carry1.repressor_maturation_half_life_min)

    def _rates(self, mrna_half_life_min, maturation_half_life_min):
        lm = np.log(2)*60/mrna_half_life_min + self.base.growth
        kmat = np.log(2)*60/maturation_half_life_min
        lu = kmat + self.e.immature_decay_h + self.base.growth
        return lm, kmat, lu

    def _expression(self, m, u, protein, source, gamma, rates):
        lm, kmat, lu = rates
        tx = source*lm*lu/(self.e.translation_h*kmat)
        return np.asarray((tx-lm*m, self.e.translation_h*m-lu*u,
                           kmat*u-(gamma+self.base.growth)*protein))

    def _precursor_ss(self, source, rates):
        _, kmat, lu = rates
        return source*lu/(self.e.translation_h*kmat), source/kmat

    @property
    def n_A1_gate_effective(self):
        """Frozen table value unless an explicit gate exponent was supplied."""
        return ZENG['n_A'][1] if self.n_A1_gate is None else float(self.n_A1_gate)

    @staticmethod
    def _carry1_sources(S1, A1):
        auto = rep(A1, ZENG['K_auto1'], ZENG['n_auto1'])
        source_A1 = ZENG['alpha_A'][1]*(1-S1)*auto
        source_F1 = ZENG['alpha_F'][1]*act(A1, ZENG['K_A'][1], ZENG['n_A'][1])
        return source_A1, source_F1

    @staticmethod
    def _embed_base43(y):
        """Convert prefix layout to model.py's osc+3bits+A0/F0/A1/F1 layout."""
        y = np.asarray(y)
        if y.ndim == 1:
            out = np.zeros(43)
            out[:28] = y[:28]
            out[28:39] = y[34:45]
            out[39:41] = y[28:30]
            out[41:43] = y[45:47]
            return out
        out = np.zeros((43,y.shape[1]))
        out[:28] = y[:28]
        out[28:39] = y[34:45]
        out[39:41] = y[28:30]
        out[41:43] = y[45:47]
        return out

    def initial_state(self, cold=False):
        y34 = self.twobit.initial_state(cold=cold)
        y43 = self.base.initial_state(cold=cold)
        bit2 = y43[28:39]
        A1, F1 = y43[41], y43[42]
        if cold:
            return np.r_[y34,bit2,A1,F1,np.zeros(4)]
        source_A1, source_F1 = self._carry1_sources(y34[27],A1)
        mA1,uA1 = self._precursor_ss(source_A1,self.rates_A1)
        mF1,uF1 = self._precursor_ss(source_F1,self.rates_F1)
        return np.r_[y34,bit2,A1,F1,mA1,uA1,mF1,uF1]

    def rhs(self,t,y):
        y=np.asarray(y)
        if y.shape!=(N_STATE_51,):
            raise ValueError(f'y must have shape ({N_STATE_51},)')
        d34=self.twobit.rhs(t,y[:34])
        y43=self._embed_base43(y)
        d43=self.base.rhs(t,y43)
        S1,A1,F1=y[27],y[45],y[46]
        mA1,uA1,mF1,uF1=y[47:51]
        source_A1,source_F1=self._carry1_sources(S1,A1)
        eA1=self._expression(mA1,uA1,A1,source_A1,ZENG['gamma_A'][1],self.rates_A1)
        eF1=self._expression(mF1,uF1,F1,source_F1,ZENG['gamma_F'][1],self.rates_F1)
        d = np.r_[d34,d43[28:39],eA1[2],eF1[2],eA1[0],eA1[1],eF1[0],eF1[1]]
        n_gate = self.n_A1_gate_effective
        if n_gate != ZENG['n_A'][1]:
            # Correct ONLY bit2's M_I derivative: the carry-2 Int source enters
            # the transcript equation of M_I and nothing else.  bit2's M_I is the
            # first entry of d43[28:39], i.e. index 34 of the 51-state layout.
            _, g1_frozen, clock = self.base.carry_promoters(y43)
            g1_gate = (act(y43[41], ZENG['K_A'][1], n_gate) *
                       rep(y43[42], ZENG['K_F'][1], ZENG['n_F'][1]) * clock)
            scale = (ZENG['alpha_Int'][1] * self.base.lm * self.base.lu /
                     (self.e.translation_h * self.base.kmat))
            d[34] += float(g1_gate - g1_frozen) * scale
        return d

    def simulate(self,hours=600.0,cold=False,sample_min=2.0,
                 rtol=2e-7,atol=2e-9,max_step_min=2.0,initial_state=None):
        if hours<=0 or sample_min<=0 or max_step_min<=0:
            raise ValueError('positive duration and time steps required')
        y0=self.initial_state(cold=cold) if initial_state is None else np.asarray(initial_state,dtype=float)
        if y0.shape!=(N_STATE_51,):
            raise ValueError(f'initial_state must have shape ({N_STATE_51},)')
        ts=np.arange(0.0,hours+sample_min/120.0,sample_min/60.0)
        sol=solve_ivp(self.rhs,(0,hours),y0,t_eval=ts,method='DOP853',
                      rtol=rtol,atol=atol,max_step=max_step_min/60.0)
        if not sol.success:raise RuntimeError(sol.message)
        return sol

    @staticmethod
    def project_twobit34(y51):
        return np.asarray(y51)[:34]

    @staticmethod
    def project_twobit34_derivative(d51):
        return np.asarray(d51)[:34]

    def flux(self,y):
        return self.twobit.flux(np.asarray(y)[:34])

    def diagnostic_signals(self,y):
        y=np.asarray(y)
        if y.ndim==1:y=y[:,None]
        def hill(x,K,n):
            x=np.maximum(np.asarray(x),0.0)
            return x**n/(K**n+x**n)
        out=self.twobit.diagnostic_signals(y[:34])
        I1,R1,C1,S1=y[19],y[25],y[26],y[27]
        vf1=(ZENG['k_fwd']*hill(I1,ZENG['K_D_int'][1],2)*
             ZENG['K_inh']/(ZENG['K_inh']+np.maximum(R1,0.0)))
        vr1=ZENG['k_rev']*hill(C1,self.base.K_complex,2)
        I2,R2,C2,S2=y[36],y[42],y[43],y[44]
        vf2=(ZENG['k_fwd']*hill(I2,ZENG['K_D_int'][2],2)*
             ZENG['K_inh']/(ZENG['K_inh']+np.maximum(R2,0.0)))
        vr2=ZENG['k_rev']*hill(C2,self.base.K_complex,2)
        A1,F1=y[45],y[46]
        clock=hill(y[8],self.e.clock_K_au,self.e.clock_n)
        g1=(hill(A1,ZENG['K_A'][1],self.n_A1_gate_effective)*
            (1-hill(F1,ZENG['K_F'][1],ZENG['n_F'][1]))*clock)
        out.update(J_fwd1=vf1*(1-S1),J_rev1=vr1*S1,
                   S2=S2,J_fwd2=vf2*(1-S2),J_rev2=vr2*S2,
                   g1=g1,clock_gate=clock,Int2_source=ZENG['alpha_Int'][1]*g1)
        return out
