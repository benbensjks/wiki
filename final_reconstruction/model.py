"""43-state continuous reconstruction. Time: h; downstream concentrations: a.u.

Six oscillator states + three original eleven-state modules + A0,F0,A1,F1.
F0/F1 are Zeng's R0/R1 carry repressors, renamed to distinguish them from RDF.
Bit0 M_I IS the upstream C31 mRNA, not a second transcript pool.
The original upstream source is imported read-only and reused at each RHS call.
"""
from __future__ import annotations

import importlib.util
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parent
WIKI = ROOT.parent
UPSTREAM_PATH = WIKI / '其他小组成员任务/Week4_振荡器-C31建模/code/Flux_Driven_Translation_Burden_Model.py'
REFERENCE_PATH = UPSTREAM_PATH.parent.parent / 'data/reference_inputs/unloaded_C31_translation_trajectory.csv'
spec = importlib.util.spec_from_file_location('han_upstream_v53d', UPSTREAM_PATH)
upstream = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = upstream
spec.loader.exec_module(upstream)

BIT_NAMES = ('M_I', 'I_u', 'I', 'M_T', 'T_u', 'T', 'M_R', 'R_u', 'R', 'C', 'S')
OSC_NAMES = ('m_TetR', 'TetR_total', 'm_CI', 'CI', 'm_LacI', 'LacI')
STATE_NAMES = OSC_NAMES + tuple(f'b{i}_{s}' for i in range(3) for s in BIT_NAMES) + ('A0', 'F0', 'A1', 'F1')
N_STATE = len(STATE_NAMES)
BIT0_STATE_NAMES = OSC_NAMES + tuple(f'b0_{s}' for s in BIT_NAMES)
TWOBIT_STATE_NAMES = (OSC_NAMES + tuple(f'b{i}_{s}' for i in range(2) for s in BIT_NAMES)
                      + ('A0', 'F0'))

# Exact published table values. a.u. is not silently relabelled as micromolar.
ZENG = dict(k_fwd=7.0, k_rev=5.0, K_D_int=(1.0, 1.0, 0.6),
            K_D_comp=1.2, K_inh=0.1, k_int=6.0,
            gamma_int=(2.0, 1.4, 2.2), alpha_rep=3.5, gamma_rep=0.6,
            K_rep=0.85, n_rep=3.9, alpha_rdf=6.0, gamma_rdf=0.8,
            alpha_A=(8.0, 16.0), gamma_A=(1.9, 2.5), K_A=(0.5, 1.2), n_A=(2.0, 4.0),
            alpha_F=(3.8, 5.0), gamma_F=(0.6, 1.1), K_F=(0.6, 0.4), n_F=(4.0, 4.0),
            alpha_Int=(18.0, 38.0), K_auto1=0.6, n_auto1=2.0)


@dataclass(frozen=True)
class Extension:
    # New parameters NOT identified by Zeng's reduced model; explicit assumptions.
    uM_per_au: float = 1.0
    mrna_half_life_min: float = 2.0  # intrinsic downstream mRNA half-life
    maturation_half_life_min: float = 5.0
    translation_h: float = 30.0
    immature_decay_h: float = 0.0
    complex_on_au_inv_h: float = 1.0
    complex_off_h: float = 10.0
    complex_decay_h: float = 1.0
    clock_K_au: float = 0.3  # no numerical K_gate/n_gate found in supplied Zeng table
    clock_n: float = 2.0
    # True interprets table gamma as intrinsic degradation, as table labels state.
    add_growth: bool = True
    def __post_init__(self):
        for k in ('uM_per_au', 'mrna_half_life_min', 'maturation_half_life_min',
                  'translation_h', 'complex_on_au_inv_h', 'complex_off_h', 'clock_K_au', 'clock_n'):
            if not np.isfinite(getattr(self, k)) or getattr(self, k) <= 0:
                raise ValueError(f'{k} must be positive and finite')
        for k in ('immature_decay_h', 'complex_decay_h'):
            if not np.isfinite(getattr(self, k)) or getattr(self, k) < 0:
                raise ValueError(f'{k} must be nonnegative and finite')


def act(x, K, n):
    z = max(float(x), 0.0) / K
    return z**n / (1.0 + z**n)


def rep(x, K, n):
    return 1.0 - act(x, K, n)


def bit_slice(i):
    return slice(6 + 11*i, 6 + 11*(i+1))


class Model:
    def __init__(self, extension=None):
        self.e = extension or Extension()
        self.a = upstream.make_reference_parameters()
        self.mu = 60*self.a.mu
        self.growth = self.mu if self.e.add_growth else 0.0
        self.copies_per_au = 602.214076*self.e.uM_per_au  # original 1 fL volume
        self.lm = np.log(2)*60/self.e.mrna_half_life_min + self.growth
        self.kmat = np.log(2)*60/self.e.maturation_half_life_min
        self.lu = self.kmat + self.e.immature_decay_h + self.growth
        # At QSS C = q I R; preserve Zeng's product-response half-saturation.
        # K_D_comp in her formula acts on I*R (thus a.u.^2 dimensionally).
        self.q = self.e.complex_on_au_inv_h/(self.e.complex_off_h+self.e.complex_decay_h+self.growth)
        self.K_complex = self.q*ZENG['K_D_comp']

    def transcript_source(self, target_mature_source):
        """QSS matching: kmat * Iu_ss = target_mature_source (a.u./h)."""
        return target_mature_source*self.lm*self.lu/(self.e.translation_h*self.kmat)

    def expression(self, m, u, protein, source, gamma):
        return np.array((self.transcript_source(source)-self.lm*m,
                         self.e.translation_h*m-self.lu*u,
                         self.kmat*u-(gamma+self.growth)*protein))

    def upstream_state(self, y):
        # The unused reference C31 protein is NOT a dynamic duplicate state.
        return np.r_[y[:6], y[6]*self.copies_per_au, 0.0]

    def flux(self, y):
        # h^-1 conversion; rho=1 at frozen unloaded baseline.
        return 60*self.a.c31_translation_per_mrna_per_min*y[6]*self.e.uM_per_au

    def carry_promoters(self, y):
        A0, F0, A1, F1 = y[39:43]
        g0 = act(A0, ZENG['K_A'][0], ZENG['n_A'][0])*rep(F0, ZENG['K_F'][0], ZENG['n_F'][0])
        clock = act(y[8], self.e.clock_K_au, self.e.clock_n)  # actual mature free Int0
        g1 = act(A1, ZENG['K_A'][1], ZENG['n_A'][1])*rep(F1, ZENG['K_F'][1], ZENG['n_F'][1])*clock
        return g0, g1, clock

    def rhs(self, t, y):
        # No clipping/projection, events, forced toggles or state reset in dynamics.
        d = np.zeros(N_STATE)
        da = 60*np.asarray(upstream.rhs(t*60, self.upstream_state(y), self.a, 1.0))
        d[:6] = da[:6]
        g0, g1, _ = self.carry_promoters(y)
        sources = (None, ZENG['alpha_Int'][0]*g0, ZENG['alpha_Int'][1]*g1)
        for i in range(3):
            b = y[bit_slice(i)]
            mI, Iu, I, mT, Tu, T, mR, Ru, R, C, S = b
            db = np.zeros(11)
            if i == 0:
                db[0] = da[6]/self.copies_per_au
                db[1] = self.flux(y)/self.e.uM_per_au-self.lu*Iu
                db[2] = self.kmat*Iu-(ZENG['gamma_int'][0]+self.growth)*I
            else:
                db[:3] = self.expression(mI, Iu, I, sources[i], ZENG['gamma_int'][i])
            # T is Zeng's PB-driven Rep/BM3R1, not oscillator TetR.
            db[3:6] = self.expression(mT, Tu, T, ZENG['alpha_rep']*(1-S), ZENG['gamma_rep'])
            db[6:9] = self.expression(mR, Ru, R, ZENG['alpha_rdf']*S*rep(T, ZENG['K_rep'], ZENG['n_rep']), ZENG['gamma_rdf'])
            binding = self.e.complex_on_au_inv_h*I*R
            unbinding = self.e.complex_off_h*C
            db[2] += -binding+unbinding
            db[8] += -binding+unbinding
            db[9] = binding-unbinding-(self.e.complex_decay_h+self.growth)*C
            # Explicit free-I depletion implements sequestration; the empirical
            # K_inh factor is retained as an additional Zeng inhibitory effect.
            vf = ZENG['k_fwd']*act(I, ZENG['K_D_int'][i], 2)*ZENG['K_inh']/(ZENG['K_inh']+max(R, 0.0))
            vr = ZENG['k_rev']*act(C, self.K_complex, 2)
            db[10] = vf*(1-S)-vr*S
            d[bit_slice(i)] = db
        for j in range(2):
            A, F = y[39+2*j:41+2*j]
            PB = 1-y[bit_slice(j)][10]
            auto = rep(A, ZENG['K_auto1'], ZENG['n_auto1']) if j == 1 else 1.0
            d[39+2*j] = ZENG['alpha_A'][j]*PB*auto-(ZENG['gamma_A'][j]+self.growth)*A
            d[40+2*j] = ZENG['alpha_F'][j]*act(A, ZENG['K_A'][j], ZENG['n_A'][j])-(ZENG['gamma_F'][j]+self.growth)*F
        return d

    def initial_state(self, cold=False):
        y = np.zeros(N_STATE)
        a0 = upstream.initial_state(self.a)
        y[:6] = a0[:6]
        y[6] = a0[6]/self.copies_per_au
        if not cold:
            # DNA initially PB; pre-equilibrated repressors prevent artificial
            # carry pulses caused solely by starting every regulator at zero.
            for i in range(3):
                b = y[bit_slice(i)]
                b[3] = self.transcript_source(ZENG['alpha_rep'])/self.lm
                b[4] = self.e.translation_h*b[3]/self.lu
                b[5] = ZENG['alpha_rep']/(ZENG['gamma_rep']+self.growth)
            for j in range(2):
                loss = ZENG['gamma_A'][j]+self.growth
                A = ZENG['alpha_A'][j]/loss if j == 0 else brentq(
                    lambda a: ZENG['alpha_A'][j]*rep(a, ZENG['K_auto1'], ZENG['n_auto1'])-loss*a,
                    0, ZENG['alpha_A'][j]/loss)
                y[39+2*j] = A
                y[40+2*j] = ZENG['alpha_F'][j]*act(A, ZENG['K_A'][j], ZENG['n_A'][j])/(ZENG['gamma_F'][j]+self.growth)
        return y

    def simulate(self, hours=150.0, cold=False, strict=False):
        if not np.isfinite(hours) or hours <= 0:
            raise ValueError('hours must be positive')
        ts = np.linspace(0, hours, int(np.ceil(hours*60))+1)
        sol = solve_ivp(self.rhs, (0, hours), self.initial_state(cold), t_eval=ts,
                        method='DOP853', rtol=2e-10 if strict else 2e-8,
                        atol=2e-12 if strict else 2e-10, max_step=1/120 if strict else 1/60)
        if not sol.success:
            raise RuntimeError(sol.message)
        return sol

    def _expand_bit0(self, y17):
        """Embed the genuinely dynamic 17-state subsystem in the full layout."""
        y = self.initial_state()
        y[:17] = y17
        return y

    def rhs_bit0(self, t, y17):
        return self.rhs(t, self._expand_bit0(y17))[:17]

    def initial_state_bit0(self, cold=False):
        return self.initial_state(cold)[:17].copy()

    def simulate_bit0(self, hours=100.0, cold=False, sample_min=1.0,
                      rtol=2e-7, atol=2e-9, max_step_min=2.0):
        if hours <= 0 or sample_min <= 0 or max_step_min <= 0:
            raise ValueError('positive duration and time steps required')
        ts = np.arange(0.0, hours + sample_min/120.0, sample_min/60.0)
        sol = solve_ivp(self.rhs_bit0, (0, hours), self.initial_state_bit0(cold),
                        t_eval=ts, method='DOP853', rtol=rtol, atol=atol,
                        max_step=max_step_min/60.0)
        if not sol.success:
            raise RuntimeError(sol.message)
        return sol

    def _expand_twobit(self, y30):
        y = self.initial_state()
        y[:28] = y30[:28]
        y[39:41] = y30[28:30]
        return y

    def rhs_twobit(self, t, y30):
        d = self.rhs(t, self._expand_twobit(y30))
        return np.r_[d[:28], d[39:41]]

    def initial_state_twobit(self, cold=False):
        y = self.initial_state(cold)
        return np.r_[y[:28], y[39:41]]

    def simulate_twobit(self, hours=100.0, cold=False, sample_min=1.0,
                        rtol=2e-7, atol=2e-9, max_step_min=2.0):
        ts = np.arange(0.0, hours + sample_min/120.0, sample_min/60.0)
        sol = solve_ivp(self.rhs_twobit, (0, hours), self.initial_state_twobit(cold),
                        t_eval=ts, method='DOP853', rtol=rtol, atol=atol,
                        max_step=max_step_min/60.0)
        if not sol.success:
            raise RuntimeError(sol.message)
        return sol
