"""Preliminary cross-model counter. Integration time: hours.

ZMH donor functions are read from a hash-pinned source, without executing its
top-level integrations/plotting. No upstream process runs at import time.
"""
from __future__ import annotations

import ast
import hashlib
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parent
WIKI = ROOT.parent
DONOR_PATH = WIKI / 'wiki任务' / '完整二级级联.py'
BASE_PATH = WIKI / 'final_reconstruction' / 'model.py'
DONOR_SHA = 'E790743F493F415E43296CD10911BB6C2A4E5A024F70E807A40267E9705A0E2B'
BASE_SHA = '256D2D104EECB4CE17F27CF4DBC08257BEEEA79E5805FD7D5614E47BBC7E543A'
DONOR_NAMES = ('pb0', 'int0', 'rep0', 'rdf0', 'a0', 'r0',
               'pb1', 'int1', 'rep1', 'rdf1')
CARRY_NAMES = ('A1', 'F1', 'M_A1', 'A1_u', 'M_F1', 'F1_u')
BIT_NAMES = ('M_I', 'I_u', 'I', 'M_T', 'T_u', 'T', 'M_R', 'R_u', 'R', 'C', 'S')
TAIL_NAMES = CARRY_NAMES + tuple('b2_' + n for n in BIT_NAMES)
STATE_NAMES = DONOR_NAMES + TAIL_NAMES
INDEX = {n: i for i, n in enumerate(STATE_NAMES)}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest().upper()


def require_hash(path, expected):
    actual = sha256(path)
    if actual != expected:
        raise RuntimeError(f'Source changed: {path}; expected {expected}, got {actual}')


def source_hashes():
    return {str(p.relative_to(WIKI)): sha256(p) for p in (DONOR_PATH, BASE_PATH)}


def load_zeng_table():
    """Read only the literal parameter table; do not import the old model."""
    require_hash(BASE_PATH, BASE_SHA)
    tree = ast.parse(BASE_PATH.read_text(encoding='utf-8-sig'))
    node = next(n for n in tree.body if isinstance(n, ast.Assign)
                and any(isinstance(t, ast.Name) and t.id == 'ZENG' for t in n.targets))
    if not isinstance(node.value, ast.Call) or not isinstance(node.value.func, ast.Name) or node.value.func.id != 'dict':
        raise RuntimeError('Unexpected ZENG table format')
    return {k.arg: ast.literal_eval(k.value) for k in node.value.keywords}


def donor_namespace():
    """Isolate source-owned globals per instance; retain exact RHS expressions."""
    require_hash(DONOR_PATH, DONOR_SHA)
    tree = ast.parse(DONOR_PATH.read_text(encoding='utf-8-sig'))
    functions = {'hill_rep', 'PLtetO', '_exact_free', 'free_TetR',
                 'repressilator_ode', 'counter_2bit_ode'}
    assignments = {'_Ttot_grid', '_Tfree_grid', 'p'}
    selected, found = [], set()
    for node in tree.body:
        name = None
        if isinstance(node, ast.ClassDef) and node.name == 'P':
            name = node.name
        elif isinstance(node, ast.FunctionDef) and node.name in functions:
            name = node.name
        elif isinstance(node, ast.Assign) and len(node.targets) == 1:
            t = node.targets[0]
            if isinstance(t, ast.Name) and t.id in assignments:
                name = t.id
        if name:
            selected.append(node)
            found.add(name)
    if found != functions | assignments | {'P'}:
        raise RuntimeError(f'Donor extraction contract mismatch: {found}')
    ns = {'np': np, 'math': math}
    exec(compile(ast.Module(body=selected, type_ignores=[]), str(DONOR_PATH), 'exec'), ns)
    return ns


class ZmhDonor:
    """Original five-state oscillator + original ten-state donor, time adapted."""
    def __init__(self, upstream_mode='original'):
        if upstream_mode not in ('original', 'tight'):
            raise ValueError('upstream_mode must be original or tight')
        self.ns = donor_namespace()
        self.p = self.ns['p'].copy()
        self.upstream_mode = upstream_mode
        options = {} if upstream_mode == 'original' else {'rtol': 1e-9, 'atol': 1e-11}
        self.upstream_options = dict(method='BDF', rtol=options.get('rtol', 1e-3),
                                     atol=options.get('atol', 1e-6))
        self.upstream = solve_ivp(self.ns['repressilator_ode'], (0, 3000),
                                  [0.5, 0.2, 1, 0, 0], t_eval=np.linspace(0, 3000, 3000),
                                  method='BDF', **options)
        if not self.upstream.success or not np.all(np.isfinite(self.upstream.y)):
            raise RuntimeError(f'Upstream failed: {self.upstream.message}')
        self.c31_max = float(np.max(self.upstream.y[4][1000:]))
        if self.c31_max <= 0:
            raise RuntimeError('Invalid normalization denominator')
        self.ns['get_u_in'] = self.input_min

    def input_min(self, t_min):
        # Reject silent extrapolation. Small endpoint tolerance is roundoff only.
        if not -1e-9 <= t_min <= 3000 + 1e-9:
            raise ValueError('Input outside the original 0..3000 minute horizon')
        val = np.interp(t_min, self.upstream.t, self.upstream.y[4])
        return float(np.clip(val / self.c31_max, 0, 2))

    def initial_state(self):
        p = self.p
        a = p['alpha_A0'] / p['gamma_A0']
        act = a**p['n_A0'] / (p['K_A0']**p['n_A0'] + a**p['n_A0'])
        r = p['alpha_R0'] * act / p['gamma_R0']
        return np.array([1., 0., p['alpha_rep0']/p['gamma_rep0'], 0., a, r,
                         1., 0., p['alpha_rep1']/p['gamma_rep1'], 0.])

    def rhs_min(self, t, y):
        return np.asarray(self.ns['counter_2bit_ode'](t, y, self.p))

    def rhs(self, t_h, y):
        return 60.0 * self.rhs_min(60.0*t_h, y)

    def gate0(self, y):
        p = self.p
        a, r = y[4], y[5]
        return (a**p['n_A0']/(p['K_A0']**p['n_A0']+a**p['n_A0']) *
                p['K_R0']**p['n_R0']/(p['K_R0']**p['n_R0']+r**p['n_R0']))

    def metadata(self):
        return dict(source_sha256=DONOR_SHA, parameters_per_minute=self.p,
                    oscillator_parameters={k: v for k, v in vars(self.ns['P']).items()
                                           if not k.startswith('_')},
                    upstream_options=self.upstream_options, c31_max=self.c31_max,
                    normalization='max(C31[1000:]) on original 3000-sample 0..3000 min grid',
                    input='normalized C31 protein concentration, NOT translation flux',
                    table_max_TetR=3000., max_TetR=float(self.upstream.y[0].max()),
                    free_TetR_table_exceeded=bool(self.upstream.y[0].max()>3000),
                    counter_time_original='min', adapter_time='h')


@dataclass(frozen=True)
class HbyConfig:
    n_A1_gate: float = 6.0
    clock_K_au: float = 0.3
    clock_n: float = 2.0
    clock_scale: float = 1.0  # new, uncalibrated donor-to-receiver interface assumption
    carry_mrna_min: float = 2.0
    carry_maturation_min: float = 32.5
    bit_mrna_min: float = 2.0
    bit_maturation_min: float = 20.0
    translation_h: float = 30.0
    kon: float = 0.1
    koff: float = 1.0
    complex_decay: float = 1.0
    receiver_uM_per_au: float = 5.75  # metadata only; does NOT calibrate donor concentrations

    def __post_init__(self):
        if any(not np.isfinite(v) or v <= 0 for v in asdict(self).values()):
            raise ValueError('All preliminary configuration values must be finite and positive')


def hill(x, k, n):
    z = max(float(x), 0.) / k
    return z**n / (1+z**n)


class HbyReceiver:
    """Six-state carry1 followed by eleven-state bit2; all rates per hour."""
    def __init__(self, config=None):
        self.c = config or HbyConfig()
        self.p = load_zeng_table()
        self.lmc = np.log(2)*60/self.c.carry_mrna_min
        self.kmc = np.log(2)*60/self.c.carry_maturation_min
        self.lmb = np.log(2)*60/self.c.bit_mrna_min
        self.kmb = np.log(2)*60/self.c.bit_maturation_min
        self.q = self.c.kon/(self.c.koff+self.c.complex_decay)
        self.k_complex = self.q*self.p['K_D_comp']

    def expression(self, m, u, protein, source, gamma, lm, km):
        # lu=km, additional growth=0, immature intrinsic loss=0.
        return np.array([source*lm/self.c.translation_h-lm*m,
                         self.c.translation_h*m-km*u, km*u-gamma*protein])

    def signals(self, y, int0):
        p, c = self.p, self.c
        a, f = y[:2]
        i, r, comp, s = y[8], y[14], y[15], y[16]
        clock = hill(c.clock_scale*int0, c.clock_K_au, c.clock_n)
        gate = hill(a, p['K_A'][1], c.n_A1_gate)*(1-hill(f,p['K_F'][1],p['n_F'][1]))*clock
        vf = p['k_fwd']*hill(i,p['K_D_int'][2],2)*p['K_inh']/(p['K_inh']+max(r,0))
        vr = p['k_rev']*hill(comp,self.k_complex,2)
        return dict(clock_gate=clock, g1=gate, u2_target_au_per_h=p['alpha_Int'][1]*gate,
                    J_fwd2_per_h=vf*(1-s), J_rev2_per_h=vr*s)

    def rhs(self, t, y, pb1, int0):
        p, c = self.p, self.c
        a, f, ma, au, mf, fu = y[:6]
        mi, iu, i, mt, tu, tr, mr, ru, r, comp, s = y[6:]
        sa = p['alpha_A'][1]*pb1*(1-hill(a,p['K_auto1'],p['n_auto1']))
        sf = p['alpha_F'][1]*hill(a,p['K_A'][1],p['n_A'][1])
        ea = self.expression(ma,au,a,sa,p['gamma_A'][1],self.lmc,self.kmc)
        ef = self.expression(mf,fu,f,sf,p['gamma_F'][1],self.lmc,self.kmc)
        sig = self.signals(y,int0)
        db = np.zeros(11)
        db[:3] = self.expression(mi,iu,i,sig['u2_target_au_per_h'],p['gamma_int'][2],self.lmb,self.kmb)
        db[3:6] = self.expression(mt,tu,tr,p['alpha_rep']*(1-s),p['gamma_rep'],self.lmb,self.kmb)
        db[6:9] = self.expression(mr,ru,r,p['alpha_rdf']*s*(1-hill(tr,p['K_rep'],p['n_rep'])),p['gamma_rdf'],self.lmb,self.kmb)
        bind, unbind = c.kon*i*r, c.koff*comp
        db[2] += -bind+unbind
        db[8] += -bind+unbind
        db[9] = bind-unbind-c.complex_decay*comp
        db[10] = sig['J_fwd2_per_h']-sig['J_rev2_per_h']
        return np.r_[ea[2],ef[2],ea[0],ea[1],ef[0],ef[1],db]

    def initial_state(self, pb1=1.):
        p = self.p
        a = brentq(lambda a: p['alpha_A'][1]*pb1*(1-hill(a,p['K_auto1'],p['n_auto1']))-p['gamma_A'][1]*a,
                   0, p['alpha_A'][1]/p['gamma_A'][1])
        sa = p['gamma_A'][1]*a
        sf = p['alpha_F'][1]*hill(a,p['K_A'][1],p['n_A'][1])
        carry = [a,sf/p['gamma_F'][1],sa/self.c.translation_h,sa/self.kmc,
                 sf/self.c.translation_h,sf/self.kmc]
        bit = np.zeros(11)
        bit[3:6] = (p['alpha_rep']/self.c.translation_h,
                    p['alpha_rep']/self.kmb,p['alpha_rep']/p['gamma_rep'])
        return np.r_[carry,bit]


class HybridModel:
    state_names = STATE_NAMES
    def __init__(self, donor, receiver=None):
        self.donor = donor
        self.receiver = receiver or HbyReceiver()

    def initial_state(self):
        d = self.donor.initial_state()
        return np.r_[d,self.receiver.initial_state(pb1=d[6])]

    def rhs(self, t, y):
        if np.shape(y) != (27,):
            raise ValueError('Expected 27 states')
        return np.r_[self.donor.rhs(t,y[:10]), self.receiver.rhs(t,y[10:],pb1=y[6],int0=y[1])]


def time_grid(hours, sample_min):
    if not np.isfinite(hours) or not np.isfinite(sample_min) or hours <= 0 or sample_min <= 0:
        raise ValueError('Positive finite duration and sampling required')
    if hours > 50:
        raise ValueError('Original upstream data ends at 50 h; no extrapolation allowed')
    step=sample_min/60
    t = np.arange(int(np.floor(hours/step))+1)*step
    t = t[t < hours-1e-12*max(1.,hours)]
    return np.r_[t,hours]


def integrate(rhs, y0, t, method='DOP853', rtol=2e-7, atol=2e-9, max_step_min=1.):
    if max_step_min <= 0:
        raise ValueError('max_step_min must be positive')
    sol = solve_ivp(rhs,(float(t[0]),float(t[-1])),y0,t_eval=t,method=method,
                    rtol=rtol,atol=atol,max_step=max_step_min/60)
    if not sol.success or not np.all(np.isfinite(sol.y)):
        raise RuntimeError(f'Integration failed: {sol.message}')
    return sol
