"""23-state HBY bit0 -> current ZMH middle -> exploratory ZMH tail.

Time is hours. The middle retains its source per-minute parameters and the
existing HZH conversion by 60. The new ZMH tail is already in hours.
"""
from pathlib import Path
import ast
import sys
import numpy as np
from scipy.optimize import brentq

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
PARENT = ROOT.parent
for directory in (PARENT, PARENT / 'hby_zmh_hby'):
    sys.path.insert(0, str(directory))
from model_hzh import HZHModel, NAMES as HZH_NAMES
from hybrid_model import require_hash

TAIL_SOURCE = PARENT.parent / 'wiki任务' / '前馈三级级联.py'
TAIL_SHA = 'C28FB9587EF12948F40ECB626B117E7C4757DB1F63001F794F8FE4FF19206CEA'
NAMES = HZH_NAMES[:17] + ('A1_zmh', 'F1_zmh', 'pb2_zmh', 'I2_zmh', 'T2_zmh', 'RDF2_zmh')
IDX = {name: i for i, name in enumerate(NAMES)}


def source_namespace():
    """Extract literal parameters and original RHS without its simulations/plots."""
    require_hash(TAIL_SOURCE, TAIL_SHA)
    tree = ast.parse(TAIL_SOURCE.read_text(encoding='utf-8-sig'))
    node = next(n for n in tree.body if isinstance(n, ast.Assign)
                and any(isinstance(t, ast.Name) and t.id == 'p' for t in n.targets))
    p = ast.literal_eval(node.value)
    fun = next(n for n in tree.body if isinstance(n, ast.FunctionDef)
               and n.name == 'counter_3bit_ode')
    ns = {'get_u_in': lambda t: 0.0}
    exec(compile(ast.Module(body=[fun], type_ignores=[]), str(TAIL_SOURCE), 'exec'), ns)
    return p, ns['counter_3bit_ode']


class HZZModel:
    names = NAMES
    clock_K = 0.4
    clock_n = 3.0

    def __init__(self, prefix, autoregulation=False):
        self.prefix = prefix
        self.autoregulation = bool(autoregulation)
        self.p, self.original_rhs = source_namespace()

    def tail_signals(self, y, int0):
        a, f, pb, i, tr, r = y
        p = self.p
        act = a**p['n_A1'] / (p['K_A1']**p['n_A1'] + a**p['n_A1'])
        repress = p['K_R1']**p['n_R1'] / (p['K_R1']**p['n_R1'] + f**p['n_R1'])
        clock = int0**self.clock_n / (self.clock_K**self.clock_n + int0**self.clock_n)
        auto = p['K_auto1']**p['n_auto1'] / (p['K_auto1']**p['n_auto1'] + a**p['n_auto1'])
        gate = act * repress * clock
        vf = p['k_fwd'] * pb * i**2 / (p['K_D_int2']**2+i**2) * p['Kinh']/(p['Kinh']+r)
        vr = p['k_rev'] * (1-pb) * (i*r)**2 / (p['K_D_comp']**2+(i*r)**2)
        return dict(g1=gate, clock=clock, auto=auto, activation=act, repression=repress,
                    int2_source=p['alpha_Int2']*gate, J_fwd2=vf, J_rev2=vr)

    def tail_rhs(self, t, y, pb1, int0):
        a, f, pb, i, tr, r = y
        p = self.p
        q = self.tail_signals(y, int0)
        auto = q['auto'] if self.autoregulation else 1.0
        return np.array([
            p['alpha_A1']*pb1*auto-p['gamma_A1']*a,
            p['alpha_R1']*q['activation']-p['gamma_R1']*f,
            -q['J_fwd2']+q['J_rev2'],
            q['int2_source']-p['gamma_int2']*i,
            p['alpha_rep']*pb-p['gamma_rep']*tr,
            p['alpha_rdf']*(1-pb)/(1+(tr/p['K_rep'])**p['n'])-p['gamma_rdf']*r])

    def tail_initial(self, pb1=1.0):
        p = self.p
        if self.autoregulation:
            a = brentq(lambda a: p['alpha_A1']*pb1*p['K_auto1']**p['n_auto1'] /
                       (p['K_auto1']**p['n_auto1']+a**p['n_auto1'])-p['gamma_A1']*a,
                       0, p['alpha_A1']/p['gamma_A1'])
        else:
            a = p['alpha_A1']*pb1/p['gamma_A1']
        act = a**p['n_A1']/(p['K_A1']**p['n_A1']+a**p['n_A1'])
        return np.array([a, p['alpha_R1']*act/p['gamma_R1'], 1., 0.,
                         p['alpha_rep']/p['gamma_rep'], 0.])

    def initial_state(self):
        low = self.prefix.initial_state()[:17].copy()
        return np.r_[low, self.tail_initial(low[13])]

    def rhs(self, t, y):
        if np.shape(y) != (23,):
            raise ValueError('HZZ requires 23 states')
        return np.r_[self.prefix.bit0_rhs(t, y[:11]),
                     self.prefix.middle_rhs(y[10], y[11:17]),
                     self.tail_rhs(t, y[17:], pb1=y[13], int0=y[2])]

    def signals(self, y):
        p, z = self.prefix.p, self.prefix.z
        i0, r0, c0, s0 = y[2], y[8], y[9], y[10]
        a0, f0, pb1, i1, tr1, r1 = y[11:17]
        act0 = a0**z['n_A0']/(z['K_A0']**z['n_A0']+a0**z['n_A0'])
        gate0 = z['K_R0']**z['n_R0']/(z['K_R0']**z['n_R0']+f0**z['n_R0'])
        q = self.tail_signals(y[17:], i0)
        q.update(int0=i0, S0=s0, S1=1-pb1, S2=1-y[19], g0=act0*gate0,
                 J_rev0=p['k_rev']*s0*c0**2/(self.prefix.tail.k_complex**2+c0**2),
                 J_rev1=60*z['k_rev']*(1-pb1)*(i1*r1)**2/(z['K_D_comp']**2+(i1*r1)**2))
        return q


def structural_checks(prefix):
    off, on = HZZModel(prefix, False), HZZModel(prefix, True)
    rng = np.random.default_rng(20260930)
    max_tail_gap = max_prefix_gap = max_toggle_gap = 0.0
    for _ in range(40):
        y = rng.uniform(.01, 5, 23)
        y[[10, 13, 19]] = rng.uniform(0, 1, 3)
        old = np.zeros(16)
        old[:4] = [1-y[10], y[2], y[5], y[8]]
        old[4:10] = y[11:17]
        old[10:] = y[17:]
        actual = off.tail_rhs(1., y[17:], y[13], y[2])
        expected = np.asarray(off.original_rhs(1., old, off.p))[10:]
        max_tail_gap = max(max_tail_gap, float(np.max(np.abs(actual-expected))))
        changed = on.tail_rhs(1., y[17:], y[13], y[2])-actual
        expected_change = off.p['alpha_A1']*y[13]*(off.tail_signals(y[17:], y[2])['auto']-1)
        max_toggle_gap = max(max_toggle_gap, abs(changed[0]-expected_change), float(np.max(np.abs(changed[1:]))))
        y34 = prefix.initial_state(); y34[:17] = y[:17]
        ref = prefix.rhs(1., y34)[:17]
        max_prefix_gap = max(max_prefix_gap, float(np.max(np.abs(off.rhs(1., y)[:17]-ref))))
        y[17:] *= 100
        max_prefix_gap = max(max_prefix_gap, float(np.max(np.abs(on.rhs(1., y)[:17]-ref))))
    initial_residuals = {name: float(np.max(np.abs(m.tail_rhs(0, m.tail_initial(), 1, 0))))
                         for name, m in [('off', off), ('on', on)]}
    checks = dict(original_tail_rhs_gap=max_tail_gap, prefix_rhs_gap=max_prefix_gap,
                  toggle_only_A1_gap=max_toggle_gap, initial_tail_residuals=initial_residuals,
                  original_tail_autoregulation=False, probes=40)
    if max_tail_gap > 1e-12 or max_prefix_gap != 0 or max_toggle_gap > 1e-12 or max(initial_residuals.values()) > 1e-10:
        raise RuntimeError(checks)
    return checks
