"""Item 3: middle-module replacement contrast.

Holds everything else fixed - same Han upstream, same HBY bit0, same HBY
A1/F1 + bit2 receiver, same uM_per_au and clock parameters - and swaps ONLY
the middle stage:

  arm 'zmh'  : ZMH's current reduced A0/F0 + bit1 module
               (A0_zmh, F0_zmh, pb1_zmh, I1_zmh, T1_zmh, RDF1_zmh; 6 states)
  arm 'm51'  : the 51-state model's carry0 expression chain + full 11-state bit1
               (M_A0, A0_u, A0, M_F0, F0_u, F0) + b1(11) = 17 states

State counts: zmh 11+6+17 = 34 (identical to hby_zmh_hby/model_hzh.py),
              m51 11+17+17 = 45.

The frozen 51-state extension and HbyConfig agree numerically on every shared
quantity (mRNA 2.0 min, bit maturation 20.0 min, carry maturation 32.5 min,
translation_h 30, kon/koff/complex_decay 0.1/1.0/1.0, uM_per_au 5.75,
clock 0.3/2), so this swap changes the module STRUCTURE and the middle-stage
PARAMETER TABLE, and nothing else. Those two are not separated here.

The criterion is applied with the same verify_hzh functions, on a signal dict
that is built for each layout, so both arms are judged identically.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from dataclasses import asdict
from datetime import datetime
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
PARENT = ROOT.parent
for _p in (PARENT, PARENT / 'hby_zmh_hby', PARENT / 'hby_zmh_hby' / 'certification',
           PARENT.parent / 'final_reconstruction'):
    sys.path.insert(0, str(_p))

import scan_oat as S                                              # noqa: E402
import verify_hzh as VH                                           # noqa: E402
from hybrid_model import HbyConfig, HbyReceiver, hill            # noqa: E402
from model_hzh import HZHModel, NAMES as NAMES_HZH                # noqa: E402
from model import ZENG as ZENG51                                  # noqa: E402

FROZEN = S.FROZEN
BIT_NAMES = ('M_I', 'I_u', 'I', 'M_T', 'T_u', 'T', 'M_R', 'R_u', 'R', 'C', 'S')
M51_NAMES = (tuple('b0_' + n for n in BIT_NAMES)
             + ('M_A0', 'A0_u', 'A0', 'M_F0', 'F0_u', 'F0')
             + tuple('b1_' + n for n in BIT_NAMES)
             + ('A1', 'F1', 'M_A1', 'A1_u', 'M_F1', 'F1_u')
             + tuple('b2_' + n for n in BIT_NAMES))
N_M51 = len(M51_NAMES)
M51_IDX = {n: i for i, n in enumerate(M51_NAMES)}
assert N_M51 == 45, N_M51


class SwapModel:
    """HBY bit0 -> middle (swappable) -> HBY A1/F1 + bit2."""

    def __init__(self, mode, han, config=None, middle='m51',
                 b1_rdf_scale=1.0):
        if mode not in ('square', 'han'):
            raise ValueError(mode)
        if mode == 'han' and han is None:
            raise ValueError('Han upstream required')
        if middle not in ('zmh', 'm51'):
            raise ValueError(middle)
        self.mode, self.han, self.middle = mode, han, middle
        # single-lever probe: scales alpha_rdf for the middle bit1 ONLY, so the
        # RDF1 pool / bit1 reverse flux can be moved without touching bit0,
        # bit2 or the carry chain.
        self.b1_rdf_scale = float(b1_rdf_scale)
        self.tail = HbyReceiver(config)
        self.c = self.tail.c
        self.rc = self.tail                      # reuse expression/rates
        self.p = self.tail.p                     # HBY/ZENG table
        self.rhs_model = HZHModel(mode, han)     # for h31 only
        self.z = self.rhs_model.z
        if han is not None:
            self.h31 = self.rhs_model.h31
            self.copies_per_au = 602.214076 * self.c.receiver_uM_per_au
        self.square = self.rhs_model.square

    # ---- bit0: same equations as the frozen hybrid, own config ----------
    def bit0_rhs(self, t, b):
        p, c = self.p, self.c
        m, iu, i, mt, tu, tr, mr, ru, r, comp, s = b
        lm, km = self.rc.lmb, self.rc.kmb
        d = np.zeros(11)
        if self.mode == 'square':
            d[0] = self.rhs_model.input_source(t) * lm / c.translation_h - lm * m
            translation = c.translation_h * m
        else:
            hp = self.han.p
            tx = (60 * hp.oscillator_plasmid_copies
                  * hp.c31_tx_per_plasmid_per_min
                  * self.han.interp(60 * t, self.h31) / self.copies_per_au)
            loss = 60 * (hp.c31_mrna_intrinsic_loss_per_min + hp.mu)
            d[0] = tx - loss * m
            translation = 60 * hp.c31_translation_per_mrna_per_min * m
        d[1] = translation - km * iu
        d[2] = km * iu - p['gamma_int'][0] * i
        d[3:6] = self.rc.expression(mt, tu, tr, p['alpha_rep'] * (1 - s),
                                    p['gamma_rep'], lm, km)
        d[6:9] = self.rc.expression(
            mr, ru, r,
            p['alpha_rdf'] * s * (1 - hill(tr, p['K_rep'], p['n_rep'])),
            p['gamma_rdf'], lm, km)
        bind, un = c.kon * i * r, c.koff * comp
        d[2] += -bind + un
        d[8] += -bind + un
        d[9] = bind - un - c.complex_decay * comp
        vf = (p['k_fwd'] * hill(i, p['K_D_int'][0], 2)
              * p['K_inh'] / (p['K_inh'] + max(r, 0)))
        vr = p['k_rev'] * hill(comp, self.rc.k_complex, 2)
        d[10] = vf * (1 - s) - vr * s
        return d

    def initial_state(self):
        b0 = self.rc.initial_state()[6:].copy()
        if self.mode == 'han':
            b0[0] = self.han.sol.y[6, 0] / self.copies_per_au
        if self.middle == 'zmh':
            p = self.z
            a = p['alpha_A0'] / p['gamma_A0']
            act = a ** p['n_A0'] / (p['K_A0'] ** p['n_A0'] + a ** p['n_A0'])
            mid = [a, p['alpha_R0'] * act / p['gamma_R0'], 1., 0.,
                   p['alpha_rep1'] / p['gamma_rep1'], 0.]
            return np.r_[b0, mid, self.rc.initial_state()]
        # m51: carry0 expression chain + full bit1, both pre-equilibrated
        P = self.p
        pb0 = 1.0 - b0[10]                       # bit0 starts PB
        srcA = ZENG51['alpha_A'][0] * pb0
        A0 = srcA / ZENG51['gamma_A'][0]
        srcF = (ZENG51['alpha_F'][0]
                * hill(A0, ZENG51['K_A'][0], ZENG51['n_A'][0]))
        c0 = [srcA / self.c.translation_h, srcA / self.rc.kmc, A0,
              srcF / self.c.translation_h, srcF / self.rc.kmc,
              srcF / ZENG51['gamma_F'][0]]
        b1 = np.zeros(11)
        b1[3:6] = (ZENG51['alpha_rep'] / self.c.translation_h,
                   ZENG51['alpha_rep'] / self.rc.kmb,
                   ZENG51['alpha_rep'] / ZENG51['gamma_rep'])
        return np.r_[b0, c0, b1, self.rc.initial_state(pb1=1.0)]

    # ---- middle-stage bit body (shared by b1 and b2 shapes) ------------
    def _bit_body(self, b, source_int, gamma_int, rdf_scale=1.0):
        P, c, Z = self.p, self.c, ZENG51
        mI, Iu, I, mT, Tu, T, mR, Ru, R, C, S = b
        db = np.zeros(11)
        db[:3] = self.rc.expression(mI, Iu, I, source_int, gamma_int,
                                    self.rc.lmb, self.rc.kmb)
        db[3:6] = self.rc.expression(mT, Tu, T, Z['alpha_rep'] * (1 - S),
                                     Z['gamma_rep'], self.rc.lmb, self.rc.kmb)
        db[6:9] = self.rc.expression(
            mR, Ru, R,
            rdf_scale * Z['alpha_rdf'] * S
            * (1 - hill(T, Z['K_rep'], Z['n_rep'])),
            Z['gamma_rdf'], self.rc.lmb, self.rc.kmb)
        bind, unbind = c.kon * I * R, c.koff * C
        db[2] += -bind + unbind
        db[8] += -bind + unbind
        db[9] = bind - unbind - c.complex_decay * C
        vf = (Z['k_fwd'] * hill(I, Z['K_D_int'][1], 2)
              * Z['K_inh'] / (Z['K_inh'] + max(R, 0.0)))
        vr = Z['k_rev'] * hill(C, self.rc.k_complex, 2)
        db[10] = vf * (1 - S) - vr * S
        return db

    def rhs(self, t, y):
        if self.middle == 'zmh':
            if np.shape(y) != (34,):
                raise ValueError('zmh arm requires 34 states')
            return np.r_[self.bit0_rhs(t, y[:11]),
                         self.rhs_model.middle_rhs(y[10], y[11:17]),
                         self.rc.rhs(t, y[17:], pb1=y[13], int0=y[2])]
        if np.shape(y) != (N_M51,):
            raise ValueError(f'm51 arm requires {N_M51} states')
        b0 = y[:11]
        i = 11
        mA, uA, A0, mF, uF, F0 = y[i:i + 6]
        b1 = y[i + 6:i + 17]
        tail = y[i + 17:]
        Z = ZENG51
        pb0 = 1.0 - b0[10]
        srcA = Z['alpha_A'][0] * pb0
        srcF = Z['alpha_F'][0] * hill(A0, Z['K_A'][0], Z['n_A'][0])
        eA = self.rc.expression(mA, uA, A0, srcA, Z['gamma_A'][0],
                                self.rc.lmc, self.rc.kmc)
        eF = self.rc.expression(mF, uF, F0, srcF, Z['gamma_F'][0],
                                self.rc.lmc, self.rc.kmc)
        g0 = (hill(A0, Z['K_A'][0], Z['n_A'][0])
              * (1 - hill(F0, Z['K_F'][0], Z['n_F'][0])))
        d = np.zeros(N_M51)
        d[:11] = self.bit0_rhs(t, b0)
        d[i:i + 6] = np.r_[eA, eF]
        d[i + 6:i + 17] = self._bit_body(
            b1, Z['alpha_Int'][0] * g0, Z['gamma_int'][1],
            rdf_scale=self.b1_rdf_scale)
        pb1 = 1.0 - b1[10]
        d[i + 17:] = self.rc.rhs(t, tail, pb1=pb1, int0=b0[2])
        return d

    # ---- signals, in the same key set verify_hzh expects ---------------
    def _tail_signals(self, tail, int0):
        """Vectorised equivalent of HbyReceiver.signals (which is scalar-only)."""
        P, c = self.p, self.c
        a, f = tail[0], tail[1]
        i, r, comp, s = tail[8], tail[14], tail[15], tail[16]
        clock = self._hill_arr(c.clock_scale * int0, c.clock_K_au, c.clock_n)
        g1 = (self._hill_arr(a, P['K_A'][1], c.n_A1_gate)
              * (1 - self._hill_arr(f, P['K_F'][1], P['n_F'][1])) * clock)
        vf = (P['k_fwd'] * self._hill_arr(i, P['K_D_int'][2], 2)
              * P['K_inh'] / (P['K_inh'] + np.maximum(r, 0.0)))
        vr = P['k_rev'] * self._hill_arr(comp, self.rc.k_complex, 2)
        return dict(clock_gate=clock, g1=g1,
                    u2_target_au_per_h=P['alpha_Int'][1] * g1,
                    J_fwd2_per_h=vf * (1 - s), J_rev2_per_h=vr * s)

    def signals(self, y):
        P, c = self.p, self.c
        if self.middle == 'zmh':
            s = VH.signals(y, self.z, self.rc)
            s['int0'] = y[2]
            # u2_target is alpha_Int[1]*g1 in both arms; VH.signals omits it
            z = self.z
            a0, f0 = y[11], y[12]
            act0 = a0 ** z['n_A0'] / (z['K_A0'] ** z['n_A0'] + a0 ** z['n_A0'])
            gate0 = act0 * (z['K_R0'] ** z['n_R0']
                            / (z['K_R0'] ** z['n_R0'] + f0 ** z['n_R0']))
            return s, dict(
                u2_target=self.p['alpha_Int'][1] * s['g1'],
                PB0=1 - y[10], PB1=y[13],
                A0=a0, F0=f0, Int1=y[14], T1=y[15], RDF1=y[16],
                Int1_source=z['alpha_Int1'] * gate0)
        Z = ZENG51
        b0 = y[:11]
        b1 = y[17:28]
        tail = y[28:]
        A0, F0 = y[13], y[16]
        g0 = (self._hill_arr(A0, Z['K_A'][0], Z['n_A'][0])
              * (1 - self._hill_arr(F0, Z['K_F'][0], Z['n_F'][0])))
        ts = self._tail_signals(tail, b0[2])
        s0 = b0[10]
        s1 = b1[10]
        c0m, c1m = b0[9], b1[9]
        rev0 = P['k_rev'] * self._hill_arr(c0m, self.rc.k_complex, 2) * s0
        rev1 = P['k_rev'] * self._hill_arr(c1m, self.rc.k_complex, 2) * s1
        out = dict(g0=g0, g1=ts['g1'], clock=ts['clock_gate'],
                   J_rev0=rev0, J_rev1=rev1,
                   J_rev2=ts['J_rev2_per_h'], J_fwd2=ts['J_fwd2_per_h'],
                   S0=s0, S1=s1, S2=tail[16], int0=b0[2])
        extra = dict(A0=A0, F0=F0, PB0=1 - b0[10], PB1=1 - b1[10],
                     Int1=b1[2], RDF1=b1[8], T1=b1[5],
                     Int1_source=Z['alpha_Int'][0] * g0,
                     u2_target=P['alpha_Int'][1] * ts['g1'])
        return out, extra

    @staticmethod
    def _hill_arr(x, k, n):
        z = np.maximum(np.asarray(x, dtype=float), 0.0) / k
        return z ** n / (1 + z ** n)


def criterion(t, sig):
    """Apply HZH_MOD8_CAUSAL_V1 to an arbitrary signal dict."""
    peaks, reads = VH.read_windows(t, sig)
    cross = {b: VH.crossings(t, sig[b]) for b in ('S0', 'S1', 'S2')}
    events = {}
    for stage, revname, gatename, bit in ((0, 'J_rev0', 'g0', 'S1'),
                                          (1, 'J_rev1', 'g1', 'S2')):
        rev = VH.segments(t, sig[revname], VH.JREV_THRESHOLD_H)
        gates = VH.segments(t, sig[gatename], VH.GATE_THRESHOLD,
                            min_dose=VH.GATE_MIN_DOSE_H)
        events[f'bit{stage}_to_bit{stage+1}'] = VH.associate(
            rev, gates, cross[bit], float(t[-1]))
    steady = VH.read_verdict(reads, VH.DROP)
    bit_margins = {b: VH.margins(reads, cross[b]) for b in cross}
    nonnull = [x for m in bit_margins.values()
               for x in (m['min_setup_h'], m['min_hold_h']) if x is not None]
    return dict(
        counting_passed=bool(steady['increments_mod8']),
        event_causality_passed=bool(all(e['passed'] for e in events.values())),
        full_certified=bool(steady['passed']
                            and all(e['passed'] for e in events.values())),
        steady=steady, reads=reads, crossings=cross, events=events,
        bit_margins=bit_margins, clock_peak_count=len(peaks),
        margin_global_h=min(nonnull) if nonnull else None)


def cycle_stats(t, y, reads, name_key):
    peaks, widths, integrals, residuals = [], [], [], []
    for r in reads:
        m = (t >= r['cycle_start_h']) & (t <= r['cycle_end_h'])
        if m.sum() < 3:
            continue
        tt, yy = t[m], np.asarray(y[m], dtype=float)
        pk = float(yy.max())
        peaks.append(pk)
        if pk > 0:
            above = tt[yy >= pk / 2]
            widths.append(float(above.max() - above.min()) if above.size > 1 else 0.0)
        integrals.append(float(np.trapezoid(np.maximum(yy, 0.0), tt)))
        k = max(1, int(0.2 * tt.size))
        residuals.append(float(np.median(yy[np.argsort(yy)[:k]])))
    return dict(peak_max=max(peaks) if peaks else None,
                width_max=max(widths) if widths else None,
                integral_min=min(integrals) if integrals else None,
                integral_max=max(integrals) if integrals else None,
                far_dwell_median=float(np.median(residuals)) if residuals else None)


def run(han, config, middle, n_gate, label, b1_rdf_scale=1.0):
    cfg = HbyConfig(**{**asdict(config), 'n_A1_gate': n_gate})
    m = SwapModel('han', han, cfg, middle, b1_rdf_scale=b1_rdf_scale)
    tic = time.perf_counter()
    t, y = S.integrate(m)
    sig, extra = m.signals(y)
    v = criterion(t, sig)
    reads = v['reads']
    inside = np.zeros(t.size, dtype=bool)
    for r in reads:
        inside |= (t >= r['start_h']) & (t <= r['end_h'])
    outside = ~inside
    row = dict(
        label=label, middle=middle, n_A1_gate=n_gate,
        b1_rdf_scale=float(b1_rdf_scale),
        n_states=int(y.shape[0]), runtime_s=time.perf_counter() - tic,
        counting_passed=v['counting_passed'],
        event_causality_passed=v['event_causality_passed'],
        full_certified=v['full_certified'],
        sequence=v['steady']['sequence'], reads=int(v['steady']['reads']),
        commitment=(None if v['steady']['minimum_commitment'] is None
                    else float(v['steady']['minimum_commitment'])),
        margin_global_h=v['margin_global_h'],
        margin_S1_h=S.min_of(v['bit_margins'].get('S1')),
        margin_S2_h=S.min_of(v['bit_margins'].get('S2')),
        clock_peak_count=v['clock_peak_count'],
        s0={k: v['events']['bit0_to_bit1'][k] for k in
            ('reverse_events', 'gate_events', 'flip_events', 'one_to_one',
             'causal_order', 'alternating_directions')},
        s1={k: v['events']['bit1_to_bit2'][k] for k in
            ('reverse_events', 'gate_events', 'flip_events', 'one_to_one',
             'causal_order', 'alternating_directions')},
        g1_peak=float(np.max(sig['g1'])),
        g1_window=cycle_stats(t, sig['g1'], reads, 'g1'),
        int2_peak=float(np.max(extra['u2_target'])),
        j_fwd2_integral=float(np.trapezoid(np.maximum(sig['J_fwd2'], 0), t)),
        j_rev2_integral=float(np.trapezoid(np.maximum(sig['J_rev2'], 0), t)),
        j_rev1_integral=float(np.trapezoid(np.maximum(sig['J_rev1'], 0), t)),
        interface=dict(
            PB0_high=float(np.median(extra['PB0'][inside])) if inside.any() else None,
            PB0_low=float(np.median(extra['PB0'][outside])) if outside.any() else None,
            PB1_high=float(np.median(extra['PB1'][inside])) if inside.any() else None,
            PB1_low=float(np.median(extra['PB1'][outside])) if outside.any() else None,
            A0=max_med(extra['A0'], inside, outside),
            F0=max_med(extra['F0'], inside, outside),
            Int1=max_med(extra['Int1'], inside, outside),
            T1=max_med(extra['T1'], inside, outside),
            RDF1=max_med(extra['RDF1'], inside, outside),
            Int1_source=max_med(extra['Int1_source'], inside, outside),
            Int1_peak=float(np.max(extra['Int1'])),
            Int1_source_peak=float(np.max(extra['Int1_source'])),
        ),
    )
    return row, t, y, sig, extra


def max_med(x, inside, outside):
    return dict(peak=float(np.max(x)),
                median_in=float(np.median(x[inside])) if inside.any() else None,
                median_out=float(np.median(x[outside])) if outside.any() else None)


def _integrate(model, y0):
    from scipy.integrate import solve_ivp
    count = round(S.HOURS * 60 / S.SAMPLE_MIN)
    t = np.linspace(0, S.HOURS, count + 1)
    sol = solve_ivp(model.rhs, (0, S.HOURS), y0, t_eval=t, method='DOP853',
                    rtol=S.RTOL, atol=S.ATOL, max_step=S.MAX_STEP_MIN / 60)
    if not sol.success or not np.all(np.isfinite(sol.y)):
        raise RuntimeError(sol.message)
    return sol.t, sol.y


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--hours', type=float, default=None)
    ap.add_argument('--no-lever', action='store_true')
    args = ap.parse_args()
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8')
    if args.hours:
        S.HOURS = args.hours

    out = ROOT / 'results' / ('swap_' + datetime.now().strftime('%Y%m%d_%H%M%S_%f'))
    out.mkdir(parents=True, exist_ok=False)
    S.dump(out / 'status.json', dict(status='RUNNING'))

    print('building Han upstream once', flush=True)
    han = S.HanInput(S.HOURS)
    assert han.p.peak_load_fraction == 0.0

    rows = []
    for middle in ('zmh', 'm51'):
        for n in (4.0, 6.0):
            label = f'{middle} n={n:g}'
            row, t, y, sig, extra = run(han, FROZEN, middle, n, label)
            rows.append(row)
            print(f'[{label:12s}] states={row["n_states"]:2d} '
                  f'cnt={int(row["counting_passed"])} '
                  f'ev={int(row["event_causality_passed"])} '
                  f'full={int(row["full_certified"])} '
                  f'reads={row["reads"]:3d} seq={row["sequence"]} '
                  f'S2m={row["margin_S2_h"]} '
                  f'g1pk={row["g1_peak"]:.4f} '
                  f'Int2pk={row["int2_peak"]:.3f} '
                  f'{row["runtime_s"]:.0f}s', flush=True)
            np.savez_compressed(out / f'traj_{middle}_n{int(n)}.npz',
                                time_h=t, states=y,
                                names=np.array(NAMES_HZH if middle == 'zmh'
                                               else M51_NAMES))
            S.dump(out / f'signals_{middle}_n{int(n)}.json',
                   {k: [float(np.min(v)), float(np.max(v))]
                    for k, v in sig.items()})

    # single-lever probe: does shrinking the middle bit1 RDF pool (and hence its
    # reverse flux) restore n=4 in the m51 arm?
    if not args.no_lever:
        for sc in (0.5, 0.25, 0.1):
            label = f'm51 n=4 rdf_x{sc:g}'
            row, t, y, sig, extra = run(han, FROZEN, 'm51', 4.0, label,
                                        b1_rdf_scale=sc)
            rows.append(row)
            print(f'[{label:20s}] cnt={int(row["counting_passed"])} '
                  f'ev={int(row["event_causality_passed"])} '
                  f'full={int(row["full_certified"])} seq={row["sequence"]} '
                  f'PB1hi={row["interface"]["PB1_high"]:.4g} '
                  f'Jrev1={row["j_rev1_integral"]:.3g} '
                  f'g1pk={row["g1_peak"]:.4f} '
                  f'S2m={row["margin_S2_h"]}', flush=True)

    # control: the zmh arm at n=6 must reproduce the frozen certification
    ctl = next(r for r in rows if r['middle'] == 'zmh' and r['n_A1_gate'] == 6.0)
    checks = dict(full=ctl['full_certified'] is True,
                  reads=ctl['reads'] == S.REF['steady_reads'],
                  sequence=ctl['sequence'] == S.REF['sequence'],
                  margin=abs(ctl['margin_global_h'] - S.REF['margin_h']) < 1e-6)
    print('control (zmh n=6 vs frozen):', json.dumps(checks), flush=True)

    S.dump(out / 'swap_all.json', dict(
        hours=S.HOURS, sample_min=S.SAMPLE_MIN, rtol=S.RTOL, atol=S.ATOL,
        frozen_config=asdict(FROZEN), reference=S.REF, control_checks=checks,
        middle_definition=dict(
            zmh='A0_zmh,F0_zmh,pb1_zmh,I1_zmh,T1_zmh,RDF1_zmh (6 states)',
            m51='M_A0,A0_u,A0,M_F0,F0_u,F0 + full 11-state b1 (17 states)'),
        confound_note=('the swap changes module structure AND the middle-stage '
                       'parameter table (ZENG vs ZMH current); they are not '
                       'separated here'),
        rows=rows))
    S.dump(out / 'status.json', dict(status='COMPLETED', rows=len(rows)))
    print('OUTPUT', out, flush=True)


if __name__ == '__main__':
    main()
