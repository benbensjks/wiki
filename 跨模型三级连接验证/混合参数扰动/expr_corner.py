"""Expression-time corner test at the interface boundaries.

Question: does expression timing move the working region, or does it only set
whether that region is dynamically reachable?

Points (the nominal point plus the four measured interface boundaries)
  1 nominal
  2 uM lower boundary near  (x0.80)
  3 uM upper boundary near  (x1.30)
  4 K_A[1] upper near       (n=6: x1.30)
  5 K_F[1] upper near       (x1.40)

Time configurations
  all_fast        all four timing parameters x0.5
  nominal         x1
  all_slow        all four timing parameters x2
  A1fast_F1slow   carry A1 maturation x0.5, carry F1 maturation x2
  A1slow_F1fast   carry A1 maturation x2,   carry F1 maturation x0.5

The four timing parameters are carry_mrna_min, carry_maturation_min,
bit_mrna_min, bit_maturation_min.

LOCAL EXTENSION, not a change to any existing file
--------------------------------------------------
HbyConfig carries a single `carry_maturation_min` used for both A1 and F1.
The mismatch configurations therefore need separate A1/F1 maturation rates, so
this script defines SplitCarryReceiver(HbyReceiver) with a copied rhs and
initial_state that use self.kmc_A / self.kmc_F. The original class in
hybrid_model.py is untouched and is still used for the nominal configuration.

Every run records J_fwd2 and J_rev2 integrals in addition to pass/fail, because
a pass/fail boolean hides the case where the count is right but the event
segmentation or the write dose has moved.
"""
from __future__ import annotations

import argparse, copy, csv, json, sys, time
from dataclasses import asdict, replace
from datetime import datetime
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

import scan_oat as S                                              # noqa: E402
import verify_hzh as VH                                           # noqa: E402
import perbit_margin as PB                                        # noqa: E402
from hybrid_model import HbyConfig, HbyReceiver, hill            # noqa: E402
from model_hzh import IDX as IDX_HZH                              # noqa: E402

FROZEN, ZENG = S.FROZEN, S.ZENG
B2_I = IDX_HZH['b2_I']

TIME_KEYS = ('carry_mrna_min', 'carry_maturation_min',
             'bit_mrna_min', 'bit_maturation_min')


class SplitCarryReceiver(HbyReceiver):
    """A1 and F1 mature at independent rates. Local subclass, read-only use."""

    def __init__(self, config, a1_mat_min, f1_mat_min):
        super().__init__(config)
        self.kmc_A = np.log(2) * 60 / a1_mat_min
        self.kmc_F = np.log(2) * 60 / f1_mat_min

    def rhs(self, t, y, pb1, int0):
        p, c = self.p, self.c
        a, f, ma, au, mf, fu = y[:6]
        mi, iu, i, mt, tu, tr, mr, ru, r, comp, s = y[6:]
        sa = p['alpha_A'][1] * pb1 * (1 - hill(a, p['K_auto1'], p['n_auto1']))
        sf = p['alpha_F'][1] * hill(a, p['K_A'][1], p['n_A'][1])
        ea = self.expression(ma, au, a, sa, p['gamma_A'][1], self.lmc, self.kmc_A)
        ef = self.expression(mf, fu, f, sf, p['gamma_F'][1], self.lmc, self.kmc_F)
        sig = self.signals(y, int0)
        db = np.zeros(11)
        db[:3] = self.expression(mi, iu, i, sig['u2_target_au_per_h'],
                                 p['gamma_int'][2], self.lmb, self.kmb)
        db[3:6] = self.expression(mt, tu, tr, p['alpha_rep'] * (1 - s),
                                  p['gamma_rep'], self.lmb, self.kmb)
        db[6:9] = self.expression(
            mr, ru, r,
            p['alpha_rdf'] * s * (1 - hill(tr, p['K_rep'], p['n_rep'])),
            p['gamma_rdf'], self.lmb, self.kmb)
        bind, unbind = c.kon * i * r, c.koff * comp
        db[2] += -bind + unbind
        db[8] += -bind + unbind
        db[9] = bind - unbind - c.complex_decay * comp
        db[10] = sig['J_fwd2_per_h'] - sig['J_rev2_per_h']
        return np.r_[ea[2], ef[2], ea[0], ea[1], ef[0], ef[1], db]

    def initial_state(self, pb1=1.):
        p = self.p
        a = brentq(lambda a: p['alpha_A'][1] * pb1
                   * (1 - hill(a, p['K_auto1'], p['n_auto1']))
                   - p['gamma_A'][1] * a, 0, p['alpha_A'][1] / p['gamma_A'][1])
        sa = p['gamma_A'][1] * a
        sf = p['alpha_F'][1] * hill(a, p['K_A'][1], p['n_A'][1])
        carry = [a, sf / p['gamma_F'][1], sa / self.c.translation_h,
                 sa / self.kmc_A, sf / self.c.translation_h, sf / self.kmc_F]
        bit = np.zeros(11)
        bit[3:6] = (p['alpha_rep'] / self.c.translation_h,
                    p['alpha_rep'] / self.kmb, p['alpha_rep'] / p['gamma_rep'])
        return np.r_[carry, bit]


# point label -> (uM fold, K_A[1] fold, K_F[1] fold, n_A1_gate)
POINTS = {
    'nominal':      (1.00, 1.00, 1.00, 6.0),
    'uM_low':       (0.80, 1.00, 1.00, 6.0),
    'uM_high':      (1.30, 1.00, 1.00, 6.0),
    'K_A_high':     (1.00, 1.30, 1.00, 6.0),
    'K_F_high':     (1.00, 1.00, 1.40, 6.0),
}

# config label -> {timing key: factor} plus optional a1/f1 maturation factors
CONFIGS = {
    'all_fast':      {k: 0.5 for k in TIME_KEYS},
    'nominal':       {k: 1.0 for k in TIME_KEYS},
    'all_slow':      {k: 2.0 for k in TIME_KEYS},
    'A1fast_F1slow': {'carry_mrna_min': 1.0, 'bit_mrna_min': 1.0,
                      'bit_maturation_min': 1.0, 'a1_mat': 0.5, 'f1_mat': 2.0},
    'A1slow_F1fast': {'carry_mrna_min': 1.0, 'bit_mrna_min': 1.0,
                      'bit_maturation_min': 1.0, 'a1_mat': 2.0, 'f1_mat': 0.5},
}


def build(prototype, point, config):
    """Build one arm, then assert every requested value by reading it back.

    BUG FIXED HERE: the mismatch configurations rebuild the receiver, and the
    previous version applied the ZENG patch with `m.tail.p[key] = val` where
    key was the TUPLE ('K_A', 1). That inserted a new tuple key into the dict
    while the RHS reads p['K_A'][1], so the K_A_high and K_F_high rows combined
    with A1fast_F1slow / A1slow_F1fast silently ran at NOMINAL thresholds --
    four rows labelled high but computed at nominal. The patch is now applied
    with the same per-index logic as scan_oat.build_receiver, and the readback
    assertion below would fail loudly if it regressed.
    """
    um, ka, kf, n = POINTS[point]
    factors = CONFIGS[config]
    kw = {k: FROZEN.__getattribute__(k) * factors[k]
          for k in TIME_KEYS if k in factors}
    cfg = HbyConfig(**{**asdict(FROZEN), **kw, 'n_A1_gate': n,
                       'receiver_uM_per_au': FROZEN.receiver_uM_per_au * um})
    want_KA = ZENG['K_A'][1] * ka
    want_KF = ZENG['K_F'][1] * kf
    patch = {('K_A', 1): want_KA, ('K_F', 1): want_KF}

    a1 = cfg.carry_maturation_min * factors.get('a1_mat', 1.0)
    f1 = cfg.carry_maturation_min * factors.get('f1_mat', 1.0)
    split = not (np.isclose(a1, cfg.carry_maturation_min)
                 and np.isclose(f1, cfg.carry_maturation_min))
    if split:
        # rebuild through the same helper S.variant uses, so the patch logic is
        # identical to the non-split path
        m = copy.copy(prototype)
        m.tail = SplitCarryReceiver(cfg, a1, f1)
        _apply_zeng_patch(m.tail, patch)
        m.p = m.tail.p
        m.copies_per_au = 602.214076 * m.tail.c.receiver_uM_per_au
    else:
        m = S.variant(prototype, cfg, patch)

    # ---- readback assertion: the requested thresholds must be the live ones
    got_KA = float(m.p['K_A'][1])
    got_KF = float(m.p['K_F'][1])
    assert abs(got_KA - want_KA) < 1e-12, (point, config, got_KA, want_KA)
    assert abs(got_KF - want_KF) < 1e-12, (point, config, got_KF, want_KF)
    assert not any(isinstance(k, tuple) for k in m.p), \
        f'tuple key leaked into the parameter dict: {[k for k in m.p if isinstance(k, tuple)]}'
    assert abs(float(m.tail.c.receiver_uM_per_au)
               - FROZEN.receiver_uM_per_au * um) < 1e-12
    assert float(m.tail.c.n_A1_gate) == n
    return m, dict(
        uM_fold=um, K_A_fold=ka, K_F_fold=kf, n_A1_gate=n,
        timing={k: cfg.__getattribute__(k) for k in TIME_KEYS},
        a1_maturation_min=a1, f1_maturation_min=f1, split_receiver=split,
        readback=dict(K_A1=got_KA, K_F1=got_KF,
                      uM_per_au=float(m.tail.c.receiver_uM_per_au),
                      n_A1_gate=float(m.tail.c.n_A1_gate),
                      kmc_A=float(getattr(m.tail, 'kmc_A', m.tail.kmc)),
                      kmc_F=float(getattr(m.tail, 'kmc_F', m.tail.kmc)),
                      carry_mrna_min=float(m.tail.c.carry_mrna_min),
                      carry_maturation_min=float(m.tail.c.carry_maturation_min),
                      bit_mrna_min=float(m.tail.c.bit_mrna_min),
                      bit_maturation_min=float(m.tail.c.bit_maturation_min),
                      n_tuple_keys=sum(1 for k in m.p if isinstance(k, tuple))))


def _apply_zeng_patch(rec, zeng_patch):
    """Identical logic to scan_oat.build_receiver's patch application."""
    rec.p = dict(rec.p)
    for (key, slot), value in zeng_patch.items():
        if slot is None:
            rec.p[key] = float(value)
        else:
            seq = list(rec.p[key])
            seq[slot] = float(value)
            rec.p[key] = tuple(seq)
    rec.k_complex = rec.q * rec.p['K_D_comp']
    return rec


def self_test_receiver_equivalence():
    """SplitCarryReceiver with equal rates must equal HbyReceiver state-by-state.

    This is the RHS-equivalence test the mismatch configurations depend on: the
    copied rhs and initial_state must be numerically identical to the parent
    when kmc_A == kmc_F == kmc.
    """
    cfg = HbyConfig(**asdict(FROZEN))
    parent = HbyReceiver(cfg)
    c = cfg.carry_maturation_min
    child = SplitCarryReceiver(cfg, c, c)
    y0p, y0c = parent.initial_state(), child.initial_state()
    gap0 = float(np.max(np.abs(y0p - y0c)))
    rng = np.random.default_rng(20260930)
    worst = 0.0
    for _ in range(24):
        y = rng.uniform(0, 3, 17)
        int0 = float(rng.uniform(0, 2))
        pb1 = float(rng.uniform(0, 1))
        dp = np.asarray(parent.rhs(0.0, y, pb1=pb1, int0=int0), dtype=float)
        dc = np.asarray(child.rhs(0.0, y, pb1=pb1, int0=int0), dtype=float)
        worst = max(worst, float(np.max(np.abs(dp - dc))))
    return dict(initial_state_gap=gap0, rhs_gap_max=worst,
                n_states=17, n_probes=24,
                equivalent=bool(gap0 < 1e-12 and worst < 1e-12))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--points', default='ALL')
    ap.add_argument('--configs', default='ALL')
    args = ap.parse_args()
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8')

    out = ROOT / 'results' / ('exprcorner_' + datetime.now().strftime('%Y%m%d_%H%M%S_%f'))
    out.mkdir(parents=True, exist_ok=True)
    S.dump(out / 'status.json', dict(status='RUNNING'))

    print('building Han upstream once', flush=True)
    han = S.HanInput(S.HOURS)
    assert han.p.peak_load_fraction == 0.0
    prototype = S.HZHModel('han', han)

    # ---- (a) RHS equivalence of the locally defined split receiver --------
    eq = self_test_receiver_equivalence()
    print('split-receiver equivalence:', json.dumps(eq), flush=True)
    if not eq['equivalent']:
        S.dump(out / 'status.json', dict(status='FAILED',
                                         error='SplitCarryReceiver is not '
                                               'equivalent to HbyReceiver at '
                                               'equal maturation rates',
                                         equivalence=eq))
        raise RuntimeError(f'SplitCarryReceiver equivalence failed: {eq}')

    # ---- (b) baseline guard ---------------------------------------------
    mb = S.variant(prototype, HbyConfig(**asdict(FROZEN)))
    tb, yb = S.integrate(mb)
    vb, _ = VH.analyse(tb, yb, mb.z, mb.tail)
    checks = dict(full=bool(vb['certified_v1']) is True,
                  reads=int(vb['steady']['reads']) == S.REF['steady_reads'],
                  sequence=str(vb['steady']['sequence']) == S.REF['sequence'],
                  margin=abs(vb['global_min_timing_margin_h']
                             - S.REF['margin_h']) < 1e-6)
    print('baseline guard:', json.dumps(checks), flush=True)
    if not all(checks.values()):
        S.dump(out / 'status.json', dict(status='FAILED', checks=checks))
        raise RuntimeError(f'Baseline guard failed: {checks}')
    S.dump(out / 'guards.json', dict(equivalence=eq, baseline=checks))

    pts = list(POINTS) if args.points == 'ALL' else args.points.split(',')
    cfs = list(CONFIGS) if args.configs == 'ALL' else args.configs.split(',')
    spec = [(p, c) for p in pts for c in cfs]
    if args.points == 'ALL' and args.configs == 'ALL':
        spec = [x for x in spec if x != ('nominal', 'nominal')]  # already known
    print(f'runs: {len(spec)}', flush=True)

    rows = []
    for i, (p, c) in enumerate(spec, 1):
        tic = time.perf_counter()
        m, meta = build(prototype, p, c)
        t, y = S.integrate(m)
        v, sig = VH.analyse(t, y, m.z, m.tail)
        sig['_t'] = t
        reads = v['read_windows']
        pb = PB.per_bit(sig, reads, t)
        bm = v['bit_margins']
        s0, s1 = v['events']['bit0_to_bit1'], v['events']['bit1_to_bit2']
        row = dict(
            label=f'{p} | {c}', point=p, config=c, **meta,
            coupling_readback=meta.pop('readback'),
            split_receiver=meta.pop('split_receiver'),
            counting_passed=bool(v['steady']['increments_mod8']),
            event_causality_passed=bool(all(e['passed'] for e in v['events'].values())),
            full_certified=bool(v['certified_v1']),
            sequence=str(v['steady']['sequence']),
            commitment_global=(None if v['steady']['minimum_commitment'] is None
                               else float(v['steady']['minimum_commitment'])),
            per_bit={b: dict(pb[b], setup_h=bm[b]['min_setup_h'],
                             hold_h=bm[b]['min_hold_h']) for b in PB.BITS},
            s0=dict(reverse=int(s0['reverse_events']), gate=int(s0['gate_events']),
                    flip=int(s0['flip_events']), one_to_one=bool(s0['one_to_one'])),
            s1=dict(reverse=int(s1['reverse_events']), gate=int(s1['gate_events']),
                    flip=int(s1['flip_events']), one_to_one=bool(s1['one_to_one'])),
            j_fwd2_integral=float(np.trapezoid(np.maximum(sig['J_fwd2'], 0), t)),
            j_rev2_integral=float(np.trapezoid(np.maximum(sig['J_rev2'], 0), t)),
            j_rev1_integral=float(np.trapezoid(np.maximum(sig['J_rev1'], 0), t)),
            j_rev0_integral=float(np.trapezoid(np.maximum(sig['J_rev0'], 0), t)),
            g1_peak=float(np.max(sig['g1'])),
            # alpha_Int[1]*g1 is the Int2 TARGET PRODUCTION SOURCE, not the
            # mature free b2_I pool. Both are recorded, under distinct names.
            int2_source_peak=float(m.p['alpha_Int'][1] * np.max(sig['g1'])),
            b2_I_peak=float(np.max(y[B2_I])),
            clock_peak=float(np.max(sig['clock'])),
            runtime_s=time.perf_counter() - tic)
        rows.append(row)
        print(f'[{i}/{len(spec)}] {p:9s} {c:14s} '
              f'cnt={int(row["counting_passed"])} '
              f'ev={int(row["event_causality_passed"])} '
              f'full={int(row["full_certified"])} '
              f'cmt={row["commitment_global"]} '
              f'Jfwd2={row["j_fwd2_integral"]:.3f} Jrev2={row["j_rev2_integral"]:.3f} '
              f'{row["runtime_s"]:.0f}s', flush=True)

    S.dump(out / 'exprcorner.json', dict(
        hours=S.HOURS, sample_min=S.SAMPLE_MIN, rtol=S.RTOL, atol=S.ATOL,
        frozen_config=asdict(FROZEN), points=POINTS, configs=CONFIGS,
        time_keys=TIME_KEYS, rows=rows))
    with (out / 'exprcorner.csv').open('w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh)
        cols = ['label', 'point', 'config', 'uM_fold', 'K_A_fold', 'K_F_fold',
                'n_A1_gate', 'counting_passed', 'event_causality_passed',
                'full_certified', 'sequence', 'commitment_global']
        for b in PB.BITS:
            cols += [f'{b}_commitment', f'{b}_band_margin', f'{b}_hold_h']
        cols += ['s0_rev', 's0_gate', 's0_flip', 's1_rev', 's1_gate', 's1_flip',
                 'j_rev0_integral', 'j_rev1_integral', 'j_fwd2_integral',
                 'j_rev2_integral', 'g1_peak', 'int2_source_peak',
                 'b2_I_peak', 'clock_peak',
                 'runtime_s']
        w.writerow(cols)
        for r in rows:
            vals = [r['label'], r['point'], r['config'], r['uM_fold'],
                    r['K_A_fold'], r['K_F_fold'], r['n_A1_gate'],
                    r['counting_passed'], r['event_causality_passed'],
                    r['full_certified'], r['sequence'], r['commitment_global']]
            for b in PB.BITS:
                p = r['per_bit'][b]
                vals += [p['commitment'], p['band_margin'], p['hold_h']]
            vals += [r['s0']['reverse'], r['s0']['gate'], r['s0']['flip'],
                     r['s1']['reverse'], r['s1']['gate'], r['s1']['flip'],
                     r['j_rev0_integral'], r['j_rev1_integral'],
                     r['j_fwd2_integral'], r['j_rev2_integral'],
                     r['g1_peak'], r['int2_source_peak'], r['b2_I_peak'],
                     r['clock_peak'],
                     round(r['runtime_s'], 3)]
            w.writerow(vals)
    S.dump(out / 'status.json', dict(status='COMPLETED', rows=len(rows)))
    print('OUTPUT', out, flush=True)


if __name__ == '__main__':
    main()
