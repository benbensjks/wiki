"""Clock-compensation grid for the upper uM_per_au boundary.

Motivation
----------
`uM_per_au` enters the model through ONE channel only: `copies_per_au`
(`model_hzh.py:49,81`) scales bit0's C31 mRNA source, hence the a.u. level of
mature free Int0, hence the second-stage clock Hill term.

If  Int0 -> Int0/f  and  clock_K -> clock_K/f  then

    H(Int0/f ; K/f, n) = H(Int0 ; K, n)

so the clock term should be restored exactly. This grid tests whether that
restores counting at the upper boundary.

IMPORTANT CAVEATS (must travel with every reading of this experiment)
--------------------------------------------------------------------
1. The identity above is exact ONLY IF I0(t) really equals I0_base(t)/f.
   It does NOT follow from scaling uM_per_au: bit0 keeps fixed K_D_int[0],
   K_inh and k_complex here, and its RDF/complex dynamics are not rescaled, so
   I0 is only approximately proportional to 1/f. Both deviations are therefore
   MEASURED, not assumed:
     * int0_scale_deviation          = max|f*I0_new(t) - I0_base(t)|
     * clock_pointwise_dev_vs_baseline = max|clock_new(t) - clock_base(t)|
   "the clock was exactly restored" must never be written as a premise.
2. Even if compensation restores counting, that shows ADJUSTING THE CLOCK GATE
   IS SUFFICIENT TO RESCUE THAT POINT. It does NOT by itself show that the
   original failure was mainly caused by the clock.
3. uM_per_au rescales bit0's whole concentration axis, so a failure of
   compensation does NOT mean the clock is unimportant; it means clock
   rescaling alone cannot cancel the whole uM effect.

Note: `clock_K_au` and `clock_scale` are exactly degenerate in this model
(verified: max|dy| = 0.000e+00), so compensating through either is equivalent.

Design
------
uM fold in {1.30, 1.35, 1.40, 1.45} x clock_K in
  { K0 (=0.3, uncompensated), K0/f (compensated),
    0.9*K0/f, 1.1*K0/f }   -> 16 runs.
The uncompensated column cross-checks the dense scan.

Reported per run: per-bit commitment and band margin (S0/S1/S2), both event
chains, clock Hill waveform statistics, g1 and Int2 waveform statistics, and
the S2 forward/reverse flux integrals.
"""
from __future__ import annotations

import argparse, csv, json, sys, time
from dataclasses import asdict
from datetime import datetime
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

import scan_oat as S                                              # noqa: E402
import verify_hzh as VH                                           # noqa: E402
import perbit_margin as PB                                        # noqa: E402
from hybrid_model import HbyConfig                                # noqa: E402
from model_hzh import IDX as IDX_HZH                              # noqa: E402

FROZEN = S.FROZEN
UMS = (1.30, 1.35, 1.40, 1.45)
B2_I = IDX_HZH['b2_I']


def wave(t, x, reads):
    """Peak, width at half max, per-cycle integral and far-dwell for a signal."""
    peaks, widths, integrals, resid = [], [], [], []
    for r in reads:
        m = (t >= r['cycle_start_h']) & (t <= r['cycle_end_h'])
        if m.sum() < 3:
            continue
        tt, yy = t[m], np.asarray(x[m], dtype=float)
        pk = float(yy.max())
        peaks.append(pk)
        if pk > 0:
            above = tt[yy >= pk / 2]
            widths.append(float(above.max() - above.min()) if above.size > 1 else 0.0)
        integrals.append(float(np.trapezoid(np.maximum(yy, 0.0), tt)))
        k = max(1, int(0.2 * tt.size))
        resid.append(float(np.median(yy[np.argsort(yy)[:k]])))
    return dict(peak_max=max(peaks) if peaks else None,
                peak_min=min(peaks) if peaks else None,
                width_max=max(widths) if widths else None,
                integral_min=min(integrals) if integrals else None,
                integral_max=max(integrals) if integrals else None,
                far_dwell_median=float(np.median(resid)) if resid else None)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--ums', default=','.join(str(u) for u in UMS))
    args = ap.parse_args()
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8')

    out = ROOT / 'results' / ('clkcomp_' + datetime.now().strftime('%Y%m%d_%H%M%S_%f'))
    out.mkdir(parents=True, exist_ok=True)
    S.dump(out / 'status.json', dict(status='RUNNING'))

    print('building Han upstream once', flush=True)
    han = S.HanInput(S.HOURS)
    assert han.p.peak_load_fraction == 0.0
    prototype = S.HZHModel('han', han)

    # ---- baseline guard, and the reference waveforms for the identity check
    mb = S.variant(prototype, HbyConfig(**asdict(FROZEN)))
    tb, yb = S.integrate(mb)
    vb, sigb = VH.analyse(tb, yb, mb.z, mb.tail)
    checks = dict(full=bool(vb['certified_v1']) is True,
                  reads=int(vb['steady']['reads']) == S.REF['steady_reads'],
                  sequence=str(vb['steady']['sequence']) == S.REF['sequence'],
                  margin=abs(vb['global_min_timing_margin_h']
                             - S.REF['margin_h']) < 1e-6)
    print('baseline guard:', json.dumps(checks), flush=True)
    if not all(checks.values()):
        S.dump(out / 'status.json', dict(status='FAILED', checks=checks))
        raise RuntimeError(f'Baseline guard failed: {checks}')
    int0_base = np.asarray(sigb['int0'], dtype=float)
    clock_base = np.asarray(sigb['clock'], dtype=float)
    S.dump(out / 'baseline_guard.json', checks)

    rows = []
    n_done = 0
    for f in [float(x) for x in args.ums.split(',')]:
        k0 = FROZEN.clock_K_au
        levels = [('none', k0), ('comp', k0 / f),
                  ('comp_x0.9', 0.9 * k0 / f), ('comp_x1.1', 1.1 * k0 / f)]
        for tag, k in levels:
            n_done += 1
            label = f'uM=x{f:.2f} clock_K={tag}'
            tic = time.perf_counter()
            cfg = HbyConfig(**{**asdict(FROZEN),
                               'receiver_uM_per_au': FROZEN.receiver_uM_per_au * f,
                               'clock_K_au': k})
            m = S.variant(prototype, cfg)
            t, y = S.integrate(m)
            v, sig = VH.analyse(t, y, m.z, m.tail)
            sig['_t'] = t
            reads = v['read_windows']
            pb = PB.per_bit(sig, reads, t)
            bm = v['bit_margins']
            s0, s1 = v['events']['bit0_to_bit1'], v['events']['bit1_to_bit2']
            int2_src = m.p['alpha_Int'][1] * sig['g1']
            b2i = np.asarray(y[B2_I], dtype=float)
            int0_new = np.asarray(sig['int0'], dtype=float)
            clock_new = np.asarray(sig['clock'], dtype=float)
            # The algebraic identity H(I0/f; K/f, n) = H(I0; K, n) is exact only
            # if I0(t) really is I0_base(t)/f. bit0 keeps fixed K_D_int[0],
            # K_inh and k_complex here, so that scaling is NOT guaranteed.
            # Both deviations are measured rather than assumed.
            scale_dev = float(np.max(np.abs(f * int0_new - int0_base)))
            clock_dev = float(np.max(np.abs(clock_new - clock_base)))
            row = dict(
                label=label, uM_fold=f, level=tag,
                clock_K_requested=k,
                clock_K_readback=float(m.tail.c.clock_K_au),
                clock_K_effective=float(m.tail.c.clock_K_au / m.tail.c.clock_scale),
                counting_passed=bool(v['steady']['increments_mod8']),
                event_causality_passed=bool(all(e['passed']
                                                for e in v['events'].values())),
                full_certified=bool(v['certified_v1']),
                sequence=str(v['steady']['sequence']),
                commitment_global=(None if v['steady']['minimum_commitment'] is None
                                   else float(v['steady']['minimum_commitment'])),
                per_bit={b: dict(pb[b], setup_h=bm[b]['min_setup_h'],
                                 hold_h=bm[b]['min_hold_h']) for b in PB.BITS},
                s0=dict(reverse=int(s0['reverse_events']),
                        gate=int(s0['gate_events']), flip=int(s0['flip_events']),
                        one_to_one=bool(s0['one_to_one'])),
                s1=dict(reverse=int(s1['reverse_events']),
                        gate=int(s1['gate_events']), flip=int(s1['flip_events']),
                        one_to_one=bool(s1['one_to_one'])),
                clock_wave=wave(t, sig['clock'], reads),
                g1_wave=wave(t, sig['g1'], reads),
                int2_source_wave=wave(t, int2_src, reads),
                b2_I_wave=wave(t, b2i, reads),
                b2_I_peak=float(np.max(b2i)),
                int2_source_peak=float(np.max(int2_src)),
                int0_scale_deviation=scale_dev,
                clock_pointwise_dev_vs_baseline=clock_dev,
                clock_range=[float(np.min(sig['clock'])), float(np.max(sig['clock']))],
                j_fwd2_integral=float(np.trapezoid(np.maximum(sig['J_fwd2'], 0), t)),
                j_rev2_integral=float(np.trapezoid(np.maximum(sig['J_rev2'], 0), t)),
                int0_range=[float(np.min(sig['int0'])), float(np.max(sig['int0']))],
                runtime_s=time.perf_counter() - tic)
            rows.append(row)
            cw = row['clock_wave']
            print(f'[{n_done}] uM={f:.2f} {tag:10s} K={k:.4f} '
                  f'cnt={int(row["counting_passed"])} '
                  f'ev={int(row["event_causality_passed"])} '
                  f'full={int(row["full_certified"])} '
                  f'clockPk={cw["peak_max"]:.4f} '
                  f'|clockD|={clock_dev:.3e} '
                  f'|fI0-I0b|={scale_dev:.3e} '
                  f'b2Ipk={row["b2_I_peak"]:.3f} '
                  f'cmt S0/S1/S2={pb["S0"]["commitment"]}/{pb["S1"]["commitment"]}/'
                  f'{pb["S2"]["commitment"]} {row["runtime_s"]:.0f}s',
                  flush=True)

    S.dump(out / 'clkcomp.json', dict(
        hours=S.HOURS, sample_min=S.SAMPLE_MIN, rtol=S.RTOL, atol=S.ATOL,
        frozen_config=asdict(FROZEN),
        design=dict(ums=list(UMS),
                    levels=['none', 'comp', 'comp_x0.9', 'comp_x1.1'],
                    compensation='clock_K_au -> clock_K_au / uM_fold'),
        identity='H(Int0/f ; K/f, n) = H(Int0 ; K, n)',
        caveat=('(1) the identity needs I0_new(t) = I0_base(t)/f, which does NOT '
                'follow from rescaling uM_per_au because bit0 keeps fixed '
                'K_D_int[0]/K_inh/k_complex; int0_scale_deviation and '
                'clock_pointwise_dev_vs_baseline measure the departure. '
                '(2) restoring counting shows clock adjustment is SUFFICIENT to '
                'rescue the point, not that the clock was the original cause. '
                '(3) a failure of compensation does not mean the clock is '
                'unimportant.'),
        rows=rows))

    with (out / 'clkcomp.csv').open('w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh)
        cols = ['label', 'uM_fold', 'level', 'clock_K_requested',
                'clock_K_readback', 'counting_passed', 'event_causality_passed',
                'full_certified', 'sequence', 'commitment_global']
        for b in PB.BITS:
            cols += [f'{b}_commitment', f'{b}_band_margin', f'{b}_hold_h']
        cols += ['s0_rev', 's0_gate', 's0_flip', 's1_rev', 's1_gate', 's1_flip',
                 'clock_peak_max', 'clock_width_max', 'clock_integral_min',
                 'clock_lo', 'clock_hi', 'int0_lo', 'int0_hi',
                 'clock_pointwise_dev_vs_baseline', 'int0_scale_deviation',
                 'g1_peak_max', 'g1_integral_min',
                 'int2_source_peak_max', 'b2_I_peak',
                 'j_fwd2_integral', 'j_rev2_integral', 'runtime_s']
        w.writerow(cols)
        for r in rows:
            vals = [r['label'], r['uM_fold'], r['level'], r['clock_K_requested'],
                    r['clock_K_readback'], r['counting_passed'],
                    r['event_causality_passed'], r['full_certified'],
                    r['sequence'], r['commitment_global']]
            for b in PB.BITS:
                p = r['per_bit'][b]
                vals += [p['commitment'], p['band_margin'], p['n_unlabelled'],
                         p['hold_h']]
            cw, gw = r['clock_wave'], r['g1_wave']
            iw, bw = r['int2_source_wave'], r['b2_I_wave']
            vals += [r['s0']['reverse'], r['s0']['gate'], r['s0']['flip'],
                     r['s1']['reverse'], r['s1']['gate'], r['s1']['flip'],
                     cw['peak_max'], cw['width_max'], cw['integral_min'],
                     r['clock_range'][0], r['clock_range'][1],
                     r['int0_range'][0], r['int0_range'][1],
                     r['clock_pointwise_dev_vs_baseline'],
                     r['int0_scale_deviation'],
                     gw['peak_max'], gw['integral_min'],
                     iw['peak_max'], bw['peak_max'],
                     r['j_fwd2_integral'], r['j_rev2_integral'],
                     round(r['runtime_s'], 3)]
            w.writerow(vals)
    S.dump(out / 'status.json', dict(status='COMPLETED', rows=len(rows)))
    print('OUTPUT', out, flush=True)


if __name__ == '__main__':
    main()
