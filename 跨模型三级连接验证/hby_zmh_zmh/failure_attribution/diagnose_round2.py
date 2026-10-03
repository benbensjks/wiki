"""Round 2: locate the clock_K operating window, not just one passing point.

Round 1/1b established that clock_K 0.40 fails and 0.20 passes, on both A1
autoregulation arms. A single passing point does not say how far it is from the
boundary, so this round brackets it.

Grid
  upper bracket   clock_K in {0.25, 0.30, 0.35} x arms {off, on}
  lower check     clock_K in {0.05, 0.02} x arm {off}
  anchors         clock_K in {0.40, 0.20} x arm {off}   (must reproduce
                  rounds 1/1b exactly, otherwise this scan is invalid)

Every run also records continuous quantities, because a pass/fail boolean hides
the case where the count survives but the margin collapses:
  * per-bit band margin  max(LOW - max, min - HIGH) over the steady read
    windows, including unlabelled ones
  * per-bit in-band occupancy and unlabelled count
  * S2 hold margin, gate dose and duration, RDF2 peak, reverse Hill, v_rev/v_fwd

No existing file is modified; everything is written under
failure_attribution/results/round2_<stamp>/.
"""
from __future__ import annotations

import csv
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np

sys.dont_write_bytecode = True

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT.parent))
sys.path.insert(0, str(ROOT.parent.parent))

import diagnose_round1 as D                                       # noqa: E402
import verify_hzh as V                                            # noqa: E402

UPPER = (0.25, 0.30, 0.35)
UPPER_ARMS = ('off', 'on')
LOWER = (0.05, 0.02)
LOWER_ARMS = ('off',)
ANCHORS = ((0.40, 'off'), (0.20, 'off'))
# previously measured anchor verdicts (rounds 1 and 1b)
ANCHOR_EXPECTED = {
    (0.40, 'off'): dict(counting=False, events=False, s2_flips=1,
                        steady='xxx4567456745674567'),
    (0.20, 'off'): dict(counting=True, events=True, s2_flips=7,
                        steady='1234567012345670123'),
}
BITS = ('S0', 'S1', 'S2')


def band_margins(res):
    """max(LOW - window max, window min - HIGH) per bit, unlabelled included.

    Positive means the whole window sits inside one band; the value is how much
    room is left. Unlabelled windows are NOT skipped -- they are exactly where
    the margin collapses.
    """
    out = {}
    for i, b in enumerate(BITS):
        worst, unlab = None, 0
        for r in res['read_windows'][V.DROP:]:
            if r['labels'][i] is None:
                unlab += 1
            lo, hi = r['extremes'][i]
            d = max(V.LOW - hi, lo - V.HIGH)
            worst = d if worst is None else min(worst, d)
        out[b] = dict(band_margin=worst, n_unlabelled=unlab,
                      n_windows=len(res['read_windows'][V.DROP:]))
    return out


def main():
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8')
    out = ROOT / 'results' / ('round2_' + datetime.now().strftime('%Y%m%d_%H%M%S'))
    out.mkdir(parents=True, exist_ok=False)
    D.dump(out / 'status.json', dict(status='running'))

    print('build shared %g h Han input' % D.HOURS, flush=True)
    han = D.HanInput(D.HOURS)
    prefix = D.HZHModel('han', han)

    tb, yb = D.integrate(prefix)
    bv, _ = V.analyse(tb, yb, prefix.z, prefix.tail)
    guard = dict(certified=bv['certified_v1'],
                 reads=bv['steady']['reads'] == 19,
                 sequence=bv['steady']['sequence'] == '1234567012345670123',
                 margin=abs(bv['global_min_timing_margin_h']
                            - 1.1784609018563117) < 1e-6)
    print('HZH prefix guard:', json.dumps(guard), flush=True)
    if not all(guard.values()):
        D.dump(out / 'status.json', dict(status='FAILED', guard=guard))
        raise RuntimeError(guard)

    # anchors FIRST, so a harness mismatch aborts before the scan is spent
    spec = ([(k, a, 'anchor') for k, a in ANCHORS]
            + [(k, a, 'upper') for k in UPPER for a in UPPER_ARMS]
            + [(k, a, 'lower') for k in LOWER for a in LOWER_ARMS])
    print('runs:', len(spec), flush=True)

    rows = []
    anchor_checks = {}
    for clock_K, arm, block in spec:
        label = f'K{clock_K:.2f}|{arm}'
        model = D.HZZModel(prefix, arm == 'on')
        model.clock_K = clock_K
        t, y = D.integrate(model)
        sig = model.signals(y)
        sig['_p'] = model.p
        if y.min() < -1e-7:
            raise RuntimeError(f'{label}: nonphysical state')
        res = D.evaluate(t, sig)
        row = dict(label=label, arm=arm, block=block, clock_K=clock_K)
        row.update(D.summarise(t, y, sig, res, clock_K, {}))
        row['band_margins'] = band_margins(res)
        row['per_bit'] = res['per_bit']
        row['prefix_gap'] = float(np.max(np.abs(y[:17] - yb[:17])))
        rows.append(row)

        if block == 'anchor':
            want = ANCHOR_EXPECTED[(clock_K, arm)]
            ok = bool(row['counting_passed'] == want['counting']
                      and row['event_causality_passed'] == want['events']
                      and row['s2_flips'] == want['s2_flips']
                      and row['steady'] == want['steady'])
            anchor_checks[label] = dict(
                ok=ok, want=want,
                got=dict(counting=row['counting_passed'],
                         events=row['event_causality_passed'],
                         s2_flips=row['s2_flips'], steady=row['steady']))
            if not ok:
                D.dump(out / 'anchor_checks.json', anchor_checks)
                D.dump(out / 'status.json', dict(
                    status='FAILED', error='anchor did not reproduce',
                    anchors=anchor_checks))
                raise RuntimeError(f'anchor reproduction failed at {label}: '
                                   f'{anchor_checks[label]}')
            print(f'  anchor {label} reproduces rounds 1/1b', flush=True)

        np.savez_compressed(out / f'traj_K{clock_K:.2f}_{arm}.npz',
                            time_h=t, states=y, state_names=np.array(D.NAMES))
        D.dump(out / f'verdict_K{clock_K:.2f}_{arm}.json', res)

        bm = row['band_margins']['S2']
        m = row['mechanism']
        dose_max = (row['gate_dose_range'] or [None, None])[1]
        dose_str = 'None' if dose_max is None else f'{dose_max:.4f}'
        hold = row['per_bit']['S2']['timing']['min_hold_h']
        hold_str = 'None' if hold is None else f'{hold:.3f}'
        print(f"[{label:10s} {block:6s}] cnt={int(row['counting_passed'])} "
              f"ev={int(row['event_causality_passed'])} "
              f"S2fl={row['s2_flips']} dose={dose_str} "
              f"meetBoth={row['gate_segments_meeting_both']} "
              f"RDF2pk={m['RDF2_peak']:.4f} hill={m['reverse_hill_value']:.4f} "
              f"S2 bandM={bm['band_margin']:+.4f} unlab={bm['n_unlabelled']} "
              f"S2hold={hold_str} seq={row['steady']}", flush=True)

    # ---- anchors already verified inline above
    D.dump(out / 'anchor_checks.json', anchor_checks)
    print('anchor checks:', json.dumps(anchor_checks), flush=True)

    # ---- boundary estimate per arm, from the pooled tested points
    pooled = {}
    for r in rows:
        pooled.setdefault(r['arm'], []).append((r['clock_K'],
                                                bool(r['counting_passed']),
                                                bool(r['event_causality_passed'])))
    for arm in ('off', 'on'):
        pooled.setdefault(arm, []).extend(
            [(0.40, False, False), (0.07, True, True)] if arm == 'off'
            else [(0.40, False, False), (0.20, True, True), (0.07, True, True)])
    bounds = {}
    for arm, pts in pooled.items():
        pts = sorted(set(pts))
        passing = sorted(k for k, c, e in pts if c and e)
        failing = sorted(k for k, c, e in pts if not (c and e))
        b = dict(tested=sorted({k for k, _, _ in pts}),
                 passing=passing, failing=failing,
                 any_passing=bool(passing))
        if passing:
            hi, lo = max(passing), min(passing)
            above = [k for k in failing if k > hi]
            below = [k for k in failing if k < lo]
            b['upper_boundary_between'] = [hi, min(above)] if above else [hi, None]
            b['lower_boundary_between'] = [max(below), lo] if below else [None, lo]
        bounds[arm] = b

    D.dump(out / 'round2_summary.json', dict(
        hours=D.HOURS, solver=D.SOLVER, prefix_guard=guard,
        grid=dict(upper=UPPER, upper_arms=UPPER_ARMS, lower=LOWER,
                  lower_arms=LOWER_ARMS,
                  anchors=[dict(clock_K=k, arm=a) for k, a in ANCHORS]),
        anchor_checks=anchor_checks, boundaries=bounds,
        note=('clock_K is the only varied quantity; both A1 autoregulation arms '
              'are covered above the lower check. A pass here means counting AND '
              'event causality; the two are also reported separately.'),
        rows=rows))

    cols = ['block', 'arm', 'clock_K', 'counting', 'event', 'full', 's2_flips',
            'steady', 'gate_dose_max', 'gate_duration_max', 'meets_both',
            'g1_peak', 'clock_peak', 'I2_peak', 'RDF2_peak', 'T2_min', 'pb2_min',
            'i_times_r', 'reverse_hill', 'v_rev', 'v_fwd', 'vr_over_vf',
            'S0_band_margin', 'S1_band_margin', 'S2_band_margin',
            'S0_commit', 'S1_commit', 'S2_commit',
            'S2_unlabelled', 'S2_hold_h', 'J_fwd2_total', 'J_rev2_total',
            'prefix_gap']
    with (out / 'round2_all.csv').open('w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for r in rows:
            m, bm = r['mechanism'], r['band_margins']
            dur = r['gate_duration_range'] or [None, None]
            dose = r['gate_dose_range'] or [None, None]
            w.writerow([
                r['block'], r['arm'], r['clock_K'], r['counting_passed'],
                r['event_causality_passed'], r['diagnostic_passed'],
                r['s2_flips'], r['steady'], dose[1], dur[1],
                r['gate_segments_meeting_both'], r['g1_peak'], r['clock_peak'],
                m['I2_peak'], m['RDF2_peak'], m['T2_min'], m['pb2_min'],
                m['i_times_r_peak'], m['reverse_hill_value'],
                m['v_reverse_at_peak_per_h'], m['v_forward_at_peak_per_h'],
                m['reverse_over_forward'],
                bm['S0']['band_margin'], bm['S1']['band_margin'],
                bm['S2']['band_margin'],
                r['per_bit']['S0']['minimum_commitment'],
                r['per_bit']['S1']['minimum_commitment'],
                r['per_bit']['S2']['minimum_commitment'],
                bm['S2']['n_unlabelled'],
                r['per_bit']['S2']['timing']['min_hold_h'],
                r['J_fwd2_total'], r['J_rev2_total'], r['prefix_gap']])

    D.dump(out / 'status.json', dict(status='COMPLETED', rows=len(rows),
                                     boundaries=bounds))
    D.dump(out / 'SHA256SUMS.json',
           {str(p.relative_to(out)): D.sha256(p)
            for p in sorted(out.rglob('*')) if p.is_file()})
    print('BOUNDARIES', json.dumps(bounds), flush=True)
    print('OUTPUT', out, flush=True)


if __name__ == '__main__':
    main()
