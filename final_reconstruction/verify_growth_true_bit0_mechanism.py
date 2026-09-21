"""Mechanism and parity audit for the robust add_growth=True bit0 candidate."""
from __future__ import annotations

import json
from dataclasses import asdict

import numpy as np

from model import Extension, Model, ROOT, ZENG, act
from rescan800_readwindow import single_bit_verdict
from verify_bit0_part2 import CAND, sim
from verify_bit0_part3 import consistent_bit0_init
from verify_twobit_causal import _read_windows


CANDIDATE = dict(uM_per_au=0.8, maturation_half_life_min=14.0,
                 complex_on_au_inv_h=1.0, complex_off_h=3.0,
                 add_growth=True)


def signals(model, sol):
    S, I, R, C = sol.y[16], sol.y[8], sol.y[14], sol.y[15]
    flux = np.asarray([model.flux(model._expand_bit0(sol.y[:, k]))
                       for k in range(sol.y.shape[1])])
    vf = ZENG['k_fwd'] * np.asarray([act(x, ZENG['K_D_int'][0], 2) for x in I]) * \
        ZENG['K_inh'] / (ZENG['K_inh'] + np.maximum(R, 0.0))
    vr = ZENG['k_rev'] * np.asarray([act(x, model.K_complex, 2) for x in C])
    return flux, vf * (1 - S), vr * S


def transition_audit(model, sol):
    t, S, R, C = sol.t, sol.y[16], sol.y[14], sol.y[15]
    flux, jf, jr = signals(model, sol)
    reads = _read_windows(t, flux, S, S)
    rows = []
    for prev, cur in zip(reads[:-1], reads[1:]):
        a = int(np.searchsorted(t, prev['trough_h']))
        b = int(np.searchsorted(t, cur['trough_h']))
        before, after = prev['bit0']['label'], cur['bit0']['label']
        if 'x' in (before, after) or before == after or t[a] < 60:
            continue
        rows.append(dict(direction=f'{before}->{after}', start_h=float(t[a]), end_h=float(t[b]),
                         R_start=float(R[a]), C_peak=float(C[a:b + 1].max()),
                         forward_peak=float(jf[a:b + 1].max()),
                         reverse_peak=float(jr[a:b + 1].max()),
                         forward_integral=float(np.trapezoid(jf[a:b + 1], t[a:b + 1])),
                         reverse_integral=float(np.trapezoid(jr[a:b + 1], t[a:b + 1])),
                         overlap_integral=float(np.trapezoid(np.minimum(jf[a:b + 1], jr[a:b + 1]),
                                                             t[a:b + 1]))))
    grouped = {}
    for direction in ('0->1', '1->0'):
        g = [r for r in rows if r['direction'] == direction]
        grouped[direction] = ({k: float(np.median([r[k] for r in g]))
                               for k in ('R_start', 'C_peak', 'forward_peak', 'reverse_peak',
                                         'forward_integral', 'reverse_integral', 'overlap_integral')}
                              if g else {})
        grouped[direction]['events'] = len(g)
    return reads, rows, grouped


def parity_audit(model, low, high, hours=180.0):
    runs = {}
    for name, state in (('high', high), ('low', low)):
        runs[name] = sim(model, consistent_bit0_init(model, state), 0.0, hours)
    out = {}
    for name, sol in runs.items():
        flux, _, _ = signals(model, sol)
        reads = _read_windows(sol.t, flux, sol.y[16], sol.y[16])
        out[name] = dict(initial_S=high if name == 'high' else low,
                         sequence=''.join(r['bit0']['label'] for r in reads),
                         steady=single_bit_verdict(reads, 4))
    gap = np.abs(runs['high'].y[16] - runs['low'].y[16])
    opposite = all((a != 'x' and b != 'x' and a != b)
                   for a, b in zip(out['high']['sequence'][4:], out['low']['sequence'][4:]))
    return dict(branches=out, opposite_after_drop4=opposite,
                max_gap_last_40h=float(gap[runs['high'].t >= hours - 40].max()))


def orbit_branch_audit(model, reference_sol, reference_reads, hours=180.0):
    """Start from the two actual period-2 orbit branches at one clock phase."""
    pair = None
    for a, b in zip(reference_reads[:-1], reference_reads[1:]):
        if (a['trough_h'] >= 180 and 'x' not in (a['bit0']['label'], b['bit0']['label'])
                and a['bit0']['label'] != b['bit0']['label']):
            pair = (a, b)
            break
    if pair is None:
        raise RuntimeError('could not locate opposite late orbit branches')
    ia = int(np.argmin(np.abs(reference_sol.t - pair[0]['trough_h'])))
    ib = int(np.argmin(np.abs(reference_sol.t - pair[1]['trough_h'])))
    ya = reference_sol.y[:, ia].copy()
    yb = reference_sol.y[:, ib].copy()
    upstream_mismatch_before_reset = float(np.max(np.abs(ya[:7] - yb[:7])))
    yb[:7] = ya[:7]
    runs = {'branch_A': sim(model, ya, 0.0, hours),
            'branch_B': sim(model, yb, 0.0, hours)}
    out = {}
    for name, sol in runs.items():
        flux, _, _ = signals(model, sol)
        reads = _read_windows(sol.t, flux, sol.y[16], sol.y[16])
        out[name] = dict(initial_S=float(sol.y[16, 0]),
                         sequence=''.join(r['bit0']['label'] for r in reads),
                         steady=single_bit_verdict(reads, 4))
    gap = np.abs(runs['branch_A'].y[16] - runs['branch_B'].y[16])
    sa, sb = out['branch_A']['sequence'][4:], out['branch_B']['sequence'][4:]
    opposite = bool(sa and len(sa) == len(sb) and
                    all(a != 'x' and b != 'x' and a != b for a, b in zip(sa, sb)))
    return dict(source_troughs_h=[pair[0]['trough_h'], pair[1]['trough_h']],
                source_labels=[pair[0]['bit0']['label'], pair[1]['bit0']['label']],
                upstream_mismatch_before_phase_reset=upstream_mismatch_before_reset,
                branches=out, opposite_after_drop4=opposite,
                max_gap_last_40h=float(gap[runs['branch_A'].t >= hours - 40].max()))


def kick_audit(model, reference, kick_h=100.0, hours=400.0):
    y0 = model.initial_state_bit0(cold=False)
    pre = sim(model, y0, 0.0, kick_h)
    ykick = pre.y[:, -1].copy(); before = float(ykick[16]); ykick[16] = 1 - ykick[16]
    post = sim(model, ykick, kick_h, hours)
    t = np.r_[pre.t, post.t[1:]]
    y = np.c_[pre.y, post.y[:, 1:]]
    class Combined: pass
    combined = Combined(); combined.t = t; combined.y = y
    flux, _, _ = signals(model, combined)
    reads = _read_windows(t, flux, y[16], y[16])
    kicked = ''.join(r['bit0']['label'] for r in reads)
    refseq = ''.join(r['bit0']['label'] for r in reference)
    tail_k, tail_r = kicked[-6:], refseq[-6:]
    flipped = all(a != 'x' and b != 'x' and a != b for a, b in zip(tail_k, tail_r))
    return dict(kick_h=kick_h, S_before=before, S_after=float(ykick[16]),
                reference_tail=tail_r, kicked_tail=tail_k,
                parity_permanently_flipped=flipped)


def main():
    cfg = asdict(Extension()); cfg.update(CAND); cfg.update(CANDIDATE)
    model = Model(Extension(**cfg))
    sol = model.simulate_bit0(hours=400.0, sample_min=2.0, max_step_min=2.0)
    reads, transitions, grouped = transition_audit(model, sol)
    steady = reads[4:]
    lows = [r['bit0']['median'] for r in steady if r['bit0']['label'] == '0']
    highs = [r['bit0']['median'] for r in steady if r['bit0']['label'] == '1']
    low, high = float(np.median(lows)), float(np.median(highs))
    result = dict(parameters=CANDIDATE,
                  readout=single_bit_verdict(reads, 4),
                  self_consistent_levels=dict(low=low, high=high),
                  transition_events=transitions, transition_medians=grouped,
                  quasi_steady_parity=parity_audit(model, low, high),
                  orbit_branch_parity=orbit_branch_audit(model, sol, reads),
                  kick=kick_audit(model, reads))
    out = ROOT / 'bit0_results' / 'verification' / 'growth_true_bit0_mechanism.json'
    out.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(dict(readout=result['readout'], levels=result['self_consistent_levels'],
                          transitions=grouped,
                          quasi_steady_parity=result['quasi_steady_parity'],
                          orbit_branch_parity=result['orbit_branch_parity'],
                          kick=result['kick']),
                     ensure_ascii=False, indent=2))
    print(f'wrote {out}')


if __name__ == '__main__':
    main()
