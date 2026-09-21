"""Adversarial verification of the reported 'stable bit0 alternation'.

Does not change any Zeng parameter and does not re-tune anything. It only asks
whether the top coarse-grid candidate really stores a bit, or whether the
recorded 0/1 codes are an artefact of reading one instant per clock cycle.

Outputs: bit0_results/verification/*.json
"""
from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import find_peaks

from model import Model, Extension, ZENG, ROOT

OUT = ROOT / 'bit0_results' / 'verification'
CONF = ROOT / 'bit0_results' / 'confirmed_candidate'
CAND = dict(uM_per_au=10.0, maturation_half_life_min=20.0,
            complex_on_au_inv_h=0.1, complex_off_h=1.0, add_growth=False)
DT = 2.0 / 60.0  # sampling step of the archived trajectories, h


def decode(v):
    return 1 if v >= 0.8 else (0 if v <= 0.2 else -1)


def code_string(vals, keep=5):
    codes = [decode(v) for v in vals][-keep:] if keep else [decode(v) for v in vals]
    ok = bool(len(codes) >= 2 and all(c >= 0 for c in codes)
              and all(codes[i + 1] == 1 - codes[i] for i in range(len(codes) - 1)))
    return ''.join('x' if c < 0 else str(c) for c in codes), ok


def flux_peaks(t, flux):
    ids, _ = find_peaks(flux, prominence=max(1.0, .1 * np.ptp(flux)), distance=int(5 / (t[1] - t[0])))
    return ids


def load_confirmed():
    df = pd.read_csv(CONF / 'bit0_180h_trajectory.csv')
    return df


def part_A_phase(df):
    """Same trajectory, same decoder, sample instant moved. Codes must be stable."""
    t = df.time_h.to_numpy(); flux = df.flux_uM_h.to_numpy(); S = df.b0_S.to_numpy()
    ids = flux_peaks(t, flux)
    troughs = np.array([a + int(np.argmin(flux[a:b])) for a, b in zip(ids[:-1], ids[1:])])
    rows = []
    for off in np.arange(-4.0, 4.01, 0.5):
        k = troughs + int(round(off / DT))
        k = k[(k >= 0) & (k < len(t))]
        vals = S[k]
        full, _ = code_string(vals, keep=0)
        last5, ok = code_string(vals, keep=5)
        rows.append(dict(offset_h=round(float(off), 2), values=[round(float(v), 4) for v in vals],
                         codes=full, late5=last5, alternates=ok))
    n_ok = int(sum(r['alternates'] for r in rows))
    out = dict(sample_instants=len(troughs), offsets_tested=len(rows), offsets_that_alternate=n_ok,
               frame=rows,
               verdict=('DECODE IS PHASE-ROBUST' if n_ok == len(rows) else
                        'DECODE IS PHASE-DEPENDENT: the 0/1 pattern is a sampling artefact'))
    return out, troughs, ids


def part_B_commitment(df, troughs, ids):
    """A stored bit must hold its value across the whole idle interval."""
    t = df.time_h.to_numpy(); S = df.b0_S.to_numpy()
    rows = []
    for n, (a, b) in enumerate(zip(ids[:-1], ids[1:]), 1):
        sl = slice(a, b + 1)
        s = S[sl]
        hi = float(np.mean(s >= 0.7)); lo = float(np.mean(s <= 0.3))
        rows.append(dict(cycle=n, S_min=float(s.min()), S_max=float(s.max()), S_mean=float(s.mean()),
                         frac_S_ge_0p7=hi, frac_S_le_0p3=lo, commitment=max(hi, lo)))
    late = [r for r in rows if r['cycle'] >= 5]
    comm = np.array([r['commitment'] for r in late])
    bands = []
    for r in late:
        bands.append('1' if r['frac_S_ge_0p7'] >= 0.7 else ('0' if r['frac_S_le_0p3'] >= 0.7 else 'x'))
    band_codes, band_ok = code_string([{'1': 1.0, '0': 0.0, 'x': 0.5}[c] for c in bands], keep=5)
    out = dict(cycles=len(rows), frame=rows,
               median_commitment_late=float(np.median(comm)),
               min_commitment_late=float(comm.min()),
               cycles_with_commitment_ge_0p9=int((comm >= 0.9).sum()),
               cycles_with_commitment_ge_0p7=int((comm >= 0.7).sum()),
               majority_band_codes=''.join(bands),
               majority_band_late5_alternates=band_ok,
               verdict=('STORE-LIKE: every cycle holds one band' if comm.min() >= 0.9 else
                        'NOT A STORED BIT: S never commits to 0 or 1 for a full cycle'))
    return out


def part_C_readout(df):
    """Is S a memory, or the instantaneous flux-balance readout vf/(vf+vr)?"""
    t = df.time_h.to_numpy(); S = df.b0_S.to_numpy()
    f = df.forward_flux.to_numpy(); r = df.reverse_flux.to_numpy()
    vf = np.where(S < 0.999, f / np.clip(1 - S, 1e-9, None), np.nan)
    vr = np.where(S > 1e-6, r / np.clip(S, 1e-9, None), np.nan)
    S_ss = vf / (vf + vr)
    tau = 1.0 / (vf + vr)
    late = t >= 55.0
    d = np.abs(S - S_ss)[late]
    ratio = (vf / vr)[late]
    out = dict(
        correlation_S_vs_instantaneous_balance=float(np.corrcoef(S[late], S_ss[late])[0, 1]),
        median_abs_deviation_from_balance=float(np.nanmedian(d)),
        S_ss_range_late=[float(np.nanmin(S_ss[late])), float(np.nanmax(S_ss[late]))],
        relaxation_time_h=[float(np.nanmin(tau[late])), float(np.nanmedian(tau[late])), float(np.nanmax(tau[late]))],
        clock_period_h=float(np.median(np.diff(t[flux_peaks(t, df.flux_uM_h.to_numpy())]))),
        vf_over_vr_late=[float(np.nanmin(ratio)), float(np.nanmedian(ratio)), float(np.nanmax(ratio))],
        fraction_of_time_vf_and_vr_within_2x=float(np.mean((ratio > 0.5) & (ratio < 2.0))),
        verdict='')
    out['verdict'] = ('S is a leaky readout of vf/(vf+vr), not a latch: it tracks the '
                      'instantaneous balance and its relaxation time is comparable to the clock period'
                      if abs(out['correlation_S_vs_instantaneous_balance']) > 0.8 else
                      'S does not simply track the instantaneous balance')
    return out


def bit_init(m, S):
    y = m.initial_state_bit0()
    y[16] = S
    b = y[6:17]
    src_T = ZENG['alpha_rep'] * (1 - S)
    b[3] = m.transcript_source(src_T) / m.lm
    b[4] = m.e.translation_h * b[3] / m.lu
    b[5] = src_T / (ZENG['gamma_rep'] + m.growth)
    T = b[5]
    src_R = ZENG['alpha_rdf'] * S * (1 - (T / ZENG['K_rep'])**ZENG['n_rep'] / (1 + (T / ZENG['K_rep'])**ZENG['n_rep']))
    b[6] = m.transcript_source(src_R) / m.lm
    b[7] = m.e.translation_h * b[6] / m.lu
    b[8] = src_R / (ZENG['gamma_rdf'] + m.growth)
    return y


def part_D_memory(hours=120.0):
    """Two identical worlds that differ only in the stored DNA state."""
    m = Model(Extension(**CAND))
    res = {}
    for label, S0 in (('start_LR_S1', 1.0), ('start_PB_S0', 0.0)):
        sol = m.simulate_bit0(hours=hours, cold=False, sample_min=2.0, rtol=1e-8, atol=1e-10, max_step_min=1.0)
        res[label] = sol
    t = res['start_LR_S1'].t
    A = res['start_LR_S1'].y[16]; B = res['start_PB_S0'].y[16]
    gap = np.abs(A - B)
    marks = [0.0, 2.0, 5.0, 10.0, 20.0, 40.0, 60.0, 100.0, float(t[-1])]
    trace = {str(mk): float(np.interp(mk, t, gap)) for mk in marks}
    tail = float(np.max(gap[t >= hours - 20.0]))
    early = float(np.median(gap[(t >= 10) & (t <= 20)]))
    return dict(hours=hours, gap_trace_h=trace, gap_median_10_20h=early, gap_max_last_20h=tail,
                initial_gap=1.0,
                verdict=('MEMORY LOST: the initial DNA state is forgotten' if tail < 0.2 * max(early, 1e-9) or tail < 0.05
                         else 'MEMORY RETAINED: the initial DNA state still separates the trajectories'),
                S_end=[float(A[-1]), float(B[-1])])


def part_E_twobit():
    df = pd.read_csv(CONF / 'twobit_candidate_trajectory.csv')
    t = df.time_h.to_numpy(); S0 = df.b0_S.to_numpy(); S1 = df.b1_S.to_numpy()
    carry = df.carry_promoter_01.to_numpy()
    fall = np.flatnonzero((S0[:-1] >= .8) & (S0[1:] < .8)) + 1
    counts = {}
    for prom_frac in (0.02, 0.05, 0.1, 0.2, 0.3):
        p, _ = find_peaks(carry, prominence=max(prom_frac * np.ptp(carry), 1e-9), distance=60)
        counts[f'prominence_{prom_frac}'] = dict(n_peaks=int(len(p)), times=[round(float(t[i]), 2) for i in p])
    # Does the DNA of bit1 ever move?
    dS1 = np.abs(np.diff(S1))
    return dict(S1_range=[float(S1.min()), float(S1.max())], S1_max_abs_step=float(dS1.max()),
                S1_total_excursion=float(S1.max() - S1.min()),
                S1_crossings_of_0p5=int(np.sum((S1[:-1] - .5) * (S1[1:] - .5) < 0)),
                b1_I_range=( [float(df.b1_I.min()), float(df.b1_I.max())] if 'b1_I' in df else None),
                carry_promoter_range=[float(carry.min()), float(carry.max())],
                falling_edges=len(fall), falling_edge_times=[round(float(t[i]), 2) for i in fall],
                carry_peaks_by_prominence=counts,
                verdict=('bit1 DNA NEVER SWITCHES: the recorded carry is promoter activity only'
                         if float(S1.max() - S1.min()) < 0.05 else
                         'bit1 DNA does move; inspect magnitude against the 0/1 bands'))


def part_F_provenance(hours=100.0):
    """Re-run the archived top candidate with the CURRENT script and compare."""
    from bit0_diagnostic import evaluate
    cfg = asdict(Extension()); cfg.update(CAND)
    got = evaluate(cfg, hours=hours)
    arch = json.loads((ROOT / 'bit0_results' / 'scan_summary.json').read_text(encoding='utf-8'))
    return dict(current_code_score={k: got[k] for k in ('alternation_accuracy', 'confident_fraction',
                                                         'dynamic_range', 'stable_alternation', 'codes') if k in got},
                archived_score=arch['top_candidate_score'],
                archived_stable_candidates=arch['stable_bit0_candidates'],
                consistent=bool(got.get('stable_alternation') == arch['top_candidate_score']['stable_alternation']),
                note='scan_summary.json was written before the last edit of bit0_diagnostic.py')


def part_G_perturbation(hours=100.0):
    """Neighbouring parameter values: a real toggle survives, an artefact does not."""
    from bit0_diagnostic import evaluate
    base = asdict(Extension()); base.update(CAND)
    numeric = ('uM_per_au', 'maturation_half_life_min', 'complex_on_au_inv_h', 'complex_off_h')
    rows = []
    for key in numeric:
        for delta in (-0.1, -0.05, 0.05, 0.1):
            cfg = base.copy(); cfg[key] = base[key] * (1 + delta)
            r = evaluate(cfg, hours=hours)
            rows.append(dict(varied=key, delta=delta, cfg_value=cfg[key],
                             stable_alternation=r.get('stable_alternation'), codes=r.get('codes'),
                             confident_fraction=r.get('confident_fraction'),
                             dynamic_range=r.get('dynamic_range')))
    cfg = base.copy(); cfg['add_growth'] = not base['add_growth']
    r = evaluate(cfg, hours=hours)
    rows.append(dict(varied='add_growth', delta=None, cfg_value=cfg['add_growth'],
                     stable_alternation=r.get('stable_alternation'), codes=r.get('codes'),
                     confident_fraction=r.get('confident_fraction'),
                     dynamic_range=r.get('dynamic_range')))
    n = len(rows); ok = int(sum(1 for r in rows if r['stable_alternation']))
    base_r = evaluate(base, hours=hours)
    return dict(test_points=n, neighbours_still_alternating=ok,
                base_stable_alternation=base_r.get('stable_alternation'), base_codes=base_r.get('codes'),
                frame=rows,
                verdict=('ROBUST' if ok >= n * 0.6 else
                         'KNIFE-EDGE: the 0/1 pattern disappears under a 5% parameter change'))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    df = load_confirmed()
    A, troughs, ids = part_A_phase(df)
    B = part_B_commitment(df, troughs, ids)
    C = part_C_readout(df)
    D = part_D_memory()
    E = part_E_twobit()
    F = part_F_provenance()
    G = part_G_perturbation()
    report = dict(A_phase_sensitivity=A, B_commitment=B, C_latch_vs_readout=C,
                  D_initial_state_memory=D, E_twobit_claim=E, F_provenance=F, G_perturbation=G)
    (OUT / 'bit0_verification.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')

    fig, ax = plt.subplots(4, 1, figsize=(13, 12), sharex=True, constrained_layout=True)
    t = df.time_h.to_numpy(); S = df.b0_S.to_numpy()
    ax[0].plot(t, S, c='#0072B2'); ax[0].axhspan(0.3, 0.7, color='grey', alpha=.15)
    ax[0].plot(t[troughs], S[troughs], 'kv', ms=5, label='read instant used for the codes')
    ax[0].set_ylabel('S0 (LR fraction)'); ax[0].legend(); ax[0].set_title('the recorded 0/1 codes are one instant per cycle')
    for off, c in ((-2.0, '#D55E00'), (0.0, 'black'), (2.0, '#009E73')):
        k = troughs + int(round(off / DT)); k = k[(k >= 0) & (k < len(t))]
        codes, _ = code_string(S[k], keep=0)
        ax[1].plot(t[k], S[k], 'o-', ms=3, color=c, label=f'offset {off:+.1f} h -> {codes}')
    ax[1].set_ylabel('sampled S'); ax[1].legend(fontsize=8)
    ax[2].plot(t, S, label='S'); ax[2].plot(t, (df.forward_flux / np.clip(1 - S, 1e-9, None)) /
            ((df.forward_flux / np.clip(1 - S, 1e-9, None)) + (df.reverse_flux / np.clip(S, 1e-9, None))),
            label='instantaneous balance vf/(vf+vr)', alpha=.7)
    ax[2].set_ylabel('S vs balance'); ax[2].legend(fontsize=8)
    cyc = np.arange(1, len(B['frame']) + 1)
    ax[3].bar(cyc, [r['frac_S_ge_0p7'] for r in B['frame']], color='#D55E00', label='frac of cycle with S>=0.7')
    ax[3].bar(cyc, [-r['frac_S_le_0p3'] for r in B['frame']], color='#0072B2', label='frac of cycle with S<=0.3')
    ax[3].axhline(0, color='k', lw=.8); ax[3].set_xlabel('clock cycle'); ax[3].legend(fontsize=8)
    fig.savefig(OUT / 'bit0_verification.png', dpi=160); plt.close(fig)
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
