"""
============================================================
rdf_outputs.py

Summary outputs for the Zhao RDF/counting system.

The GSA wrapper asks this module for one scalar feature, mirroring the
single-output pattern used in the oscillator GSA code.
============================================================
"""

import numpy as np

import zhao_core as zc


OUTPUT_DESCRIPTIONS = {
    "toggle_score": "Clean one-toggle-per-period score sampled between input pulses.",
    "H_mean": "Mean high-state LR fraction after transient.",
    "L_mean": "Mean low-state LR fraction after transient.",
    "alternation_fidelity": "Fraction of high/low samples that strictly alternate.",
    "switching_contrast": "H_mean - L_mean.",
    "lr_final": "Final LR DNA fraction.",
    "pb_final": "Final PB DNA fraction.",
    "lr_mean": "Mean LR fraction over the post-transient window.",
    "lr_amplitude": "Peak-to-peak LR fraction over the post-transient window.",
    "rdf_max": "Maximum total RDF concentration.",
    "rdf_mean": "Mean total RDF concentration after transient.",
    "rdf_auc": "Area under total RDF after transient.",
    "int_peak": "Maximum total integrase concentration.",
    "int_trough": "Minimum total integrase concentration after transient.",
    "rep_max": "Maximum repressor protein concentration.",
    "rep_mean": "Mean repressor protein concentration after transient.",
    "bm3r1_mrna_ss": "PB-driven BM3R1 mRNA steady state.",
    "bm3r1_production": "BM3R1 production strength, krep_tsl * mRNA_ss.",
    "bm3r1_ss": "BM3R1 protein steady state under PB expression.",
    "bm3r1_delay_h": "Estimated delay until BM3R1 decays below K_rep after promoter flip.",
    "bm3r1_threshold_margin": "log(BM3R1 steady state / K_rep).",
    "bm3r1_design_score": "Heuristic promoter-strength score favoring a 2-4 h BM3R1 delay and clean toggling.",
}


def _trapz(y, x):
    if hasattr(np, "trapezoid"):
        return np.trapezoid(y, x)
    return np.trapz(y, x)


def square_sample_times(P, t_end):
    """Mid-gap samples for the built-in square-pulse input."""

    period = P["period"]
    mid_gap = 0.5 * (P["ara_off"] + period + P["ara_on"])
    samples = []
    k = 0

    while True:
        sample_t = k * period + mid_gap
        if sample_t > t_end:
            break
        samples.append(sample_t)
        k += 1

    return np.array(samples)


def toggle_score(t, LRf, sample_times, skip=3):
    """Same metric as analysis_lib.toggle_score, kept local for portability."""

    samples_t = np.array(sample_times)
    samples_t = samples_t[samples_t <= t[-1]]

    if len(samples_t) < skip + 4:
        return 0.0, np.nan, np.nan, 0.0, []

    samples = np.interp(samples_t, t, LRf)[skip:]
    hi = samples > 0.6
    lo = samples < 0.4
    mid = ~(hi | lo)

    alt = 0
    n = 0
    for k in range(len(samples) - 1):
        if mid[k] or mid[k + 1]:
            continue
        n += 1
        if hi[k] != hi[k + 1]:
            alt += 1

    fidelity = alt / n if n else 0.0
    H = samples[hi].mean() if hi.any() else np.nan
    L = samples[lo].mean() if lo.any() else np.nan
    core = min(np.nan_to_num(H, nan=0.0), 1 - np.nan_to_num(L, nan=1.0))

    return core * fidelity, H, L, fidelity, list(zip(samples_t[skip:], samples))


def summarize_rdf(sol, P, sample_times=None, transient=20.0, skip=3):
    """Return all scalar RDF/counting-system outputs in one dictionary."""

    t = sol.t
    Y = sol.y.T

    LRf = zc.LR_total(Y) / P["Dtot"]
    PBf = zc.PB_total(Y) / P["Dtot"]
    INT = zc.int_total(Y)
    RDF = zc.rdf_total(Y)
    REP = Y[:, 37]

    if sample_times is None:
        sample_times = square_sample_times(P, t[-1])

    score, H, L, fidelity, samples = toggle_score(
        t,
        LRf,
        sample_times,
        skip=skip,
    )

    mask = t > transient
    if not np.any(mask):
        mask = np.ones_like(t, dtype=bool)

    mrna_ss = P["k_tscr"] * P["Dtot"] / P["k_rna"]
    bm3r1_production = P["krep_tsl"] * mrna_ss
    bm3r1_ss = bm3r1_production / P["k_dil"]
    threshold_margin = np.log(max(bm3r1_ss, 1e-12) / max(P["K_rep"], 1e-12))
    bm3r1_delay = max(0.0, threshold_margin / P["k_dil"])

    # A practical promoter design metric: prefer a delay centered near 3 h,
    # keep enough threshold margin, and include the observed toggle score.
    delay_score = np.exp(-0.5 * ((bm3r1_delay - 3.0) / 1.0) ** 2)
    margin_score = 1.0 / (1.0 + np.exp(-(threshold_margin - 3.0)))
    design_score = delay_score * margin_score * (0.5 + 0.5 * score)

    return {
        "toggle_score": score,
        "H_mean": H,
        "L_mean": L,
        "alternation_fidelity": fidelity,
        "switching_contrast": H - L if H == H and L == L else np.nan,
        "lr_final": LRf[-1],
        "pb_final": PBf[-1],
        "lr_mean": np.mean(LRf[mask]),
        "lr_amplitude": np.ptp(LRf[mask]),
        "rdf_max": np.max(RDF),
        "rdf_mean": np.mean(RDF[mask]),
        "rdf_auc": _trapz(RDF[mask], t[mask]),
        "int_peak": np.max(INT),
        "int_trough": np.min(INT[mask]),
        "rep_max": np.max(REP),
        "rep_mean": np.mean(REP[mask]),
        "bm3r1_mrna_ss": mrna_ss,
        "bm3r1_production": bm3r1_production,
        "bm3r1_ss": bm3r1_ss,
        "bm3r1_delay_h": bm3r1_delay,
        "bm3r1_threshold_margin": threshold_margin,
        "bm3r1_design_score": design_score,
        "samples": samples,
    }


def get_output(summary, output):
    if output not in OUTPUT_DESCRIPTIONS:
        raise ValueError(f"Unknown OUTPUT : {output}")

    return summary[output]
