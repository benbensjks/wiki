"""
============================================================
model_wrapper2.py

RDF system model wrapper for Sobol GSA.

According to MODE and OUTPUT defined in parameters1.py, update RDF model
parameters, run one simulation, and return one scalar feature.
============================================================
"""

import os
import sys

import numpy as np
from scipy.integrate import solve_ivp

MODEL_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "model"))
if MODEL_DIR not in sys.path:
    sys.path.append(MODEL_DIR)

import zhao_core as zc

from parameters2 import (
    BASELINE,
    DT,
    MODE,
    OUTPUT,
    SOURCE_MODE,
    T_END,
    T_END_BY_SOURCE,
    TRANSIENT_SKIP,
)
from rdf_outputs import get_output, square_sample_times, summarize_rdf


# ==========================================================
# Parameter update
# ==========================================================

def _base_params():
    P = zc.default_params()
    P.update(
        K_rep=BASELINE["K_rep"],
        n_rep=BASELINE["n_rep"],
        krep_tsl=BASELINE["krep_tsl"],
        krdf_tsl=BASELINE["krdf_tsl"],
        k_tag_int=BASELINE["k_tag_int"],
        k_int=BASELINE["k_int"],
    )
    return P


def _update_params(P, params):
    if MODE == "rdf_square":

        (
            k_int,
            pulse_width,
            kt,
            krep_tsl,
            krdf_tsl,
            k_tag_int,
            K_rep,
            n_rep,
        ) = params

        P.update(
            k_int=k_int,
            ara_on=0.5,
            ara_off=0.5 + pulse_width,
            kt=kt,
            krep_tsl=krep_tsl,
            krdf_tsl=krdf_tsl,
            k_tag_int=k_tag_int,
            K_rep=K_rep,
            n_rep=n_rep,
        )

    elif MODE == "bm3r1_promoter":

        (
            krep_tsl,
            K_rep,
            n_rep,
            krdf_tsl,
            k_tag_int,
            k_int,
            pulse_width,
        ) = params

        P.update(
            krep_tsl=krep_tsl,
            K_rep=K_rep,
            n_rep=n_rep,
            krdf_tsl=krdf_tsl,
            k_tag_int=k_tag_int,
            k_int=k_int,
            ara_on=0.5,
            ara_off=0.5 + pulse_width,
        )

    elif MODE == "rdf_waveform":

        (
            baseline_sub,
            waveform_scale,
            krep_tsl,
            krdf_tsl,
            k_tag_int,
            K_rep,
            n_rep,
        ) = params

        P.update(
            baseline_sub=baseline_sub,
            waveform_scale=waveform_scale,
            krep_tsl=krep_tsl,
            krdf_tsl=krdf_tsl,
            k_tag_int=k_tag_int,
            K_rep=K_rep,
            n_rep=n_rep,
        )

    elif MODE == "rdf_core_binding":

        for name, value in zip(
            ["Ki", "Kii", "Kir", "Ks01", "Ks02", "Kb1", "Kb2", "Dtot"],
            params,
        ):
            P[name] = value

    else:

        raise ValueError(f"Unknown MODE : {MODE}")


# ==========================================================
# Input source
# ==========================================================

def _make_source(P):
    t_end = _simulation_t_end()

    if SOURCE_MODE == "square":

        def source(t):
            return P["k_int"] * zc.square_pulse(t, P)

        return source, square_sample_times(P, t_end)

    if SOURCE_MODE == "c31_waveform":
        from analysis_lib import make_source_baseline
        from couple_oscillator import pulse_times

        source = make_source_baseline(
            P["k_dil"],
            baseline_sub=P.get("baseline_sub", BASELINE["baseline_sub"]),
            scale=P.get("waveform_scale", BASELINE["waveform_scale"]),
        )
        _, troughs = pulse_times()
        return source, troughs

    raise ValueError(f"Unknown SOURCE_MODE : {SOURCE_MODE}")


# ==========================================================
# Simulation
# ==========================================================

def _simulation_t_end():
    return T_END_BY_SOURCE.get(SOURCE_MODE, T_END)


def _initial_condition(P):
    rep_mrna = P["k_tscr"] * P["Dtot"] / P["k_rna"]
    rep = P["krep_tsl"] * rep_mrna / P["k_dil"]
    return zc.y0_PB(P, rep_mrna=rep_mrna, rep=rep)


def _simulate(P, source, t_end):
    C = zc._rate_constants(P)
    y0 = _initial_condition(P)
    t_eval = np.arange(0, t_end, DT)

    return solve_ivp(
        lambda t, y: zc.rhs(t, y, P, C, source),
        (0, t_end),
        y0,
        method="LSODA",
        rtol=1e-6,
        atol=1e-11,
        max_step=0.1,
        t_eval=t_eval,
    )


def _analytic_feature(P):
    mrna_ss = P["k_tscr"] * P["Dtot"] / P["k_rna"]
    bm3r1_production = P["krep_tsl"] * mrna_ss
    bm3r1_ss = bm3r1_production / P["k_dil"]
    threshold_margin = np.log(max(bm3r1_ss, 1e-12) / max(P["K_rep"], 1e-12))
    bm3r1_delay = max(0.0, threshold_margin / P["k_dil"])

    features = {
        "bm3r1_mrna_ss": mrna_ss,
        "bm3r1_production": bm3r1_production,
        "bm3r1_ss": bm3r1_ss,
        "bm3r1_delay_h": bm3r1_delay,
        "bm3r1_threshold_margin": threshold_margin,
    }

    return features.get(OUTPUT)


# ==========================================================
# Wrapper
# ==========================================================

def model_wrapper(params):
    """
    Parameters
    ----------
    params
        Sobol sampled RDF/counter parameters.

    Returns
    -------
    One scalar output.
    """

    P = _base_params()
    _update_params(P, params)

    feature = _analytic_feature(P)
    if feature is not None:
        return float(feature)

    t_end = _simulation_t_end()
    source, sample_times = _make_source(P)
    sol = _simulate(P, source, t_end)

    if not sol.success:
        return 0.0

    summary = summarize_rdf(
        sol,
        P,
        sample_times=sample_times,
        transient=0.2 * t_end,
        skip=TRANSIENT_SKIP,
    )

    feature = get_output(summary, OUTPUT)

    if feature is None:
        feature = 0.0
    elif np.isnan(feature):
        feature = 0.0
    elif np.isinf(feature):
        feature = 0.0

    return float(feature)

