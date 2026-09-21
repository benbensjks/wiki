"""Morris global sensitivity screening for the mechanistic repressilator-C31 ODE.

This program keeps the frozen v52a1 state equations and the v53d dynamic
translation-resource term.  It changes parameter values only; it does not add
new biochemical states and it does not fit parameters to a desired phenotype.

Run a short workflow check first:
    python Mechanistic_ODE_Global_Sensitivity_Analysis.py --mode smoke

Run the preregistered screening design:
    python Mechanistic_ODE_Global_Sensitivity_Analysis.py --mode full

Resume an interrupted run by reusing its output directory:
    python Mechanistic_ODE_Global_Sensitivity_Analysis.py --mode full \
        --output-dir "D:\\...\\outputs\\Morris_GSA_full_YYYYMMDD_HHMMSS"

The model uses copies/cell and minutes internally.  C31 production flux is
reported in uM/h using an explicitly recorded cell-volume assumption.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass, replace
from datetime import datetime
from pathlib import Path
from typing import Callable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp, trapezoid
from scipy.optimize import brentq
from scipy.signal import find_peaks


VERSION = "Mechanistic_ODE_Morris_GSA_20260819_v2_interface_leak_metrics"
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
DEFAULT_OUTPUT_ROOT = PROJECT_ROOT / "outputs"
DEFAULT_FIGURE_ROOT = PROJECT_ROOT / "figures"


@dataclass(frozen=True)
class SimulationSettings:
    duration_min: float
    sample_step_min: float
    transient_min: float
    morris_trajectories: int
    morris_levels: int
    bootstrap_replicates: int
    random_seed: int = 20260819
    minimum_relative_amplitude: float = 0.01
    minimum_peaks: int = 4
    cell_volume_fL: float = 1.0
    active_C31_fraction: float = 1.0
    molecules_per_uM_per_fL: float = 602.214076


MODE_SETTINGS = {
    # Smoke mode checks the complete workflow.  Its rankings are not evidence.
    "smoke": SimulationSettings(
        duration_min=3600.0,
        sample_step_min=2.0,
        transient_min=900.0,
        morris_trajectories=2,
        morris_levels=4,
        bootstrap_replicates=300,
    ),
    # 12 trajectories x (15 parameters + 1) = 192 model evaluations.
    "full": SimulationSettings(
        duration_min=6000.0,
        sample_step_min=1.0,
        transient_min=1500.0,
        morris_trajectories=12,
        morris_levels=4,
        bootstrap_replicates=2000,
    ),
}

# This is a reporting threshold for the A-to-B interface, not a fitted
# biochemical threshold in the Zhao B model.  It makes "residual input" an
# explicit output beside a dimensionless peak-normalized trough metric.
B_INTERFACE_RESIDUAL_FLUX_THRESHOLD_uM_h = 0.5

# A reference trajectory has an almost-zero C31 trough.  Dividing trough
# metrics by that value creates a misleadingly huge Morris mu*.  Fractions and
# duties are already bounded quantities; their scientific scale is one.
NORMALIZATION_SCALE_OVERRIDES = {
    "resolved_oscillation": 1.0,
    "C31_flux_trough_uM_h": 1.0,
    "C31_flux_trough_fraction_of_peak": 1.0,
    "C31_flux_residual_duty_above_0p5_uM_h": 1.0,
    "C31_flux_trough_fraction_of_peak_when_resolved": 1.0,
    "C31_flux_residual_duty_above_0p5_uM_h_when_resolved": 1.0,
}


@dataclass(frozen=True)
class Parameters:
    # Growth, repression chemistry and plasmid/site accounting.
    Td_min: float = 50.0
    K_copies: float = 13.0
    hill_n: float = 3.0
    oscillator_plasmid_copies: float = 10.0
    functional_tetO_per_plasmid: float = 4.0
    functional_accessible_fraction: float = 1.0
    sponge_plasmid_copies: float = 40.0
    sponge_tetO_per_plasmid: float = 2.0
    sponge_accessible_fraction: float = 1.0
    oscillator_promoter_leak: float = 0.0

    # Explicit oscillator mRNA/protein layer.  Reference rates are frozen once.
    oscillator_tx_per_plasmid_per_min: float = 0.0
    oscillator_translation_per_mrna_per_min: float = 0.5
    oscillator_mrna_intrinsic_loss_per_min: float = 0.0
    oscillator_protein_intrinsic_loss_per_min: float = 0.0

    # C31 replaces the reporter and uses the same PLtetO1 transfer function.
    c31_promoter_leak: float = 0.005
    c31_tx_per_plasmid_per_min: float = 5.0
    c31_translation_per_mrna_per_min: float = 0.5
    c31_mrna_intrinsic_loss_per_min: float = 0.0
    c31_protein_intrinsic_loss_per_min: float = 0.0

    # v53d finite shared-translation-resource boundary parameter.
    peak_load_fraction: float = 0.0

    @property
    def mu(self) -> float:
        return math.log(2.0) / self.Td_min

    @property
    def effective_functional_sites(self) -> float:
        return (
            self.functional_accessible_fraction
            * self.functional_tetO_per_plasmid
            * self.oscillator_plasmid_copies
        )

    @property
    def effective_sponge_sites(self) -> float:
        return (
            self.sponge_accessible_fraction
            * self.sponge_tetO_per_plasmid
            * self.sponge_plasmid_copies
        )


@dataclass(frozen=True)
class ParameterSpec:
    name: str
    label: str
    lower: float
    upper: float
    scale: str
    meaning: str
    evidence_class: str

    def from_unit(self, x: float) -> float:
        value = float(np.clip(x, 0.0, 1.0))
        if self.scale == "linear":
            return self.lower + value * (self.upper - self.lower)
        if self.scale == "log":
            return math.exp(math.log(self.lower) + value * math.log(self.upper / self.lower))
        raise ValueError(f"Unknown scale {self.scale!r} for {self.name}")


# These are uncertainty/design ranges, not fitted confidence intervals.
# Common K, n and oscillator leak retain the current symmetric parameterization;
# arm-specific transfer curves should replace them only after measurements exist.
PARAMETER_SPECS = (
    ParameterSpec("Td_min", "doubling time", 40.0, 60.0, "linear",
                  "MC4100 doubling time in minutes", "measured-range target"),
    ParameterSpec("oscillator_plasmid_copies", "oscillator plasmid copies", 7.0, 13.0, "linear",
                  "mean copies of the oscillator plasmid per cell", "measurement target"),
    ParameterSpec("K_copies", "common repression K", 8.0, 18.0, "linear",
                  "effective half-repression threshold in copies/cell", "effective model range"),
    ParameterSpec("hill_n", "common Hill n", 2.0, 4.0, "linear",
                  "effective repression-curve steepness", "effective model range"),
    ParameterSpec("oscillator_promoter_leak", "oscillator promoter leak", 0.0, 0.02, "linear",
                  "fractional lower asymptote shared by the three oscillator arms", "measurement target"),
    ParameterSpec("oscillator_tx_scale", "oscillator transcription scale", 0.5, 1.5, "linear",
                  "multiplicative scale on the frozen oscillator alpha", "engineering uncertainty"),
    ParameterSpec("oscillator_translation_scale", "oscillator translation scale", 0.5, 1.5, "linear",
                  "multiplicative scale on the frozen oscillator beta", "engineering uncertainty"),
    ParameterSpec("oscillator_mrna_total_half_life_min", "oscillator mRNA half-life", 1.0, 5.0, "log",
                  "total mRNA half-life including dilution", "measurement target"),
    ParameterSpec("sponge_plasmid_copies", "sponge plasmid copies Nt", 0.0, 80.0, "linear",
                  "mean sponge plasmid copies per cell, not TetO-site count", "design range"),
    ParameterSpec("sponge_accessible_fraction", "sponge accessibility eta", 0.5, 1.0, "linear",
                  "effective accessible fraction of sponge TetO sites", "uncalibrated uncertainty"),
    ParameterSpec("c31_promoter_leak", "PLtetO1-C31 leak", 0.001, 0.02, "log",
                  "fractional PLtetO1 lower asymptote for C31", "measurement target"),
    ParameterSpec("c31_tx_scale", "C31 transcription scale", 0.5, 1.5, "linear",
                  "multiplicative scale on C31 transcription", "engineering uncertainty"),
    ParameterSpec("c31_translation_scale", "C31 RBS scale", 0.3, 1.2, "linear",
                  "multiplicative scale on C31 translation", "engineering design range"),
    ParameterSpec("c31_mrna_total_half_life_min", "C31 mRNA half-life", 0.5, 5.0, "log",
                  "total C31 mRNA half-life including dilution", "measurement target"),
    ParameterSpec("peak_load_fraction", "peak shared-resource load", 0.0, 0.30, "linear",
                  "nominal translation-resource loss at reference C31 demand", "uncalibrated boundary"),
)


def make_reference_parameters() -> Parameters:
    """Create the frozen v52a1 reference using the QSS mapping only once."""
    Td0 = 50.0
    mu0 = math.log(2.0) / Td0
    total_mrna_loss0 = math.log(2.0) / 2.0
    intrinsic_mrna_loss0 = total_mrna_loss0 - mu0
    beta0 = 0.5
    No0 = 10.0
    lambda0 = 1000.0
    alpha0 = mu0 * lambda0 * total_mrna_loss0 / (beta0 * No0)
    return Parameters(
        Td_min=Td0,
        oscillator_tx_per_plasmid_per_min=alpha0,
        oscillator_translation_per_mrna_per_min=beta0,
        oscillator_mrna_intrinsic_loss_per_min=intrinsic_mrna_loss0,
        c31_mrna_intrinsic_loss_per_min=intrinsic_mrna_loss0,
    )


def repression(repressor: float, p: Parameters, leak: float) -> float:
    r = max(float(repressor), 0.0)
    regulated = p.K_copies**p.hill_n / (p.K_copies**p.hill_n + r**p.hill_n)
    return leak + (1.0 - leak) * regulated


def occupancy(free: float, p: Parameters) -> float:
    f = max(float(free), 0.0)
    if f == 0.0:
        return 0.0
    return f**p.hill_n / (p.K_copies**p.hill_n + f**p.hill_n)


def free_tetR(total: float, p: Parameters) -> float:
    total = max(float(total), 0.0)
    if total == 0.0:
        return 0.0
    sites = p.effective_functional_sites + p.effective_sponge_sites

    def balance(free: float) -> float:
        return free + sites * occupancy(free, p) - total

    free = brentq(balance, 0.0, total, xtol=1e-12, rtol=1e-12, maxiter=100)
    residual = abs(balance(free)) / max(total, 1.0)
    if residual >= 1e-6 or not 0.0 <= free <= total:
        raise RuntimeError(f"TetR conservation failed: residual={residual:g}")
    return free


def promoter_activities(tet_total: float, ci: float, laci: float, p: Parameters):
    tet_free = free_tetR(tet_total, p)
    return (
        repression(laci, p, p.oscillator_promoter_leak),
        repression(tet_free, p, p.oscillator_promoter_leak),
        repression(ci, p, p.oscillator_promoter_leak),
        repression(tet_free, p, p.c31_promoter_leak),
        tet_free,
    )


def resource_availability(m31: float, p: Parameters, reference_demand_peak: float) -> float:
    load = float(p.peak_load_fraction)
    if not 0.0 <= load < 1.0:
        raise ValueError("peak_load_fraction must be in [0, 1)")
    if load == 0.0:
        return 1.0
    q31 = p.c31_translation_per_mrna_per_min * max(float(m31), 0.0) / reference_demand_peak
    burden_ratio = load / (1.0 - load)
    return 1.0 / (1.0 + burden_ratio * q31)


def rhs(_: float, y: np.ndarray, p: Parameters, reference_demand_peak: float):
    m_tet, tet, m_ci, ci, m_laci, laci, m31, c31 = y
    h_tet, h_ci, h_laci, h31, _ = promoter_activities(tet, ci, laci, p)
    rho = resource_availability(m31, p, reference_demand_peak)
    mrna_loss = p.oscillator_mrna_intrinsic_loss_per_min + p.mu
    protein_loss = p.oscillator_protein_intrinsic_loss_per_min + p.mu
    c31_mrna_loss = p.c31_mrna_intrinsic_loss_per_min + p.mu
    c31_protein_loss = p.c31_protein_intrinsic_loss_per_min + p.mu
    No = p.oscillator_plasmid_copies
    alpha = p.oscillator_tx_per_plasmid_per_min
    beta = p.oscillator_translation_per_mrna_per_min
    return (
        No * alpha * h_tet - mrna_loss * m_tet,
        rho * beta * m_tet - protein_loss * tet,
        No * alpha * h_ci - mrna_loss * m_ci,
        rho * beta * m_ci - protein_loss * ci,
        No * alpha * h_laci - mrna_loss * m_laci,
        rho * beta * m_laci - protein_loss * laci,
        No * p.c31_tx_per_plasmid_per_min * h31 - c31_mrna_loss * m31,
        rho * p.c31_translation_per_mrna_per_min * m31 - c31_protein_loss * c31,
    )


def initial_state(p: Parameters):
    tet, ci, laci = 100.0, 20.0, 300.0
    h_tet, h_ci, h_laci, h31, _ = promoter_activities(tet, ci, laci, p)
    mrna_loss = p.oscillator_mrna_intrinsic_loss_per_min + p.mu
    c31_mrna_loss = p.c31_mrna_intrinsic_loss_per_min + p.mu
    No = p.oscillator_plasmid_copies
    alpha = p.oscillator_tx_per_plasmid_per_min
    return (
        No * alpha * h_tet / mrna_loss, tet,
        No * alpha * h_ci / mrna_loss, ci,
        No * alpha * h_laci / mrna_loss, laci,
        No * p.c31_tx_per_plasmid_per_min * h31 / c31_mrna_loss, 0.0,
    )


def apply_parameter_vector(unit_values: np.ndarray, base: Parameters) -> tuple[Parameters, dict]:
    physical = {spec.name: spec.from_unit(value) for spec, value in zip(PARAMETER_SPECS, unit_values)}
    p = replace(
        base,
        Td_min=physical["Td_min"],
        oscillator_plasmid_copies=physical["oscillator_plasmid_copies"],
        K_copies=physical["K_copies"],
        hill_n=physical["hill_n"],
        oscillator_promoter_leak=physical["oscillator_promoter_leak"],
        oscillator_tx_per_plasmid_per_min=(
            base.oscillator_tx_per_plasmid_per_min * physical["oscillator_tx_scale"]),
        oscillator_translation_per_mrna_per_min=(
            base.oscillator_translation_per_mrna_per_min
            * physical["oscillator_translation_scale"]),
        sponge_plasmid_copies=physical["sponge_plasmid_copies"],
        sponge_accessible_fraction=physical["sponge_accessible_fraction"],
        c31_promoter_leak=physical["c31_promoter_leak"],
        c31_tx_per_plasmid_per_min=(base.c31_tx_per_plasmid_per_min * physical["c31_tx_scale"]),
        c31_translation_per_mrna_per_min=(
            base.c31_translation_per_mrna_per_min * physical["c31_translation_scale"]),
        peak_load_fraction=physical["peak_load_fraction"],
    )
    # Total half-life is the experimentally interpretable input.  Intrinsic
    # degradation is derived after Td has been applied, without remapping alpha.
    osc_total_loss = math.log(2.0) / physical["oscillator_mrna_total_half_life_min"]
    c31_total_loss = math.log(2.0) / physical["c31_mrna_total_half_life_min"]
    p = replace(
        p,
        oscillator_mrna_intrinsic_loss_per_min=max(0.0, osc_total_loss - p.mu),
        c31_mrna_intrinsic_loss_per_min=max(0.0, c31_total_loss - p.mu),
    )
    return p, physical


def build_morris_design(settings: SimulationSettings):
    """Return a reproducible elementary-effects design on the unit hypercube."""
    rng = np.random.default_rng(settings.random_seed)
    dimension = len(PARAMETER_SPECS)
    levels = np.linspace(0.0, 1.0, settings.morris_levels)
    delta = settings.morris_levels / (2.0 * (settings.morris_levels - 1.0))
    rows = []
    effect_links = []
    run_id = 0
    for trajectory in range(settings.morris_trajectories):
        directions = rng.choice((-1.0, 1.0), size=dimension)
        x = np.empty(dimension, dtype=float)
        for index, direction in enumerate(directions):
            valid = levels[(levels + direction * delta >= -1e-12)
                           & (levels + direction * delta <= 1.0 + 1e-12)]
            x[index] = float(rng.choice(valid))
        order = rng.permutation(dimension)
        rows.append({
            "run_id": run_id,
            "trajectory_id": trajectory,
            "step_id": 0,
            "changed_parameter": "START",
            **{f"x__{spec.name}": x[i] for i, spec in enumerate(PARAMETER_SPECS)},
        })
        before_run = run_id
        run_id += 1
        for step, parameter_index in enumerate(order, start=1):
            before_x = float(x[parameter_index])
            x = x.copy()
            x[parameter_index] += directions[parameter_index] * delta
            after_x = float(x[parameter_index])
            rows.append({
                "run_id": run_id,
                "trajectory_id": trajectory,
                "step_id": step,
                "changed_parameter": PARAMETER_SPECS[parameter_index].name,
                **{f"x__{spec.name}": x[i] for i, spec in enumerate(PARAMETER_SPECS)},
            })
            effect_links.append({
                "trajectory_id": trajectory,
                "parameter": PARAMETER_SPECS[parameter_index].name,
                "before_run_id": before_run,
                "after_run_id": run_id,
                "delta_x": after_x - before_x,
            })
            before_run = run_id
            run_id += 1
    return pd.DataFrame(rows), pd.DataFrame(effect_links), delta


def oscillation_metrics(t: np.ndarray, signal: np.ndarray, p: Parameters,
                        settings: SimulationSettings):
    span = float(np.ptp(signal))
    scale = max(float(np.max(np.abs(signal))), 1.0)
    relative_amplitude = span / scale
    empty = np.array([], dtype=int)
    if relative_amplitude < settings.minimum_relative_amplitude:
        return "steady_state", float("nan"), float("nan"), float("nan"), empty
    minimum_distance = max(1, int(round(2.0 * p.Td_min / settings.sample_step_min)))
    peaks, _ = find_peaks(signal, prominence=max(0.05 * span, 1e-9), distance=minimum_distance)
    if len(peaks) < settings.minimum_peaks:
        return "unresolved_oscillation", float("nan"), float("nan"), relative_amplitude, peaks
    periods = np.diff(t[peaks])
    period = float(np.mean(periods))
    period_cv = float(np.std(periods, ddof=1) / period) if len(periods) > 1 else 0.0
    return "resolved_oscillation", period, period_cv, relative_amplitude, peaks


def median_peak_fwhm(t: np.ndarray, signal: np.ndarray, peaks: np.ndarray) -> float:
    widths = []
    baseline = float(np.min(signal))
    for peak in peaks[-6:]:
        level = baseline + 0.5 * (float(signal[peak]) - baseline)
        left = int(peak)
        right = int(peak)
        while left > 0 and signal[left] >= level:
            left -= 1
        while right < len(signal) - 1 and signal[right] >= level:
            right += 1
        if left > 0 and right < len(signal) - 1:
            widths.append(float(t[right] - t[left]))
    return float(np.median(widths)) if widths else float("nan")


def last_cycle_dose(t: np.ndarray, signal: np.ndarray, peaks: np.ndarray) -> float:
    if len(peaks) < 2:
        return float("nan")
    segment = slice(peaks[-2], peaks[-1] + 1)
    return float(trapezoid(signal[segment], x=t[segment]))


def phase_lag_fraction(t: np.ndarray, tet_peaks: np.ndarray, ci_signal: np.ndarray,
                       period_min: float, p: Parameters, settings: SimulationSettings) -> float:
    if len(tet_peaks) < 2 or not np.isfinite(period_min):
        return float("nan")
    ci_span = float(np.ptp(ci_signal))
    ci_peaks, _ = find_peaks(
        ci_signal,
        prominence=max(0.05 * ci_span, 1e-9),
        distance=max(1, int(round(2.0 * p.Td_min / settings.sample_step_min))),
    )
    lags = []
    for tet_peak in tet_peaks[-5:]:
        following = ci_peaks[ci_peaks > tet_peak]
        if len(following):
            lag = float(t[following[0]] - t[tet_peak])
            if 0.0 <= lag <= 1.5 * period_min:
                lags.append((lag / period_min) % 1.0)
    return float(np.mean(lags)) if lags else float("nan")


def reference_demand_peak(base: Parameters, settings: SimulationSettings) -> float:
    """Compute Jref once from the frozen no-load reference trajectory."""
    p = replace(base, peak_load_fraction=0.0)
    time = np.arange(0.0, settings.duration_min + 0.5 * settings.sample_step_min,
                     settings.sample_step_min)
    sol = solve_ivp(
        lambda t, y: rhs(t, y, p, 1.0),
        (time[0], time[-1]), initial_state(p), t_eval=time,
        method="LSODA", rtol=1e-7, atol=1e-9,
    )
    if not sol.success:
        raise RuntimeError(f"Reference simulation failed: {sol.message}")
    keep = time >= settings.transient_min
    peak = float(np.max(p.c31_translation_per_mrna_per_min * sol.y[6, keep]))
    if not np.isfinite(peak) or peak <= 0.0:
        raise RuntimeError("Reference C31 demand peak is not positive")
    return peak


def simulate_run(run_id: int, unit_values: np.ndarray, base: Parameters,
                 settings: SimulationSettings, reference_peak: float) -> dict:
    p, physical = apply_parameter_vector(unit_values, base)
    time = np.arange(0.0, settings.duration_min + 0.5 * settings.sample_step_min,
                     settings.sample_step_min)
    keep = time >= settings.transient_min
    t_ss = time[keep]
    row = {"run_id": int(run_id), **physical}
    try:
        sol = solve_ivp(
            lambda t, y: rhs(t, y, p, reference_peak),
            (time[0], time[-1]), initial_state(p), t_eval=time,
            method="LSODA", rtol=1e-7, atol=1e-9,
        )
        if not sol.success:
            raise RuntimeError(sol.message)
        tet_free = np.array([free_tetR(total, p) for total in sol.y[1]])
        activity = np.array([repression(value, p, p.c31_promoter_leak) for value in tet_free])
        resource = np.array([resource_availability(value, p, reference_peak) for value in sol.y[6]])
        demand = p.c31_translation_per_mrna_per_min * sol.y[6]
        actual_flux_uM_h = (
            60.0 * resource * demand * settings.active_C31_fraction
            / (settings.molecules_per_uM_per_fL * settings.cell_volume_fL)
        )
        status, period_min, period_cv, relative_amplitude, tet_peaks = oscillation_metrics(
            t_ss, tet_free[keep], p, settings
        )
        flux_span = float(np.ptp(actual_flux_uM_h[keep]))
        flux_peaks, _ = find_peaks(
            actual_flux_uM_h[keep],
            prominence=max(0.05 * flux_span, 1e-9),
            distance=max(1, int(round(2.0 * p.Td_min / settings.sample_step_min))),
        )
        row.update({
            "solver_success": True,
            "error_message": "",
            "status": status,
            "resolved_oscillation": float(status == "resolved_oscillation"),
            "period_min": period_min,
            "period_generations": period_min / p.Td_min if np.isfinite(period_min) else np.nan,
            "deterministic_period_cv": period_cv,
            "TetR_relative_amplitude": relative_amplitude,
            "TetR_free_peak": float(np.max(tet_free[keep])),
            "TetR_free_trough": float(np.min(tet_free[keep])),
            "CI_peak": float(np.max(sol.y[3, keep])),
            "LacI_peak": float(np.max(sol.y[5, keep])),
            "CI_phase_lag_fraction_cycle": phase_lag_fraction(
                t_ss, tet_peaks, sol.y[3, keep], period_min, p, settings),
            "PLtetO1_duty_activity_gt_0p5": float(np.mean(activity[keep] > 0.5)),
            "C31_protein_peak_copies": float(np.max(sol.y[7, keep])),
            "C31_protein_trough_copies": float(np.min(sol.y[7, keep])),
            "C31_flux_peak_uM_h": float(np.max(actual_flux_uM_h[keep])),
            "C31_flux_trough_uM_h": float(np.min(actual_flux_uM_h[keep])),
            "C31_flux_trough_fraction_of_peak": float(
                np.min(actual_flux_uM_h[keep])
                / max(np.max(actual_flux_uM_h[keep]), 1e-12)
            ),
            "C31_flux_residual_duty_above_0p5_uM_h": float(
                np.mean(actual_flux_uM_h[keep] > B_INTERFACE_RESIDUAL_FLUX_THRESHOLD_uM_h)
            ),
            "C31_flux_FWHM_h": median_peak_fwhm(t_ss, actual_flux_uM_h[keep], flux_peaks) / 60.0,
            "C31_flux_dose_uM_per_cycle": last_cycle_dose(
                t_ss, actual_flux_uM_h[keep], tet_peaks) / 60.0,
            "resource_min": float(np.min(resource[keep])),
            "resource_mean": float(np.mean(resource[keep])),
            "effective_functional_sites": p.effective_functional_sites,
            "effective_sponge_sites": p.effective_sponge_sites,
        })
        # These conditional outputs isolate C31 interface leakage from total
        # system failure.  The corresponding all-state metrics above remain
        # useful for a complete-system risk analysis.
        if status == "resolved_oscillation":
            row["C31_flux_trough_fraction_of_peak_when_resolved"] = row[
                "C31_flux_trough_fraction_of_peak"
            ]
            row["C31_flux_residual_duty_above_0p5_uM_h_when_resolved"] = row[
                "C31_flux_residual_duty_above_0p5_uM_h"
            ]
        else:
            row["C31_flux_trough_fraction_of_peak_when_resolved"] = np.nan
            row["C31_flux_residual_duty_above_0p5_uM_h_when_resolved"] = np.nan
    except Exception as exc:  # A failed parameter point is an outcome, not a lost run.
        row.update({
            "solver_success": False,
            "error_message": f"{type(exc).__name__}: {exc}",
            "status": "solver_failure",
            "resolved_oscillation": 0.0,
        })
        for metric in METRICS:
            row.setdefault(metric, np.nan)
    return row


METRICS = (
    "resolved_oscillation",
    "period_min",
    "period_generations",
    "TetR_relative_amplitude",
    "TetR_free_peak",
    "TetR_free_trough",
    "CI_phase_lag_fraction_cycle",
    "PLtetO1_duty_activity_gt_0p5",
    "C31_flux_peak_uM_h",
    "C31_flux_trough_uM_h",
    "C31_flux_trough_fraction_of_peak",
    "C31_flux_residual_duty_above_0p5_uM_h",
    "C31_flux_trough_fraction_of_peak_when_resolved",
    "C31_flux_residual_duty_above_0p5_uM_h_when_resolved",
    "C31_flux_FWHM_h",
    "C31_flux_dose_uM_per_cycle",
    "resource_min",
)


def compute_elementary_effects(metrics: pd.DataFrame, links: pd.DataFrame,
                               reference_scales: dict[str, float]) -> pd.DataFrame:
    indexed = metrics.set_index("run_id")
    rows = []
    for link in links.to_dict("records"):
        before = indexed.loc[int(link["before_run_id"])]
        after = indexed.loc[int(link["after_run_id"])]
        for metric in METRICS:
            y0 = float(before.get(metric, np.nan))
            y1 = float(after.get(metric, np.nan))
            valid = np.isfinite(y0) and np.isfinite(y1) and link["delta_x"] != 0.0
            effect = (y1 - y0) / float(link["delta_x"]) if valid else np.nan
            scale = reference_scales.get(metric, 1.0)
            normalized = effect / scale if valid and np.isfinite(scale) and scale > 0.0 else np.nan
            rows.append({
                **link,
                "metric": metric,
                "before_value": y0,
                "after_value": y1,
                "elementary_effect": effect,
                "normalized_elementary_effect": normalized,
                "valid_pair": bool(valid),
            })
    return pd.DataFrame(rows)


def bootstrap_mu_star_ci(values: np.ndarray, replicates: int, rng: np.random.Generator):
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return np.nan, np.nan
    if len(values) == 1:
        value = float(abs(values[0]))
        return value, value
    estimates = np.empty(replicates, dtype=float)
    for index in range(replicates):
        sample = rng.choice(values, size=len(values), replace=True)
        estimates[index] = float(np.mean(np.abs(sample)))
    return tuple(np.quantile(estimates, (0.025, 0.975)))


def summarize_effects(effects: pd.DataFrame, settings: SimulationSettings) -> pd.DataFrame:
    rng = np.random.default_rng(settings.random_seed + 1)
    rows = []
    for (metric, parameter), group in effects.groupby(["metric", "parameter"], sort=False):
        raw = group["elementary_effect"].to_numpy(dtype=float)
        normalized = group["normalized_elementary_effect"].to_numpy(dtype=float)
        valid_raw = raw[np.isfinite(raw)]
        valid_normalized = normalized[np.isfinite(normalized)]
        lower, upper = bootstrap_mu_star_ci(
            valid_normalized, settings.bootstrap_replicates, rng
        )
        rows.append({
            "metric": metric,
            "parameter": parameter,
            "n_effects_total": len(group),
            "n_effects_valid": len(valid_raw),
            "valid_fraction": len(valid_raw) / len(group),
            "mu_raw": float(np.mean(valid_raw)) if len(valid_raw) else np.nan,
            "mu_star_raw": float(np.mean(np.abs(valid_raw))) if len(valid_raw) else np.nan,
            "sigma_raw": float(np.std(valid_raw, ddof=1)) if len(valid_raw) > 1 else 0.0,
            "mu_normalized": (float(np.mean(valid_normalized))
                              if len(valid_normalized) else np.nan),
            "mu_star_normalized": (float(np.mean(np.abs(valid_normalized)))
                                   if len(valid_normalized) else np.nan),
            "sigma_normalized": (float(np.std(valid_normalized, ddof=1))
                                 if len(valid_normalized) > 1 else 0.0),
            "mu_star_normalized_ci95_low": lower,
            "mu_star_normalized_ci95_high": upper,
        })
    summary = pd.DataFrame(rows)
    summary["rank_within_metric"] = summary.groupby("metric")["mu_star_normalized"].rank(
        ascending=False, method="min"
    )
    return summary.sort_values(["metric", "rank_within_metric", "parameter"])


def make_plots(summary: pd.DataFrame, figure_dir: Path):
    figure_dir.mkdir(parents=True, exist_ok=True)
    selected_metrics = (
        "period_generations",
        "TetR_relative_amplitude",
        "C31_flux_peak_uM_h",
        "C31_flux_trough_fraction_of_peak_when_resolved",
        "C31_flux_residual_duty_above_0p5_uM_h_when_resolved",
        "C31_flux_FWHM_h",
        "C31_flux_dose_uM_per_cycle",
        "resolved_oscillation",
        "resource_min",
    )
    fig, axes = plt.subplots(3, 3, figsize=(21, 15), constrained_layout=True)
    for ax, metric in zip(axes.flat, selected_metrics):
        data = summary[summary["metric"] == metric].dropna(subset=["mu_star_normalized"])
        data = data.sort_values("mu_star_normalized", ascending=True).tail(10)
        ax.barh(data["parameter"], data["mu_star_normalized"], color="tab:blue", alpha=0.85)
        ax.set_xlabel("normalized Morris mu* (larger = more influential)")
        ax.set_title(metric)
        ax.grid(axis="x", alpha=0.25)
    fig.suptitle("Mechanistic repressilator-C31 ODE: Morris global sensitivity screening")
    fig.savefig(figure_dir / "Morris_GSA_top_parameters.png", dpi=180)
    plt.close(fig)

    heat = summary.pivot(index="parameter", columns="metric", values="mu_star_normalized")
    heat = heat.reindex(columns=[metric for metric in selected_metrics if metric in heat.columns])
    order = heat.fillna(0.0).max(axis=1).sort_values().index
    heat = heat.loc[order]
    fig, ax = plt.subplots(figsize=(12, 8), constrained_layout=True)
    image = ax.imshow(heat.to_numpy(dtype=float), aspect="auto", cmap="viridis")
    ax.set_yticks(np.arange(len(heat.index)), labels=heat.index)
    ax.set_xticks(np.arange(len(heat.columns)), labels=heat.columns, rotation=35, ha="right")
    ax.set_title("Normalized Morris mu* across clock and C31 outputs")
    fig.colorbar(image, ax=ax, label="normalized mu*")
    fig.savefig(figure_dir / "Morris_GSA_metric_heatmap.png", dpi=180)
    plt.close(fig)


def write_checkpoint(rows: list[dict], path: Path):
    pd.DataFrame(rows).sort_values("run_id").to_csv(path, index=False, encoding="utf-8-sig")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=tuple(MODE_SETTINGS), default="smoke",
                        help="smoke checks the workflow; full is the evidential Morris screen")
    parser.add_argument("--trajectories", type=int, default=None,
                        help="override the number of Morris trajectories")
    parser.add_argument("--seed", type=int, default=None, help="override the random seed")
    parser.add_argument("--output-dir", type=Path, default=None,
                        help="reuse an existing directory to resume a run")
    parser.add_argument("--no-plots", action="store_true")
    return parser.parse_args()


def main():
    args = parse_args()
    settings = MODE_SETTINGS[args.mode]
    if args.trajectories is not None:
        if args.trajectories < 2:
            raise ValueError("At least two Morris trajectories are required")
        settings = replace(settings, morris_trajectories=args.trajectories)
    if args.seed is not None:
        settings = replace(settings, random_seed=args.seed)

    run_stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = (args.output_dir.resolve() if args.output_dir is not None else
                  (DEFAULT_OUTPUT_ROOT / f"Morris_GSA_{args.mode}_{run_stamp}").resolve())
    figure_dir = DEFAULT_FIGURE_ROOT / output_dir.name
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    design, links, delta = build_morris_design(settings)
    design_path = output_dir / "morris_design.csv"
    links_path = output_dir / "morris_effect_links.csv"
    if design_path.exists():
        old_design = pd.read_csv(design_path)
        same_columns = list(old_design.columns) == list(design.columns)
        same_shape = old_design.shape == design.shape
        same_text = (
            same_columns and same_shape
            and old_design[["changed_parameter"]].astype(str).equals(
                design[["changed_parameter"]].astype(str)
            )
        )
        numeric_columns = [column for column in design.columns if column != "changed_parameter"]
        same_numbers = (
            same_columns and same_shape
            and np.allclose(
                old_design[numeric_columns].to_numpy(dtype=float),
                design[numeric_columns].to_numpy(dtype=float),
                rtol=0.0,
                atol=1e-12,
            )
        )
        if not (same_text and same_numbers):
            raise RuntimeError("Existing output directory contains a different Morris design")
        if not links_path.exists():
            links.to_csv(links_path, index=False, encoding="utf-8-sig")
    else:
        design.to_csv(design_path, index=False, encoding="utf-8-sig")
        links.to_csv(links_path, index=False, encoding="utf-8-sig")

    base = make_reference_parameters()
    reference_peak = reference_demand_peak(base, settings)
    reference_unit = np.array([
        (math.log(getattr(base, spec.name, spec.lower) / spec.lower)
         / math.log(spec.upper / spec.lower))
        if spec.scale == "log" and hasattr(base, spec.name)
        else 0.5
        for spec in PARAMETER_SPECS
    ], dtype=float)
    # Reference metrics are simulated at the documented baseline rather than at
    # the synthetic midpoint above.  The midpoint vector is unused otherwise.
    del reference_unit
    reference_physical = {
        "Td_min": 50.0,
        "oscillator_plasmid_copies": 10.0,
        "K_copies": 13.0,
        "hill_n": 3.0,
        "oscillator_promoter_leak": 0.0,
        "oscillator_tx_scale": 1.0,
        "oscillator_translation_scale": 1.0,
        "oscillator_mrna_total_half_life_min": 2.0,
        "sponge_plasmid_copies": 40.0,
        "sponge_accessible_fraction": 1.0,
        "c31_promoter_leak": 0.005,
        "c31_tx_scale": 1.0,
        "c31_translation_scale": 1.0,
        "c31_mrna_total_half_life_min": 2.0,
        "peak_load_fraction": 0.0,
    }
    reference_x = np.array([
        ((value - spec.lower) / (spec.upper - spec.lower)
         if spec.scale == "linear" else
         math.log(value / spec.lower) / math.log(spec.upper / spec.lower))
        for spec in PARAMETER_SPECS
        for value in [reference_physical[spec.name]]
    ])
    reference_row = simulate_run(-1, reference_x, base, settings, reference_peak)
    reference_scales = {
        metric: (max(abs(float(reference_row.get(metric, np.nan))), 1e-12)
                 if np.isfinite(float(reference_row.get(metric, np.nan))) else 1.0)
        for metric in METRICS
    }
    reference_scales.update(NORMALIZATION_SCALE_OVERRIDES)
    pd.DataFrame([reference_row]).to_csv(
        output_dir / "reference_metrics.csv", index=False, encoding="utf-8-sig"
    )

    checkpoint_path = output_dir / "run_metrics_checkpoint.csv"
    completed_rows: list[dict] = []
    completed_ids: set[int] = set()
    if checkpoint_path.exists():
        old = pd.read_csv(checkpoint_path)
        completed_rows = old.to_dict("records")
        completed_ids = set(old["run_id"].astype(int))
        print(f"Resuming with {len(completed_ids)} completed model evaluations")

    x_columns = [f"x__{spec.name}" for spec in PARAMETER_SPECS]
    total = len(design)
    for position, design_row in design.iterrows():
        run_id = int(design_row["run_id"])
        if run_id in completed_ids:
            continue
        print(f"[{position + 1:03d}/{total:03d}] trajectory={int(design_row['trajectory_id']):02d} "
              f"step={int(design_row['step_id']):02d} changed={design_row['changed_parameter']}")
        unit_values = design_row[x_columns].to_numpy(dtype=float)
        result = simulate_run(run_id, unit_values, base, settings, reference_peak)
        result.update({
            "trajectory_id": int(design_row["trajectory_id"]),
            "step_id": int(design_row["step_id"]),
            "changed_parameter": design_row["changed_parameter"],
        })
        completed_rows.append(result)
        completed_ids.add(run_id)
        write_checkpoint(completed_rows, checkpoint_path)

    metrics = pd.DataFrame(completed_rows).sort_values("run_id")
    metrics_path = output_dir / "morris_run_metrics.csv"
    metrics.to_csv(metrics_path, index=False, encoding="utf-8-sig")
    effects = compute_elementary_effects(metrics, links, reference_scales)
    effects.to_csv(output_dir / "morris_elementary_effects.csv", index=False,
                   encoding="utf-8-sig")
    summary = summarize_effects(effects, settings)
    summary.to_csv(output_dir / "morris_sensitivity_summary.csv", index=False,
                   encoding="utf-8-sig")

    top = (summary[summary["rank_within_metric"] <= 8]
           .sort_values(["metric", "rank_within_metric"]))
    top.to_csv(output_dir / "top_parameters_for_next_stage.csv", index=False,
               encoding="utf-8-sig")

    scope = {
        "version": VERSION,
        "mode": args.mode,
        "settings": asdict(settings),
        "morris_delta_on_unit_scale": delta,
        "model_evaluations": len(design),
        "reference_parameters": asdict(base),
        "reference_metrics": reference_row,
        "reference_C31_demand_peak_molecules_per_min": reference_peak,
        "C31_residual_flux_reporting_threshold_uM_h": B_INTERFACE_RESIDUAL_FLUX_THRESHOLD_uM_h,
        "Morris_output_normalization_overrides": NORMALIZATION_SCALE_OVERRIDES,
        "parameter_ranges": [asdict(spec) for spec in PARAMETER_SPECS],
        "metrics": list(METRICS),
        "method": (
            "Morris elementary-effects screening on normalized inputs; mu* ranks total "
            "importance, sigma flags nonlinearity/interactions; trajectory bootstrap gives CI"
        ),
        "resource_equation": "q31=beta31*m31/Jref; b=L/(1-L); rho=1/(1+b*q31)",
        "interpretation_boundaries": [
            "Ranges are uncertainty/design ranges, not fitted confidence intervals.",
            "A high Morris sigma suggests nonlinearity or interaction but does not identify its partner.",
            "Failed or unresolved oscillations are analyzed through resolved_oscillation; period metrics use valid pairs only.",
            "The symmetric common K, n and oscillator leak remain an approximation pending arm-specific transfer curves.",
            "peak_load_fraction is an uncalibrated boundary parameter, not the measured B burden.",
            "The 0.5 uM/h residual-flux threshold is an interface reporting convention, not a fitted B kinetic threshold.",
            "Conditional C31 residual metrics are defined only when the upstream clock is resolved; all-state metrics retain complete-system failure information.",
            "Smoke-mode rankings are not scientific results.",
        ],
    }
    with (output_dir / "scope_parameters_and_method.json").open("w", encoding="utf-8") as handle:
        json.dump(scope, handle, ensure_ascii=False, indent=2, allow_nan=True)

    if not args.no_plots:
        make_plots(summary, figure_dir)

    print("\nCompleted:", VERSION)
    print("Mode:", args.mode)
    print("Raw tables:", output_dir)
    print("Figures:", figure_dir)
    print("Resolved fraction:", f"{metrics['resolved_oscillation'].mean():.3f}")
    print("\nTop period-generation parameters:")
    period_top = summary[summary["metric"] == "period_generations"].head(8)
    for row in period_top.to_dict("records"):
        print(f"  {row['parameter']:<38s} mu*={row['mu_star_normalized']:.3g} "
              f"sigma={row['sigma_normalized']:.3g} valid={row['valid_fraction']:.2f}")
    print("\nInterpretation: use this screen to select 5-8 parameters for Sobol analysis; "
          "do not select a design from the largest or prettiest response alone.")


if __name__ == "__main__":
    main()
