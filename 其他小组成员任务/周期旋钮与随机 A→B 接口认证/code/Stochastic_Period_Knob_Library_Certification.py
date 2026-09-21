"""Stochastic certification of the deterministic period-knob library.

This program maps the explicit-mRNA mechanistic ODE to a stochastic reaction
network and tests the four deterministic A-to-B-certified knob settings.

Important scope boundary
------------------------
This is a non-negative binomial tau-leap approximation of the mechanistic ODE,
not the reduced Potvin-Trottier SI Gillespie model.  It is used to test whether
the designed 10--13 generation ordering and variability survive intrinsic
reaction and plasmid-copy noise.  The paper-aligned sponge/no-sponge phenotype
remains a separate analysis and is not re-fitted here.

Reactions represented explicitly:
* transcription and mRNA loss for TetR, CI, LacI and phiC31;
* translation and protein loss;
* mean-reverting birth/death fluctuations of oscillator and sponge plasmids;
* quasi-equilibrium TetR conservation over functional and sponge TetO sites.

No parameter is fitted by this program.  It reads the frozen certified library.

Run from the project root:
    python ./code/Stochastic_Period_Knob_Library_Certification.py --mode smoke
    python ./code/Stochastic_Period_Knob_Library_Certification.py --mode audit
    python ./code/Stochastic_Period_Knob_Library_Certification.py --mode confirmation

Use ``--fixed-plasmids`` only as a diagnostic that removes plasmid-copy noise.
"""

from __future__ import annotations

import argparse
import json
import math
import time
from dataclasses import asdict, dataclass, replace
from datetime import datetime
from functools import lru_cache
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import find_peaks, savgol_filter

import Mechanistic_ODE_Global_Sensitivity_Analysis as core


VERSION = "Stochastic_Period_Knob_Library_Certification_20260819"
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
OUTPUT_ROOT = PROJECT_ROOT / "outputs"
FIGURE_ROOT = PROJECT_ROOT / "figures"

TETR, CI, LACI = 0, 1, 2


@dataclass(frozen=True)
class RunSettings:
    trajectories_per_setting: int
    duration_generations: float
    transient_generations: float
    tau_step_min: float
    sample_step_min: float
    bootstrap_replicates: int
    random_seed: int = 20260819
    minimum_cycles: int = 4
    minimum_resolved_fraction: float = 0.80
    maximum_mean_period_relative_error: float = 0.20
    maximum_across_trajectory_cv: float = 0.30
    minimum_adjacent_order_probability: float = 0.80


MODE_SETTINGS = {
    "smoke": RunSettings(3, 90.0, 20.0, 0.20, 1.0, 500),
    "audit": RunSettings(12, 160.0, 30.0, 0.20, 1.0, 2000),
    "confirmation": RunSettings(50, 180.0, 30.0, 0.20, 1.0, 5000),
}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=tuple(MODE_SETTINGS), default="smoke")
    parser.add_argument("--library-csv", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument(
        "--fixed-plasmids",
        action="store_true",
        help="Diagnostic only: hold No and Nt fixed instead of adding copy-number noise.",
    )
    return parser.parse_args()


def locate_library(explicit: Path | None) -> Path:
    if explicit is not None:
        path = explicit.expanduser().resolve()
        if not path.exists():
            raise FileNotFoundError(path)
        return path
    candidates = sorted(
        OUTPUT_ROOT.glob(
            "Certified_Period_Knob_Library_*/period_knob_library_for_stochastic_validation.csv"
        ),
        key=lambda path: path.stat().st_mtime,
    )
    if not candidates:
        raise FileNotFoundError(
            "No certified period-knob library was found under outputs/. "
            "Run Consolidate_Certified_Period_Knob_Library.py first."
        )
    return candidates[-1]


def validate_library(frame: pd.DataFrame) -> pd.DataFrame:
    required = {
        "knob_setting_id",
        "target_period_generations",
        "realised_period_generations",
        "synchronous_oscillator_RBS_scale",
        "synchronous_oscillator_mRNA_half_life_min",
        "C31_RBS_scale",
        "C31_mRNA_total_half_life_min",
        "sponge_plasmid_copies",
        "sponge_accessible_fraction",
        "deterministic_A_pass",
        "deterministic_B_pass",
    }
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"Certified library is missing columns: {missing}")
    selected = frame.copy()
    selected = selected.loc[
        selected["deterministic_A_pass"].astype(bool)
        & selected["deterministic_B_pass"].astype(bool)
    ].sort_values("target_period_generations")
    targets = tuple(pd.to_numeric(selected["target_period_generations"]).astype(float))
    if targets != (10.0, 11.0, 12.0, 13.0):
        raise ValueError(f"Expected certified targets (10, 11, 12, 13); found {targets}")
    return selected.reset_index(drop=True)


def build_parameters(row: pd.Series) -> core.Parameters:
    base = core.make_reference_parameters()
    mu = math.log(2.0) / 50.0
    oscillator_half_life = float(row["synchronous_oscillator_mRNA_half_life_min"])
    c31_half_life = float(row["C31_mRNA_total_half_life_min"])
    return replace(
        base,
        Td_min=50.0,
        K_copies=13.0,
        hill_n=3.0,
        oscillator_plasmid_copies=10.0,
        functional_tetO_per_plasmid=4.0,
        functional_accessible_fraction=1.0,
        sponge_plasmid_copies=float(row["sponge_plasmid_copies"]),
        sponge_tetO_per_plasmid=2.0,
        sponge_accessible_fraction=float(row["sponge_accessible_fraction"]),
        # Final construct: cI and phiC31 use the same PLtetO1 curve.
        oscillator_promoter_leak=0.005,
        c31_promoter_leak=0.005,
        oscillator_translation_per_mrna_per_min=(
            base.oscillator_translation_per_mrna_per_min
            * float(row["synchronous_oscillator_RBS_scale"])
        ),
        oscillator_mrna_intrinsic_loss_per_min=max(
            0.0, math.log(2.0) / oscillator_half_life - mu
        ),
        c31_translation_per_mrna_per_min=(
            base.c31_translation_per_mrna_per_min * float(row["C31_RBS_scale"])
        ),
        c31_mrna_intrinsic_loss_per_min=max(
            0.0, math.log(2.0) / c31_half_life - mu
        ),
        peak_load_fraction=0.0,
    )


@lru_cache(maxsize=1_000_000)
def free_tetr_cached(
    total: int,
    oscillator_plasmids: int,
    sponge_plasmids: int,
    K: float,
    hill_n: float,
    functional_sites_per_plasmid: float,
    functional_accessibility: float,
    sponge_sites_per_plasmid: float,
    sponge_accessibility: float,
) -> float:
    total_value = float(max(int(total), 0))
    if total_value <= 0.0:
        return 0.0
    capacity = (
        functional_accessibility
        * functional_sites_per_plasmid
        * max(int(oscillator_plasmids), 0)
        + sponge_accessibility
        * sponge_sites_per_plasmid
        * max(int(sponge_plasmids), 0)
    )
    if capacity <= 0.0:
        return total_value
    low, high = 0.0, total_value
    K_n = K**hill_n
    free = max(total_value - 0.5 * capacity, 0.5 * total_value)
    for _ in range(35):
        free_n = max(free, 0.0) ** hill_n
        denominator = K_n + free_n
        residual = free + capacity * free_n / denominator - total_value
        if abs(residual) <= 1e-10 * max(total_value, 1.0):
            return float(free)
        if residual > 0.0:
            high = free
        else:
            low = free
        derivative = 1.0 + (
            capacity
            * hill_n
            * K_n
            * max(free, 1e-15) ** (hill_n - 1.0)
            / denominator**2
        )
        candidate = free - residual / derivative
        free = candidate if low < candidate < high else 0.5 * (low + high)
    return float(free)


def repression(repressor: float, K: float, hill_n: float, leak: float) -> float:
    value = max(float(repressor), 0.0)
    regulated = K**hill_n / (K**hill_n + value**hill_n)
    return leak + (1.0 - leak) * regulated


def initial_integer_state(p: core.Parameters) -> np.ndarray:
    continuous = np.asarray(core.initial_state(p), dtype=float)
    return np.maximum(np.rint(continuous), 0.0).astype(np.int64)


def binomial_losses(rng, counts: np.ndarray, rates: np.ndarray, dt: float):
    probabilities = 1.0 - np.exp(-np.maximum(rates, 0.0) * dt)
    probabilities = np.clip(probabilities, 0.0, 1.0)
    return rng.binomial(counts.astype(np.int64), probabilities)


def simulate_trajectory(
    row: pd.Series,
    settings: RunSettings,
    seed: int,
    fixed_plasmids: bool,
):
    p = build_parameters(row)
    dt = settings.tau_step_min
    sample_stride = int(round(settings.sample_step_min / dt))
    if sample_stride < 1 or not math.isclose(sample_stride * dt, settings.sample_step_min):
        raise ValueError("sample_step_min must be an integer multiple of tau_step_min")
    total_steps = int(round(settings.duration_generations * p.Td_min / dt))
    sample_count = total_steps // sample_stride + 1
    rng = np.random.default_rng(seed)
    state = initial_integer_state(p)
    No_mean = int(round(p.oscillator_plasmid_copies))
    Nt_mean = int(round(p.sponge_plasmid_copies))
    No, Nt = No_mean, Nt_mean

    times = np.empty(sample_count, dtype=float)
    proteins = np.empty((sample_count, 3), dtype=np.int32)
    free_tetr = np.empty(sample_count, dtype=float)
    c31_flux = np.empty(sample_count, dtype=float)
    plasmids = np.empty((sample_count, 2), dtype=np.int16)
    sample_index = 0

    mrna_loss = p.oscillator_mrna_intrinsic_loss_per_min + p.mu
    protein_loss = p.oscillator_protein_intrinsic_loss_per_min + p.mu
    c31_mrna_loss = p.c31_mrna_intrinsic_loss_per_min + p.mu
    c31_protein_loss = p.c31_protein_intrinsic_loss_per_min + p.mu
    loss_rates = np.array(
        [mrna_loss, protein_loss, mrna_loss, protein_loss,
         mrna_loss, protein_loss, c31_mrna_loss, c31_protein_loss],
        dtype=float,
    )

    for step in range(total_steps + 1):
        tet_free = free_tetr_cached(
            int(state[1]), No, Nt, p.K_copies, p.hill_n,
            p.functional_tetO_per_plasmid, p.functional_accessible_fraction,
            p.sponge_tetO_per_plasmid, p.sponge_accessible_fraction,
        )
        h_tet = repression(state[5], p.K_copies, p.hill_n, p.oscillator_promoter_leak)
        h_ci = repression(tet_free, p.K_copies, p.hill_n, p.oscillator_promoter_leak)
        h_laci = repression(state[3], p.K_copies, p.hill_n, p.oscillator_promoter_leak)
        h31 = repression(tet_free, p.K_copies, p.hill_n, p.c31_promoter_leak)

        if step % sample_stride == 0:
            times[sample_index] = step * dt
            proteins[sample_index] = (state[1], state[3], state[5])
            free_tetr[sample_index] = tet_free
            c31_flux[sample_index] = p.c31_translation_per_mrna_per_min * state[6]
            plasmids[sample_index] = (No, Nt)
            sample_index += 1
        if step == total_steps:
            break

        transcription_rates = np.array(
            [No * p.oscillator_tx_per_plasmid_per_min * h_tet,
             No * p.oscillator_tx_per_plasmid_per_min * h_ci,
             No * p.oscillator_tx_per_plasmid_per_min * h_laci,
             No * p.c31_tx_per_plasmid_per_min * h31],
            dtype=float,
        )
        translation_rates = np.array(
            [p.oscillator_translation_per_mrna_per_min * state[0],
             p.oscillator_translation_per_mrna_per_min * state[2],
             p.oscillator_translation_per_mrna_per_min * state[4],
             p.c31_translation_per_mrna_per_min * state[6]],
            dtype=float,
        )
        births_mrna = rng.poisson(np.maximum(transcription_rates * dt, 0.0))
        births_protein = rng.poisson(np.maximum(translation_rates * dt, 0.0))
        losses = binomial_losses(rng, state, loss_rates, dt)
        state -= losses
        state[[0, 2, 4, 6]] += births_mrna
        state[[1, 3, 5, 7]] += births_protein

        if not fixed_plasmids:
            No += int(rng.poisson(p.mu * No_mean * dt))
            Nt += int(rng.poisson(p.mu * Nt_mean * dt))
            if No > 0:
                No -= int(rng.binomial(No, 1.0 - math.exp(-p.mu * dt)))
            if Nt > 0:
                Nt -= int(rng.binomial(Nt, 1.0 - math.exp(-p.mu * dt)))

    return {
        "time_min": times[:sample_index],
        "proteins": proteins[:sample_index],
        "free_tetr": free_tetr[:sample_index],
        "c31_flux_molecules_per_min": c31_flux[:sample_index],
        "plasmids": plasmids[:sample_index],
        "seed": seed,
    }


def smooth(signal: np.ndarray, sample_step_min: float, window_min: float = 30.0):
    values = np.asarray(signal, dtype=float)
    window = max(5, int(round(window_min / sample_step_min)))
    if window % 2 == 0:
        window += 1
    maximum = len(values) if len(values) % 2 else len(values) - 1
    window = min(window, maximum)
    if window < 5:
        return values
    return savgol_filter(values, window, polyorder=2)


def valid_periods(periods_min: np.ndarray, deterministic_period_min: float):
    periods = np.asarray(periods_min, dtype=float)
    return periods[
        (periods >= 0.45 * deterministic_period_min)
        & (periods <= 1.80 * deterministic_period_min)
    ]


def cyclic_phase_periods(
    time_min: np.ndarray,
    proteins: np.ndarray,
    deterministic_period_min: float,
):
    x = np.asarray(proteins, dtype=float)
    x = np.column_stack([smooth(x[:, i], 1.0) for i in range(3)])
    scale = np.std(x, axis=0, ddof=1)
    scale[scale < 1e-9] = 1.0
    z = (x - np.mean(x, axis=0)) / scale
    roots = np.exp(2j * np.pi * np.arange(3) / 3.0)
    phase = np.unwrap(np.angle(z @ roots))
    if np.nanmedian(np.diff(phase)) < 0.0:
        phase = -phase
    monotone = np.maximum.accumulate(phase)
    start = 2.0 * np.pi * math.ceil(monotone[0] / (2.0 * np.pi))
    stop = 2.0 * np.pi * math.floor(monotone[-1] / (2.0 * np.pi))
    levels = np.arange(start, stop + np.pi, 2.0 * np.pi)
    crossings = []
    for level in levels:
        index = int(np.searchsorted(monotone, level, side="left"))
        if index <= 0 or index >= len(time_min):
            continue
        p0, p1 = monotone[index - 1], monotone[index]
        fraction = 0.0 if p1 <= p0 else (level - p0) / (p1 - p0)
        crossings.append(
            float(time_min[index - 1] + fraction * (time_min[index] - time_min[index - 1]))
        )
    crossings = np.asarray(crossings, dtype=float)
    return valid_periods(np.diff(crossings), deterministic_period_min), crossings, phase


def ci_peak_periods(
    time_min: np.ndarray,
    ci: np.ndarray,
    deterministic_period_min: float,
):
    signal = smooth(ci, 1.0)
    peaks, _ = find_peaks(
        signal,
        distance=max(1, int(round(0.45 * deterministic_period_min))),
        prominence=max(0.10 * float(np.ptp(signal)), 10.0),
    )
    periods = valid_periods(np.diff(time_min[peaks]), deterministic_period_min)
    return periods, time_min[peaks]


def coefficient_of_variation(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size < 2 or np.mean(values) == 0.0:
        return float("nan")
    return float(np.std(values, ddof=1) / np.mean(values))


def timing_jitter_generations(crossings_min: np.ndarray, Td_min: float):
    crossings = np.asarray(crossings_min, dtype=float)
    if crossings.size < 4:
        return float("nan")
    cycle = np.arange(crossings.size, dtype=float)
    fit = np.polyfit(cycle, crossings, deg=1)
    residual = crossings - np.polyval(fit, cycle)
    return float(np.sqrt(np.mean(residual**2)) / Td_min)


def analyse_trajectory(row: pd.Series, result: dict, settings: RunSettings):
    deterministic_gen = float(row["realised_period_generations"])
    deterministic_min = deterministic_gen * 50.0
    keep = result["time_min"] >= settings.transient_generations * 50.0
    time_min = result["time_min"][keep]
    proteins = result["proteins"][keep].astype(float)
    phase_periods, phase_crossings, phase = cyclic_phase_periods(
        time_min, proteins, deterministic_min
    )
    ci_periods, ci_peaks = ci_peak_periods(time_min, proteins[:, CI], deterministic_min)
    resolved = len(phase_periods) >= settings.minimum_cycles
    phase_gen = phase_periods / 50.0
    ci_gen = ci_periods / 50.0
    return {
        "resolved": bool(resolved),
        "phase_period_mean_generations": (
            float(np.mean(phase_gen)) if phase_gen.size else float("nan")
        ),
        "phase_period_within_trajectory_cv": coefficient_of_variation(phase_gen),
        "ci_peak_period_mean_generations": (
            float(np.mean(ci_gen)) if ci_gen.size else float("nan")
        ),
        "ci_peak_period_within_trajectory_cv": coefficient_of_variation(ci_gen),
        "phase_timing_jitter_generations": timing_jitter_generations(phase_crossings, 50.0),
        "n_phase_cycles": int(phase_gen.size),
        "n_ci_peak_cycles": int(ci_gen.size),
        "phase": phase,
        "time_min": time_min,
        "proteins": proteins,
        "phase_crossings_min": phase_crossings,
        "ci_peak_times_min": ci_peaks,
    }


def bootstrap_mean_interval(values, replicates, seed):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size < 2:
        return float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    means = np.empty(replicates, dtype=float)
    for index in range(replicates):
        means[index] = np.mean(rng.choice(values, size=values.size, replace=True))
    low, high = np.percentile(means, [2.5, 97.5])
    return float(low), float(high)


def summarize_setting(row, trajectory_rows, settings, seed):
    values = np.array(
        [item["phase_period_mean_generations"] for item in trajectory_rows], dtype=float
    )
    finite = values[np.isfinite(values)]
    mean = float(np.mean(finite)) if finite.size else float("nan")
    low, high = bootstrap_mean_interval(finite, settings.bootstrap_replicates, seed)
    resolved_fraction = float(np.mean([item["resolved"] for item in trajectory_rows]))
    deterministic = float(row["realised_period_generations"])
    relative_error = abs(mean / deterministic - 1.0) if np.isfinite(mean) else float("nan")
    across_cv = coefficient_of_variation(finite)
    individual_pass = bool(
        resolved_fraction >= settings.minimum_resolved_fraction
        and np.isfinite(relative_error)
        and relative_error <= settings.maximum_mean_period_relative_error
        and np.isfinite(across_cv)
        and across_cv <= settings.maximum_across_trajectory_cv
    )
    return {
        "knob_setting_id": row["knob_setting_id"],
        "target_period_generations": float(row["target_period_generations"]),
        "deterministic_period_generations": deterministic,
        "design_RBS_scale": float(row["synchronous_oscillator_RBS_scale"]),
        "design_oscillator_mRNA_half_life_min": float(
            row["synchronous_oscillator_mRNA_half_life_min"]
        ),
        "n_trajectories": len(trajectory_rows),
        "n_finite_trajectory_means": int(finite.size),
        "stochastic_resolved_fraction": resolved_fraction,
        "stochastic_mean_period_generations": mean,
        "stochastic_mean_period_CI_low": low,
        "stochastic_mean_period_CI_high": high,
        "stochastic_period_cv": across_cv,
        "mean_within_trajectory_period_cv": float(np.nanmean([
            item["phase_period_within_trajectory_cv"] for item in trajectory_rows
        ])),
        "mean_phase_timing_jitter_generations": float(np.nanmean([
            item["phase_timing_jitter_generations"] for item in trajectory_rows
        ])),
        "stochastic_vs_deterministic_relative_error": relative_error,
        "individual_stochastic_constraint_pass": individual_pass,
    }


def ordering_audit(summaries, trajectories_by_setting, settings):
    rows = []
    rng = np.random.default_rng(settings.random_seed + 8_000_000)
    for left, right in zip(summaries[:-1], summaries[1:]):
        left_values = np.array([
            item["phase_period_mean_generations"]
            for item in trajectories_by_setting[left["knob_setting_id"]]
        ], dtype=float)
        right_values = np.array([
            item["phase_period_mean_generations"]
            for item in trajectories_by_setting[right["knob_setting_id"]]
        ], dtype=float)
        left_values = left_values[np.isfinite(left_values)]
        right_values = right_values[np.isfinite(right_values)]
        differences = np.full(settings.bootstrap_replicates, np.nan)
        if left_values.size >= 2 and right_values.size >= 2:
            for index in range(settings.bootstrap_replicates):
                differences[index] = (
                    np.mean(rng.choice(right_values, size=right_values.size, replace=True))
                    - np.mean(rng.choice(left_values, size=left_values.size, replace=True))
                )
        finite = differences[np.isfinite(differences)]
        probability = float(np.mean(finite > 0.0)) if finite.size else float("nan")
        low, high = (
            tuple(np.percentile(finite, [2.5, 97.5]))
            if finite.size else (float("nan"), float("nan"))
        )
        rows.append({
            "lower_setting": left["knob_setting_id"],
            "upper_setting": right["knob_setting_id"],
            "mean_period_difference_generations": (
                right["stochastic_mean_period_generations"]
                - left["stochastic_mean_period_generations"]
            ),
            "difference_CI_low": float(low),
            "difference_CI_high": float(high),
            "bootstrap_probability_upper_period_is_longer": probability,
            "ordering_pass": bool(
                np.isfinite(probability)
                and probability >= settings.minimum_adjacent_order_probability
            ),
        })
    return rows


def make_figures(summary, trajectory_table, examples, figure_dir):
    figure_dir.mkdir(parents=True, exist_ok=True)
    x = np.arange(len(summary))
    deterministic = summary["deterministic_period_generations"].to_numpy(float)
    stochastic = summary["stochastic_mean_period_generations"].to_numpy(float)
    low = summary["stochastic_mean_period_CI_low"].to_numpy(float)
    high = summary["stochastic_mean_period_CI_high"].to_numpy(float)
    error = np.vstack([stochastic - low, high - stochastic])
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    axes[0].plot(x, deterministic, "o--", label="deterministic")
    axes[0].errorbar(x, stochastic, yerr=error, fmt="o-", capsize=4, label="stochastic")
    axes[0].set_xticks(x, summary["knob_setting_id"], rotation=20)
    axes[0].set_ylabel("period (generations)")
    axes[0].set_title("Period-knob ordering under stochastic reactions")
    axes[0].legend()
    axes[0].grid(alpha=0.25)
    axes[1].bar(x, summary["stochastic_period_cv"], color="tab:orange")
    axes[1].axhline(0.30, color="black", ls="--", label="provisional maximum")
    axes[1].set_xticks(x, summary["knob_setting_id"], rotation=20)
    axes[1].set_ylabel("CV across trajectory-mean periods")
    axes[1].set_title("Stochastic period variability")
    axes[1].legend()
    axes[1].grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(figure_dir / "Stochastic_period_knob_certification_summary.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(len(examples), 2, figsize=(13, 3.0 * len(examples)), squeeze=False)
    for row_index, (setting_id, example) in enumerate(examples.items()):
        time_gen = example["time_min"] / 50.0
        for protein_index, label in enumerate(("TetR", "CI", "LacI")):
            axes[row_index, 0].plot(time_gen, example["proteins"][:, protein_index], label=label)
        axes[row_index, 0].set_ylabel(f"{setting_id}\ncopies/cell")
        axes[row_index, 0].legend(fontsize=7)
        axes[row_index, 0].grid(alpha=0.2)
        axes[row_index, 1].plot(time_gen, example["phase"] / (2.0 * np.pi))
        axes[row_index, 1].set_ylabel("unwrapped phase / 2pi")
        axes[row_index, 1].grid(alpha=0.2)
    axes[-1, 0].set_xlabel("time (generations)")
    axes[-1, 1].set_xlabel("time (generations)")
    fig.suptitle("Example stochastic trajectories after transient")
    fig.tight_layout()
    fig.savefig(figure_dir / "Stochastic_period_knob_example_trajectories.png", dpi=180)
    plt.close(fig)


def main():
    args = parse_args()
    settings = MODE_SETTINGS[args.mode]
    library_path = locate_library(args.library_csv)
    library = validate_library(pd.read_csv(library_path))
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = args.output_dir or OUTPUT_ROOT / f"Stochastic_Period_Knob_{args.mode}_{stamp}"
    figure_dir = FIGURE_ROOT / output_dir.name
    output_dir.mkdir(parents=True, exist_ok=False)
    figure_dir.mkdir(parents=True, exist_ok=False)

    print(f"Loaded: {VERSION}")
    print(f"Mode: {args.mode}; settings={settings}")
    print(f"Certified library: {library_path}")
    print(f"Plasmid-copy noise: {not args.fixed_plasmids}")
    print("Primary period estimator: three-repressor cyclic phase")

    wall_start = time.perf_counter()
    trajectory_rows = []
    trajectories_by_setting = {}
    examples = {}
    summaries = []
    for setting_index, (_, row) in enumerate(library.iterrows()):
        setting_id = str(row["knob_setting_id"])
        trajectories_by_setting[setting_id] = []
        for trajectory_index in range(settings.trajectories_per_setting):
            seed = settings.random_seed + setting_index * 100_000 + trajectory_index
            print(
                f"[{setting_id}] trajectory {trajectory_index + 1}/"
                f"{settings.trajectories_per_setting}"
            )
            result = simulate_trajectory(row, settings, seed, args.fixed_plasmids)
            analysis = analyse_trajectory(row, result, settings)
            trajectories_by_setting[setting_id].append(analysis)
            if trajectory_index == 0:
                examples[setting_id] = analysis
            trajectory_rows.append({
                "knob_setting_id": setting_id,
                "target_period_generations": float(row["target_period_generations"]),
                "deterministic_period_generations": float(row["realised_period_generations"]),
                "design_RBS_scale": float(row["synchronous_oscillator_RBS_scale"]),
                "design_oscillator_mRNA_half_life_min": float(
                    row["synchronous_oscillator_mRNA_half_life_min"]
                ),
                "trajectory_index": trajectory_index,
                "seed": seed,
                **{key: value for key, value in analysis.items()
                   if key not in {"phase", "time_min", "proteins",
                                  "phase_crossings_min", "ci_peak_times_min"}},
            })
        summaries.append(summarize_setting(
            row,
            trajectories_by_setting[setting_id],
            settings,
            settings.random_seed + 9_000_000 + setting_index,
        ))

    ordering_rows = ordering_audit(summaries, trajectories_by_setting, settings)
    ordering_pass = bool(ordering_rows and all(item["ordering_pass"] for item in ordering_rows))
    for summary in summaries:
        summary["library_ordering_pass"] = ordering_pass
        summary["stochastic_validation_status"] = (
            "pass"
            if summary["individual_stochastic_constraint_pass"] and ordering_pass
            else "fail"
        )

    trajectory_table = pd.DataFrame(trajectory_rows)
    summary_table = pd.DataFrame(summaries)
    ordering_table = pd.DataFrame(ordering_rows)
    merge_table = summary_table[[
        "design_RBS_scale",
        "design_oscillator_mRNA_half_life_min",
        "stochastic_period_cv",
        "stochastic_mean_period_generations",
        "stochastic_resolved_fraction",
        "stochastic_validation_status",
    ]].copy()
    updated_library = library.merge(
        summary_table.drop(columns=["target_period_generations"]),
        on="knob_setting_id",
        how="left",
        suffixes=("", "__stochastic"),
    )

    trajectory_table.to_csv(output_dir / "stochastic_trajectory_metrics.csv", index=False)
    summary_table.to_csv(output_dir / "stochastic_setting_summary.csv", index=False)
    ordering_table.to_csv(output_dir / "stochastic_adjacent_ordering_audit.csv", index=False)
    merge_table.to_csv(output_dir / "stochastic_validation_for_designer_merge.csv", index=False)
    updated_library.to_csv(output_dir / "certified_library_with_stochastic_results.csv", index=False)
    make_figures(summary_table, trajectory_table, examples, figure_dir)

    metadata = {
        "version": VERSION,
        "mode": args.mode,
        "settings": asdict(settings),
        "library_path": str(library_path),
        "fixed_plasmids": bool(args.fixed_plasmids),
        "method": "non-negative binomial tau-leap",
        "primary_period_estimator": "three-repressor cyclic phase",
        "scope_boundary": (
            "Mechanistic period-knob stochastic certification; not the paper-SI "
            "sponge/no-sponge phenotype analysis."
        ),
        "wall_time_seconds": time.perf_counter() - wall_start,
    }
    (output_dir / "scope_and_provenance.json").write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    print(f"\nCompleted: {VERSION}")
    print(f"Raw tables: {output_dir}")
    print(f"Figures: {figure_dir}")
    print(f"Adjacent monotonic ordering pass: {ordering_pass}")
    print("\nStochastic period-knob certification:")
    for item in summaries:
        print(
            f"{item['knob_setting_id']} | "
            f"det={item['deterministic_period_generations']:.3f} gen | "
            f"stochastic={item['stochastic_mean_period_generations']:.3f} gen | "
            f"CV={item['stochastic_period_cv']:.3f} | "
            f"resolved={item['stochastic_resolved_fraction']:.2f} | "
            f"status={item['stochastic_validation_status']}"
        )


if __name__ == "__main__":
    main()
