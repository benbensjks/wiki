r"""Layer-by-layer mean-consistency audit for deterministic and stochastic A models.

Purpose
-------
The stochastic A-to-B tests produced substantially more phiC31 translation
flux than the matched deterministic ODE.  This script attributes that mean
shift without fitting, rescaling or changing any biological parameter.

For robust K4 and K5 it compares three conditions:

1. deterministic ODE reference;
2. stochastic tau-leap with fixed oscillator/sponge plasmid copy numbers;
3. stochastic tau-leap with the existing mean-reverting plasmid-copy noise.

The diagnostic records the quantities that earlier interface scripts did not
export: free TetR, instantaneous PLtetO1 activity, No*activity transcription
input, C31 mRNA and C31 translation flux.  It separates the flux ratio into
mean copy number, transfer-curve shift, nonlinear/Jensen contribution,
copy-activity covariance, mRNA balance and translation-scale contributions.

This is an attribution audit, not a parameter calibration and not an A-to-B
construct certification.  Do not use it to divide the stochastic flux by an
empirical correction factor.

Run from the project root:

    python .\code\Stochastic_Deterministic_Mean_Consistency_Audit.py --mode smoke
    python .\code\Stochastic_Deterministic_Mean_Consistency_Audit.py --mode audit

The user runs this program locally.  Creating this file does not run it.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import Mechanistic_ODE_Global_Sensitivity_Analysis as core
import Robust_Continuous_Period_Knob_AB_Certification as deterministic_a
import Robust_K4_K5_Stochastic_AB_Fidelity_Certification as library_tools
import Stochastic_Period_Knob_Library_Certification as stochastic_a


VERSION = "Stochastic_Deterministic_Mean_Consistency_Audit_20260823"
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
OUTPUT_ROOT = PROJECT_ROOT / "outputs"
FIGURE_ROOT = PROJECT_ROOT / "figures"


@dataclass(frozen=True)
class AuditSettings:
    trajectories_per_setting: int
    duration_generations: float
    transient_generations: float
    tau_step_min: float
    sample_step_min: float
    random_seed: int = 20260823


MODE_SETTINGS = {
    "smoke": AuditSettings(2, 100.0, 25.0, 0.20, 1.0),
    "audit": AuditSettings(12, 150.0, 30.0, 0.20, 1.0),
}


PARAMETER_PARITY_FIELDS = (
    "Td_min",
    "K_copies",
    "hill_n",
    "oscillator_plasmid_copies",
    "functional_tetO_per_plasmid",
    "functional_accessible_fraction",
    "sponge_plasmid_copies",
    "sponge_tetO_per_plasmid",
    "sponge_accessible_fraction",
    "oscillator_promoter_leak",
    "c31_promoter_leak",
    "oscillator_tx_per_plasmid_per_min",
    "oscillator_translation_per_mrna_per_min",
    "oscillator_mrna_intrinsic_loss_per_min",
    "c31_tx_per_plasmid_per_min",
    "c31_translation_per_mrna_per_min",
    "c31_mrna_intrinsic_loss_per_min",
    "peak_load_fraction",
)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=tuple(MODE_SETTINGS), default="smoke")
    parser.add_argument("--library-csv", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=None)
    return parser.parse_args()


def stochastic_run_settings(settings: AuditSettings) -> stochastic_a.RunSettings:
    return stochastic_a.RunSettings(
        trajectories_per_setting=settings.trajectories_per_setting,
        duration_generations=settings.duration_generations,
        transient_generations=settings.transient_generations,
        tau_step_min=settings.tau_step_min,
        sample_step_min=settings.sample_step_min,
        bootstrap_replicates=500,
        random_seed=settings.random_seed,
        minimum_cycles=4,
    )


def parameter_parity_rows(setting_id, deterministic_p, stochastic_p):
    rows = []
    for field in PARAMETER_PARITY_FIELDS:
        deterministic_value = float(getattr(deterministic_p, field))
        stochastic_value = float(getattr(stochastic_p, field))
        scale = max(abs(deterministic_value), abs(stochastic_value), 1.0)
        relative_difference = abs(deterministic_value - stochastic_value) / scale
        rows.append(
            {
                "knob_setting_id": setting_id,
                "parameter": field,
                "deterministic_value": deterministic_value,
                "stochastic_value": stochastic_value,
                "relative_difference": relative_difference,
                "parity_pass": relative_difference <= 1e-10,
            }
        )
    return rows


def simulate_stochastic_detailed(
    library_row: pd.Series,
    settings: stochastic_a.RunSettings,
    seed: int,
    fixed_plasmids: bool,
):
    """Use the frozen tau-leap network while exporting intermediate means."""
    p = stochastic_a.build_parameters(library_row)
    dt = settings.tau_step_min
    sample_stride = int(round(settings.sample_step_min / dt))
    if sample_stride < 1 or not math.isclose(
        sample_stride * dt, settings.sample_step_min, rel_tol=0.0, abs_tol=1e-12
    ):
        raise ValueError("sample_step_min must be an integer multiple of tau_step_min")

    total_steps = int(round(settings.duration_generations * p.Td_min / dt))
    sample_count = total_steps // sample_stride + 1
    rng = np.random.default_rng(seed)
    state = stochastic_a.initial_integer_state(p)
    no_mean = int(round(p.oscillator_plasmid_copies))
    nt_mean = int(round(p.sponge_plasmid_copies))
    no, nt = no_mean, nt_mean

    time_min = np.empty(sample_count, dtype=float)
    states = np.empty((sample_count, 8), dtype=float)
    free_tetr = np.empty(sample_count, dtype=float)
    h31 = np.empty(sample_count, dtype=float)
    tx_input = np.empty(sample_count, dtype=float)
    c31_flux = np.empty(sample_count, dtype=float)
    plasmids = np.empty((sample_count, 2), dtype=float)
    sample_index = 0

    mrna_loss = p.oscillator_mrna_intrinsic_loss_per_min + p.mu
    protein_loss = p.oscillator_protein_intrinsic_loss_per_min + p.mu
    c31_mrna_loss = p.c31_mrna_intrinsic_loss_per_min + p.mu
    c31_protein_loss = p.c31_protein_intrinsic_loss_per_min + p.mu
    loss_rates = np.array(
        [
            mrna_loss,
            protein_loss,
            mrna_loss,
            protein_loss,
            mrna_loss,
            protein_loss,
            c31_mrna_loss,
            c31_protein_loss,
        ],
        dtype=float,
    )

    for step in range(total_steps + 1):
        tet_free = stochastic_a.free_tetr_cached(
            int(state[1]),
            no,
            nt,
            p.K_copies,
            p.hill_n,
            p.functional_tetO_per_plasmid,
            p.functional_accessible_fraction,
            p.sponge_tetO_per_plasmid,
            p.sponge_accessible_fraction,
        )
        h_tet = stochastic_a.repression(
            state[5], p.K_copies, p.hill_n, p.oscillator_promoter_leak
        )
        h_ci = stochastic_a.repression(
            tet_free, p.K_copies, p.hill_n, p.oscillator_promoter_leak
        )
        h_laci = stochastic_a.repression(
            state[3], p.K_copies, p.hill_n, p.oscillator_promoter_leak
        )
        local_h31 = stochastic_a.repression(
            tet_free, p.K_copies, p.hill_n, p.c31_promoter_leak
        )

        if step % sample_stride == 0:
            time_min[sample_index] = step * dt
            states[sample_index] = state
            free_tetr[sample_index] = tet_free
            h31[sample_index] = local_h31
            tx_input[sample_index] = no * p.c31_tx_per_plasmid_per_min * local_h31
            c31_flux[sample_index] = p.c31_translation_per_mrna_per_min * state[6]
            plasmids[sample_index] = (no, nt)
            sample_index += 1

        if step == total_steps:
            break

        transcription_rates = np.array(
            [
                no * p.oscillator_tx_per_plasmid_per_min * h_tet,
                no * p.oscillator_tx_per_plasmid_per_min * h_ci,
                no * p.oscillator_tx_per_plasmid_per_min * h_laci,
                no * p.c31_tx_per_plasmid_per_min * local_h31,
            ],
            dtype=float,
        )
        translation_rates = np.array(
            [
                p.oscillator_translation_per_mrna_per_min * state[0],
                p.oscillator_translation_per_mrna_per_min * state[2],
                p.oscillator_translation_per_mrna_per_min * state[4],
                p.c31_translation_per_mrna_per_min * state[6],
            ],
            dtype=float,
        )
        births_mrna = rng.poisson(np.maximum(transcription_rates * dt, 0.0))
        births_protein = rng.poisson(np.maximum(translation_rates * dt, 0.0))
        losses = stochastic_a.binomial_losses(rng, state, loss_rates, dt)
        state -= losses
        state[[0, 2, 4, 6]] += births_mrna
        state[[1, 3, 5, 7]] += births_protein

        if not fixed_plasmids:
            no += int(rng.poisson(p.mu * no_mean * dt))
            nt += int(rng.poisson(p.mu * nt_mean * dt))
            if no > 0:
                no -= int(rng.binomial(no, 1.0 - math.exp(-p.mu * dt)))
            if nt > 0:
                nt -= int(rng.binomial(nt, 1.0 - math.exp(-p.mu * dt)))

    use = slice(0, sample_index)
    return {
        "parameters": p,
        "time_min": time_min[use],
        "states": states[use],
        "free_tetr": free_tetr[use],
        "h31": h31[use],
        "tx_input_mrna_per_min": tx_input[use],
        "c31_flux_molecules_per_min": c31_flux[use],
        "plasmids": plasmids[use],
        "seed": seed,
    }


def deterministic_detailed(robust_row: pd.Series, settings: AuditSettings):
    p, _ = deterministic_a.parameters_from_design(robust_row, 0.0)
    duration_h = settings.duration_generations * p.Td_min / 60.0
    reference_peak = deterministic_a.reference_demand_peak(
        core.make_reference_parameters(), duration_h
    )
    time_h, solution, tet_free, resource, _ = deterministic_a.simulate_a(
        p, reference_peak, duration_h
    )
    h31 = np.array(
        [core.repression(value, p, p.c31_promoter_leak) for value in tet_free],
        dtype=float,
    )
    no = np.full_like(time_h, p.oscillator_plasmid_copies, dtype=float)
    nt = np.full_like(time_h, p.sponge_plasmid_copies, dtype=float)
    return {
        "parameters": p,
        "time_min": time_h * 60.0,
        "states": solution.y.T,
        "free_tetr": tet_free,
        "h31": h31,
        "tx_input_mrna_per_min": no * p.c31_tx_per_plasmid_per_min * h31,
        "c31_flux_molecules_per_min": (
            resource * p.c31_translation_per_mrna_per_min * solution.y[6]
        ),
        "plasmids": np.column_stack((no, nt)),
        "seed": -1,
    }


def safe_ratio(numerator: float, denominator: float) -> float:
    return float(numerator / denominator) if abs(denominator) > 1e-15 else math.nan


def summarize_trace(
    setting_id: str,
    condition: str,
    trajectory_index: int,
    trace: dict,
    transient_generations: float,
):
    p = trace["parameters"]
    keep = trace["time_min"] >= transient_generations * p.Td_min
    if np.count_nonzero(keep) < 10:
        raise RuntimeError("Too few post-transient samples for mean audit")

    free = np.asarray(trace["free_tetr"][keep], dtype=float)
    h31 = np.asarray(trace["h31"][keep], dtype=float)
    tx_input = np.asarray(trace["tx_input_mrna_per_min"][keep], dtype=float)
    flux = np.asarray(trace["c31_flux_molecules_per_min"][keep], dtype=float)
    states = np.asarray(trace["states"][keep], dtype=float)
    plasmids = np.asarray(trace["plasmids"][keep], dtype=float)
    no = plasmids[:, 0]
    nt = plasmids[:, 1]
    m31 = states[:, 6]
    tet_total = states[:, 1]

    mean_free = float(np.mean(free))
    mean_h = float(np.mean(h31))
    h_at_mean_free = stochastic_a.repression(
        mean_free, p.K_copies, p.hill_n, p.c31_promoter_leak
    )
    mean_no = float(np.mean(no))
    mean_nt = float(np.mean(nt))
    mean_m31 = float(np.mean(m31))
    mean_tx = float(np.mean(tx_input))
    mean_flux = float(np.mean(flux))
    c31_mrna_loss = p.c31_mrna_intrinsic_loss_per_min + p.mu
    expected_m31_from_mean_input = mean_tx / c31_mrna_loss
    no_h_coupling = safe_ratio(float(np.mean(no * h31)), mean_no * mean_h)
    no_h_correlation = (
        float(np.corrcoef(no, h31)[0, 1])
        if np.std(no) > 0.0 and np.std(h31) > 0.0
        else 0.0
    )

    return {
        "knob_setting_id": setting_id,
        "condition": condition,
        "trajectory_index": trajectory_index,
        "seed": int(trace["seed"]),
        "n_post_transient_samples": int(np.count_nonzero(keep)),
        "mean_oscillator_plasmids": mean_no,
        "cv_oscillator_plasmids": safe_ratio(float(np.std(no, ddof=1)), mean_no),
        "mean_sponge_plasmids": mean_nt,
        "cv_sponge_plasmids": safe_ratio(float(np.std(nt, ddof=1)), mean_nt),
        "mean_TetR_total_copies": float(np.mean(tet_total)),
        "mean_TetR_free_copies": mean_free,
        "mean_PLtetO1_activity": mean_h,
        "PLtetO1_activity_at_mean_free_TetR": h_at_mean_free,
        "Jensen_factor_Eh_over_hEfree": safe_ratio(mean_h, h_at_mean_free),
        "No_activity_coupling_factor": no_h_coupling,
        "No_activity_correlation": no_h_correlation,
        "mean_C31_transcription_input_mrna_per_min": mean_tx,
        "mean_C31_mRNA_copies": mean_m31,
        "C31_mRNA_expected_from_mean_input": expected_m31_from_mean_input,
        "C31_mRNA_balance_factor": safe_ratio(mean_m31, expected_m31_from_mean_input),
        "mean_C31_flux_molecules_per_min": mean_flux,
        "mean_C31_flux_uM_h": float(library_tools.flux_to_uM_h(np.array([mean_flux]))[0]),
        "C31_translation_per_mrna_per_min": p.c31_translation_per_mrna_per_min,
        "C31_mRNA_total_loss_per_min": c31_mrna_loss,
    }


def aggregate_trajectory_rows(frame: pd.DataFrame) -> pd.DataFrame:
    numeric = [
        column
        for column in frame.columns
        if column not in {"knob_setting_id", "condition"}
        and pd.api.types.is_numeric_dtype(frame[column])
        and column not in {"trajectory_index", "seed"}
    ]
    rows = []
    for (setting_id, condition), group in frame.groupby(
        ["knob_setting_id", "condition"], sort=True
    ):
        row = {
            "knob_setting_id": setting_id,
            "condition": condition,
            "n_trajectories": int(group["trajectory_index"].nunique()),
        }
        for column in numeric:
            values = pd.to_numeric(group[column], errors="coerce")
            row[column] = float(values.mean())
            row[f"{column}_trajectory_sd"] = float(values.std(ddof=1)) if len(values) > 1 else 0.0
        rows.append(row)
    return pd.DataFrame(rows)


def decomposition_rows(summary: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for setting_id, local in summary.groupby("knob_setting_id", sort=True):
        lookup = local.set_index("condition")
        deterministic = lookup.loc["deterministic_reference"]
        for condition in ("stochastic_fixed_plasmids", "stochastic_fluctuating_plasmids"):
            stochastic = lookup.loc[condition]

            copy_ratio = safe_ratio(
                stochastic["mean_oscillator_plasmids"],
                deterministic["mean_oscillator_plasmids"],
            )
            transfer_at_mean_ratio = safe_ratio(
                stochastic["PLtetO1_activity_at_mean_free_TetR"],
                deterministic["PLtetO1_activity_at_mean_free_TetR"],
            )
            jensen_ratio = safe_ratio(
                stochastic["Jensen_factor_Eh_over_hEfree"],
                deterministic["Jensen_factor_Eh_over_hEfree"],
            )
            coupling_ratio = safe_ratio(
                stochastic["No_activity_coupling_factor"],
                deterministic["No_activity_coupling_factor"],
            )
            loss_ratio = safe_ratio(
                deterministic["C31_mRNA_total_loss_per_min"],
                stochastic["C31_mRNA_total_loss_per_min"],
            )
            balance_ratio = safe_ratio(
                stochastic["C31_mRNA_balance_factor"],
                deterministic["C31_mRNA_balance_factor"],
            )
            beta_ratio = safe_ratio(
                stochastic["C31_translation_per_mrna_per_min"],
                deterministic["C31_translation_per_mrna_per_min"],
            )
            observed_flux_ratio = safe_ratio(
                stochastic["mean_C31_flux_molecules_per_min"],
                deterministic["mean_C31_flux_molecules_per_min"],
            )
            predicted_flux_ratio = float(
                np.prod(
                    [
                        copy_ratio,
                        transfer_at_mean_ratio,
                        jensen_ratio,
                        coupling_ratio,
                        loss_ratio,
                        balance_ratio,
                        beta_ratio,
                    ]
                )
            )
            rows.append(
                {
                    "knob_setting_id": setting_id,
                    "condition": condition,
                    "observed_flux_ratio_vs_deterministic": observed_flux_ratio,
                    "predicted_flux_ratio_from_factors": predicted_flux_ratio,
                    "closure_relative_error": safe_ratio(
                        predicted_flux_ratio - observed_flux_ratio,
                        observed_flux_ratio,
                    ),
                    "mean_plasmid_copy_ratio": copy_ratio,
                    "transfer_at_mean_free_TetR_ratio": transfer_at_mean_ratio,
                    "nonlinear_Jensen_ratio": jensen_ratio,
                    "copy_activity_coupling_ratio": coupling_ratio,
                    "mRNA_loss_ratio": loss_ratio,
                    "mRNA_balance_ratio": balance_ratio,
                    "translation_parameter_ratio": beta_ratio,
                    "mean_free_TetR_ratio": safe_ratio(
                        stochastic["mean_TetR_free_copies"],
                        deterministic["mean_TetR_free_copies"],
                    ),
                    "mean_PLtetO1_activity_ratio": safe_ratio(
                        stochastic["mean_PLtetO1_activity"],
                        deterministic["mean_PLtetO1_activity"],
                    ),
                    "mean_C31_mRNA_ratio": safe_ratio(
                        stochastic["mean_C31_mRNA_copies"],
                        deterministic["mean_C31_mRNA_copies"],
                    ),
                }
            )
    return pd.DataFrame(rows)


def make_figures(summary, decomposition, examples, figure_dir: Path):
    setting_ids = tuple(sorted(summary["knob_setting_id"].unique()))
    condition_order = (
        "deterministic_reference",
        "stochastic_fixed_plasmids",
        "stochastic_fluctuating_plasmids",
    )
    condition_labels = ("ODE", "SSA fixed plasmids", "SSA fluctuating plasmids")

    fig, axes = plt.subplots(2, len(setting_ids), figsize=(6 * len(setting_ids), 8), squeeze=False)
    for column, setting_id in enumerate(setting_ids):
        local = summary[summary["knob_setting_id"].eq(setting_id)].set_index("condition")
        local = local.loc[list(condition_order)]
        axes[0, column].bar(condition_labels, local["mean_C31_flux_uM_h"])
        axes[0, column].set_title(f"{setting_id}: mean C31 flux")
        axes[0, column].set_ylabel("uM/h")
        axes[0, column].tick_params(axis="x", rotation=25)
        axes[1, column].bar(condition_labels, local["mean_PLtetO1_activity"])
        axes[1, column].set_title(f"{setting_id}: mean PLtetO1 activity")
        axes[1, column].set_ylabel("fraction")
        axes[1, column].tick_params(axis="x", rotation=25)
    fig.suptitle("Deterministic-stochastic mean consistency")
    fig.tight_layout()
    fig.savefig(figure_dir / "Mean_consistency_layer_summary.png", dpi=200)
    plt.close(fig)

    factor_columns = (
        "mean_plasmid_copy_ratio",
        "transfer_at_mean_free_TetR_ratio",
        "nonlinear_Jensen_ratio",
        "copy_activity_coupling_ratio",
        "mRNA_balance_ratio",
        "translation_parameter_ratio",
    )
    fig, axes = plt.subplots(1, len(setting_ids), figsize=(7 * len(setting_ids), 5), squeeze=False)
    for column, setting_id in enumerate(setting_ids):
        local = decomposition[decomposition["knob_setting_id"].eq(setting_id)].set_index(
            "condition"
        )
        x = np.arange(len(factor_columns))
        width = 0.36
        for offset, condition in enumerate(
            ("stochastic_fixed_plasmids", "stochastic_fluctuating_plasmids")
        ):
            values = np.log2(np.maximum(local.loc[condition, list(factor_columns)], 1e-12))
            axes[0, column].bar(
                x + (offset - 0.5) * width,
                values,
                width,
                label="fixed" if offset == 0 else "fluctuating",
            )
        axes[0, column].axhline(0.0, color="black", linewidth=1)
        axes[0, column].set_xticks(x, [name.replace("_ratio", "") for name in factor_columns], rotation=35, ha="right")
        axes[0, column].set_ylabel("log2 contribution ratio")
        axes[0, column].set_title(setting_id)
        axes[0, column].legend()
    fig.suptitle("Multiplicative attribution of stochastic/ODE C31 flux")
    fig.tight_layout()
    fig.savefig(figure_dir / "Mean_flux_ratio_decomposition.png", dpi=200)
    plt.close(fig)

    if examples:
        fig, axes = plt.subplots(len(setting_ids), 2, figsize=(13, 4 * len(setting_ids)), squeeze=False)
        for row, setting_id in enumerate(setting_ids):
            for condition, color in (
                ("deterministic_reference", "black"),
                ("stochastic_fixed_plasmids", "tab:blue"),
                ("stochastic_fluctuating_plasmids", "tab:orange"),
            ):
                trace = examples[(setting_id, condition)]
                p = trace["parameters"]
                show = trace["time_min"] >= max(trace["time_min"][-1] - 25.0 * 60.0, 0.0)
                time_h = trace["time_min"][show] / 60.0
                axes[row, 0].plot(time_h, trace["h31"][show], color=color, label=condition)
                flux = library_tools.flux_to_uM_h(
                    trace["c31_flux_molecules_per_min"][show]
                )
                axes[row, 1].plot(time_h, flux, color=color, label=condition)
            axes[row, 0].set_title(f"{setting_id}: PLtetO1 activity")
            axes[row, 1].set_title(f"{setting_id}: C31 translation flux")
            axes[row, 0].set_ylabel("activity")
            axes[row, 1].set_ylabel("uM/h")
            axes[row, 0].legend(fontsize=8)
            axes[row, 1].legend(fontsize=8)
        axes[-1, 0].set_xlabel("time (h)")
        axes[-1, 1].set_xlabel("time (h)")
        fig.tight_layout()
        fig.savefig(figure_dir / "Mean_consistency_representative_traces.png", dpi=200)
        plt.close(fig)


def main():
    args = parse_args()
    settings = MODE_SETTINGS[args.mode]
    stochastic_settings = stochastic_run_settings(settings)
    library_path, library = library_tools.load_k4_k5_library(args.library_csv)
    library = library_tools.validate_k4_k5_library(library)
    robust_source = deterministic_a.load_representatives(library_path)
    robust_lookup = robust_source.set_index("setting_id")

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = (
        args.output_dir.expanduser().resolve()
        if args.output_dir is not None
        else OUTPUT_ROOT / f"Mean_Consistency_Audit_{args.mode}_{stamp}"
    )
    figure_dir = FIGURE_ROOT / output_dir.name
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    print(f"{VERSION}: mode={args.mode}")
    print(f"Robust library: {library_path}")
    print("Conditions: deterministic, stochastic fixed plasmids, stochastic fluctuating plasmids")

    trajectory_rows = []
    parity_rows = []
    examples = {}

    for setting_index, (_, library_row) in enumerate(library.iterrows()):
        setting_id = str(library_row["knob_setting_id"])
        robust_row = robust_lookup.loc[str(library_row["source_setting_id"])]
        deterministic_trace = deterministic_detailed(robust_row, settings)
        stochastic_parameters = stochastic_a.build_parameters(library_row)
        parity_rows.extend(
            parameter_parity_rows(
                setting_id, deterministic_trace["parameters"], stochastic_parameters
            )
        )
        trajectory_rows.append(
            summarize_trace(
                setting_id,
                "deterministic_reference",
                -1,
                deterministic_trace,
                settings.transient_generations,
            )
        )
        examples[(setting_id, "deterministic_reference")] = deterministic_trace

        for fixed_plasmids, condition in (
            (True, "stochastic_fixed_plasmids"),
            (False, "stochastic_fluctuating_plasmids"),
        ):
            for trajectory_index in range(settings.trajectories_per_setting):
                seed = (
                    settings.random_seed
                    + setting_index * 100_000
                    + (0 if fixed_plasmids else 50_000)
                    + trajectory_index
                )
                print(
                    f"[{setting_id} | {condition} | "
                    f"{trajectory_index + 1:02d}/{settings.trajectories_per_setting:02d}]"
                )
                trace = simulate_stochastic_detailed(
                    library_row, stochastic_settings, seed, fixed_plasmids
                )
                trajectory_rows.append(
                    summarize_trace(
                        setting_id,
                        condition,
                        trajectory_index,
                        trace,
                        settings.transient_generations,
                    )
                )
                if trajectory_index == 0:
                    examples[(setting_id, condition)] = trace

    parity = pd.DataFrame(parity_rows)
    if not parity["parity_pass"].all():
        mismatch = parity.loc[~parity["parity_pass"]]
        mismatch.to_csv(output_dir / "parameter_mapping_mismatches.csv", index=False)
        raise RuntimeError(
            "Deterministic and stochastic parameter mappings differ. "
            "See parameter_mapping_mismatches.csv; attribution was stopped."
        )

    trajectories = pd.DataFrame(trajectory_rows)
    summary = aggregate_trajectory_rows(trajectories)
    decomposition = decomposition_rows(summary)

    trajectories.to_csv(output_dir / "mean_consistency_by_trajectory.csv", index=False)
    summary.to_csv(output_dir / "mean_consistency_condition_summary.csv", index=False)
    decomposition.to_csv(output_dir / "mean_flux_ratio_decomposition.csv", index=False)
    parity.to_csv(output_dir / "deterministic_stochastic_parameter_parity.csv", index=False)

    metadata = {
        "version": VERSION,
        "mode": args.mode,
        "settings": asdict(settings),
        "library_csv": str(library_path),
        "conditions": [
            "deterministic_reference",
            "stochastic_fixed_plasmids",
            "stochastic_fluctuating_plasmids",
        ],
        "primary_question": (
            "Why is mean stochastic C31 translation flux higher than the matched ODE?"
        ),
        "scope_boundary": (
            "Attribution only. No parameter fitting, empirical normalization or B certification."
        ),
        "fixed_plasmids_recorded_explicitly": True,
    }
    (output_dir / "scope_and_settings.json").write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    make_figures(summary, decomposition, examples, figure_dir)

    print(f"\nCompleted: {VERSION}")
    print(f"Raw tables: {output_dir}")
    print(f"Figures: {figure_dir}")
    print("\nMean stochastic/ODE C31 flux attribution:")
    for _, row in decomposition.iterrows():
        print(
            f"{row['knob_setting_id']:>9s} | {row['condition']:<33s} | "
            f"flux={row['observed_flux_ratio_vs_deterministic']:.3f}x | "
            f"transfer={row['transfer_at_mean_free_TetR_ratio']:.3f}x | "
            f"Jensen={row['nonlinear_Jensen_ratio']:.3f}x | "
            f"copy-coupling={row['copy_activity_coupling_ratio']:.3f}x | "
            f"mRNA-balance={row['mRNA_balance_ratio']:.3f}x"
        )
    print(
        "\nRead the decomposition before changing alpha31, beta31 or leak. "
        "This run does not certify an empirical flux correction."
    )


if __name__ == "__main__":
    main()
