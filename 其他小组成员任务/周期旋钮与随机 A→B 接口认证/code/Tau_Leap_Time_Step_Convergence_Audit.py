r"""Time-step convergence audit for the aligned stochastic repressilator-C31 model.

This program follows the deterministic/stochastic mean-consistency audit.  It
does not fit biological parameters.  It tests whether the remaining small mean
shift is the expected finite-step bias of the tau-leap birth/death update.

For robust K4 and K5, with both fixed and fluctuating plasmid copies, it compares
tau steps 0.20, 0.10 and 0.05 min while holding the one-minute output sampling
grid, duration, random-seed scheme and biological parameters fixed.

Primary readouts:

* mean stochastic/ODE C31 translation-flux ratio;
* C31 mRNA balance factor;
* theoretical finite-step immigration/death factor;
* phase-period mean, within-trajectory CV and resolved fraction;
* agreement of the two smallest time steps.

Run from the project root:

    python .\code\Tau_Leap_Time_Step_Convergence_Audit.py --mode smoke
    python .\code\Tau_Leap_Time_Step_Convergence_Audit.py --mode audit

Smoke mode checks workflow only.  Report conclusions from audit mode.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass, replace
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import Robust_Continuous_Period_Knob_AB_Certification as deterministic_a
import Robust_K4_K5_Stochastic_AB_Fidelity_Certification as library_tools
import Stochastic_Deterministic_Mean_Consistency_Audit as mean_audit
import Stochastic_Period_Knob_Library_Certification as stochastic_a


VERSION = "Tau_Leap_Time_Step_Convergence_Audit_v1p1_20260823"
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
OUTPUT_ROOT = PROJECT_ROOT / "outputs"
FIGURE_ROOT = PROJECT_ROOT / "figures"


@dataclass(frozen=True)
class ConvergenceSettings:
    trajectories_per_setting: int
    duration_generations: float
    transient_generations: float
    sample_step_min: float
    tau_steps_min: tuple[float, ...]
    random_seed: int = 20260823
    maximum_small_step_relative_change: float = 0.02
    maximum_small_step_mrna_balance_error: float = 0.015
    minimum_resolved_fraction: float = 0.80


MODE_SETTINGS = {
    "smoke": ConvergenceSettings(2, 90.0, 20.0, 1.0, (0.20, 0.10, 0.05)),
    "audit": ConvergenceSettings(8, 130.0, 25.0, 1.0, (0.20, 0.10, 0.05)),
}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=tuple(MODE_SETTINGS), default="smoke")
    parser.add_argument("--library-csv", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument(
        "--tau-steps-min",
        default=None,
        help="Optional comma-separated steps in min; default 0.20,0.10,0.05.",
    )
    return parser.parse_args()


def parse_steps(text: str | None, default: tuple[float, ...]) -> tuple[float, ...]:
    if text is None:
        values = default
    else:
        values = tuple(float(item.strip()) for item in text.split(",") if item.strip())
    values = tuple(sorted(set(values), reverse=True))
    if len(values) < 2 or any(value <= 0.0 for value in values):
        raise ValueError("At least two positive tau steps are required")
    return values


def make_run_settings(settings: ConvergenceSettings, tau_step: float):
    ratio = settings.sample_step_min / tau_step
    if not math.isclose(ratio, round(ratio), rel_tol=0.0, abs_tol=1e-10):
        raise ValueError(
            f"sample_step_min={settings.sample_step_min} must be divisible by tau={tau_step}"
        )
    return stochastic_a.RunSettings(
        trajectories_per_setting=settings.trajectories_per_setting,
        duration_generations=settings.duration_generations,
        transient_generations=settings.transient_generations,
        tau_step_min=tau_step,
        sample_step_min=settings.sample_step_min,
        bootstrap_replicates=500,
        random_seed=settings.random_seed,
        minimum_cycles=4,
    )


def theoretical_birth_death_factor(loss_rate_per_min: float, tau_step_min: float) -> float:
    value = loss_rate_per_min * tau_step_min
    if value <= 1e-12:
        return 1.0
    return float(value / (1.0 - math.exp(-value)))


def stochastic_analysis_payload(trace: dict) -> dict:
    return {
        "time_min": trace["time_min"],
        "proteins": np.asarray(trace["states"][:, [1, 3, 5]], dtype=float),
    }


def aggregate(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    keys = ["knob_setting_id", "condition", "tau_step_min"]
    for group_keys, group in frame.groupby(keys, sort=True):
        row = dict(zip(keys, group_keys))
        row["n_trajectories"] = int(len(group))
        for column in (
            "mean_C31_flux_uM_h",
            "C31_mRNA_balance_factor",
            "mean_PLtetO1_activity",
            "Jensen_factor_Eh_over_hEfree",
            "phase_period_mean_generations",
            "phase_period_within_trajectory_cv",
        ):
            values = pd.to_numeric(group[column], errors="coerce")
            finite = values[np.isfinite(values)]
            row[f"{column}_mean"] = float(finite.mean()) if len(finite) else math.nan
            row[f"{column}_sd"] = float(finite.std(ddof=1)) if len(finite) > 1 else 0.0
        row["resolved_fraction"] = float(group["resolved"].astype(bool).mean())
        row["theoretical_mRNA_discretization_factor"] = float(
            group["theoretical_mRNA_discretization_factor"].iloc[0]
        )
        rows.append(row)
    return pd.DataFrame(rows)


def add_ode_ratios(summary: pd.DataFrame, deterministic_rows: pd.DataFrame) -> pd.DataFrame:
    deterministic = deterministic_rows.set_index("knob_setting_id")
    result = summary.copy()
    flux_ratios = []
    corrected_flux_ratios = []
    period_ratios = []
    for _, row in result.iterrows():
        reference = deterministic.loc[row["knob_setting_id"]]
        flux_ratio = mean_audit.safe_ratio(
            row["mean_C31_flux_uM_h_mean"], reference["mean_C31_flux_uM_h"]
        )
        factor = row["theoretical_mRNA_discretization_factor"]
        flux_ratios.append(flux_ratio)
        corrected_flux_ratios.append(mean_audit.safe_ratio(flux_ratio, factor))
        period_ratios.append(
            mean_audit.safe_ratio(
                row["phase_period_mean_generations_mean"],
                reference["realised_period_generations"],
            )
        )
    result["mean_flux_ratio_vs_ODE"] = flux_ratios
    result["mean_flux_ratio_vs_ODE_after_theoretical_mRNA_correction"] = corrected_flux_ratios
    result["mean_period_ratio_vs_ODE"] = period_ratios
    return result


def relative_change(value_a: float, value_b: float) -> float:
    scale = max(abs(value_a), abs(value_b), 1e-12)
    return float(abs(value_a - value_b) / scale)


def convergence_decisions(summary: pd.DataFrame, settings: ConvergenceSettings):
    rows = []
    ordered_steps = sorted(summary["tau_step_min"].unique())
    smallest, next_smallest = ordered_steps[0], ordered_steps[1]
    for (setting_id, condition), group in summary.groupby(
        ["knob_setting_id", "condition"], sort=True
    ):
        lookup = group.set_index("tau_step_min")
        fine = lookup.loc[smallest]
        next_fine = lookup.loc[next_smallest]
        flux_change = relative_change(
            fine["mean_flux_ratio_vs_ODE"], next_fine["mean_flux_ratio_vs_ODE"]
        )
        period_change = relative_change(
            fine["phase_period_mean_generations_mean"],
            next_fine["phase_period_mean_generations_mean"],
        )
        cv_change = relative_change(
            fine["phase_period_within_trajectory_cv_mean"],
            next_fine["phase_period_within_trajectory_cv_mean"],
        )
        balance_error = abs(fine["C31_mRNA_balance_factor_mean"] - 1.0)
        pass_flag = bool(
            flux_change <= settings.maximum_small_step_relative_change
            and period_change <= settings.maximum_small_step_relative_change
            and balance_error <= settings.maximum_small_step_mrna_balance_error
            and fine["resolved_fraction"] >= settings.minimum_resolved_fraction
        )
        rows.append(
            {
                "knob_setting_id": setting_id,
                "condition": condition,
                "smallest_tau_step_min": smallest,
                "next_smallest_tau_step_min": next_smallest,
                "small_step_flux_relative_change": flux_change,
                "small_step_period_relative_change": period_change,
                "small_step_CV_relative_change_diagnostic": cv_change,
                "smallest_step_mRNA_balance_error": balance_error,
                "smallest_step_resolved_fraction": fine["resolved_fraction"],
                "provisional_convergence_pass": pass_flag,
            }
        )
    return pd.DataFrame(rows)


def make_figures(summary: pd.DataFrame, decisions: pd.DataFrame, figure_dir: Path):
    setting_ids = tuple(sorted(summary["knob_setting_id"].unique()))
    conditions = (
        "stochastic_fixed_plasmids",
        "stochastic_fluctuating_plasmids",
    )
    labels = {conditions[0]: "fixed plasmids", conditions[1]: "fluctuating plasmids"}
    colors = {conditions[0]: "tab:blue", conditions[1]: "tab:orange"}

    fig, axes = plt.subplots(2, len(setting_ids), figsize=(6 * len(setting_ids), 8), squeeze=False)
    for column, setting_id in enumerate(setting_ids):
        local = summary[summary["knob_setting_id"].eq(setting_id)]
        for condition in conditions:
            curve = local[local["condition"].eq(condition)].sort_values("tau_step_min")
            axes[0, column].plot(
                curve["tau_step_min"],
                curve["mean_flux_ratio_vs_ODE"],
                marker="o",
                color=colors[condition],
                label=labels[condition],
            )
            axes[0, column].plot(
                curve["tau_step_min"],
                curve["theoretical_mRNA_discretization_factor"],
                linestyle="--",
                color=colors[condition],
                alpha=0.5,
            )
            axes[1, column].plot(
                curve["tau_step_min"],
                curve["C31_mRNA_balance_factor_mean"],
                marker="s",
                color=colors[condition],
                label=labels[condition],
            )
        axes[0, column].axhline(1.0, color="black", linewidth=1)
        axes[1, column].axhline(1.0, color="black", linewidth=1)
        axes[0, column].set_title(f"{setting_id}: mean flux / ODE")
        axes[1, column].set_title(f"{setting_id}: C31 mRNA balance")
        axes[0, column].set_ylabel("ratio")
        axes[1, column].set_ylabel("ratio")
        axes[1, column].set_xlabel("tau step (min)")
        axes[0, column].invert_xaxis()
        axes[1, column].invert_xaxis()
        axes[0, column].legend()
        axes[1, column].legend()
    fig.suptitle("Tau-leap mean convergence; dashed = theoretical finite-step factor")
    fig.tight_layout()
    fig.savefig(figure_dir / "Tau_step_mean_convergence.png", dpi=200)
    plt.close(fig)

    fig, axes = plt.subplots(2, len(setting_ids), figsize=(6 * len(setting_ids), 8), squeeze=False)
    for column, setting_id in enumerate(setting_ids):
        local = summary[summary["knob_setting_id"].eq(setting_id)]
        for condition in conditions:
            curve = local[local["condition"].eq(condition)].sort_values("tau_step_min")
            axes[0, column].errorbar(
                curve["tau_step_min"],
                curve["phase_period_mean_generations_mean"],
                yerr=curve["phase_period_mean_generations_sd"],
                marker="o",
                capsize=3,
                color=colors[condition],
                label=labels[condition],
            )
            axes[1, column].plot(
                curve["tau_step_min"],
                curve["phase_period_within_trajectory_cv_mean"],
                marker="s",
                color=colors[condition],
                label=labels[condition],
            )
        axes[0, column].set_title(f"{setting_id}: period")
        axes[1, column].set_title(f"{setting_id}: within-trajectory CV")
        axes[0, column].set_ylabel("generations")
        axes[1, column].set_ylabel("CV")
        axes[1, column].set_xlabel("tau step (min)")
        axes[0, column].invert_xaxis()
        axes[1, column].invert_xaxis()
        axes[0, column].legend()
        axes[1, column].legend()
    fig.suptitle("Tau-leap clock convergence")
    fig.tight_layout()
    fig.savefig(figure_dir / "Tau_step_period_CV_convergence.png", dpi=200)
    plt.close(fig)


def main():
    args = parse_args()
    base_settings = MODE_SETTINGS[args.mode]
    tau_steps = parse_steps(args.tau_steps_min, base_settings.tau_steps_min)
    settings = replace(base_settings, tau_steps_min=tau_steps)

    library_path, library = library_tools.load_k4_k5_library(args.library_csv)
    library = library_tools.validate_k4_k5_library(library)
    robust_source = deterministic_a.load_representatives(library_path)
    robust_lookup = robust_source.set_index("setting_id")

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = (
        args.output_dir.expanduser().resolve()
        if args.output_dir is not None
        else OUTPUT_ROOT / f"Tau_Step_Convergence_{args.mode}_{stamp}"
    )
    figure_dir = FIGURE_ROOT / output_dir.name
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    print(f"{VERSION}: mode={args.mode}")
    print(f"Tau steps (min): {tau_steps}")
    print(f"Robust library: {library_path}")

    deterministic_rows = []
    trajectory_rows = []
    parity_rows = []

    for setting_index, (_, library_row) in enumerate(library.iterrows()):
        setting_id = str(library_row["knob_setting_id"])
        robust_row = robust_lookup.loc[str(library_row["source_setting_id"])]
        deterministic_trace = mean_audit.deterministic_detailed(robust_row, mean_audit.AuditSettings(
            settings.trajectories_per_setting,
            settings.duration_generations,
            settings.transient_generations,
            tau_steps[0],
            settings.sample_step_min,
            settings.random_seed,
        ))
        deterministic_summary = mean_audit.summarize_trace(
            setting_id,
            "deterministic_reference",
            -1,
            deterministic_trace,
            settings.transient_generations,
        )
        deterministic_summary["realised_period_generations"] = float(
            library_row["realised_period_generations"]
        )
        deterministic_rows.append(deterministic_summary)
        stochastic_parameters = stochastic_a.build_parameters(library_row)
        parity_rows.extend(
            mean_audit.parameter_parity_rows(
                setting_id, deterministic_trace["parameters"], stochastic_parameters
            )
        )

        for tau_index, tau_step in enumerate(tau_steps):
            run_settings = make_run_settings(settings, tau_step)
            for fixed_plasmids, condition in (
                (True, "stochastic_fixed_plasmids"),
                (False, "stochastic_fluctuating_plasmids"),
            ):
                for trajectory_index in range(settings.trajectories_per_setting):
                    # Common random-number labels are retained across tau steps.
                    # Exact event streams differ because step counts differ.
                    seed = (
                        settings.random_seed
                        + setting_index * 100_000
                        + (0 if fixed_plasmids else 50_000)
                        + trajectory_index
                    )
                    print(
                        f"[{setting_id} | tau={tau_step:g} | {condition} | "
                        f"{trajectory_index + 1:02d}/{settings.trajectories_per_setting:02d}]"
                    )
                    trace = mean_audit.simulate_stochastic_detailed(
                        library_row, run_settings, seed, fixed_plasmids
                    )
                    row = mean_audit.summarize_trace(
                        setting_id,
                        condition,
                        trajectory_index,
                        trace,
                        settings.transient_generations,
                    )
                    analysis = stochastic_a.analyse_trajectory(
                        library_row, stochastic_analysis_payload(trace), run_settings
                    )
                    loss_rate = trace["parameters"].c31_mrna_intrinsic_loss_per_min + trace[
                        "parameters"
                    ].mu
                    row.update(
                        {
                            "tau_step_min": tau_step,
                            "resolved": bool(analysis["resolved"]),
                            "phase_period_mean_generations": analysis[
                                "phase_period_mean_generations"
                            ],
                            "phase_period_within_trajectory_cv": analysis[
                                "phase_period_within_trajectory_cv"
                            ],
                            "n_phase_cycles": analysis["n_phase_cycles"],
                            "theoretical_mRNA_discretization_factor": (
                                theoretical_birth_death_factor(loss_rate, tau_step)
                            ),
                        }
                    )
                    trajectory_rows.append(row)

    parity = pd.DataFrame(parity_rows)
    if not parity["parity_pass"].all():
        mismatch = parity.loc[~parity["parity_pass"]]
        mismatch.to_csv(output_dir / "parameter_mapping_mismatches.csv", index=False)
        raise RuntimeError(
            "Deterministic/stochastic mappings differ; see parameter_mapping_mismatches.csv"
        )

    trajectory_frame = pd.DataFrame(trajectory_rows)
    deterministic_frame = pd.DataFrame(deterministic_rows)
    summary = aggregate(trajectory_frame)
    summary = add_ode_ratios(summary, deterministic_frame)
    decisions = convergence_decisions(summary, settings)

    trajectory_frame.to_csv(output_dir / "tau_step_by_trajectory.csv", index=False)
    deterministic_frame.to_csv(output_dir / "deterministic_references.csv", index=False)
    summary.to_csv(output_dir / "tau_step_condition_summary.csv", index=False)
    decisions.to_csv(output_dir / "tau_step_convergence_decisions.csv", index=False)
    parity.to_csv(output_dir / "deterministic_stochastic_parameter_parity.csv", index=False)

    metadata = {
        "version": VERSION,
        "mode": args.mode,
        "settings": asdict(settings),
        "library_csv": str(library_path),
        "scope_boundary": (
            "Numerical convergence audit only; no biological parameter fitting and no B certification."
        ),
        "decision_rule": (
            "At the two smallest steps, flux and period change <=2%, smallest-step mRNA balance "
            "error <=1.5%, and resolved fraction >=0.8. CV change is reported diagnostically."
        ),
    }
    (output_dir / "scope_and_settings.json").write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    make_figures(summary, decisions, figure_dir)

    print(f"\nCompleted: {VERSION}")
    print(f"Raw tables: {output_dir}")
    print(f"Figures: {figure_dir}")
    print("\nTime-step convergence decisions:")
    for _, row in decisions.iterrows():
        print(
            f"{row['knob_setting_id']:>9s} | {row['condition']:<33s} | "
            f"flux change={row['small_step_flux_relative_change']:.3f} | "
            f"period change={row['small_step_period_relative_change']:.3f} | "
            f"mRNA error={row['smallest_step_mRNA_balance_error']:.3f} | "
            f"pass={bool(row['provisional_convergence_pass'])}"
        )
    if args.mode == "smoke":
        print(
            "\nSmoke decisions are workflow diagnostics only. Use audit mode before choosing "
            "a production tau step or rerunning stochastic A-to-B certification."
        )
    else:
        print(
            "\nAudit decisions use strict point-estimate tolerances. Inspect trajectory "
            "variability and convergence trends before interpreting a false flag as "
            "numerical non-convergence."
        )


if __name__ == "__main__":
    main()
