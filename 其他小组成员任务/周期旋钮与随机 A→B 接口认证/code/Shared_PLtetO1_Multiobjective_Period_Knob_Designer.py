"""Build a calibrated multi-objective period-knob lookup from existing ODE results.

This script DOES NOT run the ODE or SSA.  It post-processes the completed
synchronous-RBS x oscillator-mRNA-lifetime design grid and turns it into a
reproducible period-knob map.

The output is a design library, not a fitted final construct.  For every
requested target period it reports:

* period error and clock-quality constraints;
* C31 interface metrics;
* experimental-change cost relative to the frozen reference;
* local robustness on the discrete design grid;
* RBS, mRNA-lifetime and nonlinear-interaction contributions to period change;
* optional B-module and stochastic-CV validation when matching CSV files are
  supplied.

Examples (run from the project root):

    python ./code/Shared_PLtetO1_Multiobjective_Period_Knob_Designer.py
    python ./code/Shared_PLtetO1_Multiobjective_Period_Knob_Designer.py ^
        --target-periods 8,9,10,11,12,13 --period-tolerance 0.35

Optional validation tables must contain the two design-coordinate columns
``design_RBS_scale`` and ``design_oscillator_mRNA_half_life_min``.  The B table
may additionally contain ``B_pass``; the stochastic table may contain
``stochastic_period_cv`` and ``stochastic_mean_period_generations``.
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


VERSION = "Shared_PLtetO1_Multiobjective_Period_Knob_Designer_20260819"
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
OUTPUT_ROOT = PROJECT_ROOT / "outputs"
FIGURE_ROOT = PROJECT_ROOT / "figures"

RBS_COLUMN_CANDIDATES = ("design_RBS_scale", "oscillator_RBS_scale")
MRNA_COLUMN_CANDIDATES = (
    "design_oscillator_mRNA_half_life_min",
    "oscillator_mRNA_total_half_life_min",
)

BASELINE_RBS = 1.0
BASELINE_MRNA_HALF_LIFE_MIN = 2.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-csv",
        type=Path,
        default=None,
        help="Completed RBS x mRNA grid CSV. By default the newest matching output is used.",
    )
    parser.add_argument(
        "--target-periods",
        default="8,9,10,11,12,13",
        help="Comma-separated target periods in generations.",
    )
    parser.add_argument("--period-tolerance", type=float, default=0.35)
    parser.add_argument("--minimum-relative-amplitude", type=float, default=0.90)
    parser.add_argument("--maximum-residual-duty", type=float, default=0.55)
    parser.add_argument("--top-per-target", type=int, default=5)
    parser.add_argument(
        "--b-validation-csv",
        type=Path,
        default=None,
        help="Optional B validation table keyed by the two design coordinates.",
    )
    parser.add_argument(
        "--stochastic-validation-csv",
        type=Path,
        default=None,
        help="Optional SSA table keyed by the two design coordinates.",
    )
    parser.add_argument(
        "--maximum-stochastic-cv",
        type=float,
        default=0.20,
        help="Applied only when stochastic_period_cv is supplied.",
    )
    return parser.parse_args()


def parse_targets(text: str) -> list[float]:
    values = sorted({float(item.strip()) for item in text.split(",") if item.strip()})
    if not values:
        raise ValueError("At least one target period is required.")
    if any(value <= 0 for value in values):
        raise ValueError("Target periods must be positive.")
    return values


def locate_input_csv(explicit: Path | None) -> Path:
    if explicit is not None:
        path = explicit.expanduser().resolve()
        if not path.exists():
            raise FileNotFoundError(path)
        return path

    matches = sorted(
        OUTPUT_ROOT.glob("Shared_PLtetO1_RBS_mRNA_2D_*/*RBS_mRNA_2D_all_conditions.csv"),
        key=lambda path: path.stat().st_mtime,
        reverse=True,
    )
    if not matches:
        raise FileNotFoundError(
            "No completed RBS-mRNA grid was found under the project outputs folder. "
            "Pass --input-csv explicitly."
        )
    return matches[0]


def find_column(frame: pd.DataFrame, candidates: tuple[str, ...], label: str) -> str:
    for column in candidates:
        if column in frame.columns:
            return column
    raise KeyError(f"Could not find {label}; expected one of {candidates}")


def normalise_design_columns(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.copy()
    rbs_column = find_column(frame, RBS_COLUMN_CANDIDATES, "RBS coordinate")
    mrna_column = find_column(frame, MRNA_COLUMN_CANDIDATES, "mRNA coordinate")
    frame["design_RBS_scale"] = pd.to_numeric(frame[rbs_column], errors="coerce")
    frame["design_oscillator_mRNA_half_life_min"] = pd.to_numeric(
        frame[mrna_column], errors="coerce"
    )
    return frame


def optional_merge(frame: pd.DataFrame, path: Path | None, prefix: str) -> pd.DataFrame:
    if path is None:
        return frame
    validation = normalise_design_columns(pd.read_csv(path.expanduser().resolve()))
    keys = ["design_RBS_scale", "design_oscillator_mRNA_half_life_min"]
    keep = keys + [column for column in validation.columns if column not in frame.columns and column not in keys]
    validation = validation[keep].drop_duplicates(keys)
    renamed = {
        column: f"{prefix}{column}"
        for column in validation.columns
        if column not in keys and column in frame.columns
    }
    return frame.merge(validation.rename(columns=renamed), on=keys, how="left", validate="many_to_one")


def nearest_value(values: np.ndarray, target: float) -> float:
    return float(values[np.argmin(np.abs(values - target))])


def period_at(grid: pd.DataFrame, rbs: float, mrna: float) -> float:
    rows = grid.loc[
        np.isclose(grid["design_RBS_scale"], rbs)
        & np.isclose(grid["design_oscillator_mRNA_half_life_min"], mrna),
        "period_generations",
    ]
    if rows.empty:
        return float("nan")
    return float(rows.iloc[0])


def add_period_contributions(frame: pd.DataFrame) -> pd.DataFrame:
    """Decompose the discrete period shift into two main effects + interaction."""
    frame = frame.copy()
    rbs_values = np.sort(frame["design_RBS_scale"].dropna().unique())
    mrna_values = np.sort(frame["design_oscillator_mRNA_half_life_min"].dropna().unique())
    reference_rbs = nearest_value(rbs_values, BASELINE_RBS)
    reference_mrna = nearest_value(mrna_values, BASELINE_MRNA_HALF_LIFE_MIN)
    baseline_period = period_at(frame, reference_rbs, reference_mrna)
    if not np.isfinite(baseline_period):
        raise ValueError("The RBS=1, mRNA half-life=2 min reference is missing or unresolved.")

    rows: list[dict[str, float]] = []
    for _, row in frame.iterrows():
        rbs = float(row["design_RBS_scale"])
        mrna = float(row["design_oscillator_mRNA_half_life_min"])
        period = float(row["period_generations"])
        rbs_main = period_at(frame, rbs, reference_mrna) - baseline_period
        mrna_main = period_at(frame, reference_rbs, mrna) - baseline_period
        interaction = period - baseline_period - rbs_main - mrna_main
        denominator = abs(rbs_main) + abs(mrna_main) + abs(interaction)
        rows.append(
            {
                "baseline_period_generations": baseline_period,
                "period_change_from_baseline_generations": period - baseline_period,
                "period_change_RBS_main_generations": rbs_main,
                "period_change_mRNA_main_generations": mrna_main,
                "period_change_interaction_generations": interaction,
                "absolute_contribution_RBS_fraction": abs(rbs_main) / denominator if denominator else 0.0,
                "absolute_contribution_mRNA_fraction": abs(mrna_main) / denominator if denominator else 0.0,
                "absolute_contribution_interaction_fraction": abs(interaction) / denominator if denominator else 0.0,
            }
        )
    return pd.concat([frame.reset_index(drop=True), pd.DataFrame(rows)], axis=1)


def add_local_robustness(frame: pd.DataFrame) -> pd.DataFrame:
    """Estimate local period sensitivity from immediate neighbours on the grid."""
    frame = frame.copy()
    rbs_values = np.sort(frame["design_RBS_scale"].dropna().unique())
    mrna_values = np.sort(frame["design_oscillator_mRNA_half_life_min"].dropna().unique())
    rbs_index = {value: index for index, value in enumerate(rbs_values)}
    mrna_index = {value: index for index, value in enumerate(mrna_values)}

    def slope(values: np.ndarray, index: int, fixed: float, coordinate: str) -> float:
        left = max(index - 1, 0)
        right = min(index + 1, len(values) - 1)
        if left == right:
            return float("nan")
        if coordinate == "rbs":
            p_left = period_at(frame, float(values[left]), fixed)
            p_right = period_at(frame, float(values[right]), fixed)
        else:
            p_left = period_at(frame, fixed, float(values[left]))
            p_right = period_at(frame, fixed, float(values[right]))
        return (p_right - p_left) / float(values[right] - values[left])

    rbs_slopes = []
    mrna_slopes = []
    for _, row in frame.iterrows():
        rbs = float(row["design_RBS_scale"])
        mrna = float(row["design_oscillator_mRNA_half_life_min"])
        rbs_slopes.append(slope(rbs_values, rbs_index[rbs], mrna, "rbs"))
        mrna_slopes.append(slope(mrna_values, mrna_index[mrna], rbs, "mrna"))
    frame["local_period_slope_per_RBS_scale"] = rbs_slopes
    frame["local_period_slope_per_mRNA_min"] = mrna_slopes
    frame["local_sensitivity_norm"] = np.sqrt(np.square(rbs_slopes) + np.square(mrna_slopes))
    return frame


def scale_zero_one(series: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(series, errors="coerce")
    low = numeric.min(skipna=True)
    high = numeric.max(skipna=True)
    if not np.isfinite(low) or not np.isfinite(high) or math.isclose(low, high):
        return pd.Series(0.0, index=series.index)
    return (numeric - low) / (high - low)


def pareto_mask(values: np.ndarray) -> np.ndarray:
    """Return non-dominated rows for minimisation objectives."""
    count = values.shape[0]
    efficient = np.ones(count, dtype=bool)
    for index in range(count):
        if not efficient[index]:
            continue
        dominated_by_other = np.any(
            np.all(values <= values[index], axis=1)
            & np.any(values < values[index], axis=1)
        )
        if dominated_by_other:
            efficient[index] = False
    return efficient


def build_lookup(
    frame: pd.DataFrame,
    targets: list[float],
    period_tolerance: float,
    minimum_amplitude: float,
    maximum_residual_duty: float,
    maximum_stochastic_cv: float,
    top_per_target: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    frame = frame.copy()
    resolved = frame["status"].eq("resolved_oscillation") & frame["resolved_oscillation"].astype(bool)
    amplitude_ok = pd.to_numeric(frame["TetR_relative_amplitude"], errors="coerce") >= minimum_amplitude
    duty_column = "C31_flux_residual_duty_above_0p5_uM_h_when_resolved"
    if duty_column not in frame:
        duty_column = "C31_flux_residual_duty_above_0p5_uM_h"
    interface_shape_ok = pd.to_numeric(frame[duty_column], errors="coerce") <= maximum_residual_duty
    frame["deterministic_clock_constraint_pass"] = resolved & amplitude_ok
    frame["A_interface_shape_constraint_pass"] = interface_shape_ok

    if "B_pass" in frame:
        b_known = frame["B_pass"].notna()
        b_pass = frame["B_pass"].fillna(False).astype(bool)
    else:
        b_known = pd.Series(False, index=frame.index)
        b_pass = pd.Series(False, index=frame.index)
    frame["B_validation_status"] = np.where(b_known, np.where(b_pass, "pass", "fail"), "pending")

    if "stochastic_period_cv" in frame:
        stochastic_cv = pd.to_numeric(frame["stochastic_period_cv"], errors="coerce")
        ssa_known = stochastic_cv.notna()
        ssa_pass = stochastic_cv <= maximum_stochastic_cv
    else:
        stochastic_cv = pd.Series(np.nan, index=frame.index)
        ssa_known = pd.Series(False, index=frame.index)
        ssa_pass = pd.Series(False, index=frame.index)
        frame["stochastic_period_cv"] = np.nan
    frame["stochastic_validation_status"] = np.where(
        ssa_known, np.where(ssa_pass, "pass", "fail"), "pending"
    )

    # Experimental-change cost is deliberately simple and transparent: it
    # penalises fold changes on a log2 scale rather than declaring a construct
    # biologically optimal from simulation alone.
    frame["experimental_change_cost"] = (
        np.abs(np.log2(frame["design_RBS_scale"] / BASELINE_RBS))
        + np.abs(
            np.log2(
                frame["design_oscillator_mRNA_half_life_min"]
                / BASELINE_MRNA_HALF_LIFE_MIN
            )
        )
    )
    frame["output_deviation_cost"] = (
        scale_zero_one(abs(frame["C31_flux_peak_uM_h"] - frame["C31_flux_peak_uM_h"].median()))
        + scale_zero_one(abs(frame["C31_flux_dose_uM_per_cycle"] - frame["C31_flux_dose_uM_per_cycle"].median()))
    ) / 2.0
    frame["robustness_cost"] = scale_zero_one(frame["local_sensitivity_norm"])

    lookup_rows = []
    for target in targets:
        candidates = frame.copy()
        candidates["target_period_generations"] = target
        candidates["period_error_generations"] = abs(candidates["period_generations"] - target)
        candidates["within_period_tolerance"] = candidates["period_error_generations"] <= period_tolerance
        candidates["all_available_constraints_pass"] = (
            candidates["deterministic_clock_constraint_pass"]
            & candidates["A_interface_shape_constraint_pass"]
            & (~b_known | b_pass)
            & (~ssa_known | ssa_pass)
        )

        # The score is a transparent ranking aid, never a fitted biological
        # truth. Period accuracy dominates, followed by robustness, construct
        # effort and preservation of the current C31 interface.
        candidates["provisional_design_score"] = (
            candidates["period_error_generations"] / max(period_tolerance, 1e-9)
            + 0.30 * candidates["robustness_cost"]
            + 0.20 * scale_zero_one(candidates["experimental_change_cost"])
            + 0.20 * candidates["output_deviation_cost"]
        )
        objectives = candidates[
            [
                "period_error_generations",
                "robustness_cost",
                "experimental_change_cost",
                "output_deviation_cost",
            ]
        ].to_numpy(float)
        candidates["pareto_efficient"] = pareto_mask(objectives)
        candidates = candidates.sort_values(
            [
                "all_available_constraints_pass",
                "within_period_tolerance",
                "provisional_design_score",
            ],
            ascending=[False, False, True],
        )
        candidates["rank_for_target"] = np.arange(1, len(candidates) + 1)
        lookup_rows.append(candidates.head(top_per_target))

    lookup = pd.concat(lookup_rows, ignore_index=True)
    representatives = lookup.loc[lookup["rank_for_target"].eq(1)].copy()
    representatives["selection_label"] = np.where(
        representatives["within_period_tolerance"]
        & representatives["all_available_constraints_pass"],
        "provisional_representative_pending_B_and_SSA",
        "nearest_available_design_not_certified",
    )
    return frame, lookup, representatives


def make_figures(frame: pd.DataFrame, representatives: pd.DataFrame, figure_dir: Path) -> None:
    figure_dir.mkdir(parents=True, exist_ok=True)
    rbs_values = np.sort(frame["design_RBS_scale"].unique())
    mrna_values = np.sort(frame["design_oscillator_mRNA_half_life_min"].unique())
    period_matrix = frame.pivot(
        index="design_oscillator_mRNA_half_life_min",
        columns="design_RBS_scale",
        values="period_generations",
    ).reindex(index=mrna_values, columns=rbs_values)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.4), constrained_layout=True)
    image = axes[0].imshow(period_matrix.to_numpy(float), origin="lower", aspect="auto", cmap="viridis")
    axes[0].set_xticks(range(len(rbs_values)), [f"{value:g}" for value in rbs_values])
    axes[0].set_yticks(range(len(mrna_values)), [f"{value:g}" for value in mrna_values])
    axes[0].set_xlabel("synchronous three-arm RBS scale")
    axes[0].set_ylabel("oscillator mRNA half-life (min)")
    axes[0].set_title("Calibrated genetic period-knob map")
    for row in range(period_matrix.shape[0]):
        for column in range(period_matrix.shape[1]):
            value = period_matrix.iat[row, column]
            axes[0].text(column, row, f"{value:.2f}", ha="center", va="center", fontsize=8,
                         color="white" if value < np.nanmedian(period_matrix.to_numpy(float)) else "black")
    colorbar = fig.colorbar(image, ax=axes[0], shrink=0.88)
    colorbar.set_label("period (generations)")

    axes[1].plot(
        representatives["target_period_generations"],
        representatives["period_generations"],
        "o-",
        label="provisional representative",
    )
    limits = [
        min(representatives["target_period_generations"].min(), representatives["period_generations"].min()),
        max(representatives["target_period_generations"].max(), representatives["period_generations"].max()),
    ]
    axes[1].plot(limits, limits, "--", color="0.4", label="ideal target = realised")
    for _, row in representatives.iterrows():
        axes[1].annotate(
            f"RBS={row['design_RBS_scale']:g}\nt1/2={row['design_oscillator_mRNA_half_life_min']:g} min",
            (row["target_period_generations"], row["period_generations"]),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=7,
        )
    axes[1].set_xlabel("requested period (generations)")
    axes[1].set_ylabel("model-realised period (generations)")
    axes[1].set_title("Lookup behaviour (not final construct certification)")
    axes[1].legend(fontsize=8)
    fig.savefig(figure_dir / "Period_Knob_calibration_and_lookup.png", dpi=200)
    plt.close(fig)

    contribution_columns = [
        "absolute_contribution_RBS_fraction",
        "absolute_contribution_mRNA_fraction",
        "absolute_contribution_interaction_fraction",
    ]
    labels = ["RBS main effect", "mRNA-lifetime main effect", "interaction"]
    colors = ["#4472C4", "#ED7D31", "#70AD47"]
    fig, ax = plt.subplots(figsize=(10, 5.4), constrained_layout=True)
    bottom = np.zeros(len(representatives))
    x = np.arange(len(representatives))
    for column, label, color in zip(contribution_columns, labels, colors):
        values = representatives[column].to_numpy(float)
        ax.bar(x, values, bottom=bottom, label=label, color=color)
        bottom += values
    ax.set_xticks(x, [f"target {value:g}" for value in representatives["target_period_generations"]])
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("absolute contribution fraction")
    ax.set_title("How each representative setting changes the period")
    ax.legend()
    fig.savefig(figure_dir / "Period_Knob_effect_contributions.png", dpi=200)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    targets = parse_targets(args.target_periods)
    input_csv = locate_input_csv(args.input_csv)
    frame = normalise_design_columns(pd.read_csv(input_csv))
    frame = optional_merge(frame, args.b_validation_csv, "B__")
    frame = optional_merge(frame, args.stochastic_validation_csv, "SSA__")
    frame = add_period_contributions(frame)
    frame = add_local_robustness(frame)
    frame, lookup, representatives = build_lookup(
        frame=frame,
        targets=targets,
        period_tolerance=args.period_tolerance,
        minimum_amplitude=args.minimum_relative_amplitude,
        maximum_residual_duty=args.maximum_residual_duty,
        maximum_stochastic_cv=args.maximum_stochastic_cv,
        top_per_target=args.top_per_target,
    )

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = OUTPUT_ROOT / f"Multiobjective_Period_Knob_{stamp}"
    figure_dir = FIGURE_ROOT / f"Multiobjective_Period_Knob_{stamp}"
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    frame.to_csv(output_dir / "period_knob_all_designs.csv", index=False)
    lookup.to_csv(output_dir / "period_knob_ranked_lookup.csv", index=False)
    representatives.to_csv(output_dir / "period_knob_representative_library.csv", index=False)
    make_figures(frame, representatives, figure_dir)

    scope = {
        "version": VERSION,
        "input_csv": str(input_csv),
        "target_periods_generations": targets,
        "period_tolerance_generations": args.period_tolerance,
        "minimum_relative_amplitude": args.minimum_relative_amplitude,
        "maximum_residual_duty_above_0p5_uM_h": args.maximum_residual_duty,
        "maximum_stochastic_period_cv_if_supplied": args.maximum_stochastic_cv,
        "B_validation_supplied": args.b_validation_csv is not None,
        "stochastic_validation_supplied": args.stochastic_validation_csv is not None,
        "baseline_for_effect_decomposition": {
            "RBS_scale": BASELINE_RBS,
            "oscillator_mRNA_half_life_min": BASELINE_MRNA_HALF_LIFE_MIN,
        },
        "interpretation": [
            "This is a calibrated period-knob library, not optimisation to one preferred period.",
            "Representative rows are provisional until B and stochastic validation are supplied.",
            "Deterministic period CV is not used as biological noise.",
            "The contribution fractions are design-specific; no universal RBS:mRNA ratio is assumed.",
        ],
    }
    (output_dir / "scope_and_constraints.json").write_text(
        json.dumps(scope, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    print(f"Completed: {VERSION}")
    print(f"Input grid: {input_csv}")
    print(f"Raw tables: {output_dir}")
    print(f"Figures: {figure_dir}")
    print("\nProvisional lookup (B/SSA remain pending unless validation CSVs were supplied):")
    for _, row in representatives.iterrows():
        print(
            f"target={row['target_period_generations']:.2f} gen -> "
            f"realised={row['period_generations']:.3f} gen | "
            f"RBS={row['design_RBS_scale']:.2f}, "
            f"mRNA t1/2={row['design_oscillator_mRNA_half_life_min']:.2f} min | "
            f"B={row['B_validation_status']}, SSA={row['stochastic_validation_status']}"
        )


if __name__ == "__main__":
    main()
