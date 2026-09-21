"""Post-process period-knob robustness runs to identify failure drivers.

This program DOES NOT solve the ODE and does not generate new parameter
scenarios. It reads the completed confirmation tables and separates:

1. common-mode period shifts shared by all knob settings;
2. setting-specific target-period errors;
3. adjacent-setting ordering failures;
4. loss of oscillation or amplitude.

The analysis is diagnostic rather than causal. Spearman rank associations and
high-versus-low quartile failure contrasts are reported together so that a
large effect is not inferred from one statistic alone.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
OUTPUT_ROOT = PROJECT_ROOT / "outputs"
FIGURE_ROOT = PROJECT_ROOT / "figures"

PARAMETERS = (
    "Td_min",
    "oscillator_plasmid_copies",
    "K_copies",
    "hill_n",
    "PLlacO1_leak",
    "PLtetO1_leak",
    "PR_leak",
    "oscillator_RBS_multiplier",
    "oscillator_mrna_half_life_multiplier",
    "effective_sponge_sites",
    "C31_RBS_multiplier",
    "C31_mrna_half_life_multiplier",
    "peak_load_fraction",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Attribute period-knob robustness and ordering failures."
    )
    parser.add_argument(
        "--input-dir", type=Path, default=None,
        help="Completed Period_Knob_Robustness confirmation directory.",
    )
    return parser.parse_args()


def locate_latest_input() -> Path:
    candidates = sorted(
        OUTPUT_ROOT.glob("Period_Knob_Robustness_confirmation_*"),
        key=lambda path: path.stat().st_mtime,
        reverse=True,
    )
    for candidate in candidates:
        required = (
            candidate / "period_knob_robustness_run_metrics.csv",
            candidate / "period_knob_robustness_ordering_audit.csv",
            candidate / "uncertainty_scenarios.csv",
        )
        if all(path.exists() for path in required):
            return candidate
    raise FileNotFoundError(
        "No completed Period_Knob_Robustness_confirmation directory was found."
    )


def safe_spearman(x: pd.Series, y: pd.Series) -> tuple[float, float, int]:
    frame = pd.DataFrame({"x": x, "y": y}).replace([np.inf, -np.inf], np.nan).dropna()
    if len(frame) < 5 or frame["x"].nunique() < 2 or frame["y"].nunique() < 2:
        return np.nan, np.nan, len(frame)
    result = spearmanr(frame["x"], frame["y"])
    return float(result.statistic), float(result.pvalue), len(frame)


def quartile_contrast(x: pd.Series, outcome: pd.Series) -> tuple[float, float, float, int, int]:
    frame = pd.DataFrame({"x": x, "outcome": outcome}).dropna()
    if len(frame) < 8 or frame["x"].nunique() < 4:
        return np.nan, np.nan, np.nan, 0, 0
    low_cut, high_cut = frame["x"].quantile([0.25, 0.75])
    low = frame.loc[frame["x"] <= low_cut, "outcome"]
    high = frame.loc[frame["x"] >= high_cut, "outcome"]
    low_mean = float(low.mean()) if len(low) else np.nan
    high_mean = float(high.mean()) if len(high) else np.nan
    return low_mean, high_mean, high_mean - low_mean, len(low), len(high)


def build_scenario_library_metrics(run_metrics: pd.DataFrame) -> pd.DataFrame:
    periods = run_metrics.pivot(
        index="scenario_id", columns="setting_order", values="period_generations"
    ).sort_index(axis=1)
    targets = (
        run_metrics[["setting_order", "target_period_generations"]]
        .drop_duplicates()
        .set_index("setting_order")["target_period_generations"]
        .sort_index()
    )
    signed_errors = periods.subtract(targets, axis="columns")
    gaps = periods.diff(axis=1).iloc[:, 1:]
    result = pd.DataFrame(index=periods.index)
    result["library_mean_signed_period_shift"] = signed_errors.mean(axis=1)
    result["library_mean_absolute_period_error"] = signed_errors.abs().mean(axis=1)
    result["minimum_adjacent_period_gap"] = gaps.min(axis=1)
    result["mean_adjacent_period_gap"] = gaps.mean(axis=1)
    result["all_periods_resolved"] = (
        run_metrics.pivot(
            index="scenario_id", columns="setting_order", values="resolved_oscillation"
        ).min(axis=1) >= 0.5
    )
    return result.reset_index()


def setting_attribution(run_metrics: pd.DataFrame) -> pd.DataFrame:
    rows = []
    data = run_metrics.copy()
    data["signed_period_error"] = (
        data["period_generations"] - data["target_period_generations"]
    )
    data["absolute_period_error"] = data["signed_period_error"].abs()
    data["period_failure"] = (~data["period_tolerance_pass"].astype(bool)).astype(float)
    data["amplitude_failure"] = (~data["amplitude_pass"].astype(bool)).astype(float)
    data["unresolved_failure"] = (data["resolved_oscillation"] < 0.5).astype(float)

    for setting_id, subset in data.groupby("setting_id", sort=False):
        for parameter in PARAMETERS:
            signed_rho, signed_p, n_signed = safe_spearman(
                subset[parameter], subset["signed_period_error"]
            )
            absolute_rho, absolute_p, n_absolute = safe_spearman(
                subset[parameter], subset["absolute_period_error"]
            )
            low_fail, high_fail, fail_delta, n_low, n_high = quartile_contrast(
                subset[parameter], subset["period_failure"]
            )
            _, _, unresolved_delta, _, _ = quartile_contrast(
                subset[parameter], subset["unresolved_failure"]
            )
            _, _, amplitude_delta, _, _ = quartile_contrast(
                subset[parameter], subset["amplitude_failure"]
            )
            rows.append(
                {
                    "setting_id": setting_id,
                    "target_period_generations": float(
                        subset["target_period_generations"].iloc[0]
                    ),
                    "parameter": parameter,
                    "signed_error_spearman_rho": signed_rho,
                    "signed_error_spearman_p": signed_p,
                    "absolute_error_spearman_rho": absolute_rho,
                    "absolute_error_spearman_p": absolute_p,
                    "period_failure_rate_low_quartile": low_fail,
                    "period_failure_rate_high_quartile": high_fail,
                    "period_failure_high_minus_low": fail_delta,
                    "unresolved_failure_high_minus_low": unresolved_delta,
                    "amplitude_failure_high_minus_low": amplitude_delta,
                    "n_spearman_signed": n_signed,
                    "n_spearman_absolute": n_absolute,
                    "n_low_quartile": n_low,
                    "n_high_quartile": n_high,
                }
            )
    return pd.DataFrame(rows)


def ordering_attribution(
    scenarios: pd.DataFrame,
    ordering: pd.DataFrame,
    library_metrics: pd.DataFrame,
) -> pd.DataFrame:
    data = scenarios.merge(ordering, on="scenario_id", validate="one_to_one")
    data = data.merge(library_metrics, on="scenario_id", validate="one_to_one")
    data["ordering_failure"] = (~data["all_adjacent_ordered"].astype(bool)).astype(float)
    rows = []
    for parameter in PARAMETERS:
        failure_rho, failure_p, n_failure = safe_spearman(
            data[parameter], data["ordering_failure"]
        )
        gap_rho, gap_p, n_gap = safe_spearman(
            data[parameter], data["minimum_adjacent_period_gap"]
        )
        shift_rho, shift_p, n_shift = safe_spearman(
            data[parameter], data["library_mean_signed_period_shift"]
        )
        low_fail, high_fail, failure_delta, n_low, n_high = quartile_contrast(
            data[parameter], data["ordering_failure"]
        )
        rows.append(
            {
                "parameter": parameter,
                "ordering_failure_spearman_rho": failure_rho,
                "ordering_failure_spearman_p": failure_p,
                "minimum_gap_spearman_rho": gap_rho,
                "minimum_gap_spearman_p": gap_p,
                "common_period_shift_spearman_rho": shift_rho,
                "common_period_shift_spearman_p": shift_p,
                "ordering_failure_rate_low_quartile": low_fail,
                "ordering_failure_rate_high_quartile": high_fail,
                "ordering_failure_high_minus_low": failure_delta,
                "n_failure": n_failure,
                "n_gap": n_gap,
                "n_shift": n_shift,
                "n_low_quartile": n_low,
                "n_high_quartile": n_high,
            }
        )
    return pd.DataFrame(rows)


def make_figures(
    setting_results: pd.DataFrame,
    ordering_results: pd.DataFrame,
    figure_dir: Path,
) -> None:
    setting_order = (
        setting_results[["setting_id", "target_period_generations"]]
        .drop_duplicates()
        .sort_values("target_period_generations")["setting_id"]
        .tolist()
    )
    heat = setting_results.pivot(
        index="parameter", columns="setting_id", values="absolute_error_spearman_rho"
    ).reindex(index=PARAMETERS, columns=setting_order)
    fig, ax = plt.subplots(figsize=(8, 7))
    image = ax.imshow(heat.to_numpy(float), cmap="coolwarm", vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(range(len(heat.columns)), heat.columns, rotation=25, ha="right")
    ax.set_yticks(range(len(heat.index)), heat.index)
    ax.set_title("Parameter association with absolute target-period error")
    fig.colorbar(image, ax=ax, label="Spearman rho")
    fig.tight_layout()
    fig.savefig(figure_dir / "Period_Knob_absolute_error_attribution.png", dpi=180)
    plt.close(fig)

    ranked = ordering_results.assign(
        absolute_risk_contrast=lambda frame: frame[
            "ordering_failure_high_minus_low"
        ].abs()
    ).sort_values("absolute_risk_contrast", ascending=True)
    fig, ax = plt.subplots(figsize=(9, 6))
    colors = np.where(
        ranked["ordering_failure_high_minus_low"] >= 0, "tab:red", "tab:blue"
    )
    ax.barh(
        ranked["parameter"], ranked["ordering_failure_high_minus_low"], color=colors
    )
    ax.axvline(0, color="black", lw=0.8)
    ax.set_xlabel("ordering-failure rate: high quartile minus low quartile")
    ax.set_title("Which uncertainties disrupt adjacent knob ordering?")
    fig.tight_layout()
    fig.savefig(figure_dir / "Period_Knob_ordering_failure_drivers.png", dpi=180)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    input_dir = args.input_dir.resolve() if args.input_dir else locate_latest_input()
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = OUTPUT_ROOT / f"Period_Knob_Failure_Attribution_{timestamp}"
    figure_dir = FIGURE_ROOT / output_dir.name
    output_dir.mkdir(parents=True, exist_ok=False)
    figure_dir.mkdir(parents=True, exist_ok=False)

    run_metrics = pd.read_csv(input_dir / "period_knob_robustness_run_metrics.csv")
    ordering = pd.read_csv(input_dir / "period_knob_robustness_ordering_audit.csv")
    scenarios = pd.read_csv(input_dir / "uncertainty_scenarios.csv")

    library_metrics = build_scenario_library_metrics(run_metrics)
    setting_results = setting_attribution(run_metrics)
    ordering_results = ordering_attribution(scenarios, ordering, library_metrics)

    ranking = ordering_results.copy()
    ranking["ordering_driver_score"] = (
        ranking["ordering_failure_spearman_rho"].abs().fillna(0.0)
        + ranking["ordering_failure_high_minus_low"].abs().fillna(0.0)
        + ranking["minimum_gap_spearman_rho"].abs().fillna(0.0)
    )
    ranking = ranking.sort_values("ordering_driver_score", ascending=False)

    library_metrics.to_csv(output_dir / "scenario_level_library_metrics.csv", index=False)
    setting_results.to_csv(output_dir / "setting_parameter_attribution.csv", index=False)
    ordering_results.to_csv(output_dir / "ordering_parameter_attribution.csv", index=False)
    ranking.to_csv(output_dir / "ranked_ordering_failure_drivers.csv", index=False)
    (output_dir / "scope_and_provenance.json").write_text(
        json.dumps(
            {
                "version": "Period_Knob_Robustness_Failure_Attribution_20260821",
                "input_dir": str(input_dir),
                "methods": [
                    "Spearman rank association",
                    "high-versus-low quartile failure-rate contrast",
                    "scenario-level common shift and adjacent gap decomposition",
                ],
                "limitations": (
                    "Associations diagnose the sampled uncertainty space and do not "
                    "prove molecular causation. Parameter ranges are stress-test "
                    "ranges unless experimentally calibrated."
                ),
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )
    make_figures(setting_results, ordering_results, figure_dir)

    print("Completed: Period_Knob_Robustness_Failure_Attribution_20260821")
    print(f"Input: {input_dir}")
    print(f"Raw tables: {output_dir}")
    print(f"Figures: {figure_dir}")
    print("\nTop provisional ordering-failure drivers:")
    for row in ranking.head(8).itertuples(index=False):
        print(
            f"{row.parameter:<42} score={row.ordering_driver_score:.3f} | "
            f"failure rho={row.ordering_failure_spearman_rho:+.3f} | "
            f"high-low risk={row.ordering_failure_high_minus_low:+.3f} | "
            f"minimum-gap rho={row.minimum_gap_spearman_rho:+.3f}"
        )


if __name__ == "__main__":
    main()
