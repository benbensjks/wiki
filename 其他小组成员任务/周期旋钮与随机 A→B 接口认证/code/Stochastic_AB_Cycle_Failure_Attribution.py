r"""Post-process stochastic A-to-B interface failures without rerunning simulations.

Primary inference unit
----------------------
One stochastic A cycle.  The PB-start and LR-start B simulations driven by the
same A cycle are paired.  A cycle is successful only when both initial DNA
states complete exactly one correct transition.

The program compares successful and failed cycles using peak flux, FWHM, total
dose, post-peak dose and cycle duration.  Each flux feature is also normalized
to the median deterministic-reference value for the same knob setting.
Associations are descriptive, not causal and not a parameter fit.
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
from scipy.stats import mannwhitneyu


VERSION = "Stochastic_AB_Cycle_Failure_Attribution_20260824"
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
OUTPUT_ROOT = PROJECT_ROOT / "outputs"
FIGURE_ROOT = PROJECT_ROOT / "figures"

FEATURES = (
    "flux_peak_uM_h",
    "flux_FWHM_h",
    "flux_dose_uM",
    "post_peak_dose_uM",
    "cycle_duration_generations",
)
DISPLAY_LABELS = {
    "flux_peak_uM_h": "peak flux",
    "flux_FWHM_h": "FWHM",
    "flux_dose_uM": "dose/cycle",
    "post_peak_dose_uM": "post-peak dose",
    "cycle_duration_generations": "cycle duration",
}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=None,
        help="Stochastic_AB_Interface_Audit_audit_tau0.05min_* output directory.",
    )
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--bootstrap-replicates", type=int, default=2000)
    parser.add_argument("--random-seed", type=int, default=20260824)
    return parser.parse_args()


def newest_input_dir(explicit: Path | None) -> Path:
    if explicit is not None:
        path = explicit.expanduser().resolve()
        if not path.is_dir():
            raise FileNotFoundError(path)
        return path
    candidates = sorted(
        OUTPUT_ROOT.glob("Stochastic_AB_Interface_Audit_audit_tau0.05min_*"),
        key=lambda path: path.stat().st_mtime,
    )
    if not candidates:
        raise FileNotFoundError(
            "No Stochastic_AB_Interface_Audit_audit_tau0.05min_* directory found"
        )
    return candidates[-1]


def as_bool(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    return series.astype(str).str.lower().isin({"true", "1", "yes"})


def failure_category(row: pd.Series) -> str:
    if bool(row["cycle_correct"]):
        return "correct"
    crossings = int(row["total_crossings"])
    if crossings == 0:
        return "zero crossing"
    if crossings > 1:
        return "multiple crossings"
    if not bool(row["endpoint_complete"]):
        return "incomplete endpoint"
    return "wrong direction/phase"


def cluster_bootstrap_median_difference(
    frame: pd.DataFrame,
    feature: str,
    replicates: int,
    seed: int,
) -> tuple[float, float, float]:
    """Failure minus success median, bootstrapped by A trajectory."""
    clusters = frame["cluster_id"].drop_duplicates().to_numpy()
    rng = np.random.default_rng(seed)
    values = []
    for _ in range(replicates):
        sampled = rng.choice(clusters, size=len(clusters), replace=True)
        pieces = [frame.loc[frame["cluster_id"] == cluster] for cluster in sampled]
        boot = pd.concat(pieces, ignore_index=True)
        failed = boot.loc[~boot["both_initial_states_correct"], feature].dropna()
        passed = boot.loc[boot["both_initial_states_correct"], feature].dropna()
        if len(failed) and len(passed):
            values.append(float(failed.median() - passed.median()))
    observed_failed = frame.loc[~frame["both_initial_states_correct"], feature].dropna()
    observed_passed = frame.loc[frame["both_initial_states_correct"], feature].dropna()
    observed = float(observed_failed.median() - observed_passed.median())
    if not values:
        return observed, math.nan, math.nan
    return observed, float(np.quantile(values, 0.025)), float(np.quantile(values, 0.975))


def make_cycle_table(b_cycles: pd.DataFrame, flux: pd.DataFrame) -> pd.DataFrame:
    raw_b = b_cycles.loc[b_cycles["interface_mode"] == "stochastic_raw"].copy()
    raw_flux = flux.loc[flux["interface_mode"] == "stochastic_raw"].copy()
    raw_b["cycle_correct"] = as_bool(raw_b["cycle_correct"])
    raw_b["endpoint_complete"] = as_bool(raw_b["endpoint_complete"])
    raw_b["failure_category"] = raw_b.apply(failure_category, axis=1)

    keys = ["knob_setting_id", "trajectory_index", "interface_mode", "cycle_index"]
    feature_columns = keys + [
        "flux_peak_uM_h",
        "flux_FWHM_h",
        "flux_dose_uM",
        "post_peak_dose_uM",
    ]
    merged = raw_b.merge(raw_flux[feature_columns], on=keys, how="left", validate="many_to_one")

    outcome = merged.pivot_table(
        index=keys,
        columns="initial_DNA_state",
        values="cycle_correct",
        aggfunc="first",
    ).reset_index()
    if not {"PB", "LR"}.issubset(outcome.columns):
        raise RuntimeError("Both PB and LR initial-state results are required")
    outcome = outcome.rename(columns={"PB": "PB_correct", "LR": "LR_correct"})
    outcome["both_initial_states_correct"] = outcome["PB_correct"] & outcome["LR_correct"]

    cycle_features = merged.groupby(keys, as_index=False).agg(
        flux_peak_uM_h=("flux_peak_uM_h", "first"),
        flux_FWHM_h=("flux_FWHM_h", "first"),
        flux_dose_uM=("flux_dose_uM", "first"),
        post_peak_dose_uM=("post_peak_dose_uM", "first"),
        cycle_duration_generations=("cycle_duration_generations", "mean"),
    )
    cycles = outcome.merge(cycle_features, on=keys, how="left", validate="one_to_one")
    cycles["cluster_id"] = (
        cycles["knob_setting_id"].astype(str)
        + "__"
        + cycles["trajectory_index"].astype(str)
    )

    deterministic = flux.loc[flux["interface_mode"] == "deterministic_reference"].copy()
    reference = deterministic.groupby("knob_setting_id").agg(
        flux_peak_uM_h=("flux_peak_uM_h", "median"),
        flux_FWHM_h=("flux_FWHM_h", "median"),
        flux_dose_uM=("flux_dose_uM", "median"),
        post_peak_dose_uM=("post_peak_dose_uM", "median"),
    )
    for feature in FEATURES[:-1]:
        lookup = cycles["knob_setting_id"].map(reference[feature])
        cycles[f"{feature}_ratio_to_deterministic"] = cycles[feature] / lookup
        cycles[f"{feature}_absolute_log_deviation"] = np.abs(
            np.log(cycles[f"{feature}_ratio_to_deterministic"].clip(lower=1e-12))
        )
    cycles["cycle_duration_generations_ratio_to_deterministic"] = np.nan
    return cycles


def feature_attribution(cycles: pd.DataFrame, replicates: int, seed: int) -> pd.DataFrame:
    rows = []
    analysed = list(FEATURES)
    analysed += [f"{feature}_ratio_to_deterministic" for feature in FEATURES[:-1]]
    analysed += [f"{feature}_absolute_log_deviation" for feature in FEATURES[:-1]]
    for index, feature in enumerate(analysed):
        valid = cycles[["cluster_id", "both_initial_states_correct", feature]].dropna()
        failed = valid.loc[~valid["both_initial_states_correct"], feature].to_numpy(float)
        passed = valid.loc[valid["both_initial_states_correct"], feature].to_numpy(float)
        if not len(failed) or not len(passed):
            continue
        risk_auc = float(mannwhitneyu(failed, passed, alternative="two-sided").statistic)
        risk_auc /= len(failed) * len(passed)
        difference, low, high = cluster_bootstrap_median_difference(
            cycles, feature, replicates, seed + index
        )
        rows.append(
            {
                "feature": feature,
                "n_failed_cycles": len(failed),
                "n_successful_cycles": len(passed),
                "failure_median": float(np.median(failed)),
                "success_median": float(np.median(passed)),
                "failure_minus_success_median": difference,
                "cluster_bootstrap_CI_low": low,
                "cluster_bootstrap_CI_high": high,
                "risk_AUC_P_failure_value_gt_success_value": risk_auc,
                "association_strength_abs_AUC_minus_0p5": abs(risk_auc - 0.5),
            }
        )
    return pd.DataFrame(rows).sort_values(
        "association_strength_abs_AUC_minus_0p5", ascending=False
    )


def failure_summary(b_cycles: pd.DataFrame) -> pd.DataFrame:
    raw = b_cycles.loc[b_cycles["interface_mode"] == "stochastic_raw"].copy()
    raw["cycle_correct"] = as_bool(raw["cycle_correct"])
    raw["endpoint_complete"] = as_bool(raw["endpoint_complete"])
    raw["failure_category"] = raw.apply(failure_category, axis=1)
    result = (
        raw.groupby(["knob_setting_id", "initial_DNA_state", "failure_category"])
        .size()
        .rename("n_cycles")
        .reset_index()
    )
    totals = result.groupby(["knob_setting_id", "initial_DNA_state"])["n_cycles"].transform("sum")
    result["fraction"] = result["n_cycles"] / totals
    return result


def quartile_table(cycles: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for feature in FEATURES[:-1]:
        ratio = f"{feature}_ratio_to_deterministic"
        for setting, subset in cycles.groupby("knob_setting_id"):
            work = subset[[ratio, "both_initial_states_correct"]].dropna().copy()
            try:
                work["bin"] = pd.qcut(work[ratio], q=4, duplicates="drop")
            except ValueError:
                continue
            for interval, group in work.groupby("bin", observed=True):
                rows.append(
                    {
                        "knob_setting_id": setting,
                        "feature": ratio,
                        "interval": str(interval),
                        "n_cycles": len(group),
                        "both_initial_success_fraction": float(
                            group["both_initial_states_correct"].mean()
                        ),
                    }
                )
    return pd.DataFrame(rows)


def make_figures(
    cycles: pd.DataFrame,
    failures: pd.DataFrame,
    feature_table: pd.DataFrame,
    figure_dir: Path,
):
    ratio_features = [f"{feature}_ratio_to_deterministic" for feature in FEATURES[:-1]]
    fig, axes = plt.subplots(1, len(ratio_features), figsize=(17, 4.7), sharey=False)
    for axis, feature in zip(axes, ratio_features):
        passed = cycles.loc[cycles["both_initial_states_correct"], feature].dropna()
        failed = cycles.loc[~cycles["both_initial_states_correct"], feature].dropna()
        axis.boxplot([passed, failed], labels=["both pass", "failed"], showfliers=False)
        axis.axhline(1.0, color="black", linewidth=1, linestyle="--")
        base = feature.replace("_ratio_to_deterministic", "")
        axis.set_title(DISPLAY_LABELS[base])
        axis.set_ylabel("ratio to deterministic reference")
        axis.tick_params(axis="x", rotation=20)
    fig.suptitle("Stochastic A-to-B cycle outcome versus input waveform features")
    fig.tight_layout()
    fig.savefig(figure_dir / "Stochastic_AB_failure_feature_distributions.png", dpi=200)
    plt.close(fig)

    categories = [
        "correct",
        "zero crossing",
        "multiple crossings",
        "incomplete endpoint",
        "wrong direction/phase",
    ]
    pivot = failures.pivot_table(
        index=["knob_setting_id", "initial_DNA_state"],
        columns="failure_category",
        values="fraction",
        fill_value=0.0,
    ).reindex(columns=categories, fill_value=0.0)
    fig, axis = plt.subplots(figsize=(9, 5.2))
    bottom = np.zeros(len(pivot))
    for category in categories:
        values = pivot[category].to_numpy(float)
        axis.bar(range(len(pivot)), values, bottom=bottom, label=category)
        bottom += values
    axis.set_xticks(range(len(pivot)))
    axis.set_xticklabels([f"{a}\n{b}" for a, b in pivot.index], rotation=0)
    axis.set_ylim(0, 1)
    axis.set_ylabel("fraction of audited cycles")
    axis.set_title("Failure-mode composition for raw stochastic A input")
    axis.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=3)
    fig.tight_layout()
    fig.savefig(figure_dir / "Stochastic_AB_failure_mode_composition.png", dpi=200)
    plt.close(fig)

    primary = feature_table.loc[
        feature_table["feature"].isin(ratio_features)
    ].sort_values("association_strength_abs_AUC_minus_0p5")
    fig, axis = plt.subplots(figsize=(8.5, 4.8))
    colors = [
        "#d95f02" if value >= 0.5 else "#1b9e77"
        for value in primary["risk_AUC_P_failure_value_gt_success_value"]
    ]
    axis.barh(
        [DISPLAY_LABELS[name.replace("_ratio_to_deterministic", "")] for name in primary["feature"]],
        primary["risk_AUC_P_failure_value_gt_success_value"] - 0.5,
        color=colors,
    )
    axis.axvline(0.0, color="black", linewidth=1)
    axis.set_xlabel("risk AUC - 0.5  (positive: higher values enriched in failures)")
    axis.set_title("Univariate waveform associations with paired PB/LR failure")
    fig.tight_layout()
    fig.savefig(figure_dir / "Stochastic_AB_failure_association_AUC.png", dpi=200)
    plt.close(fig)


def main():
    args = parse_args()
    if args.bootstrap_replicates < 100:
        raise ValueError("--bootstrap-replicates must be at least 100")
    input_dir = newest_input_dir(args.input_dir)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = (
        args.output_dir.expanduser().resolve()
        if args.output_dir is not None
        else OUTPUT_ROOT / f"Stochastic_AB_Failure_Attribution_{stamp}"
    )
    figure_dir = FIGURE_ROOT / output_dir.name
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    b_cycles = pd.read_csv(input_dir / "B_cycle_metrics_by_interface.csv")
    flux = pd.read_csv(input_dir / "cycle_flux_metrics.csv")
    all_cycles = make_cycle_table(b_cycles, flux)
    # Primary attribution must use the unmodified stochastic interface only.
    # Deterministic-reference and causally averaged rows answer different
    # diagnostic questions and would otherwise pseudo-replicate each A cycle.
    cycles = all_cycles.loc[
        all_cycles["interface_mode"] == "stochastic_raw"
    ].copy()
    if cycles.empty:
        raise RuntimeError("No stochastic_raw cycles were found for attribution")
    features = feature_attribution(
        cycles, args.bootstrap_replicates, args.random_seed
    )
    failures = failure_summary(b_cycles)
    quartiles = quartile_table(cycles)

    cycles.to_csv(output_dir / "paired_cycle_failure_attribution.csv", index=False)
    features.to_csv(output_dir / "feature_failure_associations.csv", index=False)
    failures.to_csv(output_dir / "failure_mode_summary.csv", index=False)
    quartiles.to_csv(output_dir / "feature_quartile_success_rates.csv", index=False)
    metadata = {
        "version": VERSION,
        "input_dir": str(input_dir),
        "primary_inference_unit": "one stochastic A cycle paired across PB and LR starts",
        "primary_success_definition": "both PB and LR initial-state B simulations are correct",
        "bootstrap_unit": "A trajectory",
        "bootstrap_replicates": args.bootstrap_replicates,
        "random_seed": args.random_seed,
        "scope": "descriptive association; no causal claim and no parameter fitting",
    }
    (output_dir / "scope_and_settings.json").write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    make_figures(cycles, failures, features, figure_dir)

    print(f"Completed: {VERSION}")
    print(f"Input: {input_dir}")
    print(f"Raw tables: {output_dir}")
    print(f"Figures: {figure_dir}")
    print("\nPrimary paired-cycle result:")
    for setting, subset in cycles.groupby("knob_setting_id"):
        print(
            f"{setting}: both-initial success="
            f"{subset['both_initial_states_correct'].mean():.3f} "
            f"({int(subset['both_initial_states_correct'].sum())}/{len(subset)})"
        )
    print("\nStrongest descriptive feature associations:")
    for _, row in features.head(8).iterrows():
        direction = "higher in failures" if row["risk_AUC_P_failure_value_gt_success_value"] > 0.5 else "lower in failures"
        print(
            f"{row['feature']:<55s} | AUC={row['risk_AUC_P_failure_value_gt_success_value']:.3f} "
            f"| {direction}"
        )


if __name__ == "__main__":
    main()
