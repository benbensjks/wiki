"""Reanalyse an existing Morris GSA result without solving the ODE again.

Use this once for a v2 result generated before conditional C31 leakage metrics
were added.  It reads existing CSV files, derives conditional metrics only for
resolved oscillator states, recomputes Morris summaries, and creates corrected
figures.  It performs zero ODE simulations.

Example:
    python .\\code\\Postprocess_Morris_Interface_Leak_Metrics.py `
      --input-dir "D:\\...\\outputs\\Morris_GSA_full_20260819_132518"
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from Mechanistic_ODE_Global_Sensitivity_Analysis import (
    B_INTERFACE_RESIDUAL_FLUX_THRESHOLD_uM_h,
    METRICS,
    NORMALIZATION_SCALE_OVERRIDES,
    PROJECT_ROOT,
    SimulationSettings,
    compute_elementary_effects,
    make_plots,
    summarize_effects,
)


VERSION = "Morris_GSA_interface_leak_postprocess_20260819"


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True,
                        help="Existing Morris_GSA output directory containing CSV files")
    return parser.parse_args()


def add_conditional_interface_metrics(metrics: pd.DataFrame) -> pd.DataFrame:
    required = {
        "status",
        "C31_flux_trough_fraction_of_peak",
        "C31_flux_residual_duty_above_0p5_uM_h",
    }
    missing = required.difference(metrics.columns)
    if missing:
        raise ValueError(
            "Input result predates v2 leakage outputs; missing columns: "
            + ", ".join(sorted(missing))
        )
    corrected = metrics.copy()
    resolved = corrected["status"].eq("resolved_oscillation")
    corrected["C31_flux_trough_fraction_of_peak_when_resolved"] = np.where(
        resolved, corrected["C31_flux_trough_fraction_of_peak"], np.nan
    )
    corrected["C31_flux_residual_duty_above_0p5_uM_h_when_resolved"] = np.where(
        resolved, corrected["C31_flux_residual_duty_above_0p5_uM_h"], np.nan
    )
    return corrected


def main():
    args = parse_args()
    input_dir = args.input_dir.resolve()
    metrics_path = input_dir / "morris_run_metrics.csv"
    links_path = input_dir / "morris_effect_links.csv"
    reference_path = input_dir / "reference_metrics.csv"
    for path in (metrics_path, links_path, reference_path):
        if not path.exists():
            raise FileNotFoundError(path)

    metrics = add_conditional_interface_metrics(pd.read_csv(metrics_path))
    links = pd.read_csv(links_path)
    reference = add_conditional_interface_metrics(pd.read_csv(reference_path))
    reference_row = reference.iloc[0].to_dict()
    reference_scales = {
        metric: (
            max(abs(float(reference_row.get(metric, np.nan))), 1e-12)
            if np.isfinite(float(reference_row.get(metric, np.nan))) else 1.0
        )
        for metric in METRICS
    }
    reference_scales.update(NORMALIZATION_SCALE_OVERRIDES)

    # Bootstrap settings match the full GSA and use the same seed as its
    # published workflow.  No dynamic model is evaluated here.
    settings = SimulationSettings(
        duration_min=6000.0,
        sample_step_min=1.0,
        transient_min=1500.0,
        morris_trajectories=12,
        morris_levels=4,
        bootstrap_replicates=2000,
        random_seed=20260819,
    )
    effects = compute_elementary_effects(metrics, links, reference_scales)
    summary = summarize_effects(effects, settings)
    output_dir = input_dir / "interface_leak_reanalysis_v2"
    output_dir.mkdir(exist_ok=False)
    figure_dir = PROJECT_ROOT / "figures" / f"{input_dir.name}_interface_leak_reanalysis_v2"
    figure_dir.mkdir(parents=True, exist_ok=False)

    metrics.to_csv(output_dir / "morris_run_metrics_with_conditional_interface_metrics.csv",
                   index=False, encoding="utf-8-sig")
    effects.to_csv(output_dir / "morris_elementary_effects_interface_corrected.csv",
                   index=False, encoding="utf-8-sig")
    summary.to_csv(output_dir / "morris_sensitivity_summary_interface_corrected.csv",
                   index=False, encoding="utf-8-sig")
    (summary[summary["rank_within_metric"] <= 8]
     .sort_values(["metric", "rank_within_metric"])
     .to_csv(output_dir / "top_parameters_interface_corrected.csv",
             index=False, encoding="utf-8-sig"))
    with (output_dir / "postprocess_scope.json").open("w", encoding="utf-8") as handle:
        json.dump({
            "version": VERSION,
            "source_result_directory": str(input_dir),
            "ODE_simulations_run": 0,
            "derived_metrics": [
                "C31_flux_trough_fraction_of_peak_when_resolved",
                "C31_flux_residual_duty_above_0p5_uM_h_when_resolved",
            ],
            "normalization_scale_overrides": NORMALIZATION_SCALE_OVERRIDES,
            "residual_flux_reporting_threshold_uM_h": B_INTERFACE_RESIDUAL_FLUX_THRESHOLD_uM_h,
            "interpretation": [
                "All-state residual metrics include clock failure and describe complete-system risk.",
                "Conditional residual metrics isolate C31 leakage only where the upstream clock remains resolved.",
                "The 0.5 uM/h threshold is a reporting convention, not a fitted B kinetic threshold.",
            ],
        }, handle, ensure_ascii=False, indent=2)
    make_plots(summary, figure_dir)

    print("Completed:", VERSION)
    print("ODE simulations run: 0")
    print("Corrected tables:", output_dir)
    print("Corrected figures:", figure_dir)
    for metric in (
        "C31_flux_trough_fraction_of_peak_when_resolved",
        "C31_flux_residual_duty_above_0p5_uM_h_when_resolved",
    ):
        print(f"\nTop parameters for {metric}:")
        top = summary[summary["metric"] == metric].head(6)
        for row in top.to_dict("records"):
            print(f"  {row['parameter']:<38s} mu*={row['mu_star_normalized']:.3g} "
                  f"valid={row['valid_fraction']:.2f}")


if __name__ == "__main__":
    main()
