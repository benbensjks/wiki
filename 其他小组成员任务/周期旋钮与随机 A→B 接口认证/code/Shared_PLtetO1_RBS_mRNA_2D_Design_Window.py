"""Two-dimensional period-design window: synchronous oscillator RBS x mRNA lifetime.

This is the local design stage following Morris, Sobol and the one-factor
period-knob map.  It retains the final shared-PLtetO1 constraint and scans the
two strongest practical genetic candidates together:

    all-three-arm oscillator RBS scale x all-three-arm oscillator mRNA half-life.

C31 RBS remains fixed at the current A-to-B candidate (0.45).  C31 flux
metrics are exported for B-model validation but are NOT labelled as B pass/fail
inside this A-only analysis.

Run from the project root, for example:
    python ./code/Shared_PLtetO1_RBS_mRNA_2D_Design_Window.py --mode smoke
    python ./code/Shared_PLtetO1_RBS_mRNA_2D_Design_Window.py --mode full --target-period-generations 12.5
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

import Mechanistic_ODE_Global_Sensitivity_Analysis as core
import Shared_PLtetO1_Period_Knob_Design_Map as local_map


VERSION = "Shared_PLtetO1_RBS_mRNA_2D_Design_Window_20260819"
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
DEFAULT_OUTPUT_ROOT = PROJECT_ROOT / "outputs"
DEFAULT_FIGURE_ROOT = PROJECT_ROOT / "figures"

RBS_VALUES = (0.60, 0.80, 1.00, 1.20, 1.40)
MRNA_HALF_LIFE_VALUES_MIN = (1.00, 1.50, 2.00, 3.00, 4.00)


def evaluate_grid(mode: str, target_period_generations: float):
    settings = local_map.DesignSettings(mode=mode)
    simulation_settings = core.MODE_SETTINGS[mode]
    base = core.make_reference_parameters()
    demand_peak = core.reference_demand_peak(base, simulation_settings)
    reference_vector, _ = local_map.build_vector({}, settings)
    reference = core.simulate_run(-1, reference_vector, base, simulation_settings, demand_peak)
    if reference["status"] != "resolved_oscillation":
        raise RuntimeError(f"Frozen reference did not resolve: {reference['status']}")

    total = len(RBS_VALUES) * len(MRNA_HALF_LIFE_VALUES_MIN)
    rows = []
    index = 0
    for half_life in MRNA_HALF_LIFE_VALUES_MIN:
        for rbs in RBS_VALUES:
            index += 1
            print(f"[{index:02d}/{total:02d}] RBS={rbs:.2f}, oscillator_mRNA_t1/2={half_life:.2f} min")
            changes = {
                "oscillator_translation_scale": rbs,
                "oscillator_mrna_total_half_life_min": half_life,
            }
            vector, physical = local_map.build_vector(changes, settings)
            result = core.simulate_run(index, vector, base, simulation_settings, demand_peak)
            period = float(result.get("period_generations", np.nan))
            clock_pass = bool(
                result.get("status") == "resolved_oscillation"
                and float(result.get("TetR_relative_amplitude", 0.0)) >= 0.10
            )
            rows.append({
                # Prefix design coordinates so they cannot collide with the
                # lower-case model parameter columns in PowerShell/Excel.
                "design_RBS_scale": rbs,
                "design_oscillator_mRNA_half_life_min": half_life,
                "target_period_generations": target_period_generations,
                "clock_pass": clock_pass,
                "period_target_absolute_error_generations": (
                    abs(period - target_period_generations) if clock_pass and np.isfinite(period) else np.nan
                ),
                "period_change_percent_vs_frozen_reference": (
                    100.0 * (period / float(reference["period_generations"]) - 1.0)
                    if clock_pass and np.isfinite(period) else np.nan
                ),
                **{f"physical__{key}": value for key, value in physical.items()},
                **result,
            })
    table = pd.DataFrame(rows)
    return table, reference, settings


def matrix(table: pd.DataFrame, column: str):
    pivot = table.pivot(index="design_oscillator_mRNA_half_life_min", columns="design_RBS_scale", values=column)
    return pivot.reindex(index=MRNA_HALF_LIFE_VALUES_MIN, columns=RBS_VALUES)


def heatmap(ax, data: pd.DataFrame, title: str, label: str, fmt: str = ".2f", cmap="viridis"):
    image = ax.imshow(data.to_numpy(float), origin="lower", aspect="auto", cmap=cmap)
    ax.set_xticks(range(len(data.columns)), [f"{x:g}" for x in data.columns])
    ax.set_yticks(range(len(data.index)), [f"{y:g}" for y in data.index])
    ax.set_xlabel("synchronous oscillator RBS scale")
    ax.set_ylabel("oscillator mRNA half-life (min)")
    ax.set_title(title)
    for row in range(data.shape[0]):
        for column in range(data.shape[1]):
            value = data.iat[row, column]
            text = "NA" if not np.isfinite(value) else format(value, fmt)
            ax.text(column, row, text, ha="center", va="center", fontsize=8,
                    color="white" if np.isfinite(value) and value > np.nanmedian(data.to_numpy(float)) else "black")
    colorbar = plt.colorbar(image, ax=ax, shrink=0.86)
    colorbar.set_label(label)


def make_plots(table: pd.DataFrame, reference: dict, figure_dir: Path, target: float):
    figure_dir.mkdir(parents=True, exist_ok=True)
    panels = (
        ("period_generations", f"period (target={target:g} gen)", "generations", ".2f", "viridis"),
        ("period_target_absolute_error_generations", "absolute target error", "generations", ".2f", "magma_r"),
        ("TetR_relative_amplitude", "TetR relative amplitude", "relative amplitude", ".2f", "viridis"),
        ("C31_flux_peak_uM_h", "C31 flux peak", "uM/h", ".2f", "plasma"),
        ("C31_flux_FWHM_h", "C31 flux FWHM", "h", ".2f", "cividis"),
        ("C31_flux_dose_uM_per_cycle", "C31 flux dose/cycle", "uM/cycle", ".2f", "magma"),
    )
    fig, axes = plt.subplots(2, 3, figsize=(15, 8), squeeze=False)
    for ax, panel in zip(axes.ravel(), panels):
        column, title, label, fmt, cmap = panel
        heatmap(ax, matrix(table, column), title, label, fmt, cmap)
    fig.suptitle(
        "Shared-PLtetO1 local design window: all-three-arm RBS x all-three-arm mRNA lifetime\n"
        f"Frozen reference = {reference['period_generations']:.2f} generations; C31 RBS scale=0.45",
        y=1.02,
    )
    fig.tight_layout()
    fig.savefig(figure_dir / "Shared_PLtetO1_RBS_mRNA_2D_design_window.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("smoke", "full"), default="smoke")
    parser.add_argument("--target-period-generations", type=float, default=12.5,
                        help="planning target used only to rank candidates; default=12.5")
    parser.add_argument("--output-dir", type=Path, default=None)
    return parser.parse_args()


def main():
    args = parse_args()
    if args.target_period_generations <= 0:
        raise ValueError("target period must be positive")
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = args.output_dir or DEFAULT_OUTPUT_ROOT / f"Shared_PLtetO1_RBS_mRNA_2D_{args.mode}_{stamp}"
    figure_dir = DEFAULT_FIGURE_ROOT / output_dir.name
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    print(f"{VERSION}: mode={args.mode}, target={args.target_period_generations:g} generations")
    table, reference, settings = evaluate_grid(args.mode, args.target_period_generations)
    ranked = table.loc[table["clock_pass"]].sort_values("period_target_absolute_error_generations")
    table.to_csv(output_dir / "RBS_mRNA_2D_all_conditions.csv", index=False)
    ranked.to_csv(output_dir / "RBS_mRNA_2D_ranked_clock_candidates.csv", index=False)
    pd.DataFrame([reference]).to_csv(output_dir / "frozen_reference_metrics.csv", index=False)
    make_plots(table, reference, figure_dir, args.target_period_generations)
    scope = {
        "version": VERSION,
        "mode": args.mode,
        "target_period_generations": args.target_period_generations,
        "grid": {"RBS_scale": RBS_VALUES, "oscillator_mRNA_total_half_life_min": MRNA_HALF_LIFE_VALUES_MIN},
        "shared_PLtetO1_constraint": {
            "K_shared": True, "hill_n_shared": True, "leak_shared": True,
            "transcription_scale_shared": True,
        },
        "frozen_C31_interface": {
            "C31_RBS_scale": settings.c31_RBS_scale,
            "C31_mRNA_total_half_life_min": settings.c31_mRNA_total_half_life_min,
            "peak_load_fraction": settings.peak_load_fraction,
        },
        "interpretation": [
            "The target-period ranking is an A-module design ranking, not B compatibility certification.",
            "All top candidates must be passed to B as beta31*m31 flux before construct selection.",
            "A design is biologically interpretable only if all three oscillator arms are changed symmetrically.",
        ],
    }
    with (output_dir / "scope_and_design_constraints.json").open("w", encoding="utf-8") as handle:
        json.dump(scope, handle, indent=2, default=float)
    print(f"\nCompleted: {VERSION}")
    print(f"Raw tables: {output_dir}")
    print(f"Figures: {figure_dir}")
    print("The ranking is by clock target only.  Revalidate C31 flux in B before selecting a construct.")


if __name__ == "__main__":
    main()
