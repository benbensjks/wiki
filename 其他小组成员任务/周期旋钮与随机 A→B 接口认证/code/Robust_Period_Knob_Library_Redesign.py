"""Uncertainty-aware redesign of the repressilator period-knob library.

The period knob is intended to provide several dosing-time settings, not to
fit one predetermined period. Therefore this program searches for a discrete
library whose settings remain ordered and well separated under the SAME
uncertainty scenarios.

This stage is conditional on a model-derived PR-leak gate. The previous
confirmation showed that PR leak above approximately 0.00449 collapses the
slow settings and dominates ordering failures. The gate is not an experimental
measurement; it is a provisional engineering requirement to be tested.

The candidate designs are the already generated synchronous oscillator
RBS/mRNA grid. C31-specific parameters, the shared PLtetO1 transfer function,
and the stochastic model are not retuned here. B-interface and SSA validation
remain subsequent certification stages.
"""

from __future__ import annotations

import argparse
import itertools
import json
from dataclasses import asdict, replace
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import qmc

import Period_Knob_Robustness_Certification as robust


VERSION = "Robust_Continuous_Period_Knob_Library_Redesign_20260821"
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
OUTPUT_ROOT = PROJECT_ROOT / "outputs"
FIGURE_ROOT = PROJECT_ROOT / "figures"

PR_LEAK_GATE = 0.0044891988797823
MINIMUM_RESOLVED_FRACTION = 0.95
MINIMUM_AMPLITUDE_FRACTION = 0.95
MINIMUM_ADJACENT_GAP_GENERATIONS = 0.50
MINIMUM_OPERATIONAL_SPAN_GENERATIONS = 4.0

MODES = {
    "smoke": robust.ModeSettings(4, 3600.0, 900.0, 2.0, random_seed=20260821),
    "audit": robust.ModeSettings(16, 6000.0, 1500.0, 2.0, random_seed=20260821),
    "confirmation": robust.ModeSettings(
        64, 6000.0, 1500.0, 1.0, random_seed=20260821
    ),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=tuple(MODES), default="smoke")
    parser.add_argument("--input-grid", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--checkpoint-every", type=int, default=20)
    return parser.parse_args()


def locate_input_grid() -> Path:
    matches = sorted(
        OUTPUT_ROOT.glob(
            "Shared_PLtetO1_RBS_mRNA_2D_full_*/RBS_mRNA_2D_all_conditions.csv"
        ),
        key=lambda path: path.stat().st_mtime,
        reverse=True,
    )
    if not matches:
        raise FileNotFoundError("RBS/mRNA 2D candidate grid was not found")
    return matches[0]


def load_candidates(path: Path) -> pd.DataFrame:
    grid = pd.read_csv(path)
    required = {
        "oscillator_RBS_scale",
        "oscillator_mRNA_total_half_life_min",
        "resolved_oscillation",
        "period_generations",
        "TetR_relative_amplitude",
        "physical__c31_mrna_total_half_life_min",
    }
    missing = required - set(grid.columns)
    if missing:
        raise ValueError(f"Input grid is missing columns: {sorted(missing)}")
    grid = grid[
        (grid["resolved_oscillation"] >= 0.5)
        & np.isfinite(grid["period_generations"])
    ].copy()
    grid = grid.sort_values(
        ["oscillator_mRNA_total_half_life_min", "oscillator_RBS_scale"]
    ).drop_duplicates(
        ["oscillator_RBS_scale", "oscillator_mRNA_total_half_life_min"]
    )
    grid["candidate_id"] = grid.apply(
        lambda row: (
            f"RBS{row['oscillator_RBS_scale']:.3g}_"
            f"mRNA{row['oscillator_mRNA_total_half_life_min']:.3g}min"
        ).replace(".", "p"),
        axis=1,
    )
    return grid.reset_index(drop=True)


def gated_uncertainties() -> tuple[robust.Uncertainty, ...]:
    values = []
    for spec in robust.UNCERTAINTIES:
        if spec.name == "PR_leak":
            values.append(
                replace(
                    spec,
                    upper=PR_LEAK_GATE,
                    evidence="model-derived provisional ordering gate",
                )
            )
        else:
            values.append(spec)
    return tuple(values)


def build_scenarios(settings: robust.ModeSettings) -> pd.DataFrame:
    specs = gated_uncertainties()
    sampler = qmc.LatinHypercube(d=len(specs), seed=settings.random_seed)
    rows = []
    for scenario_id, vector in enumerate(sampler.random(settings.n_scenarios)):
        row = {"scenario_id": scenario_id}
        for spec, value in zip(specs, vector):
            row[spec.name] = spec.from_unit(float(value))
        rows.append(row)
    return pd.DataFrame(rows)


def library_row(candidate: pd.Series) -> pd.Series:
    return pd.Series(
        {
            "knob_setting_id": candidate["candidate_id"],
            "target_period_generations": float(candidate["period_generations"]),
            "synchronous_oscillator_RBS_scale": float(
                candidate["oscillator_RBS_scale"]
            ),
            "synchronous_oscillator_mRNA_half_life_min": float(
                candidate["oscillator_mRNA_total_half_life_min"]
            ),
            "C31_mRNA_total_half_life_min": float(
                candidate["physical__c31_mrna_total_half_life_min"]
            ),
        }
    )


def write_checkpoint(rows: list[dict], path: Path) -> None:
    pd.DataFrame(rows).sort_values(["scenario_id", "candidate_id"]).to_csv(
        path, index=False, encoding="utf-8-sig"
    )


def candidate_summary(metrics: pd.DataFrame, candidates: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for candidate in candidates.itertuples(index=False):
        block = metrics[metrics["candidate_id"] == candidate.candidate_id]
        periods = block["period_generations"].dropna()
        rows.append(
            {
                "candidate_id": candidate.candidate_id,
                "RBS_scale": candidate.oscillator_RBS_scale,
                "mRNA_half_life_min": candidate.oscillator_mRNA_total_half_life_min,
                "nominal_period_generations": candidate.period_generations,
                "n_scenarios": len(block),
                "resolved_fraction": float(block["resolved_oscillation"].mean()),
                "amplitude_pass_fraction": float(
                    (block["TetR_relative_amplitude"] >= 0.20).mean()
                ),
                "period_median_generations": float(periods.median()),
                "period_mean_generations": float(periods.mean()),
                "period_sd_generations": float(periods.std(ddof=1)),
                "period_cv_across_scenarios": float(
                    periods.std(ddof=1) / periods.mean()
                ),
                "period_q05_generations": float(periods.quantile(0.05)),
                "period_q95_generations": float(periods.quantile(0.95)),
                "TetR_amplitude_min": float(block["TetR_relative_amplitude"].min()),
                "C31_flux_peak_median_uM_h": float(
                    block["C31_flux_peak_uM_h"].median()
                ),
                "C31_flux_FWHM_median_h": float(block["C31_flux_FWHM_h"].median()),
                "resource_min_worst": float(block["resource_min"].min()),
            }
        )
    summary = pd.DataFrame(rows)
    summary["A_clock_eligible"] = (
        (summary["resolved_fraction"] >= MINIMUM_RESOLVED_FRACTION)
        & (summary["amplitude_pass_fraction"] >= MINIMUM_AMPLITUDE_FRACTION)
    )
    return summary


def evaluate_libraries(
    metrics: pd.DataFrame,
    summary: pd.DataFrame,
    library_size: int,
) -> pd.DataFrame:
    eligible = summary[summary["A_clock_eligible"]].copy()
    if len(eligible) < library_size:
        return pd.DataFrame()
    period_table = metrics.pivot(
        index="scenario_id", columns="candidate_id", values="period_generations"
    )
    summary_index = summary.set_index("candidate_id")
    rows = []
    for candidate_ids in itertools.combinations(eligible["candidate_id"], library_size):
        ordered_ids = sorted(
            candidate_ids,
            key=lambda cid: summary_index.loc[cid, "period_median_generations"],
        )
        periods = period_table[list(ordered_ids)]
        gaps = periods.diff(axis=1).iloc[:, 1:]
        finite = periods.notna().all(axis=1)
        strict_order = finite & (gaps > 0.0).all(axis=1)
        separated = finite & (gaps >= MINIMUM_ADJACENT_GAP_GENERATIONS).all(axis=1)
        medians = summary_index.loc[ordered_ids, "period_median_generations"]
        median_gaps = medians.diff().dropna()
        mean_median_gap = float(median_gaps.mean())
        median_gap_cv = (
            float(median_gaps.std(ddof=0) / mean_median_gap)
            if mean_median_gap > 0.0
            else float("inf")
        )
        median_span = float(medians.iloc[-1] - medians.iloc[0])
        cvs = summary_index.loc[ordered_ids, "period_cv_across_scenarios"]
        rows.append(
            {
                "library_size": library_size,
                "setting_ids": ";".join(ordered_ids),
                "RBS_scales": ";".join(
                    f"{summary_index.loc[cid, 'RBS_scale']:.6g}" for cid in ordered_ids
                ),
                "mRNA_half_lives_min": ";".join(
                    f"{summary_index.loc[cid, 'mRNA_half_life_min']:.6g}"
                    for cid in ordered_ids
                ),
                "median_periods_generations": ";".join(
                    f"{value:.6g}" for value in medians
                ),
                "all_resolved_fraction": float(finite.mean()),
                "strict_adjacent_ordering_fraction": float(strict_order.mean()),
                "gap_at_least_0p5_fraction": float(separated.mean()),
                "median_period_span_generations": median_span,
                "operational_span_at_least_4_generations": bool(
                    median_span >= MINIMUM_OPERATIONAL_SPAN_GENERATIONS
                ),
                "minimum_median_gap_generations": float(median_gaps.min()),
                "maximum_median_gap_generations": float(median_gaps.max()),
                "median_gap_cv": median_gap_cv,
                "mean_candidate_period_cv": float(cvs.mean()),
                "maximum_candidate_period_cv": float(cvs.max()),
                "B_validation_status": "pending",
                "SSA_validation_status": "pending",
            }
        )
    result = pd.DataFrame(rows)
    return result.sort_values(
        [
            "strict_adjacent_ordering_fraction",
            "gap_at_least_0p5_fraction",
            "operational_span_at_least_4_generations",
            "median_gap_cv",
            "maximum_candidate_period_cv",
            "median_period_span_generations",
            "minimum_median_gap_generations",
        ],
        ascending=[False, False, False, True, True, False, False],
    ).reset_index(drop=True)


def make_figures(
    summary: pd.DataFrame,
    best_library: pd.Series | None,
    metrics: pd.DataFrame,
    figure_dir: Path,
) -> None:
    ordered = summary.sort_values("period_median_generations")
    fig, ax = plt.subplots(figsize=(10, 6))
    low = ordered["period_median_generations"] - ordered["period_q05_generations"]
    high = ordered["period_q95_generations"] - ordered["period_median_generations"]
    colors = np.where(ordered["A_clock_eligible"], "tab:blue", "tab:gray")
    x = np.arange(len(ordered))
    ax.vlines(
        x,
        ordered["period_q05_generations"],
        ordered["period_q95_generations"],
        colors=colors,
        alpha=0.75,
    )
    ax.scatter(x, ordered["period_median_generations"], c=colors)
    ax.set_xticks(x, ordered["candidate_id"], rotation=70, ha="right")
    ax.set_ylabel("period under gated uncertainty (generations)")
    ax.set_title("Candidate medians and 5-95% uncertainty intervals")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(figure_dir / "Robust_Knob_candidate_period_intervals.png", dpi=180)
    plt.close(fig)

    if best_library is None:
        return
    ids = best_library["setting_ids"].split(";")
    pivot = metrics.pivot(
        index="scenario_id", columns="candidate_id", values="period_generations"
    )[ids]
    fig, ax = plt.subplots(figsize=(9, 5))
    for _, row in pivot.iterrows():
        ax.plot(range(len(ids)), row.to_numpy(float), color="tab:blue", alpha=0.10)
    ax.plot(range(len(ids)), pivot.median(axis=0), "o-", color="black", lw=2, label="median")
    ax.set_xticks(range(len(ids)), ids, rotation=20, ha="right")
    ax.set_ylabel("period (generations)")
    ax.set_title("Best provisional robust period-knob library")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(figure_dir / "Robust_Knob_best_library_scenarios.png", dpi=180)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    settings = MODES[args.mode]
    input_grid = args.input_grid.resolve() if args.input_grid else locate_input_grid()
    candidates = load_candidates(input_grid)
    scenarios = build_scenarios(settings)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = (
        args.output_dir.resolve()
        if args.output_dir
        else OUTPUT_ROOT / f"Robust_Period_Knob_Redesign_{args.mode}_{timestamp}"
    )
    figure_dir = FIGURE_ROOT / output_dir.name
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    checkpoint = output_dir / "candidate_scenario_metrics_checkpoint.csv"
    manifest = output_dir / "run_manifest.json"

    expected_manifest = json.loads(
        json.dumps(
            {
                "version": VERSION,
                "mode": args.mode,
                "settings": asdict(settings),
                "input_grid": str(input_grid),
                "PR_leak_gate": PR_LEAK_GATE,
                "candidate_ids": candidates["candidate_id"].tolist(),
            }
        )
    )
    if manifest.exists():
        if json.loads(manifest.read_text(encoding="utf-8")) != expected_manifest:
            raise RuntimeError("Existing output directory has a different run manifest")
    else:
        manifest.write_text(
            json.dumps(expected_manifest, ensure_ascii=False, indent=2), encoding="utf-8"
        )
        scenarios.to_csv(output_dir / "gated_uncertainty_scenarios.csv", index=False)

    completed = pd.read_csv(checkpoint).to_dict("records") if checkpoint.exists() else []
    done = {(int(row["scenario_id"]), str(row["candidate_id"])) for row in completed}
    library_rows = {
        str(row["candidate_id"]): library_row(row)
        for _, row in candidates.iterrows()
    }
    reference_peaks = {
        cid: robust.nominal_reference_peak(lib, settings)
        for cid, lib in library_rows.items()
    }
    total = len(scenarios) * len(candidates)
    new_since_checkpoint = 0

    for _, scenario in scenarios.iterrows():
        for _, candidate in candidates.iterrows():
            candidate_id = str(candidate["candidate_id"])
            key = (int(scenario["scenario_id"]), candidate_id)
            if key in done:
                continue
            print(f"[{len(done)+1:04d}/{total:04d}] scenario={key[0]:03d} candidate={candidate_id}")
            lib = library_rows[candidate_id]
            parameters = robust.parameters_for_scenario(lib, scenario)
            result = robust.simulate_scenario(parameters, reference_peaks[candidate_id], settings)
            completed.append(
                {
                    **result,
                    "scenario_id": key[0],
                    "candidate_id": candidate_id,
                    "RBS_scale": float(candidate["oscillator_RBS_scale"]),
                    "mRNA_half_life_min": float(
                        candidate["oscillator_mRNA_total_half_life_min"]
                    ),
                    "nominal_period_generations": float(candidate["period_generations"]),
                    **{
                        spec.name: float(scenario[spec.name])
                        for spec in gated_uncertainties()
                    },
                }
            )
            done.add(key)
            new_since_checkpoint += 1
            if new_since_checkpoint >= max(args.checkpoint_every, 1):
                write_checkpoint(completed, checkpoint)
                new_since_checkpoint = 0
                print(f"Checkpoint saved: {len(done)}/{total}")

    write_checkpoint(completed, checkpoint)
    metrics = pd.DataFrame(completed).sort_values(["scenario_id", "candidate_id"])
    summary = candidate_summary(metrics, candidates)
    libraries_3 = evaluate_libraries(metrics, summary, 3)
    libraries_4 = evaluate_libraries(metrics, summary, 4)
    libraries_5 = evaluate_libraries(metrics, summary, 5)

    metrics.to_csv(output_dir / "candidate_scenario_metrics.csv", index=False)
    summary.to_csv(output_dir / "candidate_robustness_summary.csv", index=False)
    libraries_3.to_csv(output_dir / "ranked_three_setting_libraries.csv", index=False)
    libraries_4.to_csv(output_dir / "ranked_four_setting_libraries.csv", index=False)
    libraries_5.to_csv(output_dir / "ranked_five_setting_libraries.csv", index=False)

    best_library = None
    if not libraries_5.empty:
        best_library = libraries_5.iloc[0]
    elif not libraries_4.empty:
        best_library = libraries_4.iloc[0]
    elif not libraries_3.empty:
        best_library = libraries_3.iloc[0]
    make_figures(summary, best_library, metrics, figure_dir)

    provenance = {
        "version": VERSION,
        "input_grid": str(input_grid),
        "mode_settings": asdict(settings),
        "PR_leak_gate": PR_LEAK_GATE,
        "PR_leak_gate_status": "model-derived provisional requirement, not measured",
        "ranking_rule": [
            "maximize strict adjacent ordering fraction",
            "maximize fraction with every adjacent gap >= 0.5 generations",
            "require at least 4 generations of operational span",
            "minimize CV of the median adjacent gaps for a continuous knob",
            "minimize maximum candidate period CV",
            "maximize median period span after continuity criteria pass",
            "maximize minimum median gap",
        ],
        "pending": ["B-interface recertification", "SSA recertification"],
    }
    (output_dir / "scope_criteria_and_provenance.json").write_text(
        json.dumps(provenance, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    print(f"\nCompleted: {VERSION}")
    print(f"Raw tables: {output_dir}")
    print(f"Figures: {figure_dir}")
    print(f"PR leak gate: <= {PR_LEAK_GATE:.6f} (provisional model-derived)")
    if best_library is not None:
        print("\nBest provisional A-clock library:")
        print(f"settings={best_library['setting_ids']}")
        print(f"median periods={best_library['median_periods_generations']}")
        print(
            f"ordering={best_library['strict_adjacent_ordering_fraction']:.3f}, "
            f"gap>=0.5={best_library['gap_at_least_0p5_fraction']:.3f}, "
            f"span={best_library['median_period_span_generations']:.3f} gen, "
            f"gap CV={best_library['median_gap_cv']:.3f}"
        )
        print("B and SSA validation remain pending.")


if __name__ == "__main__":
    main()
