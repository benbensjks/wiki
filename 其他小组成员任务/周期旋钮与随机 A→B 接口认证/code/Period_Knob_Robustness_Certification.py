"""Robustness certification for the certified T10-T13 period-knob library.

This program performs a paired uncertainty stress test.  The SAME Latin
hypercube scenario is applied to T10, T11, T12 and T13, so adjacent ordering is
tested without adding avoidable Monte Carlo noise.  It certifies four separate
questions:

1. does the clock remain resolved?
2. is the realised period within the target tolerance?
3. is the oscillation amplitude still measurable?
4. does the ordering T10 < T11 < T12 < T13 survive?

Nominal stochastic-period evidence is merged from the existing 50-trajectory
confirmation.  Robust A-to-B evidence can be merged later with
``--b-robustness-csv``; until then the complete-system status is explicitly
``pending_B_robustness`` rather than silently passing.

No parameter range in this file is a fitted confidence interval.  The defaults
are preregistered engineering uncertainty scenarios and must be replaced by
experimental ranges when available.

Run from the project root:
    python ./code/Period_Knob_Robustness_Certification.py --mode smoke
    python ./code/Period_Knob_Robustness_Certification.py --mode audit
    python ./code/Period_Knob_Robustness_Certification.py --mode confirmation

Use ``--output-dir`` to resume an interrupted audit or confirmation run.
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
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks
from scipy.stats import qmc

import Arm_Specific_Promoter_Leak_Sensitivity as arm
import Mechanistic_ODE_Global_Sensitivity_Analysis as core


VERSION = "Period_Knob_Robustness_Certification_20260820"
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
OUTPUT_ROOT = PROJECT_ROOT / "outputs"
FIGURE_ROOT = PROJECT_ROOT / "figures"


@dataclass(frozen=True)
class ModeSettings:
    n_scenarios: int
    duration_min: float
    transient_min: float
    sample_step_min: float
    random_seed: int = 20260820


MODES = {
    "smoke": ModeSettings(8, 3600.0, 900.0, 2.0),
    "audit": ModeSettings(32, 6000.0, 1500.0, 2.0),
    "confirmation": ModeSettings(128, 6000.0, 1500.0, 1.0),
}


@dataclass(frozen=True)
class Uncertainty:
    name: str
    label: str
    lower: float
    upper: float
    scale: str = "linear"
    evidence: str = "provisional engineering scenario"

    def from_unit(self, value: float) -> float:
        x = float(np.clip(value, 0.0, 1.0))
        if self.scale == "linear":
            return self.lower + x * (self.upper - self.lower)
        if self.scale == "log":
            return math.exp(math.log(self.lower) + x * math.log(self.upper / self.lower))
        raise ValueError(self.scale)


# These ranges are deliberately explicit.  They are not experimental CIs.
UNCERTAINTIES = (
    Uncertainty("Td_min", "doubling time (min)", 45.0, 55.0),
    Uncertainty("oscillator_plasmid_copies", "oscillator plasmid copies", 8.0, 12.0),
    Uncertainty("K_copies", "common effective K (copies)", 10.0, 16.0),
    Uncertainty("hill_n", "common effective Hill n", 2.5, 3.5),
    Uncertainty("PLlacO1_leak", "PLlacO1 -> TetR leak", 0.001, 0.020, "log"),
    Uncertainty("PLtetO1_leak", "PLtetO1 -> cI/C31 leak", 0.001, 0.020, "log"),
    Uncertainty("PR_leak", "PR -> LacI leak", 0.001, 0.020, "log"),
    Uncertainty("oscillator_RBS_multiplier", "oscillator RBS implementation multiplier", 0.90, 1.10),
    Uncertainty("oscillator_mrna_half_life_multiplier", "oscillator mRNA lifetime multiplier", 0.80, 1.20),
    Uncertainty("effective_sponge_sites", "effective sponge TetO capacity", 60.0, 100.0),
    Uncertainty("C31_RBS_multiplier", "C31 RBS implementation multiplier", 0.90, 1.10),
    Uncertainty("C31_mrna_half_life_multiplier", "C31 mRNA lifetime multiplier", 0.80, 1.20),
    Uncertainty("peak_load_fraction", "peak shared translation load", 0.0, 0.20),
)


CRITERIA = {
    "minimum_resolved_fraction": 0.95,
    "minimum_period_within_tolerance_fraction": 0.95,
    "minimum_amplitude_pass_fraction": 0.95,
    "minimum_adjacent_ordering_fraction": 0.95,
    "period_tolerance_generations": 0.50,
    "minimum_relative_amplitude": 0.20,
    "minimum_B_robust_pass_fraction": 0.95,
}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=tuple(MODES), default="smoke")
    parser.add_argument("--output-dir", type=Path, default=None,
                        help="reuse an existing output folder to resume")
    parser.add_argument("--b-robustness-csv", type=Path, default=None,
                        help="optional scenario-level B result with setting_id, scenario_id, B_pass")
    parser.add_argument("--stochastic-summary-csv", type=Path, default=None,
                        help="optional nominal 50-trajectory stochastic summary")
    return parser.parse_args()


def locate_certified_library() -> Path:
    matches = sorted(
        OUTPUT_ROOT.glob("Certified_Period_Knob_Library_*/certified_period_knob_library.csv"),
        key=lambda p: p.stat().st_mtime, reverse=True,
    )
    if not matches:
        raise FileNotFoundError("Certified period-knob library was not found")
    return matches[0]


def locate_stochastic_summary(explicit: Path | None) -> Path | None:
    if explicit is not None:
        path = explicit.expanduser().resolve()
        if not path.exists():
            raise FileNotFoundError(path)
        return path
    matches = sorted(
        OUTPUT_ROOT.glob("Stochastic_Period_Knob_confirmation_*/stochastic_setting_summary.csv"),
        key=lambda p: p.stat().st_mtime, reverse=True,
    )
    return matches[0] if matches else None


def build_scenarios(settings: ModeSettings) -> pd.DataFrame:
    sampler = qmc.LatinHypercube(d=len(UNCERTAINTIES), seed=settings.random_seed)
    unit = sampler.random(settings.n_scenarios)
    rows = []
    for scenario_id, vector in enumerate(unit):
        row = {"scenario_id": scenario_id}
        for spec, value in zip(UNCERTAINTIES, vector):
            row[spec.name] = spec.from_unit(float(value))
        rows.append(row)
    return pd.DataFrame(rows)


def parameters_for_scenario(library_row: pd.Series, scenario: pd.Series) -> arm.ArmParameters:
    setting = {
        "target": float(library_row["target_period_generations"]),
        "RBS": float(library_row["synchronous_oscillator_RBS_scale"]),
        "mrna_half_life": float(library_row["synchronous_oscillator_mRNA_half_life_min"]),
    }
    nominal = arm.make_arm_parameters(setting)
    p = replace(
        nominal,
        Td_min=float(scenario["Td_min"]),
        oscillator_plasmid_copies=float(scenario["oscillator_plasmid_copies"]),
        K_copies=float(scenario["K_copies"]),
        hill_n=float(scenario["hill_n"]),
        PLlacO1_leak=float(scenario["PLlacO1_leak"]),
        PLtetO1_leak=float(scenario["PLtetO1_leak"]),
        PR_leak=float(scenario["PR_leak"]),
        oscillator_translation_per_mrna_per_min=(
            nominal.oscillator_translation_per_mrna_per_min
            * float(scenario["oscillator_RBS_multiplier"])
        ),
        sponge_plasmid_copies=float(scenario["effective_sponge_sites"]) / 2.0,
        sponge_accessible_fraction=1.0,
        c31_translation_per_mrna_per_min=(
            nominal.c31_translation_per_mrna_per_min
            * float(scenario["C31_RBS_multiplier"])
        ),
        peak_load_fraction=float(scenario["peak_load_fraction"]),
    )
    osc_half_life = (
        float(library_row["synchronous_oscillator_mRNA_half_life_min"])
        * float(scenario["oscillator_mrna_half_life_multiplier"])
    )
    c31_half_life = (
        float(library_row["C31_mRNA_total_half_life_min"])
        * float(scenario["C31_mrna_half_life_multiplier"])
    )
    return replace(
        p,
        oscillator_mrna_intrinsic_loss_per_min=max(0.0, math.log(2.0) / osc_half_life - p.mu),
        c31_mrna_intrinsic_loss_per_min=max(0.0, math.log(2.0) / c31_half_life - p.mu),
    )


def nominal_reference_peak(library_row: pd.Series, settings: ModeSettings) -> float:
    nominal_setting = {
        "target": float(library_row["target_period_generations"]),
        "RBS": float(library_row["synchronous_oscillator_RBS_scale"]),
        "mrna_half_life": float(library_row["synchronous_oscillator_mRNA_half_life_min"]),
    }
    nominal = arm.make_arm_parameters(nominal_setting)
    run_settings = arm.RunSettings(
        duration_min=settings.duration_min,
        transient_min=settings.transient_min,
        sample_step_min=settings.sample_step_min,
        trajectories=2, levels=4, bootstrap_replicates=100,
    )
    return arm.reference_demand_peak(nominal, run_settings)


def simulate_scenario(p: arm.ArmParameters, reference_peak: float,
                      settings: ModeSettings) -> dict:
    time = np.arange(0.0, settings.duration_min + 0.5 * settings.sample_step_min,
                     settings.sample_step_min)
    keep = time >= settings.transient_min
    t_ss = time[keep]
    try:
        sol = solve_ivp(
            lambda t, y: arm.rhs(t, y, p, reference_peak),
            (time[0], time[-1]), arm.initial_state(p), t_eval=time,
            method="LSODA", rtol=1e-7, atol=1e-9,
        )
        if not sol.success:
            raise RuntimeError(sol.message)
        tet_free = np.array([core.free_tetR(value, p) for value in sol.y[1]])
        resource = np.array([core.resource_availability(value, p, reference_peak)
                             for value in sol.y[6]])
        flux = (60.0 * resource * p.c31_translation_per_mrna_per_min * sol.y[6]
                / 602.214076)
        metric_settings = core.SimulationSettings(
            duration_min=settings.duration_min,
            transient_min=settings.transient_min,
            sample_step_min=settings.sample_step_min,
            morris_trajectories=2, morris_levels=4, bootstrap_replicates=100,
        )
        status, period_min, period_cv, amplitude, tet_peaks = core.oscillation_metrics(
            t_ss, tet_free[keep], p, metric_settings
        )
        span = float(np.ptp(flux[keep]))
        flux_peaks, _ = find_peaks(
            flux[keep], prominence=max(0.05 * span, 1e-9),
            distance=max(1, int(round(2.0 * p.Td_min / settings.sample_step_min))),
        )
        peak = float(np.max(flux[keep]))
        return {
            "solver_success": True,
            "status": status,
            "resolved_oscillation": float(status == "resolved_oscillation"),
            "period_generations": period_min / p.Td_min if np.isfinite(period_min) else np.nan,
            "period_h": period_min / 60.0 if np.isfinite(period_min) else np.nan,
            "deterministic_period_cv": period_cv,
            "TetR_relative_amplitude": amplitude,
            "C31_flux_peak_uM_h": peak,
            "C31_flux_trough_fraction_of_peak": (
                float(np.min(flux[keep]) / max(peak, 1e-12))
            ),
            "C31_flux_FWHM_h": core.median_peak_fwhm(t_ss, flux[keep], flux_peaks) / 60.0,
            "C31_flux_dose_uM_per_cycle": core.last_cycle_dose(
                t_ss, flux[keep], tet_peaks) / 60.0,
            "resource_min": float(np.min(resource[keep])),
        }
    except Exception as exc:
        return {"solver_success": False, "status": f"failure: {type(exc).__name__}: {exc}",
                "resolved_oscillation": 0.0, "period_generations": np.nan,
                "period_h": np.nan, "deterministic_period_cv": np.nan,
                "TetR_relative_amplitude": np.nan, "C31_flux_peak_uM_h": np.nan,
                "C31_flux_trough_fraction_of_peak": np.nan, "C31_flux_FWHM_h": np.nan,
                "C31_flux_dose_uM_per_cycle": np.nan, "resource_min": np.nan}


def write_checkpoint(rows: list[dict], path: Path):
    pd.DataFrame(rows).sort_values(["scenario_id", "setting_order"]).to_csv(
        path, index=False, encoding="utf-8-sig"
    )


def summarize(metrics: pd.DataFrame, library: pd.DataFrame,
              stochastic: pd.DataFrame | None, b_results: pd.DataFrame | None):
    rows = []
    for _, lib in library.iterrows():
        setting_id = str(lib["knob_setting_id"])
        group = metrics[metrics["setting_id"] == setting_id]
        resolved_fraction = float(group["resolved_oscillation"].mean())
        period_pass_fraction = float(group["period_tolerance_pass"].mean())
        amplitude_pass_fraction = float(group["amplitude_pass"].mean())
        clock_pass = (
            resolved_fraction >= CRITERIA["minimum_resolved_fraction"]
            and period_pass_fraction >= CRITERIA["minimum_period_within_tolerance_fraction"]
            and amplitude_pass_fraction >= CRITERIA["minimum_amplitude_pass_fraction"]
        )
        b_fraction = np.nan
        b_status = "pending"
        if b_results is not None:
            block = b_results[b_results["setting_id"].astype(str) == setting_id]
            if len(block):
                b_fraction = float(block["B_pass"].astype(float).mean())
                b_status = ("pass" if b_fraction >= CRITERIA["minimum_B_robust_pass_fraction"]
                            else "fail")
        stochastic_status = "pending"
        stochastic_cv = np.nan
        if stochastic is not None:
            block = stochastic[stochastic["knob_setting_id"].astype(str) == setting_id]
            if len(block):
                stochastic_status = str(block.iloc[0]["stochastic_validation_status"])
                stochastic_cv = float(block.iloc[0]["stochastic_period_cv"])
        rows.append({
            "setting_id": setting_id,
            "target_period_generations": float(lib["target_period_generations"]),
            "n_scenarios": len(group),
            "resolved_fraction": resolved_fraction,
            "period_within_tolerance_fraction": period_pass_fraction,
            "amplitude_pass_fraction": amplitude_pass_fraction,
            "period_mean_across_scenarios": float(group["period_generations"].mean()),
            "period_sd_across_scenarios": float(group["period_generations"].std(ddof=1)),
            "period_min_across_scenarios": float(group["period_generations"].min()),
            "period_max_across_scenarios": float(group["period_generations"].max()),
            "clock_robust_pass": clock_pass,
            "nominal_stochastic_status": stochastic_status,
            "nominal_stochastic_period_cv": stochastic_cv,
            "B_robust_pass_fraction": b_fraction,
            "B_robust_status": b_status,
            "complete_system_status": (
                "pass" if clock_pass and stochastic_status == "pass" and b_status == "pass"
                else ("fail" if (not clock_pass or stochastic_status == "fail" or b_status == "fail")
                      else "pending_B_robustness")
            ),
        })
    return pd.DataFrame(rows)


def ordering_audit(metrics: pd.DataFrame):
    rows = []
    ordered_ids = list(arm.KNOB_SETTINGS)
    for scenario_id, group in metrics.groupby("scenario_id"):
        indexed = group.set_index("setting_id")
        periods = [float(indexed.loc[s, "period_generations"]) for s in ordered_ids]
        resolved = [bool(indexed.loc[s, "resolved_oscillation"]) for s in ordered_ids]
        adjacent = [resolved[i] and resolved[i + 1] and periods[i + 1] > periods[i]
                    for i in range(len(periods) - 1)]
        rows.append({"scenario_id": int(scenario_id),
                     **{f"{ordered_ids[i+1]}_gt_{ordered_ids[i]}": adjacent[i]
                        for i in range(len(adjacent))},
                     "all_adjacent_ordered": all(adjacent)})
    frame = pd.DataFrame(rows)
    fraction = float(frame["all_adjacent_ordered"].mean())
    return frame, fraction


def make_figures(metrics: pd.DataFrame, summary: pd.DataFrame,
                 ordering_fraction: float, figure_dir: Path):
    figure_dir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)
    data = [metrics[metrics["setting_id"] == setting]["period_generations"].dropna()
            for setting in arm.KNOB_SETTINGS]
    axes[0].boxplot(data, tick_labels=list(arm.KNOB_SETTINGS), showfliers=False)
    axes[0].set_ylabel("period (generations)")
    axes[0].set_title("Period distributions under shared uncertainty scenarios")
    axes[0].grid(axis="y", alpha=0.25)
    x = np.arange(len(summary))
    width = 0.25
    axes[1].bar(x - width, summary["resolved_fraction"], width, label="resolved")
    axes[1].bar(x, summary["period_within_tolerance_fraction"], width, label="period")
    axes[1].bar(x + width, summary["amplitude_pass_fraction"], width, label="amplitude")
    axes[1].axhline(0.95, color="black", linestyle="--", linewidth=1)
    axes[1].set_xticks(x, summary["setting_id"])
    axes[1].set_ylim(0, 1.05)
    axes[1].set_title(f"Robustness pass fractions; ordering={ordering_fraction:.2f}")
    axes[1].legend()
    axes[1].grid(axis="y", alpha=0.25)
    fig.savefig(figure_dir / "Period_knob_robustness_summary.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 5), constrained_layout=True)
    for setting in arm.KNOB_SETTINGS:
        block = metrics[metrics["setting_id"] == setting]
        ax.scatter(block["period_generations"], block["TetR_relative_amplitude"],
                   s=18, alpha=0.65, label=setting)
    ax.axhline(CRITERIA["minimum_relative_amplitude"], color="black", linestyle="--")
    ax.set_xlabel("period (generations)")
    ax.set_ylabel("TetR relative amplitude")
    ax.set_title("Joint period-amplitude robustness")
    ax.legend()
    ax.grid(alpha=0.25)
    fig.savefig(figure_dir / "Period_knob_period_amplitude_robustness.png", dpi=180)
    plt.close(fig)


def main():
    args = parse_args()
    settings = MODES[args.mode]
    library_path = locate_certified_library()
    library = pd.read_csv(library_path)
    library = library.sort_values("target_period_generations").reset_index(drop=True)
    scenarios = build_scenarios(settings)

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = (args.output_dir.expanduser().resolve() if args.output_dir else
                  OUTPUT_ROOT / f"Period_Knob_Robustness_{args.mode}_{stamp}")
    figure_dir = FIGURE_ROOT / output_dir.name
    output_dir.mkdir(parents=True, exist_ok=True)
    scenarios.to_csv(output_dir / "uncertainty_scenarios.csv", index=False,
                     encoding="utf-8-sig")

    checkpoint = output_dir / "robustness_run_metrics_checkpoint.csv"
    completed = pd.read_csv(checkpoint).to_dict("records") if checkpoint.exists() else []
    done = {(int(row["scenario_id"]), str(row["setting_id"])) for row in completed}
    reference_peaks = {
        str(row["knob_setting_id"]): nominal_reference_peak(row, settings)
        for _, row in library.iterrows()
    }
    total = len(scenarios) * len(library)
    position = len(done)
    for _, scenario in scenarios.iterrows():
        for setting_order, (_, lib) in enumerate(library.iterrows()):
            key = (int(scenario["scenario_id"]), str(lib["knob_setting_id"]))
            if key in done:
                continue
            position += 1
            print(f"[{position:04d}/{total:04d}] scenario={key[0]:03d} setting={key[1]}")
            p = parameters_for_scenario(lib, scenario)
            result = simulate_scenario(p, reference_peaks[key[1]], settings)
            target = float(lib["target_period_generations"])
            result.update({
                "scenario_id": key[0], "setting_id": key[1],
                "setting_order": setting_order,
                "target_period_generations": target,
                "period_error_generations": (
                    abs(float(result["period_generations"]) - target)
                    if np.isfinite(float(result["period_generations"])) else np.nan
                ),
                **{spec.name: float(scenario[spec.name]) for spec in UNCERTAINTIES},
            })
            result["period_tolerance_pass"] = bool(
                result["resolved_oscillation"]
                and np.isfinite(result["period_error_generations"])
                and result["period_error_generations"] <= CRITERIA["period_tolerance_generations"]
            )
            result["amplitude_pass"] = bool(
                result["resolved_oscillation"]
                and np.isfinite(result["TetR_relative_amplitude"])
                and result["TetR_relative_amplitude"] >= CRITERIA["minimum_relative_amplitude"]
            )
            completed.append(result)
            done.add(key)
            write_checkpoint(completed, checkpoint)

    metrics = pd.DataFrame(completed).sort_values(["scenario_id", "setting_order"])
    metrics.to_csv(output_dir / "period_knob_robustness_run_metrics.csv", index=False,
                   encoding="utf-8-sig")
    stochastic_path = locate_stochastic_summary(args.stochastic_summary_csv)
    stochastic = pd.read_csv(stochastic_path) if stochastic_path else None
    b_results = None
    if args.b_robustness_csv is not None:
        b_path = args.b_robustness_csv.expanduser().resolve()
        if not b_path.exists():
            raise FileNotFoundError(b_path)
        b_results = pd.read_csv(b_path)
        required = {"setting_id", "scenario_id", "B_pass"}
        if not required.issubset(b_results.columns):
            raise ValueError(f"B robustness CSV requires columns: {sorted(required)}")
    setting_summary = summarize(metrics, library, stochastic, b_results)
    ordering, ordering_fraction = ordering_audit(metrics)
    setting_summary["library_adjacent_ordering_fraction"] = ordering_fraction
    setting_summary["library_ordering_pass"] = (
        ordering_fraction >= CRITERIA["minimum_adjacent_ordering_fraction"]
    )
    setting_summary.to_csv(output_dir / "period_knob_robustness_setting_summary.csv",
                           index=False, encoding="utf-8-sig")
    ordering.to_csv(output_dir / "period_knob_robustness_ordering_audit.csv",
                    index=False, encoding="utf-8-sig")
    make_figures(metrics, setting_summary, ordering_fraction, figure_dir)

    scope = {
        "version": VERSION,
        "mode": args.mode,
        "mode_settings": asdict(settings),
        "uncertainties": [asdict(item) for item in UNCERTAINTIES],
        "criteria": CRITERIA,
        "certified_library": str(library_path),
        "stochastic_summary": str(stochastic_path) if stochastic_path else None,
        "B_robustness_csv": str(args.b_robustness_csv) if args.b_robustness_csv else None,
        "interpretation_boundaries": [
            "Default uncertainty ranges are engineering scenarios, not experimental confidence intervals.",
            "Nominal stochastic CV and parameter-uncertainty spread are distinct outputs.",
            "Complete-system robustness remains pending until scenario-level B evidence is merged.",
            "K and n remain common effective parameters pending arm-specific transfer curves.",
            "Smoke mode is a workflow check; audit and confirmation are evidential stages.",
        ],
    }
    with (output_dir / "scope_parameters_criteria_and_provenance.json").open(
            "w", encoding="utf-8") as handle:
        json.dump(scope, handle, ensure_ascii=False, indent=2)

    print("\nCompleted:", VERSION)
    print("Mode:", args.mode)
    print("Raw tables:", output_dir)
    print("Figures:", figure_dir)
    print(f"Adjacent ordering fraction: {ordering_fraction:.3f}")
    print("\nRobustness summary:")
    for row in setting_summary.to_dict("records"):
        print(f"{row['setting_id']}: resolved={row['resolved_fraction']:.2f}, "
              f"period={row['period_within_tolerance_fraction']:.2f}, "
              f"amplitude={row['amplitude_pass_fraction']:.2f}, "
              f"clock_pass={row['clock_robust_pass']}, "
              f"complete={row['complete_system_status']}")


if __name__ == "__main__":
    main()
