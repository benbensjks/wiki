"""Arm-specific promoter-leak sensitivity for the certified T10-T13 clock library.

Scientific question
-------------------
The earlier mechanistic ODE used one symmetric ``oscillator_promoter_leak``
for all three repression arms.  This focused analysis separates the lower
asymptotes according to the physical promoters:

* PLlacO1 -> TetR, repressed by LacI;
* PLtetO1 -> cI, repressed by free TetR;
* PR      -> LacI, repressed by cI;
* PLtetO1 -> phiC31 uses the SAME PLtetO1 leak as the cI arm.

The program runs a three-factor Morris screen at every certified T10-T13
operating point.  It is a measurement-priority analysis, not a parameter fit.

Run from this code directory:
    python Arm_Specific_Promoter_Leak_Sensitivity.py --mode smoke
    python Arm_Specific_Promoter_Leak_Sensitivity.py --mode full

Smoke mode only checks the workflow.  Full mode produces reportable rankings.
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

import Mechanistic_ODE_Global_Sensitivity_Analysis as core


VERSION = "Arm_Specific_Promoter_Leak_Sensitivity_20260820"
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
OUTPUT_ROOT = PROJECT_ROOT / "outputs"
FIGURE_ROOT = PROJECT_ROOT / "figures"


@dataclass(frozen=True)
class LeakSpec:
    name: str
    label: str
    lower: float = 0.001
    upper: float = 0.020

    def from_unit(self, value: float) -> float:
        """Map the unit interval to a provisional log-uniform leak range."""
        x = float(np.clip(value, 0.0, 1.0))
        return math.exp(math.log(self.lower) + x * math.log(self.upper / self.lower))


LEAK_SPECS = (
    LeakSpec("PLlacO1_leak", "PLlacO1 -> TetR leak"),
    LeakSpec("PLtetO1_leak", "PLtetO1 -> cI and C31 shared leak"),
    LeakSpec("PR_leak", "PR -> LacI leak"),
)


@dataclass(frozen=True)
class ArmParameters(core.Parameters):
    """Core parameters plus physical promoter-specific lower asymptotes."""

    PLlacO1_leak: float = 0.005
    PLtetO1_leak: float = 0.005
    PR_leak: float = 0.005


@dataclass(frozen=True)
class RunSettings:
    duration_min: float
    transient_min: float
    sample_step_min: float
    trajectories: int
    levels: int
    bootstrap_replicates: int
    random_seed: int = 20260820


MODES = {
    "smoke": RunSettings(3600.0, 900.0, 2.0, 2, 4, 300),
    "full": RunSettings(6000.0, 1500.0, 1.0, 16, 4, 2000),
}


KNOB_SETTINGS = {
    "Certified_T10": {"target": 10.0, "RBS": 0.60, "mrna_half_life": 2.0},
    "Certified_T11": {"target": 11.0, "RBS": 0.80, "mrna_half_life": 2.0},
    "Certified_T12": {"target": 12.0, "RBS": 1.20, "mrna_half_life": 2.0},
    "Certified_T13": {"target": 13.0, "RBS": 1.20, "mrna_half_life": 3.0},
}


METRICS = (
    "resolved_oscillation",
    "period_generations",
    "TetR_relative_amplitude",
    "CI_relative_amplitude",
    "LacI_relative_amplitude",
    "C31_flux_peak_uM_h",
    "C31_flux_trough_fraction_of_peak_when_resolved",
    "C31_flux_FWHM_h",
    "C31_flux_dose_uM_per_cycle",
)


def make_arm_parameters(setting: dict) -> ArmParameters:
    """Build one certified deterministic operating point without refitting it."""
    original = core.make_reference_parameters()
    p = ArmParameters(**asdict(original))
    p = replace(
        p,
        Td_min=50.0,
        oscillator_plasmid_copies=10.0,
        K_copies=13.0,
        hill_n=3.0,
        oscillator_translation_per_mrna_per_min=(
            original.oscillator_translation_per_mrna_per_min * setting["RBS"]
        ),
        oscillator_mrna_intrinsic_loss_per_min=max(
            0.0, math.log(2.0) / setting["mrna_half_life"] - original.mu
        ),
        sponge_plasmid_copies=40.0,
        sponge_accessible_fraction=1.0,
        c31_translation_per_mrna_per_min=(
            original.c31_translation_per_mrna_per_min * 0.45
        ),
        c31_mrna_intrinsic_loss_per_min=max(0.0, math.log(2.0) / 2.0 - original.mu),
        peak_load_fraction=0.0,
        PLlacO1_leak=0.005,
        PLtetO1_leak=0.005,
        PR_leak=0.005,
    )
    return p


def promoter_activities(tet_total: float, ci: float, laci: float, p: ArmParameters):
    tet_free = core.free_tetR(tet_total, p)
    return (
        core.repression(laci, p, p.PLlacO1_leak),
        core.repression(tet_free, p, p.PLtetO1_leak),
        core.repression(ci, p, p.PR_leak),
        core.repression(tet_free, p, p.PLtetO1_leak),
        tet_free,
    )


def rhs(_: float, y: np.ndarray, p: ArmParameters, reference_peak: float):
    m_tet, tet, m_ci, ci, m_laci, laci, m31, c31 = y
    h_tet, h_ci, h_laci, h31, _ = promoter_activities(tet, ci, laci, p)
    rho = core.resource_availability(m31, p, reference_peak)
    mrna_loss = p.oscillator_mrna_intrinsic_loss_per_min + p.mu
    protein_loss = p.oscillator_protein_intrinsic_loss_per_min + p.mu
    c31_mrna_loss = p.c31_mrna_intrinsic_loss_per_min + p.mu
    c31_protein_loss = p.c31_protein_intrinsic_loss_per_min + p.mu
    No = p.oscillator_plasmid_copies
    alpha = p.oscillator_tx_per_plasmid_per_min
    beta = p.oscillator_translation_per_mrna_per_min
    return (
        No * alpha * h_tet - mrna_loss * m_tet,
        rho * beta * m_tet - protein_loss * tet,
        No * alpha * h_ci - mrna_loss * m_ci,
        rho * beta * m_ci - protein_loss * ci,
        No * alpha * h_laci - mrna_loss * m_laci,
        rho * beta * m_laci - protein_loss * laci,
        No * p.c31_tx_per_plasmid_per_min * h31 - c31_mrna_loss * m31,
        rho * p.c31_translation_per_mrna_per_min * m31 - c31_protein_loss * c31,
    )


def initial_state(p: ArmParameters):
    tet, ci, laci = 100.0, 20.0, 300.0
    h_tet, h_ci, h_laci, h31, _ = promoter_activities(tet, ci, laci, p)
    mrna_loss = p.oscillator_mrna_intrinsic_loss_per_min + p.mu
    c31_loss = p.c31_mrna_intrinsic_loss_per_min + p.mu
    No = p.oscillator_plasmid_copies
    alpha = p.oscillator_tx_per_plasmid_per_min
    return (
        No * alpha * h_tet / mrna_loss, tet,
        No * alpha * h_ci / mrna_loss, ci,
        No * alpha * h_laci / mrna_loss, laci,
        No * p.c31_tx_per_plasmid_per_min * h31 / c31_loss, 0.0,
    )


def reference_demand_peak(p: ArmParameters, settings: RunSettings) -> float:
    time = np.arange(0.0, settings.duration_min + 0.5 * settings.sample_step_min,
                     settings.sample_step_min)
    sol = solve_ivp(
        lambda t, y: rhs(t, y, p, 1.0),
        (time[0], time[-1]), initial_state(p), t_eval=time,
        method="LSODA", rtol=1e-7, atol=1e-9,
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    keep = time >= settings.transient_min
    peak = float(np.max(p.c31_translation_per_mrna_per_min * sol.y[6, keep]))
    if peak <= 0.0 or not np.isfinite(peak):
        raise RuntimeError("Reference C31 demand peak is invalid")
    return peak


def build_design(settings: RunSettings):
    rng = np.random.default_rng(settings.random_seed)
    dimension = len(LEAK_SPECS)
    levels = np.linspace(0.0, 1.0, settings.levels)
    delta = settings.levels / (2.0 * (settings.levels - 1.0))
    rows, links = [], []
    run_id = 0
    for trajectory in range(settings.trajectories):
        directions = rng.choice((-1.0, 1.0), size=dimension)
        x = np.empty(dimension, dtype=float)
        for index, direction in enumerate(directions):
            valid = levels[(levels + direction * delta >= -1e-12)
                           & (levels + direction * delta <= 1.0 + 1e-12)]
            x[index] = float(rng.choice(valid))
        order = rng.permutation(dimension)
        rows.append({"run_id": run_id, "trajectory_id": trajectory, "step_id": 0,
                     "changed_parameter": "START",
                     **{f"x__{s.name}": x[i] for i, s in enumerate(LEAK_SPECS)}})
        before = run_id
        run_id += 1
        for step, index in enumerate(order, start=1):
            x0 = float(x[index])
            x = x.copy()
            x[index] += directions[index] * delta
            rows.append({"run_id": run_id, "trajectory_id": trajectory, "step_id": step,
                         "changed_parameter": LEAK_SPECS[index].name,
                         **{f"x__{s.name}": x[i] for i, s in enumerate(LEAK_SPECS)}})
            links.append({"trajectory_id": trajectory,
                          "parameter": LEAK_SPECS[index].name,
                          "before_run_id": before, "after_run_id": run_id,
                          "delta_x": float(x[index] - x0)})
            before = run_id
            run_id += 1
    return pd.DataFrame(rows), pd.DataFrame(links), delta


def amplitude(signal: np.ndarray) -> float:
    return float(np.ptp(signal) / max(float(np.max(np.abs(signal))), 1.0))


def simulate(unit_values: np.ndarray, base: ArmParameters, settings: RunSettings,
             reference_peak: float) -> dict:
    physical = {spec.name: spec.from_unit(x) for spec, x in zip(LEAK_SPECS, unit_values)}
    p = replace(base, **physical)
    time = np.arange(0.0, settings.duration_min + 0.5 * settings.sample_step_min,
                     settings.sample_step_min)
    keep = time >= settings.transient_min
    t_ss = time[keep]
    try:
        sol = solve_ivp(
            lambda t, y: rhs(t, y, p, reference_peak),
            (time[0], time[-1]), initial_state(p), t_eval=time,
            method="LSODA", rtol=1e-7, atol=1e-9,
        )
        if not sol.success:
            raise RuntimeError(sol.message)
        tet_free = np.array([core.free_tetR(value, p) for value in sol.y[1]])
        resource = np.array([core.resource_availability(value, p, reference_peak)
                             for value in sol.y[6]])
        flux = (60.0 * resource * p.c31_translation_per_mrna_per_min * sol.y[6]
                / (602.214076 * 1.0))
        metric_settings = core.SimulationSettings(
            duration_min=settings.duration_min,
            sample_step_min=settings.sample_step_min,
            transient_min=settings.transient_min,
            morris_trajectories=settings.trajectories,
            morris_levels=settings.levels,
            bootstrap_replicates=settings.bootstrap_replicates,
        )
        status, period_min, _, tet_amp, tet_peaks = core.oscillation_metrics(
            t_ss, tet_free[keep], p, metric_settings
        )
        span = float(np.ptp(flux[keep]))
        flux_peaks, _ = find_peaks(
            flux[keep], prominence=max(0.05 * span, 1e-9),
            distance=max(1, int(round(2.0 * p.Td_min / settings.sample_step_min))),
        )
        peak = float(np.max(flux[keep]))
        trough_fraction = float(np.min(flux[keep]) / max(peak, 1e-12))
        return {
            **physical,
            "solver_success": True,
            "status": status,
            "resolved_oscillation": float(status == "resolved_oscillation"),
            "period_generations": period_min / p.Td_min if np.isfinite(period_min) else np.nan,
            "TetR_relative_amplitude": tet_amp,
            "CI_relative_amplitude": amplitude(sol.y[3, keep]),
            "LacI_relative_amplitude": amplitude(sol.y[5, keep]),
            "C31_flux_peak_uM_h": peak,
            "C31_flux_trough_fraction_of_peak_when_resolved": (
                trough_fraction if status == "resolved_oscillation" else np.nan
            ),
            "C31_flux_FWHM_h": core.median_peak_fwhm(t_ss, flux[keep], flux_peaks) / 60.0,
            "C31_flux_dose_uM_per_cycle": core.last_cycle_dose(
                t_ss, flux[keep], tet_peaks) / 60.0,
        }
    except Exception as exc:
        return {**physical, "solver_success": False,
                "status": f"failure: {type(exc).__name__}: {exc}",
                **{metric: (0.0 if metric == "resolved_oscillation" else np.nan)
                   for metric in METRICS}}


def bootstrap_ci(values: np.ndarray, replicates: int, rng: np.random.Generator):
    values = values[np.isfinite(values)]
    if not len(values):
        return np.nan, np.nan
    estimates = [np.mean(np.abs(rng.choice(values, len(values), replace=True)))
                 for _ in range(replicates)]
    return tuple(np.quantile(estimates, (0.025, 0.975)))


def effects_and_summary(metrics: pd.DataFrame, links: pd.DataFrame, reference: dict,
                        settings: RunSettings, setting_id: str):
    indexed = metrics.set_index("run_id")
    rows = []
    for link in links.to_dict("records"):
        before = indexed.loc[int(link["before_run_id"])]
        after = indexed.loc[int(link["after_run_id"])]
        for metric in METRICS:
            y0, y1 = float(before.get(metric, np.nan)), float(after.get(metric, np.nan))
            valid = np.isfinite(y0) and np.isfinite(y1)
            raw = ((y1 - y0) / link["delta_x"]) if valid else np.nan
            if metric in {"resolved_oscillation", "C31_flux_trough_fraction_of_peak_when_resolved"}:
                scale = 1.0
            else:
                ref = float(reference.get(metric, np.nan))
                scale = max(abs(ref), 1e-12) if np.isfinite(ref) else 1.0
            rows.append({"setting_id": setting_id, **link, "metric": metric,
                         "elementary_effect": raw,
                         "normalized_elementary_effect": raw / scale if valid else np.nan,
                         "valid_pair": valid})
    effects = pd.DataFrame(rows)
    rng = np.random.default_rng(settings.random_seed + 7)
    summary_rows = []
    for (metric, parameter), group in effects.groupby(["metric", "parameter"], sort=False):
        values = group["normalized_elementary_effect"].to_numpy(float)
        valid = values[np.isfinite(values)]
        low, high = bootstrap_ci(valid, settings.bootstrap_replicates, rng)
        summary_rows.append({
            "setting_id": setting_id, "metric": metric, "parameter": parameter,
            "n_effects": len(values), "valid_fraction": len(valid) / len(values),
            "mu_normalized": float(np.mean(valid)) if len(valid) else np.nan,
            "mu_star_normalized": float(np.mean(np.abs(valid))) if len(valid) else np.nan,
            "sigma_normalized": float(np.std(valid, ddof=1)) if len(valid) > 1 else 0.0,
            "mu_star_ci95_low": low, "mu_star_ci95_high": high,
        })
    return effects, pd.DataFrame(summary_rows)


def make_figures(summary: pd.DataFrame, pooled: pd.DataFrame, figure_dir: Path):
    figure_dir.mkdir(parents=True, exist_ok=True)
    plot_metrics = ("period_generations", "TetR_relative_amplitude",
                    "C31_flux_trough_fraction_of_peak_when_resolved",
                    "C31_flux_FWHM_h")
    labels = {spec.name: spec.label for spec in LEAK_SPECS}
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)
    for ax, metric in zip(axes.flat, plot_metrics):
        block = summary[summary["metric"] == metric]
        for parameter in labels:
            data = block[block["parameter"] == parameter].set_index("setting_id")
            data = data.reindex(KNOB_SETTINGS)
            ax.plot(list(KNOB_SETTINGS), data["mu_star_normalized"], marker="o",
                    label=labels[parameter])
        ax.set_title(metric)
        ax.set_ylabel("normalized Morris mu*")
        ax.grid(alpha=0.25)
        ax.tick_params(axis="x", rotation=20)
    axes[0, 0].legend(fontsize=8)
    fig.suptitle("Arm-specific promoter leak sensitivity across T10-T13")
    fig.savefig(figure_dir / "Arm_specific_leak_sensitivity_by_setting.png", dpi=180)
    plt.close(fig)

    period = pooled[pooled["metric"] == "period_generations"].sort_values(
        "max_mu_star_across_settings", ascending=True)
    fig, ax = plt.subplots(figsize=(9, 4.8), constrained_layout=True)
    ax.barh([labels[x] for x in period["parameter"]],
            period["max_mu_star_across_settings"], color="tab:blue")
    ax.set_xlabel("maximum normalized Morris mu* across T10-T13")
    ax.set_title("Measurement priority for promoter leak: clock period")
    ax.grid(axis="x", alpha=0.25)
    fig.savefig(figure_dir / "Arm_specific_leak_measurement_priority.png", dpi=180)
    plt.close(fig)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=tuple(MODES), default="smoke")
    return parser.parse_args()


def main():
    args = parse_args()
    settings = MODES[args.mode]
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = OUTPUT_ROOT / f"Arm_Specific_Leak_{args.mode}_{stamp}"
    figure_dir = FIGURE_ROOT / output_dir.name
    output_dir.mkdir(parents=True, exist_ok=True)
    design, links, delta = build_design(settings)
    design.to_csv(output_dir / "arm_specific_leak_morris_design.csv", index=False,
                  encoding="utf-8-sig")
    links.to_csv(output_dir / "arm_specific_leak_effect_links.csv", index=False,
                 encoding="utf-8-sig")

    all_metrics, all_effects, all_summary = [], [], []
    for setting_index, (setting_id, setting) in enumerate(KNOB_SETTINGS.items(), start=1):
        print(f"[{setting_index}/4] {setting_id}")
        base = make_arm_parameters(setting)
        reference_peak = reference_demand_peak(base, settings)
        reference = simulate(np.full(3, math.log(0.005 / 0.001) / math.log(20.0)),
                             base, settings, reference_peak)
        reference.update({"setting_id": setting_id, "run_id": -1})
        rows = []
        xcols = [f"x__{spec.name}" for spec in LEAK_SPECS]
        for _, row in design.iterrows():
            result = simulate(row[xcols].to_numpy(float), base, settings, reference_peak)
            result.update({"setting_id": setting_id, "run_id": int(row["run_id"]),
                           "trajectory_id": int(row["trajectory_id"]),
                           "step_id": int(row["step_id"]),
                           "changed_parameter": row["changed_parameter"]})
            rows.append(result)
        metrics = pd.DataFrame(rows)
        effects, summary = effects_and_summary(metrics, links, reference, settings, setting_id)
        all_metrics.extend([pd.DataFrame([reference]), metrics])
        all_effects.append(effects)
        all_summary.append(summary)

    metrics = pd.concat(all_metrics, ignore_index=True)
    effects = pd.concat(all_effects, ignore_index=True)
    summary = pd.concat(all_summary, ignore_index=True)
    pooled = (summary.groupby(["metric", "parameter"], as_index=False)
              .agg(median_mu_star_across_settings=("mu_star_normalized", "median"),
                   max_mu_star_across_settings=("mu_star_normalized", "max"),
                   minimum_valid_fraction=("valid_fraction", "min")))
    pooled["rank_by_max_mu_star"] = pooled.groupby("metric")[
        "max_mu_star_across_settings"].rank(ascending=False, method="min")

    metrics.to_csv(output_dir / "arm_specific_leak_run_metrics.csv", index=False,
                   encoding="utf-8-sig")
    effects.to_csv(output_dir / "arm_specific_leak_elementary_effects.csv", index=False,
                   encoding="utf-8-sig")
    summary.to_csv(output_dir / "arm_specific_leak_sensitivity_by_setting.csv", index=False,
                   encoding="utf-8-sig")
    pooled.to_csv(output_dir / "arm_specific_leak_measurement_priority.csv", index=False,
                  encoding="utf-8-sig")
    make_figures(summary, pooled, figure_dir)

    scope = {
        "version": VERSION,
        "mode": args.mode,
        "settings": asdict(settings),
        "leak_ranges": [asdict(spec) for spec in LEAK_SPECS],
        "certified_knob_settings": KNOB_SETTINGS,
        "model_change": {
            "PLlacO1_to_TetR": "independent PLlacO1_leak",
            "PLtetO1_to_cI": "shared PLtetO1_leak",
            "PLtetO1_to_C31": "same shared PLtetO1_leak",
            "PR_to_LacI": "independent PR_leak",
        },
        "morris_delta": delta,
        "interpretation_boundaries": [
            "The 0.001-0.02 ranges are provisional uncertainty ranges, not confidence intervals.",
            "mu* ranks influence magnitude; signed mu is needed for direction.",
            "The screen prioritizes measurements and does not fit promoter parameters.",
            "Smoke-mode rankings are not scientific results.",
            "K and n remain common effective parameters in this focused leak analysis.",
        ],
    }
    with (output_dir / "scope_parameters_and_method.json").open("w", encoding="utf-8") as fh:
        json.dump(scope, fh, ensure_ascii=False, indent=2)

    print("\nCompleted:", VERSION)
    print("Raw tables:", output_dir)
    print("Figures:", figure_dir)
    print("\nPeriod-leak measurement priority (read full mode only):")
    period = pooled[pooled["metric"] == "period_generations"].sort_values(
        "rank_by_max_mu_star")
    for row in period.to_dict("records"):
        print(f"  {row['parameter']:<16s} max mu*={row['max_mu_star_across_settings']:.3g} "
              f"median mu*={row['median_mu_star_across_settings']:.3g} "
              f"valid>={row['minimum_valid_fraction']:.2f}")


if __name__ == "__main__":
    main()
