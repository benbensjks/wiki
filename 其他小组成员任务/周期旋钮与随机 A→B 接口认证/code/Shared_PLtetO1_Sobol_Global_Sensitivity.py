"""Sobol global sensitivity analysis with a shared PLtetO1 constraint.

This is the second-stage analysis after Morris screening.  It reuses the
frozen mechanistic repressilator-C31 ODE in
``Mechanistic_ODE_Global_Sensitivity_Analysis.py`` but changes the *parameter
bookkeeping* to match the final plasmid design:

* cI and phiC31 are controlled by the same PLtetO1 transfer curve;
* therefore their leak fraction and promoter-transcription scale are tied;
* C31 RBS/translation and C31 mRNA lifetime remain C31-specific;
* sponge uncertainty is represented by effective TetO capacity
  S_eff = 2 * Nt * eta, rather than treating Nt and eta as independent
  Sobol factors.

It does not fit parameters or search for a visually preferred waveform.

Examples (run from this code folder):
    python Shared_PLtetO1_Sobol_Global_Sensitivity.py --analysis clock --mode smoke
    python Shared_PLtetO1_Sobol_Global_Sensitivity.py --analysis interface --mode full

The smoke mode is a workflow check only.  Full mode is the reportable design.
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


VERSION = "Shared_PLtetO1_Sobol_GSA_20260819"
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
DEFAULT_OUTPUT_ROOT = PROJECT_ROOT / "outputs"
DEFAULT_FIGURE_ROOT = PROJECT_ROOT / "figures"


@dataclass(frozen=True)
class SobolSettings:
    base_samples: int
    random_seed: int = 20260819


SETTINGS = {
    # 8 factors x (N * (D + 2)) = 80 solves; workflow check only.
    "smoke": SobolSettings(base_samples=8),
    # 8 factors x (N * (D + 2)) = 640 solves for each selected analysis.
    # This is intentionally a moderate first-order/total-order Sobol design;
    # increase --base-samples only after checking runtime locally.
    "full": SobolSettings(base_samples=64),
}


@dataclass(frozen=True)
class Factor:
    name: str
    label: str
    lower: float
    upper: float
    scale: str
    meaning: str

    def from_unit(self, x: float) -> float:
        x = float(np.clip(x, 0.0, 1.0))
        if self.scale == "linear":
            return self.lower + x * (self.upper - self.lower)
        if self.scale == "log":
            return math.exp(math.log(self.lower) + x * math.log(self.upper / self.lower))
        raise ValueError(f"Unknown scale {self.scale!r}")


# The two analyses deliberately each contain eight factors.  Their ranges
# inherit the Morris uncertainty/design ranges and are not confidence limits.
CLOCK_FACTORS = (
    Factor("Td_min", "doubling time", 40.0, 60.0, "linear", "MC4100 doubling time (min)"),
    Factor("K_copies", "common repression K", 8.0, 18.0, "linear", "effective half repression"),
    Factor("hill_n", "common Hill n", 2.0, 4.0, "linear", "shared transfer steepness"),
    Factor("shared_PLtetO1_leak", "shared PLtetO1 leak", 0.001, 0.02, "log",
           "same lower asymptote for cI and C31"),
    Factor("shared_PLtetO1_tx_scale", "shared PLtetO1 tx scale", 0.5, 1.5, "linear",
           "same multiplicative promoter scale for cI and C31"),
    Factor("oscillator_translation_scale", "oscillator RBS scale", 0.5, 1.5, "linear",
           "translation scale for the three oscillator repressors"),
    Factor("oscillator_mrna_total_half_life_min", "oscillator mRNA half-life", 1.0, 5.0,
           "log", "total oscillator mRNA half-life (min)"),
    Factor("effective_sponge_sites", "effective sponge TetO capacity", 0.0, 160.0, "linear",
           "2 * Nt * eta; represents functional decoy capacity, not plasmid count"),
)

INTERFACE_FACTORS = (
    Factor("Td_min", "doubling time", 40.0, 60.0, "linear", "shared A/B doubling time (min)"),
    Factor("hill_n", "common Hill n", 2.0, 4.0, "linear", "shared transfer steepness"),
    Factor("shared_PLtetO1_leak", "shared PLtetO1 leak", 0.001, 0.02, "log",
           "same lower asymptote for cI and C31"),
    Factor("shared_PLtetO1_tx_scale", "shared PLtetO1 tx scale", 0.5, 1.5, "linear",
           "same multiplicative promoter scale for cI and C31"),
    Factor("oscillator_mrna_total_half_life_min", "oscillator mRNA half-life", 1.0, 5.0,
           "log", "sets the upstream waveform recovery"),
    Factor("c31_translation_scale", "C31 RBS scale", 0.3, 1.2, "linear",
           "C31-specific translation/RBS choice"),
    Factor("c31_mrna_total_half_life_min", "C31 mRNA half-life", 0.5, 5.0, "log",
           "C31-specific mRNA persistence and flux tail"),
    Factor("peak_load_fraction", "peak shared-resource load", 0.0, 0.30, "linear",
           "uncalibrated scenario for dynamic translation burden"),
)


METRICS = (
    "resolved_oscillation",
    "period_generations",
    "TetR_relative_amplitude",
    "C31_flux_peak_uM_h",
    "C31_flux_FWHM_h",
    "C31_flux_dose_uM_per_cycle",
    "C31_flux_trough_fraction_of_peak_when_resolved",
    "C31_flux_residual_duty_above_0p5_uM_h_when_resolved",
    "resource_min",
)


def unit_for_original_spec(name: str, physical: float) -> float:
    """Invert the original Morris unit mapping so the frozen solver is reused."""
    spec = next(spec for spec in core.PARAMETER_SPECS if spec.name == name)
    if spec.scale == "linear":
        return (physical - spec.lower) / (spec.upper - spec.lower)
    if spec.scale == "log":
        return math.log(physical / spec.lower) / math.log(spec.upper / spec.lower)
    raise ValueError(f"Unknown original scale {spec.scale!r}")


def constrained_original_vector(unit_row: np.ndarray, factors: tuple[Factor, ...]) -> tuple[np.ndarray, dict]:
    """Map one constrained Sobol row into the original solver's 15 inputs."""
    values = {factor.name: factor.from_unit(x) for factor, x in zip(factors, unit_row)}
    # Start from frozen reference values for factors not in the selected group.
    reference = {
        "Td_min": 50.0,
        "oscillator_plasmid_copies": 10.0,
        "K_copies": 13.0,
        "hill_n": 3.0,
        "oscillator_promoter_leak": 0.005,
        "oscillator_tx_scale": 1.0,
        "oscillator_translation_scale": 1.0,
        "oscillator_mrna_total_half_life_min": 2.0,
        "sponge_plasmid_copies": 40.0,
        "sponge_accessible_fraction": 1.0,
        "c31_promoter_leak": 0.005,
        "c31_tx_scale": 1.0,
        "c31_translation_scale": 1.0,
        "c31_mrna_total_half_life_min": 2.0,
        "peak_load_fraction": 0.0,
    }
    reference.update(values)
    # Enforce the final plasmid topology: same PLtetO1 regulation for cI/C31.
    shared_leak = reference.pop("shared_PLtetO1_leak", reference["oscillator_promoter_leak"])
    shared_tx = reference.pop("shared_PLtetO1_tx_scale", reference["oscillator_tx_scale"])
    reference["oscillator_promoter_leak"] = shared_leak
    reference["c31_promoter_leak"] = shared_leak
    reference["oscillator_tx_scale"] = shared_tx
    reference["c31_tx_scale"] = shared_tx
    # A Sobol factor must be identifiable.  We sample the experimentally
    # relevant product 2*Nt*eta, fix eta=1 here, and report the equivalent Nt.
    if "effective_sponge_sites" in reference:
        capacity = reference.pop("effective_sponge_sites")
        reference["sponge_accessible_fraction"] = 1.0
        reference["sponge_plasmid_copies"] = capacity / 2.0
        reference["effective_sponge_sites_target"] = capacity

    vector = np.array([
        unit_for_original_spec(spec.name, reference[spec.name])
        for spec in core.PARAMETER_SPECS
    ], dtype=float)
    return np.clip(vector, 0.0, 1.0), reference


def make_saltelli_design(n: int, dimension: int, rng: np.random.Generator):
    """Return A, B and A_Bi matrices using independent uniform samples."""
    a = rng.random((n, dimension))
    b = rng.random((n, dimension))
    ab = []
    for index in range(dimension):
        matrix = a.copy()
        matrix[:, index] = b[:, index]
        ab.append(matrix)
    return a, b, ab


def run_design(a: np.ndarray, b: np.ndarray, ab: list[np.ndarray], factors, base, settings, demand_peak):
    rows = []
    blocks = [("A", None, a), ("B", None, b)] + [("AB", i, matrix) for i, matrix in enumerate(ab)]
    total = sum(len(matrix) for _, _, matrix in blocks)
    counter = 0
    for block, factor_index, matrix in blocks:
        for row_index, unit_row in enumerate(matrix):
            counter += 1
            vector, physical = constrained_original_vector(unit_row, factors)
            print(f"[{counter:04d}/{total:04d}] {block}" +
                  (f"_{factors[factor_index].name}" if factor_index is not None else "") +
                  f" sample={row_index + 1:03d}")
            result = core.simulate_run(counter - 1, vector, base, settings, demand_peak)
            rows.append({
                "block": block,
                "factor": factors[factor_index].name if factor_index is not None else "",
                "sample_id": row_index,
                **{f"x__{f.name}": float(unit_row[i]) for i, f in enumerate(factors)},
                **{f"physical__{key}": value for key, value in physical.items()},
                **result,
            })
    return pd.DataFrame(rows)


def sobol_summary(results: pd.DataFrame, factors: tuple[Factor, ...]) -> pd.DataFrame:
    """Compute Saltelli first and total order indices, reporting valid support."""
    a = results.loc[results["block"] == "A"].sort_values("sample_id").reset_index(drop=True)
    b = results.loc[results["block"] == "B"].sort_values("sample_id").reset_index(drop=True)
    rows = []
    for metric in METRICS:
        ya = a[metric].to_numpy(float)
        yb = b[metric].to_numpy(float)
        for factor in factors:
            abi = results.loc[(results["block"] == "AB") & (results["factor"] == factor.name)]
            yabi = abi.sort_values("sample_id")[metric].to_numpy(float)
            valid = np.isfinite(ya) & np.isfinite(yb) & np.isfinite(yabi)
            n_valid = int(np.sum(valid))
            if n_valid < 8:
                s1 = st = np.nan
                variance = np.nan
            else:
                ya_v, yb_v, yabi_v = ya[valid], yb[valid], yabi[valid]
                variance = float(np.var(np.r_[ya_v, yb_v], ddof=1))
                if not np.isfinite(variance) or variance <= 1e-15:
                    s1 = st = np.nan
                else:
                    s1 = float(np.mean(yb_v * (yabi_v - ya_v)) / variance)
                    st = float(0.5 * np.mean((ya_v - yabi_v) ** 2) / variance)
            rows.append({
                "metric": metric,
                "factor": factor.name,
                "factor_label": factor.label,
                "S1": s1,
                "ST": st,
                "interaction_gap_ST_minus_S1": st - s1 if np.isfinite(st) and np.isfinite(s1) else np.nan,
                "n_base_samples": len(a),
                "n_valid_triplets": n_valid,
                "valid_fraction": n_valid / len(a),
                "output_variance": variance,
            })
    summary = pd.DataFrame(rows)
    summary["rank_ST"] = summary.groupby("metric")["ST"].rank(ascending=False, method="min")
    return summary.sort_values(["metric", "rank_ST", "factor"])


def make_plots(summary: pd.DataFrame, figure_dir: Path, analysis: str):
    figure_dir.mkdir(parents=True, exist_ok=True)
    metrics = [metric for metric in METRICS if summary.loc[summary.metric == metric, "ST"].notna().any()]
    columns = 2
    rows = int(math.ceil(len(metrics) / columns))
    fig, axes = plt.subplots(rows, columns, figsize=(14, 4.6 * rows), squeeze=False)
    for ax, metric in zip(axes.ravel(), metrics):
        table = summary.loc[summary.metric == metric].sort_values("ST", ascending=True)
        ax.barh(table["factor_label"], table["ST"], color="tab:blue", label="total-order ST")
        ax.plot(table["S1"], np.arange(len(table)), "o", color="tab:orange", label="first-order S1")
        ax.set_title(metric)
        ax.set_xlabel("Sobol sensitivity index")
        ax.set_xlim(left=min(-0.1, float(np.nanmin(table["S1"])) - 0.05))
        ax.grid(axis="x", alpha=0.25)
        ax.legend(loc="lower right", fontsize=8)
    for ax in axes.ravel()[len(metrics):]:
        ax.axis("off")
    fig.suptitle(f"Shared-PLtetO1 Sobol GSA: {analysis}", y=1.01)
    fig.tight_layout()
    fig.savefig(figure_dir / "Shared_PLtetO1_Sobol_indices.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis", choices=("clock", "interface"), required=True,
                        help="clock: upstream timing; interface: C31 flux delivered to B")
    parser.add_argument("--mode", choices=tuple(SETTINGS), default="smoke")
    parser.add_argument("--base-samples", type=int, default=None,
                        help="override Saltelli base sample count N; total solves=N*(D+2)")
    parser.add_argument("--output-dir", type=Path, default=None)
    return parser.parse_args()


def main():
    args = parse_args()
    sobol_settings = SETTINGS[args.mode]
    n = args.base_samples or sobol_settings.base_samples
    if n < 8:
        raise ValueError("--base-samples must be at least 8")
    factors = CLOCK_FACTORS if args.analysis == "clock" else INTERFACE_FACTORS
    simulation_settings = core.MODE_SETTINGS["smoke" if args.mode == "smoke" else "full"]
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = args.output_dir or DEFAULT_OUTPUT_ROOT / f"Shared_PLtetO1_Sobol_{args.analysis}_{args.mode}_{stamp}"
    figure_dir = DEFAULT_FIGURE_ROOT / output_dir.name
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    base = core.make_reference_parameters()
    demand_peak = core.reference_demand_peak(base, simulation_settings)
    rng = np.random.default_rng(sobol_settings.random_seed)
    a, b, ab = make_saltelli_design(n, len(factors), rng)
    print(f"{VERSION}: analysis={args.analysis}, mode={args.mode}, N={n}, D={len(factors)}")
    print(f"Total ODE evaluations: {n * (len(factors) + 2)}")
    print("Shared constraint: leak and transcription scale are tied between cI and C31.")
    results = run_design(a, b, ab, factors, base, simulation_settings, demand_peak)
    summary = sobol_summary(results, factors)
    results.to_csv(output_dir / "sobol_run_metrics.csv", index=False)
    summary.to_csv(output_dir / "sobol_sensitivity_summary.csv", index=False)
    make_plots(summary, figure_dir, args.analysis)
    scope = {
        "version": VERSION,
        "analysis": args.analysis,
        "mode": args.mode,
        "base_samples": n,
        "dimension": len(factors),
        "total_ode_evaluations": n * (len(factors) + 2),
        "factors": [asdict(factor) for factor in factors],
        "shared_PLtetO1_constraint": {
            "c31_promoter_leak_equals_oscillator_promoter_leak": True,
            "c31_tx_scale_equals_oscillator_tx_scale": True,
            "K_and_hill_n_are_shared": True,
        },
        "effective_sponge_definition": "S_eff = 2 * Nt * eta; eta is fixed to 1 in the clock Sobol parameterization.",
        "limitations": [
            "Sobol indices are conditional on the stated uncertainty ranges, not fitted biological truths.",
            "NaN conditional interface metrics for lost oscillations reduce valid support; inspect valid_fraction.",
            "peak_load_fraction remains a scenario variable until cellular burden is measured.",
        ],
    }
    with (output_dir / "scope_parameters_and_method.json").open("w", encoding="utf-8") as handle:
        json.dump(scope, handle, indent=2)
    print(f"\nCompleted: {VERSION}")
    print(f"Raw tables: {output_dir}")
    print(f"Figures: {figure_dir}")
    print("Read ST first; a large ST-S1 gap indicates interactions/nonlinearity.")


if __name__ == "__main__":
    main()
