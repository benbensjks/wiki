"""Local period-knob design map for the shared-PLtetO1 repressilator-C31 ODE.

This script follows Morris and Sobol screening.  It is deliberately a local,
one-factor-at-a-time engineering map rather than another global fit.  Every
condition keeps the final construct constraint that cI and C31 use the same
PLtetO1 regulation: common K, n, leak and promoter-transcription scale.

It asks a practical question for each construct/culture change:
    Does the clock remain resolved, how does its period move, and does the
    resulting C31 translation-flux waveform remain available for B validation?

Included candidate knobs:
* synchronous oscillator RBS scale (same change on TetR, CI and LacI);
* oscillator mRNA total half-life;
* oscillator-plasmid copy number (replicon/copy-number scenario);
* sponge plasmid copies Nt, with eta fixed to 1;
* doubling time Td (culture condition, not a genetic knob).

The C31 RBS is frozen at 0.45, the current A-to-B candidate.  This program
does not itself run the B-model chemistry; its C31 flux outputs are the inputs
to be revalidated in B after choosing a candidate.

Run from the project root:
    python ./code/Shared_PLtetO1_Period_Knob_Design_Map.py --mode smoke
    python ./code/Shared_PLtetO1_Period_Knob_Design_Map.py --mode full
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


VERSION = "Shared_PLtetO1_Period_Knob_Design_Map_20260819"
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
DEFAULT_OUTPUT_ROOT = PROJECT_ROOT / "outputs"
DEFAULT_FIGURE_ROOT = PROJECT_ROOT / "figures"


@dataclass(frozen=True)
class DesignSettings:
    mode: str
    shared_PLtetO1_leak: float = 0.005
    shared_PLtetO1_tx_scale: float = 1.0
    K_copies: float = 13.0
    hill_n: float = 3.0
    sponge_accessible_fraction: float = 1.0
    c31_RBS_scale: float = 0.45
    c31_mRNA_total_half_life_min: float = 2.0
    peak_load_fraction: float = 0.0


# These are OAT design scenarios, not uncertainty confidence intervals.
# Every list deliberately includes the frozen reference value.
KNOB_VALUES = {
    "synchronous_oscillator_RBS_scale": (0.60, 0.80, 1.00, 1.20, 1.40),
    "oscillator_mRNA_total_half_life_min": (1.00, 1.50, 2.00, 3.00, 4.00),
    "oscillator_plasmid_copies": (7.0, 8.5, 10.0, 11.5, 13.0),
    "sponge_plasmid_copies_Nt": (0.0, 20.0, 40.0, 60.0, 80.0),
    "Td_min": (40.0, 45.0, 50.0, 55.0, 60.0),
}
REFERENCE_KNOB_VALUES = {
    "synchronous_oscillator_RBS_scale": 1.0,
    "oscillator_mRNA_total_half_life_min": 2.0,
    "oscillator_plasmid_copies": 10.0,
    "sponge_plasmid_copies_Nt": 40.0,
    "Td_min": 50.0,
}


def unit_for_original_spec(name: str, physical: float) -> float:
    spec = next(spec for spec in core.PARAMETER_SPECS if spec.name == name)
    if spec.scale == "linear":
        return (physical - spec.lower) / (spec.upper - spec.lower)
    if spec.scale == "log":
        return math.log(physical / spec.lower) / math.log(spec.upper / spec.lower)
    raise ValueError(f"Unknown scale {spec.scale!r}")


def build_vector(values: dict[str, float], settings: DesignSettings) -> tuple[np.ndarray, dict]:
    """Create the old solver input vector while enforcing shared PLtetO1."""
    physical = {
        "Td_min": 50.0,
        "oscillator_plasmid_copies": 10.0,
        "K_copies": settings.K_copies,
        "hill_n": settings.hill_n,
        "oscillator_promoter_leak": settings.shared_PLtetO1_leak,
        "oscillator_tx_scale": settings.shared_PLtetO1_tx_scale,
        "oscillator_translation_scale": 1.0,
        "oscillator_mrna_total_half_life_min": 2.0,
        "sponge_plasmid_copies": 40.0,
        "sponge_accessible_fraction": settings.sponge_accessible_fraction,
        # Enforced physical identity of the two PLtetO1 regulatory curves.
        "c31_promoter_leak": settings.shared_PLtetO1_leak,
        "c31_tx_scale": settings.shared_PLtetO1_tx_scale,
        "c31_translation_scale": settings.c31_RBS_scale,
        "c31_mrna_total_half_life_min": settings.c31_mRNA_total_half_life_min,
        "peak_load_fraction": settings.peak_load_fraction,
    }
    physical.update(values)
    vector = np.array([
        unit_for_original_spec(spec.name, physical[spec.name])
        for spec in core.PARAMETER_SPECS
    ], dtype=float)
    if not np.all((vector >= -1e-12) & (vector <= 1.0 + 1e-12)):
        raise ValueError("A design value is outside the frozen GSA range")
    return np.clip(vector, 0.0, 1.0), physical


def make_conditions():
    rows = []
    for knob, values in KNOB_VALUES.items():
        for value in values:
            changed = {
                "synchronous_oscillator_RBS_scale": {"oscillator_translation_scale": value},
                "oscillator_mRNA_total_half_life_min": {
                    "oscillator_mrna_total_half_life_min": value},
                "oscillator_plasmid_copies": {"oscillator_plasmid_copies": value},
                "sponge_plasmid_copies_Nt": {"sponge_plasmid_copies": value},
                "Td_min": {"Td_min": value},
            }[knob]
            rows.append({"knob": knob, "knob_value": value, "changes": changed})
    return rows


def analyse(mode: str, output_dir: Path, figure_dir: Path):
    settings = DesignSettings(mode=mode)
    sim_settings = core.MODE_SETTINGS[mode]
    base = core.make_reference_parameters()
    demand_peak = core.reference_demand_peak(base, sim_settings)
    reference_vector, _ = build_vector({}, settings)
    reference = core.simulate_run(-1, reference_vector, base, sim_settings, demand_peak)
    if reference["status"] != "resolved_oscillation":
        raise RuntimeError(f"Frozen reference did not resolve: {reference['status']}")

    conditions = make_conditions()
    rows = []
    for index, condition in enumerate(conditions, start=1):
        print(f"[{index:02d}/{len(conditions):02d}] {condition['knob']}={condition['knob_value']}")
        vector, physical = build_vector(condition["changes"], settings)
        result = core.simulate_run(index, vector, base, sim_settings, demand_peak)
        period = float(result.get("period_generations", np.nan))
        ref_period = float(reference["period_generations"])
        rows.append({
            "knob": condition["knob"],
            "knob_value": condition["knob_value"],
            **{f"physical__{key}": value for key, value in physical.items()},
            **result,
            "period_change_percent": (
                100.0 * (period / ref_period - 1.0) if np.isfinite(period) else np.nan),
            "clock_pass": bool(
                result.get("status") == "resolved_oscillation"
                and float(result.get("TetR_relative_amplitude", 0.0)) >= 0.10
            ),
        })
    table = pd.DataFrame(rows)
    table.to_csv(output_dir / "period_knob_design_results.csv", index=False)
    pd.DataFrame([reference]).to_csv(output_dir / "frozen_reference_metrics.csv", index=False)
    make_plots(table, reference, figure_dir)
    return table, reference, settings


def make_plots(table: pd.DataFrame, reference: dict, figure_dir: Path):
    figure_dir.mkdir(parents=True, exist_ok=True)
    knobs = list(KNOB_VALUES)
    fig, axes = plt.subplots(len(knobs), 2, figsize=(13, 3.3 * len(knobs)), squeeze=False)
    for row_index, knob in enumerate(knobs):
        subset = table.loc[table["knob"] == knob].sort_values("knob_value")
        ax_period, ax_interface = axes[row_index]
        color = np.where(subset["clock_pass"], "tab:blue", "tab:red")
        ax_period.plot(subset["knob_value"], subset["period_generations"], color="0.75", zorder=1)
        ax_period.scatter(subset["knob_value"], subset["period_generations"], c=color, s=42, zorder=2)
        ax_period.axhline(reference["period_generations"], color="black", ls="--", lw=1)
        ax_period.set_title(f"{knob}: clock period")
        ax_period.set_ylabel("period (generations)")
        ax_period.grid(alpha=0.25)
        ax_interface.plot(subset["knob_value"], subset["C31_flux_peak_uM_h"], "o-", label="flux peak")
        ax_interface.plot(subset["knob_value"], subset["C31_flux_FWHM_h"], "s--", label="FWHM")
        ax_interface.plot(subset["knob_value"], subset["C31_flux_residual_duty_above_0p5_uM_h"], "^:",
                          label="residual duty")
        ax_interface.set_title(f"{knob}: C31 interface diagnostics")
        ax_interface.set_ylabel("mixed units; see CSV")
        ax_interface.legend(fontsize=8)
        ax_interface.grid(alpha=0.25)
    axes[-1, 0].set_xlabel("knob value")
    axes[-1, 1].set_xlabel("knob value")
    fig.suptitle("Shared-PLtetO1 local period-knob map (blue=clock pass; red=lost/weak clock)", y=1.002)
    fig.tight_layout()
    fig.savefig(figure_dir / "Shared_PLtetO1_period_knob_design_map.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("smoke", "full"), default="smoke")
    parser.add_argument("--output-dir", type=Path, default=None)
    return parser.parse_args()


def main():
    args = parse_args()
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = args.output_dir or DEFAULT_OUTPUT_ROOT / f"Shared_PLtetO1_Period_Knob_{args.mode}_{stamp}"
    figure_dir = DEFAULT_FIGURE_ROOT / output_dir.name
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    table, reference, settings = analyse(args.mode, output_dir, figure_dir)
    scope = {
        "version": VERSION,
        "mode": args.mode,
        "shared_PLtetO1_constraint": {
            "K_shared": True,
            "hill_n_shared": True,
            "leak_shared": True,
            "transcription_scale_shared": True,
        },
        "frozen_C31_interface_candidate": {
            "C31_RBS_scale": settings.c31_RBS_scale,
            "C31_mRNA_total_half_life_min": settings.c31_mRNA_total_half_life_min,
            "peak_load_fraction": settings.peak_load_fraction,
            "B_Int_tag_note": "Int degradation tag belongs to B and is not part of this A-only ODE.",
        },
        "knob_values": KNOB_VALUES,
        "reference_knob_values": REFERENCE_KNOB_VALUES,
        "reference_metrics": reference,
        "interpretation": [
            "Blue means resolved clock with TetR relative amplitude >= 0.10, not B-module validation.",
            "Td is a growth-condition time-scale control, not a genetic period design knob.",
            "Sponge Nt is a plasmid-copy scenario with eta fixed at 1; physical TetO sites = 2*Nt.",
            "Any candidate selected here must next be sent to the B model using its C31 production flux.",
        ],
    }
    with (output_dir / "scope_and_design_constraints.json").open("w", encoding="utf-8") as handle:
        json.dump(scope, handle, indent=2, default=float)
    print(f"\nCompleted: {VERSION}")
    print(f"Raw tables: {output_dir}")
    print(f"Figures: {figure_dir}")
    print("Interpret a period change only when clock_pass=True; then revalidate its C31 flux in B.")


if __name__ == "__main__":
    main()
