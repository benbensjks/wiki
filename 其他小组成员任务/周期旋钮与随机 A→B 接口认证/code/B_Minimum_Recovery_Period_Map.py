"""Map the minimum recoverable input period of module B at fixed Int tag.

The script first simulates the already certified K4 A-module setting and
extracts one representative C31 translation-flux pulse.  Its peak, active
shape, FWHM and excess pulse dose are held fixed.  Only the low-input interval
between successive copies of that pulse is changed.  This separates B's
recovery/refractory time from A-side amplitude and pulse-shape changes.

Run from the project root:

    python ./code/B_Minimum_Recovery_Period_Map.py --mode smoke
    python ./code/B_Minimum_Recovery_Period_Map.py --mode full

No parameter is fitted.  The default main design uses k_tag_int=12 h^-1 and
Td=50 min through the frozen A design.  Both PB and LR initial DNA states must
pass.
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
from scipy.interpolate import CubicSpline
from scipy.integrate import trapezoid

import Mechanistic_ODE_Global_Sensitivity_Analysis as core
import Robust_Continuous_Period_Knob_AB_Certification as cert


VERSION = "B_Minimum_Recovery_Period_Map_20260821"
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
OUTPUT_ROOT = PROJECT_ROOT / "outputs"
FIGURE_ROOT = PROJECT_ROOT / "figures"

TAG_RATE_H_INV = 12.0
ACTIVE_THRESHOLD_FRACTION = 0.05
WARMUP_CYCLES = 4
MODE_SETTINGS = {
    "smoke": {
        "periods_generations": tuple(np.arange(8.0, 13.01, 1.0)),
        "audited_cycles": 5,
        "a_duration_h": 105.0,
    },
    "full": {
        "periods_generations": tuple(np.arange(7.0, 13.001, 0.25)),
        "audited_cycles": 10,
        "a_duration_h": 150.0,
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=tuple(MODE_SETTINGS), default="smoke")
    parser.add_argument("--periods-generations", default=None)
    parser.add_argument("--tag-rate", type=float, default=TAG_RATE_H_INV)
    parser.add_argument("--representative-csv", type=Path, default=None)
    parser.add_argument("--b-model-dir", type=Path, default=SCRIPT_DIR)
    return parser.parse_args()


def parse_periods(text: str | None, default: tuple[float, ...]) -> tuple[float, ...]:
    if text is None:
        return default
    values = tuple(sorted({float(x.strip()) for x in text.split(",") if x.strip()}))
    if not values or any(value <= 0.0 for value in values):
        raise ValueError("Periods must be a non-empty list of positive generation values.")
    return values


def select_k4(representatives: pd.DataFrame) -> pd.Series:
    ordered = representatives.sort_values("target_period_generations").reset_index(drop=True)
    if len(ordered) < 4:
        raise RuntimeError("At least four robust clock settings are required to select K4.")
    row = ordered.iloc[3].copy()
    row["clock_rank"] = "K4"
    return row


def extract_representative_pulse(
    time_h: np.ndarray,
    flux: np.ndarray,
    boundaries: np.ndarray,
) -> dict:
    metrics = cert.a_cycle_metrics(time_h, flux, boundaries)
    if not metrics:
        raise RuntimeError("No complete K4 cycles were available for pulse extraction.")
    median_dose = float(np.median([row["flux_dose_uM_per_cycle"] for row in metrics]))
    chosen = min(metrics, key=lambda row: abs(row["flux_dose_uM_per_cycle"] - median_dose))
    left, right = chosen["cycle_left_h"], chosen["cycle_right_h"]
    choose = (time_h >= left) & (time_h < right)
    y = np.asarray(flux[choose], float)
    dt = float(np.median(np.diff(time_h)))
    if len(y) < 10:
        raise RuntimeError("The selected K4 cycle is too short to define a pulse motif.")

    # Rotate the cycle so a trough is at phase zero.  This prevents the active
    # pulse from being split across the cycle boundary.
    trough_index = int(np.argmin(y))
    rotated = np.concatenate((y[trough_index:], y[:trough_index]))
    low = float(np.min(rotated))
    high = float(np.max(rotated))
    threshold = low + ACTIVE_THRESHOLD_FRACTION * (high - low)
    active_indices = np.flatnonzero(rotated > threshold)
    if not len(active_indices):
        raise RuntimeError("The K4 pulse never exceeded the active threshold.")
    first = max(0, int(active_indices[0]) - 1)
    last = min(len(rotated) - 1, int(active_indices[-1]) + 1)
    active_y = rotated[first:last + 1]
    active_t = np.arange(len(active_y), dtype=float) * dt
    active_duration = float(active_t[-1])
    half = low + 0.5 * (high - low)
    fwhm = float(np.sum(active_y >= half) * dt)
    excess_dose = float(trapezoid(np.maximum(active_y - low, 0.0), x=active_t))
    return {
        "active_time_h": active_t,
        "active_flux_uM_h": active_y,
        "trough_uM_h": low,
        "peak_uM_h": high,
        "active_duration_h": active_duration,
        "FWHM_h": fwhm,
        "excess_pulse_dose_uM": excess_dose,
        "source_cycle_period_h": float(right - left),
        "source_cycle_total_dose_uM": float(chosen["flux_dose_uM_per_cycle"]),
    }


def build_fixed_pulse_train(
    period_h: float,
    motif: dict,
    audited_cycles: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    active_duration = float(motif["active_duration_h"])
    if period_h <= active_duration + 0.05:
        raise ValueError(
            f"Requested period {period_h:.3f} h is not longer than the retained "
            f"active pulse ({active_duration:.3f} h)."
        )
    duration_h = (WARMUP_CYCLES + audited_cycles + 1) * period_h
    dt = min(0.01, float(np.median(np.diff(motif["active_time_h"]))))
    time_h = np.arange(0.0, duration_h + dt / 2.0, dt)
    phase = np.mod(time_h, period_h)
    pulse_start = 0.5 * (period_h - active_duration)
    local = phase - pulse_start
    flux = np.full_like(time_h, float(motif["trough_uM_h"]), dtype=float)
    inside = (local >= 0.0) & (local <= active_duration)
    flux[inside] = np.interp(
        local[inside], motif["active_time_h"], motif["active_flux_uM_h"]
    )
    boundaries = (
        np.arange(WARMUP_CYCLES, WARMUP_CYCLES + audited_cycles + 1, dtype=float)
        * period_h
    )
    if len(boundaries) != audited_cycles + 1:
        raise RuntimeError("Synthetic-cycle boundary construction failed.")
    return time_h, flux, boundaries


def make_figures(summary: pd.DataFrame, examples: dict, figure_dir: Path) -> None:
    figure_dir.mkdir(parents=True, exist_ok=True)
    ordered = summary.sort_values("period_h")
    fig, axes = plt.subplots(2, 1, figsize=(9.5, 8.0), constrained_layout=True)
    axes[0].plot(ordered["period_h"], ordered["worst_one_complete"], "o-", label="one + complete")
    axes[0].plot(ordered["period_h"], ordered["worst_zero"], "s--", label="zero")
    axes[0].plot(ordered["period_h"], ordered["worst_multiple"], "^--", label="multiple")
    axes[0].axhline(cert.SUCCESS_TARGET, color="black", ls=":", label="pass threshold")
    axes[0].set(xlabel="input period (h)", ylabel="worst PB/LR cycle fraction", ylim=(-0.03, 1.03))
    axes[0].grid(alpha=0.25)
    axes[0].legend()

    axes[1].step(ordered["period_h"], ordered["both_initial_states_pass"].astype(int), where="mid")
    axes[1].set(xlabel="input period (h)", ylabel="B pass", yticks=(0, 1), yticklabels=("fail", "pass"))
    axes[1].grid(alpha=0.25)
    fig.suptitle("B minimum recovery-period map; fixed K4 pulse, Int tag=12 h^-1")
    fig.savefig(figure_dir / "B_minimum_recovery_period_map.png", dpi=190, bbox_inches="tight")
    plt.close(fig)

    if examples:
        fig, axes = plt.subplots(len(examples), 2, figsize=(12, 3.4 * len(examples)), squeeze=False)
        for row_index, (label, data) in enumerate(examples.items()):
            lower = max(0.0, data["time"][-1] - 3.0 * data["period_h"])
            choose_a = data["time"] >= lower
            choose_b = data["b_time"] >= lower
            axes[row_index, 0].plot(data["time"][choose_a], data["flux"][choose_a], color="magenta")
            axes[row_index, 0].set(ylabel="C31 flux (uM/h)", title=label)
            axes[row_index, 1].plot(data["b_time"][choose_b], data["lr"][choose_b], color="red", label="LR")
            axes[row_index, 1].plot(data["b_time"][choose_b], data["pb"][choose_b], color="blue", label="PB")
            axes[row_index, 1].axhline(0.5, color="black", ls="--", lw=0.8)
            axes[row_index, 1].set(ylabel="DNA fraction", ylim=(-0.03, 1.03))
        axes[-1, 0].set_xlabel("time (h)")
        axes[-1, 1].set_xlabel("time (h)")
        axes[0, 1].legend()
        fig.tight_layout()
        fig.savefig(figure_dir / "B_recovery_boundary_examples.png", dpi=190, bbox_inches="tight")
        plt.close(fig)


def main() -> None:
    args = parse_args()
    mode = MODE_SETTINGS[args.mode]
    periods_gen = parse_periods(args.periods_generations, mode["periods_generations"])
    if args.tag_rate < 0.0:
        raise ValueError("The Int tag rate must be non-negative.")

    representative_csv = cert.locate_representative_csv(args.representative_csv)
    representatives = cert.load_representatives(representative_csv)
    k4 = select_k4(representatives)
    parameters, _ = cert.parameters_from_design(k4, peak_load_fraction=0.0)
    expected_period_h = float(k4["period_generations"]) * parameters.Td_min / 60.0
    source_audited_cycles = min(8, int(mode["audited_cycles"]))
    requested_a_duration_h = float(mode["a_duration_h"])
    a_duration_h = max(
        requested_a_duration_h,
        30.0 + (source_audited_cycles + 3.0) * expected_period_h,
    )
    if a_duration_h > requested_a_duration_h + 1e-9:
        print(
            f"Auto-extended K4 source simulation from {requested_a_duration_h:.1f} h "
            f"to {a_duration_h:.1f} h to obtain enough complete pulses."
        )
    base = core.make_reference_parameters()
    reference_peak = cert.reference_demand_peak(base, a_duration_h)
    a_time, _, _, _, a_flux = cert.simulate_a(
        parameters, reference_peak, a_duration_h
    )
    source_boundaries, source_period_h = cert.cycle_boundaries(
        a_time, a_flux, 30.0, source_audited_cycles, expected_period_h
    )
    motif = extract_representative_pulse(a_time, a_flux, source_boundaries)
    zc, zhao_path = cert.import_zhao_core(args.b_model_dir)

    condition_rows = []
    cycle_rows = []
    examples = {}
    total = len(periods_gen) * len(cert.DNA_INITIAL_STATES)
    counter = 0
    for period_gen in periods_gen:
        period_h = period_gen * parameters.Td_min / 60.0
        time_h, flux, boundaries = build_fixed_pulse_train(
            period_h, motif, int(mode["audited_cycles"])
        )
        spline = CubicSpline(time_h, flux, bc_type="natural")
        for initial_name in cert.DNA_INITIAL_STATES:
            counter += 1
            print(
                f"[{counter:03d}/{total:03d}] period={period_gen:.2f} gen/"
                f"{period_h:.2f} h, initial={initial_name}"
            )
            b_time, _, lr, pb, _ = cert.solve_b(
                zc, spline, float(time_h[-1]), args.tag_rate,
                parameters.Td_min, initial_name,
            )
            conservation = float(np.max(np.abs(lr + pb - 1.0)))
            cycles = cert.b_cycle_metrics(b_time, lr, pb, spline, boundaries)
            state_summary = cert.summarise_initial_state(cycles, conservation)
            condition_rows.append(
                {
                    "period_generations": period_gen,
                    "period_h": period_h,
                    "Int_tag_rate_h_inv": args.tag_rate,
                    "initial_DNA_state": initial_name,
                    **state_summary,
                }
            )
            for cycle in cycles:
                cycle_rows.append(
                    {
                        "period_generations": period_gen,
                        "period_h": period_h,
                        "initial_DNA_state": initial_name,
                        **cycle,
                    }
                )
            if initial_name == "PB" and period_gen in (periods_gen[0], periods_gen[-1]):
                examples[f"period={period_h:.2f} h"] = {
                    "time": time_h,
                    "flux": flux,
                    "b_time": b_time,
                    "lr": lr,
                    "pb": pb,
                    "period_h": period_h,
                }

    conditions = pd.DataFrame(condition_rows)
    summary_rows = []
    for (period_gen, period_h), subset in conditions.groupby(
        ["period_generations", "period_h"], sort=True
    ):
        summary_rows.append(
            {
                "period_generations": float(period_gen),
                "period_h": float(period_h),
                "both_initial_states_pass": bool(
                    len(subset) == len(cert.DNA_INITIAL_STATES)
                    and subset["pass_initial_state"].all()
                ),
                "worst_one_complete": float(subset["one_crossing_complete_fraction"].min()),
                "worst_zero": float(subset["zero_crossing_fraction"].max()),
                "worst_multiple": float(subset["multiple_crossing_fraction"].max()),
                "worst_alternation": float(subset["alternation_fidelity"].min(skipna=True)),
                "max_DNA_conservation_residual": float(subset["DNA_conservation_residual"].max()),
            }
        )
    summary = pd.DataFrame(summary_rows).sort_values("period_h").reset_index(drop=True)
    raw_passing = summary.loc[summary["both_initial_states_pass"]]
    raw_minimum = float(raw_passing["period_h"].min()) if len(raw_passing) else float("nan")
    sustained_minimum = float("nan")
    for index, row in summary.iterrows():
        if bool(row["both_initial_states_pass"]) and summary.loc[index:, "both_initial_states_pass"].all():
            sustained_minimum = float(row["period_h"])
            break

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = OUTPUT_ROOT / f"B_Minimum_Recovery_Period_{args.mode}_{stamp}"
    figure_dir = FIGURE_ROOT / output_dir.name
    output_dir.mkdir(parents=True, exist_ok=False)
    figure_dir.mkdir(parents=True, exist_ok=False)
    conditions.to_csv(output_dir / "B_recovery_initial_state_summary.csv", index=False)
    summary.to_csv(output_dir / "B_recovery_period_map.csv", index=False)
    pd.DataFrame(cycle_rows).to_csv(output_dir / "B_recovery_cycle_metrics.csv", index=False)
    make_figures(summary, examples, figure_dir)

    metadata = {
        "version": VERSION,
        "mode": args.mode,
        "source_clock": str(k4.get("setting_id", "K4")),
        "source_period_h": source_period_h,
        "Td_min": parameters.Td_min,
        "Int_tag_rate_h_inv": args.tag_rate,
        "active_threshold_fraction": ACTIVE_THRESHOLD_FRACTION,
        "fixed_pulse_metrics": {
            key: value for key, value in motif.items()
            if key not in {"active_time_h", "active_flux_uM_h"}
        },
        "raw_minimum_passing_period_h": raw_minimum,
        "sustained_minimum_passing_period_h": sustained_minimum,
        "zhao_core": str(zhao_path),
        "scope": [
            "Peak, active pulse shape, FWHM and excess pulse dose are fixed from K4.",
            "Only the inter-pulse low-input interval is changed.",
            "This is a deterministic B recovery boundary, not stochastic fidelity.",
        ],
    }
    (output_dir / "scope_and_boundary_definition.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    print(f"\nCompleted: {VERSION}")
    print(f"Raw tables: {output_dir}")
    print(f"Figures: {figure_dir}")
    print(f"Raw minimum passing period: {raw_minimum:.3f} h")
    print(f"Sustained minimum passing period: {sustained_minimum:.3f} h")
    print("Use the sustained boundary as the conservative deterministic limit.")


if __name__ == "__main__":
    main()
