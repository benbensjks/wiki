"""Certify the robust continuous period-knob library against Zhao module B.

The script reads the top five-setting library produced by
``Robust_Period_Knob_Library_Redesign.py``. For every frozen period setting it
re-solves module A, converts the actual C31 translation flux
``rho * beta31 * m31`` to uM/h, and supplies that flux directly to Zhao module
B as the Int production term.

No period-knob parameter is fitted here.  C31 and cI retain the same PLtetO1
transfer curve.  Module B owns Int binding, dilution and tag-mediated loss.
Both PB and LR DNA initial states are required to pass.

Run from the project root:

    python ./code/Robust_Continuous_Period_Knob_AB_Certification.py --mode smoke
    python ./code/Robust_Continuous_Period_Knob_AB_Certification.py --mode full

The full mode is the reportable interface certification.  Smoke mode only
checks the workflow on a reduced panel.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import asdict, replace
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp, trapezoid
from scipy.interpolate import CubicSpline
from scipy.signal import find_peaks

import Mechanistic_ODE_Global_Sensitivity_Analysis as core
import Shared_PLtetO1_Period_Knob_Design_Map as local_map


VERSION = "Robust_Continuous_Period_Knob_AB_Certification_20260821"
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
OUTPUT_ROOT = PROJECT_ROOT / "outputs"
FIGURE_ROOT = PROJECT_ROOT / "figures"


MODE_SETTINGS = {
    "smoke": {
        "duration_h": 75.0,
        "transient_h": 20.0,
        "audited_cycles": 4,
        "tag_rates_h_inv": (12.0,),
        "maximum_representatives": 3,
    },
    "full": {
        "duration_h": 160.0,
        "transient_h": 30.0,
        "audited_cycles": 10,
        "tag_rates_h_inv": (8.0, 12.0),
        "maximum_representatives": None,
    },
}

A_SAMPLE_STEP_MIN = 2.0
B_SAMPLE_STEP_H = 0.04
B_MAX_STEP_H = 0.08
DNA_INITIAL_STATES = ("PB", "LR")
DNA_MIDPOINT = 0.5
COMPLETE_STATE_THRESHOLD = 0.95
SUCCESS_TARGET = 0.95
ALTERNATION_TARGET = 0.95
DNA_CONSERVATION_TOLERANCE = 1e-8
PERIOD_TOLERANCE_GENERATIONS = 0.35


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=tuple(MODE_SETTINGS), default="smoke")
    parser.add_argument(
        "--representative-csv",
        type=Path,
        default=None,
        help="Ranked five-setting robust library; newest matching table is used by default.",
    )
    parser.add_argument(
        "--b-model-dir",
        type=Path,
        default=SCRIPT_DIR,
        help="Folder containing zhao_core.py. The delivery copies it beside this script.",
    )
    parser.add_argument(
        "--tag-rates",
        default=None,
        help="Optional comma-separated Int tag rates in h^-1, overriding the mode defaults.",
    )
    return parser.parse_args()


def parse_float_list(text: str | None, default: tuple[float, ...]) -> tuple[float, ...]:
    if text is None:
        return default
    values = tuple(float(item.strip()) for item in text.split(",") if item.strip())
    if not values or any(value < 0.0 for value in values):
        raise ValueError("Tag rates must be a non-empty list of non-negative values.")
    return values


def locate_representative_csv(explicit: Path | None) -> Path:
    if explicit is not None:
        path = explicit.expanduser().resolve()
        if not path.exists():
            raise FileNotFoundError(path)
        return path
    matches = sorted(
        OUTPUT_ROOT.glob(
            "Robust_Period_Knob_Redesign_confirmation_*/ranked_five_setting_libraries.csv"
        ),
        key=lambda path: path.stat().st_mtime,
        reverse=True,
    )
    if not matches:
        raise FileNotFoundError(
            "No confirmed robust five-setting library was found. "
            "Pass --representative-csv explicitly."
        )
    return matches[0]


def load_representatives(path: Path) -> pd.DataFrame:
    """Expand the top ranked robust-library row into five design rows."""
    ranked = pd.read_csv(path)
    legacy_required = {
        "target_period_generations",
        "period_generations",
        "period_error_generations",
        "within_period_tolerance",
        "design_RBS_scale",
        "design_oscillator_mRNA_half_life_min",
    }
    if legacy_required.issubset(ranked.columns):
        return ranked.copy()
    required = {
        "setting_ids", "RBS_scales", "mRNA_half_lives_min",
        "median_periods_generations",
    }
    missing = required.difference(ranked.columns)
    if missing or ranked.empty:
        raise KeyError(f"Robust library table lacks required columns: {sorted(missing)}")
    top = ranked.iloc[0]
    setting_ids = str(top["setting_ids"]).split(";")
    rbs = [float(value) for value in str(top["RBS_scales"]).split(";")]
    half_lives = [float(value) for value in str(top["mRNA_half_lives_min"]).split(";")]
    medians = [
        float(value) for value in str(top["median_periods_generations"]).split(";")
    ]
    if not (len(setting_ids) == len(rbs) == len(half_lives) == len(medians)):
        raise ValueError("The robust-library row contains inconsistent list lengths.")

    summary_path = path.parent / "candidate_robustness_summary.csv"
    if not summary_path.is_file():
        raise FileNotFoundError(
            f"Candidate summary required for nominal periods was not found: {summary_path}"
        )
    summary = pd.read_csv(summary_path).set_index("candidate_id")
    rows = []
    for setting_id, rbs_scale, half_life, median in zip(
        setting_ids, rbs, half_lives, medians
    ):
        if setting_id not in summary.index:
            raise KeyError(f"{setting_id} is absent from {summary_path.name}")
        nominal = float(summary.loc[setting_id, "nominal_period_generations"])
        error = nominal - median
        rows.append(
            {
                "setting_id": setting_id,
                "target_period_generations": median,
                "period_generations": nominal,
                "period_error_generations": error,
                "within_period_tolerance": abs(error)
                <= PERIOD_TOLERANCE_GENERATIONS,
                "design_RBS_scale": rbs_scale,
                "design_oscillator_mRNA_half_life_min": half_life,
            }
        )
    return pd.DataFrame(rows)


def import_zhao_core(model_dir: Path):
    directory = model_dir.expanduser().resolve()
    module_path = directory / "zhao_core.py"
    if not module_path.is_file():
        raise FileNotFoundError(
            f"zhao_core.py was not found in {directory}. "
            "Keep the delivered zhao_core.py beside this script or pass --b-model-dir."
        )
    if str(directory) not in sys.path:
        sys.path.insert(0, str(directory))
    import zhao_core as zc

    required = (
        "default_params", "_rate_constants", "rhs", "y0_PB_ss", "y0_LR_ss",
        "LR_total", "PB_total", "int_total", "rdf_total",
    )
    missing = [name for name in required if not hasattr(zc, name)]
    if missing:
        raise ImportError(f"The selected zhao_core.py lacks the required API: {missing}")
    return zc, module_path


def select_smoke_panel(frame: pd.DataFrame, maximum: int | None) -> pd.DataFrame:
    frame = frame.sort_values("target_period_generations").reset_index(drop=True)
    if maximum is None or len(frame) <= maximum:
        return frame
    indices = np.unique(np.round(np.linspace(0, len(frame) - 1, maximum)).astype(int))
    return frame.iloc[indices].reset_index(drop=True)


def parameters_from_design(row: pd.Series, peak_load_fraction: float = 0.0):
    settings = local_map.DesignSettings(mode="full")
    changes = {
        "oscillator_translation_scale": float(row["design_RBS_scale"]),
        "oscillator_mrna_total_half_life_min": float(
            row["design_oscillator_mRNA_half_life_min"]
        ),
        "peak_load_fraction": float(peak_load_fraction),
    }
    vector, physical = local_map.build_vector(changes, settings)
    base = core.make_reference_parameters()
    parameters, reconstructed = core.apply_parameter_vector(vector, base)

    # The final construct constraint is asserted, not silently assumed.
    if not math.isclose(
        parameters.c31_promoter_leak,
        parameters.oscillator_promoter_leak,
        rel_tol=0.0,
        abs_tol=1e-12,
    ):
        raise RuntimeError("C31 and oscillator PLtetO1 leak values are not shared.")
    if not math.isclose(
        reconstructed["c31_tx_scale"],
        reconstructed["oscillator_tx_scale"],
        rel_tol=0.0,
        abs_tol=1e-12,
    ):
        raise RuntimeError("C31 and oscillator PLtetO1 transcription scales are not shared.")
    return parameters, physical


def reference_demand_peak(base, duration_h: float) -> float:
    time_min = np.arange(0.0, 60.0 * duration_h + A_SAMPLE_STEP_MIN / 2.0, A_SAMPLE_STEP_MIN)
    p = replace(base, peak_load_fraction=0.0)
    solution = solve_ivp(
        lambda t, y: core.rhs(t, y, p, 1.0),
        (float(time_min[0]), float(time_min[-1])),
        core.initial_state(p),
        t_eval=time_min,
        method="LSODA",
        rtol=1e-7,
        atol=1e-9,
    )
    if not solution.success:
        raise RuntimeError("A reference solve failed: " + solution.message)
    keep = time_min >= min(30.0 * 60.0, 0.35 * time_min[-1])
    demand = p.c31_translation_per_mrna_per_min * solution.y[6]
    peak = float(np.max(demand[keep]))
    if not np.isfinite(peak) or peak <= 0.0:
        raise RuntimeError("Reference C31 demand peak is not positive.")
    return peak


def simulate_a(parameters, reference_peak: float, duration_h: float):
    time_min = np.arange(0.0, 60.0 * duration_h + A_SAMPLE_STEP_MIN / 2.0, A_SAMPLE_STEP_MIN)
    solution = solve_ivp(
        lambda t, y: core.rhs(t, y, parameters, reference_peak),
        (float(time_min[0]), float(time_min[-1])),
        core.initial_state(parameters),
        t_eval=time_min,
        method="LSODA",
        rtol=1e-7,
        atol=1e-9,
    )
    if not solution.success:
        raise RuntimeError("A solve failed: " + solution.message)
    tet_free = np.array([core.free_tetR(total, parameters) for total in solution.y[1]])
    resource = np.array(
        [core.resource_availability(m31, parameters, reference_peak) for m31 in solution.y[6]]
    )
    requested = parameters.c31_translation_per_mrna_per_min * solution.y[6]
    flux_uM_h = (
        60.0
        * resource
        * requested
        * core.MODE_SETTINGS["full"].active_C31_fraction
        / (
            core.MODE_SETTINGS["full"].molecules_per_uM_per_fL
            * core.MODE_SETTINGS["full"].cell_volume_fL
        )
    )
    return time_min / 60.0, solution, tet_free, resource, flux_uM_h


def cycle_boundaries(
    time_h: np.ndarray,
    flux: np.ndarray,
    transient_h: float,
    audited_cycles: int,
    expected_period_h: float,
) -> tuple[np.ndarray, float]:
    keep = time_h >= transient_h
    indices = np.flatnonzero(keep)
    local_signal = flux[keep]
    span = float(np.ptp(local_signal))
    dt = float(np.median(np.diff(time_h)))
    peaks, _ = find_peaks(
        local_signal,
        prominence=max(0.05 * span, 1e-9),
        distance=max(1, int(round(0.55 * expected_period_h / dt))),
    )
    peak_times = time_h[indices[peaks]]
    if len(peak_times) < audited_cycles + 2:
        raise RuntimeError(
            f"Only {len(peak_times)} post-transient input peaks were found; "
            "increase the duration or inspect the A trajectory."
        )
    period_h = float(np.mean(np.diff(peak_times[-(audited_cycles + 1):])))
    boundaries = 0.5 * (peak_times[:-1] + peak_times[1:])
    return boundaries[-(audited_cycles + 1):], period_h


def a_cycle_metrics(time_h: np.ndarray, flux: np.ndarray, boundaries: np.ndarray) -> list[dict]:
    rows = []
    dt = float(np.median(np.diff(time_h)))
    for index, (left, right) in enumerate(zip(boundaries[:-1], boundaries[1:])):
        choose = (time_h >= left) & (time_h < right)
        t = time_h[choose]
        y = flux[choose]
        if len(t) < 3:
            continue
        low = float(np.min(y))
        high = float(np.max(y))
        half_level = low + 0.5 * (high - low)
        rows.append(
            {
                "cycle_index": index,
                "cycle_left_h": float(left),
                "cycle_right_h": float(right),
                "flux_peak_uM_h": high,
                "flux_trough_uM_h": low,
                "flux_FWHM_h": float(np.sum(y >= half_level) * dt),
                "flux_dose_uM_per_cycle": float(trapezoid(y, x=t)),
            }
        )
    return rows


def crossing_events(time: np.ndarray, lr_fraction: np.ndarray, left: float, right: float):
    local = np.flatnonzero((time >= left) & (time < right))
    events = []
    for first, second in zip(local[:-1], local[1:]):
        before = lr_fraction[first] - DNA_MIDPOINT
        after = lr_fraction[second] - DNA_MIDPOINT
        if before * after < 0.0:
            fraction = abs(before) / (abs(before) + abs(after))
            event_time = time[first] + fraction * (time[second] - time[first])
            events.append((event_time, "PB_to_LR" if after > before else "LR_to_PB"))
    return events


def sample_at(time: np.ndarray, signal: np.ndarray, when: float) -> float:
    return float(np.interp(float(when), time, signal))


def b_cycle_metrics(
    time: np.ndarray,
    lr_fraction: np.ndarray,
    pb_fraction: np.ndarray,
    flux_spline: CubicSpline,
    boundaries: np.ndarray,
) -> list[dict]:
    rows = []
    for index, (left, right) in enumerate(zip(boundaries[:-1], boundaries[1:])):
        events = crossing_events(time, lr_fraction, left, right)
        start_lr = sample_at(time, lr_fraction, left)
        end_lr = sample_at(time, lr_fraction, right)
        direction = events[0][1] if len(events) == 1 else ""
        endpoint_complete = bool(
            (direction == "PB_to_LR" and end_lr >= COMPLETE_STATE_THRESHOLD)
            or (direction == "LR_to_PB" and end_lr <= 1.0 - COMPLETE_STATE_THRESHOLD)
        )
        integration_grid = np.arange(left, right + 0.002, 0.002)
        dose = float(trapezoid(np.maximum(0.0, flux_spline(integration_grid)), x=integration_grid))
        rows.append(
            {
                "cycle_index": index,
                "cycle_left_h": float(left),
                "cycle_right_h": float(right),
                "start_LR_fraction": start_lr,
                "end_LR_fraction": end_lr,
                "total_crossings": len(events),
                "single_crossing_direction": direction,
                "endpoint_complete": endpoint_complete,
                "one_crossing_and_complete": bool(len(events) == 1 and endpoint_complete),
                "cycle_input_dose_uM": dose,
            }
        )
    return rows


def alternation_fidelity(rows: list[dict]) -> float:
    eligible = []
    for previous, current in zip(rows[:-1], rows[1:]):
        if previous["one_crossing_and_complete"] and current["one_crossing_and_complete"]:
            eligible.append(
                previous["single_crossing_direction"] != current["single_crossing_direction"]
            )
    return float(np.mean(eligible)) if eligible else float("nan")


def solve_b(zc, flux_spline: CubicSpline, duration_h: float, tag_rate: float,
            Td_min: float, initial_name: str):
    parameters = zc.default_params()
    parameters.update(
        K_rep=0.0186,
        n_rep=3.4,
        krep_tsl=15.0,
        krdf_tsl=200.0,
        k_tag_int=float(tag_rate),
        k_dil=math.log(2.0) / (Td_min / 60.0),
    )
    constants = zc._rate_constants(parameters)
    if initial_name == "PB":
        initial = zc.y0_PB_ss(parameters)
    elif initial_name == "LR":
        initial = zc.y0_LR_ss(parameters)
    else:
        raise ValueError(initial_name)

    time = np.arange(0.0, duration_h + B_SAMPLE_STEP_H / 2.0, B_SAMPLE_STEP_H)

    def int_production(t):
        return max(0.0, float(flux_spline(float(t))))

    solution = solve_ivp(
        lambda t, y: zc.rhs(t, y, parameters, constants, int_production),
        (float(time[0]), float(time[-1])),
        initial,
        t_eval=time,
        method="LSODA",
        rtol=1e-7,
        atol=1e-9,
        max_step=B_MAX_STEP_H,
    )
    if not solution.success:
        raise RuntimeError("B solve failed: " + solution.message)
    states = solution.y.T
    lr = zc.LR_total(states) / parameters["Dtot"]
    pb = zc.PB_total(states) / parameters["Dtot"]
    return time, states, lr, pb, parameters


def summarise_initial_state(cycles: list[dict], conservation: float) -> dict:
    one_complete = float(np.mean([row["one_crossing_and_complete"] for row in cycles]))
    zero = float(np.mean([row["total_crossings"] == 0 for row in cycles]))
    multiple = float(np.mean([row["total_crossings"] > 1 for row in cycles]))
    alternation = alternation_fidelity(cycles)
    passed = bool(
        one_complete >= SUCCESS_TARGET
        and np.isfinite(alternation)
        and alternation >= ALTERNATION_TARGET
        and conservation <= DNA_CONSERVATION_TOLERANCE
    )
    return {
        "one_crossing_complete_fraction": one_complete,
        "zero_crossing_fraction": zero,
        "multiple_crossing_fraction": multiple,
        "alternation_fidelity": alternation,
        "DNA_conservation_residual": conservation,
        "pass_initial_state": passed,
    }


def make_figures(condition_summary: pd.DataFrame, tag_summary: pd.DataFrame,
                 examples: dict, figure_dir: Path) -> None:
    figure_dir.mkdir(parents=True, exist_ok=True)
    targets = np.sort(tag_summary["target_period_generations"].unique())
    tags = np.sort(tag_summary["Int_tag_rate_h_inv"].unique())
    matrix = tag_summary.pivot(
        index="target_period_generations",
        columns="Int_tag_rate_h_inv",
        values="worst_one_crossing_complete_fraction",
    ).reindex(index=targets, columns=tags)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2), constrained_layout=True)
    image = axes[0].imshow(matrix.to_numpy(float), origin="lower", aspect="auto", vmin=0, vmax=1)
    axes[0].set_xticks(range(len(tags)), [f"{tag:g}" for tag in tags])
    axes[0].set_yticks(range(len(targets)), [f"{target:g}" for target in targets])
    axes[0].set_xlabel("Int tag rate (h^-1)")
    axes[0].set_ylabel("requested period (generations)")
    axes[0].set_title("Worst PB/LR one-crossing completion")
    for row in range(matrix.shape[0]):
        for column in range(matrix.shape[1]):
            value = matrix.iat[row, column]
            axes[0].text(column, row, f"{value:.2f}", ha="center", va="center",
                         color="white" if value < 0.55 else "black")
    fig.colorbar(image, ax=axes[0], shrink=0.85, label="fraction")

    for tag in tags:
        subset = tag_summary[tag_summary["Int_tag_rate_h_inv"].eq(tag)].sort_values(
            "target_period_generations"
        )
        axes[1].plot(
            subset["target_period_generations"],
            subset["both_initial_states_B_pass"].astype(int),
            "o-",
            label=f"tag={tag:g}/h",
        )
    axes[1].set(
        xlabel="requested period (generations)",
        ylabel="B interface certified",
        yticks=(0, 1),
        ylim=(-0.05, 1.05),
        title="Tag-specific deterministic certification",
    )
    axes[1].legend()
    fig.savefig(figure_dir / "Period_Knob_AB_certification_summary.png", dpi=200)
    plt.close(fig)

    if examples:
        ordered = sorted(examples.items(), key=lambda item: item[0][0])
        fig, axes = plt.subplots(len(ordered), 2, figsize=(13, 3.0 * len(ordered)), squeeze=False)
        for row_index, ((target, tag), data) in enumerate(ordered):
            time, flux, lr, pb = data
            view = time >= max(20.0, time[-1] - 45.0)
            axes[row_index, 0].plot(time[view], flux[view], color="magenta")
            axes[row_index, 0].set_ylabel(f"target={target:g}\ntag={tag:g}/h\nflux (uM/h)")
            axes[row_index, 1].plot(time[view], lr[view], "r", label="LR total")
            axes[row_index, 1].plot(time[view], pb[view], "b", label="PB total")
            axes[row_index, 1].axhline(DNA_MIDPOINT, color="0.4", ls="--")
            axes[row_index, 1].set_ylabel("DNA fraction")
        axes[0, 1].legend()
        axes[-1, 0].set_xlabel("time (h)")
        axes[-1, 1].set_xlabel("time (h)")
        fig.suptitle("Period-knob A input and B response (PB initial state)")
        fig.tight_layout()
        fig.savefig(figure_dir / "Period_Knob_AB_representative_trajectories.png", dpi=180)
        plt.close(fig)


def main() -> None:
    args = parse_args()
    mode = MODE_SETTINGS[args.mode]
    tag_rates = parse_float_list(args.tag_rates, mode["tag_rates_h_inv"])
    representative_csv = locate_representative_csv(args.representative_csv)
    representatives = load_representatives(representative_csv)
    required = {
        "target_period_generations",
        "period_generations",
        "period_error_generations",
        "within_period_tolerance",
        "design_RBS_scale",
        "design_oscillator_mRNA_half_life_min",
    }
    missing = required.difference(representatives.columns)
    if missing:
        raise KeyError(f"Representative table lacks required columns: {sorted(missing)}")
    representatives = select_smoke_panel(representatives, mode["maximum_representatives"])
    zc, zhao_path = import_zhao_core(args.b_model_dir)

    transient_h = float(mode["transient_h"])
    audited_cycles = int(mode["audited_cycles"])
    # A complete audit cycle is bounded by midpoints between input peaks.  It
    # therefore needs two more detected peaks than the requested number of
    # cycles, plus allowance for the phase between the transient cutoff and the
    # first subsequent peak.  Compute the duration from the slowest selected
    # period instead of relying on a fixed smoke/full duration.
    maximum_expected_period_h = float(
        pd.to_numeric(representatives["period_generations"], errors="raise").max()
        * 50.0
        / 60.0
    )
    minimum_duration_h = (
        transient_h + (audited_cycles + 3.0) * maximum_expected_period_h
    )
    duration_h = max(float(mode["duration_h"]), minimum_duration_h)
    if duration_h > float(mode["duration_h"]) + 1e-9:
        print(
            f"Auto-extended duration from {mode['duration_h']:.1f} h to "
            f"{duration_h:.1f} h so the slowest period has enough complete cycles."
        )
    base = core.make_reference_parameters()
    reference_peak = reference_demand_peak(base, duration_h)

    condition_rows = []
    cycle_rows = []
    a_rows = []
    examples = {}
    total_b_solves = len(representatives) * len(tag_rates) * len(DNA_INITIAL_STATES)
    b_counter = 0

    for candidate_number, (_, design) in enumerate(representatives.iterrows(), start=1):
        target = float(design["target_period_generations"])
        expected_period_h = float(design["period_generations"]) * 50.0 / 60.0
        parameters, physical = parameters_from_design(design, peak_load_fraction=0.0)
        print(
            f"[A {candidate_number:02d}/{len(representatives):02d}] target={target:g} gen | "
            f"RBS={design['design_RBS_scale']:g}, "
            f"mRNA t1/2={design['design_oscillator_mRNA_half_life_min']:g} min"
        )
        time_h, a_solution, tet_free, resource, flux = simulate_a(
            parameters, reference_peak, duration_h
        )
        boundaries, realised_period_h = cycle_boundaries(
            time_h, flux, transient_h, audited_cycles, expected_period_h
        )
        a_cycles = a_cycle_metrics(time_h, flux, boundaries)
        realised_period_generations = 60.0 * realised_period_h / parameters.Td_min
        flux_spline = CubicSpline(time_h, flux, bc_type="natural")

        for row in a_cycles:
            a_rows.append(
                {
                    "target_period_generations": target,
                    "design_RBS_scale": float(design["design_RBS_scale"]),
                    "design_oscillator_mRNA_half_life_min": float(
                        design["design_oscillator_mRNA_half_life_min"]
                    ),
                    **row,
                }
            )

        for tag_rate in tag_rates:
            for initial_name in DNA_INITIAL_STATES:
                b_counter += 1
                print(
                    f"    [B {b_counter:02d}/{total_b_solves:02d}] "
                    f"tag={tag_rate:g}/h, initial={initial_name}"
                )
                b_time, states, lr, pb, b_parameters = solve_b(
                    zc, flux_spline, duration_h, tag_rate, parameters.Td_min, initial_name
                )
                conservation = float(np.max(np.abs(lr + pb - 1.0)))
                cycles = b_cycle_metrics(b_time, lr, pb, flux_spline, boundaries)
                summary = summarise_initial_state(cycles, conservation)
                condition_rows.append(
                    {
                        "target_period_generations": target,
                        "realised_A_period_generations": realised_period_generations,
                        "period_error_generations": abs(realised_period_generations - target),
                        "source_within_period_tolerance": bool(design["within_period_tolerance"]),
                        "design_RBS_scale": float(design["design_RBS_scale"]),
                        "design_oscillator_mRNA_half_life_min": float(
                            design["design_oscillator_mRNA_half_life_min"]
                        ),
                        "Int_tag_rate_h_inv": tag_rate,
                        "Int_total_half_life_min": (
                            60.0 * math.log(2.0) / (b_parameters["k_dil"] + tag_rate)
                        ),
                        "initial_DNA_state": initial_name,
                        "A_flux_peak_uM_h_median": float(
                            np.median([row["flux_peak_uM_h"] for row in a_cycles])
                        ),
                        "A_flux_FWHM_h_median": float(
                            np.median([row["flux_FWHM_h"] for row in a_cycles])
                        ),
                        "A_flux_dose_uM_cycle_median": float(
                            np.median([row["flux_dose_uM_per_cycle"] for row in a_cycles])
                        ),
                        "A_resource_min": float(np.min(resource[time_h >= transient_h])),
                        "audited_cycles": len(cycles),
                        **summary,
                    }
                )
                for cycle in cycles:
                    cycle_rows.append(
                        {
                            "target_period_generations": target,
                            "design_RBS_scale": float(design["design_RBS_scale"]),
                            "design_oscillator_mRNA_half_life_min": float(
                                design["design_oscillator_mRNA_half_life_min"]
                            ),
                            "Int_tag_rate_h_inv": tag_rate,
                            "initial_DNA_state": initial_name,
                            **cycle,
                        }
                    )
                if initial_name == "PB":
                    current = examples.get((target, tag_rate))
                    candidate_quality = summary["one_crossing_complete_fraction"]
                    if current is None or candidate_quality > current[-1]:
                        examples[(target, tag_rate)] = (
                            b_time,
                            np.maximum(0.0, flux_spline(b_time)),
                            lr,
                            pb,
                            candidate_quality,
                        )

    condition_summary = pd.DataFrame(condition_rows)
    tag_rows = []
    merge_rows = []
    for _, design in representatives.iterrows():
        target = float(design["target_period_generations"])
        candidate_tag_rows = []
        for tag_rate in tag_rates:
            subset = condition_summary.loc[
                condition_summary["target_period_generations"].eq(target)
                & condition_summary["Int_tag_rate_h_inv"].eq(tag_rate)
            ]
            both_initials = bool(
                len(subset) == len(DNA_INITIAL_STATES)
                and subset["pass_initial_state"].all()
            )
            row = {
                "target_period_generations": target,
                "realised_A_period_generations": float(
                    subset["realised_A_period_generations"].iloc[0]
                ),
                "design_RBS_scale": float(design["design_RBS_scale"]),
                "design_oscillator_mRNA_half_life_min": float(
                    design["design_oscillator_mRNA_half_life_min"]
                ),
                "Int_tag_rate_h_inv": tag_rate,
                "both_initial_states_B_pass": both_initials,
                "worst_one_crossing_complete_fraction": float(
                    subset["one_crossing_complete_fraction"].min()
                ),
                "worst_alternation_fidelity": float(
                    subset["alternation_fidelity"].min(skipna=True)
                ),
                "max_DNA_conservation_residual": float(
                    subset["DNA_conservation_residual"].max()
                ),
            }
            tag_rows.append(row)
            candidate_tag_rows.append(row)

        passing_tags = [
            row["Int_tag_rate_h_inv"]
            for row in candidate_tag_rows
            if row["both_initial_states_B_pass"]
        ]
        source_period_ok = bool(design["within_period_tolerance"])
        b_pass = bool(passing_tags)
        merge_rows.append(
            {
                "target_period_generations": target,
                "design_RBS_scale": float(design["design_RBS_scale"]),
                "design_oscillator_mRNA_half_life_min": float(
                    design["design_oscillator_mRNA_half_life_min"]
                ),
                "B_pass": b_pass,
                "B_passing_Int_tag_rates_h_inv": ";".join(f"{value:g}" for value in passing_tags),
                "source_period_tolerance_pass": source_period_ok,
                "complete_deterministic_knob_certification": bool(b_pass and source_period_ok),
                "best_worst_one_crossing_complete_fraction": max(
                    row["worst_one_crossing_complete_fraction"] for row in candidate_tag_rows
                ),
            }
        )

    tag_summary = pd.DataFrame(tag_rows)
    merge_table = pd.DataFrame(merge_rows)
    cycle_table = pd.DataFrame(cycle_rows)
    a_cycle_table = pd.DataFrame(a_rows)

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = OUTPUT_ROOT / f"Robust_Continuous_Period_Knob_AB_{args.mode}_{stamp}"
    figure_dir = FIGURE_ROOT / f"Robust_Continuous_Period_Knob_AB_{args.mode}_{stamp}"
    output_dir.mkdir(parents=True, exist_ok=False)
    figure_dir.mkdir(parents=True, exist_ok=False)
    condition_summary.to_csv(output_dir / "AB_initial_state_condition_summary.csv", index=False)
    tag_summary.to_csv(output_dir / "AB_tag_specific_certification.csv", index=False)
    merge_table.to_csv(output_dir / "period_knob_B_validation_for_merge.csv", index=False)
    cycle_table.to_csv(output_dir / "B_cycle_metrics.csv", index=False)
    a_cycle_table.to_csv(output_dir / "A_flux_cycle_metrics.csv", index=False)

    plot_examples = {
        key: value[:-1]
        for key, value in examples.items()
        if key[1] == tag_rates[0]
    }
    make_figures(condition_summary, tag_summary, plot_examples, figure_dir)

    metadata = {
        "version": VERSION,
        "mode": args.mode,
        "representative_csv": str(representative_csv),
        "zhao_core": str(zhao_path),
        "tag_rates_h_inv": tag_rates,
        "DNA_initial_states": DNA_INITIAL_STATES,
        "duration_h": duration_h,
        "transient_h": transient_h,
        "audited_cycles": audited_cycles,
        "success_target": SUCCESS_TARGET,
        "alternation_target": ALTERNATION_TARGET,
        "DNA_conservation_tolerance": DNA_CONSERVATION_TOLERANCE,
        "interface_equation": "int_production = rho * beta31 * m31 converted to uM/h",
        "shared_PLtetO1_constraint": True,
        "resource_load_fraction": 0.0,
        "interpretation": [
            "This is deterministic A-to-B interface certification, not stochastic fidelity.",
            "A period-knob setting passes B only if PB and LR initial states both pass.",
            "The five A-clock settings were selected by shared-scenario robust ordering.",
            "Resource-load robustness is a later test and is not mixed into this first certification.",
        ],
    }
    (output_dir / "scope_and_certification_rules.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    print(f"\nCompleted: {VERSION}")
    print(f"Mode: {args.mode}")
    print(f"Raw tables: {output_dir}")
    print(f"Figures: {figure_dir}")
    print("\nPeriod-knob interface certification:")
    for row in merge_rows:
        print(
            f"target={row['target_period_generations']:.2f} gen | "
            f"B_pass={row['B_pass']} | "
            f"passing tags={row['B_passing_Int_tag_rates_h_inv'] or 'none'} | "
            f"period tolerance={row['source_period_tolerance_pass']} | "
            f"complete={row['complete_deterministic_knob_certification']}"
        )


if __name__ == "__main__":
    main()
