r"""Stochastic end-to-end A-to-B certification of robust K4 and K5.

The two slow settings K4 and K5 from the confirmed robust continuous clock
library have already passed deterministic A-clock and A-to-B checks.  They are
the only settings in that five-point library above the independently mapped B
recovery boundary (10.5 generations at Td=50 min, tag=12 h^-1).

1. deterministic A-module period and waveform constraints;
2. deterministic B-module interface certification; and
3. stochastic A-module period/order certification.

This program performs the missing final test.  For every independent stochastic
A trajectory it sends the *instantaneous phiC31 translation flux* directly to
the Zhao B-module ODE, then audits both PB and LR initial DNA states.

Primary endpoints
-----------------
* fraction of cycles with exactly one midpoint crossing and a complete endpoint;
* alternation fidelity across consecutive cycles;
* probability that seven consecutive cycles are all correct;
* DNA-state conservation residual;
* stochastic A period and C31 input summaries for traceability.

Scope boundary
--------------
The A module is stochastic (non-negative binomial tau leap, inherited from
``Stochastic_Period_Knob_Library_Certification.py``); the B chemistry is the
deterministic Zhao ODE driven by each stochastic A input.  The resulting
probabilities therefore quantify propagation of A-side intrinsic noise through
the current B model.  They do not include stochastic recombination chemistry,
cell-to-cell B-parameter variation, or population selection.

No model parameter is fitted here.

Run from the project root::

    python .\code\Robust_K4_K5_Stochastic_AB_Fidelity_Certification.py --mode smoke --tau-step-min 0.05
    python .\code\Robust_K4_K5_Stochastic_AB_Fidelity_Certification.py --mode audit --tau-step-min 0.05
    python .\code\Robust_K4_K5_Stochastic_AB_Fidelity_Certification.py --mode confirmation --tau-step-min 0.05

The default physical B design uses k_tag_int = 12 h^-1, the tag under which K4
and K5 passed deterministic interface certification. Extra tag rates may be
supplied with ``--tag-rates 8,12`` as a robustness diagnostic.
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
from scipy.integrate import trapezoid

import Stochastic_Period_Knob_Library_Certification as stochastic_a
import Robust_Continuous_Period_Knob_AB_Certification as robust_cert


VERSION = "Robust_K4_K5_Stochastic_AB_Fidelity_Certification_v3_20260823"
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
OUTPUT_ROOT = PROJECT_ROOT / "outputs"
FIGURE_ROOT = PROJECT_ROOT / "figures"

MOLECULES_PER_UM_PER_FL = 602.214076
CELL_VOLUME_FL = 1.0
DNA_MIDPOINT = 0.5
COMPLETE_STATE_THRESHOLD = 0.95
DNA_CONSERVATION_TOLERANCE = 1e-6
B_SAMPLE_STEP_H = 1.0 / 60.0
B_MAX_STEP_H = 0.02
B_MINIMUM_RECOVERY_PERIOD_GENERATIONS = 10.5


@dataclass(frozen=True)
class ABRunSettings:
    trajectories_per_setting: int
    duration_generations: float
    transient_generations: float
    audited_cycles: int
    tau_step_min: float
    sample_step_min: float
    bootstrap_replicates: int
    random_seed: int = 20260819
    minimum_cycle_success_fraction: float = 0.95
    minimum_seven_cycle_success_fraction: float = 0.70
    minimum_alternation_fidelity: float = 0.95


MODE_SETTINGS = {
    # Workflow checks only.  Do not report smoke probabilities as confirmation.
    "smoke": ABRunSettings(2, 135.0, 20.0, 7, 0.05, 1.0, 500),
    "audit": ABRunSettings(12, 150.0, 25.0, 7, 0.05, 1.0, 2000),
    "confirmation": ABRunSettings(50, 165.0, 25.0, 7, 0.05, 1.0, 5000),
}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=tuple(MODE_SETTINGS), default="smoke")
    parser.add_argument(
        "--tau-step-min",
        type=float,
        default=None,
        help=(
            "Tau-leap integration step in minutes. The convergence-audited "
            "production default is 0.05 min."
        ),
    )
    parser.add_argument(
        "--library-csv",
        type=Path,
        default=None,
        help=(
            "Optional ranked_five_setting_libraries.csv from the confirmed "
            "robust redesign. By default the newest confirmed file is used."
        ),
    )
    parser.add_argument("--b-model-dir", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument(
        "--tag-rates",
        default="12",
        help="Comma-separated Int tag rates in h^-1. Default: 12.",
    )
    parser.add_argument(
        "--fixed-plasmids",
        action="store_true",
        help="Diagnostic only: remove plasmid-copy fluctuations from stochastic A.",
    )
    return parser.parse_args()


def parse_float_list(text: str) -> tuple[float, ...]:
    values = tuple(float(item.strip()) for item in text.split(",") if item.strip())
    if not values or any(value < 0.0 for value in values):
        raise ValueError("--tag-rates must contain non-negative numbers")
    return values


def load_k4_k5_library(explicit_path: Path | None) -> tuple[Path, pd.DataFrame]:
    """Convert robust K4/K5 definitions to the stochastic-A library schema."""
    path = robust_cert.locate_representative_csv(explicit_path)
    robust = robust_cert.load_representatives(path)
    robust = robust.sort_values("target_period_generations").reset_index(drop=True)
    if len(robust) != 5:
        raise ValueError(
            f"Expected a five-setting robust library, found {len(robust)} rows in {path}."
        )

    selected = robust.iloc[[3, 4]].copy().reset_index(drop=True)
    rows = []
    for index, row in selected.iterrows():
        nominal = float(row["period_generations"])
        rows.append(
            {
                "knob_setting_id": f"Robust_K{index + 4}",
                "source_setting_id": str(row["setting_id"]),
                "target_period_generations": float(row["target_period_generations"]),
                "realised_period_generations": nominal,
                "synchronous_oscillator_RBS_scale": float(row["design_RBS_scale"]),
                "synchronous_oscillator_mRNA_half_life_min": float(
                    row["design_oscillator_mRNA_half_life_min"]
                ),
                # Keep the same C31 output branch used by the deterministic
                # robust design.  Shared_PLtetO1_Period_Knob_Design_Map freezes
                # the C31 RBS/translation scale at 0.45; using 1.0 here would
                # inflate stochastic C31 flux by 1/0.45 = 2.22-fold before any
                # stochastic effect is considered.
                "C31_RBS_scale": 0.45,
                "C31_mRNA_total_half_life_min": 2.0,
                "sponge_plasmid_copies": 40.0,
                "sponge_accessible_fraction": 1.0,
                "deterministic_A_pass": True,
                "deterministic_B_pass": True,
                "B_minimum_recovery_period_generations": (
                    B_MINIMUM_RECOVERY_PERIOD_GENERATIONS
                ),
                "deterministic_recovery_margin_generations": (
                    nominal - B_MINIMUM_RECOVERY_PERIOD_GENERATIONS
                ),
            }
        )
    return path, pd.DataFrame(rows)


def validate_k4_k5_library(frame: pd.DataFrame) -> pd.DataFrame:
    """Validate the two robust settings without the legacy T10--T13 gate."""
    required = {
        "knob_setting_id",
        "target_period_generations",
        "realised_period_generations",
        "synchronous_oscillator_RBS_scale",
        "synchronous_oscillator_mRNA_half_life_min",
        "C31_RBS_scale",
        "C31_mRNA_total_half_life_min",
        "sponge_plasmid_copies",
        "sponge_accessible_fraction",
        "deterministic_A_pass",
        "deterministic_B_pass",
    }
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"K4/K5 library is missing columns: {missing}")
    selected = frame.copy().sort_values("target_period_generations").reset_index(drop=True)
    if tuple(selected["knob_setting_id"]) != ("Robust_K4", "Robust_K5"):
        raise ValueError(
            "Expected exactly Robust_K4 and Robust_K5; found "
            f"{tuple(selected['knob_setting_id'])}."
        )
    if not (
        selected["deterministic_A_pass"].astype(bool).all()
        and selected["deterministic_B_pass"].astype(bool).all()
    ):
        raise ValueError("K4 and K5 must pass deterministic A and B certification first.")
    return selected


def import_zhao_core(explicit_directory: Path | None):
    import importlib.util

    candidates = []
    if explicit_directory is not None:
        candidates.append(explicit_directory.expanduser().resolve() / "zhao_core.py")
    candidates.extend(
        [
            SCRIPT_DIR / "zhao_core.py",
            PROJECT_ROOT / "zhao_core.py",
            Path(
                r"E:\iGEM\建模工作\暑期\工作汇总\SummerWeek3 - 副本\model\zhao_core.py"
            ),
        ]
    )
    for path in candidates:
        if not path.exists():
            continue
        spec = importlib.util.spec_from_file_location("zhao_core_stochastic_ab", path)
        if spec is None or spec.loader is None:
            continue
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        required = (
            "default_params", "_rate_constants", "rhs", "y0_PB_ss", "y0_LR_ss",
            "LR_total", "PB_total",
        )
        missing = [name for name in required if not hasattr(module, name)]
        if missing:
            raise ImportError(f"{path} lacks required Zhao API: {missing}")
        return module, path
    raise FileNotFoundError(
        "zhao_core.py was not found. Copy it into this code folder or use "
        "--b-model-dir with the folder that contains it."
    )


def make_a_settings(settings: ABRunSettings) -> stochastic_a.RunSettings:
    return stochastic_a.RunSettings(
        trajectories_per_setting=settings.trajectories_per_setting,
        duration_generations=settings.duration_generations,
        transient_generations=settings.transient_generations,
        tau_step_min=settings.tau_step_min,
        sample_step_min=settings.sample_step_min,
        bootstrap_replicates=settings.bootstrap_replicates,
        random_seed=settings.random_seed,
        minimum_cycles=max(4, settings.audited_cycles),
    )


def flux_to_uM_h(flux_molecules_per_min: np.ndarray) -> np.ndarray:
    return (
        60.0
        * np.maximum(np.asarray(flux_molecules_per_min, dtype=float), 0.0)
        / (MOLECULES_PER_UM_PER_FL * CELL_VOLUME_FL)
    )


def choose_cycle_boundaries(
    phase_crossings_min: np.ndarray,
    deterministic_period_min: float,
    audited_cycles: int,
) -> np.ndarray:
    crossings = np.asarray(phase_crossings_min, dtype=float)
    crossings = crossings[np.isfinite(crossings)]
    if crossings.size < audited_cycles + 1:
        return np.array([], dtype=float)
    # Select the latest contiguous window whose periods remain within the same
    # broad validity limits used by the stochastic A certification.
    for stop in range(crossings.size, audited_cycles, -1):
        candidate = crossings[stop - audited_cycles - 1 : stop]
        intervals = np.diff(candidate)
        if np.all(
            (intervals >= 0.45 * deterministic_period_min)
            & (intervals <= 1.80 * deterministic_period_min)
        ):
            return candidate
    return np.array([], dtype=float)


def solve_b(
    zc,
    input_time_h: np.ndarray,
    input_flux_uM_h: np.ndarray,
    tag_rate_h_inv: float,
    Td_min: float,
    initial_name: str,
):
    parameters = zc.default_params()
    parameters.update(
        K_rep=0.0186,
        n_rep=3.4,
        krep_tsl=15.0,
        krdf_tsl=200.0,
        k_tag_int=float(tag_rate_h_inv),
        k_dil=math.log(2.0) / (Td_min / 60.0),
    )
    constants = zc._rate_constants(parameters)
    if initial_name == "PB":
        initial = zc.y0_PB_ss(parameters)
    elif initial_name == "LR":
        initial = zc.y0_LR_ss(parameters)
    else:
        raise ValueError(initial_name)

    duration_h = float(input_time_h[-1])
    output_time_h = np.arange(0.0, duration_h + B_SAMPLE_STEP_H / 2.0, B_SAMPLE_STEP_H)

    def int_production(t):
        return max(0.0, float(np.interp(t, input_time_h, input_flux_uM_h)))

    solution = solve_ivp(
        lambda t, y: zc.rhs(t, y, parameters, constants, int_production),
        (0.0, duration_h),
        initial,
        t_eval=output_time_h,
        method="LSODA",
        rtol=1e-7,
        atol=1e-9,
        max_step=B_MAX_STEP_H,
    )
    if not solution.success:
        raise RuntimeError("B solve failed: " + solution.message)
    states = solution.y.T
    lr = np.asarray(zc.LR_total(states), dtype=float) / parameters["Dtot"]
    pb = np.asarray(zc.PB_total(states), dtype=float) / parameters["Dtot"]
    return output_time_h, lr, pb, states, parameters


def state_label(lr_fraction: float) -> str:
    if lr_fraction >= COMPLETE_STATE_THRESHOLD:
        return "LR"
    if lr_fraction <= 1.0 - COMPLETE_STATE_THRESHOLD:
        return "PB"
    return "incomplete"


def midpoint_events(time_h, lr_fraction, left_h, right_h):
    select = np.flatnonzero((time_h >= left_h) & (time_h <= right_h))
    events = []
    for first, second in zip(select[:-1], select[1:]):
        before = lr_fraction[first] - DNA_MIDPOINT
        after = lr_fraction[second] - DNA_MIDPOINT
        if before * after < 0.0:
            fraction = abs(before) / (abs(before) + abs(after))
            crossing_time = time_h[first] + fraction * (time_h[second] - time_h[first])
            direction = "PB_to_LR" if after > before else "LR_to_PB"
            events.append((float(crossing_time), direction))
    return events


def audit_cycles(
    b_time_h: np.ndarray,
    lr_fraction: np.ndarray,
    input_time_h: np.ndarray,
    input_flux_uM_h: np.ndarray,
    boundaries_min: np.ndarray,
):
    boundaries_h = np.asarray(boundaries_min, dtype=float) / 60.0
    initial_lr = float(np.interp(boundaries_h[0], b_time_h, lr_fraction))
    initial_label = state_label(initial_lr)
    expected_label = initial_label
    rows = []
    for cycle_index, (left_h, right_h) in enumerate(
        zip(boundaries_h[:-1], boundaries_h[1:])
    ):
        expected_label = (
            "LR" if expected_label == "PB" else "PB" if expected_label == "LR" else "incomplete"
        )
        events = midpoint_events(b_time_h, lr_fraction, left_h, right_h)
        start_lr = float(np.interp(left_h, b_time_h, lr_fraction))
        end_lr = float(np.interp(right_h, b_time_h, lr_fraction))
        end_label = state_label(end_lr)
        input_select = (input_time_h >= left_h) & (input_time_h <= right_h)
        local_t = input_time_h[input_select]
        local_flux = input_flux_uM_h[input_select]
        dose = float(trapezoid(local_flux, x=local_t)) if local_t.size >= 2 else float("nan")
        peak = float(np.max(local_flux)) if local_flux.size else float("nan")
        correct = bool(
            initial_label != "incomplete"
            and len(events) == 1
            and end_label == expected_label
        )
        rows.append(
            {
                "cycle_index": cycle_index,
                "cycle_left_h": float(left_h),
                "cycle_right_h": float(right_h),
                "cycle_duration_generations": float((right_h - left_h) * 60.0 / 50.0),
                "start_LR_fraction": start_lr,
                "end_LR_fraction": end_lr,
                "start_state": state_label(start_lr),
                "expected_end_state": expected_label,
                "observed_end_state": end_label,
                "total_crossings": len(events),
                "single_crossing_direction": events[0][1] if len(events) == 1 else "",
                "endpoint_complete": end_label != "incomplete",
                "cycle_correct": correct,
                "input_peak_uM_h": peak,
                "input_dose_uM": dose,
            }
        )
    return rows, initial_label


def alternation_fidelity(cycles: list[dict]) -> float:
    directions = [row["single_crossing_direction"] for row in cycles]
    eligible = []
    for previous, current in zip(directions[:-1], directions[1:]):
        if previous and current:
            eligible.append(previous != current)
    return float(np.mean(eligible)) if eligible else float("nan")


def bootstrap_interval(values, replicates, seed):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size < 2:
        return float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    estimates = np.empty(replicates, dtype=float)
    for index in range(replicates):
        estimates[index] = np.mean(rng.choice(values, size=values.size, replace=True))
    return tuple(float(value) for value in np.percentile(estimates, [2.5, 97.5]))


def summarize_conditions(trajectory_frame: pd.DataFrame, settings: ABRunSettings):
    rows = []
    group_columns = ["knob_setting_id", "target_period_generations", "Int_tag_rate_h_inv", "initial_DNA_state"]
    for group_key, subset in trajectory_frame.groupby(group_columns, sort=True):
        setting_id, target, tag_rate, initial_state = group_key
        cycle_values = subset["cycle_success_fraction"].to_numpy(float)
        seven_values = subset["seven_cycle_sequence_success"].astype(float).to_numpy()
        cycle_low, cycle_high = bootstrap_interval(
            cycle_values, settings.bootstrap_replicates,
            settings.random_seed + int(round(target * 100)) + int(round(tag_rate * 10)),
        )
        seven_low, seven_high = bootstrap_interval(
            seven_values, settings.bootstrap_replicates,
            settings.random_seed + 1_000_000 + int(round(target * 100)) + int(round(tag_rate * 10)),
        )
        cycle_mean = float(np.mean(cycle_values))
        seven_mean = float(np.mean(seven_values))
        alternation = float(np.nanmean(subset["alternation_fidelity"]))
        conservation = float(np.max(subset["DNA_conservation_residual"]))
        condition_pass = bool(
            cycle_mean >= settings.minimum_cycle_success_fraction
            and seven_mean >= settings.minimum_seven_cycle_success_fraction
            and np.isfinite(alternation)
            and alternation >= settings.minimum_alternation_fidelity
            and conservation <= DNA_CONSERVATION_TOLERANCE
        )
        rows.append(
            {
                "knob_setting_id": setting_id,
                "target_period_generations": float(target),
                "Int_tag_rate_h_inv": float(tag_rate),
                "initial_DNA_state": initial_state,
                "n_stochastic_A_trajectories": len(subset),
                "mean_cycle_success_fraction": cycle_mean,
                "cycle_success_CI_low": cycle_low,
                "cycle_success_CI_high": cycle_high,
                "seven_cycle_sequence_success_fraction": seven_mean,
                "seven_cycle_success_CI_low": seven_low,
                "seven_cycle_success_CI_high": seven_high,
                "mean_alternation_fidelity": alternation,
                "max_DNA_conservation_residual": conservation,
                "condition_pass": condition_pass,
            }
        )
    return pd.DataFrame(rows)


def summarize_settings(condition_frame: pd.DataFrame):
    rows = []
    for (setting_id, target), subset in condition_frame.groupby(
        ["knob_setting_id", "target_period_generations"], sort=True
    ):
        passing_tags = []
        for tag_rate, tag_subset in subset.groupby("Int_tag_rate_h_inv"):
            if set(tag_subset["initial_DNA_state"]) == {"PB", "LR"} and bool(
                tag_subset["condition_pass"].all()
            ):
                passing_tags.append(float(tag_rate))
        rows.append(
            {
                "knob_setting_id": setting_id,
                "target_period_generations": float(target),
                "worst_cycle_success_fraction": float(subset["mean_cycle_success_fraction"].min()),
                "worst_seven_cycle_success_fraction": float(
                    subset["seven_cycle_sequence_success_fraction"].min()
                ),
                "worst_alternation_fidelity": float(subset["mean_alternation_fidelity"].min()),
                "stochastic_AB_passing_tag_rates_h_inv": ";".join(f"{x:g}" for x in passing_tags),
                "stochastic_AB_interface_pass": bool(passing_tags),
            }
        )
    return pd.DataFrame(rows)


def make_figures(condition_frame, setting_frame, examples, figure_dir: Path):
    figure_dir.mkdir(parents=True, exist_ok=True)
    labels = []
    cycle_values = []
    seven_values = []
    for _, row in condition_frame.sort_values(
        ["target_period_generations", "Int_tag_rate_h_inv", "initial_DNA_state"]
    ).iterrows():
        labels.append(
            f"T{row['target_period_generations']:g}\ntag={row['Int_tag_rate_h_inv']:g}\n{row['initial_DNA_state']}"
        )
        cycle_values.append(row["mean_cycle_success_fraction"])
        seven_values.append(row["seven_cycle_sequence_success_fraction"])
    x = np.arange(len(labels))
    fig, axes = plt.subplots(2, 1, figsize=(max(11, 0.75 * len(labels)), 8), constrained_layout=True)
    axes[0].bar(x, cycle_values, color="#4C78A8")
    axes[0].axhline(0.95, color="black", ls="--", label="cycle target = 0.95")
    axes[0].set(ylabel="mean correct-cycle fraction", ylim=(0, 1.05), xticks=x, xticklabels=labels)
    axes[0].legend()
    axes[1].bar(x, seven_values, color="#F58518")
    axes[1].axhline(0.70, color="black", ls="--", label="7-cycle target = 0.70")
    axes[1].set(ylabel="seven-cycle sequence success", ylim=(0, 1.05), xticks=x, xticklabels=labels)
    axes[1].legend()
    fig.suptitle("Stochastic A-to-B period-knob fidelity certification")
    fig.savefig(figure_dir / "Stochastic_AB_fidelity_summary.png", dpi=200)
    plt.close(fig)

    if examples:
        ordered = sorted(examples.items(), key=lambda item: item[0])
        fig, axes = plt.subplots(len(ordered), 2, figsize=(13, 3.0 * len(ordered)), squeeze=False)
        for row_index, (key, data) in enumerate(ordered):
            target, initial_state = key
            time_h, flux, b_time, lr, pb, boundaries_h = data
            left = boundaries_h[0]
            right = boundaries_h[-1]
            a_view = (time_h >= left) & (time_h <= right)
            b_view = (b_time >= left) & (b_time <= right)
            axes[row_index, 0].plot(time_h[a_view], flux[a_view], color="magenta", lw=1)
            axes[row_index, 0].set_ylabel(f"T{target:g} {initial_state}\nflux (uM/h)")
            axes[row_index, 1].plot(b_time[b_view], lr[b_view], "r", label="LR total")
            axes[row_index, 1].plot(b_time[b_view], pb[b_view], "b", label="PB total")
            axes[row_index, 1].axhline(0.5, color="0.4", ls="--")
            for boundary in boundaries_h:
                axes[row_index, 0].axvline(boundary, color="0.7", lw=0.5)
                axes[row_index, 1].axvline(boundary, color="0.7", lw=0.5)
            axes[row_index, 1].set_ylim(-0.03, 1.03)
        axes[0, 1].legend()
        axes[-1, 0].set_xlabel("time (h)")
        axes[-1, 1].set_xlabel("time (h)")
        fig.suptitle("Representative stochastic A inputs and B DNA responses")
        fig.tight_layout()
        fig.savefig(figure_dir / "Stochastic_AB_representative_trajectories.png", dpi=180)
        plt.close(fig)


def main():
    args = parse_args()
    settings = MODE_SETTINGS[args.mode]
    if args.tau_step_min is not None:
        if not math.isfinite(args.tau_step_min) or args.tau_step_min <= 0.0:
            raise ValueError("--tau-step-min must be a finite positive number")
        settings = replace(settings, tau_step_min=float(args.tau_step_min))
    a_settings = make_a_settings(settings)
    tag_rates = parse_float_list(args.tag_rates)
    library_path, library = load_k4_k5_library(args.library_csv)
    library = validate_k4_k5_library(library)
    zc, zhao_path = import_zhao_core(args.b_model_dir)

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = (
        args.output_dir.expanduser().resolve()
        if args.output_dir is not None
        else OUTPUT_ROOT
        / (
            f"Robust_K4_K5_Stochastic_AB_{args.mode}_"
            f"tau{settings.tau_step_min:g}min_{stamp}"
        )
    )
    figure_dir = FIGURE_ROOT / output_dir.name
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    print(f"{VERSION}: mode={args.mode}")
    print(f"Certified library: {library_path}")
    print(f"B model: {zhao_path}")
    print(
        f"A trajectories={settings.trajectories_per_setting}/setting; "
        f"audited cycles={settings.audited_cycles}; tag rates={tag_rates}; "
        f"tau step={settings.tau_step_min:g} min"
    )

    trajectory_rows = []
    cycle_rows = []
    a_rows = []
    examples = {}
    total_a = len(library) * settings.trajectories_per_setting
    total_b = total_a * len(tag_rates) * 2
    a_counter = 0
    b_counter = 0

    for setting_index, (_, row) in enumerate(library.iterrows()):
        deterministic_period_min = float(row["realised_period_generations"]) * 50.0
        for trajectory_index in range(settings.trajectories_per_setting):
            a_counter += 1
            seed = settings.random_seed + setting_index * 100_000 + trajectory_index
            print(
                f"[A {a_counter:03d}/{total_a:03d}] {row['knob_setting_id']} "
                f"trajectory={trajectory_index + 1:02d}"
            )
            result = stochastic_a.simulate_trajectory(
                row, a_settings, seed, args.fixed_plasmids
            )
            analysis = stochastic_a.analyse_trajectory(row, result, a_settings)
            boundaries_min = choose_cycle_boundaries(
                analysis["phase_crossings_min"], deterministic_period_min, settings.audited_cycles
            )
            time_h = result["time_min"] / 60.0
            flux_uM_h = flux_to_uM_h(result["c31_flux_molecules_per_min"])
            a_rows.append(
                {
                    "knob_setting_id": row["knob_setting_id"],
                    "target_period_generations": float(row["target_period_generations"]),
                    "trajectory_index": trajectory_index,
                    "seed": seed,
                    "A_resolved": bool(analysis["resolved"]),
                    "A_phase_period_mean_generations": analysis["phase_period_mean_generations"],
                    "A_within_trajectory_period_cv": analysis[
                        "phase_period_within_trajectory_cv"
                    ],
                    "n_usable_boundaries": int(boundaries_min.size),
                    "A_flux_peak_uM_h": float(np.max(flux_uM_h)),
                    "A_flux_mean_uM_h": float(np.mean(flux_uM_h)),
                }
            )
            if boundaries_min.size != settings.audited_cycles + 1:
                print("    skipped B: insufficient valid stochastic phase boundaries")
                for tag_rate in tag_rates:
                    for initial_name in ("PB", "LR"):
                        trajectory_rows.append(
                            {
                                "knob_setting_id": row["knob_setting_id"],
                                "target_period_generations": float(row["target_period_generations"]),
                                "trajectory_index": trajectory_index,
                                "seed": seed,
                                "Int_tag_rate_h_inv": tag_rate,
                                "initial_DNA_state": initial_name,
                                "audited_cycles": 0,
                                "cycle_success_fraction": 0.0,
                                "alternation_fidelity": float("nan"),
                                "seven_cycle_sequence_success": False,
                                "DNA_conservation_residual": float("nan"),
                                "A_resolved": False,
                            }
                        )
                continue

            for tag_rate in tag_rates:
                for initial_name in ("PB", "LR"):
                    b_counter += 1
                    print(
                        f"    [B {b_counter:03d}/{total_b:03d}] "
                        f"tag={tag_rate:g}/h initial={initial_name}"
                    )
                    b_time, lr, pb, states, b_parameters = solve_b(
                        zc, time_h, flux_uM_h, tag_rate, 50.0, initial_name
                    )
                    conservation = float(np.max(np.abs(lr + pb - 1.0)))
                    cycles, audit_initial_label = audit_cycles(
                        b_time, lr, time_h, flux_uM_h, boundaries_min
                    )
                    successes = np.array([cycle["cycle_correct"] for cycle in cycles], dtype=bool)
                    cycle_success = float(np.mean(successes)) if successes.size else 0.0
                    seven_success = bool(
                        successes.size >= 7 and np.all(successes[:7])
                    )
                    alternation = alternation_fidelity(cycles)
                    trajectory_rows.append(
                        {
                            "knob_setting_id": row["knob_setting_id"],
                            "target_period_generations": float(row["target_period_generations"]),
                            "trajectory_index": trajectory_index,
                            "seed": seed,
                            "Int_tag_rate_h_inv": tag_rate,
                            "Int_total_half_life_min": (
                                60.0 * math.log(2.0) / (b_parameters["k_dil"] + tag_rate)
                            ),
                            "initial_DNA_state": initial_name,
                            "audit_initial_state": audit_initial_label,
                            "audited_cycles": len(cycles),
                            "cycle_success_fraction": cycle_success,
                            "zero_crossing_fraction": float(
                                np.mean([cycle["total_crossings"] == 0 for cycle in cycles])
                            ),
                            "multiple_crossing_fraction": float(
                                np.mean([cycle["total_crossings"] > 1 for cycle in cycles])
                            ),
                            "alternation_fidelity": alternation,
                            "seven_cycle_sequence_success": seven_success,
                            "DNA_conservation_residual": conservation,
                            "A_resolved": bool(analysis["resolved"]),
                            "A_period_generations": analysis["phase_period_mean_generations"],
                        }
                    )
                    for cycle in cycles:
                        cycle_rows.append(
                            {
                                "knob_setting_id": row["knob_setting_id"],
                                "target_period_generations": float(row["target_period_generations"]),
                                "trajectory_index": trajectory_index,
                                "seed": seed,
                                "Int_tag_rate_h_inv": tag_rate,
                                "initial_DNA_state": initial_name,
                                **cycle,
                            }
                        )
                    example_key = (float(row["target_period_generations"]), initial_name)
                    if trajectory_index == 0 and tag_rate == tag_rates[-1]:
                        examples[example_key] = (
                            time_h,
                            flux_uM_h,
                            b_time,
                            lr,
                            pb,
                            boundaries_min / 60.0,
                        )

    a_frame = pd.DataFrame(a_rows)
    trajectory_frame = pd.DataFrame(trajectory_rows)
    cycle_frame = pd.DataFrame(cycle_rows)
    condition_frame = summarize_conditions(trajectory_frame, settings)
    setting_frame = summarize_settings(condition_frame)

    merged = library.merge(
        setting_frame,
        on=["knob_setting_id", "target_period_generations"],
        how="left",
    )
    merged["complete_deterministic_and_stochastic_AB_pass"] = (
        merged["deterministic_A_pass"].astype(bool)
        & merged["deterministic_B_pass"].astype(bool)
        & merged["stochastic_AB_interface_pass"].fillna(False).astype(bool)
    )

    a_frame.to_csv(output_dir / "stochastic_A_trajectory_audit.csv", index=False)
    trajectory_frame.to_csv(output_dir / "stochastic_AB_trajectory_summary.csv", index=False)
    cycle_frame.to_csv(output_dir / "stochastic_AB_cycle_metrics.csv", index=False)
    condition_frame.to_csv(output_dir / "stochastic_AB_condition_summary.csv", index=False)
    setting_frame.to_csv(output_dir / "stochastic_AB_setting_certification.csv", index=False)
    merged.to_csv(output_dir / "period_knob_library_with_stochastic_AB_results.csv", index=False)

    metadata = {
        "version": VERSION,
        "mode": args.mode,
        "settings": asdict(settings),
        "library_csv": str(library_path),
        "zhao_core": str(zhao_path),
        "tag_rates_h_inv": tag_rates,
        "cell_volume_fL": CELL_VOLUME_FL,
        "fixed_plasmids_diagnostic": bool(args.fixed_plasmids),
        "B_minimum_recovery_period_generations": B_MINIMUM_RECOVERY_PERIOD_GENERATIONS,
        "scope": (
            "stochastic A intrinsic/plasmid noise propagated into deterministic B; "
            "B intrinsic noise and parameter heterogeneity are not included"
        ),
    }
    (output_dir / "stochastic_AB_scope_and_settings.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    make_figures(condition_frame, setting_frame, examples, figure_dir)

    print(f"\nCompleted: {VERSION}")
    print(f"Mode: {args.mode}")
    print(f"Raw tables: {output_dir}")
    print(f"Figures: {figure_dir}")
    print("\nStochastic A-to-B period-knob certification:")
    for _, result in setting_frame.sort_values("target_period_generations").iterrows():
        print(
            f"{result['knob_setting_id']} | target={result['target_period_generations']:.2f} gen | "
            f"worst cycle={result['worst_cycle_success_fraction']:.3f} | "
            f"worst 7-cycle={result['worst_seven_cycle_success_fraction']:.3f} | "
            f"tags={result['stochastic_AB_passing_tag_rates_h_inv'] or 'none'} | "
            f"pass={bool(result['stochastic_AB_interface_pass'])}"
        )
    if args.mode == "smoke":
        print("\nSmoke mode is a workflow check only; use audit/confirmation for inference.")


if __name__ == "__main__":
    main()
