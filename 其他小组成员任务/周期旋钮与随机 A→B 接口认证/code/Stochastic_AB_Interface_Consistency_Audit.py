r"""Audit the stochastic-A to deterministic-B flux interface for robust K4/K5.

This is a diagnostic, not a parameter fit and not a certification run.  Each
stochastic A trajectory is coupled to B several times:

* raw one-minute translation-flux samples;
* causal moving averages over 2, 5, 10 and 20 minutes.

The corresponding deterministic A waveform is included as a reference.  The
test asks whether the poor stochastic A-to-B fidelity is robust to the temporal
coarse-graining required when a discrete tau-leap event rate is supplied to a
continuous B-module ODE.

Run from the project root::

    python .\code\Stochastic_AB_Interface_Consistency_Audit.py --mode smoke --tau-step-min 0.05
    python .\code\Stochastic_AB_Interface_Consistency_Audit.py --mode audit --tau-step-min 0.05

No parameter is selected from the best-looking curve.  The output separates
time-resolution sensitivity from biological waveform sensitivity.
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
from scipy.integrate import trapezoid

import Robust_K4_K5_Stochastic_AB_Fidelity_Certification as ab
import Robust_Continuous_Period_Knob_AB_Certification as deterministic_a
import Stochastic_Period_Knob_Library_Certification as stochastic_a


VERSION = "Stochastic_AB_Interface_Consistency_Audit_v2_20260823"
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
OUTPUT_ROOT = PROJECT_ROOT / "outputs"
FIGURE_ROOT = PROJECT_ROOT / "figures"

TAG_RATE_H_INV = 12.0
TD_MIN = 50.0
SMOOTHING_WINDOWS_MIN = (0.0, 2.0, 5.0, 10.0, 20.0)


@dataclass(frozen=True)
class AuditSettings:
    trajectories_per_setting: int
    duration_generations: float
    transient_generations: float
    audited_cycles: int
    tau_step_min: float
    sample_step_min: float
    random_seed: int = 20260822


MODE_SETTINGS = {
    "smoke": AuditSettings(2, 135.0, 20.0, 7, 0.05, 1.0),
    "audit": AuditSettings(12, 150.0, 25.0, 7, 0.05, 1.0),
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
    parser.add_argument("--library-csv", type=Path, default=None)
    parser.add_argument("--b-model-dir", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument(
        "--smoothing-windows-min",
        default=",".join(f"{value:g}" for value in SMOOTHING_WINDOWS_MIN),
        help="Comma-separated causal averaging windows in minutes; include 0 for raw.",
    )
    parser.add_argument(
        "--fixed-plasmids",
        action="store_true",
        help="Diagnostic only: remove plasmid-copy fluctuations from stochastic A.",
    )
    return parser.parse_args()


def parse_windows(text: str) -> tuple[float, ...]:
    values = tuple(sorted({float(item.strip()) for item in text.split(",") if item.strip()}))
    if not values or any(value < 0.0 for value in values):
        raise ValueError("Smoothing windows must be non-negative numbers.")
    return values


def make_stochastic_settings(settings: AuditSettings) -> stochastic_a.RunSettings:
    return stochastic_a.RunSettings(
        trajectories_per_setting=settings.trajectories_per_setting,
        duration_generations=settings.duration_generations,
        transient_generations=settings.transient_generations,
        tau_step_min=settings.tau_step_min,
        sample_step_min=settings.sample_step_min,
        bootstrap_replicates=500,
        random_seed=settings.random_seed,
        minimum_cycles=max(4, settings.audited_cycles),
    )


def causal_moving_average(
    signal: np.ndarray, sample_step_min: float, window_min: float
) -> np.ndarray:
    """Dose-preserving causal rectangular average on a uniform time grid."""
    values = np.maximum(np.asarray(signal, dtype=float), 0.0)
    if window_min <= 0.0:
        return values.copy()
    width = max(1, int(round(window_min / sample_step_min)))
    cumulative = np.concatenate(([0.0], np.cumsum(values)))
    result = np.empty_like(values)
    for index in range(values.size):
        left = max(0, index + 1 - width)
        result[index] = (cumulative[index + 1] - cumulative[left]) / (index + 1 - left)
    return result


def cycle_flux_metrics(
    time_h: np.ndarray, flux_uM_h: np.ndarray, boundaries_min: np.ndarray
) -> list[dict]:
    boundaries_h = np.asarray(boundaries_min, dtype=float) / 60.0
    rows = []
    for cycle_index, (left, right) in enumerate(zip(boundaries_h[:-1], boundaries_h[1:])):
        use = (time_h >= left) & (time_h <= right)
        local_t = time_h[use]
        local_flux = flux_uM_h[use]
        if local_t.size < 2:
            continue
        peak_index = int(np.argmax(local_flux))
        peak = float(local_flux[peak_index])
        half = 0.5 * peak
        above = local_flux >= half
        fwhm = float(np.sum(above) * np.median(np.diff(local_t)))
        dose = float(trapezoid(local_flux, x=local_t))
        post_peak_dose = float(trapezoid(local_flux[peak_index:], x=local_t[peak_index:]))
        rows.append(
            {
                "cycle_index": cycle_index,
                "flux_peak_uM_h": peak,
                "flux_FWHM_h": fwhm,
                "flux_dose_uM": dose,
                "post_peak_dose_uM": post_peak_dose,
            }
        )
    return rows


def deterministic_reference(
    library_row: pd.Series,
    robust_row: pd.Series,
    settings: AuditSettings,
):
    duration_h = settings.duration_generations * TD_MIN / 60.0
    parameters, _ = deterministic_a.parameters_from_design(robust_row, 0.0)
    reference_peak = deterministic_a.reference_demand_peak(
        deterministic_a.core.make_reference_parameters(), duration_h
    )
    time_h, _, _, _, flux_uM_h = deterministic_a.simulate_a(
        parameters, reference_peak, duration_h
    )
    boundaries_h, realised_period_h = deterministic_a.cycle_boundaries(
        time_h,
        flux_uM_h,
        settings.transient_generations * TD_MIN / 60.0,
        settings.audited_cycles,
        float(library_row["realised_period_generations"]) * TD_MIN / 60.0,
    )
    return time_h, flux_uM_h, boundaries_h * 60.0, realised_period_h


def summarize_cycle_frame(frame: pd.DataFrame) -> pd.DataFrame:
    grouping = ["knob_setting_id", "interface_mode", "initial_DNA_state"]
    rows = []
    for keys, group in frame.groupby(grouping, sort=True):
        correct = group["cycle_correct"].astype(bool).to_numpy()
        crossings = pd.to_numeric(group["total_crossings"]).to_numpy(dtype=int)
        rows.append(
            {
                "knob_setting_id": keys[0],
                "interface_mode": keys[1],
                "initial_DNA_state": keys[2],
                "n_cycles": len(group),
                "cycle_success_fraction": float(np.mean(correct)),
                "zero_crossing_fraction": float(np.mean(crossings == 0)),
                "single_crossing_fraction": float(np.mean(crossings == 1)),
                "multiple_crossing_fraction": float(np.mean(crossings > 1)),
                "endpoint_complete_fraction": float(
                    np.mean(group["endpoint_complete"].astype(bool))
                ),
                "mean_input_peak_uM_h": float(group["input_peak_uM_h"].mean()),
                "mean_input_dose_uM": float(group["input_dose_uM"].mean()),
            }
        )
    return pd.DataFrame(rows)


def summarize_flux_frame(frame: pd.DataFrame) -> pd.DataFrame:
    grouping = ["knob_setting_id", "interface_mode"]
    aggregations = {
        "flux_peak_uM_h": ["mean", "median"],
        "flux_FWHM_h": ["mean", "median"],
        "flux_dose_uM": ["mean", "median"],
        "post_peak_dose_uM": ["mean", "median"],
    }
    summary = frame.groupby(grouping, sort=True).agg(aggregations)
    summary.columns = [f"{name}_{stat}" for name, stat in summary.columns]
    return summary.reset_index()


def make_figures(
    condition_summary: pd.DataFrame,
    flux_summary: pd.DataFrame,
    examples: dict,
    figure_dir: Path,
):
    settings = tuple(sorted(condition_summary["knob_setting_id"].unique()))
    modes = list(dict.fromkeys(condition_summary["interface_mode"]))
    fig, axes = plt.subplots(2, len(settings), figsize=(6 * len(settings), 8), squeeze=False)
    for column, setting_id in enumerate(settings):
        local_b = condition_summary[condition_summary["knob_setting_id"].eq(setting_id)]
        worst = local_b.groupby("interface_mode")["cycle_success_fraction"].min()
        ordered = [mode for mode in modes if mode in worst.index]
        axes[0, column].bar(np.arange(len(ordered)), worst.loc[ordered])
        axes[0, column].axhline(0.95, color="black", linestyle="--", linewidth=1)
        axes[0, column].set_xticks(np.arange(len(ordered)), ordered, rotation=35, ha="right")
        axes[0, column].set_ylim(0, 1.05)
        axes[0, column].set_title(f"{setting_id}: worst PB/LR cycle success")
        axes[0, column].set_ylabel("fraction correct")

        local_f = flux_summary[flux_summary["knob_setting_id"].eq(setting_id)].set_index(
            "interface_mode"
        )
        ordered_f = [mode for mode in modes if mode in local_f.index]
        axes[1, column].plot(
            np.arange(len(ordered_f)),
            local_f.loc[ordered_f, "flux_peak_uM_h_median"],
            marker="o",
            label="median peak (uM/h)",
        )
        axes[1, column].plot(
            np.arange(len(ordered_f)),
            local_f.loc[ordered_f, "flux_dose_uM_median"],
            marker="s",
            label="median dose (uM)",
        )
        axes[1, column].set_xticks(
            np.arange(len(ordered_f)), ordered_f, rotation=35, ha="right"
        )
        axes[1, column].set_title(f"{setting_id}: peak versus dose")
        axes[1, column].legend()
    fig.suptitle("Stochastic A-to-B interface time-resolution audit")
    fig.tight_layout()
    fig.savefig(figure_dir / "Stochastic_AB_interface_consistency_summary.png", dpi=200)
    plt.close(fig)

    if examples:
        fig, axes = plt.subplots(len(settings), 1, figsize=(12, 4 * len(settings)), squeeze=False)
        for row_index, setting_id in enumerate(settings):
            payload = examples.get(setting_id)
            if payload is None:
                continue
            time_h, signals = payload
            show = time_h >= max(time_h[-1] - 35.0, 0.0)
            for label, signal in signals.items():
                axes[row_index, 0].plot(time_h[show], signal[show], label=label, linewidth=1.2)
            axes[row_index, 0].set_title(setting_id)
            axes[row_index, 0].set_ylabel("C31 translation flux (uM/h)")
            axes[row_index, 0].legend(ncol=3)
        axes[-1, 0].set_xlabel("time (h)")
        fig.suptitle("Representative raw and causally averaged stochastic flux")
        fig.tight_layout()
        fig.savefig(figure_dir / "Stochastic_AB_interface_flux_examples.png", dpi=200)
        plt.close(fig)


def main():
    args = parse_args()
    settings = MODE_SETTINGS[args.mode]
    if args.tau_step_min is not None:
        if not math.isfinite(args.tau_step_min) or args.tau_step_min <= 0.0:
            raise ValueError("--tau-step-min must be a finite positive number")
        settings = replace(settings, tau_step_min=float(args.tau_step_min))
    stochastic_settings = make_stochastic_settings(settings)
    windows = parse_windows(args.smoothing_windows_min)
    library_path, library = ab.load_k4_k5_library(args.library_csv)
    library = ab.validate_k4_k5_library(library)
    robust_source = deterministic_a.load_representatives(library_path)
    robust_source = robust_source.sort_values("target_period_generations").reset_index(drop=True)
    robust_lookup = robust_source.set_index("setting_id")
    zc, zhao_path = ab.import_zhao_core(args.b_model_dir)

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = (
        args.output_dir.expanduser().resolve()
        if args.output_dir is not None
        else OUTPUT_ROOT
        / (
            f"Stochastic_AB_Interface_Audit_{args.mode}_"
            f"tau{settings.tau_step_min:g}min_{stamp}"
        )
    )
    figure_dir = FIGURE_ROOT / output_dir.name
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    print(f"{VERSION}: mode={args.mode}")
    print(f"Robust library: {library_path}")
    print(f"B model: {zhao_path}")
    print(f"Tau-leap step (min): {settings.tau_step_min:g}")
    print(f"Causal averaging windows (min): {windows}")

    flux_rows = []
    b_cycle_rows = []
    trajectory_rows = []
    examples = {}

    for setting_index, (_, library_row) in enumerate(library.iterrows()):
        setting_id = str(library_row["knob_setting_id"])
        source_id = str(library_row["source_setting_id"])
        robust_row = robust_lookup.loc[source_id]

        # Deterministic reference uses the same physical K4/K5 design.
        det_time_h, det_flux, det_boundaries_min, det_period_h = deterministic_reference(
            library_row, robust_row, settings
        )
        for metric in cycle_flux_metrics(det_time_h, det_flux, det_boundaries_min):
            flux_rows.append(
                {
                    "knob_setting_id": setting_id,
                    "trajectory_index": -1,
                    "interface_mode": "deterministic_reference",
                    **metric,
                }
            )
        for initial_name in ("PB", "LR"):
            b_time, lr, _, _, _ = ab.solve_b(
                zc, det_time_h, det_flux, TAG_RATE_H_INV, TD_MIN, initial_name
            )
            cycles, _ = ab.audit_cycles(
                b_time, lr, det_time_h, det_flux, det_boundaries_min
            )
            for cycle in cycles:
                b_cycle_rows.append(
                    {
                        "knob_setting_id": setting_id,
                        "trajectory_index": -1,
                        "interface_mode": "deterministic_reference",
                        "initial_DNA_state": initial_name,
                        **cycle,
                    }
                )
        trajectory_rows.append(
            {
                "knob_setting_id": setting_id,
                "trajectory_index": -1,
                "interface_mode": "deterministic_reference",
                "A_period_generations": det_period_h * 60.0 / TD_MIN,
                "raw_flux_peak_uM_h": float(np.max(det_flux)),
            }
        )

        for trajectory_index in range(settings.trajectories_per_setting):
            seed = settings.random_seed + setting_index * 100_000 + trajectory_index
            print(
                f"[{setting_id} {trajectory_index + 1:02d}/"
                f"{settings.trajectories_per_setting:02d}] stochastic A"
            )
            result = stochastic_a.simulate_trajectory(
                library_row, stochastic_settings, seed, args.fixed_plasmids
            )
            analysis = stochastic_a.analyse_trajectory(
                library_row, result, stochastic_settings
            )
            deterministic_period_min = (
                float(library_row["realised_period_generations"]) * TD_MIN
            )
            boundaries_min = ab.choose_cycle_boundaries(
                analysis["phase_crossings_min"],
                deterministic_period_min,
                settings.audited_cycles,
            )
            if boundaries_min.size != settings.audited_cycles + 1:
                print("    skipped: insufficient valid stochastic phase boundaries")
                continue
            time_h = result["time_min"] / 60.0
            raw_flux = ab.flux_to_uM_h(result["c31_flux_molecules_per_min"])
            signals = {}
            for window in windows:
                mode = "stochastic_raw" if window == 0.0 else f"stochastic_causal_{window:g}min"
                flux = causal_moving_average(raw_flux, settings.sample_step_min, window)
                signals[mode] = flux
                metrics = cycle_flux_metrics(time_h, flux, boundaries_min)
                for metric in metrics:
                    flux_rows.append(
                        {
                            "knob_setting_id": setting_id,
                            "trajectory_index": trajectory_index,
                            "interface_mode": mode,
                            **metric,
                        }
                    )
                for initial_name in ("PB", "LR"):
                    b_time, lr, _, _, _ = ab.solve_b(
                        zc, time_h, flux, TAG_RATE_H_INV, TD_MIN, initial_name
                    )
                    cycles, _ = ab.audit_cycles(
                        b_time, lr, time_h, flux, boundaries_min
                    )
                    for cycle in cycles:
                        b_cycle_rows.append(
                            {
                                "knob_setting_id": setting_id,
                                "trajectory_index": trajectory_index,
                                "interface_mode": mode,
                                "initial_DNA_state": initial_name,
                                **cycle,
                            }
                        )
                trajectory_rows.append(
                    {
                        "knob_setting_id": setting_id,
                        "trajectory_index": trajectory_index,
                        "seed": seed,
                        "interface_mode": mode,
                        "A_resolved": bool(analysis["resolved"]),
                        "A_period_generations": analysis["phase_period_mean_generations"],
                        "A_within_trajectory_period_cv": analysis[
                            "phase_period_within_trajectory_cv"
                        ],
                        "raw_flux_peak_uM_h": float(np.max(raw_flux)),
                        "processed_flux_peak_uM_h": float(np.max(flux)),
                        "whole_trace_dose_ratio_vs_raw": float(
                            trapezoid(flux, x=time_h) / max(trapezoid(raw_flux, x=time_h), 1e-12)
                        ),
                    }
                )
            if trajectory_index == 0:
                examples[setting_id] = (time_h, signals)

    flux_frame = pd.DataFrame(flux_rows)
    cycle_frame = pd.DataFrame(b_cycle_rows)
    trajectory_frame = pd.DataFrame(trajectory_rows)
    flux_summary = summarize_flux_frame(flux_frame)
    condition_summary = summarize_cycle_frame(cycle_frame)

    flux_frame.to_csv(output_dir / "cycle_flux_metrics.csv", index=False)
    flux_summary.to_csv(output_dir / "interface_flux_summary.csv", index=False)
    cycle_frame.to_csv(output_dir / "B_cycle_metrics_by_interface.csv", index=False)
    condition_summary.to_csv(output_dir / "B_interface_condition_summary.csv", index=False)
    trajectory_frame.to_csv(output_dir / "trajectory_interface_audit.csv", index=False)

    metadata = {
        "version": VERSION,
        "mode": args.mode,
        "settings": asdict(settings),
        "library_csv": str(library_path),
        "zhao_core": str(zhao_path),
        "tag_rate_h_inv": TAG_RATE_H_INV,
        "Td_min": TD_MIN,
        "smoothing_windows_min": windows,
        "interpretation": (
            "Sensitivity to short causal averaging supports a hybrid-interface "
            "time-resolution artifact; persistent failure after averaging supports "
            "a biological waveform/noise limitation."
        ),
    }
    (output_dir / "scope_and_settings.json").write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    make_figures(condition_summary, flux_summary, examples, figure_dir)

    print(f"\nCompleted: {VERSION}")
    print(f"Raw tables: {output_dir}")
    print(f"Figures: {figure_dir}")
    print("\nWorst PB/LR cycle-success fraction by interface:")
    worst = (
        condition_summary.groupby(["knob_setting_id", "interface_mode"])[
            "cycle_success_fraction"
        ]
        .min()
        .reset_index()
    )
    for _, row in worst.iterrows():
        print(
            f"{row['knob_setting_id']:>9s} | {row['interface_mode']:<28s} | "
            f"worst cycle={row['cycle_success_fraction']:.3f}"
        )
    print(
        "\nThis diagnostic does not certify a construct. Compare raw, causally "
        "averaged and deterministic-reference rows before changing parameters."
    )


if __name__ == "__main__":
    main()
