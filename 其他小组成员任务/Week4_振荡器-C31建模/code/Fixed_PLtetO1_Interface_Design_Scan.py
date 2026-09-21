# Standalone program: fixed-PLtetO1 A/B interface engineering.
# Source version v53g. The packaged zhao_core.py is loaded from this folder.
#
# Biological constraint enforced here:
# - C31 and cI use the SAME PLtetO1 transfer function.
# - K, n, leak, period and pulse shape are NOT changed.
# - Only C31 RBS/translation gain and Int degradation-tag strength are scanned.
#
# A input to B:
#   int_production(t) = RBS_scale * v53d_no_load_C31_flux(t)
# B owns Int binding, dilution and tag-mediated loss.

from datetime import datetime
from pathlib import Path
import csv, json, math, sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp, trapezoid
from scipy.interpolate import CubicSpline
from scipy.signal import find_peaks

# Portable project paths: independent of the terminal's current directory.
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
DATA_DIR = PROJECT_ROOT / "data"
OUTPUT_ROOT = PROJECT_ROOT / "outputs"

try:
    from IPython.display import display
except ImportError:
    display = None

VERSION = "v53g_fixed_PLtetO1_RBS_tag_window_20260814"

# ============================== USER SETTINGS ==============================
SETTINGS = {
    "v53d_output_directory": None,  # None -> latest matching output folder
    # RBS scale multiplies beta31 and therefore peak/dose, not normalized width.
    "C31_RBS_scales": (0.45, 0.60, 0.70, 0.80, 1.00),
    # h^-1, representing candidate Int degradation tags in the Zhao B model.
    "Int_tag_rates_h_inv": (4.0, 8.0, 12.0, 16.0, 20.0),
    "A_transient_h": 25.0,
    "minimum_peak_distance_h": 5.0,
    "B_duration_h": 90.0,
    "B_sample_step_h": 0.02,
    "B_max_step_h": 0.05,
    "DNA_midpoint": 0.5,
    "complete_state_threshold": 0.95,
    "success_target": 0.95,
    "alternation_target": 0.95,
    "DNA_conservation_tolerance": 1e-8,
    "expected_Td_min": 50.0,
    "make_plots": True,
}
# ============================================================================


def locate_v53d_root():
    root = DATA_DIR / "reference_inputs"
    if not root.is_dir():
        raise FileNotFoundError(root)
    return root


def load_fixed_PLtetO1_flux(root):
    path = root / "unloaded_C31_translation_trajectory.csv"
    if not path.is_file():
        raise FileNotFoundError(path)
    frame = pd.read_csv(path)
    required = {"time_min", "actual_C31_translation_flux_uM_h"}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"Trajectory is missing columns: {sorted(missing)}")
    return (
        path,
        frame["time_min"].to_numpy(float)/60.0,
        np.maximum(
            frame["actual_C31_translation_flux_uM_h"].to_numpy(float), 0.0
        ),
    )


def peak_times_and_boundaries(time_h, flux):
    keep = time_h >= SETTINGS["A_transient_h"]
    indices = np.flatnonzero(keep)
    y = flux[keep]
    dt = float(np.median(np.diff(time_h)))
    span = float(np.ptp(y))
    peaks, _ = find_peaks(
        y,
        prominence=max(0.05*span, 1e-9),
        distance=max(1, int(round(SETTINGS["minimum_peak_distance_h"]/dt))),
    )
    times = time_h[indices[peaks]]
    times = times[times <= SETTINGS["B_duration_h"]]
    if len(times) < 5:
        raise RuntimeError("Fewer than five post-transient C31 pulses")
    period = float(np.mean(np.diff(times[-6:])))
    boundaries = 0.5*(times[:-1]+times[1:])
    return times, boundaries, period


def baseline_corrected_FWHM(time_h, flux, boundaries):
    widths = []
    for left, right in zip(boundaries[:-1], boundaries[1:]):
        choose = (time_h >= left) & (time_h < right)
        t = time_h[choose]; y = flux[choose]
        if len(t) < 3:
            continue
        low = float(np.min(y)); high = float(np.max(y))
        level = low+0.5*(high-low)
        widths.append(float(np.sum(y >= level)*np.median(np.diff(t))))
    return float(np.median(widths)) if widths else float("nan")


def crossing_events(t, lr, left, right):
    local = np.flatnonzero((t >= left) & (t < right))
    events = []
    for i0, i1 in zip(local[:-1], local[1:]):
        a = lr[i0]-SETTINGS["DNA_midpoint"]
        b = lr[i1]-SETTINGS["DNA_midpoint"]
        if a*b < 0.0:
            fraction = abs(a)/(abs(a)+abs(b))
            event_time = t[i0]+fraction*(t[i1]-t[i0])
            direction = "PB_to_LR" if b > a else "LR_to_PB"
            events.append((event_time, direction))
    return events


def integrate_spline(spline, left, right, step=0.002):
    grid = np.arange(left, right+0.5*step, step)
    values = np.maximum(0.0, spline(grid))
    return float(trapezoid(values, x=grid))


def alternation_fidelity(cycle_rows):
    eligible = []
    for previous, current in zip(cycle_rows[:-1], cycle_rows[1:]):
        if previous["total_crossings"] == 1 and current["total_crossings"] == 1:
            eligible.append(
                previous["single_crossing_direction"]
                != current["single_crossing_direction"]
            )
    return (
        float(np.mean(eligible)) if eligible else float("nan"), len(eligible)
    )


def write_rows(path, rows):
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def metric_matrix(rows, metric):
    tags = list(SETTINGS["Int_tag_rates_h_inv"])
    scales = list(SETTINGS["C31_RBS_scales"])
    matrix = np.full((len(tags), len(scales)), np.nan)
    for row in rows:
        i = tags.index(row["Int_tag_rate_h_inv"])
        j = scales.index(row["C31_RBS_scale"])
        matrix[i, j] = row[metric]
    return matrix


def annotate_heatmap(axis, matrix):
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = matrix[i, j]
            if np.isfinite(value):
                color = "white" if value < 0.35 or value > 0.75 else "black"
                axis.text(j, i, f"{value:.2f}", ha="center", va="center",
                          color=color, fontsize=8)


def main():
    # ---------------------------- Inputs and model ------------------------------
    MODEL_DIR = SCRIPT_DIR
    if not (MODEL_DIR/"zhao_core.py").is_file():
        raise FileNotFoundError(
            "Run this notebook from the B model folder containing zhao_core.py"
        )
    if str(MODEL_DIR) not in sys.path:
        sys.path.insert(0, str(MODEL_DIR))
    global zc
    import zhao_core as zc
    
    v53d_root = locate_v53d_root()
    source_path, A_time_h, source_flux = load_fixed_PLtetO1_flux(v53d_root)
    peak_times, cycle_boundaries, PERIOD_H = peak_times_and_boundaries(
        A_time_h, source_flux
    )
    SOURCE_PEAK = float(np.max(source_flux[A_time_h >= SETTINGS["A_transient_h"]]))
    SOURCE_FWHM_H = baseline_corrected_FWHM(
        A_time_h, source_flux, cycle_boundaries
    )
    
    probe = zc.default_params()
    expected_k_dil = math.log(2.0)/(SETTINGS["expected_Td_min"]/60.0)
    if not math.isclose(
        probe["k_dil"], expected_k_dil, rel_tol=2e-4, abs_tol=2e-4
    ):
        raise RuntimeError(
            f"A/B growth mismatch: B k_dil={probe['k_dil']:.6g} h^-1; "
            f"Td={SETTINGS['expected_Td_min']:g} min requires "
            f"{expected_k_dil:.6g} h^-1"
        )
    
    output = OUTPUT_ROOT / (
        "Fixed_PLtetO1_Interface_Design_Scan_"+datetime.now().strftime("%Y%m%d_%H%M%S")
    )
    output.mkdir(parents=True, exist_ok=False)
    
    # ------------------------------ 2D scan -------------------------------------
    t_eval = np.arange(0.0, SETTINGS["B_duration_h"], SETTINGS["B_sample_step_h"])
    conditions = [
        (float(tag), float(scale))
        for tag in SETTINGS["Int_tag_rates_h_inv"]
        for scale in SETTINGS["C31_RBS_scales"]
    ]
    summary_rows = []
    cycle_rows = []
    saved = {}
    
    for run_number, (tag, scale) in enumerate(conditions, start=1):
        label = f"tag_{tag:g}__RBS_scale_{scale:g}"
        print(f"[{run_number:02d}/{len(conditions):02d}] Running {label}")
        scaled_flux = scale*source_flux
        spline = CubicSpline(A_time_h, scaled_flux, bc_type="natural")
    
        def int_production(t):
            return max(0.0, float(spline(float(t))))
    
        P = zc.default_params()
        P.update(
            K_rep=0.0186, n_rep=3.4, krep_tsl=15.0,
            krdf_tsl=200.0, k_tag_int=tag,
        )
        C = zc._rate_constants(P)
        y0 = zc.y0_PB(
            P,
            rep_mrna=P["k_tscr"]*P["Dtot"]/P["k_rna"],
            rep=(15.0*P["k_tscr"]*P["Dtot"]
                 /P["k_rna"]/P["k_dil"]),
        )
        sol = solve_ivp(
            lambda t, y: zc.rhs(t, y, P, C, int_production),
            (0.0, SETTINGS["B_duration_h"]), y0, t_eval=t_eval,
            method="LSODA", rtol=1e-7, atol=1e-12,
            max_step=SETTINGS["B_max_step_h"],
        )
        if not sol.success:
            raise RuntimeError(f"{label}: {sol.message}")
    
        t = sol.t
        states = sol.y.T
        lr = zc.LR_total(states)/P["Dtot"]
        pb = zc.PB_total(states)/P["Dtot"]
        int_total = zc.int_total(states)
        rdf_total = zc.rdf_total(states)
        conservation = np.abs(lr+pb-1.0)
        local_rows = []
    
        for cycle_index, (left, right) in enumerate(
            zip(cycle_boundaries[:-1], cycle_boundaries[1:]), start=1
        ):
            if right > SETTINGS["B_duration_h"]:
                continue
            choose = (t >= left) & (t < right)
            indices = np.flatnonzero(choose)
            if len(indices) < 3:
                continue
            events = crossing_events(t, lr, left, right)
            end = indices[-1]
            complete_fraction = max(float(lr[end]), float(pb[end]))
            cycle_dose = integrate_spline(spline, left, right)
            residual_dose = float("nan")
            if events:
                residual_dose = integrate_spline(spline, events[0][0], right)
            row = {
                "C31_RBS_scale": scale,
                "Int_tag_rate_h_inv": tag,
                "effective_Int_total_half_life_min": (
                    60.0*math.log(2.0)/(tag+P["k_dil"])
                ),
                "cycle_index": cycle_index,
                "cycle_start_h": left,
                "cycle_end_h": right,
                "period_h": right-left,
                "total_crossings": len(events),
                "PB_to_LR_crossings": sum(
                    direction == "PB_to_LR" for _, direction in events
                ),
                "LR_to_PB_crossings": sum(
                    direction == "LR_to_PB" for _, direction in events
                ),
                "single_crossing_direction": (
                    events[0][1] if len(events) == 1 else ""
                ),
                "LR_start": float(lr[indices[0]]),
                "LR_end": float(lr[end]),
                "PB_end": float(pb[end]),
                "endpoint_complete_fraction": complete_fraction,
                "endpoint_complete": bool(
                    complete_fraction >= SETTINGS["complete_state_threshold"]
                ),
                "input_peak_uM_h": float(np.max(scaled_flux[
                    (A_time_h >= left) & (A_time_h < right)
                ])),
                "input_FWHM_h": SOURCE_FWHM_H,
                "input_dose_uM_cycle": cycle_dose,
                "residual_dose_after_first_crossing_uM": residual_dose,
                "Int_peak_uM": float(np.max(int_total[choose])),
                "RDF_peak_uM": float(np.max(rdf_total[choose])),
                "max_DNA_conservation_residual": float(
                    np.max(conservation[choose])
                ),
            }
            row["one_and_complete"] = bool(
                row["total_crossings"] == 1 and row["endpoint_complete"]
            )
            local_rows.append(row)
            cycle_rows.append(row)
    
        counts = np.array([row["total_crossings"] for row in local_rows])
        endpoints = np.array([row["endpoint_complete"] for row in local_rows])
        combined = np.array([row["one_and_complete"] for row in local_rows])
        alternation, eligible_pairs = alternation_fidelity(local_rows)
        residual = np.array([
            row["residual_dose_after_first_crossing_uM"] for row in local_rows
        ], dtype=float)
        residual = residual[np.isfinite(residual)]
        summary = {
            "C31_RBS_scale": scale,
            "Int_tag_rate_h_inv": tag,
            "effective_Int_total_half_life_min": (
                60.0*math.log(2.0)/(tag+P["k_dil"])
            ),
            "period_h": PERIOD_H,
            "fixed_PLtetO1_FWHM_h": SOURCE_FWHM_H,
            "input_peak_uM_h": scale*SOURCE_PEAK,
            "mean_input_dose_uM_cycle": float(np.mean([
                row["input_dose_uM_cycle"] for row in local_rows
            ])),
            "audited_cycles": len(local_rows),
            "fraction_exactly_one_crossing": float(np.mean(counts == 1)),
            "fraction_zero_crossings": float(np.mean(counts == 0)),
            "fraction_multiple_crossings": float(np.mean(counts > 1)),
            "fraction_complete_endpoint": float(np.mean(endpoints)),
            "fraction_one_and_complete": float(np.mean(combined)),
            "alternation_fidelity": alternation,
            "alternation_eligible_pairs": eligible_pairs,
            "mean_residual_dose_after_first_crossing_uM": (
                float(np.mean(residual)) if len(residual) else float("nan")
            ),
            "max_DNA_conservation_residual": float(np.max(conservation)),
        }
        summary["system_pass"] = bool(
            summary["fraction_one_and_complete"] >= SETTINGS["success_target"]
            and summary["fraction_complete_endpoint"] >= SETTINGS["success_target"]
            and np.isfinite(summary["alternation_fidelity"])
            and summary["alternation_fidelity"] >= SETTINGS["alternation_target"]
            and summary["max_DNA_conservation_residual"]
            <= SETTINGS["DNA_conservation_tolerance"]
        )
        summary_rows.append(summary)
        saved[(tag, scale)] = {
            "time": t, "input": np.maximum(0.0, spline(t)),
            "lr": lr, "pb": pb,
        }
    
    write_rows(output/"v53g_fixed_PLtetO1_RBS_tag_summary.csv", summary_rows)
    write_rows(output/"v53g_fixed_PLtetO1_RBS_tag_cycle_audit.csv", cycle_rows)
    
    with (output/"v53g_scope_and_parameters.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump({
            "version": VERSION,
            "settings": SETTINGS,
            "source_v53d_directory": str(v53d_root.resolve()),
            "source_trajectory": str(source_path.resolve()),
            "fixed_period_h": PERIOD_H,
            "fixed_PLtetO1_FWHM_h": SOURCE_FWHM_H,
            "fixed_source_peak_uM_h": SOURCE_PEAK,
            "enforced_construct_constraint": (
                "C31 and cI share the same PLtetO1; K, n, leak and normalized "
                "pulse shape are not scanned"
            ),
            "RBS_model": "int_production(t)=RBS_scale*source_C31_flux(t)",
            "tag_model": "k_tag_int is varied inside the Zhao B ODE",
            "acceptance": (
                "one-and-complete >=0.95, endpoint complete >=0.95, "
                "alternation >=0.95, DNA conservation <=1e-8"
            ),
            "interpretation_limits": [
                "RBS scales and tag rates are candidate designs, not measurements",
                "RBS scaling is first screened open-loop; any candidate must later be rechecked with resource feedback",
                "if no point passes, the fixed-promoter topology needs an added post-transcriptional pulse-shaping element or B redesign",
            ],
        }, handle, ensure_ascii=False, indent=2)
    
    # --------------------------------- Figures ---------------------------------
    if SETTINGS["make_plots"]:
        metrics = (
            ("fraction_one_and_complete", "one crossing + complete"),
            ("fraction_zero_crossings", "zero crossings"),
            ("fraction_multiple_crossings", "multiple crossings"),
        )
        tags = list(SETTINGS["Int_tag_rates_h_inv"])
        scales = list(SETTINGS["C31_RBS_scales"])
        fig, axes = plt.subplots(1, len(metrics), figsize=(16, 5),
                                 constrained_layout=True)
        for axis, (metric, title) in zip(axes, metrics):
            matrix = metric_matrix(summary_rows, metric)
            image = axis.imshow(
                matrix, origin="lower", vmin=0.0, vmax=1.0,
                cmap="viridis", aspect="auto",
            )
            annotate_heatmap(axis, matrix)
            axis.set_xticks(range(len(scales)), labels=[f"{x:g}" for x in scales])
            axis.set_yticks(range(len(tags)), labels=[f"{x:g}" for x in tags])
            axis.set_xlabel("C31 RBS / translation scale")
            axis.set_ylabel("k_tag_int (h^-1)")
            axis.set_title(title)
            fig.colorbar(image, ax=axis, shrink=0.8)
        fig.suptitle(
            f"v53g fixed PLtetO1: FWHM={SOURCE_FWHM_H:.2f} h; "
            "RBS x Int-tag window"
        )
        fig.savefig(output/"v53g_RBS_tag_heatmaps.png", dpi=180)
        if display is not None:
            display(fig)
        plt.close(fig)
    
        # Plot up to three passing candidates; otherwise show the three highest
        # success scores.  This is visualization, not automatic parameter choice.
        ordered = sorted(
            summary_rows,
            key=lambda row: (
                row["system_pass"], row["fraction_one_and_complete"],
                row["fraction_complete_endpoint"],
                -row["fraction_multiple_crossings"],
            ),
            reverse=True,
        )[:3]
        fig, axes = plt.subplots(len(ordered), 2, figsize=(13, 9), sharex=True,
                                 constrained_layout=True)
        if len(ordered) == 1:
            axes = np.array([axes])
        for row_index, row in enumerate(ordered):
            data = saved[(row["Int_tag_rate_h_inv"], row["C31_RBS_scale"])]
            axes[row_index, 0].plot(data["time"], data["input"], color="m")
            axes[row_index, 0].set_ylabel(
                f"RBS={row['C31_RBS_scale']:.2f}\ntag={row['Int_tag_rate_h_inv']:.0f}/h\ninput"
            )
            axes[row_index, 1].plot(data["time"], data["lr"], "r", label="LR")
            axes[row_index, 1].plot(data["time"], data["pb"], "b", label="PB")
            axes[row_index, 1].axhline(
                SETTINGS["DNA_midpoint"], color="0.4", linestyle="--"
            )
            axes[row_index, 1].set_ylim(-0.05, 1.05)
            axes[row_index, 1].set_ylabel("DNA fraction")
            axes[row_index, 1].legend()
            for axis in axes[row_index]:
                axis.set_xlim(25, SETTINGS["B_duration_h"])
                axis.grid(alpha=0.2)
        for axis in axes[-1]:
            axis.set_xlabel("time (h)")
        fig.savefig(output/"v53g_candidate_trajectories.png", dpi=180)
        if display is not None:
            display(fig)
        plt.close(fig)
    
    print("Loaded diagnostic code:", VERSION)
    print("Source fixed-PLtetO1 trajectory:", source_path.resolve())
    print(
        f"Frozen period={PERIOD_H:.3f} h; FWHM={SOURCE_FWHM_H:.3f} h; "
        f"unscaled peak={SOURCE_PEAK:.3f} uM/h"
    )
    print("Saved v53g outputs to:", output.resolve())
    for row in summary_rows:
        print(
            f"RBS={row['C31_RBS_scale']:.2f} | "
            f"tag={row['Int_tag_rate_h_inv']:>4.0f}/h "
            f"(t1/2={row['effective_Int_total_half_life_min']:.2f} min) | "
            f"peak={row['input_peak_uM_h']:.2f} | "
            f"dose={row['mean_input_dose_uM_cycle']:.2f} | "
            f"one+complete={row['fraction_one_and_complete']:.2f} | "
            f"zero={row['fraction_zero_crossings']:.2f} | "
            f"multiple={row['fraction_multiple_crossings']:.2f} | "
            f"alternation={row['alternation_fidelity']:.2f} | "
            f"pass={row['system_pass']}"
        )
    
    passing = [row for row in summary_rows if row["system_pass"]]
    print("\nPassing fixed-PLtetO1 designs:")
    if passing:
        for row in passing:
            print(
                f"  RBS scale={row['C31_RBS_scale']:.2f}, "
                f"k_tag_int={row['Int_tag_rate_h_inv']:.0f}/h, "
                f"effective total Int t1/2={row['effective_Int_total_half_life_min']:.2f} min, "
                f"peak={row['input_peak_uM_h']:.2f} uM/h, "
                f"dose={row['mean_input_dose_uM_cycle']:.2f} uM/cycle"
            )
    else:
        print("  None in the scanned RBS-tag grid.")
        print("  The unchanged 3.62-h PLtetO1 pulse cannot be rescued by these two levers alone.")
    
    print("\nInterpretation boundaries:")
    print("- PLtetO1 K, n, leak and normalized pulse shape are identical to cI and remain frozen.")
    print("- RBS changes peak/dose; the Int tag changes active-protein clearance in B.")
    print("- Any passing candidate must next be revalidated under the v53d resource-load model.")


if __name__ == "__main__":
    main()
