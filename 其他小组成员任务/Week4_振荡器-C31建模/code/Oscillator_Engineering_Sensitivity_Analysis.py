# Standalone program: one-at-a-time engineering sensitivity analysis.
# Source version v53a, built around the frozen v52a1 ODE.
# This script does NOT couple module B and does NOT run a stochastic model.

from dataclasses import dataclass, asdict, replace
from datetime import datetime
from pathlib import Path
import csv, json, math

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp, trapezoid
from scipy.optimize import brentq
from scipy.signal import find_peaks

# Portable project paths: independent of the terminal's current directory.
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
DATA_DIR = PROJECT_ROOT / "data"
OUTPUT_ROOT = PROJECT_ROOT / "outputs"

try:
    from IPython.display import display as notebook_display
except ImportError:
    notebook_display = None

VERSION = "v53a_engineering_oat_scan_20260814"

# ============================== USER SETTINGS ==============================
SETTINGS = {
    "duration_min": 6000.0,
    "sample_step_min": 1.0,
    "transient_min": 1500.0,
    "make_plots": True,
    "save_trajectories": True,
    "minimum_relative_amplitude": 0.01,
    "minimum_peaks": 4,
    # These conversion assumptions are only for displaying the C31 production
    # flux in uM/h. The ODE itself remains in copies/cell and minutes.
    "cell_volume_fL": 1.0,
    "active_C31_fraction": 1.0,
    "molecules_per_uM_per_fL": 602.214076,
}

# Reference design. Every scan changes only ONE entry from this design.
REFERENCE = {
    "Td_min": 50.0,
    "sponge_plasmid_copies": 40.0,
    "sponge_accessible_fraction": 1.0,
    "oscillator_promoter_scale": 1.0,
    "oscillator_RBS_scale": 1.0,
}

SCAN_AXES = {
    "oscillator_promoter_scale": (0.50, 0.75, 1.00, 1.25, 1.50),
    "oscillator_RBS_scale": (0.50, 0.75, 1.00, 1.25, 1.50),
    "sponge_plasmid_copies": (0.0, 20.0, 40.0, 60.0),
    "Td_min": (40.0, 50.0, 60.0),
}
# ============================================================================


@dataclass(frozen=True)
class Parameters:
    # Frozen v52a1 chemistry and plasmid accounting.
    Td_min: float = 50.0
    lambda_protein_scale: float = 1000.0
    K_copies: float = 13.0
    hill_n: float = 3.0
    oscillator_plasmid_copies: float = 10.0
    functional_tetO_per_plasmid: float = 4.0
    functional_accessible_fraction: float = 1.0
    sponge_plasmid_copies: float = 40.0
    sponge_tetO_per_plasmid: float = 2.0
    sponge_accessible_fraction: float = 1.0
    oscillator_promoter_leak: float = 0.0

    # IMPORTANT: alpha0, beta0 and intrinsic losses are calibrated ONCE at the
    # reference design. They are not remapped during the engineering scans.
    oscillator_tx_per_plasmid_per_min: float = 0.0
    oscillator_translation_per_mrna_per_min: float = 0.5
    oscillator_mrna_intrinsic_loss_per_min: float = 0.0
    oscillator_protein_intrinsic_loss_per_min: float = 0.0

    # C31 replaces the reporter on the same plasmid, so its tetO sites are
    # already included in functional_tetO_per_plasmid=4.
    c31_promoter_leak: float = 0.005
    c31_tx_per_plasmid_per_min: float = 5.0
    c31_translation_per_mrna_per_min: float = 0.5
    c31_mrna_intrinsic_loss_per_min: float = 0.0
    c31_protein_intrinsic_loss_per_min: float = 0.0

    @property
    def mu(self):
        return math.log(2.0) / self.Td_min

    @property
    def effective_functional_sites(self):
        return (self.functional_accessible_fraction
                * self.functional_tetO_per_plasmid
                * self.oscillator_plasmid_copies)

    @property
    def effective_sponge_sites(self):
        return (self.sponge_accessible_fraction
                * self.sponge_tetO_per_plasmid
                * self.sponge_plasmid_copies)


def make_reference_parameters():
    """Freeze the physical v52a1 rates at Td=50 min and mRNA t1/2=2 min."""
    Td0 = REFERENCE["Td_min"]
    mu0 = math.log(2.0) / Td0
    total_mrna_loss0 = math.log(2.0) / 2.0
    mrna_intrinsic_loss0 = total_mrna_loss0 - mu0
    beta0 = 0.5
    No = 10.0
    lambda0 = 1000.0
    protein_intrinsic_loss0 = 0.0
    # v52a1 QSS mapping, used only once to establish the reference alpha0:
    # beta*No*alpha/(delta_m+mu) = (delta_p+mu)*lambda.
    alpha0 = ((protein_intrinsic_loss0 + mu0) * lambda0 * total_mrna_loss0
              / (beta0 * No))
    return Parameters(
        Td_min=Td0,
        sponge_plasmid_copies=REFERENCE["sponge_plasmid_copies"],
        sponge_accessible_fraction=REFERENCE["sponge_accessible_fraction"],
        oscillator_tx_per_plasmid_per_min=alpha0,
        oscillator_translation_per_mrna_per_min=beta0,
        oscillator_mrna_intrinsic_loss_per_min=mrna_intrinsic_loss0,
        c31_mrna_intrinsic_loss_per_min=mrna_intrinsic_loss0,
    )


def repression(repressor, p, leak=0.0):
    r = max(float(repressor), 0.0)
    regulated = p.K_copies**p.hill_n / (p.K_copies**p.hill_n + r**p.hill_n)
    return leak + (1.0 - leak) * regulated


def occupancy(free, p):
    f = max(float(free), 0.0)
    if f == 0.0:
        return 0.0
    return f**p.hill_n / (p.K_copies**p.hill_n + f**p.hill_n)


def free_tetR(total, p):
    total = max(float(total), 0.0)
    if total == 0.0:
        return 0.0
    sites = p.effective_functional_sites + p.effective_sponge_sites
    def balance(free):
        return free + sites * occupancy(free, p) - total
    free = brentq(balance, 0.0, total, xtol=1e-12, rtol=1e-12, maxiter=100)
    residual = abs(balance(free)) / max(total, 1.0)
    if residual >= 1e-6 or not 0.0 <= free <= total:
        raise RuntimeError(f"TetR conservation failed: residual={residual:g}")
    return free


def promoter_activities(tet_total, ci, laci, p):
    tet_free = free_tetR(tet_total, p)
    return (
        repression(laci, p, p.oscillator_promoter_leak),
        repression(tet_free, p, p.oscillator_promoter_leak),
        repression(ci, p, p.oscillator_promoter_leak),
        repression(tet_free, p, p.c31_promoter_leak),
        tet_free,
    )


def rhs(_, y, p):
    m_tet, tet, m_ci, ci, m_laci, laci, m31, c31 = y
    h_tet, h_ci, h_laci, h31, _ = promoter_activities(tet, ci, laci, p)
    mrna_loss = p.oscillator_mrna_intrinsic_loss_per_min + p.mu
    protein_loss = p.oscillator_protein_intrinsic_loss_per_min + p.mu
    c31_mrna_loss = p.c31_mrna_intrinsic_loss_per_min + p.mu
    c31_protein_loss = p.c31_protein_intrinsic_loss_per_min + p.mu
    No = p.oscillator_plasmid_copies
    alpha = p.oscillator_tx_per_plasmid_per_min
    beta = p.oscillator_translation_per_mrna_per_min
    return (
        No * alpha * h_tet - mrna_loss * m_tet,
        beta * m_tet - protein_loss * tet,
        No * alpha * h_ci - mrna_loss * m_ci,
        beta * m_ci - protein_loss * ci,
        No * alpha * h_laci - mrna_loss * m_laci,
        beta * m_laci - protein_loss * laci,
        No * p.c31_tx_per_plasmid_per_min * h31 - c31_mrna_loss * m31,
        p.c31_translation_per_mrna_per_min * m31 - c31_protein_loss * c31,
    )


def initial_state(p):
    # Non-symmetric protein initial condition prevents the false symmetric state.
    tet, ci, laci = 100.0, 20.0, 300.0
    h_tet, h_ci, h_laci, h31, _ = promoter_activities(tet, ci, laci, p)
    mrna_loss = p.oscillator_mrna_intrinsic_loss_per_min + p.mu
    c31_mrna_loss = p.c31_mrna_intrinsic_loss_per_min + p.mu
    No = p.oscillator_plasmid_copies
    alpha = p.oscillator_tx_per_plasmid_per_min
    return (
        No * alpha * h_tet / mrna_loss, tet,
        No * alpha * h_ci / mrna_loss, ci,
        No * alpha * h_laci / mrna_loss, laci,
        No * p.c31_tx_per_plasmid_per_min * h31 / c31_mrna_loss, 0.0,
    )


def oscillation_metrics(t, signal, p):
    span = float(np.max(signal) - np.min(signal))
    scale = max(float(np.max(np.abs(signal))), 1.0)
    relative_amplitude = span / scale
    if relative_amplitude < SETTINGS["minimum_relative_amplitude"]:
        return {"status": "steady_state", "period_min": float("nan"),
                "period_generations": float("nan"), "period_cv": float("nan"),
                "relative_amplitude": relative_amplitude,
                "peaks": np.array([], dtype=int)}
    minimum_distance = max(1, int(round(2.0 * p.Td_min / SETTINGS["sample_step_min"])))
    peaks, _ = find_peaks(signal, prominence=max(0.05 * span, 1e-9),
                          distance=minimum_distance)
    if len(peaks) < SETTINGS["minimum_peaks"]:
        return {"status": "unresolved_oscillation", "period_min": float("nan"),
                "period_generations": float("nan"), "period_cv": float("nan"),
                "relative_amplitude": relative_amplitude, "peaks": peaks}
    periods = np.diff(t[peaks])
    period = float(np.mean(periods))
    return {"status": "resolved_oscillation", "period_min": period,
            "period_generations": period / p.Td_min,
            "period_cv": (float(np.std(periods, ddof=1) / period)
                          if len(periods) > 1 else 0.0),
            "relative_amplitude": relative_amplitude, "peaks": peaks}


def median_peak_fwhm(t, signal, peaks):
    """Median absolute half-maximum width of the last several resolved pulses."""
    widths = []
    for peak in peaks[-6:]:
        level = 0.5 * float(signal[peak])
        left = int(peak)
        right = int(peak)
        while left > 0 and signal[left] >= level:
            left -= 1
        while right < len(signal) - 1 and signal[right] >= level:
            right += 1
        if left > 0 and right < len(signal) - 1:
            widths.append(float(t[right] - t[left]))
    return float(np.median(widths)) if widths else float("nan")


def last_cycle_dose(t, signal, peaks):
    if len(peaks) < 2:
        return float("nan")
    segment = slice(peaks[-2], peaks[-1] + 1)
    return float(trapezoid(signal[segment], x=t[segment]))


def safe_name(axis, value):
    return f"{axis}__{value:g}".replace(".", "p")


def main():
    base = make_reference_parameters()
    time = np.arange(0.0, SETTINGS["duration_min"] + 0.5 * SETTINGS["sample_step_min"],
                     SETTINGS["sample_step_min"])
    keep = time >= SETTINGS["transient_min"]
    t_ss = time[keep]
    output = OUTPUT_ROOT / f"Oscillator_Engineering_Sensitivity_Analysis_{datetime.now():%Y%m%d_%H%M%S}"
    output.mkdir(parents=True, exist_ok=False)
    
    conditions = []
    for axis, values in SCAN_AXES.items():
        for value in values:
            p = base
            if axis == "oscillator_promoter_scale":
                p = replace(p, oscillator_tx_per_plasmid_per_min=(
                    base.oscillator_tx_per_plasmid_per_min * value))
            elif axis == "oscillator_RBS_scale":
                p = replace(p, oscillator_translation_per_mrna_per_min=(
                    base.oscillator_translation_per_mrna_per_min * value))
            elif axis == "sponge_plasmid_copies":
                p = replace(p, sponge_plasmid_copies=value)
            elif axis == "Td_min":
                # Only mu changes. Intrinsic mRNA/protein losses and synthesis rates
                # remain physical constants calibrated at the 50-min reference.
                p = replace(p, Td_min=value)
            else:
                raise KeyError(axis)
            conditions.append((axis, float(value), p))
    
    rows = []
    for condition_index, (axis, value, p) in enumerate(conditions, start=1):
        scenario = safe_name(axis, value)
        print(f"[{condition_index:02d}/{len(conditions):02d}] Running {scenario}")
        sol = solve_ivp(lambda t, y: rhs(t, y, p), (time[0], time[-1]), initial_state(p),
                        t_eval=time, method="LSODA", rtol=1e-7, atol=1e-9)
        if not sol.success:
            raise RuntimeError(f"{scenario}: {sol.message}")
    
        tet_free = np.array([free_tetR(value_total, p) for value_total in sol.y[1]])
        activity = np.array([repression(value_free, p, p.c31_promoter_leak)
                             for value_free in tet_free])
        osc = oscillation_metrics(t_ss, tet_free[keep], p)
        c31 = sol.y[7]
        flux_copies_per_min = p.c31_translation_per_mrna_per_min * sol.y[6]
        flux_uM_per_h = (
            60.0 * flux_copies_per_min * SETTINGS["active_C31_fraction"]
            / (SETTINGS["molecules_per_uM_per_fL"] * SETTINGS["cell_volume_fL"])
        )
        flux_peaks, _ = find_peaks(
            flux_uM_per_h[keep],
            prominence=max(0.05 * float(np.ptp(flux_uM_per_h[keep])), 1e-9),
            distance=max(1, int(round(2.0 * p.Td_min / SETTINGS["sample_step_min"]))))
    
        row = {
            "scenario": scenario,
            "scan_axis": axis,
            "scan_value": value,
            "is_reference_value_on_axis": bool(math.isclose(value, REFERENCE[axis])),
            "status": osc["status"],
            "period_min": osc["period_min"],
            "period_generations": osc["period_generations"],
            "deterministic_period_cv": osc["period_cv"],
            "relative_amplitude": osc["relative_amplitude"],
            "TetR_free_peak": float(np.max(tet_free[keep])),
            "TetR_free_trough": float(np.min(tet_free[keep])),
            "PLtetO1_duty_activity_gt_0p5": float(np.mean(activity[keep] > 0.5)),
            "C31_protein_peak_copies": float(np.max(c31[keep])),
            "C31_protein_trough_copies": float(np.min(c31[keep])),
            "C31_translation_flux_peak_uM_h": float(np.max(flux_uM_per_h[keep])),
            "C31_translation_flux_FWHM_h": (
                median_peak_fwhm(t_ss, flux_uM_per_h[keep], flux_peaks) / 60.0),
            "C31_translation_flux_dose_uM_per_cycle": (
                last_cycle_dose(t_ss, flux_uM_per_h[keep], osc["peaks"]) / 60.0),
            "Td_min": p.Td_min,
            "mu_per_min": p.mu,
            "oscillator_tx_per_plasmid_per_min": p.oscillator_tx_per_plasmid_per_min,
            "oscillator_translation_per_mrna_per_min": p.oscillator_translation_per_mrna_per_min,
            "oscillator_mrna_total_loss_per_min": (
                p.oscillator_mrna_intrinsic_loss_per_min + p.mu),
            "oscillator_mrna_total_half_life_min": (
                math.log(2.0) / (p.oscillator_mrna_intrinsic_loss_per_min + p.mu)),
            "sponge_plasmid_copies": p.sponge_plasmid_copies,
            "effective_functional_sites": p.effective_functional_sites,
            "effective_sponge_sites": p.effective_sponge_sites,
        }
        rows.append(row)
    
        if SETTINGS["save_trajectories"]:
            folder = output / scenario
            folder.mkdir()
            with (folder / "trajectory.csv").open("w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle)
                writer.writerow([
                    "time_min", "time_generations", "TetR_mRNA", "TetR_total",
                    "CI_mRNA", "CI", "LacI_mRNA", "LacI", "C31_mRNA",
                    "C31_protein", "TetR_free", "PLtetO1_activity",
                    "C31_translation_flux_copies_min", "C31_translation_flux_uM_h",
                ])
                for index, t_value in enumerate(time):
                    writer.writerow([
                        t_value, t_value / p.Td_min, *sol.y[:, index],
                        tet_free[index], activity[index], flux_copies_per_min[index],
                        flux_uM_per_h[index],
                    ])
    
    # Compute percent changes relative to the reference value WITHIN each axis.
    for axis in SCAN_AXES:
        axis_rows = [row for row in rows if row["scan_axis"] == axis]
        reference_row = next(row for row in axis_rows if row["is_reference_value_on_axis"])
        reference_period = reference_row["period_generations"]
        reference_dose = reference_row["C31_translation_flux_dose_uM_per_cycle"]
        for row in axis_rows:
            row["period_change_vs_axis_reference_percent"] = (
                100.0 * (row["period_generations"] / reference_period - 1.0)
                if np.isfinite(reference_period) and np.isfinite(row["period_generations"])
                else float("nan"))
            row["C31_dose_change_vs_axis_reference_percent"] = (
                100.0 * (row["C31_translation_flux_dose_uM_per_cycle"] / reference_dose - 1.0)
                if np.isfinite(reference_dose)
                and np.isfinite(row["C31_translation_flux_dose_uM_per_cycle"])
                else float("nan"))
    
    metrics_path = output / "v53a_engineering_scan_metrics.csv"
    with metrics_path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    
    with (output / "v53a_parameters_and_scope.json").open("w", encoding="utf-8") as handle:
        json.dump({
            "version": VERSION,
            "settings": SETTINGS,
            "reference_design": REFERENCE,
            "scan_axes": SCAN_AXES,
            "frozen_reference_parameters": asdict(base),
            "reference_alpha_mapping": (
                "beta*No*alpha/(delta_m+mu)=(delta_p+mu)*lambda, applied once at Td=50 min"),
            "scan_rule": (
                "One factor at a time; alpha and beta are not QSS-remapped during scans; "
                "Td changes dilution only; module B and stochastic noise are excluded."),
        }, handle, ensure_ascii=False, indent=2)
    
    if SETTINGS["make_plots"]:
        labels = {
            "oscillator_promoter_scale": "oscillator promoter scale",
            "oscillator_RBS_scale": "oscillator RBS scale",
            "sponge_plasmid_copies": "sponge plasmid copies Nt",
            "Td_min": "doubling time Td (min)",
        }
        fig, axes = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)
        for ax, axis in zip(axes.flat, SCAN_AXES):
            axis_rows = [row for row in rows if row["scan_axis"] == axis]
            x = np.array([row["scan_value"] for row in axis_rows], dtype=float)
            period = np.array([row["period_generations"] for row in axis_rows], dtype=float)
            width = np.array([row["C31_translation_flux_FWHM_h"] for row in axis_rows], dtype=float)
            ax.plot(x, period, "o-", color="tab:blue", label="period")
            ax.axvline(REFERENCE[axis], color="0.5", linestyle="--", linewidth=1)
            ax.set_xlabel(labels[axis])
            ax.set_ylabel("period (generations)", color="tab:blue")
            ax.tick_params(axis="y", labelcolor="tab:blue")
            ax.grid(alpha=0.25)
            twin = ax.twinx()
            twin.plot(x, width, "s--", color="tab:orange", label="C31 flux FWHM")
            twin.set_ylabel("C31 flux FWHM (h)", color="tab:orange")
            twin.tick_params(axis="y", labelcolor="tab:orange")
            ax.set_title(axis)
        fig.suptitle("v53a: one-at-a-time engineering effects on clock and C31 output")
        fig.savefig(output / "v53a_period_and_C31_width.png", dpi=180)
        if notebook_display is not None:
            notebook_display(fig)
        plt.close(fig)
    
        fig, axes = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)
        for ax, axis in zip(axes.flat, SCAN_AXES):
            axis_rows = [row for row in rows if row["scan_axis"] == axis]
            x = np.array([row["scan_value"] for row in axis_rows], dtype=float)
            dose = np.array([row["C31_translation_flux_dose_uM_per_cycle"]
                             for row in axis_rows], dtype=float)
            duty = np.array([row["PLtetO1_duty_activity_gt_0p5"]
                             for row in axis_rows], dtype=float)
            ax.plot(x, dose, "o-", color="tab:green", label="C31 dose")
            ax.axvline(REFERENCE[axis], color="0.5", linestyle="--", linewidth=1)
            ax.set_xlabel(labels[axis])
            ax.set_ylabel("C31 flux dose (uM/cycle)", color="tab:green")
            ax.tick_params(axis="y", labelcolor="tab:green")
            ax.grid(alpha=0.25)
            twin = ax.twinx()
            twin.plot(x, duty, "d--", color="tab:purple", label="PLtetO1 duty")
            twin.set_ylabel("PLtetO1 duty fraction", color="tab:purple")
            twin.tick_params(axis="y", labelcolor="tab:purple")
            ax.set_title(axis)
        fig.suptitle("v53a: output dose and promoter duty cycle")
        fig.savefig(output / "v53a_C31_dose_and_duty.png", dpi=180)
        if notebook_display is not None:
            notebook_display(fig)
        plt.close(fig)
    
    print(f"Loaded diagnostic code: {VERSION}")
    print(f"Saved v53a outputs to: {output.resolve()}")
    for row in rows:
        print(
            f"{row['scenario']:>36s} | {row['status']:>22s} | "
            f"P={row['period_min']:.1f} min/{row['period_generations']:.3f} gen | "
            f"dP={row['period_change_vs_axis_reference_percent']:+.2f}% | "
            f"C31 width={row['C31_translation_flux_FWHM_h']:.2f} h | "
            f"dose={row['C31_translation_flux_dose_uM_per_cycle']:.2f} uM/cycle"
        )
    
    print("\nInterpretation boundaries:")
    print("- This is an OAT engineering sensitivity scan, not parameter fitting or GSA.")
    print("- Select parameters from evidence/design constraints, not from the prettiest curve.")
    print("- Deterministic period CV is numerical regularity, not biological trajectory CV.")
    print("- B-module compatibility and stochastic sponge noise reduction are separate tests.")


if __name__ == "__main__":
    main()
