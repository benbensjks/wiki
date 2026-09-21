# Standalone program: flux-driven finite shared-translation-resource model.
# Source version v53d.
#
# The frozen v52a1 oscillator is retained.  Only a shared ribosome/resource
# availability factor is added to protein translation.  The factor is driven
# by the instantaneous, unthrottled C31 translation demand beta31*m31:
#
#   q31(t) = beta31*m31(t) / J31_reference_peak
#   b       = L_peak / (1-L_peak)
#   rho(t)  = 1 / (1 + b*q31(t))
#
# Thus rho=1-L_peak when q31=1.  rho is positive and saturating, unlike the
# older linear boundary rho=1-L*H31.  L_peak is still an uncertainty/boundary
# parameter until expression-burden experiments are available.

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
    from IPython.display import display
except ImportError:
    display = None

VERSION = "v53d_flux_driven_shared_translation_resource_20260814"

# ============================== USER SETTINGS ==============================
SETTINGS = {
    "duration_min": 6000.0,
    "sample_step_min": 1.0,
    "transient_min": 1500.0,
    # Peak fractional loss of translation availability at the reference C31
    # demand.  These are boundary scenarios, not fitted properties of B.
    "peak_load_fractions": (0.00, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50),
    "timebase_tolerance_percent": 5.0,
    "amplitude_tolerance_percent": 10.0,
    "minimum_relative_amplitude": 0.01,
    "minimum_peaks": 4,
    "cell_volume_fL": 1.0,
    "active_C31_fraction": 1.0,
    "molecules_per_uM_per_fL": 602.214076,
    "save_trajectories": True,
    "make_plots": True,
}
# ============================================================================


@dataclass(frozen=True)
class Parameters:
    # Frozen v52a1 design.
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
    oscillator_tx_per_plasmid_per_min: float = 0.0
    oscillator_translation_per_mrna_per_min: float = 0.5
    oscillator_mrna_intrinsic_loss_per_min: float = 0.0
    oscillator_protein_intrinsic_loss_per_min: float = 0.0
    c31_promoter_leak: float = 0.005
    c31_tx_per_plasmid_per_min: float = 5.0
    c31_translation_per_mrna_per_min: float = 0.5
    c31_mrna_intrinsic_loss_per_min: float = 0.0
    c31_protein_intrinsic_loss_per_min: float = 0.0
    peak_load_fraction: float = 0.0

    @property
    def mu(self):
        return math.log(2.0) / self.Td_min

    @property
    def effective_sites(self):
        return (
            self.functional_accessible_fraction
            * self.functional_tetO_per_plasmid
            * self.oscillator_plasmid_copies
            + self.sponge_accessible_fraction
            * self.sponge_tetO_per_plasmid
            * self.sponge_plasmid_copies
        )


def make_reference_parameters():
    mu = math.log(2.0) / 50.0
    total_mrna_loss = math.log(2.0) / 2.0
    intrinsic_mrna_loss = total_mrna_loss - mu
    beta = 0.5
    alpha = mu * 1000.0 * total_mrna_loss / (beta * 10.0)
    return Parameters(
        oscillator_tx_per_plasmid_per_min=alpha,
        oscillator_mrna_intrinsic_loss_per_min=intrinsic_mrna_loss,
        c31_mrna_intrinsic_loss_per_min=intrinsic_mrna_loss,
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

    def balance(free):
        return free + p.effective_sites * occupancy(free, p) - total

    free = brentq(balance, 0.0, total, xtol=1e-12, rtol=1e-12)
    if abs(balance(free)) / max(total, 1.0) >= 1e-6:
        raise RuntimeError("TetR conservation failed")
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


def resource_availability(m31, p, reference_C31_demand_peak):
    """Finite resource-pool approximation driven by requested C31 translation."""
    load = float(p.peak_load_fraction)
    if not 0.0 <= load < 1.0:
        raise ValueError("peak_load_fraction must be in [0,1)")
    if load == 0.0:
        return 1.0
    unthrottled_demand = p.c31_translation_per_mrna_per_min * max(float(m31), 0.0)
    q31 = unthrottled_demand / reference_C31_demand_peak
    burden_ratio = load / (1.0 - load)
    return 1.0 / (1.0 + burden_ratio * q31)


def rhs(_, y, p, reference_C31_demand_peak):
    m_tet, tet, m_ci, ci, m_laci, laci, m31, c31 = y
    h_tet, h_ci, h_laci, h31, _ = promoter_activities(tet, ci, laci, p)
    rho = resource_availability(m31, p, reference_C31_demand_peak)
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


def initial_state(p):
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


def simulate(p, time, reference_C31_demand_peak):
    sol = solve_ivp(
        lambda t, y: rhs(t, y, p, reference_C31_demand_peak),
        (float(time[0]), float(time[-1])), initial_state(p), t_eval=time,
        method="LSODA", rtol=1e-7, atol=1e-9,
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    tet_free = np.array([free_tetR(total, p) for total in sol.y[1]])
    h31 = np.array([repression(free, p, p.c31_promoter_leak) for free in tet_free])
    resource = np.array([
        resource_availability(m31, p, reference_C31_demand_peak)
        for m31 in sol.y[6]
    ])
    demand = p.c31_translation_per_mrna_per_min * sol.y[6]
    actual_flux_copies_min = resource * demand
    actual_flux_uM_h = (
        60.0 * actual_flux_copies_min * SETTINGS["active_C31_fraction"]
        / (SETTINGS["molecules_per_uM_per_fL"] * SETTINGS["cell_volume_fL"])
    )
    return sol, tet_free, h31, resource, demand, actual_flux_uM_h


def oscillation_metrics(t, signal, p):
    span = float(np.ptp(signal))
    relative_amplitude = span / max(float(np.max(np.abs(signal))), 1.0)
    distance = max(1, int(round(2.0 * p.Td_min / SETTINGS["sample_step_min"])))
    peaks, _ = find_peaks(
        signal, prominence=max(0.05 * span, 1e-9), distance=distance
    )
    if relative_amplitude < SETTINGS["minimum_relative_amplitude"]:
        return "steady_state", float("nan"), float("nan"), peaks
    if len(peaks) < SETTINGS["minimum_peaks"]:
        return "unresolved_oscillation", float("nan"), float("nan"), peaks
    periods = np.diff(t[peaks])
    period = float(np.mean(periods))
    cv = float(np.std(periods, ddof=1) / period) if len(periods) > 1 else 0.0
    return "resolved_oscillation", period, cv, peaks


def median_fwhm(t, signal, peaks):
    widths = []
    for peak in peaks[-6:]:
        level = 0.5 * float(signal[peak])
        left = right = int(peak)
        while left > 0 and signal[left] >= level:
            left -= 1
        while right < len(signal)-1 and signal[right] >= level:
            right += 1
        if left > 0 and right < len(signal)-1:
            widths.append(float(t[right]-t[left]))
    return float(np.median(widths)) if widths else float("nan")


def last_cycle_dose(t, signal, peaks):
    if len(peaks) < 2:
        return float("nan")
    segment = slice(peaks[-2], peaks[-1]+1)
    return float(trapezoid(signal[segment], x=t[segment]) / 60.0)


def write_rows(path, rows):
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def main():
    base = make_reference_parameters()
    time = np.arange(
        0.0, SETTINGS["duration_min"] + 0.5*SETTINGS["sample_step_min"],
        SETTINGS["sample_step_min"],
    )
    keep = time >= SETTINGS["transient_min"]
    t_ss = time[keep]
    
    # First obtain the physical reference peak demand without resource feedback.
    # A dummy positive denominator is safe because load=0 makes rho exactly 1.
    reference_sol, reference_tet_free, _, _, reference_demand, _ = simulate(
        base, time, 1.0
    )
    reference_C31_demand_peak = float(np.max(reference_demand[keep]))
    if reference_C31_demand_peak <= 0.0:
        raise RuntimeError("Reference C31 translation demand peak was not positive")
    
    output = OUTPUT_ROOT / (
        "Flux_Driven_Translation_Burden_Model_" + datetime.now().strftime("%Y%m%d_%H%M%S")
    )
    output.mkdir(parents=True, exist_ok=False)
    
    rows = []
    saved = {}
    loads = SETTINGS["peak_load_fractions"]
    for index, load in enumerate(loads, start=1):
        scenario = f"peak_resource_load_{100*load:g}percent".replace(".", "p")
        print(f"[{index:02d}/{len(loads):02d}] Running {scenario}")
        p = replace(base, peak_load_fraction=float(load))
        sol, tet_free, h31, resource, demand, flux = simulate(
            p, time, reference_C31_demand_peak
        )
        status, period, period_cv, osc_peaks = oscillation_metrics(
            t_ss, tet_free[keep], p
        )
        flux_peaks, _ = find_peaks(
            flux[keep], prominence=max(0.05*float(np.ptp(flux[keep])), 1e-9),
            distance=max(1, int(round(2*p.Td_min/SETTINGS["sample_step_min"]))),
        )
        row = {
            "scenario": scenario,
            "peak_load_fraction": load,
            "burden_ratio_L_over_1_minus_L": (load/(1-load) if load < 1 else float("inf")),
            "status": status,
            "period_min": period,
            "period_generations": period/p.Td_min if np.isfinite(period) else float("nan"),
            "deterministic_period_cv": period_cv,
            "resource_min": float(np.min(resource[keep])),
            "resource_mean": float(np.mean(resource[keep])),
            "unthrottled_C31_demand_peak_molecules_min": float(np.max(demand[keep])),
            "actual_C31_flux_peak_uM_h": float(np.max(flux[keep])),
            "actual_C31_flux_FWHM_h": median_fwhm(t_ss, flux[keep], flux_peaks)/60.0,
            "actual_C31_flux_dose_uM_cycle": last_cycle_dose(
                t_ss, flux[keep], osc_peaks
            ),
            "TetR_free_peak": float(np.max(tet_free[keep])),
            "TetR_free_trough": float(np.min(tet_free[keep])),
            "TetR_total_peak": float(np.max(sol.y[1, keep])),
            "CI_peak": float(np.max(sol.y[3, keep])),
            "LacI_peak": float(np.max(sol.y[5, keep])),
            "C31_peak_copies": float(np.max(sol.y[7, keep])),
            "C31_trough_copies": float(np.min(sol.y[7, keep])),
            "PLtetO1_duty_gt_0p5": float(np.mean(h31[keep] > 0.5)),
        }
        rows.append(row)
        saved[scenario] = (p, sol, tet_free, h31, resource, demand, flux)
    
        if SETTINGS["save_trajectories"]:
            folder = output / scenario
            folder.mkdir()
            with (folder / "trajectory.csv").open(
                "w", newline="", encoding="utf-8-sig"
            ) as handle:
                writer = csv.writer(handle)
                writer.writerow([
                    "time_min", "time_generations", "TetR_mRNA", "TetR_total",
                    "CI_mRNA", "CI", "LacI_mRNA", "LacI", "C31_mRNA",
                    "C31_protein", "TetR_free", "PLtetO1_activity",
                    "resource_availability", "unthrottled_C31_demand_molecules_min",
                    "actual_C31_translation_flux_uM_h",
                ])
                for j, t_value in enumerate(time):
                    writer.writerow([
                        t_value, t_value/p.Td_min, *sol.y[:, j], tet_free[j], h31[j],
                        resource[j], demand[j], flux[j],
                    ])
    
    reference_row = next(row for row in rows if row["peak_load_fraction"] == 0.0)
    for row in rows:
        row["period_change_percent"] = (
            100.0*(row["period_generations"]/reference_row["period_generations"]-1.0)
            if np.isfinite(row["period_generations"]) else float("nan")
        )
        row["TetR_peak_change_percent"] = (
            100.0*(row["TetR_free_peak"]/reference_row["TetR_free_peak"]-1.0)
        )
        row["C31_dose_change_percent"] = (
            100.0*(row["actual_C31_flux_dose_uM_cycle"]
                   / reference_row["actual_C31_flux_dose_uM_cycle"]-1.0)
        )
        row["timebase_stable"] = bool(
            row["status"] == "resolved_oscillation"
            and abs(row["period_change_percent"])
            <= SETTINGS["timebase_tolerance_percent"]
        )
        row["amplitude_stable"] = bool(
            abs(row["TetR_peak_change_percent"])
            <= SETTINGS["amplitude_tolerance_percent"]
        )
        row["engineering_pass"] = bool(
            row["timebase_stable"] and row["amplitude_stable"]
        )
    
    write_rows(output / "v53d_flux_driven_load_metrics.csv", rows)
    with (output / "v53d_scope_and_parameters.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump({
            "version": VERSION,
            "settings": SETTINGS,
            "frozen_parameters": asdict(base),
            "reference_C31_demand_peak_molecules_min": reference_C31_demand_peak,
            "resource_equation": "q31=beta31*m31/Jref; b=L/(1-L); rho=1/(1+b*q31)",
            "translation_equation": "all protein synthesis beta_i*m_i is multiplied by rho(t)",
            "not_included_to_avoid_double_counting": [
                "extra TetO pool beyond the already audited 4*No plus sponge sites",
                "Int-att sequestration, which belongs to the Zhao B ODE",
                "growth-rate burden, pending construct-specific Td measurements",
                "explicit RNAP competition and protease saturation without calibration data",
            ],
            "interpretation_limit": (
                "L_peak is a boundary parameter until resource usage is calibrated "
                "with expression and growth measurements"
            ),
        }, handle, ensure_ascii=False, indent=2)
    
    if SETTINGS["make_plots"]:
        x = 100*np.array([r["peak_load_fraction"] for r in rows])
        fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
        axes[0, 0].plot(x, [r["period_change_percent"] for r in rows], "o-")
        tol = SETTINGS["timebase_tolerance_percent"]
        axes[0, 0].axhspan(-tol, tol, color="green", alpha=0.15)
        axes[0, 0].set_ylabel("period change (%)")
        axes[0, 1].plot(x, [r["TetR_peak_change_percent"] for r in rows], "o-",
                        label="TetR peak")
        atol = SETTINGS["amplitude_tolerance_percent"]
        axes[0, 1].axhspan(-atol, atol, color="green", alpha=0.15)
        axes[0, 1].set_ylabel("TetR peak change (%)")
        axes[1, 0].plot(x, [r["resource_min"] for r in rows], "o-", label="minimum")
        axes[1, 0].plot(x, [r["resource_mean"] for r in rows], "s--", label="mean")
        axes[1, 0].set_ylabel("translation resource available")
        axes[1, 0].legend()
        axes[1, 1].plot(x, [r["C31_dose_change_percent"] for r in rows], "o-")
        axes[1, 1].set_ylabel("C31 dose change (%)")
        for axis in axes.flat:
            axis.set_xlabel("peak resource load at reference C31 demand (%)")
            axis.grid(alpha=0.25)
        fig.suptitle("v53d flux-driven shared translation-resource load")
        fig.savefig(output / "v53d_load_boundaries.png", dpi=180)
        if display is not None:
            display(fig)
        plt.close(fig)
    
        fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True,
                                 constrained_layout=True)
        for load, color in ((0.0, "black"), (0.10, "tab:orange"),
                            (0.30, "tab:red"), (0.50, "tab:purple")):
            scenario = f"peak_resource_load_{100*load:g}percent".replace(".", "p")
            if scenario not in saved:
                continue
            _, sol, tet_free, _, resource, _, _ = saved[scenario]
            axes[0].plot(time/60.0, tet_free, color=color,
                         label=f"peak load={100*load:.0f}%")
            axes[1].plot(time/60.0, resource, color=color,
                         label=f"peak load={100*load:.0f}%")
        axes[0].set_xlim(30, 70)
        axes[0].set_ylabel("free TetR (copies/cell)")
        axes[0].set_title("Clock waveform under dynamic translation burden")
        axes[1].set_ylabel("resource availability rho(t)")
        axes[1].set_xlabel("time (h)")
        for axis in axes:
            axis.legend(); axis.grid(alpha=0.25)
        fig.savefig(output / "v53d_load_trajectories.png", dpi=180)
        if display is not None:
            display(fig)
        plt.close(fig)
    
    print("Loaded diagnostic code:", VERSION)
    print("Saved v53d outputs to:", output.resolve())
    print(
        "Resource model: q31=beta31*m31/Jref, b=L/(1-L), "
        "rho=1/(1+b*q31); all protein translation terms are multiplied by rho."
    )
    for row in rows:
        print(
            f"peak load={100*row['peak_load_fraction']:>4.0f}% | "
            f"{row['status']:>22s} | "
            f"P={row['period_min']:.1f} min/{row['period_generations']:.3f} gen | "
            f"dP={row['period_change_percent']:+.2f}% | "
            f"TetR peak d={row['TetR_peak_change_percent']:+.1f}% | "
            f"resource(min/mean)={row['resource_min']:.3f}/{row['resource_mean']:.3f} | "
            f"C31 dose d={row['C31_dose_change_percent']:+.1f}% | "
            f"pass={row['engineering_pass']}"
        )
    
    passing = [r for r in rows if r["engineering_pass"]]
    if passing:
        largest = max(passing, key=lambda r: r["peak_load_fraction"])
        print(
            "Largest scanned load passing the provisional period/amplitude criteria: "
            f"{100*largest['peak_load_fraction']:.0f}%"
        )
    else:
        print("No scanned load passed the provisional period/amplitude criteria.")
    
    print("\nInterpretation boundaries:")
    print("- This is more physical than the old linear H31-gated boundary, but L_peak is still uncalibrated.")
    print("- Int-att sequestration remains in B and is not duplicated in A.")
    print("- Td stays fixed at 50 min; construct-specific growth burden must be tested separately when measured.")
    print("- Do not call a load fraction 'the real B burden' until expression/growth data calibrate it.")


if __name__ == "__main__":
    main()
