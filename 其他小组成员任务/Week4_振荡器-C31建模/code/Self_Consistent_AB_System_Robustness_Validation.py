# Standalone program: self-consistent translation burden and long A/B audit.
# Source version v53h with fixed PLtetO1.
#
# Frozen biology:
# - C31 and cI use the same PLtetO1 transfer function.
# - oscillator equations, TetO accounting, sponge and Td=50 min are unchanged.
# - C31 RBS changes beta31 only; Int tag-mediated loss is owned by module B.
#
# New validation:
# - A is re-solved for every RBS/load pair (not scaled from an RBS=1 trace).
# - resource availability is driven by the actual requested beta31*m31 flux.
# - the resulting rho*beta31*m31 flux is passed directly to Zhao module B.
# - PB and LR initial DNA states are both checked over >=10 steady cycles.

from dataclasses import dataclass, replace, asdict
from datetime import datetime
from pathlib import Path
import csv, json, math, sys

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp, trapezoid
from scipy.interpolate import CubicSpline
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

VERSION = "v53h_self_consistent_load_long_AB_validation_20260815"

# ============================== USER SETTINGS ==============================
SETTINGS = {
    # "main" = 1 candidate x 4 loads x 2 DNA initial states = 8 B solves.
    # "full" = 4 candidates x 4 loads x 2 DNA initial states = 32 B solves.
    "candidate_panel": "main",
    "main_candidate": (0.60, 12.0),       # (C31 RBS scale, k_tag_int h^-1)
    "full_candidates": (
        (0.45, 8.0), (0.60, 8.0), (0.60, 12.0), (0.70, 16.0),
    ),
    "peak_load_fractions": (0.00, 0.10, 0.20, 0.30),
    "DNA_initial_states": ("PB", "LR"),

    # Long enough for about 10-12 audited cycles after the transient.
    "duration_h": 160.0,
    "A_sample_step_min": 2.0,
    "A_transient_h": 30.0,
    "minimum_peak_distance_h": 5.0,
    "minimum_audited_cycles": 10,

    "B_sample_step_h": 0.04,
    "B_max_step_h": 0.08,
    "DNA_midpoint": 0.5,
    "complete_state_threshold": 0.95,
    "success_target": 0.95,
    "alternation_target": 0.95,
    "DNA_conservation_tolerance": 1e-8,
    "A_period_tolerance_percent": 5.0,
    "A_TetR_peak_tolerance_percent": 10.0,

    "expected_Td_min": 50.0,
    "cell_volume_fL": 1.0,
    "active_C31_fraction": 1.0,
    "molecules_per_uM_per_fL": 602.214076,
    "make_plots": True,
    "save_A_trajectories": True,
}
# ===========================================================================


@dataclass(frozen=True)
class AParameters:
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
    return AParameters(
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
    residual = abs(balance(free)) / max(total, 1.0)
    if residual >= 1e-6:
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


def resource_availability(m31, p, reference_demand_peak):
    load = float(p.peak_load_fraction)
    if not 0.0 <= load < 1.0:
        raise ValueError("peak_load_fraction must be in [0,1)")
    if load == 0.0:
        return 1.0
    requested = p.c31_translation_per_mrna_per_min * max(float(m31), 0.0)
    q31 = requested / reference_demand_peak
    burden_ratio = load / (1.0 - load)
    return 1.0 / (1.0 + burden_ratio*q31)


def a_rhs(_, y, p, reference_demand_peak):
    m_tet, tet, m_ci, ci, m_laci, laci, m31, c31 = y
    h_tet, h_ci, h_laci, h31, _ = promoter_activities(tet, ci, laci, p)
    rho = resource_availability(m31, p, reference_demand_peak)
    mrna_loss = p.oscillator_mrna_intrinsic_loss_per_min + p.mu
    protein_loss = p.oscillator_protein_intrinsic_loss_per_min + p.mu
    c31_mrna_loss = p.c31_mrna_intrinsic_loss_per_min + p.mu
    c31_protein_loss = p.c31_protein_intrinsic_loss_per_min + p.mu
    No = p.oscillator_plasmid_copies
    alpha = p.oscillator_tx_per_plasmid_per_min
    beta = p.oscillator_translation_per_mrna_per_min
    return (
        No*alpha*h_tet - mrna_loss*m_tet,
        rho*beta*m_tet - protein_loss*tet,
        No*alpha*h_ci - mrna_loss*m_ci,
        rho*beta*m_ci - protein_loss*ci,
        No*alpha*h_laci - mrna_loss*m_laci,
        rho*beta*m_laci - protein_loss*laci,
        No*p.c31_tx_per_plasmid_per_min*h31 - c31_mrna_loss*m31,
        rho*p.c31_translation_per_mrna_per_min*m31 - c31_protein_loss*c31,
    )


def a_initial_state(p):
    tet, ci, laci = 100.0, 20.0, 300.0
    h_tet, h_ci, h_laci, h31, _ = promoter_activities(tet, ci, laci, p)
    mrna_loss = p.oscillator_mrna_intrinsic_loss_per_min + p.mu
    c31_mrna_loss = p.c31_mrna_intrinsic_loss_per_min + p.mu
    No = p.oscillator_plasmid_copies
    alpha = p.oscillator_tx_per_plasmid_per_min
    return (
        No*alpha*h_tet/mrna_loss, tet,
        No*alpha*h_ci/mrna_loss, ci,
        No*alpha*h_laci/mrna_loss, laci,
        No*p.c31_tx_per_plasmid_per_min*h31/c31_mrna_loss, 0.0,
    )


def simulate_a(p, time_min, reference_demand_peak):
    sol = solve_ivp(
        lambda t, y: a_rhs(t, y, p, reference_demand_peak),
        (float(time_min[0]), float(time_min[-1])), a_initial_state(p),
        t_eval=time_min, method="LSODA", rtol=1e-7, atol=1e-9,
    )
    if not sol.success:
        raise RuntimeError("A solve failed: " + sol.message)
    tet_free = np.array([free_tetR(total, p) for total in sol.y[1]])
    resource = np.array([
        resource_availability(m31, p, reference_demand_peak) for m31 in sol.y[6]
    ])
    requested = p.c31_translation_per_mrna_per_min * sol.y[6]
    actual_copies_min = resource * requested
    actual_uM_h = (
        60.0 * actual_copies_min * SETTINGS["active_C31_fraction"]
        / (SETTINGS["molecules_per_uM_per_fL"] * SETTINGS["cell_volume_fL"])
    )
    return sol, tet_free, resource, requested, actual_uM_h


def peak_data(time_h, signal):
    keep = time_h >= SETTINGS["A_transient_h"]
    indices = np.flatnonzero(keep)
    y = signal[keep]
    dt = float(np.median(np.diff(time_h)))
    span = float(np.ptp(y))
    peaks, _ = find_peaks(
        y, prominence=max(0.05*span, 1e-9),
        distance=max(1, int(round(SETTINGS["minimum_peak_distance_h"]/dt))),
    )
    absolute_peaks = indices[peaks]
    times = time_h[absolute_peaks]
    if len(times) < SETTINGS["minimum_audited_cycles"] + 1:
        raise RuntimeError(
            f"Only {len(times)-1} post-transient cycles; increase duration_h"
        )
    period = float(np.mean(np.diff(times[-10:])))
    boundaries = 0.5*(times[:-1] + times[1:])
    return absolute_peaks, times, boundaries, period


def cycle_flux_metrics(time_h, flux, boundaries):
    rows = []
    for index, (left, right) in enumerate(zip(boundaries[:-1], boundaries[1:])):
        choose = (time_h >= left) & (time_h < right)
        t = time_h[choose]; y = flux[choose]
        if len(t) < 3:
            continue
        low, high = float(np.min(y)), float(np.max(y))
        level = low + 0.5*(high-low)
        width = float(np.sum(y >= level)*np.median(np.diff(t)))
        rows.append({
            "cycle_index": index,
            "cycle_left_h": float(left), "cycle_right_h": float(right),
            "flux_peak_uM_h": high,
            "flux_trough_uM_h": low,
            "flux_FWHM_h": width,
            "flux_dose_uM_per_cycle": float(trapezoid(y, x=t)),
        })
    if len(rows) < SETTINGS["minimum_audited_cycles"]:
        raise RuntimeError(f"Only {len(rows)} complete A cycles were audited")
    return rows[-SETTINGS["minimum_audited_cycles"]:]


def crossing_events(t, lr_fraction, left, right):
    local = np.flatnonzero((t >= left) & (t < right))
    events = []
    for i0, i1 in zip(local[:-1], local[1:]):
        a = lr_fraction[i0] - SETTINGS["DNA_midpoint"]
        b = lr_fraction[i1] - SETTINGS["DNA_midpoint"]
        if a*b < 0.0:
            fraction = abs(a)/(abs(a)+abs(b))
            event_time = t[i0] + fraction*(t[i1]-t[i0])
            events.append((event_time, "PB_to_LR" if b > a else "LR_to_PB"))
    return events


def sample_at(t, y, when):
    return float(np.interp(float(when), t, y))


def b_cycle_metrics(t, lr, pb, flux_spline, boundaries):
    rows = []
    complete = SETTINGS["complete_state_threshold"]
    audit_boundaries = boundaries[-(SETTINGS["minimum_audited_cycles"]+1):]
    for cycle_index, (left, right) in enumerate(
        zip(audit_boundaries[:-1], audit_boundaries[1:])
    ):
        events = crossing_events(t, lr, left, right)
        start_lr = sample_at(t, lr, left)
        end_lr = sample_at(t, lr, right)
        direction = events[0][1] if len(events) == 1 else ""
        endpoint_complete = (
            (direction == "PB_to_LR" and end_lr >= complete)
            or (direction == "LR_to_PB" and end_lr <= 1.0-complete)
        )
        one_complete = len(events) == 1 and endpoint_complete
        grid = np.arange(left, right+0.002, 0.002)
        dose = float(trapezoid(np.maximum(0.0, flux_spline(grid)), x=grid))
        rows.append({
            "cycle_index": cycle_index,
            "cycle_left_h": float(left), "cycle_right_h": float(right),
            "start_LR_fraction": start_lr, "end_LR_fraction": end_lr,
            "total_crossings": len(events),
            "single_crossing_direction": direction,
            "endpoint_complete": bool(endpoint_complete),
            "one_crossing_and_complete": bool(one_complete),
            "cycle_input_dose_uM": dose,
        })
    return rows


def alternation_fidelity(rows):
    eligible = []
    for previous, current in zip(rows[:-1], rows[1:]):
        if previous["one_crossing_and_complete"] and current["one_crossing_and_complete"]:
            eligible.append(
                previous["single_crossing_direction"]
                != current["single_crossing_direction"]
            )
    return float(np.mean(eligible)) if eligible else float("nan")


def write_rows(path, rows):
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def main():
    # ------------------------- Locate and validate B model ----------------------
    MODEL_DIR = SCRIPT_DIR
    if not (MODEL_DIR/"zhao_core.py").is_file():
        raise FileNotFoundError(
            "Run this notebook from the model folder containing zhao_core.py"
        )
    if str(MODEL_DIR) not in sys.path:
        sys.path.insert(0, str(MODEL_DIR))
    global zc
    import zhao_core as zc
    
    probe = zc.default_params()
    expected_k_dil = math.log(2.0)/(SETTINGS["expected_Td_min"]/60.0)
    if not math.isclose(probe["k_dil"], expected_k_dil, rel_tol=2e-4, abs_tol=2e-4):
        raise RuntimeError(
            f"A/B growth mismatch: B k_dil={probe['k_dil']:.6g} h^-1; "
            f"Td={SETTINGS['expected_Td_min']:g} min requires {expected_k_dil:.6g} h^-1"
        )
    
    panel = SETTINGS["candidate_panel"].lower()
    if panel == "main":
        candidates = (tuple(SETTINGS["main_candidate"]),)
    elif panel == "full":
        candidates = tuple(tuple(item) for item in SETTINGS["full_candidates"])
    else:
        raise ValueError("candidate_panel must be 'main' or 'full'")
    
    # ----------------------- Establish shared reference demand -----------------
    base = make_reference_parameters()
    time_min = np.arange(
        0.0, 60.0*SETTINGS["duration_h"] + 0.5*SETTINGS["A_sample_step_min"],
        SETTINGS["A_sample_step_min"],
    )
    time_h = time_min/60.0
    
    # The denominator is the RBS=1, no-load requested C31 peak, identical to v53d.
    reference_sol, _, _, reference_requested, _ = simulate_a(base, time_min, 1.0)
    reference_keep = time_h >= SETTINGS["A_transient_h"]
    reference_demand_peak = float(np.max(reference_requested[reference_keep]))
    if reference_demand_peak <= 0.0:
        raise RuntimeError("Reference C31 demand peak is not positive")
    
    output = OUTPUT_ROOT / (
        "Self_Consistent_AB_System_Robustness_Validation_" + datetime.now().strftime("%Y%m%d_%H%M%S")
    )
    output.mkdir(parents=True, exist_ok=False)
    
    summary_rows, a_cycle_rows, b_cycle_rows, design_rows = [], [], [], []
    saved_examples = {}
    conditions = [
        (float(scale), float(tag), float(load))
        for scale, tag in candidates
        for load in SETTINGS["peak_load_fractions"]
    ]
    total_b_solves = len(conditions)*len(SETTINGS["DNA_initial_states"])
    b_counter = 0
    
    for condition_number, (scale, tag, load) in enumerate(conditions, start=1):
        label = f"RBS_{scale:g}__tag_{tag:g}__load_{100*load:g}pct".replace(".", "p")
        print(f"[A {condition_number:02d}/{len(conditions):02d}] Running {label}")
        p = replace(
            base,
            c31_translation_per_mrna_per_min=(
                base.c31_translation_per_mrna_per_min*scale
            ),
            peak_load_fraction=load,
        )
        a_sol, tet_free, rho, requested, flux = simulate_a(
            p, time_min, reference_demand_peak
        )
        _, peak_times, boundaries, period_h = peak_data(time_h, flux)
        tet_peak_indices, _, _, _ = peak_data(time_h, tet_free)
        tetR_peak_copies = float(np.median(tet_free[tet_peak_indices[-10:]]))
        a_cycles = cycle_flux_metrics(time_h, flux, boundaries)
        for row in a_cycles:
            a_cycle_rows.append({"condition": label, **row})
    
        spline = CubicSpline(time_h, flux, bc_type="natural")
    
        if SETTINGS["save_A_trajectories"]:
            condition_dir = output/label
            condition_dir.mkdir(parents=True, exist_ok=False)
            a_rows = []
            for i in range(len(time_min)):
                a_rows.append({
                    "time_min": float(time_min[i]),
                    "TetR_total_copies": float(a_sol.y[1, i]),
                    "TetR_free_copies": float(tet_free[i]),
                    "CI_copies": float(a_sol.y[3, i]),
                    "LacI_copies": float(a_sol.y[5, i]),
                    "C31_mRNA_copies": float(a_sol.y[6, i]),
                    "C31_protein_shadow_copies": float(a_sol.y[7, i]),
                    "translation_resource_rho": float(rho[i]),
                    "requested_C31_translation_copies_min": float(requested[i]),
                    "actual_C31_translation_flux_uM_h": float(flux[i]),
                })
            write_rows(condition_dir/"trajectory.csv", a_rows)
    
        def int_production(t):
            return max(0.0, float(spline(float(t))))
    
        for initial_name in SETTINGS["DNA_initial_states"]:
            b_counter += 1
            print(
                f"    [B {b_counter:02d}/{total_b_solves:02d}] "
                f"initial={initial_name}"
            )
            P = zc.default_params()
            P.update(
                K_rep=0.0186, n_rep=3.4, krep_tsl=15.0,
                krdf_tsl=200.0, k_tag_int=tag,
            )
            C = zc._rate_constants(P)
            if initial_name == "PB":
                y0 = zc.y0_PB_ss(P)
            elif initial_name == "LR":
                y0 = zc.y0_LR_ss(P)
            else:
                raise ValueError(initial_name)
    
            b_time = np.arange(
                0.0, SETTINGS["duration_h"] + 0.5*SETTINGS["B_sample_step_h"],
                SETTINGS["B_sample_step_h"],
            )
            b_sol = solve_ivp(
                lambda t, y: zc.rhs(t, y, P, C, int_production),
                (float(b_time[0]), float(b_time[-1])), y0, t_eval=b_time,
                method="LSODA", rtol=1e-7, atol=1e-9,
                max_step=SETTINGS["B_max_step_h"],
            )
            if not b_sol.success:
                raise RuntimeError("B solve failed: " + b_sol.message)
    
            states = b_sol.y.T
            lr = zc.LR_total(states)/P["Dtot"]
            pb = zc.PB_total(states)/P["Dtot"]
            conservation = float(np.max(np.abs(lr+pb-1.0)))
            cycles = b_cycle_metrics(b_time, lr, pb, spline, boundaries)
            for row in cycles:
                b_cycle_rows.append({
                    "condition": label, "initial_DNA_state": initial_name, **row
                })
    
            one_complete = float(np.mean([
                row["one_crossing_and_complete"] for row in cycles
            ]))
            zero = float(np.mean([row["total_crossings"] == 0 for row in cycles]))
            multiple = float(np.mean([row["total_crossings"] > 1 for row in cycles]))
            alternation = alternation_fidelity(cycles)
            pass_initial = bool(
                one_complete >= SETTINGS["success_target"]
                and np.isfinite(alternation)
                and alternation >= SETTINGS["alternation_target"]
                and conservation <= SETTINGS["DNA_conservation_tolerance"]
            )
            summary_rows.append({
                "condition": label,
                "C31_RBS_scale": scale,
                "Int_tag_rate_h_inv": tag,
                "Int_total_half_life_min": 60.0*math.log(2.0)/(P["k_dil"]+tag),
                "nominal_peak_load_fraction": load,
                "initial_DNA_state": initial_name,
                "A_period_h": period_h,
                "A_period_generations": 60.0*period_h/p.Td_min,
                "A_TetR_free_peak_copies": tetR_peak_copies,
                "realized_resource_min": float(np.min(rho[reference_keep])),
                "realized_resource_mean": float(np.mean(rho[reference_keep])),
                "flux_peak_uM_h_median": float(np.median([
                    row["flux_peak_uM_h"] for row in a_cycles
                ])),
                "flux_FWHM_h_median": float(np.median([
                    row["flux_FWHM_h"] for row in a_cycles
                ])),
                "flux_dose_uM_cycle_median": float(np.median([
                    row["flux_dose_uM_per_cycle"] for row in a_cycles
                ])),
                "audited_cycles": len(cycles),
                "one_crossing_complete_fraction": one_complete,
                "zero_crossing_fraction": zero,
                "multiple_crossing_fraction": multiple,
                "alternation_fidelity": alternation,
                "DNA_conservation_residual": conservation,
                "pass_initial_state": pass_initial,
            })
            saved_examples[(label, initial_name)] = (
                b_time, np.maximum(0.0, spline(b_time)), lr, pb,
                zc.int_total(states), zc.rdf_total(states),
            )
    
    # A design passes only if all requested initial DNA states pass at that load.
    for scale, tag, load in conditions:
        matching = [
            row for row in summary_rows
            if row["C31_RBS_scale"] == scale
            and row["Int_tag_rate_h_inv"] == tag
            and row["nominal_peak_load_fraction"] == load
        ]
        both_initials_B_pass = bool(
            len(matching) == len(SETTINGS["DNA_initial_states"])
            and all(row["pass_initial_state"] for row in matching)
        )
        reference_matching = [
            row for row in summary_rows
            if row["C31_RBS_scale"] == scale
            and row["Int_tag_rate_h_inv"] == tag
            and row["nominal_peak_load_fraction"] == 0.0
        ]
        if not reference_matching:
            raise RuntimeError("Every candidate needs a zero-load reference")
        a_row = matching[0]
        reference_row = reference_matching[0]
        period_change_percent = 100.0*(
            a_row["A_period_h"]/reference_row["A_period_h"] - 1.0
        )
        tetR_peak_change_percent = 100.0*(
            a_row["A_TetR_free_peak_copies"]
            / reference_row["A_TetR_free_peak_copies"] - 1.0
        )
        A_clock_pass = bool(
            abs(period_change_percent) <= SETTINGS["A_period_tolerance_percent"]
            and abs(tetR_peak_change_percent)
            <= SETTINGS["A_TetR_peak_tolerance_percent"]
        )
        design_rows.append({
            "C31_RBS_scale": scale,
            "Int_tag_rate_h_inv": tag,
            "nominal_peak_load_fraction": load,
            "A_period_change_percent": period_change_percent,
            "A_TetR_peak_change_percent": tetR_peak_change_percent,
            "A_clock_pass": A_clock_pass,
            "both_initial_states_B_pass": both_initials_B_pass,
            "complete_system_pass": bool(A_clock_pass and both_initials_B_pass),
            "worst_one_crossing_complete_fraction": min(
                row["one_crossing_complete_fraction"] for row in matching
            ),
            "worst_alternation_fidelity": min(
                row["alternation_fidelity"] for row in matching
                if np.isfinite(row["alternation_fidelity"])
            ) if any(np.isfinite(row["alternation_fidelity"]) for row in matching) else float("nan"),
            "max_DNA_conservation_residual": max(
                row["DNA_conservation_residual"] for row in matching
            ),
        })
    
    write_rows(output/"v53h_condition_summary.csv", summary_rows)
    write_rows(output/"v53h_A_cycle_metrics.csv", a_cycle_rows)
    write_rows(output/"v53h_B_cycle_metrics.csv", b_cycle_rows)
    write_rows(output/"v53h_design_robustness.csv", design_rows)
    
    with (output/"v53h_run_metadata.json").open("w", encoding="utf-8") as handle:
        json.dump({
            "version": VERSION,
            "settings": SETTINGS,
            "A_reference_parameters": asdict(base),
            "reference_C31_demand_peak_copies_per_min": reference_demand_peak,
            "resource_equation": "rho=1/(1+[L/(1-L)]*(beta31*m31/Jref))",
            "interface_equation": "int_production=rho*beta31*m31 converted to uM/h",
        }, handle, indent=2, ensure_ascii=False)
    
    if SETTINGS["make_plots"]:
        # Main robustness plot: worst result across PB/LR starts at each load.
        fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
        for scale, tag in candidates:
            subset = [
                row for row in design_rows
                if row["C31_RBS_scale"] == scale
                and row["Int_tag_rate_h_inv"] == tag
            ]
            subset.sort(key=lambda row: row["nominal_peak_load_fraction"])
            x = 100*np.array([row["nominal_peak_load_fraction"] for row in subset])
            axes[0].plot(
                x, [row["worst_one_crossing_complete_fraction"] for row in subset],
                marker="o", label=f"RBS={scale:g}, tag={tag:g}/h",
            )
            axes[1].plot(
                x, [row["complete_system_pass"] for row in subset],
                marker="s", label=f"RBS={scale:g}, tag={tag:g}/h",
            )
        axes[0].axhline(SETTINGS["success_target"], color="black", ls="--")
        axes[0].set(xlabel="nominal peak resource load (%)",
                    ylabel="worst one-crossing + complete fraction", ylim=(-0.03, 1.05))
        axes[1].set(xlabel="nominal peak resource load (%)",
                    ylabel="complete A+B system passes", ylim=(-0.05, 1.05),
                    yticks=(0, 1))
        axes[0].legend(fontsize=8); axes[1].legend(fontsize=8)
        fig.suptitle("v53h long A-to-B validation under self-consistent resource load")
        fig.tight_layout()
        fig.savefig(output/"v53h_robustness_summary.png", dpi=180)
        plt.show()
    
        # Plot the main candidate at the lowest and highest load, PB initial state.
        main_scale, main_tag = map(float, SETTINGS["main_candidate"])
        available = sorted(set(
            row["nominal_peak_load_fraction"] for row in summary_rows
            if row["C31_RBS_scale"] == main_scale
            and row["Int_tag_rate_h_inv"] == main_tag
        ))
        if available:
            chosen_loads = (available[0], available[-1])
            fig, axes = plt.subplots(len(chosen_loads), 3, figsize=(15, 7), sharex=True)
            axes = np.atleast_2d(axes)
            for row_index, load in enumerate(chosen_loads):
                label = (
                    f"RBS_{main_scale:g}__tag_{main_tag:g}__load_{100*load:g}pct"
                    .replace(".", "p")
                )
                b_time, b_input, lr, pb, int_total, rdf_total = saved_examples[(label, "PB")]
                view = b_time >= SETTINGS["A_transient_h"]
                axes[row_index, 0].plot(b_time[view], b_input[view], color="magenta")
                axes[row_index, 1].plot(b_time[view], lr[view], "r", label="LR total")
                axes[row_index, 1].plot(b_time[view], pb[view], "b", label="PB total")
                axes[row_index, 2].plot(b_time[view], int_total[view], "g", label="Int total")
                axes[row_index, 2].plot(b_time[view], rdf_total[view], "k", label="RDF total")
                axes[row_index, 0].set_ylabel(f"load={100*load:g}%\ninput (uM/h)")
                axes[row_index, 1].set_ylabel("DNA fraction")
                axes[row_index, 2].set_ylabel("uM")
            for axis in axes[-1]:
                axis.set_xlabel("time (h)")
            axes[0, 1].legend(); axes[0, 2].legend()
            fig.suptitle("Main fixed-PLtetO1 candidate: long trajectories from PB")
            fig.tight_layout()
            fig.savefig(output/"v53h_main_candidate_trajectories.png", dpi=180)
            plt.show()
    
    print("\nLoaded diagnostic code:", VERSION)
    print("Saved v53h outputs to:", output)
    print(
        "Resource model: rho=1/(1+[L/(1-L)]*(beta31*m31/Jref)); "
        "A is re-solved for every RBS/load pair."
    )
    for row in design_rows:
        print(
            f"RBS={row['C31_RBS_scale']:.2f} | tag={row['Int_tag_rate_h_inv']:g}/h | "
            f"load={100*row['nominal_peak_load_fraction']:g}% | "
            f"worst one+complete={row['worst_one_crossing_complete_fraction']:.2f} | "
            f"A clock pass={row['A_clock_pass']} | "
            f"B both initials pass={row['both_initial_states_B_pass']} | "
            f"system pass={row['complete_system_pass']}"
        )
    print("\nInterpretation boundaries:")
    print("- PLtetO1 K, n, leak and topology are identical to cI and remain frozen.")
    print("- Nominal load fractions are uncertainty scenarios, not measured burden.")
    print("- RBS=0.60/tag=12 h^-1 can be frozen only if it passes both PB/LR starts across the chosen load range.")
    print("- This validates deterministic interface robustness; it is not a stochastic fidelity estimate.")
    
    if display is not None:
        try:
            import pandas as pd
            display(pd.DataFrame(design_rows))
        except Exception:
            pass


if __name__ == "__main__":
    main()
