# Standalone program: mechanistic repressilator-C31 ODE with explicit mRNA.
# Source version v52a1: symmetric oscillator mRNA-protein ODE with QSS mapping.
# It embeds the frozen v52a0.1 reduced model as the comparison baseline.

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

# ============================== USER SETTINGS ==============================
SETTINGS = {
    "duration_min": 6000.0,
    "sample_step_min": 1.0,
    "transient_min": 1500.0,
    "Td_min": 50.0,
    "make_plots": True,
    "minimum_relative_amplitude": 0.01,
    "minimum_peaks": 4,
    "fast_limit_max_period_error_percent": 10.0,
    "fast_limit_min_waveform_correlation": 0.95,
}

# Only the four requested full-model conditions. The reduced baseline is run
# once for each Nt and is not a parameter scan.
SPONGE_PLASMID_COPIES = (0.0, 40.0)
MRNA_SCENARIOS = {
    "fast_mRNA_0p2min": 0.2,   # 10x faster total loss than the 2-min case
    "mRNA_2min": 2.0,
}
# ============================================================================

@dataclass(frozen=True)
class Parameters:
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

    # Symmetric oscillator translation and protein loss.
    oscillator_translation_per_mrna_per_min: float = 0.5
    oscillator_protein_intrinsic_loss_per_min: float = 0.0

    # C31 replaces the reporter on the same plasmid. Its two tetO sites are
    # already included in functional_tetO_per_plasmid=4.
    c31_promoter_leak: float = 0.005
    c31_tx_per_plasmid_per_min: float = 5.0
    c31_translation_per_mrna_per_min: float = 0.5
    c31_mrna_total_half_life_min: float = 2.0
    c31_protein_intrinsic_loss_per_min: float = 0.0

    @property
    def mu(self):
        return math.log(2.0) / self.Td_min

    @property
    def effective_functional_sites(self):
        return (self.functional_accessible_fraction * self.functional_tetO_per_plasmid
                * self.oscillator_plasmid_copies)

    @property
    def effective_sponge_sites(self):
        return (self.sponge_accessible_fraction * self.sponge_tetO_per_plasmid
                * self.sponge_plasmid_copies)

    @property
    def c31_mrna_intrinsic_loss(self):
        return max(0.0, math.log(2.0) / self.c31_mrna_total_half_life_min - self.mu)

def repression(repressor, p, leak=0.0):
    r = max(float(repressor), 0.0)
    regulated = p.K_copies**p.hill_n / (p.K_copies**p.hill_n + r**p.hill_n)
    return leak + (1.0 - leak) * regulated

def occupancy(free, p):
    f = max(float(free), 0.0)
    if f == 0.0:
        return 0.0
    return f**p.hill_n / (p.K_copies**p.hill_n + f**p.hill_n)

def tetR_balance(free, total, p):
    sites = p.effective_functional_sites + p.effective_sponge_sites
    return free + sites * occupancy(free, p) - max(float(total), 0.0)

def free_tetR(total, p):
    total = max(float(total), 0.0)
    if total == 0.0:
        return 0.0
    free = brentq(lambda value: tetR_balance(value, total, p), 0.0, total,
                  xtol=1e-12, rtol=1e-12, maxiter=100)
    residual = abs(tetR_balance(free, total, p)) / max(total, 1.0)
    if residual >= 1e-6 or not 0.0 <= free <= total:
        raise RuntimeError(f"TetR conservation failed: residual={residual:g}")
    return free

def promoter_activities(tet_total, ci, laci, p):
    tet_free = free_tetR(tet_total, p)
    return (
        repression(laci, p, p.oscillator_promoter_leak),  # TetR production
        repression(tet_free, p, p.oscillator_promoter_leak),  # CI production
        repression(ci, p, p.oscillator_promoter_leak),  # LacI production
        repression(tet_free, p, p.c31_promoter_leak),  # C31 production
        tet_free,
    )

def reduced_rhs(_, y, p):
    tet, ci, laci, c31_mrna, c31 = y
    h_tet, h_ci, h_laci, h_c31, _ = promoter_activities(tet, ci, laci, p)
    mu = p.mu
    return (
        mu * (p.lambda_protein_scale * h_tet - tet),
        mu * (p.lambda_protein_scale * h_ci - ci),
        mu * (p.lambda_protein_scale * h_laci - laci),
        p.oscillator_plasmid_copies * p.c31_tx_per_plasmid_per_min * h_c31
        - (p.c31_mrna_intrinsic_loss + mu) * c31_mrna,
        p.c31_translation_per_mrna_per_min * c31_mrna
        - (p.c31_protein_intrinsic_loss_per_min + mu) * c31,
    )

def oscillator_mrna_rates(total_mrna_half_life_min, p):
    gamma_m_total = math.log(2.0) / total_mrna_half_life_min
    delta_m_intrinsic = max(0.0, gamma_m_total - p.mu)
    beta = p.oscillator_translation_per_mrna_per_min
    gamma_p_total = p.oscillator_protein_intrinsic_loss_per_min + p.mu
    # QSS mapping:
    # beta*No*alpha/(delta_m+mu) = (delta_p+mu)*lambda.
    alpha = gamma_p_total * p.lambda_protein_scale * gamma_m_total / (
        beta * p.oscillator_plasmid_copies
    )
    return alpha, delta_m_intrinsic, gamma_m_total

def full_rhs(_, y, p, total_mrna_half_life_min):
    m_tet, tet, m_ci, ci, m_laci, laci, c31_mrna, c31 = y
    h_tet, h_ci, h_laci, h_c31, _ = promoter_activities(tet, ci, laci, p)
    alpha, delta_m, _ = oscillator_mrna_rates(total_mrna_half_life_min, p)
    beta = p.oscillator_translation_per_mrna_per_min
    protein_loss = p.oscillator_protein_intrinsic_loss_per_min + p.mu
    mrna_loss = delta_m + p.mu
    return (
        p.oscillator_plasmid_copies * alpha * h_tet - mrna_loss * m_tet,
        beta * m_tet - protein_loss * tet,
        p.oscillator_plasmid_copies * alpha * h_ci - mrna_loss * m_ci,
        beta * m_ci - protein_loss * ci,
        p.oscillator_plasmid_copies * alpha * h_laci - mrna_loss * m_laci,
        beta * m_laci - protein_loss * laci,
        p.oscillator_plasmid_copies * p.c31_tx_per_plasmid_per_min * h_c31
        - (p.c31_mrna_intrinsic_loss + p.mu) * c31_mrna,
        p.c31_translation_per_mrna_per_min * c31_mrna
        - (p.c31_protein_intrinsic_loss_per_min + p.mu) * c31,
    )

def full_initial_state(protein_initial, p, total_mrna_half_life_min):
    tet, ci, laci = protein_initial
    h_tet, h_ci, h_laci, _, _ = promoter_activities(tet, ci, laci, p)
    alpha, _, gamma_m_total = oscillator_mrna_rates(total_mrna_half_life_min, p)
    factor = p.oscillator_plasmid_copies * alpha / gamma_m_total
    # mRNA begins at its QSS for the current repressors. Consequently the
    # initial protein derivative matches the reduced model exactly.
    return (factor * h_tet, tet, factor * h_ci, ci, factor * h_laci, laci, 0.0, 0.0)

def oscillation_metrics(t, signal, p):
    span = float(np.max(signal) - np.min(signal))
    scale = max(float(np.max(np.abs(signal))), 1.0)
    relative_amplitude = span / scale
    if relative_amplitude < SETTINGS["minimum_relative_amplitude"]:
        return {"status": "steady_state", "period_min": float("nan"),
                "period_generations": float("nan"), "period_cv": float("nan"),
                "relative_amplitude": relative_amplitude, "peaks": np.array([], dtype=int)}
    minimum_distance = max(1, int(round(2.0 * p.Td_min / SETTINGS["sample_step_min"])))
    peaks, _ = find_peaks(signal, prominence=max(0.05 * span, 1e-9), distance=minimum_distance)
    if len(peaks) < SETTINGS["minimum_peaks"]:
        return {"status": "unresolved_oscillation", "period_min": float("nan"),
                "period_generations": float("nan"), "period_cv": float("nan"),
                "relative_amplitude": relative_amplitude, "peaks": peaks}
    periods = np.diff(t[peaks])
    period = float(np.mean(periods))
    return {"status": "resolved_oscillation", "period_min": period,
            "period_generations": period / p.Td_min,
            "period_cv": float(np.std(periods, ddof=1) / period) if len(periods) > 1 else 0.0,
            "relative_amplitude": relative_amplitude, "peaks": peaks}

def last_cycle_phase(signal, peaks, points=500):
    if len(peaks) < 2:
        return None
    start, stop = peaks[-2], peaks[-1]
    original_phase = np.linspace(0.0, 1.0, stop - start + 1)
    target_phase = np.linspace(0.0, 1.0, points)
    return np.interp(target_phase, original_phase, signal[start:stop + 1])

def waveform_correlation(reference, target):
    if reference is None or target is None or np.std(reference) == 0 or np.std(target) == 0:
        return float("nan")
    return float(np.corrcoef(reference, target)[0, 1])

def cycle_dose(t, c31, peaks):
    if len(peaks) < 2:
        return float("nan")
    segment = slice(peaks[-2], peaks[-1] + 1)
    return float(trapezoid(c31[segment], x=t[segment]))


def main():
    base = Parameters(Td_min=SETTINGS["Td_min"])
    time = np.arange(0.0, SETTINGS["duration_min"] + 0.5 * SETTINGS["sample_step_min"],
                     SETTINGS["sample_step_min"])
    keep = time >= SETTINGS["transient_min"]
    t_ss = time[keep]
    protein_initial = (100.0, 20.0, 300.0)
    output = OUTPUT_ROOT / f"Mechanistic_Repressilator_C31_ODE_{datetime.now():%Y%m%d_%H%M%S}"
    output.mkdir(parents=True, exist_ok=False)
    rows, saved, reduced_by_Nt = [], {}, {}
    
    # Frozen reduced reference for Nt=0 and Nt=40.
    for Nt in SPONGE_PLASMID_COPIES:
        p = replace(base, sponge_plasmid_copies=Nt)
        sol = solve_ivp(lambda t, y: reduced_rhs(t, y, p), (time[0], time[-1]),
                        (*protein_initial, 0.0, 0.0), t_eval=time, method="LSODA",
                        rtol=1e-7, atol=1e-9)
        if not sol.success:
            raise RuntimeError(sol.message)
        tet_free = np.array([free_tetR(value, p) for value in sol.y[0]])
        activity = np.array([repression(value, p, p.c31_promoter_leak) for value in tet_free])
        osc = oscillation_metrics(t_ss, tet_free[keep], p)
        scenario = f"reduced__Nt_{Nt:g}"
        reduced_by_Nt[Nt] = {"p": p, "sol": sol, "tet_free": tet_free,
                             "activity": activity, "osc": osc}
        saved[scenario] = (p, sol.y, tet_free, activity, "reduced")
    
    # Four full-model conditions.
    for Nt in SPONGE_PLASMID_COPIES:
        p = replace(base, sponge_plasmid_copies=Nt)
        reference = reduced_by_Nt[Nt]
        ref_osc = reference["osc"]
        ref_tet_cycle = last_cycle_phase(reference["tet_free"][keep], ref_osc["peaks"])
        ref_ci_cycle = last_cycle_phase(reference["sol"].y[1, keep], ref_osc["peaks"])
        ref_laci_cycle = last_cycle_phase(reference["sol"].y[2, keep], ref_osc["peaks"])
        for mrna_name, mrna_half_life in MRNA_SCENARIOS.items():
            initial = full_initial_state(protein_initial, p, mrna_half_life)
            sol = solve_ivp(lambda t, y: full_rhs(t, y, p, mrna_half_life),
                            (time[0], time[-1]), initial, t_eval=time, method="LSODA",
                            rtol=1e-7, atol=1e-9)
            if not sol.success:
                raise RuntimeError(f"Nt={Nt}, {mrna_name}: {sol.message}")
            tet_free = np.array([free_tetR(value, p) for value in sol.y[1]])
            activity = np.array([repression(value, p, p.c31_promoter_leak) for value in tet_free])
            osc = oscillation_metrics(t_ss, tet_free[keep], p)
            tet_cycle = last_cycle_phase(tet_free[keep], osc["peaks"])
            ci_cycle = last_cycle_phase(sol.y[3, keep], osc["peaks"])
            laci_cycle = last_cycle_phase(sol.y[5, keep], osc["peaks"])
            period_error = (100.0 * (osc["period_min"] / ref_osc["period_min"] - 1.0)
                            if osc["status"] == ref_osc["status"] == "resolved_oscillation"
                            else float("nan"))
            correlations = [waveform_correlation(ref_tet_cycle, tet_cycle),
                            waveform_correlation(ref_ci_cycle, ci_cycle),
                            waveform_correlation(ref_laci_cycle, laci_cycle)]
            c31_ss = sol.y[7, keep]
            promoter_ss = activity[keep]
            alpha, delta_m, gamma_m_total = oscillator_mrna_rates(mrna_half_life, p)
            scenario = f"full__{mrna_name}__Nt_{Nt:g}"
            fast_case = mrna_name.startswith("fast_")
            qss_pass = (not fast_case or
                        (np.isfinite(period_error)
                         and abs(period_error) <= SETTINGS["fast_limit_max_period_error_percent"]
                         and np.nanmin(correlations) >= SETTINGS["fast_limit_min_waveform_correlation"]))
            row = {
                "scenario": scenario, "sponge_plasmid_copies": Nt,
                "oscillator_mrna_total_half_life_min": mrna_half_life,
                "oscillator_mrna_intrinsic_loss_per_min": delta_m,
                "oscillator_tx_per_plasmid_per_min_QSS_mapped": alpha,
                "oscillator_translation_per_mrna_per_min": p.oscillator_translation_per_mrna_per_min,
                "status": osc["status"], "period_min": osc["period_min"],
                "period_generations": osc["period_generations"],
                "deterministic_period_cv": osc["period_cv"],
                "period_error_vs_reduced_percent": period_error,
                "TetR_cycle_correlation_vs_reduced": correlations[0],
                "CI_cycle_correlation_vs_reduced": correlations[1],
                "LacI_cycle_correlation_vs_reduced": correlations[2],
                "fast_QSS_acceptance_pass": qss_pass,
                "TetR_free_peak": float(np.max(tet_free[keep])),
                "TetR_free_trough": float(np.min(tet_free[keep])),
                "C31_peak": float(np.max(c31_ss)), "C31_trough": float(np.min(c31_ss)),
                "C31_dose_copies_min_per_cycle": cycle_dose(t_ss, c31_ss, osc["peaks"]),
                "PLtetO1_duty_activity_gt_0p5": float(np.mean(promoter_ss > 0.5)),
                "effective_functional_sites": p.effective_functional_sites,
                "effective_sponge_sites": p.effective_sponge_sites,
                "QSS_mapping_lhs": (p.oscillator_translation_per_mrna_per_min
                                    * p.oscillator_plasmid_copies * alpha / gamma_m_total),
                "QSS_mapping_rhs": ((p.oscillator_protein_intrinsic_loss_per_min + p.mu)
                                    * p.lambda_protein_scale),
            }
            rows.append(row)
            saved[scenario] = (p, sol.y, tet_free, activity, "full")
    
            folder = output / scenario
            folder.mkdir()
            with (folder / "trajectory.csv").open("w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle)
                writer.writerow(["time_min", "time_generations", "TetR_mRNA", "TetR_total",
                                 "CI_mRNA", "CI", "LacI_mRNA", "LacI", "C31_mRNA",
                                 "C31_protein", "TetR_free", "PLtetO1_activity"])
                for index, t in enumerate(time):
                    writer.writerow([t, t / p.Td_min, *sol.y[:, index], tet_free[index], activity[index]])
            with (folder / "parameters.json").open("w", encoding="utf-8") as handle:
                json.dump({"parameters": asdict(p), "mRNA_scenario": mrna_name,
                           "mRNA_total_half_life_min": mrna_half_life,
                           "QSS_mapped_alpha": alpha, "initial_state": initial},
                          handle, ensure_ascii=False, indent=2)
    
    with (output / "v52a1_metrics.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    with (output / "metadata.json").open("w", encoding="utf-8") as handle:
        json.dump({"settings": SETTINGS, "base_parameters": asdict(base),
                   "mRNA_scenarios": MRNA_SCENARIOS,
                   "QSS_mapping": "beta*No*alpha/(delta_m+mu)=(delta_p+mu)*lambda",
                   "purpose": "Test whether the symmetric 6-state oscillator is a credible mechanistic extension of the frozen reduced ODE."},
                  handle, ensure_ascii=False, indent=2)
    
    if SETTINGS["make_plots"]:
        fig, axes = plt.subplots(2, 2, figsize=(13, 9), sharex=True, constrained_layout=True)
        for column, Nt in enumerate(SPONGE_PLASMID_COPIES):
            reference = reduced_by_Nt[Nt]
            axes[0, column].plot(time / base.Td_min, reference["tet_free"], color="black",
                                 linewidth=1.4, label="reduced")
            axes[1, column].plot(time / base.Td_min, reference["sol"].y[4], color="black",
                                 linewidth=1.4, label="reduced")
            for mrna_name, color in zip(MRNA_SCENARIOS, ("tab:blue", "tab:orange")):
                scenario = f"full__{mrna_name}__Nt_{Nt:g}"
                p, y, tet_free, activity, _ = saved[scenario]
                axes[0, column].plot(time / p.Td_min, tet_free, color=color, alpha=0.85, label=mrna_name)
                axes[1, column].plot(time / p.Td_min, y[7], color=color, alpha=0.85, label=mrna_name)
            axes[0, column].set_title(f"Nt={Nt:g}, eta=1")
            axes[0, column].set_ylabel("free TetR (copies/cell)")
            axes[1, column].set_ylabel("C31 protein (copies/cell)")
            axes[1, column].set_xlabel("time (generations)")
            axes[0, column].legend(); axes[1, column].legend()
        figure_path = output / "v52a1_reduced_vs_full.png"
        fig.savefig(figure_path, dpi=180)
        if notebook_display is not None:
            notebook_display(fig)
        plt.close(fig)
    
    print(f"Saved v52a1 outputs to: {output.resolve()}")
    for Nt in SPONGE_PLASMID_COPIES:
        ref = reduced_by_Nt[Nt]["osc"]
        print(f"reduced Nt={Nt:g}: {ref['status']}, P={ref['period_min']:.1f} min/{ref['period_generations']:.3f} gen")
    for row in rows:
        print(f"{row['scenario']:>30s} | {row['status']:>22s} | P={row['period_min']:.1f} min/{row['period_generations']:.3f} gen "
              f"| error={row['period_error_vs_reduced_percent']:+.2f}% "
              f"| corr(T/C/L)={row['TetR_cycle_correlation_vs_reduced']:.3f}/"
              f"{row['CI_cycle_correlation_vs_reduced']:.3f}/{row['LacI_cycle_correlation_vs_reduced']:.3f} "
              f"| QSS_pass={row['fast_QSS_acceptance_pass']} | C31={row['C31_peak']:.1f}/{row['C31_trough']:.1f}")


if __name__ == "__main__":
    main()
