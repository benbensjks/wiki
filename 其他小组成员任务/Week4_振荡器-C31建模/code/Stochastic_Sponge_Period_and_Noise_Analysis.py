"""Standalone paper-aligned stochastic period and noise analysis.

Source version: v54s0 period-estimator audit.

Purpose
-------
Keep the v47 / Potvin-Trottier SI stochastic reaction system unchanged, but
compare four period estimators on exactly the same trajectories:

1. legacy CI-total peaks;
2. PLtetO1 activation-onset intervals;
3. three-repressor cyclic phase crossings;
4. reporter-proxy ACF and PSD periods.

This version does NOT tune parameters to force sponge-mediated period
lengthening. It asks whether the apparent sign depends on how a noisy,
low-amplitude no-sponge trajectory is measured.

No simulation was run while preparing this file.
"""

from dataclasses import dataclass, replace
from datetime import datetime
from functools import lru_cache
from pathlib import Path
import csv
import json
import math
import time

import numpy as np
from scipy.signal import find_peaks, savgol_filter

# Portable project paths: independent of the terminal's current directory.
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
DATA_DIR = PROJECT_ROOT / "data"
OUTPUT_ROOT = PROJECT_ROOT / "outputs"

try:
    import matplotlib.pyplot as plt
except ModuleNotFoundError:
    plt = None


TETR, CI, LACI = 0, 1, 2


@dataclass(frozen=True)
class PaperParameters:
    # Potvin-Trottier SI 4.3.1 / Extended Data Fig. 4f.
    lambda_events_per_promoter_per_lifetime: float = 60.0
    mean_original_plasmids: int = 10
    mean_sponge_plasmids: int = 40
    sites_per_promoter: int = 2
    effective_sponge_fraction: float = 1.0
    hill_n: float = 1.5
    burst_mean: float = 10.0
    K_copies: tuple = (5.0, 10.0, 10.0)  # TetR, CI, LacI
    # Only maps dimensionless protein-lifetime time to minutes. It does not
    # change SSA dynamics or periods reported in generations. We use the
    # current A/B working value for MC4100 in LB at 37 C.
    doubling_time_min: float = 50.0

    @property
    def protein_lifetime_min(self):
        return self.doubling_time_min / math.log(2.0)

    @property
    def generations_per_lifetime(self):
        return 1.0 / math.log(2.0)


@dataclass(frozen=True)
class Settings:
    trajectories_per_condition: int
    total_time_lifetimes: float
    transient_lifetimes: float
    sample_dt_lifetimes: float
    random_seed: int = 20260815
    bootstrap_replicates: int = 3000


# Start with "audit". Use "confirmation" only after inspecting audit plots.
RUN_MODE = "audit"  # "smoke", "audit", or "confirmation"
MODE_SETTINGS = {
    "smoke": Settings(3, 90.0, 20.0, 0.10, bootstrap_replicates=500),
    "audit": Settings(12, 160.0, 30.0, 0.10, bootstrap_replicates=2000),
    "confirmation": Settings(50, 180.0, 30.0, 0.10, bootstrap_replicates=5000),
}
SETTINGS = MODE_SETTINGS[RUN_MODE]
P = PaperParameters()

# Keep one mode for the first run. Adding both doubles compute time.
# paper_si: preserve the historical v43/v47 2*N_o capacity convention.
# our_pLPT119: two PLtetO1 promoters on the oscillator plasmid => 4*N_o TetO.
CONSTRUCT_MODES = ("paper_si",)
# After the paper_si audit, optionally use:
# CONSTRUCT_MODES = ("paper_si", "our_pLPT119")

OUTPUT_ROOT = PROJECT_ROOT / "outputs"

# Pre-registered analysis ranges in protein lifetimes.
MIN_PERIOD_LIFETIMES = 4.0
MAX_PERIOD_LIFETIMES = 16.0
ACTIVITY_ON_THRESHOLD = 0.50
ACTIVITY_OFF_THRESHOLD = 0.25


def geometric_burst(rng, mean_burst):
    probability = 1.0 / (1.0 + mean_burst)
    return int(rng.geometric(probability) - 1)


@lru_cache(maxsize=700_000)
def solve_free_repressor_cached(total_integer, site_capacity_rounded, K, n):
    total = float(max(total_integer, 0))
    capacity = max(float(site_capacity_rounded), 0.0)
    if total == 0.0 or capacity == 0.0:
        return total

    low, high = 0.0, total
    free = max(total - 0.5 * capacity, 0.5 * total)
    K_n = K**n
    for _ in range(28):
        free_n = free**n
        denominator = free_n + K_n
        residual = free + capacity * free_n / denominator - total
        if abs(residual) <= 1e-11 * max(total, 1.0):
            return free
        if residual > 0.0:
            high = free
        else:
            low = free
        derivative = 1.0 + (
            capacity
            * n
            * K_n
            * max(free, 1e-15) ** (n - 1.0)
            / denominator**2
        )
        candidate = free - residual / derivative
        free = candidate if low < candidate < high else 0.5 * (low + high)
    return free


def hill_unrepressed(free, K, n):
    return K**n / (K**n + max(float(free), 0.0) ** n)


def functional_capacities(N_o, construct_mode, p=P):
    base = p.sites_per_promoter * N_o
    if construct_mode == "paper_si":
        return np.array([base, base, base], dtype=float)
    if construct_mode == "our_pLPT119":
        # cI and phiC31 each have a PLtetO1 with two TetO sites.
        return np.array([2.0 * base, base, base], dtype=float)
    raise ValueError(f"Unknown construct mode: {construct_mode}")


def free_repressors(proteins_total, N_o, N_t, with_sponge, construct_mode, p=P):
    capacities = functional_capacities(N_o, construct_mode, p)
    if with_sponge:
        capacities[TETR] += (
            p.effective_sponge_fraction * p.sites_per_promoter * N_t
        )
    return np.array(
        [
            solve_free_repressor_cached(
                int(proteins_total[i]),
                round(float(capacities[i]), 6),
                p.K_copies[i],
                p.hill_n,
            )
            for i in range(3)
        ],
        dtype=float,
    )


def promoter_activities(free, p=P):
    # LacI -> TetR, TetR -> CI, CI -> LacI.
    return np.array(
        [
            hill_unrepressed(free[LACI], p.K_copies[LACI], p.hill_n),
            hill_unrepressed(free[TETR], p.K_copies[TETR], p.hill_n),
            hill_unrepressed(free[CI], p.K_copies[CI], p.hill_n),
        ],
        dtype=float,
    )


def gillespie_trajectory(with_sponge, seed, construct_mode, p=P, settings=SETTINGS):
    rng = np.random.default_rng(seed)
    sample_times = np.arange(
        0.0,
        settings.total_time_lifetimes + 0.5 * settings.sample_dt_lifetimes,
        settings.sample_dt_lifetimes,
    )
    proteins_out = np.zeros((sample_times.size, 3), dtype=np.int32)
    free_out = np.zeros((sample_times.size, 3), dtype=np.float32)
    activities_out = np.zeros((sample_times.size, 3), dtype=np.float32)
    plasmids_out = np.zeros((sample_times.size, 2), dtype=np.int16)

    proteins = np.array([1, 1, 1], dtype=np.int64)
    N_o = int(p.mean_original_plasmids)
    N_t = int(p.mean_sponge_plasmids if with_sponge else 0)
    tau = 0.0
    sample_index = 0
    event_count = 0

    while sample_index < sample_times.size:
        free = free_repressors(
            proteins, N_o, N_t, with_sponge, construct_mode, p
        )
        activities = promoter_activities(free, p)
        production_rates = p.lambda_events_per_promoter_per_lifetime * N_o * activities
        loss_rates = proteins.astype(float)
        plasmid_rates = np.array(
            [
                float(p.mean_original_plasmids),
                float(N_o),
                float(p.mean_sponge_plasmids if with_sponge else 0.0),
                float(N_t if with_sponge else 0.0),
            ],
            dtype=float,
        )
        rates = np.concatenate([production_rates, loss_rates, plasmid_rates])
        total_rate = float(rates.sum())
        if total_rate <= 0.0:
            raise RuntimeError("All Gillespie propensities became zero")

        next_tau = tau - math.log(max(rng.random(), 1e-300)) / total_rate
        while sample_index < sample_times.size and sample_times[sample_index] < next_tau:
            proteins_out[sample_index] = proteins
            free_out[sample_index] = free
            activities_out[sample_index] = activities
            plasmids_out[sample_index] = (N_o, N_t)
            sample_index += 1
        if sample_index >= sample_times.size:
            break

        tau = next_tau
        threshold = rng.random() * total_rate
        event = int(np.searchsorted(np.cumsum(rates), threshold, side="right"))
        if event < 3:
            proteins[event] += geometric_burst(rng, p.burst_mean)
        elif event < 6:
            i = event - 3
            if proteins[i] > 0:
                proteins[i] -= 1
        elif event == 6:
            N_o += 1
        elif event == 7:
            if N_o > 0:
                N_o -= 1
        elif event == 8:
            N_t += 1
        elif event == 9:
            if N_t > 0:
                N_t -= 1
        event_count += 1

    return {
        "tau": sample_times,
        "proteins_total": proteins_out,
        "proteins_free": free_out,
        "activities": activities_out,
        "plasmids": plasmids_out,
        "event_count": event_count,
        "seed": seed,
        "with_sponge": with_sponge,
        "construct_mode": construct_mode,
    }


def smooth_signal(signal, sample_dt, window_lifetimes=0.8):
    window = max(5, int(round(window_lifetimes / sample_dt)))
    if window % 2 == 0:
        window += 1
    max_odd = len(signal) if len(signal) % 2 == 1 else len(signal) - 1
    window = min(window, max_odd)
    if window < 5:
        return np.asarray(signal, dtype=float)
    return savgol_filter(np.asarray(signal, dtype=float), window, polyorder=2)


def valid_intervals(event_times):
    periods = np.diff(np.asarray(event_times, dtype=float))
    return periods[
        (periods >= MIN_PERIOD_LIFETIMES)
        & (periods <= MAX_PERIOD_LIFETIMES)
    ]


def legacy_ci_peak_periods(tau, ci_total, settings=SETTINGS):
    signal = smooth_signal(ci_total, settings.sample_dt_lifetimes)
    peaks, _ = find_peaks(
        signal,
        distance=max(1, int(round(4.5 / settings.sample_dt_lifetimes))),
        prominence=max(0.12 * float(np.ptp(signal)), 20.0),
    )
    return valid_intervals(tau[peaks]), tau[peaks], signal


def activity_onset_periods(tau, activity):
    """Hysteretic PLtetO1 activation onsets using fixed thresholds."""
    active = False
    onsets = []
    for i, value in enumerate(activity):
        if not active and value >= ACTIVITY_ON_THRESHOLD:
            onsets.append(float(tau[i]))
            active = True
        elif active and value <= ACTIVITY_OFF_THRESHOLD:
            active = False
    return valid_intervals(onsets), np.asarray(onsets, dtype=float)


def cyclic_phase_periods(tau, proteins):
    """Cycle phase from the three standardized repressors.

    A complex three-phase order parameter is less dependent on one large peak
    than a single-channel peak detector.
    """
    x = np.asarray(proteins, dtype=float)
    scale = np.std(x, axis=0, ddof=1)
    scale[scale < 1e-9] = 1.0
    z = (x - np.mean(x, axis=0)) / scale
    roots = np.exp(2j * np.pi * np.arange(3) / 3.0)
    order = z @ roots
    phase = np.unwrap(np.angle(order))
    if np.nanmedian(np.diff(phase)) < 0.0:
        phase = -phase

    # Remove small backward jitter without fabricating extra revolutions.
    phase_monotone = np.maximum.accumulate(phase)
    start_level = 2.0 * np.pi * math.ceil(phase_monotone[0] / (2.0 * np.pi))
    stop_level = 2.0 * np.pi * math.floor(phase_monotone[-1] / (2.0 * np.pi))
    levels = np.arange(start_level, stop_level + 0.5 * np.pi, 2.0 * np.pi)
    crossings = []
    for level in levels:
        idx = int(np.searchsorted(phase_monotone, level, side="left"))
        if idx <= 0 or idx >= len(tau):
            continue
        p0, p1 = phase_monotone[idx - 1], phase_monotone[idx]
        fraction = 0.0 if p1 <= p0 else (level - p0) / (p1 - p0)
        crossings.append(float(tau[idx - 1] + fraction * (tau[idx] - tau[idx - 1])))
    return valid_intervals(crossings), np.asarray(crossings), phase


def reporter_proxy(activity, sample_dt):
    """A normalized stable reporter driven by PLtetO1.

    This is a measurement proxy, not an added stochastic circuit reaction.
    Its lifetime is one source protein lifetime in dimensionless units.
    """
    activity = np.asarray(activity, dtype=float)
    reporter = np.zeros_like(activity)
    retention = math.exp(-sample_dt)
    for i in range(1, len(activity)):
        reporter[i] = (
            reporter[i - 1] * retention
            + activity[i - 1] * (1.0 - retention)
        )
    return reporter


def acf_and_psd_periods(signal, sample_dt):
    x = np.asarray(signal, dtype=float)
    x = x - np.mean(x)
    if np.std(x) < 1e-12:
        return float("nan"), float("nan"), np.array([]), np.array([])

    n = len(x)
    fft = np.fft.rfft(x, n=2 * n)
    acf = np.fft.irfft(fft * np.conjugate(fft))[:n]
    denominator = np.arange(n, 0, -1, dtype=float)
    acf = acf / denominator
    acf = acf / max(acf[0], 1e-15)
    lags = np.arange(n, dtype=float) * sample_dt
    mask = (lags >= MIN_PERIOD_LIFETIMES) & (lags <= MAX_PERIOD_LIFETIMES)
    candidates, _ = find_peaks(acf)
    candidates = candidates[mask[candidates]]
    acf_period = (
        float(lags[candidates[np.argmax(acf[candidates])]])
        if candidates.size else float("nan")
    )

    windowed = x * np.hanning(n)
    frequencies = np.fft.rfftfreq(n, d=sample_dt)
    power = np.abs(np.fft.rfft(windowed)) ** 2
    fmask = (
        (frequencies >= 1.0 / MAX_PERIOD_LIFETIMES)
        & (frequencies <= 1.0 / MIN_PERIOD_LIFETIMES)
    )
    if np.any(fmask):
        selected = np.flatnonzero(fmask)
        best = selected[np.argmax(power[selected])]
        psd_period = float(1.0 / frequencies[best])
    else:
        psd_period = float("nan")
    return acf_period, psd_period, lags, acf


def analyze_trajectory(result, p=P, settings=SETTINGS):
    mask = result["tau"] >= settings.transient_lifetimes
    tau = result["tau"][mask]
    proteins = result["proteins_total"][mask].astype(float)
    free = result["proteins_free"][mask].astype(float)
    activities = result["activities"][mask].astype(float)
    pltet_activity = activities[:, CI]  # TetR-controlled cI/phiC31 promoter.

    legacy_periods, legacy_events, legacy_signal = legacy_ci_peak_periods(
        tau, proteins[:, CI], settings
    )
    onset_periods, onset_events = activity_onset_periods(tau, pltet_activity)
    phase_periods, phase_events, phase = cyclic_phase_periods(tau, proteins)
    reporter = reporter_proxy(pltet_activity, settings.sample_dt_lifetimes)
    acf_period, psd_period, acf_lags, acf = acf_and_psd_periods(
        reporter, settings.sample_dt_lifetimes
    )

    return {
        "tau": tau,
        "proteins": proteins,
        "free": free,
        "pltet_activity": pltet_activity,
        "reporter_proxy": reporter,
        "legacy_signal": legacy_signal,
        "phase": phase,
        "acf_lags": acf_lags,
        "acf": acf,
        "legacy_periods": legacy_periods,
        "activity_periods": onset_periods,
        "phase_periods": phase_periods,
        "acf_period": acf_period,
        "psd_period": psd_period,
        "legacy_event_times": legacy_events,
        "activity_event_times": onset_events,
        "phase_event_times": phase_events,
        "n_legacy_cycles": int(len(legacy_periods)),
        "n_activity_cycles": int(len(onset_periods)),
        "n_phase_cycles": int(len(phase_periods)),
    }


def finite_mean(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    return float(np.mean(values)) if values.size else float("nan")


def trajectory_row(condition, index, result, analysis, p=P):
    def mean_gen(periods):
        return finite_mean(periods) * p.generations_per_lifetime

    return {
        "construct_mode": result["construct_mode"],
        "condition": condition,
        "with_sponge": result["with_sponge"],
        "trajectory_index": index,
        "seed": result["seed"],
        "event_count": result["event_count"],
        "legacy_CI_peak_period_gen": mean_gen(analysis["legacy_periods"]),
        "PLtetO1_onset_period_gen": mean_gen(analysis["activity_periods"]),
        "three_repressor_phase_period_gen": mean_gen(analysis["phase_periods"]),
        "reporter_ACF_period_gen": analysis["acf_period"] * p.generations_per_lifetime,
        "reporter_PSD_period_gen": analysis["psd_period"] * p.generations_per_lifetime,
        "n_legacy_cycles": analysis["n_legacy_cycles"],
        "n_activity_cycles": analysis["n_activity_cycles"],
        "n_phase_cycles": analysis["n_phase_cycles"],
        "activity_minus_legacy_cycle_count": (
            analysis["n_activity_cycles"] - analysis["n_legacy_cycles"]
        ),
        "phase_minus_legacy_cycle_count": (
            analysis["n_phase_cycles"] - analysis["n_legacy_cycles"]
        ),
        "TetR_total_peak": float(np.max(analysis["proteins"][:, TETR])),
        "TetR_free_trough": float(np.min(analysis["free"][:, TETR])),
    }


METRIC_COLUMNS = (
    "legacy_CI_peak_period_gen",
    "PLtetO1_onset_period_gen",
    "three_repressor_phase_period_gen",
    "reporter_ACF_period_gen",
    "reporter_PSD_period_gen",
)


def bootstrap_ratio(no_values, sponge_values, replicates, seed):
    no_values = np.asarray(no_values, dtype=float)
    sponge_values = np.asarray(sponge_values, dtype=float)
    no_values = no_values[np.isfinite(no_values)]
    sponge_values = sponge_values[np.isfinite(sponge_values)]
    if no_values.size < 2 or sponge_values.size < 2:
        return float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    ratios = np.empty(replicates, dtype=float)
    for i in range(replicates):
        no_sample = rng.choice(no_values, size=no_values.size, replace=True)
        sp_sample = rng.choice(sponge_values, size=sponge_values.size, replace=True)
        ratios[i] = np.mean(sp_sample) / np.mean(no_sample)
    low, high = np.percentile(ratios, [2.5, 97.5])
    return float(low), float(high)


def condition_summary(rows, construct_mode, condition):
    selected = [
        row for row in rows
        if row["construct_mode"] == construct_mode and row["condition"] == condition
    ]
    summary = {
        "construct_mode": construct_mode,
        "condition": condition,
        "n_trajectories": len(selected),
    }
    for metric in METRIC_COLUMNS:
        values = np.array([row[metric] for row in selected], dtype=float)
        values = values[np.isfinite(values)]
        summary[f"{metric}__n"] = int(values.size)
        summary[f"{metric}__mean"] = finite_mean(values)
        summary[f"{metric}__sd"] = (
            float(np.std(values, ddof=1)) if values.size >= 2 else float("nan")
        )
        summary[f"{metric}__cv"] = (
            float(np.std(values, ddof=1) / np.mean(values))
            if values.size >= 2 and np.mean(values) != 0.0 else float("nan")
        )
    summary["mean_activity_minus_legacy_cycles"] = finite_mean(
        [row["activity_minus_legacy_cycle_count"] for row in selected]
    )
    summary["mean_phase_minus_legacy_cycles"] = finite_mean(
        [row["phase_minus_legacy_cycle_count"] for row in selected]
    )
    return summary


def write_csv(path, rows):
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def plot_diagnostics(output, all_results, all_analyses, rows, ratio_rows):
    if plt is None:
        return

    # Example trajectories: raw CI, promoter activity, reporter proxy and phase.
    fig, axes = plt.subplots(len(all_analyses), 4, figsize=(16, 3.2 * len(all_analyses)))
    if len(all_analyses) == 1:
        axes = np.asarray([axes])
    for row_index, ((mode, condition), analyses) in enumerate(all_analyses.items()):
        a = analyses[0]
        label = f"{mode} | {condition}"
        axes[row_index, 0].plot(a["tau"], a["proteins"][:, CI], lw=0.8)
        axes[row_index, 0].scatter(
            a["legacy_event_times"],
            np.interp(a["legacy_event_times"], a["tau"], a["proteins"][:, CI]),
            s=13, color="tab:red",
        )
        axes[row_index, 0].set_title(label + " | CI peaks")
        axes[row_index, 1].plot(a["tau"], a["pltet_activity"], lw=0.9)
        axes[row_index, 1].axhline(ACTIVITY_ON_THRESHOLD, ls="--", c="k", lw=0.7)
        axes[row_index, 1].set_title("PLtetO1 activity")
        axes[row_index, 2].plot(a["tau"], a["reporter_proxy"], lw=0.9)
        axes[row_index, 2].set_title("reporter proxy")
        axes[row_index, 3].plot(a["tau"], a["phase"] / (2.0 * np.pi), lw=0.9)
        axes[row_index, 3].set_title("unwrapped cyclic phase / 2pi")
        for ax in axes[row_index]:
            ax.set_xlabel("time (protein lifetimes)")
            ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(output / "v54s0_example_estimators.png", dpi=180)
    plt.show()

    # Mean periods and sponge/no-sponge ratios by estimator.
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    x = np.arange(len(METRIC_COLUMNS))
    width = 0.36
    for mode_index, mode in enumerate(CONSTRUCT_MODES):
        mode_rows = [r for r in rows if r["construct_mode"] == mode]
        means_no = [
            finite_mean([r[m] for r in mode_rows if r["condition"] == "no_sponge"])
            for m in METRIC_COLUMNS
        ]
        means_sp = [
            finite_mean([r[m] for r in mode_rows if r["condition"] == "sponge_Nt40_eta1"])
            for m in METRIC_COLUMNS
        ]
        offset = (mode_index - (len(CONSTRUCT_MODES) - 1) / 2.0) * 0.08
        axes[0].bar(x - width / 2 + offset, means_no, width, alpha=0.75, label=f"{mode} no")
        axes[0].bar(x + width / 2 + offset, means_sp, width, alpha=0.75, label=f"{mode} sponge")
        mode_ratios = [
            next(r["ratio"] for r in ratio_rows if r["construct_mode"] == mode and r["metric"] == m)
            for m in METRIC_COLUMNS
        ]
        axes[1].plot(x, mode_ratios, "o-", label=mode)
    labels = ["CI peak", "PLtet onset", "phase", "ACF", "PSD"]
    axes[0].set_xticks(x, labels, rotation=20)
    axes[0].set_ylabel("trajectory-mean period (generations)")
    axes[0].legend(fontsize=8)
    axes[0].grid(axis="y", alpha=0.2)
    axes[1].axhline(1.0, color="black", lw=1)
    axes[1].axhline(1.2, color="grey", lw=0.8, ls="--")
    axes[1].set_xticks(x, labels, rotation=20)
    axes[1].set_ylabel("sponge / no-sponge mean-period ratio")
    axes[1].legend(fontsize=8)
    axes[1].grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(output / "v54s0_period_estimator_comparison.png", dpi=180)
    plt.show()


def main():
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output = OUTPUT_ROOT / f"v54s0_{RUN_MODE}_{timestamp}"
    output.mkdir(parents=True, exist_ok=True)

    print("Loaded: v54s0_paper_aligned_period_estimator_audit_20260815")
    print(f"Run mode: {RUN_MODE}; settings={SETTINGS}")
    print(f"Construct modes: {CONSTRUCT_MODES}")
    print("No parameter is tuned to period lengthening in this version.")

    all_rows = []
    all_results = {}
    all_analyses = {}
    wall_start = time.perf_counter()

    for mode_index, construct_mode in enumerate(CONSTRUCT_MODES):
        conditions = (
            ("no_sponge", False, P),
            (
                "sponge_Nt40_eta1",
                True,
                replace(P, mean_sponge_plasmids=40, effective_sponge_fraction=1.0),
            ),
        )
        for condition_index, (condition, with_sponge, parameters) in enumerate(conditions):
            key = (construct_mode, condition)
            all_results[key] = []
            all_analyses[key] = []
            for trajectory_index in range(SETTINGS.trajectories_per_condition):
                seed = (
                    SETTINGS.random_seed
                    + mode_index * 1_000_000
                    + condition_index * 100_000
                    + trajectory_index
                )
                print(
                    f"[{construct_mode} | {condition}] "
                    f"trajectory {trajectory_index + 1}/"
                    f"{SETTINGS.trajectories_per_condition}"
                )
                result = gillespie_trajectory(
                    with_sponge, seed, construct_mode, parameters, SETTINGS
                )
                analysis = analyze_trajectory(result, parameters, SETTINGS)
                all_results[key].append(result)
                all_analyses[key].append(analysis)
                all_rows.append(
                    trajectory_row(
                        condition, trajectory_index, result, analysis, parameters
                    )
                )

    summaries = []
    ratio_rows = []
    for mode_index, construct_mode in enumerate(CONSTRUCT_MODES):
        no_summary = condition_summary(all_rows, construct_mode, "no_sponge")
        sp_summary = condition_summary(all_rows, construct_mode, "sponge_Nt40_eta1")
        summaries.extend([no_summary, sp_summary])
        for metric_index, metric in enumerate(METRIC_COLUMNS):
            no_values = [
                row[metric] for row in all_rows
                if row["construct_mode"] == construct_mode
                and row["condition"] == "no_sponge"
            ]
            sp_values = [
                row[metric] for row in all_rows
                if row["construct_mode"] == construct_mode
                and row["condition"] == "sponge_Nt40_eta1"
            ]
            ratio = finite_mean(sp_values) / finite_mean(no_values)
            low, high = bootstrap_ratio(
                no_values,
                sp_values,
                SETTINGS.bootstrap_replicates,
                SETTINGS.random_seed + 9_000_000 + mode_index * 100 + metric_index,
            )
            ratio_rows.append(
                {
                    "construct_mode": construct_mode,
                    "metric": metric,
                    "no_sponge_mean_gen": finite_mean(no_values),
                    "sponge_mean_gen": finite_mean(sp_values),
                    "ratio": ratio,
                    "ratio_ci_low": low,
                    "ratio_ci_high": high,
                    "supports_lengthening_ratio_gt_1": bool(low > 1.0),
                }
            )

    write_csv(output / "v54s0_trajectory_metrics.csv", all_rows)
    write_csv(output / "v54s0_condition_summary.csv", summaries)
    write_csv(output / "v54s0_period_ratio_bootstrap.csv", ratio_rows)

    metadata = {
        "version": "v54s0_paper_aligned_period_estimator_audit_20260815",
        "run_mode": RUN_MODE,
        "settings": SETTINGS.__dict__,
        "parameters": P.__dict__,
        "construct_modes": list(CONSTRUCT_MODES),
        "analysis": {
            "min_period_lifetimes": MIN_PERIOD_LIFETIMES,
            "max_period_lifetimes": MAX_PERIOD_LIFETIMES,
            "activity_on_threshold": ACTIVITY_ON_THRESHOLD,
            "activity_off_threshold": ACTIVITY_OFF_THRESHOLD,
        },
        "wall_time_seconds": time.perf_counter() - wall_start,
    }
    (output / "v54s0_metadata.json").write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    plot_diagnostics(output, all_results, all_analyses, all_rows, ratio_rows)

    print(f"\nSaved v54s0 outputs to: {output}")
    print("\nSponge/no-sponge period ratios (trajectory is the inference unit):")
    for row in ratio_rows:
        print(
            f"{row['construct_mode']:>14} | {row['metric']:<38} | "
            f"no={row['no_sponge_mean_gen']:.3f} gen | "
            f"sp={row['sponge_mean_gen']:.3f} gen | "
            f"ratio={row['ratio']:.3f} "
            f"[{row['ratio_ci_low']:.3f}, {row['ratio_ci_high']:.3f}] | "
            f"lengthening={row['supports_lengthening_ratio_gt_1']}"
        )
    print("\nInterpretation:")
    print("- If ACF/phase/onset lengthen but legacy CI peaks do not, the old estimator masked the effect.")
    print("- If every estimator remains near or below 1, the SI reaction model does not reproduce the paired experimental shift.")
    print("- Do not tune Nt or eta in v54s0; move next to explicit division and independently constrained thresholds.")



if __name__ == "__main__":
    main()
