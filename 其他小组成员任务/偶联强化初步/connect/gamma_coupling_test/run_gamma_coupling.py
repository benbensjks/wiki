import math
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.signal import find_peaks

ROOT = Path(__file__).resolve().parents[1]
RDF_DIR = ROOT / "rdfmodel"
sys.path.insert(0, str(RDF_DIR))

import zhao_core as zc


class OscillatorParams:
    lambd = 1000
    K = 13
    n = 3
    Td = 30
    gamma = math.log(2) / Td

    S_phys = 80
    eta = 0.5
    S_eff = S_phys * eta
    Kd = 6
    n_s = 1.5

    K_R = 6
    n_R = 2
    H_leak = 0.05

    N_C = 10
    k_tx = 5
    gamma_m = math.log(2) / 2
    beta = 0.5
    mu = math.log(2) / 30
    copies_per_uM = 602


P_OSC = OscillatorParams()


def hill_rep(x):
    return (P_OSC.K ** P_OSC.n) / (P_OSC.K ** P_OSC.n + x ** P_OSC.n + 1e-12)


def plteto_activity(tet_r_free):
    return P_OSC.H_leak + (1 - P_OSC.H_leak) / (
        1 + (tet_r_free / P_OSC.K_R) ** P_OSC.n_R
    )


def free_tetr(tet_r_total, sponge_capacity):
    if tet_r_total <= 0:
        return 0.0
    low, high = 0.0, tet_r_total
    for _ in range(100):
        tf = (low + high) / 2
        bound = sponge_capacity * (tf ** P_OSC.n_s) / (
            tf ** P_OSC.n_s + P_OSC.Kd ** P_OSC.n_s
        )
        if tf + bound > tet_r_total:
            high = tf
        else:
            low = tf
    return (low + high) / 2


def make_oscillator_model(gamma_factor, sponge_capacity):
    gamma = P_OSC.gamma * gamma_factor

    def ode(_t, y):
        tet_r, ci, lac_i, c31_mrna, c31 = y
        tet_r_free = free_tetr(tet_r, sponge_capacity)

        d_tet_r = gamma * (P_OSC.lambd * hill_rep(lac_i) - tet_r)
        d_ci = gamma * (P_OSC.lambd * hill_rep(tet_r_free) - ci)
        d_lac_i = gamma * (P_OSC.lambd * hill_rep(ci) - lac_i)

        promoter = plteto_activity(tet_r_free)
        d_c31_mrna = P_OSC.N_C * P_OSC.k_tx * promoter - P_OSC.gamma_m * c31_mrna
        d_c31 = P_OSC.beta * c31_mrna - P_OSC.mu * c31

        return [d_tet_r, d_ci, d_lac_i, d_c31_mrna, d_c31]

    return ode


def simulate_oscillator(gamma_factor):
    t_eval = np.linspace(0, 6000, 6000)
    y0 = [0.5, 0.2, 1.0, 0.0, 0.0]
    sol = solve_ivp(
        make_oscillator_model(gamma_factor, P_OSC.S_eff),
        (0, 6000),
        y0,
        t_eval=t_eval,
        method="BDF",
        rtol=1e-6,
        atol=1e-8,
    )
    if not sol.success:
        raise RuntimeError(sol.message)

    tet_r = sol.y[0]
    tet_r_free = np.array([free_tetr(x, P_OSC.S_eff) for x in tet_r])
    promoter = np.array([plteto_activity(x) for x in tet_r_free])
    buffering = (tet_r - tet_r_free) / np.maximum(tet_r, 1)

    return pd.DataFrame(
        {
            "time_min": sol.t,
            "TetR_total": tet_r,
            "TetR_free": tet_r_free,
            "CI": sol.y[1],
            "LacI": sol.y[2],
            "PLtetO1_activity": promoter,
            "C31_mRNA": sol.y[3],
            "C31_protein": sol.y[4],
            "sponge_buffering": buffering,
        }
    )


def waveform_features(df):
    t = df["time_min"].to_numpy()
    c31 = df["C31_protein"].to_numpy()
    steady = t >= 1500
    c31_steady = c31[steady]
    t_steady = t[steady]
    peaks, _ = find_peaks(c31_steady, prominence=np.ptp(c31_steady) * 0.1, distance=100)
    troughs, _ = find_peaks(-c31_steady, prominence=np.ptp(c31_steady) * 0.1, distance=100)
    peak_t = t_steady[peaks]
    trough_t = t_steady[troughs]
    period = float(np.mean(np.diff(peak_t))) if len(peak_t) >= 2 else np.nan
    duty = float(np.mean(c31_steady > 0.5 * np.nanmax(c31_steady)))
    return {
        "period_min": period,
        "period_h": period / 60 if period == period else np.nan,
        "period_generations": period / 30 if period == period else np.nan,
        "c31_peak_copies": float(np.nanmax(c31_steady)),
        "c31_trough_copies": float(np.nanmin(c31_steady)),
        "c31_peak_uM": float(np.nanmax(c31_steady) / P_OSC.copies_per_uM),
        "duty_above_half_peak": duty,
        "peak_times_h": peak_t / 60,
        "trough_times_h": trough_t / 60,
    }


def make_counter_source(df, k_dil, baseline_sub=0.26, scale=1.0):
    t_h = df["time_min"].to_numpy() / 60
    c31_uM = df["C31_protein"].to_numpy() / P_OSC.copies_per_uM
    spl = CubicSpline(t_h, c31_uM, bc_type="natural")
    dspl = spl.derivative()

    def source(t):
        c = float(spl(t))
        w = max(0.0, c - baseline_sub)
        dwdt = float(dspl(t)) if c > baseline_sub else 0.0
        return scale * max(0.0, dwdt + k_dil * w)

    return source, spl


def toggle_score(t, lr_fraction, troughs_h, skip=3):
    troughs_h = np.array(troughs_h)
    troughs_h = troughs_h[troughs_h <= t[-1]]
    if len(troughs_h) < skip + 4:
        return 0.0, np.nan, np.nan, 0.0, []
    samples = np.interp(troughs_h, t, lr_fraction)[skip:]
    high = samples > 0.6
    low = samples < 0.4
    mid = ~(high | low)
    alt = 0
    total = 0
    for i in range(len(samples) - 1):
        if mid[i] or mid[i + 1]:
            continue
        total += 1
        if high[i] != high[i + 1]:
            alt += 1
    fidelity = alt / total if total else 0.0
    h_mean = samples[high].mean() if high.any() else np.nan
    l_mean = samples[low].mean() if low.any() else np.nan
    core = min(np.nan_to_num(h_mean, nan=0.0), 1 - np.nan_to_num(l_mean, nan=1.0))
    return core * fidelity, h_mean, l_mean, fidelity, list(zip(troughs_h[skip:], samples))


def simulate_counter(df, features):
    p = zc.default_params()
    p.update(K_rep=0.0186, n_rep=3.4, krep_tsl=15.0, krdf_tsl=200.0, k_tag_int=12.0)
    c = zc._rate_constants(p)
    source, input_spline = make_counter_source(df, p["k_dil"], baseline_sub=0.26, scale=1.0)
    y0 = zc.y0_PB(
        p,
        rep_mrna=p["k_tscr"] * p["Dtot"] / p["k_rna"],
        rep=15.0 * p["k_tscr"] * p["Dtot"] / p["k_rna"] / p["k_dil"],
    )
    t_eval = np.arange(0, 100, 0.05)
    sol = solve_ivp(
        lambda t, y: zc.rhs(t, y, p, c, source),
        (0, 100),
        y0,
        method="LSODA",
        rtol=1e-6,
        atol=1e-11,
        max_step=0.1,
        t_eval=t_eval,
    )
    if not sol.success:
        raise RuntimeError(sol.message)

    y = sol.y.T
    lr_fraction = zc.LR_total(y) / p["Dtot"]
    pb_fraction = zc.PB_total(y) / p["Dtot"]
    int_total = zc.int_total(y)
    rdf_total = zc.rdf_total(y)
    rep = y[:, 37]
    input_uM = np.maximum(0.0, input_spline(sol.t) - 0.26)
    score, h_mean, l_mean, fidelity, samples = toggle_score(
        sol.t, lr_fraction, features["trough_times_h"]
    )

    out = pd.DataFrame(
        {
            "time_h": sol.t,
            "input_C31_uM_baseline_corrected": input_uM,
            "LR_fraction": lr_fraction,
            "PB_fraction": pb_fraction,
            "Int_total_uM": int_total,
            "RDF_total_uM": rdf_total,
            "BM3R1_uM": rep,
        }
    )
    metrics = {
        "toggle_score": score,
        "H_mean": h_mean,
        "L_mean": l_mean,
        "alternation_fidelity": fidelity,
        "Int_peak_uM": float(np.nanmax(int_total)),
        "Int_final_uM": float(int_total[-1]),
        "RDF_peak_uM": float(np.nanmax(rdf_total)),
        "BM3R1_peak_uM": float(np.nanmax(rep)),
        "sample_count": len(samples),
    }
    return out, metrics, samples


def plot_case(case_dir, gamma_factor, waveform, counter, metrics, features):
    fig, ax = plt.subplots(4, 1, figsize=(12, 10), sharex=False)
    ax[0].plot(waveform["time_min"], waveform["C31_protein"], color="tab:purple", lw=1)
    ax[0].set_ylabel("C31 copies")
    ax[0].set_title(f"Oscillator output, gamma = {gamma_factor:.2f}x")
    ax[1].plot(counter["time_h"], counter["input_C31_uM_baseline_corrected"], color="tab:pink", lw=1)
    ax[1].set_ylabel("Input uM")
    ax[2].plot(counter["time_h"], counter["LR_fraction"], label="LR", color="tab:red", lw=1.2)
    ax[2].plot(counter["time_h"], counter["PB_fraction"], label="PB", color="tab:blue", lw=1.2)
    ax[2].set_ylabel("DNA fraction")
    ax[2].set_ylim(-0.05, 1.05)
    ax[2].legend()
    ax[3].plot(counter["time_h"], counter["RDF_total_uM"], label="RDF", color="black", lw=1)
    ax[3].plot(counter["time_h"], counter["BM3R1_uM"], label="BM3R1", color="tab:cyan", lw=1)
    ax[3].set_ylabel("uM")
    ax[3].set_xlabel("time (h)")
    ax[3].legend()
    for axis in ax:
        axis.grid(alpha=0.3)
    fig.suptitle(
        f"period={features['period_min']:.1f} min, "
        f"score={metrics['toggle_score']:.3f}, fidelity={metrics['alternation_fidelity']:.2f}"
    )
    fig.tight_layout()
    fig.savefig(case_dir / "gamma_coupled_counter.png", dpi=180)
    plt.close(fig)


def write_case_summary(case_dir, gamma_factor, features, metrics):
    rows = {
        "gamma_factor": gamma_factor,
        "gamma_value_per_min": P_OSC.gamma * gamma_factor,
        "period_min": features["period_min"],
        "period_generations": features["period_generations"],
        "c31_peak_copies": features["c31_peak_copies"],
        "c31_peak_uM": features["c31_peak_uM"],
        "c31_trough_copies": features["c31_trough_copies"],
        "duty_above_half_peak": features["duty_above_half_peak"],
        **metrics,
    }
    pd.DataFrame([rows]).to_csv(case_dir / "summary_metrics.csv", index=False)
    with open(case_dir / "README.md", "w", encoding="utf-8") as f:
        f.write(f"# Gamma coupling test: {gamma_factor:.2f}x\n\n")
        f.write("This case changes the repressilator gamma term in the oscillator equations:\n\n")
        f.write("`dTetR/dt`, `dCI/dt`, and `dLacI/dt` use `gamma = gamma_base * factor`.\n\n")
        f.write("The generated C31 waveform is then used as the time-dependent input for the RDF counter model.\n\n")
        for key, value in rows.items():
            f.write(f"- {key}: {value}\n")


def main():
    out_root = Path(__file__).resolve().parent
    factors = [0.8, 1.0, 1.2]
    summaries = []

    for factor in factors:
        case_dir = out_root / f"gamma_{factor:.1f}x"
        case_dir.mkdir(parents=True, exist_ok=True)
        waveform = simulate_oscillator(factor)
        features = waveform_features(waveform)
        counter, metrics, _samples = simulate_counter(waveform, features)

        waveform.to_csv(case_dir / "oscillator_waveform.csv", index=False)
        counter.to_csv(case_dir / "counter_timecourse.csv", index=False)
        plot_case(case_dir, factor, waveform, counter, metrics, features)
        write_case_summary(case_dir, factor, features, metrics)

        summaries.append(
            {
                "case": f"gamma_{factor:.1f}x",
                "gamma_factor": factor,
                "period_min": features["period_min"],
                "period_generations": features["period_generations"],
                "c31_peak_uM": features["c31_peak_uM"],
                "duty_above_half_peak": features["duty_above_half_peak"],
                **metrics,
            }
        )

    summary = pd.DataFrame(summaries)
    summary.to_csv(out_root / "gamma_sweep_summary.csv", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
