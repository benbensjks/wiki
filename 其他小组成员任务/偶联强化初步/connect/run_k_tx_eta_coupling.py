import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
GAMMA_TEST_DIR = ROOT / "gamma_coupling_test"
sys.path.insert(0, str(GAMMA_TEST_DIR))

import run_gamma_coupling as base


BASE_K_TX = base.P_OSC.k_tx
BASE_ETA = base.P_OSC.eta
BASE_S_EFF = base.P_OSC.S_eff


def reset_oscillator_params():
    base.P_OSC.k_tx = BASE_K_TX
    base.P_OSC.eta = BASE_ETA
    base.P_OSC.S_eff = BASE_S_EFF


def set_k_tx_factor(factor):
    reset_oscillator_params()
    base.P_OSC.k_tx = BASE_K_TX * factor
    return base.P_OSC.k_tx


def set_eta_factor(factor):
    reset_oscillator_params()
    base.P_OSC.eta = BASE_ETA * factor
    base.P_OSC.S_eff = base.P_OSC.S_phys * base.P_OSC.eta
    return base.P_OSC.eta


def plot_case(case_dir, param_label, factor, value, waveform, counter, metrics, features):
    fig, ax = plt.subplots(4, 1, figsize=(12, 10), sharex=False)
    ax[0].plot(waveform["time_min"], waveform["C31_protein"], color="tab:purple", lw=1)
    ax[0].set_ylabel("C31 copies")
    ax[0].set_title(f"Oscillator output, {param_label} = {factor:.2f}x ({value:.3g})")

    ax[1].plot(
        counter["time_h"],
        counter["input_C31_uM_baseline_corrected"],
        color="tab:pink",
        lw=1,
    )
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
    fig.savefig(case_dir / f"{param_label}_coupled_counter.png", dpi=180)
    plt.close(fig)


def write_case_readme(case_dir, param_label, factor, value, features, metrics):
    with open(case_dir / "README.md", "w", encoding="utf-8") as f:
        f.write(f"# {param_label} coupling test: {factor:.2f}x\n\n")
        if param_label == "k_tx":
            f.write("This case changes the C31 output transcription term in the oscillator equation:\n\n")
            f.write("`dC31_mRNA/dt = N_C * k_tx * PLtetO1 - gamma_m * C31_mRNA`.\n\n")
        else:
            f.write("This case changes the TetR sponge accessibility parameter in the oscillator model:\n\n")
            f.write("`S_eff = S_phys * eta`, which changes TetR buffering and therefore PLtetO1/C31 timing.\n\n")
        f.write("The generated C31 waveform is used as the time-dependent input for the RDF counter model.\n\n")
        f.write(f"- {param_label}_factor: {factor}\n")
        f.write(f"- {param_label}_value: {value}\n")
        f.write(f"- period_min: {features['period_min']}\n")
        f.write(f"- period_generations: {features['period_generations']}\n")
        f.write(f"- c31_peak_uM: {features['c31_peak_uM']}\n")
        f.write(f"- duty_above_half_peak: {features['duty_above_half_peak']}\n")
        for key, val in metrics.items():
            f.write(f"- {key}: {val}\n")


def run_sweep(param_label, setter, out_root):
    out_root.mkdir(parents=True, exist_ok=True)
    summaries = []
    for factor in [0.8, 1.0, 1.2]:
        value = setter(factor)
        case_dir = out_root / f"{param_label}_{factor:.1f}x"
        case_dir.mkdir(parents=True, exist_ok=True)

        waveform = base.simulate_oscillator(1.0)
        features = base.waveform_features(waveform)
        counter, metrics, _samples = base.simulate_counter(waveform, features)

        waveform.to_csv(case_dir / "oscillator_waveform.csv", index=False)
        counter.to_csv(case_dir / "counter_timecourse.csv", index=False)
        plot_case(case_dir, param_label, factor, value, waveform, counter, metrics, features)

        row = {
            "case": f"{param_label}_{factor:.1f}x",
            f"{param_label}_factor": factor,
            f"{param_label}_value": value,
            "period_min": features["period_min"],
            "period_generations": features["period_generations"],
            "c31_peak_uM": features["c31_peak_uM"],
            "c31_trough_copies": features["c31_trough_copies"],
            "duty_above_half_peak": features["duty_above_half_peak"],
            **metrics,
        }
        pd.DataFrame([row]).to_csv(case_dir / "summary_metrics.csv", index=False)
        write_case_readme(case_dir, param_label, factor, value, features, metrics)
        summaries.append(row)

    summary = pd.DataFrame(summaries)
    summary.to_csv(out_root / f"{param_label}_sweep_summary.csv", index=False)
    reset_oscillator_params()
    return summary


def main():
    k_tx_summary = run_sweep("k_tx", set_k_tx_factor, ROOT / "k_tx_coupling_test")
    eta_summary = run_sweep("eta", set_eta_factor, ROOT / "eta_coupling_test")

    print("\n[k_tx]")
    print(k_tx_summary.to_string(index=False))
    print("\n[eta]")
    print(eta_summary.to_string(index=False))


if __name__ == "__main__":
    main()
