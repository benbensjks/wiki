#!/usr/bin/env python3
"""Full continuous biochemical ODE reconstruction for a three-bit cascade.

The base model follows ``前置生化背景/初步.pdf`` with eleven dynamic
variables per orthogonal recombinase module:

    (M_I, I_u, I, M_T, T_u, T, M_R, R_u, R, C, S)

Three modules therefore contribute 33 states.  Two additional continuous
carry-memory proteins E0 and E1 implement an explicit prospective falling-edge
detector.  They are biochemical ODE states, not event handlers.  The model
contains no discrete logical bits, threshold-triggered resets, event queues,
forced toggles, or direction-selective locking.

All time units are minutes and concentrations are micromolar (uM), except S,
U and promoter activities, which are dimensionless.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks


STATE_NAMES = ("M_I", "I_u", "I", "M_T", "T_u", "T", "M_R", "R_u", "R", "C", "S")
N_BASE = len(STATE_NAMES)
N_BITS = 3
EDGE_NAMES = ("E0", "E1")


@dataclass(frozen=True)
class Parameters:
    # Cell growth and transcriptional input.
    doubling_time: float = 50.0
    K_U: float = 0.35
    n_U: float = 3.0

    # Integrase expression.
    alpha_I_0: float = 1.0e-5
    alpha_I: float = 0.030
    delta_mI: float = 0.12
    k_tl_I: float = 0.80
    k_mat_I: float = 0.10
    delta_Iu: float = 0.010
    delta_I: float = 0.003

    # TetR expression.
    alpha_T_0: float = 1.0e-5
    alpha_T: float = 0.040
    delta_mT: float = 0.12
    k_tl_T: float = 0.80
    k_mat_T: float = 0.10
    delta_Tu: float = 0.010
    delta_T: float = 0.020

    # RDF expression from the LR-oriented promoter.
    alpha_R_0: float = 1.0e-6
    alpha_R: float = 0.018
    a_R_PB: float = 0.002
    a_R_LR: float = 1.0
    K_T: float = 0.20
    n_T: float = 3.0
    delta_mR: float = 0.12
    k_tl_R: float = 0.55
    k_mat_R: float = 0.06
    delta_Ru: float = 0.006
    delta_R: float = 0.004

    # Explicit Integrase-RDF complex.  nu=1 implements the PDF's initial
    # coarse-grained stoichiometric choice.
    nu: float = 1.0
    k_on_C: float = 1.20
    k_off_C: float = 0.030
    delta_C: float = 0.003

    # DNA recombination.
    k_rec_fwd: float = 0.12
    k_rec_rev: float = 0.10
    K_fwd: float = 0.25
    K_rev: float = 0.25
    n_fwd: float = 2.0
    n_rev: float = 2.0

    # Smooth external clock: Gaussian transcriptional input pulses.
    clock_period: float = 480.0
    clock_sigma: float = 24.0
    clock_amplitude: float = 1.0
    first_clock: float = 180.0
    n_clock_pulses: int = 12

    # Continuous biochemical carry extension.  E is a slowly clearing protein
    # induced by LR state.  Output requires E high and the current S low, an
    # incoherent feed-forward falling-edge detector.
    alpha_E: float = 0.055
    delta_E: float = 0.022
    K_E: float = 0.30
    n_E: float = 3.0
    K_S_edge: float = 0.25
    n_S_edge: float = 4.0
    carry_gain_01: float = 1.30
    carry_gain_12: float = 1.30

    # Numerical integration.
    rtol: float = 1.0e-8
    atol: float = 1.0e-10
    max_step: float = 1.0

    @property
    def mu(self) -> float:
        return math.log(2.0) / self.doubling_time

    @property
    def t_end(self) -> float:
        return self.first_clock + (self.n_clock_pulses + 1) * self.clock_period


def hill_activation(x: float | np.ndarray, K: float, n: float):
    x = np.maximum(x, 0.0)
    return np.power(x, n) / (np.power(K, n) + np.power(x, n) + 1.0e-30)


def hill_repression(x: float | np.ndarray, K: float, n: float):
    x = np.maximum(x, 0.0)
    return np.power(K, n) / (np.power(K, n) + np.power(x, n) + 1.0e-30)


def bit_slice(bit: int) -> slice:
    return slice(bit * N_BASE, (bit + 1) * N_BASE)


def clock_centres(p: Parameters) -> np.ndarray:
    return p.first_clock + np.arange(p.n_clock_pulses, dtype=float) * p.clock_period


def external_clock(t: float | np.ndarray, p: Parameters):
    """Smooth, event-free external input U0(t)."""
    t_arr = np.asarray(t, dtype=float)
    centres = clock_centres(p)
    z = (np.expand_dims(t_arr, axis=-1) - centres) / p.clock_sigma
    value = p.clock_amplitude * np.exp(-0.5 * z * z).sum(axis=-1)
    return float(value) if np.ndim(t) == 0 else value


def carry_output(E: float, S: float, gain: float, p: Parameters) -> float:
    """Continuous prospective falling-edge signal.

    E stores the recent LR state.  The AND of high E and low current S creates
    a finite pulse after LR->PB.  This is an explicit biochemical hypothesis,
    not part of the original 11-state equations.
    """
    memory_high = hill_activation(E, p.K_E, p.n_E)
    state_low = hill_repression(S, p.K_S_edge, p.n_S_edge)
    return float(gain * memory_high * state_low)


def module_rhs(yb: np.ndarray, U: float, p: Parameters) -> np.ndarray:
    """The eleven equations for one orthogonal module."""
    M_I, I_u, I, M_T, T_u, T, M_R, R_u, R, C, S = np.maximum(yb, 0.0)
    S = float(np.clip(S, 0.0, 1.0))
    mu = p.mu
    F_U = hill_activation(U, p.K_U, p.n_U)

    dM_I = p.alpha_I_0 + p.alpha_I * F_U - (p.delta_mI + mu) * M_I
    dI_u = p.k_tl_I * M_I - (p.k_mat_I + p.delta_Iu + mu) * I_u

    bind = p.k_on_C * I * (R ** p.nu)
    unbind = p.k_off_C * C
    dI = p.k_mat_I * I_u - (p.delta_I + mu) * I - bind + unbind

    dM_T = p.alpha_T_0 + p.alpha_T * F_U - (p.delta_mT + mu) * M_T
    dT_u = p.k_tl_T * M_T - (p.k_mat_T + p.delta_Tu + mu) * T_u
    dT = p.k_mat_T * T_u - (p.delta_T + mu) * T

    A_R = p.a_R_PB * (1.0 - S) + p.a_R_LR * S
    F_T = hill_repression(T, p.K_T, p.n_T)
    dM_R = p.alpha_R_0 + p.alpha_R * A_R * F_T - (p.delta_mR + mu) * M_R
    dR_u = p.k_tl_R * M_R - (p.k_mat_R + p.delta_Ru + mu) * R_u
    dR = p.k_mat_R * R_u - (p.delta_R + mu) * R - p.nu * bind + p.nu * unbind
    dC = bind - (p.k_off_C + p.delta_C + mu) * C

    H_fwd = hill_activation(I, p.K_fwd, p.n_fwd)
    H_rev = hill_activation(C, p.K_rev, p.n_rev)
    dS = p.k_rec_fwd * H_fwd * (1.0 - S) - p.k_rec_rev * H_rev * S

    return np.array((dM_I, dI_u, dI, dM_T, dT_u, dT, dM_R, dR_u, dR, dC, dS))


def full_rhs(t: float, y: np.ndarray, p: Parameters) -> np.ndarray:
    """33-state PDF model plus two continuous carry-memory states."""
    S0 = float(np.clip(y[bit_slice(0)][10], 0.0, 1.0))
    S1 = float(np.clip(y[bit_slice(1)][10], 0.0, 1.0))
    E0 = max(float(y[N_BASE * N_BITS]), 0.0)
    E1 = max(float(y[N_BASE * N_BITS + 1]), 0.0)

    U0 = external_clock(t, p)
    U1 = carry_output(E0, S0, p.carry_gain_01, p)
    U2 = carry_output(E1, S1, p.carry_gain_12, p)

    dy = np.zeros_like(y)
    dy[bit_slice(0)] = module_rhs(y[bit_slice(0)], U0, p)
    dy[bit_slice(1)] = module_rhs(y[bit_slice(1)], U1, p)
    dy[bit_slice(2)] = module_rhs(y[bit_slice(2)], U2, p)

    # E proteins are produced in the LR state and clear continuously.  Their
    # joint action with low S creates an analog edge pulse without events.
    dy[N_BASE * N_BITS] = p.alpha_E * S0 - (p.delta_E + p.mu) * E0
    dy[N_BASE * N_BITS + 1] = p.alpha_E * S1 - (p.delta_E + p.mu) * E1
    return dy


def initial_state() -> np.ndarray:
    return np.zeros(N_BASE * N_BITS + len(EDGE_NAMES), dtype=float)


def simulate(p: Parameters):
    t_eval = np.arange(0.0, p.t_end + 0.5, 1.0)
    sol = solve_ivp(
        lambda t, y: full_rhs(t, y, p),
        (0.0, p.t_end),
        initial_state(),
        method="BDF",
        t_eval=t_eval,
        rtol=p.rtol,
        atol=p.atol,
        max_step=p.max_step,
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol


def trajectories(sol, p: Parameters) -> pd.DataFrame:
    data: dict[str, np.ndarray] = {"t_min": sol.t, "U0": external_clock(sol.t, p)}
    for bit in range(N_BITS):
        ys = sol.y[bit_slice(bit)]
        for j, name in enumerate(STATE_NAMES):
            data[f"bit{bit}_{name}"] = ys[j]

    E0 = sol.y[N_BASE * N_BITS]
    E1 = sol.y[N_BASE * N_BITS + 1]
    S0 = data["bit0_S"]
    S1 = data["bit1_S"]
    data["E0"] = E0
    data["E1"] = E1
    data["U1"] = np.array([carry_output(e, s, p.carry_gain_01, p) for e, s in zip(E0, S0)])
    data["U2"] = np.array([carry_output(e, s, p.carry_gain_12, p) for e, s in zip(E1, S1)])
    return pd.DataFrame(data)


def hysteresis_decode(values: np.ndarray, low: float = 0.30, high: float = 0.70) -> np.ndarray:
    state = 0
    out = np.zeros(values.size, dtype=int)
    for i, value in enumerate(values):
        if state == 0 and value >= high:
            state = 1
        elif state == 1 and value <= low:
            state = 0
        out[i] = state
    return out


def analyse(df: pd.DataFrame, p: Parameters):
    for bit in range(N_BITS):
        df[f"q{bit}"] = hysteresis_decode(df[f"bit{bit}_S"].to_numpy())
    df["code"] = df["q0"] + 2 * df["q1"] + 4 * df["q2"]

    centres = clock_centres(p)
    rows = []
    for k, centre in enumerate(centres):
        sample_t = centre + 0.45 * p.clock_period
        idx = int(np.argmin(np.abs(df["t_min"].to_numpy() - sample_t)))
        rows.append(
            {
                "clock_index": k,
                "clock_centre_min": centre,
                "sample_time_min": float(df.loc[idx, "t_min"]),
                "S0": float(df.loc[idx, "bit0_S"]),
                "S1": float(df.loc[idx, "bit1_S"]),
                "S2": float(df.loc[idx, "bit2_S"]),
                "q0": int(df.loc[idx, "q0"]),
                "q1": int(df.loc[idx, "q1"]),
                "q2": int(df.loc[idx, "q2"]),
                "code": int(df.loc[idx, "code"]),
            }
        )
    samples = pd.DataFrame(rows)

    pulse_rows = []
    for name in ("U0", "U1", "U2"):
        signal = df[name].to_numpy()
        peaks, props = find_peaks(signal, height=0.05, distance=int(0.35 * p.clock_period))
        for number, peak in enumerate(peaks):
            pulse_rows.append(
                {
                    "signal": name,
                    "pulse_index": number,
                    "peak_time_min": float(df.loc[peak, "t_min"]),
                    "peak_amplitude": float(signal[peak]),
                }
            )
    pulses = pd.DataFrame(pulse_rows)

    pulse_periods = {}
    for name in ("U0", "U1", "U2"):
        times = pulses.loc[pulses["signal"] == name, "peak_time_min"].to_numpy()
        intervals = np.diff(times)
        pulse_periods[name] = {
            "n_peaks": int(times.size),
            "median_interval_min": float(np.median(intervals)) if intervals.size else None,
            "mean_interval_min": float(np.mean(intervals)) if intervals.size else None,
        }

    # Each row is sampled after a clock pulse, so the ideal counter has
    # advanced once at the first row.
    expected = (np.arange(len(samples)) + 1) % 8
    observed = samples["code"].to_numpy()
    accuracy = float(np.mean(observed == expected)) if len(observed) else float("nan")
    summary = {
        "solver_success": True,
        "n_base_states": N_BASE * N_BITS,
        "n_carry_states": len(EDGE_NAMES),
        "n_total_states": int(sol_state_count(df)),
        "uses_discrete_events": False,
        "uses_logical_state_variables": False,
        "uses_forced_toggles": False,
        "sampled_codes": observed.tolist(),
        "ideal_binary_codes": expected.tolist(),
        "binary_sequence_accuracy": accuracy,
        "carry_pulse_periods": pulse_periods,
        "minimum_dynamic_value": float(
            df[[f"bit{bit}_{name}" for bit in range(N_BITS) for name in STATE_NAMES] + list(EDGE_NAMES)]
            .min()
            .min()
        ),
        "all_values_finite": bool(np.isfinite(df.select_dtypes(include=[np.number]).to_numpy()).all()),
        "state_ranges": {
            f"S{bit}": [float(df[f"bit{bit}_S"].min()), float(df[f"bit{bit}_S"].max())]
            for bit in range(N_BITS)
        },
    }
    return samples, pulses, summary


def sol_state_count(df: pd.DataFrame) -> int:
    return N_BASE * N_BITS + len(EDGE_NAMES)


def make_figures(df: pd.DataFrame, samples: pd.DataFrame, out_dir: Path):
    t_h = df["t_min"] / 60.0

    fig, axes = plt.subplots(4, 1, figsize=(14, 11), sharex=True)
    axes[0].plot(t_h, df["U0"], label="U0 external clock", color="#111111", lw=1.2)
    axes[0].plot(t_h, df["U1"], label="U1 continuous carry", color="#0072B2")
    axes[0].plot(t_h, df["U2"], label="U2 continuous carry", color="#D55E00")
    axes[0].set_ylabel("Input activity")
    axes[0].legend(ncol=3, loc="upper right")

    colors = ("#009E73", "#0072B2", "#D55E00")
    for bit, color in enumerate(colors):
        axes[1].plot(t_h, df[f"bit{bit}_S"], label=f"S{bit}", color=color, lw=1.4)
    axes[1].axhline(0.30, color="0.6", ls="--", lw=0.8)
    axes[1].axhline(0.70, color="0.6", ls="--", lw=0.8)
    axes[1].set_ylabel("LR fraction S")
    axes[1].set_ylim(-0.03, 1.03)
    axes[1].legend(ncol=3, loc="upper right")

    for bit, color in enumerate(colors):
        axes[2].plot(t_h, df[f"bit{bit}_I"], label=f"I{bit}", color=color, lw=1.1)
        axes[2].plot(t_h, df[f"bit{bit}_C"], label=f"C{bit}", color=color, lw=0.9, ls="--")
    axes[2].set_ylabel("I and I-RDF (uM)")
    axes[2].legend(ncol=6, loc="upper right")

    axes[3].step(t_h, df["code"], where="post", color="#6A3D9A", lw=1.3, label="decoded biochemical code")
    axes[3].scatter(samples["sample_time_min"] / 60.0, samples["code"], color="#6A3D9A", s=24, zorder=3)
    axes[3].set_ylabel("Code")
    axes[3].set_xlabel("Time (h)")
    axes[3].set_yticks(range(8))
    axes[3].legend(loc="upper right")
    fig.suptitle("Full 33-state biochemical ODE with continuous carry extension")
    fig.tight_layout()
    fig.savefig(out_dir / "01_full_3bit_overview.png", dpi=200)
    plt.close(fig)

    fig, axes = plt.subplots(3, 1, figsize=(14, 9), sharex=True)
    for bit, color in enumerate(colors):
        ax = axes[bit]
        ax.plot(t_h, df[f"bit{bit}_M_I"], label="M_I", lw=0.9)
        ax.plot(t_h, df[f"bit{bit}_I_u"], label="I_u", lw=0.9)
        ax.plot(t_h, df[f"bit{bit}_I"], label="I", lw=1.1)
        ax.plot(t_h, df[f"bit{bit}_M_T"], label="M_T", lw=0.9)
        ax.plot(t_h, df[f"bit{bit}_T"], label="T", lw=1.0)
        ax.plot(t_h, df[f"bit{bit}_M_R"], label="M_R", lw=0.9)
        ax.plot(t_h, df[f"bit{bit}_R"], label="R", lw=1.0)
        ax.plot(t_h, df[f"bit{bit}_C"], label="C", lw=1.1)
        ax.set_ylabel(f"bit {bit} (uM)")
        ax.legend(ncol=8, fontsize=8, loc="upper right")
    axes[-1].set_xlabel("Time (h)")
    fig.suptitle("Explicit expression, maturation and complex dynamics")
    fig.tight_layout()
    fig.savefig(out_dir / "02_all_biochemical_states.png", dpi=200)
    plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(14, 7), sharex=True)
    axes[0].plot(t_h, df["E0"], label="E0 memory protein", color="#0072B2")
    axes[0].plot(t_h, df["E1"], label="E1 memory protein", color="#D55E00")
    axes[0].set_ylabel("Edge-memory protein (uM)")
    axes[0].legend()
    axes[1].plot(t_h, df["U1"], label="carry 0 to 1", color="#0072B2")
    axes[1].plot(t_h, df["U2"], label="carry 1 to 2", color="#D55E00")
    axes[1].set_ylabel("Carry promoter activity")
    axes[1].set_xlabel("Time (h)")
    axes[1].legend()
    fig.suptitle("Continuous prospective biochemical falling-edge detector")
    fig.tight_layout()
    fig.savefig(out_dir / "03_continuous_carry_module.png", dpi=200)
    plt.close(fig)


def write_report(out_dir: Path, summary: dict, samples: pd.DataFrame, p: Parameters):
    codes = ", ".join(map(str, summary["sampled_codes"]))
    expected = ", ".join(map(str, summary["ideal_binary_codes"]))
    ranges = summary["state_ranges"]
    report = f"""# 完整 11 状态三级生化 ODE 重建结果

## 模型范围

- 每个正交重组模块严格保留 11 个状态：`M_I, I_u, I, M_T, T_u, T, M_R, R_u, R, C, S`。
- 三级基础模型共 33 个状态。
- 为连接三级，另加入两个连续的 carry-memory 蛋白 `E0, E1`，总计 35 个 ODE 状态。
- 整个仿真不使用事件回调、逻辑位、队列、强制翻转或 `dS/dt=0` 锁存。

## 与 `初步.pdf` 的关系

单 bit 的 11 条方程逐项实现 PDF 的转录、翻译、成熟、降解、Int-RDF 显式复合以及 DNA 正反重组。PDF 没有给出可直接实现的三级进位反应，因此 carry 扩展被单独标记为 prospective hypothesis：LR 状态产生一个有寿命的记忆蛋白 E；当 E 尚高而 S 已回到 PB 时，AND 门产生有限宽度输入给下一位。

## 当前基准结果

- 采样得到的生化状态码：`[{codes}]`
- 理想 3-bit 序列：`[{expected}]`
- 逐拍一致率：{summary['binary_sequence_accuracy']:.3f}
- S0 范围：{ranges['S0'][0]:.4f} 到 {ranges['S0'][1]:.4f}
- S1 范围：{ranges['S1'][0]:.4f} 到 {ranges['S1'][1]:.4f}
- S2 范围：{ranges['S2'][0]:.4f} 到 {ranges['S2'][1]:.4f}

三级输入峰的中位间隔分别为 U0={summary['carry_pulse_periods']['U0']['median_interval_min']:.1f} min、U1={summary['carry_pulse_periods']['U1']['median_interval_min']:.1f} min、U2={summary['carry_pulse_periods']['U2']['median_interval_min']:.1f} min。三者基本同频，说明当前连续 falling-edge carry 只是逐级延迟传递，并未产生二进制计数需要的二分频。

所有保存值均为有限数：{summary['all_values_finite']}；所有动态状态的最小数值为 {summary['minimum_dynamic_value']:.3e}。

## 解释边界

这份结果首先回答“完整生化 ODE 自己实际做了什么”，不以得到漂亮的 000 到 111 为先验目标。若状态码未形成理想二进制循环，不能用外部数字逻辑修补后再宣称三级 ODE 已经跑通。后续参数扫描和结构修改必须继续保持 35 个状态的连续性，并明确区分 PDF 原方程与新增 carry 假设。

## 主要参数

- 细胞倍增时间：{p.doubling_time:.1f} min
- 外部时钟周期：{p.clock_period:.1f} min
- 时钟脉冲 sigma：{p.clock_sigma:.1f} min
- Int-RDF 结合：k_on={p.k_on_C:.3g} uM^-1 min^-1，k_off={p.k_off_C:.3g} min^-1
- 正向/逆向最大重组速率：{p.k_rec_fwd:.3g}/{p.k_rec_rev:.3g} min^-1

## 文件

- `full_11state_3bit_ode.py`：模型与分析脚本。
- `parameters.json`：完整参数。
- `trajectories.csv`：35 个状态、三级输入和解码结果。
- `clock_samples.csv`：逐拍采样。
- `pulse_summary.csv`：连续输入峰值。
- `summary.json`：机器可读结论。
- `01_full_3bit_overview.png`：总体轨迹。
- `02_all_biochemical_states.png`：各 bit 的显式生化状态。
- `03_continuous_carry_module.png`：连续 carry 模块。
"""
    (out_dir / "重建结果说明.md").write_text(report, encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parent / "results")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    p = Parameters()
    sol = simulate(p)
    df = trajectories(sol, p)
    samples, pulses, summary = analyse(df, p)

    df.to_csv(args.output / "trajectories.csv", index=False)
    samples.to_csv(args.output / "clock_samples.csv", index=False)
    pulses.to_csv(args.output / "pulse_summary.csv", index=False)
    (args.output / "parameters.json").write_text(
        json.dumps({**asdict(p), "mu": p.mu, "t_end": p.t_end}, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    (args.output / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    make_figures(df, samples, args.output)
    write_report(args.output, summary, samples, p)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
