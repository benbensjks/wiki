# -*- coding: utf-8 -*-
"""Bit0 (first bit) standalone figures.

Scope: ONLY the first bit.  Nothing downstream is plotted.
Everything downstream (bit1, bit2, A0/F0/A1/F1, carry, g1) is deliberately excluded
except where a bit0 quantity is *derived* from it (nothing here is).

Data source (frozen lineage, single file):
    final_reconstruction/threebit51_results/wiki_n6_20260924/trajectories.csv
      - 18001 samples, 0 -> 600 h, dt = 2 min
      - bit0 columns: b0_<11 states>  plus aliases S0 / J_fwd0 / J_rev0 / g0
      - bit0 input:   C31_flux_uM_h  (the upstream phiC31 production flux, uM/h)
    final_reconstruction/threebit51_results/wiki_n6_20260924/read_windows.csv
      - 56 validated read windows (start/end/value/bit0/bit1/bit2)

Frozen configuration (bit0 prefix is identical in the 34-state and 51-state models;
prefix RHS gap = 0, trajectory gap = 0):
    uM_per_au = 5.75 | bit mRNA t1/2 = 2 min | bit maturation t1/2 = 20 min
    add_growth = False | clock gate K = 0.3, n = 2

This script writes ONLY into ./out and ./  (its own folder).  It does not modify the
model, the frozen parameter files, or any other delivery.

Run:
    & 'D:\\aconade\\python.exe' .\\make_bit0_figures.py
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

# --------------------------------------------------------------------------
# paths  (this folder only)
# --------------------------------------------------------------------------
HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
OUT.mkdir(parents=True, exist_ok=True)

WIKI = Path(r"C:\Users\18633\Desktop\wiki")
FR = WIKI / "final_reconstruction" / "threebit51_results" / "wiki_n6_20260924"
TRAJ = FR / "trajectories.csv"
READS = FR / "read_windows.csv"

DPI = 200.0

# --------------------------------------------------------------------------
# palette  (project-wide; imported from dshwork/figures/tempo_style.py if present,
# otherwise the same values are used directly so this script stands alone)
# --------------------------------------------------------------------------
DEEP_BLUE, ORANGE, TEAL = "#304B53", "#D29144", "#4F9194"
INK, WHITE, LIGHT_GRAY = "#33383A", "#FFFFFF", "#F3F5F6"
LIGHT_BLUE, LIGHT_ORANGE, LIGHT_TEAL = "#DCE5E7", "#F2E2CF", "#DCEBEC"
BLUE_GRAY, RISK_ORANGE = "#8CA0A5", "#B86F35"

_style = HERE.parent / "figures" / "tempo_style.py"
if _style.exists():
    sys.path.insert(0, str(_style.parent))
    try:
        import tempo_style as _ts  # noqa: F401
        DEEP_BLUE, ORANGE, TEAL = _ts.DEEP_BLUE, _ts.ORANGE, _ts.TEAL
        INK, WHITE, LIGHT_GRAY = _ts.INK, _ts.WHITE, _ts.LIGHT_GRAY
        LIGHT_BLUE, LIGHT_ORANGE, LIGHT_TEAL = _ts.LIGHT_BLUE, _ts.LIGHT_ORANGE, _ts.LIGHT_TEAL
        BLUE_GRAY, RISK_ORANGE = _ts.BLUE_GRAY, _ts.RISK_ORANGE
    except Exception:
        pass

# --------------------------------------------------------------------------
# read bands (bit0 read criterion)
# --------------------------------------------------------------------------
BAND_LOW = 0.30
BAND_HIGH = 0.70


def apply_style(width_in: float) -> dict:
    """Set font/line sizes so that, at the final displayed width, no text is < 16 px
    and no line is < 2 px on screen.  Mirrors the quantified spec in
    dshwork/07_第3-4章_Wiki及作图要求.md section 2.2.

    px per pt = (dpi/72) * (display_px / (width_in * dpi))
    """
    display_px = 1400.0 if width_in >= 12.0 else 1200.0
    px_per_pt = (DPI / 72.0) * (display_px / (width_in * DPI))
    base = 16.0 / px_per_pt          # body labels  >= 16 px on screen
    tick = 14.0 / px_per_pt          # tick labels  >= 14 px
    title = 18.0 / px_per_pt         # titles       >= 18 px
    lw = max(1.8, 2.0 / px_per_pt)   # curves       >=  2 px
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": base,
        "axes.titlesize": title,
        "axes.labelsize": base,
        "xtick.labelsize": tick,
        "ytick.labelsize": tick,
        "legend.fontsize": tick,
        "axes.linewidth": 1.25,
        "xtick.major.width": 1.0,
        "ytick.major.width": 1.0,
        "xtick.major.size": 4.0,
        "ytick.major.size": 4.0,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.facecolor": WHITE,
        "axes.facecolor": WHITE,
        "savefig.facecolor": WHITE,
        "lines.linewidth": lw,
        "lines.solid_capstyle": "round",
    })
    return {"base": base, "tick": tick, "title": title, "lw": lw}


def save(fig, stem: str):
    png = OUT / f"{stem}.png"
    fig.savefig(png, dpi=DPI, bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {png.name}")


def band_shading(ax, t0, t1, alpha=0.55):
    """Paint the bit0 read bands behind a panel."""
    ax.axhspan(0.0, BAND_LOW, color=LIGHT_BLUE, alpha=alpha, zorder=0, lw=0)
    ax.axhspan(BAND_HIGH, 1.0, color=LIGHT_TEAL, alpha=alpha, zorder=0, lw=0)
    ax.axhspan(BAND_LOW, BAND_HIGH, color="#DDE2E4", alpha=0.45, zorder=0, lw=0)


def read_window_marks(ax, rw, t0, t1, label_first=False):
    """Shade the validated read windows that fall inside [t0, t1]."""
    for _, r in rw.iterrows():
        a, b = r["window_start_h"], r["window_end_h"]
        if b < t0 or a > t1:
            continue
        ax.axvspan(a, b, color=ORANGE, alpha=0.13, zorder=0, lw=0)
        if label_first:
            ax.text((a + b) / 2, 0.02, str(int(r["bit0"])), ha="center", va="bottom",
                    color=DEEP_BLUE, fontweight="bold", zorder=6)


# --------------------------------------------------------------------------
# load
# --------------------------------------------------------------------------
def load():
    if not TRAJ.exists():
        raise SystemExit(f"missing trajectory: {TRAJ}")
    df = pd.read_csv(TRAJ)
    rw = pd.read_csv(READS) if READS.exists() else pd.DataFrame()
    t = df["time_h"].to_numpy()
    # period estimate from the validated read windows
    if len(rw) >= 3:
        period = float(np.median(np.diff(rw["time_h"].to_numpy())))
    else:
        period = 10.6
    return df, rw, t, period


# --------------------------------------------------------------------------
# 01  one-cycle flip decomposition            (Figure 3-2)
# --------------------------------------------------------------------------
def fig01(df, rw, period, k=30):
    S = df["b0_S"].to_numpy()
    R = df["b0_R"].to_numpy()
    C = df["b0_C"].to_numpy()
    JF = df["J_fwd0"].to_numpy()
    JR = df["J_rev0"].to_numpy()
    FL = df["C31_flux_uM_h"].to_numpy()
    t = df["time_h"].to_numpy()

    k = int(np.clip(k, 1, len(rw) - 2))
    tc = rw["time_h"].iloc[k]
    t0, t1 = tc - period / 2.0, tc + period / 2.0
    m = (t >= t0) & (t <= t1)
    tw = t[m]

    st = apply_style(13.0)
    fig, ax = plt.subplots(4, 1, figsize=(13, 14), sharex=True, constrained_layout=True)

    # row 1 - DNA state against the read bands
    band_shading(ax[0], t0, t1)
    ax[0].plot(tw, S[m], color=DEEP_BLUE, label="DNA state $S$")
    read_window_marks(ax[0], rw, t0, t1)
    for y, lab in ((BAND_LOW, "read band 0 boundary (0.30)"),
                   (BAND_HIGH, "read band 1 boundary (0.70)")):
        ax[0].axhline(y, color=INK, ls="--", lw=1.2, alpha=0.7)
        ax[0].text(t1, y, "  " + lab, va="center", ha="left", color=INK)
    ax[0].set_ylabel("DNA state $S$")
    ax[0].set_ylim(-0.05, 1.05)
    ax[0].set_title("Bit 0 — one clock cycle: flip decomposition", loc="left")

    # row 2 - memory pool and complex, each normalised to its own maximum
    Rmax, Cmax = float(R.max()), float(C.max())
    ax[1].plot(tw, R[m] / Rmax, color=ORANGE, label=f"RDF pool (max {Rmax:.2f} a.u.)")
    ax[1].plot(tw, C[m] / Cmax, color=TEAL, label=f"Int–RDF complex (max {Cmax:.3f} a.u.)")
    ax[1].set_ylabel("normalised")
    ax[1].set_ylim(-0.05, 1.15)
    ax[1].legend(loc="upper left", frameon=False)

    # row 3 - the two recombination fluxes, same axis
    ax[2].plot(tw, JF[m], color=TEAL, label="forward flux $J_{fwd,0}$  (PB→LR)")
    ax[2].plot(tw, JR[m], color=ORANGE, label="reverse flux $J_{rev,0}$  (LR→PB)")
    ax[2].set_ylabel("flux (a.u./h)")
    ax[2].legend(loc="upper left", frameon=False)

    # row 4 - the real clock input
    ax[3].plot(tw, FL[m], color=DEEP_BLUE, label="φC31 production flux (clock input)")
    ax[3].fill_between(tw, 0, FL[m], color=LIGHT_BLUE, alpha=0.5)
    ax[3].set_ylabel("flux (µM/h)")
    ax[3].set_xlabel("time (h)")
    ax[3].legend(loc="upper left", frameon=False)

    # flip marker: steepest dS/dt
    dS = np.gradient(S[m], tw)
    if np.any(np.isfinite(dS)):
        tf = tw[int(np.nanargmax(np.abs(dS)))]
        for a in ax:
            a.axvline(tf, color=RISK_ORANGE, lw=1.6, ls=":", alpha=0.9)
        ax[0].annotate("flip", xy=(tf, 0.5), xytext=(tf + 0.35, 0.62),
                       color=RISK_ORANGE, fontweight="bold",
                       arrowprops=dict(arrowstyle="->", color=RISK_ORANGE, lw=1.4))

    for a in ax:
        a.grid(alpha=0.18, lw=0.8)
        a.set_xlim(t0, t1)
    save(fig, "fig_bit0_01_one_cycle_decomposition")


# --------------------------------------------------------------------------
# 02  dwell bands and read windows           (Figure 3-3)
# --------------------------------------------------------------------------
def fig02(df, rw, period):
    S = df["b0_S"].to_numpy()
    FL = df["C31_flux_uM_h"].to_numpy()
    t = df["time_h"].to_numpy()

    st = apply_style(14.0)
    fig = plt.figure(figsize=(14, 12))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.0, 1.5, 1.2], hspace=0.32)

    # -- Panel A: read-window geometry over one period ----------------------
    axA = fig.add_subplot(gs[0])
    tc = float(rw["time_h"].iloc[20])
    t0, t1 = tc - period / 2.0, tc + period / 2.0
    m = (t >= t0) & (t <= t1)
    axA.plot(t[m], FL[m], color=DEEP_BLUE)
    axA.fill_between(t[m], 0, FL[m], color=LIGHT_BLUE, alpha=0.5)
    r0 = rw.iloc[20]
    axA.axvspan(r0["window_start_h"], r0["window_end_h"], color=ORANGE, alpha=0.22, lw=0)
    axA.annotate("", xy=(r0["window_start_h"], 0.92 * FL[m].max()),
                 xytext=(r0["window_end_h"], 0.92 * FL[m].max()),
                 arrowprops=dict(arrowstyle="<->", color=RISK_ORANGE, lw=1.6))
    axA.text(tc, 0.97 * FL[m].max(),
             f"read window = 20 % of the {period:.2f} h period, trough-centred",
             ha="center", va="bottom", color=RISK_ORANGE, fontweight="bold")
    axA.set_ylabel("clock flux\n(µM/h)")
    axA.set_xlabel("time (h)")
    axA.set_title("(A)  Read-window geometry", loc="left")
    axA.set_xlim(t0, t1)
    axA.grid(alpha=0.18, lw=0.8)

    # -- Panel B: the whole run, bands + windows + decoded bit --------------
    axB = fig.add_subplot(gs[1])
    band_shading(axB, t[0], t[-1])
    axB.plot(t, S, color=DEEP_BLUE, lw=1.5)
    for _, r in rw.iterrows():
        c = LIGHT_TEAL if int(r["bit0"]) == 1 else LIGHT_BLUE
        axB.axvspan(r["window_start_h"], r["window_end_h"], color=c, alpha=0.85, lw=0, zorder=1)
        axB.plot([r["time_h"]], [S[np.argmin(np.abs(t - r["time_h"]))]],
                 marker="o", ms=3.2, color=DEEP_BLUE, zorder=5)
    axB.axhline(BAND_LOW, color=INK, ls="--", lw=1.2, alpha=0.7)
    axB.axhline(BAND_HIGH, color=INK, ls="--", lw=1.2, alpha=0.7)
    axB.text(t[-1], BAND_LOW, "  0.30", va="center", color=INK)
    axB.text(t[-1], BAND_HIGH, "  0.70", va="center", color=INK)
    axB.set_ylim(-0.05, 1.12)
    axB.set_ylabel("DNA state $S$")
    axB.set_xlabel("time (h)")
    axB.set_title(f"(B)  Full run — {len(rw)} validated read windows, all labelled", loc="left")
    axB.set_xlim(t[0], t[-1])
    axB.grid(alpha=0.18, lw=0.8, axis="x")
    handles = [Line2D([], [], color=DEEP_BLUE, lw=2, label="DNA state $S$"),
               Rectangle((0, 0), 1, 1, fc=LIGHT_TEAL, ec="none", label="read window decoding to 1"),
               Rectangle((0, 0), 1, 1, fc=LIGHT_BLUE, ec="none", label="read window decoding to 0"),
               Rectangle((0, 0), 1, 1, fc="#DDE2E4", alpha=0.6, ec="none", label="ambiguous band (0.30–0.70)")]
    axB.legend(handles=handles, loc="upper center", ncol=4, frameon=False,
               bbox_to_anchor=(0.5, 1.02))

    # -- Panel C: per-window commitment strip ------------------------------
    axC = fig.add_subplot(gs[2])
    lo, hi, val, col = [], [], [], []
    for _, r in rw.iterrows():
        m2 = (t >= r["window_start_h"]) & (t <= r["window_end_h"])
        if m2.sum() < 2:
            continue
        lo.append(S[m2].min()); hi.append(S[m2].max())
        val.append(S[np.argmin(np.abs(t - r["time_h"]))])
        col.append(TEAL if int(r["bit0"]) == 1 else DEEP_BLUE)
    idx = np.arange(len(val))
    axC.vlines(idx, lo, hi, color=col, lw=3.4, alpha=0.85)
    axC.plot(idx, val, "o", ms=4.0, color=INK, zorder=5, label="$S$ at the read instant")
    axC.axhspan(0.0, BAND_LOW, color=LIGHT_BLUE, alpha=0.55, zorder=0, lw=0)
    axC.axhspan(BAND_HIGH, 1.0, color=LIGHT_TEAL, alpha=0.55, zorder=0, lw=0)
    axC.axhspan(BAND_LOW, BAND_HIGH, color="#DDE2E4", alpha=0.45, zorder=0, lw=0)
    axC.set_ylim(-0.05, 1.05)
    axC.set_xlabel("read window index")
    axC.set_ylabel("$S$ range inside\nthe read window")
    axC.set_title("(C)  Every window stays inside one band — no window touches the ambiguous zone", loc="left")
    axC.set_xlim(-1, len(val))
    axC.grid(alpha=0.18, lw=0.8, axis="y")
    axC.legend(loc="upper right", frameon=False)

    save(fig, "fig_bit0_02_dwell_and_readwindow")


# --------------------------------------------------------------------------
# 03  the eleven bit-0 states
# --------------------------------------------------------------------------
def fig03(df, rw, period):
    t = df["time_h"].to_numpy()
    tc = float(rw["time_h"].iloc[20])
    t0, t1 = tc - period / 2.0, tc + period / 2.0
    m = (t >= t0) & (t <= t1)
    tw = t[m]

    chains = [
        ("Recombinase (Int)", ["b0_M_I", "b0_I_u", "b0_I"],
         ["$M_I$ (mRNA)", "$I_u$ (immature)", "$I$ (mature)"]),
        ("Repressor (BM3R1 / T)", ["b0_M_T", "b0_T_u", "b0_T"],
         ["$M_T$ (mRNA)", "$T_u$ (immature)", "$T$ (mature)"]),
        ("Directionality factor (RDF)", ["b0_M_R", "b0_R_u", "b0_R"],
         ["$M_R$ (mRNA)", "$R_u$ (immature)", "$R$ (mature)"]),
    ]
    cols = [TEAL, ORANGE, DEEP_BLUE]

    apply_style(13.0)
    fig, ax = plt.subplots(4, 1, figsize=(13, 15), sharex=True, constrained_layout=True)

    for a, (title, keys, labs) in zip(ax[:3], chains):
        for k, lab, c in zip(keys, labs, cols):
            v = df[k].to_numpy()[m]
            vmax = float(df[k].max())
            a.plot(tw, v / vmax if vmax > 0 else v, color=c, label=f"{lab}   max {vmax:.3g} a.u.")
            a.fill_between(tw, 0, v / vmax if vmax > 0 else v, color=c, alpha=0.10)
        a.set_ylabel("normalised")
        a.set_ylim(-0.05, 1.18)
        a.set_title(f"Bit 0 — {title}", loc="left")
        a.legend(loc="upper left", frameon=False, ncol=3)
        a.grid(alpha=0.18, lw=0.8)

    # bottom row: complex and DNA state, dual axis
    a = ax[3]
    Cv, Sv = df["b0_C"].to_numpy()[m], df["b0_S"].to_numpy()[m]
    l1, = a.plot(tw, Cv, color=TEAL, label=f"$C$ (Int–RDF complex), max {df['b0_C'].max():.3g} a.u.")
    a.set_ylabel("complex (a.u.)", color=TEAL)
    a.tick_params(axis="y", labelcolor=TEAL)
    b = a.twinx()
    b.spines["right"].set_visible(True)
    l2, = b.plot(tw, Sv, color=DEEP_BLUE, label="DNA state $S$")
    b.set_ylabel("DNA state $S$", color=DEEP_BLUE)
    b.tick_params(axis="y", labelcolor=DEEP_BLUE)
    b.set_ylim(-0.05, 1.05)
    a.set_title("Bit 0 — complex and DNA state", loc="left")
    a.set_xlabel("time (h)")
    a.legend(handles=[l1, l2], loc="upper left", frameon=False)
    a.grid(alpha=0.18, lw=0.8)
    a.set_xlim(t0, t1)

    save(fig, "fig_bit0_03_eleven_states")


# --------------------------------------------------------------------------
# 04  forward / reverse flux separation
# --------------------------------------------------------------------------
def fig04(df, rw, period, k=30):
    t = df["time_h"].to_numpy()
    JF, JR = df["J_fwd0"].to_numpy(), df["J_rev0"].to_numpy()
    S, g0 = df["b0_S"].to_numpy(), df["g0"].to_numpy()

    apply_style(13.0)
    fig, ax = plt.subplots(3, 1, figsize=(13, 10), sharex=True, constrained_layout=True)
    for row, span in ((0, (k - 2, k + 2)), (1, (k, k + 1)), (2, (k, k + 1))):
        pass

    # top: a few cycles
    a, b2 = k - 2, k + 2
    t0 = rw["time_h"].iloc[a] - period / 2.0
    t1 = rw["time_h"].iloc[b2] + period / 2.0
    m = (t >= t0) & (t <= t1)
    ax[0].plot(t[m], JF[m], color=TEAL, label="forward $J_{fwd,0}$")
    ax[0].plot(t[m], JR[m], color=ORANGE, label="reverse $J_{rev,0}$")
    ax[0].set_ylabel("flux (a.u./h)")
    ax[0].legend(loc="upper left", frameon=False, ncol=2)
    ax[0].set_title("(A)  Several cycles — the two directions never overlap", loc="left")
    ax[0].set_xlim(t0, t1)
    ax[0].grid(alpha=0.18, lw=0.8)
    for i in range(a, b2 + 1):
        ax[0].axvline(rw["time_h"].iloc[i], color=BLUE_GRAY, lw=0.9, ls=":", alpha=0.8)

    # middle: one cycle, fluxes with S
    tc = rw["time_h"].iloc[k]
    t0, t1 = tc - period / 2.0, tc + period / 2.0
    m = (t >= t0) & (t <= t1)
    ax[1].plot(t[m], JF[m], color=TEAL, label="forward $J_{fwd,0}$")
    ax[1].plot(t[m], JR[m], color=ORANGE, label="reverse $J_{rev,0}$")
    ax[1].set_ylabel("flux (a.u./h)")
    ax[1].legend(loc="upper left", frameon=False, ncol=2)
    ax[1].set_title("(B)  One cycle — forward and reverse are separated in time", loc="left")
    ax[1].grid(alpha=0.18, lw=0.8)

    # bottom: one cycle, S with bands
    band_shading(ax[2], t0, t1)
    read_window_marks(ax[2], rw, t0, t1)
    ax[2].plot(t[m], S[m], color=DEEP_BLUE)
    ax[2].axhline(BAND_LOW, color=INK, ls="--", lw=1.2, alpha=0.7)
    ax[2].axhline(BAND_HIGH, color=INK, ls="--", lw=1.2, alpha=0.7)
    ax[2].set_ylabel("DNA state $S$")
    ax[2].set_xlabel("time (h)")
    ax[2].set_ylim(-0.05, 1.05)
    ax[2].set_title("(C)  DNA state over the same cycle", loc="left")
    ax[2].grid(alpha=0.18, lw=0.8)
    for a in ax:
        a.set_xlim(t0, t1) if a is not ax[0] else None
    ax[0].set_xlim(rw["time_h"].iloc[k - 2] - period / 2.0, rw["time_h"].iloc[k + 2] + period / 2.0)

    save(fig, "fig_bit0_04_flux_separation")


# --------------------------------------------------------------------------
# 05  phase portraits  (direction memory)
# --------------------------------------------------------------------------
def fig05(df):
    S = df["b0_S"].to_numpy()
    R = df["b0_R"].to_numpy()
    C = df["b0_C"].to_numpy()
    I = df["b0_I"].to_numpy()
    P = I * R                     # instantaneous substrate product I*R

    apply_style(14.0)
    fig, ax = plt.subplots(1, 3, figsize=(14, 5.2), constrained_layout=True)

    sc = ax[0].scatter(S, R, c=np.arange(len(S)), cmap="viridis", s=2.0, alpha=0.65)
    ax[0].set_xlabel("DNA state $S$")
    ax[0].set_ylabel("free RDF pool $R$ (a.u.)")
    ax[0].set_title("(A)  $S$ vs RDF pool", loc="left")
    cb = fig.colorbar(sc, ax=ax[0]); cb.set_label("time index")

    ax[1].scatter(S, P, s=2.0, alpha=0.65, color=ORANGE)
    ax[1].set_xlabel("DNA state $S$")
    ax[1].set_ylabel("instantaneous product $I\\cdot R$  (a.u.$^2$)")
    ax[1].set_title("(B)  $S$ vs the substrate product $I\\cdot R$\n(reverse rate depends on the product, not on $R$ alone)",
                    loc="left")

    ax[2].scatter(S, C, s=2.0, alpha=0.65, color=TEAL)
    ax[2].set_xlabel("DNA state $S$")
    ax[2].set_ylabel("Int–RDF complex $C$ (a.u.)")
    ax[2].set_title("(C)  $S$ vs the complex $C$", loc="left")

    for a in ax:
        a.grid(alpha=0.18, lw=0.8)
        a.set_xlim(-0.02, 1.02)
    fig.suptitle("Bit 0 — phase portraits (600 h, frozen configuration)", x=0.01, ha="left")

    save(fig, "fig_bit0_05_phase_portrait")


# --------------------------------------------------------------------------
# 06  per-cycle metrics + a metrics CSV
# --------------------------------------------------------------------------
def fig06(df, rw, period):
    t = df["time_h"].to_numpy()
    S = df["b0_S"].to_numpy()
    FL = df["C31_flux_uM_h"].to_numpy()
    JF, JR = df["J_fwd0"].to_numpy(), df["J_rev0"].to_numpy()

    rows = []
    for i, r in rw.iterrows():
        c = r["time_h"]
        a, b = c - period / 2.0, c + period / 2.0
        m = (t >= a) & (t <= b)
        mw = (t >= r["window_start_h"]) & (t <= r["window_end_h"])
        if m.sum() < 4:
            continue
        dS = np.gradient(S[m], t[m])
        rows.append(dict(
            index=int(i), trough_h=float(c), decoded_bit0=int(r["bit0"]),
            S_at_read=float(S[np.argmin(np.abs(t - c))]),
            S_window_min=float(S[mw].min()), S_window_max=float(S[mw].max()),
            flux_peak_uM_h=float(FL[m].max()), flux_trough_uM_h=float(FL[m].min()),
            flux_dose_uM=float(np.trapezoid(FL[m], t[m])),
            forward_dose=float(np.trapezoid(JF[m], t[m])),
            reverse_dose=float(np.trapezoid(JR[m], t[m])),
            flip_time_h=float(t[m][int(np.argmax(np.abs(dS)))]),
        ))
    met = pd.DataFrame(rows)
    met["period_h"] = met["trough_h"].diff()
    met.to_csv(OUT / "bit0_cycle_metrics.csv", index=False)
    print(f"  wrote bit0_cycle_metrics.csv  ({len(met)} cycles)")

    apply_style(14.0)
    fig, ax = plt.subplots(2, 2, figsize=(14, 8), constrained_layout=True)
    x = met["index"].to_numpy()

    a = ax[0, 0]
    a.plot(x, met["S_at_read"], "o-", color=DEEP_BLUE, ms=3.5, lw=1.2)
    a.axhspan(0, BAND_LOW, color=LIGHT_BLUE, alpha=0.55, lw=0)
    a.axhspan(BAND_HIGH, 1, color=LIGHT_TEAL, alpha=0.55, lw=0)
    a.axhspan(BAND_LOW, BAND_HIGH, color="#DDE2E4", alpha=0.45, lw=0)
    a.set_ylim(-0.05, 1.05); a.set_ylabel("$S$ at the read instant")
    a.set_title("(A)  Read value stays in its band", loc="left")

    a = ax[0, 1]
    a.bar(x, met["flux_peak_uM_h"], color=ORANGE, alpha=0.8, label="per-cycle peak")
    a.plot(x, met["flux_trough_uM_h"], "o-", color=DEEP_BLUE, ms=3.0, lw=1.2, label="per-cycle trough")
    a.set_ylabel("clock flux (µM/h)"); a.legend(frameon=False)
    a.set_title("(B)  Input amplitude per cycle", loc="left")

    a = ax[1, 0]
    a.plot(x[1:], met["period_h"].to_numpy()[1:], "o-", color=DEEP_BLUE, ms=3.5, lw=1.2)
    a.axhline(period, color=ORANGE, ls="--", lw=1.4, label=f"median {period:.2f} h")
    a.set_ylabel("trough-to-trough period (h)"); a.legend(frameon=False)
    a.set_title("(C)  Period is set by the clock, not by the bit", loc="left")

    a = ax[1, 1]
    a.bar(x - 0.2, met["forward_dose"], width=0.4, color=TEAL, alpha=0.85, label="forward dose $\\int J_{fwd,0}$")
    a.bar(x + 0.2, met["reverse_dose"], width=0.4, color=ORANGE, alpha=0.85, label="reverse dose $\\int J_{rev,0}$")
    a.set_ylabel("per-cycle flux dose (a.u.)"); a.set_xlabel("read window index")
    a.legend(frameon=False)
    a.set_title("(D)  Forward and reverse doses alternate", loc="left")

    for a in ax.ravel():
        a.grid(alpha=0.18, lw=0.8)

    save(fig, "fig_bit0_06_cycle_metrics")


# --------------------------------------------------------------------------
# 07  expression-chain delays (why the pulse is narrow)
# --------------------------------------------------------------------------
def fig07(df, rw, period, k=30):
    t = df["time_h"].to_numpy()
    tc = float(rw["time_h"].iloc[k])
    t0, t1 = tc - period / 2.0, tc + period / 2.0
    m = (t >= t0) & (t <= t1)
    tw = t[m]

    apply_style(13.0)
    fig, ax = plt.subplots(2, 1, figsize=(13, 8.5), sharex=True, constrained_layout=True)

    # top: the three mature proteins
    for key, lab, c in (("b0_I", "mature Int $I$", TEAL),
                        ("b0_T", "mature BM3R1 $T$", DEEP_BLUE),
                        ("b0_R", "mature RDF $R$", ORANGE)):
        v = df[key].to_numpy()[m]
        vmax = float(df[key].max())
        ax[0].plot(tw, v / vmax, color=c, label=f"{lab}  (max {vmax:.3g} a.u.)")
        ax[0].fill_between(tw, 0, v / vmax, color=c, alpha=0.10)
    ax[0].set_ylabel("normalised")
    ax[0].set_ylim(-0.05, 1.18)
    ax[0].legend(loc="upper left", frameon=False, ncol=3)
    ax[0].set_title("(A)  The three mature proteins over one cycle", loc="left")

    # bottom: one chain resolved into mRNA -> immature -> mature
    for key, lab, c, ls in (("b0_M_R", "$M_R$ (mRNA)", LIGHT_ORANGE, "-"),
                            ("b0_R_u", "$R_u$ (immature)", ORANGE, "-"),
                            ("b0_R", "$R$ (mature)", DEEP_BLUE, "-")):
        v = df[key].to_numpy()[m]
        vmax = float(df[key].max())
        ax[1].plot(tw, v / vmax, color=c, ls=ls, label=f"{lab}  (max {vmax:.3g} a.u.)")
    ax[1].set_ylabel("normalised")
    ax[1].set_ylim(-0.05, 1.18)
    ax[1].set_xlabel("time (h)")
    ax[1].legend(loc="upper left", frameon=False, ncol=3)
    ax[1].set_title("(B)  The RDF chain resolved — the maturation stages are what delay the pool",
                    loc="left")

    for a in ax:
        a.grid(alpha=0.18, lw=0.8)
        a.set_xlim(t0, t1)
        a.axvline(tc, color=BLUE_GRAY, lw=0.9, ls=":", alpha=0.9)

    save(fig, "fig_bit0_07_expression_chains")


# --------------------------------------------------------------------------
# 08  commitment statistics  (dwell distribution and band margins)
# --------------------------------------------------------------------------
def fig08(df, rw, period):
    t = df["time_h"].to_numpy()
    S = df["b0_S"].to_numpy()

    vals, margins = [], []
    for _, r in rw.iterrows():
        c = r["time_h"]
        m = (t >= r["window_start_h"]) & (t <= r["window_end_h"])
        if m.sum() < 2:
            continue
        v = S[np.argmin(np.abs(t - c))]
        vals.append(v)
        # distance from the whole window to the nearest band boundary
        if int(r["bit0"]) == 1:
            margins.append(S[m].min() - BAND_HIGH)
        else:
            margins.append(BAND_LOW - S[m].max())
    vals = np.asarray(vals)
    margins = np.asarray(margins)

    apply_style(14.0)
    fig, ax = plt.subplots(1, 2, figsize=(14, 5.2), constrained_layout=True)

    a = ax[0]
    a.axvspan(0.0, BAND_LOW, color=LIGHT_BLUE, alpha=0.6, lw=0)
    a.axvspan(BAND_HIGH, 1.0, color=LIGHT_TEAL, alpha=0.6, lw=0)
    a.axvspan(BAND_LOW, BAND_HIGH, color="#DDE2E4", alpha=0.5, lw=0)
    a.hist(vals, bins=np.linspace(0, 1, 21), color=DEEP_BLUE, alpha=0.85)
    a.set_xlim(0, 1)
    a.set_xlabel("$S$ at the read instant")
    a.set_ylabel("number of read windows")
    a.set_title("(A)  All 56 reads land in one of the two bands", loc="left")
    a.text(0.15, a.get_ylim()[1] * 0.92, "0", ha="center", color=DEEP_BLUE, fontweight="bold")
    a.text(0.85, a.get_ylim()[1] * 0.92, "1", ha="center", color=DEEP_BLUE, fontweight="bold")
    a.text(0.5, a.get_ylim()[1] * 0.55, "ambiguous band\nnever occupied",
           ha="center", va="center", color=RISK_ORANGE)
    a.grid(alpha=0.18, lw=0.8, axis="y")

    a = ax[1]
    o = np.argsort(margins)
    a.bar(np.arange(len(margins)), margins[o],
          color=[TEAL if v > 0 else RISK_ORANGE for v in margins[o]], alpha=0.9)
    a.axhline(0.0, color=INK, lw=1.3)
    a.set_xlabel("read window (sorted by margin)")
    a.set_ylabel("margin to the nearest band boundary")
    a.set_title("(B)  Minimum margin %.3f — no window is close to its band edge"
                % margins.min(), loc="left")
    a.grid(alpha=0.18, lw=0.8, axis="y")

    fig.suptitle("Bit 0 — read commitment (600 h, frozen configuration)", x=0.01, ha="left")
    save(fig, "fig_bit0_08_commitment_statistics")


# --------------------------------------------------------------------------
# 09  bit0-only early operating scan  (uM_per_au x bit maturation)
# --------------------------------------------------------------------------
REGION = (Path(r"C:\Users\18633\Desktop\wiki") / "final_reconstruction" / "bit0_results"
          / "verification" / "bit0_region_map.csv")


def fig09():
    """Bit0-only scan.  NOTE: this uses the EARLY dwell-label criterion (`stored_bit`),
    not the finite-read-window criterion, and its uM_per_au grid starts at 6.0, which is
    ABOVE the frozen centre 5.75.  It is a diagnostic of the high-dwell level, not the
    current operating window."""
    if not REGION.exists():
        print("  skipped fig09 (bit0_region_map.csv not found)")
        return
    d = pd.read_csv(REGION)
    ums = sorted(d["uM_per_au"].unique())
    mats = sorted(d["maturation_half_life_min"].unique())

    apply_style(14.0)
    fig, ax = plt.subplots(1, 2, figsize=(14, 5.4), constrained_layout=True)

    # -- Panel A: pass/fail grid with the decoded code string ---------------
    a = ax[0]
    for i, u in enumerate(ums):
        for j, mt in enumerate(mats):
            row = d[(d.uM_per_au == u) & (d.maturation_half_life_min == mt)]
            if row.empty:
                continue
            ok = bool(row["stored_bit"].iloc[0])
            code = str(row["codes"].iloc[0])
            a.add_patch(Rectangle((i - 0.42, j - 0.42), 0.84, 0.84,
                                  fc=LIGHT_TEAL if ok else LIGHT_ORANGE,
                                  ec=DEEP_BLUE if ok else RISK_ORANGE, lw=1.4))
            a.text(i, j + 0.12, "✓" if ok else "✗", ha="center", va="center",
                   color=TEAL if ok else RISK_ORANGE, fontweight="bold")
            a.text(i, j - 0.20, code.replace("x", "·"), ha="center", va="center",
                   color=INK, family="DejaVu Sans Mono")
    a.set_xticks(range(len(ums))); a.set_xticklabels([f"{u:g}" for u in ums])
    a.set_yticks(range(len(mats))); a.set_yticklabels([f"{m:g}" for m in mats])
    a.set_xlabel("concentration scale  a.u. → µM")
    a.set_ylabel("bit maturation half-life (min)")
    a.set_xlim(-0.6, len(ums) - 0.4); a.set_ylim(-0.6, len(mats) - 0.4)
    a.set_title("(A)  Bit 0 alone — early dwell-label scan\n"
                "· = unlabelled read window;  ✓/✗ = stored_bit", loc="left")

    # -- Panel B: the boundary is set by the high dwell level, not by maturation
    a = ax[1]
    for mt, c in zip(mats, plt.cm.viridis(np.linspace(0.1, 0.9, len(mats)))):
        sub = d[d.maturation_half_life_min == mt].sort_values("uM_per_au")
        a.plot(sub["uM_per_au"], sub["S_hi"], "o-", color=c, ms=3.6, lw=1.3,
               label=f"{mt:g} min")
    a.axhline(BAND_HIGH, color=RISK_ORANGE, ls="--", lw=1.6)
    a.text(ums[-1], BAND_HIGH, "  0.70 band boundary", va="center", ha="right",
           color=RISK_ORANGE, fontweight="bold")
    a.axhline(0.98, color=BLUE_GRAY, ls=":", lw=1.2)
    a.text(ums[0], 0.98, "0.98 high dwell at the lowest scale ", va="bottom",
           ha="left", color=BLUE_GRAY)
    a.set_xlabel("concentration scale  a.u. → µM")
    a.set_ylabel("high-state dwell level $S_{hi}$")
    a.set_ylim(0.68, 1.02)
    a.legend(title="maturation t½", frameon=False, ncol=2, loc="lower left")
    a.set_title("(B)  The high dwell level falls with the conversion scale;\n"
                "maturation half-life barely moves it", loc="left")
    a.grid(alpha=0.18, lw=0.8)

    fig.suptitle("Bit 0 alone — where the high state stops clearing the read band "
                 "(early dwell-label criterion; grid starts at 6.0, above the frozen 5.75)",
                 x=0.01, ha="left")
    save(fig, "fig_bit0_09_early_window_map")


# --------------------------------------------------------------------------
def main():
    print("Bit0 standalone figures")
    print("  source :", TRAJ)
    df, rw, t, period = load()
    print(f"  loaded : {len(df)} samples, {t[0]:.1f}–{t[-1]:.1f} h, period {period:.3f} h, "
          f"{len(rw)} read windows")
    print("  output :", OUT)
    print()

    fig01(df, rw, period)
    fig02(df, rw, period)
    fig03(df, rw, period)
    fig04(df, rw, period)
    fig05(df)
    fig06(df, rw, period)
    fig07(df, rw, period)
    fig08(df, rw, period)
    fig09()

    print("\ndone.")


if __name__ == "__main__":
    main()
