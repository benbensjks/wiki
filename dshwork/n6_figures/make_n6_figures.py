# -*- coding: utf-8 -*-
"""n = 6 (frozen three-bit) parameter-analysis and robustness figures.

Four figures, all at the frozen working point
    n_A1_gate = 6 | carry1 mRNA t1/2 = 2 min | A1/F1 maturation t1/2 = 32.5 min
    uM_per_au = 5.75 | F1 production exponent = 4 (frozen table value)
    model_threebit51.py sha256 D244D3CD1914F50F994E6856084D4A012B820C13249DCBA0E963E444CE3C34CC

    01  gate-exponent selection          (why n = 6)
    02  molecular-pool perturbation     (136 runs)
    03  eight phase-consistent states   (000..111)
    04  strict tolerance                (5 points x 3 levels)

READ-ONLY on every existing artefact.  Writes only into ./out.

DATA TRAP (handled deliberately, see fig 03 panel B):
    `plausibility/eight_initial_verdict.json` per-value field `leak5pct_ratio` is the
    RETIRED single-window (aliased) metric.  Its spread across the eight initial states is
    min 0.0049 / median 0.0277 / max 34705 -- a 3.5e4 apparent difference produced by one
    and the same circuit.  The corrected paired metric lives in
    `plausibility/carry_pairing_eight_verdict.json` and has a relative spread of 0.196 %.
    Figure 03 plots BOTH so the aliasing is visible rather than hidden.

Run:
    & 'D:\\aconade\\python.exe' .\\make_n6_figures.py
"""
from __future__ import annotations

import json
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
HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
OUT.mkdir(parents=True, exist_ok=True)

WIKI = Path(r"C:\Users\18633\Desktop\wiki")
PL = WIKI / "final_reconstruction" / "plausibility"
PF = WIKI / "final_reconstruction" / "threebit51_results" / "postfreeze_peak_perturbation"

GRID_VERDICT = PL / "carry_pairing_grid_verdict.json"
EIGHT_VERDICT = PL / "carry_pairing_eight_verdict.json"      # per_point -> PAIRED metric
EIGHT_INITIAL = PL / "eight_initial_verdict.json"            # per_value -> sequences, margins, aliased leak
POOL_CSV = PL / "pool_perturbations_all_all.csv"
TOL_CSV = PL / "strict_tolerance_tol_all.csv"
SELECTED = PL / "threebit51_selected_v1.json"

DPI = 200.0
N_SELECTED = 6.0
CUTS = ("1pct", "5pct", "10pct")

# --------------------------------------------------------------------------
# palette
# --------------------------------------------------------------------------
DEEP_BLUE, ORANGE, TEAL = "#304B53", "#D29144", "#4F9194"
INK, WHITE, LIGHT_GRAY = "#33383A", "#FFFFFF", "#F3F5F6"
LIGHT_BLUE, LIGHT_ORANGE, LIGHT_TEAL = "#DCE5E7", "#F2E2CF", "#DCEBEC"
BLUE_GRAY, RISK_ORANGE = "#8CA0A5", "#B86F35"

_s = HERE.parent / "figures" / "tempo_style.py"
if _s.exists():
    sys.path.insert(0, str(_s.parent))
    try:
        import tempo_style as _ts
        DEEP_BLUE, ORANGE, TEAL = _ts.DEEP_BLUE, _ts.ORANGE, _ts.TEAL
        INK, WHITE, LIGHT_GRAY = _ts.INK, _ts.WHITE, _ts.LIGHT_GRAY
        LIGHT_BLUE, LIGHT_ORANGE, LIGHT_TEAL = _ts.LIGHT_BLUE, _ts.LIGHT_ORANGE, _ts.LIGHT_TEAL
        BLUE_GRAY, RISK_ORANGE = _ts.BLUE_GRAY, _ts.RISK_ORANGE
    except Exception:
        pass


def apply_style(width_in: float) -> dict:
    """Font/line sizes derived from the final displayed width: labels >= 16 px,
    ticks >= 14 px, titles >= 18 px, lines >= 2 px on screen."""
    display_px = 1400.0 if width_in >= 12.0 else 1200.0
    px_per_pt = (DPI / 72.0) * (display_px / (width_in * DPI))
    base, tick, title = 16.0 / px_per_pt, 14.0 / px_per_pt, 18.0 / px_per_pt
    lw = max(1.8, 2.0 / px_per_pt)
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": base,
        "axes.titlesize": title, "axes.labelsize": base,
        "xtick.labelsize": tick, "ytick.labelsize": tick, "legend.fontsize": tick,
        "axes.linewidth": 1.25, "xtick.major.width": 1.0, "ytick.major.width": 1.0,
        "xtick.major.size": 4.0, "ytick.major.size": 4.0,
        "axes.spines.top": False, "axes.spines.right": False,
        "figure.facecolor": WHITE, "axes.facecolor": WHITE, "savefig.facecolor": WHITE,
        "lines.linewidth": lw, "lines.solid_capstyle": "round",
    })
    return {"base": base, "tick": tick, "title": title, "lw": lw}


def save(fig, stem: str):
    fig.savefig(OUT / f"{stem}.png", dpi=DPI, bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {stem}.png")


def grid_ok(fig, ax, color=DEEP_BLUE, alpha=0.18):
    ax.grid(alpha=alpha, lw=0.8)


# ==========================================================================
# 01  gate-exponent selection
# ==========================================================================
def parse_point(key: str) -> dict:
    """'n=4|mRNA=2|mat=32.5' -> dict"""
    out = {}
    for part in key.split("|"):
        k, _, v = part.partition("=")
        out[k.strip()] = float(v)
    return out


def fig01():
    g = json.load(open(GRID_VERDICT, encoding="utf-8"))
    pp = g["per_point"]

    recs = []
    for key, v in pp.items():
        p = parse_point(key)
        r = dict(n=p["n"], mrna=p["mRNA"], mat=p["mat"], certified=bool(v["certified"]),
                 crossings=v.get("crossings"), gates=v.get("gates"))
        for c in CUTS:
            r[f"L_{c}"] = v.get(f"{c}_L_symmetric_median")
            r[f"L_{c}_min"] = v.get(f"{c}_L_symmetric_min")
            r[f"L_{c}_max"] = v.get(f"{c}_L_symmetric_max")
            r[f"old_{c}"] = v.get(f"{c}_old_aliased_leak_ratio")
        recs.append(r)
    d = pd.DataFrame(recs).sort_values(["n", "mrna", "mat"])
    ns = sorted(d["n"].unique())

    apply_style(14.0)
    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1.0], hspace=0.34, wspace=0.24)

    # -- A: L_symmetric vs n, three tail cuts -------------------------------
    a = fig.add_subplot(gs[0, :])
    styles = {"1pct": (TEAL, "o", "-"), "5pct": (DEEP_BLUE, "s", "-"), "10pct": (ORANGE, "^", "-")}
    for c in CUTS:
        col, mk, ls = styles[c]
        for n in ns:
            sub = d[d["n"] == n]
            a.plot([n] * len(sub), sub[f"L_{c}"], mk, color=col, ms=6, alpha=0.55,
                   markeredgecolor="none")
        med = [d[d["n"] == n][f"L_{c}"].median() for n in ns]
        a.plot(ns, med, ls, marker=mk, color=col, ms=8, lw=2.0,
               label=f"{c.replace('pct','')} % tail cut  (median of 4 grid points)")
    a.axvspan(N_SELECTED - 0.5, 8.5, color=LIGHT_TEAL, alpha=0.45, lw=0, zorder=0)
    a.axvline(4.5, color=RISK_ORANGE, ls="--", lw=1.6)
    a.axvline(N_SELECTED, color=DEEP_BLUE, ls=":", lw=2.2)
    a.set_yscale("log")
    a.set_xticks(ns)
    a.set_xlabel("gate exponent  $n_{A1,gate}$")
    a.set_ylabel("paired leak metric  $L_{symmetric}$  (log)")
    a.set_title("(A)  Why $n=6$: the paired leak metric orders $8<7<6<5\\ll4$ at every tail cut",
                loc="left")
    a.text(6.5, a.get_ylim()[1] * 0.35, "plateau\n$n=6$–$8$", ha="center", va="center",
           color=TEAL, fontweight="bold")
    a.text(4.5, a.get_ylim()[1] * 0.85, "  discrete pass boundary\n  ($n\\geq5$ passes)",
           ha="left", va="center", color=RISK_ORANGE, fontweight="bold")
    a.annotate("selected", xy=(N_SELECTED, a.get_ylim()[0] * 2.2),
               xytext=(N_SELECTED - 0.05, a.get_ylim()[0] * 1.15),
               ha="center", color=DEEP_BLUE, fontweight="bold")
    a.legend(loc="lower left", frameon=False)
    grid_ok(fig, a)

    # -- B: the 20-point certification grid ---------------------------------
    b = fig.add_subplot(gs[1, 0])
    combos = sorted({(r.mrna, r.mat) for r in d.itertuples()})
    for i, n in enumerate(ns):
        for j, (mr, mt) in enumerate(combos):
            row = d[(d["n"] == n) & (d["mrna"] == mr) & (d["mat"] == mt)]
            if row.empty:
                continue
            ok = bool(row["certified"].iloc[0])
            b.add_patch(Rectangle((i - 0.44, j - 0.44), 0.88, 0.88,
                                  fc=LIGHT_TEAL if ok else LIGHT_ORANGE,
                                  ec=DEEP_BLUE if ok else RISK_ORANGE, lw=1.4))
            b.text(i, j, "✓" if ok else "✗", ha="center", va="center",
                   color=TEAL if ok else RISK_ORANGE, fontweight="bold")
    b.set_xticks(range(len(ns))); b.set_xticklabels([f"{n:g}" for n in ns])
    b.set_yticks(range(len(combos)))
    b.set_yticklabels([f"mRNA {mr:g} min / mat {mt:g} min" for mr, mt in combos])
    b.set_xlabel("gate exponent  $n_{A1,gate}$")
    b.set_xlim(-0.6, len(ns) - 0.4); b.set_ylim(-0.6, len(combos) - 0.4)
    n4 = int(d[d["n"] == 4]["certified"].sum()); n58 = int(d[d["n"] >= 5]["certified"].sum())
    b.set_title(f"(B)  20-point grid: $n=4$ → {n4}/4 certified, $n\\geq5$ → {n58}/16", loc="left")

    # -- C: ordering stability across the three cuts ------------------------
    c = fig.add_subplot(gs[1, 1])
    ords = [g["ordering_1pct"], g["ordering_5pct"], g["ordering_10pct"]]
    for i, o in enumerate(ords):
        vals = [float(x) for x in o]
        c.plot(vals, [i] * len(vals), "-o", color=DEEP_BLUE, ms=8, lw=1.6)
        for v in vals:
            c.text(v, i + 0.16, f"{v:g}", ha="center", color=DEEP_BLUE)
    c.set_yticks([0, 1, 2]); c.set_yticklabels(["1 % cut", "5 % cut", "10 % cut"])
    c.set_xlabel("gate exponent, ordered best → worst")
    c.set_xlim(3.6, 8.4); c.set_ylim(-0.5, 2.7)
    c.set_title("(C)  Ordering identical at all three cuts\n"
                "(ordering_stable_across_cuts = %s)" % g["ordering_stable_across_cuts"], loc="left")
    grid_ok(fig, c)

    fig.suptitle("n = 6 gate-exponent selection — frozen point: $n_{A1,gate}=6$, "
                 "carry1 mRNA 2 min, maturation 32.5 min", x=0.01, ha="left")
    save(fig, "fig_n6_01_gate_exponent_selection")


# ==========================================================================
# 02  molecular-pool perturbation matrix (136 runs)
# ==========================================================================
POOL_ORDER = ["A1_full", "F1_full", "Int2_full", "RDF2_full", "C2", "S2"]
POOL_LABEL = {"A1_full": "$A_1$ (activator)", "F1_full": "$F_1$ (repressor)",
              "Int2_full": "Int2 (integrase)", "RDF2_full": "RDF2",
              "C2": "$C_2$ (complex)", "S2": "$S_2$ (DNA state)"}


def perturb_columns(d: pd.DataFrame):
    cols = []
    for mode in ("mult", "add", "flip"):
        sub = d[d["mode"] == mode]
        facs = sorted({f for f in sub["factor"].dropna().unique()})
        if mode == "flip":
            cols.append(("flip", np.nan))
        else:
            for f in facs:
                cols.append((mode, f))
    return cols


def fig02():
    d = pd.read_csv(POOL_CSV)
    cols = perturb_columns(d)

    apply_style(14.0)
    fig = plt.figure(figsize=(14, 10.5))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.35, 1.0], hspace=0.36, wspace=0.22)

    def col_label(mode, f):
        if mode == "flip":
            return "flip"
        return f"{f:g}×" if mode == "mult" else f"{f:+g}"

    # -- A: outcome matrix --------------------------------------------------
    a = fig.add_subplot(gs[0, :])
    for j, kind in enumerate(POOL_ORDER):
        pass
    for i, kind in enumerate(POOL_ORDER):
        for j, (mode, f) in enumerate(cols):
            if mode == "flip":
                sub = d[(d["kind"] == kind) & (d["mode"] == "flip")]
            else:
                sub = d[(d["kind"] == kind) & (d["mode"] == mode) & (np.isclose(d["factor"], f))]
            n = len(sub)
            if n == 0:
                a.add_patch(Rectangle((j - 0.45, i - 0.45), 0.9, 0.9,
                                      fc=WHITE, ec=BLUE_GRAY, lw=1.0, ls=":"))
                continue
            nrec = int((sub["outcome"] == "recovered_original").sum())
            nshift = int((sub["outcome"] == "legal_mod8_shift").sum())
            nlost = n - nrec - nshift
            fc = LIGHT_TEAL if nlost == 0 else LIGHT_ORANGE
            a.add_patch(Rectangle((j - 0.45, i - 0.45), 0.9, 0.9,
                                  fc=fc, ec=DEEP_BLUE, lw=1.3))
            a.text(j, i + 0.10, f"{nrec}/{n}", ha="center", va="center",
                   color=INK, fontweight="bold")
            if nshift:
                a.text(j, i - 0.22, f"+{nshift} shift", ha="center", va="center",
                       color=DEEP_BLUE)
            elif nlost:
                a.text(j, i - 0.22, f"{nlost} lost", ha="center", va="center",
                       color=RISK_ORANGE, fontweight="bold")
    a.set_xticks(range(len(cols)))
    a.set_xticklabels([col_label(m, f) for m, f in cols])
    a.set_yticks(range(len(POOL_ORDER)))
    a.set_yticklabels([POOL_LABEL[k] for k in POOL_ORDER])
    a.set_xlabel("perturbation applied to the pool")
    a.set_xlim(-0.6, len(cols) - 0.4); a.set_ylim(len(POOL_ORDER) - 0.4, -0.6)
    a.set_title("(A)  136 pool-perturbation runs — every cell keeps counting "
                "(green = all 8 initial states recovered; '+4 shift' = legal mod-8 phase shift)",
                loc="left")

    # -- B: what the S2 flip does ------------------------------------------
    b = fig.add_subplot(gs[1, 0])
    bits_lab = [f"{int(v):03d}" for v in
                sorted(d[(d["kind"] == "S2") & (d["mode"] == "flip")]["init_bits"].unique())]
    before, after, ok_shift = [], [], []
    for lab in bits_lab:
        ib = int(lab)
        base_row = d[(d["kind"] == "S2") & (d["mode"] == "mult") &
                     (np.isclose(d["factor"], 1.2)) & (d["init_bits"] == ib)]
        flip_row = d[(d["kind"] == "S2") & (d["mode"] == "flip") & (d["init_bits"] == ib)]
        if base_row.empty or flip_row.empty:
            continue
        d0 = int(str(base_row["steady_sequence"].iloc[0])[0])
        d1 = int(str(flip_row["steady_sequence"].iloc[0])[0])
        before.append(d0); after.append(d1)
        ok_shift.append((d1 - d0) % 8 == 4)
    y = np.arange(len(bits_lab))
    for i in range(len(bits_lab)):
        b.annotate("", xy=(after[i], i), xytext=(before[i], i),
                   arrowprops=dict(arrowstyle="-|>", color=RISK_ORANGE, lw=1.8,
                                   connectionstyle="arc3,rad=-0.35"))
        b.plot([before[i]], [i], "o", color=DEEP_BLUE, ms=10)
        b.plot([after[i]], [i], "o", color=ORANGE, ms=10)
        b.text(before[i], i + 0.30, f"{before[i]}", ha="center", color=DEEP_BLUE,
               fontweight="bold")
        b.text(after[i], i + 0.30, f"{after[i]}", ha="center", color=RISK_ORANGE,
               fontweight="bold")
    b.set_yticks(y); b.set_yticklabels(bits_lab)
    b.set_xticks(range(8))
    b.set_xlabel("first decoded digit after the perturbation")
    b.set_ylabel("initial state")
    b.set_xlim(-0.6, 7.6); b.set_ylim(len(bits_lab) - 0.4, -0.7)
    b.set_title("(B)  A full $S_2$ flip moves the count by exactly +4 mod 8\n"
                "in %d/%d initial states — a legal phase shift, not a failure"
                % (int(sum(ok_shift)), len(ok_shift)), loc="left")
    grid_ok(fig, b)

    # -- C: outcome totals --------------------------------------------------
    c = fig.add_subplot(gs[1, 1])
    oc = d["outcome"].value_counts()
    labels = ["recovered_original", "legal_mod8_shift"]
    vals = [int(oc.get(k, 0)) for k in labels]
    c.barh([0, 1], vals, color=[TEAL, ORANGE], alpha=0.9, height=0.55)
    for i, v in enumerate(vals):
        c.text(v + 2, i, str(v), va="center", color=INK, fontweight="bold")
    c.set_yticks([0, 1])
    c.set_yticklabels(["recovered\noriginal phase", "legal mod-8\nphase shift"])
    c.set_xlabel("number of runs")
    c.set_xlim(0, max(vals) * 1.22)
    c.set_ylim(-0.6, 1.6)
    c.set_title("(C)  0 lost lock, 0 readout-ok-but-causal-fail\n"
                "0 same-type adjacent carry in 136/136 runs", loc="left")
    grid_ok(fig, c)

    fig.suptitle("n = 6 molecular-pool perturbations — 136 runs, all eight initial states",
                 x=0.01, ha="left")
    save(fig, "fig_n6_02_pool_perturbation_matrix")


# ==========================================================================
# 03  eight phase-consistent initial states
# ==========================================================================
def fig03():
    e_init = json.load(open(EIGHT_INITIAL, encoding="utf-8"))     # per_value, spread
    e_pair = json.load(open(EIGHT_VERDICT, encoding="utf-8"))     # per_point (paired metric)

    pv = e_init["per_value"]
    pp = e_pair["per_point"]
    keys = sorted(pv.keys())                                      # '000' .. '111'

    sv = json.load(open(SELECTED, encoding="utf-8"))
    plc = sv["eight_initial_states"]["paired_leak_consistency"]

    paired = np.array([pp[f"init={k}"]["5pct_L_symmetric_median"] for k in keys], dtype=float)
    aliased = np.array([pv[k]["leak5pct_ratio"] for k in keys], dtype=float)

    apply_style(14.0)
    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1.0], hspace=0.40, wspace=0.24)

    # -- A: the eight decoded sequences as a rotation diagram ---------------
    a = fig.add_subplot(gs[0, :])
    seqs = [str(pv[k]["observed"]) for k in keys]
    base = seqs[0]
    for i, s in enumerate(seqs):
        off = next((j for j in range(len(base)) if s[:30] == base[j:j + 30]), 0)
        for j, ch in enumerate(s[:36]):
            want = base[(j + off) % len(base)]
            col = LIGHT_TEAL if ch == want else LIGHT_ORANGE
            a.add_patch(Rectangle((1 + j * 0.55, i - 0.34), 0.5, 0.68, fc=col, ec="none"))
            a.text(1 + j * 0.55 + 0.25, i, ch, ha="center", va="center", color=INK)
        a.text(-0.15, i, keys[i], ha="right", va="center", color=INK, fontweight="bold")
        a.text(1 + 36 * 0.55 + 0.25, i, f"starts at {s[0]}", ha="left", va="center", color=INK)
    a.set_xlim(-0.9, 1 + 36 * 0.55 + 5.4)
    a.set_ylim(-0.8, len(keys) - 0.2)
    a.invert_yaxis()
    a.axis("off")
    a.set_title("(A)  The eight phase-consistent initial states 000–111 — each is a pure "
                "rotation of the same strictly mod-8 sequence  (8/8 certified,\n"
                "      56 reads each, commitment 1.0 for every bit)", loc="left")

    # -- B: the aliasing trap, shown on purpose ----------------------------
    b = fig.add_subplot(gs[1, 0])
    x = np.arange(len(keys))
    b.semilogy(x, aliased, "o-", color=RISK_ORANGE, ms=8, lw=1.8,
               label="retired single-window ratio\n(aliased — must not be used)")
    b.axhspan(plc["min"], plc["max"], color=LIGHT_TEAL, alpha=0.9, zorder=0)
    b.axhline(plc["median"], color=TEAL, lw=2.6,
              label=f"paired $L_{{symmetric}}$ (5 % cut)\n"
                    f"median {plc['median']:.5f}, relative spread "
                    f"{plc['relative_spread']*100:.3f} %")
    b.set_xticks(x); b.set_xticklabels(keys, rotation=45)
    b.set_xlabel("initial state"); b.set_ylabel("leak metric (log scale)")
    b.set_title("(B)  The retired metric spans $3.5\\times10^4$ for one and the same\n"
                "circuit; the paired metric spans 0.196 %", loc="left")
    b.legend(loc="upper left", frameon=False)
    grid_ok(fig, b)

    # -- C: margins per bit, per initial state ------------------------------
    c = fig.add_subplot(gs[1, 1])
    w = 0.26
    for j, (bit, col, lab) in enumerate((("bit0", DEEP_BLUE, "bit 0"),
                                         ("bit1", TEAL, "bit 1"),
                                         ("bit2", ORANGE, "bit 2"))):
        vals = [float(pv[k][f"{bit}_margin_to_band_low"]) for k in keys]
        c.bar(x + (j - 1) * w, vals, width=w, color=col, alpha=0.9, label=lab)
    c.set_xticks(x); c.set_xticklabels(keys, rotation=45)
    c.set_xlabel("initial state")
    c.set_ylabel("margin to the low-read band")
    c.set_title("(C)  Per-bit read margin in every initial state\n"
                "(commitment = 1.0 for all bits in all eight)", loc="left")
    c.legend(frameon=False, ncol=3)
    grid_ok(fig, c)

    fig.suptitle("n = 6 eight phase-consistent initial states — one stable period-8 attractor",
                 x=0.01, ha="left")
    save(fig, "fig_n6_03_eight_initial_states")


# ==========================================================================
# 04  strict tolerance, 5 points x 3 levels
# ==========================================================================
DELTA_COLS = [
    ("5pct_L_symmetric_delta", "paired leak $L_{symmetric}$"),
    ("drop8_bit2_margin_to_band_low_delta", "bit2 margin to low band"),
    ("setup_min_h_global_delta", "setup margin (h)"),
    ("hold_min_h_global_delta", "hold margin (h)"),
    ("drop8_bit2_commitment_min_delta", "bit2 commitment"),
]


def fig04():
    d = pd.read_csv(TOL_CSV)
    levels = ["baseline", "tight", "ultra"]
    labels = list(dict.fromkeys(d["label"].tolist()))

    apply_style(14.0)
    fig = plt.figure(figsize=(14, 9.8))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0], hspace=0.42, wspace=0.26)

    # -- A: five points x three levels, with the negative control marked ----
    a = fig.add_subplot(gs[0, 0])
    for i, lab in enumerate(labels):
        sub = d[d["label"] == lab]
        nval = float(sub["n_A1_gate"].iloc[0])
        for j, lv in enumerate(levels):
            r = sub[sub["level"] == lv]
            if r.empty:
                continue
            r = r.iloc[0]
            ok = bool(r["certified"])
            a.add_patch(Rectangle((j - 0.44, i - 0.42), 0.88, 0.84,
                                  fc=LIGHT_TEAL if ok else LIGHT_ORANGE,
                                  ec=DEEP_BLUE if ok else RISK_ORANGE, lw=1.4))
            a.text(j, i + 0.10, "✓" if ok else "✗", ha="center", va="center",
                   color=TEAL if ok else RISK_ORANGE, fontweight="bold")
            a.text(j, i - 0.20, f"{int(r['crossings'])}/14", ha="center", va="center",
                   color=INK)
    a.set_xticks(range(3)); a.set_xticklabels([lv.capitalize() for lv in levels])
    a.set_yticks(range(len(labels)))
    a.set_yticklabels([f"{lab}\n($n$={d[d['label']==lab]['n_A1_gate'].iloc[0]:g})"
                       for lab in labels])
    a.set_xlabel("tolerance level")
    a.set_xlim(-0.6, 2.6); a.set_ylim(len(labels) - 0.4, -0.6)
    a.set_title("(A)  5 points × 3 levels, 600 h each.  The $n$=4 point is a deliberate\n"
                "negative control: tighter integration never rescues it", loc="left")

    # -- B: metric deltas vs baseline (passing points only) -----------------
    b = fig.add_subplot(gs[0, 1])
    passing = [l for l in labels if bool(d[d["label"] == l]["certified"].all())]
    x = np.arange(len(DELTA_COLS))
    w = 0.26
    for j, lv in enumerate(levels):
        vals, zeros = [], 0
        for c, _ in DELTA_COLS:
            sub = d[(d["level"] == lv) & (d["label"].isin(passing))][c].abs().to_numpy(float)
            v = np.nanmax(sub) if len(sub) else np.nan
            if np.isfinite(v) and v == 0.0:
                zeros += 1
            vals.append(v if (np.isfinite(v) and v > 0) else np.nan)
        b.bar(x + (j - 1) * w, vals, width=w, color=[BLUE_GRAY, TEAL, ORANGE][j],
              alpha=0.9, label=f"{lv}  ({zeros} metric(s) exactly 0)")
    b.set_yscale("log")
    b.set_xticks(x)
    b.set_xticklabels([lab for _, lab in DELTA_COLS], rotation=22, ha="right")
    b.set_ylabel("max |difference from baseline|")
    b.set_title("(B)  Numerically converged: the largest change in any reported metric\n"
                "is $\\sim10^{-10}$, far below every decision margin", loc="left")
    b.legend(frameon=False, ncol=1, loc="lower left")
    grid_ok(fig, b)

    # -- C: prefix trajectory gap ------------------------------------------
    c = fig.add_subplot(gs[1, 0])
    for lab in labels:
        sub = d[d["label"] == lab]
        ys = [max(float(sub[sub["level"] == lv]["prefix_traj_gap"].iloc[0]), 1e-12)
              for lv in levels] if len(sub) == 3 else []
        if not ys:
            continue
        ls = "--" if lab == "certified_failure" else "-"
        c.plot(range(3), ys, "o" + ls, ms=7, lw=1.6, label=lab)
    c.axhline(1e-12, color=BLUE_GRAY, lw=0.0)
    c.set_yscale("log")
    c.set_ylim(1e-11, 1e-5)
    c.set_xticks(range(3)); c.set_xticklabels([lv.capitalize() for lv in levels])
    c.set_ylabel("|prefix trajectory gap|  (51-state vs frozen 34-state)")
    c.set_xlabel("tolerance level")
    c.set_title("(C)  Baseline is the reference itself (gap 0).  Tightening the integrator\n"
                "leaves the gap at $\\sim4.2\\times10^{-7}$ — it saturates, it does not grow",
                loc="left")
    c.legend(frameon=False, fontsize=plt.rcParams["legend.fontsize"] * 0.85)
    grid_ok(fig, c)

    # -- D: identity flags --------------------------------------------------
    dd = fig.add_subplot(gs[1, 1])
    flags = [c for c in d.columns if c.endswith("_identical")]
    for i, fl in enumerate(flags):
        frac = float(d[fl].astype(bool).mean())
        dd.barh(i, frac, color=LIGHT_TEAL if frac == 1.0 else RISK_ORANGE,
                edgecolor=DEEP_BLUE, lw=1.0, height=0.6)
        dd.text(min(frac, 1.0) - 0.02, i, f"{frac*100:.0f} %", va="center", ha="right",
                color=INK, fontweight="bold")
    dd.set_yticks(range(len(flags)))
    dd.set_yticklabels([f.replace("_", " ") for f in flags])
    dd.set_xlim(0, 1.08); dd.set_xticks([0, 0.5, 1.0])
    dd.set_xlabel("fraction of the 15 runs matching the baseline")
    dd.set_title("(D)  Same code string, same events, same R/F pattern and the same\n"
                 "certification verdict in 15/15 runs", loc="left")
    grid_ok(fig, dd)

    fig.suptitle("n = 6 strict numerical tolerance — 15/15 runs stable; conclusions are "
                 "invariant to the integrator settings", x=0.01, ha="left")
    save(fig, "fig_n6_04_strict_tolerance")


# ==========================================================================
def main():
    print("n = 6 parameter-analysis and robustness figures")
    for p in (GRID_VERDICT, EIGHT_VERDICT, POOL_CSV, TOL_CSV):
        print("   ", "OK " if p.exists() else "MISSING", p.name)
    print("    output:", OUT)
    print()
    fig01()
    fig02()
    fig03()
    fig04()
    print("\ndone.")


if __name__ == "__main__":
    main()
