"""
============================================================
heatmap2.py

RDF GSA heatmap.

This mirrors GSA_Oscillator-C31/heatmap.py, but it reads Sobol result files
saved by gsa2.py instead of requiring hand-entered ST values.
============================================================
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from parameters2 import MODE, SOURCE_MODE, problem


# ==============================
# Parameters and outputs
# ==============================

params = problem["names"]

outputs = [
    "bm3r1_delay_h",
    "bm3r1_design_score",
    "bm3r1_production",
    "bm3r1_ss",
    "bm3r1_threshold_margin",
    "toggle_score",
    "switching_contrast",
    "rdf_max",
    "rdf_mean",
    "int_peak",
    "rep_max",
]


# ==============================
# Load ST data
# ==============================

base_result_dir = os.path.join(os.path.dirname(__file__), "results", MODE, SOURCE_MODE)
ST_matrix = []
available_outputs = []

for output in outputs:
    path = os.path.join(
        base_result_dir,
        output,
        f"sobol_{MODE}_{SOURCE_MODE}_{output}.npz",
    )

    if not os.path.exists(path):
        print(f"Skip missing result : {path}")
        continue

    data = np.load(path, allow_pickle=True)
    names = [str(x) for x in data["names"]]

    if names != params:
        raise RuntimeError(
            f"Parameter order mismatch in {path}.\n"
            f"Expected {params}\n"
            f"Got      {names}"
        )

    ST_matrix.append(data["ST"])
    available_outputs.append(output)

if not ST_matrix:
    raise RuntimeError(
        "No Sobol result files found.\n"
        "Run gsa2.py once for each OUTPUT you want in the heatmap."
    )

ST_matrix = np.array(ST_matrix)


# ==============================
# Heatmap
# ==============================

plt.figure(figsize=(12, 7))

sns.heatmap(
    ST_matrix,
    annot=True,
    fmt=".2f",
    cmap="viridis",
    xticklabels=params,
    yticklabels=available_outputs,
)

plt.xlabel("Parameters")
plt.ylabel("Outputs")
plt.title(
    f"RDF Global Sensitivity Analysis\n"
    f"{MODE} / {SOURCE_MODE} / Total-effect Sobol Index"
)
plt.tight_layout()

heatmap_dir = os.path.join(base_result_dir, "heatmap")
os.makedirs(heatmap_dir, exist_ok=True)
fig_path = os.path.join(heatmap_dir, f"heatmap_{MODE}_{SOURCE_MODE}.png")
plt.savefig(fig_path, dpi=150)
print(f"Saved figure : {fig_path}")
plt.show()
