"""
============================================================
parameters2.py

RDF system Global Sensitivity Analysis

Keep the same role as GSA_Oscillator-C31/parameters1.py:
select model mode, select one scalar output, and define the SALib problem.
============================================================
"""

import os
import sys

MODEL_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "model"))
if MODEL_DIR not in sys.path:
    sys.path.append(MODEL_DIR)

import zhao_core as zc


# =====================================================
# Select model
# =====================================================

MODE = "rdf_waveform"

# MODE = "rdf_waveform"

# MODE = "rdf_core_binding"


# =====================================================
# Select input source
# =====================================================

SOURCE_MODE = "c31_waveform"

# SOURCE_MODE = "c31_waveform"


# =====================================================
# Select output
# =====================================================

OUTPUT = "rdf_mean"
"""
RDF system
----------
toggle_score
H_mean
L_mean
alternation_fidelity
switching_contrast
lr_final
pb_final
lr_mean
lr_amplitude
rdf_max
rdf_mean
rdf_auc
int_peak
int_trough
rep_max
rep_mean
bm3r1_mrna_ss
bm3r1_production
bm3r1_ss
bm3r1_delay_h
bm3r1_threshold_margin
bm3r1_design_score
"""


# =====================================================
# Simulation settings
# =====================================================

T_END = 180.0
DT = 0.05
# toggle_score needs at least TRANSIENT_SKIP + 4 period samples.
# With a 24 h input period, 180 h gives enough samples after transient removal.
T_END_BY_SOURCE = {
    "square": 180.0,
    "c31_waveform": 100.0,
}
TRANSIENT_SKIP = 3

# Saltelli base sample size. The RDF model is much heavier than the oscillator.
# Raise this to 128 or 256 for final runs after the output is selected.
N = 256
N_JOBS = -1


# =====================================================
# Baseline RDF parameter set
# =====================================================

P0 = zc.default_params()

# Operating point used by the RDF/counter scans.
BASELINE = dict(
    K_rep=0.0186,
    n_rep=3.4,
    krep_tsl=15.0,
    krdf_tsl=200.0,
    k_tag_int=12.0,
    baseline_sub=0.26,
    waveform_scale=1.0,
    k_int=P0["k_int"],
    pulse_width=P0["ara_off"] - P0["ara_on"],
    kt=P0["kt"],
)


# =====================================================
# RDF system parameters
# =====================================================

if MODE == "rdf_square":

    problem = {
        "num_vars": 8,
        "names": [
            "k_int",
            "pulse_width",
            "kt",
            "krep_tsl",
            "krdf_tsl",
            "k_tag_int",
            "K_rep",
            "n_rep",
        ],
        "bounds": [
            [0.5 * BASELINE["k_int"], 1.5 * BASELINE["k_int"]],
            [0.10, 0.60],
            [0.10, 0.60],
            [0.5 * BASELINE["krep_tsl"], 1.5 * BASELINE["krep_tsl"]],
            [0.5 * BASELINE["krdf_tsl"], 1.5 * BASELINE["krdf_tsl"]],
            [0.0, 20.0],
            [0.5 * BASELINE["K_rep"], 1.5 * BASELINE["K_rep"]],
            [2.5, 4.0],
        ],
    }

elif MODE == "bm3r1_promoter":

    problem = {
        "num_vars": 7,
        "names": [
            "krep_tsl",
            "K_rep",
            "n_rep",
            "krdf_tsl",
            "k_tag_int",
            "k_int",
            "pulse_width",
        ],
        "bounds": [
            [1.0, 120.0],
            [0.005, 0.050],
            [2.5, 4.0],
            [50.0, 400.0],
            [0.0, 20.0],
            [0.5 * BASELINE["k_int"], 1.5 * BASELINE["k_int"]],
            [0.10, 0.60],
        ],
    }

elif MODE == "rdf_waveform":

    problem = {
        "num_vars": 7,
        "names": [
            "baseline_sub",
            "waveform_scale",
            "krep_tsl",
            "krdf_tsl",
            "k_tag_int",
            "K_rep",
            "n_rep",
        ],
        "bounds": [
            [0.0, 0.30],
            [0.5, 2.0],
            [0.5 * BASELINE["krep_tsl"], 1.5 * BASELINE["krep_tsl"]],
            [0.5 * BASELINE["krdf_tsl"], 1.5 * BASELINE["krdf_tsl"]],
            [0.0, 20.0],
            [0.5 * BASELINE["K_rep"], 1.5 * BASELINE["K_rep"]],
            [2.5, 4.0],
        ],
    }

elif MODE == "rdf_core_binding":

    problem = {
        "num_vars": 8,
        "names": [
            "Ki",
            "Kii",
            "Kir",
            "Ks01",
            "Ks02",
            "Kb1",
            "Kb2",
            "Dtot",
        ],
        "bounds": [
            [0.5 * P0["Ki"], 1.5 * P0["Ki"]],
            [0.5 * P0["Kii"], 1.5 * P0["Kii"]],
            [0.5 * P0["Kir"], 1.5 * P0["Kir"]],
            [0.5 * P0["Ks01"], 1.5 * P0["Ks01"]],
            [0.5 * P0["Ks02"], 1.5 * P0["Ks02"]],
            [0.5 * P0["Kb1"], 1.5 * P0["Kb1"]],
            [0.5 * P0["Kb2"], 1.5 * P0["Kb2"]],
            [0.5 * P0["Dtot"], 1.5 * P0["Dtot"]],
        ],
    }

else:

    raise ValueError("Unknown MODE")


