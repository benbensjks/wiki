"""couple_oscillator.py — couple A-module's φC31 oscillator waveform to the counter.

Loads v36_final_simulation.csv (A module, minutes, copies/cell), converts to
µM (602 copies = 1 µM in a 1 fL cell) and builds the integrase production
source S_int(t) = dC31/dt + k_dil·C31  (µM/h), so that free Int tracks the
oscillator's C31 protein trajectory when recombination consumption is small.
"""
import os
import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline

COPIES_PER_UM = 602.0  # 1 µM in 10^-15 L E. coli (A-module report, table row 浓度转换)

_HERE = os.path.dirname(os.path.abspath(__file__))
CSV_PATH = os.path.join(_HERE, '..', 'docs', '振荡器-C31_v3',
                        'v36_final_simulation.csv')


def load_c31_spline(scale=1.0):
    """Return (spline_uM(t_h), spline_deriv(t_h)) for the C31 protein in µM."""
    df = pd.read_csv(CSV_PATH)
    t_h = df['time_min'].values / 60.0
    c31_uM = df['C31_protein'].values / COPIES_PER_UM * scale
    spl = CubicSpline(t_h, c31_uM, bc_type='natural')
    return spl, spl.derivative()


def make_int_source(k_dil, scale=1.0):
    spl, dspl = load_c31_spline(scale)

    def source(t):
        return max(0.0, float(dspl(t)) + k_dil * float(spl(t)))

    return source, spl


# period / timing metadata of the waveform
def pulse_times():
    df = pd.read_csv(CSV_PATH)
    t = df['time_min'].values
    c = df['C31_protein'].values
    from scipy.signal import find_peaks
    pk, _ = find_peaks(c, prominence=200)
    tr, _ = find_peaks(-c, prominence=200)
    return t[pk] / 60.0, t[tr] / 60.0  # hours
