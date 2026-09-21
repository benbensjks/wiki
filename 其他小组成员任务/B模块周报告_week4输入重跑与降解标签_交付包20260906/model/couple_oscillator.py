"""couple_oscillator.py — couple A-module's phiC31 output to the B counter.

Week4 定稿接口（v53g 语义，与韩亚轩 Fixed_PLtetO1_Interface_Design_Scan.py 一致）:
    int_production(t) = RBS_scale * unloaded_C31_translation_flux(t)     [µM/h]

A 模块冻结参考输入:
    Week4_振荡器-C31建模.zip / data/reference_inputs/
    unloaded_C31_translation_trajectory.csv   (2026-08-14)
列 actual_C31_translation_flux_uM_h = rho(t)·beta31·m31(t)，unloaded 表示无 B 负载
回馈（resource_availability ≡ 1）。RBS_scale 乘在通量上 = C31 RBS 相对翻译强度。
Int 的结合、稀释与标签清除全部由 B 模块（zhao_core）自己处理，避免重复积分。

对照口径（画图/复核用）：同一 CSV 的 C31_protein 列（copies/cell，602 copies = 1 µM）。
旧 v36 蛋白轨迹 → 源项 的近似接口 (make_int_source / load_c31_spline) 保留，
但 Week4 §6.1 明确通量口径优先（蛋白轨迹近似会把 Int 的标签/稀释双算）。
"""
import os
import numpy as np

COPIES_PER_UM = 602.0  # 1 µM in 10^-15 L E. coli

_HERE = os.path.dirname(os.path.abspath(__file__))
CSV_PATH = os.path.join(_HERE, '..', 'docs', 'week4_unloaded',
                        'unloaded_C31_translation_trajectory.csv')

_CACHE = {}


def _load():
    if 'data' not in _CACHE:
        import pandas as pd
        df = pd.read_csv(CSV_PATH)
        t_h = df['time_min'].values / 60.0
        flux = df['actual_C31_translation_flux_uM_h'].values.astype(float)
        c31 = df['C31_protein'].values.astype(float)
        _CACHE['data'] = (t_h, flux, c31)
    return _CACHE['data']


def load_flux_spline(scale=1.0):
    """C31 产生通量 spline (µM/h)，×RBS_scale。同 scale 各级进程缓存。"""
    from scipy.interpolate import CubicSpline
    t_h, flux, _ = _load()
    key = round(float(scale), 6)
    if key not in _CACHE:
        _CACHE[key] = CubicSpline(t_h, key * flux, bc_type='natural')
    return _CACHE[key]


def make_flux_source(scale=1.0):
    """v53g 接口: int_production(t) = max(0, scale·flux(t))，单位 µM/h。"""
    spl = load_flux_spline(scale)

    def source(t):
        return max(0.0, float(spl(t)))

    return source


def load_c31_spline(scale=1.0):
    """A 蛋白轨迹 spline (µM)（对照口径）。"""
    from scipy.interpolate import CubicSpline
    t_h, _, c31 = _load()
    return CubicSpline(t_h, scale * c31 / COPIES_PER_UM, bc_type='natural')


def make_int_source(k_dil, scale=1.0):
    """旧 v36 近似口径：S(t) = dC/dt + k_dil·C（蛋白轨迹反推源项）。仅作对照。"""
    from scipy.interpolate import CubicSpline
    spl = load_c31_spline(scale)
    dspl = CubicSpline.derivative(spl)

    def source(t):
        return max(0.0, float(dspl(t)) + k_dil * float(spl(t)))

    return source


def pulse_times(mode='flux', tail_frac=0.5, prominence_frac=0.02):
    """稳态段（后 tail_frac）峰/谷时刻 (h)。mode='flux' 用通量（推荐采样口径），
    'c31' 用蛋白。返回 (peak_times_h, trough_times_h)。"""
    from scipy.signal import find_peaks
    t_h, flux, c31 = _load()
    y = flux if mode == 'flux' else c31
    lo = int((1.0 - tail_frac) * len(y))
    yt = y[lo:]
    prom = prominence_frac * (np.ptp(yt) if np.ptp(yt) > 0 else 1.0)
    pk, _ = find_peaks(yt, prominence=prom)
    tr, _ = find_peaks(-yt, prominence=prom)
    return t_h[lo:][pk], t_h[lo:][tr]


# ---- 输入波形元数据(阶段1报告用) ---------------------------------------------
def waveform_summary():
    """稳态通量/蛋白关键指标 dict。"""
    t_h, flux, c31 = _load()
    mask = t_h > 0.5 * t_h[-1]
    pk, tr = pulse_times(mode='flux')
    if len(pk) < 2:
        return {'error': 'peaks<2'}
    period = float(np.diff(pk).mean())
    flux_v = flux[mask]
    return dict(
        duration_h=float(t_h[-1]),
        period_h=period,
        flux_peak_uM_h=float(flux_v.max()),
        flux_trough_uM_h=float(flux_v.min()),
        flux_peak_trough_ratio=float(flux_v.max() / max(flux_v.min(), 1e-9)),
        c31_peak_uM=float((c31[mask] / COPIES_PER_UM).max()),
        c31_trough_uM=float((c31[mask] / COPIES_PER_UM).min()),
        n_cycles=int(len(pk)),
    )
