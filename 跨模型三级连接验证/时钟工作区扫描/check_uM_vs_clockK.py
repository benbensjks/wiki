"""uM×f 与 clock_K×f 是否**在判读所看的量上**等价？

上一支脚本（`check_uM_is_dynamical.py`）已经确认两件事：
  1. `receiver_uM_per_au` **是活参数**（`max|Δy|` 远非常数 0），
     `hybrid_model.py:160` 的 "metadata only" 注释是错的；
  2. `Int0` 峰值几乎精确按 `1/uM` 缩放（×0.8 → ×1.2499，×1.5 → ×0.6667，
     与 1/f 相符到 0.07%）。

但 `max|Δy|` 是**全 34 态**的最大值，会被前缀里按 `1/uM` 整体缩放的
mRNA/蛋白态主导（那些态改 uM 就整体平移，与判读无关）。判读真正看的是：
  * 读窗时刻（由 `Int0` 的峰/谷时刻决定）
  * 三个开关态 S0 / S1 / S2
  * 末级门信号 g1（时钟门唯一吃 Int0 的地方）

所以本脚本只量这几个量。若 uM×f 与 clock_K×f 在这些量上的差 ≪ 各自对基线的差，
则 uM 轴在判读层面**近似简并**于 clock_K 轴 —— 那么再单独扫一遍 uM 就是把
同一件事做第二遍。
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import common as C                                                  # noqa: E402

BASE_UM = 5.75
BASE_K = 0.30
HOURS = 60.0
FACTORS = (0.5, 0.8, 1.5)


def peaks(t, x):
    from scipy.signal import find_peaks
    amp = float(np.ptp(x))
    if amp <= 1e-12:
        return np.array([])
    return t[find_peaks(x, prominence=0.1 * amp)[0]]


def main():
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8')
    han = C.HanInput(HOURS)

    def run(kind, **kw):
        m, _ = C.build(kind, han, **kw)
        t, y = C.integrate(m, hours=HOURS)
        s = C.sig_of(m, y, kind)
        return t, y, s

    for kind in ('hzh', 'hzz'):
        t0, y0, s0 = run(kind)
        i0 = s0['int0']
        p0 = peaks(t0, i0)
        tail = slice(17, None)
        print('=' * 78)
        print(f'{kind}: 基线 uM={BASE_UM} clock_K={BASE_K}  '
              f'Int0峰值={np.max(i0):.5f}  峰数={len(p0)}')
        print('=' * 78)
        print(f'{"f":>5s} | {"uM×f vs 基线":>34s} | {"clock_K×f vs 基线":>34s} | '
              f'{"uM×f vs clock_K×f":>30s}')
        print(f'{"":>5s} | {"dTail":>8s} {"dSwitch":>8s} {"dg1":>8s} {"dPeak":>6s} | '
              f'{"dTail":>8s} {"dSwitch":>8s} {"dg1":>8s} {"dPeak":>6s} | '
              f'{"dTail":>8s} {"dSwitch":>8s} {"dg1":>8s} {"dPeak":>6s}')
        print('-' * 108)

        for f in FACTORS:
            tu, yu, su = run(kind, uM=BASE_UM * f)
            tk, yk, sk = run(kind, clock_K=BASE_K * f)

            def met(y, s, t):
                sw = np.max(np.abs(np.vstack([s['S0'] - s0['S0'],
                                              s['S1'] - s0['S1'],
                                              s['S2'] - s0['S2']])))
                dl = float(np.max(np.abs(s['g1'] - s0['g1'])))
                p = peaks(t, s['int0'])
                dp = (float(np.max(np.abs(p - p0))) if len(p) == len(p0) and len(p)
                      else float('nan'))
                return (float(np.max(np.abs(y[tail] - y0[tail]))), float(sw), dl, dp)

            a = met(yu, su, tu)
            b = met(yk, sk, tk)

            # 两者之差：以 uM×f 那一侧为参照
            def cross(y1, s1, t1, y2, s2, t2):
                sw = np.max(np.abs(np.vstack([s1['S0'] - s2['S0'],
                                              s1['S1'] - s2['S1'],
                                              s1['S2'] - s2['S2']])))
                p1, p2 = peaks(t1, s1['int0']), peaks(t2, s2['int0'])
                dp = (float(np.max(np.abs(p1 - p2))) if len(p1) == len(p2) and len(p1)
                      else float('nan'))
                return (float(np.max(np.abs(y1[tail] - y2[tail]))), float(sw),
                        float(np.max(np.abs(s1['g1'] - s2['g1']))), dp)

            c = cross(yu, su, tu, yk, sk, tk)
            print(f'{f:5.2f} | {a[0]:8.2e} {a[1]:8.2e} {a[2]:8.2e} {a[3]:6.3f} | '
                  f'{b[0]:8.2e} {b[1]:8.2e} {b[2]:8.2e} {b[3]:6.3f} | '
                  f'{c[0]:8.2e} {c[1]:8.2e} {c[2]:8.2e} {c[3]:6.3f}')
        print()

    print('读法：看最后一栏（uM×f vs clock_K×f）是否与左边两栏同量级。')
    print('  同量级 → 两条轴在判读层面不等价，单独扫 uM 有意义；')
    print('  小得多 → uM 只是"换了个说法的 clock_K"，uM 轴应当折算成 r = K/Int0_peak，')
    print('           不要再单列一条轴（否则同一件事会被数两遍）。')


if __name__ == '__main__':
    main()
