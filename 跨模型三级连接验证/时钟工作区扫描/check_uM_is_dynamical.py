"""`receiver_uM_per_au` 到底进不进动力学？——直接量，不看注释。

起因：`hybrid_model.py:160` 给这个字段写的行内注释是
    receiver_uM_per_au: float = 5.75  # metadata only; does NOT calibrate donor concentrations
但 `model_hzh.py` 里它被用来把上游信号换算成接收端单位：
    29:  self.copies_per_au = 602.214076 * self.tail.c.receiver_uM_per_au
    37:  return 60*han.interp(...)/self.copies_per_au
    49:  tx = ...*interp(...)/self.copies_per_au
    81:  b0[0] = han.sol.y[6,0]/self.copies_per_au
若这几行真的生效，那么注释就是错的，而且**踩上这个错注释的人会得到"uM 扰动是空操作"的
假结论**（只碰 `HbyReceiver` 而不碰 `HZHModel.copies_per_au` 的 harness 正好是这样）。

本脚本只做一件事：同一条上游、同一个时钟配置，只把 uM 换掉，看轨迹是否**逐位相同**。
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import common as C                                                  # noqa: E402

BASE_UM = 5.75
FACTORS = (0.5, 0.8, 1.0, 1.5)
HOURS = 60.0


def main():
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8')
    han = C.HanInput(HOURS)

    ref = {}
    for kind in ('hzh', 'hzz'):
        m, meta = C.build(kind, han, uM=BASE_UM)
        t, y = C.integrate(m, hours=HOURS)
        ref[kind] = (t, y, m, meta, float(np.max(C.sig_of(m, y, kind)['int0'])))

    print(f'{"kind":5s} {"uM":>7s} {"x":>5s} {"max|dy|":>12s} {"Int0峰值":>10s} '
          f'{"比值":>8s} {"copies_per_au":>14s}')
    print('-' * 70)
    for kind in ('hzh', 'hzz'):
        t, y, m, meta, i0 = ref[kind]
        for f in FACTORS:
            um = BASE_UM * f
            m2, meta2 = C.build(kind, han, uM=um)
            if f == 1.0:
                t2, y2, i02 = t, y, i0
            else:
                t2, y2 = C.integrate(m2, hours=HOURS)
                i02 = float(np.max(C.sig_of(m2, y2, kind)['int0']))
            gap = float(np.max(np.abs(y2 - y)))
            print(f'{kind:5s} {um:7.3f} {f:5.2f} {gap:12.3e} {i02:10.5f} '
                  f'{i02 / i0:8.4f} {meta2["copies_per_au"]:14.3f}')
        print()

    # 同一时间网格吗？
    for kind in ('hzh', 'hzz'):
        t, y, m, meta, i0 = ref[kind]
        for f in FACTORS[1:]:
            m2, _ = C.build(kind, han, uM=BASE_UM * f)
            t2, y2 = C.integrate(m2, hours=HOURS)
            same_grid = t.shape == t2.shape and bool(np.all(t == t2))
            print(f'{kind} x{f}: 时间网格逐点相同 = {same_grid}')
        print()

    # ------------------------------------------------------------------
    # 第二部分：uM 与 clock_K 是否**近似简并**？
    # 论证：门只通过 hill(clock_scale*Int0, clock_K, n) 吃 Int0，
    # 而 hill(s*x,K,n) == hill(x,K/s,n) 精确成立。所以只要 uM 只把 Int0
    # 的**幅度**按 1/uM 缩放、波形形状不变，uM×f 就等价于 clock_K×f。
    # 上面测到 Int0 峰值比 = 1/uM 到 0.07% —— 说明形状几乎不变但并非完全不变。
    # 这里直接量：把两组本该等价的配置跑出来比。
    # ------------------------------------------------------------------
    print('=' * 70)
    print('uM ×f  vs  clock_K ×f  （若近似简并，max|dy| 应当很小）')
    print('=' * 70)
    print(f'{"kind":5s} {"f":>5s} {"uM×f 的 max|dy|":>18s} {"clock_K×f 的 max|dy|":>20s} '
          f'{"两者之差 max|dy|":>18s}')
    print('-' * 70)
    for kind in ('hzh', 'hzz'):
        m0, _ = C.build(kind, han)
        t0, y0 = C.integrate(m0, hours=HOURS)
        for f in (0.5, 0.8, 1.5):
            mu, _ = C.build(kind, han, uM=BASE_UM * f)
            tu, yu = C.integrate(mu, hours=HOURS)
            mk, _ = C.build(kind, han, clock_K=0.3 * f)
            tk, yk = C.integrate(mk, hours=HOURS)
            du = float(np.max(np.abs(yu - y0)))
            dk = float(np.max(np.abs(yk - y0)))
            duk = float(np.max(np.abs(yu - yk)))
            print(f'{kind:5s} {f:5.2f} {du:18.3e} {dk:20.3e} {duk:18.3e}')
        print()

    print('读法：')
    print('  * 若 max|dy| 在 uM 改动下 **恰为 0**，注释"metadata only"成立，uM 轴是空的；')
    print('  * 若 max|dy| > 0 且 Int0 峰值随 uM 单调移动，则该字段是**活参数**，')
    print('    注释错误，且 uM 是一条真实（且尚未标定）的输入标尺轴；')
    print('  * 若最后一行"两者之差"远小于前两行，则 uM 与 clock_K 近似简并 ——')
    print('    那么再单独扫一遍 uM 就是把同一件事做第二遍（除非要量那 0.07% 的形状差）。')


if __name__ == '__main__':
    main()
