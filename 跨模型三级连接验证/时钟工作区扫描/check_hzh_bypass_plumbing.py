"""消融（时钟门保留 vs 去掉）能不能照 HZZ 那一套做法直接搬到 34 态 HZH 上？

背景
----
`hby_zmh_zmh` 的 round 4/4b 用 `class BypassClock(D.HZZModel)` 覆写 `tail_signals`，
把 `clock` 强制为 1。那条路在 HZZ 上**能**走通，因为 HZZ 的判读函数签名是
`diagnose_round3.evaluate(t, sig, hours)` —— **sig 由调用方提供**，
所以 `mb.signals(y)` 里被改写的 `g1` 会一路传进判读。

34 态 HZH 不是这样。它的判读入口是 `verify_hzh.analyse(t, y, zmh, receiver)`，
而 `analyse` 在内部自己调 `verify_hzh.signals(y, zmh, receiver)` **重算** clock/g1/S0/S1/S2
——**从不读模型 `signals()` 返回的字典**。于是：

    动力学被旁路（因为 rhs 走的是模型的 signals），
    判读却仍在看**没被旁路的那条 g1**。

这正是 round 4 那个 bug 的形状（"动力学已旁路、上报信号未旁路"），
只不过 round 4 是靠 `q['g1']` 修好的，而在 HZH 上**改键根本修不了**，
因为判读器根本不看那个键。

本脚本只做四件事，不产出任何科学结论：
  1. 打印认证门的常数（核对 `g1 = H(A1;1.2,6)·(1−H(F1;0.4,4))·H(Int0;K,n)` 与 0.3/2）
  2. 造一个"只替换第三个因子"的旁路末级，跑 60 h，证明**动力学确实变了**
  3. 对**同一条被旁路的轨迹**调 `verify_hzh.analyse`，证明**判读完全看不见旁路**
  4. 给出可用的正确写法（把旁路后的 sig 显式传给 `common.verdict`，并保留 parity 守卫）
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import common as C                                                  # noqa: E402
import verify_hzh as VH                                             # noqa: E402
from hybrid_model import HbyReceiver, hill                          # noqa: E402
from model_hzh import HZHModel                                      # noqa: E402

HOURS = 60.0


class BypassTail(HbyReceiver):
    """只替换 g1 的第三个因子：clock ≡ 1。其余（act、repress、重组块）原样。"""

    def signals(self, y, int0):
        q = super().signals(y, int0)
        p, c = self.p, self.c
        a, f = y[0], y[1]
        act = hill(a, p['K_A'][1], c.n_A1_gate)
        repress = 1.0 - hill(f, p['K_F'][1], p['n_F'][1])
        gate = act * repress                          # 不含 clock
        q['clock_gate'] = 1.0
        q['g1'] = gate
        q['u2_target_au_per_h'] = p['alpha_Int'][1] * gate
        return q


def main():
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8')
    han = C.HanInput(HOURS)

    m0 = HZHModel('han', han)
    p, c = m0.p, m0.tail.c
    print('=== 1. 认证门的常数（来源：hybrid_model.load_zeng_table + HbyConfig） ===')
    for k in ('K_A', 'K_F', 'n_F', 'n_A', 'alpha_Int'):
        print(f'  {k}[1] = {p[k][1]!r}')
    print(f'  n_A1_gate   = {c.n_A1_gate}   (HbyConfig 默认 6.0)')
    print(f'  clock_K_au  = {c.clock_K_au}   clock_n = {c.clock_n}   '
          f'clock_scale = {c.clock_scale}')
    print(f'  => g1 = H(A1;{p["K_A"][1]},{c.n_A1_gate}) · '
          f'(1-H(F1;{p["K_F"][1]},{p["n_F"][1]})) · H(Int0;{c.clock_K_au},{c.clock_n})')
    print(f'  alpha_Int[1] = {p["alpha_Int"][1]} a.u./h')

    # ---------------------------------------------------------------- 2
    t0, y0 = C.integrate(m0, hours=HOURS)
    mb = HZHModel('han', han)
    mb.tail = BypassTail(mb.tail.c)
    mb.p = mb.tail.p
    tb, yb = C.integrate(mb, hours=HOURS)
    print()
    print('=== 2. 旁路是否真的改了动力学 ===')
    print(f'  max|y_bypass − y_normal| = {np.max(np.abs(yb - y0)):.6e}   '
          f'(>0 说明 rhs 确实用了旁路后的 signals)')
    s0 = VH.signals(yb, mb.z, mb.tail)
    print(f'  旁路轨迹上的 clock：min={np.min(s0["clock"]):.5f} '
          f'max={np.max(s0["clock"]):.5f}  (旁路应当恒为 1)')

    # ---------------------------------------------------------------- 3
    print()
    print('=== 3. 判读器看不看得见这次旁路 ===')
    v_norm, _ = VH.analyse(t0, y0, m0.z, m0.tail)
    v_byp, _ = VH.analyse(tb, yb, m0.z, m0.tail)          # 注意：用的是**原始** receiver
    v_byp2, _ = VH.analyse(tb, yb, mb.z, mb.tail)         # 用**旁路后的** receiver
    import json
    same_z = json.dumps(v_byp, sort_keys=True, default=str) == \
        json.dumps(v_byp2, sort_keys=True, default=str)
    print(f'  analyse(旁路轨迹, 原始 receiver) == analyse(旁路轨迹, 旁路 receiver) ? {same_z}')
    print(f'  -> 判读结果只取决于 (t, y, c, p)，与"末级是怎么算 g1 的"无关')
    print(f'     旁路轨迹的判读：counting={v_byp["steady"]["passed"]} '
          f'events={all(e["passed"] for e in v_byp["events"].values())} '
          f'seq={v_byp["steady"]["sequence"]}')
    print(f'     原始轨迹的判读：counting={v_norm["steady"]["passed"]} '
          f'events={all(e["passed"] for e in v_norm["events"].values())} '
          f'seq={v_norm["steady"]["sequence"]}')

    # 门信号对比：判读用的 g1 vs 旁路后的 g1
    g_vh = np.asarray(s0['g1'], dtype=float)
    g_bp = np.asarray([mb.tail.signals(yb[17:, j], yb[2, j])['g1']
                       for j in range(yb.shape[1])], dtype=float)
    clk = np.asarray(s0['clock'], dtype=float)
    j = int(np.argmax(g_vh))
    print()
    print(f'  在门峰值处（t={tb[j]:.3f} h）：')
    print(f'    判读器用的 g1 = {g_vh[j]:.6f}   （= act·repress·clock）')
    print(f'    旁路后的  g1 = {g_bp[j]:.6f}   （= act·repress）')
    print(f'    比值 = {g_bp[j]/g_vh[j]:.4f}   该点 clock = {clk[j]:.6f}')
    print(f'  整条轨迹：判读器 g1 峰值 {g_vh.max():.6f} vs 旁路 g1 峰值 {g_bp.max():.6f}')
    print(f'  两者最大差 = {np.max(np.abs(g_vh - g_bp)):.6e}')
    print('  => **判读器把被删掉的 clock 因子又乘了回去**：'
          '它数的是"门还在"的那条脉冲链。')

    # ---------------------------------------------------------------- 4
    print()
    print('=== 4. 正确写法（本仓库已有工具） ===')
    sig_bp = C.sig_of(mb, yb, 'hzh')       # 对 hzh 也是调 verify_hzh.signals -> 仍会重算
    sig_bp['clock'] = np.ones_like(clk)    # 覆写必须发生在 sig_of **之后**
    sig_bp['g1'] = g_bp
    v_ok = C.verdict(tb, sig_bp)
    print(f'  改成"先 sig_of，再显式覆写 clock/g1，再交给 common.verdict"后：')
    print(f'    counting={v_ok["steady"]["passed"]} '
          f'events={all(e["passed"] for e in v_ok["events"].values())} '
          f'certified={v_ok["certified_v1"]} '
          f'gate_events(s1)={v_ok["events"]["bit1_to_bit2"]["gate_events"]}')
    print(f'  （对照：判读器不知道旁路时的 gate_events(s1)='
          f'{v_byp["events"]["bit1_to_bit2"]["gate_events"]}）')
    print()
    print('  要点：`C.sig_of` 对 hzh 也是调 `verify_hzh.signals`，**同样会重算**，')
    print('        所以覆写必须发生在 sig_of 之后、verdict 之前。')
    print('        parity 守卫保证 `common.verdict` 与 `analyse` 在未旁路时逐字段一致，')
    print('        因此这条路径只是"把 sig 交出来"，不会引入第二套判读逻辑。')


if __name__ == '__main__':
    main()
