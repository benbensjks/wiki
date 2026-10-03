"""冒烟预检：在花 37 min 之前，用 ~15 s 证明整条调用链在两套模型上都能走通。

只做四件事，全部是"会不会当场抛异常"级别的检查，不产出任何科学结论：
  1. `structural_guard`（40 探针）——确认 HZZ 末级与原文件逐点一致
  2. 两套模型各建一次、各积分 20 h（短），确认 `build/integrate/sig_of/verdict` 通路
  3. `per_bit` 与 `mechanism` 在**短轨迹**上能否取到全部字段（长轨迹才会用到的字段也过一遍）
  4. 报出 HZZ 两条已被 `diagnose_round*` 用过的量（S2 翻转数、g1 剂量）作横向对照
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import common as C                                                  # noqa: E402
import scan_clock_axis as A                                         # noqa: E402


def main():
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8')
    print('HAN sha =', C.sha256(Path(C.HAN_PATH)), flush=True)
    han = C.HanInput(20.0)

    sg = C.structural_guard(C.HZHModel('han', han))
    print('structural :', sg, flush=True)

    pg = C.parity_guard(han, hours=40.0)
    print('parity(40h):', pg, flush=True)

    for kind, kw in (('hzh', {}), ('hzz', dict(autoregulation=False)),
                     ('hzz', dict(autoregulation=True))):
        m, meta = C.build(kind, han, **kw)
        t, y = C.integrate(m, hours=20.0, sample_min=1.0, max_step_min=1.0)
        sig = C.sig_of(m, y, kind)
        v = C.verdict(t, sig)
        pb = C.pb_per_bit(sig, v['read_windows'], t)
        mech = C.mechanism(t, sig, v['read_windows'])
        print(f'{kind}/{meta["arm"] if "arm" in meta else kw.get("autoregulation")}: '
              f'n={y.shape[0]} t={t[-1]:.1f}h windows={len(v["read_windows"])} '
              f'seq={v["steady"]["sequence"]} s2x={len(v["crossings"]["S2"])} '
              f'int0peak={np.max(sig["int0"]):.5f} '
              f'pb={ {b: (pb[b]["commitment"], pb[b]["band_margin"]) for b in C.PB_BITS} } '
              f'mech_cycles={len(mech["per_cycle"])} '
              f'g1dose={mech["per_cycle"][0]["g1_dose"] if mech["per_cycle"] else None} '
              f'i2r2={mech["int2_times_rdf2_peak"]:.4f}', flush=True)
        # 逐字段探一遍：任何一个是 NaN/None 都会在长跑落盘时炸掉 allow_nan=False
        for k, val in list(mech.items()) + list(v.items()):
            if isinstance(val, float) and not np.isfinite(val):
                raise RuntimeError(f'non-finite field {k}={val}')
        # 折行 + 打印路径：第一版把这段内联在 main 里，撞键 bug 在护栏跑完 4 min 后才炸
        arm = (kind == 'hzz' and kw.get('autoregulation', False))
        row = A.build_row(kind, arm, A.ARM_LABEL[kind][arm], 0.30, 2.0, 0.54848,
                          meta, v, pb, mech)
        line = A.fmt_row(1, 1, row, pb, mech)
        C.dump(HERE / '_preflight_row.json', row)     # 顺带验证 allow_nan=False 落盘
        print('row ok    :', line, flush=True)
        (HERE / '_preflight_row.json').unlink()
    print('PREFLIGHT OK', flush=True)


if __name__ == '__main__':
    main()
