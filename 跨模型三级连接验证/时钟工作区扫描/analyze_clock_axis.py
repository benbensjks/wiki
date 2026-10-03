"""把时钟轴扫描的 all.csv 折成「成功区在哪、边界由什么决定」的对照表。

输入：`results/clock_axis_*/`（取最新一个，或 `--dir` 指定）
输出：stdout 的 markdown 表 + 同目录 `summary.json` / `summary.md`

回答三个问题
------------
1. 每个 (模型, 臂, n) 的**成功 K 区间**在哪，两端分别因为什么失效
   —— 用 `counting_passed`（位本身有没有按 mod8 走）与 `event_causality_passed`
   （进位事件链是否一一对应）**分开**判定，这与 round 2 记录的"两种失效模式"一致。
2. `n=2 → n=3` 把边界移动了多少（门的陡度 vs 阈值位置）。
3. 同一个 `K` 下，HZZ「自馈关」（源文件真实行为）与「自馈开」（已认证行为）
   差多少；以及 HZH 与 HZZ 在**同一个无量纲比值 r = K / Int0_peak** 上是否一致。

注意：两套模型的时钟 Hill 指数不同（HZH n=2、HZZ 硬编码 n=3），
所以 `K` 的绝对值**不能直接排名**，可迁移的量是 `r`。
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = HERE / 'results'

KIND_ORDER = ('hzh', 'hzz')
ARM_ORDER = {'hzh': ('base',), 'hzz': ('off', 'on')}


def load(dirpath: Path):
    """优先读最终 `all.json`；扫描还没跑完时退回 `partial.json`（字段类型一致）。

    两者都是原生 JSON 类型（bool/float/dict），不走 csv —— csv 那一路会把
    `mechanism` 存成字符串、把布尔存成 'True'/'False'，是纯粹的自找麻烦。
    """
    for name in ('all.json', 'partial.json'):
        p = dirpath / name
        if p.exists():
            blob = json.loads(p.read_text(encoding='utf-8'))
            if blob.get('rows'):
                return blob['rows'], blob
    raise SystemExit(f'no all.json / partial.json with rows under {dirpath}')


def classify(r):
    if r['certified_v1']:
        return 'PASS'
    if r['counting_passed'] and not r['event_causality_passed']:
        return 'events-only'          # 位在走，但进位链对不上
    if not r['counting_passed'] and r['event_causality_passed']:
        return 'counting-only'
    return 'both'


def window(rows):
    """成功 K 区间（按 certified_v1）。返回 (lo, hi, 全部 K 列表)。"""
    ks = sorted(r['K'] for r in rows)
    good = [r['K'] for r in rows if r['certified_v1']]
    return (min(good) if good else None, max(good) if good else None, ks)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dir', default=None)
    args = ap.parse_args()
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8')

    d = Path(args.dir) if args.dir else sorted(RESULTS.glob('clock_axis_*'))[-1]
    rows, top = load(d)
    int0 = top['int0_peak']
    gd = top.get('guards') or {}
    print(f'dir        = {d.name}')
    print(f'Int0_peak  = {int0:.5f}   ->  r = K / {int0:.5f}')
    print(f'rows       = {len(rows)}'
          + ('' if not gd else
             f'   guards: parity {gd["parity"]["ok"]}, '
             f'prefix_gap {gd["structural"]["prefix_rhs_gap"]}, '
             f'repro_hzh {gd["reproduction_hzh"]["ok"]}, '
             f'repro_hzz {all(v["ok"] for v in gd["reproduction_hzz"].values())}'))
    if not gd:
        print('  (partial: 护栏记录只在最终 all.json 里，本表可能不完整)')
    print()

    out = {'dir': d.name, 'int0_peak': int0, 'groups': {}, 'guards': gd}

    # ---------------- 1. 成功区 + 两端失效原因 -----------------------------
    print('| 模型 | 臂 | n | 成功 K 区间 | 对应 r 区间 | 高 K 端失效方式 | 低 K 端 |')
    print('|---|---|---|---|---|---|---|')
    for kind in KIND_ORDER:
        for arm in ARM_ORDER[kind]:
            for n in (2.0, 3.0):
                g = [r for r in rows
                     if r['kind'] == kind and r['arm'] == arm and r['n'] == n]
                if not g:
                    continue
                lo, hi, ks = window(g)
                hi_row = next((r for r in g if r['K'] == max(ks)), None)
                fail_hi = classify(hi_row) if hi_row and not hi_row['certified_v1'] else '-'
                fail_lo = None
                lo_row = next((r for r in g if r['K'] == min(ks)), None)
                if lo_row and not lo_row['certified_v1']:
                    fail_lo = classify(lo_row)
                rtxt = '-' if lo is None else f'{lo/int0:.3f} ~ {hi/int0:.3f}'
                print(f'| {kind} | {arm} | {n:g} | '
                      f'{"无" if lo is None else f"{lo:g} ~ {hi:g}"} | {rtxt} | '
                      f'{fail_hi} | {fail_lo or "通过"} |')
                out['groups'][f'{kind}|{arm}|n={n:g}'] = dict(
                    k_min=min(ks), k_max=max(ks), pass_lo=lo, pass_hi=hi,
                    n_pass=sum(r['certified_v1'] for r in g), n_rows=len(g),
                    fail_high_mode=fail_hi, fail_low_mode=fail_lo,
                    r_lo=None if lo is None else lo / int0,
                    r_hi=None if hi is None else hi / int0)

    # ---------------- 2. 逐点明细 -----------------------------------------
    print()
    print('| 模型 | 臂 | n | K | r | 判读 | counting | events | 稳态窗 | S2翻转 | '
          '边界裁剪 | 最小时序裕度 h | S2承诺 | 门剂量/周期 | 同刻 I2·RDF2 峰值 |')
    print('|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|')
    for kind in KIND_ORDER:
        for arm in ARM_ORDER[kind]:
            for n in (2.0, 3.0):
                g = sorted((r for r in rows
                            if r['kind'] == kind and r['arm'] == arm and r['n'] == n),
                           key=lambda r: r['K'])
                for r in g:
                    mech = r['mechanism']
                    dose = (mech['per_cycle'][0]['g1_dose']
                            if mech.get('per_cycle') else None)
                    marg = r['global_min_timing_margin_h']
                    print(f"| {kind} | {arm} | {n:g} | {r['K']:g} | {r['r']:.3f} | "
                          f"{classify(r)} | {int(r['counting_passed'])} | "
                          f"{int(r['event_causality_passed'])} | {r['steady_reads']} | "
                          f"{r['s2_crossings']} | {r['boundary_clips']} | "
                          f"{'-' if marg is None else f'{marg:.4f}'} | "
                          f"{r['minimum_commitment']} | "
                          f"{'-' if dose is None else f'{dose:.5f}'} | "
                          f"{mech['int2_times_rdf2_peak']:.4f} |")

    # ---------------- 3. 陡度 / 自馈 / 跨模型 -------------------------------
    print()
    print('### n 的效应（同一模型同一臂：成功区上界随 n 移动多少）')
    for kind in KIND_ORDER:
        for arm in ARM_ORDER[kind]:
            a = out['groups'].get(f'{kind}|{arm}|n=2')
            b = out['groups'].get(f'{kind}|{arm}|n=3')
            if not a or not b or a['pass_hi'] is None or b['pass_hi'] is None:
                print(f'- {kind}|{arm}: 有一侧没有成功点，无法比较')
                continue
            print(f"- {kind}|{arm}: 上界 K {a['pass_hi']:g}(n=2) → {b['pass_hi']:g}(n=3)，"
                  f"r {a['r_hi']:.3f} → {b['r_hi']:.3f}，"
                  f"相对移动 {100*(b['r_hi']-a['r_hi'])/a['r_hi']:+.1f}%")

    print()
    print('### 自馈开/关（HZZ）在同一个 K、n 下是否给出相同判读')
    diff = same = 0
    for n in (2.0, 3.0):
        for K in sorted({r['K'] for r in rows if r['kind'] == 'hzz'}):
            off = next((r for r in rows if r['kind'] == 'hzz' and r['arm'] == 'off'
                        and r['n'] == n and r['K'] == K), None)
            on = next((r for r in rows if r['kind'] == 'hzz' and r['arm'] == 'on'
                       and r['n'] == n and r['K'] == K), None)
            if not off or not on:
                continue
            if off['certified_v1'] == on['certified_v1']:
                same += 1
            else:
                diff += 1
                print(f"- n={n:g} K={K:g} (r={K/int0:.3f}): "
                      f"off={classify(off)} vs on={classify(on)}")
    print(f'- 判读相同 {same} 组，不同 {diff} 组')
    out['autoregulation_effect'] = dict(same=same, differ=diff)

    print()
    print('### 跨模型：同一个 r 上 HZH 与 HZZ 是否一致')
    for n in (2.0, 3.0):
        for K in sorted({r['K'] for r in rows if r['kind'] == 'hzh'}):
            h = next((r for r in rows if r['kind'] == 'hzh' and r['n'] == n
                      and r['K'] == K), None)
            z = next((r for r in rows if r['kind'] == 'hzz' and r['arm'] == 'on'
                      and r['n'] == n and r['K'] == K), None)
            if not h or not z:
                continue
            if h['certified_v1'] != z['certified_v1']:
                print(f"- n={n:g} r={K/int0:.3f} (K={K:g}): "
                      f"HZH={classify(h)} vs HZZ(on)={classify(z)}")

    (d / 'summary.json').write_text(
        json.dumps(out, ensure_ascii=False, indent=2, default=str), encoding='utf-8')
    print(f'\nwrote {d/"summary.json"}')


if __name__ == '__main__':
    main()
