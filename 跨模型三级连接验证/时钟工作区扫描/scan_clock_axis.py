"""第一组：时钟门本身 —— 成功区有多宽；差异来自阈值位置还是门的陡度。

轴向设计（一个能省掉一半工作的结论）
------------------------------------
HZZ 的 `NAMES = HZH_NAMES[:17] + 6 态末级`，且末级不反馈（prefix RHS gap 实测 0）。
⇒ **两套末级收到的 `Int0(t)` 是逐点相同的同一个信号**。
⇒ 用同一个 `K` 网格扫两套，就已经是在同一个无量纲比值 `r = K / Int0_peak` 上扫，
  不需要再"补测两边共同的 K/Int0_peak 比值"。
⇒ 但两者时钟 Hill 指数不同（HZH n=2、HZZ n=3，且 HZZ 是硬编码的 3），
  所以 `K` 的绝对值**不能直接排名**；本脚本同时报 `K` 与 `r`，`r` 才是可迁移的量。

网格
----
K ∈ {0.075, 0.15, 0.20, 0.22, 0.24, 0.26, 0.28, 0.30, 0.40, 0.60}
n ∈ {2, 3}
  * HZH × 1 臂（A1 自馈与成熟表达链按已认证工作点固定）
  * HZZ × 2 臂（自馈**关** = 源文件真实行为；自馈**开**）
合计 10 × 2 × 3 = 60 次，约 35 min。
K=0.30/n=2 是 HZH 冻结点；K=0.40/n=3 是 HZZ 已发布失败点 —— 两者兼作复现守卫。

守卫（任一不过则整轮作废，不产出任何结论）
------------------------------------------
1. `parity`    本文件的 verdict 必须与 `verify_hzh.analyse` 在 HZH 上**逐字段一致**
2. `structural` `model_hzz.structural_checks(prefix)`（40 随机探针）
3. `reproduction` HZH 冻结点 / HZZ 已发布失败点都必须逐字复现
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import common as C                                                  # noqa: E402

KEYS = (0.075, 0.15, 0.20, 0.22, 0.24, 0.26, 0.28, 0.30, 0.40, 0.60)
NS = (2.0, 3.0)
ARMS = {'hzh': (False,), 'hzz': (False, True)}
# HZH 没有自馈开关，标成 base —— 不写成 off，免得与 hzz|off 混为一谈
ARM_LABEL = {'hzh': {False: 'base'}, 'hzz': {False: 'off', True: 'on'}}


def run_guards(han, out):
    g = {}
    g['structural'] = C.structural_guard(C.HZHModel('han', han))

    # 复现守卫 1：HZH 冻结点（K=0.30, n=2）。这条 300 h 轨迹同时用作 parity 的输入。
    m, meta = C.build('hzh', han)
    t, y = C.integrate(m)
    v = C.verdict(t, C.sig_of(m, y, 'hzh'))
    g['parity'] = C.parity_guard(han, traj=(t, y))
    got = dict(certified_v1=bool(v['certified_v1']),
               steady=v['steady']['sequence'],
               stage0=(v['events']['bit0_to_bit1']['reverse_events'],
                       v['events']['bit0_to_bit1']['gate_events'],
                       v['events']['bit0_to_bit1']['flip_events']),
               stage1=(v['events']['bit1_to_bit2']['reverse_events'],
                       v['events']['bit1_to_bit2']['gate_events'],
                       v['events']['bit1_to_bit2']['flip_events']),
               margin_h=v['global_min_timing_margin_h'])
    want = C.FROZEN_HZH
    ok = (got['certified_v1'] == want['certified_v1'] and got['steady'] == want['steady']
          and got['stage0'] == want['stage0'] and got['stage1'] == want['stage1']
          and got['margin_h'] is not None and abs(got['margin_h'] - want['margin_h']) < 1e-6)
    g['reproduction_hzh'] = dict(ok=bool(ok), want=want, got=got)
    int0_peak = float(np.max(C.sig_of(m, y, 'hzh')['int0']))
    if not ok:
        C.dump(out / 'guards.json', g)
        raise SystemExit('REPRODUCTION GUARD FAILED (HZH) -- see guards.json')

    # 复现守卫 2：HZZ 已发布失败点（K=0.40, n=3, 两臂）
    hzz = {}
    for arm in (False, True):
        mm, _ = C.build('hzz', han, clock_K=0.40, clock_n=3.0, autoregulation=arm)
        tt, yy = C.integrate(mm)
        vv = C.verdict(tt, C.sig_of(mm, yy, 'hzz'))
        key = 'on' if arm else 'off'
        gotz = dict(steady=vv['steady']['sequence'],
                    stage1=(vv['events']['bit1_to_bit2']['reverse_events'],
                            vv['events']['bit1_to_bit2']['gate_events'],
                            vv['events']['bit1_to_bit2']['flip_events']),
                    s2_flips=len(vv['crossings']['S2']))
        okz = (gotz['steady'] == C.FROZEN_HZZ_FAIL[key]
               and gotz['stage1'] == C.FROZEN_HZZ_FAIL['stage1']
               and gotz['s2_flips'] == C.FROZEN_HZZ_FAIL['s2_flips'])
        hzz[key] = dict(ok=bool(okz), want=dict(steady=C.FROZEN_HZZ_FAIL[key],
                                                stage1=C.FROZEN_HZZ_FAIL['stage1'],
                                                s2_flips=C.FROZEN_HZZ_FAIL['s2_flips']),
                        got=gotz)
    g['reproduction_hzz'] = hzz
    C.dump(out / 'guards.json', g)
    if not all(x['ok'] for x in hzz.values()):
        raise SystemExit('REPRODUCTION GUARD FAILED (HZZ) -- see guards.json')
    return g, int0_peak


def build_row(kind, arm, arm_label, K, n, int0_peak, meta, v, pb, mech):
    """把一个点的全部读数折成一行。**单独成函数**，好让 `preflight.py` 用短轨迹
    把这条路径也走一遍——第一版把它内联在 main 里，结果 `kind` 与 `**meta` 撞键，
    在护栏跑完 4 min 之后才抛 `TypeError`，白白浪费一轮。

    `arm_label` 单独传：HZH **没有**自馈开关，把它也写成 `off` 会让人误以为
    `hzh|off` 与 `hzz|off` 是同一个开关的两个状态。HZH 一律记 `base`。
    """
    return dict(
        label=f'{kind}|{arm_label}|K={K:g}|n={n:g}',
        arm=arm_label, K=float(K), n=float(n),
        r=float(K) / int0_peak, **meta,
        certified_v1=bool(v['certified_v1']),
        counting_passed=bool(v['steady']['passed']),
        event_causality_passed=bool(all(e['passed'] for e in v['events'].values())),
        steady_sequence=v['steady']['sequence'],
        steady_reads=int(v['steady']['reads']),
        cold_reads=int(v['cold']['reads']),
        minimum_commitment=v['steady']['minimum_commitment'],
        boundary_clips=int(v['steady']['boundary_clips']),
        clock_peak_count=int(v['clock_peak_count']),
        s0=(v['events']['bit0_to_bit1']['reverse_events'],
            v['events']['bit0_to_bit1']['gate_events'],
            v['events']['bit0_to_bit1']['flip_events']),
        s1=(v['events']['bit1_to_bit2']['reverse_events'],
            v['events']['bit1_to_bit2']['gate_events'],
            v['events']['bit1_to_bit2']['flip_events']),
        s0_passed=bool(v['events']['bit0_to_bit1']['passed']),
        s1_passed=bool(v['events']['bit1_to_bit2']['passed']),
        s2_crossings=int(len(v['crossings']['S2'])),
        global_min_timing_margin_h=v['global_min_timing_margin_h'],
        per_bit={b: {k: pb[b][k] for k in pb[b]} for b in C.PB_BITS} if pb else {},
        run_timing={b: dict(setup_h=v['bit_margins'][b]['min_setup_h'],
                            hold_h=v['bit_margins'][b]['min_hold_h']) for b in v['bit_margins']},
        mechanism=mech,
        signal_ranges=v['signal_ranges'])


def fmt_row(i, total, row, pb, mech):
    per = {b: (pb[b]['commitment'], pb[b]['band_margin']) for b in C.PB_BITS} if pb else {}
    dose = mech['per_cycle'][0]['g1_dose'] if mech['per_cycle'] else float('nan')
    return (f'[{i}/{total}] {row["label"]:26s} full={int(row["certified_v1"])} '
            f'cnt={int(row["counting_passed"])} ev={int(row["event_causality_passed"])} '
            f'seq={row["steady_sequence"][:19]} '
            f'cmt={row["minimum_commitment"]} '
            f'gateDose={dose:.4f} '
            f'i2r2={mech["int2_times_rdf2_peak"]:.4f} '
            + ' '.join(f'{b}:{v2[0]:.2f}/{v2[1]:.3f}' for b, v2 in per.items()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--hours', type=float, default=C.HOURS)
    ap.add_argument('--quick', action='store_true', help='只跑 K=0.30/0.40 作冒烟')
    args = ap.parse_args()
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8')

    out = HERE / 'results' / f'clock_axis_{C.datetime.now().strftime("%Y%m%d_%H%M%S_%f")}'
    out.mkdir(parents=True, exist_ok=True)
    C.dump(out / 'status.json', dict(status='RUNNING'))

    print('building Han upstream once', flush=True)
    han = C.HanInput(args.hours)
    print(f'HAN sha = {C.sha256(Path(C.HAN_PATH))}  (expected {C.HAN_SHA})', flush=True)

    guards, int0_peak = run_guards(han, out)
    print(f'GUARDS OK  parity={guards["parity"]["ok"]}  '
          f'prefix_gap={guards["structural"]["prefix_rhs_gap"]}  '
          f'Int0_peak={int0_peak:.5f}  -> r = K / {int0_peak:.5f}', flush=True)

    keys = KEYS
    if args.quick:
        keys = tuple(k for k in KEYS if k in (0.30, 0.40))

    run = C.Run(out.parent, 'rows', explicit=out)
    total = len(keys) * len(NS) * (len(ARMS['hzh']) + len(ARMS['hzz']))
    i = 0
    for kind in ('hzh', 'hzz'):
        for K in keys:
            for n in NS:
                for arm in ARMS[kind]:
                    i += 1
                    m, meta = C.build(kind, han, clock_K=K, clock_n=n, autoregulation=arm)
                    t, y = C.integrate(m, hours=args.hours)
                    sig = C.sig_of(m, y, kind)
                    v = C.verdict(t, sig)
                    pb = C.pb_per_bit(sig, v['read_windows'], t)
                    mech = C.mechanism(t, sig, v['read_windows'])
                    row = build_row(kind, arm, ARM_LABEL[kind][arm], K, n, int0_peak,
                                    meta, v, pb, mech)
                    run.add(row)
                    print(fmt_row(i, total, row, pb, mech), flush=True)
                    C.dump(out / 'partial.json',
                           dict(hours=args.hours, int0_peak=int0_peak, rows=run.rows))

    final = run.save(dict(hours=args.hours, int0_peak=int0_peak, grid=dict(K=keys, n=NS),
                          arms=ARMS, guards=guards,
                          note=('same Int0(t) in both models (shared 17-state prefix, no feedback); '
                                'quote r = K / Int0_peak, not K; verdict dict is transcribe-identical '
                                'to verify_hzh.analyse (parity-guarded) so the adapter disclaimer lives '
                                'here and not in scope')))
    C.dump(out / 'status.json', dict(status='COMPLETED', rows=len(run.rows), final=str(final)))
    print('OUTPUT', final, flush=True)


if __name__ == '__main__':
    main()
