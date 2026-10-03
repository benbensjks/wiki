"""共用适配器：把 HZH 与 HZZ 两套末级送进**同一份判读代码**，外加护栏与机制量。

设计原则（三条，都是为了"两套模型并列可比"这件事本身成立）
------------------------------------------------------------------
1. **判读逻辑不复制。** `read_windows / read_verdict / margins / associate / segments /
   crossings` 全部直接调用 `verify_hzh` 的纯函数（它们只吃 `(t, sig)`）。
   本文件只转写 `verify_hzh.analyse` 的**编排**（源：`verify_hzh.py:162-191`）。
2. **两套模型唯一的差别是 `sig` 从哪来**：HZH 用 `verify_hzh.signals(y, m.z, m.tail)`
   （硬编码 34 态索引），HZZ 用 `model_hzz.HZZModel.signals(y)`（23 态）。
3. **`parity_guard()` 用同一条 HZH 轨迹断言本文件的 `verdict(t, sig)` 与
   `verify_hzh.analyse(t, y, m.z, m.tail)` 输出完全一致**——防止第 1 条的转写漂移。
   这道守卫不过，任何扫描结果都不得引用。

本文件只读既有模型与验证器，不修改任何既有文件。
"""
from __future__ import annotations

import csv
import hashlib
import json
import sys
import time
from dataclasses import asdict, replace
from datetime import datetime
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

WIKI = Path(r'C:\Users\18633\Desktop\wiki')
XM = WIKI / '跨模型三级连接验证'
HZH_DIR = XM / 'hby_zmh_hby'
HZZ_DIR = XM / 'hby_zmh_zmh'
CERT_DIR = HZH_DIR / 'certification'
MIX_DIR = XM / '混合参数扰动'
for _d in (XM, HZH_DIR, HZZ_DIR, CERT_DIR, MIX_DIR):
    if str(_d) not in sys.path:
        sys.path.insert(0, str(_d))

import verify_hzh as VH                                            # noqa: E402
from run_han_comparison import HAN_PATH, HAN_SHA, HanInput          # noqa: E402
from model_hzh import HZHModel, NAMES as HZH_NAMES                  # noqa: E402
from model_hzz import HZZModel, NAMES as HZZ_NAMES                  # noqa: E402
from hybrid_model import HbyConfig, HbyReceiver                     # noqa: E402
from perbit_margin import BITS as PB_BITS, per_bit as pb_per_bit     # noqa: E402

HOURS = 300.0
SAMPLE_MIN = 1.0
MAX_STEP_MIN = 1.0
RTOL = 2e-7
ATOL = 2e-9
COPIES_PER_UM_PER_FL = 602.214076

# 已发布的对照值（用于复现守卫；来源见各模型自己的 results/）
FROZEN_HZH = dict(certified_v1=True, steady='1234567012345670123',
                  stage0=(14, 14, 14), stage1=(7, 7, 7), margin_h=1.1784609018563117)
FROZEN_HZZ_FAIL = dict(off='xxx4567456745674567', on='5674567456745674567',
                       stage1=(7, 0, 1), s2_flips=1)

SOURCE_FILES = {
    'hby_zmh_hby/model_hzh.py': HZH_DIR / 'model_hzh.py',
    'hby_zmh_hby/certification/verify_hzh.py': CERT_DIR / 'verify_hzh.py',
    'hby_zmh_zmh/model_hzz.py': HZZ_DIR / 'model_hzz.py',
    'hybrid_model.py': XM / 'hybrid_model.py',
    'run_han_comparison.py': XM / 'run_han_comparison.py',
    'wiki任务/前馈三级级联.py': WIKI / 'wiki任务' / '前馈三级级联.py',
    'wiki任务/完整二级级联.py': WIKI / 'wiki任务' / '完整二级级联.py',
    'final_reconstruction/model.py': WIKI / 'final_reconstruction' / 'model.py',
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest().upper()


def dump(path: Path, obj) -> None:
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False),
                    encoding='utf-8')


def source_hashes() -> dict:
    return {k: sha256(v) for k, v in SOURCE_FILES.items() if v.exists()}


# --------------------------------------------------------------------- 构建
def build(kind: str, han, *, clock_K=None, clock_n=None, uM=None,
          alpha_int2=None, autoregulation=False, seed_model=None):
    """构造一套末级模型。`kind` = 'hzh' | 'hzz'。

    clock_K/clock_n/uM/alpha_int2 为 None 时沿用该模型自己的基线值。
    全部改动都作用在**本实例自己的**对象上（HbyConfig 用 replace 新建、
    receiver.p 是 load_zeng_table() 每次返回的新字典），不触碰全局表。
    """
    meta = dict(kind=kind, autoregulation=bool(autoregulation))
    if kind == 'hzh':
        m = seed_model if seed_model is not None else HZHModel('han', han)
        cfg = m.tail.c
        new_cfg = replace(cfg,
                          clock_K_au=cfg.clock_K_au if clock_K is None else float(clock_K),
                          clock_n=cfg.clock_n if clock_n is None else float(clock_n),
                          receiver_uM_per_au=cfg.receiver_uM_per_au if uM is None else float(uM))
        m.tail = HbyReceiver(new_cfg)
        m.copies_per_au = COPIES_PER_UM_PER_FL * m.tail.c.receiver_uM_per_au
        if alpha_int2 is not None:
            m.tail.p = dict(m.tail.p)
            seq = list(m.tail.p['alpha_Int'])
            seq[1] = float(alpha_int2)
            m.tail.p['alpha_Int'] = tuple(seq)
        meta.update(clock_K=float(m.tail.c.clock_K_au), clock_n=float(m.tail.c.clock_n),
                    clock_scale=float(m.tail.c.clock_scale),
                    uM_per_au=float(m.tail.c.receiver_uM_per_au),
                    alpha_int2=float(m.tail.p['alpha_Int'][1]),
                    n_A1_gate=float(m.tail.c.n_A1_gate))
    elif kind == 'hzz':
        prefix = seed_model if seed_model is not None else HZHModel('han', han)
        if uM is not None:
            cfg = replace(prefix.tail.c, receiver_uM_per_au=float(uM))
            prefix.tail = HbyReceiver(cfg)
            prefix.copies_per_au = COPIES_PER_UM_PER_FL * prefix.tail.c.receiver_uM_per_au
        m = HZZModel(prefix, autoregulation)
        if clock_K is not None:
            m.clock_K = float(clock_K)       # 实例属性遮蔽类属性，不改类
        if clock_n is not None:
            m.clock_n = float(clock_n)
        if alpha_int2 is not None:
            m.p = dict(m.p)
            m.p['alpha_Int2'] = float(alpha_int2)
        meta.update(clock_K=float(m.clock_K), clock_n=float(m.clock_n), clock_scale=1.0,
                    uM_per_au=float(prefix.tail.c.receiver_uM_per_au),
                    alpha_int2=float(m.p['alpha_Int2']),
                    n_A1_gate=None)
    else:
        raise ValueError(kind)
    meta['copies_per_au'] = float((m.prefix if kind == 'hzz' else m).copies_per_au)
    return m, meta


def integrate(model, hours=HOURS, sample_min=SAMPLE_MIN, max_step_min=MAX_STEP_MIN,
              rtol=RTOL, atol=ATOL):
    n = round(hours * 60 / sample_min)
    t = np.linspace(0, hours, n + 1)
    if len(t) < 2 or not np.all(np.diff(t) > 0):
        raise ValueError('bad sampling grid')
    sol = solve_ivp(model.rhs, (0, hours), model.initial_state(), t_eval=t,
                    method='DOP853', rtol=rtol, atol=atol, max_step=max_step_min / 60)
    if not sol.success or not np.all(np.isfinite(sol.y)):
        raise RuntimeError(sol.message)
    return sol.t, sol.y


# ------------------------------------------------------------- sig 与判读
def sig_of(model, y, kind):
    """两套模型唯一的差别就在这里。"""
    if kind == 'hzh':
        sig = VH.signals(y, model.z, model.tail)
        sig['int0'] = y[VH.IDX['b0_I']]
        sig['S2'] = y[VH.IDX['b2_S']]
        sig['_Int2'] = y[VH.IDX['b2_I']]
        sig['_RDF2'] = y[VH.IDX['b2_R']]
        sig['_u2'] = sig['g1'] * float(model.tail.p['alpha_Int'][1])
    else:
        sig = model.signals(y)
        sig['_Int2'] = y[20]
        sig['_RDF2'] = y[22]
        sig['_u2'] = sig['g1'] * float(model.p['alpha_Int2'])
    return sig


def verdict(t, sig):
    """`verify_hzh.analyse` 编排的忠实转写（源：verify_hzh.py:162-191）。

    与原函数的唯一差别：`sig` 由调用方提供，因此两套模型共用本函数。
    `parity_guard()` 断言它与原函数在 HZH 上逐字段一致。
    """
    peaks, reads = VH.read_windows(t, sig)
    cross = {b: VH.crossings(t, sig[b]) for b in ('S0', 'S1', 'S2')}
    events = {}
    for stage, revname, gatename, bit in ((0, 'J_rev0', 'g0', 'S1'), (1, 'J_rev1', 'g1', 'S2')):
        rev = VH.segments(t, sig[revname], VH.JREV_THRESHOLD_H)
        gates = VH.segments(t, sig[gatename], VH.GATE_THRESHOLD, min_dose=VH.GATE_MIN_DOSE_H)
        events[f'bit{stage}_to_bit{stage+1}'] = VH.associate(rev, gates, cross[bit], float(t[-1]))
    steady = VH.read_verdict(reads, VH.DROP)
    bit_margins = {b: VH.margins(reads, cross[b]) for b in cross}
    nonnull = [x for m in bit_margins.values()
               for x in (m['min_setup_h'], m['min_hold_h']) if x is not None]
    v = dict(version='HZH_MOD8_CAUSAL_V1', hours=float(t[-1]),
             rule=dict(window='Int0 peak-to-peak trough, total width 20% period',
                       bands=[VH.LOW, VH.HIGH], occupancy=VH.OCC, drop=VH.DROP,
                       min_steady_reads=VH.MIN_STEADY,
                       reverse_flux_threshold_per_h=VH.JREV_THRESHOLD_H,
                       gate_threshold=VH.GATE_THRESHOLD,
                       min_event_h=VH.PULSE_DURATION_H,
                       min_gate_dose_h=VH.GATE_MIN_DOSE_H,
                       association='nearest reverse-event peak interval; causal order separately checked',
                       initial_carry_events_dropped=2,
                       acceptance='steady mod8 + both event chains one-to-one, ordered, direction alternating'),
             cold=VH.read_verdict(reads, 0), steady=steady, clock_peak_count=len(peaks),
             read_windows=reads, crossings=cross, events=events, bit_margins=bit_margins,
             global_min_timing_margin_h=min(nonnull) if nonnull else None,
             signal_ranges={k: [float(np.min(x)), float(np.max(x))]
                            for k, x in sig.items() if not k.startswith('_')},
             certified_v1=bool(steady['passed'] and all(e['passed'] for e in events.values())),
             legacy_certified=None,
             # 逐字抄自 verify_hzh.analyse:190 —— 这一行曾经被我改成"shared adapter"
             # 的说明文字，`parity_guard` 当场报出 fields=['scope'] 不匹配。
             # 结论：判读字典里**一个字符都不许改**；本适配器的说明写在 run 的
             # `note` 字段里，不写进 verdict。
             scope='Deterministic 300 h HZH v1; no inheritance of 51-state gate-contrast predicate; no experimental reliability claim')
    return v


def parity_guard(han, traj=None, hours=40.0):
    """本文件的 `verdict(t, sig)` 必须与 `verify_hzh.analyse` 在 HZH 上逐字段一致。

    `traj=(t, y)` 可传入一条**已有的**轨迹（例如 300 h 冻结点），
    这样检查覆盖 56 个读窗与 14 次进位事件，比 40 h 冒烟强得多且不额外花算力。
    """
    m = HZHModel('han', han)
    t, y = traj if traj is not None else integrate(m, hours=hours)
    ref, _ = VH.analyse(t, y, m.z, m.tail)
    mine = verdict(t, sig_of(m, y, 'hzh'))
    if json.dumps(ref, sort_keys=True, default=str) != json.dumps(mine, sort_keys=True, default=str):
        diffs = [k for k in set(ref) | set(mine)
                 if json.dumps(ref.get(k), sort_keys=True, default=str)
                 != json.dumps(mine.get(k), sort_keys=True, default=str)]
        raise RuntimeError(f'PARITY FAILED on fields: {diffs}')
    return dict(ok=True, hours=float(t[-1]), fields=len(ref),
                read_windows=len(ref['read_windows']),
                events={k: (v['reverse_events'], v['gate_events'], v['flip_events'])
                        for k, v in ref['events'].items()})


# ------------------------------------------------------------------ 机制量
def mechanism(t, sig, cycles):
    """机制量：门的开启占比/剂量、逐进位 g1 剂量与驻留残余、
    实际成熟 Int2 波形、**同一时刻** Int2·RDF2、逐进位正反向通量积分。"""
    i2 = np.asarray(sig['_Int2'], dtype=float)
    r2 = np.asarray(sig['_RDF2'], dtype=float)
    prod = i2 * r2                       # 同刻乘积（不取峰值再相乘）
    clk = np.asarray(sig['clock'], dtype=float)
    g1 = np.asarray(sig['g1'], dtype=float)
    out = dict(
        clock_duty_gt_0p1=float(np.mean(clk > 0.1)),
        clock_duty_gt_0p5=float(np.mean(clk > 0.5)),
        clock_integral=float(np.trapezoid(clk, t)),
        clock_peak=float(clk.max()),
        int2_peak=float(i2.max()), rdf2_peak=float(r2.max()),
        int2_times_rdf2_peak=float(prod.max()),
        int2_times_rdf2_median=float(np.median(prod)),
        per_cycle=[])
    for c in cycles:
        m = (t >= c['cycle_start_h']) & (t <= c['cycle_end_h'])
        if m.sum() < 3:
            continue
        tt, gg = t[m], g1[m]
        k = max(1, int(0.2 * tt.size))
        out['per_cycle'].append(dict(
            cycle=c['cycle'],
            g1_dose=float(np.trapezoid(np.maximum(gg, 0.0), tt)),
            g1_peak=float(gg.max()),
            g1_far_dwell=float(np.median(gg[np.argsort(gg)[:k]])),
            clock_dose=float(np.trapezoid(clk[m], tt)),
            int2_dose=float(np.trapezoid(np.maximum(i2[m], 0.0), tt)),
            J_fwd2=float(np.trapezoid(np.maximum(np.asarray(sig['J_fwd2'])[m], 0.0), tt)),
            J_rev2=float(np.trapezoid(np.maximum(np.asarray(sig['J_rev2'])[m], 0.0), tt)),
            i2r2_same_time_peak=float(prod[m].max())))
    return out


def structural_guard(prefix):
    """结构护栏：直接用 model_hzz 自带的那一份（40 个随机探针，更严）。

    它同时校验：末级 RHS 与原文件 RHS 逐点一致、自馈开关只改 A1 那一项、
    前缀 RHS 对末级状态缩放完全免疫（gap 必须恰为 0）。三项任一不达标即抛错。
    """
    from model_hzz import structural_checks
    return structural_checks(prefix)


# --------------------------------------------------------------- 运行与落盘
class Run:
    def __init__(self, root: Path, tag: str, explicit: Path | None = None):
        if explicit is not None:
            self.dir = Path(explicit)
        else:
            self.dir = root / f'{tag}_{datetime.now().strftime("%Y%m%d_%H%M%S_%f")}'
        self.dir.mkdir(parents=True, exist_ok=True)
        self.rows = []

    def add(self, row):
        self.rows.append(row)

    def save(self, extra: dict) -> Path:
        dump(self.dir / 'all.json', dict(rows=self.rows, **extra))
        if self.rows:
            keys = sorted({k for r in self.rows for k in r})
            with (self.dir / 'all.csv').open('w', newline='', encoding='utf-8') as fh:
                w = csv.DictWriter(fh, fieldnames=keys, extrasaction='ignore')
                w.writeheader()
                for r in self.rows:
                    w.writerow({k: (json.dumps(v, ensure_ascii=False) if isinstance(v, (dict, list)) else v)
                                for k, v in r.items()})
        dump(self.dir / 'source_hashes.json', source_hashes())
        dump(self.dir / 'SHA256SUMS.json',
             {p.name: sha256(p) for p in sorted(self.dir.iterdir())
              if p.is_file() and p.name != 'SHA256SUMS.json'})
        dump(self.dir / 'status.json', dict(status='COMPLETED', rows=len(self.rows)))
        return self.dir
