r"""Per-item acceptance breakdown of the 21-point static tolerance scan.

Why this file exists
--------------------
The pre-audit verdict carries ONE boolean per point, `certified = steady AND one AND
order AND alt AND contrast` (verify_threebit51.py:107), and the review asked for the
components to be shown side by side instead.  A single boolean hides the distinction
that matters most here: a point can be counting perfectly and still be uncertified
because the off/on gate-peak ratio - the arm H4 already demoted as a leak proxy -
crossed its limit.

This script re-tabulates the ALREADY SAVED scan.  It does not integrate anything, it
does not define a new threshold, and it does not rename anything:

  * `certified`            copied verbatim from `preaudit_stochasticity_all.csv` and
                           cross-checked against the verdict's `certification_table`;
  * `counting_and_causality_passed`  NEW, and defined exactly as
                           arm_steady AND arm_one AND arm_order AND arm_alt.
                           It is deliberately NOT called `certified`, and no leak
                           threshold is attached to it: H4's paired `L_symmetric` is
                           reported as a number, never as a pass/fail here.
  * the off/on threshold    is READ from the artefact (`off_on_gate_peak_threshold`),
                           never re-invented.

Sub-arm definitions (from preaudit_stochasticity.py:263-267, which calls
verify_threebit51.analyse_threebit):
  steady   : verify_threebit51.counter_verdict(reads, 8).passed
             = >=16 windows left after dropping one full mod-8 period, every label
               defined, every consecutive increment == +1 (mod 8)
  one      : exactly one carry gate and exactly one bit-2 flip per late reverse event,
             and no gate / crossing left unassigned
  order    : both the gate and the flip start after the reverse event starts
  alt      : bit-2 flip directions alternate
  contrast : off/on gate-peak ratio <= 0.10

Run:  & 'D:\aconade\python.exe' .\make_acceptance_breakdown.py
Out:  preaudit_acceptance_breakdown.json / .csv / .md
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
WIKI = HERE.parents[1]
STEM = 'preaudit_acceptance_breakdown'
SCAN = HERE / 'preaudit_stochasticity_all.csv'
VERDICT = HERE / 'preaudit_stochasticity_verdict.json'

HOURS = 600.0
SAMPLE_MIN = 2.0
TAIL_CUT = 0.05
MIN_STEADY_READS = 16
ARMS = ('steady', 'one', 'order', 'alt', 'contrast')
CAUSALITY_ARMS = ('one', 'order', 'alt')


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for blk in iter(lambda: fh.read(1 << 20), b''):
            h.update(blk)
    return h.hexdigest().upper()


def factor_key(table: dict, factor: float) -> str:
    """The verdict uses stringified factors ('1.0', not '1'); match by value."""
    for k in table:
        if abs(float(k) - float(factor)) < 1e-12:
            return k
    raise KeyError(f'factor {factor} not present in {sorted(table)}')


def main() -> int:
    try:
        sys.stdout.reconfigure(encoding='utf-8', errors='replace')
    except Exception:                                                   # noqa: BLE001
        pass
    scan = pd.read_csv(SCAN)
    verdict = json.loads(VERDICT.read_text(encoding='utf-8'))
    problems: list[str] = []

    # ---------------------------------------------------------------- build rows
    rows = []
    for r in scan.itertuples():
        arms = {a: bool(getattr(r, f'arm_{a}')) for a in ARMS}
        ccp = all(arms[a] for a in ('steady',) + CAUSALITY_ARMS)
        measured = dict(
            point_id=r.point_id, knob=r.knob, factor=float(r.factor),
            realised_factor=float(r.realised_factor),
            n_A1_gate_effective=float(r.n_A1_gate_effective),
            # 1. historical certification, ORIGINAL definition, kept as-is
            certified=bool(r.certified),
            certified_definition=('steady AND one AND order AND alt AND contrast '
                                  '(verify_threebit51.py:107)'),
            # 2. digital counting
            counting_sequence=r.sequence,
            counting_unlabelled_windows=int(r.unlabelled),
            counting_unlabelled_scope='over ALL 56 read windows (includes the 8 dropped '
                                      'for the steady verdict)',
            counting_steady_reads=int(r.steady_reads),
            counting_steady_reads_sufficient=bool(r.steady_reads_sufficient),
            counting_mod8_increments=arms['steady'],
            counting_cold_reads=int(r.cold_reads),
            # 3. carry causality
            causality_one_gate_and_flip_per_reverse=arms['one'],
            causality_order_after_reverse_start=arms['order'],
            causality_directions_alternate=arms['alt'],
            causality_unassigned_gates=int(r.unassigned_gates),
            causality_unassigned_crossings=int(r.unassigned_crossings),
            causality_crossings=int(r.crossings),
            causality_gates=int(r.gates),
            causality_reverse_events=int(r.reverse_events),
            # 4. gate contrast (raw ratio, no re-thresholding)
            contrast_off_on_ratio=float(r.off_on_gate_peak_ratio),
            contrast_threshold=float(r.off_on_gate_peak_threshold),
            contrast_passed=arms['contrast'],
            # 5. paired leak - a NUMBER, deliberately not a pass/fail
            leak_symmetric_median=float(r.leak_symmetric),
            leak_n_pairs=int(r.leak_pairs),
            leak_tail_cut='5 %',
            leak_output_grid=f'{SAMPLE_MIN:g} min, {HOURS:g} h',
            # 6. timing and read windows
            timing_setup_min_h=float(r.setup_min_h),
            timing_hold_min_h=float(r.hold_min_h),
            timing_min_commitment=float(r.min_commitment),
            timing_gate_peak_median=float(r.gate_peak_median),
            timing_gate_dose_median=float(r.gate_dose_median),
            timing_gate_width_median_h=float(r.gate_width_median),
            # the derived label
            counting_and_causality_passed=ccp,
            counting_and_causality_definition=('counting_mod8_increments AND '
                                               'causality_one_gate_and_flip_per_reverse '
                                               'AND causality_order_after_reverse_start '
                                               'AND causality_directions_alternate'),
        )
        rows.append(measured)

    # ---------------------------------------------------------------- cross-checks
    # `leak_table` in the verdict is a DISPLAY table rounded to 6 dp; the shard CSV keeps
    # full precision.  Compare at the verdict's own displayed precision and report the
    # residual, rather than pretending the two are bit-identical.
    leak_display_dp = 6
    leak_max_residual = 0.0
    for row in rows:
        fk = factor_key(verdict['certification_table'], row['factor'])
        knob = row['knob']
        v_cert = bool(verdict['certification_table'][fk][knob])
        if v_cert != row['certified']:
            problems.append(f"{row['point_id']}: certified {row['certified']} != verdict "
                            f"{v_cert}")
        v_leak = float(verdict['leak_table'][fk][knob])
        residual = abs(round(row['leak_symmetric_median'], leak_display_dp) - v_leak)
        leak_max_residual = max(leak_max_residual, residual)
        if residual > 10 ** (-leak_display_dp):
            problems.append(f"{row['point_id']}: leak {row['leak_symmetric_median']} != "
                            f"verdict {v_leak} even at {leak_display_dp} dp")
        v_ratio = float(verdict['off_on_gate_peak_ratio_table'][fk][knob])
        if abs(v_ratio - row['contrast_off_on_ratio']) > 1e-12:
            problems.append(f"{row['point_id']}: off/on {row['contrast_off_on_ratio']} != "
                            f"verdict {v_ratio}")

    if len({r['contrast_threshold'] for r in rows}) != 1:
        problems.append('the off/on threshold is not constant across points; this script '
                        'must not pick one')
    threshold = rows[0]['contrast_threshold']

    def arm_ok(row: dict, arm: str) -> bool:
        return {
            'steady': row['counting_mod8_increments'],
            'one': row['causality_one_gate_and_flip_per_reverse'],
            'order': row['causality_order_after_reverse_start'],
            'alt': row['causality_directions_alternate'],
            'contrast': row['contrast_passed'],
        }[arm]

    by_arm = {a: sum(1 for r in rows if not arm_ok(r, a)) for a in ARMS}
    contrast_only = sorted(r['point_id'] for r in rows
                           if not r['contrast_passed']
                           and all(arm_ok(r, a) for a in ('steady', 'one', 'order', 'alt')))
    perfect_but_uncertified = sorted(
        r['point_id'] for r in rows
        if r['counting_and_causality_passed'] and not r['certified'])
    # The v1 report asserted "3 points have a perfect digital readout yet are
    # uncertified".  That count is NOT reproducible from the artefacts under any single
    # crisp definition, and saying so is more useful than quietly picking one.  Two
    # defensible definitions both give 4, and they are DIFFERENT SETS.
    scope_a = perfect_but_uncertified
    scope_b = sorted(
        r['point_id'] for r in rows
        if r['counting_mod8_increments'] and r['counting_unlabelled_windows'] == 0
        and not r['certified'])
    v1_recheck = dict(
        v1_claim='3 个点的数字读出完美却未认证',
        reproducible_from_artefacts=False,
        definitions=[
            dict(name='counting_and_causality_passed AND NOT certified '
                      '(steady-window scope)',
                 members=scope_a, n=len(scope_a)),
            dict(name='steady mod-8 increments AND zero unlabelled windows over all 56 '
                      'AND NOT certified',
                 members=scope_b, n=len(scope_b)),
        ],
        note=('both crisp definitions give 4, not 3, and they are different sets: the '
              'first contains conc_scale|f=1.1 (2 unlabelled windows in the dropped cold '
              'prefix) and the second contains conc_scale|f=0.8 (perfect sequence but '
              'fails the three causality arms). Over ALL 56 windows there are 4 points '
              'whose sequence is exact AND fully labelled, so v1\'s 3 cannot be '
              'reconstructed; the reproducible statements are the two above.'))
    if by_arm != verdict['certification_arms']['failures_by_arm']:
        problems.append(f'failures by arm {by_arm} != verdict '
                        f"{verdict['certification_arms']['failures_by_arm']}")
    v_contrast_only = verdict['certification_arms']['contrast_only_failures']
    v_contrast_only = sorted(v_contrast_only) if isinstance(v_contrast_only, list) \
        else int(v_contrast_only)
    if isinstance(v_contrast_only, list):
        if contrast_only != v_contrast_only:
            problems.append(f'contrast-only {contrast_only} != verdict {v_contrast_only}')
    elif len(contrast_only) != v_contrast_only:
        problems.append(f'contrast-only count {len(contrast_only)} != verdict '
                        f'{v_contrast_only}')

    # ------------------------------------------------------- discrete passing sets
    def contiguous(factors: list, tested: list) -> bool:
        """True only if every TESTED factor between min and max is in the set."""
        if not factors:
            return True
        lo, hi = min(factors), max(factors)
        inside = [f for f in tested if lo - 1e-12 <= f <= hi + 1e-12]
        return len(inside) == len(factors)

    passing_sets = {}
    for knob in sorted({r['knob'] for r in rows}):
        tested = sorted(r['factor'] for r in rows if r['knob'] == knob)
        mine = [r for r in rows if r['knob'] == knob]
        cert = sorted(r['factor'] for r in mine if r['certified'])
        ccp = sorted(r['factor'] for r in mine if r['counting_and_causality_passed'])
        failed_by_arm = {a: sorted(r['factor'] for r in mine if not arm_ok(r, a))
                         for a in ARMS}
        failed_by_arm = {a: f for a, f in failed_by_arm.items() if f}
        passing_sets[knob] = dict(
            tested_factors=tested,
            certified_factors=cert,
            counting_and_causality_factors=ccp,
            certified_is_contiguous=contiguous(cert, tested),
            counting_and_causality_is_contiguous=contiguous(ccp, tested),
            failed_factors_by_arm=failed_by_arm,
        )
    for knob, s in passing_sets.items():
        why = '；'.join(f'{a}: {f}' for a, f in s['failed_factors_by_arm'].items()) \
            or '无失败'
        s['reading'] = (
            f"已测因子 {s['tested_factors']}。历史谓词通过的是 "
            f"{s['certified_factors']}；计数+因果通过的是 "
            f"{s['counting_and_causality_factors']}。失败原因按因子分布：{why}。"
            f"**两组都只是「已测网格点」的集合** —— 未采样的因子从未运行，"
            f"因此都不能当作连续可靠范围引用。")

    doc = dict(
        kind='preaudit_acceptance_breakdown',
        source_scan=str(SCAN.relative_to(WIKI)),
        source_scan_sha256=sha256(SCAN),
        source_verdict=str(VERDICT.relative_to(WIKI)),
        source_verdict_sha256=sha256(VERDICT),
        no_reintegration=('purely a re-tabulation of the saved 21-point scan; no model '
                          'was built and no trajectory integrated'),
        definition_note=('counting_and_causality_passed is a NEW label and is NOT a '
                         'renaming of certified; certified keeps its original '
                         'five-arm definition. No leak threshold is introduced here: '
                         'H4 paired L_symmetric is reported as a number only.'),
        scan_settings=dict(hours=HOURS, sample_min=SAMPLE_MIN,
                           tail_cut=TAIL_CUT, min_steady_reads=MIN_STEADY_READS,
                           factor_grid=sorted({r['factor'] for r in rows})),
        off_on_gate_peak_threshold=threshold,
        leak_verdict_display_precision_dp=leak_display_dp,
        leak_verdict_max_rounding_residual=leak_max_residual,
        leak_rounding_note=('the verdict\'s leak_table is a DISPLAY table rounded to '
                            f'{leak_display_dp} dp; the shard CSV keeps full precision. '
                            'The cross-check compares at the displayed precision.'),
        totals=dict(points=len(rows),
                    certified=sum(r['certified'] for r in rows),
                    counting_and_causality_passed=sum(
                        r['counting_and_causality_passed'] for r in rows)),
        failures_by_arm=by_arm,
        contrast_only_failures=contrast_only,
        perfect_readout_but_uncertified=perfect_but_uncertified,
        v1_claim_recheck=v1_recheck,
        discrete_passing_points=passing_sets,
        points=rows,
        cross_check_problems=problems,
        model_sha256=verdict['model_sha256'],
        source_sha256=verdict['source_sha256'],
    )
    (HERE / f'{STEM}.json').write_text(
        json.dumps(doc, ensure_ascii=False, indent=2), encoding='utf-8')
    pd.DataFrame(rows).to_csv(HERE / f'{STEM}.csv', index=False, encoding='utf-8')

    # ------------------------------------------------------------------ markdown
    L = []
    L.append('# 静态容差扫描：21 点分项验收表\n')
    L.append(f'来源：`{doc["source_scan"]}`（sha256 `{doc["source_scan_sha256"]}`）'
             f'与 `{doc["source_verdict"]}`。'
             f'**纯后处理**：未构建模型、未积分任何轨迹。\n')
    L.append(f'- 设置：{HOURS:g} h、输出网格 {SAMPLE_MIN:g} min、配对泄漏切尾 '
             f'{TAIL_CUT:.0%}、稳态读数下限 {MIN_STEADY_READS}\n'
             f'- 因子网格：{doc["scan_settings"]["factor_grid"]}\n'
             f'- 门对比度限值：**{threshold}**（从产物读出，未重新设定）\n')
    L.append('## 0. 两个标签的区别（重要）\n')
    L.append('- `certified`：**历史值，原定义不变** —— '
             '`steady ∧ one ∧ order ∧ alt ∧ contrast`（`verify_threebit51.py:107`）。\n'
             '- `counting_and_causality_passed`：**本轮新增**，定义恰为 '
             '`steady ∧ one ∧ order ∧ alt`。**它不是 `certified` 的改名**，'
             '也不附带任何泄漏阈值 —— H4 的配对 `L_symmetric` 在这里只作为**数字**列出。\n')
    L.append(f'## 1. 总计\n')
    L.append(f'| 项 | 值 |\n|---|---|\n'
             f'| 点数 | {doc["totals"]["points"]} |\n'
             f'| 历史 `certified` 通过 | **{doc["totals"]["certified"]}** / '
             f'{doc["totals"]["points"]} |\n'
             f'| `counting_and_causality_passed` | **'
             f'{doc["totals"]["counting_and_causality_passed"]}** / '
             f'{doc["totals"]["points"]} |\n'
             f'| 逐臂失败次数 | {doc["failures_by_arm"]} |\n'
             f'| 仅 contrast 失败 | {len(contrast_only)} 个：{contrast_only} |\n'
             f'| 数字与因果全过、但未认证 | {len(perfect_but_uncertified)} 个：'
             f'{perfect_but_uncertified} |\n')
    L.append(f'> **v1 报告里「3 个点的数字读出完美却未认证」不可复现**：两个清晰定义都给出 4，'
             f'而且成员不同 ——\n'
             f'> - A = `counting_and_causality_passed ∧ ¬certified` → `{scope_a}`\n'
             f'> - B = `稳态模 8 递增 ∧ 全 56 窗零未标注 ∧ ¬certified` → `{scope_b}`\n'
             f'> 差别在 `conc_scale|f=1.1`（2 个未标注窗落在被丢弃的起始段，故只进 A）'
             f'与 `conc_scale|f=0.8`（码串完美但三条因果臂全败，故只进 B）。'
             f'本表只报这两个可复现的说法，不去凑 v1 的 3。\n')
    L.append(f'> 泄漏列在此处只作**数字**列出：`preaudit_stochasticity_verdict.json` 的 '
             f'`leak_table` 是**显示用表、四舍五入到 {leak_display_dp} 位**，'
             f'分片 CSV 才是全精度；交叉核对按显示精度比，'
             f'最大残差 `{leak_max_residual:.2e}`。\n')
    L.append('## 2. 逐点明细\n')
    hdr = ('| 点 | `certified` | 计数(模8递增) | 未标注窗 | 稳态码串 | 一一对应 | 因果顺序 | '
           '方向交替 | off/on | 超限? | `L_symmetric`(5%) | setup/hold (h) | commitment | '
           '**计数+因果** |')
    L.append(hdr)
    L.append('|' + '---|' * 14)
    for r in rows:
        leak = ('n/a（无有效配对）' if not np.isfinite(r['leak_symmetric_median'])
                else f"{r['leak_symmetric_median']:.6f}")
        L.append(
            f"| `{r['point_id']}` | "
            f"{'✅' if r['certified'] else '❌'} | "
            f"{'✅' if r['counting_mod8_increments'] else '❌'} | "
            f"{r['counting_unlabelled_windows']} | "
            f"`{r['counting_sequence']}` | "
            f"{'✅' if r['causality_one_gate_and_flip_per_reverse'] else '❌'} | "
            f"{'✅' if r['causality_order_after_reverse_start'] else '❌'} | "
            f"{'✅' if r['causality_directions_alternate'] else '❌'} | "
            f"{r['contrast_off_on_ratio']:.4f} | "
            f"{'**是**' if not r['contrast_passed'] else '否'} | "
            f"{leak} | "
            f"{r['timing_setup_min_h']:.2f} / {r['timing_hold_min_h']:.2f} | "
            f"{r['timing_min_commitment']:.2f} | "
            f"{'✅' if r['counting_and_causality_passed'] else '❌'} |")
    L.append('\n注：`未标注窗` 是对**全部 56 个读窗**计的；`计数(模8递增)` 用的是丢掉一个完整'
             '模 8 周期后的 48 个稳态窗。两者口径不同，因此 `conc_scale|f=1.1` 会出现'
             '「未标注 2 个、但稳态模 8 递增成立」——那 2 个落在被丢弃的起始段。\n')
    L.append('## 3. 离散通过点（**不是**连续可靠范围）\n')
    for knob, s in sorted(passing_sets.items()):
        L.append(f'### `{knob}`\n')
        L.append(f'- 已测因子：`{s["tested_factors"]}`\n'
                 f'- 历史 `certified` 通过：`{s["certified_factors"]}`'
                 f'（连续={s["certified_is_contiguous"]}）\n'
                 f'- `counting_and_causality_passed` 通过：'
                 f'`{s["counting_and_causality_factors"]}`'
                 f'（连续={s["counting_and_causality_is_contiguous"]}）\n')
        L.append(f'{s["reading"]}\n')
    if problems:
        L.append('## 4. 交叉核对问题\n')
        L += [f'- {p}\n' for p in problems]
    else:
        L.append('## 4. 交叉核对\n')
        L.append('与 `preaudit_stochasticity_verdict.json` 的 `certification_table`、'
                 '`leak_table`、`off_on_gate_peak_ratio_table`、`failures_by_arm` '
                 '**逐点一致，0 个问题**。\n')
    (HERE / f'{STEM}.md').write_text('\n'.join(L), encoding='utf-8')

    print(f'points                       : {doc["totals"]["points"]}')
    print(f'certified (historical)       : {doc["totals"]["certified"]}')
    print(f'counting_and_causality_pass  : {doc["totals"]["counting_and_causality_passed"]}')
    print(f'failures by arm              : {by_arm}')
    print(f'contrast-only failures       : {contrast_only}')
    print(f'perfect readout, uncertified : {perfect_but_uncertified}')
    print(f'cross-check problems         : {len(problems)}')
    for p in problems:
        print('  PROBLEM:', p)
    print(f'wrote {STEM}.json / .csv / .md')
    return 1 if problems else 0


if __name__ == '__main__':
    sys.exit(main())
