# CURRENT_BASELINE — TEMPO 建模权威口径

> **用途**：Wiki 写作与对外汇报的唯一取数口径。任何数字若不在本文中，必须先落盘为产物再引用。
> **状态**：只读核对完成，未修改任何模型或结果文件。
> **建立日期**：2026-09-21
> **核对范围**：`threebit51_selected_v1.json`、`threebit51_provisional_n6.json`、`threebit51_leak_correction.json`、
> `strict_tolerance_tol_verdict.json`、`strict_tolerance_grid_verdict.json`、`pool_perturbations_all_verdict.json`、
> `eight_initial_verdict.json`、`FROZEN_34_STATE.json`、`selected_profile.json`、`postfreeze_verdict.json`、
> `absolute_scale.json`、`local_identifiability_verdict.json`、`projection_equivalence.json`、
> `scan_carry_pairing.py`、`check_absolute_scale.py`

---

## 0. 一句话权威口径

> **51 状态 n=6 确定性三级模型已冻结，并在规定参数、初态、读窗和有限扰动条件下通过模 8 计数验证。
> 后冻结附件增加了 30 次有效扰动验证，但不同扰动时相的证据强度应分别说明。
> 生长条件统一、浓度标定、独立门指数的生物实现、随机性及可辨识性仍属开放问题。**

### 0.1 三个状态必须分开陈述，不得互相顶替

| 状态 | 结论 |
|---|---|
| 模型冻结 | ✅ 已完成 |
| Wiki 章节定稿 | ❌ 未完成 |
| 同一细胞条件下的生物可行性确认 | ❌ 未完成 |

**禁止表述**：因模型冻结而把后两者写成已完成。
**禁止表述**：`冻结 n=6 完整电路占细胞翻译预算 5.33%，已证明表达可行`（见 §6）。

---

## 1. 冻结件（可直接引用）

### 1.1 两比特：`twobit34_selected_v1`

来源 `twobit34_results/FROZEN_34_STATE.json`，`status = "frozen lower-stage reference for three-bit development"`。

| 项 | 值 |
|---|---|
| `uM_per_au` | 5.75 |
| `carry0_mrna_half_life_min` | 2.0 |
| `carry0_maturation_half_life_min` | 32.5 |
| `add_growth` | false |

冻结规则（原文）：三级开发期间不得改动冻结源与参数档；任何 51 状态候选的前 34 投影状态必须保持动力学等价；未来两比特改动须发新版本与新哈希，不得覆盖本次冻结。

有效指标：`minimum_timing_margin_h = 2.9137227753952324`、`bit1_setup_h = 4.637884499952506`、`minimum_commitment = 1.0`、冷启动/稳态/因果认证/四初态全部通过。

### 1.2 三级：`threebit51_selected_v1`

`status = "FROZEN"`、`name = "threebit51_selected_v1"`、`supersedes = null`。

| 项 | 值 |
|---|---|
| `n_A1_gate` | 6.0 |
| `carry1_mrna_half_life_min` | 2.0 |
| `a1_f1_maturation_half_life_min` | 32.5 |
| `f1_production_exponent` | 4.0 |

`model.zeng_table = "unmodified"`；`f1_production_hill` 保持冻结的 `ZENG["n_A"][1]`，**只有门臂携带 `n_A1_gate`**（解耦由构造保证）。

**读出**（`readout`）：

| 指标 | 值 |
|---|---|
| `clock_period_h` | 10.60 |
| `read_window_h_median` | 2.1333 |
| `reads_total` | 56 |
| `reads_steady_drop8` | 48 |
| `read_windows_clipped` | 0 |
| `unlabelled_reads` | 0 |
| `commitment_min_all_bits` | 1.0 |
| `bit2_margin_to_band_low` | 0.282823 |
| `bit2_margin_to_band_high` | 0.278372 |
| `bit2_low_window_max_max` | 0.0171771 |
| `bit2_high_window_min_min` | 0.978372 |
| `setup_min_h_global` | 4.03039 |
| `hold_min_h_global` | 2.91372 |
| R/F pattern | `FRFRFRFRFRFRFR` |

**结构护栏**：`prefix_guard_passed = true`（`rhs_gap_n=4..8` 全为 `0.0`）、`eight_state_prefix_guard_passed = true`（`eight_state_no_feedback_gap = 0.0`）、`decoupling_verification_passed = true`。

**八初态**：`8/8 strict mod-8`，`expected sequence == observed`。配对泄漏 n=8：min `0.0142342614514237`、max `0.0142621389888085`、median `0.01424080981304355`、spread `2.7877537384798923e-05`、relative_spread `0.0019575809066184648`（0.196%）、`all_certified = true`。
初态构造：晚期轨道导出的完整 51 状态，六个上游状态与共享 C31 mRNA（index 6, `b0_M_I`）重置到同一时钟相位，**S0/S1/S2 从未被编辑**。

**严格容差**：`strict_tolerance_tol_verdict.json` `rows = 15, ok = 15`（5 点 × baseline/tight/ultra），`all_points_stable = true`；输出网格轴 `strict_tolerance_grid_verdict.json` `rows = 2, ok = 2`。

**池扰动**：`pool_perturbations_all_verdict.json` `rows = 136, ok = 136`，`any_failure = false`，`all_eight_intact_under_medium = true`。
`overall = {lost_lock: 0, recovered_original: 128, legal_mod8_shift: 8, readout_ok_no_causal: 0, alternation_broken: 0, errors: 0}`。
S2 全翻转：八初态一致地把计数平移 **+4 mod 8** 并继续计数，按合法相移吸收，不计失败。

**投影等价**：`projected_state_max_abs_error = 1.0345502232667059e-08`、`projected_rhs_max_abs_error = 0.0`，51 状态投影码串与冻结两比特均为 `1230123012301230123012301230`，两者 `certified` 均 `true`。

---

## 2. 版本链与证据哈希

```
provisional_n6 (09-20 23:16)  →  leak_correction (09-21 00:19)  →  selected_v1 (09-21 16:02)  →  postfreeze 附件 (09-21 17:29)
   7CE9D5EB…D2973                 484745AF…C77B                     6BFD7FF6…3A36                  （不改冻结件）
```

| 对象 | SHA256 |
|---|---|
| `threebit51_selected_v1.json` | `6BFD7FF6F0A4D897753D02789E5B22879D1A9DE94D53C4C8253F63D52E953A36` |
| `threebit51_provisional_n6.json` | `7CE9D5EB256045D837569AB106DEA6743B821C9BEDCD4FCCAD0CE7459E9D2973` |
| `threebit51_leak_correction.json` | `484745AFF81E34374F14B383EB4D1640A407A6C652D30426BC569ACC6479C77B` |
| `model_threebit51.py` | `D244D3CD1914F50F994E6856084D4A012B820C13249DCBA0E963E444CE3C34CC` |
| `model.py` | `256D2D104EECB4CE17F27CF4DBC08257BEEEA79E5805FD7D5614E47BBC7E543A` |
| `model_twobit34.py` | `4385E3B8C9B8AE39D09061A9DF7D5B1E7F645AC5CF56B22E4C7183B421501062` |
| `selected_profile.json` | `0571047D63AD0FABCD67A0B0E833C716A01C8E923AAC0F1860DEB2A302B2A768` |

`selected_v1` 通过 `cites` 用哈希钉住 provisional 与 correction；`evidence_sha256` 含 24 项证据；另有 `frozen_source_sha256` 与 `hash_manifest`（指向 `plausibility/SHA256SUMS.json`）。

---

## 3. 统计口径（易错，必须逐字照抄）

### 3.1 泄漏指标：三种口径不可混用

冻结中心点为 n=6、5% 尾巴切分、**2 min 输出网格**。同一名义条件下的三个数：

| 口径 | 值 |
|---|---:|
| 「n=6」四个参数组合的**组内中位数** | **0.0168830** |
| **正式中心单点** | **0.0142177** |
| **八初态**结果的中位数 | **0.0142408** |

按 n 的组内中位数（均 4 点/n，5% 切分）：

| n | `L_symmetric_5pct_median` | certified |
|---:|---:|---:|
| 4 | 0.4488681943757129 | 0/4 |
| 5 | 0.0520612629019367 | 4/4 |
| 6 | 0.01688300301854215 | 4/4 |
| 7 | 0.01385250868912815 | 4/4 |
| 8 | 0.0130503595937485 | 4/4 |

排序 `8 < 7 < 6 < 5 << 4` 在 1% / 5% / 10% 三档切分下完全一致（`ordering_stable_across_cuts = true`）；`n6_beats_n5 = true`；n=5 是离散通过边界，n=6..8 构成平台。

**两个正交的轴，不得混为一谈**：

- **尾巴切分轴**（1% / 5% / 10%）：n=6 组内中位数 = `0.01301566 / 0.01688300 / 0.02383767`
- **中心点在同一三档下的值** = `0.01241157 / 0.01421772 / 0.01764609`

把「n=6 组内中位数的三档值」写成「中心点的三档值」是错误口径。

指标定义（`carry_pairing_eight_verdict.metric_definition`，逐字）：

```
classification : R if S2 at carry start >= 0.5 else F
pairing        : one R + one F per 8-read super-period; end singles dropped
L_rev          : far_off_rev(F) / gate_on_rev(R)
L_fwd          : far_off_fwd(R) / gate_on_fwd(F)
L_symmetric    : [far_off_rev(F) + far_off_fwd(R)] / [gate_on_rev(R) + gate_on_fwd(F)]
note           : all sums inside a pair first, statistics over pairs second
```

**报告规则**：一个 R 加一个 F 窗口组成一个超周期；先在对内求和，再跨对统计；尾巴单独保留；三档切分都要报。退休的旧单窗口 `leak_ratio` **不得与配对指标并列**。

### 3.2 输出网格敏感性

中心点收敛检查：2 min → `0.014217717356355263`；1 min → `0.014162346932809777`；0.5 min → `0.01415636280346843`。
相对变化：2→1 min = `0.0038945`；1→0.5 min = `0.00042243`；2→0.5 min = `0.0043154`。每次减半变化约降 9 倍。码串在三个网格上完全相同。
积分器收紧 5000 倍带来的变化为 `3e-10`。

**正确表述**：

> 跨网格比较必须注明网格和离散误差，不能把不同网格的结果直接混用。

**禁止表述**：「跨网格不可比」——本次 2→1→0.5 min 本身就是有效的网格收敛检查。
**同时注意**：不能由该中心点的 0.43% 偏差推断**所有参数点**的排序都对网格不敏感。

### 3.3 `anomaly` / `dropped` 语义（`scan_carry_pairing.py:137 pair_rows`）

`pair_rows` 只在两种情况下记 `anomalies`：

| reason | 含义 |
|---|---|
| `same type adjacent` | 相邻 carry 同型（`a['type'] == b['type']`） |
| `not evaluable` | 配对中存在 `far_off_evaluable = False` 的窗口（有限仿真终点缺少完整 far-off 区间） |

`dropped` **只**记录循环末尾剩下的单个窗口（`if i == n - 1`，`reason = 'incomplete single window at the end - dropped, not padded'`）。

因此：

- `n_dropped = 1` **不等于**「总共只排除了一个窗口」；因不可评估而跳过的窗口记在 `anomalies` 里。
- 八初态：6 组各有 1 个末端 `not evaluable`，其余 2 组为 0；**没有同型相邻异常**。`any_anomalies = true` 来自有限仿真终点，**不是隐藏失锁**。
- 网格 16 个成功点：各有 1 个末端不可评估配对。
- n=4 的四个失败点：7 / 11 / 9 / 11 个异常，主要为同型相邻。

**禁止表述**：把「0–11 个异常」整体归给冻结的 n=6 模型。

---

## 4. 扰动证据的时相分层（不得笼统称为「补齐」）

`postfreeze_verdict.json`：`total_runs = 36`，其中 `baselines = 6`、`perturbed = 30`、`degenerate = 0`、`effective = 30`、`recovered_original = 30`、`legal_mod8_shift = 0`、`lost_lock = 0`、`readout_ok_but_causal_fail = 0`，`min_steady_reads_observed = 33`，时限 450 h（= 40 窗口 × 10.60 h，为满足 drop-8 后 ≥32 个稳态读窗）。

**正确表述**：36 条轨迹全部通过，其中 **30 条扰动、6 条对照**。
**禁止表述**：「36 次扰动全部通过」。

扰动时相（`pool_phase_structure.lead_time_before_gate_h`）：

| 对象 \ 方向 | F 型 carry | R 型 carry |
|---|---:|---:|
| `Int2_full` | −1.067 h（门内附近） | −1.067 h（门内附近） |
| `C2` | +40.70 h（提前约 41 h） | −1.600 h（门内附近） |
| `RDF2_full` | +40.77 h（提前约 41 h） | +29.0 h（提前 29 h） |

`degenerate_at_spec_named_instant = true` 出现在 `RDF2_full|F`、`C2|F`。

因此证据强度必须分层：

- **真正的门内近时相检验**：`Int2_full` 双向、`C2` 反向；
- **经过驻留期后恢复计数的检验**：`RDF2_full` 正向与反向、`C2` 正向。

**禁止表述**：「峰相位把所有弱检验都补齐了」，或把 `RDF2_full|R`（提前 29 h）称为脉冲内检验。

---

## 5. 生长条件：物理解释阻塞，不是确定性计算阻塞

`gamma` 已由曾墨涵确认 = **本征降解 + 生长稀释**的总清除率，故正式模型 `add_growth = false`；`add_growth = true` 会重复计入稀释，仅作历史诊断。

若要求同一生长条件且本征降解非负：

```
gamma = delta + mu,  delta >= 0   =>   mu <= 0.6 /h
mu = 60*ln2/Td                     =>   Td >= 69.3 min
```

- 允许 `delta = 0`（边界）：`Td >= 69.3 min`
- 要求严格正降解：`Td > 69.3 min`
- 上游振荡器：`Td = 50 min` → `mu = 0.8318 /h`，**不满足**

**正确表述**：该不一致阻止我们声称「当前上下游参数已对应同一细胞条件」，但**不阻止**冻结或报告现有 ODE 的数学结果。即使问到曾墨涵的 μ，也只是**明确差异**；要实现共同生长条件仍需重新定标并重跑受影响模块。

---

## 6. 绝对量级：计算对象与遗漏范围

`absolute_scale.json`：`hours = 300`、`uM_per_au = 5.75`、`max_free_pool_uM = 40.70689509998698`、`total_synthesis_copies_per_cell_h = 213169.58161990505`、`burden_fraction = 0.05329239540497626`、`burden_flagged = false`。
`flagged_pools = [bit0.T, bit0.R, bit1.T, bit1.R, bit2.R, F0]`。
`reference = {free_tf_uM: 30, burden_fraction: 0.1, protein_synthesis_budget_per_h: 4e6}`。

**引用前必须交代的四点**（`check_absolute_scale.py`）：

1. 脚本调用默认 `build_threebit()`（`:53`），**未设置独立 `n_A1_gate = 6`**，对应默认 n=4 配置 —— 版本与冻结中心不一致；
2. 5.33% 由各物种**各自峰值之和**得到（`:120` `copies_per_cell_h_max.sum()`），**未必发生在同一时刻**；
3. 该求和**不包含六个上游振荡器状态**对应的蛋白合成负荷；
4. 30 µM、10%、400 万 copies/cell/h 都是**检查脚本的参考假设**，不是实验测得的安全边界。

**禁止表述**：「冻结 n=6 完整电路占细胞翻译预算 5.33%，已证明表达可行」。
**同时**：六个浓度 FLAG 不是「已证明不可能」，**同样不能省略不报**。正确用法是早期量级筛查，并连同 1 fL 等换算假设、统计方式与遗漏范围一并说明。

---

## 7. 待审计结果（不得背书）

### 7.1 局部可辨识性

`plausibility/local_identifiability/local_identifiability_verdict.json` 已落盘：`hours = 350`、`analysis_start_h = 200`、`sample_min = 2`（实验观测 30）、`steps = [0.01, 0.02, 0.05]`、`rank_tol = 0.001`、9 个参数、17 个理论观测量 / 5 个实验观测量；step 0.01 下 `effective_rank = 8`、`condition_number = 1236.665`。

**状态**：

> 局部可辨识性计算已完成初版，**差分收敛与解释待审计**。

差分步长 1% / 2% / 5% 改变时灵敏度矩阵变化约 46%–68%，有效秩在 8 与 9 之间变动，**尚未展示稳定的局部差分估计**。

另：文件中「`kon`/`koff`/`decay` 只通过 `q` 进入模型，因此只能识别 `K_complex`」的解释**过强**——在显式 ODE 中它们分别进入结合、解离与复合物损失，不能把准稳态 Hill 映射的性质直接当作完整动态模型的结构不可辨识证明。

### 7.2 该包的自我限定（`conditional_on`，逐字要点）

理想无噪声观测、无测量误差模型、状态按模型 a.u. 绝对定量、初值为冻结默认初态（非拟合）。这些是**实验的**假设，不是模型的假设；可辨识性只在其内成立。

---

## 8. 历史结果 vs 过时表述（区分保留与更新）

### 8.1 必须保留为设计迭代证据（不得删除、不得改写）

**第 3 章的原 800 点统计**是有效的历史实验结果，**没有被 `selected_v1` 取代**：

```
原粗网格 drop0 : 0/800
原粗网格 drop4 : 22/800
drop6          : 46 点
最长烧入        : 48 点
```

并且当前第 3 章已写明**不能外推到全部参数或正式中心**。新旧判据同轨迹比较：两者都通过 2 点、仅新判据通过 20 点、仅旧判据通过 0 点。`2.8%` 只是**旧网格覆盖率**，不可解释为真实可行区域宽度。

### 8.2 必须更新（不是删除）

`dshwork/03_第4章_多比特级联.md` 仍写「三级模 8 未完成」「独立门指数尚待实施」——这部分应更新为 n=6 三级冻结结果，而**不是**删除全部旧扫描数字。

### 8.3 两类失败是两件不同的事，不得全部归因于后处理

1. **n=4 的三级模型真实计数失败**：低态写入不足、错误方向作用、时序错配；
2. **泄漏统计的度量错误**：尾巴混入、FR 奇偶中位数混叠。

n=6 成功来自**独立门指数改变了动力学行为**；修正统计量只是让泄漏比较变得可靠。**两者都要保留**，不得把早期失败改写成「只是泄漏指标错了」。

### 8.4 历史 PPT 的处理

`组会9.15.pptx`、`ZUHUI.pptx` 应标为**当时的设计状态**，或由新结果**补充**，不宜把历史问题改写成从未存在。

---

## 9. 开放问题清单

| # | 问题 | 性质 |
|---|---|---|
| 1 | 生长条件未统一（`mu <= 0.6/h` vs 上游 `0.8318/h`） | 物理解释阻塞；需曾墨涵 μ 或倍增时间；共同条件需重新定标重跑 |
| 2 | `uM_per_au = 5.75` | 未标定接口假设，非实测换算 |
| 3 | `alpha_rep / gamma_rep` = 3.5 / 0.6（PDF）vs 3.2 / 0.7（方波代码） | 未决，待曾墨涵确认；现按 PDF 保留 |
| 4 | `absolute_scale` 六个浓度 FLAG + 版本/口径/遗漏 | 见 §6 |
| 5 | 显式复合物 | 相对曾墨涵约化模型的结构性添加，螯合部分游离 Int 池 |
| 6 | 谷相位 basin 证据对 Int2 / C2 偏弱 | 初态采在通量谷，池处于极小值；峰相位附件只部分覆盖（§4） |
| 7 | 随机性 | 纯确定性；无 Gillespie / 质粒分离 / 单细胞成功率 |
| 8 | 局部可辨识性 | 初版已算，差分收敛与解释待审计（§7.1） |
| 9 | 独立门指数的生物实现 | `n_A1_gate` 是模型接口假设；实际构建需选择传感 TF / 杂合启动子 / split 方案 |

---

## 10. 引用规则（写 Wiki 时逐条遵守）

1. **数字必须带口径**：泄漏值必须同时写 n、尾巴切分、输出网格、以及是「组内中位数 / 中心单点 / 八初态中位数」中的哪一个。
2. **不得把诊断配置写成冻结配置**：`add_growth = true`、`uM_per_au ≈ 1`、`k_int = 6`、默认 n=4 的 absolute-scale 结果均为诊断或历史，不是候选工作点。
3. **不得把模型结果写成实验结果**：全部为确定性模型预测；`uM_per_au` 未标定；无湿实验闭环。
4. **失败与限制随文保留**：n=4 真实失败、六个浓度 FLAG、生长条件不一致、可辨识性未审计，都是设计证据，不删不藏。
5. **证据强度分层**：`[certified]` 产物内认证 / `[mechanism-inferred]` 机制与扫描支持 / `[pending audit]` 待审计 / `[pending calibration]` 待实验标定。
6. **跨章引用不改写历史**：第 3 章 800 点统计与第 4 章 n=6 冻结结果并列呈现，标明各自的条件与适用范围。
