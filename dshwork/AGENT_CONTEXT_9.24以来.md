# 交接上下文：2026-09-24 之后做的全部工作

**读者**：接手本项目的另一个 agent。
**写作口径**：只写**已被本次会话实测/核对过**的事实；每个数字都能在文末"产物清单"里找到出处。凡是我自己犯过的错，单列 §7，不要重复踩。

---

## 0. 三十秒速览

1. 9/24 当天发现一个**静默工作点事故**（`n_A1_gate` 回落），据此作废 18 个产物、立下 H1–H8 硬约束、并把 51 状态模型冻结在 `n_A1_gate = 6`（附正式论证：泄漏膝点）。
2. 随后完成：翻译不可观测性证明、可辨识性收尾、**随机性预审**（这是最可能改变结论的一项）、以及 10 张图。
3. **最重要的未决问题**：随机性预审的预登记规则字面输出是"必须做随机性研究"，但我逐臂归因后发现，7 个未认证点里 **4 个只是 contrast 臂单独失败**（那条指标已被项目自己降级）。**在这个归属决定之前，不要引用那句结论。**
4. **仍未闭合**：RDF 承重性、相位扩散/快噪声/质粒丢失、生长条件、`uM_per_au` 标定。
5. **写完本文档后复查时新查出的第二条工作点陷阱**：同一个 51 状态模型有**两条构造路径**给出不同 `clock`（认证链 `0.3/2.0` vs 模型默认 `0.4/3.0`），而 `0.4/3.0` 从未进入任何扰动/认证网格。**动手构造模型前先读 §1 的警告块。**

---

## 1. 冻结件与工作点（先看这一节，动任何东西前对照）

**冻结工作点**（＝认证链 `plausibility_common.frozen_extension()`＋`build_threebit()` 的实测口径）：

| 量 | 冻结值 | 出处 |
|---|---|---|
| `n_A1_gate` | `6.0` | constructor 参数（`SELECTED_N_A1_GATE`） |
| `n_A[1]` | `4.0` | ZENG 表值，**仍用于 F1 产生**（与门指数解耦，H1） |
| `uM_per_au` | `5.75` | `twobit34_results/selected_profile.json` |
| `carry0/1_mrna_half_life_min` | `2.0` | 同上 |
| `a1_f1_maturation_half_life_min` | `32.5` | 同上（carry0 = carry1） |
| `clock_K_au` | **`0.3`** | 同上（`extension.clock_K_au`） |
| `clock_n` | **`2.0`** | 同上（`extension.clock_n`） |
| `add_growth` | `False`（**只是不额外加生长稀释项**，见 §4 第 4 条；**不要写成 μ=0**） | 同上 |

> ⚠️ **`clock` 是第二个"静默工作点"陷阱，与当初的 `n_A1_gate` 同类，务必先读完这一块再构造模型。**
> - **认证链**走 `frozen_extension()`（读 `selected_profile.json` 的 `extension`）⇒ **`clock_K_au = 0.3`、`clock_n = 2.0`**。600 h 认证网格、随机性预审、冻结点泄漏 `0.014218` 全部出自这一条路径（`scan_carry_pairing.py:99`、`preaudit_stochasticity.py:223` 都显式 `extension=frozen_extension()`）。
> - **`ThreeBit51Model()` 不传 `extension=`** ⇒ 落到 `selected_threebit_extension()`（`model_threebit51.py:34`）⇒ **`0.4 / 3.0`**；`test_threebit51.py:46-47` 断言的也是 0.4/3.0。`scan_threebit51_carry1_dsh.py:123` 正是 `ThreeBit51Model(carry=carry)` 这种写法，所以 `threebit51_results/carry1_scan_dsh/meta_shard*.json` 与 `threebit51_results/diagnostic/parameters.json` 记录的是 **0.4**。
> - **后果（本次实测，`n_A1_gate = 6`）**：`b2_I`（state 36）峰值 **0.3/2.0 → `1.7784`**，**0.4/3.0 → `1.1711`（−34.2 %）**。复现配方：`solve_ivp(model.rhs, (0, 200), model.initial_state(cold=False), method='LSODA', rtol=1e-8, atol=1e-10, max_step=0.05)`，**自变量单位是小时**（`model.py:126` 在 RHS 内部才 `upstream.rhs(t*60, ...)` 换算成分钟；`max_step=0.05` 即 3 min），取 `sol.y[36].max()`。也就是说，照本节重建工作点却漏传 `extension=`，轨迹会被静默换掉——这正是 §1 必须钉住 clock 的原因。
> - **0.4/3.0 不在预审扰动族里**：`preaudit_stochasticity.py:230` 只缩放 K（`replace(frozen, clock_K_au=frozen.clock_K_au * f)`），**n 固定为 2**；认证因子区间 `[0.8, 1.2]` ⇒ K ∈ `[0.24, 0.36]`，且 `touches_grid_edge = True`。所以 K = 0.4（1.333×）与 **n: 2→3 从未被扰动**，两者都**未经认证**。
> - **没有任何断言检查 clock**：H1 的断言只覆盖门指数（`plausibility_common.py:147` 的 `n_A1_gate_effective`）。下一个 agent 若想让代码自己钉住 clock，得自己加断言——现在没有。
> - **想改代码消除这条分歧要先读这句**：`model_threebit51.py` 与 `test_threebit51.py` **都在 `FROZEN_FILES` 里**（`plausibility_common.py:63-65`），其哈希进入 `source_hashes()`，而 `load_valid_decoupling()` 会比对 `art['source_sha256'] != source_hashes()`（`plausibility_common.py:668`）⇒ 改完必须重跑 `verify-decoupling`，否则门会拒绝。

**只读冻结件（哈希已核对，勿改）**

| 文件 | SHA256 |
|---|---|
| `final_reconstruction/model.py` | `256D2D104EECB4CE17F27CF4DBC08257BEEEA79E5805FD7D5614E47BBC7E543A` |
| `final_reconstruction/model_twobit34.py` | `4385E3B8C9B8AE39D09061A9DF7D5B1E7F645AC5CF56B22E4C7183B421501062` |
| `final_reconstruction/model_threebit51.py` | `D244D3CD1914F50F994E6856084D4A012B820C13249DCBA0E963E444CE3C34CC` |
| `final_reconstruction/verify_threebit51.py` | `9AC2B6A9AD83BBEC41A87891A620EFE4EA050ACDAD8BF39757FD5F14281B8827` |
| `final_reconstruction/verify_twobit_causal.py` | `C8678A03EC6201E65EBBE2BD182AB89BC2AB8DAAD502A4E302CAFFBC361BCFF0` |
| `plausibility/threebit51_selected_v1.json` | `6BFD7FF6F0A4D897753D02789E5B22879D1A9DE94D53C4C8253F63D52E953A36` |
| `plausibility/threebit51_provisional_n6.json` | `7CE9D5EB256045D837569AB106DEA6743B821C9BEDCD4FCCAD0CE7459E9D2973` |
| `plausibility/threebit51_leak_correction.json` | `484745AFF81E34374F14B383EB4D1640A407A6C652D30426BC569ACC6479C77B` |

**清单**（2026-09-26 收尾轮之后的重生成值）：
- 内层 `plausibility/SHA256SUMS.json` = `9143A7A7303819C4858B39DFDD868BCA82E78A993BFF23C941B8BFE6E84855D4`（**943 条**；收尾轮前为 932 条 / `016CDB86…`，旧清单已归档在 `plausibility/preaudit_report_archive/ARCHIVED_plausibility_SHA256SUMS.v1.json`）。
- 根级 `SHA256SUMS_authority.json` = `1E12E1361CC2D7F2F044E3102756D313F8D1350893D46CD2967293F2D89CA039`（**30 条**；收尾轮前为 27 条 / `800B8DB9…`；新增的是预审报告、分项验收表、预审归档哈希三件）。
- **更新顺序不能颠倒**：`make_authority_manifest.py` 本身就在 `plausibility/` 内，所以动它之后必须**先**重生成内层清单（`plausibility_common.write_manifest()`）、**再** `make_authority_manifest.py write`、最后 `check`。搞反会让内层清单立刻过期。
自检命令：`D:\aconade\python.exe plausibility/make_authority_manifest.py check`（会**重算** `plausibility/SHA256SUMS.json` 每一条；本目录任何改动都会让它响亮失败——这是有意的）。**别用裸 `python`**，见 §7 第 6 条。
解耦守卫：`decoupling_verification.json` 的 `implementation_sha256 = 7F4658472384B730892591133B2FC6AC37954914D557C5A867FE504057C271DA`，改 `model_threebit51.py` 或 `plausibility_common.py` 后必须重跑 `verify-decoupling`。

---

## 2. 按主题的全部工作

### A. 工作点事故与整改（9/24，一切的起点）

- **事故**：多个检查脚本调用 `build_threebit()` 时未显式给 `n_A1_gate`，于是模型静默回落到 ZENG 表值 **4.0**；49/55 次局部可辨识性运行、`check_threshold_margin`/`check_absolute_scale`/`check_mass_balance`/`check_rdf_collapse`/旧 `check_identifiability`/`scan_sharpen_vs_lengthen` arm A 全部工作点错误。
- **处置**：18 个路径/5 个脚本组的产物标记 `INVALID_AT_FROZEN_POINT`（`plausibility/INVALIDATED_ARTIFACTS.json`），n=4 产物归档到 `plausibility/invalidated_n4/`（含 `ARCHIVED_HASHES.json`），并逐项重跑。**已核实等价**：constructor `n_A1_gate=6` ≡ 全局补丁 6 + `n_f1_drive=4`（前缀 RHS gap 恰 0.0，尾部 4.44e-16／相对 8.84e-19）。
- **制度化**：`plausibility/HARD_CONSTRAINTS.md` H1–H8；`CURRENT_BASELINE.md` 新增 **§11「判据与验证纪律」**（H1–H8 的交付版）；`plausibility_common` 的 `build_threebit` 现在**断言**实测指数，`working_point_block(model)` 读回工作点而不复述意图。

### B. 模型性质的三项收尾

1. **`translation_h` 结构性不可观测**（`check_translation_unobservable.py`）：13 个 mRNA 状态全部列出后，**12 个移动**（`M(T)/M(30) = 30/T` 精确：T=45 → 0.6667，T=15 → 2.0000），**`b0_M_I` 不动**（gap 9.09e-16/4.44e-16，比值 1.0000，因为它是共享的上游 C31 mRNA）；可观测状态差 1.08e-13、通量差 1.84e-13、读出 `1234567012` 三档一致。⇒ `translation_h` 只能按假设固定，**不能称为"标定"**。
2. **等剂量成熟时间三态**：只有参考点 n=6 落在网格内；n=5/5.5 为 `below_grid`（13.98／23.88 min，下界 25），n=6.5/7 为 `above_grid`（48.39／62.71 min，上界 45）。**补偿量有两个不能混用的口径**：切线型 23.638–23.812 min/单位 n（仅无穷小步长成立） vs 割线型 17.237–31.789；把 n 从 6 挪到 6.5 买回冻结剂量，切线给 44.347 min、网格实测 48.395 min，**差 4.047 min**。
3. **剂量交易不传导到功能**：固定 n 扫遍成熟时间，泄漏最多动 **2.65 %**（绝对 3.3e-4），而半步 n 动 0.00059–0.0282（1.8×–85×）。网格下界保留意见：5 行里 4 行最小值落在下界 25 min，线性外推到 mat=0 的预算 4.3e-4 仍小于最小半步 5.89e-4（比值 **0.73**）。

### C. n=6 冻结的正式论证（`plausibility/N6_FREEZE_JUSTIFICATION.md`）

600 h 认证网格（n=4..8 × mRNA × mat，20 点，`carry_pairing_grid_verdict.json`）逐 n 中位数（5 % 切尾）：n=4 不认证（0/4，crossings 1–6/14）／n=5 `0.052061`／n=6 `0.016883`／n=7 `0.013853`／n=8 `0.013050`。
- 排序 `8<7<6<5<4` 在 **1 %/5 %/10 % 三档完全一致**；
- 逐步比值有**两个不能混用的口径**（源文件 §3 把两者放在两张分列表格里，§5 明文写"两者都报，不混用"）：
  - **整数步 / 冻结单点**（600 h，冻结点 5 % 切尾）：5→6 = **3.496×**、6→7 = 1.232×、7→8 = 1.044×；对应的 `dlnL/dn` = **1.252 → 0.209 → 0.043**（跨 n=6 后边际收益先降 6.0×、再降 4.9×）。**这三个 `dlnL/dn` 只配这三个比值**（`ln 3.496 = 1.2517`、`ln 1.232 = 0.2086`、`ln 1.044 = 0.0431`）。
  - **4 点中位数**（600 h，5 % 切尾）：5→6 = **3.084×**（= `0.052061/0.016883`）、6→7 = 1.219×（= `0.016883/0.013853`）、7→8 = 1.062×（= `0.013853/0.013050`）。
  - **3.084× 与 3.496× 不是同一个量**（前者 4 点中位数、后者冻结单点）。引用时必须写明口径，**尤其不要拿中位数比值去配整数步的 `dlnL/dn`**（`ln 3.084 = 1.1262 ≠ 1.252`）；
- 半步口径膝点落在 n ≈ 5.5–6，n=6 是**位于膝点上且已认证的最小整数**；
- 文档同时写明**不声称**的事：n=7/8 泄漏更低（这是简约性选择不是极值选择）；膝点读成 5.5/6.5 不改变决策，读成 8 才会。
- **必须保留的边界（复核意见）**：`n_A1_gate = 6` 是模型里的**有效响应陡峭度**，**不自动对应"六个结合位点"或某个可直接搭建的启动子**。实际实现仍需实验的剂量—响应曲线支撑。因此 n=6 是**工程折中**（已认证 + 收益开始变缓 + 采用较小指数），不是已验证的生物实现。

### D. 随机性预审（`preaudit_stochasticity.py`，最可能改变结论的一项）

- **`h2check` 通过**：三种冻结点表述（冻结 uM 的 extension／冻结 clock_K 的 extension／把 K_A1 patch 回自身冻结值）100 h 轨迹 **`max|Δstate| = 0.0`**（逐位相同）⇒ 无 ZENG 跨模型渗漏。
- **静态容差扫描**（21 点 × 600 h，**跑两遍**：run1 `D724B361…`、run2 `4B393567…` 十个可比字段**逐位相同**）：认证区间宽度 `conc_scale **0.00**`（只有 f=0.95/1.00 认证）、`K_A1 0.05`、`clock_K 0.20`；预登记规则字面输出 **SMALL →「随机性研究是必需的」**。
- **但归因改变了这句话**：`certified = steady ∧ one ∧ order ∧ alt ∧ contrast`，7 个未认证点里 **4 个是 `contrast` 臂单独失败**（off/on 门峰值比 > 0.10），它们的**数字读出完全正确**（14/14、序列精确、0 未标定），其中 `conc_scale|1.05` 的配对泄漏 **0.01069 还优于冻结点 0.01422**。冻结点自身在该臂只有 **0.0545/0.10 = 54.5 %，余量 1.83×**。而这条 off/on 峰值比正是 **H4 已为泄漏量化降级**的指标。
- **小数字（独立成立）**：`copies_per_au = 3462.7309`（其中 `cell_volume_fl = 1.0` 是**假设**）。`b2_I` 的门相分布**极不均匀**：最小 **20.0**（CV ≈ 22 %）、p05 **38.6**（CV ≈ 16 %）、**中位数 2159.8**（CV ≈ 2.2 %）；全周期中位数只有 **1.97**（CV ≈ 71 %）。`F1` 全周期谷值 **40.2**（CV ≈ 16 %）；显式复合物 `b2_C` 门相中位数 **1.48**（CV ≈ 82 %）；`A1` 门相最小 2368。
  - ⚠️ **不要把"门相最小 20 个分子"说成"主要写入发生在 20 个分子时"**：门相**中位数是 2160**，20 只是极低尾。要下这个判断，必须把**低拷贝时段**与**实际正/反向重组通量**对应起来 —— 这一步**还没做**。
  - ⚠️ `CV ≈ 1/√N` 是**泊松计数的参考量级**，不是这个反馈网络已经测得的噪声；`cell_volume_fl = 1.0` 与 `uM_per_au = 5.75` 都是未标定假设，拷贝数线性依赖后者。
  - ⚠️ "把 `uM_per_au` 提高五倍就把最小结合点从 20 抬到 100" **只对"保持 a.u. 轨迹不变"的换算估计成立**（`uM_per_au_for_100_copies / 5.75 = 4.99`）。代码里这个参数**同时改变动力学输入**，真去改参数必须重算。
  - ⚠️ **`conc_scale` 那一维实际改的就是 `uM_per_au`**（`replace(frozen, uM_per_au=frozen.uM_per_au*f)`，`preaudit_stochasticity.py:226`），把它统称为"全局浓度增益噪声"会误导。
- **判据三次修订**（全部留在产物里）：Rule A 命中 `b2_C = 5.1e-9`（那是"该时刻不存在"）、Rule B 又把全周期最小值折回敏感性得出 `uM_per_au` 需 9.4e4×（荒谬）、Rule C **拒绝给单池判定**，改报 per-pool/per-phase 全表。**一次会话两次修订说明"最小池"这个统计量本身是错的工具。**
- 报告：`plausibility/PREAUDIT_STOCHASTICITY_REPORT.md` **已修订为 v2**（v1 逐字归档在 `plausibility/preaudit_report_archive/`）。v2 相对 v1 改了 6 处：小数字改为「门相 最小 20.03 / p05 38.56 / **中位 2159.8**」并声明**尚未确定主要写入是否发生在低拷贝时段**；`1/√N` 降级为**泊松参考量级**而非已算出的 CV；补「改 `uM_per_au` 会改变动力学，不能只按比例换算后沿用原轨迹」；结论由「随机性研究是必需的」改为「**值得进一步评估随机效应，目前尚未确定其对计数可靠性的影响**」；旋钮名改为「**`uM_per_au` 接口参数扰动**」；指向新的分项验收表。
- **分项验收表（必须用它，别只引单一布尔）**：`plausibility/preaudit_acceptance_breakdown.{md,json,csv}`，由 `make_acceptance_breakdown.py` **纯后处理**生成（不建模型、不积分），与 verdict 的 `certification_table`/`leak_table`/`off_on_gate_peak_ratio_table`/`failures_by_arm` 逐点交叉核对 0 问题。
  - 新增标签 `counting_and_causality_passed := steady ∧ one ∧ order ∧ alt`。**不是 `certified` 的改名，也不附任何泄漏阈值**。实测历史 `certified` **14/21**、新标签 **18/21**；4 个点是 **contrast 臂单独失败**（计数与因果全过）。
  - **离散通过点**：`uM_per_au` 维在历史谓词下通过 {0.95, 1.00}、在计数+因果下通过 {0.95, 1.00, 1.05, 1.10}；`clock_K` 全 7 点；`K_A1` 为前 5 点 / 全 7 点。**都只是已测网格点，未采样因子从未运行，不得写成连续可靠范围。** 该维**失败原因随因子变化**（<0.95 败在 steady/causality，>1.00 只败在 contrast）。
  - **一处如实报告的不一致**：v1 的「3 个点的数字读出完美却未认证」**不可复现**；两个清晰定义都给 4 且成员不同（表里列了 A/B 两式与差异点）。**不要再引用那个 3。**
- **v2 里 `uM_per_au` 认证集合若按 H4 重定**会变成 {0.95, 1.00, 1.05, 1.10}，对称宽度 **0.05** —— 而原判据是 `SMALL if width < 0.05`，**0.05 正好在边界上**，所以"SMALL"这个判定本身要重审，**报告里没有代为决定**。
- **结论边界**：本预审**不能**回答快噪声（周期 ≲ 进位周期 42.4 h）、时钟相位扩散、质粒分离/低拷贝灭绝、单细胞成功率；也**没有**证明冻结电路在真实随机动力学下仍可靠。

### E. RDF 承重性（**OPEN**）

`plausibility/RDF_PROBE_DESIGN.md` 是**预登记设计**（尚未实现）：四个旋钮（bit2-only 的 `K_rep/n_rep`、`alpha_rdf`、`gamma_rdf`、`alpha_rep`）、非线性 vs 纯剂量对照、三种结局 A/B/C、以及反证条件。关键教训来自两次失败：v1 全局 patch `ZENG['n_rep']/['K_rep']` **毁掉冻结前缀**（`reverse_events` 14→0，判决已撤回）；v2 改成 bit2-only 后**判决空洞**（n=6 本来就计数，"能恢复"恒真），现字段为 `null` 并打印 `vacuity_text()`。**正确方向是反过来：从能计数的工作点出发把池往下压，测边界与余量。**

### F. 图（全部在 `dshwork/figures/`，脚本自带自检）

| 图 | 文件 | 自检项数 |
|---|---|---|
| 扫描图 4 张（2bit/3bit 单参数 + 鲁棒性） | `final_reconstruction/scan_figures/` | 47 |
| 3-1 单比特记忆 | `fig03_1_single_bit_memory` | 40 |
| 3-2 翻转与读窗 | `fig03_2_switching_read_windows` | 43 |
| 4-1 carry 架构 | `fig04_1_carry_architecture` | 45 |
| 4-2 mod4 → mod8 | `fig04_2_mod4_mod8` | 33 |
| 4-3 更陡的门 | `fig04_3_sharper_gate` | 30 |
| 4-4 验证摘要 | `fig04_4_validation_summary` | 33 |

纪律要点：机制图与数据图都**导出 SVG 且保持真 `<text>`**（`svg.fonttype='none'`），因此标签可被逐字校验、也方便美工改字；LEGACY（`n_A[1]` 共享槽位）与 FROZEN 面板**必须分开标注**；泄漏一律写清 n／切尾档／输出网格。旧图（`fig03_single_bit_mechanism` 等 6 张）**未被改动**。
**4-3 需要一次新运行**：`dshwork/figures/data/n4_matched_600h.{csv,json}`（由 `make_n4_matched_trajectory.py` 生成，同一 extension/carry、只有门指数取 4，脚本断言实测指数）。结果：`certified=False`、**14 次进位里只有 5 次翻转 bit2**、**17 个读窗位标签从未确定**。该 JSON 自述为"**新运行、非既有产物**"。
- 🔴 **本图曾有一个真实 bug，复查处已修（2026-09-26）**：生成器原先把列写成 `S0 = sol.y[44], S1 = sol.y[45], S2 = sol.y[46]`，即把 **`b2_S` / `A1` / `F1`** 挂在了 S0/S1/S2 名下（真值 **S0 = 16、S1 = 27、S2 = 44**）。于是 Panel A 的 **n=4 面板画的是 F1、n=6 面板画的是真 S2** —— 同一个轴标签下是两个不同物理量；原图注「S2 峰值 2.25、超出物理区间被裁掉」描述的其实是 F1（以 a.u. 计的浓度，本来就不受 1 约束）。**这句话作废，不要再引用。**
- **判据数字完全不受影响**：`crossings = 5`、`gates = 14`、17 个未确定窗、`certified = False` 都来自 `verify_threebit51.analyse_threebit`，它读 `states = (y[16], y[27], y[44])`，一直是对的。
- **修复是三件事**：(a) 生成器改为**按 `model.state_names` 名字取索引并断言**、写盘后自检，列名改为显式的 `b0_S/b1_S/b2_S/A1/F1`；(b) 新增独立检查器 `dshwork/figures/check_state_column_semantics.py`（见 §7 第 8 条）；(c) 图脚本加「曲线必须落在 [0,1] 内」闸。新旧数据逐列比对**逐位相同**（新 `b2_S` ≡ 旧 `S0`、`A1` ≡ 旧 `S1`、`F1` ≡ 旧 `S2`）⇒ 修复是**纯改标签**，数值未变。
- **顺带修正的两处图设计**：共用时间窗从 `[200,320]` 改为 **`[30,150] h`**（原窗口里 n=4 的 S2 最低只到 0.309，看不到它真实的翻上去那一刻）；n=4 的图注改为数据驱动的真实读数 —— 同一窗口内 **9 个 1、0 个 0、2 个未确定**（bit2 在 t≈38 h 翻上去后再没回低态），n=6 为 **7 个 1、4 个 0、0 个未确定**（干净的 mod-8 交替）。
**4-4 的计数从产物重算**：冻结前 136 = 128 恢复 + 8 合法 mod-8 相移 + 0 失锁（`plausibility/pool_perturbations_all_all.csv` 与 `threebit51_selected_v1.json` 的 `pool_perturbations.overall` 双向核对）；冻结后 30 = 30 + 0 + 0（`postfreeze_all.csv`，另含 6 条无扰动对照）。8 次相移来自 S2 翻转；提前时间实测 1.6–40.8 h。

### G. 公式与参数体系统一（以代码为准）

把 `wiki任务/初步.pdf`（正交重组酶，每模块 11 状态：5 条辅助函数 + 12 条 ODE）与 `其他小组成员任务/前馈环脉冲进位级联模型.pdf`（I1-FFL carry：A/F 两条臂、AND 门、方案 A 负自馈、方案 B 时钟门）**统一到代码命名**。全部差异按"代码正确、PDF 侧改述"处理；下面五条是**承重的五条**（完整清单的出处见本节末）：

1. `初步.pdf` 的 TetR 由**输入**诱导；我们 `b*_T` 由 **PB 状态**驱动（`model.py:141`）；
2. `初步.pdf` 的整合酶由输入诱导；我们的 **bit0 用真实 C31 通量、bit1/2 用 carry 门输出**，且**没有** `α_{I,0}` 泄漏项；
3. RDF 泄漏 `α_{R,0}` 我们取 **0**，构象活性隐式取 `a_LR=1, a_PB=0`；
4. 正/逆向 Hill 指数我们**固定为 2**，阈值分别是 `K_D_int[i]` 与 `K_complex = q·K_D_comp`（`q = k_on/(k_off+δ_C+μ)`）；
5. 门指数拆成 **门臂 `n_A1_gate = 6` / F1 产生 `n_A[1] = 4`**，并用"只改 `d[b2_M_I]`"的精确修正实现（H1）。

**索引语义（最容易搞错）**：`alpha_Int / K_A / n_A / K_F / n_F / alpha_A / gamma_A / alpha_F / gamma_F` 按 **carry 级 j** 索引（j=0 由 bit0 触发、产出 bit1 的整合酶；j=1 产出 bit2）；`gamma_int / K_D_int` 按 **位 i** 索引。

**本节的两点边界（别越读）**：

1. 上面五条是**承重**的五条，**不是全部**：本轮比对得到的完整清单**尚未落盘**，只存在于会话记录里。盘上现有的 `方程与实现对照.md`（wiki 根，9/19）是**更早的**一份对照，只覆盖 `初步.pdf` 的 11 状态方程，且它的 §8 carry 写法（记忆蛋白 `dE_i/dt = α_E S_i − (δ_E+μ)E_i` ＋ `U_{i+1} = g_i·…`）**已被后来的 I1-FFL carry 取代**，不要拿它当本节的完整版。要用完整清单，先把统一文档写到盘上。
2. 本节谈的是 **PDF ↔ 代码** 的差异。还有一条 **代码内部**的分歧——同一条认证链有两条构造路径给出不同 `clock`（0.3/2.0 vs 0.4/3.0）——在 **§1 的警告块**里，别混进本节。

---

## 3. 关键数字速查

| 量 | 值 |
|---|---|
| 时钟周期 / 进位周期 | 10.60 h / 42.4 h（= 4 拍） |
| 600 h 读窗数 / 稳态（drop 8） | 56 / 48 |
| 读带 / 读窗几何 | `BAND_LOW=0.30`、`BAND_HIGH=0.70`；读窗 = 周期 **20 %（谷值 ±10 %）**，实测 0.2013 |
| contrast 限值 | `MAX_OFF_ON_GATE_RATIO = 0.10`；冻结点 0.0545 |
| 冻结点配对泄漏（5 % 切尾） | 0.014218（1 %: 0.012412，10 %: 0.017646） |
| 冻结点开关余量 | setup 4.03 h、hold 3.70 h |
| 拷贝数换算 | `copies_per_au = 3462.7309` |
| 写入池 | `b2_I` 门相 最小 20 / p05 38.6 / **中位数 2159.8** 个分子（分布极不均匀；绝对拷贝数线性依赖未标定的 `uM_per_au = 5.75`，见 §2D、§4 第 5 条） |

---

## 4. 未决问题与建议的下一步

1. **（零算力，最优先）`contrast` 臂的归属**：若认证谓词改按现行 H4 配对泄漏标准，则 `K_A1|1.10`、`|1.20` 与 `uM_per_au|1.05`、`|1.10` 都应通过（**计数与因果全过**，且后两者的泄漏 0.010691 / 0.009785 **优于冻结点 0.014218**）。后果要说准：`uM_per_au` 维的认证集合会从 **{0.95, 1.00}（对称宽度 0，原判据判为 SMALL）** 变成 **{0.95, 1.00, 1.05, 1.10}（对称宽度 0.05）**，而原判据是 `SMALL if width < 0.05` —— **0.05 正好落在边界上**，"SMALL" 本身要重审。若保留 off/on 峰值比，则必须把"±5 % 的 `uM_per_au` 变动"写成设计硬边界。**调整必须版本化并列重评，不能删掉失败项后沿用原认证名称。在此之前不要引用"随机性研究是必需的"。**
2. **随机性研究**：理由应来自 §2D 的小数字（**门相 p05 只有 38.6 个分子、且门相分布从 20 跨到 2160**），而不是那条规则。需实现真正的随机模型；预审已列明它不能覆盖的四类问题。**预审支持的是"值得进一步研究离散分子效应"，尚未证明冻结电路会因随机性失效，也没有估计单细胞成功率。**
3. **RDF 承重性**：按 `RDF_PROBE_DESIGN.md` 实现（设计密集型，可与随机性并行写）。
3b. **（建议的下一项计算工作）低拷贝时段究竟承担了多少写入通量**：先确认轨迹是否含 `I2`、RDF2(`b2_R`)、`C2`(`b2_C`) 与正反向通量 —— **n=6 的 `wiki_n6_20260924/trajectories.csv` 本来就有**；n=4 的 `dshwork/figures/data/n4_matched_600h.csv` 已在收尾轮补齐 `b2_I/b2_R/b2_C/J_fwd2/J_rev2`，所以**两侧都可只做后处理**。按 F/R 两类进位分别统计：`I2` 或 `C2` 低于 10 / 30 / 100 个分子时各承担了多少**预期方向**的重组通量积分；主要写入发生在哪个分子数量级；低拷贝时段位于有效写入段、脉冲尾部还是远驻留期。**这些档位只用于描述，不是失效阈值。先弄清实际写入是否依赖低拷贝时段，再设计随机模型。**
4. **生长条件（表述已改，不要再写"μ=0"）**：`model.py:88-89` 是
   `self.mu = 60*self.a.mu`、`self.growth = self.mu if self.e.add_growth else 0.0`。
   所以 `add_growth=False` 的准确含义是「**下游不额外添加生长稀释项**，使用曾同学确认的、
   已包含降解与稀释的总清除率 γ」——**不是**"细胞不生长"，也**不是**"整个模型 μ=0"
   （上游仍用自己的生长参数）。真正未解决的是**物理一致性**：上游 μ≈0.832/h **大于**下游
   某些总清除率 0.6/h，两者无法在相同生长条件下解释为非负内禀降解。这条应继续写在 Wiki 限制里。
5. **`uM_per_au = 5.75` 未标定**：所有绝对浓度线性依赖；小数字结论可被它 ×5 翻转。
6. **`check_identifiability.py`（旧）**从未执行、已被 `check_local_identifiability.py` 取代；`scan_sharpen_vs_lengthen` **不得重跑**（arm A 钉在 4.0）。
7. **drop 数不是"不一致"，是按各自计数器周期取的（此处此前误判）**：两套计数器的模数不同，规则却一致——都是"丢掉第一个完整计数周期"。
   - `verify_twobit_causal._counter_verdict` 的判据是 `(v[i+1]-v[i]) % 4 == 1`（**模 4**），故 `STEADY_DROP_READS = 4`（`verify_twobit_causal.py:32`）；
   - `verify_threebit51.counter_verdict` 的判据是 `% 8 == 1`（**模 8**），故稳态调用 `counter_verdict(reads, 8)`（`verify_threebit51.py:31,69`）；
   - `check_read_commitment.py` 已把两口径并列（`PRIMARY_DROP = 8` 是产生认证的那个、`ALT_DROP = 4` 仅供对照，见该文件 docstring 第 35-37 行）。
   **不要把 4 和 8 统一**——那会破坏其中一个计数器；引用时写明是哪套计数器即可。

---

## 5. 纪律（改任何东西之前读一遍）

- **H1** 构造 51 状态模型必须显式给 `n_A1_gate` 并断言实测值；
- **H2** 绝不同时持有两个语义依赖全局 `ZENG` 的模型实例（`_carry1_sources` 是 staticmethod，运行时读全局）；
- **H3/H7** 失效产物**先归档再重跑**，不要就地覆盖、不要删目录；
- **H4** 泄漏只能用配对 R/F 超周期 `L_symmetric`，并注明切尾档与输出网格；
- **H5** 结构性检查用 **RHS 级相对比较**，不要用轨迹逐位相等；
- **H6** 秩/条件数必须连切点和最小奇异值余量一起报；
- **H8** 不要在 workspace 里用 `tempfile.mkdtemp()`（目录删不掉；若落在 `plausibility/` 会让清单永久不一致）。

**禁止**：把 LEGACY 面板当冻结结果；引用 `INVALIDATED_ARTIFACTS.json` 里的产物；把 `certified=False` 直接写成"计数器失效"（先看是哪条臂失败）；把 `best_maturation` 当推荐值（它常落在网格下界）。

**认证规则的分层与版本（复核意见要求，尚未实施）**：后续报告应**分项并列**展示 ——
(a) 数字读出是否通过、(b) 进位因果是否通过、(c) 门峰值对比度、(d) 配对泄漏 / 时序 / 读窗裕量 ——
而不是只给一个 `certified` 布尔。**历史认证结果保留原定义**；若要调整认证规则（例如按 H4 的配对泄漏
口径重定 contrast 臂），必须**建立新版本、并列重评**，**不能删掉失败项后沿用原来的认证名称**。
H4 把峰值比降级为"不是泄漏的代理"，但这本身**还没有决定**它是否应作为一条独立设计要求保留。

---

## 6. 产物清单（按目录）

- `final_reconstruction/plausibility/`：本包主体。关键新增：`preaudit_stochasticity.py` + `_verdict.json`/`_inventory.json`/`_h2check.json`、`PREAUDIT_STOCHASTICITY_REPORT.md`、`RDF_PROBE_DESIGN.md`、`N6_FREEZE_JUSTIFICATION.md`、`HARD_CONSTRAINTS.md`（H1–H8）、`make_authority_manifest.py`、`INVALIDATED_ARTIFACTS.json`、`invalidated_n4/`、`preaudit_stochasticity_run1/`（第一批归档：8 分片 + 三次 inventory 快照）。
- `final_reconstruction/scan_figures/`：4 张扫描图 + `plot_scan_figures.py` + `selfcheck_scan_figures.py` + `figure_data_manifest.json`（输入/输出哈希）。
- `dshwork/figures/`：6 张 Wiki 图 + 6 个作图脚本 + `tempo_style.py` + `make_n4_matched_trajectory.py`（4-3 的取数脚本，已改为按 `state_names` 取名定位）+ **`check_state_column_semantics.py`（列语义检查器，新增；见 §7 第 8 条）** + `data/n4_matched_600h.{csv,json}` + `README.md`（含可直接用的中文图注与逐项验证表）。
- wiki 根：`CURRENT_BASELINE.md`（新增 §11 判据与验证纪律）、`SHA256SUMS_authority.json/.txt`。

---

## 7. 本次会话我犯过的错（请勿重复）

1. **删目录丢证据**：`Remove-Item -Recurse` 删掉 `local_identifiability`，原始撤回记录与 55 份 n=4 产物未先归档；`RETRACTION_n4_run.md` 是**重建版**，不能声称字节一致 → 立 H3/H7。
2. **归档不全**：tolerance 的 run1 归档了，但 **inventory 的第一次产物被第二次运行覆盖前没归档**；只能如实说明"共享字段逐位相同、但无字节级可追溯"。
3. **判据两次修订**（§2D）：说明聚合统计量选错了，不是阈值需要调。
4. **校验器自己会错**（至少 6 次）：`fig.findobj` 会重复返回同一 Text（行名计数虚高）；`Bbox.intersection` 在这版 matplotlib 是静态方法；条带用 `aspect='equal'` 会把刻度标签挤出画布；`per_point` 与 `by_n_A1_gate` 的键序相反（`5pct_L_symmetric_median` vs `L_symmetric_5pct_median`）；把 30 h 探针的"0 次穿越"当成 n=4 的结论（600 h 实际是 5 次）；`axhspan` 的宽度是轴分数而非数据坐标。**结论：写完检查脚本后要单独验它。**
5. **环境陷阱**：`Wait-Process` 被禁；前台 pwsh 退出会收割子进程（分片扫描必须用常驻驱动）；`ProcessPoolExecutor` 不可用；控制台 GBK 会让含 `−`(U+2212)/`×` 的输出崩（脚本里用 `sys.stdout.reconfigure(encoding='utf-8')`）。
6. **裸 `python` 是 Python 2.7.11，不是本项目的解释器**：`where.exe python` 首个命中是 `D:\python.exe` = **Python 2.7.11**（`python -V` 在本轮实测如此），后面才是 `miniconda3\python.exe`(3.9) 与 `Python38-32`。本轮写检查脚本时就因此吃到 `AttributeError: 'file' object has no attribute 'reconfigure'`。**本项目解释器是 `D:\aconade\python.exe`（Python 3.13.9）；所有命令、复现步骤里的 `python` 都应当写成这个绝对路径。** 顺带：`MPLBACKEND=Agg`，读 PDF 用 `C:\texlive\2025\bin\windows\pdftotext.EXE -q -enc UTF-8 <pdf> -`。
7. **这份交接文档自己的缺陷（写完后复核才发现，已在本版修掉，记录以免再犯）**：
   - §2C 曾把**两个口径混进一句话**（4 点中位数比值 `3.084×` 并排配上整数步的 `dlnL/dn` 1.252），而源文件 §5 明文禁止混用 —— **违反自己引用的纪律，且位置在最承重的 n=6 论证上**。凡"同一现象有两种口径"，务必分开写、各自标注；
   - §1 曾漏钉 `clock`，把当初 `n_A1_gate` 那类静默回落的坑原样复制了一遍 —— 教训是：**列工作点时，凡是"模型有默认值、而默认值不等于冻结值"的字段都要显式列出并注明两条路径**（`clock` 就是第二个）；
   - §4 曾把 `drop 4` vs `drop 8` 误判为"不一致、要么统一"，实际上那是两套不同模数的计数器各丢一个周期，真去统一会弄坏其中一个 —— **看到两个不同的常数，先确认它们是不是各自单位的同一个量**；
   - §2G 曾声称"10 条差异"却只列 5 条且无落盘产物 —— **计数型断言必须配可核对的产物，否则只报承重的那几条**。
8. **"列保真"不等于"列语义正确"（Figure 4-3 的真实 bug，见 §2F）**：图件自检原本只验"画出来的曲线 == CSV 里名叫 `S2` 的列"，而那个列里装的其实是 `F1` —— 于是 n=4 面板画 F1、n=6 面板画真 S2，**两条不同物理量被放在同一个轴标签下**，还配了"S2 峰值 2.25、超出物理区间"的说明。**教训：验证必须分两层——(a) 曲线 == 列（保真）、(b) 列 == 模型状态（语义）。** 第二层现在由 `dshwork/figures/check_state_column_semantics.py` 承担：按模型 `state_names` 定位索引、要求"声称是分数"的列必须落在 `[0,1]`、用 `g1` 恒等式做**判别性**检验（把 A1/F1 对调必须失败，实测差 0.389）、并要求该列能**重放出读窗的位标签**（n=4 重放 56 个窗 0 不符）。**最便宜的那道闸（DNA 构象分数必须在 [0,1]）本来一行就能抓住它，而且不需要重新积分。**
9. **图注里的机理断言必须机检**（收尾轮）：Figure 4-3 新图注说「首次进入高态后反向写入不充分 / 未重新进入低态带」。原来的 n=4 数据文件**没有** `J_rev2`，这句话当时**不可验证**。修法是让取数脚本把 `b2_I/b2_R/b2_C/J_fwd2/J_rev2` 与 `excursion` 块一起写出（首次进高带时刻、其后最小值、是否回低带、正/反向通量积分、穿越与门时段表），并**在断言不成立时直接报错退出**。实测：首次进高带 **t = 33.13 h**，其后 17 007 个采样点最小值 **0.3087 > 0.30**，确实从不回低带。**写进图的话，必须有一条代码能让它为假。**
10. **交付物的数值必须与产物一致（v1 的"3"）**：v1 报告写「3 个点的数字读出完美却未认证」，而产物在两个清晰定义下都给 **4** 且成员不同。**不可复现的计数要么给出定义、要么不要写**；收尾轮选择如实报告这个矛盾（列出 A/B 两式），而不是去凑那个 3。
