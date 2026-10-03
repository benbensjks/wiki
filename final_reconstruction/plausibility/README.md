# 参数合理性检验包（`plausibility/`）

> 状态（2026-09-20 更新）：**四项 checks 与主扫描（44 点 × 600 h）已跑完**。见"实测结论"一节。
> 这里只放"除了继续扫参数之外，还需要检验什么"的那部分；单点机理验证与探针由另一条线负责。
> 所有脚本只写本目录，**不碰** `model.py`、`model_twobit34.py`、`model_threebit51.py`、
> `verify_threebit51.py`、`verify_twobit_causal.py`、`selected_profile.json` 以及任何其他 agent 的产物。

## 实测结论（2026-09-20）

**加长 vs 加陡：加陡胜出（R3）。** 44 点 600 h 同判据扫描：

- **臂 B（独立 `n_A1_gate`）20 点里 16 点通过认证**，`n_A1_gate ∈ {5,6,7,8}` 全部给出
  `crossings = gates = reverse_events = 14`、冷启动读窗 `1234567012345670…`（56 拍 / 7 个完整模 8 周期）；
  **在完全冻结的成熟时间 32.5 min 上就成立**。
- **臂 A（加长成熟时间）24 点全部失败**；`mat ≥ 90` 时穿越数冲到 22–27（远大于 14）——剂量过大的振铃；
  泄漏剂量在整条轴上恒为 0.62–0.85，与成熟时间无关。
- 泄漏随门指数单调改善但**饱和**：n=5 约 0.22、n=8 约 0.17。因此**单纯加陡到不了 0.10 那条线**，
  残余泄漏另有来源。

**R1 的实际阻塞点是泄漏线，不是事件计数。** 16 个认证点的 `leak_worst_Jrev2` 在 0.172–0.365，
没有一个 ≤ 0.10。`crossings == gates == 14` 在 600 h 里是可达的（本包对事件收集没有烧入限制；
烧入只用于因果配对的后段窗口）。

## 为什么需要这一包

当前 51 状态三级模型的状况是：前 34 状态冻结且已认证；bit2 在第 1 轮 35 点扫描中 **0/35 认证**，
最优 12/14 穿越，且最优解落在成熟时间网格边界。此时有两种可能：

1. **参数没调够**（剂量不足）→ 继续扫；
2. **参数本身不合理或结构上受限**（地板、简并、量级、崩塌边界）→ 再扫也白跑。

这一包的五个脚本就是为了在花掉下一轮 30–60 分钟扫描之前，先把第 2 类问题排除或坐实。

## 脚本清单

| 脚本 | 回答的问题 | 决定什么 | 代价 |
|---|---|---|---|
| `check_threshold_margin.py` | 模型里每一个 Hill 项 `H(x;K,n)`，工作轨迹是否压在 `x/K≈1` 上？ | 找出还有哪些"刀刃上的阈值"。已知 `A1/K_A1` 就踩在阈值上（残余 0.303、峰值 1.24、K=1.2），这正是上一轮改指数能奏效的原因 | 1 次 300 h 仿真 |
| `check_absolute_scale.py` | `1 a.u. = uM_per_au µM` 下，各池的实际浓度、copies/cell 是多少？总翻译负荷占细胞预算多少？ | 决定 `uM_per_au` 能不能被称作"标定值"，以及这个设计在表达上是否可行 | 1 次 300 h 仿真 |
| `check_mass_balance.py` | 显式复合物对游离 I/RDF 的消耗是否实现正确、量级多大？ | 这是相对曾墨涵约化模型的**结构性添加**，必须量化并写进报告 | 1 次 300 h 仿真 |
| `check_identifiability.py` | (a) `kon/koff/δ_C` 是否只有 `q` 可辨识？(b) 表达链的 (mRNA, 成熟) 是否只有"有效延迟"一个自由度？ | 决定报告里能声称标定了几个参数；避免把简并自由度当成独立标定 | 6 + 4 次 100 h 仿真 |
| `check_rdf_collapse.py` | bit2 停滞是不是因为 RDF 池被自己的 Hill 崩塌掏空（结构性），而不是剂量不足？ | **直接决定"继续延长成熟时间"这条路要不要走** | 3 次 600 h 仿真 |
| `scan_sharpen_vs_lengthen.py` | **加长脉冲** vs **加陡门**，同一判据下谁有效？ | 下一轮扫描的主决策；两条臂都失败则转结构 | 44 点 × 600 h（分片） |

## 预先写死的判据（不得事后更改）

所有阈值都是脚本顶部的常量：

- 阈值余量（v2 规则）：**主指标** `frac_time_near_threshold ≥ 0.25`；**次指标** 仅对 `role` 以 `gate` 开头的项、
  `ratio_median > 0.10`。`straddles_threshold` 在 v2 里**只作描述、不参与判定**——开关型项按定义就会跨越阈值，
  v1 把它放进 OR 导致 16/18 全中、聚合标签失效。可用 `--rescore-from-csv` 对已存 CSV 重打分，不必重跑仿真。
- **绝对量级**：单个游离池 > `30 µM` → FLAG；总合成需求 > 细胞翻译预算 `4×10⁶ copies/cell/h` 的 `10%` → FLAG。
  （两个参考值都是数量级假设，输出里会连同结果一起打印，便于读者重新标定。）
- **质量守恒**：`d(I+C)/dt` 与 `d(R+C)/dt` 的最大残差 > 该项最大量级的 `1e-3` → 实现错误。
- **可辨识性**：同 `q` 不同 `(kon,koff,δ_C)` 之间 `max|ΔS| < 0.02` → 简并；
  同时跑一个**故意错配 q** 的对照，若对照位移 < 0.10 则判定"该测试不灵敏"，不下结论。
- **RDF 崩塌**：若某个探针（软化抑制 / 移除抑制）能把崩塌点移到停滞点以下**并使 bit2 恢复认证**，
  则崩塌是承重结构，停止剂量扫描；反之剂量解释成立。
- **加长 vs 加陡**：
  - `R1 成功`：存在点满足 `certified=True` 且 `crossings == gates == 14` 且 `leak_worst_Jrev2 ≤ 0.10`；
  - `R2 饱和`：臂 A 固定 mRNA 时 `crossings(mat=120) ≤ crossings(mat=60)` → 延长时间不再有用；
  - `R3 加陡胜出`：臂 B 在**冻结的 32.5 min** 上就能认证，而臂 A 在任何成熟时间都不能；
  - `R4 双臂皆败`：两臂都没有认证点 → 杠杆不在参数上，转结构。

## 新增的泄漏指标（替换 off/on 峰值比）

`plausibility_common.leak_budget()` 计算**两次 carry 之间**的反向剂量积分 `∫J_rev2`。
理由是复合物 Hill 在低驱动下是**平方律**，真正侵蚀存储态的是积分剂量而不是峰值比。
参考量级：要把 21 h 驻留期的侵蚀压到 0.1 以下，需要 `g1_off ≲ 3e-4`，
对应 off/on 峰值比约 `2e-3` —— 比现在沿用的 `0.10` 严约 50 倍。
本包的输出同时给出两种口径，便于对照，但**判定用积分剂量**。

## 解耦要求（臂 B 的硬前置条件）

`model_threebit51` 里 `ZENG['n_A'][1]` 被**两处**使用：

1. `model.Model.carry_promoters` → `g1 = H(A1;K_A1,n)·G(F1;K_F1,n_F1)·clock`，这是我们要**加陡的门臂**；
2. `ThreeBit51Model._carry1_sources` → `source_F1 = alpha_F1·H(A1;K_A1,n)`，这是**不相干臂的驱动，必须保持冻结值**。

所以直接 patch `ZENG['n_A'][1]` 会同时移动两者，测出来的效应无法归因。

`plausibility_common.DecoupledThreeBit51Model` 用一个**子类覆写**把两者分开：门臂用 patch 后的值，
`_carry1_sources` 与 `initial_state` 里的 F1 产生/预载锁定在 `FROZEN_N_A1 = 4.0`。**不修改任何冻结文件。**

`scan_sharpen_vs_lengthen.py` 在 `--arm B` 时会先读 `decoupling_verification.json`；
若该产物缺失、未通过、或与当前实现哈希不符（`implementation_sha256` / `source_sha256`），
**直接 `SystemExit` 拒绝启动**，不会跑出一个不可归因的结果。

先运行验证、再扫描：

```powershell
& 'D:\aconade\python.exe' .\scan_sharpen_vs_lengthen.py verify-decoupling --hours 100
# 通过后才会写出 decoupling_verification.json，臂 B 才放行
```

验证做两件事：

- **测试 1**：门指数取冻结值 4.0 时，解耦模型必须与冻结模型**逐位相同**（初值 gap 应为 0，轨迹 gap < 1e-9）——证明解耦没有引入偏差；
- **测试 2**：把门指数移到 6 / 8 后，在同一条冻结轨迹上，**F1 产生项必须与冻结值完全相等**，
  而**共享槽位的旧做法必须使它改变**——证明解耦确实把两个用途分开了。

两项都通过才写产物。

## 运行方法

```powershell
cd final_reconstruction\plausibility

# 1) 五个便宜检查（每个 1–3 次仿真，逐个串行，日志写到 log_*.txt）
powershell -ExecutionPolicy Bypass -File .\run_plausibility.ps1 checks

# 2) 主扫描（分片，文件系统完成检测；本沙箱禁用 Wait-Process）
powershell -ExecutionPolicy Bypass -File .\run_plausibility.ps1 scan -Nshards 12 -Arm AB -Hours 600

# 3) 汇总 + 预登记判据 + 投影护栏
powershell -ExecutionPolicy Bypass -File .\run_plausibility.ps1 merge -Nshards 12

# 4) 刷新 SHA256 清单
powershell -ExecutionPolicy Bypass -File .\run_plausibility.ps1 manifest
```

也可以单独跑：

```powershell
& 'D:\aconade\python.exe' .\check_threshold_margin.py --hours 300 --no-write
& 'D:\aconade\python.exe' .\scan_sharpen_vs_lengthen.py scan --arm B --shard 0 --nshards 12 --hours 600
& 'D:\aconade\python.exe' .\scan_sharpen_vs_lengthen.py merge --nshards 12 --projection-guard
```

## 关于时长的一条预登记说明

32 状态验收要求"丢弃前 8 个读窗后仍有 ≥16 个有效读窗"，而一个读窗 ≈ 10.58 h，
所以 **稳态判定至少需要 24 个读窗 ≈ 254 h**。

- `--hours 300` 可以给出稳态判定，可用于筛选；
- **`R1 成功` 只承认 600 h**（56 个读窗 ≈ 7 个完整模 8 周期），与既定的验收长度一致；
- `--hours 200` 只能看穿越数与泄漏剂量，**不能**作为认证依据。

## 两条诊断性说明（务必连同结果一起读）

1. **臂 B 的实现在生产模型里必须改。** 本包用运行时 monkey-patch `model.ZENG['n_A'][1]` 来诊断，
   而这个槽位同时被"门里的 A1 臂"和"驱动 F1 的 Hill"使用。在工作点上第二处是惰性的
   （n=4 时 F1 稳态 0.018、n=6 时 0.0012，都远低于 `K_F1=0.4`），所以观察到的效应应归于门里的 A1 臂；
   但正式实现必须把它暴露成一个**独立的新参数**（例如 `n_A1_gate`），并在生产模型里重跑同一比较。
   这属于新增接口假设，不是修改曾墨涵的表。
2. **这些检查回答的是"参数是否合理"，不是"计数器是否成立"。**
   它们不替代四初态、严格容差、扰动扫描，也不构成 γ 语义或 a.u. 标定的实验依据。

## 随机性预审的收尾（2026-09-26）

复核要求把单一 `certified` 布尔拆成分项，并收紧几处过度解读。本轮**没有重新积分任何轨迹**。

| 文件 | 作用 |
|---|---|
| `PREAUDIT_STOCHASTICITY_REPORT.md` | 预审报告 **v2**。改了 6 处：小数字改为「门相 最小 20 / p05 38.6 / **中位 2160**」并声明**尚未确定**主要写入是否发生在低拷贝时段；`1/√N` 降级为**泊松参考量级**而非已算出的 CV；补上「改 `uM_per_au` 会改变动力学，不能只按比例换算」；结论由「随机性研究是必需的」改为「**值得进一步评估随机效应**」；旋钮名改为「`uM_per_au` 接口参数扰动」；指向新的分项验收表 |
| `make_acceptance_breakdown.py` | 生成器。**纯后处理**：只读 `preaudit_stochasticity_all.csv`，不建模型、不积分。与 verdict 的 `certification_table`/`leak_table`/`off_on_gate_peak_ratio_table`/`failures_by_arm` **逐点交叉核对** |
| `preaudit_acceptance_breakdown.{md,json,csv}` | 21 点分项验收表：历史 `certified`（原值原定义）、数字计数（码串 / 未标注窗 / 模 8 递增）、进位因果（一一对应 / 因果顺序 / 方向交替）、门对比度（原始 off/on 比值 + 是否超 0.10）、配对泄漏（`L_symmetric` + 切尾档 + 输出网格）、时序与读窗（setup/hold、commitment） |
| `archive_preaudit_report.py` + `preaudit_report_archive/` | v1 报告与**两套旧清单**的逐字归档 + `ARCHIVED_HASHES.json`（存档前后哈希必须相等，否则报错退出） |

**新增标签 `counting_and_causality_passed`**，定义恰为 `steady ∧ one ∧ order ∧ alt`。
它**不是** `certified` 的改名，也**不附带任何泄漏阈值** —— H4 的配对 `L_symmetric` 只作数字列出。
实测：历史 `certified` 通过 **14/21**，新标签通过 **18/21**；4 个点是 **contrast 臂单独失败**
（计数与因果全过），其中 `conc_scale|1.05`/`|1.10` 的泄漏（0.010691 / 0.009785）**优于冻结点 0.014218**。

**离散通过点，不是连续可靠范围**：`uM_per_au` 那一维在历史谓词下通过 {0.95, 1.00}、
在计数+因果下通过 {0.95, 1.00, 1.05, 1.10}；`clock_K` 全 7 点；`K_A1` 分别为前 5 点 / 全 7 点。
两者都**只是已测网格点**，未采样因子从未运行。另：`uM_per_au` 的**失败原因随因子变化**
（< 0.95 失败在 steady/causality，> 1.00 只失败在 contrast），读成"某区间可靠"是不对的。

**一处如实报告的不一致**：v1 写的「3 个点的数字读出完美却未认证」**无法从产物复现**；
两个清晰定义都给 4 且成员不同（分项表 §1 有 A/B 两式与差异点）。分项表只报这两个可复现说法。

**认证规则的版本纪律**：历史结果**保留原定义**；要调整认证谓词（例如按 H4 配对泄漏口径重定
contrast 臂），必须**建立新版本、并列重评**，**不能删掉失败项后沿用原认证名称**。

## 交付物

运行后本目录会出现：各检查的 CSV/JSON、`sharpen_vs_lengthen_all.csv` 与 `_verdict.json`、
分片 CSV 与 meta（含源文件哈希）、`SHA256SUMS.json/.txt`、以及各次运行的日志。

## 两份必读文档（不在本目录的清单里，各有自己的清单）

| 文档 | 内容 | 谁在管它的哈希 |
|---|---|---|
| `plausibility/HARD_CONSTRAINTS.md` | H1–H8 硬约束全文与故障记录 | 本目录 `SHA256SUMS.json` |
| `plausibility/N6_FREEZE_JUSTIFICATION.md` | 为什么冻结在 n=6 的正式论证（含它不声称什么） | 本目录 `SHA256SUMS.json` |
| `plausibility/RDF_PROBE_DESIGN.md` | RDF 承重性探针的预登记设计（尚未实现） | 本目录 `SHA256SUMS.json` |
| `plausibility/PREAUDIT_STOCHASTICITY_REPORT.md` | 随机性预审报告 **v2**（含 v1 的归档指针与 6 处修订） | 本目录 `SHA256SUMS.json` |
| `plausibility/preaudit_acceptance_breakdown.md` | 21 点分项验收表（纯后处理；唯一可引用的分项口径） | 本目录 `SHA256SUMS.json` |
| `../CURRENT_BASELINE.md` | 项目权威基线；**§11「判据与验证纪律」**是 H1–H8 的交付版 | wiki 根 `SHA256SUMS_authority.json` |

`SHA256SUMS.json` 只覆盖 `plausibility/**` 与少数具名源文件，**不覆盖 wiki 根目录的权威文档**。
后者由根级清单覆盖，用下面这条命令生成并校验：

```powershell
& 'D:\aconade\python.exe' .\make_authority_manifest.py write   # 重新生成根级清单
& 'D:\aconade\python.exe' .\make_authority_manifest.py check   # 只读校验（含本目录清单是否过期）
```

`check` 会重算本目录 `SHA256SUMS.json` 里的每一条，因此**任何本目录内的改动都会让根级清单的
校验失败**，必须重新 `write`。这是有意的：清单过期必须报错，不能静默通过。
