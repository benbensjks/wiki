# 扫描图：图注与参数（`scan_figures/`）

本目录把**已经跑完**的扫描结果画成四张图。**没有做任何新的仿真**：每张图都只是读取既有产物，
输入清单（含每个文件的 SHA256 与行数）写在 `figure_data_manifest.json` 里，
产物本身的 SHA256 / 像素尺寸也记在同一份清单里。

```powershell
& 'D:\aconade\python.exe' .\plot_scan_figures.py        # 生成四张图（可重复：两次运行逐字节相同）
& 'D:\aconade\python.exe' .\selfcheck_scan_figures.py    # 47 项数值断言，全部通过
```

## 0. 必须先知道的三件事

1. **认证 ≠ 计数正确。** `certified` 是冻结验证器的裁决：
   `steady ∧ one ∧ order ∧ alt ∧ contrast`。其中 `contrast` 用的是 off/on 门峰值比
   （限值 `MAX_OFF_ON_GATE_RATIO = 0.10`），而 H4 早已把这条指标从**泄漏量化**里降级。
   图上绿/红只表示"该裁决通过与否"，不表示"数得对不对"。
2. **LEGACY 面板不是冻结结果。** 51 状态的早期扫描（carry1 首轮、A1/F1 分开成熟）
   跑在 `ZENG['n_A'][1] = 4.0` 的**共享槽位**上，早于 `n_A1_gate` 参数出现，
   **不是**冻结的 `n_A1_gate = 6.0`。图里这些面板标题带 LEGACY，禁止当作冻结配置引用。
3. **泄漏一律是配对 R/F 超周期 `L_symmetric`**（H4），2 min 输出网格，每张图注明切尾档。

## 1. `fig_2bit_single_param.png` — 2bit 单参数扫描

**数据**：`twobit34_results/robustness/robustness_all.csv`（75 点 = 3 mRNA × 5 成熟 × 5 uM_per_au，600 h）。

| 面板 | 内容 |
|---|---|
| (a) | 固定 成熟 32.5 / uM 5.75，扫 carry mRNA（1/2/4 min）：bit1 setup 与 hold 裕量 |
| (b) | 固定 mRNA 2 / uM 5.75，扫成熟（25–35 min） |
| (c) | 固定 mRNA 2 / 成熟 32.5，扫 uM_per_au（5.5–6.5） |
| (d) | **边际稳健性**：每个参数取值的认证比例（另两维遍历） |

**可直接使用的图注**：
> 图 X. 两比特计数器（34 状态）的单参数切片与边际稳健性。75 点联合网格（3 种 carry mRNA
> 半衰期 × 5 种成熟半衰期 × 5 种 uM 换算，各 600 h，有限读窗判据），共 53/75 认证。
> 三个切片均过选定的中心点（mRNA = 2 min，成熟 = 32.5 min，uM_per_au = 5.75）；
> 圆点 = 认证，红叉 = 未认证。

**读图保留意见**：
- 失败原因**只有一种**：22 个未认证点全部是 `readout_unlabelled_window`（读窗内有未定值样本）。
- (d) 的 uM 依赖**非单调**：6.00 → 15/15、6.25 → **4/15**、6.50 → 15/15。
  这与 `strict_tolerance` 里 `pass_notch_lip` / `fail_notch_mid` / `fail_notch_low` 的命名一致，
  是真实存在的"凹口"，**不要**读成单调趋势，也不要用单调插值去外推。

## 2. `fig_3bit_single_param.png` — 3bit 单参数扫描（LEGACY 与 FROZEN 并列）

**数据**：`threebit51_results/carry1_scan_dsh/carry1_all.csv`（35 点，LEGACY）、
`carry1_split_scan/split_all.csv`（63 点，LEGACY）、`plausibility/nA1_vs_maturation_all.csv`（25 点，FROZEN）、
`plausibility/preaudit_stochasticity_all.csv`（21 点，FROZEN）。

| 面板 | 工作点 | 内容 |
|---|---|---|
| (a) | **LEGACY n_A[1]=4.0** | 固定成熟 60 min，扫 carry1 mRNA（0.5–8 min）： crossings/cycles 与门峰值 |
| (b) | **LEGACY n_A[1]=4.0** | 固定 mRNA 4 min，扫成熟（5–60 min）：crossings/cycles 与门剂量 |
| (c) | FROZEN n_A1_gate=6 | 固定 A1 成熟 32.5，扫 `n_A1_gate`（5–7）：进位次数 + 配对泄漏 |
| (d) | FROZEN n_A1_gate=6 | 固定 n=6，扫 A1 成熟（25–45）：进位次数 + 门剂量 |
| (e) | FROZEN | 全局浓度增益 f × uM_per_au（7 档）：认证 + off/on 峰值比（限值 0.10） |
| (f) | FROZEN | 时钟阈值 `clock_K` 与门阈值 `K_A1`（各 7 档） |

**可直接使用的图注**：
> 图 Y. 三比特计数器（51 状态）的单参数扫描。上排为**早期共享槽位指数**（ZENG n_A[1] = 4.0）
> 的 carry1 mRNA 与成熟时间扫描：35 点与 63 点**无一认证**，bit2 从不计数（63 点全部
> `gate_not_formed`），这是"停滞"期证据，不是冻结配置的结果。下排为**冻结工作点**
> （`n_A1_gate` = 6.0，构造器路径）：门指数 n=4 直接失败、n≥5 计数；A1 成熟时间 5 档全部认证
> 且剂量变化 <5 %、配对泄漏变化 ≤2.7 %；全局浓度增益只在 f = 0.95 与 1.00 认证。

**读图保留意见**：
- (b) 的门剂量曲线在成熟 ≤ 5 min 处**断开**：该点 `gate_events = 0`（门根本没形成），
  剂量不可测。断口是"未测到"，不是零——**不要**把折线连过去。
- (c)(d) 的圆点/红叉是冻结验证器的裁决；泄漏画在右侧对数轴，**5 % 切尾**。
- (e)(f) 的 off/on 峰值比与限值 0.10 的关系解释见 `PREAUDIT_STOCHASTICITY_REPORT.md`：
  f=1.05 与 K_A1 f=1.10 是**只有 contrast 臂失败**，其数字读出完全正确。

## 3. `fig_2bit_robustness.png` — 2bit 鲁棒性分析

**数据**：75 点联合网格、`strict_tolerance_all.csv`（9 点 × 3 组求解器设置）、
`perturbations_{s,rdf,affl,int1}_all.csv`（共 140 行 = 4 组扰动 × 4 初态）。

| 面板 | 内容 |
|---|---|
| (a)(b)(c) | 认证地图：uM × 成熟，分别固定 mRNA = 1 / 2 / 4；金星 = 选定中心 |
| (d) | mRNA = 2 的 `timing_margin_h` 热图，白圈 = 认证、红叉 = 未认证 |
| (e) | 严格求解器容差复核（9 点名 × 3 设置） |
| (f) | 分子池扰动的**结局分布**（5 类，堆叠） |

**可直接使用的图注**：
> 图 Z. 两比特模型的鲁棒性。(a–c) 三种 carry mRNA 下的认证地图；(d) 时序裕量热图；
> (e) 九个代表点在 rtol/atol/max_step 三组严格设置下的复核——**没有任何裁决翻转**
> （4 个通过点保持通过、4 个失败点保持失败）；(f) 分子池扰动（S0/S1、RDF0/RDF1、A0/F0、Int1，
> 各 4 个初态共 140 行）的结局分布：同相位 125、相位偏移 +1/+2/+3 共 8、计数丢失 7。

**读图保留意见（重要）**：
- (f) 里**橙色不是失败**：`stable_phase_shift_*` 的行 `certified = True`，计数器仍在正常计数，
  只是相对无扰动运行**相位偏了 1–3 个读窗**。把这五类压成"恢复/未恢复"会把这一点抹掉
  （本目录第一版正是这么画的，已改）。
- 计数丢失集中在 S 态扰动（`flip` 4 行里 1 行、`delta_-0.20` 8 行里 4 行）与 A0/F0 池 0.5 倍
  （8 行里 2 行）；**RDF 池与 Int1 池扰动 60 行全部同相位恢复**。

## 4. `fig_3bit_robustness.png` — 3bit 鲁棒性分析

**数据**：`nA1_vs_maturation_all.csv`（25 点）、`carry_pairing_grid_verdict.json`（20 点，600 h）、
`preaudit_stochasticity_verdict.json`、`preaudit_stochasticity_inventory.json`、
`postfreeze_peak_perturbation/postfreeze_all.csv`（36 行）。

| 面板 | 内容 |
|---|---|
| (a) | 冻结 25 点网格（n × A1 成熟）认证地图，金星 = 冻结工作点 |
| (b) | **泄漏膝点**：逐 n 中位数，1 %/5 %/10 % 三档切尾 |
| (c) | 20 点网格里四种 (mRNA, 成熟) 组合的逐 n 曲线 |
| (d) | 静态容差对称宽度（预登记规则） |
| (e) | 未认证点**失败在哪条臂** |
| (f) | 拷贝数盘点：决策相关池（门相最小 / 全周期谷值），含 100 拷贝小数字线 |

**可直接使用的图注**：
> 图 W. 冻结三比特工作点（n_A1_gate = 6.0）的鲁棒性。(a) 25 点（n × A1 成熟）全部认证，
> 说明决定因素是门指数而非成熟时间；(b) 配对 R/F 泄漏随 n 单调下降且排序 8<7<6<5<4 在三档
> 切尾下一致，5→6 改善 3.08 倍、6→7 降到 1.22 倍、7→8 仅 1.06 倍（膝点在 6）；
> (c) 该排序不是单一设定的产物，四种 (mRNA, 成熟) 组合都成立；(d) 静态容差最差的是全局浓度增益
> （宽度 0.00，只有 f = 0.95/1.00 认证）；(e) 7 个未认证点中 4 个是 **contrast 臂单独失败**；
> (f) 执行写入的 Int2 池在门相最低约 20 个拷贝，F1 谷值约 40 个。

**读图保留意见**：
- (d) 的红色阈值线是**预登记**的 SMALL 门限 0.05；但 (e) 说明该结论的成因是认证谓词里的
  contrast 臂，不是计数失效——引用前务必读 `PREAUDIT_STOCHASTICITY_REPORT.md` §1。
- (f) 的**灰色斜纹 = 结构零**（该相位 < 1 拷贝，是"没有"而不是"很少"），
  为落在对数轴上统一截到 1e-3；蓝色才是可用于噪声估计的水平。
  拷贝数换算线性依赖未标定的 `uM_per_au = 5.75`。

## 5. 验证方式（本文档的诚实边界）

**这些图我没有目视检查过** —— 当前模型不支持图像输入。替代方案是数值验证：
`selfcheck_scan_figures.py` 在进程内重新渲染每张图，然后断言：面板数与系列数、
**抽样数值与源产物逐位一致**（例如 (b) 的 5 % 中位数、b2_I 的 20.028 拷贝、
2bit 中心切片的 setup 裕量）、缺失数据的**位置**与源文件一致、LEGACY/FROZEN 标注齐全、
无缺字警告、PNG 真实且尺寸合理、以及清单里每个输入哈希仍然匹配。47 项全部通过，
且两次运行产出的 PNG **逐字节相同**。

**但数值正确不等于排版好看**：如果渲染出来的图有标签重叠、文字压线之类的问题，
需要看图的人（你）指出，我按反馈调布局。

## 6. 不允许的用法

- 把 LEGACY 面板当成冻结工作点的结果（它与 `n_A1_gate = 6.0` 无关）。
- 把 `certified = False` 直接写成"计数器失效"（4 个点是 contrast 臂单独失败）。
- 把 (f) 的灰色斜纹条当成"测得 1e-3 个拷贝"。
- 把不同切尾档或不同输出网格的泄漏值放在一起比较（H4）。
