# 机制图 / 流程图的代码生成

对应 `04_第3章_作图说明.md` 与 `05_第4章_作图说明.md` 中**可以用代码绘制**的示意图。
数据图（读数条带、热图、剂量—响应）不在本目录，它们需要仿真输出。

## 运行

```powershell
$env:MPLBACKEND='Agg'
& 'D:\aconade\python.exe' .\make_schematics.py --out .\out
# 只画其中一张：
& 'D:\aconade\python.exe' .\make_schematics.py --out .\out --only fig04_ffl_carry
```

依赖只有 `matplotlib`（`Agg` 后端，无需图形界面）与 `numpy`。

## 生成的文件

| 图 | 文件名 | 类型 | 来源 |
|---|---|---|---|
| Figure 1 | `fig01_modeling_workflow.png` | 流程图 | `01` 图 1：Design → Build → Test → Learn → Decision 闭环 |
| Figure 3-1 (旧版构图) | `fig03_single_bit_mechanism.png` | 机制图 | `04` §3：Int / RDF / BM3R1 开关 |
| **Figure 3-1（按精简版要求）** | `fig03_1_single_bit_memory.png` / `.svg` | 机制图 | `精简交付版/02` §Figure 3-1：左 PB(0) / 右 LR(1) + 下方一组调控关系 |
| **Figure 3-2（按精简版要求）** | `fig03_2_switching_read_windows.png` / `.svg` | 数据图 | `精简交付版/02` §Figure 3-2：两个连续时钟周期、三行共用时间轴、真实轨迹 |
| **Figure 4-1（按精简版要求）** | `fig04_1_carry_architecture.png` / `.svg` | 机制图 | `精简交付版/02` §Figure 4-1：整体链 + 两级 carry 并排（含时钟 AND 门） |
| **Figure 4-2（按精简版要求）** | `fig04_2_mod4_mod8.png` / `.svg` | 数据图 | `精简交付版/02` §Figure 4-2：两位/三位读数条带（28 / 56 reads） |
| **Figure 4-2a（新增，2026-09-27）** | `fig04_2a_twobit_under_upstream.png` / `.svg` / `.pdf` | 数据图 | **真实韩氏振荡器在振荡 + 两位计数器在完整 mod 4 计数**（同一时间轴三面板）；见下节 |
| **Figure 4-3（按精简版要求）** | `fig04_3_sharper_gate.png` / `.svg` | 数据图 | `精简交付版/02` §Figure 4-3：n=4 vs n=6 的 S2；五个指数的配对泄漏 |
| **Figure 4-4（按精简版要求）** | `fig04_4_validation_summary.png` / `.svg` | 验证摘要图 | `精简交付版/02` §Figure 4-4：八初态条带 + 扰动战役汇总 |
| Figure 4-1 | `fig04_ffl_carry.png` | 机制图 | `05` §3：I1-FFL 进位门（含两个时间过程小图） |
| Figure 4-6 | `fig04_6_bm3r1_shutdown_interface.png` | 机制图 | `05` §8：BM3R1 共享接口 |
| Figure 7 | `fig07_full_system_coupling.png` | 框架图 | `01` 图 7：全系统耦合 |
| Figure 8 | `fig08_design_panel.png` | 面板 | `01` 图 8：Modeling-Guided Design |

## 两张 Figure 3-1 的关系（不要混用）

`fig03_single_bit_mechanism.png` 与 `fig03_1_single_bit_memory.png` **不是同一张图的修订版**，
而是两套构图：

| | 旧版 `fig03_single_bit_mechanism` | 新版 `fig03_1_single_bit_memory` |
|---|---|---|
| 依据 | `04_第3章_作图说明.md` §3 | `精简交付版/02_作图要求精简版.md` §Figure 3-1 |
| DNA 构型 | 一个合并框 `DNA state / PB (1−S)  LR (S)` | **左右两张卡** `PB / 0`、`LR / 1` |
| 往返箭头 | `Forward recombination (PB→LR)` / `Reverse recombination (LR→PB)` | `Int` + `Forward switching` / `Int–RDF complex` + `Reverse switching` |
| 下方调控关系 | 分散在四周 | **集中成一个虚线块**（PB→BM3R1、BM3R1 ─┤ RDF expression、LR→RDF expression、Int + RDF ⇌ Int–RDF complex） |
| RDF 池命名 | `RDF pool (memory)` | `RDF`（**不写 memory**，按要求） |
| 图上文字 | 含脚注 `Memory is stored in the RDF pool, not in Int` | 无任何 memory 字样 |

**新版没有改动旧版**，旧版文件与 `make_schematics.py` 原样保留，便于对照。
新图的生成脚本是 `make_fig03_1_single_bit_memory.py`（独立文件，只复用 `tempo_style.py` 的配色与线型原语）。

## 新版 Figure 3-1 的自检（因为出图的人看不到图）

`make_fig03_1_single_bit_memory.py` 画完即自检，**40 项**，覆盖：

1. 要求里 10 个字符串**逐字**作为独立文本元素存在（含 en dash 的 `Int–RDF complex`）；
2. 禁止措辞（`memor*`、`epigenetic`）在图内与导出文件里都不存在；
3. PB 卡在 LR 卡左侧、两卡同高、往返箭头方向正确（从实际创建的 artist 读回几何，不复述输入）；
4. `BM3R1` 在它所抑制的 `RDF expression` 箭头的左侧，`RDF` 框在该箭头正下方；
5. 所有文本的渲染包围盒**两两不重叠**且不越出画布（这是"看不到图"时唯一能替代目视的检查）；
6. 导出的 SVG 用真 `<text>`（`svg.fonttype='none'`），因此**标签可编辑**且可被逐字校验。

SVG 保持文本（而不是把字转成轮廓）是刻意的：既方便美工改字，也让上面的逐字校验成为可能；
代价是 SVG 依赖 `DejaVu Sans`；若交付机器没有该字体，文字会重排 —— 以 PNG 为准。

## Figure 3-2 — Switching and read windows（数据图）

脚本 `make_fig03_2_switching_read_windows.py`，**不含任何仿真**：全部读取既有冻结产物。

**数据来源（只读，含哈希）**

| 文件 | 大小 | SHA256（前 16 位） | 提供什么 |
|---|---|---|---|
| `final_reconstruction/twobit34_results/trajectories.csv` | 7 553 517 B | `805D31CF82D5D005` | 300 h、2 min 网格的 34 状态轨迹；`S0`、`b0_R`、`b0_C`、`J_fwd0`、`J_rev0`、`C31_flux_au_h` |
| `final_reconstruction/twobit34_results/read_windows.csv` | 2 944 B | `C2B1A09B28CC7FF3` | 验证器自己算出的 28 个读窗（起始/结束/谷值/位值） |

**截取哪两个周期，为什么**：取**最后两个完整周期**（cycle 25、26，谷值 271.833 与 282.433 h），
时间轴 `[266.533, 287.733] h` = **恰好 2 个周期**（10.600 h × 2）。选这一对是因为：
两次翻转（S0 = 0.5）都**落在画面内部**（276.775、286.631 h），两个读窗都在范围内，
而且 cycle 25 读 **0**、cycle 26 读 **1**，因此低带与高带各出现一次。
两端都落在数据 `[0, 300] h` 之内 —— **末端不需要拼接或补窗**。

**画对了什么（逐项已验证）**

| 要求 | 实现与验证 |
|---|---|
| 真实轨迹，不是理想方波 | 6 条曲线与 CSV 列**逐元素相同**（637 点，`np.array_equal`）；S0 有 **42 个采样**落在 0.30–0.70 之间；10–90 % 翻转宽度 **2.53 h / 1.23 h**（周期 10.60 h）；相邻采样最大跳变 0.0304 |
| 读窗总宽 20 % | 用 `read_windows.csv` 的原值（未重算）；实测 **2.1333 h = 周期的 0.2013**，且**以谷值为中心**（偏差 0.0000 h） |
| 末端不补齐 | 画面范围与两个读窗全部落在数据内；全 300 h 中**没有任何读窗越界**（0 个） |
| 读带 0–0.30 / 0.70–1.0 | 断言绘制值与轴上的矩形补丁都是 `[0, 0.30]`、`[0.70, 1.0]` |
| 中行量级不同 | `RDF0` 与 `Complex C0` **分轴**，并在图上写明 `RDF0 and Complex C0 on separate axes - not normalised`（RDF0 峰值约 7 a.u.，C0 峰值约 0.085 a.u.，相差约 80 倍） |
| 共用时间轴 | 三个面板 `sharex`，只有底行有刻度与 `Time (h)` |

**布局也是被检查过的**：所有标签两两不重叠、不越出画布、**不压在曲线上**
（逐标签把包围盒换算到该轴数据坐标，检查是否有采样点落入）。
后一项当场抓到两个真实缺陷：`High band` 原本落在 Int 输入脉冲（峰值 6.64）上、
`Low band` 原本压在 S0 曲线上。现在两类标签的位置是**由数据搜索出来**的
（`emptiest_x`：在若干候选 x 中取"S0 距离该电平"与"输入距该高度余量"的较小者最大处），
所以换数据窗口也能自动落到空处。

**图上文字清单（13 个语义文本 + 刻度数字）**：`PB/LR` 无、本图含
`LR fraction S0`（轴标签 + 图例）、`Int input`（图例）、`Int input (C31 flux, a.u./h)`（右轴）、
`Low band`、`High band`、`Read window`、读值 `0` 与 `1`、`RDF0`（图例）、`RDF0 (a.u.)`（左轴）、
`Complex C0`（图例）、`Complex C0 (a.u.)`（右轴）、
`RDF0 and Complex C0 on separate axes - not normalised`、
`recombination flux (a.u./h)`（左轴）、`Forward flux`、`Reverse flux`（图例）、`Time (h)`。

**超出字面标签清单的三处（如需删除请说）**：① 输入轨迹标签 `Int input`（要求只写了"Int输入"）；
② 两个读窗内的读值 `0` / `1`（总览表说 3-2 要讲清"在哪个时间窗读 0 或 1"）；
③ 面板 B 的 "not normalised" 说明（要求允许"不同量级分轴**或**明确标注归一化"，两者都做了）。

## Figure 4-1 — Carry architecture（机制图）

脚本 `make_fig04_1_carry_architecture.py`，**不含任何仿真**：架构全部从模型源码读出来，
并在运行时用冻结产物复核，所以图不可能与模型脱节。

**依据（模型源码，逐条写进脚本 docstring）**

| 源码 | 内容 | 图上的体现 |
|---|---|---|
| `model.py:118` | `g0 = act(A0,K_A[0],n_A[0]) * rep(F0,K_F[0],n_F[0])` | **Carry 0 没有时钟因子** |
| `model.py:119` | `clock = act(y[8], clock_K_au, clock_n)`（y[8] = 成熟游离 Int0） | 时钟来自 **Int0** |
| `model.py:120` | `g1 = act(A1,...) * rep(F1,...) * clock` | **只有 Carry 1 有时钟 AND 门** |
| `model.py:129` | `sources = (None, alpha_Int[0]*g0, alpha_Int[1]*g1)` | 门乘的是**整合酶表达**（输出） |
| `model.py:157` | `auto = rep(A,K_auto1,n_auto1) if j == 1 else 1.0` | **负自调控只存在于 A1**（A0 没有） |
| `model.py:158` | `d[A_j] = alpha_A[j]*PB_j*auto - loss*A` | **PB 驱动 A** |
| `model.py:159` | `d[F_j] = alpha_F[j]*act(A_j,K_A[j],n_A[j]) - loss*F` | **A 驱动 F**；F 的生产指数是 `n_A[j]` |

**运行时的数值复核（不是硬编码）**：`n_A1_gate = 6` 读自
`plausibility/threebit51_selected_v1.json` 的 `frozen_parameters`（SHA256 `6BFD7FF6F0A4D897…`），
`ZENG['n_A'][1] = 4` 由 `import model` 现读（`model.py` SHA256 `256D2D104EECB4CE…`）。
两者都写进自检：若哪天配置变了，图会先失败而不是继续印 6/4。

**你标的"必须画对"逐条验证（45 项自检全过）**

| 要求 | 验证方式 |
|---|---|
| 只有第二级带时钟 AND 门 | `Clock gate` 恰好出现 1 次且在 Carry 1 框内；`Int0 (clock)` 与 `Clock gate` 都不在 Carry 0 框内；`×` AND 节点恰好 1 个且在 Carry 1 内 |
| F 抑制的是**输出** | 从实际画出的 T 形线（`repression()` 返回 None，所以脚本改为截获它添加的 Line2D）读回横杠中心：F0 落在 `Int1 expression` 框边上、F1 落在 `Int2 expression` 框边上，且**都不落在 F 自己的框上** |
| 门乘的是响应，不是浓度 | 图上直接印出模型里真实的乘积：`Carry 0 gate = act(A0) × (1 − act(F0))`、`Carry 1 gate = act(A1) × (1 − act(F1)) × clock`，并写明 "the gate multiplies RESPONSES, not A and F concentrations"；**没有**任何把 A 与 F 浓度相乘的节点 |
| `n_A1_gate=6` 只标在 A1→Int2 门 | 该标签唯一，且 x 坐标在 F1 列右侧（即门臂上）；`n=4` 唯一，且在 A1 下方（即 F1 产生臂上） |
| 方向不能反 | PB→A 箭头断言从左到右且端点贴两框；A→F 箭头断言向下且从 A 框底到 F 框顶 |
| 上方整体链 | 断言六个框的文字按 x 排序恰为 `Clock, Bit 0, Carry 0, Bit 1, Carry 1, Bit 2` |
| 表达链只画一次 | 底部单独一个 `Expression chain` 插图（`mRNA → immature → mature`）+ 一行说明；**没有**在任一箭头旁重复 |

**布局同样被检查**：标签两两不重叠、不压在任何**实心框**上（`scope()` 的虚线框是容器，允许内部放注释——
这条区分是检查器第一版漏掉的）、不压在**任何直线连接件或箭头**上（沿每条线段采样 40 点判定；
唯一的斜线——A1 自调控的弧——不计入并如实报出"1 diagonal connector not sampled"）、以及全部在画布内。

**图上的文字清单（40 个文本元素）**：链 6 项、`Carry 0`/`Carry 1` 框标题 ×2、`PB0/A0/F0/Int1 expression`、
`PB1/A1/F1/Int2 expression`、`Int0 (clock)`、`×`、`mRNA/immature/mature`、
`Delayed repression` ×2、`Negative autoregulation`、`Clock gate`、`Next-bit integrase` ×2、
`n_A1_gate=6`、`n=4`、两条 gate 公式、两段说明、以及
`Modelled interface; molecular implementation pending`。

**两处需要你确认**：
1. `Modelled interface; molecular implementation pending` 我放在**右下角**，语义上指向
   Carry 1 的门（时钟 AND 与 `n_A1_gate=6` 都是模型接口）。若它本意是指别的元件（例如只指时钟门、
   或指 A1 自调控），我把它挪到那一步旁边。
2. 图上额外印了**两条 gate 公式**与"响应而非浓度"一句。这是为了把"必须画对"的第二条钉死；
   若嫌字多，可以删掉公式、只留拓扑（但那样"不是 A×F 浓度"就只能靠箭头语义表达）。

## Figure 4-2 — From modulo 4 to modulo 8（数据图）

脚本 `make_fig04_2_mod4_mod8.py`，**不含任何仿真**。数据：34 状态
`twobit34_results/{trajectories,read_windows}.csv`（300 h，28 窗）与 51 状态
`threebit51_results/wiki_n6_20260924/{trajectories,read_windows}.csv`（600 h，56 窗）。

**图注（可直接用）**：
> **Figure 4-2. 从 mod 4 到 mod 8。** (A) 选定 34 状态模型的 300 h 轨迹（S0、S1），
> 下方为 28 格读数条带，每格一个读窗，环序 1-2-3-0。（B）冻结 n_A1_gate = 6 的 51 状态模型
> 600 h 轨迹（S0、S1、S2）与 56 格条带，环序 1-2-3-4-5-6-7-0。条带数字是**读窗结果**；
> 曲线是状态轨迹，连线不代表读窗之外的每个瞬间都能读成整数。drop 定义为丢弃前 8 个读窗
> （冷启动），实际初态与读窗规则见 Methods。56 个读数含 **7 段完整的 8 值序列**，
> 不宣称观察到 56 次完整状态转移。

**验证**：28 / 56 格与产物窗数一致；每格等于 `read_windows.value` 原值；未定义/失败窗**不被丢弃**
（本图两份产物恰好没有未定义窗，检查会如实报出数量）；曲线与 CSV 列逐元素相同；
`28 reads` / `56 reads` 是图上仅有的两个计数；无 `100%`／`56 transitions` 之类过度声称。

> ⚠️ **Figure 4-2 里没有上游。** Panel A 只画 `S0`/`S1` 与读窗条带 —— 脚本里
> `flux|C31|upstream|TetR|CI|LacI|oscillator` **零命中**。所以它回答的是"读数对不对"，
> 不回答"这个计数是不是真的架在韩的振荡器上跑出来的"。后者由下面的 **4-2a** 承担。

## Figure 4-2a — Two-bit counting under the real Han oscillator（数据图，2026-09-27 新增）

脚本 `make_fig04_2a_twobit_under_upstream.py`；取数脚本
`make_twobit_under_upstream_run.py`（**一次 300 h 仿真**，约 5 min）。数据：
`data/twobit_selected_300h.{csv,json}`。

**为什么新增**：交付图集里**没有任何一张**把"上游振荡器在振荡"和"两位计数器在完整 mod 4 计数"
放在一起。`twobit34_results/figures/01_counter_overview` 只有一条 C31 通量线，
`03_state_groups` 的振荡器是归一化定性图，4-2 干脆没有上游。

**三个面板**（共用 300 h 时间轴；精简版要求"单张最多 3 个面板"）：

| 面板 | 内容 |
|---|---|
| A | **韩氏振荡器本体**：`TetR total` / `CI` / `LacI`（各自归一化到自身最大值）＋ 右轴 **C31 翻译通量 (a.u./h)** —— 驱动 bit 0 的那一路 |
| B | **两位 DNA 状态** `b0_S` / `b1_S`，0–0.30 低带与 0.70–1 高带底纹，28 个读窗阴影 |
| C | **解码读数**：0→1→2→3→0 阶梯，逐窗标值，附 28 位码串 |

**图注（可直接用）**：
> **Figure 4-2a. 两位计数器在真实振荡器输出下完成 mod 4 计数。**
> (A) 韩氏振荡器的三个蛋白（各自归一化到自身最大）与其 C31 翻译通量（右轴，a.u./h）；
> 这就是驱动 bit 0 的真实上游。(B) 两位的 DNA 构象状态 `S0`、`S1`，阴影为 0–0.30 低带与
> 0.70–1 高带，浅蓝为 28 个读窗。(C) 解码读数成阶梯上升，**7 个完整 mod 4 循环**，
> 码串 `1230123012301230123012301230`，无未标注窗。运行 **认证通过**，
> bit 1 setup **4.64 h** / hold **2.91 h**，最小承诺度 **1.00**。
> 工作点：`uM_per_au = 5.75`、进位 A/F 成熟 **32.5 min**、进位 mRNA **2 min**、
> clock **0.3 / 2.0**、额外稀释项 = 0（`selected_profile.json`）。
> **耦合是单向的，且运行在上游无负载极限（rho = 1）**；`uM_per_au` 未经实验标定。

**工作点声明（重要）**：本图跑在**选定档 `5.75 / 32.5 / 2`**，即 51 状态链、140 点扰动扫描、
`initial_states_selected` 共用的那一个；**不是** `run_twobit34.py` 裸 `TwoBit34Model()` 落到的
`nominal_extension()` 默认档 `6.0 / 30.0`。两者实测差别：

| | 本图（选定档 5.75/32.5） | 主交付 `twobit34_results/`（6.0/30.0） |
|---|---|---|
| bit 1 setup | **4.637884 h** | 4.430958 h |
| bit 1 hold | **2.913723 h** | 2.799375 h |
| reads / 码串 / certified | 28 / `1230`×7 / True | 28 / `1230`×7 / True |

取数脚本因此**逐字段断言**模型等于 `selected_profile.json`，并**断言自己没有落在裸默认值上**
（`uM_per_au ≠ 6.0`、成熟 ≠ 30.0）—— 这正是当初缺的那道闸。它还把结果与 75 点网格的同一行
逐位交叉核对（`grid_row_cross_check.matches = True`）。

**验证（42 项自检，`FAILURES: 0`）**：7 个必需标签逐字存在；三条振荡器曲线各自等于归一化后的
原始列、通量曲线等于原始 `C31_flux_au_h` 列、S0/S1 曲线等于原始 `b0_S`/`b1_S` 列（逐元素）；
`b0_S`/`b1_S` 落在 [0,1] 且都访问过低带与高带；阶梯等于全部 28 个 `read_windows.value`（未改一个）；
码串恰为 7 个完整 mod 4 循环且相邻值 `+1 mod 4`；工作点与 provenance 断言（含网格行交叉核对）；
**页脚三要素 `Model:` / `Input:` / `Parameters:` 齐备**；**图内不含任何泄漏数字**
（两位认证没有泄漏臂，所以没有可引用的数）；无标签重叠、全部在画布内；SVG 保留真文本。

## Figure 4-3 — A sharper gate enables three-bit counting（数据图）

脚本 `make_fig04_3_sharper_gate.py`。Panel A 的两条 S2 轨迹来自**同一表达设置**：
n=6 用既有冻结产物 `wiki_n6_20260924/trajectories.csv`；n=4 用
`data/n4_matched_600h.{csv,json}` —— 由 `make_n4_matched_trajectory.py` 生成的
**新匹配运行**（同一 extension/carry，只有门指数取 4；`build_threebit` 会断言实测指数）。
Panel B 用 `plausibility/carry_pairing_grid_verdict.json` 的逐点与逐 n 中位数。
两条曲线都读**名为 `b2_S` 的列**（模型自己对 bit2 DNA 构象态的命名），不读裸 `S2`。

> ⚠️ **已修正的事故（2026-09-26）**：`make_n4_matched_trajectory.py` 原先把列写成
> `S0 = y[44], S1 = y[45], S2 = y[46]`，即把 **b2_S、A1、F1** 挂在了 S0、S1、S2 名下
> （真实索引是 S0 = 16、S1 = 27、S2 = 44）。于是 Panel A 的 n=4 面板画的是 bit2 的
> **F1 蛋白**，而 n=6 面板画的是 bit2 真正的 S2 —— **同一个轴标签下是两个不同物理量**；
> 原图注「S2 峰值 2.25、超出物理区间被裁掉」描述的其实是 F1（以 a.u. 计的浓度，本来就不受 1 约束）。
>
> **原错误数据已归档**（H3/H7，不是就地覆盖）：`archive/` 下有
> `n4_matched_600h_gen1_mislabelled_S0S1S2.csv`（sha256 `FBE88CF7…`，就是产出错误图的那一份）
> 与 `..._gen2_relabelled_9col.csv`（sha256 `F2DB0AA4…`，改对列名但尚未补齐 bit2 池与通量的那一版），
> 逐件哈希记在 `archive/ARCHIVED_HASHES.json`（由 `archive_n4_generations.py` 生成，
> 归档前后哈希必须相等，否则脚本报错退出）。
>
> **新旧列对应验证（"修复只改标签、数值未变"由此确立）**：逐列逐元素比对 ——
> `gen2[b2_S] ≡ gen1[S0]`、`gen2[A1] ≡ gen1[S1]`、`gen2[F1] ≡ gen1[S2]`，**`max|Δ| = 0`**；
> `time_h / g1 / clock_gate / C31_flux` 四列亦逐位不变。随后补齐 `b2_I / b2_R / b2_C /
> J_fwd2 / J_rev2` 五列时，原 9 列再次逐位不变。
>
> **判据数字不受影响**：`crossings = 5`、`gates = 14`、17 个未确定读窗、`certified = False`
> 都来自 `verify_threebit51.analyse_threebit`（它读 `states = (y[16], y[27], y[44])`，一直是对的）。

**图注（可直接用）**：
> **Figure 4-3. 更陡的门使三位计数成立。** (A) 同一表达设置下 n=4 与 n=6 的 S2 轨迹
> （[30, 150] h，读窗阴影；**这一段是启动阶段示例**，全程证据见下）。n=4 **未认证**：
> **首次进入高态后反向写入不充分**；本次 600 h 运行中，**首次进入高态（t ≈ 33.1 h）后未重新
> 进入低态带**（全程最小值 0.31），**无法连续完成模 8 计数** —— 14 次进位里只有 5 次翻转 bit 2，
> 17 个读窗的位标签从未确定；本窗口内读数为 9 个 1、**0 个 0**、2 个未确定。
> n=6 **已认证**：14/14 进位各翻转一次，同一窗口内 7 个 1、4 个 0、0 个未确定，是干净的
> mod-8 交替。
> (B) 五个门指数下的配对反向通量比 `L_symmetric`：每个指数 **4 个实际配置点**（空心圆）与其中位数
> （横杠），星号标出选定的 n=6。**5% tail cut; 2 min output grid; 4 model configurations per exponent。**
> 这 4 个点是 **4 个参数组合，不是实验重复**；n=4 的失败是真实动力学失败，不是评分规则造成的；
> n=5 只是**已测试整数中首个通过的值**，不是已定位的精确连续分界。

**「反向写入不充分」这句话的证据（全程 600 h，不是那 120 h 窗口）**：`data/n4_matched_600h.json`
的 `excursion` 块给出 —— 首次进入高带 `t_first_high_h = 33.13`；此后 `b2_S` **从不**回到
`≤ 0.30`（`returns_to_low_band_after_first_high = false`，17007 个采样点最小值 **0.3087**）；
反向通量积分 `∫J_rev2 dt`、正向 `∫J_fwd2 dt`、以及两者的全程值、峰值、5 次穿越的时刻表与
14 次进位门的时段表都记在同一块里。**生成脚本会在该断言不成立时直接报错退出**，
所以图注这句话不可能与数据脱钩。

**验证**：两条曲线与各自的 `b2_S` 列逐元素相同、共用同一时间轴 [30, 150] h、且都保留低态与真实过渡；
**两条曲线都落在 [0,1] 内**（错列如 F1 会立刻越过 1，这一行就是当初能立刻抓住该 bug 的闸）；
n=6 产物的派生 `S2` 列与 `b2_S` 逐位相同；n=4 文件通过 `check_state_column_semantics.py`
（按模型 `state_names` 定位索引、重放 56 个读窗的 bit2 标签得 0 不符、并用 `g1` 恒等式做**判别性**检验
—— 把 A1/F1 对调会差 0.389）；`excursion` 的两个断言被单独检查（含 `∫J_rev2` 是否记录在案）；
n=4 的运行在其 provenance 里自述为"新运行、非既有产物"且实测指数 = 4；四个中位数
0.4489/0.0521/0.0169/0.0139/0.0131 与冻结值一致；20 个点分 5 组每组 4 个；星号在 n=6；
两条必须保留的措辞（"4 个参数组合不是重复"、"n=5 是首个通过的被测整数"）在图上；
图上另有一行页脚注明「[30, 150] h 是启动阶段示例」。当前脚本输出 **67 项通过 / `FAILURES: 0`**。

## Figure 4-4 — Validation summary（验证摘要图）

脚本 `make_fig04_4_validation_summary.py`。Panel A 用 `plausibility/eight_initial_all.csv`
的 8 条 `cold_sequence`（每条 56 读）；Panel B 用 `plausibility/pool_perturbations_all_all.csv`
（136 次）与其在 `threebit51_selected_v1.json` 里的 `pool_perturbations.overall`
（128 恢复 / 8 合法相移 / 0 失锁），加 `postfreeze_peak_perturbation/postfreeze_all.csv`
（30 次扰动 + 6 条无扰动对照）。

**图注（可直接用）**：
> **Figure 4-4. 验证摘要。** (A) 八个轨道一致初态（行标 000–111）各 56 个读窗的条带，
> 颜色为解码值 0–7（右侧数值色标）；可见不同初态之间只是**相位偏移**，8/8 序列正确。
> (B) 两组扰动运行：冻结前 136 次 = 128 恢复原相位 + 8 次合法 mod-8 相移 + 0 失锁；
> 冻结后 30 次 = 30 恢复 + 0 + 0。8 次合法相移来自 S2 翻转；136 次中部分近零池扰动本身很弱
> （例如 C2×0.8：最大变化 3e-12、2 条退化运行）；冻结后另有 6 条无扰动对照，且部分扰动提前
> 1.6–40.8 h 施加。严格容差 5 点 × 3 档结论一致（不另设热图）。**仅限所列扰动集合**，
> 不构成"任意噪声下鲁棒"或实验成功率的声明。

**验证**：8×56 格全部等于产物序列；8 个初态全部认证；色标存在且为 0–7 数值刻度；
128/8/0 与 136、30/0/0 与 30 都是**从产物重算**而非手写；禁用措辞
（`100%`/`arbitrary noise`/`any noise`/`experimental success`/`guaranteed`）在图与 SVG 中都不存在；
`Specified perturbations only` 作为限定语在图上。

## 排版对应（按"只保留三件事"）

三段结构对每张图都成立：**一句话结论 → 主图 → 简短解释**。
参数表、ODE、旧扫描与复现信息放 Methods／补充材料；**shutdown 接口只留文字**，
本轮不新增第 7 张主图。上面每张图的"图注（可直接用）"就是那段简短解释的候选文案。





## 两个文件的分工

- **`tempo_style.py`**：调色板、线型约定与绘图原语。所有颜色与线型与 Oscillator 章节一致，
  这样同一个 Wiki 上的机制图属于同一视觉家族。
  - 原语：`box` / `arrow` / `repression`（T 形端点）/ `binding`（双向箭头）/
    `and_node`（× 号 AND 节点）/ `hline`（阈值虚线）/ `scope`（模块边界虚线框）/
    `hollow_arrow`（"尚未建模"的空心虚线箭头）/ `legend_lines`（线型图例）
- **`make_schematics.py`**：每张图一个函数，坐标全部用数据坐标显式写出，便于微调。

## 约定

1. **全部标签只用英文**。作图说明要求"英文栏可直接复制到图片中"，
   而且英文彻底避开中文字体依赖，headless 渲染不会出现缺字方框。
2. **不要改颜色**。若要调整，改 `tempo_style.py` 里的常量，所有图一起变。
3. **AND 必须画成 `and_node`（× 号圆节点）**，不要用加号或并联箭头——
   "相乘"是进位机制成立的前提。
4. **抑制必须用 `repression`（T 形端点）**，与激活箭头区分开。
5. 需要新增图时，在 `make_schematics.py` 里写一个返回 `save(fig, ...)` 的函数，
   并注册到末尾的 `FIGURES` 字典即可。

## 尚未由代码生成的图

| 图 | 原因 |
|---|---|
| 3-2 单周期翻转分解 | 需要冻结模型的真实轨迹 |
| 3-3 读窗几何与驻留带 | Panel 1–2 是定义、可以画；Panel 3 需要真实轨迹做新旧判据对照 |
| 3-4 工作窗口 | 需要 `uM × 进位成熟时间` 二维细化扫描结果 |
| 3-5 扰动与容差 | 需要四组分子池扰动扫描结果 |
| 4-2 门对比度 | 需要单周期的 $g_0$ 曲线 |
| 4-3 进位因果链 | 需要 $J_{rev,0}$ / $g_0$ / 重组酶 / $S_1$ 四行同步轨迹 |
| 4-4 300 h 两位读数 | 需要 300 h 轨迹 |
| 4-5 工作窗口热图 | 同上，需要二维扫描结果 |
| 4-7 三级机制与悬崖 | 需要 63 点扫描的逐点结果 |

这些可以由同一套 `tempo_style.py` 原语加上仿真输出来生成，配色自动保持一致。
