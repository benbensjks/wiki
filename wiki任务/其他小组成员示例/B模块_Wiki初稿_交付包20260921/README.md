# B 模块 9.20 交付物：第 3 章 Wiki 初稿、机制图需求与 Engineering 构想

**2026-09-21（第二版）｜TEMPO 建模组 B 模块（王厚骅）｜对应组长 2026-09-17 22:18 通知**

组长通知要求每个文件夹包含三件东西：① 自己那一块 wiki 的初稿 md（大致小标题即可）；
② 机制图的绘图要求说明（最好附简单草图）；③ Engineering 构想。本文件夹按此组织。

**语言约定：Wiki 正文与结构以英文为准（`01/02/03` 三个英文主稿）；中文稿仅作辅助阅读，
放在 `zh_aux/`。** 交付范围只覆盖 **B 模块（第 3 章单比特计数器）**；
对其他章节只保留边界指向（如"多比特进位见第 4 章"），不写其内容。

| 文件 | 对应要求 | 内容 |
|---|---|---|
| `01_Chapter3_Wiki_Draft.md`（EN 主稿） | ① wiki 初稿 | Chapter 3 "Coupling the Clock to a Reliable Single-bit Counter"：page flow、modeling status、mechanism、criteria & operating window、final-clock re-validation、BM3R1–RDF delay circuit、Int degradation tag、promoter、robustness、engineering summary、scope、conclusions、figure list |
| `02_Figure_Brief.md`（EN 主稿） | ② 配图 | Figure 3-1（状态依赖的开关）与 Figure 3-2（设计窗口与失效模式）的构图、禁止画错点、英文标签与图注草稿；另含参数表图与可选生成式概念图说明 |
| `03_Engineering.md`（EN 主稿） | ③ Engineering | three Design → Build → Test → Learn → Redesign rounds + **cross-module A↔B feedback/rework iteration**、design-rule table、Learn list、experiment interfaces、**judging-criteria alignment (2026 Judge Handbook / Best Model)**、**external Engineering benchmarks (IZJU-China 2025 / Heidelberg 2025 / USTC 2025)** |
| `zh_aux/01–03_中文辅助.md` | — | 上述三稿的中文辅助版（便于本人与网页组阅读理解；内容与英文主稿一致） |
| `figures/` | — | 概念草图（fig03_1 机制、fig03_2 设计窗口）、参数表图（fig03_7）、可选生成式概念图（细胞情境）+ 10 张候选数据图；来源对照见 `figures/README.md` |
| `scripts/` | — | 各图生成脚本（`make_mechanism_v2.py`、`make_design_window_v1.py`、`make_parameter_table.py`）与池量复算脚本（`verify_pool_levels.py`），可在 srv2026 `igem-tempo-2026` 环境重跑 |
| `B模块_Wiki初稿_交付包20260921.zip` | — | 本文件夹完整打包（发群用，UTF-8 文件名，附 SHA-256） |

---

## 1. 本稿基于哪些已完成工作

本稿不是新计算，而是把 B 模块已交付并核验过的工作整理成 wiki 叙事：

| 来源 | 日期 | 本稿引用内容 |
|---|---|---|
| 《Week 3 · B 模块：ϕC31 二进制计数器单级开关》 | 2026-07-30（v3 08-01） | BM3R1 替换、失效机制、接口泄漏需求（≤1.5%）、初版工作点 |
| 《week4 定稿输入下的 RDF 计数重验证与 Int 降解标签设计》 | 2026-09-06 | 通量接口语义、新工作点、K 敏感带消失、标签窗口与选型 |
| 《非零 RDF 启动子泄漏下的工作方案》 | 2026-09-19 | 非零泄漏方案（β_B=2/β_R=8/scale=0.30/k_I=12）、可行域、启动条件、η 不确定性 |
| 《启动子审查与实验交接》 | 2026-09-19 | Cello pBM3R1 66 bp 推荐、36 bp 序列核验 |

数字均可在上述交付包中找到；本稿只做转写、组织与图表规划，未新增仿真结论。
复现环境：srv2026 conda `igem-tempo-2026`（Python 3.12）；Week 4 扫描脚本在
`~/work/Week4_RDF_REV/model/`，非零泄漏交付包在 `~/work/B_nonzero_leak_20260919/`。

---

## 2. 与 D 已交付第 3 章草稿的关系（需协调，请组长/组员确认）

胡犇岩同学 9/20 交付的文件夹中已包含一份第 3 章草稿
（`02_第3章_单比特计数器.md`，基于其仓库中的约化 single-bit 模型：
11 个连续状态（mRNA–未成熟–成熟链 + 复合物 + DNA 比例）、浓度标尺换算与 M1–M6 读窗判据），
并生成了单比特机制图 `fig03_single_bit_mechanism.png`。两份草稿描述同一模块，但模型口径不同：

| | 本稿（B 模块） | D 的草稿 |
|---|---|---|
| 模型 | 38 维 Zhao 2019 核心 + BM3R1 + 显式泄漏 | 11 状态约化核心（含成熟链） |
| 用途 | 单比特模块自身的机制与实验设计参数 | 其进位 / 级联建模的比较基线（据 ZUHUI 讨论） |
| 判据 | 端点 $Q$ + 中点穿越 + 计数奇偶性 | M1–M6 读窗判据 |
| 数据 | Week 3 / Week 4 / 9-19 交付包 | D 仓库冻结模型 |

**建议在 9.30 冻结前明确其中一种处理**：

1. **统一口径**：选定一个模型作为第 3 章权威口径（需要两人对照后决定），
   另一份作为补充/交叉验证；
2. **分工合并**：机制与判据（D 草稿） + 设计参数与实验接口（本稿）合成一章，明确每节负责人；
3. 无论哪种处理，**图号需要统一**（两份草稿都使用 Figure 3-x，存在重号）。

本稿以"B 模块能直接交付的实验设计信息"为重：启动子、表达配比、标签窗口、
启动条件、η 判别——这些内容不依赖上述模型口径选择。

---

## 3. 证据边界（务必随文保留）

- 全部结论为**确定性模型预测**，不是 MC4100 湿实验成功率；
- β 是模型翻译系数，不是启动子 RPU，也不唯一对应 RBS 序列；
- η=1 是工作假设（Int 清除时 RDF 完整返回），需实验判别；
- 判据（端点 $Q$、中点穿越）是工程验收判据，不是实验标准；
- 输出时长 / 自杀效应器不在本章范围，BM3R1 池是登记的共享接口。

---

## 4. 待确认 / 下一步

1. 第 3 章模型口径与图号（见第 2 节）；
2. 最终 caption 与全文译校（10.5 前随完整 Markdown 交付）；
3. 9.30 冻结前补：η 判别实验设计、表达标定换算表更新；
4. 网页组件与交互（见 `03_Engineering.md` §9）；
5. 图集状态：Figure 3-1/3-2/3-7 草图已交；one-pulse-one-flip 时序图待细化
   （本轮不作为美工交接图，后续单独提供）；生成式细胞情境图为可选资产。

---

## 5. 交付前自检记录

- 数字核对：$Q=0.9966$、$\beta$ 参数、泄漏容差（≤0.01）、标签窗口（3–24 / 6–18）、
  启动通过率（37/40、200/200、199/200）等均与 9-19 交付包
  `data/final_summary.json`、`design_targets.json` 一致；
- 第 3 章正文与 Figure 3-1 的"脉冲到达前 RDF 池相差约 2500 倍"由非零泄漏交付包
  `point_B2.0_R8.0_s0.3_I12.0_f1.0_PB.csv` 逐脉冲复算得到
  （LR 态 4.457 µM vs PB 态 0.00175 µM，比值 ≈2550×；复算脚本 `scripts/verify_pool_levels.py`）；
- 英文主稿与中文辅助稿内容一致；md 结构为 `#/##/###` 三级，无跳级；
- 各图脚本（`make_mechanism_v2.py`、`make_design_window_v1.py`、`make_parameter_table.py`）可在 srv2026 `igem-tempo-2026` 环境重跑；
- 外部 Wiki 对标（`03_Engineering.md` §7）统一使用**浏览器渲染**核对
  （headless Chrome）：iGEM Wiki 存在客户端渲染页面，直接 HTTP 抓取会误判为空页；
  USTC 2025 的 Engineering/Model 页经浏览器渲染后确认内容完整（Iteration 1–4 + 五个模型），
  本版已按实际内容改写分析；
- 候选数据图均为 B 模块既往交付中的图，未修改数据。
