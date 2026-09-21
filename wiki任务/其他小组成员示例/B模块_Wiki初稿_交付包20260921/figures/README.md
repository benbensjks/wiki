# figures 目录说明

本目录含本轮概念草图、参数表图与候选数据图。候选图均来自 B 模块既往交付，**未修改数据**；
文件名前缀 `candidate_` 表示可直接改造为 Wiki 正式图（重新配色、统一字号、拆分面板）。

## 本轮概念图 / 表图

| 文件 | 内容 | 生成方式 |
|---|---|---|
| `fig03_1_single_bit_mechanism_sketch.png` | Figure 3-1 状态依赖的开关（PB/LR 双卡 + 生产分配 + RDF 池量 + 两条脉冲弧线） | `../scripts/make_mechanism_v2.py`（matplotlib） |
| `fig03_2_design_window_sketch.png` | Figure 3-2 设计窗口与失效模式（tag 窗口三段 + 联合表达窗口示意） | `../scripts/make_design_window_v1.py`（matplotlib） |
| `fig03_7_parameter_table.png` | Figure 3-7 参数表图（模型条件 / 主要方案 / 设计目标 / 工作窗口，含证据标签） | `../scripts/make_parameter_table.py`（matplotlib） |
| `optional_concept_cell_context.png` | 可选生成式概念图（细胞情境，供 hero / 情境插图；生成式草图，需美工统一风格） | 千问生图（提示词存档：`../scripts/prompts/`） |

## 候选数据图

| 文件 | 原始文件 | 来源交付包 | 建议用途 |
|---|---|---|---|
| `candidate_clock_input_flux.png` | `fig1_input.png` | Week 4 交付包（2026-09-06） | Figure 3-3 上半：A 输入通量与 C31 蛋白 |
| `candidate_steady_state_trajectory.png` | `primary_mechanism.png` | 非零泄漏交付包（2026-09-19） | Figure 3-3 下半：Int 脉冲 / $S$ / RDF / BM3R1 稳态轨迹 |
| `candidate_rbs_tag_workspace.png` | `fig7_rbs_tag_region.pdf`（已转 PNG） | Week 4 交付包 | Figure 3-4(A)：RBS × tag 联合工作区 |
| `candidate_feasible_region.png` | `feasible_region.png` | 非零泄漏交付包 | Figure 3-4(B)(C)：$(\beta_B,\beta_R)$ 可行域与泄漏容差 |
| `candidate_tag_mechanism.png` | `fig5_tag_mechanism.png` | Week 4 交付包 | Figure 3-5(A)：间隙 / 峰值 Int vs $k_{tag}$ |
| `candidate_heterogeneity_mc.png` | `fig8_stochastic_fidelity.pdf`（已转 PNG） | Week 4 交付包 | Figure 3-5(B)：异质性蒙特卡洛保真度 |
| `candidate_robustness_and_initialization.png` | `robustness_and_initialization.png` | 非零泄漏交付包 | Figure 3-6(A)：启动相位与扰动通过数 |
| `candidate_long_run_comparison.png` | `long_run_comparison.png` | 非零泄漏交付包 | Figure 3-6(B)：旧点失效 / 两套方案恢复的长程对比 |
| `candidate_input_clearance_pairing.png` | `input_clearance_family.png` | 非零泄漏交付包 | Figure 3-6(C)：$k_I$ 与输入 scale 的配对规则 |
| `candidate_week4_operating_point.png` | `fig4_operating_point.png` | Week 4 交付包 | 备用：Week 4 工作点全程（可替换 3-3 下半） |

## 备注

- Figure 3-1 中池量（≈0.002 / ≈4.5 µM）的复算脚本：`../scripts/verify_pool_levels.py`；
- one-pulse-one-flip 时序图作为配套概念图**待日后细化**（本轮不作为美工交接图）；
- 生成式概念图的提示词存档：`../scripts/prompts/`（可复现）；
- 图内坐标、单位与图注口径以原交付包报告为准；
- 若第 3 章最终采用 D 的统一 core 模型口径，本目录候选图需相应替换。
