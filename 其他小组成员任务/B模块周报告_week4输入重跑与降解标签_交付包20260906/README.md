# Week4_RDF_REV — B 模块修订交付

本目录为 iGEM 2026 TEMPO B 模块周报的修订版本（审查后重建），与原始交付 `../Week4_RDF` 平级。所有数值均可由下述脚本复现。

## 文档

| 文件 | 说明 |
|---|---|
| `报告.tex` / `报告.pdf` | 主报告：week4 定稿输入下的 RDF 计数重验证 + Int 降解标签设计（14 页）。编译：`xelatex 报告.tex` 两遍 |
| `Int降解标签选型报告.tex` / `.pdf` | 标签选型专项报告（6 页）。编译同上 |
| `调研_降解标签/` | 7 篇一手材料（5 篇全文 + Jadhav SI 全表 + Lies 摘要级核验存档）与两份调研 md |
| `docs/week4_unloaded/` | A 模块 week4 冻结参考输入（2026-08-14，unloaded，100 h） |
| `figures/` | fig1–fig8（矢量 PDF） |

## 主要修订点（相对审查前版本）

数值口径与科学结论层面的改动，均来自复核与补充实验：

1. **输入剂量**：§1 表“每周期通量剂量”按通量-时间积分实测为 22.8 µM（峰间谷逐周期积分，σ<0.001）；另列“蛋白×周期积分 ≈27.5 µM·h”作为对照口径。原 28.5 µM 为蛋白×周期积分值，系口径混淆。
2. **判据统一**：主判据 P = banded toggle score（与 Week 3 一致）；辅判据 M = 二值化匹配率，用于边界交叉对照，凡引用 M 处均注明。联合工作区以 P 判据重扫（`scan_joint_v2.py`，tag 1–28 × RBS 0.30–0.60，392×2 模拟）。
3. **最优平台修正**：原“tag 4–18 全域可行”不成立（RBS0.6×tag4/5/6 为 0.0/0.6/0.2）。联合区（P 判据）实测全稳平台 = **tag 7–18 × RBS 0.30–0.60**（≥0.99；tag 6 需 RBS≤0.55）；低速档 k≤4 需 RBS≤0.45–0.55；高速档 19–28 需逐步上调 RBS（0.35–0.55）。
4. **上下缘精扫**（`edge_scan.py`）：下缘双初态同时通过 = 2.8（2.4–2.7 初态依赖的易碎带）；上缘 RBS0.45 = 24 双初态 0.997，25/27 碎裂（0.24）、26 边缘（含中间态样本）、28+ 卡滞失效。原“26 起失效、上界≈25”的表述过粗，且未覆盖 30–36 的“部分复原”假象（实为回翻失败的碎裂退化带，非单调）。
5. **上界机制修正**：`exp1_upper_diag.py` 显示高 tag 失效形态为 LR→PB 回翻失败/周期卡滞，而非单纯“峰被削平”。
6. **蒙特卡洛重做**（`rand_robust_v2.py`）：修正相位采样错位（按偏移谷时刻采样）；12 档 ×30 实例（RBS0.45 基点 ±0.04 扰动）。随机全稳 = **6–18**（0/30 失败），5 与 22 边缘，3 与 25–26 退化。原“窗口 5–22”上端过宽；且原 MC 的 M 判据对相位错位不敏感,对易碎档给出过高得分。
7. **LAA-LAA 定量**：SI 高诱导残留为原生的 13.7% vs 37.2%（比例 0.37，即 ≈1/2.7；速率比 3.74×），非“约为 1/4”。换算公式 `k=k_dil(1/R−1)` 显式写出，全诱导点换算表见 `si_k_convert.py`（显示降解机强饱和）。低诱导速率比（4.8×）锚定文献原生 3–8 → LAA-LAA ≈14–38。
8. **文献核验与修正**：Lies 2008 全文非 OA（PMC2516991），所引数字（≥90%、≥1.4 min⁻¹、clpX⁻ 0.10 min⁻¹、SspB⁻ 3×）经 PubMed 摘要逐项核验并归档摘要；其“clpX⁻≈5 min”改为 ≈7 min。Andersen E. coli 档 LAA/LVA 40、AAV 60、ASV 110 经全文复核；P. putida 中 LVA≈60 min、其余≈190 min（表中区分）。Jadhav “ClpP 过表达无效”表述修正为“ClpP 中等效果（机制为增加 X/P 相遇）、SspB 无效”。Bernabé-Orts 文献改动为直接引用其“degradation tags to regulate the half-life of the integrase”原文语句（原“26S 体系”无来源支持）。
9. **“11 周期”修正**：稳态谷采样为 5 点（100 h ≈9.4 驱动周期），原数字来自 Week 3 v36 语境（11 次脉冲），不再沿用。
10. **复现指引补全**：新增 `edge_scan`、`scan_joint_v2`、`rand_robust_v2`、`si_k_convert`、`exp1_upper_diag`、`make_fig6/7/8_v2`，并声明重扫描策略（本地 12 核 / srv2026 32 进程；两机结果 392 项中仅 2 项 3 位小数差 0.001）。

## 环境与复现

- conda env `igem-w3`（本机与 srv2026 各一份；numpy 2.5.1/2.5.2、scipy 1.18.0、matplotlib 3.11.1）
- 输入 CSV 在本目录 `docs/week4_unloaded/`（脚本 `couple_oscillator.py` 相对定位）
- **重扫描**（`scan_joint_v2`、`rand_robust_v2`、`edge_scan`）：`NPROC=32 python xxx.py`，srv2026 全量 ~1 min；本地 12 核即可复现，结果一致
- 判据、初值、积分器设定见 `报告.tex §2`
- 绘图：`python make_figs.py`（fig1–5）、`make_fig6_v2.py`（fig6）、`make_fig7_v2.py`（fig7）、`make_fig8_v2.py`（fig8）

## 文件约定

- `model/` 中脚本分三组：最后管线（上文引用各脚本）、复核/探针（`verify_*`、`diag_default`、`exp1_upper_diag`、`mech_tag`、`smoke0`）、遗留（Week 3 v36 时代及过程脚本：`scan_*`（旧）、`diagnose*`、`fig_baseline_req.py`、`verify_audit.py`、`final_run.py`、`leak_rate.py`、`run_bm3r1.py`、`zhao_core.v0730.py.bak` 等）。遗留脚本未被引用报告，导入路径为旧的绝对路径（Week4_RDF/model），不作为复现组成部分。
- `*.json` 为各扫描输出，含 `cons`（DNA 守恒）字段。
