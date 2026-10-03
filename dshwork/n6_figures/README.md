# n = 6 参数分析与鲁棒性图（51 状态三级冻结点）

> **范围声明**：只画 **n = 6 冻结工作点**的参数分析与鲁棒性证据。
> **未改动任何现有文件**——本目录为新建目录，脚本只读现有产物、只向 `./out` 写入。

---

## 1. 数据来源（全部只读）

| 用途 | 文件 |
|---|---|
| 门指数选型 | `final_reconstruction/plausibility/carry_pairing_grid_verdict.json`（20 点 × 3 档尾巴切分） |
| 八初态（序列/裕量） | `final_reconstruction/plausibility/eight_initial_verdict.json`（`per_value`） |
| 八初态（**配对**泄漏） | `final_reconstruction/plausibility/carry_pairing_eight_verdict.json`（`per_point`） |
| 池扰动 136 条 | `final_reconstruction/plausibility/pool_perturbations_all_all.csv` |
| 严格容差 15 点 | `final_reconstruction/plausibility/strict_tolerance_tol_all.csv` |
| 冻结参数与汇总 | `final_reconstruction/plausibility/threebit51_selected_v1.json` |

**冻结配置**（四张图全部基于此）：

```text
n_A1_gate            = 6
carry1 mRNA 半衰期    = 2 min
A1 / F1 成熟半衰期     = 32.5 min
F1 产生 Hill 指数      = 4      （冻结表值，未被 n_A1_gate 改动）
uM_per_au            = 5.75
model_threebit51.py  sha256 = D244D3CD1914F50F994E6856084D4A012B820C13249DCBA0E963E444CE3C34CC
```

---

## 2. 四张图

### `fig_n6_01_gate_exponent_selection.png` — 门指数选型（**为什么选 n=6**）

| 面板 | 内容 |
|---|---|
| **(A)** | `L_symmetric` vs `n_A1_gate`（对数纵轴）。20 个网格点全部画出，三档尾巴切分（1 % / 5 % / 10 %）各一条中位数线 |
| **(B)** | 20 点认证网格：`n_A1_gate` × (mRNA, 成熟) 四组合 |
| **(C)** | 三档切分下的排序，验证 `8 < 7 < 6 < 5 ≪ 4` 在三档下完全一致 |

**核对结果**：`n=4 → 0/4`；`n=5,6,7,8 → 各 4/4（共 16/16）`；`ordering_stable_across_cuts = True`；`n6_beats_n5 = True`。
**关键读法**：n=4 与 n≥5 之间是**离散通过边界**（跨了一个数量级），n=6–8 是**平台**（差异 < 2 倍）。所以 n=6 不是"最优"，而是**平台上的第一个点**——取它是因为它离边界最远、且已经在平台上。

> **Caption draft.** Gate-exponent selection at the frozen working point. (A) The paired leak
> metric $L_{symmetric}$ against the gate exponent for all 20 grid points, with the median at
> each of the three tail cuts. (B) Certification over the 20-point grid: $n=4$ gives 0/4 and
> $n\ge5$ gives 16/16. (C) The ordering $8<7<6<5\ll4$ holds at every tail cut, so $n=6$ is
> chosen as the first point of the $n=6$–$8$ plateau rather than as the minimum.

### `fig_n6_02_pool_perturbation_matrix.png` — 池扰动矩阵（136 条）

| 面板 | 内容 |
|---|---|
| **(A)** | 6 个扰动对象 × 7 种扰动（0.8×/1.2×、±0.1/±0.2、全翻转）的矩阵，格内为"恢复原相位的条数 / 总数" |
| **(B)** | $S_2$ 全翻转的相移：每个初态扰动前 → 扰动后的首位数字 |
| **(C)** | 结局统计 |

**核对结果**：136 条 = **128 恢复原相位 + 8 合法模 8 相移**，**0 失锁**、**0 读出通过但因果失败**、`same_type_adjacent` 全为 0、`certified` 全为 True。
**$S_2$ 全翻转**：8 个初态**逐一无误**地平移 **+4 mod 8**（1→5、2→6、3→7、4→0、5→1、6→2、7→3、0→4），按合法相移吸收，不计失败。

> **Caption draft.** Molecular-pool perturbations at the frozen point. (A) 136 runs over six
> pools and seven perturbation types; every cell keeps counting. (B) A full $S_2$ flip shifts
> the count by exactly $+4 \bmod 8$ in all eight initial states. 128 runs recover the original
> phase, 8 shift by a legal multiple of the period, and none loses lock.

### `fig_n6_03_eight_initial_states.png` — 八初态一致性

| 面板 | 内容 |
|---|---|
| **(A)** | 000–111 八个初态的解码序列，逐格显示，可见彼此是**同一序列的纯旋转** |
| **(B)** | ⚠️ **配对指标 vs 已退休的混叠指标**（见下） |
| **(C)** | 每个初态下三个 bit 的读取裕量 |

**核对结果**：八初态**全部 certified**；配对指标 median `0.0142408`、相对极差 **0.1958 %**；每格 56 个读窗、commitment 全为 1.0。

> ⚠️ **面板 (B) 故意把陷阱画出来**：`eight_initial_verdict.json` 里每个初态的
> `leak5pct_ratio` 是**已退休的单窗口（混叠）指标**，它的八初态范围是
> min `0.0049` / median `0.0277` / **max `34705.5`** —— 同一个电路显示出 **3.5×10⁴ 倍的假差异**。
> 若照这个字段作图，会得出"八初态极不一致"的**完全错误**的结论。
> 正确的配对指标在 `carry_pairing_eight_verdict.json`，相对极差 **0.196 %**。

> **Caption draft.** The eight phase-consistent initial states. (A) Each decodes to a pure
> rotation of the same strictly mod-8 sequence; all eight are certified. (B) The retired
> single-window leak ratio would suggest a $3.5\times10^4$-fold spread within one and the same
> circuit; the paired metric spans 0.196 %. (C) Per-bit read margins, all far from the band edges.

### `fig_n6_04_strict_tolerance.png` — 严格容差（5 点 × 3 档）

| 面板 | 内容 |
|---|---|
| **(A)** | 5 点 × 3 档的认证与穿越数；**明确标出 n=4 是故意设置的阴性对照** |
| **(B)** | 各报告指标相对 baseline 的最大变化（对数轴），只统计通过点 |
| **(C)** | `prefix_traj_gap`（51 状态 vs 冻结 34 状态前缀）随容差档位 |
| **(D)** | 七个 `*_identical` 标志的满足比例 |

**核对结果**：

- **15/15 全部完成且稳定**（`ok = 15/15`，`all_points_stable = true`）；
- **`certified_failure` 是 n=4 的阴性对照**，三档下**始终失败**（crossings 5/14、17 个未标读窗）。这正是该实验的价值：**收紧积分器既救不活失败点，也弄不坏通过点**；
- 指标最大变化 `|5pct_L_symmetric_delta| ≤ 3.7×10⁻¹⁰`；
- `prefix_traj_gap`：baseline 为 0（它就是参考），tight `4.19×10⁻⁷`、ultra `4.21×10⁻⁷` —— **饱和而非随容差增长**，说明它不是积分误差主导；
- 七个 `*_identical` 标志**全部 15/15**：码串、稳态码串、前缀码串、bit2 码串、事件、R/F pattern、认证判定**全一致**。

> **Caption draft.** Strict numerical tolerance at the frozen point. (A) Five parameter points
> at three tolerance levels; the $n=4$ point is a deliberate negative control and fails at every
> level, so tighter integration neither rescues a failure nor breaks a pass. (B) The largest
> change in any reported metric is about $10^{-10}$. (C) The frozen 34-state prefix gap
> saturates at $4.2\times10^{-7}$ rather than growing with tolerance. (D) Code string, events,
> R/F pattern and certification verdict are identical in 15/15 runs.

---

## 3. 运行

```powershell
$env:MPLBACKEND='Agg'
& 'D:\aconade\python.exe' .\make_n6_figures.py
```

- 只依赖 `numpy / pandas / matplotlib`；
- 配色优先从 `../figures/tempo_style.py` 读取，读不到时用同一组内置常量；
- 字号与线宽按**最终网页显示尺寸**反推（标签 ≥ 16 px、刻度 ≥ 14 px、标题 ≥ 18 px、线条 ≥ 2 px），与 `07_第3-4章_Wiki及作图要求.md` §2.2 一致；
- 每张图同时输出 **PNG (dpi 200)** 与 **PDF**。

---

## 4. 本目录**不**包含的 n=6 图（已在别处，无需重画）

| 内容 | 现成位置 |
|---|---|
| 600 h 全程 / 90–190 h 放大 / 三位数字条带 | `threebit51_results/wiki_n6_20260924/01–03` |
| R 型与 F 型进位并列 / 读窗 vs 读带 / R-F 超周期与 far-off 泄漏 | `threebit51_results/flip_diagnostics/fig1–fig4` |
| 后冻结峰相位补充（36 条 = 30 扰动 + 6 对照） | `threebit51_results/postfreeze_peak_perturbation/`（只有 CSV/JSON，未作图） |

---

## 5. 与两比特结果的边界

本目录全部是**三级（51 状态）**的 n=6 证据。两比特的参数分析与鲁棒性图在
`final_reconstruction/twobit34_results/`（75 点工作窗口、四组分子池扰动、四初态、严格容差），
请不要与这里的图混排成同一阶段——两边的判据与网格不同。
