# Bit0 单独作图（第一个 Bit）

> **范围声明**：本目录**只画第一个 Bit**。bit1 / bit2 / A0 / F0 / A1 / F1 / carry / g1 一律不画。
> 唯一例外是 bit0 的输入通量 `C31_flux_uM_h`——它是 bit0 的接口，不是下游量。
>
> **未改动任何现有文件**。本目录为新目录，脚本只向 `./out` 写入。

---

## 1. 数据来源与配置

| 项 | 值 |
|---|---|
| 轨迹 | `final_reconstruction/threebit51_results/wiki_n6_20260924/trajectories.csv` |
| 规模 | 18001 采样 × 67 列，0 → 600 h，dt = 2 min |
| 读窗 | `.../wiki_n6_20260924/read_windows.csv`，**56 个已验证读窗** |
| 哈希清单 | `.../wiki_n6_20260924/SHA256SUMS.json` |
| 浓度换算 | `uM_per_au = 5.75` |
| bit mRNA 半衰期 | 2 min |
| bit 成熟半衰期 | 20 min |
| `add_growth` | `False`（gamma 已含稀释，不再叠加） |
| 时钟门 | **K = 0.3, n = 2**（实际成功扫描的配置，**不是** 0.4/3） |
| 上游 | 真实 ϕC31 production flux，µM/h |

**为什么用这一份**：bit0 的前 34 状态在 34 状态与 51 状态模型间是严格前缀（RHS 误差 0、轨迹误差 0），所以这份 600 h 轨迹里的 bit0 与冻结 34 状态完全一致，同时提供更长的统计样本（56 拍 vs 300 h 的 28 拍）。

**bit0 用到的列**：

```text
输入      C31_flux_uM_h
DNA       b0_S          (= S0)
通量       J_fwd0, J_rev0
重组酶     b0_M_I, b0_I_u, b0_I
BM3R1/T   b0_M_T, b0_T_u, b0_T
RDF       b0_M_R, b0_R_u, b0_R
复合物/状态 b0_C, b0_S
```

---

## 2. 数据自检（跑完即验，全部通过）

| 项 | 结果 |
|---|---|
| 解码交替 | 28 个 1 / 28 个 0，完美交替 |
| 周期 | **10.591 ± 0.015 h** |
| 读窗内 S 范围（高态） | 0.9714 – 0.9830 |
| 读窗内 S 范围（低态） | 0.0744 – 0.0760 |
| 到读带边界的裕量 | min **0.2240** / median 0.2482 / max 0.2785；**0 个窗口贴近边界** |
| 通量峰 / 谷 | 6.52–6.64 / 0.03597 µM/h |
| 每周期通量剂量 | 22.47 – 22.90 µM |

> 注意驻留电平是 **0.075 / 0.980**，与早期 bit0 诊断报的 0.18 / 0.84 不同——那是
> `uM_per_au = 10` 这一**旧配置**下的结果。冻结的 5.75 给出远干净的两条带。

---

## 3. 九张图

### `fig_bit0_01_one_cycle_decomposition.png` — 单周期翻转分解
对应 `04_第3章_作图说明.md` 的 **Figure 3-2**。四行同步：DNA 状态（带读带底色）／RDF 池与复合物（各自归一化，图例给真实峰值）／正反向通量／时钟输入。橙点线标出翻转时刻。
**用途**：第 3 章核心机制图——一眼看到"正反两个方向在时间上分开"。
> **Caption draft.** Bit 0 over one clock cycle (period 10.59 h). The DNA state is shown
> against the read bands; the RDF pool and the Int–RDF complex are each normalised to their own
> maximum; forward and reverse recombination fluxes are plotted on a common axis; the bottom
> panel is the ϕC31 production flux that drives the bit. Forward and reverse events are
> separated in time, and the flip occurs on the pulse. Model output, frozen configuration.

### `fig_bit0_02_dwell_and_readwindow.png` — 驻留带与读窗（三联）
对应 **Figure 3-3**。**(A)** 读窗几何：谷值居中、占周期 20%；**(B)** 全程 600 h 的 S，叠加 56 个读窗（按解码值着色）与解码点；**(C)** 每个读窗内 S 的 min–max 条（蓝=0，绿=1），显示**没有一个窗口碰到中间灰区**。
> **Caption draft.** (A) Each clock period carries one read window centred on the flux trough,
> spanning 20 % of the period. (B) The full 600 h trajectory with all 56 validated read windows;
> every window is labelled. (C) The range of S inside each read window: the two sets never
> overlap and no window enters the ambiguous 0.30–0.70 band.

### `fig_bit0_03_eleven_states.png` — 十一个 bit0 状态
四行：重组酶链 `M_I→I_u→I`／BM3R1 链 `M_T→T_u→T`／RDF 链 `M_R→R_u→R`（各自归一化，图例给真实峰值）／复合物 C 与 DNA 状态 S（双轴）。
**用途**：状态总览，也可作 Methods 折叠栏的配图。
> **Caption draft.** All eleven continuous states of bit 0 over one cycle. Each expression
> chain is resolved into mRNA, immature and mature protein and normalised to its own maximum;
> the true maxima are given in the legend. The bottom panel shows the Int–RDF complex and the
> DNA state on separate axes.

### `fig_bit0_04_flux_separation.png` — 正反向通量分离
三联：**(A)** 几个周期的正反通量（竖虚线为读窗时刻）；**(B)** 单周期放大；**(C)** 同一周期的 DNA 状态 + 读带。
**用途**：机制论证的核心证据——**每个周期内两个方向不重叠**。
> **Caption draft.** Forward and reverse recombination fluxes of bit 0. (A) Several cycles:
> the two directions never overlap. (B) One cycle enlarged. (C) The DNA state over the same
> cycle, with the read bands shaded.

### `fig_bit0_05_phase_portrait.png` — 相图
三联：**(A)** `S` vs RDF 池（按时间着色）；**(B)** `S` vs **同一时刻的乘积** `I·R`；**(C)** `S` vs 复合物 `C`。
**用途**：展示"方向记忆寄存在 RDF 池"以及"反向速率取决于乘积而非 R 单值"。**(B)** 直接对应第 3 章必须讲准的那句话。
> **Caption draft.** Phase portraits of bit 0 over 600 h. (B) plots the DNA state against the
> instantaneous product $I\cdot R$: the reverse reaction depends on the product rather than on
> the free RDF pool alone.

### `fig_bit0_06_cycle_metrics.png` — 逐周期指标（四联）
**(A)** 每个读窗的 S 值（始终在带内）；**(B)** 每周期通量峰/谷；**(C)** 逐周期周期长度（中位 10.60 h）——**周期由时钟决定，不由 bit 决定**；**(D)** 每周期正/反通量积分剂量交替。
> **Caption draft.** Per-cycle metrics. (C) shows the trough-to-trough period stays at
> 10.60 h: the bit follows the clock, it does not set the period. (D) forward and reverse flux
> doses alternate cycle by cycle.

### `fig_bit0_07_expression_chains.png` — 表达链延迟
**(A)** 三个成熟蛋白在同一周期内的归一化轨迹；**(B)** RDF 链展开成 mRNA → 未成熟 → 成熟，显示成熟级带来的延迟。
**用途**：解释 bit0 的记忆与延迟为什么存在。
> **Caption draft.** (B) The RDF expression chain resolved into mRNA, immature and mature
> protein within one cycle; the maturation stages are what delay the pool relative to
> transcription.

### `fig_bit0_08_commitment_statistics.png` — 读取承诺统计
**(A)** 56 个读点 S 的直方图（双峰，**中间灰区从未被占据**）；**(B)** 每个读窗到最近读带边界的裕量，排序后展示——min 0.224。
> **Caption draft.** (A) Histogram of S at the 56 read instants: strictly bimodal, the
> ambiguous band is never occupied. (B) Margin from each read window to its nearest band
> boundary, sorted; the minimum is 0.224, so no read is close to a decision edge.

### `fig_bit0_09_early_window_map.png` — bit0 单独早期窗口扫描
**(A)** `uM_per_au × bit 成熟半衰期` 的通过/未通过网格，格内标出解码码串（`·` = 未标读窗）；**(B)** 高态驻留电平 `S_hi` 随换算标尺下降（0.980 → 0.727），**成熟半衰期几乎不影响**。

> ⚠️ **这张图的判据与网格必须随图一起说明**：
> 1. 用的是**早期驻留标签判据 `stored_bit`**，**不是**有限读窗的 commitment 判据；
> 2. 网格 `uM_per_au` 从 **6.0 起**，而冻结中心 **5.75 在这个网格之外**；
> 3. 因此它是"高态电平为什么会在某个换算标尺上失守"的**机制诊断**，不能当作当前工作窗口。
> 4. 它来自 `bit0_results/verification/bit0_region_map.csv`（56 点，8 × 7）。

> **Caption draft.** Bit 0 alone over the concentration-scale factor and the maturation
> half-life, scored with the early dwell-label criterion. The high-state dwell level falls
> monotonically from 0.980 at scale 6 to 0.727 at scale 18, while the maturation half-life
> barely moves it; above scale 10 the read windows become unlabelled or the bit latches high.
> Diagnostic scan only: the grid starts at 6.0, above the frozen centre of 5.75.

---

## 4. 附带的指标表

`out/bit0_cycle_metrics.csv`（56 行 × 13 列），逐拍可复算：

```text
index, trough_h, decoded_bit0, S_at_read, S_window_min, S_window_max,
flux_peak_uM_h, flux_trough_uM_h, flux_dose_uM,
forward_dose, reverse_dose, flip_time_h, period_h
```

---

## 5. 运行

```powershell
$env:MPLBACKEND='Agg'
& 'D:\aconade\python.exe' .\make_bit0_figures.py
```

- 脚本自包含，只依赖 `numpy / pandas / matplotlib`；
- 配色优先从 `../figures/tempo_style.py` 读取；该文件不存在时使用同一组内置常量，因此不会因路径问题失败；
- **字号与线宽按最终网页显示尺寸反推**（标签 ≥ 16 px、刻度 ≥ 14 px、标题 ≥ 18 px、线条 ≥ 2 px），
  与 `07_第3-4章_Wiki及作图要求.md` §2.2 一致，不是 matplotlib 默认值；
- 每张图同时输出 **PNG (dpi 200)** 与 **PDF**。

---

## 6. 未包含的内容（有意排除）

| 内容 | 为什么不在本目录 |
|---|---|
| Figure 3-5 扰动与容差 | 扫描网格定义在两比特系统上（S0/S1 成对），不是 bit0 单独的量 |
| Figure 3-4 工作窗口（正式版） | 数据未冻结，需二维细化扫描；本目录只给了早期诊断版（fig09） |
| Figure 3-1 机制图 | 不依赖数据，交美工重绘 |
| `bit0_results/` 血统的轨迹 | 该线的"确认候选"在 `uM_per_au = 10`，其两比特测试 `bit1_S_range = [0, 0.0057]`（bit1 不翻），且 kon/koff 与 conversion 扫描的 `stored_by_commitment` 均为 **0**。只能作判据演进的历史证据，不能作当前 bit0 的机制图数据源 |
