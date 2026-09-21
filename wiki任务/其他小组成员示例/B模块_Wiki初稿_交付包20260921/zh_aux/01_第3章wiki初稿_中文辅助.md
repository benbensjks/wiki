# 第 3 章 Coupling the Clock to a Reliable Single-bit Counter
## 时钟耦合到可靠的单比特计数器（B 模块 wiki 初稿，2026-09-21）

> 对应组会 9.15 目录框架的第 3 章，回答该 PPT 第 3 页列出的 B 模块问题：
> BM3R1 启动子强度、BM3R1/RDF 关系、降解标签加在哪里、A 最新波形下旧结论是否成立、
> one-pulse-one-flip 工作窗口如何定义。**单比特 / 耦合 robustness 按组会要求放在本章。**
>
> **文档状态：wiki 初稿（未冻结）**。用于网页组先搭结构、美工先画机制图；
> 9.30 冻结前会随最终模型口径与实验标定持续修订。
> 数值来源：B 模块 Week 3 / Week 4 交付包与 2026-09-19 两份交付包（非零泄漏方案、启动子审查）。
> 未实验标定的一律按第 10 节的证据标签标注；本稿不把模型结论写成实验结论。

---

## 目录

1. [页面主线](#1-页面主线)
2. [建模完成状态](#2-建模完成状态)
3. [The Problem — 本章要回答什么](#3-the-problem--本章要回答什么)
4. [How the Switch Works — 单比特如何工作](#4-how-the-switch-works--单比特如何工作)
5. [One Pulse, One Flip — 判据与工作窗口](#5-one-pulse-one-flip--判据与工作窗口)
6. [Re-validation under the Final Clock — A 最新波形后的重验证](#6-re-validation-under-the-final-clock--a-最新波形后的重验证)
7. [The BM3R1–RDF Delay Circuit — 关系与表达比例](#7-the-bm3r1rdf-delay-circuit--关系与表达比例)
8. [Int Degradation Tag — 加在哪里、什么条件下需要](#8-int-degradation-tag--加在哪里什么条件下需要)
9. [Promoter Recommendation — RDF 端启动子](#9-promoter-recommendation--rdf-端启动子)
10. [Robustness and Boundaries — 鲁棒性与边界](#10-robustness-and-boundaries--鲁棒性与边界)
11. [Engineering — Model-guided Decisions（简版）](#11-engineering--model-guided-decisions简版)
12. [给实验组的指引（Guidance to Wet Lab）](#12-给实验组的指引guidance-to-wet-lab)
13. [Scope, Limitations and Reproducibility](#13-scope-limitations-and-reproducibility)
14. [本章结论](#14-本章结论)
15. [主图清单与图注草稿](#15-主图清单与图注草稿)
16. [参考文献](#16-参考文献)

---

## 1. 页面主线

```text
真实上游 ϕC31 波形（周期 10.59 h，产生通量接口，µM/h）
          ↓ 驱动
单比特开关：Int / RDF / BM3R1 + 一个 PB ⇄ LR 的 DNA 双状态
          ↓ 机制
一次脉冲 = 一次翻转：记忆寄存在 RDF 池，方向由脉冲到达时的 RDF 存量决定
          ↓ 判据
面向状态而非瞬间的严格计数判据（端点 Q + 每周期一次中点穿越 + 计数奇偶性）
          ↓ 设计
表达比例窗口（BM3R1 / RDF）× Int 降解标签窗口 × RDF 启动子泄漏容差
          ↓ 压力测试
启动相位、参数扰动、异质性蒙特卡洛、Int 清除时 RDF 命运
          ↓ 设计决策
启动子 66 bp 推荐、表达配比目标、标签梯度标尺、初始化条件
```

---

## 2. 建模完成状态

| 内容 | 状态 | 当前证据边界 |
|---|---|---|
| ϕC31–RDF–BM3R1 单比特模型（38 维 ODE） | 已完成 | 重组核心逐行移植自 Zhao 2019 原文（35 式），BM3R1 延迟电路 3 式替换；DNA 守恒至 5 位有效数字 |
| 一次脉冲一次翻转机制 | 已完成 | 失效诊断 + 工作点扫描；翻转方向由脉冲时刻的 RDF 存量决定 |
| A 最新波形（Week 4）下的重验证 | 已完成 | 旧默认参数失效、新工作点成立；K 敏感带消失 |
| Int 降解标签功能窗口与机制 | 已完成 | 确定性窗口、随机全稳窗口、上界碎裂带、文献速率核对 |
| RDF 启动子非零泄漏（0.8%）下的工作方案 | 已完成（模型） | 双初态 430 h、35 个稳态周期；固定相位扰动通过；任意相位仍有失败 |
| 可行域（表达比例 × 泄漏） | 已完成（离散网格） | β_R 4–8 × β_B 0.75–6 宽平台；泄漏容差约 ≤1% |
| 启动子序列与实验交接 | 已完成 | Cello pBM3R1 66 bp 推荐；36 bp 版本核验为 AI 生成且缺碱基 |
| Int 清除时 RDF 是否共丢失（η） | 未完成（实验判别） | η=1 与 η=0 两套参数不能混用；η≥0.9 通过、η≤0.75 失败 |
| 表达强度实验标定（α、β、拷贝数） | 未完成 | 现有 β 是模型翻译系数，不是启动子 RPU，不唯一对应 RBS 序列 |
| 分子级噪声 / 单细胞成功率 | 未完成 | 现有 MC 是参数情景扰动，不等于 Gillespie 单细胞模拟 |

---

## 3. The Problem — 本章要回答什么

> 一个用位点特异性重组实现的双稳态开关，在**真实上游波形**（而不是理想方波）驱动下，
> 能不能做到"每个时钟脉冲恰好翻转一次"？判定它"能"的标准是什么？工作窗口有多宽？
> 哪些设计量真正决定成败？

四个子问题：

1. **机制**：为什么它不是"输入来了就跟着走"，而是每个脉冲只翻一次？
2. **接口**：从 A 模块传过来的到底是什么量？旧波形换成最新波形后，原来的结论还成立吗？
3. **判据与窗口**：$S$ 是连续量，凭什么读成 0 / 1？one-pulse-one-flip 的工作窗口如何定义？
4. **设计量**：BM3R1/RDF 表达比例、Int 降解标签、RDF 启动子泄漏各自要求什么范围？

---

## 4. How the Switch Works — 单比特如何工作

### 4.1 状态变量与 DNA 双状态

计数器基本单元是一段两侧为 ϕC31 `attP`/`attB` 位点的可翻转 DNA：

- **PB 态（bit 0）**：`attP × attB`；LR 比例记作 $S$，则 PB 比例为 $1-S$；
- **LR 态（bit 1）**：`attL × attR`。

$S$ 是连续比例而不是逻辑位：模型中没有阈值事件、没有事件队列、没有"空闲期强制冻结"，
所有数字行为都必须是连续动力学的结果。单比特模型共 38 个 ODE 状态：
Int / BM3R1 / RDF 各自的 mRNA–未成熟–成熟链、游离 Int 与二聚体、RDF、
Int–RDF 复合物、以及全部 DNA–Int(±RDF) 中间复合物与突触复合物。

### 4.2 两个方向的重组

$$v_f = k_{fwd}\,H(I;K_{D,int},2)\,\frac{K_{inh}}{K_{inh}+R},
\qquad
v_r = k_{rev}\,H(C;K_C,2)$$

$$\frac{dS}{dt} = v_f\,(1-S) - v_r\,S$$

其中 $H(x;K,n)=x^n/(K^n+x^n)$；$I$、$R$ 分别为游离 Int 与游离 RDF，
$C$ 为 Int–RDF 复合物。**正向（PB→LR）由游离 Int 驱动；反向（LR→PB）由 Int–RDF 复合物驱动。**

> 关键机制：反向重组的自变量是**同一时刻的乘积** $I\cdot R$（复合物准稳态下
> $C \approx q\,I\,R$，等价于 $H(I\cdot R;K_{D,comp},2)$），不是 $R$ 本身；
> 且反向重组**不消耗** RDF——自由 RDF 的下降来自结合隔离与复合物清除。

### 4.3 BM3R1 延迟电路

$$ \dot m_B = \alpha_B D_{PB} - \gamma_B m_B, \qquad \dot B = \beta_B m_B - \mu B $$

$$ \dot m_R = \alpha_R D_{LR}\left[\ell_R + \frac{1-\ell_R}{1+(B/K_B)^n}\right] - \gamma_R m_R $$

- BM3R1（记作 $T$）**只在 PB 态表达**；翻转后其 mRNA 与蛋白逐步消退；
- RDF **只在 LR 态转录**，并被 BM3R1 以 Hill 形式抑制；$\ell_R$ 是不可被无限强阻遏压低的泄漏底限。

延迟时长近似 $\tau \approx \ln(B_{ss}/K_B)/\mu$（仅生长稀释清除时）。
**延迟必须盖过整个脉冲宽度**，否则 RDF 会在同一个脉冲内解锁，DNA 来回双翻，计数失败——
这是本模块最核心的设计约束。

### 4.4 为什么能"一次脉冲一次翻转"

1. **RDF 池就是记忆载体。** LR 驻留期 RDF 持续累积；PB 期 BM3R1 把 RDF 产生压到近零，
   池被清空。在非零泄漏主要方案（430 h、双初态）中，**脉冲到达前**的 RDF 游离池为：
   LR 态约 **4.46 µM**、PB 态约 **0.0018 µM**，相差约 **2500 倍**（三个数量级）。
2. **脉冲到达时，RDF 存量决定方向。** 上一状态为 LR 时复合物大、反向重组胜出；
   上一状态为 PB 时 RDF 池已空，正向重组没有对手。
3. **翻转是单方向的单调过渡**：LR→PB 与 PB→LR 各自在约 2 h 内单调完成
   （LR 从 ~1.0 单调降到 ~0.003，或反向），同一脉冲内没有来回双翻——
   方向在脉冲到达前就由 RDF 池存量决定。

> 模型中的 $dS/dt$ 只依赖当下的 $I,R,C$，没有延迟项或双稳项——
> 记忆全部寄存在 RDF 池里，而 RDF 池在两次时钟脉冲之间被保留。

### Figure 3-1（机制图，交美工）

状态依赖的开关：PB / LR 两态、各自的生产分配、脉冲到达时的 RDF 池量与两个翻转方向。
构图、箭头语义、禁止画错的点与英文标签见配套《02_机制图需求与草图》；
本轮草图 `figures/fig03_1_single_bit_mechanism_sketch.png`。
One-pulse-one-flip 时序示意图作为配套概念图保留，待日后细化。

---

## 5. One Pulse, One Flip — 判据与工作窗口

### 5.1 严格计数判据

$S$ 是连续量，必须公开"怎么读"。本章采用的判据：

1. **端点标准**：高低端点分别 $\ge 0.95$、$\le 0.05$；逐点取最差正确状态比例 $Q$，
   不能跳过中间态；
2. **每周期恰一次中点穿越**：用积分器内部步长上的**求根事件**检测 $LR=0.5$ 穿越，
   而不是只看稀疏作图点；
3. **计数奇偶性**：在明确初态与启动相位时，核对从第一个完整输入开始的计数奇偶性与全部穿越次数。

$Q$ 是 DNA 端点状态指标，**不是"单细胞计数成功率"**。

> 声明：这是工程验收判据，不是实验标准。判据本身会改变结论，因此必须公开定义。

### 5.2 工作窗口（A 最新波形、零 RDF 泄漏假设下）

以 Week 4 定稿输入（周期 10.59 h、产生通量接口、PLtetO1 泄漏实测 0.5%）驱动，
未调参的原文默认参数与 BM3R1 未调参组合全部失效；重扫后的工作点：

| 元件 / 参数 | 推荐值 | 窗口 |
|---|---|---|
| ϕC31 RBS scale（A 输入接口） | 0.45 | — |
| $k_{BM3,tsl}$（BM3R1 翻译） | 15 h⁻¹ | 6–30 h⁻¹ |
| $k_{rdf,tsl}$（RDF 翻译） | 200 h⁻¹ | ≥200 h⁻¹ 稳健（100 可用） |
| $K_{BM3}$ / $n$ | 18.6 nM / 3.4 | 12–50 nM × 2.9–3.4 全过 |
| $k_{tag,int}$（Int 降解标签） | 8–16 h⁻¹ | 确定性 [3, 24]，随机全稳 6–18 |

双初态 score 0.995–0.998，稳态谷采样序列 0.00/1.00 干净交替、fidelity 1.00。

> **K 敏感带消失**是本次重跑最重要的鲁棒性结论之一：v36 时代工作带仅 18–35 nM
> （$K\le12$ 崩塌、$K\ge50$ 失效），新波形（周期 10.59 h、脉冲间隙约 8 h）下
> 12–50 nM 任意 $K$、两种 Hill 构型均可靠计数——实验组对 $K_{BM3}$ 精确取值不再敏感。

### 5.3 非零 RDF 泄漏（0.8%）下的当前方案

启动子审查发现 RDF 启动子存在不可忽略的泄漏底限（Cello B1 门 $y_{min}/y_{max}=0.8\%$）。
在显式非零泄漏下，Week 4 原工作点失效；重新匹配表达 / 输入 / 清除后找到方案：

| 量 | Week 4 原工作点 | 非零泄漏主要方案 | 属性 |
|---|---:|---:|---|
| BM3R1 翻译系数 $\beta_B$ | 15 h⁻¹ | **2 h⁻¹** | 模型目标 |
| RDF 翻译系数 $\beta_R$ | 200 h⁻¹ | **8 h⁻¹** | 模型目标 |
| Int 合成通量 scale | 0.45 | **0.30** | 相对 A 冻结输入 |
| Int 额外清除 $k_I$ | 12 h⁻¹ | **12 h⁻¹** | 有效模型参数 |
| RDF 启动子泄漏底限 $\ell_R$ | 0（假设） | **0.008** | Cello B1 文献先验情景 |
| Int 清除后完整 RDF 返回比例 $\eta$ | 0（代码隐含） | **1** | 工作假设，需实验核验 |

**结果**：双初态 430 h、35 个稳态完整周期严格通过，最差 $Q=0.9966$；
已知初态、低谷启动下 40 个完整脉冲全部正确计数；固定相位扰动 200/200（较小）、
199/200（较宽）；任意相位 37/40——**初始化条件不能省**。

**可行域**（主要方案条件、固定低谷启动、430 h）：

- $(\beta_B,\beta_R)$ 宽平台：$\beta_R$ 4–8 × $\beta_B$ 0.75–6 全部通过（60 格）；
  含平台外扩展共 76/130 格通过；$\beta_R\ge16$ 在本次条件下无解；
- 泄漏容差：$\beta_R=8$ 时 $\ell_R\le0.01$ 通过（0.0125 起失败）；
  $\ell_R\le0.006$ 时 $\beta_R$ 4–16 全部通过；泄漏越高必须配越低的 $\beta_R$；
- 边界：网格为离散分辨率；固定低谷启动，随机相位会收窄；scale、tag 与机制 $\eta$ 固定，
  变更后需重扫。

**Figure 3-2**（概念草图）解释窗口"为什么有边"：两侧各自的失效机制
（多翻 / 回翻失败 / 泄漏放大的自发翻转），以及泄漏升高时窗口向低 RDF 强度移动。

### 5.4 启动条件（初始化是设计条件）

| 情景 | 通过 / 测试数 | 说明 |
|---|---:|---|
| 20 个启动相位 × PB/LR | **37/40** | 任意相位接入仍可失败 |
| 较小扰动，随机相位 | 193/200 | ±10% 情景集 |
| 较小扰动，固定低谷启动 | **200/200** | 最差 $Q=0.9933$ |
| 较宽扰动，随机相位 | 194/200 | ±20% 情景集 |
| 较宽扰动，固定低谷启动 | 199/200 | 余下一例 $Q=0.9249$ |

结论：**已知初态 + 从通量低谷接入第一个完整脉冲**是本方案的启动要求；
"完全不受接入相位影响"尚未达到。

---

## 6. Re-validation under the Final Clock — A 最新波形后的重验证

### 6.1 接口语义：传的是产生通量，不是浓度

A→B 的接口是 **ϕC31 产生通量** $J_I(t)=s_I \cdot v_A(t)$（µM/h，$s_I$ 为 RBS 相对翻译强度）：

$$ \dot I_{tot} = J_I(t) - (\mu + k_I)\,I_{tot} $$

Int 的结合、稀释与标签清除由 B 模块自行处理，**不再把蛋白浓度当作外部状态强迫**
（避免标签/稀释被重复积分）。这是全项目接口统一要求（A 给 B 的究竟是浓度、通量还是归一化波形）
在 B 侧的落实。

### 6.2 波形规格（Week 4 定稿输入，RBS=1 基准）

| 指标 | v36（旧，Week 3 用） | Week 4 unloaded（新） |
|---|---|---|
| 周期 | 6.77 h | **10.59 h** |
| C31 蛋白峰 | 3.99 µM | 7.01 µM |
| C31 蛋白谷 | 0.28 µM | 0.104 µM |
| 通量峰 / 谷（RBS=1） | —（蛋白反推） | 6.64 / 0.0360 µM/h（峰谷比 185） |
| 每周期通量剂量（RBS=1） | — | 22.8 µM（稳态谷–谷积分） |
| PLtetO1 leak（模型） | 5%（假设） | **0.5%**（谷值实测） |

非零泄漏方案取 scale=0.30：通量峰约 1.992 µM/h、谷约 0.0108 µM/h、
每周期通量积分约 6.853 µM。

### 6.3 重验证结论

1. **旧结论在旧输入下不成立、在新输入下重新成立**：未调参组合在新输入下依然失效；
   重扫工作点（5.2）双初态 score 0.995–0.998。
2. **失败机制可分解**（Week 3 已定位，新输入下缓解）：
   旧波形基线泄漏导致脉冲间隙自发翻转；脉冲太宽、延迟太短导致同脉冲回翻。
3. **接口需求**：旧波形要求 PLtetO1 泄漏 $\le1.5\%$；新波形实测 0.5%，已满足。
4. **但非零 RDF 泄漏又把 Week 4 点推翻**——这是当前方案（5.3）产生的原因。
   两轮"失效 → 重新匹配"说明：单比特的工作点不是一组固定常数，
   而是随输入波形与启动子泄漏共同移动的联合窗口。

---

## 7. The BM3R1–RDF Delay Circuit — 关系与表达比例

### 7.1 BM3R1 / RDF 的关系到底是什么

- BM3R1 的作用是**产生合适的 RDF 表达延迟**，其表达**不是越强越好**：
  - 太强 → RDF 解锁过迟 → 回翻失败（LR 卡高态）；
  - 太弱 → 延迟不足 → 同一脉冲内来回双翻。
- RDF 既可能干扰下一次正向翻转，也必须足量支持反向翻转；
  **不能把它当作单纯的"抑制剂"清除**。
- 延迟由 BM3R1 的清除速度决定：$\tau \approx \ln(B_{ss}/K_B)/\mu$。

### 7.2 表达比例窗口（转成实验可理解的目标）

非零泄漏主要方案在模型条件（$\alpha_B=\alpha_R=120\ \mathrm{h^{-1}}$、
$\gamma=4\ \mathrm{h^{-1}}$、$\mu=0.832\ \mathrm{h^{-1}}$、$D_{tot}=0.017\ \mu M$、
1 fL 假设）下换算：

| 指标 | 模型目标 | 可测形式（建议） |
|---|---:|---|
| $\alpha_B\beta_B$ | 240 h⁻² | — |
| $\alpha_R\beta_R$ | 960 h⁻² | — |
| BM3R1 最大合成通量 | 1.02 µM/h | PB 稳态 1.226 µM ≈ **740 拷贝/细胞** |
| RDF 最大合成通量 | 4.08 µM/h | 去阻遏稳态 4.905 µM ≈ **2950 拷贝/细胞** |
| BM3R1 : RDF 最大产能比 | 1 : 4 | 同一宿主/载体/生长条件比较 |
| RDF 相对泄漏 | ≤0.8% | 饱和阻遏与完全去阻遏荧光之比 |

**边界声明**：$\beta$ 是翻译系数，不是启动子强度（RPU），也不唯一对应一条 RBS 序列；
"15→2、200→8"是模型工作点之间的比较，不是对实验组现有构件实测表达量的倍数判断。
实际 $\alpha$、$\gamma$ 或拷贝数不同时，应以测得的产能与动态重新拟合。

### 7.3 启动子与 K/n 的证据级别

- 参考启动子为已表征的 Cello **pBM3R1**（66 bp），见第 9 节；
- $K_{BM3}=18.6$ nM 是**锚定估计**：由 Cello 实测泄漏 $y_{min}/y_{max}=0.008$
  与原文模型标度反推（$S/K=4.12$，取 $S=0.0765\ \mu M$），落在 TetR 家族
  体内半抑制浓度典型区间；Cello 的 $K=0.04$ RPU 不能直接换成 nM；
- 工作带对 $K$ 的宽度已由 5.2 / 5.3 给出；若实验实现为 B2 构型（$n=2.9$），
  最适 $K$ 需相应重选。

---

## 8. Int Degradation Tag — 加在哪里、什么条件下需要

### 8.1 对象与作用

**标签加在 Int 上**（所有含 Int 的物种：游离 Int、二聚体、Int–RDF 复合物及一切 DNA–Int 复合物）；
本方案**不需要**给 BM3R1 或 RDF 加降解标签。

标签的两种作用：

1. **压缩脉冲间隙游离 Int**：tag 0→4 使间隙游离 Int 降约 56 倍
   （0.0560 → 0.0010 µM），防止 RDF 解锁后的自发回翻与多翻——对应 tag 过低的失效；
2. **削脉冲 Int 峰**：峰 Int 随 tag 单调下降；削得太狠则 LR→PB 回翻失败——对应 tag 过高的失效。

### 8.2 功能窗口与失效形态

以 RBS 0.45 为基点：

- **确定性干净窗口 $k_{tag,int}\in[3,24]$ h⁻¹**（双初态 + 无中间态样本口径；
  下缘 2.8 起双初态通过，25 起碎裂）；
- **随机全稳窗口 6–18 h⁻¹**（30 实例/档的异质性蒙特卡洛：tag 6–18 全稳 30/30；
  tag 5 与 22 边缘；tag 3 与 25–26 明显退化）；
- **RBS × tag 联合工作区**呈"平台 + 两侧收窄"形：全稳平台 tag 7–18 在 RBS 0.30–0.60 全域通过；
  低速标签（$k\le4$）不能配高 RBS，高速标签（≥19）需相应上调 RBS 补偿；
- **上界不是单调失效**：tag 25–28 是碎裂的退化带（相邻档位实绩可跳变，
  失效形态以 LR→PB 回翻失败、周期卡滞为主），tag ≥30 后多数实例卡在 LR 高态。

### 8.3 标签速率必须实测

文献中同一 ssrA 标签的速率数据因测量口径不同相差很大：

- 原生 ssrA：E. coli 直测表观 $t_{1/2}\approx6$–14 min（$k\approx3$–8 h⁻¹），
  覆盖窗口下缘；
- LAA-LAA 等高速档：按 SI 残留比例换算约 14–38 h⁻¹，部分越过 RBS 0.45 上缘
  （可借 RBS 上调补偿，但 ≥22 后异质性下已不稳定）；
- 若 Int 融合使速率降到 1–2.7 h⁻¹（悲观情形），系统失效——
  此时换标签序列无济于事，必须考虑 ClpX 共表达。

**实验建议（梯度标尺法）**：构建同启动子/同 RBS 的
sfGFP–{无标签对照, 原生 ssrA, LAA+4, SsrA2X, LAA-LAA} 五个梯度；
在 MC4100、37 °C 指数期、低诱导浓度下用去诱导剂与翻译抑制双法交叉测 $k$；
按实测 $k$ 对窗口落档选型。

### 8.4 与输入强度的配对

标签速率不必恰好等于 12 h⁻¹。非零泄漏方案测试了近似维持 Int 有效幅度的接口族：

$$ s_I(k_I) = 0.30\,\frac{\mu+k_I}{\mu+12} $$

| 实测后回代的 $k_I$ / h⁻¹ | 对应输入 scale 候选 |
|---:|---:|
| 2 | 0.0662 |
| 4 | 0.1130 |
| 8 | 0.2065 |
| 12 | 0.3000 |
| 20 | 0.4870 |

这是**若干已测试离散点的配对规则，不是连续区间的数学证明**；
标签需根据 Int 融合蛋白实测速率匹配输入强度，不能把某个 ssrA 名称固定等同于 12 h⁻¹。

---

## 9. Promoter Recommendation — RDF 端启动子

### 9.1 推荐元件

可追溯的参考元件是 Cello **pBM3R1（66 bp）**，沿转录方向 5′→3′：

```text
AATCCGCGTGATAGGTCTGATTCGTTACCAATTGACGGAATGAACGTTCATTCCGATAATGCTAGC
```

| 参考门 | 去阻遏输出 | 受抑制输出 | 相对底限 |
|---|---:|---:|---:|
| Cello B1_BM3R1 | 0.5 RPU | 0.004 RPU | **0.8%** |

以上为 Cello 特定宿主/表达盒条件下、相对 J23101 参考盒的测量标尺，
**不是 MC4100 的绝对强度**。B1/B2/B3 使用同一启动子、不同 RBS，
不能把它们当作三条独立的低泄漏 promoter。

### 9.2 36 bp 版本的核验结论

实验组提供的 36 bp 序列经核验：

```text
DeepSeek 原始输出 (37 bp)：TTGACACGGAATGAACGTTCATTCCGATAATGCTAGC
实验组转述     (36 bp)：TTGACACGGAATGAACGTTCATCCGATAATGCTAGC
按同设计原则修正 (38 bp)：TTGACACGGAATGAACGTTCATTCCGTATAATGCTAGC
```

- 来源为 AI 生成（"J23119 + BM3R1 operator"），**不是文献或 Registry 中的现成元件**；
- 36 bp 版本把 20 bp 完美回文 operator 破坏为 19 bp 非回文；
  36/37 bp 两个版本都相对 J23119 的 -10 区（`TATAAT`）少一个 T；
- 不能套用 Cello/Stanton 的强度、泄漏或 $K/n$ 参数。

**建议：不要按 36/37 bp 原样合成。** 首选 Cello pBM3R1 66 bp；
若走 J23119 路线，用修正后的 38 bp 并作为未表征元件实测。

> 方向性提醒：非零泄漏方案要求 RDF 表达**中等且需标定**；
> "选最强启动子"来自零泄漏旧工作点，与当前结论方向相反。

---

## 10. Robustness and Boundaries — 鲁棒性与边界

### 10.1 数值与模型边界

- DNA 总量守恒：全程保持 0.01700 µM（5 位有效数字）；
- 收紧积分容差与最大步长后结论一致（详见交付包 `data/nominal_comparison.csv`）；
- 现有全部结论是**确定性平均场**：没有模拟分子级噪声、质粒分离或单细胞实测成功率；
  异质性 MC 是参数情景扰动，不是 Gillespie 单细胞模拟。

### 10.2 机制不确定性：η（Int 清除时 RDF 是否共丢失）

主要方案假设 Int 被标签清除时，非共价结合的 RDF 可以完整返回（$\eta=1$）；
Week 4 代码隐含 RDF 随含 Int 复合物共同丢失（$\eta=0$）。两套处理不能混用参数：

- $\eta=1$：主要方案通过；
- $\eta=0.9$：仍通过（$Q\approx0.983$）；
- $\eta\le0.75$：已测档位不通过；
- 保留原 $\eta=0$ 处理也存在解：$\beta_B=1.5$、$\beta_R=20$、scale=0.085、
  $k_I=2.2$，最差 $Q=0.9879$。

**选择哪种处理，应回到"Int 清除时是否额外消耗 RDF"的实测证据**
（同条件 RDF 稳定性/含量变化）。

### 10.3 输入对照与 pBAD-INT

针对主要 B 参数，另测试了"合成通量峰 2 µM/h、相对底限 0.5%"的平滑方波：
周期 10.6 h、宽 2 或 3 h 时双初态通过；宽 0.25/0.5/1 h 未达到严格端点标准。
这是用来规划 **pBAD-INT 分步表征**的通量条件对照，
**不是阿拉伯糖浓度或加糖时长的直接操作规格**。

### 10.4 未完成 / 未建模清单

1. 表达强度标定（$\alpha$、$\beta$、拷贝数）与 Int 融合蛋白清除速率；
2. RDF 共丢失判别（$\eta$）；
3. 分子噪声与单细胞成功率（Gillespie / 质粒分离）；
4. 高 tag 档（≥22）在 RBS 补偿配置下的随机稳健性复测；
5. A+B 联合的资源负载自洽审计（本模块结论基于 A 的 unloaded 输入）。

---

## 11. Engineering — Model-guided Decisions（简版）

| 环节 | 内容 |
|---|---|
| **Design** | 用 Int / RDF / BM3R1 三元件与一个 DNA 双状态实现一位存储；BM3R1 提供延迟，RDF 池承载记忆；不引入逻辑位或事件式翻转 |
| **Build** | 38 维显式 ODE（原文重组核心 + BM3R1 延迟电路 + 显式泄漏 + Int 标签清除）；上游真实通量作为唯一驱动 |
| **Test** | 严格计数判据、双初态长程（430 h）、启动相位与参数扰动、异质性 MC、启动子序列审计 |
| **Learn** | 记忆在 RDF 池；BM3R1 不是越强越好；标签窗口两端都有边界且上界非单调；初始化与表达标定是最关键的两个实验前置 |
| **Decision** | 采用 66 bp 参考启动子；给出 BM3R1/RDF 表达配比与 Int 清除的配对规则；把"低谷启动"写进实验条件；将 η 判别列为下一步实验 |

> 完整的 Design → Build → Test → Learn → Redesign 叙事见配套《03_Engineering构想》。

---

## 12. 给实验组的指引（Guidance to Wet Lab）

本章"模型 → 实验"决策，按实验组应执行的顺序排列：

| # | 决策 | 实验动作 | 状态 |
|---|---|---|---|
| 1 | **RDF 启动子** | 使用已表征的 Cello pBM3R1 66 bp；**不要**按 36/37 bp 原样合成；若走 J23119 路线，用修正的 38 bp 并作为未表征元件实测 | 已由审计决定 |
| 2 | **表达标定** | 在同宿主/载体中测当前 BM3R1 / RDF 构件的稳态表达与泄漏；换算成 µM 与拷贝数；对照目标（BM3R1 ≈ 740 拷贝/细胞、RDF ≈ 2950 拷贝/细胞、比例 ≈ 1:4）调 RBS 或启动子强度 | 实验第一优先级 |
| 3 | **Int 降解标签** | 构建 sfGFP 梯度标尺（无标签 / 原生 ssrA / LAA+4 / SsrA2X / LAA-LAA）；在 MC4100、37 °C 指数期、低诱导下用去诱导与翻译抑制双法测 $k$；按 $k$ 落档选型；若 $k<3\ \mathrm{h^{-1}}$ 考虑 ClpX 共表达 | 标签必需；速率必须实测 |
| 4 | **η（Int 清除时 RDF 命运）** | 在配平条件下比较有/无 Int 标签清除时的 RDF 稳定性/含量；若确有共丢失，切换到 η=0 备选参数组并重扫 | 最大机制不确定性 |
| 5 | **启动安排** | 提供已知初态，从通量低谷接入第一个完整脉冲；不要假设任意相位接入（37/40） | 当前方案的一部分 |
| 6 | **pBAD-INT 分步表征** | 先用方波通量对照（峰 2 µM/h、宽 2–3 h 通过）表征开关，再接入时钟；模型通量**不等于**阿拉伯糖浓度或加糖时长 | 已列入实验安排 |
| 7 | **回代重校** | 把实测值回代参数化扫描脚本，冻结构建选择前重核工作窗口 | 持续进行 |

## 13. Scope, Limitations and Reproducibility

### In scope

- Int / RDF / BM3R1 与 PB ⇄ LR 双状态的连续动力学；
- 真实上游 ϕC31 产生通量作为驱动；
- 严格计数判据、工作窗口、启动条件与扰动；
- Int 降解标签、BM3R1/RDF 表达比例、RDF 启动子泄漏与序列。

### Model boundary

- 本模块只负责**一位**；多比特进位见第 4 章（FFL 与游离重组酶架构）；
- 输出时长的编程由输出时长模块负责，本模块登记共享接口（BM3R1 池）；
- 细胞间差异、分子噪声与群体相位不在本章范围内。

### Parameter evidence labels

| 标签 | 在本章的对象 |
|---|---|
| literature value | 重组核心速率常数（Zhao 2019 / Pokhilko 2016）；Cello B1 门的 ymax/ymin |
| literature-informed range | 原生 ssrA 等标签的速率区间 |
| model-effective parameter | $\beta_B$、$\beta_R$、$K_B$、$n$、$\gamma$、$k_I$、$\ell_R$、$\eta$ |
| design variable | Int 输入 scale、RBS、启动子档位 |
| measurement pending | 表达标定、Int 融合蛋白清除速率、η、单细胞噪声 |

### Reproducibility

- 模型与扫描脚本随交付包提供（Week 4 交付包 + 非零泄漏交付包）；
- 计算环境：srv2026 `igem-tempo-2026`（Python 3.12；numpy/scipy/matplotlib）；
- 输出保存 CSV / JSON 与作图数据；扫描保留参数网格与元数据；
- 运行分级：`smoke` 流程检查、`scan` 统计、`confirmation` 报告结果。

---

## 14. 本章结论

1. 单比特在真实上游波形下可以成为**有记忆的数字存储元件**：
   记忆寄存在 RDF 池中，Int–RDF 复合物在脉冲到达时决定翻转方向；
   "一次脉冲一次翻转"可以通过设计实现，但需要满足延迟 > 脉宽。
2. 判定"能"的标准必须公开且面向状态（端点 $Q$ + 中点穿越 + 计数奇偶性）；
   工作窗口不是一组固定常数，而是随输入波形与启动子泄漏共同移动的联合窗口。
3. 当前候选方案（非零泄漏 $\ell_R=0.008$）：$\beta_B=2$、$\beta_R=8$、scale=0.30、
   $k_I=12\ \mathrm{h^{-1}}$，双初态 430 h / 35 稳态周期通过，最差 $Q=0.9966$；
   启动要求为"已知初态 + 低谷接入第一个完整脉冲"。
4. 两个直接实验接口已经给出：**BM3R1/RDF 表达配比（1:4 最大产能）**与
   **Int 标签梯度标尺（原生 ssrA / LAA+4 / SsrA2X / LAA-LAA）**；
   RDF 端启动子推荐 Cello pBM3R1 66 bp，36 bp 版本不建议原样合成。
5. 主要未决项：表达标定、$\eta$（RDF 共丢失）与单细胞噪声；
   在这些完成前，本方案是**明确模型条件下的候选工作点**，不是实验成功率。

---

## 15. 主图清单与图注草稿

| Figure | 类型 | 内容 | 现状 |
|---|---|---|---|
| Figure 3-1 | 机制图（概念） | 状态依赖的开关：PB/LR 两态、生产分配、脉冲到达时的 RDF 池量、两个翻转方向 | 草图已交（`figures/fig03_1_single_bit_mechanism_sketch.png`） |
| Figure 3-2 | 设计图（概念） | 窗口为何有边：tag 窗口三段（多翻 / 通过 / 卡滞）；联合表达窗口与各侧失效机制 | 草图已交（`figures/fig03_2_design_window_sketch.png`） |
| Figure 3-3 | 数据图 | 时钟输入通量 + 稳态轨迹（Int 脉冲 / $S$ / RDF / BM3R1） | 候选图已有 |
| Figure 3-4 | 数据图 | 工作窗口三联：(A) RBS×tag；(B) $(\beta_B,\beta_R)$ 可行域；(C) 泄漏×$\beta_R$ | 候选图已有 |
| Figure 3-5 | 数据图 | 标签机制（间隙/峰值 Int vs $k_{tag}$）+ 异质性 MC | 候选图已有 |
| Figure 3-6 | 数据图 | 启动与扰动通过数 + 长程对比 + input–clearance 配对 | 候选图已有 |
| Figure 3-7 | 表格图 | 参数表：模型条件、主要方案、设计目标、工作窗口 | 已交（`figures/fig03_7_parameter_table.png`） |

> One-pulse-one-flip 时序示意图作为配套概念图保留（待细化）；另附一张生成式细胞情境概念图
> （`figures/optional_concept_cell_context.png`，供美工参考）。

### 图注草稿

**Figure 3-1. 状态依赖的开关（概念草图）。**
DNA 元件有 PB / LR 两种构型。PB 态下 BM3R1 表达并阻遏 RDF 产生，脉冲到达时 RDF 池近乎空
（≈0.002 µM），因此时钟脉冲驱动正向反应（PB→LR）；LR 态下 BM3R1 消退、RDF 累积
（脉冲到达时 ≈4.5 µM，约 2500× 对比），下一个脉冲形成 Int–RDF 复合物驱动反向反应（LR→PB）。
记忆寄存在 RDF 池；延迟必须盖过脉冲宽度。

**Figure 3-2. 设计窗口与失效模式（概念草图）。**
（A）Int 降解标签窗口的两侧：太弱 → 间隙 Int 残留 → 多翻；太强 → 脉冲峰被削 → 回翻失败
（卡在 LR）；确定性通过带 [3, 24] h⁻¹、随机全稳带 [6, 18] h⁻¹。
（B）联合表达窗口是一个平台（非零泄漏方案实测：β_R 4–8 × β_B 0.75–6）；各侧失效机制不同，
泄漏升高使窗口向低 RDF 强度移动。示意版；实测通过/失败图见可行域数据图。

**Figure 3-3. 时钟输入与稳态计数。**
（A）时钟模块交付的 ϕC31 产生通量（Week 4 冻结输入，RBS scale 1）；阴影带为启动参考的低谷。
（B）非零泄漏主要方案（$\beta_B=2$、$\beta_R=8$、scale 0.30）的稳态轨迹：
每个 Int 脉冲（上）恰好触发一次 LR 翻转（中），BM3R1 与 RDF 池（下）决定翻转方向。
模型输出；无实验重复。

**Figure 3-4. 工作窗口。**
（A）Week 4 零泄漏点的 RBS × tag 联合工作区；圆圈为双初态严格通过；全稳平台为 tag 7–18 h⁻¹。
（B）0.8% RDF 泄漏下的 $(\beta_B,\beta_R)$ 可行域（430 h、固定低谷启动）；
红圈为严格通过；宽平台为 $\beta_R$ 4–8 × $\beta_B$ 0.75–6。
（C）泄漏 × $\beta_R$ 容差：泄漏越高，要求的 $\beta_R$ 越低。颜色 = 双初态最差 $Q$。

**Figure 3-5. 降解标签窗口与异质性。**
（A）间隙与峰值游离 Int 随 tag 速率的变化（对数轴）；灰点线为 ~0.1 µM 翻转阈值；
绿带为 RBS 0.45 工作窗口。
（B）异质性蒙特卡洛（每档 30 实例）：tag 6–18 全稳，5 与 22 边缘，3 与 25–26 明显退化。
这是模型参数扰动，不是实验噪声分布。

**Figure 3-6. 启动、长程与输入–清除配对。**
（A）启动相位与参数扰动的严格判据通过数（固定低谷 vs 随机相位）。
（B）长程对比：Week 4 点在非零泄漏下失效，RDF 返回方案与重新调参的原处理均可恢复计数。
（C）用于把实测 Int 标签速率匹配到输入 scale 的配对规则。

**Figure 3-7. 参数表。** 模型条件、非零泄漏主要方案、可测形式的设计目标与工作窗口；
每行标注证据等级。

> 候选图见 `figures/`（`candidate_*.png`），来源与文件名对照见 `figures/README.md`。
> 图号与 D 已交付的第 3 章草稿存在重叠，合并时需统一编号（见交付 README）。

---

## 16. 参考文献

1. Zhao, J., Pokhilko, A., Ebenhöh, O., Rosser, S. J., & Colloms, S. D. (2019). A single-input
   binary counting module based on serine integrase site-specific recombination.
   *Nucleic Acids Research, 47*(9), 4896–4909. DOI: 10.1093/nar/gkz245.
2. Pokhilko, A., Zhao, J., Ebenhöh, O., Colloms, S. D., et al. (2016). The mechanism of ϕC31
   integrase directionality: experimental analysis and computational modelling.
   *Nucleic Acids Research, 44*(15), 7360–7372. DOI: 10.1093/nar/gkw616.
3. Potvin-Trottier, L., Lord, N. D., Vinnicombe, G., & Paulsson, J. (2016). Synchronous
   long-term oscillations in a synthetic gene circuit. *Nature, 538*, 514–517.
   DOI: 10.1038/nature19841.
4. Stanton, B. C., Nielsen, A. A. K., Tamsir, A., Clancy, K., Peterson, T., & Voigt, C. A.
   (2014). Genomic mining of prokaryotic repressors for orthogonal logic gates.
   *Nature Chemical Biology, 10*, 99–105. DOI: 10.1038/nchembio.1411.
5. Nielsen, A. A. K., Der, B. S., Shin, J., Vaidyanathan, P., Paralanov, V., Strychalski, E. A.,
   Ross, D., Densmore, D., & Voigt, C. A. (2016). Genetic circuit design automation.
   *Science, 352*, aac7341. DOI: 10.1126/science.aac7341.
6. Liang, Q., & Fulco, A. J. (1995). Transcriptional regulation of the genes encoding
   cytochromes P450BM-1 and P450BM-3 in *Bacillus megaterium* by the binding of Bm3R1 repressor
   to Barbie box elements and operator sites. *Journal of Biological Chemistry, 270*(31),
   18606–18614. DOI: 10.1074/jbc.270.31.18606.
7. Andersen, J. B., Sternberg, C., Poulsen, L. K., Bjorn, S. P., Givskov, M., & Molin, S.
   (1998). New unstable variants of green fluorescent protein for studies of transient gene
   expression in bacteria. *Applied and Environmental Microbiology, 64*(6), 2240–2246.
   DOI: 10.1128/AEM.64.6.2240-2246.1998.
8. Lies, M., & Maurizi, M. R. (2008). Turnover of endogenous SsrA-tagged proteins mediated by
   ATP-dependent proteases in *Escherichia coli*. *Journal of Biological Chemistry, 283*(32),
   22918–22929. DOI: 10.1074/jbc.M801692200.
9. Szydlo, K., Ignatova, Z., & Gorochowski, T. E. (2022). Improving the robustness of
   engineered bacteria to nutrient stress using programmed proteolysis.
   *ACS Synthetic Biology, 11*(3), 1049–1059. DOI: 10.1021/acssynbio.1c00490.
10. Jadhav, P., Roy, S., Butzin, X. Y., & Butzin, N. C. (2025). Engineering a new SsrA-based
    degradation tag (LAA-LAA) and a bacterial synthetic oscillator.
    *ACS Synthetic Biology, 14*(4), 1062–1071. DOI: 10.1021/acssynbio.4c00612.
11. Klimecka, M. M., Antosiewicz, A., Izert, M. A., et al. (2021). A uniform benchmark for
    testing SsrA-derived degrons in *Escherichia coli*. *Molecules, 26*(19), 5936.
    DOI: 10.3390/molecules26195936.

> 引用格式：完整作者列表、年份、题名、期刊、卷(期)、页码、DOI。
> 本表为文献值或文献区间；模型有效参数与待标定量在正文中单独标注。
