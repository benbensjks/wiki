# scillator–ϕC31 Model wiki
## 1. 页面主线

```text
应用需求
  ↓
振荡器如何产生 ϕC31 脉冲
  ↓
sponge 如何改变周期
  ↓
哪些参数适合用于工程调节
  ↓
如何建立周期旋钮
  ↓
噪声、负载和群体相位带来的使用边界
  ↓
模型指导的构建与实验决策
```

确定性模型、文献随机模型和项目扩展随机模型作为证据来源，写入图注和 Methods 折叠栏。

### 建模完成状态

| 内容                                   | 状态     | 当前证据边界                                               |
| ---------------------------------------- | ---------- | ------------------------------------------------------------ |
| 核心 repressilator–sponge–ϕC31 模型 | 已完成   | 已有确定性轨迹、 ϕC31 波形和公式审计                      |
| Sponge 周期效应                        | 已完成   | 已完成有/无 sponge 的轨迹级随机确认                        |
| Morris/Sobol 与工程变量筛选            | 已完成   | 已识别周期、波形和排序相关参数                             |
| RBS × mRNA lifetime 周期旋钮          | 部分完成 | 已有模型档位；真实 RBS 序列与稳定性元件尚待标定            |
| 文献一致随机模型                       | 已完成   | 已完成 SI §4.3.1 结构与统计复现                           |
| 随机 ϕC31–下游接口优化               | 部分完成 | 已提高成功率，尚未达到冻结认证阈值                         |
| Effective-load 情景分析                | 已完成   | 已分析生产能力和 growth/dilution 情景                      |
| 显式宿主资源机制                       | 未完成   | RNAP、ribosome、energy 和资源复合物动力学尚未建立          |
| 独立细胞群体相位扩散                   | 已完成   | 已有 phase-aligned stochastic population audit             |
| 群体同步控制                           | 未完成   | reset、entrainment、cell–cell coupling 和空间模型尚未建立 |

---

带 `🟨` 的标题对应本轮需要美工组绘制的机制图、流程图或框架图。双语标签和具体任务见

`Oscillator–ϕC31 Model wiki作图`。

---

## 2. Hero — A Tunable Intracellular Clock for Timed Delivery

### 三个设计目标

- **Stable timing — 稳定计时：**     形成持续且可辨识的振荡；
- **Tunable period — 周期调节：**     根据给药间隔选择周期档位；
- **Readable output — 输出可读：**     形成可供下游模块读取的 ϕC31 脉冲。

### 🟨 Figure 1 — Graphical abstract / 项目总览图

```text
Timed-delivery schedule
          ↓ selects
Period setting → Repressilator + TetO sponge → ϕC31 pulse → downstream reader
                         │
                         ├── Stable timing
                         ├── Tunable period
                         └── Readable output
```

用细虚线标出模块边界：本模块生成并评估 ϕC31 时间信号；下游计数器由相应模块负责。

---

## 3. How the Clock Works — 回路如何计时并输出

### 🟨 Figure 2A — Repressilator–sponge–ϕC31 mechanism / 核心机制动图

动画分为五步：

1. LacI 较低时 PLlacO1 驱动 TetR 积累；游离 TetR 抑制 PLtetO1，使 CI 与 ϕC31 的新生转录保持较低；
2. CI 下降后 PR 恢复，LacI 上升并抑制 PLlacO1，TetR 合成关闭并进入以稀释为主的衰减相；
3. TetR 在游离池与 sponge 结合池之间快速、可逆分配；PLtetO1 读取游离 TetR，sponge 使相同游离阈值对应更高的总 TetR 有效阈值；
4. 游离 TetR 跨过阈值后 PLtetO1 去抑制，同时驱动 CI 与 ϕC31；ϕC31 mRNA 和蛋白依次累积并形成滞后脉冲；
5. CI 上升并抑制 PR，LacI 下降后 PLlacO1 恢复；TetR 重新积累并逐渐压低 PLtetO1，结束 ϕC31 输出脉冲。

### Figure 2B — Baseline waveform and metrics / 代表性波形与指标

在同一时间轴上展示 TetR、CI、LacI 和 ϕC31，并标出：

- period；
- ϕC31 peak；
- pulse width / FWHM；
- cycle dose / AUC；
- residual level。

完整方程和参数定义放入 `Mathematical formulation` 折叠栏。

---

## 4. Sponge as a Timing-Precision Element — Sponge 如何改变阈值、精度与周期

### 核心问题

TetO sponge 如何通过游离 TetR 提高有效去抑制阈值、降低计时噪声，并在项目模型中改变周期

### 🟨 Figure 3A — Free-TetR threshold mechanism / Sponge 阈值机制图

采用三联机制图：

```text
Without sponge
very low TetR threshold → noisy last few loss events → broad crossing-time distribution

Reversible titration
total TetR = free TetR + DNA-bound TetR
PLtetO1 senses free TetR → higher effective total-TetR threshold

Timing consequence
avoid low-copy decay tail → narrower crossing-time distribution → lower period CV and phase drift
```

### Figure 3B–C — Stochastic sponge confirmation / 随机确认结果

- Panel B：有/无 sponge 的代表性随机轨迹；
- Panel C：轨迹级 period ratio 与 bootstrap confidence interval。

主指标结果：

- 无 sponge：9.896 generations；
- 有 sponge：12.287 generations；
- period ratio：1.242；
- 95% CI：[1.232, 1.251]。

---

## 5. Find Engineering Handles — 识别可调参数

### 核心问题

哪些参数同时具有周期影响力、实验可测性和构建可行性？

### 🟨 Figure 4A — Sensitivity-to-design funnel / 敏感性到工程变量

```text
Morris screening
      ↓
Sobol interaction analysis
      ↓
biological measurability
      ↓
constructability
      ↓
primary knob axes and calibration parameters
```

参数分为三类：

| 类型                  | 参数                                            | 用途           |
| ----------------------- | ------------------------------------------------- | ---------------- |
| Knob variables        | oscillator RBS、mRNA lifetime                   | 形成周期档位   |
| Calibration variables | promoter leak、growth/dilution、sponge capacity | 校准周期和排序 |
| Output variables      | ϕC31 RBS、ϕC31 degradation、active fraction   | 调整输出可读性 |

主页显示 5–6 个关键参数。完整 Morris 热图和 Sobol `S1/ST` 进入补充内容。

---

## 6. Build the Period Knob — 构建周期旋钮

### 核心问题

如何将 RBS 和 mRNA lifetime 映射为可构建的周期档位？

### 🟨 Figure 4B — Biological period-knob mechanism / 周期旋钮机制图

```text
Oscillator cassette
promoter → RBS → TetR / CI / LacI → mRNA stability element
             ↑                            ↑
  translation efficiency           mRNA lifetime
             └──────────────┬─────────────┘
                            ↓
           protein accumulation and threshold timing
                            ↓
                selectable period settings K1–K5
```

图中区分三个层级：

1. **Selected discrete library：**     已完成模型筛选的离散周期档位，并分别标注各项检查状态；
2. **Local deterministic exploration：**     网页中的 promoter/RBS scale 局部扫描

网页中超出 frozen conditions 的组合标记为 exploratory。

**未完成：**   RBS scale 到真实 RBS 序列、mRNA lifetime 到具体稳定性元件的实验映射。该层保留在正文的下一阶段计划中，暂不交给美工组绘制。

### Figure 4C — Period-library performance / 周期档位性能

- 推荐参数组合；
- realized period；
- 相邻档位间距和单调性；
- 随机周期分布；
- ϕC31 输出与负载检查状态。

### 实验映射

1. 用候选 RBS 库或 RBS Calculator 形成强度梯度；
2. 用 reporter 标定相对翻译强度；
3. 测量候选 mRNA 稳定性元件的实际半衰期；
4. 将实测值回填模型；
5. 重新检查周期、ϕC31 波形和负载边界。

---

## 7. Implementation Boundaries — 实现边界

## 7.1 Molecular noise / 分子噪声

文献一致的 exact SSA 用于核对随机反应、参数映射、周期和 CV。项目扩展模型用于评估 sponge、周期档位和 ϕC31 输出。

### Figure 5A — Stochastic evidence / 随机证据

- literature-reference reproduction；
- project stochastic period distributions；
- cycle CV 和 resolved fraction。

图注明确模型层、轨迹数、随机算法和置信区间方法。

## 7.2 Resource burden / 资源负载

### Figure 5B — Effective-load scenario analysis / 有效负载情景分析

**已完成：**     用 effective load factor 与 growth/dilution response 表示负载情景，并通过配对轨迹比较 period change、ϕC31 waveform change 和 provisional load window。所有变化均相对于同一构建、同一 sponge 条件下的 no-load baseline。

**未完成：**     RNAP、ribosome、energy、mRNA–ribosome complex 和 protease competition 的显式动力学；ϕC31 或下游表达经共享资源反馈到振荡器的完整闭环。

主页当前只放负载剂量—响应数据图，并标注 `phenomenological load scenario`。完整资源反馈机制图暂不进入美工作图任务单。

## 7.3 Population phase drift / 群体相位扩散

### Figure 5C — Population phase mechanism / 群体相位机制图

**已完成范围：**     无细胞间耦合的独立随机细胞、初始相位对齐、内在相位扩散、群体平均幅度和相干度。

图中展示：

1. 多个细胞从接近的初始相位开始；
2. 独立周期噪声使细胞相位逐渐分散；
3. 单细胞仍振荡；
4. 群体平均幅度和相干度下降；

正文报告 no-offset phase diffusion 与 initial-delay 的额外影响。

**未完成**    **：**     reset、entrainment、cell–cell coupling、群落空间结构和信号扩散。这些内容只列入后续工作，不出现在当前机制图中。

## 7.4 ϕC31–downstream compatibility / 跨模块可读性（可放在后续模块之后）

用状态卡呈现：

- 输入：period、peak、width、dose、residual；
- 判据：recovery margin 和 cycle-level success；
- 结果：确定性可行区域与随机成功率；
- 边界：B 模块内部重设计由对应模块负责。

---

## 8. Engineering — Model-guided Decisions

## 8.1 Construct-aware clock model / 可映射到构建的时钟模型

| 环节     | 内容                                                                         |
| ---------- | ------------------------------------------------------------------------------ |
| Design   | 连接 promoter、RBS、mRNA lifetime、protein degradation、sponge 和 ϕC31 输出 |
| Build    | 建立显式确定性核心与文献一致的 exact SSA 基线                                |
| Test     | 公式审计、参数映射、文献周期和 CV 复现                                       |
| Learn    | 明确各参数的定义、量纲、证据等级和适用模型                                   |
| Decision | 文献基线与项目扩展分别管理；优先使用可测量和可构建参数                       |

## 8.2 Period-library design / 周期档位设计

| 环节     | 内容                                                          |
| ---------- | --------------------------------------------------------------- |
| Design   | 用 RBS × mRNA lifetime 建立候选档位                          |
| Build    | GSA、二维扫描、随机排序与 ϕC31 输出检查                      |
| Test     | 检查周期误差、档位间距、排序概率、CV、leak 和恢复时间         |
| Learn    | promoter leak 和下游恢复窗口决定部分档位的可用性              |
| Decision | 引入 leak gate、minimum gap 和 output constraints，更新候选库 |

## 8.3 Host and population context / 宿主与群体环境

| 环节     | 内容                                                               |
| ---------- | -------------------------------------------------------------------- |
| Design   | 评估 ϕC31/downstream 负载、随机脉冲和群体相位                     |
| Build    | 资源负载情景、随机接口评价和多细胞相位模型                         |
| Test     | period shift、ϕC31 waveform、cycle success、phase coherence       |
| Learn    | 可用性同时受到下游恢复、资源状态和群体相位的限制                   |
| Decision | 分开设计 timing control 与 output conditioning，并规划群体重置方案 |

### Figure 6 — Engineering decision framework / 工程决策框架图（内容待定）

```text
Biological design
      ↓
Mechanistic prediction
      ↓
Quantitative test
      ↓
Design decision
```

图旁放三个实例卡：parameter mapping、period library、host/population context。

实验测量和模型回填列为下一阶段计划，暂不画成已经闭合的反馈环。

Engineering 的最终叙事尚未确定，本轮不交付美工组绘图。

### 下一阶段实验闭环

1. 测量 PLtetO1、PLlacO1 和 PR 的 ON/OFF 与 leak；
2. 标定候选 RBS 的相对翻译强度；
3. 测量 mRNA/protein lifetime；
4. 测量 sponge copy number 和有效 TetO 容量；
5. 同步记录 growth、resource reporter、ϕC31 pulse 和单细胞周期；
6. 用实测参数更新候选档位并选择代表构建。

---

## 9. Scope, Limitations and Reproducibility

### In scope

- TetR–CI–LacI 振荡器；
- TetO sponge 与游离 TetR；
- PLtetO1–ϕC31 输出；
- GSA 与周期旋钮；
- 分子噪声、资源负载和群体相位。

### Model boundary

- 下游模块用于评价 ϕC31 输入兼容性；
- 染色体架构作为补充探索；
- 体内药代动力学与完整给药效果属于更高层级模型。

### Parameter evidence labels

1. literature value；
2. literature-informed range；
3. model-effective parameter；
4. design variable；
5. measurement pending。

### Reproducibility

- 提供核心代码、输入、参数来源和运行环境；
- 保存 CSV、作图数据、随机种子和版本号；
- 标注 ODE、exact SSA 和 tau-leap；
- `smoke` 用于流程检查，`audit` 用于初步统计，`confirmation` 用于报告结果。

---

## 10. 主图清单

| Figure     | 类型            | 内容                                              |
| ------------ | ----------------- | --------------------------------------------------- |
| Fig. 1     | 框架图          | 项目总览与三个设计目标                            |
| Fig. 2A    | 机制动图        | repressilator–sponge–ϕC31                      |
| Fig. 2B    | 数据图          | 代表性波形与指标                                  |
| Fig. 3A    | 机制图          | sponge 与 free-TetR 阈值                          |
| Fig. 3B–C | 数据图          | 有/无 sponge 轨迹与效应量                         |
| Fig. 4A    | 流程图          | GSA 到工程变量                                    |
| Fig. 4B    | 机制图          | RBS × mRNA lifetime 周期旋钮                     |
| Fig. 4C    | 数据图          | 周期档位性能                                      |
| Fig. 5A    | 数据图          | 随机模型证据                                      |
| Fig. 5B    | 数据图          | effective-load 情景与响应边界；显式资源机制待完成 |
| Fig. 5C    | 机制图 + 数据图 | 已完成范围内的独立细胞群体相位扩散                |
| Fig. 6     | 待定            | Engineering 叙事确定后再决定是否绘制              |