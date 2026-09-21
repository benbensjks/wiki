# scillator–ϕC31 Model wiki作图
## 1. 使用说明

正文使用中文说明构图。每个需要放入图中的标签均给出英文版本，英文栏可直接复制到图片中。

本任务单包含 6 项，均对应已经完成的建模内容，还有部分建模没有优化完成，不放入此次作图：

1. Figure 1：项目总览框架图 / Project overview；
2. Figure 2A：核心回路机制动图 / Core circuit animation；
3. Figure 3A：sponge 阈值机制图 / Sponge threshold mechanism；
4. Figure 4A：敏感性到工程变量流程图 / Sensitivity-to-design funnel；
5. Figure 4B：周期旋钮机制图 / Biological period-knob mechanism；

---

## 2. 统一视觉规范

我截了一张主页的图让ai帮我提取了一下主要的颜色以及让它帮我设计了一下颜色，仅供参考

![屏幕截图 2026-09-20 095655](assets/屏幕截图%202026-09-20%20095655-20260920120409-15qd19s.png)

### 主色

以下色值由主页截图近似取样，正式作图前可用网站 CSS 中的最终变量替换：

| 主色                            | Hex | 主要用途                                 |
| --------------------------------- | ----- | ------------------------------------------ |
| TEMPO 深蓝 / TEMPO deep blue    | `#304B53`    | 振荡器主体、结构框、主文字、稳定状态     |
| TEMPO 橙 / TEMPO orange         | `#D29144`    | 时间变化、周期旋钮、阈值、风险与强调     |
| Dry Lab 绿 / Dry Lab teal-green | `#4F9194`    | ϕC31 输出、模型结果、已验证或已选择路径 |
| 白色 / White                    | `#FFFFFF`    | 主背景                                   |
| 浅灰 / Light gray               | `#F3F5F6`    | 次级卡片和补充信息背景                   |
| 深灰 / Charcoal                 | `#33383A`    | 正文与坐标文字                           |

### 浅色辅助色

| 辅助色                  | Hex | 用途                              |
| ------------------------- | ----- | ----------------------------------- |
| 浅蓝 / Light blue       | `#DCE5E7`    | 振荡器模块浅色填充                |
| 浅橙 / Light orange     | `#F2E2CF`    | 周期档位、阈值区和探索状态背景    |
| 浅绿 / Light teal-green | `#DCEBEC`    | ϕC31、Dry Lab 结果和通过状态背景 |
| 蓝灰 / Blue-gray        | `#8CA0A5`    | 待测参数、边界和次级线条          |

### 元件配色

| 中文含义     | 图中英文 | 颜色建议                    |
| -------------- | ---------- | ----------------------------- |
| TetR         | `TetR`         | TEMPO 深蓝 `#304B53`                 |
| CI           | `CI`         | TEMPO 橙 `#D29144`                   |
| LacI         | `LacI`         | Dry Lab 绿 `#4F9194`                 |
| TetO sponge  | `TetO sponge`         | 浅橙填充 `#F2E2CF` + 橙色位点 `#D29144`       |
| PLtetO1      | `PLtetO1`         | 深蓝轮廓 `#304B53` + 橙色 ON 指示 `#D29144`   |
| ϕC31 mRNA   | `ϕC31 mRNA`         | 浅绿 `#DCEBEC`                       |
| ϕC31 蛋白   | `ϕC31 protein`         | Dry Lab 绿 `#4F9194`                 |
| 下游模块     | `Downstream module`         | 深蓝 `#304B53`，用浅蓝填充区分振荡器 |
| 周期控制     | `Timing control`         | TEMPO 橙 `#D29144`                   |
| 已选择或通过 | `Selected / passed`         | Dry Lab 绿 `#4F9194` + 勾号          |
| 风险或边界   | `Design boundary`         | 深橙 `#B86F35` + 三角图标            |
| 探索性设置   | `Exploratory setting`         | 橙色虚线框 `#D29144`                 |
| 待实验标定   | `Pending calibration`         | 蓝灰空心框 `#8CA0A5`                 |

### 线条

| 中文含义 | 可直接使用的英文标签 | 画法               |
| ---------- | ---------------------- | -------------------- |
| 抑制     | `Repression`                     | 实线 + T 形端点    |
| 激活     | `Activation`                     | 实线箭头           |
| 可逆结合 | `Reversible binding`                     | 双向箭头           |
| 模块边界 | `Model scope boundary`                     | 蓝灰虚线框         |
| 阈值     | `Derepression threshold`                     | TEMPO 橙色水平虚线 |
| 后续方案 | `Future design`                     | 深蓝或蓝灰虚线     |

---

## 3. Figure 1 — 项目总览图 / Graphical abstract

### 构图 / Layout

从左到右排列：

```text
给药时间需求 / Timed-delivery schedule
            ↓
可选择的周期档位 / Selectable period setting
            ↓
可调细胞内时钟 / Tunable intracellular clock
            ↓
ϕC31 脉冲序列 / ϕC31 pulse train
            ↓
下游读取模块 / Downstream pulse reader
```

中央时钟内部画：

- 三基因抑制环 / `TetR–CI–LacI repressilator`；
- TetO sponge / `TetO sponge`；
- PLtetO1–ϕC31 支路 / `PLtetO1-driven ϕC31 output`。

细胞下方放三个并列标签：

- 稳定计时 / `Stable timing`；
- 周期可调 / `Tunable period`；
- 输出可读 / `Readable output`。

ϕC31 与下游模块之间画虚线边界：

- 模型范围边界 / `Model scope boundary`。

### 可直接粘贴的英文

```text
Timed-delivery schedule
Selectable period setting
Tunable intracellular clock
TetR–CI–LacI repressilator
TetO sponge
PLtetO1-driven ϕC31 output
ϕC31 pulse train
Downstream pulse reader
Stable timing
Tunable period
Readable output
Model scope boundary
```

---

## 4. Figure 2A — 核心机制动图 / Core mechanism animation

### 固定元件 / Fixed elements

- 大肠杆菌细胞 / `E. coli cell`；
- 三基因抑制环 / `TetR–CI–LacI repressilator`；
- PLlacO1 驱动 TetR / `PLlacO1-driven TetR expression`；
- PLtetO1 驱动 CI 与 ϕC31 / `PLtetO1-driven CI and ϕC31 expression`；
- PR 驱动 LacI / `PR-driven LacI expression`；
- TetO sponge 位点 / `TetO decoy sites`；
- ϕC31 转录和翻译支路 / `ϕC31 transcription and translation`；
- 四条同步波形 / `TetR`, `CI`, `LacI`, `ϕC31`。

动画主线采用文献中的弛豫振荡解释：每个抑制蛋白经历近似恒定产生的积累阶段，随后进入以生长稀释为主的衰减阶段；一个周期可近似理解为三个连续衰减阶段。

### Frame 1 — TetR 积累并抑制 PLtetO1 / TetR accumulation and PLtetO1 repression

- LacI 较低，PLlacO1 活性较高 / `Low LacI permits PLlacO1 activity`；
- TetR 以近似稳定速率产生并积累 / `TetR accumulates during the production phase`；
- 游离 TetR 抑制 PLtetO1，CI 与 ϕC31 的新生转录较低 / `Free TetR represses PLtetO1, reducing new CI and ϕC31 transcription`；
- 已有 CI 继续下降，PR 逐渐解除抑制 / `Existing CI declines and PR is progressively derepressed`。

### Frame 2 — LacI 上升并启动 TetR 衰减相 / LacI rise initiates the TetR decay phase

- CI 降低后 PR 活性恢复，LacI 开始上升 / `PR activity recovers as CI declines, and LacI begins to rise`；
- LacI 抑制 PLlacO1 / `LacI represses PLlacO1`；
- TetR 新生合成关闭，进入以生长稀释为主的衰减阶段 / `TetR production shuts off and a dilution-dominated decay phase begins`；
- 若模型包含基础降解，可与生长稀释共同标注。

### Frame 3 — Sponge 滴定与有效阈值提高 / Sponge titration raises the effective threshold

- TetR 与高拷贝 TetO decoy sites 快速、可逆结合 / `TetR binds rapidly and reversibly to high-copy TetO decoy sites`；
- TetR 在游离池和 DNA 结合池之间分配 / `TetR partitions between free and DNA-bound pools`；
- TetR 总量 / `Total TetR`；
- 游离 TetR / `Free TetR`；
- sponge 结合 TetR / `Sponge-bound TetR`。

可放一个短公式标签：

```text
Total TetR = free TetR + sponge-bound TetR
```

- PLtetO1 读取的是游离 TetR，而不是总 TetR / `PLtetO1 responds to free TetR, not total TetR`；
- sponge 使同一游离 TetR 阈值对应更高的总 TetR，即提高有效去抑制阈值 / `Sponge binding makes the same free-TetR threshold correspond to a higher total-TetR level`；
- 这样可避开低分子数衰减末端中噪声最大的最后几步 / `This avoids the noisiest last few low-copy decay steps`。

### Frame 4 — PLtetO1 去抑制并形成 ϕC31 脉冲 / PLtetO1 derepression and ϕC31 pulse formation

- 游离 TetR 降至去抑制阈值以下后，PLtetO1 活性升高 / `As free TetR falls below the derepression threshold, PLtetO1 activity increases.`；
- CI 与 ϕC31 转录同步增强 / `CI and ϕC31 transcription increase together`；
- ϕC31 转录 / `ϕC31 transcription`；
- ϕC31 mRNA 累积 / `ϕC31 mRNA accumulation`；
- ϕC31 翻译 / `ϕC31 translation`；
- ϕC31 蛋白脉冲滞后于启动子活性 / `The ϕC31 protein pulse lags behind promoter activity`。

### Frame 5 — CI 介导复位并结束输出脉冲 / CI-mediated reset and pulse termination

- CI 蛋白上升并抑制 PR，LacI 产生下降 / `CI rises and represses PR, reducing LacI production`；
- LacI 降低后 PLlacO1 恢复，TetR 重新进入产生相 / `As LacI declines, PLlacO1 recovers and TetR re-enters its production phase`；
- TetR 增加时有相当一部分分配到 decoy sites，从而缓冲游离 TetR 的上升 / `Decoy binding buffers the rise in free TetR as total TetR increases`；
- PLtetO1 活性逐渐下降 / `PLtetO1 activity decreases`；
- ϕC31 转录停止，mRNA 与蛋白随后下降 / `ϕC31 transcription stops, followed by mRNA and protein decline`；
- 回到 Frame 1 / `The cycle returns to Frame 1`。

### 可直接粘贴的英文 / Copy-ready labels

```text
E. coli cell
TetR–CI–LacI repressilator
TetO decoy sites
PLlacO1-driven TetR expression
PLtetO1-driven CI and ϕC31 expression
PR-driven LacI expression
TetR accumulation and PLtetO1 repression
Low LacI permits PLlacO1 activity
TetR accumulates during the production phase
Free TetR represses PLtetO1
Existing CI declines and PR is progressively derepressed
LacI rise initiates the TetR decay phase
PR activity recovers as CI declines
LacI represses PLlacO1
TetR production shuts off
Dilution-dominated decay phase
Sponge titration raises the effective threshold
Rapid reversible TetR–decoy binding
TetR partitions between free and DNA-bound pools
Total TetR
Free TetR
Sponge-bound TetR
Total TetR = free TetR + sponge-bound TetR
PLtetO1 responds to free TetR, not total TetR
Higher effective total-TetR threshold
PLtetO1 derepression and ϕC31 pulse formation
Intrinsic free-TetR threshold
CI and ϕC31 transcription increase together
ϕC31 transcription
ϕC31 mRNA accumulation
ϕC31 translation
The ϕC31 protein pulse lags behind promoter activity
CI-mediated reset and pulse termination
CI rises and represses PR
LacI production decreases
TetR re-enters its production phase
Decoy binding buffers the rise in free TetR
PLtetO1 activity decreases
ϕC31 transcription stops
The cycle returns to Frame 1
```

### 动画参数 / Animation settings

- 循环时长 / `Loop duration: 8–12 s`；
- 平滑相位切换 / `Smooth phase transitions`。

---

## 5. Figure 3A — Sponge 滴定与计时精度机制 / Sponge titration and timing-precision mechanism

### Panel 1：无 sponge 的低阈值衰减 / Low-threshold decay without sponge

```text
TetR 总量约等于游离 TetR / Total TetR ≈ free TetR
        ↓
PLtetO1 只在很低 TetR 水平下去抑制 / PLtetO1 derepresses only at very low TetR
        ↓
最后几个随机清除事件主导等待时间 / The last few stochastic loss events dominate timing
        ↓
跨阈值时间分布较宽 / Broad threshold-crossing-time distribution
```

画法：用多条 TetR 衰减轨迹表示它们在高拷贝区接近、进入低拷贝尾部后明显分散；低位虚线标 `Intrinsic free-TetR threshold K`。

### Panel 2：可逆 DNA 滴定 / Reversible DNA titration

```text
High-copy TetO decoy sites + PLtetO1 operator sites
        ↓
TetR 在 free 与 DNA-bound pools 间快速分配
Rapid partitioning between free and DNA-bound pools
        ↓
PLtetO1 senses free TetR, not total TetR
        ↓
同一 free-TetR 阈值对应更高的 total-TetR 有效阈值
The same free-TetR threshold maps to a higher effective total-TetR threshold
```

中央保留简式：

```text
Total TetR = free TetR + DNA-bound TetR
```

### Panel 3：减少衰减末端的计时噪声 / Reduced timing noise in the decay tail

```text
提高有效阈值 / Higher effective threshold
        ↓
避开低分子数的最后几步 / Avoid the last few low-copy decay steps
        ↓
更窄的跨阈值时间分布 / Narrower threshold-crossing-time distribution
        ↓
较低的周期变异和相位漂移 / Lower period variability and phase drift
```

### 项目结果条带 / Project-result strip

在机制图底部标记模拟结果：

```text
Project-specific result: sponge increased the measured three-repressor phase period
9.896 → 12.287 generations; ratio = 1.242 [1.232, 1.251]
```

### 曲线与结构标签

- TetR 总量 / `Total TetR`；
- 游离 TetR / `Free TetR`；
- DNA 结合 TetR / `DNA-bound TetR`；
- TetO decoy sites / `TetO decoy sites`；
- 内在游离 TetR 阈值 / `Intrinsic free-TetR threshold K`；
- 有效总 TetR 阈值 / `Effective total-TetR threshold`；
- 跨阈值时间分布 / `Threshold-crossing-time distribution`。

### 可直接粘贴的英文

```text
Without sponge
With sponge
Total TetR ≈ free TetR
PLtetO1 derepresses only at very low TetR
The last few stochastic loss events dominate timing
Broad threshold-crossing-time distribution
High-copy TetO decoy sites
Rapid reversible TetR–decoy binding
Free TetR
DNA-bound TetR
Total TetR
Total TetR = free TetR + DNA-bound TetR
PLtetO1 senses free TetR, not total TetR
Intrinsic free-TetR threshold K
Effective total-TetR threshold
Higher effective threshold
Avoid the last few low-copy decay steps
Narrower threshold-crossing-time distribution
Lower period variability and phase drift
Project-specific period effect
```

---

## 6. Figure 4A — 敏感性到工程变量 / Sensitivity-to-design funnel

### 流程层级 / Funnel stages

```text
全部模型参数 / All model parameters
        ↓
Morris 初筛 / Morris screening
        ↓
关键参数 / Influential parameters
        ↓
Sobol 交互分析 / Sobol interaction analysis
        ↓
主效应与交互效应 / Main effects and interactions
        ↓
生物学可测性 / Biological measurability
        ↓
构建可行性 / Constructability
        ↓
工程变量分类 / Engineering variable classes
```

### 最终三类

#### 周期旋钮变量 / Knob variables

- 振荡器 RBS / `Oscillator RBS strength`；
- mRNA 寿命 / `mRNA lifetime`。

#### 校准参数 / Calibration variables

- 启动子漏表达 / `Promoter leak`；
- 生长与稀释 / `Growth and dilution`；
- sponge 容量 / `Effective sponge capacity`。

#### 输出参数 / Output variables

- ϕC31 RBS / `ϕC31 RBS strength`；
- ϕC31 降解 / `ϕC31 degradation`；
- ϕC31 活性比例 / `Active ϕC31 fraction`。

### 可直接粘贴的英文

```text
All model parameters
Morris screening
Influential parameters
Sobol interaction analysis
Main effects and interactions
Biological measurability
Constructability
Engineering variable classes
Knob variables
Oscillator RBS strength
mRNA lifetime
Calibration variables
Promoter leak
Growth and dilution
Effective sponge capacity
Output variables
ϕC31 RBS strength
ϕC31 degradation
Active ϕC31 fraction
```

---

## 7. Figure 4B — 周期旋钮机制 / Biological period-knob mechanism

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

### 上部：元件位置 / Component locations

```text
启动子 / Promoter
→ 核糖体结合位点 /                 Ribosome-binding site (RBS)
→ TetR / CI / LacI 编码序列 /     TetR / CI / LacI coding sequence
→ mRNA 稳定性元件 /                mRNA stability element
```

### RBS 分支 / RBS branch

```text
RBS 序列 / RBS sequence
→ 翻译起始效率 / Translation initiation efficiency
→ 每条 mRNA 的蛋白产量 / Protein output per mRNA
→ 抑制蛋白积累 / Repressor accumulation
→ 阈值跨越时间 / Threshold-crossing time
```

### mRNA 寿命分支 / mRNA-lifetime branch

```text
mRNA 稳定性元件 / mRNA stability element
→ 转录本寿命 / Transcript lifetime
→ 蛋白产生窗口 / Protein-production window
→ 相位持续时间 / Phase duration
```

两条分支汇合至：

- 振荡周期 / `Oscillation period`。

### 档位 / Period settings

```text
短周期 / Short period
K1 — K2 — K3 — K4 — K5
长周期 / Long period
```

### 可直接粘贴的英文

```text
Promoter
Ribosome-binding site (RBS)
TetR / CI / LacI coding sequence
mRNA stability element
RBS sequence
Translation initiation efficiency
Protein output per mRNA
Repressor accumulation
Threshold-crossing time
Transcript lifetime
Protein-production window
Phase duration
Oscillation period
Short period
Long period
Selected discrete setting
Local deterministic exploration
Oscillator RBS — timing control
ϕC31 RBS — output conditioning
```