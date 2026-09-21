# Modeling Overview — 第 3 / 4 章接口摘要
## Interface summary for the Overview chapter

> **本文件不接管 Overview 整章**，只提供本模块（单比特计数器 + 多比特进位）可以并入 Overview 的接口文字。
> **文档状态：接口摘要，供 Overview 负责人取用**；本模块只提供文字，不修改 Overview 的结构与叙事。
> 供 Overview 负责人直接取用；本模块正文见 `02_第3章_单比特计数器.md` 与 `03_第4章_多比特级联.md`。

---

## 1. 供「From Biological Functions to Modeling Questions」表格使用

Overview 现有表格已经包含 Single-bit counter 与 Multi-bit counter 两行。建议按下表补足第三列，
使"为什么重要"落到我们实际测到的东西上：

| Biological function | Modeling question | Why it matters（本模块的实测依据） |
|---|---|---|
| Single-bit counter | Under what conditions does one pulse produce exactly one PB↔LR transition? | 翻转方向由**脉冲时刻的 RDF 池量**决定，而不是由两个方向的瞬时竞争决定；因此"一个脉冲一次翻转"要求 RDF 池在两次脉冲之间被正确清空或正确保留。读出需要面向状态的判据：同一批轨迹在瞬间阈值判据下会被误判。 |
| Multi-bit counter | How can a completed transition in one bit generate a carry signal for the next bit? | 直接耦合传递的是**电平**而非**边沿**，下一位会每个周期都翻转。进位需要把状态平台整形成一次窄脉冲，而可靠性的瓶颈不是"门是否打开"，而是**门的关态残留**在下一级大 RDF 池背景下积累出的错误方向剂量。 |
| （建议新增一行）Cross-module timing | How much of the timing behaviour is shared between the counter and the downstream output module? | 两个模块共享 **BM3R1 池**：counter 用它门控 RDF 产生，输出模块用它的残余量当计时器。因此两章的时间尺度不能各自独立选定。 |

---

## 2. 供「Connecting Models to Engineering Decisions」使用

Overview 的信息流句（`oscillator dynamics → Integrase input → DNA-state switching → carry propagation → timed downstream expression`）
在本模块这一段的展开可以写成：

> 时钟模块以真实 C31 翻译通量驱动下游，而不是理想方波。单比特模块检验这个输入能否产生可靠的
> PB↔LR 跃迁，以及重组酶、RDF 与调控蛋白能否在下一次脉冲到来前恢复。
> 由此得到的时序约束再成为多比特进位的要求：上一位的反向重组必须被整形为下一位的**一次**重组酶脉冲。
> 最后，稳定 DNA 状态与输出时长被拆成两个问题——状态决定输出**何时开始**，
> 而共享的 BM3R1 池决定输出**持续多久**。

本模块为 Overview 提供的三条可引用结论：

1. **判据会被误用**：同一批参数点，用"瞬间阈值采样"与"面向状态的有限读窗"两套判据评分，
   后者多认出 20 个真实交替点，而旧判据一个都没多认。因此项目统一采用面向状态的判据。
2. **记忆不在 DNA 上**：DNA 状态只是连续比例 $S$，其走向由 RDF 池与复合物决定；
   记忆寄存在 RDF 池里（两种状态间相差约三个数量级）。
3. **加一级不会自动成功**：两级已经过完整因果认证；三级在骨架与投影等价通过之后，
   模 8 验收仍然失败，失败机制已定位为"次阈值泄漏 + 下一级 RDF 存量被提前消耗"。

---

## 3. 供「Modeling as an Engineering Cycle」使用

本模块贡献的 engineering decisions（可并入 Overview 的清单）：

- 用**真实生成的重组酶波形**替换理想化的计数输入；
- 定义 one-pulse-one-flip 的工作窗口与统一读出判据；
- 约束 RDF / BM3R1 / 清除速率可用的表达区间；
- 比较不同的进位机制，并选择**非相干前馈环**做边沿整形；
- 识别级联扩展时的泄漏与时序限制，并给出**门对比度 + 积分泄漏**的验收口径；
- 识别出"计数器决定何时开始、输出模块决定持续多久"这一功能分离，
  并把它落实到两章共享的 BM3R1 接口上。

---

## 4. 供「How to Read the Modeling Section」使用：本模块对四问的回答

Overview 规定每个模型小节回答四个问题。本模块的回答如下。

### 第 3 章 Single-Bit Counter

| 问题 | 回答 |
|---|---|
| What engineering problem were we trying to solve? | 一个重组酶开关能否在**真实**振荡器波形下，每个周期产生恰好一次可靠的 0↔1 跃迁，而不是跟着输入电平走。 |
| How did we represent the biology mathematically? | 每个比特 11 个连续状态（重组酶 / RDF / BM3R1 各自的 mRNA、未成熟、成熟，加显式 Int–RDF 复合物与 DNA 比例 $S$）；没有逻辑位、没有事件式翻转。 |
| What did the simulations reveal? | 记忆寄存在 RDF 池（两种状态相差约 1400 倍）；翻转方向由脉冲时刻的复合物决定；存在有限工作窗口，其中最关键的未标定量是浓度标尺换算；单比特需要一个建立过程，而不是上电即有效。 |
| How did those results affect the design of TEMPO? | 采用面向状态的读出判据；把浓度标尺列为必须标定的量；把"建立时间"写成设计指标；把同一套判据沿用到多比特。 |

### 第 4 章 Multi-Bit Counting

| 问题 | 回答 |
|---|---|
| What engineering problem were we trying to solve? | 如何把上一位的状态变成下一位的**一次**进位脉冲，并判断这套进位在多级扩展下能走多远。 |
| How did we represent the biology mathematically? | 用非相干前馈环（相干激活臂 $A$ + 滞后阻遏臂 $F$ + 时钟门相乘）；进位脉冲按**事件链**定义：$J_{rev}\to g$ 窗口 → 下一位重组酶 → $S$ 翻转。 |
| What did the simulations reveal? | 两级在 300 h 内产生 7 个完整模 4 循环且四种轨道自洽初态全部通过（是状态机，不是受迫振荡）；可靠性的瓶颈是门的**关态残留**；三级在投影等价通过后仍在模 8 验收失败，机制为次阈值泄漏与下一级 RDF 存量的相位错配。 |
| How did those results affect the design of TEMPO? | 采用事件链式 carry 定义并公开验收口径；为进位调控蛋白增加完整表达链以把平台整形成窄脉冲；把三级如实写成"机制与限制"，不宣称成功；把与输出模块共享的 BM3R1 接口登记为跨章约束。 |

---

## 5. 边界声明（Overview 中不要写成我们已完成的事）

| 事项 | 正确写法 |
|---|---|
| 进位架构 | 是 **abstract / prospective carry architecture**，目前只有模型层面的定义与验证，尚无实验实现 |
| 三级计数 | **未通过验收**。可以写"机制已定位"，不能写"三级计数器已实现" |
| 进位调控蛋白的表达层 | 属新增未标定参数（mRNA 半衰期、成熟半衰期均为 measurement pending） |
| 下游接口 | 两章共享 **BM3R1**；其交付包内没有源码，尚不能确认两边使用同一个 BM3R1 子模型 |
| 生长条件 | 上下游来自未对齐的生长条件，存在三条互斥路线且尚无依据选择；不能写成"同一细胞条件" |
| 分子噪声 | 本模块的扰动是**模型扰动**，不等同实验噪声分布 |

---

## 6. 供 Design Explorer 使用（若纳入本模块）

建议暴露的控件与反馈：

| 控件 | 反馈 |
|---|---|
| 浓度标尺换算 | 单比特与两级是否仍能计数；工作窗口图上当前点是否落在认证集内 |
| 进位成熟半衰期 | 门的开态宽度与最小时间裕量；工作窗口图上的位置 |
| 进位 mRNA 半衰期 | 同上 |
| 读出判据切换（瞬间采样 / 有限读窗） | 同一轨迹下两种判据的读数差异（用于说明判据的必要性） |
| 门对比度显示 | $L_{peak}$ 与积分泄漏剂量，区分"主窗口"与"次阈值开口" |

不建议在 Explorer 中暴露三级控件，除非三级通过验收；当前只建议以"机制演示"形式展示其失败模式。

---

## 7. 与 Overview 的接口约定

- Overview 只引用本模块的结论与边界声明，不复制正文公式；
- 引用具体数字时以 `selected_profile` 与两份章节草稿为准，避免引用已被取代的工作点；
- 交叉引用请使用**页面名或锚点**，不要写"见下一章第几页"；
- 资产路径由网页组统一确认，本模块只提供最终文件名。
