# Scaling from One Bit to Multi-bit Counting wiki 大纲

## 1. 页面主线与建模状态 (Page Main Line & Status)

### 1.1 页面主线 (Page Flow)

```
       Carry 的接口要求
              ↓
          候选架构比较
              ↓
 I1-FFL 非相干前馈环脉冲进位 (二级级联)
              ↓
3-Bit 级联计数逻辑与静态等待期泄漏 ( Early Leakage )

两种泄漏抑制设计：负自馈 ( Negative Autoregulation ) 与 杂合门控 ( Split System )
              ↓
 完整的 8 阶段二进制循环计数 ( P1~P8 ) 
              ↓
     鲁棒性、实现边界和架构选择
```

### 1.2 建模状态与分层证据链 (Status & Evidence Stages)

|**阶段 / 模块**|**状态**|**证据级别 (Evidence Level)**|**建模边界与干/湿实验声明**|
| --| ------------| ----------------------------| ------------------------------------------------|
|**前提继承**|已完成|文献 / Single-Bit 页面链接|继承单 Bit 拓扑开关模型与基本生物物理参数|
|**接口与候选比较**|框架建立|7 维理论对比矩阵|定义 Carry 触发要求与各架构 Scaling 瓶颈|
|**二级级联 (阶段一)**|已完成|**Proof of Principle**|基于理想化 ODE 实现单次进位触发与脉冲自关断|
|**三级失效诊断 (阶段二)**|问题已识别|**Failure Diagnosis**|发现静态等待期底噪积累导致 Bit 2 提前翻转|
|**泄漏抑制 (方案 A & B)**|方程已验证|**In-Silico Verification**|在方程层面补齐负自馈与 Split System 双因子门控|
|**P1–P8 计数 (阶段三)**|逻辑闭环|**In-Silico Verification**|理想参数下完整 8 阶段二进制循环计数波形|
|**实体构建与湿实验**|待完成|**Future Work**|$A_i/R_i$仍为抽象因子，实体元件筛选与基因构筑未完成|

## 2. 系统规范与继承前提 (System Definitions & Prerequisites)

### 2.1 继承的接口与前提

本页面建立在 **单 Bit 的拓扑动力学** 基础之上。关于整合酶 $Int$ 与方向性因子 $RDF$ 驱动的 $PB \leftrightarrow LR$ 拓扑重组速率方程、热力学构象修饰参数及 $BM3R1$ 时间延迟网络。

### 2.2 状态与编码规范 (Binary State & Bit Order)

为消除位序混淆，本项目统一规范如下：

- **Bit 编号与高低位**：Bit 0 为最低有效位 (LSB)，Bit 1 为中间位，Bit 2 为最高有效位 (MSB)。
- **二进制向量表达**：定义状态向量为 $(Bit 2, Bit 1, Bit 0)$。

  - 例如状态 $(0,0,1)$ 表示：$Bit 2 = 0$ (PB 态), $Bit 1 = 0$ (PB 态), $Bit 0 = 1$ (LR 态)。
- **8 阶段计数映射**：P1 $(0,0,1) \rightarrow$ P2 $(0,1,0) \rightarrow$ P3 $(0,1,1) \rightarrow$ P4 $(1,0,0) \rightarrow$ P5 $(1,0,1) \rightarrow$ P6 $(1,1,0) \rightarrow$ P7 $(1,1,1) \rightarrow$ P8 $(0,0,0)$。

### 2.3 进位事件定义 (Carry Trigger Event)

- **Carry Trigger Event**：当 Bit $i$ 发生 **$LR \rightarrow PB$**  **(即从 1 复位归零至 0)**  的下跳沿时，触发进位信号产生；而当 Bit $i$ 发生 $PB \rightarrow LR$ ($0 \rightarrow 1$) 时，不触发进位。

## 3. Carry 接口要求与候选架构比较 (Carry Interface & Trade-offs)

### 3.1 Carry 信号的接口要求

- **单次脉冲触发 (Single-Pulse Triggering)** ：必须在触发沿吐出单一窄脉冲，不可产生连续多重脉冲。
- **脉冲窗口对齐 (Pulse Window)** ：脉冲峰值浓度需超过下级重组阈值 $K_{int}$，且脉冲持续时间 $\tau_{pulse}$ 必须严格限制在重组时间与复位时间之间。
- **底噪漏失控制 (Low Basal Leakage)** ：非进位期间（静态等待期）进位蛋白的表达量需低于重组引发阈值。

### 3.2 候选架构 7 维对比矩阵 (Candidate Architectures Comparison)

|**比较维度 (Metrics)**|**I1-FFL (非相干前馈环)**|**Cascade Delay Switch**|**Split-Integrase Gate**|
| --| -------------------------| ------------------------| ---------------------------|
|**Carry Trigger**|复位归零沿激活$A_i$|启动子方向切替直连|级联 AND 逻辑协同触发|
|**Pulse Width (**​**$\tau_{pulse}$**​ **)**|由$A_i/R_i$抑制延迟决定|由阻遏蛋白降解速率决定|取决于双因子重叠时间窗口|
|**False Carry Risk**|中 (受静态底噪积累影响)|高 (受启动子泄漏影响)|极低 (需双因子同时存在)|
|**Delay Accumulation**|低 (局部前馈快速关断)|高 (多级连续表达延迟)|低 (受时钟脉冲严格对齐)|
|**Leakage Tolerance**|需额外泄漏抑制设计|极差|优秀|
|**Construct Complexity**|中等 (需额外$A_i/R_i$单元)|较低|较高 (需拆分蛋白质粒构建)|
|**Scaling Limit**|预计 3\~4 Bit|1\~2 Bit|预计 \>4 Bit|

## 4. 阶段一：二级级联逻辑验证 (Level-2 Carry Cascade: Proof of Principle)

### 4.1 I1-FFL 脉冲生成逻辑

- **双分支时序差**：

  - **快速分支 (Direct Activation)** ：$A_0$ 迅速拉高下一级整合酶 $Int_1$ 的转录。
  - **慢速分支 (Delayed Repression)** ：$A_0$ 延迟激活阻遏因子 $R_0$，随着 $R_0$ 超过阈值关断 $Int_1$。
- **脉冲形成**：产生非对称 $Int_1$ 窄脉冲（上升沿由 $A_0$ 主导，下降沿由 $R_0$ 关断）。

### 4.2 驱动源与二级级联 ODE 模型

- **输入驱动**：Repressilator 驱动产生 $Int_0$ 方波/振荡脉冲。
- **动态方程**：

  $$
  \frac{dA_0}{dt} = \alpha_A \cdot [Bit0_{PB}] - \gamma_A A_0
  $$

$$
\frac{dR_0}{dt} = \alpha_R \frac{A_0^{n_1}}{K_1^{n_1} + A_0^{n_1}} - \gamma_R R_0
$$

$$
Pulse(Int_1) = \alpha_{int} \frac{A_0^{n_2}}{K_2^{n_2} + A_0^{n_2}} \cdot \frac{K_3^{n_3}}{K_3^{n_3} + R_0^{n_3}}
$$

### 4.3 二级级联仿真结果 (Proof of Principle)

- **结果阐述**：在理想参数下，验证 Bit 0 归零成功挤出 $Int_1$ 单次脉冲并驱动 Bit 1精准翻转一次，证明前馈环进位机制逻辑通畅。

## 5. 阶段二：三级级联失效诊断 (Level-3 Cascade: Failure Diagnosis)

### 5.1 3-Bit 理想逻辑递增与拓扑链路

- **级联链路**：$Bit 0 \xrightarrow{Int_1} Bit 1 \xrightarrow{Int_2} Bit 2$。

### 5.2 核心瓶颈：静态等待期的提前泄漏 (Early Leakage)

- **失效诊断现象**：在 P2 $(0,1,0)$ 到 P3 $(0,1,1)$ 的静态等待阶段，上游激活因子 $A_1$ 的稳态表达持续存在，导致 $Int_2$ 出现底噪积累。
- **系统故障模式**：

  1. **False Carry (误进位)** ：$Int_2$ 底噪超过重组阈值，导致 Bit 2 在未收到进位信号时提前触发翻转。
  2. **Back-tracking (状态回退)** ：过量 $Int_1/Int_2$ 底噪引发 Bit 1 非预期的反向重组。

## 6. 两种泄漏抑制设计 (Leakage Mitigation Designs)

### 6.1 方案 A：$A_1$ 蛋白负自馈调控 (Negative Autoregulation)

- **机制与方程修饰**：将激活因子 $A_1$ 的启动子改造为带负反馈关断的启动子：

  $$
  \frac{dA_1}{dt} = \alpha_A \cdot [Bit1_{PB}] \cdot \frac{K_{auto}^{m}}{K_{auto}^{m} + A_1^{m}} - \gamma_A A_1
  $$
- **抑漏效果**：极大地压缩稳态下的 $A_1$ 浓度，使进位过程仅产生陡峭的瞬态脉冲，脉冲结束后驱动力迅速压低至阈值以下。

### 6.2 方案 B：$Int_0$ 与进位信号的双因子杂合门控 (Split System / Hybrid Promoter)

- **机制与方程修饰**：强制使 $Int_2$ 的表达依赖于外源主时钟脉冲 $Int_0$ 与上游进位信号 $A_1$ 的 AND 逻辑：

  $$
  Pulse(Int_2) = f(A_1, R_1) \cdot \frac{Int_0^{p}}{K_{clock}^{p} + Int_0^{p}}
  $$
- **抑漏效果**：即使 $A_1$ 在静态等待期有泄漏，若缺乏 $Int_0$ 主时钟脉冲的协同，下级 $Int_2$ 仍被锁定在关断状态。

## 7. 阶段三：优化后 3-Bit 计数验证 (In-Silico P1–P8 Validation)

>  **边界声明**：本章展示的 8 阶段二进制循环计数均属于 **理想化 ODE 数学模型 (In-Silico)**  验证结果。

### 7.1 完整 8 阶段二进制循环计数 (P1\~P8)

- **仿真波形验证**​：展示方案 A 与方案 B 作用下，3-Bit 逻辑计数器在 P1\~P8 阶段无误进位、无回退的完整波形。

### 7.2 动态脉冲对齐与状态转换时序

- **时序拆解**：展示 $Int_0$ 时钟脉冲、$Int_1$ 及 $Int_2$ 进位窄脉冲在 8 个阶段中的精准时间对齐轨迹。

## 8. 鲁棒性、实现边界和架构选择 (Benchmarking & Boundaries)

### 8.1 模型假设与参数化 (Assumptions & Parameterization)

- **假设条件 (Assumptions)** ：假设细胞内资源（核糖体/RNA 聚合酶）无极化竞争；假设质粒拷贝数均匀分布。
- **参数来源分类 (Parameter Classification)** ：

  1. **文献测量参数**：整合酶结合速率 $k_{int}$、重组速率 $k_{rec}$ 。
  2. **拟合参数**：$BM3R1$ 的 Hill 抑制系数与解离常数 $K_D$。
  3. **纯设计变量 (Design Variables)** ：激活/阻遏因子降解速率 $\gamma_A, \gamma_R$、前馈环转录强度比 $\alpha_A/\alpha_R$。

### 8.2 成功评估指标与鲁棒性测试 (Benchmarking Metrics)

定义计数器运行成功与否的量化评估标准：

- **False Carry Rate (误进位率)** ：静态等待期未发生归零事件时下级 Bit 翻转的概率（目标 $< 1\%$）。
- **Missing Carry Rate (漏进位率)** ：发生归零事件时下级未能在时间窗口内翻转的概率（目标 $< 1\%$）。
- **Carry Delay Window (**​**$\tau_{delay}$**​ **)** ：进位脉冲峰值与上级归零沿之间的时间差。
- **周期通过率 (Pass Rate)** ​：计数器无故障连续完成 P1\~P8 完整 8 阶段循环的参数空间比例。

### 8.3 全参数敏感性排行 (GSA)

- 结合 Sobol 敏感性分析，列出对“周期通过率”影响最大的 Top 5 生物物理参数。

### 8.4 物理实现边界与实体构建展望

- **正交蛋白系统可确定**：

  - **低位（Bit 0）：采用** **$\Phi\text{C31}$** **重组酶系统，阻遏蛋白为 BM3R1。**
  - **中位（Bit 1）：采用** **$\text{Bxb1}$** **重组酶系统，阻遏蛋白为正交的 PhlF。**
  - **高位（Bit 2）：采用** **$\text{TP901-1}$** **重组酶系统，阻遏蛋白为正交的 PsaR。**

    > 文献来源于：Weinberg, B. H., Pham, N. H., Caraballo, L. D., Lozanoski, T., Engel, A., Bhatia, S., & Wong, W. W. (2017). Large-scale design of robust genetic circuits with multiple inputs and outputs for mammalian cells. *Nature biotechnology*, *35*(5), 453-462.
    >
- **干/湿实验脱节澄清**：目前 $A_i/R_i$ 仅为逻辑功能因子，湿实验中尚未完成特定正交蛋白（如 TetR/LacI 衍生物或 SigC/Anti-SigC）的实体映射。
- **架构选择建议**：基于干实验模型预测，方案 B (Split System) 在抑制静态泄漏方面具有更高的参数容忍度，建议湿实验组优先考虑杂合启动子或拆分整合酶构建方案。
