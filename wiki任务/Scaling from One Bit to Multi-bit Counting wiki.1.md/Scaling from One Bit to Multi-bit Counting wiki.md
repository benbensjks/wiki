# Scaling from One Bit to Multi-bit Counting wiki

## 1. 页面主线与建模状态 (Page Main Line & Status)

### 页面主线 (Page Flow)

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

### 建模完成状态与证据表 (Status Table)

|**内容/模块**|**状态**|**当前证据与模型边界**|
| --| ------------| -----------------------------------------------------|
|**单 Bit 拓扑重组与 ODE 动力学**|已完成|包含 Pokhilko 构象修饰模型与 13 个生物物理参数标定|
|**BM3R1 时间延迟网络**|已完成|实现高协同性 Hill 抑制与反转门控|
|**二级级联 (I1-FFL 脉冲进位)**|已完成|已实现 Repressilator/方波驱动与$Int_1$进位脉冲关断|
|**三级级联静态泄漏 (Early Leakage)**|问题已识别|发现 P2/P3 阶段 Bit2 提前翻转及未预期回退|
|**泄漏优化 (方案 A & B)**|方程已验证|$A_1$负自馈 +$Int_0$杂合门控在逻辑上可闭合 P1\~P8|
|**实体蛋白标定与参数敏感性**|待完成|具体$A_i/R_i$蛋白基因实体未确定，高 RDF 细胞毒性待实验组评估|

## 2. 第一部分：Carry 接口要求与候选架构比较 (Carry Interface & Candidate Comparison)

> **目标**：明确多级计数器中进位信号 (Carry Signal) 的逻辑与生物学接口规范，评估并筛选可行的进位架构。

### 2.1 Carry 信号的接口要求

- **生物学与逻辑规范**：

  - **单次触发性 (Single-Pulse Triggering)** ：上级归零时仅能吐出单个有效脉冲，避免下级连续翻转。
  - **脉冲宽/高窗口 (Pulse Width & Amplitude Window)** ：脉冲强度必须超过下级整合酶重组阈值，但持续时间必须小于重组完成时间，防止反复翻转。
  - **低静态底噪 (Low Basal Leakage)** ：非进位阶段（静态等待期）转录底噪必须极低，防止累积触发。

### 2.2 候选架构比较 (Candidate Architecture Comparison)

>  *（提示：本小节暂时留空，等和hby交流后补充）*

## 3. 第二部分：I1-FFL 脉冲进位与二级级联 (Level-2 Carry Cascade)

> **目标**：基于非相干前馈环 (I1-FFL) 构造进位脉冲，展示 Bit 0 到 Bit 1 的单次精准进位触发过程（干实验 ODE 仿真）。

### 3.1 I1-FFL 脉冲生成逻辑

- **进位触发条件**：当 Bit $i$ 完成 $LR \rightarrow PB$ (复位归零) 时，恢复组成型启动子表达激活因子 $A_i$。
- **双分支时序差**：

  - **快速分支 (Direct Activation)** ：$A_i$ 迅速拉高下一级整合酶 $Int_{i+1}$ 的转录。
  - **慢速分支 (Delayed Repression)** ：$A_i$ 延迟激活阻遏因子 $R_i$，随着 $R_i$ 超过阈值后关断 $Int_{i+1}$。
- **脉冲形成**：挤出一个有限宽度的窄脉冲 $Int_{i+1}$，确保下一级仅精准翻转一次。

### 3.2 振荡器驱动与二级计数 ODE 方程

- **驱动源输入**：Repressilator 驱动产生 $Int_0$ 脉冲信号。
- **级联方程表达**：

  - 激活因子 $A_0$ 与阻遏因子 $R_0$ 动态方程。
  - AND 逻辑门控下的 $Int_1$ 脉冲生成方程：$f_A(A_0) \cdot g_R(R_0)$。

### 3.3 二级级联仿真结果与三子图拆解

- **Panel 1 (Bit State)** ：Bit 1 发生 $LR \rightarrow PB$ 触发，Bit 2 接收进位信号完成单次翻转。
- **Panel 2 (Concentration)** ：$A_0$ 快速上升与 $R_0$ 缓升的时间差，表现明显的时间延迟效应。
- **Panel 3 (Int Level)** ：生成的非对称 $Int_1$ 窄脉冲 (上升沿由 $A_0$ 主导，下降沿由 $R_0$ 关断)。

### 3.4 二级模型的稳健性小结

- 验证前馈环脉冲进位在二级结构中的可行性，证实“前级归零 $\rightarrow$ 产生单次脉冲 $\rightarrow$ 驱动下级翻转 $\rightarrow$ 脉冲自关断”链路通畅。

## 4. 第三部分：理想化三级级联与泄漏抑制优化 (Idealized Level-3 Cascade)

> **边界声明**：本章展示的所有 3-Bit 计数仿真与泄漏抑制逻辑，均建立在**理想化 ODE 数学模型 (In-Silico Simulation)**  之上，仅作理论可行性验证。调控因子 $A_i/R_i$ 尚未映射至具体的生物实体蛋白，亦未开展湿实验电路构建。

### 4.1 三级级联拓扑架构与二进制递增逻辑

- **级联链路**：$Bit 0 \xrightarrow{Int_1} Bit 1 \xrightarrow{Int_2} Bit 2$。
- **理想二进制状态映射**：P1 $(0,0,1)=1 \rightarrow \dots \rightarrow$ P7 $(1,1,1)=7 \rightarrow$ P8 $(0,0,0)=0$ 循环。

### 4.2 核心瓶颈：静态等待期的提前泄漏 (Early Leakage)

- **问题现象**：在 P2/P3 等静态等待阶段，上游高稳态表达导致 $Int_2$ 出现底噪积累，引发 Bit 2 提前触发翻转或 Bit 1 非预期回退。

### 4.3 两种泄漏抑制设计 (方程已验证)

- **方案 A：**​**$A_1$** **蛋白负自馈调控 (Negative Autoregulation)**

  - **机制与方程**：将 $A_1$ 启动子改造为负反馈启动子。
  - **效果**：压缩稳态浓度，生成陡峭瞬态脉冲后迅速压低驱动力，削弱静态期泄露。
- **方案 B：**​**$Int_0$** **与进位信号的双因子杂合门控 (Split System / Hybrid Promoter)**

  - **机制与方程**：$Int_2$ 的表达采用 AND 双重逻辑，强制依赖外源时钟脉冲 $Int_0$。
  - **理论生物学对应**：杂合启动子 (Hybrid Promoter) 或拆分整合酶 (Split-Integrase)。

### 4.4 完整的 8 阶段二进制循环计数 (P1\~P8)

- **理想化 ODE 仿真验证**：

  - **图 1 (DNA States)** ​：Bit 0, Bit 1, Bit 2 在理想参数下的陡峭翻转与 P1\~P8 完整逻辑闭环。
  - **图 2 (Factors)** ：负自馈作用下 $A_1$ 的低稳态锁死与 $R_1$ 的高位关断门控。
  - **图 3 (Pulses)** ：$Int_0$ 时钟脉冲与 $Int_1, Int_2$ 进位脉冲的时序精准对齐。

### 4.5 鲁棒性、实现边界与架构选择

- **参数鲁棒性分析**：评估理想模型在参数波动下的计数成功率与敏感度排行。
- **现实实现边界与瓶颈**：

  - **正交蛋白系统可确定**：

    - **低位（Bit 0）：采用** **$\Phi\text{C31}$** **重组酶系统，阻遏蛋白为 BM3R1。**
    - **中位（Bit 1）：采用** **$\text{Bxb1}$** **重组酶系统，阻遏蛋白为正交的 PhlF。**
    - **高位（Bit 2）：采用** **$\text{TP901-1}$** **重组酶系统，阻遏蛋白为正交的 PsaR。**

      > 文献来源于：Weinberg, B. H., Pham, N. H., Caraballo, L. D., Lozanoski, T., Engel, A., Bhatia, S., & Wong, W. W. (2017). Large-scale design of robust genetic circuits with multiple inputs and outputs for mammalian cells. *Nature biotechnology*, *35*(5), 453-462.
      >
  - **部分实体蛋白缺失**：真正的湿实验需要寻找具体且无交叉干扰的 $A_i/R_i$ 转录因子。
  - **宿主负担与毒性**：高浓度 RDF 的毒性以及多级级联带来的代谢负载限制。
- **架构选择建议**：基于干实验模型的理论预测，为后续湿实验元件筛选与电路构建提供设计边界参考。

‍
