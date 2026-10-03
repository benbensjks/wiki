# Engineering——C部分迭代：优化FFL模型的级联逻辑门

## 1.循环一

**Design（设计）**

在 3-Bit 基因逻辑计数器的设计初期，我们试图利用非相干前馈环（I1-FFL）生成的暂态整合酶脉冲，实现多级重组状态的精准递增。然而，理论分析表明，系统在静态等待阶段（如 P2/P3 阶段）极易出现高位（Bit 2）因上游基底表达积累而提前触发翻转的“早期泄露”（Early Leakage）问题。为此，我们针对泄露抑制提出了两种改进设计策略：一是引入 $A_1$ 蛋白负自馈调控（Negative Autoregulation），利用其“加速响应、压缩稳态”的特性将 $A_1$ 稳态浓度压制在低位，从源头削弱静态等待期的驱动力；二是设计双因子杂合门控，如双输入杂合启动子或拆分整合酶 Split-Integrase 系统，构建逻辑 AND 门，使下一级整合酶 $Int_2$ 的表达同时依赖上游进位信号与外源时钟信号 $Int_0$。

**Build（构建）**

我们将 Pokhilko 重组动力学与 I1-FFL 级联网络相结合，建立了描述 DNA 拓扑状态（$PB_i/LR_i$）与转录调控因子动态的非线性常微分方程（ODE）模型。为了验证两种改进策略的效果，我们在方程构建中分别予以实施：

1. **方案 A（负自馈调控）** ：在 $A_1$ 动态方程中加入负自馈抑制项 $\frac{K_{auto1}^{n_{auto1}}}{K_{auto1}^{n_{auto1}}+A_1^{n_{auto1}}}$，成功使 $A_1$ 在 Bit 1 翻转后快速达到稳态并平稳锁定在极低水平，有效降低了高稳态对下游引发泄露的隐患；
2. **方案 B（双因子杂合门控）** ：将 $Int_2$ 的表达控制改写为外源脉冲门控 $gate_{Int0}$ 与上游激活/阻遏因子（$A_1, R_1$）的乘积逻辑 AND 门，强制要求上游进位信号与外源时钟信号同时存在，成功阻断了静态等待期间因单一上游积累导致的非预期表达。

两套方案在数学模型层面的分别单独引入，结果显示方案A未能明显抑制提前泄露问题，而方案B能够完美解决当前的泄露问题。

**Test（测试）**

数值仿真结果显示，在经过参数调优后，系统成功实现了从 P1 $(0,0,1)$ 到 P8 $(0,0,0)$ 的 8 阶段完整二进制递增计数与循环复位，并在特定参数下平稳锁定了 $A_1$ 浓度，缓解了 P2/P3 阶段 Bit 2 的提前泄露。但在进一步的测试中也暴露出了关键缺陷：在参数波动下，部分 Bit 的提前翻转仍未能被彻底避免。模型虽然在数学逻辑上可行，但对参数窗口的依赖度较高。

**Learn（学习）**

从本次迭代中我们认识到，不能拘泥于模型最初设计的单一或固有拓扑构造。面对复杂的生物系统和动态泄露问题，仅依靠原有的简单前馈环微调难以彻底消除系统的不稳定性。我们需要打破原有的架构限制，主动引入其他生物学上合理且高效的逻辑门来重构与优化现有模型。通过拓展拓扑自由度并融合多种逻辑门控机制，从网络结构层面增强抗干扰能力，才能构建出更加稳健、适应性更强的基因逻辑计数器。

‍

# 英文版

### **CYCLE 1**

**Design**

In designing the 3-bit genetic counter, we used I1-FFL-generated integrase pulses for stepwise state transitions. However, theoretical analysis revealed "early leakage" in higher-order bits during static waiting phases  due to basal expression accumulation. To suppress leakage, we proposed two strategies: (1) introducing negative autoregulation to $A_1$ to limit its steady-state concentration, and (2) implementing a dual-input hybrid AND gate for $Int_2$ that requires both the upstream carry signal and the clock signal $Int_0$.

**Build**

We integrated Pokhilko recombination dynamics with the I1-FFL cascade network to construct a non-linear ordinary differential equation model that describes DNA topological states ($PB_i/LR_i$) and transcription factor dynamics. To evaluate the effectiveness of both strategies, we implemented them individually within the model:

1. **Strategy A (Negative Autoregulation):**  By incorporating a negative autoregulatory term $\frac{K_{\text{auto}1}^{n_{\text{auto}1}}}{K_{\text{auto}1}^{n_{\text{auto}1}} + A_1^{n_{\text{auto}1}}}$ into the $A_1$ dynamic equation, $A_1$ rapidly reached steady state after Bit 1 flipped and was maintained at a low concentration .
2. **Strategy B (Dual-Input Hybrid Gating):**  By modeling $Int_2$ expression as a multiplicative AND gate combining external pulse gating $gate_{Int0}$ with upstream activator/repressor dynamics , we enforced strict co-dependency on both the carry signal and the clock signal, successfully blocking unintended expression during static waiting periods caused by single-factor upstream accumulation.

When evaluated separately in the mathematical model, Strategy A failed to significantly suppress early leakage, whereas Strategy B perfectly resolved the current leakage issue.

**Test**

> 这里要插入优化前后的对比图

Numerical simulations showed that after parameter optimization, the system successfully performed an 8-stage complete binary increment from P1 $(0,0,1)$ to P8 $(0,0,0)$ with cyclic resetting. However, further testing revealed a key limitation: under parameter fluctuations, premature flipping in certain bits could not be entirely eliminated. Although the model proved logically viable, it exhibited a strong dependency on a narrow parameter window.

**Learn**

From this iteration, we realized that model design should not be constrained by the initial or default network topology. When addressing complex biological systems and dynamic leakage, relying solely on minor adjustments to simple feed-forward loops is insufficient to resolve system instability. We must break away from original structural limitations and proactively incorporate other biologically plausible and effective logic gates, such as AND, NOR, or hybrid composite gates  to optimize and reconstruct the existing model. By expanding topological flexibility and integrating multi-gate control mechanisms, we can enhance system robustness at the network level, laying the foundation for a more reliable and adaptable genetic logic counter.
