<!-- 状态：本稿与中文初稿、Hybrid_Three_Bit_Cascade_Wiki_EN.md同步，取代wikiformal.md第5节混合模型正文；其他章节不在本次替换范围。中英逐段对照，公式和数值表共用。图位M-25、M-26尚未成图，仅用编辑注释预留。 -->

## 混合三级级联：让低位状态产生一次有效进位
## A Three-Bit Hybrid Cascade: Generating a Single Effective Carry from a Lower-Order State Transition

**中文**

多位计数需要同时完成两件事：每个位保存当前的0或1；低位回零时，向高位传递一次有效写入。为检验两种模型能否完成这一过程，我们把游离重组酶生化动力学模型中的bit0、bit2，与前馈环约化模型中的第一级进位网络和bit1连接，在同一次ODE积分中运行。

**English**

Multi-bit counting requires both state retention within each bit and a single effective write signal to the next bit when the lower-order bit resets. To examine whether the two model formulations can jointly satisfy these requirements, we coupled bit0 and bit2 from the free-integrase biochemical kinetics model to the first carry network and bit1 from the reduced feed-forward-loop model. The connected modules were integrated as a single system of ordinary differential equations (ODEs).

**中文**

韩亚轩上游驱动下，这个混合模型在300 h内获得27个连续可判定读数。丢弃第一个完整模8周期后，剩余19个读数仍逐拍加1，完成三级计数。

**English**

With the upstream oscillator model developed by Han as the driver, the hybrid system produced 27 consecutive, unambiguous readouts over 300 h. After excluding the first complete modulo-8 counting cycle, the remaining 19 readouts continued to increment by one modulo 8 at each clock cycle, demonstrating continuous three-bit counting under the specified conditions.

### 三个存储位，两处进位接口 / Three State-Holding Bits and Two Carry Interfaces

**中文**

主结构为 **HBY bit0 → ZMH bit1 → HBY bit2**，下文简称HZH。bit0与bit2区分mRNA、未成熟蛋白、成熟游离蛋白和Int–RDF复合物；中间位沿用曾同学现行《完整二级级联.py》的前馈环及约化方程。三部分各自保留参数来源。

**English**

The principal architecture is **HBY bit0 → ZMH bit1 → HBY bit2**, hereafter referred to as HZH. The bit0 and bit2 modules distinguish mRNA, immature protein, mature free protein, and integrase–recombination directionality factor (Int–RDF) complexes. The middle module retains the feed-forward network and reduced equations from Zeng’s current two-bit implementation, `完整二级级联.py`. Each module preserves its documented parameter provenance.

| 模块 / Module | 在混合模型中的作用 / Role in the hybrid model | 连续状态数 / Number of continuous states |
|---|---|---:|
| HBY bit0 | 接收上游C31表达输入，保存低位DNA状态 / Receives the upstream C31 expression input and stores the least significant DNA state | 11 |
| ZMH A0/F0与bit1 / ZMH A0/F0 and bit1 | 产生Int1进位脉冲，保存中位DNA状态 / Generates the Int1 carry pulse and stores the intermediate DNA state | 2＋4 |
| HBY A1/F1与bit2 / HBY A1/F1 and bit2 | 产生Int2输入，保存高位DNA状态 / Generates the Int2 input and stores the most significant DNA state | 6＋11 |

**中文**

接收模型共34态。韩上游先单独求解，再单向驱动bit0，其振荡器状态不计入这34态。bit0、bit1、bit2分别是最低位、中位和最高位，数字由三个DNA状态读出。

**English**

The downstream system comprises 34 continuous states. The upstream oscillator is solved separately and supplies a prescribed, unidirectional input to bit0; its states are excluded from this count. Bits 0, 1, and 2 represent the least significant, intermediate, and most significant bits, respectively. The numerical count is decoded from their DNA configurations.

<!-- 图位 M-25 / Figure placeholder M-25: 完整混合机制总览 / Complete hybrid mechanism overview. Caption: Module coupling in the HZH three-bit hybrid model. PB0 and Int1 connect the lower two bits; the PB state of bit1 drives the A1/F1 network, while mature free Int0 from the same bit0 supplies the second-stage clock factor. DNA inversion and RNA transcription are represented separately. -->

### 第一处连接：把bit0回零转换为Int1脉冲 / Interface 1: Converting the Reset of Bit0 into an Int1 Pulse

**中文**

bit0由韩上游的C31转录与翻译驱动。其C31 mRNA只建模一次，经未成熟态形成游离Int0，再参与DNA正反向重组。这里使用实际表达通量，没有再乘方波模型中的输入增益。

**English**

Bit0 is driven by the C31 transcription and translation fluxes computed from Han’s upstream model. The C31 transcript is represented once and supplies the immature Int0 pool, which subsequently matures into free Int0 and participates in DNA recombination. No additional gain from the square-wave input model is applied to this expression pathway.

**中文**

bit0从LR回到PB时，$PB_0=1-S_0$增加，驱动A0产生。A0迅速激活Int1产生，同时建立较慢的F0抑制臂；随着F0积累，Int1产生被关断。这样，持续的DNA电平能够形成一次暂态写入输入。

**English**

When bit0 undergoes the LR-to-PB reset, the increasing PB fraction, $PB_0=1-S_0$, drives A0 production. A0 activates Int1 synthesis through the direct arm and induces the slower F0 inhibitory arm. As F0 accumulates, it suppresses further Int1 synthesis. This incoherent feed-forward structure therefore converts a sustained DNA-state signal into a transient write input.

**中文**

定义激活函数$H(x;K,n)=x^n/(K^n+x^n)$，抑制函数$G=1-H$，第一级门为：

**English**

Using the activation function $H(x;K,n)=x^n/(K^n+x^n)$ and its inhibitory counterpart $G=1-H$, the first-stage gate and Int1 dynamics are:

$$
g_0=H(A_0;1.0,2)\,G(F_0;2.6454,5.4332),
\qquad
\frac{dI_1}{dt}=22.728\,g_0-4.8I_1.
$$

**中文**

该式按小时计，$I_1$是bit1的成熟整合酶。此段沿用直接蛋白表达方程，没有额外增加mRNA或成熟延迟。F0是进位抑制因子，在《完整二级级联.py》中记为r0／R0，在混合代码中记为F0_zmh；它与位内RDF是不同分子池。

**English**

Time is expressed in hours, and $I_1$ denotes the mature integrase in bit1. The production coefficient is 22.728 a.u. h$^{-1}$, and the clearance rate is 4.8 h$^{-1}$. This module retains direct protein-level dynamics without introducing additional mRNA or maturation states. The carry inhibitor F0 is denoted r0/R0 in `完整二级级联.py` and F0_zmh in the hybrid implementation; it is distinct from the RDF pool within each bit.

<!-- 图位 M-26a / Figure placeholder M-26a: PB0, A0/F0, and g0/Int1 on a common time axis, with an inset showing the direct activation arm and delayed inhibitory arm. A0 need not itself be a narrow pulse. -->

### 第二处连接：由bit1状态与同源时钟控制Int2表达 / Interface 2: Controlling Int2 Expression with the Bit1 State and a Shared Clock

**中文**

曾同学bit1的PB比例直接驱动A1表达。A1带负自馈，并激活F1产生；A1激活、F1抑制与bit0的成熟游离Int0共同决定第二级门：

**English**

The PB fraction of the ZMH bit1 module drives A1 expression. A1 is negatively autoregulated and activates F1 production. The second-stage gate combines A1 activation, F1 inhibition, and a clock factor derived from mature free Int0 in bit0:

$$
g_1=H(A_1;1.2,6)\,G(F_1;0.4,4)\,H(I_0;0.3,2),
\qquad u_{I2}=38g_1\;\mathrm{a.u./h}.
$$

**中文**

门输出先换算为bit2的mRNA产生源，再经过$mRNA\rightarrow$未成熟Int2$\rightarrow$成熟游离Int2。游离Int2支持正向重组，与RDF2形成的复合物支持反向重组，bit2据此交替写入和回零。

**English**

The gate specifies a target mature-protein production rate, which is mapped to the mRNA synthesis source in bit2. The expression pathway then proceeds from mRNA through immature Int2 to mature free Int2. Free Int2 supports forward recombination, whereas its complex with RDF2 supports reverse recombination, allowing bit2 to alternate between writing and resetting.

**中文**

时钟读取的是这条混合电路自身的Int0。参考配置中的激活门指数6只作用于Int2产生门，A1驱动F1产生的指数仍为4；6是有效响应陡度，不表示已经设计出六个结合位点。

**English**

The clock factor reads Int0 from the same hybrid circuit. In the reference configuration, the activation exponent of 6 applies only to the Int2 production gate; the exponent governing A1-induced F1 synthesis remains 4. The value 6 represents effective response steepness and does not imply an experimentally implemented promoter with six binding sites.

<!-- 图位 M-26b / Figure placeholder M-26b: A local view of PB1, A1/F1, g1, mature Int2, and S2. The inset should distinguish the inhibitory F1 input from the activating A1 and clock inputs and show the Int2 expression chain. -->

### 300 h内连续完成三级计数 / Continuous Three-Bit Counting over 300 h

**中文**

我们在Int0相邻峰之间的谷值附近读取DNA状态，读窗总宽约为一个时钟周期的20%。LR比例不高于0.30读为0，不低于0.70读为1；至少80%的窗内样本需要落在相应带中，才能确定该位标签。

**English**

DNA states were decoded within windows centered on the Int0 trough between consecutive peaks. Each window spanned approximately 20% of the corresponding clock period. A bit was assigned 0 if at least 80% of the samples had an LR fraction no greater than 0.30, or 1 if at least 80% had an LR fraction no less than 0.70. Windows that satisfied neither condition were left unlabelled.

![300 h counting trajectory of the Han-driven HZH hybrid model](wiki_submission_assets/hybrid_threebit_300h.png)

**中文图注**：同一次积分得到的输入、三位DNA状态、数字读出和整合酶波形。27个读窗均可判定，读数依次为1、2、3、4、5、6、7、0并重复；三个位的窗内最低承诺度均为1.0。结果对应HZH的指定参数与韩上游输入。

**English caption**: Input, DNA-state trajectories, decoded counts, and integrase waveforms from a single simulation of the HZH system. All 27 read windows were unambiguous, yielding repeated sequences of 1, 2, 3, 4, 5, 6, 7, and 0. The minimum within-window commitment across all three bits was 1.0. These results apply to the specified HZH parameter set and Han upstream input.

| 检查项 / Assessment | HZH运行结果 / HZH result |
|---|---|
| 连续仿真 / Simulation duration | 300 h |
| 全部／稳态读窗 / Total / steady-state read windows | 27 / 19；丢弃前8窗后仍逐拍模8递增 / Modulo-8 increments were retained after excluding the first eight windows |
| 未标注／边界裁剪读窗 / Unlabelled / boundary-clipped windows | 0 / 0 |
| bit0／bit1／bit2穿越0.5次数 / Number of crossings of 0.5 by bit0 / bit1 / bit2 | 29 / 14 / 7 |

**中文**

低位最频繁翻转，高位逐级减少约一半。三个DNA状态因此组成连续的模8计数，而每一级进位只需要在低位回零时作用一次。

**English**

The least significant bit switched most frequently, with successive higher-order bits switching at approximately half the frequency of the preceding bit. Together, the three DNA states encoded continuous modulo-8 counting, with one effective carry associated with each lower-order reset.

### 进位事件与完整生化初态的复核 / Verification of Carry Events and Complete Biochemical Initial States

**中文**

除数字读出外，我们还独立检测反向重组、进位门和高位翻转，核对事件关联。两级分别得到：

**English**

In addition to decoding the numerical sequence, we independently detected reverse-recombination episodes, carry-gate events, and higher-order DNA transitions to evaluate their associations:

| 进位级 / Carry stage | 低位反向重组 / Lower-order reverse-recombination events | 门事件 / Gate events | 高位翻转 / Higher-order transitions |
|---|---:|---:|---:|
| bit0 → bit1 | 14 | 14 | 14 |
| bit1 → bit2 | 7 | 7 | 7 |

**中文**

最初两次进位作为启动阶段单列；之后的主要事件满足既定的一一对应、时间筛选与翻转方向交替规则。读窗与相邻DNA翻转之间的最小时序间隔约1.18 h，受限项是bit0的hold。

**English**

The first two carry events at each stage were treated separately as startup events. Subsequent events satisfied the predefined one-to-one association and temporal criteria, with alternating higher-order transition directions. The minimum timing margin between a read window and an adjacent DNA transition was approximately 1.18 h; the limiting quantity was the hold margin of bit0.

**中文**

收紧积分容差，并把输出间隔从1 min缩短到0.5 min后，码串和两级事件结论保持一致。我们还从连续八个稳态读窗的谷值抽取完整34态初值，跨一个完整模8周期，依次覆盖数字1–7、0。各组延续自己的原上游时刻再运行300 h，均从预期下一数字开始，并通过相同判据。对应的一手记录见[八数字初态验证摘要](hby_zmh_hby/certification/eight_phase_results/20260928_144314_590882/summary.json)。

**English**

The decoded sequences and event-association verdicts were unchanged when integration tolerances were tightened and the output sampling interval was reduced from 1 min to 0.5 min. Eight complete 34-state initial conditions were extracted at the troughs of eight consecutive steady-state read windows, spanning one complete modulo-8 cycle and representing counts 1–7 followed by 0. Each was continued for an additional 300 h while preserving its corresponding absolute upstream time. All eight simulations began with the expected successor count and passed the same evaluation criteria, as documented in the [primary eight-digital-state verification record](hby_zmh_hby/certification/eight_phase_results/20260928_144314_590882/summary.json).

**中文**

这些初态来自同一条确定性轨道，说明计数相位可以被保存。它们不代表任意初态或噪声条件下的成功率。详细规则与复算记录见[混合模型验收报告](hby_zmh_hby/certification/认证结果说明.md)。

**English**

These initial conditions are phase-consistent states on a single deterministic trajectory and support retention of the counting phase. They do not quantify performance under arbitrary initial conditions or stochastic fluctuations. The complete criteria and numerical checks are documented in the [hybrid-model verification report](hby_zmh_hby/certification/认证结果说明.md).

<!-- 八数字初态图建议采用紧凑矩阵 / Suggested initial-state figure: a compact matrix of initial count, first new count, and evaluation outcome. Avoid repeating eight full-length trajectories. -->

### 时钟门改变了写入波形和落态深度 / The Clock Gate Modifies the Write Waveform and State-Settling Depth

**中文**

为检验Int0时钟因子的作用，我们保留其余方程，将该因子置为1。标称点的模8读出及两级事件关联仍通过，但g1增强，bit2驻留状态更接近读出带边。时钟门在这个工作点不是维持数字计数的必要条件，却会改变写入过程。

**English**

To assess the contribution of the Int0 clock factor, we set this factor to 1 while retaining the remaining equations. Modulo-8 decoding and the event-association criteria at both stages remained satisfied at the nominal operating point. However, the carry output g1 increased, and the bit2 plateau states moved closer to the read-band boundaries. Thus, the clock gate was not required for nominal digital counting in this test, although it altered the write dynamics.

![HZH comparison with the Int0 clock factor retained or bypassed](M21_时钟门消融/out/fig_M21_clock_gate_bypass.png)

**中文图注**：两组数字标签均保持模8，连续DNA状态并不相同。移除时钟因子后，bit2最小读窗带边距离由0.2821降至0.1857；这一距离描述落态深度，不是实测抗噪声能力。

**English caption**: Both conditions retained modulo-8 digital labels, but their continuous DNA-state trajectories differed. Bypassing the clock factor reduced the minimum bit2 read-window distance to a band boundary from 0.2821 to 0.1857. This metric describes state-settling depth and does not constitute an experimental measurement of noise tolerance.

**中文**

进一步把无时钟组在一个完整模8周期内的Int2源总积分配平，带边距离恢复到0.2524，仍与保留时钟组不同：

**English**

We subsequently rescaled the Int2 synthesis source in the clock-bypassed condition to match its time integral to that of the clock-retained condition over one complete modulo-8 cycle. The minimum band-edge distance increased to 0.2524 but remained below the clock-retained value:

| 配置 / Condition | 模8读出及事件判据 / Modulo-8 decoding and event criteria | bit2最小读窗带边距离 / Minimum bit2 read-window band-edge distance |
|---|---|---:|
| 保留时钟门 / Clock retained | 通过 / Passed | 0.2821 |
| 时钟因子置1 / Clock bypassed | 通过 / Passed | 0.1857 |
| 时钟因子置1，并配平Int2源总积分 / Clock bypassed, Int2 source integral matched | 通过 / Passed | 0.2524 |

**中文**

该对照说明总输入量能够解释部分变化，输入的时间分布仍需单独考察。相同的源积分不保证成熟Int2波形和DNA驻留状态相同，不能据此把各项差额拆成独立、可相加的生化作用。

**English**

This comparison is consistent with a contribution from total input magnitude, while also indicating that the temporal distribution of the input requires separate consideration. Equal source integrals do not ensure identical mature Int2 waveforms or DNA plateau states. Differences between these conditions cannot be interpreted as independent, additive biochemical causal effects.

<!-- M-22完整三臂图 / Full three-condition figure: M21_时钟门消融/out/fig_M21_three_arm.png. Use as supplementary material; reduce long in-image titles before publication. -->

### 工作区取决于接口和接收模块 / The Operating Region Depends on the Interface and Receiving Module

**中文**

混合模型的参数扫描显示，调高第二级A1门指数能够增加bit2的时序余量，并提高对A1激活阈值偏移的容忍度；在已测网格内，它没有移动浓度换算系数（uM_per_au）的计数边界。上侧浓度边界首先表现为bit1无法稳定落入读出带。由于后级门不反馈至bit1，改变它不能修复这个低位状态。

**English**

Parameter scans of the hybrid model showed that increasing the second-stage A1 activation exponent improved bit2 timing margins and extended tolerance to shifts in the A1 activation threshold. Within the tested grid, it did not shift the counting boundary along the concentration-conversion factor (uM_per_au) axis. Failure near the upper boundary of this factor first manifested as inadequate read-band commitment in bit1. Because the downstream gate does not feed back to bit1, changing that gate cannot restore the lower-order state.

**中文**

时钟阈值也需要与实际Int0波形配合。同一批配对扫描包含10个K值、2个时钟指数和3个末级配置，共60次运行。其中HZH的20个点在时钟指数2和3下、所测10个K格点（0.075–0.60）均通过计数及事件判据；其余40个点属于HZZ负自馈开、关两臂。这是离散采样结果，格点之间的连续可行区和界外行为仍未证明。

**English**

The clock threshold must also be evaluated relative to the Int0 waveform supplied by the coupled model. The paired scan comprised 60 runs: 10 K values, two clock exponents, and three receiving-stage configurations, all retaining the same lower two bits. All 20 HZH runs passed both counting and event-association criteria at the 10 tested K values spanning 0.075–0.60 and clock exponents of 2 and 3. The remaining 40 runs evaluated HZZ with A1 autoregulation enabled or disabled. These results establish performance at the sampled points; they do not demonstrate a continuous feasible interval or characterize behavior outside the tested range.

<!-- 图位 M-18 / Figure placeholder M-18: plot discrete operating points against K and K/Int0_peak. Distinguish counting failure from event-criterion failure explicitly, without implying continuity between sampled points. -->

### 更换末级：HBY–ZMH–ZMH的补充连接 / An Alternative Receiving Stage: The HBY–ZMH–ZMH Cascade

**中文**

我们保留相同的前17态，把末级改为曾同学《前馈三级级联.py》的探索性A1/F1与bit2方程，得到23态HZZ模型。该末级使用另一套参数，原时钟门为$H(I_0;0.4,3)$，不能视作现行二级参数的简单换算。

**English**

As a supplementary coupling test, we retained the same first 17 states and replaced the receiving stage with the exploratory A1/F1 and bit2 equations from Zeng’s three-bit script, `前馈三级级联.py`. The resulting HZZ system comprises 23 states. Its receiving stage uses a distinct parameter set and originally employs the clock factor $H(I_0;0.4,3)$; these parameters are not a simple unit conversion of the current two-bit implementation.

**中文**

保留末级方程与参数、关闭负自馈，并设置与方程一致的初值后，直接接入HBY Int0未形成模8。将时钟阈值重新设为0.10后，600 h运行获得48个稳态模8读数；从数字0–7抽取完整23态并保持各自上游相位续算300 h，八组均通过计数及事件判据。

**English**

With A1 negative autoregulation disabled, the receiving-stage equations and parameters retained, and initial conditions consistent with the active equations, direct coupling to HBY Int0 did not produce modulo-8 counting. After the clock threshold was reset to 0.10, a 600 h simulation yielded 48 steady-state modulo-8 readouts. Eight complete 23-state initial conditions representing counts 0–7 were then continued for 300 h each, preserving their corresponding upstream phases. All eight continuations passed the counting and event-association criteria.

**中文**

这一补充说明，另一种末级也可完成连接，但需要重新检查时钟与进位输入的匹配。它不证明原探索性三级参数直接成功，也不能把两种末级的K数值当成相同的物理阈值。[HZZ长时程结果](hby_zmh_zmh/failure_attribution/results/round3_20261001_135824/round3_summary.json)和[修正后的八数字初态结果](hby_zmh_zmh/failure_attribution/results/round5b_20261001_173719/eight_states.json)分别记录了两项验证。

**English**

This supplementary test demonstrates that an alternative receiving stage can support the cascade following reassessment of the clock–carry interface. It does not establish successful direct transfer of the original exploratory three-bit parameter set, nor does it justify treating the K values of the two receiving stages as identical physical thresholds. The [HZZ long-duration results](hby_zmh_zmh/failure_attribution/results/round3_20261001_135824/round3_summary.json) and [corrected eight-state continuations](hby_zmh_zmh/failure_attribution/results/round5b_20261001_173719/eight_states.json) document these evaluations separately.

<!-- HZZ可作为页面末尾补充 / HZZ may be placed in a supplementary section. If a figure is included, use the corrected comparison or eight-digital-state figure and identify the exploratory receiving-stage parameters, K=0.10, and disabled A1 autoregulation. -->

### 对实际构建的意义与当前边界 / Implications for Implementation and Current Limitations

**中文**

混合连接给实验表征提出了具体任务：先确认低位能够稳定落态，再测回零后产生的进位脉冲是否能完成高位写入，并在下一次读取前关断。需要共同比较脉冲宽度、实际成熟整合酶、RDF时序和正反向写入通量，单独的峰值不能替代完整周期。

**English**

The hybrid coupling identifies specific requirements for experimental characterization. Lower-order bits must first exhibit well-resolved plateau states. The carry pulse generated after a reset must then complete the higher-order write operation and cease before the next read window. Pulse width, mature integrase dynamics, RDF timing, and forward and reverse recombination fluxes should therefore be assessed over complete cycles rather than inferred from peak amplitudes alone.

**中文**

这些结果支持指定方程、参数及输入条件下的确定性三级计数。浓度换算系数与新增表达参数尚未实验标定；上游50 min倍增时间与部分下游总清除率的生长条件仍未对齐。实际正交酶与调控因子也需要各自表征，当前数值不能作为三个酶系统的实测参数直接移植。

**English**

The results support deterministic three-bit counting under the specified equations, parameter sets, and input conditions. The concentration-conversion factor and added expression parameters remain experimentally uncalibrated. The growth conditions underlying the 50 min upstream doubling time and some downstream total clearance rates have not yet been reconciled. The candidate orthogonal recombinases and regulatory factors also require individual characterization; the current numerical values should not be treated as experimentally measured parameters for three distinct enzyme systems.

<details>
<summary>模型来源、完整参数与复现记录 / Model Provenance, Complete Parameters, and Reproducibility Records</summary>

**中文**

完整连接方程见[HZH核心方程与符号附录](混合三级级联_核心方程与符号.md)，其中逐位区分正反重组、表达源到mRNA的换算，以及显式复合物的准稳态映射。HBY bit0／bit2的正向Hill指数为2，ZMH bit1的正向指数为4。

**English**

The [HZH core-equation and notation appendix](混合三级级联_核心方程与符号.md) specifies the bit-dependent forward and reverse recombination laws, the conversion from protein production targets to mRNA sources, and the quasi-steady-state mapping for explicit complexes. The forward Hill exponent is 2 for HBY bit0/bit2 and 4 for ZMH bit1.

**中文**

HZH中间模块采用《完整二级级联.py》的现行参数；HBY两端采用独立游离酶模型的早期基础表及选定扩展配置。源文件和单位转换见[主结果完整参数与来源](hby_zmh_hby/主结果_完整参数与来源.md)，连接代码见[model_hzh.py](hby_zmh_hby/model_hzh.py)，验收规则见[verify_hzh.py](hby_zmh_hby/certification/verify_hzh.py)。

**English**

The HZH middle module uses the current parameters from `完整二级级联.py`. The two HBY bit modules retain the early baseline parameter table and selected extensions of the independent free-integrase model. Source files and unit conversions are documented in the [complete parameter and provenance record](hby_zmh_hby/主结果_完整参数与来源.md). Module coupling is implemented in [model_hzh.py](hby_zmh_hby/model_hzh.py), and the evaluation criteria are defined in [verify_hzh.py](hby_zmh_hby/certification/verify_hzh.py).

**中文**

八初态结果虽保存于eight_phase_results目录，构造方式是连续八个稳态读窗覆盖八个不同数字；它与HZZ中已撤回的“同一时钟周期内八子相位、重启上游”的旧实验不同，不引用旧fig4_eight_phases或fig6_phase_events。

**English**

Although the valid HZH records are stored in a directory named eight_phase_results, their initial conditions span eight consecutive steady-state read windows and eight distinct digital values. This construction is distinct from the superseded HZZ experiment that sampled eight subphases within a single clock cycle and restarted the upstream driver. The superseded fig4_eight_phases and fig6_phase_events are not evidence for the present claim.

**中文**

读窗采用0.30／0.70带与80%占用率。主要事件阈值及其适用范围记录于验收报告；这些是分析约定，不能当作DNA写入的生化硬阈值。参数工作区数据见[60点时钟扫描](时钟工作区扫描/results/clock_axis_20261001_180035_369646/summary.json)。主结果轨迹的图源见[M-6](wiki_submission_assets/hybrid_threebit_300h.provenance.json)，时钟消融与输入配平的图源见[M-21](M21_时钟门消融/out/fig_M21_provenance.json)、[M-22](M21_时钟门消融/out/fig_M21_three_arm_provenance.json)。

**English**

Readout uses the 0.30/0.70 bands and an 80% occupancy requirement. Event-detection thresholds and their scope are specified in the verification report. These thresholds are analysis conventions rather than biochemical write thresholds. Operating-point data are available in the [60-point clock scan](时钟工作区扫描/results/clock_axis_20261001_180035_369646/summary.json). Figure provenance is provided in [M-6](wiki_submission_assets/hybrid_threebit_300h.provenance.json) for the principal trajectory, and in [M-21](M21_时钟门消融/out/fig_M21_provenance.json) and [M-22](M21_时钟门消融/out/fig_M21_three_arm_provenance.json) for clock ablation and input-integral matching.

</details>
