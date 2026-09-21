# 最终结构重建版：真实上游 + 11 状态 bit + I1-FFL

> **2026-09-19 追加：单比特诊断结论见
> [`bit0_results/verification/bit0_诊断报告.md`](bit0_results/verification/bit0_诊断报告.md)。**
> 要点：`bit0` 在 `uM_per_au≈6–10`、`add_growth=False` 时是**有驻留、有记忆的 0↔1 交替**，
> 原 `bit0_diagnostic.score_cycles` 的“单点采样 0.8/0.2”判据会把真比特读成 `x`，应换成驻留判据
> （`verify_bit0_part4.dual_code`）；归档的 `bit0_results/scan_summary.json` 用仓内当前代码无法复现。
> 在 `uM_per_au=6` 上，两比特两三百小时读数连续得到 `1230` 循环。以下正文仍是结构说明。

本目录落实 2026-09-19 的五项要求。结构与代码可运行；是否形成 000→001→…→111 的连续计数，以 `results/summary.json` 的实测结果为准。这里“最终”指本次指定结构的交付，不表示参数已实验标定或计数已经通过。

## 五项要求与代码对应

| 要求 | 实现 |
|---|---|
| 韩亚轩的真实上游 | `model.py` 直接导入交付包 v53d，复用其 `rhs`、`initial_state`、`make_reference_parameters`；没有人工周期函数 |
| 曾墨涵参数 | `ZENG` 保留其参数分析 PDF §1–3 的数值；新引入参数只放在 `Extension` |
| 删除 E0/E1，采用 A0,R0,A1,R1 | 新模型仅含 A0,F0,A1,F1 四个 carry 状态；F0/F1 就是原文 R0/R1，避免与 bit 内 RDF 的 R 重名 |
| 下一级 M_I,I_u,I 承接 carry | I1-FFL 的启动子活性调控 M_I 转录，然后翻译、成熟；没有额外 Int1/Int2 状态或无量纲 U 整合酶 |
| A1 负自馈、外部时钟 AND、真实 C31 通量 | A1 产生项乘负自馈 Hill；bit2 M_I 产生项乘成熟游离 Int0 的 Hill；bit0 接收真实 C31 翻译通量 |

## 一次运行

```powershell
cd C:\Users\18633\Desktop\wiki\final_reconstruction
& 'D:\aconade\python.exe' .\run.py --hours 150
& 'D:\aconade\python.exe' -m unittest test_model -v
```

依赖：numpy、scipy、pandas、matplotlib。必须连同 `wiki/其他小组成员任务/Week4_振荡器-C31建模` 一起保留，代码通过相对目录找到真实上游。Python 3.10+；该机默认 `python` 没有 scipy，所以上述命令使用已经具备依赖的环境。

`run.py --cold --output results_cold` 可检查零调控蛋白初态；`--strict --output results_strict` 可收紧精度。默认 150 h、每分钟输出一次、DOP853、rtol=2e-8、atol=2e-10、最大步长 1 min。严格配置最大步长 0.5 min、容差缩小 100 倍。输出选项应设在本 wiki 文件夹内。

## 状态与单位

总计 **43 个连续 ODE 状态**：

1. 6 个上游状态：mTetR、TetR_total、mCI、CI、mLacI、LacI，单位 copies/cell。
2. 每个 bit 均保留 `(M_I,I_u,I,M_T,T_u,T,M_R,R_u,R,C,S)`，共 33 个状态。S 为 LR 比例，PB=1−S；其余为 a.u.。
3. A0、F0、A1、F1 共 4 个调控蛋白，单位 a.u.。

上游原来的 mC31 与 bit0 的 M_I 是同一个物理量，通过 `602.214076 × uM_per_au` 换算，不各积分一份。上游原来的 C31 蛋白积累方程由 bit0 的 I_u/I/C 反应替代。给原上游 rhs 的 C31 占位值不作为动态状态，且不影响原方程其他七个导数。

全模型时间单位 h。上游源文件按 min 计算，rhs 统一乘 60。细胞体积沿用上游 1 fL。

**a.u. 到 µM 的转换尚未标定**：默认 `uM_per_au=1` 是显式接口假设，不是实验事实；源模型 a.u. 数字全部保留。mRNA 的 a.u. 也作为数量尺度使用，便于表达链计算，不代表已测定的转录本浓度。

## 真实上游与翻译通量

从源代码保留 TetR→CI→LacI→TetR 抑制环、tetO sponge 的准平衡占据和 TetR 总量守恒。参数保留 Td=50 min、lambda=1000、K=13 copies、n=3、No=10、Nt=40、mRNA 总半衰期=2 min、翻译速率=0.5/min、C31 leak=0.005。

使用冻结 unloaded 基准（rho=1）：

`J31_uM_h = 60 * beta31_per_min * M_I0_au * uM_per_au`

`dI_u0/dt = J31_uM_h/uM_per_au - (kmat + delta_u + mu)*I_u0`

这个 J31 是翻译产生通量。M_I0 的转录来自原上游 PLtetO1 方程，不对 J31 再套一遍 Hill/转录滤波。I_u→I 成熟属于此次 11 状态扩展，因此成熟游离 Int0 不必与源文件中忽略成熟/结合的 C31 蛋白浓度一致。冻结 CSV 只用于验证，运行时不插值或重放 CSV。

王厚骅相关筛选交付的 `GSA_RDF/GSA_RDF/parameters2.py` 采用 38 状态 Zhao B 模型，工作点含 K_rep=0.0186、n_rep=3.4、krep_tsl=15、krdf_tsl=200、k_tag_int=12。其方程/单位与曾墨涵模型不同，不能把这些值混入本次要求的曾墨涵基准。后续周期旋钮报告也提示：确定性接口通过不能替代随机接口认证。本版未加入随机噪声。

## I1-FFL 与下游表达链

令 H(x;K,n)=x^n/(K^n+x^n)，G=1−H。F 是原文 carry 阻遏因子，R 是 RDF：

```
dA0/dt = 8*(1-S0) - (1.9+mu)*A0
dF0/dt = 3.8*H(A0;0.5,2) - (0.6+mu)*F0
dA1/dt = 16*(1-S1)*G(A1;0.6,2) - (2.5+mu)*A1
dF1/dt = 5*H(A1;1.2,4) - (1.1+mu)*F1

q1 = 18*H(A0;0.5,2)*G(F0;0.6,4)
q2 = 38*H(A1;1.2,4)*G(F1;0.4,4)*H(I0;K_gate,n_gate)
```

q1/q2 是目标成熟蛋白产生速率（a.u./h），不是新蛋白状态。两者分别通过下式进入现有 M_I1/M_I2：

```
lambda_m = delta_m + mu
lambda_u = kmat + delta_u + mu
alpha_tx(q) = q * lambda_m * lambda_u / (beta * kmat)
dM/dt = alpha_tx(q) - lambda_m*M
dP_u/dt = beta*M - lambda_u*P_u
dP/dt = kmat*P_u - (gamma_P+mu)*P + binding terms
```

如此在恒定 q 的转录/成熟准稳态下 `kmat*P_u=q`，保留原 18/38 等蛋白产生速率的含义。有限响应时间会改变动态，不声称与原粗粒化曲线完全相同。

Bit 内 T 表示曾墨涵的 Rep/BM3R1，与振荡器 TetR 是不同调控池：其目标产生速率 `qT=3.5*(1-S)`；RDF 目标产生速率 `qR=6*S*G(T;0.85,3.9)`。两者同样通过原 M_T/T_u/T、M_R/R_u/R 展开。这是对《初步.pdf》中 T 跟随共同输入项的有意结构替换，以满足本次曾墨涵模型要求。

Int0 在原文 AND 门中是时钟信号；这里取 bit0 成熟游离 I。原文只给了 gate 形式，未在所提供数值表中找到 K_gate/n_gate，默认 0.3/2 明确属于新增假设。该 Hill 门仍是有效调控描述；φC31 本身不是已指定的转录激活因子，实际构建还需选择传感 TF/杂合启动子或 split-integrase 方案。A1 的负自馈用于限制表达和泄漏，单靠负自馈并不能保证产生脉冲或计数成功。

## 生化复合物与 DNA 动力学

保留《初步.pdf》的 1:1 显式结合：

```
binding = kon*I*R; unbinding = koff*C
dI += -binding+unbinding
dR += -binding+unbinding
dC = binding-unbinding-(delta_C+mu)*C
dS = 7*H(I;KD_int,2)*0.1/(0.1+R)*(1-S) - 5*H(C;KC,2)*S
```

KD_int 对 bit0/1/2 分别为 1、1、0.6。保留曾墨涵正向反应的经验 RDF 抑制项；这与显式游离 Int 减少共同起作用，是待验证的组合假设，不把它冒充独立测得的第二种抑制机制。

曾墨涵反向速率的自变量是 **I×R**，《初步.pdf》是显式 C，二者的阈值不能直接照搬。令 `qC=kon/(koff+delta_C+mu)`，则 C 的准稳态为 qC*I*R，故 `KC=qC*1.2`。这样反向 Hill 在该准稳态下等于她的原速率函数。原文 K_D_comp 标为 a.u.，但用于 I×R 的表达式要求 a.u.^2；这里明确作量纲修正，而非将 1.2 当作结合解离常数。保留原数值与导出 KC 同时输出。

结合本身严格满足 I+C、R+C 收支；复合物损失按 I、R 两者共同移除处理。该降解命运也属新增假设，不是 DNA 损失；S+PB 恒为 1。

## 参数一致性的边界

原文给的是粗粒化蛋白模型，并未给出完整 11 状态所需的 mRNA/成熟/复合物动力学参数。当前目录内未找到可明确署名为曾墨涵的完整原始 Python 文件；本版参数以她的两份 PDF 为直接依据，旧 `mixed_2bit.py` 仅作交叉核对。因此可核实的是“与提供的数值表一致”，不能冒称已经逐行复现她的原始代码。

| 参数来源 | 处理 |
|---|---|
| 曾墨涵表内全部核心/两级 FFL 参数 | 原值保留在 ZENG，包括 k_int=6 |
| k_int=6 的方波输入 | 曾墨涵已确认只是为得到正常Int0方波脉冲的专用调参量；真实C31通量不可与其直接相乘，故弃用该因子 |
| gamma_I=(2,1.4,2.2)、gamma_T=0.6、gamma_R=0.8 | 曾墨涵已确认是“降解+稀释”的总清除率；不得再叠加 growth |
| mRNA 半衰期 2 min、成熟半衰期 5 min、beta=30/h | 新增展开参数；无下游实验标定 |
| kon=1/(a.u. h)、koff=10/h、delta_C=1/h | 新增显式复合物参数；不由 K_D_comp=1.2 唯一决定 |
| uM_per_au=1、K_gate=0.3、n_gate=2 | 接口/门控假设，需实验或原始代码补充 |

曾墨涵已确认 gamma 是“本征降解+生长稀释”的总清除率，因此正式模型采用 `add_growth=False`，避免重复计入稀释。`add_growth=True` 仅作为历史诊断模式保留，物理上会重复加一次稀释，不能作为候选工作点。尚未闭合的是两段模型的生长条件：韩亚轩上游 `mu≈0.8318/h` 大于曾表最小总清除率 `0.6/h`，所以两者不能在未知条件下直接分解成同一组非负本征降解率；需要曾墨涵补充其假设的 mu 或倍增时间。

## 初值和判据

默认三个位 PB；Rep、A、F 按 PB 时的调控稳态预装载，I/RDF 从零开始，mC31 沿用源振荡器初值。此初值需要实验预培养准备，并不是整个网络的稳态；保留源振荡器启动相位，报告所有启动段。冷启动可另运行。

计数读点为相邻真实通量峰之间的谷时刻，S≤0.3 解码为 0，S≥0.7 解码为 1，中间区记为 −1。解码、找峰、模 8 比较均只做事后分析，不进入积分。除完整顺序正确率，也输出逐拍模 8 增量正确率。至少 8 个完整读点且全部匹配才标记 binary_counter_pass。

## 结果与验证文件

- `results/trajectories.csv`：全部 43 个状态、真实通量、carry 启动子活性、各级成熟产生通量。
- `results/parameters.json`：源参数、新增假设、导出量、源文件 SHA256。
- `results/summary.json`、`clock_samples.csv`、`结果说明.md`：计数结果。
- `results/overview.png`、`all_33_states.png`：总体和全部 bit 状态图。
- `validation/numerical_validation.json`：100 h 冻结上游对照、严格容差对照的逐状态误差。
- `test_model.py`：状态契约、QSS 映射、AND/负自馈、物料收支、非负边界、上游隔离和数值回归。

## 可追溯原始资料

- [韩亚轩上游 v53d 代码](../其他小组成员任务/Week4_振荡器-C31建模/code/Flux_Driven_Translation_Burden_Model.py)
- [曾墨涵前馈环模型](../其他小组成员任务/前馈环脉冲进位级联模型.pdf)：第 5–6 页单 bit/I1-FFL/时钟 AND，第 4 页 A1 自反馈。
- [曾墨涵数值表](../其他小组成员任务/前馈环脉冲进位级联模型参数分析%20(1).pdf)：第 1–2 页。
- 原 11 状态结构：`C:/Users/18633/Desktop/final/前置生化背景/初步.pdf`。

本版只读其他成员原件，旧 35 状态诊断程序和旧结果仍在 wiki 根目录留作追溯；运行入口以本目录为准。

## 2026-09-19 两比特因果时序复核

新增 `verify_twobit_causal.py`，以 bit0 实际反向重组通量、绝对 carry 门剂量和有限读出窗口统一验证两比特。基线 `uM_per_au=6` 的同步读点序列正确，但 bit1 翻转早于主要反向重组事件约 0.73–1.55 h，故不再标为完整认证。

新增 `verify_carry_delay.py` 作为尚未并入基线模型的生化延迟变体：为 A0/F0 各追加 mRNA 和未成熟蛋白状态，不修改任何 `ZENG` 数值。`uM_per_au=6`、carry mRNA 半衰期 2 min、成熟半衰期 25–35 min 时，300 h 冷启动与稳态均通过；暂以 30 min 为标称候选。

`verify_carry_phase.py` 进一步将稳定上游平移 0、2、4、6、8、10 h，并把下游重置为相同 PB 初态；六组 300 h 均通过，排除了唯一上游初相位造成的偶然通过。

## 正式34状态两比特模型

`model_twobit34.py` 将验证通过的A0/F0表达延迟固化为正式34状态ODE：6个上游状态、两个原始11状态bit、
A0/F0成熟蛋白以及四个A0/F0转录/成熟状态。`verify_twobit34_equivalence.py` 证明其与临时验证变体在300 h
内逐状态、逐信号和逐事件完全一致；`run_twobit34.py` 生成完整轨迹、读窗、参数和五组报告插图。

正式结果与使用说明见 `twobit34_results/正式34状态模型说明.md`。

`verify_twobit34_initial_states.py` 从稳定周期轨道提取四个完整生化初态，并把上游统一到同一时钟相位。
`00/01/10/11` 四种初态分别续跑300 h后均严格按模4计数并通过完整因果认证。结果和图见
`twobit34_results/initial_states/四初态验证报告.md`。

75点联合扫描得到53个完整认证点。严格要求三个扫描轴均不在边界且六个轴向邻居全部通过后，正式中心选为
`uM_per_au=5.75`、A0/F0成熟半衰期32.5 min、carry mRNA半衰期2 min；四初态在该中心重新验证后全部通过。
选点依据见 `twobit34_results/selected_profile.md`，复核热图见 `twobit34_results/robustness/review_summary.md`。

完整方法、边界和结果见 `bit0_results/verification/两比特因果时序与延迟层报告.md`。曾墨涵已确认 gamma 包含生长稀释，因此 `add_growth=False` 的结构解释成立；仍需补充她使用的具体 mu，以解决与韩亚轩50 min倍增条件之间的不一致。

## 2026-09-20 显式生长稀释工作区

原800点重扫发现的 `add_growth=True,uM_per_au≈1` 区域在数学上确实能形成单bit0周期二分频，但在gamma已确认包含稀释后，
该模式等于重复计入稀释，只作为诊断区域保留，不再作为物理候选。
但它在曾墨涵I1-FFL表内参数不变时无法驱动bit1：直接接口及A0/F0分别成熟的34个组合均失败，bit1最高约0.549。
其数学与机制结果见 `bit0_results/verification/growth_true后续验证报告.md`，解释口径以本节最新确认优先。

## 51状态三级扩展

`model_threebit51.py` 在冻结34状态前缀之后追加bit2的11状态和完整A1/F1表达链，总计51状态。
数值Jacobian证明下游对前34状态无反馈；300 h投影轨迹与冻结两比特模型在 `1.03e-8` 内一致。
首轮600 h未调参诊断显示14个bit1反向事件产生14个干净g1窗口，但bit2仅完成3次阈值穿越，之后因正反重组通量抵消停在高态周期轨道。
详细结果见 `threebit51_results/51状态骨架与首轮诊断.md`。
