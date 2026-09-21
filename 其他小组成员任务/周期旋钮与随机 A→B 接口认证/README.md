# iGEM A 模块周期旋钮与随机接口：工作汇总

本文件夹汇总 2026-08-24 阶段交付内容，重点包括：

- 机制化 ODE 的 Morris/Sobol 全局敏感性分析；
- RBS–mRNA 半衰期周期旋钮设计；
- 周期档位的确定性、随机和联合参数鲁棒性认证；
- B 模块最小恢复周期与确定性 A→B 接口边界；
- 随机 A→B 接口一致性和逐周期失败归因；
- sponge 500 条/条件随机轨迹确认。

## 文件夹结构

```text
code/       正式分析脚本及运行所需的本地依赖
figures/    报告引用的关键图片
outputs/    关键汇总 CSV/JSON；不包含大体积逐点轨迹
report/     Markdown 阶段报告
README.md   本说明文件
requirements.txt
```

## 直接查看结果

建议先阅读：

```text
report/阶段报告_周期旋钮与随机接口认证_20260824.md
```

报告中的图片使用相对路径，移动整个提交文件夹后仍可显示。

## 运行环境

- Python 3.12（项目当前环境）
- NumPy、Pandas、SciPy、Matplotlib

安装依赖：

```powershell
python -m pip install -r requirements.txt
```

所有脚本均应从本提交文件夹根目录运行，例如：

```powershell
python .\code\Mechanistic_ODE_Global_Sensitivity_Analysis.py --mode smoke
```

`smoke` 仅检查流程；正式统计应使用对应脚本支持的 `full`、`audit` 或 `confirmation` 模式。随机确认和 GSA 计算时间较长。

## 主分析顺序

1. `Mechanistic_ODE_Global_Sensitivity_Analysis.py`：Morris 筛选；
2. `Postprocess_Morris_Interface_Leak_Metrics.py`：统一 leak 与接口指标口径；
3. `Shared_PLtetO1_Sobol_Global_Sensitivity.py`：Sobol 主效应/总效应；
4. `Shared_PLtetO1_RBS_mRNA_2D_Design_Window.py`：生成二维周期设计图；
5. `Shared_PLtetO1_Multiobjective_Period_Knob_Designer.py`：候选周期档位选择；
6. `Stochastic_Period_Knob_Library_Certification.py`：随机 A 时钟认证；
7. `Period_Knob_Robustness_Certification.py`、`Period_Knob_Robustness_Failure_Attribution.py`：联合不确定性认证与失败归因；
8. `Robust_Period_Knob_Library_Redesign.py`：鲁棒周期库重设计；
9. `Robust_Continuous_Period_Knob_AB_Certification.py`、`B_Minimum_Recovery_Period_Map.py`：确定性 A→B 与 B 恢复边界；
10. `Stochastic_Deterministic_Mean_Consistency_Audit.py`、`Tau_Leap_Time_Step_Convergence_Audit.py`：随机数值一致性；
11. `Stochastic_AB_Interface_Consistency_Audit.py`、`Stochastic_AB_Cycle_Failure_Attribution.py`：随机接口审计与失败归因；
12. `Stochastic_Sponge_500_Trajectory_Confirmation.py`：sponge 500 条/条件确认。

## 支持文件

以下文件主要提供函数或数据结构，通常无需单独运行：

- `Shared_PLtetO1_Period_Knob_Design_Map.py`
- `Arm_Specific_Promoter_Leak_Sensitivity.py`
- `Robust_K4_K5_Stochastic_AB_Fidelity_Certification.py`
- `zhao_core.py`

`Period_Knob_AB_Interface_Certification.py` 和 `Consolidate_Certified_Period_Knob_Library.py` 保留用于复现早期 T10–T13 库及兼容后续脚本，但其结果已由鲁棒连续档位流程更新。

## 当前主要结论

- 标称 T10–T13 随机 A 时钟保持单调，周期 CV 约为 3.4%–3.7%；
- 联合不确定性下，$P_R$ leak 是周期排序的主要风险参数；
- 重设计 K4/K5 在确定性 A→B 接口中通过；
- 随机 raw A→B 严格 PB/LR 配对周期成功率为 K4 50.0%、K5 60.7%；
- 接口失败以方向/相位错误为主，单一峰值、FWHM、剂量或周期指标的区分能力有限；
- sponge 500 条/条件确认中，phase 周期比为 1.242，PSD 周期比为 1.113，并伴随周期 CV 降低。

## 数据范围

`outputs/` 仅保留支持报告结论的汇总 CSV/JSON。大体积原始逐点轨迹仍保存在原项目：

```text
D:\CodexWork\iGEM_A模块_Week1模型审计\41_Mechanistic_ODE_Global_Sensitivity_Analysis_20260819
```

