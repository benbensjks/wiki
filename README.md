# 完整 11 状态三级生化 ODE 重建

> 2026-09-19 更新：当前实现入口为 [final_reconstruction/README.md](final_reconstruction/README.md)。新版使用真实上游、曾墨涵表内参数和 A0/R0/A1/R1 I1-FFL，共 43 个 ODE 状态。下文及根目录 `full_11state_3bit_ode.py`、`results/` 是旧 35 状态诊断版的历史说明，请勿与新版结果混用。

本文件夹重建 `前置生化背景/初步.pdf` 中每个正交重组酶模块的 11 个连续状态，并把三个模块连接为完全连续的三级模型。

## 模型层级

- 每个 bit：11 个状态。
- 三个 bit：33 个 PDF 基础状态。
- 连续 carry 扩展：2 个记忆蛋白 `E0`、`E1`。
- 总数：35 条 ODE。

模型没有以下机制：

- Python 逻辑 bit；
- 阈值事件产生脉冲；
- 事件队列；
- 直接执行 `q = 1-q`；
- 空闲期强制 `dS/dt = 0`；
- 根据目标方向人为关闭正向或逆向重组。

## 运行

本机可使用：

```powershell
$env:MPLBACKEND='Agg'
$env:PYTHONIOENCODING='utf-8'
& 'D:\aconade\python.exe' '.\full_11state_3bit_ode.py'
```

也可以指定输出目录：

```powershell
& 'D:\aconade\python.exe' '.\full_11state_3bit_ode.py' --output '.\results'
```

## 当前结果的正确读法

当前基准下，完整模型能够产生三级连续响应，所有三个 DNA 构象变量都能在 PB/LR 之间变化。但 `U0`、`U1`、`U2` 基本同频，只出现逐级相位延迟，没有产生 1:2:4 的频率除法。因此当前模型是一个连续级联的一次脉冲传播链，不是已经成功的 3-bit 二进制计数器。

这个失败结果不能用外部数字变量替换。下一步若要得到真正的二进制计数，必须在连续生化层加入能够保留奇偶状态的双稳或方向性交替机制，并重新验证，而不能在代码事件处理器中直接翻转逻辑位。

## 文件结构

- `full_11state_3bit_ode.py`：完整模型、分析和绘图。
- `方程与实现对照.md`：PDF 方程到代码的逐项映射。
- `results/parameters.json`：全部参数和单位。
- `results/trajectories.csv`：35 个状态、输入和解码结果。
- `results/clock_samples.csv`：逐拍采样。
- `results/pulse_summary.csv`：三级输入峰。
- `results/summary.json`：机器可读总结。
- `results/重建结果说明.md`：结果解释。
- `results/01_full_3bit_overview.png`：总体结果。
- `results/02_all_biochemical_states.png`：显式表达、成熟与复合物状态。
- `results/03_continuous_carry_module.png`：连续 carry 模块。

## 参数状态

目前参数用于结构验证和暴露机制问题，并非经过实验标定的最终参数集。特别是三个正交重组酶模块暂时共享同一套参数；这避免把未经证实的差异伪装成正交系统的实测差异，但后续必须按具体 Int/RDF 对分别标定。
