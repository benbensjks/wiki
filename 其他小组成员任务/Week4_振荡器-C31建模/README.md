# iGEM A 模块可复现建模交付包

本目录整理了 repressilator–TetR sponge–PLtetO1–phiC31 项目的阶段性正式代码、报告、代表性图片和结果摘要。代码由探索期 Notebook 脚本重构而来；本次整理只改变运行入口、文件路径和输出组织，不改变模型方程、参数或分析判据。

## 目录结构

```text
.
├── README.md
├── requirements.txt
├── iGEM_A模块建模阶段报告.md
├── code/                     # 六个可独立运行的正式程序及 B 模块依赖
├── data/reference_inputs/    # 固定的跨程序参考输入及来源说明
├── figures/                  # 报告使用的冻结图片，不会被新运行覆盖
├── results/                  # 已完成分析的冻结 CSV 摘要
└── outputs/                  # 本地重新运行后产生的时间戳输出目录
```

## 环境准备

推荐 Python 3.10–3.12。在项目根目录运行：

```powershell
python -m pip install -r requirements.txt
python code\verify_setup.py
```

`verify_setup.py` 只检查依赖、代码语法和必要输入文件，不运行模拟。

## 六个正式程序

| 程序 | 来源版本 | 回答的问题 | 主要输出 |
|---|---|---|---|
| `Mechanistic_Repressilator_C31_ODE.py` | v52a1 | 显式 mRNA–蛋白 ODE 是否能回到约化极限，并输出 C31 | 周期、QSS 误差、波形相关、C31 峰谷 |
| `Oscillator_Engineering_Sensitivity_Analysis.py` | v53a | 启动子、RBS、sponge 和倍增时间如何影响时钟 | OAT 周期、脉宽和剂量 |
| `Flux_Driven_Translation_Burden_Model.py` | v53d | 通量驱动的共享翻译资源负载会不会拖动时钟 | 负载–周期、幅度和剂量边界 |
| `Fixed_PLtetO1_Interface_Design_Scan.py` | v53g | C31 与 cI 使用同一 PLtetO1 时，RBS 与 Int 标签怎样匹配 B | RBS–标签成功窗口 |
| `Self_Consistent_AB_System_Robustness_Validation.py` | v53h | 候选设计在自洽负载、PB/LR 初态和长时间运行下是否可靠 | A/B 系统通过率和守恒审计 |
| `Stochastic_Sponge_Period_and_Noise_Analysis.py` | v54s0 | 随机模型中 sponge 是否产生可观周期延长 | 多周期估计器、bootstrap 区间、轨迹图 |

这些名称是交付名称；源版本号仍保存在代码头部和输出元数据中，便于追溯。

## 推荐运行顺序

每个程序都可直接从任意当前目录启动，因为路径由程序文件自身的位置确定：

```powershell
python code\Mechanistic_Repressilator_C31_ODE.py
python code\Oscillator_Engineering_Sensitivity_Analysis.py
python code\Flux_Driven_Translation_Burden_Model.py
python code\Fixed_PLtetO1_Interface_Design_Scan.py
python code\Self_Consistent_AB_System_Robustness_Validation.py
python code\Stochastic_Sponge_Period_and_Noise_Analysis.py
```

建议先单独运行需要查看的程序，不建议第一次使用时“一键运行全部”。A/B 长期验证和随机确认模式的计算量明显更大。

如果仍希望在 Notebook 中查看图像，可在项目根目录新建 Notebook，并运行：

```python
%run code/Mechanistic_Repressilator_C31_ODE.py
```

其余程序同理。正式交付和复现以 `.py` 入口为准，Notebook 仅用于交互式查看。

## 输入与输出规则

- 所有新运行结果写入 `outputs/<程序名_时间戳>/`，不会覆盖冻结结果。
- `figures/` 是报告展示用图片；它们对应已经完成的运行，不会随重新运行自动改变。
- `results/` 是阶段报告中结论对应的 CSV 摘要，可直接交给其他组员作图或审计。
- `Fixed_PLtetO1_Interface_Design_Scan.py` 使用 `data/reference_inputs/unloaded_C31_translation_trajectory.csv`，不再搜索某台电脑上的旧 v53d 文件夹。
- `zhao_core.py` 是 A/B 接口分析的本地依赖，放在 `code/` 中；请勿单独移动两个 A/B 程序而遗漏该文件。

## 结果解释边界

- 参数扫描是工程敏感性/边界分析，不是从实验数据拟合“最好看的参数”。
- 名义资源负载比例仍是未标定情景，不代表真实细胞负载。
- 随机周期比较以独立轨迹为统计单位；不同周期估计器对应不同可观测量，应同时报告而非只选择单一有利读数。
- `T_d=50 min`、PLtetO1 leak、C31 活性比例、RBS 比例和 Int 标签损失率中仍含工作假设，最终需由实验构建和测量确认。

## 当前阶段推荐候选

在固定使用与 cI 相同 PLtetO1 的条件下，当前最稳健的确定性工程候选为：

- `C31_RBS_scale = 0.45`
- `k_tag_int = 8 h^-1`

该候选在已扫描的 0–30% 名义动态资源负载、PB/LR 两种初态和长期确定性测试中均通过；这是一项待实验校准的设计建议，不是已验证的生物学常数。

## 本交付包的静态审计

- 六个正式程序均有 `if __name__ == "__main__":` 入口。
- 所有程序均使用相对项目路径，不依赖 `Path.cwd()`。
- 六个程序无 D:/E: 绝对路径。
- 必要本地依赖与参考输入已随包提供。
- 已完成 Python AST 语法检查；整理过程中未运行任何模拟。

