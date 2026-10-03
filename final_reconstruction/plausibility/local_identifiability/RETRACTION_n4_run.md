# 撤回记录（重建版）：局部可辨识性第一次运行

> **这份文件是重建的，不是原件。**
> 原始 `ANALYSIS_STATUS.json` 与 `RETRACTION.md` 写在
> `plausibility/local_identifiability/` 里；我在重跑前用
> `Remove-Item local_identifiability -Recurse -Force` **清空了该目录，把这两份记录
> 连同 55 次 n=4 运行的 npz/json 一起删掉了**，归档目录 `invalidated_n4/` 里也没有副本。
> 本文件依据本次会话记录重建，内容忠实于原件，但**不能声称与原件逐字节一致**。
> 该删除本身是一次流程失误，已记入下方 `process_failures`。

**状态：`FAILED_ACCEPTANCE` — 第一次运行的分析结果不得被引用。**

## 1. 已确认的错误：工作点混用

`check_local_identifiability.build_with`（旧版）写成：

```python
ext, carry1, n_gate = frozen_extension(), frozen_carry0(), None
for name, value in overrides.items():
    ...
    elif k == 'ctor':
        n_gate = float(value)      # 只有 n_A1_gate 被扰动时才赋值
```

而 `model_threebit51.py`：

```python
return ZENG['n_A'][1] if self.n_A1_gate is None else float(self.n_A1_gate)
```

`NOMINAL['n_A1_gate'] = 6` 只是字典常量，**从未进入模型构造**。实测复核：

| 运行 | `overrides` | `n_A1_gate_effective` |
|---|---|---:|
| baseline | `{}` | **4.0** |
| 另外 8 个参数各 ±1 %/±2 %/±5 %，共 48 次 | 仅该参数 | **4.0** |
| `n_A1_gate` 自身 6 次 | 含 `n_A1_gate` | 5.70 – 6.30 |

**49 次在 n=4，6 次在 n≈6。** 后果：拼出的 9 列矩阵**不是同一参数点上的 Jacobian**，
其秩、条件数、最弱方向、全部相关系数**一律无效**。

这与我在审核旧 `check_identifiability.py` 时指出的缺陷是同一个——我把同样的错误在新脚本里又犯了一次。

## 2. 撤回的结论

1. **"相位伪影是已定位根因"** —— 撤回。固定时刻差分中翻转时刻变化是
   `∂S/∂p = ∂s/∂p − (∂s/∂t)(∂τ/∂p)` 的第二项，是**真实时序灵敏度**。σ₁ 随步长缩小上升
   不足以证明 `1/h` 发散；不同观测集给出不同相关矩阵是常态；且存在已确认的工作点混用必须先排除。
   正确表述：*在当前实现与扰动尺度下未得到可信的局部灵敏度估计；工作点混用是已确认错误，
   时序灵敏度与有限幅度非线性是待区分的因素。* 不得据此宣布"局部灵敏度方法对本系统不可用"。
2. **"`q·K_D_comp` 精确结构退化"** —— 撤回。`K_complex = q·K_D_comp` 使 q 从反向 Hill
   约掉**只在准稳态 `C ≈ q·I·R` 下成立**；模型显式积分 C，`kon/koff/δ_C` 还决定复合物弛豫时间、
   游离 Int/RDF 隔离量、复合物损失。判据不能只看阈值定义。
3. **"仍然站得住的三条"** —— 全部撤回：`n_A1_gate` 与 `uM_per_au` 不混淆（两列来自不同工作点）、
   mRNA 与 A1 成熟强负相关（对应 n=4 配置）、q↔K_D 精确退化（见上）。
4. "有效秩 8" 是阈值 `σ/σ_max > 1e-3` 下的有效秩，不是数学秩。
5. 收敛检查口径有误：比较的是按时间压缩后的 **RMS 摘要表**，不是完整的
   "时间 × 观测量 × 参数"矩阵；报出的 0.507 / 0.680 应称"摘要差异"。

## 3. 仍然成立的部分（当时保留）

- 收敛检查确实报警：三个步长下奇异谱与条件数相差约 15 倍，那批运行**没有给出收敛的局部估计**。
- 冻结与冻结后结果**完全不受影响**：其余脚本都显式传 `n_A1_gate`。
- 旧 `check_identifiability.py` 同样跑在 ZENG 默认值上，但**从未执行**，无产物。

## 4. 当时已实施的修正

`resolved_params` 先解析完整参数向量；`n_A1_gate` 每次显式传入并断言；每个 job json 记录
`resolved_params` 与 `n_A1_gate_effective`；`load_jobs` 加硬守卫；收敛判据改为完整矩阵 + RMS 分别报告。

## 5. process_failures（本次补记）

1. **记录被误删**：重跑前清空目录，销毁了原始撤回记录与 55 份 n=4 产物。正确做法是把它们
   移入归档目录（`invalidated_n4/`）后再重跑，如我处理另外 4 个检查时那样。本文件是补建的。
2. 该轮结束时 `INVALIDATED_ARTIFACTS.json` 的 `already_marked` 仍指向已不存在的两份文件，
   已一并更正（见该索引的 `already_marked` 字段）。

## 6. 结论

本目录当前的分析结果来自**重跑**（统一 n=6，审计通过），其验收状态见
`ANALYSIS_STATUS.json`。**本文件只描述已被取代的第一次运行，不适用于当前结果。**
