# 硬约束（改代码前必读）

这些是本项目用真实故障换来的约束。违反任何一条都会**静默**产生看起来正常的结果，
所以它们不是风格建议，而是必须先检查的前置条件。

---

## H1. 构造 51 状态模型时，`n_A1_gate` 必须显式给出

```python
ThreeBit51Model(..., n_A1_gate=None)   # → n_A1_gate_effective 回落到 ZENG['n_A'][1] = 4.0
```

`model.py` 的 `ZENG['n_A'][1]` 同时被用在**两处**：carry 门（`model.py:120`）与 F1 产生项
（`model_threebit51.py:96`）。`n_A1_gate` 只覆盖**门臂**，它是相对已发表参数用法的**实质性
修改**，取值 6.0（见 `threebit51_selected_v1.json`）。

**已因此作废的分析**：`check_threshold_margin`、`check_absolute_scale`、`check_mass_balance`、
`check_rdf_collapse`、`check_identifiability`(旧)、`scan_sharpen_vs_lengthen` arm A、
以及局部可辨识性第一次运行的 55 次中的 49 次。详见 `INVALIDATED_ARTIFACTS.json`。

**检查方法**：构造后立刻断言

```python
assert float(model.n_A1_gate_effective) == 你想用的值
```

或使用 `working_point.working_point_block(model)`，它**读回**实测值而不复述意图。

## H2. 绝不同时持有两个语义依赖 `ZENG` 的模型实例

`ThreeBit51Model._carry1_sources` 是 **staticmethod，在调用时读全局 `ZENG['n_A'][1]`**。
因此一次全局补丁会**渗漏到所有存活的模型实例**，包括那些本应使用表值的实例。

**症状**：等价性测试出现只集中在某一个索引（如 `M_F1`/49）的假差异。
**同源故障**：`check_rdf_collapse` v1 全局改写 `ZENG['n_rep']/['K_rep']` → 冻结前缀被毁
（`reverse_events` 14 → 0）。

**规则**：构造、积分、比较必须**严格串行**，并在相位之间 `restore(patch)`。
`build_threebit(n_A1_gate=...)` 走构造函数、不带此风险；`build_decoupled` 走全局补丁、
带此风险（两者已被验证语义等价到舍入）。

## H3. 只有 `invalidated_n4/` 里的东西是经过归档的失效产物

任何被判定失效的产物，**先把原件移入归档目录再重跑**，不要就地覆盖，更不要删目录。

**已发生的损失**：局部可辨识性第一次运行的 55 份 n=4 产物与原始撤回记录被 `Remove-Item`
删掉，未先归档；现存的 `RETRACTION_n4_run.md` 是**重建版**，不能声称字节一致。

## H4. 泄漏指标只能用配对 R/F 超周期版本

单窗 `median(far_off_rev)/median(gate_on_rev)` 已停用（`threebit51_leak_correction.json`）：
进位按 R/F 型交替，奇数窗数时中位数落在单一物理类型上，实测离散可达 2.87e6 倍。

现行定义：按 S2 在进位开始时的值分类；一 R + 一 F 组成一个 8 拍超周期；**先在对内求和，
再跨配对统计**；tail 独立；报 1 %/5 %/10 % 三档。
**配对指标依赖输出网格**（2 min → 1 min 变化 3.9e-3），引用绝对数值必须声明网格。

## H5. 结构性检查用 RHS 级比较，不用轨迹相等

自适应步长控制器让同一模型的两次积分出现 ~1e-13 的舍入差。判据要写成相对量：

```
rhs_relative_gap = max|Δrhs| / max|rhs|     阈值 1e-12
```

**注意**：代数等价的两条实现路径**不会**给出逐位相同的 RHS（运算顺序不同 → 末位舍入不同）。
`decoupling_verification.json` 的 Test 3 就是这种情形：前缀（状态 0..33）gap **恰为 0.0**，
而状态 34..50 的 gap 是 4.4e-16（相对 8.8e-19）。**"恰好 0" 是不该用的判据。**

## H6. 秩、条件数这类阈值量必须带阈值上报

`effective_rank` 是阈值量。同一个矩阵：

```
1e-2 切点 → rank 7（理论） / 6（实验）
1e-3 切点 → rank 9 / 9
```

只报一个 rank 而不报切点与最小奇异值的余量，等于把结论建立在一个任意选择上。

## H7. 不要用 `Remove-Item -Recurse` 清理分析目录

见 H3。清理前先归档；`write_manifest()` 用 `rglob`，子目录会被纳入清单，归档不会破坏可追溯性。

## H8. 不要把 `tempfile.mkdtemp()` 的临时目录建在工作区内

实测（2026-09-24）：在工作区里 `tempfile.mkdtemp()` 建的目录带 owner-only 权限，沙箱随后
**连读都拒绝**（`scandir` 与 `rmdir` 都是 WinError 5），`os.chmod` 也修不回来 ——
留下一块**永远删不掉**的垃圾目录。

如果建在 `plausibility/` 里，后果不只是垃圾：`write_manifest()` 用 `rglob` 收录子目录，
而清单**永远无法与一个删不掉的目录达成一致**，`make_authority_manifest.py check` 会永久失败。

**规则**：自检需要临时树时

1. 用**固定路径** + `mkdir(parents=True, exist_ok=True)`，**不要**用 `mkdtemp`；
2. 放在 `final_reconstruction/` **之外**（根级清单是具名列表，不会误收），
   或在系统临时区（`$env:TEMP` 下的深度 3 建目录实测可行）；
3. 用完自己 `shutil.rmtree` 删掉，并**打印删除是否成功**；
4. 脚本本体不要留在项目里（放在临时区），避免污染交付物。

**已发生的损失**：`wiki/_selftest_scratch/mam_selftest_iub9mba2` 因第 1 条被违反而无法删除，
只能留在原地；该目录不在任何清单内，无害但已记录在此。
