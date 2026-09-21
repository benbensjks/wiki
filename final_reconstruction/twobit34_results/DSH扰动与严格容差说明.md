# dsh：分子池扰动与严格容差扫描

## 冻结输入

不得修改：

- `selected_profile.json`；
- `initial_states_selected/four_initial_states.csv`；
- `model.py`、`model_twobit34.py`；
- `verify_twobit_causal.py`；
- `scan_twobit34_perturbations.py`；
- `scan_twobit34_strict_tolerance.py`。

建议继续使用 `D:\aconade\python.exe`。若 `Wait-Process` 被拒绝，可沿用文件系统完成检测：每片CSV达到预期行数且对应meta JSON出现后视为完成。

## A. 分子池扰动：严格按顺序运行

### A1. S0/S1

40点，10片，每片4点：

```powershell
python scan_twobit34_perturbations.py scan --group S --shard <0..9> --nshards 10 --hours 200
python scan_twobit34_perturbations.py merge --group S --nshards 10
```

### A2. RDF0/RDF1

40点，10片：

```powershell
python scan_twobit34_perturbations.py scan --group RDF --shard <0..9> --nshards 10 --hours 200
python scan_twobit34_perturbations.py merge --group RDF --nshards 10
```

### A3. A0/F0完整表达池

40点，10片：

```powershell
python scan_twobit34_perturbations.py scan --group AFFL --shard <0..9> --nshards 10 --hours 200
python scan_twobit34_perturbations.py merge --group AFFL --nshards 10
```

### A4. Int1完整表达池

20点，5片：

```powershell
python scan_twobit34_perturbations.py scan --group INT1 --shard <0..4> --nshards 5 --hours 200
python scan_twobit34_perturbations.py merge --group INT1 --nshards 5
```

扰动定义：

- S0/S1：`-0.20,-0.10,+0.10,+0.20,flip`，限制在0–1；
- RDF0/RDF1：对应 `M_R/R_u/R` 同时乘 `0.5,0.8,1.2,1.5,2.0`；
- A0/F0：对应mRNA、未成熟蛋白、成熟蛋白整池同时乘上述因子；
- Int1：`M_I1/I1_u/I1` 整池同时乘上述因子。

结果分类：

- `recovered_same_phase`：完整认证且恢复原计数相位；
- `stable_phase_shift_1/2/3`：完整认证但永久转入另一合法数字相位；
- `readout_pass_causal_fail`：码串能读但因果验收失败；
- `lost_counting`：失去稳定模4计数；
- `integration_failed`：数值失败。

S翻转导致合法相位改变不应自动视为故障；小幅扰动的首要目标才是恢复原相位。

## B. 中心与边界严格容差

9个参数点×3种数值设置=27条轨迹。使用9片，每片3条：

```powershell
python scan_twobit34_strict_tolerance.py scan --shard <0..8> --nshards 9 --hours 300
python scan_twobit34_strict_tolerance.py merge --nshards 9
```

三种设置：

1. baseline：`rtol=2e-7, atol=2e-9, max_step=2 min`；
2. tight：`rtol=1e-9, atol=1e-11, max_step=1 min`；
3. ultra：`rtol=1e-11, atol=1e-13, max_step=0.5 min`。

包含正式中心、4个通过边界和4个失败边界。必须比较完整冷启动/稳态码串、认证结论、失败原因及连续指标，而不只比较布尔值。

## dsh总结要求

### 扰动

1. 各组完整性、积分失败和非有限轨迹；
2. 每个目标在四种初态下的结果分类；
3. 保持原相位、合法相移、失锁分别多少；
4. 最大可恢复扰动幅度；
5. 哪个分子池最敏感；
6. S0/S1强制flip是否产生预期永久相移；
7. 是否出现读出正确但carry因果失败的危险点；
8. 所有分片的环境和源文件哈希一致性。

### 严格容差

1. 27/27完整性和积分健康；
2. 9个点是否全部保持原通过/失败结论；
3. 边界失败是否在严格容差下仍失败；
4. 每点baseline相对ultra的最大连续指标差；
5. 若出现判据变化，列出完整码串和具体差异；
6. 生成并复核结果文件SHA256。

## 解释限制

这些实验验证的是当前34状态ODE在选定参数和有限扰动下的吸引盆及数值稳定性，不替代实验测得的扰动幅度、细胞噪声或参数分布。
