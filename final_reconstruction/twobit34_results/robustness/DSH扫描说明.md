# dsh：正式34状态模型75点联合鲁棒性扫描

## 任务边界

- 只运行 `scan_twobit34_robustness.py`；
- 不修改 `model.py`、`model_twobit34.py`、`verify_twobit_causal.py`、扫描脚本或 `ZENG`；
- 15个独立OS进程，每片5点，不使用 `ProcessPoolExecutor`；
- 每点300 h；
- 扫描完成后必须调用同一脚本的 `merge` 模式生成汇总和热图。

## 参数网格

```text
uM_per_au = 5.5, 5.75, 6.0, 6.25, 6.5
A0/F0共同成熟半衰期 = 25, 27.5, 30, 32.5, 35 min
carry mRNA半衰期 = 1, 2, 4 min
```

共75点。只改变这三个新增参数；其余采用 `model_twobit34.nominal_extension()`，曾墨涵表内参数保持不动。

## PowerShell启动参考

```powershell
$root = 'C:\Users\18633\Desktop\wiki\final_reconstruction'
$py = 'D:\aconade\python.exe'
$out = Join-Path $root 'twobit34_results\robustness'
$processes = @()

0..14 | ForEach-Object {
    $shard = $_
    $stdout = Join-Path $out ("shard{0:D2}.out.log" -f $shard)
    $stderr = Join-Path $out ("shard{0:D2}.err.log" -f $shard)
    $processes += Start-Process -FilePath $py -WorkingDirectory $root -WindowStyle Hidden -PassThru `
        -ArgumentList @('scan_twobit34_robustness.py', 'scan', '--shard', $shard,
                        '--nshards', 15, '--hours', 300) `
        -RedirectStandardOutput $stdout -RedirectStandardError $stderr
}

$processes | Wait-Process
$processes | Select-Object Id, ExitCode
& $py (Join-Path $root 'scan_twobit34_robustness.py') merge --nshards 15
```

## 自动交付物

`merge` 成功后应产生：

- `robustness_all.csv`：75点总表；
- `robustness_summary.json` / `.md`：完整性、失败分类、通过率和最佳点；
- `robustness_heatmaps.png` / `.pdf`：三种mRNA半衰期下的认证和时间裕量；
- `SHA256SUMS.json` / `.txt`：分片、汇总、模型和验证器哈希。

## dsh需要总结的内容

1. 75行/75唯一组合、重复/缺失/积分失败/非有限值；
2. 冷启动通过、稳态通过、完整因果认证数量；
3. 各失败原因数量；
4. carry mRNA半衰期1/2/4 min各自通过率；
5. 通过区是否连续，以及其边界；
6. 按 `min(bit1 setup, bit1 hold)` 排名前5的工作点；
7. 标称点 `6/30/2` 是否复现，裕量是多少；
8. 15份元数据的Python/依赖版本和四个源文件哈希是否完全一致；
9. SHA256清单重新计算是否0处不符；
10. 不把本扫描解释成gamma已含稀释的实验依据；它只验证该假设下的工程鲁棒性。
