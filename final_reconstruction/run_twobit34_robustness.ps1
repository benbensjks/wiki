# dsh driver: 75-point robustness scan of the formal 34-state two-bit model.
# 15 independent OS processes, 5 points each, no ProcessPoolExecutor.
# Follows twobit34_results/robustness/DSH扫描说明.md.
#
# This sandbox denies Wait-Process on child handles, so completion is detected
# from the filesystem (each shard writes its 5-row CSV after every point and its
# meta JSON as the last act).  The driver therefore stays alive for the whole
# scan, which also keeps the children alive.
#
#   powershell -ExecutionPolicy Bypass -File .\run_twobit34_robustness.ps1 scan
#   powershell -ExecutionPolicy Bypass -File .\run_twobit34_robustness.ps1 merge
#
param([ValidateSet('scan', 'merge')][string]$Mode = 'scan')

$ErrorActionPreference = 'Continue'
$env:MPLBACKEND = 'Agg'
$env:PYTHONIOENCODING = 'utf-8'
$env:OMP_NUM_THREADS = '1'
$env:OPENBLAS_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'

$root = 'C:\Users\18633\Desktop\wiki\final_reconstruction'
$py = 'D:\aconade\python.exe'
$out = Join-Path $root 'twobit34_results\robustness'
$nshards = 15
$points_per_shard = 5
$hours = 300
$deadline = (Get-Date).AddHours(6)
$started = Get-Date

New-Item -ItemType Directory -Force -Path $out | Out-Null

if ($Mode -eq 'merge') {
    & $py (Join-Path $root 'scan_twobit34_robustness.py') merge --nshards $nshards
    Write-Host "merge exit=$LASTEXITCODE"
    exit $LASTEXITCODE
}

Write-Host "robustness scan start $($started.ToString('s')) shards=$nshards hours=$hours"
$launched = @()
0..($nshards - 1) | ForEach-Object {
    $shard = $_
    $stdout = Join-Path $out ("shard{0:D2}.out.log" -f $shard)
    $stderr = Join-Path $out ("shard{0:D2}.err.log" -f $shard)
    $p = Start-Process -FilePath $py -WorkingDirectory $root -NoNewWindow -PassThru `
        -ArgumentList @('scan_twobit34_robustness.py', 'scan', '--shard', $shard,
                        '--nshards', $nshards, '--hours', $hours) `
        -RedirectStandardOutput $stdout -RedirectStandardError $stderr
    $launched += [pscustomobject]@{ shard = $shard; pid = $p.Id }
    Start-Sleep -Milliseconds 400
}
Write-Host ("launched {0} shard processes: {1}" -f $launched.Count,
            (($launched | ForEach-Object { "s$($_.shard)=$($_.pid)" }) -join ' '))

function Get-ShardState {
    $rows = 0; $done = 0; $stalled = @()
    for ($s = 0; $s -lt $nshards; $s++) {
        $tag = '{0:D2}of{1:D2}' -f $s, $nshards
        $csv = Join-Path $out "robustness_shard$tag.csv"
        $meta = Join-Path $out "meta_shard$tag.json"
        if (Test-Path $csv) {
            $n = (Import-Csv $csv).Count
            $rows += $n
            $age = ((Get-Date) - (Get-Item $csv).LastWriteTime).TotalMinutes
            if ($n -lt $points_per_shard -and $age -gt 25) { $stalled += "s$s(rows=$n,idle=$([int]$age)m)" }
        }
        if (Test-Path $meta) { $done++ }
    }
    return [pscustomobject]@{ rows = $rows; done = $done; stalled = $stalled }
}

$lastRows = -1
while ($true) {
    Start-Sleep -Seconds 60
    $st = Get-ShardState
    $elapsed = [int]((Get-Date) - $started).TotalMinutes
    if ($st.rows -ne $lastRows -or $st.done -ge $nshards) {
        Write-Host ("t+{0,3} min  metas={1}/{2}  rows={3}/{4}{5}" -f $elapsed, $st.done, $nshards,
                    $st.rows, ($nshards * $points_per_shard),
                    $(if ($st.stalled.Count) { "  STALLED: " + ($st.stalled -join ',') } else { '' }))
        $lastRows = $st.rows
    }
    if ($st.done -ge $nshards) { break }
    if ((Get-Date) -gt $deadline) { Write-Host 'DEADLINE reached'; break }
}

$final = Get-ShardState
Write-Host ("scan finished: metas={0}/{1} rows={2}/{3}" -f $final.done, $nshards, $final.rows,
            ($nshards * $points_per_shard))
for ($s = 0; $s -lt $nshards; $s++) {
    $tag = '{0:D2}of{1:D2}' -f $s, $nshards
    $csv = Join-Path $out "robustness_shard$tag.csv"
    $n = if (Test-Path $csv) { (Import-Csv $csv).Count } else { -1 }
    Write-Host ("shard {0:D2} rows={1}" -f $s, $n)
}
$badErr = Get-ChildItem $out -Filter 'shard*.err.log' | Where-Object { $_.Length -gt 0 }
Write-Host ("non-empty stderr logs: {0}" -f $(if ($badErr) { ($badErr.Name -join ',') } else { 'none' }))
Write-Host "robustness scan done, elapsed $([int]((Get-Date) - $started).TotalMinutes) min"
