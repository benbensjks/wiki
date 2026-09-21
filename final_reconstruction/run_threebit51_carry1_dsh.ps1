# dsh driver: stage-1 carry1 scan of the 51-state three-bit model.
# 7 independent OS processes, 5 points each, 600 h per point.
#
# This sandbox denies Wait-Process on child handles, so completion is detected
# from the filesystem: every shard writes its 5-row CSV after each point and a
# meta JSON as its last act.  The driver stays alive for the whole scan, which
# also keeps the children alive.
#
#   powershell -ExecutionPolicy Bypass -File .\run_threebit51_carry1_dsh.ps1 scan
#   powershell -ExecutionPolicy Bypass -File .\run_threebit51_carry1_dsh.ps1 merge
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
$script = Join-Path $root 'scan_threebit51_carry1_dsh.py'
$out = Join-Path $root 'threebit51_results\carry1_scan_dsh'
$nshards = 7
$points_per_shard = 5
$hours = 600
$deadline = (Get-Date).AddHours(4)
$started = Get-Date

New-Item -ItemType Directory -Force -Path $out | Out-Null

if ($Mode -eq 'merge') {
    & $py $script merge --nshards $nshards --hours $hours
    Write-Host "merge exit=$LASTEXITCODE"
    exit $LASTEXITCODE
}

# record the scanner hash before starting, so it can be compared afterwards
$hashBefore = (Get-FileHash -Algorithm SHA256 $script).Hash
Write-Host "carry1 scan start $($started.ToString('s')) shards=$nshards hours=$hours"
Write-Host "scanner sha256 before: $hashBefore"

$launched = @()
0..($nshards - 1) | ForEach-Object {
    $s = $_
    $p = Start-Process -FilePath $py -WorkingDirectory $root -NoNewWindow -PassThru `
        -ArgumentList @('scan_threebit51_carry1_dsh.py', 'scan', '--shard', $s,
                        '--nshards', $nshards, '--hours', $hours) `
        -RedirectStandardOutput (Join-Path $out ("shard{0:D2}.out.log" -f $s)) `
        -RedirectStandardError (Join-Path $out ("shard{0:D2}.err.log" -f $s))
    $launched += [pscustomobject]@{ shard = $s; pid = $p.Id }
    Start-Sleep -Milliseconds 400
}
Write-Host ("launched {0} shard processes: {1}" -f $launched.Count,
            (($launched | ForEach-Object { "s$($_.shard)=$($_.pid)" }) -join ' '))

function Get-ShardState {
    $rows = 0; $done = 0; $stalled = @()
    for ($s = 0; $s -lt $nshards; $s++) {
        $tag = '{0:D2}of{1:D2}' -f $s, $nshards
        $csv = Join-Path $out "carry1_shard$tag.csv"
        $meta = Join-Path $out "meta_shard$tag.json"
        if (Test-Path $csv) {
            $n = (Import-Csv $csv).Count
            $rows += $n
            $age = ((Get-Date) - (Get-Item $csv).LastWriteTime).TotalMinutes
            if ($n -lt $points_per_shard -and $age -gt 45) { $stalled += "s$s(rows=$n,idle=$([int]$age)m)" }
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
$hashAfter = (Get-FileHash -Algorithm SHA256 $script).Hash
Write-Host ("scan finished: metas={0}/{1} rows={2}/{3}" -f $final.done, $nshards, $final.rows,
            ($nshards * $points_per_shard))
for ($s = 0; $s -lt $nshards; $s++) {
    $tag = '{0:D2}of{1:D2}' -f $s, $nshards
    $csv = Join-Path $out "carry1_shard$tag.csv"
    $n = if (Test-Path $csv) { (Import-Csv $csv).Count } else { -1 }
    Write-Host ("shard {0:D2} rows={1}" -f $s, $n)
}
$badErr = Get-ChildItem $out -Filter 'shard*.err.log' | Where-Object { $_.Length -gt 0 }
Write-Host ("non-empty stderr logs: {0}" -f $(if ($badErr) { ($badErr.Name -join ',') } else { 'none' }))
Write-Host "scanner sha256 after : $hashAfter"
Write-Host ("scanner unchanged    : {0}" -f ($hashBefore -eq $hashAfter))
Write-Host "carry1 scan done, elapsed $([int]((Get-Date) - $started).TotalMinutes) min"
