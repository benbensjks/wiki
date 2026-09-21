# dsh driver: stage-3 split-maturation scan (63 points, 600 h, 9 shards).
# Children are launched with Start-Process; completion is detected from the
# filesystem (each shard writes its CSV after every point and a meta JSON last),
# so the driver stays alive for the whole run and the children are not reaped.
#
#   powershell -ExecutionPolicy Bypass -File .\run_carry1_split_dsh.ps1 scan
#   powershell -ExecutionPolicy Bypass -File .\run_carry1_split_dsh.ps1 merge
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
$script = Join-Path $root 'scan_carry1_split_dsh.py'
$out = Join-Path $root 'threebit51_results\carry1_split_scan'
$nshards = 9
$points_per_shard = 7
$hours = 600
$deadline = (Get-Date).AddHours(4)
$started = Get-Date

New-Item -ItemType Directory -Force -Path $out | Out-Null

if ($Mode -eq 'merge') {
    & $py $script merge --nshards $nshards --hours $hours
    Write-Host "merge exit=$LASTEXITCODE"
    exit $LASTEXITCODE
}

$hashBefore = (Get-FileHash -Algorithm SHA256 $script).Hash
Write-Host "split scan start $($started.ToString('s')) shards=$nshards hours=$hours"
Write-Host "scanner sha256 before: $hashBefore"

0..($nshards - 1) | ForEach-Object {
    $s = $_
    Start-Process -FilePath $py -WorkingDirectory $root -NoNewWindow -PassThru `
        -ArgumentList @('scan_carry1_split_dsh.py', 'scan', '--shard', $s,
                        '--nshards', $nshards, '--hours', $hours) `
        -RedirectStandardOutput (Join-Path $out ("shard{0:D2}.out.log" -f $s)) `
        -RedirectStandardError (Join-Path $out ("shard{0:D2}.err.log" -f $s)) | Out-Null
    Start-Sleep -Milliseconds 400
}
Write-Host "launched $nshards shard processes"

function Get-State {
    $rows = 0; $done = 0
    for ($s = 0; $s -lt $nshards; $s++) {
        $tag = '{0:D2}of{1:D2}' -f $s, $nshards
        $csv = Join-Path $out "split_shard$tag.csv"
        if (Test-Path $csv) { $rows += (Import-Csv $csv).Count }
        if (Test-Path (Join-Path $out "meta_shard$tag.json")) { $done++ }
    }
    return [pscustomobject]@{ rows = $rows; done = $done }
}

$last = -1
while ($true) {
    Start-Sleep -Seconds 60
    $st = Get-State
    if ($st.rows -ne $last -or $st.done -ge $nshards) {
        Write-Host ("t+{0,3} min  metas={1}/{2}  rows={3}/{4}" -f
                    [int]((Get-Date) - $started).TotalMinutes, $st.done, $nshards, $st.rows,
                    ($nshards * $points_per_shard))
        $last = $st.rows
    }
    if ($st.done -ge $nshards) { break }
    if ((Get-Date) -gt $deadline) { Write-Host 'DEADLINE reached'; break }
}

$hashAfter = (Get-FileHash -Algorithm SHA256 $script).Hash
Write-Host ("scan finished: metas={0}/{1} rows={2}/{3}" -f $st.done, $nshards, $st.rows,
            ($nshards * $points_per_shard))
$badErr = Get-ChildItem $out -Filter 'shard*.err.log' | Where-Object { $_.Length -gt 0 }
Write-Host ("non-empty stderr logs: {0}" -f $(if ($badErr) { ($badErr.Name -join ',') } else { 'none' }))
Write-Host "scanner sha256 after : $hashAfter"
Write-Host ("scanner unchanged    : {0}" -f ($hashBefore -eq $hashAfter))
Write-Host "elapsed $([int]((Get-Date) - $started).TotalMinutes) min"
