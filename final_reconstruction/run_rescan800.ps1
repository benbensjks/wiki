# Serial-equivalent sharded re-scan of the original 800-point Extension grid.
# Each shard is an independent OS process launched with Start-Process: no
# multiprocessing pipes are used, so this works under the current file sandbox
# and on Windows PowerShell 5.1 (which has no ForEach-Object -Parallel).
#
#   powershell -ExecutionPolicy Bypass -File .\run_rescan800.ps1
#
$ErrorActionPreference = 'Continue'
$env:MPLBACKEND = 'Agg'
$env:PYTHONIOENCODING = 'utf-8'
$env:OMP_NUM_THREADS = '1'
$env:OPENBLAS_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'

$here = Split-Path -Parent $MyInvocation.MyCommand.Path
$python = 'D:\aconade\python.exe'
$script = Join-Path $here 'rescan800_readwindow.py'
$out = Join-Path $here 'bit0_results\rescan800'
$nshards = 16
$hours = 300

New-Item -ItemType Directory -Force -Path $out | Out-Null
$started = Get-Date
Write-Host "rescan800 start $($started.ToString('s')) shards=$nshards hours=$hours"

$procs = @()
for ($s = 0; $s -lt $nshards; $s++) {
    $tag = '{0:d2}' -f $s
    $a = @($script, '--shard', $s, '--nshards', $nshards, '--hours', $hours, '--out-dir', $out)
    $p = Start-Process -FilePath $python -ArgumentList $a -NoNewWindow -PassThru `
        -RedirectStandardOutput (Join-Path $out "shard$tag.out.log") `
        -RedirectStandardError  (Join-Path $out "shard$tag.err.log")
    $procs += [pscustomobject]@{ shard = $s; proc = $p }
}
Write-Host "launched $($procs.Count) shard processes"

foreach ($entry in $procs) {
    $entry.proc.WaitForExit()
    $entry.proc.Refresh()
    Write-Host ("shard {0:d2} exit={1}" -f $entry.shard, $entry.proc.ExitCode)
}

$rows = 0
$missing = @()
for ($s = 0; $s -lt $nshards; $s++) {
    $tag = '{0:d2}' -f $s
    $csv = Join-Path $out "rescan800_shard${tag}of${nshards}.csv"
    if (Test-Path $csv) {
        $n = (Import-Csv $csv).Count
        $rows += $n
        Write-Host ("shard {0:d2} rows={1}" -f $s, $n)
    } else { $missing += $tag }
}
Write-Host "rescan800 total rows=$rows expected=$((800 / $nshards) * $nshards) missing=$($missing -join ',')"
Write-Host "rescan800 done, elapsed $([int]((Get-Date) - $started).TotalMinutes) min"
