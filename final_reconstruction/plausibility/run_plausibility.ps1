# Driver for the parameter-plausibility package.
#
#   powershell -ExecutionPolicy Bypass -File .\run_plausibility.ps1 checks
#   powershell -ExecutionPolicy Bypass -File .\run_plausibility.ps1 scan -Nshards 12 -Arm AB
#   powershell -ExecutionPolicy Bypass -File .\run_plausibility.ps1 merge -Nshards 12
#   powershell -ExecutionPolicy Bypass -File .\run_plausibility.ps1 manifest
#
# NOTE ON THIS SANDBOX: Wait-Process on child handles is denied here, so shard
# completion is detected from the filesystem (each shard writes its CSV after
# every point and its meta JSON as its last act).  The driver therefore stays
# alive for the whole scan, which is also what keeps the children alive.
param(
    [Parameter(Position = 0)][ValidateSet('checks', 'scan', 'merge', 'manifest')][string]$Mode = 'checks',
    [int]$Nshards = 12,
    [ValidateSet('A', 'B', 'AB')][string]$Arm = 'AB',
    [double]$Hours = 600.0,
    [double]$SampleMin = 2.0,
    [double]$MaxStepMin = 2.0
)

$ErrorActionPreference = 'Continue'
$env:MPLBACKEND = 'Agg'
$env:PYTHONIOENCODING = 'utf-8'
$env:OMP_NUM_THREADS = '1'
$env:OPENBLAS_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'

$here = Split-Path -Parent $MyInvocation.MyCommand.Path
$root = Split-Path -Parent $here
$py = 'D:\aconade\python.exe'
$started = Get-Date

if ($Mode -eq 'checks') {
    # Five independent processes; completion is detected from their output JSON
    # files (Wait-Process is denied in this sandbox).
    $checkScripts = @('check_threshold_margin.py', 'check_absolute_scale.py',
                      'check_mass_balance.py', 'check_identifiability.py',
                      'check_rdf_collapse.py')
    $sentinels = @('threshold_margin\threshold_margin_baseline.json', 'absolute_scale.json',
                   'mass_balance.json', 'identifiability.json', 'rdf_collapse.json')
    $procs = @()
    for ($i = 0; $i -lt $checkScripts.Count; $i++) {
        $name = [IO.Path]::GetFileNameWithoutExtension($checkScripts[$i])
        $procs += Start-Process -FilePath $py -WorkingDirectory $here -NoNewWindow -PassThru `
            -ArgumentList @($checkScripts[$i]) `
            -RedirectStandardOutput (Join-Path $here "log_$name.txt") `
            -RedirectStandardError (Join-Path $here "log_$name.err.txt")
        Start-Sleep -Milliseconds 500
    }
    Write-Host "launched $($procs.Count) check processes"
    $deadline = (Get-Date).AddHours(6)
    $last = -1
    while ($true) {
        Start-Sleep -Seconds 30
        $done = 0
        foreach ($s in $sentinels) { if (Test-Path (Join-Path $here $s)) { $done++ } }
        if ($done -ne $last) {
            Write-Host ("t+{0,4} min  checks finished={1}/{2}" -f
                        [int]((Get-Date) - $started).TotalMinutes, $done, $sentinels.Count)
            $last = $done
        }
        if ($done -ge $sentinels.Count) { break }
        if ((Get-Date) -gt $deadline) { Write-Host 'DEADLINE reached'; break }
    }
    $bad = Get-ChildItem $here -Filter 'log_*.err.txt' | Where-Object { $_.Length -gt 0 }
    Write-Host ("non-empty stderr: {0}" -f $(if ($bad) { ($bad.Name -join ',') } else { 'none' }))
    & $py -c "import sys; sys.path.insert(0,r'$here'); import plausibility_common as C; print(len(C.write_manifest()), 'manifest entries')"
    Write-Host "checks done, elapsed $([int]((Get-Date) - $started).TotalMinutes) min"
    exit 0
}

if ($Mode -eq 'manifest') {
    & $py -c "import sys; sys.path.insert(0,r'$here'); import plausibility_common as C; print(len(C.write_manifest()), 'entries')"
    exit $LASTEXITCODE
}

if ($Mode -eq 'merge') {
    & $py (Join-Path $here 'scan_sharpen_vs_lengthen.py') merge --nshards $Nshards --projection-guard
    exit $LASTEXITCODE
}

# ------------------------------------------------------------------- scan
Write-Host "plausibility scan start $($started.ToString('s')) arm=$Arm shards=$Nshards hours=$Hours"
$procs = @()
for ($s = 0; $s -lt $Nshards; $s++) {
    $tag = '{0:d2}' -f $s
    $p = Start-Process -FilePath $py -WorkingDirectory $here -NoNewWindow -PassThru `
        -ArgumentList @('scan_sharpen_vs_lengthen.py', 'scan', '--arm', $Arm,
                        '--shard', $s, '--nshards', $Nshards, '--hours', $Hours,
                        '--sample-min', $SampleMin, '--max-step-min', $MaxStepMin) `
        -RedirectStandardOutput (Join-Path $here "scan$tag.out.log") `
        -RedirectStandardError (Join-Path $here "scan$tag.err.log")
    $procs += [pscustomobject]@{ shard = $s; pid = $p.Id }
    Start-Sleep -Milliseconds 400
}
Write-Host ("launched {0} shard processes" -f $procs.Count)

$deadline = (Get-Date).AddHours(12)
$last = -1
while ($true) {
    Start-Sleep -Seconds 60
    $done = 0; $rows = 0
    for ($s = 0; $s -lt $Nshards; $s++) {
        $tag = '{0:d2}of{1:d2}' -f $s, $Nshards
        $meta = Join-Path $here "meta_sharpen_vs_lengthen_$tag.json"
        $csv = Join-Path $here "sharpen_vs_lengthen_$tag.csv"
        if (Test-Path $meta) { $done++ }
        if (Test-Path $csv) { $rows += (Import-Csv $csv).Count }
    }
    if ($rows -ne $last -or $done -ge $Nshards) {
        Write-Host ("t+{0,4} min  metas={1}/{2}  rows={3}" -f [int]((Get-Date) - $started).TotalMinutes,
                    $done, $Nshards, $rows)
        $last = $rows
    }
    if ($done -ge $Nshards) { break }
    if ((Get-Date) -gt $deadline) { Write-Host 'DEADLINE reached'; break }
}
$bad = Get-ChildItem $here -Filter 'scan*.err.log' | Where-Object { $_.Length -gt 0 }
Write-Host ("non-empty stderr: {0}" -f $(if ($bad) { ($bad.Name -join ',') } else { 'none' }))
Write-Host "plausibility scan done, elapsed $([int]((Get-Date) - $started).TotalMinutes) min"
