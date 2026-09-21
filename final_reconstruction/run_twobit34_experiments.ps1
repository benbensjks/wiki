# dsh driver for the frozen 34-state perturbation and strict-tolerance runs.
# Follows twobit34_results/DSH扰动与严格容差说明.md.
#
# This sandbox denies Wait-Process on child handles, so completion is detected
# from the filesystem: a shard is done when its CSV holds the expected number of
# rows and its meta JSON exists (the meta is the shard's last act).  The driver
# therefore stays alive for the whole run, which also keeps the children alive.
#
#   powershell -ExecutionPolicy Bypass -File .\run_twobit34_experiments.ps1 -Kind S -Mode scan
#   powershell -ExecutionPolicy Bypass -File .\run_twobit34_experiments.ps1 -Kind S -Mode merge
#
param(
    [ValidateSet('S', 'RDF', 'AFFL', 'INT1', 'STRICT')][string]$Kind,
    [ValidateSet('scan', 'merge')][string]$Mode = 'scan'
)

$ErrorActionPreference = 'Continue'
$env:MPLBACKEND = 'Agg'
$env:PYTHONIOENCODING = 'utf-8'
$env:OMP_NUM_THREADS = '1'
$env:OPENBLAS_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'

$root = 'C:\Users\18633\Desktop\wiki\final_reconstruction'
$py = 'D:\aconade\python.exe'
$res = Join-Path $root 'twobit34_results'

$cfg = @{
    S      = @{ script = 'scan_twobit34_perturbations.py';   out = 'perturbations';    nshards = 10; hours = 200; per = 4;
                csv = 'perturb_s_shard{0:D2}of{1:D2}.csv';    meta = 'meta_s_shard{0:D2}of{1:D2}.json';    merge = @('merge', '--group', 'S', '--nshards', '10') }
    RDF    = @{ script = 'scan_twobit34_perturbations.py';   out = 'perturbations';    nshards = 10; hours = 200; per = 4;
                csv = 'perturb_rdf_shard{0:D2}of{1:D2}.csv';  meta = 'meta_rdf_shard{0:D2}of{1:D2}.json';  merge = @('merge', '--group', 'RDF', '--nshards', '10') }
    AFFL   = @{ script = 'scan_twobit34_perturbations.py';   out = 'perturbations';    nshards = 10; hours = 200; per = 4;
                csv = 'perturb_affl_shard{0:D2}of{1:D2}.csv'; meta = 'meta_affl_shard{0:D2}of{1:D2}.json'; merge = @('merge', '--group', 'AFFL', '--nshards', '10') }
    INT1   = @{ script = 'scan_twobit34_perturbations.py';   out = 'perturbations';    nshards = 5;  hours = 200; per = 4;
                csv = 'perturb_int1_shard{0:D2}of{1:D2}.csv'; meta = 'meta_int1_shard{0:D2}of{1:D2}.json'; merge = @('merge', '--group', 'INT1', '--nshards', '5') }
    STRICT = @{ script = 'scan_twobit34_strict_tolerance.py'; out = 'strict_tolerance'; nshards = 9; hours = 300; per = 3;
                csv = 'strict_shard{0:D2}of{1:D2}.csv';       meta = 'meta_shard{0:D2}of{1:D2}.json';      merge = @('merge', '--nshards', '9') }
}
$c = $cfg[$Kind]
$out = Join-Path $res $c.out
$scriptPath = Join-Path $root $c.script
$n = [int]$c.nshards
$per = [int]$c.per
$hours = [double]$c.hours
$started = Get-Date

New-Item -ItemType Directory -Force -Path $out | Out-Null

if ($Mode -eq 'merge') {
    $argv = @($scriptPath) + $c.merge
    & $py @argv
    Write-Host "$Kind merge exit=$LASTEXITCODE"
    exit $LASTEXITCODE
}

function Shard-Csv($s) { Join-Path $out ($c.csv -f $s, $n) }
function Shard-Meta($s) { Join-Path $out ($c.meta -f $s, $n) }

function Get-State {
    $rows = 0; $done = 0; $stalled = @()
    for ($s = 0; $s -lt $n; $s++) {
        $csv = Shard-Csv $s
        if (Test-Path $csv) {
            $k = (Import-Csv $csv).Count
            $rows += $k
            $age = ((Get-Date) - (Get-Item $csv).LastWriteTime).TotalMinutes
            if ($k -lt $per -and $age -gt 25) { $stalled += "s$s(rows=$k,idle=$([int]$age)m)" }
        }
        if (Test-Path (Shard-Meta $s)) { $done++ }
    }
    [pscustomobject]@{ rows = $rows; done = $done; stalled = $stalled }
}

Write-Host "$Kind scan start $($started.ToString('s')) shards=$n hours=$hours per_shard=$per"
$launched = @()
for ($s = 0; $s -lt $n; $s++) {
    $stdout = Join-Path $out ("{0}_shard{1:D2}.out.log" -f $Kind.ToLower(), $s)
    $stderr = Join-Path $out ("{0}_shard{1:D2}.err.log" -f $Kind.ToLower(), $s)
    $argv = @($c.script, 'scan', '--shard', $s, '--nshards', $n, '--hours', $hours)
    if ($Kind -ne 'STRICT') { $argv += @('--group', $Kind) }
    $p = Start-Process -FilePath $py -WorkingDirectory $root -NoNewWindow -PassThru `
        -ArgumentList $argv -RedirectStandardOutput $stdout -RedirectStandardError $stderr
    $launched += $p.Id
    Start-Sleep -Milliseconds 400
}
Write-Host ("launched $($launched.Count) processes: " + ($launched -join ','))

$lastRows = -1
$deadline = $started.AddHours(6)
while ($true) {
    Start-Sleep -Seconds 60
    $st = Get-State
    $elapsed = [int]((Get-Date) - $started).TotalMinutes
    if ($st.rows -ne $lastRows -or $st.done -ge $n) {
        Write-Host ("t+{0,3} min metas={1}/{2} rows={3}/{4}{5}" -f $elapsed, $st.done, $n, $st.rows,
                    ($n * $per), $(if ($st.stalled.Count) { "  STALLED: " + ($st.stalled -join ',') } else { '' }))
        $lastRows = $st.rows
    }
    if ($st.done -ge $n) { break }
    if ((Get-Date) -gt $deadline) { Write-Host 'DEADLINE reached'; break }
}

$final = Get-State
Write-Host ("$Kind scan finished: metas={0}/{1} rows={2}/{3}" -f $final.done, $n, $final.rows, ($n * $per))
for ($s = 0; $s -lt $n; $s++) {
    $csv = Shard-Csv $s
    $k = if (Test-Path $csv) { (Import-Csv $csv).Count } else { -1 }
    Write-Host ("  shard {0:D2} rows={1}" -f $s, $k)
}
$bad = Get-ChildItem $out -Filter ("{0}_shard*.err.log" -f $Kind.ToLower()) | Where-Object { $_.Length -gt 0 }
Write-Host ("non-empty stderr: " + $(if ($bad) { ($bad.Name -join ',') } else { 'none' }))
Write-Host "$Kind scan done, elapsed $([int]((Get-Date) - $started).TotalMinutes) min"
