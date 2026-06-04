$ErrorActionPreference = "Stop"
Set-Location "D:\STAPpp-master\STAPpp-master\Bridge-1"

$exe = "D:\STAPpp-master\STAPpp-master\stap++.exe"
$dat = "D:\STAPpp-master\STAPpp-master\Bridge-1\Bridge-1.dat"
$out = "D:\STAPpp-master\STAPpp-master\Bridge-1\Bridge-1-test8.out"

Write-Host "Running: $exe $dat"

try {
    $result = & $exe $dat 2>&1
    $result | Out-File -FilePath $out -Encoding UTF8
    Write-Host "Exit code: $LASTEXITCODE"
    Write-Host "Output lines: $($result.Count)"
} catch {
    Write-Host "Exception: $_"
    $_.Exception.Message | Out-File -FilePath $out -Encoding UTF8
}
