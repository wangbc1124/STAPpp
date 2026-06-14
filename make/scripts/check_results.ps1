param(
    [string]$ResultDir = "results"
)

$ErrorActionPreference = "Stop"
$root = Resolve-Path -LiteralPath (Join-Path $PSScriptRoot "..")
Set-Location $root

$outFiles = @(
    "Bridge-1\Bridge-1.out",
    "Bridge-2\Bridge-2.out",
    "Bridge-3\Bridge-3.out",
    "Bridge-4.out"
)

Write-Host "=== Solver summaries ==="
foreach ($outFile in $outFiles) {
    if (-not (Test-Path -LiteralPath $outFile)) {
        Write-Host "Missing output: $outFile"
        continue
    }
    Write-Host "--- $outFile ---"
    Select-String -LiteralPath $outFile -Pattern "NUMBER OF EQUATIONS|NUMBER OF CSR NONZEROS|ACTUAL SOLVER|ITERATIVE SOLVER CONVERGED|CHECKED RELATIVE RESIDUAL|PARDISO_FACT_NNZ|TOTAL_TIME|DISPLACEMENT CSV"
}

Write-Host "=== CSV files ==="
$csvFiles = @(
    "$ResultDir\Bridge-1.displacements.csv",
    "$ResultDir\Bridge-2.displacements.csv",
    "$ResultDir\Bridge-3.displacements.csv",
    "$ResultDir\Bridge-4.displacements.csv"
)
foreach ($csvFile in $csvFiles) {
    if (Test-Path -LiteralPath $csvFile) {
        Get-Item -LiteralPath $csvFile | Select-Object Length,LastWriteTime,FullName
    } else {
        Write-Host "Missing CSV: $csvFile"
    }
}
