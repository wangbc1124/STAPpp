param(
    [string]$BuildDir = "build-mkl",
    [string]$Exe = "",
    [string]$ResultDir = "results",
    [string]$MklRoot = $(if ($env:STAPPP_MKL_ROOT) { $env:STAPPP_MKL_ROOT } else { $env:MKLROOT }),
    [string]$MklRuntimeDir = $env:STAPPP_MKL_RUNTIME_DIR,
    [string]$CompilerRuntimeDir = $env:STAPPP_COMPILER_RUNTIME_DIR
)

$ErrorActionPreference = "Stop"
$root = Resolve-Path -LiteralPath (Join-Path $PSScriptRoot "..")
Set-Location $root

if ([string]::IsNullOrWhiteSpace($MklRuntimeDir) -and -not [string]::IsNullOrWhiteSpace($MklRoot)) {
    $redist = Join-Path $MklRoot "redist\intel64"
    $bin = Join-Path $MklRoot "bin"
    if (Test-Path -LiteralPath $redist) {
        $MklRuntimeDir = $redist
    } elseif (Test-Path -LiteralPath $bin) {
        $MklRuntimeDir = $bin
    }
}
if (-not [string]::IsNullOrWhiteSpace($MklRuntimeDir) -and (Test-Path -LiteralPath $MklRuntimeDir)) {
    $env:PATH = "$MklRuntimeDir;$env:PATH"
}
if (-not [string]::IsNullOrWhiteSpace($CompilerRuntimeDir) -and (Test-Path -LiteralPath $CompilerRuntimeDir)) {
    $env:PATH = "$CompilerRuntimeDir;$env:PATH"
}

if ([string]::IsNullOrWhiteSpace($Exe)) {
    $candidates = @(
        (Join-Path $BuildDir "stap++.exe"),
        (Join-Path $BuildDir "Release\stap++.exe"),
        (Join-Path $BuildDir "Debug\stap++.exe")
    )
    foreach ($candidate in $candidates) {
        if (Test-Path -LiteralPath $candidate) {
            $Exe = $candidate
            break
        }
    }
}

if ([string]::IsNullOrWhiteSpace($Exe) -or -not (Test-Path -LiteralPath $Exe)) {
    throw "stap++.exe was not found. Run scripts\configure_mkl.ps1 first, or pass -Exe explicitly."
}

New-Item -ItemType Directory -Force -Path $ResultDir | Out-Null

$cases = @(
    @{ Name = "Bridge-1"; Input = "Bridge-1\Bridge-1.dat"; Csv = "$ResultDir\Bridge-1.displacements.csv"; Mtype = "auto" },
    @{ Name = "Bridge-2"; Input = "Bridge-2\Bridge-2.dat"; Csv = "$ResultDir\Bridge-2.displacements.csv"; Mtype = "sym-indef" },
    @{ Name = "Bridge-3"; Input = "Bridge-3\Bridge-3.dat"; Csv = "$ResultDir\Bridge-3.displacements.csv"; Mtype = "auto" },
    @{ Name = "Bridge-4"; Input = "Bridge-4.dat"; Csv = "$ResultDir\Bridge-4.displacements.csv"; Mtype = "auto" }
)

foreach ($case in $cases) {
    if (-not (Test-Path -LiteralPath $case.Input)) {
        throw "Input file not found: $($case.Input)"
    }

    Write-Host "=== Running $($case.Name) ==="
    $args = @(
        $case.Input,
        "--solver", "sparse-auto",
        "--backend", "pardiso",
        "--pardiso-mtype", $case.Mtype,
        "--output", "summary",
        "--tol", "1e-6",
        "--max-iter", "5000",
        "--csv", $case.Csv
    )
    Write-Host "+ $Exe $($args -join ' ')"
    & $Exe @args
    if ($LASTEXITCODE -ne 0) {
        exit $LASTEXITCODE
    }
}

Write-Host "All bridge cases completed."
& (Join-Path $PSScriptRoot "check_results.ps1") -ResultDir $ResultDir
