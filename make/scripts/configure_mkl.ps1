param(
    [string]$BuildDir = "build-mkl",
    [string]$Generator = "",
    [string]$CMakeExe = $env:STAPPP_CMAKE_EXE,
    [string]$MklRoot = $(if ($env:STAPPP_MKL_ROOT) { $env:STAPPP_MKL_ROOT } else { $env:MKLROOT }),
    [string]$MklIncludeDir = $env:STAPPP_MKL_INCLUDE_DIR,
    [string]$MklLibrary = $env:STAPPP_MKL_LIBRARY,
    [string]$MklRuntimeDir = $env:STAPPP_MKL_RUNTIME_DIR,
    [string]$CompilerRuntimeDir = $env:STAPPP_COMPILER_RUNTIME_DIR,
    [string]$ExtraRuntimeDirs = $env:STAPPP_EXTRA_RUNTIME_DIRS,
    [switch]$UseVisualStudio
)

$ErrorActionPreference = "Stop"
$root = Resolve-Path -LiteralPath (Join-Path $PSScriptRoot "..")
Set-Location $root

function ConvertTo-CMakePath {
    param([string]$Path)
    if ([string]::IsNullOrWhiteSpace($Path)) {
        return ""
    }
    return $Path.Replace("\", "/")
}

if ([string]::IsNullOrWhiteSpace($CMakeExe)) {
    $cmakeCommand = Get-Command cmake -ErrorAction SilentlyContinue
    if ($cmakeCommand) {
        $CMakeExe = $cmakeCommand.Source
    }
}

if ([string]::IsNullOrWhiteSpace($CMakeExe)) {
    $cmakeCandidates = @(
        "C:\Program Files\CMake\bin\cmake.exe",
        "C:\Program Files (x86)\CMake\bin\cmake.exe",
        "C:\Program Files\Microsoft Visual Studio\2022\Community\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe",
        "C:\Program Files\Microsoft Visual Studio\2022\Professional\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe",
        "C:\Program Files\Microsoft Visual Studio\2022\Enterprise\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe",
        "C:\Program Files\Microsoft Visual Studio\2022\BuildTools\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe",
        "C:\Program Files (x86)\Microsoft Visual Studio\2019\BuildTools\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe"
    )
    foreach ($candidate in $cmakeCandidates) {
        if (Test-Path -LiteralPath $candidate) {
            $CMakeExe = $candidate
            break
        }
    }
}

if ([string]::IsNullOrWhiteSpace($CMakeExe) -or -not (Test-Path -LiteralPath $CMakeExe)) {
    throw "cmake.exe was not found. Add CMake to PATH, install CMake, or pass -CMakeExe explicitly."
}

if ([string]::IsNullOrWhiteSpace($MklRoot) -and
    ([string]::IsNullOrWhiteSpace($MklIncludeDir) -or [string]::IsNullOrWhiteSpace($MklLibrary))) {
    $defaultMkl = "C:\Program Files (x86)\Intel\oneAPI\mkl\latest"
    if (Test-Path -LiteralPath $defaultMkl) {
        $MklRoot = $defaultMkl
    }
}

if ([string]::IsNullOrWhiteSpace($MklIncludeDir)) {
    $MklIncludeDir = Join-Path $MklRoot "include"
}
if ([string]::IsNullOrWhiteSpace($MklLibrary)) {
    $MklLibrary = Join-Path $MklRoot "lib\mkl_rt.lib"
}
if ([string]::IsNullOrWhiteSpace($MklRuntimeDir)) {
    $redistCandidate = Join-Path $MklRoot "redist\intel64"
    $binCandidate = Join-Path $MklRoot "bin"
    if (Test-Path -LiteralPath $redistCandidate) {
        $MklRuntimeDir = $redistCandidate
    } elseif (Test-Path -LiteralPath $binCandidate) {
        $MklRuntimeDir = $binCandidate
    } else {
        $MklRuntimeDir = Split-Path -Parent $MklLibrary
    }
}
if ([string]::IsNullOrWhiteSpace($CompilerRuntimeDir)) {
    $gxxCommand = Get-Command g++ -ErrorAction SilentlyContinue
    if ($gxxCommand) {
        $CompilerRuntimeDir = Split-Path -Parent $gxxCommand.Source
    }
}

if (-not (Test-Path -LiteralPath $MklIncludeDir)) {
    throw "MKL include directory not found: $MklIncludeDir"
}
if (-not (Test-Path -LiteralPath $MklLibrary)) {
    throw "MKL library not found: $MklLibrary"
}
if (-not (Test-Path -LiteralPath $MklRuntimeDir)) {
    throw "MKL runtime directory not found: $MklRuntimeDir"
}

if (Test-Path -LiteralPath $MklRuntimeDir) {
    $env:PATH = "$MklRuntimeDir;$env:PATH"
}
if (-not [string]::IsNullOrWhiteSpace($CompilerRuntimeDir) -and (Test-Path -LiteralPath $CompilerRuntimeDir)) {
    $env:PATH = "$CompilerRuntimeDir;$env:PATH"
}

$cmakeArgs = @(
    "-S", ".",
    "-B", $BuildDir,
    "-DSTAPPP_ENABLE_MKL_PARDISO=ON",
    "-DSTAPPP_MKL_ROOT=$(ConvertTo-CMakePath $MklRoot)",
    "-DSTAPPP_MKL_INCLUDE_DIR=$(ConvertTo-CMakePath $MklIncludeDir)",
    "-DSTAPPP_MKL_LIBRARY=$(ConvertTo-CMakePath $MklLibrary)",
    "-DSTAPPP_MKL_RUNTIME_DIR=$(ConvertTo-CMakePath $MklRuntimeDir)"
)
if (-not [string]::IsNullOrWhiteSpace($CompilerRuntimeDir)) {
    $cmakeArgs += "-DSTAPPP_COMPILER_RUNTIME_DIR=$(ConvertTo-CMakePath $CompilerRuntimeDir)"
}
if (-not [string]::IsNullOrWhiteSpace($ExtraRuntimeDirs)) {
    $cmakeArgs += "-DSTAPPP_EXTRA_RUNTIME_DIRS=$(ConvertTo-CMakePath $ExtraRuntimeDirs)"
}

if ($UseVisualStudio) {
    $cmakeArgs += @("-G", "Visual Studio 17 2022", "-A", "x64")
} elseif (-not [string]::IsNullOrWhiteSpace($Generator)) {
    $cmakeArgs += @("-G", $Generator)
} else {
    $cmakeArgs += @("-DCMAKE_BUILD_TYPE=Release")
}

Write-Host "+ $CMakeExe $($cmakeArgs -join ' ')"
& $CMakeExe @cmakeArgs
if ($LASTEXITCODE -ne 0) {
    exit $LASTEXITCODE
}

$buildArgs = @("--build", $BuildDir, "--config", "Release")
Write-Host "+ $CMakeExe $($buildArgs -join ' ')"
& $CMakeExe @buildArgs
if ($LASTEXITCODE -ne 0) {
    exit $LASTEXITCODE
}

Write-Host "Build completed."
Write-Host "Candidate executables:"
Get-ChildItem -LiteralPath $BuildDir -Recurse -Filter "stap++.exe" | Select-Object FullName,Length,LastWriteTime
