param(
    [string]$ProjectDrive = "X:",
    [string]$MklDrive = "Y:",
    [string]$ProjectPath = "",
    [string]$MklPath = "",
    [string]$CMakeExe = "",
    [string]$CompilerBin = "",
    [string]$CCompiler = "",
    [string]$CxxCompiler = ""
)

$ErrorActionPreference = "Stop"

function Find-CommandPath {
    param([string]$Name)
    $command = Get-Command $Name -ErrorAction SilentlyContinue
    if ($command) {
        return $command.Source
    }
    return ""
}

function Find-FirstExistingPath {
    param([string[]]$Candidates)
    foreach ($candidate in $Candidates) {
        if (-not [string]::IsNullOrWhiteSpace($candidate) -and (Test-Path -LiteralPath $candidate)) {
            return (Resolve-Path -LiteralPath $candidate).Path
        }
    }
    return ""
}

function ConvertTo-CMakePath {
    param([string]$Path)
    if ([string]::IsNullOrWhiteSpace($Path)) {
        return ""
    }
    return $Path.Replace("\", "/")
}

if ([string]::IsNullOrWhiteSpace($ProjectPath)) {
    $ProjectPath = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot "..")).Path
}

if ([string]::IsNullOrWhiteSpace($MklPath)) {
    if (-not [string]::IsNullOrWhiteSpace($env:MKLROOT) -and (Test-Path -LiteralPath $env:MKLROOT)) {
        $MklPath = (Resolve-Path -LiteralPath $env:MKLROOT).Path
    } elseif (-not [string]::IsNullOrWhiteSpace($env:STAPPP_MKL_ROOT) -and (Test-Path -LiteralPath $env:STAPPP_MKL_ROOT)) {
        $MklPath = (Resolve-Path -LiteralPath $env:STAPPP_MKL_ROOT).Path
    }
}

if (-not (Test-Path -LiteralPath $ProjectPath)) {
    throw "Project path not found: $ProjectPath"
}
if ([string]::IsNullOrWhiteSpace($MklPath) -or -not (Test-Path -LiteralPath $MklPath)) {
    throw "MKL path not found. Set STAPPP_MKL_ROOT, set MKLROOT, or pass -MklPath."
}

subst $ProjectDrive $ProjectPath
subst $MklDrive $MklPath

if ([string]::IsNullOrWhiteSpace($CMakeExe)) {
    if (-not [string]::IsNullOrWhiteSpace($env:STAPPP_CMAKE_EXE)) {
        $CMakeExe = $env:STAPPP_CMAKE_EXE
    } else {
        $CMakeExe = Find-CommandPath "cmake.exe"
    }
}
if ([string]::IsNullOrWhiteSpace($CMakeExe)) {
    $CMakeExe = Find-FirstExistingPath @(
        "C:\Program Files\CMake\bin\cmake.exe",
        "C:\Program Files (x86)\CMake\bin\cmake.exe",
        "C:\Program Files\Microsoft Visual Studio\2022\Community\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe",
        "C:\Program Files\Microsoft Visual Studio\2022\Professional\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe",
        "C:\Program Files\Microsoft Visual Studio\2022\Enterprise\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe",
        "C:\Program Files\Microsoft Visual Studio\2022\BuildTools\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe",
        "C:\Program Files (x86)\Microsoft Visual Studio\2019\BuildTools\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe"
    )
}

if ([string]::IsNullOrWhiteSpace($CompilerBin)) {
    if (-not [string]::IsNullOrWhiteSpace($env:STAPPP_COMPILER_RUNTIME_DIR)) {
        $CompilerBin = $env:STAPPP_COMPILER_RUNTIME_DIR
    } elseif (-not [string]::IsNullOrWhiteSpace($env:STAPPP_MINGW_BIN)) {
        $CompilerBin = $env:STAPPP_MINGW_BIN
    }
}
if ([string]::IsNullOrWhiteSpace($CompilerBin)) {
    $gxx = Find-CommandPath "g++.exe"
    if ($gxx) {
        $CompilerBin = Split-Path -Parent $gxx
    }
}

if ([string]::IsNullOrWhiteSpace($CMakeExe) -or -not (Test-Path -LiteralPath $CMakeExe)) {
    throw "cmake.exe not found. Set STAPPP_CMAKE_EXE, add cmake to PATH, or pass -CMakeExe."
}
if ([string]::IsNullOrWhiteSpace($CompilerBin) -or -not (Test-Path -LiteralPath $CompilerBin)) {
    throw "Compiler bin directory not found. Set STAPPP_COMPILER_RUNTIME_DIR or STAPPP_MINGW_BIN, add g++ to PATH, or pass -CompilerBin."
}

$cmakeBin = Split-Path -Parent (Resolve-Path -LiteralPath $CMakeExe).Path
$compilerBinResolved = (Resolve-Path -LiteralPath $CompilerBin).Path

if ([string]::IsNullOrWhiteSpace($CCompiler)) {
    if (-not [string]::IsNullOrWhiteSpace($env:STAPPP_C_COMPILER)) {
        $CCompiler = $env:STAPPP_C_COMPILER
    } else {
        $CCompiler = Join-Path $compilerBinResolved "gcc.exe"
    }
}
if ([string]::IsNullOrWhiteSpace($CxxCompiler)) {
    if (-not [string]::IsNullOrWhiteSpace($env:STAPPP_CXX_COMPILER)) {
        $CxxCompiler = $env:STAPPP_CXX_COMPILER
    } else {
        $CxxCompiler = Join-Path $compilerBinResolved "g++.exe"
    }
}
if (-not (Test-Path -LiteralPath $CxxCompiler)) {
    throw "C++ compiler not found: $CxxCompiler"
}

$env:MKLROOT = "$MklDrive\"
$env:STAPPP_MKL_ROOT = "$MklDrive\"
$env:STAPPP_MKL_INCLUDE_DIR = ConvertTo-CMakePath "$MklDrive\include"
$env:STAPPP_MKL_LIBRARY = ConvertTo-CMakePath "$MklDrive\lib\mkl_rt.lib"
$env:STAPPP_MKL_RUNTIME_DIR = ConvertTo-CMakePath "$MklDrive\bin"
$env:STAPPP_COMPILER_RUNTIME_DIR = ConvertTo-CMakePath $compilerBinResolved
$env:STAPPP_C_COMPILER = ConvertTo-CMakePath $CCompiler
$env:STAPPP_CXX_COMPILER = ConvertTo-CMakePath $CxxCompiler
$env:STAPPP_CMAKE_EXE = ConvertTo-CMakePath (Resolve-Path -LiteralPath $CMakeExe).Path
$env:PATH = "$cmakeBin;$compilerBinResolved;$MklDrive\bin;$env:PATH"

Write-Host "Project mapped to $ProjectDrive\"
Write-Host "MKL mapped to $MklDrive\"
Write-Host "CMake: $env:STAPPP_CMAKE_EXE"
Write-Host "C compiler: $env:STAPPP_C_COMPILER"
Write-Host "C++ compiler: $env:STAPPP_CXX_COMPILER"
Write-Host "Opening VS Code..."

code "$ProjectDrive\"
