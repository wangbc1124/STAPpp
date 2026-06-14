# STAP++ Evaluation Clean Package

This directory is a minimal evaluation package for the finite element course project.

## Included

- `src/`: STAP++ source code.
- `CMakeLists.txt`: root CMake build entry.
- `Bridge-1/Bridge-1.dat`
- `Bridge-2/Bridge-2.dat`
- `Bridge-3/Bridge-3.dat`
- `Bridge-4.dat`
- `tools/inp2dat/`: optional `.inp` to `.dat` converter.
- `scripts/configure_mkl.ps1`: configure and build with Intel oneMKL PARDISO.
- `scripts/run_all_bridges.ps1`: run Bridge-1 to Bridge-4 in sequence.
- `scripts/check_results.ps1`: extract residuals, timings, and CSV status.
- `scripts/open_vscode_env.ps1`: map local ASCII drives and open VS Code for CMake Tools.
- `.vscode/`: VS Code CMake Tools settings driven by `STAPPP_*` environment variables.
- `docs/command_line_run_manual.md`: detailed Chinese command-line manual.
- `docs/environment_variables.md`: portable environment variable reference.

## Excluded

The package intentionally excludes `.git`, `.claude`, `tmp`, build directories, old executables, Doxygen HTML, previous `.out`, `.csv`, and `.vtk` generated files.

## Quick Start

Open PowerShell in this directory:

```powershell
$env:STAPPP_MKL_ROOT = "<your-mkl-root>"
$env:STAPPP_CMAKE_EXE = "<path-to-cmake.exe>"
$env:STAPPP_COMPILER_RUNTIME_DIR = "<path-to-mingw-bin>"
$env:STAPPP_C_COMPILER = "$env:STAPPP_COMPILER_RUNTIME_DIR\gcc.exe"
$env:STAPPP_CXX_COMPILER = "$env:STAPPP_COMPILER_RUNTIME_DIR\g++.exe"
.\scripts\configure_mkl.ps1
.\scripts\run_all_bridges.ps1
```

If using Visual Studio generator:

```powershell
.\scripts\configure_mkl.ps1 -UseVisualStudio
.\scripts\run_all_bridges.ps1 -BuildDir build-mkl
```

Check results:

```powershell
.\scripts\check_results.ps1
```

## VS Code CMake Tools

Run from PowerShell:

```powershell
.\scripts\open_vscode_env.ps1
```

Then in VS Code:

1. Install `CMake Tools` and `C/C++`.
2. Run `CMake: Select a Kit`.
3. Choose `STAP++ MinGW from environment`.
4. Run `CMake: Configure`.
5. Run `CMake: Build`.
