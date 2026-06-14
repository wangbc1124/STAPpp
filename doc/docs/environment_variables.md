# STAP++ 可迁移环境变量清单

生成时间：2026-06-13

本文列出测评环境、命令行脚本、VS Code CMake Tools 和 CMake 运行目标可能需要设置或覆盖的全部环境变量。迁移到新机器时，优先改这些变量或脚本参数，不需要修改源码。

## 1. 推荐最小配置

如果使用 MinGW + MKL/PARDISO，建议先在 PowerShell 中设置：

```powershell
$env:STAPPP_MKL_ROOT = "<your-mkl-root>"
$env:STAPPP_CMAKE_EXE = "<path-to-cmake.exe>"
$env:STAPPP_COMPILER_RUNTIME_DIR = "<path-to-mingw-bin>"
$env:STAPPP_C_COMPILER = "$env:STAPPP_COMPILER_RUNTIME_DIR\gcc.exe"
$env:STAPPP_CXX_COMPILER = "$env:STAPPP_COMPILER_RUNTIME_DIR\g++.exe"
```

其中 `<your-mkl-root>` 应该是包含 `include`、`lib` 和运行时 DLL 目录的 MKL 根目录。正式 oneAPI 安装通常类似：

```text
<oneapi-root>\mkl\latest
```

封装的临时 MKL 运行目录通常类似：

```text
<package-or-cache-root>\Library
```

## 2. 变量总表

| 变量 | 必需性 | 使用位置 | 作用 |
| --- | --- | --- | --- |
| `STAPPP_MKL_ROOT` | 推荐 | PowerShell 脚本、CMake、VS Code | MKL 根目录。CMake 可由它推导 include、lib 和运行时目录。 |
| `MKLROOT` | 兼容可选 | PowerShell 脚本、CMake、VS Code | Intel oneAPI 常用变量。未设置 `STAPPP_MKL_ROOT` 时作为后备。 |
| `STAPPP_MKL_INCLUDE_DIR` | 可选 | CMake | MKL 头文件目录，通常是 `<mkl-root>\include`。 |
| `STAPPP_MKL_LIBRARY` | 可选 | CMake | MKL 导入库，Windows 通常是 `<mkl-root>\lib\mkl_rt.lib`。 |
| `STAPPP_MKL_RUNTIME_DIR` | 推荐 | CMake 运行目标、VS Code、脚本 | MKL DLL 所在目录。oneAPI 安装常用 `<mkl-root>\redist\intel64`；封装目录常用 `<mkl-root>\bin`。 |
| `STAPPP_COMPILER_RUNTIME_DIR` | 推荐 | CMake 运行目标、VS Code、脚本 | 编译器运行时 DLL 所在目录。MinGW 通常是包含 `g++.exe` 的 `bin` 目录。 |
| `STAPPP_EXTRA_RUNTIME_DIRS` | 可选 | CMake 运行目标、脚本 | 额外运行时目录列表；Windows 用分号分隔。 |
| `STAPPP_CMAKE_EXE` | VS Code 推荐 | `open_vscode_env.ps1`、VS Code | 指定 `cmake.exe` 完整路径。未设置时脚本会从 `PATH` 和常见安装路径查找。 |
| `STAPPP_MINGW_BIN` | 可选 | `open_vscode_env.ps1` | MinGW `bin` 目录后备变量；优先级低于 `STAPPP_COMPILER_RUNTIME_DIR`。 |
| `STAPPP_C_COMPILER` | VS Code/MinGW 推荐 | VS Code Kit、CMake Tools | C 编译器完整路径，例如 `<mingw-bin>\gcc.exe`。 |
| `STAPPP_CXX_COMPILER` | VS Code/MinGW 推荐 | VS Code Kit、CMake Tools | C++ 编译器完整路径，例如 `<mingw-bin>\g++.exe`。 |
| `PATH` | 运行时必需 | PowerShell、VS Code、程序运行 | 运行 `stap++.exe` 时必须能找到 MKL DLL 和编译器运行时 DLL。CMake `run-bridge-*` 目标会自动注入配置过的运行时目录。 |

## 3. CMake 缓存变量

这些不是必须作为环境变量存在，也可以通过 `cmake -D...` 显式传入：

```powershell
cmake -S . -B build-mkl -G "MinGW Makefiles" `
  -DCMAKE_BUILD_TYPE=Release `
  -DCMAKE_C_COMPILER="$env:STAPPP_C_COMPILER" `
  -DCMAKE_CXX_COMPILER="$env:STAPPP_CXX_COMPILER" `
  -DSTAPPP_ENABLE_MKL_PARDISO=ON `
  -DSTAPPP_MKL_ROOT="$env:STAPPP_MKL_ROOT" `
  -DSTAPPP_MKL_INCLUDE_DIR="$env:STAPPP_MKL_INCLUDE_DIR" `
  -DSTAPPP_MKL_LIBRARY="$env:STAPPP_MKL_LIBRARY" `
  -DSTAPPP_MKL_RUNTIME_DIR="$env:STAPPP_MKL_RUNTIME_DIR" `
  -DSTAPPP_COMPILER_RUNTIME_DIR="$env:STAPPP_COMPILER_RUNTIME_DIR" `
  -DSTAPPP_EXTRA_RUNTIME_DIRS="$env:STAPPP_EXTRA_RUNTIME_DIRS"
```

如果只传：

```powershell
-DSTAPPP_MKL_ROOT="$env:STAPPP_MKL_ROOT"
```

CMake 会尝试自动推导：

```text
STAPPP_MKL_INCLUDE_DIR = <mkl-root>/include
STAPPP_MKL_LIBRARY = <mkl-root>/lib/mkl_rt.lib
STAPPP_MKL_RUNTIME_DIR = <mkl-root>/redist/intel64 或 <mkl-root>/bin
```

## 4. 脚本参数与环境变量对应关系

`scripts\configure_mkl.ps1` 支持：

```powershell
.\scripts\configure_mkl.ps1 `
  -BuildDir build-mkl `
  -Generator "MinGW Makefiles" `
  -CMakeExe "$env:STAPPP_CMAKE_EXE" `
  -MklRoot "$env:STAPPP_MKL_ROOT" `
  -MklIncludeDir "$env:STAPPP_MKL_INCLUDE_DIR" `
  -MklLibrary "$env:STAPPP_MKL_LIBRARY" `
  -MklRuntimeDir "$env:STAPPP_MKL_RUNTIME_DIR" `
  -CompilerRuntimeDir "$env:STAPPP_COMPILER_RUNTIME_DIR" `
  -ExtraRuntimeDirs "$env:STAPPP_EXTRA_RUNTIME_DIRS"
```

`scripts\open_vscode_env.ps1` 支持：

```powershell
.\scripts\open_vscode_env.ps1 `
  -ProjectDrive X: `
  -MklDrive Y: `
  -ProjectPath "<path-to-STAPpp_evaluation_clean>" `
  -MklPath "$env:STAPPP_MKL_ROOT" `
  -CMakeExe "$env:STAPPP_CMAKE_EXE" `
  -CompilerBin "$env:STAPPP_COMPILER_RUNTIME_DIR" `
  -CCompiler "$env:STAPPP_C_COMPILER" `
  -CxxCompiler "$env:STAPPP_CXX_COMPILER"
```

参数优先级高于环境变量。

## 5. VS Code 使用规则

`.vscode\settings.json` 和 `.vscode\cmake-kits.json` 不再写死本机路径，而是读取 VS Code 启动进程中的 `STAPPP_*` 变量。

推荐用：

```powershell
.\scripts\open_vscode_env.ps1
```

启动 VS Code。该脚本会设置：

```text
STAPPP_MKL_ROOT
STAPPP_MKL_INCLUDE_DIR
STAPPP_MKL_LIBRARY
STAPPP_MKL_RUNTIME_DIR
STAPPP_COMPILER_RUNTIME_DIR
STAPPP_C_COMPILER
STAPPP_CXX_COMPILER
STAPPP_CMAKE_EXE
MKLROOT
PATH
```

如果不用脚本启动 VS Code，则需要在启动 VS Code 前由用户自己设置上述变量。

