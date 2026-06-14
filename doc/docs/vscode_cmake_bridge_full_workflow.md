# VS Code + CMake Tools 运行 Bridge 完整流程

生成时间：2026-06-13

本文档用于使用 VS Code、CMake Tools、MinGW 和 MKL/PARDISO 调用 `STAP++` 程序，完成 `Bridge-1` 至 `Bridge-4` 算例运行。路径配置通过环境变量或脚本参数传入，不依赖某一台机器的固定安装路径。

## 1. 前提

需要准备的组件：

```text
MinGW-w64 或 MSVC
CMake
Intel oneAPI MKL 或封装的 MKL 运行目录
VS Code
CMake Tools
C/C++
```

建议先阅读环境变量清单：

```text
docs\environment_variables.md
```

如果工程路径或 MKL 路径包含中文，建议通过 `subst` 映射 ASCII 盘符运行；脚本默认使用 `X:` 映射工程目录、`Y:` 映射 MKL 目录，这两个盘符可通过参数修改。

## 2. 打开 VS Code

打开 PowerShell：

```powershell
cd <path-to-STAPpp_evaluation_clean>
.\scripts\open_vscode_env.ps1
```

迁移到新机器时，建议显式传入关键路径：

```powershell
.\scripts\open_vscode_env.ps1 `
  -ProjectPath "<path-to-STAPpp_evaluation_clean>" `
  -MklPath "<path-to-mkl-root>" `
  -CMakeExe "<path-to-cmake.exe>" `
  -CompilerBin "<path-to-mingw-bin>"
```

该脚本会设置 `STAPPP_*` 环境变量，并以映射后的工程目录打开 VS Code。

如果脚本成功，会看到类似：

```text
Project mapped to X:\
MKL mapped to Y:\
CMake: <path-to-cmake.exe>
C compiler: <path-to-gcc.exe>
C++ compiler: <path-to-g++.exe>
Opening VS Code...
```

## 3. 安装 VS Code 插件

确认安装以下插件：

```text
CMake Tools
C/C++
```

对应插件通常是：

```text
ms-vscode.cmake-tools
ms-vscode.cpptools
```

## 4. 选择 Kit

在 VS Code 中：

```text
Ctrl + Shift + P
CMake: Select a Kit
```

选择：

```text
STAP++ MinGW from environment
```

如果没有看到，执行：

```text
CMake: Scan for Kits
```

然后重新执行：

```text
CMake: Select a Kit
```

说明：

- 本仓库已提供 `.vscode\cmake-kits.json`。
- CMake Tools 应能从当前工程的 `.vscode\cmake-kits.json` 读取该 Kit。
- 旧的 `VisualStudio.11.0` warning 可以忽略，那是用户全局 CMake Tools 缓存里的旧 Kit。

## 5. 配置 CMake

执行：

```text
Ctrl + Shift + P
CMake: Configure
```

成功标志：

```text
Configuring done
Generating done
Build files have been written to: X:/build-mkl
```

配置中应包含：

```text
-DSTAPPP_ENABLE_MKL_PARDISO=ON
-DSTAPPP_MKL_ROOT=<path-to-mkl-root>
-DSTAPPP_MKL_RUNTIME_DIR=<path-to-mkl-runtime-dir>
-G "MinGW Makefiles"
```

## 6. 编译

执行：

```text
Ctrl + Shift + P
CMake: Build
```

成功后应生成：

```text
build-mkl\stap++.exe
```

可以在 VS Code Terminal 中检查：

```powershell
Test-Path .\build-mkl\stap++.exe
```

输出应为：

```text
True
```

## 7. 准备运行

推荐直接使用 CMake 运行目标。CMake 已配置运行期环境变量，会在运行目标中自动设置：

```text
PATH = <mkl-runtime-dir>;<compiler-runtime-dir>
```

因此不需要在 PowerShell 中手动设置 `$env:PATH`。

在 VS Code 中执行：

```text
Ctrl + Shift + P
CMake: Build Target
```

然后选择：

```text
run-bridge-1
run-bridge-2
run-bridge-3
run-bridge-4
```

其中：

- `run-bridge-2` 已固定使用 `--pardiso-mtype sym-indef`。
- 其他 Bridge 算例使用 `--pardiso-mtype auto`。

## 8. 备用：在 Terminal 中单独运行 Bridge-1

如果不用 CMake 运行目标，也可以在 VS Code Terminal 中手动执行：

```powershell
$env:PATH = "$env:STAPPP_MKL_RUNTIME_DIR;$env:STAPPP_COMPILER_RUNTIME_DIR;$env:PATH"
$exe = ".\build-mkl\stap++.exe"
New-Item -ItemType Directory -Force -Path results | Out-Null
```

```powershell
& $exe Bridge-1\Bridge-1.dat --solver sparse-auto --backend pardiso --pardiso-mtype auto --output summary --tol 1e-6 --max-iter 5000 --csv results\Bridge-1.displacements.csv
```

检查：

```powershell
Select-String -LiteralPath Bridge-1\Bridge-1.out -Pattern "ACTUAL SOLVER|ITERATIVE SOLVER CONVERGED|CHECKED RELATIVE RESIDUAL|TOTAL_TIME|DISPLACEMENT CSV"
```

通过标准：

```text
ACTUAL SOLVER = pardiso
ITERATIVE SOLVER CONVERGED = YES
CHECKED RELATIVE RESIDUAL <= 1.0e-6
DISPLACEMENT CSV = results\Bridge-1.displacements.csv
```

## 9. 备用：在 Terminal 中单独运行 Bridge-2

```powershell
& $exe Bridge-2\Bridge-2.dat --solver sparse-auto --backend pardiso --pardiso-mtype sym-indef --output summary --tol 1e-6 --max-iter 5000 --csv results\Bridge-2.displacements.csv
```

检查：

```powershell
Select-String -LiteralPath Bridge-2\Bridge-2.out -Pattern "ACTUAL SOLVER|PARDISO_SELECTED_MTYPE|PARDISO_ATTEMPT_COUNT|PARDISO_RETRY_FROM_SPD_TO_SYM_INDEF|ITERATIVE SOLVER CONVERGED|CHECKED RELATIVE RESIDUAL|TOTAL_TIME|DISPLACEMENT CSV"
```

说明：

```text
Bridge-2 已知按 spd 首次分解会失败，因此这里直接指定 sym-indef，避免 auto 模式下先失败再重试。
如果最终 ITERATIVE SOLVER CONVERGED = YES 且残差 <= 1.0e-6，则通过。
```

## 10. 备用：在 Terminal 中单独运行 Bridge-3

```powershell
& $exe Bridge-3\Bridge-3.dat --solver sparse-auto --backend pardiso --pardiso-mtype auto --output summary --tol 1e-6 --max-iter 5000 --csv results\Bridge-3.displacements.csv
```

检查：

```powershell
Select-String -LiteralPath Bridge-3\Bridge-3.out -Pattern "ACTUAL SOLVER|ITERATIVE SOLVER CONVERGED|CHECKED RELATIVE RESIDUAL|TOTAL_TIME|DISPLACEMENT CSV"
```

## 11. 备用：在 Terminal 中单独运行 Bridge-4

```powershell
& $exe Bridge-4.dat --solver sparse-auto --backend pardiso --pardiso-mtype auto --output summary --tol 1e-6 --max-iter 5000 --csv results\Bridge-4.displacements.csv
```

检查：

```powershell
Select-String -LiteralPath Bridge-4.out -Pattern "ACTUAL SOLVER|ITERATIVE SOLVER CONVERGED|CHECKED RELATIVE RESIDUAL|TOTAL_TIME|DISPLACEMENT CSV"
```

说明：

```text
Bridge-4 耗时较长，本机实测约 1441 s。
```

## 12. 一次性检查全部结果

```powershell
.\scripts\check_results.ps1
```

成功标准：

```text
ACTUAL SOLVER = pardiso
ITERATIVE SOLVER CONVERGED = YES
CHECKED RELATIVE RESIDUAL <= 1.0e-6
DISPLACEMENT CSV = results\*.displacements.csv
```

## 13. 查看计算时间

`.out` 文件中有：

```text
S O L U T I O N   T I M E   L O G   I N   S E C
```

关键字段：

```text
READ_TIME
CSR_ASSEMBLY_TIME
EXPORT_UPPER_CSR_TIME
ITER_SOLVE_TIME
RESIDUAL_CHECK_TIME
CSV_WRITE_TIME
TOTAL_TIME
```

如果要查看除去读入和 CSV 输出后的时间：

```powershell
$out = Get-Content Bridge-4.out
$total = [double](($out | Select-String "TOTAL_TIME").ToString().Split("=")[1])
$read  = [double](($out | Select-String "READ_TIME").ToString().Split("=")[1])
$csv   = [double](($out | Select-String "CSV_WRITE_TIME").ToString().Split("=")[1])
$total - $read - $csv
```

四个算例一起计算：

```powershell
$files = @(
  "Bridge-1\Bridge-1.out",
  "Bridge-2\Bridge-2.out",
  "Bridge-3\Bridge-3.out",
  "Bridge-4.out"
)

foreach ($f in $files) {
  $out = Get-Content $f
  $total = [double](($out | Select-String "TOTAL_TIME").ToString().Split("=")[1])
  $read  = [double](($out | Select-String "READ_TIME").ToString().Split("=")[1])
  $csv   = [double](($out | Select-String "CSV_WRITE_TIME").ToString().Split("=")[1])

  [pscustomobject]@{
    Case = $f
    TOTAL_TIME = $total
    READ_TIME = $read
    CSV_WRITE_TIME = $csv
    COMPUTE_TIME = $total - $read - $csv
  }
}
```

## 14. 常见问题

### 14.1 Select a Kit 看不到 MinGW

先确认是通过脚本打开的 VS Code：

```powershell
.\scripts\open_vscode_env.ps1
```

然后在 VS Code 中执行：

```text
CMake: Scan for Kits
CMake: Select a Kit
```

选择：

```text
STAP++ MinGW from environment
```

### 14.2 Configure 成功但手动运行时报找不到 MKL DLL

优先使用 CMake 运行目标 `run-bridge-*`，它会自动携带 `PATH`。如果手动运行，在 VS Code Terminal 中执行：

```powershell
$env:PATH = "$env:STAPPP_MKL_RUNTIME_DIR;$env:STAPPP_COMPILER_RUNTIME_DIR;$env:PATH"
```

再重新运行 Bridge 命令。

### 14.3 中文路径导致 CMake 或编译失败

不要直接打开中文路径下的目录。使用：

```powershell
.\scripts\open_vscode_env.ps1
```

确认 VS Code 打开的是：

```text
映射后的工程目录，例如 X:\
```

### 14.4 旧 VisualStudio.11.0 Kit warning

类似 warning：

```text
VS installation instance not found for kit "VisualStudio.11.0 ..."
```

可以忽略。它不影响 `STAP++ MinGW from environment`。
