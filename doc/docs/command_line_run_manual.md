# STAP++ 课程项目命令行运行手册

生成时间：2026-06-13

本文档用于课程项目现场测评或迁移到新机器后复现实验环境。范围从干净测评仓库、`cmake` 配置、编译开始，到 `Bridge-1` 至 `Bridge-4` 四个算例运行结束。

优先使用干净测评仓库：

```text
STAPpp_evaluation_clean
```

该目录只包含测评运行必需的源码、输入文件、构建脚本、运行脚本和本文档，不包含 `.git`、`tmp`、历史构建目录、旧可执行文件、历史 `.out/.csv/.vtk` 结果和 Doxygen HTML。

## 1. 适用环境

测评说明中与运行直接相关的要求包括：

- 使用 `cmake` 作为编译入口。
- 运行 `Bridge-1` 至 `Bridge-4`。
- 程序运行不能崩溃，内存不超过物理内存限制。
- 初始向量为零向量。
- 结果需要能用于 Tecplot 或 ParaView 后处理。
- 需要记录 walltime 和内存占用情况。

当前程序的大算例主路线使用 Intel oneMKL PARDISO 后端。原因是 `Bridge-3` 和 `Bridge-4` 的自由度规模较大，标准库自实现迭代求解器可以作为保留路径，但当前通过测评时间约束的路线是 `--backend pardiso`。

## 2. 目录约定

以下命令均假定当前目录为干净测评仓库根目录：

```powershell
cd <path-to-STAPpp_evaluation_clean>
```

迁移到其他电脑时，将上面的路径替换为新的 `STAPpp_evaluation_clean` 目录即可。本文后续命令均使用相对路径，不依赖本机绝对路径。

### 2.1 干净仓库内容

```text
STAPpp_evaluation_clean
|-- CMakeLists.txt
|-- README.md
|-- README_EVALUATION.md
|-- src
|   |-- cpp
|   `-- h
|-- Bridge-1
|   `-- Bridge-1.dat
|-- Bridge-2
|   `-- Bridge-2.dat
|-- Bridge-3
|   `-- Bridge-3.dat
|-- Bridge-4.dat
|-- .vscode
|   |-- settings.json
|   `-- cmake-kits.json
|-- tools
|   `-- inp2dat
|-- scripts
|   |-- configure_mkl.ps1
|   |-- run_all_bridges.ps1
|   `-- check_results.ps1
`-- docs
    |-- command_line_run_manual.md
    |-- environment_variables.md
    `-- vscode_cmake_bridge_full_workflow.md
```

封装时有意排除：

```text
.git
.claude
tmp
build
results
*.exe
*.out
*.csv
*.vtk
Doxygen HTML
历史调试脚本和历史验证输出
```

## 3. 软件依赖

最低依赖：

- Windows 10 或 Windows 11。
- CMake 3.10 或更高。
- C++ 编译器：Microsoft Visual Studio 2022 MSVC，或 MinGW-w64 GCC。
- Python 3，用于从 `.inp` 重新转换 `.dat` 时运行 `tools\inp2dat\inp2dat.py`。
- Intel oneAPI MKL，用于启用 PARDISO 后端。

迁移到新机器时，不要直接修改源码或 `.vscode` 文件中的路径。优先设置环境变量，完整清单见：

```text
docs\environment_variables.md
```

推荐最小 PowerShell 配置模板如下，尖括号内容按本机实际路径替换：

```powershell
$env:STAPPP_MKL_ROOT = "<your-mkl-root>"
$env:STAPPP_CMAKE_EXE = "<path-to-cmake.exe>"
$env:STAPPP_COMPILER_RUNTIME_DIR = "<path-to-mingw-bin>"
$env:STAPPP_C_COMPILER = "$env:STAPPP_COMPILER_RUNTIME_DIR\gcc.exe"
$env:STAPPP_CXX_COMPILER = "$env:STAPPP_COMPILER_RUNTIME_DIR\g++.exe"
```

如果使用 Intel oneAPI 正式安装，`STAPPP_MKL_ROOT` 通常指向 `mkl\latest`。如果使用封装的 MKL 运行目录，`STAPPP_MKL_ROOT` 应指向同时包含 `include`、`lib`、`bin` 的目录。

检查工具是否可用：

```powershell
& $env:STAPPP_CMAKE_EXE --version
python --version
where python
```

如果使用 MinGW：

```powershell
& $env:STAPPP_CXX_COMPILER --version
```

如果使用 MSVC：

```powershell
cl
```

## 4. 从 CMake 开始配置和编译

### 4.0 推荐一键配置脚本

在干净仓库根目录执行：

```powershell
.\scripts\configure_mkl.ps1 `
  -BuildDir build-mkl `
  -Generator "MinGW Makefiles"
```

该脚本会：

- 检查 `STAPPP_MKL_ROOT`、`STAPPP_MKL_INCLUDE_DIR` 或脚本参数指定的 MKL include 目录。
- 检查 `STAPPP_MKL_LIBRARY` 或 `<mkl-root>\lib\mkl_rt.lib`。
- 将 `STAPPP_MKL_RUNTIME_DIR` 和 `STAPPP_COMPILER_RUNTIME_DIR` 加入当前 PowerShell 进程的 `PATH`。
- 使用 `cmake` 配置 `STAPPP_ENABLE_MKL_PARDISO=ON`。
- 编译生成 `stap++.exe`。

如果现场使用 Visual Studio 2022 生成器：

```powershell
.\scripts\configure_mkl.ps1 -UseVisualStudio
```

如果需要显式指定生成器：

```powershell
.\scripts\configure_mkl.ps1 -Generator "MinGW Makefiles"
```

如果当前终端找不到 `cmake`，但 Visual Studio BuildTools 自带 CMake，可显式传入：

```powershell
.\scripts\configure_mkl.ps1 `
  -Generator "MinGW Makefiles" `
  -CMakeExe "$env:STAPPP_CMAKE_EXE"
```

如果工程路径或 MKL 路径包含中文，旧版 CMake/MinGW 组合可能在生成文件或编译参数中出现路径编码问题。实测可用的处理方式是将工程目录和 MKL 目录临时映射到 ASCII 盘符：

```powershell
subst X: "<path-to-STAPpp_evaluation_clean>"
subst Y: "<path-to-mkl-root>"

X:
$env:STAPPP_MKL_ROOT = "Y:\"
.\scripts\configure_mkl.ps1 `
  -BuildDir build-mkl `
  -Generator "MinGW Makefiles" `
  -MklRoot "$env:STAPPP_MKL_ROOT"
```

如果使用正式 oneAPI 安装目录，通常不需要 `Y:` 映射，只需要设置 `MKLROOT`。

### 4.1 MinGW-w64 推荐命令

```powershell
cmake -S . -B build-mkl -G "MinGW Makefiles" `
  -DCMAKE_BUILD_TYPE=Release `
  -DCMAKE_C_COMPILER="$env:STAPPP_C_COMPILER" `
  -DCMAKE_CXX_COMPILER="$env:STAPPP_CXX_COMPILER" `
  -DSTAPPP_ENABLE_MKL_PARDISO=ON `
  -DSTAPPP_MKL_ROOT="$env:STAPPP_MKL_ROOT" `
  -DSTAPPP_MKL_RUNTIME_DIR="$env:STAPPP_MKL_RUNTIME_DIR" `
  -DSTAPPP_COMPILER_RUNTIME_DIR="$env:STAPPP_COMPILER_RUNTIME_DIR"

cmake --build build-mkl --config Release
```

生成的程序通常位于：

```text
build-mkl\stap++.exe
```

### 4.2 VS Code CMake Tools 可迁移配置

本仓库已补充 VS Code CMake Tools 配置文件：

```text
.vscode\settings.json
.vscode\cmake-kits.json
scripts\open_vscode_env.ps1
```

使用方式：

```powershell
cd <path-to-STAPpp_evaluation_clean>
.\scripts\open_vscode_env.ps1
```

迁移到新机器时可显式传入路径：

```powershell
cd <path-to-STAPpp_evaluation_clean>
.\scripts\open_vscode_env.ps1 `
  -ProjectPath "<path-to-STAPpp_evaluation_clean>" `
  -MklPath "$env:STAPPP_MKL_ROOT" `
  -CMakeExe "$env:STAPPP_CMAKE_EXE" `
  -CompilerBin "$env:STAPPP_COMPILER_RUNTIME_DIR"
```

脚本会将以下变量加入 VS Code 启动进程的环境：

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

在 VS Code 中操作：

1. 安装 `CMake Tools` 和 `C/C++` 插件。
2. `Ctrl + Shift + P`。
3. 执行 `CMake: Select a Kit`。
4. 选择 `STAP++ MinGW from environment`。
5. 执行 `CMake: Configure`。
6. 执行 `CMake: Build`。
7. 如需直接运行算例，执行 `CMake: Build Target`，选择 `run-bridge-1`、`run-bridge-2`、`run-bridge-3` 或 `run-bridge-4`。

CMake 已为这些运行目标注入 `STAPPP_MKL_RUNTIME_DIR`、`STAPPP_COMPILER_RUNTIME_DIR` 和 `STAPPP_EXTRA_RUNTIME_DIRS` 指定的运行期路径。因此通过 CMake 目标运行时，不需要手动设置 `$env:PATH`。

如果没有看到 `STAP++ MinGW from environment`，执行：

```text
CMake: Scan for Kits
```

或者确认当前 VS Code 是通过 `scripts\open_vscode_env.ps1` 启动的。

### 4.3 Visual Studio 2022 可选命令

```powershell
cmake -S . -B build-mkl-vs -G "Visual Studio 17 2022" -A x64 `
  -DSTAPPP_ENABLE_MKL_PARDISO=ON `
  -DSTAPPP_MKL_ROOT="$env:STAPPP_MKL_ROOT" `
  -DSTAPPP_MKL_RUNTIME_DIR="$env:STAPPP_MKL_RUNTIME_DIR"

cmake --build build-mkl-vs --config Release
```

生成的程序通常位于：

```text
build-mkl-vs\Release\stap++.exe
```

后续命令统一用 `$exe` 指向可执行文件。按实际构建方式选择其中一个：

```powershell
$exe = ".\build-mkl\stap++.exe"
```

或：

```powershell
$exe = ".\build-mkl-vs\Release\stap++.exe"
```

检查程序参数：

```powershell
& $exe
```

正常情况下会输出：

```text
Usage: stap++ InputFileName [--solver skyline|sparse-cg|sparse-bicgstab|sparse-gmres|sparse-auto] ...
```

## 5. 输入文件准备

干净测评仓库已经包含可直接运行的 `.dat` 文件：

```text
Bridge-1\Bridge-1.dat
Bridge-2\Bridge-2.dat
Bridge-3\Bridge-3.dat
Bridge-4.dat
```

说明：

- `Bridge-1\Bridge-1.dat` 是当前已验证的 Bridge-1 输入，包含 4107 个节点和 66 条可化简 MPC。
- 历史文件 `Bridge-1\Bridge-1-mpc.dat` 只保留为旧实验记录，不作为当前默认运行入口。
- `Bridge-3\Bridge-3.dat` 和 `Bridge-4.dat` 是已经转换好的大算例输入。

如果需要从 Abaqus `.inp` 重新生成 `.dat`，可使用：

```powershell
python tools\inp2dat\inp2dat.py Bridge-3\Bridge-3.inp Bridge-3\Bridge-3.dat `
  --solid-type H8RPIER `
  --pier-instances Part-Pier `
  --tie-mode auto `
  --node-order rcm
```

`Bridge-4.dat` 当前已经存在，通常不建议在现场重新转换，除非需要从原始 `.inp` 完整复现数据准备流程。

## 6. 通用运行参数

PARDISO 主路线使用以下公共参数：

```text
--solver sparse-auto
--backend pardiso
--pardiso-mtype auto
--output summary
--tol 1e-6
--max-iter 5000
```

说明：

- `--solver sparse-auto`：进入稀疏求解路径。
- `--backend pardiso`：使用 Intel oneMKL PARDISO 直接求解。
- `--pardiso-mtype auto`：先按 SPD 矩阵尝试，必要时程序可切换。
- `--output summary`：只输出关键规模、残差、计时和 CSV 路径，避免大算例输出巨型 `.out`。
- `--tol 1e-6`：复核相对残差必须不大于该阈值。
- `--max-iter 5000`：保留接口参数。PARDISO 是直接法，成功时输出 `ITERATIONS = 0`。

程序内部会在求解后复核残差。若 `CHECKED RELATIVE RESIDUAL > --tol`，程序退出，不写位移 CSV。

## 7. Bridge-1 到 Bridge-4 运行命令

### 7.1 Bridge-1

```powershell
& $exe Bridge-1\Bridge-1.dat `
  --solver sparse-auto `
  --backend pardiso `
  --pardiso-mtype auto `
  --output summary `
  --tol 1e-6 `
  --max-iter 5000 `
  --csv results\Bridge-1.displacements.csv
```

输出文件：

```text
Bridge-1\Bridge-1.out
results\Bridge-1.displacements.csv
```

当前已验证结果：

```text
NEQ = 15122
NNZ = 349613
ACTUAL SOLVER = pardiso
ITERATIVE SOLVER CONVERGED = YES
CHECKED RELATIVE RESIDUAL = 6.61471e-13
TOTAL_TIME = 3.52000e-01 s
```

### 7.2 Bridge-2

```powershell
& $exe Bridge-2\Bridge-2.dat `
  --solver sparse-auto `
  --backend pardiso `
  --pardiso-mtype sym-indef `
  --output summary `
  --tol 1e-6 `
  --max-iter 5000 `
  --csv results\Bridge-2.displacements.csv
```

输出文件：

```text
Bridge-2\Bridge-2.out
results\Bridge-2.displacements.csv
```

`Bridge-2` 已知按 `spd` 首次分解会失败，因此建议直接指定 `--pardiso-mtype sym-indef`，避免 `auto` 模式下先失败再重试。若重新运行，成功判据以新生成的 `Bridge-2\Bridge-2.out` 为准：

```text
ITERATIVE SOLVER CONVERGED = YES
CHECKED RELATIVE RESIDUAL <= 1.0e-6
DISPLACEMENT CSV = results\Bridge-2.displacements.csv
TOTAL_TIME = ...
```

### 7.3 Bridge-3

```powershell
& $exe Bridge-3\Bridge-3.dat `
  --solver sparse-auto `
  --backend pardiso `
  --pardiso-mtype auto `
  --output summary `
  --tol 1e-6 `
  --max-iter 5000 `
  --csv results\Bridge-3.displacements.csv
```

输出文件：

```text
Bridge-3\Bridge-3.out
results\Bridge-3.displacements.csv
```

当前已验证结果：

```text
NEQ = 809526
NNZ = 29646042
ACTUAL SOLVER = pardiso
ITERATIVE SOLVER CONVERGED = YES
CHECKED RELATIVE RESIDUAL = 1.71609e-11
PARDISO_FACT_NNZ = 617189494
TOTAL_TIME = 2.16730e+01 s
```

### 7.4 Bridge-4

```powershell
& $exe Bridge-4.dat `
  --solver sparse-auto `
  --backend pardiso `
  --pardiso-mtype auto `
  --output summary `
  --tol 1e-6 `
  --max-iter 5000 `
  --csv results\Bridge-4.displacements.csv
```

输出文件：

```text
Bridge-4.out
results\Bridge-4.displacements.csv
```

当前已验证结果：

```text
NEQ = 5837346
NNZ = 226464642
ACTUAL SOLVER = pardiso
ITERATIVE SOLVER CONVERGED = YES
CHECKED RELATIVE RESIDUAL = 4.72514e-11
PARDISO_FACT_NNZ = 1552010006
TOTAL_TIME = 1.44876e+03 s
```

## 8. 单独运行四个算例

测评默认建议各算例单独运行，便于记录每个算例的退出码、walltime、残差和输出文件。

如果使用 VS Code CMake Tools，推荐：

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

CMake 会自动携带 MKL 和 MinGW 运行期路径。

如果在 PowerShell 里手动运行，以下命令假定已设置：

```powershell
$env:PATH = "$env:STAPPP_MKL_RUNTIME_DIR;$env:STAPPP_COMPILER_RUNTIME_DIR;$env:PATH"
$exe = ".\build-mkl\stap++.exe"
New-Item -ItemType Directory -Force -Path results | Out-Null
```

### 8.1 Bridge-1 单独运行

```powershell
& $exe Bridge-1\Bridge-1.dat --solver sparse-auto --backend pardiso --pardiso-mtype auto --output summary --tol 1e-6 --max-iter 5000 --csv results\Bridge-1.displacements.csv
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
```

### 8.2 Bridge-2 单独运行

```powershell
& $exe Bridge-2\Bridge-2.dat --solver sparse-auto --backend pardiso --pardiso-mtype sym-indef --output summary --tol 1e-6 --max-iter 5000 --csv results\Bridge-2.displacements.csv
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
```

说明：`Bridge-2` 实测按 `spd` 首次分解会失败，因此单独运行命令直接使用 `sym-indef`。只要最终输出 `ITERATIVE SOLVER CONVERGED = YES` 且 `CHECKED RELATIVE RESIDUAL <= 1e-6`，该算例通过。

### 8.3 Bridge-3 单独运行

```powershell
& $exe Bridge-3\Bridge-3.dat --solver sparse-auto --backend pardiso --pardiso-mtype auto --output summary --tol 1e-6 --max-iter 5000 --csv results\Bridge-3.displacements.csv
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
```

### 8.4 Bridge-4 单独运行

```powershell
& $exe Bridge-4.dat --solver sparse-auto --backend pardiso --pardiso-mtype auto --output summary --tol 1e-6 --max-iter 5000 --csv results\Bridge-4.displacements.csv
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
```

## 9. 可选一键顺序运行脚本

一键脚本只作为辅助入口，不作为默认测评流程。

完成第 4 节编译后，在干净仓库根目录执行：

```powershell
.\scripts\run_all_bridges.ps1
```

该脚本会顺序运行：

```text
Bridge-1\Bridge-1.dat
Bridge-2\Bridge-2.dat
Bridge-3\Bridge-3.dat
Bridge-4.dat
```

并将位移 CSV 写入：

```text
results\Bridge-1.displacements.csv
results\Bridge-2.displacements.csv
results\Bridge-3.displacements.csv
results\Bridge-4.displacements.csv
```

脚本会在任一算例返回非零退出码时停止。

运行结束后检查：

```powershell
.\scripts\check_results.ps1
```

### 9.1 手动连续运行命令

如果需要从命令行连续运行四个算例，可直接执行以下 PowerShell 命令。该命令不会重新编译，只运行已配置好的 `$exe`。

```powershell
& $exe Bridge-1\Bridge-1.dat --solver sparse-auto --backend pardiso --pardiso-mtype auto --output summary --tol 1e-6 --max-iter 5000 --csv results\Bridge-1.displacements.csv
& $exe Bridge-2\Bridge-2.dat     --solver sparse-auto --backend pardiso --pardiso-mtype sym-indef --output summary --tol 1e-6 --max-iter 5000 --csv results\Bridge-2.displacements.csv
& $exe Bridge-3\Bridge-3.dat     --solver sparse-auto --backend pardiso --pardiso-mtype auto --output summary --tol 1e-6 --max-iter 5000 --csv results\Bridge-3.displacements.csv
& $exe Bridge-4.dat              --solver sparse-auto --backend pardiso --pardiso-mtype auto --output summary --tol 1e-6 --max-iter 5000 --csv results\Bridge-4.displacements.csv
```

如果希望某个算例失败后立即停止，使用：

```powershell
$ErrorActionPreference = "Stop"
& $exe Bridge-1\Bridge-1.dat --solver sparse-auto --backend pardiso --pardiso-mtype auto --output summary --tol 1e-6 --max-iter 5000 --csv results\Bridge-1.displacements.csv
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
& $exe Bridge-2\Bridge-2.dat --solver sparse-auto --backend pardiso --pardiso-mtype sym-indef --output summary --tol 1e-6 --max-iter 5000 --csv results\Bridge-2.displacements.csv
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
& $exe Bridge-3\Bridge-3.dat --solver sparse-auto --backend pardiso --pardiso-mtype auto --output summary --tol 1e-6 --max-iter 5000 --csv results\Bridge-3.displacements.csv
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
& $exe Bridge-4.dat --solver sparse-auto --backend pardiso --pardiso-mtype auto --output summary --tol 1e-6 --max-iter 5000 --csv results\Bridge-4.displacements.csv
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
```

## 10. 结果检查命令

运行结束后，用以下命令提取四个输出文件中的关键结果：

```powershell
Select-String -LiteralPath `
  Bridge-1\Bridge-1.out, `
  Bridge-2\Bridge-2.out, `
  Bridge-3\Bridge-3.out, `
  Bridge-4.out `
  -Pattern "NUMBER OF EQUATIONS|NUMBER OF CSR NONZEROS|ACTUAL SOLVER|ITERATIVE SOLVER CONVERGED|CHECKED RELATIVE RESIDUAL|PARDISO_FACT_NNZ|TOTAL_TIME|DISPLACEMENT CSV"
```

成功标准：

```text
ACTUAL SOLVER = pardiso
ITERATIVE SOLVER CONVERGED = YES
CHECKED RELATIVE RESIDUAL <= 1.0e-6
DISPLACEMENT CSV = ...
TOTAL_TIME = ...
```

检查 CSV 文件是否生成：

```powershell
Get-Item `
  results\Bridge-1.displacements.csv, `
  results\Bridge-2.displacements.csv, `
  results\Bridge-3.displacements.csv, `
  results\Bridge-4.displacements.csv |
  Select-Object Length,LastWriteTime,FullName
```

## 11. 本机封装后验证结果

以下结果来自 2026-06-13 在封装后的 `STAPpp_evaluation_clean` 中逐个单独运行四个算例，构建路径使用 `subst X:` 映射到 ASCII 盘符，MKL 路径使用 `subst Y:` 映射到 ASCII 盘符。

| 算例 | NEQ | NNZ | PARDISO 类型 | 尝试次数 | 复核相对残差 | TOTAL_TIME |
| --- | ---: | ---: | --- | ---: | ---: | ---: |
| Bridge-1 | 15122 | 349613 | spd | 1 | 6.61471e-13 | 4.11000e-01 s |
| Bridge-2 | 120410 | 3914399 | sym-indef | 1 | 4.57451e-12 | 4.02800e+00 s |
| Bridge-3 | 809526 | 29646042 | spd | 1 | 1.72596e-11 | 2.28980e+01 s |
| Bridge-4 | 5837346 | 226464642 | spd | 1 | 4.66155e-11 | 1.44114e+03 s |

对应 CSV 输出：

```text
results\Bridge-1.displacements.csv
results\Bridge-2.displacements.csv
results\Bridge-3.displacements.csv
results\Bridge-4.displacements.csv
```

## 12. 后处理说明

当前程序直接写出位移 CSV，文件包含节点位移结果，可作为误差分析或可视化前处理输入：

```text
results\Bridge-1.displacements.csv
results\Bridge-2.displacements.csv
results\Bridge-3.displacements.csv
results\Bridge-4.displacements.csv
```

干净测评仓库默认不携带历史 `.vtk` 结果文件。若需要把新的 CSV 位移合并回 VTK 或 Tecplot 格式，应单独执行后处理脚本，并记录转换命令、输入 CSV、输出可视化文件和误差分析结果。

## 13. 常见问题

### 13.1 CMake 找不到 MKL

现象：

```text
STAPPP_ENABLE_MKL_PARDISO=ON requires STAPPP_MKL_INCLUDE_DIR ...
```

处理：

```powershell
$env:STAPPP_MKL_ROOT = "<your-mkl-root>"
Test-Path "$env:STAPPP_MKL_ROOT\include"
Test-Path "$env:STAPPP_MKL_ROOT\lib\mkl_rt.lib"
```

确认路径存在后重新运行 `cmake` 配置命令。

### 13.2 运行时找不到 `mkl_rt.dll`

处理：

```powershell
$env:PATH = "$env:STAPPP_MKL_RUNTIME_DIR;$env:STAPPP_COMPILER_RUNTIME_DIR;$env:PATH"
```

然后重新运行算例。

### 13.3 PARDISO 后端不可用

如果输出：

```text
BACKEND = pardiso
BACKEND_AVAILABLE = NO
```

说明当前可执行文件没有启用 `STAPPP_ENABLE_MKL_PARDISO=ON`，或运行环境未正确加载 MKL。应重新执行第 4 节的 CMake 配置和编译命令。

### 13.4 Bridge-4 内存压力过大

`Bridge-4` 当前规模：

```text
NEQ = 5837346
NNZ = 226464642
PARDISO_FACT_NNZ = 1552010006
```

该算例对内存和连续可用内存要求较高。现场运行前建议关闭不必要程序，并使用 `--output summary` 避免输出完整节点、单元和应力表。
