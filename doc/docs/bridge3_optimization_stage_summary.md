# Bridge-3 优化阶段性总结

生成时间：2026-06-12

## 1. 当前目标

本阶段目标是让 STAP++ 在不依赖第三方数值库的前提下，尽可能快地跑通 `Bridge-3` 算例，并为后续 `Bridge-4` 和误差分析保留可迁移、可解释的实现基础。

当前目标已经收窄为：

- 暂不处理 `Bridge-4`。
- 暂不做 Abaqus 误差分析。
- 先保证 `Bridge-3` 能在可接受精度下完成求解。
- 使用 `Bridge-1` 和 `Bridge-2` 作为中小规模验证算例。
- 保留原 skyline LDLT 路径作为小算例参考，不删除原实现。

## 2. 项目架构概览

STAP++ 的核心仍然是原有有限元程序结构：

- `src/cpp/main.cpp`
  - 程序入口。
  - 解析命令行参数。
  - 选择 skyline 或 sparse 迭代求解路径。
  - 输出计时、残差和位移 CSV。

- `src/h/Domain.h` 与 `src/cpp/Domain.cpp`
  - 管理整体有限元模型。
  - 读取节点、荷载、MPC、单元组。
  - 计算方程编号。
  - 调用单元刚度矩阵并装配全局矩阵。
  - 写出位移 CSV。

- `src/h/ElementGroup.h` 与 `src/cpp/ElementGroup.cpp`
  - 读取不同单元组。
  - 根据单元类型分配具体单元对象和材料对象。

- `src/h/Solver.h` 与 `src/cpp/Solver.cpp`
  - 原有 skyline LDLT 求解器。
  - 新增 CSR 稀疏矩阵和迭代求解器。

- `src/h/Outputter.h` 与 `src/cpp/Outputter.cpp`
  - 原有详细输出。
  - 新增 summary 输出模式，避免大算例输出巨大 `.out` 文件。

- `tools/run_bridge3_fast.py`
  - Bridge-3 构建、转换和运行入口。
  - 使用相对路径，不写死本机目录。

- `tools/run_sparse_validation.py`
  - Bridge-1 和 Bridge-2 的稀疏迭代验证入口。

- `CMakeLists.txt`
  - 根目录构建入口。
  - 使用 C++11。
  - 支持 MSVC 和 MinGW。

## 3. 已实现的主要优化思路

### 3.1 从 skyline 直接法转向 CSR 稀疏迭代

原 skyline LDLT 对大规模桥梁算例不适合。`Bridge-3` 的 skyline 存储估计已经达到不可接受规模，且存在索引溢出风险。

因此新增了标准库自实现的 CSR 稀疏矩阵：

- `row_ptr`
- `col_ind`
- `values`
- full CSR 存储
- `matvec`
- diagonal 提取
- Jacobi 预条件器
- SSOR 风格预条件器

当前支持的稀疏求解器：

- `sparse-cg`
- `sparse-bicgstab`
- `sparse-gmres`

其中 `sparse-bicgstab` 当前是默认候选，因为它在 `Bridge-1-mpc` 上比 GMRES 更快。

### 3.2 保留 skyline 作为参考路径

skyline LDLT 没有删除，仍可通过：

```powershell
.\stap++.exe input.dat --solver skyline
```

使用。它适合小算例校验和结果对照，但不建议作为 `Bridge-3` 主路线。

### 3.3 MPC 简单约束转 DOF alias

原先 MPC 约束若用 penalty 形式，会显著恶化条件数。当前实现会在方程编号之前识别简单零右端二项 MPC：

```text
slave_dof - master_dof = 0
```

并转换为 DOF alias，使 slave DOF 直接使用 master DOF 的方程号。

已验证：

- `Bridge-1-mpc` 转换 66 条简单 MPC，剩余 penalty MPC 为 0。
- `Bridge-3` 转换 312 条简单 MPC，剩余 penalty MPC 为 0。

这一步减少了方程数，也避免了 penalty 刚度带来的病态问题。

### 3.4 Summary 输出模式

大算例不再默认输出完整节点表、方程号表、单元表和应力表。summary 模式只输出：

- `NUMNP`
- `NUMEG`
- `NEQ`
- CSR `NNZ`
- 求解器类型
- tolerance
- max iterations
- preconditioner
- residual
- iteration count
- 计时
- 位移 CSV 路径

运行示例：

```powershell
.\stap++.exe Bridge-3.dat --solver sparse-bicgstab --output summary --precond ssor --tol 1e-6 --max-iter 5000 --csv runs\Bridge-3\displacements.csv
```

### 3.5 收敛约束加强

当前迭代求解不允许“失败但继续输出结果”。必须同时满足：

- 求解器内部 `converged = YES`
- 复核残差 `CHECKED RELATIVE RESIDUAL <= --tol`

否则程序退出，不写位移 CSV。

这样做的目的是避免把未收敛结果误用于后续误差分析。

## 4. 条件限制

### 4.1 课程项目限制

当前阶段按最保守规则处理：

- 不使用 Eigen。
- 不使用 MKL。
- 不使用第三方矩阵库或第三方求解器库。
- 迭代求解器初始向量固定为零向量。
- 稀疏矩阵、预条件器和迭代算法均在项目内自实现。

### 4.2 迁移限制

为了保证迁移到另一台电脑仍可用：

- 使用 CMake 作为正式构建入口。
- 使用 C++11 和标准库。
- 脚本使用相对路径。
- 不写死本机路径。
- 不依赖本机专有环境变量。

推荐构建命令：

```powershell
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
```

### 4.3 当前未处理内容

以下内容本阶段没有完成：

- `Bridge-4` 求解。
- Abaqus 误差分析。
- 多线程并行。
- 第三方库求解器。
- 完整 `Bridge-3` 收敛运行。

## 5. 当前验证结果

### 5.1 Bridge-1-mpc

运行设置：

```powershell
.\tmp\verify\stappp_sparse.exe Bridge-1\Bridge-1-mpc.dat --solver sparse-bicgstab --output summary --precond ssor --tol 1e-6 --max-iter 5000 --csv tmp\verify\Bridge-1-mpc.valid.csv
```

结果：

- `NUMNP = 4107`
- `NEQ = 15122`
- `NNZ = 684104`
- simple MPC alias：66 条
- 收敛：是
- 迭代次数：1551
- `CHECKED RELATIVE RESIDUAL = 9.68897e-07`
- 总时间约 `5.126s`

结论：`Bridge-1-mpc` 能通过当前稀疏迭代路径。

### 5.2 Bridge-2

运行设置：

```powershell
.\tmp\verify\stappp_sparse.exe Bridge-2\Bridge-2.dat --solver sparse-bicgstab --output summary --precond ssor --tol 1e-6 --max-iter 5000 --csv tmp\verify\Bridge-2.valid.csv
```

结果：

- `NUMNP = 37065`
- `NEQ = 120410`
- `NNZ = 7708388`
- 无 MPC
- CSR 装配成功
- 5000 次后未收敛
- `CHECKED RELATIVE RESIDUAL = 4.43047e-01`
- 实际运行约 `220s`

结论：读取和装配问题已修复，但当前预条件器不足，`Bridge-2` 尚不能作为通过算例。

### 5.3 Bridge-3 短跑

运行设置：

```powershell
.\tmp\verify\stappp_sparse.exe tmp\baseline\Bridge-3.dat --solver sparse-bicgstab --output summary --precond ssor --tol 1e-6 --max-iter 20 --csv tmp\verify\Bridge-3.short.csv
```

结果：

- `NUMNP = 259567`
- simple MPC alias：312 条
- 剩余 penalty MPC：0
- `NEQ = 809526`
- `NNZ = 58482558`
- 20 次短跑未收敛
- `CHECKED RELATIVE RESIDUAL = 4.18561e-01`
- 短跑耗时约 `128.6s`

结论：`Bridge-3` 能读入并完成 CSR 装配，MPC alias 生效；但尚未完成收敛求解。

## 6. 为什么不建议回到 skyline 主路线

skyline LDLT 的问题不是单纯慢，而是规模上接近不可用：

- `Bridge-3` 的 skyline 存储估计曾达到约 `11.9GB`。
- 转换日志中也出现过 `NWK-est 33731547015 -> 11231123616`，即优化后仍是约 `11.2e9` 个 skyline 条目。
- 若按 double 存储，`11.2e9` 个条目约 `89.8GB`，还不含其他数组和分解过程开销。
- 当前 skyline 存储部分索引使用 `unsigned int`，超过 `4.29e9` 会有溢出风险。

因此 skyline 可用于小算例对照，但不适合 `Bridge-3` 主路线。

## 7. 当前主要瓶颈

### 7.1 收敛性是第一瓶颈

当前最主要瓶颈不是 CSR matvec 的单次速度，而是预条件器不够强。

证据：

- `Bridge-1-mpc` 能收敛，但需要 1551 次。
- `Bridge-2` 在 5000 次后仍未收敛。
- `Bridge-3` 20 次短跑残差仍在 `O(1e-1)` 量级。

这说明当前 Jacobi/SSOR 预条件对桥梁大规模混合单元系统不足。

### 7.2 Bridge-2 是下一阶段门槛

`Bridge-2` 的规模明显小于 `Bridge-3`：

- `Bridge-2`: `NEQ = 120410`, `NNZ = 7708388`
- `Bridge-3`: `NEQ = 809526`, `NNZ = 58482558`

如果 `Bridge-2` 不能稳定收敛，直接长跑 `Bridge-3` 风险很高。

因此下一阶段应先让 `Bridge-2` 在 `tol = 1e-6` 下通过，再跑完整 `Bridge-3`。

### 7.3 GMRES 稳但慢

已实现 `sparse-gmres` 作为备选诊断求解器，但当前测试表明：

- GMRES 残差下降较平滑。
- 但在 `Bridge-1-mpc` 上 3000 次仍未达到 `1e-6`。
- 因此不适合作为当前默认主求解器。

## 8. 下一阶段建议

推荐继续稀疏迭代路线，但重点从“是否使用迭代”转为“改进预条件器”。

优先级建议：

1. 实现 block Jacobi 预条件器
   - 按节点自由度块构造。
   - 对每个节点局部自由度做小矩阵求逆或分解。
   - 标准库可实现，迁移风险低。

2. 尝试 ILU(0) 或 IC(0)
   - 若课程规则允许自实现不完整分解，这是更强的方向。
   - 实现复杂度高于 block Jacobi。
   - 需要严格检查零主元、负主元和数值稳定性。

3. 继续使用 `Bridge-2` 作为收敛门槛
   - 验收条件：`CHECKED RELATIVE RESIDUAL <= 1e-6`。
   - 通过后再完整运行 `Bridge-3`。

4. 保留 skyline 小算例对照
   - 用于验证小算例位移一致性。
   - 不建议用于 `Bridge-3` 主求解。

## 9. 当前命令入口

构建：

```powershell
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
```

Bridge-1 / Bridge-2 验证：

```powershell
python tools\run_sparse_validation.py
```

Bridge-3 快速入口：

```powershell
python tools\run_bridge3_fast.py
```

手动运行 Bridge-3：

```powershell
.\build\Release\stap++.exe Bridge-3\Bridge-3.dat --solver sparse-bicgstab --output summary --precond ssor --tol 1e-6 --max-iter 5000 --csv runs\Bridge-3\displacements.csv
```

## 10. 阶段结论

当前阶段已经完成：

- 大算例 summary 输出。
- CSR 稀疏矩阵。
- 稀疏迭代求解器入口。
- MPC alias 转换。
- Bridge-2 无 MPC 读取修复。
- Bridge-1-mpc 收敛验证。
- Bridge-3 读入、alias 和 CSR 装配短跑验证。

当前尚未完成：

- Bridge-2 稳定收敛。
- Bridge-3 完整收敛。

下一阶段的核心问题是预条件器，而不是放宽收敛标准或回退到 skyline。
