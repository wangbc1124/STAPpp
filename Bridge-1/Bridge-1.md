# Bridge-1: Cable-Stayed Bridge — STAP++ 3D Analysis

## Model Overview

| Property | Value |
|----------|-------|
| Input file | `Bridge-1.inp` (Abaqus format) |
| Global nodes | 4163 |
| Global elements | 2884 |
| Equations (STAP++) | 9348 |
| Equations (Abaqus) | 15846 |
| Skyline matrix size | 6,573,903 entries (sparse: 498,250 non-zeros) |

## Coordinate Verification

All 4163 node coordinates have been verified against Abaqus `all_instances_nodes.csv` — **100% match** (tolerance 1e-4).

### Key Fix: Instance Rotation Convention

The Abaqus rotation convention was corrected: Abaqus rotates instances about the **global origin (0,0,0)**, not about the axis point. The two axis points only define the rotation axis direction.

```
# Correct: rotate about origin, then translate
v_rot = R @ v        # R = Rodrigues rotation about origin
global = v_rot + translation

# Incorrect (previous): rotate about axis_a
p = v - axis_a
v_rot = R @ p + axis_a   # Wrong — axis_a is NOT the rotation center
```

## Element Mapping

| Abaqus Type | STAP++ Type | Count | Description |
|------------|------------|-------|-------------|
| T3D2 | Bar (1) | 20 | Steel cables (3D truss) |
| S4R | Q4 (2) | 400 | Deck (plane stress membrane, XY plane) |
| C3D8R | H8 (5) | 1760 | Pier + RiverBank (3D hexahedral) |
| B31 | Bar (1) | 704 | Support beams (converted from 3D beam to truss) |

> **Note**: B31 beams were mapped to Bar (truss) because STAP++ Beam (type 6) only supports 2D frame analysis in the XY plane.

## Materials

| Material | E (GPa) | ν | ρ (kg/m³) | Used in |
|----------|---------|---|-----------|---------|
| Concrete | 25.0 | 0.30 | 2320 | Pier, RiverBank |
| Steel | 117.0 | 0.266 | 7860 | Cables |
| Aluminum | 70.0 | 0.346 | 2710 | SupportBeam |
| Granite | 60.0 | 0.27 | 2770 | RiverBank |

## Boundary Conditions

| Category | Count | Description |
|----------|-------|-------------|
| Abaqus BC Set-102 | 152 nodes | Pier base + RiverBank foundation (all 3 DOFs fixed) |
| Orphan (tie-merged) | 72 nodes | Fully fixed (no element connectivity) |
| Q4-only Z-fixed | 477 nodes | Deck nodes with Z direction fixed |
| Bar-plane fixed | 664 nodes | SupportBeam internal nodes (fully fixed) |
| **Total fully fixed** | **888 nodes** | bc = (1,1,1) |
| **Fully free** | **2798 nodes** | bc = (0,0,0) |

## Tie Constraints

48 tie constraints connect:
- Cable ends → Deck floor nodes (42 ties)
- Cable ends → Pier tower nodes (4 ties)
- SupportBeam ends → RiverBank nodes (4 ties)
- Deck floor ends → RiverBank nodes (10 ties)

72 slave nodes were merged into 34 unique master nodes via transitive closure.

## Geometry (Corrected)

| Feature | Value |
|---------|-------|
| Deck span (X) | -250 ~ 250 m (500m main span) |
| Deck width (Y) | -10 ~ 10 m (20m) |
| Deck elevation | Z = 0 |
| Tower height | 200 m (Z = -50 to 150) |
| Tower-1 center | Y = 10 ~ 20, X = -18 ~ 18 |
| Tower-2 center | Y = -20 ~ -10, X = -18 ~ 18 |
| Cable arrangement | 4 groups × 5 positions (50, 100, 150, 200, 250m) |
| Cable tower end | Z = 150 (tower top) |
| Cable deck end | Z = 0 (deck level) |

## Load

Gravity: magnitude = 10 m/s², direction = (0, 0, -1)
Applied as consistent nodal forces computed from element mass × gravity.

Total gravity force (Z): -5.742 × 10⁹ N (≈ 574,000 tonnes)

## Displacement Results

### Comparison with Abaqus

| Metric | STAP++ | Abaqus | Ratio |
|--------|--------|--------|-------|
| Max |Z| displacement | 0.0159 m | 0.557 m | 35× |
| Location | Node 513 (tower top) | Node 198 (Part-Floor-1) | — |
| Max |Y| displacement | 0.0093 m | — | — |
| Max |X| displacement | 0.0009 m | — | — |
| Equation count | 9348 | 15846 | 0.59× |
| Skyline matrix size (NWK) | 6,573,903 | — | — |
| Half bandwidth (MK) | 7304 | — | — |

### STAP++ Detailed Displacements — Top 10 |dz|

| Node | dx (m) | dy (m) | dz (m) |
|------|--------|--------|--------|
| 513 | +5.84e-5 | +9.32e-3 | **-1.588e-2** |
| 1371 | +5.84e-5 | -9.32e-3 | -1.588e-2 |
| 692 | +5.66e-5 | +9.27e-3 | -1.567e-2 |
| 1555 | +5.66e-5 | -9.27e-3 | -1.567e-2 |
| 528 | +5.97e-5 | +9.27e-3 | -1.566e-2 |
| 1387 | +5.97e-5 | -9.27e-3 | -1.566e-2 |
| 569 | +5.54e-5 | +8.71e-3 | -1.566e-2 |
| 1546 | +5.54e-5 | -8.71e-3 | -1.566e-2 |
| 1170 | +8.46e-5 | +8.72e-3 | -1.565e-2 |
| 1953 | +8.46e-5 | -8.72e-3 | -1.565e-2 |

### Pier Tower Displacements (typical)

| Node | dx (m) | dy (m) | dz (m) |
|------|--------|--------|--------|
| 506 (base) | 0.0 | 0.0 | 0.0 |
| 804 (mid) | -1.56e-4 | +1.04e-3 | -8.40e-3 |
| 510 (top) | +5.59e-5 | +9.24e-3 | -1.46e-2 |

### Analysis of Displacement Differences

The STAP++ max Z displacement (1.59 cm) is approximately **35× smaller** than Abaqus (55.7 cm). This large difference is primarily due to element formulation differences:

1. **Deck (S4R → Q4)**: Abaqus S4R is a 4-node **shell** with full bending stiffness (6 DOFs/node). STAP++ Q4 is a **plane stress membrane** (2 DOFs/node, no bending). Under gravity (out-of-plane loading), the shell bends significantly while the membrane is stiff in-plane but has no out-of-plane deformation mode. 477 deck nodes have Z fixed, further stiffening the deck.

2. **Support beams (B31 → Bar)**: Abaqus B31 is a 3D Timoshenko beam with bending, shear, and axial stiffness (6 DOFs/node). STAP++ Bar is a truss with only axial stiffness (3 DOFs/node). 664 Bar-only nodes are fully fixed to prevent mechanism formation, removing their flexibility contribution entirely.

3. **H8 full integration locking**: STAP++ H8 uses full 2×2×2 Gauss integration, which can exhibit shear and volumetric locking for thin sections. Abaqus C3D8R uses reduced integration with hourglass control, making it more flexible and accurate for bending-dominated problems.

4. **DOF count**: Abaqus has 15846 equations vs STAP++ 9348, reflecting the richer element formulations (shells have 6 DOFs, beams have 6 DOFs, vs 2 for Q4 and 3 for Bar).

## Stress Results

### STAP++ Element Stresses

| Group | Type | Count | Max Force | Max Stress | Notes |
|-------|------|-------|-----------|------------|-------|
| 1 | Bar (SupportBeam) | 704 | 4.21 MN | 5.54 MPa | 46 of 704 elements have non-zero result |
| 2 | Bar (Cables) | 20 | 2.70 MN | 10.8 MPa | All 20 cables in tension |
| 3 | Q4 (Deck) | 400 | — | 0.66 MPa (von Mises) | Membrane stress, low |
| 4 | H8 (Pier-2) | 960 | — | 3.34 MPa (von Mises) | σ_z dominant (compression) |
| 5 | H8 (Pier-1+RiverBank) | 800 | — | 1.23 MPa (von Mises) | Low stress, σ_z dominant |

### Stress Analysis

- **Cables**: Maximum force of 2.70 MN in steel cables. At A=0.25 m², stress = 10.8 MPa, well below steel yield strength (~250 MPa). All 20 cables are in tension as expected for a cable-stayed bridge under gravity.

- **Deck**: Maximum von Mises stress of 0.66 MPa. Very low — the Q4 membrane doesn't capture bending stresses. This is expected given the plane stress formulation.

- **Pier (H8)**: Maximum von Mises stress of ~3.34 MPa in Pier-2 and ~1.23 MPa in Pier-1+RiverBank, well within concrete capacity (25–40 MPa). Stresses are predominantly compressive in the Z direction (gravity-aligned).

- **SupportBeam bars**: 658 of 704 elements (93%) show zero force because their nodes are fully fixed to prevent mechanism formation (B31→Bar conversion loses bending stiffness). The 46 active elements carry forces up to 4.21 MN.

## Abaqus 对比分析

从 `Bridge-1.odb` 中提取了完整场输出数据（使用 Abaqus Python API），与 STAP++ C++ 求解器结果进行逐项对比。ODB 数据包含 4163 个节点的位移（U）和 5396 个单元积分点的应力（S, TENSOR_3D_SURFACE），以及 4163 个节点的反力（RF）。

### 数据来源

| 数据 | 来源 | 条目数 |
|------|------|--------|
| Abaqus 位移 | `Bridge-1.odb` → Step-1, Frame 1, Field U | 4163 |
| Abaqus 应力 | `Bridge-1.odb` → Step-1, Frame 1, Field S | 5396 |
| Abaqus 反力 | `Bridge-1.odb` → Step-1, Frame 1, Field RF | 4163 |
| STAP++ 位移 | `Bridge-1-cpp.out` | 4163 |
| STAP++ 应力 | `Bridge-1-cpp.out` (5 个单元组) | 2884 |

节点映射通过坐标匹配完成（容差 0.005m）。在 4163 个 Abaqus 节点中，**3235 个可唯一映射**到 STAP++ 全局节点，**928 个因绑定约束合并而共享坐标**（Abaqus 保留为独立节点，STAP++ 合并为同一节点）。

### 整体位移对比

| 指标 | Abaqus | STAP++ | 比值 S/A |
|------|--------|--------|----------|
| Max \|DZ\| (m) | 2.973×10⁻¹ | 1.588×10⁻² | 0.053 |
| Max \|DY\| (m) | 2.428×10⁻¹ | 9.320×10⁻³ | 0.038 |
| Max \|DX\| (m) | 3.294×10⁻² | 7.344×10⁻⁴ | 0.022 |
| 方程数 | 15,846 | 9,348 | 0.590 |
| 自由度/节点 | 1~6 (取决于单元) | 2~3 (取决于单元) | — |

> **核心发现**：STAP++ 整体位移幅值约为 Abaqus 的 **2%–5%**。最大 Z 向位移 STAP++ 仅 1.59 cm（塔顶），而 Abaqus 达 29.7 cm（支撑梁）。

### 各部件位移对比

| 部件（实例） | 对比节点数 | Abaqus Max\|DZ\| | STAP++ Max\|DZ\| | 比值 | 评估 |
|-------------|-----------|-----------------|-----------------|------|------|
| PART-CABLE50 | 8 | 0.244 m | 0.0159 m | 0.065 | 显著偏小 |
| PART-CABLE100 | 8 | 0.235 m | 0.0159 m | 0.068 | 显著偏小 |
| PART-CABLE150 | 8 | 0.216 m | 0.0159 m | 0.073 | 显著偏小 |
| PART-CABLE200 | 8 | 0.138 m | 0.0159 m | 0.115 | 偏小 |
| PART-CABLE250 | 8 | 0.087 m | 0.0159 m | 0.183 | 偏小 |
| **PART-FLOOR-1（甲板）** | 28 | **0.294 m** | **0.0063 m** | **0.021** | **显著偏小** |
| PART-PIER-1（桥塔1） | 861 | 0.087 m | 0.0159 m | 0.183 | 偏小 |
| PART-PIER-2（桥塔2） | 410 | 0.087 m | 0.0156 m | 0.180 | 偏小 |
| PART-RIVERBANK-1（河岸1） | 605 | 0.0023 m | 0.0006 m | 0.271 | 偏小 |
| PART-RIVERBANK-2（河岸2） | 605 | 0.0023 m | 0.0006 m | 0.271 | 偏小 |
| PART-SUPPORTBEAM-1（支撑梁1） | 343 | 0.297 m | 0.0003 m | 0.001 | 极大偏小 |
| PART-SUPPORTBEAM-2（支撑梁2） | 343 | 0.297 m | 0.0063 m | 0.021 | 显著偏小 |

> **关键说明**：上表中"对比节点数"是坐标可唯一匹配的节点数，不等于部件实际节点数。Abaqus 中多个实例在相同位置有独立节点（如索端部与甲板/塔节点重合），STAP++ 通过绑定约束将其合并。绳索节点 Z 向位移为 0 是因为其端点与甲板或塔的固定节点绑定。

### 甲板位移详细对比

甲板是差异最大的部件。Abaqus 壳单元（S4R）具有完整的弯曲刚度，在重力（面外载荷）作用下产生 0–0.29 m 的竖向挠度。STAP++ 的 Q4 膜单元只有面内刚度，Z 向必须固定以防止奇异性。

| Abaqus 节点 | 坐标 (X, Y) | Abaqus dz (m) | STAP++ dz (m) | 差值 (m) |
|------------|-------------|---------------|---------------|----------|
| 9 | (0, 10) | -0.2944 | -0.0063 | -0.2881 |
| 10 | (0, -10) | -0.2944 | 0.0 | -0.2944 |
| 7 | (50, 10) | -0.2443 | 0.0 | -0.2443 |
| 8 | (50, -10) | -0.2443 | 0.0 | -0.2443 |
| 11 | (-50, 10) | -0.2443 | 0.0 | -0.2443 |
| 12 | (-50, -10) | -0.2443 | 0.0 | -0.2443 |
| 5 | (100, 10) | -0.2350 | 0.0 | -0.2350 |
| 6 | (100, -10) | -0.2350 | 0.0 | -0.2350 |
| 13 | (-100, 10) | -0.2350 | 0.0 | -0.2350 |
| 14 | (-100, -10) | -0.2350 | 0.0 | -0.2350 |

> **结论**：甲板柔度完全缺失——Abaqus 中跨中挠度 0.29 m（L/1700），符合斜拉桥刚度预期。STAP++ 甲板基本无竖向变形。

### 桥塔位移对比

选取桥塔三个典型标高进行对比：

| 位置 | 标高 | Abaqus dz | STAP++ dz | 比值 S/A | Abaqus dy | STAP++ dy |
|------|------|-----------|-----------|----------|-----------|-----------|
| 塔底 | Z ≈ -50 m | ~0 | 0 | — | ~0 | 0 |
| 塔中 | Z ≈ 20 m | -0.0054 | -0.0080 | **1.47** | -0.0238 | -0.0007 |
| 塔顶 | Z ≈ 150 m | -0.0868 | -0.0159 | 0.183 | -0.2428 | -0.0093 |

> **关键发现**：塔顶 Y 向位移差距达 26 倍（Abaqus -0.24 m vs STAP++ -0.009 m），反映索力通过甲板弯曲传递的侧向分量在 STAP++ 中严重不足。

### 支撑梁位移对比

支撑梁是差异最大的部件：Abaqus B31 梁单元在重力下有 **0.30 m** 的竖向挠度，而 STAP++ Bar 桁架单元节点全部固定（防止机构），**位移为零**。以下是位移差最大的节点：

| STAP++ 节点 | 坐标 | Abaqus dz (m) | STAP++ dz (m) |
|------------|------|---------------|---------------|
| 3581 | (-13.3, -10.0, -8.0) | -0.2973 | 0.0 |
| 3595 | (13.3, -10.0, -8.0) | -0.2973 | 0.0 |
| 3924 | (-13.3, 10.0, -8.0) | -0.2973 | 0.0 |
| 3938 | (13.3, 10.0, -8.0) | -0.2973 | 0.0 |

### 应力对比

#### Abaqus 各部件应力（ODB 提取）

| 部件 | 应力记录数 | Max von Mises (MPa) | Avg von Mises (MPa) |
|------|-----------|--------------------|--------------------|
| PART-CABLE50 | 4 | 109.1 | 109.1 |
| PART-CABLE100 | 4 | 77.2 | 77.2 |
| PART-CABLE150 | 4 | 46.1 | 46.1 |
| PART-CABLE200 | 4 | 6.84 | 6.84 |
| PART-CABLE250 | 4 | 17.7 | 17.7 |
| **PART-FLOOR-1（甲板）** | 800 | **39.9** | 15.5 |
| PART-PIER-1（桥塔1） | 480 | 9.16 | 2.47 |
| PART-PIER-2（桥塔2） | 480 | 9.16 | 2.47 |
| PART-RIVERBANK-1 | 400 | 3.39 | 0.83 |
| PART-RIVERBANK-2 | 400 | 3.39 | 0.83 |
| PART-SUPPORTBEAM-1 | 1408 | 69.4 | 12.5 |
| PART-SUPPORTBEAM-2 | 1408 | 69.4 | 12.5 |

#### STAP++ 各单元组应力

| 单元组 | 类型 | 单元数 | Max von Mises (MPa) | Avg von Mises (MPa) | 说明 |
|--------|------|--------|--------------------|--------------------|------|
| 1 | Bar（支撑梁） | 704 | 5.54 | — | 46/704 非零，其余全固定 |
| 2 | Bar（索） | 20 | 10.8 | — | 全部 20 根受拉 |
| 3 | Q4（甲板） | 400 | 0.66 | 0.22 | 仅膜应力，无弯曲分量 |
| 4 | H8（桥塔2） | 960 | 3.34 | 1.88 | 全积分，可能锁死 |
| 5 | H8（桥塔1+河岸） | 800 | 1.23 | 0.66 | 全积分，可能锁死 |

#### 应力差异分析

| 对比项 | Abaqus | STAP++ | 比值 | 原因 |
|--------|--------|--------|------|------|
| 甲板 Max VM | 39.9 MPa | 0.66 MPa | **0.017** | Q4无弯曲应力，仅膜应力 |
| 桥塔 Max VM | 9.16 MPa | 3.34 MPa | 0.365 | H8全积分锁死 vs C3D8R减积分 |
| 索 Max 应力 | 109.1 MPa | 10.8 MPa | 0.099 | 索端与甲板/塔绑定，甲板无弯曲导致索力偏低 |
| 支撑梁 Max VM | 69.4 MPa | 5.54 MPa | 0.080 | Bar无弯曲，节点固定 |

> **核心发现**：所有部件的 STAP++ 应力均低于 Abaqus。甲板应力差异最大（60 倍），因为 Q4 膜单元完全丢失了弯曲应力分量——这是重力作用下甲板的主导应力形式。

### 索力对比

STAP++ 索力（Group 2，20 根 Bar 单元）全部为拉力，符合斜拉桥受力特征：

- 最大索力：**2.70 MN**（Elem 1, 最短索 CABLE50）
- 最小索力：**0.77 MN**（Elem 18, 最长索 CABLE250）
- 索应力范围：3.1–10.8 MPa（远低于钢索屈服强度 ~250 MPa）

Abaqus ODB 中索的应力数据每个实例仅 1 个记录（积分点外推值），参考价值有限。但 Abaqus 索 von Mises 应力（CABLE50: 109 MPa）明显高于 STAP++（10.8 MPa），原因同上——Abaqus 中甲板弯曲带动索端产生更大位移差，从而增大索力。

### 差异根因总结

STAP++ 与 Abaqus 之间的定量差异源于**单元列式的本质限制**，并非模型或求解器错误：

1. **Q4 膜单元 ≠ S4R 壳单元**（最大影响）：
   - Q4 仅 2 个平动自由度/节点（X, Y），无面外刚度，无弯曲
   - S4R 有 6 个自由度/节点，含旋转自由度，完全捕捉弯曲
   - 重力是面外载荷 → 弯曲是甲板主导变形模式 → Q4 完全丢失

2. **Bar 桁架单元 ≠ B31 铁木辛柯梁**：
   - Bar 仅有轴向刚度（3 DOFs），无弯曲/剪切/扭转
   - B31 含全部 6 个 DOFs
   - 664 个 Bar 独有节点被全固定以防止机构，丧失全部柔度贡献

3. **H8 全积分 ≠ C3D8R 减积分**：
   - 全积分（2×2×2）在薄壁截面可产生剪切/体积锁死
   - 减积分（1 个积分点 + 沙漏控制）更柔、更准确
   - 桥塔应力偏低 2.7 倍反映了锁死效应

4. **绑定约束节点合并**：
   - Abaqus 保留独立节点（通过约束方程耦合）
   - STAP++ 直接合并为同一节点
   - 928 个节点无法唯一映射

5. **自由度总数差距**：
   - Abaqus 15,846 vs STAP++ 9,348（1.7 倍）
   - 反映了壳和梁的旋转自由度以及保留的独立节点

### 数据可用性说明

当前 ODB 中应力类型为 `TENSOR_3D_SURFACE`（3 分量面应力），这意味着 ODB 仅包含壳/膜单元的应力输出。对于实体单元（C3D8R），ODB 中可能存储为单独的输出请求或默认未输出完整 6 分量应力。更全面的对比需要重新运行 Abaqus 作业并输出 `TENSOR_3D_FULL` 应力分量。

## 求解器性能

| 指标 | 数值 |
|------|------|
| 矩阵维度 (NEQ) | 9348 |
| 天际线存储条目 (NWK) | 6,573,903 |
| 最大半带宽 (MK) | 7304 |
| 最大初始对角元 | 5.43 × 10¹¹ |
| 最小主元（分解后） | 3.55 × 10⁷ (eq. 6047) |
| 主元比率 | 6.53 × 10⁻⁵ |
| 刚度矩阵组装时间 | 0.046 s |
| LDLT 分解时间 | 4.386 s |
| 输入解析时间 | 0.151 s |
| **总求解时间** | **4.583 s** |
| 主元容差 | 1.0 × 10⁻¹⁷ × max_diag = 5.43 × 10⁻⁶ |
| 结果 | 矩阵正定，求解收敛 |

## 求解器实现

STAP++ C++ 求解器从源码（`src/cpp/*.cpp`）使用 clang++ 编译：

```
clang++ -O2 -std=c++11 -Isrc/h -o Bridge-1/stap.exe src/cpp/*.cpp
```

**编译详情：**
- 编译器：clang++ 22.1.4 (MSYS2 UCRT64)
- 优化级别：`-O2`
- C++ 标准：C++11
- 二进制大小：264 KB
- 注：g++ 16.1.0 在当前 MSYS2 环境中存在问题（cc1plus.exe 崩溃），选用 clang++ 作为替代方案

**求解流程：**
1. `convert_inp_to_dat.py` — 解析 Abaqus `.inp`，施加三维实例变换、绑定约束、边界条件和重力载荷，输出 STAP++ `.dat` 格式
2. `stap.exe` — C++ 天际线 LDLT 求解器，读取 `.dat`，组装全局刚度矩阵，求解 Ku=F
3. 输出：全部 4163 个节点的位移和全部 2884 个单元的应力（24,336 行）

**支持的单元类型：** Bar(1)、Q4(2)、H8(5)，使用完全 Gauss 积分（Q4: 2×2, H8: 2×2×2）

## 关键观测

1. **坐标修正至关重要**：Abaqus 的旋转约定是绕全局原点旋转再平移，而非绕轴点旋转。全部 4163 个节点现已与 Abaqus 精确匹配。

2. **C++ 求解器运行正确**：NEQ=9348, NWK=6,573,903, MK=7304。LDLT 分解耗时 4.386s，主元比率良好（6.53×10⁻⁵），矩阵正定。总耗时 4.583s。

3. **单元列式限制是位移差异的根本原因**：STAP++ 整体位移约为 Abaqus 的 2%–5%，详见 Abaqus 对比章节。核心差异为：
   - Q4 膜单元（无弯曲） vs S4R 壳单元（含弯曲）
   - Bar 桁架（无弯曲/剪切） vs B31 梁（完整 6-DOF）
   - H8 全积分（可能锁死） vs C3D8R 减积分

4. **修正后的几何模型更刚硬**：塔顶 Z=150m（此前错误为 Z=244m），索更短，柔度更小。

5. **绑定约束处理正确**：传递闭包正确合并了约束链。72 个从节点 → 34 个主节点。但这也导致 928 个 Abaqus 节点无法唯一映射到 STAP++ 节点（共享坐标但 Abaqus 保留为独立节点）。

6. **Bar 节点固定是保真度损失主因**：B31→Bar 转换导致 664 个节点必须全固定以避免刚度奇异。这是与 Abaqus 结果差异的最主要原因，支撑梁在 STAP++ 中基本无贡献。

## 局限性

1. **无壳弯曲**：Q4 膜单元无法表示桥面板弯曲。Abaqus S4R 壳单元捕捉弯曲——这是重力载荷下的主导变形模式。

2. **无梁弯曲**：B31→Bar 转换消除了支撑梁的弯曲和剪切刚度。

3. **H8 全积分**：薄壁桥塔截面可能出现剪切/体积锁死。Abaqus C3D8R 使用减积分 + 沙漏控制。

4. **无预应力**：未模拟索的预应力效应（实际斜拉桥的关键受力特征）。

5. **仅线性静力**：无几何或材料非线性。

6. **全固定的 Bar 节点**：704 个支撑梁单元中 658 个贡献为零（节点全固定以防止刚度奇异）。

7. **Abaqus 应力数据不完整**：ODB 中应力为 TENSOR_3D_SURFACE（3 分量），缺少实体单元的完整 6 分量应力输出。更全面的对比需要重新运行 Abaqus 作业。

---

*生成日期：2026-06-09 — STAP++ C++ 求解器 (clang++ 22.1.4), Bridge-1 三维分析, 含 Abaqus ODB 全量对比*
