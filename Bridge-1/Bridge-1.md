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

## Abaqus Comparison

Abaqus reference results were extracted from `abaqus.rpt` (generated from `Bridge-1.odb`). The report contains displacement (U) and von Mises stress (S, Mises) for 8 sample nodes on the deck (Part-Floor-1). These nodes exist in both models with identical IDs and coordinates.

### Deck Node Displacement Comparison

All STAP++ deck (Q4) nodes have Z direction fixed (`bc = (0,0,1)`) because the Q4 plane stress membrane has only 2 DOFs per node (X, Y). The Z DOF must be fixed to prevent a singular stiffness matrix. Abaqus S4R shells have 6 DOFs per node and capture the full out-of-plane bending deformation.

| Abaqus Node | Position (X, Y) | Abaqus dz (m) | STAP++ dz (m) | Ratio |
|-------------|-----------------|---------------|---------------|-------|
| 271 | (135, 5) | -0.413 | 0.000 | 0 |
| 249 | (175, 0) | -0.347 | 0.000 | 0 |
| 492 | (225, 0) | -0.536 | 0.000 | 0 |
| 297 | (85, 0) | -0.403 | 0.000 | 0 |
| 79 | (50, -5) | -0.296 | 0.000 | 0 |
| 329 | (25, -5) | -0.499 | 0.000 | 0 |
| 341 | (5, -5) | -0.362 | 0.000 | 0 |
| 346 | (-5, 5) | -0.362 | 0.000 | 0 |

> **Key finding**: STAP++ Q4 deck nodes have exactly zero Z displacement because their Z DOF is fixed. Abaqus S4R shells show 0.3–0.5 m of gravity-induced bending deformation, which the Q4 membrane fundamentally cannot represent. The deck's primary deformation mode under vertical gravity load is out-of-plane bending — a mode the Q4 element cannot capture.

### Deck Node Stress Comparison

Abaqus reports von Mises stress at the same 8 deck nodes (extrapolated from integration points). STAP++ stresses are computed at Gauss points and averaged to nodes from connected Q4 elements.

| Abaqus Node | Connected Q4 Elements | STAP++ VM (MPa) | Abaqus VM (MPa) | Ratio STAP++/Abaqus |
|-------------|----------------------|-----------------|-----------------|---------------------|
| 271 | 4 | 0.241 | 9 | 0.027 |
| 249 | 4 | 0.284 | 11 | 0.026 |
| 492 | 4 | 0.077 | 29 | 0.003 |
| 297 | 4 | 0.151 | 7 | 0.022 |
| 79 | 4 | 0.575 | 28 | 0.021 |
| 329 | 4 | 0.275 | 12 | 0.023 |
| 341 | 4 | 0.273 | 21 | 0.013 |
| 346 | 4 | 0.272 | 21 | 0.013 |

> **Key finding**: STAP++ Q4 von Mises stresses are **13× to 370× smaller** than Abaqus S4R stresses. The Abaqus shell captures both membrane and bending stress components, while the Q4 membrane only captures in-plane (membrane) stresses. For a deck under gravity (out-of-plane load), bending stress dominates — typically 10–50× larger than membrane stress.

### Overall Comparison Summary

| Metric | STAP++ | Abaqus | Ratio |
|--------|--------|--------|-------|
| Equation count | 9,348 | 15,846 | 0.59× |
| Max deck |Z| displacement | 0 m (fixed) | 0.54 m | 0 |
| Max tower top |Z| displacement | 0.016 m | — | — |
| Deck von Mises stress range | 0.08–0.66 MPa | 7–29 MPa | ~0.02× |
| Cable forces | 0.8–2.7 MN (all tension) | — | — |
| Pier von Mises stress | 1.2–3.3 MPa | — | — |

### Root Cause of Differences

The large quantitative difference between STAP++ and Abaqus results stems from **element formulation limitations**, not from errors in the model or solver:

1. **Q4 ≠ S4R** (dominant factor): The deck under gravity is a bending-dominated problem. Abaqus S4R shells capture bending through rotational DOFs and curved element formulation. STAP++ Q4 is a plane stress membrane — it only resists in-plane stretching, not bending. This is the single largest source of discrepancy.

2. **Z-fixed deck nodes**: Because Q4 has no Z-stiffness (only 2 DOFs: X, Y), 477 deck nodes have their Z DOF fixed to prevent singular stiffness. This entirely eliminates deck vertical displacement in STAP++, while Abaqus shows the deck bending downward by 0.3–0.5 m.

3. **Bar ≠ B31**: B31 Timoshenko beams include bending, shear, and torsional stiffness (6 DOFs). Bar truss elements have only axial stiffness (3 DOFs). 664 Bar-only nodes are fully fixed to prevent mechanisms.

4. **H8 vs C3D8R**: Full integration H8 may exhibit locking in thin sections, while C3D8R uses reduced integration with hourglass control.

### Data Availability

The Abaqus `.rpt` contains data for only 8 deck nodes. No Abaqus data is available for:
- Pier (H8) displacements or stresses
- Cable (T3D2) forces or stresses
- SupportBeam (B31) displacements or stresses
- Tower top (node 513) displacement

A more complete comparison would require extracting additional field output from the Abaqus `.odb` file.

## Solver Performance

| Metric | Value |
|--------|-------|
| Matrix dimension (NEQ) | 9348 |
| Skyline entries (NWK) | 6,573,903 |
| Maximum half bandwidth (MK) | 7304 |
| Max initial diagonal | 5.43 × 10¹¹ |
| Min pivot (after factorization) | 3.55 × 10⁷ (at eq. 6047) |
| Pivot ratio | 6.53 × 10⁻⁵ |
| Stiffness assembly time | 0.046 s |
| LDLT factorization time | 4.386 s |
| Input phase time | 0.151 s |
| **Total solution time** | **4.583 s** |
| Pivot tolerance | 1.0 × 10⁻¹⁷ × max_diag = 5.43 × 10⁻⁶ |
| Result | Matrix is positive definite, solution converged |

## Solver Implementation

The STAP++ C++ solver was compiled from source (`src/cpp/*.cpp`) using clang++:

```
clang++ -O2 -std=c++11 -Isrc/h -o Bridge-1/stap.exe src/cpp/*.cpp
```

**Compilation details:**
- Compiler: clang++ 22.1.4 (MSYS2 UCRT64)
- Optimization: `-O2`
- Standard: C++11
- Binary size: 264 KB
- Note: g++ 16.1.0 was broken in the current MSYS2 environment (cc1plus.exe crash), so clang++ was used as the workaround.

**Solver pipeline:**
1. `convert_inp_to_dat.py` — Parses Abaqus `.inp`, applies 3D instance transformations, tie constraints, BCs, and gravity loads, outputs STAP++ `.dat` format
2. `stap.exe` — C++ skyline LDLT solver, reads `.dat`, assembles global stiffness, solves Ku=F
3. Output: displacements for all 4163 nodes and stresses for all 2884 elements (24,336 lines)

**Supported element types:** Bar(1), Q4(2), H8(5) with full Gauss integration (2×2 for Q4, 2×2×2 for H8).

## Key Observations

1. **Coordinate fix is critical**: The Abaqus rotation convention (rotate about origin, then translate) is different from the intuitive "rotate about axis point" interpretation. All 4163 nodes now match Abaqus exactly.

2. **C++ solver works correctly**: NEQ=9348, NWK=6,573,903, MK=7304. LDLT factorization completes in 4.386s. Pivot ratio is healthy (6.53×10⁻⁵). Matrix is positive definite with minimum pivot 3.55×10⁷ at equation 6047. Converged solution in 4.583s total. Compiled with clang++ 22.1.4 (g++ 16.1.0 was broken in the current MSYS2 environment).

3. **Element limitations drive the displacement difference**: The 35× ratio vs Abaqus is expected given:
   - Q4 membrane (no bending) vs S4R shell (with bending)
   - Bar truss (no bending/shear) vs B31 beam (full 6-DOF)
   - H8 full integration (potential locking) vs C3D8R reduced integration

4. **The corrected geometry produces a stiffer structure**: With the tower top at Z=150 (vs previously Z=244), the cables are shorter and provide less flexibility. The previous wrong coordinates accidentally created a more flexible structure closer to Abaqus results.

5. **Tie constraint handling**: The transitive closure correctly merges chains of tied nodes. 72 slave nodes → 34 master nodes.

6. **Bar-only node stabilization**: SupportBeam B31→Bar conversion creates mechanisms. 664 Bar-only nodes are fully fixed. This is the primary loss of fidelity compared to Abaqus.

## Limitations

1. **No shell bending**: Q4 membrane cannot represent deck plate bending. Abaqus S4R shells capture this, which is the dominant deformation mode under gravity.

2. **No beam bending**: B31→Bar conversion eliminates support beam bending and shear stiffness.

3. **Full integration H8**: May exhibit shear/volumetric locking for thin pier sections. Abaqus C3D8R uses reduced integration with hourglass control.

4. **No prestress**: Cable prestress effects are not modeled (significant in real cable-stayed bridges).

5. **Linear static only**: No geometric or material nonlinearity.

6. **Fully fixed Bar nodes**: 658/704 SupportBeam elements contribute zero to the solution because their nodes are fixed to prevent singular stiffness.

---

*Generated: June 3, 2026 — STAP++ C++ Solver (clang++ 22.1.4), Bridge-1 3D Analysis*
