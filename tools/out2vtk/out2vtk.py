#!/usr/bin/env python3
"""
STAP++ .out to ParaView .vtk converter

Converts finite element output from STAP++ (.out) to legacy VTK unstructured
grid format for visualization in ParaView.

Usage: python out2vtk.py <input.out> [output.vtk]

Supports all element types:
  1  Truss  (2-node line)
  2  Q4     (4-node quad, full integration)
  3  Q4R    (4-node quad, reduced integration)
  4  T3     (3-node triangle, CST)
  5  H8     (8-node hexahedron)
  6  Beam   (2-node Euler-Bernoulli frame)
  7  Plate  (4-node Mindlin-Reissner plate bending)
  8  Shell  (4-node shell)
  9  Beam3D (2-node 3D beam)
  10 Shell4 (4-node flat shell, membrane + plate bending)
  11 Beam3DTimoshenko (2-node B31-equivalent spatial beam)
  12 H8R    (8-node reduced-integration hexahedron)
  14 H8RPier (8-node pier-local reduced-integration hexahedron)
"""

import re
import sys
import os
import math

# ── VTK cell type constants ───────────────────────────────────────────
VTK_LINE       = 3
VTK_TRIANGLE   = 5
VTK_QUAD       = 9
VTK_HEXAHEDRON = 12

# ── Element type definitions ──────────────────────────────────────────
# elem_type -> (vtk_cell_type, nodes_per_elem, stress_names, data_lines_per_elem)
ELEM_DEF = {
    1:  (VTK_LINE,        2, ['Force',  'Stress'],  1),
    2:  (VTK_QUAD,        4, ['Sigma_X','Sigma_Y','Tau_XY'],  1),
    3:  (VTK_QUAD,        4, ['Sigma_X','Sigma_Y','Tau_XY'],  1),
    4:  (VTK_TRIANGLE,    3, ['Sigma_X','Sigma_Y','Tau_XY'],  1),
    5:  (VTK_HEXAHEDRON,  8, ['Sigma_X','Sigma_Y','Sigma_Z',
                              'Tau_XY', 'Tau_YZ', 'Tau_ZX'],  2),
    6:  (VTK_LINE,        2, ['Axial_Stress','Bend_Stress1','Bend_Stress2'], 1),
    7:  (VTK_QUAD,        4, ['Mx','My','Mxy','Qx','Qy'], 2),
    8:  (VTK_QUAD,        4, ['Sigma_X','Sigma_Y','Tau_XY',
                              'Mx','My','Mxy','Qx','Qy'],  3),
    9:  (VTK_LINE,        2, ['Axial_Force','Moment_Y1','Moment_Z1',
                              'Torque1','Axial_Stress','Moment_Y2',
                              'Moment_Z2','Torque2'],  3),
    10: (VTK_QUAD,        4, ['Sigma_X','Sigma_Y','Tau_XY',
                              'Mx','My','Mxy','Qx','Qy'],  3),
    11: (VTK_LINE,        2, ['Axial_Force','Moment_Y1','Moment_Z1',
                              'Torque1','Axial_Stress','Moment_Y2',
                              'Moment_Z2','Torque2'],  3),
    12: (VTK_HEXAHEDRON,  8, ['Sigma_X','Sigma_Y','Sigma_Z',
                              'Tau_XY', 'Tau_YZ', 'Tau_ZX'],  2),
    14: (VTK_HEXAHEDRON,  8, ['Sigma_X','Sigma_Y','Sigma_Z',
                              'Tau_XY', 'Tau_YZ', 'Tau_ZX'],  2),
}


def _find_numbers(line):
    """Return all floating-point numbers found in a line as float list."""
    return [float(v) for v in re.findall(r'[-+]?\d+\.?\d*(?:[eE][-+]?\d+)?', line)]


def _has_kw(line, *keywords):
    """Check if line contains all space-separated keyword groups,
       ignoring the exact number of spaces between letters/words."""
    nospace = line.replace(' ', '')
    for kw in keywords:
        if kw.replace(' ', '') not in nospace:
            return False
    return True


def parse_out(filepath):
    """Parse a STAP++ .out file and return structured data dict."""
    with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
        lines = f.readlines()

    result = {
        'title': '',
        'num_nodes': 0,
        'num_elem_groups': 0,
        'num_load_cases': 0,
        'nodes': [],
        'elements': [],
        'elem_group_types': [],
        'elem_type_by_group': [],  # one entry per group (may have duplicates)
        'elem_count_by_group': [],
        'displacements': [],
        'stresses': [],
    }

    current_elem_type = None
    i = 0
    while i < len(lines):
        line = lines[i]

        if line.startswith('TITLE :'):
            result['title'] = line.split(':', 1)[1].strip()

        elif 'NUMNP' in line and '=' in line:
            m = re.search(r'=\s*(\d+)', line)
            if m:
                result['num_nodes'] = int(m.group(1))

        elif 'NUMEG' in line and '=' in line:
            m = re.search(r'=\s*(\d+)', line)
            if m:
                result['num_elem_groups'] = int(m.group(1))

        elif 'NLCASE' in line and '=' in line:
            m = re.search(r'=\s*(\d+)', line)
            if m:
                result['num_load_cases'] = int(m.group(1))

        elif _has_kw(line, 'NODAL', 'POINT', 'DATA'):
            i = _parse_nodal_data(lines, i, result)

        elif _has_kw(line, 'ELEMENT', 'GROUP', 'DATA'):
            i, current_elem_type = _parse_element_group(lines, i, result, current_elem_type)

        elif 'ELEMENT TYPE' in line and 'NPAR(1)' in line:
            # Files that skip the "ELEMENT GROUP DATA" header line
            i, current_elem_type = _parse_element_group(lines, i - 1, result, current_elem_type)

        elif _has_kw(line, 'DISPLACEMENTS'):
            i = _parse_displacements(lines, i, result)

        elif _has_kw(line, 'STRESS', 'CALCULATIONS'):
            i = _parse_stresses(lines, i, result)

        i += 1

    return result


def _parse_nodal_data(lines, start, result):
    """Parse node coordinates from NODAL POINT DATA section."""
    i = start + 1
    while i < len(lines):
        line = lines[i]
        if 'EQUATION NUMBERS' in line:
            break
        nums = _find_numbers(line)
        if len(nums) >= 7:
            result['nodes'].append((nums[-3], nums[-2], nums[-1]))
        i += 1
    return i


def _parse_element_group(lines, start, result, prev_elem_type):
    """Parse a single element group: type, count, materials, connectivity.
    Returns (next_i, elem_type).  If next_i points to a line that was not
    consumed (next-group header), the caller's i += 1 will land back on it."""
    i = start + 1
    elem_type = None
    num_elems = 0

    while i < len(lines):
        line = lines[i]

        if 'TOTAL SYSTEM DATA' in line:
            break
        if _has_kw(line, 'ELEMENT', 'GROUP', 'DATA') and i > start + 1:
            break
        # Multi-group files may skip GROUP DATA header between groups.
        # Rewind i by 1 so main-loop's i += 1 lands on this line again.
        if 'ELEMENT TYPE' in line and 'NPAR(1)' in line and elem_type is not None:
            i -= 1
            break

        if 'NPAR(1)' in line:
            m = re.search(r'=\s*(\d+)', line)
            if m:
                elem_type = int(m.group(1))
                if elem_type not in result['elem_group_types']:
                    result['elem_group_types'].append(elem_type)
                result['elem_type_by_group'].append(elem_type)

        if 'NPAR(2)' in line:
            m = re.search(r'=\s*(\d+)', line)
            if m:
                num_elems = int(m.group(1))

        if _has_kw(line, 'ELEMENT', 'INFORMATION') and elem_type is not None:
            i, parsed_count = _parse_element_info(lines, i, elem_type, num_elems, result)
            result['elem_count_by_group'].append(parsed_count)

        i += 1

    return i, (elem_type if elem_type is not None else prev_elem_type)


def _parse_element_info(lines, start, elem_type, num_elems, result):
    """Parse element connectivity table."""
    if elem_type not in ELEM_DEF:
        return start, 0

    _, nodes_per_elem, _, _ = ELEM_DEF[elem_type]
    expected_nums = 1 + nodes_per_elem + 1

    i = start + 1
    parsed = 0
    while i < len(lines) and parsed < num_elems:
        nums = _find_numbers(lines[i])
        if nums and int(nums[0]) == parsed + 1:
            while len(nums) < expected_nums and i + 1 < len(lines):
                i += 1
                nums += _find_numbers(lines[i])
        if len(nums) >= expected_nums and int(nums[0]) == parsed + 1:
            node_ids = [int(v) - 1 for v in nums[1:1 + nodes_per_elem]]
            mat_id = int(nums[1 + nodes_per_elem])
            result['elements'].append({
                'type': elem_type,
                'connectivity': node_ids,
                'material': mat_id,
            })
            parsed += 1
        i += 1
    return i, parsed


def _parse_displacements(lines, start, result):
    """Parse nodal displacements for one load case.

    Supports both NDF=3 (3-component: X,Y,Z) and NDF=6
    (6-component: X,Y,Z,RX,RY,RZ) formats.
    Rotations are stored separately in result['rotations'].
    """
    i = start + 1
    disps = []
    rotations = []
    target_count = result['num_nodes']

    while i < len(lines):
        line = lines[i]

        if _has_kw(line, 'STRESS', 'CALCULATIONS'):
            break
        if _has_kw(line, 'SOLUTION', 'TIME'):
            break
        if 'LOAD CASE' in line and 'LOAD CASE NUMBER' not in line:
            break

        nums = _find_numbers(line)
        if len(nums) >= 7:
            # NDF=6: node_id, dx, dy, dz, rx, ry, rz
            disps.append((nums[-6], nums[-5], nums[-4]))
            rotations.append((nums[-3], nums[-2], nums[-1]))
        elif len(nums) >= 4:
            disps.append((nums[-3], nums[-2], nums[-1]))

        if len(disps) >= target_count:
            break
        i += 1

    result.setdefault('rotations', [])
    result['displacements'].append(disps)
    if rotations:
        result['rotations'].append(rotations)
    return i


def _parse_stresses(lines, start, result):
    """Parse element stresses for one element group / load case."""
    i = start + 1
    stresses = {}

    # Figure out which element group this stress block belongs to.
    # The .out file prints stresses group by group, matching element order.
    # Use the per-group list (may have duplicates, e.g. two H8 groups).
    n_stress_blocks = len(result['stresses'])
    etypes = result.get('elem_type_by_group', [])
    if n_stress_blocks < len(etypes):
        elem_type = etypes[n_stress_blocks]
    elif etypes:
        elem_type = etypes[-1]
    else:
        return i

    if elem_type not in ELEM_DEF:
        return i

    _, _, stress_names, data_lines_per_elem = ELEM_DEF[elem_type]
    ncomp = len(stress_names)

    # Advance past headers to first data line
    while i < len(lines):
        line = lines[i]

        if _has_kw(line, 'SOLUTION', 'TIME'):
            result['stresses'].append(stresses)
            return i
        if 'LOAD CASE' in line and 'LOAD CASE NUMBER' not in line:
            result['stresses'].append(stresses)
            return i
        if _has_kw(line, 'STRESS', 'CALCULATIONS', 'ELEMENT', 'GROUP'):
            i -= 1  # rewind so main-loop i+=1 lands back on this header
            break

        stripped = line.strip()
        if stripped and stripped[0].isdigit():
            nums = _find_numbers(line)
            if len(nums) >= 2 and int(nums[0]) >= 1:
                break
        i += 1

    # Parse element stress values
    while i < len(lines):
        line = lines[i]

        if _has_kw(line, 'SOLUTION', 'TIME'):
            break
        if 'LOAD CASE' in line and 'LOAD CASE NUMBER' not in line:
            break
        if _has_kw(line, 'STRESS', 'CALCULATIONS', 'ELEMENT', 'GROUP'):
            i -= 1  # rewind so main-loop i+=1 lands back on this header
            break

        stripped = line.strip()
        if not stripped:
            i += 1
            continue

        nums = _find_numbers(line)
        if len(nums) >= 2:
            try:
                elem_idx = int(nums[0]) - 1

                if data_lines_per_elem == 1:
                    stresses[elem_idx] = nums[1:1 + ncomp]
                elif data_lines_per_elem == 2:
                    row1 = nums[1:]
                    i += 1
                    if i < len(lines):
                        row2 = _find_numbers(lines[i])
                        stresses[elem_idx] = (row1 + row2)[:ncomp]
                else:  # 3 lines
                    row1 = nums[1:]
                    rows = row1
                    for _ in range(data_lines_per_elem - 1):
                        i += 1
                        if i < len(lines):
                            rows += _find_numbers(lines[i])
                    stresses[elem_idx] = rows[:ncomp]
            except (ValueError, OverflowError):
                pass

        i += 1

    result['stresses'].append(stresses)
    return i


def write_vtk(result, outpath):
    """Write parsed data as legacy VTK unstructured grid."""
    nodes = result['nodes']
    elements = result['elements']
    num_nodes = len(nodes)
    num_cells = len(elements)

    if num_nodes == 0:
        raise ValueError("No nodes found in .out file")
    if num_cells == 0:
        raise ValueError("No elements found in .out file")

    with open(outpath, 'w', encoding='utf-8') as f:
        f.write('# vtk DataFile Version 2.0\n')
        f.write(f"{result['title']}\n")
        f.write('ASCII\n')
        f.write('DATASET UNSTRUCTURED_GRID\n')

        # ── Points ──
        f.write(f'\nPOINTS {num_nodes} float\n')
        for x, y, z in nodes:
            f.write(f'  {x:14.6e} {y:14.6e} {z:14.6e}\n')

        # ── Cells ──
        total_size = sum(1 + len(e['connectivity']) for e in elements)
        f.write(f'\nCELLS {num_cells} {total_size}\n')
        for e in elements:
            conn = e['connectivity']
            f.write(f'  {len(conn)} {" ".join(str(n) for n in conn)}\n')

        # ── Cell types ──
        f.write(f'\nCELL_TYPES {num_cells}\n')
        for e in elements:
            vtk_type = ELEM_DEF.get(e['type'], (0,))[0]
            f.write(f'  {vtk_type}\n')

        # ── Point data: displacements ──
        has_point_data = False
        if result['displacements']:
            disps = result['displacements'][0]
            if len(disps) == num_nodes:
                f.write(f'\nPOINT_DATA {num_nodes}\n')
                has_point_data = True
                f.write('VECTORS Displacement float\n')
                for dx, dy, dz in disps:
                    f.write(f'  {dx:14.6e} {dy:14.6e} {dz:14.6e}\n')

        # ── Point data: rotations (NDF=6) ──
        if result.get('rotations'):
            rots = result['rotations'][0]
            if len(rots) == num_nodes:
                if not has_point_data:
                    f.write(f'\nPOINT_DATA {num_nodes}\n')
                    has_point_data = True
                f.write('VECTORS Rotation float\n')
                for rx, ry, rz in rots:
                    f.write(f'  {rx:14.6e} {ry:14.6e} {rz:14.6e}\n')
        elem_types = result.get('elem_group_types', [])
        stress_data = _collect_stress_data(result, num_cells)
        has_stress = bool(stress_data[0])
        if has_stress:
            if not has_point_data:
                f.write(f'\nPOINT_DATA {num_nodes}\n')
                has_point_data = True
            _write_recovered_point_stress_data(f, nodes, elements, stress_data)
        if len(elem_types) > 1 or has_stress:
            f.write(f'\nCELL_DATA {num_cells}\n')
        if len(elem_types) > 1:
            f.write('SCALARS ElementType int 1\n')
            f.write('LOOKUP_TABLE default\n')
            for e in elements:
                f.write(f'  {e["type"]}\n')
            f.write('SCALARS Material int 1\n')
            f.write('LOOKUP_TABLE default\n')
            for e in elements:
                f.write(f'  {e["material"]}\n')

        # ── Cell data: stresses (per-group) ──
        if has_stress:
            _write_stress_data(f, stress_data, num_cells, header_already_written=(len(elem_types) > 1))


def _collect_stress_data(result, num_cells):
    """Collect mixed element-group stresses into per-cell dictionaries."""
    etypes = result.get('elem_type_by_group', [])
    counts = result.get('elem_count_by_group', [])

    # Build a map: cell_index -> stress_values_list (by group)
    # and collect all unique stress names in order
    all_names = []
    seen = set()

    # First pass: collect all stress names in order of first appearance
    for etype, stress_dict in zip(etypes, result['stresses']):
        if etype not in ELEM_DEF:
            continue
        names = ELEM_DEF[etype][2]
        for n in names:
            if n not in seen:
                all_names.append(n)
                seen.add(n)

    # Build per-cell stress map: cell_idx -> {name: value}
    cell_stress = [{} for _ in range(num_cells)]

    # Map each group's stress dict to correct cell indices.
    # Elements are stored in group order — offset tracks the start of each group
    offset = 0
    for gidx, (etype, stress_dict) in enumerate(zip(etypes, result['stresses'])):
        if etype not in ELEM_DEF:
            continue
        names = ELEM_DEF[etype][2]

        for local_idx, values in stress_dict.items():
            global_idx = offset + local_idx
            if global_idx < num_cells:
                for j, name in enumerate(names):
                    if j < len(values):
                        cell_stress[global_idx][name] = values[j]

        if gidx < len(counts):
            offset += counts[gidx]
        else:
            offset += max(stress_dict.keys()) + 1 if stress_dict else 0

    return all_names, cell_stress


def _write_stress_data(f, stress_data, num_cells, header_already_written=False):
    """Write original per-cell stress data to VTK."""
    all_names, cell_stress = stress_data
    if not all_names:
        return

    if not header_already_written:
        f.write(f'CELL_DATA {num_cells}\n')

    for name in all_names:
        f.write(f'SCALARS {name} float 1\n')
        f.write('LOOKUP_TABLE default\n')
        for ci in range(num_cells):
            v = cell_stress[ci].get(name, 0.0)
            f.write(f'  {v:14.6e}\n')

    f.write('SCALARS Von_Mises float 1\n')
    f.write('LOOKUP_TABLE default\n')
    for ci in range(num_cells):
        f.write(f'  {_von_mises(cell_stress[ci]):14.6e}\n')

    sigma_components = [
        ('Sigma-Magnitude', lambda s: _von_mises(s)),
        ('Sigma-s11', lambda s: _stress_tensor_components(s)[0]),
        ('Sigma-s22', lambda s: _stress_tensor_components(s)[1]),
        ('Sigma-s33', lambda s: _stress_tensor_components(s)[2]),
        ('Sigma-s12', lambda s: _stress_tensor_components(s)[3]),
        ('Sigma-s23', lambda s: _stress_tensor_components(s)[4]),
        ('Sigma-s13', lambda s: _stress_tensor_components(s)[5]),
    ]
    for name, getter in sigma_components:
        f.write(f'SCALARS {name} float 1\n')
        f.write('LOOKUP_TABLE default\n')
        for ci in range(num_cells):
            f.write(f'  {getter(cell_stress[ci]):14.6e}\n')

    f.write('TENSORS Sigma float\n')
    for ci in range(num_cells):
        sx, sy, sz, txy, tyz, tzx = _stress_tensor_components(cell_stress[ci])
        f.write(f'  {sx:14.6e} {txy:14.6e} {tzx:14.6e}\n')
        f.write(f'  {txy:14.6e} {sy:14.6e} {tyz:14.6e}\n')
        f.write(f'  {tzx:14.6e} {tyz:14.6e} {sz:14.6e}\n')

    f.write('TENSORS StressTensor float\n')
    for ci in range(num_cells):
        sx, sy, sz, txy, tyz, tzx = _stress_tensor_components(cell_stress[ci])
        f.write(f'  {sx:14.6e} {txy:14.6e} {tzx:14.6e}\n')
        f.write(f'  {txy:14.6e} {sy:14.6e} {tyz:14.6e}\n')
        f.write(f'  {tzx:14.6e} {tyz:14.6e} {sz:14.6e}\n')


def _write_recovered_point_stress_data(f, nodes, elements, stress_data):
    """Write superconvergent patch recovered nodal stresses."""
    all_names, cell_stress = stress_data
    if not all_names:
        return

    recovered = _recover_nodal_stress_patch(nodes, elements, all_names, cell_stress)

    for name in all_names:
        f.write(f'SCALARS SPR_{name} float 1\n')
        f.write('LOOKUP_TABLE default\n')
        for stress in recovered:
            f.write(f'  {stress.get(name, 0.0):14.6e}\n')

    f.write('SCALARS SPR_Von_Mises float 1\n')
    f.write('LOOKUP_TABLE default\n')
    for stress in recovered:
        f.write(f'  {_von_mises(stress):14.6e}\n')

    sigma_components = [
        ('SPR_Sigma-Magnitude', lambda s: _von_mises(s)),
        ('SPR_Sigma-s11', lambda s: _stress_tensor_components(s)[0]),
        ('SPR_Sigma-s22', lambda s: _stress_tensor_components(s)[1]),
        ('SPR_Sigma-s33', lambda s: _stress_tensor_components(s)[2]),
        ('SPR_Sigma-s12', lambda s: _stress_tensor_components(s)[3]),
        ('SPR_Sigma-s23', lambda s: _stress_tensor_components(s)[4]),
        ('SPR_Sigma-s13', lambda s: _stress_tensor_components(s)[5]),
    ]
    for name, getter in sigma_components:
        f.write(f'SCALARS {name} float 1\n')
        f.write('LOOKUP_TABLE default\n')
        for stress in recovered:
            f.write(f'  {getter(stress):14.6e}\n')

    f.write('TENSORS SPR_Sigma float\n')
    for stress in recovered:
        sx, sy, sz, txy, tyz, tzx = _stress_tensor_components(stress)
        f.write(f'  {sx:14.6e} {txy:14.6e} {tzx:14.6e}\n')
        f.write(f'  {txy:14.6e} {sy:14.6e} {tyz:14.6e}\n')
        f.write(f'  {tzx:14.6e} {tyz:14.6e} {sz:14.6e}\n')


def _recover_nodal_stress_patch(nodes, elements, all_names, cell_stress):
    centers = _cell_centers(nodes, elements)
    node_to_cells = [[] for _ in nodes]
    for ci, elem in enumerate(elements):
        for ni in elem['connectivity']:
            if 0 <= ni < len(node_to_cells):
                node_to_cells[ni].append(ci)

    recovered = []
    for ni, point in enumerate(nodes):
        patch = node_to_cells[ni]
        node_stress = {}
        for name in all_names:
            samples = [
                (centers[ci], cell_stress[ci][name])
                for ci in patch
                if ci < len(cell_stress) and name in cell_stress[ci]
            ]
            if samples:
                node_stress[name] = _fit_patch_value(point, samples)
        recovered.append(node_stress)
    return recovered


def _cell_centers(nodes, elements):
    centers = []
    for elem in elements:
        conn = [ni for ni in elem['connectivity'] if 0 <= ni < len(nodes)]
        if not conn:
            centers.append((0.0, 0.0, 0.0))
            continue
        inv = 1.0 / len(conn)
        centers.append((
            sum(nodes[ni][0] for ni in conn) * inv,
            sum(nodes[ni][1] for ni in conn) * inv,
            sum(nodes[ni][2] for ni in conn) * inv,
        ))
    return centers


def _fit_patch_value(point, samples):
    if len(samples) < 4:
        return _weighted_average(point, samples)

    x0 = sum(p[0] for p, _ in samples) / len(samples)
    y0 = sum(p[1] for p, _ in samples) / len(samples)
    z0 = sum(p[2] for p, _ in samples) / len(samples)
    scale = max(
        math.sqrt((p[0] - x0) ** 2 + (p[1] - y0) ** 2 + (p[2] - z0) ** 2)
        for p, _ in samples
    )
    if scale <= 1.0e-20:
        return _weighted_average(point, samples)

    normal = [[0.0] * 4 for _ in range(4)]
    rhs = [0.0] * 4
    for p, value in samples:
        basis = [
            1.0,
            (p[0] - x0) / scale,
            (p[1] - y0) / scale,
            (p[2] - z0) / scale,
        ]
        distance = math.sqrt((p[0] - point[0]) ** 2 + (p[1] - point[1]) ** 2 + (p[2] - point[2]) ** 2)
        weight = 1.0 / max(distance, scale * 1.0e-6)
        for i in range(4):
            rhs[i] += weight * basis[i] * value
            for j in range(4):
                normal[i][j] += weight * basis[i] * basis[j]

    coeff = _solve_linear_system(normal, rhs)
    if coeff is None:
        return _weighted_average(point, samples)

    basis = [
        1.0,
        (point[0] - x0) / scale,
        (point[1] - y0) / scale,
        (point[2] - z0) / scale,
    ]
    return sum(coeff[i] * basis[i] for i in range(4))


def _weighted_average(point, samples):
    weighted_sum = 0.0
    weight_total = 0.0
    for p, value in samples:
        distance = math.sqrt((p[0] - point[0]) ** 2 + (p[1] - point[1]) ** 2 + (p[2] - point[2]) ** 2)
        weight = 1.0 / max(distance, 1.0e-12)
        weighted_sum += weight * value
        weight_total += weight
    if weight_total <= 0.0:
        return samples[0][1]
    return weighted_sum / weight_total


def _solve_linear_system(matrix, rhs):
    n = len(rhs)
    a = [row[:] + [rhs[i]] for i, row in enumerate(matrix)]
    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(a[r][col]))
        if abs(a[pivot][col]) < 1.0e-18:
            return None
        if pivot != col:
            a[col], a[pivot] = a[pivot], a[col]
        inv_pivot = 1.0 / a[col][col]
        for j in range(col, n + 1):
            a[col][j] *= inv_pivot
        for r in range(n):
            if r == col:
                continue
            factor = a[r][col]
            if factor == 0.0:
                continue
            for j in range(col, n + 1):
                a[r][j] -= factor * a[col][j]
    return [a[i][n] for i in range(n)]


def _stress_tensor_components(stress):
    """Return sx, sy, sz, txy, tyz, tzx with sensible fallbacks."""
    if 'Sigma_X' in stress or 'Sigma_Y' in stress or 'Sigma_Z' in stress:
        return (
            stress.get('Sigma_X', 0.0),
            stress.get('Sigma_Y', 0.0),
            stress.get('Sigma_Z', 0.0),
            stress.get('Tau_XY', 0.0),
            stress.get('Tau_YZ', 0.0),
            stress.get('Tau_ZX', 0.0),
        )
    axial = stress.get('Axial_Stress', stress.get('Stress', 0.0))
    return axial, 0.0, 0.0, 0.0, 0.0, 0.0


def _von_mises(stress):
    sx, sy, sz, txy, tyz, tzx = _stress_tensor_components(stress)
    return (0.5 * ((sx - sy) ** 2 + (sy - sz) ** 2 + (sz - sx) ** 2)
            + 3.0 * (txy ** 2 + tyz ** 2 + tzx ** 2)) ** 0.5


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    inpath = sys.argv[1]
    if not os.path.exists(inpath):
        print(f"Error: file not found: {inpath}")
        sys.exit(1)

    outpath = sys.argv[2] if len(sys.argv) >= 3 else os.path.splitext(inpath)[0] + '.vtk'

    print(f"Converting: {inpath}")
    print(f"       -> : {outpath}")

    result = parse_out(inpath)

    from collections import Counter
    tc = Counter(e['type'] for e in result['elements'])

    print(f"  Title:        {result['title']}")
    print(f"  Nodes:        {len(result['nodes'])}")
    print(f"  Elements:     {len(result['elements'])}")
    print(f"                   " + ", ".join(
        f"type {t} x {c}" for t, c in sorted(tc.items())))
    print(f"  Elem types:   {result['elem_group_types']}")
    print(f"  Load cases:   {len(result['displacements'])}")
    print(f"  Stress sets:  {len(result['stresses'])}")

    write_vtk(result, outpath)
    print("Done.")


if __name__ == '__main__':
    main()
