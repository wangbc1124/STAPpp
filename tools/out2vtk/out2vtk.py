#!/usr/bin/env python3
"""
STAP++ .out to ParaView .vtk converter

Converts finite element output from STAP++ (.out) to legacy VTK unstructured
grid format for visualization in ParaView.

Usage: python out2vtk.py <input.out> [output.vtk]

Supports all element types:
  1 - Truss (2-node line)
  2 - Q4   (4-node quad, full integration)
  3 - Q4R  (4-node quad, reduced integration)
  4 - T3   (3-node triangle, CST)
  5 - H8   (8-node hexahedron)
  6 - Beam (2-node Euler-Bernoulli frame)
  7 - Plate (4-node Mindlin-Reissner plate bending)
"""

import re
import sys
import os

# ── VTK cell type constants ───────────────────────────────────────────
VTK_LINE       = 3
VTK_TRIANGLE   = 5
VTK_QUAD       = 9
VTK_HEXAHEDRON = 12

# ── Element type definitions ──────────────────────────────────────────
# elem_type -> (vtk_cell_type, nodes_per_elem, stress_names, stress_data_lines)
ELEM_DEF = {
    1: (VTK_LINE,        2, ['Force',  'Stress'],  1),
    2: (VTK_QUAD,        4, ['Sigma_X','Sigma_Y','Tau_XY'],  1),
    3: (VTK_QUAD,        4, ['Sigma_X','Sigma_Y','Tau_XY'],  1),
    4: (VTK_TRIANGLE,    3, ['Sigma_X','Sigma_Y','Tau_XY'],  1),
    5: (VTK_HEXAHEDRON,  8, ['Sigma_X','Sigma_Y','Sigma_Z',
                              'Tau_XY', 'Tau_YZ', 'Tau_ZX'],  2),
    6: (VTK_LINE,        2, ['Axial_Stress','Bend_Stress1','Bend_Stress2'], 1),
    7: (VTK_QUAD,        4, ['Mx','My','Mxy','Qx','Qy'], 2),
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
    """Parse a single element group: type, count, materials, connectivity."""
    i = start + 1
    elem_type = None
    num_elems = 0

    while i < len(lines):
        line = lines[i]

        if 'TOTAL SYSTEM DATA' in line:
            break
        if _has_kw(line, 'ELEMENT', 'GROUP', 'DATA') and i > start + 1:
            break

        if 'NPAR(1)' in line:
            m = re.search(r'=\s*(\d+)', line)
            if m:
                elem_type = int(m.group(1))
                if elem_type not in result['elem_group_types']:
                    result['elem_group_types'].append(elem_type)

        if 'NPAR(2)' in line:
            m = re.search(r'=\s*(\d+)', line)
            if m:
                num_elems = int(m.group(1))

        if _has_kw(line, 'ELEMENT', 'INFORMATION') and elem_type is not None:
            i = _parse_element_info(lines, i, elem_type, num_elems, result)

        i += 1

    return i, (elem_type if elem_type is not None else prev_elem_type)


def _parse_element_info(lines, start, elem_type, num_elems, result):
    """Parse element connectivity table."""
    if elem_type not in ELEM_DEF:
        return start

    _, nodes_per_elem, _, _ = ELEM_DEF[elem_type]
    expected_nums = 1 + nodes_per_elem + 1  # elem# + nodes + material#

    i = start + 1
    parsed = 0
    while i < len(lines) and parsed < num_elems:
        nums = _find_numbers(lines[i])
        if len(nums) >= expected_nums:
            elem_num = int(nums[0])
            node_ids = [int(v) - 1 for v in nums[1:1 + nodes_per_elem]]
            mat_id = int(nums[1 + nodes_per_elem])
            result['elements'].append({
                'type': elem_type,
                'connectivity': node_ids,
                'material': mat_id,
            })
            parsed += 1
        i += 1
    return i


def _parse_displacements(lines, start, result):
    """Parse nodal displacements for one load case."""
    i = start + 1
    disps = []
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
        if len(nums) >= 4:
            disps.append((nums[-3], nums[-2], nums[-1]))

        if len(disps) >= target_count:
            break
        i += 1

    result['displacements'].append(disps)
    return i


def _parse_stresses(lines, start, result):
    """Parse element stresses for one element group / load case."""
    i = start + 1
    stresses = {}

    elem_types = result.get('elem_group_types', [])
    if not elem_types:
        return i
    elem_type = elem_types[-1]  # current group

    if elem_type not in ELEM_DEF:
        return i

    _, _, stress_names, data_lines_per_elem = ELEM_DEF[elem_type]
    ncomp = len(stress_names)

    # Advance past the stress-component header to the first data line.
    # Data lines begin with an integer (the element number).
    while i < len(lines):
        line = lines[i]

        if _has_kw(line, 'SOLUTION', 'TIME'):
            result['stresses'].append(stresses)
            return i
        if 'LOAD CASE' in line and 'LOAD CASE NUMBER' not in line:
            result['stresses'].append(stresses)
            return i
        if _has_kw(line, 'STRESS', 'CALCULATIONS', 'ELEMENT', 'GROUP'):
            break

        stripped = line.strip()
        if stripped and stripped[0].isdigit():
            nums = _find_numbers(line)
            if len(nums) >= 2:
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
                else:
                    row1 = nums[1:]
                    i += 1
                    if i < len(lines):
                        row2 = _find_numbers(lines[i])
                        combined = row1 + row2
                        stresses[elem_idx] = combined[:ncomp]
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
        f.write(f'POINTS {num_nodes} float\n')
        for x, y, z in nodes:
            f.write(f'{x:14.6e} {y:14.6e} {z:14.6e}\n')

        # ── Cells (connectivity) ──
        total_size = sum(1 + len(e['connectivity']) for e in elements)
        f.write(f'\nCELLS {num_cells} {total_size}\n')
        for e in elements:
            conn = e['connectivity']
            f.write(str(len(conn)) + ' ' + ' '.join(str(n) for n in conn) + '\n')

        # ── Cell types ──
        f.write(f'\nCELL_TYPES {num_cells}\n')
        for e in elements:
            vtk_type = ELEM_DEF.get(e['type'], (0,))[0]
            f.write(f'{vtk_type}\n')

        # ── Point data: displacements ──
        if result['displacements']:
            disps = result['displacements'][0]
            if len(disps) == num_nodes:
                f.write(f'\nPOINT_DATA {num_nodes}\n')
                f.write('VECTORS Displacement float\n')
                for dx, dy, dz in disps:
                    f.write(f'{dx:14.6e} {dy:14.6e} {dz:14.6e}\n')

        # ── Cell data: stresses ──
        if result['stresses']:
            all_stress = {}
            for stress_dict in result['stresses']:
                all_stress.update(stress_dict)

            if all_stress and result.get('elem_group_types'):
                elem_type = result['elem_group_types'][0]
                if elem_type in ELEM_DEF:
                    stress_names = ELEM_DEF[elem_type][2]
                    f.write(f'\nCELL_DATA {num_cells}\n')
                    for comp_idx, name in enumerate(stress_names):
                        f.write(f'SCALARS {name} float 1\n')
                        f.write('LOOKUP_TABLE default\n')
                        for cell_idx in range(num_cells):
                            vals = all_stress.get(cell_idx, [0.0] * len(stress_names))
                            v = vals[comp_idx] if comp_idx < len(vals) else 0.0
                            f.write(f'{v:14.6e}\n')


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

    print(f"  Title:        {result['title']}")
    print(f"  Nodes:        {len(result['nodes'])}")
    print(f"  Elements:     {len(result['elements'])}")
    print(f"  Elem types:   {result['elem_group_types']}")
    print(f"  Load cases:   {len(result['displacements'])}")
    print(f"  Stress sets:  {len(result['stresses'])}")

    write_vtk(result, outpath)
    print(f"Done.")


if __name__ == '__main__':
    main()
