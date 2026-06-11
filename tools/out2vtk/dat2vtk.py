#!/usr/bin/env python3
"""
STAP++ .dat to ParaView .vtk converter (mesh-only)

Reads a STAP++ input data file (.dat) and writes a legacy VTK unstructured
grid for mesh visualization in ParaView.  Only geometric information (nodes +
elements) is extracted — boundary conditions, loads, and material properties
are skipped because they are not rendered directly.

Usage:  python dat2vtk.py <input.dat> [output.vtk]

Supported element types:
  0  Shell4  (4-node flat shell)
  1  Truss   (2-node line)
  2  Q4      (4-node quad, full integration)
  3  Q4R     (4-node quad, reduced integration)
  4  T3      (3-node triangle, CST)
  5  H8      (8-node hexahedron)
  6  Beam    (2-node Euler-Bernoulli frame)
  7  Plate   (4-node Mindlin-Reissner plate bending)
  9  Beam3D  (2-node 3D beam, NDF=6)
  10 Shell4  (4-node flat shell, alt)
  11 Beam3DTimoshenko (2-node B31-equivalent beam)
  12 H8R     (8-node reduced-integration hexahedron)
  14 H8RPier (8-node pier-local reduced-integration hexahedron)

Both NDF=3 and NDF=6 per-node formats are auto-detected.
"""

import sys
import os


# ── VTK / element constants ─────────────────────────────────────────────
VTK_LINE       = 3
VTK_TRIANGLE   = 5
VTK_QUAD       = 9
VTK_HEXAHEDRON = 12

# elem_type -> (vtk_cell_type, nodes_per_elem, material_params_per_set)
ELEM_DEF = {
    0:  (VTK_QUAD,        4, 4),    # Shell4: mat#, E, nu, t
    1:  (VTK_LINE,        2, 3),    # Truss:  mat#, E, A
    2:  (VTK_QUAD,        4, 4),    # Q4:     mat#, E, nu, t
    3:  (VTK_QUAD,        4, 4),    # Q4R:    mat#, E, nu, t
    4:  (VTK_TRIANGLE,    3, 4),    # T3:     mat#, E, nu, t
    5:  (VTK_HEXAHEDRON,  8, 3),    # H8:     mat#, E, nu
    6:  (VTK_LINE,        2, 4),    # Beam:   mat#, E, A, I
    7:  (VTK_QUAD,        4, 4),    # Plate:  mat#, E, nu, t
    9:  (VTK_LINE,        2, 10),   # Beam3D: mat#, E, G, A, Iy, Iz, J, ox, oy, oz
    10: (VTK_QUAD,        4, 4),    # Shell4: mat#, E, nu, t  (alt)
    11: (VTK_LINE,        2, 12),   # Beam3DTimoshenko: mat#, E, nu, A, Iy, Iz, J, Asy, Asz, ox, oy, oz
    12: (VTK_HEXAHEDRON,  8, 3),    # H8R:    mat#, E, nu
    14: (VTK_HEXAHEDRON,  8, 3),    # H8RPier: mat#, E, nu
}


def _detect_ndf(tokens, pos, numnp):
    """Detect BC codes per node (3 or 6) by checking where load case data lands."""
    for ndf in (3, 6):
        check = pos + numnp * (1 + ndf + 3)
        if check < len(tokens):
            try:
                lc = int(tokens[check])
                if 1 <= lc <= 100:
                    return ndf
            except (ValueError, OverflowError):
                pass

    # If both fail, prefer 6 when file is large enough
    if pos + numnp * (1 + 6 + 3) <= len(tokens):
        return 6
    return 3


def _scan_material_boundary(tokens, pos, nmat, numnp):
    """Find where material section ends by looking for elem#1 in connectivity.

    Returns the number of material-parameter tokens per material set
    (skipping the leading mat# token, which is already counted in nmat).
    """
    best = 3
    for candidate in range(1, 20):
        end = pos + nmat * candidate
        if end >= len(tokens):
            break
        try:
            eid = int(tokens[end])
            if eid == 1:
                # Verify a few connectivity rows look valid
                all_ok = True
                for trial in range(3):
                    base = end + trial * 4  # assume 2-node elements for trial
                    if base + 3 >= len(tokens):
                        all_ok = False
                        break
                    n1 = int(tokens[base + 1])
                    n2 = int(tokens[base + 2])
                    m  = int(tokens[base + 3])
                    if not (1 <= n1 <= numnp and 1 <= n2 <= numnp and m >= 1):
                        all_ok = False
                        break
                if all_ok:
                    return candidate
        except (ValueError, OverflowError):
            pass
    return best


def _scan_nodes_per_elem(tokens, pos, nelem, numnp):
    """For unknown element types, auto-detect nodes_per_elem from connectivity.

    Scans from pos, finding the first few rows of connectivity to infer
    how many node IDs follow each element number.
    """
    # Try scanning forward to find element 1
    for offset in range(min(20, len(tokens) - pos)):
        try:
            if int(tokens[pos + offset]) == 1:
                # Found elem 1 at pos+offset
                start = pos + offset
                # Now find how many consecutive valid node IDs follow
                # before the next elem number (2) or material number
                node_ids = []
                for j in range(1, 21):  # max 20 nodes per elem
                    idx = start + j
                    if idx >= len(tokens):
                        break
                    try:
                        nid = int(tokens[idx])
                        if 1 <= nid <= numnp:
                            node_ids.append(nid)
                        else:
                            break
                    except ValueError:
                        break

                # Verify with element 2
                e2_start = start + 1 + len(node_ids) + 1  # elem# + nodes + mat#
                if e2_start + 2 < len(tokens):
                    try:
                        e2 = int(tokens[e2_start])
                        if e2 == 2:
                            ok = True
                            for j in range(1, len(node_ids) + 1):
                                if e2_start + j >= len(tokens):
                                    ok = False
                                    break
                                nid = int(tokens[e2_start + j])
                                if not (1 <= nid <= numnp):
                                    ok = False
                                    break
                            if ok:
                                mat_pos = e2_start + len(node_ids) + 1
                                if mat_pos < len(tokens):
                                    mat_id = int(tokens[mat_pos])
                                    if mat_id >= 1:
                                        return len(node_ids)
                    except (ValueError, OverflowError):
                        pass
                if node_ids:
                    return len(node_ids)
        except ValueError:
            pass
    return 2  # default fallback


def parse_dat(filepath):
    """Parse a STAP++ .dat file.  Returns (title, nodes, elements, elem_types, ndf)."""

    with open(filepath, 'r', encoding='utf-8', errors='replace') as fh:
        title = fh.readline().strip()
        token_lines = []
        for raw in fh:
            token_lines.extend(raw.split())

    if len(token_lines) < 4:
        raise ValueError("File too short: missing control line")

    tokens = [t for t in token_lines if t]   # strip any empty wrappers

    pos = 0
    numnp  = int(tokens[pos]); pos += 1
    numeg  = int(tokens[pos]); pos += 1
    nlcase = int(tokens[pos]); pos += 1
    pos += 1   # skip MODEX

    ndf = _detect_ndf(tokens, pos, numnp)

    # ── Nodes ─────────────────────────────────────────────────────────
    nodes = []
    node_end = pos + numnp * (1 + ndf + 3)

    if node_end > len(tokens):
        raise ValueError(
            f"File too short: need {node_end} tokens for {numnp} nodes "
            f"(NDF={ndf}), have {len(tokens)}"
        )

    for _ in range(numnp):
        pos += 1                        # skip node id
        pos += ndf                      # skip BC codes
        x = float(tokens[pos]); pos += 1
        y = float(tokens[pos]); pos += 1
        z = float(tokens[pos]); pos += 1
        nodes.append((x, y, z))

    # ── Load cases ────────────────────────────────────────────────────
    for _ in range(nlcase):
        pos += 1                       # skip load case number
        nload = int(tokens[pos]); pos += 1
        pos += nload * 3               # skip (node, dir, mag) triplets

    # ── Element groups ────────────────────────────────────────────────
    elements = []
    elem_types = set()

    for _g in range(numeg):
        if pos + 2 >= len(tokens):
            print(f"  [warning] Truncated at group {_g + 1}/{numeg} — stopping")
            break

        try:
            etype = int(tokens[pos])
        except (ValueError, OverflowError):
            print(f"  [warning] Non-numeric element type \"{tokens[pos]}\" "
                  f"at group {_g + 1}/{numeg} — stopping")
            break
        pos += 1
        nelem = int(tokens[pos]); pos += 1
        nmat  = int(tokens[pos]); pos += 1

        if nmat > 1000 or nelem == 0 or nelem > 100000:
            raise ValueError(
                f"Implausible group header: type={etype} nelem={nelem} nmat={nmat}"
            )

        elem_types.add(etype)

        if etype in ELEM_DEF:
            vtk_cell_type, nodes_per_elem, mat_params = ELEM_DEF[etype]
        else:
            # ── Auto-detect unknown element type ─────────────────────
            mat_params = _scan_material_boundary(tokens, pos, nmat, numnp)
            pos += nmat * mat_params
            nodes_per_elem = _scan_nodes_per_elem(tokens, pos, nelem, numnp)
            vtk_cell_type = {2: VTK_LINE, 3: VTK_TRIANGLE, 4: VTK_QUAD, 8: VTK_HEXAHEDRON}.get(
                nodes_per_elem, VTK_LINE)
            ELEM_DEF[etype] = (vtk_cell_type, nodes_per_elem, mat_params)
            print(f"  [auto-detected] type {etype}: {nodes_per_elem} nodes/elem, "
                  f"{mat_params} material params, VTK cell {vtk_cell_type}")

        # Skip material properties
        pos += nmat * mat_params

        # Read connectivity
        conn_width = 1 + nodes_per_elem + 1   # elem# + nodes + mat#

        for _ in range(nelem):
            if pos + conn_width > len(tokens):
                raise ValueError("Truncated connectivity data")
            pos += 1                    # skip element number
            node_ids = [int(tokens[pos + i]) - 1 for i in range(nodes_per_elem)]
            pos += nodes_per_elem
            pos += 1                    # skip material set number
            elements.append({
                'type':          etype,
                'connectivity':  node_ids,
            })

    return title, nodes, elements, sorted(elem_types), ndf


def write_vtk(title, nodes, elements, outpath):
    """Write parsed mesh as legacy VTK unstructured grid."""

    num_nodes = len(nodes)
    num_cells = len(elements)

    if num_nodes == 0:
        raise ValueError("No nodes found")
    if num_cells == 0:
        raise ValueError("No elements found")

    with open(outpath, 'w', encoding='utf-8') as f:
        f.write('# vtk DataFile Version 2.0\n')
        f.write(f'{title}\n')
        f.write('ASCII\n')
        f.write('DATASET UNSTRUCTURED_GRID\n')

        f.write(f'\nPOINTS {num_nodes} float\n')
        for x, y, z in nodes:
            f.write(f'  {x:14.6e} {y:14.6e} {z:14.6e}\n')

        total_size = sum(1 + len(e['connectivity']) for e in elements)
        f.write(f'\nCELLS {num_cells} {total_size}\n')
        for e in elements:
            conn = e['connectivity']
            f.write(f'  {len(conn)} {" ".join(str(n) for n in conn)}\n')

        f.write(f'\nCELL_TYPES {num_cells}\n')
        for e in elements:
            vtk_type = ELEM_DEF[e['type']][0]
            f.write(f'  {vtk_type}\n')

        f.write(f'\nCELL_DATA {num_cells}\n')
        f.write('SCALARS ElementType int 1\n')
        f.write('LOOKUP_TABLE default\n')
        for e in elements:
            f.write(f'  {e["type"]}\n')


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    inpath = sys.argv[1]
    if not os.path.exists(inpath):
        print(f"Error: file not found: {inpath}", file=sys.stderr)
        sys.exit(1)

    outpath = sys.argv[2] if len(sys.argv) >= 3 else os.path.splitext(inpath)[0] + '.vtk'

    print(f"Converting: {inpath}")
    print(f"        -> : {outpath}")

    title, nodes, elements, elem_types, ndf = parse_dat(inpath)

    from collections import Counter
    tc = Counter(e['type'] for e in elements)

    print(f"  Title:        {title}")
    print(f"  Nodes:        {len(nodes)}  (NDF={ndf})")
    print(f"  Elements:     {len(elements)}")
    print(f"                   " + ", ".join(
        f"type {t} × {c}" for t, c in sorted(tc.items())))
    print(f"  Elem types:   {elem_types}")

    write_vtk(title, nodes, elements, outpath)
    print("Done.")


if __name__ == '__main__':
    main()
