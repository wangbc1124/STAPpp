#!/usr/bin/env python3
"""
Convert Abaqus .inp file to STAP++ .dat format — FULL 3D VERSION.

Preserves 3D coordinates, element types, and instance orientations.
Element type mappings:
  - T3D2 (2-node 3D truss)  -> Bar  (type 1) — fully 3D
  - S4R  (4-node shell)      -> Q4   (type 2) — plane stress in XY plane
  - C3D8R (8-node brick)     -> H8   (type 5) — fully 3D
  - B31  (2-node 3D beam)    -> Beam (type 6) — 2D frame in XY plane

Gravity direction: Abaqus (0,0,-1) -> STAP++ DOF 3 (Z).
"""

import math
import numpy as np
from collections import defaultdict

# ============================================================
# Rodrigues rotation formula
# ============================================================

def rotation_matrix(axis_point_a, axis_point_b, angle_deg):
    """Compute 3x3 rotation matrix from two points on axis and angle in degrees.

    Abaqus format: xa, ya, za, xb, yb, zb, angle
    Axis direction = (xb-xa, yb-ya, zb-za), angle in degrees.
    """
    ax, ay, az = axis_point_a
    bx, by, bz = axis_point_b
    ux = bx - ax
    uy = by - ay
    uz = bz - az
    norm = math.sqrt(ux*ux + uy*uy + uz*uz)
    if norm < 1e-15:
        return np.eye(3)
    ux /= norm
    uy /= norm
    uz /= norm
    theta = math.radians(angle_deg)
    c = math.cos(theta)
    s = math.sin(theta)
    R = np.array([
        [c + ux*ux*(1-c),      ux*uy*(1-c) - uz*s,   ux*uz*(1-c) + uy*s],
        [uy*ux*(1-c) + uz*s,   c + uy*uy*(1-c),       uy*uz*(1-c) - ux*s],
        [uz*ux*(1-c) - uy*s,   uz*uy*(1-c) + ux*s,    c + uz*uz*(1-c)]
    ])
    return R


# ============================================================
# 1. Parse Abaqus .inp file
# ============================================================

def parse_abaqus_inp(filepath):
    """Parse an Abaqus .inp file."""
    with open(filepath, 'r') as f:
        lines = f.readlines()

    clean = []
    for line in lines:
        s = line.strip()
        if s.startswith('**'):
            continue
        clean.append(s)

    data = {
        'parts': {},
        'assembly': {'instances': [], 'nsets': [], 'elsets': [], 'ties': [], 'surfaces': {}},
        'materials': {},
        'step': {'boundary': [], 'loads': [], 'static': None}
    }

    i = 0
    while i < len(clean):
        line = clean[i]

        if line.startswith('*Part,'):
            name = line.split('name=')[1].split(',')[0].strip()
            part = parse_part(clean, i + 1)
            data['parts'][name] = part
            i = part['_end']

        elif line.startswith('*Assembly,'):
            asm = parse_assembly(clean, i + 1)
            data['assembly'] = asm
            i = asm['_end']

        elif line.startswith('*Material,'):
            name = line.split('name=')[1].strip()
            mat = {}
            i += 1
            while i < len(clean):
                cl = clean[i]
                if cl.startswith('**') or cl.startswith('*Material,') or \
                   cl.startswith('*Part,') or cl.startswith('*Assembly,') or \
                   cl.startswith('*Step,') or cl.startswith('*End Assembly'):
                    break
                if cl.startswith('*Density'):
                    i += 1
                    if i < len(clean):
                        mat['density'] = float(clean[i].split(',')[0])
                elif cl.startswith('*Elastic'):
                    i += 1
                    if i < len(clean):
                        parts = clean[i].split(',')
                        mat['E'] = float(parts[0])
                        mat['nu'] = float(parts[1])
                i += 1
            data['materials'][name] = mat
            continue

        elif line.startswith('*Step,'):
            step = parse_step(clean, i + 1)
            data['step'] = step
            i = step['_end']

        elif line.startswith('*End Assembly'):
            i += 1
        else:
            i += 1

    return data


def parse_part(lines, start):
    """Parse a *Part section."""
    nodes = {}
    elements = []
    sections = {}
    nsets = {}
    elsets = {}

    i = start
    while i < len(lines):
        line = lines[i]
        if line.startswith('*End Part'):
            break

        if line.startswith('*Node'):
            i += 1
            while i < len(lines):
                s = lines[i]
                if s.startswith('*'):
                    break
                parts = s.split(',')
                if len(parts) >= 4:
                    nid = int(parts[0])
                    nodes[nid] = (float(parts[1]), float(parts[2]), float(parts[3]))
                i += 1
            continue

        elif line.startswith('*Element'):
            etype = 'UNKNOWN'
            if 'type=' in line:
                etype = line.split('type=')[1].split(',')[0].strip()
            i += 1
            while i < len(lines):
                s = lines[i]
                if s.startswith('*'):
                    break
                parts = s.split(',')
                if len(parts) >= 2:
                    eid = int(parts[0])
                    conn = [int(x.strip()) for x in parts[1:]]
                    elements.append({'id': eid, 'type': etype, 'nodes': conn})
                i += 1
            continue

        elif line.startswith('*Nset'):
            nset_name = line.split('nset=')[1].split(',')[0].strip()
            is_gen = 'generate' in line
            ns = parse_set_lines(lines, i + 1, is_gen)
            nsets[nset_name] = ns
            i = ns['_end']
            continue

        elif line.startswith('*Elset'):
            eset_name = line.split('elset=')[1].split(',')[0].strip()
            is_gen = 'generate' in line
            es = parse_set_lines(lines, i + 1, is_gen)
            elsets[eset_name] = es
            i = es['_end']
            continue

        elif line.startswith('*Solid Section') or line.startswith('*Shell Section') or \
             line.startswith('*Beam Section'):
            sec = {'type': line.split('*')[1].split(',')[0].strip()}
            if 'material=' in line:
                sec['material'] = line.split('material=')[1].split(',')[0].strip()
            if 'elset=' in line:
                sec['elset'] = line.split('elset=')[1].split(',')[0].strip()
            i += 1
            if i < len(lines) and not lines[i].startswith('*'):
                sec['data'] = [float(x.strip()) for x in lines[i].split(',') if x.strip()]
                i += 1
                # For Beam Section, also read n1 direction vector on next line
                if 'Beam Section' in sec['type']:
                    if i < len(lines) and not lines[i].startswith('*'):
                        try:
                            n1_vals = [float(x.strip()) for x in lines[i].split(',') if x.strip()]
                            if len(n1_vals) >= 3:
                                sec['n1'] = n1_vals[:3]
                            i += 1
                        except ValueError:
                            pass  # not a float line, skip
            sections[sec.get('elset', '')] = sec
            continue

        else:
            i += 1

    return {
        'nodes': nodes,
        'elements': elements,
        'sections': sections,
        'nsets': nsets,
        'elsets': elsets,
        '_end': i
    }


def parse_set_lines(lines, start, is_generate=False):
    """Parse node/element set data (comma-separated numbers, possibly multi-line).

    If is_generate is True and exactly 3 numbers found, expands as range(start, end+1, step).
    """
    numbers = []
    i = start
    while i < len(lines):
        s = lines[i].strip()
        if s.startswith('*'):
            break
        if s:
            for tok in s.split(','):
                tok = tok.strip()
                if tok:
                    try:
                        numbers.append(int(tok))
                    except ValueError:
                        pass
        i += 1

    if is_generate and len(numbers) == 3:
        start_val, end_val, step = numbers
        numbers = list(range(start_val, end_val + 1, step))

    return {'numbers': numbers, '_end': i}


def parse_assembly(lines, start):
    """Parse *Assembly section including instances, nsets, ties."""
    instances = []
    nsets = []
    elsets = []
    ties = []
    surfaces = {}

    i = start
    while i < len(lines):
        line = lines[i]
        if line.startswith('*End Assembly'):
            break

        if line.startswith('*Instance,'):
            name = line.split('name=')[1].split(',')[0].strip()
            part = line.split('part=')[1].split(',')[0].strip()
            translation = (0.0, 0.0, 0.0)
            rotation_data = None  # (axis_a, axis_b, angle_deg)

            i += 1
            # Read translation
            if i < len(lines) and not lines[i].startswith('*'):
                parts = [float(x.strip()) for x in lines[i].split(',') if x.strip()]
                if len(parts) == 3:
                    translation = tuple(parts)
                    i += 1
            # Read rotation
            if i < len(lines) and not lines[i].startswith('*'):
                parts = [float(x.strip()) for x in lines[i].split(',') if x.strip()]
                if len(parts) == 7:
                    rotation_data = (tuple(parts[:3]), tuple(parts[3:6]), parts[6])
                    i += 1

            instances.append({
                'name': name, 'part_name': part,
                'translation': translation,
                'rotation': rotation_data
            })

        elif line.startswith('*Nset,'):
            nsets.append(parse_assembly_set(lines, i))
            i = nsets[-1]['_end']
            continue

        elif line.startswith('*Elset,'):
            elsets.append(parse_assembly_set(lines, i))
            i = elsets[-1]['_end']
            continue

        elif line.startswith('*Tie,'):
            tie = {'name': line.split('name=')[1].split(',')[0].strip()}
            i += 1
            if i < len(lines) and not lines[i].startswith('*'):
                parts = lines[i].split(',')
                if len(parts) >= 2:
                    tie['slave'] = parts[0].strip()
                    tie['master'] = parts[1].strip()
            ties.append(tie)
            i += 1
            continue

        elif line.startswith('*Surface,'):
            # *Surface, type=NODE, name=m_Set-6_CNS_, internal
            surf_name = line.split('name=')[1].split(',')[0].strip()
            i += 1
            if i < len(lines) and not lines[i].startswith('*'):
                surf_parts = [x.strip() for x in lines[i].split(',') if x.strip()]
                # First part is the nset name, second is tolerance
                if surf_parts:
                    surfaces[surf_name] = surf_parts[0]  # nset name
            i += 1
            continue

        else:
            i += 1

    if not lines[i].startswith('*End Assembly'):
        i += 1  # skip past *End Assembly

    return {
        'instances': instances, 'nsets': nsets, 'elsets': elsets, 'ties': ties,
        'surfaces': surfaces, '_end': i
    }


def parse_assembly_set(lines, start):
    """Parse assembly-level nset or elset.

    Formats:
      *Nset, nset=NAME, instance=INSTANCE
      *Nset, nset=NAME, instance=INSTANCE, generate
      *Elset, elset=NAME, instance=INSTANCE
      *Elset, elset=NAME, instance=INSTANCE, generate
    """
    header = lines[start]
    name = ''
    instance = ''
    is_generate = False

    for part in header.split(','):
        p = part.strip()
        if p.startswith('nset='):
            name = p.split('=', 1)[1]
        elif p.startswith('elset='):
            name = p.split('=', 1)[1]
        elif p.startswith('instance='):
            instance = p.split('=', 1)[1]
        elif p.strip() == 'generate':
            is_generate = True

    numbers = parse_set_lines(lines, start + 1, is_generate)
    return {
        'name': name, 'instance': instance,
        'numbers': numbers['numbers'],
        '_end': numbers['_end']
    }


def parse_step(lines, start):
    """Parse *Step section for BCs, loads."""
    boundary = []
    loads = []
    static = None

    i = start
    while i < len(lines):
        line = lines[i]
        if line.startswith('*End Step'):
            break

        if line.startswith('*Static'):
            i += 1
            if i < len(lines):
                static = lines[i].strip()
            i += 1
            continue

        elif line.startswith('*Boundary'):
            i += 1
            while i < len(lines):
                s = lines[i]
                if s.startswith('*'):
                    break
                parts = s.split(',')
                if len(parts) >= 3:
                    set_name = parts[0].strip()
                    dof_first = int(parts[1].strip())
                    dof_last = int(parts[2].strip())
                    boundary.append({
                        'set': set_name, 'dof_first': dof_first, 'dof_last': dof_last
                    })
                i += 1
            continue

        elif line.startswith('*Dload'):
            i += 1
            while i < len(lines):
                s = lines[i]
                if s.startswith('*'):
                    break
                if 'GRAV' in s:
                    parts = [x.strip() for x in s.split(',')]
                    grav_idx = parts.index('GRAV')
                    mag = float(parts[grav_idx + 1])
                    gx = float(parts[grav_idx + 2])
                    gy = float(parts[grav_idx + 3])
                    gz = float(parts[grav_idx + 4])
                    loads.append({
                        'type': 'gravity', 'magnitude': mag,
                        'direction': (gx, gy, gz)
                    })
                i += 1
            continue

        else:
            i += 1

    return {'boundary': boundary, 'loads': loads, 'static': static, '_end': i}


# ============================================================
# 2. Flatten assembly: apply 3D transforms, merge coincident nodes
# ============================================================

def flatten_assembly(data, merge_tol=1e-3):
    """Apply 3D instance transformations and build global node/element tables.

    Nodes from each instance are kept separate (no XYZ merging).
    Only tie constraints merge slave nodes into master nodes.
    This preserves the Abaqus node count.
    """

    part_data = data['parts']
    assembly = data['assembly']
    materials = data['materials']

    global_nodes = {}    # gid -> (x, y, z)
    node_map = {}         # (instance_name, local_nid) -> global_nid
    global_elements = []  # [{id, type, nodes, material, section_data}]
    element_groups = defaultdict(list)

    next_gid = 1
    next_eid = 1

    def create_node(x, y, z):
        nonlocal next_gid
        gid = next_gid
        next_gid += 1
        global_nodes[gid] = (x, y, z)
        return gid

    for inst in assembly['instances']:
        inst_name = inst['name']
        part_name = inst['part_name']
        part = part_data.get(part_name)
        if part is None:
            print(f"  WARNING: Part '{part_name}' not found for instance '{inst_name}'")
            continue

        tx, ty, tz = inst['translation']
        rot_data = inst['rotation']

        # Transform local nodes to global (with offset axis handling)
        local_to_global = {}
        for lnid, (lx, ly, lz) in part['nodes'].items():
            v = np.array([lx, ly, lz])

            if rot_data:
                R = rotation_matrix(rot_data[0], rot_data[1], rot_data[2])
                # Abaqus rotates about the GLOBAL ORIGIN, not about axis_a.
                # The two axis points only define the direction vector.
                # Rotation order: 1) rotate about origin, 2) translate.
                v_rot = R @ v
            else:
                v_rot = v

            gx = v_rot[0] + tx
            gy = v_rot[1] + ty
            gz = v_rot[2] + tz
            gid = create_node(gx, gy, gz)
            local_to_global[lnid] = gid
            node_map[(inst_name, lnid)] = gid

        # Map elements
        section_data = None
        for elem in part['elements']:
            etype = elem['type']
            # Find section/material for this element
            mat_name = None
            sec_info = None
            for eset_name, eset in part['elsets'].items():
                if elem['id'] in eset['numbers']:
                    for sec_eset, sec in part['sections'].items():
                        if sec_eset.split('.')[-1] == eset_name or sec_eset == eset_name:
                            mat_name = sec.get('material', '')
                            sec_info = {'data': sec.get('data', None),
                                        'n1': sec.get('n1', None)}
                            break
                    break

            global_conn = [local_to_global[lnid] for lnid in elem['nodes']]
            global_elements.append({
                'id': next_eid,
                'type': etype,
                'nodes': global_conn,
                'material': mat_name or '',
                'section': sec_info,
                'instance': inst_name
            })
            element_groups[(etype, mat_name or '')].append(next_eid)
            next_eid += 1

    # === Process tie constraints: merge slave nodes into master nodes ===
    slave_to_master = {}
    surfaces = assembly.get('surfaces', {})

    for tie in assembly.get('ties', []):
        slave_surf = tie.get('slave', '')
        master_surf = tie.get('master', '')

        # Map surface name to nset name (surface naming convention: {nset}_CNS_)
        slave_nset = surfaces.get(slave_surf, slave_surf.replace('_CNS_', ''))
        master_nset = surfaces.get(master_surf, master_surf.replace('_CNS_', ''))

        # Find the assembly nsets for slave and master
        slave_nodes = []
        master_nodes = []

        for nset in assembly.get('nsets', []):
            if nset['name'] == slave_nset:
                inst_name = nset.get('instance', '')
                for lnid in nset['numbers']:
                    gid = node_map.get((inst_name, lnid))
                    if gid:
                        slave_nodes.append(gid)
            if nset['name'] == master_nset:
                inst_name = nset.get('instance', '')
                for lnid in nset['numbers']:
                    gid = node_map.get((inst_name, lnid))
                    if gid:
                        master_nodes.append(gid)

        if not slave_nodes or not master_nodes:
            print(f"  WARNING: Tie '{tie.get('name','')}' has empty node set")
        elif len(slave_nodes) != len(master_nodes):
            # Mismatched sizes: merge all slaves to first master
            target = master_nodes[0]
            for s_gid in slave_nodes:
                slave_to_master[s_gid] = target
        else:
            for s_gid, m_gid in zip(slave_nodes, master_nodes):
                slave_to_master[s_gid] = m_gid

    if slave_to_master:
        print(f"  Tie constraints: merging {len(slave_to_master)} slave nodes into masters")

        # Compute transitive closure: if A->B and B->C, then A->C
        # This handles chained ties (e.g., Cable->Pier, Pier->RiverBank)
        def resolve_transitive(n):
            """Follow the chain of slave->master mappings transitively."""
            visited = set()
            original = n
            while n in slave_to_master:
                if n in visited:
                    # Cycle detected — break cycle by keeping last valid mapping
                    print(f"  WARNING: Cycle detected in tie constraints for node {original}")
                    break
                visited.add(n)
                n = slave_to_master[n]
            return n

        # Build transitive mapping
        for s_gid in list(slave_to_master.keys()):
            slave_to_master[s_gid] = resolve_transitive(s_gid)

        # Remap element connectivities using transitive mapping
        for elem in global_elements:
            elem['nodes'] = [slave_to_master.get(n, n) for n in elem['nodes']]

    # Renumber nodes sequentially (1..N) — STAP++ requires this
    old_to_new = {}
    new_nodes = {}
    new_id = 1
    for old_nid in sorted(global_nodes.keys()):
        old_to_new[old_nid] = new_id
        new_nodes[new_id] = global_nodes[old_nid]
        new_id += 1
    global_nodes = new_nodes

    # Update element connectivities
    for elem in global_elements:
        elem['nodes'] = [old_to_new[n] for n in elem['nodes']]

    # Update node_map
    new_node_map = {}
    for key, gid in node_map.items():
        if gid in old_to_new:
            new_node_map[key] = old_to_new[gid]
    node_map = new_node_map

    # Store element_groups info for later use
    return global_nodes, global_elements, element_groups, node_map, materials


# ============================================================
# 3. Resolve boundary conditions
# ============================================================

def resolve_fixed_nodes(data, node_map):
    """Find global node IDs that should have fixed BCs based on assembly nsets."""
    assembly = data['assembly']
    step = data['step']

    # Build (instance_name, local_nid) -> {BC DOFs}
    bc_dofs = {}  # key -> set of DOF numbers (1,2,3)

    for bc in step['boundary']:
        set_name = bc['set']
        dofs = set(range(bc['dof_first'], bc['dof_last'] + 1))

        # Find this set in assembly nsets
        for nset in assembly['nsets']:
            if nset['name'] == set_name:
                inst_name = nset['instance']
                for lnid in nset['numbers']:
                    key = (inst_name, lnid)
                    if key not in bc_dofs:
                        bc_dofs[key] = set()
                    bc_dofs[key] |= dofs

    # Convert to global node -> (bc_x, bc_y, bc_z)
    fixed = {}
    for (inst_name, lnid), dofs in bc_dofs.items():
        gid = node_map.get((inst_name, lnid))
        if gid:
            if gid not in fixed:
                fixed[gid] = [0, 0, 0, 1, 1, 1]
            for d in dofs:
                fixed[gid][d - 1] = 1  # Abaqus DOF 1,2,3 -> bc[0],[1],[2]

    return fixed


# ============================================================
# 4. Write STAP++ .dat file
# ============================================================


def box_section_props(a, b, t1, t2, t3, t4):
    """Compute cross-section properties for a rectangular box section."""
    tw = t1 + t2
    tf = t3 + t4
    area = a * b - (a - tw) * (b - tf)
    Iy = (b * a**3 - (b - tf) * (a - tw)**3) / 12.0
    Iz = (a * b**3 - (a - tw) * (b - tf)**3) / 12.0
    a0 = a - tw/2.0
    b0 = b - tf/2.0
    A0 = a0 * b0
    sum_st = a0/t1 + a0/t2 + b0/t3 + b0/t4
    J = 4.0 * A0**2 / sum_st if sum_st > 1e-15 else 0.0
    return area, Iy, Iz, J

STAP_TYPE_MAP = {
    'T3D2': 1,     # Bar (3D truss, steel cables)
    'S4R': 10,     # Shell4 (flat shell: membrane + plate bending)
    'C3D8R': 5,    # H8 (3D hexahedral)
    'B31': 9,      # Beam3D (3D Euler-Bernoulli beam)
}

# 3D element type numbers (have stiffness in all 3 directions)
THREE_D_TYPES = {1, 5, 9, 10}

def polygon_area(coords):
    """Shoelace formula for polygon area."""
    n = len(coords)
    area = 0.0
    for i in range(n):
        x1, y1 = coords[i]
        x2, y2 = coords[(i + 1) % n]
        area += x1 * y2 - x2 * y1
    return abs(area) / 2.0


def write_stap_dat(output_path, data, global_nodes, global_elements,
                   element_groups, node_map):
    """Write STAP++ .dat file in 3D format."""

    materials = data['materials']
    step = data['step']

    # Resolve BCs
    bc_map = resolve_fixed_nodes(data, node_map)
    print(f"  Fixed nodes (from BCs): {len(bc_map)}")

    # Filter to valid element types
    valid_elems = []
    skipped = defaultdict(int)
    for elem in global_elements:
        st = STAP_TYPE_MAP.get(elem['type'])
        if st is None:
            skipped[elem['type']] += 1
            continue
        valid_elems.append(elem)

    for etype, count in skipped.items():
        print(f"  WARNING: Unknown type {etype}, skipping {count} elements")

    # Organize by STAP++ type + abaqus type + material
    # (separate abaqus type is needed for Bar to differentiate cable vs support beam)
    stap_groups = defaultdict(list)
    for elem in valid_elems:
        st = STAP_TYPE_MAP[elem['type']]
        mat = elem['material']
        abaqus_type = elem['type']
        stap_groups[(st, abaqus_type, mat)].append(elem['id'])

    # Identify connected nodes (nodes referenced by at least one valid element)
    connected = set()
    node_elem_types = defaultdict(set)  # node -> set of STAP++ type numbers
    for elem in valid_elems:
        st = STAP_TYPE_MAP[elem['type']]
        for n in elem['nodes']:
            connected.add(n)
            node_elem_types[n].add(st)

    orphan_nodes = set(global_nodes.keys()) - connected
    if orphan_nodes:
        print(f"  Orphan nodes (will be fully fixed): {len(orphan_nodes)}")

    # Nodes only connected to 2D elements (Q4) need Z fixed
    z_fix_nodes = set()
    for nid in connected:
        types = node_elem_types.get(nid, set())
        if types and not (types & THREE_D_TYPES):
            z_fix_nodes.add(nid)
    print(f"  Z-fixed nodes (2D elements only): {len(z_fix_nodes)}")

    # Nodes connected to Shell4 (type 10) at mesh boundaries have insufficient
    # rotational stiffness from drilling DOF stabilization alone. Fix rx,ry,rz
    # for corner/edge nodes (≤ 2 Shell4 elements) to prevent spurious large rotations.
    shell4_boundary_fix = set()  # node -> set of DOF indices {3,4,5} to fix
    node_shell4_count = defaultdict(int)
    for elem in valid_elems:
        st = STAP_TYPE_MAP[elem['type']]
        if st == 10:  # Shell4
            for n in elem['nodes']:
                node_shell4_count[n] += 1

    for nid, count in node_shell4_count.items():
        if count <= 2 and 10 in node_elem_types.get(nid, set()):
            shell4_boundary_fix.add(nid)
    print(f"  Shell4 boundary nodes (rotation-fixed): {len(shell4_boundary_fix)}")

    # Nodes connected only to Bar elements (type 1) with coplanar bars need
    # out-of-plane DOFs fixed. SupportBeam (B31) bars are in XZ plane -> fix Y.
    # Cable (T3D2) bars are 3D -> check direction vectors.
    bar_fix_nodes = defaultdict(set)  # nid -> set of DOF indices to fix

    # Build node -> [connected bar elements] mapping
    node_bar_elems = defaultdict(list)
    for elem in valid_elems:
        st = STAP_TYPE_MAP[elem['type']]
        if st == 1:  # Bar
            for n in elem['nodes']:
                node_bar_elems[n].append(elem)

    for nid, bar_elems in node_bar_elems.items():
        # Only check nodes that are connected EXCLUSIVELY to Bar elements
        types = node_elem_types.get(nid, set())
        if types != {1}:  # Has other types (Q4, H8) -> skip
            continue

        # Collect unit direction vectors of all connected Bar elements
        directions = []
        for elem in bar_elems:
            if len(elem['nodes']) != 2:
                continue
            n1, n2 = elem['nodes'][0], elem['nodes'][1]
            other = n2 if n1 == nid else n1
            p_self = np.array(global_nodes[nid])
            p_other = np.array(global_nodes[other])
            vec = p_other - p_self
            norm = np.linalg.norm(vec)
            if norm > 1e-12:
                directions.append(vec / norm)

        if not directions:
            continue

        # Stack direction vectors into a matrix (N x 3)
        dirs = np.array(directions)  # shape: (N, 3)
        N, D = dirs.shape  # N = number of bars, D = 3

        # Compute SVD to find the nullspace of the direction matrix
        U, s, Vt = np.linalg.svd(dirs, full_matrices=False)

        # Use RELATIVE threshold for singular values: if a singular value is
        # less than 1e-3 of the maximum, the corresponding direction has
        # negligible stiffness. Absolute threshold 1e-6 was too tight —
        # numerically collinear bars can have sv2 ~ 1.5e-6.
        max_sv = s[0] if len(s) > 0 else 1.0
        sv_threshold = max_sv * 1.0e-3

        # Identify nullspace directions from near-zero singular values
        null_dirs = []
        for i, sv in enumerate(s):
            if sv < sv_threshold:
                null_dirs.append(Vt[i])

        # For N < 3, SVD gives fewer than 3 singular values.
        # The remaining nullspace directions are perpendicular to ALL bar directions.
        if N < 3:
            if N == 1:
                # One bar: two perpendicular directions needed
                bar_dir = dirs[0]
                v1 = np.array([-bar_dir[1], bar_dir[0], 0.0])
                if np.linalg.norm(v1) < 1e-10:
                    v1 = np.array([-bar_dir[2], 0.0, bar_dir[0]])
                v1 /= np.linalg.norm(v1)
                v2 = np.cross(bar_dir, v1)
                v2 /= np.linalg.norm(v2)
                if len(null_dirs) == 0:
                    null_dirs = [v1, v2]
                elif len(null_dirs) == 1:
                    # One from SVD, need the second perpendicular
                    null_dirs.append(v1 if abs(np.dot(v1, null_dirs[0])) < 0.5 else v2)
            elif N == 2:
                # Two bars: check if they are effectively collinear
                # Use both cross product norm AND SVD condition number
                n = np.cross(dirs[0], dirs[1])
                n_norm = np.linalg.norm(n)
                # Collinear if: cross product near zero, OR sv ratio indicates
                # effectively 1D span (sv[1] < sv_threshold means second
                # direction has < 0.1% of primary stiffness)
                is_collinear = (n_norm < 1.0e-6) or (len(s) >= 2 and s[1] < sv_threshold)

                if is_collinear:
                    # Effectively collinear => 1D span => need 2 nullspace dirs
                    bar_dir = dirs[0]
                    v1 = np.array([-bar_dir[1], bar_dir[0], 0.0])
                    if np.linalg.norm(v1) < 1e-10:
                        v1 = np.array([-bar_dir[2], 0.0, bar_dir[0]])
                    v1 /= np.linalg.norm(v1)
                    v2 = np.cross(bar_dir, v1)
                    v2 /= np.linalg.norm(v2)
                    if len(null_dirs) == 0:
                        null_dirs = [v1, v2]
                    elif len(null_dirs) == 1:
                        null_dirs.append(v1 if abs(np.dot(v1, null_dirs[0])) < 0.5 else v2)
                else:
                    # Genuinely non-collinear (span 2D plane): cross product
                    # gives the plane normal (1 nullspace direction)
                    n /= n_norm
                    found = False
                    for nd in null_dirs:
                        if abs(np.dot(nd, n)) > 0.99:
                            found = True
                            break
                    if not found and len(null_dirs) < 2:
                        null_dirs.append(n)

        # Fix DOFs that have significant component in any nullspace direction
        for null_vec in null_dirs:
            for dof in range(3):
                if abs(null_vec[dof]) > 0.5:
                    bar_fix_nodes[nid].add(dof)

        # Bar-only nodes (no Q4, no H8) are internal discretization points
        # of the SupportBeam. Bar elements can only carry axial loads, so
        # nodes with N <= 2 connected bars create mechanisms in the
        # perpendicular directions. Even for N >= 3, if the bars are
        # coplanar, there's a nullspace direction.
        #
        # Strategy: fully fix ALL Bar-only nodes unless they have N >= 3
        # bars that genuinely span 3D space (no nullspace). This prevents
        # chain mechanisms and ensures a solvable stiffness matrix.
        if N <= 2:
            # N=1 or N=2: at most 2D span, fully fix to prevent mechanisms
            bar_fix_nodes[nid] = {0, 1, 2}
        elif len(null_dirs) > 0:
            # N >= 3 but coplanar: fix nullspace directions, AND fully
            # fix to prevent chain mechanisms propagating through the
            # in-plane DOFs
            bar_fix_nodes[nid] = {0, 1, 2}

    # Count bar-fixed nodes
    bar_fix_count = sum(1 for nid, dofs in bar_fix_nodes.items() if dofs)
    print(f"  Bar-plane-fixed nodes (coplanar bars): {bar_fix_count}")
    if bar_fix_count:
        # Show a few examples
        examples = [(nid, dofs) for nid, dofs in bar_fix_nodes.items() if dofs][:5]
        dof_names = {0:'X',1:'Y',2:'Z'}
        for nid, dofs in examples:
            print(f"    Node {nid}: fix DOFs {[dof_names[d] for d in dofs]}")

    # Nodes connected to Q4 + Bar (no H8) may have insufficient Z stiffness
    # if the Bar elements are nearly horizontal (small dz/L).
    # Q4 is plane stress (XY only, no Z stiffness), so Z support comes
    # entirely from Bar elements.
    q4_bar_z_fix = set()
    for nid in connected:
        types = node_elem_types.get(nid, set())
        # Only check nodes with Q4 + Bar, but no H8
        if 2 not in types or 1 not in types or 5 in types:
            continue
        # Check max |dz/L| across all connected Bar elements
        max_dz_ratio = 0.0
        for elem in node_bar_elems.get(nid, []):
            if len(elem['nodes']) != 2:
                continue
            n1, n2 = elem['nodes'][0], elem['nodes'][1]
            other = n2 if n1 == nid else n1
            p_self = np.array(global_nodes[nid])
            p_other = np.array(global_nodes[other])
            vec = p_other - p_self
            L = np.linalg.norm(vec)
            if L > 1e-12:
                dz_ratio = abs(vec[2]) / L
                if dz_ratio > max_dz_ratio:
                    max_dz_ratio = dz_ratio
        # If ALL bars have |dz/L| < 0.05, Z has essentially no stiffness
        if max_dz_ratio < 0.05:
            q4_bar_z_fix.add(nid)

    if q4_bar_z_fix:
        print(f"  Q4+Bar Z-fixed nodes (bars nearly horizontal): {len(q4_bar_z_fix)}")
        for nid in list(q4_bar_z_fix)[:5]:
            print(f"    Node {nid}: max |dz/L| < 0.05")

    num_nodes = len(global_nodes)
    num_eg = len(stap_groups)
    num_loads = 1

    # --- Compute gravity loads ---
    # Abaqus gravity: magnitude * direction
    grav = None
    for ld in step['loads']:
        if ld['type'] == 'gravity':
            grav = ld
            break

    nodal_force = defaultdict(lambda: np.zeros(3))

    if grav:
        mag = grav['magnitude']
        gx, gy, gz = grav['direction']

        for elem in valid_elems:
            etype = elem['type']
            econn = elem['nodes']
            mat_name = elem['material']
            mat = materials.get(mat_name, {})
            rho = mat.get('density', 2320)

            if etype == 'S4R':
                # Shell -> plane stress: compute area and thickness
                coords_xy = [(global_nodes[n][0], global_nodes[n][1]) for n in econn]
                area = polygon_area(coords_xy)
                t = 0.2  # Default floor thickness
                elem_mass = area * t * rho
            elif etype == 'T3D2':
                p1 = np.array(global_nodes[econn[0]])
                p2 = np.array(global_nodes[econn[1]])
                length = np.linalg.norm(p2 - p1)
                sec = elem.get('section')
                sec_data = sec.get('data', []) if sec else []
                A = sec_data[0] if len(sec_data) > 0 else 0.25
                elem_mass = A * length * rho
            elif etype == 'C3D8R':
                # H8 brick: compute volume from 8-node hex
                coords = [np.array(global_nodes[n]) for n in econn]
                # Jacobian at center approximates volume/8
                vol = hex_volume(coords)
                elem_mass = vol * rho
            elif etype == 'B31':
                p1 = np.array(global_nodes[econn[0]])
                p2 = np.array(global_nodes[econn[1]])
                length = np.linalg.norm(p2 - p1)
                sec = elem.get('section')
                sec_data = sec.get('data', []) if sec else []
                if len(sec_data) >= 6:
                    area_beam, _, _, _ = box_section_props(sec_data[0], sec_data[1], sec_data[2], sec_data[3], sec_data[4], sec_data[5])
                else:
                    area_beam = 0.76
                elem_mass = area_beam * length * rho
            else:
                elem_mass = 0.0

            weight = elem_mass * mag
            force_per_node = weight / len(econn)

            for n in econn:
                nodal_force[n][0] += force_per_node * gx
                nodal_force[n][1] += force_per_node * gy
                nodal_force[n][2] += force_per_node * gz

    # Collect non-zero loads
    load_entries = []
    for nid in sorted(nodal_force.keys()):
        fx, fy, fz = nodal_force[nid]
        if abs(fx) > 1e-10:
            load_entries.append((nid, 1, fx))  # DOF 1 = X
        if abs(fy) > 1e-10:
            load_entries.append((nid, 2, fy))  # DOF 2 = Y
        if abs(fz) > 1e-10:
            load_entries.append((nid, 3, fz))  # DOF 3 = Z

    # ============ WRITE ============
    with open(output_path, 'w') as f:
        f.write("Bridge-1: Cable-stayed bridge (3D conversion)\n")
        f.write(f"{num_nodes}  {num_eg}  {num_loads}  1\n")

        # --- Nodal data ---
        for nid in sorted(global_nodes.keys()):
            x, y, z = global_nodes[nid]

            if nid in bc_map:
                bc = bc_map[nid]
                bc_x, bc_y, bc_z, bc_rx, bc_ry, bc_rz = bc[0], bc[1], bc[2], bc[3], bc[4], bc[5]
            elif nid in orphan_nodes:
                bc_x, bc_y, bc_z, bc_rx, bc_ry, bc_rz = 1, 1, 1, 1, 1, 1
            else:
                bc_x, bc_y, bc_z, bc_rx, bc_ry, bc_rz = 0, 0, 0, 1, 1, 1  # rotations fixed by default

            # Nodes connected to Shell4(10) or Beam3D(9) need rotations active
            types = node_elem_types.get(nid, set())
            if 10 in types or 9 in types:
                # Shell4 boundary nodes: fix rotations to prevent spurious modes
                if nid in shell4_boundary_fix:
                    bc_rx = bc_ry = bc_rz = 1
                else:
                    bc_rx = bc_ry = bc_rz = 0

            # Fix Z for nodes only connected to 2D elements (Q4)
            if nid in z_fix_nodes:
                bc_z = 1

            # Fix DOFs for Bar-only nodes with coplanar bar directions
            if nid in bar_fix_nodes:
                dofs = bar_fix_nodes[nid]
                if 0 in dofs: bc_x = 1
                if 1 in dofs: bc_y = 1
                if 2 in dofs: bc_z = 1

            # Fix Z for Q4+Bar nodes where bars are nearly horizontal
            if nid in q4_bar_z_fix:
                bc_z = 1

            f.write(f"  {nid}  {bc_x}  {bc_y}  {bc_z}  {bc_rx}  {bc_ry}  {bc_rz}  {x:.6f}  {y:.6f}  {z:.6f}\n")
        # --- Load data (ALL on ONE line) ---
        nload = len(load_entries)
        f.write(f"  1  {nload}")
        for nid, dof, val in load_entries:
            f.write(f"  {nid}  {dof}  {val:.6e}")
        f.write("\n")

        # --- Element groups ---
        mat_index = {}  # material_name -> local mat set number
        mat_set_counter = 1

        for (stap_type, abaqus_type, mat_name), eid_list in sorted(stap_groups.items()):
            eids = sorted(eid_list)
            nume = len(eids)

            if mat_name not in mat_index:
                mat_index[mat_name] = mat_set_counter
                mat_set_counter += 1

            f.write(f"  {stap_type}  {nume}  1\n")

            # Write material line
            mat = materials.get(mat_name, {})
            E = mat.get('E', 2e11)
            nu = mat.get('nu', 0.3)

            if stap_type in (1,):  # Bar (truss)
                area = 0.25  # Steel cable
                f.write(f"    1  {E:.6e}  {area:.6f} ")


            elif stap_type in (2,):  # Q4 (membrane)
                t = 0.2
                f.write(f"    1  {E:.6e}  {nu:.6f}  {t:.6f} ")


            elif stap_type in (5,):  # H8 (solid)
                f.write(f"    1  {E:.6e}  {nu:.6f} ")


            elif stap_type in (9,):  # Beam3D
                sec_info = None
                for eid in eids:
                    e = eid_to_elem.get(eid)
                    if e and e.get("section"):
                        sec_info = e["section"]
                        break
                if sec_info:
                    sec_data = sec_info.get('data', [])
                    n1_vec = sec_info.get('n1', [0.0, 0.0, -1.0])
                else:
                    sec_data = []
                    n1_vec = [0.0, 0.0, -1.0]
                if len(sec_data) >= 6:
                    a, b, t1, t2, t3, t4 = sec_data[0:6]
                    area, Iy, Iz, J = box_section_props(a, b, t1, t2, t3, t4)
                else:
                    area, Iy, Iz, J = 0.76, 0.4585, 0.4585, 0.6859
                # Validate n1: if zero or parallel to beam axis, fall back
                n1_norm = sum(v*v for v in n1_vec)**0.5
                if n1_norm < 1e-10:
                    n1_vec = [0.0, 0.0, -1.0]
                    area, Iy, Iz, J = 0.76, 0.4585, 0.4585, 0.6859
                f.write(f"    1  {E:.6e}  {nu:.6f}  {area:.6f}  {Iy:.6f}  {Iz:.6f}  {J:.6f}  {n1_vec[0]:.6f}  {n1_vec[1]:.6f}  {n1_vec[2]:.6f} ")


            elif stap_type in (10,):  # Shell4 (flat shell)
                t = 0.2
                f.write(f"    1  {E:.6e}  {nu:.6f}  {t:.6f} ")


            elif stap_type in (2,):  # Q4
                f.write(f"    1  {E:.6e}  {nu:.6f}  0.200000\n")
            elif stap_type in (5,):  # H8
                f.write(f"    1  {E:.6e}  {nu:.6f}\n")

            # Write element connectivity (renumber from 1 within group)
            eid_to_elem = {elem['id']: elem for elem in valid_elems}
            for local_eid, gid in enumerate(eids, start=1):
                elem = eid_to_elem[gid]
                conn = elem['nodes']
                f.write(f"    {local_eid}")
                for n in conn:
                    f.write(f"  {n}")
                f.write(f"  1\n")

        # --- End marker ---
        f.write("  0\n")

    return output_path


def hex_volume(coords):
    """Approximate volume of an 8-node hex element using Jacobian at center."""
    # Shape function derivatives at origin
    dN = [
        [-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
        [-1, -1,  1], [1, -1,  1], [1, 1,  1], [-1, 1,  1]
    ]
    dN = [[x/8.0 for x in d] for d in dN]

    J = np.zeros((3, 3))
    for a in range(8):
        x, y, z = coords[a]
        for i in range(3):
            J[i][0] += dN[a][0] * (x if i == 0 else y if i == 1 else z)
            J[i][1] += dN[a][1] * (x if i == 0 else y if i == 1 else z)
            J[i][2] += dN[a][2] * (x if i == 0 else y if i == 1 else z)

    detJ = abs(np.linalg.det(J))
    return detJ * 8.0  # Gauss quadrature weight sum = 8


# ============================================================
# Main
# ============================================================

def main():
    import sys
    inp_file = sys.argv[1] if len(sys.argv) > 1 else 'Bridge-1.inp'
    out_file = sys.argv[2] if len(sys.argv) > 2 else 'Bridge-1.dat'

    print("=" * 60)
    print("Abaqus -> STAP++ 3D Converter")
    print("=" * 60)

    print(f"\nParsing {inp_file}...")
    data = parse_abaqus_inp(inp_file)

    print(f"  Parts: {list(data['parts'].keys())}")
    print(f"  Instances: {len(data['assembly']['instances'])}")
    print(f"  Materials: {list(data['materials'].keys())}")
    print(f"  BC sets: {len(data['step']['boundary'])}")
    print(f"  Loads: {len(data['step']['loads'])}")

    print(f"\nFlattening assembly (3D transforms, no flattening)...")
    global_nodes, global_elements, element_groups, node_map, materials = \
        flatten_assembly(data)

    # Count element types
    type_counts = defaultdict(int)
    for elem in global_elements:
        type_counts[elem['type']] += 1
    print(f"  Global nodes: {len(global_nodes)}")
    print(f"  Global elements: {len(global_elements)}")
    print(f"  Element types: {dict(type_counts)}")

    print(f"\nWriting STAP++ .dat to {out_file}...")
    write_stap_dat(out_file, data, global_nodes, global_elements,
                   element_groups, node_map)

    print(f"Written {out_file}")
    print("Done!")


if __name__ == '__main__':
    main()
