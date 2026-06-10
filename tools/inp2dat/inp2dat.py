#!/usr/bin/env python3
import argparse
import math
from collections import defaultdict


def trim(s):
    return s.strip()


def split_comma(s):
    return [x.strip() for x in s.split(",")]


def header_value(header, key):
    prefix = key.lower() + "="
    for part in split_comma(header):
        if part.lower().startswith(prefix):
            return part[len(prefix):].strip()
    return ""


def header_has_flag(header, flag):
    return any(part.strip().lower() == flag.lower() for part in split_comma(header))


def parse_floats(line):
    return [float(x) for x in split_comma(line) if x]


def parse_set_lines(lines, i, generate=False):
    nums = []
    while i < len(lines) and not lines[i].startswith("*"):
        for tok in split_comma(lines[i]):
            if tok:
                nums.append(int(tok))
        i += 1
    if generate and len(nums) == 3:
        a, b, step = nums
        nums = list(range(a, b + 1, step))
    return nums, i


def vec_add(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def vec_sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def vec_dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def vec_cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def vec_norm(a):
    return math.sqrt(vec_dot(a, a))


def vec_scale(a, s):
    return (a[0] * s, a[1] * s, a[2] * s)


def vec_normalize(a):
    n = vec_norm(a)
    if n < 1.0e-15:
        return (0.0, 0.0, 0.0)
    return vec_scale(a, 1.0 / n)


def rotate_rodrigues(v, axis_a, axis_b, angle_deg):
    u = vec_normalize(vec_sub(axis_b, axis_a))
    if vec_norm(u) < 1.0e-15:
        return v
    theta = math.radians(angle_deg)
    c = math.cos(theta)
    s = math.sin(theta)
    uxv = vec_cross(u, v)
    udotv = vec_dot(u, v)
    return vec_add(vec_add(vec_scale(v, c), vec_scale(uxv, s)), vec_scale(u, udotv * (1.0 - c)))


def polygon_area_xy(conn, nodes):
    area = 0.0
    for i, nid in enumerate(conn):
        a = nodes[nid]
        b = nodes[conn[(i + 1) % len(conn)]]
        area += a[0] * b[1] - b[0] * a[1]
    return abs(area) * 0.5


def det3(m):
    return (
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    )


def hex_volume(conn, nodes):
    signs = [
        (-1, -1, -1), (1, -1, -1), (1, 1, -1), (-1, 1, -1),
        (-1, -1, 1), (1, -1, 1), (1, 1, 1), (-1, 1, 1),
    ]
    j = [[0.0, 0.0, 0.0] for _ in range(3)]
    for a, nid in enumerate(conn[:8]):
        x, y, z = nodes[nid]
        for col in range(3):
            d = signs[a][col] / 8.0
            j[0][col] += d * x
            j[1][col] += d * y
            j[2][col] += d * z
    return abs(det3(j)) * 8.0


def parse_part(lines, i):
    part = {"nodes": {}, "elements": [], "sections": {}, "nsets": {}, "elsets": {}}
    while i < len(lines):
        line = lines[i]
        if line.startswith("*End Part"):
            return part, i + 1
        if line.startswith("*Node"):
            i += 1
            while i < len(lines) and not lines[i].startswith("*"):
                p = split_comma(lines[i])
                if len(p) >= 4:
                    part["nodes"][int(p[0])] = (float(p[1]), float(p[2]), float(p[3]))
                i += 1
            continue
        if line.startswith("*Element"):
            etype = header_value(line, "type")
            i += 1
            while i < len(lines) and not lines[i].startswith("*"):
                p = split_comma(lines[i])
                if len(p) >= 2:
                    part["elements"].append({
                        "id": int(p[0]),
                        "type": etype,
                        "nodes": [int(x) for x in p[1:] if x],
                    })
                i += 1
            continue
        if line.startswith("*Nset"):
            name = header_value(line, "nset")
            nums, i = parse_set_lines(lines, i + 1, header_has_flag(line, "generate"))
            part["nsets"][name] = nums
            continue
        if line.startswith("*Elset"):
            name = header_value(line, "elset")
            nums, i = parse_set_lines(lines, i + 1, header_has_flag(line, "generate"))
            part["elsets"][name] = nums
            continue
        if line.startswith("*Solid Section") or line.startswith("*Shell Section") or line.startswith("*Beam Section"):
            elset = header_value(line, "elset")
            section = {"material": header_value(line, "material"), "data": [], "n1": None}
            is_beam = line.startswith("*Beam Section")
            i += 1
            if i < len(lines) and not lines[i].startswith("*"):
                section["data"] = parse_floats(lines[i])
                i += 1
                if is_beam and i < len(lines) and not lines[i].startswith("*"):
                    try:
                        n1_vals = parse_floats(lines[i])
                        if len(n1_vals) >= 3:
                            section["n1"] = tuple(n1_vals[:3])
                        i += 1
                    except ValueError:
                        pass
            part["sections"][elset] = section
            continue
        i += 1
    return part, i


def parse_assembly_set(lines, i, key):
    header = lines[i]
    name = header_value(header, key)
    instance = header_value(header, "instance")
    nums, i = parse_set_lines(lines, i + 1, header_has_flag(header, "generate"))
    return {"name": name, "instance": instance, "numbers": nums}, i


def parse_assembly(lines, i):
    assembly = {"instances": [], "nsets": [], "elsets": [], "ties": [], "surfaces": {}}
    while i < len(lines):
        line = lines[i]
        if line.startswith("*End Assembly"):
            return assembly, i + 1
        if line.startswith("*Instance"):
            inst = {
                "name": header_value(line, "name"),
                "part": header_value(line, "part"),
                "translation": (0.0, 0.0, 0.0),
                "rotation": None,
            }
            i += 1
            if i < len(lines) and not lines[i].startswith("*"):
                vals = parse_floats(lines[i])
                if len(vals) == 3:
                    inst["translation"] = tuple(vals)
                    i += 1
            if i < len(lines) and not lines[i].startswith("*"):
                vals = parse_floats(lines[i])
                if len(vals) == 7:
                    inst["rotation"] = (tuple(vals[:3]), tuple(vals[3:6]), vals[6])
                    i += 1
            assembly["instances"].append(inst)
            continue
        if line.startswith("*Nset"):
            item, i = parse_assembly_set(lines, i, "nset")
            assembly["nsets"].append(item)
            continue
        if line.startswith("*Elset"):
            item, i = parse_assembly_set(lines, i, "elset")
            assembly["elsets"].append(item)
            continue
        if line.startswith("*Surface"):
            sname = header_value(line, "name")
            i += 1
            if i < len(lines) and not lines[i].startswith("*"):
                assembly["surfaces"][sname] = lines[i].strip().split(",")[0].strip()
                i += 1
            continue
        if line.startswith("*Tie"):
            tname = header_value(line, "name")
            i += 1
            tie = {"name": tname}
            while i < len(lines) and not lines[i].startswith("*"):
                p = split_comma(lines[i])
                if len(p) >= 2:
                    tie["slave"] = p[0].strip()
                    tie["master"] = p[1].strip()
                i += 1
            assembly["ties"].append(tie)
            continue
        i += 1
    return assembly, i


def parse_material(lines, i):
    name = header_value(lines[i], "name")
    mat = {"density": 2320.0, "E": 2.0e11, "nu": 0.3}
    i += 1
    while i < len(lines):
        line = lines[i]
        if line.startswith("*Material") or line.startswith("*Part") or line.startswith("*Assembly") or line.startswith("*Step"):
            break
        if line.startswith("*Density"):
            i += 1
            if i < len(lines):
                vals = parse_floats(lines[i])
                if vals:
                    mat["density"] = vals[0]
            i += 1
            continue
        if line.startswith("*Elastic"):
            i += 1
            if i < len(lines):
                vals = parse_floats(lines[i])
                if len(vals) >= 2:
                    mat["E"] = vals[0]
                    mat["nu"] = vals[1]
            i += 1
            continue
        i += 1
    return name, mat, i


def parse_step(lines, i):
    boundaries = []
    gravity = []
    while i < len(lines):
        line = lines[i]
        if line.startswith("*End Step"):
            return boundaries, gravity, i + 1
        if line.startswith("*Boundary"):
            i += 1
            while i < len(lines) and not lines[i].startswith("*"):
                p = split_comma(lines[i])
                if len(p) >= 3:
                    boundaries.append({"set": p[0], "first": int(p[1]), "last": int(p[2])})
                i += 1
            continue
        if line.startswith("*Dload"):
            i += 1
            while i < len(lines) and not lines[i].startswith("*"):
                p = split_comma(lines[i])
                if "GRAV" in p:
                    k = p.index("GRAV")
                    gravity.append({
                        "magnitude": float(p[k + 1]),
                        "direction": (float(p[k + 2]), float(p[k + 3]), float(p[k + 4])),
                    })
                i += 1
            continue
        i += 1
    return boundaries, gravity, i


def parse_inp(path):
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        lines = [trim(line) for line in f if not trim(line).startswith("**")]

    data = {"parts": {}, "assembly": {}, "materials": {}, "boundaries": [], "gravity": []}
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith("*Part"):
            name = header_value(line, "name")
            part, i = parse_part(lines, i + 1)
            data["parts"][name] = part
            continue
        if line.startswith("*Assembly"):
            data["assembly"], i = parse_assembly(lines, i + 1)
            continue
        if line.startswith("*Material"):
            name, mat, i = parse_material(lines, i)
            data["materials"][name] = mat
            continue
        if line.startswith("*Step"):
            bcs, grav, i = parse_step(lines, i + 1)
            data["boundaries"].extend(bcs)
            data["gravity"].extend(grav)
            continue
        i += 1
    return data


def find_section(part, elem_id):
    for elset_name, elem_ids in part["elsets"].items():
        if elem_id not in elem_ids:
            continue
        for sec_name, sec in part["sections"].items():
            short_name = sec_name.rsplit(".", 1)[-1]
            if sec_name == elset_name or short_name == elset_name:
                return sec
    return {"material": "", "data": [], "n1": None}


def flatten_assembly(data, merge_tol=1.0e-6):
    raw_nodes = {}
    raw_elements = []
    node_map_raw = {}
    next_node = 1
    next_elem = 1

    for inst in data["assembly"].get("instances", []):
        part = data["parts"].get(inst["part"])
        if part is None:
            print(f"  WARNING: part '{inst['part']}' not found for instance '{inst['name']}'")
            continue

        local_to_global = {}
        for local_id, coord in part["nodes"].items():
            p = coord
            if inst["rotation"]:
                axis_a, axis_b, angle = inst["rotation"]
                p = rotate_rodrigues(p, axis_a, axis_b, angle)
            p = vec_add(p, inst["translation"])
            gid = next_node
            next_node += 1
            raw_nodes[gid] = p
            local_to_global[local_id] = gid
            node_map_raw[(inst["name"], local_id)] = gid

        for src in part["elements"]:
            sec = find_section(part, src["id"])
            raw_elements.append({
                "id": next_elem,
                "type": src["type"],
                "nodes": [local_to_global[n] for n in src["nodes"]],
                "material": sec.get("material", ""),
                "section": sec.get("data", []),
                "n1": sec.get("n1"),
            })
            next_elem += 1

    coord_to_node = {}
    old_to_new = {}
    nodes = {}
    next_new = 1
    for old_id, p in sorted(raw_nodes.items()):
        key = tuple(round(x / merge_tol) for x in p)
        if key not in coord_to_node:
            coord_to_node[key] = next_new
            old_to_new[old_id] = next_new
            nodes[next_new] = p
            next_new += 1
        else:
            old_to_new[old_id] = coord_to_node[key]

    if len(nodes) != len(raw_nodes):
        print(f"  Coincident nodes merged: {len(raw_nodes) - len(nodes)}")

    for elem in raw_elements:
        elem["nodes"] = [old_to_new[n] for n in elem["nodes"]]

    node_map = {key: old_to_new[old] for key, old in node_map_raw.items()}

    # Process tie constraints: merge slave nodes into master nodes
    slave_to_master = {}
    for tie in data["assembly"].get("ties", []):
        slave_surf = tie.get("slave", "")
        master_surf = tie.get("master", "")
        slave_nset = data["assembly"]["surfaces"].get(slave_surf, slave_surf.replace("_CNS_", ""))
        master_nset = data["assembly"]["surfaces"].get(master_surf, master_surf.replace("_CNS_", ""))
        slave_nodes, master_nodes = [], []
        for nset in data["assembly"].get("nsets", []):
            if nset["name"] == slave_nset:
                for lnid in nset["numbers"]:
                    gid = node_map.get((nset["instance"], lnid))
                    if gid: slave_nodes.append(gid)
            if nset["name"] == master_nset:
                for lnid in nset["numbers"]:
                    gid = node_map.get((nset["instance"], lnid))
                    if gid: master_nodes.append(gid)
        if not slave_nodes or not master_nodes:
            continue
        import math
        used = set()
        for sg in slave_nodes:
            sp = raw_nodes[sg]
            best_mg, best_d = None, float("inf")
            for mg in master_nodes:
                if mg in used:
                    continue
                mp = raw_nodes[mg]
                d = math.sqrt((sp[0]-mp[0])**2 + (sp[1]-mp[1])**2 + (sp[2]-mp[2])**2)
                if d < best_d:
                    best_d = d
                    best_mg = mg
            if best_mg and best_d < 100.0:
                used.add(best_mg)
                slave_to_master[sg] = best_mg
    if slave_to_master:
        print(f"  Tie constraints: merging {len(slave_to_master)} slave nodes")
        for elem in raw_elements:
            elem["nodes"] = [slave_to_master.get(n, n) for n in elem["nodes"]]

    # Fix cable elements: ensure both end nodes have the same Y sign.
    # Tie constraint processing or coordinate merging can occasionally
    # connect a cable tower-end to the wrong tower (Y sign mismatch).
    cable_bad = 0
    for elem in raw_elements:
        if elem["type"] != "T3D2":
            continue
        if len(elem["nodes"]) != 2:
            continue
        n1, n2 = elem["nodes"]
        p1 = nodes.get(n1, (0, 0, 0))
        p2 = nodes.get(n2, (0, 0, 0))
        if abs(p1[1] - p2[1]) < 0.01:
            continue  # Y matches, OK
        # Find correct tower node: same Z (≈150), sign of Y matches deck end
        deck_n = n1 if abs(p1[2]) < 1 else n2
        tower_n = n2 if deck_n == n1 else n1
        deck_p = nodes.get(deck_n, (0, 0, 0))
        target_y = deck_p[1]  # desired Y sign for tower end
        best_tower, best_dist = None, 1e9
        for tid, tp in nodes.items():
            if abs(tp[2] - abs(deck_p[2]) - 150) > 5:  # only consider tower-top-level nodes
                continue
            if abs(tp[0]) > 5:  # tower top X ≈ 0
                continue
            # Must have same Y sign and close to expected position
            dy = abs(tp[1] - target_y)
            if dy < best_dist:
                best_dist = dy
                best_tower = tid
        if best_tower and best_dist < 0.1:
            if deck_n == n1:
                elem["nodes"] = [n1, best_tower]
            else:
                elem["nodes"] = [best_tower, n2]
            cable_bad += 1
    if cable_bad:
        print(f"  Cable Y-consistency fix: {cable_bad} elements corrected")
    return nodes, raw_elements, node_map


def stap_type(etype):
    return {"T3D2": 1, "S4R": 10, "C3D8R": 5, "B31": 9}.get(etype, 0)


THREE_D_TYPES = {1, 5, 9, 10}


def box_section_props(a, b, t1, t2, t3, t4):
    """Compute cross-section properties for rectangular box."""
    tw, tf = t1 + t2, t3 + t4
    area = a * b - (a - tw) * (b - tf)
    Iy = (b * a**3 - (b - tf) * (a - tw)**3) / 12.0
    Iz = (a * b**3 - (a - tw) * (b - tf)**3) / 12.0
    a0, b0 = a - tw / 2.0, b - tf / 2.0
    A0 = a0 * b0
    sum_st = a0 / t1 + a0 / t2 + b0 / t3 + b0 / t4
    J = 4.0 * A0**2 / sum_st if sum_st > 1e-15 else 0.0
    return area, Iy, Iz, J


def resolve_fixed_nodes(data, node_map):
    fixed = defaultdict(lambda: [0, 0, 0, 1, 1, 1])
    for bc in data["boundaries"]:
        for nset in data["assembly"].get("nsets", []):
            if nset["name"] != bc["set"]:
                continue
            for local_id in nset["numbers"]:
                gid = node_map.get((nset["instance"], local_id))
                if not gid:
                    continue
                for dof in range(bc["first"], min(bc["last"], 3) + 1):
                    fixed[gid][dof - 1] = 1
    return dict(fixed)


def write_dat(path, data, nodes, elements, node_map):
    valid = [e for e in elements if stap_type(e["type"])]
    skipped = defaultdict(int)
    for e in elements:
        if not stap_type(e["type"]):
            skipped[e["type"]] += 1
    for etype, count in skipped.items():
        print(f"  WARNING: unknown element type {etype}, skipping {count} elements")

    fixed = resolve_fixed_nodes(data, node_map)
    print(f"  Fixed nodes (from BCs): {len(fixed)}")

    connected = set()
    node_types = defaultdict(set)
    node_bar_elems = defaultdict(list)
    node_shell4_count = defaultdict(int)
    node_beam3d_count = defaultdict(int)
    for e in valid:
        st = stap_type(e["type"])
        for n in e["nodes"]:
            connected.add(n)
            node_types[n].add(st)
            if st == 1:
                node_bar_elems[n].append(e)
            if st == 10:  # Shell4
                node_shell4_count[n] += 1
            if st == 9:  # Beam3D
                node_beam3d_count[n] += 1

    orphan_nodes = set(nodes) - connected
    if orphan_nodes:
        print(f"  Orphan nodes (will be fully fixed): {len(orphan_nodes)}")

    # Shell4 boundary: nodes with <=2 connected Shell4 OR at deck ends (|X| > 220)
    # get rotations fixed. This prevents cantilever-like behavior at the
    # unsupported deck extremities where cables do not reach.
    shell4_boundary_fix = set()
    for nid, count in node_shell4_count.items():
        if 10 not in node_types.get(nid, set()):
            continue
        x, y, z = nodes.get(nid, (0, 0, 0))
        if count <= 2 or abs(x) > 220:
            shell4_boundary_fix.add(nid)
    print(f"  Shell4 boundary nodes (rotation-fixed): {len(shell4_boundary_fix)}")

    # H8-Shell4/Beam3D interface: H8 has only 3 translational DOFs, so at shared
    # nodes the rotational DOFs have no H8 stiffness. Fix drilling rotation to
    # prevent artificial hinge behavior at pier-deck and pier-beam connections.
    h8_interface_rz_fix = set()
    h8_interface_rxyz_fix = set()
    for nid, types in node_types.items():
        if 5 not in types:  # Must involve H8
            continue
        if 10 in types or 9 in types:  # H8 + Shell4/Beam3D
            if 10 in types:
                h8_interface_rz_fix.add(nid)  # Shell4: fix rz (drilling)
            if 9 in types:
                h8_interface_rxyz_fix.add(nid)  # Beam3D: fix all rotations at H8 end
    print(f"  H8-Shell4/Beam3D interface nodes (rz-fixed): {len(h8_interface_rz_fix)}")
    print(f"  H8-Beam3D interface nodes (rxyz-fixed): {len(h8_interface_rxyz_fix)}")

    # Bar-only nodes: fix DOFs for mechanism prevention
    # (Beam3D nodes excluded - they have full bending/torsional stiffness)
    bar_fix_nodes = {}
    for nid, elems in node_bar_elems.items():
        types = node_types[nid]
        if types != {1}:  # Only Bar, no other types
            continue
        if len(elems) <= 2:
            bar_fix_nodes[nid] = {0, 1, 2}
            continue
        dirs = []
        for e in elems:
            if len(e["nodes"]) != 2:
                continue
            other = e["nodes"][1] if e["nodes"][0] == nid else e["nodes"][0]
            d = vec_normalize(vec_sub(nodes[other], nodes[nid]))
            if vec_norm(d) > 1e-12:
                dirs.append(d)
        has_non_coplanar = False
        for a in range(len(dirs)):
            if has_non_coplanar:
                break
            for b in range(a + 1, len(dirs)):
                normal = vec_cross(dirs[a], dirs[b])
                if vec_norm(normal) < 1.0e-6:
                    continue
                normal = vec_normalize(normal)
                for c in range(b + 1, len(dirs)):
                    if abs(vec_dot(normal, dirs[c])) > 1.0e-3:
                        has_non_coplanar = True
                        break
                if has_non_coplanar:
                    break
        if not has_non_coplanar:
            bar_fix_nodes[nid] = {0, 1, 2}
    print(f"  Bar-plane-fixed nodes: {len(bar_fix_nodes)}")

    # Gravity loads
    nodal_force = defaultdict(lambda: [0.0, 0.0, 0.0])
    if data["gravity"]:
        grav = data["gravity"][0]
        mag = grav["magnitude"]
        gx, gy, gz = grav["direction"]
        for e in valid:
            mat = data["materials"].get(e["material"], {"density": 2320.0})
            rho = mat.get("density", 2320.0)
            mass = 0.0
            if e["type"] == "S4R":
                sec = e.get("section", [])
                t_shell = sec[0] if sec else 0.2
                mass = polygon_area_xy(e["nodes"], nodes) * t_shell * rho
            elif e["type"] == "T3D2":
                length = vec_norm(vec_sub(nodes[e["nodes"][1]], nodes[e["nodes"][0]]))
                area = e["section"][0] if e["section"] else 0.25
                mass = area * length * rho
            elif e["type"] == "B31":
                length = vec_norm(vec_sub(nodes[e["nodes"][1]], nodes[e["nodes"][0]]))
                sec = e.get("section", [])
                if len(sec) >= 6:
                    area_b, _, _, _ = box_section_props(sec[0], sec[1], sec[2], sec[3], sec[4], sec[5])
                else:
                    area_b = 0.76
                mass = area_b * length * rho
            elif e["type"] == "C3D8R":
                mass = hex_volume(e["nodes"], nodes) * rho
            per_node = mass * mag / max(1, len(e["nodes"]))
            for n in e["nodes"]:
                nodal_force[n][0] += per_node * gx
                nodal_force[n][1] += per_node * gy
                nodal_force[n][2] += per_node * gz

    load_entries = []
    for nid in sorted(nodal_force):
        for dof, val in enumerate(nodal_force[nid], start=1):
            if abs(val) > 1.0e-10:
                load_entries.append((nid, dof, val))

    groups = defaultdict(list)
    for e in valid:
        groups[(stap_type(e["type"]), e["type"], e["material"])].append(e)

    with open(path, "w", encoding="utf-8") as f:
        f.write("Bridge-1: Cable-stayed bridge (3D conversion)\n")
        f.write(f"{len(nodes)}  {len(groups)}  1  1\n")
        for nid in sorted(nodes):
            # Start with NDF=6 bc: translations from fixed set, rotations default=1
            bc = list(fixed.get(nid, [0, 0, 0, 1, 1, 1]))
            if nid in orphan_nodes:
                bc = [1, 1, 1, 1, 1, 1]
            if nid in bar_fix_nodes:
                for dof in bar_fix_nodes[nid]:
                    bc[dof] = 1
            # Shell4/Beam3D nodes: activate rotations
            types = node_types.get(nid, set())
            if 10 in types or 9 in types:
                if nid in shell4_boundary_fix and 9 not in types:
                    # Shell4-only boundary: fix only drilling (rz), keep bending (rx,ry) active
                    bc[3] = 0
                    bc[4] = 0
                    bc[5] = 1
                elif nid in h8_interface_rxyz_fix:
                    # Beam3D at H8 interface: fix all rotations
                    bc[3] = bc[4] = bc[5] = 1
                elif nid in h8_interface_rz_fix:
                    # Shell4 at H8 interface: fix only drilling (rz)
                    bc[3] = 0
                    bc[4] = 0
                    bc[5] = 1
                else:
                    # Interior Shell4 or any Beam3D node: all rotations active
                    bc[3] = bc[4] = bc[5] = 0
            x, y, z = nodes[nid]
            f.write(f"  {nid}  {bc[0]}  {bc[1]}  {bc[2]}  {bc[3]}  {bc[4]}  {bc[5]}  {x:.6f}  {y:.6f}  {z:.6f}\n")

        f.write(f"  1  {len(load_entries)}")
        for nid, dof, val in load_entries:
            f.write(f"  {nid}  {dof}  {val:.6e}")
        f.write("\n")

        for (st, abaqus_type, mat_name), elems in sorted(groups.items()):
            mat = data["materials"].get(mat_name, {"E": 2.0e11, "nu": 0.3})
            f.write(f"  {st}  {len(elems)}  1\n")
            if st == 1:  # Bar
                area = 0.25
                if abaqus_type == "B31":
                    area = 0.76
                elif abaqus_type == "T3D2" and elems:
                    sec = elems[0].get("section", [])
                    area = sec[0] if sec else 0.25
                f.write(f"    1  {mat.get('E', 2.0e11):.6e}  {area:.6f}\n")
            elif st == 2:  # Q4
                f.write(f"    1  {mat.get('E', 2.0e11):.6e}  {mat.get('nu', 0.3):.6f}  0.200000\n")
            elif st == 5:  # H8
                f.write(f"    1  {mat.get('E', 2.0e11):.6e}  {mat.get('nu', 0.3):.6f}\n")
            elif st == 9:  # Beam3D
                E_val = mat.get("E", 2.0e11)
                nu_val = mat.get("nu", 0.3)
                # Get section data from first element
                sec_info = None
                for e in elems:
                    s = e.get("section", [])
                    if len(s) >= 6:
                        sec_info = s
                        break
                if sec_info:
                    a, b, t1, t2, t3, t4 = sec_info[:6]
                    area, Iy, Iz, J = box_section_props(a, b, t1, t2, t3, t4)
                else:
                    area, Iy, Iz, J = 0.76, 0.4585, 0.4585, 0.6859
                # Get n1 vector
                n1_vec = elems[0].get("n1") if elems else None
                if n1_vec is None:
                    n1_vec = (0.0, 0.0, -1.0)
                f.write(f"    1  {E_val:.6e}  {nu_val:.6f}  {area:.6f}  {Iy:.6f}  {Iz:.6f}  {J:.6f}  {n1_vec[0]:.6f}  {n1_vec[1]:.6f}  {n1_vec[2]:.6f}\n")
            elif st == 10:  # Shell4
                sec_data = elems[0].get("section", [])
                t = sec_data[0] if sec_data else 0.2
                f.write(f"    1  {mat.get('E', 2.0e11):.6e}  {mat.get('nu', 0.3):.6f}  {t:.6f}\n")
            for local_id, elem in enumerate(elems, start=1):
                conn = "  ".join(str(n) for n in elem["nodes"])
                f.write(f"    {local_id}  {conn}  1\n")
        f.write("  0\n")


def main():
    parser = argparse.ArgumentParser(description="Convert Abaqus .inp to STAP++ .dat")
    parser.add_argument("inp", nargs="?", default="Bridge-1.inp")
    parser.add_argument("dat", nargs="?", default="Bridge-1.dat")
    args = parser.parse_args()

    print("Abaqus .inp -> STAP++ .dat converter")
    print(f"  Input : {args.inp}")
    print(f"  Output: {args.dat}")
    data = parse_inp(args.inp)
    print(f"  Parts: {len(data['parts'])}")
    print(f"  Instances: {len(data['assembly'].get('instances', []))}")
    print(f"  Materials: {len(data['materials'])}")
    nodes, elements, node_map = flatten_assembly(data)
    print(f"  Global nodes: {len(nodes)}")
    print(f"  Global elements: {len(elements)}")
    write_dat(args.dat, data, nodes, elements, node_map)
    print("Done.")


if __name__ == "__main__":
    main()
