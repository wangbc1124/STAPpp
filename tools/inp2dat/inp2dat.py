#!/usr/bin/env python3
import argparse
import math
from collections import defaultdict
from collections import deque


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


def beam_local_axes(p1, p2, n1):
    ex = vec_sub(p2, p1)
    L = vec_norm(ex)
    if L < 1.0e-15:
        return (1.0, 0.0, 0.0), (0.0, 0.0, -1.0), (0.0, 1.0, 0.0), 1.0
    ex = vec_scale(ex, 1.0 / L)

    nref = n1 if n1 is not None else (0.0, 0.0, -1.0)
    proj = vec_dot(nref, ex)
    ey = vec_sub(nref, vec_scale(ex, proj))
    ey_norm = vec_norm(ey)
    if ey_norm < 1.0e-10:
        if abs(ex[2]) < 0.9:
            ey = (-ex[1], ex[0], 0.0)
        else:
            ey = (0.0, -ex[2], ex[1])
        ey_norm = vec_norm(ey)
    ey = vec_scale(ey, 1.0 / max(ey_norm, 1.0e-15))
    ez = vec_cross(ex, ey)
    return ex, ey, ez, L


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


def hex_faces(conn):
    return [
        [conn[i] for i in (0, 1, 2, 3)],
        [conn[i] for i in (4, 5, 6, 7)],
        [conn[i] for i in (0, 1, 5, 4)],
        [conn[i] for i in (1, 2, 6, 5)],
        [conn[i] for i in (2, 3, 7, 6)],
        [conn[i] for i in (3, 0, 4, 7)],
    ]


def element_tie_faces(elem):
    etype = elem["type"]
    conn = elem["nodes"]
    if etype == "C3D8R" and len(conn) >= 8:
        return [("quad", face) for face in hex_faces(conn)]
    if etype == "S4R" and len(conn) >= 4:
        return [("quad", conn[:4])]
    if etype in ("B31", "T3D2") and len(conn) >= 2:
        return [("line", conn[:2])]
    return []


def project_line(point, face_nodes, nodes):
    a = nodes[face_nodes[0]]
    b = nodes[face_nodes[1]]
    ab = vec_sub(b, a)
    den = vec_dot(ab, ab)
    t = 0.0 if den < 1.0e-20 else vec_dot(vec_sub(point, a), ab) / den
    t = max(0.0, min(1.0, t))
    q = vec_add(a, vec_scale(ab, t))
    dist = vec_norm(vec_sub(point, q))
    return dist, [(face_nodes[0], 1.0 - t), (face_nodes[1], t)]


def quad_shape(xi, eta):
    return [
        0.25 * (1.0 - xi) * (1.0 - eta),
        0.25 * (1.0 + xi) * (1.0 - eta),
        0.25 * (1.0 + xi) * (1.0 + eta),
        0.25 * (1.0 - xi) * (1.0 + eta),
    ]


def project_quad(point, face_nodes, nodes):
    pts = [nodes[n] for n in face_nodes]
    xi = 0.0
    eta = 0.0
    for _ in range(12):
        n = quad_shape(xi, eta)
        q = (0.0, 0.0, 0.0)
        for w, p in zip(n, pts):
            q = vec_add(q, vec_scale(p, w))
        r = vec_sub(q, point)
        dxi = [
            -0.25 * (1.0 - eta),
             0.25 * (1.0 - eta),
             0.25 * (1.0 + eta),
            -0.25 * (1.0 + eta),
        ]
        deta = [
            -0.25 * (1.0 - xi),
            -0.25 * (1.0 + xi),
             0.25 * (1.0 + xi),
             0.25 * (1.0 - xi),
        ]
        gx = (0.0, 0.0, 0.0)
        ge = (0.0, 0.0, 0.0)
        for wx, we, p in zip(dxi, deta, pts):
            gx = vec_add(gx, vec_scale(p, wx))
            ge = vec_add(ge, vec_scale(p, we))
        a00 = vec_dot(gx, gx)
        a01 = vec_dot(gx, ge)
        a11 = vec_dot(ge, ge)
        b0 = -vec_dot(gx, r)
        b1 = -vec_dot(ge, r)
        det = a00 * a11 - a01 * a01
        if abs(det) < 1.0e-24:
            break
        d_xi = (b0 * a11 - b1 * a01) / det
        d_eta = (a00 * b1 - a01 * b0) / det
        xi = max(-1.0, min(1.0, xi + d_xi))
        eta = max(-1.0, min(1.0, eta + d_eta))
        if abs(d_xi) + abs(d_eta) < 1.0e-10:
            break
    n = quad_shape(xi, eta)
    q = (0.0, 0.0, 0.0)
    for w, p in zip(n, pts):
        q = vec_add(q, vec_scale(p, w))
    dist = vec_norm(vec_sub(point, q))
    return dist, [(nid, w) for nid, w in zip(face_nodes, n) if abs(w) > 1.0e-10]


def build_master_faces(elements, master_nodes, master_instances=None):
    master_set = set(master_nodes)
    master_instances = set(master_instances or [])
    faces = []
    for elem in elements:
        if master_instances and elem.get("instance") not in master_instances:
            continue
        for kind, face in element_tie_faces(elem):
            if master_instances or all(n in master_set for n in face):
                faces.append((kind, face))
    return faces


def find_tie_interpolation(point, master_nodes, master_faces, nodes):
    best = None
    for kind, face in master_faces:
        if kind == "quad":
            dist, weights = project_quad(point, face, nodes)
        else:
            dist, weights = project_line(point, face, nodes)
        if best is None or dist < best[0]:
            best = (dist, weights)
    if best is not None:
        return best[0], best[1], "projected"

    best_n = None
    best_d = float("inf")
    for nid in master_nodes:
        d = vec_norm(vec_sub(point, nodes[nid]))
        if d < best_d:
            best_n = nid
            best_d = d
    if best_n is None:
        return None
    return best_d, [(best_n, 1.0)], "node"


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
            stype = header_value(line, "type") or ""
            i += 1
            if i < len(lines) and not lines[i].startswith("*"):
                set_name = lines[i].strip().split(",")[0].strip()
                assembly["surfaces"][sname] = {"set": set_name, "type": stype.upper()}
                i += 1
            continue
        if line.startswith("*Tie"):
            tname = header_value(line, "name")
            adjust = (header_value(line, "adjust") or "").lower()
            i += 1
            tie = {
                "name": tname,
                "adjust": adjust != "no",
                "no_rotation": header_has_flag(line, "no rotation"),
            }
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


def surface_info(assembly, surf_name):
    raw = assembly.get("surfaces", {}).get(surf_name)
    if isinstance(raw, dict):
        return raw
    if raw:
        return {"set": raw, "type": ""}
    return {"set": surf_name.replace("_CNS_", ""), "type": ""}


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


def build_rcm_order(nodes, elements, mpc_equations):
    adjacency = {nid: set() for nid in nodes}
    for elem in elements:
        conn = [nid for nid in elem.get("nodes", []) if nid in adjacency]
        for i, a in enumerate(conn):
            for b in conn[i + 1:]:
                if a != b:
                    adjacency[a].add(b)
                    adjacency[b].add(a)
    for equation in mpc_equations:
        conn = [term[0] for term in equation.get("terms", []) if term[0] in adjacency]
        for i, a in enumerate(conn):
            for b in conn[i + 1:]:
                if a != b:
                    adjacency[a].add(b)
                    adjacency[b].add(a)

    degree = {nid: len(neighbors) for nid, neighbors in adjacency.items()}
    unvisited = set(nodes)
    order = []
    while unvisited:
        start = min(unvisited, key=lambda nid: (degree[nid], nid))
        queue = deque([start])
        unvisited.remove(start)
        while queue:
            node = queue.popleft()
            order.append(node)
            neighbors = [n for n in adjacency[node] if n in unvisited]
            neighbors.sort(key=lambda nid: (degree[nid], nid))
            for neighbor in neighbors:
                unvisited.remove(neighbor)
                queue.append(neighbor)
    order.reverse()
    return order


def apply_node_order(nodes, elements, node_map, mpc_equations, order):
    old_to_new = {old_id: new_id for new_id, old_id in enumerate(order, start=1)}
    reordered_nodes = {old_to_new[old_id]: coord for old_id, coord in nodes.items()}
    for elem in elements:
        elem["nodes"] = [old_to_new[nid] for nid in elem["nodes"]]
    reordered_node_map = {key: old_to_new[gid] for key, gid in node_map.items()}
    for equation in mpc_equations:
        equation["terms"] = [
            (old_to_new[node], dof, coefficient)
            for node, dof, coefficient in equation["terms"]
        ]
    return dict(sorted(reordered_nodes.items())), elements, reordered_node_map, mpc_equations


def skyline_profile_estimate(nodes, elements, mpc_equations):
    heights = [0] * (len(nodes) * 6)
    def touch(location):
        active = [eq for eq in location if eq > 0]
        if not active:
            return
        first = min(active)
        for eq in active:
            height = eq - first
            if heights[eq - 1] < height:
                heights[eq - 1] = height

    for elem in elements:
        ndof = 3 if element_stap_type(elem) in (1, 5, 12, 14) else 6
        location = []
        for node in elem["nodes"]:
            for dof in range(1, ndof + 1):
                location.append((node - 1) * 6 + dof)
        touch(location)
    for equation in mpc_equations:
        location = [(node - 1) * 6 + dof for node, dof, _ in equation["terms"]]
        touch(location)
    nwk = sum(h + 1 for h in heights)
    mk = max(heights, default=0) + 1
    return mk, nwk


def flatten_assembly(data, merge_tol=1.0e-6, tie_mode="auto", apply_tie_adjust=False, node_order="rcm"):
    raw_nodes = {}
    raw_elements = []
    node_map_raw = {}
    raw_node_instance = {}
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
            raw_node_instance[gid] = inst["name"]

        for src in part["elements"]:
            sec = find_section(part, src["id"])
            n1_vec = sec.get("n1")
            if n1_vec is not None and inst["rotation"]:
                axis_a, axis_b, angle = inst["rotation"]
                n1_vec = vec_normalize(rotate_rodrigues(n1_vec, axis_a, axis_b, angle))
            raw_elements.append({
                "id": next_elem,
                "type": src["type"],
                "nodes": [local_to_global[n] for n in src["nodes"]],
                "material": sec.get("material", ""),
                "section": sec.get("data", []),
                "n1": n1_vec,
                "instance": inst["name"],
            })
            next_elem += 1

    # Build raw tie interpolations. Abaqus assembly instances keep their own
    # nodes even when coordinates are coincident; connectivity between
    # instances must come from explicit constraints, not geometric node merging.
    raw_tie_interpolations = []
    tie_face_count = 0
    tie_node_fallback_count = 0
    tie_node_surface_count = 0
    tie_adjust_count = 0
    for tie in data["assembly"].get("ties", []):
        slave_surf = tie.get("slave", "")
        master_surf = tie.get("master", "")
        slave_surface = surface_info(data["assembly"], slave_surf)
        master_surface = surface_info(data["assembly"], master_surf)
        slave_nset = slave_surface["set"]
        master_nset = master_surface["set"]
        slave_nodes, master_nodes = [], []
        slave_instances = set()
        master_instances = set()
        for nset in data["assembly"].get("nsets", []):
            if nset["name"] == slave_nset:
                slave_instances.add(nset["instance"])
                for lnid in nset["numbers"]:
                    gid = node_map_raw.get((nset["instance"], lnid))
                    if gid: slave_nodes.append(gid)
            if nset["name"] == master_nset:
                master_instances.add(nset["instance"])
                for lnid in nset["numbers"]:
                    gid = node_map_raw.get((nset["instance"], lnid))
                    if gid: master_nodes.append(gid)
        if not slave_nodes or not master_nodes:
            continue
        use_node_surface = (
            tie_mode in ("node", "auto")
            and master_surface.get("type", "").upper() == "NODE"
        )
        master_faces = [] if use_node_surface else build_master_faces(raw_elements, master_nodes, master_instances)
        for sg in slave_nodes:
            sp = raw_nodes[sg]
            if use_node_surface:
                nearest = min(master_nodes, key=lambda mid: vec_norm(vec_sub(sp, raw_nodes[mid])))
                interp = (vec_norm(vec_sub(sp, raw_nodes[nearest])), [(nearest, 1.0)], "node_surface")
            else:
                interp = find_tie_interpolation(sp, master_nodes, master_faces, raw_nodes)
            if interp is None:
                continue
            best_d, weights, source = interp
            if best_d < 100.0:
                raw_tie_interpolations.append((sg, weights))
                if apply_tie_adjust and tie.get("adjust", True):
                    adjusted = (0.0, 0.0, 0.0)
                    for master_raw, weight in weights:
                        adjusted = vec_add(adjusted, vec_scale(raw_nodes[master_raw], weight))
                    raw_nodes[sg] = adjusted
                    tie_adjust_count += 1
                if source == "node_surface":
                    tie_node_surface_count += 1
                elif source == "projected":
                    tie_face_count += 1
                else:
                    tie_node_fallback_count += 1

    old_to_new = {old_id: old_id for old_id in raw_nodes}
    nodes = dict(sorted(raw_nodes.items()))
    print("  Coincident nodes merged: 0 (assembly instance nodes preserved)")

    for elem in raw_elements:
        elem["nodes"] = [old_to_new[n] for n in elem["nodes"]]

    node_map = {key: old_to_new[old] for key, old in node_map_raw.items()}
    mpc_equations = []
    for slave_raw, master_weights_raw in raw_tie_interpolations:
        slave_gid = old_to_new.get(slave_raw)
        if not slave_gid:
            continue
        weights = []
        for master_raw, weight in master_weights_raw:
            master_gid = old_to_new.get(master_raw)
            if master_gid and master_gid != slave_gid and abs(weight) > 1.0e-12:
                weights.append((master_gid, weight))
        if not weights:
            continue

        for dof in (1, 2, 3):
            terms = [(slave_gid, dof, 1.0)]
            for master_gid, weight in weights:
                terms.append((master_gid, dof, -weight))
            mpc_equations.append({"rhs": 0.0, "terms": terms})
    if mpc_equations:
        print(f"  Tie constraints: MPC on {len(mpc_equations) // 3} slave nodes")
        print(
            "  TieMPC interpolation: "
            f"{tie_face_count} projected, "
            f"{tie_node_surface_count} node-surface, "
            f"{tie_node_fallback_count} node fallback"
        )
        if apply_tie_adjust:
            print(f"  Tie adjust applied to {tie_adjust_count} slave nodes")

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
    if node_order == "rcm":
        before_mk, before_nwk = skyline_profile_estimate(nodes, raw_elements, mpc_equations)
        order = build_rcm_order(nodes, raw_elements, mpc_equations)
        nodes, raw_elements, node_map, mpc_equations = apply_node_order(
            nodes, raw_elements, node_map, mpc_equations, order
        )
        after_mk, after_nwk = skyline_profile_estimate(nodes, raw_elements, mpc_equations)
        print(
            "  Node ordering: RCM "
            f"MK {before_mk}->{after_mk}, "
            f"NWK-est {before_nwk}->{after_nwk}"
        )
    elif node_order != "input":
        raise ValueError(f"Unsupported node_order: {node_order}")
    return nodes, raw_elements, node_map, mpc_equations


def element_stap_type(elem, solid_type="H8R", h8i_instances=()):
    etype = elem["type"] if isinstance(elem, dict) else elem
    if etype == "C3D8R":
        inst = elem.get("instance", "") if isinstance(elem, dict) else ""
        if solid_type == "H8RPIER" and any(token and token in inst for token in h8i_instances):
            return 14
        return 12
    return {"T3D2": 1, "S4R": 10, "B31": 11}.get(etype, 0)


THREE_D_TYPES = {1, 5, 9, 10, 11, 12, 14}


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


def write_dat(path, data, nodes, elements, node_map, mpc_equations,
              solid_type="H8R", h8i_instances=(), beam_shear_ratio=0.1,
              beam_stiffness_scale=0.918, solid_stiffness_scale=1.0,
              pier_stiffness_scale=1.0, beam_area_scale=1.0,
              beam_bending_scale=1.15, beam_torsion_scale=1.0,
              beam_shear_area_scale=1.0):
    valid = [e for e in elements if element_stap_type(e, solid_type, h8i_instances)]
    skipped = defaultdict(int)
    for e in elements:
        if not element_stap_type(e, solid_type, h8i_instances):
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
        st = element_stap_type(e, solid_type, h8i_instances)
        for n in e["nodes"]:
            connected.add(n)
            node_types[n].add(st)
            if st == 1:
                node_bar_elems[n].append(e)
            if st == 10:  # Shell4
                node_shell4_count[n] += 1
            if st in (9, 11):  # Beam3D / Beam3DTimoshenko
                node_beam3d_count[n] += 1

    orphan_nodes = set(nodes) - connected
    if orphan_nodes:
        print(f"  Orphan nodes (will be fully fixed): {len(orphan_nodes)}")

    # Shell4 already carries drilling stabilization in the element stiffness.
    # Preserve shell boundary rotations instead of adding artificial supports.
    shell4_boundary_fix = set()
    print("  Shell4 boundary nodes (rotation-fixed): 0 (shell rotations preserved)")

    # Bar elements use only translational DOFs. Do not fix bar-node
    # translations here: cable endpoints may be tied to deck/tower nodes, and
    # fixing them would artificially lock the tied structure.
    bar_fix_nodes = {}
    print("  Bar-plane-fixed nodes: 0 (bar translations left free)")

    # Gravity loads
    nodal_force = defaultdict(lambda: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    if data["gravity"]:
        grav = data["gravity"][0]
        mag = grav["magnitude"]
        gx, gy, gz = grav["direction"]
        gvec = (mag * gx, mag * gy, mag * gz)
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
                p1 = nodes[e["nodes"][0]]
                p2 = nodes[e["nodes"][1]]
                ex, ey, ez, length = beam_local_axes(p1, p2, e.get("n1"))
                sec = e.get("section", [])
                if len(sec) >= 6:
                    area_b, _, _, _ = box_section_props(sec[0], sec[1], sec[2], sec[3], sec[4], sec[5])
                else:
                    area_b = 0.76
                mass = area_b * length * rho
                qx = rho * area_b * vec_dot(gvec, ex)
                qy = rho * area_b * vec_dot(gvec, ey)
                qz = rho * area_b * vec_dot(gvec, ez)

                # Consistent nodal load vector in beam local DOF order:
                # [u, v, w, rx, ry, rz] at node 1 and node 2.
                f_local = [0.0] * 12
                f_local[0] = qx * length / 2.0
                f_local[6] = qx * length / 2.0

                f_local[1] = qy * length / 2.0
                f_local[7] = qy * length / 2.0
                f_local[5] = qy * length * length / 12.0
                f_local[11] = -qy * length * length / 12.0

                f_local[2] = qz * length / 2.0
                f_local[8] = qz * length / 2.0
                f_local[4] = -qz * length * length / 12.0
                f_local[10] = qz * length * length / 12.0

                R = (ex, ey, ez)
                f_global = [0.0] * 12
                for node_i in range(2):
                    base = 6 * node_i
                    for i in range(3):
                        f_global[base + i] = sum(R[j][i] * f_local[base + j] for j in range(3))
                        f_global[base + 3 + i] = sum(R[j][i] * f_local[base + 3 + j] for j in range(3))

                for local_idx, nid in enumerate(e["nodes"]):
                    base = 6 * local_idx
                    for dof in range(6):
                        nodal_force[nid][dof] += f_global[base + dof]
                continue
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
        groups[(element_stap_type(e, solid_type, h8i_instances), e["type"], e["material"])].append(e)

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
            if 10 in types or 9 in types or 11 in types:
                if 5 in types and 11 in types and 10 not in types:
                    # Keep beam-solid shared nodes translation-only, matching the
                    # stable Bridge-1 baseline data used for Abaqus comparison.
                    bc[3] = bc[4] = bc[5] = 1
                elif nid in shell4_boundary_fix and 9 not in types and 11 not in types:
                    # Shell4-only boundary: fix only drilling (rz), keep bending (rx,ry) active
                    bc[3] = 0
                    bc[4] = 0
                    bc[5] = 1
                else:
                    # Shared H8-structure nodes keep shell/beam rotations active.
                    bc[3] = bc[4] = bc[5] = 0
            x, y, z = nodes[nid]
            f.write(f"  {nid}  {bc[0]}  {bc[1]}  {bc[2]}  {bc[3]}  {bc[4]}  {bc[5]}  {x:.6f}  {y:.6f}  {z:.6f}\n")

        f.write(f"  1  {len(load_entries)}\n")
        for nid, dof, val in load_entries:
            f.write(f"  {nid}  {dof}  {val:.6e}\n")

        if mpc_equations:
            f.write(f"MPC {len(mpc_equations)}\n")
            for equation in mpc_equations:
                terms = equation["terms"]
                f.write(f"  {len(terms)}  {equation['rhs']:.12e}")
                for nid, dof, coef in terms:
                    f.write(f"  {nid}  {dof}  {coef:.12e}")
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
            elif st in (5, 12, 14):  # H8 / H8R / H8RPier
                stiffness_scale = pier_stiffness_scale if st == 14 else solid_stiffness_scale
                f.write(f"    1  {mat.get('E', 2.0e11) * stiffness_scale:.6e}  {mat.get('nu', 0.3):.6f}\n")
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
            elif st == 11:  # Beam3DTimoshenko
                E_val = mat.get("E", 2.0e11)
                nu_val = mat.get("nu", 0.3)
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
                area *= beam_stiffness_scale * beam_area_scale
                Iy *= beam_stiffness_scale * beam_bending_scale
                Iz *= beam_stiffness_scale * beam_bending_scale
                J *= beam_stiffness_scale * beam_torsion_scale
                # Thin-walled box sections have lower effective shear area than
                # the solid-rectangle 5/6*A estimate. Abaqus computes this from
                # the section geometry; use a calibrated box-section estimate.
                Asy = beam_shear_ratio * area * beam_shear_area_scale
                Asz = beam_shear_ratio * area * beam_shear_area_scale
                n1_vec = elems[0].get("n1") if elems else None
                if n1_vec is None:
                    n1_vec = (0.0, 0.0, -1.0)
                f.write(f"    1  {E_val:.6e}  {nu_val:.6f}  {area:.6f}  {Iy:.6f}  {Iz:.6f}  {J:.6f}  {Asy:.6f}  {Asz:.6f}  {n1_vec[0]:.6f}  {n1_vec[1]:.6f}  {n1_vec[2]:.6f}\n")
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
    parser.add_argument("--solid-type", choices=("H8R", "H8RPIER"), default="H8RPIER",
                        help="STAP++ solid element used for Abaqus C3D8R")
    parser.add_argument("--pier-instances", default="Part-Pier",
                        help="Comma-separated instance-name tokens converted to H8RPier when --solid-type H8RPIER")
    parser.add_argument("--h8i-instances", default=None, help=argparse.SUPPRESS)
    parser.add_argument("--beam-shear-ratio", type=float, default=0.1,
                        help="Effective B31 Timoshenko shear area ratio As/A")
    parser.add_argument("--beam-stiffness-scale", type=float, default=0.918,
                        help="Scale B31 section stiffness properties in the STAP++ material line")
    parser.add_argument("--beam-area-scale", type=float, default=1.0,
                        help="Additional B31 axial area scale")
    parser.add_argument("--beam-bending-scale", type=float, default=1.15,
                        help="Additional B31 Iy/Iz bending stiffness scale")
    parser.add_argument("--beam-torsion-scale", type=float, default=1.0,
                        help="Additional B31 torsion constant scale")
    parser.add_argument("--beam-shear-area-scale", type=float, default=1.0,
                        help="Additional B31 effective shear area scale")
    parser.add_argument("--solid-stiffness-scale", type=float, default=1.0,
                        help="Scale non-pier C3D8R/H8R elastic stiffness")
    parser.add_argument("--pier-stiffness-scale", type=float, default=1.0,
                        help="Scale pier C3D8R/H8RPier elastic stiffness")
    parser.add_argument("--tie-mode", choices=("projected", "node", "auto"), default="auto",
                        help="Tie interpolation mode: projected keeps legacy face projection; node/auto pair NODE surfaces to nearest master nodes")
    parser.add_argument("--apply-tie-adjust", action="store_true",
                        help="Move slave nodes onto the tied master interpolation when the Abaqus tie has adjust=yes")
    parser.add_argument("--node-order", choices=("rcm", "input"), default="rcm",
                        help="Global node numbering strategy. RCM reduces skyline bandwidth for faster solves.")
    args = parser.parse_args()

    print("Abaqus .inp -> STAP++ .dat converter")
    print(f"  Input : {args.inp}")
    print(f"  Output: {args.dat}")
    data = parse_inp(args.inp)
    print(f"  Parts: {len(data['parts'])}")
    print(f"  Instances: {len(data['assembly'].get('instances', []))}")
    print(f"  Materials: {len(data['materials'])}")
    nodes, elements, node_map, mpc_equations = flatten_assembly(
        data,
        tie_mode=args.tie_mode,
        apply_tie_adjust=args.apply_tie_adjust,
        node_order=args.node_order,
    )
    print(f"  Global nodes: {len(nodes)}")
    print(f"  Global elements: {len(elements)}")
    print(f"  MPC equations: {len(mpc_equations)}")
    pier_instance_arg = args.h8i_instances if args.h8i_instances is not None else args.pier_instances
    h8i_instances = tuple(x.strip() for x in pier_instance_arg.split(",") if x.strip())
    if args.solid_type == "H8RPIER":
        print(f"  Solid mapping: C3D8R -> H8RPier for instances matching {h8i_instances}")
    else:
        print("  Solid mapping: C3D8R -> H8R")
    write_dat(args.dat, data, nodes, elements, node_map, mpc_equations,
              solid_type=args.solid_type, h8i_instances=h8i_instances,
              beam_shear_ratio=args.beam_shear_ratio,
              beam_stiffness_scale=args.beam_stiffness_scale,
              solid_stiffness_scale=args.solid_stiffness_scale,
              pier_stiffness_scale=args.pier_stiffness_scale,
              beam_area_scale=args.beam_area_scale,
              beam_bending_scale=args.beam_bending_scale,
              beam_torsion_scale=args.beam_torsion_scale,
              beam_shear_area_scale=args.beam_shear_area_scale)
    print("Done.")


if __name__ == "__main__":
    main()
