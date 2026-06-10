"""Apply all patches to inp2dat.py to support NDF=6, Beam3D, Shell4."""
with open('D:/STAPpp-master/STAPpp-master/tools/inp2dat/inp2dat.py', 'r') as f:
    content = f.read()

# Patch 1: stap_type mapping
content = content.replace(
    'return {"T3D2": 1, "S4R": 2, "C3D8R": 5, "B31": 1}.get(etype, 0)',
    'return {"T3D2": 1, "S4R": 10, "C3D8R": 5, "B31": 9}.get(etype, 0)'
)

# Patch 2: add THREE_D_TYPES and box_section_props
marker = 'def resolve_fixed_nodes(data, node_map):'
insert = '''THREE_D_TYPES = {1, 5, 9, 10}


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


'''
content = content.replace(marker, insert + marker)

# Patch 3: parse_assembly dict init - add ties and surfaces
content = content.replace(
    'assembly = {"instances": [], "nsets": [], "elsets": []}',
    'assembly = {"instances": [], "nsets": [], "elsets": [], "ties": [], "surfaces": {}}'
)

# Patch 4: add surface/tie parsing in parse_assembly
old_asm = '''        if line.startswith("*Elset"):
            item, i = parse_assembly_set(lines, i, "elset")
            assembly["elsets"].append(item)
            continue
        i += 1
    return assembly, i'''
new_asm = '''        if line.startswith("*Elset"):
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
    return assembly, i'''
content = content.replace(old_asm, new_asm)

# Patch 5: Beam Section n1 parsing
old_sec = '''        if line.startswith("*Solid Section") or line.startswith("*Shell Section") or line.startswith("*Beam Section"):
            elset = header_value(line, "elset")
            section = {"material": header_value(line, "material"), "data": []}
            i += 1
            if i < len(lines) and not lines[i].startswith("*"):
                section["data"] = parse_floats(lines[i])
                i += 1
            part["sections"][elset] = section
            continue'''
new_sec = '''        if line.startswith("*Solid Section") or line.startswith("*Shell Section") or line.startswith("*Beam Section"):
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
            continue'''
content = content.replace(old_sec, new_sec)

# Patch 6: find_section default return
content = content.replace(
    'return {"material": "", "data": []}',
    'return {"material": "", "data": [], "n1": None}'
)

# Patch 7: flatten_assembly - add n1 to element data
old_flatten = '''            sec = find_section(part, src["id"])
            raw_elements.append({
                "id": next_elem,
                "type": src["type"],
                "nodes": [local_to_global[n] for n in src["nodes"]],
                "material": sec.get("material", ""),
                "section": sec.get("data", []),
            })'''
new_flatten = '''            sec = find_section(part, src["id"])
            raw_elements.append({
                "id": next_elem,
                "type": src["type"],
                "nodes": [local_to_global[n] for n in src["nodes"]],
                "material": sec.get("material", ""),
                "section": sec.get("data", []),
                "n1": sec.get("n1"),
            })'''
content = content.replace(old_flatten, new_flatten)

# Patch 8: add tie merging in flatten_assembly
old_ret = '''    node_map = {key: old_to_new[old] for key, old in node_map_raw.items()}
    return nodes, raw_elements, node_map'''
new_ret = '''    node_map = {key: old_to_new[old] for key, old in node_map_raw.items()}

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
        for s, m in zip(slave_nodes, master_nodes):
            slave_to_master[s] = m
    if slave_to_master:
        print(f"  Tie constraints: merging {len(slave_to_master)} slave nodes")
        for elem in raw_elements:
            elem["nodes"] = [slave_to_master.get(n, n) for n in elem["nodes"]]
    return nodes, raw_elements, node_map'''
content = content.replace(old_ret, new_ret)

# Patch 9: resolve_fixed_nodes - NDF=6
old_fix = "    fixed = defaultdict(lambda: [0, 0, 0])\n    for bc in data[\"boundaries\"]:"
new_fix = "    fixed = defaultdict(lambda: [0, 0, 0, 1, 1, 1])\n    for bc in data[\"boundaries\"]:"
content = content.replace(old_fix, new_fix)

# Now replace the entire write_dat function
# Find start and end of write_dat
wstart = content.find('def write_dat(path, data, nodes, elements, node_map):')
wend = content.find('\n\ndef main():')
if wstart < 0 or wend < 0:
    print(f"ERROR: write_dat start={wstart}, end={wend}")
else:
    new_write_dat = '''def write_dat(path, data, nodes, elements, node_map):
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

    # Shell4 boundary: nodes with <=2 connected Shell4 get rotations fixed
    shell4_boundary_fix = set()
    for nid, count in node_shell4_count.items():
        if count <= 2 and 10 in node_types.get(nid, set()):
            shell4_boundary_fix.add(nid)
    print(f"  Shell4 boundary nodes (rotation-fixed): {len(shell4_boundary_fix)}")

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
                mass = polygon_area_xy(e["nodes"], nodes) * 0.2 * rho
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
        f.write("Bridge-1: Cable-stayed bridge (3D conversion)\\n")
        f.write(f"{len(nodes)}  {len(groups)}  1  1\\n")
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
                if nid in shell4_boundary_fix:
                    bc[3] = bc[4] = bc[5] = 1  # fix rotations on boundary
                else:
                    bc[3] = bc[4] = bc[5] = 0  # rotations active for shell/beam
            x, y, z = nodes[nid]
            f.write(f"  {nid}  {bc[0]}  {bc[1]}  {bc[2]}  {bc[3]}  {bc[4]}  {bc[5]}  {x:.6f}  {y:.6f}  {z:.6f}\\n")

        f.write(f"  1  {len(load_entries)}")
        for nid, dof, val in load_entries:
            f.write(f"  {nid}  {dof}  {val:.6e}")
        f.write("\\n")

        for (st, abaqus_type, mat_name), elems in sorted(groups.items()):
            mat = data["materials"].get(mat_name, {"E": 2.0e11, "nu": 0.3})
            f.write(f"  {st}  {len(elems)}  1\\n")
            if st == 1:  # Bar
                area = 0.25
                if abaqus_type == "B31":
                    area = 0.76
                elif abaqus_type == "T3D2" and elems:
                    sec = elems[0].get("section", [])
                    area = sec[0] if sec else 0.25
                f.write(f"    1  {mat.get('E', 2.0e11):.6e}  {area:.6f}\\n")
            elif st == 2:  # Q4
                f.write(f"    1  {mat.get('E', 2.0e11):.6e}  {mat.get('nu', 0.3):.6f}  0.200000\\n")
            elif st == 5:  # H8
                f.write(f"    1  {mat.get('E', 2.0e11):.6e}  {mat.get('nu', 0.3):.6f}\\n")
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
                f.write(f"    1  {E_val:.6e}  {nu_val:.6f}  {area:.6f}  {Iy:.6f}  {Iz:.6f}  {J:.6f}  {n1_vec[0]:.6f}  {n1_vec[1]:.6f}  {n1_vec[2]:.6f}\\n")
            elif st == 10:  # Shell4
                t = 0.2
                f.write(f"    1  {mat.get('E', 2.0e11):.6e}  {mat.get('nu', 0.3):.6f}  {t:.6f}\\n")
            for local_id, elem in enumerate(elems, start=1):
                conn = "  ".join(str(n) for n in elem["nodes"])
                f.write(f"    {local_id}  {conn}  1\\n")
        f.write("  0\\n")
'''
    content = content[:wstart] + new_write_dat + content[wend:]

with open('D:/STAPpp-master/STAPpp-master/tools/inp2dat/inp2dat.py', 'w') as f:
    f.write(content)
print("All patches applied successfully!")
print(f"New file length: {len(content)} chars")
