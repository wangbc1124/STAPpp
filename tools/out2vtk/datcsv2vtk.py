#!/usr/bin/env python3
"""
Merge a STAP++ .dat mesh and displacement CSV into a legacy VTK file.

Usage:
  python datcsv2vtk.py Bridge-3.dat Bridge-3.displacements.csv Bridge-3.vtk
"""

import csv
import math
import os
import sys
import time

from dat2vtk import ELEM_DEF


def read_displacements(path):
    displacements = {}
    with open(path, "r", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.DictReader(f)
        required = ("node", "U1", "U2", "U3", "R1", "R2", "R3")
        if reader.fieldnames is None:
            raise ValueError("CSV is missing header row")
        missing = [name for name in required if name not in reader.fieldnames]
        if missing:
            raise ValueError(
                "CSV is missing required columns: {0}".format(", ".join(missing))
            )
        for row in reader:
            node_id = int(row["node"])
            displacements[node_id] = (
                float(row["U1"]),
                float(row["U2"]),
                float(row["U3"]),
                float(row["R1"]),
                float(row["R2"]),
                float(row["R3"]),
            )
    return displacements


def parse_dat_with_mpc(path):
    title = ""
    nodes = []
    elements = []
    elem_types = set()

    with open(path, "r", encoding="utf-8", errors="replace") as f:
        title = f.readline().strip()
        token_lines = []
        for raw in f:
            token_lines.extend(raw.split())

    if len(token_lines) < 4:
        raise ValueError("File too short: missing control line")

    tokens = [tok for tok in token_lines if tok]
    pos = 0
    numnp = int(tokens[pos]); pos += 1
    numeg = int(tokens[pos]); pos += 1
    nlcase = int(tokens[pos]); pos += 1
    pos += 1  # MODEX

    last_error = None
    for ndf in (3, 6):
        try:
            trial_pos = pos
            trial_nodes = []
            trial_elements = []
            trial_elem_types = set()

            for _ in range(numnp):
                trial_pos += 1
                trial_pos += ndf
                x = float(tokens[trial_pos]); trial_pos += 1
                y = float(tokens[trial_pos]); trial_pos += 1
                z = float(tokens[trial_pos]); trial_pos += 1
                trial_nodes.append((x, y, z))

            for _ in range(nlcase):
                trial_pos += 1  # load case id
                nload = int(tokens[trial_pos]); trial_pos += 1
                trial_pos += nload * 3

            if trial_pos < len(tokens) and tokens[trial_pos].upper() == "MPC":
                nmpc = int(tokens[trial_pos + 1])
                trial_pos += 2
                for _ in range(nmpc):
                    nterms = int(tokens[trial_pos]); trial_pos += 1
                    trial_pos += 1  # rhs
                    trial_pos += nterms * 3

            group_index = 0
            while trial_pos < len(tokens):
                if tokens[trial_pos] == "0":
                    trial_pos += 1
                    break

                if trial_pos + 2 >= len(tokens):
                    raise ValueError("Truncated element group header")

                group_index += 1
                etype = int(tokens[trial_pos]); trial_pos += 1
                nume = int(tokens[trial_pos]); trial_pos += 1
                nmat = int(tokens[trial_pos]); trial_pos += 1

                if etype not in ELEM_DEF:
                    raise ValueError(
                        "Unsupported element type {0} in group {1}".format(
                            etype, group_index
                        )
                    )

                trial_elem_types.add(etype)
                _, nodes_per_elem, mat_params = ELEM_DEF[etype]
                trial_pos += nmat * mat_params

                for _ in range(nume):
                    if trial_pos + 1 + nodes_per_elem >= len(tokens):
                        raise ValueError("Truncated connectivity data")
                    trial_pos += 1  # element id
                    conn = [int(tokens[trial_pos + i]) - 1 for i in range(nodes_per_elem)]
                    trial_pos += nodes_per_elem
                    trial_pos += 1  # material id
                    trial_elements.append({
                        "type": etype,
                        "connectivity": conn,
                    })

            if group_index != numeg:
                raise ValueError(
                    "Parsed group count {0} does not match control NUMEG {1}".format(
                        group_index, numeg
                    )
                )

            return title, trial_nodes, trial_elements, sorted(trial_elem_types), ndf
        except (ValueError, OverflowError, IndexError) as exc:
            last_error = exc

    raise ValueError(
        "Failed to detect node DOF count from .dat: {0}".format(last_error)
    )


def write_vtk(path, title, nodes, elements, displacements_by_node):
    num_nodes = len(nodes)
    num_cells = len(elements)

    if num_nodes == 0:
        raise ValueError("No nodes found in .dat")
    if num_cells == 0:
        raise ValueError("No elements found in .dat")

    missing_nodes = [nid for nid in range(1, num_nodes + 1) if nid not in displacements_by_node]
    extra_nodes = [nid for nid in displacements_by_node if nid < 1 or nid > num_nodes]
    if missing_nodes:
        preview = ", ".join(str(nid) for nid in missing_nodes[:10])
        raise ValueError(
            "Displacement CSV is missing {0} nodes; first missing ids: {1}".format(
                len(missing_nodes), preview
            )
        )
    if extra_nodes:
        preview = ", ".join(str(nid) for nid in sorted(extra_nodes)[:10])
        raise ValueError(
            "Displacement CSV contains node ids outside mesh range; first ids: {0}".format(
                preview
            )
        )

    out_dir = os.path.dirname(path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(path, "w", encoding="utf-8") as f:
        f.write("# vtk DataFile Version 2.0\n")
        f.write("{0}\n".format(title))
        f.write("ASCII\n")
        f.write("DATASET UNSTRUCTURED_GRID\n")

        f.write("\nPOINTS {0} float\n".format(num_nodes))
        for x, y, z in nodes:
            f.write("  {0:14.6e} {1:14.6e} {2:14.6e}\n".format(x, y, z))

        total_size = sum(1 + len(elem["connectivity"]) for elem in elements)
        f.write("\nCELLS {0} {1}\n".format(num_cells, total_size))
        for elem in elements:
            conn = elem["connectivity"]
            f.write("  {0} {1}\n".format(len(conn), " ".join(str(nid) for nid in conn)))

        f.write("\nCELL_TYPES {0}\n".format(num_cells))
        for elem in elements:
            vtk_type = ELEM_DEF[elem["type"]][0]
            f.write("  {0}\n".format(vtk_type))

        f.write("\nPOINT_DATA {0}\n".format(num_nodes))
        f.write("VECTORS Displacement float\n")
        for node_id in range(1, num_nodes + 1):
            u1, u2, u3, _, _, _ = displacements_by_node[node_id]
            f.write("  {0:14.6e} {1:14.6e} {2:14.6e}\n".format(u1, u2, u3))

        f.write("VECTORS Rotation float\n")
        for node_id in range(1, num_nodes + 1):
            _, _, _, r1, r2, r3 = displacements_by_node[node_id]
            f.write("  {0:14.6e} {1:14.6e} {2:14.6e}\n".format(r1, r2, r3))

        f.write("SCALARS DisplacementMagnitude float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for node_id in range(1, num_nodes + 1):
            u1, u2, u3, _, _, _ = displacements_by_node[node_id]
            mag = math.sqrt(u1 * u1 + u2 * u2 + u3 * u3)
            f.write("  {0:14.6e}\n".format(mag))

        f.write("SCALARS RotationMagnitude float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for node_id in range(1, num_nodes + 1):
            _, _, _, r1, r2, r3 = displacements_by_node[node_id]
            mag = math.sqrt(r1 * r1 + r2 * r2 + r3 * r3)
            f.write("  {0:14.6e}\n".format(mag))

        f.write("\nCELL_DATA {0}\n".format(num_cells))
        f.write("SCALARS ElementType int 1\n")
        f.write("LOOKUP_TABLE default\n")
        for elem in elements:
            f.write("  {0}\n".format(elem["type"]))


def main(argv):
    if len(argv) < 3:
        print(__doc__)
        return 1

    dat_path = argv[1]
    csv_path = argv[2]
    vtk_path = argv[3] if len(argv) >= 4 else os.path.splitext(csv_path)[0] + ".vtk"

    if not os.path.exists(dat_path):
        raise FileNotFoundError("Missing .dat: {0}".format(dat_path))
    if not os.path.exists(csv_path):
        raise FileNotFoundError("Missing displacement CSV: {0}".format(csv_path))

    print("Converting mesh+csv:")
    print("  dat: {0}".format(dat_path))
    print("  csv: {0}".format(csv_path))
    print("  vtk: {0}".format(vtk_path))

    t0 = time.perf_counter()
    title, nodes, elements, elem_types, ndf = parse_dat_with_mpc(dat_path)
    t1 = time.perf_counter()
    displacements = read_displacements(csv_path)
    t2 = time.perf_counter()

    print("  nodes: {0} (NDF={1})".format(len(nodes), ndf))
    print("  elements: {0}".format(len(elements)))
    print("  elem types: {0}".format(elem_types))
    print("  csv rows: {0}".format(len(displacements)))
    print("  parse dat time: {0:.3f} s".format(t1 - t0))
    print("  read csv time: {0:.3f} s".format(t2 - t1))

    write_vtk(vtk_path, title, nodes, elements, displacements)
    t3 = time.perf_counter()
    print("  write vtk time: {0:.3f} s".format(t3 - t2))
    print("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
