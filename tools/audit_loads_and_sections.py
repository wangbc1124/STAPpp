#!/usr/bin/env python3
"""Audit STAP++ load and section equivalence against the Abaqus input model."""

from __future__ import annotations

import argparse
import contextlib
import csv
import importlib.util
import io
import json
import math
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
INP2DAT_PATH = ROOT / "tools" / "inp2dat" / "inp2dat.py"


def load_inp2dat():
    spec = importlib.util.spec_from_file_location("stappp_inp2dat", str(INP2DAT_PATH))
    if spec is None or spec.loader is None:
        raise RuntimeError("Cannot load inp2dat from {}".format(INP2DAT_PATH))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def add_vec(a, b):
    for i in range(3):
        a[i] += b[i]


def moment_from_force(point, force):
    x, y, z = point
    fx, fy, fz = force
    return [y * fz - z * fy, z * fx - x * fz, x * fy - y * fx]


def risk_for_value(value, expected=None, scaled=False):
    if value is None:
        return "missing"
    if scaled:
        return "scaled"
    if expected is not None and abs(value - expected) > max(1.0e-9, abs(expected) * 1.0e-8):
        return "suspicious"
    return "exact"


def parse_dat_loads(dat_path):
    tokens = []
    with dat_path.open("r", encoding="utf-8", errors="replace") as f:
        _title = f.readline()
        for line in f:
            tokens.extend(line.split())
    pos = 0
    numnp = int(tokens[pos]); pos += 1
    _numeg = int(tokens[pos]); pos += 1
    nlcase = int(tokens[pos]); pos += 1
    pos += 1

    # Current bridge files are NDF=6. Keep a fallback for older 3-DOF data.
    for ndf in (6, 3):
        trial = pos
        nodes = {}
        try:
            for _ in range(numnp):
                nid = int(tokens[trial]); trial += 1
                trial += ndf
                nodes[nid] = (float(tokens[trial]), float(tokens[trial + 1]), float(tokens[trial + 2]))
                trial += 3
            loads = []
            for _ in range(nlcase):
                lcase = int(tokens[trial]); trial += 1
                nload = int(tokens[trial]); trial += 1
                for _ in range(nload):
                    loads.append((lcase, int(tokens[trial]), int(tokens[trial + 1]), float(tokens[trial + 2])))
                    trial += 3
            return nodes, loads
        except (ValueError, IndexError):
            continue
    raise RuntimeError("Cannot parse nodal loads from {}".format(dat_path))


def summarize_force(rows):
    total = {
        "mass": 0.0,
        "gravity_force_x": 0.0,
        "gravity_force_y": 0.0,
        "gravity_force_z": 0.0,
        "gravity_moment_x": 0.0,
        "gravity_moment_y": 0.0,
        "gravity_moment_z": 0.0,
        "element_count": 0,
    }
    for row in rows:
        for key in total:
            total[key] += row.get(key, 0.0)
    return total


def compute_element_load_rows(inp2dat, elements, nodes, data):
    if not data["gravity"]:
        return []
    grav = data["gravity"][0]
    mag = grav["magnitude"]
    gx, gy, gz = grav["direction"]
    gvec = (mag * gx, mag * gy, mag * gz)
    rows = []
    for elem in elements:
        mat = data["materials"].get(elem["material"], {"density": 2320.0})
        rho = mat.get("density", 2320.0)
        mass = 0.0
        forces = []
        if elem["type"] == "S4R":
            sec = elem.get("section", [])
            thickness = sec[0] if sec else 0.2
            mass = inp2dat.polygon_area_xy(elem["nodes"], nodes) * thickness * rho
        elif elem["type"] == "T3D2":
            length = inp2dat.vec_norm(inp2dat.vec_sub(nodes[elem["nodes"][1]], nodes[elem["nodes"][0]]))
            area = elem["section"][0] if elem["section"] else 0.25
            mass = area * length * rho
        elif elem["type"] == "B31":
            p1 = nodes[elem["nodes"][0]]
            p2 = nodes[elem["nodes"][1]]
            ex, ey, ez, length = inp2dat.beam_local_axes(p1, p2, elem.get("n1"))
            sec = elem.get("section", [])
            if len(sec) >= 6:
                area_b, _, _, _ = inp2dat.box_section_props(sec[0], sec[1], sec[2], sec[3], sec[4], sec[5])
            else:
                area_b = 0.76
            mass = area_b * length * rho
            qx = rho * area_b * inp2dat.vec_dot(gvec, ex)
            qy = rho * area_b * inp2dat.vec_dot(gvec, ey)
            qz = rho * area_b * inp2dat.vec_dot(gvec, ez)
            f_local = [0.0] * 12
            f_local[0] = f_local[6] = qx * length / 2.0
            f_local[1] = f_local[7] = qy * length / 2.0
            f_local[5] = qy * length * length / 12.0
            f_local[11] = -qy * length * length / 12.0
            f_local[2] = f_local[8] = qz * length / 2.0
            f_local[4] = -qz * length * length / 12.0
            f_local[10] = qz * length * length / 12.0
            axes = (ex, ey, ez)
            for node_i, nid in enumerate(elem["nodes"]):
                base = 6 * node_i
                force = [
                    sum(axes[j][i] * f_local[base + j] for j in range(3))
                    for i in range(3)
                ]
                forces.append((nid, force))
        elif elem["type"] == "C3D8R":
            mass = inp2dat.hex_volume(elem["nodes"], nodes) * rho

        if not forces:
            per_node = mass * mag / max(1, len(elem["nodes"]))
            force = [per_node * gx, per_node * gy, per_node * gz]
            forces = [(nid, force) for nid in elem["nodes"]]

        force_total = [0.0, 0.0, 0.0]
        moment_total = [0.0, 0.0, 0.0]
        for nid, force in forces:
            add_vec(force_total, force)
            add_vec(moment_total, moment_from_force(nodes[nid], force))
        rows.append({
            "instance": elem.get("instance", ""),
            "abaqus_type": elem["type"],
            "material": elem.get("material", ""),
            "mass": mass,
            "gravity_force_x": force_total[0],
            "gravity_force_y": force_total[1],
            "gravity_force_z": force_total[2],
            "gravity_moment_x": moment_total[0],
            "gravity_moment_y": moment_total[1],
            "gravity_moment_z": moment_total[2],
            "element_count": 1,
        })
    return rows


def group_rows(rows, keys):
    grouped = defaultdict(list)
    for row in rows:
        grouped[tuple(row[k] for k in keys)].append(row)
    out = []
    for key, items in sorted(grouped.items()):
        summary = summarize_force(items)
        for name, value in zip(keys, key):
            summary[name] = value
        out.append(summary)
    return out


def write_csv(path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def read_reactions(path):
    totals = defaultdict(lambda: [0.0, 0.0, 0.0])
    if not path or not path.exists():
        return None
    with path.open("r", encoding="utf-8", errors="replace", newline="") as f:
        for row in csv.DictReader(f):
            inst = row["instance"]
            totals[inst][0] += float(row["RF1"])
            totals[inst][1] += float(row["RF2"])
            totals[inst][2] += float(row["RF3"])
            totals["__TOTAL__"][0] += float(row["RF1"])
            totals["__TOTAL__"][1] += float(row["RF2"])
            totals["__TOTAL__"][2] += float(row["RF3"])
    return {k: {"RF1": v[0], "RF2": v[1], "RF3": v[2]} for k, v in totals.items()}


def section_rows(inp2dat, elements, data, solid_type, h8i_instances,
                 beam_shear_ratio, beam_stiffness_scale, beam_area_scale,
                 beam_bending_scale, beam_torsion_scale, beam_shear_area_scale,
                 solid_stiffness_scale, pier_stiffness_scale):
    rows = []
    seen = set()
    for elem in elements:
        st = inp2dat.element_stap_type(elem, solid_type, h8i_instances)
        if not st:
            continue
        key = (elem.get("instance", ""), elem["type"], elem.get("material", ""), tuple(elem.get("section", [])), st)
        if key in seen:
            continue
        seen.add(key)
        mat = data["materials"].get(elem["material"], {})
        sec = elem.get("section", [])
        row = {
            "instance": elem.get("instance", ""),
            "abaqus_type": elem["type"],
            "stappp_type": st,
            "material": elem.get("material", ""),
            "E_inp": mat.get("E"),
            "nu_inp": mat.get("nu"),
            "density": mat.get("density"),
            "A_inp": None,
            "Iy_inp": None,
            "Iz_inp": None,
            "J_inp": None,
            "Asy_stappp": None,
            "Asz_stappp": None,
            "E_stappp": mat.get("E"),
            "risk": "exact",
            "note": "",
        }
        if elem["type"] == "B31":
            if len(sec) >= 6:
                area, iy, iz, j = inp2dat.box_section_props(sec[0], sec[1], sec[2], sec[3], sec[4], sec[5])
                row.update({
                    "A_inp": area,
                    "Iy_inp": iy,
                    "Iz_inp": iz,
                    "J_inp": j,
                    "A_stappp": area * beam_stiffness_scale * beam_area_scale,
                    "Iy_stappp": iy * beam_stiffness_scale * beam_bending_scale,
                    "Iz_stappp": iz * beam_stiffness_scale * beam_bending_scale,
                    "J_stappp": j * beam_stiffness_scale * beam_torsion_scale,
                })
                row["Asy_stappp"] = beam_shear_ratio * row["A_stappp"] * beam_shear_area_scale
                row["Asz_stappp"] = beam_shear_ratio * row["A_stappp"] * beam_shear_area_scale
                row["risk"] = "scaled" if (
                    beam_stiffness_scale != 1.0 or beam_bending_scale != 1.0
                    or beam_shear_ratio != 1.0
                ) else "exact"
                row["note"] = "B31 box section converted to Beam3DTimoshenko"
            else:
                row["risk"] = "defaulted"
                row["note"] = "B31 section data missing; converter defaults are used"
        elif elem["type"] == "S4R":
            row["thickness_inp"] = sec[0] if sec else None
            row["thickness_stappp"] = sec[0] if sec else 0.2
            row["risk"] = risk_for_value(row["thickness_inp"])
            row["note"] = "S4R mapped to Shell4"
        elif elem["type"] == "C3D8R":
            scale = pier_stiffness_scale if st == 14 else solid_stiffness_scale
            row["E_stappp"] = mat.get("E") * scale if mat.get("E") is not None else None
            row["risk"] = risk_for_value(row["E_stappp"], mat.get("E"), scaled=(scale != 1.0 or st == 14))
            row["note"] = "C3D8R mapped to H8RPier" if st == 14 else "C3D8R mapped to H8R"
        elif elem["type"] == "T3D2":
            row["A_inp"] = sec[0] if sec else None
            row["A_stappp"] = sec[0] if sec else 0.25
            row["risk"] = risk_for_value(row["A_inp"])
            row["note"] = "T3D2 mapped to Bar"
        rows.append(row)
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", required=True)
    parser.add_argument("--inp", required=True, type=Path)
    parser.add_argument("--dat", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--reactions-csv", type=Path, default=None)
    parser.add_argument("--solid-type", choices=("H8R", "H8RPIER"), default="H8RPIER")
    parser.add_argument("--pier-instances", default="Part-Pier")
    parser.add_argument("--tie-mode", choices=("projected", "node", "auto"), default="auto")
    parser.add_argument("--apply-tie-adjust", action="store_true")
    parser.add_argument("--node-order", choices=("rcm", "input"), default="rcm")
    parser.add_argument("--beam-shear-ratio", type=float, default=0.1)
    parser.add_argument("--beam-stiffness-scale", type=float, default=0.918)
    parser.add_argument("--beam-area-scale", type=float, default=1.0)
    parser.add_argument("--beam-bending-scale", type=float, default=1.15)
    parser.add_argument("--beam-torsion-scale", type=float, default=1.0)
    parser.add_argument("--beam-shear-area-scale", type=float, default=1.0)
    parser.add_argument("--solid-stiffness-scale", type=float, default=1.0)
    parser.add_argument("--pier-stiffness-scale", type=float, default=1.0)
    args = parser.parse_args()

    inp2dat = load_inp2dat()
    data = inp2dat.parse_inp(str(args.inp))
    with contextlib.redirect_stdout(io.StringIO()):
        nodes, elements, _node_map, _mpc = inp2dat.flatten_assembly(
            data,
            tie_mode=args.tie_mode,
            apply_tie_adjust=args.apply_tie_adjust,
            node_order=args.node_order,
        )

    load_rows = compute_element_load_rows(inp2dat, elements, nodes, data)
    by_instance = group_rows(load_rows, ["instance"])
    by_type = group_rows(load_rows, ["abaqus_type", "material"])
    total = summarize_force(load_rows)

    dat_nodes, dat_loads = parse_dat_loads(args.dat)
    dat_total_force = [0.0, 0.0, 0.0]
    dat_total_moment = [0.0, 0.0, 0.0]
    for _lcase, nid, dof, value in dat_loads:
        if 1 <= dof <= 3:
            dat_total_force[dof - 1] += value
            force = [0.0, 0.0, 0.0]
            force[dof - 1] = value
            add_vec(dat_total_moment, moment_from_force(dat_nodes[nid], force))

    h8i_instances = tuple(x.strip() for x in args.pier_instances.split(",") if x.strip())
    sections = section_rows(
        inp2dat, elements, data, args.solid_type, h8i_instances,
        args.beam_shear_ratio, args.beam_stiffness_scale, args.beam_area_scale,
        args.beam_bending_scale, args.beam_torsion_scale, args.beam_shear_area_scale,
        args.solid_stiffness_scale, args.pier_stiffness_scale,
    )

    args.out_dir.mkdir(parents=True, exist_ok=True)
    load_fields = [
        "instance", "abaqus_type", "material", "mass",
        "gravity_force_x", "gravity_force_y", "gravity_force_z",
        "gravity_moment_x", "gravity_moment_y", "gravity_moment_z", "element_count",
    ]
    write_csv(args.out_dir / "load_summary_by_instance.csv", by_instance, [f for f in load_fields if f != "abaqus_type" and f != "material"])
    write_csv(args.out_dir / "load_summary_by_element_type.csv", by_type, [f for f in load_fields if f != "instance"])
    section_fields = sorted({key for row in sections for key in row.keys()})
    write_csv(args.out_dir / "section_stiffness_audit.csv", sections, section_fields)

    reactions = read_reactions(args.reactions_csv)
    balance = {
        "case": args.case,
        "converter_total": total,
        "dat_total_force": {
            "Fx": dat_total_force[0],
            "Fy": dat_total_force[1],
            "Fz": dat_total_force[2],
            "Mx": dat_total_moment[0],
            "My": dat_total_moment[1],
            "Mz": dat_total_moment[2],
            "load_entries": len(dat_loads),
        },
        "abaqus_reactions": reactions,
        "notes": [],
    }
    if reactions and "__TOTAL__" in reactions:
        rf = reactions["__TOTAL__"]
        balance["force_plus_reaction"] = {
            "x": dat_total_force[0] + rf["RF1"],
            "y": dat_total_force[1] + rf["RF2"],
            "z": dat_total_force[2] + rf["RF3"],
        }
    else:
        balance["notes"].append("Abaqus reaction CSV missing; reaction balance was not evaluated.")
    with (args.out_dir / "load_balance.json").open("w", encoding="utf-8") as f:
        json.dump(balance, f, indent=2, ensure_ascii=False)

    print(json.dumps({
        "case": args.case,
        "out_dir": str(args.out_dir),
        "dat_total_force": balance["dat_total_force"],
        "abaqus_total_reaction": reactions.get("__TOTAL__") if reactions else None,
        "force_plus_reaction": balance.get("force_plus_reaction"),
        "section_risk_counts": dict((risk, sum(1 for row in sections if row.get("risk") == risk))
                                    for risk in sorted({row.get("risk") for row in sections})),
    }, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    raise SystemExit(main())
