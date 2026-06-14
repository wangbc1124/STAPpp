#!/usr/bin/env python3
"""Compare Abaqus ODB-exported displacements with STAP++ displacement CSV.

The primary alignment is Abaqus (instance, local node label) to the STAP++
global node id produced by tools/inp2dat/inp2dat.py.  Coordinates are used
only as an explicit fallback diagnostic.
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import importlib.util
import io
import json
import math
import sys
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
INP2DAT_PATH = ROOT / "tools" / "inp2dat" / "inp2dat.py"


def load_inp2dat():
    spec = importlib.util.spec_from_file_location("stappp_inp2dat", str(INP2DAT_PATH))
    if spec is None or spec.loader is None:
        raise RuntimeError("Cannot load inp2dat module from {}".format(INP2DAT_PATH))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def coord_key(x, y, z, ndigits=6):
    return (round(float(x), ndigits), round(float(y), ndigits), round(float(z), ndigits))


def vec_norm3(v):
    return math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])


def read_stappp_displacements(path):
    rows = {}
    with path.open("r", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.DictReader(f)
        required = ("node", "U1", "U2", "U3")
        missing = [name for name in required if name not in (reader.fieldnames or [])]
        if missing:
            raise ValueError("{} missing columns: {}".format(path, ", ".join(missing)))
        for row in reader:
            rows[int(row["node"])] = (
                float(row["U1"]),
                float(row["U2"]),
                float(row["U3"]),
            )
    return rows


def read_abaqus_coords_displacements(path):
    rows = []
    with path.open("r", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.DictReader(f)
        required = ("node", "instance", "x", "y", "z", "U1", "U2", "U3")
        missing = [name for name in required if name not in (reader.fieldnames or [])]
        if missing:
            raise ValueError("{} missing columns: {}".format(path, ", ".join(missing)))
        for row in reader:
            rows.append({
                "node": int(row["node"]),
                "instance": row["instance"],
                "coord": (float(row["x"]), float(row["y"]), float(row["z"])),
                "u": (float(row["U1"]), float(row["U2"]), float(row["U3"])),
            })
    return rows


def build_mapping(inp_path, tie_mode, apply_tie_adjust, node_order):
    inp2dat = load_inp2dat()
    data = inp2dat.parse_inp(str(inp_path))
    # Silence converter diagnostics; the comparison report carries its own stats.
    with contextlib.redirect_stdout(io.StringIO()):
        nodes, _elements, node_map, _mpc = inp2dat.flatten_assembly(
            data,
            tie_mode=tie_mode,
            apply_tie_adjust=apply_tie_adjust,
            node_order=node_order,
        )
    normalized_map = {(inst.upper(), int(local)): int(gid) for (inst, local), gid in node_map.items()}
    return nodes, normalized_map


def build_node_type_categories(inp_path, tie_mode, apply_tie_adjust, node_order):
    inp2dat = load_inp2dat()
    data = inp2dat.parse_inp(str(inp_path))
    with contextlib.redirect_stdout(io.StringIO()):
        _nodes, elements, _node_map, _mpc = inp2dat.flatten_assembly(
            data,
            tie_mode=tie_mode,
            apply_tie_adjust=apply_tie_adjust,
            node_order=node_order,
        )
    node_types = defaultdict(set)
    for elem in elements:
        st = inp2dat.element_stap_type(elem)
        if not st:
            continue
        for node in elem["nodes"]:
            node_types[node].add(st)

    categories = {}
    for node, types in node_types.items():
        has_h8 = bool(types & {5, 12, 14})
        has_shell = 10 in types
        has_beam = bool(types & {9, 11})
        if has_h8 and has_shell and has_beam:
            category = "H8+Shell+Beam"
        elif has_h8 and has_shell:
            category = "H8+Shell"
        elif has_h8 and has_beam:
            category = "H8+Beam"
        elif has_shell and has_beam:
            category = "Shell+Beam"
        elif has_h8:
            category = "H8 only"
        elif has_shell:
            category = "Shell only"
        elif has_beam:
            category = "Beam only"
        else:
            category = "Other"
        categories[node] = category
    return categories


def aggregate_coord_displacements(items):
    groups = defaultdict(list)
    for coord, disp in items:
        groups[coord_key(*coord)].append(disp)
    out = {}
    for key, vals in groups.items():
        n = float(len(vals))
        out[key] = (
            sum(v[0] for v in vals) / n,
            sum(v[1] for v in vals) / n,
            sum(v[2] for v in vals) / n,
        )
    return out


def update_component_stats(stats, comp_names, a, s, meta):
    diffs = {
        "U1": abs(s[0] - a[0]),
        "U2": abs(s[1] - a[1]),
        "U3": abs(s[2] - a[2]),
        "UMag": abs(vec_norm3(s) - vec_norm3(a)),
    }
    refs = {
        "U1": abs(a[0]),
        "U2": abs(a[1]),
        "U3": abs(a[2]),
        "UMag": vec_norm3(a),
    }
    for name in comp_names:
        stats[name]["sq_sum"] += diffs[name] * diffs[name]
        stats[name]["ref_sq_sum"] += refs[name] * refs[name]
        if diffs[name] > stats[name]["max_abs_error"]:
            stats[name]["max_abs_error"] = diffs[name]
            stats[name]["max_abs_error_meta"] = dict(meta)


def finalize_stats(stats, comp_names, count):
    for name in comp_names:
        sq = stats[name].pop("sq_sum")
        ref_sq = stats[name].pop("ref_sq_sum")
        stats[name]["rms_error"] = math.sqrt(sq / float(count)) if count else None
        stats[name]["relative_l2_error"] = math.sqrt(sq / ref_sq) if ref_sq > 0.0 else None


def compare(args):
    nodes, node_map = build_mapping(
        args.inp,
        tie_mode=args.tie_mode,
        apply_tie_adjust=args.apply_tie_adjust,
        node_order=args.node_order,
    )
    node_categories = build_node_type_categories(
        args.inp,
        tie_mode=args.tie_mode,
        apply_tie_adjust=args.apply_tie_adjust,
        node_order=args.node_order,
    )
    stappp = read_stappp_displacements(args.stappp_csv)
    abaqus_rows = read_abaqus_coords_displacements(args.abaqus_csv)

    comp_names = ("U1", "U2", "U3", "UMag")
    stats = {
        name: {
            "max_abs_error": -1.0,
            "max_abs_error_meta": {},
            "sq_sum": 0.0,
            "ref_sq_sum": 0.0,
        }
        for name in comp_names
    }
    worst = []
    by_instance = defaultdict(lambda: {
        name: {
            "max_abs_error": -1.0,
            "max_abs_error_meta": {},
            "sq_sum": 0.0,
            "ref_sq_sum": 0.0,
        }
        for name in comp_names
    })
    by_instance_count = defaultdict(int)
    worst_by_instance = defaultdict(list)
    by_shared_category = defaultdict(lambda: {
        "node_count": 0,
        "stappp_zero_umag_count": 0,
        "large_umag_error_count": 0,
        "umag_abs_error_sum": 0.0,
        "umag_ref_sum": 0.0,
        "umag_error_sq_sum": 0.0,
        "umag_ref_sq_sum": 0.0,
    })
    mapped_count = 0
    missing_mapping = []
    missing_stappp = []
    mapped_abaqus_items = []
    mapped_stappp_items = []

    for row in abaqus_rows:
        map_key = (row["instance"].upper(), row["node"])
        global_node = node_map.get(map_key)
        if global_node is None:
            missing_mapping.append(row)
            continue
        s = stappp.get(global_node)
        if s is None:
            missing_stappp.append((row, global_node))
            continue
        mapped_count += 1
        a = row["u"]
        mapped_abaqus_items.append((row["coord"], a))
        mapped_stappp_items.append((nodes[global_node], s))
        meta = {
            "alignment": "instance_node",
            "instance": row["instance"],
            "abaqus_node": row["node"],
            "stappp_node": global_node,
            "coord": row["coord"],
            "abaqus_u": a,
            "stappp_u": s,
            "umag_abs_error": abs(vec_norm3(s) - vec_norm3(a)),
        }
        update_component_stats(stats, comp_names, a, s, meta)
        inst_key = row["instance"]
        update_component_stats(by_instance[inst_key], comp_names, a, s, meta)
        by_instance_count[inst_key] += 1
        worst_by_instance[inst_key].append((meta["umag_abs_error"], meta))
        worst.append((meta["umag_abs_error"], meta))
        category = node_categories.get(global_node, "Unconnected")
        cat_stats = by_shared_category[category]
        a_umag = vec_norm3(a)
        s_umag = vec_norm3(s)
        err = meta["umag_abs_error"]
        cat_stats["node_count"] += 1
        cat_stats["umag_abs_error_sum"] += err
        cat_stats["umag_ref_sum"] += a_umag
        cat_stats["umag_error_sq_sum"] += err * err
        cat_stats["umag_ref_sq_sum"] += a_umag * a_umag
        if s_umag <= 1.0e-14 and a_umag > 1.0e-10:
            cat_stats["stappp_zero_umag_count"] += 1
        if err > max(1.0e-8, 0.5 * a_umag):
            cat_stats["large_umag_error_count"] += 1

    if mapped_count == 0:
        raise RuntimeError("No nodes matched by instance/local node mapping.")
    finalize_stats(stats, comp_names, mapped_count)
    instance_reports = {}
    for inst_key, inst_stats in by_instance.items():
        finalize_stats(inst_stats, comp_names, by_instance_count[inst_key])
        instance_reports[inst_key] = {
            "node_count": by_instance_count[inst_key],
            "U1": inst_stats["U1"],
            "U2": inst_stats["U2"],
            "U3": inst_stats["U3"],
            "UMag": inst_stats["UMag"],
        }
        worst_by_instance[inst_key].sort(key=lambda item: item[0], reverse=True)
    worst.sort(key=lambda item: item[0], reverse=True)

    coord_abaqus = aggregate_coord_displacements([(r["coord"], r["u"]) for r in abaqus_rows])
    coord_stappp = aggregate_coord_displacements([(nodes[nid], disp) for nid, disp in stappp.items() if nid in nodes])
    coord_common = sorted(set(coord_abaqus) & set(coord_stappp))
    coord_stats = {
        name: {
            "max_abs_error": -1.0,
            "max_abs_error_meta": {},
            "sq_sum": 0.0,
            "ref_sq_sum": 0.0,
        }
        for name in comp_names
    }
    for key in coord_common:
        a = coord_abaqus[key]
        s = coord_stappp[key]
        update_component_stats(
            coord_stats,
            comp_names,
            a,
            s,
            {
                "alignment": "coordinate_fallback",
                "coord": key,
                "abaqus_u": a,
                "stappp_u": s,
                "umag_abs_error": abs(vec_norm3(s) - vec_norm3(a)),
            },
        )
    finalize_stats(coord_stats, comp_names, len(coord_common))

    report = {
        "case": args.case,
        "alignment": {
            "primary": "instance_node",
            "fallback": "coordinate_aggregate",
            "abaqus_rows": len(abaqus_rows),
            "stappp_rows": len(stappp),
            "mapped_common_nodes": mapped_count,
            "missing_abaqus_mapping": len(missing_mapping),
            "missing_stappp_nodes": len(missing_stappp),
            "coordinate_common_points": len(coord_common),
            "coordinate_only_abaqus_points": len(set(coord_abaqus) - set(coord_stappp)),
            "coordinate_only_stappp_points": len(set(coord_stappp) - set(coord_abaqus)),
        },
        "instance_node": {name: stats[name] for name in comp_names},
        "coordinate_fallback": {name: coord_stats[name] for name in comp_names},
        "worst_umag_points": [meta for _err, meta in worst[: args.top]],
        "instance_node_by_instance": instance_reports,
        "worst_umag_points_by_instance": {
            inst: [meta for _err, meta in items[: args.top_per_instance]]
            for inst, items in worst_by_instance.items()
        },
        "shared_node_category_stats": {
            category: {
                "node_count": values["node_count"],
                "stappp_zero_umag_count": values["stappp_zero_umag_count"],
                "large_umag_error_count": values["large_umag_error_count"],
                "mean_umag_abs_error": values["umag_abs_error_sum"] / values["node_count"] if values["node_count"] else None,
                "mean_abaqus_umag": values["umag_ref_sum"] / values["node_count"] if values["node_count"] else None,
                "UMag_relative_l2": math.sqrt(values["umag_error_sq_sum"] / values["umag_ref_sq_sum"]) if values["umag_ref_sq_sum"] > 0.0 else None,
            }
            for category, values in sorted(by_shared_category.items())
        },
        "missing_mapping_samples": [
            {
                "instance": row["instance"],
                "abaqus_node": row["node"],
                "coord": row["coord"],
            }
            for row in missing_mapping[: args.top]
        ],
        "missing_stappp_samples": [
            {
                "instance": row["instance"],
                "abaqus_node": row["node"],
                "stappp_node": global_node,
                "coord": row["coord"],
            }
            for row, global_node in missing_stappp[: args.top]
        ],
    }
    return report


def write_keypoints(path, report):
    rows = []
    for idx, item in enumerate(report["worst_umag_points"], start=1):
        rows.append({
            "rank": idx,
            "label": "worst_umag_{:02d}".format(idx),
            "alignment": item["alignment"],
            "instance": item.get("instance", ""),
            "abaqus_node": item.get("abaqus_node", ""),
            "stappp_node": item.get("stappp_node", ""),
            "x": item["coord"][0],
            "y": item["coord"][1],
            "z": item["coord"][2],
            "abaqus_U1": item["abaqus_u"][0],
            "abaqus_U2": item["abaqus_u"][1],
            "abaqus_U3": item["abaqus_u"][2],
            "stappp_U1": item["stappp_u"][0],
            "stappp_U2": item["stappp_u"][1],
            "stappp_U3": item["stappp_u"][2],
            "UMag_abs_error": item["umag_abs_error"],
        })
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        fieldnames = [
            "rank", "label", "alignment", "instance", "abaqus_node", "stappp_node",
            "x", "y", "z", "abaqus_U1", "abaqus_U2", "abaqus_U3",
            "stappp_U1", "stappp_U2", "stappp_U3", "UMag_abs_error",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_instance_reports(error_path, worst_path, report):
    error_path.parent.mkdir(parents=True, exist_ok=True)
    with error_path.open("w", encoding="utf-8", newline="") as f:
        fieldnames = [
            "instance", "node_count",
            "U1_relative_l2", "U2_relative_l2", "U3_relative_l2", "UMag_relative_l2",
            "U1_max_abs", "U2_max_abs", "U3_max_abs", "UMag_max_abs",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for instance, stats in sorted(report["instance_node_by_instance"].items()):
            writer.writerow({
                "instance": instance,
                "node_count": stats["node_count"],
                "U1_relative_l2": stats["U1"]["relative_l2_error"],
                "U2_relative_l2": stats["U2"]["relative_l2_error"],
                "U3_relative_l2": stats["U3"]["relative_l2_error"],
                "UMag_relative_l2": stats["UMag"]["relative_l2_error"],
                "U1_max_abs": stats["U1"]["max_abs_error"],
                "U2_max_abs": stats["U2"]["max_abs_error"],
                "U3_max_abs": stats["U3"]["max_abs_error"],
                "UMag_max_abs": stats["UMag"]["max_abs_error"],
            })

    with worst_path.open("w", encoding="utf-8", newline="") as f:
        fieldnames = [
            "instance", "rank", "abaqus_node", "stappp_node", "x", "y", "z",
            "abaqus_U1", "abaqus_U2", "abaqus_U3",
            "stappp_U1", "stappp_U2", "stappp_U3", "UMag_abs_error",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for instance, items in sorted(report["worst_umag_points_by_instance"].items()):
            for rank, item in enumerate(items, start=1):
                writer.writerow({
                    "instance": instance,
                    "rank": rank,
                    "abaqus_node": item.get("abaqus_node", ""),
                    "stappp_node": item.get("stappp_node", ""),
                    "x": item["coord"][0],
                    "y": item["coord"][1],
                    "z": item["coord"][2],
                    "abaqus_U1": item["abaqus_u"][0],
                    "abaqus_U2": item["abaqus_u"][1],
                    "abaqus_U3": item["abaqus_u"][2],
                    "stappp_U1": item["stappp_u"][0],
                    "stappp_U2": item["stappp_u"][1],
                    "stappp_U3": item["stappp_u"][2],
                    "UMag_abs_error": item["umag_abs_error"],
                })


def write_shared_category_report(path, report):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        fieldnames = [
            "category", "node_count", "stappp_zero_umag_count",
            "large_umag_error_count", "mean_umag_abs_error",
            "mean_abaqus_umag", "UMag_relative_l2",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for category, stats in sorted(report["shared_node_category_stats"].items()):
            row = {"category": category}
            row.update(stats)
            writer.writerow(row)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", required=True)
    parser.add_argument("--inp", required=True, type=Path)
    parser.add_argument("--abaqus-csv", required=True, type=Path)
    parser.add_argument("--stappp-csv", required=True, type=Path)
    parser.add_argument("--out-json", required=True, type=Path)
    parser.add_argument("--keypoints-csv", required=True, type=Path)
    parser.add_argument("--error-by-instance-csv", type=Path, default=None)
    parser.add_argument("--worst-by-instance-csv", type=Path, default=None)
    parser.add_argument("--shared-category-csv", type=Path, default=None)
    parser.add_argument("--tie-mode", choices=("projected", "node", "auto"), default="auto")
    parser.add_argument("--apply-tie-adjust", action="store_true")
    parser.add_argument("--node-order", choices=("rcm", "input"), default="rcm")
    parser.add_argument("--top", type=int, default=20)
    parser.add_argument("--top-per-instance", type=int, default=5)
    args = parser.parse_args()

    report = compare(args)
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    with args.out_json.open("w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, ensure_ascii=False)
    write_keypoints(args.keypoints_csv, report)
    if args.error_by_instance_csv or args.worst_by_instance_csv:
        error_path = args.error_by_instance_csv or (args.out_json.parent / "error_by_instance.csv")
        worst_path = args.worst_by_instance_csv or (args.out_json.parent / "worst_points_by_instance.csv")
        write_instance_reports(error_path, worst_path, report)
    if args.shared_category_csv:
        write_shared_category_report(args.shared_category_csv, report)
    print(json.dumps({
        "case": report["case"],
        "mapped_common_nodes": report["alignment"]["mapped_common_nodes"],
        "missing_abaqus_mapping": report["alignment"]["missing_abaqus_mapping"],
        "missing_stappp_nodes": report["alignment"]["missing_stappp_nodes"],
        "UMag_relative_l2": report["instance_node"]["UMag"]["relative_l2_error"],
        "U3_relative_l2": report["instance_node"]["U3"]["relative_l2_error"],
        "out_json": str(args.out_json),
        "keypoints_csv": str(args.keypoints_csv),
    }, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    raise SystemExit(main())
