#!/usr/bin/env python3
"""Compare Bridge-1 STAP++ displacements with Abaqus CSV data by coordinates."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


def coord_key(x, y, z, ndigits):
    return (round(float(x), ndigits), round(float(y), ndigits), round(float(z), ndigits))


def norm3(v):
    return math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])


def component(instance):
    name = instance.upper()
    if "SUPPORTBEAM" in name:
        return "SupportBeam/B31"
    if "FLOOR" in name:
        return "Floor/S4R"
    if "PIER" in name:
        return "Pier/C3D8R"
    if "RIVERBANK" in name:
        return "RiverBank/C3D8R"
    if "CABLE" in name:
        return "Cable/T3D2"
    return "Other"


def read_dat_nodes(path):
    with Path(path).open("r", encoding="utf-8", errors="replace") as f:
        lines = f.readlines()
    if len(lines) < 2:
        raise ValueError("{} is too short to be a STAP++ dat file".format(path))
    numnp = int(lines[1].split()[0])
    nodes = {}
    for line in lines[2:2 + numnp]:
        parts = line.split()
        if len(parts) < 10:
            continue
        nodes[int(parts[0])] = (float(parts[7]), float(parts[8]), float(parts[9]))
    if len(nodes) != numnp:
        raise ValueError("Read {} nodes from {}, expected {}".format(len(nodes), path, numnp))
    return nodes


def read_stappp_csv(path):
    rows = {}
    with Path(path).open("r", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows[int(row["node"])] = (float(row["U1"]), float(row["U2"]), float(row["U3"]))
    return rows


def read_abaqus_nodes(path, ndigits):
    coords = {}
    with Path(path).open("r", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            row = {str(k).strip(): v for k, v in row.items()}
            instance = (row.get("Instance_Name") or row.get("instance") or row.get("Instance") or "").strip()
            node = int(row.get("Node_ID") or row.get("node") or row.get("NodeID"))
            x = row.get("X") or row.get("x")
            y = row.get("Y") or row.get("y")
            z = row.get("Z") or row.get("z")
            coords[(instance.upper(), node)] = coord_key(x, y, z, ndigits)
    return coords


def read_abaqus_displacements(path):
    disp = {}
    with Path(path).open("r", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            instance = (row.get("Instance") or row.get("instance") or "").strip()
            node = int(row.get("NodeID") or row.get("node"))
            u1 = row.get("DX") or row.get("U1")
            u2 = row.get("DY") or row.get("U2")
            u3 = row.get("DZ") or row.get("U3")
            disp[(instance.upper(), node)] = (float(u1), float(u2), float(u3))
    return disp


def average_vectors(items):
    groups = defaultdict(list)
    for key, vec in items:
        groups[key].append(vec)
    out = {}
    for key, values in groups.items():
        n = float(len(values))
        out[key] = (
            sum(v[0] for v in values) / n,
            sum(v[1] for v in values) / n,
            sum(v[2] for v in values) / n,
        )
    return out


def percentile(values, fraction):
    if not values:
        return None
    values = sorted(values)
    index = min(len(values) - 1, max(0, int(round((len(values) - 1) * fraction))))
    return values[index]


def compare(args):
    dat_nodes = read_dat_nodes(args.dat)
    stap = read_stappp_csv(args.stappp_csv)
    abaqus_nodes = read_abaqus_nodes(args.abaqus_nodes, args.coord_digits)
    abaqus_disp = read_abaqus_displacements(args.abaqus_displacements)

    stap_by_coord = average_vectors(
        (coord_key(*dat_nodes[node], args.coord_digits), disp)
        for node, disp in stap.items()
        if node in dat_nodes
    )

    rows = []
    missing_coord = 0
    for key, a_disp in abaqus_disp.items():
        coord = abaqus_nodes.get(key)
        if coord is None:
            missing_coord += 1
            continue
        s_disp = stap_by_coord.get(coord)
        if s_disp is None:
            continue
        diff = (s_disp[0] - a_disp[0], s_disp[1] - a_disp[1], s_disp[2] - a_disp[2])
        a_mag = norm3(a_disp)
        s_mag = norm3(s_disp)
        abs_vec = norm3(diff)
        rel_vec = abs_vec / a_mag if a_mag > args.min_reference else None
        rel_mag = abs(s_mag - a_mag) / a_mag if a_mag > args.min_reference else None
        rows.append({
            "instance": key[0],
            "node": key[1],
            "component": component(key[0]),
            "coord": coord,
            "aba": a_disp,
            "stap": s_disp,
            "abs_vec": abs_vec,
            "rel_vec": rel_vec,
            "rel_mag": rel_mag,
        })

    if not rows:
        raise RuntimeError("No common coordinate-displacement rows were found.")

    rel_rows = [r for r in rows if r["rel_vec"] is not None]
    thresholds = {}
    for threshold in args.thresholds:
        thresholds[str(threshold)] = 100.0 * sum(
            1 for row in rel_rows if row["rel_vec"] < threshold
        ) / float(len(rel_rows) or 1)

    by_component = defaultdict(list)
    for row in rel_rows:
        by_component[row["component"]].append(row)

    summary_rows = []
    for comp, items in sorted(by_component.items()):
        values = [r["rel_vec"] for r in items if r["rel_vec"] is not None]
        summary_rows.append({
            "component": comp,
            "count": len(items),
            "mean_rel_percent": 100.0 * sum(values) / len(values) if values else None,
            "p50_rel_percent": 100.0 * percentile(values, 0.50) if values else None,
            "p90_rel_percent": 100.0 * percentile(values, 0.90) if values else None,
            "under_10_percent": 100.0 * sum(1 for v in values if v < 0.10) / float(len(values) or 1),
        })

    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.out_csv.open("w", encoding="utf-8", newline="") as f:
        fieldnames = [
            "instance", "node", "component", "x", "y", "z",
            "aba_U1", "aba_U2", "aba_U3", "stap_U1", "stap_U2", "stap_U3",
            "abs_vec", "rel_vec_percent", "rel_mag_percent",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in sorted(rows, key=lambda r: r["abs_vec"], reverse=True):
            x, y, z = row["coord"]
            writer.writerow({
                "instance": row["instance"],
                "node": row["node"],
                "component": row["component"],
                "x": x,
                "y": y,
                "z": z,
                "aba_U1": row["aba"][0],
                "aba_U2": row["aba"][1],
                "aba_U3": row["aba"][2],
                "stap_U1": row["stap"][0],
                "stap_U2": row["stap"][1],
                "stap_U3": row["stap"][2],
                "abs_vec": row["abs_vec"],
                "rel_vec_percent": None if row["rel_vec"] is None else 100.0 * row["rel_vec"],
                "rel_mag_percent": None if row["rel_mag"] is None else 100.0 * row["rel_mag"],
            })

    if args.component_csv:
        args.component_csv.parent.mkdir(parents=True, exist_ok=True)
        with args.component_csv.open("w", encoding="utf-8", newline="") as f:
            fieldnames = ["component", "count", "mean_rel_percent", "p50_rel_percent", "p90_rel_percent", "under_10_percent"]
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(summary_rows)

    worst = sorted(rows, key=lambda r: r["abs_vec"], reverse=True)[:args.top]
    return {
        "abaqus_rows": len(abaqus_disp),
        "stappp_rows": len(stap),
        "common_rows": len(rows),
        "relative_rows": len(rel_rows),
        "missing_abaqus_coords": missing_coord,
        "threshold_percentages": thresholds,
        "mean_rel_vec_percent": 100.0 * sum(r["rel_vec"] for r in rel_rows) / len(rel_rows),
        "p50_rel_vec_percent": 100.0 * percentile([r["rel_vec"] for r in rel_rows], 0.50),
        "p90_rel_vec_percent": 100.0 * percentile([r["rel_vec"] for r in rel_rows], 0.90),
        "out_csv": str(args.out_csv),
        "component_csv": str(args.component_csv) if args.component_csv else "",
        "worst": worst,
        "by_component": summary_rows,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dat", required=True, type=Path)
    parser.add_argument("--stappp-csv", required=True, type=Path)
    parser.add_argument("--abaqus-nodes", required=True, type=Path)
    parser.add_argument("--abaqus-displacements", required=True, type=Path)
    parser.add_argument("--out-csv", required=True, type=Path)
    parser.add_argument("--component-csv", type=Path)
    parser.add_argument("--coord-digits", type=int, default=5)
    parser.add_argument("--min-reference", type=float, default=1.0e-8)
    parser.add_argument("--thresholds", type=float, nargs="+", default=[0.10, 0.50, 1.00])
    parser.add_argument("--top", type=int, default=8)
    args = parser.parse_args()

    report = compare(args)
    print("Abaqus rows: {}".format(report["abaqus_rows"]))
    print("STAP++ rows: {}".format(report["stappp_rows"]))
    print("Common coordinate rows: {}".format(report["common_rows"]))
    print("Mean relative vector error: {:.2f}%".format(report["mean_rel_vec_percent"]))
    print("P50 relative vector error: {:.2f}%".format(report["p50_rel_vec_percent"]))
    print("P90 relative vector error: {:.2f}%".format(report["p90_rel_vec_percent"]))
    for threshold, value in report["threshold_percentages"].items():
        print("Nodes with error < {:.0f}%: {:.1f}%".format(float(threshold) * 100.0, value))
    print("\nBy component:")
    for row in report["by_component"]:
        print("  {component}: n={count}, p50={p50_rel_percent:.2f}%, p90={p90_rel_percent:.2f}%, <10%={under_10_percent:.1f}%".format(**row))
    print("\nWorst nodes:")
    for row in report["worst"]:
        print("  {instance}:{node} {component} coord={coord} abs={abs_vec:.6e} rel={rel:.2f}%".format(
            instance=row["instance"],
            node=row["node"],
            component=row["component"],
            coord=row["coord"],
            abs_vec=row["abs_vec"],
            rel=100.0 * row["rel_vec"] if row["rel_vec"] is not None else 0.0,
        ))
    print("\nWrote {}".format(report["out_csv"]))
    if report["component_csv"]:
        print("Wrote {}".format(report["component_csv"]))


if __name__ == "__main__":
    raise SystemExit(main())
