#!/usr/bin/env python3
"""Compare STAP++ nodal displacements with Abaqus ODB-exported displacements."""

import argparse
import csv
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR / "inp2dat"))

import inp2dat  # noqa: E402


def read_stap_displacements(path):
    text = Path(path).read_text(encoding="utf-8", errors="ignore")
    displacements = {}
    in_disp = False
    for line in text.splitlines():
        if "D I S P L A C E M E N T S" in line:
            in_disp = True
            continue
        if not in_disp:
            continue
        match = re.match(
            r"\s*(\d+)\s+([-+\d.eE]+)\s+([-+\d.eE]+)\s+([-+\d.eE]+)"
            r"\s+([-+\d.eE]+)\s+([-+\d.eE]+)\s+([-+\d.eE]+)",
            line,
        )
        if match:
            displacements[int(match.group(1))] = tuple(float(match.group(i)) for i in range(2, 8))
        elif "S T R E S S" in line or "T O T A L" in line:
            break
    return displacements


def read_abaqus_displacements(path):
    displacements = {}
    with Path(path).open(newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            key = (row["instance"].upper(), int(row["node"]))
            displacements[key] = (
                float(row["U1"]),
                float(row["U2"]),
                float(row["U3"]),
            )
    return displacements


def instance_kind(instance):
    name = instance.upper()
    if "FLOOR" in name:
        return "Floor"
    if "SUPPORTBEAM" in name:
        return "SupportBeam"
    if "PIER" in name:
        return "Pier"
    if "CABLE" in name:
        return "Cable"
    if "RIVERBANK" in name:
        return "RiverBank"
    return "Other"


def rel_error(abs_err, ref):
    return abs_err / abs(ref) * 100.0 if abs(ref) > 1.0e-12 else None


def percentile(values, fraction):
    if not values:
        return 0.0
    values = sorted(values)
    idx = min(len(values) - 1, max(0, int(round((len(values) - 1) * fraction))))
    return values[idx]


def summarize(rows, top):
    print("Mapped nodes: {0}".format(len(rows)))
    if not rows:
        return

    components = ("U1", "U2", "U3", "UMag")
    for comp in components:
        abs_key = "abs_" + comp
        rel_key = "rel_" + comp
        abs_vals = [r[abs_key] for r in rows]
        rel_vals = [r[rel_key] for r in rows if r[rel_key] is not None]
        rms = math.sqrt(sum(v * v for v in abs_vals) / len(abs_vals))
        print(
            "{0}: max_abs={1:.6e} rms_abs={2:.6e} mean_rel={3:.3f}% "
            "p95_rel={4:.3f}% max_rel={5:.3f}%".format(
                comp,
                max(abs_vals),
                rms,
                sum(rel_vals) / len(rel_vals) if rel_vals else 0.0,
                percentile(rel_vals, 0.95),
                max(rel_vals) if rel_vals else 0.0,
            )
        )

    by_kind = defaultdict(list)
    for row in rows:
        by_kind[row["kind"]].append(row)

    print("\nBy part group, using displacement magnitude:")
    for kind in ("Floor", "SupportBeam", "Pier", "Cable", "RiverBank", "Other"):
        sub = by_kind.get(kind, [])
        if not sub:
            continue
        rel_vals = [r["rel_UMag"] for r in sub if r["rel_UMag"] is not None]
        rms = math.sqrt(sum(r["abs_UMag"] * r["abs_UMag"] for r in sub) / len(sub))
        print(
            "{0}: n={1} max_abs={2:.6e} rms_abs={3:.6e} mean_rel={4:.3f}% "
            "p95_rel={5:.3f}% max_rel={6:.3f}%".format(
                kind,
                len(sub),
                max(r["abs_UMag"] for r in sub),
                rms,
                sum(rel_vals) / len(rel_vals) if rel_vals else 0.0,
                percentile(rel_vals, 0.95),
                max(rel_vals) if rel_vals else 0.0,
            )
        )

    print("\nTop {0} nodes by displacement-magnitude absolute error:".format(top))
    for row in sorted(rows, key=lambda r: r["abs_UMag"], reverse=True)[:top]:
        x, y, z = row["xyz"]
        print(
            "{instance} node {local_node} gid {gid} kind={kind} "
            "xyz=({x:.3f},{y:.3f},{z:.3f}) "
            "aba=({au1:.6e},{au2:.6e},{au3:.6e}) stap=({su1:.6e},{su2:.6e},{su3:.6e}) "
            "absUMag={abs_umag:.6e} relUMag={rel_umag:.3f}%".format(
                instance=row["instance"],
                local_node=row["local_node"],
                gid=row["gid"],
                kind=row["kind"],
                x=x,
                y=y,
                z=z,
                au1=row["aba"][0],
                au2=row["aba"][1],
                au3=row["aba"][2],
                su1=row["stap"][0],
                su2=row["stap"][1],
                su3=row["stap"][2],
                abs_umag=row["abs_UMag"],
                rel_umag=row["rel_UMag"] if row["rel_UMag"] is not None else 0.0,
            )
        )


def write_error_csv(path, rows):
    with Path(path).open("w", newline="") as f:
        fieldnames = [
            "instance",
            "local_node",
            "gid",
            "kind",
            "x",
            "y",
            "z",
            "aba_U1",
            "aba_U2",
            "aba_U3",
            "stap_U1",
            "stap_U2",
            "stap_U3",
            "abs_U1",
            "abs_U2",
            "abs_U3",
            "abs_UMag",
            "rel_U1_percent",
            "rel_U2_percent",
            "rel_U3_percent",
            "rel_UMag_percent",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in sorted(rows, key=lambda r: r["gid"]):
            x, y, z = row["xyz"]
            writer.writerow({
                "instance": row["instance"],
                "local_node": row["local_node"],
                "gid": row["gid"],
                "kind": row["kind"],
                "x": x,
                "y": y,
                "z": z,
                "aba_U1": row["aba"][0],
                "aba_U2": row["aba"][1],
                "aba_U3": row["aba"][2],
                "stap_U1": row["stap"][0],
                "stap_U2": row["stap"][1],
                "stap_U3": row["stap"][2],
                "abs_U1": row["abs_U1"],
                "abs_U2": row["abs_U2"],
                "abs_U3": row["abs_U3"],
                "abs_UMag": row["abs_UMag"],
                "rel_U1_percent": row["rel_U1"],
                "rel_U2_percent": row["rel_U2"],
                "rel_U3_percent": row["rel_U3"],
                "rel_UMag_percent": row["rel_UMag"],
            })


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inp", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--odb-displacements", required=True)
    parser.add_argument("--csv-errors", default=None)
    parser.add_argument("--top", type=int, default=5)
    args = parser.parse_args()

    data = inp2dat.parse_inp(args.inp)
    nodes, _elements, node_map, _mpc = inp2dat.flatten_assembly(data)
    normalized_node_map = {(inst.upper(), nid): gid for (inst, nid), gid in node_map.items()}

    stap_displacements = read_stap_displacements(args.out)
    abaqus_displacements = read_abaqus_displacements(args.odb_displacements)

    rows = []
    for key, gid in normalized_node_map.items():
        if key not in abaqus_displacements or gid not in stap_displacements:
            continue
        aba = abaqus_displacements[key]
        stap6 = stap_displacements[gid]
        stap = stap6[:3]
        aba_mag = math.sqrt(sum(v * v for v in aba))
        stap_mag = math.sqrt(sum(v * v for v in stap))
        diff = tuple(stap[i] - aba[i] for i in range(3))
        diff_mag = abs(stap_mag - aba_mag)
        row = {
            "instance": key[0],
            "local_node": key[1],
            "gid": gid,
            "kind": instance_kind(key[0]),
            "xyz": nodes[gid],
            "aba": aba,
            "stap": stap,
            "abs_U1": abs(diff[0]),
            "abs_U2": abs(diff[1]),
            "abs_U3": abs(diff[2]),
            "abs_UMag": diff_mag,
            "rel_U1": rel_error(abs(diff[0]), aba[0]),
            "rel_U2": rel_error(abs(diff[1]), aba[1]),
            "rel_U3": rel_error(abs(diff[2]), aba[2]),
            "rel_UMag": rel_error(diff_mag, aba_mag),
        }
        rows.append(row)

    print("STAP nodes with displacement: {0}".format(len(stap_displacements)))
    print("Abaqus ODB displacement rows: {0}".format(len(abaqus_displacements)))
    print("Flattened INP nodes: {0}".format(len(normalized_node_map)))
    summarize(rows, args.top)
    if args.csv_errors:
        write_error_csv(args.csv_errors, rows)
        print("\nWrote node error CSV: {0}".format(args.csv_errors))


if __name__ == "__main__":
    main()
