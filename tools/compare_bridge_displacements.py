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
            instance = row.get("instance") or row.get("Instance")
            node = row.get("node") or row.get("NodeID")
            u1 = row.get("U1") or row.get("DX")
            u2 = row.get("U2") or row.get("DY")
            u3 = row.get("U3") or row.get("DZ")
            if instance is None or node is None or u1 is None or u2 is None or u3 is None:
                raise ValueError(
                    "Unsupported Abaqus displacement CSV header. Expected either "
                    "node,instance,U1,U2,U3 or NodeID,Instance,DX,DY,DZ."
                )
            key = (instance.upper(), int(node))
            displacements[key] = (
                float(u1),
                float(u2),
                float(u3),
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


def tolerant_rel_error(abs_err, ref, floor):
    denom = max(abs(ref), floor)
    return abs_err / denom * 100.0 if denom > 0.0 else None


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

    components = ("U1", "U2", "U3", "UMag", "UVec")
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

    print("\nBy part group:")
    for kind in ("Floor", "SupportBeam", "Pier", "Cable", "RiverBank", "Other"):
        sub = by_kind.get(kind, [])
        if not sub:
            continue
        rel_vals = [r["rel_UVec"] for r in sub if r["rel_UVec"] is not None]
        mag_rel_vals = [r["rel_UMag"] for r in sub if r["rel_UMag"] is not None]
        tol_rel_vals = [r["tol_rel_UVec"] for r in sub if r["tol_rel_UVec"] is not None]
        rms = math.sqrt(sum(r["abs_UVec"] * r["abs_UVec"] for r in sub) / len(sub))
        print(
            "{0}: n={1} max_vec_abs={2:.6e} rms_vec_abs={3:.6e} "
            "vec_p50={4:.3f}% vec_p95={5:.3f}% mag_p95={6:.3f}% tol_vec_p95={7:.3f}%".format(
                kind,
                len(sub),
                max(r["abs_UVec"] for r in sub),
                rms,
                percentile(rel_vals, 0.50),
                percentile(rel_vals, 0.95),
                percentile(mag_rel_vals, 0.95),
                percentile(tol_rel_vals, 0.95),
            )
        )

    print("\nTop {0} nodes by displacement-vector absolute error:".format(top))
    for row in sorted(rows, key=lambda r: r["abs_UVec"], reverse=True)[:top]:
        x, y, z = row["xyz"]
        print(
            "{instance} node {local_node} gid {gid} kind={kind} "
            "xyz=({x:.3f},{y:.3f},{z:.3f}) "
            "aba=({au1:.6e},{au2:.6e},{au3:.6e}) stap=({su1:.6e},{su2:.6e},{su3:.6e}) "
            "absUVec={abs_uvec:.6e} relUVec={rel_uvec:.3f}%".format(
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
                abs_uvec=row["abs_UVec"],
                rel_uvec=row["rel_UVec"] if row["rel_UVec"] is not None else 0.0,
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
            "abs_UVec",
            "rel_U1_percent",
            "rel_U2_percent",
            "rel_U3_percent",
            "rel_UMag_percent",
            "rel_UVec_percent",
            "tol_rel_UVec_percent",
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
                "abs_UVec": row["abs_UVec"],
                "rel_U1_percent": row["rel_U1"],
                "rel_U2_percent": row["rel_U2"],
                "rel_U3_percent": row["rel_U3"],
                "rel_UMag_percent": row["rel_UMag"],
                "rel_UVec_percent": row["rel_UVec"],
                "tol_rel_UVec_percent": row["tol_rel_UVec"],
            })


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inp", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--odb-displacements", required=True)
    parser.add_argument("--csv-errors", default=None)
    parser.add_argument("--top", type=int, default=5)
    parser.add_argument(
        "--relative-floor",
        type=float,
        default=None,
        help="Reference displacement floor for tolerant relative vector error. "
             "Defaults to 0.5%% of max Abaqus displacement magnitude.",
    )
    args = parser.parse_args()

    data = inp2dat.parse_inp(args.inp)
    nodes, _elements, node_map, _mpc = inp2dat.flatten_assembly(data)
    normalized_node_map = {(inst.upper(), nid): gid for (inst, nid), gid in node_map.items()}

    stap_displacements = read_stap_displacements(args.out)
    abaqus_displacements = read_abaqus_displacements(args.odb_displacements)
    if len(stap_displacements) != len(normalized_node_map):
        print(
            "WARNING: STAP displacement count ({0}) differs from flattened INP nodes ({1}). "
            "The .out file may not match the .inp/.dat version.".format(
                len(stap_displacements), len(normalized_node_map)
            )
        )

    max_aba_mag = max((math.sqrt(sum(v * v for v in aba)) for aba in abaqus_displacements.values()), default=0.0)
    relative_floor = args.relative_floor
    if relative_floor is None:
        relative_floor = 0.005 * max_aba_mag

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
        diff_vec = math.sqrt(sum(v * v for v in diff))
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
            "abs_UVec": diff_vec,
            "rel_U1": rel_error(abs(diff[0]), aba[0]),
            "rel_U2": rel_error(abs(diff[1]), aba[1]),
            "rel_U3": rel_error(abs(diff[2]), aba[2]),
            "rel_UMag": rel_error(diff_mag, aba_mag),
            "rel_UVec": rel_error(diff_vec, aba_mag),
            "tol_rel_UVec": tolerant_rel_error(diff_vec, aba_mag, relative_floor),
        }
        rows.append(row)

    print("STAP nodes with displacement: {0}".format(len(stap_displacements)))
    print("Abaqus ODB displacement rows: {0}".format(len(abaqus_displacements)))
    print("Flattened INP nodes: {0}".format(len(normalized_node_map)))
    print("Relative error floor: {0:.6e}".format(relative_floor))
    summarize(rows, args.top)
    if args.csv_errors:
        write_error_csv(args.csv_errors, rows)
        print("\nWrote node error CSV: {0}".format(args.csv_errors))


if __name__ == "__main__":
    main()
