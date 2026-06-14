#!/usr/bin/env python3
import argparse
import csv
import math
import re
from collections import defaultdict
from pathlib import Path

import inp2dat


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
            r"\s+(\d+)\s+([-+\d.eE]+)\s+([-+\d.eE]+)\s+([-+\d.eE]+)"
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
        reader = csv.reader(f)
        next(reader)
        for row in reader:
            displacements[(row[1].upper(), int(row[0]))] = (float(row[2]), float(row[3]), float(row[4]))
    return displacements


def instance_kind(instance):
    if "FLOOR" in instance:
        return "Floor"
    if "SUPPORTBEAM" in instance:
        return "SupportBeam"
    if "PIER" in instance:
        return "Pier"
    if "CABLE" in instance:
        return "Cable"
    if "RIVERBANK" in instance:
        return "RiverBank"
    return "Other"


def summarize(rows, top_n):
    by_kind = defaultdict(list)
    for row in rows:
        by_kind[instance_kind(row["instance"])].append(row)

    print(f"Mapped nodes: {len(rows)}")
    for kind in ("Floor", "SupportBeam", "Pier", "Cable", "RiverBank", "Other"):
        sub = by_kind.get(kind, [])
        if not sub:
            continue
        sum_aba = sum(abs(r["aba_u3"]) for r in sub)
        sum_stap = sum(abs(r["stap_u3"]) for r in sub)
        rms = math.sqrt(sum(r["diff_u3"] ** 2 for r in sub) / len(sub))
        rel = [
            abs(r["diff_u3"]) / abs(r["aba_u3"]) * 100.0
            for r in sub
            if abs(r["aba_u3"]) > 1.0e-8
        ]
        rel.sort()
        mean_rel = sum(rel) / len(rel) if rel else 0.0
        median_rel = rel[len(rel) // 2] if rel else 0.0
        ratio = sum_stap / sum_aba if sum_aba > 0.0 else 0.0
        print(
            f"\n{kind}: n={len(sub)} ratio_sum_abs_u3={ratio:.6f} "
            f"rms_u3={rms:.6e} mean_rel={mean_rel:.2f}% median_rel={median_rel:.2f}%"
        )
        for r in sorted(sub, key=lambda item: abs(item["diff_u3"]), reverse=True)[:top_n]:
            x, y, z = r["xyz"]
            rel_u3 = abs(r["diff_u3"]) / abs(r["aba_u3"]) * 100.0 if abs(r["aba_u3"]) > 1.0e-30 else 0.0
            print(
                f"  {r['instance']} node {r['local_node']} gid {r['gid']} "
                f"xyz=({x:.3f},{y:.3f},{z:.3f}) "
                f"aba={r['aba_u3']:.6e} stap={r['stap_u3']:.6e} "
                f"diff={r['diff_u3']:+.6e} rel={rel_u3:.2f}%"
            )


def main():
    parser = argparse.ArgumentParser(description="Compare Bridge-1 Abaqus and STAP++ U3 by instance node map.")
    parser.add_argument("--inp", default="abaqus/Bridge-1.inp")
    parser.add_argument("--out", default="Bridge-1/Bridge-1.out")
    parser.add_argument("--odb-displacements", default="abaqus/odb_displacements.csv")
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
        stap = stap_displacements[gid]
        rows.append({
            "instance": key[0],
            "local_node": key[1],
            "gid": gid,
            "xyz": nodes[gid],
            "aba_u3": aba[2],
            "stap_u3": stap[2],
            "diff_u3": stap[2] - aba[2],
        })

    summarize(rows, args.top)


if __name__ == "__main__":
    main()
