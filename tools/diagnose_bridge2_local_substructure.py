#!/usr/bin/env python3
"""Extract local substructure diagnostics around worst Bridge-2 nodes."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

from diagnose_bridge_aliases import (
    build_alias_rows,
    build_flattened_context,
    classify_instance,
    parse_dat_mpc_aliases,
)


def load_csv(path: Path):
    with path.open("r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def int_or_none(value):
    if value is None or value == "":
        return None
    return int(float(value))


def float_or_none(value):
    if value is None or value == "":
        return None
    return float(value)


def build_node_adjacency(elements):
    node_to_elements = defaultdict(list)
    for idx, elem in enumerate(elements):
        elem_id = elem.get("id", idx + 1)
        elem["_diag_id"] = elem_id
        for node in elem["nodes"]:
            node_to_elements[node].append(elem)
    return node_to_elements


def build_alias_adjacency(alias_rows):
    node_to_alias = defaultdict(list)
    for row in alias_rows:
        slave = int(row["slave_node"])
        master = int(row["master_node"])
        node_to_alias[slave].append(row)
        node_to_alias[master].append(row)
    return node_to_alias


def summarize_elements(elements, focus_node):
    element_types = Counter()
    instance_types = Counter()
    neighbor_nodes = set()
    summary_rows = []
    for elem in elements:
        etype = elem.get("type", "")
        inst = elem.get("instance", "")
        element_types[etype] += 1
        instance_types[inst] += 1
        for nid in elem["nodes"]:
            if nid != focus_node:
                neighbor_nodes.add(nid)
        summary_rows.append({
            "element_id": elem.get("_diag_id"),
            "element_type": etype,
            "instance": inst,
            "material": elem.get("material", ""),
            "node_count": len(elem["nodes"]),
            "nodes": ";".join(str(n) for n in elem["nodes"]),
        })
    return {
        "element_type_counts": dict(element_types),
        "instance_counts": dict(instance_types),
        "neighbor_node_count": len(neighbor_nodes),
        "neighbor_nodes": sorted(neighbor_nodes),
        "rows": summary_rows,
    }


def summarize_aliases(alias_rows, focus_node):
    pair_counter = Counter()
    role_counter = Counter()
    rows = []
    for row in alias_rows:
        slave = int(row["slave_node"])
        master = int(row["master_node"])
        role = "slave" if slave == focus_node else "master" if master == focus_node else "adjacent"
        role_counter[role] += 1
        pair_name = "{}->{}".format(row["slave_kind"], row["master_kind"])
        pair_counter[pair_name] += 1
        rows.append({
            "role": role,
            "tie_name": row.get("tie_name", ""),
            "tie_no_rotation": row.get("tie_no_rotation", ""),
            "slave_node": slave,
            "slave_instance": row.get("slave_instance", ""),
            "slave_local_node": row.get("slave_local_node", ""),
            "master_node": master,
            "master_instance": row.get("master_instance", ""),
            "master_local_node": row.get("master_local_node", ""),
            "pair_name": pair_name,
            "same_coord": row.get("same_coord", ""),
        })
    return {
        "role_counts": dict(role_counter),
        "pair_counts": dict(pair_counter),
        "rows": rows,
    }


def diagnose_targets(target_rows, context, alias_rows):
    node_to_elements = build_node_adjacency(context["elements"])
    node_to_alias = build_alias_adjacency(alias_rows)
    inverse_node_map = {gid: key for key, gid in context["node_map"].items()}
    results = []
    for row in target_rows:
        focus_node = int_or_none(row.get("stappp_node"))
        if focus_node is None:
            continue
        coord = context["nodes"].get(focus_node, (None, None, None))
        local_info = inverse_node_map.get(focus_node, ("", ""))
        attached_elements = node_to_elements.get(focus_node, [])
        attached_aliases = node_to_alias.get(focus_node, [])
        elem_summary = summarize_elements(attached_elements, focus_node)
        alias_summary = summarize_aliases(attached_aliases, focus_node)
        results.append({
            "instance": row.get("instance", ""),
            "rank": int_or_none(row.get("rank")),
            "abaqus_node": row.get("abaqus_node", ""),
            "stappp_node": focus_node,
            "x": coord[0],
            "y": coord[1],
            "z": coord[2],
            "umag_abs_error": float_or_none(row.get("UMag_abs_error")),
            "instance_kind": classify_instance(row.get("instance", "")),
            "local_node": local_info[1],
            "element_type_counts": elem_summary["element_type_counts"],
            "element_instance_counts": elem_summary["instance_counts"],
            "neighbor_node_count": elem_summary["neighbor_node_count"],
            "alias_role_counts": alias_summary["role_counts"],
            "alias_pair_counts": alias_summary["pair_counts"],
            "elements": elem_summary["rows"],
            "aliases": alias_summary["rows"],
        })
    return results


def write_focus_csv(path: Path, results):
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for item in results:
        rows.append({
            "instance": item["instance"],
            "rank": item["rank"],
            "stappp_node": item["stappp_node"],
            "local_node": item["local_node"],
            "x": item["x"],
            "y": item["y"],
            "z": item["z"],
            "umag_abs_error": item["umag_abs_error"],
            "element_type_counts": json.dumps(item["element_type_counts"], ensure_ascii=False),
            "element_instance_counts": json.dumps(item["element_instance_counts"], ensure_ascii=False),
            "neighbor_node_count": item["neighbor_node_count"],
            "alias_role_counts": json.dumps(item["alias_role_counts"], ensure_ascii=False),
            "alias_pair_counts": json.dumps(item["alias_pair_counts"], ensure_ascii=False),
        })
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else [
            "instance", "rank", "stappp_node", "local_node", "x", "y", "z",
            "umag_abs_error", "element_type_counts", "element_instance_counts",
            "neighbor_node_count", "alias_role_counts", "alias_pair_counts",
        ])
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inp", required=True, type=Path)
    parser.add_argument("--dat", required=True, type=Path)
    parser.add_argument("--worst-csv", required=True, type=Path)
    parser.add_argument("--instances", nargs="+", required=True)
    parser.add_argument("--top", type=int, default=5)
    parser.add_argument("--out-csv", required=True, type=Path)
    parser.add_argument("--out-json", required=True, type=Path)
    parser.add_argument("--tie-mode", choices=("projected", "node", "auto"), default="auto")
    parser.add_argument("--apply-tie-adjust", action="store_true")
    parser.add_argument("--node-order", choices=("rcm", "input"), default="rcm")
    args = parser.parse_args()

    worst_rows = load_csv(args.worst_csv)
    selected = []
    for instance in args.instances:
        matches = [row for row in worst_rows if row.get("instance") == instance]
        matches.sort(key=lambda item: float_or_none(item.get("UMag_abs_error")) or -1.0, reverse=True)
        selected.extend(matches[: args.top])

    context = build_flattened_context(args.inp, args.tie_mode, args.apply_tie_adjust, args.node_order)
    aliases = parse_dat_mpc_aliases(args.dat)
    alias_rows = build_alias_rows(aliases, context)
    results = diagnose_targets(selected, context, alias_rows)

    write_focus_csv(args.out_csv, results)
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    with args.out_json.open("w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)

    print(json.dumps({
        "target_count": len(results),
        "instances": args.instances,
        "top_per_instance": args.top,
    }, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
