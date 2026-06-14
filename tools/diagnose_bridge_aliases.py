#!/usr/bin/env python3
"""Diagnose Bridge tie/MPC alias relationships against flattened INP assembly."""

from __future__ import annotations

import argparse
import contextlib
import csv
import importlib.util
import io
import json
from collections import Counter, defaultdict
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


def parse_dat_mpc_aliases(dat_path: Path):
    aliases = {}
    with dat_path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.strip().upper().startswith("MPC "):
                mpc_count = int(line.split()[1])
                for _ in range(mpc_count):
                    tokens = f.readline().split()
                    if not tokens:
                        continue
                    nterm = int(tokens[0])
                    rhs = float(tokens[1])
                    terms = []
                    pos = 2
                    for _term in range(nterm):
                        node = int(tokens[pos]); pos += 1
                        dof = int(tokens[pos]); pos += 1
                        coef = float(tokens[pos]); pos += 1
                        terms.append((node, dof, coef))
                    if rhs != 0.0 or nterm != 2:
                        continue
                    a, b = terms
                    if abs(a[2] - 1.0) <= 1.0e-10 and abs(b[2] + 1.0) <= 1.0e-10:
                        aliases[(a[0], a[1])] = (b[0], b[1])
                    elif abs(b[2] - 1.0) <= 1.0e-10 and abs(a[2] + 1.0) <= 1.0e-10:
                        aliases[(b[0], b[1])] = (a[0], a[1])
                break
    return aliases


def build_flattened_context(inp_path: Path, tie_mode: str, apply_tie_adjust: bool, node_order: str):
    inp2dat = load_inp2dat()
    data = inp2dat.parse_inp(str(inp_path))
    with contextlib.redirect_stdout(io.StringIO()):
        nodes, elements, node_map, _mpc = inp2dat.flatten_assembly(
            data,
            tie_mode=tie_mode,
            apply_tie_adjust=apply_tie_adjust,
            node_order=node_order,
        )

    node_instances = {}
    for (instance, _local_id), global_id in node_map.items():
        node_instances[global_id] = instance

    node_element_types = defaultdict(Counter)
    node_stap_types = defaultdict(Counter)
    for elem in elements:
        stap_type = inp2dat.element_stap_type(elem)
        for node in elem["nodes"]:
            node_element_types[node][elem["type"]] += 1
            node_stap_types[node][str(stap_type)] += 1

    return {
        "data": data,
        "nodes": nodes,
        "elements": elements,
        "node_map": node_map,
        "node_instances": node_instances,
        "node_element_types": node_element_types,
        "node_stap_types": node_stap_types,
    }


def classify_instance(instance: str) -> str:
    name = instance.upper()
    if "SUPPORTBEAM" in name:
        return "supportbeam"
    if "FLOOR" in name:
        return "floor"
    if "CABLE" in name:
        return "cable"
    if "PIER" in name:
        return "pier"
    if "RIVERBANK" in name:
        return "riverbank"
    return "other"


def summarize_counter(counter_obj: Counter):
    return ";".join("{}:{}".format(key, counter_obj[key]) for key in sorted(counter_obj))


def build_nset_lookup(data):
    lookup = defaultdict(set)
    for item in data.get("assembly", {}).get("nsets", []):
        instance = item.get("instance") or ""
        for local_node in item.get("numbers", []):
            lookup[(instance, int(local_node))].add(item["name"])
    return lookup


def build_tie_lookup(data):
    tie_map = []
    ties = data.get("assembly", {}).get("ties", [])
    for tie in ties:
        slave_surface = tie.get("slave", "")
        master_surface = tie.get("master", "")
        tie_map.append({
            "tie_name": tie.get("name", ""),
            "slave_surface": slave_surface,
            "master_surface": master_surface,
            "adjust": tie.get("adjust", True),
            "no_rotation": tie.get("no_rotation", False),
        })
    return tie_map


def build_inverse_node_map(node_map):
    inverse = {}
    for key, value in node_map.items():
        inverse[value] = key
    return inverse


def match_tie_info(slave_nsets, master_nsets, tie_infos):
    slave_surfaces = {"{}_CNS_".format(name) for name in slave_nsets}
    master_surfaces = {"{}_CNS_".format(name) for name in master_nsets}
    for item in tie_infos:
        if item["slave_surface"] in slave_surfaces and item["master_surface"] in master_surfaces:
            return item
    return {}


def build_alias_rows(aliases, context):
    rows = []
    nodes = context["nodes"]
    instances = context["node_instances"]
    element_types = context["node_element_types"]
    stap_types = context["node_stap_types"]
    inverse_node_map = build_inverse_node_map(context["node_map"])
    nset_lookup = build_nset_lookup(context["data"])
    tie_lookup = build_tie_lookup(context["data"])
    for (slave_node, slave_dof), (master_node, master_dof) in sorted(aliases.items()):
        slave_instance = instances.get(slave_node, "")
        master_instance = instances.get(master_node, "")
        slave_coord = nodes.get(slave_node, (None, None, None))
        master_coord = nodes.get(master_node, (None, None, None))
        slave_local = inverse_node_map.get(slave_node, ("", ""))[1]
        master_local = inverse_node_map.get(master_node, ("", ""))[1]
        slave_nsets = sorted(nset_lookup.get((slave_instance, int(slave_local)), set())) if slave_local != "" else []
        master_nsets = sorted(nset_lookup.get((master_instance, int(master_local)), set())) if master_local != "" else []
        tie_info = match_tie_info(slave_nsets, master_nsets, tie_lookup)
        rows.append({
            "slave_node": slave_node,
            "slave_dof": slave_dof,
            "slave_instance": slave_instance,
            "slave_local_node": slave_local,
            "slave_nsets": ";".join(slave_nsets),
            "slave_kind": classify_instance(slave_instance),
            "slave_x": slave_coord[0],
            "slave_y": slave_coord[1],
            "slave_z": slave_coord[2],
            "slave_element_types": summarize_counter(element_types.get(slave_node, Counter())),
            "slave_stap_types": summarize_counter(stap_types.get(slave_node, Counter())),
            "master_node": master_node,
            "master_dof": master_dof,
            "master_instance": master_instance,
            "master_local_node": master_local,
            "master_nsets": ";".join(master_nsets),
            "master_kind": classify_instance(master_instance),
            "master_x": master_coord[0],
            "master_y": master_coord[1],
            "master_z": master_coord[2],
            "master_element_types": summarize_counter(element_types.get(master_node, Counter())),
            "master_stap_types": summarize_counter(stap_types.get(master_node, Counter())),
            "tie_name": tie_info.get("tie_name", ""),
            "tie_slave_surface": tie_info.get("slave_surface", ""),
            "tie_master_surface": tie_info.get("master_surface", ""),
            "tie_adjust": tie_info.get("adjust", ""),
            "tie_no_rotation": tie_info.get("no_rotation", ""),
            "same_instance": int(slave_instance == master_instance),
            "same_coord": int(slave_coord == master_coord),
        })
    return rows


def write_csv(path: Path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        with path.open("w", encoding="utf-8", newline="") as f:
            writer = csv.writer(f)
            writer.writerow([
                "slave_node", "slave_dof", "slave_instance", "slave_kind",
                "slave_x", "slave_y", "slave_z", "slave_element_types", "slave_stap_types",
                "master_node", "master_dof", "master_instance", "master_kind",
                "master_x", "master_y", "master_z", "master_element_types", "master_stap_types",
                "same_instance", "same_coord",
            ])
        return
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inp", required=True, type=Path)
    parser.add_argument("--dat", required=True, type=Path)
    parser.add_argument("--out-csv", required=True, type=Path)
    parser.add_argument("--out-json", required=True, type=Path)
    parser.add_argument("--tie-mode", choices=("projected", "node", "auto"), default="auto")
    parser.add_argument("--apply-tie-adjust", action="store_true")
    parser.add_argument("--node-order", choices=("rcm", "input"), default="rcm")
    args = parser.parse_args()

    aliases = parse_dat_mpc_aliases(args.dat)
    context = build_flattened_context(args.inp, args.tie_mode, args.apply_tie_adjust, args.node_order)
    rows = build_alias_rows(aliases, context)
    write_csv(args.out_csv, rows)

    summary = {
        "alias_count": len(rows),
        "by_slave_kind": Counter(row["slave_kind"] for row in rows),
        "by_master_kind": Counter(row["master_kind"] for row in rows),
        "cross_instance_aliases": sum(1 for row in rows if not row["same_instance"]),
        "same_coord_aliases": sum(1 for row in rows if row["same_coord"]),
        "supportbeam_floor_aliases": sum(
            1 for row in rows
            if row["slave_kind"] in ("supportbeam", "floor", "cable")
            or row["master_kind"] in ("supportbeam", "floor", "cable")
        ),
    }
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    with args.out_json.open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    raise SystemExit(main())
