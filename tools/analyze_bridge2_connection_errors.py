#!/usr/bin/env python3
"""Cross-reference Bridge-2 worst displacement points with alias/tie connections."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path


FOCUS_PAIRS = {
    ("supportbeam", "floor"): "supportbeam-floor",
    ("cable", "floor"): "cable-floor",
    ("cable", "pier"): "cable-pier",
}


def load_csv(path: Path):
    with path.open("r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def float_or_none(value):
    if value is None or value == "":
        return None
    return float(value)


def int_or_none(value):
    if value is None or value == "":
        return None
    return int(float(value))


def build_alias_maps(rows):
    slave_map = defaultdict(list)
    master_map = defaultdict(list)
    for row in rows:
        pair_name = FOCUS_PAIRS.get((row.get("slave_kind", ""), row.get("master_kind", "")))
        if not pair_name:
            continue
        row = dict(row)
        row["pair_name"] = pair_name
        slave_map[int(row["slave_node"])].append(row)
        master_map[int(row["master_node"])].append(row)
    return slave_map, master_map


def enrich_worst_rows(worst_rows, slave_map, master_map):
    enriched = []
    for row in worst_rows:
        stappp_node = int_or_none(row.get("stappp_node"))
        if stappp_node is None:
            continue
        slave_hits = slave_map.get(stappp_node, [])
        master_hits = master_map.get(stappp_node, [])
        hit_types = []
        pair_names = []
        tie_names = []
        if slave_hits:
            hit_types.append("slave")
            pair_names.extend(hit["pair_name"] for hit in slave_hits)
            tie_names.extend(hit.get("tie_name", "") for hit in slave_hits if hit.get("tie_name"))
        if master_hits:
            hit_types.append("master")
            pair_names.extend(hit["pair_name"] for hit in master_hits)
            tie_names.extend(hit.get("tie_name", "") for hit in master_hits if hit.get("tie_name"))
        enriched.append({
            "instance": row.get("instance", ""),
            "rank": int_or_none(row.get("rank")),
            "stappp_node": stappp_node,
            "abaqus_node": row.get("abaqus_node", ""),
            "x": float_or_none(row.get("x")),
            "y": float_or_none(row.get("y")),
            "z": float_or_none(row.get("z")),
            "UMag_abs_error": float_or_none(row.get("UMag_abs_error")),
            "abaqus_U1": float_or_none(row.get("abaqus_U1")),
            "abaqus_U2": float_or_none(row.get("abaqus_U2")),
            "abaqus_U3": float_or_none(row.get("abaqus_U3")),
            "stappp_U1": float_or_none(row.get("stappp_U1")),
            "stappp_U2": float_or_none(row.get("stappp_U2")),
            "stappp_U3": float_or_none(row.get("stappp_U3")),
            "connection_hit_type": ";".join(hit_types),
            "connection_pairs": ";".join(sorted(set(pair_names))),
            "tie_names": ";".join(sorted(set(tie_names))),
            "slave_hit_count": len(slave_hits),
            "master_hit_count": len(master_hits),
        })
    return enriched


def summarize(enriched_rows):
    pair_counter = Counter()
    hit_counter = Counter()
    instance_pair_counter = defaultdict(Counter)
    for row in enriched_rows:
        hit_type = row["connection_hit_type"] or "none"
        hit_counter[hit_type] += 1
        pairs = [item for item in row["connection_pairs"].split(";") if item]
        if not pairs:
            pair_counter["none"] += 1
            instance_pair_counter[row["instance"]]["none"] += 1
            continue
        for pair in pairs:
            pair_counter[pair] += 1
            instance_pair_counter[row["instance"]][pair] += 1
    return {
        "rows": len(enriched_rows),
        "hit_type_counts": dict(hit_counter),
        "pair_counts": dict(pair_counter),
        "instance_pair_counts": {
            key: dict(counter)
            for key, counter in sorted(instance_pair_counter.items())
        },
    }


def write_csv(path: Path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        with path.open("w", encoding="utf-8", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["instance", "rank", "stappp_node", "connection_pairs", "connection_hit_type"])
        return
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--alias-csv", required=True, type=Path)
    parser.add_argument("--worst-csv", required=True, type=Path)
    parser.add_argument("--out-csv", required=True, type=Path)
    parser.add_argument("--out-json", required=True, type=Path)
    args = parser.parse_args()

    alias_rows = load_csv(args.alias_csv)
    worst_rows = load_csv(args.worst_csv)
    slave_map, master_map = build_alias_maps(alias_rows)
    enriched_rows = enrich_worst_rows(worst_rows, slave_map, master_map)
    summary = summarize(enriched_rows)

    write_csv(args.out_csv, enriched_rows)
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    with args.out_json.open("w", encoding="utf-8") as f:
        json.dump({
            "summary": summary,
            "top_rows": sorted(
                enriched_rows,
                key=lambda item: (item["UMag_abs_error"] is not None, item["UMag_abs_error"]),
                reverse=True,
            )[:40],
        }, f, indent=2, ensure_ascii=False)

    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
