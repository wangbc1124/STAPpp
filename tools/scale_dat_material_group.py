#!/usr/bin/env python3
"""Scale material Young's modulus for selected STAP++ element groups."""

from __future__ import annotations

import argparse
from pathlib import Path


def scale_groups(src, dst, element_type, occurrence, scale):
    lines = Path(src).read_text(encoding="utf-8", errors="replace").splitlines()
    numnp = int(lines[1].split()[0])
    i = 2 + numnp + 1
    if i < len(lines) and lines[i].strip().startswith("MPC"):
        i += 1 + int(lines[i].split()[1])

    hit = 0
    while i < len(lines):
        parts = lines[i].split()
        if not parts:
            i += 1
            continue
        if parts[0] == "0":
            break
        st, ne, nm = map(int, parts[:3])
        if st == element_type:
            hit += 1
            if occurrence == 0 or occurrence == hit:
                for mat_line in range(i + 1, i + 1 + nm):
                    mat_parts = lines[mat_line].split()
                    if len(mat_parts) < 3:
                        raise ValueError("Unsupported material line: {}".format(lines[mat_line]))
                    mat_parts[1] = "{:.6e}".format(float(mat_parts[1]) * scale)
                    lines[mat_line] = "    " + "  ".join(mat_parts)
        i += 1 + nm + ne

    if hit == 0:
        raise ValueError("No element group of type {} found in {}".format(element_type, src))
    if occurrence and occurrence > hit:
        raise ValueError("Requested occurrence {}, but only {} groups of type {} exist".format(
            occurrence, hit, element_type
        ))

    Path(dst).parent.mkdir(parents=True, exist_ok=True)
    Path(dst).write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("Wrote {} (type {}, occurrence {}, scale {})".format(dst, element_type, occurrence, scale))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("src", type=Path)
    parser.add_argument("dst", type=Path)
    parser.add_argument("--element-type", type=int, required=True)
    parser.add_argument("--occurrence", type=int, default=0,
                        help="1-based occurrence of the element type; 0 scales all matching groups.")
    parser.add_argument("--scale", type=float, required=True)
    args = parser.parse_args()
    scale_groups(args.src, args.dst, args.element_type, args.occurrence, args.scale)


if __name__ == "__main__":
    raise SystemExit(main())
