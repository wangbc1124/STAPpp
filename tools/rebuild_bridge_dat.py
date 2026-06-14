#!/usr/bin/env python3
"""Rebuild Bridge-1 to Bridge-4 .dat files from the repository .inp inputs."""

import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
CONVERTER = ROOT / "tools" / "inp2dat" / "inp2dat.py"

CASES = [
    (ROOT / "Bridge-1" / "Bridge-1.inp", ROOT / "Bridge-1" / "Bridge-1.dat"),
    (ROOT / "Bridge-2" / "Bridge-2.inp", ROOT / "Bridge-2" / "Bridge-2.dat"),
    (ROOT / "Bridge-3" / "Bridge-3.inp", ROOT / "Bridge-3.dat"),
    (ROOT / "Bridge-4.inp", ROOT / "Bridge-4.dat"),
]


def run_case(inp_path, dat_path):
    cmd = [sys.executable, str(CONVERTER), str(inp_path), str(dat_path)]
    print("+ " + " ".join(str(x) for x in cmd))
    subprocess.run(cmd, cwd=ROOT, check=True)


def main():
    for inp_path, dat_path in CASES:
        if not inp_path.exists():
            raise FileNotFoundError(inp_path)
        run_case(inp_path, dat_path)
    print("Rebuilt Bridge-1 to Bridge-4 .dat files.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
