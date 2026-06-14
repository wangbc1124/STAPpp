#!/usr/bin/env python3
"""Rebuild the Bridge-1 VTK post-processing output from Bridge-1.out."""

import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
SCRIPT = ROOT / "tools" / "out2vtk" / "out2vtk.py"
OUT_PATH = ROOT / "data" / "Bridge-1" / "Bridge-1.out"
VTK_PATH = ROOT / "data" / "Bridge-1" / "Bridge-1.vtk"


def main():
    print("[postprocess] Bridge-1")
    cmd = [sys.executable, str(SCRIPT), str(OUT_PATH), str(VTK_PATH)]
    subprocess.run(cmd, check=True)
    print("  -> {0}".format(VTK_PATH))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
