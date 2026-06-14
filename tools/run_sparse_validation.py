#!/usr/bin/env python3
"""Run Bridge-1 and Bridge-2 with the standard-library sparse iterative solver."""

import argparse
import os
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def run(cmd, cwd=ROOT, env=None, stdout=None):
    print("+ " + " ".join(str(x) for x in cmd))
    subprocess.run(cmd, cwd=cwd, env=env, stdout=stdout, stderr=subprocess.STDOUT, check=True)


def find_exe(build_dir):
    candidates = [
        build_dir / "stap++.exe",
        build_dir / "Release" / "stap++.exe",
        build_dir / "Debug" / "stap++.exe",
        build_dir / "stap++",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise FileNotFoundError("stap++ executable was not found under {0}".format(build_dir))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-dir", default="build")
    parser.add_argument("--run-dir", default="runs/validation")
    parser.add_argument("--solver", choices=("sparse-auto", "sparse-bicgstab", "sparse-gmres", "sparse-cg"), default="sparse-auto")
    parser.add_argument("--precond", choices=("none", "jacobi", "ssor", "block-jacobi", "ilu0"), default="block-jacobi")
    parser.add_argument("--scale", choices=("none", "diag"), default="diag")
    parser.add_argument("--tol", default="1e-6")
    parser.add_argument("--max-iter", default="5000")
    parser.add_argument("--skip-build", action="store_true")
    args = parser.parse_args()

    build_dir = ROOT / args.build_dir
    run_dir = ROOT / args.run_dir
    run_dir.mkdir(parents=True, exist_ok=True)

    if not args.skip_build:
        run(["cmake", "-S", str(ROOT), "-B", str(build_dir), "-DCMAKE_BUILD_TYPE=Release"])
        run(["cmake", "--build", str(build_dir), "--config", "Release"])

    exe = find_exe(build_dir)
    env = os.environ.copy()
    env["STAP_SOLVER"] = args.solver
    env["STAP_OUTPUT_MODE"] = "summary"
    env["STAP_PCG_PRECOND"] = args.precond
    env["STAP_PCG_SCALE"] = args.scale

    cases = [
        ("Bridge-1", ROOT / "Bridge-1" / "Bridge-1.dat"),
        ("Bridge-2", ROOT / "Bridge-2" / "Bridge-2.dat"),
    ]
    for name, dat_path in cases:
        if not dat_path.exists():
            raise FileNotFoundError(dat_path)
        csv_path = run_dir / (name + ".displacements.csv")
        log_path = run_dir / (name + ".stdout.log")
        with log_path.open("w", encoding="utf-8") as log:
            run([
                str(exe),
                str(dat_path),
                "--solver",
                args.solver,
                "--output",
                "summary",
                "--precond",
                args.precond,
                "--scale",
                args.scale,
                "--tol",
                str(args.tol),
                "--max-iter",
                str(args.max_iter),
                "--csv",
                str(csv_path),
            ], env=env, stdout=log)
        print("{0}: log={1} csv={2}".format(name, log_path, csv_path))


if __name__ == "__main__":
    main()
