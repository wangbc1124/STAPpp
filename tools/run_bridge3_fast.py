#!/usr/bin/env python3
"""Build and run Bridge-3 with the PARDISO mainline by default."""

import argparse
import os
import shutil
import subprocess
import sys
import time
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
    parser.add_argument("--run-dir", default="runs/Bridge-3")
    parser.add_argument("--inp", default="Bridge-3/Bridge-3.inp")
    parser.add_argument("--dat", default=None)
    parser.add_argument("--tol", default="1e-6")
    parser.add_argument("--max-iter", default="5000")
    parser.add_argument("--precond", choices=("none", "jacobi", "ssor", "block-jacobi", "ilu0"), default="block-jacobi")
    parser.add_argument("--scale", choices=("none", "diag"), default="diag")
    parser.add_argument("--solver", choices=("sparse-auto", "sparse-bicgstab", "sparse-gmres", "sparse-cg"), default="sparse-auto")
    parser.add_argument("--backend", choices=("pardiso", "standard"), default="pardiso")
    parser.add_argument("--pardiso-mtype", choices=("auto", "spd", "sym-indef", "unsym"), default="auto")
    parser.add_argument("--mpc-penalty-scale", default="1e4")
    parser.add_argument("--skip-build", action="store_true")
    parser.add_argument("--skip-convert", action="store_true")
    args = parser.parse_args()

    build_dir = ROOT / args.build_dir
    run_dir = ROOT / args.run_dir
    run_dir.mkdir(parents=True, exist_ok=True)

    dat_path = Path(args.dat) if args.dat else run_dir / "Bridge-3.dat"
    if not dat_path.is_absolute():
        dat_path = ROOT / dat_path

    if not args.skip_build:
        run(["cmake", "-S", str(ROOT), "-B", str(build_dir), "-DCMAKE_BUILD_TYPE=Release"])
        run(["cmake", "--build", str(build_dir), "--config", "Release"])

    if not args.skip_convert and not dat_path.exists():
        convert_log = run_dir / "convert.log"
        with convert_log.open("w", encoding="utf-8") as log:
            run([
                sys.executable,
                str(ROOT / "tools" / "inp2dat" / "inp2dat.py"),
                str(ROOT / args.inp),
                str(dat_path),
                "--solid-type",
                "H8RPIER",
                "--pier-instances",
                "Part-Pier",
                "--tie-mode",
                "auto",
                "--node-order",
                "rcm",
            ], stdout=log)

    if not dat_path.exists():
        raise FileNotFoundError(dat_path)

    exe = find_exe(build_dir)
    out_log = run_dir / "solve.stdout.log"
    csv_path = run_dir / "displacements.csv"
    env = os.environ.copy()
    env["STAP_SOLVER"] = args.solver
    env["STAP_BACKEND"] = args.backend
    env["STAP_OUTPUT_MODE"] = "summary"
    env["STAP_PCG_PRECOND"] = args.precond
    env["STAP_PCG_SCALE"] = args.scale
    env["STAP_MPC_PENALTY_SCALE"] = str(args.mpc_penalty_scale)

    start = time.perf_counter()
    with out_log.open("w", encoding="utf-8") as log:
        run([
            str(exe),
            str(dat_path),
            "--solver",
            str(args.solver),
            "--backend",
            str(args.backend),
            "--pardiso-mtype",
            str(args.pardiso_mtype),
            "--output",
            "summary",
            "--precond",
            str(args.precond),
            "--scale",
            str(args.scale),
            "--tol",
            str(args.tol),
            "--max-iter",
            str(args.max_iter),
            "--csv",
            str(csv_path),
        ], env=env, stdout=log)
    elapsed = time.perf_counter() - start

    print("Bridge-3 run completed.")
    print("  seconds: {0:.3f}".format(elapsed))
    print("  backend: {0}".format(args.backend))
    print("  dat: {0}".format(dat_path))
    print("  stdout log: {0}".format(out_log))
    print("  displacement csv: {0}".format(csv_path))


if __name__ == "__main__":
    main()
