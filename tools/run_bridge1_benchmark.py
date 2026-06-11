#!/usr/bin/env python3
"""Run the Bridge-1 STAP++ benchmark against Abaqus ODB-exported CSV data."""

import argparse
import os
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def run(cmd, env=None, capture=False):
    print("+ " + " ".join(str(x) for x in cmd))
    if capture:
        result = subprocess.run(cmd, cwd=ROOT, env=env, text=True, capture_output=True)
        if result.returncode != 0:
            if result.stdout:
                print(result.stdout)
            if result.stderr:
                print(result.stderr, file=sys.stderr)
            result.check_returncode()
        return
    subprocess.run(cmd, cwd=ROOT, env=env, check=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case-name", default="Bridge-1-benchmark")
    parser.add_argument("--inp", default="abaqus/Bridge-1.inp")
    parser.add_argument("--odb-displacements", default="abaqus/odb_displacements.csv")
    parser.add_argument("--solid-type", default="H8RPIER", choices=("H8R", "H8RPIER"))
    parser.add_argument("--pier-instances", default="Part-Pier")
    parser.add_argument("--beam-shear-ratio", type=float, default=0.1)
    parser.add_argument("--beam-stiffness-scale", type=float, default=0.918)
    parser.add_argument("--beam-area-scale", type=float, default=1.0)
    parser.add_argument("--beam-bending-scale", type=float, default=1.15)
    parser.add_argument("--beam-torsion-scale", type=float, default=1.0)
    parser.add_argument("--beam-shear-area-scale", type=float, default=1.0)
    parser.add_argument("--solid-stiffness-scale", type=float, default=1.0)
    parser.add_argument("--pier-stiffness-scale", type=float, default=1.0)
    parser.add_argument("--pier-h8-alpha", type=float, default=None)
    parser.add_argument("--pier-h8-alpha-min", type=float, default=None)
    parser.add_argument("--h8-alpha", type=float, default=None)
    parser.add_argument("--h8-alpha-min", type=float, default=None)
    parser.add_argument("--pier-global-y-scale", type=float, default=None)
    parser.add_argument("--h8-global-y-scale", type=float, default=None)
    parser.add_argument("--pier-h8-sri", default=None)
    parser.add_argument("--pier-h8-sri-hg-blend", type=float, default=None)
    parser.add_argument("--pier-h8-fb-hg-scale", type=float, default=None)
    parser.add_argument("--h8-fb-hg-scale", type=float, default=None)
    parser.add_argument("--pier-h8-ortho-hg", default=None)
    parser.add_argument("--pier-h8-ortho-hg-scale", type=float, default=None)
    parser.add_argument("--pier-h8-ortho-hg-x-scale", type=float, default=None)
    parser.add_argument("--pier-h8-ortho-hg-y-scale", type=float, default=None)
    parser.add_argument("--pier-h8-ortho-hg-z-scale", type=float, default=None)
    parser.add_argument("--h8-ortho-hg", default=None)
    parser.add_argument("--h8-ortho-hg-scale", type=float, default=None)
    parser.add_argument("--shell-membrane-scale", type=float, default=None)
    parser.add_argument("--shell-bending-scale", type=float, default=None)
    parser.add_argument("--shell-shear-scale", type=float, default=None)
    parser.add_argument("--shell-hg-alpha", type=float, default=None)
    parser.add_argument("--tie-mode", choices=("projected", "node", "auto"), default="auto")
    parser.add_argument("--apply-tie-adjust", action="store_true")
    parser.add_argument("--node-order", choices=("rcm", "input"), default="rcm")
    parser.add_argument("--exe", default=None, help="Path to stap++.exe. Defaults to build/stap++.exe if present.")
    parser.add_argument("--top", type=int, default=20)
    args = parser.parse_args()

    case_base = ROOT / "Bridge-1" / args.case_name
    dat_path = case_base.with_suffix(".dat")
    out_path = case_base.with_suffix(".out")
    csv_path = case_base.with_name(case_base.name + "_node_errors.csv")

    convert_cmd = [
        sys.executable,
        str(ROOT / "tools" / "inp2dat" / "inp2dat.py"),
        str(ROOT / args.inp),
        str(dat_path),
        "--solid-type",
        args.solid_type,
        "--pier-instances",
        args.pier_instances,
        "--beam-shear-ratio",
        str(args.beam_shear_ratio),
        "--beam-stiffness-scale",
        str(args.beam_stiffness_scale),
        "--beam-area-scale",
        str(args.beam_area_scale),
        "--beam-bending-scale",
        str(args.beam_bending_scale),
        "--beam-torsion-scale",
        str(args.beam_torsion_scale),
        "--beam-shear-area-scale",
        str(args.beam_shear_area_scale),
        "--solid-stiffness-scale",
        str(args.solid_stiffness_scale),
        "--pier-stiffness-scale",
        str(args.pier_stiffness_scale),
        "--tie-mode",
        args.tie_mode,
        "--node-order",
        args.node_order,
    ]
    if args.apply_tie_adjust:
        convert_cmd.append("--apply-tie-adjust")
    run(convert_cmd)

    env = os.environ.copy()
    if args.pier_h8_alpha is not None:
        env["STAP_PIER_H8_ALPHA"] = str(args.pier_h8_alpha)
    if args.pier_h8_alpha_min is not None:
        env["STAP_PIER_H8_ALPHA_MIN"] = str(args.pier_h8_alpha_min)
    if args.h8_alpha is not None:
        env["STAP_H8_ALPHA"] = str(args.h8_alpha)
    if args.h8_alpha_min is not None:
        env["STAP_H8_ALPHA_MIN"] = str(args.h8_alpha_min)
    if args.pier_global_y_scale is not None:
        env["STAP_PIER_GLOBAL_Y_SCALE"] = str(args.pier_global_y_scale)
    if args.h8_global_y_scale is not None:
        env["STAP_H8_GLOBAL_Y_SCALE"] = str(args.h8_global_y_scale)
    if args.pier_h8_sri is not None:
        env["STAP_PIER_H8_SRI"] = str(args.pier_h8_sri)
    if args.pier_h8_sri_hg_blend is not None:
        env["STAP_PIER_H8_SRI_HG_BLEND"] = str(args.pier_h8_sri_hg_blend)
    if args.pier_h8_fb_hg_scale is not None:
        env["STAP_PIER_H8_FB_HG_SCALE"] = str(args.pier_h8_fb_hg_scale)
    if args.h8_fb_hg_scale is not None:
        env["STAP_H8_FB_HG_SCALE"] = str(args.h8_fb_hg_scale)
    if args.pier_h8_ortho_hg is not None:
        env["STAP_PIER_H8_ORTHO_HG"] = str(args.pier_h8_ortho_hg)
    if args.pier_h8_ortho_hg_scale is not None:
        env["STAP_PIER_H8_ORTHO_HG_SCALE"] = str(args.pier_h8_ortho_hg_scale)
    if args.pier_h8_ortho_hg_x_scale is not None:
        env["STAP_PIER_H8_ORTHO_HG_X_SCALE"] = str(args.pier_h8_ortho_hg_x_scale)
    if args.pier_h8_ortho_hg_y_scale is not None:
        env["STAP_PIER_H8_ORTHO_HG_Y_SCALE"] = str(args.pier_h8_ortho_hg_y_scale)
    if args.pier_h8_ortho_hg_z_scale is not None:
        env["STAP_PIER_H8_ORTHO_HG_Z_SCALE"] = str(args.pier_h8_ortho_hg_z_scale)
    if args.h8_ortho_hg is not None:
        env["STAP_H8_ORTHO_HG"] = str(args.h8_ortho_hg)
    if args.h8_ortho_hg_scale is not None:
        env["STAP_H8_ORTHO_HG_SCALE"] = str(args.h8_ortho_hg_scale)
    if args.shell_membrane_scale is not None:
        env["STAP_SHELL4_MEMBRANE_SCALE"] = str(args.shell_membrane_scale)
    if args.shell_bending_scale is not None:
        env["STAP_SHELL4_BENDING_SCALE"] = str(args.shell_bending_scale)
    if args.shell_shear_scale is not None:
        env["STAP_SHELL4_SHEAR_SCALE"] = str(args.shell_shear_scale)
    if args.shell_hg_alpha is not None:
        env["STAP_SHELL4_HG_ALPHA"] = str(args.shell_hg_alpha)

    exe = Path(args.exe) if args.exe else ROOT / "build" / "stap++.exe"
    if not exe.exists():
        exe = ROOT / "stap++.exe"
    run([str(exe), str(dat_path)], env=env, capture=True)
    if not out_path.exists():
        raise FileNotFoundError(out_path)

    compare_cmd = [
        sys.executable,
        str(ROOT / "tools" / "compare_bridge_displacements.py"),
        "--inp",
        str(ROOT / args.inp),
        "--out",
        str(out_path),
        "--odb-displacements",
        str(ROOT / args.odb_displacements),
        "--csv-errors",
        str(csv_path),
        "--top",
        str(args.top),
    ]
    run(compare_cmd)


if __name__ == "__main__":
    main()
