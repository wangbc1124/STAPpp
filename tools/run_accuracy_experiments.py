#!/usr/bin/env python3
"""Run accuracy experiments against the Abaqus displacement baseline."""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXE_DEFAULT = ROOT / "build-vscode-ninja" / "stap++.exe"
CONVERTER = ROOT / "tools" / "inp2dat" / "inp2dat.py"
COMPARE = ROOT / "tools" / "compare_abaqus_stappp.py"
DIAGNOSE_ALIASES = ROOT / "tools" / "diagnose_bridge_aliases.py"
ANALYZE_CONNECTIONS = ROOT / "tools" / "analyze_bridge2_connection_errors.py"
DIAGNOSE_LOCAL = ROOT / "tools" / "diagnose_bridge2_local_substructure.py"
AUDIT_LOADS_AND_SECTIONS = ROOT / "tools" / "audit_loads_and_sections.py"
MKL_RUNTIME_DIR = ROOT / "tmp" / "mkl-venv" / "Library" / "bin"
MSVC_RUNTIME_DIR = Path("C:/Program Files (x86)/Microsoft Visual Studio/2019/BuildTools/VC/Tools/MSVC/14.29.30133/bin/Hostx64/x64")

CASES = {
    "Bridge-1": {
        "inp": ROOT / "Bridge-1" / "Bridge-1.inp",
        "abaqus_csv": ROOT / "tmp" / "abaqus_baseline_20260613" / "Bridge-1" / "Bridge-1.coords_displacements.csv",
    },
    "Bridge-2": {
        "inp": ROOT / "Bridge-2" / "Bridge-2.inp",
        "abaqus_csv": ROOT / "tmp" / "abaqus_baseline_20260613" / "Bridge-2" / "Bridge-2.coords_displacements.csv",
    },
    "Bridge-3": {
        "inp": ROOT / "Bridge-3" / "Bridge-3.inp",
        "abaqus_csv": ROOT / "tmp" / "abaqus_baseline_20260613" / "Bridge-3" / "Bridge-3.coords_displacements.csv",
    },
}

TIE_EXPERIMENTS = [
    {"name": "tie_auto_noadjust", "tie_mode": "auto", "apply_tie_adjust": False},
    {"name": "tie_node_noadjust", "tie_mode": "node", "apply_tie_adjust": False},
    {"name": "tie_projected_noadjust", "tie_mode": "projected", "apply_tie_adjust": False},
    {"name": "tie_auto_adjust", "tie_mode": "auto", "apply_tie_adjust": True},
    {"name": "tie_node_adjust", "tie_mode": "node", "apply_tie_adjust": True},
]

UNIT_EXPERIMENTS = [
    {
        "name": "unit_current",
        "tie_mode": "auto",
        "apply_tie_adjust": False,
        "converter_args": [],
        "env": {},
    },
    {
        "name": "unit_solid_h8r",
        "tie_mode": "auto",
        "apply_tie_adjust": False,
        "converter_args": ["--solid-type", "H8R"],
        "env": {},
    },
    {
        "name": "unit_hg_lower",
        "tie_mode": "auto",
        "apply_tie_adjust": False,
        "converter_args": [],
        "env": {
            "STAP_H8_ALPHA": "0.003",
            "STAP_H8_ALPHA_MIN": "0.001",
            "STAP_PIER_H8_ALPHA": "0.003",
            "STAP_PIER_H8_ALPHA_MIN": "0.001",
            "STAP_PIER_H8_ORTHO_HG_SCALE": "0.005",
            "STAP_PIER_H8_ORTHO_HG_Z_SCALE": "0.1",
        },
    },
    {
        "name": "unit_no_pier_ortho_hg",
        "tie_mode": "auto",
        "apply_tie_adjust": False,
        "converter_args": [],
        "env": {
            "STAP_PIER_H8_ORTHO_HG": "0",
            "STAP_H8_ORTHO_HG": "0",
        },
    },
    {
        "name": "unit_beam_euler",
        "tie_mode": "auto",
        "apply_tie_adjust": False,
        "converter_args": ["--beam-shear-ratio", "1000000"],
        "env": {},
    },
]

BEAM_EXPERIMENTS = [
    {
        "name": "beam_unscaled_sr_0p1",
        "tie_mode": "auto",
        "apply_tie_adjust": False,
        "converter_args": [
            "--beam-stiffness-scale", "1.0",
            "--beam-bending-scale", "1.0",
            "--beam-torsion-scale", "1.0",
            "--beam-shear-ratio", "0.1",
        ],
        "env": {},
    },
    {
        "name": "beam_unscaled_sr_0p25",
        "tie_mode": "auto",
        "apply_tie_adjust": False,
        "converter_args": [
            "--beam-stiffness-scale", "1.0",
            "--beam-bending-scale", "1.0",
            "--beam-torsion-scale", "1.0",
            "--beam-shear-ratio", "0.25",
        ],
        "env": {},
    },
    {
        "name": "beam_unscaled_sr_0p5",
        "tie_mode": "auto",
        "apply_tie_adjust": False,
        "converter_args": [
            "--beam-stiffness-scale", "1.0",
            "--beam-bending-scale", "1.0",
            "--beam-torsion-scale", "1.0",
            "--beam-shear-ratio", "0.5",
        ],
        "env": {},
    },
    {
        "name": "beam_unscaled_sr_0p833",
        "tie_mode": "auto",
        "apply_tie_adjust": False,
        "converter_args": [
            "--beam-stiffness-scale", "1.0",
            "--beam-bending-scale", "1.0",
            "--beam-torsion-scale", "1.0",
            "--beam-shear-ratio", "0.833333333333",
        ],
        "env": {},
    },
    {
        "name": "beam_unscaled_sr_1p0",
        "tie_mode": "auto",
        "apply_tie_adjust": False,
        "converter_args": [
            "--beam-stiffness-scale", "1.0",
            "--beam-bending-scale", "1.0",
            "--beam-torsion-scale", "1.0",
            "--beam-shear-ratio", "1.0",
        ],
        "env": {},
    },
]

ROTATION_EXPERIMENTS = [
    {
        "name": "rotation_current",
        "tie_mode": "auto",
        "apply_tie_adjust": False,
        "converter_args": [],
        "env": {},
    },
    {
        "name": "rotation_free_beam_solid",
        "tie_mode": "auto",
        "apply_tie_adjust": False,
        "converter_args": ["--shared-rotation-mode", "free-beam-solid"],
        "env": {},
    },
]

BEAM_DEEP_SCAN_EXPERIMENTS = [
    {
        "name": "bridge2_rotation_current_baseline",
        "tie_mode": "auto",
        "apply_tie_adjust": False,
        "converter_args": [],
        "env": {},
        "supportbeam_floor_rotation_mode": "current",
    },
    {
        "name": "bridge2_rotation_fix_r123_baseline",
        "tie_mode": "auto",
        "apply_tie_adjust": False,
        "converter_args": ["--supportbeam-floor-rotation-mode", "fix-r123"],
        "env": {},
        "supportbeam_floor_rotation_mode": "fix-r123",
    },
    {
        "name": "beam_unscaled_sr_0p1_current",
        "tie_mode": "auto",
        "apply_tie_adjust": False,
        "converter_args": [
            "--beam-stiffness-scale", "1.0",
            "--beam-bending-scale", "1.0",
            "--beam-torsion-scale", "1.0",
            "--beam-shear-ratio", "0.1",
        ],
        "env": {},
        "supportbeam_floor_rotation_mode": "current",
    },
    {
        "name": "beam_unscaled_sr_0p05_current",
        "tie_mode": "auto",
        "apply_tie_adjust": False,
        "converter_args": [
            "--beam-stiffness-scale", "1.0",
            "--beam-bending-scale", "1.0",
            "--beam-torsion-scale", "1.0",
            "--beam-shear-ratio", "0.05",
        ],
        "env": {},
        "supportbeam_floor_rotation_mode": "current",
    },
    {
        "name": "beam_unscaled_sr_0p25_current",
        "tie_mode": "auto",
        "apply_tie_adjust": False,
        "converter_args": [
            "--beam-stiffness-scale", "1.0",
            "--beam-bending-scale", "1.0",
            "--beam-torsion-scale", "1.0",
            "--beam-shear-ratio", "0.25",
        ],
        "env": {},
        "supportbeam_floor_rotation_mode": "current",
    },
    {
        "name": "beam_unscaled_sr_0p5_current",
        "tie_mode": "auto",
        "apply_tie_adjust": False,
        "converter_args": [
            "--beam-stiffness-scale", "1.0",
            "--beam-bending-scale", "1.0",
            "--beam-torsion-scale", "1.0",
            "--beam-shear-ratio", "0.5",
        ],
        "env": {},
        "supportbeam_floor_rotation_mode": "current",
    },
    {
        "name": "beam_unscaled_sr_0p833_current",
        "tie_mode": "auto",
        "apply_tie_adjust": False,
        "converter_args": [
            "--beam-stiffness-scale", "1.0",
            "--beam-bending-scale", "1.0",
            "--beam-torsion-scale", "1.0",
            "--beam-shear-ratio", "0.833333333333",
        ],
        "env": {},
        "supportbeam_floor_rotation_mode": "current",
    },
    {
        "name": "beam_unscaled_sr_1p0_current",
        "tie_mode": "auto",
        "apply_tie_adjust": False,
        "converter_args": [
            "--beam-stiffness-scale", "1.0",
            "--beam-bending-scale", "1.0",
            "--beam-torsion-scale", "1.0",
            "--beam-shear-ratio", "1.0",
        ],
        "env": {},
        "supportbeam_floor_rotation_mode": "current",
    },
    {
        "name": "beam_unscaled_sr_0p1_sheararea_0p5_current",
        "tie_mode": "auto",
        "apply_tie_adjust": False,
        "converter_args": [
            "--beam-stiffness-scale", "1.0",
            "--beam-bending-scale", "1.0",
            "--beam-torsion-scale", "1.0",
            "--beam-shear-ratio", "0.1",
            "--beam-shear-area-scale", "0.5",
        ],
        "env": {},
        "supportbeam_floor_rotation_mode": "current",
    },
    {
        "name": "beam_unscaled_sr_0p1_sheararea_2p0_current",
        "tie_mode": "auto",
        "apply_tie_adjust": False,
        "converter_args": [
            "--beam-stiffness-scale", "1.0",
            "--beam-bending-scale", "1.0",
            "--beam-torsion-scale", "1.0",
            "--beam-shear-ratio", "0.1",
            "--beam-shear-area-scale", "2.0",
        ],
        "env": {},
        "supportbeam_floor_rotation_mode": "current",
    },
    {
        "name": "beam_unscaled_euler_current",
        "tie_mode": "auto",
        "apply_tie_adjust": False,
        "converter_args": [
            "--beam-stiffness-scale", "1.0",
            "--beam-bending-scale", "1.0",
            "--beam-torsion-scale", "1.0",
            "--beam-shear-ratio", "1000000",
        ],
        "env": {},
        "supportbeam_floor_rotation_mode": "current",
    },
]

BEAM_CONNECTION_SCAN_EXPERIMENTS = [
    {
        "name": "beam_unscaled_sr_0p1_current",
        "tie_mode": "auto",
        "apply_tie_adjust": False,
        "converter_args": [
            "--beam-stiffness-scale", "1.0",
            "--beam-bending-scale", "1.0",
            "--beam-torsion-scale", "1.0",
            "--beam-shear-ratio", "0.1",
        ],
        "env": {},
        "supportbeam_floor_rotation_mode": "current",
    },
    {
        "name": "beam_unscaled_sr_0p1_fix_r123",
        "tie_mode": "auto",
        "apply_tie_adjust": False,
        "converter_args": [
            "--beam-stiffness-scale", "1.0",
            "--beam-bending-scale", "1.0",
            "--beam-torsion-scale", "1.0",
            "--beam-shear-ratio", "0.1",
            "--supportbeam-floor-rotation-mode", "fix-r123",
        ],
        "env": {},
        "supportbeam_floor_rotation_mode": "fix-r123",
    },
    {
        "name": "beam_unscaled_sr_0p25_current",
        "tie_mode": "auto",
        "apply_tie_adjust": False,
        "converter_args": [
            "--beam-stiffness-scale", "1.0",
            "--beam-bending-scale", "1.0",
            "--beam-torsion-scale", "1.0",
            "--beam-shear-ratio", "0.25",
        ],
        "env": {},
        "supportbeam_floor_rotation_mode": "current",
    },
    {
        "name": "beam_unscaled_sr_0p25_fix_r123",
        "tie_mode": "auto",
        "apply_tie_adjust": False,
        "converter_args": [
            "--beam-stiffness-scale", "1.0",
            "--beam-bending-scale", "1.0",
            "--beam-torsion-scale", "1.0",
            "--beam-shear-ratio", "0.25",
            "--supportbeam-floor-rotation-mode", "fix-r123",
        ],
        "env": {},
        "supportbeam_floor_rotation_mode": "fix-r123",
    },
]

SHELL_SOLID_SCAN_EXPERIMENTS = [
    {
        "name": "beam_unscaled_sr_0p1_current_h8r",
        "tie_mode": "auto",
        "apply_tie_adjust": False,
        "converter_args": [
            "--solid-type", "H8R",
            "--beam-stiffness-scale", "1.0",
            "--beam-bending-scale", "1.0",
            "--beam-torsion-scale", "1.0",
            "--beam-shear-ratio", "0.1",
        ],
        "env": {},
        "supportbeam_floor_rotation_mode": "current",
    },
    {
        "name": "beam_unscaled_sr_0p1_current_h8rpiersoft",
        "tie_mode": "auto",
        "apply_tie_adjust": False,
        "converter_args": [
            "--beam-stiffness-scale", "1.0",
            "--beam-bending-scale", "1.0",
            "--beam-torsion-scale", "1.0",
            "--beam-shear-ratio", "0.1",
        ],
        "env": {
            "STAP_H8_ALPHA": "0.003",
            "STAP_H8_ALPHA_MIN": "0.001",
            "STAP_PIER_H8_ALPHA": "0.003",
            "STAP_PIER_H8_ALPHA_MIN": "0.001",
            "STAP_PIER_H8_ORTHO_HG_SCALE": "0.005",
            "STAP_PIER_H8_ORTHO_HG_Z_SCALE": "0.1",
        },
        "supportbeam_floor_rotation_mode": "current",
    },
]


def run(cmd, cwd=ROOT, env=None, stdout_path=None):
    print("+ " + " ".join(str(x) for x in cmd))
    start = time.perf_counter()
    if stdout_path:
        stdout_path.parent.mkdir(parents=True, exist_ok=True)
        with stdout_path.open("w", encoding="utf-8", errors="replace") as f:
            subprocess.run(cmd, cwd=cwd, env=env, stdout=f, stderr=subprocess.STDOUT, check=True)
    else:
        subprocess.run(cmd, cwd=cwd, env=env, check=True)
    return time.perf_counter() - start


def parse_solver_out(path):
    fields = {
        "checked_residual": None,
        "total_time": None,
        "csr_assembly_time": None,
        "iter_solve_time": None,
        "converted_simple_mpc": None,
        "remaining_penalty_mpc": None,
    }
    patterns = {
        "checked_residual": re.compile(r"CHECKED RELATIVE RESIDUAL\s*=\s*([0-9eE+\-.]+)"),
        "total_time": re.compile(r"TOTAL_TIME\s*=\s*([0-9eE+\-.]+)"),
        "csr_assembly_time": re.compile(r"CSR_ASSEMBLY_TIME\s*=\s*([0-9eE+\-.]+)"),
        "iter_solve_time": re.compile(r"ITER_SOLVE_TIME\s*=\s*([0-9eE+\-.]+)"),
        "converted_simple_mpc": re.compile(r"Converted simple MPC equations to DOF aliases:\s*(\d+)"),
        "remaining_penalty_mpc": re.compile(r"Remaining penalty MPC equations:\s*(\d+)"),
    }
    text = path.read_text(encoding="utf-8", errors="replace")
    for key, pattern in patterns.items():
        m = pattern.search(text)
        if m:
            fields[key] = float(m.group(1)) if key not in ("converted_simple_mpc", "remaining_penalty_mpc") else int(m.group(1))
    return fields


def load_instance_metrics(path):
    metrics = {}
    with path.open("r", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            metrics[row["instance"]] = row
    return metrics


def mean_instance_metric(instance_metrics, prefix, metric_name):
    values = []
    for instance, row in instance_metrics.items():
        if instance.startswith(prefix):
            value = row.get(metric_name)
            if value not in (None, ""):
                values.append(float(value))
    if not values:
        return None
    return sum(values) / float(len(values))


def compute_gain(value, baseline):
    if value is None or baseline in (None, 0.0):
        return None
    return (baseline - value) / baseline


def filter_args_for_audit(extra_args):
    allowed_flags = {
        "--solid-type",
        "--pier-instances",
        "--beam-shear-ratio",
        "--beam-stiffness-scale",
        "--beam-area-scale",
        "--beam-bending-scale",
        "--beam-torsion-scale",
        "--beam-shear-area-scale",
        "--solid-stiffness-scale",
        "--pier-stiffness-scale",
    }
    filtered = []
    args = list(extra_args or [])
    idx = 0
    while idx < len(args):
        flag = args[idx]
        if flag in allowed_flags and idx + 1 < len(args):
            filtered.extend([flag, args[idx + 1]])
            idx += 2
        else:
            idx += 1
    return filtered


def run_case(case, experiment, run_root, exe, extra_converter_args=None, env_overrides=None):
    cfg = CASES[case]
    case_dir = run_root / experiment["name"] / case
    case_dir.mkdir(parents=True, exist_ok=True)
    dat_path = case_dir / (case + ".dat")
    csv_path = case_dir / (case + ".displacements.csv")
    out_base = case_dir / case
    solver_stdout = case_dir / "solve.stdout.log"
    convert_stdout = case_dir / "convert.stdout.log"
    compare_json = case_dir / "compare.json"
    keypoints_csv = case_dir / "keypoints.csv"
    error_by_instance_csv = case_dir / "error_by_instance.csv"
    worst_by_instance_csv = case_dir / "worst_points_by_instance.csv"
    shared_category_csv = case_dir / "shared_rotation_audit.csv"
    alias_csv = case_dir / "bridge_alias_diagnosis.csv"
    alias_json = case_dir / "bridge_alias_diagnosis.json"
    connection_crossref_csv = case_dir / "connection_error_crossref.csv"
    connection_crossref_json = case_dir / "connection_error_crossref.json"
    local_focus_csv = case_dir / "local_substructure_focus.csv"
    local_focus_json = case_dir / "local_substructure_focus.json"
    section_audit_dir = case_dir / "section_audit"
    section_audit_csv = section_audit_dir / "section_stiffness_audit.csv"

    solid_type = "H8RPIER"
    filtered_extra_args = list(extra_converter_args or [])
    if "--solid-type" in filtered_extra_args:
        idx = filtered_extra_args.index("--solid-type")
        if idx + 1 >= len(filtered_extra_args):
            raise ValueError("--solid-type requires a value in extra_converter_args")
        solid_type = filtered_extra_args[idx + 1]
        del filtered_extra_args[idx:idx + 2]

    convert_cmd = [
        sys.executable,
        str(CONVERTER),
        str(cfg["inp"]),
        str(dat_path),
        "--solid-type",
        solid_type,
        "--pier-instances",
        "Part-Pier",
        "--tie-mode",
        experiment["tie_mode"],
        "--node-order",
        "rcm",
    ]
    if experiment["apply_tie_adjust"]:
        convert_cmd.append("--apply-tie-adjust")
    convert_cmd.extend(filtered_extra_args)
    convert_elapsed = run(convert_cmd, stdout_path=convert_stdout)

    env = os.environ.copy()
    runtime_dirs = [p for p in (MKL_RUNTIME_DIR, MSVC_RUNTIME_DIR) if p.exists()]
    if runtime_dirs:
        env["PATH"] = os.pathsep.join(str(p) for p in runtime_dirs) + os.pathsep + env.get("PATH", "")
    env.update({
        "STAP_BACKEND": "pardiso",
        "STAP_SOLVER": "sparse-auto",
        "STAP_OUTPUT_MODE": "summary",
    })
    if env_overrides:
        env.update(env_overrides)
    solve_elapsed = run([
        str(exe),
        str(dat_path),
        "--solver",
        "sparse-auto",
        "--backend",
        "pardiso",
        "--pardiso-mtype",
        "auto",
        "--output",
        "summary",
        "--csv",
        str(csv_path),
    ], env=env, stdout_path=solver_stdout)

    produced_out = dat_path.with_suffix(".out")
    final_out = out_base.with_suffix(".out")
    if produced_out.exists() and produced_out != final_out:
        shutil.move(str(produced_out), str(final_out))
    elif not final_out.exists():
        final_out = produced_out

    compare_elapsed = run([
        sys.executable,
        str(COMPARE),
        "--case",
        case,
        "--inp",
        str(cfg["inp"]),
        "--abaqus-csv",
        str(cfg["abaqus_csv"]),
        "--stappp-csv",
        str(csv_path),
        "--out-json",
        str(compare_json),
        "--keypoints-csv",
        str(keypoints_csv),
        "--error-by-instance-csv",
        str(error_by_instance_csv),
        "--worst-by-instance-csv",
        str(worst_by_instance_csv),
        "--shared-category-csv",
        str(shared_category_csv),
        "--tie-mode",
        experiment["tie_mode"],
        "--node-order",
        "rcm",
    ] + (["--apply-tie-adjust"] if experiment["apply_tie_adjust"] else []), stdout_path=case_dir / "compare.stdout.log")

    diag_elapsed = 0.0
    connection_elapsed = 0.0
    local_elapsed = 0.0
    section_elapsed = 0.0

    if case == "Bridge-2":
        diag_elapsed += run([
            sys.executable,
            str(DIAGNOSE_ALIASES),
            "--inp",
            str(cfg["inp"]),
            "--dat",
            str(dat_path),
            "--out-csv",
            str(alias_csv),
            "--out-json",
            str(alias_json),
            "--tie-mode",
            experiment["tie_mode"],
            "--node-order",
            "rcm",
        ] + (["--apply-tie-adjust"] if experiment["apply_tie_adjust"] else []), stdout_path=case_dir / "diagnose_aliases.stdout.log")

        connection_elapsed += run([
            sys.executable,
            str(ANALYZE_CONNECTIONS),
            "--alias-csv",
            str(alias_csv),
            "--worst-csv",
            str(worst_by_instance_csv),
            "--out-csv",
            str(connection_crossref_csv),
            "--out-json",
            str(connection_crossref_json),
        ], stdout_path=case_dir / "connection_crossref.stdout.log")

        local_elapsed += run([
            sys.executable,
            str(DIAGNOSE_LOCAL),
            "--inp",
            str(cfg["inp"]),
            "--dat",
            str(dat_path),
            "--worst-csv",
            str(worst_by_instance_csv),
            "--instances",
            "PART-SUPPORTBEAM-2",
            "PART-FLOOR-1",
            "--top",
            "5",
            "--out-csv",
            str(local_focus_csv),
            "--out-json",
            str(local_focus_json),
            "--tie-mode",
            experiment["tie_mode"],
            "--node-order",
            "rcm",
        ] + (["--apply-tie-adjust"] if experiment["apply_tie_adjust"] else []), stdout_path=case_dir / "local_focus.stdout.log")

        section_cmd = [
            sys.executable,
            str(AUDIT_LOADS_AND_SECTIONS),
            "--case",
            case,
            "--inp",
            str(cfg["inp"]),
            "--dat",
            str(dat_path),
            "--out-dir",
            str(section_audit_dir),
            "--tie-mode",
            experiment["tie_mode"],
            "--node-order",
            "rcm",
        ]
        if experiment["apply_tie_adjust"]:
            section_cmd.append("--apply-tie-adjust")
        section_cmd.extend(filter_args_for_audit(filtered_extra_args))
        section_elapsed += run(section_cmd, stdout_path=case_dir / "section_audit.stdout.log")

    solver = parse_solver_out(final_out)
    report = json.loads(compare_json.read_text(encoding="utf-8"))
    instance_metrics = load_instance_metrics(error_by_instance_csv)
    supportbeam2_umag = float(instance_metrics["PART-SUPPORTBEAM-2"]["UMag_relative_l2"]) if "PART-SUPPORTBEAM-2" in instance_metrics else None
    floor1_umag = float(instance_metrics["PART-FLOOR-1"]["UMag_relative_l2"]) if "PART-FLOOR-1" in instance_metrics else None
    supportbeam1_umag = float(instance_metrics["PART-SUPPORTBEAM-1"]["UMag_relative_l2"]) if "PART-SUPPORTBEAM-1" in instance_metrics else None
    cable_avg_umag = mean_instance_metric(instance_metrics, "PART-CABLE", "UMag_relative_l2")
    pier_avg_umag = mean_instance_metric(instance_metrics, "PART-PIER", "UMag_relative_l2")
    return {
        "experiment": experiment["name"],
        "case": case,
        "tie_mode": experiment["tie_mode"],
        "apply_tie_adjust": experiment["apply_tie_adjust"],
        "supportbeam_floor_rotation_mode": experiment.get("supportbeam_floor_rotation_mode", "current"),
        "convert_elapsed": convert_elapsed,
        "solve_elapsed": solve_elapsed,
        "compare_elapsed": compare_elapsed,
        "diagnose_alias_elapsed": diag_elapsed,
        "connection_crossref_elapsed": connection_elapsed,
        "local_focus_elapsed": local_elapsed,
        "section_audit_elapsed": section_elapsed,
        "checked_residual": solver["checked_residual"],
        "total_time": solver["total_time"],
        "csr_assembly_time": solver["csr_assembly_time"],
        "iter_solve_time": solver["iter_solve_time"],
        "converted_simple_mpc": solver["converted_simple_mpc"],
        "remaining_penalty_mpc": solver["remaining_penalty_mpc"],
        "mapped_common_nodes": report["alignment"]["mapped_common_nodes"],
        "missing_abaqus_mapping": report["alignment"]["missing_abaqus_mapping"],
        "missing_stappp_nodes": report["alignment"]["missing_stappp_nodes"],
        "umag_relative_l2": report["instance_node"]["UMag"]["relative_l2_error"],
        "u3_relative_l2": report["instance_node"]["U3"]["relative_l2_error"],
        "supportbeam1_umag_relative_l2": supportbeam1_umag,
        "supportbeam2_umag_relative_l2": supportbeam2_umag,
        "floor1_umag_relative_l2": floor1_umag,
        "cable_avg_umag_relative_l2": cable_avg_umag,
        "pier_avg_umag_relative_l2": pier_avg_umag,
        "compare_json": str(compare_json),
        "keypoints_csv": str(keypoints_csv),
        "error_by_instance_csv": str(error_by_instance_csv),
        "worst_by_instance_csv": str(worst_by_instance_csv),
        "shared_category_csv": str(shared_category_csv),
        "alias_csv": str(alias_csv) if alias_csv.exists() else "",
        "alias_json": str(alias_json) if alias_json.exists() else "",
        "connection_crossref_csv": str(connection_crossref_csv) if connection_crossref_csv.exists() else "",
        "connection_crossref_json": str(connection_crossref_json) if connection_crossref_json.exists() else "",
        "local_focus_csv": str(local_focus_csv) if local_focus_csv.exists() else "",
        "local_focus_json": str(local_focus_json) if local_focus_json.exists() else "",
        "section_audit_csv": str(section_audit_csv) if section_audit_csv.exists() else "",
    }


def write_summary(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "experiment", "case", "tie_mode", "apply_tie_adjust", "supportbeam_floor_rotation_mode",
        "checked_residual", "total_time", "csr_assembly_time", "iter_solve_time",
        "converted_simple_mpc", "remaining_penalty_mpc",
        "mapped_common_nodes", "missing_abaqus_mapping", "missing_stappp_nodes",
        "umag_relative_l2", "u3_relative_l2",
        "supportbeam1_umag_relative_l2", "supportbeam2_umag_relative_l2",
        "floor1_umag_relative_l2", "cable_avg_umag_relative_l2", "pier_avg_umag_relative_l2",
        "global_umag_gain_vs_baseline", "supportbeam2_umag_gain_vs_baseline",
        "floor1_umag_gain_vs_baseline", "cable_instances_avg_umag_gain",
        "convert_elapsed", "solve_elapsed", "compare_elapsed",
        "diagnose_alias_elapsed", "connection_crossref_elapsed",
        "local_focus_elapsed", "section_audit_elapsed",
        "compare_json", "keypoints_csv", "error_by_instance_csv", "worst_by_instance_csv",
        "shared_category_csv", "alias_csv", "alias_json",
        "connection_crossref_csv", "connection_crossref_json",
        "local_focus_csv", "local_focus_json", "section_audit_csv",
    ]
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def annotate_bridge2_gains(rows, baseline_experiment):
    baseline = None
    for row in rows:
        if row.get("case") == "Bridge-2" and row.get("experiment") == baseline_experiment:
            baseline = row
            break
    if baseline is None:
        return
    for row in rows:
        if row.get("case") != "Bridge-2":
            continue
        row["global_umag_gain_vs_baseline"] = compute_gain(
            row.get("umag_relative_l2"), baseline.get("umag_relative_l2")
        )
        row["supportbeam2_umag_gain_vs_baseline"] = compute_gain(
            row.get("supportbeam2_umag_relative_l2"), baseline.get("supportbeam2_umag_relative_l2")
        )
        row["floor1_umag_gain_vs_baseline"] = compute_gain(
            row.get("floor1_umag_relative_l2"), baseline.get("floor1_umag_relative_l2")
        )
        row["cable_instances_avg_umag_gain"] = compute_gain(
            row.get("cable_avg_umag_relative_l2"), baseline.get("cable_avg_umag_relative_l2")
        )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase", choices=(
        "phase-a",
        "tie-matrix",
        "unit-matrix",
        "beam-matrix",
        "rotation-matrix",
        "beam-deep-scan",
        "beam-connection-scan",
        "shell-solid-scan",
    ), default="phase-a")
    parser.add_argument("--run-root", default="tmp/accuracy_experiments_20260613")
    parser.add_argument("--exe", default=str(EXE_DEFAULT))
    parser.add_argument("--include-bridge3", action="store_true")
    args = parser.parse_args()

    exe = Path(args.exe)
    if not exe.exists():
        raise FileNotFoundError(exe)
    run_root = ROOT / args.run_root
    rows = []

    if args.phase == "phase-a":
        rows.append(run_case("Bridge-1", TIE_EXPERIMENTS[0], run_root, exe))
    elif args.phase == "tie-matrix":
        for experiment in TIE_EXPERIMENTS:
            b1 = run_case("Bridge-1", experiment, run_root, exe)
            b2 = run_case("Bridge-2", experiment, run_root, exe)
            rows.extend([b1, b2])
            if args.include_bridge3:
                ok = (
                    b1["checked_residual"] is not None and b1["checked_residual"] <= 1.0e-8
                    and b2["checked_residual"] is not None and b2["checked_residual"] <= 1.0e-8
                )
                if ok:
                    rows.append(run_case("Bridge-3", experiment, run_root, exe))
    elif args.phase == "unit-matrix":
        for experiment in UNIT_EXPERIMENTS:
            b1 = run_case(
                "Bridge-1",
                experiment,
                run_root,
                exe,
                extra_converter_args=experiment.get("converter_args"),
                env_overrides=experiment.get("env"),
            )
            b2 = run_case(
                "Bridge-2",
                experiment,
                run_root,
                exe,
                extra_converter_args=experiment.get("converter_args"),
                env_overrides=experiment.get("env"),
            )
            rows.extend([b1, b2])
            if args.include_bridge3:
                ok = (
                    b1["checked_residual"] is not None and b1["checked_residual"] <= 1.0e-8
                    and b2["checked_residual"] is not None and b2["checked_residual"] <= 1.0e-8
                )
                if ok:
                    rows.append(run_case(
                        "Bridge-3",
                        experiment,
                        run_root,
                        exe,
                        extra_converter_args=experiment.get("converter_args"),
                        env_overrides=experiment.get("env"),
                    ))
    elif args.phase == "beam-matrix":
        for experiment in BEAM_EXPERIMENTS:
            b1 = run_case(
                "Bridge-1",
                experiment,
                run_root,
                exe,
                extra_converter_args=experiment.get("converter_args"),
                env_overrides=experiment.get("env"),
            )
            b2 = run_case(
                "Bridge-2",
                experiment,
                run_root,
                exe,
                extra_converter_args=experiment.get("converter_args"),
                env_overrides=experiment.get("env"),
            )
            rows.extend([b1, b2])
            if args.include_bridge3:
                ok = (
                    b1["checked_residual"] is not None and b1["checked_residual"] <= 1.0e-8
                    and b2["checked_residual"] is not None and b2["checked_residual"] <= 1.0e-8
                    and b1["umag_relative_l2"] is not None and b1["umag_relative_l2"] <= 0.7 * 2.535876884286347
                    and b2["umag_relative_l2"] is not None and b2["umag_relative_l2"] <= 0.7 * 2.2174526174001286
                )
                if ok:
                    rows.append(run_case(
                        "Bridge-3",
                        experiment,
                        run_root,
                        exe,
                        extra_converter_args=experiment.get("converter_args"),
                        env_overrides=experiment.get("env"),
                    ))
    elif args.phase == "rotation-matrix":
        for experiment in ROTATION_EXPERIMENTS:
            b1 = run_case(
                "Bridge-1",
                experiment,
                run_root,
                exe,
                extra_converter_args=experiment.get("converter_args"),
                env_overrides=experiment.get("env"),
            )
            b2 = run_case(
                "Bridge-2",
                experiment,
                run_root,
                exe,
                extra_converter_args=experiment.get("converter_args"),
                env_overrides=experiment.get("env"),
            )
            rows.extend([b1, b2])
            if args.include_bridge3:
                ok = (
                    b1["checked_residual"] is not None and b1["checked_residual"] <= 1.0e-8
                    and b2["checked_residual"] is not None and b2["checked_residual"] <= 1.0e-8
                    and b1["umag_relative_l2"] is not None and b1["umag_relative_l2"] <= 0.7 * 2.535876884286347
                    and b2["umag_relative_l2"] is not None and b2["umag_relative_l2"] <= 0.7 * 2.2174526174001286
                )
                if ok:
                    rows.append(run_case(
                        "Bridge-3",
                        experiment,
                        run_root,
                        exe,
                        extra_converter_args=experiment.get("converter_args"),
                        env_overrides=experiment.get("env"),
                    ))
    elif args.phase == "beam-deep-scan":
        for experiment in BEAM_DEEP_SCAN_EXPERIMENTS:
            rows.append(run_case(
                "Bridge-2",
                experiment,
                run_root,
                exe,
                extra_converter_args=experiment.get("converter_args"),
                env_overrides=experiment.get("env"),
            ))
        annotate_bridge2_gains(rows, "bridge2_rotation_current_baseline")
    elif args.phase == "beam-connection-scan":
        for experiment in BEAM_CONNECTION_SCAN_EXPERIMENTS:
            rows.append(run_case(
                "Bridge-2",
                experiment,
                run_root,
                exe,
                extra_converter_args=experiment.get("converter_args"),
                env_overrides=experiment.get("env"),
            ))
        annotate_bridge2_gains(rows, "beam_unscaled_sr_0p1_current")
    elif args.phase == "shell-solid-scan":
        for experiment in SHELL_SOLID_SCAN_EXPERIMENTS:
            rows.append(run_case(
                "Bridge-2",
                experiment,
                run_root,
                exe,
                extra_converter_args=experiment.get("converter_args"),
                env_overrides=experiment.get("env"),
            ))
        annotate_bridge2_gains(rows, "beam_unscaled_sr_0p1_current_h8r")

    summary = run_root / ("summary_" + args.phase + ".csv")
    write_summary(summary, rows)
    print("Wrote summary: {}".format(summary))
    print(json.dumps(rows, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    raise SystemExit(main())
