# STAPpp Script Classification

This document tracks the script entrypoints that are currently worth maintaining.

It does not treat these as formal entrypoints:

- `tmp/**`
- `__pycache__/`
- `*.pyc`
- generated `.vtk`
- local mirror directories
- virtual environments

## 1. Primary experiment entry

These are the current top-level scripts for active development.

| Path | Role | Current use |
| --- | --- | --- |
| `tools/run_accuracy_experiments.py` | Main experiment entry | Batch `Bridge-2` accuracy scans such as `beam-deep-scan`, `beam-connection-scan`, and `shell-solid-scan` |
| `tools/run_bridge3_fast.py` | Performance run entry | `Bridge-3` fast run / summary workflow |
| `tools/run_sparse_validation.py` | Validation entry | Sparse-path regression and fallback validation |

## 2. Comparison and diagnostic tools

These scripts support the current accuracy-attribution workflow.

| Path | Role | Current use |
| --- | --- | --- |
| `tools/compare_abaqus_stappp.py` | Displacement comparison | Abaqus vs STAP++ displacement error report |
| `tools/diagnose_bridge_aliases.py` | Alias diagnosis | Inspect slave/master alias relationships from `.dat` and `.inp` |
| `tools/analyze_bridge2_connection_errors.py` | Connection cross-reference | Link worst `Bridge-2` points to `cable-floor`, `cable-pier`, and `supportbeam-floor` connections |
| `tools/diagnose_bridge2_local_substructure.py` | Local substructure diagnosis | Inspect local elements and adjacent connections around worst `Bridge-2` points |
| `tools/audit_loads_and_sections.py` | Load/section audit | Check gravity equivalence and section/material mapping |

## 3. Conversion and preprocessing tools

| Path | Role | Current use |
| --- | --- | --- |
| `tools/inp2dat/inp2dat.py` | Main converter | Abaqus `.inp` to STAP++ `.dat` |
| `tools/rebuild_bridge_dat.py` | Rebuild helper | Rebuild bridge `.dat` files with current converter settings |

## 4. Post-processing and export tools

| Path | Role | Current use |
| --- | --- | --- |
| `tools/run_bridge_postprocess.py` | Post-process entry | Bridge result post-processing workflow |
| `tools/out2vtk/out2vtk.py` | Export tool | `.out` to VTK |
| `tools/out2vtk/dat2vtk.py` | Export tool | `.dat` to VTK |
| `tools/out2vtk/datcsv2vtk.py` | Export tool | `.dat` + displacement CSV to VTK |
| `tools/export_odb_reactions.py` | Abaqus export | Export support reactions from ODB-derived data |

## 5. Current maintenance guidance

Maintain first:

1. `tools/run_accuracy_experiments.py`
2. `tools/compare_abaqus_stappp.py`
3. `tools/diagnose_bridge_aliases.py`
4. `tools/analyze_bridge2_connection_errors.py`
5. `tools/diagnose_bridge2_local_substructure.py`
6. `tools/audit_loads_and_sections.py`

Maintain second:

- `tools/run_bridge3_fast.py`
- `tools/run_sparse_validation.py`
- `tools/rebuild_bridge_dat.py`
- `tools/out2vtk/*.py`

Not part of the current script mainline:

- ad hoc case-specific scripts under `Bridge-*` unless explicitly promoted
- historical artifacts under `tmp/`
- generated or cached files
