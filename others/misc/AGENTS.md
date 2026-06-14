# STAPpp Current Handoff Guide

## Goal

Current mainline objective:

- Use `PARDISO` as the default large-case solver path.
- Treat `Bridge-2` accuracy attribution as the current primary R&D task.
- Use `Bridge-2` to select the next recommended `.inp -> .dat` conversion candidates before expanding back to `Bridge-1` and `Bridge-3`.

Current secondary objective:

- Keep `Bridge-1/2/3` runnable and reproducible.
- Preserve the standard-library iterative path as a fallback and diagnostic branch, but do not optimize it first.

Out of scope for the current stage:

- `Bridge-4`
- New Abaqus baseline exports
- Reworking the C++ solver main path
- Returning to a full tie-mode matrix sweep

## Hard Constraints

- Large-case default remains `PARDISO`.
- Do not change the C++ solver default path during the current accuracy stage.
- Do not hard-code local absolute paths in scripts.
- Do not treat temporary experiment winners as converter defaults until follow-up scans and regressions are complete.
- Do not move the current focus back to standard-library iterative solvers.
- Keep dirty worktree files unless the user explicitly asks to clean them.

## Current Facts

Already implemented:

- `PARDISO` backend is available and used for large bridge cases.
- Alias CSV backfill is already fixed in source and verified via `stappp_aliasfix.exe`.
- `supportbeam-floor` rotation experiments exist:
  - `current`
  - `fix-r123`
  - `fix-r12`
  - `fix-r3`
- Accuracy diagnostics exist:
  - `tools/compare_abaqus_stappp.py`
  - `tools/diagnose_bridge_aliases.py`
  - `tools/analyze_bridge2_connection_errors.py`
  - `tools/diagnose_bridge2_local_substructure.py`
  - `tools/audit_loads_and_sections.py`
- Batch scan entry exists:
  - `tools/run_accuracy_experiments.py`
  - phases:
    - `beam-deep-scan`
    - `beam-connection-scan`
    - `shell-solid-scan`

Verified status:

- `Bridge-1/2/3` can run with `PARDISO`.
- `Bridge-2` accuracy work is no longer blocked by solver residual.
- `beam-deep-scan` has been run successfully and archived under:
  - `tmp/accuracy_experiments_20260614_beam_scan/`

## Current Conclusions

Most important conclusions so far:

- `Part-SupportBeam-2` is the current primary accuracy suspect region.
- `Part-Floor-1` is a secondary linked region.
- The worst `Bridge-2` support beam errors are not tie nodes themselves; they are internal `B31` beam-band nodes.
- `Part-Floor-1` shell parameter mapping currently looks correct at the section/material level.
- Beam-parameter experiments are currently more valuable than further tie-rotation micro-sweeps.

Current evidence from `beam-deep-scan`:

- Removing the current mixed beam scaling significantly reduces error versus:
  - `bridge2_rotation_current_baseline`
  - `bridge2_rotation_fix_r123_baseline`
- Better candidates are currently near high `beam_shear_ratio` / near-Euler behavior:
  - `beam_unscaled_euler_current`
  - `beam_unscaled_sr_1p0_current`
  - `beam_unscaled_sr_0p833_current`
- `fix-r123` helps, but its effect is currently smaller than the beam-parameter effect.

## Mainline Policy

- Large-case mainline:
  - `solver = sparse-auto`
  - `backend = pardiso`
- Current accuracy work:
  - do not modify the C++ solver path first
  - work at the conversion / equivalence / diagnostics layer first
- Fallback branch:
  - standard backend and iterative solvers remain available, but are not the current optimization target

## Execution Order

Current fixed execution order:

1. Run `beam-connection-scan`
2. Run `shell-solid-scan`
3. Pick the best `Bridge-2` candidates
4. Regress on `Bridge-1`
5. Expand to `Bridge-3`

Do not do these first:

- another full tie-mode sweep
- `Bridge-4`
- standard iterative solver tuning
- direct default changes in `inp2dat.py`

## Key Files

- `tools/run_accuracy_experiments.py`
- `tools/compare_abaqus_stappp.py`
- `tools/diagnose_bridge_aliases.py`
- `tools/analyze_bridge2_connection_errors.py`
- `tools/diagnose_bridge2_local_substructure.py`
- `tools/audit_loads_and_sections.py`
- `docs/bridge2_accuracy_stage_summary.md`

## Commands

Recommended experiment entrypoints:

```powershell
python tools\run_accuracy_experiments.py --phase beam-deep-scan --run-root tmp/accuracy_experiments_YYYYMMDD_beam_scan --exe build-vscode-ninja\stappp_aliasfix.exe
python tools\run_accuracy_experiments.py --phase beam-connection-scan --run-root tmp/accuracy_experiments_YYYYMMDD_beam_connection --exe build-vscode-ninja\stappp_aliasfix.exe
python tools\run_accuracy_experiments.py --phase shell-solid-scan --run-root tmp/accuracy_experiments_YYYYMMDD_shell_solid --exe build-vscode-ninja\stappp_aliasfix.exe
```

Existing key result roots:

- `tmp/accuracy_experiments_20260613/`
- `tmp/accuracy_experiments_20260614/`
- `tmp/accuracy_experiments_20260614_beam_scan/`

## Do Not Do

- Do not switch the current priority back to `Bridge-3` assembly performance.
- Do not treat `fix-r12` and `fix-r3` as separate future scan dimensions unless new evidence appears.
- Do not write experimental beam winners into converter defaults before:
  - `beam-connection-scan`
  - `Bridge-1` regression
- Do not clean the dirty tree as part of this stage unless the user explicitly asks.

## Background

Current stage summary:

- [docs/bridge2_accuracy_stage_summary.md](D:/Desktop/2026春/有限元/大作业/STAPpp/docs/bridge2_accuracy_stage_summary.md)

Older performance-stage background:

- [docs/bridge3_optimization_stage_summary.md](D:/Desktop/2026春/有限元/大作业/STAPpp/docs/bridge3_optimization_stage_summary.md)
