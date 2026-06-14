# Bridge-2 Accuracy Attribution Stage Summary

## Summary

This document captures the current accuracy-attribution stage after the earlier
performance-oriented `Bridge-3` work.

Current stage goal:

- explain and reduce the displacement error of `Bridge-2` against the Abaqus baseline
- keep `PARDISO` as the large-case mainline
- avoid changing the C++ solver path until the conversion / equivalence layer is better understood

Current scope:

- `Bridge-2` is the primary gate
- `Bridge-1` and `Bridge-3` are follow-up confirmation cases
- `Bridge-4` is out of scope for this stage

## Current Repository State

Large-case solver status:

- `PARDISO` backend is available and already used as the mainline for bridge cases.
- `Bridge-1/2/3` can run under the current `PARDISO` path.

Accuracy tooling status:

- `tools/compare_abaqus_stappp.py`
- `tools/diagnose_bridge_aliases.py`
- `tools/analyze_bridge2_connection_errors.py`
- `tools/diagnose_bridge2_local_substructure.py`
- `tools/audit_loads_and_sections.py`
- `tools/run_accuracy_experiments.py`

The batch experiment script now supports:

- `beam-deep-scan`
- `beam-connection-scan`
- `shell-solid-scan`

## Current Evidence Chain

### 1. Alias / tie diagnosis

Current evidence does not support “tied to the wrong object” as the dominant explanation.

Confirmed patterns include:

- `supportbeam -> floor`
- `supportbeam -> riverbank`
- `cable -> floor`
- `cable -> pier`

The major zero-displacement symptom at some cable endpoints was partly a CSV/export issue and has already been traced through the alias-displacement backfill fix.

### 2. Connection cross-reference

The current worst `Bridge-2` errors are not explained only by tie nodes.

Important result:

- many high-error cable endpoints are directly associated with `cable-floor` or `cable-pier` links
- the largest `Part-SupportBeam-2` and `Part-Floor-1` errors are mostly not alias nodes themselves

This pushed the investigation away from “wrong tie target” and toward local stiffness / equivalence behavior.

### 3. Local substructure diagnosis

`Part-SupportBeam-2`:

- the current worst points are internal `B31` beam-band nodes
- they are not tie nodes themselves

`Part-Floor-1`:

- the current worst points are internal `S4R` shell-band nodes
- they are not tie nodes themselves

Current interpretation:

- `Part-SupportBeam-2` is the primary suspect region
- `Part-Floor-1` is a linked secondary region

### 4. Section stiffness audit

`Part-Floor-1`:

- current `S4R -> Shell4` mapping looks correct at the section/material level
- thickness, `E`, and `nu` currently appear consistent with the Abaqus input

`Part-SupportBeam-2`:

- current `B31 -> Beam3DTimoshenko` mapping uses mixed beam scaling
- this made beam parameters the next main scan surface

## Beam Deep Scan

Main result root:

- [tmp/accuracy_experiments_20260614_beam_scan](D:/Desktop/2026春/有限元/大作业/STAPpp/tmp/accuracy_experiments_20260614_beam_scan)

Main summary:

- [summary_beam-deep-scan.csv](D:/Desktop/2026春/有限元/大作业/STAPpp/tmp/accuracy_experiments_20260614_beam_scan/summary_beam-deep-scan.csv)

Baseline comparison:

- `bridge2_rotation_current_baseline`
  - `UMag relative_l2 = 2.2211814683790703`
- `bridge2_rotation_fix_r123_baseline`
  - `UMag relative_l2 = 2.1885429065135527`

Best current beam candidates:

- `beam_unscaled_euler_current`
  - `UMag relative_l2 = 2.070537214685912`
  - `supportbeam2_umag_relative_l2 = 4.263118277072766`
  - `floor1_umag_relative_l2 = 1.984402782051219`
- `beam_unscaled_sr_1p0_current`
  - `UMag relative_l2 = 2.0707160814210313`
- `beam_unscaled_sr_0p833_current`
  - `UMag relative_l2 = 2.0707516436022253`

Current trustworthy conclusion:

- removing the current mixed beam scaling gives a larger benefit than using `fix-r123` alone
- high `beam_shear_ratio` / near-Euler candidates are currently the leading direction

Current non-final conclusion:

- `fix-r123` still appears useful, but likely as a secondary correction after beam parameters are chosen

## Script Map

Primary experiment entry:

- `tools/run_accuracy_experiments.py`

Comparison:

- `tools/compare_abaqus_stappp.py`

Connection diagnostics:

- `tools/diagnose_bridge_aliases.py`
- `tools/analyze_bridge2_connection_errors.py`
- `tools/diagnose_bridge2_local_substructure.py`

Audit:

- `tools/audit_loads_and_sections.py`

## Result Directory Convention

Current active result roots:

- `tmp/accuracy_experiments_20260613/`
  - older baseline and earlier iteration results
- `tmp/accuracy_experiments_20260614/`
  - single-purpose diagnostics and local experiments
- `tmp/accuracy_experiments_20260614_beam_scan/`
  - current beam deep-scan main results

Each experiment case directory is expected to contain:

- `Bridge-2.dat`
- `Bridge-2.displacements.csv`
- `Bridge-2.out`
- `compare.json`
- `error_by_instance.csv`
- `worst_points_by_instance.csv`
- `bridge_alias_diagnosis.csv`
- `bridge_alias_diagnosis.json`
- `connection_error_crossref.csv`
- `connection_error_crossref.json`
- `local_substructure_focus.csv`
- `local_substructure_focus.json`
- `section_audit/section_stiffness_audit.csv`

## Recommended Next Order

The next fixed order is:

1. run `beam-connection-scan`
2. run `shell-solid-scan`
3. select the best `Bridge-2` candidates
4. regress on `Bridge-1`
5. expand to `Bridge-3`

## Confidence Level

Conclusions currently treated as trustworthy:

- `PARDISO` is not the main accuracy blocker for `Bridge-2`
- `Part-SupportBeam-2` is the current primary suspect region
- `Part-Floor-1` is a linked secondary region
- beam-parameter scans are more valuable right now than more tie-rotation micro-scans
- `beam`-side changes currently outperform `fix-r123` alone

Conclusions currently treated as provisional:

- the final recommended `beam_shear_ratio` is not frozen yet
- `fix-r123` may still provide meaningful combined benefit after beam candidates are narrowed
- `Bridge-2` winners may not transfer directly to `Bridge-3`
