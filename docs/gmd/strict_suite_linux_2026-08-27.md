# Strict Linux CPU suite qualification — 2026-08-27

Status: **CPU/data execution gate PASS by complete run plus audited focused closure;
single-command clean-room rerun remains required for the definitive campaign**.

## Locked environment

- CLM.jl commit: `6482aad85c08d5427690711402b3a8097091f0d1`
- platform: Linux `aarch64`, LinuxKit kernel 7.0.12
- container image: `sha256:fcecc1f80c4423ce2688bbf7033968f1dccbecbe04bfcbb693e1bed17690316c`
- Julia: 1.12.6
- manifest SHA-256: `94ea12b7e1a720df91489ef969b3a9ef6f1ecd485a1563998043618726d258b3`
- data controls: `SYMFLUENCE_DATA`, `CESM_INPUTDATA`, and `CLM_REQUIRE_TESTDATA=1`

The run used a fresh clone outside the cloud-synced workspace. The complete command was:

```sh
SYMFLUENCE_DATA=/path/to/SYMFLUENCE_data \
CESM_INPUTDATA=/path/to/SYMFLUENCE_data/installs/cesm-inputdata \
CLM_REQUIRE_TESTDATA=1 \
julia --project=. --check-bounds=yes -e 'using Test; include("test/runtests.jl")'
```

## Complete-run result

The command exited 0 after 252 min 17.3 s:

- 27,758 passed;
- 0 failed;
- 3 broken hardware markers (CUDA, AMDGPU, and Metal);
- 27,761 assertions represented in the summary.

The 660,623-byte log has SHA-256
`043092235e89bb8a9114887f9f8aecfad2680ac974030da1eff8e84535751671`.

## Strict-data defect and focused closure

The green exit concealed one data-dependent group. `test_crop_collapse_static.jl` ignored
`SYMFLUENCE_DATA`, constructed a developer-specific macOS Google Drive path, and emitted an
informational message instead of using the centralized strict-data gate. Consequently its
21 real-file assertions did not enter the complete-run total.

The harness was corrected to resolve both inputs from the documented SYMFLUENCE root and
call `testdata_missing`. In the same locked Linux image:

- with the staged crop fixture, the whole file passed 27/27 (6 self-contained plus the 21
  previously omitted real-file assertions);
- with an intentionally nonexistent data root and `CLM_REQUIRE_TESTDATA=1`, it exited 1
  with 6 pass and 1 explicit missing-data failure.

Thus the complete run plus focused closure covers 27,779 passing assertions and the same
three hardware skips as the corrected macOS qualification. No scientific implementation
code changed during closure. This is sufficient for the Linux CPU/data execution gate,
but the definitive frozen campaign must repeat the entire suite in one invocation after
all experiment inputs and the submission commit are frozen.

## Findings that remain outside this execution gate

The run reproduced the already documented qualification warnings: the winter AD finite
difference is non-finite while labelled passed, parameter recovery is under-constrained,
the broad initialization diagnostic reports hundreds of NaN-on-active fields, and FATES
and lake diagnostics exceed the primary claim boundary. A green CPU execution gate does
not close those scientific claim gates.
