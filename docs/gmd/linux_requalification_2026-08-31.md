# Locked-Linux requalification of C03/C05 — 2026-08-31

Status: **all three repaired gates PASS in the locked Linux ARM64 environment; full
strict-suite rerun at the same commit in progress.**

## Environment

- CLM.jl commit: `46ef5cb` (the `gmd/submission` tip carrying the 12 requalification
  commits: calibration override fixes, dual-copy state isolation, device-loop fixes,
  harness repairs, and the Metal E8/E9 records).
- Container image digest `sha256:fcecc1f80c4423ce2688bbf7033968f1dccbecbe04bfcbb693e1bed17690316c`
  (`clm-jl-gmd-linux-arm64:2026-08-27` — identical to the 2026-08-27 qualification image).
- Julia 1.12.6; Manifest SHA-256
  `94ea12b7e1a720df91489ef969b3a9ef6f1ecd485a1563998043618726d258b3` (identical to the
  2026-08-27 locked manifest).
- Fresh non-synced clone at `/private/tmp/clm-jl-linux-qual`; cached locked depot;
  `SYMFLUENCE_DATA`, `CESM_INPUTDATA`, `CLM_REQUIRE_TESTDATA=1`;
  `julia --project=. --check-bounds=yes` throughout.

## Gate results

| Gate | 2026-08-27 (pre-fix, Linux) | 2026-08-31 (post-fix, Linux) |
|---|---|---|
| `test_parameter_recovery.jl` | 15 pass / 4 fail (joint Medlyn/Vcmax zero-signal) | **19 / 19 pass** (10m01.6s) |
| `test_ad_robustness.jl` | 135 pass / 9 fail (winter non-finite references) | **150 / 150 pass** (6m21.6s) |
| `test_calibration.jl` | — | **30 / 30 pass** (6m04.8s) |

The gate test files are unchanged from the strict forms introduced in `96609ab`/`da38f8e`;
only the model-side defects they exposed were repaired (see
`calibration_requalification_2026-08-27.md` and `ad_requalification_2026-08-27.md`).

## Claim consequence

- **C03** (verified derivatives): PASS-QUALIFICATION on macOS **and** locked Linux —
  all six forcing regimes, winter included. The seasonal-cycle differentiability
  exclusion is lifted at qualification level.
- **C05** (synthetic parameter recovery): PASS-QUALIFICATION on macOS **and** locked
  Linux — scalar and joint recovery.
- Both remain subject to the frozen definitive campaign like every other claim; this
  record retires their "locked-Linux rerun pending" condition.

## Full strict-suite single-invocation rerun (same day)

The complete suite at the same commit, image, manifest, and data controls, in one
invocation (`julia --project=. --check-bounds=yes -e 'using Test; include("test/runtests.jl")'`):

- **27,785 pass / 0 fail / 3 hardware-broken (27,788 total), exit 0, 69m29.1s**;
- +6 assertions over the 2026-08-27 composite total (27,779): the winter AD guarded
  comparisons that previously skipped on non-finite references now execute;
- log SHA-256 `5cab008a8b02ac26aefb4f40702a63680c40f25a765c94d24fbd199c46adad69`.

This closes the 2026-08-27 record's "complete run in one invocation" caveat at the
current tip; the definitive campaign still reruns it through `run_campaign.py` after
the protocol freeze.
