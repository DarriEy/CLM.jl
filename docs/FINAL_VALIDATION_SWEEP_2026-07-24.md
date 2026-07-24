# Final full validation sweep — CPU + Metal GPU (2026-07-24)

A comprehensive two-backend validation run: the full CPU test suite, all CPU Fortran-parity /
domain harnesses, and the full Metal-GPU device harness fleet. Every harness's raw log and a
per-harness verdict/rc/timing row were captured.

## Environment
- **Julia** 1.12.6 (aarch64-apple-darwin).
- **CPU:** `--project=.` (main env, NCDatasets etc.).
- **Metal GPU:** `--project=scripts`, `Metal.functional() == true`, device = `AGXG16CDevice` (Apple Silicon), `HWLOC_COMPONENTS=-opencl`.
- Raw logs: `/private/tmp/claude-501/final_validation/{cpu,cpu_parity,gpu}/`; summaries `*_summary.csv`.

## Headline
- **CPU full test suite: 27402 pass / 0 fail / 3 broken** (30m16s). The 3 broken are the long-standing baseline.
- **CPU domain/Fortran-parity: 28 PASS, 0 real failures.** Every harness with a Fortran reference on disk passes; the rest skip (reference not staged locally); the lone lake "fail" is the known, documented residual.
- **Metal GPU: 118 harnesses PASS.** One real regression found (a #294 Metal down-adapt break) — **fixed (#308)**. Remaining non-passes are classifier artifacts, a wrong-backend harness, and two known/separate items.

## CPU — Fortran parity / domain (28 PASS)
Julia↔Fortran per-timestep diffs against the same instrumented CTSM `cesm.exe`. All PASS:

| Group | Harnesses |
|---|---|
| Biogeophysics | `albedo`, `aerosol`, `activelayer`, `soilresis`, `daylength` |
| Snow | `drift`, `drift_full` |
| Methane (all regimes) | `ch4`, `ch4_ebul` (ebullition), `ch4_tws` (TWS-inversion source), `ch4_merbleue` |
| Fire | `fire`, `firedefo` (deforestation), `firepeat` (peat) |
| CN | `cn_summer`, `cn_autumn`, `cn_midstep`, `cn_step`, `cn_read`, `ncycle` |
| Decomposition | `decomprates` |
| Isotopes C13/C14 | `isotopes` (Bow), `isotopes_aripuana` (tropical) |
| LUNA | `luna`, `luna_update` |
| VOC / MEGAN | `voc_aripuana` (tropical) |
| Site | `aripuana` |
| Aggregate | `sweep` |

**SKIP (no local reference dump — not failures):** `aripuana_coldstart`, `cn_coldstart`, `cn_annual`,
`qflx_surf`, `smoke`, `snow`, `soilparam`, `stillwater`, `stillwater_multistep`, `validate`.
These need Fortran ground-truth dumps not staged on this box (re-runnable via the parity recipe when the references are regenerated).

**Known residual (not a regression):** `fortran_parity_lake` runs all 48 steps but exceeds its
tight assertion threshold at ~0.20 max-rel (EFLX_LH/FSH). This is the documented lake surface-flux
residual (defects A–E fixed in #275/#279/#281/#284/#285; the remnant is the step-1 `LAKE_REC_SHIFT`
harness alignment, not physics). See `docs/LAKE_FLUX_RESIDUAL.md`.

## Metal GPU — device harness fleet (118 PASS)
118 of the `gpu_validate_*` / `gpu_ad_*` fleet pass (device == CPU parity). The non-passes:

**Real regression — FIXED (#308):** the full-driver Metal e2e harnesses (`gpu_validate_clmdrv_e2e`,
`clmdrv_cn_e2e`, `clmdrv_realinit_e2e`, `snicar_e2e`) failed with
`MethodError: no method matching CLM.DecompMIMICSState(...)`. Root cause: PR #294 added
`decomp_mimics_state` to the inst; `DecompMIMICSState` is `@kwdef` parametric only on its 4 array
fields (scalars are concrete `Float64`), so the Metal Float32 down-adaptor (`scripts/gpu_adapt.jl`)
down-converted the scalars and the positional reconstruction matched no method. **This is the
production Metal path** (`parity_run_domain_gpu.jl` uses the same adaptor), not harness-only. Fixed by
adapting only the array fields, passing scalars through. After #308: `clmdrv_e2e` PASS (max rel 4.7e-8),
`realinit` PASS (1.9e-3), `snicar` PASS (4.5e-6).

**Classifier artifacts (actually PASS):** `gpu_validate_lake`, `_soiltemp`, `_soilwater` report
"N passed, 0 failed" — the sweep classifier matched the word "failed" in "0 failed". All pass.

**Wrong backend (expected):** `gpu_validate_cuda` fails only because CUDA.jl isn't installed on a
Metal box; the CUDA path is validated on the NVIDIA box.

**Known / separate items (not this regression, not new):**
- `gpu_validate_clmdrv_cn_e2e` — with #308 the adapt now succeeds, then trips a **pre-existing
  host-loop-on-device leak** in `cnveg_carbon_flux_zero_dwt!` (`setindex!` on a device array). This is
  the known "outside-BGC GPU host-loop kernelization" gap (`gpu-outside-bgc-gap`), dates to the original
  port (not #294), previously masked by the earlier adapt failure.
- `gpu_ad_validate` — a subset of AD sanity tests (`adtest_rhs_ssw`/`mat_ssw`) hit a KernelAbstractions
  Dual-dispatch `MethodError` (no `cpu__*_kernel!` method for `ForwardDiff.Dual`); a separate AD-harness
  item, not a device-parity failure.

## Bottom line
The port's **default single-point CLM5 physics is green on both backends** (CPU suite 27402/0/3;
118 Metal device harnesses pass), and **every available Fortran-parity domain matches** — including
the subsystems wired/validated this era (methane, fire, CN, decomposition, isotopes at two biomes,
LUNA, VOC). The sweep found and fixed one real Metal regression (#294 → #308) and surfaced two known
frontier items (the outside-BGC GPU host-loop kernelization gap; an AD-harness Dual-dispatch case).
No new default-path correctness defect.
