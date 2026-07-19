# Metal-GPU re-validation — 2026-07-18

**Why:** PR #242 was the last GPU validation. Eight `src/` files changed after it,
several carrying device kernels, and none of those paths had been checked on
Metal. The Apple-Silicon dev box is returned ~2026-07-24; after that this exact
hardware is gone (the next box is CUDA), so Metal validation is the one piece of
project work that is genuinely deadline-bound. Everything Fortran-related is
re-clonable anywhere.

Precedent for treating this as urgent: **#232 found 7 GPU-adaptation regressions**
when validation-era code had gone unchecked for a stretch (CNBalanceData adaptor,
dribbler kernels, RH kernel, cn_driver guards, `is_fates_col`). Those were general
Float32/GPU bugs — they would have bitten CUDA too.

## Changed since #242

| file | change | device coverage |
|---|---|---|
| `biogeochem/methane.jl` | #243 aerenchyma grass/crop porosity; #245 hardened 5 flux floors | `gpu_validate_ch4aere_e2e`, `gpu_validate_ch4_e2e` |
| `biogeochem/cn_driver.jl` | #246 c13/c14 `c_state_update` cascade | `gpu_validate_c_iso_flux_e2e`, `gpu_validate_carbon_isotopes_e2e` |
| `biogeophys/hydrology_no_drainage.jl` | #248 **new `compute_wf2!` kernel** | **none — added here** |
| `driver/clm_driver.jl` | #254 `transient_landcover` gate | `gpu_validate_clmdrv_cn_e2e` (composite); gate itself host-side, see below |
| `driver/clm_initialize.jl`, `clm_run.jl` | #252 `use_bedrock` conditional, #225 `h2osfcflag` | host-only (init/config, no kernel) |
| `infrastructure/dyn_file_io.jl`, `dyn_pft_crop_file.jl` | #254 real-CTSM-file reading | host-only (NetCDF readers, no kernel) |

## Results (M4 Max, Float32 device vs Float64 host)

All harnesses exit 0.

| harness | worst rel(device, host) |
|---|---|
| `gpu_validate_ch4aere_e2e` | pass |
| `gpu_validate_ch4_e2e` | `conc_ch4_unsat` 1.7e-12, `finundated` 2.5e-09, `grnd_ch4_cond` 2.2e-10, `totcolch4_grc` 5.6e-12 |
| `gpu_validate_c_iso_flux_e2e` | `harvest_c_to_litr_c_col` 4.5e-08, `gap_mortality_c_to_cwdc_col` 3.3e-09 |
| `gpu_validate_carbon_isotopes_e2e` | decay/photo terms ≤ 6.9e-08; `c14_psnsun` exact |
| `gpu_validate_clmdrv_cn_e2e` | `t_soisno` 2.2e-08; `decomp_cpools_vr`, `leafc` exact |
| **`gpu_validate_computewf2_e2e`** (new) | **1.8e-08** |

**No regressions found; no `src/` change was required.**

Known non-issue: `gpu_validate_clmdrv_cn_e2e` emits a NaN gridcell water-balance
warning. That is the documented synthetic-fixture artifact (also noted in #242),
not a parity failure — the fixture does not close a water budget.

## The gap that was closed

`compute_wf2!` (#248) is a real KernelAbstractions kernel
(`_hydrond_compute_wf2_kernel!`) and had **zero** GPU coverage. It exists to
reproduce a CTSM accumulation-order quirk — `HydrologyNoDrainageMod.F90:640-699`
computes wf (0.05 m) then wf2 (0.17 m) without resetting `rwat/swat/rz`, so the
top layers are counted twice — and getting that wrong biases wf2 by O(0.01),
enough to flip it across the `max(wf2,0)` kink that `baf_peatf` sits on (it ran
7.7% high until #248). An accumulation-order subtlety is precisely the kind of
thing a device port can get silently wrong, so it now has a pin.

`scripts/gpu_validate_computewf2_e2e.jl` is written to fail loudly if it ever
becomes vacuous: it asserts the fixture actually straddles both depths (inner set
non-empty, outer strictly larger), that per-column answers differ, that a masked
column is left untouched, and that at least one column lands in the negative-wf2
kink regime.

## Host-only, by design (documented, not defects)

- `dyn_file_io.jl` / `dyn_pft_crop_file.jl` — NetCDF readers, no kernel.
- `clm_initialize.jl` / `clm_run.jl` config defaults — init-time scalars.
- The #254 `transient_landcover` gate is a host-side predicate; what it *selects*
  (CNFireFluxes' burned-fraction field) runs inside the CN composite, which is
  covered by `gpu_validate_clmdrv_cn_e2e`. A dedicated device pin for the gated
  fire branch is the one piece of coverage still worth adding.
