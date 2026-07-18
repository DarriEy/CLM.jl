# Crop / irrigation / fire-crop Fortran parity — status (2026-07-18)

Backlog items **A3 (crops + crop N)**, **A4 (irrigation)** and the **A2 fire crop
branch** share a single Fortran ground-truth run: a crop-bearing spin-up with
`use_crop = .true.`, `use_cn = .true.`, `irrigate = .true.`. This documents an
attempt to generate that reference at the candidate site
`domain_Cropland_Mead_USA`, the **precise blocker** that stopped it, and the
verified state of the Julia crop/irrigation/fire-crop wiring (which side of each
gap is code vs. reference).

## TL;DR

- **Fortran crop reference: BLOCKED — not generated.** No crop-resolved surfdata
  (a `cft` dimension carrying the crop CFTs) exists for Mead or anywhere on disk,
  and none can be built: there is no `inputdata` tree and no raw crop datasets to
  drive `mksurfdata_esmf`, and no global crop surfdata to `subset_data` from. The
  migrated Mead assets are an **SP run** (`use_crop/use_cn/irrigate = .false.`,
  `cft = 2`), whose dumps are snow/hydrology-only.
- **`n_fert!` / `n_soyfix!`: NOT wired** — correctly refused (#218), still
  awaiting the reference. Not changed here.
- **No source-verifiable mis-port fix was made.** Every remaining crop gap
  (drive `plant_crop!`/`vernalization!`; wire crop N) is a *behavioural* change to
  gated crop physics whose correctness cannot be established without the Fortran
  run — blind-wiring it is exactly the trap #218 / conservation-is-not-accuracy
  warn against.

## What is on disk (Mead) — and why it is not a crop reference

The Mead domain lives only on the Drive migration copy
(`My Drive/data/SYMFLUENCE_data/domain_Cropland_Mead_USA`), not in the local
`compHydro/SYMFLUENCE_data`. Its single completed CTSM run
(`simulations/clm_crop/CLM`) is **satellite-phenology**, confirmed three ways:

| Evidence | Value |
|---|---|
| `lnd_in` | `use_crop = .false.`, `use_cn = .false.`, `irrigate = .false.`, `create_crop_landunit = .true.`, `finidat = ''` |
| build `env_run.xml` | `CLM_BLDNML_OPTS = "-bgc sp"` |
| `surfdata_clm.nc` | `cft = 2` (generic rainfed + irrigated crop — the SP layout), `natpft = 15` |
| `pdump_*_n134xx.nc` | 3 pfts, `column = 1`; only `T_SOISNO/H2OSOI_*/T_GRND/ZWT/SNOW_DEPTH/elai/tlai/T_VEG/…` — **no crop, CN, or irrigation fields** |

The prognostic crop model (`use_crop = .true.`) requires the surfdata to carry the
**crop CFTs** (the `cft = 64` "78-pft" layout with corn/soy/wheat/… — the layout
that `create_crop_landunit` + `-bgc bgc -crop` bldnml produces). A `cft = 2`
surfdata is the SP layout and is **incompatible** with `use_crop = .true.` (CTSM
aborts on the CFT-dimension mismatch). Flipping the namelist alone cannot fix
this.

### Why a crop-resolved surfdata cannot be produced here

- **No `inputdata` tree** exists in the CLM install; there is **no global crop
  surfdata** to extract a Mead single point from via
  `tools/site_and_regional/subset_data`.
- **No raw crop datasets** (`mksrf_*` PCT_CFT / land-use) are on disk, so
  `tools/mksurfdata_esmf` cannot build one either. Acquiring them is a large,
  network-heavy, uncertain download on a shared box.
- The only crop-resolved surfdatas present are CTSM **test inputs** for the wrong
  places (`surfdata_5x5_amazon_hist_78pfts_…_with_crop.nc`,
  `surfdata_1x1_mexicocityMEX_…`), which would require standing up an entirely
  different NUOPC domain / mesh / datm forcing — a fresh crop configuration, not a
  Mead run.

**To unblock A2/A3/A4, the missing prerequisite is a crop-resolved surfdata (crop
CFTs) for a crop-bearing point, plus a CN-crop cold-start spin-up long enough for
the crop to plant, accumulate GDD/HUI, grow and harvest, and crop pdump/BGCDUMP
instrumentation.** The site need not be Mead specifically — single-step parity
validates the crop *code path* against the same inputs — but it must be a
crop-CFT surfdata.

## Verified Julia wiring state (no run required)

Checked against current `main` (source, not comments) and CTSM
`CNPhenologyMod.F90` / `CropType.F90` / `irrigationMod.F90` / `CNFireLi2016Mod`.

### Crop accumulators — WIRED (the old SPVAL "instantly mature" bug is fixed)
- `crop_update_acc_vars!` (`src/types/crop.jl:308`) does real degree-day
  accumulation of `hui_patch` / `gddaccum_patch` / `gddtsoi_patch`, reset to 0
  when `!croplive`, runaccum-clamped `[0, 99999]`. Live call site
  `src/driver/clm_driver.jl:2956`, gated `config.use_crop`.
- Planting GDDs `gdd0/8/10 → gdd020/820/1020` via
  `temperature_update_acc_vars_crop_gdds!` (`src/types/temperature.jl:1104`),
  driven from `clm_driver.jl:2874`.

### Crop phenology — driver WIRED, planting NOT driven
- `cn_phenology!` is wired (`src/biogeochem/cn_driver.jl:790`, both phases);
  reaches `crop_phenology!` (`src/biogeochem/phenology.jl:936`) and crop harvest
  `cn_crop_harvest_to_product_pools!` (`phenology.jl:964`, gated `use_crop`).
- **`crop_phenology!` (`phenology.jl:2344`) is a deliberately SIMPLIFIED kernel** —
  it only sets `cphase` (planted → grainfill when `hui >= huigrain`) and the
  grain-fill `bglfr`. In CTSM, `CropPhenology` (`CNPhenologyMod.F90:2006`) also
  makes the **GDD-based planting decision, vernalization, the
  planted→leafemerge→grainfill transitions, and harvest** (sets `croplive`,
  `harvdate`).
- **`plant_crop!` (`phenology.jl:2394`) and `vernalization!` (`phenology.jl:2438`)
  are DEFINED but have no call site** — so the model never sets `croplive` on its
  own; the sowing lifecycle is not driven. Wiring this correctly needs the Fortran
  reference (planting-window thresholds, phase transitions, harvest triggers).

### Crop N (`n_fert!` / `n_soyfix!`) — UNWIRED (correctly)
- Defined + unit-tested (`src/biogeochem/n_dynamics.jl:332` / `:453`); the only
  driver "call site" is the explanatory non-call comment at
  `src/biogeochem/cn_driver.jl:403-407`. Real callers are tests / gpu-validate
  scripts only. Per #218, not wired blind. See `docs/N_CYCLE_PARITY.md` (rows
  246–247).

### Irrigation — driver calls are REAL (contrary to the old "empty stub" note)
- `set_irrig_method!` (`clm_driver.jl:781`), `calc_irrigation_needed!`
  (`clm_driver.jl:1599`, sets `irrig_rate` + `n_irrig_steps_left`),
  `calc_irrigation_fluxes!` (`clm_driver.jl:1096`), plus
  `withdraw_groundwater_irrigation!` (`hydrology_no_drainage.jl:722`). All gated
  `config.irrigate && !isempty(irr.irrig_method_patch) && !isempty(pftcon.irrigated)`.
- The A4 note "driver call sites were empty stubs — verify wired" is **resolved:
  they are genuine calls.** They remain a no-op only until an irrigated crop CFT
  patch sets `n_irrig_steps_left` — i.e. exactly what the crop reference would
  exercise. Not validatable against Fortran without it.

### Fire crop branch (`baf_crop`) — reachable, gated on crop CFT
- Computed by `_firea2016_crop_fire_kernel!` (`src/biogeochem/fire_li2016.jl:210`),
  gated `itype[p] > nc4_grass` (managed-crop CFT) `&& kmo == abm_lf[c]`; consumed
  in the crop/non-crop split at `src/biogeochem/fire_base.jl:630`. Non-zero only on
  crop-CFT patches, so `≡ 0` at Bow — needs a crop CFT run to exercise.

## Bottom line per backlog row

| Row | Status |
|---|---|
| **A3 crop pools / GDD-HUI / harvest / allocation** | BLOCKED (no crop-CFT surfdata). Accumulators verified live; `crop_phenology!` planting/vernalization/harvest not driven — needs the reference to wire+validate. |
| **A3 crop N (`n_fert!` / `n_soyfix!`)** | NOT wired — correctly refused (#218), blocked on reference. |
| **A4 irrigation** | Driver calls verified REAL + gated (A4's "empty stub" concern resolved). Behaviour un-validatable without a crop CFT run. |
| **A2 fire crop branch (`baf_crop`)** | Crop-CFT branch reachable + gated; `≡ 0` without a crop CFT — blocked on the same reference. |

Ground-truth recipe once a crop-CFT surfdata exists: reuse the
`scripts/validation/` BGCDUMP/PDUMP nstep-window instrumentation (run `cesm.exe`
plainly on macOS-26 per `macos26-lldb-xzone-brk-artifact`), then add a
`scripts/fortran_parity_crop.jl` on the `fortran_parity_common.jl` single-step
pattern — verifying the harness datm year/hour + injected inputs match the
Fortran run first (the #233/#240/#243 lesson).
