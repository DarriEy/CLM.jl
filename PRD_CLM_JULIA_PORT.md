# PRD: CLM → Julia Differentiable GPU Port

## Goal
Port CTSM/CLM (Community Land Model) from Fortran 90 to Julia, producing a functionally identical model that is:
1. **Numerically equivalent** to Fortran CLM (RMSE < 1e-12 per field per timestep)
2. **Differentiable** via Enzyme.jl (reverse-mode AD through full timestep)
3. **GPU-executable** via KernelAbstractions.jl (CUDA/ROCm/Metal)

## Scope
- **In scope:** All physics in `src/main/`, `src/biogeophys/`, `src/biogeochem/`, `src/soilbiogeochem/`
- **Deferred:** FATES (`src/fates/`, ~52K lines) — port as Phase 2 after core CLM works
- **Out of scope:** CIME build infrastructure, I/O (use Julia NetCDF), MPI (replace with Julia MPI.jl)
- **Lines to port:** ~150K (excluding FATES, I/O infrastructure, build system)

## Architecture Decisions (Pre-settled)

### Language & Stack
- **Julia 1.11+** with Enzyme.jl, KernelAbstractions.jl, CUDA.jl, NetCDF.jl, MPI.jl
- Package name: `CLM.jl`
- Module structure mirrors Fortran directory layout

### Data Layout
- **Keep SoA (Structure of Arrays)** — Fortran CLM already uses this pattern
- Fortran: `col%dz(begc:endc, 1:nlevgrnd)` → Julia: `dz::Matrix{Float64}` in a struct
- All state structs are plain Julia structs (not mutable where possible for AD)

### Subgrid Hierarchy
```
GridcellData → LandunitData → ColumnData → PatchData
```
- Each is a Julia struct holding SoA arrays
- Topology encoded as index arrays (same as Fortran: `col.landunit[c]` → `column_data.landunit[c]`)
- Bounds passed as `UnitRange{Int}` instead of Fortran `bounds_type`

### Filters → Masks
- **Replace filter index arrays with boolean masks** for GPU compatibility
- Fortran: `do fc = 1,num_soilc; c = filter_soilc(fc); ...`
- Julia: `@kernel function soil_kernel(mask, ...)  i = @index(Global); mask[i] || return; ...`
- Masks are static-sized (padded to max), no dynamic rebuild needed
- This is the key GPU enabler

### Snow Layers
- **Fixed max depth, masked** — pad all columns to `nlevsno` layers
- Active layer count stored as integer per column (same as Fortran `snl`)
- Operations on inactive layers produce zero (masked out)
- No dynamic array resizing

### AD Strategy
- **Phase 1:** Port 1:1, don't worry about AD
- **Phase 2:** GPU kernels with KernelAbstractions
- **Phase 3:** Smooth discontinuities for AD (separate PRD)

---

## Module Port Sequence

Each module below is **one Ralph loop iteration**. The agent:
1. Reads the Fortran source
2. Writes the Julia equivalent
3. Writes a test comparing Julia output to Fortran reference data
4. Commits on green

### Tier 0: Project Skeleton (1 iteration)

**Task 0.1:** Initialize `CLM.jl` package structure
```
CLM.jl/
├── Project.toml
├── src/
│   ├── CLM.jl              # top-level module
│   ├── constants/           # Tier 1
│   ├── types/               # Tier 2
│   ├── infrastructure/      # Tier 3
│   ├── biogeophys/          # Tier 4-6
│   ├── biogeochem/          # Tier 7-8
│   └── driver/              # Tier 9
├── test/
│   ├── runtests.jl
│   └── reference_data/      # Fortran output for validation
└── benchmarks/
```
- Generate `Project.toml` with deps: Enzyme, KernelAbstractions, CUDA, NetCDF, MPI, Test
- Create module skeleton with `include()` chain
- **Acceptance:** `using CLM` loads without error

---

### Tier 1: Constants & Pure Math (8 iterations, ~2 days)

These have ZERO CLM dependencies. Pure leaf modules.

| # | Fortran Source | Julia Target | Lines | Notes |
|---|---------------|-------------|-------|-------|
| 1.1 | `clm_varctl.F90` | `constants/varctl.jl` | 380 | Control flags → `Base.@kwdef` struct |
| 1.2 | `clm_varpar.F90` | `constants/varpar.jl` | 290 | Parameter dimensions |
| 1.3 | `clm_varcon.F90` | `constants/varcon.jl` | 460 | Physical constants (gravity, Boltzmann, etc.) |
| 1.4 | `landunit_varcon.F90` | `constants/landunit_varcon.jl` | 85 | Landunit type enums → `@enum LandunitType` |
| 1.5 | `column_varcon.F90` | `constants/column_varcon.jl` | 120 | Column type enums |
| 1.6 | `QSatMod.F90` | `biogeophys/qsat.jl` | 200 | Saturation vapor pressure. Pure math. |
| 1.7 | `TridiagonalMod.F90` | `infrastructure/tridiagonal.jl` | 150 | Tridiagonal solver. GPU-critical. |
| 1.8 | `BandDiagonalMod.F90` | `infrastructure/band_diagonal.jl` | 180 | Band diagonal solver. GPU-critical. |

**Acceptance per module:** Unit tests pass comparing Julia output to known Fortran values for 10+ test cases.

---

### Tier 2: Core Spatial Types (6 iterations, ~2 days)

The subgrid hierarchy. Everything else depends on these.

| # | Fortran Source | Julia Target | Lines | Notes |
|---|---------------|-------------|-------|-------|
| 2.1 | `GridcellType.F90` | `types/gridcell.jl` | 130 | `struct GridcellData` with SoA fields |
| 2.2 | `LandunitType.F90` | `types/landunit.jl` | 162 | `struct LandunitData` — includes urban canyon params |
| 2.3 | `ColumnType.F90` | `types/column.jl` | 252 | `struct ColumnData` — vertical levels, hillslope vars |
| 2.4 | `PatchType.F90` | `types/patch.jl` | 210 | `struct PatchData` — vegetation/FATES flags |
| 2.5 | `decompMod.F90` | `infrastructure/decomp.jl` | 350 | Domain decomposition → Julia ranges |
| 2.6 | `filterMod.F90` + `filterColMod.F90` | `infrastructure/filters.jl` | 400 | **Convert filter arrays to boolean mask arrays** |

**Key design decision for 2.6:** Filters become `BitVector` masks sized to max bounds. Example:
```julia
struct ColumnFilters
    soil::BitVector      # true for soil columns
    snow::BitVector      # true for snow-covered columns
    lake::BitVector      # true for lake columns
    nolake::BitVector    # true for non-lake columns
    urban::BitVector     # true for urban columns
    active::BitVector    # true for active columns
end
```

**Acceptance:** Can construct full subgrid hierarchy from test surface dataset. Topology indices are correct.

---

### Tier 3: Biogeophysics State Types (18 iterations, ~5 days)

All the `*Type.F90` files that hold state arrays. Mechanical translation — ideal for Ralph loop.

| # | Fortran Source | Julia Target | Lines | Notes |
|---|---------------|-------------|-------|-------|
| 3.1 | `TemperatureType.F90` | `types/temperature.jl` | 1819 | Large — Init + restart + state arrays |
| 3.2 | `EnergyFluxType.F90` | `types/energy_flux.jl` | 1053 | Energy flux state |
| 3.3 | `SoilStateType.F90` | `types/soil_state.jl` | ~700 | Soil properties |
| 3.4 | `SoilHydrologyType.F90` | `types/soil_hydrology.jl` | ~600 | Soil water state |
| 3.5 | `CanopyStateType.F90` | `types/canopy_state.jl` | ~600 | Canopy state |
| 3.6 | `LakeStateType.F90` | `types/lake_state.jl` | ~400 | Lake state |
| 3.7 | `SurfaceAlbedoType.F90` | `types/surface_albedo.jl` | ~500 | Albedo state |
| 3.8 | `SolarAbsorbedType.F90` | `types/solar_absorbed.jl` | ~300 | Solar absorption |
| 3.9 | `UrbanParamsType.F90` | `types/urban_params.jl` | ~400 | Urban parameters |
| 3.10 | `WaterInfoBaseType.F90` | `types/water_info_base.jl` | ~100 | Water info base |
| 3.11 | `WaterStateType.F90` | `types/water_state.jl` | ~500 | Generic water state |
| 3.12 | `WaterFluxType.F90` | `types/water_flux.jl` | ~500 | Generic water flux |
| 3.13 | `WaterDiagnosticBulkType.F90` | `types/water_diagnostic_bulk.jl` | 1257 | Bulk water diagnostics |
| 3.14 | `WaterStateBulkType.F90` | `types/water_state_bulk.jl` | ~800 | Bulk water state |
| 3.15 | `WaterFluxBulkType.F90` | `types/water_flux_bulk.jl` | ~800 | Bulk water flux |
| 3.16 | `WaterBalanceType.F90` | `types/water_balance.jl` | ~400 | Water balance tracking |
| 3.17 | `WaterType.F90` | `types/water.jl` | 1082 | Master water container |
| 3.18 | `FrictionVelocityMod.F90` | `types/friction_velocity.jl` | 1200 | Type + computation |

**Pattern for each:** Fortran `pointer` arrays → Julia `Vector{Float64}` / `Matrix{Float64}` fields. Init subroutines → inner constructors or factory functions. Restart read/write → separate I/O module (deferred).

**Acceptance:** Can instantiate each type with correct dimensions. Memory layout matches Fortran.

---

### Tier 4: Biogeochemistry State Types (12 iterations, ~4 days)

| # | Fortran Source | Julia Target | Lines |
|---|---------------|-------------|-------|
| 4.1 | `CNVegStateType.F90` | `types/cn_veg_state.jl` | ~600 |
| 4.2 | `CNVegCarbonStateType.F90` | `types/cn_veg_carbon_state.jl` | 4876 |
| 4.3 | `CNVegNitrogenStateType.F90` | `types/cn_veg_nitrogen_state.jl` | 2511 |
| 4.4 | `CNVegCarbonFluxType.F90` | `types/cn_veg_carbon_flux.jl` | 5566 |
| 4.5 | `CNVegNitrogenFluxType.F90` | `types/cn_veg_nitrogen_flux.jl` | 2755 |
| 4.6 | `SoilBiogeochemCarbonStateType.F90` | `types/soil_bgc_carbon_state.jl` | 1726 |
| 4.7 | `SoilBiogeochemNitrogenStateType.F90` | `types/soil_bgc_nitrogen_state.jl` | 1426 |
| 4.8 | `SoilBiogeochemCarbonFluxType.F90` | `types/soil_bgc_carbon_flux.jl` | 1022 |
| 4.9 | `SoilBiogeochemNitrogenFluxType.F90` | `types/soil_bgc_nitrogen_flux.jl` | 1253 |
| 4.10 | `SoilBiogeochemStateType.F90` | `types/soil_bgc_state.jl` | ~500 |
| 4.11 | `CropType.F90` | `types/crop.jl` | 1017 |
| 4.12 | `CNSharedParamsMod.F90` | `types/cn_shared_params.jl` | ~300 |

**Note:** 4.2 and 4.4 are the largest files in the entire biogeochem directory (4876 and 5566 lines). These are mostly array declarations and init code — mechanical but long. May need 2 iterations each.

---

### Tier 5: Radiation & Surface Physics (10 iterations, ~4 days)

First real physics. Each module gets a GPU kernel version.

| # | Fortran Source | Julia Target | Lines | GPU? |
|---|---------------|-------------|-------|------|
| 5.1 | `DaylengthMod.F90` | `biogeophys/daylength.jl` | ~200 | Trivial |
| 5.2 | `SurfaceAlbedoMod.F90` | `biogeophys/surface_albedo.jl` | 1781 | Yes — per-patch |
| 5.3 | `UrbanAlbedoMod.F90` | `biogeophys/urban_albedo.jl` | 1305 | Yes — per-column |
| 5.4 | `SurfaceRadiationMod.F90` | `biogeophys/surface_radiation.jl` | 1027 | Yes — per-patch |
| 5.5 | `UrbanRadiationMod.F90` | `biogeophys/urban_radiation.jl` | ~500 | Yes |
| 5.6 | `SnowSnicarMod.F90` | `biogeophys/snow_snicar.jl` | 2264 | Yes — per-column |
| 5.7 | `AerosolMod.F90` | `biogeophys/aerosol.jl` | ~600 | Yes |
| 5.8 | `SurfaceHumidityMod.F90` | `biogeophys/surface_humidity.jl` | ~200 | Trivial |
| 5.9 | `SurfaceResistanceMod.F90` | `biogeophys/surface_resistance.jl` | ~400 | Yes |
| 5.10 | `SoilMoistStressMod.F90` | `biogeophys/soil_moist_stress.jl` | ~400 | Yes |

**Pattern:** Each Fortran subroutine becomes a Julia function AND a `@kernel` function. CPU version calls kernel on CPU backend. GPU version dispatches to CUDA.

```julia
@kernel function surface_radiation_kernel!(
    absorbed_sw, @Const(albedo), @Const(incoming_sw), @Const(mask)
)
    i = @index(Global)
    mask[i] || return
    absorbed_sw[i] = incoming_sw[i] * (1.0 - albedo[i])
end
```

**Acceptance:** Per-field max absolute error < 1e-12 vs Fortran reference for a standard test case (1-month global 2° run).

---

### Tier 6: Temperature & Hydrology Core (12 iterations, ~5 days)

The expensive physics. These are the GPU payoff modules.

| # | Fortran Source | Julia Target | Lines | GPU Priority |
|---|---------------|-------------|-------|-------------|
| 6.1 | `SoilTemperatureMod.F90` | `biogeophys/soil_temperature.jl` | 2974 | **HIGH** — tridiagonal solve per column |
| 6.2 | `LakeTemperatureMod.F90` | `biogeophys/lake_temperature.jl` | 1482 | Medium |
| 6.3 | `SoilWaterRetentionCurveMod.F90` | `biogeophys/swrc_base.jl` | ~150 | Abstract type |
| 6.4 | `SoilWaterRetentionCurveClappHornberg1978Mod.F90` | `biogeophys/swrc_clapp_hornberg.jl` | ~200 | Per-column |
| 6.5 | `SoilWaterRetentionCurveVanGenuchten1980Mod.F90` | `biogeophys/swrc_van_genuchten.jl` | ~200 | Per-column |
| 6.6 | `SoilWaterMovementMod.F90` | `biogeophys/soil_water_movement.jl` | 2210 | **HIGH** |
| 6.7 | `SoilHydrologyMod.F90` | `biogeophys/soil_hydrology.jl` | 2910 | **HIGH** |
| 6.8 | `SnowHydrologyMod.F90` | `biogeophys/snow_hydrology.jl` | 4070 | **HIGH** — layer dynamics |
| 6.9 | `HydrologyNoDrainageMod.F90` | `biogeophys/hydrology_no_drainage.jl` | ~800 | Driver |
| 6.10 | `HydrologyDrainageMod.F90` | `biogeophys/hydrology_drainage.jl` | ~700 | Driver |
| 6.11 | `LakeHydrologyMod.F90` | `biogeophys/lake_hydrology.jl` | ~800 | Medium |
| 6.12 | `HillslopeHydrologyMod.F90` | `biogeophys/hillslope_hydrology.jl` | 1148 | Medium |

**Critical note for 6.8 (SnowHydrologyMod):** This is 4070 lines and contains the AD-hostile snow layer merge/split logic. For Phase 1, port it faithfully (including the discontinuities). The smoothing happens in Phase 3. But DO convert the dynamic layer indexing to fixed-max-depth masked arrays now.

**Acceptance:** Column-level regression tests. Run a single column for 1 year, compare all state variables to Fortran at every timestep.

---

### Tier 7: Canopy & Flux Physics (8 iterations, ~3 days)

| # | Fortran Source | Julia Target | Lines | GPU? |
|---|---------------|-------------|-------|------|
| 7.1 | `PhotosynthesisMod.F90` | `biogeophys/photosynthesis.jl` | 5209 | **HIGH** — per-patch, compute-heavy |
| 7.2 | `CanopyFluxesMod.F90` | `biogeophys/canopy_fluxes.jl` | 1756 | Yes |
| 7.3 | `CanopyHydrologyMod.F90` | `biogeophys/canopy_hydrology.jl` | 1191 | Yes |
| 7.4 | `BareGroundFluxesMod.F90` | `biogeophys/bareground_fluxes.jl` | ~500 | Yes |
| 7.5 | `SoilFluxesMod.F90` | `biogeophys/soil_fluxes.jl` | ~400 | Yes |
| 7.6 | `UrbanFluxesMod.F90` | `biogeophys/urban_fluxes.jl` | 1143 | Yes |
| 7.7 | `LunaMod.F90` | `biogeophys/luna.jl` | 1410 | Per-patch |
| 7.8 | `OzoneMod.F90` | `biogeophys/ozone.jl` | ~400 | Per-patch |

**Note:** 7.1 (Photosynthesis) is the single most compute-intensive per-patch routine and the biggest GPU payoff target. It's also 5209 lines. Budget 2-3 iterations.

---

### Tier 8: Biogeochemistry Physics (20 iterations, ~7 days)

| # | Fortran Source | Julia Target | Lines |
|---|---------------|-------------|-------|
| 8.1 | `CNPhenologyMod.F90` | `biogeochem/phenology.jl` | 4517 |
| 8.2 | `CNAllocationMod.F90` | `biogeochem/allocation.jl` | ~800 |
| 8.3 | `CNGRespMod.F90` | `biogeochem/growth_resp.jl` | ~400 |
| 8.4 | `CNMRespMod.F90` | `biogeochem/maint_resp.jl` | ~500 |
| 8.5 | `CNFUNMod.F90` | `biogeochem/fun.jl` | 1800 |
| 8.6 | `CNGapMortalityMod.F90` | `biogeochem/gap_mortality.jl` | ~400 |
| 8.7 | `SoilBiogeochemDecompCascadeBGCMod.F90` | `biogeochem/decomp_bgc.jl` | ~800 |
| 8.8 | `SoilBiogeochemDecompCascadeMIMICSMod.F90` | `biogeochem/decomp_mimics.jl` | 1367 |
| 8.9 | `SoilBiogeochemDecompMod.F90` | `biogeochem/decomp.jl` | ~600 |
| 8.10 | `SoilBiogeochemNitrifDenitrifMod.F90` | `biogeochem/nitrif_denitrif.jl` | ~600 |
| 8.11 | `SoilBiogeochemNLeachingMod.F90` | `biogeochem/n_leaching.jl` | ~300 |
| 8.12 | `CNNDynamicsMod.F90` | `biogeochem/n_dynamics.jl` | ~500 |
| 8.13 | `CNCStateUpdate1Mod.F90` | `biogeochem/c_state_update1.jl` | ~400 |
| 8.14 | `CNCStateUpdate2Mod.F90` | `biogeochem/c_state_update2.jl` | ~400 |
| 8.15 | `CNCStateUpdate3Mod.F90` | `biogeochem/c_state_update3.jl` | ~400 |
| 8.16 | `CNNStateUpdate1Mod.F90` | `biogeochem/n_state_update1.jl` | ~400 |
| 8.17 | `CNNStateUpdate2Mod.F90` | `biogeochem/n_state_update2.jl` | ~400 |
| 8.18 | `CNNStateUpdate3Mod.F90` | `biogeochem/n_state_update3.jl` | ~400 |
| 8.19 | `CNDriverMod.F90` | `biogeochem/cn_driver.jl` | 1337 |
| 8.20 | `CNVegetationFacade.F90` | `biogeochem/vegetation_facade.jl` | 1665 |

---

### Tier 9: Specialized Subsystems (12 iterations, ~4 days)

| # | Fortran Source | Julia Target | Lines |
|---|---------------|-------------|-------|
| 9.1 | `CNFireBaseMod.F90` | `biogeochem/fire_base.jl` | 1339 |
| 9.2 | `CNFireLi2014Mod.F90` | `biogeochem/fire_li2014.jl` | 1500 |
| 9.3 | `ch4Mod.F90` | `biogeochem/methane.jl` | 4450 |
| 9.4 | `VOCEmissionMod.F90` | `biogeochem/voc_emission.jl` | 1086 |
| 9.5 | `DustEmisBase.F90` + variants | `biogeochem/dust_emission.jl` | ~600 |
| 9.6 | `IrrigationMod.F90` | `biogeophys/irrigation.jl` | 1786 |
| 9.7 | `CNDVDriverMod.F90` | `biogeochem/cndv_driver.jl` | ~500 |
| 9.8 | `NutrientCompetitionFlexibleCNMod.F90` | `biogeochem/nutrient_competition.jl` | 1685 |
| 9.9 | `SatellitePhenologyMod.F90` | `biogeochem/satellite_phenology.jl` | ~400 |
| 9.10 | `CNCIsoFluxMod.F90` | `biogeochem/c_iso_flux.jl` | 1936 |
| 9.11 | `CNBalanceCheckMod.F90` | `biogeochem/cn_balance_check.jl` | ~600 |
| 9.12 | `BalanceCheckMod.F90` | `biogeophys/balance_check.jl` | 1053 |

---

### Tier 10: Infrastructure & Coupling (8 iterations, ~3 days)

| # | Fortran Source | Julia Target | Lines |
|---|---------------|-------------|-------|
| 10.1 | `atm2lndType.F90` | `types/atm2lnd.jl` | 1149 |
| 10.2 | `lnd2atmType.F90` | `types/lnd2atm.jl` | ~500 |
| 10.3 | `TopoMod.F90` | `infrastructure/topo.jl` | ~400 |
| 10.4 | `subgridAveMod.F90` | `infrastructure/subgrid_ave.jl` | 1098 |
| 10.5 | `accumulMod.F90` | `infrastructure/accumul.jl` | ~400 |
| 10.6 | `pftconMod.F90` | `constants/pftcon.jl` | 1602 |
| 10.7 | `controlMod.F90` | `infrastructure/control.jl` | 1350 |
| 10.8 | `clm_instMod.F90` | `infrastructure/instances.jl` | ~800 |

---

### Tier 11: Driver & Integration (3 iterations, ~2 days)

| # | Fortran Source | Julia Target | Lines |
|---|---------------|-------------|-------|
| 11.1 | `clm_driver.F90` | `driver/clm_driver.jl` | 1711 |
| 11.2 | — | `driver/timestep.jl` | new | Single-timestep entry point |
| 11.3 | — | `driver/run.jl` | new | Time loop, I/O orchestration |

**11.1 is the capstone.** It wires all Tier 5-9 modules together in the correct physics ordering. This is where the mask-based GPU dispatch replaces the OpenMP clump loop:

```julia
function clm_timestep!(state, forcing, params, backend)
    # All columns processed in parallel via GPU kernels
    surface_radiation!(state, forcing, params, backend)
    canopy_fluxes!(state, forcing, params, backend)
    soil_temperature!(state, params, backend)  # batched tridiagonal
    hydrology_no_drainage!(state, params, backend)
    # ... etc, same order as Fortran driver
end
```

---

## Validation Strategy

### Per-Module (runs in Ralph loop)
- Generate reference data: run Fortran CLM single-column for 1 year, dump all state/flux fields at every timestep to NetCDF
- Julia test loads reference, runs same forcing, compares field-by-field
- Pass criterion: max |error| < 1e-12 (double precision roundoff)

### Integration (manual checkpoints)
- After Tier 6: run 10-year single-point, compare annual means
- After Tier 8: run 1-year global 2°, compare spatial patterns
- After Tier 11: run standard CESM test suite (ERS, ERI tests)

### GPU Validation
- CPU and GPU backends must produce identical results (bitwise)
- Use KernelAbstractions CPU backend as reference

---

## Size Estimates

| Tier | Modules | ~Lines of Fortran | ~Iterations | ~Days |
|------|---------|-------------------|-------------|-------|
| 0 | 1 | 0 | 1 | 0.5 |
| 1 | 8 | 1,865 | 8 | 2 |
| 2 | 6 | 1,504 | 6 | 2 |
| 3 | 18 | ~12,000 | 18 | 5 |
| 4 | 12 | ~22,000 | 14 | 4 |
| 5 | 10 | ~9,000 | 10 | 4 |
| 6 | 12 | ~18,000 | 14 | 5 |
| 7 | 8 | ~12,500 | 10 | 3 |
| 8 | 20 | ~15,000 | 20 | 7 |
| 9 | 12 | ~15,000 | 12 | 4 |
| 10 | 8 | ~7,700 | 8 | 3 |
| 11 | 3 | ~3,000 | 3 | 2 |
| **Total** | **118** | **~117,000** | **~124** | **~41** |

---

## Ralph Loop Agent Configuration

### Per-Iteration Prompt Template
```
You are porting CLM from Fortran to Julia.

Current task: Port {FORTRAN_FILE} to {JULIA_FILE}

Context:
- Read the Fortran source at {FORTRAN_PATH}
- The Julia package is at CLM.jl/
- Already-ported dependencies: {LIST_OF_COMPLETED_MODULES}

Rules:
1. Translate ALL subroutines and functions. Do not skip any.
2. Preserve variable names exactly (for traceability to Fortran).
3. Replace Fortran `pointer` arrays with Julia typed fields.
4. Replace filter-based loops with mask-based iteration.
5. Use `@kernel` from KernelAbstractions for any loop over columns/patches.
6. Do NOT smooth discontinuities (min/max, phase change thresholds).
   Keep them identical to Fortran for now.
7. Write tests in test/{module_name}_test.jl that compare against
   reference data in test/reference_data/.
8. Run tests. Fix until green. Commit.

Do not modify any previously ported module.
Do not port modules not listed in the current task.
```

### Exit Conditions
- All tests pass
- `using CLM` loads without error
- One clean commit with message: "Port {module_name} from Fortran to Julia"

### Failure Recovery
- If tests fail after 3 attempts, commit what exists with `[WIP]` prefix
- Log the failure reason in `CLM.jl/PORTING_LOG.md`
- Move to next module (may unblock later)

---

## Phase 2 Preview: FATES Port (~52K lines, separate PRD)

After core CLM works in Julia:
- Port `src/fates/` following same tier strategy
- FATES has its own type system (cohorts, patches) that maps well to Julia
- `FatesHistoryInterfaceMod.F90` (9332 lines) is the single largest file — needs special handling

## Phase 3 Preview: AD Smoothing (separate PRD, requires scientist)

For each discontinuity class:
- **Phase changes:** Replace `if (t > tfrz) t = tfrz` with sigmoid relaxation
- **Min/max clamps:** Replace `min(x, 1.0)` with `softclamp(x, 1.0)`
- **Snow layer dynamics:** Replace discrete merge/split with continuous blending
- **Filter membership:** Already handled by mask conversion in Phase 1
- Validate that gradients are physically meaningful (finite, correct sign, reasonable magnitude)
- This is research, not engineering. Budget accordingly.
