# Subroutine-level port audit (Phase A2)

Fortran-vs-Julia subroutine coverage across the **138 ported physics modules**
(biogeophys + biogeochem + soilbiogeochem), via 16 parallel read-only agents.
Complements the module-level audit in the chat / `COMPLETION_ROADMAP.md`.

## Headline
| Status | Count |
|---|---|
| **full** (all substantive routines ported) | **102** |
| partial (some routines missing) | 23 |
| stub | 13 |

But most "stub" are **pure type-definition modules** (no physics — `LakeStateType`,
`SolarAbsorbedType`, `SurfaceAlbedoType`, `CNVegStateType`, `SoilBiogeochem*Type`,
`Species*Type`, `SoilStateType`, `SoilHydrologyType`) — correctly have nothing to
port beyond the struct, which exists in `types/`. So the real signal is in the
*partials with substantive missing routines*, triaged below.

## A. Genuine substantive physics gaps  → real Phase-A tasks
- **`SoilHydrologyMod`** (9/16) — the biggest. Missing **lateral flow** (`PerchedLateralFlow`, `SubsurfaceLateralFlow`), **infiltration-excess routing** (`RouteInfiltrationExcess`, `SetFloodc`), and **groundwater irrigation withdrawal** (`WithdrawGroundwaterIrrigation`). Core water-table/drainage routines ARE ported; lateral redistribution + irrigation aren't.
- **`GlacierSurfaceMassBalanceMod`** (1/3) — missing **`HandleIceMelt`**, **`ComputeSurfaceMassBalance`** (only runoff-term adjustment ported). Real glacier SMB gap.
- **`SnowHydrologyMod`** (13/14) — missing **`SnowCapping`** (excess-snow → ice/runoff capping).
- **`SoilStateInitTimeConstMod`** (3/7) — missing dust emission **threshold functions** (`ThresholdSoilMoistZender2003/Kok2014`) + the main time-const init.
- **`SoilTemperatureMod`** (11/14) — missing **`BuildingHAC`** = the prognostic urban building energy balance = **already Phase A1** (`UrbBuildTempOleson2015`). Confirms that gap.
- **`SurfaceWaterMod`** (3/6) — missing h2osfc-fraction state updates (`UpdateFracH2oSfc`, `BulkDiag_FracH2oSfc`, `UpdateState_TooSmallH2osfcToSoil`). **VERIFY** — may be partly folded into the hydrology driver.

## B. Alternative methods (not gaps — the port chose one)
- **`SoilWaterMovementMod`** (7/10) — missing `soilwater_zengdecker2009`, `soilwater_moisture_form` (alternative Richards-eq integration forms; one form is ported).
- **`SoilHydrologyInitTimeConstMod`** (2/5) — missing VIC-specific init (`initSoilParVIC`, `initCLMVICMap`) — VIC is an alternate hydrology option.

## C. Parameter readers / init (minor — params loaded elsewhere in the port)
`CNAllocationMod`, `CNFUNMod`, `SoilBiogeochemNLeachingMod`, `InfiltrationExcessRunoffMod` (`readParams`); `CNFireBaseMod` (`CNFireInit`); `SnowSnicarMod` (`SnowOptics_init`/`SnowAge_init` — optics/aging file loaders; the RT physics is ported). Low priority — the Julia port loads these params via its own `read_params.jl` path.

## D. Dynamic-subgrid coupling (= Phase C, dyn_subgrid)
- **`CNVegetationFacade`** (4/7) — missing `UpdateSubgridWeights`, `DynamicAreaConservation`, `InterpFileInputs` — the transient land-use hooks (deferred subsystem).
- **`CNDVDriverMod`** (1/3) — missing `CNDVHist` (history), `BuildNatVegFilter`.

## E. Likely false-negatives (ARE ported, agent missed structural renaming) — verify
- **`UrbanParamsType`** (2/3, "missing `UrbanInput`") — **ported this session** as `src/infrastructure/urban_input.jl` (`read_urban_input!`).
- **`SaturatedExcessRunoffMod`** ("kernel-only") — likely ported as a kernel.
- **`MEGANFactorsMod`** (`megan_factors_get`), **`CNSpeciesMod`** (`species_from_string`), **`TemperatureType`** (`UpdateAccVars_CropGDDs`), **`TotalWaterAndHeatMod`** (`AccumulateLiquidWaterHeat`) — small utilities likely inlined in callers.

## Verification + remediation (Phase A, this round)
Parallel agents verified the flagged partials and ported the genuine gaps. **Key correction: the subroutine audit over-counted gaps** — it checked the Julia files mapped via "Ported from" headers but missed sibling files, producing *file-mapping false-negatives*.

**Were already ported (false-negatives — added missing test coverage):**
- **Soil lateral flow + irrigation** — all 5 routines (`PerchedLateralFlow`, `SubsurfaceLateralFlow`, `RouteInfiltrationExcess`, `SetFloodc`, `WithdrawGroundwaterIrrigation`) live in `soil_hydrology.jl` (the audit checked `soil_hydrology_mod.jl`). Added conservation tests (36).
- **Snow capping** — full chain in `snow_hydrology.jl`, wired in the driver. Added conservation tests (26).
- **h2osfc fraction updates, saturated-excess runoff, liquid-water-heat accumulation** — all present (`surface_water.jl`, `sat_excess_runoff.jl`, `total_water_heat.jl`).

**Genuinely missing → PORTED this round (real ports, FD/conservation-tested):**
- **Prognostic urban building temperature** (`UrbBuildTempOleson2015` + `BuildingHAC`) → `urb_build_temp_oleson2015.jl` (20 tests: energy-balance closure + ForwardDiff check).
- **Glacier surface mass balance** (`HandleIceMelt`, `ComputeSurfaceMassBalance`, real `AdjustRunoffTerms`) → `glacier_surface_mass_balance.jl` (27 tests).
- **Dust soil-moisture thresholds** (`ThresholdSoilMoist{Zender2003,Kok2014}`) + wired `gwc_thr_col`/`mss_frc_cly_vld_col` (were NaN).

**Genuinely missing → still TODO (small):** `megan_factors_get` (MEGAN compound lookup; VOC emission not wired into the driver), `species_from_string` (C12/C13/C14/N parser; numeric IDs hardcoded), `UpdateAccVars_CropGDDs` (crop-GDD accumulation is a stub, blocked on the unported `accumulMod`).

## Bottom line
The subroutine audit's 102/138 "full" was **pessimistic** — most flagged "partials" were file-mapping false-negatives (already ported elsewhere). After this round's three real ports (urban prognostic building, glacier SMB, dust thresholds), the remaining single-point CLM5 physics gaps are a **very short list**: MEGAN/VOC driver wiring, `species_from_string`, and crop-GDD accumulation (the last gated on porting `accumulMod`). Everything else outstanding is alternative methods (Phase E), I/O (Phase B), transient subgrid (Phase C), AD/GPU (Phase D), or FATES (Phase F).
