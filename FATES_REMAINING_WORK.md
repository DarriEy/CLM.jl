# FATES — Remaining Work Scoping (post Tier-F module port)

**Status as of 2026-06-23:** all 19 FATES batches (0–18) are ported — 62 modules in
`src/fates/`, ~77k Fortran lines, full suite 20,642 green (PRs #62–#80). Every module is
**standalone**: none is referenced from `src/driver/`, none is in `CLMInstances`. The FATES
*physics* is done; what remains is **integration + a small set of deferred fill routines**,
not new science.

This document scopes the three remaining work areas, derived from a read-only audit of both
the Julia port and the Fortran sources (`installs/clm/src/fates` + `src/utils/clmfates_interfaceMod.F90`).

---

## ✅ Progress (updated 2026-06-23) — FATES now RUNS in the live driver

Done since this doc was written (suite 20,642 → 20,772; PRs #81–#84):
- **W1+W2 (PR #81):** `CLMInstances.fates::Union{AbstractFatesInterface,Nothing}` ownership model (skipped from AD dual-copy + GPU adapt, like `surfdata`); `clm_fates_init!` cold-starts a finite carbon-only single FATES site. Proven by test.
- **W3+W4 (PR #84):** all four per-timestep hooks wired & gated behind `use_fates` — `FatesSunShadeFracs` (sun/shade), `FatesNormalizedCanopyRadiation` (canopy albedo), `btran_ed!` (water stress), `FatesPlantRespPhotosynthDrive` (photosynthesis); bc_in packed from CLM column state, bc_out unpacked into canopystate/surfalb/photosyns/energyflux. Init wiring (`clm_fates_init!` into `clm_initialize.jl`, gated). **FATES drives radiation/btran/photosynthesis through `clm_driver` for a `use_fates` column; default path byte-identical.**
- **Restart fills R1–R5 (PR #82):** the full deferred restart pack/unpack — cohort diagnostics/mortality, patch fuel/albedo/litter, site demographic/mortality/damage arrays, `update_3dpatch_radiation!`. Demographic round-trip complete.
- **Plant-hydraulics Tier A (PR #83):** 13 init/lifecycle routines ported; all `hlm_use_planthydro`-gated stub call-sites closed; the no-op `AccumulateMortalityWaterStorage` shadow removed.

**Update 2 (PR #85) — real parameters + the daily demographic step:**
- **Real FATES param-file reader** — `read_fates_params!` (CDL parser → the port's existing register/receive infra; `data/fates/fates_params_default.cdl` in-repo); 297/297 params, **numpft now 14** (real). `clm_fates_init!` uses it; synthetic `_fates_spike_setup_pft!` superseded.
- **W5 daily-dynamics hook** — `ed_ecosystem_dynamics` + `ed_update_site` + `TotalBalanceCheck` gated into the daily driver branch; **runs mass-conserving (|err|≤1e-5) on the real 14-PFT params.** FATES now advances its cohort/patch population on real parameters. Suite → 20,969.

**Remaining tail (the table at the bottom is updated with ✅ markers):**
1. **W4b full *in-solve* photosynthesis coupling** — currently `FatesPlantRespPhotosynthDrive` is called *adjacent to* `canopy_fluxes_core!` (post-solve `t_veg`), not inside its Enzyme-compilable iterative solve. Future: move inside for two-way leaf-temp↔flux coupling.
2. **Plant-hydraulics Tier B** — `Hydraulics_BC` transpiration solve (~544 F90 lines) + `hydraulics_drive`/`FillDrainRhizShells`/`RecruitWUptake` + a NEW gated `hydraulics_drive` callsite in the canopy-flux path. This is what fills `ftc_*`/`btran`/`leaf_psi` (consumer branches already ported, reading zeros until then).
3. **History fills H1–H6** — NOW UNBLOCKED (per-step state is live under `use_fates`): `update_history_dyn1!` finish, `dyn2`, `nutrflux`, `hifrq1/2`, `hydraulics`. Plus the size/age/coage/height class-index helpers shared by `dyn2`/`hifrq2`.
4. **Restart R6/R7/R8** — disturbance running-means (needs SetRMean/GetRMean restart accessors); cohort PRT pools (needs `InitPRTObject` on restart cohorts); cohort hydraulics (needs `InitHydrCohort`, now ported via Tier A).
5. ✅ ~~Real FATES param-file reader~~ — DONE (PR #85): `read_fates_params!` from `data/fates/fates_params_default.cdl`.
6. **Multi-veg-patch + multi-site** — current column↔patch map is `p = col.patchi[c]+1` (single veg patch MVP); the full `ifp` walk is needed for real FATES columns. AD/GPU for FATES columns remains explicitly out of scope.
7. **Multi-day driver-integrated run / stability** — W5 is currently *helper-tested* (drives `fates_daily_dynamics_step!` directly); a full multi-day `clm_drv!` run with `use_fates=true` exercising the whole loop (radiation→photosynthesis→daily-dynamics→history) end-to-end, plus a spin-up/stability check, is the natural next validation.

---

## TL;DR — the dependency ordering that drives sequencing

The three areas are **not independent**; one unblocks the others:

```
            ┌─────────────────────────────────────────────┐
            │  LIVE-DRIVER WIRING  (the keystone, W1–W6)   │
            │  makes FATES actually run in clm_driver       │
            └───────────────┬───────────────┬──────────────┘
                            │ unblocks      │ unblocks
                            ▼               ▼
        ┌───────────────────────────┐   ┌──────────────────────────────┐
        │ HISTORY buffer-fills (H1–H6)│  │ Plant-hydraulics Tier B (BC   │
        │ per-timestep state is NaN   │  │ solve) — needs a new gated    │
        │ until the driver populates  │  │ hydraulics_drive callsite     │
        │ it; do AFTER W4/W5          │  │ in the canopy-flux path       │
        └───────────────────────────┘   └──────────────────────────────┘

   Standalone NOW (no driver dependency, do anytime, round-trip / unit testable):
     • RESTART fills R1,R2,R3,R4,R5,R6   • Plant-hydraulics Tier A (init/lifecycle)
   Blocked on cohort object init (InitPRTObject / InitHydrCohort on restart cohorts):
     • RESTART fills R7 (PRT pools), R8 (hydro)
```

**Recommended order:** (1) **Live-driver wiring W1–W4** to get a runnable single-site
carbon-only FATES — this is the highest-value milestone and validates the whole port
end-to-end. (2) Interleave the **standalone restart fills** and **plant-hydraulics Tier A**
(independent, cheap, can run in parallel waves). (3) **W5–W6 + history fills** once per-step
state is live. (4) **Plant-hydraulics Tier B + PRT/hydro restart fills** last, when
`hlm_use_planthydro` / restart-parity are actually wanted.

**Rough total:** ~15–18 agent-batches. Live-driver wiring ≈6, B18-followup fills ≈8 (split
by dimension), plant-hydraulics ≈3–4. None is blocked by a missing type — the Julia type
port is structurally complete.

---

## 1. Live-driver wiring (the keystone) — ~6 sub-batches

**Goal:** turn the standalone modules into a FATES that runs inside `clm_driver`.

**Good news:** the CLM.jl driver is already **FATES-aware scaffolding** — a
`FATESMode <: AbstractBGCMode`, plumbed `use_fates*` config flags, `col.is_fates`/`pch.is_fates`
tags, an `all_soil_patches` filter ("for FATES SP"), and **~7 named placeholder comment-sites
in `clm_drv_core!` at the exact lines where the Fortran host hooks live**. This is wiring into
existing hooks, not greenfield.

**Fortran spec:** host glue is `src/utils/clmfates_interfaceMod.F90` (3896 lines,
`type hlm_fates_interface_type`, global `clm_fates`). It maps FATES sites 1:1 onto `istsoil`
columns via `f2hmap`. One `fates(nc)` per clump.

### The one hard decision — CLMInstances ownership (linked-list vs SoA/AD)
`ed_site_type`/`fates_patch_type`/`fates_cohort_type` are heap-allocated, pointer-linked
mutable structs (`taller/shorter/older/younger`, abstract `prt`), fundamentally incompatible
with the SoA/Metal/Enzyme design of `CLMInstances` and its AD dual-copy loop
(`calibration.jl:326` blindly converts every field `Float64→Dual`).

**Recommended resolution:** attach FATES as a single
`fates::Union{fates_interface_type,Nothing}=nothing` field, treated exactly like `surfdata` —
**excluded from the dual-copy skip-list (`calibration.jl:332`) and from `@adapt_structure`
GPU movement.** Documented v1 boundary: **FATES columns run CPU-only and are not
AD-differentiable.** The ported `fates_interface_type{sites::Vector, bc_in/out::Vector,
bc_pconst}` (`FatesInterfaceMod.jl:83`) already has the right container shape.

### Insertion points (CLM.jl placeholders already present in `clm_drv_core!`)
| Hook | CLM.jl line | Fortran source |
|---|---|---|
| SP-LAI / InterpFileInputs | 615–651 | clm_drv top |
| `wrap_sunfrac` → `FatesSunShadeFracs` | **831** | clm_drv:608 |
| `wrap_btran` → `btran_ed!` | 896 | CanopyFluxesMod:865 |
| `wrap_photosynthesis` → `FatesPlantRespPhotosynthDrive` (inside `canopy_fluxes_core!`) | 932 | CanopyFluxesMod:1117 |
| BGC/litter coupling | 1450–1472 | dynamics_driv Part I/III |
| FATES-SP phenology | 1531 | SP path |
| **daily dynamics** (`UpdateFatesRmean`, hifrq hist, `if is_beg_curr_day`→`dynamics_driv`, `setFilters!`) | **1624–1631** | clm_drv:1143 |
| `wrap_canopy_radiation` → `FatesNormalizedCanopyRadiation` | ~1687 (in `surface_albedo!`) | SurfaceAlbedoMod:1048 |
| seed dispersal | 1703 | clm_drv:1277 |
| `TransferZ0mDisp` | in `biogeophys_pre_flux_calcs!` :851 | BiogeophysPreFluxCalcs:162 |

**Note:** btran/photosynthesis/accumulate run *inside* the CanopyFluxes iterative solve, not
as standalone driver phases — `canopy_fluxes_core!` needs an internal FATES branch.

### bc_in/bc_out mapping — CLM.jl already computes ~all carbon-only sources
Spot-checked OK: `t_soisno`, `smp_sl` (reuse the PHS smp_l path), `h2osoi_liqvol`, `solad/solai`,
`pbot`, `cair` (367ppm·pbot), `coszen`, `albgrd/albgri`, canopy-solve locals (esat/eair/rb/t_veg).
⚠️ only-if-needed: SPITFIRE fire-forcing streams (lightning/pop/RH24), SP-LAI streams.
bc_out sinks all exist (`canopystate.*`, `surfalb.*`, `photosyns.rssun/rssha`).

### MVP path (single-site, carbon-only, no hydro/fire/CNP/LUH)
- **W1** CLMInstances `fates` field + dual-copy/adapt skip + init globals/alloc + column→site map. *(High risk — the ownership decision; do as a spike first.)*
- **W2** Cold-start chain `init_site_vars!→set_site_properties!→init_patches!→init_cohorts!→ed_update_site` — finite single-site build.
- **W3** Radiation hooks (831 + ~1687); pack solad/coszen, unpack fsun/laisun/albd/fab*.
- **W4** BTRAN + photosynthesis hooks (896, 932) inside the canopy solve → GPP/respiration.
- **W5** Daily dynamics (1624–1631) + `setFilters!` + `TotalBalanceCheck` (fire/hydro gated off).
- **W6** History `update_history_dyn1!` wiring + single-site carbon-only validation harness vs Fortran.

**Top risks:** (1) the AD/GPU exclusion must be a conscious documented boundary; (2) patch
indexing `p = ifp + col%patchi(c)` + bare-ground patch + `setFilters!` rebuild must mirror
Fortran exactly or column averages break; (3) hooks-inside-CanopyFluxes; (4) delete the legacy
orphan `src/biogeochem/fates_interface.jl` (1313-line no-op stub, not `include`d) to avoid confusion.

---

## 2. B18-followup buffer-fill routines — ~8 sub-batches

All registries (472 history vars / 571 handles; 184 restart vars) are **100% complete** — only
the fill bodies are deferred (`# TODO Batch 18-followup`). No deferral is blocked by a missing
type/field; the Julia type port is structurally complete.

### 2a. RESTART fills (`set_restart_vectors!`/`get_restart_vectors!`) — ~138 vars / ~12k flat slots
The demographic skeleton round-trip + counts already work. Deferred = diagnostic/flux/PRT/hydro
copies, symmetric on pack+unpack.

**Standalone NOW (pure round-trip tests, fields already exist):**
- **R1 (S)** cohort diagnostic scalars (gpp/npp/resp_acc…) + mortality fluxes (bmort/hmort/…) — ~19 fields.
- **R2 (S)** patch fuel/scorch + radiation/albedo (`gnd_alb_dir/dif`) + site hydraulics scalars.
- **R3 (M)** `update_3dpatch_radiation!` (jl:970; F90 ~132) — recompute per-patch 3D albedo after restart; **solvers already ported** (`PatchNormanRadiation`/two-stream); only needs R2's gnd-albedo round-trip.
- **R4 (L)** patch litter blocks (~456 fields; nested element×pft×dcmpy×cwd×soillayer) — `litter_type` fully exists.
- **R5 (L)** site demographic/mortality/damage/flux-diag arrays (~11k of the 12k slots) — regular nested scpf/scag/cd/element loops; mechanical but bulky; split the (gated) tree-damage cross-tab into its own sub-batch.
- **R6 (M)** site + patch running means — needs `SetRMeanRestartVar`/`GetRMeanRestartVar` accessors ported first.

**Blocked on cohort object init** (`FatesRestartInterfaceMod.jl:821` — restart cohorts have
`prt=nothing`/`co_hydr=nothing`):
- **R7 (L)** cohort PRT pool pack/unpack + run `InitPRTObject` on restart cohorts (+ `InitPRTBoundaryConditions`/`UpdateCohortBioPhysRates`).
- **R8 (M)** cohort hydraulics + `InitHydrCohort` (planthydro-gated; depends on §3 Tier A).

### 2b. HISTORY fills — ~3,600 Fortran lines; **all blocked on live-driver wiring**
Standalone, the cohort `gpp_tstep`/`resp_*`/PRT/`co_hydr` and site `flux_diags` fields are
NaN/unset — nothing meaningful to write until the driver populates per-step state. Do AFTER W4/W5.
- **H1 (M)** `update_history_nutrflux!` (220) — CNP fluxes (CNP-mode only).
- **H2 (M)** `update_history_hifrq1!` (258) — site GPP/AR/NEP/resp.
- **H3 (L)** `update_history_hifrq2!` (448) — scpf/scls/cnlf disaggregation; needs the size/age/coage/height **class-index helpers**.
- **H4 (L)** `update_history_hydraulics!` (374) — planthydro-gated.
- **H5 (L)** finish `update_history_dyn1!` (~60 remaining fills of the 672-line source).
- **H6 (XL)** `update_history_dyn2!` (1853) — split by output dimension (pft / scpf / scls+scag / age / cnlf / cwd+fuel+landuse) into ~5 sub-batches; shares the class-index helpers with H3.

Cosmetic field nuances (not blockers): `ar`=`resp_m`+`resp_g`; fire mort=`crownfiremort`+`cambialfiremort`;
no `lmort_logging` (use `lmort_direct`); running means live at patch not cohort level.

---

## 3. Plant-hydraulics stubs — ~3–4 sub-batches

The hard math is **already ported** (1D Taylor solver core ~1,300 lines, all plant-level helpers,
WTF/WKF retention+conductivity, mem types). ~2,700 Fortran lines / 19 routines remain (the two
2D solvers `MatSolve2D`/`PicardSolve2D` are correctly `fates_endrun`-stubbed and **not** ported).
All call-sites are gated behind `hlm_use_planthydro` (default off) — inert today.

The field contract is already satisfied: `ed_cohort_hydr_type` exposes
`ftc_ag/ftc_troot/ftc_aroot/btran/is_newly_recruited/errh2o`, and the **consumer branches are
already ported** (`EDMortalityFunctionsMod.jl:136` reads `ftc_*`; `FatesPlantRespPhotosynthMod.jl:671/898`
reads `co_hydr.btran`/`leaf_psi`). They read zeros until the BC solver runs.

**Tier A — init/lifecycle (standalone, closes ~80% of stub call-sites):** `InitHydrCohort`(15),
`DeallocateHydrCohort`(14), `SavePreviousRhizVolumes`(13), `constrain_water_contents`(23),
`AccumulateMortalityWaterStorage`(36, replaces the misleading no-op shadow in
`EDLoggingMortalityMod.jl:67`), `RecruitWaterStorage`(53), `UpdateSizeDepRhizVolLenCon`(109)+
`UpdateSizeDepRhizHydProps`(26), `UpdateSizeDepPlantHydStates`(79), `FuseCohortHydraulics`(100),
`ConstrainRecruitNumber`(97), `InitHydrSites`(170), `HydrSiteColdStart`(186). ~860 lines.
Testable by flipping `hlm_use_planthydro=itrue` + constructing a hydro site, no transpiration solve.

**Tier B — the transpiration solve (needs a NEW gated driver callsite):** `SetMaxCondConnections`(86),
`Report1DError`(136), `RecruitWUptake`(83), `FillDrainRhizShells`(148),
**`Hydraulics_BC`(~544 — the core per-cohort uptake orchestrator that fills `ftc_*`/`btran`/`leaf_psi`)**,
`hydraulics_drive`(~27 dispatcher). **There is currently no `hydraulics_drive` callsite anywhere
in Julia** — it must be added into the FATES canopy-flux timestep path (Fortran calls it from the
photosynthesis path). Only meaningful end-to-end once §1 wiring exists.

**Tier C — I/O (deferrable):** `RestartHydrStates`(215) — restart parity only; gated hooks already
exist at `FatesRestartInterfaceMod.jl:496,821` (overlaps R8).

---

## Consolidated batch table

| # | Batch | Area | Effort | Depends on | Status |
|---|---|---|---|---|---|
| 1 | W1 CLMInstances + init globals/alloc + col→site map | wiring | M (High risk) | — | ✅ #81 |
| 2 | W2 cold-start chain → finite site | wiring | M | W1 | ✅ #81 |
| 3 | W3 radiation hooks (sunfrac/albedo) | wiring | M | W2 | ✅ #84 |
| 4 | W4 btran + photosynthesis hooks | wiring | M | W2 | ✅ #84 (W4b adjacent, not in-solve) |
| 5 | W5 daily dynamics + balance check (helper-tested) | wiring | L | W2–W4 | ✅ #85 (setFilters/multi-day deferred) |
| – | Real FATES param-file reader (CDL → register/receive) | params | L | — | ✅ #85 (numpft=14) |
| 6 | W6 history dyn1 wiring + validation harness | wiring | M | W5 | ⬜ (init-wiring done #84; dyn1 fills pending) |
| – | R1 cohort diag+mort restart fills | restart | S | — | ✅ #82 |
| – | R2 patch fuel/rad/site-hydro-scalar restart | restart | S | — | ✅ #82 |
| – | R3 `update_3dpatch_radiation!` | restart | M | R2 | ✅ #82 |
| – | R4 patch litter restart (~456) | restart | L | — | ✅ #82 |
| – | R5 site demog/mort/damage/flux arrays (~11k) | restart | L | — | ✅ #82 |
| – | R6 running-means restart (+accessors) | restart | M | — | ⬜ not started |
| – | R7 cohort PRT restart + InitPRTObject | restart | L | InitPRTObject | ⬜ not started |
| – | R8 cohort hydro restart + InitHydrCohort | restart | M | §3 Tier A (done) | ⬜ not started |
| – | H1 nutrflux / H2 hifrq1 / H3 hifrq2 / H4 hydraulics / H5 dyn1-finish / H6 dyn2(×5) | history | M–XL | W4 (done) → live state | ⬜ UNBLOCKED, not started |
| – | PH-A plant-hydraulics init/lifecycle | hydro | M (×1–1.5) | — | ✅ #83 |
| – | PH-B `Hydraulics_BC` transpiration solve + new driver callsite | hydro | L (×1.5–2) | §1 wiring (done) | ⬜ not started |
| – | PH-C `RestartHydrStates` | hydro | M | PH-A | overlaps R8 |

**Note on AD/GPU:** every FATES column is CPU-only / non-differentiable in this plan (the
`Union{…,Nothing}` ownership decision). Making FATES AD/GPU-capable is a separate, much larger
effort (would require an SoA re-expression of the cohort/patch hierarchy) and is explicitly out
of scope here.
