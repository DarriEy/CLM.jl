# Fortran-dependent validation backlog

**What this is.** The list of CLM.jl subsystems whose Fortran parity is *not yet
established* because it requires **generating new Fortran (CTSM) ground truth** — a new
site, a non-soil landunit, or a config-flag rerun — that does not exist on disk today.

**When to do it.** The CTSM source is re-clonable and buildable on any box
(`git clone -b ctsm5.3.012 https://github.com/ESCOMP/CTSM.git`, commit
`ab466d6f9789ca3df2c72bda46cf7afed2d04102`), so this is **anytime** work — it is *not*
gated by any particular machine. (Contrast: Metal-GPU validation, which is gated by the
Apple-Silicon dev box and must happen before it is returned.)

**What is deliberately NOT here** (doable anytime, needs no new Fortran run):
- Fire `NFIRE`/`FAREA_BURNED` — blocked only on porting accumulator-restart I/O; scores
  against the existing `bgc_ref_firech4` dump. (In progress.)
- Diffing the **existing-but-unread** mid-CN-driver dumps `after_competition` /
  `after_ecosysdyn_predrain` (ground truth already on disk, no harness reads it).
- Albedo bands (`albd/albi/fabd/fabi`), LUNA-recomputed `vcmx25_z/jmx25_z` (harness
  currently *injects* it), aerosol/dust/ozone — all present in existing dumps, just not
  diffed. Pure harness work.

**Run-verified harness status (2026-07-17):** see
[`docs/PARITY_COVERAGE_2026-07.md`](PARITY_COVERAGE_2026-07.md) for the refreshed
scorecard — every `scripts/fortran_parity_*.jl` added since the 2026-06-17 audit was
executed against the surviving local dumps. Headlines: active-layer / soilresis /
albedo (no-injection) are now **VALIDATED**; aerosol / daylength are green-by-skip
(need a snowy / `dayl` dump); **LUNA's EOD vcmax update DIVERGES ~5.3–5.6% high** (a
documented Rubisco-N optimum residual, jmax matches exactly); lake was validated
historically but its dump domain is gone; isotopes still need a c13/c14 spinup.

Cross-references: `docs/HARNESS_COVERAGE_AUDIT.md` (per-module coverage, tables A/B/C),
`docs/CH4_FIRE_PARITY.md` (§7 "not validated"), `docs/N_CYCLE_PARITY.md`,
`docs/PARITY_STATUS.md`. Ground-truth recipe: the `restFile_write_dump` /
`bgcdumpMod.F90` SourceMods instrumentation under `scripts/validation/`, run an
instrumented `cesm.exe` **plainly** (never under lldb — macOS xzone `brk` guards) with
`BGCDUMP_NSTEP_LO/HI` + `PDUMP_NSTEP_LO/HI` set to the diff window.

---

## A. Needs a new SITE (surfdata + forcing)

| # | Item | Julia module(s) | What a Fortran run must exercise | Candidate site (migrated) |
|---|------|-----------------|----------------------------------|---------------------------|
| A1 | **Methane as a SOURCE** (transport / aerenchyma / ebullition / `finundated>0`) — the regime CH4 actually matters for | `methane.jl` | `use_lch4=.true.` at a site with `frac_h2osfc>0` so `finundated>0`; diff `CH4_EBUL_*`, `CH4_AERE_*`, `CH4_SURF_*`, `CONC_CH4/O2_*`, `TOTCOLCH4`. At Bow these are vacuous (`finundated≡0`). **STILL BLOCKED 2026-07-17 (round 2): MerBleue assets RESTORED and a full Fortran CH4 reference was GENERATED (`clm_bgc_spinup/merbleue_ref_ch4`, 169-step window), but the site does NOT inundate — the water table equilibrates ~1.9 m deep even in mid-July, `H2OSFC≡0`, so `finundated≡0` across the whole window (same degenerate regime as Bow). `CH4_EBUL_*` and `CH4_DFSAT_FLUX` are identically zero (vacuous). CH4-as-a-source is a physical-hydrology blocker, not an assets blocker.** The MerBleue restart is also SP-only (no BGC restart exists), so CN+CH4 are `use_init_interp` cold-started ⇒ a production/O2-demand confound. See docs/CH4_FIRE_PARITY.md §9. | `domain_Peatland_MerBleue_Canada` (present, **does not inundate**) |
| A2 | **Fire crop/peat/tropical/land-use branches** | `fire_li2016.jl`, `fire_base.jl` | `baf_crop` (crop CFT), `baf_peatf` (peatf>0), `dtrotr`/`TROTR` (tropical tree), `lfc` (land-use transition) — all `≡0` at Bow | Mead (crop), MerBleue (peat), Aripuana/Leticia (tropical) |
| A3 | **Crops** (`crop.jl` + crop phenology/allocation) + **crop N** (`n_fert!`, `n_soyfix!` — currently UNWIRED, correctly refused blind-wiring per #218) | `crop.jl`, crop phen/alloc, `n_dynamics.jl` | `use_crop=.true.` + surfdata with a crop CFT + crop spinup; diff crop C/N pools, GDD, harvest, fertilizer | `domain_Cropland_Mead_USA` |
| A4 | **Irrigation** | `irrigation.jl` | crop-bearing site + `irrigate=.true.`; at Bow `n_irrig_steps_left==0`, `irrig_rate` NaN (no crop CFT). Driver call sites were empty stubs — verify wired. | `domain_Cropland_Mead_USA` |

## B. Needs a non-soil LANDUNIT run

| # | Item | Julia module(s) | Fortran run | Notes |
|---|------|-----------------|-------------|-------|
| B1 | **Lake** | `lake_*.jl` (con/fluxes/hydrology/temperature) | lake-landunit column | `fortran_parity_lake.jl` was **VALIDATED historically** (thermodynamics + surface turbulent-flux solve, through 3 real fixes ending commit `3db0c1a` "closes the surface-flux residual"). **The `clm_lake_run` h0 domain is now ABSENT on disk (2026-07-17)** → needs regenerating to reproduce, but is not un-validated. See `docs/PARITY_COVERAGE_2026-07.md`. |
| B2 | **Urban** | `urban_*.jl` (albedo/fluxes/radiation) | urban-landunit column | Metal-validated; `urban_radiation` is called but no-ops in Bow |
| B3 | **Glacier (glcmec)** | glacier path | glcmec-landunit column | only surfdata/init coverage today |
| B4 | **Hillslope** | `hillslope_hydrology.jl` | hillslope routing config | unit-tested only |

## C. Needs a config-flag rerun

| # | Item | Julia module(s) | Fortran run | Notes |
|---|------|-----------------|-------------|-------|
| C1 | **Isotopes** (C13/C14) | `carbon_isotopes.jl`, `c_iso_flux.jl` | `use_c13/use_c14=.true.` rerun | `fortran_parity_isotopes.jl` RAN (2026-07-17): confirmed the Fortran `use_c13/use_c14=.true.` branch **crashes at restart read** (the `tier1_iso` branch restart carries no c13/c14 vars) → **still needs a c13/c14-enabled spinup**. The irrigate=.true. half runs clean but is inactive at Bow (no crop CFT). See `docs/PARITY_COVERAGE_2026-07.md`. |
| C1a | **LUNA vcmax EOD update** (Rubisco-N optimum) | `luna.jl` `update_photosynthesis_capacity!` / `nitrogen_allocation!` | Fortran `Allocation` **internal** instrumentation (emit `Ncb`/`vcmx25_opt`, not just post `vcmx25_z`) | **KNOWN DIVERGE (2026-07-17):** `fortran_parity_luna_update.jl` shows Julia's post-EOD `vcmx25_z` runs **+5.2–5.6% high** (jmx25_z matches exactly). Constants byte-identical to `LunaMod.F90`; the residual is in the joint N-allocation optimizer's carboxylation pool `Ncb`. Not tuned. See `docs/PARITY_COVERAGE_2026-07.md` §"LUNA vcmax residual". |
| C2 | **CNDV** (dynamic vegetation) | `cndv.jl` | `use_cndv` rerun | Metal-validated only |
| C3 | **VOC / MEGAN** | `voc_emission.jl` | non-empty MEGAN compound list | not called in Bow (empty compound list) |
| C4 | **MIMICS decomposition** | `decomp_mimics.jl` | `decomp_method='MIMICS'` rerun | Bow uses CENTURY/BGC |

## D. FATES deeper validation

| # | Item | What a Fortran run must exercise |
|---|------|----------------------------------|
| D1 | Cold-start (#213, 27/27) + time-stepped dynamics (#217, phase-attributed) are DONE. Remaining: a **multi-year / disturbance** FATES-enabled Fortran run, and validating the **boreal cold-start fix (#227)** and the **fixed-biogeog screen (#197)** against Fortran (not just internal). Needs FATES-enabled CTSM reference runs (the `fates_parity_1pt` / `fates_parity_ne3` cases are the scaffolding; ground truth `fates_pdump_fortran.txt.gz` for the current scenario is committed). |

---

## Notes on prioritization

- **A1 (methane-as-source at MerBleue)** is the single highest-value scientific gap — it
  validates the half of the CH4 model that is the point of the CH4 model.
- **A2/A3/A4 (crop + fire crop/peat)** can share Fortran runs (a crop-bearing and a
  peatland run exercise several rows at once).
- **B (landunits)** are robustness/coverage rather than a specific science claim; batch
  them if a multi-landunit case is convenient.
- Several rows have an **existing `fortran_parity_*.jl` harness** (lake, isotopes) — the
  first action for those is to *run the existing harness* and see whether a dump already
  exists, before generating new ground truth.
