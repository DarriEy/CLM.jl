# CLM.jl — Full-Platform Fortran-Parity Coverage Plan

_Written 2026-06-17. Companion to `docs/HARNESS_COVERAGE_AUDIT.md`. This is the
plan to take CLM.jl from "biogeophysics + CN/BGC pools validated in one Bow soil
column" to "every active process module validated against Fortran ground truth"._

The reference Fortran run is the **Bow-at-Banff single column** instrumented
`cesm.exe` (case `cases/symfluence_build/`), dumping CLM5 restart-format
snapshots `pdump_<boundary>_n<nstep>.nc` via SourceMods `restFile_write_dump`.

## Ground-truth facts established by this audit (verified, not assumed)

The Bow `lnd_in` (verified from `CaseDocs/lnd_in`) actually has:
`use_lch4=.true.`, `fire_method='li2016crufrc'`, `use_nitrif_denitrif=.true.`,
`use_fun=.true.`, `use_luna=.true.`, `use_hydrstress=.true.`,
`soil_decomp_method='CENTURYKoven2013'`, **`use_crop=.false.`**,
`use_fates=.false.`, `use_c13=.false.`, `use_c14=.false.`, `irrigate=.false.`.

Important corrections to the audit's Table C assumptions:
- **Methane and fire are ON in the Fortran run** (not gated off). But the
  CLM5 **restart dump does not persist the ch4 prognostics** (verified: 0 ch4
  vars in `pdump_before_step_n1757852.nc`; only `burndate`, `lfc` are present
  for fire). So although the processes run in Fortran, their state/fluxes are
  **not in the existing dumps** — a new dump boundary is required to validate
  them. The Julia driver also still has `ch4!` as a **placeholder**
  (`clm_driver.jl:1563`) and no fire call.
- Bow surfdata has **PCT_CROP = 0.0** and `PCT_NAT_PFT` summing to 100% (15
  natural PFTs, `PCT_CFT (2,1,1)` but zero crop area). Crop validation
  therefore needs a **different surface dataset with crop CFTs**, not just a
  flag flip.
- `use_c13/use_c14` are OFF → isotopes need a flag flip + rerun.

## Free-win results delivered with EXISTING data (Part 1, done)

1. **CN-driver internals** (`scripts/fortran_parity_cn_midstep.jl`): the unused
   `after_ecosysdyn_predrain` boundary is **byte-identical** to the already-read
   `after_hydrologydrainage` for all CN fields (`after_competition` differs only
   by the `annsum_counter` +3600 s tick — no CN pool). The new probe validates
   the mid-driver allocation storage/xfer pools **and the FUN N-acquisition flux
   split** (`Nuptake/Nactive/Nfix/Necm/Nnonmyc/Npassive`, `availc`, `gpp_pepv`,
   `storage_c/ndemand`) — none of which are injected, so the result is a genuine
   recompute. **Result: 0.000e+00 on every veg group; mineral N (sminn_vr)
   1.9e-4** (small post-drainage N-leaching delta, expected). This is real new
   coverage of the CN allocation + FUN internals beyond the end-of-step pools.

2. **SurfaceAlbedo bands without injection** (`scripts/fortran_parity_albedo.jl`):
   lets Julia COMPUTE albd/albi/fabd/fabi/fsun and the soil/ground bands, then
   diffs vs the dump. **Result: fabd/fabi/fsun match exactly (0.0), but albd/albi
   AND albgrd/albgri/albsod/albsoi are off by 6.5% (max|abs| 0.076).** Root cause
   localized: the soil-albedo moisture blend `inc = max(0.11 − 0.40·h2osoi_vol[1],
   0)` reads `h2osoi_vol_col`, which `read_fortran_restart!` does **not** inject —
   it sits at the cold-start default `0.75·watsat = 0.410` (wet) while Fortran's
   top soil was dry (~0.085, backed out from `albsod=0.166`). Same stale-
   `h2osoi_vol` mechanism MEMORY records for grass-stress/smp_l, now shown to also
   bias soil albedo — **the SABV/SABG check masked it because the harness injects
   the bands.** Fix: inject or recompute `h2osoi_vol` from the injected liq/ice
   before `surface_albedo!`. **Actionable, src-only change.**

---

## Module-by-module roadmap to Fortran parity

Effort: **S** ≤1 day, **M** a few days, **L** ≥1 week / needs new run or data.

### Tier 0 — FREE with existing dumps (no Fortran rerun)

| Module group | Flag (already ON) | What to dump | Julia entry | Effort | Status / blocker |
|---|---|---|---|---|---|
| CN-driver internals (alloc storage/xfer, FUN N-split) | `use_cn/use_fun` | already in `after_ecosysdyn_predrain` | `cn_vegetation_ecosystem_pre_drainage!` outputs | **S** | **DONE** (`fortran_parity_cn_midstep.jl`, 0.0) |
| SurfaceAlbedo bands | always-on | already in every dump (albd/albi/fabd/fabi/fsun/albgr*/albso*) | `surface_albedo!` (recompute) | **S** | **DONE + FIXED 2026-06-17** — 6.5% soil-albedo bug closed to 0.0; `read_fortran_restart!` now recomputes `h2osoi_vol_col` from injected liq/ice (src/infrastructure/fortran_restart.jl). SP parity 30/30, CN parity 11/11 unchanged. |
| `aerosol.jl` deposition | always-on | `mss_bcphi/bcpho/dst1-4/ocphi/ocpho` in dumps | `aerosol_masses!` | **S** | free; diff snow aerosol mass directly |
| `active_layer.jl` | always-on | `altmax/altmax_indx/altmax_lastyear` in dumps | `alt_calc!` | **S** | **DONE + FIXED 2026-06-17** — found `active_layer_init_cold!` never called → altmax stuck at SPVAL; fixed (cold_start.jl) + inject altmax/lastyear (fortran_restart.jl). Now 0.0; asserted in test_fortran_parity_freewins.jl (4/4). |
| `daylength.jl` | always-on | `dayl/prev_dayl` in dumps | `init_daylength!` | **S** | **CHECKED** — `dayl/prev_dayl` NOT in any pdump (only `coszen`); no Fortran ground truth. Self-consistency vs astronomical formula rel 1.2e-9. Not assertable. |
| `surface_resistance.jl` soilresis | always-on | `SOILRESIS` in dumps | `soil_resistance!` | **S** | **DONE** — rel 1.7e-8 (genuine recompute, overwrite-proven); asserted. |
| Decomp cascade rates / N-cycling fluxes | `use_cn` | `f_nit_vr/gross_nmin/potential_immob/decomp_*_rate` in predrain dump | decomp/nitrif outputs | **M** | **PARTIAL** — the `*_VR_P` rate prints are DEAD (identically 0; instrumentation hook never populated) → direct rate parity needs a Fortran dump fix. Integrated effect validated via mineral-N pool deltas (NO3 0.985, NH4/sminn the documented residual). |
| LUNA recompute (don't inject vcmx25_z) | `use_luna` | `vcmx25_z/jmx25_z` after step in dumps | `luna!` | **M** | **PARTIAL + GAP FOUND** — cadence proven (EOD-only, updates exactly at n1757856); but (a) `update_photosynthesis_capacity!` is NEVER called by `clm_drv!` (injection is the only source — wire it in), (b) no-inject recompute blocked: dump omits LUNA accumulators (`t_veg10_day/night`, `fpsn24`, `dayl`) → needs a new dump boundary. Cadence + invariant asserted (9/9). |

### Tier 1 — ONE Fortran rerun with a flag flipped (or a new dump boundary)

| Module group | Flag to flip / dump to add | What to dump | Julia entry | Effort | Blocker |
|---|---|---|---|---|---|
| **Methane (`methane.jl`/ch4)** | `use_lch4` already ON, but ch4 state NOT in restart dump → **add a `restFile_write_dump`-style boundary AFTER the `ch4` call** writing `conc_ch4/conc_o2/ch4_prod/ch4_oxid/ch4_aere/ch4_ebul/ch4_surf_flux/o2stress/finundated` | new boundary `after_ch4` | `ch4!` — but **Julia driver has ch4! as a PLACEHOLDER (`clm_driver.jl:1563`)** | **L** | needs (a) Julia ch4! wired into driver, (b) new Fortran dump boundary + rerun; module already CPU↔Metal validated |
| **Fire (`fire_li2014.jl`)** | `fire_method='li2016crufrc'` ON; only `burndate/lfc` dumped → add a boundary writing `farea_burned/baf_crop/nfire/fire C/N fluxes/cwdc_to_fire` | new boundary `after_fire` (inside EcosysDynPostDrainage) | `cnfire_area_li2014!`, `cnfire_fluxes_li2014!` — **not called in Julia driver** | **L** | needs Julia fire wired into CN driver + new dump boundary + rerun |
| **C13/C14 isotopes** | flip `use_c13=.true.`/`use_c14=.true.`, rerun | `c13_*`/`c14_*` veg+soil pools | `carbon_isotopes.jl` (exists, **never driver-called**) | **L** | **RE-TIERED TO TIER 2 (2026-06-17, tried it):** the 2202 restart has ZERO c13/c14 fields → branch read CRASHES (rc=133); the cold-init escape only works under init_interp, gated off. Needs an isotope-ENABLED spinup, not a flag flip. Also wire `carbon_isotopes.jl` into `clm_drv!`. |
| **MIMICS decomposition** | flip `soil_decomp_method='MIMICSWieder2015'`, rerun | `MIM_*` pools / decomp rates | `decomp_mimics.jl` | **M** | one rerun; **needs a separate MIMICS spinup** (different equilibrium than CENTURY) — so effectively L |
| **VOC / MEGAN** | set `megan_specifier` (non-empty compound list), rerun | VOC emission fluxes (history, not restart → add a dump) | `voc_emission!` (Julia driver L998 placeholder, empty compound list) | **M** | needs MEGAN namelist + a flux dump boundary + Julia call wired |
| **Dust / dry-dep** | `n_drydep>0` + dust always-on | dust emission flux, `depvel` (history fields → add dump) | `dust_dry_dep!`, `depvel_compute!` | **M** | fluxes are history not restart → add a flux dump boundary |
| **Irrigation** | flip `irrigate=.true.`, rerun | `irrig_rate/irrig_rate_demand/n_irrig_steps_left` | `irrigation.jl` (ports exist, **driver calls are empty stubs** ~clm_driver.jl:752/982) | **L** | **RAN OK but NO-OP at Bow (2026-06-17):** `IrrigationMod` gates on `pftcon%irrigated(m)==1` (crop only); Bow has no crop CFT → `n_irrig_steps_left==0`, `irrig_rate` NaN. Needs a crop-bearing site to be a real test. New irrigate=.true. dumps in `tier1_iso/`. |

### Tier 2 — needs new spinup and/or surface data (real NCAR-scale work)

| Module group | What's needed | Julia entry | Effort | Blocker |
|---|---|---|---|---|
| **Crops (`crop.jl` + crop phen/alloc)** | new **surfdata with crop CFTs** (Bow PCT_CROP=0) + `use_crop=.true.` + crop spinup + rerun | `crop.jl`, crop phenology/allocation paths | **L** | needs a crop-bearing gridcell; Bow column can't activate crop code at all |
| **CNDV (dynamic veg)** | `use_cndv=.true.` + a long enough run for biogeography to evolve + rerun | `cndv.jl` | **L** | needs multi-year CNDV spinup; single-step parity won't exercise the establishment/mortality biogeography |
| **Lake landunit** | a **lake-containing gridcell** (Bow has none) + dumps at lake filters | `lake_temperature/fluxes/hydrology/con.jl` | **L** | different domain; lake filter empty in Bow |
| **Urban landunit** | an **urban gridcell** + urban dumps | `urban_albedo/fluxes/radiation.jl` | **L** | different domain; urban filters empty in Bow |
| **Glacier (glcmec)** | a **glaciated gridcell** + glc2lnd coupling | glacier SMB path | **L** | different domain |
| **Hillslope hydrology** | a **hillslope-configured surfdata** + routing wired in driver | `hillslope_hydrology.jl` | **L** | not invoked in Bow driver; needs hillslope columns |
| **FATES** | `use_fates=.true.` (mutually exclusive with the validated CN config) + FATES boundary dumps | `fates_interface.jl` (stub) | **L** | a whole separate validation track; FATES has its own restart format |

---

## Prioritized sprint plan

**Sprint 1 — harvest the free wins (1–2 days, no Fortran rerun).**
Land the two delivered probes as guarded regression tests, then add the
cheap diff-only probes: aerosol snow-mass, active-layer, daylength,
soilresis, and the decomp-rate / N-flux source terms from the predrain dump.
Fix the `h2osoi_vol` injection so the soil-albedo 6.5% closes, then assert
albd/albi parity. Add a no-inject LUNA recompute check (mirror the albedo
pattern) so LUNA's own vcmx25_z update is validated, not bypassed.

**Sprint 2 — one rerun unlocks isotopes + irrigation (M).**
Flip `use_c13/use_c14/irrigate` and rerun the instrumented exe from the
existing 2202 spinup IC; the dump infra already namespaces c13/c14 pools.
Validate isotope pools and irrigation fluxes against the new dumps. No new
spinup needed (isotopes ride existing C; irrigation demand exists on natveg).

**Sprint 3 — wire + dump methane and fire (L).**
These run in the Fortran reference but aren't restart-persisted and aren't in
the Julia driver. Wire `ch4!` and `cnfire_*!` into `clm_drv!` (both modules
are already CPU↔Metal validated), add `after_ch4`/`after_fire`
`restFile_write_dump` boundaries to the SourceMods driver, rerun, and diff.

**Sprint 4 — MIMICS + VOC/dust flux dumps (L).** Separate MIMICS spinup +
rerun for the alternate decomposition; add history-flux dump boundaries for
VOC/dust/dry-dep and wire the Julia calls.

**Sprint 5+ — multi-landunit / multi-PFT tracks (L each).** Crop, CNDV, lake,
urban, glacier, hillslope, FATES each need a different domain or surfdata and
a dedicated spinup; these are the genuine "new ground-truth generation"
sprints to scope with NCAR, not single-column flag flips.

### One-line honesty summary
Tier 0 (free) closes the CN-internals and albedo-band gaps today and already
**found a real 6.5% soil-albedo bug** (stale `h2osoi_vol`); Tier 1 is a handful
of single reruns/dump-boundary edits; Tier 2 is where "full platform parity"
actually requires generating new Fortran ground truth on new domains.
