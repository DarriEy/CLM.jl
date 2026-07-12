# CLM.jl Fortran-Parity Harness Coverage Audit

> ## 🕰️ HISTORICAL SNAPSHOT — 2026-06-17. Do not read this as current status.
>
> **Banner added 2026-07-12.** Preserved for its per-module harness inventory.
> Its *wiring* column has since gone stale in at least these rows (verified
> in-tree 2026-07-12): fire IS called from the CN driver (dispatched on
> `cnfire_method`, `src/biogeochem/cn_driver.jl:53`, default `:nofire`);
> methane, VOC, dust, irrigation and the C13/C14 isotope path are all called
> from `clm_drv!`; and FATES is ported and runs live under `use_fates` (its
> Fortran bit-parity is still not established). See the banner in
> [`docs/FULL_PLATFORM_COVERAGE_PLAN.md`](FULL_PLATFORM_COVERAGE_PLAN.md) for the
> line-referenced table, and the README for current parity coverage.

_Read-only audit, 2026-06-17. Scope: which CLM modules the parity harness actually
exercises against Fortran ground truth, vs. which run unvalidated, vs. which are
config-gated off in the validated Bow-at-Banff configuration._

The validated configuration is the **Bow-at-Banff single-column** case:
`use_cn=true`, `use_hydrstress=true` (PHS), `use_luna=true`, `use_aquifer_layer=false`,
soil column (not lake/urban/glacier/crop). Fortran ground truth = per-boundary
restart dumps `pdump_<boundary>_n<nstep>.nc` emitted by an instrumented `cesm.exe`
(SourceMods `restFile_write_dump`) at 9 boundaries per step:
`before_step`, `after_canopyfluxes`, `after_soiltemperature`, `after_soilfluxes`,
`after_hydrologynodrainage`, `after_hydrologydrainage`, `after_ecosysdyn_predrain`,
`after_competition`. Dumps are CLM5 restart-format and very rich (~400 variables),
so the *dump content* is rarely the limiting factor — the limiting factor is which
variables each harness actually reads and diffs.

---

## Honest scope statement

The parity claim is **strong for the biogeophysics + hydrology + snow core, and for
the CN/BGC prognostic pools**, and **weak-to-absent for everything gated off in the
single soil column** (crops, fire, methane, isotopes, VOC, lake, urban, glacier,
CNDV, irrigation, hillslope).

Of the *active* process surface in the Bow config:

- The **biogeophysics/hydrology/snow chain is the most thoroughly validated**: a
  gated CI regression test (`test/test_fortran_parity.jl`) asserts ~15 prognostic
  fields (T_GRND, T_SOISNO, T_VEG, T_STEM, SABV, SABG, EFLX_GNET, ZWT, ZWT_PERCH,
  H2OSOI_LIQ/ICE, WA, H2OSFC, SNOW_DEPTH, frac_sno) within tolerances on a daytime
  and a night step, and a manual `fortran_parity_sweep.jl` diffs ~45 science fields
  at the boundary each is last written. This covers roughly the full always-on
  biogeophysics call tree by *output state*, though it validates **end-of-step
  prognostics**, not every intermediate per-module variable.
- The **CN/BGC pools are validated by single-step and multi-step (drift) scripts**
  (`fortran_parity_cn_summer.jl`, `fortran_parity_drift{,_full}.jl`), reported in
  MEMORY at ~1× agreement for tree/grass availc and the C/N pools, with bounded
  drift (buffer pools cpool/xsmrpool ~3% over a day, mineral N ~1% then plateau).
  **These CN scripts are NOT in the test suite** — they are run-by-hand diagnostics,
  so there is no automated regression guard on CN parity.
- A large class of fields appears in the dumps but is **never diffed by any harness**
  (e.g. surface-albedo bands albd/albi/fabd/fabi, the LUNA vcmax/jmax outputs after
  the step, FUN sub-fluxes, decomposition cascade rates, ozone/dust/aerosol mass) —
  these are validated only transitively (if they were wrong, a downstream diffed
  field would drift), or asserted-by-construction.
- **No process gated off in the soil column is Fortran-validated at the driver
  level at all** (Tables C). Several of these modules *are* internally unit-tested
  and/or GPU-validated (CPU↔Metal), but unit/GPU tests check Julia-vs-Julia, not
  Julia-vs-Fortran.

**Bottom line:** the project's "parity achieved" claims are well-supported for the
*end-of-step prognostic state of the always-on biogeophysics core and the CN/BGC
pools in the Bow column*. They do **not** extend to (a) every intermediate
per-module quantity, (b) any multi-column / non-soil landunit, or (c) the ~15
config-gated process modules, which are validated only against Julia-internal
references, if at all.

> Note on `scripts/gpu_validate_*.jl` (~100 files): these are **CPU-vs-GPU (Metal)**
> equivalence checks, **not** Fortran parity. They prove a kernel runs on the GPU at
> the same answer as its own CPU version; they say nothing about Fortran fidelity.
> Do not count them as parity coverage.

---

## (A) VALIDATED against Fortran ground truth

| Module(s) | Harness | Quantities diffed vs Fortran dump | Agreement (per MEMORY / test tol) |
|---|---|---|---|
| `surface_radiation.jl`, `surface_albedo.jl` (RT outputs) | `test_fortran_parity.jl` (CI, gated); `fortran_parity_sweep.jl` | SABV_P, SABG_P (also fsun, elai, esai in sweep) | tol 5 W/m² in CI; MEMORY: SABV/SABG "exact" after RT injection |
| `canopy_fluxes.jl` / `canopy_fluxes_core!` (energy balance, aero) | `test_fortran_parity.jl`; `fortran_parity_validate.jl`; `sweep` | T_VEG, T_STEM, TAF_P, QAF_P, RAH1_P, RAH2_P, OBU, FV_P, RAM1_P, EFLX_GNET | T_VEG within ~0.02–0.21 K (MEMORY); CI tol 1.2 K |
| `photosynthesis.jl` + `luna.jl` injection + PHS | `fortran_parity_cn_summer.jl`; `probe_phs/availc.jl` | availc, vcmax_z, vegwp, psn (via CN pools downstream) | tree availc ~0.998×, grass ~0.985–0.995× (MEMORY) |
| `soil_temperature.jl` | `test_fortran_parity.jl`; `sweep` | T_SOISNO, T_GRND, T_GRND_R, THK_C | T_GRND ~0.27 K, T_SOISNO tol 0.20 K |
| `soil_fluxes.jl` | `fortran_parity_validate.jl`; `sweep` | EFLX_SHG, EFLX_LH, EFLX_SOIG, EFLX_GNET, EFLX_LWNET, CGRNDS, CGRNDL, T_REF2M | diffed in validate/sweep (no CI tol on most) |
| `bareground_fluxes.jl` | (transitively via T_GRND/fluxes); `sweep` patches | ground SH/LH | indirectly validated |
| `surface_humidity.jl`, `surface_resistance.jl`, `soil_moist_stress.jl` | `fortran_parity_soilparam.jl` (SMP from shared water); BTRAN in `validate.jl` | SMP per-layer, BTRAN | SMP matches to ~4 digits (MEMORY) |
| `snow_hydrology.jl`, `snow_cover_fraction.jl`, `snow_snicar.jl` | `fortran_parity_snow.jl`; `test_fortran_parity.jl`; `sweep` | INT_SNOW, SNOW_DEPTH, frac_sno, frac_sno_eff, snw_rds | INT_SNOW "exact" (MEMORY); CI tol 1e-3 |
| `hydrology_no_drainage.jl`, `canopy_hydrology.jl` | `sweep`; `validate` | LIQCAN, SNOCAN, FWET, H2OSFC | diffed |
| `soil_water_movement.jl`, `soil_hydrology.jl`, `hydrology_drainage.jl` | `test_fortran_parity.jl`; `validate` per-layer ΔH2OSOI; `sweep` | H2OSOI_LIQ/ICE, ZWT, ZWT_PERCH, WA, SMP | CI tol 0.05 (liq/ice); ZWT_PERCH 6.3e-4 (MEMORY) |
| `pre_flux_calcs.jl`, surface water / runoff (`sat_excess_runoff`, `infilt_excess_runoff`, `surface_water`) | (transitive via H2OSFC/H2OSOI in sweep) | column water state | indirectly validated |
| **CN veg C/N pools** (`allocation`, `growth_resp`, `maint_resp`, `phenology`, `cn_driver` veg path) | `fortran_parity_cn_summer.jl` (single step); `fortran_parity_drift{,_full}.jl` (multi-step) | leafc, frootc, livestemc, deadstemc, live/deadcrootc, cpool, xsmrpool, gresp_storage; leafn, frootn, deadstemn, retransn, npool | global CN max\|rel\| ~1×; drift bounded (MEMORY) |
| **Soil BGC decomposition** (`decomp_bgc`, `decomp_cascade`, `litter_vert_transp`) | `cn_summer`; `drift` | litr1/2/3c_vr, soil1/2/3c_vr, cwdc_vr, litr1/soil1/cwdn_vr | leafc/decomp ≤1e-4 drift (MEMORY) |
| **N cycling** (`nitrif_denitrif`, `n_leaching`, `n_dynamics`, `nutrient_competition`) | `cn_summer`; `drift`; `probe_sminbalance/nh4no3/sminabs.jl` | sminn_vr, smin_no3_vr, smin_nh4_vr | mineral N ~1% drift then plateau (MEMORY) |
| **FUN** (`fun.jl`) N uptake | `probe_fun/funuptake/ndemand.jl`; `cn_summer` (via Nuptake/Nretrans) | Nuptake, Nretrans, avail_retransn, plant_ndemand | tree Nuptake →1.00× (MEMORY; commit bd52c1a) |
| **Restart injection** (`fortran_restart.jl`) | `fortran_parity_smoke.jl`; `fortran_parity_cn_read.jl` | round-trip of injected fields ≈ 0 | reader plumbing validated |
| **Soil pedotransfer params** (`cold_start.jl` watsat/bsw/sucsat) | `fortran_parity_soilparam.jl` | SMP recomputed from shared water vs Fortran SMP | matches (forcing-free isolation) |

CI/automation status: **only `test/test_fortran_parity.jl` is a real automated test**
(included in `test/runtests.jl`, but *skips* when the machine-local dumps/inputs are
absent — i.e. green-by-skip on CI). It validates the **SP biogeophysics path only**
(`run_one_parity_step!` with `use_cn=false`). All CN/BGC, PHS, LUNA, FUN parity lives
in **hand-run scripts** with no assertion/regression guard.

---

## (B) RUNS in the Bow config but NOT validated against Fortran

These modules execute every step (or under the Bow flags) and **affect** the
validated outputs, but **no harness diffs their own outputs** against a dump — they
are covered only transitively, or asserted-by-construction.

| Module | Runs because | Why unchecked |
|---|---|---|
| `surface_albedo.jl` (albd/albi/albgrd/albsod/fabd/fabi bands) | always-on radiation | dumps contain albd/albi/fabd/fabi but **no harness diffs them**; only the derived SABV/SABG are checked. RT injection in the harness can mask an albedo bug. |
| `daylength.jl` | always-on | dayl seeded from declination in the harness, not diffed vs a dump field |
| `active_layer.jl` (`alt_calc!`) | always-on | altmax/altmax_indx in dumps, never diffed |
| `aerosol.jl` (`aerosol_masses!`, deposition into snow) | always-on snow path | mss_bcphi/bcpho/dst1-4/ocphi/ocpho in dumps, never diffed; only snw_rds (downstream) is in sweep |
| `dry_dep_velocity.jl` (`depvel_compute!`) | `n_drydep>0` if set; monthly veg read | no dump field diffed |
| `dust_emission.jl` (`dust_dry_dep!`) | always-on | no dust flux diffed |
| `ozone.jl` (`calc_ozone_stress!`) | always-on | ozone stress factor not diffed (would show only via BTRAN/photosynthesis) |
| `surface_resistance.jl` soilresis | always-on | SOILRESIS in dumps, not explicitly diffed |
| `total_water_heat.jl`, `balance_check.jl` | always-on | internal closure check, not a Fortran diff (DYNBAL_* in dumps unused) |
| `lnd2atm_mod.jl` | always-on | gridcell→atm fluxes not diffed |
| `cn_precision_control.jl`, `cn_balance_check.jl` | use_cn | internal C/N conservation, not a Fortran diff |
| `cn_annual_update.jl`, `veg_struct_update.jl`, `veg_compute_seed.jl` | use_cn | annual/structural updates; only the resulting pools (Table A) are diffed, not the per-routine intermediates |
| `gap_mortality.jl` | use_cn | mortality fluxes not separately diffed (folded into pools) |
| `c_state_update{1,2,3}.jl`, `n_state_update{1,2,3}.jl` | use_cn | state-update cascade; validated only via final pools |
| `decomp_potential.jl`, `decomp_competition.jl`, `decomp_precision_control.jl`, `decomp_vertical_profile.jl` | use_cn decomp | rates/cascade not diffed; only resulting decomp pools (Table A). `after_competition`/`after_ecosysdyn_predrain` dumps **exist** but no harness reads them. |
| LUNA outputs (`luna.jl` vcmx25_z/jmx25_z **after** the step) | use_luna | the harness **injects** vcmx25_z/jmx25_z from the dump as input; it does not diff the Julia-recomputed LUNA output (vcmx25_z is input-injected, so LUNA's own update is effectively bypassed/unvalidated) |
| PHS solver internals (`soil_water_plant_sink.jl`, `root_biophys.jl`) | use_hydrstress | vegwp diffed at injectable steps (MEMORY <0.5%), but k_soil_root and the per-segment Newton path are not field-diffed |

The unused boundary dumps `after_competition` and `after_ecosysdyn_predrain` are the
clearest gap: Fortran ground truth for the mid-CN-driver state exists on disk but **no
harness compares against it**, so CN errors are only caught at end-of-step.

---

## (C) NOT exercised at all in the validated config (gated OFF)

These modules have a Julia port and (mostly) internal unit tests + GPU validation,
but the **Bow soil-column driver never enters their code path**, so there is **zero
Fortran-driver validation** of them. The gating condition is the disabling flag/type.

| Module | Gating flag / condition (OFF in Bow) | Other coverage that exists |
|---|---|---|
| `methane.jl` (ch4) | `use_lch4=false` | `test_methane.jl` (unit); `gpu_validate_ch4*` (CPU↔Metal) |
| `carbon_isotopes.jl`, `c_iso_flux.jl` | `use_c13=false`, `use_c14=false` | `test_c_iso_flux.jl`; `gpu_validate_c_iso_flux` |
| `voc_emission.jl` | not called (empty MEGAN compound list, driver L997) | `test_voc_emission.jl`; `gpu_validate_voc_emission` |
| `fire_base.jl`, `fire_li2014.jl` | no fire call in `clm_drv!` (would need use_cn + a fire gate; not wired) | `test_fire_*.jl` (incl. FD check); `gpu_validate_fire_*` |
| `crop.jl` + crop phenology/alloc paths | `use_crop=false`; no crop-itype patches | `test_crop.jl`; alloc pipeline GPU validate |
| `irrigation.jl` | `irrigate=false` | `test_irrigation.jl`; `gpu_validate_irrigation` |
| `cndv.jl` (dynamic veg) | `use_cndv` not active in driver | `gpu_validate_cndv_light` (CPU↔Metal) |
| `lake_*.jl` (con/fluxes/hydrology/temperature) | landunit ≠ lake (filter masks empty) | `test_lake_*.jl`; `gpu_validate_lake*` |
| `urban_*.jl` (albedo/fluxes/radiation) | landunit ≠ urban (urban filters empty); urban_radiation is *called* but no-ops | `test_urban_*.jl`; `gpu_validate_urban*` |
| glacier (glcmec path) | landunit ≠ glcmec | none Fortran; covered by surfdata/init |
| `hillslope_hydrology.jl` | no hillslope routing invoked in driver | `test_hillslope_hydrology.jl` (unit) |
| `fates_interface.jl` | `use_fates=false` (CN active instead) | interface stub |
| `cn_products_mod.jl` (wood/crop product pools) | partition only fires for crop/harvest; ~0 in Bow | `gpu_validate_cn_products` |
| `satellite_phenology.jl` (prescribed LAI) | suppressed when `use_cn=true` | `test_satellite_phenology.jl`; GPU validate |
| `decomp_mimics.jl` | `decomp_method` ≠ MIMICS (Bow uses century/BGC) | `test_decomp_mimics.jl`; `gpu_validate_decompmimics` |

For every Table-C module, "tested" means **Julia-internal** (unit invariants and/or
CPU-vs-GPU equivalence). None has been run through `clm_drv!` against a Fortran dump,
because the Bow soil column never activates them.

---

## Files referenced (all verified to exist)

Harnesses: `scripts/fortran_parity_common.jl` (build_bow_inst, run_one_parity_step!,
compare_inst_to_dump, _parity_registry), `fortran_parity_validate.jl`,
`fortran_parity_sweep.jl` (sweep_registry, ~45 fields), `fortran_parity_smoke.jl`,
`fortran_parity_snow.jl`, `fortran_parity_soilparam.jl`, `fortran_parity_cn_read.jl`,
`fortran_parity_cn_summer.jl` (file header reads `cn_step`), `fortran_parity_drift.jl`,
`fortran_parity_drift_full.jl`. Probes: `scripts/probe_{phs,availc,frootmr,fun,
funuptake,ndemand,nh4no3,sminabs,sminbalance,allometry,1757873}.jl`.

Tests: `test/test_fortran_parity.jl` (the **only** Fortran-diff test in
`test/runtests.jl`, gated-skip when dumps absent; SP path only).
`test/verify_vs_fortran.jl` exists but is **not** in `runtests.jl`.

Driver: `src/driver/clm_driver.jl` (`clm_drv_core!`). Reference dumps:
`/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup/bgc_ref_summer/`
and `.../clm_parity_run/` (9 per-boundary dump types, ~400 vars each).
