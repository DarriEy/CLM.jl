# GPU port — the remaining gap (Phase 1 census result)

Authoritative burn-down list for getting **all of CLM outside BGC** (and the last BGC
diagnostics) onto the GPU. Produced by `scripts/gpu_host_loop_census.jl` (AST host-loop
detector) + a classification pass over every flagged function.

## Status update (post-census burn-down)

The clusters below have since been kernelized + validated on Metal (per-cluster
`scripts/gpu_validate_*` harnesses):

- **Cluster A (urban, ~34 loops)** — DONE (`urban_albedo!` + canyon-RT + `building_temperature!` + `wasteheat!`).
- **Cluster B (BGC summaries, ~87 loops)** — DONE (`soil_bgc_{carbon,nitrogen}_state_summary!`, `cnveg_*_summary!`).
- **Cluster C (~15 misc)** — DONE (`dust_emission_zender2003!`, `cn_veg_struct_update!`, `nitrif_denitrif`, glacier SMB, dry-dep, `_fun_p2c!`).
- **Feature-gated** — DONE (irrigation, LUNA `update_photosynthesis_capacity!`, flexibleCN, `dyn_cnbal_*`).
- **Forcing downscaling (`downscale_forcings.jl`, 7 loops)** — DONE. All per-column
  loops of `downscale_forcings!` / `partition_precip!` / `downscale_longwave!` are now
  KernelAbstractions kernels (`_downscale_baseline/temp/pco2/…_kernel!`). Host KA.CPU
  path proven **byte-identical** to the original loop; Metal Float32 parity ~1e-7 over
  all 12 output fields. Harness: `scripts/gpu_validate_downscale_forcings.jl`.

Remaining flagged loops in the physics dirs are host-reference fallbacks (device variant
already exists), `*_init_cold!` construction, parsing, accessors, reverse-AD, or the
config-gated `hillslope_hydrology.jl` (non-default) — no genuine default-hot un-kernelized
physics remains in `biogeophys`/`biogeochem`. The original per-cluster tables below are
kept for historical reference.

## Headline

The raw census flags **1,511** `for`-loops outside `@kernel`. After classifying each as
gap vs. host-reference / init / I/O / dead, the **real remainder is ~136 default-hot loops
in a handful of clusters**, plus ~67 behind non-default feature flags. Everything else is:
a host-reference path kept beside an existing device kernel (`isa Array ? host : device`),
once-per-run construction (`*_init_cold!`, decomposition/subgrid build), file I/O
(history/restart/surfdata/forcing), reverse-AD, or dead/unwired code.

| bucket | loops | in scope? |
|---|---|---|
| **Default per-step, host-only (real gap)** | **~136** | **yes — Phases 3–4** |
| Feature-gated (crop/irrigation, LUNA, flexibleCN, transient land-use) | ~67 | yes, lower priority |
| Host-reference (device variant already exists) | ~700 | no — already ported |
| Init / setup / topology (once per run) | ~400 | no |
| I/O / restart / history | ~150 | no (stays host) |
| Reverse-AD / dead / unwired | ~60 | no |

## The real gap — default per-step, no device variant

Ordered by priority. Template = the proven kernelization pattern each needs.

### A. Urban biogeophysics — the one substantive physics cluster (~34 loops)
The default CLM5 urban path is entirely un-kernelized; blocks a whole-driver GPU run with urban landunits.

| file:line | function | loops | template |
|---|---|---|---|
| `urban_albedo.jl:658` | `urban_albedo!` | 12 | grouped-struct (per-landunit) |
| `urban_albedo.jl:254` | `net_solar!` | 7 | grouped-struct (30+ matrix args) |
| `urban_albedo.jl:140` | `incident_direct!` | 3 | grouped-struct |
| `urban_albedo.jl:210` | `incident_diffuse!` | 2 | grouped-struct |
| `urb_build_temp_oleson2015.jl:87` | `building_temperature!` (+ `_solve5x5!` inner) | 9 | grouped-struct + device-callable 5×5 GE inner (like the PHS `_*_core!`) — **default** `is_prog_build_temp()` |
| `urban_fluxes.jl:63` | `wasteheat!` | 1 | elementwise (per-landunit) |

### B. BGC state summary reductions (~87 loops, mostly one template)
Per-step diagnostic pool-total aggregations reached each step via `clm_drv!`. Reductions, not new physics.

| file:line | function | loops | template |
|---|---|---|---|
| `soil_bgc_nitrogen_state.jl:381` | `soil_bgc_nitrogen_state_summary!` | 43 | per-column reduction |
| `soil_bgc_carbon_state.jl:355` | `soil_bgc_carbon_state_summary!` | 37 | per-column reduction |
| `cn_veg_carbon_flux.jl:1304` | `cnveg_carbon_flux_summary!` | 3 | per-column reduction |
| `cn_veg_carbon_state.jl:608` | `cnveg_carbon_state_summary!` | 2 | per-column reduction |
| `cn_veg_nitrogen_state.jl:643` | `cnveg_nitrogen_state_summary!` | 2 | per-column reduction |

### C. Misc default-hot (~15 loops)
| file:line | function | loops | template | note |
|---|---|---|---|---|
| `dust_emission.jl:681` | `dust_emission_zender2003!` | 3 | elementwise per-patch | **default** dust scheme (Leung2023 sibling already kernelized) |
| `veg_struct_update.jl:44` | `cn_veg_struct_update!` | 1 | grouped-struct | zero-kernel file; feeds radiation + photosynthesis |
| `nitrif_denitrif.jl:552` | `soilbiogeochem_n_state_update1!` | 1–2 | grouped-struct (~15 flux fields) | `use_nitrif_denitrif` default true |
| `glacier_surface_mass_balance.jl:29/96/172` | ice-melt / SMB / runoff-adjust | 6 | elementwise | glacier landunit |
| `dry_dep_velocity.jl:505` | `depvel_compute!` | 2 | grouped-struct | per-patch dry-dep velocity |
| `fun.jl:86` | `_fun_p2c!` | 2 | patch2col-scatter | minor; `use_fun` |
| `soil_water_movement.jl:1644` | `soil_water!` flexibleCN fixup | 5 | elementwise (sequential vertical) | config-gated |

## Feature-gated (real, but not on the default path — lower priority)

| file / function | loops | gate |
|---|---|---|
| `irrigation.jl` (6 fns: need/apply/withdrawals/p2c/c2g) | 15 | crop + irrigate |
| `luna.jl` (6 fns → ~4 kernel entry points; `update_photosynthesis_capacity!` is the core) | 12 | `use_luna` |
| `nutrient_competition_flexiblecn.jl` | 14 | `use_flexiblecn` (default false; `nutrient_competition.jl` IS kernelized) |
| `dyn_cons_biogeochem.jl` (`dyn_cnbal_patch!/col!`) | 26 | transient land-use, via `dyn_subgrid_driver` |

## Explicitly out of scope
- **FATES** (65 files) — separate ecosystem-demography subsystem; own track.
- **Restart / history / surfdata / forcing I/O** — stays host (GPU doesn't accelerate I/O; only gather-at-write needed).
- **`*_init_cold!` / `*_set_values!` / topology build** — once-per-run construction.
- **Reverse-AD** (`canopy_fluxes_reverse.jl`, `driver_reverse.jl`) — off the forward device path.
- **Host-reference paths** — kept intentionally beside device kernels for byte-identity + AD.

## Suggested sequencing (feeds Phases 3–4 of the plan)
1. **Phase 2 first** — validate `clm_drv!` `use_cn=true` (biogeophys + BGC composite) on Metal; that decides whether module-boundary precision holds before more kernelizing.
2. **Cluster B** (BGC summaries, ~87 loops, one reduction template) — highest loop-count, lowest risk; unblocks a clean BGC-integrated device timestep.
3. **Cluster A** (urban, ~34 loops) — the one substantive physics cluster; needed for urban-landunit runs.
4. **Cluster C** (~15 misc) + glacier/dry-dep for full default-landunit coverage.
5. Feature-gated (irrigation, LUNA, flexibleCN, dyn_cnbal) as those configs are targeted.
