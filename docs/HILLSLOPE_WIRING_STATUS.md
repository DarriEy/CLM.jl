# Hillslope Hydrology (`use_hillslope`) — Wiring Status

Status as of this branch. Scope: what it would take to run CTSM's hillslope
hydrology (`use_hillslope=.true.`) end-to-end in CLM.jl, what is already ported,
and the concrete blockers that make a *real* end-to-end run infeasible today.

**Bottom line (updated): the reader + multi-column subgrid builder + init wire
are now DONE and validated on a SYNTHETIC catena; a full CTSM-parity run is
still blocked only by the absence of a real hillslope `fsurdat`.** The physics
core (lateral saturated subsurface flow *and* stream routing math) was already
ported. This branch closes the two *plumbing* blockers that were open:
1. **Surfdata reader** — `surfrd_hillslope!` (`src/infrastructure/surfdata.jl`)
   now reads the CTSM `InitHillslope` variables (`nhillcolumns`, `pct_hillslope`,
   `hillslope_index`/`column_index`/`downhill_column_index`,
   `hillslope_slope/aspect/area/distance/width/elevation`, `hillslope_pftndx`,
   and routing `hillslope_stream_*`) into new `SurfaceInputData` fields — GATED
   on `varctl.use_hillslope` and the file actually carrying the variables.
2. **Multi-column subgrid builder** — `count_subgrid_elements` and
   `set_landunit_veg_compete!` build `max(nhillcolumns,1)` columns per ISTSOIL
   landunit (CTSM `subgrid_get_info_natveg`), and `init_hillslope_columns!`
   (`init_gridcells.jl`) wires the `cold`/`colu` catena + per-column geometry
   via the already-ported `init_hillslope!`. Wired into `clm_initialize.jl`
   Step 10b.

All three changes are **strictly gated** on `varctl.use_hillslope` (default
`false`): when off, the natveg build is a single column of weight 1.0 on the
identical code path, so the default run is **byte-identical** (proven by the
full suite + `test/test_hillslope_e2e.jl`'s single-column assertion). The
init-consistency **guard** (`check_aquifer_layer!`) remains wired.

**Still open (honest):** no CTSM-produced hillslope `fsurdat` exists locally, so
there is no Fortran hillslope run to compare against — validation is
port-self-consistency (conservation) + an independent scalar oracle on a
clearly-labeled SYNTHETIC catena, NOT CTSM parity. The full *driver* run on a
multi-column catena (vertical/soil/canopy init per column through `clm_run!`,
plus threading `use_hillslope_routing` + a real `qflx_streamflow_grc` through
`CLMDriverConfig`/`hydrology_drainage.jl:293-297`/`balance_check!`) is not yet
exercised end-to-end and is deferred.

---

## 1. What is ported vs missing

### Ported — physics core (usable once fed real geometry)

- **Lateral saturated subsurface flow between catena columns** — the heart of
  hillslope hydrology — IS ported, including the upland→lowland downhill/uphill
  scatter:
  - `subsurface_lateral_flow!` — `src/biogeophys/soil_hydrology.jl:2705`. Full
    hillslope branch at `soil_hydrology.jl:2893` onward: kinematic/darcy head
    gradient (`HEAD_GRADIENT_KINEMATIC`/`_DARCY`), LayerSum/Uniform
    transmissivity, source/destination column selection via `col.cold`, inflow
    scatter `qflx_latflow_in[col_cold[c]] += …` (`soil_hydrology.jl:3007`), and
    volumetric discharge from the lowest column. (CTSM `SoilHydrologyMod.
    LateralFlowHillslope` / the hillslope additions in `SoilHydrologyMod.F90`.)
  - `perched_lateral_flow!` — `soil_hydrology.jl:2174`. Hillslope branch at
    `soil_hydrology.jl:2300` with the downhill perched-drainage scatter
    (`qflx_drain_perched[cold[c]] -= …`). (CTSM `SoilHydrologyMod.PerchedLateralFlow`.)
  - Both are **already wired into the runtime** — `hydrology_drainage!` calls
    them in the `use_aquifer_layer=false` branch and threads
    `use_hillslope_routing` (`src/biogeophys/hydrology_drainage.jl:275-291`).
  - Device (Metal/CUDA) kernels exist but deliberately restrict to the
    non-hillslope, per-column-independent case (the catena scatter is not a
    race-free per-column kernel); the CPU scalar body handles both paths.
    See `soil_hydrology.jl:2201-2203`, `2736`.

- **Stream channel routing math** — ported and unit-tested, in
  `src/biogeophys/hillslope_hydrology.jl`:
  - `hillslope_stream_outflow!` (`:588`) — Manning discharge + overbank methods.
  - `hillslope_update_stream_water!` (`:683`) — stream volume/depth bookkeeping.
  (CTSM `HillslopeHydrologyMod.HillslopeStreamOutflow` / `…UpdateStreamWater`.)

- **Geomorphology / init routines** — ported (`hillslope_hydrology.jl`):
  `init_hillslope!` (`:146`), `hillslope_properties_init!` (`:63`),
  `set_hillslope_soil_thickness!` (`:287`), `hillslope_soil_thickness_profile!`
  / `…_linear!`, `hillslope_set_lowland_upland_pfts!`, `hillslope_pft_from_file!`,
  and the init guard `check_aquifer_layer!` (`:112`). These accept pre-loaded
  arrays rather than reading NetCDF (the reader is the missing piece — §2).

- **Meteorological downscaling hook** — the topo path threads `use_hillslope`
  through `src/infrastructure/topo.jl` (default `false`); actual hillslope
  meteorology downscaling is explicitly skipped in
  `src/infrastructure/downscale_forcings.jl:20`.

- **Column state** — `ColumnData` carries the full hillslope field set
  (`colu`, `cold`, `hillslope_ndx`, `hill_elev/slope/area/width/distance/aspect`,
  `is_hillslope_column`) — `src/types/column.jl:67-80`. `LandunitData` carries
  `stream_channel_depth/width/length/slope/number` and `WaterStateData` carries
  `stream_water_volume_lun`.

### Missing / no-op

- `hillslope_dominant_lowland_pft!` — **no-op body** (`hillslope_hydrology.jl:526`).
  Not a blocked dependency: `pftcon` is fully ported; only the `find_k_max_indices`
  helper and this routine's own body are unported (PFT distribution methods 2/3).
- The NetCDF **reader** for hillslope surfdata variables — not ported (see §2).
- The **subgrid catena builder** — not ported (see §3).

### `balance_check.jl` expectations

`balance_check.jl` already has `use_hillslope_routing` hooks at 10 sites
(e.g. `:214,:231,:304,:339,:780,:798,:981,:1169-1178`). They only read the
streamflow array *under* `use_hillslope_routing`; when off they pass a dummy, so
the default balance path is unaffected. They are ready to receive a real
`qflx_streamflow_grc` once routing runs.

---

## 2. Surfdata blocker (HARD)

Hillslope needs an `fsurdat` with dims `nhillslope`/`nmaxhillcol` and the
variables CTSM reads in `HillslopeHydrologyMod.InitHillslope`
(`installs/clm/src/biogeophys/HillslopeHydrologyMod.F90:249-399`):
`nhillcolumns`, `pct_hillslope`, `hillslope_index`, `column_index`,
`downhill_column_index`, `hillslope_slope`, `hillslope_aspect`,
`hillslope_area`, `hillslope_distance`, `hillslope_width`,
`hillslope_elevation`, `hillslope_pftndx`, and (routing) `hillslope_stream_*`,
plus `hillslope_bedrock_depth` for the FromFile soil profile.

**None exists locally.** A scan of every `*surfdata*.nc` under
`/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/installs/clm/` for
`nhillslope`/`pct_hillslope` returned zero hits. The port's reader
(`src/infrastructure/surfdata.jl`) has **no** hillslope variable handling.

Fabricating a surfdata is explicitly out of bounds (it would manufacture
geometry rather than test the port), so a real end-to-end run is blocked here
regardless of code. Hillslope surfdata is produced by CTSM's offline
`tools/` hillslope preprocessing from a DEM — not present in this environment.

---

## 3. Subgrid-structure blocker (HARD)

Hillslope represents a **catena**: multiple hydrologically-connected soil
columns per natural-veg landunit, linked by `col.cold`/`col.colu` neighbor
pointers. The port builds the natural-veg (`ISTSOIL`) landunit with a **single**
column:

- `src/infrastructure/init_gridcells.jl:108` —
  `add_column!(col, lun, ci, li[], 1, 1.0)` (one column, weight 1.0).
- Column counts flow from `ncols_per_g` (default `fill(1, numg)`) in
  `src/infrastructure/decomp_init.jl:101`.

There is no code path that creates `nhillcolumns` columns per landunit, assigns
`hillslope_index`/`column_index`/`downhill_column_index`, or wires the
`cold`/`colu` neighbor graph from surfdata. CTSM does this in
`subgridMod`/`initSubgridMod` (hillslope column counting) + `InitHillslope`
(neighbor wiring). Until the builder produces the catena, the lateral solve —
though ported — has no multi-column landunit to run on.

---

## 4. What was wired on this branch (safe, gated, byte-identical)

Only the **init-consistency guard** — faithful to CTSM and provably a no-op for
the default run:

- `src/driver/clm_initialize.jl` — after the aquifer/bedrock flags resolve, call
  `check_aquifer_layer!(varctl.use_hillslope, use_aquifer_layer)`. CTSM calls
  this from `clm_instInit` because hillslope replaces the unconfined-aquifer
  lower boundary with an explicit lateral catena, so the two are mutually
  exclusive. `varctl.use_hillslope` is `false` on every current init path (no
  `clm_initialize!` keyword sets it) → the guard returns `nothing` → **default
  run byte-identical**. It fatals early if a future caller enables both.

- Test: `test/test_driver_defaults_audit.jl` — new testset
  “check_aquifer_layer! guard is WIRED into clm_initialize!” source-scans the
  init module for the call (so the wire can't be silently dropped) and re-checks
  the guard's functional invariant.

**Deliberately NOT wired** (would read absent/NaN geometry or add unreachable
code, violating “don't half-wire something that reads zeros/NaN”):
`init_hillslope!` + `set_hillslope_soil_thickness!` on the init path (no
surfdata, no catena), and the `hillslope_stream_outflow!`/`…update_stream_water!`
runtime calls (the `use_hillslope_routing` stub at
`hydrology_drainage.jl:293-297`) — these can only run atop a catena that does
not exist yet.

---

## 5. Staged plan to finish (in dependency order)

1. **Surfdata**: obtain or generate (via CTSM `tools/` from a DEM) a
   hillslope-enabled `fsurdat` for a single test point (e.g. Bow/Aripuanã) with
   the §2 variables. Prerequisite for everything below; do not fabricate.
2. **Reader**: extend `src/infrastructure/surfdata.jl` to read the §2 hillslope
   variables + dims into host arrays (mirror `InitHillslope`'s reads).
3. **Subgrid catena builder**: extend `decomp_init.jl`/`init_gridcells.jl`/
   `init_subgrid.jl` to (a) count `nhillcolumns` columns per `ISTSOIL` landunit
   and (b) call `init_hillslope!` to set `cold`/`colu`, per-column geometry, and
   recompute column weights. Port CTSM `subgridMod`/`initSubgridMod` hillslope
   column counting.
4. **Init wire**: on the live init path, gated on `use_hillslope`, call
   `hillslope_properties_init!` → subgrid catena build → `init_hillslope!` →
   `set_hillslope_soil_thickness!`. Add a `use_hillslope` keyword to
   `clm_initialize!`/`clm_run!` (currently only a VarCtl field).
5. **Runtime wire**: replace the `hydrology_drainage.jl:293-297` stub with the
   `hillslope_stream_outflow!` + `hillslope_update_stream_water!` calls, gated on
   `use_hillslope_routing`; feed `qflx_streamflow_grc` to the balance-check hooks.
6. **Finish PFT distribution**: implement `find_k_max_indices` +
   `hillslope_dominant_lowland_pft!` body (only needed for PFT methods 2/3).
7. **Validate**: against a Fortran `use_hillslope=.true.` run, or an independent
   oracle of the lateral-flow + stream-outflow math on the catena.

Steps 1–3 are the real work (data + subgrid surgery). Steps 4–6 are then
mechanical wiring of already-ported routines.
