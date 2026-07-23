# CTSM Hillslope `hillslope_file` — Variable Spec + Parity-Test Surfdata

Source of truth (Fortran CTSM at `installs/clm`):
- `src/biogeophys/HillslopeHydrologyMod.F90::InitHillslope` (reads the geomorphology)
- `src/main/surfrdMod.F90::surfrd_hillslope` (reads dims + `nhillcolumns`, sets column count)
- `src/main/surfrdMod.F90::check_domain_attributes` (lon/lat must match the land domain)
- `src/main/controlMod.F90`, `clm_varctl.F90`, `clm_initializeMod.F90` (namelist wiring)

## 1. How CTSM enables hillslopes

Hillslope geomorphology is **NOT** read from `fsurdat`. It is a **separate netCDF file**
given by the `hillslope_file` namelist entry (`&clm_inparm`, `char*256`, `input_pathname="abs"`).
`surfrd_get_data(begg,endg,ldomain,fsurdat,hillslope_file,...)` opens both; the hillslope
file is validated against the land domain (`check_domain_attributes`, needs `LONGXY/LATIXY`
or `xc/yc`), then `surfrd_hillslope` reads its dims and `InitHillslope` reads the geometry.

Namelist flags (all `&clm_inparm`, default `.false.` unless noted):
| flag | group | effect |
|---|---|---|
| `use_hillslope` | clm_inparm | master toggle; requires `hillslope_file` set |
| `use_hillslope_routing` | clm_inparm | enables stream routing → reads the 3 `hillslope_stream_*` vars |
| `downscale_hillslope_meteorology` | clm_inparm | met downscaling across columns |
| `hillslope_fsat_equals_zero` | clm_inparm | fsat=0 on hillslope columns |
| `hillslope_head_gradient_method` | hillslope_hydrology_inparm | `Kinematic` \| `Darcy` |
| `hillslope_transmissivity_method` | hillslope_hydrology_inparm | `LayerSum` \| `Uniform` |
| `hillslope_pft_distribution_method` | hillslope_properties_inparm | `Standard`\|`FromFile`\|`DominantPftUniform`\|`DominantPftLowland`\|`PftLowlandUpland` |
| `hillslope_soil_profile_method` | hillslope_properties_inparm | `Uniform`\|`FromFile`\|`SetLowlandUpland`\|`Linear` |

Hard constraint: `use_hillslope` and `use_aquifer_layer()` **cannot both be true**
(`check_aquifer_layer`). `use_aquifer_layer()` is true iff the soil lower boundary condition
is `aquifer` or `watertable`; so hillslope requires a flux/zero-flux lower BC.
`n_dom_pfts` is **not** required. Hillslope is pure hydrology → runs under SP (no BGC/FATES needed)
and is compiled unconditionally (no CPP guard; `HillslopeHydrologyMod.o` present in the SP build).

## 2. Exact variable spec

Dimensions on `hillslope_file`:
- grid dimension = **`gridcell`** (unstructured/SE grids) **or** **`lsmlat`,`lsmlon`** (structured; e.g. single point 1×1). CTSM reads with `dim1name=grlnd`.
- `nhillslope` — number of representative hillslopes per landunit (sets `clm_varctl::nhillslope`).
- `nmaxhillcol` — max hillslope columns per landunit (sets `clm_varctl::max_columns_hillslope`).

CDL dimension order = **extra dim first, then grid dims** (same convention as
`PCT_NAT_PFT(natpft, gridcell)` / `(natpft, lsmlat, lsmlon)`).

| variable | dims (CDL, gridcell grid) | type | units | notes |
|---|---|---|---|---|
| `LONGXY` | `(gridcell)` / `(lsmlat,lsmlon)` | double | degrees east | domain check |
| `LATIXY` | `(gridcell)` / `(lsmlat,lsmlon)` | double | degrees north | domain check |
| `nhillcolumns` | `(gridcell)` | int | unitless | # columns in gridcell; 0 ⇒ no hillslope (uses standard patches) |
| `pct_hillslope` | `(nhillslope, gridcell)` | double | percent | % of landunit occupied by each hillslope; per-gridcell sums to 100 over active hillslopes |
| `hillslope_index` | `(nmaxhillcol, gridcell)` | int | unitless | which hillslope each column belongs to (1..nhillslope) |
| `column_index` | `(nmaxhillcol, gridcell)` | int | unitless | external column id, spans 1..nhillcolumns |
| `downhill_column_index` | `(nmaxhillcol, gridcell)` | int | unitless | downstream neighbor's `column_index`; **≤ -999 ⇒ lowland, discharges to stream** (→ `col%cold=ispval`) |
| `hillslope_slope` | `(nmaxhillcol, gridcell)` | double | m/m | along-hill slope |
| `hillslope_aspect` | `(nmaxhillcol, gridcell)` | double | radians | azimuth (0..2π) |
| `hillslope_area` | `(nmaxhillcol, gridcell)` | double | m² | column surface area; column `wtlunit` recomputed = (area/Σarea)·(pct·0.01) |
| `hillslope_distance` | `(nmaxhillcol, gridcell)` | double | m | distance of column lower edge from hill bottom |
| `hillslope_width` | `(nmaxhillcol, gridcell)` | double | m | width of column lower edge |
| `hillslope_elevation` | `(nmaxhillcol, gridcell)` | double | m | mean column elevation rel. to gridcell mean |
| `hillslope_pftndx` | `(nmaxhillcol, gridcell)` | int | unitless | **optional**; used only by `FromFile`/`PftLowlandUpland` pft methods |
| `hillslope_stream_depth` | `(gridcell)` | double | m | **routing only** → `lun%stream_channel_depth` |
| `hillslope_stream_width` | `(gridcell)` | double | m | **routing only** → `lun%stream_channel_width` |
| `hillslope_stream_slope` | `(gridcell)` | double | m/m | **routing only** → `lun%stream_channel_slope` |

`InitHillslope` maps external→internal columns by relative offset within the landunit's
column range: `col%cold(c) = c + (downhill_column_index - column_index)`; uphill neighbor
`col%colu` is found by inverse scan. It also `endrun`s if `nhillcolumns>0` but Σ`hillslope_area==0`.
(`hill_bedrock` is declared but only read by `SetHillslopeSoilThickness` for the
`FromFile` soil profile method — not part of the core geometry read.)

## 3. The parity-test surfdata (real base, synthetic geometry)

Generator: `scripts/make_hillslope_surfdata.jl <base_surfdata.nc> <out.nc>`
(auto-detects `gridcell` vs `lsmlat/lsmlon`). Two files were produced:

- `hillslope_ne3np4_synthetic_c.nc` — base = present ne3np4 surfdata (`gridcell=488`).
- `hillslope_aripuana_synthetic_c.nc` — base = SYMFLUENCE single-point Aripuanã surfdata (`lsmlat=lsmlon=1`).

**Real base / synthetic geometry**: `LONGXY`/`LATIXY` are copied verbatim from the real base
surfdata (so `check_domain_attributes` passes against that grid's mesh). The catena itself is an
**idealized, mass-consistent synthetic geometry** — it is *not* from a DEM. One gridcell (the most
vegetated, `PCT_NATVEG=100`) hosts the catena; all others get `nhillcolumns=0`.

Synthetic catena: `nhillslope=1`, `nmaxhillcol=4`, columns ordered `ci=1` (lowland, at stream) →
`ci=4` (upland):
- `column_index=[1,2,3,4]`, `hillslope_index=[1,1,1,1]`
- `downhill_column_index=[-999,1,2,3]` (ci1→stream, chain upward)
- `hillslope_distance=[0,25,50,75] m`, `hillslope_width=[100,90,80,70] m`
- `hillslope_area=[2500,2250,2000,1750] m²` (= width·25 m, Σ=8500)
- `hillslope_elevation=[1,6,11,16] m` (monotonic up), `hillslope_slope=[0.05,0.08,0.10,0.12] m/m`
- `hillslope_aspect=π` (south), `hillslope_pftndx=13` (c3 grass)
- `pct_hillslope=100`; stream: depth 2 m, width 5 m, slope 1e-3
- recomputed `col%wtlunit` = [0.2941,0.2647,0.2353,0.2059], **Σ=1.0** ✔

Suggested namelist for a Fortran parity run (SP, no aquifer):
```
use_hillslope = .true.
use_hillslope_routing = .true.
hillslope_file = '<abs path to hillslope_*_synthetic_c.nc>'
hillslope_pft_distribution_method = 'Standard'
hillslope_soil_profile_method     = 'Uniform'
```

## 4. Port reader status

There is **no** `surfrd_hillslope!` netCDF reader on `origin/main` (the task's premise was
inaccurate). The port's hillslope entry point is `init_hillslope!(...)` in
`src/biogeophys/hillslope_hydrology.jl`, which takes the arrays **already read**. Verified by
`scripts/verify_hillslope_read.jl`: it reads the file with NCDatasets (exact CTSM var names/dims)
and feeds `init_hillslope!`, which reproduces the catena (downhill/uphill chains, areas,
`wtlunit` Σ=1, stream props). A production port run still needs a small netCDF→arrays reader
(surfrd-equivalent) wired into the driver.
