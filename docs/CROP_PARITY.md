# Crop / irrigation / fire-crop Fortran parity — status (2026-07-18)

> **SUPERSEDED IN PART (2026-07-19).** The crop lifecycle is now wired and
> validated against this document's reference — see **`docs/CROP_LIFECYCLE_MAP.md`**,
> which is authoritative for the lifecycle. In particular:
> - The § "Wiring `n_fert!` — it would be a NO-OP today" analysis was correct and
>   its precondition is now MET: `crop_phenology!` computes `fert_patch`, whose
>   value matches the May reference to relative error 0.0, so `n_fert!` is wired.
>   **`n_soyfix!` is still NOT wired** — `SOYFIXN` ≡ 0 in both windows. The cause
>   is now fully diagnosed and is **not** a window/phase issue: these runs have
>   `use_fun = .true.`, and CTSM calls `CNSoyfix` only under `.not. use_fun`, so
>   it never executed. Non-FUN runs were generated and hit a second wall,
>   `fpg ≡ 1` (the site is N-replete). See `docs/CROP_LIFECYCLE_MAP.md` §3c-§3g.
> - The § "Verified Julia wiring state" rows for `crop_phenology!` /
>   `plant_crop!` / `vernalization!` are out of date (all three now have live call
>   sites).
> - **Path correction:** the reference lives at
>   `SYMFLUENCE_data/clm_bgc_spinup/crop_ref_usplains/` — directly under
>   `SYMFLUENCE_data`, *not* under `installs/` as written below.
> - The **restart file** (`Crop_USplains.clm2.r.2020-05-30-14400.nc`) is a richer
>   reference than the bgcdumps and was unused by this document: it carries
>   `idop`, `gddmaturity`, `huileaf`, `huigrain`, `fert_counter`, `onset_counter`,
>   letting each lifecycle step be validated in isolation.

Backlog items **A3 (crops + crop N)**, **A4 (irrigation)** and the **A2 fire crop
branch** share a single Fortran ground-truth run: a crop-bearing spin-up with
`use_crop = .true.`, `use_cn = .true.`, `irrigate = .true.`. This documents an
attempt to generate that reference at the candidate site
`domain_Cropland_Mead_USA`, the **precise blocker** that stopped it, and the
verified state of the Julia crop/irrigation/fire-crop wiring (which side of each
gap is code vs. reference).

## TL;DR

- **UPDATE 2026-07-18 — a crop-CFT surfdata now EXISTS (the "no surfdata"
  blocker is retired).** The earlier "no `inputdata` tree / no global crop
  surfdata" conclusion was drawn from the Drive-only Mead migration copy. The
  **local box does have an inputdata tree** —
  `installs/cesm-inputdata/` — and it ships a **global crop-resolved surfdata**:
  `lnd/clm2/surfdata_esmf/ctsm5.4.0/surfdata_ne3np4_hist_2000_78pfts_c251022.nc`
  (`cft = 64`, `natpft = 15`, `lsmpft = 79`; 488 gridcells on the ne3np4
  cubed-sphere) plus its ESMF mesh + domain. A **single-point crop surfdata** has
  been extracted from it (US Great Plains corn/soybean/wheat cell, `PCT_CROP=45%`)
  and validated CLM.jl-readable. See § "Crop-CFT surfdata OBTAINED" below. What
  remains for the Fortran reference is only multi-year datm forcing at the point +
  a CN-crop spin-up (the serialized Fortran-box step), **not** a missing surfdata.
- **SYMFLUENCE cannot build a crop surfdata itself.** Its
  `CLMSurfaceGenerator.generate_surface_data` is hardwired to the SP layout
  (`cft=2, natpft=15, lsmpft=17`) and `CLMNuopcGenerator` writes `use_crop=.false.`
  into `lnd_in` with no crop config surface. So the crop surfdata had to be
  extracted from the shipped global 78-pft file (option 2), not produced via
  SYMFLUENCE (option 1).
- **Original blocker (now superseded), for the record:** No crop-resolved surfdata
  was thought to exist; the migrated Mead assets are an **SP run**
  (`use_crop/use_cn/irrigate = .false.`, `cft = 2`), whose dumps are
  snow/hydrology-only. The stock `subset_data` tool still cannot slice the ne3np4
  file (it expects a structured `lsmlat/lsmlon` global surfdata, not the
  unstructured `gridcell`-dimensioned one) — hence the direct-reshape extract.
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

### Why the *Mead* assets are not a crop reference (original finding)

- The migrated Mead assets on the Drive are the SP `cft=2` run above.
- The stock `tools/site_and_regional/subset_data` still cannot be pointed at the
  ne3np4 global crop surfdata: that tool's `SinglePointCase` reads a **structured**
  (`lsmlat/lsmlon`) global surfdata (`create_1d_coord(..., "lsmlon","lsmlat")`),
  but the shipped 78-pft file is **unstructured** (`gridcell`-dimensioned). Hence
  the direct-reshape extract below rather than a `subset_data` call.
- `tools/mksurfdata_esmf` would need the raw `mksrf_*` PCT_CFT / land-use rasters,
  which are not on disk — a large, network-heavy download, and unnecessary now
  that a global crop surfdata is available to slice.

## Crop-CFT surfdata OBTAINED (2026-07-18)

The local inputdata tree ships a **global crop-resolved surfdata**:

    installs/cesm-inputdata/lnd/clm2/surfdata_esmf/ctsm5.4.0/surfdata_ne3np4_hist_2000_78pfts_c251022.nc

with `cft = 64`, `natpft = 15`, `lsmpft = 79` over 488 ne3np4 cubed-sphere
gridcells (real crop area — many cells `PCT_CROP > 40%`), plus its ESMF mesh
(`share/meshes/ne3np4_ESMFmesh_c230714_cdf5.nc`) and domain
(`share/domains/domain.lnd.ne3np4_gx3v7.230718.nc`). This is a complete,
known-good global crop configuration.

**Single-point crop surfdata built from it** (the new asset), via
`scripts/build_crop_cft_surfdata.py` → stored on the Drive at
`data/SYMFLUENCE_data/crop_cft_surfdata/`:

| File | What |
|---|---|
| `surfdata_cropCFT_USplains_1pt.nc` | `lsmlat=1,lsmlon=1,natpft=15,cft=64,lsmpft=79`; `PCT_CROP=45.0`, `PCT_CFT` sums 100; lon 263.29 (−96.71) / lat 44.80. Dominant CFTs: temperate soybean (35%), temperate corn (31%), spring wheat (12%). |
| `esmf_mesh_croppt.nc` | matching 1-element ESMF mesh, `centerCoords == LONGXY/LATIXY`. |

Built by extracting **gridcell 352** and applying the same
`gridcell → (lsmlat=1, lsmlon=1)` reshape that
`subset_data.create_surfdata_at_point` performs (done directly on the
unstructured file). Validated **without the Fortran box**: all CLM init-required
vars present (`AREA` excepted — mesh-supplied in NUOPC), and CLM.jl's own reader
(`surfrd_get_num_patches`/`surfrd_get_grid_dims`) reads it as
`numpft=15, numcft=64, grid=(1,1)`.

### Runnable crop config (recipe)

Mirror the working single-point MerBleue CN case
(`domain_Peatland_MerBleue_Canada/settings/CLM/`), swapping in the crop assets:

1. `fsurdat = …/surfdata_cropCFT_USplains_1pt.nc`; land + datm `meshfile =
   …/esmf_mesh_croppt.nc`. (Alternatively, create a SYMFLUENCE domain whose
   centroid is 263.29/44.80 so its generated mesh matches, then overwrite the
   generated `cft=2` surfdata with the crop one.)
2. `lnd_in`: `use_cn=.true.  use_crop=.true.  irrigate=.true.
   use_fertilizer=.true.  use_grainproduct=.true.  create_crop_landunit=.true.
   finidat=''` (cold start). SYMFLUENCE writes `lnd_in` directly (bypassing
   `bldnml`), so these are set by editing the generated `lnd_in`.
3. Multi-year datm forcing at the point — **the one remaining acquisition**
   (GSWP3 or ERA5-derived, as for the other domains) — plus a CN-crop spin-up long
   enough for a crop to plant / accumulate GDD-HUI / grow / harvest.
4. **Open question for the run step:** whether the shared `cesm.exe` engages crops
   from the namelist flags alone. CTSM compiles the crop code unconditionally
   (`use_crop` is a runtime flag, not an ifdef), so it *should*; if the build
   rejects `use_crop=.true.` (e.g. it was configured `-bgc sp`), rebuild the case
   `-bgc bgc -crop`. This needs the Fortran box and is deferred as the serialized
   stretch (a sibling agent holds the box).

**To finish A2/A3/A4, the only remaining prerequisite is the datm forcing + a
CN-crop cold-start spin-up** long enough for the crop to plant, accumulate
GDD/HUI, grow and harvest, plus crop pdump/BGCDUMP instrumentation. The site need
not be Mead — single-step parity validates the crop *code path* against the same
inputs.

## Fortran-box survey (2026-07-18) — the run step is now scripted

The remaining prerequisites above were surveyed on the box and are resolved.
`scripts/setup_crop_ref_case.py` builds the whole case in one command; the
findings behind it:

### OPEN QUESTION ANSWERED — **no rebuild is needed**

The shared `cesm.exe` engages crops **from the namelist flags alone**. Two
independent lines of evidence:

1. **No CPP gate exists.** CTSM has no `#ifdef` on crop anywhere in `src/`;
   `use_crop` is a runtime `clm_varctl` flag branched on in normal code (e.g.
   `clm_varpar.F90:329`). `CLM_BLDNML_OPTS="-bgc sp"` steers `bldnml`'s
   *namelist defaults* — it is not a compile-time configuration.
2. **The crop code is already in the binary.** `nm` on the built
   `cases/symfluence_build/bld/cesm.exe` (the `-bgc sp` one every existing
   reference uses) finds the full crop path compiled in:

   | Module | Symbols | Notable |
   |---|---|---|
   | `CropType` | 15 | `cropupdateaccvars`, `initcold`, `latbaset`, `restart` |
   | `irrigationMod` | 33 | the whole A4 path |
   | `CNPhenologyMod` | — | `cropphenology`, `cropphase`, `cncropharvesttoproductpools` |
   | `CNNDynamicsMod` | — | **`cnnfert`, `cnsoyfix`** — the two unwired routines |
   | `CNFire*` | 68 | the A2 `baf_crop` branch |

   `PlantCrop`/`Vernalization` are `contains`-scoped inside `CropPhenology` and
   so carry no separate external symbols; their enclosing routine is present.

This retires the last deferred unknown from the § "Runnable crop config" recipe.

### Parameter file — covers the crop CFTs unchanged

The Bow CN reference paramfile
(`domain_Bow_at_Banff_lumped/.../clm5_params.nc`) is dimensioned `pft = 79` and
carries every crop parameter (`gddmin`, `hybgdd`, `mxmat`, `baset`, `fertnitro`,
`declfact`, …). No crop-specific paramfile has to be built.

### Forcing at the point — acquisition path + the substitution actually used

There is **no datm forcing at or near the crop point**: the box's `inputdata`
tree ships no `atm/datm7`, and every SYMFLUENCE domain with `clmforc.YYYY.nc`
is elsewhere (MerBleue 45.41N/-75.52 ×2yr; Bow ×11yr; the rest have none).
ERA5 is reachable **without CDS credentials** — SYMFLUENCE pulls it from the
public ARCO-ERA5 zarr on GCS (`gcp-public-data-arco-era5/ar/…`, anonymous), and
`xarray`/`zarr`/`gcsfs` are all installed. `scripts/acquire_crop_forcing.py`
does that pull directly (a single-cell extraction does not need a whole
SYMFLUENCE domain build).

**Cost caveat:** the ARCO store chunks *one global field per (hour, variable)*,
so a 2-year hourly 8-variable pull is ~140k chunk reads and runs for hours — it
is the slowest step by far, and it does not parallelise with the serialized
Fortran box.

**What the generated case therefore uses by default:** the MerBleue 2016–2017
forcing, **relabeled onto the crop point's coordinates** (`--relabel`;
meteorological values copied through unaltered, only coordinate labels change).
This is sound *for a parity test and only for that* — CTSM and CLM.jl read the
identical file, so the diff isolates the crop code path exactly as site-native
met would. It is **not** a climatology of the surfdata's location, so absolute
yields from this run carry no site meaning. The donor (45.41N) and the crop
point (44.80N) are 0.6° apart and share a humid-continental regime that
genuinely grows corn and soybean, so a crop lifecycle remains physically
plausible. Pass `--real-era5` to `setup_crop_ref_case.py` to substitute a
site-native pull once someone can absorb the transfer time.

> **The `--real-era5` path is UNVALIDATED.** `acquire_crop_forcing.py` carries a
> `--validate-merbleue` self-check that re-derives forcing at the MerBleue point
> and correlates it against that domain's independently-produced
> `clmforc.2016.nc`, specifically to prove the ERA5→CLM conversions (the
> accumulated-flux `/3600`, the dewpoint→specific-humidity inversion) before
> they are trusted as a parity input. **It was started and never completed** —
> the ARCO transfer outran the session. So the conversion arithmetic in
> `to_clm()` has *not* been checked against ground truth and must not be
> assumed correct.
>
> This does **not** affect the reference this document reports: the case was
> built with `--relabel`, which copies every meteorological value through
> unaltered and touches only coordinate labels — no conversion is involved.
> Anyone taking the `--real-era5` path must run `--validate-merbleue` to
> completion first.

### The case is a self-contained run dir, not a CIME case

Every working CN reference on this box (`clm_bgc_spinup/bow_ref_*`,
`merbleue_ref_ch4`, `bgc_ref_*`) is run by writing `lnd_in`/`drv_in`/`datm_in`/
`nuopc.*` directly and invoking the shared `cesm.exe` in place — SYMFLUENCE
bypasses `bldnml` and `case.submit` entirely. The crop case mirrors that, so it
stays consistent with the references the parity harnesses already consume.
Instrumentation is already baked into the exe and gated by the
`PDUMP_NSTEP_LO/HI` and `BGCDUMP_NSTEP_LO/HI` environment variables — the
spin-up runs with the window closed, then the diff run re-runs a short bracket
with it open.

**Every flag is set explicitly** in `LND_FLAGS`, including ones that merely
match the donor. Per #251, an unset flag is not a matched flag.

### Three harness bugs caught by verifying the GENERATED namelists

The case generator was not trusted; every field it wrote was read back. That
found three defects, each of which would have produced a "port bug" that was
really a harness bug — the #233/#240/#243/#248/#251 pattern, a sixth and
seventh time:

1. **Mesh path truncated by a space.** The crop assets live under
   `.../My Drive/...`, and `nuopc.runconfig` fields are **unquoted**
   (`mesh_lnd = /path/…`). The space would have truncated the path and dropped
   the run onto a different grid. Fixed by staging the assets to a space-free
   directory on the box; asserted in the script.
2. **`start_ymd` clobbered `restart_ymd`.** An unanchored regex matched the
   substring, overwriting the `-999` sentinel with `20160101`. Fixed by
   anchoring to line start and asserting each key matches exactly once.
3. **The topo stream still pointed at the donor CN case's Bow file — the worst
   of the three.** `topo_forcing.nc` supplies the elevation the atm forcing is
   lapse-rate downscaled *from*: **Bow 2138.8 m vs MerBleue 71.5 m**, a 2067 m
   error, i.e. **~12.4 K** of spurious temperature shift at
   `lapse_rate = 0.006`. That feeds directly into GDD accumulation and the crop
   planting decision — it would have silently prevented or grossly distorted
   the crop lifecycle this whole reference exists to exercise. The topo stream
   must match whichever site the met came from.

Note that none of the three would have raised an error; all three would have
run to completion and produced confident, wrong numbers.

### Reproducing it

```bash
python3 scripts/setup_crop_ref_case.py          # build case (default 4 yr spin-up)
bash /Users/.../clm_bgc_spinup/crop_ref_usplains/run_crop_ref.sh
```

## THE FORTRAN CN-CROP RUN SUCCEEDED (2026-07-18)

`FINAL_OK=1` on the first attempt — 4 model years (2016→2020, nstep 35064),
cold start, **no rebuild**, no macOS-26 xzone retry needed. That confirms the
symbol-table analysis above empirically: `use_crop=.true.` from the namelist
alone engages the crop model in the stock `-bgc sp` `cesm.exe`. CTSM's own
`cropcal_stream DERIVED settings` block appears in `lnd.log`, which only prints
on the crop path.

**A crop plants, grows, and harvests every year.** Annual means from
`Crop_USplains.clm2.h0.2016-01-01-14400.nc`:

| Field | 2016 | 2016-12 | 2017 | 2018 | 2019 |
|---|---|---|---|---|---|
| `CPHASE` | nan | nan | 2.848 | 3.469 | 3.394 |
| `HUI` | 0 | 0 | 129.45 | 107.85 | 124.38 |
| `GDDACCUM` | 0 | 0 | 123.81 | 107.85 | 124.38 |
| `TLAI` | 0 | 0.522 | 1.451 | 1.518 | 1.534 |
| `LEAFC` | 2.56 | 15.36 | 38.92 | 40.27 | 40.39 |
| `TOTVEGC` | 102.1 | 189.2 | 391.7 | 448.6 | 450.6 |
| `TOTVEGN` | 2.99 | 5.00 | 10.02 | 11.86 | 11.49 |
| **`NFERTILIZATION`** | 0 | 0 | 1.360e-7 | 1.360e-7 | 1.360e-7 |
| **`QIRRIG_DEMAND`** | 0 | 0 | 1.750e-8 | 3.466e-8 | 1.504e-8 |
| **`BAF_CROP`** | 0 | 0 | 1.178e-11 | 8.566e-12 | 1.178e-11 |

`HUI_PERHARV` = **1469 / 1481 / 1425** in 2017/2018/2019 (second slot `-1`, the
no-second-harvest sentinel) — i.e. exactly **one completed harvest per year**,
each after ~1400–1500 accumulated heat units. `CPHASE ≈ 3` is grain fill. The
lifecycle is complete: plant → GDD/HUI accumulation → grain fill → harvest.

### Every blocked backlog row now has a live, non-zero reference

| Row | Was | Now |
|---|---|---|
| **A3** crop pools / GDD-HUI / harvest | BLOCKED — no crop-CFT surfdata | `CPHASE`/`HUI`/`GDDACCUM`/`*_PERHARV` all populated, 3 harvests |
| **A3** crop N (`n_fert!`/`n_soyfix!`) | unwired, no reference | `NFERTILIZATION` = 1.360e-7, non-zero and steady |
| **A4** irrigation | un-validatable (no irrigated CFT patch) | `QIRRIG_DEMAND` non-zero, varies year to year |
| **A2** fire crop branch | `≡ 0` without a crop CFT | `BAF_CROP` non-zero (1.18e-11) |

Note `NFIRE` is non-zero only in 2016 while `BAF_CROP` is non-zero only from
2017 — the crop-CFT fire branch and the natural-veg branch are separately
exercised, which is what the A2 crop/non-crop split needs.

**Caveat on interpretation:** these are 4 cold-start years, not a spun-up CN
equilibrium; `TOTVEGC` is still rising (102 → 451). They establish that the
crop *code path* runs and produces non-degenerate state — which is what
single-step parity needs. They are not a claim about equilibrium pool sizes,
and (per the forcing substitution above) not a claim about this location's
agronomy.

### Two more generator defects the run itself surfaced

Both aborted *after* initialization, so only a real run could find them:

4. **`fd.yaml` (the NUOPC field dictionary) was never copied** — aborts inside
   `ESMF_IO_YAMLRead` with a bare "Caught exception reading content from file",
   naming no file.
5. **A global `<datafiles>` edit injected the met files into `topo.observed`**,
   which then tried to read `TOPO` out of a `clmforc` file. datm reported only
   "NetCDF: Variable not found" — no file, no variable name; the PIO stack
   frame was the only handle. Stream patching is now per-`<stream_info>`.

Plus a sixth, from the history stream: `GDDPLANT`, `GRAINC`, `GRAINN`,
`CROPPROD1C`, `GRAINC_TO_FOOD` and `QIRRIG` are **not** CTSM history field
names (`QIRRIG_DEMAND` is; the grain-C fields are not registered under those
literals). An unknown name aborts at `histFileMod.F90:887` after
initialization. The generator now validates `hist_fincl1` against the names
CTSM actually registers and drops unknowns with a warning.

## Instrumented crop reference dumps GENERATED (two windows)

`bgcdumpMod.F90` was extended with the crop lifecycle + crop N fields (see
`scripts/validation/fortran_sourcemods/bgcdumpMod.F90`, all behind
`if (use_crop)` so non-crop references stay byte-identical) and the case
rebuilt — an incremental `./case.build`, ~9 s total; verify it took by checking
`strings bld/cesm.exe | grep CROPLIVE`, since CIME reports SUCCESS even when
nothing recompiled.

Dumps are produced by the proven two-stage pattern: advance to the window with
`PDUMP/BGCDUMP_NSTEP_*` closed, write a restart, then run a short bracket from
it with the window wide open.

| Ref dir | Date | nsteps | Files | What it captures |
|---|---|---|---|---|
| `refs_crop_may/` | 2020-05-30 | 38665–38676 | 48 | **the fertilization window** |
| `refs_crop_july/` | 2020-07-19 | 39865–39876 | 48 | mid-season growth / grain fill |

(Both under `clm_bgc_spinup/crop_ref_usplains/`; 12 `bgcdump_after_fire` +
36 `pdump` per window.)

### What the reference actually contains

Patch indices **5–12 are the eight live crop patches** (`CROPLIVE = 1`,
`NYRS_CROP_ACTIVE = 4`); everything else is `1e36`/0 fill, which is the correct
non-crop signature.

**May window (n38665)** — CPHASE = 2 (leaf emergence) on all eight:

| Field | Values across the 8 crop patches |
|---|---|
| `HUI` | 131.99, 131.85, 342.10, 342.03, 64.93, 63.19, 162.45, 162.29 |
| **`FERT`** | **9.693e-6, 9.693e-6, 5.727e-6, 5.779e-6, 1.567e-6, 1.567e-6, 9.693e-6, 9.693e-6** |
| `FERT_TO_SMINN` | 8 non-zero columns, max 9.693e-6 |
| `SOYFIXN` | 0 |

**July window (n39865)** — CPHASE 2 and 3 (grain fill), `HUI` 552–1274,
`FERTNITRO` 0.71–14.75, but `FERT = 0`: **the fertilization window has already
closed by mid-July.** That is why the May window exists, and it is a concrete
reminder that "the crop is growing" and "the crop-N flux is live" are different
windows — a mid-season dump alone would have shown `FERT ≡ 0` and looked like a
dead path.

`SOYFIXN` is 0 in **both** windows. CTSM's `CNSoyfix` requires a later growth
phase than either window sampled, so the soybean-fixation reference is still
outstanding — a third window (or a phase-targeted bracket) is needed before
`n_soyfix!` can be validated. Recorded as a known gap rather than papered over.

## Wiring `n_fert!` — it would be a NO-OP today. Do not wire it yet.

The task brief scoped "wire `n_fert!`/`n_soyfix!` now that a reference exists".
Reading the two sides against each other says **the reference is not the
binding constraint for `n_fert!`, and wiring it would create a green-looking
dead path** — the exact bug class in `dead-initcold-systemic` ("assert the
wiring **and a non-no-op body**").

**What `n_fert!` actually does.** `src/biogeochem/n_dynamics.jl:332` is a pure
aggregation: it zeroes `soilbgc_nf.fert_to_sminn_col` and scatters
`cnveg_nf.fert_patch` into it, weighted by `patch.wtcol`. It computes no
fertilizer. It is a faithful port of `CNNFert`, which is also only a `p2c`.

**`fert_patch` is never computed in CLM.jl.** Every occurrence in `src/` is a
declaration, an allocation, or a zeroing:

| Site | What |
|---|---|
| `types/cn_veg_nitrogen_flux.jl:172` | field declaration |
| `types/cn_veg_nitrogen_flux.jl:548` | `nanvec(np)` allocation |
| `types/cn_veg_nitrogen_flux.jl:1111` | `= 0.0` (zero-flux reset) |
| `biogeochem/n_dynamics.jl:341` | **read** by `n_fert!` |

**Where CTSM computes it:** `CNPhenologyMod.F90:2529`, *inside*
`CropPhenology`:

```fortran
fert(p) = (manunitro(ivt(p)) * 1000._r8 + fertnitro(p)) / fert_counter(p)
```

i.e. fertilizer is applied over a counted window after planting, from the
`fertnitro_patch` surfdata/param input and the manure term. That is precisely
the part of `CropPhenology` that `crop_phenology!` (`phenology.jl:2344`) does
**not** port — it only sets `cphase` and the grain-fill `bglfr`.

**So the dependency order is:** port the fertilization window inside crop
phenology (which needs `plant_crop!` driven, so a planting date exists to count
from) → *then* `fert_patch` becomes non-zero → *then* wiring `n_fert!` moves
real nitrogen. Wiring `n_fert!` first would scatter an all-zero array, pass
any "is it wired?" test, and change nothing — a vacuous green.

`n_soyfix!` is different: it computes `soyfixn_patch` itself from soil
moisture / growth phase / mineral N, and takes `crop::CropData` directly. It
is not blocked on `fert_patch`. It *is* blocked on the crop lifecycle state
(`cphase`, `croplive`) being driven, which is the same `plant_crop!` gap.

**Conclusion:** the honest next step is not "wire the two N routines" but
"finish `crop_phenology!` against the Fortran reference — planting decision,
vernalization, phase transitions, the fertilization window, harvest — and wire
the N routines behind it". The reference this task produced is what makes that
verifiable; the reference alone does not make the N wiring meaningful. Not
wired here, deliberately, and for a sharper reason than #218's original one.

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
