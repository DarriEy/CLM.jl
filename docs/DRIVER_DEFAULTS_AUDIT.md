# Driver defaults audit — CLM.jl entry points vs CTSM namelist defaults

Systematic sweep of every control flag / parameter default exposed by CLM.jl's
driver entry points (`clm_initialize!`, `clm_run!`, `CLMDriverConfig`) and the
global control state (`src/constants/varctl.jl`), checked against CTSM's
authoritative defaults.

## Why this audit exists

Two independently-found, independently-merged bugs had the *same shape*:

| PR | Flag | CLM.jl default | CTSM default | Consequence |
|----|------|----------------|--------------|-------------|
| #252 | `use_bedrock` | `true` unconditionally | **conditional**: `.false.` under `use_fates`/`vichydro`/`clm4_5`, `.true.` otherwise | `clamp_zwt_to_bedrock!` dragged zwt 8.6→2.28 m; `qflx_drain ~ exp(-zwt/hkdepth)` ignited 3–17 mm/day drainage; top-layer soil moisture 0.12 vs 0.32, starving FATES 22% by day 19 |
| #225 | `h2osfcflag` | `0` | `1` (`SoilHydrologyType.F90`, `clm_soilhydrology_inparm`) | surface-water store silently disabled; infiltration excess went straight to surface runoff instead of ponding + re-infiltrating; `frac_h2osfc` pinned at 0 |

Both were invisible for the same reason: **the harnesses that mattered passed
explicit values, so the wrong default only bit callers who trusted it.** A wrong
default is not caught by any conservation check, any parity scorecard that
overrides it, or any test that constructs its own config.

The diagnostic signature of this bug class is therefore not a failing test — it
is **a workaround**: a harness passing a value explicitly that it should have
been able to inherit.

## Authoritative sources

CTSM tree: `/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/installs/clm/`

1. `bld/namelist_files/namelist_defaults_ctsm.xml` — the real namelist defaults.
   **Many are CONDITIONAL** on `phys` / `use_fates` / `vichydro` /
   `configuration` / `structure`. Capturing the conditions is the whole point:
   that is exactly what #252 got wrong.
2. `bld/namelist_files/namelist_definition_ctsm.xml` — types and valid values.
3. `*_inparm` blocks and module-level initialisers in `src/` (e.g.
   `clm_varctl.F90`, `controlMod.F90`, `SoilHydrologyType.F90`,
   `SoilWaterMovementMod.F90`) — the **code fallback**, which can *disagree*
   with the namelist default. That disagreement was the subtlety in #252
   (code fallback `.false.`, namelist default `.true.` for CLM5).

**Resolution conditions for this port.** CLM.jl targets `phys="clm5_0"`,
`use_fates=.false.`, `configuration="clm"`, `structure="standard"`,
`vichydro=0`. All "CTSM default" values below are resolved under those
conditions, with the conditional variants recorded.

---

## Headline finding: the default driver configuration is CLM4.5-shaped, not CLM5

The single most consequential result of this sweep is not one flag but a
**coherent pattern across four of them**. CLM.jl's driver defaults select the
CLM4.5 soil-hydrology and photosynthesis configuration, while CTSM's `clm5_0`
defaults select a different one. The multi-biome parity harness
(`scripts/parity_run_domain.jl`) already hand-rebuilds the CLM5 configuration
flag by flag — and says so in its own header comment:

```
# All three reference runs share the Bow lnd_in physics config (use_hydrstress +
# use_luna, Jackson1996 rooting, Vionnet2012 + wind-dependent snow density,
# baseflow_scalar/int_snow_max per lnd_in), so the only per-domain differences are
...
    use_aquifer_layer = get(cfg, :aquifer, false), use_hydrstress = true, use_luna = true,
    h2osfcflag = 1,  # CLM5 default (surface water active); matters for wet sites
```

That block is four workarounds in four consecutive lines. `h2osfcflag` was
already recognised as one and fixed in #225. The other three are the same bug,
unfixed.

### M1 — `use_aquifer_layer` — **CLOSED** (was: MISMATCH, live, highest severity)

> **Resolved.** Default flipped to `false` in `clm_initialize!`, `clm_run!` and
> `CLMDriverConfig`, matching CTSM's derivation for `clm5_0`. Verified: the flip
> reaches the solver switch (`clm_driver.jl` `swm_cfg`) and now selects
> `MOISTURE_FORM` + `BC_ZERO_FLUX`, so `use_aquifer_layer(cfg)` derives `false` —
> exactly CTSM's chain. Full suite `--check-bounds=yes` is **byte-identical to
> baseline: 26150 passed / 3 broken / 0 failed / 0 errored** (baseline confirmed
> by two independent runs). That is a NO-REGRESSION result, not proof the physics
> is unchanged — it cannot be, the solver differs. Nothing moved because the
> harnesses sensitive to the lower boundary condition pass the flag explicitly
> (`longhorizon_conservation.jl` via a const; `parity_run_domain.jl:195` via
> `get(cfg, :aquifer, false)` — note the scorecard *already* defaulted to
> `false`, independent corroboration), and the suites that do inherit the default
> assert on canopy/energy quantities over 1–4 timesteps. Pinned in
> `test/test_driver_defaults_audit.jl`.


CTSM does not expose `use_aquifer_layer` as a namelist flag at all. It is a
**derived function** (`src/biogeophys/SoilWaterMovementMod.F90:221`):

```fortran
function use_aquifer_layer() result(lres)
   if(lower_boundary_condition == bc_aquifer .or. lower_boundary_condition == bc_watertable)then
      lres=.true.
   else
      lres=.false.
   endif
```

with `bc_zero_flux = 2`, `bc_waterTable = 3`, `bc_aquifer = 4`
(`SoilWaterMovementMod.F90:55-59`). Resolving the chain for `clm5_0`:

```xml
<soilwater_movement_method               >1</soilwater_movement_method>   <!-- clm5_0 -->
<soilwater_movement_method phys="clm4_5" >0</soilwater_movement_method>
<lower_boundary_condition soilwater_movement_method="0"                      >4</lower_boundary_condition>
<lower_boundary_condition soilwater_movement_method="1" use_bedrock=".true." >2</lower_boundary_condition>
```

`clm5_0` → method 1 → `lower_boundary_condition = 2 = bc_zero_flux` →
**`use_aquifer_layer() = .false.`**

CLM.jl defaults it to **`true`** in all three places
(`clm_initialize!`, `clm_run!`, `CLMDriverConfig`). That flag is fully live —
it selects the entire soil-water solver (`clm_driver.jl:1884-1886`):

```
#   use_aquifer_layer=true  → Zeng-Decker 2009 with BC_AQUIFER (default)
#   use_aquifer_layer=false → Moisture form with BC_ZERO_FLUX
```

So the default driver runs **the CLM4.5 Zeng-Decker-2009 solver with an aquifer
lower boundary**, plus the `use_aquifer_layer`-gated aquifer-recharge branch at
`clm_driver.jl:1908` and the `QCHARGE` diagnostic — where CTSM `clm5_0` runs the
moisture-form solver against a zero-flux bedrock boundary with no aquifer.

**CTSM would refuse to start in CLM.jl's default configuration.** After #252,
CLM.jl defaults `use_bedrock=true`; combined with `use_aquifer_layer=true` that
is `use_bedrock .and. lbc == bc_aquifer`, which `SoilWaterMovementMod.F90:181`
explicitly aborts on:

```fortran
if((use_bedrock) .and. (lower_boundary_condition /= bc_zero_flux)) then
   call endrun(subname // ':: ERROR inconsistent soilwater_movement namelist: use_bedrock requires bc_zero_flux lbc')
```

This is the strongest possible evidence that the default pair is wrong: it is a
combination the reference model treats as a fatal namelist inconsistency.

**Who is already working around it.** Essentially every harness that runs real
physics passes `use_aquifer_layer=false` explicitly — `parity_run_domain.jl`,
`parity_run_domain_gpu.jl`, `parity_run_spinup.jl`, `fortran_parity_validate.jl`,
`fortran_parity_snow.jl`, `fortran_parity_cn_coldstart.jl`,
`initcold_bow_bitcheck.jl`, `multisite_smoke.jl`, `lake_water_balance.jl`,
`spring_spinup_test.jl`, `validate_luna_coupled.jl`, `probe_deep_soil.jl`,
`probe_h2osfc_subdaily.jl`, `probe_veg_pool_audit.jl`,
`probe_cbalance_decompose.jl`, `maritime_snow_step_dump.jl`. The default is
believed by nobody who has checked it against Fortran.

#### M1 RE-CONFIRMED on the fixed tridiagonal solver

#259 measured the `false` flip as byte-identical across the full suite. That
result was valid but rested on a compromised foundation: **both sides of it were
solving soil water with a broken solver.** #263 subsequently found that
`_tridiag_multi_kernel!` wrote `cp[col,jj]` and re-read `cp[col,jj-1]` on the
next iteration, a store→load recurrence LLVM reorders at `-O2`, so
`tridiagonal_multi!` returned wrong values on the CPU for any column with
`nlevs >= 5`. `soil_water_movement.jl:1709` solves ~21 levels — and it is
precisely the solver `use_aquifer_layer` selects between. The comparison has
therefore been re-run on the corrected solver.

**It still holds, and the evidence is now stronger than "nothing moved".**
Re-running the 240-day Bow trajectory on the fixed solver with
`use_aquifer_layer=true` against the `false` default:

| | `false` (default) | `true` |
|---|---|---|
| `qflx_drain` (mm / 240 d) | 0.000000 | **8 748 338.42** |
| total runoff (mm) | 377.032388 | **8 748 715.40** |
| `qflx_surf` (mm) | 377.032388 | 376.959705 |
| ET (mm) | 114.162495 | 112.046123 |
| soil water (mm) | 284.044305 | 286.216041 |
| precipitation (mm) | 404.066486 | 404.066486 |

`true` returns **8.7 million mm of drainage from a site that received 404 mm of
precipitation** — a runaway of about four orders of magnitude beyond the entire
water input. So the two settings are emphatically *not* equivalent physics; the
byte-identical suite is a statement about test coverage, nothing else.

Note what did *not* happen: the run completed. Soil water, ET and surface runoff
stayed plausible, and no balance check fired, because the runaway lives in a
drainage term that is not debited against column storage. A conservation check
cannot see this — same lesson as M2's partition shift, one register louder.

Two honest caveats. First, this arm also carries `use_bedrock=true`, which is
exactly the pair `SoilWaterMovementMod.F90:181` calls a fatal namelist
inconsistency — so this measures the forbidden combination the port used to
ship, which is the point, rather than the ZD09 solver in isolation. Second,
CLM.jl does **not** currently mirror that `endrun`; `clm_initialize!` accepts
`use_bedrock=true, use_aquifer_layer=true` silently. Porting the guard would
have turned this runaway into a clear error message and is worth doing, but it
is a behaviour change beyond this PR's scope.

**Verdict: `use_aquifer_layer=false` is re-confirmed on the corrected solver.**

### M2 — `baseflow_scalar` (MISMATCH, live, 10×)

```xml
<baseflow_scalar                lower_boundary_condition="2">0.001d00</baseflow_scalar>
<baseflow_scalar                lower_boundary_condition="1">1.d-2</baseflow_scalar>
<baseflow_scalar phys="clm4_5"  lower_boundary_condition="2">1.d-2</baseflow_scalar>
```

For `clm5_0` the resolved `lower_boundary_condition` is `2` (see M1), so the
CTSM default is **`0.001`**. `clm_run!` defaults **`1.0e-2`** — the value that
belongs to `lbc=1`, or to `clm4_5`. A **10× multiplier on the baseflow /
drainage rate**, consistent with the same CLM4.5-shaped drift as M1.

Note CLM.jl already disagrees with *itself* here: the calibration default in
`src/calibration/parameters.jl:72` is `0.001` (the CTSM value), while the driver
entry point and `soil_hydrology.jl:66,77` use `1.0e-2`.

Live: `BASEFLOW_SCALAR[]` feeds the drainage computation in `soil_hydrology.jl`.

### M3 — `use_hydrstress` (CONDITIONAL-MISMATCH, live)

```xml
<use_hydrstress                                                       >.false.</use_hydrstress>
<use_hydrstress                use_fates=".false." configuration="clm">.true.</use_hydrstress>
<use_hydrstress phys="clm4_5"  use_fates=".false." configuration="clm">.false.</use_hydrstress>
```

CTSM `clm5_0` / non-FATES / `configuration="clm"` → **`.true.`** (plant hydraulic
stress is the CLM5 default photosynthesis path). CLM.jl defaulted **`false`** in
`clm_initialize!`, `clm_run!` and `CLMDriverConfig` — CTSM's *code fallback*
(`clm_varctl.F90`), not its namelist default: the same root cause as
#252/#259/#267. Consequence: the default driver used the Soil-Moisture-Stress
(BTRAN) path rather than PHS, changing stomatal conductance, transpiration and
GPP under water stress.

**FIXED in the PHS campaign** — `use_hydrstress` is now `Union{Bool,Nothing}`
with `nothing` → `!use_fates`, matching `use_bedrock` (#252) and `use_luna`
(#267). The default has to be CONDITIONAL rather than a bare `true` because
`control.jl:102` endruns on PHS+FATES, mirroring CTSM: the `.true.` namelist
default applies to the non-FATES configuration only.

`varctl.use_hydrstress` keeps its `false` **struct** default, which mirrors
CTSM's code fallback; it is the DERIVATION that carries the namelist default
(the #265 principle, as for `create_crop_landunit`).

See `docs/PHS_DEFAULT_BLOCKERS.md` for the blocker history and the measured
suite diff.

### M4 — `use_luna` (CONDITIONAL-MISMATCH, live)

```xml
<use_luna                                   >.true.</use_luna>
<use_luna                use_fates=".true." >.false.</use_luna>
<use_luna phys="clm4_5"  use_fates=".false.">.false.</use_luna>
```

CTSM `clm5_0` / non-FATES → **`.true.`**; the default is CONDITIONAL on
`use_fates` (CTSM `endrun`s on LUNA+FATES, `controlMod.F90:505`). CLM.jl
defaulted **`false`** — CTSM's *code fallback* (`clm_varctl.F90:371`), not its
namelist default: the same root cause as #252/#259.

**FIXED in the LUNA campaign** — `use_luna` is now `Union{Bool,Nothing}` with
`nothing` → `!use_fates`, matching `use_bedrock`. Full suite unchanged
(26351/0/0/3 before and after; zero movers).

**But the flip is currently INERT on the default path.** See
`docs/LUNA_DEFAULT_CAMPAIGN.md`: CTSM applies the LUNA-acclimated
`vcmx25_z`/`jmx25_z` in **both** photosynthesis routines —
`Photosynthesis` (non-PHS) at `PhotosynthesisMod.F90:1721,1748` and
`PhotosynthesisHydraulicStress` at `:3335,3377`. CLM.jl ported only the PHS
pair (`photosynthesis.jl:2100,2132`). The default, non-PHS `photosynthesis!`
accepts a `use_luna` kwarg and **never reads it**, so LUNA acclimates
`vcmx25_z` (measured: 30.0 → 38.1/43.9 over 16 days at Bow) and the result is
discarded. That is why the flip moved zero tests, and it is why closing M4
does *not* by itself deliver CLM5 photosynthesis.

**The LUNA consumption guard (verified against CTSM source, 2026-07-19).** CTSM
gates every LUNA read with *three* conditions, not two — the same expression at
all four consumption sites:

```fortran
if(use_luna.and.c3flag(p).and.crop(patch%itype(p))== 0)
```

- `src/biogeophys/PhotosynthesisMod.F90:1721` — non-PHS, `lmr25` (leaf resp).
- `src/biogeophys/PhotosynthesisMod.F90:1748` — non-PHS, `vcmax25`/`jmax25`/`tpu25`.
- `src/biogeophys/PhotosynthesisMod.F90:3335` — PHS, `lmr25_sun`/`lmr25_sha`.
- `src/biogeophys/PhotosynthesisMod.F90:3377` — PHS, `vcmax25_sun`/`_sha` etc.

The two PHS sites guard **identically** to the two non-PHS sites (byte-identical
condition text apart from whitespace). `crop` is the `pftcon` crop flag indexed
by `patch%itype(p)` — a 0/1-valued per-PFT parameter, so `crop(ivt)==0` is
"this PFT is not a crop". LUNA therefore never acclimates a crop patch: crops
stay on the static `vcmax25top*nscaler` profile in CTSM.

`crop(patch%itype(p))` also appears at `:3498`/`:3515`/`:3718`/`:3753`, but those
are the separate `modifyphoto_and_lmr_forcrop` branches and are unrelated to LUNA.

**Status:** #268 ported the non-PHS pair (`:1721`/`:1748`) CTSM-exact, crop guard
included. The PHS pair was ported earlier gating on `use_luna && c3flag_p` only
— the `crop(itype)==0` condition was **missing**, so crop patches under PHS
received LUNA-acclimated `vcmax25`/`jmax25`/`lmr25`. Fixed in the PR that added
this note.

### M5 — `create_crop_landunit` (FIXED — was MISMATCH, live)

```xml
<create_crop_landunit     use_fates=".false.">.true.</create_crop_landunit>
<create_crop_landunit     use_fates=".true." >.false.</create_crop_landunit>
```

CTSM keys this on **`use_fates`**, not on `use_crop`: for any non-FATES run it is
`.true.`. `clm_initialize.jl:138` used to set
`varctl.create_crop_landunit = use_crop`, so a non-FATES, non-crop run got
`.false.` where CTSM gives `.true.`. Live (4 reads in `src/`) — it affects
subgrid landunit construction.

**Fixed:** `clm_initialize.jl` now derives `!use_fates`, and `varpar_init!`
(`varpar.jl`) branches on `varctl.create_crop_landunit` rather than
`varctl.use_crop`, mirroring `clm_varpar.F90:209`.

The `VarCtl` struct field default stays `false` — that is *correct*, mirroring
CTSM's code fallback `clm_varctl.F90:154`. Only the derivation (CTSM's
`CLMBuildNamelist` step) and the `varpar_init!` branch were wrong. A consequence
worth stating: tests that call `varpar_init!` directly without setting the flag
keep the old `cft_size = 0` behaviour, so the change is confined to the driver
path.

**Verified against CTSM source (2026-07-19).** The default is not merely a
namelist preference — for a non-FATES run CTSM makes `.false.` a *fatal build
error*:

- `bld/namelist_files/namelist_defaults_ctsm.xml:2377-2378` — keyed on
  `use_fates` alone; `use_crop` is not a selector.
- `bld/CLMBuildNamelist.pm:2248-2250` —
  `"$var is false which is ONLY allowed when FATES is being used"` (fatal).
  So no non-FATES CTSM run can ever have `create_crop_landunit=.false.`
- `src/main/controlMod.F90:460-461` — the runtime check is only the
  *converse* (`use_crop .and. .not. create_crop_landunit` → `endrun`); it does
  **not** permit the non-crop case to turn the landunit off.
- `src/main/clm_varpar.F90:209` — the `cft_size` branch keys on
  `create_crop_landunit`, and its else-branch is commented
  `"only true when FATES is active"`.

The port's mirror of that branch is `src/constants/varpar.jl:178`, which keys on
`varctl.use_crop`. Consequence for the **default** (non-crop, non-FATES) run:
CTSM takes the `create_crop_landunit=.true.` branch (`cft_size = surf_numcft`,
crop landunit built), CLM.jl takes the else-branch (`cft_size = 0`, crop area
folded into natural veg). Two different subgrid structures from one surfdata.

Note the downstream consumers (`init_gridcells.jl:144`,
`surfdata.jl:442,734`) already gate on `create_crop_landunit` correctly — only
the *derivation* and the `varpar_init!` branch are wrong.

#### Measured blast radius

`scripts/probe_crop_landunit_default.jl` builds the real `initGridCells!`
subgrid both ways for a default non-crop, non-FATES run. The answer depends
entirely on whether the surfdata has crop area:

| surfdata | `PCT_CROP` | subgrid change |
|---|---|---|
| `domain_Bow_at_Banff_lumped` (and every scorecard domain checked: Stillwater, Donga, HubbardBrook, Tagus, Cropland_Mead) | `0` | **none** — `nl/nc/np = 2/2/4` both ways, identical landunit types and weights |
| `crop_cft_surfdata/surfdata_cropCFT_USplains_1pt.nc` | `45.01%` | **structural** — `nl/nc/np` `4/12/16` → `5/16/20`; `wt_lunit[ISTSOIL]` `0.9628` → `0.5127` with a new `ISTCROP` landunit at `0.4501`, carrying 4 crop patches `[17, 19, 23, 75]` (corn / spring wheat / soybean / generic, rainfed-only at `irrigate=false`) |

So the divergence is **latent on every domain this port currently validates**:
`set_landunit_crop_noncompete!` returns early on
`wt_lunit[g,ISTCROP] > 0.0`, so with `PCT_CROP = 0` no crop landunit is built
either way. It becomes structural the moment a surfdata with crop area is used —
where today's default silently folds 45% of the gridcell into natural veg.

What *does* change on every non-FATES run is the varpar bounds, since
`cft_size` goes `0 → surf_numcft` (2 on the generic-crop Bow file, 64 on the
full-CFT file), moving `cft_ub` and `maxveg` (`15 → 16` on Bow).

#### Latent inconsistency found while probing

`count_subgrid_elements` and `set_landunit_crop_noncompete!` disagree when
`create_crop_landunit=true` **and** `cft_size==0`: the counter's `n_cft == 0`
fallback reserves a column and patch (`surfdata.jl:740-747`) while the builder
returns early (`init_gridcells.jl:147`), so `initGridCells!` dies with
`landunit count mismatch: 4 != 5`. That state is unreachable once the flag and
`varpar_init!` move together (which is the fix), but it is a live landmine for
anyone who sets the flag alone. CTSM forbids the state outright —
`surfrdMod.F90:958-961` `endrun`s on `cft_size == 0 .and. any(PCT_CROP > 0)`.

---

## Full sweep

Legend: **MATCH** — agrees with CTSM under this port's conditions.
**MISMATCH** — disagrees. **COND-MISMATCH** — agrees under some conditions,
disagrees under this port's. **INERT** — the field is never read in `src/`, so a
mismatch is documentation-only until something starts reading it (recorded
anyway: an inert wrong default is a landmine for the first consumer).

### Driver entry-point keywords

| Flag | CLM.jl default | CTSM default (+ conditions) | Verdict | Consequence if wrong | Explicit passers (workarounds) |
|---|---|---|---|---|---|
| `use_bedrock` | `nothing` → `!use_fates` | `.false.` if `use_fates`/`vichydro`/`clm4_5`; else `.true.` | **MATCH** (fixed #252) | zwt clamped to bedrock → runaway drainage | — |
| `h2osfcflag` | `1` | `1` (code fallback `SoilHydrologyType.F90:339,347`) | **MATCH** (fixed #225) | surface-water store disabled | `parity_run_domain.jl` (now redundant) |
| `use_aquifer_layer` | ~~`true`~~ → **`false`** | **`.false.`** (derived: `clm5_0`→method 1→lbc 2 = `bc_zero_flux`) | **CLOSED** (#259; re-confirmed on the fixed tridiagonal solver — see "M1 RE-CONFIRMED") | selected CLM4.5 Zeng-Decker-2009 + aquifer solver instead of CLM5 moisture-form + zero-flux; CTSM `endrun`s on this combined with `use_bedrock=true`. Measured on the corrected solver: `true` returns 8.7e6 mm of drainage from a 404 mm-of-precipitation site | 16 scripts (see M1) |
| `baseflow_scalar` | ~~`1.0e-2`~~ → **`0.001`** | **`0.001`** (`lbc=2`); `1.d-2` only for `lbc=1` or `clm4_5` | **CLOSED** (see "M2 CLOSED" below) | was 10× on the drainage/baseflow rate; measured −7.6 % baseflow / +0.57 % surface runoff at MerBleue, exactly zero at Bow (branch never entered, in CTSM either) | `parity_run_domain.jl` (per-domain; 18/20 already `0.001`) |
| `use_hydrstress` | `nothing` → `!use_fates` (**FIXED**) | **`.true.`** (`clm5_0`, non-FATES, `configuration="clm"`) | **COND-MISMATCH — CLOSED** | was: PHS off → BTRAN path; different gs/transpiration/GPP under stress | `parity_run_domain.jl`, `parity_run_domain_gpu.jl`, `probe_h2osfc_subdaily.jl` |
| `use_luna` | `nothing` → `!use_fates` | **`.true.`** (`clm5_0`, non-FATES) | **FIXED** (this campaign) | — but see the LUNA-consumption gap below: the flip is currently INERT on the default non-PHS path | `parity_run_domain.jl`, `parity_run_domain_gpu.jl`, `fortran_parity_cn_coldstart.jl` (now redundant) |
| `int_snow_max` | `2000.0` | `2000.` (`1.e30` for `clm4_5`) | **MATCH** | — | — |
| `dtime` | `1800` | `1800` | **MATCH** | — | — |
| `use_cn` | `false` | `.true.` for BGC compsets, `.false.` for SP | **MATCH** (compset-level, not a physics default) | — | — |
| `use_crop` | `false` | compset-level | **MATCH** | — | — |
| `use_fates` | `false` | compset-level | **MATCH** | — | — |
| `use_lch4` | `false` | `.true.` unless `soil_decomp_method="None"` | **COND-MISMATCH** (CN runs only) | methane biogeochem off by default in CN runs | — |
| `use_c13` / `use_c14` | `false` | `.false.` | **MATCH** | — | — |
| `cnfire_method` | `:nofire` | `li2016crufrc` for CN | **COND-MISMATCH** (documented deliberate: fire needs stream files) | — | — |
| `soil_layerstruct` | `"20SL_8.5m"` | `20SL_8.5m` (`clm5_0`) | **MATCH** | — | — |
| `all_active` | `false` | `.false.` | **MATCH** | — | — |
| `ncopies` / `read_full_grid` / `verbose` / file paths | — | CLM.jl-only plumbing | **NOT-IN-CTSM** | — | — |

### `varctl` global control state

79 fields swept. "reads" = `grep -rn "varctl.<field>" src/ | wc -l`; a field with
0 reads is **INERT** (mismatch is documentation-only *today*, but is a landmine
for the first consumer — recorded rather than dismissed).

A structural observation that explains most of the mismatches below: **where
CTSM's namelist default and its Fortran code fallback disagree, CLM.jl copied
the code fallback every single time.** That is precisely the #252 failure mode,
and it recurs in 14 fields: `co2_ppmv`, `convert_ocean_to_land`,
`n_dom_landunits`, `n_dom_pfts`, the six `toosmall_*`,
`downscale_hillslope_meteorology`, `nyr_forcing`, `nsegspc`,
`glc_snow_persistence_max_days`. The port faithfully mirrored
`clm_varctl.F90`'s initialisers — but those initialisers are the value you get
when build-namelist emits *nothing*, which never happens in a real CTSM run.

| Field | CLM.jl | CTSM (clm5_0, non-FATES) + conditions / code fallback | Verdict | Reads |
|---|---|---|---|---|
| `convert_ocean_to_land` | `false` | **`.true.`** (`defaults:2382`, unconditional). Code fallback `.false.` — **disagrees** | **MISMATCH (live)** | 1 |
| `co2_ppmv` | `355.0` | **`379.0`** at default `sim_year=2000` (`defaults:25`). 1850/PtVg 284.7, 1979 336.6, 2010 388.8, 2015 397.5, 2018 408.83. Code fallback `355.0` — **disagrees** | **MISMATCH** | 0 (but threaded as a plain arg through 19 sites) |
| `glc_snow_persistence_max_days` | `7300` | **`0`**; `7300` is the **clm4_5** variant (`defaults:540-541`). Code fallback 7300 — disagrees | **MISMATCH (live)** | 0 in varctl, but the live kwarg at `glacier_surface_mass_balance.jl:146` also defaults 7300 |
| `nsegspc` | `20` | **`35`** (`defaults:2404`). Code fallback 20; CTSM's own *definition* docs say 20 (stale) | MISMATCH (INERT) | 0 |
| `nyr_forcing` | `10` | **`1`** (`defaults:688`); 20 under matrix spinup. Code fallback 10 | MISMATCH (INERT) | 0 |
| `h2osno_max` | `-999.0` | **`10000.0`** (`defaults:465`); fast 5000, clm4_5 1000. No code fallback — `controlMod.F90:557` **hard-errors if `<= 0`** | MISMATCH (INERT — physics correctly uses `varcon.H2OSNO_MAX = 10000.0`; `varctl` copy is a dead duplicate) | 0 |
| `n_dom_landunits` | `-1` | **`0`** (`defaults:2387`); fast 1. Code fallback −1 | MISMATCH (inert: guards are `> 0`, so −1 ≡ 0) | 2 |
| `n_dom_pfts` | `-1` | **`0`** (`defaults:2390`) | MISMATCH (INERT) | 0 |
| `toosmall_soil`/`_crop`/`_glacier`/`_lake`/`_wetland`/`_urban` | `-1.0` | **`0.d00`** (`defaults:2393-2398`). Code fallback −1 | MISMATCH ×6 (inert: guards `> 0`) | 2 each |
| `downscale_hillslope_meteorology` | `false` | **`.true.`** (`defaults:665`, unconditional). Code fallback `.false.` | MISMATCH (INERT — hillslope off) | 0 |
| `z0param_method` | `""` | `ZengWang2007`; `Meier2022` for clm5_1/clm6_0. No code fallback (namelist mandatory) | MISMATCH (cosmetic — `friction_velocity.jl:322` branches only on `== "Meier2022"`, so `""` *is* ZengWang2007) | 0 |
| `irrigate` | `false` | `.false.` at `sim_year_range=constant`; **`.true.`** for clm5_0 transient `1850-2100` non-CNDV | MATCH (constant); COND for transient | 1 |
| `snow_cover_fraction_method` | `""` | `SwensonLawrence2012` | MATCH in effect (`clm_initialize.jl:382-383` fills it) | 2 |
| `use_biomass_heat_storage` | `false` | `.false.` for **clm5_0** (`.true.` generic/clm5_1/clm6_0) | MATCH | — |
| `spinup_state` | `0` | `0` at `clm_accelerated_spinup=off` | MATCH | — |
| `use_c13`, `use_c14`, `use_cndv`, `use_cn`, `use_noio`, `use_vichydro`, `use_extralakelayers`, `use_excess_ice`, `use_snicar_frc`, `use_SSRE`, `use_soil_moisture_streams`, `use_lai_streams`, `use_z0m_snowmelt`, `hist_wrtch4diag`, `glc_do_dynglacier`, `use_hillslope`, `use_hillslope_routing`, `hillslope_fsat_equals_zero`, `run_zero_weight_urban`, `collapse_urban`, `crop_fsat_equals_zero`, `crop_residue_removal_frac`, `all_active`, `use_subgrid_fluxes`, `spinup_matrixcn`, `hist_wrt_matrixcn_diag`, `nyr_SASU`, `iloop_avg`, `o3_veg_stress_method`, `co2_type`, `reduce_dayl_factor` | as declared | agree under this port's conditions | **MATCH** ×31 | — |
| `snicar_solarspec`, `snicar_dust_optics`, `snicar_use_aerosol`, `snicar_snobc_intmix`, `snicar_snodst_intmix`, `do_sno_oc`, `snicar_numrad_snw`, `snicar_snw_shape` | as declared | agree for clm5_0 (`snw_shape=sphere` is the clm5_0/clm4_5 variant; `hexagonal_plate` is clm5_1+ — already documented in-source) | **MATCH** ×8 | — |
| `outnc_large_files`, `ndep_from_cpl`, `bound_h2osoi`, `nhillslope`, `max_columns_hillslope`, `o3_ppbv`, `anoxia`, `nfix_timeconst` | as declared | **not namelist variables in CTSM** — code fallback is the only authority, and CLM.jl matches all | **MATCH** (NOT-IN-CTSM) ×8 | — |

### M6 — the `use_flexibleCN` cascade (CONDITIONAL-MISMATCH cluster, applies to every CN run)

The largest *cluster* found. Under `phys=clm5_0` **with `use_cn=.true.`**, CTSM
turns on `use_flexibleCN`:

```xml
<use_flexibleCN                    >.false.</use_flexibleCN>
<use_flexibleCN use_cn=".true."    >.true.</use_flexibleCN>   <!-- clm4_5+CN → .false. -->
```

which then cascades into seven further defaults, all of which CLM.jl hardcodes
at their SP-mode (`false`/`0`) values:

| Flag | CLM.jl | CTSM under `clm5_0` + `use_cn` | Reads |
|---|---|---|---|
| `use_flexibleCN` | `false` | **`.true.`** | 0 |
| `MM_Nuptake_opt` | `false` | **`.true.`** | 0 |
| `CNratio_floating` | `false` | **`.true.`** | 3 (live) |
| `lnc_opt` | `false` | **`.true.`** | 0 |
| `vcmax_opt` | `0` | **`3`** | 0 |
| `CN_evergreen_phenology_opt` | `0` | **`1`** | 1 (live) |
| `carbon_resp_opt` | `0` | **`1`** (`0` if FUN on) | live branch at `cn_veg_carbon_flux.jl:1376` |
| `use_nguardrail` | `false` | **`.true.`** | 0 |

These are **not hypothetical**: CLM.jl's BGC work runs with `use_cn=true`, so
the port is running CLM5 BGC with the flexible-CN package off. This is flagged
here rather than fixed — flipping flexible-CN on is a model-configuration change
of a different magnitude to a single wrong boolean, and it needs its own
validation campaign against a Fortran reference (see "Not fixed here" below).

---

---

## Follow-up campaign #2 — closing the INERT half (this PR)

#256 deliberately changed nothing and pinned all 27 mismatches at their current
(wrong) values. Those pins were never meant to calcify into "expected". This
section records the first tranche closed: **the 12 INERT `varctl` mismatches**.

### Re-verification of the INERT classification

Every one was re-checked independently of the #256 sweep, by grepping *each*
read of the field and confirming the value cannot reach a branch under the
default configuration. `varctl reads` counts `varctl.<field>` in `src/`
(excluding the declaration and the `control.jl` echo-print).

| Field | −> CTSM | `varctl.` reads | Why the flip cannot move anything |
|---|---|---|---|
| `nsegspc` | 20 → **35** | 0 | Only consumer is `decomp_init.jl`, which takes its own `nsegspc::Integer = 20` kwarg; nothing wires `varctl` to it. |
| `nyr_forcing` | 10 → **1** | 0 | No read anywhere in `src/`. |
| `h2osno_max` | −999.0 → **10000.0** | 0 | Snow capping uses `varcon.H2OSNO_MAX` (already 10000.0). The `varctl` copy is a dead duplicate. |
| `n_dom_landunits` | −1 → **0** | 2 | Both guards are `> 0` (`surfdata.jl:324`, `control.jl:126`); −1 and 0 are indistinguishable. |
| `n_dom_pfts` | −1 → **0** | 0 | Guard `> 0` at `control.jl:126` only. |
| `toosmall_soil/_crop/_glacier/_lake/_wetland/_urban` | −1.0 → **0.0** (×6) | 2 each | `surfdata.jl:311` guards with `any(x -> x > 0.0, …)`; `surfrd_utils.jl:135` guards with `toosmall_any > 0.0`. Both reject −1.0 and 0.0 alike. |
| `downscale_hillslope_meteorology` | false → **true** | 0 | No read. The physics consumers (`topo.jl`, `surface_albedo.jl`) take their own kwarg and are additionally gated on `col.is_hillslope_column`, and `use_hillslope` defaults false. |

All twelve CTSM values were re-read from
`namelist_defaults_ctsm.xml` under `phys="clm5_0"`, `structure="standard"`:
lines 2404, 688, 465, 2387, 2390, 2393–2398, 665 respectively.

**Byte-identity is the proof of inertness.** The full suite was run to
completion on unmodified `main` and again after the twelve flips, Julia 1.12,
`--check-bounds=yes`. See "Verification" below. A flip that moved the suite
would have been evidence the flag was *not* inert, and would have been reverted
rather than accommodated.

### Also corrected: `surfrd_utils.jl` kwarg defaults

`collapse_individual_lunits!` declared `toosmall_* = -1.0` while
`dyn_subgrid_control.jl` already declared the same six at `0.0`. Aligned on
`0.0` (CTSM). Inert for the same reason as above — the `toosmall_any > 0.0`
early-return rejects both.

### NOT closed, and why — the internal defaults that are *not* provably inert

Three internal function-level defaults disagree with CTSM but are **reachable**,
so they are recorded rather than flipped:

- `decomp_init.jl:39,74` — `nsegspc::Integer = 20` (CTSM 35). Live whenever
  `numg/nclumps` falls between 20 and 35: it flips `seglen1` and therefore
  selects segment-based instead of round-robin gridcell decomposition. The real
  defect is that **nothing wires `varctl.nsegspc` into `decompInit!` at all**;
  fixing the wiring is the correct change, not moving a second literal.
- `topo.jl:49,112,189` and `surface_albedo.jl:1297` —
  `downscale_hillslope_meteorology::Bool = false` (CTSM `.true.`). Gated on
  `col.is_hillslope_column`, so unreachable while `use_hillslope=false`, but a
  caller that enables hillslopes reaches it. Belongs with the hillslope campaign.

## M2 revisited — `baseflow_scalar`: the contradiction is real, and worse than a 10×

#256 recorded that CLM.jl holds two values for one constant: `1.0e-2` on the
driver path and `0.001` (the CTSM value) in `src/calibration/parameters.jl:72`.
Chasing which one is live produced a **separate, larger finding**.

**Which is live: only `1.0e-2`.** The chain is
`clm_run.jl:96 (baseflow_scalar = 1.0e-2)` →
`init_soil_hydrology_config` → `BASEFLOW_SCALAR[]` (`soil_hydrology.jl:44,87`) →
`qflx_latflow_out` (`soil_hydrology.jl:2750, 3002`). `soil_hydrology.jl:66,77`
agree at `1.0e-2`. That is the CTSM value for `lower_boundary_condition=1`, or
for `clm4_5` — **not** for `clm5_0`, where `lbc=2` gives `0.001`
(`namelist_defaults_ctsm.xml:195-197`, re-read and confirmed).

**The calibration `0.001` is unreachable — its applier crashes.**
`_apply_baseflow_scalar!` (`parameters.jl:16-22`) does:

```julia
inst.overrides.baseflow_scalar = val          # correct: clm_run.jl:145 consumes this
for c in eachindex(inst.column.baseflow_scalar)   # ColumnData HAS NO SUCH FIELD
    inst.column.baseflow_scalar[c] = val
end
```

Verified live on Julia 1.12:

```
ColumnData.baseflow_scalar? false
ColumnData.fff?             false
param baseflow_scalar default=0.001 -> APPLIER THREW: FieldError: type CLM.ColumnData has no field `baseflow_scalar`
param fff             default=0.5   -> APPLIER THREW: FieldError: type CLM.ColumnData has no field `fff`
```

So **both** hydrology calibration parameters are dead: calibrating
`baseflow_scalar` or `fff` throws on the first `apply!`. The `inst.overrides.…`
assignment on the line above is the real, working injection point (consumed at
`clm_run.jl:145-146` and `:150-151`); the array loop is unreachable-by-design
code that could never have run.

**Why nobody noticed** — `test/test_calibration.jl:59-62` asserts only
`hydro[1].name == "baseflow_scalar"`. It never calls `apply!`. A textbook
instance of the vacuous-check class: the test proves the parameter is *listed*,
not that it *works*.

**What is fixed here** — the bogus array loops are deleted (they can only
throw), and `test_calibration.jl` now actually invokes `apply!` and asserts the
override lands. This is byte-identical for every non-calibration caller.

**What is deliberately NOT fixed here** — the value. Repairing the applier means
`0.001` can now reach physics through the calibration path while the driver path
still starts at `1.0e-2`, i.e. the 10× disagreement becomes genuinely live
rather than merely written down. Resolving it is a **physics change to the
drainage/baseflow rate** and needs the dedicated `baseflow_scalar` campaign with
Fortran parity evidence — exactly the treatment `use_hydrstress` and `use_luna`
are getting. Flipping it inside a sweep labelled "inert" would be the very bug
this audit exists to catch.

## M2 CLOSED — `baseflow_scalar` = `0.001`, and where the 10× actually bites

The dedicated campaign M2-revisited asked for. Default changed `1.0e-2` → `0.001`
at all four literals that carry it (`soil_hydrology.jl:44` `BASEFLOW_SCALAR[]`,
`:66` the `SoilHydrologyConfig` field, `:77` the `init_soil_hydrology_config`
kwarg, `clm_run.jl` the `clm_run!` kwarg). All four must move together: only
`:44`/`:87` is read by physics, but `test_soil_hydrology_mod.jl:59` calls bare
`init_soil_hydrology_config()`, so a partial change makes the live global
depend on test *ordering*.

### CTSM's conditional, verified three ways

```xml
<baseflow_scalar                lower_boundary_condition="2">0.001d00</baseflow_scalar>
<baseflow_scalar                lower_boundary_condition="1">1.d-2</baseflow_scalar>
<baseflow_scalar phys="clm4_5"  lower_boundary_condition="2">1.d-2</baseflow_scalar>
```

1. **The XML chain.** For `clm5_0`: `soilwater_movement_method=1`
   (`defaults:419`) + `use_bedrock=.true.` (`defaults:180`) ⇒
   `lower_boundary_condition=2` (`defaults:426`) ⇒ row 195 ⇒ **`0.001`**.
   Row 196 needs `lbc=1`; row 197 needs `phys=clm4_5`. Neither matches, so the
   resolution is unambiguous. (Under `clm4_5` proper, `swmm=0` ⇒ `lbc=4`, and
   `baseflow_scalar` is not used at all — the definition entry says "ONLY used
   if lower_boundary_condition is not aquifer or table".)
2. **CTSM-emitted namelists on this machine.** Two independent
   build-namelist outputs, both `phys=clm5_0`:
   `domain_Bow_at_Banff_lumped/settings/CLM/lnd_in` and
   `domain_Peatland_MerBleue_Canada/settings/CLM/lnd_in` each contain
   `soilwater_movement_method = 1`, `use_bedrock = .true.`,
   `lower_boundary_condition = 2`, and **`baseflow_scalar = 0.001d00`**.
   This is the derivation confirmed end-to-end, not merely read off the XML.
3. **CLM.jl already agreed with itself everywhere except the driver.**
   `scripts/parity_run_domain.jl` sets `baseflow = 0.001` for **18 of its 20**
   domains; the two exceptions are DDS-*calibrated* (Bow `0.0022119554`,
   Stillwater `0.016035343`). `src/calibration/parameters.jl` already used
   `0.001`. The driver default was the lone `1.0e-2` holdout — the same
   "code fallback copied instead of the namelist default" pattern as
   #252/#259/#273. `SoilHydrologyMod.F90:71`'s `= 1.e-2_r8` is exactly that
   pre-namelist initialiser.

### Measured effect of the 10× (this is a real physics change)

`qflx_latflow_out` is guarded by `if zwt <= zi(nbedrock)` — the water table must
reach *into* the bedrock layer for the power-law baseflow term to exist at all.
That guard, not the coefficient, dominates:

**Bow at Banff (240 d from the Fortran spun-up 2003 restart): exactly zero effect.**

| | `1.0e-2` | `0.001` |
|---|---|---|
| runoff (mm) | 377.032388 | 377.032388 |
| `qflx_latflow_out` (mm) | 0.000000 | 0.000000 |
| soil water (mm) | 284.044305 | 284.044305 |
| mean `zwt` (m) | 2.280000 | 2.280000 |

Bit-identical, because `nbedrock=12` ⇒ `zi_bedrock = 1.88 m` while `zwt` sits at
`2.28 m` — *below* bedrock — so `zwt <= zi_bedrock` is false at every one of the
11 520 steps and the branch never executes. **This is not a port artifact:** the
Fortran reference for the same site
(`domain_Bow_at_Banff_lumped/simulations/clm_dds_calibration/CLM/…h0.2004-12-31…`)
reports `QDRAI = 0.0` and `ZWT = 2.28 m` for all 365 days, with
`QRUNOFF ≡ QOVER`. CLM.jl reproduces CTSM's own zero. At Bow the coefficient
multiplies a branch the reference model never enters either.

**MerBleue peat bog (240 d from its 2017 restart): the branch is live, and the
10× moves it.**

| | `1.0e-2` | `0.001` | Δ |
|---|---|---|---|
| `qflx_drain` / `latflow_out` (mm) | 63.690862 | 58.843222 | **−7.6 %** |
| `qflx_surf` (mm) | 637.480694 | 641.131278 | **+0.57 %** |
| total runoff (mm) | 701.171557 | 699.974499 | −0.17 % |
| soil water (mm) | 745.978369 | 747.278950 | +0.17 % |
| mean `zwt` (m) | 2.273917 | 2.223477 | 5.0 cm higher |
| end `zwt` (m) | 2.280000 | 2.064378 | 21.6 cm higher |
| ET (mm) | 260.082967 | 259.979444 | −0.04 % |

Two things worth stating plainly:

* **A 10× cut in the coefficient produces only a 7.6 % cut in baseflow.** The
  term is strongly self-limiting: less drainage ⇒ water accumulates ⇒ `zwt`
  rises ⇒ `(zi_bedrock − zwt)^n_baseflow` grows ⇒ drainage partly recovers. The
  water table, not the coefficient, sets the flux.
* **Total runoff barely moves (−0.17 %) while the partition moves much more**
  (baseflow −7.6 %, surface runoff +0.57 %). A water-balance check sees almost
  nothing here — the [CONSERVATION IS NOT ACCURACY] pattern again. The number
  that changed is the one no balance check constrains.

MerBleue's own CTSM-emitted `lnd_in` specifies `0.001`, so the direction is
toward the reference, not away from it.

### Why the test suite is nearly blind to this

`test/test_hydrology_drainage.jl:281` calls `hydrology_drainage!` with
`mask_hydrology = BitVector([false, false])`, so the drainage body never runs.
`test_soil_lateral_flow.jl:274` and `test_soil_hydrology_mod.jl:465` call
`subsurface_lateral_flow!` directly but assert only finiteness and
`qflx_latflow_out ≥ 0` — one-sided bounds that hold for *any* non-negative
coefficient. Every Fortran-parity test passes `baseflow_scalar` explicitly
(`scripts/fortran_parity_common.jl:41` = the calibrated `0.0022119554`). The
only assertion on the value itself was `test_soil_hydrology_mod.jl:45`, updated
here, plus the new pin in `test_driver_defaults_audit.jl`.

## M6 revisited — classifying the 8 `use_cn` cascade flags

Asked whether any of the cascade is inert under the default `use_cn=false` and
therefore closeable now. Result: **none should be closed now, but for two
different reasons.**

| Flag | `varctl.` reads | Status |
|---|---|---|
| `use_flexibleCN` | 0 | **Unwired.** `fun.jl` takes its own `use_flexiblecn` argument; `soil_water_movement.jl:1746` takes its own kwarg defaulting `false`. Neither reads `varctl`. |
| `MM_Nuptake_opt` | 0 | **Unwired** — zero occurrences in `src/` outside the declaration. |
| `lnc_opt` | 0 | **Unwired** — zero occurrences. |
| `vcmax_opt` | 0 | **Unwired** — appears only in two `history_writer.jl` comments. |
| `CN_evergreen_phenology_opt` | 1 | **Unwired, and a known divergence.** The single "read" is inside a comment block (`phenology.jl:1028-1042`); the code sets `tranr = 0.0002` *unconditionally*, where Fortran applies it only when the flag is 1. Documented in-source as deliberate: the 16-biome scorecard was validated against the always-on behaviour. |
| `CNratio_floating` | 3 | **Live under `use_cn=true`** — `phenology.jl:952,958,962`, inside `cn_phenology!`. Inert under the default `use_cn=false`. |
| `carbon_resp_opt` | 0 | **Live under `use_cn=true`**, but through `CNDriverConfig.carbon_resp_opt` (`cn_driver.jl:44`, branch at `cn_veg_carbon_flux.jl:1376`), not through `varctl`. |
| `use_nguardrail` | 0 | **Live under `use_cn=true`**, likewise through `CNDriverConfig` (`cn_driver.jl:48` → `cn_precision_control.jl:82,92`). |

**Verdict.** Four flags (`use_flexibleCN`, `MM_Nuptake_opt`, `lnc_opt`,
`vcmax_opt`) plus `CN_evergreen_phenology_opt` have zero live reads and are
technically inert *even under `use_cn=true`* — flipping them would pass a
byte-identity suite. They are **still not closed here**, deliberately: setting
`varctl.use_flexibleCN = true` while every consumer independently defaults
`false` would make the control state *assert* that the port runs flexible-CN
when it demonstrably does not. That is a worse landmine than the current honest
`false`, not a smaller one. The correct fix for all five is **wiring plus
validation**, which is the CN campaign.

The remaining three (`CNratio_floating`, `carbon_resp_opt`, `use_nguardrail`)
are live the moment `use_cn=true`, which is how all BGC work runs. They need the
CN campaign outright.

## Scorecard

| Category | Count |
|---|---|
| Flags audited | **~110** (16 driver-entry keywords + 79 `varctl` fields + `CLMDriverConfig` switches + the `use_flexibleCN` cascade) |
| MATCH | ~85 |
| MISMATCH — live physics | **1 remaining** (`convert_ocean_to_land`) + `glc_snow_persistence_max_days` (live via kwarg). Originally 6: `create_crop_landunit` **fixed** (M5), `use_luna` **fixed** (#267, M4), `use_hydrstress` **fixed** (PHS campaign, M3), `use_aquifer_layer` **fixed** (#259, re-confirmed here on the corrected tridiagonal solver), `baseflow_scalar` **fixed** (here, M2) |
| MISMATCH — inert (0 reads / guarded) | 12 (`nsegspc`, `nyr_forcing`, `h2osno_max`, `n_dom_landunits`, `n_dom_pfts`, 6× `toosmall_*`, `downscale_hillslope_meteorology`, `z0param_method`) |
| CONDITIONAL-MISMATCH under `use_cn=true` | 8 (the `use_flexibleCN` cascade, M6) |
| Namelist-vs-code-fallback disagreements found | **14** — CLM.jl copied the code fallback in *every* case |

## Not fixed here, and why

Three of the live mismatches are deliberately **documented rather than
changed**, because each is a model-configuration decision rather than a
mechanical correction, and flipping them silently would be the mirror image of
the bug this audit exists to catch:

- ~~**`use_hydrstress` / `use_luna` (M3, M4).**~~ **Both are now CLOSED**, each as
  its own change with its own evidence, exactly as this section asked: `use_luna`
  in #267 and `use_hydrstress` in the PHS campaign
  (`docs/PHS_DEFAULT_BLOCKERS.md`). Both landed as CONDITIONAL defaults
  (`nothing` → `!use_fates`), not bare flips.
- **The `use_flexibleCN` cascade (M6).** Eight coupled flags; enabling
  flexible-CN needs a validation campaign against a Fortran CN reference, not a
  one-line default change.

These are recorded as **open, evidenced findings**. The value of writing them
down without fixing them is that the next person to see
`use_hydrstress = true` in a harness will know it is a workaround for a known
wrong default, not a deliberate experimental setting.

## Verification

- `test/test_driver_defaults_audit.jl` standalone, `--check-bounds=yes`: **70/70 pass**.
- Full suite, Julia 1.12, `--check-bounds=yes`:
  **26141 passed, 0 failed, 0 errored, 3 broken** (40m54s, exit 0).
  Subtracting this file's 70 new assertions gives a 26071 baseline, i.e. the
  suite is unchanged apart from the additions — as expected, since **this audit
  changes no default**. The three broken tests are the pre-existing ones.

A first full-suite attempt was killed by SIGTERM (`EXIT=143`) mid-run on the
loaded shared box; it was rerun to completion rather than reported as green.

### M2 close + M1 re-confirm (baseflow campaign)

Julia 1.12, full suite, baseline established on the same merge-base (`37db702`):

| Run | Pass | Fail | Broken | Total | Time |
|---|---|---|---|---|---|
| Baseline (unmodified `main`) | 26578 | 0 | 3 | 26581 | 27m15.6s |
| After (`baseflow_scalar = 0.001` + pins) | 26585 | 0 | 3 | 26588 | 27m01.1s |
| `use_aquifer_layer=true` (measurement only, not merged) | 26576 | **2** | 3 | 26581 | 27m03.7s |

Baseline → After is `+7` total / `+7` pass, exactly the seven `@test` lines of
the new `baseflow_scalar` testset; **no pre-existing assertion changed status.**

The `use_aquifer_layer=true` arm is the sharper result. Flipping the whole
soil-water lower boundary condition — a change the Bow probe shows moves
drainage by seven orders of magnitude — produced exactly **two** failures out
of 26 581, and both were `test_driver_defaults_audit.jl:100-101`, the assertions
that read the default out of the AST. Not one physics assertion noticed. The
only thing in the suite that saw it was the file written specifically to see
defaults, which is this file's own thesis demonstrated rather than argued.

## Method note

For each flag the check is three-layered, because agreeing with only one layer
is how #252 happened:

1. the **namelist default** with all its conditional variants,
2. the **code fallback** in the `*_inparm` initialiser,
3. whether CLM.jl's value is **actually read** anywhere (an unread field cannot
   cause a physics bug today but will silently mislead the first consumer).

A flag is only marked MATCH when the resolved namelist default under this port's
conditions agrees, *and* any disagreement with the code fallback is understood.
