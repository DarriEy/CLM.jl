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
(campaign #3 used the ctsm5.3.012 checkout at
`/Users/darri/repos/SYMFLUENCE_data/installs/clm/ctsm5.3.012/`; line numbers below
from that campaign refer to its `namelist_defaults_ctsm.xml`)

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
| `use_lch4` | `false` | `.true.` unless `soil_decomp_method="None"` | **COND-MISMATCH — deliberate** (campaign #3: both BGC references override it to `.false.`; TWS_inversion needs the ch4finundated stream) | methane biogeochem off by default in CN runs | — |
| `use_c13` / `use_c14` | `false` | `.false.` | **MATCH** | — | — |
| `use_fun` | `nothing` → `use_cn` (campaign #3; was `false`) | `.true.` for `use_cn .and. use_nitrif_denitrif` (`:2020-2022`) | **CLOSED** | FUN N uptake/fixation off in every default CN run | `fortran_parity_common.jl` (now redundant) |
| `use_flexibleCN` | `nothing` → `use_cn` (campaign #3; was not a kwarg, `false`) | `.true.` for `use_cn` (`:573-575`) | **CLOSED** (cascade: see campaign #3) | flexible leaf C:N off | `fortran_parity_common.jl` (now redundant) |
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
| `convert_ocean_to_land` | ~~`false`~~ → **`true`** | **`.true.`** (`defaults:2382`, unconditional). Code fallback `.false.` | **CLOSED** (campaign #3; inert on land-only surfdata) | 1 |
| `co2_ppmv` | ~~`355.0`~~ → **`367.0`** | I2000 compsets pass `CCSM_CO2_PPMV = 367.0` (every reference build line has `-co2_ppmv 367.0`); `defaults:25` 379.0 only applies when none is given; code fallback 355.0 | **CLOSED** (campaign #3; unread — physics is `forcing_reader.jl`'s `367e-6*pbot`) | 0 |
| `glc_snow_persistence_max_days` | ~~`7300`~~ → **`0`** | **`0`**; `7300` is the **clm4_5** variant (`defaults:540-541`) | **CLOSED** (campaign #3; the `clm_driver.jl` SMB call now passes `varctl`'s value) | 1 |
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

### M6 — the `use_flexibleCN` cascade (CONDITIONAL-MISMATCH cluster, applies to every CN run) — **CLOSED in campaign #3**

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

## Campaign #3 — the defaults of a standard CTSM 5.3.012 case

Goal (approved by the project owner, knowing it changes default results):
`clm_initialize!` / `clm_run!` called with default arguments run the configuration a
standard CTSM 5.3.012 case runs. Every value below was read from the ctsm5.3.012 tree
(`installs/clm/ctsm5.3.012/bld/namelist_files/namelist_defaults_ctsm.xml`; `:NNN` is a
line number in that file) and cross-checked against what build-namelist actually
wrote for the reference cases on this machine: the SP Bow case
(`domain_Bow_at_Banff_lumped/settings/CLM/lnd_in`, the strict-gate configuration) and
the two BGC references (`runs/cn_annual_ref/lnd_in` and the cold-start 2-year reference
of #362, both `-bgc bgc -use_case 2000_control`, `phys=clm5_0`).

### The pattern, again

Every changed row is the #252 pattern: CTSM has a namelist default and a module code
fallback that disagree, and the port had copied the code fallback. Two were worse than
a wrong value, because the code said they *were* CTSM's default:
`cn_veg_carbon_state.jl` ("its CTSM default is 20 gC/m2" — 20 is the fallback; the
clm5_0 BGC default is 100) and `snow_hydrology.jl` (Anderson1976 + no wind "for AD" —
the clm4_5 values; the AD tests pass with the clm5_0 values, see below).

**How the values are applied.** Following the #265 principle, per-instance structs
(`PhotosynthesisData`, `CNSharedParamsData`, the CN facade and driver configs,
`Atm2LndParamsData`, `Lnd2AtmParamsData`, `VarCtl`) keep CTSM's code fallbacks, and a
new step in `clm_initialize!` (`_clm5_namelist_defaults!`, just before the cold BGC init) writes the
namelist defaults onto the instance, the way build-namelist writes them into `lnd_in`.
Callers still override any of them on `inst` after `clm_initialize!` returns.
Module-global switches that callers set *before* `clm_initialize!` (the snow namelist
`Ref`s, `rooting_profile_config`, `cnvegcstate_const`) take the namelist default as
their initial value instead, so `clm_initialize!` never clobbers a caller's setting.

`clm_run!` used to apply three of these itself (`light_inhibit`, the atm2lnd rain/snow
repartitioning, the lapse rate) *after* `clm_initialize!`. Callers that drive
`clm_initialize!` + `clm_drv!` directly (every Fortran-parity harness, the AD tests,
the FATES harnesses) did not get them. They now live only in `clm_initialize!`.

### Defaults changed

| # | Setting | Old port default | New default | CTSM 5.3.012 source | Scope |
|---|---|---|---|---|---|
| 1 | `use_fun` | `false` | `true` when `use_cn` (`nothing` kwarg on `clm_initialize!`) | `:2020-2022` (`use_cn .and. use_nitrif_denitrif`); Bow BGC lnd_in `use_fun = .true.` | CN |
| 2 | `use_flexibleCN` | `false` | `true` when `use_cn` (`nothing` kwarg) | `:573-575`; lnd_in `use_flexiblecn = .true.` | CN |
| 3 | `lnc_opt` | `false` | `true` when `use_cn` | `:596-597`; lnd_in `lnc_opt = .true.` | CN |
| 4 | `vcmax_opt` | `0` | `3` when `use_flexibleCN` | `:586`; lnd_in `vcmax_opt = 3` | CN |
| 5 | `MM_Nuptake_opt` | `false` | `true` when `use_flexibleCN` | `:583` | CN (varctl only: the port's cold init already hard-codes the `.true.` branch of `CNVeg{Carbon,Nitrogen}StateType` InitCold) |
| 6 | `CNratio_floating` | `false` | `true` when `use_flexibleCN` | `:584` | CN (live: `phenology.jl` litterfall C:N) |
| 7 | `CN_evergreen_phenology_opt` | `0` | `1` when `use_flexibleCN` | `:587` | CN (live since the evergreen follow-up below: gates `cn_evergreen_phenology!`'s storage → transfer and transfer → display) |
| 8 | `carbon_resp_opt` | `0` | `0` with FUN, `1` for flexibleCN without FUN | `:588-589` | CN (unchanged in the default case) |
| 9 | `use_nguardrail` | `false` | `true` when `use_cn` | `:642-643` | CN |
| 10 | `leafresp_method` | `1` (Ryan 1991) | `2` (Atkin 2015) when `use_cn` | `:102-104`; lnd_in `leafresp_method = 2` | CN (SP never reads it) |
| 11 | `br_root` | `SPVAL` → `br` = 2.525e-6 | `0.83e-6` when `use_cn` | `:2027`; lnd_in `&cnmresp_inparm br_root = 0.83d-06` | CN |
| 12 | `constrain_stress_deciduous_onset` | `false` | `true` when `use_cn` | `:635-636` | CN |
| 13 | `initial_vegC` | `20` | `100` | `:2033-2035` (100 for `use_cn .and. mm_nuptake_opt`; 20 = code fallback `CNVegCarbonStateType.F90:232`) | CN cold start |
| 14 | CN `ncrit` / `nnegcrit` | `1e-8` / `-7` | `1e-9` / `-6` | `:646`, `:652` (code fallback `CNPrecisionControlMod.F90:27-28`); lnd_in `&cnprecision_inparm` | CN |
| 15 | `modifyphoto_and_lmr_forcrop` | `false` | `true` | `:107-108`; both lnd_in | all (PHS canopy lmr bsun/bsha weighting) |
| 16 | `light_inhibit` | `true` only via `clm_run!` | `true` from `clm_initialize!` | `:98-99` | all (unchanged for `clm_run!` callers) |
| 17 | atm2lnd `repartition_rain_snow`, `lapse_rate`, `lapse_rate_longwave`, `precip_repartition_*` | set only by `clm_run!` | set by `clm_initialize!` | `:209-234` | all (unchanged for `clm_run!` callers) |
| 18 | `melt_non_icesheet_ice_runoff` | `false` | `true` | `:237-238`; both lnd_in | all (ice runoff from non-ice-sheet columns) |
| 19 | `rooting_profile_method_water` / `_carbon` | Zeng 2001 (0) | Jackson 1996 (1) | `:184`, `:188` (code fallback `RootBiophysMod.F90:64-65`) | all |
| 20 | `snow_overburden_compaction_method` | Anderson1976 | Vionnet2012 | `:486-487` | all |
| 21 | `wind_dependent_snow_density` | `false` | `true` | `:483-484` | all |
| 22 | `convert_ocean_to_land` | `false` | `true` | `:2382` | all (inert on land-only surfdata) |
| 23 | `glc_snow_persistence_max_days` | `7300` | `0` | `:540-541` (7300 is clm4_5) | glacier SMB diagnostic |
| 24 | `varctl.co2_ppmv` | `355.0` | `367.0` | I2000 `CCSM_CO2_PPMV` (`-co2_ppmv 367.0` on every reference build line); unread — physics uses `forcing_reader.jl`'s `367e-6*pbot` | none (bookkeeping) |
| 25 | History `c2l_scale_type` / `l2g_scale_type` | `unity` on every field | the Fortran `hist_addfld` value, 47 fields (36 `urbanf`, 1 `urbans`, 10 `veg`) | `EnergyFluxType.F90`, `SolarAbsorbedType.F90`, `WaterFluxType.F90`, `WaterDiagnostic*Type.F90`, `TemperatureType.F90`, `SoilHydrologyType.F90`, `SaturatedExcessRunoffMod.F90` | output only, multi-landunit gridcells |
| 26 | History `BTRAN` | daily-min snapshot | hourly-averaged `btran_patch` (`l2g='veg'`); new `BTRANMN` = the daily-min field | `EnergyFluxType.F90:661-669` | output only |

Also: `clm_initialize!` now warns once when `fsnowoptics`/`fsnowaging` are empty (the
SNICAR tables stay zero; CTSM's defaults are the 5-band inputdata files, `:2039,2044`).

### Deliberately left different

| Setting | CTSM default | Port default | Why |
|---|---|---|---|
| `use_lch4` | `.true.` (`:2538`, any BGC) | `false` | Both BGC references override it to `.false.`, so nothing validates the CH4 default on the multi-year path; and the clm5_0 `finundation_method='TWS_inversion'` reads its inversion parameters from the `ch4finundated` stream file, which the port does not read. |
| `fire_method` | `li2016crufrc` (`:242`) | `:nofire` | Needs the lightning and population-density streams (input data, not configuration). Both references run `nofire`. |
| `glcmec_downscale_longwave` / `longwave_downscaling_limit` | `.true.` / `0.5` (`:212`, `:220`) | `false` / unset | CTSM downscales longwave only over TopoMod's downscale filter (glacier-MEC and hillslope columns, `needs_downscaling_col`), which the port never populates in the driver (`topo_update!` is not called; `clm_run!` sets `topo_col = forc_topo`), and `downscale_longwave!` applies to every column without CTSM's `build_normalization`. Turning it on would be less faithful, not more. |
| `reseed_dead_plants` | `.false.` (`:60`; `.true.` only with `clm_accelerated_spinup=on`) | `false` | Matches the standard case. The BGC references were built with `-clm_accelerated_spinup on` and so carry `.true.`; that is a reference-case setting. |
| `fsnowoptics` / `fsnowaging` | inputdata files | `""` (with a warning) | Input paths, like `fsurdat`; the port has no inputdata root. |
| `glacier_region_ice_runoff_behavior` | `'melted'` for region 0 | ice-sheet columns never melted | Not wired (`ice_runoff_melted_grc` is never set); matters only for ISTICE columns. |
| hist `set_urb` / `set_lake` / `set_nourb` | per field | not ported | Only the scale types were in scope; these spval masks are a separate history gap. |
| `dtime` | 1800 | 1800 | Matches; the reference cases run 3600 by case setting. |

### The FUN / flexibleCN set on the CN path

The Bow Fortran-parity harnesses (`fortran_parity_common.jl` `build_bow_inst`) already
set items 1-4, 10, 11, 15 and 16 by hand, which is why the single-step CN parity tests
did not move. The two multi-year CN harnesses did **not**: `fortran_parity_cn_annual.jl`
and the #362 cold-start harness inherited `use_fun=false`, flexibleCN off,
`lnc_opt=false`, Ryan leaf respiration, `br_root=SPVAL`, the onset constraint off, the
Anderson snow pair and (cold start) `initial_vegC=20`, against Fortran references that
ran with the CTSM values. Both harnesses were run unmodified before and after, with
`CESM_INPUTDATA` pointing at the local inputdata. (Without it they silently run with
zeroed SNICAR tables — `fortran_parity_common.jl` falls back to `""` when the file is
missing. That is how #362 ran them; the new `clm_initialize!` warning is what exposed
it. The numbers below are with the real SNICAR files; the zeroed-SNICAR rerun moves
in the same direction for every field.)

**cn_annual** (§8 reference: Fortran 2202 BGC restart injected, 365 days, Bow).
Annual mean of the monthly means; rel = |J−F| / (1+|F|).

| field | Fortran | before | after | rel before | rel after |
|---|---|---|---|---|---|
| TOTVEGC (gC/m2) | 374.37 | 366.01 | 374.22 | 2.2e-2 | **4.0e-4** |
| TOTECOSYSC | 41672 | 41666 | 41672 | 1.6e-4 | **<1e-5** |
| SMINN (gN/m2) | 51.175 | 51.479 | 51.172 | 5.8e-3 | **6.4e-5** |
| GPP (gC/m2/s) | 2.454e-6 | 2.278e-6 | 2.415e-6 | -7.2 % | **-1.6 %** |
| NPP (gC/m2/s) | 1.690e-6 | 7.77e-7 | 1.676e-6 | **-54 %** | **-0.8 %** |
| EFLX_LH_TOT (W/m2) | 8.317 | 8.248 | 8.265 | 7.4e-3 | 5.6e-3 |
| TLAI | 0.06611 | 0.06387 | 0.06263 | 2.1e-3 | 3.3e-3 (away) |
| LEAFC (gC/m2) | 9.891 | 9.559 | 9.374 | 3.1e-2 | 4.8e-2 (away) |
| TG (K) | 273.62 | 273.63 | 273.64 | 3e-5 | 6e-5 (away, 0.01 K) |
| FSH (W/m2) | 39.644 | 39.696 | 39.587 | 1.3e-3 | 1.4e-3 (away, 0.005 W/m2) |

NPP was the big error: with FUN off and Ryan leaf respiration the port's Bow
vegetation had **negative** NPP from August to December where Fortran's is positive.
It is now within 1 % annually and within 0.3e-6 in every month.

**The two that moved away were investigated before shipping.** Each default was
reverted alone on top of the new configuration (`CN_TOGGLE` runs of the same
harness, zeroed-SNICAR configuration):

| reverted alone | TLAI rel | LEAFC rel | TOTVEGC rel | NPP rel |
|---|---|---|---|---|
| (none, new defaults) | 3.4e-3 | 4.9e-2 | 1.4e-3 | 4.7e-8 |
| `constrain_stress_deciduous_onset` | 3.4e-3 | 4.9e-2 | 1.4e-3 | 4.7e-8 |
| `CNratio_floating` | 3.4e-3 | 4.9e-2 | 1.5e-3 | 5.1e-8 |
| `leafresp_method` (Ryan) | 3.4e-3 | 5.0e-2 | 0 | 9e-10 |
| `lnc_opt` + `vcmax_opt` | 3.4e-3 | 4.9e-2 | 6.7e-3 | 1.6e-7 |
| `use_fun` | 3.4e-3 | 4.9e-2 | 1.0e-3 | 3.8e-8 |
| snow pair (Anderson, no wind) | 3.5e-3 | 5.0e-2 | 1.6e-3 | 5.0e-8 |
| **`use_flexibleCN`** | **3.7e-4** | **5.4e-3** | 1.9e-3 | 5.6e-8 |

The whole LEAFC/TLAI move comes from `use_flexibleCN`, and the leaf deficit opens in
June and grows through December (LEAFC Dec: Fortran 4.45, before 3.68, after 3.37), so
it is in the growing-season allocation, not in onset. Fortran runs `use_flexibleCN =
.true.` (both BGC `lnd_in`), so the CTSM value is the right default, and the move
exposes a residual in the port's flexible-CN leaf allocation that the FUN-off,
flexible-off configuration had been offsetting. It is recorded as an open finding, not
fixed here. (Ryan leaf respiration happens to give a better TOTVEGC/NPP here, but
Fortran's lnd_in has `leafresp_method = 2`, so that would be compensation, not parity.)
The TG/FSH moves are 0.01 K and 0.005 W/m2.

**Cold-start 2-year CN** (#362's harness, Bow cold start 2002-01-01, against the
`cn-alt-roll/fortran_cold2y` ctsm5.3.012 reference; annual means).

| field | year | Fortran | before | after | rel before | rel after |
|---|---|---|---|---|---|---|
| TOTVEGC (gC/m2) | 1 | 175.04 | 40.29 (**0.23×**) | 177.20 (**1.01×**) | 0.77 | 1.2e-2 |
| TOTVEGC | 2 | 136.25 | 44.69 (0.33×) | 139.84 (1.03×) | 0.67 | 2.6e-2 |
| TOTECOSYSC | 2 | 1129.6 | 968.4 | 1132.6 | 0.14 | 2.7e-3 |
| TLAI | 1 | 0.402 | 0.082 | 0.399 | 0.23 | 2.2e-3 |
| FROOTC | 1 | 58.2 | 13.5 | 57.0 | 0.76 | 2.1e-2 |
| SMINN (gN/m2) | 2 | 0.0180 | 0.498 | 0.0191 | 0.47 | 1.1e-3 |
| GPP | 2 | 2.28e-6 | 1.24e-6 | 2.17e-6 | −46 % | −4.9 % |
| NPP | 2 | 8.75e-7 | 6.34e-7 | 8.69e-7 | −28 % | −0.6 % |
| NFIX (gN/m2/s) | 2 | 1.84e-9 | 0 | 1.55e-9 | — | — |
| EFLX_LH_TOT | 2 | 8.57 | 7.55 | 8.51 | 0.11 | 6.3e-3 |
| TG (K) | 1 / 2 | 273.10 / 273.78 | 273.27 / 274.03 | 273.31 / 274.27 | — | away (+0.04 / +0.24 K) |

(The harness still passes its own `h2osfcflag = 0`, `baseflow_scalar = 0.0022`
and `int_snow_max = 3113` where the reference lnd_in has 1 / 0.001 / 2000; those are
harness choices, unchanged before and after.) Every carbon, nitrogen and water field moves toward Fortran; most by one to two
orders of magnitude. The "TOTVEGC about 0.3×" of #362 was mostly `initial_vegC`
(20 vs 100: the cold vegetation started with a fifth of Fortran's carbon) plus FUN
being off (zero N fixation). TG is the one field that moves away. Reverting the snow
pair alone on top of the new defaults changes TG by < 0.003 K (`CN_TOGGLE=snow`), so
it is not the snow namelist; it comes with the 4× larger, now Fortran-sized canopy
(TLAI 0.08 → 0.40), i.e. a residual in the canopy/ground energy partition that the
undersized canopy was hiding. Recorded as a follow-up, not reverted.

### Snow namelist and AD

The `snow_hydrology.jl` comment said the wind term's `(·)^8.8` "breaks the strict
FD-vs-AD derivative check at deep-winter (250 K) conditions (test_ad_robustness.jl)".
Re-measured with the clm5_0 values as the defaults:

| `test_ad_robustness.jl` "Winter cold (250 K)" | old pair (Anderson1976, no wind) | clm5_0 pair (Vionnet2012 + wind) |
|---|---|---|
| d(LH)/dT  AD vs FD | 0.0505 vs **NaN** | 0.0495 vs 0.0497 (0.28 %) |
| d(SH)/dT  AD vs FD | -14.97 vs **NaN** | -15.08 vs -14.16 (6.5 %, tolerance 10 %) |
| d(T_grnd)/dT AD vs FD | 0.382 vs **NaN** | 0.375 vs 0.376 (0.13 %) |
| assertions run | 90 | 96 |

So the claim was backwards. Under the old pair the deep-winter finite-difference
perturbation runs return NaN, and the strict FD-vs-AD comparison was skipped without
failing (`abs(NaN) > 1e-4` is false), a vacuous pass. Under the clm5_0 pair the FD is
finite and agrees with AD, and the six extra assertions run and pass. The same result
holds with every other new default in place and only the snow pair reverted
(`snow_hydrology_ad_smooth_namelist!()` before the test: FD = NaN, 90 assertions), so
the NaN belongs to the Anderson/no-wind path, not to another change. The NaN itself is
not chased here; it is a finding for the AD campaign.

With the new defaults: `test_ad_robustness.jl` 96/96, `test_ad_e2e.jl` 122/122,
`test_ad_smoke.jl` 28/28, `test_phase_snow_smoothing.jl` 48/48.

`docs/FORWARD_AD_SMOOTHNESS_BOUNDARY.md`'s snow analysis is about the discrete layer
events (`snl`, `imelt`, combine/divide), not about the compaction law or the new-snow
density. The old pair stays available as a documented opt-in,
`CLM.snow_hydrology_ad_smooth_namelist!()`, for an AD application that wants the snow
path that document characterised; nothing in the test suite needs it.

### Validation

**Multisite strict gate** (`parity_run_domain.jl`, 1 year each, scored with the
`parity_config.py` gate: annual |Δ| ≤ 1 % / 0.05 K and daily nRMSE ≤ 0.10 / 0.20 K).
"Before" is unmodified `origin/main` (42409bb); "after" is this branch, where
`parity_run_domain.jl` no longer sets the rooting profile, the snow namelist,
`use_hydrstress`, `use_luna` or `h2osfcflag` by hand (it inherits them).

| site | before | after | failing cells (after) |
|---|---|---|---|
| Bow | 68/69 | 68/69 | BTRANMN |
| Stillwater | 68/69 | 68/69 | SNOW_DEPTH |
| Abisko | 69/69 | 69/69 | — |
| Baltimore | 68/69 | 68/69 | SNOW_DEPTH (documented) |
| Savanna | 69/69 | 69/69 | — |
| Iceland (extra) | 69/69 | 69/69 | — |
| Massa (extra) | 64/69 | 64/69 | FGEV, QOVER, QSOIL, QFLX_INFL, SNOWLIQ |

Every one of the 102 history variables the two runs share is **bit-identical** at all
seven sites, except Baltimore and Iceland `QRUNOFF` (max |Δ| 1.4e-20 and 7.9e-20
mm/s, the order of a floating-point rounding) and `BTRAN`, which is now the hourly
field (`BTRANMN` after is bit-identical to `BTRAN` before, so the gate cell is
unchanged). That is the expected result: this harness already set the snow pair,
the rooting profile and PHS/LUNA by hand, `modifyphoto_and_lmr_forcrop` only changes
the canopy leaf-respiration diagnostic on the PHS SP path, and the CN changes do not
run with `use_cn = false`. What it proves is that the inherited defaults are exactly
the configuration the harness used to rebuild by hand.

**The scale types are inert on every reference domain.** All seven SYMFLUENCE
surface datasets are 100 % natural vegetation (`PCT_URBAN = PCT_GLACIER = PCT_LAKE
= 0`, including the "urban" Baltimore and the "glacier" Iceland/Massa domains), so
no reference gridcell has a second landunit for `urbanf`/`urbans`/`veg` to act on.
The aggregation itself is pinned by a unit test in `test_driver_defaults_audit.jl`
(one soil + one urban landunit with roof, sunwall and road columns; `unity`,
`urbanf`, `urbans` and `veg` against hand-computed `c2g` values).

**Fortran-parity and unit tests.** Full suite (Julia 1.12, `SYMFLUENCE_DATA` and
`CESM_INPUTDATA` set): 33217 passed, 31 failed, 4 broken (the `@test_broken`
markers). All 31 failures were assertions pinning the old defaults in three files:
`test_constants.jl` (`co2_ppmv`), `test_cn_veg_carbon_state.jl` (cold-start seeds
of `INITIAL_VEGC`), `test_cn_precision_control.jl` (truncation inputs between the
new and old `ncrit`). Updated; re-run standalone: 70/70, 198/198, 162/162, and
`test_driver_defaults_audit.jl` 157/157. No Fortran-parity test (`test_fortran_parity*`,
the FATES parity tests, `test_fun_coldstart_finite.jl`, `test_cn_active_layer_roll.jl`)
moved.

### Open findings from this campaign
1. ~~**flexibleCN leaf allocation residual**~~ — closed, see "Follow-up: evergreen
   transfer → display" below.
2. ~~**Cold-start TG +0.24 K in year 2**~~ — closed by the same fix.
3. **The old snow pair NaNs the 250 K FD perturbation** in `test_ad_robustness.jl`;
   that test silently skipped its strict check under the old defaults.
4. `fortran_parity_common.jl` falls back to empty SNICAR paths when
   `CESM_INPUTDATA` is unset, so CN harness runs without it use zeroed snow optics.

### Follow-up: evergreen transfer → display (closes findings 1 and 2)

**Root cause.** `CNEvergreenPhenology` with `CN_evergreen_phenology_opt == 1` does two
things for evergreen patches every step (`CNPhenologyMod.F90:736-831`): storage →
transfer at the fixed `tranr = 0.0002` (`:739-776`), then transfer → display with
`t1 = 1/dt`, i.e. the whole transfer pool (`:778-829`). The port
(`phenology.jl` `cn_evergreen_phenology!`) did the first unconditionally and never the
second, so `leafc_xfer` / `frootc_xfer` (and the wood and N transfer pools) of an
evergreen patch only ever filled: `LEAFC_XFER_TO_LEAFC` was 0 all year on the Bow
needleleaf tree where Fortran's is ~6e-8 gC/m2/s, and by December the tree's
`leafc_xfer` was ~8000× Fortran's.

This was hidden while `use_flexibleCN` was off. The CLM4.5 allocation uses the PFT's
`fcur` for evergreens, so their growth went straight to display. FlexibleCN sets
`fcur = 0` for evergreens (`NutrientCompetitionFlexibleCNMod.F90:412-414`, ported
faithfully at `nutrient_competition_flexiblecn.jl:206`), so all new evergreen leaf and
fine-root C and N go to storage and reach display only through the missing block. The
flexibleCN allocation, demand and FUN branches themselves were checked line by line
against Fortran and match. The difference was found with a second, patch-level
daily-mean history tape on the pristine ctsm5.3.012 exe (no rebuild: `hist_fincl2`
with ~90 CN fields, `hist_dov2xy = .false.`) diffed against the same fields dumped
per patch from the Julia harness. The first field that is wrong is
`LEAFC_XFER_TO_LEAFC` (−100 % from day 1). Everything downstream on the tree patch
(`LEAFC`, `TLAI`, `FROOTC`, `GPP`) follows it.

**Fix.** `cn_evergreen_phenology!` now runs both halves when
`CN_evergreen_phenology_opt == 1` (new `phen_evergreen_xfer!` kernel for the
transfer → display fluxes) and neither when it is 0, as Fortran does. The value comes
from `varctl`, which `clm_initialize!` sets to 1 whenever `use_flexibleCN`
(`namelist_defaults_ctsm.xml:587`). So the previously deliberate "always-on tranr"
divergence is gone too. The SP strict gate does not run CN phenology and is
unaffected.

**cn_annual** (Fortran 2202 restart, 365 days, real SNICAR; rel = |J−F|/(1+|F|),
annual mean of monthly means):

| field | Fortran | before (#371) | after | rel before | rel after |
|---|---|---|---|---|---|
| LEAFC (gC/m2) | 9.8907 | 9.3737 (−4.8 %) | 9.8877 (−0.03 %) | 4.8e-2 | **2.7e-4** |
| TLAI | 0.066109 | 0.062631 | 0.066081 | 3.3e-3 | **2.7e-5** |
| GPP (gC/m2/s) | 2.4542e-6 | 2.4151e-6 (−1.6 %) | 2.4346e-6 (−0.8 %) | 3.9e-8 | 2.0e-8 |
| NPP (gC/m2/s) | 1.6903e-6 | 1.6760e-6 (−0.8 %) | 1.6728e-6 (−1.0 %) | 1.4e-8 | 1.8e-8 |
| TOTVEGC (gC/m2) | 374.37 | 374.22 | 374.14 | 4.0e-4 | 6.0e-4 |
| EFLX_LH_TOT (W/m2) | 8.3172 | 8.2654 | 8.3148 | 5.6e-3 | 2.6e-4 |
| FSH (W/m2) | 39.644 | 39.587 | 39.627 | 1.4e-3 | 4.0e-4 |
| TG (K) | 273.62 | 273.64 | 273.62 | 6e-5 | <1e-5 |

On the tree patch alone, LEAFC is now within 0.3 % in every month (it was −25 % by
December). The remaining grass-patch differences are onset/offset timing in April and
October, which predate this change.

**Cold-start 2-year** (#362 harness vs `fortran_cold2y`, annual means):

| field | year | Fortran | before (#371) | after |
|---|---|---|---|---|
| TG (K) | 1 | 273.10 | 273.31 (+0.21) | 273.04 (**−0.06**) |
| TG (K) | 2 | 273.78 | 274.27 (+0.49) | 273.75 (**−0.03**) |
| TOTVEGC (gC/m2) | 1 | 175.04 | 177.20 (1.2e-2) | 177.57 (1.4e-2) |
| TOTVEGC (gC/m2) | 2 | 136.25 | 139.84 (2.6e-2) | 140.06 (2.8e-2) |
| TLAI | 1 / 2 | 0.4024 / 0.2982 | 0.3993 / 0.2871 | 0.4032 / 0.3030 |
| FROOTC (gC/m2) | 2 | 46.61 | 41.62 (−11 %) | 47.88 (+2.7 %) |
| GPP (gC/m2/s) | 2 | 2.2825e-6 | 2.1701e-6 (−4.9 %) | 2.3640e-6 (+3.6 %) |
| EFLX_LH_TOT (W/m2) | 2 | 8.5665 | 8.5061 | 8.6009 |

So the cold-start TG bias came from the evergreen canopy too. With leaf C stranded in
transfer, the tree's displayed canopy and fine roots were too small, and they shaded
and transpired too little. TOTVEGC moves slightly away (by 0.2 percentage points)
because carbon that used to sit in the transfer pool, which counts in TOTVEGC, is now
displayed leaf that turns over.

### Follow-up: grass phenology timing and the cold-start winter GPP

Two residuals were left after the evergreen fix: the grass patch's LEAFC in cn_annual
(+20 % in April, −5 % in October, read as onset/offset timing), and the cold start's
year-1 GPP (+15.8 %). The first was not in the model. The second was two cold-start
bugs, one of them a deliberate deviation from Fortran.

**1. The grass "timing" residual was a record offset in the comparison.** A startup
run writes a zero-length nstep-0 snapshot as history record 1 (time 0, bounds [0, 0];
`lnd.log`: "Writing current time sample ... at nstep = 0"). The mean over model day d
is the record with time == d. `fortran_parity_cn_annual.jl` (and the per-patch
comparison scripts derived from it) indexed the records by day, so Julia day d was
compared with Fortran day d−1. Every field that changes fast looked one day early;
the grass onset jump in late April and the leaf drop in October are the largest such
changes. Aligned by the time stamp, the grass onset fires on the same day and the same
step as Fortran. The harness also took the orbit from the pre-advance calday, which is
one step behind Fortran's (`lnd_comp_nuopc.F90:890-891` takes calday and declin before
`clm_drv` at :903 and advances the clock after it at :921, with nstep 0 at the start
time, so the step a Julia loop step matches sees calday = start + k·dtime; `clm_run!`
already does this). That put the daylength-triggered grass offset one hour late. Both
are fixed in the harness. There is no phenology code change.

Grass patch (itype 12) LEAFC, monthly means, J/F − 1:

| comparison | Jan | Mar | Apr | May | Jun | Sep | Oct | Nov | annual |
|---|---|---|---|---|---|---|---|---|---|
| index-aligned (old harness) | −0.01 % | −0.01 % | **+19.82 %** | +4.55 % | −0.01 % | −0.01 % | **−4.67 %** | −0.01 % | +0.02 % |
| time-aligned, pre-advance orbit | 0.00 % | 0.00 % | 0.00 % | 0.00 % | 0.00 % | 0.00 % | +0.21 % | 0.00 % | +0.02 % |
| time-aligned, post-advance orbit (this branch) | 0.00 % | 0.00 % | 0.00 % | 0.00 % | 0.00 % | 0.00 % | 0.00 % | 0.00 % | 0.00 % |

The gridcell numbers in the cn_annual tables above are annual means of index-aligned
monthly means. The one-day shift barely moves those annual means, but the monthly
"timing" rows were artifacts.

**2. The cold-start winter GPP: frozen-soil PHS on a cold-start canopy.** Per-patch
daily and hourly history from the pristine ctsm5.3.012 exe (130-day cold start, extra
`hist_fincl2` tape with VEGWP, KROOT, HK, SMP, TSOI, SOILLIQ, SOILICE, FPSN, FCTR,
TSAI, ...) against the same fields dumped from the Julia harness. The whole year-1
excess is in January to April. Fortran's needleleaf GPP is about 0 from day 13 to day
96, because the whole root zone is frozen (SOILLIQ 0.01 kg/m2 per layer). Every
k_soil_root then sits at the 1e-16 conductance floor, the plant water potential falls to
−4e5 to −2e6 mm, and the plc stress closes the stomata. In Julia, VEGWP stayed at the
soil potential (−1.29e5 mm), btran stayed near 1, and the tree photosynthesised all
winter (March GPP 1.8e-6 against Fortran's 1.9e-9). Two causes:

- **`getvegwp` threshold** (`photosynthesis.jl` `_getvegwp` and the host `getvegwp!`;
  Fortran `PhotosynthesisMod.F90:5018-5023`). Fortran takes the mean of
  (smp − grav2) only when `sum(k_soil_root) == 0`. Otherwise it solves
  x(root) = (Σk(smp − grav2) − qflx)/Σk. The port used `abs(sum_ksr) <= 1e-12` and,
  below that, a k-weighted mean of (smp − grav2). That is the qflx → 0 limit, so
  vegwp was pinned at the soil potential. The threshold was added in d0ad46a, when
  the below-threshold branch was an *unweighted* mean that the dry bedrock dragged to
  −4e7 mm. #228 later replaced that branch with the weighted mean and made soilflux
  the per-layer sum. The threshold's own reason is gone with both of those: the
  floors are Fortran's own `max(…, 1e-16)` conductances (`:3086-3087`), not an AD
  artifact, and the cn_annual harness header already records that the Fortran-literal
  `== 0` produces no NaN. Now Fortran-literal: soilflux telescopes to the demand, the
  HydStress sink debits the same flux, and the water balance closes.
- **SP canopy seeded on a CN cold start** (`cold_start.jl` `init_satellite_phenology!`).
  CTSM fills the canopy state from the surface-data monthly LAI/SAI only for SP,
  FATES-SP or dry deposition (`clm_initializeMod.F90:461-486,697-707`). Under `use_cn`,
  tlai/tsai/elai/esai/htop/hbot keep `CanopyStateType::InitCold`'s zeros (`:529-536`).
  The port seeded them on every cold start. CNVegStructUpdate then started from the
  surface SAI (2.77 for Bow's needleleaf; Fortran's CN value is 0.5) as `tsai_old`, and
  relaxed it only over months. Through the first winter the tree carried 5.5× the stem
  area, 2.9× the root conductance (through the root-area index) and 5.5× the xylem
  conductance. `cold_start_initialize!` now takes `use_cn` and seeds only the SP
  interpolation arrays when it is set.

**AD consequence.** With vegwp now supply-limited on frozen soil, the stomata shut
completely (bsun = bsha = 0, net assimilation exactly 0). The Medlyn gs quadratic then
has an exact double root (b² = 4ac), and `quadratic_solve`'s `sqrt(0)` gave an
infinite Dual partial, so Inf·0 = NaN in every derivative. `test_ad_robustness.jl`
"Spring (285K)" and the scenarios after it failed (AD NaN, FD finite). Negative
round-off was already clamped to 0 and hit the same `sqrt(Dual(0))`. `quadratic_solve`
now takes sqrt = 0 as a constant when the discriminant is ≤ 0. The value is identical,
and every Float64 result is bit-identical.

Applied first, the SAI fix brings TSAI and SNOCAN to Fortran's values, but winter GPP
stays high. The getvegwp fix closes it. With both, the needleleaf
winter GPP and transpiration track Fortran day by day (day 9: GPP 1.166e-7 vs 1.158e-7;
FCTR 1.09e-3 W/m2 vs QVEGT 4.34e-10 mm/s, the same flux).

**Cold-start 2-year** (#362 harness vs `fortran_cold2y`, annual means, both runs with real
SNICAR; "before" is main after #372):

| field | year | Fortran | before | after |
|---|---|---|---|---|
| GPP (gC/m2/s) | 1 | 1.7916e-6 | 2.0742e-6 (+15.8 %) | 1.8084e-6 (**+0.9 %**) |
| GPP | 2 | 2.2825e-6 | 2.3640e-6 (+3.6 %) | 2.3473e-6 (+2.8 %) |
| NPP (gC/m2/s) | 1 | 3.1971e-7 | 4.4299e-7 (+38.6 %) | 3.3436e-7 (**+4.6 %**) |
| NPP | 2 | 8.7454e-7 | 8.8176e-7 (+0.8 %) | 8.8417e-7 (+1.1 %) |
| TG (K) | 1 / 2 | 273.10 / 273.78 | 273.04 / 273.75 | 273.09 / 273.77 |
| TOTVEGC (gC/m2) | 1 | 175.04 | 177.57 (+1.45 %) | 175.35 (**+0.18 %**) |
| TOTVEGC | 2 | 136.25 | 140.06 (+2.80 %) | 137.08 (**+0.61 %**) |
| TLAI | 2 | 0.29817 | 0.30298 (+1.6 %) | 0.30083 (+0.9 %) |
| FROOTC (gC/m2) | 2 | 46.608 | 47.877 (+2.7 %) | 47.035 (+0.9 %) |
| EFLX_LH_TOT (W/m2) | 1 | 8.8868 | 9.2147 (+3.7 %) | 8.9201 (+0.4 %) |
| SMINN (gN/m2) | 1 | 0.0047905 | 0.003975 (−17 %) | 0.0047737 (−0.4 %) |
| TOTECOSYSC (gC/m2) | 2 | 1129.6 | 1134.6 | 1131.2 |

(The "after" run also uses the post-advance orbit in the cold harness.)

**cn_annual** (Fortran 2202 restart, time-aligned; annual mean of monthly means,
J/F − 1; "before" is main after #372 run through the same aligned comparison):

| field | Fortran | before | after |
|---|---|---|---|
| LEAFC | 9.927 | −0.03 % | −0.05 % |
| TLAI | 0.06635 | −0.03 % | −0.05 % |
| GPP | 2.463e-6 | −0.82 % | −0.98 % |
| NPP | 1.696e-6 | −1.08 % | −1.27 % |
| TOTVEGC | 374.5 | −0.06 % | −0.07 % |
| TG | 273.7 | +0.01 % | +0.00 % |

On the restart run the root zone is never fully frozen, so the threshold rarely
binds. The small move (winter GPP −10 to −14 % of a ~5e-8 flux, against −2 to −7.5 %
before) comes from the hours when it did.

**Strict gate** (5 sites, scored without the old BTRAN→BTRANMN alias): Bow 68/69,
Stillwater 68/69, Abisko 69/69, Baltimore 68/69, Savanna 69/69, the same as before.
Stillwater, Abisko, Baltimore and Savanna are bit-identical. The SP path never
touches the CN canopy seed, and their root zones never sit at the floor. Bow's cells
move by hundredths of a percent (FPSN annual −2.73 % → −2.87 %, BTRANMN +3.92 % →
+4.11 %, still its one failing cell; FCEV and SOILICE move toward Fortran).

**Still open** (both items are closed in the next follow-up). The October grass and gridcell GPP in cn_annual (−14 % / −9.7 %)
follows a soil-water BTRAN difference after a rain/snow event on days 279-294
(Fortran's daily-min BTRAN rises to 0.18-0.73 there, Julia's lags). That is
hydrology, not CN, and it is the same cell that fails Bow's SP gate (BTRANMN). Year-2
cold-start GPP is still +2.8 %, mostly in April-May (+36 % / +8 %, against −13 to −19 %
in Jan-Mar), and these months are essentially unchanged by this fix. Not traced yet;
a year-2 per-patch diff is the next step. Year-1 monthly GPP is now within 6 % in every
month (it was +64 % to +96000 % in Jan-Apr). The cold harness passes `doalb = true`
on every step, while a Fortran cold start skips the albedo (and so CNVegStructUpdate)
on its first steps. So Julia's step-1 TLAI is 0.71 where Fortran's is 0, a one-step
difference.

### Follow-up: the October BTRAN lag, Bow's BTRANMN cell and the year-2 spring GPP

The two residuals left open above were four separate problems: one in the
cn_annual harness and three in the model. Closing the Bow strict gate and the AD check
afterwards took three more (5-7, at the end of this section). None was in the soil-water physics itself,
although the October residual looked like it.

**Method.** A second history tape on the pristine ctsm5.3.012 exe, with one record per
step for a full year (`hist_nhtfrq = 1`, `hist_dov2xy = .false.`). It holds per-layer
H2OSOI, SOILLIQ, SOILICE, SMP, HK and TSOI, the infiltration and runoff chain, the
snow layers (SNO_Z, SNO_ICE, SNO_LIQH2O), INT_SNOW, VEGWP, KROOT, Vcmx25Z/Jmx25Z and
more. It was written for three runs:
- cn_annual, from the 2202 restart.
- The SP Bow gate year, a startup from the 2003-01-01 restart. Its h0 matches the
  gate reference to 0.3 % in BTRANMN.
- The 2-year cold start.

The same fields were dumped per step from the Julia harnesses and aligned by step:
Julia step k is Fortran record k+1, because record 1 is the nstep-0 snapshot.

**1. The cn_annual harness used the SP namelist calibration on a CN reference.**
`fortran_parity_cn_annual.jl` passed `int_snow_max = 3113.2227` and
`baseflow_scalar = 0.0022119554`, which are the Bow SP calibration. The cn_annual
reference is a `-bgc bgc -use_case 2000_control` case. Its `lnd_in` has the CTSM
defaults: `int_snow_max = 2000`, `baseflow_scalar = 0.001`, and no `h2osfcflag`, so the
default 1 (the harness passed 0).

`int_snow_max` caps `int_snow` in `FracSnowDuringMelt`
(`SnowCoverFractionSwensonLawrence2012Mod.F90:263`). That function sets the melt term
of the snow compaction (`SnowHydrologyMod.F90:2003-2035`). Bow's `int_snow` is about
9000 mm, so the cap binds. The chain of effects:
- The first quantity to diverge was snow-layer thickness on day 262 (step 6283), at a
  melt step, while SWE was still equal.
- The late-September snowpack then did not melt out on day 280 as Fortran's does.
  About 5 mm SWE survived into the October snowfalls.
- The top soil stayed about 0.03 m³/m³ drier through the rain on days 290-296.
- The grass daily-minimum BTRAN and GPP followed. This was the "lag".

The harness now uses the reference's values. This was not a model bug, but it was the
whole October residual.

**2. `e_ice` never reached the soil-water solve.** CTSM reads the ice-impedance factor
`e_ice` from the params file into three modules:
- `SoilHydrologyMod.F90:196`.
- `SoilWaterMovementMod.F90:104`, used by `IceImpedance` (`:2144`) for `hk_l` and the
  Richards solve.
- `InfiltrationExcessRunoffMod.F90:96`, used for `qinmax` (`:299`).

The port read it only into `soilhydrology_params`. `SoilWaterMovementConfig` had its
own `e_ice = 6.0` (`soil_water_movement.jl`), and `clm_drv!` built the config without
passing the params value. `infilt_excess_params.e_ice` was never set either. The stock
clm5 params file has `e_ice = 6`, so the bug was invisible everywhere except Bow, whose
calibrated file has `e_ice = 3.417`.

How it showed on Bow's frozen soil:
- Liquid, ice and temperature were identical step by step, but Julia's `hk_l` was
  30-37 % low in every layer.
- The PHS soil conductance, `min(hksat, hk_l)/(1e3·r_soil)`, followed it down.
- The leaf water potential fell further each afternoon (−5.94e5 mm against
  −5.73e5 mm).
- Transpiration, photosynthesis and the daily-minimum BTRAN fell 10-20 % in Jan-Mar.

That is the winter half of Bow's BTRANMN cell and of the cn_annual winter GPP deficit.
`clm_drv!` and the reverse driver now pass `e_ice = soilhydrology_params.e_ice`, and
`read_params!` now sets `infilt_excess_params.e_ice`. After the fix, winter `hk_l`,
`vegwp`, transpiration and FPSN agree hour by hour to 4 digits.

**3. The LUNA restart was incomplete.** With the winter fixed, Bow's BTRANMN got
*worse*: +7.3 % for the year, with April, May and June at +13 %, +52 % and +14 %. The
winter bias had been cancelling a spring one. The cause was the grass vcmax25.
Fortran's held at the restart's 125.5. Julia's dropped to the 30.0 seed within two
weeks.

LUNA skips a day when a patch has no leaf area or no leaf N, and falls back to
`vcmx25_z_last_valid_patch` (`LunaMod.F90:491-495`). CTSM restarts that pair
(`PhotosynthesisMod.F90:1078-1087`), but `read_fortran_restart!` read only `vcmx25_z`
and `jmx25_z`. So the dormant grass reverted to the `InitAllocate` seed (30/60). As a
result:
- Its spring photosynthesis was a third of Fortran's.
- Its transpiration was low, so its leaf water potential and BTRAN were high.

The reader now also loads `vcmx25_z_last_valid_patch`, `jmx25_z_last_valid_patch`,
`pnlc_z`, `enzs_z` and `gpp24`, the rest of the `use_luna` restart set. This affects
every restart-driven run with a patch that is leafless for part of the year.

**4. `tempavg_t2m` was a one-day relaxation, not the annual mean.** In Fortran,
`CNPhenologyClimate` accumulates `tempavg_t2m += t_ref2m · fracday/dayspyr`
(`CNPhenologyMod.F90:579`). That running sum reaches the annual mean of the 2-m
temperature after one year. `CNAnnualUpdate` copies it to `annavg_t2m` and resets it
(`CNAnnualUpdateMod.F90:86-87`). `annavg_t2m` sets the seasonal-deciduous onset
threshold, `crit_onset_gdd = exp(4.8 + 0.13·(annavg_t2m − tfrz))` (`:1019`).

The port had `tempavg += fracday·(t_ref2m − tempavg)`, an exponential relaxation with
a one-day time constant. On the cold start:
- Julia's year-2 `annavg_t2m` was the late-December temperature, 258.7 K, where
  Fortran's was about 270 K.
- The grass onset threshold fell from about 80 GDD to about 18 GDD.
- The grass leafed out on day 468 instead of day 485.

That was the whole year-2 April (+36 %) and May (+8 %) GPP excess. Restart runs were
unaffected in their first year, because `annavg_t2m` comes from the restart.
`phen_climate!` now accumulates the Fortran sum.

The cold harness passes `doalb = true` on step 1. That does not matter at this level of
agreement: with the four fixes, year-1 and year-2 monthly GPP are within 0.4 % in every
month except year-1 January (−3.9 %, a winter flux of 2e-8).

**cn_annual** (2202 restart, time-aligned monthly means, J/F − 1; "before" is
`origin/cn/grass-phenology-timing`, "after" is this branch with the corrected harness and all seven fixes):

| field | | Jan | Feb | Mar | Apr | May | Jun | Jul | Aug | Sep | Oct | Nov | Dec | annual |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| grass GPP | before | −22.9 % | −36.0 % | −24.3 % | −0.65 % | −0.32 % | −0.27 % | −1.01 % | −0.37 % | +0.05 % | **−14.1 %** | +0.14 % | −0.27 % | −0.88 % |
| | after | +0.21 % | +0.16 % | −0.01 % | −0.01 % | +0.00 % | −0.01 % | +0.01 % | +0.10 % | −0.01 % | **−0.01 %** | −0.01 % | +0.01 % | +0.00 % |
| tree GPP | before | −9.25 % | −13.9 % | −10.6 % | −1.01 % | +0.02 % | −2.02 % | −1.31 % | −0.60 % | −0.80 % | −1.25 % | −0.90 % | −1.58 % | −1.34 % |
| | after | +0.07 % | +0.08 % | +0.02 % | −0.01 % | −0.01 % | −0.02 % | −0.00 % | +0.05 % | −0.01 % | −0.01 % | −0.02 % | −0.03 % | −0.00 % |
| gridcell GPP | before | −9.64 % | −14.3 % | −10.8 % | −1.00 % | −0.25 % | −0.51 % | −1.05 % | −0.46 % | −0.21 % | **−9.65 %** | −0.86 % | −1.55 % | −0.98 % |
| | after | +0.08 % | +0.08 % | +0.02 % | −0.01 % | −0.00 % | −0.01 % | +0.01 % | +0.08 % | −0.01 % | **−0.01 %** | −0.02 % | −0.03 % | +0.00 % |
| gridcell NPP | before | +5.6 % | +10.3 % | +29.7 % | −1.29 % | −0.27 % | −0.63 % | −1.37 % | −2.18 % | −0.58 % | −17.6 % | +0.40 % | +0.21 % | −1.27 % |
| | after | −0.29 % | −0.08 % | −0.05 % | −0.02 % | −0.01 % | −0.01 % | +0.00 % | +0.47 % | −0.03 % | −0.03 % | +0.00 % | −0.00 % | +0.00 % |

Annual means, before → after:
- TG: +0.014 K → +0.002 K.
- EFLX_LH_TOT: −0.12 % → +0.00 %.
- TOTVEGC: −0.07 % → −0.00 %.
- LEAFC: −0.05 % → −0.00 %.

**Cold-start 2-year** (#362 harness vs `fortran_cold2y`, annual means; "before" is
`origin/cn/grass-phenology-timing`, "after" is this branch, all seven fixes, with the
cn_annual namelist values in the harness):

| field | year | Fortran | before | after |
|---|---|---|---|---|
| GPP (gC/m2/s) | 1 | 1.7916e-6 | +0.94 % | **−0.02 %** |
| GPP | 2 | 2.2825e-6 | +2.84 % | **+0.02 %** |
| NPP | 1 | 3.1971e-7 | +4.58 % | +1.53 % |
| NPP | 2 | 8.7454e-7 | +1.10 % | +2.57 % |
| TLAI | 2 | 0.29817 | +0.89 % | +0.04 % |
| TOTVEGC | 2 | 136.25 | +0.61 % | +0.41 % |
| TG | 1 / 2 | 273.10 / 273.78 K | −0.01 / −0.01 K | +0.00 / −0.00 K |

Year-2 monthly GPP, before → after:

| | Jan | Feb | Mar | Apr | May | Jun |
|---|---|---|---|---|---|---|
| before | −13.1 % | −19.3 % | −14.4 % | +36.2 % | +8.1 % | −4.2 % |
| after | +0.06 % | +0.03 % | +0.09 % | +0.10 % | −0.02 % | +0.00 % |

The harness namelist change alone moved year-2 GPP only from +2.84 % to +2.97 %. The
e_ice fix closed Jan-Mar, and the `tempavg_t2m` fix closed April-May.

**Strict gate after the four fixes alone** (scored with `scripts/parity_gate_diff.py`,
without the old BTRAN → BTRANMN alias): Bow 66/69. BTRANMN passed (+4.11 % → +0.37 %),
but QRUNOFF, QOVER and QFLX_INFL failed daily nRMSE (0.0045 → 0.13), because the
June melt-out pulse ran off one hour earlier, across midnight. The other four sites
were unchanged. Two more fixes, below, close this.

#### Closing the Bow runoff cells and the AD near-freezing check

**Use the branch run as the reference, not a startup.** The hourly tapes above came from a
startup-from-restart Fortran run. A startup executes a real physics step at nstep 0:
its first record already holds a stepped state, T_GRND 259.90 → 257.79 K. So the
startup run is one physics step ahead of the gate reference. The gate reference itself
is a continuous run. A `branch` run from the same 2003-01-01 restart
(`start_type = branch`, `nrevsn`, `brnch_retain_casename = .true.`) reproduces the
gate h0 bit for bit, and its record k is Julia step k. All hourly tracing below uses
that branch run.

**5. Photosynthesis ran on patches that had already converged** (`canopy_fluxes.jl`).
- CTSM calls PhotosynthesisHydraulicStress with the shrinking `filterp`
  (`CanopyFluxesMod.F90:1125-1137`; the filter is rebuilt at `:1435-1444`).
- The port called it on every exposed-veg patch on every iteration. A patch that had
  converged kept getting its `qflx_tran_veg` rewritten while another patch iterated.
- Its final transpiration then no longer matched the transpiration its `ecidif` cap
  and leaf energy balance were computed with.

On Bow on 21 April (step 2657) the effect ran as follows:
- The tree should evaporate its whole canopy store. Instead the rewritten transpiration
  left 1.7e-7 mm of `snocan`; Fortran's is exactly 0.
- The canopy-snow fraction goes as `snocan^0.15`, so that residue made a 0.6 % error in
  absorbed solar for the next four hours.
- That left a drift of about 1e-4 in the ground and canopy-air state.
- On 26 May 18:00 (step 3498) the drift tipped a bistable branch. The ground-evaporation
  conductance switches on the sign of `delq` (`:1253-1260`). `delq` was within a few
  1e-6 of zero, and its sign flipped from one iteration to the next:
  −2.3e-4, −1.4e-5, +2.9e-5, +2.2e-5, −2.5e-6, +3.0e-7, +2.7e-6.
- Fortran ended the iteration on the dew side, without the soil resistance (ground
  evaporation under the grass 9.8e-6 mm/s). Julia ended on the other side
  (0.7e-6 mm/s). That was the +0.2 mm SWE drift that moved the June pulse.

The in-loop photosynthesis calls now take `mask_exposedvegp .& active`.

**6. The SP LAI interpolation was one step late** (`clm_driver.jl`,
`satellite_phenology.jl`). interpMonthlyVeg uses `get_curr_date(..., offset=dtime)`
(`SatellitePhenologyMod.F90:234`). The port passed the current date, so every daily
change in LAI/SAI landed one step late, at 24:00 instead of 23:00. New
`sp_interp_date()` gives the date on the no-leap calendar.

**7. The Dual copy shared Int and Bool state with its Float64 source**
(`calibration.jl` `_calib_dual_copy`).
- This was the `test_ad_robustness` "Near-freezing (272K)" failure. Neither the AD
  derivative nor the kink was wrong.
- The copy converted Float64 arrays and passed everything else by reference, so the
  Dual instance and `inst_f64` shared `column.snl`, `imelt` and the other integer and
  boolean state.
- The first (AD-seeded) run combined two snow layers and wrote the new `snl` into
  `inst_f64.column.snl`. The two FD runs then started from that `snl` with the original
  two-layer arrays: an inconsistent state.
- So repeated runs from "the same" state gave SH = 6.80 W/m² (first run) and 3.04 W/m²
  (every later run). The T − 0.01 K run was not on a plateau; it was the second run.
- With the Bow e_ice, the warm-up now ended where a snow-layer combine happens. Before,
  it did not, which is why the bug did not show.
- Bits arrays (Int, Bool) and BitArrays are now copied. Every scenario now matches FD to
  0.0 %; Winter cold's SH had been at 6.47 %.

**Native restart.** `restart_io.jl` now writes and reads `vcmx25_z`, `jmx25_z`, the
last-valid pair, `pnlc_z`, `enzs_z` and `gpp24`, so a Julia-to-Julia restart keeps
the LUNA fallback state. A round-trip test was added.

**Strict gate** (5 sites, final; "before" is `origin/cn/grass-phenology-timing`):

| site | before | after | cells that changed status |
|---|---|---|---|
| Bow | 68/69 | **69/69** | BTRANMN FAIL → PASS (+4.11 % → +0.03 %, nRMSE 0.21 → 0.001) |
| Stillwater | 68/69 | 68/69 | none (SNOW_DEPTH still fails: annual −0.04 % → −2.49 %, nRMSE 0.161 → 0.131) |
| Abisko | 69/69 | 69/69 | none |
| Baltimore | 68/69 | **69/69** | SNOW_DEPTH FAIL → PASS (−7.19 % → +0.42 %, nRMSE 0.152 → 0.056) |
| Savanna | 69/69 | 69/69 | none |

Bow, before → after:

| cell | annual | daily nRMSE |
|---|---|---|
| FPSN | −2.87 % → −0.003 % | 0.066 → 0.0001 |
| QVEGT | −0.57 % → −0.001 % | 0.022 → 0.000 |
| H2OSNO | −0.04 % → −0.000 % | 0.0012 → 0.000 |
| QOVER / QRUNOFF / QFLX_INFL | −0.02 % → −0.00 % | 0.0045 → 0.0000 |

Stillwater: FPSN −0.24 % → −0.01 % (nRMSE 0.033 → 0.002), QVEGT −0.21 % → −0.01 %.
Abisko FPSN/QVEGT move toward Fortran by 0.04 percentage points.

**Still open.**
- Cold-start NPP is +1.5 % (year 1) / +2.6 % (year 2) while GPP matches to 0.02 %, so the
  difference is in autotrophic respiration. cn_annual (restart) NPP matches to 0.00 %, so
  it is cold-start specific.
- Stillwater's SNOW_DEPTH cell.
- `fortran_parity_cn_coldstart.jl` (single step) still uses the SP `INT_SNOW_MAX`. It
  has no effect at nstep 0 with no snow.

## Scorecard

| Category | Count |
|---|---|
| Flags audited | **~110** (16 driver-entry keywords + 79 `varctl` fields + `CLMDriverConfig` switches + the `use_flexibleCN` cascade) |
| MATCH | ~85 |
| MISMATCH — live physics | **0 remaining** after campaign #3 (`convert_ocean_to_land`, `glc_snow_persistence_max_days` closed there; see that section for the deliberate differences). Originally 6: `create_crop_landunit` **fixed** (M5), `use_luna` **fixed** (#267, M4), `use_hydrstress` **fixed** (PHS campaign, M3), `use_aquifer_layer` **fixed** (#259, re-confirmed here on the corrected tridiagonal solver), `baseflow_scalar` **fixed** (here, M2) |
| MISMATCH — inert (0 reads / guarded) | 12 (`nsegspc`, `nyr_forcing`, `h2osno_max`, `n_dom_landunits`, `n_dom_pfts`, 6× `toosmall_*`, `downscale_hillslope_meteorology`, `z0param_method`) |
| CONDITIONAL-MISMATCH under `use_cn=true` | 8 (the `use_flexibleCN` cascade, M6) — **closed in campaign #3**, together with `use_fun`, `leafresp_method`, `br_root`, `constrain_stress_deciduous_onset`, `initial_vegC`, `ncrit`/`nnegcrit` |
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
