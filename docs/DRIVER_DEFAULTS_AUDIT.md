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

### M1 — `use_aquifer_layer` (MISMATCH, live, highest severity)

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
stress is the CLM5 default photosynthesis path). CLM.jl defaults **`false`** in
`clm_initialize!`, `clm_run!` and `CLMDriverConfig`. Consequence: the default
driver uses the Soil-Moisture-Stress (BTRAN) path rather than PHS, changing
stomatal conductance, transpiration and GPP under water stress.

### M4 — `use_luna` (CONDITIONAL-MISMATCH, live)

```xml
<use_luna                                   >.true.</use_luna>
<use_luna                use_fates=".true." >.false.</use_luna>
<use_luna phys="clm4_5"  use_fates=".false.">.false.</use_luna>
```

CTSM `clm5_0` / non-FATES → **`.true.`**. CLM.jl defaults **`false`**.
Consequence: no LUNA photosynthetic-N acclimation, so `vcmax25` is the static
PFT value rather than the acclimated `vcmx25_z`.

### M5 — `create_crop_landunit` (MISMATCH, live)

```xml
<create_crop_landunit     use_fates=".false.">.true.</create_crop_landunit>
<create_crop_landunit     use_fates=".true." >.false.</create_crop_landunit>
```

CTSM keys this on **`use_fates`**, not on `use_crop`: for any non-FATES run it is
`.true.`. `clm_initialize.jl:127` sets `varctl.create_crop_landunit = use_crop`,
so a non-FATES, non-crop run gets `.false.` where CTSM gives `.true.`. Live
(4 reads in `src/`) — it affects subgrid landunit construction.

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
| `use_aquifer_layer` | `true` | **`.false.`** (derived: `clm5_0`→method 1→lbc 2 = `bc_zero_flux`) | **MISMATCH** | selects CLM4.5 Zeng-Decker-2009 + aquifer solver instead of CLM5 moisture-form + zero-flux; CTSM `endrun`s on this combined with `use_bedrock=true` | 16 scripts (see M1) |
| `baseflow_scalar` | `1.0e-2` | **`0.001`** (`lbc=2`); `1.d-2` only for `lbc=1` or `clm4_5` | **MISMATCH** | 10× drainage/baseflow rate | `parity_run_domain.jl` (per-domain) |
| `use_hydrstress` | `false` | **`.true.`** (`clm5_0`, non-FATES, `configuration="clm"`) | **COND-MISMATCH** | PHS off → BTRAN path; different gs/transpiration/GPP under stress | `parity_run_domain.jl`, `parity_run_domain_gpu.jl`, `probe_h2osfc_subdaily.jl` |
| `use_luna` | `false` | **`.true.`** (`clm5_0`, non-FATES) | **COND-MISMATCH** | no photosynthetic-N acclimation; static `vcmax25` | `parity_run_domain.jl`, `parity_run_domain_gpu.jl`, `fortran_parity_cn_coldstart.jl` |
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

## Scorecard

| Category | Count |
|---|---|
| Flags audited | **~110** (16 driver-entry keywords + 79 `varctl` fields + `CLMDriverConfig` switches + the `use_flexibleCN` cascade) |
| MATCH | ~85 |
| MISMATCH — live physics | **6** (`use_aquifer_layer`, `baseflow_scalar`, `use_hydrstress`, `use_luna`, `create_crop_landunit`, `convert_ocean_to_land`) + `glc_snow_persistence_max_days` (live via kwarg) |
| MISMATCH — inert (0 reads / guarded) | 12 (`nsegspc`, `nyr_forcing`, `h2osno_max`, `n_dom_landunits`, `n_dom_pfts`, 6× `toosmall_*`, `downscale_hillslope_meteorology`, `z0param_method`) |
| CONDITIONAL-MISMATCH under `use_cn=true` | 8 (the `use_flexibleCN` cascade, M6) |
| Namelist-vs-code-fallback disagreements found | **14** — CLM.jl copied the code fallback in *every* case |

## Not fixed here, and why

Three of the live mismatches are deliberately **documented rather than
changed**, because each is a model-configuration decision rather than a
mechanical correction, and flipping them silently would be the mirror image of
the bug this audit exists to catch:

- **`use_hydrstress` / `use_luna` (M3, M4).** Turning PHS and LUNA on by default
  switches the photosynthesis pathway for every caller. Both are already
  validated *when enabled* (see the PHS and LUNA parity work), so the physics is
  not in question — but the default flip changes every number in every test that
  does not pass them explicitly, and it should land as its own change with its
  own parity evidence, not buried in a defaults sweep.
- **The `use_flexibleCN` cascade (M6).** Eight coupled flags; enabling
  flexible-CN needs a validation campaign against a Fortran CN reference, not a
  one-line default change.

These are recorded as **open, evidenced findings**. The value of writing them
down without fixing them is that the next person to see
`use_hydrstress = true` in a harness will know it is a workaround for a known
wrong default, not a deliberate experimental setting.

## Method note

For each flag the check is three-layered, because agreeing with only one layer
is how #252 happened:

1. the **namelist default** with all its conditional variants,
2. the **code fallback** in the `*_inparm` initialiser,
3. whether CLM.jl's value is **actually read** anywhere (an unread field cannot
   cause a physics bug today but will silently mislead the first consumer).

A flag is only marked MATCH when the resolved namelist default under this port's
conditions agrees, *and* any disagreement with the code fallback is understood.
