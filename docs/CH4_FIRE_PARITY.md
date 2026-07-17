# Methane (CH4) and Fire (Li2016) — the first Fortran parity scorecards

**Status: both subsystems have now been diffed against Fortran CTSM for the first time.
Neither passes. Both were broken, in ways nothing in the repo could have detected.**

Nothing was tuned to match Fortran. Every fix below is a correction of a mis-ported
constant, a missing call, or a missing computation, verified line-by-line against the
CTSM source. The divergences that remain are reported, not hidden.

---

## 1. Why these had never been validated (the premise was wrong)

The task brief assumed the Bow Fortran reference ran with `use_lch4=.true.` and
`fire_method='li2016crufrc'`, and that the blocker was that the CLM5 restart does not
persist their prognostics.

**That is not what was happening.** The actual reference run's `lnd_in` had:

```
use_lch4    = .false.
fire_method = 'nofire'
```

So there were **0 CH4 variables in the dumps** and only `burndate`/`lfc` for fire, not
because the restart format drops them, but because **neither subsystem ever ran in the
Fortran reference at all**. There was nothing to diff against. The "GPU-validated"
claim in the README was internal-consistency only, and the README's
"Implemented, Not Yet Validated" was, if anything, generous.

A second, independent blocker sat on the Julia side (§4).

---

## 2. The new Fortran reference

New instrumentation (`scripts/validation/fortran_pdump/`, mirroring the `restFile_write_dump`
recipe from #216/#221):

* **`bgcdumpMod.F90`** — a new SourceMods module writing the CH4 + fire **diagnostics**,
  which are *not* restart variables and therefore appear in no existing dump:
  `conc_ch4/conc_o2` (sat+unsat), `ch4_prod/oxid/aere/ebul/tran` per layer, `o2stress`,
  `ch4stress`, the surface-flux decomposition, `finundated`, `totcolch4`; and
  `farea_burned`, `nfire`, `baf_crop`, `baf_peatf`, `fuelc`, `lgdp/lgdp1/lpop`, `fsr`, `fd`,
  the full fire C/N flux set, `somc_fire`, `burndate`, `lfc`.
* **Two new dump boundaries** in `clm_driver.F90`:
  * `after_ch4` — immediately after the `ch4()` call.
  * `after_fire` — immediately after `EcosystemDynamicsPreDrainage`.
    **NB: fire runs inside `CNDriverNoLeaching`, i.e. inside PreDrainage — NOT in
    PostDrainage**, contrary to the task brief. That is the first point at which the fire
    fluxes exist.
* `ch4Mod.F90` copied to SourceMods with its component pointers made `public` (a pure
  visibility change; no physics touched) so the dump module can read them.

Gated by `BGCDUMP_NSTEP_LO/HI`, so an unconfigured run is unaffected.

**Run** (`SYMFLUENCE_data/clm_bgc_spinup/bgc_ref_firech4`):

| | |
|---|---|
| config | `use_lch4=.true.`, `fire_method='li2016crufrc'`, `finundation_method='h2osfc'`, `use_cn/use_fun/use_nitrif_denitrif/use_hydrstress/use_luna=.true.`, `use_crop=.false.`, `soil_decomp_method='CENTURYKoven2013'` |
| start | **startup** from the converged 2202-01-01 BGC-spinup restart. A **branch is impossible**: it aborts with `Field missing from restart file: tempavg_agnpp` — the `use_lch4=.false.` spinup never wrote the CH4 restart fields. Startup + `finidat` cold-starts exactly those, and nothing else. |
| spin-up | 196 days with CH4 + fire ON *before* the window, so the CH4 state and the fire accumulators (`prec10`/`prec60`/`rh30`) are physically converged, not cold |
| window | nstep 4720..4888 (2202-07-16 16:00 → +7 days), `dtime = 3600 s` |

### `finundation_method`

The reference's namelist default was `'TWS_inversion'`, whose Prigent satellite-regression
stream file **does not exist** in the local inputdata tree. Switched to `'h2osfc'` — CTSM's
*prognostic* method (`finundated = frac_h2osfc`), which needs no stream file and is the
only one the Julia port implements. Both codes therefore run the same method.

### Lightning and population density — a synthetic stream

Li-family fire needs two streams. The real files
(`clmforc.Li_2012_...lnfm_Total...`, `clmforc.Li_2017_HYDEv3.2...hdm...`) are **absent**
from the local inputdata tree (fire was never on, so they were never downloaded) and could
not be re-fetched. Resolution — **the same approach PR #221 used for the missing CMIP6 ndep
file**: a synthetic stream, uniform in space and constant in time, driving **both** codes
(`scripts/validation/make_fire_streams.jl`):

```
lnfm = 4.57e-4 counts/km^2/hr    (= 4.0 flashes/km^2/yr, annual-mean total IC+CG density)
hdm  = 4.2     counts/km^2       (~9,300 people over the ~2,210 km^2 Bow-at-Banff basin)
```

Because the field is uniform and constant, ESMF's bilinear regrid and `shr_strdata`'s time
interpolation both return the stream value **exactly**, so Fortran's `forc_lnfm`/`forc_hdm`
are known analytically. Verified: Julia reads back `4.57e-4` / `4.2`, matching. The diff
therefore measures **fire physics**, not regrid fidelity. It reuses the existing
`fv0.9x1.25` ESMF mesh, so no mesh had to be generated.

> The stream is **synthetic**. It is not a real lightning/population dataset and this
> document does not claim it is. It is physically defensible for the site, and — because
> both codes see the identical file — it cannot bias the comparison.

### Non-vacuity — does Bow actually burn? Does it inundate?

**Fire: YES.** Fortran produces `FAREA_BURNED = 1.72e-8 s^-1`, `NFIRE = 9.90e-9`, and
non-zero per-patch fire C fluxes (`M_LEAFC_TO_FIRE = [0, 9.6e-8, 4.8e-7]`,
`M_DEADSTEMC_TO_FIRE = 1.6e-6`). A 0-vs-0 pass is **not possible** here.

> Honest caveat: `1.7e-8 s^-1` ≈ 0.54 yr^-1 — roughly a 2-year fire return interval, which is
> ~10× more fire than Banff really sees. That is a property of applying an **annual-mean**
> lightning rate year-round instead of concentrating it in summer storms. It is not a
> statement about either model, both codes see the identical forcing, and it makes the test
> **stronger** (a large, unambiguous fire signal to diff).

**Methane: partially.** Bow is a dry mountain site: `frac_h2osfc ≈ 0`, so
**`finundated = 0`** and the column is a net CH4 **sink** (`CH4_SURF_FLUX_TOT = -2.89e-12`),
which is the physically correct answer for an upland forest. This is *not* vacuous — CTSM
integrates **both** the saturated and unsaturated columns every step regardless of
`finundated`, so production (sat: `CH4_PROD_SAT = 5.7e-8`), oxidation (unsat:
`CH4_OXID_UNSAT = 5.0e-10`), aerenchyma, ebullition and both tridiagonal transport solves
all run with non-zero values and are all compared.

**What is NOT exercised at Bow: the high-inundation WETLAND regime** — the regime methane
actually matters for. `finundated ≡ 0` means the sat/unsat surface-flux *weighting* is
degenerate, and `ch4_dfsat_flux`, `CH4_EBUL_*` and `BAF_PEATF` are identically zero in
Fortran (marked *vacuous* below — they prove nothing). **Validating CH4 as a source
requires a peatland/wetland site** (e.g. MerBleue, already in the multi-biome scorecard).
That remains **unvalidated**.

---

## 3. Bugs found and FIXED — methane

### 3a. `CH4VarCon` defaults: three were the **opposite** of CTSM's

Checked line-by-line against `ch4varcon.F90`:

| flag | CTSM | CLM.jl (was) | effect |
|---|---|---|---|
| `anoxicmicrosites` | `.false.` (`:63`) | **`true`** | CH4 production **above the water table** switched on. Fortran's `CH4_PROD_UNSAT ≡ 0`; Julia's was not. |
| `ch4rmcnlim` | `.false.` (`:58`) | **`true`** | the CN / low-moisture limitation on SOM HR was being **removed** |
| `use_aereoxid_prog` | `.true.` (`:17`) | **`false`** | aerenchyma oxidation prescribed instead of prognostic |

These made the Julia methane model **a different model**. The Fortran reference confirms the
fix: `CH4_PROD_UNSAT` is **identically zero** in Fortran, exactly as `anoxicmicrosites=.false.`
requires. Fixed in `methane.jl`, with each default annotated with its `ch4varcon.F90` line.

### 3b. `rgasm` was **1000× too small** — every atmospheric boundary concentration was 1000× too large

`ch4Mod.F90:1843`: `rgasm = rgas / 1000` — where Fortran's `rgas` is `SHR_CONST_RGAS`
= **8314.47 J/K/kmol**. The `/1000` is the **kmol → mol conversion**.

`methane.jl:1570,2361`: `rgasm = RGAS / 1000.0` — but Julia's `RGAS` (`varcon.jl:29`) is
**already per-mol** (8.3145 J/K/mol). The `/1000` was applied **twice**.

Consequence: `c_atm = p / (rgasm · T)` was **1000× too large** for all three species.
Measured before the fix: `c_atm = [0.0563, 6921, 12.15]` (CH4, O2, CO2). After:
`[5.63e-5, 6.92, 0.0122]` — matching Fortran. `rgasm` also scales the ebullition gas volume,
so ebullition was simultaneously 1000× too small. Fixed at both sites.

### 3c. `finundated` was a frozen constant 0.1 — `CalcFinundated` was never ported

`ch4Mod.F90:1892-1920` is **two** steps: (1) the `CalcFinundated` call
(`ch4FInundatedStreamType.F90:295-317`), which for `finundation_mtd_h2osfc` is simply
`finundated(c) = frac_h2osfc(c)`; then (2) the snow-season hold + redox lag.

**Only step (2) was ported.** The Julia kernel clamped `finundated` against *itself* and
lagged it toward *itself* — a no-op fixed point that left every column pinned at its
cold-start `0.1` **forever**. The model ran "10 % inundated, always, everywhere".
`CH4VarCon.finundation_mtd_h2osfc` was a declared-and-never-read field.

Fixed: `CalcFinundated` ported (h2osfc route), `finundation_mtd` is now a real setting, and
the two unported inversion methods `error()` explicitly rather than silently doing nothing.

### 3d. `soil_bgc_carbon_flux_summary!` was a DEAD PORT — **all heterotrophic respiration was zero**

The single largest find. `soil_bgc_carbon_flux_summary!` (`types/soil_bgc_carbon_flux.jl`)
existed, was unit-tested, and was **never called from anywhere in `src/`** — both "call
sites" (`cn_driver.jl:1520`, `:1543`) were **comments**, and the docstring openly said
*"ported, not wired"*.

It is what vertically integrates `decomp_cascade_hr_vr → decomp_cascade_hr_col` and then
sums that into `somhr` / `lithr` / `cwdhr` / `michr` / `hr_col` / `hr_vr`. Without it:

```
decomp_cascade_hr_vr_col  sum = 9.28e-6      (nonzero — decomposition works)
decomp_cascade_hr_col     sum = 0.0          (never integrated)
somhr = lithr = hr_col = hr_vr = 0.0         (everything downstream)
```

**Every heterotrophic-respiration diagnostic in CLM.jl was identically zero** — and since
`ch4_prod` reads `somhr`/`lithr`/`hr_vr`, **CH4 production was identically zero too**. This
is a bug well beyond methane; it silently zeroed soil respiration everywhere.

Wired at the Fortran position (inside `CNDriverSummarizeFluxes`). After the fix:
`somhr = 7.75e-7`, `lithr = 4.67e-7`, and `CH4_PROD_SAT` goes from `0.0` to `2.90e-9`
against Fortran's `2.23e-9`.

---

## 4. Bugs found and FIXED — fire

### 4a. The entire Li fire chain was a DEAD PORT — `_fire_active` could **never** be true

Every Li2014/2016/2021/2024 kernel exists and is faithful. But on the production path
(`clm_drv!` → `vegetation_facade.jl` → `cn_driver_no_leaching!`) the gate
`_fire_active` (`cn_driver.jl:1109-1116`) was **structurally false**, for two independent
reasons, either of which alone was fatal:

1. `CNVegetationConfig` had **no fire field at all**, and `_sync_driver_config!` copied 12
   flags, none of them `cnfire_method` — so `config.cnfire_method` could never leave
   `:nofire`.
2. The seven fire data structs (`CNFireLi2014Data`, `CNFireBaseData`, `CNFireConstData`,
   `CNFireParams`, `PftConFireBase`, `PftConFireLi2014`, `DgvsFireData`) were **never
   instantiated anywhere in `src/`** — no `CLMInstances` member, no init-cold. Only
   `test/test_cnfire_wiring.jl` built them, by hand.

And even with the gate forced open, fire had **no ignition source**: `forc_lnfm`,
`forc_hdm`, `gdp_lf_col`, `peatf_lf_col`, `abm_lf_col` were all zero-length `Float64[]`
with **no writer anywhere in `src/`**.

**Now wired**, gated so the default (`:nofire`) is bit-identical:
`cnfire_method` threaded through `CNVegetationConfig`/`CLMDriverConfig`; the seven structs
instantiated and init-cold in `CLMInstances`; `readParams_CNFire!` wired (verified:
`ignition_efficiency = 0.22`, `prh30 = 0.7` now read from the param file); the fire bundle
plus `fsat`/`wf`/`wf2`/`forc_t`/`forc_rain`/`forc_snow`/`kmo`/`kda`/`nstep` threaded through
the facade; a new `src/infrastructure/fire_streams.jl` reading lnfm/hdm; and `gdp`/`peatf`/
`abm` read from the surface dataset (verified: `gdp_lf = 40.0`, matching Fortran).

**The Li2016 scheme question, answered:** `fire_li2016.jl` **does** exist, and
`cnfire_method_symbol` maps `"li2016crufrc" → :li2016` correctly (`fire_factory.jl:33`).
The port is **not** Li2014-only. Confirmed by the scorecard: `nfire`/`farea_burned` in
`fire_li2016.jl:340-355` are **line-for-line identical** to `CNFireLi2016Mod.F90`.

### 4b. Five more real port bugs, found by *wiring* the chain

Turning fire on immediately exposed five further genuine bugs. Two of them are **not
fire-gated** — they were corrupting **every** CN run:

1. **`forc_q_not_downscaled_grc` was never written by the forcing reader.** Its only
   consumer is `atm2lnd_update_rh!`, which therefore computed **`forc_rh_grc ≡ 0` in every
   single CLM.jl run, ever.** Relative humidity was identically zero: `RH30`/`RH24` were dead,
   and every RH-dependent term (all of Li fire's fuel-combustibility factors) ran permanently
   bone-dry. Now set from `QBOT`; Bow's RH goes **0 → 44.9 %**.
2. **`prec10` / `prec30` / `prec60` / `rh30_patch` were allocate-NaN.** Fortran registers them
   with `init_value=0` and `InitAccVars` fills them before step 1. Julia's step-1 fire — *and*
   CNPhenology's rain-triggered stress-deciduous onset — read **NaN**. Now 0-initialised
   (`accum_runmean` still returns `val` verbatim at `nstep<=1`, so cold-start accumulation is
   unchanged).
3. **`cnveg_carbon_state_summary!`'s column `p2c` was stubbed** ("pending subgridAveMod"), so
   `totvegc_col` / `totc_p2c_col` stayed **allocate-NaN forever**. `totvegc_col` is a live
   input to the fire fuel load → **NaN `FUELC` → NaN burned area**. This alone would have made
   any fire run NaN.
4. **`fd_pft` off-by-one** (`fire_li2014.jl:414`, `fire_li2016.jl:161`): indexed
   `fd_pft[itype[p]]` where Fortran's array is **0-based** — while `fsr_pft`, on the adjacent
   line, correctly used `[itype[p]+1]`. Wrong PFT's fire duration, and index 0 (out-of-bounds
   under `@inbounds`) for bare soil. `FD_COL` 30240 → **82080**, which is now **exact** against
   Fortran.
5. **`forc_rh_grc` was derived only in the END-of-step accumulator block**, so fire read the
   *previous* step's RH (and 0 on step 1). Now also derived at the top of the step; the
   end-of-step call remains, so `RH30`/`RH24` accumulate bit-identically.

> Bugs 1 and 2 are **not gated by `cnfire_method`**. They are genuine fixes that change CN
> runs in general. That is the one behavioural change in this PR that reviewers should weigh
> deliberately — it is a bug-fix, not tuning, but "default output is bit-identical" holds for
> the **fire/CH4 wiring**, not for these two dead-input repairs (nor for the HR summary, §3d).

### 4c. Reference-run namelist: `&lifire_inparm`

The Bow case is built `-lnd_tuning_mode clm5_0_GSWP3v1`, whose `&lifire_inparm` **overrides
10 of the CTSM module defaults** that `CNFireConstData` carries — e.g.
`occur_hi_gdp_tree = 0.33` (not 0.39), `lfuel = 105` (not 75), `ufuel = 1050` (not 650),
`bt_min/bt_max = 0.85/0.98` (not 0.3/0.7), `pot_hmn_ign_counts_alpha = 0.010` (not 0.0035).

This is **configuration, not tuning**: both codes must be driven by the same namelist or the
comparison is meaningless (exactly as `br_root` / `leafresp_method` / the rooting profile are
already handled in `build_bow_inst`). Applying it took `LGDP` from `relerr 1.5e-1` to
**exactly 0.0**, which independently *validates* the `lgdp` port.

---

## 4d. Bonus: the CN gridcell balance check was reading GARBAGE MEMORY

Running the suite under `--check-bounds=yes` (as CI does) surfaced an eleventh bug, latent
since #221 and unrelated to fire/CH4 except that it blocked a green suite.

`_cnbal_begin_grc_kernel!` (`cn_balance_check.jl:266`) indexes three "dribbler" arrays
(`hrv_xsmrpool_amount_left`, `gru_conv_cflux_amount_left`, `dwt_conv_cflux_amount_left`)
**unconditionally** on the `!use_fates_bgc` branch. They default to `Float64[]`, and the
production caller (`vegetation_facade.jl:472`) **never supplies them** — the dribblers are not
ported. So the kernel indexed a **0-element array**, inside `@inbounds`:

* with `@inbounds` active (a normal run), this **silently read out-of-bounds memory** straight
  into `begcb_grc` — the *beginning* carbon mass of the CN gridcell balance check;
* only `--check-bounds=yes` turns it into the `BoundsError` that exposed it.

Fixed by substituting explicit zeros (via `similar`, so the GPU backends are preserved) —
"no dribbler" contributes `0.0`, which is exactly what Fortran's `InitGridcellBalance` sums
when the dribblers are empty. This is the `check-bounds` trap the repo has hit before, and it
is a reminder that `@inbounds` + an unported-and-therefore-empty array is a silent-corruption
pattern worth grepping for.

---

## 5. SCORECARD — Fire (Li2016), boundary `after_fire`, single-step oracle

Injected Fortran state → one Julia step → diff the same-step dump.
`relerr` = worst |J−F| over the profile, scaled by the profile's own |F| max.

| field | worst rel.err | \|F\| scale | verdict |
|---|---|---|---|
| `LGDP` | **0.0** | 2.33e-1 | **EXACT** |
| `LGDP1` | **0.0** | 4.42e-1 | **EXACT** |
| `LPOP` | **0.0** | 8.41e-1 | **EXACT** |
| `FSR_COL` (fire spread rate) | **0.0** | 2.72e-1 | **EXACT** |
| `FD_COL` (fire duration) | **0.0** | 8.21e+4 | **EXACT** |
| `LFWT` | **0.0** | 9.50e-1 | **EXACT** |
| `WTLF` | **0.0** | 9.50e-1 | **EXACT** |
| **`NFIRE`** | **6.1e-16** | 1.06e-8 | **EXACT** (was 1.24) |
| **`FAREA_BURNED`** | **1.3e-15** | 1.61e-8 | **EXACT** (was 4.76) |
| `FUELC` | 1.7e-06 | 2.39e+3 | near-exact |
| `ROOTC_COL` | 2.5e-07 | 1.65e+2 | near-exact |
| `M_LEAFC_TO_FIRE` | 8.5e-06 | 4.46e-7 | near-exact (was 4.758) |
| `M_LIVESTEMC_TO_FIRE` | 9.9e-05 | 3.84e-9 | near-exact (was 4.757) |
| `M_DEADSTEMC_TO_FIRE` | 2.3e-06 | 1.50e-6 | near-exact (was 4.758) |
| `M_LEAFC_TO_LITTER_FIRE` | 8.5e-06 | 8.92e-8 | near-exact (was 4.758) |
| `FIRE_MORTALITY_C_TO_CWDC` | 1.4e-03 | 6.60e-6 | near-exact (was 4.749) |
| `FIRE_MORTALITY_N_TO_CWDN` | 1.4e-03 | 1.54e-8 | near-exact (was 4.749) |
| `BAF_CROP`, `BAF_PEATF`, `LFC`, `FBAC`, `FBAC1`, `CROPF`, `TROTR1/2`, `DTROTR`, `FUELC_CROP`, `SOMC_FIRE`, `M_FROOTC_TO_FIRE` | — | **0** | *vacuous* (F ≡ 0 — proves nothing; no crop, no peat, no tropical tree, no land-use transition at Bow) |

**21 / 29 fields within 1e-9** (up from 19). `NFIRE` and `FAREA_BURNED` are now
**bit-exact**; the fire C/N fluxes, no longer inheriting a burned-area error, fall
to their own near-exact injected-state residual (≤1.4e-3, same class as `FUELC`).

### Root cause — a scorecard forcing-alignment bug, NOT a fire-physics defect (RESOLVED)

The original ~1.6–4.8× `NFIRE`/`FAREA_BURNED` divergence was chased to a single term —
**`fire_m`** (the climate/fuel-moisture factor) in `nfire = ig/secsphr · fb · fire_m · lgdp` —
because `ig`, `fb` and `lgdp` are all bit-exact (synthetic uniform lnfm/hdm, exact `lgdp`,
saturated `fb`). Instrumenting `fire_m`'s three inputs at Bow:

* **`afuel = 0`** — the fuel load (`fuelc ≈ 2385`) sits below the 2500 gC/m² threshold, so
  the `rh30`-weighted branch carries **zero weight**. (This is why restoring the `rh30`
  accumulator from the restart, PR #231, did *not* move `NFIRE` — a genuine finding that
  correctly **ruled `rh30` out**.)
* **`btran_term = 1`** — the injected soil moisture (`H2OSOI_LIQ` → `h2osoi_vol ≈ 0.08–0.16`,
  `s ≈ 0.3`) makes the leaf-weighted `btran2 ≈ 0.5`, well below `bt_min = 0.85`, so the root-
  wetness factor saturates at 1 in **both** codes. (An earlier note blamed `btran2`; that was
  wrong — one hydrology step cannot lift `s` from 0.3 to near-saturation, and the injected
  and fire-time `h2osoi_vol` are identical.)

So `fire_m` reduces to **`arh^1.5`**, a pure function of the *instantaneous* `forc_rh`. Julia
computed `forc_rh = 41.8 %`; back-deriving Fortran's from its dumped `NFIRE` gave the value
Fortran's `fire_m` must have used: **`forc_rh_F ≈ 57.7 %`**. The `forc_rh` formula
(`100·q/qsat`) and the forcing file (`clmforc.2003.nc`) are identical between codes, so the
difference could only be *which forcing record was read*. It was:

> **The Fortran BGC-spinup datm recycles model-year 2202 to forcing-year `2002` (not `2003`)
> and delivers each record one hour behind the model timestamp.** Julia's fire scorecard read
> `clmforc.2003` at the model hour.

Verified across all 5 probe steps — the implied Fortran `forc_rh` matches **`clmforc.2002` at
(model_date − 1 h)** to within 0.04 % (residual = the QSat-vs-approx formula), and with that
forcing `fire_m` matches Fortran to ~0.3 %:

| step | model date | implied `forc_rh_F` | `clmforc.2002` (−1 h) |
|---|---|---|---|
| 4721 | 07-16 17:00 | 57.68 | **57.73** |
| 4763 | 07-18 11:00 | 66.81 | **66.85** |
| 4804 | 07-20 04:00 | 41.34 | **41.37** |
| 4846 | 07-21 22:00 | 33.66 | **33.69** |
| 4888 | 07-23 16:00 | 45.93 | **45.97** |

**Fix (this PR):** point only the *forcing read* at `clmforc.2002` at `(model_hour − 1)` while
leaving the model clock on the 2003 day-of-year proxy that drives orbital/daylength (added a
default-preserving `forcing_date` override to `run_one_parity_step!`). This is configuration —
matching the reference's datm stream, exactly like the `&lifire_inparm` namelist in §4c — not
tuning. `NFIRE`/`FAREA_BURNED` go from 1.24 / 4.76 to **6e-16 / 1e-15** (bit-exact).

**The Li2016 fire scheme is now validated end-to-end against Fortran, `farea_burned` included.**
The only forcing-dependent quantity in the fire dump was `NFIRE`/`FAREA_BURNED` via `forc_rh`;
every other field is injected state, which is why the misalignment hid until `fire_m` was
instrumented. The residual ≤1.4e-3 on the fire C/N *mortality* fluxes is a separate near-exact
injected-carbon-state artefact of the single-step oracle (same class as `FUELC`'s 1.7e-6), not
a fire-formula defect.

---

## 6. SCORECARD — Methane, boundary `after_ch4`, single-step oracle

After all four §3 fixes. CH4 prognostics injected from the dump.

| field | worst rel.err | \|F\| scale | verdict |
|---|---|---|---|
| `O2STRESS_UNSAT` | **0.0** | 1.0 | **EXACT** |
| `CH4_PROD_SAT` | 7.7e-02 | 5.67e-8 | close (was **1.0** — Julia produced *nothing*) |
| `CONC_O2_UNSAT` | 1.0e-01 | 3.33e+0 | close (was 1.7e+3) |
| `ZWT_CH4_UNSAT` | 1.3e-01 | 8.60e+0 | DIVERGES |
| `CH4_OXID_SAT` | 5.0e-01 | 1.04e-7 | DIVERGES |
| `QFLX_SURF_LAG` | 8.5e-04 | 3.48e-5 | DIVERGES |
| `CH4STRESS_UNSAT` | 8.5e-01 | 1.0 | DIVERGES |
| `CH4STRESS_SAT` | 9.4e-01 | 1.0 | DIVERGES |
| `O2STRESS_SAT` | 9.99e-01 | 1.0 | DIVERGES |
| `CONC_CH4_SAT` | 1.00e+00 | 4.91e-2 | DIVERGES (Julia ≈ 0) |
| `CONC_O2_SAT` | 1.00e+00 | 1.60e-1 | DIVERGES (Julia ≈ 0) |
| `CONC_CH4_UNSAT` | 1.33e+00 | 2.70e-5 | DIVERGES |
| `CH4_AERE_SAT` / `_UNSAT` | ~1.0 | 6.29e-8 / 8.81e-14 | DIVERGES (Julia ≈ 0) |
| `CH4_TRAN_SAT` / `_UNSAT` | ~1.0 | 6.48e-11 / 8.81e-14 | DIVERGES (Julia ≈ 0) |
| `CH4_SURF_AERE_SAT` / `_UNSAT` | ~1.0 | 1.67e-8 / 1.72e-14 | DIVERGES (Julia ≈ 0) |
| `CH4_OXID_UNSAT` | 1.8e+01 | 4.96e-10 | DIVERGES |
| `CH4_SURF_DIFF_UNSAT` | 3.5e+02 | 2.41e-10 | DIVERGES |
| `CH4_SURF_FLUX_TOT` | 3.5e+02 | 2.89e-12 | DIVERGES |
| `CH4_SURF_DIFF_SAT` | 1.5e+07 | 1.36e-12 | DIVERGES |
| `TOTCOLCH4` | 5.8e+00 | 5.23e-4 | DIVERGES |
| `CH4_PROD_UNSAT`, `CH4_EBUL_*`, `CH4_SURF_EBUL_*`, `CH4_EBUL_TOTAL_*`, `CH4_DFSAT_FLUX`, `FINUNDATED`, `FINUNDATED_LAG` | — | **0** | *vacuous* (F ≡ 0 — see §2 non-vacuity; `CH4_PROD_UNSAT ≡ 0` is the expected consequence of `anoxicmicrosites=.false.`) |

**11 / 33 fields within 1e-9** (8 of those are vacuous).

**Methane is NOT validated.** The four fixes moved it from *catastrophically* wrong
(concentrations 1000× out, production identically zero) to *the right order of magnitude*
for production and O2, but the transport/aerenchyma half is still badly wrong.

### Lead (a) — `grnd_ch4_cond_patch` frozen at cold-start `0.01` — FIXED (PR #234)

`bareground_fluxes.jl` and `canopy_fluxes.jl` each have a guarded write
(`if use_lch4 && length(grnd_ch4_cond_patch) >= p: grnd_ch4_cond_patch[p] = 1/raw`),
mirroring `BareGroundFluxesMod.F90:401` and `CanopyFluxesMod.F90:1097`, which in Fortran
write `1/raw` **directly into `ch4_inst%grnd_ch4_cond_patch`**. But the Julia driver called
**both** flux routines with `use_lch4=false` and an **empty** `grnd_ch4_cond_patch`
(`bareground_fluxes!` omitted the kwargs entirely; `canopy_fluxes_core!` passed `false, …,
Float64[]` in those positional slots). So the write never fired and `ch4.grnd_ch4_cond_patch`
stayed pinned at its cold-start `0.01` forever — starving the aerenchyma boundary conductance.

**Fixed:** the driver now threads `inst.ch4.grnd_ch4_cond_patch` into both calls under
`use_lch4` (Float64[] otherwise, so the non-CH4 and reverse-AD paths are byte-identical; the
kernel write itself is already `use_lch4`-guarded). Verified live: at Bow the field goes from
`[0.01, 0.01, 0.01, 0.01]` (frozen) to physical aerodynamic values `[0.0082, 0.0167, 0.0106,
0.01]` after one step. `test_methane.jl` 124/124; the ch4 path passes `--check-bounds=yes`
(the newly-active `@inbounds` write does not OOB).

> **Its effect on the Bow scorecard is masked** — `CH4_AERE_*` did *not* move, because at Bow
> `CONC_CH4_SAT` itself collapses to ≈ 0 in Julia (see the dominant lead below) and the
> aerenchyma flux is proportional to concentration. The fix is a **verified correctness
> repair**, not a validation of aerenchyma — that still needs a `finundated>0` wetland
> reference, which is **blocked** (see §7).

### Dominant remaining lead — the SATURATED column collapses to ≈ 0 in Julia (NEW, isolated this pass)

The single biggest signal in the §6 table: the **saturated** column's state collapses after
one step while the **unsaturated** column is near-right. The harness *injects* Fortran's
`CONC_CH4_SAT` / `CONC_O2_SAT` / `O2STRESS_SAT` as the IC, yet one Julia step drives them to
≈ 0 (`CONC_CH4_SAT` relerr 1.0, `O2STRESS_SAT` relerr 0.999 → Julia ≈ 0), whereas
`O2STRESS_UNSAT` is **exact** and `CONC_O2_UNSAT` is close (0.17). Because every `*_SAT`
transport/aerenchyma diagnostic (`CH4_AERE_SAT`, `CH4_TRAN_SAT`, `CH4_SURF_*_SAT`) is
proportional to the sat concentration, they are all ≈ 0 as a *consequence*. This — not
`grnd_ch4_cond` — is what dominates the sat-column divergence, and it is exactly the regime a
wetland (`finundated>0`, sat column heavily weighted) would stress hardest. **Prime suspect:
the saturated-column diffusion/transport solve (`ch4_tran` `sat=1` path).** Not fixed here:
isolating it responsibly needs the wetland ground truth (blocked), and CTSM integrates the sat
column every step *regardless* of `finundated`, so it cannot be chased from Bow alone.

### Other known remaining leads (reported, not fixed)
* **`fphr_col` is NaN** at the `ch4_prod!` read: `soil_biogeochem_decomp!` (`cn_driver.jl:699`)
  is called **without** `use_lch4`, so the `fphr` kernel never runs. Silently inert
  (guarded by `> 0`, and `NaN > 0` is false), which makes the `ch4rmcnlim` CN-limit removal
  a no-op.
* **The CH4 → decomp anoxia feedback is open-loop**: the `ch4!` call in `clm_driver.jl`
  (the `config.use_cn, false, false` positional args, currently ~line 2608) hardcodes
  `use_nitrif_denitrif=false, anoxia=false` — even though the reference config has
  `use_nitrif_denitrif=.true.`. So `o_scalar` is a constant 1.0 and `o2stress_unsat` never
  feeds back into decomposition. (Not a clean flag-flip: closing the loop needs the decomp
  step to read the ch4-written `o_scalar`, and validation needs the wetland reference.)
* **`site_ox_aere!`** (`methane.jl:1190`) is defined and never called from `src/` (its body
  was inlined into `_ch4aere_patch_kernel!`); **`qflx_surf_lag_col`** is computed and read by
  nothing.
* **No CH4 restart I/O** in `read_fortran_restart!` — CLM.jl cannot restart a methane run
  from a Fortran restart (the harness injects it manually).
* **No consumer of any CH4 output**: nothing writes them to history and nothing returns
  `ch4_surf_flux_tot` to the atmosphere.
* `cnveg_nitrogen_flux_summary!` and `soil_bgc_nitrogen_flux_summary!` remain
  **ported-but-unwired** (N-flux *diagnostics* only — no downstream physics reads them).

---

## 7. What is now true, and what is still not

**Validated:**
* The **Li2016 fire scheme is correctly ported** — the burned-area *formula*, the spread
  rate, the fire duration, all three suppression terms, **AND `NFIRE` / `FAREA_BURNED`
  themselves** are **bit-exact** against Fortran (once the scorecard's forcing was aligned to
  the reference's recycle-year `2002` + 1-hour datm lag; see §5). The fire C/N combustion/
  partitioning chain is exact-up-to a ≤1.4e-3 injected-carbon-state residual.
* The methane **atmospheric boundary condition** and **production** terms now agree with
  Fortran to the right order of magnitude, from a state where they were 1000× out and zero
  respectively.
* **Heterotrophic respiration is alive** for the first time.

**NOT validated (do not claim otherwise):**
* The whole methane **transport / aerenchyma / surface-flux** half. The `grnd_ch4_cond`
  wiring (lead a) is now correct and live, but the dominant divergence is the
  **saturated-column collapse** (above), still open.
* **Methane as a SOURCE.** Bow has `finundated ≡ 0`. The wetland regime — the one that
  matters — is untested and needs a peatland site.
* Fire's crop / peat / deforestation / tropical-tree branches: `baf_crop`, `baf_peatf`,
  `lfc`, `fbac`, `dtrotr` are all identically zero at Bow. **Vacuous.**

### CH4-as-source at MerBleue — attempted, BLOCKED on missing site assets (2026-07-17)

The plan was to generate the first `finundated>0` Fortran reference at MerBleue (peatland,
`finundation_method='h2osfc'` ⇒ `finundated = frac_h2osfc > 0`) and diff the Julia transport/
aerenchyma/ebullition half against it. **Phase 1 could not start**: the MerBleue domain assets
referenced by `scripts/parity_config.py:36` and `scripts/parity_run_domain.jl:121-125`
(`$SYMFLUENCE_DATA/domain_Peatland_MerBleue_Canada/` — surfdata, `clmforc.2017.nc`, the
converged `clm_peatland` restart) are **absent from this machine's `SYMFLUENCE_data`** (only
`domain_Bow_at_Banff_*`, `domain_Iceland_*`, `domain_Logan_River_*` are present; `ls`/`du`/
`stat` all confirm no MerBleue/Peatland tree — not FS flakiness). No other wetland site is
available, and with `finundation_method='h2osfc'` a `finundated>0` reference requires a
site whose surfdata + wet forcing produce `frac_h2osfc>0` (Bow's harness even runs
`h2osfcflag=0`). So the CTSM case could not be configured, let alone built/run.

**To unblock:** restore `domain_Peatland_MerBleue_Canada/` (surfdata + `data/forcing` +
`simulations/clm_peatland/CLM/*.clm2.r.*`), then follow §2/§8 (startup + `finidat`, spin CH4
ON ~196 days, dump `after_ch4`, **verify `finundated>0`** before diffing). The dominant
sat-column lead above is the first thing that reference would expose.

## 8. Reproducing

```bash
# Fortran reference (already generated; regenerate only if the instrumentation changes)
julia --project=. scripts/validation/make_fire_streams.jl <inputdata>/lnd/clm2/firedata
#   then: SourceMods bgcdumpMod.F90 + ch4Mod.F90 + clm_driver.F90 hooks; ./case.build;
#   run cesm.exe PLAINLY (never under lldb) with
#   BGCDUMP_NSTEP_LO=4720 BGCDUMP_NSTEP_HI=4888 PDUMP_NSTEP_LO=4720 PDUMP_NSTEP_HI=4888

julia --project=. scripts/fortran_parity_ch4.jl  5    # methane scorecard
julia --project=. scripts/fortran_parity_fire.jl 5    # fire scorecard
```
