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

> **UPDATE (§9): the assets were RESTORED and the reference was generated — but the site does
> not inundate (`finundated≡0`, water table ~1.9 m deep). The blocker is now physical
> hydrology, not missing assets. See §9 for the full round-2 result.**

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

---

## 9. MerBleue peatland — CH4-as-a-SOURCE, round 2: the site does NOT inundate (2026-07-17)

The plan (§7): now that `domain_Peatland_MerBleue_Canada` has been **restored**, generate the
first `finundated>0` Fortran CH4 reference at a peatland and finally diff the transport /
aerenchyma / **ebullition** half of the model — the regime Bow cannot reach.

**Phase 1 was executed end-to-end. The Fortran reference was generated. The GATE FAILED:
`finundated ≡ 0` across the whole window. MerBleue, as configured, does not inundate.**

### What was built and run
* Run dir `SYMFLUENCE_data/clm_bgc_spinup/merbleue_ref_ch4` (adapted from §2's Bow firech4
  recipe): the existing `symfluence_build` `cesm.exe` (already carries the `bgcdumpMod` /
  `ch4Mod` public-pointer / `clm_driver` `after_ch4` SourceMods — **no rebuild needed**),
  pointed at the MerBleue surfdata + `clmforc.2016/2017.nc` + the converged SP restart.
* Config = the §2 CH4 namelist (`use_lch4=.true.`, `finundation_method='h2osfc'`,
  `use_cn/use_fun/use_nitrif_denitrif/use_hydrstress/use_luna=.true.`,
  `soil_decomp_method='CENTURYKoven2013'`), **`fire_method='nofire'`** (fire is orthogonal to
  CH4-as-source and drops the synthetic lnfm/hdm stream dependency).
* **Startup + `use_init_interp` from the SP restart.** The restored MerBleue restart is
  **SP-only** — 114 vars, **no CN and no CH4 fields** (a `use_cn=.false.` spinup). No BGC
  restart for MerBleue exists anywhere in the tree, so CN + CH4 are `init_interp`
  cold-started, then spun **4720 steps** (start 2017-01-01 → mid-July 2017), dumping the
  `after_ch4` + `pdump` window nstep 4720..4888.

### The macOS-26 launch recipe (differs from §8 — banked)
The bare `./cesm.exe` run **aborts with `SIGTRAP` (exit 133)** — the macOS-26 xzone `brk #0x1`
heap-guard, caught fatally by gfortran's own signal handler. Three findings, in order:
1. `export MallocNanoZone=0` clears the **init-phase** trap (routes small allocs off the nano
   zone). Necessary but not sufficient.
2. `init_interp` writes `./init_generated_files/finidat_interp_dest.nc` — that subdirectory
   must be **`mkdir`'d** in the run dir first (else PIO `createfile` → "No such file or
   directory").
3. A *second* `brk` guard still fires in the first `clm_drv!` when run plainly — but **running
   the same binary under `lldb -b` (batch) passes the guards and the run completes to status
   0**, producing the full dump window. So for this config the working launch is `lldb -b`,
   the **opposite** of the §8/`macos26-lldb-xzone` note (which was for a warm-CN Bow config).
   Also needs the `fd.yaml` NUOPC field dictionary copied into the run dir.

### GATE RESULT — `finundated ≡ 0` (BLOCKER, physical this time)

| nstep | date | FINUNDATED | H2OSFC | ZWT | CH4_SURF_FLUX_TOT |
|---|---|---|---|---|---|
| 4720 | 2017-07-16 16:00 | **0.0** | 0.0 | 1.90 m | −4.4e-12 (sink) |
| 4800 | 2017-07-20 | **0.0** | 0.0 | 1.92 m | −4.3e-12 (sink) |
| 4888 | 2017-07-23 16:00 | **0.0** | 0.0 | 1.94 m | −5.0e-12 (sink) |

**Why it does not inundate:** the MerBleue surfdata models the bog as `PCT_NATVEG=100` (a
normal vegetated column, **not** a CLM `wetland` landunit — `PCT_WETLAND=0`), `FMAX=0.5`,
`20SL_8.5m` + `use_bedrock`, `baseflow_scalar=0.001`. That column equilibrates with the water
table ~1.9 m deep even in mid-July; the `h2osfc` surface-water store (which `finundated =
frac_h2osfc` reads under the `h2osfc` method) requires the water table to reach the **surface**
(saturation excess), which never happens. This is a genuine property of the SYMFLUENCE MerBleue
CLM configuration — the multi-biome scorecard's "MerBleue 69/69" parity agrees on this (dry)
hydrology between the two codes. Forcing inundation would mean re-engineering the site
(wetland landunit / zeroed drainage / raised WT) — a different, un-converged model, not "point
the reference at MerBleue." **`CH4_EBUL_*` and `CH4_DFSAT_FLUX` are identically zero — vacuous,
exactly as at Bow. CH4-as-a-source is NOT unlocked.**

### Phase 2 — the sat column reproduces the divergence, but MerBleue cannot isolate it

`scripts/fortran_parity_ch4_merbleue.jl` (single-step injection oracle, same pattern as the Bow
harness, parameterized for the MerBleue site/dump/date; disables the fatal CN veg balance check
since injecting partial state is not mass-conservative). Scorecard (worst rel.err over 3 window
steps): the SAT column is badly wrong — `CONC_O2_SAT` 44×, `CONC_CH4_SAT` 0.88, `O2STRESS_SAT`
0.38, `CH4_OXID_SAT` 1.68, `CH4_SURF_DIFF_SAT` 1e5 — while the UNSAT path is largely right
(`O2STRESS_UNSAT`, `CH4STRESS_SAT` **exact**). So the "saturated-column collapse" reproduces at
a **second, organic-soil site** — it is not Bow-specific.

**But MerBleue cannot cleanly attribute it to a transport bug**, for two reasons found here:
1. **Cold-CN confound.** No BGC restart exists ⇒ Julia's CN (and O₂ demand / production)
   cold-starts and diverges from the Fortran reference. `CH4_PROD_SAT` differs 27% for this
   reason alone; less decomposition ⇒ less O₂ consumed ⇒ `CONC_O2_SAT` accumulates — most of
   the 44× is this, not transport.
2. **At `finundated=0` the sat surface fluxes are zero-weighted and physically negligible**
   (`|F|` scale ~1e-10), so their huge rel.err (`CH4_SURF_DIFF_SAT` 1e5) is noise on a
   vanishing quantity, not a diagnostic signal.

### Hypotheses tested against the source (methane.jl ↔ ch4Mod.F90)
* **`jwt` not zeroed for the sat pass** (prime suspect from §6) — **REFUTED.** Fortran sets
  `jwt(c)=0` for `sat==1` (`ch4Mod.F90:2060`); Julia does the same via
  `_ch4o_jwt_zero_kernel!(jwt, mask_s)` (`methane.jl:2438`), over the correct soil mask. The
  `epsilon_t` / diffusivity branch (`j<=jwt` gas-phase vs `j>jwt` liquid) is therefore driven
  identically ⇒ the sat column correctly uses liquid diffusion.
* **AD-smoothing (`smooth_min/max`, the dimensional-`k` bug class) polluting the Float64 parity
  run** — **REFUTED.** `smooth_ad.jl` makes the Float64 methods **exact** (they smooth only for
  `ForwardDiff.Dual`); the parity run uses exact min/max.

### Verdict (honest)
* **CH4-as-a-source (`finundated>0` wetland regime): STILL BLOCKED**, now for a
  physical-hydrology reason — MerBleue as configured is a dry column. Needs a site whose
  surfdata + forcing actually produce `frac_h2osfc>0` (a true wetland landunit, or a
  water-table-shallow config), which is not available.
* **The saturated-column divergence is real and site-independent** (reproduced at MerBleue),
  but **isolating the transport-solve bug is not possible from either dry site** while the CN
  state is cold-started. The clean path — **inject Fortran's per-layer flux terms**
  (`CH4_PROD/OXID/AERE/EBUL_SAT`, all present in the `bgcdump`) so only the `ch4_tran` solve is
  under test, removing the cold-CN confound — is the documented next step. **No `methane.jl`
  fix was made: the honest "isolated to the sat column, jwt+smoothing refuted, clean isolation
  needs flux-injection or `finundated>0`" beats a tuned match.** `test_methane.jl` 124/124;
  `src/` untouched (default path byte-identical).

### Reproducing
```bash
# Fortran reference (in SYMFLUENCE_data/clm_bgc_spinup/merbleue_ref_ch4/):
#   cp fd.yaml + MerBleue datm/drv/runconfig; lnd_in = §2 CH4 namelist + MerBleue paths
#   + use_init_interp=.true. + fire_method='nofire'; runconfig start_ymd=20170101,
#   stop_option=nsteps stop_n=4900; mkdir init_generated_files
#   export MallocNanoZone=0 HWLOC_COMPONENTS=-opencl \
#          PDUMP_NSTEP_LO=4720 PDUMP_NSTEP_HI=4888 BGCDUMP_NSTEP_LO=4720 BGCDUMP_NSTEP_HI=4888
#   lldb -b -o run -o quit ./cesm.exe        # NB: lldb, not plain — see §9 launch recipe
julia --project=. scripts/fortran_parity_ch4_merbleue.jl 3   # MerBleue methane scorecard
```

---

## 10. The CH4 SOURCE-term residual — four mis-ported constants (oxidation kinetics + Henry's law)

After #237 fixed the saturated-column TRANSPORT (the `1.0e3` diffusivity placeholder), the
Bow scorecard still showed `CONC_CH4_SAT ~12%`, `O2STRESS_SAT ~1`, and huge `CH4_OXID_*`.
#237 isolated the remaining lead to the SOURCE recompute (production + oxidation). This pass
found it: **the oxidation kinetics were running 1000× (sat) / 100× (unsat) too fast, plus a
mangled Henry's-law table** — five constants, all verified line-by-line against `ch4Mod.F90`
`readParams` and `clm_varcon.F90`.

### The clean isolation — `scripts/ch4_oxid_source_injection.jl`
Oxidation (`ch4_oxid`) is a pure function of INJECTABLE before-step state (conc_ch4, conc_o2,
watsat, h2osoi_vol, smp_l, t_soisno, jwt) — nothing from the CN driver enters it, unlike
production (which reads warm-CN `somhr`/`hr_vr`, not dumped). So a clean oracle: inject
Fortran's before-step state, call `ch4_oxid!` directly, and diff.

**Key subtlety:** the dumped `CH4_OXID` is POST-competition — `ch4_tran` overwrites
`ch4_oxid_depth *= ch4stress` (or `*= o2stress`) (`ch4Mod.F90:3609,3623`). At Bow's sat column
`O2STRESS_SAT≈0`, scaling the raw rate down ~1000×. The isolation applies Fortran's own dumped
`min(ch4stress,o2stress)` to Julia's raw rate to compare against the dump. The residual was a
**uniform** ×1000 (sat) / ×100 (unsat) offset — a constant/param bug, not a per-layer formula
error. Deep layers (where the MM factor is linear) pinned the factors to *exactly* 1000 and 100.

### The bugs (all HARDCODED in `ch4Mod.F90 readParams`, file read commented out)
| param | Fortran | CLM.jl (was) | factor | affects |
|---|---|---|---|---|
| `vmax_ch4_oxid`   | `45e-6*1000/3600`      = 1.25e-5 | **0.0125**  | 1000× | sat + unsat oxidation |
| `vmax_oxid_unsat` | `45e-6*1000/3600/10`   = 1.25e-6 | **1.25e-3** | 1000× | unsat oxidation |
| `k_m_unsat`       | `5e-6*1000/10`         = 5.0e-4  | **5.0e-3**  | 10×   | unsat oxidation MM |

The defaults had dropped the `/3600` (per-hour→per-sec) and `/10` factors. `k_m`=5.0e-3 and
`k_m_o2`=2.0e-2 were already correct. In the linear MM regime `oxid ∝ vmax/k_m`, so the unsat
error = 1000/10 = **exactly 100×**, and the sat error = **exactly 1000×** — matching the probe.

### Plus the Henry's-law table (`clm_varcon.F90`), mangled — `src/constants/varcon.jl`
`k_h_inv = exp(-c_h_inv·(1/T − 1/kh_tbase) + log(kh_theta))` drives every CH4/O2 gas–liquid
partition (oxidation above-WT, aerenchyma, ebullition, and the `ch4_tran` diffusivity).
| const (CH4,O2,CO2) | Fortran | CLM.jl (was) |
|---|---|---|
| `C_H_INV`  | `[1600, 1500, 2400]`     | **[600, 1.3, 36]** (garbage) |
| `KH_THETA` | `[714.29, 769.23, 29.4]` | **[1600, 1500, 2400]** (= Fortran's c_h_inv!) |

### Before → after (Bow scorecard, worst rel.err over 3 window steps)
| field | baseline (post-#237) | + Henry | + oxid params (this PR) |
|---|---|---|---|
| `CH4_OXID_UNSAT` | 18.2   | 18.2   | **0.098** |
| `CH4_OXID_SAT`   | 1.67   | 1.68   | **0.167** |
| `CONC_CH4_SAT`   | 0.125  | 0.0085 | **0.0033** |
| `O2STRESS_SAT`   | 0.997  | 0.997  | **0.558** |
| `CH4STRESS_SAT`  | 0.938  | 0.938  | **0.238** |
| `TOTCOLCH4`      | 0.171  | 0.170  | **0.0045** |
| `CH4_SURF_FLUX_TOT` | 2.45 | 2.50  | **0.256** |
| `CONC_O2_SAT`    | 0.120  | **0.069** | 0.069 |
| `CH4_SURF_DIFF_SAT` | 0.596 | **0.065** | 0.065 |

The oxidation-kernel isolation (competition-compensated) went from **987 / 92** (sat/unsat) to
**0.095 / 0.143**, with deep layers now bit-matching Fortran (e.g. lev 20: 1.9857e-54 in both).

### What still remains (reported honestly, NOT tuned)
* **`CH4_PROD_SAT` = 7.7%, unchanged.** Production reads the warm-CN `somhr`/`lithr`/`hr_vr`,
  which are Julia-recomputed and NOT dumped, so this residual is a CN-state confound, not a
  clear production-kernel bug. The production kernel matches `ch4Mod.F90` line-by-line.
* **The saturated-column transport half** (`CH4_TRAN_SAT ~1`, `O2STRESS_SAT 0.56`,
  `CONC_O2_SAT/UNSAT ~6%`) is much improved but not closed. This is the transport/competition
  solve, still needing a `finundated>0` wetland ground truth (blocked, §7/§9) to isolate.
* `CH4_AERE_SAT` moved 0.14→0.61: a single-step sensitivity on a Bow-degenerate
  (`finundated≡0`, tiny `|F|`) diagnostic — not load-bearing.

`test_methane.jl` + `test_constants.jl` 198/198 under `--check-bounds=yes`; the five corrected
values are now asserted against their exact Fortran expressions. Default (non-CH4) path is
byte-identical — `C_H_INV`/`KH_THETA` and the `CH4Params` oxidation constants are used ONLY in
`methane.jl` under `use_lch4`.

### Reproducing
```bash
julia --project=. scripts/ch4_oxid_source_injection.jl 4721   # clean oxidation-kernel oracle
julia --project=. scripts/fortran_parity_ch4.jl 3             # full Bow scorecard
```

---

## 11. CH4-AS-A-SOURCE, UNBLOCKED — `finundation_mtd=TWS_inversion` at Bow (2026-07-18)

**The three-round blocker is gone: `finundated>0` was achieved, the saturated source
regime was diffed against Fortran for the first time, and a real mis-port was found and
fixed.** The sat-column divergence that §6/§9/§10 chased through jwt / smoothing /
diffusivity (#237) / oxidation constants (#239) had TWO remaining causes — a harness
forcing misalignment and a genuine aerenchyma porosity bug — both resolved here.

### The unblock — a wetter finundation method, not a wetter site
Rounds 1–2 (§7/§9) assumed CH4-as-a-source needed a site that PONDS (`frac_h2osfc>0`).
It does not. #238 ported `finundation_mtd=TWS_inversion`: `finundated = FWS_TWS_A*TWS +
FWS_TWS_B`, the Prigent **satellite regression**, which is DECOUPLED from the column's own
h2osfc pond. Evaluating the fitted coefficients (`finundated_inversiondata_0.9x1.25`) at
every migrated site's gridcell showed **several inundate**, including Bow itself: at Bow's
`(51.4,244.0)` cell `FWS_TWS_A=-3.957e-5, FWS_TWS_B=+0.0678`, so `finundated ≈ 0.045 > 0`
for the site's `TWS≈586 kg/m²`. (MerBleue does NOT: its cell needs `TWS>2242` — the h2osfc
result of §9 was not a fluke, the coarse-grid regression there is genuinely dry.)

**Why Bow, not a peatland:** Bow is the ONLY site with a converged BGC restart, so the CN
state (production, O₂ demand) is warm — avoiding the cold-CN confound that made the
MerBleue-as-source diff un-attributable (§9 Phase 2). finundated≈0.045 is modest but makes
the whole source regime NON-vacuous: the sat/unsat surface-flux partition is no longer
degenerate and `CH4_SURF_AERE_SAT` (saturated aerenchyma surface flux) is non-zero.

### The reference + harness
* Reference `SYMFLUENCE_data/clm_bgc_spinup/bow_ref_ch4tws`: IDENTICAL to §2's
  `bgc_ref_firech4` except `finundation_method='TWS_inversion'`. The `&ch4finundated` stream
  file was already in the namelist; the missing `_ESMFmesh_cdf5_130621` mesh was substituted
  with the `fv0.9x1.25` ESMFmesh (same 0.9x1.25 grid, exactly as the ndep stream does). No
  rebuild (the shared `cesm.exe` carries the bgcdump SourceMods). Launch = plain (the
  macOS-26 xzone `brk` guard is nondeterministic; it fired once, a retry ran clean to
  EXIT_CODE=0 — do NOT reach for `lldb -b`, which `-o quit`s at the first guard). Verified
  **`FINUNDATED=0.0446` across the whole nstep 4720..4888 window.**
* Harness `scripts/fortran_parity_ch4_tws.jl`: sets `finundation_mtd=TWS_inversion`, reads
  the REAL stream coefficients at Bow's cell, and back-derives Fortran's per-step `tws` from
  the dumped `FINUNDATED` (the bgcdump carries no TWS; adding it needs a Fortran rebuild) so
  `CalcFinundated` reproduces `FINUNDATED` exactly (validated: `4e-5`) and the transport diff
  is clean.

### Bug 1 (harness) — the datm misalignment, exactly as #233
grnd_ch4_cond (the aerodynamic boundary conductance = `1/raw`, which gates the aerenchyma
O₂/CH4 supply) is RECOMPUTED from the forcing wind/T each step; it diverged **11–78%
step-dependent**. The CH4 harness read `clmforc.2003` at the model hour, but the Bow
BGC-spinup datm recycles model-year 2202 to **forcing-year 2002, one hour behind** (§5).
Pointing the forcing at `clmforc.2002` at `(model_hour−1)` made grnd_ch4_cond **bit-exact**
(0.0123960 vs 0.0123960) and collapsed the whole UNSAT column to parity
(`CONC_CH4/O2_UNSAT`, `CH4_OXID_UNSAT`, `CH4_TRAN_UNSAT` all `≤1e-5`). CH4 prognostics are
injected so they barely feel the forcing — but the aerodynamic conductance does not, and it
gates aerenchyma. **The lesson of #233 repeats for methane: align the datm before blaming
the port.**

### Bug 2 (real mis-port, FIXED) — grass patches lost 2/3 of their aerenchyma
With grnd_ch4_cond exact, `CH4_AERE_SAT` still diverged **61%** (Julia ~2× low, layer-varying
factor). Isolation: grnd exact, `CONC_CH4_SAT` exact (0.4%), so the error was in the plant
tiller term `area_tiller·rootfr/z`. Bow's patches are 5% bare / 60% needleleaf tree / **35%
C3 grass**. `ch4Mod.F90 ch4_aere` gives **grass and crop** patches the full root porosity
(`poros_tiller=0.3`, Colmer 2003) but every other PFT only `0.3*nongrassporosratio`
(=0.099). **`methane.jl:_ch4aere_patch_kernel!` applied `nongrassporosratio`
UNCONDITIONALLY**, so the grass patch — the dominant CH4 conduit — had `1/nongrassporosratio
≈ 3×` too little aerenchyma. Fixed: `poros_tiller` is full for `itype >= nc3_arctic_grass`
(grasses 12/13/14 and crops ≥15 are exactly that set, reproducing CTSM's OR-condition),
reduced otherwise. `use_lch4`-gated (the kernel runs only under `use_lch4`); default path
byte-identical; guarded by a new `test_methane.jl` assertion (grass/tree aere ratio →
`1/nongrassporosratio`).

### Scorecard — Bow TWS_inversion, `after_ch4`, converged CN, aligned forcing
Both fixes applied (worst rel.err over 5 window steps; `|F|` = Fortran profile scale):

| field | before (this section) | after |
|---|---|---|
| `CH4_AERE_SAT`      | 0.610 | **4.7e-3** |
| `CH4_SURF_AERE_SAT` | 0.569 | **1.1e-3** |
| `CONC_O2_SAT`       | 0.069 | **2.5e-4** |
| `O2STRESS_SAT`      | 0.550 | **1.6e-2** |
| `CH4_OXID_SAT`      | 0.383 | **4.1e-2** |
| `CONC_CH4_SAT`      | 4.4e-3| **3.7e-4** |
| `CH4_SURF_FLUX_TOT` | 9.09  | **9.9e-3** |
| `TOTCOLCH4`         | 2.8e-3| **4.4e-4** |
| `FINUNDATED`        | —     | **4.2e-5** (validates the `CalcFinundated` TWS_inversion port) |
| `GRND_CH4_COND`     | 0.78  | **1.5e-5** |
| UNSAT column (`CONC_CH4/O2_UNSAT`, `CH4_OXID/TRAN/AERE_UNSAT`, stresses) | ~1 | **≤1e-5–1e-3 / EXACT** |

**The saturated CH4 source regime now agrees with Fortran to the single-step-oracle floor
(~1e-3–1e-4), the same near-exact class as the fire C/N fluxes (§5).** The "saturated-column
collapse" (Julia≈0) of §6 is GONE — it was the diffusivity (#237) + oxidation constants
(#239) + the datm alignment + the grass porosity, in that order of discovery.

### What remains (reported, NOT tuned)
* **`CH4_PROD_SAT` = 7.9%, unchanged** — production reads the warm-CN `somhr`/`lithr`/`hr_vr`,
  which are Julia-recomputed and NOT in the dump; a CN-state confound, not a production-kernel
  bug (the kernel matches line-by-line; same residual as §10).
* **Ebullition is still physically ZERO at Bow** (`CH4_EBUL_*`, `CH4_DFSAT_FLUX` ≡ 0): the
  code path runs, but Bow's saturated CH4 concentration stays below the bubble threshold (a
  dry, low-production column). Strong ebullition needs a wet, high-production column WITH a
  converged BGC restart — which does not exist in the migrated set. This is the one remaining
  source-regime term the current assets cannot exercise (honest, not vacuous-hidden).
* `CH4_DFSAT_FLUX` relerr ~1 on a `1e-14` scale — finundated is quasi-steady (0.0446±1e-4),
  so the dfsat redistribution is negligible; not load-bearing.

### Reproducing
```bash
# Fortran reference (already generated): SYMFLUENCE_data/clm_bgc_spinup/bow_ref_ch4tws/
#   lnd_in = bgc_ref_firech4 lnd_in with finundation_method='TWS_inversion' + the fv0.9x1.25
#   mesh for &ch4finundated; run ./run_bow_ch4tws.sh (plain; retry if the xzone brk fires).
julia --project=. scripts/fortran_parity_ch4_tws.jl 5    # Bow TWS_inversion source-regime scorecard
```

---

## 12. FIRE PEATLAND branch (`baf_peatf`) VALIDATED — a real wf2 double-count bug (2026-07-18)

**The peatland-fire term `baf_peatf` is now diffed against Fortran for the first
time, is BIT-EXACT, and finding it exposed & fixed a real port bug in the soil
`wf2` diagnostic.** (Backlog A2, peat row.)

### The unblock — a synthetic `peatf`, not a new site
`baf_peatf` is `F≡0` vacuous at Bow (and at *every* migrated site) because the
SYMFLUENCE surfdata pipeline never populates `peatf` — it is **0 in all 21
domain surfdatas**, MerBleue included. So the peat branch cannot be reached by
picking a peatland site. The unblock: `peatf` feeds **only** the CNFire modules
in CTSM (verified — `CNFire*Mod` + `CNVegStateType`/`FireDataBaseType`, nothing
else), so setting `peatf=0.5` in a COPY of Bow's surfdata is surgical — the
warm-CN Bow state is otherwise byte-unperturbed, and both CLM.jl and CTSM read
the SAME edited surfdata, so the diff measures the `_fire_peatland_kernel!`
physics, not the surfdata value (exactly like the synthetic uniform lnfm/hdm
streams of §2). Bow at 51.4°N > `borealat`=40 with unfrozen summer soil hits the
**boreal** peat branch.

### The reference + harness
* Reference `SYMFLUENCE_data/clm_bgc_spinup/bow_ref_firepeat`: identical to
  `bow_ref_ch4tws` (warm-CN, `fire_method='li2016crufrc'`) except `fsurdat`
  points at `surfdata_peatf.nc` (peatf=0.5). Ran **plain**, EXIT_CODE=0, the full
  nstep 4720..4888 window. `BAF_PEATF=1.67382e-10 /s` across the window (nonzero).
* Harness `scripts/fortran_parity_firepeat.jl`: single-step oracle; the pre-step
  hook injects the fire accumulators + `apply_bow_lifire_inparm!` + sets
  `cnfire_li2014.peatf_lf_col .= 0.5` to match. `baf_peatf` is a pure per-column
  function of injected before-step state (prec60, wf2, tsoi17, fsat) + peatf, so
  no cold-CN confound.

### The bug — `wf2` omitted the CTSM carryover double-count (FIXED)
First diff showed `BAF_PEATF` **7.7% high**, `NFIRE` bit-exact. Attribution probe
(`scripts/probe_firepeat_factors.jl`): the residual sat **exactly on the
`max(wf2,0)` kink** in the boreal formula `exp(-π·max(wf2,0)/0.3)`. At two of
three steps Julia's `wf2` was negative → factor 1.0 → EXACT; at nstep 4720
Julia's `wf2=+0.008` gave factor 0.919 while Fortran's implied `wf2≤0`.

Root cause: **CTSM `HydrologyNoDrainageMod.F90:640-699` computes `wf` (top 0.05m)
then `wf2` (top 0.17m) in one routine and does NOT reset `rwat/swat/rz` between
the two depth loops** — so the top-0.05m layers are accumulated TWICE into `wf2`.
Julia's `compute_wf!` computed `wf2` as an independent single accumulation to
0.17m, omitting the double-count. At Bow's dry summer top-soil the shallow layers
are drier than watdry (negative contribution), so double-counting them pushes
`wf2` from +0.008 to ≤0 — flipping the kink. Fix: new `compute_wf2!`
(hydrology_no_drainage.jl) replicates the carryover exactly (accumulate ≤0.05,
then continue ≤0.17 without reset). `wf2_col` is consumed **only** by the peat
fire branch (not history/restart), so the default (`:nofire`) path is
byte-identical. Guarded by a `compute_wf2! carryover` test (proves it differs
from the single accumulation) under `--check-bounds=yes`.

### Scorecard — Bow peat branch, `after_fire`, converged CN
| field | before fix | after fix |
|---|---|---|
| `BAF_PEATF`    | 7.70e-2 | **0.0 (bit-exact)** |
| `FAREA_BURNED` | 6.47e-3 | **1.3e-15** |
| `SOMC_FIRE`    | 7.70e-2 | **0.0** |
| `NFIRE`        | 6.1e-16 | 6.1e-16 |

**5/5 within 1e-9.** The peatland-fire branch is validated to the
single-step-oracle floor, same class as the fire C/N fluxes (§5).

### Still blocked — the tropical DEFORESTATION-fire term (`dtrotr`)
`baf_crop` (crop) and the tropical **deforestation** fire term remain
un-generatable on current assets:
* `dtrotr_col` accumulates `-dwt_smoothed[p]` **only** when `transient_landcover`
  is true (a `flanduse_timeseries` with `do_transient_pfts/lu`). **No migrated
  domain has a flanduse_timeseries** — so `dwt≡0` ⇒ `dtrotr≡0` even at
  Aripuanã/Leticia, where `trotr1+trotr2 > 0.6` (tropical BET ≈ 0.61/0.70) DOES
  meet the tree-cover gate. Plus those tropical sites have **no BGC restart**
  (cold-CN confound). Generating `dtrotr>0` needs BOTH a transient-land-use
  config AND a BGC spinup at a tropical site — a precise, un-generatable blocker.

### Reproducing
```bash
# Fortran reference (already generated): SYMFLUENCE_data/clm_bgc_spinup/bow_ref_firepeat/
#   cp bow_ref_ch4tws config; fsurdat -> a peatf=0.5 copy of Bow surfdata (surfdata_peatf.nc);
#   run ./run_bow_firepeat.sh (plain).
julia +1.12 --project=. scripts/fortran_parity_firepeat.jl 5   # Bow peat-branch scorecard
julia +1.12 --project=. scripts/probe_firepeat_factors.jl      # per-factor attribution
```

---

## 13. CH4 EBULLITION (`CH4_EBUL_*`) — the original blocker (2026-07-18)

**Ebullition (`CH4_EBUL_TOTAL_SAT/UNSAT`) could not be validated numerically on any
on-disk site — the Fortran reference itself was `≡0`.** (Backlog A1 residual.)
**This was RESOLVED — see §14.**

* The `_meth_ebul_kernel!` port is verified **line-by-line** against `ch4Mod.F90`
  `ch4_ebullition`: `bubble_f=0.57`, `vgc_max=vgc_min=0.15`, threshold
  `vgc > vgc_max·bubble_f = 0.0855`, `vgc = conc_ch4/watsat/k_h_cc·rgasm·T/pressure`.
* The Bow-TWS_inversion reference (§11, the ONLY site with `finundated>0` AND a
  converged BGC restart) gives **`CH4_EBUL_TOTAL_SAT ≡ 0`**: its peak
  `CONC_CH4_SAT ≈ 0.029 mol/m³` yields `vgc` well below the 0.0855 bubble
  threshold — Bow is a dry, low-CH4-production column. There is nothing to diff.
* Crossing the threshold needs a column that is SIMULTANEOUSLY (a) inundated
  (`finundated>0`), (b) high-CH4-production, and (c) has a converged BGC restart.
  **No migrated site provides all three:** Bow-TWS has (a)+(c) but not (b);
  MerBleue has the peat/production potential but neither inundates under
  TWS_inversion (its cell needs `TWS>2242` vs ~627, §11) nor has a BGC restart
  (SP-only, cold-CN). A full MerBleue BGC spinup (centuries of accelerated
  decomposition) is un-generatable in-session — but that is NOT required (§14).

---

## 14. CH4 EBULLITION — VALIDATED via a synthetic conc perturbation (2026-07-18)

**Ebullition is now diffed against Fortran for the first time and matches to the
single-step-oracle floor (`CH4_EBUL_TOTAL_SAT` rel.err 3.0e-5). No port bug — the
line-by-line-verified kernel is numerically confirmed.**

### The unblock — a threshold-crossing CONCENTRATION, not a new site
§13 assumed ebullition needs a physically high-production wetland (a MerBleue BGC
spinup). It does not — for a *kernel parity diff* it needs only that BOTH codes
compute ebullition from the **same** saturated CH4 concentration that happens to
cross the bubble threshold. Bow's peak `CONC_CH4_SAT ≈ 0.029` is only **~1.9×
below** threshold (vgc ≈ 0.045 vs 0.0855) — a *modest* gap, not orders of
magnitude. So the reference is Bow-TWS (which already has `finundated>0` + a warm
BGC restart, §11) with `CONC_CH4_SAT/UNSAT` scaled **×5** in a copy of the warm
restart. Both CLM.jl and CTSM read the identical scaled before-step concentration
from the `pdump`, so the diff isolates the `_meth_ebul_kernel!` physics — the exact
**synthetic-shared-input** technique of the peatf (§12, edited surfdata) and
lnfm/hdm (§2, synthetic streams) references, applied to a state variable.

### The reference + harness
* Reference `SYMFLUENCE_data/clm_bgc_spinup/bow_ref_ch4ebul`: a **short (6-step)
  startup continuation** from the ch4tws warm restart at `2202-07-24-14400`, with
  `finidat = finidat_ch4ebul_x5.nc` (`CONC_CH4_SAT/UNSAT ×5` via
  `scripts/probe_ch4_ebul_threshold.jl`'s edit). `start_ymd=22020724
  start_tod=14400`, `stop_n=6`, dump window wide (`get_nstep` resets to 0 for a
  startup → dumps `n0..n6`). Ran **plain**, EXIT_CODE=0. `CH4_EBUL_TOTAL_SAT =
  1.208e-6` at n0, decaying (7.9e-7, 5.7e-7, …) as transport dissipates the
  injected conc — n0 (pristine scaled conc) is the highest-signal step.
* The conc scale needed was verified in Julia FIRST (`probe_ch4_ebul_threshold.jl`,
  no Fortran box): with the model's real `watsat`, scale 1 → `CH4_EBUL_TOTAL_SAT=0`
  (Bow reality), scale 5 → `1.07e-6` (fires). This de-risked the box run.
* Harness `scripts/fortran_parity_ch4_ebul.jl`: identical to the §11 TWS harness
  (injects the scaled before-step conc + TWS_inversion finundated) but pointed at
  the ebul reference, `N0=0`, boundary `after_ch4`.

### Forcing alignment (the §11/§233 lesson, re-checked)
Unlike the recycled LONG run (datm 1 h behind, §5), this **fresh startup** aligns
the datm record to the model hour EXACTLY (offset **0**): at n0
`GRND_CH4_COND` is **bit-exact** (J=F=1.03706e-2, rel 1.7e-11) and `CH4_AERE_SAT`
is bit-exact at every step. (At n1+ the fresh-startup datm record drifts from a
constant-offset read for `GRND_CH4_COND` specifically, but this does NOT propagate
to the aerenchyma or ebullition fluxes — n0 is the fully-aligned reference step.)

### Scorecard — Bow ebul reference, n0, `after_ch4`, aligned forcing
| field | worst rel.err (n0) | \|F\| scale | note |
|---|---|---|---|
| `GRND_CH4_COND`     | **1.7e-11** | 1.04e-2 | forcing-sensitive path aligned (bit-exact) |
| `CH4_EBUL_SAT`      | **2.7e-5**  | 6.53e-6 | per-layer ebullition rate |
| `CH4_EBUL_TOTAL_SAT`| **3.0e-5**  | 1.21e-6 | **the term unlocked** |
| `CH4_SURF_EBUL_SAT` | **3.0e-5**  | 1.21e-6 | surface ebullition flux |
| `CONC_CH4_SAT`      | 4.3e-5      | 1.21e-1 | the injected (scaled) conc |
| `CH4_AERE_SAT`      | 1.5e-5      | 1.44e-7 | |

**Ebullition agrees with Fortran to ~3e-5** — the same single-step-oracle floor as
the §11 CH4 source fluxes and the §5 fire C/N fluxes. The residual traces to
`CONC_CH4_SAT` (4.3e-5), the one-step-recomputed warm-CN production confound
(§10/§11), NOT the ebullition kernel. The bubble threshold, vgc formula, and
release rate are validated. **No `methane.jl` change** — the `use_lch4` path is
byte-identical; guarded by a new `test_methane.jl` analytic-oracle assertion
(threshold gating + exact rate) under `--check-bounds=yes`.

### What is now validated end-to-end (the last CH4 source term)
Combined with §11 (production/oxidation/aerenchyma/diffusive surface flux at
`finundated>0`), **every CH4 source-regime term is now diffed against Fortran.**
Ebullition was the last one that no on-disk asset could exercise.

### Reproducing
```bash
# 1. edit a warm restart's CONC_CH4 (×5) and generate the Fortran reference:
#    julia +1.12 --project=. scripts/edit_ch4_restart.jl <raw.r.nc> finidat_ch4ebul_x5.nc 5.0
#    (already applied; run dir SYMFLUENCE_data/clm_bgc_spinup/bow_ref_ch4ebul/)
#    bash SYMFLUENCE_data/clm_bgc_spinup/bow_ref_ch4ebul/run_bow_ch4ebul.sh   # plain, ~6 steps
julia +1.12 --project=. scripts/probe_ch4_ebul_threshold.jl 5   # confirm scale crosses threshold (no box)
julia +1.12 --project=. scripts/fortran_parity_ch4_ebul.jl 1    # n0 ebullition scorecard
```
