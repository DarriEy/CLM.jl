# Crop lifecycle — CTSM → CLM.jl dependency map (2026-07-19)

Prerequisite reading: `docs/CROP_PARITY.md` (the reference), #253 (Fortran CN-crop
run), #264 (patch index alignment).

This document maps CTSM's `CropPhenology` (`CNPhenologyMod.F90:2006-2638`) and its
callees onto CLM.jl, establishes the **true dependency order**, and records exactly
which reference field validates which step. It is written **before** implementation
deliberately: the ordering is the part that previous attempts got wrong, and it is
useful on its own.

## 0. The reference — where it is and what it can prove

`SYMFLUENCE_data/clm_bgc_spinup/crop_ref_usplains/` (note: **directly under
`SYMFLUENCE_data`, not under `installs/`** — `docs/CROP_PARITY.md` gives the wrong
path).

| Asset | Content |
|---|---|
| `refs_crop_may/bgcdump_after_fire_n38665..76.nc` | fertilization window OPEN |
| `refs_crop_july/bgcdump_after_fire_n39865..76.nc` | mid-season, fertilization window CLOSED |
| `Crop_USplains.clm2.r.2020-05-30-14400.nc` | **restart — carries every intermediate** |

The bgcdumps carry only lifecycle *outputs* (`CROPLIVE`, `CPHASE`, `HUI`,
`GDDACCUM`, `GDDTSOI`, `VF`, `FERTNITRO`, `FERT`, `SOYFIXN`, `HARVDATE`,
`NYRS_CROP_ACTIVE`, `FERT_TO_SMINN`, `SOYFIXN_TO_SMINN`). **The restart file is the
more valuable reference** and was not used by prior work: it additionally carries
`idop`, `gddmaturity`, `huileaf`, `huigrain`, `peaklai`, `fert_counter`,
`onset_counter`, `cumvd`, `hdidx` — i.e. the intermediates that let each step be
validated *in isolation* rather than only end-to-end.

### The 8 crop patches (0-based 5..12; #264 made Julia build these same 8)

| p0 | itype | idop | gddmaturity | huileaf | huigrain | fertnitro | FERT (May) |
|---|---|---|---|---|---|---|---|
| 5 | 17 tmp_corn | 136 | 1477 | 44.319 | 850.077 | 14.7488 | 9.69261e-6 |
| 6 | 18 irr_tmp_corn | 136 | 1477 | 44.303 | 849.893 | 14.7488 | 9.69261e-6 |
| 7 | 19 swheat | 127 | 1700 | 85.000 | 1020.000 | 7.8962 | 5.72700e-6 |
| 8 | 20 irr_swheat | 127 | 1700 | 85.000 | 1020.000 | 7.9868 | 5.77940e-6 |
| 9 | 23 tmp_soybean | 145 | 1259 | 37.778 | 629.631 | 0.7075 | 1.56687e-6 |
| 10 | 24 irr_tmp_soybean | 145 | 1257 | 37.721 | 628.680 | 0.7075 | 1.56687e-6 |
| 11 | 75 (corn-family) | 106 | 1335 | 40.050 | 595.314 | 14.7488 | 9.69261e-6 |
| 12 | 76 (corn-family) | 106 | 1334 | 40.028 | 595.130 | 14.7488 | 9.69261e-6 |

`cumvd`/`hdidx` are **NaN on every patch** and no winter-wheat CFT (itype 21/22) is
present. **Vernalization is NOT exercised by this reference** — it can be wired but
NOT validated here. Recorded as a gap, not papered over.

`SOYFIXN = 0` in *both* windows (CTSM's `CNSoyfix` needs a later phase than either
samples), so **`n_soyfix!` remains unvalidatable** — the same standard #253 applied
to `n_fert!`.

### Three formulas the reference pins EXACTLY (checked, not assumed)

These were back-solved from the restart + dumps and agree to 5-6 digits on all
8 patches, so they are genuine acceptance targets rather than "is it wired" tests:

1. **`huileaf = lfemerg(ivt) * gddmaturity`** — implied `lfemerg` = 0.03 (corn,
   soybean, 75/76) and 0.05 (spring wheat). Exact round numbers.
2. **`huigrain`** — the corn-family branch (`crmcorn` regression) vs the plain
   branch is *distinguishable from the data*: itypes **17, 18, 75, 76 take the
   corn-family branch** (implied `grnfill` 0.65, 0.65, 0.50, 0.50) and **19, 20,
   23, 24 take the plain branch** (0.60, 0.60, 0.50, 0.50). Taking the wrong branch
   changes `huigrain` by 7-11% and is immediately visible.
3. **`fert = (manunitro(ivt)*1000 + fertnitro) / fert_counter`** with
   `fert_counter = ndays_on * secspday = 20 * 86400`. Back-solving `manunitro*1000`
   from the May dump gives **2.000000 on all 8 patches**. This is the single most
   precise target available and it is what makes `n_fert!` non-vacuous.

## 1. Dependency chain (the actual ordering)

CTSM's `CropPhenology` is one loop, but the data dependencies inside it are strictly
ordered. Wiring out of order is what produces green-but-dead paths.

```
                    [ crop_phenology_init! ]
                     inhemi, minplantjday/maxplantjday   <-- from mn/mxNH/SHplantdate
                              |
                              v
       [ accumulators — ALREADY LIVE, verified #218 ]
        gdd020/gdd820/gdd1020 (20yr means)  t_a10, t_a5min, t_a10min
        hui, gddaccum, gddtsoi  (reset to 0 when !croplive)
                              |
                              v
  (1) SOWING WINDOW      get_swindow -> is_in_sowing_window, is_end_sowing_window
                              |
                              v
  (2) PLANTING DECISION  do_plant_normal / do_plant_lastchance / do_plant_prescribed
        winter cereal:  a5tmin <= minplanttemp  AND gdd020 >= gddmin
        everything else: t10 > planttemp AND a10tmin > minplanttemp AND gdd820 >= gddmin
                              |
                              v
  (3) plant_crop!  <-- THE ONLY WRITER OF croplive=true
        croplive=T, idop=jday, iyop=kyr, harvdate=NOT_Harvested, sowing_count++,
        leafc_xfer=seed, crop_seedc_to_leaf, aleaf/astem/aroot,
        gddmaturity  (per-crop-family from gdd20 + hybgdd)
                              |
                              v
  (4) PHASE THRESHOLDS   huileaf = lfemerg*gddmaturity
                         huigrain = grnfill*gddmaturity  (corn family: crmcorn adj)
                              |
                              v
  (5) if croplive:  cphase = planted
        (5a) vernalization!   [winter cereal only; modifies vf -> scales hui/gddtsoi]
        (5b) idpp = DaysPastPlanting(idop, jday)
        (5c) onset_counter -= dt
        (5d) HARVEST DECISION  do_harvest = hui >= gddmaturity  OR  idpp >= mxmat
                              |
              +---------------+---------------+
              |               |               |
              v               v               v
      leafemerge         do_harvest       grainfill
      (gddtsoi>=huileaf  croplive=F       (hui>=huigrain)
       AND hui<huigrain) cphase=harvest   cphase=grainfill
       cphase=leafemerge harvdate=jday    bglfr=1/(leaf_long*...)
       onset_flag=1      offset_flag=1
       onset_counter=dt
   >>> fert_counter = ndays_on*secspday
   >>> fert = (manunitro*1000+fertnitro)/fert_counter   <-- (6) FERTILIZER
              |
              v
  (7) every step while live:  if fert_counter <= 0 then fert = 0
                              else fert_counter -= dt
                              |
                              v
  (8) n_fert!   p2c scatter of fert_patch -> fert_to_sminn_col   [NOW non-vacuous]
      n_soyfix! needs cphase+croplive (not fert_patch) -- but has NO reference
```

**The load-bearing fact:** `fert` is assigned **once**, at the single timestep of
the phase-1→phase-2 onset transition (step 5, `onset_flag` branch), and then *held*
every subsequent step until `fert_counter` decrements to 0 — it is not recomputed.
So `fert` is non-zero for exactly `ndays_on = 20` days after leaf emergence. That is
why the May window shows `FERT` = 9.69e-6 and the July window shows `FERT ≡ 0`:
**the window closed, the path is not dead.** A mid-season dump alone would have read
as a dead path — this is the trap `docs/CROP_PARITY.md` flags, and it is structural,
not incidental.

## 2. Julia-side gap census (verified against source, not comments)

The **data structures are complete** — every state field, param and constant the
port needs already exists. Nothing needs to be added to the types. The gaps are
wiring and logic:

| # | Gap | Location | Severity |
|---|---|---|---|
| G1 | `crop_phenology!` body is a SIMPLIFIED stub — sets only `cphase` and grain-fill `bglfr`. No planting decision, no `croplive` transition, no harvest, no fertilizer, no onset. | `src/biogeochem/phenology.jl:2344` | **the port gap** |
| G2 | `plant_crop!` — the only writer of `croplive=true` — has ZERO callers | `phenology.jl:2394` | **dead** |
| G3 | `vernalization!` has ZERO callers | `phenology.jl:2438` | dead (unvalidatable here) |
| G4 | `jday`/`kyr`/`dayspyr`/`use_fertilizer` are kwargs of `crop_phenology!` that **nothing ever supplies** — they silently take defaults `1/1/365.0/false`. `cn_phenology!` has no such parameters at all, so there is no plumbing path from the driver. | `phenology.jl:2358`, `cn_driver.jl:791` | **silent-default** |
| G5 | `PftConPhenology` is constructed **without** `mnNHplantdate`/`mxNHplantdate`/`mnSHplantdate`/`mxSHplantdate` — so those four arrive **empty** | `cn_driver.jl:772-787` | **dead** |
| G6 | `cn_phenology_init!` is called without pftcon/patch/gridcell, so `crop_phenology_init!` never fires → `inhemi`, `minplantjday`, `maxplantjday` stay **empty** | `cn_driver.jl:771` | **dead** |
| G7 | `n_fert!`/`n_soyfix!` unwired (correctly — blocked on G1 producing a real `fert_patch`) | `cn_driver.jl:402-407` | blocked by G1 |

**G5 and G6 are the deepest and were not previously identified.** They are the
`dead-initcold-systemic` bug class again: even if `crop_phenology!` were fully
ported, `minplantjday`/`maxplantjday` would be **empty vectors**, so the sowing
window could not be computed at all and planting could never fire. G4 compounds it —
`jday` defaulting to `1` means the planting-window test would be evaluated on
January 1 every single step. **Any attempt to wire the lifecycle without fixing G4-G6
first would have produced a crop that never plants, and the failure would have looked
like a physics bug.**

## 3. Implementation order (each step validated before the next)

| Step | What | Validated against |
|---|---|---|
| **S0** | Plumbing: fix G4, G5, G6. Thread `jday`/`kyr`/`dayspyr`/`use_fertilizer`; add the 4 plantdate params; fire `crop_phenology_init!`. | `minplantjday`/`maxplantjday` non-empty and matching `get_calday(mn/mxNHplantdate)`; assert non-no-op |
| **S1** | Planting decision + `plant_crop!` + `gddmaturity` | restart `idop` (136/136/127/127/145/145/106/106), `gddmaturity` (1477…1334), `croplive` |
| **S2** | Phase thresholds `huileaf`/`huigrain` incl. the corn-family `crmcorn` branch | restart `huileaf`, `huigrain` — **exact**, branch-distinguishing |
| **S3** | Phase transitions `cphase` + onset/offset flags | dump `CPHASE` = 2 (May) / 2,2,3,3,2,2,3 (July) |
| **S4** | Fertilizer window: `fert_counter`, `fert` | dump `FERT` (May, 8 values) — **exact to 6 digits**; `FERT ≡ 0` in July |
| **S5** | `n_fert!` wiring | dump `FERT_TO_SMINN` (8 non-zero cols, max 9.693e-6) |
| **S6** | Harvest: `croplive=false`, `harvdate`, `offset_flag` | `HUI_PERHARV` 1425-1481, one harvest/yr |
| — | `vernalization!` | **NO REFERENCE** — wire behind `croplive`, mark unvalidated |
| — | `n_soyfix!` | **NO REFERENCE** (`SOYFIXN ≡ 0` both windows) — do not wire |

Everything gated on `use_crop` so the default path stays byte-identical.

## 3b. WIRED AND VALIDATED (2026-07-19)

S0-S6 are implemented and validated; `test/test_crop_lifecycle.jl` (106 assertions)
carries the reference numbers. Tolerances are at **round-off** (`rtol` 1e-12/1e-13),
not a loose band — the full-precision reference values were taken from the files
rather than the rounded ones printed in `docs/CROP_PARITY.md`, which is what makes a
tight tolerance possible.

| Step | Field | Julia | Reference | Agreement |
|---|---|---|---|---|
| S1 | `croplive` | `true` ×8 | `CROPLIVE` = 1 ×8 | flips (was structurally impossible) |
| S1 | `idop` / `iyop` | 136 / 2020 | — | set (date itself not under test, see below) |
| S1 | `gddmaturity` | 1477.2873853187507 … 1334.2627384206983 | identical | `rtol` 1e-12 |
| S2 | `huileaf` | 44.318621559562516 … 40.02788215262095 | identical | `rtol` 1e-12 |
| S2 | `huigrain` | 850.0767538833009 … 595.1300578305277 | identical | `rtol` 1e-12 |
| S3 | `cphase` | 2 (leafemerge) | `CPHASE` = 2 (May) | exact |
| S4 | `fert` | 9.69260860321513e-6, 5.726995403717921e-6, 1.5668663836927799e-6, … | identical | **relative error 0.0** |
| S4b | `fert` after window | 0.0, `croplive` still true | July `FERT` ≡ 0 | exact |
| S5 | `n_fert!` | wired | `FERT_TO_SMINN` 8 non-zero cols | scatter now non-vacuous |
| S6 | harvest | `croplive`→false, `harvdate`, `harvest_reason`=MATURE | one harvest/yr | exact |

The `huigrain` branch is asserted **discriminatively**: the test requires that itypes
17/18/75/76 do NOT equal `grnfill*gddmaturity` and that 19/20/23/24 DO. Taking the
wrong branch is a 7-18% error and fails.

Confirmation the parameters are real, not fitted: `clm5_params.nc` (pft=79) gives
`lfemerg` = 0.03/0.05, `grnfill` = 0.65/0.60/0.50, `manunitro` = 0.002 — exactly the
values back-solved from the reference before the paramfile was opened. itypes 75/76
are `tropical_corn` / `irrigated_tropical_corn`, which is why they take the
corn-family branch.

### Not validated — stated so a green suite is not over-read

- **Planting DATE selection.** The reference `idop` (136/127/145/106) depends on the
  run's 20-year GDD climatology (`gdd020/gdd820/gdd1020` are 20-yr running means),
  which is not reconstructible from the dumps. The planting *conditions* are
  synthesized in the test; the date is an input, not an output under test. The
  decision logic is ported but its date output is unverified.
- **`vernalization!`** — wired (it had no caller, which is what kept it dead) but
  **UNVALIDATED**: the reference has no winter-wheat CFT and its `cumvd`/`hdidx` are
  NaN on every patch.
- **`n_soyfix!`** — still NOT wired, deliberately. `SOYFIXN` ≡ 0 in *both* reference
  windows, so wiring it would be precisely the blind wiring #218/#253 refused.

### Remaining work, in order

1. A third reference window (or phase-targeted bracket) that reaches the growth phase
   `CNSoyfix` requires → then wire and validate `n_soyfix!`.
2. A winter-wheat CFT reference → validate `vernalization!`.
3. An end-to-end `scripts/fortran_parity_crop.jl` single-step diff on the
   `fortran_parity_common.jl` pattern, which would also close the planting-date gap
   by running from the reference restart (verify the harness datm year/hour first —
   the #233/#240/#243 lesson).
4. Crop allocation (`arepr`/grain pools) and `CNCropHarvestToProductPools` value
   parity — reachable now that the lifecycle drives `cphase`.

## 3c. `CNSoyfix` activation criterion (2026-07-19, read from source)

Source: `CNNDynamicsMod.F90:319-463`, called from `CNDriverMod.F90:325`.

**Correction to a standing assumption.** `docs/CROP_LIFECYCLE_MAP.md` §3 and
`docs/CROP_PARITY.md` both say `n_soyfix!` "needs `cphase`" and that
`CNSoyfix` "fires only in a later growth phase". **Neither is true.**
`CNSoyfix` never reads `cphase`. It is not phase-gated. The real gate is
`fpg`, and that changes what a valid reference window is.

The routine writes `soyfixn_patch` on **every** patch every step (it has an
explicit `else → soyfixn = 0` on both branches), then `p2c`-scatters to
`soyfixn_to_sminn_col`. So it is never "skipped" — a zero is a computed zero.

### The three gates, in order

```
(A) croplive(p) == .true.
(B) patch%itype(p) ∈ {ntmp_soybean, nirrig_tmp_soybean,
                      ntrp_soybean, nirrig_trp_soybean}
(C) fpg(c) < 1.0            <-- THE LOAD-BEARING GATE
```

If (A) or (B) fails → `soyfixn = 0`. If (C) fails → `soyfixn = 0` explicitly
("if nitrogen demand met, no fixation"). Only when all three hold is anything
computed.

In the 8-patch collapsed set, (B) selects **exactly patches p0 = 9 (itype 23,
`temperate_soybean`) and p0 = 10 (itype 24, `irrigated_temperate_soybean`)**.
`ntrp_soybean`/`nirrig_trp_soybean` are the *tropical* soybeans (not 75/76 —
those are `tropical_corn`/`irrigated_tropical_corn`, per §3b) and are absent
from this surfdata. **So at most 2 of 8 patches can ever be non-zero here**,
and a test asserting non-zero on all 8 would be wrong.

### The computation when all three gates hold

```
soy_ndemand = plant_ndemand(p) - plant_ndemand(p)*fpg(c)      ! = pnd*(1-fpg)

fxw = wf(c) / 0.85                                            ! wf = soil water
                                                              !   frac of whc, top 0.5 m
fxn = 0                       if sminn(c) >  30
    = 1.5 - 0.005*(sminn*10)  if 10 < sminn(c) <= 30
    = 1                       if sminn(c) <= 10

GDDfrac = hui(p) / gddmaturity(p)
fxg = 0                    if GDDfrac <= 0.15
    = 6.67*GDDfrac - 1     if 0.15 < GDDfrac <= 0.30
    = 1                    if 0.30 < GDDfrac <= 0.55
    = 3.75 - 5*GDDfrac     if 0.55 < GDDfrac <= 0.75
    = 0                    if GDDfrac >  0.75

fxr = max(0, min(1, fxw, fxn) * fxg)
soyfixn(p) = min(fxr * soy_ndemand, soy_ndemand)
```

Note `min(1, fxw, fxn)` is a 3-way min including the literal `1`, so the water
and N factors are capped at 1 but `fxg` is applied *after* the cap and is not
itself capped by it. The final `min(soyfixn, soy_ndemand)` is redundant given
`fxr <= 1` but is kept for fidelity.

### Where to bracket a dump — and the trap

The **GDD window** is `0.15 < hui/gddmaturity < 0.75`. With soybean
`gddmaturity` ≈ 1259 / 1257 (p0 9/10, from the restart), that is
`hui ∈ (189, 944)`, peaking at `fxg = 1` for `hui ∈ (378, 692)`. A dump
bracketing the `fxg = 1` plateau is the most discriminative target, because
there `soyfixn = min(1, fxw, fxn) * plant_ndemand * (1 - fpg)` exactly.

**Why the two existing windows are zero is NOT the same reason in each:**

- **May (n38665, 2020-05-30):** soybean `idop` = 145, so at jday 151 `hui` is
  barely accumulating → `GDDfrac <= 0.15` → `fxg = 0`. A *timing* zero. A
  later window fixes this.
- **July (n39865, 2020-07-19):** `GDDfrac` is plausibly inside the window here,
  so if this dump is also zero the cause is likely **gate (C)**: `fpg = 1`,
  i.e. the soil supplies all the N the plants demand and there is no deficit to
  fix. That is a *state* zero, and **no choice of window fixes it** — it is a
  property of the run's N budget (a fertilized, spun-up soil column), not of
  when you look.

**This distinction is the whole risk of the task.** Before generating any new
window, `fpg` and `GDDfrac` must be read out of the July reference to determine
which of the two zeros is in play. If `fpg == 1` at July, a new dump window is
the wrong instrument and the blocker is an N-budget one.

## 3d. Why `SOYFIXN` ≡ 0 — ROOT CAUSE FOUND: `use_fun = .true.`

**The existing reference cannot produce a non-zero `SOYFIXN` at any window, and
the reason is not timing.** `CNSoyfix` is never called in it.

`CNDriverMod.F90:319-329`:

```fortran
if (use_crop) then
   call CNNFert(...)
   if (.not. use_fun) then  ! if FUN is active, then soy fixation handled by FUN
      call CNSoyfix (...)
   end if
end if
```

and `crop_ref_usplains/lnd_in` has:

```
 use_cn   = .true.
 use_crop = .true.
 use_fun  = .true.        <-- CNSoyfix is compiled in but NEVER CALLED
 use_nitrif_denitrif = .true.
```

With FUN active, soybean fixation is handled inside `CNFUNMod` and
`soyfixn_patch` is simply never written — it holds its initialized value, 0.
`SOYFIXN` in the dumps is therefore **not a computed zero at all**; it is an
untouched field. The `p2c` scatter to `SOYFIXN_TO_SMINN` likewise never runs.

### The evidence that rules out the timing hypothesis

The natural assumption (recorded in prior docs) was that the windows sample the
wrong growth phase. **The July window disproves that**: at `n39865` the two
soybean patches sit at the *most favourable possible* point of the criterion.

| Quantity | p0 = 9 (itype 23) | p0 = 10 (itype 24) | Gate |
|---|---|---|---|
| `CROPLIVE` | 1 | 1 | (A) PASS |
| itype | 23 `temperate_soybean` | 24 `irrigated_temperate_soybean` | (B) PASS |
| `HUI` | 552.196 | 549.678 | — |
| `gddmaturity` | 1259.261 | 1257.360 | — |
| `GDDfrac` | **0.43851** | **0.43717** | in (0.30, 0.55] → **`fxg` = 1**, the plateau |
| `SOYFIXN` | **0.0** | **0.0** | — |

Every gate that a *window* controls is satisfied, and the growth-stage factor is
at its maximum. A zero here cannot be explained by phase. (For completeness the
May window *is* a genuine timing zero: `HUI` = 64.93 → `GDDfrac` = 0.0516 →
`fxg` = 0. The two windows are zero for two different reasons, and neither is
fixable by choosing a third window.)

Corroborating: the July restart has `fpg` = **NaN** on all soil columns and
`plant_ndemand` = **0.0 on all 24 patches** — both are FUN-path artifacts. Under
`use_fun = .true.` the `fpg`/`plant_ndemand` pair that `CNSoyfix` consumes is not
maintained by the non-FUN nutrient-competition path at all, so even a forced call
would have read NaN/zero inputs.

### Consequence — `n_soyfix!` stays UNWIRED

`SOYFIXN ≡ 0` is **not** evidence about a window. It is evidence that this
reference is configured for a code path that excludes `CNSoyfix` entirely.
Wiring `n_soyfix!` against it would scatter zeros and pass any "is it wired?"
test while moving nothing — the exact failure #218, #253 and #266 each refused.
**The refusal stands.**

### What would actually unblock it (precise, in order)

1. **A reference run with `use_fun = .false.`** — this is the whole blocker. It
   is a *runtime* namelist flag (no rebuild; same binary), so the existing case
   from `scripts/setup_crop_ref_case.py` can be reconfigured. Note this is not a
   cosmetic switch: it swaps FUN for `NutrientCompetitionFlexibleCNMod`, which is
   what computes the `fpg`/`plant_ndemand` pair `CNSoyfix` needs
   (`NutrientCompetitionFlexibleCNMod.F90:442`).
2. **Restart it from `Crop_USplains.clm2.r.2020-07-19-14400.nc`** and dump the
   same `n39865..76` bracket. That restart already sits on the `fxg = 1` plateau
   (table above), so the growth-stage gate is pre-satisfied and only a handful of
   steps are needed — a warm restart, not a re-spin-up. Whether a FUN restart
   loads cleanly into a non-FUN run is the open risk and must be verified, not
   assumed.
3. **Extend the BGCDUMP instrumentation to carry `CNSoyfix`'s inputs**, which the
   current dump does not: `fpg_col`, `plant_ndemand_patch`, `sminn_col`, `wf_col`
   (plus the existing `HUI`/`gddmaturity`/`CROPLIVE`). Without them the only
   possible test is end-to-end; with them each of `fxw`/`fxn`/`fxg`/`fxr` can be
   validated in isolation, which is what makes the branch structure
   discriminative rather than a single number that many wrong formulas reproduce.
4. Only then wire `n_soyfix!` and assert values.

Until (1) exists, **`n_soyfix!` is blocked on a reference-configuration gap, not
on Julia code.** The port side is ready; the ground truth is not.

## 4. Standing rules for this work

- Do **not** wire a call whose upstream input is still zero (#253's refusal of
  `n_fert!` was correct and the standard is retained).
- Assert the **value** against the reference, not that the call happened
  (`conservation-is-not-accuracy`, `vacuous-checks-bug-class`).
- Assert the **wiring AND a non-no-op body** (`dead-initcold-systemic`).
- Pick the window where the process is **active** — `FERT` is only alive in May.
