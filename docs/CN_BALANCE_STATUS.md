# CN mass-conservation check — status

**The carbon check is ON by default and FATAL.** The nitrogen check is ON and fatal at
Fortran's `nerror`, which it passes, but it has one documented open residual (below).

| | threshold | Bow `use_cn`, worst over summer (28 steps) + autumn (480-step free run) |
|---|---|---|
| `errcb_col` | `cerror` = 1e-7 | **1.4e-11** ✅ |
| `errcb_grc` | `cerror` = 1e-7 | **1.4e-11** ✅ |
| `errnb_col` | `nerror` = 1e-3 (fatal) | **3.1e-4** ✅ passes |
| `errnb_col` | `nwarning` = 1e-7 (warn) | **3.1e-4** ❌ **exceeds — see "Open" below** |

Harness: `scripts/cn_conservation_audit.jl [summer|autumn|both]`.

---

## Why the carbon check "could not run" — and what it found when it could

The C half was default-OFF, and the reason given was a supposed wood/crop product-pool
gap. That reason was wrong. The real reason was that **every term the check reads was
structurally dead**, so it could not disagree with anything:

* `gpp_col`, `er_col`, `ar_col`, `rr_col`, `npp_col`, `nep_col`, `fire_closs_col`,
  `hrv_xsmrpool_to_atm_col`, `gru_conv_cflux_col`, `nbp_grc`, `nee_grc`,
  `landuseflux_grc` — the **entire column/gridcell half of
  `CNVegCarbonFluxType::Summary` was a one-line stub** ("Column-level p2c aggregation is
  stubbed pending subgridAveMod port") even though `subgridAveMod` *is* ported.
* `hr_col`, `somhr_col`, `lithr_col`, `cwdhr_col`, `som_c_leached_col`,
  `denit_col`, `f_n2o_nit_col`, `smin_no3_leached_col`, `som_n_leached_col` — the two
  **SoilBiogeochem flux summaries were ported but never called** from
  `CNDriverSummarizeFluxes`.
* `totc_col` / `totn_col` — the begin/end masses — **contained no vegetation carbon or
  nitrogen at all**. Fortran builds `totc_col` from `totc_p2c_col` (the p2c of
  `totc_patch`, which *includes* cpool/xsmrpool/ctrunc) and `totecosysc_col` from
  `totvegc_col` (which *excludes* cpool). The port passed the cpool-excluding aggregate
  for both and never supplied `totc_p2c_col`, so it defaulted to `zeros()`.

Verified against the live model before the fix: `ar_patch` = 1.5e-7 while `ar_col` = 0.0.

Once wired, the check ran — **and failed, by 4.5e-2 gC/m2/step against a 1e-7
threshold.** Chasing that down produced four independent, real carbon-destroying or
carbon-creating bugs:

### 1. `lf_f` / `fr_f` litter fractions were never populated → all leaf + fine-root litter destroyed
`lf_flab/lf_fcel/lf_flig` (and the fine-root trio) are read from the parameter file as
vectors, but the 2-D `lf_f[pft, i]` / `fr_f[pft, i]` arrays that every litter-routing
kernel indexes were allocated as `zeros(npft, ndecomp_pools)` and **never filled**
(Fortran packs them in `pftconMod.F90:822-834`). Every leaf and fine-root litter flux is
formed as `flux * lf_f[ivt,i] * ...`, so **the entire leaf + fine-root litterfall — from
phenology, from gap mortality, and from fire; C and N alike — was multiplied by zero and
destroyed.** It left the vegetation pools and never arrived in the soil. Only the
woody→CWD path (which does not use `lf_f`) survived, which is why the soil still received
~38% of the litter flux and the leak was not obvious.
Fixed: `_pftcon_derive_litter_fractions!` (`src/infrastructure/read_params.jl`).

### 2. Phenological litterfall was routed into ONE soil layer → 86% destroyed
`cn_litter_to_column!` was called from `cn_phenology!` **with no `nlevdecomp` argument**,
so it silently took the kwarg default `nlevdecomp = 1`. The kernel loops
`for j in 1:nlevdecomp` and deposits `flux * lf_f[ivt,i] * wtcol[p] * leaf_prof[p,j]` into
layer `j`. The profiles are normalised over the **full** decomposition column
(`sum_j prof*dzsoi_decomp == 1`), so restricting the loop to `j = 1` delivered only
`prof[p,1]*dz[1]` ≈ **14%** — and the other ~86% of every phenological leaf and fine-root
litterfall flux (C and N) was destroyed.
Fixed: `nlevdecomp`/`i_litr_*`/`npcropmin` threaded through `cn_phenology!`.

### 3. `totc_col` carried no vegetation carbon (see above)
Symptom: a **constant** ~391.7 gC/m2 offset in every step's `errcb` at Bow (= `totvegc`),
because the BEGIN summary (`cn_driver_summarize_states!`) and the END summary
(`clm_driver.jl`) disagreed on the basis. Fixed in both.

### 4. A PFT off-by-one made every TREE non-woody in the C/N state update
`c_state_update1.jl` and `n_state_update1.jl` computed
```julia
is_woody = woody[ivt[p]] == 1     # WRONG
```
`ivt[p]` is the **raw 0-based Fortran PFT index** (bare ground == 0 — exactly what the
`ivt[p] >= 1` guard on the line above skips). pftcon arrays are 1-based, so the PFT's row
is `ivt[p] + 1` — which is what `allocation.jl` (`ivt = itype[p] + 1`), the gap-mortality
kernel and `c_iso_flux.jl` all do. Indexing with the raw value read the **previous PFT's**
woody flag: for patch type 1 (needleleaf evergreen temperate **tree**) it read `woody[1]`
— bare ground — and got 0.0.

**Every tree in the model was treated as non-woody by the C and N state updates.** The
whole `if is_woody` branch was skipped, so cpool was never debited for livestem/livecroot
maintenance respiration, for allocation to livestemc/deadstemc/livecrootc/deadcrootc
(+ storage), or for woody growth respiration — while `allocation.jl` (which indexes
*correctly*) still put all of those fluxes into `mr`, `gr` and hence `ar`. Carbon was
counted as respired but never removed from any pool: **~4.4e-4 gC/m2/step created from
nothing in the Bow tree patch, every step.** The grass patch closed to 1.9e-14 precisely
*because* grass is genuinely non-woody, so the wrong flag happened to be harmless there —
which is what made the bug pinpointable.

### 5. Gridcell check disabled by a `NaN * 0` special column
Special (non-soil) columns have `wtgcell == 0` and were left at the NaN allocation default
for `totc_col`/`totn_col`. The gridcell c2g then computed `NaN * 0.0 == NaN`, and
`abs(NaN) > cerror` is **false** — so the gridcell half of the check could never fail no
matter how badly it was violated. Fortran zeroes special columns
(`SoilBiogeochemCarbonStateType.F90:686`); so do we now.

### Also fixed (same class, found en route)
* `lag_npp_col` was blanket-zeroed every step by `_zero_cnveg_flux_arrays!`. Fortran's
  `SetValues` does **not** zero it — it is spval-init, a restart variable, and written
  only by the Summary relaxation. Zeroing it destroyed the exponential memory every step,
  collapsing the 10-day relaxation that CNNFixation reads to a step-local fraction.

**Result:** carbon now conserves at **|errcb| ≤ 1.4e-11** (column and gridcell) across the
summer window, a 480-step autumn free run through the leaf-offset ramp, and cold start —
four orders of margin under `cerror` = 1e-7. The check is **ON by default and FATAL**.

---

## Open: nitrogen is created in the vegetation pools (~3.3e-4 gN/m2/step)

The N budget clears Fortran's `nerror` (1e-3) but not its `nwarning` (1e-7).
**Neither threshold has been retuned and no term has been re-zeroed.** The warning is
capped with `maxlog` only so that one documented, non-fatal residual does not emit an
identical line on all 480 steps of a free run.

Localised (`scripts/probe_veg_pool_audit.jl`, summer, Bow column 1, one step):

```
d(totn_col) = +3.25e-4      <- the N pool GROWS
d(soil N)   = -6.70e-4      <- soil supplies 6.70e-4
d(veg N)    = +9.95e-4      <- vegetation gains 9.95e-4
                                => +3.25e-4 gN/m2/step created from nothing, in VEG
```

So nitrogen is created inside the **vegetation** N pools, not lost in the veg→soil
transfer (that transfer conserves exactly now — it shares the `lf_f` and `nlevdecomp`
fixes with carbon, and the N side of `n_state_update1` shares the `is_woody` fix). It is
concentrated in the **growing season**: in the autumn window most steps sit at ~1e-13 and
the worst is 1.07e-5, versus 3.1e-4 in summer — consistent with an error in the active
N-uptake / allocation / retranslocation path rather than in turnover.

Two adjacent dead inputs were also observed and are *not* the cause (they would make the
pool grow *less*, not more), but should be chased with it:

* `ndep_to_sminn_col` reads **0.0** — atmospheric N deposition reaches the mineral pool
  nowhere.
* `nfix_to_sminn_col` reads **0.0** — biological N fixation likewise.

**Next step:** run the N analogue of `scripts/probe_veg_pool_audit.jl` — snapshot every
patch-level N pool (`npool`, `retransn`, leaf/froot/livestem/... N and their
storage/xfer) across one step and compare each pool's delta against the fluxes that are
supposed to move it. Carbon's culprit fell out of exactly that audit the moment the grass
patch was seen to close to machine precision while the tree patch did not.
