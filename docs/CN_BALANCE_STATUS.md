# CN mass-conservation check — status

**Both the carbon and nitrogen checks are ON by default and FATAL, and both now
conserve to the machine-precision floor.** The one previously-documented nitrogen
residual (net mineralization created in the soil) is **RESOLVED** — see "Resolved" below.

| | threshold | Bow `use_cn`, worst over summer (28 steps) + autumn (480-step free run) |
|---|---|---|
| `errcb_col` | `cerror` = 1e-7 | **1.4e-11** ✅ |
| `errcb_grc` | `cerror` = 1e-7 | **1.4e-11** ✅ |
| `errnb_col` | `nerror` = 1e-3 (fatal) | **1.1e-12** ✅ passes |
| `errnb_col` | `nwarning` = 1e-7 (warn) | **1.1e-12** ✅ (was 3.1e-4) |

Harness: `scripts/cn_conservation_audit.jl [summer|autumn|both]`. N-side pool
audits: `scripts/probe_veg_pool_audit_n.jl`, `scripts/probe_n_split_vegsoil.jl`,
`scripts/probe_n_column_balance.jl`.

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

## Resolved: net mineralization was created in the SOIL, not the vegetation

The residual was **~3.05e-4 gN/m2/step of net mineralization added to the soil mineral
pool but never removed from the soil organic (decomp) pools** — the soil decomposition
cascade's N-transfer step into `decomp_npools_vr` was **entirely missing** from the port.

### How the earlier "it's in the vegetation" reading was a red herring
The old localisation compared `d(veg N)` against a `soil N` total and concluded veg gains
more than soil supplies. But splitting the column N into veg / soil-mineral / soil-decomp
(`scripts/probe_n_split_vegsoil.jl`) shows the veg side conserves exactly:

```
veg (p2c totn_patch)  d = +9.95e-4   = plant uptake 1.02e-3 − litterfall 2.5e-5   ✔
soil mineral (sminn)  d = −6.95e-4                                                  ← too small
soil decomp pools     d = +2.52e-5   = litterfall only                             ← should also LOSE ~3e-4
TOTN_COL              d = +3.25e-4   ← created
```

`scripts/probe_n_column_balance.jl` confirms the plant-uptake path is airtight under FUN:
`p2c(sminn_to_plant_fun) − vr_soil_debit = 0.0` (and `Nfix ≈ 0`, so the FUN-fixation
`nfix_to_sminn` p2c is not involved here). The mineral pool is debited exactly what the
plant takes up. The leak is that decomposition adds **gross_nmin − actual_immob**
(net mineralization) to `smin_nh4_vr` in `soilbiogeochem_n_state_update1!`
(`nitrif_denitrif.jl`) while the matching **organic-side** debit was never applied.

### The bug
`SoilBiogeochemNStateUpdate1Mod.F90:118-160` builds `decomp_npools_sourcesink` from the
cascade — donor pool loses `decomp_cascade_ntransfer_vr`, receiver gains
`ntransfer + decomp_cascade_sminn_flux_vr`, terminal transitions debit the donor by the
mineral-N flux — and that sourcesink is applied to `decomp_npools_vr` (via litter vertical
transport). The **carbon** column kernel `_csu1_col_kernel!` (`c_state_update1.jl:367-382`)
ports this exactly. The **nitrogen** column kernel `_nsu1_col_kernel!`
(`n_state_update1.jl`) did **not**: it only SET `decomp_npools_sourcesink_col = phenology
litterfall` and stopped. So the organic pools received litterfall but never lost the N
that decomposition mineralized into the smin pool → N created = net mineralization, every
step. It concentrated in the growing season because that is when decomposition (warm, wet
soil) and net mineralization peak. Carbon could not catch it: carbon has no
mineral-pool analogue — its cascade just respires to the atmosphere and was ported.

### The fix
Add the N-cascade donor/receiver/terminal block to `_nsu1_col_kernel!`, a line-for-line
mirror of the carbon kernel and of the Fortran N formulas
(`decomp_cascade_ntransfer_vr_col` + `decomp_cascade_sminn_flux_vr_col`), gated on
`!use_soil_matrixcn` and threaded through `n_state_update1!` from `cn_driver.jl`
(`cascade_donor_pool` / `cascade_receiver_pool` / `ndecomp_cascade_transitions`, the same
three the C update already receives). The new kwargs default to a no-op
(`ndecomp_cascade_transitions = 0`) so existing unit-test callers stay bit-identical.

**Result:** `errnb_col` fell from **3.1e-4 → 1.1e-12** (summer) and holds at the ~1e-12
floor across the 480-step autumn free run — four orders under `nwarning` = 1e-7, matching
the carbon side. Both checks are **ON by default and FATAL**.

### Adjacent observations (not the cause; left as-is)
* `ndep_to_sminn_col` reads 0.0 (no ndep stream fed to this harness) and FUN `Nfix ≈ 0`
  in this window — neither contributes to the residual, and both are correctly wired
  (deposition via `n_deposition!`, FUN fixation p2c'd to `nfix_to_sminn_col` in `fun.jl`).
