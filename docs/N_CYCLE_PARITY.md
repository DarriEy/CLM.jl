# The nitrogen input/loss cycle — wiring + Fortran parity

**Status: nitrogen now enters and leaves the ecosystem.** Before this work it did
neither. This document says exactly what was wired, what it was validated
against, and what is deliberately still unwired.

Nothing was tuned to match Fortran.

---

## 1. The gap

`n_deposition!`, `n_fixation!`, `n_free_living_fixation!`, `n_leaching!`,
`n_fert!`, `n_soyfix!`, `cn_annual_update!`, `cn_products_update!` and the CN
balance checks were all **dead ports**: fully written, unit-tested, and never
called. Every "call site" in `cn_driver.jl` / `vegetation_facade.jl` was a
comment. Confirmed by grep for live call sites (a call in a comment does not
count; a unit test calling a routine in isolation does not make it alive).

The consequence was not subtle:

* `atm2lnd.forc_ndep_grc` was initialised to `0.0` and **never assigned
  anywhere**, so `ndep_to_sminn ≡ 0`.
* `cnveg_carbonflux.lag_npp_col` was initialised to `SPVAL` and **never
  updated or read from restart**, so even if `n_fixation!` had been called it
  would have taken the SPVAL branch and returned exactly zero.
* `cnveg_carbonflux.npp_col` was **never computed** (the patch→column `p2c` was
  missing), so the input to `lag_npp` did not exist either.
* `n_leaching!` was never called, so the mineral-N loss term was identically
  zero.

So mineral N had a **budget with no inputs and no losses**. Per step that is
small — which is why it hid under the ~2e-3 CN residual — but over a spinup it
sets the equilibrium C:N and hence the whole soil-carbon state.

The CN mass-conservation check, which exists precisely to catch this, was
**itself dead** — the same pattern.

---

## 2. The ndep stream problem, and how it was solved

CTSM reads N deposition from an `fndep_*` NetCDF stream (`&ndepdyn_nml`). The
CMIP6 file had been deleted from the local inputdata tree and could not be
re-fetched (the server's certificate does not cover the inputdata hostname;
disabling TLS verification was declined). PR #216's Fortran reference therefore
ran with a **synthetic zero-ndep stream**, which conveniently matched Julia's
(broken) zero — so it could not have detected this.

**Resolution: option (b).** A synthetic **non-zero** stream was created and
**both** codes were driven with the same file:

```
fndep_clm_UNIFORM1e-8_0.9x1.25_yr2000_CLMjl-parity.nc
NDEP_month = 1.0e-8 g(N)/m2/s  (= 0.3154 gN/m2/yr), uniform in space, constant in time
```

`0.32 gN/m2/yr` is in the observed range for the Canadian Rockies (Banff is a
relatively remote site, ~0.2–0.5 gN/m2/yr).

**Why uniform + constant is the *right* choice, not a cop-out.** CTSM regrids the
stream with ESMF (`ndepmapalgo='bilinear'`) and time-interpolates it with
`shr_strdata` (`tintalgo='linear'`, `taxmode='cycle'`). Reproducing ESMF's
regridding bit-for-bit would require the ESMF weight file. But for a field that
is uniform in space and constant in time, **bilinear weights sum to 1 and the
linear time weights sum to 1**, so both codes return the stream value *exactly*.
Fortran's `forc_ndep_grc` is therefore known analytically to be `1.0e-8`, and the
parity test measures the **N-cycle physics wiring** rather than forcing-file
interpolation fidelity — which is what we actually needed to test.

> The stream is **synthetic**. It is not a real N-deposition dataset, and this
> document does not claim it is. A spatially-varying real (CMIP6) file will be
> handled physically correctly by `ndep_stream.jl` (bilinear + linear-in-time,
> the same scheme), but a *bit-for-bit* match with CTSM on such a file is **not**
> claimed — see the fidelity note at the top of `src/infrastructure/ndep_stream.jl`.

---

## 3. The Fortran reference

Generated with the restored `restFile_write_dump` CN instrumentation (#216),
plus a **new dump boundary** added for this work:

| boundary | when |
|---|---|
| `before_step` | injected as the shared IC |
| `after_hydrologydrainage` | after `EcosystemDynamicsPreDrainage` — ndep+fixation in, decomposition applied, **before** leaching |
| `after_ecosysdyn_postdrain` | **NEW** — after `EcosystemDynamicsPostDrainage`, i.e. **after** leaching |

The new boundary matters: leaching happens in `PostDrainage`, so its effect was
invisible in every previously existing dump. Because Fortran now dumps **both
sides** of the leaching call, the leaching flux itself has an *exact per-level
Fortran ground truth*:

```
(leached + runoff)[j] == (no3_pre[j] - no3_post[j]) / dt
```

Two windows at Bow-at-Banff (model year 2202, `dtime = 3600 s`), both run with
the synthetic non-zero stream:

| window | nsteps | dates | why |
|---|---|---|---|
| autumn | 1760377–1760560 (184) | 2202-10-29 → | dormant |
| **summer** | **1757873–1758060 (188)** | **2202-07-16 →** | **active drainage — the only window that exercises leaching** |

Reference config (Bow `lnd_in`): `use_cn`, **`use_fun = .true.`**,
`use_nitrif_denitrif = .true.`, `use_crop = .false.`.
So Fortran runs `CNNDeposition` + **`CNFreeLivingFixation` (ffix)**, *not*
`CNNFixation`, and puts `ndep + ffix` into `smin_nh4` (the FUN branch of
`SoilBiogeochemNStateUpdate1`).

Harness: `scripts/fortran_parity_ncycle.jl [summer|autumn] [nprobe]`.

---

## 4. Scorecard

### 4a. N INPUT fluxes — **EXACT**

| flux | Julia vs Fortran |
|---|---|
| `ndep_to_sminn` | **EXACT** (rel ≤ 1e-12) at every probe, in both windows |
| `ffix_to_sminn` | **EXACT** (rel ≤ 1e-12) at every probe, in both windows |

`ffix` is checked against Fortran's own formula evaluated on the `AnnET` the dump
carries; Julia computes its own from the injected `AnnET`, so this is a real
check of the formula, not a tautology. Typical magnitudes at Bow:

* `ndep` = 1.0e-8 gN/m²/s → **3.6e-5 gN/m² per step**
* `ffix` ≈ 2.8e-9 gN/m²/s → **1.0e-5 gN/m² per step** (≈ 0.087 gN/m²/yr)
* total N input ≈ **0.40 gN/m²/yr**, previously **0.00**

### 4b. Mineral-N pools (vs `after_ecosysdyn_postdrain`), summer window

| field | worst \|rel\| over probes |
|---|---|
| `smin_nh4_vr` | **2.0e-05** |
| `smin_no3_vr` | **2.5e-07** |
| `sminn_vr`    | **6.5e-07** |

(Relative to the *profile's* magnitude — see the note in the harness on why
per-element scaling is meaningless for the deep levels, where Fortran holds
exact zeros.)

### 4c. Direction of the change — **TOWARD Fortran**

This is the number that matters: wiring the N cycle changes `use_cn` output, and
it must move Julia *toward* the reference.

| nstep | `sminn_vr` err ON | `sminn_vr` err OFF (pre-PR) | `smin_nh4_vr` ON | OFF |
|---|---|---|---|---|
| 1757875 | 6.52e-07 | 4.65e-06 | 2.01e-05 | 1.80e-04 |
| 1757921 | 1.19e-07 | 5.05e-06 | 4.33e-06 | 1.96e-04 |
| 1757966 | 3.97e-07 | 5.05e-06 | 1.31e-05 | 1.88e-04 |
| autumn 1760378 | 7.69e-09 | 5.17e-06 | — | — |

**Every field, every probe, improves.** Mineral-N error drops by roughly
**7–670×**. Direction: unambiguously toward Fortran.

### 4d. N LOSS (leaching) — **CLOSED**: the totals now match to 4.4e-07

The autumn window turned out to be **unable to test leaching at all**: Fortran's
own `no3` pre→post change there is *identically zero* (frozen/dry column, no
drainage). That is itself a finding — it is why the summer window was generated.

In the summer window, `n_leaching!` reproduces Fortran's leaching operator
exactly. The *totals* originally disagreed (rel 1.0–16.4), which this document
attributed to a `qflx_surf` hydrology bug. **That diagnosis was wrong, and the
follow-up PR proved it: the surface-runoff port is faithful.**

> **The real cause was in THIS harness.** Bow's datm streams cycle
> `year_first=2002, year_last=2009, year_align=2002, taxmode='cycle'`, so the
> BGC-spinup windows — at Fortran model year **2202** — are forced by data year
> `2002 + mod(200, 8)` = **2002**. `fortran_parity_ncycle.jl` drove Julia from
> `clmforc.2003.nc` on a 2003 step date. Julia and Fortran were therefore seeing
> **different years' weather**, which decorrelates (corr ≈ 0, not a clean offset)
> everything downstream of precipitation. `qflx_surf` is `fsat · (rain+snowmelt)`,
> so it was the first place the wrong weather became visible, and `n_leaching!`
> the first consumer that multiplied it.
>
> `scripts/fortran_parity_cn_summer.jl` had already documented this exact rule;
> the N-cycle harness simply did not apply it. The mapping now lives in
> `parity_forcing()` (`scripts/fortran_parity_common.jl`) so it cannot be
> forgotten again.

With the correct forcing year, every N-loss number falls into line:

| metric | with 2003 forcing (wrong) | with 2002 forcing (correct) |
|---|---|---|
| leaching total, worst \|rel\| | **1.0 – 16.4** | **4.37e-07** |
| `smin_nh4_vr` | 2.0e-05 | **7.26e-07** |
| `smin_no3_vr` | 2.5e-07 | **5.55e-09** |
| `sminn_vr`    | 6.5e-07 | **1.88e-08** |

`ndep_to_sminn` and `ffix_to_sminn` remain EXACT, and the directional test still
improves at every probe. The leaching **operator** result quoted before (predict
levels 2…20 from one scalar calibrated on level 1, to ≤ 4.4e-07) was correct —
it was measuring the operator around the wrong water flux, and that same 4.4e-07
is now the error of the **absolute totals**.

**Surface-runoff (`qflx_surf`) parity: CLOSED.** Scored directly against a
200-step instantaneous Fortran history reference
(`scripts/fortran_parity_qflx_surf.jl`): `QOVER` max \|rel\| = **6.9e-09**,
`FSAT` = 5.8e-09, `RAIN` exact. See §7.

---

## 5. The CN balance check — turned on, and what it says

Both `begin_cn_column_balance!` / `begin_cn_gridcell_balance!` and
`c_balance_check!` / `n_balance_check!` are now wired into the facade, with the
`CNBalanceData` + `CNProductsData` instances the facade never had.

Wiring them exposed a **second dead link**: `cn_driver_summarize_states!` did not
call the SoilBiogeochem state summaries, so `totc_col` / `totn_col` — the check's
*begin* and *end* masses — were never computed at the point the balance is
seeded. `begcb`/`begnb` could only ever be zero, and the check would have
reported the entire column mass as an error. Fortran's `CNDriverSummarizeStates`
does call them; Julia's now does too.

### NITROGEN half: **RUNS and PASSES**

At every probe, with the N inputs visible in the check's own accounting for the
first time:

```
inputs_ndep = 3.6e-05      (= 1.0e-8 * 3600  ✓)
inputs_ffix = 9.9e-06      (✓)
inputs_nfix = 0.0          (✓ correct — under FUN, Fortran uses ffix, not nfix)
col_errnb  ≈ 1.3e-04 … 2.0e-04 gN/m2/step
```

`|err| ≈ 1.5e-4 gN/m²` is **below the error threshold (`nerror` = 1e-3)** — the
check passes — but **above the warning threshold (`nwarning` = 1e-7)**, so it
warns.

> **FINDING (reported, not gated): the N budget does not close to machine
> precision.** The residual is ~1.5e-4 gN/m²/step. The check's *output* side
> still reads `outputs_dnit = 0.0`, i.e. the column-level N loss aggregates
> (notably `f_n2o_nit_col`, and the vertical-transport N term) are not populated —
> the same dead-`p2c` bug class. That is the most likely home of the residual.
> It needs its own parity-validated PR; it is **not** fixed here and **not**
> hidden.

### CARBON half: **cannot run yet**

The check reads `gpp_col`, `er_col`, `col_fire_closs`, … — the column-level C
flux aggregates. **Every one of them is 0.0**: like `npp_col`, the patch→column
`p2c` that would fill them is never performed. So the C check sees zero fluxes
and reports the whole net storage change (0.0024 gC/m² per step) as an error.

This is a **real finding of the same bug class**, and it is why the check is
**default OFF** (`CLM.cn_balance_check_enabled!(true)` to turn it on): switching
it on globally today would abort every `use_cn` run for a *carbon* reason
unrelated to this PR. The N half is validated and passing; the C half needs the
column C-flux aggregates wired in a separate PR.

---

## 6. What is deliberately NOT wired, and why

| routine | status | why |
|---|---|---|
| `n_fert!`, `n_soyfix!` | **UNWIRED** | Crop N. Bow's reference has `use_crop = .false.` and **no crop CFT**, so there is no Fortran reference here to validate against. Both are ported + unit-tested. Wiring them blind is exactly what #218 correctly refused to do. **TODO:** needs a crop-active Fortran reference run (a site with a crop CFT). |
| `plant_crop!`, `vernalization!` | **UNWIRED** | Same reason — crop-gated, no Bow reference. |
| `cn_products_update!` | **UNWIRED** | Needs the harvest / gross-unrepresented-disturbance flux inputs, which are themselves not on the live path at Bow (`do_harvest = false`). The product **pools** are now allocated (the balance check needs them as a sink) but are not updated. **TODO:** needs a transient-landcover Fortran reference. |
| non-`use_nitrif_denitrif` branch of `SoilBiogeochemNStateUpdate1` | **NOT PORTED** | Julia's mineral-N state update implements only the nitrif/denitrif path, and `clm_initialize!` sets `use_nitrif_denitrif = use_cn`, so the other branch is unreachable in practice. Documented rather than written blind. |
| column-level C flux aggregates (`gpp_col`, `er_col`, …) | **UNWIRED** | Blocks the carbon half of the balance check — see §5. |
| `qflx_surf` parity | **CLOSED** | Was never a hydrology bug — it was this harness's forcing year. `QOVER` max \|rel\| = 6.9e-09 against a 200-step instantaneous Fortran history. See §4d and `scripts/fortran_parity_qflx_surf.jl`. |

---

## 7. Reproducing

Fortran reference (run **plainly** — never under lldb; macOS 26's xzone allocator
uses inline `brk #0x1` guards that halt lldb but run fine directly):

The case (`clm_bgc_spinup`) is left in its **default** state: the ZERO stream, so
every pre-existing reference stays reproducible. To regenerate the N-cycle
references, swap the stream and re-run:

```bash
cd .../SYMFLUENCE_data/clm_bgc_spinup

# 1. point lnd_in at the SYNTHETIC stream (the only edit needed; nrevsn and the
#    three rpointers are already consistent at the 2202-07-16-57600 restart)
sed -i '' 's/fndep_clm_ZERO_/fndep_clm_UNIFORM1e-8_/' lnd_in

# 2. shorten the run (nuopc.runconfig: stop_n=200, restart_n=100000 so no
#    restart is written and the rpointers are not clobbered)

# 3. run PLAINLY — never under lldb (macOS 26's xzone allocator uses inline
#    `brk #0x1` guards that halt lldb but run fine directly)
source .../cases/symfluence_build/.env_mach_specific.sh
export DYLD_INSERT_LIBRARIES=/opt/homebrew/lib/libmimalloc.dylib
export HWLOC_COMPONENTS=-opencl          # libhwloc -> OpenCL -> Metal SIGILL at MPI startup
env PDUMP_NSTEP_LO=1757873 PDUMP_NSTEP_HI=1758060 .../bld/cesm.exe

# 4. restore lnd_in + nuopc.runconfig afterwards, and collect the dumps into
#    bgc_ref_ndep_summer/ (each reference dir keeps the lnd_in.used that made it)
```

CIME's `case.build` needs a Python with `distutils` (removed in 3.12+); build with
`/usr/bin/python3` (3.9) on `PATH`.

The new `after_ecosysdyn_postdrain` dump boundary lives in the case SourceMods
(`clm_driver.F90`, immediately after `EcosystemDynamicsPostDrainage`).

Julia:

```bash
julia --project=. scripts/fortran_parity_ncycle.jl summer 5
julia --project=. scripts/fortran_parity_ncycle.jl autumn 5
```
