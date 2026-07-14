# Surface runoff (`qflx_surf`) — Fortran parity

**Status: CLOSED. The surface-runoff port is faithful.** `qflx_surf` agrees with
CTSM to **6.9e-09** relative — the Float64 cross-compiler floor.

The "surface-runoff bug" that PR #221 reported did not exist. It was an artifact
of #221's own parity harness driving Julia with **the wrong forcing year**.

Nothing was tuned to match Fortran.

---

## 1. The reported bug, and what it actually was

PR #221 validated `n_leaching!` and found the leaching *operator* exact but the
leaching *totals* off by 1.0–16.4 relative. It correctly localized the cause to
the water flux — `qflx_surf` — and (reasonably) concluded there was a hydrology
bug, naming the dead `compute_fsat_topmodel!` / `compute_fsat_vic!` as prime
suspects.

Both leads were wrong:

**(a) The dead `compute_fsat_*` functions are a red herring.** They *are* dead —
nothing calls them. But `saturated_excess_runoff!` **inlines** the identical
TOPMODEL/VIC logic into its GPU kernel (`_sat_excess_runoff_kernel!`,
`src/biogeophys/sat_excess_runoff.jl`). `fsat` is computed, and it is *correct*:
it agrees with Fortran to **5.8e-09**. The two standalone functions are simply an
unused duplicate of live, correct code. **A dead function is not automatically a
dead code path** — check whether its body was inlined somewhere before concluding
the physics is missing.

**(b) The real cause was the harness's forcing year.**

Bow's `datm.streams.xml` — all three CLMNCEP streams:

```
year_first = 2002   year_last = 2009   year_align = 2002   taxmode = cycle
```

so a model year `Y` is driven by **data year**
`year_first + mod(Y - year_align, year_last - year_first + 1)`.

The BGC-spinup reference windows sit at Fortran model year **2202**:

```
data year = 2002 + mod(2202 - 2002, 8) = 2002 + mod(200, 8) = 2002
```

`scripts/fortran_parity_ncycle.jl` drove Julia from **`clmforc.2003.nc`** on a
2003 step date. Julia and Fortran were seeing **different years' weather**.

This is a nasty failure mode because it does *not* look like a forcing bug. The
two precipitation series are not offset, they are **uncorrelated** (corr ≈ −0.03
at every lag from −3 to +3 steps), so the resulting flux mismatch looks like
broken physics: at one step Julia's `qflx_surf` was 17× Fortran's, at another
1/36500 of it.

Feeding the *correct* year (`clmforc.2002.nc`) makes the two precipitation series
agree at **corr = 0.99999, lag 0**, value-for-value.

`scripts/fortran_parity_cn_summer.jl` had **already documented this exact rule**.
The N-cycle harness just did not apply it. The mapping is now a shared, mandatory
helper — `parity_forcing(model_date)` in `scripts/fortran_parity_common.jl` — so
it cannot be forgotten again.

---

## 2. The Fortran reference

The restart-format `pdump` snapshots carry *state*, not *fluxes*, so they cannot
score `qflx_surf` at all. Instead this uses CTSM's **history** mechanism at
per-step instantaneous resolution:

```
hist_nhtfrq          = -8760, 1     ! tape 2 = every step
hist_mfilt           = 20, 400
hist_avgflag_pertape = 'A', 'I'     ! tape 2 = INSTANTANEOUS, not a time mean
```

| item | value |
|---|---|
| run dir | `SYMFLUENCE_data/clm_ref_water_summer` |
| file | `Bow_at_Banff_lumped.clm2.h1.2202-07-16-61200.nc` |
| records | 200 contiguous steps, `dtime = 3600 s` |
| record ↔ nstep | record `i` ⇔ `nstep = 1757872 + i` (record 1 ends 2202-07-16 17:00 = the step starting 16:00 = nstep 1757873) |

Eight of the scored fields are **not** stock CTSM history fields. They were added
as `hist_addfld1d` registrations in a SourceMods copy of `WaterFluxBulkType.F90`
so the runoff partition could be bisected term by term:

`QSATEXCS` `QINFLEXC` `QINFLEXCS` `QINSOIL` `QINSOILL` `QTS2SFC` `QINH2OSFC` `QSFCDRAIN`

Harness: `scripts/fortran_parity_qflx_surf.jl [nprobe] [h2osfcflag]`.

> **Boundary discipline.** The h1 record is written at the END of the step, so it
> is the SAME-step value of everything computed during that step. `before_step` is
> used only as the injected IC — never as a comparison target.

---

## 3. Scorecard — the whole surface-water chain

Julia single-step oracle (inject Fortran `before_step`, run one `clm_drv!` step)
vs the same-step Fortran h1 record, 8 probes across the summer window:

| field | max \|abs\| | max \|rel\| | verdict |
|---|---|---|---|
| `FSAT` | 5.0e-09 | **5.8e-09** | exact |
| `FCOV` | 5.0e-09 | 5.8e-09 | exact |
| `RAIN` (forcing) | **0.0** | **0.0** | bit-exact |
| `QRAINSNOM` (`qflx_rain_plus_snomelt`) | 7.3e-15 | 3.5e-09 | exact |
| `QSATEXCS` (`qflx_sat_excess_surf`) | 1.3e-14 | 6.9e-09 | exact |
| `QTOPSOIL` (`qflx_top_soil`) | 2.5e-14 | 1.2e-08 | exact |
| **`QOVER` (`qflx_surf`)** | **1.3e-14** | **6.9e-09** | **exact** |
| `QINFLEXC` / `QINFLEXCS` | 0.0 | 0.0 | zero in both |
| `QH2OSFC` / `FH2OSFC` / `H2OSFC` | 0.0 | 0.0 | zero in both |
| `QTS2SFC` / `QINH2OSFC` / `QSFCDRAIN` | 0.0 | 0.0 | zero in both |
| `QDRAI` / `QRGWL` | 0.0 | 0.0 | zero in both |
| `ZWT` | 2.9e-08 | 1.3e-08 | exact |
| `QINSOIL` / `QINSOILL` / `QINFL` | 4.3e-10 | 5.4e-05 | see below |

Both codes satisfy the partition identity exactly:

```
QOVER == QSATEXCS + QINFLEXCS + QH2OSFC        (both codes, every step)
QSATEXCS == FSAT * qflx_rain_plus_snomelt      (both codes, every step)
```

Before the forcing-year fix, `QOVER` max \|rel\| was **0.9987**, with per-step
ratios of 0.000, 0.001, Inf and 1.239. After: **6.9e-09**, ratio 1.000 at every
probe.

### Downstream: `n_leaching!` totals

#221 predicted the leaching totals would come into line once the water flux did.
They do (`scripts/fortran_parity_ncycle.jl summer`):

| metric | 2003 forcing (wrong) | 2002 forcing (correct) |
|---|---|---|
| leaching total, worst \|rel\| | **1.0 – 16.4** | **4.37e-07** |
| `smin_nh4_vr` | 2.0e-05 | **7.26e-07** |
| `smin_no3_vr` | 2.5e-07 | **5.55e-09** |
| `sminn_vr` | 6.5e-07 | **1.88e-08** |

`ndep_to_sminn` / `ffix_to_sminn` remain EXACT.

### The one honest residual: `QINSOIL` (5.4e-05 rel, 4.3e-10 abs)

`qflx_in_soil` (and its dependents `qflx_in_soil_limited`, `qflx_infl`) carry a
~4e-10 mm/s residual. It does **not** enter `qflx_surf` — the surface-runoff
partition takes `qflx_sat_excess_surf` and `qflx_infl_excess_surf`, both of which
are exact. `qflx_in_soil` subtracts an evaporation term
(`(1 - fsno - frac_h2osfc)·qflx_ev_soil`), so this is the tail of the known
canopy/ground energy-balance residual (the same coupled-solve floor that leaves
`T_VEG` at 2.8e-04 in the SP anchor), not a hydrology defect.

---

## 4. The real bug this work DID find: `h2osfcflag` default

CTSM's namelist default is **`h2osfcflag = 1`** (`SoilHydrologyType.F90`,
`clm_soilhydrology_inparm`) — the surface-water (h2osfc) store is ACTIVE unless a
case turns it off. Bow's `lnd_in` never sets it, so **the Fortran reference runs
with h2osfc active.**

In CLM.jl, `SoilHydrologyData` itself already defaulted to 1 — but **both driver
entry points overrode it to 0**:

* `src/driver/clm_initialize.jl` — `h2osfcflag::Int = 0`
* `src/driver/clm_run.jl` — `h2osfcflag::Int = 0`

So every caller that did not explicitly pass the flag ran with the surface-water
store **silently disabled**: infiltration excess went straight to surface runoff
instead of being ponded and re-infiltrated, and `frac_h2osfc` was pinned at 0.
`scripts/parity_run_domain.jl` (the multi-biome scorecard) had already been forced
to pass `h2osfcflag = 1` explicitly to work around this, with the comment
*"CLM5 default (surface water active); matters for wet sites"* — the library
default was the outlier.

**Both entry points now default to 1**, matching CTSM.

### Why this did NOT fix `qflx_surf` at Bow (and why that is not a let-down)

Measured, not assumed: at Bow the h2osfc pond **never activates in either code**.
Fortran's own history says so — `H2OSFC`, `FH2OSFC`, `QH2OSFC`, `QINFLEXC` and
`QINH2OSFC` are **identically zero at all 200 steps**, because `qinmax` at Bow is
large enough that there is never any infiltration excess to pond. So at Bow the
flag is inert, and flipping it is a **byte-for-byte no-op**:

* SP parity anchor (`fortran_parity_validate.jl 13461`): global max \|rel\|
  **2.831e-04 → 2.831e-04**, identical.
* Bow scorecard run (`parity_run_domain.jl`): output identical to the digit.

The flag fix is therefore a **correctness fix that Bow cannot exercise**. It will
change results at wet / low-relief sites, where `qinmax` is small enough to
generate infiltration excess and the pond does fill. That is the intended CLM5
behaviour and matches CTSM.

---

## 5. Independently confirmed (fixed upstream by #223): `begcb_grc` OOB read

While validating this work, running the suite under `--check-bounds=yes` (which CI
uses) surfaced a `BoundsError` in `_cnbal_begin_grc_kernel!`
(`src/biogeochem/cn_balance_check.jl`) at **6 separate tests**. Root cause: the
three OPTIONAL dribble-flux terms (`hrv_xsmrpool_ / gru_conv_cflux_ /
dwt_conv_cflux_amount_left`) default to `Float64[]`, the only live caller
(`vegetation_facade.jl`) never passes them, yet the `!use_fates_bgc` branch read
all three unconditionally. Under `@inbounds` that does not fail — it silently pulls
three words of **out-of-bounds garbage** into `begcb_grc`, the carbon balance
check's *begin* mass.

**This was fixed upstream in #223** (via a `has_dribbler` guard passed into the
kernel) while this branch was in flight, so no fix is carried here — the two
diagnoses agree. Recorded because it is a clean instance of the banked
`check-bounds CI trap`: `@inbounds` OOB reads pass in a default build and only fail
under `--check-bounds=yes`. **Test new kernel paths with `--check-bounds=yes`.**

---

## 5b. Is this the same residual as #224's top-soil-layer melt plateau?

**No — and the measurement rules it out for the surface-runoff chain.**

#224 reports that with FATES's root profile fixed, its BTRAN reads almost entirely
off CLM soil layer 1, and Julia's top layer leaves the 273.15 K melt plateau ~4 h
early, holding `h2ovol1` 0.31 vs Fortran's 0.35. Both that and `qflx_surf` are
"top-of-column water/energy", so a shared cause is a reasonable hypothesis.

It is not shared, at least not anywhere `qflx_surf` can see it:

* In the scored window `qflx_surf` and every term feeding it (`FSAT`, `ZWT`,
  `QRAINSNOM`, `QSATEXCS`, `QTOPSOIL`) agree with Fortran to **≤ 1.2e-08**. A top-
  layer thermal/phase-change error large enough to move `h2ovol1` by 0.04 does not
  propagate into the surface-runoff partition here.
* The window is **snow-free and melt-free**: Fortran's own `FSNO` and `QSNOMELT`
  are zero (and Julia matches them exactly, 0.0 abs). So this window **cannot
  exercise** a melt-plateau residual at all.

So the two are separable: the melt-plateau residual is a *soil-temperature /
phase-change* defect that this (summer, snow-free) reference has no power to test,
and it demonstrably does not contaminate `qflx_surf`. Reproducing it needs a
**melt-season** Fortran reference (spring, `QSNOMELT > 0`), which is a different
window and a different PR.

The one honest residual this window *does* show — `QINSOIL` at 4.3e-10 mm/s — is an
**evaporation**-term residual (`qflx_in_soil` subtracts
`(1 - fsno - frac_h2osfc)·qflx_ev_soil`), not a phase-change one, and it does not
enter `qflx_surf`.

---

## 6. Pre-existing, NOT introduced here: the Bow water-balance failure

`DOMAIN=Bow scripts/parity_run_domain.jl` **aborts on the live+fatal water-balance
check** at `nstep = 2041` with `errh2o = 1.066869312643659e-6`.

This reproduces **bit-identically on `main` (21eb9ca)** with these changes stashed,
so it is pre-existing and is *not* a regression from this work — but it is real
and it is open. It is not a `qflx_surf` defect (`qflx_surf` there is 6.5e-7 mm/s,
and is now proven exact against Fortran); it is a separate balance-closure bug,
most likely surfaced when #211 made the check live and fatal. **It needs its own
PR.**

---

## 7. Reproducing

Fortran reference (run **plainly** — never under lldb; macOS 26's xzone allocator
uses inline `brk #0x1` guards that halt lldb but run fine directly):

```bash
CASE=.../installs/clm/cases/symfluence_build

# 1. History registrations for the runoff-chain internals
cp .../src/biogeophys/WaterFluxBulkType.F90 "$CASE/SourceMods/src.clm/"
#    add hist_addfld1d for QSATEXCS/QINFLEXC/QINFLEXCS/QINSOIL/QINSOILL/
#    QTS2SFC/QINH2OSFC/QSFCDRAIN at the end of InitBulkHistory

# 2. Build. NOTE: case.build needs a python WITH distutils (3.11), not 3.12+.
SHIM=$(mktemp -d); ln -sf "$(command -v python3.11)" "$SHIM/python3"; export PATH="$SHIM:$PATH"
rm -f "$CASE"/bld/*/*/*/*/*/*/clm/obj/waterfluxbulktype.o
cd "$CASE" && ./case.build

# 3. Run into an ISOLATED rundir (never clobber the validated dumps)
RUN=.../clm_ref_water_summer; mkdir -p "$RUN"; cd "$RUN"
cp .../clm_bgc_spinup/{lnd_in,datm_in,datm.streams.xml,drv_in,drv_flds_in,fd.yaml,\
   nuopc.runconfig,nuopc.runseq,rpointer.*} .
cp .../clm_bgc_spinup/Bow_at_Banff_lumped.{clm2.r,clm2.rh0,cpl.r,datm.r}.2202-07-16-57600.nc .
#    ^^^ the .rh0 (restart-history) file is REQUIRED for a branch run; without it
#    CTSM dies with a bare SIGTRAP and no error message.
#    then set hist_nhtfrq/-mfilt/-fincl2/-avgflag_pertape as in §2, stop_n = 200
source "$CASE/.env_mach_specific.sh"; export NETCDF_PATH=/usr/local
export HWLOC_COMPONENTS=-opencl
"$CASE/bld/cesm.exe"          # runs plainly, exit 0
```

Julia side:

```bash
julia +1.12 --project=. scripts/fortran_parity_qflx_surf.jl 8 1   # qflx_surf chain
julia +1.12 --project=. scripts/fortran_parity_ncycle.jl summer 6 # leaching totals
```

### Gotchas banked

* **The `.rh0` file is mandatory for a CTSM branch run.** Omitting it produces a
  bare `SIGTRAP` right after the decomposition printout, with **no CLM error
  message** — it looks exactly like a compiler/allocator crash. It is not.
* **`case.build` needs Python ≤ 3.11** (CIME imports `distutils`). Under 3.12+ it
  dies in `standard_script_setup`; if you pipe the output you will see a spurious
  exit 0 and silently run a stale exe.
* **Do not add `hist_addfld1d` to the BASE `WaterFluxType.F90`.** That type is
  instantiated per water tracer; a registration there traps during init. The
  bulk-only `WaterFluxBulkType.F90` is the correct place.
* **Always route parity forcing through `parity_forcing()`.** Hand-picking a
  forcing file is how #221 got a phantom hydrology bug.
