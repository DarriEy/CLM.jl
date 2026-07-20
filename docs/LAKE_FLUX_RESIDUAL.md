# Lake surface turbulent-flux residual (backlog B1)

Open residual: on the cold-start PCT_LAKE=100 SP reference
(`scripts/validation/gen_lake_ref.sh` → `SYMFLUENCE_data/clm_lake_run`), the lake
THERMODYNAMICS track Fortran (TLAKE / LAKEICEFRAC ~1%) but the SURFACE TURBULENT-FLUX
solve does not: FSH / EFLX_LH ~10%, TSA ~7-11%. (The FSA ~13x spikes are relative-diff
artifacts on near-zero dawn/dusk solar — not real.)

The same residual was chased through three real fixes ending at `3db0c1a`
("closes the surface-flux residual"), so it has recurred or was never fully closed.

## STEP 1 — namelist flag audit (Fortran `lnd_in` vs what the harness passes)

Rationale: a ~10% flux residual sitting on top of VALIDATED thermodynamics is the
classic signature of a mismatched INPUT, not a mis-ported kernel. In this repo an
unset flag has been the actual cause seven times (#233, #240, #243, #248, #251, the
crop data-locality case, and the qflx_surf case in #225). **An UNSET flag is not a
matched flag** — if the harness does not pass it, the port takes its own default,
and several driver defaults changed in the last week (#252, #259, #225, #265, #267).

Sources:
- Fortran: `/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_lake_run/lnd_in`
  (regenerated 2026-07-18 by `scripts/validation/gen_lake_ref.sh`, rc=0)
- Harness: `scripts/fortran_parity_lake.jl` — the `CLM.clm_initialize!` call and the
  `CLM.CLMDriverConfig` construction
- Port defaults: `src/driver/clm_initialize.jl:65-130`, `src/constants/varpar.jl:148`,
  `src/types/friction_velocity.jl:187`, `src/biogeophys/lake_con.jl:79-93`

Legend: **SET** = harness passes it explicitly. **DEFAULT** = harness is silent and the
port's own default applies (the dangerous case).

### Flags that reach the lake surface-flux path

| Flag | Fortran `lnd_in` | Harness | Port default | Verdict |
|---|---|---|---|---|
| `zetamaxstable` | `0.5d00` | DEFAULT | `0.5` (friction_velocity.jl:187) | MATCH |
| `wind_min` | paramfile (`clm5_params.nc`) | from paramfile | read via read_params.jl:104 | MATCH (same paramfile) |
| `z0param_method` | `'ZengWang2007'` | DEFAULT | `"ZengWang2007"` (varctl.jl:163) | MATCH — see note Z |
| `nlevsno` | `12` | DEFAULT | `12` (varpar.jl:148) | MATCH |
| `soil_layerstruct_predefined` | `'20SL_8.5m'` | DEFAULT | `"20SL_8.5m"` | MATCH |
| `use_subgrid_fluxes` | `.true.` | DEFAULT | n/a (single gridcell) | inert here |
| `lakepuddling` | absent → CTSM default `.false.` | DEFAULT | `false` (lake_con.jl:93) | MATCH |
| `lake_use_old_fcrit_minz0` | absent → `.false.` | DEFAULT | `false` (lake_con.jl:79) | MATCH |
| `deepmixing_depthcrit` | absent → `25.0` | DEFAULT | `25.0` (lake_con.jl:82) | MATCH |
| `lake_melt_icealb` | absent → paramfile/default | DEFAULT | — | MATCH |

Note Z: the whole `&clm_lake_inparm` group is ABSENT from the reference `lnd_in`, so the
Fortran run takes CTSM's lake defaults — which the port's `lake_con.jl` defaults
reproduce. The lake namelist is therefore **not** the source of the residual.

**`z0param_method` is a DEAD END — corrected after first writing this table.** An
earlier revision of this document flagged it as UNSET. That was wrong, for two
independent reasons, and it is recorded here so nobody chases it again:
1. `varctl.z0param_method` already defaults to `"ZengWang2007"` (`src/constants/varctl.jl:163`,
   set by #260), matching the reference. The `""` defaults visible in
   `src/types/friction_velocity.jl:303,566` are per-function argument defaults, not the
   driver default — reading those first is what produced the false flag.
2. More decisively, **`lake_fluxes.jl` never consults `z0param_method` or `z0method` at
   all** (zero occurrences). It applies the ZengWang2007 form unconditionally at
   `lake_fluxes.jl:259`. The flag cannot move the lake result in either setting.

### Flags that differ

| Flag | Fortran `lnd_in` | Harness | Verdict |
|---|---|---|---|
| `int_snow_max` | `3113.2227d00` (`&scf_swenson_lawrence_2012_inparm`) | **SET `2000.0`** | **MISMATCH** |
| `h2osfcflag` | absent → CTSM default `1` | **SET `0`** | **MISMATCH** (port default is now 1, #225) |
| `use_luna` | `.true.` | SET `false` | MISMATCH (expected inert: no veg patch on a 100% lake) |
| `use_hydrstress` | `.true.` | SET `false` | MISMATCH (expected inert: same) |
| `n_melt_glcmec` | `10.0d00` | DEFAULT | to confirm |

### Flags verified matched

| Flag | Fortran `lnd_in` | Harness / port | Verdict |
|---|---|---|---|
| `use_bedrock` | `.true.` | SET `true` | MATCH (survives the #252 conditional-default change because the harness sets it) |
| `use_aquifer_layer` | derived `.false.` (soilwater_movement_method=1 → lower_boundary_condition=2) | SET `false` | MATCH |
| `create_crop_landunit` | `.true.` | derived `!use_fates` = `true` (#265) | MATCH |
| `use_cn` | `.false.` | SET `false` | MATCH |
| `use_crop` | `.false.` | DEFAULT `false` | MATCH |
| `use_lch4` | `.false.` | DEFAULT `false` | MATCH |
| `use_fates` | `.false.` | DEFAULT `false` | MATCH |
| `use_nitrif_denitrif` | `.false.` | DEFAULT (use_cn=false) | MATCH |
| `use_excess_ice` | `.false.` | DEFAULT `false` | MATCH |
| `use_hillslope` | `.false.` | DEFAULT `false` | MATCH |
| `irrigate` | `.false.` | DEFAULT `false` | MATCH |
| `finidat` / `nrevsn` | `''` / `''` (cold start) | cold start | MATCH |
| `dtime` | 3600 (hist_nhtfrq=1, 48 steps) | SET `3600` | MATCH |
| `fsurdat` | `surfdata_lake100.nc` | same file (repo `test_inputs/lake/`) | MATCH |
| `paramfile` | `clm5_params.nc` (Bow dds_run_1) | same path | MATCH |
| `fsnowoptics` / `fsnowaging` | snicar 5bnd / drdt_bst_60 | SET, same paths | MATCH |
| `co2_ppmv` | `367.0` | to confirm | — |

### Non-namelist inputs (the harness supplies these directly)

| Input | Fortran | Harness | Verdict |
|---|---|---|---|
| `forc_hgt_u/t/q` | `Sa_z = 30.0` | hard-set `30.0` | **MATCH — verified** |
| TBOT/WIND/QBOT/PSRF/FLDS time interp | `tintalgo = linear` | exact record read | MATCH at exact stream times |
| FSDS time interp | `tintalgo = coszen` | — | differs off-record (dt == stream interval, so inert) |
| PRECT time interp | `tintalgo = nearest` | — | inert |
| stream `<offset>` | `0` (all streams) | n/a | MATCH |
| forcing hour mapping | datm advances with the model clock | `base + Hour(max(s-1, 1))` | **SUSPECT — see below** |

`Sa_z` verified at `components/cdeps/datm/datm_datamode_clmncep_mod.F90:416` —
`if (.not. associated(strm_z)) Sa_z(n) = 30.0_r8`. The Bow forcing file
(`clmforc.2003.nc`) carries no `ZBOT`/height variable, so this fallback fires and the
Fortran reference genuinely runs at 30 m. The harness's hard-set 30.0 is correct.

Stream settings read from `clm_lake_run/datm.streams.xml` (`CLMNCEP.TPQW` block).
Because the model timestep (3600 s) equals the forcing interval (hourly), `linear`
interpolation evaluated at a stream time returns that record exactly — so the
interpolation algorithm is not a source of error, but **which hour is evaluated** is.

**Forcing-hour mapping.** The harness maps step `s` to forcing hour `max(s-1, 1)`:
steps 1 and 2 BOTH read hour 1, and every step from 2 onward is therefore offset one
hour behind the Fortran datm. The same expression appears in
`scripts/fortran_parity_stillwater_multistep.jl:133`, so there is precedent — but that
precedent was established for a SINGLE-step probe (`fortran_parity_stillwater.jl:48`
"Fixed via force_date=base+Hour(max(N-1,1))"), where `max(N-1,1)` only guards `N=1`.
In a LOOP the same expression duplicates hour 1 and lags the rest. A one-hour forcing
lag is exactly the failure that produced the phantom fire residual in #233 (NFIRE
4.76 → 1.3e-15 once the harness year/hour mismatch was fixed; the port was correct).

A one-hour lag in wind / air temperature / humidity would bias FSH, EFLX_LH and TSA —
the three fields that are off — while leaving the deep-lake thermodynamics (a slow
integrator with a ~day-plus time constant) visibly unaffected. That matches the
observed signature better than a mis-ported profile function does.

## STEP 2 — what the residual actually is

Two findings, both in the HARNESS rather than in the lake physics. Neither is a
mis-ported kernel.

### Finding A — the h0 record index was off by one (the FSH / EFLX_LH residual)

The lake reference's history file indexes its own records unambiguously:

```
nstep       = [0, 1, 2, 3, ...]
time        = [0.0, 0.0416667, 0.0833333, 0.125]        days since 2003-01-01
time_bounds = [[0,0], [0,1/24], [1/24,2/24], [2/24,3/24]]
```

**Record 1 is the `nstep=0` cold-start dump over a ZERO-WIDTH interval `[0,0]`** — an
instantaneous initial-state write, not a timestep average. Record 2 (`nstep=1`,
bounds `[0, 1h]`) is the first actual model step. So **model step `s` is h0 record
`s+1`**, and the harness's `fds[fn][1, s]` indexing compared every Julia step against
the Fortran state ONE STEP BEHIND it.

That single off-by-one accounts for the bulk of the reported FSH / EFLX_LH
divergence. Cross-check on the raw values (`LAKE_DUMP=1`), before any change:

| | Julia step 1 | h0 record 1 (`nstep=0`) | h0 record 2 (`nstep=1`) |
|---|---|---|---|
| FSH | 69.999 | 84.913 | **70.255** |
| EFLX_LH | 44.276 | 53.756 | (see run output) |

Julia's step 1 FSH matches record 2 to **0.4%**, against the 17.6% "residual" that the
record-1 comparison reported. The port was right; the comparison was misaligned.

This is the same class as the trap already recorded in MEMORY as
*parity-before-vs-after-dump-trap* ("wrong dump manufactured a phantom bug"), and the
eighth instance of the harness-input trap in this repo.

**On the `LAKE_FORC_OFFSET = −1` measurement.** Shifting the FORCING back one hour
appeared to cut `max|rel|` from 1.30e+01 to 9.6e-01. That number is **not** evidence
and was not acted on: `max|rel|` over the run is dominated by the FSA dawn/dusk
relative-diff artifact on near-zero solar (a known artifact, see the header of
`gen_lake_ref.sh`), so shifting solar by an hour changes which near-zero denominator
the maximum lands on. A one-hour forcing shift was ALSO partially compensating for
Finding A by corrupting a second input — two wrongs moving the aggregate metric.
The forcing mapping itself is fine: step `s` spans `(s-1)h → s·h` and reads its
forcing accordingly; the only genuine defect there is the `max(s-1, 1)` clamp, which
makes steps 1 and 2 both read hour 1 (`LAKE_FORC_EXACT=1` drops the clamp).

### Finding B — `t_ref2m_patch` is never written on a lake patch (the TSA residual)

`TSA` is not a 7-11% physics divergence. It is an **unwritten field**:

```
[dump s= 1] TSA J= 283.000 F= 259.463
[dump s=24] TSA J= 283.000 F= 264.752
[dump s=48] TSA J= 283.000 F= 270.141
```

Julia's `t_ref2m_patch` is **exactly 283.000 at all 48 steps** — the cold-start
initialisation constant — while Fortran's TSA ranges over 254.98 – 270.14 K.

`grep -rln 't_ref2m_patch\[' src/biogeophys/` returns only `bareground_fluxes.jl` and
`urban_fluxes.jl`. Neither runs on a lake patch, so nothing ever writes the field
there. CTSM's `LakeFluxesMod.F90` does compute it, carrying the full stability-
dependent profile relations out of the Monin-Obukhov iteration.

This is a genuine port gap of the **"ported then never called"** master class
(`dead-initcold-systemic` in MEMORY) — a dead write, not a physics error. It is also
the one place where the dead agents' 2-m-diagnostic hypothesis was pointing at
something real, though for the wrong reason: the diagnostic is not computed from a
single-regime profile, it is not computed at all.

A fix is lake-landunit-gated by construction: nothing writes `t_ref2m_patch` on a lake
patch today, so writing it there cannot change any non-lake result.

### Measured effect of the record fix

`LAKE_DUMP=1`, 47 steps (step `s` needs record `s+1`, so 47 not 48). `LAKE_REC_SHIFT`
and `LAKE_FORC_OFFSET` / `LAKE_FORC_EXACT` are env knobs on the harness.

| config | record map | forcing map | max\|rel\| |
|---|---|---|---|
| historical | `s` | `max(s-1,1)` | **1.298e+01** |
| forcing shifted −1h (rejected) | `s` | `max(s-1,1) − 1` | 9.617e-01 |
| **record fix (adopted)** | **`s+1`** | **`max(s-1,1)`** | **3.108e-01** |
| record fix + forcing `s` | `s+1` | `s` | 2.012e+00 |

**42× reduction** from the record fix alone, with the forcing mapping left untouched.
Note the last row: "correcting" the forcing on top of the record fix makes it *worse*,
which independently confirms the historical forcing mapping was already right and that
the −1h shift was compensation, not a fix.

Sample raw values after the record fix (Julia vs Fortran):

| step | FSH J / F | EFLX_LH J / F | TG J / F |
|---|---|---|---|
| 1 | 69.999 / 70.255 | 44.276 / 41.737 | 272.899 / 269.497 |
| 23 | 40.004 / 36.954 | 32.646 / 32.794 | 270.844 / 271.078 |
| 44 | 18.275 / 18.329 | 20.953 / 22.809 | 272.016 / 272.353 |
| 47 | 10.488 / 10.675 | 19.226 / 20.893 | 272.268 / 272.456 |

### Finding C — the lake MO kernel implements ONE of Fortran's FOUR stability regimes

A real residual SURVIVES correct alignment: FSH ~6-12% (31% at step 2), EFLX_LH ~6-10%,
plus H2OSNO ~10-20% in the last few steps. Part of that is Finding B (TSA is dead and
still contributes its ~7-11%), but the flux bias is separate, and there is a concrete
kernel-level cause.

`FrictionVelocityMod.F90:960-1050` computes the temperature/humidity profile relations
under a **four-way branch on `zeta`**:

1. `zeta < -zetat` (very unstable) — log capped at `-zetat*obu/z0h`, plus a convective
   correction `0.8*((zetat)^-0.333 - (-zeta)^-0.333)`
2. `zeta < 0` (unstable) — plain `log(zldis/z0h) - psi2(zeta) + psi2(z0h/obu)`
3. `0 <= zeta <= 1` (stable) — `log(zldis/z0h) + 5*zeta - 5*z0h/obu`
4. `zeta > 1` (very stable) — `log(obu/z0h) + 5 - 5*z0h/obu + (5*log(zeta)+zeta-1)`

`lake_fluxes.jl:336-337` implements **regime 2 only**:

```julia
temp1 = VKC / (log(zldis_t / z0hg) - stability_func2(zeta_t) + stability_func2(zeta0h))
```

and `stability_func2` (`friction_velocity.jl:357-364`) clamps `z = min(zeta, 0)`, so it
is a no-op above neutral — it is explicitly documented "only valid for zeta <= 0".
`zeta_t` is clamped to `[-100, ZETAMAX_LAKE=0.5]`, so it genuinely can be positive.

Consequence for THIS reference: a 273 K lake under 255 K air with light wind is
strongly unstable, so **regime 1 is the one being missed** — the port applies the plain
unstable form at large negative `zeta` where CTSM caps the log and adds the convective
term. That changes `temp1` → `rah` → FSH/EFLX_LH directly, in the right magnitude for a
~10% flux bias. Regimes 3/4 are additionally unreachable in the port but do not fire in
this particular winter-lake case.

## STEP 3 — B and C fixed; C exposed a FOURTH defect (D)

Both port defects are now fixed. Each was measured on its own run so neither is
tuned against the other. Harness: `scripts/fortran_parity_lake.jl 47`, per-field
max|rel| over the run.

| field | baseline | after B | after C (= B+C) |
|---|---|---|---|
| FSH | 3.10e-01 | 3.10e-01 | **2.70e-01** |
| EFLX_LH | 2.80e-01 | 2.80e-01 | **2.70e-01** |
| TSA | 1.10e-01 | **8.30e-03** | 7.40e-03 |
| TLAKE | 6.50e-03 | 6.50e-03 | 1.30e-02 |
| LAKEICE / LAKEICE_SRF | 1.80e-02 | 1.80e-02 | **3.80e-01** |
| TG | 1.30e-02 | 1.30e-02 | 3.30e-02 |
| FIRE | 5.00e-02 | 5.00e-02 | 1.40e-01 |
| run max\|rel\| | 3.108e-01 | 3.108e-01 | 4.966e-01 |

### B — `t_ref2m` wired (LakeFluxesMod.F90:707-717)

CTSM computes the lake 2-m diagnostics in the loop AFTER the stability
iteration, from the profile coefficients the iteration converged on:

```
t_ref2m = thm + temp1*dth*(1/temp12m - 1/temp1)
q_ref2m = forc_q + temp2*dqh*(1/temp22m - 1/temp2)
rh_ref2m = min(100, q_ref2m/QSat(t_ref2m)*100)
```

`temp12m`/`temp22m` are the same profile relations evaluated at `zldis = 2 + z0h`
(`FrictionVelocityMod.F90:1010-1050`), and `dth`/`dqh` are as last set INSIDE the
iteration (F90:535-536) — i.e. before the Phase-3 freeze corrections touch
`t_grnd`, so they must be carried out of the loop rather than recomputed.

Effect: **TSA 1.10e-01 -> 8.30e-03** (a flat 283.000 K against a 255-270 K
reference, i.e. 13-28 K, becomes max ~1.9 K on the cold-start transient and
~0.1 K by the end of the run). Every other field is BYTE-IDENTICAL, which is the
expected signature of a pure dead-write fix: nothing else reads the field.

Note lake `z0` is a function of `ustar` (F90:570-582 Subin/Charnock for open
water, ZengWang2007 off `ust_lake` when frozen), unlike the land path's fixed
`zlnd`/`zsno` — so `temp12m` had to be evaluated inside the iteration where the
converged `z0hg` lives, not from a stored roughness afterwards.

### C — all four stability regimes ported

`FrictionVelocityMod.F90` branches four ways on `zeta`, identically for the
momentum profile (`:847-861`, transition `zetam = 1.574`) and the
temperature/humidity profile (`:946-960`, transition `zetat = 0.465`):

1. `zeta < -zeta*` very unstable — log capped at the transition + a convective
   correction (`1.14*((-zeta)^0.333 - zetam^0.333)` momentum,
   `0.8*(zetat^-0.333 - (-zeta)^-0.333)` heat)
2. `zeta < 0` unstable — plain log + `StabilityFunc`
3. `0 <= zeta <= 1` stable — `log + 5*zeta - 5*z0/obu`
4. `zeta > 1` very stable — log capped at `obu` + `5*log(zeta)+zeta-1`

The port implemented **regime 2 only**, and `stability_func1/2` clamp their
argument at 0, so above neutral the correction silently vanished and regimes 3/4
were unreachable. Now a single scalar, GPU-safe pair —
`mo_profile_denom_m` / `mo_profile_denom_h` in `src/types/friction_velocity.jl` —
carries all four branches for every caller. The kernel's pre-existing
`clamp(zeta, -100, ZETAMAX_LAKE)` before the profile evaluation was also removed:
that is Fortran's clamp on the zeta UPDATE (F90:544-551, still applied), not a
clamp inside `FrictionVelocity`, and it truncated exactly the convective range
the four-way branch exists to handle.

**Which regime does the reference exercise?** All 47 steps are in **regime 1**.
`zeta` stays in `[-100.000, -15.729]` and never rises above the `-zetat = -0.465`
transition (`frictionvel.zeta_patch` is now written by the lake kernel; the
harness prints the histogram). A 273 K lake under 255 K air is strongly
convective, so the one branch the port had was the one branch the reference never
needed. Regimes 3 and 4 are unreachable on this reference and are covered by the
unit tests in `test/test_lake_ref2m.jl` instead, which assert each branch is
reached and matches the independently-ported array form of the same Fortran
(`friction_velocity!`) to 1e-12.

**How much the regime-1 branch actually changes the answer — and where.** A first
version of that test asserted ">5% different from the old regime-2 expression" at
`zeta = -3.0` and FAILED at **0.52%**. The failure was correct and is worth
recording: regimes 1 and 2 JOIN CONTINUOUSLY at `-zeta*` by construction, so just
past the transition they agree to a fraction of a percent. Measured at the
reference's own condition (`zeta = -100`, `obu = -0.3 m`) the split is:

| profile | old regime-2-only | regime 1 | change |
|---|---|---|---|
| momentum `denom_m` | 5.96 | 8.78 | **+32%** |
| heat `denom_h` | 5.89 | 6.07 | +3.0% |

That asymmetry IS the parity result: `ustar` moved a long way (0.075 -> 0.058,
onto Fortran's ~0.061) while `temp1` barely moved, which is why C corrected the
momentum profile but left `rah` dominated by defect D below. Testing the branch
near its transition would have made this fix look cosmetic.

### D (NEW, open) — the lake never freezes, so it is on the WRONG roughness branch

C is faithful to CTSM and fixes the momentum profile measurably, but the run
metric got WORSE (3.11e-01 -> 4.97e-01, driven by LAKEICE 1.8e-02 -> 3.8e-01).
That is a compensating-error pair coming apart, and the second error is now
localised exactly.

Fortran does not write `ustar`/`obu`/`rah` to h0, but the flux definitions invert
from fields that ARE there (`LAKE_AERO_PROBE=1` does this):

```
rah = rho*cp*(TG - thm)/FSH      (F90:526)
ram = -rho*u/TAUX                (F90:699)
```

| step | zeta J | ustar J | z0mg J | z0hg J | rah J | rah F | ram F | TG J | TG F |
|---|---|---|---|---|---|---|---|---|---|
| 1 | -100.0 | 0.0534 | 3.39e-05 | 1.00e-04 | 316.8 | 159.7 | 372.9 | 273.25 | 269.50 |
| 8 | -100.0 | 0.0578 | 3.11e-05 | 9.20e-05 | 296.5 | 145.7 | 353.5 | 276.17 | 267.29 |
| 16 | -100.0 | 0.0559 | 3.19e-05 | 9.44e-05 | 305.4 | 148.7 | 357.7 | 275.23 | 268.29 |

Reading it:

* **Momentum is now right.** `ram_F = 353` implies `ustar_F = sqrt(um/ram) ~ 0.061`.
  The port is at 0.0578 after C, against ~0.075 before it — the regime-1 momentum
  cap plus convective term moved `ustar` onto the reference. This is C working.
* **Heat is not**, and it is not the regime. At `zeta = -100` with `z0h ~ 1e-4`
  the regime-1 and regime-2 heat denominators nearly coincide (6.84 vs 6.66), so
  C barely moved `temp1`. The reference needs `denom_h = 3.64`; the port has 6.84.
  Solving regime 1 for the roughness that closes that gap gives
  **`z0hg ~ 2.3e-03` against the port's 9.2e-05 — 25x too small.**
* **Why**: `TG_J = 273-276 K` while `TG_F = 267-269 K`. The port's lake surface
  never drops below freezing, so `lake_fluxes.jl` stays on the UNFROZEN roughness
  branch (`z0mg = max(minz0lake, cus*kva/ustar, cur*ustar^2/g) ~ 3e-05`) while
  CTSM is on the FROZEN branch (`z0mg = z0frzlake = 1e-03`, LakeCon.F90:50). It is
  self-reinforcing: too-small roughness -> `rah` 2x too high -> too little
  sensible heat lost -> surface stays warm -> never freezes.

So D is upstream of the surface-flux kernel: the lake surface energy/ice state,
not the Monin-Obukhov solve. Do not chase it inside `lake_fluxes.jl`'s profile
code again — that is now correct.

### Also found: `MINZ0LAKE` is a MIXED PAIR

`LakeCon.F90:143-158` sets `fcrit` and `minz0lake` together from one flag:

| `lake_use_old_fcrit_minz0` | `fcrit` | `minz0lake` |
|---|---|---|
| `.true.` (Subin 2011) | 22 | 1e-05 |
| `.false.` (**the default, and the reference's setting**) | 100 | **1e-10** |

`lake_fluxes.jl` has `FCRIT = 100.0` (the `.false.` value) and
`MINZ0LAKE = 1.0e-5` (the `.true.` value) — one from each branch, which cannot
both be right. Not on the critical path for D (the port's `z0hg ~ 9e-05` sits
above the 1e-05 floor, so it is currently inert), but it is a latent 1e5x floor
error that will bite the moment D is fixed and the roughnesses move. Left
unchanged here so it lands with its own evidence.

## STEP 4 — D is CLOSED: `ks` (the eddy depth-decay rate) was hardwired to zero

### The mechanism, with values

`LakeFluxesMod.F90:755` computes the depth-decay rate of the lake's wind-driven
eddy diffusivity, and `LakeTemperatureMod.F90:373` uses it:

```
ks = 6.6 * sqrt(|sin(lat)|) * u2m^-1.84                       [1/m], lat in RADIANS
ke = vkc*ws*z_lake/p0 * exp(-ks*z_lake) / (1 + 37*ri^2)       [m2/s]
```

`lake_fluxes.jl` wrote **`ks_col[c] = zero(T)`** — a literal zero, because the
gridcell latitude had never been wired into the lake kernel. The comment there
justified the literal on AD grounds (`sqrt(0)` has an `Inf` ForwardDiff
derivative), which is true at the equator and irrelevant everywhere else.

`exp(-0*z) == 1` at every depth. So the surface eddy diffusivity **never decayed
downward**: `ke` grew linearly with `z` all the way to the lake bed and the whole
50 m column stayed turbulently coupled to the skin. At Bow (lat 51.17 deg,
`u2m ~1.6 m/s`) the correct `ks` is **~2.4 /m**, so at 5 m depth the port's eddy
diffusivity was `exp(2.4*5) = 1.6e5` times too large.

The consequence is directly visible in `savedtke1 = kme(1)*cwat`, the top-layer
eddy conductivity the surface solve conducts through:

| | port (before) | port (after) | Fortran regime |
|---|---|---|---|
| `savedtke1`, steps 1-3 | 150-160 W/m/K | 131-133 | unfrozen eddy branch |
| `savedtke1`, steps 7-16 | 148-160 W/m/K | **0.61-0.67** | frozen branch (`km + fangkm`) |

With `tksur/dzsur = 150/0.05 = 3000 W/m2/K`, `t_grnd = ax/bx` is welded to
`t_lake(1)`, and `t_lake(1)` is welded to 50 m of 277 K water. Every W/m2 the
surface lost was replenished from the deep lake instead of from the top ~0.3 m,
so `t_lake(1)` fell only ~0.12 K/step against a `gnet` of -280 W/m2 — which over
0.1 m of water is 2.4 K/step. In the reference, layers 1+2 alone account for the
entire `gnet` (-207 and -34 W/m2 against `gnet = -243`): **the reference barely
mixes at all.**

The chain, end to end:

```
ks = 0  ->  eddy diffusivity never decays with depth
        ->  the whole 50 m column feeds the skin  (savedtke1 150 vs 0.6)
        ->  t_grnd welded to t_lake(1) ~276 K, never reaches 273.15 K
        ->  the UNFROZEN roughness branch stays selected (z0mg ~3.2e-5)
        ->  z0hg ~9.5e-5 instead of ~7.6e-4
        ->  rah ~305 s/m instead of ~145
        ->  too little sensible+latent heat lost  ->  surface stays warm.
```

It is self-reinforcing, which is why three previous attempts inside the
Monin-Obukhov solve could not move it.

### The fix

Wire `grc%lat` into the lake flux kernel (a new `grc_lat` field on `LakeFluxDV`,
passed from `clm_drv_core!`) and restore the Fortran expression, branching on
`|sin(lat)| > 0` so the equatorial `sqrt(0)` AD hazard the literal zero was
guarding against is still handled — and handled without deleting the term.

This is another instance of the master bug class in MEMORY
(`dead-initcold-systemic`, "ported then never called"): the formula was ported
correctly, its input never was, and the placeholder was documented as a
simplification rather than as a defect.

### Measured: the four numbers that decide whether D is closed

`LAKE_AERO_PROBE=1`, steps 1-16. `rah_F` is inverted from the reference's own
`FSH`/`TG`/`TBOT` (LakeFluxesMod.F90:526), so it compares like with like.

| step | z0mg J before | z0mg J after | rah J before | rah J after | rah F | TG J before | TG J after | TG F |
|---|---|---|---|---|---|---|---|---|
| 1 | 3.39e-05 | 3.39e-05 | 316.8 | 316.8 | 159.7 | 273.25 | 273.25 | 269.50 |
| 5 | 3.19e-05 | 3.19e-05 | 305.4 | 305.4 | 146.6 | 274.90 | 273.99 | 267.57 |
| 7 | 3.19e-05 | **1.000e-03** | 296.5 | **142.0** | 144.8 | 276.29 | **267.50** | 267.40 |
| 10 | 3.11e-05 | **1.000e-03** | 296.5 | **144.4** | 144.9 | 275.96 | **267.28** | 267.25 |
| 16 | 3.19e-05 | **1.000e-03** | 305.4 | **148.1** | 148.7 | 275.23 | **268.34** | 268.29 |

* the surface reaches freezing (`t_lake(1)` hits exactly 273.150 at step 6 and
  stays; before, `min(TG_J)` over 47 steps was 273.246 — it NEVER went below);
* it picks up `z0frzlake = 1.000e-03` exactly, 31x the unfrozen Charnock value;
* `rah` moves 305 -> ~145, onto the reference's 145-149 (**within 0.5%**);
* `TG` tracks the reference to **0.04 K** from step 7 on, against 8-9 K before.

Lake ice now forms: `LAKEICEFRAC_SURF` 0.000 at every one of 47 steps before,
0.277 at step 16 after, against the reference's 0.361.

### Per-field parity, before -> after (47 steps, max|rel|)

| field | #279 head | after D | note |
|---|---|---|---|
| LAKEICE / LAKEICE_SRF | 3.80e-01 | **8.50e-02** | 4.5x |
| FSA | 5.00e-01 | **3.80e-02** | dawn/dusk near-zero-solar artifact, moved off |
| EFLX_LH | 2.70e-01 | 2.40e-01 | max now in the freeze-onset window |
| FIRE | 1.40e-01 | 1.30e-01 | |
| TG | 3.30e-02 | 3.10e-02 | |
| TLAKE | 1.30e-02 | 1.20e-02 | |
| TSA | 7.40e-03 | 7.40e-03 | unchanged (finding B) |
| H2OSNO | 2.00e-01 | 2.00e-01 | unchanged |
| **FSH** | 2.70e-01 | **3.30e-01** | **WORSE — see below** |
| **run max\|rel\|** | 4.966e-01 | **3.308e-01** | |

**FSH gets worse and that is reported, not tuned.** The FSH maximum moves INTO
the freeze-onset window and grows there, while collapsing everywhere else:

| step | 1 | 4 | 5 | 6 | 7 | 8 | 12 | 20 | 30 | 40 |
|---|---|---|---|---|---|---|---|---|---|---|
| FSH rel before | 0.220 | 0.160 | 0.180 | 0.190 | 0.190 | 0.180 | 0.230 | 0.210 | 0.130 | 0.220 |
| FSH rel after | 0.220 | 0.230 | 0.290 | 0.330 | **0.026** | **0.006** | **0.007** | **0.032** | **0.009** | **0.009** |

From step 7 on FSH improves by **20-30x**. The surviving maximum is a
**cold-start phase lag**, not an aerodynamic bias: the port cold-starts from a
uniform 277 K lake (exactly CTSM's `TemperatureType.F90:833-836`) while the
reference's first h0 record already shows a partly-cooled top layer
(`TLAKE(1) = 274.81`, i.e. the reference has already shed ~255 W/m2 for one hour
before the harness's step 1). The port therefore reaches freezing at step 6
where the reference is frozen from its first record, and across steps 4-6 the
two models are in DIFFERENT SURFACE PHASES — one on ice, one on open water.
That is a residual worth its own investigation (it is a harness/reference
alignment question about what the `nstep=0` record contains, not a kernel), and
it is deliberately left open rather than absorbed into a looser tolerance.

### Also fixed: the `MINZ0LAKE` / `FCRIT` mixed pair

`MINZ0LAKE` is now `1.0e-10`, matching `FCRIT = 100.0` on the
`lake_use_old_fcrit_minz0 = .false.` side of `LakeCon.F90:143-158` (CTSM's
default, `LakeCon.F90:98`, and the reference's setting). This is **coupled to the
D fix and had to land with it**: while the lake never froze, `z0hg ~9.5e-05` sat
above the old `1e-05` floor and the mismatch was inert. Once D lets the surface
freeze, the roughnesses move, and a floor 1e5x above CTSM's is a live error.

## Status

STEP 1 complete (table above; `z0param_method` retracted as a dead end).
STEP 2 complete. STEP 3 complete (B and C fixed; C exposed D). The "surface
turbulent-flux residual" is **four** distinct defects,
none of which is the Monin-Obukhov *iteration* that three previous fixes chased:

- **A (harness, dominant):** h0 record off-by-one — 1.30e+01 → 3.11e-01 max|rel|. Fixed.
- **B (port, dead write): FIXED.** `t_ref2m`/`q_ref2m`/`rh_ref2m` are now written
  on the lake path (LakeFluxesMod.F90:707-717). TSA 1.10e-01 -> 8.30e-03; every
  other field byte-identical.
- **C (port, physics): FIXED.** All four `zeta` stability regimes ported as one
  scalar pair shared with the array form. The reference sits in regime 1 at every
  one of its 47 steps. FSH 3.10e-01 -> 2.70e-01, EFLX_LH 2.80e-01 -> 2.70e-01,
  and `ustar` moves onto the reference (0.075 -> 0.058 vs Fortran's ~0.061).
- **D (port, physics): FIXED.** `ks`, the depth-decay rate of the lake eddy
  diffusivity (`LakeFluxesMod.F90:755`), was hardwired to zero because the
  gridcell latitude had never been wired into the lake kernel — so `exp(-ks*z)`
  was 1 at every depth and the whole 50 m column mixed to the skin. The surface
  therefore could not reach freezing and stayed on the unfrozen roughness
  branch. Fixed: lake now freezes at step 6, picks up `z0frzlake = 1e-3`,
  `rah` 305 -> 145 (reference 145-149), `TG` within 0.04 K, `LAKEICE`
  3.80e-01 -> 8.50e-02, run max|rel| 4.97e-01 -> 3.31e-01. FSH's maximum moves
  into (and grows in) the freeze-onset window while improving 20-30x elsewhere;
  see STEP 4.
- **E (open, NEW):** the port and the reference are not aligned at step 1. The
  port cold-starts at a uniform 277 K lake (CTSM's own cold start) while the
  reference's first h0 record already shows `TLAKE(1) = 274.81`, one hour of
  cooling ahead. The port reaches freezing ~6 steps late as a result, and that
  lag is now the largest single contributor to the remaining metric.

B and C are both lake-landunit-gated by construction (`lake_fluxes.jl` runs only
over `mask_lakep`; nothing wrote `t_ref2m_patch` on a lake patch before), so
neither can perturb any non-lake result.
