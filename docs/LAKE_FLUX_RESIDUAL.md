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

## Status

STEP 1 complete (table above; `z0param_method` retracted as a dead end).
STEP 2 complete: the residual is two harness/port defects, not a mis-ported
Monin-Obukhov kernel. Instrumenting the MO iteration was NOT needed and would have
chased a phantom — which is what the record-1 comparison had been manufacturing.
