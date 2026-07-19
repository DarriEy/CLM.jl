# `use_luna` default flip — campaign record (row M4)

Campaign #4 of the defaults programme. Precedent: `docs/PHS_DEFAULT_BLOCKERS.md`
(row M3, measured and REVERTED), `use_bedrock` (#252), `use_aquifer_layer` (#259).

## What CTSM says — namelist default vs code fallback

The two disagree, which is the systemic root cause this whole programme exists
to remove (14 such disagreements found in the audit; CLM.jl copied the code
fallback in every one).

| Source | Value | File:line |
|---|---|---|
| **Namelist default** (authority) | **`.true.`** | `bld/namelist_files/namelist_defaults_ctsm.xml:578` |
| — under `use_fates=".true."` | `.false.` | `namelist_defaults_ctsm.xml:579` |
| — under `phys="clm4_5"`, non-FATES | `.false.` | `namelist_defaults_ctsm.xml:580` |
| Namelist *definition* entry fallback | `.false.` | `namelist_definition_ctsm.xml:874-875` |
| **Code fallback** (NOT authority) | `.false.` | `clm_varctl.F90:371` |

So for `phys=clm5_0`, non-FATES — the configuration this port targets — CTSM
defaults **`use_luna = .true.`**, and the default is **CONDITIONAL on
`use_fates`**, exactly like `use_bedrock`. CLM.jl shipped `false`, i.e. the
Fortran code fallback.

CTSM enforces the condition rather than merely defaulting it:
`controlMod.F90:505-508` `endrun`s with *"luna is not compatible with FATES"*.
CLM.jl already mirrors this at `src/infrastructure/control.jl:104`, so a bare
`true` default would make every FATES run throw — the conditional is required,
not cosmetic.

## Baseline (this campaign's own measurement)

Full suite, Julia 1.12, `--check-bounds=yes`, at `f6ac7f5`:

**26351 passed / 0 failed / 0 errored / 3 broken** (26354 total, 40m54s).

Recorded before any change. Note this is NOT the 26150 of the `use_hydrstress`
campaign — the baseline moved under #262/#263/#264, which is why older numbers
must not be trusted.

## Blast radius

- 72 explicit `use_luna=true` call sites, 42 explicit `false`.
- **32 test files** call `clm_initialize!` / `clm_run!` / `CLMDriverConfig(`
  without passing `use_luna`, i.e. they inherit the default and are exposed to
  the flip.

## Result of the flip

`use_luna` became `Union{Bool,Nothing}`, `nothing` → `!use_fates`, in
`clm_initialize!`, `clm_run!` and the `CLMDriverConfig` constructor — the same
shape as `use_bedrock` (#252).

| | baseline | `use_luna` conditional-on |
|---|---|---|
| passed | 26351 | **26351** |
| failed | 0 | **0** |
| errored | 0 | **0** |
| broken | 3 | **3** |

**Zero movers.** Per-testset tallies are identical; the only diffs between the
two logs are testset *ordering* and wall-clock times. No test needed
recalibration and nothing threw. The conditional resolves correctly:
`CLMDriverConfig().use_luna == true`, `…(use_fates=true).use_luna == false`.

## Why zero movers is NOT a clean bill of health

A green suite here would be a *vacuous* green, so it was checked rather than
trusted. Two independent reasons the flip cannot move tests:

1. **LUNA only fires at `is_end_curr_day`** (`clm_driver.jl:2923`). Most unit
   tests run a handful of steps and never cross a day boundary, so LUNA is a
   no-op in them regardless of the flag. The harnesses that *do* run for days
   already pass `use_luna` explicitly.

2. **The LUNA→vcmax feedback is missing from the default photosynthesis path.**
   This is a real port gap, found by this campaign.

## THE FINDING — LUNA is computed and then discarded (non-PHS path)

Measured at Bow (`scripts`-free probe, 4 patches, hourly, 16 days, so the
240-hour means fill), default config vs explicit `use_luna=false`:

| quantity | LUNA on | LUNA off |
|---|---|---|
| `par240d_z` allocated | `(4, 1)` | `(0, 0)` |
| `vcmx25_z` | `[30.0, 38.14, 43.92, 30.0]` | *(unallocated)* |
| `jmx25_z` | `[60.0, 76.27, 87.84, 60.0]` | *(unallocated)* |
| `vcmax_z` (the capacity actually used) | `0.35120209710203465` | `0.35120209710203465` |
| `psnsun` | `0.12483706075381139` | `0.12483706075381139` |
| `fpsn` | max abs diff **0.0** | — |

LUNA runs and genuinely acclimates (30.0 → 38.1/43.9 on the two C3 patches,
which carry `elai` 4.03 and 1.21). **And not one digit of GPP changes.**

Direct instrumentation of the branch confirms it: a temporary print inside
`if use_luna && c3flag_p` (`photosynthesis.jl:2132`) fired **0 times** across
the whole 16-day run.

The cause is structural. That branch — and the leaf-respiration branch at
`photosynthesis.jl:2100` — live inside `_psn_phs_pass2_body!`, the **PHS**
(plant-hydraulic-stress) Pass-2 kernel, reached only via
`photosynthesis_phs!` when `use_hydrstress=true`. CTSM has the LUNA branch in
**both** routines:

| CTSM routine | starts | LUNA branches | ported to CLM.jl? |
|---|---|---|---|
| `Photosynthesis` (non-PHS) | `PhotosynthesisMod.F90:1220` | **1721** (lmr), **1748** (vcmax25/jmax25) | **NO** |
| `PhotosynthesisHydraulicStress` | `PhotosynthesisMod.F90:2694` | 3335 (lmr), 3377 (vcmax25/jmax25) | yes → `photosynthesis.jl:2100, 2132` |

The default non-PHS `photosynthesis!` (`photosynthesis.jl:2305`) declares
`use_luna::Bool=false` in its signature and **never references it again** —
an accepted-and-dropped kwarg. This is the *"ported then never called"* bug
class (`dead-initcold-systemic`), in its subtler form: the flag is threaded all
the way from `CLMDriverConfig` through `canopy_fluxes` into `photosynthesis!`,
so every wiring assertion along the chain passes; only the *body* is missing.

Note the coupling with row M3: because CLM.jl's only LUNA consumption sits in
the PHS kernel, and the `use_hydrstress` default flip was measured and reverted
(`docs/PHS_DEFAULT_BLOCKERS.md`), LUNA has **no effect on any default run** —
neither flag alone is sufficient.

## Verdict

`use_luna` **can** default on, and now does:

- it is CTSM's verified namelist default for `clm5_0` non-FATES;
- it is correctly conditional, so FATES runs are unaffected and CLM.jl's
  existing LUNA+FATES guard (`control.jl:104`) never trips;
- it costs zero test failures and zero errors.

But **do not read this row's closure as "the port now has CLM5 photosynthesis."**
It does not. Until `PhotosynthesisMod.F90:1721` and `:1748` are ported into the
non-PHS path, enabling LUNA buys the acclimation diagnostic (`VCMX25LUNA`,
`history_writer.jl:610`) and the per-step accumulator cost, with no effect on
GPP, transpiration or leaf respiration.

A reviewer who would rather not carry an inert default may prefer to hold this
flip until the consumption branches land. Both options are defensible; the
measurements above are the input to that decision, and the flip is a
three-line revert either way.

## Order of work to make the default meaningful

1. Port `PhotosynthesisMod.F90:1748` (vcmax25/jmax25 from `vcmx25_z`/`jmx25_z`)
   and `:1721` (LUNA leaf-respiration base rate) into the non-PHS
   `photosynthesis!` body, mirroring `photosynthesis.jl:2132` / `:2100`.
2. Add a value-asserting regression that LUNA on vs off **changes GPP** on a
   multi-day run — the check this campaign had to perform by hand, and whose
   absence is why the gap survived. Assert the VALUE, not the wiring.
3. Re-run the M3 (`use_hydrstress`) experiment once its own blockers are fixed;
   M3 and M4 together are what actually constitute CLM5 photosynthesis.
