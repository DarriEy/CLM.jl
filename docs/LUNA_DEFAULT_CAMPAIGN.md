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

Results of the flip are recorded below.
