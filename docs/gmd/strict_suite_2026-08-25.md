# Strict external-data suite qualification — 2026-08-25

Status: **FAIL**.

## Execution

- CLM.jl base: `e5c3c183ebb318bdaaa89cf16fa1f458a39fbd20`
- Submission branch commit at execution: `cf7abdd`
- Julia: 1.12.6
- Platform: macOS 26.5, Apple M1
- Command:

```sh
SYMFLUENCE_DATA=/path/to/SYMFLUENCE_data \
CLM_REQUIRE_TESTDATA=1 \
julia --project=. --check-bounds=yes -e 'using Test; include("test/runtests.jl")'
```

Julia reported **27,766 passed, 1 failed, 3 broken (27,770 total)** in 49 min
06.9 s and exited with status 1. The strict release gate therefore did not pass.

## Hard failure

The failure was under `Hillslope FULL-DRIVER multi-column timestep (E2E)`. Its asset-gated
path requires SNICAR optics and aging files that resolved to a developer-specific
`CESM_INPUTDATA` location rather than the staged qualification archive. A focused invocation
confirmed that strict mode converts these missing files into a failing assertion. This is
an input-discovery/reproducibility failure; it must not be relabelled as a successful
hillslope result.

## Broken tests despite strict mode

Three tests still recorded `Broken` rather than failing when their inputs were unavailable:

- FUN cold-start finiteness;
- lake water-balance closure;
- standalone `run_clm` (expects the Bow `clmforc.2002_2004.nc` fixture).

These custom gates do not all route through the shared `testdata_missing` helper. Thus
`CLM_REQUIRE_TESTDATA=1` is not yet a comprehensive enforcement mechanism. Every external-
data skip must use one gate before the definitive clean-room run.

## Additional qualification observations

- The primary Fortran per-step, CN/BGC, LUNA/decomposition, multisite, glacier, lake,
  snow, urban, and FATES real-parameter groups that executed all passed.
- The AD multi-scenario group passed 90/90, but the winter finite-difference comparison
  produced non-finite values while the scenario was labelled passed. This requires a
  tighter acceptance definition before supporting a derivative-verification claim.
- The synthetic parameter-recovery experiment labelled a case passed after returning
  Medlyn 6 versus truth 8 and Vcmax multiplier 1 versus truth 1.3, with a zero objective
  after one iteration. This is evidence of an under-constrained test or acceptance rule,
  not evidence that both parameters are identifiable and recoverable.
- The full initialization diagnostic reported many NaNs on nominally active fields while
  its selected-field assertions passed. The paper must specify which state subset is
  required finite rather than describing the whole initialized state as finite.

## Required closure

1. Make all required data resolve from a single documented archive and checksum manifest.
2. Route every data-gated test through the strict helper; demonstrate zero broken tests.
3. Re-run the entire suite on Linux in the locked CPU environment and require exit 0.
4. Prespecify derivative tolerances and make non-finite finite differences fail or be an
   explicitly justified exclusion.
5. Redesign parameter recovery for identifiability and assert parameter error, not merely
   optimizer termination or objective value.
