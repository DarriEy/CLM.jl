# Fortran reference-generation audit

Date: 2026-08-25

## Result

The exact stock CTSM source is identifiable and fetchable, but the current repository does
**not yet provide a clean-room, end-to-end generator for all Fortran results intended for
the paper**. Existing material is valuable reconstruction evidence, not yet a submission-
grade reference pipeline.

## Verified source

- Official remote: `https://github.com/ESCOMP/CTSM.git`
- Annotated tag: `ctsm5.3.012` (`ae99861566e6aa8fdbdfff0f2bbe3fdad60ba7bb`)
- Peeled commit: `ab466d6f9789ca3df2c72bda46cf7afed2d04102`
- A shallow clean clone at that tag resolves to the expected commit and exact tag.
- External repositories are pinned by the CTSM checkout but are not populated by a normal
  shallow clone; the definitive source manifest must include their hashes.

## Available instrumentation

The repository contains:

- `scripts/validation/fortran_pdump/pdumpMod.F90`;
- `clm_driver.hooks.diff` and additional process-specific diffs;
- BGC dump modules and hooks;
- FATES dump modules and an instrumentation helper;
- several reference-generation scripts and extensive methodological notes.

This is enough to reconstruct important parts of the oracle, but the relationship between
each historic reference file, its exact SourceMods set, case configuration, executable,
input files, and source commit is not represented by one manifest.

## Blocking findings

1. `regen_parity_references.sh` describes itself as documentation-as-code, contains
   historical absolute paths and macOS-specific commands, and relies on a previously built
   executable that is no longer present in the staged data tree.
2. `gen_fortran_refs.sh` requires an already instrumented case and references a “banked
   recipe” and machine state outside the repository. It is not a clean-source build recipe.
3. The prior `installs/clm` data location contains case directories but no CTSM source or
   built executable in the inspected tree.
4. Several reference modes use different instrumentation files or additional history-field
   registrations. There is no per-reference instrumentation checksum ledger.
5. Reference inputs and restarts live in an external archive whose file-level licences,
   checksums, and redistribution status have not yet been audited.
6. The existing recipes document flaky startup retries and old Python/compiler constraints;
   these must be re-evaluated on the release platform rather than normalized as expected
   scientific workflow behaviour.

## Required remediation

Build one parameterized reference pipeline that starts from the verified CTSM tag and:

1. checks out all pinned externals and records their hashes;
2. applies a named, checksummed instrumentation patch set;
3. proves instrumentation state-neutrality on a non-dump control run;
4. creates/builds a fully specified case from committed configuration;
5. acquires inputs from cited locations and verifies checksums;
6. runs a declared window with deterministic success/failure handling;
7. writes reference NetCDF with source/configuration metadata;
8. produces a run record and artifact manifest;
9. can be executed from a clean Linux environment by a second person.

Historic references may be used for debugging, but definitive manuscript claims should use
newly generated references from this pipeline or provide a documented equivalence audit
that is equally strong.

## Clean-source execution update — 2026-08-26

- `./bin/git-fleximod update` completed successfully in the clean tagged checkout.
- The recursive source hashes are frozen in
  `repro/manifests/ctsm5.3.012-externals.txt`.
- The Bow BTRAN oracle specification and checksummed staged inputs are frozen in
  `repro/configs/bow_btran_oracle.toml`.
- CIME fails under the system Python because `distutils` is unavailable. Running under an
  explicit uv-managed Python 3.11 creates the case successfully.
- A clean `I2000Clm50SpRs`, `f09_g17`, one-task stock case now builds successfully in an
  isolated MPICH environment. The executable SHA-256 is
  `9719074f2300508c722bf6326c420f8e19339b737ce28d5c268b81094db0a57d`.
- The frozen toolchain is Python 3.11.16, CMake 3.31.8, ESMF 8.9.1, MPICH 4.3.2,
  netCDF-C 4.9.3, and netCDF-Fortran 4.6.2. Required Darwin build-only settings are
  recorded in `repro/patches/ctsm5.3.012-darwin-build-compat.patch`; no CLM scientific
  source was changed for the stock build.
- The staged archive has a January 2003 Bow restart and complete run namelists, but not the
  April 29 shared-state dump cited by historical notes. The new pipeline must run the
  January trajectory through the target window, then emit the shared-state oracle.
- A relocated `startup + finidat` replay from the January restart completed for six model
  steps with the stock executable. The same replay with the established `pdumpMod` and two
  driver hooks also completed and exercised all snapshot writes. Every selected history
  variable and value was identical; only the netCDF root name and creation timestamp
  differed. Checksums and executable identity are frozen in the oracle specification.

No BTRAN exception has been added. The annual consumers support the daily-minimum
amplification hypothesis, but only the new identical-state oracle can close the claim.
Base instrumentation state-neutrality, toolchain provisioning, and the target trajectory
are closed; first-divergence localization remains active.

The target trajectory has now also completed from 2003-01-01 through the full
2003-04-29--2003-05-07 oracle window. It emitted paired prognostic and BTRAN/LUNA
sidecars at 193 hourly steps. The lowest vegetated-patch value is BTRAN
`0.05617311246487857` at 2003-05-01 21:00 (step 2901). Re-injecting the registered
17-field initial state into CLM.jl gives zero initial relative error, but the next Julia
step produces BTRAN `1.0` for both vegetated patches versus CTSM values `0.0561731`
and `0.127648`. This is a reproduced discrepancy, not grounds for an exception: the
next localization must freeze the remaining instantaneous forcing and PHS conductance
inputs before assigning the first divergent routine.

The first expansion is complete. At step 2901, Julia matches CTSM's post-call
`K_SOIL_ROOT`, root conductance, and soil conductance exactly (zero max absolute and
relative error), excluding conductance construction as the source of the BTRAN gap.
CTSM's final sunlit/shaded stress factors for the tree patch are `0.0492946` and
`0.0598674`; its weighted BTRAN is therefore `0.0561731`. The Julia one-step harness
initially returned BTRAN `1.0` and did not reproduce CTSM's sunlit/shaded LAI because
the replay omitted canopy-layer albedo state carried from the preceding step.

Adding `NRAD`, `TLAI_Z`, `FSUN_Z`, and the direct/diffuse sunlit/shaded absorbed-light
profiles makes `LAISUN` and `LAISHA` exact. The same replay now gives Julia BTRAN
`0.7974531` and `0.4016276` versus CTSM `0.0561731` and `0.1276478`; all three
soil/root conductance arrays remain exact. Read-only call traces place the next
divergence in the coupled photosynthesis/hydraulic demand: the final sunlit stomatal
conductance is about `7.40e5` versus CTSM `6.61e5` for the tree patch, and `8.38e5`
versus `1.72e6` for the grass patch. Conductance construction and sun/shade geometry
are therefore excluded, but the photosynthesis-demand/calcstress boundary is not yet
closed. This remains an open localization result, not a scientific exception.

The Julia trace is now reproducible rather than console-only. For the last recorded
coupled call it reconstructs unstressed sunlit/shaded transpiration demands of
`3.21064e-7` / `5.95820e-7 kg m-2 s-1` for the tree and `3.46312e-6` /
`6.11948e-7 kg m-2 s-1` for the grass, alongside unstressed stomatal conductance,
`BSUN`/`BSHA`, and all four plant water potentials. The previous CTSM call trace was
not retained as an artifact, so its console transcription is not being promoted to
submission evidence. The next oracle run must persist that trace and compare the call
sequence mechanically; only then can demand generation be separated from the hydraulic
Newton solve.

The CTSM trace is now persisted and checksummed. It contains 33 tree and 21 grass
`calcstress` calls, whereas Julia's current trace contains 17 pass-level records per
patch; the two files therefore must not be joined by line number. At the earliest
recorded outer state, `qsatl`, LAI/SAI, dry fraction, density, pressure, and reference
temperature agree to displayed precision; `qaf` differs by only `1.36e-7`, and boundary
conductance is exact for grass and within 0.36% for the tree. Despite that, the first
recorded unstressed shaded stomatal conductance is `8.05e5` in Julia versus `3.02e5` in
CTSM for the tree, and `2.07e6` versus `7.95e5` for grass. This moves the first-divergence
search upstream of `getqflx` and the hydraulic Newton solve into the sun/shade
photosynthesis/CI coupling. An inner-call Julia capture is still required before naming
the exact expression, because the current traces have different call granularity.

A matching CTSM trace at the `hybrid_PHS` return boundary removes the inner-call
granularity ambiguity but reveals a second control-flow difference: CTSM completes 11
tree and 7 grass outer calls during the canopy-temperature solve, while Julia completes
17 for each. The first outer result is already close but unequal. Tree `AN_SUN` is
`0.73845` in CTSM versus `0.60670` in Julia and grass is `3.82239` versus `3.72144`;
their corresponding stress and stomatal-conductance values differ in the same direction.
The trajectories then separate sharply. This evidence supports amplification through
the coupled convergence path, but does not yet justify changing its iteration limits.
A historical replay immediately before the fused PHS kernel was attempted and rejected
as a comparator because that revision predates the required LUNA and current driver
configuration. The next defensible boundary is therefore the first-call photosynthesis
input/output tuple on the current code, not a cross-version result.

That tuple exposed two additional carried fields rather than a missing science branch.
CTSM scales shaded LUNA `VCMAX`, `JMAX`, and `TPU` by
`VCMAXCINTSHA/VCMAXCINTSUN`; Julia already contains the same branch, but the replay had
not carried those integration coefficients. Adding them restores the shaded capacities
exactly. The first tree `hybrid_PHS` result consequently moves from `AN_SUN=0.60670`
to `0.73843694`, versus CTSM `0.73844529`; stressed sunlit conductance is `48985.65`
versus `48986.83`. The initial leaf photosynthesis/PHS tuple is therefore close. The
remaining large BTRAN difference develops later as the canopy energy solve follows a
longer convergence trajectory. The complete replay now yields BTRAN `0.8455352` and
`0.4414279` versus CTSM `0.0561731` and `0.1276478`; first-call agreement is not
whole-step parity.

End-of-iteration instrumentation in both canopy solvers then located the first large
energy-path discrepancy before convergence can amplify it. At iteration one, CTSM
enters the tree and grass energy balances with `SABV=772.8918277` and
`190.2960005 W m-2`; Julia has zero for both. The resulting bounded Newton increment
is `+1 K` in CTSM and `-1 K` in Julia. The remaining conductance and photosynthesis
quantities at that boundary are comparatively close, so convergence settings are not
the first cause. `SABV` is now part of the diagnostic oracle and replay registry. The
current replay hook runs before `surface_radiation!`, however, and that routine
recomputes the field before `canopy_fluxes_core!`. The next isolation step is a
post-surface-radiation/pre-canopy injection hook; changing iteration limits or adding a
scientific exception would be premature.

The boundary replay confirms that diagnosis. Applying the CTSM `SABV` values only at
entry to `canopy_fluxes_core!` changes tree BTRAN from `0.8455352` to `0.05617330`
against CTSM's `0.05617311`, and grass from `0.4414279` to `0.12763999` against
`0.12764776`. Both solvers now emit 18 patch-iteration records (11 tree and 7 grass),
instead of Julia retaining both patches for 17 iterations. At the first iteration the
largest remaining relative differences in the traced flux/conductance quantities are
approximately `2e-5`–`7e-5`; final absolute BTRAN errors are `1.86e-7` and `7.77e-6`.
This is causal isolation, not a production fix: the override is explicitly gated and
off by default. Inspection of `surface_radiation!` further narrows the replay omission:
it constructs `SABV` from the bulk two-band `FABD`/`FABI` canopy absorptances, while the
oracle currently carries only the canopy-layer visible-band sun/shade arrays used by
photosynthesis. The replay's cold-started bulk arrays are therefore zero. Those bulk
absorptances must be added to the next oracle revision and injected before
`surface_radiation!`; this evidence does not indicate a production radiation defect.

That upstream test is now complete. The oracle carries both `FABD` and `FABI`, and
injecting them before the ordinary `surface_radiation!` call reproduces the same
`SABV`, iteration counts, and near-parity BTRAN obtained in the causal override test.
The temporary direct-`SABV` hook has therefore been removed from production source.
The remaining discrepancy is small but still above the frozen strict threshold, so the
comparison remains formally failing and no exception is claimed.

The first remaining photosynthesis residual is now localized to canopy humidity and
boundary conductance. Julia enters the tree call with `QAF=0.0065372678704` and
`gbmol=802654.6136`, versus CTSM `0.0065374042231` and `802655.6602`. A direct
`QAF` carry experiment leaves the result unchanged because `CanopyFluxes` overwrites
that state before photosynthesis using `(forc_q + qg)/2` in both codes. Consequently,
`QAF` was not retained as an oracle field. The next comparison must trace `qg` and its
ground saturation inputs at canopy entry, plus the independent `rb`-to-`gbmol` path.

An expanded end-of-iteration trace separates those paths. Atmospheric `forc_q` and
leaf `qsatl` are exact, while column `qg` is lower in Julia by `2.72705e-7`; this
fully accounts for the initial `QAF` direction because `QAF=(forc_q+qg)/2`.
Separately, Julia's first tree `uaf` is lower by `5.88790e-7 m s-1` and `rb` is
higher by `5.48786e-5 s m-1`; the grass differences have the same sign. The residual
therefore has two concrete upstream branches: `SurfaceHumidity` (`qg_snow`, `qg_soil`,
surface fractions and their temperatures) and the canopy `ram1 -> uaf -> rb -> gbmol`
chain. Neither is evidence for changing the hydraulic or convergence algorithms.

Tracing the surface mixture shows `qg_snow`, `qg_soil`, `qg_h2osfc`, and
`frac_h2osfc` are identical. The full `2.72705e-7` humidity difference is produced by
`frac_sno_eff`: Julia `0.2749586811`, CTSM `0.2749147620`. The humidity investigation
therefore moves upstream to the snow-cover-fraction update. For aerodynamic resistance,
`ram1=um/ustar^2` and `uaf=um*sqrt(1/(ram1*um))`; the traced pre-iteration `um` agrees,
so `uaf` reduces to the already-different incoming `ustar`. The next aerodynamic
boundary is the initial Monin–Obukhov solve and its roughness/displacement inputs.
