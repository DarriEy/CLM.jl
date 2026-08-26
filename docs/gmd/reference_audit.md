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
The active task is now adding the BTRAN/LUNA diagnostic fields and executing the target
trajectory; base instrumentation state-neutrality and toolchain provisioning are closed.
