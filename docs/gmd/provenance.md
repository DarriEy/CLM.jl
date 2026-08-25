# GMD release provenance

## Release workspace

- Created: 2026-08-25
- Repository: `https://github.com/DarriEy/CLM.jl.git`
- Branch: `gmd/submission`
- Base commit: `e5c3c183ebb318bdaaa89cf16fa1f458a39fbd20`
- Base commit date: 2026-08-05 04:53:59 -0600
- Base subject: `Merge pull request #329: harden the last unguarded in-kernel Thomas sweep`
- Base commit count: 1292
- Julia used for initial qualification: 1.12.6
- Workspace: a fresh clone under `/private/tmp`, outside Google Drive

The previous cloud-synced working copy was not used as a release base. Its HEAD was
`3f5a7063552f2a9242eeeaf4cf6d04b49cc7901e` (2026-07-24), 22 commits behind the current
public base. Its `paper/` directory and historical draft were absent from the public
repository and are intentionally excluded from the new manuscript baseline.

## Environment policy

The package supports Julia 1.10 and later. Development leaves the root `Manifest.toml`
untracked so CI resolves supported dependency sets. The paper artifact will additionally
provide a locked CPU environment for exact reproduction. Supported-package CI and exact-
artifact reproduction serve different purposes and both are required.

Initial instantiation on Julia 1.12.6 resolved successfully. The generated root
`Manifest.toml` is ignored and did not modify tracked files. Its resolved dependency set
must not yet be treated as the final paper lockfile.

## Source-model provenance: partially resolved gate

Repository documentation consistently identifies the reference source as:

- Repository: `https://github.com/ESCOMP/CTSM.git`
- Tag: `ctsm5.3.012`
- Annotated tag object: `ae99861566e6aa8fdbdfff0f2bbe3fdad60ba7bb`
- Commit: `ab466d6f9789ca3df2c72bda46cf7afed2d04102`

This identity is recorded in `docs/FORTRAN_VALIDATION_BACKLOG.md`; the upstream CTSM change
summary identifies `ctsm5.3.012` as the 2024-11-13 release point. It was independently
verified against the official remote with `git ls-remote` on 2026-08-25. The definitive
campaign must archive that verification and distinguish stock source from every
instrumentation patch.

Before experiments are frozen, record:

- upstream licence and citation;
- compiler and build flags;
- case XML/namelist settings and enabled options;
- changes used solely for reference dumps;
- proof that dump instrumentation does not change model state;
- source and checksums of restart, surface, parameter, aerosol, and forcing data.

The current developer documentation points to a machine-specific Fortran installation,
which is insufficient for publication.

## Licensing inventory

| Material | Current indication | Required action |
|---|---|---|
| Original Julia contributions | MIT (`LICENSE`) | Confirm all release files are covered |
| CTSM-derived implementation | BSD 3-Clause attribution in `NOTICE` | Legal/provenance review against exact upstream source |
| Fortran reference source | Not redistributed in this repository | Cite/archive lawful patch and acquisition recipe |
| CTSM input and parameter files | Mixed/unknown | File-level redistribution audit |
| Site forcing and observations | Mixed/unknown | File-level licence, citation, and redistribution audit |
| Generated model output | Author-generated | Confirm upstream/input restrictions do not attach |
| Plotting/analysis code | Repository code | Include in software or experiment archive |

## Validation ledger

The initial clean-clone full CPU test used this command:

```sh
julia --project=. --check-bounds=yes -e 'using Test; include("test/runtests.jl")'
```

Result on Julia 1.12.6 / macOS 26.5 / Apple M1:

- exit status: 0;
- elapsed test time reported by Julia: 12 min 14.2 s;
- 27,503 passed, 0 failed, 16 broken, 27,519 total;
- all 16 broken results were external-data gates because `SYMFLUENCE_DATA` was not set.

This is the clean-source/no-data baseline, not strict release qualification. The second run
with `SYMFLUENCE_DATA` set and `CLM_REQUIRE_TESTDATA=1` completed with 27,766 passes, one
failure, three broken tests, and exit status 1 after 49 min 06.9 s. Strict mode did not make
every missing fixture fatal. The complete interpretation and closure requirements are in
`strict_suite_2026-08-25.md`. The generated Julia 1.12.6 Manifest used for this qualification has SHA-256
`94ea12b7e1a720df91489ef969b3a9ef6f1ecd485a1563998043618726d258b3`.

The exact base commit previously completed GitHub Actions run
[`30999270755`](https://github.com/DarriEy/CLM.jl/actions/runs/30999270755) successfully
on 2026-08-05. Its four jobs were Julia 1.10 on Ubuntu, latest stable Julia on Ubuntu,
two-rank MPI bit-identity, and a default-bounds conservation gate. This establishes a
remote baseline for the source commit, but it does not replace the frozen paper campaign.
