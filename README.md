# CLM.jl

[![CI](https://github.com/DarriEy/CLM.jl/actions/workflows/test.yml/badge.svg)](https://github.com/DarriEy/CLM.jl/actions/workflows/test.yml)
[![Julia 1.10+](https://img.shields.io/badge/Julia-1.10%2B-9558B2.svg)](https://julialang.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

CLM.jl is a differentiable Julia port of the Community Land Model in
[CTSM](https://github.com/ESCOMP/CTSM). It implements the CLM5 land-surface
physics, hydrology, snow, carbon–nitrogen biogeochemistry, and a broad set of
optional process models in a structure-of-arrays design that supports automatic
differentiation and accelerator execution.

The project is currently an **unreleased research artifact** (`v0.1.0`), not a
drop-in replacement for CTSM and not yet registered in Julia's General registry.
Use it for research and development only after independently validating the
configuration and outputs that matter to you.

> [!IMPORTANT]
> CLM.jl was produced primarily by AI coding agents, with human supervision at
> the goal and validation level. The test and parity results below are real, but
> they do not substitute for independent scientific review. Expect latent defects
> outside the validated surface and do not use the model for consequential
> scientific, operational, or policy decisions without additional verification.

## Project status

Status at commit `32fd5ef` (2026-08-03):

| Layer | Result | Scope |
|---|---:|---|
| GitHub CI | 3/3 jobs passed | Julia 1.10, latest stable Julia, and a two-rank MPI bit-identity smoke test |
| Strict-data CPU suite | 27,537 passed, 0 failed, 3 hardware-skipped | 27,540 runtime checks; all external data gates enabled |
| Annual domain runs | Passed | Aripuanã: 8,784/8,784 steps; Stillwater: 8,760/8,760 steps with strict balance checks |
| Multi-biome annual parity | 1,102/1,104 variables in tolerance | 16 biomes × 69 variables against Fortran histories |
| Per-timestep Fortran parity | 35 harnesses passed | Instrumented CTSM references across physics and biogeochemistry; documented exceptions remain |
| Apple Metal fleet | 128 device harnesses passed | Device-versus-CPU parity on Apple Silicon; Metal uses `Float32` |

The three `Broken` results in the strict CPU suite are explicit optional-device
skips for CUDA, AMDGPU, and Metal. They are not failed tests or missing fixtures.
The complete validation record and reproduction notes are in
[`docs/HANDOFF_2026-07-24.md`](docs/HANDOFF_2026-07-24.md).

### What these results establish

- The tested single-point and single-column CLM5 configurations run through full
  annual cycles without fatal water-balance failures.
- The validated CPU implementation agrees closely with instrumented Fortran CTSM
  at matched timesteps and across annual multi-biome history variables.
- The tested GPU kernels and composite driver paths agree with same-precision CPU
  references on real hardware.
- Forward-mode and reverse-mode AD paths are exercised through substantial model
  components and driver workflows.

### What they do not establish

- Universal parity for every flag combination, land unit, forcing product, or
  transient multi-century simulation.
- Scientific validity of configurations that have not been compared against CTSM
  or observations.
- Bit identity between Julia and Fortran. Compiler, solver, precision, and
  operation-order differences make that an inappropriate general target.
- Production readiness of every optional FATES, crop, dynamic-land-use, methane,
  or coupled ice-sheet path.

## Capabilities

| Area | Implemented scope | Validation posture |
|---|---|---|
| Surface physics | Canopy and bare-ground radiation, turbulence, photosynthesis, LUNA, plant hydraulics | Broad unit, annual, and Fortran-parity coverage |
| Hydrology | Infiltration, runoff, drainage, water table, aquifer/bedrock options, hillslope routing | Annual and process-level parity coverage; some coupled modes remain opt-in |
| Snow and ice | Multilayer snow hydrology, compaction, phase change, SNICAR, glaciers | Seasonal and per-timestep parity coverage |
| Lakes and urban | Lake thermodynamics/hydrology and CLMU canyon energy balance | Dedicated strict-data fixtures and robustness gates |
| Biogeochemistry | CN pools and fluxes, CENTURY and MIMICS decomposition, N cycle, isotopes | CN and subsystem parity harnesses; coverage varies by process and site |
| Disturbance and trace gases | Li fire families, methane, VOC/MEGAN, dry deposition, dust | Major branches exercised against Fortran; see subsystem reports for residuals |
| Vegetation | Satellite phenology, crop lifecycle, CNDV, transient subgrid, FATES | Core paths run; deep optional and long-timescale configurations are less mature |
| Differentiation | ForwardDiff and Enzyme support, calibration framework, smoothed physics alternatives | Gradient and driver-level validation in tested configurations |
| Parallel execution | KernelAbstractions CPU/GPU kernels, MPI domain split, multi-GPU plumbing | CPU/MPI CI plus real-hardware Metal validation; CUDA/AMDGPU are optional extensions |

The source preserves CTSM variable names where practical, making routines such as
`h2osoi_liq`, `t_soisno`, and `qflx_*` traceable back to the Fortran model.

## Installation

CLM.jl supports Julia 1.10 and later. Until the package is registered, install it
directly from GitHub:

```julia
using Pkg
Pkg.add(url = "https://github.com/DarriEy/CLM.jl")
```

For development:

```bash
git clone https://github.com/DarriEy/CLM.jl.git
cd CLM.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

The committed `Project.toml` declares the supported dependency ranges. Manifests
are intentionally not committed, so each Julia version resolves a compatible
environment.

## Input data

An offline run requires CTSM-compatible NetCDF inputs:

- surface data (`fsurdat`);
- CLM parameter data (`paramfile`);
- atmospheric forcing (`fforcing`);
- optional snow optics, snow aging, restart, fire, methane, crop, or transient
  land-use inputs for the corresponding configurations.

These large upstream datasets are not distributed in the Git repository. The
test harnesses locate the project data archive through `SYMFLUENCE_DATA`:

```bash
export SYMFLUENCE_DATA=/absolute/path/to/SYMFLUENCE_data
```

Without the external archive, self-contained tests still run and data-dependent
tests report their gates explicitly. For strict release validation, stage the
fixtures described in [`docs/HANDOFF_2026-07-24.md`](docs/HANDOFF_2026-07-24.md).

## Running an offline simulation

The high-level API is currently namespaced rather than exported:

```julia
using CLM
using Dates

state = CLM.clm_run!(
    fsurdat = "/path/to/surfdata.nc",
    paramfile = "/path/to/clm5_params.nc",
    fforcing = "/path/to/clmforc.2003.nc",
    fhistory = "clm_history.nc",
    start_date = DateTime(2003, 1, 1),
    end_date = DateTime(2004, 1, 1),
    dtime = 1800,
    use_cn = false,
)
```

`CLM.clm_run!` returns the final `CLMInstances` state and writes the requested
history file. Important defaults follow the CLM5 configuration where available:
LUNA and plant hydraulic stress resolve conditionally, the CLM5 zero-flux lower
boundary is used, and balance violations are fatal by default.

For lower-level experiments, initialize the state tree directly:

```julia
inst, bounds, filters, time_manager = CLM.clm_initialize!(
    fsurdat = "/path/to/surfdata.nc",
    paramfile = "/path/to/clm5_params.nc",
    use_cn = true,
)
```

The driver and initialization keyword documentation in
[`src/driver/clm_run.jl`](src/driver/clm_run.jl) and
[`src/driver/clm_initialize.jl`](src/driver/clm_initialize.jl) is the authoritative
API reference while a Documenter site is being prepared.

## GPU execution

The base package has no hard GPU dependency. CUDA and AMDGPU register through
package extensions; Metal is carried by the dedicated script environment and
uses `Float32` because Apple GPUs do not support `Float64` compute.

To validate the available GPU on a machine:

```bash
julia --project=scripts -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=scripts scripts/gpu_validate.jl
```

On Apple Silicon, do **not** apply the CTSM/MPI
`HWLOC_COMPONENTS=-opencl` workaround to Julia Metal processes; it can make
`Metal.functional()` return `false`.

The repository contains component and end-to-end harnesses named
`scripts/gpu_validate_*`. GPU speedups require enough independent columns to
amortize kernel-launch and host/device transfer overhead; single-point workloads
usually belong on the CPU.

## Automatic differentiation and calibration

State and parameter structures are generally parameterized by `FT<:Real`, and
the codebase provides ForwardDiff and Enzyme paths plus a calibration framework.
A typical workflow constructs `CalibrationParameter` and `CalibrationTarget`
objects, then calls the namespaced optimizer:

```julia
using CLM

problem = CLM.CalibrationProblem(
    params = parameters,
    targets = targets,
    fsurdat = "/path/to/surfdata.nc",
    paramfile = "/path/to/clm5_params.nc",
)

result = CLM.calibrate(problem; maxiter = 20)
```

The concrete setters, getters, bounds, and observation data are application
specific. See `src/calibration/` and the `scripts/calibrate_*` examples before
building a new workflow. Differentiability is configuration-dependent: discrete
events and optional process branches may require the supplied smoothing or
compositional reverse-mode paths.

## Testing and validation

Run the ordinary package test suite with:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

The direct form used during release validation is:

```bash
julia --project=. --check-bounds=yes \
    -e 'using Test; include("test/runtests.jl")'
```

With the external data archive:

```bash
SYMFLUENCE_DATA=/absolute/path/to/SYMFLUENCE_data \
    julia --project=. --check-bounds=yes \
    -e 'using Test; include("test/runtests.jl")'
```

Validation is layered rather than represented by a single coverage number:

1. Unit and analytic-oracle tests exercise translated routines and invariants.
2. Integration tests drive initialization, complete timesteps, restarts, AD, and
   optional process configurations.
3. Strict-data tests exercise official urban, lake, snow, glacier, and multisite
   fixtures.
4. Fortran-parity harnesses compare against instrumented CTSM state boundaries.
5. Annual domain runs compare free-running Julia and Fortran history variables.
6. Real-device GPU harnesses compare device and same-precision CPU results.

CI intentionally resolves fresh environments for Julia 1.10 and latest stable
Julia, then runs the package tests and a two-rank MPI bit-identity smoke test.

## Repository layout

```text
src/
  constants/       Physical constants, control flags, and precision policy
  types/           Structure-of-arrays model state and parameter containers
  infrastructure/  I/O, solvers, kernels, backends, MPI, and subgrid machinery
  biogeophys/      Radiation, turbulence, hydrology, snow, soil, and photosynthesis
  biogeochem/      Carbon, nitrogen, decomposition, fire, methane, and vegetation
  fates/            FATES ecosystem-demography port and CLM coupling
  driver/           Initialization, timestep orchestration, and offline runs
  calibration/      AD rules, parameter injection, objectives, and optimization
test/               Unit, integration, strict-data, AD, backend, and MPI tests
scripts/            Parity, diagnostics, calibration, GPU, and validation tools
docs/               Validation reports, audits, residual analyses, and handoff notes
test_inputs/        Small redistributable fixtures
```

## Known limitations

- The package is unreleased and its public API is not yet stabilized.
- Large CTSM inputs and instrumented Fortran references must be staged separately.
- Validation breadth is uneven: biogeophysics and core CN have the strongest
  coverage; deep FATES toggles and some crop/transient/coupled paths have less.
- Apple Metal uses `Float32`; CPU, CUDA, and AMDGPU configurations generally use
  `Float64`. Cross-precision results should be assessed scientifically, not only
  by device-versus-CPU parity.
- Some historical reports under `docs/` are explicitly marked as snapshots and
  may describe issues fixed by later commits. Start with the current handoff.
- This is a source port of CTSM, not an official CESM or NCAR product.

The detailed remaining validation tail is tracked in
[`docs/HANDOFF_2026-07-24.md`](docs/HANDOFF_2026-07-24.md) and
[`docs/FORTRAN_VALIDATION_BACKLOG.md`](docs/FORTRAN_VALIDATION_BACKLOG.md).

## Provenance

The port was developed through a long-running, agentic workflow: coding agents
read the upstream Fortran, translated modules, added tests, generated parity
harnesses, and iterated against CTSM references. Human direction focused on scope,
validation targets, and release decisions rather than line-by-line authorship.

For traceability, the repository preserves upstream naming, porting notes,
instrumented-reference recipes, and a detailed commit history. Useful starting
points are:

- [`PORTING_LOG.md`](PORTING_LOG.md) — translation history and conventions;
- [`PRD_CLM_JULIA_PORT.md`](PRD_CLM_JULIA_PORT.md) — original project scope;
- [`docs/HANDOFF_2026-07-24.md`](docs/HANDOFF_2026-07-24.md) — current validation record;
- [`docs/MULTISITE_PARITY_MATRIX.md`](docs/MULTISITE_PARITY_MATRIX.md) — annual biome matrix;
- [`docs/CH4_FIRE_PARITY.md`](docs/CH4_FIRE_PARITY.md) — methane and fire oracle work;
- [`docs/GPU_PORT_GAP.md`](docs/GPU_PORT_GAP.md) — accelerator scope and limitations.

## Contributing

Issues and pull requests are welcome. A useful contribution should:

- name the CTSM configuration and upstream routine being matched;
- include a focused regression test or parity oracle;
- distinguish self-consistency from comparison against external ground truth;
- preserve default-path results unless an intentional scientific correction is
  documented;
- pass Julia 1.10, latest stable Julia, and the MPI smoke test.

Run `git diff --check` and the relevant targeted tests before the full suite. Do
not commit machine-specific manifests, generated parity output, or proprietary/
large upstream datasets.

## Citation

Until a versioned release and `CITATION.cff` are published, cite the repository:

> Eythorsson, D. (2026). *CLM.jl: A differentiable Julia port of the Community
> Land Model* (version 0.1.0, research software). GitHub.
> https://github.com/DarriEy/CLM.jl

For scientific descriptions of CLM5, also cite the appropriate CTSM/CLM papers
for the configuration and processes used.

## License and upstream attribution

The original Julia contributions in this repository are licensed under the
[MIT License](LICENSE). CLM.jl is a derivative port of CTSM/CLM5; the upstream
copyright and BSD 3-Clause terms retained in [`NOTICE`](NOTICE) continue to apply
to derived portions of the work.

CLM.jl is based on the Community Land Model developed by NCAR and the broader
CESM community. It is not affiliated with or endorsed by NCAR, CESM, or the CTSM
development team.
