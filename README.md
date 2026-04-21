# CLM.jl

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Julia](https://img.shields.io/badge/Julia-1.10%2B-purple.svg)](https://julialang.org)
[![Tests](https://img.shields.io/badge/Tests-15%2C528_passing-brightgreen.svg)]()

A differentiable Julia port of the [Community Land Model version 5](https://www.cesm.ucar.edu/models/clm/) (CLM5/CTSM) — the land surface component of the Community Earth System Model (CESM).

CLM.jl reproduces the full CLM5 biogeophysics in pure Julia, enabling **automatic differentiation** for gradient-based parameter calibration, **GPU-ready architecture** with BitVector masks, and **composability** with the Julia scientific ecosystem (ForwardDiff.jl, Enzyme.jl, Flux.jl, Optim.jl).

## Validation Against Fortran CLM5

CLM.jl has been validated in SP mode against Fortran CLM5 (CTSM release-clm5.0) for a full-year simulation at the Bow at Banff site (51.2°N, 115.6°W, 2003, restart-initialized from a 2002 cold-start spinup). The figure below shows 7-day running averages for 12 key output variables:

![CLM.jl vs Fortran CLM5 — Bow at Banff 2003 (restart init, 7-day avg)](scripts/comparison.png)

| Variable | Mean Bias | RMSE | Correlation |
|---|---|---|---|
| 2-m Air Temperature | +0.01 K | 0.15 K | 0.999 |
| Absorbed Solar Radiation | +0.5% | 1.2 W/m² | 0.999 |
| Latent Heat Flux | +3% | 1.8 W/m² | 0.98 |
| Canopy Transpiration | +12% | 2.4 W/m² | 0.97 |
| Sensible Heat, Ground Temp, Snow, Runoff | < 5% | — | > 0.98 |

## Features

### Biogeophysics (Validated)
- Surface albedo with two-stream canopy radiative transfer
- Snow radiative transfer (SNICAR) with aerosol-snow interactions
- Turbulent fluxes via Monin-Obukhov similarity (bare ground, canopy, lake, urban)
- Photosynthesis (Farquhar/Collatz) with Ball-Berry and Medlyn stomatal conductance
- LUNA photosynthetic optimization (Vcmax/Jmax acclimation)
- Plant hydraulic stress (PHS)
- Snow hydrology (12 layers, compaction, grain-size evolution)
- Soil hydrology (Richards equation, Clapp-Hornberger / van Genuchten)
- Soil temperature with phase change
- Lake temperature and hydrology
- Urban canyon energy balance (CLMU)
- Irrigation, hillslope lateral flow

### Biogeochemistry (Implemented, Not Yet Validated)
- Carbon-nitrogen cycling (CN mode) with allocation, respiration
- Decomposition (Century cascade and MIMICS)
- Nitrification-denitrification, nitrogen leaching
- Fire (Li 2014), gap mortality, dynamic vegetation (CNDV)
- Methane (CH4), VOC emissions, dust emission

### AD & Calibration
- All 50+ data structs parameterized on `{FT<:Real}` for dual-number propagation
- ~520 smoothed discontinuities (`smooth_max`, `smooth_min`, `smooth_heaviside`)
- Pure-Julia LU fallback for banded matrix solves with Dual numbers
- Built-in `CalibrationProblem` framework with ForwardDiff gradients
- 7 tunable parameters: `vcmax25_scale`, `medlyn_slope`, `csoilc`, `jmax25top_sf`, `baseflow_scalar`, `fff`, `ksat_scale`

## Quick Start

```julia
using CLM

# Run a full simulation
clm_run!(
    fsurdat  = "path/to/surfdata.nc",
    paramfile = "path/to/params.nc",
    fforcing  = "path/to/forcing/",
    fhistory  = "output/history.nc"
)

# Or initialize and step manually
inst, bounds, filt, tm = clm_initialize!(
    fsurdat  = "path/to/surfdata.nc",
    paramfile = "path/to/params.nc"
)
clm_drv!(config, inst, filt, filt_ia, bounds, ...)
```

## Gradient-Based Calibration

```julia
using CLM, ForwardDiff

prob = CalibrationProblem(
    params = [
        CalibrationParameter("vcmax25_scale", 1.0, setter!, (0.5, 2.0)),
        CalibrationParameter("medlyn_slope", 4.1, setter!, (1.0, 10.0)),
    ],
    targets = [CalibrationTarget("LH", getter, obs_LH, 1.0)],
    fsurdat = "surfdata.nc",
    paramfile = "params.nc",
    fforcing = "forcing/"
)

# AD gradients + Armijo line search
result = calibrate(prob; maxiter=20)
```

## Architecture

```
src/
  constants/       Physical constants, control flags, PFT parameters
  types/           50+ mutable structs (SoA layout, {FT<:Real} parameterized)
  infrastructure/  Solvers, I/O, initialization, filters, subgrid
  biogeophys/      Radiation, turbulence, hydrology, snow, soil, photosynthesis
  biogeochem/      Phenology, decomposition, nutrient cycling, fire, CH4, VOC
  driver/          Timestep driver, initialization, top-level run
  calibration/     AD-based calibration framework
```

**Key design decisions:**
- **SoA (Structure of Arrays)** layout matching Fortran CLM for GPU compatibility
- **BitVector masks** replace Fortran integer filter arrays — no dynamic rebuild, GPU-ready
- **Fortran variable names preserved** for line-by-line traceability (`h2osoi_liq`, `t_soisno`, etc.)
- **Fixed-size snow arrays** padded to `nlevsno=12` — no dynamic resizing

## Performance

| Configuration | Wall-clock (1 year, single-point) |
|---|---|
| Fortran CLM5 (ifort -O2) | ~2 s |
| CLM.jl (Float64) | ~4 s |
| CLM.jl + ForwardDiff (1 param) | ~8 s |
| CLM.jl + ForwardDiff (7 params) | ~28 s |

## Testing

```bash
julia --project=. -e 'using Test; include("test/runtests.jl")'
```

15,528 tests: unit tests, end-to-end SP simulation, AD gradient verification (6 climate scenarios), calibration framework, parameter recovery, CN integration, and Enzyme feasibility.

## Requirements

- Julia 1.10+
- NCDatasets.jl, ForwardDiff.jl, JSON.jl
- CLM5 surface data and parameter files (NetCDF)

## Citation

If you use CLM.jl, please cite:

> Eythorsson, D. (2026). CLM.jl: A differentiable Julia port of the Community Land Model. *Geoscientific Model Development* (in preparation).

## Contributing

Issues and pull requests are welcome. Please run the full test suite before submitting changes:

```bash
julia --project=. -e 'using Test; include("test/runtests.jl")'
```

## License

[MIT](LICENSE) — free for academic and commercial use.

## Acknowledgements

CLM.jl is based on the [Community Land Model](https://github.com/ESCOMP/CTSM) developed by the National Center for Atmospheric Research (NCAR) and the broader CESM community. The original Fortran CLM5 is described in Lawrence et al. (2019).
