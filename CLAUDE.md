# CLM.jl — Fortran CLM to Julia Port

## Project Overview
This is a Julia port of the Community Land Model (CLM/CTSM) from Fortran 90.
The goal is full process fidelity, GPU execution, and differentiability via AD.

## Fortran Source Location
The original Fortran CLM source is at:
`/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/installs/clm/`

## Architecture & Conventions

### Module Structure
All Julia source files are `include()`d into the top-level `CLM` module in `src/CLM.jl`.
Files are organized by category:
- `src/constants/` — Physical constants, control flags, parameters
- `src/types/` — Data structure definitions (state types)
- `src/infrastructure/` — Solvers, decomposition, filters, utilities
- `src/biogeophys/` — Biogeophysics process modules
- `src/biogeochem/` — Biogeochemistry process modules
- `src/driver/` — Main timestep driver

### Translation Patterns

**Fortran types → Julia structs (mutable, @kwdef):**
```julia
Base.@kwdef mutable struct TemperatureData
    t_soisno::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    t_grnd::Vector{Float64} = Float64[]
end
```

**Fortran filter loops → mask-based loops:**
```julia
for c in eachindex(mask_soil)
    mask_soil[c] || continue
    # physics for column c
end
```

**Fortran subroutines → Julia bang functions:**
```julia
function soil_temperature!(temp::TemperatureData, mask::BitVector, bounds::UnitRange{Int})
```

### Key Decisions
- **Preserve Fortran variable names** for traceability back to the original code
- **SoA (Structure of Arrays)** layout — matches Fortran and is GPU-friendly
- **BitVector masks** replace integer filter arrays — GPU-compatible, no dynamic rebuild
- **Snow layers padded to fixed max** — no dynamic array resizing
- **Discontinuities kept as-is** in Phase 1 — smoothing for AD happens in Phase 3
- **Float64 throughout** — parametric types for Float32/dual numbers come later

### Adding a New Module
1. Read the Fortran source completely
2. Create the Julia file in the appropriate `src/` subdirectory
3. Add `include("path/to/file.jl")` in `src/CLM.jl`
4. Add `include("test_modulename.jl")` in `test/runtests.jl`
5. Write tests that verify correctness (finite-difference derivative checks where applicable)
6. Ensure ALL existing tests still pass before committing

## Build & Test
```bash
julia --project=. -e 'using Test; include("test/runtests.jl")'
```

## PRD
See `PRD_CLM_JULIA_PORT.md` for the full porting plan, module dependency graph, and tier ordering.
