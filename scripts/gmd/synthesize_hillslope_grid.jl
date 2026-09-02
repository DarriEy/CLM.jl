# ==========================================================================
# synthesize_hillslope_grid.jl — replace the ONLY upstream-derived content in
# the tracked hillslope fixtures (redistribution audit D02/D03) with analytic
# synthetic coordinates, and strip the leaked provenance attributes.
#
# The hillslope catena geometry in these files was always synthetic
# (make_hillslope_surfdata.jl); the LONGXY/LATIXY arrays were copied from the
# Aripuana surfdata (D02, S02) and the CTSM ne3np4 surface file (D03, S03),
# and a developer-absolute base path leaked into a global attribute. After
# this pass the files contain no upstream values:
#
#  - single-point file: a generic synthetic tropical point (-10, 300) —
#    deliberately nobody's study site;
#  - unstructured file: a Fibonacci sphere lattice of the same gridcell
#    count — uniform, deterministic, and obviously analytic.
#
# Real-grid variants for Fortran-side parity remain locally regenerable via
# scripts/make_hillslope_surfdata.jl from SYMFLUENCE_DATA / CESM_INPUTDATA;
# those outputs are not tracked.
#
#   julia --project=. scripts/gmd/synthesize_hillslope_grid.jl
# ==========================================================================

using NCDatasets

const REPO = joinpath(@__DIR__, "..", "..")

# Golden-angle (Fibonacci) sphere lattice — analytic, deterministic.
fib_lat(i, n) = asind(clamp(1 - 2 * (i - 0.5) / n, -1.0, 1.0))
fib_lon(i)    = mod(i * 137.50776405003785, 360.0)

function synthesize!(path::AbstractString)
    NCDataset(path, "a") do ds
        lat = ds["LATIXY"]; lon = ds["LONGXY"]
        n = length(lat)
        if n == 1
            lat[1] = -10.0
            lon[1] = 300.0
        else
            for i in 1:n
                lat[i] = fib_lat(i, n)
                lon[i] = fib_lon(i)
            end
        end
        haskey(ds.attrib, "base_surfdata") && delete!(ds.attrib, "base_surfdata")
        ds.attrib["note"] = "FULLY SYNTHETIC: idealized hillslope catena geometry " *
            "(1 hillslope, 4 cols, upland->lowland) on an analytic synthetic grid " *
            "(single point: generic (-10,300); multi-cell: Fibonacci sphere lattice). " *
            "No upstream data values. Real-grid variants for Fortran parity are " *
            "regenerated locally via scripts/make_hillslope_surfdata.jl."
        ds.attrib["grid_coordinates"] = "analytic synthetic (scripts/gmd/synthesize_hillslope_grid.jl)"
        println(path, ": ", n, " gridcell coordinate(s) replaced, provenance attrs sanitized")
    end
end

synthesize!(joinpath(REPO, "data", "hillslope", "hillslope_aripuana_synthetic_c.nc"))
synthesize!(joinpath(REPO, "data", "hillslope", "hillslope_ne3np4_synthetic_c.nc"))
