# ==========================================================================
# gen_lake_surfdata.jl — Generate LAKE-enabled surfdata for broad-config parity
#
# The Bow surfdata is 100% natural vegetation (PCT_NATVEG=100, all other
# landunit fractions 0). SYMFLUENCE can only synthesize that single-landunit
# natveg column. To exercise the (driver-wired, cold-start-ready) LAKE path
# we clone the Bow surfdata and rewrite the landunit fractions.
#
# This produces TWO variants under test_inputs/lake/:
#   surfdata_lake100.nc   — PCT_LAKE=100 (pure deep-lake column; isolates lake
#                           physics: lake_fluxes!/lake_temperature!/lake_hydrology!)
#   surfdata_mixed.nc     — PCT_NATVEG=50, PCT_LAKE=50 (exercises the
#                           multi-landunit subgrid build: one soil column + one
#                           lake column in the same gridcell)
#
# Both keep Bow's lat/lon, soil texture, LAKEDEPTH, and monthly phenology so a
# Fortran reference run is a one-flag change (set PCT_LAKE on the same site).
#
# Usage:
#   julia --project=. scripts/gen_lake_surfdata.jl [path/to/source/surfdata_clm.nc]
# Default source = the Bow-at-Banff surfdata in SYMFLUENCE_data.
# ==========================================================================

using NCDatasets

const DEFAULT_SRC = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/" *
                    "domain_Bow_at_Banff_lumped/settings/CLM/parameters/surfdata_clm.nc"

const OUTDIR = joinpath(@__DIR__, "..", "test_inputs", "lake")

"""
    clone_surfdata(src, dst; pct_natveg, pct_lake, lakedepth)

Byte-copy `src` to `dst`, then overwrite the landunit-fraction scalars
(PCT_NATVEG, PCT_LAKE and the remaining PCT_* set to 0) and LAKEDEPTH.
Everything else (soil texture, PFT distribution, monthly LAI/SAI/heights,
lat/lon) is preserved so the only varied axis is the landunit structure.
"""
function clone_surfdata(src::String, dst::String; pct_natveg::Float64,
                        pct_lake::Float64, lakedepth::Float64 = 10.0)
    isfile(src) || error("source surfdata not found: $src")
    cp(src, dst; force = true)
    NCDataset(dst, "a") do ds
        set_scalar!(ds, "PCT_NATVEG", pct_natveg)
        set_scalar!(ds, "PCT_LAKE",   pct_lake)
        set_scalar!(ds, "PCT_CROP",    0.0)
        set_scalar!(ds, "PCT_WETLAND", 0.0)
        set_scalar!(ds, "PCT_GLACIER", 0.0)
        set_scalar!(ds, "PCT_URBAN",   0.0)
        set_scalar!(ds, "PCT_OCEAN",   0.0)
        if haskey(ds, "LAKEDEPTH")
            set_scalar!(ds, "LAKEDEPTH", lakedepth)
        end
    end
    return dst
end

# Write a single value into a (lsmlat, lsmlon) scalar field, preserving shape.
function set_scalar!(ds::NCDataset, name::String, val::Real)
    haskey(ds, name) || return
    v = ds[name]
    v[:] = fill(Float64(val), size(v))
    return
end

function main()
    src = length(ARGS) >= 1 ? ARGS[1] : DEFAULT_SRC
    mkpath(OUTDIR)

    f_lake100 = joinpath(OUTDIR, "surfdata_lake100.nc")
    f_mixed   = joinpath(OUTDIR, "surfdata_mixed.nc")

    clone_surfdata(src, f_lake100; pct_natveg = 0.0,  pct_lake = 100.0)
    clone_surfdata(src, f_mixed;   pct_natveg = 50.0, pct_lake = 50.0)

    println("Source surfdata: ", src)
    println("Wrote: ", f_lake100, "  (PCT_LAKE=100)")
    println("Wrote: ", f_mixed,   "  (PCT_NATVEG=50, PCT_LAKE=50)")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
