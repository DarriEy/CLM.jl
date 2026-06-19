# ==========================================================================
# gen_glacier_surfdata.jl — Generate GLACIER-enabled surfdata for broad-config parity
#
# The Bow surfdata is 100% natural vegetation (PCT_NATVEG=100, all other
# landunit fractions 0). To exercise the GLACIER (istice) land-unit path we
# clone the Bow surfdata and rewrite the landunit fractions so the gridcell is
# a pure glacier column.
#
# GLC_MEC mode is active in the Bow namelists (maxpatch_glc=10), so the glacier
# landunit is built from PCT_GLACIER together with the elevation-class
# distribution PCT_GLC_MEC. The Bow surfdata already carries
# PCT_GLC_MEC = (100,0,...,0) — all glacier in elevation class 0 (the
# lowest/bare class) — and GLACIER_REGION=0 (single_at_atm_topo behavior).
# We keep that distribution and just turn the landunit on via PCT_GLACIER=100.
#
# Produces under test_inputs/glacier/:
#   surfdata_glacier100.nc — PCT_GLACIER=100 (pure glacier column; isolates the
#                            istice land-unit path: generic soil/snow column
#                            with glacier albedo/runoff behavior)
#
# Keeps Bow's lat/lon, TOPO, soil texture and monthly phenology so a Fortran
# reference run is a one-flag change (set PCT_GLACIER on the same site).
#
# Usage:
#   julia --project=. scripts/gen_glacier_surfdata.jl [path/to/source/surfdata_clm.nc]
# Default source = the Bow-at-Banff surfdata in SYMFLUENCE_data.
# ==========================================================================

using NCDatasets

const DEFAULT_SRC = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/" *
                    "domain_Bow_at_Banff_lumped/settings/CLM/parameters/surfdata_clm.nc"

const OUTDIR = joinpath(@__DIR__, "..", "test_inputs", "glacier")

# Write a single value into a (lsmlat, lsmlon) scalar field, preserving shape.
function set_scalar!(ds::NCDataset, name::String, val::Real)
    haskey(ds, name) || return
    v = ds[name]
    v[:] = fill(Float64(val), size(v))
    return
end

"""
    clone_surfdata(src, dst; pct_natveg, pct_glacier)

Byte-copy `src` to `dst`, then overwrite the landunit-fraction scalars so the
gridcell is a pure glacier column. PCT_GLC_MEC / TOPO_GLC_MEC / GLACIER_REGION
are left untouched (Bow already has all glacier in MEC class 0). Everything
else (soil texture, PFT distribution, monthly LAI/SAI/heights, lat/lon, TOPO)
is preserved so the only varied axis is the landunit structure.
"""
function clone_surfdata(src::String, dst::String; pct_natveg::Float64,
                        pct_glacier::Float64)
    isfile(src) || error("source surfdata not found: $src")
    cp(src, dst; force = true)
    NCDataset(dst, "a") do ds
        set_scalar!(ds, "PCT_NATVEG", pct_natveg)
        set_scalar!(ds, "PCT_GLACIER", pct_glacier)
        set_scalar!(ds, "PCT_CROP",    0.0)
        set_scalar!(ds, "PCT_LAKE",    0.0)
        set_scalar!(ds, "PCT_WETLAND", 0.0)
        set_scalar!(ds, "PCT_URBAN",   0.0)
        set_scalar!(ds, "PCT_OCEAN",   0.0)
    end
    return dst
end

function main()
    src = length(ARGS) >= 1 ? ARGS[1] : DEFAULT_SRC
    mkpath(OUTDIR)

    f_glac = joinpath(OUTDIR, "surfdata_glacier100.nc")
    clone_surfdata(f_glac == src ? error("dst==src") : src, f_glac;
                   pct_natveg = 0.0, pct_glacier = 100.0)

    println("Source surfdata: ", src)
    println("Wrote: ", f_glac, "  (PCT_GLACIER=100, PCT_GLC_MEC class-0=100)")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
