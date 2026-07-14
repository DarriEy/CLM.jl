# =============================================================================
# Synthetic lightning (lnfm) + population-density (hdm) streams for the
# Fortran CTSM fire reference run.
#
# WHY THIS EXISTS
# ---------------
# CTSM's Li-family fire methods (`fire_method='li2016crufrc'`) need two data
# streams:
#   * lnfm -- lightning flash density  (counts/km^2/hr), stream var 'lnfm'
#   * hdm  -- human population density (counts/km^2),    stream var 'hdm'
# The real inputdata files
#   atm/datm7/NASA_LIS/clmforc.Li_2012_climo1995-2011.T62.lnfm_Total_*.nc
#   lnd/clm2/firedata/clmforc.Li_2017_HYDEv3.2_CMIP6_hdm_0.5x0.5_*.nc
# are NOT present in the local inputdata tree (fire was never switched on here,
# so they were never downloaded) and could not be re-fetched.
#
# Resolution -- the SAME approach PR #221 used for the missing CMIP6 ndep file:
# build a synthetic stream that is UNIFORM IN SPACE and CONSTANT IN TIME, and
# drive BOTH codes with it.
#
# Why uniform+constant is the right choice and not a cop-out:
#   CTSM regrids these streams with ESMF (`mapalgo='bilinear'`) and
#   time-interpolates with shr_strdata (lnfm: tintalgo='linear', taxmode='cycle';
#   hdm: tintalgo='nearest', taxmode='extend').  For a field that is uniform in
#   space and constant in time, the bilinear weights sum to 1 and the temporal
#   weights sum to 1, so BOTH codes recover the stream value EXACTLY.  Fortran's
#   forc_lnfm/forc_hdm are therefore known analytically, and the parity test
#   measures the FIRE PHYSICS rather than ESMF interpolation fidelity -- which is
#   what we actually need to test.
#
# These are SYNTHETIC values.  They are not a real lightning/population dataset
# and this script does not claim they are.  They are chosen to be physically
# defensible for the Bow-at-Banff basin so the reference run is non-vacuous:
#
#   LNFM = 4.57e-4 counts/km^2/hr
#     = 4.0 flashes/km^2/yr, the annual-mean TOTAL (IC+CG) flash density.
#       Mid-latitude continental North America is ~2-6 flashes/km^2/yr in the
#       NASA LIS/OTD climatology; the Canadian Rockies sit in the lower-middle of
#       that range.  CLM applies its own cloud-to-ground fraction internally
#       (the 0.22 `ignition_efficiency` in the Li ignition term), so this is the
#       TOTAL flash density, not the CG density.
#
#   HDM = 4.2 counts/km^2
#     = ~9,300 people (Banff townsite ~8,300 + Lake Louise ~1,000) over the
#       ~2,210 km^2 Bow-at-Banff drainage.
#
# The grid is fv0.9x1.25 (288 lon x 192 lat), reusing the EXISTING ESMF mesh
#   share/meshes/fv0.9x1.25_141008_polemod_ESMFmesh.nc
# so no new mesh file has to be generated.  The lon/lat axes are copied verbatim
# from the ndep stream on the same grid, guaranteeing mesh consistency.
#
# Usage:
#   julia --project=. scripts/validation/make_fire_streams.jl [outdir]
# =============================================================================

using NCDatasets

const NDEP_REF = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/ndepdata/" *
                 "fndep_clm_ZERO_0.9x1.25_yr2000_CLMjl-parity.nc"

# --- the two synthetic constants (see header for provenance) ---
const LNFM_VALUE = 4.57e-4   # counts/km^2/hr  (= 4.0 flashes/km^2/yr)
const HDM_VALUE  = 4.2       # counts/km^2

function main(outdir::AbstractString)
    mkpath(outdir)

    # `.var` on the read side too: ds["time"] would CF-decode to DateTime, which then
    # cannot be written back into a Float64 axis (it lands as fill).
    lon, lat, tvals = NCDataset(NDEP_REF, "r") do ds
        (Array(ds["lon"]), Array(ds["lat"]), Array(ds["time"].var))
    end
    nlon, nlat, nt = length(lon), length(lat), length(tvals)
    @info "grid from ndep reference" nlon nlat nt

    write_stream(joinpath(outdir, "clmforc.UNIFORM_lnfm_0.9x1.25_CLMjl-parity.nc"),
                 lon, lat, tvals, "lnfm", LNFM_VALUE, "counts/km^2/hr",
                 "Lightning flash density (UNIFORM synthetic - CLM.jl fire parity)")

    write_stream(joinpath(outdir, "clmforc.UNIFORM_hdm_0.9x1.25_CLMjl-parity.nc"),
                 lon, lat, tvals, "hdm", HDM_VALUE, "counts/km^2",
                 "Human population density (UNIFORM synthetic - CLM.jl fire parity)")

    @info "done" LNFM_VALUE HDM_VALUE outdir
end

function write_stream(path, lon, lat, tvals, varname, value, units, longname)
    isfile(path) && rm(path)
    NCDataset(path, "c") do ds
        ds.dim["lon"]  = length(lon)
        ds.dim["lat"]  = length(lat)
        ds.dim["time"] = length(tvals)

        v = defVar(ds, "lon", Float64, ("lon",))
        v[:] = lon
        v.attrib["units"] = "degrees_east"
        v.attrib["long_name"] = "longitude"

        v = defVar(ds, "lat", Float64, ("lat",))
        v[:] = lat
        v.attrib["units"] = "degrees_north"
        v.attrib["long_name"] = "latitude"

        # NB: set the CF time attributes FIRST, then write through `.var` to bypass
        # NCDatasets' CF date encoding -- writing Float64 offsets to a variable that
        # NCDatasets considers a CF time axis otherwise silently stores fill values,
        # which CTSM rejects with "shr_stream_verifyTCoord: elapsed seconds must be
        # increasing".
        v = defVar(ds, "time", Float64, ("time",))
        v.attrib["units"] = "days since 2000-01-01 00:00:00"
        v.attrib["calendar"] = "gregorian"
        v.attrib["long_name"] = "time"
        v.var[:] = tvals

        v = defVar(ds, varname, Float64, ("lon", "lat", "time"))
        v[:, :, :] = fill(value, length(lon), length(lat), length(tvals))
        v.attrib["units"] = units
        v.attrib["long_name"] = longname

        ds.attrib["title"] = longname
        ds.attrib["note"] = "SYNTHETIC uniform-in-space, constant-in-time stream generated by " *
                            "CLM.jl scripts/validation/make_fire_streams.jl. NOT a real dataset. " *
                            "Used to drive BOTH Fortran CTSM and CLM.jl identically so that fire " *
                            "PHYSICS can be compared without ESMF-regrid fidelity confounding it."
        ds.attrib["constant_value"] = value
    end
    @info "wrote" path varname value
end

main(length(ARGS) >= 1 ? ARGS[1] : pwd())
