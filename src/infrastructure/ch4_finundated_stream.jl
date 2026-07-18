# ==========================================================================
# Ported from: src/cpl/share_esmf/ch4FInundatedStreamType.F90
#   (ch4finundatedstream_type: Init / InitAllocate / UseStreams / CalcFinundated)
#
# Reads the satellite-derived inundated-fraction regression stream used by the
# two non-default finundation methods in the CH4 model:
#   * TWS_inversion : finundated = FWS_TWS_A * TWS + FWS_TWS_B          (Swenson)
#   * ZWT_inversion : finundated = F0 * exp(-ZWT/ZWT0) + P3*qflx_surf_lag (Koven)
# Coefficients are per-gridcell fields on the CTSM 0.9x1.25 grid
# (`finundated_inversiondata_0.9x1.25_c170706.nc`). The `finundated` regression
# DECOUPLES inundation from the column's own fill-and-spill h2osfc pond, so a
# wetland gridcell can inundate directly from the observed regression even when
# the single-point soil column drains (h2osfc method gives finundated == 0).
#
# FIDELITY NOTE
# -------------
# CTSM regrids this stream onto the model mesh with ESMF (mapalgo='bilinear')
# and time-interpolates with shr_strdata (a single 1996 slice, hardcoded to
# mcdate 1996-12-31 in Init -> effectively a constant field). This reader does
# NEAREST-NEIGHBOUR spatial selection from the file's own LATIXY/LONGXY grid,
# which is exact for a single-point run whose gridcell coincides with a file
# cell, and is the physically-correct constant-in-time field. A bit-for-bit
# match with ESMF's bilinear regridder on a general (multi-cell) grid is NOT
# claimed (same caveat as ndep_stream.jl / fire_streams.jl).
#
# Public:
#   CH4FInundatedStream            — per-gridcell regression coefficient arrays
#   read_ch4_finundated_stream!    — populate coefficients for given gridcells
# ==========================================================================

"""
    CH4FInundatedStream

Per-gridcell regression coefficients for the satellite-inversion finundation
methods, mirroring the private pointers of `ch4finundatedstream_type`:

* `zwt0_gdc`, `f0_gdc`, `p3_gdc`         — ZWT_inversion  (from ZWT0/F0/P3)
* `fws_slope_gdc`, `fws_intercept_gdc`   — TWS_inversion  (from FWS_TWS_A/FWS_TWS_B)

Each array is gridcell-length (indexed by `col_gridcell[c]`). Unpopulated for the
default `FINUNDATION_MTD_H2OSFC` method, in which case the CH4 driver never reads
these (the h2osfc branch of `CalcFinundated` uses `frac_h2osfc` only).
"""
Base.@kwdef mutable struct CH4FInundatedStream{V<:AbstractVector{<:Real}}
    active::Bool = false
    filename::String = ""
    zwt0_gdc::V          = Float64[]
    f0_gdc::V            = Float64[]
    p3_gdc::V            = Float64[]
    fws_slope_gdc::V     = Float64[]
    fws_intercept_gdc::V = Float64[]
end

# Cyclic great-circle-ish nearest index into the file's 2-D LATIXY/LONGXY grid.
# LONGXY is degrees_east in [0,360); we compare with a wrapped longitude delta so
# a gridcell near the 0/360 seam still selects its true neighbour.
function _ch4finund_nearest(latixy::AbstractMatrix, longxy::AbstractMatrix,
                            xlat::Real, xlon::Real)
    xlon360 = mod(Float64(xlon), 360.0)
    best = CartesianIndex(1, 1)
    bestd = Inf
    @inbounds for I in CartesianIndices(latixy)
        dlat = Float64(latixy[I]) - Float64(xlat)
        dlon = mod(Float64(longxy[I]) - xlon360 + 180.0, 360.0) - 180.0
        d = dlat * dlat + dlon * dlon
        if d < bestd
            bestd = d
            best = I
        end
    end
    return best
end

"""
    read_ch4_finundated_stream!(stream, filename, lats_deg, lons_deg; method)

Read the finundation-inversion regression coefficients for the gridcells at
`lats_deg` / `lons_deg` (degrees north / degrees east, one entry per gridcell)
from `filename`, into `stream`.

`method` selects which coefficient set is read (matching the Fortran `Init`,
which only allocates+reads the fields the active method needs):
* `FINUNDATION_MTD_TWS_INVERSION` -> `fws_slope_gdc` (FWS_TWS_A), `fws_intercept_gdc` (FWS_TWS_B)
* `FINUNDATION_MTD_ZWT_INVERSION` -> `zwt0_gdc` (ZWT0), `f0_gdc` (F0), `p3_gdc` (P3)

For `FINUNDATION_MTD_H2OSFC` this is a no-op (the h2osfc branch needs no stream).
"""
function read_ch4_finundated_stream!(stream::CH4FInundatedStream,
                                     filename::AbstractString,
                                     lats_deg::AbstractVector{<:Real},
                                     lons_deg::AbstractVector{<:Real};
                                     method::Int = FINUNDATION_MTD_TWS_INVERSION)
    if method == FINUNDATION_MTD_H2OSFC
        return stream
    end
    isfile(filename) ||
        error("read_ch4_finundated_stream!: stream file not found: $filename")
    length(lats_deg) == length(lons_deg) ||
        error("read_ch4_finundated_stream!: lats/lons length mismatch")
    ng = length(lats_deg)

    NCDataset(filename, "r") do ds
        # NOTE: read via `.var` (raw, no CF decoding). The coefficient variables
        # carry DESCRIPTIVE-STRING `scale_factor`/`add_offset` attributes (e.g.
        # "TWS", "exp(-ZWT/ZWT0)") that document the regression -- they are NOT
        # numeric CF transforms, and NCDatasets' CF path errors trying to apply
        # them (`Float64 * String`). The stored values are already the fitted
        # coefficients, so raw access is correct.
        latixy = Float64.(Array(ds["LATIXY"].var[:, :]))
        longxy = Float64.(Array(ds["LONGXY"].var[:, :]))

        # Nearest file cell for each model gridcell (constant across the single
        # time slice; take the first time level to mirror the hardcoded mcdate).
        idx = [_ch4finund_nearest(latixy, longxy, lats_deg[g], lons_deg[g]) for g in 1:ng]

        _sel(varname) = begin
            # Vars are (time, lsmlat, lsmlon) in the file -> (lsmlon, lsmlat, time)
            # in NCDatasets; take time level 1.
            slab = Float64.(Array(ds[varname].var[:, :, 1]))
            [slab[idx[g]] for g in 1:ng]
        end

        if method == FINUNDATION_MTD_TWS_INVERSION
            stream.fws_slope_gdc     = _sel("FWS_TWS_A")
            stream.fws_intercept_gdc = _sel("FWS_TWS_B")
        elseif method == FINUNDATION_MTD_ZWT_INVERSION
            stream.zwt0_gdc = _sel("ZWT0")
            stream.f0_gdc   = _sel("F0")
            stream.p3_gdc   = _sel("P3")
        else
            error("read_ch4_finundated_stream!: unknown finundation method $method")
        end
    end

    stream.active   = true
    stream.filename = String(filename)
    return stream
end
