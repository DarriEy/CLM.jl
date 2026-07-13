# ==========================================================================
# Ported from: src/cpl/share_esmf/ndepStreamMod.F90
#
# Atmospheric nitrogen-deposition stream: reads the `fndep_*` NetCDF stream
# file and fills `atm2lnd.forc_ndep_grc` (gN/m2/s) — the single entry point of
# nitrogen into the CLM ecosystem.
#
# Public functions:
#   ndep_stream_init!   — read the stream file (CNNDeposition's `ndep_init`)
#   ndep_interp!        — time+space interpolate to the model grid (`ndep_interp`)
#
# FIDELITY NOTE (read before trusting this for a new dataset)
# ----------------------------------------------------------
# CTSM regrids the stream with ESMF (`ndepmapalgo = 'bilinear'`) and
# time-interpolates with `shr_strdata` (`ndep_tintalgo = 'linear'`,
# `ndep_taxmode = 'cycle'`). This module reimplements *the same scheme* —
# bilinear in space, linear in time with cyclic wraparound — but it is NOT
# bit-validated against ESMF's regridder for a general (spatially varying)
# field, because doing so needs the ESMF weight file.
#
# What IS validated against Fortran (see scripts/fortran_parity_ndep.jl) is the
# reduction that matters for the port: for a stream that is uniform in space and
# constant in time, bilinear weights sum to 1 and the linear time weights sum to
# 1, so BOTH codes return the stream value exactly. The Bow N-cycle parity
# reference is generated with exactly such a stream, which isolates the N-cycle
# *physics* wiring from forcing-interpolation fidelity.
#
# Consequence: a spatially varying real (e.g. CMIP6) ndep file will be handled
# physically correctly here, but a *bit-for-bit* match with CTSM on such a file
# is not claimed.
# ==========================================================================

"""
    NDepStreamData

Holds the raw atmospheric N-deposition stream read from a CTSM `fndep_*` file,
plus the unit-conversion flag derived from the file's `units` attribute.

`ndep` is stored in the file's native units; `ndep_interp!` applies the
conversion to gN/m2/s (mirroring `check_units` + `ndep_interp` in
`ndepStreamMod.F90`).

Ported from the `sdat_ndep` stream state in `ndepStreamMod.F90`.
"""
Base.@kwdef mutable struct NDepStreamData
    active::Bool                 = false
    lon::Vector{Float64}         = Float64[]        # degrees_east, ascending
    lat::Vector{Float64}         = Float64[]        # degrees_north, ascending
    doy::Vector{Float64}         = Float64[]        # stream time as day-of-year
    ndep::Array{Float64,3}       = zeros(0, 0, 0)   # (lon, lat, time), file units
    # true when the file's units are g(N)/m2/yr -> divide by secspday*dayspyr.
    divide_by_secs_per_yr::Bool  = false
    varname::String              = "NDEP_month"
    filename::String             = ""
end

"""
    ndep_stream_init!(s::NDepStreamData, filename; varname="NDEP_month")

Read the N-deposition stream file and check its units.

Mirrors `ndep_init` + `check_units` in `ndepStreamMod.F90`: CTSM accepts exactly
two unit strings and errors on anything else. We do the same rather than
guessing — a silently mis-scaled N input is precisely the class of bug this
module exists to fix.
"""
function ndep_stream_init!(s::NDepStreamData, filename::AbstractString;
                           varname::AbstractString = "NDEP_month")
    isfile(filename) || error("ndep_stream_init!: stream file not found: $filename")

    NCDataset(filename, "r") do ds
        haskey(ds, varname) ||
            error("ndep_stream_init!: variable '$varname' not found in $filename")

        # --- units (check_units) ---
        units = get(ds[varname].attrib, "units", "")
        if units == "g(N)/m2/s"
            s.divide_by_secs_per_yr = false
        elseif units == "g(N)/m2/yr"
            s.divide_by_secs_per_yr = true
        else
            error("ndep_stream_init!: unexpected units for nitrogen deposition: " *
                  "'$units' (expected 'g(N)/m2/s' or 'g(N)/m2/yr')")
        end

        s.lon  = Float64.(Array(ds["lon"][:]))
        s.lat  = Float64.(Array(ds["lat"][:]))
        s.ndep = Float64.(Array(ds[varname][:, :, :]))

        # Stream time -> day-of-year. The CTSM ndep files carry a monthly time
        # axis as "days since <ref>"; taxmode='cycle' means only the position
        # within the year matters. `.var` bypasses CF decoding so we get the raw
        # numeric offsets regardless of the calendar attribute.
        traw = Float64.(Array(ds["time"].var[:]))
        s.doy = _ndep_time_to_doy(traw, get(ds["time"].attrib, "units", ""))
    end

    s.varname  = String(varname)
    s.filename = String(filename)
    s.active   = true
    return nothing
end

# Convert the stream's numeric time axis to day-of-year in [0, dayspyr).
# For "days since <date>" the offsets are day counts from the reference date, so
# the day-of-year is the offset modulo the year length (the ndep streams are a
# single climatological year, which is what taxmode='cycle' assumes).
function _ndep_time_to_doy(traw::Vector{Float64}, units::AbstractString)
    if occursin("days since", units)
        # Offset 0 == the reference date == day-of-year 1.
        return @. mod(traw, 365.0) + 1.0
    elseif occursin("months since", units)
        # Mid-month day-of-year for a 365-day year.
        mid = [15.5, 45.0, 74.5, 105.0, 135.5, 166.0,
               196.5, 227.5, 258.0, 288.5, 319.0, 349.5]
        return [mid[Int(mod(t, 12.0)) + 1] for t in traw]
    else
        error("_ndep_time_to_doy: unsupported time units '$units'")
    end
end

# Bilinear interpolation of a regular lon/lat field at (xlon, xlat).
# Longitude wraps periodically; latitude is clamped at the poles. Weights sum to
# 1 exactly, so a uniform field returns its value exactly.
function _ndep_bilinear(fld::AbstractMatrix{Float64},
                        lon::Vector{Float64}, lat::Vector{Float64},
                        xlon::Real, xlat::Real)
    nlon = length(lon); nlat = length(lat)
    nlon == 1 && nlat == 1 && return fld[1, 1]

    # --- longitude (periodic) ---
    lon0 = lon[1]
    dlon = nlon > 1 ? (lon[2] - lon[1]) : 360.0
    xl   = mod(xlon - lon0, 360.0) / dlon        # fractional index from lon[1]
    i0   = clamp(floor(Int, xl), 0, nlon - 1)
    wlon = xl - i0
    ia   = i0 + 1                                 # 1-based
    ib   = (i0 + 1) % nlon + 1                    # wraps to 1 past the last cell

    # --- latitude (clamped) ---
    j = searchsortedlast(lat, xlat)
    ja = clamp(j, 1, nlat - 1)
    jb = ja + 1
    wlat = (xlat - lat[ja]) / (lat[jb] - lat[ja])
    wlat = clamp(wlat, 0.0, 1.0)

    return (1 - wlon) * (1 - wlat) * fld[ia, ja] +
           wlon       * (1 - wlat) * fld[ib, ja] +
           (1 - wlon) * wlat       * fld[ia, jb] +
           wlon       * wlat       * fld[ib, jb]
end

# Linear-in-time weights on a cyclic (climatological) day-of-year axis.
# Returns (k1, k2, w) so that value = (1-w)*f[k1] + w*f[k2].
function _ndep_time_weights(doy::Vector{Float64}, calday::Real, dayspyr::Real)
    nt = length(doy)
    nt == 1 && return (1, 1, 0.0)

    d = mod(Float64(calday) - 1.0, Float64(dayspyr)) + 1.0
    k1 = findlast(t -> t <= d, doy)
    if k1 === nothing
        # Before the first sample: wrap to the last sample of the previous cycle.
        k1 = nt
        k2 = 1
        span = (doy[1] + dayspyr) - doy[nt]
        w    = (d + dayspyr - doy[nt]) / span
    elseif k1 == nt
        k2   = 1
        span = (doy[1] + dayspyr) - doy[nt]
        w    = (d - doy[nt]) / span
    else
        k2   = k1 + 1
        span = doy[k2] - doy[k1]
        w    = (d - doy[k1]) / span
    end
    return (k1, k2, clamp(w, 0.0, 1.0))
end

"""
    ndep_interp!(a2l, s, grc; calday, dayspyr, bounds_grc)

Time- and space-interpolate the N-deposition stream onto the model gridcells and
fill `a2l.forc_ndep_grc` in **gN/m2/s**.

Mirrors `ndep_interp` in `ndepStreamMod.F90`, including its unit conversion
(`dataptr1d / (secspday * dayspyr)` when the file is in g(N)/m2/yr).

No-op (leaving `forc_ndep_grc` untouched) when the stream is inactive, so the
default non-CN path is unaffected.
"""
function ndep_interp!(a2l::Atm2LndData,
                      s::NDepStreamData,
                      grc::GridcellData;
                      calday::Real,
                      dayspyr::Real,
                      bounds_grc::UnitRange{Int})
    s.active || return nothing
    isempty(bounds_grc) && return nothing

    (k1, k2, w) = _ndep_time_weights(s.doy, calday, dayspyr)
    conv = s.divide_by_secs_per_yr ? 1.0 / (SECSPDAY * dayspyr) : 1.0

    f1 = @view s.ndep[:, :, k1]
    f2 = @view s.ndep[:, :, k2]

    for g in bounds_grc
        xlon = grc.londeg[g]
        xlat = grc.latdeg[g]
        v1 = _ndep_bilinear(f1, s.lon, s.lat, xlon, xlat)
        v2 = _ndep_bilinear(f2, s.lon, s.lat, xlon, xlat)
        a2l.forc_ndep_grc[g] = ((1 - w) * v1 + w * v2) * conv
    end

    return nothing
end
