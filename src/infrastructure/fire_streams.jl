# ==========================================================================
# Ported from: src/cpl/share_esmf/FireDataBaseType.F90
#   (hdm_init / hdm_interp / lnfm_init / lnfm_interp / FireInterp)
#
# The two data streams the Li-family fire methods need:
#   * lnfm — lightning flash density  (counts/km2/hr) -> fire_li2014.forc_lnfm (gridcell)
#   * hdm  — human population density (counts/km2)    -> fire_li2014.forc_hdm  (gridcell)
# Without them the Li ignition term
#   ig = (lh + forc_lnfm[g]/(...)*ignition_efficiency) * (1-fs) * (1-cropf)
# has no lightning source at all, so this module is what makes the ported fire
# chain runnable.
#
# FIDELITY NOTE (mirrors the one on infrastructure/ndep_stream.jl — read it too)
# ----------------------------------------------------------------------------
# CTSM regrids both streams with ESMF (`popdensmapalgo`/`lightngmapalgo` =
# 'bilinear') and time-interpolates with shr_strdata:
#   lnfm : lightng_tintalgo = 'linear' , taxmode = 'cycle'
#   hdm  : popdens_tintalgo = 'nearest', taxmode = 'extend'
# This module reimplements *the same scheme* — bilinear in space (reusing the
# ndep bilinear kernel: same regular lon/lat weights), linear-in-time-cyclic for
# lnfm and nearest-in-time for hdm — but it is NOT bit-validated against ESMF's
# regridder for a general (spatially varying) field, because that needs the ESMF
# weight file.
#
# What the Bow fire parity reference DOES rest on: the reference streams
# (scripts/validation/make_fire_streams.jl) are UNIFORM IN SPACE and CONSTANT IN
# TIME. Bilinear weights sum to 1 and the temporal weights sum to 1, so BOTH
# codes recover the stream value exactly — Fortran's forc_lnfm/forc_hdm are known
# analytically. That isolates the fire PHYSICS from interpolation fidelity.
#
# Consequence: a real (spatially varying) NASA-LIS lnfm or HYDE hdm file is
# handled physically correctly here, but a *bit-for-bit* match with CTSM on such
# a file is not claimed.
#
# Public functions:
#   fire_stream_init!    — read the lnfm/hdm stream files (hdm_init + lnfm_init)
#   fire_stream_interp!  — time+space interpolate onto the gridcells (FireInterp)
# ==========================================================================

"""
    FireStreamField

One raw fire stream (lnfm or hdm) read from its NetCDF file, in the file's
native units (CTSM applies no unit conversion to either — the Li equations
consume counts/km2/hr and counts/km2 directly).

Ported from the `sdat_lnfm` / `sdat_hdm` stream states in `FireDataBaseType.F90`.
"""
Base.@kwdef mutable struct FireStreamField
    active::Bool           = false
    lon::Vector{Float64}   = Float64[]        # degrees_east, ascending
    lat::Vector{Float64}   = Float64[]        # degrees_north, ascending
    doy::Vector{Float64}   = Float64[]        # stream time as day-of-year
    tabs::Vector{Float64}  = Float64[]        # stream time, noleap days since 0000-01-01
    # stream_year_first / _last / model_year_align (hdm: taxmode='extend'). 0 = no
    # year window: nearest by day-of-year over all records (synthetic one-year streams).
    year_first::Int        = 0
    year_last::Int         = 0
    year_align::Int        = 0
    fld::Array{Float64,3}  = zeros(0, 0, 0)   # (lon, lat, time), file units
    varname::String        = ""
    filename::String       = ""
end

"""
    FireStreamData

The pair of fire streams (`lnfm`, `hdm`). Inactive by default, so a model that
never supplies the files behaves exactly as before (the fire path itself is
additionally gated on `cnfire_method !== :nofire`).
"""
Base.@kwdef mutable struct FireStreamData
    lnfm::FireStreamField = FireStreamField()
    hdm::FireStreamField  = FireStreamField()
end

# Read one stream file into a FireStreamField. Same file layout as the ndep
# stream (lon, lat, time + a 3-D (lon,lat,time) field), so the ndep time-axis
# helper is reused verbatim.
function _fire_stream_read!(f::FireStreamField, filename::AbstractString,
                            varname::AbstractString, expected_units::AbstractString)
    isfile(filename) || error("fire_stream_init!: stream file not found: $filename")

    NCDataset(filename, "r") do ds
        haskey(ds, varname) ||
            error("fire_stream_init!: variable '$varname' not found in $filename")

        units = get(ds[varname].attrib, "units", "")
        # CTSM does not convert units for these two streams — it consumes the
        # stream values as-is. We therefore do not convert either, but we DO warn
        # on an unexpected unit string rather than silently mis-scaling ignition.
        if !isempty(units) && units != expected_units
            @warn "fire_stream_init!: unexpected units for '$varname'" units expected_units file=filename
        end

        f.lon = Float64.(Array(ds["lon"][:]))
        f.lat = Float64.(Array(ds["lat"][:]))
        f.fld = Float64.(Array(ds[varname][:, :, :]))

        # `.var` bypasses CF decoding so we get the raw numeric offsets.
        traw  = Float64.(Array(ds["time"].var[:]))
        tunits = get(ds["time"].attrib, "units", "")
        f.doy  = _ndep_time_to_doy(traw, tunits)
        f.tabs = _fire_time_to_abs(traw, tunits)
    end

    f.varname  = String(varname)
    f.filename = String(filename)
    f.active   = true
    return nothing
end

# Stream time axis as noleap days since 0000-01-01 (both CTSM fire streams are on the
# noleap calendar), so a record's year is fld(t, 365). Used for the hdm year window.
function _fire_time_to_abs(traw::Vector{Float64}, units::AbstractString)
    m = match(r"days since\s+(-?\d+)-(\d+)-(\d+)", units)
    m === nothing && return Float64[]
    y0, mo, d0 = parse(Int, m[1]), parse(Int, m[2]), parse(Int, m[3])
    cum = (0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334)
    return @. 365.0 * y0 + cum[mo] + (d0 - 1) + traw
end

"""
    fire_stream_init!(s::FireStreamData; flnfm="", fhdm="",
                      lnfm_varname="lnfm", hdm_varname="hdm",
                      stream_year_first_popdens=0, stream_year_last_popdens=0,
                      model_year_align_popdens=0)

Read the lightning and population-density stream files.

Mirrors `lnfm_init` + `hdm_init` in `FireDataBaseType.F90` (which are called from
`BaseFireInit` only when `need_lightning_and_popdens()` is true — i.e. for the
Li-family methods). An empty filename leaves that stream inactive.

`stream_year_first_popdens` / `stream_year_last_popdens` / `model_year_align_popdens`
are the `popd_streams` namelist window (`hdm_init`, FireDataBaseType.F90). When set
(> 0), only records dated in [first-01-01, (last+1)-01-01) are used and the model
year maps to data year `first + year - align` (CDEPS `shr_stream_findBounds`,
taxmode 'extend'), as needed for the multi-year HYDE file (1850-2016: CTSM runs
use one year of it). Left at 0, hdm is nearest by day of year over every record,
which suits a single-year stream only.
"""
function fire_stream_init!(s::FireStreamData;
                           flnfm::AbstractString = "",
                           fhdm::AbstractString  = "",
                           lnfm_varname::AbstractString = "lnfm",
                           hdm_varname::AbstractString  = "hdm",
                           stream_year_first_popdens::Int = 0,
                           stream_year_last_popdens::Int  = 0,
                           model_year_align_popdens::Int  = 0)
    isempty(flnfm) || _fire_stream_read!(s.lnfm, flnfm, lnfm_varname, "counts/km^2/hr")
    isempty(fhdm)  || _fire_stream_read!(s.hdm,  fhdm,  hdm_varname,  "counts/km^2")
    if !isempty(fhdm) && stream_year_first_popdens > 0
        last_ = max(stream_year_last_popdens, stream_year_first_popdens)
        align = model_year_align_popdens > 0 ? model_year_align_popdens : stream_year_first_popdens
        h = s.hdm
        isempty(h.tabs) && error("fire_stream_init!: hdm year window needs a 'days since' time axis")
        any(t -> stream_year_first_popdens <= fld(t, 365.0) <= last_, h.tabs) ||
            error("fire_stream_init!: no hdm record in years $(stream_year_first_popdens)-$(last_)")
        h.year_first, h.year_last, h.year_align = stream_year_first_popdens, last_, align
    end
    return nothing
end

# hdm record for model (year, calday) under CDEPS taxmode='extend', tintalgo='nearest'
# (dshr_stream_mod.F90 shr_stream_findBounds; dshr_tinterp_mod.F90:108-119). Valid
# records are those dated in [year_first-01-01, (year_last+1)-01-01); the model year
# maps to data year year_first + year - year_align; before the first valid record
# (LVD) or at/after the last (GVD) the end record is held; in between, the nearer
# bound wins and a tie goes to the lower bound (abs(sint1) <= abs(sint2)).
function _fire_extend_nearest_index(f::FireStreamField, year::Integer, calday::Real)
    lo = 365.0 * f.year_first
    hi = 365.0 * (f.year_last + 1)
    lvd = findfirst(t -> t >= lo, f.tabs)
    gvd = findlast(t -> t < hi, f.tabs)
    dyear = f.year_first + Int(year) - f.year_align
    tm = 365.0 * dyear + (Float64(calday) - 1.0)
    tm < f.tabs[lvd]  && return lvd
    tm >= f.tabs[gvd] && return gvd
    k1 = findlast(t -> t <= tm, f.tabs)
    k2 = k1 + 1
    return abs(tm - f.tabs[k1]) <= abs(f.tabs[k2] - tm) ? k1 : k2
end

# Nearest-in-time sample index on a cyclic day-of-year axis (hdm's
# popdens_tintalgo='nearest'). With taxmode='extend' and a single-year stream the
# nearest sample is simply the closest day-of-year, wrapping across the year end.
function _fire_nearest_time_index(doy::Vector{Float64}, calday::Real, dayspyr::Real)
    nt = length(doy)
    nt <= 1 && return 1
    d = mod(Float64(calday) - 1.0, Float64(dayspyr)) + 1.0
    best = 1
    bestdist = Inf
    for k in 1:nt
        raw = abs(d - doy[k])
        dist = min(raw, Float64(dayspyr) - raw)   # cyclic distance
        if dist < bestdist
            bestdist = dist
            best = k
        end
    end
    return best
end

# Space-interpolate one field slice at a gridcell (bilinear, weights sum to 1).
_fire_space_interp(f::FireStreamField, k::Int, xlon::Real, xlat::Real) =
    _ndep_bilinear(@view(f.fld[:, :, k]), f.lon, f.lat, xlon, xlat)

"""
    fire_stream_interp!(fire_li2014, s, grc; calday, dayspyr, bounds_grc)

Time- and space-interpolate the lightning / population-density streams onto the
model gridcells, filling `fire_li2014.forc_lnfm` (counts/km2/hr) and
`fire_li2014.forc_hdm` (counts/km2).

Mirrors `FireInterp` (`lnfm_interp` + `hdm_interp`) in `FireDataBaseType.F90`:
lnfm is linear-in-time on a cyclic axis, hdm is nearest-in-time. No unit
conversion is applied (CTSM applies none).

No-op for a stream that was never read, leaving its `forc_*` array untouched.
"""
function fire_stream_interp!(fire_li2014::CNFireLi2014Data,
                             s::FireStreamData,
                             grc::GridcellData;
                             calday::Real,
                             dayspyr::Real,
                             bounds_grc::UnitRange{Int},
                             year::Union{Nothing,Integer} = nothing)
    _fire_stream_interp_into!(fire_li2014.forc_lnfm, fire_li2014.forc_hdm, s, grc;
                              calday = calday, dayspyr = dayspyr, bounds_grc = bounds_grc,
                              year = year)
end

# FireInterp core shared by the CN (Li) and FATES (FATESFireDataMod) hosts: fill the
# gridcell lightning [counts/km2/hr] and population-density [counts/km2] arrays.
function _fire_stream_interp_into!(forc_lnfm::AbstractVector, forc_hdm::AbstractVector,
                                   s::FireStreamData, grc::GridcellData;
                                   calday::Real, dayspyr::Real,
                                   bounds_grc::UnitRange{Int},
                                   year::Union{Nothing,Integer} = nothing)
    isempty(bounds_grc) && return nothing

    if s.lnfm.active
        (k1, k2, w) = _ndep_time_weights(s.lnfm.doy, calday, dayspyr)
        for g in bounds_grc
            xlon = grc.londeg[g]
            xlat = grc.latdeg[g]
            v1 = _fire_space_interp(s.lnfm, k1, xlon, xlat)
            v2 = _fire_space_interp(s.lnfm, k2, xlon, xlat)
            forc_lnfm[g] = (1 - w) * v1 + w * v2
        end
    end

    if s.hdm.active
        k = (s.hdm.year_first > 0 && year !== nothing) ?
            _fire_extend_nearest_index(s.hdm, year, calday) :
            _fire_nearest_time_index(s.hdm.doy, calday, dayspyr)
        for g in bounds_grc
            forc_hdm[g] = _fire_space_interp(s.hdm, k, grc.londeg[g], grc.latdeg[g])
        end
    end

    return nothing
end

# ---------------------------------------------------------------------------
# FATES fire data (FATESFireDataMod.F90) — the host side of SPITFIRE modes 2-5
# ---------------------------------------------------------------------------
"""
    FATESFireData

Host state for SPITFIRE `fates_spitfire_mode > 1` (`FATESFireDataMod.F90`):
the lightning and population-density streams interpolated to the gridcell
(`FireInterp`), and the 24-hour running mean of lightning (`lnfm24`, a `runmean`
accumulator with period -1 day and init_value 0, updated every step in
`UpdateAccVars`). `dynamics_driv` hands `lnfm24 * 24` [#/km2/day] to
`bc_in%lightning24` and `forc_hdm` to `bc_in%pop_density` (modes >= 3).
Zero-length (inactive) unless a FATES run with `fates_spitfire_mode > 1` is
initialized with stream files.
"""
Base.@kwdef mutable struct FATESFireData
    stream::FireStreamData          = FireStreamData()
    forc_lnfm::Vector{Float64}      = Float64[]   # lightning [counts/km2/hr]
    forc_hdm::Vector{Float64}       = Float64[]   # population density [counts/km2]
    lnfm24::Vector{Float64}         = Float64[]   # 24-h running mean of forc_lnfm
    lnfm24_nsteps::Vector{Float64}  = Float64[]   # accumulMod per-point nsteps
    # Seconds added to the driver's (jday, secs) before the streams are advanced.
    # Fortran advances them at get_curr_date (lnfm_interp/hdm_interp,
    # FireDataBaseType.F90), the END of the step, which is what clm_run passes. The
    # FATES parity harness passes the START of the step (scripts/fates_fortran_parity.jl)
    # and sets this to dtime. That matters for the 3-hourly LIS lightning
    # climatology: a one-step lag moves the day-mean forc_lnfm by 1.5 % on average.
    tod_offset::Int                 = 0
end

fates_fire_data_active(d::FATESFireData) = !isempty(d.lnfm24)

"""
    fates_fire_data_init!(d, ng; flnfm="", fhdm="")

`FATESFireDataMod` Init + InitAccBuffer/InitAccVars: read the streams and allocate
the gridcell arrays (lnfm24 = its init_value 0).
"""
function fates_fire_data_init!(d::FATESFireData, ng::Int;
                               flnfm::AbstractString = "", fhdm::AbstractString = "",
                               stream_kwargs...)
    fire_stream_init!(d.stream; flnfm = flnfm, fhdm = fhdm, stream_kwargs...)
    d.forc_lnfm     = zeros(ng)
    d.forc_hdm      = zeros(ng)
    d.lnfm24        = zeros(ng)
    d.lnfm24_nsteps = zeros(ng)
    return d
end

"`FATESFireDataMod::FireInterp` (called from clm_fates%InterpFileInputs each step)."
fates_fire_interp!(d::FATESFireData, grc::GridcellData; calday::Real, dayspyr::Real,
                   bounds_grc::UnitRange{Int}, year::Union{Nothing,Integer} = nothing) =
    _fire_stream_interp_into!(d.forc_lnfm, d.forc_hdm, d.stream, grc;
                              calday = calday, dayspyr = dayspyr, bounds_grc = bounds_grc,
                              year = year)

"""
`FATESFireDataMod::UpdateAccVars`: lnfm24 runmean (1 day) of forc_lnfm over the
gridcells, with accumulMod's per-point step counter.
"""
function fates_fire_update_acc_vars!(d::FATESFireData, bounds_grc::UnitRange{Int};
                                     dtime::Int)
    w1 = Float64(accum_window_steps(1, dtime))
    for g in bounds_grc
        n = min(d.lnfm24_nsteps[g] + 1.0, w1)
        d.lnfm24_nsteps[g] = n
        d.lnfm24[g] = ((n - 1.0) * d.lnfm24[g] + d.forc_lnfm[g]) / n
    end
    return nothing
end
