# ==========================================================================
# NetCDF forcing data reader
# No Fortran equivalent (CESM coupler provides forcings)
#
# Public functions:
#   forcing_reader_init!   — Open forcing file
#   read_forcing_step!     — Read one timestep of forcing data
#   forcing_reader_close!  — Close forcing file
# ==========================================================================

"""
    ForcingReader

Stateful reader for NetCDF atmospheric forcing files. Tracks the current
time index and provides step-by-step reading of forcing variables.

Variable mapping from forcing file → Atm2LndData:
  TBOT     → forc_t_not_downscaled_grc
  PSRF     → forc_pbot_not_downscaled_grc
  WIND     → forc_wind_grc, forc_u_grc
  FLDS     → forc_lwrad_not_downscaled_grc
  FSDS     → forc_solad_not_downscaled_grc, forc_solai_grc
  PRECTmms → forc_rain_not_downscaled_grc / forc_snow_not_downscaled_grc
  QBOT     → forc_vp_grc (via q → e conversion)
  RH       → forc_vp_grc (via RH → e conversion, alternative)
"""
Base.@kwdef mutable struct ForcingReader
    ds::Union{NCDataset, Nothing} = nothing
    time_index::Int = 0
    ntimes::Int = 0
    times::Vector{DateTime} = DateTime[]
    filepath::String = ""

    # --- Atm spatial grid (multi-gridcell support) ---
    # When the forcing file carries a 2D (lon,lat) atm grid, these hold the
    # native atm-grid coordinates (degrees). For a single-point forcing file
    # (no horizontal dims) `n_atm` stays 1 and the legacy broadcast path is
    # used (byte-identical to the historical single-gridcell behavior).
    atm_lon::Vector{Float64} = Float64[]   # atm grid longitudes (degrees), length nlon
    atm_lat::Vector{Float64} = Float64[]   # atm grid latitudes  (degrees), length nlat
    nlon::Int = 0
    nlat::Int = 0
    n_atm::Int = 1                         # number of distinct atm points (nlon*nlat, or 1 for a point)
    is_point::Bool = true                  # true => single atm point, broadcast to all land cells

    # atm→land index map: for each land gridcell g, the linear index (1-based,
    # column-major over (lon,lat)) of its mapped atm cell. Built once, cached.
    map_g2atm::Vector{Int} = Int[]
    map_ng::Int = 0                        # ng the map was built for (rebuild if it changes)

    # Linear time-interpolation of the forcing to the model step (datm `tintalgo=
    # linear`). When the forcing is coarser than the model step (e.g. hourly forcing
    # driving a 30-min run), nearest-neighbour selection makes a sub-step temperature
    # error that can flip the rain/snow partition at a marginal event. Off by default
    # (nearest-neighbour, byte-identical to the legacy reader); enabled per-run.
    interp_time::Bool = false
end

"""
    forcing_reader_init!(fr, filepath)

Open a NetCDF forcing file and read the time coordinate.
"""
function forcing_reader_init!(fr::ForcingReader, filepath::String)
    fr.filepath = filepath
    fr.ds = NCDataset(filepath, "r")
    fr.time_index = 0

    # Read time coordinate (may return DateTimeNoLeap from NCDatasets)
    if haskey(fr.ds, "time")
        raw_times = collect(fr.ds["time"][:])
        # Convert any CF-time types to DateTime
        fr.times = map(raw_times) do t
            if t isa DateTime
                t
            else
                # Handle DateTimeNoLeap, DateTimeAllLeap, etc.
                DateTime(Dates.year(t), Dates.month(t), Dates.day(t),
                         Dates.hour(t), Dates.minute(t), Dates.second(t))
            end
        end
        fr.ntimes = length(fr.times)
    else
        fr.ntimes = 0
        fr.times = DateTime[]
    end

    # --- Detect / read the atm spatial grid ---
    _read_atm_grid!(fr)

    # Reset the cached atm→land map (rebuilt lazily on the first read with coords)
    fr.map_g2atm = Int[]
    fr.map_ng = 0

    return nothing
end

"""
    _read_atm_grid!(fr)

Read the atmospheric forcing grid's horizontal coordinates, if present.

Recognized coordinate variable names (case-insensitive aliases, in order of
preference): longitude → `LONGXY`, `lon`, `longitude`, `LON`; latitude →
`LATIXY`, `lat`, `latitude`, `LAT`. Both 1D axes (`lon(nlon)`, `lat(nlat)`) and
2D meshes (`LONGXY(nlon,nlat)`, `LATIXY(nlon,nlat)`, as in CLM domain/surfdata)
are supported; a 2D mesh is reduced to its per-cell coordinate lists.

If no horizontal coordinates are found, the file is treated as a single atm
point (`is_point = true`) and the legacy broadcast path is used.
"""
function _read_atm_grid!(fr::ForcingReader)
    ds = fr.ds
    fr.atm_lon = Float64[]
    fr.atm_lat = Float64[]
    fr.nlon = 0
    fr.nlat = 0
    fr.n_atm = 1
    fr.is_point = true

    _find_coord(names) = begin
        for nm in names
            haskey(ds, nm) && return Array(ds[nm])
        end
        return nothing
    end

    lonv = _find_coord(("LONGXY", "lon", "longitude", "LON", "Longitude"))
    latv = _find_coord(("LATIXY", "lat", "latitude", "LAT", "Latitude"))

    (lonv === nothing || latv === nothing) && return nothing

    if ndims(lonv) == 1 && ndims(latv) == 1
        # 1D axes: full lon×lat tensor grid
        fr.atm_lon = Float64.(vec(lonv))
        fr.atm_lat = Float64.(vec(latv))
        fr.nlon = length(fr.atm_lon)
        fr.nlat = length(fr.atm_lat)
    elseif ndims(lonv) == 2 && ndims(latv) == 2
        # 2D mesh (nlon,nlat): collapse to per-axis coordinate lists. CLM
        # domain/surfdata meshes are regular, so a row/column slice suffices.
        nlon, nlat = size(lonv)
        fr.nlon = nlon
        fr.nlat = nlat
        fr.atm_lon = Float64.(vec(lonv[:, 1]))   # lon varies along dim 1
        fr.atm_lat = Float64.(vec(latv[1, :]))   # lat varies along dim 2
    else
        return nothing
    end

    fr.n_atm = fr.nlon * fr.nlat
    # A 1×1 grid is still effectively a single point → keep the broadcast path
    # (byte-identical to single-gridcell forcing).
    fr.is_point = (fr.n_atm <= 1)

    return nothing
end

"""
    _build_atm2land_map!(fr, gridcell_latdeg, gridcell_londeg)

Build (and cache) the nearest-neighbour atm→land index map: for each land
gridcell, the linear index (column-major over the atm `(lon,lat)` grid) of the
nearest atm cell by great-circle distance. Called once per forcing-grid /
gridcell-set; cached in `fr.map_g2atm`.
"""
function _build_atm2land_map!(fr::ForcingReader,
                              gridcell_latdeg::AbstractVector,
                              gridcell_londeg::AbstractVector)
    ng = length(gridcell_latdeg)
    fr.map_g2atm = Vector{Int}(undef, ng)
    fr.map_ng = ng

    for g in 1:ng
        latg = Float64(gridcell_latdeg[g])
        long = Float64(gridcell_londeg[g])
        fr.map_g2atm[g] = _nearest_atm_index(fr.atm_lon, fr.atm_lat,
                                             fr.nlon, fr.nlat, latg, long)
    end
    return nothing
end

"""
Great-circle (haversine) nearest-neighbour lookup over a tensor `(lon,lat)`
grid. Returns the column-major linear index `ilon + (ilat-1)*nlon`. Longitudes
are compared on the circle so the ±180/360 wrap is handled correctly.
"""
function _nearest_atm_index(atm_lon::AbstractVector, atm_lat::AbstractVector,
                            nlon::Int, nlat::Int, latg::Float64, long::Float64)
    deg2rad = π / 180.0
    φg = latg * deg2rad
    λg = long * deg2rad
    sinφg = sin(φg)
    cosφg = cos(φg)

    best_idx = 1
    best_d = Inf
    for jlat in 1:nlat
        φa = atm_lat[jlat] * deg2rad
        sinφa = sin(φa)
        cosφa = cos(φa)
        for ilon in 1:nlon
            λa = atm_lon[ilon] * deg2rad
            # central angle via spherical law of cosines (monotone in distance)
            c = sinφg * sinφa + cosφg * cosφa * cos(λa - λg)
            c = clamp(c, -1.0, 1.0)
            d = acos(c)
            if d < best_d
                best_d = d
                best_idx = ilon + (jlat - 1) * nlon
            end
        end
    end
    return best_idx
end

"""
    read_forcing_step!(fr, a2l, target_time, ng, nc; gridcell_latdeg, gridcell_londeg)

Read the forcing timestep closest to `target_time` and populate `a2l`.

For a single-point forcing file (no horizontal grid) each variable's scalar
value is broadcast to all `ng` land gridcells — byte-identical to the historical
single-gridcell behavior. For a 2D `(lon,lat)` atm grid, each land gridcell `g`
reads the forcing of its nearest atm cell, using the atm→land map built from the
gridcell coordinates (`gridcell_latdeg`/`gridcell_londeg`, in degrees). The map
is built once and cached on `fr`. If the coordinates are omitted, the legacy
broadcast (atm cell 1) is used.
"""
function read_forcing_step!(fr::ForcingReader, a2l::Atm2LndData,
                            target_time::DateTime, ng::Int, nc::Int;
                            gridcell_latdeg::Union{AbstractVector, Nothing} = nothing,
                            gridcell_londeg::Union{AbstractVector, Nothing} = nothing,
                            dtime::Int = 1800)
    ds = fr.ds
    ds === nothing && error("ForcingReader not initialized")

    # Time selection. Default: nearest-neighbour (ti0==ti1, w=0 → byte-identical to
    # the legacy reader). With interp_time: linearly interpolate between the two
    # bracketing forcing times (datm tintalgo='linear'), which matters when the
    # forcing is coarser than the model step.
    local ti, ti0, ti1, tw
    if fr.interp_time
        (ti0, ti1, tw) = _find_bracket_time(fr.times, target_time)
        ti = ti0
    else
        ti = _find_closest_time(fr.times, target_time, fr.time_index)
        ti0 = ti; ti1 = ti; tw = 0.0
    end
    fr.time_index = ti

    # Decide whether to spatially map. We map only when the file has a real 2D
    # atm grid AND gridcell coordinates were supplied; otherwise we broadcast a
    # single atm point (the legacy, byte-identical path).
    do_map = (!fr.is_point) && (gridcell_latdeg !== nothing) &&
             (gridcell_londeg !== nothing) && fr.n_atm > 1

    if do_map
        # (Re)build the atm→land index map if needed (static; built once).
        if isempty(fr.map_g2atm) || fr.map_ng != ng
            _build_atm2land_map!(fr, gridcell_latdeg, gridcell_londeg)
        end
    end

    # Read a variable into a per-gridcell vector (length ng).
    #   - point / no-map path: the single broadcast scalar is filled into all g,
    #     identical to the historical reader (data[1,ti] / data[1,1,ti]).
    #   - mapped path: each g gets the value at its mapped atm cell.
    # The variable's data array is read from disk exactly once per call.
    function _read_var(varname::String, default::Float64; interp::Bool = true)
        # Per-field interpolation control: datm linearly interpolates the STATE
        # fields (T, P, wind, humidity, LW) but holds PRECIP constant over the
        # source interval and uses coszen for solar — so precip/solar pass interp=
        # false and keep nearest-neighbour even when interp_time is on. (Smearing
        # precip linearly destroys the snowfall-event timing → wrecks snow_depth.)
        w_eff = interp ? tw : 0.0
        out = Vector{Float64}(undef, ng)
        if haskey(ds, varname)
            data = Array(ds[varname])
            nd = ndims(data)
            # Value at a single time index for the broadcast (no-map) path.
            scalar_at(tx) = if nd == 1
                Float64(data[tx])
            elseif nd == 2
                sz = size(data)
                sz[end] >= tx ? Float64(data[1, tx]) : Float64(data[tx, 1])
            elseif nd == 3
                Float64(data[1, 1, tx])
            else
                default
            end
            # Value at a single time index for a mapped atm cell.
            mapped_at(aidx, tx) = if nd == 1
                Float64(data[tx])
            elseif nd == 2
                Float64(data[aidx, tx])
            elseif nd == 3
                ilon = ((aidx - 1) % fr.nlon) + 1
                jlat = ((aidx - 1) ÷ fr.nlon) + 1
                Float64(data[ilon, jlat, tx])
            else
                default
            end
            if !do_map
                # w_eff==0 → byte-identical to the legacy nearest read.
                val = w_eff == 0.0 ? scalar_at(ti0) :
                      (1.0 - w_eff) * scalar_at(ti0) + w_eff * scalar_at(ti1)
                @inbounds for g in 1:ng
                    out[g] = val
                end
            else
                @inbounds for g in 1:ng
                    aidx = fr.map_g2atm[g]
                    out[g] = w_eff == 0.0 ? mapped_at(aidx, ti0) :
                             (1.0 - w_eff) * mapped_at(aidx, ti0) + w_eff * mapped_at(aidx, ti1)
                end
            end
        else
            @inbounds for g in 1:ng
                out[g] = default
            end
        end
        return out
    end

    # --- Read each variable (per-gridcell vectors) and populate forcings ---

    # Temperature [K]
    tbot = _read_var("TBOT", 270.0)
    for g in 1:ng
        a2l.forc_t_not_downscaled_grc[g] = tbot[g]
    end

    # Surface pressure [Pa]
    psrf = _read_var("PSRF", 85000.0)
    for g in 1:ng
        a2l.forc_pbot_not_downscaled_grc[g] = psrf[g]
    end

    # Wind speed [m/s]
    wind = _read_var("WIND", 3.0)
    for g in 1:ng
        a2l.forc_wind_grc[g] = wind[g]
        # Decompose the scalar forcing wind speed equally into east/north
        # components, matching the CLMNCEP datm (datm_datamode_clmncep_mod.F90:435
        # Sa_u = strm_wind/sqrt(2); Sa_v = Sa_u). The wind SPEED sqrt(u²+v²)=wind
        # is unchanged (so aerodynamics/physics are identical), but the taux/tauy
        # momentum-flux diagnostics need this split — the old u=wind,v=0 put all
        # stress in taux (×sqrt(2) too large) and left tauy=0.
        a2l.forc_u_grc[g] = wind[g] / sqrt(2.0)
        a2l.forc_v_grc[g] = wind[g] / sqrt(2.0)
    end

    # Longwave radiation [W/m2]
    flds = _read_var("FLDS", 250.0)
    for g in 1:ng
        a2l.forc_lwrad_not_downscaled_grc[g] = flds[g]
    end

    # Shortwave radiation [W/m2] — split into VIS/NIR (50/50) then
    # direct/diffuse using CAM-derived polynomial (Fortran DATM CLMNCEP)
    # FSDS uses datm tintalgo='coszen': take the LOWER-bound interval value and
    # scale by cosz(t)/avg_cosz([LB,UB]). interp=false reads the floor (ti0=LB) value.
    fsds = _read_var("FSDS", 0.0; interp=false)
    if fr.interp_time && ti1 > ti0 && gridcell_latdeg !== nothing && gridcell_londeg !== nothing
        cday = _calday_of(target_time)
        @inbounds for g in 1:ng
            latr = Float64(gridcell_latdeg[g]) * _DEG2RAD
            lonr = Float64(gridcell_londeg[g]) * _DEG2RAD
            czm = _cosz_at(cday, latr, lonr)
            if czm > _SOLZENMIN
                avgcz = _avg_cosz(fr.times[ti0], fr.times[ti1], latr, lonr, dtime)
                fsds[g] = fsds[g] * czm / avgcz
            else
                fsds[g] = 0.0
            end
        end
    end
    for g in 1:ng
        fsds_g = fsds[g]
        swndr = fsds_g * 0.50  # NIR half
        ratio_nir = clamp(0.29548 + 0.00504*swndr - 1.4957e-5*swndr^2 + 1.4881e-8*swndr^3, 0.01, 0.99)
        swvdr = fsds_g * 0.50  # VIS half
        ratio_vis = clamp(0.17639 + 0.00380*swvdr - 9.0039e-6*swvdr^2 + 8.1351e-9*swvdr^3, 0.01, 0.99)
        a2l.forc_solad_not_downscaled_grc[g, 1] = ratio_vis * swvdr        # VIS direct
        a2l.forc_solad_not_downscaled_grc[g, 2] = ratio_nir * swndr        # NIR direct
        a2l.forc_solai_grc[g, 1] = (1.0 - ratio_vis) * swvdr               # VIS diffuse
        a2l.forc_solai_grc[g, 2] = (1.0 - ratio_nir) * swndr               # NIR diffuse
        a2l.forc_solar_not_downscaled_grc[g] = fsds_g
    end

    # Precipitation [mm/s] — total, partition by temperature. DATM holds
    # precipitation constant over the source interval, and the corresponding
    # RAIN/SNOW forcing split uses the lower-bound interval temperature rather
    # than the linearly interpolated state temperature used for surface physics.
    precip = _read_var("PRECTmms", 0.0; interp=false)  # datm holds precip constant over the interval
    tbot_precip = _read_var("TBOT", 270.0; interp=false)
    # Partition precipitation using a linear ramp matching Fortran DATM
    # (shr_precip_mod.F90): frac_rain = (T - TFRZ) * 0.5, clamped to [0,1]
    # All snow at T <= 0°C, all rain at T >= +2°C
    for g in 1:ng
        precip_g = precip[g] < 0.0 ? 0.0 : precip[g]
        frac_rain = clamp((tbot_precip[g] - TFRZ) * 0.5, 0.0, 1.0)
        a2l.forc_rain_not_downscaled_grc[g] = precip_g * frac_rain
        a2l.forc_snow_not_downscaled_grc[g] = precip_g * (1.0 - frac_rain)
    end

    # Specific humidity → vapor pressure
    # Try QBOT first, then RH
    if haskey(ds, "QBOT")
        qbot = _read_var("QBOT", 0.003)
        for g in 1:ng
            # e = q * p / (0.622 + 0.378 * q)
            a2l.forc_vp_grc[g] = qbot[g] * psrf[g] / (0.622 + 0.378 * qbot[g])
            # Fortran carries Sa_shum straight through to forc_q_not_downscaled_grc.
            # This port was leaving that field at its init value, which made
            # atm2lnd_update_rh! (its ONLY consumer) compute forc_rh_grc = 0 for every
            # run — i.e. the RH30 accumulator and the Li fire RH terms saw a
            # permanently bone-dry atmosphere. Set it from QBOT.
            a2l.forc_q_not_downscaled_grc[g] = qbot[g]
        end
    elseif haskey(ds, "RH")
        rh = _read_var("RH", 70.0)
        for g in 1:ng
            rh_g = rh[g] / 100.0  # convert % to fraction
            # Saturation vapor pressure (Tetens formula)
            tc = tbot[g] - TFRZ
            esat = 611.0 * exp(17.27 * tc / (tc + 237.3))
            a2l.forc_vp_grc[g] = rh_g * esat
            # q from e (inverse of the QBOT branch above), so forc_rh_grc is
            # reconstructible on the RH-forced path too.
            e = a2l.forc_vp_grc[g]
            a2l.forc_q_not_downscaled_grc[g] = 0.622 * e / max(psrf[g] - 0.378 * e, 1.0)
        end
    end

    # Potential temperature. The Fortran datm for this observed single-point
    # forcing delivers forc_th = forc_t (Sa_ptem == TBOT; no reference-pressure
    # adjustment), so thv/thvstar in the surface layer match. Applying the
    # standard (100000/pbot)^kappa factor here (pbot≈79000 at this altitude)
    # inflated forc_th/thv by ~7% (307 vs 287 K) and biased the Monin-Obukhov
    # stability solve.
    for g in 1:ng
        t = a2l.forc_t_not_downscaled_grc[g]
        a2l.forc_th_not_downscaled_grc[g] = t
    end

    # Density from equation of state
    for g in 1:ng
        pbot = a2l.forc_pbot_not_downscaled_grc[g]
        t = a2l.forc_t_not_downscaled_grc[g]
        vp = a2l.forc_vp_grc[g]
        # ρ = (p - 0.378*e) / (Rd * T)
        a2l.forc_rho_not_downscaled_grc[g] = (pbot - 0.378 * vp) / (RAIR * t)
    end

    # Reference heights (set defaults if zero)
    for g in 1:ng
        if a2l.forc_hgt_grc[g] <= 0.0
            a2l.forc_hgt_grc[g] = 30.0  # default 30m
        end
        if a2l.forc_hgt_u_grc[g] <= 0.0
            a2l.forc_hgt_u_grc[g] = a2l.forc_hgt_grc[g]
        end
        if a2l.forc_hgt_t_grc[g] <= 0.0
            a2l.forc_hgt_t_grc[g] = a2l.forc_hgt_grc[g]
        end
        if a2l.forc_hgt_q_grc[g] <= 0.0
            a2l.forc_hgt_q_grc[g] = a2l.forc_hgt_grc[g]
        end
    end

    # CO2 / O2 defaults. CLM forms the partial pressures from the molar fractions
    # times the ACTUAL surface pressure (forc_pco2 = co2_ppmv*1e-6*pbot,
    # forc_po2 = 0.209*pbot), so they must scale with pbot — a fixed sea-level value
    # over-estimates both at high elevation (Bow pbot≈79 kPa → 40 Pa would be 506 ppm,
    # vs Fortran's 367 ppm = 29 Pa). co2_ppmv = 367 matches the I2000 (year-2000) CO2.
    # Seed here from the NOT-downscaled grc pbot; downscale_forcings! then rescales
    # these partial pressures to the elevation-corrected (downscaled) column pbot the
    # photosynthesis solve actually uses for `cair`/cs (matching Fortran, which reads
    # forc_pco2 built from the same surface pbot). Recompute each step (NOT gated on
    # <=0) so the value tracks pbot and the downscale rescale stays non-compounding —
    # the old <=0 freeze pinned forc_pco2 at the not-downscaled step-1 value, leaving
    # cair ~2.7% low at elevation-corrected columns → stomata too open → +2% transp.
    for g in 1:ng
        a2l.forc_pco2_grc[g] = 367.0e-6 * a2l.forc_pbot_not_downscaled_grc[g]
        a2l.forc_po2_grc[g]  = 0.209    * a2l.forc_pbot_not_downscaled_grc[g]
    end

    return nothing
end

"""
    forcing_reader_close!(fr)

Close the forcing NetCDF file.
"""
function forcing_reader_close!(fr::ForcingReader)
    if fr.ds !== nothing
        close(fr.ds)
        fr.ds = nothing
    end
    return nothing
end

# ---- Internal helpers ----

"""
Find the time index closest to target_time, starting search from hint.
"""
# Bracket `target` between two forcing times and return (i0, i1, w) such that the
# linearly-interpolated value is (1-w)*v[i0] + w*v[i1], with w∈[0,1]. Clamps to the
# endpoints outside the forcing span. Mirrors datm `tintalgo='linear'`.
function _find_bracket_time(times::Vector{DateTime}, target::DateTime)
    n = length(times)
    n == 0 && return (1, 1, 0.0)
    n == 1 && return (1, 1, 0.0)
    if target <= times[1]
        return (1, 1, 0.0)
    elseif target >= times[n]
        return (n, n, 0.0)
    end
    # LB = largest forcing time <= target, UB = next. Find the first time STRICTLY
    # after target (UB); LB is the one before it. The strict `>` is essential: a model
    # step landing exactly on a forcing time tk must bracket [tk, tk+1] (LB=tk), not
    # [tk-1, tk] — otherwise the coszen LB-value pick is off by one interval.
    i1 = n
    @inbounds for i in 2:n
        if times[i] > target
            i1 = i
            break
        end
    end
    i0 = i1 - 1
    span = Dates.value(times[i1] - times[i0])
    w = span > 0 ? Dates.value(target - times[i0]) / span : 0.0
    return (i0, i1, clamp(w, 0.0, 1.0))
end

# --- coszen solar time-interpolation (datm tintalgo='coszen') ---
# Solar is interpolated by the cosine of the solar zenith angle, not linearly:
#   FSDS(t) = FSDS_LB · cosz(t) / avg_cosz([LB,UB])   (0 when cosz ≤ solZenMin)
# This shapes the diurnal cycle while conserving the interval-mean flux, matching
# CDEPS dshr_strdata_mod / dshr_tinterp_mod. Ported faithfully (shr_orb_cosz).
const _SOLZENMIN = 0.001        # min solar zenith cosine (dshr_tinterp_mod)
const _DEG2RAD   = π / 180.0

# Calendar day (day-of-year + day fraction), matching get_curr_calday.
_calday_of(dt::DateTime) =
    Float64(Dates.dayofyear(dt)) +
    (Dates.hour(dt) * 3600 + Dates.minute(dt) * 60 + Dates.second(dt)) / SECSPDAY

# Cosine of the solar zenith angle (shr_orb_cosz), lat/lon in radians.
function _cosz_at(calday::Float64, latr::Float64, lonr::Float64)
    (declin, _) = compute_orbital(calday)
    return sin(latr) * sin(declin) -
           cos(latr) * cos(declin) * cos((calday - floor(calday)) * 2.0 * π + lonr)
end

# Time-average of max(solZenMin, cosz) over [t0, t1), stepping at the model dt
# (adjusted so the interval divides evenly). Mirrors shr_tInterp_getAvgCosz.
function _avg_cosz(t0::DateTime, t1::DateTime, latr::Float64, lonr::Float64, dtime::Int)
    dtsec = Dates.value(Dates.Second(t1 - t0))
    dtsec <= 0 && return max(_SOLZENMIN, _cosz_at(_calday_of(t0), latr, lonr))
    ldt = dtime
    if dtsec % ldt != 0
        ldt = dtsec ÷ (dtsec ÷ ldt + 1)
    end
    total = 0.0
    n = 0
    t = t0
    @inbounds while t < t1
        total += max(_SOLZENMIN, _cosz_at(_calday_of(t), latr, lonr)) * ldt
        n += ldt
        t += Dates.Second(ldt)
    end
    return n > 0 ? total / n : _SOLZENMIN
end

function _find_closest_time(times::Vector{DateTime}, target::DateTime, hint::Int)
    if isempty(times)
        return 1
    end

    best_idx = max(1, hint)
    best_diff = abs(Dates.value(times[best_idx] - target))

    for i in eachindex(times)
        d = abs(Dates.value(times[i] - target))
        if d < best_diff
            best_diff = d
            best_idx = i
        end
    end

    return best_idx
end
