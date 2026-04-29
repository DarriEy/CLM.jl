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

    return nothing
end

"""
    read_forcing_step!(fr, a2l, target_time, ng, nc)

Read the forcing timestep closest to `target_time` and populate `a2l`.
"""
function read_forcing_step!(fr::ForcingReader, a2l::Atm2LndData,
                            target_time::DateTime, ng::Int, nc::Int)
    ds = fr.ds
    ds === nothing && error("ForcingReader not initialized")

    # Find closest time index
    ti = _find_closest_time(fr.times, target_time, fr.time_index)
    fr.time_index = ti

    # Helper to read a variable (handles spatial dims gracefully)
    function _read_var(varname::String, default::Float64)
        if haskey(ds, varname)
            data = Array(ds[varname])
            if ndims(data) == 1  # (time,) only
                return data[ti]
            elseif ndims(data) == 2  # (spatial, time) or (time, spatial)
                # Try to get scalar for single gridcell
                sz = size(data)
                if sz[end] >= ti
                    return Float64(data[1, ti])
                else
                    return Float64(data[ti, 1])
                end
            elseif ndims(data) == 3  # (lon, lat, time)
                return Float64(data[1, 1, ti])
            end
        end
        return default
    end

    # --- Read each variable and populate gridcell-level forcings ---

    # Temperature [K]
    tbot = _read_var("TBOT", 270.0)
    for g in 1:ng
        a2l.forc_t_not_downscaled_grc[g] = tbot
    end

    # Surface pressure [Pa]
    psrf = _read_var("PSRF", 85000.0)
    for g in 1:ng
        a2l.forc_pbot_not_downscaled_grc[g] = psrf
    end

    # Wind speed [m/s]
    wind = _read_var("WIND", 3.0)
    for g in 1:ng
        a2l.forc_wind_grc[g] = wind
        a2l.forc_u_grc[g] = wind
        a2l.forc_v_grc[g] = 0.0
    end

    # Longwave radiation [W/m2]
    flds = _read_var("FLDS", 250.0)
    for g in 1:ng
        a2l.forc_lwrad_not_downscaled_grc[g] = flds
    end

    # Shortwave radiation [W/m2] — split into VIS/NIR (50/50) then
    # direct/diffuse using CAM-derived polynomial (Fortran DATM CLMNCEP)
    fsds = _read_var("FSDS", 0.0)
    for g in 1:ng
        swndr = fsds * 0.50  # NIR half
        ratio_nir = clamp(0.29548 + 0.00504*swndr - 1.4957e-5*swndr^2 + 1.4881e-8*swndr^3, 0.01, 0.99)
        swvdr = fsds * 0.50  # VIS half
        ratio_vis = clamp(0.17639 + 0.00380*swvdr - 9.0039e-6*swvdr^2 + 8.1351e-9*swvdr^3, 0.01, 0.99)
        a2l.forc_solad_not_downscaled_grc[g, 1] = ratio_vis * swvdr        # VIS direct
        a2l.forc_solad_not_downscaled_grc[g, 2] = ratio_nir * swndr        # NIR direct
        a2l.forc_solai_grc[g, 1] = (1.0 - ratio_vis) * swvdr               # VIS diffuse
        a2l.forc_solai_grc[g, 2] = (1.0 - ratio_nir) * swndr               # NIR diffuse
        a2l.forc_solar_not_downscaled_grc[g] = fsds
    end

    # Precipitation [mm/s] — total, partition by temperature
    precip = _read_var("PRECTmms", 0.0)
    if precip < 0.0
        precip = 0.0
    end
    # Partition precipitation using a linear ramp matching Fortran DATM
    # (shr_precip_mod.F90): frac_rain = (T - TFRZ) * 0.5, clamped to [0,1]
    # All snow at T <= 0°C, all rain at T >= +2°C
    for g in 1:ng
        frac_rain = clamp((tbot - TFRZ) * 0.5, 0.0, 1.0)
        a2l.forc_rain_not_downscaled_grc[g] = precip * frac_rain
        a2l.forc_snow_not_downscaled_grc[g] = precip * (1.0 - frac_rain)
    end

    # Specific humidity → vapor pressure
    # Try QBOT first, then RH
    if haskey(ds, "QBOT")
        qbot = _read_var("QBOT", 0.003)
        for g in 1:ng
            # e = q * p / (0.622 + 0.378 * q)
            a2l.forc_vp_grc[g] = qbot * psrf / (0.622 + 0.378 * qbot)
        end
    elseif haskey(ds, "RH")
        rh = _read_var("RH", 70.0) / 100.0  # convert % to fraction
        for g in 1:ng
            # Saturation vapor pressure (Tetens formula)
            tc = tbot - TFRZ
            esat = 611.0 * exp(17.27 * tc / (tc + 237.3))
            a2l.forc_vp_grc[g] = rh * esat
        end
    end

    # Potential temperature
    for g in 1:ng
        pbot = a2l.forc_pbot_not_downscaled_grc[g]
        t = a2l.forc_t_not_downscaled_grc[g]
        a2l.forc_th_not_downscaled_grc[g] = t * (100000.0 / pbot)^(RAIR / CPAIR)
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

    # CO2 / O2 defaults
    for g in 1:ng
        if a2l.forc_pco2_grc[g] <= 0.0
            a2l.forc_pco2_grc[g] = 40.0  # ~400 ppm at sea level
        end
        if a2l.forc_po2_grc[g] <= 0.0
            a2l.forc_po2_grc[g] = 20900.0  # ~20.9% of 100 kPa
        end
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
