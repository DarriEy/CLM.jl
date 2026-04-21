# ==========================================================================
# Ported from: src/main/atm2lndType.F90
# Atmosphere-to-land forcing data types and initialization
# ==========================================================================

# --------------------------------------------------------------------------
# Parameters type
# --------------------------------------------------------------------------

"""
    Atm2LndParamsData

Parameters controlling atmosphere-to-land forcing adjustments:
rain/snow repartitioning, longwave downscaling, and lapse rates.

Ported from `atm2lnd_params_type` in `atm2lndType.F90`.
"""
Base.@kwdef mutable struct Atm2LndParamsData{FT<:Real}
    # true => repartition rain/snow from atm based on temperature
    repartition_rain_snow::Bool = false

    # true => downscale longwave radiation
    glcmec_downscale_longwave::Bool = false

    # Surface temperature lapse rate (K m-1)
    lapse_rate::FT = NaN

    # longwave radiation lapse rate (W m-2 m-1)
    lapse_rate_longwave::FT = NaN

    # Relative limit for how much longwave downscaling can be done (unitless)
    longwave_downscaling_limit::FT = NaN

    # Rain-snow ramp for glacier landunits
    # frac_rain = (temp - all_snow_t) * frac_rain_slope  (all_snow_t in K)
    precip_repartition_glc_all_snow_t::FT = NaN
    precip_repartition_glc_frac_rain_slope::FT = NaN

    # Rain-snow ramp for non-glacier landunits
    precip_repartition_nonglc_all_snow_t::FT = NaN
    precip_repartition_nonglc_frac_rain_slope::FT = NaN
end

# --------------------------------------------------------------------------
# Main atm2lnd data type
# --------------------------------------------------------------------------

"""
    Atm2LndData

Atmosphere-to-land forcing data. Contains gridcell-level atmospheric forcing
fields (not downscaled), column-level downscaled fields, and patch-level
time-averaged quantities.

Ported from `atm2lnd_type` in `atm2lndType.F90`.
"""
Base.@kwdef mutable struct Atm2LndData{FT<:Real}
    params::Atm2LndParamsData = Atm2LndParamsData()

    # --- atm->lnd not downscaled (gridcell-level) ---
    forc_u_grc                    ::Vector{FT} = Float64[]   # atm wind speed, east direction (m/s)
    forc_v_grc                    ::Vector{FT} = Float64[]   # atm wind speed, north direction (m/s)
    forc_wind_grc                 ::Vector{FT} = Float64[]   # atmospheric wind speed
    forc_hgt_grc                  ::Vector{FT} = Float64[]   # atmospheric reference height (m)
    forc_topo_grc                 ::Vector{FT} = Float64[]   # atmospheric surface height (m)
    forc_hgt_u_grc                ::Vector{FT} = Float64[]   # obs height of wind [m]
    forc_hgt_t_grc                ::Vector{FT} = Float64[]   # obs height of temperature [m]
    forc_hgt_q_grc                ::Vector{FT} = Float64[]   # obs height of humidity [m]
    forc_vp_grc                   ::Vector{FT} = Float64[]   # atmospheric vapor pressure (Pa)
    forc_pco2_grc                 ::Vector{FT} = Float64[]   # CO2 partial pressure (Pa)
    forc_pco2_240_patch           ::Vector{FT} = Float64[]   # 10-day mean CO2 partial pressure (Pa)
    forc_solad_not_downscaled_grc ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # direct beam radiation (numrad)
    forc_solai_grc                ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # diffuse radiation (numrad)
    forc_solar_not_downscaled_grc ::Vector{FT} = Float64[]   # incident solar radiation
    forc_solar_downscaled_col     ::Vector{FT} = Float64[]   # incident solar radiation (downscaled)
    forc_ndep_grc                 ::Vector{FT} = Float64[]   # nitrogen deposition rate (gN/m2/s)
    forc_pc13o2_grc               ::Vector{FT} = Float64[]   # C13O2 partial pressure (Pa)
    forc_po2_grc                  ::Vector{FT} = Float64[]   # O2 partial pressure (Pa)
    forc_po2_240_patch            ::Vector{FT} = Float64[]   # 10-day mean O2 partial pressure (Pa)
    forc_aer_grc                  ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # aerosol deposition array
    forc_pch4_grc                 ::Vector{FT} = Float64[]   # CH4 partial pressure (Pa)
    forc_o3_grc                   ::Vector{FT} = Float64[]   # ozone partial pressure (mol/mol)

    forc_t_not_downscaled_grc     ::Vector{FT} = Float64[]   # not downscaled atm temperature (K)
    forc_th_not_downscaled_grc    ::Vector{FT} = Float64[]   # not downscaled atm potential temperature (K)
    forc_pbot_not_downscaled_grc  ::Vector{FT} = Float64[]   # not downscaled atm pressure (Pa)
    forc_pbot240_downscaled_patch ::Vector{FT} = Float64[]   # 10-day mean downscaled atm pressure (Pa)
    forc_rho_not_downscaled_grc   ::Vector{FT} = Float64[]   # not downscaled atm density (kg/m**3)
    forc_lwrad_not_downscaled_grc ::Vector{FT} = Float64[]   # not downscaled atm downwrd IR longwave radiation (W/m**2)

    # --- atm->lnd downscaled (column-level) ---
    forc_t_downscaled_col         ::Vector{FT} = Float64[]   # downscaled atm temperature (K)
    forc_th_downscaled_col        ::Vector{FT} = Float64[]   # downscaled atm potential temperature (K)
    forc_pbot_downscaled_col      ::Vector{FT} = Float64[]   # downscaled atm pressure (Pa)
    forc_rho_downscaled_col       ::Vector{FT} = Float64[]   # downscaled atm density (kg/m**3)
    forc_lwrad_downscaled_col     ::Vector{FT} = Float64[]   # downscaled atm downwrd IR longwave radiation (W/m**2)
    forc_solad_downscaled_col     ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # direct beam radiation downscaled (numrad)

    # --- atm->lnd precipitation/humidity (gridcell-level, not downscaled) ---
    forc_rain_not_downscaled_grc      ::Vector{FT} = Float64[]   # rain rate (mm H2O/s)
    forc_snow_not_downscaled_grc      ::Vector{FT} = Float64[]   # snow rate (mm H2O/s)
    forc_q_not_downscaled_grc         ::Vector{FT} = Float64[]   # specific humidity (kg/kg)

    # --- atm->lnd precipitation/humidity (column-level, downscaled) ---
    forc_rain_downscaled_col          ::Vector{FT} = Float64[]   # rain rate (mm H2O/s)
    forc_snow_downscaled_col          ::Vector{FT} = Float64[]   # snow rate (mm H2O/s)
    forc_q_downscaled_col             ::Vector{FT} = Float64[]   # specific humidity (kg/kg)

    # --- time averaged quantities (patch-level) ---
    fsd24_patch                   ::Vector{FT} = Float64[]   # patch 24hr average of direct beam radiation
    fsd240_patch                  ::Vector{FT} = Float64[]   # patch 240hr average of direct beam radiation
    fsi24_patch                   ::Vector{FT} = Float64[]   # patch 24hr average of diffuse beam radiation
    fsi240_patch                  ::Vector{FT} = Float64[]   # patch 240hr average of diffuse beam radiation
    wind24_patch                  ::Vector{FT} = Float64[]   # patch 24-hour running mean of wind
    t_mo_patch                    ::Vector{FT} = Float64[]   # patch 30-day average temperature (K)
    t_mo_min_patch                ::Vector{FT} = Float64[]   # patch annual min of t_mo (K)
end

# --------------------------------------------------------------------------
# Helper: compute ramp parameters from snow/rain temperature endpoints
# --------------------------------------------------------------------------

"""
    compute_ramp_params(all_snow_t_c, all_rain_t_c) -> (all_snow_t_k, frac_rain_slope)

Convert rain/snow ramp temperature endpoints (in Celsius) into the all-snow
temperature in Kelvin and the fractional rain slope.
"""
function compute_ramp_params(all_snow_t_c::Real, all_rain_t_c::Real)
    frac_rain_slope = 1.0 / (all_rain_t_c - all_snow_t_c)
    all_snow_t_k = all_snow_t_c + TFRZ
    return (all_snow_t_k, frac_rain_slope)
end

# --------------------------------------------------------------------------
# Parameters constructor (mirrors atm2lnd_params_constructor)
# --------------------------------------------------------------------------

"""
    atm2lnd_params_init!(params; repartition_rain_snow, glcmec_downscale_longwave,
                          lapse_rate, [lapse_rate_longwave, longwave_downscaling_limit,
                          precip_repartition_glc_all_snow_t, precip_repartition_glc_all_rain_t,
                          precip_repartition_nonglc_all_snow_t, precip_repartition_nonglc_all_rain_t])

Initialize atm2lnd parameters. Validates inputs and converts rain/snow ramp
endpoints from Celsius to the internal representation.

Ported from `atm2lnd_params_constructor` in `atm2lndType.F90`.
"""
function atm2lnd_params_init!(params::Atm2LndParamsData;
        repartition_rain_snow::Bool,
        glcmec_downscale_longwave::Bool,
        lapse_rate::Real,
        lapse_rate_longwave::Real = NaN,
        longwave_downscaling_limit::Real = NaN,
        precip_repartition_glc_all_snow_t::Real = NaN,
        precip_repartition_glc_all_rain_t::Real = NaN,
        precip_repartition_nonglc_all_snow_t::Real = NaN,
        precip_repartition_nonglc_all_rain_t::Real = NaN)

    params.repartition_rain_snow = repartition_rain_snow
    params.glcmec_downscale_longwave = glcmec_downscale_longwave
    params.lapse_rate = lapse_rate

    if glcmec_downscale_longwave
        if isnan(lapse_rate_longwave)
            error("atm2lnd_params_init!: For glcmec_downscale_longwave true, lapse_rate_longwave must be provided")
        end
        if isnan(longwave_downscaling_limit)
            error("atm2lnd_params_init!: For glcmec_downscale_longwave true, longwave_downscaling_limit must be provided")
        end
        if longwave_downscaling_limit < 0.0 || longwave_downscaling_limit > 1.0
            error("atm2lnd_params_init!: longwave_downscaling_limit must be between 0 and 1")
        end
        params.lapse_rate_longwave = lapse_rate_longwave
        params.longwave_downscaling_limit = longwave_downscaling_limit
    else
        params.lapse_rate_longwave = NaN
        params.longwave_downscaling_limit = NaN
    end

    if repartition_rain_snow
        if isnan(precip_repartition_glc_all_snow_t)
            error("atm2lnd_params_init!: For repartition_rain_snow true, precip_repartition_glc_all_snow_t must be provided")
        end
        if isnan(precip_repartition_glc_all_rain_t)
            error("atm2lnd_params_init!: For repartition_rain_snow true, precip_repartition_glc_all_rain_t must be provided")
        end
        if isnan(precip_repartition_nonglc_all_snow_t)
            error("atm2lnd_params_init!: For repartition_rain_snow true, precip_repartition_nonglc_all_snow_t must be provided")
        end
        if isnan(precip_repartition_nonglc_all_rain_t)
            error("atm2lnd_params_init!: For repartition_rain_snow true, precip_repartition_nonglc_all_rain_t must be provided")
        end
        if precip_repartition_glc_all_rain_t <= precip_repartition_glc_all_snow_t
            error("atm2lnd_params_init!: Must have precip_repartition_glc_all_snow_t < precip_repartition_glc_all_rain_t")
        end
        if precip_repartition_nonglc_all_rain_t <= precip_repartition_nonglc_all_snow_t
            error("atm2lnd_params_init!: Must have precip_repartition_nonglc_all_snow_t < precip_repartition_nonglc_all_rain_t")
        end

        (glc_snow_k, glc_slope) = compute_ramp_params(
            precip_repartition_glc_all_snow_t, precip_repartition_glc_all_rain_t)
        params.precip_repartition_glc_all_snow_t = glc_snow_k
        params.precip_repartition_glc_frac_rain_slope = glc_slope

        (nonglc_snow_k, nonglc_slope) = compute_ramp_params(
            precip_repartition_nonglc_all_snow_t, precip_repartition_nonglc_all_rain_t)
        params.precip_repartition_nonglc_all_snow_t = nonglc_snow_k
        params.precip_repartition_nonglc_frac_rain_slope = nonglc_slope
    else
        params.precip_repartition_glc_all_snow_t = NaN
        params.precip_repartition_glc_frac_rain_slope = NaN
        params.precip_repartition_nonglc_all_snow_t = NaN
        params.precip_repartition_nonglc_frac_rain_slope = NaN
    end

    nothing
end

# --------------------------------------------------------------------------
# Allocation (mirrors InitAllocate)
# --------------------------------------------------------------------------

"""
    atm2lnd_init!(a2l, ng, nc, np)

Allocate and zero-initialize all arrays in `Atm2LndData`.
- `ng`: number of gridcells
- `nc`: number of columns
- `np`: number of patches

Ported from `InitAllocate` in `atm2lndType.F90`.
"""
function atm2lnd_init!(a2l::Atm2LndData{FT}, ng::Int, nc::Int, np::Int) where {FT}
    ival = 0.0

    # atm->lnd (gridcell-level)
    a2l.forc_u_grc                    = fill(ival, ng)
    a2l.forc_v_grc                    = fill(ival, ng)
    a2l.forc_wind_grc                 = fill(ival, ng)
    a2l.forc_hgt_grc                  = fill(ival, ng)
    a2l.forc_topo_grc                 = fill(ival, ng)
    a2l.forc_hgt_u_grc                = fill(ival, ng)
    a2l.forc_hgt_t_grc                = fill(ival, ng)
    a2l.forc_hgt_q_grc                = fill(ival, ng)
    a2l.forc_vp_grc                   = fill(ival, ng)
    a2l.forc_pco2_grc                 = fill(ival, ng)
    a2l.forc_solad_not_downscaled_grc = fill(ival, ng, NUMRAD)
    a2l.forc_solai_grc                = fill(ival, ng, NUMRAD)
    a2l.forc_solar_not_downscaled_grc = fill(ival, ng)
    a2l.forc_ndep_grc                 = fill(ival, ng)
    a2l.forc_pc13o2_grc               = fill(ival, ng)
    a2l.forc_po2_grc                  = fill(ival, ng)
    a2l.forc_aer_grc                  = fill(ival, ng, 14)
    a2l.forc_pch4_grc                 = fill(ival, ng)
    a2l.forc_o3_grc                   = fill(ival, ng)

    if varctl.use_luna
        a2l.forc_pco2_240_patch           = fill(ival, np)
        a2l.forc_po2_240_patch            = fill(ival, np)
        a2l.forc_pbot240_downscaled_patch = fill(ival, np)
    end

    # atm->lnd not downscaled (gridcell-level)
    a2l.forc_t_not_downscaled_grc     = fill(ival, ng)
    a2l.forc_pbot_not_downscaled_grc  = fill(ival, ng)
    a2l.forc_th_not_downscaled_grc    = fill(ival, ng)
    a2l.forc_rho_not_downscaled_grc   = fill(ival, ng)
    a2l.forc_lwrad_not_downscaled_grc = fill(ival, ng)

    # atm->lnd precipitation/humidity (gridcell-level, not downscaled)
    a2l.forc_rain_not_downscaled_grc  = fill(ival, ng)
    a2l.forc_snow_not_downscaled_grc  = fill(ival, ng)
    a2l.forc_q_not_downscaled_grc     = fill(ival, ng)

    # atm->lnd downscaled (column-level)
    a2l.forc_t_downscaled_col         = fill(ival, nc)
    a2l.forc_pbot_downscaled_col      = fill(ival, nc)
    a2l.forc_th_downscaled_col        = fill(ival, nc)
    a2l.forc_rho_downscaled_col       = fill(ival, nc)
    a2l.forc_lwrad_downscaled_col     = fill(ival, nc)
    a2l.forc_solad_downscaled_col     = fill(ival, nc, NUMRAD)
    a2l.forc_solar_downscaled_col     = fill(ival, nc)

    # atm->lnd precipitation/humidity (column-level, downscaled)
    a2l.forc_rain_downscaled_col      = fill(ival, nc)
    a2l.forc_snow_downscaled_col      = fill(ival, nc)
    a2l.forc_q_downscaled_col         = fill(ival, nc)

    # time-averaged (patch-level)
    a2l.fsd24_patch                   = fill(FT(NaN), np)
    a2l.fsd240_patch                  = fill(FT(NaN), np)
    a2l.fsi24_patch                   = fill(FT(NaN), np)
    a2l.fsi240_patch                  = fill(FT(NaN), np)
    if varctl.use_fates
        a2l.wind24_patch              = fill(FT(NaN), np)
    end
    a2l.t_mo_patch                    = fill(FT(NaN), np)
    a2l.t_mo_min_patch                = fill(FT(SPVAL), np)

    nothing
end

# --------------------------------------------------------------------------
# InitForTesting (mirrors Fortran InitForTesting)
# --------------------------------------------------------------------------

"""
    atm2lnd_init_for_testing!(a2l, ng, nc, np; params=nothing)

Initialize for unit testing. Allocates arrays and optionally sets params.

Ported from `InitForTesting` in `atm2lndType.F90`.
"""
function atm2lnd_init_for_testing!(a2l::Atm2LndData, ng::Int, nc::Int, np::Int;
        params::Union{Atm2LndParamsData, Nothing} = nothing)
    atm2lnd_init!(a2l, ng, nc, np)
    if params !== nothing
        a2l.params = params
    else
        a2l.params = Atm2LndParamsData(
            repartition_rain_snow = false,
            glcmec_downscale_longwave = false,
            lapse_rate = 0.01)
    end
    nothing
end

# --------------------------------------------------------------------------
# ReadNamelist stub (Fortran I/O not applicable)
# --------------------------------------------------------------------------

"""
    atm2lnd_read_namelist!(a2l; kwargs...)

Set parameters from keyword arguments (replaces Fortran namelist read).
"""
function atm2lnd_read_namelist!(a2l::Atm2LndData;
        repartition_rain_snow::Bool = false,
        glcmec_downscale_longwave::Bool = false,
        lapse_rate::Real = NaN,
        lapse_rate_longwave::Real = NaN,
        longwave_downscaling_limit::Real = NaN,
        precip_repartition_glc_all_snow_t::Real = NaN,
        precip_repartition_glc_all_rain_t::Real = NaN,
        precip_repartition_nonglc_all_snow_t::Real = NaN,
        precip_repartition_nonglc_all_rain_t::Real = NaN)

    atm2lnd_params_init!(a2l.params;
        repartition_rain_snow = repartition_rain_snow,
        glcmec_downscale_longwave = glcmec_downscale_longwave,
        lapse_rate = lapse_rate,
        lapse_rate_longwave = lapse_rate_longwave,
        longwave_downscaling_limit = longwave_downscaling_limit,
        precip_repartition_glc_all_snow_t = precip_repartition_glc_all_snow_t,
        precip_repartition_glc_all_rain_t = precip_repartition_glc_all_rain_t,
        precip_repartition_nonglc_all_snow_t = precip_repartition_nonglc_all_snow_t,
        precip_repartition_nonglc_all_rain_t = precip_repartition_nonglc_all_rain_t)
    nothing
end

# --------------------------------------------------------------------------
# InitHistory stub (Fortran history output not applicable in Julia)
# --------------------------------------------------------------------------

"""
    atm2lnd_init_history!(a2l, ng, nc, np)

Stub for history field registration. No-op in Julia port.

Ported from `InitHistory` in `atm2lndType.F90`.
"""
function atm2lnd_init_history!(a2l::Atm2LndData{FT}, ng::Int, nc::Int, np::Int) where {FT}
    nothing
end

# --------------------------------------------------------------------------
# InitAccBuffer stub
# --------------------------------------------------------------------------

"""
    atm2lnd_init_acc_buffer!(a2l)

Stub for accumulation buffer initialization. No-op in Julia port.

Ported from `InitAccBuffer` in `atm2lndType.F90`.
"""
function atm2lnd_init_acc_buffer!(a2l::Atm2LndData{FT}) where {FT}
    nothing
end

# --------------------------------------------------------------------------
# InitAccVars stub
# --------------------------------------------------------------------------

"""
    atm2lnd_init_acc_vars!(a2l, bounds)

Stub for accumulation variable initialization. No-op in Julia port.

Ported from `InitAccVars` in `atm2lndType.F90`.
"""
function atm2lnd_init_acc_vars!(a2l::Atm2LndData{FT}, bounds::UnitRange{Int}) where {FT}
    nothing
end

# --------------------------------------------------------------------------
# UpdateAccVars (mirrors Fortran UpdateAccVars)
# --------------------------------------------------------------------------

"""
    atm2lnd_update_acc_vars!(a2l, bounds_p, bounds_c;
        patch_gridcell, patch_column,
        fsd24_update!, fsd240_update!, fsi24_update!, fsi240_update!,
        [wind24_update!, tda_update!, pco2_240_update!, po2_240_update!, pbot240_update!])

Update time-accumulated forcing variables. The actual accumulation/extraction
is delegated to caller-provided update functions (replacing Fortran accumulMod).

Ported from `UpdateAccVars` in `atm2lndType.F90`.
"""
function atm2lnd_update_acc_vars!(a2l::Atm2LndData,
        bounds_p::UnitRange{Int},
        patch_gridcell::AbstractVector{Int},
        patch_column::AbstractVector{Int})

    # Accumulate direct beam radiation: FSD24, FSD240
    for p in bounds_p
        g = patch_gridcell[p]
        val = a2l.forc_solad_not_downscaled_grc[g, 1]
        a2l.fsd24_patch[p] = val
        a2l.fsd240_patch[p] = val
    end

    # Accumulate diffuse radiation: FSI24, FSI240
    for p in bounds_p
        g = patch_gridcell[p]
        val = a2l.forc_solai_grc[g, 1]
        a2l.fsi24_patch[p] = val
        a2l.fsi240_patch[p] = val
    end

    # CNDV: accumulate temperature (TDA)
    if varctl.use_cndv
        for p in bounds_p
            c = patch_column[p]
            t_val = a2l.forc_t_downscaled_col[c]
            a2l.t_mo_patch[p] = t_val
            a2l.t_mo_min_patch[p] = min(a2l.t_mo_min_patch[p], t_val)
        end
    end

    # FATES: accumulate wind
    if varctl.use_fates
        for p in bounds_p
            g = patch_gridcell[p]
            a2l.wind24_patch[p] = a2l.forc_wind_grc[g]
        end
    end

    # LUNA: accumulate pco2, po2, pbot
    if varctl.use_luna
        for p in bounds_p
            g = patch_gridcell[p]
            a2l.forc_pco2_240_patch[p] = a2l.forc_pco2_grc[g]
        end
        for p in bounds_p
            g = patch_gridcell[p]
            a2l.forc_po2_240_patch[p] = a2l.forc_po2_grc[g]
        end
        for p in bounds_p
            c = patch_column[p]
            a2l.forc_pbot240_downscaled_patch[p] = a2l.forc_pbot_downscaled_col[c]
        end
    end

    nothing
end

# --------------------------------------------------------------------------
# Restart stub
# --------------------------------------------------------------------------

"""
    atm2lnd_restart!(a2l, bounds)

Stub for restart variable read/write. No-op in Julia port.

Ported from `Restart` in `atm2lndType.F90`.
"""
function atm2lnd_restart!(a2l::Atm2LndData{FT}, bounds::UnitRange{Int}) where {FT}
    nothing
end

# --------------------------------------------------------------------------
# Clean (mirrors Fortran Clean)
# --------------------------------------------------------------------------

"""
    atm2lnd_clean!(a2l)

Reset all array fields to empty. Mirrors Fortran `Clean` subroutine.

Ported from `Clean` in `atm2lndType.F90`.
"""
function atm2lnd_clean!(a2l::Atm2LndData{FT}) where {FT}
    # atm->lnd gridcell
    a2l.forc_u_grc                    = FT[]
    a2l.forc_v_grc                    = FT[]
    a2l.forc_wind_grc                 = FT[]
    a2l.forc_hgt_grc                  = FT[]
    a2l.forc_topo_grc                 = FT[]
    a2l.forc_hgt_u_grc                = FT[]
    a2l.forc_hgt_t_grc                = FT[]
    a2l.forc_hgt_q_grc                = FT[]
    a2l.forc_vp_grc                   = FT[]
    a2l.forc_pco2_grc                 = FT[]
    a2l.forc_solad_not_downscaled_grc = Matrix{FT}(undef, 0, 0)
    a2l.forc_solai_grc                = Matrix{FT}(undef, 0, 0)
    a2l.forc_solar_not_downscaled_grc = FT[]
    a2l.forc_ndep_grc                 = FT[]
    a2l.forc_pc13o2_grc               = FT[]
    a2l.forc_po2_grc                  = FT[]
    a2l.forc_aer_grc                  = Matrix{FT}(undef, 0, 0)
    a2l.forc_pch4_grc                 = FT[]
    a2l.forc_o3_grc                   = FT[]

    # atm->lnd not downscaled
    a2l.forc_t_not_downscaled_grc     = FT[]
    a2l.forc_pbot_not_downscaled_grc  = FT[]
    a2l.forc_th_not_downscaled_grc    = FT[]
    a2l.forc_rho_not_downscaled_grc   = FT[]
    a2l.forc_lwrad_not_downscaled_grc = FT[]

    # atm->lnd precipitation/humidity (gridcell-level, not downscaled)
    a2l.forc_rain_not_downscaled_grc  = FT[]
    a2l.forc_snow_not_downscaled_grc  = FT[]
    a2l.forc_q_not_downscaled_grc     = FT[]

    # atm->lnd downscaled
    a2l.forc_t_downscaled_col         = FT[]
    a2l.forc_pbot_downscaled_col      = FT[]
    a2l.forc_th_downscaled_col        = FT[]
    a2l.forc_rho_downscaled_col       = FT[]
    a2l.forc_lwrad_downscaled_col     = FT[]
    a2l.forc_solad_downscaled_col     = Matrix{FT}(undef, 0, 0)
    a2l.forc_solar_downscaled_col     = FT[]

    # atm->lnd precipitation/humidity (column-level, downscaled)
    a2l.forc_rain_downscaled_col      = FT[]
    a2l.forc_snow_downscaled_col      = FT[]
    a2l.forc_q_downscaled_col         = FT[]

    # time-averaged
    a2l.fsd24_patch                   = FT[]
    a2l.fsd240_patch                  = FT[]
    a2l.fsi24_patch                   = FT[]
    a2l.fsi240_patch                  = FT[]
    if varctl.use_fates
        a2l.wind24_patch              = FT[]
    end
    a2l.t_mo_patch                    = FT[]
    a2l.t_mo_min_patch                = FT[]

    nothing
end
