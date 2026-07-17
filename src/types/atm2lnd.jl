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
Base.@kwdef mutable struct Atm2LndData{FT<:Real,
                           V<:AbstractVector{FT},
                           M<:AbstractMatrix{FT}}
    params::Atm2LndParamsData = Atm2LndParamsData()

    # --- atm->lnd not downscaled (gridcell-level) ---
    forc_u_grc                    ::V = Float64[]   # atm wind speed, east direction (m/s)
    forc_v_grc                    ::V = Float64[]   # atm wind speed, north direction (m/s)
    forc_wind_grc                 ::V = Float64[]   # atmospheric wind speed
    forc_hgt_grc                  ::V = Float64[]   # atmospheric reference height (m)
    forc_topo_grc                 ::V = Float64[]   # atmospheric surface height (m)
    forc_hgt_u_grc                ::V = Float64[]   # obs height of wind [m]
    forc_hgt_t_grc                ::V = Float64[]   # obs height of temperature [m]
    forc_hgt_q_grc                ::V = Float64[]   # obs height of humidity [m]
    forc_vp_grc                   ::V = Float64[]   # atmospheric vapor pressure (Pa)
    forc_pco2_grc                 ::V = Float64[]   # CO2 partial pressure (Pa)
    forc_pco2_240_patch           ::V = Float64[]   # 10-day mean CO2 partial pressure (Pa)
    forc_solad_not_downscaled_grc ::M = Matrix{Float64}(undef, 0, 0)  # direct beam radiation (numrad)
    forc_solai_grc                ::M = Matrix{Float64}(undef, 0, 0)  # diffuse radiation (numrad)
    forc_solar_not_downscaled_grc ::V = Float64[]   # incident solar radiation
    forc_solar_downscaled_col     ::V = Float64[]   # incident solar radiation (downscaled)
    forc_ndep_grc                 ::V = Float64[]   # nitrogen deposition rate (gN/m2/s)
    forc_pc13o2_grc               ::V = Float64[]   # C13O2 partial pressure (Pa)
    forc_po2_grc                  ::V = Float64[]   # O2 partial pressure (Pa)
    forc_po2_240_patch            ::V = Float64[]   # 10-day mean O2 partial pressure (Pa)
    forc_aer_grc                  ::M = Matrix{Float64}(undef, 0, 0)  # aerosol deposition array
    forc_pch4_grc                 ::V = Float64[]   # CH4 partial pressure (Pa)
    forc_o3_grc                   ::V = Float64[]   # ozone partial pressure (mol/mol)

    forc_t_not_downscaled_grc     ::V = Float64[]   # not downscaled atm temperature (K)
    forc_th_not_downscaled_grc    ::V = Float64[]   # not downscaled atm potential temperature (K)
    forc_pbot_not_downscaled_grc  ::V = Float64[]   # not downscaled atm pressure (Pa)
    forc_pbot240_downscaled_patch ::V = Float64[]   # 10-day mean downscaled atm pressure (Pa)
    forc_rho_not_downscaled_grc   ::V = Float64[]   # not downscaled atm density (kg/m**3)
    forc_lwrad_not_downscaled_grc ::V = Float64[]   # not downscaled atm downwrd IR longwave radiation (W/m**2)

    # --- atm->lnd downscaled (column-level) ---
    forc_t_downscaled_col         ::V = Float64[]   # downscaled atm temperature (K)
    forc_th_downscaled_col        ::V = Float64[]   # downscaled atm potential temperature (K)
    forc_pbot_downscaled_col      ::V = Float64[]   # downscaled atm pressure (Pa)
    forc_rho_downscaled_col       ::V = Float64[]   # downscaled atm density (kg/m**3)
    forc_lwrad_downscaled_col     ::V = Float64[]   # downscaled atm downwrd IR longwave radiation (W/m**2)
    forc_solad_downscaled_col     ::M = Matrix{Float64}(undef, 0, 0)  # direct beam radiation downscaled (numrad)

    # --- atm->lnd precipitation/humidity (gridcell-level, not downscaled) ---
    forc_rain_not_downscaled_grc      ::V = Float64[]   # rain rate (mm H2O/s)
    forc_snow_not_downscaled_grc      ::V = Float64[]   # snow rate (mm H2O/s)
    forc_q_not_downscaled_grc         ::V = Float64[]   # specific humidity (kg/kg)

    # --- atm->lnd precipitation/humidity (column-level, downscaled) ---
    forc_rain_downscaled_col          ::V = Float64[]   # rain rate (mm H2O/s)
    forc_snow_downscaled_col          ::V = Float64[]   # snow rate (mm H2O/s)
    forc_q_downscaled_col             ::V = Float64[]   # specific humidity (kg/kg)

    # --- relative humidity forcing (gridcell-level) ---
    # Fortran carries forc_rh as an atm->lnd import field on wateratm2lndbulk_type.
    # This port has no coupler import, so it is DERIVED from the (already present)
    # q / T / pbot forcing by `atm2lnd_update_rh!` each step. It is the source term
    # for the RH30 (use_cn) and RH24 (use_fates) accumulators below.
    forc_rh_grc                   ::V = Float64[]   # relative humidity (%)

    # --- time averaged quantities (patch-level) ---
    fsd24_patch                   ::V = Float64[]   # patch 24hr average of direct beam radiation
    fsd240_patch                  ::V = Float64[]   # patch 240hr average of direct beam radiation
    fsi24_patch                   ::V = Float64[]   # patch 24hr average of diffuse beam radiation
    fsi240_patch                  ::V = Float64[]   # patch 240hr average of diffuse beam radiation
    wind24_patch                  ::V = Float64[]   # patch 24-hour running mean of wind
    t_mo_patch                    ::V = Float64[]   # patch 30-day average temperature (K)
    t_mo_min_patch                ::V = Float64[]   # patch annual min of t_mo (K)

    # TDA is a `timeavg` accumulator (accum_period = -30 days), NOT a runmean:
    # `extract_accum_field` returns the 30-day mean only at a period boundary and
    # SPVAL in between, so `t_mo_min` tracks the minimum of the 30-day MEANS, not
    # of the instantaneous temperature. These two scratch arrays are the in-array
    # equivalent of accumulMod's per-field `val`/`nsteps` buffers for TDA.
    tda_accum_patch               ::V = Float64[]   # running sum over the current 30-day window
    tda_naccum_patch              ::V = Float64[]   # steps accumulated in the current window

    # --- Fortran `wateratm2lndbulk_type` accumulators ---------------------------
    # Ported from `Wateratm2lndBulkType.F90` (InitAccBuffer / UpdateAccVars).
    # CLM.jl has no separate wateratm2lndbulk_type; these members are hosted on
    # Atm2LndData because that is already the home of the forcing accumulators
    # (fsd24/fsi24/wind24/forc_pco2_240) and because their source terms
    # (forc_rain_downscaled_col, forc_snow_downscaled_col, forc_rh_grc) live here.
    prec10_patch                  ::V = Float64[]   # 10-day running mean of total precip (mm H2O/s) [use_cn]
    prec30_patch                  ::V = Float64[]   # 30-day running mean of total precip (mm H2O/s) [use_cn]
    prec60_patch                  ::V = Float64[]   # 60-day running mean of total precip (mm H2O/s) [use_cn]
    rh30_patch                    ::V = Float64[]   # 30-day running mean of relative humidity (%)   [use_cn]
    prec365_col                   ::V = Float64[]   # 365-day running mean of total precip (mm H2O/s) [use_cndv]
    prec24_patch                  ::V = Float64[]   # 24-hr running mean of total precip (mm H2O/s)  [use_fates]
    rh24_patch                    ::V = Float64[]   # 24-hr running mean of relative humidity (%)    [use_fates]
end

Atm2LndData{FT}(; kwargs...) where {FT<:Real} =
    Atm2LndData{FT, Vector{FT}, Matrix{FT}}(; kwargs...)
Adapt.@adapt_structure Atm2LndData


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

    # relative humidity (derived from q/T/pbot each step by atm2lnd_update_rh!)
    a2l.forc_rh_grc                   = fill(ival, ng)

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
    a2l.tda_accum_patch               = fill(FT(0), np)
    a2l.tda_naccum_patch              = fill(FT(0), np)

    # wateratm2lndbulk accumulators — allocated under the same gates Fortran
    # registers them under (Wateratm2lndBulkType.F90::InitAccBuffer), so an
    # ungated run neither allocates nor accumulates them.
    if varctl.use_cn
        # Zero, NOT NaN: Fortran registers these with `init_value=0._r8`
        # (Wateratm2lndBulkType.F90::InitAccBuffer) and `InitAccVars` copies the
        # accumulator into the patch field BEFORE the first timestep, so they are
        # finite (0 on a cold start, the restart value on a warm start) by the time
        # any consumer runs. Allocating them NaN meant the step-1 consumers — the Li
        # fire fuel-moisture/ignition terms (prec10/prec60/rh30) and CNPhenology's
        # rain-triggered stress-deciduous onset — read NaN on the first step.
        # accum_runmean still returns `val` verbatim at nstep<=1, so a cold-start run
        # accumulates exactly as before.
        a2l.prec10_patch              = fill(FT(0), np)
        a2l.prec30_patch              = fill(FT(0), np)
        a2l.prec60_patch              = fill(FT(0), np)
        a2l.rh30_patch                = fill(FT(0), np)
    end
    if varctl.use_cndv
        a2l.prec365_col               = fill(FT(NaN), nc)
    end
    if varctl.use_fates
        a2l.prec24_patch              = fill(FT(NaN), np)
        a2l.rh24_patch                = fill(FT(NaN), np)
    end

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
    atm2lnd_update_rh!(a2l, bounds_grc)

Derive `forc_rh_grc` (relative humidity, %) from the specific-humidity,
temperature and pressure forcing: `rh = 100 * q / qsat(T, pbot)`, clamped to
[0, 100].

Fortran imports `forc_rh` from the coupler onto `wateratm2lndbulk_type`; this
port has no coupler, so RH is reconstructed from forcing that IS present. It is
the source term for the RH30 (`use_cn`) and RH24 (`use_fates`) accumulators.
"""
@kernel function _atm2lnd_rh_kernel!(forc_rh, @Const(forc_t), @Const(forc_pbot),
                                     @Const(forc_q), gmin::Int, gmax::Int)
    g = @index(Global)
    @inbounds if gmin <= g <= gmax
        t = forc_t[g]; p = forc_pbot[g]; q = forc_q[g]
        if isfinite(t) && isfinite(p) && isfinite(q) && p > zero(p)
            qs, _es = qsat_no_derivs(t, p)
            T = eltype(forc_rh)
            forc_rh[g] = qs > zero(qs) ? min(T(100), max(zero(T), T(100) * q / qs)) : zero(T)
        end
    end
end

function atm2lnd_update_rh!(a2l::Atm2LndData, bounds_grc::UnitRange{Int})
    isempty(bounds_grc) && return nothing
    length(a2l.forc_rh_grc) >= last(bounds_grc) || return nothing
    # Kernelized (was a host @inbounds loop that scalar-indexed device arrays and
    # broke the use_cn Metal composite). _launch! runs the KA CPU kernel on host
    # (byte-identical: each g writes its own forc_rh_grc[g]) and the device kernel on GPU.
    _launch!(_atm2lnd_rh_kernel!, a2l.forc_rh_grc,
        a2l.forc_t_not_downscaled_grc, a2l.forc_pbot_not_downscaled_grc,
        a2l.forc_q_not_downscaled_grc, first(bounds_grc), last(bounds_grc);
        ndrange = length(a2l.forc_rh_grc))
    return nothing
end

"""
    atm2lnd_update_acc_vars!(a2l, bounds_p, bounds_c, patch_gridcell, patch_column;
                             nstep=0, dtime=1800)

Update the time-accumulated forcing variables — every one a genuine `accumulMod`
running mean, with its window derived from Fortran's accumulation period in DAYS
and the model timestep (see [`accum_runmean`](@ref) / [`accum_window_steps`](@ref)):

| field                                | acctype | Fortran period | gate      |
|--------------------------------------|---------|----------------|-----------|
| `fsd24` / `fsi24`                    | runmean | -1 day         | always    |
| `fsd240` / `fsi240`                  | runmean | -10 days       | always    |
| `prec10` / `prec30` / `prec60`       | runmean | -10/-30/-60 d  | use_cn    |
| `rh30`                               | runmean | -30 days       | use_cn    |
| `prec365`                            | runmean | -365 days      | use_cndv  |
| `t_mo` (TDA) → `t_mo_min`            | timeavg | -30 days       | use_cndv  |
| `wind24` / `prec24` / `rh24`         | runmean | -1 day         | use_fates |
| `forc_pco2_240` / `_po2_240` / `_pbot240` | runmean | -10 days  | use_luna  |

BUG HISTORY: every one of these was previously a PASS-THROUGH — the kernels
assigned the INSTANTANEOUS forcing value to the "24hr"/"240hr"/"10-day mean"
field. The routine was called every step, so it looked wired, but no averaging
happened at all. Consequences: LUNA acclimated to the instantaneous
CO2/O2/pressure rather than their 10-day means (`LunaMod` `CO2_p240`/`O2_p240`/
`forc_pbot10`); MEGAN/VOC saw instantaneous rather than 24h/240h radiation; and
CNDV's `t_mo_min` tracked the coldest INSTANT of the year instead of the coldest
30-day MEAN, which is a completely different (and far colder) bioclimatic limit.

Ported from `UpdateAccVars` in `atm2lndType.F90` and
`Wateratm2lndBulkType.F90`.
"""
function atm2lnd_update_acc_vars!(a2l::Atm2LndData,
        bounds_p::UnitRange{Int},
        patch_gridcell::AbstractVector{Int},
        patch_column::AbstractVector{Int};
        bounds_c::UnitRange{Int} = 1:0,
        nstep::Int = 0,
        dtime::Int = 1800)

    isempty(bounds_p) && return nothing
    pmin = first(bounds_p)
    pmax = last(bounds_p)

    w1   = accum_window_steps(1,   dtime)
    w10  = accum_window_steps(10,  dtime)
    w30  = accum_window_steps(30,  dtime)
    w60  = accum_window_steps(60,  dtime)
    w365 = accum_window_steps(365, dtime)

    # FSD24 (-1 day) / FSD240 (-10 days) — direct-beam radiation.
    _launch!(_a2l_solad_kernel!, a2l.fsd24_patch, a2l.fsd240_patch,
        patch_gridcell, a2l.forc_solad_not_downscaled_grc, nstep, w1, w10, pmin, pmax;
        ndrange = length(a2l.fsd24_patch))

    # FSI24 (-1 day) / FSI240 (-10 days) — diffuse radiation.
    _launch!(_a2l_solai_kernel!, a2l.fsi24_patch, a2l.fsi240_patch,
        patch_gridcell, a2l.forc_solai_grc, nstep, w1, w10, pmin, pmax;
        ndrange = length(a2l.fsi24_patch))

    # PREC10 / PREC30 / PREC60 / RH30 (use_cn). Precip is the column-downscaled
    # rain + snow mapped onto the patch; RH is the gridcell forcing RH.
    if varctl.use_cn && length(a2l.prec10_patch) >= pmax
        _launch!(_a2l_preccn_kernel!, a2l.prec10_patch, a2l.prec30_patch,
            a2l.prec60_patch, a2l.rh30_patch, patch_column, patch_gridcell,
            a2l.forc_rain_downscaled_col, a2l.forc_snow_downscaled_col,
            a2l.forc_rh_grc, nstep, w10, w30, w60, pmin, pmax;
            ndrange = length(a2l.prec10_patch))
    end

    # CNDV: TDA (30-day TIMEAVG) → t_mo, and the running min of those 30-day
    # means → t_mo_min. PREC365 is COLUMN-level (Fortran: "we cannot use a
    # patch-level accumulator for CNDV because this is used for establishment,
    # so must be available for inactive patches").
    if varctl.use_cndv
        _launch!(_a2l_cndv_kernel!, a2l.t_mo_patch, a2l.t_mo_min_patch,
            a2l.tda_accum_patch, a2l.tda_naccum_patch,
            patch_column, a2l.forc_t_downscaled_col, nstep, w30, pmin, pmax;
            ndrange = length(a2l.t_mo_patch))

        if !isempty(bounds_c) && length(a2l.prec365_col) >= last(bounds_c)
            _launch!(_a2l_prec365_kernel!, a2l.prec365_col,
                a2l.forc_rain_downscaled_col, a2l.forc_snow_downscaled_col,
                nstep, w365, first(bounds_c), last(bounds_c);
                ndrange = length(a2l.prec365_col))
        end
    end

    # FATES: WIND24 / PREC24 / RH24 (all -1 day).
    if varctl.use_fates
        _launch!(_a2l_fates_kernel!, a2l.wind24_patch, a2l.prec24_patch,
            a2l.rh24_patch, patch_gridcell, patch_column, a2l.forc_wind_grc,
            a2l.forc_rain_downscaled_col, a2l.forc_snow_downscaled_col,
            a2l.forc_rh_grc, nstep, w1, pmin, pmax;
            ndrange = length(a2l.wind24_patch))
    end

    # LUNA: pco2_240 / po2_240 / pbot240 (all -10 days).
    if varctl.use_luna
        _launch!(_a2l_luna_kernel!, a2l.forc_pco2_240_patch,
            a2l.forc_po2_240_patch, a2l.forc_pbot240_downscaled_patch,
            patch_gridcell, patch_column, a2l.forc_pco2_grc, a2l.forc_po2_grc,
            a2l.forc_pbot_downscaled_col, nstep, w10, pmin, pmax;
            ndrange = length(a2l.forc_pco2_240_patch))
    end

    nothing
end

# Per-patch running means of direct-beam radiation: FSD24 (1 day), FSD240 (10 days).
@kernel function _a2l_solad_kernel!(fsd24_patch, fsd240_patch,
        @Const(patch_gridcell), @Const(forc_solad_not_downscaled_grc),
        nstep::Int, w1::Int, w10::Int, pmin::Int, pmax::Int)
    p = @index(Global)
    @inbounds if pmin <= p <= pmax
        g = patch_gridcell[p]
        val = forc_solad_not_downscaled_grc[g, 1]
        if isfinite(val)
            fsd24_patch[p]  = accum_runmean(fsd24_patch[p],  val, nstep, w1)
            fsd240_patch[p] = accum_runmean(fsd240_patch[p], val, nstep, w10)
        end
    end
end

# Per-patch running means of diffuse radiation: FSI24 (1 day), FSI240 (10 days).
@kernel function _a2l_solai_kernel!(fsi24_patch, fsi240_patch,
        @Const(patch_gridcell), @Const(forc_solai_grc),
        nstep::Int, w1::Int, w10::Int, pmin::Int, pmax::Int)
    p = @index(Global)
    @inbounds if pmin <= p <= pmax
        g = patch_gridcell[p]
        val = forc_solai_grc[g, 1]
        if isfinite(val)
            fsi24_patch[p]  = accum_runmean(fsi24_patch[p],  val, nstep, w1)
            fsi240_patch[p] = accum_runmean(fsi240_patch[p], val, nstep, w10)
        end
    end
end

# use_cn: PREC10/PREC30/PREC60 (patch running means of column rain+snow) and
# RH30 (patch running mean of gridcell RH). Feeds the Li fire schemes' fuel
# moisture / ignition terms and CNPhenology's rain-triggered stress-decid onset.
@kernel function _a2l_preccn_kernel!(prec10_patch, prec30_patch, prec60_patch,
        rh30_patch, @Const(patch_column), @Const(patch_gridcell),
        @Const(forc_rain_downscaled_col), @Const(forc_snow_downscaled_col),
        @Const(forc_rh_grc), nstep::Int, w10::Int, w30::Int, w60::Int,
        pmin::Int, pmax::Int)
    p = @index(Global)
    @inbounds if pmin <= p <= pmax
        c = patch_column[p]
        g = patch_gridcell[p]
        prec = forc_rain_downscaled_col[c] + forc_snow_downscaled_col[c]
        if isfinite(prec)
            prec10_patch[p] = accum_runmean(prec10_patch[p], prec, nstep, w10)
            prec30_patch[p] = accum_runmean(prec30_patch[p], prec, nstep, w30)
            prec60_patch[p] = accum_runmean(prec60_patch[p], prec, nstep, w60)
        end
        rh = forc_rh_grc[g]
        if isfinite(rh)
            rh30_patch[p] = accum_runmean(rh30_patch[p], rh, nstep, w30)
        end
    end
end

# use_cndv: TDA is a 30-day TIMEAVG (not a runmean). Accumulate the sum over the
# window; at a window boundary emit the mean into t_mo and fold it into the
# running minimum t_mo_min, then reset. This mirrors accumulMod's timeavg
# extract, which returns the mean only at `mod(nstep, period) == 0` and SPVAL
# otherwise — so `min(t_mo_min, SPVAL)` is a no-op on non-boundary steps and
# t_mo_min ends up tracking the coldest 30-day MEAN.
@kernel function _a2l_cndv_kernel!(t_mo_patch, t_mo_min_patch,
        tda_accum_patch, tda_naccum_patch,
        @Const(patch_column), @Const(forc_t_downscaled_col),
        nstep::Int, w30::Int, pmin::Int, pmax::Int)
    p = @index(Global)
    @inbounds if pmin <= p <= pmax
        T = eltype(t_mo_patch)
        c = patch_column[p]
        t_val = forc_t_downscaled_col[c]
        if isfinite(t_val)
            tda_accum_patch[p]  += t_val
            tda_naccum_patch[p] += one(T)
        end
        if nstep % w30 == 0
            n = tda_naccum_patch[p]
            if n > zero(T)
                avg = tda_accum_patch[p] / n
                t_mo_patch[p] = avg
                cur = t_mo_min_patch[p]
                t_mo_min_patch[p] = (cur == T(SPVAL) || !isfinite(cur)) ? avg : min(cur, avg)
            end
            tda_accum_patch[p]  = zero(T)
            tda_naccum_patch[p] = zero(T)
        end
    end
end

# use_cndv: PREC365 — COLUMN-level 365-day running mean of rain+snow.
@kernel function _a2l_prec365_kernel!(prec365_col,
        @Const(forc_rain_downscaled_col), @Const(forc_snow_downscaled_col),
        nstep::Int, w365::Int, cmin::Int, cmax::Int)
    c = @index(Global)
    @inbounds if cmin <= c <= cmax
        prec = forc_rain_downscaled_col[c] + forc_snow_downscaled_col[c]
        if isfinite(prec)
            prec365_col[c] = accum_runmean(prec365_col[c], prec, nstep, w365)
        end
    end
end

# use_fates: WIND24 / PREC24 / RH24 — all 1-day running means (SPITFIRE inputs).
@kernel function _a2l_fates_kernel!(wind24_patch, prec24_patch, rh24_patch,
        @Const(patch_gridcell), @Const(patch_column), @Const(forc_wind_grc),
        @Const(forc_rain_downscaled_col), @Const(forc_snow_downscaled_col),
        @Const(forc_rh_grc), nstep::Int, w1::Int, pmin::Int, pmax::Int)
    p = @index(Global)
    @inbounds if pmin <= p <= pmax
        g = patch_gridcell[p]
        c = patch_column[p]
        wnd = forc_wind_grc[g]
        if isfinite(wnd)
            wind24_patch[p] = accum_runmean(wind24_patch[p], wnd, nstep, w1)
        end
        prec = forc_rain_downscaled_col[c] + forc_snow_downscaled_col[c]
        if isfinite(prec)
            prec24_patch[p] = accum_runmean(prec24_patch[p], prec, nstep, w1)
        end
        rh = forc_rh_grc[g]
        if isfinite(rh)
            rh24_patch[p] = accum_runmean(rh24_patch[p], rh, nstep, w1)
        end
    end
end

# use_luna: 10-day running means of CO2/O2 partial pressure and air pressure.
# These are LunaMod's CO2_p240 / O2_p240 / forc_pbot10 — the acclimation climate.
@kernel function _a2l_luna_kernel!(forc_pco2_240_patch, forc_po2_240_patch,
        forc_pbot240_downscaled_patch, @Const(patch_gridcell), @Const(patch_column),
        @Const(forc_pco2_grc), @Const(forc_po2_grc), @Const(forc_pbot_downscaled_col),
        nstep::Int, w10::Int, pmin::Int, pmax::Int)
    p = @index(Global)
    @inbounds if pmin <= p <= pmax
        g = patch_gridcell[p]
        c = patch_column[p]
        vco2 = forc_pco2_grc[g]
        if isfinite(vco2)
            forc_pco2_240_patch[p] = accum_runmean(forc_pco2_240_patch[p], vco2, nstep, w10)
        end
        vo2 = forc_po2_grc[g]
        if isfinite(vo2)
            forc_po2_240_patch[p] = accum_runmean(forc_po2_240_patch[p], vo2, nstep, w10)
        end
        vpb = forc_pbot_downscaled_col[c]
        if isfinite(vpb)
            forc_pbot240_downscaled_patch[p] =
                accum_runmean(forc_pbot240_downscaled_patch[p], vpb, nstep, w10)
        end
    end
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

    a2l.forc_rh_grc                   = FT[]

    # time-averaged
    a2l.fsd24_patch                   = FT[]
    a2l.fsd240_patch                  = FT[]
    a2l.fsi24_patch                   = FT[]
    a2l.fsi240_patch                  = FT[]
    if varctl.use_fates
        a2l.wind24_patch              = FT[]
        a2l.prec24_patch              = FT[]
        a2l.rh24_patch                = FT[]
    end
    a2l.t_mo_patch                    = FT[]
    a2l.t_mo_min_patch                = FT[]
    a2l.tda_accum_patch               = FT[]
    a2l.tda_naccum_patch              = FT[]
    if varctl.use_cn
        a2l.prec10_patch              = FT[]
        a2l.prec30_patch              = FT[]
        a2l.prec60_patch              = FT[]
        a2l.rh30_patch                = FT[]
    end
    if varctl.use_cndv
        a2l.prec365_col               = FT[]
    end

    nothing
end
