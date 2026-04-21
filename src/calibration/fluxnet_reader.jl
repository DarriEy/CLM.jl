# ==========================================================================
# FLUXNET forcing reader for CLM.jl calibration
#
# Maps FLUXNET half-hourly CSV/NetCDF data to CLM atmospheric forcing fields.
# Standard FLUXNET variable names → CLM atm2lnd fields.
#
# Public API:
#   read_fluxnet_forcing!  — Read FLUXNET data and fill CLM forcing
#   FluxnetSiteInfo        — Site metadata
# ==========================================================================

"""
    FluxnetSiteInfo

Metadata for a FLUXNET site.

- `site_id`: FLUXNET site identifier (e.g., "US-Ha1")
- `lat`: latitude [degrees]
- `lon`: longitude [degrees]
- `elev`: elevation [m]
- `pft`: dominant PFT index
- `forcing_file`: path to FLUXNET forcing data (CSV or NetCDF)
- `obs_file`: path to observed fluxes (CSV or NetCDF)
"""
Base.@kwdef struct FluxnetSiteInfo
    site_id::String
    lat::Float64
    lon::Float64
    elev::Float64 = 0.0
    pft::Int = 1
    forcing_file::String = ""
    obs_file::String = ""
end

"""
    FluxnetForcing

Parsed FLUXNET forcing data for a single timestep.
All values in CLM-compatible units.
"""
Base.@kwdef struct FluxnetForcing
    ta::Float64 = 285.0        # Air temperature [K]
    pa::Float64 = 85000.0      # Air pressure [Pa]
    ws::Float64 = 3.0          # Wind speed [m/s]
    lw_in::Float64 = 250.0     # Incoming longwave [W/m2]
    sw_in::Float64 = 300.0     # Incoming shortwave [W/m2]
    precip::Float64 = 0.0      # Precipitation [mm/s → kg/m2/s]
    vpd::Float64 = 500.0       # Vapor pressure deficit [Pa]
    co2::Float64 = 400.0       # CO2 mixing ratio [ppm]
    rh::Float64 = 0.5          # Relative humidity [fraction]
end

"""
    FluxnetTarget

Observed FLUXNET flux data for calibration targets.
"""
Base.@kwdef struct FluxnetTarget
    gpp::Float64 = NaN         # Gross primary productivity [umol CO2/m2/s]
    et::Float64 = NaN          # Evapotranspiration [mm/day or W/m2]
    nee::Float64 = NaN         # Net ecosystem exchange [umol CO2/m2/s]
    sh::Float64 = NaN          # Sensible heat [W/m2]
    lh::Float64 = NaN          # Latent heat [W/m2]
    swc::Float64 = NaN         # Soil water content [m3/m3]
end

"""
    setup_fluxnet_forcing!(a2l, forcing::FluxnetForcing, ng::Int)

Fill CLM atmospheric forcing fields from parsed FLUXNET data.
Converts FLUXNET variables to CLM field naming and units.
"""
function setup_fluxnet_forcing!(a2l, forcing::FluxnetForcing, ng::Int)
    T = forcing.ta
    P = forcing.pa
    SW = forcing.sw_in
    LW = forcing.lw_in
    WS = forcing.ws

    # Compute vapor pressure from VPD and saturation
    es, _, _, _ = qsat(T, P)
    es_pa = es  # qsat returns Pa
    vp = max(es_pa - forcing.vpd, 10.0)  # ensure positive

    for g in 1:ng
        a2l.forc_t_not_downscaled_grc[g] = T
        a2l.forc_pbot_not_downscaled_grc[g] = P
        a2l.forc_th_not_downscaled_grc[g] = T * (100000.0 / P)^(RAIR / CPAIR)
        a2l.forc_rho_not_downscaled_grc[g] = P / (RAIR * T)
        a2l.forc_lwrad_not_downscaled_grc[g] = LW
        a2l.forc_vp_grc[g] = vp
        a2l.forc_hgt_grc[g] = 30.0
        a2l.forc_hgt_u_grc[g] = 30.0
        a2l.forc_hgt_t_grc[g] = 30.0
        a2l.forc_hgt_q_grc[g] = 30.0
        a2l.forc_topo_grc[g] = 0.0
        a2l.forc_wind_grc[g] = WS
        a2l.forc_u_grc[g] = WS
        a2l.forc_v_grc[g] = 0.0

        # Split shortwave into direct/diffuse (70/30 split)
        for b in 1:NUMRAD
            a2l.forc_solad_not_downscaled_grc[g, b] = SW * 0.35  # half to VIS, half to NIR
            a2l.forc_solai_grc[g, b] = SW * 0.15
        end
        a2l.forc_solar_not_downscaled_grc[g] = SW

        # Precipitation: split into rain/snow based on temperature
        precip_rate = forcing.precip  # kg/m2/s
        if T > TFRZ + 2.0
            a2l.forc_rain_not_downscaled_grc[g] = precip_rate
            a2l.forc_snow_not_downscaled_grc[g] = 0.0
        elseif T < TFRZ - 2.0
            a2l.forc_rain_not_downscaled_grc[g] = 0.0
            a2l.forc_snow_not_downscaled_grc[g] = precip_rate
        else
            frac_rain = (T - (TFRZ - 2.0)) / 4.0
            a2l.forc_rain_not_downscaled_grc[g] = precip_rate * frac_rain
            a2l.forc_snow_not_downscaled_grc[g] = precip_rate * (1.0 - frac_rain)
        end

        # CO2 and O2
        a2l.forc_pco2_grc[g] = forcing.co2 * 1e-6 * P  # ppm → partial pressure
        a2l.forc_po2_grc[g] = 0.209 * P  # ~21% of atmospheric pressure
    end
    return nothing
end
