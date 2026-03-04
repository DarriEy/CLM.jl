# ==========================================================================
# Ported from: src/biogeophys/EnergyFluxType.F90
# Energy flux data type allocation and initialization
# ==========================================================================

"""
    EnergyFluxData

Energy flux data structure. Holds all energy flux variables at patch,
column, landunit, and gridcell levels, including sensible/latent heat fluxes,
radiation, wind stress, transpiration, conductances, and balance checks.

Ported from `energyflux_type` in `EnergyFluxType.F90`.
"""
Base.@kwdef mutable struct EnergyFluxData
    # --- Sensible heat fluxes (patch-level) ---
    eflx_sh_stem_patch           ::Vector{Float64} = Float64[]   # patch sensible heat flux from stem (W/m**2) [+ to atm]
    eflx_sh_grnd_patch           ::Vector{Float64} = Float64[]   # patch sensible heat flux from ground (W/m**2) [+ to atm]
    eflx_sh_veg_patch            ::Vector{Float64} = Float64[]   # patch sensible heat flux from leaves (W/m**2) [+ to atm]
    eflx_sh_snow_patch           ::Vector{Float64} = Float64[]   # patch sensible heat flux from snow (W/m**2) [+ to atm]
    eflx_sh_soil_patch           ::Vector{Float64} = Float64[]   # patch sensible heat flux from soil (W/m**2) [+ to atm]
    eflx_sh_h2osfc_patch         ::Vector{Float64} = Float64[]   # patch sensible heat flux from surface water (W/m**2) [+ to atm]
    eflx_sh_tot_patch            ::Vector{Float64} = Float64[]   # patch total sensible heat flux (W/m**2) [+ to atm]
    eflx_sh_tot_u_patch          ::Vector{Float64} = Float64[]   # patch urban total sensible heat flux (W/m**2) [+ to atm]
    eflx_sh_tot_r_patch          ::Vector{Float64} = Float64[]   # patch rural total sensible heat flux (W/m**2) [+ to atm]

    # --- Sensible heat fluxes (column-level) ---
    eflx_sh_precip_conversion_col ::Vector{Float64} = Float64[]  # col sensible heat flux from precipitation conversion (W/m**2) [+ to atm]

    # --- Latent heat fluxes (patch-level) ---
    eflx_lh_tot_patch            ::Vector{Float64} = Float64[]   # patch total latent heat flux (W/m**2) [+ to atm]
    eflx_lh_tot_u_patch          ::Vector{Float64} = Float64[]   # patch urban total latent heat flux (W/m**2) [+ to atm]
    eflx_lh_tot_r_patch          ::Vector{Float64} = Float64[]   # patch rural total latent heat flux (W/m**2) [+ to atm]
    eflx_lh_vegt_patch           ::Vector{Float64} = Float64[]   # patch transpiration heat flux from veg (W/m**2) [+ to atm]
    eflx_lh_vege_patch           ::Vector{Float64} = Float64[]   # patch evaporation heat flux from veg (W/m**2) [+ to atm]
    eflx_lh_grnd_patch           ::Vector{Float64} = Float64[]   # patch evaporation heat flux from ground (W/m**2) [+ to atm]

    # --- Soil/ground heat fluxes (patch-level) ---
    eflx_soil_grnd_patch         ::Vector{Float64} = Float64[]   # patch soil heat flux (W/m**2) [+ = into soil]
    eflx_soil_grnd_u_patch       ::Vector{Float64} = Float64[]   # patch urban soil heat flux (W/m**2) [+ = into soil]
    eflx_soil_grnd_r_patch       ::Vector{Float64} = Float64[]   # patch rural soil heat flux (W/m**2) [+ = into soil]

    # --- Longwave radiation (patch-level) ---
    eflx_lwrad_net_patch         ::Vector{Float64} = Float64[]   # patch net infrared (longwave) rad (W/m**2) [+ = to atm]
    eflx_lwrad_net_r_patch       ::Vector{Float64} = Float64[]   # patch rural net infrared (longwave) rad (W/m**2) [+ = to atm]
    eflx_lwrad_net_u_patch       ::Vector{Float64} = Float64[]   # patch urban net infrared (longwave) rad (W/m**2) [+ = to atm]
    eflx_lwrad_out_patch         ::Vector{Float64} = Float64[]   # patch emitted infrared (longwave) radiation (W/m**2)
    eflx_lwrad_out_r_patch       ::Vector{Float64} = Float64[]   # patch rural emitted infrared (longwave) rad (W/m**2)
    eflx_lwrad_out_u_patch       ::Vector{Float64} = Float64[]   # patch urban emitted infrared (longwave) rad (W/m**2)

    # --- Snow melt fluxes (column-level) ---
    eflx_snomelt_col             ::Vector{Float64} = Float64[]   # col snow melt heat flux (W/m**2)
    eflx_snomelt_r_col           ::Vector{Float64} = Float64[]   # col rural snow melt heat flux (W/m**2)
    eflx_snomelt_u_col           ::Vector{Float64} = Float64[]   # col urban snow melt heat flux (W/m**2)
    eflx_h2osfc_to_snow_col      ::Vector{Float64} = Float64[]   # col snow melt to h2osfc heat flux (W/m**2)

    # --- Ground/net heat fluxes (patch-level) ---
    eflx_gnet_patch              ::Vector{Float64} = Float64[]   # patch net heat flux into ground (W/m**2)
    eflx_grnd_lake_patch         ::Vector{Float64} = Float64[]   # patch net heat flux into lake/snow surface, excluding light transmission (W/m**2)

    # --- Dynamic balance flux (gridcell-level) ---
    eflx_dynbal_grc              ::Vector{Float64} = Float64[]   # grc dynamic land cover change conversion energy flux (W/m**2)

    # --- Bottom/layer heat fluxes (column-level) ---
    eflx_bot_col                 ::Vector{Float64} = Float64[]   # col heat flux from beneath the soil or ice column (W/m**2)
    eflx_fgr12_col               ::Vector{Float64} = Float64[]   # col ground heat flux between soil layers 1 and 2 (W/m**2)
    eflx_fgr_col                 ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col (rural) soil downward heat flux (W/m2) (ncols, nlevgrnd)
    eflx_building_heat_errsoi_col ::Vector{Float64} = Float64[]  # col heat flux to interior surface of walls and roof for errsoi check (W/m**2)

    # --- Urban energy fluxes (column-level) ---
    eflx_urban_ac_col            ::Vector{Float64} = Float64[]   # col urban air conditioning flux (W/m**2)
    eflx_urban_heat_col          ::Vector{Float64} = Float64[]   # col urban heating flux (W/m**2)

    # --- Urban energy fluxes (patch-level) ---
    eflx_anthro_patch            ::Vector{Float64} = Float64[]   # patch total anthropogenic heat flux (W/m**2)
    eflx_traffic_patch           ::Vector{Float64} = Float64[]   # patch traffic sensible heat flux (W/m**2)
    eflx_wasteheat_patch         ::Vector{Float64} = Float64[]   # patch sensible heat flux from domestic heating/cooling sources of waste heat (W/m**2)
    eflx_ventilation_patch       ::Vector{Float64} = Float64[]   # patch sensible heat flux from building ventilation (W/m**2)
    eflx_heat_from_ac_patch      ::Vector{Float64} = Float64[]   # patch sensible heat flux put back into canyon due to removal by AC (W/m**2)

    # --- Urban energy fluxes (landunit-level) ---
    eflx_traffic_lun             ::Vector{Float64} = Float64[]   # lun traffic sensible heat flux (W/m**2)
    eflx_wasteheat_lun           ::Vector{Float64} = Float64[]   # lun sensible heat flux from domestic heating/cooling sources of waste heat (W/m**2)
    eflx_ventilation_lun         ::Vector{Float64} = Float64[]   # lun sensible heat flux from building ventilation (W/m**2)
    eflx_heat_from_ac_lun        ::Vector{Float64} = Float64[]   # lun sensible heat flux to be put back into canyon due to removal by AC (W/m**2)
    eflx_building_lun            ::Vector{Float64} = Float64[]   # lun building heat flux from change in interior building air temperature (W/m**2)
    eflx_urban_ac_lun            ::Vector{Float64} = Float64[]   # lun urban air conditioning flux (W/m**2)
    eflx_urban_heat_lun          ::Vector{Float64} = Float64[]   # lun urban heating flux (W/m**2)

    # --- Derivatives of energy fluxes (patch-level) ---
    dgnetdT_patch                ::Vector{Float64} = Float64[]   # patch derivative of net ground heat flux wrt soil temp (W/m**2 K)
    netrad_patch                 ::Vector{Float64} = Float64[]   # patch net radiation (W/m**2) [+ = to sfc]
    cgrnd_patch                  ::Vector{Float64} = Float64[]   # patch deriv. of soil energy flux wrt to soil temp [W/m2/K]
    cgrndl_patch                 ::Vector{Float64} = Float64[]   # patch deriv. of soil latent heat flux wrt soil temp [W/m**2/K]
    cgrnds_patch                 ::Vector{Float64} = Float64[]   # patch deriv. of soil sensible heat flux wrt soil temp [W/m2/K]

    # --- Canopy radiation (patch-level) ---
    dlrad_patch                  ::Vector{Float64} = Float64[]   # patch downward longwave radiation below the canopy [W/m2]
    ulrad_patch                  ::Vector{Float64} = Float64[]   # patch upward longwave radiation above the canopy [W/m2]

    # --- Wind stress (patch-level) ---
    taux_patch                   ::Vector{Float64} = Float64[]   # patch wind (shear) stress: e-w (kg/m/s**2)
    tauy_patch                   ::Vector{Float64} = Float64[]   # patch wind (shear) stress: n-s (kg/m/s**2)

    # --- Conductance (patch-level) ---
    canopy_cond_patch            ::Vector{Float64} = Float64[]   # patch tracer conductance for canopy [m/s]

    # --- Transpiration (patch-level) ---
    btran_patch                  ::Vector{Float64} = Float64[]   # patch transpiration wetness factor (0 to 1)
    btran_min_patch              ::Vector{Float64} = Float64[]   # patch daily minimum transpiration wetness factor (0 to 1)
    btran_min_inst_patch         ::Vector{Float64} = Float64[]   # patch instantaneous daily minimum transpiration wetness factor (0 to 1)
    bsun_patch                   ::Vector{Float64} = Float64[]   # patch sunlit canopy transpiration wetness factor (0 to 1)
    bsha_patch                   ::Vector{Float64} = Float64[]   # patch shaded canopy transpiration wetness factor (0 to 1)

    # --- Roots (patch-level, 2D) ---
    rresis_patch                 ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # patch root resistance by layer (0-1) (npatches, nlevgrnd)

    # --- Latent heat (column-level) ---
    htvp_col                     ::Vector{Float64} = Float64[]   # latent heat of vapor of water (or sublimation) [J/kg]

    # --- Canopy heat (patch-level) ---
    dhsdt_canopy_patch           ::Vector{Float64} = Float64[]   # patch change in heat content of canopy (leaf+stem) (W/m**2) [+ to atm]

    # --- Balance checks (patch-level) ---
    errsoi_patch                 ::Vector{Float64} = Float64[]   # soil/lake energy conservation error (W/m**2)
    errseb_patch                 ::Vector{Float64} = Float64[]   # surface energy conservation error (W/m**2)
    errsol_patch                 ::Vector{Float64} = Float64[]   # solar radiation conservation error (W/m**2)
    errlon_patch                 ::Vector{Float64} = Float64[]   # longwave radiation conservation error (W/m**2)

    # --- Balance checks (column-level) ---
    errsoi_col                   ::Vector{Float64} = Float64[]   # soil/lake energy conservation error (W/m**2)
    errseb_col                   ::Vector{Float64} = Float64[]   # surface energy conservation error (W/m**2)
    errsol_col                   ::Vector{Float64} = Float64[]   # solar radiation conservation error (W/m**2)
    errlon_col                   ::Vector{Float64} = Float64[]   # longwave radiation conservation error (W/m**2)
end

"""
    energyflux_init!(ef::EnergyFluxData, np::Int, nc::Int, nl::Int, ng::Int)

Allocate and initialize all fields of an `EnergyFluxData` instance for
`np` patches, `nc` columns, `nl` landunits, and `ng` gridcells.
Real fields are initialized to `NaN`.

Ported from `energyflux_type%InitAllocate` in `EnergyFluxType.F90`.
"""
function energyflux_init!(ef::EnergyFluxData, np::Int, nc::Int, nl::Int, ng::Int)
    nlevgrnd = varpar.nlevgrnd

    # --- Sensible heat fluxes (patch-level) ---
    ef.eflx_sh_stem_patch            = fill(NaN, np)
    ef.eflx_sh_grnd_patch            = fill(NaN, np)
    ef.eflx_sh_veg_patch             = fill(NaN, np)
    ef.eflx_sh_snow_patch            = fill(NaN, np)
    ef.eflx_sh_soil_patch            = fill(NaN, np)
    ef.eflx_sh_h2osfc_patch          = fill(NaN, np)
    ef.eflx_sh_tot_patch             = fill(NaN, np)
    ef.eflx_sh_tot_u_patch           = fill(NaN, np)
    ef.eflx_sh_tot_r_patch           = fill(NaN, np)

    # --- Sensible heat fluxes (column-level) ---
    ef.eflx_sh_precip_conversion_col = fill(NaN, nc)

    # --- Latent heat fluxes (patch-level) ---
    ef.eflx_lh_tot_patch             = fill(NaN, np)
    ef.eflx_lh_tot_u_patch           = fill(NaN, np)
    ef.eflx_lh_tot_r_patch           = fill(NaN, np)
    ef.eflx_lh_vegt_patch            = fill(NaN, np)
    ef.eflx_lh_vege_patch            = fill(NaN, np)
    ef.eflx_lh_grnd_patch            = fill(NaN, np)

    # --- Soil/ground heat fluxes (patch-level) ---
    ef.eflx_soil_grnd_patch          = fill(NaN, np)
    ef.eflx_soil_grnd_u_patch        = fill(NaN, np)
    ef.eflx_soil_grnd_r_patch        = fill(NaN, np)

    # --- Longwave radiation (patch-level) ---
    ef.eflx_lwrad_net_patch          = fill(NaN, np)
    ef.eflx_lwrad_net_r_patch        = fill(NaN, np)
    ef.eflx_lwrad_net_u_patch        = fill(NaN, np)
    ef.eflx_lwrad_out_patch          = fill(NaN, np)
    ef.eflx_lwrad_out_r_patch        = fill(NaN, np)
    ef.eflx_lwrad_out_u_patch        = fill(NaN, np)

    # --- Snow melt fluxes (column-level) ---
    ef.eflx_snomelt_col              = fill(NaN, nc)
    ef.eflx_snomelt_r_col            = fill(NaN, nc)
    ef.eflx_snomelt_u_col            = fill(NaN, nc)
    ef.eflx_h2osfc_to_snow_col       = fill(NaN, nc)

    # --- Ground/net heat fluxes (patch-level) ---
    ef.eflx_gnet_patch               = fill(NaN, np)
    ef.eflx_grnd_lake_patch          = fill(NaN, np)

    # --- Dynamic balance flux (gridcell-level) ---
    ef.eflx_dynbal_grc               = fill(NaN, ng)

    # --- Bottom/layer heat fluxes (column-level) ---
    ef.eflx_bot_col                  = fill(NaN, nc)
    ef.eflx_fgr12_col                = fill(NaN, nc)
    ef.eflx_fgr_col                  = fill(NaN, nc, nlevgrnd)
    ef.eflx_building_heat_errsoi_col = fill(NaN, nc)

    # --- Urban energy fluxes (column-level) ---
    ef.eflx_urban_ac_col             = fill(NaN, nc)
    ef.eflx_urban_heat_col           = fill(NaN, nc)

    # --- Urban energy fluxes (patch-level) ---
    ef.eflx_anthro_patch             = fill(NaN, np)
    ef.eflx_traffic_patch            = fill(NaN, np)
    ef.eflx_wasteheat_patch          = fill(NaN, np)
    ef.eflx_ventilation_patch        = fill(NaN, np)
    ef.eflx_heat_from_ac_patch       = fill(NaN, np)

    # --- Urban energy fluxes (landunit-level) ---
    ef.eflx_traffic_lun              = fill(NaN, nl)
    ef.eflx_wasteheat_lun            = fill(NaN, nl)
    ef.eflx_ventilation_lun          = fill(NaN, nl)
    ef.eflx_heat_from_ac_lun         = fill(NaN, nl)
    ef.eflx_building_lun             = fill(NaN, nl)
    ef.eflx_urban_ac_lun             = fill(NaN, nl)
    ef.eflx_urban_heat_lun           = fill(NaN, nl)

    # --- Derivatives of energy fluxes (patch-level) ---
    ef.dgnetdT_patch                 = fill(NaN, np)
    ef.netrad_patch                  = fill(NaN, np)
    ef.cgrnd_patch                   = fill(NaN, np)
    ef.cgrndl_patch                  = fill(NaN, np)
    ef.cgrnds_patch                  = fill(NaN, np)

    # --- Canopy radiation (patch-level) ---
    ef.dlrad_patch                   = fill(NaN, np)
    ef.ulrad_patch                   = fill(NaN, np)

    # --- Wind stress (patch-level) ---
    ef.taux_patch                    = fill(NaN, np)
    ef.tauy_patch                    = fill(NaN, np)

    # --- Conductance (patch-level) ---
    ef.canopy_cond_patch             = fill(NaN, np)

    # --- Transpiration (patch-level) ---
    ef.btran_patch                   = fill(NaN, np)
    ef.btran_min_patch               = fill(NaN, np)
    ef.btran_min_inst_patch          = fill(NaN, np)
    ef.bsun_patch                    = fill(NaN, np)
    ef.bsha_patch                    = fill(NaN, np)

    # --- Roots (patch-level, 2D) ---
    ef.rresis_patch                  = fill(NaN, np, nlevgrnd)

    # --- Latent heat (column-level) ---
    ef.htvp_col                      = fill(NaN, nc)

    # --- Canopy heat (patch-level) ---
    ef.dhsdt_canopy_patch            = fill(NaN, np)

    # --- Balance checks (patch-level) ---
    ef.errsoi_patch                  = fill(NaN, np)
    ef.errseb_patch                  = fill(NaN, np)
    ef.errsol_patch                  = fill(NaN, np)
    ef.errlon_patch                  = fill(NaN, np)

    # --- Balance checks (column-level) ---
    ef.errsoi_col                    = fill(NaN, nc)
    ef.errseb_col                    = fill(NaN, nc)
    ef.errsol_col                    = fill(NaN, nc)
    ef.errlon_col                    = fill(NaN, nc)

    return nothing
end

"""
    energyflux_clean!(ef::EnergyFluxData)

Deallocate (reset to empty) all fields of an `EnergyFluxData` instance.

Ported from deallocation logic in `EnergyFluxType.F90`.
"""
function energyflux_clean!(ef::EnergyFluxData)
    # Patch-level 1D
    ef.eflx_sh_stem_patch            = Float64[]
    ef.eflx_sh_grnd_patch            = Float64[]
    ef.eflx_sh_veg_patch             = Float64[]
    ef.eflx_sh_snow_patch            = Float64[]
    ef.eflx_sh_soil_patch            = Float64[]
    ef.eflx_sh_h2osfc_patch          = Float64[]
    ef.eflx_sh_tot_patch             = Float64[]
    ef.eflx_sh_tot_u_patch           = Float64[]
    ef.eflx_sh_tot_r_patch           = Float64[]
    ef.eflx_lh_tot_patch             = Float64[]
    ef.eflx_lh_tot_u_patch           = Float64[]
    ef.eflx_lh_tot_r_patch           = Float64[]
    ef.eflx_lh_vegt_patch            = Float64[]
    ef.eflx_lh_vege_patch            = Float64[]
    ef.eflx_lh_grnd_patch            = Float64[]
    ef.eflx_soil_grnd_patch          = Float64[]
    ef.eflx_soil_grnd_u_patch        = Float64[]
    ef.eflx_soil_grnd_r_patch        = Float64[]
    ef.eflx_lwrad_net_patch          = Float64[]
    ef.eflx_lwrad_net_r_patch        = Float64[]
    ef.eflx_lwrad_net_u_patch        = Float64[]
    ef.eflx_lwrad_out_patch          = Float64[]
    ef.eflx_lwrad_out_r_patch        = Float64[]
    ef.eflx_lwrad_out_u_patch        = Float64[]
    ef.eflx_gnet_patch               = Float64[]
    ef.eflx_grnd_lake_patch          = Float64[]
    ef.eflx_anthro_patch             = Float64[]
    ef.eflx_traffic_patch            = Float64[]
    ef.eflx_wasteheat_patch          = Float64[]
    ef.eflx_ventilation_patch        = Float64[]
    ef.eflx_heat_from_ac_patch       = Float64[]
    ef.dgnetdT_patch                 = Float64[]
    ef.netrad_patch                  = Float64[]
    ef.cgrnd_patch                   = Float64[]
    ef.cgrndl_patch                  = Float64[]
    ef.cgrnds_patch                  = Float64[]
    ef.dlrad_patch                   = Float64[]
    ef.ulrad_patch                   = Float64[]
    ef.taux_patch                    = Float64[]
    ef.tauy_patch                    = Float64[]
    ef.canopy_cond_patch             = Float64[]
    ef.btran_patch                   = Float64[]
    ef.btran_min_patch               = Float64[]
    ef.btran_min_inst_patch          = Float64[]
    ef.bsun_patch                    = Float64[]
    ef.bsha_patch                    = Float64[]
    ef.dhsdt_canopy_patch            = Float64[]
    ef.errsoi_patch                  = Float64[]
    ef.errseb_patch                  = Float64[]
    ef.errsol_patch                  = Float64[]
    ef.errlon_patch                  = Float64[]

    # Patch-level 2D
    ef.rresis_patch                  = Matrix{Float64}(undef, 0, 0)

    # Column-level 1D
    ef.eflx_sh_precip_conversion_col = Float64[]
    ef.eflx_snomelt_col              = Float64[]
    ef.eflx_snomelt_r_col            = Float64[]
    ef.eflx_snomelt_u_col            = Float64[]
    ef.eflx_h2osfc_to_snow_col       = Float64[]
    ef.eflx_bot_col                  = Float64[]
    ef.eflx_fgr12_col                = Float64[]
    ef.eflx_building_heat_errsoi_col = Float64[]
    ef.eflx_urban_ac_col             = Float64[]
    ef.eflx_urban_heat_col           = Float64[]
    ef.htvp_col                      = Float64[]
    ef.errsoi_col                    = Float64[]
    ef.errseb_col                    = Float64[]
    ef.errsol_col                    = Float64[]
    ef.errlon_col                    = Float64[]

    # Column-level 2D
    ef.eflx_fgr_col                  = Matrix{Float64}(undef, 0, 0)

    # Landunit-level
    ef.eflx_traffic_lun              = Float64[]
    ef.eflx_wasteheat_lun            = Float64[]
    ef.eflx_ventilation_lun          = Float64[]
    ef.eflx_heat_from_ac_lun         = Float64[]
    ef.eflx_building_lun             = Float64[]
    ef.eflx_urban_ac_lun             = Float64[]
    ef.eflx_urban_heat_lun           = Float64[]

    # Gridcell-level
    ef.eflx_dynbal_grc               = Float64[]

    return nothing
end

"""
    energyflux_init_cold!(ef, col, lun, patch_data, t_grnd_col,
                          bounds_col, bounds_patch, bounds_lun;
                          is_simple_buildtemp=false, is_prog_buildtemp=false)

Initialize cold-start conditions for energy flux state variables.
Sets emitted longwave radiation from ground temperature, initializes urban
energy fluxes, and sets root resistance to zero.

Arguments:
- `ef`: EnergyFluxData instance
- `col`: ColumnData instance
- `lun`: LandunitData instance
- `patch_data`: PatchData instance
- `t_grnd_col`: ground temperature per column (K)
- `bounds_col`: UnitRange{Int} for columns
- `bounds_patch`: UnitRange{Int} for patches
- `bounds_lun`: UnitRange{Int} for landunits
- `is_simple_buildtemp`: whether using simple building temperature method
- `is_prog_buildtemp`: whether using prognostic building temperature method

Ported from `energyflux_type%InitCold` in `EnergyFluxType.F90`.
"""
function energyflux_init_cold!(ef::EnergyFluxData,
                                col::ColumnData,
                                lun::LandunitData,
                                patch_data::PatchData,
                                t_grnd_col::Vector{Float64},
                                bounds_col::UnitRange{Int},
                                bounds_patch::UnitRange{Int},
                                bounds_lun::UnitRange{Int};
                                is_simple_buildtemp::Bool = false,
                                is_prog_buildtemp::Bool = false)

    nlevgrnd = varpar.nlevgrnd

    # Columns — initialize urban building fluxes for simple building temp method
    if is_simple_buildtemp
        for c in bounds_col
            ef.eflx_building_heat_errsoi_col[c] = 0.0
            ef.eflx_urban_ac_col[c]             = 0.0
            ef.eflx_urban_heat_col[c]           = 0.0
        end
    end

    # Patches — set non-urban spval for urban-only fields, compute lwrad_out
    for p in bounds_patch
        c = patch_data.column[p]
        l = patch_data.landunit[p]

        if !lun.urbpoi[l]  # non-urban
            ef.eflx_lwrad_net_u_patch[p] = SPVAL
            ef.eflx_lwrad_out_u_patch[p] = SPVAL
            ef.eflx_lh_tot_u_patch[p]    = SPVAL
            ef.eflx_sh_tot_u_patch[p]    = SPVAL
            ef.eflx_soil_grnd_u_patch[p] = SPVAL
        end

        ef.eflx_lwrad_out_patch[p] = SB * (t_grnd_col[c])^4
    end

    # Patches — set urban/non-urban waste heat, traffic, ventilation, etc.
    for p in bounds_patch
        l = patch_data.landunit[p]

        if !lun.urbpoi[l]
            ef.eflx_traffic_lun[l]        = SPVAL
            ef.eflx_wasteheat_lun[l]      = SPVAL
            ef.eflx_ventilation_lun[l]    = SPVAL
            if is_prog_buildtemp
                ef.eflx_building_lun[l]   = 0.0
                ef.eflx_urban_ac_lun[l]   = 0.0
                ef.eflx_urban_heat_lun[l] = 0.0
            end

            ef.eflx_wasteheat_patch[p]    = 0.0
            ef.eflx_ventilation_patch[p]  = 0.0
            ef.eflx_heat_from_ac_patch[p] = 0.0
            ef.eflx_traffic_patch[p]      = 0.0
            if is_simple_buildtemp
                ef.eflx_anthro_patch[p]   = 0.0
            end
        else
            if is_prog_buildtemp
                ef.eflx_building_lun[l]   = 0.0
                ef.eflx_urban_ac_lun[l]   = 0.0
                ef.eflx_urban_heat_lun[l] = 0.0
                ef.eflx_ventilation_lun[l] = 0.0
            end
        end
    end

    # Initialize rresis, for use in ecosystemdyn
    for p in bounds_patch
        for lev in 1:nlevgrnd
            ef.rresis_patch[p, lev] = 0.0
        end
    end

    return nothing
end

# ==========================================================================
# The following subroutines depend on infrastructure modules that are not yet
# ported (history, restart/IO, accumulator). They are provided as stubs that
# document the Fortran interface and can be filled in when those modules
# become available.
# ==========================================================================

"""
    energyflux_init_history!(ef, bounds_col, bounds_patch, bounds_lun, bounds_grc;
                              is_simple_buildtemp=false, is_prog_buildtemp=false)

Register energy flux fields for history file output.

Ported from `energyflux_type%InitHistory` in `EnergyFluxType.F90`.
Requires history infrastructure (histFileMod) — stub until that module is ported.
"""
function energyflux_init_history!(ef::EnergyFluxData,
                                   bounds_col::UnitRange{Int},
                                   bounds_patch::UnitRange{Int},
                                   bounds_lun::UnitRange{Int},
                                   bounds_grc::UnitRange{Int};
                                   is_simple_buildtemp::Bool = false,
                                   is_prog_buildtemp::Bool = false)
    # Stub: history field registration will be added when histFileMod is ported.
    # All fields that would be registered:
    #   EFLX_DYNBAL, FSM, FSM_ICE, FSM_R, FSM_U, FIRA, FIRA_ICE, FIRA_R,
    #   FIRE, FIRE_ICE, FIRE_R, LWup, FCTR, FCEV, FGEV, FSH, FSH_ICE, FSH_R, Qh,
    #   Qle, EFLX_LH_TOT, EFLX_LH_TOT_ICE, EFLX_LH_TOT_R, Qstor, FSH_V,
    #   FSH_STEM, DHSDT_CANOPY, FSH_G, FGR, FGR_ICE, FGR_R, FIRA_U,
    #   EFLX_SOIL_GRND, FIRE_U, FSH_U, FSH_PRECIP_CONVERSION,
    #   EFLX_LH_TOT_U, FGR_U, Rnet, DLRAD, ULRAD, CGRND, CGRNDL, CGRNDS,
    #   EFLX_GNET, EFLX_GRND_LAKE, BUILDHEAT, URBAN_AC, URBAN_HEAT,
    #   EFLXBUILD, DGNETDT, FGR12, FGR_SOIL_R, TRAFFICFLUX, WASTEHEAT,
    #   VENTILATION, HEAT_FROM_AC, Qanth, TAUX, Qtau, TAUY,
    #   BTRAN, BTRANMN, RRESIS, ERRSOI, ERRSEB, ERRSOL
    return nothing
end

"""
    energyflux_restart!(ef, bounds_col, bounds_patch, bounds_lun;
                         flag="read", is_simple_buildtemp=false, is_prog_buildtemp=false)

Read/write energy flux state from/to restart file.

Ported from `energyflux_type%Restart` in `EnergyFluxType.F90`.
Requires NetCDF/restart infrastructure — stub until that module is ported.
"""
function energyflux_restart!(ef::EnergyFluxData,
                              bounds_col::UnitRange{Int},
                              bounds_patch::UnitRange{Int},
                              bounds_lun::UnitRange{Int};
                              flag::String = "read",
                              is_simple_buildtemp::Bool = false,
                              is_prog_buildtemp::Bool = false)
    # Stub: restart variable I/O will be added when restUtilMod/ncdio_pio is ported.
    # Variables that would be read/written:
    #   EFLX_LWRAD_OUT, URBAN_AC_L (landunit), URBAN_HEAT_L (landunit),
    #   EFLX_VENTILATION (landunit), URBAN_AC (column), URBAN_HEAT (column),
    #   BTRAN_MIN, BTRAN_MIN_INST, eflx_dynbal_dribbler
    return nothing
end

"""
    energyflux_init_acc_buffer!(ef, bounds_patch)

Initialize accumulation buffer for energy flux accumulated fields.

Ported from `energyflux_type%InitAccBuffer` in `EnergyFluxType.F90`.
Requires accumulation infrastructure (accumulMod) — stub until that module is ported.
"""
function energyflux_init_acc_buffer!(ef::EnergyFluxData,
                                     bounds_patch::UnitRange{Int})
    # Stub: accumulation field definitions will be added when accumulMod is ported.
    # Fields that would be initialized:
    #   BTRANAV (timeavg, 1 hour)
    return nothing
end

"""
    energyflux_init_acc_vars!(ef, bounds_patch; is_startup=false)

Initialize accumulated variables from the accumulation buffer.
Called for both initial and restart runs.

Ported from `energyflux_type%InitAccVars` in `EnergyFluxType.F90`.
Requires accumulation infrastructure (accumulMod) — stub until that module is ported.
"""
function energyflux_init_acc_vars!(ef::EnergyFluxData,
                                   bounds_patch::UnitRange{Int};
                                   is_startup::Bool = false)
    if is_startup
        for p in bounds_patch
            ef.btran_min_patch[p]      = SPVAL
            ef.btran_min_inst_patch[p] = SPVAL
        end
    end

    return nothing
end

"""
    energyflux_update_acc_vars!(ef, bounds_patch;
                                 end_cd=false, secs=0, dtime=0)

Update accumulated energy flux variables each timestep.
Handles hourly average of btran and daily min tracking.

Ported from `energyflux_type%UpdateAccVars` in `EnergyFluxType.F90`.
Core logic is ported; accumulator calls are stubs until accumulMod is ported.
"""
function energyflux_update_acc_vars!(ef::EnergyFluxData,
                                     bounds_patch::UnitRange{Int};
                                     end_cd::Bool = false,
                                     secs::Int = 0,
                                     dtime::Int = 0)
    # BTRANAV: track daily min btran
    for p in bounds_patch
        bt = ef.btran_patch[p]
        isfinite(bt) || continue

        # Track instantaneous minimum for the current day
        if hasfield(typeof(ef), :btran_min_inst_patch) && length(ef.btran_min_inst_patch) >= p
            if ef.btran_min_inst_patch[p] == SPVAL
                ef.btran_min_inst_patch[p] = bt
            else
                ef.btran_min_inst_patch[p] = min(bt, ef.btran_min_inst_patch[p])
            end

            # At end of day, copy to daily min and reset
            if end_cd
                if hasfield(typeof(ef), :btran_min_patch) && length(ef.btran_min_patch) >= p
                    ef.btran_min_patch[p] = ef.btran_min_inst_patch[p]
                end
                ef.btran_min_inst_patch[p] = SPVAL
            end
        end
    end

    return nothing
end
