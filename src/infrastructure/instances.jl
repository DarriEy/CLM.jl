# ============================================================================
# Instance factory — creates and initializes all CLM data instances
#
# Ported from: src/main/clm_instMod.F90
#
# The Fortran module declares global instances and calls Init() on each.
# Here we bundle them into a single CLMInstances container and provide
# clm_instInit!() to allocate/initialise everything in one call.
#
# Skipped from Fortran: file I/O, restart, history, namelist reading,
# FATES, ozone factory, dust factory, snow-cover-fraction factory,
# excess-ice streams, soil-water-retention-curve factory.
# ============================================================================

"""
    CLMInstances

Container holding all CLM data type instances.  Corresponds to the module-level
variables declared in `clm_instMod.F90`.
"""
Base.@kwdef mutable struct CLMInstances
    # --- Hierarchy / grid ---
    gridcell::GridcellData           = GridcellData()
    landunit::LandunitData           = LandunitData()
    column::ColumnData               = ColumnData()
    patch::PatchData                 = PatchData()

    # --- Topography ---
    topo::TopoData                   = TopoData()

    # --- Urban ---
    urbanparams::UrbanParamsData     = UrbanParamsData()

    # --- Physics ---
    temperature::TemperatureData     = TemperatureData()
    energyflux::EnergyFluxData       = EnergyFluxData()
    canopystate::CanopyStateData     = CanopyStateData()
    soilstate::SoilStateData         = SoilStateData()
    soilhydrology::SoilHydrologyData = SoilHydrologyData()
    water::WaterData                 = WaterData()
    frictionvel::FrictionVelocityData = FrictionVelocityData()
    lakestate::LakeStateData         = LakeStateData()
    solarabs::SolarAbsorbedData      = SolarAbsorbedData()
    surfalb::SurfaceAlbedoData       = SurfaceAlbedoData()
    surfalb_con::SurfaceAlbedoConstants = SurfaceAlbedoConstants()

    # --- Radiation diagnostics ---
    surfrad::SurfaceRadiationData    = SurfaceRadiationData()

    # --- Photosynthesis ---
    photosyns::PhotosynthesisData    = PhotosynthesisData()

    # --- Ozone ---
    ozone::OzoneData                 = OzoneData()

    # --- Aerosol deposition ---
    aerosol::AerosolData             = AerosolData()

    # --- Irrigation ---
    irrigation::IrrigationData       = IrrigationData()

    # --- Atmosphere / land coupling ---
    atm2lnd::Atm2LndData             = Atm2LndData()
    lnd2atm::Lnd2AtmData             = Lnd2AtmData()

    # --- Crop ---
    crop::CropData                   = CropData()

    # --- Soil biogeochemistry ---
    soilbiogeochem_state::SoilBiogeochemStateData           = SoilBiogeochemStateData()
    soilbiogeochem_carbonstate::SoilBiogeochemCarbonStateData = SoilBiogeochemCarbonStateData()
    soilbiogeochem_carbonflux::SoilBiogeochemCarbonFluxData   = SoilBiogeochemCarbonFluxData()
    soilbiogeochem_nitrogenstate::SoilBiogeochemNitrogenStateData = SoilBiogeochemNitrogenStateData()
    soilbiogeochem_nitrogenflux::SoilBiogeochemNitrogenFluxData   = SoilBiogeochemNitrogenFluxData()

    # --- Emissions ---
    dust_emis::DustEmisBaseData      = DustEmisBaseData()
    vocemis::VOCEmisData             = VOCEmisData()

    # --- Dry deposition velocity ---
    drydep::DryDepVelocityData = DryDepVelocityData()

    # --- Satellite phenology ---
    satellite_phenology::SatellitePhenologyData = SatellitePhenologyData()

    # --- Active layer ---
    active_layer::ActiveLayerData = ActiveLayerData()

    # --- Balance check ---
    balcheck::BalanceCheckData = BalanceCheckData()

    # --- CN vegetation ---
    bgc_vegetation::CNVegetationData = CNVegetationData()
    cn_products::CNProductsData      = CNProductsData()

    # --- Decomposition cascade ---
    decomp_cascade::DecompCascadeConData = DecompCascadeConData()

    # --- BGC decomposition state and params ---
    decomp_bgc_state::DecompBGCState = DecompBGCState()
    decomp_bgc_params::DecompBGCParams = DecompBGCParams()
    cn_shared_params::CNSharedParamsData = CNSharedParamsData()
    decomp_params::DecompParams = DecompParams()
    competition_state::SoilBGCCompetitionState = SoilBGCCompetitionState()
    competition_params::SoilBGCCompetitionParams = SoilBGCCompetitionParams()
    litter_params::LitterVertTranspParams = LitterVertTranspParams()

    # --- Saturated / infiltration excess runoff ---
    sat_excess_runoff::SaturatedExcessRunoffData = SaturatedExcessRunoffData()
    infilt_excess_runoff::InfiltrationExcessRunoffData = InfiltrationExcessRunoffData()

    # --- Snow cover fraction ---
    scf_method::SnowCoverFractionBase = SnowCoverFractionSwensonLawrence2012()

    # --- Calibration overrides (NaN = use default) ---
    overrides::CalibrationOverrides = CalibrationOverrides()

    # --- Surface input data (for monthly phenology re-reads) ---
    surfdata::Union{SurfaceInputData, Nothing} = nothing
end

"""
    clm_instInit!(inst::CLMInstances;
                  ng, nl, nc, np,
                  nlevdecomp_full, ndecomp_pools, ndecomp_cascade_transitions)

Allocate and cold-initialise every data instance inside `inst`.

Dimension arguments
- `ng` — number of gridcells
- `nl` — number of landunits
- `nc` — number of columns
- `np` — number of patches

Biogeochemistry dimensions (keyword, with defaults)
- `nlevdecomp_full`            — full number of decomposition levels (default 10)
- `ndecomp_pools`              — number of decomposition pools      (default 7)
- `ndecomp_cascade_transitions`— number of cascade transitions      (default 5)
"""
function clm_instInit!(inst::CLMInstances;
                       ng::Int,
                       nl::Int,
                       nc::Int,
                       np::Int,
                       nlevdecomp_full::Int = 10,
                       ndecomp_pools::Int = 7,
                       ndecomp_cascade_transitions::Int = 5)

    # --- Grid hierarchy ---
    gridcell_init!(inst.gridcell, ng)
    landunit_init!(inst.landunit, nl)
    column_init!(inst.column, nc)
    patch_init!(inst.patch, np)

    # --- Urban ---
    urbanparams_init!(inst.urbanparams, nl)

    # --- Topography (allocate only; cold-start requires valid col/lun linkage) ---
    topo_init_allocate!(inst.topo, nc)

    # --- Core physics ---
    temperature_init!(inst.temperature, np, nc, nl, ng)
    energyflux_init!(inst.energyflux, np, nc, nl, ng)
    canopystate_init!(inst.canopystate, np)
    soilstate_init!(inst.soilstate, np, nc)
    soilhydrology_init!(inst.soilhydrology, nc)
    water_init!(inst.water, nc, np, nl, ng)
    frictionvel_init!(inst.frictionvel, np, nc)
    frictionvel_read_nml!(inst.frictionvel)
    frictionvel_read_params!(inst.frictionvel)
    lakestate_init!(inst.lakestate, nc, np)
    solarabs_init!(inst.solarabs, np, nl)
    surfalb_init!(inst.surfalb, np, nc, ng)

    # --- Radiation diagnostics ---
    surfrad_init!(inst.surfrad, np)

    # --- Photosynthesis ---
    photosynthesis_data_init!(inst.photosyns, np)

    # --- Ozone ---
    ozone_init!(inst.ozone, np)

    # --- Aerosol deposition ---
    aerosol_init!(inst.aerosol, nc)

    # --- Irrigation ---
    irrigation_init_allocate!(inst.irrigation, np, nc, varpar.nlevsoi)

    # --- Atmosphere / land coupling ---
    atm2lnd_init!(inst.atm2lnd, ng, nc, np)
    lnd2atm_init!(inst.lnd2atm, ng, nc)

    # --- Crop ---
    crop_init!(inst.crop, np)

    # --- Soil biogeochemistry ---
    soil_bgc_state_init!(inst.soilbiogeochem_state, nc, np,
                         nlevdecomp_full, ndecomp_cascade_transitions)
    soil_bgc_carbon_state_init!(inst.soilbiogeochem_carbonstate, nc, ng,
                                nlevdecomp_full, ndecomp_pools)
    soil_bgc_carbon_flux_init!(inst.soilbiogeochem_carbonflux, nc,
                               nlevdecomp_full, ndecomp_pools,
                               ndecomp_cascade_transitions)
    soil_bgc_nitrogen_state_init!(inst.soilbiogeochem_nitrogenstate, nc, ng,
                                  nlevdecomp_full, ndecomp_pools)
    soil_bgc_nitrogen_flux_init!(inst.soilbiogeochem_nitrogenflux, nc,
                                 nlevdecomp_full, ndecomp_pools,
                                 ndecomp_cascade_transitions)

    # --- Emissions ---
    dust_emis_init!(inst.dust_emis, np)
    vocemis_init!(inst.vocemis, np, ng, 20, 20)

    # --- Dry deposition velocity ---
    drydep_init!(inst.drydep, np, 0)  # n_drydep=0 by default; set >0 to activate

    # --- Satellite phenology ---
    satellite_phenology_init!(inst.satellite_phenology, np)

    # --- Active layer ---
    active_layer_init!(inst.active_layer, nc)

    # --- Balance check ---
    balance_check_init!(inst.balcheck, 1800.0)

    # --- CN vegetation ---
    cn_vegetation_init!(inst.bgc_vegetation, np, nc, ng;
                        nlevdecomp=nlevdecomp_full,
                        ndecomp_pools=ndecomp_pools,
                        ndecomp_cascade_transitions=ndecomp_cascade_transitions)
    cn_products_init!(inst.cn_products, ng)

    # --- Saturated / infiltration excess runoff ---
    inst.sat_excess_runoff.fsat_col = fill(0.0, nc)
    inst.sat_excess_runoff.fcov_col = fill(0.0, nc)
    inst.infilt_excess_runoff.qinmax_col = fill(0.0, nc)

    # --- Decomposition cascade (allocate with defaults) ---
    decomp_cascade_con_init!(inst.decomp_cascade,
                              ndecomp_pools, ndecomp_cascade_transitions)

    return nothing
end

"""
    decomp_cascade_con_init!(cascade, ndecomp_pools, ndecomp_cascade_transitions)

Allocate decomposition cascade configuration arrays with default sizes.
The actual cascade parameters are set by `init_decomp_cascade_bgc!` or
`init_decomp_cascade_mimics!` during model initialization.
"""
function decomp_cascade_con_init!(cascade::DecompCascadeConData,
                                   ndecomp_pools::Int,
                                   ndecomp_cascade_transitions::Int)
    cascade.cascade_donor_pool             = zeros(Int, ndecomp_cascade_transitions)
    cascade.cascade_receiver_pool          = zeros(Int, ndecomp_cascade_transitions)
    cascade.floating_cn_ratio_decomp_pools = falses(ndecomp_pools)
    cascade.is_litter                      = falses(ndecomp_pools)
    cascade.is_soil                        = falses(ndecomp_pools)
    cascade.is_cwd                         = falses(ndecomp_pools)
    cascade.initial_cn_ratio               = zeros(ndecomp_pools)
    cascade.initial_stock                  = zeros(ndecomp_pools)
    cascade.is_metabolic                   = falses(ndecomp_pools)
    cascade.is_cellulose                   = falses(ndecomp_pools)
    cascade.is_lignin                      = falses(ndecomp_pools)
    cascade.spinup_factor                  = ones(ndecomp_pools)
    cascade.decomp_pool_name_restart       = fill("", ndecomp_pools)
    cascade.decomp_pool_name_history       = fill("", ndecomp_pools)
    cascade.decomp_pool_name_long          = fill("", ndecomp_pools)
    cascade.decomp_pool_name_short         = fill("", ndecomp_pools)
    cascade.cascade_step_name              = fill("", ndecomp_cascade_transitions)
    return nothing
end
