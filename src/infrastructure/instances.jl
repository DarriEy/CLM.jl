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

    # --- CN vegetation ---
    bgc_vegetation::CNVegetationData = CNVegetationData()
    cn_products::CNProductsData      = CNProductsData()
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
    lakestate_init!(inst.lakestate, nc, np)
    solarabs_init!(inst.solarabs, np, nl)
    surfalb_init!(inst.surfalb, np, nc, ng)

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

    # --- CN vegetation ---
    cn_vegetation_init!(inst.bgc_vegetation, np, nc, ng;
                        nlevdecomp=nlevdecomp_full,
                        ndecomp_pools=ndecomp_pools,
                        ndecomp_cascade_transitions=ndecomp_cascade_transitions)
    cn_products_init!(inst.cn_products, ng)

    return nothing
end
