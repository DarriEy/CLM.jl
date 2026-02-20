#!/bin/bash
# ==========================================================================
# Ralph Loop Module List for CLM → Julia Port
# Format: TIER:FORTRAN_PATH:JULIA_PATH:DESCRIPTION:MAX_LINES
# Modules are in strict dependency order within each tier.
# ==========================================================================

MODULES=(
    # --- TIER 2: Core Spatial Types ---
    "2:src/main/GridcellType.F90:src/types/gridcell.jl:Gridcell type definition:130"
    "2:src/main/LandunitType.F90:src/types/landunit.jl:Landunit type definition:162"
    "2:src/main/ColumnType.F90:src/types/column.jl:Column type definition:252"
    "2:src/main/PatchType.F90:src/types/patch.jl:Patch type definition:210"
    "2:src/main/decompMod.F90:src/infrastructure/decomp.jl:Domain decomposition:350"
    "2:src/main/filterMod.F90:src/infrastructure/filters.jl:Filter to mask conversion:400"

    # --- TIER 3: Biogeophysics State Types ---
    "3:src/biogeophys/TemperatureType.F90:src/types/temperature.jl:Temperature state:1819"
    "3:src/biogeophys/EnergyFluxType.F90:src/types/energy_flux.jl:Energy flux state:1053"
    "3:src/biogeophys/SoilStateType.F90:src/types/soil_state.jl:Soil state:700"
    "3:src/biogeophys/SoilHydrologyType.F90:src/types/soil_hydrology.jl:Soil hydrology state:600"
    "3:src/biogeophys/CanopyStateType.F90:src/types/canopy_state.jl:Canopy state:600"
    "3:src/biogeophys/LakeStateType.F90:src/types/lake_state.jl:Lake state:400"
    "3:src/biogeophys/SurfaceAlbedoType.F90:src/types/surface_albedo.jl:Albedo state:500"
    "3:src/biogeophys/SolarAbsorbedType.F90:src/types/solar_absorbed.jl:Solar absorption:300"
    "3:src/biogeophys/UrbanParamsType.F90:src/types/urban_params.jl:Urban parameters:400"
    "3:src/biogeophys/WaterInfoBaseType.F90:src/types/water_info_base.jl:Water info base:100"
    "3:src/biogeophys/WaterStateType.F90:src/types/water_state.jl:Water state:500"
    "3:src/biogeophys/WaterFluxType.F90:src/types/water_flux.jl:Water flux:500"
    "3:src/biogeophys/WaterDiagnosticBulkType.F90:src/types/water_diagnostic_bulk.jl:Water diagnostic bulk:1257"
    "3:src/biogeophys/WaterStateBulkType.F90:src/types/water_state_bulk.jl:Water state bulk:800"
    "3:src/biogeophys/WaterFluxBulkType.F90:src/types/water_flux_bulk.jl:Water flux bulk:800"
    "3:src/biogeophys/WaterBalanceType.F90:src/types/water_balance.jl:Water balance:400"
    "3:src/biogeophys/WaterType.F90:src/types/water.jl:Master water container:1082"
    "3:src/biogeophys/FrictionVelocityMod.F90:src/types/friction_velocity.jl:Friction velocity:1200"

    # --- TIER 4: Biogeochemistry State Types ---
    "4:src/biogeochem/CNVegStateType.F90:src/types/cn_veg_state.jl:Veg CN state:600"
    "4:src/biogeochem/CNVegCarbonStateType.F90:src/types/cn_veg_carbon_state.jl:Veg carbon state:4876"
    "4:src/biogeochem/CNVegNitrogenStateType.F90:src/types/cn_veg_nitrogen_state.jl:Veg nitrogen state:2511"
    "4:src/biogeochem/CNVegCarbonFluxType.F90:src/types/cn_veg_carbon_flux.jl:Veg carbon flux:5566"
    "4:src/biogeochem/CNVegNitrogenFluxType.F90:src/types/cn_veg_nitrogen_flux.jl:Veg nitrogen flux:2755"
    "4:src/soilbiogeochem/SoilBiogeochemCarbonStateType.F90:src/types/soil_bgc_carbon_state.jl:Soil C state:1726"
    "4:src/soilbiogeochem/SoilBiogeochemNitrogenStateType.F90:src/types/soil_bgc_nitrogen_state.jl:Soil N state:1426"
    "4:src/soilbiogeochem/SoilBiogeochemCarbonFluxType.F90:src/types/soil_bgc_carbon_flux.jl:Soil C flux:1022"
    "4:src/soilbiogeochem/SoilBiogeochemNitrogenFluxType.F90:src/types/soil_bgc_nitrogen_flux.jl:Soil N flux:1253"
    "4:src/soilbiogeochem/SoilBiogeochemStateType.F90:src/types/soil_bgc_state.jl:Soil BGC state:500"
    "4:src/biogeochem/CropType.F90:src/types/crop.jl:Crop type:1017"
    "4:src/biogeochem/CNSharedParamsMod.F90:src/types/cn_shared_params.jl:CN shared params:300"

    # --- TIER 5: Radiation & Surface Physics ---
    "5:src/biogeophys/DaylengthMod.F90:src/biogeophys/daylength.jl:Day length:200"
    "5:src/biogeophys/SurfaceAlbedoMod.F90:src/biogeophys/surface_albedo.jl:Surface albedo:1781"
    "5:src/biogeophys/UrbanAlbedoMod.F90:src/biogeophys/urban_albedo.jl:Urban albedo:1305"
    "5:src/biogeophys/SurfaceRadiationMod.F90:src/biogeophys/surface_radiation.jl:Surface radiation:1027"
    "5:src/biogeophys/UrbanRadiationMod.F90:src/biogeophys/urban_radiation.jl:Urban radiation:500"
    "5:src/biogeophys/SnowSnicarMod.F90:src/biogeophys/snow_snicar.jl:Snow SNICAR:2264"
    "5:src/biogeophys/AerosolMod.F90:src/biogeophys/aerosol.jl:Aerosol:600"
    "5:src/biogeophys/SurfaceHumidityMod.F90:src/biogeophys/surface_humidity.jl:Surface humidity:200"
    "5:src/biogeophys/SurfaceResistanceMod.F90:src/biogeophys/surface_resistance.jl:Surface resistance:400"
    "5:src/biogeophys/SoilMoistStressMod.F90:src/biogeophys/soil_moist_stress.jl:Soil moisture stress:400"

    # --- TIER 6: Temperature & Hydrology Core ---
    "6:src/biogeophys/SoilTemperatureMod.F90:src/biogeophys/soil_temperature.jl:Soil temperature:2974"
    "6:src/biogeophys/LakeTemperatureMod.F90:src/biogeophys/lake_temperature.jl:Lake temperature:1482"
    "6:src/biogeophys/SoilWaterRetentionCurveMod.F90:src/biogeophys/swrc_base.jl:SWRC base:150"
    "6:src/biogeophys/SoilWaterRetentionCurveClappHornberg1978Mod.F90:src/biogeophys/swrc_clapp_hornberg.jl:SWRC Clapp-Hornberg:200"
    "6:src/biogeophys/SoilWaterRetentionCurveVanGenuchten1980Mod.F90:src/biogeophys/swrc_van_genuchten.jl:SWRC Van Genuchten:200"
    "6:src/biogeophys/SoilWaterMovementMod.F90:src/biogeophys/soil_water_movement.jl:Soil water movement:2210"
    "6:src/biogeophys/SoilHydrologyMod.F90:src/biogeophys/soil_hydrology.jl:Soil hydrology:2910"
    "6:src/biogeophys/SnowHydrologyMod.F90:src/biogeophys/snow_hydrology.jl:Snow hydrology:4070"
    "6:src/biogeophys/HydrologyNoDrainageMod.F90:src/biogeophys/hydrology_no_drainage.jl:Hydrology no drain:800"
    "6:src/biogeophys/HydrologyDrainageMod.F90:src/biogeophys/hydrology_drainage.jl:Hydrology drainage:700"
    "6:src/biogeophys/LakeHydrologyMod.F90:src/biogeophys/lake_hydrology.jl:Lake hydrology:800"
    "6:src/biogeophys/HillslopeHydrologyMod.F90:src/biogeophys/hillslope_hydrology.jl:Hillslope hydrology:1148"

    # --- TIER 7: Canopy & Flux Physics ---
    "7:src/biogeophys/PhotosynthesisMod.F90:src/biogeophys/photosynthesis.jl:Photosynthesis:5209"
    "7:src/biogeophys/CanopyFluxesMod.F90:src/biogeophys/canopy_fluxes.jl:Canopy fluxes:1756"
    "7:src/biogeophys/CanopyHydrologyMod.F90:src/biogeophys/canopy_hydrology.jl:Canopy hydrology:1191"
    "7:src/biogeophys/BareGroundFluxesMod.F90:src/biogeophys/bareground_fluxes.jl:Bare ground fluxes:500"
    "7:src/biogeophys/SoilFluxesMod.F90:src/biogeophys/soil_fluxes.jl:Soil fluxes:400"
    "7:src/biogeophys/UrbanFluxesMod.F90:src/biogeophys/urban_fluxes.jl:Urban fluxes:1143"
    "7:src/biogeophys/LunaMod.F90:src/biogeophys/luna.jl:LUNA:1410"
    "7:src/biogeophys/OzoneMod.F90:src/biogeophys/ozone.jl:Ozone damage:400"

    # --- TIER 8: Biogeochemistry Physics ---
    "8:src/biogeochem/CNPhenologyMod.F90:src/biogeochem/phenology.jl:Phenology:4517"
    "8:src/biogeochem/CNAllocationMod.F90:src/biogeochem/allocation.jl:CN allocation:800"
    "8:src/biogeochem/CNGRespMod.F90:src/biogeochem/growth_resp.jl:Growth respiration:400"
    "8:src/biogeochem/CNMRespMod.F90:src/biogeochem/maint_resp.jl:Maintenance respiration:500"
    "8:src/biogeochem/CNFUNMod.F90:src/biogeochem/fun.jl:FUN model:1800"
    "8:src/biogeochem/CNGapMortalityMod.F90:src/biogeochem/gap_mortality.jl:Gap mortality:400"
    "8:src/soilbiogeochem/SoilBiogeochemDecompCascadeBGCMod.F90:src/biogeochem/decomp_bgc.jl:Decomp BGC:800"
    "8:src/soilbiogeochem/SoilBiogeochemDecompCascadeMIMICSMod.F90:src/biogeochem/decomp_mimics.jl:Decomp MIMICS:1367"
    "8:src/soilbiogeochem/SoilBiogeochemDecompMod.F90:src/biogeochem/decomp.jl:Decomposition:600"
    "8:src/soilbiogeochem/SoilBiogeochemNitrifDenitrifMod.F90:src/biogeochem/nitrif_denitrif.jl:Nitrif-Denitrif:600"
    "8:src/soilbiogeochem/SoilBiogeochemNLeachingMod.F90:src/biogeochem/n_leaching.jl:N leaching:300"
    "8:src/biogeochem/CNNDynamicsMod.F90:src/biogeochem/n_dynamics.jl:N dynamics:500"
    "8:src/biogeochem/CNCStateUpdate1Mod.F90:src/biogeochem/c_state_update1.jl:C state update 1:400"
    "8:src/biogeochem/CNCStateUpdate2Mod.F90:src/biogeochem/c_state_update2.jl:C state update 2:400"
    "8:src/biogeochem/CNCStateUpdate3Mod.F90:src/biogeochem/c_state_update3.jl:C state update 3:400"
    "8:src/biogeochem/CNNStateUpdate1Mod.F90:src/biogeochem/n_state_update1.jl:N state update 1:400"
    "8:src/biogeochem/CNNStateUpdate2Mod.F90:src/biogeochem/n_state_update2.jl:N state update 2:400"
    "8:src/biogeochem/CNNStateUpdate3Mod.F90:src/biogeochem/n_state_update3.jl:N state update 3:400"
    "8:src/biogeochem/CNDriverMod.F90:src/biogeochem/cn_driver.jl:CN driver:1337"
    "8:src/biogeochem/CNVegetationFacade.F90:src/biogeochem/vegetation_facade.jl:Vegetation facade:1665"

    # --- TIER 9: Specialized Subsystems ---
    "9:src/biogeochem/CNFireBaseMod.F90:src/biogeochem/fire_base.jl:Fire base:1339"
    "9:src/biogeochem/CNFireLi2014Mod.F90:src/biogeochem/fire_li2014.jl:Fire Li2014:1500"
    "9:src/biogeochem/ch4Mod.F90:src/biogeochem/methane.jl:Methane:4450"
    "9:src/biogeochem/VOCEmissionMod.F90:src/biogeochem/voc_emission.jl:VOC emissions:1086"
    "9:src/biogeochem/DustEmisBase.F90:src/biogeochem/dust_emission.jl:Dust emission:600"
    "9:src/biogeophys/IrrigationMod.F90:src/biogeophys/irrigation.jl:Irrigation:1786"
    "9:src/biogeochem/SatellitePhenologyMod.F90:src/biogeochem/satellite_phenology.jl:Satellite phenology:400"
    "9:src/biogeochem/CNCIsoFluxMod.F90:src/biogeochem/c_iso_flux.jl:C isotope flux:1936"
    "9:src/biogeochem/CNBalanceCheckMod.F90:src/biogeochem/cn_balance_check.jl:CN balance check:600"
    "9:src/biogeophys/BalanceCheckMod.F90:src/biogeophys/balance_check.jl:Balance check:1053"

    # --- TIER 10: Infrastructure & Coupling ---
    "10:src/main/atm2lndType.F90:src/types/atm2lnd.jl:Atm-to-land type:1149"
    "10:src/main/lnd2atmType.F90:src/types/lnd2atm.jl:Land-to-atm type:500"
    "10:src/main/TopoMod.F90:src/infrastructure/topo.jl:Topography:400"
    "10:src/main/subgridAveMod.F90:src/infrastructure/subgrid_ave.jl:Subgrid averaging:1098"
    "10:src/main/accumulMod.F90:src/infrastructure/accumul.jl:Accumulation:400"
    "10:src/main/pftconMod.F90:src/constants/pftcon.jl:PFT constants:1602"
    "10:src/main/controlMod.F90:src/infrastructure/control.jl:Control module:1350"
    "10:src/main/clm_instMod.F90:src/infrastructure/instances.jl:Instance factory:800"

    # --- TIER 11: Driver ---
    "11:src/main/clm_driver.F90:src/driver/clm_driver.jl:Main CLM driver:1711"
)
