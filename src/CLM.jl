module CLM

using LinearAlgebra
using Dates
using NCDatasets

# ===========================================================================
# Tier 1: Constants & Parameters
# ===========================================================================
include("constants/varctl.jl")
include("constants/landunit_varcon.jl")
include("constants/column_varcon.jl")
include("constants/varpar.jl")
include("constants/varcon.jl")
include("constants/pftcon.jl")

# ===========================================================================
# Tier 1: Types (data structures)
# ===========================================================================
include("types/gridcell.jl")
include("types/landunit.jl")
include("types/column.jl")
include("types/patch.jl")
include("types/temperature.jl")
include("types/energy_flux.jl")
include("types/soil_state.jl")
include("types/soil_hydrology.jl")
include("types/canopy_state.jl")
include("types/lake_state.jl")
include("types/surface_albedo.jl")
include("types/solar_absorbed.jl")
include("types/urban_params.jl")
include("types/water_info_base.jl")
include("types/water_state.jl")
include("types/water_flux.jl")
include("types/water_diagnostic_bulk.jl")
include("types/water_state_bulk.jl")
include("types/water_flux_bulk.jl")
include("types/water_balance.jl")
include("types/water.jl")
include("types/friction_velocity.jl")
include("types/cn_veg_state.jl")
include("types/cn_veg_carbon_state.jl")
include("types/cn_veg_nitrogen_state.jl")
include("types/cn_veg_carbon_flux.jl")
include("types/cn_veg_nitrogen_flux.jl")
include("types/soil_bgc_carbon_state.jl")
include("types/soil_bgc_nitrogen_state.jl")
include("types/soil_bgc_carbon_flux.jl")
include("types/soil_bgc_nitrogen_flux.jl")
include("types/soil_bgc_state.jl")
include("types/crop.jl")
include("types/cn_shared_params.jl")
include("types/cn_products.jl")
include("types/atm2lnd.jl")
include("types/lnd2atm.jl")

# ===========================================================================
# Tier 1: Infrastructure (solvers, utilities, decomposition)
# ===========================================================================
include("infrastructure/decomp.jl")
include("infrastructure/filters.jl")
include("infrastructure/tridiagonal.jl")
include("infrastructure/band_diagonal.jl")
include("infrastructure/smooth_ad.jl")

# ===========================================================================
# Calibration overrides struct (needed by physics modules for override kwargs)
# ===========================================================================
include("calibration/overrides.jl")
include("infrastructure/topo.jl")
include("infrastructure/subgrid_ave.jl")
include("infrastructure/accumul.jl")
include("infrastructure/control.jl")

# ===========================================================================
# Initialization pipeline
# ===========================================================================
include("infrastructure/surfrd_utils.jl")
include("infrastructure/init_subgrid.jl")
include("infrastructure/subgrid_weights.jl")
include("infrastructure/time_manager.jl")
include("infrastructure/surfdata.jl")
include("infrastructure/init_gridcells.jl")
include("infrastructure/read_params.jl")
include("infrastructure/init_vertical.jl")
include("infrastructure/orbital.jl")

# ===========================================================================
# Tier 1: Biogeophysics (pure math)
# ===========================================================================
include("biogeophys/qsat.jl")
include("biogeophys/daylength.jl")
include("biogeophys/snow_snicar.jl")
include("biogeophys/aerosol.jl")
include("biogeophys/surface_albedo.jl")
include("biogeophys/urban_albedo.jl")
include("biogeophys/surface_radiation.jl")
include("biogeophys/urban_radiation.jl")
include("biogeophys/surface_humidity.jl")
include("biogeophys/surface_resistance.jl")
include("biogeophys/soil_moist_stress.jl")
include("biogeophys/soil_temperature.jl")
include("biogeophys/lake_temperature.jl")
include("biogeophys/swrc_base.jl")
include("biogeophys/swrc_clapp_hornberg.jl")
include("biogeophys/swrc_van_genuchten.jl")
include("biogeophys/soil_water_movement.jl")
include("biogeophys/soil_hydrology.jl")
include("biogeophys/snow_cover_fraction.jl")
include("biogeophys/snow_hydrology.jl")
include("biogeophys/hydrology_no_drainage.jl")
include("biogeophys/hydrology_drainage.jl")
include("biogeophys/lake_hydrology.jl")
include("biogeophys/hillslope_hydrology.jl")
include("biogeophys/photosynthesis.jl")
include("biogeophys/canopy_fluxes.jl")
include("biogeophys/canopy_hydrology.jl")
include("biogeophys/bareground_fluxes.jl")
include("biogeophys/soil_fluxes.jl")
include("biogeophys/urban_fluxes.jl")
include("biogeophys/luna.jl")
include("biogeophys/ozone.jl")
include("biogeophys/irrigation.jl")
include("biogeophys/balance_check.jl")
include("biogeophys/pre_flux_calcs.jl")
include("biogeophys/surface_water.jl")
include("biogeophys/lake_con.jl")
include("biogeophys/lake_fluxes.jl")
include("biogeophys/active_layer.jl")
include("biogeophys/root_biophys.jl")
include("biogeophys/soil_water_plant_sink.jl")
include("biogeophys/sat_excess_runoff.jl")
include("biogeophys/infilt_excess_runoff.jl")
include("biogeophys/total_water_heat.jl")
include("biogeophys/dry_dep_velocity.jl")

# ===========================================================================
# Tier 2: Biogeochemistry
# ===========================================================================
include("biogeochem/phenology.jl")
include("biogeochem/allocation.jl")
include("biogeochem/growth_resp.jl")
include("biogeochem/maint_resp.jl")
include("biogeochem/fun.jl")
include("biogeochem/gap_mortality.jl")
include("biogeochem/decomp_bgc.jl")
include("biogeochem/decomp_mimics.jl")
include("biogeochem/decomp.jl")
include("biogeochem/decomp_competition.jl")
include("biogeochem/decomp_vertical_profile.jl")
include("biogeochem/decomp_potential.jl")
include("biogeochem/decomp_precision_control.jl")
include("biogeochem/litter_vert_transp.jl")
include("biogeochem/nitrif_denitrif.jl")
include("biogeochem/n_leaching.jl")
include("biogeochem/n_dynamics.jl")
include("biogeochem/c_state_update1.jl")
include("biogeochem/c_state_update2.jl")
include("biogeochem/c_state_update3.jl")
include("biogeochem/n_state_update1.jl")
include("biogeochem/n_state_update2.jl")
include("biogeochem/n_state_update3.jl")
include("biogeochem/fire_base.jl")
include("biogeochem/fire_li2014.jl")
include("biogeochem/methane.jl")
include("biogeochem/voc_emission.jl")
include("biogeochem/dust_emission.jl")
include("biogeochem/satellite_phenology.jl")
include("biogeochem/c_iso_flux.jl")
include("biogeochem/cn_balance_check.jl")
include("biogeochem/cn_driver.jl")
include("biogeochem/cn_precision_control.jl")
include("biogeochem/veg_struct_update.jl")
include("biogeochem/nutrient_competition.jl")
include("biogeochem/veg_compute_seed.jl")
include("biogeochem/cn_annual_update.jl")
include("biogeochem/cn_products_mod.jl")
include("biogeochem/vegetation_facade.jl")

# ===========================================================================
# Instance factory (depends on all types above)
# ===========================================================================
include("infrastructure/instances.jl")

# ===========================================================================
# Infrastructure (depends on instances + all types)
# ===========================================================================
include("infrastructure/downscale_forcings.jl")
include("infrastructure/lnd2atm_mod.jl")
include("infrastructure/forcing_reader.jl")
include("infrastructure/history_writer.jl")
include("infrastructure/restart_io.jl")
include("infrastructure/fortran_restart.jl")

# ===========================================================================
# Initialization (depends on instances + all types)
# ===========================================================================
include("infrastructure/cold_start.jl")

# ===========================================================================
# Driver (depends on all modules above)
# ===========================================================================
include("driver/clm_driver.jl")
include("driver/clm_initialize.jl")
include("driver/clm_run.jl")

# ===========================================================================
# Calibration framework (depends on driver + all modules)
# ===========================================================================
include("calibration/calibration.jl")
include("calibration/parameters.jl")
include("calibration/optimize.jl")
include("calibration/enzyme_ad.jl")
include("calibration/fluxnet_reader.jl")
include("calibration/site_calibration.jl")
include("calibration/param_injection.jl")

end # module CLM
