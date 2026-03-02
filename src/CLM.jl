module CLM

using LinearAlgebra

# ===========================================================================
# Tier 1: Constants & Parameters
# ===========================================================================
include("constants/varctl.jl")
include("constants/landunit_varcon.jl")
include("constants/column_varcon.jl")
include("constants/varpar.jl")
include("constants/varcon.jl")

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

# ===========================================================================
# Tier 1: Infrastructure (solvers, utilities, decomposition)
# ===========================================================================
include("infrastructure/decomp.jl")
include("infrastructure/filters.jl")
include("infrastructure/tridiagonal.jl")
include("infrastructure/band_diagonal.jl")

# ===========================================================================
# Tier 1: Biogeophysics (pure math)
# ===========================================================================
include("biogeophys/qsat.jl")
include("biogeophys/daylength.jl")
include("biogeophys/surface_albedo.jl")
include("biogeophys/urban_albedo.jl")
include("biogeophys/surface_radiation.jl")
include("biogeophys/urban_radiation.jl")
include("biogeophys/snow_snicar.jl")
include("biogeophys/aerosol.jl")
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
include("biogeophys/snow_hydrology.jl")
include("biogeophys/hydrology_no_drainage.jl")
include("biogeophys/hydrology_drainage.jl")
include("biogeophys/lake_hydrology.jl")
include("biogeophys/hillslope_hydrology.jl")
include("biogeophys/photosynthesis.jl")
include("biogeophys/canopy_fluxes.jl")
include("biogeophys/canopy_hydrology.jl")

end # module CLM
