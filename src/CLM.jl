module CLM

using LinearAlgebra
using Dates
using NCDatasets
import Adapt  # device-movable state structs (Adapt.@adapt_structure)
import KernelAbstractions  # backend-agnostic physics kernels (CPU/GPU)
import Atomix  # atomic scatter (patch→column accumulation) in kernels
# Bring the KA macros into module scope HERE (before any include) so kernels can be
# defined in early-included files too (e.g. types/friction_velocity.jl), not just in
# infrastructure/kernels.jl. Macros must be in scope at each file's include/lower time.
using KernelAbstractions: @kernel, @index, @Const

# ===========================================================================
# Tier 1: Constants & Parameters
# ===========================================================================
include("constants/precision.jl")
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
include("types/dgvs.jl")
# accumul.jl defines AccumManager, needed by temperature.jl's crop-GDD routine
# type signature, so it must be included before the types that reference it.
include("infrastructure/accumul.jl")
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
include("infrastructure/kernels.jl")  # KernelAbstractions physics kernels (Phase 4)

# ===========================================================================
# Calibration overrides struct (needed by physics modules for override kwargs)
# ===========================================================================
include("calibration/overrides.jl")
include("infrastructure/topo.jl")
include("infrastructure/subgrid_ave.jl")
include("infrastructure/control.jl")
include("infrastructure/dyn_subgrid_control.jl")

# ===========================================================================
# Initialization pipeline
# ===========================================================================
include("infrastructure/surfrd_utils.jl")
include("infrastructure/init_subgrid.jl")
include("infrastructure/subgrid_weights.jl")
include("infrastructure/dyn_landunit_area.jl")
include("infrastructure/dyn_patch_state_updater.jl")  # conservative patch-state update on weight change (dyn_subgrid)
include("infrastructure/dyn_init_columns.jl")
include("infrastructure/dyn_column_state_updater.jl")  # conservative column-state updates on transient column-area change (dyn_subgrid)
include("infrastructure/time_manager.jl")
include("infrastructure/surfdata.jl")
include("infrastructure/urban_input.jl")
include("infrastructure/init_gridcells.jl")
include("infrastructure/read_params.jl")
include("infrastructure/init_vertical.jl")
include("infrastructure/orbital.jl")
include("infrastructure/dyn_file_io.jl")  # transient land-use file time-stepping (dyn_subgrid)
include("infrastructure/dyn_pft_crop_file.jl")  # transient PCT_NAT_PFT / PCT_CROP readers (dyn_subgrid)
include("infrastructure/dyn_lake_urban_file.jl")  # transient lake & urban land-use readers (dyn_subgrid)

# ===========================================================================
# Tier 1: Biogeophysics (pure math)
# ===========================================================================
include("biogeophys/qsat.jl")
include("biogeophys/daylength.jl")
include("biogeophys/snow_snicar.jl")
include("biogeophys/snow_snicar_device.jl")
include("biogeophys/aerosol.jl")
include("biogeophys/surface_albedo.jl")
include("biogeophys/urban_albedo.jl")
include("biogeophys/surface_radiation.jl")
include("biogeophys/urban_radiation.jl")
include("biogeophys/surface_humidity.jl")
include("biogeophys/surface_resistance.jl")
include("biogeophys/soil_moist_stress.jl")
include("biogeophys/soil_temperature.jl")
include("biogeophys/urb_build_temp_oleson2015.jl")
include("biogeophys/lake_temperature.jl")
include("biogeophys/swrc_base.jl")
include("biogeophys/swrc_clapp_hornberg.jl")
include("biogeophys/swrc_van_genuchten.jl")
include("biogeophys/soil_water_movement.jl")
include("biogeophys/soil_hydrology.jl")
include("biogeophys/snow_cover_fraction.jl")
include("biogeophys/snow_hydrology.jl")
include("biogeophys/hydrology_no_drainage.jl")
include("biogeophys/glacier_surface_mass_balance.jl")
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
# dyn_cons_biogeophys must come AFTER dyn_subgrid_control.jl (line 87) and
# total_water_heat.jl (above), on which it depends.
include("infrastructure/dyn_cons_biogeophys.jl")

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
# dyn_cons_biogeochem (C/N conservation across land-cover change) depends on the
# CN veg/soil C+N state+flux types, the dyn patch/column state updaters
# (infrastructure, included above), and veg_compute_seed.jl.
include("biogeochem/dyn_cons_biogeochem.jl")
include("biogeochem/cn_annual_update.jl")
include("biogeochem/cn_products_mod.jl")
# Transient wood-harvest mortality (dyn_subgrid/dynHarvestMod). Needs dyn_file_io.jl
# (line 106) for the DynFile/DynVarTimeUninterp readers and the CN veg flux/state
# types above for the harvest C/N flux fields.
include("biogeochem/dyn_harvest.jl")
include("biogeochem/cndv.jl")
include("biogeochem/dyn_gross_unrep.jl")  # gross unrepresented LULCC reader + CN disturbance fluxes (dyn_subgrid)
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
include("infrastructure/history_io.jl")
include("infrastructure/restart_io.jl")
include("infrastructure/init_interp.jl")
include("infrastructure/fortran_restart.jl")

# ===========================================================================
# Initialization (depends on instances + all types)
# ===========================================================================
include("infrastructure/cold_start.jl")

# ===========================================================================
# FATES (Tier F) — Functionally Assembled Terrestrial Ecosystem Simulator
# ===========================================================================
# Batch 1: the foundation layer (bedrock for all other FATES modules). Ported
# from src/fates/main/. Standalone FATES types — do NOT add any of these to
# CLMInstances or any ForwardDiff-dual-copied struct. Included in dependency
# order: Constants -> Globals -> Integrators/Utils -> RunningMean ->
# IODimensions -> IOVariableKind -> ParametersInterface -> SynchronizedParams.
include("fates/FatesConstantsMod.jl")
include("fates/FatesGlobals.jl")
include("fates/FatesIntegratorsMod.jl")
include("fates/FatesUtilsMod.jl")
include("fates/FatesRunningMeanMod.jl")
include("fates/FatesIODimensionsMod.jl")
include("fates/FatesIOVariableKindMod.jl")
include("fates/FatesParametersInterface.jl")
include("fates/FatesSynchronizedParamsMod.jl")
include("fates/EDParamsMod.jl")

# Batch 2 (depends only on the foundation above): PARTEH/allometry parameter
# storage and the SPITFIRE fire-weather base type. Standalone — NOT added to
# CLMInstances or any dual-copied struct.
include("fates/PRTParametersMod.jl")
include("fates/SFFireWeatherMod.jl")

# Batch 2 (PARTEH generic framework): the abstract prt_vartypes/prt_global_type
# machinery, organ/element indexing, the generic prt_vartype state container, and
# the generic getters/setters/initializers. Depends only on the foundation +
# PRTParametersMod. Standalone — NOT added to CLMInstances or any dual-copied struct.
include("fates/PRTGenericMod.jl")

# Batch 3 (PARTEH loss fluxes): the loss-flux routines that operate on the
# generic prt_vartypes state — event turnover (deciduous drop / burn / damage /
# repro release), maintenance/background turnover with retranslocation to
# storage, and bud-burst flush from storage. Depends on the foundation +
# PRTGenericMod + PRTParametersMod. Standalone — NOT added to CLMInstances or
# any dual-copied struct.
include("fates/PRTLossFluxesMod.jl")

# Batch 1 (cont.): standalone leaf modules depending only on the foundation.
# FatesLitterMod (litter_type: CWD + fine-litter pools by element) and
# FatesRadiationMemMod (solar-band indices/params). Standalone — do NOT add to
# CLMInstances or any dual-copied struct.
include("fates/FatesLitterMod.jl")
include("fates/FatesRadiationMemMod.jl")
include("fates/FatesInterfaceTypesMod.jl")

# Batch 2: pure index-mapping helpers for the FATES history/diagnostic
# multiplexed dimensions (size/age/height/coage/fuel/damage class lookups +
# flat 1-D index combiners). Depends only on the foundation + EDParamsMod
# (nclmax, bin-edge vectors on ed_params()) + FatesInterfaceTypesMod (nlev*).
# Standalone — NOT added to CLMInstances or any dual-copied struct.
include("fates/FatesSizeAgeTypeIndicesMod.jl")

# Batch 1 (biogeophys): plant-hydraulics water transfer functions (WTFs). Depends
# only on the foundation (Constants/Globals). Standalone — not in CLMInstances.
include("fates/FatesHydroWTFMod.jl")

# Batch 2 (biogeophys): plant-hydraulics memory/state types (ed_site_hydr_type,
# ed_cohort_hydr_type) + compartment/shell indexing constants. Depends on the
# foundation + FatesHydroWTFMod (WRFType/WKFType). Standalone — not in CLMInstances.
include("fates/FatesHydraulicsMemMod.jl")

# Batch 1 (radiation): the multi-layer, multi-PFT two-stream canopy radiative
# transfer solver. Self-contained — depends only on the foundation above.
include("fates/TwoStreamMLPEMod.jl")

# Batch 2 (fire): SPITFIRE fuel-class enumeration (FatesFuelClassesMod, depends
# on FatesLitterMod) and the Nesterov fire-weather-index concrete subtype
# (SFNesterovMod, extends SFFireWeatherMod's fire_weather). Standalone — NOT
# added to CLMInstances or any dual-copied struct.
include("fates/FatesFuelClassesMod.jl")
include("fates/SFNesterovMod.jl")
# Batch 3: SPITFIRE fuel state (FatesFuelMod) — the `fuel_type` holding per-fuel-
# class loading + derived non-trunk loading/moisture/bulk-density/SAV/MEF.
# Depends on FatesFuelClassesMod, SFFireWeatherMod, SFNesterovMod, FatesConstantsMod.
include("fates/FatesFuelMod.jl")

# Batch 3: the FATES PFT-indexed parameter container (EDPftvarcon). The large
# per-PFT trait/allometry/allocation/mortality/phenology/fire/hydraulics
# parameter table + its registration/retrieval with the FATES parameter reader,
# consistency checks (FatesCheckParams), and litter-decomposability lookup
# (GetDecompyFrac). FATES's OWN PFT param table (EDPftvarcon_inst) — distinct
# from the host CLM `pftcon`. Standalone — NOT added to CLMInstances or any
# dual-copied struct. Depends on the foundation + Batches 1-2 above.
include("fates/EDPftvarcon.jl")
# Batch 3 (fire): SPITFIRE fire-model parameter holder (SFParamsMod) — fuel
# energy / mineral fractions / moisture-of-extinction + drying ratios /
# fire-spread + intensity coefficients by fuel class, with registration /
# retrieval against the FATES parameter reader. Depends on the foundation +
# FatesFuelClassesMod (num_fuel_classes) + FatesLitterMod (ncwd). Standalone —
# NOT added to CLMInstances or any dual-copied struct.
include("fates/SFParamsMod.jl")

# Batch 4: statically-derived FATES parameters + seed dispersal. Both standalone
# and independent of each other; build on the foundation + Batches 1-3.
# FatesParameterDerivedMod: parameters DERIVED at init from the raw
# EDPftvarcon/EDParams/SFParams (per-PFT photosynthesis constants from
# vcmax25top, branch fraction from CWD fractions, damage-class transition
# matrices). FatesDispersalMod: seed dispersal across grid cells — neighbor
# linked-list topology, the dispersal-kernel probability densities
# (exponential / exppower / logsech), and the calendar-driven IsItDispersalTime.
# Depend on the merged accessors edpftvarcon_inst()/ed_params()/sf_params() and
# the FatesInterfaceTypesMod module globals. Standalone — NOT added to
# CLMInstances or any dual-copied struct.
include("fates/FatesParameterDerivedMod.jl")
include("fates/FatesDispersalMod.jl")

# ===========================================================================
# Driver (depends on all modules above)
# ===========================================================================
# Transient land-use top-level orchestrator (dyn_subgrid/dynSubgridDriverMod).
# Must come AFTER all dyn_subgrid modules above (control/file_io/pft_crop/lake_urban/
# landunit_area/state-updaters/init_columns/cons_biogeophys + harvest/gross_unrep).
include("driver/dyn_subgrid_driver.jl")
include("driver/clm_driver.jl")
include("driver/clm_initialize.jl")
include("driver/clm_run.jl")
include("driver/run_clm.jl")

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

# Enzyme reverse-mode rules (band_solve! adjoint + inactive param-struct types).
# Included last so all referenced types (param/control containers) are defined.
include("infrastructure/enzyme_rules.jl")

# Compositional reverse-mode AD for canopy_fluxes_core! (per-sub-phase Enzyme calls,
# checkpointed Newton loop). Needs the canopy/photosynthesis kernels + Enzyme above.
include("biogeophys/canopy_fluxes_reverse.jl")

# Productionized clm_drv! reverse phases (soil_temperature!/soil_water!/water_table!)
# for the compositional engine — the building blocks of a whole-timestep reverse.
include("driver/driver_reverse.jl")

end # module CLM
