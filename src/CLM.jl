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

# FATES (Tier F) Batch 5 — tree-damage module (DamageMainMod). Calendar-driven
# damage-event test (IsItDamageTime!), damage-class transition fraction lookup
# (GetDamageFrac via param_derived), crown-loss fraction (GetCrownReduction), and
# damage-class (crown-loss) dependent mortality (GetDamageMortality). Depends on
# the merged accessors ed_params()/edpftvarcon_inst()/param_derived() and the
# FatesInterfaceTypesMod calendar globals. Standalone — NOT added to CLMInstances
# or any dual-copied struct.
include("fates/DamageMainMod.jl")

# FATES (Tier F) Batch 6 — the allometry engine (FatesAllometryMod). The
# size->structure allometry library: height(dbh), crown area, AGB/BGB/sapwood/
# leaf/fineroot/storage/structure biomass, and each relationship's analytic
# derivative w.r.t. dbh. Dispatches among per-PFT allometry-mode selectors (h,
# bagw, blmax, bsap, bdead, bbgw, bfineroot, carea, store). Includes the
# diameter<->height inverses, root-profile builders, ForceDBH/size2dbh Newton+
# Regula-Falsi solvers, VegAreaLayer, and CheckIntegratedAllometries. Depends on
# the merged accessors prt_params/param_derived()/ed_params(), GetCrownReduction
# (DamageMainMod), and the EDParamsMod VAI bins. Standalone — NOT added to
# CLMInstances or any dual-copied struct.
include("fates/FatesAllometryMod.jl")

# Batch 7 — PARTEH carbon-only allometric allocation hypothesis
# (PRTAllometricCarbonMod). Concrete subtype of AbstractPRTVartypes that
# implements the deferred DailyPRT!/FastPRT! generics; depends on the allometry
# engine + integrators above. Standalone — NOT added to CLMInstances.
include("fates/PRTAllometricCarbonMod.jl")
# FATES (Tier F) Batch 7 — the Carbon-Nitrogen-Phosphorus (CNP) prioritized
# allometric allocation hypothesis (PRTAllometricCNPMod). The most complex PARTEH
# allocation strategy: a concrete `cnp_allom_prt_vartypes <: AbstractPRTVartypes`
# implementing the deferred generics DailyPRT!/FastPRT!/GetNutrientTarget/
# GetCoordVal plus the CNP-specific allocation chain (prioritized replacement →
# stature growth via Euler integration of allometric derivatives → remainder
# allocation/efflux), the multi-element stoichiometric targets, nutrient-limited
# downregulation, the PID l2fr controller, and storage overflow/efflux/burn.
# Depends on the merged PRTGenericMod (abstract type + generics + registration),
# FatesAllometryMod ((value,deriv) routines), prt_params, ed_params(), and the
# Euler/RKF45 integrators. Standalone — NOT added to CLMInstances.
include("fates/PRTAllometricCNPMod.jl")

# FATES (Tier F) Batch 8 — the COHORT type (FatesCohortMod). The fundamental
# FATES demographic unit: integrates the PRT allocation state, the per-cohort
# plant-hydraulics state, and all per-cohort vegetation/flux/mortality/phenology
# variables, plus the taller/shorter linked-list pointers. Depends on the PRT
# generic + concrete carbon/CNP hypotheses, hydraulics memory, allometry,
# size/age indices, the FATES PFT/derived parameter tables, and the interface
# control flags above. Standalone — NOT added to CLMInstances.
include("fates/FatesCohortMod.jl")

# FATES (Tier F) Batch 9 — the PATCH type (FatesPatchMod). A FATES patch is a
# collection of cohorts sharing a disturbance history/age. It owns the cohort
# linked list (tallest/shortest), the older/younger patch linked list, the
# canopy-layer/PFT/leaf-layer leaf-area + radiation profile arrays, the per-element
# litter pools, the SPITFIRE fuel object, the two-stream radiation object, and the
# per-patch running means. Depends on the cohort type (Batch 8) + litter/fuel/
# two-stream/running-mean types + PRT generic + EDParams/Radiation/Interface
# constants. Standalone — NOT added to CLMInstances.
include("fates/FatesPatchMod.jl")

# FATES (Tier F) Batch 10 — the SITE type (EDTypesMod), the TOP of the FATES
# type system and the "type-system-complete" checkpoint. ed_site_type is the
# FATES site/column structure that owns the age-ordered patch linked list
# (oldest_patch/youngest_patch), the site plant-hydraulics object, the
# fire-weather object, the per-element mass-balance + integrated-flux-balance
# arrays, the flux diagnostics, the full phenology state, soil layering, and the
# termination/recruitment/demotion/promotion/disturbance accumulators on
# size x pft history arrays. Also defines the site-scoped helper types
# (ed_resources_management_type, elem_diag_type, site_ifluxbal_type,
# site_fluxdiags_type, site_massbal_type) + the module-level FATES constants.
# Depends on the patch type (Batch 9) + cohort type + ed_site_hydr_type +
# fire_weather + litter/running-mean/PRT-generic types + EDParams/Interface
# constants. Standalone — NOT added to CLMInstances.
include("fates/EDTypesMod.jl")

# FATES (Tier F) Batch 11 — the PARTEH parameter initialization / registration
# layer (PRTParamsFATESMod, Fortran module PRTInitParamsFatesMod). Registers and
# receives all the PARTEH / PFT allocation / stoichiometry / turnover parameters
# into the merged `prt_params` via the FATES parameter reader, plus the PARTEH
# consistency checks (PRTCheckParams) + derived-parameter setup (PRTDerivedParams)
# + the new-recruit total stoichiometry helper. Depends on the merged param
# storage (PRTParametersMod), the parameter-reader API (FatesParametersInterface),
# the organ/element index consts + hypothesis consts + StorageNutrientTarget
# (PRTGenericMod), the allometry chain (FatesAllometryMod), EDPftvarcon_inst, and
# init_recruit_trim (EDTypesMod). Standalone — nothing added to CLMInstances.
include("fates/PRTParamsFATESMod.jl")
# FATES (Tier F) Batch 11 — two standalone diagnostic / accumulation modules
# that build on the complete merged FATES type system above.
#   * ChecksBalancesMod — site-level carbon/mass-balance checking + diagnostics:
#     SiteMassStock/PatchMassStock sum the biomass+seed+litter stocks per element,
#     CheckLitterPools sanity-checks for negative litter, and
#     CheckIntegratedMassPools accumulates the time-integrated net fluxes in/out
#     of vegetation + litter and compares against the instantaneous state.
#   * EDAccumulateFluxesMod — AccumulateFluxes_ED accumulates/averages the
#     per-cohort photosynthesis/respiration fluxes (npp/gpp/resp/sym_nfix/c13/
#     year_net_uptake) over the FATES daily timestep.
# Both depend on ed_site_type/patch/cohort + the EDTypes mass-balance helper
# types + PRTGenericMod element/organ consts + FatesLitter. Standalone — NOT
# added to CLMInstances or any dual-copied struct.
include("fates/ChecksBalancesMod.jl")
include("fates/EDAccumulateFluxesMod.jl")
# FATES (Tier F) Batch 11 — the FATES <-> host soil-biogeochem flux coupling
# (FatesSoilBGCFluxMod). Prepares the litter/CWD fragmentation + root-exudate
# fluxes FATES passes to the host soil-BGC (partitioned into cellulose/lignin/
# labile chemical fractions over decomposition layers), the per-competitor fine-
# root carbon + decomposer microbial biomass nutrient-acquisition boundary
# conditions, and the CH4-model boundary conditions; and unpacks the host's
# returned per-competitor N (NH4/NO3) and P uptake fluxes back to each cohort's
# daily demand/uptake. Operates on the merged site/patch/cohort types + the
# per-element litter_type + bc_in_type/bc_out_type. Standalone — NOT added to
# CLMInstances.
include("fates/FatesSoilBGCFluxMod.jl")
# FATES plant-hydraulics solver (FatesPlantHydraulicsMod). Plant water transport
# (Richards-style implicit 1D Taylor solver), stomatal water-stress (btran), and
# the soil->root->stem->leaf flow network + rhizosphere shells + water mass
# balance. Default solver path (hydr_solver_1DTaylor) ported fully; the two
# non-default 2D solvers (Newton/Picard) are stubbed. Depends on the FATES type
# system + WRF/WKF functions + allometry + PRT + params (all included above).
# Standalone — NOT added to CLMInstances.
include("fates/FatesPlantHydraulicsMod.jl")
# FATES (Tier F) Batch 11 — the SPITFIRE main fire driver (SFMainMod). The daily
# fire model: per-patch fuel characterization -> fire weather / danger index ->
# rate of spread (Rothermel 1972 / Thonicke et al. 2010) -> ground fuel
# consumption -> fire intensity & area burnt -> fire effects (crown scorch /
# damage, cambial damage / kill, post-fire tree mortality). Operates on the
# merged ed_site / patch / cohort + fuel_type + fire_weather objects. Standalone
# — NOT added to CLMInstances or any dual-copied struct. Must come AFTER
# EDTypesMod (site type) + FatesFuelMod + SFParamsMod + EDPftvarcon + the
# allometry engine (CrownDepth), all included above.
include("fates/SFMainMod.jl")
# FATES (Tier F) Batch 11 — two independent modules built on the complete FATES
# type system. (1) FatesTwoStreamUtilsMod: the glue between the FATES canopy
# (patch/cohort canopy-layer x PFT structure) and the merged TwoStreamMLPEMod
# solver — builds the scattering-element inputs from the patch canopy profiles,
# calls the solver, and maps absorbed radiation back onto cohorts/leaf-layers
# (FatesConstructRadElements/FatesGetCohortAbsRad/FatesPatchFSun/
# CheckPatchRadiationBalance/TransferRadParams). (2) FatesLandUseChangeMod: FATES
# land-use-change transitions — aggregates the LUH2 transition/state vectors to
# the 5 FATES land use categories, builds the [m2/m2/day] transition matrix +
# clearing-rules matrix, and the spin-up→land-use initialization transitions.
# Both depend on the patch/cohort/site types + two-stream solver / interface
# constants. Standalone — NOT added to CLMInstances.
include("fates/FatesTwoStreamUtilsMod.jl")
include("fates/FatesLandUseChangeMod.jl")
# FATES (Tier F) Batch 12 — three demographic/biogeophys modules on the complete
# type system. (1) EDBtranMod: the FATES water-stress factor btran + root-soil
# water uptake distribution (get_active_suction_layers / btran_ed). (2)
# EDCohortDynamicsMod: the demographic cohort engine — create/insert/sort/count/
# fuse/terminate of the self-referential cohort linked list, EvaluateAndCorrectDBH
# + DamageRecovery, routing terminated-cohort biomass to litter. (3)
# EDLoggingMortalityMod: logging/harvest disturbance mortality — direct/collateral/
# mechanical mortality fractions per PFT/size/scenario, killed-tree biomass → CWD/
# litter + exported wood product, harvest-rate/harvest-debt logic. All depend on
# the cohort/patch/site types, litter pools, allometry, PRT, EDParams/EDPftvarcon,
# and (logging) FatesLandUseChangeMod above. Standalone — NOT added to CLMInstances.
include("fates/EDBtranMod.jl")
include("fates/EDCohortDynamicsMod.jl")
include("fates/EDLoggingMortalityMod.jl")
# FATES (Tier F) Batch 13 — the physiology checkpoint. (1) FatesBstressMod:
# salinity transpiration-stress factor (bstress_sal_ft) — sibling of EDBtranMod.
# (2) EDMortalityFunctionsMod: per-cohort mortality rates (carbon-starvation /
# hydraulic-failure / background / cold-freezing / size-senescence / damage) +
# Mortality_Derivative feeding the cohort number-density ODE + logging coupling.
# (3) EDPhysiologyMod: the physiology hub — cold/drought-deciduous phenology +
# leaf on/off, trim_canopy (LAI→carbon-balance trimming), recruitment / seed
# rain / germination / decay, CWD & litter input/fragmentation/output, satellite
# (SP) phenology, recruit L2FR/stoich. All depend on the cohort/patch/site types,
# allometry, PRT, litter pools, EDParams/EDPftvarcon above. Plant-hydraulics
# init paths stubbed behind hlm_use_planthydro (off by default). Standalone —
# NOT added to CLMInstances.
include("fates/FatesBstressMod.jl")
include("fates/EDMortalityFunctionsMod.jl")
include("fates/EDPhysiologyMod.jl")
# FATES (Tier F) Batch 14 — the patch disturbance engine (EDPatchDynamicsMod):
# per-patch disturbance_rates (fire/treefall-mortality/logging/land-use) →
# spawn_patches/split_patch (move surviving & killed cohorts into new disturbed
# patches, conserving area + cohort number + mass) → the three *_litter_fluxes
# routines (fire/mortality/landusechange — killed/disturbed biomass → litter/CWD
# + burn flux) → fuse_patches/fuse_2_patches (area-weighted profile merge) →
# terminate_patches → DistributeSeeds + patch census (set_patchno/countPatches/
# check_patch_area/patch_pft_size_profile/GetPseudoPatchAge). Reuses the Batch
# 12–13 cohort/mortality/logging/physiology helpers. Plant-hydraulics init paths
# stubbed behind hlm_use_planthydro (off by default). Standalone — NOT added to
# CLMInstances.
include("fates/EDPatchDynamicsMod.jl")
# FATES (Tier F) Batch 15 — the canopy STRUCTURE engine (EDCanopyStructureMod):
# canopy_structure! arranges cohorts into discrete canopy layers (PPA) via
# DemoteFromLayer!/PromoteIntoLayer! so each layer's crown area fits the patch;
# canopy_spread! (crown-area scaling); canopy_summarization!/leaf_area_profile!
# build the per-patch canopy-layer × PFT × leaf-layer LAI/SAI profiles;
# update_hlm_dynamics! packs canopy state (LAI/SAI, fractions, z0m/displa/dleaf)
# to the host. Reuses Batch 12–14 cohort/patch/allometry helpers. Plant-hydraulics
# paths stubbed behind hlm_use_planthydro (off by default). Standalone — NOT added
# to CLMInstances.
include("fates/EDCanopyStructureMod.jl")
# FATES (Tier F) Batch 16 — radiation + photosynthesis + inventory-init + the
# ecosystem-dynamics DRIVER. (1) FatesNormanRadMod: the Norman (1979) multi-layer
# canopy radiative-transfer solver (PatchNormanRadiation — per canopy-layer × PFT
# × leaf-layer absorbed/reflected/transmitted direct+diffuse rad, albedo,
# sunlit/shaded fractions). (2) FatesPlantRespPhotosynthMod: leaf-scale
# Farquhar/Collatz photosynthesis (C3/C4) + Ball-Berry/Medlyn stomatal
# conductance + Ryan-1991/Atkin-2017 maintenance respiration + leaf→cohort flux
# scaling (FatesPlantRespPhotosynthDrive). (3) FatesInventoryInitMod: PSS/CSS
# inventory-based patch/cohort site initialization (initialize_sites_by_inventory).
# (4) EDMainMod: the daily FATES driver (ed_ecosystem_dynamics) orchestrating
# phenology→growth→mortality→disturbance→canopy structure→recruitment, plus
# ed_update_site/TotalBalanceCheck/bypass_dynamics — included LAST as it reuses
# all the above. Plant-hydraulics paths stubbed behind hlm_use_planthydro (off by
# default). Standalone — NOT added to CLMInstances.
include("fates/FatesNormanRadMod.jl")
include("fates/FatesPlantRespPhotosynthMod.jl")
include("fates/FatesInventoryInitMod.jl")
include("fates/EDMainMod.jl")
# FATES (Tier F) Batch 17 — the init checkpoint. (1) FatesRadiationDriveMod: the
# host-facing canopy radiation driver — FatesNormalizedCanopyRadiation (selects
# Norman vs two-stream, calls the per-patch solver, packs albd/albi/fabd/fabi/
# ftdd/ftid/ftii into bc_out) + FatesSunShadeFracs (sunlit/shaded LAI + absorbed
# PAR for the photosynthesis driver); reuses PatchNormanRadiation + the two-stream
# solver. (2) EDInitMod: FATES default cold-start state init — init_site_vars!/
# zero_site!/set_site_properties! + init_patches!/init_cohorts! (near-bareground
# demographic seeding, the counterpart to FatesInventoryInitMod). Reuse the Batch
# 12–16 constructors/solvers. Plant-hydraulics + LUH paths stubbed behind
# hlm_use_planthydro/hlm_use_luh (off by default). Standalone — NOT added to
# CLMInstances.
include("fates/FatesRadiationDriveMod.jl")
include("fates/EDInitMod.jl")
# FATES (Tier F) Batch 18 — host coupling + I/O (the final FATES wave). (1) The
# HISTORY pair: FatesHistoryVariableType (one output variable: metadata + bound
# data) + FatesHistoryInterfaceMod (the registry of all 472 FATES history vars +
# dimension bookkeeping + the update_history_* aggregation of site/patch/cohort
# state into output buffers; update_history_dyn1! ported, the remaining buffer-
# fill routines registry-complete but fill-logic deferred to a B18 followup —
# unexercised until FATES is driver-wired). (2) The RESTART pair:
# FatesRestartVariableType + FatesRestartInterfaceMod (the registry + the
# site→patch→cohort ↔ flat-vector serialization; the demographic round-trip is
# complete + tested, per-cohort diagnostic copies deferred). (3) FatesInterfaceMod:
# the CLM↔FATES coupling seam (set_fates_ctrlparms, SetFatesGlobalElements,
# allocate/zero/set bc_in/bc_out/bc_pconst, SetFatesTime, DetermineGridCellNeighbors).
# Var-types precede their interfaces; FatesInterfaceMod last. Standalone — NOT yet
# added to CLMInstances (live-driver wiring is the next milestone).
include("fates/FatesHistoryVariableType.jl")
include("fates/FatesHistoryInterfaceMod.jl")
include("fates/FatesRestartVariableType.jl")
include("fates/FatesRestartInterfaceMod.jl")
include("fates/FatesInterfaceMod.jl")
# W1+W2 live-driver wiring: clm_fates_init! — bootstrap + cold-start a single
# carbon-only FATES site and attach it to CLMInstances.fates. Depends on
# CLMInstances (instances.jl) + the whole FATES module stack above.
include("fates/fates_driver_init.jl")

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
