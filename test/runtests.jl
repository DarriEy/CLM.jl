using Test
using CLM

@testset "CLM.jl" begin
    include("test_constants.jl")
    include("test_driver_defaults_audit.jl")  # driver/varctl defaults pinned vs CTSM (docs/DRIVER_DEFAULTS_AUDIT.md)
    include("test_qsat.jl")
    include("test_tridiagonal.jl")
    include("test_decomp.jl")
    include("test_decomp_init.jl")
    include("test_threaded_clumps.jl")
    include("test_spmd.jl")
    include("test_distributed_driver.jl")  # MPI driver single-rank no-op path (2-rank smoke: test/mpi/run_mpi_smoke.jl under mpiexec)
    include("test_gridcell.jl")
    include("test_landunit.jl")
    include("test_column.jl")
    include("test_patch.jl")
    include("test_filters.jl")
    include("test_temperature.jl")
    include("test_energy_flux.jl")
    include("test_soil_state.jl")
    include("test_soil_hydrology.jl")
    include("test_canopy_state.jl")
    include("test_lake_state.jl")
    include("test_surface_albedo.jl")
    include("test_solar_absorbed.jl")
    include("test_urban_params.jl")
    include("test_read_params.jl")
    include("test_water_info_base.jl")
    include("test_water_state.jl")
    include("test_water_flux.jl")
    include("test_water_diagnostic_bulk.jl")
    include("test_water_state_bulk.jl")
    include("test_water_flux_bulk.jl")
    include("test_water_balance.jl")
    include("test_water.jl")
    include("test_friction_velocity.jl")
    include("test_cn_veg_state.jl")
    include("test_cn_veg_carbon_state.jl")
    include("test_cn_veg_nitrogen_state.jl")
    include("test_cn_veg_carbon_flux.jl")
    include("test_cn_offset_coldstart.jl")
    include("test_cn_veg_nitrogen_flux.jl")
    include("test_soil_bgc_carbon_state.jl")
    include("test_soil_bgc_nitrogen_state.jl")
    include("test_soil_bgc_carbon_flux.jl")
    include("test_soil_bgc_nitrogen_flux.jl")
    include("test_soil_bgc_state.jl")
    include("test_crop.jl")
    include("test_cn_shared_params.jl")
    include("test_daylength.jl")
    include("test_surface_albedo_mod.jl")
    include("test_urban_albedo.jl")
    include("test_surface_radiation.jl")
    include("test_urban_radiation.jl")
    include("test_snow_snicar.jl")
    include("test_aerosol.jl")
    include("test_surface_humidity.jl")
    include("test_surface_resistance.jl")
    include("test_soil_moist_stress.jl")
    include("test_btran_smoothing.jl")
    include("test_soil_temperature.jl")
    include("test_urban_building_temp.jl")
    include("test_lake_temperature.jl")
    include("test_lake_suboptions.jl")
    include("test_swrc_base.jl")
    include("test_swrc_clapp_hornberg.jl")
    include("test_swrc_van_genuchten.jl")
    include("test_soil_water_movement.jl")
    include("test_soil_hydrology_mod.jl")
    include("test_soil_lateral_flow.jl")
    include("test_snow_hydrology.jl")
    include("test_phase_snow_smoothing.jl")
    include("test_smooth_conservation.jl")
    include("test_smooth_axis_guard.jl")
    include("test_snow_capping.jl")
    include("test_snow_cover_fraction.jl")
    include("test_hydrology_no_drainage.jl")
    include("test_hydrology_drainage.jl")
    include("test_glacier_smb.jl")
    include("test_lake_hydrology.jl")
    include("test_hillslope_hydrology.jl")
    include("test_photosynthesis.jl")
    include("test_photosynthesis_ad_pft_params.jl")
    include("test_canopy_fluxes.jl")
    include("test_canopy_reverse.jl")
    include("test_canopy_hydrology.jl")
    include("test_bareground_fluxes.jl")
    include("test_soil_fluxes.jl")
    include("test_urban_fluxes.jl")
    include("test_luna.jl")
    include("test_ozone.jl")
    include("test_phenology.jl")
    include("test_allocation.jl")
    include("test_growth_resp.jl")
    include("test_maint_resp.jl")
    include("test_fun.jl")
    include("test_gap_mortality.jl")
    include("test_sparse_matrix_multiply.jl")
    include("test_matrixcn_gpu_ready.jl")
    include("test_cn_soil_matrix.jl")
    include("test_cn_soil_matrix_input.jl")
    include("test_cn_soil_matrix_advance.jl")
    include("test_cn_soil_matrix_bgc_topology.jl")
    include("test_cn_soil_matrix_iso.jl")
    include("test_dwt_matrix_input.jl")
    include("test_decomp_bgc.jl")
    include("test_decomp_mimics.jl")
    include("test_soil_biogeochem_decomp.jl")
    include("test_decomp_competition.jl")
    include("test_decomp_vertical_profile.jl")
    include("test_decomp_potential.jl")
    include("test_decomp_precision_control.jl")
    include("test_litter_vert_transp.jl")
    include("test_lvt_tri_ma_vr.jl")
    include("test_nitrif_denitrif.jl")
    include("test_n_leaching.jl")
    include("test_n_dynamics.jl")
    include("test_ndep_stream.jl")
    include("test_c_state_update1.jl")
    include("test_c_state_update2.jl")
    include("test_c_state_update3.jl")
    include("test_cn_veg_matrix.jl")
    include("test_cn_veg_matrix_wiring.jl")
    include("test_cn_veg_matrix_wiring_n.jl")
    include("test_cn_veg_matrix_wiring_crop.jl")
    include("test_cn_veg_matrix_wiring_crop_n.jl")
    include("test_matrix_cn_tails.jl")
    include("test_n_state_update1.jl")
    include("test_n_state_update2.jl")
    include("test_n_state_update3.jl")
    include("test_fire_base.jl")
    include("test_fire_li2014.jl")
    include("test_fire_methods.jl")
    include("test_cnfire_wiring.jl")
    include("test_methane.jl")
    include("test_voc_emission.jl")
    include("test_dust_emission.jl")
    include("test_satellite_phenology.jl")
    include("test_c_iso_flux.jl")
    include("test_cisoflux.jl")
    include("test_carbon_isotopes.jl")
    include("test_cn_balance_check.jl")
    include("test_balance_check.jl")
    include("test_irrigation.jl")
    include("test_cn_driver.jl")
    include("test_cn_precision_control.jl")
    include("test_veg_struct_update.jl")
    include("test_nutrient_competition.jl")
    include("test_flexiblecn.jl")
    include("test_flexiblecn_matrixcn_agsys.jl")
    include("test_vegetation_facade.jl")
    include("test_atm2lnd.jl")
    include("test_multigridcell_forcing.jl")
    include("test_lnd2atm.jl")
    include("test_topo.jl")
    include("test_subgrid_ave.jl")
    include("test_accumul.jl")
    include("test_accumul_restart.jl")
    include("test_crop_gdd.jl")
    include("test_pftcon.jl")
    include("test_control.jl")
    include("test_dyn_subgrid_control.jl")
    include("test_instances.jl")
    include("test_restart_io.jl")
    include("test_checkpoint.jl")
    include("test_init_interp.jl")
    include("test_band_diagonal.jl")
    include("test_clm_driver.jl")
    include("test_clm_wiring.jl")
    include("test_root_biophys.jl")
    include("test_soil_water_plant_sink.jl")
    include("test_sat_excess_runoff.jl")
    include("test_infilt_excess_runoff.jl")
    include("test_veg_compute_seed.jl")
    include("test_cn_annual_update.jl")
    include("test_cn_products_mod.jl")
    include("test_cndv.jl")
    include("test_cndv_wiring.jl")
    include("test_total_water_heat.jl")
    include("test_history_io.jl")
    include("test_multigridcell_history.jl")
    include("test_io_metadata.jl")
    include("test_dyn_file_io.jl")
    include("test_dyn_pft_crop_file.jl")
    include("test_dyn_harvest.jl")
    include("test_dyn_lake_urban_file.jl")
    include("test_dyn_landunit_area.jl")
    include("test_dyn_patch_state_updater.jl")
    include("test_dyn_cons_biogeophys.jl")
    include("test_dyn_cons_biogeochem.jl")
    include("test_dry_dep_velocity.jl")
    include("test_subgrid_build.jl")
    include("test_dyn_init_columns.jl")
    include("test_dyn_column_state_updater.jl")
    include("test_dyn_gross_unrep.jl")
    include("test_dyn_subgrid_driver.jl")
    include("test_dyn_subgrid_cnbal_wiring.jl")
    include("test_dyn_subgrid_e2e.jl")  # transient dyn_subgrid driven through the real clm_drv! loop
    include("test_initialization.jl")
    include("test_multigridcell_surfdata.jl")
    include("test_crop_collapse_static.jl")
    include("test_crop_lifecycle.jl")  # crop lifecycle vs the CTSM CN-crop Fortran reference
    include("test_multicolumn.jl")
    include("test_device_adapt.jl")
    include("test_backend.jl")  # centralized backend-selection API + multi-GPU-over-MPI plumbing (CPU proxy)
    include("test_validation_harness.jl")  # validation-matrix schema (runner exercised locally vs domain data)
    include("test_clm_run.jl")
    include("test_ad_smoke.jl")
    include("test_ad_e2e.jl")
    # test_ad_robustness.jl runs below in an isolated subprocess (flakes in-process
    # near AD phase-change discontinuities; passes standalone).
    include("test_calibration.jl")
    include("test_parameter_recovery.jl")
    include("test_enzyme_feasibility.jl")
    include("test_enzyme_smoke.jl")
    # Enzyme reverse-AD of the photosynthesis sub-phase (cf_rev_psn!) vs finite differences
    # — the psn counterpart of the NO-Enzyme forward guard test_canopy_reverse.jl.
    include("test_photosynthesis_reverse.jl")
    # Forward-parity guard for the productionized clm_drv! reverse phases
    # (src/driver/driver_reverse.jl). Builds a real inst → placed with the other
    # clm_initialize!-using integration tests, not the fresh-global unit block.
    include("test_driver_reverse.jl")
    include("test_multisite_calibration.jl")
    include("test_cn_integration.jl")

    # FATES (Tier F) — foundation-layer unit tests.
    include("test_fates_foundation.jl")
    include("test_fates_edparams.jl")
    include("test_fates_prtparams_fireweather.jl")
    include("test_fates_prtgeneric.jl")
    # FATES (Tier F) — PARTEH loss fluxes (turnover / drop / burn / damage / flush).
    include("test_fates_prtlossfluxes.jl")
    include("test_fates_litter_radmem.jl")
    include("test_fates_interfacetypes.jl")
    include("test_fates_sizeageindices.jl")
    # FATES (Tier F) — plant-hydraulics water transfer functions.
    include("test_fates_hydrowtf.jl")
    # FATES (Tier F) — plant-hydraulics memory/state types.
    include("test_fates_hydraulicsmem.jl")
    include("test_fates_twostream.jl")
    # FATES (Tier F) — fire: fuel-class enumeration + Nesterov fire-weather index.
    include("test_fates_fuelclasses_nesterov.jl")
    # FATES (Tier F) — PFT-indexed parameter container (EDPftvarcon, Batch 3).
    include("test_fates_edpftvarcon.jl")
    # FATES (Tier F) — fire: SPITFIRE parameter holder (SFParamsMod).
    include("test_fates_sfparams.jl")
    # FATES (Tier F) Batch 3 — SPITFIRE fuel state (FatesFuelMod).
    include("test_fates_fuel.jl")
    # FATES (Tier F) Batch 4 — statically-derived params + seed dispersal.
    include("test_fates_paramderived_dispersal.jl")
    # FATES (Tier F) Batch 5 — tree-damage module (DamageMainMod).
    include("test_fates_damage.jl")
    # FATES (Tier F) Batch 6 — allometry engine (FatesAllometryMod).
    include("test_fates_allometry.jl")
    # FATES (Tier F) Batch 7 — PARTEH carbon-only allometric allocation
    # hypothesis (PRTAllometricCarbonMod).
    include("test_fates_prtcarbon.jl")
    # FATES (Tier F) Batch 7 — CNP prioritized allometric allocation hypothesis.
    include("test_fates_prtcnp.jl")
    # FATES (Tier F) Batch 8 — the COHORT type (FatesCohortMod).
    include("test_fates_cohort.jl")
    # FATES (Tier F) Batch 9 — the PATCH type (FatesPatchMod).
    include("test_fates_patch.jl")
    # FATES (Tier F) Batch 10 — the SITE type (EDTypesMod, ed_site_type).
    include("test_fates_edtypes.jl")
    # FATES (Tier F) Batch 11 — PARTEH parameter init/registration (PRTParamsFATESMod).
    include("test_fates_prtparamsfates.jl")

    # FATES (Tier F) Batch 11 — ChecksBalancesMod (site mass-stock summation +
    # integrated flux/state balance check) + EDAccumulateFluxesMod
    # (daily per-cohort flux accumulation).
    include("test_fates_checksbalances_accfluxes.jl")
    # FATES (Tier F) Batch 11 — FATES <-> host soil-BGC flux coupling
    # (FatesSoilBGCFluxMod): litter/CWD fragmentation -> soil litter pools,
    # nutrient-acquisition BC pack/unpack, root-exudate efflux.
    include("test_fates_soilbgcflux.jl")
    # FATES (Tier F) Batch 11 — plant-hydraulics solver (FatesPlantHydraulicsMod).
    include("test_fates_planthydraulics.jl")
    # FATES plant-hydraulics Tier A — init/lifecycle routines.
    include("test_fates_planthydraulics_tierA.jl")
    # FATES plant-hydraulics Tier B — transpiration/uptake solve + gated driver callsite.
    include("test_fates_planthydraulics_tierB.jl")
    # FATES (Tier F) Batch 11 — SPITFIRE main fire driver (SFMainMod).
    include("test_fates_sfmain.jl")
    # FATES (Tier F) Batch 11 — two-stream radiation glue (FatesTwoStreamUtilsMod)
    # + land-use-change transitions (FatesLandUseChangeMod).
    include("test_fates_twostreamutils_landuse.jl")

    # FATES (Tier F) Batch 12 — water-stress btran + root uptake (EDBtranMod).
    include("test_fates_edbtran.jl")
    # FATES (Tier F) Batch 12 — cohort demographic engine (EDCohortDynamicsMod):
    # create/insert/sort/count/fuse/terminate of the cohort linked list.
    include("test_fates_edcohortdynamics.jl")
    # FATES (Tier F) Batch 12 — logging/harvest disturbance mortality
    # (EDLoggingMortalityMod): mortality fractions + killed-biomass litter routing.
    include("test_fates_edloggingmortality.jl")

    # FATES (Tier F) Batch 13 — salinity transpiration-stress (FatesBstressMod).
    include("test_fates_bstress.jl")
    # FATES (Tier F) Batch 13 — per-cohort mortality rates + Mortality_Derivative
    # (EDMortalityFunctionsMod): carbon-starvation/hydraulic/cold/background/damage.
    include("test_fates_edmortality.jl")
    # FATES (Tier F) Batch 13 — physiology hub (EDPhysiologyMod): phenology leaf
    # on/off, trim_canopy, recruitment/seed, CWD & litter input/fragmentation.
    include("test_fates_edphysiology.jl")
    # FATES (Tier F) Batch 14 — patch disturbance engine (EDPatchDynamicsMod):
    # disturbance_rates/spawn/split/fuse/terminate patches + *_litter_fluxes.
    include("test_fates_edpatchdynamics.jl")
    # FATES (Tier F) Batch 15 — canopy structure engine (EDCanopyStructureMod):
    # layer arrangement (demote/promote) + leaf-area profiles + host BC packing.
    include("test_fates_edcanopystructure.jl")
    # FATES (Tier F) Batch 16 — Norman (1979) two-stream canopy radiation
    # (FatesNormanRadMod), leaf photosynthesis+respiration (FatesPlantRespPhoto
    # synthMod), PSS/CSS inventory init (FatesInventoryInitMod), and the daily
    # ecosystem-dynamics driver (EDMainMod).
    include("test_fates_norman_rad.jl")
    include("test_fates_plantrespphotosynth.jl")
    include("test_fates_inventoryinit.jl")
    include("test_fates_edmain.jl")
    # FATES (Tier F) Batch 17 — canopy radiation driver (FatesRadiationDriveMod)
    # + FATES cold-start state init (EDInitMod) → init checkpoint.
    include("test_fates_radiationdrive.jl")
    include("test_fates_edinit.jl")
    # FATES (Tier F) Batch 18 — host coupling + I/O: history registry/aggregation
    # (FatesHistoryInterfaceMod), restart registry + demographic pack/unpack
    # round-trip (FatesRestartInterfaceMod), and the CLM↔FATES coupling seam
    # (FatesInterfaceMod).
    include("test_fates_history.jl")
    include("test_fates_restart.jl")
    include("test_fates_interface.jl")
    # Real FATES param-file reader: read_fates_params! parses the official FATES
    # default parameter file (data/fates/fates_params_default.cdl) and populates
    # every FATES param global; sample values match the CDL; cold-start finite.
    include("test_fates_params_reader.jl")
    # W1+W2 live-driver wiring spike: cold-start a single FATES site + attach to
    # CLMInstances + prove the AD-dual-copy / GPU-adapt ownership-model skip.
    include("test_fates_coldstart_spike.jl")

    # W3+W4 live-driver hooks: bc pack -> FATES driver -> bc unpack round-trip for
    # radiation (sun/shade + canopy albedo), btran, and photosynthesis, plus the
    # gated driver-branch column<->site mapping.
    include("test_fates_driver_hooks.jl")

    # FATES live-driver wiring W5 — the daily demographic step. Cold-starts a
    # carbon-only site, drives ed_ecosystem_dynamics + ed_update_site + the final
    # TotalBalanceCheck through fates_daily_dynamics_step! (the gated driver hook),
    # and asserts mass conservation + finite/physical cohort/patch/canopy state.
    include("test_fates_daily_dynamics.jl")

    # FATES history buffer-fills Wave 1 (B18-followup): cold-starts a use_fates
    # site (clm_fates_init! builds + Init's inst.fates.hist), drives update_history_dyn1!
    # via the daily hook + update_history_hifrq1! via fates_hifrq_history_step!, and
    # asserts the dyn1-core + hifrq1 history buffers populate finite, physical values.
    include("test_fates_history_dyn1.jl")

    # FATES history buffer-fills Wave 2: drives the level-2 (dim-level > 1) daily
    # update_history_dyn2! + per-timestep update_history_hifrq2! and asserts the
    # size/PFT/age-class DISAGGREGATED handles populate finite/physical values with
    # the cohort's mass landing in the correct size/pft/age class bin.
    include("test_fates_history_dyn2.jl")

    # FATES history buffer-fills Wave 3 (FINAL): constructs CNP-flux / fire (fuel +
    # fire-weather + fire-mort) / tree-damage / plant-hydraulics (si_hydr+co_hydr)
    # state on the cold-started site and asserts update_history_{nutrflux,dyn1,dyn2,
    # hydraulics}! fill the remaining nutrflux / fire-dimensioned / damage-cross-tab /
    # hydraulics handles with finite/physical values + correct size/pft/damage binning.
    include("test_fates_history_w3.jl")

    # FATES first REAL multi-day clm_drv!-integrated run: cold-starts a 14-PFT
    # carbon-only site on a FATES-tagged single column and loops the actual driver
    # for 2 days (96 steps), firing the four per-timestep FATES hooks + the daily
    # demographic step end-to-end; asserts finiteness, mass conservation, and that
    # FATES demonstrably drove the column (laisha/rssun populated by the FATES path).
    include("test_fates_driver_run.jl")

    # W4b in-solve FATES photosynthesis coupling: asserts FatesPlantRespPhotosynthDrive
    # now runs INSIDE canopy_fluxes_core!'s leaf-temperature Newton iteration (two-way
    # coupling) — the FATES bc_in t_veg_pa matches the converged in-loop leaf temp, and
    # rssun/rssha are finite/positive. Default (!use_fates) path stays byte-identical.
    include("test_fates_w4b_insolve.jl")

    # FATES live MODE exercises through the real clm_drv! loop: the SPITFIRE fire,
    # CNP nutrient-cycle, and tree-damage modes each cold-started + run >=2 days,
    # asserting finite + mass-conserving (C, and N/P for CNP) and that the mode
    # actually did something live (fire-weather/ignitions populated, CNP N/P demand
    # + 3 elements, damage class state active). Defines the shared
    # _fates_save_globals!/_fates_restore_globals/_fates_census helpers used by
    # test_fates_spinup.jl, so it MUST be included before it.
    include("test_fates_live_modes.jl")

    # FATES carbon-only longer-horizon spin-up / stability: a 15-day clm_drv! loop
    # asserting no blow-up/NaN, the census stays physical, carbon is conserved EVERY
    # day, and the demographic state evolves sanely (cohorts persist, patch census
    # ages/splits, carbon stays in a sane band). Reuses test_fates_live_modes.jl
    # helpers (included above).
    include("test_fates_spinup.jl")

    # FATES multi-veg-patch / multi-site coupling: the generalized ifp patch-walk
    # (p = ifp + col.patchi[c]) + setFilters-equivalent weight rebuild. PART A
    # splits a cold-start site into 2 vegetated patches and asserts the walk +
    # bc pack/unpack + fates_set_filters! cover both patch slots; PART B cold-starts
    # nsites=2 on two FATES columns and loops clm_drv! >= 1 day, asserting both run
    # finite + the daily step advances both sites.
    include("test_fates_multipatch.jl")

    # These tests each pass STANDALONE but flake when run in-process after the
    # full suite — a cumulative global-state effect (precompile / method-
    # invalidation / NCDataset-handle / float-state near AD discontinuities),
    # not a single resettable global (verified: no subset reproduces it). The
    # gated Fortran-parity files additionally depend on external dumps (skipped
    # cleanly when absent). Running each in a FRESH subprocess gives it a clean
    # global state. CLM precompilation is shared via the depot, so only the
    # first subprocess pays the precompile cost.
    @testset "Isolated subprocesses (global-state-sensitive)" begin
        isolated_files = [
            "test_fates_photosynthesis.jl",      # FATES daytime GPP guard (real Aripuana stand; gated on data)
            "test_fortran_parity.jl",            # per-step Julia↔Fortran (gated)
            "test_fortran_parity_cn.jl",         # CN/BGC/PHS/LUNA/FUN (gated)
            "test_fortran_parity_luna_decomp.jl",# LUNA cadence + decomp N-source (gated)
            "test_fortran_parity_freewins.jl",   # soilresis, aerosol free-wins (gated)
            "test_fortran_parity_luna_update.jl",# LUNA EOD vcmax/jmax update (gated)
            "test_multisite_robustness.jl",      # off-Bow driver robustness (gated)
            "test_glacier_robustness.jl",        # glacier istice driver robustness (gated)
            "test_lake_water_balance.jl",        # lake istdlak column water balance closes (gated)
            "test_lake_ref2m.jl",                # lake 2-m diagnostics vs Fortran + the four MO stability regimes (gated)
            "test_lake_eddy_ks.jl",              # lake eddy-diffusivity depth decay (ks) + the frozen roughness branch (gated)
            "test_snow_robustness.jl",           # deep-snowpack / cold-soil robustness (gated)
            "test_urban_robustness.jl",          # urban isturb roof/wall/road robustness (gated)
            "test_init_cold_wiring.jl",          # every *_init_cold! is CALLED on the live init path (+ NaN sweep, gated)
            "test_accum_wiring.jl",              # every accumulator is CALLED and actually AVERAGES (no hollow no-op bodies)
            "test_run_clm.jl",                   # standalone run harness: run→h0+restart→continue (gated)
            "test_ad_robustness.jl",             # AD near phase-change discontinuities
        ]
        for f in isolated_files
            path = joinpath(@__DIR__, f)
            cmd = `$(Base.julia_cmd()) --project=$(Base.active_project()) -e "using Test, CLM; include(raw\"$(path)\")"`
            # Under `Pkg.test` (CI), the parent sets JULIA_LOAD_PATH to the test
            # env only — no @stdlib — so a subprocess `using Printf` (stdlib, not a
            # test dep) fails with "Package Printf not found". Force the default
            # project+stdlib path so any stdlib `using` in the isolated file
            # resolves. @ resolves to --project; @stdlib provides Printf/etc.
            cmd = addenv(cmd, "JULIA_LOAD_PATH" => "@:@stdlib")
            @testset "$f (subprocess)" begin
                @test success(pipeline(cmd; stdout = stdout, stderr = stderr))
            end
        end
    end
end
