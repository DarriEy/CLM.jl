# ==========================================================================
# Ported from: src/main/clm_driver.F90
# Main CLM driver — physics calling sequence
#
# This module provides the main CLM driver physics calling sequence.
# Computations occur over "clumps" of gridcells (and associated subgrid
# scale entities). In the Julia port, OMP parallel loops are replaced
# by sequential loops (GPU parallelism will be added later).
#
# Public functions:
#   clm_drv!            — Main CLM driver
#
# Private functions:
#   clm_drv_init!       — Initialize variables from previous timestep
#   clm_drv_patch2col!  — Average patch-level fluxes to column level
#   write_diagnostic    — Write diagnostic information
# ==========================================================================

# ---------------------------------------------------------------------------
# CLMDriverConfig — control flags for the driver
# ---------------------------------------------------------------------------

"""
    CLMDriverConfig

Configuration flags for the CLM driver. Aggregates control variables
from `clm_varctl` and other modules that determine which code paths
are active in the driver loop.

Ported from module-level `use` statements in `clm_driver.F90`.
"""
Base.@kwdef mutable struct CLMDriverConfig
    use_fates::Bool = false
    use_fates_sp::Bool = false
    use_fates_bgc::Bool = false
    use_cn::Bool = false
    use_lch4::Bool = false
    use_c13::Bool = false
    use_c14::Bool = false
    use_crop::Bool = false
    irrigate::Bool = false
    use_noio::Bool = false
    use_soil_moisture_streams::Bool = false
    use_cropcal_streams::Bool = false
    use_lai_streams::Bool = false
    ndep_from_cpl::Bool = false
    use_matrixcn::Bool = false
    fates_spitfire_mode::Int = 0
    fates_seeddisp_cadence::Int = 0
    n_drydep::Int = 0
    decomp_method::Int = 1
    no_soil_decomp::Int = 0
end

# ---------------------------------------------------------------------------
# CLMDriverState — mutable state needed across the driver timestep
# ---------------------------------------------------------------------------

"""
    CLMDriverState

Mutable state for tracking driver-level variables across the timestep,
such as time step number and date components.

Ported from local variables in `clm_drv` subroutine.
"""
Base.@kwdef mutable struct CLMDriverState
    nstep::Int = 0
    yr::Int = 0
    mon::Int = 0
    day::Int = 0
    sec::Int = 0
end

# ---------------------------------------------------------------------------
# clm_drv_init! — Initialize variables from previous timestep
# Ported from clm_drv_init in clm_driver.F90
# ---------------------------------------------------------------------------

"""
    clm_drv_init!(bounds, canopystate, waterstatebulk, waterdiagnosticbulk,
                  energyflux, photosyns, col_data, pch_data,
                  mask_nolakec, mask_nolakep, mask_soilp)

Initialization of CLM driver variables needed from previous timestep.

- Initializes intracellular CO2 parameters for VOCEmission
- Resets bottom heat flux
- Sets frac_veg_nosno from albedo calculation
- Computes ice fraction of snow at previous timestep

Ported from `clm_drv_init` in `clm_driver.F90`.
"""
function clm_drv_init!(bounds::BoundsType,
                        canopystate::CanopyStateData,
                        waterstatebulk::WaterStateBulkData,
                        waterdiagnosticbulk::WaterDiagnosticBulkData,
                        energyflux::EnergyFluxData,
                        photosyns::PhotosynthesisData,
                        col_data::ColumnData,
                        pch_data::PatchData,
                        mask_nolakec::BitVector,
                        mask_nolakep::BitVector,
                        mask_soilp::BitVector)

    nlevsno_val = varpar.nlevsno

    # Initialize intracellular CO2 (Pa) parameters each timestep for VOCEmission
    for p in bounds.begp:bounds.endp
        for j in axes(photosyns.cisun_z_patch, 2)
            photosyns.cisun_z_patch[p, j] = -999.0
            photosyns.cisha_z_patch[p, j] = -999.0
        end
    end

    # Reset flux from beneath soil/ice column
    for c in bounds.begc:bounds.endc
        energyflux.eflx_bot_col[c] = 0.0
    end

    # Initialize fraction of vegetation not covered by snow
    for p in bounds.begp:bounds.endp
        if pch_data.active[p]
            canopystate.frac_veg_nosno_patch[p] = canopystate.frac_veg_nosno_alb_patch[p]
        else
            canopystate.frac_veg_nosno_patch[p] = 0
        end
    end

    # Initialize set of previous time-step variables
    # Ice fraction of snow at previous time step
    # Fortran indices: j = -nlevsno+1 : 0
    # Julia 1-based:   j_jl = j + nlevsno = 1 : nlevsno
    for j in (-nlevsno_val + 1):0
        j_jl = j + nlevsno_val  # convert to 1-based Julia index
        for c in bounds.begc:bounds.endc
            mask_nolakec[c] || continue
            if j >= col_data.snl[c] + 1
                h2o_liq = waterstatebulk.ws.h2osoi_liq_col[c, j_jl]
                h2o_ice = waterstatebulk.ws.h2osoi_ice_col[c, j_jl]
                total = h2o_liq + h2o_ice
                if total > 0.0
                    waterdiagnosticbulk.frac_iceold_col[c, j_jl] = h2o_ice / total
                end
            end
        end
    end

    return nothing
end

# ---------------------------------------------------------------------------
# clm_drv_patch2col! — Average patch-level fluxes to column level
# Ported from clm_drv_patch2col in clm_driver.F90
# ---------------------------------------------------------------------------

"""
    clm_drv_patch2col!(bounds, mask_allc, mask_nolakec,
                        energyflux, waterfluxbulk,
                        col_data, pch_data)

Average over all patches for variables defined over both soil and lake to
provide the column-level averages of flux variables defined at the patch level.

Note: lake points are excluded from many of the averages because the
corresponding fields are computed in LakeHydrology, which is called after
this routine.

Ported from `clm_drv_patch2col` in `clm_driver.F90`.
"""
function clm_drv_patch2col!(bounds::BoundsType,
                             mask_allc::BitVector,
                             mask_nolakec::BitVector,
                             energyflux::EnergyFluxData,
                             waterfluxbulk::WaterFluxBulkData,
                             col_data::ColumnData,
                             pch_data::PatchData)

    # Averaging for patch evaporative flux variables (non-lake columns)
    p2c_1d_filter!(waterfluxbulk.qflx_ev_snow_col,
                   waterfluxbulk.qflx_ev_snow_patch,
                   mask_nolakec, col_data, pch_data)

    p2c_1d_filter!(waterfluxbulk.qflx_ev_soil_col,
                   waterfluxbulk.qflx_ev_soil_patch,
                   mask_nolakec, col_data, pch_data)

    p2c_1d_filter!(waterfluxbulk.qflx_ev_h2osfc_col,
                   waterfluxbulk.qflx_ev_h2osfc_patch,
                   mask_nolakec, col_data, pch_data)

    # Averaging for patch water flux variables (non-lake columns)
    p2c_1d_filter!(waterfluxbulk.wf.qflx_evap_soi_col,
                   waterfluxbulk.wf.qflx_evap_soi_patch,
                   mask_nolakec, col_data, pch_data)

    p2c_1d_filter!(waterfluxbulk.wf.qflx_evap_tot_col,
                   waterfluxbulk.wf.qflx_evap_tot_patch,
                   mask_nolakec, col_data, pch_data)

    p2c_1d_filter!(waterfluxbulk.wf.qflx_tran_veg_col,
                   waterfluxbulk.wf.qflx_tran_veg_patch,
                   mask_nolakec, col_data, pch_data)

    p2c_1d_filter!(waterfluxbulk.wf.qflx_liqevap_from_top_layer_col,
                   waterfluxbulk.wf.qflx_liqevap_from_top_layer_patch,
                   mask_nolakec, col_data, pch_data)

    # qflx_evap_soi averaged over all columns (including lake)
    p2c_1d_filter!(waterfluxbulk.wf.qflx_evap_soi_col,
                   waterfluxbulk.wf.qflx_evap_soi_patch,
                   mask_allc, col_data, pch_data)

    p2c_1d_filter!(waterfluxbulk.wf.qflx_liqdew_to_top_layer_col,
                   waterfluxbulk.wf.qflx_liqdew_to_top_layer_patch,
                   mask_nolakec, col_data, pch_data)

    p2c_1d_filter!(waterfluxbulk.wf.qflx_solidevap_from_top_layer_col,
                   waterfluxbulk.wf.qflx_solidevap_from_top_layer_patch,
                   mask_nolakec, col_data, pch_data)

    p2c_1d_filter!(waterfluxbulk.wf.qflx_soliddew_to_top_layer_col,
                   waterfluxbulk.wf.qflx_soliddew_to_top_layer_patch,
                   mask_nolakec, col_data, pch_data)

    return nothing
end

# ---------------------------------------------------------------------------
# write_diagnostic — Write diagnostic output
# Ported from write_diagnostic in clm_driver.F90
# ---------------------------------------------------------------------------

"""
    write_diagnostic(nstep; io=stdout)

Write diagnostic surface temperature output each timestep.

Ported from `write_diagnostic` in `clm_driver.F90`.
"""
function write_diagnostic(nstep::Int; io::IO=stdout)
    println(io, "clm: completed timestep ", nstep)
    return nothing
end

# ---------------------------------------------------------------------------
# clm_drv! — Main CLM driver
# Ported from clm_drv in clm_driver.F90
# ---------------------------------------------------------------------------

"""
    clm_drv!(config, inst, filt, filt_inactive_and_active, bounds_proc,
             doalb, nextsw_cday, declinp1, declin, obliqr,
             rstwr, nlend, rdate, rof_prognostic;
             nstep, is_first_step, is_beg_curr_day, is_beg_curr_year)

Main CLM driver — first phase of the CLM driver calling the CLM physics.

## Arguments
- `config::CLMDriverConfig`    — control flags
- `inst::CLMInstances`         — all CLM data instances
- `filt::ClumpFilter`          — active-only filter masks
- `filt_inactive_and_active::ClumpFilter` — inactive+active filter masks
- `bounds_proc::BoundsType`    — processor bounds
- `doalb::Bool`                — true if time for surface albedo calc
- `nextsw_cday::Float64`       — calendar day for nstep+1
- `declinp1::Float64`          — declination angle for next time step
- `declin::Float64`            — declination angle for current time step
- `obliqr::Float64`            — obliquity
- `rstwr::Bool`                — true => write restart file this step
- `nlend::Bool`                — true => end of run on this step
- `rdate::String`              — restart file time stamp for name
- `rof_prognostic::Bool`       — whether running with prognostic ROF component

## Keyword Arguments
- `nstep::Int`                 — current time step number
- `is_first_step::Bool`        — true if first time step
- `is_beg_curr_day::Bool`      — true if beginning of current day
- `is_beg_curr_year::Bool`     — true if beginning of current year
- `photosyns::PhotosynthesisData` — photosynthesis state (not in CLMInstances)

Ported from `clm_drv` in `clm_driver.F90`.

Note: This is a high-level orchestrator. Many called subroutines are
already ported (surface radiation, canopy fluxes, soil temperature, etc.).
In this port, each call is represented; if the called function is not yet
wired up, it is documented as a placeholder comment with the Fortran
subroutine name.
"""
function clm_drv!(config::CLMDriverConfig,
                   inst::CLMInstances,
                   filt::ClumpFilter,
                   filt_inactive_and_active::ClumpFilter,
                   bounds_proc::BoundsType,
                   doalb::Bool,
                   nextsw_cday::Float64,
                   declinp1::Float64,
                   declin::Float64,
                   obliqr::Float64,
                   rstwr::Bool,
                   nlend::Bool,
                   rdate::String,
                   rof_prognostic::Bool;
                   nstep::Int = 0,
                   is_first_step::Bool = false,
                   is_beg_curr_day::Bool = false,
                   is_beg_curr_year::Bool = false,
                   photosyns::PhotosynthesisData = PhotosynthesisData())

    bounds_clump = bounds_proc  # single-clump mode

    # ========================================================================
    # Glacier area initialization on first timestep
    # (In Fortran: glc2lnd_inst%update_glc2lnd_fracs + dynSubgrid_wrapup_weight_changes)
    # ========================================================================
    if is_first_step
        # Placeholder: update_glc2lnd_fracs!(bounds_clump)
        # Placeholder: dynSubgrid_wrapup_weight_changes!(bounds_clump)
    end

    # ========================================================================
    # Specified phenology
    # ========================================================================
    if config.use_cn
        if config.n_drydep > 0
            # Placeholder: interpMonthlyVeg!(bounds_proc, inst.canopystate)
        end
    elseif config.use_fates
        if config.use_fates_sp
            # Placeholder: interpMonthlyVeg!(bounds_proc, inst.canopystate)
        end
    else
        if doalb || config.n_drydep > 0
            # Placeholder: interpMonthlyVeg!(bounds_proc, inst.canopystate)
        end
    end

    # ========================================================================
    # Decomp vertical profiles
    # (alt_calc and SoilBiogeochemVerticalProfile)
    # ========================================================================
    # Placeholder: active_layer_inst%alt_calc(...)
    if (config.use_cn || config.use_fates_bgc) && config.decomp_method != config.no_soil_decomp
        # Placeholder: SoilBiogeochemVerticalProfile!(bounds_clump, ...)
    end

    # ========================================================================
    # Initialize mass balance checks for C and N
    # ========================================================================
    if config.use_cn
        # Placeholder: bgc_vegetation_inst%InitEachTimeStep(bounds_clump, ...)
    end

    # Gridcell-level balance check init
    if config.use_cn || config.use_fates_bgc
        # Placeholder: bgc_vegetation_inst%InitGridcellBalance(bounds_clump, ...)
    end
    if config.use_lch4
        # Placeholder: ch4_init_gridcell_balance_check!(bounds_clump, ...)
    end

    # Begin water balance
    # Placeholder: water_gridcell_balance!(bounds_clump, ..., flag="begwb")

    # ========================================================================
    # Dynamic subgrid weights
    # ========================================================================
    # Placeholder: dynSubgrid_driver!(bounds_proc, ...)

    # ========================================================================
    # Prescribed soil moisture from streams
    # ========================================================================
    if config.use_soil_moisture_streams
        # Placeholder: PrescribedSoilMoistureAdvance!(bounds_proc)
    end

    # ========================================================================
    # Column-level mass balance init
    # ========================================================================
    if config.use_soil_moisture_streams
        # Placeholder: PrescribedSoilMoistureInterp!(bounds_clump, ...)
    end
    # Placeholder: BeginWaterColumnBalance!(bounds_clump, ...)

    if config.use_cn || config.use_fates_bgc
        # Placeholder: bgc_vegetation_inst%InitColumnBalance(bounds_clump, ...)
    end
    if config.use_lch4
        # Placeholder: ch4_init_column_balance_check!(bounds_clump, ...)
    end

    # ========================================================================
    # Update dynamic N deposition
    # ========================================================================
    if config.use_cn && !config.ndep_from_cpl
        # Placeholder: ndep_interp!(bounds_proc, inst.atm2lnd)
    end

    if config.use_cn
        # Placeholder: bgc_vegetation_inst%InterpFileInputs(bounds_proc)
    end

    # Placeholder: urbantv_interp!(bounds_proc)

    # LAI streams
    if doalb && config.use_lai_streams
        # Placeholder: lai_advance!(bounds_proc)
    end

    # Crop calendar streams
    if config.use_crop && config.use_cropcal_streams && is_beg_curr_year
        # Placeholder: cropcal_advance!(bounds_proc)
    end

    # ========================================================================
    # Driver initialization, downscaling, canopy interception
    # ========================================================================

    # Update daylength
    # Placeholder: update_daylength!(bounds_clump, declin=declin, obliquity=obliqr)

    # Initialize variables needed for new driver time step
    clm_drv_init!(bounds_clump,
                  inst.canopystate,
                  inst.water.waterstatebulk_inst,
                  inst.water.waterdiagnosticbulk_inst,
                  inst.energyflux,
                  photosyns,
                  inst.column,
                  inst.patch,
                  filt.nolakec,
                  filt.nolakep,
                  filt.soilp)

    # Placeholder: topo_inst%UpdateTopo!(bounds_clump, ...)
    # Placeholder: downscale_forcings!(bounds_clump, ...)
    # Placeholder: set_atm2lnd_water_tracers!(bounds_clump, ...)

    # Update exposed vegetation filter
    set_exposedvegp_filter!(filt, bounds_clump,
                            inst.canopystate.frac_veg_nosno_patch[bounds_clump.begp:bounds_clump.endp])

    # ========================================================================
    # Irrigation withdrawal
    # ========================================================================
    if config.irrigate
        # Placeholder: CalcAndWithdrawIrrigationFluxes!(bounds_clump, ...)
    end

    # ========================================================================
    # First stage of hydrology
    # ========================================================================
    # Placeholder: CanopyInterceptionAndThroughfall!(bounds_clump, ...)
    # Placeholder: HandleNewSnow!(bounds_clump, ...)
    # Placeholder: UpdateFracH2oSfc!(bounds_clump, ...)

    # ========================================================================
    # Surface radiation
    # ========================================================================
    if config.use_fates
        # Placeholder: clm_fates%wrap_sunfrac(nc, atm2lnd_inst, canopystate_inst)
    else
        # Placeholder: CanopySunShadeFracs!(...)
    end
    # Placeholder: SurfaceRadiation!(bounds_clump, ...)
    # Placeholder: UrbanRadiation!(bounds_clump, ...)

    # ========================================================================
    # Biogeophysics pre-flux calculations
    # ========================================================================
    # Placeholder: BiogeophysPreFluxCalcs!(bounds_clump, ...)
    # Placeholder: ozone_inst%CalcOzoneStress!(bounds_clump, ...)
    # Placeholder: CalculateSurfaceHumidity!(bounds_clump, ...)

    # ========================================================================
    # Determine fluxes
    # ========================================================================
    # Placeholder: BareGroundFluxes!(bounds_clump, ...)
    # Placeholder: CanopyFluxes!(bounds_clump, ...)
    # Placeholder: UrbanFluxes!(bounds_clump, ...)
    # Placeholder: LakeFluxes!(bounds_clump, ...)
    # Placeholder: frictionvel_inst%SetActualRoughnessLengths!(bounds_clump, ...)

    # ========================================================================
    # Irrigation needed
    # ========================================================================
    if config.irrigate
        # Placeholder: irrigation_inst%CalcIrrigationNeeded!(bounds_clump, ...)
    end

    # ========================================================================
    # Dust and VOC emissions
    # ========================================================================
    # Placeholder: dust_emis_inst%DustEmission!(bounds_clump, ...)
    # Placeholder: dust_emis_inst%DustDryDep!(bounds_clump, ...)
    # Placeholder: VOCEmission!(bounds_clump, ...)

    # ========================================================================
    # Determine temperatures
    # ========================================================================
    # Placeholder: LakeTemperature!(bounds_clump, ...)
    # Placeholder: SoilTemperature!(bounds_clump, ...)
    # Placeholder: glacier_smb_inst%HandleIceMelt!(bounds_clump, ...)

    # ========================================================================
    # Update surface fluxes for new ground temperature
    # ========================================================================
    # Placeholder: SoilFluxes!(bounds_clump, ...)

    # ========================================================================
    # Patch to column averaging
    # ========================================================================
    clm_drv_patch2col!(bounds_clump,
                       filt.allc,
                       filt.nolakec,
                       inst.energyflux,
                       inst.water.waterfluxbulk_inst,
                       inst.column,
                       inst.patch)

    # ========================================================================
    # Vertical soil and surface hydrology
    # ========================================================================
    # Placeholder: HydrologyNoDrainage!(bounds_clump, ...)
    # Placeholder: glacier_smb_inst%ComputeSurfaceMassBalance!(bounds_clump, ...)
    # Placeholder: AerosolMasses!(bounds_clump, ...) [snow filter]

    # ========================================================================
    # Lake hydrology
    # ========================================================================
    # Placeholder: LakeHydrology!(bounds_clump, ...)
    # Placeholder: AerosolMasses!(bounds_clump, ...) [lake snow filter]
    # Placeholder: SnowAge_grain!(bounds_clump, ...) [lake snow filter]

    # ========================================================================
    # Urban snow fraction
    # ========================================================================
    for c in bounds_clump.begc:bounds_clump.endc
        l = inst.column.landunit[c]
        if inst.landunit.urbpoi[l]
            snow_depth = inst.water.waterdiagnosticbulk_inst.snow_depth_col[c]
            inst.water.waterdiagnosticbulk_inst.frac_sno_col[c] = min(snow_depth / 0.05, 1.0)
        end
    end

    # ========================================================================
    # Snow aging
    # ========================================================================
    # Placeholder: SnowAge_grain!(bounds_clump, ...) [non-lake snow filter]

    # ========================================================================
    # Ecosystem dynamics
    # ========================================================================
    if config.use_cn || config.use_fates_bgc
        # Placeholder: bgc_vegetation_inst%EcosystemDynamicsPreDrainage!(bounds_clump, ...)
    end

    # Satellite phenology (SP mode)
    if !config.use_cn && !config.use_fates && doalb
        # Placeholder: SatellitePhenology!(bounds_clump, ...)
    end
    if config.use_fates_sp && doalb
        # Placeholder: SatellitePhenology!(bounds_clump, ...) [FATES SP uses inactive+active filters]
    end

    # Dry deposition
    # Placeholder: depvel_compute!(bounds_clump, ...)

    # Crop calendar
    if config.use_crop && config.use_cropcal_streams && is_beg_curr_year
        # Placeholder: cropcal_interp!(bounds_clump, ...)
    end

    # ========================================================================
    # Hydrology with drainage
    # ========================================================================
    # Placeholder: HydrologyDrainage!(bounds_clump, ...)

    # ========================================================================
    # Ecosystem dynamics post drainage
    # ========================================================================
    if config.use_cn || config.use_fates_bgc
        # Placeholder: bgc_vegetation_inst%EcosystemDynamicsPostDrainage!(bounds_clump, ...)
    end

    # ========================================================================
    # FATES dynamics
    # ========================================================================
    if config.use_fates
        # Placeholder: clm_fates%WrapUpdateFatesRmean(nc, temperature_inst)
        # Placeholder: clm_fates%wrap_update_hifrq_hist(bounds_clump, ...)
        if is_beg_curr_day
            # Placeholder: clm_fates%dynamics_driv(nc, bounds_clump, ...)
            # Placeholder: setFilters!(bounds_clump, ...)
        end
    end

    # ========================================================================
    # Water diagnostic summaries
    # ========================================================================
    # Placeholder: water_inst%Summary!(bounds_clump, ...)

    # ========================================================================
    # Carbon and nitrogen balance check
    # ========================================================================
    if config.use_cn || config.use_fates_bgc
        # Placeholder: bgc_vegetation_inst%BalanceCheck!(bounds_clump, ...)
    end

    # ========================================================================
    # Methane fluxes
    # ========================================================================
    if config.use_lch4
        # Placeholder: ch4!(bounds_clump, ...)
    end

    # ========================================================================
    # Albedos for next time step
    # ========================================================================
    if doalb
        # Placeholder: SurfaceAlbedo!(bounds_clump, ...)
        # Placeholder: UrbanAlbedo!(bounds_clump, ...)
    end

    # ========================================================================
    # FATES seed dispersal
    # ========================================================================
    if config.use_fates
        # Placeholder: clm_fates%WrapGlobalSeedDispersal()
    end

    # ========================================================================
    # Land to atmosphere
    # ========================================================================
    # Placeholder: lnd2atm!(bounds_proc, ...)

    # ========================================================================
    # Land to GLC
    # ========================================================================
    # Placeholder: lnd2glc_inst%update_lnd2glc!(bounds_clump, ...)

    # ========================================================================
    # Energy and water balance check
    # ========================================================================
    # Placeholder: water_gridcell_balance!(bounds_clump, ..., flag="endwb")
    # Placeholder: BalanceCheck!(bounds_clump, ...)

    # ========================================================================
    # Diagnostics
    # ========================================================================
    write_diagnostic(nstep)

    # ========================================================================
    # Update accumulators
    # ========================================================================
    if nstep > 0
        # Placeholder: atm2lnd_inst%UpdateAccVars!(bounds_proc)
        # Placeholder: temperature_inst%UpdateAccVars!(bounds_proc, crop_inst)
        # Placeholder: canopystate_inst%UpdateAccVars!(bounds_proc)
        # Placeholder: water_inst%UpdateAccVars!(bounds_proc)
        # Placeholder: energyflux_inst%UpdateAccVars!(bounds_proc)
        # Placeholder: bgc_vegetation_inst%UpdateAccVars!(bounds_proc, ...)

        if config.use_crop
            # Placeholder: crop_inst%CropUpdateAccVars!(bounds_proc, ...)
        end

        if config.use_fates
            # Placeholder: clm_fates%UpdateAccVars!(bounds_proc)
        end
    end

    # ========================================================================
    # History buffer
    # ========================================================================
    # Placeholder: hist_update_hbuf!(bounds_proc)

    # ========================================================================
    # Dynamic vegetation (CNDV)
    # ========================================================================
    if config.use_cn
        # Placeholder: bgc_vegetation_inst%EndOfTimeStepVegDynamics!(bounds_clump, ...)
    end

    # ========================================================================
    # History / Restart output
    # ========================================================================
    if !config.use_noio
        # Placeholder: hist_htapes_wrapup!(rstwr, nlend, bounds_proc, ...)
        if config.use_cn
            # Placeholder: bgc_vegetation_inst%WriteHistory!(bounds_proc)
        end
        if rstwr
            # Placeholder: restFile_write!(bounds_proc, filer, ...)
        end
    end

    return nothing
end
