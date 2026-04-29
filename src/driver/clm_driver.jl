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
# ---------------------------------------------------------------------------
# BGC mode type hierarchy — enables dispatch-based specialization
# ---------------------------------------------------------------------------

"""Abstract supertype for biogeochemistry mode. Subtypes enable
the compiler to specialize and dead-code-eliminate unused branches."""
abstract type AbstractBGCMode end

"""Satellite Phenology mode — prescribed LAI, no biogeochemistry."""
struct SPMode <: AbstractBGCMode end

"""Carbon-Nitrogen biogeochemistry mode."""
Base.@kwdef struct CNMode <: AbstractBGCMode
    use_lch4::Bool = false
    use_c13::Bool = false
    use_c14::Bool = false
    use_crop::Bool = false
    use_cropcal_streams::Bool = false
    use_matrixcn::Bool = false
    ndep_from_cpl::Bool = false
    decomp_method::Int = 1
    no_soil_decomp::Int = 0
    ndecomp_pools::Int = 7
    ndecomp_cascade_transitions::Int = 10
    i_litr_min::Int = 1
    i_litr_max::Int = 3
    i_cwd::Int = 7
    npcropmin::Int = 17
    nrepr::Int = 1
end

"""FATES (Functionally Assembled Terrestrial Ecosystem Simulator) mode."""
Base.@kwdef struct FATESMode <: AbstractBGCMode
    use_fates_sp::Bool = false
    use_fates_bgc::Bool = false
    fates_spitfire_mode::Int = 0
    fates_seeddisp_cadence::Int = 0
end

# Convenience queries for dispatch
is_bgc_active(::SPMode) = false
is_bgc_active(::CNMode) = true
is_bgc_active(m::FATESMode) = m.use_fates_bgc

has_cn(::AbstractBGCMode) = false
has_cn(::CNMode) = true

has_fates(::AbstractBGCMode) = false
has_fates(::FATESMode) = true

has_lch4(::AbstractBGCMode) = false
has_lch4(m::CNMode) = m.use_lch4

"""Check if decomposition cascade has been properly initialized (non-zero donor pools)."""
_decomp_initialized(cascade::DecompCascadeConData) =
    !isempty(cascade.cascade_donor_pool) && any(x -> x != 0, cascade.cascade_donor_pool)

mutable struct CLMDriverConfig{M <: AbstractBGCMode}
    bgc_mode::M
    irrigate::Bool
    use_noio::Bool
    use_aquifer_layer::Bool
    use_soil_moisture_streams::Bool
    use_lai_streams::Bool
    n_drydep::Int
end

# Backward-compatible property accessors so existing code like config.use_cn still works
function Base.getproperty(c::CLMDriverConfig, s::Symbol)
    m = getfield(c, :bgc_mode)
    if s === :use_cn
        return m isa CNMode
    elseif s === :use_fates
        return m isa FATESMode
    elseif s === :use_fates_sp
        return m isa FATESMode && m.use_fates_sp
    elseif s === :use_fates_bgc
        return m isa FATESMode && m.use_fates_bgc
    elseif s === :use_lch4
        return m isa CNMode && m.use_lch4
    elseif s === :use_c13
        return m isa CNMode && m.use_c13
    elseif s === :use_c14
        return m isa CNMode && m.use_c14
    elseif s === :use_crop
        return m isa CNMode && m.use_crop
    elseif s === :use_cropcal_streams
        return m isa CNMode && m.use_cropcal_streams
    elseif s === :use_matrixcn
        return m isa CNMode && m.use_matrixcn
    elseif s === :ndep_from_cpl
        return m isa CNMode && m.ndep_from_cpl
    elseif s === :fates_spitfire_mode
        return m isa FATESMode ? m.fates_spitfire_mode : 0
    elseif s === :fates_seeddisp_cadence
        return m isa FATESMode ? m.fates_seeddisp_cadence : 0
    elseif s === :decomp_method
        return m isa CNMode ? m.decomp_method : 1
    elseif s === :no_soil_decomp
        return m isa CNMode ? m.no_soil_decomp : 0
    elseif s === :ndecomp_pools
        return m isa CNMode ? m.ndecomp_pools : 7
    elseif s === :ndecomp_cascade_transitions
        return m isa CNMode ? m.ndecomp_cascade_transitions : 10
    elseif s === :i_litr_min
        return m isa CNMode ? m.i_litr_min : 1
    elseif s === :i_litr_max
        return m isa CNMode ? m.i_litr_max : 3
    elseif s === :i_cwd
        return m isa CNMode ? m.i_cwd : 7
    elseif s === :npcropmin
        return m isa CNMode ? m.npcropmin : 17
    elseif s === :nrepr
        return m isa CNMode ? m.nrepr : 1
    else
        return getfield(c, s)
    end
end

# Backward-compatible keyword constructor: CLMDriverConfig(use_cn=true, ...)
function CLMDriverConfig(; use_cn::Bool=false, use_fates::Bool=false,
                          use_fates_sp::Bool=false, use_fates_bgc::Bool=false,
                          use_lch4::Bool=false, use_c13::Bool=false, use_c14::Bool=false,
                          use_crop::Bool=false, use_cropcal_streams::Bool=false,
                          use_matrixcn::Bool=false, ndep_from_cpl::Bool=false,
                          fates_spitfire_mode::Int=0, fates_seeddisp_cadence::Int=0,
                          irrigate::Bool=false, use_noio::Bool=false,
                          use_aquifer_layer::Bool=true,
                          use_soil_moisture_streams::Bool=false, use_lai_streams::Bool=false,
                          n_drydep::Int=0,
                          decomp_method::Int=1, no_soil_decomp::Int=0,
                          ndecomp_pools::Int=7, ndecomp_cascade_transitions::Int=10,
                          i_litr_min::Int=1, i_litr_max::Int=3, i_cwd::Int=7,
                          npcropmin::Int=17, nrepr::Int=1)
    if use_fates
        mode = FATESMode(; use_fates_sp, use_fates_bgc,
                          fates_spitfire_mode, fates_seeddisp_cadence)
    elseif use_cn
        mode = CNMode(; use_lch4, use_c13, use_c14, use_crop,
                       use_cropcal_streams, use_matrixcn, ndep_from_cpl,
                       decomp_method, no_soil_decomp,
                       ndecomp_pools, ndecomp_cascade_transitions,
                       i_litr_min, i_litr_max, i_cwd, npcropmin, nrepr)
    else
        mode = SPMode()
    end
    CLMDriverConfig(mode, irrigate, use_noio, use_aquifer_layer,
                    use_soil_moisture_streams, use_lai_streams, n_drydep)
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
# bitvec_to_filter — Convert BitVector mask to (count, indices) for Fortran-style APIs
# ---------------------------------------------------------------------------

"""
    bitvec_to_filter(mask::BitVector) -> (num::Int, filter::Vector{Int})

Convert a BitVector mask to the (count, index-array) form used by some
ported Fortran routines (e.g. urban_fluxes!) that have not yet been
converted to the BitVector interface.
"""
function bitvec_to_filter(mask::BitVector)
    indices = findall(mask)
    return (length(indices), indices)
end

# ---------------------------------------------------------------------------
# clm_drv! — Main CLM driver
# Ported from clm_drv in clm_driver.F90
# ---------------------------------------------------------------------------

"""
    clm_drv!(config, inst, filt, filt_inactive_and_active, bounds_proc,
             doalb, nextsw_cday, declinp1, declin, obliqr,
             rstwr, nlend, rdate, rof_prognostic;
             nstep, is_first_step, is_beg_curr_day, is_end_curr_day, is_beg_curr_year)

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
- `is_end_curr_day::Bool`      — true if end of current day
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
                   nextsw_cday::Real,
                   declinp1::Real,
                   declin::Real,
                   obliqr::Real,
                   rstwr::Bool,
                   nlend::Bool,
                   rdate::String,
                   rof_prognostic::Bool;
                   nstep::Int = 0,
                   is_first_step::Bool = false,
                   is_beg_curr_day::Bool = false,
                   is_end_curr_day::Bool = false,
                   is_beg_curr_year::Bool = false,
                   dtime::Real = 1800.0,
                   mon::Int = 1,
                   day::Int = 1,
                   photosyns::PhotosynthesisData = PhotosynthesisData())

    bounds_clump = bounds_proc  # single-clump mode

    # Shorthand aliases for frequently-used instances
    bc = bounds_clump
    col = inst.column
    lun = inst.landunit
    pch = inst.patch
    grc = inst.gridcell
    temp = inst.temperature
    ef = inst.energyflux
    cs = inst.canopystate
    ss = inst.soilstate
    sh = inst.soilhydrology
    fv = inst.frictionvel
    ls = inst.lakestate
    sa = inst.solarabs
    alb = inst.surfalb
    sr = inst.surfrad
    up = inst.urbanparams
    a2l = inst.atm2lnd
    wsb = inst.water.waterstatebulk_inst
    wdb = inst.water.waterdiagnosticbulk_inst
    wfb = inst.water.waterfluxbulk_inst
    oz = inst.ozone
    aer = inst.aerosol
    irr = inst.irrigation
    ps = photosyns

    # Bound ranges
    bc_col = bc.begc:bc.endc
    bc_patch = bc.begp:bc.endp
    bc_grc = bc.begg:bc.endg
    bc_lun = bc.begl:bc.endl

    # Compute column-level specific humidity from gridcell vapor pressure
    # q = 0.622 * e / (p - 0.378 * e), where e = forc_vp, p = forc_pbot
    nc = length(a2l.forc_pbot_downscaled_col)
    FT = eltype(a2l.forc_pbot_downscaled_col)
    forc_q_col = zeros(FT, nc)
    for c in bc_col
        g = col.gridcell[c]
        vp = a2l.forc_vp_grc[g]
        pbot = a2l.forc_pbot_downscaled_col[c]
        forc_q_col[c] = 0.622 * vp / max(pbot - 0.378 * vp, 1.0)
    end

    # ========================================================================
    # Glacier area initialization on first timestep
    # ========================================================================
    if is_first_step
        # Placeholder: update_glc2lnd_fracs!(bounds_clump) [glacier coupling]
        # Placeholder: dynSubgrid_wrapup_weight_changes!(bounds_clump) [dynamic subgrid]
    end

    # ========================================================================
    # Specified phenology
    # ========================================================================
    if config.use_cn
        if config.n_drydep > 0
            months, needs_read = interp_monthly_veg!(inst.satellite_phenology; kmo=mon, kda=day)
            if needs_read && inst.surfdata !== nothing
                read_monthly_vegetation!(inst.satellite_phenology, cs, pch, bc_patch;
                    monthly_lai=inst.surfdata.monthly_lai,
                    monthly_sai=inst.surfdata.monthly_sai,
                    monthly_height_top=inst.surfdata.monthly_htop,
                    monthly_height_bot=inst.surfdata.monthly_hbot,
                    months=months)
            end
        end
    elseif config.use_fates
        if config.use_fates_sp
            months, needs_read = interp_monthly_veg!(inst.satellite_phenology; kmo=mon, kda=day)
            if needs_read && inst.surfdata !== nothing
                read_monthly_vegetation!(inst.satellite_phenology, cs, pch, bc_patch;
                    monthly_lai=inst.surfdata.monthly_lai,
                    monthly_sai=inst.surfdata.monthly_sai,
                    monthly_height_top=inst.surfdata.monthly_htop,
                    monthly_height_bot=inst.surfdata.monthly_hbot,
                    months=months)
            end
        end
    else
        if doalb || config.n_drydep > 0
            months, needs_read = interp_monthly_veg!(inst.satellite_phenology; kmo=mon, kda=day)
            if needs_read && inst.surfdata !== nothing
                read_monthly_vegetation!(inst.satellite_phenology, cs, pch, bc_patch;
                    monthly_lai=inst.surfdata.monthly_lai,
                    monthly_sai=inst.surfdata.monthly_sai,
                    monthly_height_top=inst.surfdata.monthly_htop,
                    monthly_height_bot=inst.surfdata.monthly_hbot,
                    months=months)
            end
        end
    end

    # ========================================================================
    # Decomp vertical profiles
    # ========================================================================
    # ActiveLayer — WIRED
    alt_calc!(inst.active_layer, filt.soilc, temp, col, grc;
              mon=mon, day=day, sec=0, dtime=Int(dtime))
    if (config.use_cn || config.use_fates_bgc) && config.decomp_method != config.no_soil_decomp
        # SoilBiogeochemVerticalProfile — WIRED
        soil_bgc_vertical_profile!(
            inst.soilbiogeochem_state, inst.active_layer, ss, col, pch;
            mask_bgc_soilc=filt.bgc_soilc,
            mask_bgc_vegp=filt.bgc_vegp,
            nlevdecomp=varpar.nlevdecomp,
            nlevdecomp_full=varpar.nlevdecomp_full,
            dzsoi_decomp=dzsoi_decomp[],
            zsoi=zsoi[])
    end

    # ========================================================================
    # Initialize mass balance checks for C and N
    # ========================================================================
    if config.use_cn
        # InitEachTimeStep — WIRED
        cn_vegetation_init_each_timestep!(inst.bgc_vegetation;
            mask_soilc=filt.soilc,
            mask_soilp=filt.soilp,
            bounds_col=bc_col,
            bounds_patch=bc_patch)
    end
    if config.use_cn || config.use_fates_bgc
        # InitGridcellBalance — WIRED
        cn_vegetation_init_gridcell_balance!(inst.bgc_vegetation;
            mask_bgc_soilc=filt.bgc_soilc,
            mask_bgc_vegp=filt.bgc_vegp,
            mask_allc=filt.allc,
            bounds_col=bc_col,
            bounds_patch=bc_patch,
            soilbgc_cs=inst.soilbiogeochem_carbonstate,
            soilbgc_ns=inst.soilbiogeochem_nitrogenstate)
    end
    if config.use_lch4
        # Placeholder: ch4_init_gridcell_balance_check!(bc, ...) [CH4 gridcell balance]
    end

    # Begin water balance — WIRED
    water_gridcell_balance!(inst.water, ls, col, lun, grc,
                            filt.nolakec, filt.lakec,
                            bc_col, bc_lun, bc_grc, "begwb")

    # ========================================================================
    # Dynamic subgrid weights
    # ========================================================================
    # Placeholder: dynSubgrid_driver!(bounds_proc, ...) [dynamic vegetation area]

    # ========================================================================
    # Prescribed soil moisture / column balance init
    # ========================================================================
    if config.use_soil_moisture_streams
        # Placeholder: PrescribedSoilMoistureAdvance!(bounds_proc) [soil moisture streams]
        # Placeholder: PrescribedSoilMoistureInterp!(bc, ...) [soil moisture interp]
    end
    # BeginWaterColumnBalance — WIRED
    begin_water_column_balance!(inst.water, sh, ls, col, lun,
                                filt.nolakec, filt.lakec, bc_col)

    if config.use_cn || config.use_fates_bgc
        # InitColumnBalance — WIRED
        cn_vegetation_init_column_balance!(inst.bgc_vegetation;
            mask_bgc_soilc=filt.bgc_soilc,
            mask_bgc_vegp=filt.bgc_vegp,
            mask_allc=filt.allc,
            bounds_col=bc_col,
            bounds_patch=bc_patch,
            soilbgc_cs=inst.soilbiogeochem_carbonstate,
            soilbgc_ns=inst.soilbiogeochem_nitrogenstate)
    end
    if config.use_lch4
        # Placeholder: ch4_init_column_balance_check!(bc, ...) [CH4 column balance]
    end

    # ========================================================================
    # Update dynamic N deposition / forcing interpolation
    # ========================================================================
    if config.use_cn && !config.ndep_from_cpl
        # Placeholder: ndep_interp!(bounds_proc, a2l) [N deposition]
    end
    if config.use_cn
        # Placeholder: bgc_vegetation_inst%InterpFileInputs(bounds_proc) [CN file inputs]
    end
    # Placeholder: urbantv_interp!(bounds_proc) [urban TV]

    if doalb && config.use_lai_streams
        # Placeholder: lai_advance!(bounds_proc) [LAI streams]
    end
    if config.use_crop && config.use_cropcal_streams && is_beg_curr_year
        # Placeholder: cropcal_advance!(bounds_proc) [crop calendar]
    end

    # ========================================================================
    # DAYLENGTH
    # ========================================================================
    update_daylength!(grc, declin, obliqr, is_first_step, bc_grc)

    # ========================================================================
    # Driver initialization
    # ========================================================================
    clm_drv_init!(bc, cs, wsb, wdb, ef, ps,
                  col, pch, filt.nolakec, filt.nolakep, filt.soilp)

    # Placeholder: topo_inst%UpdateTopo!(bc, ...) [topography update]
    # Placeholder: downscale_forcings!(bc, ...) [atmospheric downscaling]

    # Update exposed vegetation filter
    set_exposedvegp_filter!(filt, bc,
                            cs.frac_veg_nosno_patch[bc_patch])

    # ========================================================================
    # Irrigation withdrawal
    # ========================================================================
    if config.irrigate
        # Placeholder: CalcAndWithdrawIrrigationFluxes!(bc, ...) [irrigation]
    end

    # ========================================================================
    # HYDROLOGY STAGE 1: Canopy interception, new snow, surface water
    # ========================================================================
    # Rain/snow/flood forcings (from atmosphere, filled by downscale_forcings! or forcing reader)
    np = length(pch.column)
    forc_rain_col = isdefined(Main, :nothing) ? zeros(FT, nc) : zeros(FT, nc)  # will use a2l fields
    forc_snow_col = zeros(FT, nc)
    qflx_floodg = zeros(FT, length(grc.lat))
    # Use Atm2LndData rain/snow if available (populated by downscale_forcings!/forcing_reader)
    if length(a2l.forc_rain_downscaled_col) == nc
        forc_rain_col = a2l.forc_rain_downscaled_col
    end
    if length(a2l.forc_snow_downscaled_col) == nc
        forc_snow_col = a2l.forc_snow_downscaled_col
    end

    # CanopyInterceptionAndThroughfall — WIRED
    canopy_interception_and_throughfall!(
        pch, col, cs, inst.water, dtime,
        filt.soilp, filt.nolakep, filt.nolakec,
        bc_patch, bc_col, bc_grc,
        forc_rain_col, forc_snow_col, a2l.forc_t_downscaled_col,
        a2l.forc_wind_grc,
        zeros(FT, np), zeros(FT, np))  # irrigation sprinkler/drip placeholders

    # HandleNewSnow — WIRED
    handle_new_snow!(temp, wsb, wdb, col, lun,
                     filt.nolakec, bc_col, dtime, varpar.nlevsno;
                     forc_t=a2l.forc_t_downscaled_col,
                     forc_wind=a2l.forc_wind_grc,
                     qflx_snow_grnd=wfb.wf.qflx_snow_grnd_col,
                     qflx_snow_drain=wfb.wf.qflx_snow_drain_col,
                     int_snow=wsb.int_snow_col,
                     scf_method=inst.scf_method)



    # UpdateFracH2oSfc — WIRED
    update_frac_h2osfc!(inst.water, col, filt.soilc, bc_col; dtime=dtime)

    # ========================================================================
    # SURFACE RADIATION
    # ========================================================================
    if !config.use_fates
        # CanopySunShadeFracs — WIRED
        canopy_sun_shade_fracs!(alb, cs, sa,
                                a2l.forc_solad_downscaled_col, a2l.forc_solai_grc,
                                pch, filt.nourbanp, bc_patch)
    else
        # Placeholder: clm_fates%wrap_sunfrac(...) [FATES sun/shade]
    end

    # SurfaceRadiation — WIRED
    surface_radiation!(alb, cs, sa, sr, wdb,
                       col, lun, grc, pch,
                       a2l.forc_solad_downscaled_col, a2l.forc_solai_grc,
                       filt.nourbanp, filt.urbanp, bc_patch;
                       dtime=dtime)

    # UrbanRadiation — WIRED
    urban_radiation!(sa, ef, col, lun, pch, up, temp, wdb,
                     a2l.forc_lwrad_downscaled_col,
                     a2l.forc_solad_downscaled_col, a2l.forc_solai_grc,
                     filt.nourbanl, filt.urbanl, filt.urbanc, filt.urbanp,
                     bc_lun, bc_patch)

    # ========================================================================
    # BIOGEOPHYSICS PRE-FLUX CALCULATIONS — WIRED
    # ========================================================================
    biogeophys_pre_flux_calcs!(cs, ef, fv, temp, ss, wsb, wdb, wfb,
                                col, lun, pch,
                                a2l.forc_t_downscaled_col, a2l.forc_th_downscaled_col, forc_q_col,
                                a2l.forc_hgt_u_grc, a2l.forc_hgt_t_grc, a2l.forc_hgt_q_grc,
                                filt.nolakec, filt.nolakep, filt.urbanc,
                                bc_col, bc_patch)

    # Soil evaporation resistance (soilbeta or soilresis)
    calc_soilevap_resis!(col, lun, ss, wsb, wdb, temp,
                         filt.nolakec, bc_col)

    # CalcOzoneStress — WIRED
    calc_ozone_stress!(oz, filt.exposedvegp, filt.noexposedvegp,
                       bc_patch, pch, pftcon.woody;
                       is_time_to_run_luna=false)

    # CalculateSurfaceHumidity — WIRED
    calculate_surface_humidity!(col, lun, temp, ss, wsb, wdb,
                                a2l.forc_pbot_downscaled_col, forc_q_col,
                                filt.nolakec, bc_col)

    # ========================================================================
    # DETERMINE FLUXES
    # ========================================================================

    # BareGroundFluxes — WIRED
    bareground_fluxes!(cs, ef, fv, temp, ss, wfb, wsb, wdb, ps,
                       pch, col, lun,
                       filt.noexposedvegp, bc_patch,
                       forc_q_col, a2l.forc_pbot_downscaled_col,
                       a2l.forc_th_downscaled_col, a2l.forc_rho_downscaled_col,
                       a2l.forc_t_downscaled_col,
                       a2l.forc_u_grc, a2l.forc_v_grc,
                       a2l.forc_hgt_t_grc, a2l.forc_hgt_u_grc, a2l.forc_hgt_q_grc)

    # Compute volumetric liquid water content for root moisture stress
    joff = varpar.nlevsno
    for j in 1:varpar.nlevgrnd
        for c in bc_col
            filt.nolakec[c] || continue
            dz_cj = col.dz[c, j + joff]
            if dz_cj > 0.0
                # Root stress expects volumetric liquid water bounded by the
                # effective pore space (matches SoilMoistStressMod semantics).
                liqvol = wsb.ws.h2osoi_liq_col[c, j + joff] / (dz_cj * DENH2O)
                wdb.h2osoi_liqvol_col[c, j + joff] =
                    min(max(liqvol, 0.0), ss.eff_porosity_col[c, j])
            else
                wdb.h2osoi_liqvol_col[c, j + joff] = 0.0
            end
        end
    end

    # SoilMoistStress (root moisture stress / BTRAN) — WIRED
    al = inst.active_layer
    calc_root_moist_stress!(ss, ef, temp, wsb, wdb, col, pch,
                             pftcon.smpso, pftcon.smpsc,
                             Float64.(al.altmax_lastyear_indx_col),
                             Float64.(al.altmax_indx_col),
                             filt.exposedvegp, bc_patch,
                             varpar.nlevgrnd, varpar.nlevsno)

    # CanopyFluxes — WIRED
    # Get downreg/leafn from CN vegetation facade when active, else zeros for SP mode
    if config.use_cn
        downreg_patch = get_downreg_patch(inst.bgc_vegetation, bc_patch)
        leafn_patch = get_leafn_patch(inst.bgc_vegetation, bc_patch)
    else
        downreg_patch = zeros(FT, np)
        leafn_patch = zeros(FT, np)
    end
    canopy_fluxes!(cs, ef, fv, temp, sa, ss, wfb, wsb, wdb, ps,
                   pch, col, grc,
                   filt.exposedvegp, bc_patch, bc_col,
                   a2l.forc_lwrad_downscaled_col,
                   forc_q_col, a2l.forc_pbot_downscaled_col,
                   a2l.forc_th_downscaled_col, a2l.forc_rho_downscaled_col,
                   a2l.forc_t_downscaled_col,
                   a2l.forc_u_grc, a2l.forc_v_grc,
                   a2l.forc_pco2_grc, a2l.forc_po2_grc,
                   a2l.forc_hgt_t_grc, a2l.forc_hgt_u_grc, a2l.forc_hgt_q_grc,
                   grc.dayl, grc.max_dayl,
                   downreg_patch, leafn_patch,
                   dtime;
                   t10_patch=temp.t_a10_patch,
                   nrad_patch=alb.nrad_patch,
                   tlai_z_patch=alb.tlai_z_patch,
                   vcmaxcint_sun_patch=alb.vcmaxcintsun_patch,
                   vcmaxcint_sha_patch=alb.vcmaxcintsha_patch,
                   parsun_z_patch=sa.parsun_z_patch,
                   parsha_z_patch=sa.parsha_z_patch,
                   laisun_z_patch=cs.laisun_z_patch,
                   laisha_z_patch=cs.laisha_z_patch,
                   o3coefv_patch=ones(FT, np),
                   o3coefg_patch=ones(FT, np),
                   dleaf_pft=pftcon.dleaf,
                   slatop_pft=pftcon.slatop,
                   leafcn_pft=pftcon.leafcn,
                   flnr_pft=pftcon.flnr,
                   fnitr_pft=pftcon.fnitr,
                   mbbopt_pft=pftcon.mbbopt,
                   c3psn_pft=pftcon.c3psn,
                   woody_pft=Float64.(pftcon.woody),
                   overrides=inst.overrides)

    # UrbanFluxes — WIRED (uses integer-filter API via bitvec_to_filter)
    (num_nourbanl, filter_nourbanl) = bitvec_to_filter(filt.nourbanl)
    (num_urbanl, filter_urbanl) = bitvec_to_filter(filt.urbanl)
    (num_urbanc, filter_urbanc) = bitvec_to_filter(filt.urbanc)
    (num_urbanp, filter_urbanp) = bitvec_to_filter(filt.urbanp)
    urban_fluxes!(ef, fv, temp, ss, up, wfb, wsb, wdb,
                  pch, col, lun,
                  num_nourbanl, filter_nourbanl,
                  num_urbanl, filter_urbanl,
                  num_urbanc, filter_urbanc,
                  num_urbanp, filter_urbanp,
                  bc_lun, bc_col, bc_patch,
                  a2l.forc_t_downscaled_col, a2l.forc_th_downscaled_col,
                  a2l.forc_rho_downscaled_col, forc_q_col,
                  a2l.forc_pbot_downscaled_col,
                  a2l.forc_u_grc, a2l.forc_v_grc;
                  dtime=dtime, nstep=nstep)

    # LakeFluxes — WIRED
    lake_fluxes!(temp, ef, fv, sa, ls, wsb, wdb, wfb,
                 col, pch, lun,
                 a2l.forc_t_downscaled_col, a2l.forc_th_downscaled_col, forc_q_col,
                 a2l.forc_pbot_downscaled_col, a2l.forc_rho_downscaled_col, a2l.forc_lwrad_downscaled_col,
                 a2l.forc_u_grc, a2l.forc_v_grc,
                 a2l.forc_hgt_u_grc, a2l.forc_hgt_t_grc, a2l.forc_hgt_q_grc,
                 filt.lakec, filt.lakep,
                 bc_col, bc_patch;
                 dtime=dtime)

    # SetActualRoughnessLengths — WIRED
    set_actual_roughness_lengths!(fv,
                                  filt.exposedvegp, filt.noexposedvegp,
                                  filt.urbanp, filt.lakep,
                                  bc_patch,
                                  pch.column, pch.landunit,
                                  lun.z_0_town)

    # ========================================================================
    # Irrigation needed
    # ========================================================================
    if config.irrigate
        # Placeholder: irrigation_inst%CalcIrrigationNeeded!(bc, ...) [irrigation needed]
    end

    # ========================================================================
    # EMISSIONS
    # ========================================================================
    # Placeholder: DustEmission!(bc, ...) [dust emission — deferred: Zender2003 mobilization]

    # DustDryDep — WIRED
    dust_dry_dep!(inst.dust_emis, BitVector(pch.active), pch.column, bc_patch,
                  a2l.forc_pbot_downscaled_col, a2l.forc_rho_downscaled_col,
                  a2l.forc_t_downscaled_col,
                  fv.ram1_patch, fv.fv_patch)

    # VOCEmission — no-op (MEGAN compounds not initialized; empty compound list)
    # To enable: initialize MEGANCompound/MEGANMechComp arrays and call voc_emission!()

    # ========================================================================
    # DETERMINE TEMPERATURES
    # ========================================================================

    # LakeTemperature — WIRED
    grnd_ch4_cond = zeros(FT, nc)  # placeholder (CH4 module provides when active)
    lake_temperature!(col, pch, sa, ss, wsb, wdb, wfb, ef, temp, ls,
                      grnd_ch4_cond,
                      filt.lakec, filt.lakep,
                      bc_col, bc_patch, dtime)

    # SoilTemperature — WIRED
    # Placeholder for urban building temperature max
    nl = length(lun.itype)
    urbantv_t_building_max = fill(FT(323.15), nl)  # placeholder ~50°C cap
    soil_temperature!(col, lun, pch, temp, ef, ss, wsb, wdb, wfb, sa, cs, up,
                      urbantv_t_building_max,
                      a2l.forc_lwrad_downscaled_col,
                      filt.nolakec, filt.nolakep,
                      filt.urbanl, filt.urbanc,
                      bc_col, bc_lun, bc_patch,
                      dtime)



    # Placeholder: glacier_smb_inst%HandleIceMelt!(bc, ...) [glacier ice melt]

    # ========================================================================
    # SURFACE FLUXES for new ground temperature — WIRED
    # ========================================================================
    soil_fluxes!(ef, temp, cs, wsb, wdb, wfb, sa,
                 pch, col, lun,
                 filt.nolakec, filt.nolakep, filt.urbanp,
                 bc_col, bc_patch,
                 a2l.forc_lwrad_downscaled_col,
                 dtime)

    # ========================================================================
    # Patch to column averaging
    # ========================================================================
    clm_drv_patch2col!(bc, filt.allc, filt.nolakec,
                       ef, wfb, col, pch)

    # ========================================================================
    # HYDROLOGY: Water balance physics (HydrologyNoDrainage)
    # Ported calling sequence from HydrologyNoDrainage in
    # HydrologyNoDrainageMod.F90 — ~20 physics calls that move water
    # through snow, surface, soil, and water table.
    # ========================================================================

    nlevsno = varpar.nlevsno
    nlevsoi = varpar.nlevsoi

    # --- 1. Build snow/no-snow filter ---
    build_snow_filter!(filt.snowc, filt.nosnowc, col.snl,
                       filt.nolakec, bc_col)
    # --- 2. SnowWater: snow percolation → qflx_rain_plus_snomelt ---
    # 2a. Apply top-layer fluxes (sublimation, dew, rain on snow)
    update_state_top_layer_fluxes!(
        wsb.ws.h2osoi_ice_col, wsb.ws.h2osoi_liq_col,
        dtime, col.snl, wdb.frac_sno_eff_col,
        wfb.wf.qflx_soliddew_to_top_layer_col,
        wfb.wf.qflx_solidevap_from_top_layer_col,
        wfb.wf.qflx_liq_grnd_col,
        wfb.wf.qflx_liqdew_to_top_layer_col,
        wfb.wf.qflx_liqevap_from_top_layer_col,
        filt.snowc, bc_col, nlevsno)

    # 2b. Liquid percolation through snow
    bulk_flux_snow_percolation!(
        wfb.wf.qflx_snow_percolation_col,
        dtime, col.snl, col.dz, wdb.frac_sno_eff_col,
        wsb.ws.h2osoi_ice_col, wsb.ws.h2osoi_liq_col,
        filt.snowc, bc_col, nlevsno)



    # 2c. Update snow layer liquid after percolation
    update_state_snow_percolation!(
        wsb.ws.h2osoi_liq_col,
        dtime, col.snl,
        wfb.wf.qflx_snow_percolation_col,
        filt.snowc, bc_col, nlevsno)



    # 2d. Aerosol fluxes through snow layers
    calc_and_apply_aerosol_fluxes!(
        aer, dtime, col.snl,
        wsb.ws.h2osoi_ice_col, wsb.ws.h2osoi_liq_col,
        wfb.wf.qflx_snow_percolation_col,
        filt.snowc, bc_col, nlevsno)



    # 2e. Adjust layer thicknesses after percolation
    post_percolation_adjust_layer_thicknesses!(
        col.dz, col.snl,
        wsb.ws.h2osoi_ice_col, wsb.ws.h2osoi_liq_col,
        filt.snowc, bc_col, nlevsno)



    # 2f. Sum percolation out of bottom snow layer → qflx_rain_plus_snomelt
    # For snow columns: drainage from bottom + non-snow-covered rain
    # For no-snow columns: rain + snowmelt
    # Extract bottom-layer percolation (layer 0 in Fortran = index nlevsno in Julia)
    qflx_snow_perc_bottom = zeros(FT, nc)
    for c in bc_col
        if filt.snowc[c]
            qflx_snow_perc_bottom[c] = wfb.wf.qflx_snow_percolation_col[c, nlevsno]
        end
    end
    sum_flux_add_snow_percolation!(
        wfb.wf.qflx_snow_drain_col,
        wfb.wf.qflx_rain_plus_snomelt_col,
        wdb.frac_sno_eff_col,
        qflx_snow_perc_bottom,
        wfb.wf.qflx_liq_grnd_col,
        wfb.wf.qflx_snomelt_col,
        filt.snowc, filt.nosnowc, bc_col)

    # --- 3. Soil water fractions (eff_porosity, icefrac) ---
    set_soil_water_fractions!(
        sh, ss, wsb.ws,
        col.dz, filt.hydrologyc, bc_col, nlevsoi, nlevsno)

    # --- 3b. Set flood water flux (gridcell → column) ---
    set_floodc!(
        wfb.wf.qflx_floodc_col, qflx_floodg,
        col.gridcell, col.itype,
        filt.nolakec, bc_col)

    # --- 4. Saturated excess runoff ---
    saturated_excess_runoff!(
        inst.sat_excess_runoff,
        filt.hydrologyc, bc_col,
        col, lun, sh, ss, wfb)

    # --- 5. Set qflx inputs for infiltration calculation ---
    set_qflx_inputs!(wfb, wdb, col.snl, filt.hydrologyc, bc_col)

    # --- 6. Infiltration excess runoff ---
    infiltration_excess_runoff!(
        inst.infilt_excess_runoff,
        sh, ss,
        inst.sat_excess_runoff.fsat_col,
        wfb, wdb,
        filt.hydrologyc, bc_col)

    # --- 7. Route infiltration excess ---
    route_infiltration_excess!(
        wfb, sh, col.landunit, lun.itype,
        filt.hydrologyc, bc_col)

    # --- 8. Update surface water (h2osfc) ---
    update_h2osfc!(col, sh, ef, wfb, wsb, wdb,
                   filt.hydrologyc, bc_col; dtime=dtime)

    # --- 9. Infiltration ---
    infiltration!(wfb, filt.hydrologyc, bc_col)

    # --- 10. Total surface runoff ---
    total_surface_runoff!(
        wfb, sh, wsb.ws,
        col.snl, col.itype, col.landunit, lun.urbpoi,
        filt.hydrologyc, filt.urbanc,
        bc_col, dtime)

    # --- 11. Root water uptake (transpiration sink) ---
    compute_effec_rootfrac_and_vert_tran_sink!(
        bc_col, nlevsoi,
        filt.hydrologyc,
        ss, cs, wfb, ef,
        col, lun, pch)

    # --- 12. Soil water movement (Richards equation) ---
    # Configure solver to match Fortran namelist:
    #   use_aquifer_layer=true  → Zeng-Decker 2009 with BC_AQUIFER (default)
    #   use_aquifer_layer=false → Moisture form with BC_ZERO_FLUX
    swm_cfg = if config.use_aquifer_layer
        SoilWaterMovementConfig()  # ZD09 + BC_AQUIFER
    else
        SoilWaterMovementConfig(soilwater_movement_method=MOISTURE_FORM,
                                lower_boundary_condition=BC_ZERO_FLUX)
    end
    soil_water!(col,
                filt.hydrologyc, filt.urbanc,
                sh, ss, wfb, wsb,
                temp, cs, ef,
                SoilWaterRetentionCurveClappHornberg1978(), swm_cfg,
                dtime)

    # --- 13. Water table ---
    water_table!(sh, ss,
                 temp.t_soisno_col, wsb.ws, wfb,
                 col.dz, col.z, col.zi,
                 filt.hydrologyc, bc_col,
                 nlevsoi, dtime)

    # Bedrock clipping: Fortran ThetaBasedWaterTable clips ZWT at bedrock depth
    # when soil above bedrock is not saturated. Without this, ZWT drops to 80m
    # in cold-start runs because the aquifer drains without replenishment.
    nlevsno_l = varpar.nlevsno
    for c in bc_col
        filt.hydrologyc[c] || continue
        nbr = col.nbedrock[c]
        zi_bedrock = col.zi[c, nbr + nlevsno_l + 1]
        if sh.zwt_col[c] > zi_bedrock
            sh.zwt_col[c] = zi_bedrock
        end
    end

    # --- 14. Condensation renewal ---
    renew_condensation!(
        wsb.ws, wdb, wfb,
        col.snl, col.itype,
        filt.hydrologyc, filt.urbanc,
        bc_col, dtime)

    # --- 15. Snow capping excess ---
    h2osno_total = zeros(FT, nc)
    waterstate_calculate_total_h2osno!(wsb.ws, filt.snowc, bc_col,
                                       col.snl, h2osno_total)
    # Init snow capping fluxes
    init_flux_snow_capping!(
        wfb.wf.qflx_snwcp_ice_col, wfb.wf.qflx_snwcp_liq_col,
        wfb.wf.qflx_snwcp_discarded_ice_col, wfb.wf.qflx_snwcp_discarded_liq_col,
        filt.snowc, bc_col)
    # Extract bottom snow layer copies (layer 0 in Fortran = nlevsno in Julia)
    jj_bottom = nlevsno
    dz_bottom_vec = col.dz[:, jj_bottom]
    ice_bottom_vec = wsb.ws.h2osoi_ice_col[:, jj_bottom]
    liq_bottom_vec = wsb.ws.h2osoi_liq_col[:, jj_bottom]
    mask_capping = falses(nc)
    rho_orig_bottom = zeros(FT, nc)
    frac_adjust = zeros(FT, nc)
    # Compute capping fluxes
    bulk_flux_snow_capping_fluxes!(
        mask_capping, rho_orig_bottom, frac_adjust,
        wfb.wf.qflx_snwcp_ice_col, wfb.wf.qflx_snwcp_liq_col,
        wfb.wf.qflx_snwcp_discarded_ice_col, wfb.wf.qflx_snwcp_discarded_liq_col,
        dtime, dz_bottom_vec, inst.topo.topo_col, h2osno_total,
        ice_bottom_vec, liq_bottom_vec,
        col.landunit, lun.itype,
        filt.snowc, bc_col, nlevsno, nstep)
    # Remove capping mass from bottom layer
    update_state_remove_snow_capping_fluxes!(
        ice_bottom_vec, liq_bottom_vec,
        dtime,
        wfb.wf.qflx_snwcp_ice_col, wfb.wf.qflx_snwcp_liq_col,
        wfb.wf.qflx_snwcp_discarded_ice_col, wfb.wf.qflx_snwcp_discarded_liq_col,
        mask_capping, bc_col)
    # Write back bottom layer state
    for c in bc_col
        if mask_capping[c]
            wsb.ws.h2osoi_ice_col[c, jj_bottom] = ice_bottom_vec[c]
            wsb.ws.h2osoi_liq_col[c, jj_bottom] = liq_bottom_vec[c]
        end
    end
    # Update dz and aerosols after capping
    snow_capping_update_dz_and_aerosols!(
        dz_bottom_vec, aer, jj_bottom,
        rho_orig_bottom, ice_bottom_vec, frac_adjust,
        mask_capping, bc_col)
    for c in bc_col
        if mask_capping[c]
            col.dz[c, jj_bottom] = dz_bottom_vec[c]
        end
    end

    # --- 16. Snow compaction ---
    snow_compaction!(
        col.dz, dtime, col.snl,
        temp.t_soisno_col,
        wsb.ws.h2osoi_ice_col, wsb.ws.h2osoi_liq_col,
        temp.imelt_col,
        wdb.frac_sno_col, wdb.frac_h2osfc_col,
        wdb.swe_old_col, wsb.int_snow_col, wdb.frac_iceold_col,
        a2l.forc_wind_grc, col.gridcell, col.landunit,
        lun.lakpoi, lun.urbpoi,
        filt.snowc, bc_col, nlevsno)



    # --- 17. Combine thin snow layers ---
    combine_snow_layers!(
        col.snl, col.dz, col.zi, col.z,
        temp.t_soisno_col,
        wsb.ws.h2osoi_ice_col, wsb.ws.h2osoi_liq_col,
        wsb.ws.h2osno_no_layers_col,
        wdb.snow_depth_col, wdb.frac_sno_col, wdb.frac_sno_eff_col,
        wsb.int_snow_col, wdb.snw_rds_col,
        aer, lun.itype, lun.urbpoi, col.landunit,
        filt.snowc, bc_col, nlevsno)
    # --- 18. Divide thick snow layers ---
    divide_snow_layers!(
        col.snl, col.dz, col.zi, col.z,
        temp.t_soisno_col,
        wsb.ws.h2osoi_ice_col, wsb.ws.h2osoi_liq_col,
        wdb.frac_sno_col, wdb.snw_rds_col,
        aer, false,  # is_lake=false
        filt.snowc, bc_col, nlevsno)


    # --- 19. Zero empty snow layers ---
    zero_empty_snow_layers!(
        col.snl, col.dz, col.z, col.zi,
        temp.t_soisno_col,
        wsb.ws.h2osoi_ice_col, wsb.ws.h2osoi_liq_col,
        filt.snowc, bc_col, nlevsno)

    # --- 20. Rebuild snow/no-snow filter after layer changes ---
    build_snow_filter!(filt.snowc, filt.nosnowc, col.snl,
                       filt.nolakec, bc_col)

    # --- Diagnostics (hydrology_no_drainage!) ---
    hydrology_no_drainage!(temp, ss, wsb, wdb,
                           col, lun,
                           filt.nolakec, filt.hydrologyc, filt.urbanc,
                           filt.snowc, filt.nosnowc,
                           bc_col,
                           dtime,
                           nlevsno, nlevsoi,
                           varpar.nlevgrnd, varpar.nlevurb)

    # Placeholder: glacier_smb_inst%ComputeSurfaceMassBalance!(bc, ...) [glacier SMB]

    # AerosolMasses (non-lake) — WIRED
    aerosol_masses!(aer,
                    filt.snowc, filt.nosnowc, bc_col,
                    col.snl, wsb.ws.h2osoi_ice_col, wsb.ws.h2osoi_liq_col,
                    wdb.h2osno_top_col, wdb.snw_rds_col)

    # ========================================================================
    # Lake hydrology
    # ========================================================================
    # LakeHydrology — WIRED
    lake_hydrology!(temp, ef, ls, ss, wsb, wdb,
                    inst.water.waterbalancebulk_inst, wfb,
                    col, pch,
                    filt.lakec, filt.lakep,
                    forc_rain_col, forc_snow_col, qflx_floodg,
                    bc_col, bc_patch, dtime,
                    varpar.nlevsno, varpar.nlevsoi, varpar.nlevgrnd)

    # AerosolMasses (lake) — WIRED
    aerosol_masses!(aer,
                    filt.lakesnowc, filt.lakenosnowc, bc_col,
                    col.snl, wsb.ws.h2osoi_ice_col, wsb.ws.h2osoi_liq_col,
                    wdb.h2osno_top_col, wdb.snw_rds_col)

    # SnowAge_grain (lake) — WIRED
    snowage_grain!(col.snl, col.dz, wdb.frac_sno_col,
                   wsb.ws.h2osoi_liq_col, wsb.ws.h2osoi_ice_col,
                   temp.t_soisno_col, temp.t_grnd_col,
                   wfb.wf.qflx_snow_grnd_col, wfb.wf.qflx_snofrz_lyr_col,
                   wsb.ws.h2osno_no_layers_col, a2l.forc_t_downscaled_col,
                   wdb.snw_rds_col, wdb.snw_rds_top_col,
                   wdb.sno_liq_top_col, temp.snot_top_col, temp.dTdz_top_col,
                   varpar.nlevsno, dtime;
                   mask_snowc=filt.lakesnowc, mask_nosnowc=filt.lakenosnowc)

    # ========================================================================
    # Urban snow fraction
    # ========================================================================
    for c in bc_col
        l = col.landunit[c]
        if lun.urbpoi[l]
            snow_depth = wdb.snow_depth_col[c]
            if !isfinite(snow_depth)
                snow_depth = 0.0
            end
            wdb.frac_sno_col[c] = min(snow_depth / 0.05, 1.0)
        end
    end

    # ========================================================================
    # Snow aging
    # ========================================================================
    # SnowAge_grain (non-lake) — WIRED
    snowage_grain!(col.snl, col.dz, wdb.frac_sno_col,
                   wsb.ws.h2osoi_liq_col, wsb.ws.h2osoi_ice_col,
                   temp.t_soisno_col, temp.t_grnd_col,
                   wfb.wf.qflx_snow_grnd_col, wfb.wf.qflx_snofrz_lyr_col,
                   wsb.ws.h2osno_no_layers_col, a2l.forc_t_downscaled_col,
                   wdb.snw_rds_col, wdb.snw_rds_top_col,
                   wdb.sno_liq_top_col, temp.snot_top_col, temp.dTdz_top_col,
                   varpar.nlevsno, dtime;
                   mask_snowc=filt.snowc, mask_nosnowc=filt.nosnowc)

    # ========================================================================
    # Ecosystem dynamics
    # ========================================================================
    if config.use_cn || config.use_fates_bgc
        # EcosystemDynamicsPreDrainage — WIRED
        cn_vegetation_ecosystem_pre_drainage!(inst.bgc_vegetation;
            mask_bgc_soilc=filt.bgc_soilc,
            mask_bgc_vegp=filt.bgc_vegp,
            mask_pcropp=filt.pcropp,
            mask_soilnopcropp=filt.soilnopcropp,
            mask_exposedvegp=filt.exposedvegp,
            mask_noexposedvegp=filt.noexposedvegp,
            bounds_col=bc_col,
            bounds_patch=bc_patch,
            nlevdecomp=varpar.nlevdecomp,
            nlevdecomp_full=varpar.nlevdecomp_full,
            ndecomp_pools=config.ndecomp_pools,
            ndecomp_cascade_transitions=config.ndecomp_cascade_transitions,
            i_litr_min=config.i_litr_min,
            i_litr_max=config.i_litr_max,
            i_cwd=config.i_cwd,
            npcropmin=config.npcropmin,
            nrepr=config.nrepr,
            patch_column=pch.column,
            ivt=pch.itype,
            woody=pftcon.woody,
            harvdate=inst.crop.harvdate_patch,
            col_is_fates=col.is_fates,
            cascade_donor_pool=inst.decomp_cascade.cascade_donor_pool,
            cascade_receiver_pool=inst.decomp_cascade.cascade_receiver_pool,
            dt=dtime,
            soilbgc_cs=inst.soilbiogeochem_carbonstate,
            soilbgc_cf=inst.soilbiogeochem_carbonflux,
            soilbgc_ns=inst.soilbiogeochem_nitrogenstate,
            soilbgc_nf=inst.soilbiogeochem_nitrogenflux,
            soilbgc_state=inst.soilbiogeochem_state,
            # Decomposition infrastructure (only pass if cascade is initialized)
            cascade_con=(_decomp_initialized(inst.decomp_cascade) ? inst.decomp_cascade : nothing),
            decomp_bgc_state=(_decomp_initialized(inst.decomp_cascade) ? inst.decomp_bgc_state : nothing),
            decomp_bgc_params=(_decomp_initialized(inst.decomp_cascade) ? inst.decomp_bgc_params : nothing),
            cn_shared_params=(_decomp_initialized(inst.decomp_cascade) ? inst.cn_shared_params : nothing),
            decomp_params=(_decomp_initialized(inst.decomp_cascade) ? inst.decomp_params : nothing),
            competition_state=(_decomp_initialized(inst.decomp_cascade) ? inst.competition_state : nothing),
            competition_params=(_decomp_initialized(inst.decomp_cascade) ? inst.competition_params : nothing),
            litter_params=(_decomp_initialized(inst.decomp_cascade) ? inst.litter_params : nothing),
            t_soisno=(_decomp_initialized(inst.decomp_cascade) ? @view(temp.t_soisno_col[:, (varpar.nlevsno+1):end]) : nothing),
            soilpsi=(_decomp_initialized(inst.decomp_cascade) ? ss.soilpsi_col : nothing),
            col=col,
            grc=grc,
            active_layer=(_decomp_initialized(inst.decomp_cascade) ? inst.active_layer : nothing),
            dzsoi_decomp=dzsoi_decomp[],
            zsoi_vals=zsoi[],
            zisoi_vals=zisoi[],
            mask_actfirec=filt.actfirec,
            mask_actfirep=filt.actfirep)
    end

    # Satellite phenology (SP mode) — WIRED
    if !config.use_cn && !config.use_fates && doalb
        satellite_phenology!(inst.satellite_phenology, cs, wdb,
                             pch, filt.nolakep, bc_patch)
    end
    if config.use_fates_sp && doalb
        # Placeholder: SatellitePhenology!(bc, ...) [FATES-SP satellite phenology]
    end

    # Dry deposition velocity — WIRED (only if n_drydep > 0)
    if config.n_drydep > 0
        depvel_compute!(inst.drydep, filt.nolakep, bc_patch,
                        pch.gridcell, pch.column, pch.landunit, pch.itype,
                        fv.ram1_patch, fv.rb1_patch, fv.fv_patch,
                        cs.elai_patch,
                        a2l.forc_t_downscaled_col, a2l.forc_solar_downscaled_col,
                        wdb.frac_sno_col, grc.lat, 1)
    end

    if config.use_crop && config.use_cropcal_streams && is_beg_curr_year
        # Placeholder: cropcal_interp!(bc, ...) [crop calendar interp]
    end

    # ========================================================================
    # HYDROLOGY: Drainage — WIRED
    # ========================================================================
    hydrology_drainage!(temp, sh, ss, wsb, wdb,
                        inst.water.waterbalancebulk_inst, wfb,
                        col, lun,
                        filt.nolakec, filt.hydrologyc, filt.urbanc, filt.do_smb_c,
                        forc_rain_col, forc_snow_col, qflx_floodg,
                        bc_col,
                        dtime,
                        varpar.nlevsno, varpar.nlevsoi,
                        varpar.nlevgrnd, varpar.nlevurb;
                        use_aquifer_layer=config.use_aquifer_layer)

    # Bedrock clipping (post-drainage): prevent ZWT from exceeding bedrock
    for c in bc_col
        filt.hydrologyc[c] || continue
        nbr = col.nbedrock[c]
        zi_bedrock = col.zi[c, nbr + varpar.nlevsno + 1]
        if sh.zwt_col[c] > zi_bedrock
            sh.zwt_col[c] = zi_bedrock
        end
    end

    # ========================================================================
    # Ecosystem dynamics post drainage
    # ========================================================================
    if config.use_cn || config.use_fates_bgc
        # EcosystemDynamicsPostDrainage — WIRED
        cn_vegetation_ecosystem_post_drainage!(inst.bgc_vegetation;
            mask_bgc_soilc=filt.bgc_soilc,
            mask_bgc_vegp=filt.bgc_vegp,
            mask_allc=filt.allc,
            mask_actfirec=filt.actfirec,
            mask_actfirep=filt.actfirep,
            bounds_col=bc_col,
            bounds_patch=bc_patch,
            nlevdecomp=varpar.nlevdecomp,
            ndecomp_pools=config.ndecomp_pools,
            ndecomp_cascade_transitions=config.ndecomp_cascade_transitions,
            i_litr_min=config.i_litr_min,
            i_litr_max=config.i_litr_max,
            i_cwd=config.i_cwd,
            dt=dtime,
            doalb=doalb,
            soilbgc_cs=inst.soilbiogeochem_carbonstate,
            soilbgc_cf=inst.soilbiogeochem_carbonflux,
            soilbgc_ns=inst.soilbiogeochem_nitrogenstate,
            soilbgc_nf=inst.soilbiogeochem_nitrogenflux)
    end

    # ========================================================================
    # FATES dynamics
    # ========================================================================
    if config.use_fates
        # Placeholder: clm_fates%WrapUpdateFatesRmean(...) [FATES running mean]
        # Placeholder: clm_fates%wrap_update_hifrq_hist(...) [FATES history]
        if is_beg_curr_day
            # Placeholder: clm_fates%dynamics_driv(...) [FATES dynamics]
            # Placeholder: setFilters!(...) [FATES filter update]
        end
    end

    # ========================================================================
    # Water diagnostic summaries
    # ========================================================================
    h2osno_total_col = zeros(FT, nc)
    waterstate_calculate_total_h2osno!(wsb.ws, filt.allc, bc_col, col.snl, h2osno_total_col)

    # WaterSummary — WIRED
    water_summary!(inst.water, bc_col, bc_patch;
                   mask_soilp=filt.soilp,
                   mask_allc=filt.allc,
                   mask_nolakec=filt.nolakec,
                   h2osno_total_col=h2osno_total_col,
                   dz_col=col.dz,
                   zi_col=col.zi,
                   landunit_col=col.landunit,
                   urbpoi=BitVector(lun.urbpoi),
                   lun_itype=lun.itype)

    # ========================================================================
    # Carbon and nitrogen balance check
    # ========================================================================
    if config.use_cn || config.use_fates_bgc
        # CN BalanceCheck — WIRED
        cn_vegetation_balance_check!(inst.bgc_vegetation;
            mask_bgc_soilc=filt.bgc_soilc,
            bounds_col=bc_col,
            nstep_since_startup=nstep,
            soilbgc_cf=inst.soilbiogeochem_carbonflux,
            soilbgc_nf=inst.soilbiogeochem_nitrogenflux,
            soilbgc_cs=inst.soilbiogeochem_carbonstate,
            soilbgc_ns=inst.soilbiogeochem_nitrogenstate)
    end

    # ========================================================================
    # Methane fluxes
    # ========================================================================
    if config.use_lch4
        # Placeholder: ch4!(bc, ...) [methane fluxes]
    end

    # ========================================================================
    # ALBEDOS for next time step — WIRED
    # ========================================================================
    if doalb
        # SurfaceAlbedoConstants — initialized in clm_initialize!, fallback to empty
        alb_con = inst.surfalb_con
        # Simple coszen function: cos(solar zenith) from calendar day, lat, lon, declination
        # NOTE: grc.lat and grc.lon are already in radians (set in init_gridcells.jl)
        # Matches Fortran shr_orb_cosz: cosz = sin(lat)*sin(decl) - cos(lat)*cos(decl)*cos(jday*2π + lon)
        coszen_cday = (cday, lat, lon, decl) -> begin
            hour_angle = 2.0 * π * mod(cday, 1.0) + lon - π
            max(sin(lat) * sin(decl) + cos(lat) * cos(decl) * cos(hour_angle), 0.0)
        end

        surface_albedo!(alb, alb_con, grc, col, lun, pch, cs, temp, wsb, wdb, ls, aer,
                        filt_inactive_and_active.nourbanc,
                        filt_inactive_and_active.nourbanp,
                        nextsw_cday, declinp1,
                        bc_grc, bc_col, bc_patch,
                        pftcon.rhol, pftcon.rhos, pftcon.taul, pftcon.taus, pftcon.xl,
                        coszen_cday)

        # UrbanAlbedo — WIRED
        urban_albedo!(filt.urbanl, filt.urbanc, filt.urbanp,
                      lun, col, pch, wsb, wdb, up, sa, alb)
    end

    # ========================================================================
    # FATES seed dispersal
    # ========================================================================
    if config.use_fates
        # Placeholder: clm_fates%WrapGlobalSeedDispersal() [FATES seeds]
    end

    # ========================================================================
    # Land to atmosphere
    # ========================================================================
    # Placeholder: lnd2atm!(bounds_proc, ...) [land-atmosphere coupling]

    # ========================================================================
    # Land to GLC
    # ========================================================================
    # Placeholder: lnd2glc_inst%update_lnd2glc!(bc, ...) [land-glacier coupling]

    # ========================================================================
    # Energy and water balance check
    # ========================================================================
    # End water balance — WIRED
    water_gridcell_balance!(inst.water, ls, col, lun, grc,
                            filt.nolakec, filt.lakec,
                            bc_col, bc_lun, bc_grc, "endwb")

    # BalanceCheck — WIRED
    balance_check!(inst.balcheck, wfb, wsb,
                   inst.water.waterbalancebulk_inst, wdb, ef, sa, cs, alb,
                   col, lun, pch, grc,
                   filt.allc, bc_col, bc_patch, bc_grc,
                   nstep, 0, dtime;
                   forc_rain_col=forc_rain_col,
                   forc_snow_col=forc_snow_col,
                   forc_rain_grc=length(a2l.forc_rain_not_downscaled_grc) > 0 ? a2l.forc_rain_not_downscaled_grc : zeros(FT, length(grc.lat)),
                   forc_snow_grc=length(a2l.forc_snow_not_downscaled_grc) > 0 ? a2l.forc_snow_not_downscaled_grc : zeros(FT, length(grc.lat)),
                   forc_solad_col=a2l.forc_solad_downscaled_col,
                   forc_solai_grc=a2l.forc_solai_grc,
                   forc_lwrad_col=a2l.forc_lwrad_downscaled_col,
                   forc_flood_grc=zeros(FT, length(grc.lat)),
                   qflx_ice_runoff_col=zeros(FT, nc),
                   qflx_evap_tot_grc=zeros(FT, length(grc.lat)),
                   qflx_surf_grc=zeros(FT, length(grc.lat)),
                   qflx_qrgwl_grc=zeros(FT, length(grc.lat)),
                   qflx_drain_grc=zeros(FT, length(grc.lat)),
                   qflx_drain_perched_grc=zeros(FT, length(grc.lat)),
                   qflx_ice_runoff_grc=zeros(FT, length(grc.lat)),
                   qflx_sfc_irrig_grc=zeros(FT, length(grc.lat)),
                   qflx_streamflow_grc=zeros(FT, length(grc.lat)))

    # ========================================================================
    # Diagnostics
    # ========================================================================
    write_diagnostic(nstep)

    # ========================================================================
    # Update accumulators
    # ========================================================================
    if nstep > 0
        # Atm2Lnd accumulator update — WIRED
        atm2lnd_update_acc_vars!(a2l, bc_patch, pch.gridcell, pch.column)

        # Temperature accumulators — WIRED
        temperature_update_acc_vars!(temp, bc_col, bc_patch, lun, pch;
            nstep=nstep, dtime=Int(dtime))

        # Canopy state accumulators — WIRED
        canopystate_update_acc_vars!(cs, bc_patch; nstep=nstep)

        # Water accumulators — WIRED
        water_update_acc_vars!(inst.water, bc_col)

        # Energy flux accumulators — WIRED
        energyflux_update_acc_vars!(ef, bc_patch;
            end_cd=is_end_curr_day, dtime=Int(dtime))

        # Placeholder: bgc_vegetation_inst%UpdateAccVars!(bounds_proc, ...) [BGC accum]

        if config.use_crop
            # Placeholder: crop_inst%CropUpdateAccVars!(bounds_proc, ...) [crop accum]
        end
        if config.use_fates
            # Placeholder: clm_fates%UpdateAccVars!(bounds_proc) [FATES accum]
        end
    end

    # ========================================================================
    # History buffer
    # ========================================================================
    # Placeholder: hist_update_hbuf!(bounds_proc) [history buffer update]

    # ========================================================================
    # Dynamic vegetation (CNDV)
    # ========================================================================
    if config.use_cn
        # EndOfTimeStepVegDynamics — WIRED
        cn_vegetation_end_of_timestep!(inst.bgc_vegetation;
            bounds_patch=bc_patch,
            is_end_curr_year=false,
            is_first_step=is_first_step)
    end

    # ========================================================================
    # History / Restart output
    # ========================================================================
    if !config.use_noio
        # Placeholder: hist_htapes_wrapup!(rstwr, nlend, bounds_proc, ...) [history write]
        if config.use_cn
            # Placeholder: bgc_vegetation_inst%WriteHistory!(bounds_proc) [BGC history]
        end
        if rstwr
            # Placeholder: restFile_write!(bounds_proc, ...) [restart write]
        end
    end

    return nothing
end
