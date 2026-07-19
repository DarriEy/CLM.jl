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
    use_cndv::Bool = false   # dynamic global vegetation (CNDV / DGVM) coupling
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
    # CN fire method (CTSM namelist `fire_method`; use cnfire_method_symbol to map the
    # string). Default :nofire keeps the driver byte-identical to the historical
    # no-fire path. Set to :li2014/:li2016/:li2021/:li2024 (and initialize with
    # clm_initialize!(; cnfire_method=…, flnfm=…, fhdm=…)) to run the Li fire chain.
    cnfire_method::Symbol = :nofire
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
    use_hydrstress::Bool   # plant hydraulic stress (PHS) photosynthesis path
    use_luna::Bool         # LUNA photosynthetic-N acclimation (vcmax25 from vcmx25_z)
    use_voc::Bool          # MEGAN VOC emissions (gated, like shr_megan_mechcomps_n>0)
    use_ozone::Bool        # ozone uptake/stress on photosynthesis (CalcOzoneUptake → CalcOzoneStress)
    megan::MEGANConfig     # MEGAN static descriptors (factors table + compound mappings)
    # Crop growing-degree-day accumulation manager. Like `megan`, this lives on
    # the (non-dual-copied) config rather than on `CLMInstances` so the
    # ForwardDiff dual-copy of the instances never traverses its non-numeric
    # AccumField vectors (String/Matrix{Bool}/Matrix{Int}). `nothing` until a
    # use_crop run lazily registers the GDD fields on the first accumulator step.
    accum_mgr::Union{AccumManager, Nothing}
    # Dust-emission scheme selector (DustEmisFactory equivalent). Defaults to
    # inactive so the default driver path runs no dust mobilization (only the
    # historically-wired dry deposition) and is byte-identical.
    dust::DustEmisConfig
    # Transient land-use (dynamic subgrid) state. Like `megan`/`accum_mgr`, this
    # non-numeric driver-state object lives on the config (NOT on the dual-copied
    # CLMInstances). `nothing` unless clm_initialize! built a DynSubgridState for a
    # transient-landcover run; the driver's per-timestep dynSubgrid hooks no-op when
    # it is `nothing`, keeping the default (non-transient) path byte-identical.
    dyn_subgrid::Any
    # On-node intra-process multi-clump parallelism (CTSM's OpenMP-over-clumps
    # model). When `false` (default) the driver runs a single clump = the whole proc
    # domain, byte-identical to the historical single-clump path. When `true` AND the
    # decomposition has been split into `get_proc_clumps() > 1` disjoint clumps, the
    # *clump-safe* per-column physics is run threaded over clumps (each clump owns a
    # disjoint, contiguous begc:endc / begp:endp / begg:endg slice of the shared
    # proc-level state arrays — CLM has no halo, so the clumps need zero
    # communication). See the thread-safety notes on `clm_drv_core!`.
    use_threaded_clumps::Bool
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
    elseif s === :use_cndv
        return m isa CNMode && m.use_cndv
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
    elseif s === :cnfire_method
        return m isa CNMode ? m.cnfire_method : :nofire
    else
        return getfield(c, s)
    end
end

# Backward-compatible keyword constructor: CLMDriverConfig(use_cn=true, ...)
function CLMDriverConfig(; use_cn::Bool=false, use_fates::Bool=false,
                          use_fates_sp::Bool=false, use_fates_bgc::Bool=false,
                          use_lch4::Bool=false, use_c13::Bool=false, use_c14::Bool=false,
                          use_crop::Bool=false, use_cropcal_streams::Bool=false,
                          use_matrixcn::Bool=false, use_cndv::Bool=false,
                          ndep_from_cpl::Bool=false,
                          fates_spitfire_mode::Int=0, fates_seeddisp_cadence::Int=0,
                          irrigate::Bool=false, use_noio::Bool=false,
                          # CTSM-derived default: .false. for clm5_0 (see clm_initialize!).
                          use_aquifer_layer::Bool=false,
                          use_soil_moisture_streams::Bool=false, use_lai_streams::Bool=false,
                          n_drydep::Int=0, use_hydrstress::Bool=false,
                          # CTSM-conditional default (namelist_defaults_ctsm.xml:578-580):
                          # .false. under FATES, .true. otherwise. See clm_initialize!.
                          use_luna::Union{Bool,Nothing}=nothing,
                          use_voc::Bool=false, use_ozone::Bool=false,
                          megan::Union{MEGANConfig,Nothing}=nothing,
                          dust::Union{DustEmisConfig,Nothing}=nothing,
                          dyn_subgrid=nothing,
                          decomp_method::Int=1, no_soil_decomp::Int=0,
                          ndecomp_pools::Int=7, ndecomp_cascade_transitions::Int=10,
                          i_litr_min::Int=1, i_litr_max::Int=3, i_cwd::Int=7,
                          npcropmin::Int=17, nrepr::Int=1,
                          cnfire_method::Symbol=:nofire,
                          use_threaded_clumps::Bool=false)
    if use_fates
        mode = FATESMode(; use_fates_sp, use_fates_bgc,
                          fates_spitfire_mode, fates_seeddisp_cadence)
    elseif use_cn
        mode = CNMode(; use_lch4, use_c13, use_c14, use_crop,
                       use_cropcal_streams, use_matrixcn, use_cndv, ndep_from_cpl,
                       decomp_method, no_soil_decomp,
                       ndecomp_pools, ndecomp_cascade_transitions,
                       i_litr_min, i_litr_max, i_cwd, npcropmin, nrepr,
                       cnfire_method)
    else
        mode = SPMode()
    end
    if megan === nothing
        megan = MEGANConfig()
        megan_factors_init!(megan.megan_factors, 20)
    end
    if dust === nothing
        dust = DustEmisConfig()
    end
    # Resolve the conditional default: FATES runs get .false. (CTSM endrun's on
    # LUNA+FATES, controlMod.F90:505); everything else gets CLM5's .true.
    _use_luna = use_luna === nothing ? !use_fates : use_luna
    CLMDriverConfig(mode, irrigate, use_noio, use_aquifer_layer,
                    use_soil_moisture_streams, use_lai_streams, n_drydep, use_hydrstress, _use_luna,
                    use_voc, use_ozone, megan, nothing, dust, dyn_subgrid, use_threaded_clumps)
end

# ---------------------------------------------------------------------------
# clm_run_clump_physics! — on-node multi-clump parallel driver loop
# (CTSM's OpenMP-over-clumps `do nc = 1, nclumps` model)
# ---------------------------------------------------------------------------

"""
    clm_run_clump_physics!(phys!, clump_bounds; threaded=true)

Run a clump-safe per-column physics step `phys!(bounds_clump)` over a list of
disjoint clump bounds, optionally threaded over clumps with `Threads.@threads`.

This is the Julia analogue of CTSM's `do nc = 1, nclumps ... end do`
(OpenMP-over-clumps, NOT MPI): each `bounds_clump` owns a disjoint, contiguous
slice of gridcells/columns/patches of the shared proc-level state arrays, and
CLM has no halo, so the clumps are embarrassingly parallel with zero
communication. `phys!` must only read/write the (disjoint) state slice indexed by
its `bounds_clump` (`bounds_clump.begc:bounds_clump.endc`, etc.) — proc-wide
filter masks are absolute-indexed and therefore correct for any clump sub-range
with no slicing.

Determinism: because the clumps own disjoint slices, the result is independent of
clump count and of whether the loop is threaded — `phys!` applied over `nclumps`
clumps produces exactly the same state as `phys!` applied over a single
whole-proc clump (provided `phys!` itself contains no cross-clump reductions or
shared-global writes). The caller is responsible for that contract; the orchestrator
`clm_drv_core!` keeps non-clump-safe phases (full-proc FATES site walks, lazy
accumulator init, gridcell scratch) OUT of `phys!` and runs them serially.

`threaded=false` forces a deterministic serial walk (used to validate
single-clump byte-identity and as the fallback when `Threads.nthreads() == 1`).
"""
function clm_run_clump_physics!(phys!, clump_bounds::AbstractVector{BoundsType};
                                threaded::Bool = true)
    if threaded && Threads.nthreads() > 1 && length(clump_bounds) > 1
        Threads.@threads for nc in eachindex(clump_bounds)
            phys!(clump_bounds[nc])
        end
    else
        for nc in eachindex(clump_bounds)
            phys!(clump_bounds[nc])
        end
    end
    return nothing
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

# ---- clm_drv_init! per-element kernels (GPU/CPU backend-agnostic) ----

# Init intracellular CO2 to -999 over patches [begp..endp] × canopy-levels.
@kernel function _drvinit_ci_kernel!(cisun, cisha, begp::Int, endp::Int, nlevcan::Int)
    I = @index(Global, Cartesian)
    p = I[1]; j = I[2]
    @inbounds if begp <= p <= endp && j <= nlevcan
        cisun[p, j] = -oftype(cisun[p, j], 999)
        cisha[p, j] = -oftype(cisha[p, j], 999)
    end
end

# Reset flux from beneath soil/ice column to 0 over [begc..endc].
@kernel function _drvinit_eflxbot_kernel!(eflx_bot, begc::Int, endc::Int)
    c = @index(Global)
    @inbounds if begc <= c <= endc
        eflx_bot[c] = zero(eltype(eflx_bot))
    end
end

# frac_veg_nosno = active ? frac_veg_nosno_alb : 0  over [begp..endp].
@kernel function _drvinit_fracvegnosno_kernel!(fvns, @Const(fvns_alb), @Const(active),
        begp::Int, endp::Int)
    p = @index(Global)
    @inbounds if begp <= p <= endp
        fvns[p] = active[p] ? fvns_alb[p] : zero(eltype(fvns))
    end
end

# Ice fraction of snow at previous time step over [begc..endc] × snow layers.
# Fortran j = -nlevsno+1 : 0 ; Julia jj = j + nlevsno = 1 : nlevsno.
@kernel function _drvinit_fraciceold_kernel!(frac_iceold, @Const(mask), @Const(snl),
        @Const(h2osoi_liq), @Const(h2osoi_ice), begc::Int, endc::Int, nlevsno::Int)
    c = @index(Global)
    @inbounds if begc <= c <= endc && mask[c]
        for jj in 1:nlevsno
            j = jj - nlevsno            # Fortran index in [-nlevsno+1 .. 0]
            if j >= snl[c] + 1
                h2o_liq = h2osoi_liq[c, jj]
                h2o_ice = h2osoi_ice[c, jj]
                total = h2o_liq + h2o_ice
                if total > zero(total)
                    frac_iceold[c, jj] = h2o_ice / total
                end
            end
        end
    end
end

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
                        mask_nolakec::AbstractVector{Bool},
                        mask_nolakep::AbstractVector{Bool},
                        mask_soilp::AbstractVector{Bool})

    nlevsno_val = varpar.nlevsno

    # Initialize intracellular CO2 (Pa) parameters each timestep for VOCEmission
    nlevcan = size(photosyns.cisun_z_patch, 2)
    _launch!(_drvinit_ci_kernel!, photosyns.cisun_z_patch, photosyns.cisha_z_patch,
             bounds.begp, bounds.endp, nlevcan;
             ndrange = (size(photosyns.cisun_z_patch, 1), nlevcan))

    # Reset flux from beneath soil/ice column
    _launch!(_drvinit_eflxbot_kernel!, energyflux.eflx_bot_col, bounds.begc, bounds.endc)

    # Initialize fraction of vegetation not covered by snow
    _launch!(_drvinit_fracvegnosno_kernel!, canopystate.frac_veg_nosno_patch,
             canopystate.frac_veg_nosno_alb_patch, pch_data.active,
             bounds.begp, bounds.endp)

    # Initialize set of previous time-step variables
    # Ice fraction of snow at previous time step (per-column kernel with an
    # internal sequential snow-layer loop over the padded slice).
    _launch!(_drvinit_fraciceold_kernel!, waterdiagnosticbulk.frac_iceold_col,
             mask_nolakec, col_data.snl,
             waterstatebulk.ws.h2osoi_liq_col, waterstatebulk.ws.h2osoi_ice_col,
             bounds.begc, bounds.endc, nlevsno_val;
             ndrange = size(waterdiagnosticbulk.frac_iceold_col, 1))

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
                             mask_allc::AbstractVector{Bool},
                             mask_nolakec::AbstractVector{Bool},
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

    # Canopy evaporation: aggregate patch→column so the QVEGE diagnostic is populated
    # (the energy form FCEV=eflx_lh_vege was fine, but qflx_evap_veg_col stayed 0).
    p2c_1d_filter!(waterfluxbulk.wf.qflx_evap_veg_col,
                   waterfluxbulk.wf.qflx_evap_veg_patch,
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
function bitvec_to_filter(mask::AbstractVector{Bool})
    # findall on the host (Array() a device mask first) — this builds an integer
    # filter for the host-side urban integer-filter API.
    m = mask isa Array ? mask : Array(mask)
    indices = findall(m)
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
# Keyword-argument interface used throughout the codebase and tests. Thin
# wrapper that forwards to the fully-positional clm_drv_core! below.
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
                   is_end_curr_year::Bool = false,
                   dtime::Real = 1800.0,
                   mon::Int = 1,
                   day::Int = 1,
                   secs::Int = 0,
                   jday::Int = 1,
                   year::Int = 0,
                   photosyns::PhotosynthesisData = PhotosynthesisData())
    # FATES is CPU-only (host-fallback): its driver hooks scalar-index the per-column
    # CLM arrays, which are device arrays when the inst lives on a GPU. Under use_fates,
    # run the step inside GPUArraysCore.@allowscalar so those few per-column/patch scalar
    # accesses are permitted seamlessly. This is a no-op on a CPU inst (and on any CPU-only
    # run), and the non-FATES biogeophysics uses kernels (never scalar indexing), so it
    # hides nothing there — shared scalar-indexing bugs still surface on non-fates GPU runs
    # (which take the plain branch below). The default (!use_fates) path is byte-identical.
    # config.use_fates is a compile-time property of CLMDriverConfig{Mode}, so the branch
    # folds away and the differentiable !use_fates entry never sees the wrapper.
    if config.use_fates
        return GPUArraysCore.@allowscalar clm_drv_core!(config, inst, filt,
            filt_inactive_and_active, bounds_proc, doalb, nextsw_cday, declinp1, declin,
            obliqr, rstwr, nlend, rdate, rof_prognostic, nstep, is_first_step,
            is_beg_curr_day, is_end_curr_day, is_beg_curr_year, dtime, mon, day, photosyns,
            is_end_curr_year, secs, jday, year)
    end
    return clm_drv_core!(config, inst, filt, filt_inactive_and_active, bounds_proc,
        doalb, nextsw_cday, declinp1, declin, obliqr, rstwr, nlend, rdate, rof_prognostic,
        nstep, is_first_step, is_beg_curr_day, is_end_curr_day, is_beg_curr_year,
        dtime, mon, day, photosyns, is_end_curr_year, secs, jday, year)
end

# Fully-positional core: the differentiable entry point. Enzyme cannot build the
# augmented kwarg NamedTuple when it contains an active struct (photosyns), so
# reverse-mode differentiates this positional method instead.
function clm_drv_core!(config::CLMDriverConfig,
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
                   rof_prognostic::Bool,
                   nstep::Int,
                   is_first_step::Bool,
                   is_beg_curr_day::Bool,
                   is_end_curr_day::Bool,
                   is_beg_curr_year::Bool,
                   dtime::Real,
                   mon::Int,
                   day::Int,
                   photosyns::PhotosynthesisData,
                   is_end_curr_year::Bool = false,
                   secs::Int = 0,
                   jday::Int = 1,
                   year::Int = 0)

    # ------------------------------------------------------------------------
    # On-node multi-clump decomposition (CTSM's `do nc = 1, nclumps` model).
    #
    # CTSM runs the per-timestep physics inside a `do nc = 1, nclumps` loop,
    # OpenMP-threaded over CLUMPS (NOT MPI): each clump owns a disjoint,
    # contiguous slice of gridcells/columns/patches and there is no halo, so the
    # clumps are embarrassingly parallel with zero communication.
    #
    # `get_clump_bounds(nc)` returns PROC-RELATIVE bounds (offset by procinfo.beg*)
    # so a clump's begc:endc indexes directly into the proc-level state arrays, and
    # the proc-wide filter masks (filt.soilc, ...) are absolute-indexed — so the
    # SAME `filt` is correct for any clump sub-range with no slicing.
    #
    # `clm_drv_per_clump_bounds` is the list of clump bounds to walk. With a single
    # clump (the default, `nclumps == 1`) this is exactly `[bounds_proc]`, so the
    # body below runs once over the whole proc domain — BYTE-IDENTICAL to the
    # historical single-clump path. With `nclumps > 1` it partitions the proc into
    # disjoint clumps; `clm_run_clump_physics!` runs the clump-safe per-column
    # physics over them (threaded when `config.use_threaded_clumps`).
    #
    # THREAD-SAFETY NOTE: the orchestrator body below is NOT itself threaded over
    # clumps. It contains full-proc-domain operations (`for c in 1:length(...)`
    # FATES site walks with a running site counter `s`, lazy `config.accum_mgr`
    # initialization, gridcell-scratch sized to the whole proc) that are not
    # clump-decomposable without a broader refactor and would race under naive
    # whole-body threading. The orchestrator therefore runs ONCE over the proc
    # (bounds_clump = bounds_proc), exactly as before. Genuinely embarrassingly-
    # parallel per-column physics is exposed via `clm_run_clump_physics!`, which is
    # threaded over clumps and proven deterministic (multi-clump result ==
    # single-clump result regardless of clump count / thread count).
    clm_drv_per_clump_bounds = if config.use_threaded_clumps && get_proc_clumps() > 1
        cb = [get_clump_bounds(nc) for nc in 1:get_proc_clumps()]
        # Only adopt the clump partition if it exactly tiles the proc domain passed
        # in (disjoint + complete column coverage) — otherwise the global `decomp`
        # was split inconsistently with this `bounds_proc`, so fall back to the safe
        # single whole-proc clump rather than mis-slice the state arrays.
        proc_ncol = bounds_proc.endc - bounds_proc.begc + 1
        sum(b -> b.endc - b.begc + 1, cb) == proc_ncol ? cb : [bounds_proc]
    else
        [bounds_proc]
    end

    bounds_clump = bounds_proc  # orchestrator runs proc-level (single clump = whole proc)

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
    l2a = inst.lnd2atm
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
    # Backend-aware scratch/parameter helpers. On CPU `_bref` is a plain Array so
    # these are ordinary host allocations; on a GPU device run `_bref` is a device
    # array so `similar` lands the scratch on-device and `_onbk` copies global
    # (host) parameter arrays — e.g. pftcon.* — onto the working backend at FT.
    _bref = a2l.forc_pbot_downscaled_col
    _zlike(dims...)       = fill!(similar(_bref, FT, dims...), zero(FT))
    _olike(dims...)       = fill!(similar(_bref, FT, dims...), one(FT))
    _fulllike(v, dims...) = fill!(similar(_bref, FT, dims...), FT(v))
    _blike(dims...)       = fill!(similar(_bref, Bool, dims...), false)
    _onbk(arr)            = _bref isa Array ? arr : Adapt.adapt(typeof(_bref).name.wrapper, FT.(arr))
    # Match the timestep scalar to the working precision (Float64 on CPU, Float32
    # on device) so it doesn't leak a `double` into the device kernels it's passed to.
    dtime = FT(dtime)
    # Column specific humidity from vapor pressure — KernelAbstractions kernel
    # (runs on CPU or GPU depending on where the arrays live). See kernels.jl.
    forc_q_col = _zlike(nc)
    compute_forc_q!(forc_q_col, col.gridcell, a2l.forc_vp_grc, a2l.forc_pbot_downscaled_col)

    # ========================================================================
    # Glacier area initialization on first timestep
    # ========================================================================
    if is_first_step
        # Placeholder: update_glc2lnd_fracs!(bounds_clump) [glacier coupling]

        # dynSubgrid_wrapup_weight_changes! — WIRED (gated; no-op unless a transient
        # run built config.dyn_subgrid in clm_initialize!). Reconciles the cold-start
        # subgrid weights up the hierarchy before the first physics step.
        if config.dyn_subgrid !== nothing
            dynSubgrid_wrapup_weight_changes!(bc, grc, lun, col, pch)
        end

        # Irrigation cold-start — WIRED (gated on config.irrigate). The relsat
        # wilting-point/target columns + per-patch surface method + nsteps-per-day are
        # static, so initialize them once on the first step (rather than threading a
        # full IrrigationInitCold through clm_initialize!). Uses the default
        # Clapp-Hornberger SWRC and the PFT irrigated flag; surface method left UNSET
        # so set_irrig_method! falls back to the default (drip). No-op when irrigate off.
        if config.irrigate && !isempty(irr.irrig_method_patch) && !isempty(pftcon.irrigated) && irr.irrig_nsteps_per_day == 0
            _irr_swrc = SoilWaterRetentionCurveClappHornberg1978()
            ng_irr = isempty(bc_grc) ? 0 : length(bc_grc)
            npft_irr = length(pftcon.irrigated)
            _irr_method_surf = fill(IRRIG_METHOD_UNSET, max(ng_irr, 1), max(npft_irr, 1))
            irrigation_init_cold!(irr, ss, _irr_swrc, col,
                                  pftcon.irrigated, _irr_method_surf,
                                  pch, grc, bc_col, bc_patch,
                                  varpar.nlevsoi, round(Int, dtime))
        end
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
            soilbgc_ns=inst.soilbiogeochem_nitrogenstate,
            # c2g + BeginCNGridcellBalance (seeds begcb_grc/begnb_grc).
            bounds_grc=bc_grc,
            cascade_con=(_decomp_initialized(inst.decomp_cascade) ? inst.decomp_cascade : nothing),
            col=col, patch=pch,
            nlevdecomp=varpar.nlevdecomp, ndecomp_pools=config.ndecomp_pools,
            dzsoi_decomp=dzsoi_decomp[], zisoi_vals=zisoi[])
    end
    if config.use_lch4 && !isempty(inst.ch4.ch4_surf_flux_tot_col)
        # ch4_init_gridcell_balance_check! — WIRED. dz is soil-indexed (drop the
        # nlevsno snow-layer prefix that col.dz/z/zi/t_soisno carry). Seeds
        # totcolch4_bef_grc for the end-of-step CH4 conservation check.
        _ch4_nsno = varpar.nlevsno; _ch4_nsoi = varpar.nlevsoi
        _ch4_dz_soil = Matrix(col.dz[:, (_ch4_nsno + 1):(_ch4_nsno + _ch4_nsoi)])
        ch4_init_gridcell_balance_check!(inst.ch4, filt.nolakec, filt.lakec,
                                         Array(col.gridcell), Array(col.wtgcell),
                                         _ch4_dz_soil, _ch4_nsoi, bc.endg,
                                         inst.ch4_varcon.allowlakeprod)
    end

    # Begin water balance — WIRED
    water_gridcell_balance!(inst.water, ls, col, lun, grc,
                            filt.nolakec, filt.lakec,
                            bc_col, bc_lun, bc_grc, "begwb")
    # Canopy water belongs in begwb_grc, evaluated HERE with the BEGIN-of-step
    # canopy state (Fortran ComputeWaterMassNonLake p2c's canopy into the column
    # total that WaterGridcellBalance then c2g's; the Julia storage fn stubs it to
    # zero). This used to be done at the END of the step for begwb_grc AND
    # endwb_grc from the same end-of-step liqcan/snocan, so the canopy term
    # cancelled out of (endwb_grc - begwb_grc) — the gridcell balance omitted the
    # canopy storage change while qflx_evap_tot/forc_rain still carried the canopy
    # exchanges. That left a real ~0.006 mm/step residual in errh2o_grc which only
    # ever @warn'd because DAnstep was hardcoded to 0. Mirrors begwb_col below.
    add_canopy_water_to_grc_storage!(inst.water.waterbalancebulk_inst.begwb_grc,
                                     wsb.ws.liqcan_patch, wsb.ws.snocan_patch,
                                     filt.nolakec, col, lun, pch, bc_col, bc_grc)

    # ========================================================================
    # Dynamic subgrid weights — WIRED (gated; no-op unless a transient run built
    # config.dyn_subgrid). Applies the per-year land-cover change and the
    # biogeophysics (water + energy) conservation dynbal via the cons_bgp bundle.
    # Transient datasets advance on year boundaries, so only run at is_beg_curr_year
    # (matching the Fortran do_transient_*'s annual cadence). The CN-half
    # (dyn_cnbal_patch!/col!) is WIRED via the cons_bgc bundle when config.use_cn —
    # mirroring the Fortran `if (use_cn) call DynamicAreaConservation`.
    # ========================================================================
    if config.dyn_subgrid !== nothing && is_beg_curr_year && !is_first_step
        _cons_bgp = (urbanparams = up, soilstate = ss,
                     waterstatebulk = wsb, waterdiagnosticbulk = wdb,
                     temperature = temp, lakestate = ls)
        # The dyn_subgrid driver asserts PROC-level bounds; in single-clump mode the
        # clump bounds ARE the proc bounds, so present a PROC-level view of bc.
        _bc_proc = BoundsType(begg = bc.begg, endg = bc.endg, begl = bc.begl,
                              endl = bc.endl, begc = bc.begc, endc = bc.endc,
                              begp = bc.begp, endp = bc.endp,
                              level = BOUNDS_LEVEL_PROC)

        # CN-balance bundle (use_cn only). Mirrors CNVegetationFacade's
        # DynamicAreaConservation arg list. The soilp/soilc masks must include
        # inactive points (Fortran filter_inactive_and_active) so columns that
        # just shrank to 0 area are still updated.
        _cons_bgc = nothing
        if config.use_cn
            _veg = inst.bgc_vegetation
            _cons_bgc = (
                dynbal = DynConsBiogeochemState(),
                pftcon = pftcon,
                canopystate = inst.canopystate,
                cnveg_state = _veg.cnveg_state_inst,
                cnveg_carbonstate = _veg.cnveg_carbonstate_inst,
                cnveg_carbonflux = _veg.cnveg_carbonflux_inst,
                cnveg_nitrogenstate = _veg.cnveg_nitrogenstate_inst,
                cnveg_nitrogenflux = _veg.cnveg_nitrogenflux_inst,
                soilbiogeochem_state = inst.soilbiogeochem_state,
                soilbiogeochem_carbonstate = inst.soilbiogeochem_carbonstate,
                soilbiogeochem_nitrogenstate = inst.soilbiogeochem_nitrogenstate,
                mask_soilp_with_inactive = filt_inactive_and_active.soilp,
                mask_soilc_with_inactive = filt_inactive_and_active.soilc,
                dt = dtime,
                nlevdecomp = varpar.nlevdecomp,
                ndecomp_pools = config.ndecomp_pools,
                i_litr_min = config.i_litr_min,
                i_litr_max = config.i_litr_max,
                i_cwd = config.i_cwd,
                use_crop = config.use_crop,
                nrepr = config.nrepr,
                use_nitrif_denitrif = _veg.config.use_nitrif_denitrif,
            )
        end

        dynSubgrid_driver!(config.dyn_subgrid, _bc_proc, grc, lun, col, pch;
                           year = year,
                           temp = temp, ws = wsb.ws, cons_bgp = _cons_bgp,
                           cons_bgc = _cons_bgc,
                           mask_nolakec = filt.nolakec, mask_lakec = filt.lakec)
    end

    # ------------------------------------------------------------------------
    # CNDV: time-interpolate dynamic-vegetation PFT weights to this time step.
    # Mirrors CNVegetationFacade::UpdateSubgridWeights → dynCNDV_interp (called
    # every time step when use_cndv); see dynCNDVMod.F90:dynCNDV_interp. The
    # interpolation weight is wt1 = 1 - yearfrac(this step); fpcgridold is rolled
    # forward on the first step of the year (is_beg_curr_year && nstep>0).
    # Gated off by default → no-op on the empty/unsized dgvs state.
    # ------------------------------------------------------------------------
    if config.use_cndv && !isempty(inst.dgvs.fpcgrid_patch)
        _cndv_yearfrac = (jday - 1 + secs / SECSPDAY) / 365.0
        dyn_cndv_interp!(inst.dgvs, pch, lun, bc_patch;
                         wt1 = 1.0 - _cndv_yearfrac,
                         is_beg_curr_year = is_beg_curr_year && nstep > 0)
    end

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
    # Include p2c canopy water in begwb (Fortran ComputeWaterMassNonLake p2c's it
    # in; the Julia storage fn stubs canopy to zero). Required so errh2o_col treats
    # canopy interception/unloading as internal, not a phantom storage change.
    add_canopy_water_to_storage!(inst.water.waterbalancebulk_inst.begwb_col,
                                 wsb.ws.liqcan_patch, wsb.ws.snocan_patch,
                                 filt.nolakec, col, pch)

    if config.use_cn || config.use_fates_bgc
        # InitColumnBalance — WIRED
        cn_vegetation_init_column_balance!(inst.bgc_vegetation;
            mask_bgc_soilc=filt.bgc_soilc,
            mask_bgc_vegp=filt.bgc_vegp,
            mask_allc=filt.allc,
            bounds_col=bc_col,
            bounds_patch=bc_patch,
            soilbgc_cs=inst.soilbiogeochem_carbonstate,
            soilbgc_ns=inst.soilbiogeochem_nitrogenstate,
            cascade_con=(_decomp_initialized(inst.decomp_cascade) ? inst.decomp_cascade : nothing),
            col=col, patch=pch,
            nlevdecomp=varpar.nlevdecomp, ndecomp_pools=config.ndecomp_pools,
            dzsoi_decomp=dzsoi_decomp[], zisoi_vals=zisoi[])
    end
    if config.use_lch4 && !isempty(inst.ch4.ch4_surf_flux_tot_col)
        # ch4_init_column_balance_check! — WIRED. Seeds totcolch4_bef_col for the
        # column-level CH4 conservation check (dz soil-indexed as above).
        _ch4c_nsno = varpar.nlevsno; _ch4c_nsoi = varpar.nlevsoi
        _ch4_dz_soil_c = Matrix(col.dz[:, (_ch4c_nsno + 1):(_ch4c_nsno + _ch4c_nsoi)])
        ch4_init_column_balance_check!(inst.ch4, filt.nolakec, filt.lakec,
                                       _ch4_dz_soil_c, _ch4c_nsoi,
                                       inst.ch4_varcon.allowlakeprod)
    end

    # ========================================================================
    # Update dynamic N deposition / forcing interpolation
    # ========================================================================
    if config.use_cn && !config.ndep_from_cpl
        # ndep_interp! — WIRED. Fills atm2lnd.forc_ndep_grc (gN/m2/s) from the
        # N-deposition stream. Fortran clm_drv: `if (use_cn .and. .not. ndep_from_cpl)
        # call ndep_interp(bounds_proc, atm2lnd_inst)` (ndepStreamMod.F90).
        # No-op when no stream file was supplied (ndep_stream.active == false), which
        # leaves forc_ndep_grc at its initialised value exactly as before.
        # dayspyr = 365 (NO_LEAP), the convention used throughout the CN driver.
        ndep_interp!(inst.atm2lnd, inst.ndep_stream, grc;
                     calday     = jday + secs / SECSPDAY,
                     dayspyr    = 365.0,
                     bounds_grc = bc_grc)
    end
    if config.use_cn
        # bgc_vegetation_inst%InterpFileInputs → CNVegetationFacade::InterpFileInputs
        # → cnfire_method%FireInterp(bounds) (FireDataBaseType.F90). WIRED: fills
        # cnfire_li2014.forc_lnfm (counts/km2/hr) + forc_hdm (counts/km2) — the
        # lightning + population-density ignition sources of the Li fire schemes.
        # Fortran runs this only when need_lightning_and_popdens() (i.e. a Li method);
        # no-op when the streams were never read (no file supplied), so the default
        # (:nofire) path is untouched.
        if need_lightning_and_popdens(config.cnfire_method)
            fire_stream_interp!(inst.cnfire_li2014, inst.cnfire_stream, grc;
                                calday     = jday + secs / SECSPDAY,
                                dayspyr    = 365.0,
                                bounds_grc = bc_grc)
        end
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

    # Derive forc_rh_grc from the CURRENT step's forcing (t/pbot/q). Fortran gets
    # forc_rh from the coupler in the atm->lnd import at the TOP of the step; this
    # port derives it, and used to do so only in the end-of-step accumulator block
    # (it is still called there, immediately before atm2lnd_update_acc_vars!, so the
    # accumulated RH30/RH24 are bit-identical — same value, just computed earlier).
    # Deriving it here as well is what lets the Li fire schemes see THIS step's RH
    # instead of the previous step's (0 on step 1, i.e. maximum flammability).
    atm2lnd_update_rh!(a2l, bc_grc)

    # Update exposed vegetation filter
    set_exposedvegp_filter!(filt, bc,
                            cs.frac_veg_nosno_patch[bc_patch])

    # ========================================================================
    # Irrigation withdrawal — WIRED (gated on config.irrigate)
    # ========================================================================
    # CalcAndWithdrawIrrigationFluxes: route the demand computed at the END of the
    # previous step (calc_irrigation_needed!) into surface/groundwater withdrawal +
    # the per-patch drip/sprinkler application fluxes consumed by canopy interception
    # downstream. calc_irrigation_fluxes! internally runs the bulk-withdrawal and
    # (guard: !isempty(irrig_method_patch) — no-op on an unsized instance, like ch4)
    # application-flux sub-steps. No-op until n_irrig_steps_left_patch is set.
    if config.irrigate && !isempty(irr.irrig_method_patch) && !isempty(pftcon.irrigated)
        calc_irrigation_fluxes!(irr, sh, ss, wfb, col, pch,
                                filt.soilc, filt.soilp,
                                bc_col, bc_patch, varpar.nlevsoi, dtime)
    end

    # ========================================================================
    # HYDROLOGY STAGE 1: Canopy interception, new snow, surface water
    # ========================================================================
    # Rain/snow/flood forcings (from atmosphere, filled by downscale_forcings! or forcing reader)
    np = length(pch.column)
    forc_rain_col = isdefined(Main, :nothing) ? _zlike(nc) : _zlike(nc)  # will use a2l fields
    forc_snow_col = _zlike(nc)
    qflx_floodg = _zlike(length(grc.lat))
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
        _zlike(np), _zlike(np))  # irrigation sprinkler/drip placeholders

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
        # W3a sun/shade — FATES.  Pack the radiation bc_in for every FATES site,
        # run FatesSunShadeFracs, and unpack fsun/laisun/laisha back to canopystate.
        # Maps FATES site s 1:1 onto the s-th FATES column; each of the site's
        # vegetated patches maps to HLM patch p = ifp + col.patchi[c] (bare-ground
        # at +0). The fates_veg_patches walk yields the (ifp, p) pairs.
        if inst.fates !== nothing
            s = 0
            for c in 1:length(col.is_fates)
                col.is_fates[c] || continue
                s += 1
                s <= inst.fates.nsites || break
                for (ifp, p) in fates_veg_patches(inst.fates.sites[s], c, col)
                    fates_pack_bcin_radiation!(inst; s=s, c=c, p=p, ifp=ifp,
                                               coszen=alb.coszen_col[c])
                end
            end
            FatesSunShadeFracs(inst.fates.nsites, inst.fates.sites,
                               inst.fates.bc_in, inst.fates.bc_out)
            s = 0
            for c in 1:length(col.is_fates)
                col.is_fates[c] || continue
                s += 1
                s <= inst.fates.nsites || break
                for (ifp, p) in fates_veg_patches(inst.fates.sites[s], c, col)
                    fates_unpack_bcout_sunfrac!(inst; s=s, c=c, p=p, ifp=ifp)
                end
            end
        end
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
                       bc_patch, pch, _onbk(pftcon.woody);
                       is_time_to_run_luna=false)

    # CalculateSurfaceHumidity — WIRED
    calculate_surface_humidity!(col, lun, temp, ss, wsb, wdb,
                                a2l.forc_pbot_downscaled_col, forc_q_col,
                                filt.nolakec, bc_col)

    # ========================================================================
    # DETERMINE FLUXES
    # ========================================================================


    # BareGroundFluxes — WIRED
    # CH4 aerenchyma boundary conductance (grnd_ch4_cond_patch): BareGroundFluxesMod.F90:401
    # writes `1/raw` DIRECTLY into ch4_inst%grnd_ch4_cond_patch. The Julia driver used to
    # omit it, leaving ch4.grnd_ch4_cond_patch frozen at its cold-start 0.01 (which starved
    # the aerenchyma boundary conductance in ch4!). Thread the ch4 struct field so the bare
    # fraction is populated, exactly as Fortran does. Gated on use_lch4 (the write in
    # bareground_fluxes! is itself guarded by use_lch4), so the non-CH4 path is unchanged.
    _bare_ch4_cond = config.use_lch4 ? inst.ch4.grnd_ch4_cond_patch : Float64[]
    bareground_fluxes!(cs, ef, fv, temp, ss, wfb, wsb, wdb, ps,
                       pch, col, lun,
                       filt.noexposedvegp, bc_patch,
                       forc_q_col, a2l.forc_pbot_downscaled_col,
                       a2l.forc_th_downscaled_col, a2l.forc_rho_downscaled_col,
                       a2l.forc_t_downscaled_col,
                       a2l.forc_u_grc, a2l.forc_v_grc,
                       a2l.forc_hgt_t_grc, a2l.forc_hgt_u_grc, a2l.forc_hgt_q_grc;
                       use_lch4 = config.use_lch4,
                       grnd_ch4_cond_patch = _bare_ch4_cond)

    # Volumetric liquid water content for root moisture stress — 2D kernel
    # over (column, soil layer). See kernels.jl.
    joff = varpar.nlevsno
    # Refresh effective soil porosity (watsat − vol_ice) from THIS timestep's ice
    # before it is consumed below. Fortran SoilMoistStressMod calls
    # calc_effective_soilporosity at the top of the BTRAN routine, so eff_porosity is
    # always self-consistent with the current ice. The Julia port otherwise only
    # refreshes eff_porosity_col later in the hydrology block (set_soil_water_fractions!,
    # after this flux block), leaving a ONE-STEP-STALE value here that biased the
    # freeze-up/thaw shoulder-season chain liqvol-clamp → smp_l → vegwp → bsun/btran
    # (and transpiration) high. (calc_effective_soilporosity! was dead code.)
    calc_effective_soilporosity!(ss.watsat_col, wsb.ws.h2osoi_ice_col, col.dz,
                                 ss.eff_porosity_col, filt.nolakec, bc_col,
                                 varpar.nlevgrnd, varpar.nlevsno)
    compute_h2osoi_liqvol!(wdb.h2osoi_liqvol_col, filt.nolakec, col.dz,
                           wsb.ws.h2osoi_liq_col, ss.eff_porosity_col,
                           joff, varpar.nlevgrnd; denh2o=DENH2O)

    # SoilMoistStress (root moisture stress / BTRAN) — WIRED
    al = inst.active_layer
    calc_root_moist_stress!(ss, ef, temp, wsb, wdb, col, pch,
                             _onbk(pftcon.smpso), _onbk(pftcon.smpsc),
                             FT.(al.altmax_lastyear_indx_col),
                             FT.(al.altmax_indx_col),
                             filt.exposedvegp, bc_patch,
                             varpar.nlevgrnd, varpar.nlevsno)

    # W4a btran — FATES.  For FATES columns, pack the soil-hydrology bc_in and run
    # btran_ed! (the FATES water-stress factor + root-soil resistance), then unpack
    # btran/rootr back into energyflux/soilstate (overriding the column average).
    if config.use_fates && inst.fates !== nothing
        nlevsoil_f = varpar.nlevsoi
        s = 0
        for c in 1:length(col.is_fates)
            col.is_fates[c] || continue
            s += 1
            s <= inst.fates.nsites || break
            # Fortran CanopyFluxes hands wrap_btran the COLUMNS of the exposed-veg
            # patch filter (`filterc_tmp(f) = patch%column(filterp(f))`); columns
            # outside it get `filter_btran = .false.` and a -999 soil state, because
            # h2osoi_liqvol is not computed for them.
            exposed = false
            for p in eachindex(filt.exposedvegp)
                if filt.exposedvegp[p] && pch.column[p] == c
                    exposed = true
                    break
                end
            end
            fates_pack_bcin_btran!(inst; s=s, c=c, nlevsoil=nlevsoil_f,
                                   filter_btran=exposed)
        end
        btran_ed!(inst.fates.sites, inst.fates.bc_in, inst.fates.bc_out)
        s = 0
        for c in 1:length(col.is_fates)
            col.is_fates[c] || continue
            s += 1
            s <= inst.fates.nsites || break
            for (ifp, p) in fates_veg_patches(inst.fates.sites[s], c, col)
                fates_unpack_bcout_btran!(inst; s=s, c=c, p=p, ifp=ifp,
                                          nlevsoil=nlevsoil_f)
            end
        end
    end

    # CanopyFluxes — WIRED
    # Get downreg/leafn from CN vegetation facade when active, else zeros for SP mode
    if config.use_cn
        downreg_patch = get_downreg_patch(inst.bgc_vegetation, bc_patch)
        leafn_patch = get_leafn_patch(inst.bgc_vegetation, bc_patch)
    else
        downreg_patch = _zlike(np)
        leafn_patch = _zlike(np)
    end
    # PHS (plant hydraulic stress) needs fine-root carbon for the soil-to-root
    # conductance (root_biomass_density → root_length_density → r_soil → soil_cond).
    # Fortran computes it for BOTH bgc modes (CNVegetationFacade::get_froot_carbon_patch):
    # CN uses frootc_patch; SP derives a diagnostic from prescribed LAI
    # (froot_c = tlai/slatop * froot_leaf). Gating this on use_cn left SP+PHS with
    # froot_carbon=0 → root_biomass_density floored to c_to_b → r_soil ~70x too large
    # → soil-to-root conductance ~48x too small → over-stressed stomata (rss pinned to
    # rsmax). Compute it whenever PHS is on, matching Fortran for either mode.
    phs_froot_c = if config.use_hydrstress
        config.use_cn ?
            get_froot_carbon_patch(inst.bgc_vegetation, bc_patch) :
            get_froot_carbon_patch(inst.bgc_vegetation, bc_patch;
                tlai=cs.tlai_patch, slatop=pftcon.slatop,
                froot_leaf=pftcon.froot_leaf, ivt=(pch.itype .+ 1))
    else
        FT[]
    end
    # Calibrated stomatal params (per-PFT medlynslope/intercept) come from pftcon in
    # CN mode AND whenever plant-hydraulic stress (PHS) is active. Fortran ALWAYS reads
    # pftcon%medlynslope (regardless of use_cn), so the SP+PHS Bow reference run uses the
    # calibrated per-PFT slope (e.g. 11.15 for the Bow tree/grass PFTs), not the
    # canopy-core default 6.0. Hardcoding 6.0 here made Julia's unstressed stomatal
    # conductance ~1.75x too low → too little transpiration demand → bsun too high
    # (under-stressed) → the residual PHS T_VEG divergence. Gating the read on
    # use_cn||use_hydrstress restores parity while keeping the pure non-PHS SP/AD path
    # (neither flag) byte-identical at medlynslope=6.0, preserving its validated energy
    # balance. (Analogous to the SP-mode froot_carbon #138 and LUNA lmr #139 param fixes.)
    _use_calib_med = config.use_cn || config.use_hydrstress
    _medint_pft = _use_calib_med ? _onbk(pftcon.medlynintercept) : fill(100.0, MXPFT + 1)
    _medslp_pft = _use_calib_med ? _onbk(pftcon.medlynslope)     : fill(6.0, MXPFT + 1)
    _crop_pft   = config.use_cn ? _onbk(pftcon.crop)            : Float64[]
    # FATES: reset the photosynthesis filter to 1 ("do nothing") at the START of the
    # canopy-flux/photosynthesis step (mirrors Fortran filter_photo_pa(:)=1 in the
    # CanopyFluxes photosynthesis setup). The in-solve pack re-flags computing patches to
    # 2; the post-solve accumulate transitions them 2→3. Resetting here (not after the
    # accumulate) means a step that skips photosynthesis leaves the filter at 1 instead of
    # a stale 3, so AccumulateFluxes_ED does not re-accumulate stale per-step fluxes.
    if config.use_fates && inst.fates !== nothing
        for s in 1:inst.fates.nsites
            fill!(inst.fates.bc_in[s].filter_photo_pa, 1)
        end
    end

    # PhotosynthesisMod::TimeStepInit — zero psnsun/psnsha/fpsn (+ the wc/wj/wp
    # splits and the C13/C14 twins) on EVERY NON-LAKE patch at the start of the
    # canopy-flux step. Fortran calls this at the top of CanopyFluxes
    # (CanopyFluxesMod.F90:661).
    #
    # `photosynthesis_timestep_init!` was PORTED (photosynthesis.jl) and then
    # never called — the same dead-port class as the `*_init_cold!` sweep in
    # src/infrastructure/init_cold.jl. The consequence: a BARE-GROUND patch never
    # enters the photosynthesis solve, so nothing ever wrote its psnsun/psnsha/
    # fpsn and they sat at the allocator's NaN for the whole run. That NaN was
    # being papered over downstream by an `isfinite(v) ? v : 0.0` fallback in
    # history_fpsn_patch (history_writer.jl) — a guard added to fix the +43% FPSN
    # aggregation bug whose ROOT CAUSE was this missing call. Zeroing at the
    # source is what Fortran does; the history fallback stays only to map the
    # LAKE patches' deliberate spval to 0.
    photosynthesis_timestep_init!(ps, filt.nolakep, bc_patch;
                                  use_c13 = config.use_c13,
                                  use_c14 = config.use_c14)

    # Positional call into canopy_fluxes_core! (no kwarg NamedTuple) so Enzyme
    # reverse-mode can compile the differentiated driver. The trailing args after
    # `overrides` are the "default-only back group". We thread the ones that affect
    # the stomatal/photosynthesis solve — the calibrated Medlyn params, crop flags,
    # and config.use_cn (which selects the CN leaf-respiration / GPP path). Without
    # this the canopy ran with medlynslope=6.0 (Bow-calibrated 11.15) and use_cn=false
    # even in a CN run. Passing them positionally keeps the differentiated path
    # compilable. The woody-structure (dbh/nstem/is_tree/...) and z0v params keep
    # their core defaults — they only matter under use_biomass_heat_storage (off) or
    # the non-default z0param_method.
    canopy_fluxes_core!(cs, ef, fv, temp, sa, ss, wfb, wsb, wdb, ps,
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
                   dtime,
                   temp.t_a10_patch, alb.nrad_patch, alb.tlai_z_patch,
                   alb.vcmaxcintsun_patch, alb.vcmaxcintsha_patch,
                   sa.parsun_z_patch, sa.parsha_z_patch,
                   cs.laisun_z_patch, cs.laisha_z_patch,
                   _olike(np), _olike(np),
                   _onbk(pftcon.dleaf), _onbk(pftcon.slatop), _onbk(pftcon.leafcn),
                   _onbk(pftcon.flnr), _onbk(pftcon.fnitr), _onbk(pftcon.mbbopt),
                   _onbk(pftcon.c3psn), _onbk(pftcon.woody),
                   inst.overrides,
                   # (3) default-only back group (positional): woody-structure +
                   # z0v keep core defaults; Medlyn/crop/use_cn are threaded.
                   fill(0.1, MXPFT + 1), fill(0.1, MXPFT + 1), fill(1.0, MXPFT + 1),
                   fill(0.0, MXPFT + 1), fill(500.0, MXPFT + 1),
                   fill(false, MXPFT + 1), fill(false, MXPFT + 1),
                   _medint_pft, _medslp_pft, _crop_pft,
                   fill(0.35, MXPFT + 1), fill(0.003, MXPFT + 1), fill(0.25, MXPFT + 1),
                   fill(2.0, MXPFT + 1), fill(8.0, MXPFT + 1),
                   config.use_cn,
                   config.use_lch4, false, config.use_hydrstress, phs_froot_c,
                   config.use_fates, config.use_luna,
                   # Re-state the four default-only positionals between use_luna and
                   # the FATES handle. z0param_method / forc_pc13o2_grc / leaf_mr_vcm keep
                   # their canopy_fluxes_core! defaults. grnd_ch4_cond_patch is threaded
                   # from the ch4 struct under use_lch4: CanopyFluxesMod.F90:1097 writes
                   # `1/(raw_above+raw_below)` DIRECTLY into ch4_inst%grnd_ch4_cond_patch
                   # (the exposed-veg counterpart of the bareground write above). Under
                   # use_lch4=false it is Float64[], so the default/AD path is byte-identical.
                   "ZengWang2007",
                   (config.use_lch4 ? inst.ch4.grnd_ch4_cond_patch : Float64[]),
                   Float64[], 0.015,
                   # W4b: FATES handle — threaded as the trailing optional positional.
                   # In Fortran, FATES photosynthesis is called from INSIDE the
                   # CanopyFluxes iterative solve (CanopyFluxesMod:1117). The gated
                   # `use_fates` branch in canopy_fluxes_core! now runs
                   # FatesPlantRespPhotosynthDrive each Newton iteration with the
                   # current in-loop leaf temperature and unpacks rssun/rssha back
                   # into photosyns (two-way coupling). The handle is `nothing` on the
                   # non-FATES path, so the differentiable default signature is
                   # untouched. This supersedes the former adjacent post-solve call.
                   (config.use_fates ? inst.fates : nothing))

    # Total canopy water diagnostic h2ocan = snocan + liqcan, synced HERE (after
    # canopy_fluxes applied the canopy-evaporation decrement to snocan/liqcan) so it
    # matches Fortran's end-of-step H2OCAN sampling. Syncing it at the end of
    # canopy_hydrology instead (pre-evaporation) reads the store ~15% high.
    canhyd_update_h2ocan!(wdb.h2ocan_patch, wsb.ws.snocan_patch, wsb.ws.liqcan_patch,
                          filt.nolakep, 1, np)

    # C13/C14 photosynthesis fractionation — WIRED (gated on config.use_c13/use_c14).
    # Runs after the photosynthesis solve has filled ps.psnsun/psnsha + the
    # fractionation factors ps.alphapsnsun/sha, mirroring the Fortran call from
    # PhotosynthesisTotal. With both flags off (default) the block is skipped, so the
    # default path is unchanged; non-veg patches (psn=0) contribute zero isotope flux.
    if config.use_c13 || config.use_c14
        c13_c14_photosynthesis!(ps, a2l.forc_pco2_grc, a2l.forc_pc13o2_grc,
            pch.gridcell, grc.latdeg, trues(np), 1:np;
            use_c13 = config.use_c13, use_c14 = config.use_c14)
    end

    # CalcOzoneUptake — WIRED (gated on config.use_ozone)
    # Accumulate ozone dose now that the canopy/photosynthesis solve has filled the
    # stomatal resistances (ps.rssun/rssha) and the leaf-boundary/aerodynamic
    # resistances (fv.rb1/ram1). The matching CalcOzoneStress runs at the TOP of the
    # next timestep's flux block (it intentionally uses the previous step's dose, per
    # OzoneMod's design). With use_ozone off, o3uptake* stays at its cold-start 0 and
    # CalcOzoneStress returns coef 1.0 → the default path is byte-identical.
    if config.use_ozone
        calc_ozone_uptake!(oz, pch, filt.exposedvegp, bc_patch,
                           a2l.forc_pbot_downscaled_col, a2l.forc_th_downscaled_col,
                           ps.rssun_patch, ps.rssha_patch,
                           fv.rb1_patch, fv.ram1_patch,
                           cs.tlai_patch, a2l.forc_o3_grc,
                           _onbk(pftcon.evergreen), _onbk(pftcon.leaf_long),
                           dtime)
    end

    # W4b photosynthesis — FATES: the per-timestep ("hifrq") history fill +
    # plant-hydraulics step that used to follow the (now-removed) adjacent
    # photosynthesis call. The photosynthesis solve itself is done IN-LOOP inside
    # canopy_fluxes_core! above (use_fates branch); only the post-solve
    # diagnostics + hydraulics remain here.
    if config.use_fates && inst.fates !== nothing
        # (W4b: the FATES photosynthesis solve is done IN-LOOP inside
        # canopy_fluxes_core! above — its use_fates branch runs
        # FatesPlantRespPhotosynthDrive for each canopy-solve patch, so the
        # multi-patch ifp-walk is handled there by the per-patch canopy loop.)

        # W4b-accumulate — FATES: mirror the Fortran wrap_accumulatefluxes, called at
        # the END of CanopyFluxes (CanopyFluxesMod:1634). The in-solve photosynthesis
        # re-packs filter_photo_pa=2 each Newton iteration and computes FRESH per-step
        # cohort fluxes (gpp/npp/resp_tstep; 0 at night). Transition the patches that
        # computed this step (==2 → 3, the Fortran wrap_photosynthesis 2→3 step) and run
        # AccumulateFluxes_ED (which gates on ==3) so the per-step fluxes roll up into the
        # DAILY gpp/npp/resp_acc that the daily growth+allocation step (EDMainMod) reads.
        # WITHOUT this accumulation npp_acc ≡ 0 → the allocation grows nothing (dbh frozen,
        # stand can only lose carbon). The 1-reset happens at the START of the next step's
        # canopy solve (mirrors Fortran filter_photo_pa(:)=1), so a step that skips
        # photosynthesis leaves the filter at 1 and is not re-accumulated with stale fluxes.
        let s = 0
            for c in 1:length(col.is_fates)
                col.is_fates[c] || continue
                s += 1
                s <= inst.fates.nsites || break
                fbc = inst.fates.bc_in[s]
                for (ifp, _p) in fates_veg_patches(inst.fates.sites[s], c, col)
                    fbc.filter_photo_pa[ifp] == 2 && (fbc.filter_photo_pa[ifp] = 3)
                end
            end
        end
        AccumulateFluxes_ED(inst.fates.nsites, inst.fates.sites,
                            inst.fates.bc_in, inst.fates.bc_out, dtime)

        # W4b-runningmeans — FATES: mirror the Fortran UpdateFatesRMeansTStep, called on
        # every model timestep after the photosynthesis solve. It rolls the veg-temperature
        # running means (tveg24 / tveg_lpa / tveg_longterm) from bc_in.t_veg_pa, the seedling
        # PAR/SMP/MDD means, and the smoothed site NPP (bc_out.ema_npp). WITHOUT this the
        # temperature means stay FROZEN at their cold-start init (temp_init_veg = 288.15 K =
        # 15 °C): the FATES cold-deciduous phenology then reads a constant, wrong 15 °C veg
        # temperature and eventually fires a spurious synchronized leaf-off that never re-
        # flushes → carbon starvation → a phantom boom/bust die-back on multi-year runs.
        UpdateFatesRMeansTStep(inst.fates.sites, inst.fates.bc_in, inst.fates.bc_out)

        # Per-timestep ("hifrq") FATES history fill — site GPP/AR/NEP/resp +
        # stomatal/bl conductance + veg-temp + radiation-error diagnostics from
        # the just-computed cohort per-step fluxes. Write-only; no-op until the
        # history interface is built (clm_fates_init!). Mirrors the Fortran host's
        # update_history_hifrq in the timestep path.
        fates_hifrq_history_step!(inst; dt_tstep=dtime)

        # W4c plant-hydraulics — FATES.  When plant hydraulics is enabled
        # (hlm_use_planthydro==itrue), run hydraulics_drive AFTER the photosynthesis
        # solve (mirroring the Fortran call from the photosynthesis path): it
        # reconciles the rhizosphere shells with the host soil column and runs the
        # transpiration/uptake solve, populating co_hydr.ftc_*/btran/leaf_psi (the
        # fields the already-ported mortality + photosynthesis branches consume) and
        # bc_out.qflx_soil2root_sisl / plant_stored_h2o_si. Default path
        # (planthydro OFF) is byte-identical — the helper returns immediately.
        if hlm_use_planthydro[] == itrue
            fates_hydraulics_step!(inst; nlevsoil=varpar.nlevsoi, dtime=dtime)
        end
    end


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
    # Irrigation needed — WIRED (gated on config.irrigate)
    # ========================================================================
    # CalcIrrigationNeeded: at the end of the timestep, evaluate the soil-water
    # deficit on exposed irrigated patches and set the per-patch irrigation rate +
    # n_irrig_steps_left for the next step's withdrawal. local_time_sec_patch is the
    # local solar time (seconds since midnight) at each patch's longitude, mirroring
    # the Fortran get_local_time(londeg, offset=0). volr (river volume) is only used
    # when limit_irrigation_if_rof_enabled — passed as zeros (ROF limiting off).
    if config.irrigate && !isempty(irr.irrig_method_patch) && !isempty(pftcon.irrigated)
        local_time_sec_patch = zeros(Int, length(pch.gridcell))
        for p in bc_patch
            g = pch.gridcell[p]
            lon = grc.londeg[g]
            lon < 0.0 && (lon += 360.0)
            local_time_sec_patch[p] = mod(secs + round(Int, lon / DEGPSEC), ISECSPDAY)
        end
        volr_irr = zeros(eltype(ss.eff_porosity_col), length(grc.londeg))
        calc_irrigation_needed!(irr, cs.elai_patch, temp.t_soisno_col,
                                ss.eff_porosity_col, wsb.ws.h2osoi_liq_col,
                                volr_irr, rof_prognostic,
                                pftcon.irrigated, local_time_sec_patch,
                                col, grc, pch, filt.exposedvegp,
                                bc_col, bc_patch, bc_grc, varpar.nlevsoi)
    end

    # ========================================================================
    # EMISSIONS
    # ========================================================================
    # DustEmission — WIRED via the DustEmisFactory dispatcher (config.dust).
    # No-op unless config.dust.active (default false → byte-identical default
    # path). Fills inst.dust_emis.flx_mss_vrt_dst_patch from wind/soil state.
    dust_emission!(config.dust, inst.dust_emis, filt.nolakep, pch, lun,
                   bc_patch, bc.endl - bc.begl + 1,
                   a2l.forc_rho_downscaled_col, ss, cs, wsb.ws, wdb, fv)

    # DustDryDep — WIRED
    dust_dry_dep!(inst.dust_emis, pch.active, pch.column, bc_patch,
                  a2l.forc_pbot_downscaled_col, a2l.forc_rho_downscaled_col,
                  a2l.forc_t_downscaled_col,
                  fv.ram1_patch, fv.fv_patch)

    # VOCEmission — MEGAN VOC emissions (gated on use_voc + populated compounds).
    # Mirrors the Fortran `if (shr_megan_mechcomps_n < 1) return` gate: the call is
    # a no-op unless MEGAN compound/mechanism descriptors have been supplied on
    # `config.megan` (e.g. from a namelist read + megan_factors_get lookup).
    if config.use_voc && !isempty(config.megan.mech_comps) && !isempty(config.megan.meg_compounds)
        voc_emission!(
            inst.vocemis, config.megan.meg_compounds, config.megan.mech_comps, config.megan.megan_factors,
            pch, 1:length(filt.soilp), filt.soilp,
            a2l.forc_solad_downscaled_col, a2l.forc_solai_grc,
            a2l.forc_pbot_downscaled_col, a2l.forc_pco2_grc,
            a2l.fsd24_patch, a2l.fsd240_patch, a2l.fsi24_patch, a2l.fsi240_patch,
            cs.fsun_patch, cs.fsun24_patch, cs.fsun240_patch,
            cs.elai_patch, cs.elai240_patch,
            ps.cisun_z_patch, ps.cisha_z_patch,
            temp.t_veg_patch, temp.t_veg24_patch, temp.t_veg240_patch,
            ef.btran_patch)
    end

    # ========================================================================
    # DETERMINE TEMPERATURES
    # ========================================================================

    # LakeTemperature — WIRED
    grnd_ch4_cond = _zlike(nc)  # placeholder (CH4 module provides when active)
    lake_temperature!(col, pch, sa, ss, wsb, wdb, wfb, ef, temp, ls,
                      grnd_ch4_cond,
                      filt.lakec, filt.lakep,
                      bc_col, bc_patch, dtime)

    # SoilTemperature — WIRED
    # Placeholders for urban building-temperature time-varying inputs (urbantv).
    nl = length(lun.itype)
    urbantv_t_building_max = _fulllike(323.15, nl)  # placeholder ~50°C cap
    # AC adoption rate (0..1). Used only by the prognostic (CLM5) building path
    # with urban_explicit_ac=true; 1.0 = full adoption (placeholder until urbantv
    # streams are read).
    urbantv_p_ac = _fulllike(1.0, nl)
    soil_temperature!(col, lun, pch, temp, ef, ss, wsb, wdb, wfb, sa, cs, up,
                      urbantv_t_building_max,
                      a2l.forc_lwrad_downscaled_col,
                      filt.nolakec, filt.nolakep,
                      filt.urbanl, filt.urbanc,
                      bc_col, bc_lun, bc_patch,
                      dtime;
                      urbantv_p_ac = urbantv_p_ac)



    # Glacier ice melt: convert layer meltwater back to ice over do_smb istice
    # columns, accumulating qflx_glcice_melt. (GlacierSurfaceMassBalanceMod:HandleIceMelt)
    handle_ice_melt!(wsb.ws.h2osoi_liq_col, wsb.ws.h2osoi_ice_col,
                     wfb.wf.qflx_glcice_melt_col,
                     col.landunit, lun.itype,
                     filt.do_smb_c, bc_col, dtime,
                     varpar.nlevsno, varpar.nlevgrnd)

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



    # 2e'. Update int_snow, and RESET frac_sno / snow_depth / int_snow when no snow
    # present (Fortran BulkDiag_SnowWaterAccumulatedSnow, called at the end of
    # SnowWater in SnowHydrologyMod.F90:1135). Without this, a no-layer (trace) snow
    # event that fully melts leaves frac_sno / snow_depth nonzero forever → residual
    # snow-cover fraction on zero-snow days (Tagus FSNO annual mean ~5× Fortran; ~50%
    # of it from these stale-frac_sno days). Uses the same snowc/nosnowc filters.
    bulkdiag_snow_water_accumulated_snow!(
        wsb.int_snow_col,
        wdb.frac_sno_col,
        wdb.snow_depth_col,
        dtime,
        wdb.frac_sno_eff_col,
        wfb.wf.qflx_soliddew_to_top_layer_col,
        wfb.wf.qflx_liqdew_to_top_layer_col,
        wfb.wf.qflx_liq_grnd_col,
        wsb.ws.h2osno_no_layers_col,
        filt.snowc, filt.nosnowc, bc_col)



    # 2f. Sum percolation out of bottom snow layer → qflx_rain_plus_snomelt
    # For snow columns: drainage from bottom + non-snow-covered rain
    # For no-snow columns: rain + snowmelt
    # Extract bottom-layer percolation (layer 0 in Fortran = index nlevsno in Julia)
    qflx_snow_perc_bottom = _zlike(nc)
    # Masked bottom-layer gather as a device-safe broadcast (snowc ? perc : 0; the
    # array is zero-initialized so `Bool * val` is byte-identical to the masked set).
    @views qflx_snow_perc_bottom[bc_col] .=
        filt.snowc[bc_col] .* wfb.wf.qflx_snow_percolation_col[bc_col, nlevsno]
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
    # Pass the qinmax populated by infiltration_excess_runoff! (step 6) so ponded
    # surface water can infiltrate, matching Fortran UpdateH2osfc (uses
    # infiltration_excess_runoff_inst%qinmax_col).
    update_h2osfc!(col, sh, ef, wfb, wsb, wdb,
                   filt.hydrologyc, bc_col; dtime=dtime,
                   qinmax=inst.infilt_excess_runoff.qinmax_col)

    # --- 9. Infiltration ---
    infiltration!(wfb, filt.hydrologyc, bc_col)

    # --- 10. Total surface runoff ---
    total_surface_runoff!(
        wfb, sh, wsb.ws,
        col.snl, col.itype, col.landunit, lun.urbpoi,
        filt.hydrologyc, filt.urbanc,
        bc_col, dtime)

    # --- 10b. Urban ponding state (roof / impervious road) ---
    # Fortran calls UpdateUrbanPonding IMMEDIATELY after TotalSurfaceRunoff
    # (HydrologyNoDrainageMod.F90:336). `update_urban_ponding!` was PORTED and
    # never called — another dead port, the same class as the *_init_cold! sweep.
    #
    # Roof and impervious-road columns are NOT hydrologically active: they have no
    # infiltration and their h2osfc is excluded from the column water mass
    # (TotalWaterAndHeatMod: "Nothing more to add in this case"). Their ONLY water
    # store is the pond held in the top soil layer, and UpdateUrbanPonding is the
    # ONLY routine that advances it:
    #     h2osoi_liq(c,1) += (qflx_rain_plus_snomelt - qflx_liqevap_from_top_layer)*dt
    # (capped at pondmx_urban, floored at 0).
    #
    # Without it, rain landing on a roof was neither STORED nor RUN OFF, and the
    # evaporation charged against those columns was never DEBITED from any store —
    # so `endwb - begwb` could not match the flux sum. That is the 0.077 mm
    # impervious-road water-balance error, which was invisible only because
    # `errh2o` was NaN there (waterflux InitCold was dead) and `abs(NaN) > tol` is
    # false, so the check silently passed.
    update_urban_ponding!(
        wsb.ws, sh, wfb,
        col.snl, col.itype, filt.urbanc,
        bc_col, dtime)

    # --- 11. Root water uptake (transpiration sink) ---
    # Fortran SoilWaterPlantSinkMod reads use_hydrstress from clm_varctl and, when
    # true, distributes the per-layer sink from k_soil_root × (smp − vegwp − grav2)
    # (the HydStress branch) rather than the default rootr×qflx_tran profile. The
    # flag MUST be threaded through here in PHS mode, else qflx_rootsoi_col collapses
    # to the static root profile (wrong per-layer H2OSOI_LIQ). use_hydrstress=false
    # ≡ omitting it, so the non-PHS path stays byte-identical.
    compute_effec_rootfrac_and_vert_tran_sink!(
        bc_col, nlevsoi,
        filt.hydrologyc,
        ss, cs, wfb, ef,
        col, lun, pch;
        use_hydrstress=config.use_hydrstress)

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
    # Fortran HydrologyNoDrainageMod.F90:356 branches on use_aquifer_layer:
    #   true  → WaterTable (Zeng-Decker aquifer; RECOMPUTES frost_table/zwt_perched)
    #   false → ThetaBasedWaterTable (zwt from the moisture profile; does NOT touch
    #           frost_table — that stays at its restart value, and the perched
    #           drainage is handled by perched_lateral_flow! in HydrologyDrainage).
    # Calling the aquifer water_table! unconditionally wrongly recomputed
    # frost_table = col_z[k_frz] (node) every step, where Fortran leaves it at the
    # restart value; that inflated zwt_perched and qflx_drain_perched.
    if config.use_aquifer_layer
        water_table!(sh, ss,
                     temp.t_soisno_col, wsb.ws, wfb,
                     col.dz, col.z, col.zi,
                     filt.hydrologyc, bc_col,
                     nlevsoi, dtime;
                     recompute_frost_table=true)
    else
        # CLM5 default (HydrologyNoDrainageMod.F90:360-368): PerchedWaterTable +
        # ThetaBasedWaterTable diagnose the (perched) water table from the moisture
        # profile each step; zwt is NOT advanced by aquifer qcharge. Without this,
        # zwt stays pinned at its restart value, so the topographic baseflow
        # (zi_bedrock - zwt)^n_baseflow in subsurface_lateral_flow! → 0 and the soil
        # never drains — saturating the surface, collapsing the dry-surface layer,
        # and blowing up ground evaporation (Aripuana FGEV +425%).
        perched_water_table!(sh, ss, temp.t_soisno_col, wsb.ws,
                             col.dz, col.z, col.zi,
                             filt.hydrologyc, bc_col, nlevsoi)
        theta_based_water_table!(sh, ss, wsb.ws,
                             col.dz, col.z, col.zi, col.nbedrock,
                             filt.hydrologyc, bc_col, nlevsoi)
    end

    # Bedrock clipping: Fortran ThetaBasedWaterTable clips ZWT at bedrock depth
    # when soil above bedrock is not saturated. Without this, ZWT drops to 80m
    # in cold-start runs because the aquifer drains without replenishment.
    nlevsno_l = varpar.nlevsno
    clamp_zwt_to_bedrock!(sh.zwt_col, filt.hydrologyc, col.nbedrock, col.zi, nlevsno_l)

    # --- 14. Condensation renewal ---
    renew_condensation!(
        wsb.ws, wdb, wfb,
        col.snl, col.itype,
        filt.hydrologyc, filt.urbanc,
        bc_col, dtime)

    # --- 15. Snow capping excess ---
    h2osno_total = _zlike(nc)
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
    mask_capping = _blike(nc)
    rho_orig_bottom = _zlike(nc)
    frac_adjust = _zlike(nc)
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
    # Write back bottom layer state (device-safe broadcast; non-capping entries are
    # the unchanged extracted values, so writing all back is byte-identical).
    @views wsb.ws.h2osoi_ice_col[bc_col, jj_bottom] .= ice_bottom_vec[bc_col]
    @views wsb.ws.h2osoi_liq_col[bc_col, jj_bottom] .= liq_bottom_vec[bc_col]
    # Update dz and aerosols after capping
    snow_capping_update_dz_and_aerosols!(
        dz_bottom_vec, aer, jj_bottom,
        rho_orig_bottom, ice_bottom_vec, frac_adjust,
        mask_capping, bc_col)
    # Device-safe broadcast write-back (non-capping entries unchanged).
    @views col.dz[bc_col, jj_bottom] .= dz_bottom_vec[bc_col]

    # --- 16. Snow compaction ---
    # The melt-compaction term needs the SwensonLawrence2012 SCA shape parameter
    # (n_melt) and int_snow_max to compute fsno_melt = FracSnowDuringMelt. Extract
    # them from the active SCF parameterization (NiuYang2007 has neither → fall back
    # to defaults that never take the SwensonLawrence branch materially).
    _scf = inst.scf_method
    _n_melt = hasproperty(_scf, :n_melt) && length(_scf.n_melt) == length(col.snl) ?
        _scf.n_melt : fill(20.0, length(col.snl))
    _int_snow_max = hasproperty(_scf, :int_snow_max) ? _scf.int_snow_max : 2000.0
    snow_compaction!(
        col.dz, dtime, col.snl,
        temp.t_soisno_col,
        wsb.ws.h2osoi_ice_col, wsb.ws.h2osoi_liq_col,
        temp.imelt_col,
        wdb.frac_sno_col, wdb.frac_h2osfc_col,
        wdb.swe_old_col, wsb.int_snow_col, wdb.frac_iceold_col,
        a2l.forc_wind_grc, col.gridcell, col.landunit,
        lun.lakpoi, lun.urbpoi,
        filt.snowc, bc_col, nlevsno;
        use_subgrid_fluxes = varctl.use_subgrid_fluxes,
        n_melt = _n_melt, int_snow_max = _int_snow_max)



    # --- 17. Combine thin snow layers ---
    combine_snow_layers!(
        col.snl, col.dz, col.zi, col.z,
        temp.t_soisno_col,
        wsb.ws.h2osoi_ice_col, wsb.ws.h2osoi_liq_col,
        wsb.ws.h2osno_no_layers_col,
        wdb.snow_depth_col, wdb.frac_sno_col, wdb.frac_sno_eff_col,
        wsb.int_snow_col, wdb.snw_rds_col,
        aer, lun.itype, lun.urbpoi, col.landunit,
        filt.snowc, bc_col, nlevsno,
        wfb.wf.qflx_sl_top_soil_col, dtime)
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

    # Glacier surface mass balance: ice growth (frz), net glcice, dyn water-flux
    # balance term. (GlacierSurfaceMassBalanceMod:ComputeSurfaceMassBalance)
    # glc_dyn_runoff_routing is not yet ported (no glc2lnd) → 0 (standalone CLM).
    let ng_smb = isempty(bc_col) ? 0 : maximum(Array(col.gridcell)[bc_col])
        glc_dyn_runoff_routing_grc =
            fill!(similar(wfb.wf.qflx_glcice_col, FT, ng_smb), zero(FT))
        compute_surface_mass_balance!(
            wfb.wf.qflx_glcice_col,
            wfb.wf.qflx_glcice_frz_col,
            wfb.wf.qflx_glcice_dyn_water_flux_col,
            wfb.wf.qflx_snwcp_ice_col,
            wfb.wf.qflx_glcice_melt_col,
            wsb.snow_persistence_col,
            glc_dyn_runoff_routing_grc,
            col.landunit, col.gridcell, lun.itype,
            filt.allc, filt.do_smb_c, bc_col)  # default glc_snow_persistence_max_days=7300
    end

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
    update_urban_frac_sno!(wdb.frac_sno_col, col.landunit, lun.urbpoi, wdb.snow_depth_col)

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
    # CN fire (Li family) active? Needs the method selected AND the fire bundle
    # sized by clm_instInit!(; cnfire_method=…) — see the fire-bundle kwargs below.
    _fire_on = config.cnfire_method !== :nofire &&
               !isempty(inst.cnfire_li2014.forc_lnfm) &&
               !isempty(inst.cnfire_base.btran2_patch)
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
            # C13/C14 soil-BGC state/flux for the CIsoFlux* cascade — passed only
            # when the tracer is on (else nothing → cascade skipped). The facade
            # additionally guards on the parallel veg + soil state being allocated.
            c13_soilbgc_cs = config.use_c13 ? inst.c13_soilbiogeochem_carbonstate : nothing,
            c13_soilbgc_cf = config.use_c13 ? inst.c13_soilbiogeochem_carbonflux : nothing,
            c14_soilbgc_cs = config.use_c14 ? inst.c14_soilbiogeochem_carbonstate : nothing,
            c14_soilbgc_cf = config.use_c14 ? inst.c14_soilbiogeochem_carbonflux : nothing,
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
            # Vegetation-flux inputs for the veg-CN flux chain (mresp, gpp/alloc, …)
            patch=pch,
            pftcon_main=pftcon,
            crop=inst.crop,
            photosyns=inst.photosyns,
            canopystate=inst.canopystate,
            soilstate=ss,
            temperature=temp,
            water_diag=wdb,
            gridcell=grc,
            is_first_step=is_first_step,
            # Accumulated forcing (Fortran wateratm2lndbulk_type). Consumed by the
            # Li fire schemes and by CNPhenology's rain-triggered stress-decid onset.
            # These are only allocated under use_cn (as in Fortran's InitAccBuffer).
            prec10_patch=a2l.prec10_patch,
            prec30_patch=a2l.prec30_patch,
            prec60_patch=a2l.prec60_patch,
            rh30_patch=a2l.rh30_patch,
            forc_rh_grc=a2l.forc_rh_grc,
            forc_wind_grc=a2l.forc_wind_grc,
            # ---- Mineral-N inputs (deposition + fixation) ----
            # forc_ndep_grc is filled by ndep_interp! above (zero when no stream is
            # configured, which reproduces the historical no-N-input behaviour).
            # nfix_timeconst comes from clm_varctl (control.jl derives it as 10 d
            # when use_nitrif_denitrif, else 0) and selects CNNFixation's lagged-NPP
            # vs annual-mean-NPP branch. AnnET drives the FUN free-living fixation
            # (a 365-day running mean, now maintained by the accumulator subsystem).
            forc_ndep=a2l.forc_ndep_grc,
            AnnET=wfb.AnnET,
            nfix_timeconst=varctl.nfix_timeconst,
            dayspyr=365.0,
            h2osoi_vol=(_decomp_initialized(inst.decomp_cascade) ? wsb.ws.h2osoi_vol_col : nothing),
            h2osoi_liq=(_decomp_initialized(inst.decomp_cascade) ? wsb.ws.h2osoi_liq_col[:, (varpar.nlevsno+1):end] : nothing),
            # ---- CN fire (Li family) bundle ----
            # Supplied ONLY when config.cnfire_method !== :nofire; otherwise every
            # entry is `nothing`/empty and cn_driver's `_fire_active` gate is false,
            # so the default path is byte-identical. `_fire_on` also guards on the
            # bundle actually having been allocated by clm_instInit! (a caller that
            # asked for fire but never initialized it gets no fire rather than an
            # out-of-bounds read on a zero-length forc_lnfm).
            fire_data          = _fire_on ? inst.cnfire_base : nothing,
            fire_li2014        = _fire_on ? inst.cnfire_li2014 : nothing,
            cnfire_const       = _fire_on ? inst.cnfire_const : nothing,
            cnfire_params      = _fire_on ? inst.cnfire_params : nothing,
            pftcon_fire        = _fire_on ? inst.pftcon_fire : nothing,
            pftcon_fire_li2014 = _fire_on ? inst.pftcon_fire_li2014 : nothing,
            dgvs_fire          = _fire_on ? inst.dgvs_fire : nothing,
            mask_exposedvegp_fire   = _fire_on ? filt.exposedvegp : nothing,
            mask_noexposedvegp_fire = _fire_on ? filt.noexposedvegp : nothing,
            # Fire forcing + water-state inputs (all column-level).
            forc_t_fire_col    = _fire_on ? a2l.forc_t_downscaled_col : Float64[],
            forc_rain_fire_col = _fire_on ? forc_rain_col : Float64[],
            forc_snow_fire_col = _fire_on ? forc_snow_col : Float64[],
            fsat_fire_col      = _fire_on ? inst.sat_excess_runoff.fsat_col : Float64[],
            wf_fire_col        = _fire_on ? wdb.wf_col : Float64[],
            wf2_fire_col       = _fire_on ? wdb.wf2_col : Float64[],
            fire_kmo = mon, fire_kda = day, fire_mcsec = secs, fire_nstep = nstep,
            nlevgrnd_fire = varpar.nlevgrnd,
            # CTSM: `transient_landcover = run_has_transient_landcover()`, which is
            # do_transient_pfts .OR. do_transient_crops .OR. do_transient_urban —
            # NOT "a dyn_subgrid state exists". A DynSubgridState is constructed with
            # ALL transient aspects off (a near-identity reweight), and one built for
            # do_transient_LAKES only must still report false (lakes are deliberately
            # excluded from the Fortran predicate). Testing `!== nothing` turned both
            # of those into `true`, which switches CNFireFluxes' per-patch burned
            # fraction from `farea_burned` to `fbac` and enables the lfc/lfc2
            # deforestation bookkeeping in runs where CTSM leaves them off.
            transient_landcover = (config.dyn_subgrid !== nothing &&
                                   run_has_transient_landcover(config.dyn_subgrid.ctl)),
            mask_actfirec=filt.actfirec,
            mask_actfirep=filt.actfirep)
    end

    # CN vegetation structure update (leafc -> tlai/tsai/htop/elai) — WIRED.
    # Runs on the radiation step after the CN dynamics/phenology so leaf-out (leafc
    # growth from onset) is reflected in LAI for the next step's albedo + photosynthesis.
    # Mirrors CTSM's CNVegStructUpdate call in the EcosystemDynamics path.
    if config.use_cn && !config.use_fates && doalb
        cn_veg_struct_update!(filt.bgc_vegp, bc_patch, pch, cs,
            inst.bgc_vegetation.cnveg_carbonstate_inst, wdb, inst.frictionvel,
            inst.bgc_vegetation.cnveg_state_inst, inst.crop, pftcon;
            dt=dtime, npcropmin=config.npcropmin)
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
    clamp_zwt_to_bedrock!(sh.zwt_col, filt.hydrologyc, col.nbedrock, col.zi, varpar.nlevsno)

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
            soilbgc_nf=inst.soilbiogeochem_nitrogenflux,
            patch_itype=pch.itype,
            # NB: the SoilBiogeochem carbon-flux Summary (somhr/lithr/hr_vr — CH4's
            # substrate) is wired by #223, which supplies `decomp`/`dzsoi_decomp_vals`
            # further down this same call. It was independently found dead by both that
            # work and this one.
            # ---- N leaching (the mineral-N LOSS term) ----
            # Post-drainage, so qflx_drain/qflx_surf are this step's final values —
            # which is exactly why Fortran calls SoilBiogeochemNLeaching here and not
            # in the pre-drainage half. Soil-layer slices (drop the snow levels).
            h2osoi_liq=(_decomp_initialized(inst.decomp_cascade) ?
                        wsb.ws.h2osoi_liq_col[:, (varpar.nlevsno+1):end] : nothing),
            qflx_drain=(_decomp_initialized(inst.decomp_cascade) ? wfb.wf.qflx_drain_col : nothing),
            qflx_surf=(_decomp_initialized(inst.decomp_cascade) ? wfb.wf.qflx_surf_col : nothing),
            col_dz=(_decomp_initialized(inst.decomp_cascade) ?
                    col.dz[:, (varpar.nlevsno+1):end] : nothing),
            zisoi=(_decomp_initialized(inst.decomp_cascade) ? zisoi[] : nothing),
            nlevsoi=(_decomp_initialized(inst.decomp_cascade) ? varpar.nlevsoi : 0),
            # Column NPP + lagged NPP (consumed by NEXT step's n_fixation!).
            col=col,
            patch=pch,
            nfix_timeconst=varctl.nfix_timeconst,
            # CNDriverSummarizeFluxes: the soil-BGC C/N flux summaries and the
            # column/gridcell halves of the CNVeg C/N flux summaries. These produce
            # hr_col / som_c_leached_col / denit_col / f_n2o_nit_col / gpp_col /
            # er_col / nbp_grc — i.e. EVERY term the CN balance check reads. Before
            # this they were structurally dead (see cnveg_carbon_flux_summary_col!).
            bounds_grc=bc_grc,
            decomp=(_decomp_initialized(inst.decomp_cascade) ? inst.decomp_cascade : nothing),
            dzsoi_decomp_vals=(_decomp_initialized(inst.decomp_cascade) ? dzsoi_decomp[] : nothing))

        # C14 radioactive decay — WIRED (gated on config.use_c14). Mirrors the
        # Fortran C14Decay call in CNDriver: after the C/N state updates, decay
        # every C14 pool (gridcell seedc / soil decomp vr / veg patch pools). Runs
        # only when the parallel C14 state is stood up (the facade C14 veg state +
        # inst C14 soil state, both sized by use_c14); the isempty guards keep it a
        # no-op if the flag is set without the state allocated. Default off.
        let _c14cs  = inst.bgc_vegetation.c14_cnveg_carbonstate_inst,
            _c14scs = inst.c14_soilbiogeochem_carbonstate
            if config.use_c14 && !isempty(_c14cs.cpool_patch) &&
               !isempty(_c14scs.decomp_cpools_vr_col)
                c14_decay!(_c14cs, inst.bgc_vegetation.c14_cnveg_carbonflux_inst,
                    _c14scs, inst.c14_soilbiogeochem_carbonflux;
                    mask_soilc=filt.soilc, mask_soilp=filt.soilp,
                    bounds_col=bc_col, bounds_patch=bc_patch,
                    bounds_gridcell=bc.begg:bc.endg,
                    dt=dtime, nlevdecomp=varpar.nlevdecomp,
                    ndecomp_pools=config.ndecomp_pools,
                    col_gridcell=inst.column.gridcell, latdeg_grc=inst.gridcell.latdeg,
                    ivt=pch.itype, npcropmin=config.npcropmin,
                    use_matrixcn=config.use_matrixcn)
            end
        end

        # Soil-biogeochem C/N STATE summaries (totsomc_col/totsomn_col/totecosysc_col)
        # — WIRED. cn_driver_summarize_states! only does the veg-state summary; the soil
        # summaries need the decomp infrastructure, which is in scope here. Pure
        # diagnostics (no downstream physics reads them), so this is an output-correctness
        # fix. totvegc_col/totvegn_col (for the ecosystem totals) are aggregated from the
        # veg patch totals via p2c_1d! (unity scale).
        if config.use_cn && _decomp_initialized(inst.decomp_cascade)
            _scs = inst.soilbiogeochem_carbonstate
            _sns = inst.soilbiogeochem_nitrogenstate
            _dc  = inst.decomp_cascade
            _ccs = inst.bgc_vegetation.cnveg_carbonstate_inst
            _cns = inst.bgc_vegetation.cnveg_nitrogenstate_inst
            # TWO DIFFERENT veg aggregates — Fortran does not use one for both
            # (SoilBiogeochemCarbonStateType.F90:1627-1644):
            #
            #   totecosysc_col (TOTECOSYSC) += totvegc_col  <- p2c of totvegc_patch,
            #                                                  EXCLUDES cpool
            #   totc_col       (TOTCOLC)    += totc_p2c_col <- p2c of totc_patch,
            #                                                  INCLUDES cpool/xsmrpool/ctrunc
            #
            # `totc_col` is the pool the CN balance check integrates (begcb/endcb), so
            # it MUST be the cpool-INCLUSIVE one. The port passed the p2c of
            # totvegc_patch for both and never supplied totc_p2c_col at all — it
            # defaulted to zeros(), so totc_col carried *no vegetation carbon whatsoever*.
            #
            # GPP lands in cpool before it is allocated to the tissue pools. With cpool
            # outside the budget, every gC of GPP looked like carbon created from
            # nothing: the column C imbalance was ~gpp*dt (0.045 gC/m2/step at Bow in
            # summer, vs a cerror of 1e-7). This is the bug the carbon balance check
            # was built to catch, and it could not run to catch it because gpp_col/er_col
            # were themselves dead. Same story verbatim on the N side (totn_p2c_col).
            _tvc  = zeros(bc.endc); _tvn  = zeros(bc.endc)
            _tp2c = zeros(bc.endc); _tn2c = zeros(bc.endc)
            p2c_1d!(_tvc,  _ccs.totvegc_patch, bc, "unity", pch)
            p2c_1d!(_tvn,  _cns.totvegn_patch, bc, "unity", pch)
            p2c_1d!(_tp2c, _ccs.totc_patch,    bc, "unity", pch)
            p2c_1d!(_tn2c, _cns.totn_patch,    bc, "unity", pch)
            # Keep the canonical fields in sync (they are history/restart outputs).
            if !isempty(_ccs.totc_p2c_col); copyto!(_ccs.totc_p2c_col, _tp2c); end
            if !isempty(_cns.totn_p2c_col); copyto!(_cns.totn_p2c_col, _tn2c); end
            soil_bgc_carbon_state_summary!(_scs, filt.allc, bc_col;
                nlevdecomp=varpar.nlevdecomp, ndecomp_pools=config.ndecomp_pools,
                dzsoi_decomp_vals=dzsoi_decomp[], zisoi_vals=zisoi[],
                is_litter=_dc.is_litter, is_soil=_dc.is_soil, is_cwd=_dc.is_cwd,
                totvegc_col=_tvc, totc_p2c_col=_tp2c)
            soil_bgc_nitrogen_state_summary!(_sns, filt.allc, bc_col;
                nlevdecomp=varpar.nlevdecomp, ndecomp_pools=config.ndecomp_pools,
                dzsoi_decomp_vals=dzsoi_decomp[], zisoi_vals=zisoi[],
                is_litter=_dc.is_litter, is_soil=_dc.is_soil, is_cwd=_dc.is_cwd,
                totvegn_col=_tvn, totn_p2c_col=_tn2c)

            # SPECIAL (non-soil) columns carry NO carbon/nitrogen budget. Fortran zeroes
            # their C/N state outright — SoilBiogeochemCarbonStateType.F90:686
            # (`this%totc_col(c) = 0._r8` for the special-column filter). The port left
            # them at the NaN allocation default, and the gridcell balance then computed
            #     totc_grc += totc_col[c] * wtgcell[c]   ==   NaN * 0.0   ==   NaN
            # — the NaN*0 trap. Every gridcell C and N imbalance was NaN, so the
            # GRIDCELL half of the CN balance check could never fail no matter how badly
            # it was violated (NaN > cerror is false). Zeroing them, as Fortran does,
            # makes the gridcell check able to fail.
            for c in bc_col
                filt.soilc[c] && continue
                _scs.totc_col[c] = 0.0
                _scs.totecosysc_col[c] = 0.0
                _sns.totn_col[c] = 0.0
                _sns.totecosysn_col[c] = 0.0
            end
        end
    end

    # ========================================================================
    # FATES dynamics
    # ========================================================================
    if config.use_fates
        # Placeholder: clm_fates%WrapUpdateFatesRmean(...) [FATES running mean]
        # Placeholder: clm_fates%wrap_update_hifrq_hist(...) [FATES history]

        # FATES parity instrumentation — "fast" (sub-daily) phase. Fires on EVERY
        # timestep at exactly the point the Fortran clm_driver does (after the
        # photosynthesis/flux solve, the FATES running-mean update and the hifrq
        # history fill; BEFORE the daily dynamics branch), so the fast thread can be
        # scored on its own. No-op unless a harness installs FATES_PARITY_HOOK.
        if inst.fates !== nothing && FATES_PARITY_HOOK[] !== nothing
            for s in 1:inst.fates.nsites
                fates_parity_hook(inst.fates.sites[s], inst.fates.bc_in[s], 1)
            end
        end

        if is_beg_curr_day && inst.fates !== nothing
            # SPITFIRE fire-weather (gated): when SPITFIRE is on, derive the per-day
            # representative precip / relative-humidity / wind from the FATES column's
            # current downscaled atmospheric forcing and pass them to the daily step
            # (the bc daily pack fills the per-patch fire-weather bc_in the fire model
            # reads). A single-step stand-in for the true 24-hour running means (the
            # bc-pack docstring notes the running-mean nuance); finite + representative.
            # OFF (no_fire) => fire_weather=nothing and the fields stay zero (the fire
            # chain never runs), so the default path is byte-identical.
            fire_weather = nothing
            if hlm_spitfire_mode[] > hlm_sf_nofire_def[]
                fc = _fates_first_fates_column(col)
                if fc != 0
                    g = col.gridcell[fc]
                    a2l = inst.atm2lnd
                    rain = a2l.forc_rain_downscaled_col[fc]
                    snow = a2l.forc_snow_downscaled_col[fc]
                    # Relative humidity [%] from specific humidity / saturation at
                    # forcing temperature & pressure (Tetens; clamped to [1,100]).
                    qair = a2l.forc_q_downscaled_col[fc]
                    tair = a2l.forc_t_downscaled_col[fc]
                    pbot = a2l.forc_pbot_downscaled_col[fc]
                    esat_fw = 610.78 * exp(17.27 * (tair - 273.15) / (tair - 35.86))  # [Pa]
                    qsat_fw = 0.622 * esat_fw / max(pbot - 0.378 * esat_fw, 1.0)
                    rh = clamp(100.0 * qair / max(qsat_fw, 1.0e-9), 1.0, 100.0)
                    # Wind speed from the east/north components (forc_wind_grc may be
                    # an unset forcing slot; u/v are always packed by the harness/driver).
                    wind = sqrt(a2l.forc_u_grc[g]^2 + a2l.forc_v_grc[g]^2)
                    fire_weather = (precip24 = rain + snow, relhumid24 = rh,
                                    wind24 = wind)
                end
            end

            # W5 daily demographic step — FATES `dynamics_driv`. Runs once per FATES
            # site at the day boundary: pack the daily soil/decomp bc_in (+ fire
            # weather when SPITFIRE is on), advance the cohort/patch population
            # (ed_ecosystem_dynamics — includes the fire model, CNP nutrient
            # acquisition, and tree-damage when those modes are enabled), recompute
            # site diagnostics + canopy structure (ed_update_site), run the final
            # TotalBalanceCheck mass-conservation audit, and unpack the new canopy
            # structure (elai/htop/...) back to canopystate. The daily step walks
            # ALL the column's vegetated patches (p = ifp + col.patchi[c]) to unpack
            # canopy structure, and rebuilds the HLM per-patch weights
            # (fates_set_filters!, the setFilters-equivalent) from the freshly-advanced
            # FATES patch areas — so a disturbance that splits/merges patches
            # propagates into the per-column averaging. Fire weather is packed when
            # SPITFIRE is on.
            fates_daily_dynamics_step!(inst; nlevsoil=varpar.nlevsoi,
                                       nlevdecomp=varpar.nlevdecomp,
                                       use_fates_bgc=config.use_fates_bgc,
                                       fire_weather=fire_weather)
        end
    end

    # ========================================================================
    # Water diagnostic summaries
    # ========================================================================
    h2osno_total_col = _zlike(nc)
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
                   urbpoi=lun.urbpoi,
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
            soilbgc_ns=inst.soilbiogeochem_nitrogenstate,
            col=col, grc=grc, bounds_grc=bc_grc, dt=dtime,
            # The CN mass-conservation check. OFF by default (it error()s on
            # failure and the C side has known residuals from the not-yet-wired
            # product pools); flip on with `CLM.cn_balance_check_enabled!(true)`
            # or per-instance via the facade config. See docs/N_CYCLE_PARITY.md.
            run_check=cn_balance_check_enabled())
    end

    # ========================================================================
    # Methane fluxes — WIRED (gated on config.use_lch4)
    # ========================================================================
    # ch4! consumes the just-computed BGC respiration (somhr/lithr/hr_vr/o_scalar/
    # fphr/pot_f_nit_vr) + soil temp/moisture + vegetation roots/NPP and produces the
    # column CH4 surface flux + concentration profiles. The soil-physics matrices
    # col.dz/z/zi and temp/water carry an nlevsno snow-layer prefix, so they are
    # soil-sliced to [c, 1:nlevsoi]; the BGC fluxes are nlevdecomp_full-resolved and
    # sliced to the top nlevsoi layers. Gridcell forcings (forc_t/pbot/po2/pco2/pch4)
    # are read by gridcell index inside ch4!.
    if config.use_lch4 && !isempty(inst.ch4.ch4_surf_flux_tot_col)
        _sl = (nlevsno + 1):(nlevsno + nlevsoi)         # snow-offset soil slice
        _dl = 1:nlevsoi                                  # decomp-resolved soil slice
        cnf_ch4 = inst.soilbiogeochem_carbonflux
        nnf_ch4 = inst.soilbiogeochem_nitrogenflux
        cveg_cf = inst.bgc_vegetation.cnveg_carbonflux_inst
        # lake_icefrac_col is nlevlak-resolved (nlevlak may be < nlevsoi); pad to the
        # nlevsoi soil grid ch4! expects (lake path only reads lake columns, deeper
        # layers default frozen-free). Replaces NaN cold-start with 0.
        _ch4_lakeicefr = zeros(eltype(temp.t_soisno_col), size(ls.lake_icefrac_col, 1), nlevsoi)
        let _nlk = min(size(ls.lake_icefrac_col, 2), nlevsoi)
            @views _ch4_lakeicefr[:, 1:_nlk] .= ifelse.(isnan.(ls.lake_icefrac_col[:, 1:_nlk]),
                                                        zero(eltype(_ch4_lakeicefr)),
                                                        ls.lake_icefrac_col[:, 1:_nlk])
        end
        ch4!(inst.ch4, inst.ch4_params, inst.ch4_varcon,
             filt.soilc, filt.soilp, filt.lakec, filt.nolakec,
             Array(col.gridcell), Array(col.wtgcell),
             Array(pch.column), Array(pch.itype), Array(pch.wtcol),
             Array(col.is_fates), grc.latdeg,
             a2l.forc_pbot_not_downscaled_grc, a2l.forc_t_not_downscaled_grc,
             a2l.forc_po2_grc, a2l.forc_pco2_grc, a2l.forc_pch4_grc,
             Matrix(ss.watsat_col[:, _dl]),
             Matrix(wsb.ws.h2osoi_vol_col[:, _dl]),
             Matrix(wsb.ws.h2osoi_liq_col[:, _sl]),
             Matrix(wsb.ws.h2osoi_ice_col[:, _sl]),
             wsb.ws.h2osfc_col,
             Matrix(ss.bsw_col[:, _dl]),
             Matrix(ss.cellorg_col[:, _dl]),
             Matrix(ss.smp_l_col[:, _dl]),
             Matrix(temp.t_soisno_col[:, _sl]),
             temp.t_grnd_col, temp.t_h2osfc_col, wdb.frac_h2osfc_col,
             wdb.snow_depth_col, Array(col.snl), wfb.wf.qflx_surf_col,
             Matrix(ss.rootfr_patch[:, _dl]), Matrix(ss.rootfr_col[:, _dl]),
             Matrix(ss.crootfr_patch[:, _dl]), Matrix(ss.rootr_patch[:, _dl]),
             cs.elai_patch, wfb.wf.qflx_tran_veg_patch,
             cveg_cf.annsum_npp_patch, cveg_cf.rr_patch,
             cnf_ch4.somhr_col, cnf_ch4.lithr_col,
             Matrix(cnf_ch4.hr_vr_col[:, _dl]), Matrix(cnf_ch4.o_scalar_col[:, _dl]),
             Matrix(cnf_ch4.fphr_col[:, _dl]), Matrix(nnf_ch4.pot_f_nit_vr_col[:, _dl]),
             _ch4_lakeicefr, col.lakedepth,
             Matrix(col.z[:, _sl]), Matrix(col.dz[:, _sl]),
             Matrix(col.zi[:, (nlevsno + 1):(nlevsno + nlevsoi + 1)]),
             nlevsoi, nlevsno, varpar.nlevdecomp, varpar.nlevdecomp_full,
             5, 0.01, 130.0,
             bc.endg, noveg, dtime,
             config.use_cn, false, false,
             cveg_cf.agnpp_patch, cveg_cf.bgnpp_patch, 365.0 * SECSPDAY;
             finundated_stream = (inst.ch4_varcon.finundation_mtd == FINUNDATION_MTD_H2OSFC ?
                                  nothing : inst.ch4_finundated_stream),
             zwt = sh.zwt_col, zwt_perched = sh.zwt_perched_col,
             tws = inst.water.waterdiagnosticbulk_inst.tws_grc)
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
                        _onbk(pftcon.rhol), _onbk(pftcon.rhos), _onbk(pftcon.taul), _onbk(pftcon.taus), _onbk(pftcon.xl),
                        coszen_cday)

        # UrbanAlbedo — WIRED
        urban_albedo!(filt.urbanl, filt.urbanc, filt.urbanp,
                      lun, col, pch, wsb, wdb, up, sa, alb)

        # W3b canopy radiation — FATES.  After surface_albedo! has computed the
        # ground albedo (albgrd/albgri) + coszen for each column, pack the FATES
        # radiation bc_in and run FatesNormalizedCanopyRadiation to fill the
        # per-patch normalized albedo / absorbed / transmitted fractions, then
        # unpack albd/albi/fab*/ft* into surfalb for the next radiation step.
        if config.use_fates && inst.fates !== nothing
            s = 0
            for c in 1:length(col.is_fates)
                col.is_fates[c] || continue
                s += 1
                s <= inst.fates.nsites || break
                for (ifp, p) in fates_veg_patches(inst.fates.sites[s], c, col)
                    fates_pack_bcin_radiation!(inst; s=s, c=c, p=p, ifp=ifp,
                                               coszen=alb.coszen_col[c])
                end
            end
            FatesNormalizedCanopyRadiation(inst.fates.nsites, inst.fates.sites,
                                           inst.fates.bc_in, inst.fates.bc_out)
            s = 0
            for c in 1:length(col.is_fates)
                col.is_fates[c] || continue
                s += 1
                s <= inst.fates.nsites || break
                for (ifp, p) in fates_veg_patches(inst.fates.sites[s], c, col)
                    fates_unpack_bcout_canopy_radiation!(inst; s=s, c=c, p=p, ifp=ifp)
                end
            end
        end
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
    # Placeholder: lnd2atm!(bounds_proc, ...) [gridcell averaging for the coupler]
    #
    # The ice-runoff part of lnd2atm IS wired: qflx_ice_runoff_col (= snow-capping
    # solid runoff + excess-soil-ice solid runoff) is a sink in the column water
    # balance below (BalanceCheckMod.F90:650), so it must be computed before the
    # balance check — Fortran likewise calls lnd2atm before BalanceCheck.
    #
    # clm_instInit! allocates the outputs; tolerate hand-built instances (and match
    # the working backend/precision) by allocating on first use.
    if length(l2a.qflx_ice_runoff_col) != nc
        l2a.qflx_ice_runoff_col   = _zlike(nc)
        l2a.qflx_liq_from_ice_col = _zlike(nc)
    end
    if length(l2a.eflx_sh_ice_to_liq_col) != nc
        l2a.eflx_sh_ice_to_liq_col = _zlike(nc)
    end
    handle_ice_runoff!(l2a, wfb, col, lun, bc_col)
    # Melted ice runoff (only when melt_non_icesheet_ice_runoff; off by default)
    # leaves via qflx_qrgwl instead, which the balance check already accounts for.
    add_liq_from_ice_to_runoff!(wfb, l2a, col, bc_col)

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
    # Column-level end water balance: the port never set endwb_col (it stayed NaN,
    # so the errh2o_col check silently passed on NaN). Mirror begin (non-lake +
    # lake stores) and add p2c canopy water, consistent with begwb above.
    end_water_column_balance!(inst.water, ls, col, filt.nolakec, filt.lakec, bc_col)
    add_canopy_water_to_storage!(inst.water.waterbalancebulk_inst.endwb_col,
                                 wsb.ws.liqcan_patch, wsb.ws.snocan_patch,
                                 filt.nolakec, col, pch)

    # Gridcell water balance: mirror the column fix so errh2o_grc is meaningful.
    # (1) add canopy to endwb_grc with the END-of-step canopy state (begwb_grc got
    # its canopy at the top of the step, with the BEGIN-of-step state — adding both
    # here from the same liqcan/snocan cancelled the canopy storage change out of
    # the balance); (2) aggregate column output fluxes to gridcell (the call below
    # previously passed _zlike → errh2o_grc compared ΔS to precip alone). Both are
    # required: once the column fluxes are de-NaN'd, an unfixed grc check spams.
    add_canopy_water_to_grc_storage!(inst.water.waterbalancebulk_inst.endwb_grc,
                                     wsb.ws.liqcan_patch, wsb.ws.snocan_patch,
                                     filt.nolakec, col, lun, pch, bc_col, bc_grc)

    # Total water storage (tws_grc) for the CH4 TWS_inversion finundation method.
    # lnd2atmMod.F90:470-483 sets tws_grc = c2g(endwb_col) + river-storage depth;
    # with no active river model the river term is 0, so tws is the gridcell-total
    # end-of-step column water storage. Written here (end of step) so the NEXT
    # step's ch4! reads it — matching CTSM, where lnd2atm runs after the step and
    # ch4 consumes the persisted value. Only consulted by the inversion methods;
    # the default h2osfc finundation path never reads it.
    let wdbulk = inst.water.waterdiagnosticbulk_inst
        if !isempty(wdbulk.tws_grc)
            copyto!(wdbulk.tws_grc, inst.water.waterbalancebulk_inst.endwb_grc)
        end
    end
    _g_evap_tot     = _zlike(length(grc.lat)); _g_surf      = _zlike(length(grc.lat))
    _g_qrgwl        = _zlike(length(grc.lat)); _g_drain     = _zlike(length(grc.lat))
    _g_drain_perch  = _zlike(length(grc.lat)); _g_sfc_irrig = _zlike(length(grc.lat))
    _g_ice_runoff   = _zlike(length(grc.lat))
    # 'urbanf' c2l scaling, NOT unity — Fortran BalanceCheckMod aggregates every
    # water-balance column term with c2l_scale_type='urbanf' (walls x 3*canyon_hwr,
    # roads x 3), because those fluxes are per-m2 of WALL / ROAD area, not of
    # GROUND area. With unity the urban gridcell could not close even when every
    # column closed to 1e-13. Identical to c2g_unity! with no urban landunit.
    c2g_urbanf!(_g_evap_tot,    wfb.wf.qflx_evap_tot_col,      col, lun, bc_col, bc_grc)
    c2g_urbanf!(_g_surf,        wfb.wf.qflx_surf_col,          col, lun, bc_col, bc_grc)
    c2g_urbanf!(_g_qrgwl,       wfb.wf.qflx_qrgwl_col,         col, lun, bc_col, bc_grc)
    c2g_urbanf!(_g_drain,       wfb.wf.qflx_drain_col,         col, lun, bc_col, bc_grc)
    c2g_urbanf!(_g_drain_perch, wfb.wf.qflx_drain_perched_col, col, lun, bc_col, bc_grc)
    c2g_urbanf!(_g_sfc_irrig,   wfb.wf.qflx_sfc_irrig_col,     col, lun, bc_col, bc_grc)
    # qflx_rofice_grc = c2g(qflx_ice_runoff_col) (lnd2atmMod.F90:459). The dynbal
    # correction Fortran subtracts there is zero in the port (no dynamic landunits
    # coupling), so the plain aggregate is the whole term.
    c2g_urbanf!(_g_ice_runoff,  l2a.qflx_ice_runoff_col,      col, lun, bc_col, bc_grc)

    # BalanceCheck — WIRED
    # DAnstep = steps since run start / restart / last DA state jump (Fortran
    # get_nstep_since_startup_or_lastDA_restart_or_pause = get_nstep() - DA_nstep).
    # This used to be a hardcoded 0, which made every hard-error branch in
    # balance_check! unreachable (they all require DAnstep > skip_steps), so the
    # top-level water/energy conservation check could only ever @warn. Four real
    # balance bugs hid behind that. Pass the real counter.
    balance_check!(inst.balcheck, wfb, wsb,
                   inst.water.waterbalancebulk_inst, wdb, ef, sa, cs, alb,
                   col, lun, pch, grc,
                   filt.allc, bc_col, bc_patch, bc_grc,
                   nstep, get_nstep_since_startup_or_last_da(inst.balcheck, nstep), dtime;
                   forc_rain_col=forc_rain_col,
                   forc_snow_col=forc_snow_col,
                   forc_rain_grc=length(a2l.forc_rain_not_downscaled_grc) > 0 ? a2l.forc_rain_not_downscaled_grc : _zlike(length(grc.lat)),
                   forc_snow_grc=length(a2l.forc_snow_not_downscaled_grc) > 0 ? a2l.forc_snow_not_downscaled_grc : _zlike(length(grc.lat)),
                   forc_solad_col=a2l.forc_solad_downscaled_col,
                   forc_solai_grc=a2l.forc_solai_grc,
                   forc_lwrad_col=a2l.forc_lwrad_downscaled_col,
                   forc_flood_grc=_zlike(length(grc.lat)),
                   qflx_ice_runoff_col=l2a.qflx_ice_runoff_col,
                   qflx_evap_tot_grc=_g_evap_tot,
                   qflx_surf_grc=_g_surf,
                   qflx_qrgwl_grc=_g_qrgwl,
                   qflx_drain_grc=_g_drain,
                   qflx_drain_perched_grc=_g_drain_perch,
                   qflx_ice_runoff_grc=_g_ice_runoff,
                   qflx_sfc_irrig_grc=_g_sfc_irrig,
                   qflx_streamflow_grc=_zlike(length(grc.lat)),
                   # FATES trips the (now live) check on its own real coupling gaps
                   # — fsa/fsr unset on FATES patches, ~142 W/m2 solar residual. The
                   # errors are still computed and warned; only the abort is held
                   # back. See the TODO at the top of balance_check!.
                   use_fates=config.use_fates)

    # ========================================================================
    # Diagnostics
    # ========================================================================
    write_diagnostic(nstep)

    # ========================================================================
    # Update accumulators
    # ========================================================================
    if nstep > 0
        # Fortran ORDER (clm_driver.F90): atm2lnd → temperature → canopystate →
        # water → energyflux → bgc_vegetation(dgvs) → crop → fates. Preserved below.
        #
        # forc_rh_grc is DERIVED here (this port has no coupler import of RH); it is
        # the source term for the RH30 / RH24 accumulators.
        atm2lnd_update_rh!(a2l, bc_grc)

        # AnnET is registered under use_fun in Fortran (WaterFluxBulkType imports it
        # from CNSharedParamsMod). In this port use_fun lives on the CN driver config.
        _use_fun_acc = config.use_cn && inst.bgc_vegetation.driver_config.use_fun

        # Atm2Lnd accumulator update — WIRED. Also carries the Fortran
        # `wateratm2lndbulk_type` accumulators (PREC10/30/60/365/24, RH30/24),
        # which live on Atm2LndData in this port.
        atm2lnd_update_acc_vars!(a2l, bc_patch, pch.gridcell, pch.column;
            bounds_c = bc_col, nstep = nstep, dtime = Int(dtime))

        # Temperature accumulators — WIRED. The crop growing-degree-day path
        # (GDD0/8/10 runaccum + GDD020/820/1020 20-yr runmeans) activates only
        # when use_crop. Its AccumManager lives on the (non-dual) config so the
        # ForwardDiff dual-copy of `inst` never sees its non-numeric AccumField
        # vectors; the GDD fields are registered lazily on the first crop step
        # (Fortran temperature_type%InitAccBuffer, called once at init).
        #
        # SOIL10 needs `upper_soil_layer` — Fortran's CNSharedParams layer
        # containing 0.12 m (first layer whose LOWER interface reaches that depth).
        # zisoi is 1-based with zisoi[1] = 0 and zisoi[i+1] = lower interface of
        # layer i; empty before the vertical grid is built.
        _usl = let zi = zisoi[]
            if isempty(zi)
                0
            else
                l = varpar.nlevgrnd
                for i in 1:varpar.nlevgrnd
                    if (i + 1) <= length(zi) && 0.12 <= zi[i + 1]
                        l = i
                        break
                    end
                end
                l
            end
        end
        if config.use_crop
            crop = inst.crop
            if getfield(config, :accum_mgr) === nothing
                mgr = AccumManager()
                temperature_init_acc_buffer!(mgr, length(bc_patch);
                    active = collect(Bool, @view pch.active[bc_patch]),
                    use_crop = true, step_size = Int(dtime))
                config.accum_mgr = mgr
            end
            temperature_update_acc_vars!(temp, bc_col, bc_patch, lun, pch;
                nstep=nstep, dtime=Int(dtime),
                end_cd=is_end_curr_day, upper_soil_layer=_usl,
                mgr=getfield(config, :accum_mgr), use_crop=true,
                month=mon, day=day, jday=jday,
                secs=secs, is_end_curr_year=is_end_curr_year,
                latdeg=grc.latdeg, gridcell=pch.gridcell,
                active=collect(Bool, pch.active), itype=pch.itype,
                npcropmin=config.npcropmin,
                gdd20_season_start=crop.gdd20_season_start_patch,
                gdd20_season_end=crop.gdd20_season_end_patch)
        else
            temperature_update_acc_vars!(temp, bc_col, bc_patch, lun, pch;
                nstep=nstep, dtime=Int(dtime), secs=secs,
                end_cd=is_end_curr_day, upper_soil_layer=_usl)
        end

        # Canopy state accumulators — WIRED (FSUN24 / FSUN240 / LAI240; the
        # windows are derived from dtime, not hardcoded step counts).
        canopystate_update_acc_vars!(cs, bc_patch; nstep=nstep, dtime=Int(dtime))

        # Water accumulators — WIRED (AnnET [use_fun], SNOW_5D).
        water_update_acc_vars!(inst.water, bc_col;
            nstep=nstep, dtime=Int(dtime), use_fun=_use_fun_acc)

        # Energy flux accumulators — WIRED. BTRANMN = daily min of the HOURLY-AVERAGED
        # btran (BTRANAV), so nstep is needed to detect hour boundaries.
        energyflux_update_acc_vars!(ef, bc_patch;
            end_cd=is_end_curr_day, nstep=nstep, dtime=Int(dtime))

        # LUNA photosynthetic-N acclimation — WIRED. Accumulate 24h climate each
        # step; at end-of-day compute the 10-day means, recompute vcmx25_z/jmx25_z,
        # and clear the 24h buffers. LUNA self-skips while the day accumulators are
        # SPVAL (Fortran first-day bootstrap), so a restart-injected run where
        # vcmx25_z is supplied stays unaffected until the accumulators fill. The
        # forc_pco2/po2_240 + pbot240 10-day means are maintained above by
        # atm2lnd_update_acc_vars!; the LUNA PAR arrays are allocated when use_luna
        # (clm_instInit! → solarabs_init!). o3coefjmax=1: ozone–LUNA coupling off.
        if config.use_luna
            acc24_climate_luna!(cs, ps, alb, sa, temp, pch, filt.exposedvegp, bc_patch, dtime)
            if is_end_curr_day
                # Gather the (few) inputs to host so this control-flow scalar loop is
                # device-safe; no-op on the host path (Array-of-Array is identity).
                _lpg  = pch.gridcell isa Array ? pch.gridcell : Array(pch.gridcell)
                _ldl  = grc.dayl isa Array ? grc.dayl : Array{Float64}(grc.dayl)
                _lmdl = grc.max_dayl isa Array ? grc.max_dayl : Array{Float64}(grc.max_dayl)
                _luna_dayl = zeros(np)
                @inbounds for p in 1:np
                    g = _lpg[p]
                    _luna_dayl[p] = smooth_clamp((_ldl[g] * _ldl[g]) /
                                                 (_lmdl[g] * _lmdl[g]), 0.01, 1.0)
                end
                acc240_climate_luna!(temp, ps, alb, sa, wdb, fv, pch, filt.exposedvegp, bc_patch,
                    a2l.forc_po2_240_patch, a2l.forc_pco2_240_patch,
                    fv.rb1_patch, wdb.rh_af_patch, dtime)
                update_photosynthesis_capacity!(ps, temp, cs, alb, sa, wdb, fv, pch, grc,
                    filt.exposedvegp, bc_patch, _luna_dayl,
                    a2l.forc_pbot240_downscaled_patch, a2l.forc_pco2_240_patch, a2l.forc_po2_240_patch,
                    pftcon.c3psn, pftcon.slatop, pftcon.leafcn, pftcon.rhol, pftcon.taul,
                    fill(1.0, np), luna_params_inst, dtime, NLEVCAN)
                clear24_climate_luna!(sa, ps, temp, pch, filt.exposedvegp, bc_patch)
            end
        end

        # Placeholder: bgc_vegetation_inst%UpdateAccVars!(bounds_proc, ...) [BGC accum]

        # CNDV: accumulate AGDD / AGDDTW each time step (reset annually on Jan 1).
        # Mirrors CNVegetationFacade::UpdateAccVars → dgvs_inst%UpdateAccVars; see
        # CNDVType.F90:UpdateAccVars. Only runs when use_cndv and the dgvs state is
        # sized → default path untouched.
        if config.use_cndv && !isempty(inst.dgvs.agdd_patch)
            cndv_update_acc_vars!(inst.dgvs, inst.dgv_ecophyscon,
                temp.t_a10_patch, temp.t_ref2m_patch, bc_patch, dtime;
                month=mon, day=day, secs=secs)
        end

        # Crop accumulators — WIRED. HUI / GDDACCUM / GDDTSOI are `runaccum` fields
        # reset per-patch when the crop is not live. They were allocated to SPVAL and
        # NEVER written (crop_update_acc_vars! was a no-op stub with no call site), so
        # crop phenology read HUI = 1e36 and treated every crop as instantly mature.
        # Fortran calls this every step (clm_driver.F90, after bgc_vegetation).
        if config.use_crop
            crop_update_acc_vars!(inst.crop, bc_patch, temp.t_ref2m_patch,
                temp.t_soisno_col;
                patch = pch, col = col, pftcon = pftcon, dtime = dtime)
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
            is_end_curr_year=is_end_curr_year,
            is_first_step=is_first_step)

        # CNDV: run the annual dynamic-vegetation driver at the last time step of
        # the year. Mirrors CNVegetationFacade::EndOfTimeStepVegDynamics →
        # CNDVDriver (gated on use_cndv .and. is_end_curr_year .and. .not.
        # is_first_step); see CNDVDriverMod.F90:CNDVDriver. Climate20 + light
        # competition + establishment update fpcgrid_patch, which the per-step
        # dyn_cndv_interp! then maps onto patch.wtcol. Gated off → byte-identical.
        if config.use_cndv && is_end_curr_year && !is_first_step &&
           !isempty(inst.dgvs.fpcgrid_patch)
            _cveg_cf = inst.bgc_vegetation.cnveg_carbonflux_inst
            _cveg_cs = inst.bgc_vegetation.cnveg_carbonstate_inst
            # prec365_col — the 365-day running mean of precipitation, now a real
            # accumulator (atm2lnd_update_acc_vars!, gated on use_cndv). It used to be
            # a ZERO BUFFER, which made CNDV establishment's `prec365 >= prec_min_estab`
            # test false forever: no PFT could ever ESTABLISH. Fall back to zeros only
            # if the accumulator was never allocated.
            _cndv_prec365 = length(a2l.prec365_col) >= nc ? a2l.prec365_col : _zlike(nc)
            # cndv_driver! rebuilds this natveg patch mask in-place from
            # present_patch, so any same-length Bool buffer works as scratch.
            _cndv_natvegp = collect(Bool, filt.natvegp)
            cndv_driver!(inst.dgvs, inst.dgv_ecophyscon,
                a2l.t_mo_min_patch, _cndv_prec365,
                _cveg_cf.annsum_npp_patch, _cveg_cf.annsum_litfall_patch,
                _cveg_cs.deadstemc_patch, _cveg_cs.leafcmax_patch,
                pftcon, pch, lun, _cndv_natvegp, bc_patch, year;
                bounds_gridcell = bc_grc)
        end
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
