# ==========================================================================
# Ported from: src/dyn_subgrid/dynSubgridControlMod.F90 (481 lines)
# Dynamic subgrid control — stores & queries control flags for transient
# (dynamic) land use operation.
#
# Skipped: namelist file I/O (read_namelist open/read/close, getavu/relavu,
#           find_nlgroup_name) and MPI broadcast (shr_mpi_bcast). The namelist
#           values are instead passed as keyword arguments to
#           dyn_subgrid_control_init.
# Ported:   the dyn_subgrid_control_type struct, dynSubgridControl_init
#           (default values + check_namelist_consistency), and all query
#           getters. The Fortran module-private singleton becomes an explicit
#           struct passed as the first argument to each getter.
#
# Note: this is static configuration, kept standalone — it is NOT added to any
#       CLMInstances / dual-copied struct.
# ==========================================================================

"""
    DynSubgridControl

Control flags related to dynamic subgrid (transient land use) operation.

Mirrors Fortran `dyn_subgrid_control_type`. Field names preserved from the
Fortran source for traceability.
"""
Base.@kwdef mutable struct DynSubgridControl
    flanduse_timeseries::String = ""        # transient landuse dataset ("" = none)
    do_transient_pfts::Bool   = false       # apply transient natural PFTs from dataset
    do_transient_crops::Bool  = false       # apply transient crops from dataset
    do_transient_lakes::Bool  = false       # apply transient lakes from dataset
    do_transient_urban::Bool  = false       # apply transient urban from dataset
    do_harvest::Bool          = false       # apply harvest from dataset
    do_grossunrep::Bool       = false       # apply gross unrepresented landcover change

    # whether to reset baseline values of total column water and energy in the
    # first step of the run
    reset_dynbal_baselines::Bool = false

    # Testing-only: whether area changes are allowed at times other than the
    # year boundary. Only used for error-checking, not for controlling model
    # behavior.
    for_testing_allow_non_annual_changes::Bool = false

    # Testing-only: if true, set the dynbal water and energy fluxes to zero.
    # Needed in some tests with daily (rather than annual) glacier dynamics.
    # NOTE: setting this true breaks water and energy conservation.
    for_testing_zero_dynbal_fluxes::Bool = false

    initialized::Bool = false               # whether this object has been initialized
end

"""
    dyn_subgrid_control_init(; kwargs...) -> DynSubgridControl

Initialize the dynamic-subgrid control settings. Mirrors Fortran
`dynSubgridControl_init` (plus the non-I/O parts of `read_namelist`): builds the
struct from the namelist-derived values (defaulted to the Fortran defaults) and
runs `check_namelist_consistency`.

The namelist flags are taken as keyword arguments. The consistency check also
depends on a set of other run-configuration flags (from `clm_varctl`), which are
likewise taken as keyword arguments defaulted to their Fortran defaults:
`use_cndv`, `use_fates`, `use_cn`, `use_crop`, `collapse_urban`, `n_dom_pfts`,
`n_dom_landunits`, and the `toosmall_*` thresholds.

`masterproc` controls whether the consistency check runs (matching the Fortran,
where `check_namelist_consistency` is only called on the master process).
"""
function dyn_subgrid_control_init(;
        flanduse_timeseries::String = "",
        do_transient_pfts::Bool   = false,
        do_transient_crops::Bool  = false,
        do_transient_lakes::Bool  = false,
        do_transient_urban::Bool  = false,
        do_harvest::Bool          = false,
        do_grossunrep::Bool       = false,
        reset_dynbal_baselines::Bool = false,
        for_testing_allow_non_annual_changes::Bool = false,
        for_testing_zero_dynbal_fluxes::Bool = false,
        # --- run-configuration flags consulted by check_namelist_consistency ---
        masterproc::Bool = true,
        use_cndv::Bool   = false,
        use_fates::Bool  = false,
        use_cn::Bool     = false,
        use_crop::Bool   = false,
        collapse_urban::Bool = false,
        n_dom_pfts::Int      = 0,
        n_dom_landunits::Int = 0,
        toosmall_soil::Float64    = 0.0,
        toosmall_crop::Float64    = 0.0,
        toosmall_glacier::Float64 = 0.0,
        toosmall_lake::Float64    = 0.0,
        toosmall_wetland::Float64 = 0.0,
        toosmall_urban::Float64   = 0.0,
    )

    ctl = DynSubgridControl(
        flanduse_timeseries = flanduse_timeseries,
        do_transient_pfts   = do_transient_pfts,
        do_transient_crops  = do_transient_crops,
        do_transient_lakes  = do_transient_lakes,
        do_transient_urban  = do_transient_urban,
        do_harvest          = do_harvest,
        do_grossunrep       = do_grossunrep,
        reset_dynbal_baselines = reset_dynbal_baselines,
        for_testing_allow_non_annual_changes = for_testing_allow_non_annual_changes,
        for_testing_zero_dynbal_fluxes       = for_testing_zero_dynbal_fluxes,
    )

    if masterproc
        check_namelist_consistency(ctl;
            use_cndv = use_cndv,
            use_fates = use_fates,
            use_cn = use_cn,
            use_crop = use_crop,
            collapse_urban = collapse_urban,
            n_dom_pfts = n_dom_pfts,
            n_dom_landunits = n_dom_landunits,
            toosmall_soil = toosmall_soil,
            toosmall_crop = toosmall_crop,
            toosmall_glacier = toosmall_glacier,
            toosmall_lake = toosmall_lake,
            toosmall_wetland = toosmall_wetland,
            toosmall_urban = toosmall_urban,
        )
    end

    ctl.initialized = true

    return ctl
end

"""
    check_namelist_consistency(ctl::DynSubgridControl; kwargs...)

Check consistency of the dynamic-subgrid namelist settings. Mirrors Fortran
`check_namelist_consistency`: throws an `error` (Fortran `endrun`) on each
incompatible flag combination, with messages matching the Fortran source.
"""
function check_namelist_consistency(ctl::DynSubgridControl;
        use_cndv::Bool   = false,
        use_fates::Bool  = false,
        use_cn::Bool     = false,
        use_crop::Bool   = false,
        collapse_urban::Bool = false,
        n_dom_pfts::Int      = 0,
        n_dom_landunits::Int = 0,
        toosmall_soil::Float64    = 0.0,
        toosmall_crop::Float64    = 0.0,
        toosmall_glacier::Float64 = 0.0,
        toosmall_lake::Float64    = 0.0,
        toosmall_wetland::Float64 = 0.0,
        toosmall_urban::Float64   = 0.0,
    )

    # All transient/harvest/gross-unrep options require a flanduse_timeseries file.
    if ctl.flanduse_timeseries == ""
        if ctl.do_transient_pfts
            error("ERROR: do_transient_pfts can only be true if you are running with " *
                  "a flanduse_timeseries file (currently flanduse_timeseries is blank)")
        end
        if ctl.do_transient_crops
            error("ERROR: do_transient_crops can only be true if you are running with " *
                  "a flanduse_timeseries file (currently flanduse_timeseries is blank)")
        end
        if ctl.do_transient_lakes
            error("ERROR: do_transient_lakes can only be true if you are running with " *
                  "a flanduse_timeseries file (currently flanduse_timeseries is blank)")
        end
        if ctl.do_transient_urban
            error("ERROR: do_transient_urban can only be true if you are running with " *
                  "a flanduse_timeseries file (currently flanduse_timeseries is blank)")
        end
        if ctl.do_harvest
            error("ERROR: do_harvest can only be true if you are running with " *
                  "a flanduse_timeseries file (currently flanduse_timeseries is blank)")
        end
        if ctl.do_grossunrep
            error("ERROR: do_grossunrep can only be true if you are running with " *
                  "a flanduse_timeseries file (currently flanduse_timeseries is blank)")
        end
    end

    if ctl.do_transient_pfts
        if use_cndv
            error("ERROR: do_transient_pfts is incompatible with use_cndv")
        end
        if use_fates
            error("ERROR: do_transient_pfts is incompatible with use_fates")
        end
    end

    # NOTE(wjs, 2020-08-23) do_transient_lakes is treated like do_transient_pfts
    # and do_transient_crops here, to keep transient lakes consistent with the
    # other transient areas.
    if (ctl.do_transient_pfts || ctl.do_transient_crops ||
        ctl.do_transient_lakes || ctl.do_transient_urban)
        if collapse_urban
            error("ERROR: do_transient_pfts, do_transient_crops, do_transient_lakes and " *
                  "do_transient_urban are incompatible with collapse_urban = .true.")
        end
        if (n_dom_pfts > 0 || n_dom_landunits > 0 ||
            toosmall_soil > 0.0 || toosmall_crop > 0.0 ||
            toosmall_glacier > 0.0 || toosmall_lake > 0.0 ||
            toosmall_wetland > 0.0 || toosmall_urban > 0.0)
            error("ERROR: do_transient_pfts, do_transient_crops and do_transient_lakes and " *
                  "do_transient_urban are incompatible with any of the following set to > 0: " *
                  "n_dom_pfts > 0, n_dom_landunits > 0, " *
                  "toosmall_soil > 0., toosmall_crop > 0., " *
                  "toosmall_glacier > 0., toosmall_lake > 0., " *
                  "toosmall_wetland > 0., toosmall_urban > 0.")
        end
    end

    if ctl.do_transient_crops
        if use_fates
            # ED / FATES does not currently have a mechanism for changing its
            # column areas. See https://github.com/NGEET/ed-clm/issues/173
            error("ERROR: do_transient_crops does not currently work with use_fates")
        end
    end

    if ctl.do_harvest
        if !(use_cn || use_fates)
            error("ERROR: do_harvest can only be true if either use_cn or use_fates are true")
        end
    end

    if ctl.do_grossunrep
        # Check use_fates first; the .not. use_cn error appears only if
        # .not. use_fates and .not. use_cn.
        if use_fates
            error("ERROR: do_grossunrep currently does not work with use_fates")
        end
        if !use_cn
            error("ERROR: do_grossunrep can only be true if use_cn is true")
        end
    end

    return nothing
end

# --------------------------------------------------------------------------
# Query functions (Fortran getters). In the Fortran these read a module-private
# singleton; here they take the struct as the first argument.
#
# The Fortran getters SHR_ASSERT that the object is initialized; mirror that.
# --------------------------------------------------------------------------

function get_flanduse_timeseries(ctl::DynSubgridControl)
    @assert ctl.initialized "dyn_subgrid_control is not initialized"
    return ctl.flanduse_timeseries
end

function get_do_transient_pfts(ctl::DynSubgridControl)
    @assert ctl.initialized "dyn_subgrid_control is not initialized"
    return ctl.do_transient_pfts
end

function get_do_transient_crops(ctl::DynSubgridControl)
    @assert ctl.initialized "dyn_subgrid_control is not initialized"
    return ctl.do_transient_crops
end

function get_do_transient_lakes(ctl::DynSubgridControl)
    @assert ctl.initialized "dyn_subgrid_control is not initialized"
    return ctl.do_transient_lakes
end

function get_do_transient_urban(ctl::DynSubgridControl)
    @assert ctl.initialized "dyn_subgrid_control is not initialized"
    return ctl.do_transient_urban
end

function get_do_harvest(ctl::DynSubgridControl)
    @assert ctl.initialized "dyn_subgrid_control is not initialized"
    return ctl.do_harvest
end

function get_do_grossunrep(ctl::DynSubgridControl)
    @assert ctl.initialized "dyn_subgrid_control is not initialized"
    return ctl.do_grossunrep
end

function get_reset_dynbal_baselines(ctl::DynSubgridControl)
    @assert ctl.initialized "dyn_subgrid_control is not initialized"
    return ctl.reset_dynbal_baselines
end

function get_for_testing_allow_non_annual_changes(ctl::DynSubgridControl)
    @assert ctl.initialized "dyn_subgrid_control is not initialized"
    return ctl.for_testing_allow_non_annual_changes
end

function get_for_testing_zero_dynbal_fluxes(ctl::DynSubgridControl)
    @assert ctl.initialized "dyn_subgrid_control is not initialized"
    return ctl.for_testing_zero_dynbal_fluxes
end

"""
    run_has_transient_landcover(ctl::DynSubgridControl) -> Bool

Returns true if any aspects of prescribed transient landcover are enabled.
Mirrors Fortran `run_has_transient_landcover` (note: this checks pfts, crops,
and urban — not lakes — matching the Fortran).
"""
function run_has_transient_landcover(ctl::DynSubgridControl)
    return (get_do_transient_pfts(ctl) ||
            get_do_transient_crops(ctl) ||
            get_do_transient_urban(ctl))
end
