# ==========================================================================
# Ported from: src/biogeophys/WaterType.F90
# Master water container for bulk water and water tracers.
#
# Variables pertaining to bulk water can be accessed directly:
#     water.waterfluxbulk_inst
#     water.waterstatebulk_inst
#     water.waterdiagnosticbulk_inst
#     water.waterbalancebulk_inst
#
# Or as one of the indices in water.bulk_and_tracers:
#     water.bulk_and_tracers[water.i_bulk].waterflux
#
# To loop through bulk and all tracers:
#     for i in water.bulk_and_tracers_beg:water.bulk_and_tracers_end
#         bt = water.bulk_and_tracers[i]
#         # use bt.waterflux, bt.waterstate, etc.
#     end
#
# To loop through just tracers (not bulk):
#     for i in water.tracers_beg:water.tracers_end
#         bt = water.bulk_and_tracers[i]
#         # use bt.waterflux, bt.waterstate, etc.
#     end
#
# Note: Julia uses 1-based indexing. Bulk is at index i_bulk (=1).
# Tracers start at tracers_beg (=2). Fortran uses 0-based with bulk at 0.
# ==========================================================================

# ---------------------------------------------------------------------------
# Concrete info types (minimal implementations of WaterInfoBulkType.F90 and
# WaterInfoIsotopeType.F90, needed by the water container)
# ---------------------------------------------------------------------------

"""
    WaterInfoBulkData <: WaterInfoBaseType

Info type for bulk water. Implements the WaterInfoBaseType interface.

Ported from `water_info_bulk_type` in `WaterInfoBulkType.F90`.
"""
struct WaterInfoBulkData <: WaterInfoBaseType
    ratio::Float64
end

WaterInfoBulkData() = WaterInfoBulkData(1.0)

get_name(::WaterInfoBulkData) = "bulk"
fname(::WaterInfoBulkData, basename::String) = basename
lname(::WaterInfoBulkData, basename::String) = basename
is_communicated_with_coupler(::WaterInfoBulkData) = true
is_included_in_consistency_check(::WaterInfoBulkData) = false

"""
    WaterInfoIsotopeData <: WaterInfoBaseType

Info type for isotope/tracer water. Implements the WaterInfoBaseType interface.

Ported from `water_info_isotope_type` in `WaterInfoIsotopeType.F90`.
"""
struct WaterInfoIsotopeData <: WaterInfoBaseType
    tracer_name::String
    ratio::Float64
    included_in_consistency_check::Bool
    communicated_with_coupler::Bool
end

get_name(info::WaterInfoIsotopeData) = info.tracer_name
fname(info::WaterInfoIsotopeData, basename::String) = basename * "_" * info.tracer_name
lname(info::WaterInfoIsotopeData, basename::String) = basename * " " * info.tracer_name
is_communicated_with_coupler(info::WaterInfoIsotopeData) = info.communicated_with_coupler
is_included_in_consistency_check(info::WaterInfoIsotopeData) = info.included_in_consistency_check

# ---------------------------------------------------------------------------
# WaterParams — parameters controlling tracer configuration
# ---------------------------------------------------------------------------

"""
    WaterParams

Parameters controlling water tracer configuration.

Ported from `water_params_type` in `WaterType.F90`.
"""
Base.@kwdef mutable struct WaterParams
    enable_consistency_checks::Bool = false
    enable_isotopes::Bool = false
end

# ---------------------------------------------------------------------------
# BulkOrTracerData — container for one bulk or tracer instance
# ---------------------------------------------------------------------------

"""
    BulkOrTracerData

Container holding water sub-instances for either bulk water or a single tracer.

For bulk water, the fields reference the Bulk-specific data types
(WaterFluxBulkData, WaterStateBulkData, etc.).
For tracers, the fields hold base data types (WaterFluxData, WaterStateData, etc.).

Ported from `bulk_or_tracer_type` in `WaterType.F90`.
"""
Base.@kwdef mutable struct BulkOrTracerData
    waterflux::Union{WaterFluxData, WaterFluxBulkData, Nothing} = nothing
    waterstate::Union{WaterStateData, WaterStateBulkData, Nothing} = nothing
    waterdiagnostic::Union{WaterDiagnosticBulkData, Nothing} = nothing
    waterbalance::Union{WaterBalanceData, Nothing} = nothing
    # waterlnd2atm and wateratm2lnd: will be added when those types are ported

    is_isotope::Bool = false
    info::Union{WaterInfoBulkData, WaterInfoIsotopeData, Nothing} = nothing
    # vars::WaterTracerContainerType — not yet ported
end

# ---------------------------------------------------------------------------
# WaterData — master water container
# ---------------------------------------------------------------------------

"""
    WaterData

Master water container. Holds all water-related data for bulk water and tracers.

Ported from `water_type` in `WaterType.F90`.
"""
# Parametric on working precision FT. FT is a *construction-precision tag*: the keyword
# constructor builds the bulk sub-instances at FT, but the bulk fields stay LOOSE
# (un-pinned UnionAll) types so the AD path can swap in `Dual`-typed sub-instances on a
# Float64 container (and so device adapt can swap array storage). NOT @kwdef: a single
# type parameter can't keep @kwdef's synthesised constructors and also build children at
# FT via `WaterData{FT}()` (same signature), so the keyword constructor is explicit below.
mutable struct WaterData{FT<:Real}
    # --- Public index bounds (1-based; Fortran used 0-based) ---
    bulk_and_tracers_beg::Int   # first index for bulk & tracers
    bulk_and_tracers_end::Int   # last index for bulk & tracers
    tracers_beg::Int            # first index for just tracers
    tracers_end::Int            # last index for just tracers (empty range = no tracers)
    i_bulk::Int                 # index of bulk in bulk_and_tracers

    # --- Direct bulk instance access (loose types: see note above) ---
    waterfluxbulk_inst::WaterFluxBulkData
    waterstatebulk_inst::WaterStateBulkData
    waterdiagnosticbulk_inst::WaterDiagnosticBulkData
    waterbalancebulk_inst::WaterBalanceData
    # waterlnd2atmbulk_inst  — not yet ported
    # wateratm2lndbulk_inst  — not yet ported

    # --- Iteration array for bulk + tracers ---
    bulk_and_tracers::Vector{BulkOrTracerData}

    # --- Private ---
    params::WaterParams
    bulk_tracer_index::Int      # index of tracer replicating bulk (-1 if none)
end

# Keyword constructor (field-name keyword API used by the Adapt rule and callers). Bulk
# sub-instances default to working precision FT.
function WaterData{FT}(;
        bulk_and_tracers_beg::Int = 1,
        bulk_and_tracers_end::Int = 1,
        tracers_beg::Int = 2,
        tracers_end::Int = 1,
        i_bulk::Int = 1,
        waterfluxbulk_inst::WaterFluxBulkData       = WaterFluxBulkData{FT}(),
        waterstatebulk_inst::WaterStateBulkData     = WaterStateBulkData{FT}(),
        waterdiagnosticbulk_inst::WaterDiagnosticBulkData = WaterDiagnosticBulkData{FT}(),
        waterbalancebulk_inst::WaterBalanceData     = WaterBalanceData{FT}(),
        bulk_and_tracers::Vector{BulkOrTracerData}  = BulkOrTracerData[],
        params::WaterParams                         = WaterParams(),
        bulk_tracer_index::Int = -1) where {FT<:Real}
    WaterData{FT}(bulk_and_tracers_beg, bulk_and_tracers_end, tracers_beg, tracers_end,
                  i_bulk, waterfluxbulk_inst, waterstatebulk_inst, waterdiagnosticbulk_inst,
                  waterbalancebulk_inst, bulk_and_tracers, params, bulk_tracer_index)
end

# Default precision is Float64; legacy callers that don't specify FT get a Float64 tree.
WaterData(; kwargs...) = WaterData{Float64}(; kwargs...)

# Custom Adapt rule for WaterData. The default @adapt_structure can't be used:
# `bulk_and_tracers` holds BulkOrTracerData elements that *alias* the direct bulk
# sub-struct fields (waterstatebulk_inst, ...). We move the bulk sub-structs to the
# device, then rebuild the iteration vector so its references point at the moved
# instances — preserving the aliasing that the bulk/tracer loops rely on.
function Adapt.adapt_structure(to, w::WaterData{FT}) where {FT}
    wf = Adapt.adapt(to, w.waterfluxbulk_inst)
    ws = Adapt.adapt(to, w.waterstatebulk_inst)
    wd = Adapt.adapt(to, w.waterdiagnosticbulk_inst)
    wb = Adapt.adapt(to, w.waterbalancebulk_inst)
    # substitute moved bulk insts where an element aliases them; adapt other
    # (tracer) references directly; leave `nothing` as-is.
    sub(x, orig, moved) = x === orig ? moved : (x === nothing ? nothing : Adapt.adapt(to, x))
    bt = [BulkOrTracerData(
              waterflux       = sub(e.waterflux, w.waterfluxbulk_inst, wf),
              waterstate      = sub(e.waterstate, w.waterstatebulk_inst, ws),
              waterdiagnostic = sub(e.waterdiagnostic, w.waterdiagnosticbulk_inst, wd),
              waterbalance    = sub(e.waterbalance, w.waterbalancebulk_inst, wb),
              is_isotope      = e.is_isotope,
              info            = e.info)
          for e in w.bulk_and_tracers]
    return WaterData{FT}(
        bulk_and_tracers_beg = w.bulk_and_tracers_beg,
        bulk_and_tracers_end = w.bulk_and_tracers_end,
        tracers_beg = w.tracers_beg,
        tracers_end = w.tracers_end,
        i_bulk = w.i_bulk,
        waterfluxbulk_inst = wf,
        waterstatebulk_inst = ws,
        waterdiagnosticbulk_inst = wd,
        waterbalancebulk_inst = wb,
        bulk_and_tracers = bt,
        params = w.params,
        bulk_tracer_index = w.bulk_tracer_index)
end

# -----------------------------------------------------------------------
# water_init! — Main initialization
# Ported from Init in WaterType.F90
# Note: ReadNamelist is omitted (I/O infrastructure not ported);
# tracer configuration is passed directly via keyword arguments.
# -----------------------------------------------------------------------

"""
    water_init!(water, nc, np, nl, ng; enable_consistency_checks, enable_isotopes)

Initialize all water variables. Allocates all sub-instances for `nc` columns,
`np` patches, `nl` landunits, and `ng` gridcells.

Ported from `Init` in `WaterType.F90`. ReadNamelist is replaced by keyword args.
"""
function water_init!(water::WaterData, nc::Int, np::Int, nl::Int, ng::Int;
                     enable_consistency_checks::Bool = false,
                     enable_isotopes::Bool = false)
    water.params = WaterParams(
        enable_consistency_checks = enable_consistency_checks,
        enable_isotopes = enable_isotopes)

    water_do_init!(water, nc, np, nl, ng)
    return nothing
end

# -----------------------------------------------------------------------
# water_init_for_testing! — Init for unit tests
# Ported from InitForTesting in WaterType.F90
# -----------------------------------------------------------------------

"""
    water_init_for_testing!(water, nc, np, nl, ng; params)

Version of Init routine just for unit tests.
Params are passed in directly instead of reading from namelist.

Ported from `InitForTesting` in `WaterType.F90`.
"""
function water_init_for_testing!(water::WaterData, nc::Int, np::Int, nl::Int, ng::Int;
                                  params::WaterParams = WaterParams())
    water.params = params
    water_do_init!(water, nc, np, nl, ng)
    return nothing
end

# -----------------------------------------------------------------------
# water_do_init! — Shared initialization (private)
# Ported from DoInit in WaterType.F90
# -----------------------------------------------------------------------

function water_do_init!(water::WaterData, nc::Int, np::Int, nl::Int, ng::Int)
    # Set up tracer info and allocate the bulk_and_tracers array
    water_setup_tracer_info!(water)

    # Allocate bulk instances and link into bulk_and_tracers
    water_allocate_bulk!(water)

    # Initialize (allocate arrays in) bulk sub-instances
    waterfluxbulk_init!(water.waterfluxbulk_inst, nc, np, nl, ng)
    waterstatebulk_init!(water.waterstatebulk_inst, nc, np, nl, ng)
    waterdiagnosticbulk_init!(water.waterdiagnosticbulk_inst, nc, np, nl, ng)
    waterbalance_init!(water.waterbalancebulk_inst, nc, np, ng)

    # Initialize tracer sub-instances
    for i in water.tracers_beg:water.tracers_end
        water_allocate_tracer!(water, i)
        bt = water.bulk_and_tracers[i]
        if bt.waterflux isa WaterFluxData
            waterflux_init!(bt.waterflux, nc, np, nl, ng)
        end
        if bt.waterstate isa WaterStateData
            waterstate_init!(bt.waterstate, nc, np, nl, ng)
        end
        if bt.waterbalance isa WaterBalanceData
            waterbalance_init!(bt.waterbalance, nc, np, ng)
        end
        # WaterDiagnosticData (base) not yet ported — tracer diagnostics unavailable
    end

    return nothing
end

# -----------------------------------------------------------------------
# water_setup_tracer_info! — Setup information on each water tracer
# Ported from SetupTracerInfo in WaterType.F90
# -----------------------------------------------------------------------

function water_setup_tracer_info!(water::WaterData)
    water.bulk_tracer_index = -1

    num_tracers = 0
    enable_bulk_tracer = false

    if water.params.enable_consistency_checks || water.params.enable_isotopes
        enable_bulk_tracer = true
    end

    if enable_bulk_tracer
        num_tracers += 1
    end
    if water.params.enable_isotopes
        num_tracers += 2
    end
    if water.params.enable_consistency_checks
        num_tracers += 3
    end

    # Set index bounds (Julia 1-based; Fortran 0-based)
    water.bulk_and_tracers_beg = 1
    water.tracers_beg = 2
    water.bulk_and_tracers_end = num_tracers + 1
    water.tracers_end = num_tracers + 1
    water.i_bulk = 1

    # Allocate the array (index 1 = bulk, indices 2:end = tracers)
    water.bulk_and_tracers = [BulkOrTracerData() for _ in 1:(num_tracers + 1)]

    # Bulk info at index i_bulk
    water.bulk_and_tracers[water.i_bulk].info = WaterInfoBulkData()

    # Tracers start at index 2 (Fortran: tracer_num starts at 1)
    tracer_idx = 2

    if enable_bulk_tracer
        water.bulk_and_tracers[tracer_idx].info = WaterInfoIsotopeData(
            "H2OTR", 1.0, true, false)
        water.bulk_and_tracers[tracer_idx].is_isotope = true
        water.bulk_tracer_index = tracer_idx
        tracer_idx += 1
    end

    if water.params.enable_isotopes
        water.bulk_and_tracers[tracer_idx].info = WaterInfoIsotopeData(
            "HDO", 0.9, false, false)
        water.bulk_and_tracers[tracer_idx].is_isotope = true
        tracer_idx += 1

        water.bulk_and_tracers[tracer_idx].info = WaterInfoIsotopeData(
            "H218O", 0.5, false, false)
        water.bulk_and_tracers[tracer_idx].is_isotope = true
        tracer_idx += 1
    end

    if water.params.enable_consistency_checks
        water.bulk_and_tracers[tracer_idx].info = WaterInfoIsotopeData(
            "TESTMED", 0.1, true, false)
        water.bulk_and_tracers[tracer_idx].is_isotope = true
        tracer_idx += 1

        water.bulk_and_tracers[tracer_idx].info = WaterInfoIsotopeData(
            "TESTSMALL", 1.0e-10, true, false)
        water.bulk_and_tracers[tracer_idx].is_isotope = true
        tracer_idx += 1

        water.bulk_and_tracers[tracer_idx].info = WaterInfoIsotopeData(
            "TESTBIG", 10.0, true, false)
        water.bulk_and_tracers[tracer_idx].is_isotope = true
        tracer_idx += 1
    end

    # Verify count (tracer_idx started at 2, so tracers added = tracer_idx - 2)
    if tracer_idx - 2 != num_tracers
        error("water_setup_tracer_info!: tracer count discrepancy: " *
              "expected $num_tracers, got $(tracer_idx - 2)")
    end

    return nothing
end

# -----------------------------------------------------------------------
# water_allocate_bulk! — Allocate bulk objects and link into bulk_and_tracers
# Ported from AllocateBulk in WaterType.F90
# -----------------------------------------------------------------------

function water_allocate_bulk!(water::WaterData)
    bt = water.bulk_and_tracers[water.i_bulk]

    # Link the bulk instances into the iteration array
    bt.waterflux = water.waterfluxbulk_inst
    bt.waterstate = water.waterstatebulk_inst
    bt.waterdiagnostic = water.waterdiagnosticbulk_inst
    bt.waterbalance = water.waterbalancebulk_inst
    # bt.waterlnd2atm = water.waterlnd2atmbulk_inst  — not yet ported
    # bt.wateratm2lnd = water.wateratm2lndbulk_inst  — not yet ported

    return nothing
end

# -----------------------------------------------------------------------
# water_allocate_tracer! — Allocate tracer objects for tracer i
# Ported from AllocateTracer in WaterType.F90
# -----------------------------------------------------------------------

function water_allocate_tracer!(water::WaterData{FT}, i::Int) where {FT}
    bt = water.bulk_and_tracers[i]

    bt.waterflux = WaterFluxData{FT}()
    bt.waterstate = WaterStateData{FT}()
    # bt.waterdiagnostic — requires base WaterDiagnosticData (not yet ported)
    bt.waterbalance = WaterBalanceData{FT}()
    # bt.waterlnd2atm — not yet ported
    # bt.wateratm2lnd — not yet ported

    return nothing
end

# -----------------------------------------------------------------------
# water_clean! — Deallocate all sub-instances (Julia-specific)
# -----------------------------------------------------------------------

"""
    water_clean!(water)

Deallocate (reset to empty) all water sub-instances.
"""
function water_clean!(water::WaterData)
    waterfluxbulk_clean!(water.waterfluxbulk_inst)
    waterstatebulk_clean!(water.waterstatebulk_inst)
    waterdiagnosticbulk_clean!(water.waterdiagnosticbulk_inst)
    waterbalance_clean!(water.waterbalancebulk_inst)

    for i in water.tracers_beg:water.tracers_end
        bt = water.bulk_and_tracers[i]
        if bt.waterflux isa WaterFluxData
            waterflux_clean!(bt.waterflux)
        end
        if bt.waterstate isa WaterStateData
            waterstate_clean!(bt.waterstate)
        end
        if bt.waterbalance isa WaterBalanceData
            waterbalance_clean!(bt.waterbalance)
        end
    end

    water.bulk_and_tracers = BulkOrTracerData[]
    water.bulk_tracer_index = -1

    return nothing
end

# -----------------------------------------------------------------------
# Accumulation buffer / variable functions
# Ported from InitAccBuffer, InitAccVars, UpdateAccVars in WaterType.F90
# -----------------------------------------------------------------------

"""
    water_init_acc_buffer!(water, bounds_col)

Initialize accumulation buffer for all water variables.

Ported from `InitAccBuffer` in `WaterType.F90`.
"""
function water_init_acc_buffer!(water::WaterData, bounds_col::UnitRange{Int})
    waterfluxbulk_init_acc_buffer!(water.waterfluxbulk_inst, bounds_col)
    # wateratm2lndbulk_inst%InitAccBuffer — not yet ported
    waterdiagnosticbulk_init_acc_buffer!(water.waterdiagnosticbulk_inst, bounds_col)
    return nothing
end

"""
    water_init_acc_vars!(water, bounds_col)

Initialize variables associated with accumulated fields.

Ported from `InitAccVars` in `WaterType.F90`.
"""
function water_init_acc_vars!(water::WaterData, bounds_col::UnitRange{Int})
    waterfluxbulk_init_acc_vars!(water.waterfluxbulk_inst, bounds_col)
    # wateratm2lndbulk_inst%initAccVars — not yet ported
    waterdiagnosticbulk_init_acc_vars!(water.waterdiagnosticbulk_inst, bounds_col)
    return nothing
end

"""
    water_update_acc_vars!(water, bounds_col; nstep=0, dtime=1800, use_fun=false)

Update accumulated variables. Called every time step, in Fortran's sub-order
(fluxbulk → atm2lndbulk → diagnosticbulk).

BUG HISTORY: this dispatcher was called every step but did LITERALLY NOTHING —
both of its sub-calls were no-op stubs and the third (`wateratm2lndbulk`) was
simply absent. Three accumulator families were silently missing: `AnnET`
(fluxbulk), `SNOW_5D` (diagnosticbulk), and `PREC10/30/60/365/24` + `RH30/24`
(atm2lndbulk). See the individual routines for the downstream consequences.

The `wateratm2lndbulk` accumulators are driven from
[`atm2lnd_update_acc_vars!`](@ref) rather than here: CLM.jl has no separate
`wateratm2lndbulk_type`, and those fields (plus their rain/snow/RH source terms)
live on `Atm2LndData`. The ORDER is preserved — the driver calls
`atm2lnd_update_acc_vars!` before this — and no accumulator reads another's
output within a step, so the split is behaviourally identical.

Ported from `UpdateAccVars` in `WaterType.F90`.
"""
function water_update_acc_vars!(water::WaterData, bounds_col::UnitRange{Int};
                                nstep::Int = 0,
                                dtime::Int = 1800,
                                use_fun::Bool = false)
    waterfluxbulk_update_acc_vars!(water.waterfluxbulk_inst, bounds_col;
        nstep = nstep, dtime = dtime, use_fun = use_fun)
    waterdiagnosticbulk_update_acc_vars!(water.waterdiagnosticbulk_inst, bounds_col;
        nstep = nstep, dtime = dtime)
    return nothing
end

# -----------------------------------------------------------------------
# water_restart! — Read/write restart information
# Ported from Restart in WaterType.F90
# -----------------------------------------------------------------------

"""
    water_restart!(water, bounds_col; flag)

Read/write information to/from restart file for all water variables.

Ported from `Restart` in `WaterType.F90`.
"""
function water_restart!(water::WaterData, bounds_col::UnitRange{Int};
                         flag::String = "read")
    waterfluxbulk_restart!(water.waterfluxbulk_inst, bounds_col; flag=flag)
    waterstatebulk_restart!(water.waterstatebulk_inst, bounds_col; flag=flag)
    waterdiagnosticbulk_restart!(water.waterdiagnosticbulk_inst, bounds_col; flag=flag)

    for i in water.tracers_beg:water.tracers_end
        bt = water.bulk_and_tracers[i]
        if bt.waterflux isa WaterFluxData
            waterflux_restart!(bt.waterflux, bounds_col; flag=flag)
        end
        if bt.waterstate isa WaterStateData
            waterstate_restart!(bt.waterstate, bounds_col; flag=flag)
        end
        # waterdiagnostic for tracers — requires base WaterDiagnosticData (not ported)
    end

    return nothing
end

# -----------------------------------------------------------------------
# Accessor functions
# -----------------------------------------------------------------------

"""
    water_get_bulk_or_tracer_name(water, i) -> String

Get name of the given tracer (or bulk). `i` must be in
`bulk_and_tracers_beg:bulk_and_tracers_end`.

Ported from `GetBulkOrTracerName` in `WaterType.F90`.
"""
function water_get_bulk_or_tracer_name(water::WaterData, i::Int)
    @assert i >= water.bulk_and_tracers_beg "Index $i below bulk_and_tracers_beg"
    @assert i <= water.bulk_and_tracers_end "Index $i above bulk_and_tracers_end"
    return get_name(water.bulk_and_tracers[i].info)
end

"""
    water_is_isotope(water, i) -> Bool

Returns true if tracer `i` is an isotope. `i` must be in
`tracers_beg:tracers_end`.

Ported from `IsIsotope` in `WaterType.F90`.
"""
function water_is_isotope(water::WaterData, i::Int)
    @assert i >= water.tracers_beg "Index $i below tracers_beg"
    @assert i <= water.tracers_end "Index $i above tracers_end"
    return water.bulk_and_tracers[i].is_isotope
end

"""
    water_get_isotope_info(water, i) -> WaterInfoIsotopeData

Get the isotope info object for tracer `i`. `i` must be in
`tracers_beg:tracers_end` and `water_is_isotope(water, i)` must be true.

Ported from `GetIsotopeInfo` in `WaterType.F90`.
"""
function water_get_isotope_info(water::WaterData, i::Int)
    @assert i >= water.tracers_beg "Index $i below tracers_beg"
    @assert i <= water.tracers_end "Index $i above tracers_end"

    info = water.bulk_and_tracers[i].info
    if info isa WaterInfoIsotopeData
        return info
    else
        error("water_get_isotope_info: tracer $i is not an isotope")
    end
end

"""
    water_get_bulk_tracer_index(water) -> Int

Get the index of the tracer that replicates bulk water.
Returns -1 if there is no such tracer in this run.

Ported from `GetBulkTracerIndex` in `WaterType.F90`.
"""
function water_get_bulk_tracer_index(water::WaterData)
    return water.bulk_tracer_index
end

"""
    water_do_consistency_check(water) -> Bool

Returns true if TracerConsistencyCheck should be called in this run.

Ported from `DoConsistencyCheck` in `WaterType.F90`.
"""
function water_do_consistency_check(water::WaterData)
    return water.params.enable_consistency_checks
end

# -----------------------------------------------------------------------
# water_tracer_consistency_check! — Check tracer consistency with bulk
# Ported from TracerConsistencyCheck in WaterType.F90
# Requires WaterTracerContainerType — stub until ported
# -----------------------------------------------------------------------

"""
    water_tracer_consistency_check!(water, bounds_col; caller_location)

Check consistency of water tracers with bulk water.

Should only be called if `water_do_consistency_check(water)` returns true.

Ported from `TracerConsistencyCheck` in `WaterType.F90`.
Stub: requires WaterTracerContainerType for variable iteration.
"""
function water_tracer_consistency_check!(water::WaterData, bounds_col::UnitRange{Int};
                                          caller_location::String = "")
    # Stub: requires WaterTracerContainerType (vars field) for variable iteration
    # and WaterTracerUtils (CompareBulkToTracer)
    return nothing
end

# -----------------------------------------------------------------------
# water_reset_checked_tracers! — Reset checked tracers to bulk * ratio
# Ported from ResetCheckedTracers in WaterType.F90
# Requires WaterTracerContainerType — stub until ported
# -----------------------------------------------------------------------

"""
    water_reset_checked_tracers!(water, bounds_col)

For tracers included in consistency checks, reset all values to bulk * ratio.

Ported from `ResetCheckedTracers` in `WaterType.F90`.
Stub: requires WaterTracerContainerType for variable iteration.
"""
function water_reset_checked_tracers!(water::WaterData, bounds_col::UnitRange{Int})
    # Stub: requires WaterTracerContainerType (vars field) for variable iteration
    # and WaterTracerUtils (SetTracerToBulkTimesRatio)
    return nothing
end

# -----------------------------------------------------------------------
# water_summary! — Compute end-of-timestep summaries
# Ported from Summary in WaterType.F90
# -----------------------------------------------------------------------

"""
    water_summary!(water, bounds_col, bounds_patch; mask_soilp, mask_allc,
                    mask_nolakec, h2osno_total_col, dz_col, zi_col,
                    landunit_col, urbpoi, lun_itype)

Compute end-of-timestep summaries of water diagnostic terms.

Iterates over bulk and all tracers, calling the diagnostic Summary on each.
Waterstate and waterflux fields are extracted from the container's own
sub-instances; external data (column/landunit geometry) must be passed in.

Ported from `Summary` in `WaterType.F90`.
"""
function water_summary!(water::WaterData,
                         bounds_col::UnitRange{Int},
                         bounds_patch::UnitRange{Int};
                         mask_soilp::AbstractVector{Bool},
                         mask_allc::AbstractVector{Bool},
                         mask_nolakec::AbstractVector{Bool},
                         h2osno_total_col::AbstractVector{<:Real},
                         dz_col::AbstractMatrix{<:Real},
                         zi_col::AbstractMatrix{<:Real},
                         landunit_col::AbstractVector{<:Integer},
                         urbpoi::AbstractVector{Bool},
                         lun_itype::AbstractVector{<:Integer})
    # Bulk summary: extract waterstate/waterflux fields from bulk instances
    ws = water.waterstatebulk_inst.ws
    wf = water.waterfluxbulk_inst.wf

    waterdiagnosticbulk_summary!(water.waterdiagnosticbulk_inst,
        bounds_col, bounds_patch;
        mask_soilp = mask_soilp,
        mask_allc = mask_allc,
        mask_nolakec = mask_nolakec,
        h2osoi_ice_col = ws.h2osoi_ice_col,
        h2osoi_liq_col = ws.h2osoi_liq_col,
        excess_ice_col = ws.excess_ice_col,
        qflx_intercepted_liq_patch = wf.qflx_intercepted_liq_patch,
        qflx_intercepted_snow_patch = wf.qflx_intercepted_snow_patch,
        qflx_liq_grnd_col = wf.qflx_liq_grnd_col,
        qflx_snow_grnd_col = wf.qflx_snow_grnd_col,
        h2osno_total_col = h2osno_total_col,
        dz_col = dz_col,
        zi_col = zi_col,
        landunit_col = landunit_col,
        urbpoi = urbpoi,
        lun_itype = lun_itype)

    # Tracer summaries will be added when base WaterDiagnosticData is ported.
    # In Fortran, this iterates over all bulk_and_tracers and calls
    # waterdiagnostic_inst%Summary with the matching waterstate_inst and waterflux_inst.

    return nothing
end
