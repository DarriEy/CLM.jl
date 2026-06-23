# FatesRestartVariableType.jl
# Julia port of FATES src/fates/main/FatesRestartVariableType.F90
#
# The `fates_restart_variable_type` describes ONE restart variable: its
# metadata (name/units/long-name/vtype/flushval) plus a 1D data payload that is
# EITHER a real (r81d) OR an integer (int1d) vector — only one is allocated per
# variable, governed by `vtype`.
#
# Deps: FatesConstantsMod (fates_r8 == Float64), FatesIODimensionsMod
#       (fates_io_dimension_type), FatesIOVariableKindMod
#       (fates_io_variable_kind_type, the site_*/cohort_* kind constants,
#        iotype_index, set_active!), FatesGlobals (fates_log, fates_endrun).
#
# NOTE: Fortran arrays here use explicit lower bounds (lb1:ub1). FATES IO dims
# start their lower bound at 1, so Julia 1-based vectors of length (ub1-lb1+1)
# map directly. We keep lb1/ub1 around for fidelity but allocate length-based.

# ---------------------------------------------------------------------------
# fates_restart_variable_type — one restart variable (metadata + payload).
# ---------------------------------------------------------------------------
Base.@kwdef mutable struct fates_restart_variable_type
    vname::String = ""
    units::String = ""
    long::String = ""
    vtype::String = ""
    flushval::Float64 = 0.0      # DONT THINK THIS IS NEEDED IN RESTARTS
                                 # (RGK-11-2016) — kept for fidelity.
    dim_kinds_index::Int = 0
    # Pointers (only one of these is allocated per variable). Empty == nullified.
    r81d::Vector{Float64} = Float64[]
    int1d::Vector{Int} = Int[]
    # Bookkeeping: whether r81d / int1d is the active payload, and the
    # Fortran-style lower/upper bounds (for traceability).
    is_int::Bool = false
    lb1::Int = 0
    ub1::Int = 0
    lb2::Int = 0
    ub2::Int = 0
end

# small logging helper — FatesGlobals.fates_log() returns a unit; we just print.
_restart_log_println(args...) = (println(stderr, args...); nothing)

# =====================================================================================

"""
    GetBounds(this, thread, dim_bounds, dim_kinds) -> (lb1, ub1, lb2, ub2)

Return the lower/upper bounds for this variable's payload. `thread==0` returns
the whole-proc bounds; otherwise the per-clump bounds for `thread` (1-based).
"""
function GetBounds(this::fates_restart_variable_type, thread::Integer,
                   dim_bounds::AbstractVector{fates_io_dimension_type},
                   dim_kinds::AbstractVector{fates_io_variable_kind_type})
    lb1 = 0; ub1 = 0; lb2 = 0; ub2 = 0

    ndims = dim_kinds[this.dim_kinds_index].ndims

    if thread == 0
        d_index = dim_kinds[this.dim_kinds_index].dim1_index
        lb1 = dim_bounds[d_index].lower_bound
        ub1 = dim_bounds[d_index].upper_bound
        if ndims > 1
            d_index = dim_kinds[this.dim_kinds_index].dim2_index
            lb2 = dim_bounds[d_index].lower_bound
            ub2 = dim_bounds[d_index].upper_bound
        end
    else
        d_index = dim_kinds[this.dim_kinds_index].dim1_index
        lb1 = dim_bounds[d_index].clump_lower_bound[thread]
        ub1 = dim_bounds[d_index].clump_upper_bound[thread]
        if ndims > 1
            d_index = dim_kinds[this.dim_kinds_index].dim2_index
            lb2 = dim_bounds[d_index].clump_lower_bound[thread]
            ub2 = dim_bounds[d_index].clump_upper_bound[thread]
        end
    end

    return (lb1, ub1, lb2, ub2)
end

# =====================================================================================

"""
    Init!(this, vname, units, long, vtype, flushval, num_dim_kinds, dim_kinds, dim_bounds)

Initialize a restart variable: store metadata, resolve the dim-kind index,
mark that kind active, then allocate EITHER the real or integer payload
(per `vtype`) filled with `flushval`.
"""
function Init!(this::fates_restart_variable_type, vname::AbstractString,
               units::AbstractString, long::AbstractString, vtype::AbstractString,
               flushval::Real, num_dim_kinds::Integer,
               dim_kinds::AbstractVector{fates_io_variable_kind_type},
               dim_bounds::AbstractVector{fates_io_dimension_type})
    this.vname = vname
    this.units = units
    this.long = long
    this.vtype = vtype
    this.flushval = Float64(flushval)

    # nullify both payloads
    this.r81d = Float64[]
    this.int1d = Int[]

    dk_index = iotype_index(strip(vtype), num_dim_kinds, dim_kinds)
    this.dim_kinds_index = dk_index
    set_active!(dim_kinds[dk_index])

    (lb1, ub1, lb2, ub2) = GetBounds(this, 0, dim_bounds, dim_kinds)
    this.lb1 = lb1; this.ub1 = ub1; this.lb2 = lb2; this.ub2 = ub2

    n = ub1 - lb1 + 1
    n < 0 && (n = 0)

    vt = strip(vtype)
    if vt == cohort_r8 || vt == site_r8
        this.is_int = false
        this.r81d = fill(Float64(flushval), n)
    elseif vt == cohort_int || vt == site_int
        this.is_int = true
        this.int1d = fill(round(Int, flushval), n)
    else
        _restart_log_println("Incompatible vtype passed to set_restart_var")
        _restart_log_println("vtype = ", strip(vtype), " ?")
        fates_endrun("FatesRestartVariableType.Init!: incompatible vtype $vtype")
    end

    return nothing
end

# =====================================================================================

"""
    Flush!(this, thread, dim_bounds, dim_kinds)

Reset this variable's payload (over its bounds) to `flushval`.
"""
function Flush!(this::fates_restart_variable_type, thread::Integer,
                dim_bounds::AbstractVector{fates_io_dimension_type},
                dim_kinds::AbstractVector{fates_io_variable_kind_type})
    (lb1, ub1, lb2, ub2) = GetBounds(this, thread, dim_bounds, dim_kinds)

    name = strip(dim_kinds[this.dim_kinds_index].name)
    # Translate Fortran lb1:ub1 (1-based, full-array on thread 0) into Julia
    # vector indices. For thread 0, fill the whole payload.
    if name == site_r8 || name == cohort_r8
        if thread == 0
            fill!(this.r81d, this.flushval)
        else
            @inbounds for i in lb1:ub1
                this.r81d[i] = this.flushval
            end
        end
    elseif name == site_int || name == cohort_int
        fv = round(Int, this.flushval)
        if thread == 0
            fill!(this.int1d, fv)
        else
            @inbounds for i in lb1:ub1
                this.int1d[i] = fv
            end
        end
    else
        _restart_log_println("fates restart variable type undefined while flushing restart variables")
        fates_endrun("FatesRestartVariableType.Flush!: undefined vtype")
    end

    return nothing
end
