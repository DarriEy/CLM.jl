# FatesHistoryVariableType.jl
# Julia port of FATES src/fates/main/FatesHistoryVariableType.F90
#
# The `fates_history_variable_type` derived type: one FATES history output
# variable, carrying its metadata (name/units/long-name/use_default/vtype/
# avgflag/upfreq/flushval) plus the bound real/integer data buffer. The Fortran
# uses six distinct allocatable pointers (r81d/r82d/r83d/int1d/int2d/int3d) and
# allocates exactly one per variable depending on its dimension-kind. We mirror
# that with six fields, only one of which is populated.
#
# Deps: FatesIODimensionsMod (fates_io_dimension_type), FatesIOVariableKindMod
#       (fates_io_variable_kind_type, iotype_index, set_active!, the site_*_r8
#       variable-kind name constants), FatesGlobals (fates_log/fates_endrun).
#
# Upstream-Fortran quirks preserved:
#   * The Init `select case(vtype)` only allocates r81d for `site_r8`; every
#     other listed kind is 2D (r82d). There is no native 3D registration path
#     here (r83d stays unbound) — matches Fortran.
#   * Bounds use Fortran lower:upper; in Julia we allocate a vector/matrix of
#     length (ub-lb+1) and remember `lb1`/`lb2` so site h_gid indices (which are
#     1-based and within [lb,ub]) map onto the buffer. For all FATES history
#     dims lb==1, so the offset is 0, but we keep the offset general.

# ---------------------------------------------------------------------------
# fates_history_variable_type
# ---------------------------------------------------------------------------
Base.@kwdef mutable struct fates_history_variable_type
    vname::String = ""
    units::String = ""
    long::String = ""
    use_default::String = ""   # "active" / "inactive"
    vtype::String = ""
    avgflag::String = ""       # single char, e.g. "A"
    upfreq::Int = 0            # update-frequency group (group_dyna_simple, ...)
    flushval::Float64 = 0.0
    dim_kinds_index::Int = 0
    # Data buffers — exactly one is allocated per variable (per its vtype).
    r81d::Vector{Float64} = Float64[]
    r82d::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    r83d::Array{Float64,3} = Array{Float64,3}(undef, 0, 0, 0)
    int1d::Vector{Int} = Int[]
    int2d::Matrix{Int} = Matrix{Int}(undef, 0, 0)
    int3d::Array{Int,3} = Array{Int,3}(undef, 0, 0, 0)
    # Stored lower bounds (Fortran lb1/lb2) so callers can map a 1-based site/
    # dim index onto the (possibly offset) buffer. allocated:: which buffer.
    lb1::Int = 0
    lb2::Int = 0
    allocated::Symbol = :none  # :r81d | :r82d | :r83d | :none
end

# The list of 2D variable kinds (everything except site_r8) handled by Init/HFlush.
const _HVAR_2D_KINDS = Set{String}([
    site_soil_r8, site_size_pft_r8, site_size_r8, site_coage_r8, site_coage_pft_r8,
    site_pft_r8, site_age_r8, site_height_r8, site_fuel_r8, site_cwdsc_r8,
    site_can_r8, site_cnlf_r8, site_cnlfpft_r8, site_cdsc_r8, site_cdpf_r8,
    site_scag_r8, site_scagpft_r8, site_agepft_r8, site_elem_r8, site_elpft_r8,
    site_elcwd_r8, site_elage_r8, site_agefuel_r8, site_landuse_r8, site_lulu_r8,
    site_clscpf_r8,
])

# =====================================================================================

"""
    GetBounds(this, thread, dim_bounds, dim_kinds) -> (lb1, ub1, lb2, ub2)

Return the (whole-proc when `thread==0`, else per-thread) lower/upper bounds for
the variable's dimension-kind. Mirrors the Fortran `GetBounds` type-bound method.
"""
function GetBounds(this::fates_history_variable_type, thread::Integer,
                   dim_bounds::AbstractVector{fates_io_dimension_type},
                   dim_kinds::AbstractVector{fates_io_variable_kind_type})
    lb1 = 0; ub1 = 0; lb2 = 0; ub2 = 0
    dk = dim_kinds[this.dim_kinds_index]
    ndims = dk.ndims
    if thread == 0
        d_index = dk.dim1_index
        lb1 = dim_bounds[d_index].lower_bound
        ub1 = dim_bounds[d_index].upper_bound
        if ndims > 1
            d_index = dk.dim2_index
            lb2 = dim_bounds[d_index].lower_bound
            ub2 = dim_bounds[d_index].upper_bound
        end
    else
        d_index = dk.dim1_index
        lb1 = dim_bounds[d_index].clump_lower_bound[thread]
        ub1 = dim_bounds[d_index].clump_upper_bound[thread]
        if ndims > 1
            d_index = dk.dim2_index
            lb2 = dim_bounds[d_index].clump_lower_bound[thread]
            ub2 = dim_bounds[d_index].clump_upper_bound[thread]
        end
    end
    return (lb1, ub1, lb2, ub2)
end

# =====================================================================================

"""
    Init!(this, vname, units, long, use_default, vtype, avgflag, flushval, upfreq,
          num_dim_kinds, dim_kinds, dim_bounds)

Initialize one history variable: copy metadata, resolve its dimension-kind
index, mark that kind active, look up its bounds, and allocate + flush the
single data buffer appropriate to `vtype`. Mirrors the Fortran `Init` method.
"""
function Init!(this::fates_history_variable_type,
               vname::AbstractString, units::AbstractString, long::AbstractString,
               use_default::AbstractString, vtype::AbstractString,
               avgflag::AbstractString, flushval::Real, upfreq::Integer,
               num_dim_kinds::Integer,
               dim_kinds::AbstractVector{fates_io_variable_kind_type},
               dim_bounds::AbstractVector{fates_io_dimension_type})

    this.vname = String(vname)
    this.units = String(units)
    this.long = String(long)
    this.use_default = String(use_default)
    this.vtype = String(vtype)
    this.avgflag = String(avgflag)
    this.flushval = Float64(flushval)
    this.upfreq = Int(upfreq)

    # nullify buffers (already empty by default kwdef)
    this.allocated = :none

    dk_index = iotype_index(strip(vtype), num_dim_kinds, dim_kinds)
    this.dim_kinds_index = dk_index
    set_active!(dim_kinds[dk_index])

    (lb1, ub1, lb2, ub2) = GetBounds(this, 0, dim_bounds, dim_kinds)
    this.lb1 = lb1
    this.lb2 = lb2

    vt = strip(vtype)
    if vt == site_r8
        n1 = ub1 - lb1 + 1
        this.r81d = fill(Float64(flushval), max(n1, 0))
        this.allocated = :r81d
    elseif vt in _HVAR_2D_KINDS
        n1 = ub1 - lb1 + 1
        n2 = ub2 - lb2 + 1
        this.r82d = fill(Float64(flushval), max(n1, 0), max(n2, 0))
        this.allocated = :r82d
    else
        @warn "Incompatible vtype passed to set_history_var: vtype = $(vt)"
        fates_endrun("FatesHistoryVariableType.Init!: incompatible vtype '$(vt)'")
    end
    return nothing
end

# =====================================================================================

"""
    HFlush!(this, thread, dim_bounds, dim_kinds)

Flush the variable's data buffer back to its `flushval`. Mirrors the Fortran
`HFlush` method (whole-buffer flush, not per-thread slab, since the Julia
buffers are proc-local arrays indexed from 1).
"""
function HFlush!(this::fates_history_variable_type, thread::Integer,
                 dim_bounds::AbstractVector{fates_io_dimension_type},
                 dim_kinds::AbstractVector{fates_io_variable_kind_type})
    name = strip(dim_kinds[this.dim_kinds_index].name)
    if name == site_r8
        fill!(this.r81d, this.flushval)
    elseif name in _HVAR_2D_KINDS
        fill!(this.r82d, this.flushval)
    else
        @warn "fates history variable type undefined while flushing history variables"
        fates_endrun("FatesHistoryVariableType.HFlush!: undefined kind '$(name)'")
    end
    return nothing
end
