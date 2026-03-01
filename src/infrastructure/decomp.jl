# ==========================================================================
# Ported from: src/main/decompMod.F90
# Domain decomposition types, constants, and functions
# ==========================================================================

# --- Subgrid level constants ---
const SUBGRID_LEVEL_UNSPECIFIED = -1
const SUBGRID_LEVEL_LNDGRID    = 0
const SUBGRID_LEVEL_GRIDCELL   = 1
const SUBGRID_LEVEL_LANDUNIT   = 2
const SUBGRID_LEVEL_COLUMN     = 3
const SUBGRID_LEVEL_PATCH      = 4
const SUBGRID_LEVEL_COHORT     = 5

# --- Bounds level constants ---
const BOUNDS_LEVEL_PROC  = 1
const BOUNDS_LEVEL_CLUMP = 2

# =========================================================================
# Type definitions
# =========================================================================

"""
    BoundsType

Subgrid bounds for processor or clump decomposition.
Ported from `bounds_type` in `decompMod.F90`.
"""
Base.@kwdef mutable struct BoundsType
    begg::Int = 0         # beginning gridcell index
    endg::Int = 0         # ending gridcell index
    begl::Int = 0         # beginning landunit index
    endl::Int = 0         # ending landunit index
    begc::Int = 0         # beginning column index
    endc::Int = 0         # ending column index
    begp::Int = 0         # beginning patch index
    endp::Int = 0         # ending patch index
    begCohort::Int = 0    # beginning cohort index
    endCohort::Int = 0    # ending cohort index
    level::Int = 0        # whether defined on proc or clump level
    clump_index::Int = -1 # clump index (if clump level)
end

"""
    ProcessorType

Per-processor decomposition information.
Ported from `processor_type` in `decompMod.F90`.
"""
Base.@kwdef mutable struct ProcessorType
    nclumps::Int = 0          # number of clumps for this processor
    cid::Vector{Int} = Int[]  # clump indices
    ncells::Int = 0           # number of gridcells in proc
    nlunits::Int = 0          # number of landunits in proc
    ncols::Int = 0            # number of columns in proc
    npatches::Int = 0         # number of patches in proc
    nCohorts::Int = 0         # number of cohorts in proc
    begg::Int = 0             # beginning gridcell index
    endg::Int = 0             # ending gridcell index
    begl::Int = 0             # beginning landunit index
    endl::Int = 0             # ending landunit index
    begc::Int = 0             # beginning column index
    endc::Int = 0             # ending column index
    begp::Int = 0             # beginning patch index
    endp::Int = 0             # ending patch index
    begCohort::Int = 0        # beginning cohort index
    endCohort::Int = 0        # ending cohort index
end

"""
    ClumpType

Per-clump decomposition information.
Ported from `clump_type` in `decompMod.F90`.
"""
Base.@kwdef mutable struct ClumpType
    owner::Int = 0        # process id owning clump
    ncells::Int = 0       # number of gridcells in clump
    nlunits::Int = 0      # number of landunits in clump
    ncols::Int = 0        # number of columns in clump
    npatches::Int = 0     # number of patches in clump
    nCohorts::Int = 0     # number of cohorts in clump
    begg::Int = 0         # beginning gridcell index
    endg::Int = 0         # ending gridcell index
    begl::Int = 0         # beginning landunit index
    endl::Int = 0         # ending landunit index
    begc::Int = 0         # beginning column index
    endc::Int = 0         # ending column index
    begp::Int = 0         # beginning patch index
    endp::Int = 0         # ending patch index
    begCohort::Int = 0    # beginning cohort index
    endCohort::Int = 0    # ending cohort index
end

"""
    DecompData

Global decomposition state, wrapping all module-level variables from `decompMod.F90`.
"""
Base.@kwdef mutable struct DecompData
    procinfo::ProcessorType = ProcessorType()
    clumps::Vector{ClumpType} = ClumpType[]
    nclumps::Int = 0             # total number of clumps across all processors
    numg::Int = 0                # total number of gridcells on all procs
    numl::Int = 0                # total number of landunits on all procs
    numc::Int = 0                # total number of columns on all procs
    nump::Int = 0                # total number of patches on all procs
    numCohort::Int = 0           # total number of cohorts on all procs
    gindex_global::Vector{Int} = Int[]   # global index (includes ocean points)
    gindex_grc::Vector{Int} = Int[]      # gridcell global index
    gindex_lun::Vector{Int} = Int[]      # landunit global index
    gindex_col::Vector{Int} = Int[]      # column global index
    gindex_patch::Vector{Int} = Int[]    # patch global index
    gindex_cohort::Vector{Int} = Int[]   # cohort global index
    ldomain_ns::Int = 0          # domain size (from domainMod, set externally)
end

# Global decomposition state instance
const decomp = DecompData()

# =========================================================================
# Functions
# =========================================================================

"""
    get_beg(bounds::BoundsType, subgrid_level::Int) -> Int

Get beginning bound for a given subgrid level. Returns -1 for invalid level.
Ported from `get_beg` in `decompMod.F90`.
"""
function get_beg(bounds::BoundsType, subgrid_level::Int)
    if subgrid_level == SUBGRID_LEVEL_GRIDCELL
        return bounds.begg
    elseif subgrid_level == SUBGRID_LEVEL_LANDUNIT
        return bounds.begl
    elseif subgrid_level == SUBGRID_LEVEL_COLUMN
        return bounds.begc
    elseif subgrid_level == SUBGRID_LEVEL_PATCH
        return bounds.begp
    elseif subgrid_level == SUBGRID_LEVEL_COHORT
        return bounds.begCohort
    else
        return -1
    end
end

"""
    get_end(bounds::BoundsType, subgrid_level::Int) -> Int

Get ending bound for a given subgrid level. Returns -1 for invalid level.
Ported from `get_end` in `decompMod.F90`.
"""
function get_end(bounds::BoundsType, subgrid_level::Int)
    if subgrid_level == SUBGRID_LEVEL_GRIDCELL
        return bounds.endg
    elseif subgrid_level == SUBGRID_LEVEL_LANDUNIT
        return bounds.endl
    elseif subgrid_level == SUBGRID_LEVEL_COLUMN
        return bounds.endc
    elseif subgrid_level == SUBGRID_LEVEL_PATCH
        return bounds.endp
    elseif subgrid_level == SUBGRID_LEVEL_COHORT
        return bounds.endCohort
    else
        return -1
    end
end

"""
    get_clump_bounds(n::Int; decomp_data::DecompData=decomp) -> BoundsType

Determine clump bounds for processor clump index `n`.
Ported from `get_clump_bounds` in `decompMod.F90`.
"""
function get_clump_bounds(n::Int; decomp_data::DecompData=decomp)
    cid = decomp_data.procinfo.cid[n]
    clump = decomp_data.clumps[cid]
    pi = decomp_data.procinfo

    bounds = BoundsType()
    bounds.begp      = clump.begp      - pi.begp      + 1
    bounds.endp      = clump.endp      - pi.begp      + 1
    bounds.begc      = clump.begc      - pi.begc      + 1
    bounds.endc      = clump.endc      - pi.begc      + 1
    bounds.begl      = clump.begl      - pi.begl      + 1
    bounds.endl      = clump.endl      - pi.begl      + 1
    bounds.begg      = clump.begg      - pi.begg      + 1
    bounds.endg      = clump.endg      - pi.begg      + 1
    bounds.begCohort = clump.begCohort - pi.begCohort + 1
    bounds.endCohort = clump.endCohort - pi.begCohort + 1

    bounds.level = BOUNDS_LEVEL_CLUMP
    bounds.clump_index = n

    return bounds
end

"""
    get_proc_bounds(; decomp_data::DecompData=decomp) -> BoundsType

Retrieve processor bounds (1-based, relative to processor).
Ported from `get_proc_bounds` in `decompMod.F90`.
"""
function get_proc_bounds(; decomp_data::DecompData=decomp)
    pi = decomp_data.procinfo

    bounds = BoundsType()
    bounds.begp      = 1
    bounds.endp      = pi.endp      - pi.begp      + 1
    bounds.begc      = 1
    bounds.endc      = pi.endc      - pi.begc      + 1
    bounds.begl      = 1
    bounds.endl      = pi.endl      - pi.begl      + 1
    bounds.begg      = 1
    bounds.endg      = pi.endg      - pi.begg      + 1
    bounds.begCohort = 1
    bounds.endCohort = pi.endCohort - pi.begCohort + 1

    bounds.level = BOUNDS_LEVEL_PROC
    bounds.clump_index = -1

    return bounds
end

"""
    get_proc_total(pid::Int; decomp_data::DecompData=decomp)

Count up gridcells, landunits, columns, patches, and cohorts owned by process `pid`.
Returns a NamedTuple `(ncells, nlunits, ncols, npatches, nCohorts)`.
Ported from `get_proc_total` in `decompMod.F90`.
"""
function get_proc_total(pid::Int; decomp_data::DecompData=decomp)
    ncells = 0
    nlunits = 0
    ncols = 0
    npatches = 0
    nCohorts = 0
    for cid in 1:decomp_data.nclumps
        if decomp_data.clumps[cid].owner == pid
            ncells   += decomp_data.clumps[cid].ncells
            nlunits  += decomp_data.clumps[cid].nlunits
            ncols    += decomp_data.clumps[cid].ncols
            npatches += decomp_data.clumps[cid].npatches
            nCohorts += decomp_data.clumps[cid].nCohorts
        end
    end
    return (ncells=ncells, nlunits=nlunits, ncols=ncols, npatches=npatches, nCohorts=nCohorts)
end

"""
    get_proc_global(; decomp_data::DecompData=decomp)

Return global counts of gridcells, landunits, columns, patches, and cohorts.
Returns a NamedTuple `(ng, nl, nc, np, nCohorts)`.
Ported from `get_proc_global` in `decompMod.F90`.
"""
function get_proc_global(; decomp_data::DecompData=decomp)
    return (ng=decomp_data.numg, nl=decomp_data.numl, nc=decomp_data.numc,
            np=decomp_data.nump, nCohorts=decomp_data.numCohort)
end

"""
    get_proc_clumps(; decomp_data::DecompData=decomp) -> Int

Return the number of clumps for this processor.
Ported from `get_proc_clumps` in `decompMod.F90`.
"""
function get_proc_clumps(; decomp_data::DecompData=decomp)
    return decomp_data.procinfo.nclumps
end

"""
    get_global_index(subgrid_index::Int, subgrid_level::Int; decomp_data::DecompData=decomp) -> Int

Determine global index space value for target point at given subgrid level.
Ported from `get_global_index` in `decompMod.F90`.
"""
function get_global_index(subgrid_index::Int, subgrid_level::Int; decomp_data::DecompData=decomp)
    bounds_proc = get_proc_bounds(decomp_data=decomp_data)
    beg_index = get_beg(bounds_proc, subgrid_level)
    if beg_index == -1
        error("get_global_index: subgrid_level not supported: $subgrid_level")
    end

    gindex = get_subgrid_level_gindex(subgrid_level, decomp_data=decomp_data)
    return gindex[subgrid_index - beg_index + 1]
end

"""
    get_global_index_array(subgrid_index::AbstractVector{Int}, subgrid_level::Int; decomp_data::DecompData=decomp) -> Vector{Int}

Determine global index space values for target array at given subgrid level.
Ported from `get_global_index_array` in `decompMod.F90`.
"""
function get_global_index_array(subgrid_index::AbstractVector{Int}, subgrid_level::Int; decomp_data::DecompData=decomp)
    bounds_proc = get_proc_bounds(decomp_data=decomp_data)
    beg_index = get_beg(bounds_proc, subgrid_level)
    if beg_index == -1
        error("get_global_index_array: subgrid_level not supported: $subgrid_level")
    end

    gindex = get_subgrid_level_gindex(subgrid_level, decomp_data=decomp_data)
    result = similar(subgrid_index)
    for i in eachindex(subgrid_index)
        result[i] = gindex[subgrid_index[i] - beg_index + 1]
    end
    return result
end

"""
    get_subgrid_level_from_name(subgrid_level_name::AbstractString) -> Int

Given a name like "gridcell", return a subgrid level index.
Ported from `get_subgrid_level_from_name` in `decompMod.F90`.
"""
function get_subgrid_level_from_name(subgrid_level_name::AbstractString)
    if subgrid_level_name == GRLND
        return SUBGRID_LEVEL_LNDGRID
    elseif subgrid_level_name == NAMEG
        return SUBGRID_LEVEL_GRIDCELL
    elseif subgrid_level_name == NAMEL
        return SUBGRID_LEVEL_LANDUNIT
    elseif subgrid_level_name == NAMEC
        return SUBGRID_LEVEL_COLUMN
    elseif subgrid_level_name == NAMEP
        return SUBGRID_LEVEL_PATCH
    elseif subgrid_level_name == NAMECOHORT
        return SUBGRID_LEVEL_COHORT
    else
        error("get_subgrid_level_from_name: unknown subgrid_level_name: $subgrid_level_name")
    end
end

"""
    get_subgrid_level_gsize(subgrid_level::Int; decomp_data::DecompData=decomp) -> Int

Determine 1D global size from subgrid level.
Ported from `get_subgrid_level_gsize` in `decompMod.F90`.

Note: The `SUBGRID_LEVEL_LNDGRID` case uses `decomp_data.ldomain_ns` which
corresponds to `ldomain%ns` from `domainMod` in Fortran. Set `ldomain_ns`
in `DecompData` when the domain module is initialized.
"""
function get_subgrid_level_gsize(subgrid_level::Int; decomp_data::DecompData=decomp)
    if subgrid_level == SUBGRID_LEVEL_LNDGRID
        return decomp_data.ldomain_ns
    elseif subgrid_level == SUBGRID_LEVEL_GRIDCELL
        return decomp_data.numg
    elseif subgrid_level == SUBGRID_LEVEL_LANDUNIT
        return decomp_data.numl
    elseif subgrid_level == SUBGRID_LEVEL_COLUMN
        return decomp_data.numc
    elseif subgrid_level == SUBGRID_LEVEL_PATCH
        return decomp_data.nump
    elseif subgrid_level == SUBGRID_LEVEL_COHORT
        return decomp_data.numCohort
    else
        error("get_subgrid_level_gsize: unknown subgrid_level: $subgrid_level")
    end
end

"""
    get_subgrid_level_gindex(subgrid_level::Int; decomp_data::DecompData=decomp) -> Vector{Int}

Get global index array associated with subgrid level.
Ported from `get_subgrid_level_gindex` in `decompMod.F90`.
"""
function get_subgrid_level_gindex(subgrid_level::Int; decomp_data::DecompData=decomp)
    if subgrid_level == SUBGRID_LEVEL_LNDGRID
        return decomp_data.gindex_global
    elseif subgrid_level == SUBGRID_LEVEL_GRIDCELL
        return decomp_data.gindex_grc
    elseif subgrid_level == SUBGRID_LEVEL_LANDUNIT
        return decomp_data.gindex_lun
    elseif subgrid_level == SUBGRID_LEVEL_COLUMN
        return decomp_data.gindex_col
    elseif subgrid_level == SUBGRID_LEVEL_PATCH
        return decomp_data.gindex_patch
    elseif subgrid_level == SUBGRID_LEVEL_COHORT
        return decomp_data.gindex_cohort
    else
        error("get_subgrid_level_gindex: unknown subgrid_level: $subgrid_level")
    end
end
