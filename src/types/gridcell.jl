# ==========================================================================
# Ported from: src/main/GridcellType.F90
# Gridcell data type allocation and initialization
# ==========================================================================
#
# gridcell types can have values of:
#   1 => default

"""
    GridcellData

Gridcell-level data structure. Holds topological mapping, daylength, and
landunit index information for each gridcell.

Ported from `gridcell_type` in `GridcellType.F90`.
"""
Base.@kwdef mutable struct GridcellData{FT<:AbstractFloat}
    # topological mapping functionality, local 1d gdc arrays
    area     ::Vector{FT} = Float64[]   # total land area, gridcell (km^2)
    lat      ::Vector{FT} = Float64[]   # latitude (radians)
    lon      ::Vector{FT} = Float64[]   # longitude (radians)
    latdeg   ::Vector{FT} = Float64[]   # latitude (degrees)
    londeg   ::Vector{FT} = Float64[]   # longitude (degrees)
    active   ::Vector{Bool}    = Bool[]      # just needed for symmetry with other subgrid types

    nbedrock ::Vector{Int}     = Int[]       # index of uppermost bedrock layer

    # Daylength
    max_dayl  ::Vector{FT} = Float64[]  # maximum daylength for this grid cell (s)
    dayl      ::Vector{FT} = Float64[]  # daylength (seconds)
    prev_dayl ::Vector{FT} = Float64[]  # daylength from previous timestep (seconds)

    # indices into landunit-level arrays for landunits in this grid cell
    # (ispval implies this landunit doesn't exist on this grid cell)
    # [1:max_lunit, begg:endg]
    # (note: spatial dimension is last for efficiency — outer loop over g,
    #  inner loop over landunit type)
    landunit_indices ::Matrix{Int} = Matrix{Int}(undef, 0, 0)
end

"""
    gridcell_init!(grc::GridcellData, ngridcells::Int)

Allocate and initialize all fields of a `GridcellData` instance for
`ngridcells` gridcells.

Ported from `gridcell_type%Init` in `GridcellType.F90`.
"""
function gridcell_init!(grc::GridcellData{FT}, ngridcells::Int) where {FT}
    # The following is set in InitGridCells
    grc.area     = fill(FT(NaN), ngridcells)
    grc.lat      = fill(FT(NaN), ngridcells)
    grc.lon      = fill(FT(NaN), ngridcells)
    grc.latdeg   = fill(FT(NaN), ngridcells)
    grc.londeg   = fill(FT(NaN), ngridcells)
    grc.active   = fill(true, ngridcells)
    grc.nbedrock = fill(ISPVAL, ngridcells)

    # Initialized in module DayLength
    grc.max_dayl  = fill(FT(NaN), ngridcells)
    grc.dayl      = fill(FT(NaN), ngridcells)
    grc.prev_dayl = fill(FT(NaN), ngridcells)

    grc.landunit_indices = fill(ISPVAL, MAX_LUNIT, ngridcells)

    return nothing
end

"""
    gridcell_clean!(grc::GridcellData)

Deallocate (reset to empty) all fields of a `GridcellData` instance.

Ported from `gridcell_type%Clean` in `GridcellType.F90`.
"""
function gridcell_clean!(grc::GridcellData{FT}) where {FT}
    grc.area     = FT[]
    grc.lat      = FT[]
    grc.lon      = FT[]
    grc.latdeg   = FT[]
    grc.londeg   = FT[]
    grc.active   = Bool[]
    grc.nbedrock = Int[]
    grc.max_dayl  = FT[]
    grc.dayl      = FT[]
    grc.prev_dayl = FT[]
    grc.landunit_indices = Matrix{Int}(undef, 0, 0)

    return nothing
end
