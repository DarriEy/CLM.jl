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
Base.@kwdef mutable struct GridcellData
    # topological mapping functionality, local 1d gdc arrays
    area     ::Vector{Float64} = Float64[]   # total land area, gridcell (km^2)
    lat      ::Vector{Float64} = Float64[]   # latitude (radians)
    lon      ::Vector{Float64} = Float64[]   # longitude (radians)
    latdeg   ::Vector{Float64} = Float64[]   # latitude (degrees)
    londeg   ::Vector{Float64} = Float64[]   # longitude (degrees)
    active   ::Vector{Bool}    = Bool[]      # just needed for symmetry with other subgrid types

    nbedrock ::Vector{Int}     = Int[]       # index of uppermost bedrock layer

    # Daylength
    max_dayl  ::Vector{Float64} = Float64[]  # maximum daylength for this grid cell (s)
    dayl      ::Vector{Float64} = Float64[]  # daylength (seconds)
    prev_dayl ::Vector{Float64} = Float64[]  # daylength from previous timestep (seconds)

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
function gridcell_init!(grc::GridcellData, ngridcells::Int)
    # The following is set in InitGridCells
    grc.area     = fill(NaN, ngridcells)
    grc.lat      = fill(NaN, ngridcells)
    grc.lon      = fill(NaN, ngridcells)
    grc.latdeg   = fill(NaN, ngridcells)
    grc.londeg   = fill(NaN, ngridcells)
    grc.active   = fill(true, ngridcells)
    grc.nbedrock = fill(ISPVAL, ngridcells)

    # Initialized in module DayLength
    grc.max_dayl  = fill(NaN, ngridcells)
    grc.dayl      = fill(NaN, ngridcells)
    grc.prev_dayl = fill(NaN, ngridcells)

    grc.landunit_indices = fill(ISPVAL, MAX_LUNIT, ngridcells)

    return nothing
end

"""
    gridcell_clean!(grc::GridcellData)

Deallocate (reset to empty) all fields of a `GridcellData` instance.

Ported from `gridcell_type%Clean` in `GridcellType.F90`.
"""
function gridcell_clean!(grc::GridcellData)
    grc.area     = Float64[]
    grc.lat      = Float64[]
    grc.lon      = Float64[]
    grc.latdeg   = Float64[]
    grc.londeg   = Float64[]
    grc.active   = Bool[]
    grc.nbedrock = Int[]
    grc.max_dayl  = Float64[]
    grc.dayl      = Float64[]
    grc.prev_dayl = Float64[]
    grc.landunit_indices = Matrix{Int}(undef, 0, 0)

    return nothing
end
