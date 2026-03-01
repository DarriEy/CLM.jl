# ==========================================================================
# Ported from: src/main/ColumnType.F90
# Column data type allocation and initialization
# ==========================================================================
#
# column types can have values of:
#   1  => (istsoil)          soil (vegetated or bare soil)
#   2  => (istcrop)          crop (only for crop configuration)
#   3  => (UNUSED)           (formerly non-multiple elevation class land ice)
#   4  => (istice)           land ice
#   5  => (istdlak)          deep lake
#   6  => (istwet)           wetland
#   71 => (icol_roof)        urban roof
#   72 => (icol_sunwall)     urban sunwall
#   73 => (icol_shadewall)   urban shadewall
#   74 => (icol_road_imperv) urban impervious road
#   75 => (icol_road_perv)   urban pervious road

"""
    ColumnData

Column-level data structure. Holds g/l/c/p hierarchy, topological mapping,
topography, vertical levels, hillslope hydrology, and other column
characteristics.

Ported from `column_type` in `ColumnType.F90`.
"""
Base.@kwdef mutable struct ColumnData
    # g/l/c/p hierarchy, local g/l/c/p cells only
    landunit         ::Vector{Int}     = Int[]       # index into landunit level quantities
    wtlunit          ::Vector{Float64} = Float64[]   # weight (relative to landunit)
    gridcell         ::Vector{Int}     = Int[]       # index into gridcell level quantities
    wtgcell          ::Vector{Float64} = Float64[]   # weight (relative to gridcell)
    patchi           ::Vector{Int}     = Int[]       # beginning patch index for each column
    patchf           ::Vector{Int}     = Int[]       # ending patch index for each column
    npatches         ::Vector{Int}     = Int[]       # number of patches for each column

    # topological mapping functionality
    itype            ::Vector{Int}     = Int[]       # column type
    lun_itype        ::Vector{Int}     = Int[]       # landunit type (convenience copy)
    active           ::Vector{Bool}    = Bool[]      # true => do computations on this column
    type_is_dynamic  ::Vector{Bool}    = Bool[]      # true => itype can change throughout the run

    is_fates         ::Vector{Bool}    = Bool[]      # true => this is a fates column

    # topography
    micro_sigma      ::Vector{Float64} = Float64[]   # microtopography pdf sigma (m)
    topo_slope       ::Vector{Float64} = Float64[]   # gridcell topographic slope
    topo_std         ::Vector{Float64} = Float64[]   # gridcell elevation standard deviation

    # vertical levels
    snl              ::Vector{Int}     = Int[]       # number of snow layers
    dz               ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # layer thickness (m)
    z                ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # layer depth (m)
    zi               ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # interface level below a "z" level (m)
    zii              ::Vector{Float64} = Float64[]   # convective boundary height (m)
    dz_lake          ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # lake layer thickness (m)
    z_lake           ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # layer depth for lake (m)
    lakedepth        ::Vector{Float64} = Float64[]   # variable lake depth (m)
    nbedrock         ::Vector{Int}     = Int[]       # variable depth to bedrock index

    # hillslope hydrology variables
    col_ndx          ::Vector{Int}     = Int[]       # column index of column
    colu             ::Vector{Int}     = Int[]       # column index of uphill column
    cold             ::Vector{Int}     = Int[]       # column index of downhill column
    hillslope_ndx    ::Vector{Int}     = Int[]       # hillslope identifier
    hill_elev        ::Vector{Float64} = Float64[]   # mean elevation relative to stream channel (m)
    hill_slope       ::Vector{Float64} = Float64[]   # mean along-hill slope (m/m)
    hill_area        ::Vector{Float64} = Float64[]   # mean surface area (m2)
    hill_width       ::Vector{Float64} = Float64[]   # across-hill width of bottom boundary (m)
    hill_distance    ::Vector{Float64} = Float64[]   # along-hill distance from bottom (m)
    hill_aspect      ::Vector{Float64} = Float64[]   # azimuth angle wrt north (radians)

    # other column characteristics
    is_hillslope_column   ::Vector{Bool}    = Bool[]  # true if hillslope element
    hydrologically_active ::Vector{Bool}    = Bool[]  # true if hydrologically active type
    urbpoi                ::Vector{Bool}    = Bool[]  # true => urban point

    # levgrnd_class: class in which each layer falls (e.g., soil vs. bedrock).
    # ispval indicates the layer is completely unused for this column.
    levgrnd_class    ::Matrix{Int}     = Matrix{Int}(undef, 0, 0)  # (ncols, nlevmaxurbgrnd)
end

"""
    column_init!(col::ColumnData, ncols::Int)

Allocate and initialize all fields of a `ColumnData` instance for
`ncols` columns. Integer fields are set to `ISPVAL`, real fields to `NaN`
(or `SPVAL` where Fortran uses `spval`), and logical fields to `false`.

Ported from `column_type%Init` in `ColumnType.F90`.
"""
function column_init!(col::ColumnData, ncols::Int)
    nlevsno        = varpar.nlevsno
    nlevmaxurbgrnd = varpar.nlevmaxurbgrnd
    nlevlak        = varpar.nlevlak

    # Vertical dimension sizes (Fortran offset indices → Julia 1-based)
    # dz, z: -nlevsno+1 : nlevmaxurbgrnd  →  nlevsno + nlevmaxurbgrnd levels
    # zi:    -nlevsno   : nlevmaxurbgrnd  →  nlevsno + nlevmaxurbgrnd + 1 levels
    nlev_dz = nlevsno + nlevmaxurbgrnd
    nlev_zi = nlevsno + nlevmaxurbgrnd + 1

    # --- g/l/c/p hierarchy (set in initGridCellsMod) ---
    col.gridcell    = fill(ISPVAL, ncols)
    col.wtgcell     = fill(NaN, ncols)
    col.landunit    = fill(ISPVAL, ncols)
    col.wtlunit     = fill(NaN, ncols)
    col.patchi      = fill(ISPVAL, ncols)
    col.patchf      = fill(ISPVAL, ncols)
    col.npatches    = fill(ISPVAL, ncols)
    col.itype       = fill(ISPVAL, ncols)
    col.lun_itype   = fill(ISPVAL, ncols)
    col.active      = fill(false, ncols)
    col.type_is_dynamic = fill(false, ncols)

    col.is_fates    = fill(false, ncols)

    # --- vertical levels (set in initVerticalMod) ---
    col.snl         = fill(ISPVAL, ncols)
    col.dz          = fill(NaN, ncols, nlev_dz)
    col.z           = fill(NaN, ncols, nlev_dz)
    col.zi          = fill(NaN, ncols, nlev_zi)
    col.zii         = fill(NaN, ncols)
    col.lakedepth   = fill(SPVAL, ncols)
    col.dz_lake     = fill(NaN, ncols, nlevlak)
    col.z_lake      = fill(NaN, ncols, nlevlak)
    col.col_ndx     = fill(ISPVAL, ncols)
    col.colu        = fill(ISPVAL, ncols)
    col.cold        = fill(ISPVAL, ncols)
    col.hillslope_ndx = fill(ISPVAL, ncols)
    col.hill_elev   = fill(SPVAL, ncols)
    col.hill_slope  = fill(SPVAL, ncols)
    col.hill_area   = fill(SPVAL, ncols)
    col.hill_width  = fill(SPVAL, ncols)
    col.hill_distance = fill(SPVAL, ncols)
    col.hill_aspect = fill(SPVAL, ncols)
    col.nbedrock    = fill(ISPVAL, ncols)
    col.levgrnd_class = fill(ISPVAL, ncols, nlevmaxurbgrnd)
    col.micro_sigma = fill(NaN, ncols)
    col.topo_slope  = fill(NaN, ncols)
    col.topo_std    = fill(NaN, ncols)
    col.is_hillslope_column   = fill(false, ncols)
    col.hydrologically_active = fill(false, ncols)
    col.urbpoi      = fill(false, ncols)

    return nothing
end

"""
    column_clean!(col::ColumnData)

Deallocate (reset to empty) all fields of a `ColumnData` instance.

Ported from `column_type%Clean` in `ColumnType.F90`.
"""
function column_clean!(col::ColumnData)
    col.gridcell    = Int[]
    col.wtgcell     = Float64[]
    col.landunit    = Int[]
    col.wtlunit     = Float64[]
    col.patchi      = Int[]
    col.patchf      = Int[]
    col.npatches    = Int[]
    col.itype       = Int[]
    col.lun_itype   = Int[]
    col.active      = Bool[]
    col.is_fates    = Bool[]
    col.type_is_dynamic = Bool[]
    col.snl         = Int[]
    col.dz          = Matrix{Float64}(undef, 0, 0)
    col.z           = Matrix{Float64}(undef, 0, 0)
    col.zi          = Matrix{Float64}(undef, 0, 0)
    col.zii         = Float64[]
    col.lakedepth   = Float64[]
    col.dz_lake     = Matrix{Float64}(undef, 0, 0)
    col.z_lake      = Matrix{Float64}(undef, 0, 0)
    col.micro_sigma = Float64[]
    col.topo_slope  = Float64[]
    col.topo_std    = Float64[]
    col.nbedrock    = Int[]
    col.levgrnd_class = Matrix{Int}(undef, 0, 0)
    col.is_hillslope_column   = Bool[]
    col.hydrologically_active = Bool[]
    col.col_ndx     = Int[]
    col.colu        = Int[]
    col.cold        = Int[]
    col.hillslope_ndx = Int[]
    col.hill_elev   = Float64[]
    col.hill_slope  = Float64[]
    col.hill_area   = Float64[]
    col.hill_width  = Float64[]
    col.hill_distance = Float64[]
    col.hill_aspect = Float64[]
    col.urbpoi      = Bool[]

    return nothing
end

"""
    column_update_itype!(col::ColumnData, c::Int, itype::Int)

Update the column type for column `c`. Any updates to `col.itype` after
initialization should be made via this function.

This can NOT be used to change the landunit type: it can only be used to
change the column type within a fixed landunit.

Ported from `column_type%update_itype` in `ColumnType.F90`.
"""
function column_update_itype!(col::ColumnData, c::Int, itype::Int)
    if col.type_is_dynamic[c]
        col.itype[c] = itype
        col.hydrologically_active[c] = is_hydrologically_active(
            itype, col.lun_itype[c])
    else
        error("column_update_itype! ERROR: attempt to update itype when " *
              "type_is_dynamic is false (c=$c, col.itype[c]=$(col.itype[c]), itype=$itype)")
    end
    return nothing
end
