# ==========================================================================
# Ported from: src/main/LandunitType.F90
# Landunit data type allocation and initialization
# ==========================================================================
#
# landunit types can have values of (see landunit_varcon.jl):
#   1 => (istsoil)    soil (vegetated or bare soil landunit)
#   2 => (istcrop)    crop (only for crop configuration)
#   3 => (UNUSED)     (formerly non-multiple elevation class land ice)
#   4 => (istice)     land ice
#   5 => (istdlak)    deep lake
#   6 => (istwet)     wetland
#   7 => (isturb_tbd) urban tbd
#   8 => (isturb_hd)  urban hd
#   9 => (isturb_md)  urban md

"""
    LandunitData

Landunit-level data structure. Holds g/l/c/p hierarchy, topological mapping,
urban properties, and hillslope variables for each landunit.

Ported from `landunit_type` in `LandunitType.F90`.
"""
Base.@kwdef mutable struct LandunitData
    # g/l/c/p hierarchy, local g/l/c/p cells only
    gridcell    ::Vector{Int}     = Int[]       # index into gridcell level quantities
    wtgcell     ::Vector{Float64} = Float64[]   # weight (relative to gridcell)
    coli        ::Vector{Int}     = Int[]       # beginning column index per landunit
    colf        ::Vector{Int}     = Int[]       # ending column index for each landunit
    ncolumns    ::Vector{Int}     = Int[]       # number of columns for each landunit
    nhillslopes ::Vector{Int}     = Int[]       # number of hillslopes for each landunit
    patchi      ::Vector{Int}     = Int[]       # beginning patch index for each landunit
    patchf      ::Vector{Int}     = Int[]       # ending patch index for each landunit
    npatches    ::Vector{Int}     = Int[]       # number of patches for each landunit

    # topological mapping functionality
    itype       ::Vector{Int}     = Int[]       # landunit type
    ifspecial   ::Vector{Bool}    = Bool[]      # true => landunit is not vegetated
    lakpoi      ::Vector{Bool}    = Bool[]      # true => lake point
    urbpoi      ::Vector{Bool}    = Bool[]      # true => urban point
    glcpoi      ::Vector{Bool}    = Bool[]      # true => glacier point
    active      ::Vector{Bool}    = Bool[]      # true => do computations on this landunit

    # urban properties
    canyon_hwr   ::Vector{Float64} = Float64[]  # canyon height to width ratio (-)
    wtroad_perv  ::Vector{Float64} = Float64[]  # weight of pervious road column to total road (-)
    wtlunit_roof ::Vector{Float64} = Float64[]  # weight of roof with respect to urban landunit (-)
    ht_roof      ::Vector{Float64} = Float64[]  # height of urban roof (m)
    z_0_town     ::Vector{Float64} = Float64[]  # urban momentum roughness length (m)
    z_d_town     ::Vector{Float64} = Float64[]  # urban displacement height (m)

    # hillslope variables
    stream_channel_depth  ::Vector{Float64} = Float64[]  # stream channel bankfull depth (m)
    stream_channel_width  ::Vector{Float64} = Float64[]  # stream channel bankfull width (m)
    stream_channel_length ::Vector{Float64} = Float64[]  # stream channel length (m)
    stream_channel_slope  ::Vector{Float64} = Float64[]  # stream channel slope (m/m)
    stream_channel_number ::Vector{Float64} = Float64[]  # number of channels in landunit
end

"""
    landunit_init!(lun::LandunitData, nlandunits::Int)

Allocate and initialize all fields of a `LandunitData` instance for
`nlandunits` landunits. Integer fields are set to `ISPVAL`, real fields
to `NaN`, and logical fields to `false`.

Ported from `landunit_type%Init` in `LandunitType.F90`.
"""
function landunit_init!(lun::LandunitData, nlandunits::Int)
    # g/l/c/p hierarchy — set in InitGridCellsMod
    lun.gridcell    = fill(ISPVAL, nlandunits)
    lun.wtgcell     = fill(NaN, nlandunits)
    lun.coli        = fill(ISPVAL, nlandunits)
    lun.colf        = fill(ISPVAL, nlandunits)
    lun.ncolumns    = fill(ISPVAL, nlandunits)
    lun.nhillslopes = fill(ISPVAL, nlandunits)
    lun.patchi      = fill(ISPVAL, nlandunits)
    lun.patchf      = fill(ISPVAL, nlandunits)
    lun.npatches    = fill(ISPVAL, nlandunits)

    # topological mapping
    lun.itype     = fill(ISPVAL, nlandunits)
    lun.ifspecial = fill(false, nlandunits)
    lun.lakpoi    = fill(false, nlandunits)
    lun.urbpoi    = fill(false, nlandunits)
    lun.glcpoi    = fill(false, nlandunits)

    # active — initialized in routine setActive in module reweightMod
    lun.active = Vector{Bool}(undef, nlandunits)

    # urban properties — set in urbanparams_inst%Init in UrbanParamsType
    lun.canyon_hwr   = fill(NaN, nlandunits)
    lun.wtroad_perv  = fill(NaN, nlandunits)
    lun.wtlunit_roof = fill(NaN, nlandunits)
    lun.ht_roof      = fill(NaN, nlandunits)
    lun.z_0_town     = fill(NaN, nlandunits)
    lun.z_d_town     = fill(NaN, nlandunits)

    # hillslope variables — initialized in HillslopeHydrologyMod
    lun.stream_channel_depth  = fill(NaN, nlandunits)
    lun.stream_channel_width  = fill(NaN, nlandunits)
    lun.stream_channel_length = fill(NaN, nlandunits)
    lun.stream_channel_slope  = fill(NaN, nlandunits)
    lun.stream_channel_number = fill(NaN, nlandunits)

    return nothing
end

"""
    landunit_clean!(lun::LandunitData)

Deallocate (reset to empty) all fields of a `LandunitData` instance.

Ported from `landunit_type%Clean` in `LandunitType.F90`.
"""
function landunit_clean!(lun::LandunitData)
    lun.gridcell    = Int[]
    lun.wtgcell     = Float64[]
    lun.coli        = Int[]
    lun.colf        = Int[]
    lun.ncolumns    = Int[]
    lun.nhillslopes = Int[]
    lun.patchi      = Int[]
    lun.patchf      = Int[]
    lun.npatches    = Int[]
    lun.itype       = Int[]
    lun.ifspecial   = Bool[]
    lun.lakpoi      = Bool[]
    lun.urbpoi      = Bool[]
    lun.glcpoi      = Bool[]
    lun.active      = Bool[]
    lun.canyon_hwr   = Float64[]
    lun.wtroad_perv  = Float64[]
    lun.wtlunit_roof = Float64[]
    lun.ht_roof      = Float64[]
    lun.z_0_town     = Float64[]
    lun.z_d_town     = Float64[]
    lun.stream_channel_depth  = Float64[]
    lun.stream_channel_width  = Float64[]
    lun.stream_channel_length = Float64[]
    lun.stream_channel_slope  = Float64[]
    lun.stream_channel_number = Float64[]

    return nothing
end
