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
Base.@kwdef mutable struct LandunitData{FT<:Real}
    # g/l/c/p hierarchy, local g/l/c/p cells only
    gridcell    ::Vector{Int}     = Int[]       # index into gridcell level quantities
    wtgcell     ::Vector{FT} = Float64[]   # weight (relative to gridcell)
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
    canyon_hwr   ::Vector{FT} = Float64[]  # canyon height to width ratio (-)
    wtroad_perv  ::Vector{FT} = Float64[]  # weight of pervious road column to total road (-)
    wtlunit_roof ::Vector{FT} = Float64[]  # weight of roof with respect to urban landunit (-)
    ht_roof      ::Vector{FT} = Float64[]  # height of urban roof (m)
    z_0_town     ::Vector{FT} = Float64[]  # urban momentum roughness length (m)
    z_d_town     ::Vector{FT} = Float64[]  # urban displacement height (m)

    # hillslope variables
    stream_channel_depth  ::Vector{FT} = Float64[]  # stream channel bankfull depth (m)
    stream_channel_width  ::Vector{FT} = Float64[]  # stream channel bankfull width (m)
    stream_channel_length ::Vector{FT} = Float64[]  # stream channel length (m)
    stream_channel_slope  ::Vector{FT} = Float64[]  # stream channel slope (m/m)
    stream_channel_number ::Vector{FT} = Float64[]  # number of channels in landunit
end

"""
    landunit_init!(lun::LandunitData, nlandunits::Int)

Allocate and initialize all fields of a `LandunitData` instance for
`nlandunits` landunits. Integer fields are set to `ISPVAL`, real fields
to `NaN`, and logical fields to `false`.

Ported from `landunit_type%Init` in `LandunitType.F90`.
"""
function landunit_init!(lun::LandunitData{FT}, nlandunits::Int) where {FT}
    # g/l/c/p hierarchy — set in InitGridCellsMod
    lun.gridcell    = fill(ISPVAL, nlandunits)
    lun.wtgcell     = fill(FT(NaN), nlandunits)
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
    lun.canyon_hwr   = fill(FT(NaN), nlandunits)
    lun.wtroad_perv  = fill(FT(NaN), nlandunits)
    lun.wtlunit_roof = fill(FT(NaN), nlandunits)
    lun.ht_roof      = fill(FT(NaN), nlandunits)
    lun.z_0_town     = fill(FT(NaN), nlandunits)
    lun.z_d_town     = fill(FT(NaN), nlandunits)

    # hillslope variables — initialized in HillslopeHydrologyMod
    lun.stream_channel_depth  = fill(FT(NaN), nlandunits)
    lun.stream_channel_width  = fill(FT(NaN), nlandunits)
    lun.stream_channel_length = fill(FT(NaN), nlandunits)
    lun.stream_channel_slope  = fill(FT(NaN), nlandunits)
    lun.stream_channel_number = fill(FT(NaN), nlandunits)

    return nothing
end

"""
    landunit_clean!(lun::LandunitData)

Deallocate (reset to empty) all fields of a `LandunitData` instance.

Ported from `landunit_type%Clean` in `LandunitType.F90`.
"""
function landunit_clean!(lun::LandunitData{FT}) where {FT}
    lun.gridcell    = Int[]
    lun.wtgcell     = FT[]
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
    lun.canyon_hwr   = FT[]
    lun.wtroad_perv  = FT[]
    lun.wtlunit_roof = FT[]
    lun.ht_roof      = FT[]
    lun.z_0_town     = FT[]
    lun.z_d_town     = FT[]
    lun.stream_channel_depth  = FT[]
    lun.stream_channel_width  = FT[]
    lun.stream_channel_length = FT[]
    lun.stream_channel_slope  = FT[]
    lun.stream_channel_number = FT[]

    return nothing
end
