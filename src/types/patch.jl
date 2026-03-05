# ==========================================================================
# Ported from: src/main/PatchType.F90
# Patch data type allocation and initialization
# ==========================================================================
#
# patch types can have values of:
#   0  => not_vegetated
#   1  => needleleaf_evergreen_temperate_tree
#   2  => needleleaf_evergreen_boreal_tree
#   3  => needleleaf_deciduous_boreal_tree
#   4  => broadleaf_evergreen_tropical_tree
#   5  => broadleaf_evergreen_temperate_tree
#   6  => broadleaf_deciduous_tropical_tree
#   7  => broadleaf_deciduous_temperate_tree
#   8  => broadleaf_deciduous_boreal_tree
#   9  => broadleaf_evergreen_shrub
#   10 => broadleaf_deciduous_temperate_shrub
#   11 => broadleaf_deciduous_boreal_shrub
#   12 => c3_arctic_grass
#   13 => c3_non-arctic_grass
#   14 => c4_grass
#   15 => c3_crop
#   16 => c3_irrigated
#   17 => temperate_corn
#   18 => irrigated_temperate_corn
#   19 => spring_wheat
#   20 => irrigated_spring_wheat
#   21 => winter_wheat
#   22 => irrigated_winter_wheat
#   23 => temperate_soybean
#   24 => irrigated_temperate_soybean
#   25 => barley
#   26 => irrigated_barley
#   27 => winter_barley
#   28 => irrigated_winter_barley
#   29 => rye
#   30 => irrigated_rye
#   31 => winter_rye
#   32 => irrigated_winter_rye
#   33 => cassava
#   34 => irrigated_cassava
#   35 => citrus
#   36 => irrigated_citrus
#   37 => cocoa
#   38 => irrigated_cocoa
#   39 => coffee
#   40 => irrigated_coffee
#   41 => cotton
#   42 => irrigated_cotton
#   43 => datepalm
#   44 => irrigated_datepalm
#   45 => foddergrass
#   46 => irrigated_foddergrass
#   47 => grapes
#   48 => irrigated_grapes
#   49 => groundnuts
#   50 => irrigated_groundnuts
#   51 => millet
#   52 => irrigated_millet
#   53 => oilpalm
#   54 => irrigated_oilpalm
#   55 => potatoes
#   56 => irrigated_potatoes
#   57 => pulses
#   58 => irrigated_pulses
#   59 => rapeseed
#   60 => irrigated_rapeseed
#   61 => rice
#   62 => irrigated_rice
#   63 => sorghum
#   64 => irrigated_sorghum
#   65 => sugarbeet
#   66 => irrigated_sugarbeet
#   67 => sugarcane
#   68 => irrigated_sugarcane
#   69 => sunflower
#   70 => irrigated_sunflower
#   71 => miscanthus
#   72 => irrigated_miscanthus
#   73 => switchgrass
#   74 => irrigated_switchgrass
#   75 => tropical_corn
#   76 => irrigated_tropical_corn
#   77 => tropical_soybean
#   78 => irrigated_tropical_soybean

"""
    PatchData

Patch-level data structure. Holds g/l/c/p hierarchy, vegetation type,
and FATES-specific fields for each patch.

Ported from `patch_type` in `PatchType.F90`.
"""
Base.@kwdef mutable struct PatchData{FT<:AbstractFloat}
    # g/l/c/p hierarchy, local g/l/c/p cells only
    column           ::Vector{Int}     = Int[]       # index into column level quantities
    wtcol            ::Vector{FT} = Float64[]   # weight (relative to column)
    landunit         ::Vector{Int}     = Int[]       # index into landunit level quantities
    wtlunit          ::Vector{FT} = Float64[]   # weight (relative to landunit)
    gridcell         ::Vector{Int}     = Int[]       # index into gridcell level quantities
    wtgcell          ::Vector{FT} = Float64[]   # weight (relative to gridcell)

    # Non-ED only
    itype            ::Vector{Int}     = Int[]       # patch vegetation type
    mxy              ::Vector{Int}     = Int[]       # m index for laixy(i,j,m),etc. (undefined for special landunits)
    active           ::Vector{Bool}    = Bool[]      # true => do computations on this patch

    # FATES only
    is_veg           ::Vector{Bool}    = Bool[]      # this is an ACTIVE fates patch
    is_bareground    ::Vector{Bool}    = Bool[]      # true => bareground fates patch
    wt_ed            ::Vector{FT} = Float64[]   # FATES patch weight
    sp_pftorder_index ::Vector{FT} = Float64[]  # index to map 'p' onto the order of ED patches in SP mode

    is_fates         ::Vector{Bool}    = Bool[]      # true for patch vector space reserved for FATES
end

"""
    patch_init!(pch::PatchData, npatches::Int)

Allocate and initialize all fields of a `PatchData` instance for
`npatches` patches. Integer fields are set to `ISPVAL`, real fields to `NaN`,
and logical fields to `false`.

FATES-only fields (`is_veg`, `is_bareground`, `wt_ed`, `sp_pftorder_index`)
are only allocated when `varctl.use_fates` is `true`.

Ported from `patch_type%Init` in `PatchType.F90`.
"""
function patch_init!(pch::PatchData{FT}, npatches::Int) where {FT}
    # --- g/l/c/p hierarchy (set in InitGridCellsMod) ---
    pch.gridcell    = fill(ISPVAL, npatches)
    pch.wtgcell     = fill(FT(NaN), npatches)

    pch.landunit    = fill(ISPVAL, npatches)
    pch.wtlunit     = fill(FT(NaN), npatches)

    pch.column      = fill(ISPVAL, npatches)
    pch.wtcol       = fill(FT(NaN), npatches)

    pch.mxy         = fill(ISPVAL, npatches)
    pch.active      = fill(false, npatches)

    pch.itype       = fill(ISPVAL, npatches)

    pch.is_fates    = fill(false, npatches)

    # FATES-only fields
    if varctl.use_fates
        pch.is_veg           = fill(false, npatches)
        pch.is_bareground    = fill(false, npatches)
        pch.wt_ed            = fill(FT(NaN), npatches)
        pch.sp_pftorder_index = fill(FT(NaN), npatches)
    end

    return nothing
end

"""
    patch_clean!(pch::PatchData)

Deallocate (reset to empty) all fields of a `PatchData` instance.

Ported from `patch_type%Clean` in `PatchType.F90`.
"""
function patch_clean!(pch::PatchData{FT}) where {FT}
    pch.gridcell    = Int[]
    pch.wtgcell     = FT[]
    pch.landunit    = Int[]
    pch.wtlunit     = FT[]
    pch.column      = Int[]
    pch.wtcol       = FT[]
    pch.itype       = Int[]
    pch.mxy         = Int[]
    pch.active      = Bool[]
    pch.is_fates    = Bool[]

    if varctl.use_fates
        pch.is_veg           = Bool[]
        pch.is_bareground    = Bool[]
        pch.wt_ed            = FT[]
        pch.sp_pftorder_index = FT[]
    end

    return nothing
end
