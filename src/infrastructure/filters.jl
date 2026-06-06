# ==========================================================================
# Ported from: src/main/filterMod.F90
# Filter to mask conversion — replaces Fortran integer filter arrays
# with BitVector masks for GPU-compatible, allocation-free iteration.
# ==========================================================================

# npcropmin: minimum PFT index for prognostic crops
# (from pftconMod.F90, not yet ported — defined here as constant)
const NPCROPMIN = 17

"""
    ClumpFilter

Filter masks for processing columns, patches, and landunits of particular
types, including lake, non-lake, urban, soil, snow, non-snow, and
naturally-vegetated patches.

Replaces Fortran `clumpfilter` type: integer filter arrays + counts become
BitVectors. Usage:

    for c in bounds.begc:bounds.endc
        filt.soilc[c] || continue
        # physics for column c
    end

Ported from `clumpfilter` in `filterMod.F90`.
"""
Base.@kwdef mutable struct ClumpFilter{B<:AbstractVector{Bool}}
    # Column-level masks
    allc::B = BitVector()                    # all columns
    lakec::B = BitVector()          # lake filter (columns)
    nolakec::B = BitVector()        # non-lake filter (columns)
    soilc::B = BitVector()          # soil filter (columns)
    bgc_soilc::B = BitVector()      # soil with biogeochemistry active (columns)
    snowc::B = BitVector()          # snow filter (columns) — set in SnowHydrology
    nosnowc::B = BitVector()        # non-snow filter (columns) — set in SnowHydrology
    lakesnowc::B = BitVector()      # lake snow filter (columns) — set in LakeHydrology
    lakenosnowc::B = BitVector()    # lake non-snow filter (columns) — set in LakeHydrology
    hydrologyc::B = BitVector()     # hydrology filter (columns)
    urbanc::B = BitVector()         # urban filter (columns)
    nourbanc::B = BitVector()       # non-urban filter (columns)
    icec::B = BitVector()           # glacier filter (columns)
    do_smb_c::B = BitVector()       # glacier+bareland SMB filter (columns)
    actfirec::B = BitVector()       # active fire filter (columns) — set elsewhere

    # Patch-level masks
    natvegp::B = BitVector()        # CNDV nat-vegetated filter (patches) — set by CNDV
    pcropp::B = BitVector()         # prognostic crop filter (patches)
    soilnopcropp::B = BitVector()   # soil w/o prog. crops (patches)
    all_soil_patches::B = BitVector() # all soil or crop patches (for FATES SP)
    lakep::B = BitVector()          # lake filter (patches)
    nolakep::B = BitVector()        # non-lake filter (patches)
    bgc_vegp::B = BitVector()       # vegetation biochemistry active (patches)
    soilp::B = BitVector()          # soil filter (patches)
    exposedvegp::B = BitVector()    # exposed vegetation (patches) — set by set_exposedvegp_filter!
    noexposedvegp::B = BitVector()  # non-exposed vegetation (patches)
    urbanp::B = BitVector()         # urban filter (patches)
    nourbanp::B = BitVector()       # non-urban filter (patches)
    nolakeurbanp::B = BitVector()   # non-lake, non-urban filter (patches)
    actfirep::B = BitVector()       # active fire filter (patches) — set elsewhere

    # Landunit-level masks
    urbanl::B = BitVector()         # urban filter (landunits)
    nourbanl::B = BitVector()       # non-urban filter (landunits)
end

# Device-movable: adapt swaps every BitVector mask to the target backend array
# type (e.g. MtlVector{Bool}) so the whole filter set can ride along to the GPU.
Adapt.@adapt_structure ClumpFilter

# Global filter instances
# `clump_filter` is the standard set of filters (active points only).
const clump_filter = ClumpFilter()

# `clump_filter_inactive_and_active` includes both inactive and active points.
# Rarely appropriate to use — see Fortran comments in filterMod.F90.
const clump_filter_inactive_and_active = ClumpFilter()

# =========================================================================
# Allocation
# =========================================================================

"""
    alloc_filters!(filt::ClumpFilter, ncols::Int, npatches::Int, nlandunits::Int)

Allocate BitVector masks for a `ClumpFilter` instance.
All masks are initialized to `false`.

Ported from `allocFiltersOneGroup` in `filterMod.F90`.
"""
function alloc_filters!(filt::ClumpFilter, ncols::Int, npatches::Int, nlandunits::Int)
    # Column-level masks
    filt.allc = falses(ncols)
    filt.lakec = falses(ncols)
    filt.nolakec = falses(ncols)
    filt.soilc = falses(ncols)
    filt.bgc_soilc = falses(ncols)
    filt.snowc = falses(ncols)
    filt.nosnowc = falses(ncols)
    filt.lakesnowc = falses(ncols)
    filt.lakenosnowc = falses(ncols)
    filt.hydrologyc = falses(ncols)
    filt.urbanc = falses(ncols)
    filt.nourbanc = falses(ncols)
    filt.icec = falses(ncols)
    filt.do_smb_c = falses(ncols)
    filt.actfirec = falses(ncols)

    # Patch-level masks
    filt.natvegp = falses(npatches)
    filt.pcropp = falses(npatches)
    filt.soilnopcropp = falses(npatches)
    filt.all_soil_patches = falses(npatches)
    filt.lakep = falses(npatches)
    filt.nolakep = falses(npatches)
    filt.bgc_vegp = falses(npatches)
    filt.soilp = falses(npatches)
    filt.exposedvegp = falses(npatches)
    filt.noexposedvegp = falses(npatches)
    filt.urbanp = falses(npatches)
    filt.nourbanp = falses(npatches)
    filt.nolakeurbanp = falses(npatches)
    filt.actfirep = falses(npatches)

    # Landunit-level masks
    filt.urbanl = falses(nlandunits)
    filt.nourbanl = falses(nlandunits)

    return nothing
end

"""
    alloc_all_filters!(ncols::Int, npatches::Int, nlandunits::Int)

Allocate both the standard (active-only) and inactive+active filter sets.

Ported from `allocFilters` in `filterMod.F90`.
"""
function alloc_all_filters!(ncols::Int, npatches::Int, nlandunits::Int)
    alloc_filters!(clump_filter, ncols, npatches, nlandunits)
    alloc_filters!(clump_filter_inactive_and_active, ncols, npatches, nlandunits)
    return nothing
end

# =========================================================================
# Setting filters
# =========================================================================

"""
    set_filters!(bounds::BoundsType,
                 col::ColumnData, lun::LandunitData,
                 pch::PatchData, grc::GridcellData;
                 melt_replaced_by_ice::Vector{Bool}=Bool[])

Set both active-only and inactive+active filter masks.

Ported from `setFilters` in `filterMod.F90`.
"""
function set_filters!(bounds::BoundsType,
                      col::ColumnData, lun::LandunitData,
                      pch::PatchData, grc::GridcellData;
                      melt_replaced_by_ice::Vector{Bool}=Bool[])
    @assert bounds.level == BOUNDS_LEVEL_CLUMP

    set_filters_one_group!(clump_filter, bounds, false,
                           col, lun, pch, grc;
                           melt_replaced_by_ice=melt_replaced_by_ice)

    set_filters_one_group!(clump_filter_inactive_and_active, bounds, true,
                           col, lun, pch, grc;
                           melt_replaced_by_ice=melt_replaced_by_ice)

    return nothing
end

"""
    set_filters_one_group!(filt::ClumpFilter, bounds::BoundsType,
                           include_inactive::Bool,
                           col::ColumnData, lun::LandunitData,
                           pch::PatchData, grc::GridcellData;
                           melt_replaced_by_ice::Vector{Bool}=Bool[])

Set filter masks for one group (active-only or inactive+active).

"Standard" filters only include active points. Setting `include_inactive=true`
creates alternative filters that also apply over inactive points.

This routine sets filters determined by subgrid type and "active" status.
Filters based on model state (e.g., snow cover) should be set elsewhere.

Ported from `setFiltersOneGroup` in `filterMod.F90`.
"""
function set_filters_one_group!(filt::ClumpFilter, bounds::BoundsType,
                                include_inactive::Bool,
                                col::ColumnData, lun::LandunitData,
                                pch::PatchData, grc::GridcellData;
                                melt_replaced_by_ice::Vector{Bool}=Bool[])
    @assert bounds.level == BOUNDS_LEVEL_CLUMP

    # Clear column-level masks in bounds range
    for c in bounds.begc:bounds.endc
        filt.allc[c] = false
        filt.lakec[c] = false
        filt.nolakec[c] = false
        filt.soilc[c] = false
        filt.bgc_soilc[c] = false
        filt.hydrologyc[c] = false
        filt.urbanc[c] = false
        filt.nourbanc[c] = false
        filt.icec[c] = false
        filt.do_smb_c[c] = false
    end

    # Clear patch-level masks in bounds range
    for p in bounds.begp:bounds.endp
        filt.pcropp[p] = false
        filt.soilnopcropp[p] = false
        filt.lakep[p] = false
        filt.nolakep[p] = false
        filt.bgc_vegp[p] = false
        filt.soilp[p] = false
        filt.urbanp[p] = false
        filt.nourbanp[p] = false
        filt.nolakeurbanp[p] = false
    end

    # Clear landunit-level masks in bounds range
    for l in bounds.begl:bounds.endl
        filt.urbanl[l] = false
        filt.nourbanl[l] = false
    end

    # --- All columns filter ---
    for c in bounds.begc:bounds.endc
        if col.active[c] || include_inactive
            filt.allc[c] = true
        end
    end

    # --- Lake / non-lake column filter ---
    for c in bounds.begc:bounds.endc
        if col.active[c] || include_inactive
            l = col.landunit[c]
            if lun.lakpoi[l]
                filt.lakec[c] = true
            else
                filt.nolakec[c] = true
            end
        end
    end

    # --- Lake / non-lake / nolakeurban patch filter ---
    for p in bounds.begp:bounds.endp
        if pch.active[p] || include_inactive
            l = pch.landunit[p]
            if lun.lakpoi[l]
                filt.lakep[p] = true
            else
                filt.nolakep[p] = true
                if !lun.urbpoi[l]
                    filt.nolakeurbanp[p] = true
                end
            end
        end
    end

    # --- BGC soil column filter (use_cn or use_fates_bgc) ---
    if varctl.use_cn || varctl.use_fates_bgc
        for c in bounds.begc:bounds.endc
            if col.active[c] || include_inactive
                l = col.landunit[c]
                if lun.itype[l] == ISTSOIL || lun.itype[l] == ISTCROP
                    filt.bgc_soilc[c] = true
                end
            end
        end
    end

    # --- BGC veg patch filter (use_cn only, NOT FATES) ---
    if varctl.use_cn
        for p in bounds.begp:bounds.endp
            if pch.active[p] || include_inactive
                l = pch.landunit[p]
                if lun.itype[l] == ISTSOIL || lun.itype[l] == ISTCROP
                    filt.bgc_vegp[p] = true
                end
            end
        end
    end

    # --- Soil column filter ---
    for c in bounds.begc:bounds.endc
        if col.active[c] || include_inactive
            l = col.landunit[c]
            if lun.itype[l] == ISTSOIL || lun.itype[l] == ISTCROP
                filt.soilc[c] = true
            end
        end
    end

    # --- Soil patch filter ---
    for p in bounds.begp:bounds.endp
        if pch.active[p] || include_inactive
            l = pch.landunit[p]
            if lun.itype[l] == ISTSOIL || lun.itype[l] == ISTCROP
                filt.soilp[p] = true
            end
        end
    end

    # --- Hydrology column filter (soil and urban pervious road cols) ---
    # Matches Fortran filterMod.F90: excludes lake columns (they use LakeMod)
    for c in bounds.begc:bounds.endc
        if col.active[c] || include_inactive
            li = col.landunit[c]
            if col.hydrologically_active[c] && lun.itype[li] != ISTDLAK
                filt.hydrologyc[c] = true
            end
        end
    end

    # --- Prognostic crop / soil-no-pcrop patch filter ---
    if !varctl.use_fates
        for p in bounds.begp:bounds.endp
            if pch.active[p] || include_inactive
                if pch.itype[p] >= NPCROPMIN
                    filt.pcropp[p] = true
                else
                    l = pch.landunit[p]
                    if lun.itype[l] == ISTSOIL || lun.itype[l] == ISTCROP
                        filt.soilnopcropp[p] = true
                    end
                end
            end
        end
    end

    # --- Urban / non-urban landunit filter ---
    for l in bounds.begl:bounds.endl
        if lun.active[l] || include_inactive
            if lun.urbpoi[l]
                filt.urbanl[l] = true
            else
                filt.nourbanl[l] = true
            end
        end
    end

    # --- Urban / non-urban column filter ---
    for c in bounds.begc:bounds.endc
        if col.active[c] || include_inactive
            l = col.landunit[c]
            if lun.urbpoi[l]
                filt.urbanc[c] = true
            else
                filt.nourbanc[c] = true
            end
        end
    end

    # --- Urban / non-urban patch filter ---
    for p in bounds.begp:bounds.endp
        if pch.active[p] || include_inactive
            l = pch.landunit[p]
            if lun.urbpoi[l]
                filt.urbanp[p] = true
            else
                filt.nourbanp[p] = true
            end
        end
    end

    # --- Ice column filter ---
    for c in bounds.begc:bounds.endc
        if col.active[c] || include_inactive
            l = col.landunit[c]
            if lun.itype[l] == ISTICE
                filt.icec[c] = true
            end
        end
    end

    # --- SMB column filter ---
    if !isempty(melt_replaced_by_ice)
        for c in bounds.begc:bounds.endc
            if col.active[c] || include_inactive
                l = col.landunit[c]
                g = col.gridcell[c]
                if melt_replaced_by_ice[g] &&
                   (lun.itype[l] == ISTICE || lun.itype[l] == ISTSOIL)
                    filt.do_smb_c[c] = true
                end
            end
        end
    end

    # Note: snow filters (snowc, nosnowc, lakesnowc, lakenosnowc) are
    # reconstructed each time step in LakeHydrology and SnowHydrology.
    # Note: CNDV natvegp filter is reconstructed each time CNDV is run.

    return nothing
end

"""
    set_exposedvegp_filter!(filt::ClumpFilter, bounds::BoundsType,
                            frac_veg_nosno::Vector{Int})

Set the exposedvegp and noexposedvegp filter masks for one clump.

The exposedvegp filter includes nolakeurban patches where `frac_veg_nosno > 0`.
noexposedvegp includes nolakeurban patches where `frac_veg_nosno <= 0`.
Neither filter includes urban or lake points.

Only sets filters in the provided `filt`, NOT in
`clump_filter_inactive_and_active`.

Ported from `setExposedvegpFilter` in `filterMod.F90`.
"""
@kernel function _set_exposedvegp_kernel!(exposedvegp, noexposedvegp,
        @Const(nolakeurbanp), @Const(frac_veg_nosno), pmin::Int, pmax::Int)
    p = @index(Global)
    @inbounds if pmin <= p <= pmax
        exposedvegp[p] = false
        noexposedvegp[p] = false
        if nolakeurbanp[p]
            if frac_veg_nosno[p] > 0
                exposedvegp[p] = true
            else
                noexposedvegp[p] = true
            end
        end
    end
end

function set_exposedvegp_filter!(filt::ClumpFilter, bounds::BoundsType,
                                 frac_veg_nosno::AbstractVector{<:Integer})
    @assert bounds.level == BOUNDS_LEVEL_CLUMP
    @assert length(frac_veg_nosno) >= bounds.endp

    _launch!(_set_exposedvegp_kernel!, filt.exposedvegp, filt.noexposedvegp,
             filt.nolakeurbanp, frac_veg_nosno, bounds.begp, bounds.endp;
             ndrange = bounds.endp)

    return nothing
end
