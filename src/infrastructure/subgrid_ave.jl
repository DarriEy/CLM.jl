# ==========================================================================
# Ported from: src/main/subgridAveMod.F90
# Utilities to perform subgrid averaging
# ==========================================================================
#
# Note about the urban scaling types used for c2l_scale_type (urbanf / urbans):
# 'urbanf' should be used for variables expressed as something-per-m^2
# ('extensive' state or flux variables), whereas 'urbans' should be used for
# variables not expressed as per-m^2 ('intensive' state variables, e.g.
# temperature). The urbanf scaling converts from per-m^2 of vertical wall area
# to per-m^2 of ground area.

# -----------------------------------------------------------------------
# Private helper: create_scale_l2g_lookup
# -----------------------------------------------------------------------

"""
    create_scale_l2g_lookup(l2g_scale_type::String) -> Vector{Float64}

Create a lookup array of length `MAX_LUNIT` giving the scale factor for each
landunit type depending on `l2g_scale_type`.

Supported scale types:
- "unity"          : all landunits included with scale 1
- "natveg"         : natural vegetation (soil) only
- "veg"            : vegetation (soil + crop)
- "ice"            : ice only
- "nonurb"         : all non-urban landunits
- "lake"           : lake only
- "veg_plus_lake"  : vegetation + lake

Ported from `create_scale_l2g_lookup` in `subgridAveMod.F90`.
"""
function create_scale_l2g_lookup(l2g_scale_type::String)
    # Initialize to spval for all landunits — any landunit that keeps the
    # default value will be excluded from grid cell averages.
    scale_lookup = fill(SPVAL, MAX_LUNIT)

    if l2g_scale_type == "unity"
        fill!(scale_lookup, 1.0)
    elseif l2g_scale_type == "natveg"
        scale_lookup[ISTSOIL] = 1.0
    elseif l2g_scale_type == "veg"
        scale_lookup[ISTSOIL] = 1.0
        scale_lookup[ISTCROP] = 1.0
    elseif l2g_scale_type == "ice"
        scale_lookup[ISTICE] = 1.0
    elseif l2g_scale_type == "nonurb"
        fill!(scale_lookup, 1.0)
        for i in ISTURB_MIN:ISTURB_MAX
            scale_lookup[i] = SPVAL
        end
    elseif l2g_scale_type == "lake"
        scale_lookup[ISTDLAK] = 1.0
    elseif l2g_scale_type == "veg_plus_lake"
        scale_lookup[ISTSOIL] = 1.0
        scale_lookup[ISTCROP] = 1.0
        scale_lookup[ISTDLAK] = 1.0
    else
        error("create_scale_l2g_lookup: scale type '$l2g_scale_type' not supported")
    end

    return scale_lookup
end

# -----------------------------------------------------------------------
# Private helper: build_scale_l2g!
# -----------------------------------------------------------------------

"""
    build_scale_l2g!(scale_l2g, bounds, l2g_scale_type, lun)

Fill `scale_l2g[bounds.begl:bounds.endl]` with appropriate values for the
given `l2g_scale_type`. This array can later be used to scale each landunit
in forming grid cell averages.

Ported from `build_scale_l2g` in `subgridAveMod.F90`.
"""
function build_scale_l2g!(scale_l2g::Vector{Float64}, bounds::BoundsType,
                          l2g_scale_type::String, lun::LandunitData)
    scale_lookup = create_scale_l2g_lookup(l2g_scale_type)
    for l in bounds.begl:bounds.endl
        scale_l2g[l] = scale_lookup[lun.itype[l]]
    end
    return nothing
end

# -----------------------------------------------------------------------
# Private helper: set_c2l_scale!
# -----------------------------------------------------------------------

"""
    set_c2l_scale!(scale_c2l, bounds, c2l_scale_type, col, lun)

Set `scale_c2l` for different `c2l_scale_type` values.

Supported scale types:
- "unity"  : scale factor of 1.0 for all columns
- "urbanf" : urban flux scaling (extensive variables, per-m^2)
- "urbans" : urban state scaling (intensive variables, e.g. temperature)

Ported from `set_c2l_scale` in `subgridAveMod.F90`.
"""
function set_c2l_scale!(scale_c2l::Vector{Float64}, bounds::BoundsType,
                        c2l_scale_type::String,
                        col::ColumnData, lun::LandunitData)
    if c2l_scale_type == "unity"
        for c in bounds.begc:bounds.endc
            scale_c2l[c] = 1.0
        end
    elseif c2l_scale_type == "urbanf"
        for c in bounds.begc:bounds.endc
            l = col.landunit[c]
            if lun.urbpoi[l]
                if col.itype[c] == ICOL_SUNWALL
                    scale_c2l[c] = 3.0 * lun.canyon_hwr[l]
                elseif col.itype[c] == ICOL_SHADEWALL
                    scale_c2l[c] = 3.0 * lun.canyon_hwr[l]
                elseif col.itype[c] == ICOL_ROAD_PERV || col.itype[c] == ICOL_ROAD_IMPERV
                    scale_c2l[c] = 3.0
                elseif col.itype[c] == ICOL_ROOF
                    scale_c2l[c] = 1.0
                end
            else
                scale_c2l[c] = 1.0
            end
        end
    elseif c2l_scale_type == "urbans"
        for c in bounds.begc:bounds.endc
            l = col.landunit[c]
            if lun.urbpoi[l]
                if col.itype[c] == ICOL_SUNWALL
                    scale_c2l[c] = (3.0 * lun.canyon_hwr[l]) / (2.0 * lun.canyon_hwr[l] + 1.0)
                elseif col.itype[c] == ICOL_SHADEWALL
                    scale_c2l[c] = (3.0 * lun.canyon_hwr[l]) / (2.0 * lun.canyon_hwr[l] + 1.0)
                elseif col.itype[c] == ICOL_ROAD_PERV || col.itype[c] == ICOL_ROAD_IMPERV
                    scale_c2l[c] = 3.0 / (2.0 * lun.canyon_hwr[l] + 1.0)
                elseif col.itype[c] == ICOL_ROOF
                    scale_c2l[c] = 1.0
                end
            else
                scale_c2l[c] = 1.0
            end
        end
    else
        error("set_c2l_scale!: scale type '$c2l_scale_type' not supported")
    end
    return nothing
end

# =======================================================================
# p2c: Patch to Column averaging
# =======================================================================

"""
    p2c_1d!(carr, parr, bounds, p2c_scale_type, pch)

Perform subgrid-average from patches to columns (1D).
Averaging is only done for points that are not equal to `SPVAL`.

Ported from `p2c_1d` in `subgridAveMod.F90`.
"""
function p2c_1d!(carr::Vector{Float64}, parr::Vector{Float64},
                 bounds::BoundsType, p2c_scale_type::String,
                 pch::PatchData)
    # Set scale_p2c
    np = bounds.endp - bounds.begp + 1
    scale_p2c = Vector{Float64}(undef, bounds.endp)
    if p2c_scale_type == "unity"
        for p in bounds.begp:bounds.endp
            scale_p2c[p] = 1.0
        end
    else
        error("p2c_1d! error: scale type '$p2c_scale_type' not supported")
    end

    nc = bounds.endc - bounds.begc + 1
    sumwt = zeros(bounds.endc)

    for c in bounds.begc:bounds.endc
        carr[c] = SPVAL
    end

    for p in bounds.begp:bounds.endp
        if pch.active[p] && pch.wtcol[p] != 0.0
            if parr[p] != SPVAL
                c = pch.column[p]
                if sumwt[c] == 0.0
                    carr[c] = 0.0
                end
                carr[c] = carr[c] + parr[p] * scale_p2c[p] * pch.wtcol[p]
                sumwt[c] = sumwt[c] + pch.wtcol[p]
            end
        end
    end

    found = false
    index = 0
    for c in bounds.begc:bounds.endc
        if sumwt[c] > 1.0 + 1.0e-6
            found = true
            index = c
        elseif sumwt[c] != 0.0
            carr[c] = carr[c] / sumwt[c]
        end
    end
    if found
        error("p2c_1d! error: sumwt is greater than 1.0 at c=$index")
    end

    return nothing
end

"""
    p2c_2d!(carr, parr, bounds, num2d, p2c_scale_type, pch)

Perform subgrid-average from patches to columns (2D).
Averaging is only done for points that are not equal to `SPVAL`.

Ported from `p2c_2d` in `subgridAveMod.F90`.
"""
function p2c_2d!(carr::Matrix{Float64}, parr::Matrix{Float64},
                 bounds::BoundsType, num2d::Int, p2c_scale_type::String,
                 pch::PatchData)
    # Set scale_p2c
    scale_p2c = Vector{Float64}(undef, bounds.endp)
    if p2c_scale_type == "unity"
        for p in bounds.begp:bounds.endp
            scale_p2c[p] = 1.0
        end
    else
        error("p2c_2d! error: scale type '$p2c_scale_type' not supported")
    end

    sumwt = zeros(bounds.endc)

    for c in bounds.begc:bounds.endc
        for j in 1:num2d
            carr[c, j] = SPVAL
        end
    end

    for j in 1:num2d
        fill!(sumwt, 0.0)
        for p in bounds.begp:bounds.endp
            if pch.active[p] && pch.wtcol[p] != 0.0
                if parr[p, j] != SPVAL
                    c = pch.column[p]
                    if sumwt[c] == 0.0
                        carr[c, j] = 0.0
                    end
                    carr[c, j] = carr[c, j] + parr[p, j] * scale_p2c[p] * pch.wtcol[p]
                    sumwt[c] = sumwt[c] + pch.wtcol[p]
                end
            end
        end
        found = false
        index = 0
        for c in bounds.begc:bounds.endc
            if sumwt[c] > 1.0 + 1.0e-6
                found = true
                index = c
            elseif sumwt[c] != 0.0
                carr[c, j] = carr[c, j] / sumwt[c]
            end
        end
        if found
            error("p2c_2d! error: sumwt is greater than 1.0 at c=$index, lev=$j")
        end
    end

    return nothing
end

"""
    p2c_1d_filter!(colarr, patcharr, mask_c, col, pch)

Perform patch to column averaging for single level patch arrays using a
column mask. For each active column in the mask, averages over its patches.

Ported from `p2c_1d_filter` in `subgridAveMod.F90`.
"""
function p2c_1d_filter!(colarr::Vector{Float64}, patcharr::Vector{Float64},
                        mask_c::BitVector, col::ColumnData, pch::PatchData)
    for c in eachindex(mask_c)
        mask_c[c] || continue
        colarr[c] = 0.0
        for p in col.patchi[c]:col.patchf[c]
            if pch.active[p]
                colarr[c] = colarr[c] + patcharr[p] * pch.wtcol[p]
            end
        end
    end
    return nothing
end

"""
    p2c_2d_filter!(colarr, patcharr, lev, mask_c, col, pch)

Perform patch to column averaging for multi-level patch arrays using a
column mask. For each active column in the mask, averages over its patches
for each level.

Ported from `p2c_2d_filter` in `subgridAveMod.F90`.
"""
function p2c_2d_filter!(colarr::Matrix{Float64}, patcharr::Matrix{Float64},
                        lev::Int, mask_c::BitVector,
                        col::ColumnData, pch::PatchData)
    for j in 1:lev
        for c in eachindex(mask_c)
            mask_c[c] || continue
            colarr[c, j] = 0.0
            for p in col.patchi[c]:col.patchf[c]
                if pch.active[p]
                    colarr[c, j] = colarr[c, j] + patcharr[p, j] * pch.wtcol[p]
                end
            end
        end
    end
    return nothing
end

# =======================================================================
# p2l: Patch to Landunit averaging
# =======================================================================

"""
    p2l_1d!(larr, parr, bounds, p2c_scale_type, c2l_scale_type, pch, col, lun)

Perform subgrid-average from patches to landunits (1D).
Averaging is only done for points that are not equal to `SPVAL`.

Ported from `p2l_1d` in `subgridAveMod.F90`.
"""
function p2l_1d!(larr::Vector{Float64}, parr::Vector{Float64},
                 bounds::BoundsType,
                 p2c_scale_type::String, c2l_scale_type::String,
                 pch::PatchData, col::ColumnData, lun::LandunitData)
    # Set scale_c2l
    scale_c2l = Vector{Float64}(undef, bounds.endc)
    set_c2l_scale!(scale_c2l, bounds, c2l_scale_type, col, lun)

    # Set scale_p2c
    scale_p2c = Vector{Float64}(undef, bounds.endp)
    if p2c_scale_type == "unity"
        for p in bounds.begp:bounds.endp
            scale_p2c[p] = 1.0
        end
    else
        error("p2l_1d! error: scale type '$p2c_scale_type' not supported")
    end

    sumwt = zeros(bounds.endl)
    for l in bounds.begl:bounds.endl
        larr[l] = SPVAL
    end

    for p in bounds.begp:bounds.endp
        if pch.active[p] && pch.wtlunit[p] != 0.0
            c = pch.column[p]
            if parr[p] != SPVAL && scale_c2l[c] != SPVAL
                l = pch.landunit[p]
                if sumwt[l] == 0.0
                    larr[l] = 0.0
                end
                larr[l] = larr[l] + parr[p] * scale_p2c[p] * scale_c2l[c] * pch.wtlunit[p]
                sumwt[l] = sumwt[l] + pch.wtlunit[p]
            end
        end
    end

    found = false
    index = 0
    for l in bounds.begl:bounds.endl
        if sumwt[l] > 1.0 + 1.0e-6
            found = true
            index = l
        elseif sumwt[l] != 0.0
            larr[l] = larr[l] / sumwt[l]
        end
    end
    if found
        error("p2l_1d! error: sumwt is greater than 1.0 at l=$index")
    end

    return nothing
end

"""
    p2l_2d!(larr, parr, bounds, num2d, p2c_scale_type, c2l_scale_type, pch, col, lun)

Perform subgrid-average from patches to landunits (2D).
Averaging is only done for points that are not equal to `SPVAL`.

Ported from `p2l_2d` in `subgridAveMod.F90`.
"""
function p2l_2d!(larr::Matrix{Float64}, parr::Matrix{Float64},
                 bounds::BoundsType, num2d::Int,
                 p2c_scale_type::String, c2l_scale_type::String,
                 pch::PatchData, col::ColumnData, lun::LandunitData)
    # Set scale_c2l
    scale_c2l = Vector{Float64}(undef, bounds.endc)
    set_c2l_scale!(scale_c2l, bounds, c2l_scale_type, col, lun)

    # Set scale_p2c
    scale_p2c = Vector{Float64}(undef, bounds.endp)
    if p2c_scale_type == "unity"
        for p in bounds.begp:bounds.endp
            scale_p2c[p] = 1.0
        end
    else
        error("p2l_2d! error: scale type '$p2c_scale_type' not supported")
    end

    for l in bounds.begl:bounds.endl
        for j in 1:num2d
            larr[l, j] = SPVAL
        end
    end

    sumwt = zeros(bounds.endl)
    for j in 1:num2d
        fill!(sumwt, 0.0)
        for p in bounds.begp:bounds.endp
            if pch.active[p] && pch.wtlunit[p] != 0.0
                c = pch.column[p]
                if parr[p, j] != SPVAL && scale_c2l[c] != SPVAL
                    l = pch.landunit[p]
                    if sumwt[l] == 0.0
                        larr[l, j] = 0.0
                    end
                    larr[l, j] = larr[l, j] + parr[p, j] * scale_p2c[p] * scale_c2l[c] * pch.wtlunit[p]
                    sumwt[l] = sumwt[l] + pch.wtlunit[p]
                end
            end
        end
        found = false
        index = 0
        for l in bounds.begl:bounds.endl
            if sumwt[l] > 1.0 + 1.0e-6
                found = true
                index = l
            elseif sumwt[l] != 0.0
                larr[l, j] = larr[l, j] / sumwt[l]
            end
        end
        if found
            error("p2l_2d! error: sumwt is greater than 1.0 at l=$index, j=$j")
        end
    end

    return nothing
end

# =======================================================================
# p2g: Patch to Gridcell averaging
# =======================================================================

"""
    p2g_1d!(garr, parr, bounds, p2c_scale_type, c2l_scale_type, l2g_scale_type, pch, col, lun)

Perform subgrid-average from patches to gridcells (1D).
Averaging is only done for points that are not equal to `SPVAL`.

Ported from `p2g_1d` in `subgridAveMod.F90`.
"""
function p2g_1d!(garr::Vector{Float64}, parr::Vector{Float64},
                 bounds::BoundsType,
                 p2c_scale_type::String, c2l_scale_type::String,
                 l2g_scale_type::String,
                 pch::PatchData, col::ColumnData, lun::LandunitData)
    # Build scale_l2g
    scale_l2g = Vector{Float64}(undef, bounds.endl)
    build_scale_l2g!(scale_l2g, bounds, l2g_scale_type, lun)

    # Set scale_c2l
    scale_c2l = Vector{Float64}(undef, bounds.endc)
    set_c2l_scale!(scale_c2l, bounds, c2l_scale_type, col, lun)

    # Set scale_p2c
    scale_p2c = Vector{Float64}(undef, bounds.endp)
    if p2c_scale_type == "unity"
        for p in bounds.begp:bounds.endp
            scale_p2c[p] = 1.0
        end
    else
        error("p2g_1d! error: scale type '$p2c_scale_type' not supported")
    end

    sumwt = zeros(bounds.endg)
    for g in bounds.begg:bounds.endg
        garr[g] = SPVAL
    end

    for p in bounds.begp:bounds.endp
        if pch.active[p] && pch.wtgcell[p] != 0.0
            c = pch.column[p]
            l = pch.landunit[p]
            if parr[p] != SPVAL && scale_c2l[c] != SPVAL && scale_l2g[l] != SPVAL
                g = pch.gridcell[p]
                if sumwt[g] == 0.0
                    garr[g] = 0.0
                end
                garr[g] = garr[g] + parr[p] * scale_p2c[p] * scale_c2l[c] * scale_l2g[l] * pch.wtgcell[p]
                sumwt[g] = sumwt[g] + pch.wtgcell[p]
            end
        end
    end

    found = false
    index = 0
    for g in bounds.begg:bounds.endg
        if sumwt[g] > 1.0 + 1.0e-6
            found = true
            index = g
        elseif sumwt[g] != 0.0
            garr[g] = garr[g] / sumwt[g]
        end
    end
    if found
        error("p2g_1d! error: sumwt is greater than 1.0 at g=$index")
    end

    return nothing
end

"""
    p2g_2d!(garr, parr, bounds, num2d, p2c_scale_type, c2l_scale_type, l2g_scale_type, pch, col, lun)

Perform subgrid-average from patches to gridcells (2D).
Averaging is only done for points that are not equal to `SPVAL`.

Ported from `p2g_2d` in `subgridAveMod.F90`.
"""
function p2g_2d!(garr::Matrix{Float64}, parr::Matrix{Float64},
                 bounds::BoundsType, num2d::Int,
                 p2c_scale_type::String, c2l_scale_type::String,
                 l2g_scale_type::String,
                 pch::PatchData, col::ColumnData, lun::LandunitData)
    # Build scale_l2g
    scale_l2g = Vector{Float64}(undef, bounds.endl)
    build_scale_l2g!(scale_l2g, bounds, l2g_scale_type, lun)

    # Set scale_c2l
    scale_c2l = Vector{Float64}(undef, bounds.endc)
    set_c2l_scale!(scale_c2l, bounds, c2l_scale_type, col, lun)

    # Set scale_p2c
    scale_p2c = Vector{Float64}(undef, bounds.endp)
    if p2c_scale_type == "unity"
        for p in bounds.begp:bounds.endp
            scale_p2c[p] = 1.0
        end
    else
        error("p2g_2d! error: scale type '$p2c_scale_type' not supported")
    end

    for g in bounds.begg:bounds.endg
        for j in 1:num2d
            garr[g, j] = SPVAL
        end
    end

    sumwt = zeros(bounds.endg)
    for j in 1:num2d
        fill!(sumwt, 0.0)
        for p in bounds.begp:bounds.endp
            if pch.active[p] && pch.wtgcell[p] != 0.0
                c = pch.column[p]
                l = pch.landunit[p]
                if parr[p, j] != SPVAL && scale_c2l[c] != SPVAL && scale_l2g[l] != SPVAL
                    g = pch.gridcell[p]
                    if sumwt[g] == 0.0
                        garr[g, j] = 0.0
                    end
                    garr[g, j] = garr[g, j] + parr[p, j] * scale_p2c[p] * scale_c2l[c] * scale_l2g[l] * pch.wtgcell[p]
                    sumwt[g] = sumwt[g] + pch.wtgcell[p]
                end
            end
        end
        found = false
        index = 0
        for g in bounds.begg:bounds.endg
            if sumwt[g] > 1.0 + 1.0e-6
                found = true
                index = g
            elseif sumwt[g] != 0.0
                garr[g, j] = garr[g, j] / sumwt[g]
            end
        end
        if found
            error("p2g_2d! error: sumwt is greater than 1.0 at g=$index, j=$j")
        end
    end

    return nothing
end

# =======================================================================
# c2l: Column to Landunit averaging
# =======================================================================

"""
    c2l_1d!(larr, carr, bounds, c2l_scale_type, col, lun; include_inactive=false)

Perform subgrid-average from columns to landunits (1D).
Averaging is only done for points that are not equal to `SPVAL`.

If `include_inactive` is `true`, include inactive as well as active columns
in the averages. The purpose is to produce valid landunit-level output for
inactive landunits. This should only be set if `carr` has no NaN values,
even for inactive columns.

Ported from `c2l_1d` in `subgridAveMod.F90`.
"""
function c2l_1d!(larr::Vector{Float64}, carr::Vector{Float64},
                 bounds::BoundsType, c2l_scale_type::String,
                 col::ColumnData, lun::LandunitData;
                 include_inactive::Bool=false)
    # Set scale_c2l
    scale_c2l = Vector{Float64}(undef, bounds.endc)
    set_c2l_scale!(scale_c2l, bounds, c2l_scale_type, col, lun)

    sumwt = zeros(bounds.endl)
    for l in bounds.begl:bounds.endl
        larr[l] = SPVAL
    end

    for c in bounds.begc:bounds.endc
        if (col.active[c] || include_inactive) && col.wtlunit[c] != 0.0
            if carr[c] != SPVAL && scale_c2l[c] != SPVAL
                l = col.landunit[c]
                if sumwt[l] == 0.0
                    larr[l] = 0.0
                end
                larr[l] = larr[l] + carr[c] * scale_c2l[c] * col.wtlunit[c]
                sumwt[l] = sumwt[l] + col.wtlunit[c]
            end
        end
    end

    found = false
    index = 0
    for l in bounds.begl:bounds.endl
        if sumwt[l] > 1.0 + 1.0e-6
            found = true
            index = l
        elseif sumwt[l] != 0.0
            larr[l] = larr[l] / sumwt[l]
        end
    end
    if found
        error("c2l_1d! error: sumwt is greater than 1.0 at l=$index")
    end

    return nothing
end

"""
    c2l_2d!(larr, carr, bounds, num2d, c2l_scale_type, col, lun)

Perform subgrid-average from columns to landunits (2D).
Averaging is only done for points that are not equal to `SPVAL`.

Ported from `c2l_2d` in `subgridAveMod.F90`.
"""
function c2l_2d!(larr::Matrix{Float64}, carr::Matrix{Float64},
                 bounds::BoundsType, num2d::Int, c2l_scale_type::String,
                 col::ColumnData, lun::LandunitData)
    # Set scale_c2l
    scale_c2l = Vector{Float64}(undef, bounds.endc)
    set_c2l_scale!(scale_c2l, bounds, c2l_scale_type, col, lun)

    for l in bounds.begl:bounds.endl
        for j in 1:num2d
            larr[l, j] = SPVAL
        end
    end

    sumwt = zeros(bounds.endl)
    for j in 1:num2d
        fill!(sumwt, 0.0)
        for c in bounds.begc:bounds.endc
            if col.active[c] && col.wtlunit[c] != 0.0
                if carr[c, j] != SPVAL && scale_c2l[c] != SPVAL
                    l = col.landunit[c]
                    if sumwt[l] == 0.0
                        larr[l, j] = 0.0
                    end
                    larr[l, j] = larr[l, j] + carr[c, j] * scale_c2l[c] * col.wtlunit[c]
                    sumwt[l] = sumwt[l] + col.wtlunit[c]
                end
            end
        end
        found = false
        index = 0
        for l in bounds.begl:bounds.endl
            if sumwt[l] > 1.0 + 1.0e-6
                found = true
                index = l
            elseif sumwt[l] != 0.0
                larr[l, j] = larr[l, j] / sumwt[l]
            end
        end
        if found
            error("c2l_2d! error: sumwt is greater than 1.0 at l=$index, lev=$j")
        end
    end

    return nothing
end

# =======================================================================
# c2g: Column to Gridcell averaging
# =======================================================================

"""
    c2g_1d!(garr, carr, bounds, c2l_scale_type, l2g_scale_type, col, lun)

Perform subgrid-average from columns to gridcells (1D).
Averaging is only done for points that are not equal to `SPVAL`.

Ported from `c2g_1d` in `subgridAveMod.F90`.
"""
function c2g_1d!(garr::Vector{Float64}, carr::Vector{Float64},
                 bounds::BoundsType,
                 c2l_scale_type::String, l2g_scale_type::String,
                 col::ColumnData, lun::LandunitData)
    # Build scale_l2g
    scale_l2g = Vector{Float64}(undef, bounds.endl)
    build_scale_l2g!(scale_l2g, bounds, l2g_scale_type, lun)

    # Set scale_c2l
    scale_c2l = Vector{Float64}(undef, bounds.endc)
    set_c2l_scale!(scale_c2l, bounds, c2l_scale_type, col, lun)

    sumwt = zeros(bounds.endg)
    for g in bounds.begg:bounds.endg
        garr[g] = SPVAL
    end

    for c in bounds.begc:bounds.endc
        if col.active[c] && col.wtgcell[c] != 0.0
            l = col.landunit[c]
            if carr[c] != SPVAL && scale_c2l[c] != SPVAL && scale_l2g[l] != SPVAL
                g = col.gridcell[c]
                if sumwt[g] == 0.0
                    garr[g] = 0.0
                end
                garr[g] = garr[g] + carr[c] * scale_c2l[c] * scale_l2g[l] * col.wtgcell[c]
                sumwt[g] = sumwt[g] + col.wtgcell[c]
            end
        end
    end

    found = false
    index = 0
    for g in bounds.begg:bounds.endg
        if sumwt[g] > 1.0 + 1.0e-6
            found = true
            index = g
        elseif sumwt[g] != 0.0
            garr[g] = garr[g] / sumwt[g]
        end
    end
    if found
        error("c2g_1d! error: sumwt is greater than 1.0 at g=$index")
    end

    return nothing
end

"""
    c2g_2d!(garr, carr, bounds, num2d, c2l_scale_type, l2g_scale_type, col, lun)

Perform subgrid-average from columns to gridcells (2D).
Averaging is only done for points that are not equal to `SPVAL`.

Ported from `c2g_2d` in `subgridAveMod.F90`.
"""
function c2g_2d!(garr::Matrix{Float64}, carr::Matrix{Float64},
                 bounds::BoundsType, num2d::Int,
                 c2l_scale_type::String, l2g_scale_type::String,
                 col::ColumnData, lun::LandunitData)
    # Build scale_l2g
    scale_l2g = Vector{Float64}(undef, bounds.endl)
    build_scale_l2g!(scale_l2g, bounds, l2g_scale_type, lun)

    # Set scale_c2l
    scale_c2l = Vector{Float64}(undef, bounds.endc)
    set_c2l_scale!(scale_c2l, bounds, c2l_scale_type, col, lun)

    for g in bounds.begg:bounds.endg
        for j in 1:num2d
            garr[g, j] = SPVAL
        end
    end

    sumwt = zeros(bounds.endg)
    for j in 1:num2d
        fill!(sumwt, 0.0)
        for c in bounds.begc:bounds.endc
            if col.active[c] && col.wtgcell[c] != 0.0
                l = col.landunit[c]
                if carr[c, j] != SPVAL && scale_c2l[c] != SPVAL && scale_l2g[l] != SPVAL
                    g = col.gridcell[c]
                    if sumwt[g] == 0.0
                        garr[g, j] = 0.0
                    end
                    garr[g, j] = garr[g, j] + carr[c, j] * scale_c2l[c] * scale_l2g[l] * col.wtgcell[c]
                    sumwt[g] = sumwt[g] + col.wtgcell[c]
                end
            end
        end
        found = false
        index = 0
        for g in bounds.begg:bounds.endg
            if sumwt[g] > 1.0 + 1.0e-6
                found = true
                index = g
            elseif sumwt[g] != 0.0
                garr[g, j] = garr[g, j] / sumwt[g]
            end
        end
        if found
            error("c2g_2d! error: sumwt is greater than 1.0 at g=$index, lev=$j")
        end
    end

    return nothing
end

# =======================================================================
# l2g: Landunit to Gridcell averaging
# =======================================================================

"""
    l2g_1d!(garr, larr, bounds, l2g_scale_type, lun)

Perform subgrid-average from landunits to gridcells (1D).
Averaging is only done for points that are not equal to `SPVAL`.

Ported from `l2g_1d` in `subgridAveMod.F90`.
"""
function l2g_1d!(garr::Vector{Float64}, larr::Vector{Float64},
                 bounds::BoundsType, l2g_scale_type::String,
                 lun::LandunitData)
    # Build scale_l2g
    scale_l2g = Vector{Float64}(undef, bounds.endl)
    build_scale_l2g!(scale_l2g, bounds, l2g_scale_type, lun)

    sumwt = zeros(bounds.endg)
    for g in bounds.begg:bounds.endg
        garr[g] = SPVAL
    end

    for l in bounds.begl:bounds.endl
        if lun.active[l] && lun.wtgcell[l] != 0.0
            if larr[l] != SPVAL && scale_l2g[l] != SPVAL
                g = lun.gridcell[l]
                if sumwt[g] == 0.0
                    garr[g] = 0.0
                end
                garr[g] = garr[g] + larr[l] * scale_l2g[l] * lun.wtgcell[l]
                sumwt[g] = sumwt[g] + lun.wtgcell[l]
            end
        end
    end

    found = false
    index = 0
    for g in bounds.begg:bounds.endg
        if sumwt[g] > 1.0 + 1.0e-6
            found = true
            index = g
        elseif sumwt[g] != 0.0
            garr[g] = garr[g] / sumwt[g]
        end
    end
    if found
        error("l2g_1d! error: sumwt is greater than 1.0 at g=$index")
    end

    return nothing
end

"""
    l2g_2d!(garr, larr, bounds, num2d, l2g_scale_type, lun)

Perform subgrid-average from landunits to gridcells (2D).
Averaging is only done for points that are not equal to `SPVAL`.

Ported from `l2g_2d` in `subgridAveMod.F90`.
"""
function l2g_2d!(garr::Matrix{Float64}, larr::Matrix{Float64},
                 bounds::BoundsType, num2d::Int, l2g_scale_type::String,
                 lun::LandunitData)
    # Build scale_l2g
    scale_l2g = Vector{Float64}(undef, bounds.endl)
    build_scale_l2g!(scale_l2g, bounds, l2g_scale_type, lun)

    for g in bounds.begg:bounds.endg
        for j in 1:num2d
            garr[g, j] = SPVAL
        end
    end

    sumwt = zeros(bounds.endg)
    for j in 1:num2d
        fill!(sumwt, 0.0)
        for l in bounds.begl:bounds.endl
            if lun.active[l] && lun.wtgcell[l] != 0.0
                if larr[l, j] != SPVAL && scale_l2g[l] != SPVAL
                    g = lun.gridcell[l]
                    if sumwt[g] == 0.0
                        garr[g, j] = 0.0
                    end
                    garr[g, j] = garr[g, j] + larr[l, j] * scale_l2g[l] * lun.wtgcell[l]
                    sumwt[g] = sumwt[g] + lun.wtgcell[l]
                end
            end
        end
        found = false
        index = 0
        for g in bounds.begg:bounds.endg
            if sumwt[g] > 1.0 + 1.0e-6
                found = true
                index = g
            elseif sumwt[g] != 0.0
                garr[g, j] = garr[g, j] / sumwt[g]
            end
        end
        if found
            error("l2g_2d! error: sumwt is greater than 1.0 at g=$index, lev=$j")
        end
    end

    return nothing
end
