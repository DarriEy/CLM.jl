# ==========================================================================
# Ported from: src/main/column_varcon.F90
# Column type constants and utility functions
# ==========================================================================

# Urban column type indices (derived from ISTURB_MIN)
const ICOL_ROOF        = ISTURB_MIN * 10 + 1  # 71
const ICOL_SUNWALL     = ISTURB_MIN * 10 + 2  # 72
const ICOL_SHADEWALL   = ISTURB_MIN * 10 + 3  # 73
const ICOL_ROAD_IMPERV = ISTURB_MIN * 10 + 4  # 74
const ICOL_ROAD_PERV   = ISTURB_MIN * 10 + 5  # 75

"""
    is_hydrologically_active(col_itype::Int, lun_itype::Int) -> Bool

Determine whether a column is hydrologically active.
Soil, crop, lake, wetland, and pervious road columns are active.
Urban walls, roofs, and impervious roads are not.
"""
function is_hydrologically_active(col_itype::Int, lun_itype::Int)
    if lun_itype == ISTSOIL || lun_itype == ISTCROP
        return true
    elseif lun_itype == ISTICE
        return true
    elseif lun_itype == ISTDLAK
        return true
    elseif lun_itype == ISTWET
        return true
    elseif ISTURB_MIN <= lun_itype <= ISTURB_MAX
        return col_itype == ICOL_ROAD_PERV
    else
        return false
    end
end

"""
    ice_class_to_col_itype(ice_class::Int) -> Int

Convert a glacier ice class (1..maxpatch_glc) to a column itype.
Ice class 1 maps to col_itype 1, etc.
"""
function ice_class_to_col_itype(ice_class::Int)
    return ice_class
end

"""
    col_itype_to_ice_class(col_itype::Int) -> Int

Convert a column itype back to a glacier ice class.
"""
function col_itype_to_ice_class(col_itype::Int)
    return col_itype
end
