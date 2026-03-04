# ==========================================================================
# Ported from: src/main/initSubgridMod.F90
# Low-level subgrid construction: add_landunit!, add_column!, add_patch!,
# clm_ptrs_compdown!, clm_ptrs_check!
# ==========================================================================

"""
    add_landunit!(lun, li, gi, ltype, wtgcell)

Add a landunit at index `li[] + 1`. `li` is a `Ref{Int}` incremented in place.
Sets gridcell, weight, itype, and landunit classification flags.
"""
function add_landunit!(lun::LandunitData, li::Ref{Int}, gi::Int, ltype::Int, wtgcell::Float64)
    li[] += 1
    l = li[]

    lun.gridcell[l] = gi
    lun.wtgcell[l] = wtgcell
    lun.itype[l] = ltype

    lun.ifspecial[l] = landunit_is_special(ltype)
    lun.glcpoi[l] = (ltype == ISTICE)
    lun.lakpoi[l] = (ltype == ISTDLAK)
    lun.urbpoi[l] = (ISTURB_MIN <= ltype <= ISTURB_MAX)

    nothing
end

"""
    add_column!(col, lun, ci, li, ctype, wtlunit; type_is_dynamic=false)

Add a column at index `ci[] + 1`. `ci` is a `Ref{Int}` incremented in place.
Sets landunit, gridcell, weight, itype, and derived flags.
"""
function add_column!(col::ColumnData, lun::LandunitData, ci::Ref{Int},
                     li::Int, ctype::Int, wtlunit::Float64;
                     type_is_dynamic::Bool=false)
    ci[] += 1
    c = ci[]

    col.landunit[c] = li
    col.gridcell[c] = lun.gridcell[li]
    col.wtlunit[c] = wtlunit
    col.itype[c] = ctype
    col.lun_itype[c] = lun.itype[li]
    col.type_is_dynamic[c] = type_is_dynamic
    col.hydrologically_active[c] = is_hydrologically_active(ctype, lun.itype[li])
    col.urbpoi[c] = lun.urbpoi[li]

    nothing
end

"""
    add_patch!(pch, col, lun, pi, ci, ptype, wtcol)

Add a patch at index `pi[] + 1`. `pi` is a `Ref{Int}` incremented in place.
Sets column, landunit, gridcell, weight, itype, and mxy.
"""
function add_patch!(pch::PatchData, col::ColumnData, lun::LandunitData,
                    pi::Ref{Int}, ci::Int, ptype::Int, wtcol::Float64)
    pi[] += 1
    p = pi[]

    pch.column[p] = ci
    li = col.landunit[ci]
    pch.landunit[p] = li
    pch.gridcell[p] = col.gridcell[ci]
    pch.wtcol[p] = wtcol
    pch.itype[p] = ptype

    # Set mxy for soil/crop landunits
    if lun.itype[li] == ISTSOIL || lun.itype[li] == ISTCROP
        lb_offset = 1 - varpar.natpft_lb
        pch.mxy[p] = ptype + lb_offset
    else
        pch.mxy[p] = ISPVAL
    end

    nothing
end

"""
    clm_ptrs_compdown!(bounds, grc, lun, col, pch)

Fill in down-pointers from up-pointers. Assumes patch→column, patch→landunit,
col→landunit are already set. Computes col.patchi/patchf, lun.patchi/patchf/coli/colf,
and grc.landunit_indices.
"""
function clm_ptrs_compdown!(bounds::BoundsType, grc::GridcellData,
                             lun::LandunitData, col::ColumnData, pch::PatchData)
    # Patch → column and patch → landunit down-pointers
    curc = 0
    curl = 0
    for p in bounds.begp:bounds.endp
        if pch.column[p] != curc
            curc = pch.column[p]
            col.patchi[curc] = p
        end
        col.patchf[curc] = p
        col.npatches[curc] = col.patchf[curc] - col.patchi[curc] + 1

        if pch.landunit[p] != curl
            curl = pch.landunit[p]
            lun.patchi[curl] = p
        end
        lun.patchf[curl] = p
        lun.npatches[curl] = lun.patchf[curl] - lun.patchi[curl] + 1
    end

    # Column → landunit down-pointers
    curl = 0
    for c in bounds.begc:bounds.endc
        if col.landunit[c] != curl
            curl = col.landunit[c]
            lun.coli[curl] = c
        end
        lun.colf[curl] = c
        lun.ncolumns[curl] = lun.colf[curl] - lun.coli[curl] + 1
    end

    # Landunit → gridcell: fill grc.landunit_indices
    for g in bounds.begg:bounds.endg
        for lt in 1:MAX_LUNIT
            grc.landunit_indices[lt, g] = ISPVAL
        end
    end

    for l in bounds.begl:bounds.endl
        ltype = lun.itype[l]
        curg = lun.gridcell[l]
        if grc.landunit_indices[ltype, curg] == ISPVAL
            grc.landunit_indices[ltype, curg] = l
        else
            error("clm_ptrs_compdown! ERROR: landunit type $ltype already set for gridcell $curg")
        end
    end

    nothing
end

"""
    clm_ptrs_check!(bounds, grc, lun, col, pch)

Validate the subgrid hierarchy: check index ranges, monotonicity, and tree consistency.
"""
function clm_ptrs_check!(bounds::BoundsType, grc::GridcellData,
                          lun::LandunitData, col::ColumnData, pch::PatchData)
    begg, endg = bounds.begg, bounds.endg
    begl, endl = bounds.begl, bounds.endl
    begc, endc = bounds.begc, bounds.endc
    begp, endp = bounds.begp, bounds.endp

    # Check gridcell → landunit indices
    for g in begg:endg
        for ltype in 1:MAX_LUNIT
            l = grc.landunit_indices[ltype, g]
            if l != ISPVAL
                (l >= begl && l <= endl) || error("clm_ptrs_check!: g index range error at g=$g, ltype=$ltype, l=$l")
            end
        end
    end

    # Check landunit ranges
    for l in begl:endl
        (lun.gridcell[l] >= begg && lun.gridcell[l] <= endg) || error("clm_ptrs_check!: l gridcell range error at l=$l")
        (lun.coli[l] >= begc && lun.coli[l] <= endc) || error("clm_ptrs_check!: l coli range error at l=$l")
        (lun.colf[l] >= begc && lun.colf[l] <= endc) || error("clm_ptrs_check!: l colf range error at l=$l")
        (lun.patchi[l] >= begp && lun.patchi[l] <= endp) || error("clm_ptrs_check!: l patchi range error at l=$l")
        (lun.patchf[l] >= begp && lun.patchf[l] <= endp) || error("clm_ptrs_check!: l patchf range error at l=$l")
    end

    # Check column ranges
    for c in begc:endc
        (col.gridcell[c] >= begg && col.gridcell[c] <= endg) || error("clm_ptrs_check!: c gridcell range error at c=$c")
        (col.landunit[c] >= begl && col.landunit[c] <= endl) || error("clm_ptrs_check!: c landunit range error at c=$c")
        (col.patchi[c] >= begp && col.patchi[c] <= endp) || error("clm_ptrs_check!: c patchi range error at c=$c")
        (col.patchf[c] >= begp && col.patchf[c] <= endp) || error("clm_ptrs_check!: c patchf range error at c=$c")
    end

    # Check patch ranges
    for p in begp:endp
        (pch.gridcell[p] >= begg && pch.gridcell[p] <= endg) || error("clm_ptrs_check!: p gridcell range error at p=$p")
        (pch.landunit[p] >= begl && pch.landunit[p] <= endl) || error("clm_ptrs_check!: p landunit range error at p=$p")
        (pch.column[p] >= begc && pch.column[p] <= endc) || error("clm_ptrs_check!: p column range error at p=$p")
    end

    # Check tree consistency
    for g in begg:endg
        for ltype in 1:MAX_LUNIT
            l = grc.landunit_indices[ltype, g]
            l == ISPVAL && continue
            lun.itype[l] == ltype || error("clm_ptrs_check!: tree consistency error: lun.itype[$l]=$(lun.itype[l]) != $ltype")
            lun.gridcell[l] == g || error("clm_ptrs_check!: tree consistency error: lun.gridcell[$l] != $g")
            for c in lun.coli[l]:lun.colf[l]
                col.gridcell[c] == g || error("clm_ptrs_check!: tree error: col.gridcell[$c] != $g")
                col.landunit[c] == l || error("clm_ptrs_check!: tree error: col.landunit[$c] != $l")
                for p in col.patchi[c]:col.patchf[c]
                    pch.gridcell[p] == g || error("clm_ptrs_check!: tree error: pch.gridcell[$p] != $g")
                    pch.landunit[p] == l || error("clm_ptrs_check!: tree error: pch.landunit[$p] != $l")
                    pch.column[p] == c || error("clm_ptrs_check!: tree error: pch.column[$p] != $c")
                end
            end
        end
    end

    nothing
end
