# ==========================================================================
# Ported from: src/main/subgridWeightsMod.F90
# Weight computation, active flag setting, and weight validation
# ==========================================================================

"""
    compute_higher_order_weights!(bounds, col, lun, pch)

Given patch%wtcol, col%wtlunit, and lun%wtgcell, compute:
- col%wtgcell = col%wtlunit * lun%wtgcell
- patch%wtlunit = patch%wtcol * col%wtlunit
- patch%wtgcell = patch%wtcol * col%wtgcell
"""
function compute_higher_order_weights!(bounds::BoundsType,
                                        col::ColumnData, lun::LandunitData,
                                        pch::PatchData)
    for c in bounds.begc:bounds.endc
        l = col.landunit[c]
        col.wtgcell[c] = col.wtlunit[c] * lun.wtgcell[l]
    end

    for p in bounds.begp:bounds.endp
        c = pch.column[p]
        pch.wtlunit[p] = pch.wtcol[p] * col.wtlunit[c]
        pch.wtgcell[p] = pch.wtcol[p] * col.wtgcell[c]
    end

    nothing
end

"""
    set_active!(bounds, lun, col, pch)

Set `active` flags based on weights and landunit types.
Simplified from Fortran: no glc_behavior dependency.
- Landunits: active if wtgcell > 0, or if soil/urban/lake type
- Columns: active if parent landunit active AND wtlunit > 0, or urban
- Patches: active if parent column active AND wtcol > 0
If `varctl.all_active`, everything is active.
"""
function set_active!(bounds::BoundsType, lun::LandunitData,
                     col::ColumnData, pch::PatchData)
    # Landunits
    for l in bounds.begl:bounds.endl
        if varctl.all_active
            lun.active[l] = true
        else
            lun.active[l] = false
            # Active if non-zero weight
            if lun.wtgcell[l] > 0.0
                lun.active[l] = true
            end
            # Soil, urban, lake landunits always active
            if lun.itype[l] == ISTSOIL
                lun.active[l] = true
            end
            if lun.itype[l] == ISTDLAK
                lun.active[l] = true
            end
            if ISTURB_MIN <= lun.itype[l] <= ISTURB_MAX
                lun.active[l] = true
            end
        end
    end

    # Columns
    for c in bounds.begc:bounds.endc
        if varctl.all_active
            col.active[c] = true
        else
            l = col.landunit[c]
            col.active[c] = false
            if lun.active[l] && col.wtlunit[c] > 0.0
                col.active[c] = true
            end
            # All urban columns on an active urban landunit are active
            if lun.active[l] && (ISTURB_MIN <= lun.itype[l] <= ISTURB_MAX)
                col.active[c] = true
            end
        end
    end

    # Patches
    for p in bounds.begp:bounds.endp
        if varctl.all_active
            pch.active[p] = true
        else
            c = pch.column[p]
            pch.active[p] = false
            if col.active[c] && pch.wtcol[p] > 0.0
                pch.active[p] = true
            end
        end
    end

    nothing
end

"""
    check_weights!(bounds, grc, lun, col, pch; active_only::Bool)

Check that subgrid weights satisfy requirements:
- If active_only=false: sum of ALL children weights == 1
- If active_only=true: sum of ACTIVE children weights == 1 for active parents,
  0 or 1 for inactive parents
"""
function check_weights!(bounds::BoundsType, grc::GridcellData,
                         lun::LandunitData, col::ColumnData, pch::PatchData;
                         active_only::Bool=false)
    tol = 1.0e-12

    # Check patch weights on columns
    sumwtcol = zeros(bounds.endc)
    for p in bounds.begp:bounds.endp
        if !active_only || pch.active[p]
            c = pch.column[p]
            sumwtcol[c] += pch.wtcol[p]
        end
    end
    for c in bounds.begc:bounds.endc
        _check_wt_ok(sumwtcol[c], active_only, col.active[c], tol,
                      "patch-on-col at c=$c")
    end

    # Check patch weights on landunits
    sumwtlunit = zeros(bounds.endl)
    for p in bounds.begp:bounds.endp
        if !active_only || pch.active[p]
            l = pch.landunit[p]
            sumwtlunit[l] += pch.wtlunit[p]
        end
    end
    for l in bounds.begl:bounds.endl
        _check_wt_ok(sumwtlunit[l], active_only, lun.active[l], tol,
                      "patch-on-lunit at l=$l")
    end

    # Check patch weights on gridcells
    sumwtgcell_p = zeros(bounds.endg)
    for p in bounds.begp:bounds.endp
        if !active_only || pch.active[p]
            g = pch.gridcell[p]
            sumwtgcell_p[g] += pch.wtgcell[p]
        end
    end
    for g in bounds.begg:bounds.endg
        _check_wt_ok(sumwtgcell_p[g], active_only, true, tol,
                      "patch-on-gcell at g=$g")
    end

    # Check column weights on landunits
    sumwtlunit_c = zeros(bounds.endl)
    for c in bounds.begc:bounds.endc
        if !active_only || col.active[c]
            l = col.landunit[c]
            sumwtlunit_c[l] += col.wtlunit[c]
        end
    end
    for l in bounds.begl:bounds.endl
        _check_wt_ok(sumwtlunit_c[l], active_only, lun.active[l], tol,
                      "col-on-lunit at l=$l")
    end

    # Check column weights on gridcells
    sumwtgcell_c = zeros(bounds.endg)
    for c in bounds.begc:bounds.endc
        if !active_only || col.active[c]
            g = col.gridcell[c]
            sumwtgcell_c[g] += col.wtgcell[c]
        end
    end
    for g in bounds.begg:bounds.endg
        _check_wt_ok(sumwtgcell_c[g], active_only, true, tol,
                      "col-on-gcell at g=$g")
    end

    # Check landunit weights on gridcells
    sumwtgcell_l = zeros(bounds.endg)
    for l in bounds.begl:bounds.endl
        if !active_only || lun.active[l]
            g = lun.gridcell[l]
            sumwtgcell_l[g] += lun.wtgcell[l]
        end
    end
    for g in bounds.begg:bounds.endg
        _check_wt_ok(sumwtgcell_l[g], active_only, true, tol,
                      "lunit-on-gcell at g=$g")
    end

    nothing
end

"""
    _check_wt_ok(sumwts, active_only, i_am_active, tol, context)

Helper: check if sumwts satisfies weight requirements.
"""
function _check_wt_ok(sumwts::Float64, active_only::Bool, i_am_active::Bool,
                       tol::Float64, context::String)
    wt_eq_1 = abs(sumwts - 1.0) <= tol
    if active_only
        if i_am_active
            wt_eq_1 || error("check_weights! ERROR: $context, sum=$sumwts (expected 1.0, active_only=true)")
        else
            (sumwts == 0.0 || wt_eq_1) || error("check_weights! ERROR: $context, sum=$sumwts (expected 0 or 1, inactive parent)")
        end
    else
        wt_eq_1 || error("check_weights! ERROR: $context, sum=$sumwts (expected 1.0, active_only=false)")
    end
    nothing
end

"""
    get_landunit_weight(grc, lun, g, ltype) -> Float64

Get weight of landunit type `ltype` on gridcell `g`. Returns 0 if not present.
"""
function get_landunit_weight(grc::GridcellData, lun::LandunitData, g::Int, ltype::Int)
    l = grc.landunit_indices[ltype, g]
    return l == ISPVAL ? 0.0 : lun.wtgcell[l]
end
