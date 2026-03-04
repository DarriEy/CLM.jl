# ==========================================================================
# Ported from: src/main/initGridCellsMod.F90
# High-level grid construction: builds g/l/c/p hierarchy from surface data
# ==========================================================================

"""
    initGridCells!(bounds, surf, grc, lun, col, pch; lat=0.0, lon=0.0, area=1.0)

Master builder for the subgrid hierarchy. Creates landunits, columns, and patches
for each gridcell based on surface data weights.

Landunit creation order: ISTSOIL → ISTCROP → ISTURB_TBD → ISTURB_HD → ISTURB_MD
→ ISTDLAK → ISTWET → ISTICE
"""
function initGridCells!(bounds::BoundsType, surf::SurfaceInputData,
                         grc::GridcellData, lun::LandunitData,
                         col::ColumnData, pch::PatchData;
                         lat::Float64=0.0, lon::Float64=0.0, area::Float64=1.0)
    li = Ref(bounds.begl - 1)
    ci = Ref(bounds.begc - 1)
    pi = Ref(bounds.begp - 1)

    # Natural vegetation
    for g in bounds.begg:bounds.endg
        set_landunit_veg_compete!(surf, g, li, ci, pi, lun, col, pch)
    end

    # Crop
    for g in bounds.begg:bounds.endg
        set_landunit_crop_noncompete!(surf, g, li, ci, pi, lun, col, pch)
    end

    # Urban TBD, HD, MD
    for ltype in [ISTURB_TBD, ISTURB_HD, ISTURB_MD]
        for g in bounds.begg:bounds.endg
            set_landunit_urban!(surf, ltype, g, li, ci, pi, lun, col, pch)
        end
    end

    # Lake
    for g in bounds.begg:bounds.endg
        set_landunit_wet_lake!(surf, ISTDLAK, g, li, ci, pi, lun, col, pch)
    end

    # Wetland
    for g in bounds.begg:bounds.endg
        set_landunit_wet_lake!(surf, ISTWET, g, li, ci, pi, lun, col, pch)
    end

    # Glacier
    for g in bounds.begg:bounds.endg
        set_landunit_ice!(surf, g, li, ci, pi, lun, col, pch)
    end

    # Verify counts
    li[] == bounds.endl || error("initGridCells!: landunit count mismatch: $(li[]) != $(bounds.endl)")
    ci[] == bounds.endc || error("initGridCells!: column count mismatch: $(ci[]) != $(bounds.endc)")
    pi[] == bounds.endp || error("initGridCells!: patch count mismatch: $(pi[]) != $(bounds.endp)")

    # Set gridcell properties
    for g in bounds.begg:bounds.endg
        grc.area[g] = area
        grc.latdeg[g] = lat
        grc.londeg[g] = lon
        grc.lat[g] = lat * RPI / 180.0
        grc.lon[g] = lon * RPI / 180.0
    end

    # Fill down-pointers
    clm_ptrs_compdown!(bounds, grc, lun, col, pch)

    # Validate
    clm_ptrs_check!(bounds, grc, lun, col, pch)

    # Compute higher-order weights
    compute_higher_order_weights!(bounds, col, lun, pch)

    nothing
end

# --- Private per-landunit builders ---

function set_landunit_veg_compete!(surf::SurfaceInputData, gi::Int,
                                    li::Ref{Int}, ci::Ref{Int}, pi::Ref{Int},
                                    lun::LandunitData, col::ColumnData, pch::PatchData)
    wt = surf.wt_lunit[gi, ISTSOIL]
    natpft_size = size(surf.wt_nat_patch, 2)

    # Always create soil landunit
    add_landunit!(lun, li, gi, ISTSOIL, wt)
    add_column!(col, lun, ci, li[], 1, 1.0)

    # Add patches for each PFT with non-zero weight
    n_added = 0
    for m in 1:natpft_size
        p_wt = surf.wt_nat_patch[gi, m]
        if p_wt > 0.0
            # PFT type: 0-based (m=1 → ptype=0, m=2 → ptype=1, ...)
            ptype = varpar.natpft_lb + (m - 1)
            add_patch!(pch, col, lun, pi, ci[], ptype, p_wt)
            n_added += 1
        end
    end
    # Must have at least one patch (bare ground)
    if n_added == 0
        add_patch!(pch, col, lun, pi, ci[], noveg, 1.0)
    end

    nothing
end

function set_landunit_crop_noncompete!(surf::SurfaceInputData, gi::Int,
                                        li::Ref{Int}, ci::Ref{Int}, pi::Ref{Int},
                                        lun::LandunitData, col::ColumnData, pch::PatchData)
    varctl.create_crop_landunit || return nothing
    surf.wt_lunit[gi, ISTCROP] > 0.0 || return nothing

    cft_size = size(surf.wt_cft, 2)
    cft_size > 0 || return nothing

    wt = surf.wt_lunit[gi, ISTCROP]
    add_landunit!(lun, li, gi, ISTCROP, wt)

    n_added = 0
    for m in 1:cft_size
        cft_wt = surf.wt_cft[gi, m]
        if cft_wt > 0.0
            cft_type = varpar.cft_lb + (m - 1)
            ctype = ISTCROP * 100 + cft_type
            add_column!(col, lun, ci, li[], ctype, cft_wt)
            add_patch!(pch, col, lun, pi, ci[], cft_type, 1.0)
            n_added += 1
        end
    end

    # Fallback: at least one crop column
    if n_added == 0
        cft_type = varpar.cft_lb
        ctype = ISTCROP * 100 + cft_type
        add_column!(col, lun, ci, li[], ctype, 1.0)
        add_patch!(pch, col, lun, pi, ci[], cft_type, 1.0)
    end

    nothing
end

function set_landunit_urban!(surf::SurfaceInputData, ltype::Int, gi::Int,
                              li::Ref{Int}, ci::Ref{Int}, pi::Ref{Int},
                              lun::LandunitData, col::ColumnData, pch::PatchData)
    wt = surf.wt_lunit[gi, ltype]
    (wt > 0.0 || varctl.run_zero_weight_urban) || return nothing

    add_landunit!(lun, li, gi, ltype, wt)

    # Default urban column weights
    wtlunit_roof = 0.5
    wtroad_perv = 0.5

    for m in 1:varpar.maxpatch_urb
        if m == 1
            ctype = ICOL_ROOF
            wtcol = wtlunit_roof
        elseif m == 2
            ctype = ICOL_SUNWALL
            wtcol = (1.0 - wtlunit_roof) / 3.0
        elseif m == 3
            ctype = ICOL_SHADEWALL
            wtcol = (1.0 - wtlunit_roof) / 3.0
        elseif m == 4
            ctype = ICOL_ROAD_IMPERV
            wtcol = ((1.0 - wtlunit_roof) / 3.0) * (1.0 - wtroad_perv)
        else
            ctype = ICOL_ROAD_PERV
            wtcol = ((1.0 - wtlunit_roof) / 3.0) * wtroad_perv
        end

        add_column!(col, lun, ci, li[], ctype, wtcol)
        add_patch!(pch, col, lun, pi, ci[], noveg, 1.0)
    end

    nothing
end

function set_landunit_wet_lake!(surf::SurfaceInputData, ltype::Int, gi::Int,
                                 li::Ref{Int}, ci::Ref{Int}, pi::Ref{Int},
                                 lun::LandunitData, col::ColumnData, pch::PatchData)
    wt = surf.wt_lunit[gi, ltype]

    # Always create lake landunit, skip wetland if zero weight
    if ltype == ISTDLAK
        # Always create lake
    elseif wt <= 0.0
        return nothing
    end

    add_landunit!(lun, li, gi, ltype, wt)
    add_column!(col, lun, ci, li[], ltype, 1.0)
    add_patch!(pch, col, lun, pi, ci[], noveg, 1.0)

    nothing
end

function set_landunit_ice!(surf::SurfaceInputData, gi::Int,
                            li::Ref{Int}, ci::Ref{Int}, pi::Ref{Int},
                            lun::LandunitData, col::ColumnData, pch::PatchData)
    wt = surf.wt_lunit[gi, ISTICE]
    wt > 0.0 || return nothing

    add_landunit!(lun, li, gi, ISTICE, wt)

    maxpatch_glc = varpar.maxpatch_glc
    n_added = 0
    if size(surf.wt_glc_mec, 2) > 0
        for m in 1:min(maxpatch_glc, size(surf.wt_glc_mec, 2))
            glc_wt = surf.wt_glc_mec[gi, m]
            if glc_wt > 0.0
                add_column!(col, lun, ci, li[], ice_class_to_col_itype(m), glc_wt)
                add_patch!(pch, col, lun, pi, ci[], noveg, 1.0)
                n_added += 1
            end
        end
    end
    if n_added == 0
        add_column!(col, lun, ci, li[], ice_class_to_col_itype(1), 1.0)
        add_patch!(pch, col, lun, pi, ci[], noveg, 1.0)
    end

    nothing
end
