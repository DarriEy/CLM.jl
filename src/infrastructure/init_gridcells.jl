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

# Per-gridcell lat/lon
For a multi-gridcell run each gridcell carries its own latitude/longitude. These
are taken (in priority order) from:
  1. `lats`/`lons` keyword vectors (length ng), if supplied;
  2. `surf.grid_latdeg`/`surf.grid_londeg` (populated by the 2D surfdata read);
  3. the scalar `lat`/`lon` fallback (legacy single-gridcell / tiled path),
applied identically to every gridcell.
"""
function initGridCells!(bounds::BoundsType, surf::SurfaceInputData,
                         grc::GridcellData, lun::LandunitData,
                         col::ColumnData, pch::PatchData;
                         lat::Float64=0.0, lon::Float64=0.0, area::Float64=1.0,
                         lats::Union{Nothing,AbstractVector}=nothing,
                         lons::Union{Nothing,AbstractVector}=nothing)
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

    # Set gridcell properties. Per-gridcell lat/lon priority:
    #   lats/lons kwargs > surf.grid_latdeg/londeg > scalar lat/lon fallback.
    have_surf_ll = length(surf.grid_latdeg) >= bounds.endg &&
                   length(surf.grid_londeg) >= bounds.endg
    for g in bounds.begg:bounds.endg
        glat = lats !== nothing ? Float64(lats[g]) :
               (have_surf_ll ? surf.grid_latdeg[g] : lat)
        glon = lons !== nothing ? Float64(lons[g]) :
               (have_surf_ll ? surf.grid_londeg[g] : lon)
        grc.area[g] = area
        grc.latdeg[g] = glat
        grc.londeg[g] = glon
        grc.lat[g] = glat * RPI / 180.0
        grc.lon[g] = glon * RPI / 180.0
    end

    # Fill down-pointers
    clm_ptrs_compdown!(bounds, grc, lun, col, pch)

    # Validate
    clm_ptrs_check!(bounds, grc, lun, col, pch)

    # Compute higher-order weights
    compute_higher_order_weights!(bounds, col, lun, pch)

    nothing
end

"""
    init_hillslope_columns!(bounds, surf, grc, lun, col, pch)

Populate the hillslope catena geometry on the subgrid hierarchy AFTER
`initGridCells!` has built the multi-column natural-veg landunits. Bridges the
gridcell-indexed hillslope arrays read into `SurfaceInputData` (§`surfrd_hillslope!`)
to the landunit-indexed arrays `init_hillslope!` expects, then calls
`init_hillslope!` (sets `col.cold`/`col.colu`/`col.hill_*`/`is_hillslope_column`
and — under routing — the landunit `stream_channel_*`) and
`set_hillslope_soil_thickness!`. Finally recomputes higher-order weights so
`col.wtgcell`/patch weights reflect the per-column hillslope areas that
`init_hillslope!` wrote into `col.wtlunit`.

Mirrors CTSM's `clm_instInit` → `InitHillslope` call. Strictly gated: a no-op
unless `varctl.use_hillslope` is set AND `surf` carries hillslope geometry, so
the default (single-column) path is untouched.
"""
function init_hillslope_columns!(bounds::BoundsType, surf::SurfaceInputData,
                                 grc::GridcellData, lun::LandunitData,
                                 col::ColumnData, pch::PatchData)
    varctl.use_hillslope || return nothing
    isempty(surf.nhillcolumns) && return nothing

    nhill = surf.nhillslope_dim
    nmax  = surf.nmaxhillcol_dim
    (nhill > 0 && nmax > 0) || return nothing

    nl = bounds.endl
    ncolumns_hillslope = zeros(Int, nl)
    pct_hillslope = zeros(nl, nhill)
    hill_ndx   = zeros(Int, nl, nmax)
    col_ndx_in = zeros(Int, nl, nmax)
    col_dndx   = fill(-999, nl, nmax)
    hill_slope   = zeros(nl, nmax)
    hill_aspect  = zeros(nl, nmax)
    hill_area    = zeros(nl, nmax)
    hill_dist    = zeros(nl, nmax)
    hill_width   = zeros(nl, nmax)
    hill_elev    = zeros(nl, nmax)

    use_routing = varctl.use_hillslope_routing && !isempty(surf.hillslope_stream_depth)
    stream_depth = use_routing ? zeros(nl) : nothing
    stream_width = use_routing ? zeros(nl) : nothing
    stream_slope = use_routing ? zeros(nl) : nothing

    for l in bounds.begl:bounds.endl
        lun.itype[l] == ISTSOIL || continue
        g = lun.gridcell[l]
        (g >= 1 && g <= length(surf.nhillcolumns)) || continue
        ncolumns_hillslope[l] = max(surf.nhillcolumns[g], 0)
        for nh in 1:nhill
            pct_hillslope[l, nh] = surf.pct_hillslope[g, nh]
        end
        for ci in 1:nmax
            hill_ndx[l, ci]   = surf.hillslope_index[g, ci]
            col_ndx_in[l, ci] = surf.column_index[g, ci]
            col_dndx[l, ci]   = surf.downhill_column_index[g, ci]
            hill_slope[l, ci] = surf.hillslope_slope[g, ci]
            hill_aspect[l, ci]= surf.hillslope_aspect[g, ci]
            hill_area[l, ci]  = surf.hillslope_area[g, ci]
            hill_dist[l, ci]  = surf.hillslope_distance[g, ci]
            hill_width[l, ci] = surf.hillslope_width[g, ci]
            hill_elev[l, ci]  = surf.hillslope_elevation[g, ci]
        end
        if use_routing
            stream_depth[l] = surf.hillslope_stream_depth[g]
            stream_width[l] = surf.hillslope_stream_width[g]
            stream_slope[l] = surf.hillslope_stream_slope[g]
        end
    end

    init_hillslope!(col, lun, grc, bounds.begl:bounds.endl, bounds.begc:bounds.endc;
        nhillslope = nhill, max_columns_hillslope = nmax,
        ncolumns_hillslope = ncolumns_hillslope, pct_hillslope = pct_hillslope,
        hill_ndx = hill_ndx, col_ndx_in = col_ndx_in, col_dndx = col_dndx,
        hill_slope_in = hill_slope, hill_aspect_in = hill_aspect,
        hill_area_in = hill_area, hill_dist_in = hill_dist,
        hill_width_in = hill_width, hill_elev_in = hill_elev,
        stream_channel_depth = stream_depth, stream_channel_width = stream_width,
        stream_channel_slope = stream_slope,
        use_hillslope_routing = use_routing)

    set_hillslope_soil_thickness!(col, lun, bounds.begc:bounds.endc;
        soil_profile_method = hillslope_config[].soil_profile_method,
        nlevsoi = varpar.nlevsoi)

    # InitHillslope rewrote col.wtlunit from the per-column hillslope areas; push
    # that down into col.wtgcell and the patch weights.
    compute_higher_order_weights!(bounds, col, lun, pch)

    return nothing
end

# --- Private per-landunit builders ---

function set_landunit_veg_compete!(surf::SurfaceInputData, gi::Int,
                                    li::Ref{Int}, ci::Ref{Int}, pi::Ref{Int},
                                    lun::LandunitData, col::ColumnData, pch::PatchData)
    wt = surf.wt_lunit[gi, ISTSOIL]
    natpft_size = size(surf.wt_nat_patch, 2)

    # Always create soil landunit
    add_landunit!(lun, li, gi, ISTSOIL, wt)

    # Number of columns for this natural-veg landunit. Default is a single column;
    # under hillslope hydrology the landunit is a CATENA of `ncolumns_hillslope`
    # connected columns (CTSM subgrid_get_info_natveg / set_landunit_veg_compete:
    # ncols = max(ncolumns_hillslope(gi),1), each with the SAME PFT breakdown, and
    # column weight 1/ncols set here — later overwritten by InitHillslope from the
    # per-column hillslope areas). Gated so the default (use_hillslope=false or no
    # hillslope geometry) builds exactly one column of weight 1.0, byte-identical.
    ncols_hill = 1
    have_hill = varctl.use_hillslope && !isempty(surf.nhillcolumns) &&
                gi <= length(surf.nhillcolumns)
    if have_hill
        ncols_hill = max(surf.nhillcolumns[gi], 1)
    end
    wtcol2lunit = 1.0 / ncols_hill

    for ci2 in 1:ncols_hill
        add_column!(col, lun, ci, li[], 1, wtcol2lunit)
        # Record the external (surfdata) column index so InitHillslope can map the
        # downhill-neighbour pointers. For the single-column default this is 1 and
        # never read (is_hillslope_column stays false).
        if have_hill && !isempty(surf.column_index) && ci2 <= size(surf.column_index, 2)
            col.col_ndx[ci[]] = surf.column_index[gi, ci2]
        end

        # Add patches for each PFT with non-zero weight — IDENTICAL for the non-FATES
        # and FATES paths, so the initial stand + patch weights (and hence the
        # dynamics) are unchanged. For FATES we then PAD the column with extra inactive
        # vegetated patch slots so its disturbance patches have HLM slots to map onto
        # (up to fates_maxpatch + 1 total) instead of being dropped from the coupling.
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
            n_added = 1
        end
        if varctl.use_fates && varctl.fates_maxpatch > 0
            # Pad with inactive vegetated slots (weight 0 — FATES sets the real
            # weight/type from its patch areas when disturbance populates them).
            # Matches count_subgrid_elements.
            for _ in (n_added + 1):(varctl.fates_maxpatch + 1)
                add_patch!(pch, col, lun, pi, ci[], varpar.natpft_lb + 1, 0.0)
            end
        end
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
