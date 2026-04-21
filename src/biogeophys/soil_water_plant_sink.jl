# ==========================================================================
# Ported from: src/biogeophys/SoilWaterPlantSinkMod.F90
# Calculates effective root fraction and vertical transpiration sink.
#
# This module computes the soil water sink due to plant transpiration,
# distributing it vertically over the soil column. There are two main
# pathways:
#   1) Default method: weighted column-level root fractions from patch data
#   2) Hydraulic stress method: flux driven by soil-root conductance and
#      water potential gradients
#
# Two output quantities:
#   waterfluxbulk_inst.qflx_rootsoi_col(c,j)  -- root-soil water exchange
#   soilstate_inst.rootr_col(c,j)              -- column effective root fraction
# ==========================================================================

# --------------------------------------------------------------------------
# compute_effec_rootfrac_and_vert_tran_sink_default!
# --------------------------------------------------------------------------

"""
    compute_effec_rootfrac_and_vert_tran_sink_default!(
        bounds_col, nlevsoi, filterc,
        soilstate_inst, waterfluxbulk_inst,
        col, patch_data)

Default method: compute column-level effective root fraction and
vertical transpiration sink from patch-level root fractions weighted
by transpiration.

Applicable to any Richards solver not coupled to plant hydraulics.

The column-level effective root fraction is NOT a simple area-weighted
average of patch-level rootr. Instead, the weighting depends on both
the per-unit-area transpiration of the patch and the patch's area
relative to all patches.

Ported from `Compute_EffecRootFrac_And_VertTranSink_Default`
in `SoilWaterPlantSinkMod.F90`.
"""
function compute_effec_rootfrac_and_vert_tran_sink_default!(
        bounds_col::UnitRange{Int},
        nlevsoi::Int,
        filterc::Vector{Int},
        soilstate_inst::SoilStateData,
        waterfluxbulk_inst::WaterFluxBulkData,
        col::ColumnData,
        patch_data::PatchData)

    # Aliases (following Fortran associate block)
    qflx_rootsoi_col    = waterfluxbulk_inst.qflx_rootsoi_col
    qflx_tran_veg_patch = waterfluxbulk_inst.wf.qflx_tran_veg_patch
    qflx_tran_veg_col   = waterfluxbulk_inst.wf.qflx_tran_veg_col
    rootr_patch          = soilstate_inst.rootr_patch
    rootr_col            = soilstate_inst.rootr_col

    # Accumulator for transpiration-weighted denominator
    FT = eltype(rootr_patch)
    temp = zeros(FT, length(bounds_col))
    # Map column index to temp array index
    c_offset = first(bounds_col) - 1

    # Zero out rootr_col for filter columns
    for j in 1:nlevsoi
        for c in filterc
            rootr_col[c, j] = 0.0
        end
    end

    # Accumulate transpiration-weighted root fraction
    for j in 1:nlevsoi
        for c in filterc
            for p in col.patchi[c]:(col.patchi[c] + col.npatches[c] - 1)
                if patch_data.active[p]
                    rootr_col[c, j] += rootr_patch[p, j] *
                        qflx_tran_veg_patch[p] * patch_data.wtcol[p]
                end
            end
        end
    end

    # Accumulate denominator: sum of (transpiration * weight) over active patches
    for c in filterc
        for p in col.patchi[c]:(col.patchi[c] + col.npatches[c] - 1)
            if patch_data.active[p]
                temp[c - c_offset] += qflx_tran_veg_patch[p] * patch_data.wtcol[p]
            end
        end
    end

    # Normalize rootr_col and compute qflx_rootsoi_col
    for j in 1:nlevsoi
        for c in filterc
            if temp[c - c_offset] != 0.0
                rootr_col[c, j] /= temp[c - c_offset]
            end
            qflx_rootsoi_col[c, j] = rootr_col[c, j] * qflx_tran_veg_col[c]
        end
    end

    return nothing
end

# --------------------------------------------------------------------------
# compute_effec_rootfrac_and_vert_tran_sink_hydstress_roads!
# --------------------------------------------------------------------------

"""
    compute_effec_rootfrac_and_vert_tran_sink_hydstress_roads!(
        bounds_col, nlevsoi, filterc,
        soilstate_inst, waterfluxbulk_inst,
        col, patch_data)

Hydraulic stress method for pervious road columns. Identical logic to
the default method (transpiration-weighted column root fractions).

Ported from `Compute_EffecRootFrac_And_VertTranSink_HydStress_Roads`
in `SoilWaterPlantSinkMod.F90`.
"""
function compute_effec_rootfrac_and_vert_tran_sink_hydstress_roads!(
        bounds_col::UnitRange{Int},
        nlevsoi::Int,
        filterc::Vector{Int},
        soilstate_inst::SoilStateData,
        waterfluxbulk_inst::WaterFluxBulkData,
        col::ColumnData,
        patch_data::PatchData)

    # Aliases
    qflx_rootsoi_col    = waterfluxbulk_inst.qflx_rootsoi_col
    qflx_tran_veg_patch = waterfluxbulk_inst.wf.qflx_tran_veg_patch
    qflx_tran_veg_col   = waterfluxbulk_inst.wf.qflx_tran_veg_col
    rootr_patch          = soilstate_inst.rootr_patch
    rootr_col            = soilstate_inst.rootr_col

    # Accumulator for transpiration-weighted denominator
    FT = eltype(rootr_patch)
    temp = zeros(FT, length(bounds_col))
    c_offset = first(bounds_col) - 1

    # Zero out rootr_col for filter columns
    for j in 1:nlevsoi
        for c in filterc
            rootr_col[c, j] = 0.0
        end
    end

    # Accumulate transpiration-weighted root fraction
    for j in 1:nlevsoi
        for c in filterc
            for p in col.patchi[c]:(col.patchi[c] + col.npatches[c] - 1)
                if patch_data.active[p]
                    rootr_col[c, j] += rootr_patch[p, j] *
                        qflx_tran_veg_patch[p] * patch_data.wtcol[p]
                end
            end
        end
    end

    # Accumulate denominator
    for c in filterc
        for p in col.patchi[c]:(col.patchi[c] + col.npatches[c] - 1)
            if patch_data.active[p]
                temp[c - c_offset] += qflx_tran_veg_patch[p] * patch_data.wtcol[p]
            end
        end
    end

    # Normalize rootr_col and compute qflx_rootsoi_col
    for j in 1:nlevsoi
        for c in filterc
            if temp[c - c_offset] != 0.0
                rootr_col[c, j] /= temp[c - c_offset]
            end
            qflx_rootsoi_col[c, j] = rootr_col[c, j] * qflx_tran_veg_col[c]
        end
    end

    return nothing
end

# --------------------------------------------------------------------------
# compute_effec_rootfrac_and_vert_tran_sink_hydstress!
# --------------------------------------------------------------------------

"""
    compute_effec_rootfrac_and_vert_tran_sink_hydstress!(
        bounds_col, nlevsoi, filterc,
        waterfluxbulk_inst, soilstate_inst,
        canopystate_inst, col, patch_data)

Hydraulic stress method for vegetation columns. Computes root-soil water
exchange from the soil-root interface conductance and the water potential
gradient between soil and root.

    patchflux = k_soil_root(p,j) * (smp(c,j) - vegwp(p,root) - grav2)

where `grav2 = z(c,j) * 1000` is the gravitational potential in mm H2O.

Negative patchflux (water flowing from root to soil) is accumulated into
`qflx_hydr_redist_patch` and `qflx_phs_neg_col`.

Ported from `Compute_EffecRootFrac_And_VertTranSink_HydStress`
in `SoilWaterPlantSinkMod.F90`.
"""
function compute_effec_rootfrac_and_vert_tran_sink_hydstress!(
        bounds_col::UnitRange{Int},
        nlevsoi::Int,
        filterc::Vector{Int},
        waterfluxbulk_inst::WaterFluxBulkData,
        soilstate_inst::SoilStateData,
        canopystate_inst::CanopyStateData,
        col::ColumnData,
        patch_data::PatchData)

    # Index for root water potential segment (Fortran: root=4)
    root_idx = 4

    # Aliases
    k_soil_root            = soilstate_inst.k_soil_root_patch
    qflx_phs_neg_col       = waterfluxbulk_inst.qflx_phs_neg_col
    qflx_hydr_redist_patch = waterfluxbulk_inst.qflx_hydr_redist_patch
    qflx_tran_veg_col      = waterfluxbulk_inst.wf.qflx_tran_veg_col
    qflx_tran_veg_patch    = waterfluxbulk_inst.wf.qflx_tran_veg_patch
    qflx_rootsoi_col       = waterfluxbulk_inst.qflx_rootsoi_col
    smp                    = soilstate_inst.smp_l_col
    frac_veg_nosno         = canopystate_inst.frac_veg_nosno_patch
    z                      = col.z
    vegwp                  = canopystate_inst.vegwp_patch

    # Snow layer offset for z indexing.
    # col.z is stored as (ncols, nlevsno + nlevmaxurbgrnd). Soil layer j
    # corresponds to index j + nlevsno in the combined array.
    nlevsno_off = varpar.nlevsno

    for c in filterc
        qflx_phs_neg_col[c] = 0.0

        for j in 1:nlevsoi
            grav2 = z[c, j + nlevsno_off] * 1000.0
            temp_c = 0.0

            for p in col.patchi[c]:(col.patchi[c] + col.npatches[c] - 1)
                if j == 1
                    qflx_hydr_redist_patch[p] = 0.0
                end
                if patch_data.active[p] && frac_veg_nosno[p] > 0
                    if patch_data.wtcol[p] > 0.0
                        patchflux = k_soil_root[p, j] * (smp[c, j] - vegwp[p, root_idx] - grav2)
                        if patchflux < 0.0
                            qflx_hydr_redist_patch[p] += patchflux
                        end
                        temp_c += patchflux * patch_data.wtcol[p]
                    end
                end
            end
            qflx_rootsoi_col[c, j] = temp_c

            if temp_c < 0.0
                qflx_phs_neg_col[c] += temp_c
            end
        end
    end

    return nothing
end

# --------------------------------------------------------------------------
# compute_effec_rootfrac_and_vert_tran_sink!
# --------------------------------------------------------------------------

"""
    compute_effec_rootfrac_and_vert_tran_sink!(
        bounds_col, nlevsoi, mask_hydrology,
        soilstate_inst, canopystate_inst,
        waterfluxbulk_inst, energyflux_inst,
        col, lun, patch_data;
        use_hydrstress=false)

Wrapper that dispatches to the correct method for computing the
effective root fraction and soil water sink due to plant transpiration.

Three groups of columns are processed:
  1. Pervious roads (icol_road_perv)
  2. Non-natural vegetation (not pervious road and not istsoil)
  3. Natural vegetation (istsoil)

For each group, the appropriate method is called based on `use_hydrstress`.

Ported from `Compute_EffecRootFrac_And_VertTranSink`
in `SoilWaterPlantSinkMod.F90`.
"""
function compute_effec_rootfrac_and_vert_tran_sink!(
        bounds_col::UnitRange{Int},
        nlevsoi::Int,
        mask_hydrology::BitVector,
        soilstate_inst::SoilStateData,
        canopystate_inst::CanopyStateData,
        waterfluxbulk_inst::WaterFluxBulkData,
        energyflux_inst::EnergyFluxData,
        col::ColumnData,
        lun::LandunitData,
        patch_data::PatchData;
        use_hydrstress::Bool = false)

    num_filterc_tot = 0

    # ---------------------------------------------------------------
    # 1) Pervious roads
    # ---------------------------------------------------------------
    filterc_roads = Int[]
    for c in bounds_col
        mask_hydrology[c] || continue
        if col.itype[c] == ICOL_ROAD_PERV
            push!(filterc_roads, c)
        end
    end
    num_filterc_tot += length(filterc_roads)

    if use_hydrstress
        compute_effec_rootfrac_and_vert_tran_sink_hydstress_roads!(
            bounds_col, nlevsoi, filterc_roads,
            soilstate_inst, waterfluxbulk_inst,
            col, patch_data)
    else
        compute_effec_rootfrac_and_vert_tran_sink_default!(
            bounds_col, nlevsoi, filterc_roads,
            soilstate_inst, waterfluxbulk_inst,
            col, patch_data)
    end

    # ---------------------------------------------------------------
    # 2) Not pervious road and not natural vegetation (everything else)
    # ---------------------------------------------------------------
    filterc_other = Int[]
    for c in bounds_col
        mask_hydrology[c] || continue
        l = col.landunit[c]
        if col.itype[c] != ICOL_ROAD_PERV && lun.itype[l] != ISTSOIL
            push!(filterc_other, c)
        end
    end
    num_filterc_tot += length(filterc_other)

    if use_hydrstress
        compute_effec_rootfrac_and_vert_tran_sink_hydstress!(
            bounds_col, nlevsoi, filterc_other,
            waterfluxbulk_inst, soilstate_inst,
            canopystate_inst, col, patch_data)
    else
        compute_effec_rootfrac_and_vert_tran_sink_default!(
            bounds_col, nlevsoi, filterc_other,
            soilstate_inst, waterfluxbulk_inst,
            col, patch_data)
    end

    # ---------------------------------------------------------------
    # 3) Natural vegetation (istsoil)
    # ---------------------------------------------------------------
    filterc_natveg = Int[]
    for c in bounds_col
        mask_hydrology[c] || continue
        l = col.landunit[c]
        if lun.itype[l] == ISTSOIL
            push!(filterc_natveg, c)
        end
    end
    num_filterc_tot += length(filterc_natveg)

    if use_hydrstress
        compute_effec_rootfrac_and_vert_tran_sink_hydstress!(
            bounds_col, nlevsoi, filterc_natveg,
            waterfluxbulk_inst, soilstate_inst,
            canopystate_inst, col, patch_data)
    else
        compute_effec_rootfrac_and_vert_tran_sink_default!(
            bounds_col, nlevsoi, filterc_natveg,
            soilstate_inst, waterfluxbulk_inst,
            col, patch_data)
    end

    # Sanity check: all hydrology columns must have been processed
    num_hydrology = count(mask_hydrology[bounds_col])
    if num_hydrology != num_filterc_tot
        error("SoilWaterPlantSink: total columns flagged for root water uptake " *
              "($num_filterc_tot) did not match the total hydrology filter count " *
              "($num_hydrology). This is likely a problem with column/landunit filters.")
    end

    return nothing
end
