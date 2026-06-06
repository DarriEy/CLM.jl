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
# GPU helper: move a host filter (Vector{Int}) onto the backend of a device
# prototype array so kernels indexing device arrays can read it. On the plain
# CPU backend this is a no-op identity (byte-identical Vector{Int}); on a GPU
# backend it allocates a device Int vector and copyto!s the host filter in.
# --------------------------------------------------------------------------
@inline _to_device_filter(proto, filterc::Vector{Int}) =
    _to_device_filter(KA.get_backend(proto), proto, filterc)
@inline _to_device_filter(::KA.CPU, proto, filterc::Vector{Int}) = filterc
@inline function _to_device_filter(backend, proto, filterc::Vector{Int})
    d = similar(proto, Int, length(filterc))
    copyto!(d, filterc)
    return d
end

# --------------------------------------------------------------------------
# Default / hydstress-roads per-column kernel.
#
# Whole-function-on-device version of the transpiration-weighted column root
# fraction + vertical sink. Launched over `length(filterc)` columns; `c =
# filterc[i]`. Each column owns a disjoint, contiguous patch range
# [patchi[c], patchi[c]+npatches[c]-1], so the per-patch accumulation needs no
# atomics — it is an internal SEQUENTIAL reduction owned by this thread. The
# internal layer loop (j = 1..nlevsoi) and patch loop mirror the Fortran exactly:
#   rootr_col(c,j) = sum_p active(p) * rootr_patch(p,j)*qtran(p)*wtcol(p)
#   temp(c)        = sum_p active(p) * qtran(p)*wtcol(p)
#   rootr_col(c,j) /= temp(c)   (iff temp != 0)
#   qflx_rootsoi_col(c,j) = rootr_col(c,j) * qtran_col(c)
# Literals carry the working eltype (`zero(T)`) so no Float64 reaches a Float32
# backend; byte-identical on Float64 CPU.
# --------------------------------------------------------------------------
@kernel function _plantsink_default_kernel!(rootr_col, qflx_rootsoi_col,
        @Const(filterc), @Const(rootr_patch), @Const(qflx_tran_veg_patch),
        @Const(qflx_tran_veg_col), @Const(patchi), @Const(npatches),
        @Const(active), @Const(wtcol), nlevsoi::Int)
    i = @index(Global)
    @inbounds begin
        c = filterc[i]
        T = eltype(rootr_col)
        p_lo = patchi[c]
        p_hi = patchi[c] + npatches[c] - 1

        # Denominator: transpiration-weighted area sum over active patches.
        temp_c = zero(T)
        for p in p_lo:p_hi
            if active[p]
                temp_c += qflx_tran_veg_patch[p] * wtcol[p]
            end
        end

        for j in 1:nlevsoi
            acc = zero(T)
            for p in p_lo:p_hi
                if active[p]
                    acc += rootr_patch[p, j] * qflx_tran_veg_patch[p] * wtcol[p]
                end
            end
            if temp_c != zero(T)
                acc /= temp_c
            end
            rootr_col[c, j] = acc
            qflx_rootsoi_col[c, j] = acc * qflx_tran_veg_col[c]
        end
    end
end

function _plantsink_default_launch!(rootr_col, qflx_rootsoi_col, filterc, rootr_patch,
        qflx_tran_veg_patch, qflx_tran_veg_col, patchi, npatches, active, wtcol,
        nlevsoi::Int)
    _launch!(_plantsink_default_kernel!, rootr_col, qflx_rootsoi_col, filterc,
             rootr_patch, qflx_tran_veg_patch, qflx_tran_veg_col, patchi, npatches,
             active, wtcol, nlevsoi; ndrange = length(filterc))
end

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

    # Move the host filter onto the backend of the written matrix (no-op on CPU).
    filterc_d = _to_device_filter(rootr_col, filterc)

    # Whole-column kernel: zero + weighted accumulation + normalize + sink, with
    # internal sequential layer/patch loops (each column owns its patch range).
    _plantsink_default_launch!(rootr_col, qflx_rootsoi_col, filterc_d, rootr_patch,
        qflx_tran_veg_patch, qflx_tran_veg_col, col.patchi, col.npatches,
        patch_data.active, patch_data.wtcol, nlevsoi)

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

    # Identical algorithm to the default method (see _plantsink_default_kernel!).
    filterc_d = _to_device_filter(rootr_col, filterc)
    _plantsink_default_launch!(rootr_col, qflx_rootsoi_col, filterc_d, rootr_patch,
        qflx_tran_veg_patch, qflx_tran_veg_col, col.patchi, col.npatches,
        patch_data.active, patch_data.wtcol, nlevsoi)

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

# --------------------------------------------------------------------------
# Hydraulic-stress per-column kernel.
#
# Whole-function-on-device version. Launched over `length(filterc)` columns;
# `c = filterc[i]`. Each column owns a disjoint patch range, so the per-patch
# reductions (column flux temp_c, and the per-patch hydraulic-redistribution
# accumulation qflx_hydr_redist_patch[p] across layers) are SEQUENTIAL within
# this thread — no atomics. Internal sequential loops over layers j and patches
# p mirror the Fortran exactly. `nlevsno_off` is the snow-layer offset for z
# indexing (col.z is (ncols, nlevsno + nlevmaxurbgrnd)). Literals carry the
# working eltype (`zero(T)`, `T(1000)`) so no Float64 reaches a Float32 backend;
# byte-identical on Float64 CPU.
# --------------------------------------------------------------------------
@kernel function _plantsink_hydstress_kernel!(qflx_rootsoi_col, qflx_phs_neg_col,
        qflx_hydr_redist_patch, @Const(filterc), @Const(k_soil_root), @Const(smp),
        @Const(vegwp), @Const(frac_veg_nosno), @Const(z), @Const(patchi),
        @Const(npatches), @Const(active), @Const(wtcol), nlevsno_off::Int,
        nlevsoi::Int, root_idx::Int)
    i = @index(Global)
    @inbounds begin
        c = filterc[i]
        T = eltype(qflx_rootsoi_col)
        p_lo = patchi[c]
        p_hi = patchi[c] + npatches[c] - 1

        qflx_phs_neg_col[c] = zero(T)

        for j in 1:nlevsoi
            grav2 = z[c, j + nlevsno_off] * T(1000)
            temp_c = zero(T)

            for p in p_lo:p_hi
                if j == 1
                    qflx_hydr_redist_patch[p] = zero(T)
                end
                if active[p] && frac_veg_nosno[p] > 0
                    if wtcol[p] > zero(T)
                        patchflux = k_soil_root[p, j] * (smp[c, j] - vegwp[p, root_idx] - grav2)
                        if patchflux < zero(T)
                            qflx_hydr_redist_patch[p] += patchflux
                        end
                        temp_c += patchflux * wtcol[p]
                    end
                end
            end
            qflx_rootsoi_col[c, j] = temp_c

            if temp_c < zero(T)
                qflx_phs_neg_col[c] += temp_c
            end
        end
    end
end

function _plantsink_hydstress_launch!(qflx_rootsoi_col, qflx_phs_neg_col,
        qflx_hydr_redist_patch, filterc, k_soil_root, smp, vegwp, frac_veg_nosno,
        z, patchi, npatches, active, wtcol, nlevsno_off::Int, nlevsoi::Int,
        root_idx::Int)
    _launch!(_plantsink_hydstress_kernel!, qflx_rootsoi_col, qflx_phs_neg_col,
             qflx_hydr_redist_patch, filterc, k_soil_root, smp, vegwp, frac_veg_nosno,
             z, patchi, npatches, active, wtcol, nlevsno_off, nlevsoi, root_idx;
             ndrange = length(filterc))
end

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

    # Move the host filter onto the backend of the written matrix (no-op on CPU).
    filterc_d = _to_device_filter(qflx_rootsoi_col, filterc)

    # Per-column kernel: internal sequential layer/patch loops, each column owns
    # its patch range (the hydr-redist accumulation is per-patch, no atomics).
    _plantsink_hydstress_launch!(qflx_rootsoi_col, qflx_phs_neg_col,
        qflx_hydr_redist_patch, filterc_d, k_soil_root, smp, vegwp, frac_veg_nosno,
        z, col.patchi, col.npatches, patch_data.active, patch_data.wtcol,
        nlevsno_off, nlevsoi, root_idx)

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
        mask_hydrology::AbstractVector{Bool},
        soilstate_inst::SoilStateData,
        canopystate_inst::CanopyStateData,
        waterfluxbulk_inst::WaterFluxBulkData,
        energyflux_inst::EnergyFluxData,
        col::ColumnData,
        lun::LandunitData,
        patch_data::PatchData;
        use_hydrstress::Bool = false)

    num_filterc_tot = 0

    # Host copies of the small integer metadata the column-group filters are
    # built from. On the CPU these are identity (the arrays are already host
    # Arrays); on a GPU backend this pulls the column/landunit type vectors to
    # the host once so filter construction stays scalar-index-free on the device.
    # The resulting Int filters are then moved back to the device per sub-call.
    col_itype    = Array(col.itype)
    col_landunit = Array(col.landunit)
    lun_itype    = Array(lun.itype)
    # Pull the hydrology mask to the host too so the column-group filter loops
    # below stay scalar-index-free on a device backend (no-op on CPU).
    mask_host    = Array(mask_hydrology)

    # ---------------------------------------------------------------
    # 1) Pervious roads
    # ---------------------------------------------------------------
    filterc_roads = Int[]
    for c in bounds_col
        mask_host[c] || continue
        if col_itype[c] == ICOL_ROAD_PERV
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
        mask_host[c] || continue
        l = col_landunit[c]
        if col_itype[c] != ICOL_ROAD_PERV && lun_itype[l] != ISTSOIL
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
        mask_host[c] || continue
        l = col_landunit[c]
        if lun_itype[l] == ISTSOIL
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
    num_hydrology = count(mask_host[bounds_col])
    if num_hydrology != num_filterc_tot
        error("SoilWaterPlantSink: total columns flagged for root water uptake " *
              "($num_filterc_tot) did not match the total hydrology filter count " *
              "($num_hydrology). This is likely a problem with column/landunit filters.")
    end

    return nothing
end
