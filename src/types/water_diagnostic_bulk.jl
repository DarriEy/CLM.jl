# ==========================================================================
# Ported from: src/biogeophys/WaterDiagnosticBulkType.F90
# Water diagnostic variables that apply only to bulk water.
# Also includes fields from the parent WaterDiagnosticType.F90.
# ==========================================================================

# --- Module-level parameters (from params_type in Fortran) ---

"""
    WaterDiagnosticBulkParams

Parameters read from parameter file for water diagnostic bulk type.

Ported from `params_type` in `WaterDiagnosticBulkType.F90`.
"""
Base.@kwdef mutable struct WaterDiagnosticBulkParams{FT<:Real}
    zlnd::FT = 0.01           # Momentum roughness length for soil, glacier, wetland (m)
    snw_rds_min::FT = 54.526  # minimum allowed snow effective radius (also cold "fresh snow" value) [microns]
end

const waterdiagbulk_params = WaterDiagnosticBulkParams()

"""
    WaterDiagnosticBulkData

Water diagnostic bulk data structure. Holds diagnostic water variables at the
column, patch, landunit, and gridcell levels for bulk water.

Includes fields from both the parent `waterdiagnostic_type` and the child
`waterdiagnosticbulk_type` in the Fortran source.

Ported from `waterdiagnosticbulk_type` in `WaterDiagnosticBulkType.F90` and
`waterdiagnostic_type` in `WaterDiagnosticType.F90`.
"""
Base.@kwdef mutable struct WaterDiagnosticBulkData{FT<:Real,
                                  V<:AbstractVector{FT}, M<:AbstractMatrix{FT}}
    # =====================================================================
    # Fields from parent waterdiagnostic_type (WaterDiagnosticType.F90)
    # =====================================================================

    # --- Column-level 1D (parent) ---
    snowice_col                 ::V = Float64[]   # col average snow ice lens
    snowliq_col                 ::V = Float64[]   # col average snow liquid water
    total_plant_stored_h2o_col  ::V = Float64[]   # col water bound in plants (kg/m2 H2O)
    h2osoi_liqice_10cm_col      ::V = Float64[]   # col liquid water + ice lens in top 10cm of soil (kg/m2)
    qg_snow_col                 ::V = Float64[]   # col ground specific humidity [kg/kg]
    qg_soil_col                 ::V = Float64[]   # col ground specific humidity [kg/kg]
    qg_h2osfc_col               ::V = Float64[]   # col ground specific humidity [kg/kg]
    qg_col                      ::V = Float64[]   # col ground specific humidity [kg/kg]

    # --- Patch-level 1D (parent) ---
    h2ocan_patch                ::V = Float64[]   # patch total canopy water (liq+ice) (mm H2O)
    q_ref2m_patch               ::V = Float64[]   # patch 2 m height surface specific humidity (kg/kg)

    # --- Landunit-level 1D (parent) ---
    qaf_lun                     ::V = Float64[]   # lun urban canopy air specific humidity (kg/kg)

    # --- Gridcell-level 1D (parent) ---
    tws_grc                     ::V = Float64[]   # grc total water storage (mm H2O)

    # =====================================================================
    # Fields from waterdiagnosticbulk_type (WaterDiagnosticBulkType.F90)
    # =====================================================================

    # --- Column-level 1D ---
    h2osno_total_col            ::V = Float64[]   # col total snow water (mm H2O)
    snow_depth_col              ::V = Float64[]   # col snow height of snow covered area (m)
    snow_5day_col               ::V = Float64[]   # col snow height 5 day avg (m)
    snowdp_col                  ::V = Float64[]   # col area-averaged snow height (m)
    snomelt_accum_col           ::V = Float64[]   # col accumulated snow melt for z0m calculation (m H2O)
    h2osoi_liq_tot_col          ::V = Float64[]   # col vertically summed liquid water (kg/m2)
    h2osoi_ice_tot_col          ::V = Float64[]   # col vertically summed ice lens (kg/m2)
    exice_subs_tot_col          ::V = Float64[]   # col total subsidence due to excess ice melt (m)
    exice_vol_tot_col           ::V = Float64[]   # col averaged volumetric excess ice content (m3/m3)
    snw_rds_top_col             ::V = Float64[]   # col snow grain radius (top layer) [microns]
    h2osno_top_col              ::V = Float64[]   # col top-layer mass of snow [kg]
    sno_liq_top_col             ::V = Float64[]   # col snow liquid water fraction (mass), top layer [fraction]
    dqgdT_col                   ::V = Float64[]   # col d(qg)/dT

    # Fractions (column 1D)
    frac_sno_col                ::V = Float64[]   # col fraction of ground covered by snow (0 to 1)
    frac_sno_eff_col            ::V = Float64[]   # col effective fraction of ground covered by snow (0 to 1)
    frac_h2osfc_col             ::V = Float64[]   # col fractional area with surface water > 0
    frac_h2osfc_nosnow_col      ::V = Float64[]   # col fractional area with surface water > 0 (if no snow)
    wf_col                      ::V = Float64[]   # col soil water as frac. of whc for top 0.05 m (0-1)
    wf2_col                     ::V = Float64[]   # col soil water as frac. of whc for top 0.17 m (0-1)

    # Summed fluxes (column 1D)
    qflx_prec_grnd_col          ::V = Float64[]   # col water onto ground including canopy runoff (mm H2O/s)

    # --- Column-level 2D ---
    snow_layer_unity_col        ::M = Matrix{Float64}(undef, 0, 0)  # col unity values for snow layers (-nlevsno+1:0)
    bw_col                      ::M = Matrix{Float64}(undef, 0, 0)  # col partial density of water in snow (kg/m3) (-nlevsno+1:0)
    air_vol_col                 ::M = Matrix{Float64}(undef, 0, 0)  # col air filled porosity (1:nlevgrnd)
    h2osoi_liqvol_col           ::M = Matrix{Float64}(undef, 0, 0)  # col volumetric liquid water content (v/v) (-nlevsno+1:nlevgrnd)
    swe_old_col                 ::M = Matrix{Float64}(undef, 0, 0)  # col initial snow water (-nlevsno+1:0)
    exice_subs_col              ::M = Matrix{Float64}(undef, 0, 0)  # col per layer subsidence due to excess ice melt (m) (1:nlevmaxurbgrnd)
    exice_vol_col               ::M = Matrix{Float64}(undef, 0, 0)  # col per layer volumetric excess ice content (m3/m3) (1:nlevsoi)
    snw_rds_col                 ::M = Matrix{Float64}(undef, 0, 0)  # col snow grain radius (col,lyr) [microns] (-nlevsno+1:0)
    frac_iceold_col             ::M = Matrix{Float64}(undef, 0, 0)  # col fraction of ice relative to tot water (-nlevsno+1:nlevgrnd)

    # --- Patch-level 1D ---
    iwue_ln_patch               ::V = Float64[]   # patch intrinsic water use efficiency near local noon (umolCO2/molH2O)
    vpd_ref2m_patch             ::V = Float64[]   # patch 2 m height surface vapor pressure deficit (Pa)
    rh_ref2m_patch              ::V = Float64[]   # patch 2 m height surface relative humidity (%)
    rh_ref2m_r_patch            ::V = Float64[]   # patch 2 m height surface relative humidity - rural (%)
    rh_ref2m_u_patch            ::V = Float64[]   # patch 2 m height surface relative humidity - urban (%)
    rh_af_patch                 ::V = Float64[]   # patch fractional humidity of canopy air (dimensionless)
    rh10_af_patch               ::V = Float64[]   # patch 10-day mean fractional humidity of canopy air (dimensionless)
    fwet_patch                  ::V = Float64[]   # patch canopy fraction that is wet (0 to 1)
    fcansno_patch               ::V = Float64[]   # patch canopy fraction that is snow covered (0 to 1)
    fdry_patch                  ::V = Float64[]   # patch canopy fraction of foliage that is green and dry [-]
    qflx_prec_intr_patch        ::V = Float64[]   # patch interception of precipitation (mm H2O/s)

    # --- Landunit-level 1D ---
    stream_water_depth_lun      ::V = Float64[]   # lun depth of water in streams (m)
end

# Array-type params V/M stay loose (un-pinned) so `adapt(MtlArray/CuArray, wd)` can swap the
# storage to a device array (and AD can swap Dual-typed arrays in). The `{FT}` convenience
# ctor defaults them to Vector{FT}/Matrix{FT} at the working precision; @kwdef still gives
# the no-FT `WaterDiagnosticBulkData()` (infers Float64 from the defaults). Same pattern as
# WaterStateBulkData. Adapt.@adapt_structure rebuilds the struct with each field adapted.
WaterDiagnosticBulkData{FT}(; kwargs...) where {FT<:Real} =
    WaterDiagnosticBulkData{FT, Vector{FT}, Matrix{FT}}(; kwargs...)
Adapt.@adapt_structure WaterDiagnosticBulkData

"""
    waterdiagnosticbulk_init!(wd, nc, np, nl, ng)

Allocate and initialize all fields of a `WaterDiagnosticBulkData` instance for
`nc` columns, `np` patches, `nl` landunits, and `ng` gridcells.

Ported from `InitBulkAllocate` + parent `InitAllocate` in Fortran.
"""
function waterdiagnosticbulk_init!(wd::WaterDiagnosticBulkData{FT}, nc::Int, np::Int, nl::Int, ng::Int) where {FT}
    nlevgrnd       = varpar.nlevgrnd
    nlevsno        = varpar.nlevsno
    nlevsoi        = varpar.nlevsoi
    nlevmaxurbgrnd = varpar.nlevmaxurbgrnd
    nlevtot_snow_soil = nlevsno + nlevgrnd  # for combined snow+soil arrays

    # =================================================================
    # Parent fields (waterdiagnostic_type)
    # =================================================================

    # Column 1D (parent)
    wd.snowice_col                 = fill(FT(NaN), nc)
    wd.snowliq_col                 = fill(FT(NaN), nc)
    wd.total_plant_stored_h2o_col  = fill(FT(NaN), nc)
    wd.h2osoi_liqice_10cm_col      = fill(FT(NaN), nc)
    wd.qg_snow_col                 = fill(FT(NaN), nc)
    wd.qg_soil_col                 = fill(FT(NaN), nc)
    wd.qg_h2osfc_col               = fill(FT(NaN), nc)
    wd.qg_col                      = fill(FT(NaN), nc)

    # Patch 1D (parent)
    wd.h2ocan_patch                = fill(FT(NaN), np)
    wd.q_ref2m_patch               = fill(FT(NaN), np)

    # Landunit 1D (parent)
    wd.qaf_lun                     = fill(FT(NaN), nl)

    # Gridcell 1D (parent)
    wd.tws_grc                     = fill(FT(NaN), ng)

    # =================================================================
    # Bulk-specific fields (waterdiagnosticbulk_type)
    # =================================================================

    # Column 1D
    wd.h2osno_total_col            = fill(FT(NaN), nc)
    wd.snow_depth_col              = fill(FT(NaN), nc)
    wd.snow_5day_col               = fill(FT(NaN), nc)
    wd.snowdp_col                  = fill(FT(NaN), nc)
    wd.snomelt_accum_col           = fill(zero(FT), nc)
    wd.h2osoi_liq_tot_col          = fill(FT(NaN), nc)
    wd.h2osoi_ice_tot_col          = fill(FT(NaN), nc)
    wd.exice_subs_tot_col          = fill(zero(FT), nc)   # initialized to 0 in Fortran
    wd.exice_vol_tot_col           = fill(zero(FT), nc)   # initialized to 0 in Fortran
    wd.snw_rds_top_col             = fill(FT(NaN), nc)
    wd.h2osno_top_col              = fill(FT(NaN), nc)
    wd.sno_liq_top_col             = fill(FT(NaN), nc)
    wd.dqgdT_col                   = fill(FT(NaN), nc)

    # Fractions (column 1D)
    wd.frac_sno_col                = fill(FT(NaN), nc)
    wd.frac_sno_eff_col            = fill(FT(NaN), nc)
    wd.frac_h2osfc_col             = fill(FT(NaN), nc)
    wd.frac_h2osfc_nosnow_col      = fill(FT(NaN), nc)
    wd.wf_col                      = fill(FT(NaN), nc)
    wd.wf2_col                     = fill(FT(NaN), nc)

    # Summed fluxes (column 1D)
    wd.qflx_prec_grnd_col          = fill(FT(NaN), nc)

    # Column 2D
    wd.snow_layer_unity_col        = fill(FT(NaN), nc, nlevsno)              # (-nlevsno+1:0)
    wd.bw_col                      = fill(FT(NaN), nc, nlevsno)              # (-nlevsno+1:0)
    wd.air_vol_col                 = fill(FT(NaN), nc, nlevgrnd)             # (1:nlevgrnd)
    wd.h2osoi_liqvol_col           = fill(FT(NaN), nc, nlevtot_snow_soil)    # (-nlevsno+1:nlevgrnd)
    wd.swe_old_col                 = fill(zero(FT), nc, nlevsno)              # (-nlevsno+1:0)
    wd.exice_subs_col              = fill(zero(FT), nc, nlevmaxurbgrnd)       # (1:nlevmaxurbgrnd), init 0
    wd.exice_vol_col               = fill(zero(FT), nc, nlevsoi)              # (1:nlevsoi), init 0
    wd.snw_rds_col                 = fill(FT(NaN), nc, nlevsno)              # (-nlevsno+1:0)
    wd.frac_iceold_col             = fill(FT(NaN), nc, nlevtot_snow_soil)    # (-nlevsno+1:nlevgrnd)

    # Patch 1D
    wd.iwue_ln_patch               = fill(FT(NaN), np)
    wd.vpd_ref2m_patch             = fill(FT(NaN), np)
    wd.rh_ref2m_patch              = fill(FT(NaN), np)
    wd.rh_ref2m_r_patch            = fill(FT(NaN), np)
    wd.rh_ref2m_u_patch            = fill(FT(NaN), np)
    wd.rh_af_patch                 = fill(FT(NaN), np)
    wd.rh10_af_patch               = fill(FT(SPVAL), np)  # initialized to spval in Fortran
    wd.fwet_patch                  = fill(FT(NaN), np)
    wd.fcansno_patch               = fill(FT(NaN), np)
    wd.fdry_patch                  = fill(FT(NaN), np)
    wd.qflx_prec_intr_patch        = fill(FT(NaN), np)

    # Landunit 1D
    wd.stream_water_depth_lun      = fill(FT(NaN), nl)

    return nothing
end

"""
    waterdiagnosticbulk_clean!(wd)

Deallocate (reset to empty) all fields of a `WaterDiagnosticBulkData` instance.
"""
function waterdiagnosticbulk_clean!(wd::WaterDiagnosticBulkData{FT}) where {FT}
    # Parent fields
    wd.snowice_col                 = FT[]
    wd.snowliq_col                 = FT[]
    wd.total_plant_stored_h2o_col  = FT[]
    wd.h2osoi_liqice_10cm_col      = FT[]
    wd.qg_snow_col                 = FT[]
    wd.qg_soil_col                 = FT[]
    wd.qg_h2osfc_col               = FT[]
    wd.qg_col                      = FT[]
    wd.h2ocan_patch                = FT[]
    wd.q_ref2m_patch               = FT[]
    wd.qaf_lun                     = FT[]
    wd.tws_grc                     = FT[]

    # Bulk column 1D
    wd.h2osno_total_col            = FT[]
    wd.snow_depth_col              = FT[]
    wd.snow_5day_col               = FT[]
    wd.snowdp_col                  = FT[]
    wd.snomelt_accum_col           = FT[]
    wd.h2osoi_liq_tot_col          = FT[]
    wd.h2osoi_ice_tot_col          = FT[]
    wd.exice_subs_tot_col          = FT[]
    wd.exice_vol_tot_col           = FT[]
    wd.snw_rds_top_col             = FT[]
    wd.h2osno_top_col              = FT[]
    wd.sno_liq_top_col             = FT[]
    wd.dqgdT_col                   = FT[]
    wd.frac_sno_col                = FT[]
    wd.frac_sno_eff_col            = FT[]
    wd.frac_h2osfc_col             = FT[]
    wd.frac_h2osfc_nosnow_col      = FT[]
    wd.wf_col                      = FT[]
    wd.wf2_col                     = FT[]
    wd.qflx_prec_grnd_col          = FT[]

    # Bulk column 2D
    wd.snow_layer_unity_col        = Matrix{FT}(undef, 0, 0)
    wd.bw_col                      = Matrix{FT}(undef, 0, 0)
    wd.air_vol_col                 = Matrix{FT}(undef, 0, 0)
    wd.h2osoi_liqvol_col           = Matrix{FT}(undef, 0, 0)
    wd.swe_old_col                 = Matrix{FT}(undef, 0, 0)
    wd.exice_subs_col              = Matrix{FT}(undef, 0, 0)
    wd.exice_vol_col               = Matrix{FT}(undef, 0, 0)
    wd.snw_rds_col                 = Matrix{FT}(undef, 0, 0)
    wd.frac_iceold_col             = Matrix{FT}(undef, 0, 0)

    # Bulk patch 1D
    wd.iwue_ln_patch               = FT[]
    wd.vpd_ref2m_patch             = FT[]
    wd.rh_ref2m_patch              = FT[]
    wd.rh_ref2m_r_patch            = FT[]
    wd.rh_ref2m_u_patch            = FT[]
    wd.rh_af_patch                 = FT[]
    wd.rh10_af_patch               = FT[]
    wd.fwet_patch                  = FT[]
    wd.fcansno_patch               = FT[]
    wd.fdry_patch                  = FT[]
    wd.qflx_prec_intr_patch        = FT[]

    # Bulk landunit 1D
    wd.stream_water_depth_lun      = FT[]

    return nothing
end

"""
    waterdiagnosticbulk_init_cold!(wd, bounds_col, bounds_patch;
        snow_depth_input_col, h2osno_input_col, snl_col,
        landunit_col, urbpoi, snw_rds_min, zlnd)

Initialize cold-start conditions for water diagnostic bulk variables.

Ported from `InitBulkCold` in `WaterDiagnosticBulkType.F90`.
"""
function waterdiagnosticbulk_init_cold!(wd::WaterDiagnosticBulkData,
                                         bounds_col::UnitRange{Int},
                                         bounds_patch::UnitRange{Int};
                                         snow_depth_input_col::Vector{<:Real},
                                         h2osno_input_col::Vector{<:Real},
                                         snl_col::Vector{Int},
                                         landunit_col::Vector{Int},
                                         urbpoi::BitVector,
                                         snw_rds_min::Real = waterdiagbulk_params.snw_rds_min,
                                         zlnd::Real = waterdiagbulk_params.zlnd)
    nlevsno = varpar.nlevsno

    # Set snow_depth and snow_layer_unity
    for c in bounds_col
        wd.snow_depth_col[c] = snow_depth_input_col[c]
        for j in 1:nlevsno
            wd.snow_layer_unity_col[c, j] = 1.0
        end
    end

    # wf/wf2 to spval
    for c in bounds_col
        wd.wf_col[c] = SPVAL
        wd.wf2_col[c] = SPVAL
    end

    # frac_h2osfc to zero
    for c in bounds_col
        wd.frac_h2osfc_col[c] = 0.0
    end

    # Patch zeros
    for p in bounds_patch
        wd.fwet_patch[p] = 0.0
        wd.fdry_patch[p] = 0.0
        wd.fcansno_patch[p] = 0.0
        wd.qflx_prec_intr_patch[p] = 0.0
    end

    # Set snow cover fraction
    for c in bounds_col
        l = landunit_col[c]
        if urbpoi[l]
            # From Bonan 1996 (LSM technical note)
            wd.frac_sno_col[c] = min(wd.snow_depth_col[c] / 0.05, 1.0)
        else
            wd.frac_sno_col[c] = 0.0
            # snow cover fraction as in Niu and Yang 2007
            if wd.snow_depth_col[c] > 0.0
                snowbd = min(400.0, h2osno_input_col[c] / wd.snow_depth_col[c])  # bulk density (kg/m3)
                fmelt = (snowbd / 100.0)^1.0
                # 100 is assumed fresh snow density; 1 is melting factor
                wd.frac_sno_col[c] = tanh(wd.snow_depth_col[c] / (2.5 * zlnd * fmelt))
            end
        end
    end

    # Snow grain radius initialization
    for c in bounds_col
        if snl_col[c] < 0
            # Active snow layers: set to snw_rds_min
            for j in (snl_col[c] + 1):0
                jj = j + nlevsno  # convert to 1-based
                wd.snw_rds_col[c, jj] = snw_rds_min
            end
            # Inactive snow layers: zero
            for j in (-nlevsno + 1):snl_col[c]
                jj = j + nlevsno
                wd.snw_rds_col[c, jj] = 0.0
            end
            wd.snw_rds_top_col[c] = snw_rds_min
        elseif h2osno_input_col[c] > 0.0
            # No resolved layers but snow exists: set top layer
            jj_zero = 0 + nlevsno  # layer index for j=0
            wd.snw_rds_col[c, jj_zero] = snw_rds_min
            for j in (-nlevsno + 1):-1
                jj = j + nlevsno
                wd.snw_rds_col[c, jj] = 0.0
            end
            wd.snw_rds_top_col[c] = SPVAL
            wd.sno_liq_top_col[c] = SPVAL
        else
            # No snow at all
            for j in 1:nlevsno
                wd.snw_rds_col[c, j] = 0.0
            end
            wd.snw_rds_top_col[c] = SPVAL
            wd.sno_liq_top_col[c] = SPVAL
        end
    end

    return nothing
end

# ------------------------------------------------------------------
# waterdiagnosticbulk_summary! kernels
# Structs are not isbits, so each kernel takes the INDIVIDUAL wd.<field>
# output arrays + the input arrays (one thread per patch or per column).
# T(...) keeps literals/consts in the device eltype (Metal: Float32-only);
# on CPU T==Float64 so all arithmetic is byte-identical to the original loop.
# Loop-carried sums (`+=` over j) accumulate into thread-locals and write the
# output element ONCE; per-(c,j) writes (exice_vol_col) are own-index direct.
# ------------------------------------------------------------------

# Per-patch: qflx_prec_intr = intercepted liq + intercepted snow
@kernel function _wdb_prec_intr_kernel!(qflx_prec_intr_patch,
        @Const(mask_soilp), @Const(qflx_intercepted_liq_patch),
        @Const(qflx_intercepted_snow_patch), pmin::Int, pmax::Int)
    p = @index(Global)
    @inbounds if pmin <= p <= pmax && mask_soilp[p]
        qflx_prec_intr_patch[p] = qflx_intercepted_liq_patch[p] + qflx_intercepted_snow_patch[p]
    end
end

# Per-column (mask_allc): copy h2osno_total and sum qflx_prec_grnd.
@kernel function _wdb_allc_kernel!(h2osno_total_out, qflx_prec_grnd_out,
        @Const(mask_allc), @Const(h2osno_total_col),
        @Const(qflx_liq_grnd_col), @Const(qflx_snow_grnd_col), cmin::Int, cmax::Int)
    c = @index(Global)
    @inbounds if cmin <= c <= cmax && mask_allc[c]
        h2osno_total_out[c] = h2osno_total_col[c]
        qflx_prec_grnd_out[c] = qflx_liq_grnd_col[c] + qflx_snow_grnd_col[c]
    end
end

# Per-column (mask_nolakec, non-urban): init totals, sum over soil layers
# (top-10cm + liq/ice totals + excess-ice subsidence), then excess-ice volume.
@kernel function _wdb_nolakec_kernel!(h2osoi_liqice_10cm_col, h2osoi_liq_tot_col,
        h2osoi_ice_tot_col, exice_subs_tot_col, exice_vol_tot_col, exice_vol_col,
        @Const(mask_nolakec), @Const(landunit_col), @Const(urbpoi), @Const(lun_itype),
        @Const(h2osoi_liq_col), @Const(h2osoi_ice_col), @Const(zi_col), @Const(dz_col),
        @Const(exice_subs_col), @Const(excess_ice_col),
        nlevsoi::Int, nlevsno::Int, cmin::Int, cmax::Int)
    c = @index(Global)
    @inbounds if cmin <= c <= cmax && mask_nolakec[c]
        T = eltype(h2osoi_liq_tot_col)
        l = landunit_col[c]
        if !urbpoi[l]
            # Sum over soil layers (combined snow+soil index jj = j + nlevsno).
            acc_10cm = zero(T)
            acc_liq  = zero(T)
            acc_ice  = zero(T)
            acc_subs = zero(T)
            is_soilcrop = (lun_itype[l] == ISTSOIL || lun_itype[l] == ISTCROP)
            for j in 1:nlevsoi
                jj = j + nlevsno
                if zi_col[c, j] <= T(0.1)
                    fracl = one(T)
                    acc_10cm += (h2osoi_liq_col[c, jj] + h2osoi_ice_col[c, jj]) * fracl
                else
                    if zi_col[c, j] > T(0.1) && zi_col[c, j-1] < T(0.1)
                        fracl = (T(0.1) - zi_col[c, j-1]) / dz_col[c, j]
                        acc_10cm += (h2osoi_liq_col[c, jj] + h2osoi_ice_col[c, jj]) * fracl
                    end
                end
                acc_liq += h2osoi_liq_col[c, jj]
                acc_ice += h2osoi_ice_col[c, jj]
                if is_soilcrop
                    acc_subs += exice_subs_col[c, j]
                end
            end
            h2osoi_liqice_10cm_col[c] = acc_10cm
            h2osoi_liq_tot_col[c] = acc_liq
            h2osoi_ice_tot_col[c] = acc_ice
            exice_subs_tot_col[c] = acc_subs

            # Excess ice volume calculations.
            if is_soilcrop
                dz_tot = zero(T)
                acc_voltot = zero(T)
                for j in 1:nlevsoi
                    dz_ext = dz_col[c, j] + excess_ice_col[c, j] / T(DENICE)
                    exice_vol_col[c, j] = excess_ice_col[c, j] / (T(DENICE) * dz_ext)
                    dz_tot += dz_ext
                    acc_voltot += exice_vol_col[c, j] * dz_ext  # (m)
                end
                exice_vol_tot_col[c] = acc_voltot / dz_tot  # (m3/m3)
            else
                exice_vol_tot_col[c] = zero(T)
            end
        end
    end
end

"""
    waterdiagnosticbulk_summary!(wd, bounds_col;
        mask_soilp, mask_allc, mask_nolakec,
        waterstate, waterflux,
        dz_col, zi_col, landunit_col, urbpoi, lun_itype)

Compute end-of-timestep summaries of water diagnostic terms.

Ported from `Summary` in `WaterDiagnosticBulkType.F90`.
"""
function waterdiagnosticbulk_summary!(wd::WaterDiagnosticBulkData,
                                       bounds_col::UnitRange{Int},
                                       bounds_patch::UnitRange{Int};
                                       mask_soilp::AbstractVector{Bool},
                                       mask_allc::AbstractVector{Bool},
                                       mask_nolakec::AbstractVector{Bool},
                                       h2osoi_ice_col::AbstractMatrix{<:Real},
                                       h2osoi_liq_col::AbstractMatrix{<:Real},
                                       excess_ice_col::AbstractMatrix{<:Real},
                                       qflx_intercepted_liq_patch::AbstractVector{<:Real},
                                       qflx_intercepted_snow_patch::AbstractVector{<:Real},
                                       qflx_liq_grnd_col::AbstractVector{<:Real},
                                       qflx_snow_grnd_col::AbstractVector{<:Real},
                                       h2osno_total_col::AbstractVector{<:Real},
                                       dz_col::AbstractMatrix{<:Real},
                                       zi_col::AbstractMatrix{<:Real},
                                       landunit_col::AbstractVector{<:Integer},
                                       urbpoi::AbstractVector{Bool},
                                       lun_itype::AbstractVector{<:Integer})
    nlevsoi = varpar.nlevsoi
    nlevsno = varpar.nlevsno

    # qflx_prec_intr = intercepted liq + intercepted snow (per patch)
    if !isempty(bounds_patch)
        _launch!(_wdb_prec_intr_kernel!, wd.qflx_prec_intr_patch,
            mask_soilp, qflx_intercepted_liq_patch, qflx_intercepted_snow_patch,
            first(bounds_patch), last(bounds_patch);
            ndrange = length(wd.qflx_prec_intr_patch))
    end

    # Copy h2osno_total and sum qflx_prec_grnd (per column, mask_allc).
    if !isempty(bounds_col)
        _launch!(_wdb_allc_kernel!, wd.h2osno_total_col, wd.qflx_prec_grnd_col,
            mask_allc, h2osno_total_col, qflx_liq_grnd_col, qflx_snow_grnd_col,
            first(bounds_col), last(bounds_col);
            ndrange = length(wd.h2osno_total_col))
    end

    # Soil-layer sums + excess-ice volume (per column, mask_nolakec, non-urban).
    # Note: h2osoi_liq/ice_col use combined snow+soil indexing;
    # soil layer j maps to column index j + nlevsno in Julia arrays.
    if !isempty(bounds_col)
        _launch!(_wdb_nolakec_kernel!, wd.h2osoi_liqice_10cm_col, wd.h2osoi_liq_tot_col,
            wd.h2osoi_ice_tot_col, wd.exice_subs_tot_col, wd.exice_vol_tot_col, wd.exice_vol_col,
            mask_nolakec, landunit_col, urbpoi, lun_itype,
            h2osoi_liq_col, h2osoi_ice_col, zi_col, dz_col, wd.exice_subs_col, excess_ice_col,
            nlevsoi, nlevsno, first(bounds_col), last(bounds_col);
            ndrange = length(wd.h2osoi_liq_tot_col))
    end

    return nothing
end

"""
    waterdiagnosticbulk_reset_filter!(wd, mask, bounds_col)

Initialize SNICAR variables for fresh snow columns matching mask.

Ported from `ResetBulkFilter` in `WaterDiagnosticBulkType.F90`.
"""
function waterdiagnosticbulk_reset_filter!(wd::WaterDiagnosticBulkData,
                                            mask::BitVector,
                                            bounds_col::UnitRange{Int};
                                            snw_rds_min::Real = waterdiagbulk_params.snw_rds_min)
    for c in bounds_col
        mask[c] || continue
        waterdiagnosticbulk_reset!(wd, c; snw_rds_min=snw_rds_min)
    end
    return nothing
end

"""
    waterdiagnosticbulk_reset!(wd, c)

Initialize SNICAR variables for a fresh snow column.

Ported from `ResetBulk` in `WaterDiagnosticBulkType.F90`.
"""
function waterdiagnosticbulk_reset!(wd::WaterDiagnosticBulkData, c::Int;
                                     snw_rds_min::Real = waterdiagbulk_params.snw_rds_min)
    nlevsno = varpar.nlevsno
    # j=0 maps to index nlevsno in our 1-based arrays
    wd.snw_rds_col[c, nlevsno] = snw_rds_min
    return nothing
end

# ==========================================================================
# Stubs for infrastructure-dependent subroutines
# ==========================================================================

"""
    waterdiagnosticbulk_init_history!(wd, bounds_col)

Register water diagnostic bulk fields for history file output.

Ported from `InitBulkHistory` in `WaterDiagnosticBulkType.F90`.
Not implemented (no-op stub). History I/O IS ported
(`src/infrastructure/history_io.jl`); fields are registered in a central
registry rather than per-type methods.
"""
function waterdiagnosticbulk_init_history!(wd::WaterDiagnosticBulkData,
                                            bounds_col::UnitRange{Int})
    return nothing
end

"""
    waterdiagnosticbulk_restart!(wd, bounds_col; flag="read")

Read/write water diagnostic bulk from/to restart file.

Ported from `RestartBulk` in `WaterDiagnosticBulkType.F90`.
Not implemented (no-op stub). Restart I/O IS ported
(`src/infrastructure/restart_io.jl`, `fortran_restart.jl`); restart variables
are declared in a central registry rather than per-type methods.
"""
function waterdiagnosticbulk_restart!(wd::WaterDiagnosticBulkData,
                                       bounds_col::UnitRange{Int};
                                       flag::String = "read")
    return nothing
end

"""
    waterdiagnosticbulk_init_acc_buffer!(wd, bounds_col)

Initialize accumulation buffer for snow 5-day average.

Ported from `InitAccBuffer` in `WaterDiagnosticBulkType.F90`.
Not implemented (no-op stub). The accumulator IS ported
(`src/infrastructure/accumul.jl`, `AccumManager`) and is used by the live
driver (crop GDD, `t_mo_min`); these particular fields are simply not
registered with it.
"""
function waterdiagnosticbulk_init_acc_buffer!(wd::WaterDiagnosticBulkData,
                                               bounds_col::UnitRange{Int})
    return nothing
end

"""
    waterdiagnosticbulk_init_acc_vars!(wd, bounds_col)

Initialize accumulation variables from restart.

Ported from `InitAccVars` in `WaterDiagnosticBulkType.F90`.
Not implemented (no-op stub). The accumulator IS ported
(`src/infrastructure/accumul.jl`, `AccumManager`) and is used by the live
driver (crop GDD, `t_mo_min`); these particular fields are simply not
registered with it.
"""
function waterdiagnosticbulk_init_acc_vars!(wd::WaterDiagnosticBulkData,
                                             bounds_col::UnitRange{Int})
    return nothing
end

"""
    waterdiagnosticbulk_update_acc_vars!(wd, bounds_col)

Update accumulation variables (snow 5-day average).

Ported from `UpdateAccVars` in `WaterDiagnosticBulkType.F90`.
Not implemented (no-op stub). The accumulator IS ported
(`src/infrastructure/accumul.jl`, `AccumManager`) and is used by the live
driver (crop GDD, `t_mo_min`); these particular fields are simply not
registered with it.
"""
function waterdiagnosticbulk_update_acc_vars!(wd::WaterDiagnosticBulkData,
                                               bounds_col::UnitRange{Int})
    return nothing
end
