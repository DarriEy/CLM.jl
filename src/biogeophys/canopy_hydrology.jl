# ==========================================================================
# Ported from: src/biogeophys/CanopyHydrologyMod.F90
# Calculation of:
# (1) water storage of intercepted precipitation
# (2) direct throughfall and canopy drainage of precipitation
# (3) the fraction of foliage covered by water and the fraction
#     of foliage that is dry and transpiring.
# (4) snow layer initialization if the snow accumulation exceeds 10 mm.
# ==========================================================================

# --- Parameters type (from params_type in Fortran) ---

Base.@kwdef mutable struct CanopyHydrologyParamsData
    liq_canopy_storage_scalar  ::Float64 = 0.04   # Canopy-storage-of-liquid-water parameter (kg/m2)
    snow_canopy_storage_scalar ::Float64 = 6.0    # Canopy-storage-of-snow parameter (kg/m2)
    snowcan_unload_temp_fact   ::Float64 = 1.87e5  # Temperature canopy snow unload scaling (C2 in Eq. 14, Roesch et al. 2001) (K*s)
    snowcan_unload_wind_fact   ::Float64 = 1.56e5  # Wind canopy snow unload scaling (modifies 1.56e5, C3 in Eq. 15, Roesch et al. 2001) (-)
    interception_fraction      ::Float64 = 1.0    # Fraction of intercepted precipitation (-)
    maximum_leaf_wetted_fraction::Float64 = 1.0   # Maximum fraction of leaf that may be wet (-)
end

const canopy_hydrology_params = CanopyHydrologyParamsData()

# --- Module control state ---

Base.@kwdef mutable struct CanopyHydrologyControl
    use_clm5_fpi::Bool = false   # use clm5 fpi equation
end

const canopy_hydrology_ctrl = CanopyHydrologyControl()

# --- Namelist reader ---

function canopy_hydrology_read_nml!(; use_clm5_fpi::Bool = false)
    canopy_hydrology_ctrl.use_clm5_fpi = use_clm5_fpi
    return nothing
end

# --- Parameter reader ---

function canopy_hydrology_read_params!(;
        liq_canopy_storage_scalar::Real = 0.04,
        snow_canopy_storage_scalar::Real = 6.0,
        snowcan_unload_temp_fact::Real = 1.87e5,
        snowcan_unload_wind_fact::Real = 1.56e5,
        interception_fraction::Real = 1.0,
        maximum_leaf_wetted_fraction::Real = 1.0)
    canopy_hydrology_params.liq_canopy_storage_scalar   = liq_canopy_storage_scalar
    canopy_hydrology_params.snow_canopy_storage_scalar  = snow_canopy_storage_scalar
    canopy_hydrology_params.snowcan_unload_temp_fact    = snowcan_unload_temp_fact
    canopy_hydrology_params.snowcan_unload_wind_fact    = snowcan_unload_wind_fact
    canopy_hydrology_params.interception_fraction       = interception_fraction
    canopy_hydrology_params.maximum_leaf_wetted_fraction = maximum_leaf_wetted_fraction
    return nothing
end

# -----------------------------------------------------------------------
# SumFlux_TopOfCanopyInputs
# Compute patch-level precipitation inputs for bulk water or one tracer
# Ported from SumFlux_TopOfCanopyInputs in CanopyHydrologyMod.F90
# -----------------------------------------------------------------------

@kernel function _canhyd_top_of_canopy_kernel!(qflx_liq_above_canopy, forc_snow_patch,
        @Const(mask_nolakep), @Const(column), @Const(forc_rain),
        @Const(forc_snow_col), @Const(qflx_irrig_sprinkler), @Const(pstart))
    i = @index(Global)
    @inbounds begin
        p = pstart + i - 1
        if mask_nolakep[p]
            c = column[p]
            qflx_liq_above_canopy[p] = forc_rain[c] + qflx_irrig_sprinkler[p]
            forc_snow_patch[p] = forc_snow_col[c]
        end
    end
end

canhyd_top_of_canopy!(qflx_liq_above_canopy, forc_snow_patch, mask_nolakep, column,
        forc_rain, forc_snow_col, qflx_irrig_sprinkler, pstart, np) =
    _launch!(_canhyd_top_of_canopy_kernel!, qflx_liq_above_canopy, forc_snow_patch,
        mask_nolakep, column, forc_rain, forc_snow_col, qflx_irrig_sprinkler, pstart;
        ndrange=np)

function sum_flux_top_of_canopy_inputs!(
        patch::PatchData,
        mask_nolakep::AbstractVector{Bool},
        bounds_p::UnitRange{Int},
        # Inputs (column-level)
        forc_rain::AbstractVector{<:Real},
        forc_snow_col::AbstractVector{<:Real},
        # Inputs (patch-level)
        qflx_irrig_sprinkler::AbstractVector{<:Real},
        # Outputs (patch-level)
        qflx_liq_above_canopy::AbstractVector{<:Real},
        forc_snow_patch::AbstractVector{<:Real})

    canhyd_top_of_canopy!(qflx_liq_above_canopy, forc_snow_patch, mask_nolakep,
        patch.column, forc_rain, forc_snow_col, qflx_irrig_sprinkler,
        first(bounds_p), length(bounds_p))

    return nothing
end

# -----------------------------------------------------------------------
# BulkFlux_CanopyInterceptionAndThroughfall
# Compute canopy interception and throughfall for bulk water
# Ported from BulkFlux_CanopyInterceptionAndThroughfall in CanopyHydrologyMod.F90
# -----------------------------------------------------------------------

@kernel function _canhyd_intercept_throughfall_kernel!(qflx_through_snow, qflx_through_liq,
        qflx_intercepted_snow, qflx_intercepted_liq, check_point,
        @Const(mask_nolakep), @Const(column), @Const(itype), @Const(frac_veg_nosno),
        @Const(elai), @Const(esai), @Const(forc_snow), @Const(qflx_liq_above_canopy),
        @Const(use_clm5_fpi), @Const(interception_fraction),
        @Const(icol_sunwall), @Const(icol_shadewall), @Const(pstart))
    i = @index(Global)
    @inbounds begin
        p = pstart + i - 1
        if mask_nolakep[p]
            c = column[p]
            T = eltype(forc_snow)
            lai = elai[p] + esai[p]

            chk = (frac_veg_nosno[p] == 1 && (forc_snow[p] + qflx_liq_above_canopy[p]) > zero(T))
            check_point[p] = chk

            if chk
                # Coefficient of interception
                if use_clm5_fpi
                    fpiliq = oftype(lai, interception_fraction) * tanh(lai)
                else
                    fpiliq = oftype(lai, 0.25) * (one(lai) - exp(oftype(lai, -0.5) * lai))
                end

                fpisnow = (one(lai) - exp(oftype(lai, -0.5) * lai))  # max interception of 1

                # Direct throughfall
                qflx_through_snow[p] = forc_snow[p] * (one(lai) - fpisnow)
                qflx_through_liq[p]  = qflx_liq_above_canopy[p] * (one(lai) - fpiliq)

                # Canopy interception
                qflx_intercepted_snow[p] = forc_snow[p] * fpisnow
                qflx_intercepted_liq[p]  = qflx_liq_above_canopy[p] * fpiliq
            else
                # Note that special landunits will be handled here, in addition to soil points
                # with frac_veg_nosno == 0.
                qflx_intercepted_snow[p] = zero(T)
                qflx_intercepted_liq[p]  = zero(T)
                if itype[c] == icol_sunwall || itype[c] == icol_shadewall
                    qflx_through_snow[p] = zero(T)
                    qflx_through_liq[p]  = zero(T)
                else
                    qflx_through_snow[p] = forc_snow[p]
                    qflx_through_liq[p]  = qflx_liq_above_canopy[p]
                end
            end
        end
    end
end

canhyd_intercept_throughfall!(qflx_through_snow, qflx_through_liq, qflx_intercepted_snow,
        qflx_intercepted_liq, check_point, mask_nolakep, column, itype, frac_veg_nosno,
        elai, esai, forc_snow, qflx_liq_above_canopy, use_clm5_fpi, interception_fraction,
        icol_sunwall, icol_shadewall, pstart, np) =
    _launch!(_canhyd_intercept_throughfall_kernel!, qflx_through_snow, qflx_through_liq,
        qflx_intercepted_snow, qflx_intercepted_liq, check_point, mask_nolakep, column,
        itype, frac_veg_nosno, elai, esai, forc_snow, qflx_liq_above_canopy, use_clm5_fpi,
        convert(eltype(forc_snow), interception_fraction),
        icol_sunwall, icol_shadewall, pstart; ndrange=np)

function bulk_flux_canopy_interception_and_throughfall!(
        patch::PatchData,
        col::ColumnData,
        mask_nolakep::AbstractVector{Bool},
        bounds_p::UnitRange{Int},
        # Inputs
        frac_veg_nosno::AbstractVector{<:Integer},
        elai::AbstractVector{<:Real},
        esai::AbstractVector{<:Real},
        forc_snow::AbstractVector{<:Real},
        qflx_liq_above_canopy::AbstractVector{<:Real},
        # Outputs
        qflx_through_snow::AbstractVector{<:Real},
        qflx_through_liq::AbstractVector{<:Real},
        qflx_intercepted_snow::AbstractVector{<:Real},
        qflx_intercepted_liq::AbstractVector{<:Real},
        check_point_for_interception_and_excess::AbstractVector{Bool})

    canhyd_intercept_throughfall!(qflx_through_snow, qflx_through_liq, qflx_intercepted_snow,
        qflx_intercepted_liq, check_point_for_interception_and_excess,
        mask_nolakep, patch.column, col.itype, frac_veg_nosno, elai, esai, forc_snow,
        qflx_liq_above_canopy, canopy_hydrology_ctrl.use_clm5_fpi,
        canopy_hydrology_params.interception_fraction, ICOL_SUNWALL, ICOL_SHADEWALL,
        first(bounds_p), length(bounds_p))

    return nothing
end

# -----------------------------------------------------------------------
# TracerFlux_CanopyInterceptionAndThroughfall
# Calculate canopy interception and throughfall for one tracer
# Ported from TracerFlux_CanopyInterceptionAndThroughfall in CanopyHydrologyMod.F90
#
# In the Fortran code this calls CalcTracerFromBulk. We inline the
# proportional scaling: tracer_val = (tracer_source / bulk_source) * bulk_val,
# with tracer_val = 0 when bulk_source == 0.
# -----------------------------------------------------------------------

function tracer_flux_canopy_interception_and_throughfall!(
        mask_nolakep::AbstractVector{Bool},
        bounds_p::UnitRange{Int},
        # Bulk inputs
        bulk_forc_snow::AbstractVector{<:Real},
        bulk_qflx_liq_above_canopy::AbstractVector{<:Real},
        bulk_qflx_through_snow::AbstractVector{<:Real},
        bulk_qflx_intercepted_snow::AbstractVector{<:Real},
        bulk_qflx_through_liq::AbstractVector{<:Real},
        bulk_qflx_intercepted_liq::AbstractVector{<:Real},
        # Tracer inputs
        trac_forc_snow::AbstractVector{<:Real},
        trac_qflx_liq_above_canopy::AbstractVector{<:Real},
        # Tracer outputs
        trac_qflx_through_snow::AbstractVector{<:Real},
        trac_qflx_intercepted_snow::AbstractVector{<:Real},
        trac_qflx_through_liq::AbstractVector{<:Real},
        trac_qflx_intercepted_liq::AbstractVector{<:Real})

    canhyd_tracer_intercept_throughfall!(
        trac_qflx_through_snow, trac_qflx_intercepted_snow,
        trac_qflx_through_liq, trac_qflx_intercepted_liq,
        mask_nolakep, bulk_forc_snow, bulk_qflx_liq_above_canopy,
        bulk_qflx_through_snow, bulk_qflx_intercepted_snow,
        bulk_qflx_through_liq, bulk_qflx_intercepted_liq,
        trac_forc_snow, trac_qflx_liq_above_canopy,
        first(bounds_p), length(bounds_p))

    return nothing
end

@kernel function _canhyd_tracer_intercept_throughfall_kernel!(
        trac_qflx_through_snow, trac_qflx_intercepted_snow,
        trac_qflx_through_liq, trac_qflx_intercepted_liq,
        @Const(mask_nolakep), @Const(bulk_forc_snow), @Const(bulk_qflx_liq_above_canopy),
        @Const(bulk_qflx_through_snow), @Const(bulk_qflx_intercepted_snow),
        @Const(bulk_qflx_through_liq), @Const(bulk_qflx_intercepted_liq),
        @Const(trac_forc_snow), @Const(trac_qflx_liq_above_canopy), @Const(pstart))
    i = @index(Global)
    @inbounds begin
        p = pstart + i - 1
        if mask_nolakep[p]
            # Snow throughfall tracer
            if bulk_forc_snow[p] != 0.0
                ratio = trac_forc_snow[p] / bulk_forc_snow[p]
                trac_qflx_through_snow[p] = ratio * bulk_qflx_through_snow[p]
                trac_qflx_intercepted_snow[p] = ratio * bulk_qflx_intercepted_snow[p]
            else
                trac_qflx_through_snow[p] = 0.0
                trac_qflx_intercepted_snow[p] = 0.0
            end

            # Liquid throughfall tracer
            if bulk_qflx_liq_above_canopy[p] != 0.0
                ratio = trac_qflx_liq_above_canopy[p] / bulk_qflx_liq_above_canopy[p]
                trac_qflx_through_liq[p] = ratio * bulk_qflx_through_liq[p]
                trac_qflx_intercepted_liq[p] = ratio * bulk_qflx_intercepted_liq[p]
            else
                trac_qflx_through_liq[p] = 0.0
                trac_qflx_intercepted_liq[p] = 0.0
            end
        end
    end
end

canhyd_tracer_intercept_throughfall!(
        trac_qflx_through_snow, trac_qflx_intercepted_snow,
        trac_qflx_through_liq, trac_qflx_intercepted_liq,
        mask_nolakep, bulk_forc_snow, bulk_qflx_liq_above_canopy,
        bulk_qflx_through_snow, bulk_qflx_intercepted_snow,
        bulk_qflx_through_liq, bulk_qflx_intercepted_liq,
        trac_forc_snow, trac_qflx_liq_above_canopy, pstart, np) =
    _launch!(_canhyd_tracer_intercept_throughfall_kernel!,
        trac_qflx_through_snow, trac_qflx_intercepted_snow,
        trac_qflx_through_liq, trac_qflx_intercepted_liq,
        mask_nolakep, bulk_forc_snow, bulk_qflx_liq_above_canopy,
        bulk_qflx_through_snow, bulk_qflx_intercepted_snow,
        bulk_qflx_through_liq, bulk_qflx_intercepted_liq,
        trac_forc_snow, trac_qflx_liq_above_canopy, pstart; ndrange=np)

# -----------------------------------------------------------------------
# UpdateState_AddInterceptionToCanopy
# Update snocan and liqcan based on interception, for bulk or one tracer
# Ported from UpdateState_AddInterceptionToCanopy in CanopyHydrologyMod.F90
# -----------------------------------------------------------------------

function update_state_add_interception_to_canopy!(
        mask_soilp::AbstractVector{Bool},
        bounds_p::UnitRange{Int},
        dtime::Real,
        # Inputs
        qflx_intercepted_snow::AbstractVector{<:Real},
        qflx_intercepted_liq::AbstractVector{<:Real},
        # In/Out
        snocan::AbstractVector{<:Real},
        liqcan::AbstractVector{<:Real})

    canhyd_add_interception!(snocan, liqcan, mask_soilp, dtime,
        qflx_intercepted_snow, qflx_intercepted_liq,
        first(bounds_p), length(bounds_p))

    return nothing
end

@kernel function _canhyd_add_interception_kernel!(snocan, liqcan, @Const(mask_soilp),
        @Const(dtime), @Const(qflx_intercepted_snow), @Const(qflx_intercepted_liq),
        @Const(pstart))
    i = @index(Global)
    @inbounds begin
        p = pstart + i - 1
        if mask_soilp[p]
            dt = oftype(snocan[p], dtime)
            snocan[p] = max(zero(snocan[p]), snocan[p] + dt * qflx_intercepted_snow[p])
            liqcan[p] = max(zero(liqcan[p]), liqcan[p] + dt * qflx_intercepted_liq[p])
        end
    end
end

canhyd_add_interception!(snocan, liqcan, mask_soilp, dtime,
        qflx_intercepted_snow, qflx_intercepted_liq, pstart, np) =
    _launch!(_canhyd_add_interception_kernel!, snocan, liqcan, mask_soilp,
        convert(eltype(snocan), dtime),
        qflx_intercepted_snow, qflx_intercepted_liq, pstart; ndrange=np)

# -----------------------------------------------------------------------
# BulkFlux_CanopyExcess
# Compute runoff from canopy due to exceeding maximum storage, for bulk
# Ported from BulkFlux_CanopyExcess in CanopyHydrologyMod.F90
# -----------------------------------------------------------------------

function bulk_flux_canopy_excess!(
        mask_soilp::AbstractVector{Bool},
        bounds_p::UnitRange{Int},
        dtime::Real,
        # Inputs
        elai::AbstractVector{<:Real},
        esai::AbstractVector{<:Real},
        snocan::AbstractVector{<:Real},
        liqcan::AbstractVector{<:Real},
        check_point_for_interception_and_excess::AbstractVector{Bool},
        # Outputs
        qflx_snocanfall::AbstractVector{<:Real},
        qflx_liqcanfall::AbstractVector{<:Real})

    canhyd_canopy_excess!(qflx_snocanfall, qflx_liqcanfall, mask_soilp, dtime,
        elai, esai, snocan, liqcan, check_point_for_interception_and_excess,
        canopy_hydrology_params.liq_canopy_storage_scalar,
        canopy_hydrology_params.snow_canopy_storage_scalar,
        first(bounds_p), length(bounds_p))

    return nothing
end

@kernel function _canhyd_canopy_excess_kernel!(qflx_snocanfall, qflx_liqcanfall,
        @Const(mask_soilp), @Const(dtime), @Const(elai), @Const(esai),
        @Const(snocan), @Const(liqcan), @Const(check_point),
        @Const(liq_canopy_storage_scalar), @Const(snow_canopy_storage_scalar),
        @Const(pstart))
    i = @index(Global)
    @inbounds begin
        p = pstart + i - 1
        if mask_soilp[p]
            T = eltype(qflx_liqcanfall)
            qflx_liqcanfall[p] = zero(T)
            qflx_snocanfall[p] = zero(T)

            if check_point[p]
                lai = elai[p] + esai[p]
                dt = oftype(lai, dtime)
                liqcanmx = oftype(lai, liq_canopy_storage_scalar) * lai
                qflx_liqcanfall[p] = max((liqcan[p] - liqcanmx) / dt, zero(lai))
                snocanmx = oftype(lai, snow_canopy_storage_scalar) * lai
                qflx_snocanfall[p] = max((snocan[p] - snocanmx) / dt, zero(lai))
            end
        end
    end
end

canhyd_canopy_excess!(qflx_snocanfall, qflx_liqcanfall, mask_soilp, dtime,
        elai, esai, snocan, liqcan, check_point,
        liq_canopy_storage_scalar, snow_canopy_storage_scalar, pstart, np) =
    _launch!(_canhyd_canopy_excess_kernel!, qflx_snocanfall, qflx_liqcanfall, mask_soilp,
        convert(eltype(qflx_snocanfall), dtime), elai, esai, snocan, liqcan, check_point,
        convert(eltype(qflx_snocanfall), liq_canopy_storage_scalar),
        convert(eltype(qflx_snocanfall), snow_canopy_storage_scalar), pstart; ndrange=np)

# -----------------------------------------------------------------------
# TracerFlux_CanopyExcess
# Calculate runoff from canopy due to exceeding maximum storage, for one tracer
# Ported from TracerFlux_CanopyExcess in CanopyHydrologyMod.F90
# -----------------------------------------------------------------------

function tracer_flux_canopy_excess!(
        mask_soilp::AbstractVector{Bool},
        bounds_p::UnitRange{Int},
        # Bulk inputs
        bulk_liqcan::AbstractVector{<:Real},
        bulk_snocan::AbstractVector{<:Real},
        bulk_qflx_liqcanfall::AbstractVector{<:Real},
        bulk_qflx_snocanfall::AbstractVector{<:Real},
        # Tracer inputs
        trac_liqcan::AbstractVector{<:Real},
        trac_snocan::AbstractVector{<:Real},
        # Tracer outputs
        trac_qflx_liqcanfall::AbstractVector{<:Real},
        trac_qflx_snocanfall::AbstractVector{<:Real})

    canhyd_tracer_canopy_excess!(trac_qflx_liqcanfall, trac_qflx_snocanfall,
        mask_soilp, bulk_liqcan, bulk_snocan, bulk_qflx_liqcanfall, bulk_qflx_snocanfall,
        trac_liqcan, trac_snocan, first(bounds_p), length(bounds_p))

    return nothing
end

@kernel function _canhyd_tracer_canopy_excess_kernel!(trac_qflx_liqcanfall,
        trac_qflx_snocanfall, @Const(mask_soilp), @Const(bulk_liqcan), @Const(bulk_snocan),
        @Const(bulk_qflx_liqcanfall), @Const(bulk_qflx_snocanfall), @Const(trac_liqcan),
        @Const(trac_snocan), @Const(pstart))
    i = @index(Global)
    @inbounds begin
        p = pstart + i - 1
        if mask_soilp[p]
            # Liquid canopy fall: tracer proportional to bulk
            if bulk_liqcan[p] != 0.0
                trac_qflx_liqcanfall[p] = (trac_liqcan[p] / bulk_liqcan[p]) * bulk_qflx_liqcanfall[p]
            else
                trac_qflx_liqcanfall[p] = 0.0
            end

            # Snow canopy fall: tracer proportional to bulk
            if bulk_snocan[p] != 0.0
                trac_qflx_snocanfall[p] = (trac_snocan[p] / bulk_snocan[p]) * bulk_qflx_snocanfall[p]
            else
                trac_qflx_snocanfall[p] = 0.0
            end
        end
    end
end

canhyd_tracer_canopy_excess!(trac_qflx_liqcanfall, trac_qflx_snocanfall,
        mask_soilp, bulk_liqcan, bulk_snocan, bulk_qflx_liqcanfall, bulk_qflx_snocanfall,
        trac_liqcan, trac_snocan, pstart, np) =
    _launch!(_canhyd_tracer_canopy_excess_kernel!, trac_qflx_liqcanfall, trac_qflx_snocanfall,
        mask_soilp, bulk_liqcan, bulk_snocan, bulk_qflx_liqcanfall, bulk_qflx_snocanfall,
        trac_liqcan, trac_snocan, pstart; ndrange=np)

# -----------------------------------------------------------------------
# UpdateState_RemoveCanfallFromCanopy
# Update snocan and liqcan based on canfall, for bulk or one tracer
# Ported from UpdateState_RemoveCanfallFromCanopy in CanopyHydrologyMod.F90
# -----------------------------------------------------------------------

function update_state_remove_canfall_from_canopy!(
        mask_soilp::AbstractVector{Bool},
        bounds_p::UnitRange{Int},
        dtime::Real,
        # Inputs
        qflx_liqcanfall::AbstractVector{<:Real},
        qflx_snocanfall::AbstractVector{<:Real},
        # In/Out
        liqcan::AbstractVector{<:Real},
        snocan::AbstractVector{<:Real})

    canhyd_remove_canfall!(liqcan, snocan, mask_soilp, dtime,
        qflx_liqcanfall, qflx_snocanfall, first(bounds_p), length(bounds_p))

    return nothing
end

@kernel function _canhyd_remove_canfall_kernel!(liqcan, snocan, @Const(mask_soilp),
        @Const(dtime), @Const(qflx_liqcanfall), @Const(qflx_snocanfall), @Const(pstart))
    i = @index(Global)
    @inbounds begin
        p = pstart + i - 1
        if mask_soilp[p]
            dt = oftype(liqcan[p], dtime)
            liqcan[p] = liqcan[p] - dt * qflx_liqcanfall[p]
            snocan[p] = snocan[p] - dt * qflx_snocanfall[p]
        end
    end
end

canhyd_remove_canfall!(liqcan, snocan, mask_soilp, dtime,
        qflx_liqcanfall, qflx_snocanfall, pstart, np) =
    _launch!(_canhyd_remove_canfall_kernel!, liqcan, snocan, mask_soilp,
        convert(eltype(liqcan), dtime),
        qflx_liqcanfall, qflx_snocanfall, pstart; ndrange=np)

# -----------------------------------------------------------------------
# BulkFlux_SnowUnloading
# Compute snow unloading for bulk
# Ported from BulkFlux_SnowUnloading in CanopyHydrologyMod.F90
# -----------------------------------------------------------------------

function bulk_flux_snow_unloading!(
        patch::PatchData,
        mask_soilp::AbstractVector{Bool},
        bounds_p::UnitRange{Int},
        dtime::Real,
        # Inputs
        frac_veg_nosno::AbstractVector{<:Integer},
        forc_t::AbstractVector{<:Real},        # column-level
        forc_wind::AbstractVector{<:Real},     # gridcell-level
        snocan::AbstractVector{<:Real},
        # Outputs
        qflx_snotempunload::AbstractVector{<:Real},
        qflx_snowindunload::AbstractVector{<:Real},
        qflx_snow_unload::AbstractVector{<:Real})

    canhyd_snow_unloading!(qflx_snotempunload, qflx_snowindunload, qflx_snow_unload,
        mask_soilp, patch.column, patch.gridcell, dtime, frac_veg_nosno,
        forc_t, forc_wind, snocan,
        canopy_hydrology_params.snowcan_unload_temp_fact,
        canopy_hydrology_params.snowcan_unload_wind_fact,
        first(bounds_p), length(bounds_p))

    return nothing
end

@kernel function _canhyd_snow_unloading_kernel!(qflx_snotempunload, qflx_snowindunload,
        qflx_snow_unload, @Const(mask_soilp), @Const(column), @Const(gridcell),
        @Const(dtime), @Const(frac_veg_nosno), @Const(forc_t), @Const(forc_wind),
        @Const(snocan), @Const(snowcan_unload_temp_fact), @Const(snowcan_unload_wind_fact),
        @Const(pstart))
    i = @index(Global)
    @inbounds begin
        p = pstart + i - 1
        if mask_soilp[p]
            T = eltype(snocan)
            if frac_veg_nosno[p] == 1 && snocan[p] > zero(T)
                c = column[p]
                g = gridcell[p]
                dt = oftype(snocan[p], dtime)

                qflx_snotempunload[p] = max(zero(T), snocan[p] * (forc_t[c] - oftype(snocan[p], 270.15)) /
                                            oftype(snocan[p], snowcan_unload_temp_fact))
                qflx_snowindunload[p] = oftype(snocan[p], snowcan_unload_wind_fact) *
                                        snocan[p] * forc_wind[g] / oftype(snocan[p], 1.56e5)
                qflx_snow_unload[p] = min(qflx_snotempunload[p] + qflx_snowindunload[p],
                                          snocan[p] / dt)
            else
                qflx_snotempunload[p] = zero(T)
                qflx_snowindunload[p] = zero(T)
                qflx_snow_unload[p] = zero(T)
            end
        end
    end
end

canhyd_snow_unloading!(qflx_snotempunload, qflx_snowindunload, qflx_snow_unload,
        mask_soilp, column, gridcell, dtime, frac_veg_nosno, forc_t, forc_wind, snocan,
        snowcan_unload_temp_fact, snowcan_unload_wind_fact, pstart, np) =
    _launch!(_canhyd_snow_unloading_kernel!, qflx_snotempunload, qflx_snowindunload,
        qflx_snow_unload, mask_soilp, column, gridcell,
        convert(eltype(snocan), dtime), frac_veg_nosno,
        forc_t, forc_wind, snocan,
        convert(eltype(snocan), snowcan_unload_temp_fact),
        convert(eltype(snocan), snowcan_unload_wind_fact),
        pstart; ndrange=np)

# -----------------------------------------------------------------------
# TracerFlux_SnowUnloading
# Compute snow unloading for one tracer
# Ported from TracerFlux_SnowUnloading in CanopyHydrologyMod.F90
# -----------------------------------------------------------------------

function tracer_flux_snow_unloading!(
        mask_soilp::AbstractVector{Bool},
        bounds_p::UnitRange{Int},
        # Bulk inputs
        bulk_snocan::AbstractVector{<:Real},
        bulk_qflx_snow_unload::AbstractVector{<:Real},
        # Tracer inputs
        trac_snocan::AbstractVector{<:Real},
        # Tracer outputs
        trac_qflx_snow_unload::AbstractVector{<:Real})

    canhyd_tracer_snow_unloading!(trac_qflx_snow_unload, mask_soilp,
        bulk_snocan, bulk_qflx_snow_unload, trac_snocan,
        first(bounds_p), length(bounds_p))

    return nothing
end

@kernel function _canhyd_tracer_snow_unloading_kernel!(trac_qflx_snow_unload,
        @Const(mask_soilp), @Const(bulk_snocan), @Const(bulk_qflx_snow_unload),
        @Const(trac_snocan), @Const(pstart))
    i = @index(Global)
    @inbounds begin
        p = pstart + i - 1
        if mask_soilp[p]
            if bulk_snocan[p] != 0.0
                trac_qflx_snow_unload[p] = (trac_snocan[p] / bulk_snocan[p]) * bulk_qflx_snow_unload[p]
            else
                trac_qflx_snow_unload[p] = 0.0
            end
        end
    end
end

canhyd_tracer_snow_unloading!(trac_qflx_snow_unload, mask_soilp,
        bulk_snocan, bulk_qflx_snow_unload, trac_snocan, pstart, np) =
    _launch!(_canhyd_tracer_snow_unloading_kernel!, trac_qflx_snow_unload, mask_soilp,
        bulk_snocan, bulk_qflx_snow_unload, trac_snocan, pstart; ndrange=np)

# -----------------------------------------------------------------------
# UpdateState_RemoveSnowUnloading
# Update snocan based on snow unloading, for bulk or one tracer
# Ported from UpdateState_RemoveSnowUnloading in CanopyHydrologyMod.F90
# -----------------------------------------------------------------------

function update_state_remove_snow_unloading!(
        mask_soilp::AbstractVector{Bool},
        bounds_p::UnitRange{Int},
        dtime::Real,
        # Inputs
        qflx_snow_unload::AbstractVector{<:Real},
        # In/Out
        snocan::AbstractVector{<:Real})

    canhyd_remove_snow_unloading!(snocan, mask_soilp, dtime, qflx_snow_unload,
        first(bounds_p), length(bounds_p))

    return nothing
end

@kernel function _canhyd_remove_snow_unloading_kernel!(snocan, @Const(mask_soilp),
        @Const(dtime), @Const(qflx_snow_unload), @Const(pstart))
    i = @index(Global)
    @inbounds begin
        p = pstart + i - 1
        if mask_soilp[p]
            # NOTE(wjs, 2019-05-09) The following check sets snocan to exactly 0 if it would
            # be set to close to 0. The use of exact equality of these floating point values
            # in the conditional is kept for bit-for-bit answers during the refactor.
            dt = oftype(snocan[p], dtime)
            if qflx_snow_unload[p] == snocan[p] / dt
                snocan[p] = zero(snocan[p])
            else
                snocan[p] = snocan[p] - qflx_snow_unload[p] * dt
            end
        end
    end
end

canhyd_remove_snow_unloading!(snocan, mask_soilp, dtime, qflx_snow_unload, pstart, np) =
    _launch!(_canhyd_remove_snow_unloading_kernel!, snocan, mask_soilp,
        convert(eltype(snocan), dtime), qflx_snow_unload, pstart; ndrange=np)

# -----------------------------------------------------------------------
# SumFlux_FluxesOntoGround
# Compute summed fluxes onto ground, for bulk or one tracer
# Ported from SumFlux_FluxesOntoGround in CanopyHydrologyMod.F90
#
# Note: The Fortran version calls p2c (patch-to-column averaging). Here
# we implement a simplified p2c inline using patch weights (wtcol).
# -----------------------------------------------------------------------

function sum_flux_fluxes_onto_ground!(
        patch::PatchData,
        col::ColumnData,
        mask_nolakep::AbstractVector{Bool},
        mask_nolakec::AbstractVector{Bool},
        bounds_p::UnitRange{Int},
        bounds_c::UnitRange{Int},
        # Inputs (patch-level)
        qflx_through_snow::AbstractVector{<:Real},
        qflx_snocanfall::AbstractVector{<:Real},
        qflx_snow_unload::AbstractVector{<:Real},
        qflx_through_liq::AbstractVector{<:Real},
        qflx_liqcanfall::AbstractVector{<:Real},
        qflx_irrig_drip::AbstractVector{<:Real},
        # Outputs (column-level)
        qflx_snow_grnd_col::AbstractVector{<:Real},
        qflx_liq_grnd_col::AbstractVector{<:Real},
        qflx_snow_h2osfc::AbstractVector{<:Real})

    # Initialize column-level accumulators
    canhyd_init_grnd_col!(qflx_snow_grnd_col, qflx_liq_grnd_col, mask_nolakec,
        first(bounds_c), length(bounds_c))

    # Compute patch-level sums and accumulate to column via p2c.
    # This is a patch→column scatter (weighted +=): multiple patches may share a
    # column, so the per-patch kernel uses _scatter_add! (atomic on GPU, plain +=
    # on the sequential CPU/Dual backend) to avoid races on the column accumulators.
    canhyd_fluxes_onto_ground!(qflx_snow_grnd_col, qflx_liq_grnd_col,
        mask_nolakep, patch.column, patch.wtcol,
        qflx_through_snow, qflx_snocanfall, qflx_snow_unload,
        qflx_through_liq, qflx_liqcanfall, qflx_irrig_drip,
        first(bounds_p), length(bounds_p))

    # For now, no snow on surface water
    canhyd_zero_h2osfc!(qflx_snow_h2osfc, mask_nolakec,
        first(bounds_c), length(bounds_c))

    return nothing
end

# Patch→column weighted scatter of snow/liquid fluxes onto the ground.
# One thread per patch; multiple patches write the same column, so use _scatter_add!
# (atomic on a hardware GPU, plain += on the sequential CPU/Dual backend).
@kernel function _canhyd_fluxes_onto_ground_kernel!(qflx_snow_grnd_col, qflx_liq_grnd_col,
        @Const(mask_nolakep), @Const(column), @Const(wtcol),
        @Const(qflx_through_snow), @Const(qflx_snocanfall), @Const(qflx_snow_unload),
        @Const(qflx_through_liq), @Const(qflx_liqcanfall), @Const(qflx_irrig_drip),
        @Const(pstart))
    i = @index(Global)
    @inbounds begin
        p = pstart + i - 1
        if mask_nolakep[p]
            c = column[p]
            w = wtcol[p]
            qflx_snow_grnd_patch = qflx_through_snow[p] + qflx_snocanfall[p] +
                                   qflx_snow_unload[p]
            qflx_liq_grnd_patch = qflx_through_liq[p] + qflx_liqcanfall[p] +
                                  qflx_irrig_drip[p]
            _scatter_add!(qflx_snow_grnd_col, c, qflx_snow_grnd_patch * w)
            _scatter_add!(qflx_liq_grnd_col, c, qflx_liq_grnd_patch * w)
        end
    end
end

canhyd_fluxes_onto_ground!(qflx_snow_grnd_col, qflx_liq_grnd_col, mask_nolakep, column,
        wtcol, qflx_through_snow, qflx_snocanfall, qflx_snow_unload, qflx_through_liq,
        qflx_liqcanfall, qflx_irrig_drip, pstart, np) =
    _launch!(_canhyd_fluxes_onto_ground_kernel!, qflx_snow_grnd_col, qflx_liq_grnd_col,
        mask_nolakep, column, wtcol, qflx_through_snow, qflx_snocanfall, qflx_snow_unload,
        qflx_through_liq, qflx_liqcanfall, qflx_irrig_drip, pstart; ndrange=np)

@kernel function _canhyd_init_grnd_col_kernel!(qflx_snow_grnd_col, qflx_liq_grnd_col,
        @Const(mask_nolakec), @Const(cstart))
    i = @index(Global)
    @inbounds begin
        c = cstart + i - 1
        if mask_nolakec[c]
            qflx_snow_grnd_col[c] = zero(eltype(qflx_snow_grnd_col))
            qflx_liq_grnd_col[c] = zero(eltype(qflx_liq_grnd_col))
        end
    end
end

canhyd_init_grnd_col!(qflx_snow_grnd_col, qflx_liq_grnd_col, mask_nolakec, cstart, nc) =
    _launch!(_canhyd_init_grnd_col_kernel!, qflx_snow_grnd_col, qflx_liq_grnd_col,
        mask_nolakec, cstart; ndrange=nc)

@kernel function _canhyd_zero_h2osfc_kernel!(qflx_snow_h2osfc, @Const(mask_nolakec),
        @Const(cstart))
    i = @index(Global)
    @inbounds begin
        c = cstart + i - 1
        if mask_nolakec[c]
            qflx_snow_h2osfc[c] = zero(eltype(qflx_snow_h2osfc))
        end
    end
end

canhyd_zero_h2osfc!(qflx_snow_h2osfc, mask_nolakec, cstart, nc) =
    _launch!(_canhyd_zero_h2osfc_kernel!, qflx_snow_h2osfc, mask_nolakec, cstart; ndrange=nc)

# -----------------------------------------------------------------------
# BulkDiag_FracWet
# Determine fraction of vegetated surface that is wet
# Ported from BulkDiag_FracWet in CanopyHydrologyMod.F90
# -----------------------------------------------------------------------

function bulk_diag_frac_wet!(
        mask_soilp::AbstractVector{Bool},
        bounds_p::UnitRange{Int},
        # Inputs
        frac_veg_nosno::AbstractVector{<:Integer},
        elai::AbstractVector{<:Real},
        esai::AbstractVector{<:Real},
        snocan::AbstractVector{<:Real},
        liqcan::AbstractVector{<:Real},
        # Outputs
        fwet::AbstractVector{<:Real},
        fdry::AbstractVector{<:Real},
        fcansno::AbstractVector{<:Real})

    canhyd_frac_wet!(fwet, fdry, fcansno, mask_soilp, frac_veg_nosno,
        elai, esai, snocan, liqcan,
        canopy_hydrology_params.liq_canopy_storage_scalar,
        canopy_hydrology_params.snow_canopy_storage_scalar,
        canopy_hydrology_params.maximum_leaf_wetted_fraction,
        first(bounds_p), length(bounds_p))

    return nothing
end

@kernel function _canhyd_frac_wet_kernel!(fwet, fdry, fcansno, @Const(mask_soilp),
        @Const(frac_veg_nosno), @Const(elai), @Const(esai), @Const(snocan), @Const(liqcan),
        @Const(liq_canopy_storage_scalar), @Const(snow_canopy_storage_scalar),
        @Const(maximum_leaf_wetted_fraction), @Const(pstart))
    i = @index(Global)
    @inbounds begin
        p = pstart + i - 1
        if mask_soilp[p]
            T = eltype(fwet)
            if frac_veg_nosno[p] == 1
                h2ocan = snocan[p] + liqcan[p]

                if h2ocan > zero(h2ocan)
                    vegt = frac_veg_nosno[p] * (elai[p] + esai[p])
                    fwet[p] = (h2ocan / (vegt * oftype(h2ocan, liq_canopy_storage_scalar)))^oftype(h2ocan, 0.666666666666)
                    fwet[p] = min(fwet[p], oftype(fwet[p], maximum_leaf_wetted_fraction))
                    if snocan[p] > zero(T)
                        fcansno[p] = (snocan[p] / (vegt * oftype(h2ocan, snow_canopy_storage_scalar)))^oftype(h2ocan, 0.15)
                        fcansno[p] = min(fcansno[p], one(fcansno[p]))
                    else
                        fcansno[p] = zero(T)
                    end
                else
                    fwet[p] = zero(T)
                    fcansno[p] = zero(T)
                end
                fdry[p] = (one(T) - fwet[p]) * elai[p] / (elai[p] + esai[p])
            else
                fwet[p] = zero(T)
                fdry[p] = zero(T)
            end
        end
    end
end

canhyd_frac_wet!(fwet, fdry, fcansno, mask_soilp, frac_veg_nosno,
        elai, esai, snocan, liqcan,
        liq_canopy_storage_scalar, snow_canopy_storage_scalar,
        maximum_leaf_wetted_fraction, pstart, np) =
    _launch!(_canhyd_frac_wet_kernel!, fwet, fdry, fcansno, mask_soilp, frac_veg_nosno,
        elai, esai, snocan, liqcan,
        convert(eltype(fwet), liq_canopy_storage_scalar),
        convert(eltype(fwet), snow_canopy_storage_scalar),
        convert(eltype(fwet), maximum_leaf_wetted_fraction), pstart; ndrange=np)

# -----------------------------------------------------------------------
# Total canopy water diagnostic: h2ocan = snocan + liqcan (Fortran
# CanopyHydrologyMod stores h2ocan_patch each step; the Julia port only
# computed it locally for fwet, leaving the persistent diagnostic at 0).
# -----------------------------------------------------------------------
@kernel function _canhyd_h2ocan_kernel!(h2ocan, @Const(snocan), @Const(liqcan),
        @Const(mask), @Const(pstart))
    i = @index(Global)
    @inbounds begin
        p = pstart + i - 1
        if mask[p]
            h2ocan[p] = snocan[p] + liqcan[p]
        end
    end
end
canhyd_update_h2ocan!(h2ocan, snocan, liqcan, mask, pstart, np) =
    _launch!(_canhyd_h2ocan_kernel!, h2ocan, snocan, liqcan, mask, pstart; ndrange=np)

# -----------------------------------------------------------------------
# CanopyInterceptionAndThroughfall (top-level coordinator)
# Ported from CanopyInterceptionAndThroughfall in CanopyHydrologyMod.F90
#
# Coordinates all sub-routines for canopy hydrology including
# interception, throughfall, canopy excess, snow unloading, and
# summation of fluxes onto the ground.
# -----------------------------------------------------------------------

function canopy_interception_and_throughfall!(
        patch::PatchData,
        col::ColumnData,
        canopystate::CanopyStateData,
        water::WaterData,
        dtime::Real,
        # Masks
        mask_soilp::AbstractVector{Bool},
        mask_nolakep::AbstractVector{Bool},
        mask_nolakec::AbstractVector{Bool},
        # Bounds
        bounds_p::UnitRange{Int},
        bounds_c::UnitRange{Int},
        bounds_g::UnitRange{Int},
        # Atmospheric forcing (column-level)
        forc_rain::AbstractVector{<:Real},
        forc_snow_col::AbstractVector{<:Real},
        forc_t::AbstractVector{<:Real},
        # Atmospheric forcing (gridcell-level)
        forc_wind::AbstractVector{<:Real},
        # Irrigation (patch-level)
        qflx_irrig_sprinkler::AbstractVector{<:Real},
        qflx_irrig_drip::AbstractVector{<:Real})

    np = length(bounds_p)
    nc = length(bounds_c)

    b_wf = water.waterfluxbulk_inst
    b_ws = water.waterstatebulk_inst
    b_wd = water.waterdiagnosticbulk_inst

    # Local scratch — allocate via similar() of a device-resident state array so the
    # scratch lives on the same backend (CPU/Metal) as the rest of the call. FT is the
    # working element type (Float64 on CPU, Float32 on Metal, Dual under AD).
    FT = eltype(forc_rain)
    proto = b_wf.wf.qflx_through_snow_patch
    qflx_liq_above_canopy_patch = fill!(similar(proto, FT, last(bounds_p)), zero(FT))
    forc_snow_patch             = fill!(similar(proto, FT, last(bounds_p)), zero(FT))
    check_point_for_interception_and_excess =
        fill!(similar(proto, Bool, last(bounds_p)), false)

    # The canopy-fall fluxes (snocanfall/liqcanfall/snow_unload) are computed only over
    # mask_soilp, but the fluxes-onto-ground scatter reads them over mask_nolakep. Special
    # land units (e.g. glacier/istice) are in nolakep but not soilp, so their entries are
    # never set -> NaN -> poisons qflx_snow_grnd -> int_snow/h2osno. Zero them here so
    # non-soil patches contribute 0 canopy fall (they have no canopy); the soilp kernels
    # overwrite the soil patches.
    fill!(b_wf.wf.qflx_snocanfall_patch, zero(FT))
    fill!(b_wf.wf.qflx_liqcanfall_patch, zero(FT))
    fill!(b_wf.wf.qflx_snow_unload_patch, zero(FT))

    # --- Step 1: Top-of-canopy inputs for bulk ---
    sum_flux_top_of_canopy_inputs!(
        patch, mask_nolakep, bounds_p,
        forc_rain, forc_snow_col, qflx_irrig_sprinkler,
        qflx_liq_above_canopy_patch, forc_snow_patch)

    # --- Step 2: Bulk interception and throughfall ---
    bulk_flux_canopy_interception_and_throughfall!(
        patch, col, mask_nolakep, bounds_p,
        canopystate.frac_veg_nosno_patch,
        canopystate.elai_patch, canopystate.esai_patch,
        forc_snow_patch, qflx_liq_above_canopy_patch,
        b_wf.wf.qflx_through_snow_patch,
        b_wf.wf.qflx_through_liq_patch,
        b_wf.wf.qflx_intercepted_snow_patch,
        b_wf.wf.qflx_intercepted_liq_patch,
        check_point_for_interception_and_excess)

    # --- Step 3: Add interception to canopy for bulk ---
    update_state_add_interception_to_canopy!(
        mask_soilp, bounds_p, dtime,
        b_wf.wf.qflx_intercepted_snow_patch,
        b_wf.wf.qflx_intercepted_liq_patch,
        b_ws.ws.snocan_patch,
        b_ws.ws.liqcan_patch)

    # --- Step 4: Canopy excess for bulk ---
    bulk_flux_canopy_excess!(
        mask_soilp, bounds_p, dtime,
        canopystate.elai_patch, canopystate.esai_patch,
        b_ws.ws.snocan_patch, b_ws.ws.liqcan_patch,
        check_point_for_interception_and_excess,
        b_wf.wf.qflx_snocanfall_patch,
        b_wf.wf.qflx_liqcanfall_patch)

    # --- Step 5: Remove canfall from canopy for bulk ---
    update_state_remove_canfall_from_canopy!(
        mask_soilp, bounds_p, dtime,
        b_wf.wf.qflx_liqcanfall_patch,
        b_wf.wf.qflx_snocanfall_patch,
        b_ws.ws.liqcan_patch,
        b_ws.ws.snocan_patch)

    # --- Step 6: Snow unloading for bulk ---
    bulk_flux_snow_unloading!(
        patch, mask_soilp, bounds_p, dtime,
        canopystate.frac_veg_nosno_patch,
        forc_t, forc_wind,
        b_ws.ws.snocan_patch,
        b_wf.qflx_snotempunload_patch,
        b_wf.qflx_snowindunload_patch,
        b_wf.wf.qflx_snow_unload_patch)

    # --- Step 7: Remove snow unloading from canopy for bulk ---
    update_state_remove_snow_unloading!(
        mask_soilp, bounds_p, dtime,
        b_wf.wf.qflx_snow_unload_patch,
        b_ws.ws.snocan_patch)

    # --- Step 8: Sum fluxes onto ground for bulk ---
    sum_flux_fluxes_onto_ground!(
        patch, col, mask_nolakep, mask_nolakec,
        bounds_p, bounds_c,
        b_wf.wf.qflx_through_snow_patch,
        b_wf.wf.qflx_snocanfall_patch,
        b_wf.wf.qflx_snow_unload_patch,
        b_wf.wf.qflx_through_liq_patch,
        b_wf.wf.qflx_liqcanfall_patch,
        qflx_irrig_drip,
        b_wf.wf.qflx_snow_grnd_col,
        b_wf.wf.qflx_liq_grnd_col,
        b_wf.wf.qflx_snow_h2osfc_col)

    # --- Step 9: Fraction wet/dry for bulk ---
    bulk_diag_frac_wet!(
        mask_soilp, bounds_p,
        canopystate.frac_veg_nosno_patch,
        canopystate.elai_patch, canopystate.esai_patch,
        b_ws.ws.snocan_patch, b_ws.ws.liqcan_patch,
        b_wd.fwet_patch, b_wd.fdry_patch, b_wd.fcansno_patch)

    # NOTE: the total canopy water diagnostic h2ocan = snocan + liqcan is synced
    # AFTER canopy_fluxes! (see canhyd_update_h2ocan! call in clm_driver), because
    # canopy_fluxes applies the canopy-evaporation decrement to snocan/liqcan; the
    # end-of-hydrology values here are pre-evaporation and would read ~15% high.

    return nothing
end
