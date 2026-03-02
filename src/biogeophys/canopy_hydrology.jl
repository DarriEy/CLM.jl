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
        liq_canopy_storage_scalar::Float64 = 0.04,
        snow_canopy_storage_scalar::Float64 = 6.0,
        snowcan_unload_temp_fact::Float64 = 1.87e5,
        snowcan_unload_wind_fact::Float64 = 1.56e5,
        interception_fraction::Float64 = 1.0,
        maximum_leaf_wetted_fraction::Float64 = 1.0)
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

function sum_flux_top_of_canopy_inputs!(
        patch::PatchData,
        mask_nolakep::BitVector,
        bounds_p::UnitRange{Int},
        # Inputs (column-level)
        forc_rain::Vector{Float64},
        forc_snow_col::Vector{Float64},
        # Inputs (patch-level)
        qflx_irrig_sprinkler::Vector{Float64},
        # Outputs (patch-level)
        qflx_liq_above_canopy::Vector{Float64},
        forc_snow_patch::Vector{Float64})

    for p in bounds_p
        mask_nolakep[p] || continue
        c = patch.column[p]

        qflx_liq_above_canopy[p] = forc_rain[c] + qflx_irrig_sprinkler[p]
        forc_snow_patch[p] = forc_snow_col[c]
    end

    return nothing
end

# -----------------------------------------------------------------------
# BulkFlux_CanopyInterceptionAndThroughfall
# Compute canopy interception and throughfall for bulk water
# Ported from BulkFlux_CanopyInterceptionAndThroughfall in CanopyHydrologyMod.F90
# -----------------------------------------------------------------------

function bulk_flux_canopy_interception_and_throughfall!(
        patch::PatchData,
        col::ColumnData,
        mask_nolakep::BitVector,
        bounds_p::UnitRange{Int},
        # Inputs
        frac_veg_nosno::Vector{Int},
        elai::Vector{Float64},
        esai::Vector{Float64},
        forc_snow::Vector{Float64},
        qflx_liq_above_canopy::Vector{Float64},
        # Outputs
        qflx_through_snow::Vector{Float64},
        qflx_through_liq::Vector{Float64},
        qflx_intercepted_snow::Vector{Float64},
        qflx_intercepted_liq::Vector{Float64},
        check_point_for_interception_and_excess::Vector{Bool})

    for p in bounds_p
        mask_nolakep[p] || continue
        c = patch.column[p]

        check_point_for_interception_and_excess[p] =
            (frac_veg_nosno[p] == 1 && (forc_snow[p] + qflx_liq_above_canopy[p]) > 0.0)

        if check_point_for_interception_and_excess[p]
            # Coefficient of interception
            if canopy_hydrology_ctrl.use_clm5_fpi
                fpiliq = canopy_hydrology_params.interception_fraction * tanh(elai[p] + esai[p])
            else
                fpiliq = 0.25 * (1.0 - exp(-0.5 * (elai[p] + esai[p])))
            end

            fpisnow = (1.0 - exp(-0.5 * (elai[p] + esai[p])))  # max interception of 1

            # Direct throughfall
            qflx_through_snow[p] = forc_snow[p] * (1.0 - fpisnow)
            qflx_through_liq[p]  = qflx_liq_above_canopy[p] * (1.0 - fpiliq)

            # Canopy interception
            qflx_intercepted_snow[p] = forc_snow[p] * fpisnow
            qflx_intercepted_liq[p]  = qflx_liq_above_canopy[p] * fpiliq
        else
            # Note that special landunits will be handled here, in addition to soil points
            # with frac_veg_nosno == 0.
            qflx_intercepted_snow[p] = 0.0
            qflx_intercepted_liq[p]  = 0.0
            if col.itype[c] == ICOL_SUNWALL || col.itype[c] == ICOL_SHADEWALL
                qflx_through_snow[p] = 0.0
                qflx_through_liq[p]  = 0.0
            else
                qflx_through_snow[p] = forc_snow[p]
                qflx_through_liq[p]  = qflx_liq_above_canopy[p]
            end
        end
    end

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
        mask_nolakep::BitVector,
        bounds_p::UnitRange{Int},
        # Bulk inputs
        bulk_forc_snow::Vector{Float64},
        bulk_qflx_liq_above_canopy::Vector{Float64},
        bulk_qflx_through_snow::Vector{Float64},
        bulk_qflx_intercepted_snow::Vector{Float64},
        bulk_qflx_through_liq::Vector{Float64},
        bulk_qflx_intercepted_liq::Vector{Float64},
        # Tracer inputs
        trac_forc_snow::Vector{Float64},
        trac_qflx_liq_above_canopy::Vector{Float64},
        # Tracer outputs
        trac_qflx_through_snow::Vector{Float64},
        trac_qflx_intercepted_snow::Vector{Float64},
        trac_qflx_through_liq::Vector{Float64},
        trac_qflx_intercepted_liq::Vector{Float64})

    for p in bounds_p
        mask_nolakep[p] || continue

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

    return nothing
end

# -----------------------------------------------------------------------
# UpdateState_AddInterceptionToCanopy
# Update snocan and liqcan based on interception, for bulk or one tracer
# Ported from UpdateState_AddInterceptionToCanopy in CanopyHydrologyMod.F90
# -----------------------------------------------------------------------

function update_state_add_interception_to_canopy!(
        mask_soilp::BitVector,
        bounds_p::UnitRange{Int},
        dtime::Float64,
        # Inputs
        qflx_intercepted_snow::Vector{Float64},
        qflx_intercepted_liq::Vector{Float64},
        # In/Out
        snocan::Vector{Float64},
        liqcan::Vector{Float64})

    for p in bounds_p
        mask_soilp[p] || continue

        snocan[p] = max(0.0, snocan[p] + dtime * qflx_intercepted_snow[p])
        liqcan[p] = max(0.0, liqcan[p] + dtime * qflx_intercepted_liq[p])
    end

    return nothing
end

# -----------------------------------------------------------------------
# BulkFlux_CanopyExcess
# Compute runoff from canopy due to exceeding maximum storage, for bulk
# Ported from BulkFlux_CanopyExcess in CanopyHydrologyMod.F90
# -----------------------------------------------------------------------

function bulk_flux_canopy_excess!(
        mask_soilp::BitVector,
        bounds_p::UnitRange{Int},
        dtime::Float64,
        # Inputs
        elai::Vector{Float64},
        esai::Vector{Float64},
        snocan::Vector{Float64},
        liqcan::Vector{Float64},
        check_point_for_interception_and_excess::Vector{Bool},
        # Outputs
        qflx_snocanfall::Vector{Float64},
        qflx_liqcanfall::Vector{Float64})

    for p in bounds_p
        mask_soilp[p] || continue

        qflx_liqcanfall[p] = 0.0
        qflx_snocanfall[p] = 0.0

        if check_point_for_interception_and_excess[p]
            liqcanmx = canopy_hydrology_params.liq_canopy_storage_scalar * (elai[p] + esai[p])
            qflx_liqcanfall[p] = max((liqcan[p] - liqcanmx) / dtime, 0.0)
            snocanmx = canopy_hydrology_params.snow_canopy_storage_scalar * (elai[p] + esai[p])
            qflx_snocanfall[p] = max((snocan[p] - snocanmx) / dtime, 0.0)
        end
    end

    return nothing
end

# -----------------------------------------------------------------------
# TracerFlux_CanopyExcess
# Calculate runoff from canopy due to exceeding maximum storage, for one tracer
# Ported from TracerFlux_CanopyExcess in CanopyHydrologyMod.F90
# -----------------------------------------------------------------------

function tracer_flux_canopy_excess!(
        mask_soilp::BitVector,
        bounds_p::UnitRange{Int},
        # Bulk inputs
        bulk_liqcan::Vector{Float64},
        bulk_snocan::Vector{Float64},
        bulk_qflx_liqcanfall::Vector{Float64},
        bulk_qflx_snocanfall::Vector{Float64},
        # Tracer inputs
        trac_liqcan::Vector{Float64},
        trac_snocan::Vector{Float64},
        # Tracer outputs
        trac_qflx_liqcanfall::Vector{Float64},
        trac_qflx_snocanfall::Vector{Float64})

    for p in bounds_p
        mask_soilp[p] || continue

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

    return nothing
end

# -----------------------------------------------------------------------
# UpdateState_RemoveCanfallFromCanopy
# Update snocan and liqcan based on canfall, for bulk or one tracer
# Ported from UpdateState_RemoveCanfallFromCanopy in CanopyHydrologyMod.F90
# -----------------------------------------------------------------------

function update_state_remove_canfall_from_canopy!(
        mask_soilp::BitVector,
        bounds_p::UnitRange{Int},
        dtime::Float64,
        # Inputs
        qflx_liqcanfall::Vector{Float64},
        qflx_snocanfall::Vector{Float64},
        # In/Out
        liqcan::Vector{Float64},
        snocan::Vector{Float64})

    for p in bounds_p
        mask_soilp[p] || continue

        liqcan[p] = liqcan[p] - dtime * qflx_liqcanfall[p]
        snocan[p] = snocan[p] - dtime * qflx_snocanfall[p]
    end

    return nothing
end

# -----------------------------------------------------------------------
# BulkFlux_SnowUnloading
# Compute snow unloading for bulk
# Ported from BulkFlux_SnowUnloading in CanopyHydrologyMod.F90
# -----------------------------------------------------------------------

function bulk_flux_snow_unloading!(
        patch::PatchData,
        mask_soilp::BitVector,
        bounds_p::UnitRange{Int},
        dtime::Float64,
        # Inputs
        frac_veg_nosno::Vector{Int},
        forc_t::Vector{Float64},        # column-level
        forc_wind::Vector{Float64},     # gridcell-level
        snocan::Vector{Float64},
        # Outputs
        qflx_snotempunload::Vector{Float64},
        qflx_snowindunload::Vector{Float64},
        qflx_snow_unload::Vector{Float64})

    for p in bounds_p
        mask_soilp[p] || continue

        if frac_veg_nosno[p] == 1 && snocan[p] > 0.0
            c = patch.column[p]
            g = patch.gridcell[p]

            qflx_snotempunload[p] = max(0.0, snocan[p] * (forc_t[c] - 270.15) /
                                        canopy_hydrology_params.snowcan_unload_temp_fact)
            qflx_snowindunload[p] = canopy_hydrology_params.snowcan_unload_wind_fact *
                                    snocan[p] * forc_wind[g] / 1.56e5
            qflx_snow_unload[p] = min(qflx_snotempunload[p] + qflx_snowindunload[p],
                                      snocan[p] / dtime)
        else
            qflx_snotempunload[p] = 0.0
            qflx_snowindunload[p] = 0.0
            qflx_snow_unload[p] = 0.0
        end
    end

    return nothing
end

# -----------------------------------------------------------------------
# TracerFlux_SnowUnloading
# Compute snow unloading for one tracer
# Ported from TracerFlux_SnowUnloading in CanopyHydrologyMod.F90
# -----------------------------------------------------------------------

function tracer_flux_snow_unloading!(
        mask_soilp::BitVector,
        bounds_p::UnitRange{Int},
        # Bulk inputs
        bulk_snocan::Vector{Float64},
        bulk_qflx_snow_unload::Vector{Float64},
        # Tracer inputs
        trac_snocan::Vector{Float64},
        # Tracer outputs
        trac_qflx_snow_unload::Vector{Float64})

    for p in bounds_p
        mask_soilp[p] || continue

        if bulk_snocan[p] != 0.0
            trac_qflx_snow_unload[p] = (trac_snocan[p] / bulk_snocan[p]) * bulk_qflx_snow_unload[p]
        else
            trac_qflx_snow_unload[p] = 0.0
        end
    end

    return nothing
end

# -----------------------------------------------------------------------
# UpdateState_RemoveSnowUnloading
# Update snocan based on snow unloading, for bulk or one tracer
# Ported from UpdateState_RemoveSnowUnloading in CanopyHydrologyMod.F90
# -----------------------------------------------------------------------

function update_state_remove_snow_unloading!(
        mask_soilp::BitVector,
        bounds_p::UnitRange{Int},
        dtime::Float64,
        # Inputs
        qflx_snow_unload::Vector{Float64},
        # In/Out
        snocan::Vector{Float64})

    for p in bounds_p
        mask_soilp[p] || continue

        # NOTE(wjs, 2019-05-09) The following check sets snocan to exactly 0 if it would
        # be set to close to 0. The use of exact equality of these floating point values
        # in the conditional is kept for bit-for-bit answers during the refactor.
        if qflx_snow_unload[p] == snocan[p] / dtime
            snocan[p] = 0.0
        else
            snocan[p] = snocan[p] - qflx_snow_unload[p] * dtime
        end
    end

    return nothing
end

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
        mask_nolakep::BitVector,
        mask_nolakec::BitVector,
        bounds_p::UnitRange{Int},
        bounds_c::UnitRange{Int},
        # Inputs (patch-level)
        qflx_through_snow::Vector{Float64},
        qflx_snocanfall::Vector{Float64},
        qflx_snow_unload::Vector{Float64},
        qflx_through_liq::Vector{Float64},
        qflx_liqcanfall::Vector{Float64},
        qflx_irrig_drip::Vector{Float64},
        # Outputs (column-level)
        qflx_snow_grnd_col::Vector{Float64},
        qflx_liq_grnd_col::Vector{Float64},
        qflx_snow_h2osfc::Vector{Float64})

    # Initialize column-level accumulators
    for c in bounds_c
        mask_nolakec[c] || continue
        qflx_snow_grnd_col[c] = 0.0
        qflx_liq_grnd_col[c] = 0.0
    end

    # Compute patch-level sums and accumulate to column via p2c
    for p in bounds_p
        mask_nolakep[p] || continue

        qflx_snow_grnd_patch = qflx_through_snow[p] +
                               qflx_snocanfall[p] +
                               qflx_snow_unload[p]

        qflx_liq_grnd_patch = qflx_through_liq[p] +
                              qflx_liqcanfall[p] +
                              qflx_irrig_drip[p]

        # Patch-to-column weighted average (inline p2c)
        c = patch.column[p]
        qflx_snow_grnd_col[c] += qflx_snow_grnd_patch * patch.wtcol[p]
        qflx_liq_grnd_col[c]  += qflx_liq_grnd_patch * patch.wtcol[p]
    end

    # For now, no snow on surface water
    for c in bounds_c
        mask_nolakec[c] || continue
        qflx_snow_h2osfc[c] = 0.0
    end

    return nothing
end

# -----------------------------------------------------------------------
# BulkDiag_FracWet
# Determine fraction of vegetated surface that is wet
# Ported from BulkDiag_FracWet in CanopyHydrologyMod.F90
# -----------------------------------------------------------------------

function bulk_diag_frac_wet!(
        mask_soilp::BitVector,
        bounds_p::UnitRange{Int},
        # Inputs
        frac_veg_nosno::Vector{Int},
        elai::Vector{Float64},
        esai::Vector{Float64},
        snocan::Vector{Float64},
        liqcan::Vector{Float64},
        # Outputs
        fwet::Vector{Float64},
        fdry::Vector{Float64},
        fcansno::Vector{Float64})

    for p in bounds_p
        mask_soilp[p] || continue

        if frac_veg_nosno[p] == 1
            h2ocan = snocan[p] + liqcan[p]

            if h2ocan > 0.0
                vegt = frac_veg_nosno[p] * (elai[p] + esai[p])
                fwet[p] = (h2ocan / (vegt * canopy_hydrology_params.liq_canopy_storage_scalar))^0.666666666666
                fwet[p] = min(fwet[p], canopy_hydrology_params.maximum_leaf_wetted_fraction)
                if snocan[p] > 0.0
                    fcansno[p] = (snocan[p] / (vegt * canopy_hydrology_params.snow_canopy_storage_scalar))^0.15
                    fcansno[p] = min(fcansno[p], 1.0)
                else
                    fcansno[p] = 0.0
                end
            else
                fwet[p] = 0.0
                fcansno[p] = 0.0
            end
            fdry[p] = (1.0 - fwet[p]) * elai[p] / (elai[p] + esai[p])
        else
            fwet[p] = 0.0
            fdry[p] = 0.0
        end
    end

    return nothing
end

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
        dtime::Float64,
        # Masks
        mask_soilp::BitVector,
        mask_nolakep::BitVector,
        mask_nolakec::BitVector,
        # Bounds
        bounds_p::UnitRange{Int},
        bounds_c::UnitRange{Int},
        bounds_g::UnitRange{Int},
        # Atmospheric forcing (column-level)
        forc_rain::Vector{Float64},
        forc_snow_col::Vector{Float64},
        forc_t::Vector{Float64},
        # Atmospheric forcing (gridcell-level)
        forc_wind::Vector{Float64},
        # Irrigation (patch-level)
        qflx_irrig_sprinkler::Vector{Float64},
        qflx_irrig_drip::Vector{Float64})

    np = length(bounds_p)
    nc = length(bounds_c)

    # Local arrays
    qflx_liq_above_canopy_patch = zeros(Float64, last(bounds_p))
    forc_snow_patch = zeros(Float64, last(bounds_p))
    check_point_for_interception_and_excess = fill(false, last(bounds_p))

    b_wf = water.waterfluxbulk_inst
    b_ws = water.waterstatebulk_inst
    b_wd = water.waterdiagnosticbulk_inst

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

    return nothing
end
