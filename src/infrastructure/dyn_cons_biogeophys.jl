# ==========================================================================
# Ported from: src/dyn_subgrid/dynConsBiogeophysMod.F90 (794 lines)
#
# Conserves WATER and ENERGY (heat) across land-cover area changes. Computes
# grid-level total heat & water content BEFORE and AFTER a land-cover change;
# the difference (an otherwise-spurious source/sink) becomes a "dynbal"
# adjustment flux (qflx_liq_dynbal / qflx_ice_dynbal / eflx_dynbal) so the
# budget closes.
#
# This module does NOT update column/patch state — it READS column state arrays
# (h2osoi_liq/ice, t_soisno, snow, h2osfc, lake, ...) and integrates content
# using the already-ported TotalWaterAndHeatMod helpers (total_water_heat.jl):
#   compute_liq_ice_mass_non_lake! / compute_liq_ice_mass_lake!
#   compute_heat_non_lake!         / compute_heat_lake!
#   accumulate_soil_liq_ice_mass_non_lake! / accumulate_soil_heat_non_lake!
#   accumulate_liq_ice_mass_lake!  / accumulate_heat_lake!
#   adjust_delta_heat_for_delta_liq! / heat_base_temp
# then aggregates column -> gridcell via c2g_1d! (subgrid_ave.jl).
#
# Translation notes vs. Fortran:
# - Fortran integer column filters become BitVector masks (mask_nolakec /
#   mask_lakec / mask_icec), per CLAUDE.md conventions.
# - The Fortran code loops over water_inst%bulk_and_tracers. This port handles
#   bulk water only (the single waterstatebulk_inst / waterdiagnosticbulk_inst),
#   matching the rest of the CLM.jl port (no isotopic tracers ported yet).
# - The "dribbler" (eflx_dynbal_dribbler / qflx_*_dynbal_dribbler) in CTSM is a
#   pass-through smoother. With the default annual glacier dynamics and no
#   time-dribbling configured, set_curr_delta + get_curr_flux returns the delta
#   directly as the flux for the step. We therefore set the dynbal flux equal to
#   the gridcell content delta (zeroed when get_for_testing_zero_dynbal_fluxes).
# - The before/after totals and the produced dynbal fluxes are stored in a
#   standalone `DynConsBiogeophysState` struct — NOT added to CLMInstances or any
#   ForwardDiff-dual-copied struct (per the task constraints).
#
# Skipped: water tracers (bulk_and_tracers loop), the dribbler time-smoothing
#   object, and the urban-params heat path is delegated to the already-ported
#   compute_heat_non_lake! helper.
# ==========================================================================

"""
    DynConsBiogeophysState

Holds the grid-level "before" / "after" total heat & water contents and the
resulting dynbal conservation fluxes for one bulk-water accounting.

Mirrors the relevant Fortran fields that live on `waterbalance_type`,
`temperature_type`, `waterflux_type` and `energyflux_type`:
- `liq1_grc` / `ice1_grc`     : gridcell liquid/ice content BEFORE land cover change (kg/m^2)
- `liq2_grc` / `ice2_grc`     : gridcell liquid/ice content AFTER  land cover change (kg/m^2)
- `heat1_grc` / `heat2_grc`   : gridcell heat content before/after (J/m^2)
- `liquid_water_temp1_grc` / `liquid_water_temp2_grc` : gridcell weighted-avg liquid water temperature before/after (K)
- `qflx_liq_dynbal_grc` / `qflx_ice_dynbal_grc` : dynbal liquid/ice water flux (kg/m^2 = mm H2O)
- `eflx_dynbal_grc`           : dynbal energy (heat) flux (J/m^2)

This object is standalone and is NOT part of `CLMInstances`.
"""
Base.@kwdef mutable struct DynConsBiogeophysState{V<:AbstractVector{<:Real}}
    liq1_grc::V = Float64[]
    ice1_grc::V = Float64[]
    liq2_grc::V = Float64[]
    ice2_grc::V = Float64[]
    heat1_grc::V = Float64[]
    heat2_grc::V = Float64[]
    liquid_water_temp1_grc::V = Float64[]
    liquid_water_temp2_grc::V = Float64[]
    qflx_liq_dynbal_grc::V = Float64[]
    qflx_ice_dynbal_grc::V = Float64[]
    eflx_dynbal_grc::V = Float64[]
end

"""
    dyn_cons_biogeophys_state_init(ng::Int) -> DynConsBiogeophysState

Allocate a `DynConsBiogeophysState` with `ng` gridcells, all fields zeroed.
"""
function dyn_cons_biogeophys_state_init(ng::Int)
    z() = zeros(Float64, ng)
    return DynConsBiogeophysState(
        liq1_grc = z(), ice1_grc = z(),
        liq2_grc = z(), ice2_grc = z(),
        heat1_grc = z(), heat2_grc = z(),
        liquid_water_temp1_grc = z(), liquid_water_temp2_grc = z(),
        qflx_liq_dynbal_grc = z(), qflx_ice_dynbal_grc = z(),
        eflx_dynbal_grc = z(),
    )
end

# ---------------------------------------------------------------------------
# Private: dyn_water_content
# ---------------------------------------------------------------------------

"""
    dyn_water_content!(liquid_mass, ice_mass, bounds, mask_nolakec, mask_lakec,
                       col, lun, waterstate, waterdiagnostic, lakestate)

Compute gridcell total liquid and ice water contents (kg/m^2), aggregating
non-lake + lake column content to the gridcell via `c2g_1d!`.

Outputs are written into `liquid_mass[begg:endg]` and `ice_mass[begg:endg]`.

Ported from `dyn_water_content` in `dynConsBiogeophysMod.F90`.
"""
function dyn_water_content!(liquid_mass::AbstractVector{<:Real},
                            ice_mass::AbstractVector{<:Real},
                            bounds::BoundsType,
                            mask_nolakec::AbstractVector{Bool},
                            mask_lakec::AbstractVector{Bool},
                            col::ColumnData, lun::LandunitData,
                            waterstate::WaterStateData,
                            waterdiagnostic::WaterDiagnosticBulkData,
                            lakestate::LakeStateData)
    nc = bounds.endc
    liquid_mass_col = zeros(Float64, nc)
    ice_mass_col    = zeros(Float64, nc)

    # Non-lake column content (subtract dynbal baselines)
    compute_liq_ice_mass_non_lake!(mask_nolakec, col, waterstate, waterdiagnostic,
                                   true, liquid_mass_col, ice_mass_col)

    # Lake column content (add lake water and subtract dynbal baselines)
    compute_liq_ice_mass_lake!(mask_lakec, col, waterstate, lakestate,
                               true, liquid_mass_col, ice_mass_col)

    # Aggregate column -> gridcell
    c2g_1d!(liquid_mass, liquid_mass_col, bounds, "urbanf", "unity", col, lun)
    c2g_1d!(ice_mass,    ice_mass_col,    bounds, "urbanf", "unity", col, lun)

    return nothing
end

# ---------------------------------------------------------------------------
# Private: dyn_heat_content
# ---------------------------------------------------------------------------

"""
    dyn_heat_content!(heat_grc, liquid_water_temp_grc, bounds, mask_nolakec,
                      mask_lakec, col, lun, urbanparams, soilstate, temperature,
                      waterstatebulk, waterdiagnosticbulk, lakestate)

Compute grid-level total heat content (J/m^2) and the weighted-average liquid
water temperature (K). Heat content is relative to a baseline of 0 C
(`heat_base_temp = tfrz`).

Outputs written into `heat_grc[begg:endg]` and `liquid_water_temp_grc[begg:endg]`.

Ported from `dyn_heat_content` in `dynConsBiogeophysMod.F90`.
"""
function dyn_heat_content!(heat_grc::AbstractVector{<:Real},
                           liquid_water_temp_grc::AbstractVector{<:Real},
                           bounds::BoundsType,
                           mask_nolakec::AbstractVector{Bool},
                           mask_lakec::AbstractVector{Bool},
                           col::ColumnData, lun::LandunitData,
                           urbanparams::UrbanParamsData,
                           soilstate::SoilStateData,
                           temperature::TemperatureData,
                           waterstatebulk::WaterStateBulkData,
                           waterdiagnosticbulk::WaterDiagnosticBulkData,
                           lakestate::LakeStateData)
    nc = bounds.endc
    heat_col        = zeros(Float64, nc)
    heat_liquid_col = zeros(Float64, nc)
    cv_liquid_col   = zeros(Float64, nc)

    compute_heat_non_lake!(mask_nolakec, col, lun, urbanparams, soilstate,
                           temperature, waterstatebulk, waterdiagnosticbulk,
                           heat_col, heat_liquid_col, cv_liquid_col)

    compute_heat_lake!(mask_lakec, col, soilstate, temperature,
                       waterstatebulk, lakestate,
                       heat_col, heat_liquid_col, cv_liquid_col)

    # Aggregate column -> gridcell
    heat_liquid_grc = zeros(Float64, bounds.endg)
    cv_liquid_grc   = zeros(Float64, bounds.endg)
    c2g_1d!(heat_grc,        heat_col,        bounds, "urbanf", "unity", col, lun)
    c2g_1d!(heat_liquid_grc, heat_liquid_col, bounds, "urbanf", "unity", col, lun)
    c2g_1d!(cv_liquid_grc,   cv_liquid_col,   bounds, "urbanf", "unity", col, lun)

    # Weighted-average liquid water temperature
    for g in bounds.begg:bounds.endg
        if cv_liquid_grc[g] > 0.0
            liquid_water_temp_grc[g] =
                (heat_liquid_grc[g] / cv_liquid_grc[g]) + heat_base_temp
        else
            # 0 or negative water mass: arbitrary temperature
            liquid_water_temp_grc[g] = TFRZ
        end
    end

    return nothing
end

# ---------------------------------------------------------------------------
# Public: dyn_hwcontent_init
# ---------------------------------------------------------------------------

"""
    dyn_hwcontent_init!(dynbal, bounds, mask_nolakec, mask_lakec, col, lun,
                        urbanparams, soilstate, water, temperature, lakestate)

Compute grid cell-level heat and water content BEFORE land cover change and
store them as the "before" (`*1_grc`) totals on `dynbal`.

Should be called BEFORE any subgrid weight updates this time step.

Ported from `dyn_hwcontent_init` in `dynConsBiogeophysMod.F90`.
"""
function dyn_hwcontent_init!(dynbal::DynConsBiogeophysState,
                             bounds::BoundsType,
                             mask_nolakec::AbstractVector{Bool},
                             mask_lakec::AbstractVector{Bool},
                             col::ColumnData, lun::LandunitData,
                             urbanparams::UrbanParamsData,
                             soilstate::SoilStateData,
                             waterstatebulk::WaterStateBulkData,
                             waterdiagnosticbulk::WaterDiagnosticBulkData,
                             temperature::TemperatureData,
                             lakestate::LakeStateData)
    dyn_water_content!(dynbal.liq1_grc, dynbal.ice1_grc, bounds,
                       mask_nolakec, mask_lakec, col, lun,
                       waterstatebulk.ws, waterdiagnosticbulk, lakestate)

    dyn_heat_content!(dynbal.heat1_grc, dynbal.liquid_water_temp1_grc, bounds,
                      mask_nolakec, mask_lakec, col, lun, urbanparams, soilstate,
                      temperature, waterstatebulk, waterdiagnosticbulk, lakestate)

    return nothing
end

# ---------------------------------------------------------------------------
# Private: dyn_water_content_final
# ---------------------------------------------------------------------------

"""
    dyn_water_content_final!(delta_liq, dynbal, bounds, mask_nolakec, mask_lakec,
                             col, lun, waterstate, waterdiagnostic, lakestate, ctl)

Compute grid cell-level water content AFTER land cover change (`*2_grc`) and
the dynbal water fluxes from (after - before). `delta_liq` is returned (the
change in gridcell liquid content), to be passed to the heat adjustment.

When `get_for_testing_zero_dynbal_fluxes(ctl)` is true, the deltas and fluxes
are zeroed.

Ported from `dyn_water_content_final` in `dynConsBiogeophysMod.F90`.
"""
function dyn_water_content_final!(delta_liq::AbstractVector{<:Real},
                                  dynbal::DynConsBiogeophysState,
                                  bounds::BoundsType,
                                  mask_nolakec::AbstractVector{Bool},
                                  mask_lakec::AbstractVector{Bool},
                                  col::ColumnData, lun::LandunitData,
                                  waterstate::WaterStateData,
                                  waterdiagnostic::WaterDiagnosticBulkData,
                                  lakestate::LakeStateData,
                                  ctl::DynSubgridControl)
    dyn_water_content!(dynbal.liq2_grc, dynbal.ice2_grc, bounds,
                       mask_nolakec, mask_lakec, col, lun,
                       waterstate, waterdiagnostic, lakestate)

    delta_ice = zeros(Float64, bounds.endg)
    if get_for_testing_zero_dynbal_fluxes(ctl)
        for g in bounds.begg:bounds.endg
            delta_liq[g] = 0.0
            delta_ice[g] = 0.0
        end
    else
        for g in bounds.begg:bounds.endg
            delta_liq[g] = dynbal.liq2_grc[g] - dynbal.liq1_grc[g]
            delta_ice[g] = dynbal.ice2_grc[g] - dynbal.ice1_grc[g]
        end
    end

    # Dribbler pass-through: the dynbal flux equals the gridcell content delta.
    for g in bounds.begg:bounds.endg
        dynbal.qflx_liq_dynbal_grc[g] = delta_liq[g]
        dynbal.qflx_ice_dynbal_grc[g] = delta_ice[g]
    end

    return nothing
end

# ---------------------------------------------------------------------------
# Public: dyn_hwcontent_final
# ---------------------------------------------------------------------------

"""
    dyn_hwcontent_final!(dynbal, bounds, mask_nolakec, mask_lakec, col, lun,
                         urbanparams, soilstate, water, temperature, lakestate, ctl)

Compute grid cell-level heat and water content AFTER land cover change, and
from (after - before) produce the dynbal water flux (`qflx_liq_dynbal_grc` /
`qflx_ice_dynbal_grc`) and energy flux (`eflx_dynbal_grc`). Honors
`get_for_testing_zero_dynbal_fluxes(ctl)` (when true, all dynbal fluxes are
zeroed).

Should be called AFTER all subgrid weight updates this time step.

Ported from `dyn_hwcontent_final` in `dynConsBiogeophysMod.F90`.
"""
function dyn_hwcontent_final!(dynbal::DynConsBiogeophysState,
                              bounds::BoundsType,
                              mask_nolakec::AbstractVector{Bool},
                              mask_lakec::AbstractVector{Bool},
                              col::ColumnData, lun::LandunitData,
                              urbanparams::UrbanParamsData,
                              soilstate::SoilStateData,
                              waterstatebulk::WaterStateBulkData,
                              waterdiagnosticbulk::WaterDiagnosticBulkData,
                              temperature::TemperatureData,
                              lakestate::LakeStateData,
                              ctl::DynSubgridControl)
    delta_liq_bulk = zeros(Float64, bounds.endg)

    # Water content (bulk water) -> dynbal water fluxes + delta_liq
    dyn_water_content_final!(delta_liq_bulk, dynbal, bounds,
                             mask_nolakec, mask_lakec, col, lun,
                             waterstatebulk.ws, waterdiagnosticbulk, lakestate, ctl)

    # Heat content (after)
    dyn_heat_content!(dynbal.heat2_grc, dynbal.liquid_water_temp2_grc, bounds,
                      mask_nolakec, mask_lakec, col, lun, urbanparams, soilstate,
                      temperature, waterstatebulk, waterdiagnosticbulk, lakestate)

    delta_heat = zeros(Float64, bounds.endg)
    if get_for_testing_zero_dynbal_fluxes(ctl)
        for g in bounds.begg:bounds.endg
            delta_heat[g] = 0.0
        end
    else
        for g in bounds.begg:bounds.endg
            delta_heat[g] = dynbal.heat2_grc[g] - dynbal.heat1_grc[g]
        end
    end

    # Adjust delta_heat for the implicit heat carried by the dynbal liquid flux.
    bounds_grc = bounds.begg:bounds.endg
    adjust_delta_heat_for_delta_liq!(bounds_grc, delta_liq_bulk,
                                     dynbal.liquid_water_temp1_grc,
                                     dynbal.liquid_water_temp2_grc,
                                     delta_heat)

    # Dribbler pass-through: the dynbal energy flux equals the heat delta.
    for g in bounds.begg:bounds.endg
        dynbal.eflx_dynbal_grc[g] = delta_heat[g]
    end

    return nothing
end

# ---------------------------------------------------------------------------
# Private: set_glacier_baselines
# ---------------------------------------------------------------------------

"""
    set_glacier_baselines!(baselines_col, vals_col, bounds, mask_icec, col, lun, grc)

Compute start-of-run baseline values for each glacier (ice) column, for some
dynbal term. The baseline for a glacier column is the value in that glacier
column minus the natural-vegetation-landunit average value on its gridcell
(representing the missing soil-under-glacier stock).

`vals_col` must be set for at least the natural-veg and glacier columns.
`baselines_col` is set for all columns in `mask_icec`.

Ported from `set_glacier_baselines` in `dynConsBiogeophysMod.F90`.
"""
function set_glacier_baselines!(baselines_col::AbstractVector{<:Real},
                                vals_col::AbstractVector{<:Real},
                                bounds::BoundsType,
                                mask_icec::AbstractVector{Bool},
                                col::ColumnData, lun::LandunitData,
                                grc::GridcellData)
    # Average column values to landunit (natural veg), including inactive points
    vals_lun = zeros(Float64, bounds.endl)
    c2l_1d!(vals_lun, collect(Float64, vals_col), bounds, "urbanf", col, lun;
            include_inactive = true)

    for c in eachindex(mask_icec)
        mask_icec[c] || continue
        g = col.gridcell[c]

        # Start with the value in this glacier column.
        baselines_col[c] = vals_col[c]

        # Subtract the natural-veg landunit value on this gridcell (missing stock),
        # if such a landunit exists.
        l_natveg = grc.landunit_indices[ISTSOIL, g]
        if l_natveg != ISPVAL
            baselines_col[c] = baselines_col[c] - vals_lun[l_natveg]
        end
    end

    return nothing
end

# ---------------------------------------------------------------------------
# Private: dyn_water_content_set_baselines
# ---------------------------------------------------------------------------

"""
    dyn_water_content_set_baselines!(bounds, mask_natveg_and_glc, mask_icec,
                                     mask_lakec, col, lun, grc, waterstate,
                                     lakestate; reset_all_baselines,
                                     reset_lake_baselines)

Set start-of-run baseline values for water content (`dynbal_baseline_liq_col`,
`dynbal_baseline_ice_col`) for glacier (ice) columns (when
`reset_all_baselines`) and lake columns (when either flag is set).

Ported from `dyn_water_content_set_baselines` in `dynConsBiogeophysMod.F90`.
"""
function dyn_water_content_set_baselines!(bounds::BoundsType,
                                          mask_natveg_and_glc::AbstractVector{Bool},
                                          mask_icec::AbstractVector{Bool},
                                          mask_lakec::AbstractVector{Bool},
                                          col::ColumnData, lun::LandunitData,
                                          grc::GridcellData,
                                          waterstate::WaterStateData,
                                          lakestate::LakeStateData;
                                          reset_all_baselines::Bool,
                                          reset_lake_baselines::Bool)
    nc = bounds.endc

    if reset_all_baselines
        soil_liquid_mass_col = zeros(Float64, nc)
        soil_ice_mass_col    = zeros(Float64, nc)

        accumulate_soil_liq_ice_mass_non_lake!(mask_natveg_and_glc, col, waterstate,
                                               soil_liquid_mass_col, soil_ice_mass_col)

        set_glacier_baselines!(waterstate.dynbal_baseline_liq_col,
                               soil_liquid_mass_col, bounds, mask_icec, col, lun, grc)
        set_glacier_baselines!(waterstate.dynbal_baseline_ice_col,
                               soil_ice_mass_col, bounds, mask_icec, col, lun, grc)
    end

    if reset_all_baselines || reset_lake_baselines
        lake_liquid_mass_col = zeros(Float64, nc)
        lake_ice_mass_col    = zeros(Float64, nc)

        # Total water volume of the lake column (tracer_ratio = 1 for bulk water)
        accumulate_liq_ice_mass_lake!(mask_lakec, col, lakestate, 1.0,
                                      lake_liquid_mass_col, lake_ice_mass_col)

        for c in eachindex(mask_lakec)
            mask_lakec[c] || continue
            waterstate.dynbal_baseline_liq_col[c] = lake_liquid_mass_col[c]
            waterstate.dynbal_baseline_ice_col[c] = lake_ice_mass_col[c]
        end
    end

    return nothing
end

# ---------------------------------------------------------------------------
# Private: dyn_heat_content_set_baselines
# ---------------------------------------------------------------------------

"""
    dyn_heat_content_set_baselines!(bounds, mask_natveg_and_glc, mask_icec,
                                    mask_lakec, col, lun, grc, urbanparams,
                                    soilstate, lakestate, waterstatebulk,
                                    temperature; reset_all_baselines,
                                    reset_lake_baselines)

Set start-of-run baseline values for heat content
(`dynbal_baseline_heat_col`) for glacier columns (when `reset_all_baselines`)
and lake columns (when either flag is set).

Ported from `dyn_heat_content_set_baselines` in `dynConsBiogeophysMod.F90`.
"""
function dyn_heat_content_set_baselines!(bounds::BoundsType,
                                         mask_natveg_and_glc::AbstractVector{Bool},
                                         mask_icec::AbstractVector{Bool},
                                         mask_lakec::AbstractVector{Bool},
                                         col::ColumnData, lun::LandunitData,
                                         grc::GridcellData,
                                         urbanparams::UrbanParamsData,
                                         soilstate::SoilStateData,
                                         lakestate::LakeStateData,
                                         waterstatebulk::WaterStateBulkData,
                                         waterdiagnosticbulk::WaterDiagnosticBulkData,
                                         temperature::TemperatureData;
                                         reset_all_baselines::Bool,
                                         reset_lake_baselines::Bool)
    nc = bounds.endc

    if reset_all_baselines
        soil_heat_col        = zeros(Float64, nc)
        # heat_liquid / cv_liquid are needed only for the helper interface
        soil_heat_liquid_col = zeros(Float64, nc)
        soil_cv_liquid_col   = zeros(Float64, nc)

        accumulate_soil_heat_non_lake!(mask_natveg_and_glc, col, lun, urbanparams,
                                       soilstate, temperature, waterstatebulk,
                                       soil_heat_col, soil_heat_liquid_col,
                                       soil_cv_liquid_col)

        set_glacier_baselines!(temperature.dynbal_baseline_heat_col,
                               soil_heat_col, bounds, mask_icec, col, lun, grc)
    end

    if reset_all_baselines || reset_lake_baselines
        lake_heat_col = zeros(Float64, nc)
        accumulate_heat_lake!(mask_lakec, col, temperature, lakestate, lake_heat_col)

        for c in eachindex(mask_lakec)
            mask_lakec[c] || continue
            temperature.dynbal_baseline_heat_col[c] = lake_heat_col[c]
        end
    end

    return nothing
end

# ---------------------------------------------------------------------------
# Public: dyn_hwcontent_set_baselines
# ---------------------------------------------------------------------------

"""
    dyn_hwcontent_set_baselines!(bounds, mask_icec, mask_lakec, col, lun, grc,
                                 urbanparams, soilstate, lakestate,
                                 waterstatebulk, waterdiagnosticbulk, temperature;
                                 reset_all_baselines, reset_lake_baselines)

Set start-of-run baseline values for heat and water content in glacier and lake
columns. These baselines are subtracted from each column's total heat/water
calculation each time step, reducing the fictitious dynbal conservation fluxes.

The natural-veg + glacier (istsoil + istice) column mask used for soil
accumulation includes inactive points (so an inactive point that later becomes
active has a baseline). This mask is built internally from `lun.itype` /
`col.landunit`.

Ported from `dyn_hwcontent_set_baselines` in `dynConsBiogeophysMod.F90`.
"""
function dyn_hwcontent_set_baselines!(bounds::BoundsType,
                                      mask_icec::AbstractVector{Bool},
                                      mask_lakec::AbstractVector{Bool},
                                      col::ColumnData, lun::LandunitData,
                                      grc::GridcellData,
                                      urbanparams::UrbanParamsData,
                                      soilstate::SoilStateData,
                                      lakestate::LakeStateData,
                                      waterstatebulk::WaterStateBulkData,
                                      waterdiagnosticbulk::WaterDiagnosticBulkData,
                                      temperature::TemperatureData;
                                      reset_all_baselines::Bool,
                                      reset_lake_baselines::Bool)
    # Natural veg + glacier column mask (istsoil + istice), including inactive
    # points (mirrors col_filter_from_ltypes with include_inactive = .true.).
    mask_natveg_and_glc = falses(bounds.endc)
    for c in bounds.begc:bounds.endc
        lt = lun.itype[col.landunit[c]]
        if lt == ISTSOIL || lt == ISTICE
            mask_natveg_and_glc[c] = true
        end
    end

    # Water baselines (bulk water only)
    dyn_water_content_set_baselines!(bounds, mask_natveg_and_glc, mask_icec,
                                     mask_lakec, col, lun, grc,
                                     waterstatebulk.ws, lakestate;
                                     reset_all_baselines = reset_all_baselines,
                                     reset_lake_baselines = reset_lake_baselines)

    # Heat baselines
    dyn_heat_content_set_baselines!(bounds, mask_natveg_and_glc, mask_icec,
                                    mask_lakec, col, lun, grc, urbanparams,
                                    soilstate, lakestate, waterstatebulk,
                                    waterdiagnosticbulk, temperature;
                                    reset_all_baselines = reset_all_baselines,
                                    reset_lake_baselines = reset_lake_baselines)

    return nothing
end
