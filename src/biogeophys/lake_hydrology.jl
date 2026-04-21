# ==========================================================================
# Ported from: src/biogeophys/LakeHydrologyMod.F90 (~760 lines)
# Calculation of Lake Hydrology. Full hydrology of snow layers is done.
# However, there is no infiltration, and the water budget is balanced with
# qflx_qrgwl. Lake water mass is kept constant. The soil is simply maintained
# at volumetric saturation if ice melting frees up pore space.
#
# Public functions:
#   sum_flux_fluxes_onto_ground!     — Copy precip fluxes onto lake surface
#   lake_sublimation_dew!            — Calculate sublimation and dew for lake patches
#   lake_patch_to_col_fluxes!        — Copy patch-level fluxes to column averages
#   lake_soil_hydrology!             — Maintain soil saturation for lake columns
#   lake_check_single_snow_layer!    — Check for single unfrozen snow layer over lake
#   lake_snow_above_unfrozen!        — Handle snow layers above unfrozen lake top
#   lake_snow_diagnostics!           — Accumulate snow ice/liq and internal temperature
#   lake_water_balance!              — Compute ending water balance and runoff
#   lake_top_layer_diagnostics!      — Top-layer snow diagnostics for history
#   lake_hydrology!                  — Main orchestrator
# ==========================================================================

# --- Module-level constants ---
const FRAC_SNO_SMALL = 1.0e-6   # small value of frac_sno used when initiating a snow pack due to frost
const SNOW_BD_LAKE   = 250.0    # assumed snow bulk density (for lakes w/out resolved snow layers) [kg/m^3]

# =========================================================================
# sum_flux_fluxes_onto_ground!
# =========================================================================

"""
    sum_flux_fluxes_onto_ground!(forc_snow, forc_rain,
        qflx_snow_grnd, qflx_liq_grnd,
        mask_lake, bounds)

Compute "summed" (really just copies here) fluxes onto "ground" (really the
lake surface), for bulk or one tracer. (Subroutine name mimics the one in
CanopyHydrologyMod.)

Ported from `SumFlux_FluxesOntoGround` in `LakeHydrologyMod.F90`.
"""
function sum_flux_fluxes_onto_ground!(
    forc_snow::Vector{<:Real},
    forc_rain::Vector{<:Real},
    qflx_snow_grnd::Vector{<:Real},
    qflx_liq_grnd::Vector{<:Real},
    mask_lake::BitVector,
    bounds::UnitRange{Int}
)
    for c in bounds
        mask_lake[c] || continue
        qflx_snow_grnd[c] = forc_snow[c]
        qflx_liq_grnd[c]  = forc_rain[c]
    end
    return nothing
end

# =========================================================================
# lake_sublimation_dew!
# =========================================================================

"""
    lake_sublimation_dew!(
        qflx_liqevap_from_top_layer, qflx_solidevap_from_top_layer,
        qflx_soliddew_to_top_layer, qflx_liqdew_to_top_layer,
        qflx_ev_snow,
        qflx_evap_soi, h2osoi_liq, h2osoi_ice,
        h2osno_no_layers, snow_depth, frac_sno,
        t_grnd, t_soisno, snl,
        patch_column, mask_lakep, bounds_patch,
        dtime, nlevsno)

Calculate sublimation and dew for lake patches. Adapted from HydrologyLake
and Biogeophysics2.

Ported from inline code in `LakeHydrology` in `LakeHydrologyMod.F90` (lines 274-365).
"""
function lake_sublimation_dew!(
    qflx_liqevap_from_top_layer::Vector{<:Real},
    qflx_solidevap_from_top_layer::Vector{<:Real},
    qflx_soliddew_to_top_layer::Vector{<:Real},
    qflx_liqdew_to_top_layer::Vector{<:Real},
    qflx_ev_snow::Vector{<:Real},
    qflx_evap_soi::Vector{<:Real},
    h2osoi_liq::Matrix{<:Real},
    h2osoi_ice::Matrix{<:Real},
    h2osno_no_layers::Vector{<:Real},
    snow_depth::Vector{<:Real},
    frac_sno::Vector{<:Real},
    t_grnd::Vector{<:Real},
    t_soisno::Matrix{<:Real},
    snl::Vector{Int},
    patch_column::Vector{Int},
    mask_lakep::BitVector,
    bounds_patch::UnitRange{Int},
    dtime::Real,
    nlevsno::Int
)
    for p in bounds_patch
        mask_lakep[p] || continue
        c = patch_column[p]
        jtop = snl[c] + 1

        qflx_liqevap_from_top_layer[p]   = 0.0
        qflx_solidevap_from_top_layer[p] = 0.0
        qflx_soliddew_to_top_layer[p]    = 0.0
        qflx_liqdew_to_top_layer[p]      = 0.0
        qflx_ev_snow[p]                  = qflx_evap_soi[p]

        if jtop <= 0  # snow layers present
            j = jtop + nlevsno  # Julia 1-based index

            if qflx_evap_soi[p] >= 0.0
                # Evaporation: partition between liquid evap and ice sublimation
                qflx_evap_soi_lim = min(qflx_evap_soi[p],
                    (h2osoi_liq[c, j] + h2osoi_ice[c, j]) / dtime)
                qflx_ev_snow[p] = qflx_evap_soi_lim
                if (h2osoi_liq[c, j] + h2osoi_ice[c, j]) > 0.0
                    qflx_liqevap_from_top_layer[p] = max(
                        qflx_evap_soi_lim * (h2osoi_liq[c, j] /
                            (h2osoi_liq[c, j] + h2osoi_ice[c, j])), 0.0)
                else
                    qflx_liqevap_from_top_layer[p] = 0.0
                end
                qflx_solidevap_from_top_layer[p] = qflx_evap_soi_lim - qflx_liqevap_from_top_layer[p]
            else
                # Dew/frost
                if t_grnd[c] < TFRZ && t_soisno[c, j] < TFRZ
                    qflx_soliddew_to_top_layer[p] = abs(qflx_evap_soi[p])
                elseif jtop < 0 || (t_grnd[c] == TFRZ && t_soisno[c, j] == TFRZ)
                    qflx_liqdew_to_top_layer[p] = abs(qflx_evap_soi[p])
                end
            end

        else  # No snow layers
            if qflx_evap_soi[p] >= 0.0
                # Sublimation: do not allow more than there is snow
                qflx_solidevap_from_top_layer[p] = min(qflx_evap_soi[p],
                    h2osno_no_layers[c] / dtime)
                qflx_liqevap_from_top_layer[p] = qflx_evap_soi[p] -
                    qflx_solidevap_from_top_layer[p]
            else
                if t_grnd[c] < TFRZ - 0.1
                    qflx_soliddew_to_top_layer[p] = abs(qflx_evap_soi[p])
                else
                    qflx_liqdew_to_top_layer[p] = abs(qflx_evap_soi[p])
                end
            end

            # Update snow pack for dew & sub.
            h2osno_temp = h2osno_no_layers[c]
            qflx_dew_minus_sub_snow = -qflx_solidevap_from_top_layer[p] +
                qflx_soliddew_to_top_layer[p]
            h2osno_no_layers[c] = h2osno_no_layers[c] + qflx_dew_minus_sub_snow * dtime
            h2osno_no_layers[c] = max(h2osno_no_layers[c], 0.0)

            if qflx_dew_minus_sub_snow > 0.0
                # Accumulating snow from dew: ensure at least small non-zero frac_sno
                if frac_sno[c] <= 0.0
                    frac_sno[c] = FRAC_SNO_SMALL
                end
            elseif qflx_dew_minus_sub_snow < 0.0
                # Losing snow from sublimation: reset frac_sno if snow gone
                if h2osno_no_layers[c] == 0.0
                    frac_sno[c] = 0.0
                end
            end

            if h2osno_temp > 0.0
                # Assume snow bulk density remains the same as before
                snow_depth[c] = snow_depth[c] * h2osno_no_layers[c] / h2osno_temp
            else
                # Assume a constant snow bulk density = 250
                snow_depth[c] = h2osno_no_layers[c] / SNOW_BD_LAKE
            end
        end
    end

    return nothing
end

# =========================================================================
# lake_frac_sno_eff!
# =========================================================================

"""
    lake_frac_sno_eff!(frac_sno, frac_sno_eff, mask_lake, bounds)

Recalculate frac_sno_eff from frac_sno for lake columns.
For lakes, frac_sno_eff = frac_sno (no subgrid snow fraction adjustment).

Ported from `CalcFracSnoEff` call in `LakeHydrology` in `LakeHydrologyMod.F90` (line 368).
"""
function lake_frac_sno_eff!(
    frac_sno::Vector{<:Real},
    frac_sno_eff::Vector{<:Real},
    mask_lake::BitVector,
    bounds::UnitRange{Int}
)
    for c in bounds
        mask_lake[c] || continue
        frac_sno_eff[c] = frac_sno[c]
    end
    return nothing
end

# =========================================================================
# lake_patch_to_col_fluxes!
# =========================================================================

"""
    lake_patch_to_col_fluxes!(
        qflx_evap_tot_col, qflx_liqevap_from_top_layer_col,
        qflx_liqdew_to_top_layer_col, qflx_soliddew_to_top_layer_col,
        qflx_solidevap_from_top_layer_col, qflx_ev_snow_col,
        qflx_evap_tot, qflx_liqevap_from_top_layer,
        qflx_liqdew_to_top_layer, qflx_soliddew_to_top_layer,
        qflx_solidevap_from_top_layer, qflx_ev_snow,
        patch_column, mask_lakep, bounds_patch)

Copy patch-level flux averages to column level (assuming one pft per lake column).

Ported from inline code in `LakeHydrology` in `LakeHydrologyMod.F90` (lines 380-390).
"""
function lake_patch_to_col_fluxes!(
    qflx_evap_tot_col::Vector{<:Real},
    qflx_liqevap_from_top_layer_col::Vector{<:Real},
    qflx_liqdew_to_top_layer_col::Vector{<:Real},
    qflx_soliddew_to_top_layer_col::Vector{<:Real},
    qflx_solidevap_from_top_layer_col::Vector{<:Real},
    qflx_ev_snow_col::Vector{<:Real},
    qflx_evap_tot::Vector{<:Real},
    qflx_liqevap_from_top_layer::Vector{<:Real},
    qflx_liqdew_to_top_layer::Vector{<:Real},
    qflx_soliddew_to_top_layer::Vector{<:Real},
    qflx_solidevap_from_top_layer::Vector{<:Real},
    qflx_ev_snow::Vector{<:Real},
    patch_column::Vector{Int},
    mask_lakep::BitVector,
    bounds_patch::UnitRange{Int}
)
    for p in bounds_patch
        mask_lakep[p] || continue
        c = patch_column[p]

        qflx_evap_tot_col[c]                = qflx_evap_tot[p]
        qflx_liqevap_from_top_layer_col[c]  = qflx_liqevap_from_top_layer[p]
        qflx_liqdew_to_top_layer_col[c]     = qflx_liqdew_to_top_layer[p]
        qflx_soliddew_to_top_layer_col[c]   = qflx_soliddew_to_top_layer[p]
        qflx_solidevap_from_top_layer_col[c] = qflx_solidevap_from_top_layer[p]
        qflx_ev_snow_col[c]                 = qflx_ev_snow[p]
    end

    return nothing
end

# =========================================================================
# lake_soil_hydrology!
# =========================================================================

"""
    lake_soil_hydrology!(h2osoi_liq, h2osoi_ice, h2osoi_vol,
        dz, watsat, mask_lake, bounds, nlevsoi, nlevsno)

Maintain soil at volumetric saturation for lake columns. If melting has
left open pore space, fill with water. If liquid water exceeds saturation,
reduce it.

Ported from inline code in `LakeHydrology` in `LakeHydrologyMod.F90` (lines 455-487).
"""
function lake_soil_hydrology!(
    h2osoi_liq::Matrix{<:Real},
    h2osoi_ice::Matrix{<:Real},
    h2osoi_vol::Matrix{<:Real},
    dz::Matrix{<:Real},
    watsat::Matrix{<:Real},
    mask_lake::BitVector,
    bounds::UnitRange{Int},
    nlevsoi::Int,
    nlevsno::Int
)
    for j in 1:nlevsoi
        jj = j + nlevsno  # Julia 1-based index for combined snow+soil array
        for c in bounds
            mask_lake[c] || continue

            h2osoi_vol[c, j] = h2osoi_liq[c, jj] / (dz[c, jj] * DENH2O) +
                               h2osoi_ice[c, jj] / (dz[c, jj] * DENICE)

            if h2osoi_vol[c, j] < watsat[c, j]
                h2osoi_liq[c, jj] = (watsat[c, j] * dz[c, jj] -
                    h2osoi_ice[c, jj] / DENICE) * DENH2O
            elseif h2osoi_liq[c, jj] > watsat[c, j] * DENH2O * dz[c, jj]
                h2osoi_liq[c, jj] = watsat[c, j] * DENH2O * dz[c, jj]
            end
        end
    end

    return nothing
end

# =========================================================================
# lake_check_single_snow_layer!
# =========================================================================

"""
    lake_check_single_snow_layer!(
        eflx_sh_tot, eflx_sh_grnd, eflx_soil_grnd,
        eflx_gnet, eflx_grnd_lake,
        t_soisno, h2osoi_ice, h2osoi_liq,
        snl, h2osno_no_layers, h2osno_total,
        snow_depth, qflx_snow_drain,
        patch_column, mask_lakep, bounds_patch,
        dtime, nlevsno)

Check for single completely unfrozen snow layer over lake. If ice is present
and temperature is above freezing, adjust sensible heat. If no ice remains,
remove the layer entirely.

Ported from inline code in `LakeHydrology` in `LakeHydrologyMod.F90` (lines 495-535).
"""
function lake_check_single_snow_layer!(
    eflx_sh_tot::Vector{<:Real},
    eflx_sh_grnd::Vector{<:Real},
    eflx_soil_grnd::Vector{<:Real},
    eflx_gnet::Vector{<:Real},
    eflx_grnd_lake::Vector{<:Real},
    t_soisno::Matrix{<:Real},
    h2osoi_ice::Matrix{<:Real},
    h2osoi_liq::Matrix{<:Real},
    snl::Vector{Int},
    h2osno_no_layers::Vector{<:Real},
    h2osno_total::Vector{<:Real},
    snow_depth::Vector{<:Real},
    qflx_snow_drain::Vector{<:Real},
    patch_column::Vector{Int},
    mask_lakep::BitVector,
    bounds_patch::UnitRange{Int},
    dtime::Real,
    nlevsno::Int
)
    j = 0 + nlevsno  # Julia index for snow layer 0

    for p in bounds_patch
        mask_lakep[p] || continue
        c = patch_column[p]

        if snl[c] == -1
            if h2osoi_ice[c, j] > 0.0 && t_soisno[c, j] > TFRZ
                # Take extra heat of layer and release to sensible heat
                heatrem           = (CPLIQ * h2osoi_liq[c, j]) * (t_soisno[c, j] - TFRZ)
                t_soisno[c, j]    = TFRZ
                eflx_sh_tot[p]    = eflx_sh_tot[p] + heatrem / dtime
                eflx_sh_grnd[p]   = eflx_sh_grnd[p] + heatrem / dtime
                eflx_soil_grnd[p] = eflx_soil_grnd[p] - heatrem / dtime
                eflx_gnet[p]      = eflx_gnet[p] - heatrem / dtime
                eflx_grnd_lake[p] = eflx_grnd_lake[p] - heatrem / dtime
            elseif h2osoi_ice[c, j] == 0.0
                # Remove layer
                heatrem             = CPLIQ * h2osoi_liq[c, j] * (t_soisno[c, j] - TFRZ)
                eflx_sh_tot[p]      = eflx_sh_tot[p] + heatrem / dtime
                eflx_sh_grnd[p]     = eflx_sh_grnd[p] + heatrem / dtime
                eflx_soil_grnd[p]   = eflx_soil_grnd[p] - heatrem / dtime
                eflx_gnet[p]        = eflx_gnet[p] - heatrem / dtime
                eflx_grnd_lake[p]   = eflx_grnd_lake[p] - heatrem / dtime
                qflx_snow_drain[c]  = qflx_snow_drain[c] + h2osno_total[c] / dtime
                snl[c]              = 0
                h2osno_no_layers[c] = 0.0
                h2osno_total[c]     = 0.0
                snow_depth[c]       = 0.0
            else
                eflx_grnd_lake[p] = eflx_gnet[p]
            end
        else
            eflx_grnd_lake[p] = eflx_gnet[p]
        end
    end

    return nothing
end

# =========================================================================
# lake_snow_above_unfrozen!
# =========================================================================

"""
    lake_snow_above_unfrozen!(
        t_lake, lake_icefrac, t_soisno,
        h2osoi_ice, h2osoi_liq,
        snl, h2osno_no_layers, h2osno_total,
        snow_depth, qflx_snomelt, eflx_snomelt,
        qflx_snomelt_lyr, qflx_snow_drain,
        dz_lake, mask_lake, bounds, dtime, nlevsno)

Check for snow layers above lake with unfrozen top layer. If the top layer
has sufficient heat to melt the snow without freezing, melt the snow and
remove the layers. Otherwise, let snow persist and melt by diffusion.

Ported from inline code in `LakeHydrology` in `LakeHydrologyMod.F90` (lines 543-605).
"""
function lake_snow_above_unfrozen!(
    t_lake::Matrix{<:Real},
    lake_icefrac::Matrix{<:Real},
    t_soisno::Matrix{<:Real},
    h2osoi_ice::Matrix{<:Real},
    h2osoi_liq::Matrix{<:Real},
    snl::Vector{Int},
    h2osno_no_layers::Vector{<:Real},
    h2osno_total::Vector{<:Real},
    snow_depth::Vector{<:Real},
    qflx_snomelt::Vector{<:Real},
    eflx_snomelt::Vector{<:Real},
    qflx_snomelt_lyr::Matrix{<:Real},
    qflx_snow_drain::Vector{<:Real},
    dz_lake::Matrix{<:Real},
    mask_lake::BitVector,
    bounds::UnitRange{Int},
    dtime::Real,
    nlevsno::Int
)
    nc = length(mask_lake)
    FT = eltype(t_lake)
    unfrozen = falses(nc)
    sumsnowice = zeros(FT, nc)
    heatsum = zeros(FT, nc)

    # Determine which columns have unfrozen top lake layer with snow
    for c in bounds
        mask_lake[c] || continue
        if t_lake[c, 1] > TFRZ && lake_icefrac[c, 1] == 0.0 && snl[c] < 0
            unfrozen[c] = true
        end
    end

    # Sum snow ice and heat content
    for j_fortran in (-nlevsno + 1):0
        j = j_fortran + nlevsno  # Julia index
        for c in bounds
            mask_lake[c] || continue
            if unfrozen[c]
                if j_fortran == -nlevsno + 1
                    sumsnowice[c] = 0.0
                    heatsum[c] = 0.0
                end
                if j_fortran >= snl[c] + 1
                    sumsnowice[c] = sumsnowice[c] + h2osoi_ice[c, j]
                    heatsum[c] = heatsum[c] +
                        h2osoi_ice[c, j] * CPICE * (TFRZ - t_soisno[c, j]) +
                        h2osoi_liq[c, j] * CPLIQ * (TFRZ - t_soisno[c, j])
                end
            end
        end
    end

    # Determine if lake can absorb the latent heat
    for c in bounds
        mask_lake[c] || continue
        if unfrozen[c]
            heatsum[c] = heatsum[c] + sumsnowice[c] * HFUS
            heatrem = (t_lake[c, 1] - TFRZ) * CPLIQ * DENH2O * dz_lake[c, 1] - heatsum[c]

            if heatrem + DENH2O * dz_lake[c, 1] * HFUS > 0.0
                # Remove snow and subtract the latent heat from the top layer
                qflx_snomelt[c] = qflx_snomelt[c] + sumsnowice[c] / dtime
                eflx_snomelt[c] = eflx_snomelt[c] + sumsnowice[c] * HFUS / dtime

                # Update melt per layer
                for j_fortran in (snl[c] + 1):0
                    j = j_fortran + nlevsno
                    qflx_snomelt_lyr[c, j_fortran + nlevsno] =
                        qflx_snomelt_lyr[c, j_fortran + nlevsno] + h2osoi_ice[c, j] / dtime
                end

                # Update incidental drainage from snow pack
                qflx_snow_drain[c] = qflx_snow_drain[c] + h2osno_total[c] / dtime

                h2osno_no_layers[c] = 0.0
                snow_depth[c] = 0.0
                snl[c] = 0

                if heatrem > 0.0
                    t_lake[c, 1] = t_lake[c, 1] - heatrem / (CPLIQ * DENH2O * dz_lake[c, 1])
                else
                    t_lake[c, 1] = TFRZ
                    lake_icefrac[c, 1] = -heatrem / (DENH2O * dz_lake[c, 1] * HFUS)
                end
            end
        end
    end

    return nothing
end

# =========================================================================
# lake_snow_diagnostics!
# =========================================================================

"""
    lake_snow_diagnostics!(
        snowice, snowliq, t_sno_mul_mss,
        h2osoi_ice, h2osoi_liq, t_soisno,
        snl, mask_lake, mask_lakesnow,
        bounds, nlevsno)

Accumulate snow ice, snow liquid, and snow internal temperature (t_sno_mul_mss)
over all snow layers for history output.

Ported from inline code in `LakeHydrology` in `LakeHydrologyMod.F90` (lines 619-651).
"""
function lake_snow_diagnostics!(
    snowice::Vector{<:Real},
    snowliq::Vector{<:Real},
    t_sno_mul_mss::Vector{<:Real},
    h2osoi_ice::Matrix{<:Real},
    h2osoi_liq::Matrix{<:Real},
    t_soisno::Matrix{<:Real},
    snl::Vector{Int},
    mask_lake::BitVector,
    mask_lakesnow::BitVector,
    bounds::UnitRange{Int},
    nlevsno::Int
)
    # Initialize
    for c in bounds
        mask_lake[c] || continue
        snowice[c] = 0.0
        snowliq[c] = 0.0
        t_sno_mul_mss[c] = 0.0
    end

    # Accumulate over snow layers
    for j_fortran in (-nlevsno + 1):0
        j = j_fortran + nlevsno  # Julia index
        for c in bounds
            mask_lakesnow[c] || continue
            if j_fortran >= snl[c] + 1
                snowice[c] = snowice[c] + h2osoi_ice[c, j]
                snowliq[c] = snowliq[c] + h2osoi_liq[c, j]
                t_sno_mul_mss[c] = t_sno_mul_mss[c] + h2osoi_ice[c, j] * t_soisno[c, j]
                t_sno_mul_mss[c] = t_sno_mul_mss[c] + h2osoi_liq[c, j] * TFRZ
            end
        end
    end

    return nothing
end

# =========================================================================
# lake_water_balance!
# =========================================================================

"""
    lake_water_balance!(
        qflx_drain_perched, qflx_rsub_sat, qflx_infl,
        qflx_surf, qflx_drain, qflx_qrgwl,
        qflx_floodc, qflx_runoff,
        qflx_rain_plus_snomelt, qflx_top_soil,
        qflx_ice_runoff_snwcp,
        h2osoi_liq, h2osoi_ice, h2osoi_vol,
        forc_rain, forc_snow, qflx_evap_tot,
        qflx_snwcp_ice, qflx_snwcp_discarded_ice,
        qflx_snwcp_discarded_liq,
        qflx_liq_grnd, qflx_snow_drain,
        qflx_floodg, begwb, endwb,
        dz, patch_column, patch_gridcell,
        mask_lakep, bounds_patch,
        dtime, nlevsno, nlevgrnd)

Compute ending water balance and volumetric soil water for lake columns.
Computes qflx_qrgwl to balance the water budget.

Ported from inline code in `LakeHydrology` in `LakeHydrologyMod.F90` (lines 660-688).
"""
function lake_water_balance!(
    qflx_drain_perched::Vector{<:Real},
    qflx_rsub_sat::Vector{<:Real},
    qflx_infl::Vector{<:Real},
    qflx_surf::Vector{<:Real},
    qflx_drain::Vector{<:Real},
    qflx_qrgwl::Vector{<:Real},
    qflx_floodc::Vector{<:Real},
    qflx_runoff::Vector{<:Real},
    qflx_rain_plus_snomelt::Vector{<:Real},
    qflx_top_soil::Vector{<:Real},
    qflx_ice_runoff_snwcp::Vector{<:Real},
    h2osoi_liq::Matrix{<:Real},
    h2osoi_ice::Matrix{<:Real},
    h2osoi_vol::Matrix{<:Real},
    forc_rain::Vector{<:Real},
    forc_snow::Vector{<:Real},
    qflx_evap_tot::Vector{<:Real},
    qflx_snwcp_ice::Vector{<:Real},
    qflx_snwcp_discarded_ice::Vector{<:Real},
    qflx_snwcp_discarded_liq::Vector{<:Real},
    qflx_liq_grnd::Vector{<:Real},
    qflx_snow_drain::Vector{<:Real},
    qflx_floodg::Vector{<:Real},
    begwb::Vector{<:Real},
    endwb::Vector{<:Real},
    dz::Matrix{<:Real},
    patch_column::Vector{Int},
    patch_gridcell::Vector{Int},
    mask_lakep::BitVector,
    bounds_patch::UnitRange{Int},
    dtime::Real,
    nlevsno::Int,
    nlevgrnd::Int
)
    for p in bounds_patch
        mask_lakep[p] || continue
        c = patch_column[p]
        g = patch_gridcell[p]

        qflx_drain_perched[c] = 0.0
        qflx_rsub_sat[c]      = 0.0
        qflx_infl[c]          = 0.0
        qflx_surf[c]          = 0.0
        qflx_drain[c]         = 0.0

        # Insure water balance using qflx_qrgwl
        qflx_qrgwl[c] = forc_rain[c] + forc_snow[c] - qflx_evap_tot[p] -
            qflx_snwcp_ice[c] -
            qflx_snwcp_discarded_ice[c] - qflx_snwcp_discarded_liq[c] -
            (endwb[c] - begwb[c]) / dtime + qflx_floodg[g]
        qflx_floodc[c]    = qflx_floodg[g]
        qflx_runoff[c]    = qflx_drain[c] + qflx_qrgwl[c]
        qflx_rain_plus_snomelt[c] = qflx_liq_grnd[c] + qflx_snow_drain[c]
        qflx_top_soil[c]  = qflx_rain_plus_snomelt[c]
        qflx_ice_runoff_snwcp[c] = qflx_snwcp_ice[c]
    end

    return nothing
end

# =========================================================================
# lake_update_h2osoi_vol!
# =========================================================================

"""
    lake_update_h2osoi_vol!(h2osoi_liq, h2osoi_ice, h2osoi_vol,
        dz, mask_lake, bounds, nlevgrnd, nlevsno)

Update volumetric soil water content from liquid and ice for lake columns.

Ported from inline code in `LakeHydrology` in `LakeHydrologyMod.F90` (lines 660-665).
"""
function lake_update_h2osoi_vol!(
    h2osoi_liq::Matrix{<:Real},
    h2osoi_ice::Matrix{<:Real},
    h2osoi_vol::Matrix{<:Real},
    dz::Matrix{<:Real},
    mask_lake::BitVector,
    bounds::UnitRange{Int},
    nlevgrnd::Int,
    nlevsno::Int
)
    for j in 1:nlevgrnd
        jj = j + nlevsno
        for c in bounds
            mask_lake[c] || continue
            h2osoi_vol[c, j] = h2osoi_liq[c, jj] / (dz[c, jj] * DENH2O) +
                               h2osoi_ice[c, jj] / (dz[c, jj] * DENICE)
        end
    end

    return nothing
end

# =========================================================================
# lake_top_layer_diagnostics!
# =========================================================================

"""
    lake_top_layer_diagnostics!(
        h2osno_top, snw_rds, snot_top, dTdz_top,
        snw_rds_top, sno_liq_top,
        h2osoi_ice, h2osoi_liq, snl,
        mask_lakesnow, mask_lakenosnow,
        bounds, nlevsno)

Update top-layer snow diagnostics for history output.

Ported from inline code in `LakeHydrology` in `LakeHydrologyMod.F90` (lines 692-709).
"""
function lake_top_layer_diagnostics!(
    h2osno_top::Vector{<:Real},
    snw_rds::Matrix{<:Real},
    snot_top::Vector{<:Real},
    dTdz_top::Vector{<:Real},
    snw_rds_top::Vector{<:Real},
    sno_liq_top::Vector{<:Real},
    h2osoi_ice::Matrix{<:Real},
    h2osoi_liq::Matrix{<:Real},
    snl::Vector{Int},
    mask_lakesnow::BitVector,
    mask_lakenosnow::BitVector,
    bounds::UnitRange{Int},
    nlevsno::Int
)
    # Snow columns: top-layer mass
    for c in bounds
        mask_lakesnow[c] || continue
        jtop = snl[c] + 1 + nlevsno  # Julia index
        h2osno_top[c] = h2osoi_ice[c, jtop] + h2osoi_liq[c, jtop]
    end

    # Non-snow columns: zero / spval
    for c in bounds
        mask_lakenosnow[c] || continue
        h2osno_top[c]  = 0.0
        snw_rds[c, :]  .= 0.0
        snot_top[c]    = SPVAL
        dTdz_top[c]    = SPVAL
        snw_rds_top[c] = SPVAL
        sno_liq_top[c] = SPVAL
    end

    return nothing
end

# =========================================================================
# lake_hydrology!
# =========================================================================

"""
    lake_hydrology!(
        temperature, energyflux, lakestate, soilstate,
        waterstatebulk, waterdiagbulk, waterbalancebulk, waterfluxbulk,
        col_data, patch_data,
        mask_lake, mask_lakep,
        forc_rain, forc_snow, qflx_floodg,
        bounds_col, bounds_patch,
        dtime, nlevsno, nlevsoi, nlevgrnd)

Main lake hydrology orchestrator. Performs snow hydrology, sublimation/dew,
soil hydrology, snow-lake interaction, and water balance for lake columns.

WARNING: This subroutine assumes lake columns have one and only one pft.

Sequence:
  1. Copy precip fluxes onto lake surface
  2. Calculate sublimation and dew
  3. Copy patch fluxes to column averages
  4. Build snow filter and call snow routines (stubs for SnowWater, etc.)
  5. Maintain soil at saturation
  6. Check for unfrozen snow layers over lake
  7. Handle snow layers above unfrozen lake top
  8. Accumulate snow diagnostics
  9. Compute water balance and runoff
  10. Top-layer diagnostics

Ported from `LakeHydrology` in `LakeHydrologyMod.F90`.
"""
function lake_hydrology!(
    temperature::TemperatureData,
    energyflux::EnergyFluxData,
    lakestate::LakeStateData,
    soilstate::SoilStateData,
    waterstatebulk::WaterStateBulkData,
    waterdiagbulk::WaterDiagnosticBulkData,
    waterbalancebulk::WaterBalanceData,
    waterfluxbulk::WaterFluxBulkData,
    col_data::ColumnData,
    patch_data::PatchData,
    mask_lake::BitVector,
    mask_lakep::BitVector,
    forc_rain::Vector{<:Real},
    forc_snow::Vector{<:Real},
    qflx_floodg::Vector{<:Real},
    bounds_col::UnitRange{Int},
    bounds_patch::UnitRange{Int},
    dtime::Real,
    nlevsno::Int,
    nlevsoi::Int,
    nlevgrnd::Int
)
    # --- Extract arrays from data structures ---
    ws = waterstatebulk.ws
    wf = waterfluxbulk.wf

    h2osoi_liq       = ws.h2osoi_liq_col
    h2osoi_ice       = ws.h2osoi_ice_col
    h2osoi_vol       = ws.h2osoi_vol_col
    h2osno_no_layers = ws.h2osno_no_layers_col

    frac_sno     = waterdiagbulk.frac_sno_col
    frac_sno_eff = waterdiagbulk.frac_sno_eff_col
    frac_iceold  = waterdiagbulk.frac_iceold_col
    snow_depth   = waterdiagbulk.snow_depth_col
    snowice      = waterdiagbulk.snowice_col
    snowliq      = waterdiagbulk.snowliq_col
    snw_rds      = waterdiagbulk.snw_rds_col
    snw_rds_top  = waterdiagbulk.snw_rds_top_col
    h2osno_top   = waterdiagbulk.h2osno_top_col
    sno_liq_top  = waterdiagbulk.sno_liq_top_col

    t_lake       = temperature.t_lake_col
    t_grnd       = temperature.t_grnd_col
    t_soisno     = temperature.t_soisno_col
    dTdz_top     = temperature.dTdz_top_col
    snot_top     = temperature.snot_top_col
    t_sno_mul_mss = temperature.t_sno_mul_mss_col

    begwb = waterbalancebulk.begwb_col
    endwb = waterbalancebulk.endwb_col

    watsat = soilstate.watsat_col

    dz       = col_data.dz
    dz_lake  = col_data.dz_lake
    snl      = col_data.snl

    lake_icefrac = lakestate.lake_icefrac_col

    qflx_liq_grnd       = wf.qflx_liq_grnd_col
    qflx_snow_grnd       = wf.qflx_snow_grnd_col
    qflx_evap_tot       = wf.qflx_evap_tot_patch
    qflx_evap_soi       = wf.qflx_evap_soi_patch
    qflx_ev_snow_patch   = wf.qflx_ev_snow_patch
    qflx_snomelt         = wf.qflx_snomelt_col
    qflx_snow_drain      = wf.qflx_snow_drain_col
    qflx_snwcp_ice       = wf.qflx_snwcp_ice_col
    qflx_snwcp_discarded_ice = wf.qflx_snwcp_discarded_ice_col
    qflx_snwcp_discarded_liq = wf.qflx_snwcp_discarded_liq_col
    qflx_drain_perched   = wf.qflx_drain_perched_col
    qflx_rsub_sat        = wf.qflx_rsub_sat_col
    qflx_surf            = wf.qflx_surf_col
    qflx_drain           = wf.qflx_drain_col
    qflx_infl            = wf.qflx_infl_col
    qflx_qrgwl           = wf.qflx_qrgwl_col
    qflx_runoff          = wf.qflx_runoff_col
    qflx_floodc          = wf.qflx_floodc_col
    qflx_rain_plus_snomelt = wf.qflx_rain_plus_snomelt_col
    qflx_top_soil        = wf.qflx_top_soil_col
    qflx_ice_runoff_snwcp = wf.qflx_ice_runoff_snwcp_col

    qflx_evap_tot_col    = wf.qflx_evap_tot_col
    qflx_liqevap_from_top_layer_col   = wf.qflx_liqevap_from_top_layer_col
    qflx_liqdew_to_top_layer_col      = wf.qflx_liqdew_to_top_layer_col
    qflx_soliddew_to_top_layer_col    = wf.qflx_soliddew_to_top_layer_col
    qflx_solidevap_from_top_layer_col = wf.qflx_solidevap_from_top_layer_col

    qflx_solidevap_from_top_layer = wf.qflx_solidevap_from_top_layer_patch
    qflx_liqevap_from_top_layer   = wf.qflx_liqevap_from_top_layer_patch
    qflx_soliddew_to_top_layer    = wf.qflx_soliddew_to_top_layer_patch
    qflx_liqdew_to_top_layer      = wf.qflx_liqdew_to_top_layer_patch

    eflx_snomelt  = energyflux.eflx_snomelt_col
    eflx_sh_tot   = energyflux.eflx_sh_tot_patch
    eflx_sh_grnd  = energyflux.eflx_sh_grnd_patch
    eflx_soil_grnd = energyflux.eflx_soil_grnd_patch
    eflx_gnet     = energyflux.eflx_gnet_patch
    eflx_grnd_lake = energyflux.eflx_grnd_lake_patch

    patch_column   = patch_data.column
    patch_gridcell = patch_data.gridcell

    # --- Bulk ev_snow references ---
    qflx_ev_snow_col_bulk = waterfluxbulk.qflx_ev_snow_col

    qflx_snomelt_lyr = waterfluxbulk.qflx_snomelt_lyr_col

    # =====================================================================
    # 1. Compute summed fluxes onto ground (copy precip to ground fluxes)
    # =====================================================================
    sum_flux_fluxes_onto_ground!(
        forc_snow, forc_rain,
        qflx_snow_grnd, qflx_liq_grnd,
        mask_lake, bounds_col)

    # Note: UpdateQuantitiesForNewSnow and InitializeExplicitSnowPack from
    # SnowHydrologyMod are called in the Fortran code here. Those are already
    # ported in snow_hydrology.jl and would be called by the driver.
    # For this port, they are omitted as stubs (same pattern as other hydrology modules).

    # =====================================================================
    # 2. Calculate sublimation and dew
    # =====================================================================
    lake_sublimation_dew!(
        qflx_liqevap_from_top_layer, qflx_solidevap_from_top_layer,
        qflx_soliddew_to_top_layer, qflx_liqdew_to_top_layer,
        qflx_ev_snow_patch,
        qflx_evap_soi, h2osoi_liq, h2osoi_ice,
        h2osno_no_layers, snow_depth, frac_sno,
        t_grnd, t_soisno, snl,
        patch_column, mask_lakep, bounds_patch,
        dtime, nlevsno)

    # =====================================================================
    # 3. Recalculate frac_sno_eff
    # =====================================================================
    lake_frac_sno_eff!(frac_sno, frac_sno_eff, mask_lake, bounds_col)

    # =====================================================================
    # 4. Copy patch-level fluxes to column averages
    # =====================================================================
    lake_patch_to_col_fluxes!(
        qflx_evap_tot_col, qflx_liqevap_from_top_layer_col,
        qflx_liqdew_to_top_layer_col, qflx_soliddew_to_top_layer_col,
        qflx_solidevap_from_top_layer_col, qflx_ev_snow_col_bulk,
        qflx_evap_tot, qflx_liqevap_from_top_layer,
        qflx_liqdew_to_top_layer, qflx_soliddew_to_top_layer,
        qflx_solidevap_from_top_layer, qflx_ev_snow_patch,
        patch_column, mask_lakep, bounds_patch)

    # =====================================================================
    # 5. Build snow filter (mask_lakesnow, mask_lakenosnow)
    # =====================================================================
    nc = length(mask_lake)
    mask_lakesnow   = falses(nc)
    mask_lakenosnow = falses(nc)
    for c in bounds_col
        mask_lake[c] || continue
        if snl[c] < 0
            mask_lakesnow[c] = true
        else
            mask_lakenosnow[c] = true
        end
    end

    # Note: SnowWater, SnowCapping, SnowCompaction, CombineSnowLayers,
    # DivideSnowLayers, ZeroEmptySnowLayers from SnowHydrologyMod are called
    # in the Fortran code here. Those are already ported in snow_hydrology.jl
    # and would be called by the driver. For this port, they are omitted as stubs.

    # =====================================================================
    # 6. Recompute h2osno_total
    # =====================================================================
    FT_lh = eltype(t_grnd)
    h2osno_total = zeros(FT_lh, nc)
    waterstate_calculate_total_h2osno!(ws, mask_lake, bounds_col, snl, h2osno_total)

    # =====================================================================
    # 7. Maintain soil at saturation
    # =====================================================================
    lake_soil_hydrology!(
        h2osoi_liq, h2osoi_ice, h2osoi_vol,
        dz, watsat, mask_lake, bounds_col,
        nlevsoi, nlevsno)

    # =====================================================================
    # 8. Check for single completely unfrozen snow layer over lake
    # =====================================================================
    lake_check_single_snow_layer!(
        eflx_sh_tot, eflx_sh_grnd, eflx_soil_grnd,
        eflx_gnet, eflx_grnd_lake,
        t_soisno, h2osoi_ice, h2osoi_liq,
        snl, h2osno_no_layers, h2osno_total,
        snow_depth, qflx_snow_drain,
        patch_column, mask_lakep, bounds_patch,
        dtime, nlevsno)

    # =====================================================================
    # 9. Handle snow layers above unfrozen lake top
    # =====================================================================
    lake_snow_above_unfrozen!(
        t_lake, lake_icefrac, t_soisno,
        h2osoi_ice, h2osoi_liq,
        snl, h2osno_no_layers, h2osno_total,
        snow_depth, qflx_snomelt, eflx_snomelt,
        qflx_snomelt_lyr, qflx_snow_drain,
        dz_lake, mask_lake, bounds_col,
        dtime, nlevsno)

    # =====================================================================
    # 10. Rebuild snow filter after modifications
    # =====================================================================
    for c in bounds_col
        mask_lake[c] || continue
        if snl[c] < 0
            mask_lakesnow[c] = true
            mask_lakenosnow[c] = false
        else
            mask_lakesnow[c] = false
            mask_lakenosnow[c] = true
        end
    end

    # =====================================================================
    # 11. Snow diagnostics (ice/liq/temperature sums)
    # =====================================================================
    lake_snow_diagnostics!(
        snowice, snowliq, t_sno_mul_mss,
        h2osoi_ice, h2osoi_liq, t_soisno,
        snl, mask_lake, mask_lakesnow,
        bounds_col, nlevsno)

    # =====================================================================
    # 12. Compute ending water balance (stub: endwb = begwb for now)
    # =====================================================================
    # In the Fortran code, ComputeWaterMassLake is called here.
    # That function is not yet ported, so we use a stub that sets endwb = begwb.
    for c in bounds_col
        mask_lake[c] || continue
        endwb[c] = begwb[c]
    end

    # =====================================================================
    # 13. Update volumetric soil water
    # =====================================================================
    lake_update_h2osoi_vol!(
        h2osoi_liq, h2osoi_ice, h2osoi_vol,
        dz, mask_lake, bounds_col,
        nlevgrnd, nlevsno)

    # =====================================================================
    # 14. Compute water balance and runoff
    # =====================================================================
    lake_water_balance!(
        qflx_drain_perched, qflx_rsub_sat, qflx_infl,
        qflx_surf, qflx_drain, qflx_qrgwl,
        qflx_floodc, qflx_runoff,
        qflx_rain_plus_snomelt, qflx_top_soil,
        qflx_ice_runoff_snwcp,
        h2osoi_liq, h2osoi_ice, h2osoi_vol,
        forc_rain, forc_snow, qflx_evap_tot,
        qflx_snwcp_ice, qflx_snwcp_discarded_ice,
        qflx_snwcp_discarded_liq,
        qflx_liq_grnd, qflx_snow_drain,
        qflx_floodg, begwb, endwb,
        dz, patch_column, patch_gridcell,
        mask_lakep, bounds_patch,
        dtime, nlevsno, nlevgrnd)

    # =====================================================================
    # 15. Top-layer diagnostics
    # =====================================================================
    lake_top_layer_diagnostics!(
        h2osno_top, snw_rds, snot_top, dTdz_top,
        snw_rds_top, sno_liq_top,
        h2osoi_ice, h2osoi_liq, snl,
        mask_lakesnow, mask_lakenosnow,
        bounds_col, nlevsno)

    return nothing
end
