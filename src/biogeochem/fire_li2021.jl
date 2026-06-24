# ==========================================================================
# Ported from: src/biogeochem/CNFireLi2021Mod.F90
# Fire dynamics module for Li et al. (2021) version.
#
# Incremental refinement of Li2016. Deltas vs Li2016:
#   * Root wetness uses the Li2021 method (CNFire_calc_fire_root_wetness_Li2021).
#   * btran_col accumulation rescales btran2 by per-PFT (rswf_min, rswf_max).
#   * Fire-occurrence moisture term drops the bt_min/bt_max clamp: the btran
#     factor is simply (1 - btran_col/wtlf)^0.5.
#   * Deforestation-fire climate term `cli` uses a revised dtrotr formula, and
#     the resulting farea includes the fuel factor fb (fb*cli).
#
# Public functions:
#   need_lightning_and_popdens_li2021  — Returns true
#   cnfire_area_li2021!                — Compute column-level burned area
#   cnfire_fluxes_li2021!              — Fire C/N fluxes (delegates to base)
#
# NOTE: Non-default fire method — straightforward mask-based host loops.
# ==========================================================================

"""
    need_lightning_and_popdens_li2021()

Returns `true`. Ported from `need_lightning_and_popdens` in `CNFireLi2021Mod.F90`.
"""
need_lightning_and_popdens_li2021() = true

"""
    cnfire_area_li2021!(...)

Compute column-level burned-area fraction for the Li et al. (2021) fire model.
Ported from `CNFireArea` in `CNFireLi2021Mod.F90`.
"""
function cnfire_area_li2021!(
    fire_li2014::CNFireLi2014Data,
    pftcon_li2014::PftConFireLi2014,
    fire_data::CNFireBaseData,
    cnfire_const::CNFireConstData,
    cnfire_params::CNFireParams,
    pftcon::PftConFireBase,
    mask_soilc::AbstractVector{Bool},
    mask_soilp::AbstractVector{Bool},
    mask_exposedveg::AbstractVector{Bool},
    mask_noexposedveg::AbstractVector{Bool},
    bounds_c::UnitRange{Int},
    bounds_p::UnitRange{Int},
    patch::PatchData,
    col::ColumnData,
    grc::GridcellData,
    soilstate::SoilStateData,
    h2osoi_vol_col::AbstractMatrix{<:Real},
    cnveg_state::CNVegStateData,
    cnveg_cs::CNVegCarbonStateData,
    decomp_cascade_con::DecompCascadeConData,
    totlitc_col::AbstractVector{<:Real},
    decomp_cpools_vr_col::AbstractArray{<:Real,3},
    t_soi17cm_col::AbstractVector{<:Real};
    forc_rh_grc::AbstractVector{<:Real} = Float64[],
    forc_wind_grc::AbstractVector{<:Real} = Float64[],
    forc_t_col::AbstractVector{<:Real} = Float64[],
    forc_rain_col::AbstractVector{<:Real} = Float64[],
    forc_snow_col::AbstractVector{<:Real} = Float64[],
    prec60_patch::AbstractVector{<:Real} = Float64[],
    prec10_patch::AbstractVector{<:Real} = Float64[],
    rh30_patch::AbstractVector{<:Real} = Float64[],
    fsat_col::AbstractVector{<:Real} = Float64[],
    wf2_col::AbstractVector{<:Real} = Float64[],
    dt::Real = 1800.0,
    dayspyr::Real = 365.0,
    kmo::Int = 1,
    kda::Int = 1,
    mcsec::Int = 0,
    nstep::Int = 1,
    nlevgrnd::Int = 10,
    nlevdecomp::Int = 1,
    ndecomp_pools::Int = 7,
    nc4_grass::Int = 14,
    nc3crop::Int = 15,
    ndllf_evr_tmp_tree::Int = 1,
    nbrdlf_evr_trp_tree::Int = 4,
    nbrdlf_dcd_trp_tree::Int = 6,
    nbrdlf_evr_shrub::Int = 9,
    noveg::Int = 0,
    transient_landcover::Bool = false,
    spinup_state::Int = 0,
    spinup_factor_deadwood::Real = SPINUP_FACTOR_DEADWOOD_DEFAULT,
    soil_suction_fn::Function = default_soil_suction,
    _ignored...,  # absorb extra kwargs forwarded by the cnfire_area! dispatcher
)
    totlitc          = totlitc_col
    decomp_cpools_vr = decomp_cpools_vr_col
    tsoi17           = t_soi17cm_col
    lfuel            = cnfire_const.lfuel
    ufuel            = cnfire_const.ufuel
    rh_hgh           = cnfire_const.rh_hgh
    rh_low           = cnfire_const.rh_low
    cli_scale        = cnfire_const.cli_scale
    cropfire_a1      = cnfire_const.cropfire_a1
    non_boreal_peatfire_c    = cnfire_const.non_boreal_peatfire_c
    pot_hmn_ign_counts_alpha = cnfire_const.pot_hmn_ign_counts_alpha
    boreal_peatfire_c        = cnfire_const.boreal_peatfire_c
    max_rh30_affecting_fuel  = cnfire_const.max_rh30_affecting_fuel
    defo_fire_precip_thresh_bet  = cnfire_const.defo_fire_precip_thresh_bet
    defo_fire_precip_thresh_bdt  = cnfire_const.defo_fire_precip_thresh_bdt
    borpeat_fire_soilmoist_denom = cnfire_const.borpeat_fire_soilmoist_denom
    nonborpeat_fire_precip_denom = cnfire_const.nonborpeat_fire_precip_denom
    occur_hi_gdp_tree            = cnfire_const.occur_hi_gdp_tree
    g0                           = cnfire_const.g0
    prh30                        = cnfire_params.prh30
    ignition_efficiency          = cnfire_params.ignition_efficiency

    fsr_pft  = pftcon_li2014.fsr_pft
    fd_pft   = pftcon_li2014.fd_pft
    rswf_min = pftcon.rswf_min
    rswf_max = pftcon.rswf_max

    btran2   = fire_data.btran2_patch
    gdp_lf   = fire_li2014.gdp_lf_col
    peatf_lf = fire_li2014.peatf_lf_col
    abm_lf   = fire_li2014.abm_lf_col
    forc_hdm = fire_li2014.forc_hdm
    forc_lnfm = fire_li2014.forc_lnfm
    is_cwd   = decomp_cascade_con.is_cwd
    spinup_factor = decomp_cascade_con.spinup_factor

    forc_rh   = forc_rh_grc
    forc_wind = forc_wind_grc
    forc_t    = forc_t_col
    forc_rain = forc_rain_col
    forc_snow = forc_snow_col
    fsat      = fsat_col
    wf2       = wf2_col

    dwt_smoothed = cnveg_state.dwt_smoothed_patch
    cropf_col    = cnveg_state.cropf_col
    baf_crop     = cnveg_state.baf_crop_col
    baf_peatf    = cnveg_state.baf_peatf_col
    burndate     = cnveg_state.burndate_patch
    fbac         = cnveg_state.fbac_col
    fbac1        = cnveg_state.fbac1_col
    farea_burned = cnveg_state.farea_burned_col
    nfire        = cnveg_state.nfire_col
    fsr_col      = cnveg_state.fsr_col
    fd_col       = cnveg_state.fd_col
    lgdp_col     = cnveg_state.lgdp_col
    lgdp1_col    = cnveg_state.lgdp1_col
    lpop_col     = cnveg_state.lpop_col
    lfwt         = cnveg_state.lfwt_col
    trotr1_col   = cnveg_state.trotr1_col
    trotr2_col   = cnveg_state.trotr2_col
    dtrotr_col   = cnveg_state.dtrotr_col
    lfc          = cnveg_state.lfc_col
    wtlf         = cnveg_state.wtlf_col

    totvegc       = cnveg_cs.totvegc_col
    rootc_col     = cnveg_cs.rootc_col
    leafc_col     = cnveg_cs.leafc_col
    deadstemc_col = cnveg_cs.deadstemc_col
    fuelc         = cnveg_cs.fuelc_col
    fuelc_crop    = cnveg_cs.fuelc_crop_col

    deadcrootc         = cnveg_cs.deadcrootc_patch
    deadcrootc_storage = cnveg_cs.deadcrootc_storage_patch
    deadcrootc_xfer    = cnveg_cs.deadcrootc_xfer_patch
    deadstemc          = cnveg_cs.deadstemc_patch
    frootc             = cnveg_cs.frootc_patch
    frootc_storage     = cnveg_cs.frootc_storage_patch
    frootc_xfer        = cnveg_cs.frootc_xfer_patch
    livecrootc         = cnveg_cs.livecrootc_patch
    livecrootc_storage = cnveg_cs.livecrootc_storage_patch
    livecrootc_xfer    = cnveg_cs.livecrootc_xfer_patch
    leafc              = cnveg_cs.leafc_patch
    leafc_storage      = cnveg_cs.leafc_storage_patch
    leafc_xfer         = cnveg_cs.leafc_xfer_patch

    secsphr  = SECSPHR
    secspday = SECSPDAY

    nc = last(bounds_c)
    FT = eltype(totvegc)
    btran_col  = zeros(FT, nc)
    prec60_col = zeros(FT, nc)
    prec10_col = zeros(FT, nc)
    rh30_col   = zeros(FT, nc)

    p2c!(prec10_col, prec10_patch, patch, mask_soilc, bounds_c, bounds_p)
    p2c!(prec60_col, prec60_patch, patch, mask_soilc, bounds_c, bounds_p)
    p2c!(rh30_col,   rh30_patch,   patch, mask_soilc, bounds_c, bounds_p)
    p2c!(leafc_col,  leafc,        patch, mask_soilc, bounds_c, bounds_p)
    p2c!(deadstemc_col, deadstemc, patch, mask_soilc, bounds_c, bounds_p)

    if nstep == 0
        for c in bounds_c
            mask_soilc[c] || continue
            farea_burned[c] = 0.0
            baf_crop[c]     = 0.0
            baf_peatf[c]    = 0.0
            fbac[c]         = 0.0
            fbac1[c]        = 0.0
            cropf_col[c]    = 0.0
        end
        return nothing
    end

    for c in bounds_c
        mask_soilc[c] || continue
        cropf_col[c] = 0.0
        lfwt[c]      = 0.0
    end
    for p in bounds_p
        mask_soilp[p] || continue
        c = patch.column[p]
        if patch.itype[p] > nc4_grass
            cropf_col[c] += patch.wtcol[p]
        end
        if patch.itype[p] >= ndllf_evr_tmp_tree && patch.itype[p] <= nc4_grass
            lfwt[c] += patch.wtcol[p]
        end
    end

    for c in bounds_c
        mask_soilc[c] || continue
        fuelc_crop[c] = 0.0
    end
    for p in bounds_p
        mask_soilp[p] || continue
        c = patch.column[p]
        if patch.itype[p] > nc4_grass && patch.wtcol[p] > 0.0 && leafc_col[c] > 0.0
            fuelc_crop[c] += (leafc[p] + leafc_storage[p] + leafc_xfer[p]) *
                             patch.wtcol[p] / cropf_col[c] +
                             totlitc[c] * leafc[p] / leafc_col[c] *
                             patch.wtcol[p] / cropf_col[c]
        end
    end

    for c in bounds_c
        mask_soilc[c] || continue
        fsr_col[c]    = 0.0
        fd_col[c]     = 0.0
        rootc_col[c]  = 0.0
        lgdp_col[c]   = 0.0
        lgdp1_col[c]  = 0.0
        lpop_col[c]   = 0.0
        btran_col[c]  = 0.0
        wtlf[c]       = 0.0
        trotr1_col[c] = 0.0
        trotr2_col[c] = 0.0
        if transient_landcover
            dtrotr_col[c] = 0.0
        end
    end

    # NOTE(Li2021): root wetness uses the Li2021 method.
    cnfire_calc_fire_root_wetness_li2021!(
        fire_data, mask_exposedveg, mask_noexposedveg, bounds_p,
        patch, soilstate, h2osoi_vol_col, nlevgrnd)

    for p in bounds_p
        mask_exposedveg[p] || continue
        c = patch.column[p]
        if patch.itype[p] < nc3crop && cropf_col[c] < 1.0
            # NOTE(Li2021): rescale root wetness by per-PFT (rswf_min, rswf_max).
            it = patch.itype[p]
            btran_col[c] += max(0.0, min(1.0,
                (btran2[p] - rswf_min[it + 1]) / (rswf_max[it + 1] - rswf_min[it + 1]))) *
                patch.wtcol[p]
            wtlf[c] += patch.wtcol[p]
        end
    end

    for p in bounds_p
        mask_soilp[p] || continue
        c = patch.column[p]
        g = col.gridcell[c]
        if patch.itype[p] < nc3crop && cropf_col[c] < 1.0
            if patch.itype[p] == nbrdlf_evr_trp_tree && patch.wtcol[p] > 0.0
                trotr1_col[c] += patch.wtcol[p]
            end
            if patch.itype[p] == nbrdlf_dcd_trp_tree && patch.wtcol[p] > 0.0
                trotr2_col[c] += patch.wtcol[p]
            end
            if transient_landcover
                if patch.itype[p] == nbrdlf_evr_trp_tree || patch.itype[p] == nbrdlf_dcd_trp_tree
                    if dwt_smoothed[p] < 0.0
                        dtrotr_col[c] += -dwt_smoothed[p]
                    end
                end
            end
            rootc_col[c] += (frootc[p] + frootc_storage[p] + frootc_xfer[p] +
                deadcrootc[p] * spinup_factor_deadwood +
                deadcrootc_storage[p] + deadcrootc_xfer[p] +
                livecrootc[p] + livecrootc_storage[p] +
                livecrootc_xfer[p]) * patch.wtcol[p]

            fsr_col[c] += fsr_pft[patch.itype[p] + 1] * patch.wtcol[p] / (1.0 - cropf_col[c])

            hdmlf = forc_hdm[g]
            if hdmlf > 0.1
                if patch.itype[p] != noveg
                    if patch.itype[p] >= nbrdlf_evr_shrub
                        lgdp_col[c]  += (0.1 + 0.9 * exp(-1.0 * RPI * (gdp_lf[c] / 8.0)^0.5)) *
                                        patch.wtcol[p] / (1.0 - cropf_col[c])
                        lgdp1_col[c] += (0.2 + 0.8 * exp(-1.0 * RPI * (gdp_lf[c] / 7.0))) *
                                        patch.wtcol[p] / (1.0 - cropf_col[c])
                        lpop_col[c]  += (0.2 + 0.8 * exp(-1.0 * RPI * (hdmlf / 450.0)^0.5)) *
                                        patch.wtcol[p] / (1.0 - cropf_col[c])
                    else
                        if gdp_lf[c] > 20.0
                            lgdp_col[c]  += occur_hi_gdp_tree * patch.wtcol[p] / (1.0 - cropf_col[c])
                            lgdp1_col[c] += 0.62 * patch.wtcol[p] / (1.0 - cropf_col[c])
                        else
                            if gdp_lf[c] > 8.0
                                lgdp_col[c]  += 0.79 * patch.wtcol[p] / (1.0 - cropf_col[c])
                                lgdp1_col[c] += 0.83 * patch.wtcol[p] / (1.0 - cropf_col[c])
                            else
                                lgdp_col[c]  += patch.wtcol[p] / (1.0 - cropf_col[c])
                                lgdp1_col[c] += patch.wtcol[p] / (1.0 - cropf_col[c])
                            end
                        end
                        lpop_col[c] += (0.4 + 0.6 * exp(-1.0 * RPI * (hdmlf / 125.0))) *
                                       patch.wtcol[p] / (1.0 - cropf_col[c])
                    end
                end
            else
                lgdp_col[c]  += patch.wtcol[p] / (1.0 - cropf_col[c])
                lgdp1_col[c] += patch.wtcol[p] / (1.0 - cropf_col[c])
                lpop_col[c]  += patch.wtcol[p] / (1.0 - cropf_col[c])
            end

            fd_col[c] += fd_pft[patch.itype[p]] * patch.wtcol[p] * secsphr / (1.0 - cropf_col[c])
        end
    end

    if transient_landcover
        for c in bounds_c
            mask_soilc[c] || continue
            if dtrotr_col[c] > 0.0
                if kmo == 1 && kda == 1 && mcsec == 0
                    lfc[c] = 0.0
                end
                if kmo == 1 && kda == 1 && mcsec == dt
                    lfc[c] = dtrotr_col[c] * dayspyr * secspday / dt
                end
            else
                lfc[c] = 0.0
            end
        end
    end

    for c in bounds_c
        mask_soilc[c] || continue
        baf_crop[c] = 0.0
    end
    for p in bounds_p
        mask_soilp[p] || continue
        if kmo == 1 && kda == 1 && mcsec == 0
            burndate[p] = 10000
        end
    end
    for p in bounds_p
        mask_soilp[p] || continue
        c = patch.column[p]
        g = col.gridcell[c]
        if forc_t[c] >= TFRZ && patch.itype[p] > nc4_grass &&
           kmo == abm_lf[c] && burndate[p] >= 999 && patch.wtcol[p] > 0.0
            hdmlf = forc_hdm[g]
            fhd  = 0.04 + 0.96 * exp(-1.0 * RPI * (hdmlf / 350.0)^0.5)
            fgdp = 0.01 + 0.99 * exp(-1.0 * RPI * (gdp_lf[c] / 10.0))
            fb   = max(0.0, min(1.0, (fuelc_crop[c] - lfuel) / (ufuel - lfuel)))
            baf_crop[c] += cropfire_a1 / secsphr * fhd * fgdp * patch.wtcol[p]
            if fb * fhd * fgdp * patch.wtcol[p] > 0.0
                burndate[p] = kda
            end
        end
    end

    for c in bounds_c
        mask_soilc[c] || continue
        g = col.gridcell[c]
        if grc.latdeg[g] < cnfire_const.borealat
            baf_peatf[c] = non_boreal_peatfire_c / secsphr *
                max(0.0, min(1.0, (4.0 - prec60_col[c] * secspday / nonborpeat_fire_precip_denom) / 4.0))^2 *
                peatf_lf[c] * (1.0 - fsat[c])
        else
            baf_peatf[c] = boreal_peatfire_c / secsphr *
                exp(-RPI * (max(wf2[c], 0.0) / borpeat_fire_soilmoist_denom)) *
                max(0.0, min(1.0, (tsoi17[c] - TFRZ) / 10.0)) * peatf_lf[c] *
                (1.0 - fsat[c])
        end
    end

    i_cwd = 0
    for l in 1:ndecomp_pools
        if is_cwd[l]
            i_cwd = l
        end
    end
    dzsoi = dzsoi_decomp[]

    for c in bounds_c
        mask_soilc[c] || continue
        g = col.gridcell[c]
        hdmlf = forc_hdm[g]
        nfire[c] = 0.0
        if cropf_col[c] < 1.0
            fuelc[c] = totlitc[c] + totvegc[c] - rootc_col[c] - fuelc_crop[c] * cropf_col[c]
            if spinup_state == 2
                fuelc[c] += (spinup_factor_deadwood - 1.0) * deadstemc_col[c]
                for j in 1:nlevdecomp
                    fuelc[c] += decomp_cpools_vr[c, j, i_cwd] * dzsoi[j] * spinup_factor[i_cwd] *
                                get_spinup_latitude_term(grc.latdeg[col.gridcell[c]])
                end
            else
                for j in 1:nlevdecomp
                    fuelc[c] += decomp_cpools_vr[c, j, i_cwd] * dzsoi[j]
                end
            end
            fuelc[c] /= (1.0 - cropf_col[c])
            fb = max(0.0, min(1.0, (fuelc[c] - lfuel) / (ufuel - lfuel)))
            if trotr1_col[c] + trotr2_col[c] <= 0.6
                afuel = min(1.0, max(0.0, (fuelc[c] - 2500.0) / (5000.0 - 2500.0)))
                arh   = 1.0 - max(0.0, min(1.0, (forc_rh[g] - rh_low) / (rh_hgh - rh_low)))
                arh30 = 1.0 - max(prh30, min(1.0, rh30_col[c] / max_rh30_affecting_fuel))
                if forc_rh[g] < rh_hgh && wtlf[c] > 0.0 && tsoi17[c] > TFRZ
                    # NOTE(Li2021): btran term has no bt_min/bt_max clamp.
                    fire_m = ((afuel * arh30 + (1.0 - afuel) * arh)^1.5) *
                             ((1.0 - btran_col[c] / wtlf[c])^0.5)
                else
                    fire_m = 0.0
                end
                lh = pot_hmn_ign_counts_alpha * 6.8 * hdmlf^0.43 / 30.0 / 24.0
                fs = 1.0 - (0.01 + 0.98 * exp(-0.025 * hdmlf))
                ig = (lh + forc_lnfm[g] /
                      (5.16 + 2.16 * cos(RPI / 180.0 * 3 * min(60.0, abs(grc.latdeg[g])))) *
                      ignition_efficiency) * (1.0 - fs) * (1.0 - cropf_col[c])
                nfire[c] = ig / secsphr * fb * fire_m * lgdp_col[c]
                Lb_lf = 1.0 + 10.0 * (1.0 - exp(-0.06 * forc_wind[g]))
                spread_m = fire_m^0.5
                farea_burned[c] = min(1.0, (g0 * spread_m * fsr_col[c] *
                    fd_col[c] / 1000.0)^2 * lgdp1_col[c] * lpop_col[c] *
                    nfire[c] * RPI * Lb_lf + baf_crop[c] + baf_peatf[c])
            else
                farea_burned[c] = min(1.0, baf_crop[c] + baf_peatf[c])
            end

            if transient_landcover
                if trotr1_col[c] + trotr2_col[c] > 0.6
                    if (kmo == 1 && kda == 1 && mcsec == 0) || dtrotr_col[c] <= 0.0
                        fbac1[c]        = 0.0
                        farea_burned[c] = baf_crop[c] + baf_peatf[c]
                    else
                        cri = (defo_fire_precip_thresh_bet * trotr1_col[c] +
                               defo_fire_precip_thresh_bdt * trotr2_col[c]) /
                              (trotr1_col[c] + trotr2_col[c])
                        # NOTE(Li2021): revised deforestation climate term.
                        cli = (max(0.0, min(1.0, (cri - prec60_col[c] * secspday) / cri))^0.5) *
                              (max(0.0, min(1.0, (cri - prec10_col[c] * secspday) / cri))^0.5) *
                              (15.0 * min(0.0016, dtrotr_col[c] / dt * dayspyr * secspday) + 0.009) *
                              max(0.0, min(1.0, (0.25 - (forc_rain[c] + forc_snow[c]) * secsphr) / 0.25))
                        farea_burned[c] = fb * cli * (cli_scale / secspday) + baf_crop[c] + baf_peatf[c]
                        fbac1[c] = max(0.0, fb * cli * (cli_scale / secspday) - 2.0 * lfc[c] / dt)
                    end
                    fbac[c] = fbac1[c] + baf_crop[c] + baf_peatf[c]
                else
                    fbac[c] = farea_burned[c]
                end
            end
        else
            farea_burned[c] = min(1.0, baf_crop[c] + baf_peatf[c])
        end
    end

    return nothing
end

"""
    cnfire_fluxes_li2021!(...)

Fire effects routine for the Li et al. (2021) fire model. Inherits the base
`CNFireFluxes` (litter 0.5 / CWD 0.25); delegates to `cnfire_fluxes_li2014!`.
"""
function cnfire_fluxes_li2021!(args...; kwargs...)
    return cnfire_fluxes_li2014!(args...; kwargs...)
end
