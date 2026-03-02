# ==========================================================================
# Ported from: src/biogeochem/CNFireLi2014Mod.F90
# Fire dynamics module for Li et al. (2014) version.
# Based on Li et al. (2012a,b; 2013; 2014).
#
# Public types:
#   CNFireLi2014Data    — Li2014-specific fire data (forcing, GDP, peatland)
#   PftConFireLi2014    — PFT-level fire parameters specific to Li2014
#
# Public functions:
#   need_lightning_and_popdens_li2014  — Returns true
#   p2c!                               — Patch to column area-weighted average
#   cnfire_area_li2014!                — Compute column-level burned area
# ==========================================================================

# ---------------------------------------------------------------------------
# Li2014-specific data
# ---------------------------------------------------------------------------

"""
    CNFireLi2014Data

Li2014-specific fire data: population density, lightning frequency,
GDP, peatland fraction, and crop fire timing.
Ported from member data of `cnfire_li2014_type` in `CNFireLi2014Mod.F90`.
"""
Base.@kwdef mutable struct CNFireLi2014Data
    forc_hdm     ::Vector{Float64} = Float64[]  # human population density (#/km2) (gridcell)
    forc_lnfm    ::Vector{Float64} = Float64[]  # lightning frequency (#/km2/hr) (gridcell)
    gdp_lf_col   ::Vector{Float64} = Float64[]  # GDP data (k$/capita) (column)
    peatf_lf_col ::Vector{Float64} = Float64[]  # peatland fraction data (0-1) (column)
    abm_lf_col   ::Vector{Int}     = Int[]      # prescribed crop fire month (1-12) (column)
end

# ---------------------------------------------------------------------------
# PFT-level fire parameters specific to Li2014
# ---------------------------------------------------------------------------

"""
    PftConFireLi2014

PFT-level fire parameters specific to the Li2014 fire module.
Contains fire spread rate and fire duration by PFT.
"""
Base.@kwdef mutable struct PftConFireLi2014
    fsr_pft ::Vector{Float64} = Float64[]  # fire spread rate by PFT (km/hr)
    fd_pft  ::Vector{Float64} = Float64[]  # fire duration by PFT (hr)
end

# ---------------------------------------------------------------------------
# need_lightning_and_popdens_li2014
# ---------------------------------------------------------------------------

"""
    need_lightning_and_popdens_li2014()

Returns `true` — the Li2014 fire model requires lightning and population
density data.
Ported from `need_lightning_and_popdens` in `CNFireLi2014Mod.F90`.
"""
need_lightning_and_popdens_li2014() = true

# ---------------------------------------------------------------------------
# p2c! — patch to column area-weighted average
# ---------------------------------------------------------------------------

"""
    p2c!(col_out, patch_in, patch, mask_soilc, bounds_c, bounds_p)

Patch-to-column area-weighted average.
Computes: col_out[c] = Σ_p (patch_in[p] * patch.wtcol[p]) for all p on column c.
Ported from `p2c` in `subgridAveMod.F90`.
"""
function p2c!(
    col_out::Vector{Float64},
    patch_in::Vector{Float64},
    patch_d::PatchData,
    mask_soilc::BitVector,
    bounds_c::UnitRange{Int},
    bounds_p::UnitRange{Int}
)
    for c in bounds_c
        mask_soilc[c] || continue
        col_out[c] = 0.0
    end
    for p in bounds_p
        c = patch_d.column[p]
        (c in bounds_c && mask_soilc[c]) || continue
        col_out[c] += patch_in[p] * patch_d.wtcol[p]
    end
    return nothing
end

# ---------------------------------------------------------------------------
# cnfire_area_li2014! — Compute column-level burned area
# ---------------------------------------------------------------------------

"""
    cnfire_area_li2014!(fire_li2014, pftcon_li2014, fire_data, cnfire_const,
        cnfire_params, pftcon, mask_soilc, mask_soilp, mask_exposedveg,
        mask_noexposedveg, bounds_c, bounds_p, patch, col, grc,
        soilstate, h2osoi_vol_col, cnveg_state, cnveg_cs,
        decomp_cascade_con, totlitc_col, decomp_cpools_vr_col,
        t_soi17cm_col; kwargs...)

Compute column-level burned area fraction for the Li et al. (2014)
fire model. Includes cropland fire, peatland fire, deforestation fire,
and natural fire components.

Ported from `CNFireArea` in `CNFireLi2014Mod.F90`.
"""
function cnfire_area_li2014!(
    # Li2014-specific data
    fire_li2014::CNFireLi2014Data,
    pftcon_li2014::PftConFireLi2014,
    # Fire base data & constants
    fire_data::CNFireBaseData,
    cnfire_const::CNFireConstData,
    cnfire_params::CNFireParams,
    pftcon::PftConFireBase,
    # Masks and bounds
    mask_soilc::BitVector,
    mask_soilp::BitVector,
    mask_exposedveg::BitVector,
    mask_noexposedveg::BitVector,
    bounds_c::UnitRange{Int},
    bounds_p::UnitRange{Int},
    # Core type data
    patch::PatchData,
    col::ColumnData,
    grc::GridcellData,
    # Soil state (for root wetness)
    soilstate::SoilStateData,
    h2osoi_vol_col::Matrix{Float64},
    # CN Veg State and Carbon State
    cnveg_state::CNVegStateData,
    cnveg_cs::CNVegCarbonStateData,
    # Decomposition data
    decomp_cascade_con::DecompCascadeConData,
    # Input arrays
    totlitc_col::Vector{Float64},
    decomp_cpools_vr_col::Array{Float64,3},
    t_soi17cm_col::Vector{Float64};
    # Atmospheric forcing (keyword args)
    forc_rh_grc::Vector{Float64} = Float64[],
    forc_wind_grc::Vector{Float64} = Float64[],
    forc_t_col::Vector{Float64} = Float64[],
    forc_rain_col::Vector{Float64} = Float64[],
    forc_snow_col::Vector{Float64} = Float64[],
    prec60_patch::Vector{Float64} = Float64[],
    prec10_patch::Vector{Float64} = Float64[],
    # Saturated runoff
    fsat_col::Vector{Float64} = Float64[],
    # Water diagnostic
    wf_col::Vector{Float64} = Float64[],
    wf2_col::Vector{Float64} = Float64[],
    # Time/date
    dt::Float64 = 1800.0,
    dayspyr::Float64 = 365.0,
    kmo::Int = 1,
    kda::Int = 1,
    mcsec::Int = 0,
    nstep::Int = 1,
    nlevgrnd::Int = 10,
    nlevdecomp::Int = 1,
    ndecomp_pools::Int = 7,
    # PFT indices
    nc4_grass::Int = 14,
    nc3crop::Int = 15,
    ndllf_evr_tmp_tree::Int = 1,
    nbrdlf_evr_trp_tree::Int = 4,
    nbrdlf_dcd_trp_tree::Int = 6,
    nbrdlf_evr_shrub::Int = 9,
    noveg::Int = 0,
    # Control flags
    transient_landcover::Bool = false,
    # Soil suction function
    soil_suction_fn::Function = default_soil_suction
)
    # --- Aliases (matching Fortran associate block) ---
    totlitc          = totlitc_col
    decomp_cpools_vr = decomp_cpools_vr_col
    tsoi17           = t_soi17cm_col
    lfuel            = cnfire_const.lfuel
    ufuel            = cnfire_const.ufuel
    rh_hgh           = cnfire_const.rh_hgh
    rh_low           = cnfire_const.rh_low
    bt_min           = cnfire_const.bt_min
    bt_max           = cnfire_const.bt_max
    cli_scale        = cnfire_const.cli_scale
    cropfire_a1      = cnfire_const.cropfire_a1
    non_boreal_peatfire_c    = cnfire_const.non_boreal_peatfire_c
    pot_hmn_ign_counts_alpha = cnfire_const.pot_hmn_ign_counts_alpha
    boreal_peatfire_c        = cnfire_const.boreal_peatfire_c
    defo_fire_precip_thresh_bet  = cnfire_const.defo_fire_precip_thresh_bet
    defo_fire_precip_thresh_bdt  = cnfire_const.defo_fire_precip_thresh_bdt
    borpeat_fire_soilmoist_denom = cnfire_const.borpeat_fire_soilmoist_denom
    nonborpeat_fire_precip_denom = cnfire_const.nonborpeat_fire_precip_denom

    fsr_pft          = pftcon_li2014.fsr_pft
    fd_pft           = pftcon_li2014.fd_pft

    btran2           = fire_data.btran2_patch
    gdp_lf           = fire_li2014.gdp_lf_col
    peatf_lf         = fire_li2014.peatf_lf_col
    abm_lf           = fire_li2014.abm_lf_col
    is_cwd           = decomp_cascade_con.is_cwd

    forc_rh          = forc_rh_grc
    forc_wind        = forc_wind_grc
    forc_t           = forc_t_col
    forc_rain        = forc_rain_col
    forc_snow        = forc_snow_col
    prec60           = prec60_patch
    prec10           = prec10_patch
    fsat             = fsat_col
    wf               = wf_col
    wf2              = wf2_col

    dwt_smoothed     = cnveg_state.dwt_smoothed_patch
    cropf_col        = cnveg_state.cropf_col
    baf_crop         = cnveg_state.baf_crop_col
    baf_peatf        = cnveg_state.baf_peatf_col
    burndate         = cnveg_state.burndate_patch
    fbac             = cnveg_state.fbac_col
    fbac1            = cnveg_state.fbac1_col
    farea_burned     = cnveg_state.farea_burned_col
    nfire            = cnveg_state.nfire_col
    fsr_col          = cnveg_state.fsr_col
    fd_col           = cnveg_state.fd_col
    lgdp_col         = cnveg_state.lgdp_col
    lgdp1_col        = cnveg_state.lgdp1_col
    lpop_col         = cnveg_state.lpop_col
    lfwt             = cnveg_state.lfwt_col
    trotr1_col       = cnveg_state.trotr1_col
    trotr2_col       = cnveg_state.trotr2_col
    dtrotr_col       = cnveg_state.dtrotr_col
    lfc              = cnveg_state.lfc_col
    wtlf             = cnveg_state.wtlf_col

    totvegc          = cnveg_cs.totvegc_col
    rootc_col        = cnveg_cs.rootc_col
    leafc_col        = cnveg_cs.leafc_col
    fuelc            = cnveg_cs.fuelc_col
    fuelc_crop       = cnveg_cs.fuelc_crop_col

    deadcrootc       = cnveg_cs.deadcrootc_patch
    deadcrootc_storage = cnveg_cs.deadcrootc_storage_patch
    deadcrootc_xfer  = cnveg_cs.deadcrootc_xfer_patch
    frootc           = cnveg_cs.frootc_patch
    frootc_storage   = cnveg_cs.frootc_storage_patch
    frootc_xfer      = cnveg_cs.frootc_xfer_patch
    livecrootc       = cnveg_cs.livecrootc_patch
    livecrootc_storage = cnveg_cs.livecrootc_storage_patch
    livecrootc_xfer  = cnveg_cs.livecrootc_xfer_patch
    leafc            = cnveg_cs.leafc_patch
    leafc_storage    = cnveg_cs.leafc_storage_patch
    leafc_xfer       = cnveg_cs.leafc_xfer_patch

    # --- Local arrays ---
    nc = length(bounds_c)
    btran_col = zeros(last(bounds_c))

    # Temporary arrays for p2c results
    prec60_col = zeros(last(bounds_c))
    prec10_col = zeros(last(bounds_c))

    # --- Patch to column averaging ---
    p2c!(prec10_col, prec10, patch, mask_soilc, bounds_c, bounds_p)
    p2c!(prec60_col, prec60, patch, mask_soilc, bounds_c, bounds_p)
    p2c!(leafc_col,  leafc,  patch, mask_soilc, bounds_c, bounds_p)

    # --- On first time-step, set area burned to zero and exit ---
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

    # --- Calculate cropf_col and lfwt ---
    for c in bounds_c
        mask_soilc[c] || continue
        cropf_col[c] = 0.0
        lfwt[c]      = 0.0
    end
    for p in bounds_p
        mask_soilp[p] || continue
        c = patch.column[p]
        if patch.itype[p] > nc4_grass
            cropf_col[c] = cropf_col[c] + patch.wtcol[p]
        end
        if patch.itype[p] >= ndllf_evr_tmp_tree && patch.itype[p] <= nc4_grass
            lfwt[c] = lfwt[c] + patch.wtcol[p]
        end
    end

    # --- Calculate crop fuel ---
    for c in bounds_c
        mask_soilc[c] || continue
        fuelc_crop[c] = 0.0
    end
    for p in bounds_p
        mask_soilp[p] || continue
        c = patch.column[p]
        if patch.itype[p] > nc4_grass && patch.wtcol[p] > 0.0 && leafc_col[c] > 0.0
            fuelc_crop[c] = fuelc_crop[c] +
                (leafc[p] + leafc_storage[p] + leafc_xfer[p]) *
                    patch.wtcol[p] / cropf_col[c] +
                totlitc[c] * leafc[p] / leafc_col[c] *
                    patch.wtcol[p] / cropf_col[c]
        end
    end

    # --- Initialize noncrop column variables ---
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

    # --- Calculate root wetness (btran2) ---
    cnfire_calc_fire_root_wetness_li2014!(
        fire_data,
        mask_exposedveg,
        mask_noexposedveg,
        bounds_p,
        pftcon,
        patch,
        soilstate,
        h2osoi_vol_col,
        nlevgrnd;
        soil_suction_fn = soil_suction_fn
    )

    # --- Accumulate btran_col and wtlf from exposed vegetation ---
    for p in bounds_p
        mask_exposedveg[p] || continue
        c = patch.column[p]
        if patch.itype[p] < nc3crop && cropf_col[c] < 1.0
            btran_col[c] = btran_col[c] + btran2[p] * patch.wtcol[p]
            wtlf[c]      = wtlf[c] + patch.wtcol[p]
        end
    end

    # --- Main patch loop for noncrop column variables ---
    for p in bounds_p
        mask_soilp[p] || continue
        c = patch.column[p]
        g = col.gridcell[c]

        # For non-crop -- natural vegetation and bare-soil
        if patch.itype[p] < nc3crop && cropf_col[c] < 1.0

            # Tropical tree fractions
            if patch.itype[p] == nbrdlf_evr_trp_tree && patch.wtcol[p] > 0.0
                trotr1_col[c] = trotr1_col[c] + patch.wtcol[p]
            end
            if patch.itype[p] == nbrdlf_dcd_trp_tree && patch.wtcol[p] > 0.0
                trotr2_col[c] = trotr2_col[c] + patch.wtcol[p]
            end

            # Transient landcover: decreased tropical tree coverage
            if transient_landcover
                if patch.itype[p] == nbrdlf_evr_trp_tree || patch.itype[p] == nbrdlf_dcd_trp_tree
                    if dwt_smoothed[p] < 0.0
                        dtrotr_col[c] = dtrotr_col[c] - dwt_smoothed[p]
                    end
                end
            end

            # Root carbon accumulation
            rootc_col[c] = rootc_col[c] + (frootc[p] + frootc_storage[p] +
                frootc_xfer[p] + deadcrootc[p] +
                deadcrootc_storage[p] + deadcrootc_xfer[p] +
                livecrootc[p] + livecrootc_storage[p] +
                livecrootc_xfer[p]) * patch.wtcol[p]

            # Fire spread rate
            fsr_col[c] = fsr_col[c] + fsr_pft[patch.itype[p]] * patch.wtcol[p] / (1.0 - cropf_col[c])

            if lfwt[c] != 0.0
                hdmlf = fire_li2014.forc_hdm[g]

                if hdmlf > 0.1
                    # For NOT bare-soil
                    if patch.itype[p] != noveg
                        # For shrub and grass
                        if patch.itype[p] >= nbrdlf_evr_shrub
                            lgdp_col[c]  = lgdp_col[c] + (0.1 + 0.9 *
                                exp(-1.0 * RPI *
                                (gdp_lf[c] / 8.0)^0.5)) * patch.wtcol[p] /
                                (1.0 - cropf_col[c])
                            lgdp1_col[c] = lgdp1_col[c] + (0.2 + 0.8 *
                                exp(-1.0 * RPI *
                                (gdp_lf[c] / 7.0))) * patch.wtcol[p] / lfwt[c]
                            lpop_col[c]  = lpop_col[c] + (0.2 + 0.8 *
                                exp(-1.0 * RPI *
                                (hdmlf / 450.0)^0.5)) * patch.wtcol[p] / lfwt[c]
                        else   # for trees
                            if gdp_lf[c] > 20.0
                                lgdp_col[c] = lgdp_col[c] + cnfire_const.occur_hi_gdp_tree *
                                    patch.wtcol[p] / (1.0 - cropf_col[c])
                            else
                                lgdp_col[c] = lgdp_col[c] + patch.wtcol[p] / (1.0 - cropf_col[c])
                            end
                            if gdp_lf[c] > 20.0
                                lgdp1_col[c] = lgdp1_col[c] + 0.62 * patch.wtcol[p] / lfwt[c]
                            else
                                if gdp_lf[c] > 8.0
                                    lgdp1_col[c] = lgdp1_col[c] + 0.83 * patch.wtcol[p] / lfwt[c]
                                else
                                    lgdp1_col[c] = lgdp1_col[c] + patch.wtcol[p] / lfwt[c]
                                end
                            end
                            lpop_col[c] = lpop_col[c] + (0.4 + 0.6 *
                                exp(-1.0 * RPI *
                                (hdmlf / 125.0))) * patch.wtcol[p] / lfwt[c]
                        end
                    end
                else
                    lgdp_col[c]  = lgdp_col[c] + patch.wtcol[p] / (1.0 - cropf_col[c])
                    lgdp1_col[c] = lgdp1_col[c] + patch.wtcol[p] / lfwt[c]
                    lpop_col[c]  = lpop_col[c] + patch.wtcol[p] / lfwt[c]
                end
            end

            # Fire duration
            fd_col[c] = fd_col[c] + fd_pft[patch.itype[p]] * patch.wtcol[p] * SECSPHR / (1.0 - cropf_col[c])
        end
    end

    # --- Estimate annual decreased fractional coverage of BET+BDT ---
    if transient_landcover
        for c in bounds_c
            mask_soilc[c] || continue
            if dtrotr_col[c] > 0.0
                if kmo == 1 && kda == 1 && mcsec == 0
                    lfc[c] = 0.0
                end
                if kmo == 1 && kda == 1 && mcsec == dt
                    lfc[c] = dtrotr_col[c] * dayspyr * SECSPDAY / dt
                end
            else
                lfc[c] = 0.0
            end
        end
    end

    # --- Calculate burned area fraction in cropland ---
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
        # For crop
        if forc_t[c] >= TFRZ && patch.itype[p] > nc4_grass &&
           kmo == abm_lf[c] && forc_rain[c] + forc_snow[c] == 0.0 &&
           burndate[p] >= 999 && patch.wtcol[p] > 0.0

            hdmlf = fire_li2014.forc_hdm[g]

            # human density impact on agricultural fire
            fhd = 0.04 + 0.96 * exp(-1.0 * RPI * (hdmlf / 350.0)^0.5)

            # impact of GDP on agricultural fire
            fgdp = 0.01 + 0.99 * exp(-1.0 * RPI * (gdp_lf[c] / 10.0))

            # burned area
            fb = max(0.0, min(1.0, (fuelc_crop[c] - lfuel) / (ufuel - lfuel)))

            baf_crop[c] = baf_crop[c] + cropfire_a1 / SECSPHR * fb * fhd * fgdp * patch.wtcol[p]
            if fb * fhd * fgdp * patch.wtcol[p] > 0.0
                burndate[p] = kda
            end
        end
    end

    # --- Calculate peatland fire ---
    for c in bounds_c
        mask_soilc[c] || continue
        g = col.gridcell[c]
        if grc.latdeg[g] < cnfire_const.borealat
            baf_peatf[c] = non_boreal_peatfire_c / SECSPHR *
                max(0.0, min(1.0,
                    (4.0 - prec60_col[c] * SECSPDAY / nonborpeat_fire_precip_denom) /
                    4.0))^2 * peatf_lf[c] * (1.0 - fsat[c])
        else
            baf_peatf[c] = boreal_peatfire_c / SECSPHR *
                exp(-RPI * (max(wf2[c], 0.0) / borpeat_fire_soilmoist_denom)) *
                max(0.0, min(1.0, (tsoi17[c] - TFRZ) / 10.0)) * peatf_lf[c] *
                (1.0 - fsat[c])
        end
    end

    # --- Find which pool is the CWD pool ---
    i_cwd = 0
    for l in 1:ndecomp_pools
        if is_cwd[l]
            i_cwd = l
        end
    end

    # --- Main column loop: fractional area affected by fire ---
    for c in bounds_c
        mask_soilc[c] || continue
        g = col.gridcell[c]
        hdmlf = fire_li2014.forc_hdm[g]
        nfire[c] = 0.0

        if cropf_col[c] < 1.0
            if trotr1_col[c] + trotr2_col[c] > 0.6
                farea_burned[c] = min(1.0, baf_crop[c] + baf_peatf[c])
            else
                fuelc[c] = totlitc[c] + totvegc[c] - rootc_col[c] - fuelc_crop[c] * cropf_col[c]
                for j in 1:nlevdecomp
                    fuelc[c] = fuelc[c] + decomp_cpools_vr[c, j, i_cwd] * dzsoi_decomp[][j]
                end
                fuelc[c] = fuelc[c] / (1.0 - cropf_col[c])
                fb       = max(0.0, min(1.0, (fuelc[c] - lfuel) / (ufuel - lfuel)))
                m        = max(0.0, wf[c])
                fire_m   = exp(-RPI * (m / 0.69)^2) * (1.0 - max(0.0,
                    min(1.0, (forc_rh[g] - rh_low) / (rh_hgh - rh_low)))) *
                    min(1.0, exp(RPI * (forc_t[c] - TFRZ) / 10.0))
                lh       = pot_hmn_ign_counts_alpha * 6.8 * hdmlf^0.43 / 30.0 / 24.0
                fs       = 1.0 - (0.01 + 0.98 * exp(-0.025 * hdmlf))
                ig       = (lh + fire_li2014.forc_lnfm[g] /
                    (5.16 + 2.16 * cos(3.0 * grc.lat[g])) *
                    cnfire_params.ignition_efficiency) * (1.0 - fs) * (1.0 - cropf_col[c])
                nfire[c] = ig / SECSPHR * fb * fire_m * lgdp_col[c]
                Lb_lf    = 1.0 + 10.0 * (1.0 - exp(-0.06 * forc_wind[g]))
                if wtlf[c] > 0.0
                    spread_m = (1.0 - max(0.0, min(1.0,
                        (btran_col[c] / wtlf[c] - bt_min) /
                        (bt_max - bt_min)))) * (1.0 - max(0.0,
                        min(1.0, (forc_rh[g] - rh_low) / (rh_hgh - rh_low))))
                else
                    spread_m = 0.0
                end
                farea_burned[c] = min(1.0,
                    (cnfire_const.g0 * spread_m * fsr_col[c] *
                    fd_col[c] / 1000.0)^2 * lgdp1_col[c] *
                    lpop_col[c] * nfire[c] * RPI * Lb_lf +
                    baf_crop[c] + baf_peatf[c])
            end

            # Deforestation fires (transient landcover)
            if transient_landcover
                if trotr1_col[c] + trotr2_col[c] > 0.6
                    if (kmo == 1 && kda == 1 && mcsec == 0) ||
                       dtrotr_col[c] <= 0.0
                        fbac1[c]        = 0.0
                        farea_burned[c] = baf_crop[c] + baf_peatf[c]
                    else
                        # Calculate precip threshold as area-weighted mean of BET and BDT
                        cri = (defo_fire_precip_thresh_bet * trotr1_col[c] +
                               defo_fire_precip_thresh_bdt * trotr2_col[c]) /
                              (trotr1_col[c] + trotr2_col[c])

                        cli = (max(0.0, min(1.0, (cri - prec60_col[c] * SECSPDAY) / cri))^0.5) *
                              (max(0.0, min(1.0, (cri - prec10_col[c] * SECSPDAY) / cri))^0.5) *
                              max(0.0005, min(1.0, 19.0 * dtrotr_col[c] * dayspyr * SECSPDAY / dt - 0.001)) *
                              max(0.0, min(1.0, (0.25 - (forc_rain[c] + forc_snow[c]) * SECSPHR) / 0.25))
                        farea_burned[c] = cli * (cli_scale / SECSPDAY) + baf_crop[c] + baf_peatf[c]
                        # burned area out of conversion region due to land use fire
                        fbac1[c] = max(0.0, cli * (cli_scale / SECSPDAY) - 2.0 * lfc[c] / dt)
                    end
                    # total burned area out of conversion
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
