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
# GPU: cnfire_area_li2021! runs entirely through KernelAbstractions kernels on
# the _launch! path (byte-identical on CPU, Metal-validated). Reuses the generic
# Li2014 sections + the Li2016 noncrop/crop-fire kernels (identical formulas),
# and adds Li2021-specific kernels for the rswf-rescaled btran accumulation and
# the fire-spread loop. Off-by-default (non-default method).
# ==========================================================================

"""
    need_lightning_and_popdens_li2021()

Returns `true`. Ported from `need_lightning_and_popdens` in `CNFireLi2021Mod.F90`.
"""
need_lightning_and_popdens_li2021() = true

# --- Section: btran_col / wtlf accumulation with per-PFT rswf rescale -------
# Li2021 rescales btran2 by the per-PFT (rswf_min, rswf_max) window before the
# area-weighted column sum (Li2016/2014 used the raw btran2).
@kernel function _firea2021_btran_wtlf_kernel!(btran_col, wtlf, @Const(mask_exposedveg),
        @Const(column), @Const(itype), @Const(wtcol), @Const(btran2), @Const(cropf_col),
        @Const(rswf_min), @Const(rswf_max), nc3crop::Int)
    p = @index(Global)
    T = eltype(btran_col)
    @inbounds if mask_exposedveg[p]
        c = column[p]
        if itype[p] < nc3crop && cropf_col[c] < one(T)
            it = itype[p]
            _scatter_add!(btran_col, c, smooth_max(zero(T), smooth_min(one(T),
                (btran2[p] - rswf_min[it + 1]) / (rswf_max[it + 1] - rswf_min[it + 1]))) * wtcol[p])
            _scatter_add!(wtlf, c, wtcol[p])
        end
    end
end

# --- Section: main per-column fire-spread loop (Li2021) --------------------
# Reuses the Li2016 device-view bundles (_Fire2016SpreadOut/In1/In2/Mat/Scalars);
# the array set is identical. Differs from Li2016 only in the fire-occurrence
# btran term (no bt_min/bt_max clamp) and the deforestation climate term `cli`
# (revised dtrotr formula; farea/fbac1 use fb*cli). bt_min/bt_max are unused.
@kernel function _firea2021_fire_spread_kernel!(
        out::_Fire2016SpreadOut, in1::_Fire2016SpreadIn1, in2::_Fire2016SpreadIn2,
        mat::_Fire2016SpreadMat,
        @Const(mask_soilc), @Const(gridcell),
        sc::_Fire2016SpreadScalars,
        i_cwd::Int, nlevdecomp::Int, spinup_state::Int,
        date_is_jan1_0::Bool, transient_landcover::Bool)
    c = @index(Global)
    @inbounds if mask_soilc[c]
        T = eltype(out.nfire)

        nfire = out.nfire; farea_burned = out.farea_burned; fuelc = out.fuelc
        fbac  = out.fbac;  fbac1 = out.fbac1

        forc_hdm   = in1.forc_hdm;   cropf_col  = in1.cropf_col
        trotr1_col = in1.trotr1_col; trotr2_col = in1.trotr2_col
        baf_crop   = in1.baf_crop;   baf_peatf  = in1.baf_peatf
        totlitc    = in1.totlitc;    totvegc    = in1.totvegc
        rootc_col  = in1.rootc_col;  fuelc_crop = in1.fuelc_crop
        dzsoi_decomp = in1.dzsoi_decomp; rh30_col = in1.rh30_col
        forc_rh    = in1.forc_rh;    wtlf = in1.wtlf; forc_lnfm = in1.forc_lnfm

        latdeg    = in2.latdeg;    forc_wind = in2.forc_wind
        lgdp_col  = in2.lgdp_col;  lgdp1_col = in2.lgdp1_col; lpop_col = in2.lpop_col
        fsr_col   = in2.fsr_col;   fd_col    = in2.fd_col;    btran_col = in2.btran_col
        tsoi17    = in2.tsoi17;    dtrotr_col = in2.dtrotr_col
        prec60_col = in2.prec60_col; prec10_col = in2.prec10_col
        forc_rain = in2.forc_rain; forc_snow = in2.forc_snow; lfc = in2.lfc
        deadstemc_col = in2.deadstemc_col; spinup_lat_term = in2.spinup_lat_term

        decomp_cpools_vr = mat.decomp_cpools_vr

        lfuel = sc.lfuel; ufuel = sc.ufuel; rh_low = sc.rh_low; rh_hgh = sc.rh_hgh
        prh30 = sc.prh30; max_rh30_affecting_fuel = sc.max_rh30_affecting_fuel
        pot_hmn_ign_counts_alpha = sc.pot_hmn_ign_counts_alpha
        ignition_efficiency = sc.ignition_efficiency
        g0 = sc.g0; cli_scale = sc.cli_scale
        defo_fire_precip_thresh_bet = sc.defo_fire_precip_thresh_bet
        defo_fire_precip_thresh_bdt = sc.defo_fire_precip_thresh_bdt
        spinup_factor_deadwood = sc.spinup_factor_deadwood
        spinup_factor_cwd = sc.spinup_factor_cwd
        tfrz = sc.tfrz; rpi = sc.rpi; secsphr = sc.secsphr; secspday = sc.secspday
        dayspyr = sc.dayspyr; dt = sc.dt

        g = gridcell[c]
        hdmlf = forc_hdm[g]
        nfire[c] = zero(T)
        if cropf_col[c] < one(T)
            fc = totlitc[c] + totvegc[c] - rootc_col[c] - fuelc_crop[c] * cropf_col[c]
            if spinup_state == 2
                fc = fc + (spinup_factor_deadwood - one(T)) * deadstemc_col[c]
                for j in 1:nlevdecomp
                    fc = fc + decomp_cpools_vr[c, j, i_cwd] * dzsoi_decomp[j] *
                              spinup_factor_cwd * spinup_lat_term[c]
                end
            else
                for j in 1:nlevdecomp
                    fc = fc + decomp_cpools_vr[c, j, i_cwd] * dzsoi_decomp[j]
                end
            end
            fc = fc / (one(T) - cropf_col[c])
            fuelc[c] = fc
            fb = smooth_max(zero(T), smooth_min(one(T), (fc - lfuel) / (ufuel - lfuel)))
            if trotr1_col[c] + trotr2_col[c] <= T(0.6)
                afuel = smooth_min(one(T), smooth_max(zero(T), (fc - T(2500.0)) / (T(5000.0) - T(2500.0))))
                arh   = one(T) - smooth_max(zero(T), smooth_min(one(T), (forc_rh[g] - rh_low) / (rh_hgh - rh_low)))
                arh30 = one(T) - smooth_max(prh30, smooth_min(one(T), rh30_col[c] / max_rh30_affecting_fuel))
                if forc_rh[g] < rh_hgh && wtlf[c] > zero(T) && tsoi17[c] > tfrz
                    # Li2021: btran term has no bt_min/bt_max clamp.
                    fire_m = ((afuel * arh30 + (one(T) - afuel) * arh)^T(1.5)) *
                             ((one(T) - btran_col[c] / wtlf[c])^T(0.5))
                else
                    fire_m = zero(T)
                end
                lh = pot_hmn_ign_counts_alpha * T(6.8) * hdmlf^T(0.43) / T(30.0) / T(24.0)
                fs = one(T) - (T(0.01) + T(0.98) * exp(T(-0.025) * hdmlf))
                ig = (lh + forc_lnfm[g] /
                     (T(5.16) + T(2.16) * cos(rpi / T(180.0) * T(3) * smooth_min(T(60.0), abs(latdeg[g])))) *
                     ignition_efficiency) * (one(T) - fs) * (one(T) - cropf_col[c])
                nfire[c] = ig / secsphr * fb * fire_m * lgdp_col[c]
                Lb_lf = one(T) + T(10.0) * (one(T) - exp(T(-0.06) * forc_wind[g]))
                spread_m = fire_m^T(0.5)
                farea_burned[c] = smooth_min(one(T), (g0 * spread_m * fsr_col[c] *
                    fd_col[c] / T(1000.0))^2 * lgdp1_col[c] * lpop_col[c] *
                    nfire[c] * rpi * Lb_lf + baf_crop[c] + baf_peatf[c])
            else
                farea_burned[c] = smooth_min(one(T), baf_crop[c] + baf_peatf[c])
            end

            if transient_landcover
                if trotr1_col[c] + trotr2_col[c] > T(0.6)
                    if date_is_jan1_0 || dtrotr_col[c] <= zero(T)
                        fbac1[c]        = zero(T)
                        farea_burned[c] = baf_crop[c] + baf_peatf[c]
                    else
                        cri = (defo_fire_precip_thresh_bet * trotr1_col[c] +
                               defo_fire_precip_thresh_bdt * trotr2_col[c]) /
                              (trotr1_col[c] + trotr2_col[c])
                        # Li2021: revised deforestation climate term.
                        cli = (smooth_max(zero(T), smooth_min(one(T), (cri - prec60_col[c] * secspday) / cri))^T(0.5)) *
                              (smooth_max(zero(T), smooth_min(one(T), (cri - prec10_col[c] * secspday) / cri))^T(0.5)) *
                              (T(15.0) * smooth_min(T(0.0016), dtrotr_col[c] / dt * dayspyr * secspday) + T(0.009)) *
                              smooth_max(zero(T), smooth_min(one(T), (T(0.25) - (forc_rain[c] + forc_snow[c]) * secsphr) / T(0.25)))
                        farea_burned[c] = fb * cli * (cli_scale / secspday) + baf_crop[c] + baf_peatf[c]
                        fbac1[c] = max(zero(T), fb * cli * (cli_scale / secspday) - T(2) * lfc[c] / dt)
                    end
                    fbac[c] = fbac1[c] + baf_crop[c] + baf_peatf[c]
                else
                    fbac[c] = farea_burned[c]
                end
            end
        else
            farea_burned[c] = smooth_min(one(T), baf_crop[c] + baf_peatf[c])
        end
    end
end

function _firea2021_fire_spread!(
        nfire, farea_burned, fuelc, fbac, fbac1,
        forc_hdm, cropf_col, trotr1_col, trotr2_col, baf_crop, baf_peatf,
        totlitc, totvegc, rootc_col, fuelc_crop, dzsoi_decomp, rh30_col,
        forc_rh, wtlf, forc_lnfm,
        latdeg, forc_wind, lgdp_col, lgdp1_col, lpop_col, fsr_col, fd_col,
        btran_col, tsoi17, dtrotr_col, prec60_col, prec10_col, forc_rain, forc_snow,
        lfc, deadstemc_col, spinup_lat_term,
        decomp_cpools_vr, mask_soilc, gridcell,
        lfuel, ufuel, rh_low, rh_hgh, prh30, max_rh30_affecting_fuel,
        pot_hmn_ign_counts_alpha, ignition_efficiency, g0, cli_scale,
        defo_fire_precip_thresh_bet, defo_fire_precip_thresh_bdt,
        spinup_factor_deadwood, spinup_factor_cwd,
        tfrz, rpi, secsphr, secspday, dayspyr, dt,
        i_cwd::Int, nlevdecomp::Int, spinup_state::Int,
        date_is_jan1_0::Bool, transient_landcover::Bool; ndrange = length(farea_burned))

    (ndrange isa Integer ? ndrange == 0 : prod(ndrange) == 0) && return nothing

    T = eltype(farea_burned)
    out = _Fire2016SpreadOut(; nfire, farea_burned, fuelc, fbac, fbac1)
    in1 = _Fire2016SpreadIn1(; forc_hdm, cropf_col, trotr1_col, trotr2_col, baf_crop,
        baf_peatf, totlitc, totvegc, rootc_col, fuelc_crop, dzsoi_decomp, rh30_col,
        forc_rh, wtlf, forc_lnfm)
    in2 = _Fire2016SpreadIn2(; latdeg, forc_wind, lgdp_col, lgdp1_col, lpop_col,
        fsr_col, fd_col, btran_col, tsoi17, dtrotr_col, prec60_col, prec10_col,
        forc_rain, forc_snow, lfc, deadstemc_col, spinup_lat_term)
    mat = _Fire2016SpreadMat(; decomp_cpools_vr)
    sc = _Fire2016SpreadScalars(;
        lfuel=T(lfuel), ufuel=T(ufuel), rh_low=T(rh_low), rh_hgh=T(rh_hgh),
        bt_min=zero(T), bt_max=zero(T), prh30=T(prh30),
        max_rh30_affecting_fuel=T(max_rh30_affecting_fuel),
        pot_hmn_ign_counts_alpha=T(pot_hmn_ign_counts_alpha),
        ignition_efficiency=T(ignition_efficiency), g0=T(g0), cli_scale=T(cli_scale),
        defo_fire_precip_thresh_bet=T(defo_fire_precip_thresh_bet),
        defo_fire_precip_thresh_bdt=T(defo_fire_precip_thresh_bdt),
        spinup_factor_deadwood=T(spinup_factor_deadwood),
        spinup_factor_cwd=T(spinup_factor_cwd),
        tfrz=T(tfrz), rpi=T(rpi), secsphr=T(secsphr), secspday=T(secspday),
        dayspyr=T(dayspyr), dt=T(dt))

    backend = _kernel_backend(out.nfire)
    _firea2021_fire_spread_kernel!(backend)(out, in1, in2, mat,
        mask_soilc, gridcell, sc,
        i_cwd, nlevdecomp, spinup_state, date_is_jan1_0, transient_landcover;
        ndrange = ndrange)
    KA.synchronize(backend)
    return nothing
end

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
    btran_col  = fill!(similar(totvegc, FT, nc), zero(FT))
    prec60_col = fill!(similar(totvegc, FT, nc), zero(FT))
    prec10_col = fill!(similar(totvegc, FT, nc), zero(FT))
    rh30_col   = fill!(similar(totvegc, FT, nc), zero(FT))

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

    date_is_jan1_0  = (kmo == 1 && kda == 1 && mcsec == 0)
    date_is_jan1_dt = (kmo == 1 && kda == 1 && mcsec == dt)

    # --- cropf_col + lfwt ---
    _launch!(_firea_zero_cropf_lfwt_kernel!, cropf_col, lfwt, mask_soilc)
    _launch!(_firea_cropf_lfwt_kernel!, cropf_col, lfwt, mask_soilp,
             patch.column, patch.itype, patch.wtcol,
             nc4_grass, ndllf_evr_tmp_tree; ndrange = length(mask_soilp))

    # --- crop fuel ---
    _launch!(_firea_zero_fuelc_crop_kernel!, fuelc_crop, mask_soilc)
    _launch!(_firea_fuelc_crop_kernel!, fuelc_crop, mask_soilp,
             patch.column, patch.itype, patch.wtcol,
             leafc, leafc_storage, leafc_xfer, leafc_col, totlitc, cropf_col,
             nc4_grass; ndrange = length(mask_soilp))

    # --- init noncrop column variables ---
    _launch!(_firea_init_noncrop_kernel!, fsr_col, fd_col, rootc_col, lgdp_col,
             lgdp1_col, lpop_col, btran_col, wtlf, trotr1_col, trotr2_col,
             dtrotr_col, mask_soilc, transient_landcover)

    # --- root wetness (btran2), Li2021 method ---
    cnfire_calc_fire_root_wetness_li2021!(
        fire_data, mask_exposedveg, mask_noexposedveg, bounds_p,
        patch, soilstate, h2osoi_vol_col, nlevgrnd)

    # --- accumulate btran_col / wtlf with per-PFT rswf rescale (Li2021) ---
    _launch!(_firea2021_btran_wtlf_kernel!, btran_col, wtlf, mask_exposedveg,
             patch.column, patch.itype, patch.wtcol, btran2, cropf_col,
             rswf_min, rswf_max, nc3crop; ndrange = length(mask_exposedveg))

    # --- main noncrop per-patch column-var accumulation (shared Li2016 form) ---
    _firea2016_noncrop_main!(
        trotr1_col, trotr2_col, dtrotr_col, rootc_col, fsr_col,
        lgdp_col, lgdp1_col, lpop_col, fd_col,
        mask_soilp, patch.column, col.gridcell, patch.itype, patch.wtcol,
        cropf_col, lfwt, dwt_smoothed,
        frootc, frootc_storage, frootc_xfer,
        deadcrootc, deadcrootc_storage, deadcrootc_xfer,
        livecrootc, livecrootc_storage, livecrootc_xfer,
        fsr_pft, fd_pft, gdp_lf, forc_hdm,
        nc3crop, nbrdlf_evr_trp_tree, nbrdlf_dcd_trp_tree, nbrdlf_evr_shrub, noveg,
        occur_hi_gdp_tree, RPI, secsphr, spinup_factor_deadwood,
        transient_landcover; ndrange = length(mask_soilp))

    # --- annual decreased fractional coverage of BET+BDT ---
    if transient_landcover
        _launch!(_firea_lfc_kernel!, lfc, mask_soilc, dtrotr_col,
                 date_is_jan1_0, date_is_jan1_dt,
                 eltype(lfc)(dayspyr), eltype(lfc)(secspday), eltype(lfc)(dt))
    end

    # --- burned-area fraction in cropland (shared Li2016 form) ---
    _launch!(_firea_zero_baf_crop_kernel!, baf_crop, mask_soilc)
    _launch!(_firea_burndate_init_kernel!, burndate, mask_soilp,
             date_is_jan1_0; ndrange = length(mask_soilp))
    _launch!(_firea2016_crop_fire_kernel!, baf_crop, burndate, mask_soilp,
             patch.column, col.gridcell, patch.itype, patch.wtcol,
             forc_t, abm_lf, fuelc_crop, gdp_lf, forc_hdm,
             nc4_grass, kmo, kda,
             eltype(baf_crop)(TFRZ), eltype(baf_crop)(RPI), eltype(baf_crop)(lfuel),
             eltype(baf_crop)(ufuel), eltype(baf_crop)(cropfire_a1), eltype(baf_crop)(secsphr);
             ndrange = length(mask_soilp))

    # --- peatland fire (per-column independent; shared Li2014 kernel) ---
    let TFp = eltype(baf_peatf)
        fire_peatland!(
            baf_peatf, mask_soilc, col.gridcell, grc.latdeg, prec60_col,
            peatf_lf, fsat, wf2, tsoi17,
            TFp(cnfire_const.borealat), TFp(non_boreal_peatfire_c), TFp(boreal_peatfire_c),
            TFp(nonborpeat_fire_precip_denom), TFp(borpeat_fire_soilmoist_denom),
            TFp(secsphr), TFp(secspday), TFp(TFRZ), TFp(RPI))
    end

    # --- find CWD pool (host-resolved) ---
    i_cwd = 0
    for l in 1:ndecomp_pools
        if is_cwd[l]
            i_cwd = l
        end
    end

    dzsoi_decomp_arr = copyto!(similar(totvegc, FT, length(dzsoi_decomp[])), dzsoi_decomp[])

    spinup_lat_term_col = fill!(similar(totvegc, FT, nc), zero(FT))
    fire_spinup_latterm!(spinup_lat_term_col, mask_soilc, col.gridcell, grc.latdeg)
    spinup_factor_cwd = (i_cwd >= 1 && i_cwd <= length(spinup_factor)) ?
                        spinup_factor[i_cwd] : one(FT)

    # --- main column loop: fractional area affected by fire (Li2021 form) ---
    _firea2021_fire_spread!(
        nfire, farea_burned, fuelc, fbac, fbac1,
        forc_hdm, cropf_col, trotr1_col, trotr2_col, baf_crop, baf_peatf,
        totlitc, totvegc, rootc_col, fuelc_crop, dzsoi_decomp_arr, rh30_col,
        forc_rh, wtlf, forc_lnfm,
        grc.latdeg, forc_wind, lgdp_col, lgdp1_col, lpop_col, fsr_col, fd_col,
        btran_col, tsoi17, dtrotr_col, prec60_col, prec10_col, forc_rain, forc_snow,
        lfc, deadstemc_col, spinup_lat_term_col,
        decomp_cpools_vr, mask_soilc, col.gridcell,
        lfuel, ufuel, rh_low, rh_hgh, prh30, max_rh30_affecting_fuel,
        pot_hmn_ign_counts_alpha, ignition_efficiency, g0, cli_scale,
        defo_fire_precip_thresh_bet, defo_fire_precip_thresh_bdt,
        spinup_factor_deadwood, spinup_factor_cwd,
        TFRZ, RPI, secsphr, secspday, dayspyr, dt,
        i_cwd, nlevdecomp, spinup_state, date_is_jan1_0, transient_landcover;
        ndrange = length(farea_burned))

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
