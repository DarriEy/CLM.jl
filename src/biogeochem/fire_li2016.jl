# ==========================================================================
# Ported from: src/biogeochem/CNFireLi2016Mod.F90
# Fire dynamics module for Li et al. (2016) version (calibrated 2015 against
# 20th-century transient runs with CRUNCEPv5 + climatological lightning).
#
# Incremental refinement of Li2014: same CNFireBase machinery + Li2014 root
# wetness, but a revised fire-occurrence formulation (afuel/arh/arh30 +
# bt_min/bt_max btran term), revised GDP/population limitation factors, a
# revised agricultural-fire term, and an AD-spinup fuel branch.
#
# Public functions:
#   need_lightning_and_popdens_li2016  — Returns true
#   cnfire_area_li2016!                — Compute column-level burned area
#   cnfire_fluxes_li2016!              — Fire C/N fluxes (delegates to base)
#
# GPU: cnfire_area_li2016! runs entirely through KernelAbstractions kernels on
# the _launch! path (byte-identical on CPU, Metal-validated). It reuses the
# generic Li2014 sections (p2c!, the zero/init/btran/lfc/burndate kernels, and
# fire_peatland!) and adds Li2016-specific kernels for the sections whose
# formulas differ (noncrop column-var accumulation, crop fire, fire spread).
# The default li2014 path is untouched. Off-by-default here (non-default method).
# ==========================================================================

"""
    need_lightning_and_popdens_li2016()

Returns `true` — the Li2016 fire model requires lightning + population density.
Ported from `need_lightning_and_popdens` in `CNFireLi2016Mod.F90`.
"""
need_lightning_and_popdens_li2016() = true

# ---------------------------------------------------------------------------
# Shared spinup latitude-term kernel (also used by li2021/li2024, included
# after this file). Equivalent to get_spinup_latitude_term(latdeg[gridcell[c]])
# but device-safe (no host scalar-indexing of a device latdeg array) and
# byte-identical on Float64. Only consumed by the spinup_state==2 fuel branch.
# ---------------------------------------------------------------------------
@kernel function _fire_spinup_latterm_kernel!(out, @Const(mask_soilc),
                                              @Const(gridcell), @Const(latdeg))
    c = @index(Global)
    T = eltype(out)
    @inbounds if mask_soilc[c]
        g = gridcell[c]
        out[c] = one(T) + T(50.0) / (one(T) + exp(T(-0.15) * (abs(latdeg[g]) - T(60.0))))
    end
end

function fire_spinup_latterm!(out, mask_soilc, gridcell, latdeg)
    _launch!(_fire_spinup_latterm_kernel!, out, mask_soilc, gridcell, latdeg;
             ndrange = length(out))
    return nothing
end

# ===========================================================================
# Li2016-specific kernels
# ===========================================================================

# --- Section: main noncrop per-patch column-var accumulation (Li2016) ------
# Reuses the Li2014 device-view bundles _FireNCOut / _FireNCIn (defined in
# fire_li2014.jl, included first). Li2016 differs from Li2014: deadcrootc is
# weighted by spinup_factor_deadwood, the GDP/pop limitation factors all divide
# by (1-cropf_col) (Li2014 used lfwt for lgdp1/lpop and gated on lfwt), and the
# tree-branch GDP coefficients differ (0.79/0.83 mid-tier).
Base.@kwdef struct _FireNC2016Scalars{T}
    occur_hi_gdp_tree::T
    rpi::T
    secsphr::T
    spinup_factor_deadwood::T
end

@kernel function _firea2016_noncrop_main_kernel!(out::_FireNCOut, inp::_FireNCIn,
        @Const(mask_soilp), @Const(column), @Const(gridcell), @Const(itype),
        sc::_FireNC2016Scalars,
        nc3crop::Int, nbrdlf_evr_trp_tree::Int, nbrdlf_dcd_trp_tree::Int,
        nbrdlf_evr_shrub::Int, noveg::Int, transient_landcover::Bool)
    p = @index(Global)
    @inbounds if mask_soilp[p]
        T = eltype(out.trotr1_col)

        trotr1_col = out.trotr1_col; trotr2_col = out.trotr2_col; dtrotr_col = out.dtrotr_col
        rootc_col  = out.rootc_col;  fsr_col    = out.fsr_col
        lgdp_col   = out.lgdp_col;   lgdp1_col  = out.lgdp1_col; lpop_col = out.lpop_col
        fd_col     = out.fd_col

        wtcol = inp.wtcol; cropf_col = inp.cropf_col; dwt_smoothed = inp.dwt_smoothed
        frootc = inp.frootc; frootc_storage = inp.frootc_storage; frootc_xfer = inp.frootc_xfer
        deadcrootc = inp.deadcrootc; deadcrootc_storage = inp.deadcrootc_storage
        deadcrootc_xfer = inp.deadcrootc_xfer
        livecrootc = inp.livecrootc; livecrootc_storage = inp.livecrootc_storage
        livecrootc_xfer = inp.livecrootc_xfer
        fsr_pft = inp.fsr_pft; fd_pft = inp.fd_pft; gdp_lf = inp.gdp_lf; forc_hdm = inp.forc_hdm

        occur_hi_gdp_tree = sc.occur_hi_gdp_tree; rpi = sc.rpi; secsphr = sc.secsphr
        spinup_factor_deadwood = sc.spinup_factor_deadwood

        c = column[p]
        g = gridcell[c]
        if itype[p] < nc3crop && cropf_col[c] < one(T)
            if itype[p] == nbrdlf_evr_trp_tree && wtcol[p] > zero(T)
                _scatter_add!(trotr1_col, c, wtcol[p])
            end
            if itype[p] == nbrdlf_dcd_trp_tree && wtcol[p] > zero(T)
                _scatter_add!(trotr2_col, c, wtcol[p])
            end
            if transient_landcover
                if itype[p] == nbrdlf_evr_trp_tree || itype[p] == nbrdlf_dcd_trp_tree
                    if dwt_smoothed[p] < zero(T)
                        _scatter_add!(dtrotr_col, c, -dwt_smoothed[p])
                    end
                end
            end
            # Li2016: deadcrootc weighted by spinup_factor_deadwood (unlike Li2014).
            _scatter_add!(rootc_col, c, (frootc[p] + frootc_storage[p] + frootc_xfer[p] +
                deadcrootc[p] * spinup_factor_deadwood +
                deadcrootc_storage[p] + deadcrootc_xfer[p] +
                livecrootc[p] + livecrootc_storage[p] + livecrootc_xfer[p]) * wtcol[p])

            _scatter_add!(fsr_col, c, fsr_pft[itype[p] + 1] * wtcol[p] / (one(T) - cropf_col[c]))

            hdmlf = forc_hdm[g]
            if hdmlf > T(0.1)
                if itype[p] != noveg
                    if itype[p] >= nbrdlf_evr_shrub  # shrub + grass
                        _scatter_add!(lgdp_col, c, (T(0.1) + T(0.9) *
                            exp(-one(T) * rpi * (gdp_lf[c] / T(8.0))^T(0.5))) * wtcol[p] /
                            (one(T) - cropf_col[c]))
                        _scatter_add!(lgdp1_col, c, (T(0.2) + T(0.8) *
                            exp(-one(T) * rpi * (gdp_lf[c] / T(7.0)))) * wtcol[p] /
                            (one(T) - cropf_col[c]))
                        _scatter_add!(lpop_col, c, (T(0.2) + T(0.8) *
                            exp(-one(T) * rpi * (hdmlf / T(450.0))^T(0.5))) * wtcol[p] /
                            (one(T) - cropf_col[c]))
                    else  # trees
                        if gdp_lf[c] > T(20.0)
                            _scatter_add!(lgdp_col, c, occur_hi_gdp_tree * wtcol[p] /
                                (one(T) - cropf_col[c]))
                            _scatter_add!(lgdp1_col, c, T(0.62) * wtcol[p] /
                                (one(T) - cropf_col[c]))
                        else
                            if gdp_lf[c] > T(8.0)
                                _scatter_add!(lgdp_col, c, T(0.79) * wtcol[p] /
                                    (one(T) - cropf_col[c]))
                                _scatter_add!(lgdp1_col, c, T(0.83) * wtcol[p] /
                                    (one(T) - cropf_col[c]))
                            else
                                _scatter_add!(lgdp_col, c, wtcol[p] / (one(T) - cropf_col[c]))
                                _scatter_add!(lgdp1_col, c, wtcol[p] / (one(T) - cropf_col[c]))
                            end
                        end
                        _scatter_add!(lpop_col, c, (T(0.4) + T(0.6) *
                            exp(-one(T) * rpi * (hdmlf / T(125.0)))) * wtcol[p] /
                            (one(T) - cropf_col[c]))
                    end
                end
            else
                _scatter_add!(lgdp_col, c, wtcol[p] / (one(T) - cropf_col[c]))
                _scatter_add!(lgdp1_col, c, wtcol[p] / (one(T) - cropf_col[c]))
                _scatter_add!(lpop_col, c, wtcol[p] / (one(T) - cropf_col[c]))
            end

            # fd_pft is a 0-based pftcon array in Fortran (`fd_pft(patch%itype(p))`),
            # so the Julia (1-based) index is itype+1 — same as fsr_pft above. Without
            # the +1 this read the WRONG PFT's fire duration (and, for a bare-soil
            # patch itype=0, indexed element 0 under @inbounds).
            _scatter_add!(fd_col, c, fd_pft[itype[p] + 1] * wtcol[p] * secsphr / (one(T) - cropf_col[c]))
        end
    end
end

function _firea2016_noncrop_main!(trotr1_col, trotr2_col, dtrotr_col,
        rootc_col, fsr_col, lgdp_col, lgdp1_col, lpop_col, fd_col,
        mask_soilp, column, gridcell, itype, wtcol,
        cropf_col, lfwt, dwt_smoothed,
        frootc, frootc_storage, frootc_xfer,
        deadcrootc, deadcrootc_storage, deadcrootc_xfer,
        livecrootc, livecrootc_storage, livecrootc_xfer,
        fsr_pft, fd_pft, gdp_lf, forc_hdm,
        nc3crop::Int, nbrdlf_evr_trp_tree::Int, nbrdlf_dcd_trp_tree::Int,
        nbrdlf_evr_shrub::Int, noveg::Int, occur_hi_gdp_tree, rpi, secsphr,
        spinup_factor_deadwood, transient_landcover::Bool; ndrange = length(mask_soilp))

    (ndrange isa Integer ? ndrange == 0 : prod(ndrange) == 0) && return nothing

    out = _FireNCOut(; trotr1_col, trotr2_col, dtrotr_col,
                       rootc_col, fsr_col, lgdp_col, lgdp1_col, lpop_col, fd_col)
    inp = _FireNCIn(; wtcol, cropf_col, lfwt, dwt_smoothed,
                      frootc, frootc_storage, frootc_xfer,
                      deadcrootc, deadcrootc_storage, deadcrootc_xfer,
                      livecrootc, livecrootc_storage, livecrootc_xfer,
                      fsr_pft, fd_pft, gdp_lf, forc_hdm)

    T = eltype(out.trotr1_col)
    sc = _FireNC2016Scalars(; occur_hi_gdp_tree = T(occur_hi_gdp_tree),
                              rpi = T(rpi), secsphr = T(secsphr),
                              spinup_factor_deadwood = T(spinup_factor_deadwood))

    backend = _kernel_backend(out.trotr1_col)
    _firea2016_noncrop_main_kernel!(backend)(out, inp,
        mask_soilp, column, gridcell, itype, sc,
        nc3crop, nbrdlf_evr_trp_tree, nbrdlf_dcd_trp_tree,
        nbrdlf_evr_shrub, noveg, transient_landcover; ndrange = ndrange)
    KA.synchronize(backend)
    return nothing
end

# --- Section: crop burned-area fraction (Li2016 per-patch scatter + burndate) -
# Li2016 differs from Li2014: gate omits the (rain+snow==0) condition and
# instead keys on burndate>=999; the ag-fire baf uses only fhd*fgdp (no fuel
# factor fb); fb still gates burndate.
@kernel function _firea2016_crop_fire_kernel!(baf_crop, burndate, @Const(mask_soilp),
        @Const(column), @Const(gridcell), @Const(itype), @Const(wtcol),
        @Const(forc_t), @Const(abm_lf), @Const(fuelc_crop), @Const(gdp_lf), @Const(forc_hdm),
        nc4_grass::Int, kmo::Int, kda::Int, tfrz, rpi, lfuel, ufuel,
        cropfire_a1, secsphr)
    p = @index(Global)
    T = eltype(baf_crop)
    @inbounds if mask_soilp[p]
        c = column[p]
        g = gridcell[c]
        if forc_t[c] >= tfrz && itype[p] > nc4_grass &&
           kmo == abm_lf[c] && burndate[p] >= 999 && wtcol[p] > zero(T)
            hdmlf = forc_hdm[g]
            fhd  = T(0.04) + T(0.96) * exp(-one(T) * rpi * (hdmlf / T(350.0))^T(0.5))
            fgdp = T(0.01) + T(0.99) * exp(-one(T) * rpi * (gdp_lf[c] / T(10.0)))
            fb   = smooth_max(zero(T), smooth_min(one(T), (fuelc_crop[c] - lfuel) / (ufuel - lfuel)))
            _scatter_add!(baf_crop, c, cropfire_a1 / secsphr * fhd * fgdp * wtcol[p])
            if fb * fhd * fgdp * wtcol[p] > zero(T)
                burndate[p] = kda
            end
        end
    end
end

# --- Section: main per-column fire-spread loop (Li2016) --------------------
Base.@kwdef struct _Fire2016SpreadOut{V}
    nfire::V; farea_burned::V; fuelc::V; fbac::V; fbac1::V
end
Adapt.@adapt_structure _Fire2016SpreadOut

Base.@kwdef struct _Fire2016SpreadIn1{V}
    forc_hdm::V; cropf_col::V; trotr1_col::V; trotr2_col::V; baf_crop::V
    baf_peatf::V; totlitc::V; totvegc::V; rootc_col::V; fuelc_crop::V
    dzsoi_decomp::V; rh30_col::V; forc_rh::V; wtlf::V; forc_lnfm::V
end
Adapt.@adapt_structure _Fire2016SpreadIn1

Base.@kwdef struct _Fire2016SpreadIn2{V}
    latdeg::V; forc_wind::V; lgdp_col::V; lgdp1_col::V; lpop_col::V
    fsr_col::V; fd_col::V; btran_col::V; tsoi17::V; dtrotr_col::V
    prec60_col::V; prec10_col::V; forc_rain::V; forc_snow::V; lfc::V
    deadstemc_col::V; spinup_lat_term::V
end
Adapt.@adapt_structure _Fire2016SpreadIn2

Base.@kwdef struct _Fire2016SpreadMat{M}
    decomp_cpools_vr::M
end
Adapt.@adapt_structure _Fire2016SpreadMat

Base.@kwdef struct _Fire2016SpreadScalars{T}
    lfuel::T; ufuel::T; rh_low::T; rh_hgh::T; bt_min::T; bt_max::T
    prh30::T; max_rh30_affecting_fuel::T
    pot_hmn_ign_counts_alpha::T; ignition_efficiency::T; g0::T; cli_scale::T
    defo_fire_precip_thresh_bet::T; defo_fire_precip_thresh_bdt::T
    spinup_factor_deadwood::T; spinup_factor_cwd::T
    tfrz::T; rpi::T; secsphr::T; secspday::T; dayspyr::T; dt::T
end

@kernel function _firea2016_fire_spread_kernel!(
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
        bt_min = sc.bt_min; bt_max = sc.bt_max
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
                    fire_m = ((afuel * arh30 + (one(T) - afuel) * arh)^T(1.5)) *
                             ((one(T) - smooth_max(zero(T), smooth_min(one(T),
                                (btran_col[c] / wtlf[c] - bt_min) / (bt_max - bt_min))))^T(0.5))
                else
                    fire_m = zero(T)
                end
                lh = pot_hmn_ign_counts_alpha * T(6.8) * hdmlf^T(0.43) / T(30.0) / T(24.0)
                fs = one(T) - (T(0.01) + T(0.98) * exp(T(-0.025) * hdmlf))
                # Li2016: lightning normalization uses latitude in DEGREES.
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
                        cli = (smooth_max(zero(T), smooth_min(one(T), (cri - prec60_col[c] * secspday) / cri))^T(0.5)) *
                              (smooth_max(zero(T), smooth_min(one(T), (cri - prec10_col[c] * secspday) / cri))^T(0.5)) *
                              smooth_max(T(0.0005), smooth_min(one(T), T(19.0) * dtrotr_col[c] * dayspyr * secspday / dt - T(0.001))) *
                              smooth_max(zero(T), smooth_min(one(T), (T(0.25) - (forc_rain[c] + forc_snow[c]) * secsphr) / T(0.25)))
                        farea_burned[c] = cli * (cli_scale / secspday) + baf_crop[c] + baf_peatf[c]
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

function _firea2016_fire_spread!(
        nfire, farea_burned, fuelc, fbac, fbac1,
        forc_hdm, cropf_col, trotr1_col, trotr2_col, baf_crop, baf_peatf,
        totlitc, totvegc, rootc_col, fuelc_crop, dzsoi_decomp, rh30_col,
        forc_rh, wtlf, forc_lnfm,
        latdeg, forc_wind, lgdp_col, lgdp1_col, lpop_col, fsr_col, fd_col,
        btran_col, tsoi17, dtrotr_col, prec60_col, prec10_col, forc_rain, forc_snow,
        lfc, deadstemc_col, spinup_lat_term,
        decomp_cpools_vr, mask_soilc, gridcell,
        lfuel, ufuel, rh_low, rh_hgh, bt_min, bt_max, prh30, max_rh30_affecting_fuel,
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
        bt_min=T(bt_min), bt_max=T(bt_max), prh30=T(prh30),
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
    _firea2016_fire_spread_kernel!(backend)(out, in1, in2, mat,
        mask_soilc, gridcell, sc,
        i_cwd, nlevdecomp, spinup_state, date_is_jan1_0, transient_landcover;
        ndrange = ndrange)
    KA.synchronize(backend)
    return nothing
end

"""
    cnfire_area_li2016!(fire_li2014, pftcon_li2014, fire_data, cnfire_const,
        cnfire_params, pftcon, masks..., bounds..., patch, col, grc,
        soilstate, h2osoi_vol_col, cnveg_state, cnveg_cs, decomp_cascade_con,
        totlitc_col, decomp_cpools_vr_col, t_soi17cm_col; kwargs...)

Compute column-level burned-area fraction for the Li et al. (2016) fire model.

Ported from `CNFireArea` in `CNFireLi2016Mod.F90`. Re-uses the Li2014-specific
data holder (`CNFireLi2014Data`: forc_hdm/forc_lnfm/gdp/peatf/abm) and PFT fire
parameters (`PftConFireLi2014`: fsr_pft/fd_pft), as the Fortran Li2016 type
reads the same `this%`/`pftcon` members.
"""
function cnfire_area_li2016!(
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
    # --- Aliases (Fortran associate block) ---
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
    # device-resident scratch (similar(totvegc,…) lands on the input's backend).
    btran_col  = fill!(similar(totvegc, FT, nc), zero(FT))
    prec60_col = fill!(similar(totvegc, FT, nc), zero(FT))
    prec10_col = fill!(similar(totvegc, FT, nc), zero(FT))
    rh30_col   = fill!(similar(totvegc, FT, nc), zero(FT))

    # --- pft -> column averaging ---
    p2c!(prec10_col, prec10_patch, patch, mask_soilc, bounds_c, bounds_p)
    p2c!(prec60_col, prec60_patch, patch, mask_soilc, bounds_c, bounds_p)
    p2c!(rh30_col,   rh30_patch,   patch, mask_soilc, bounds_c, bounds_p)
    p2c!(leafc_col,  leafc,        patch, mask_soilc, bounds_c, bounds_p)
    p2c!(deadstemc_col, deadstemc, patch, mask_soilc, bounds_c, bounds_p)

    # --- First time-step: zero area burned and exit ---
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

    # Host-resolved uniform date predicates.
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

    # --- root wetness (btran2), Li2014 method ---
    cnfire_calc_fire_root_wetness_li2014!(
        fire_data, mask_exposedveg, mask_noexposedveg, bounds_p,
        pftcon, patch, soilstate, h2osoi_vol_col, nlevgrnd;
        soil_suction_fn = soil_suction_fn)

    # --- accumulate btran_col / wtlf ---
    _launch!(_firea_btran_wtlf_kernel!, btran_col, wtlf, mask_exposedveg,
             patch.column, patch.itype, patch.wtcol, btran2, cropf_col,
             nc3crop; ndrange = length(mask_exposedveg))

    # --- main noncrop per-patch column-var accumulation (Li2016 form) ---
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

    # --- burned-area fraction in cropland (Li2016 form) ---
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

    # dzsoi_decomp global -> working backend/precision.
    dzsoi_decomp_arr = copyto!(similar(totvegc, FT, length(dzsoi_decomp[])), dzsoi_decomp[])

    # spinup latitude term per column (only consumed by the spinup_state==2 fuel
    # branch; built device-safe & byte-identical to get_spinup_latitude_term).
    spinup_lat_term_col = fill!(similar(totvegc, FT, nc), zero(FT))
    fire_spinup_latterm!(spinup_lat_term_col, mask_soilc, col.gridcell, grc.latdeg)
    spinup_factor_cwd = (i_cwd >= 1 && i_cwd <= length(spinup_factor)) ?
                        spinup_factor[i_cwd] : one(FT)

    # --- main column loop: fractional area affected by fire (Li2016 form) ---
    _firea2016_fire_spread!(
        nfire, farea_burned, fuelc, fbac, fbac1,
        forc_hdm, cropf_col, trotr1_col, trotr2_col, baf_crop, baf_peatf,
        totlitc, totvegc, rootc_col, fuelc_crop, dzsoi_decomp_arr, rh30_col,
        forc_rh, wtlf, forc_lnfm,
        grc.latdeg, forc_wind, lgdp_col, lgdp1_col, lpop_col, fsr_col, fd_col,
        btran_col, tsoi17, dtrotr_col, prec60_col, prec10_col, forc_rain, forc_snow,
        lfc, deadstemc_col, spinup_lat_term_col,
        decomp_cpools_vr, mask_soilc, col.gridcell,
        lfuel, ufuel, rh_low, rh_hgh, bt_min, bt_max, prh30, max_rh30_affecting_fuel,
        pot_hmn_ign_counts_alpha, ignition_efficiency, g0, cli_scale,
        defo_fire_precip_thresh_bet, defo_fire_precip_thresh_bdt,
        spinup_factor_deadwood, spinup_factor_cwd,
        TFRZ, RPI, secsphr, secspday, dayspyr, dt,
        i_cwd, nlevdecomp, spinup_state, date_is_jan1_0, transient_landcover;
        ndrange = length(farea_burned))

    return nothing
end

"""
    cnfire_fluxes_li2016!(...)

Fire effects routine for the Li et al. (2016) fire model. Like all Li-family
methods this inherits the base `CNFireFluxes`; combustion-completeness factors
are 0.5 (litter) / 0.25 (CWD). Delegates to `cnfire_fluxes_li2014!`, which sets
exactly those factors.

Ported from the inherited `CNFireFluxes` under `CNFireLi2016Mod.F90`.
"""
function cnfire_fluxes_li2016!(args...; kwargs...)
    return cnfire_fluxes_li2014!(args...; kwargs...)
end
