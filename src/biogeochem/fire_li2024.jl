# ==========================================================================
# Ported from: src/biogeochem/CNFireLi2024Mod.F90
# Fire dynamics module for Li et al. (2024) version (GSWP/CRUJRA forcing).
#
# Incremental refinement of Li2021. Deltas vs Li2016/2021:
#   * Uses a single 30-day running-mean precip (prec30) instead of prec60+prec10.
#   * lfwt sums ALL natural veg (itype <= nc4_grass) weighted by wtgcell.
#   * Revised agricultural-fire human-density / GDP terms, with crop-live gating
#     (no ag fire while a managed crop is alive when use_crop is on).
#   * Peatland fire: non-boreal branch gated on tropical-tree closed-forest
#     fraction; revised (linear, no fsat) non-boreal formula.
#   * Fire-occurrence ignition `ig` uses lfwt^0.5 (no 1-cropf), and fd_col is
#     rescaled by (lfwt*lgdp1*lpop)^0.5 so the spread term drops the explicit
#     lgdp1/lpop factors.
#   * Deforestation fire accumulates additively onto farea_burned (fb*cli).
#
# Public functions:
#   need_lightning_and_popdens_li2024  — Returns true
#   cnfire_area_li2024!                — Compute column-level burned area
#   cnfire_fluxes_li2024!              — Fire C/N fluxes (delegates to base)
#
# GPU: cnfire_area_li2024! runs entirely through KA kernels on the _launch!
# path (byte-identical on CPU, Metal-validated). Reuses the generic Li2014
# sections, the Li2016 noncrop kernel, and the Li2021 rswf-btran kernel; adds
# Li2024-specific kernels for the wtgcell lfwt sum, the crop-live-gated ag fire,
# the tropical-tree-gated non-boreal peatland, and the fire-spread loop.
# Off-by-default (non-default method).
# ==========================================================================

"""
    need_lightning_and_popdens_li2024()

Returns `true`. Ported from `need_lightning_and_popdens` in `CNFireLi2024Mod.F90`.
"""
need_lightning_and_popdens_li2024() = true

# --- Section: cropf_col + lfwt (Li2024: lfwt sums all natural veg by wtgcell) -
@kernel function _firea2024_cropf_lfwt_kernel!(cropf_col, lfwt, @Const(mask_soilp),
        @Const(column), @Const(itype), @Const(wtcol), @Const(wtgcell), nc4_grass::Int)
    p = @index(Global)
    @inbounds if mask_soilp[p]
        c = column[p]
        if itype[p] > nc4_grass
            _scatter_add!(cropf_col, c, wtcol[p])
        end
        # Li2024: all natural veg (incl. bare soil, itype 0) weighted by wtgcell.
        if itype[p] <= nc4_grass
            _scatter_add!(lfwt, c, wtgcell[p])
        end
    end
end

# --- Section: crop burned-area fraction (Li2024 ag-fire + crop-live gate) ---
@kernel function _firea2024_crop_fire_kernel!(baf_crop, burndate, @Const(mask_soilp),
        @Const(column), @Const(gridcell), @Const(itype), @Const(wtcol),
        @Const(forc_t), @Const(abm_lf), @Const(gdp_lf), @Const(forc_hdm),
        @Const(croplive_patch),
        nc4_grass::Int, kmo::Int, kda::Int, tfrz, rpi, cropfire_a1, secsphr,
        use_crop::Bool, has_croplive::Bool)
    p = @index(Global)
    T = eltype(baf_crop)
    @inbounds if mask_soilp[p]
        c = column[p]
        g = gridcell[c]
        if forc_t[c] >= tfrz && itype[p] > nc4_grass &&
           kmo == abm_lf[c] && burndate[p] >= 999 && wtcol[p] > zero(T)
            hdmlf = forc_hdm[g]
            fhd  = T(0.2) + T(0.8) * exp(-one(T) * rpi * (hdmlf / T(400.0)))
            fgdp = T(0.05) + T(0.95) * exp(-one(T) * rpi * (gdp_lf[c] / T(20.0)))
            croplive = has_croplive ? croplive_patch[p] : false
            if (use_crop && !croplive) || (!use_crop)
                burndate[p] = kda
                _scatter_add!(baf_crop, c, cropfire_a1 / secsphr * fhd * fgdp * wtcol[p])
            end
        end
    end
end

# --- Section: peatland fire (Li2024: tropical-tree-gated non-boreal branch) -
@kernel function _firea2024_peatland_kernel!(baf_peatf, @Const(mask_soilc),
        @Const(gridcell), @Const(latdeg), @Const(wtgcell_col),
        @Const(trotr1_col), @Const(trotr2_col), @Const(prec30_col),
        @Const(peatf_lf), @Const(fsat), @Const(wf2), @Const(tsoi17),
        borealat, non_boreal_peatfire_c, boreal_peatfire_c,
        nonborpeat_fire_precip_denom, borpeat_fire_soilmoist_denom,
        secsphr, secspday, tfrz, rpi)
    c = @index(Global)
    T = eltype(baf_peatf)
    @inbounds if mask_soilc[c]
        g = gridcell[c]
        if latdeg[g] < borealat
            if (trotr1_col[c] + trotr2_col[c]) * wtgcell_col[c] <= T(0.8) &&
               trotr1_col[c] + trotr2_col[c] > zero(T)
                baf_peatf[c] = non_boreal_peatfire_c / secsphr *
                    smooth_max(zero(T), smooth_min(one(T),
                        one(T) - prec30_col[c] * secspday / nonborpeat_fire_precip_denom)) *
                    peatf_lf[c]
            else
                baf_peatf[c] = zero(T)
            end
        else
            baf_peatf[c] = boreal_peatfire_c / secsphr *
                exp(-rpi * (smooth_max(wf2[c], zero(T)) / borpeat_fire_soilmoist_denom)) *
                smooth_max(zero(T), smooth_min(one(T), (tsoi17[c] - tfrz) / T(10.0))) *
                peatf_lf[c] * (one(T) - fsat[c])
        end
    end
end

# --- Section: main per-column fire-spread loop (Li2024) --------------------
Base.@kwdef struct _Fire2024SpreadOut{V}
    nfire::V; farea_burned::V; fuelc::V; fbac::V; fbac1::V; fd_col::V
end
Adapt.@adapt_structure _Fire2024SpreadOut

Base.@kwdef struct _Fire2024SpreadIn1{V}
    forc_hdm::V; cropf_col::V; trotr1_col::V; trotr2_col::V; baf_crop::V
    baf_peatf::V; totlitc::V; totvegc::V; rootc_col::V; fuelc_crop::V
    dzsoi_decomp::V; rh30_col::V; forc_rh::V; wtlf::V; forc_lnfm::V
end
Adapt.@adapt_structure _Fire2024SpreadIn1

Base.@kwdef struct _Fire2024SpreadIn2{V}
    latdeg::V; forc_wind::V; lgdp_col::V; lgdp1_col::V; lpop_col::V
    fsr_col::V; btran_col::V; tsoi17::V; dtrotr_col::V; prec30_col::V
    forc_rain::V; forc_snow::V; lfc::V; lfwt::V; deadstemc_col::V; spinup_lat_term::V
end
Adapt.@adapt_structure _Fire2024SpreadIn2

Base.@kwdef struct _Fire2024SpreadMat{M}
    decomp_cpools_vr::M
end
Adapt.@adapt_structure _Fire2024SpreadMat

Base.@kwdef struct _Fire2024SpreadScalars{T}
    lfuel::T; ufuel::T; rh_low::T; rh_hgh::T; prh30::T; max_rh30_affecting_fuel::T
    pot_hmn_ign_counts_alpha::T; ignition_efficiency::T; g0::T; cli_scale::T
    defo_fire_precip_thresh_bet::T; defo_fire_precip_thresh_bdt::T
    spinup_factor_deadwood::T; spinup_factor_cwd::T
    tfrz::T; rpi::T; secsphr::T; secspday::T; dayspyr::T; dt::T
end

@kernel function _firea2024_fire_spread_kernel!(
        out::_Fire2024SpreadOut, in1::_Fire2024SpreadIn1, in2::_Fire2024SpreadIn2,
        mat::_Fire2024SpreadMat,
        @Const(mask_soilc), @Const(gridcell),
        sc::_Fire2024SpreadScalars,
        i_cwd::Int, nlevdecomp::Int, spinup_state::Int,
        date_is_jan1_0::Bool, transient_landcover::Bool)
    c = @index(Global)
    @inbounds if mask_soilc[c]
        T = eltype(out.nfire)

        nfire = out.nfire; farea_burned = out.farea_burned; fuelc = out.fuelc
        fbac  = out.fbac;  fbac1 = out.fbac1; fd_col = out.fd_col

        forc_hdm   = in1.forc_hdm;   cropf_col  = in1.cropf_col
        trotr1_col = in1.trotr1_col; trotr2_col = in1.trotr2_col
        baf_crop   = in1.baf_crop;   baf_peatf  = in1.baf_peatf
        totlitc    = in1.totlitc;    totvegc    = in1.totvegc
        rootc_col  = in1.rootc_col;  fuelc_crop = in1.fuelc_crop
        dzsoi_decomp = in1.dzsoi_decomp; rh30_col = in1.rh30_col
        forc_rh    = in1.forc_rh;    wtlf = in1.wtlf; forc_lnfm = in1.forc_lnfm

        latdeg    = in2.latdeg;    forc_wind = in2.forc_wind
        lgdp_col  = in2.lgdp_col;  lgdp1_col = in2.lgdp1_col; lpop_col = in2.lpop_col
        fsr_col   = in2.fsr_col;   btran_col = in2.btran_col; tsoi17 = in2.tsoi17
        dtrotr_col = in2.dtrotr_col; prec30_col = in2.prec30_col
        forc_rain = in2.forc_rain; forc_snow = in2.forc_snow; lfc = in2.lfc
        lfwt = in2.lfwt; deadstemc_col = in2.deadstemc_col; spinup_lat_term = in2.spinup_lat_term

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
            # Li2024: fire_m terms computed unconditionally (not gated on trotr).
            afuel = smooth_min(one(T), smooth_max(zero(T), (fc - T(2500.0)) / (T(5000.0) - T(2500.0))))
            arh   = one(T) - smooth_max(zero(T), smooth_min(one(T), (forc_rh[g] - rh_low) / (rh_hgh - rh_low)))
            arh30 = one(T) - smooth_max(prh30, smooth_min(one(T), rh30_col[c] / max_rh30_affecting_fuel))
            if forc_rh[g] < rh_hgh && wtlf[c] > zero(T) && tsoi17[c] > tfrz
                fire_m = ((afuel * arh30 + (one(T) - afuel) * arh)^T(1.5)) *
                         ((one(T) - btran_col[c] / wtlf[c])^T(0.5))
            else
                fire_m = zero(T)
            end
            lh = pot_hmn_ign_counts_alpha * T(6.8) * hdmlf^T(0.43) / T(30.0) / T(24.0)
            fs = one(T) - (T(0.01) + T(0.98) * exp(T(-0.025) * hdmlf))
            # Li2024: ignition uses lfwt^0.5; anthropogenic term only for
            # non-closed-forest columns; lightning normalization in degrees.
            lnorm = forc_lnfm[g] /
                (T(5.16) + T(2.16) * cos(rpi / T(180.0) * T(3) * smooth_min(T(60.0), abs(latdeg[g])))) *
                ignition_efficiency
            if trotr1_col[c] + trotr2_col[c] <= T(0.6)
                ig = (lh + lnorm) * (one(T) - fs) * (lfwt[c]^T(0.5))
            else
                ig = lnorm * (one(T) - fs) * (lfwt[c]^T(0.5))
            end
            nfire[c] = ig / secsphr * fb * fire_m * lgdp_col[c]
            Lb_lf = one(T) + T(10.0) * (one(T) - exp(T(-0.06) * forc_wind[g]))
            spread_m = fire_m^T(0.5)
            # Li2024: fd_col folds in (lfwt*lgdp1*lpop)^0.5; farea drops explicit
            # lgdp1/lpop factors of earlier versions.
            fd_col[c] = (lfwt[c] * lgdp1_col[c] * lpop_col[c])^T(0.5) * fd_col[c]
            farea_burned[c] = smooth_min(one(T), (g0 * spread_m * fsr_col[c] *
                fd_col[c] / T(1000.0))^2 * nfire[c] * rpi * Lb_lf +
                baf_crop[c] + baf_peatf[c])

            if transient_landcover
                if trotr1_col[c] + trotr2_col[c] > T(0.6)
                    if date_is_jan1_0 || dtrotr_col[c] <= zero(T)
                        fbac1[c] = zero(T)
                    else
                        cri = (defo_fire_precip_thresh_bet * trotr1_col[c] +
                               defo_fire_precip_thresh_bdt * trotr2_col[c]) /
                              (trotr1_col[c] + trotr2_col[c])
                        cli = smooth_max(zero(T), smooth_min(one(T), one(T) - prec30_col[c] * secspday / cri)) *
                              (T(15.0) * smooth_min(T(0.0016), dtrotr_col[c] / dt * dayspyr * secspday) + T(0.009)) *
                              smooth_max(zero(T), smooth_min(one(T), (T(0.25) - (forc_rain[c] + forc_snow[c]) * secsphr) / T(0.25)))
                        farea_burned[c] = farea_burned[c] + fb * cli * (cli_scale / secspday)
                        fbac1[c] = smooth_max(zero(T), fb * cli * (cli_scale / secspday) - T(2) * lfc[c] / dt)
                    end
                    fbac[c] = farea_burned[c] + fbac1[c]
                else
                    fbac[c] = farea_burned[c]
                end
            end
        else
            farea_burned[c] = smooth_min(one(T), baf_crop[c] + baf_peatf[c])
        end
    end
end

function _firea2024_fire_spread!(
        nfire, farea_burned, fuelc, fbac, fbac1, fd_col,
        forc_hdm, cropf_col, trotr1_col, trotr2_col, baf_crop, baf_peatf,
        totlitc, totvegc, rootc_col, fuelc_crop, dzsoi_decomp, rh30_col,
        forc_rh, wtlf, forc_lnfm,
        latdeg, forc_wind, lgdp_col, lgdp1_col, lpop_col, fsr_col,
        btran_col, tsoi17, dtrotr_col, prec30_col, forc_rain, forc_snow,
        lfc, lfwt, deadstemc_col, spinup_lat_term,
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
    out = _Fire2024SpreadOut(; nfire, farea_burned, fuelc, fbac, fbac1, fd_col)
    in1 = _Fire2024SpreadIn1(; forc_hdm, cropf_col, trotr1_col, trotr2_col, baf_crop,
        baf_peatf, totlitc, totvegc, rootc_col, fuelc_crop, dzsoi_decomp, rh30_col,
        forc_rh, wtlf, forc_lnfm)
    in2 = _Fire2024SpreadIn2(; latdeg, forc_wind, lgdp_col, lgdp1_col, lpop_col,
        fsr_col, btran_col, tsoi17, dtrotr_col, prec30_col, forc_rain, forc_snow,
        lfc, lfwt, deadstemc_col, spinup_lat_term)
    mat = _Fire2024SpreadMat(; decomp_cpools_vr)
    sc = _Fire2024SpreadScalars(;
        lfuel=T(lfuel), ufuel=T(ufuel), rh_low=T(rh_low), rh_hgh=T(rh_hgh),
        prh30=T(prh30), max_rh30_affecting_fuel=T(max_rh30_affecting_fuel),
        pot_hmn_ign_counts_alpha=T(pot_hmn_ign_counts_alpha),
        ignition_efficiency=T(ignition_efficiency), g0=T(g0), cli_scale=T(cli_scale),
        defo_fire_precip_thresh_bet=T(defo_fire_precip_thresh_bet),
        defo_fire_precip_thresh_bdt=T(defo_fire_precip_thresh_bdt),
        spinup_factor_deadwood=T(spinup_factor_deadwood),
        spinup_factor_cwd=T(spinup_factor_cwd),
        tfrz=T(tfrz), rpi=T(rpi), secsphr=T(secsphr), secspday=T(secspday),
        dayspyr=T(dayspyr), dt=T(dt))

    backend = _kernel_backend(out.nfire)
    _firea2024_fire_spread_kernel!(backend)(out, in1, in2, mat,
        mask_soilc, gridcell, sc,
        i_cwd, nlevdecomp, spinup_state, date_is_jan1_0, transient_landcover;
        ndrange = ndrange)
    KA.synchronize(backend)
    return nothing
end

"""
    cnfire_area_li2024!(...)

Compute column-level burned-area fraction for the Li et al. (2024) fire model.
Ported from `CNFireArea` in `CNFireLi2024Mod.F90`.
"""
function cnfire_area_li2024!(
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
    prec30_patch::AbstractVector{<:Real} = Float64[],
    rh30_patch::AbstractVector{<:Real} = Float64[],
    fsat_col::AbstractVector{<:Real} = Float64[],
    wf2_col::AbstractVector{<:Real} = Float64[],
    croplive_patch::AbstractVector{Bool} = Bool[],
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
    use_crop::Bool = false,
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
    prec30_col = fill!(similar(totvegc, FT, nc), zero(FT))
    rh30_col   = fill!(similar(totvegc, FT, nc), zero(FT))

    # Li2024: a single 30-day precip running mean replaces prec60+prec10.
    p2c!(prec30_col, prec30_patch, patch, mask_soilc, bounds_c, bounds_p)
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

    # --- cropf_col + lfwt (Li2024: lfwt sums all natural veg by wtgcell) ---
    _launch!(_firea_zero_cropf_lfwt_kernel!, cropf_col, lfwt, mask_soilc)
    _launch!(_firea2024_cropf_lfwt_kernel!, cropf_col, lfwt, mask_soilp,
             patch.column, patch.itype, patch.wtcol, patch.wtgcell,
             nc4_grass; ndrange = length(mask_soilp))

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

    # --- accumulate btran_col / wtlf with per-PFT rswf rescale (shared Li2021) ---
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

    # --- burned-area fraction in cropland (Li2024 form) ---
    _launch!(_firea_zero_baf_crop_kernel!, baf_crop, mask_soilc)
    _launch!(_firea_burndate_init_kernel!, burndate, mask_soilp,
             date_is_jan1_0; ndrange = length(mask_soilp))
    _launch!(_firea2024_crop_fire_kernel!, baf_crop, burndate, mask_soilp,
             patch.column, col.gridcell, patch.itype, patch.wtcol,
             forc_t, abm_lf, gdp_lf, forc_hdm, croplive_patch,
             nc4_grass, kmo, kda,
             eltype(baf_crop)(TFRZ), eltype(baf_crop)(RPI),
             eltype(baf_crop)(cropfire_a1), eltype(baf_crop)(secsphr),
             use_crop, !isempty(croplive_patch); ndrange = length(mask_soilp))

    # --- peatland fire (Li2024: tropical-tree-gated non-boreal) ---
    let TFp = eltype(baf_peatf)
        _launch!(_firea2024_peatland_kernel!, baf_peatf, mask_soilc,
             col.gridcell, grc.latdeg, col.wtgcell,
             trotr1_col, trotr2_col, prec30_col, peatf_lf, fsat, wf2, tsoi17,
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

    # --- main column loop: fractional area affected by fire (Li2024 form) ---
    _firea2024_fire_spread!(
        nfire, farea_burned, fuelc, fbac, fbac1, fd_col,
        forc_hdm, cropf_col, trotr1_col, trotr2_col, baf_crop, baf_peatf,
        totlitc, totvegc, rootc_col, fuelc_crop, dzsoi_decomp_arr, rh30_col,
        forc_rh, wtlf, forc_lnfm,
        grc.latdeg, forc_wind, lgdp_col, lgdp1_col, lpop_col, fsr_col,
        btran_col, tsoi17, dtrotr_col, prec30_col, forc_rain, forc_snow,
        lfc, lfwt, deadstemc_col, spinup_lat_term_col,
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
    cnfire_fluxes_li2024!(...)

Fire effects routine for the Li et al. (2024) fire model. Inherits the base
`CNFireFluxes` (litter 0.5 / CWD 0.25); delegates to `cnfire_fluxes_li2014!`.
"""
function cnfire_fluxes_li2024!(args...; kwargs...)
    return cnfire_fluxes_li2014!(args...; kwargs...)
end
