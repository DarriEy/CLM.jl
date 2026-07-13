# ==========================================================================
# snow_snicar_device.jl — on-device (KernelAbstractions/Metal) port of
# `snicar_rt!` from snow_snicar.jl.
#
# The host `snicar_rt!` is a single per-column loop with ~40 reused scratch
# arrays, lookup-table interpolation and a full Adding-Doubling RT solver. It is
# kept intact for the CPU path and for AD. Here we provide `snicar_rt_device!`,
# which runs the same physics one-thread-per-column on whatever backend the
# output arrays live on.
#
# Decomposition (three kernels sharing per-column scratch in global memory):
#   _snicar_setup_kernel!    — zero flx_abs output, build per-column snow-layer
#                              local arrays (handles the no-snow-layer case)
#   _snicar_band_kernel!     — launched once per SNICAR spectral band: snow +
#                              aerosol optics, delta scaling, Adding-Doubling RT
#   _snicar_finalize_kernel! — band-weight albout/flx_abs into VIS/NIR, SZA fix
#
# Scratch persists across the launches because it is allocated once (per column)
# in device global memory. The spectral-band loop runs on the host (5 launches).
#
# Differences from the host that are byte-identical for valid inputs:
#   - `sno_shp` (a per-layer String) collapses to one Int shape code (the config
#     is uniform across layers); `sno_fs`/`sno_AR` are never set in this path so
#     the shape defaults always apply.
#   - The grain-size bounds `error()`/`@warn` is dropped (surface_albedo clamps
#     snw_rds into [SNW_RDS_MIN_TBL, SNW_RDS_MAX_TBL] before the call, so the
#     check never fires); dead host locals (albsfc_lcl, flx_slr*_lcl, F_sfc_pls,
#     frac_sno) are omitted.
#   - Float literals are written `T(...)` so a Float32 backend carries no Float64.
# ==========================================================================

# Device-inline 1-D piecewise linear interpolation over a const vector `xd`
# (length nd) and the `c`-th row of scratch matrix `ym`. Mirrors
# piecewise_linear_interp1d but takes a (matrix,row) instead of a slice so no
# view is allocated inside the kernel.
@inline function _snicar_pw_interp(xd, ym, c, nd, xi)
    T = eltype(ym)
    if nd == 1
        return ym[c, 1]
    end
    if xi < xd[1]
        t = (xi - xd[1]) / (xd[2] - xd[1])
        return (one(T) - t) * ym[c, 1] + t * ym[c, 2]
    elseif xi > xd[nd]
        t = (xi - xd[nd-1]) / (xd[nd] - xd[nd-1])
        return (one(T) - t) * ym[c, nd-1] + t * ym[c, nd]
    else
        res = zero(T)
        for k in 2:nd
            if xd[k-1] <= xi && xi <= xd[k]
                t = (xi - xd[k-1]) / (xd[k] - xd[k-1])
                res = (one(T) - t) * ym[c, k-1] + t * ym[c, k]
                break
            end
        end
        return res
    end
end

# ---- Device-resident bundles (grouped to keep the kernel arg lists small) ----

# Snow optics lookup tables, copied onto the working backend at element type T.
Base.@kwdef struct _SnicarOptics{V,M}
    ss_alb_snw_drc::M; asm_prm_snw_drc::M; ext_cff_mss_snw_drc::M
    ss_alb_snw_dfs::M; asm_prm_snw_dfs::M; ext_cff_mss_snw_dfs::M
    ss_alb_bc_hphil::V; asm_prm_bc_hphil::V; ext_cff_mss_bc_hphil::V
    ss_alb_bc_hphob::V; asm_prm_bc_hphob::V; ext_cff_mss_bc_hphob::V
    ss_alb_oc_hphil::V; asm_prm_oc_hphil::V; ext_cff_mss_oc_hphil::V
    ss_alb_oc_hphob::V; asm_prm_oc_hphob::V; ext_cff_mss_oc_hphob::V
    ss_alb_dst1::V; asm_prm_dst1::V; ext_cff_mss_dst1::V
    ss_alb_dst2::V; asm_prm_dst2::V; ext_cff_mss_dst2::V
    ss_alb_dst3::V; asm_prm_dst3::V; ext_cff_mss_dst3::V
    ss_alb_dst4::V; asm_prm_dst4::V; ext_cff_mss_dst4::V
    flx_wgt_dir::V; flx_wgt_dif::V
end
Adapt.@adapt_structure _SnicarOptics

# Small constant lookup vectors (band centers, Gaussian quadrature, nonspherical
# / internal-mixing coefficients), all on the working backend at element type T.
Base.@kwdef struct _SnicarConst{V}
    wvl_ct::V; g_wvl_ct::V; dstint_wvl_ct::V; bcint_wvl_ct::V
    difgauspt::V; difgauswt::V
    g_b0::V; g_b1::V; g_b2::V
    g_f07_c0::V; g_f07_c1::V; g_f07_c2::V
    g_f07_p0::V; g_f07_p1::V; g_f07_p2::V
    bcint_d0::V; bcint_d1::V; bcint_d2::V; bcint_m::V; bcint_n::V
    dstint_a1::V; dstint_a2::V; dstint_a3::V
end
Adapt.@adapt_structure _SnicarConst

# Per-column scratch (one row per column; reused across bands within a column).
Base.@kwdef struct _SnicarScratch{M,A3,MI}
    h2osno_ice_lcl::M; h2osno_liq_lcl::M; snw_rds_lcl::MI
    mss_cnc_aer_lcl::A3
    ss_alb_snw_lcl::M; asm_prm_snw_lcl::M; ext_cff_mss_snw_lcl::M
    ss_alb_aer_lcl::M; asm_prm_aer_lcl::M; ext_cff_mss_aer_lcl::M
    g_ice_Cg_tmp::M; gg_ice_F07_tmp::M
    L_snw::M; tau_snw::M; L_aer::A3; tau_aer::A3
    tau_arr::M; omega_arr::M; g_arr::M
    enh_omg_bcint_tmp::M; enh_omg_bcint_tmp2::M
    enh_omg_dstint_tmp::M; enh_omg_dstint_tmp2::M
    tau_star::M; omega_star::M; g_star::M
    trndir::M; trntdr::M; trndif::M; rupdir::M; rupdif::M; rdndif::M
    dfdir::M; dfdif::M; dftmp::M
    rdir::M; rdif_a::M; rdif_b::M; tdir_arr::M; tdif_a::M; tdif_b::M
    trnlay::M; F_abs::M
    albout_lcl::M; flx_abs_lcl::A3
end
Adapt.@adapt_structure _SnicarScratch

# State inputs read by the kernels.
Base.@kwdef struct _SnicarIn{V,M,A3,VI,MI}
    coszen::V; h2osno_liq::M; h2osno_ice::M; h2osno_total::V
    snw_rds::MI; mss_cnc_aer_in::A3; albsfc::M; snl::VI
end
Adapt.@adapt_structure _SnicarIn

# Outputs: snow albedo (col,bnd) and absorbed flux (col,lyr,bnd).
Base.@kwdef struct _SnicarOut{M,A3}
    albout::M; flx_abs::A3
end
Adapt.@adapt_structure _SnicarOut

# --------------------------------------------------------------------------
# Setup kernel — zero flx_abs output and fill per-column snow-layer locals.
# --------------------------------------------------------------------------
@kernel function _snicar_setup_kernel!(o::_SnicarOut, s::_SnicarScratch, in::_SnicarIn,
        @Const(mask), nlevsno::Int, snw_rds_min_int::Int)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(s.h2osno_ice_lcl)
        joff = nlevsno
        for i in 1:(nlevsno+1)
            for bnd in 1:NUMRAD
                o.flx_abs[c, i, bnd] = zero(T)
            end
        end
        h2osno_lcl = in.h2osno_total[c]
        if in.coszen[c] > zero(T) && h2osno_lcl > T(MIN_SNW)
            for i in 1:nlevsno
                s.h2osno_ice_lcl[c, i] = zero(T)
                s.h2osno_liq_lcl[c, i] = zero(T)
                s.snw_rds_lcl[c, i] = 0
            end
            if in.snl[c] > -1
                # No resolved snow layers: lump all mass into layer `joff`.
                s.h2osno_ice_lcl[c, joff] = h2osno_lcl
                s.h2osno_liq_lcl[c, joff] = zero(T)
                s.snw_rds_lcl[c, joff] = snw_rds_min_int
            else
                for i in 1:nlevsno
                    s.h2osno_liq_lcl[c, i] = in.h2osno_liq[c, i]
                    s.h2osno_ice_lcl[c, i] = in.h2osno_ice[c, i]
                    s.snw_rds_lcl[c, i] = in.snw_rds[c, i]
                end
            end
        end
    end
end

# --------------------------------------------------------------------------
# Band kernel — snow & aerosol optics + Adding-Doubling RT for one band.
# --------------------------------------------------------------------------
@kernel function _snicar_band_kernel!(o::_SnicarOut, s::_SnicarScratch,
        opt::_SnicarOptics, cst::_SnicarConst, in::_SnicarIn, @Const(mask),
        flg_slr_in::Int, bnd_idx::Int, numrad_snw::Int, nir_bnd_bgn::Int,
        nlevsno::Int, sno_shp_code::Int, snobc::Bool, snodst::Bool)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(s.tau_arr)
        joff = nlevsno
        h2osno_lcl = in.h2osno_total[c]
        if in.coszen[c] > zero(T) && h2osno_lcl > T(MIN_SNW)
            flg_nosnl = in.snl[c] > -1 ? 1 : 0
            snl_lcl = in.snl[c] > -1 ? -1 : Int(in.snl[c])
            snl_top_j = (snl_lcl + 1) + joff
            snl_btm_j = 0 + joff
            snl_btm_itf = snl_btm_j + 1

            c0 = zero(T); c1 = one(T); c3 = T(3.0); c4 = T(4.0)
            cp01 = T(0.01); cp5 = T(0.5); cp75 = T(0.75); c1p5 = T(1.5)
            trmin = T(0.001); exp_min = exp(-T(10.0))

            mu_not = smooth_max(in.coszen[c], cp01)

            # Restore aerosol concentrations (zeroed last band / no-layer case).
            if flg_nosnl == 0
                for j in 1:SNO_NBR_AER
                    for i in 1:nlevsno
                        s.mss_cnc_aer_lcl[c, i, j] = in.mss_cnc_aer_in[c, i, j]
                    end
                end
            else
                for j in 1:SNO_NBR_AER
                    for i in 1:nlevsno
                        s.mss_cnc_aer_lcl[c, i, j] = zero(T)
                    end
                end
            end
            # Zero aerosols for highly absorptive bands.
            if numrad_snw == DEFAULT_NUMBER_BANDS &&
               (bnd_idx == HIGHEST_DEFAULT_BAND || bnd_idx == SEC_HIGHEST_DEFAULT_BAND)
                for j in 1:SNO_NBR_AER
                    for i in 1:nlevsno
                        s.mss_cnc_aer_lcl[c, i, j] = zero(T)
                    end
                end
            end
            if numrad_snw == HIGH_NUMBER_BANDS && bnd_idx > 100
                for j in 1:SNO_NBR_AER
                    for i in 1:nlevsno
                        s.mss_cnc_aer_lcl[c, i, j] = zero(T)
                    end
                end
            end

            # ---- Snow optical properties from lookup table ----
            for i in snl_top_j:snl_btm_j
                rds_idx = Int(s.snw_rds_lcl[c, i]) - SNW_RDS_MIN_TBL + 1
                if flg_slr_in == 1
                    s.ss_alb_snw_lcl[c, i] = opt.ss_alb_snw_drc[rds_idx, bnd_idx]
                    s.ext_cff_mss_snw_lcl[c, i] = opt.ext_cff_mss_snw_drc[rds_idx, bnd_idx]
                    if sno_shp_code == 0
                        s.asm_prm_snw_lcl[c, i] = opt.asm_prm_snw_drc[rds_idx, bnd_idx]
                    end
                else
                    s.ss_alb_snw_lcl[c, i] = opt.ss_alb_snw_dfs[rds_idx, bnd_idx]
                    s.ext_cff_mss_snw_lcl[c, i] = opt.ext_cff_mss_snw_dfs[rds_idx, bnd_idx]
                    if sno_shp_code == 0
                        s.asm_prm_snw_lcl[c, i] = opt.asm_prm_snw_dfs[rds_idx, bnd_idx]
                    end
                end
            end

            # ---- Nonspherical snow asymmetry factor ----
            wvl_b = cst.wvl_ct[bnd_idx]
            for i in snl_top_j:snl_btm_j
                if sno_shp_code != 0
                    rds = T(s.snw_rds_lcl[c, i])
                    if sno_shp_code == 1          # spheroid
                        diam_ice = T(2.0) * rds
                        fs_ratio = T(FS_SPHD_DEFAULT) / T(FS_HEX_REF)
                        AR_tmp = T(AR_TMP_DEFAULT_1)
                        for igb in 1:SNICAR_SEVEN_BANDS
                            s.g_ice_Cg_tmp[c, igb] = cst.g_b0[igb] *
                                (fs_ratio^cst.g_b1[igb]) * (diam_ice^cst.g_b2[igb])
                            s.gg_ice_F07_tmp[c, igb] = cst.g_f07_c0[igb] +
                                cst.g_f07_c1[igb] * AR_tmp + cst.g_f07_c2[igb] * (AR_tmp * AR_tmp)
                        end
                    elseif sno_shp_code == 2      # hexagonal_plate
                        diam_ice = T(2.0) * rds
                        fs_ratio = T(FS_HEX_REF) / T(FS_HEX_REF)
                        AR_tmp = T(AR_TMP_DEFAULT_2)
                        for igb in 1:SNICAR_SEVEN_BANDS
                            s.g_ice_Cg_tmp[c, igb] = cst.g_b0[igb] *
                                (fs_ratio^cst.g_b1[igb]) * (diam_ice^cst.g_b2[igb])
                            s.gg_ice_F07_tmp[c, igb] = cst.g_f07_p0[igb] +
                                cst.g_f07_p1[igb] * log(AR_tmp) + cst.g_f07_p2[igb] * (log(AR_tmp) * log(AR_tmp))
                        end
                    else                          # koch_snowflake
                        diam_ice = T(2.0) * rds / T(0.544)
                        fs_ratio = T(FS_KOCH_DEFAULT) / T(FS_HEX_REF)
                        AR_tmp = T(AR_TMP_DEFAULT_2)
                        for igb in 1:SNICAR_SEVEN_BANDS
                            s.g_ice_Cg_tmp[c, igb] = cst.g_b0[igb] *
                                (fs_ratio^cst.g_b1[igb]) * (diam_ice^cst.g_b2[igb])
                            s.gg_ice_F07_tmp[c, igb] = cst.g_f07_p0[igb] +
                                cst.g_f07_p1[igb] * log(AR_tmp) + cst.g_f07_p2[igb] * (log(AR_tmp) * log(AR_tmp))
                        end
                    end
                    g_Cg_intp = _snicar_pw_interp(cst.g_wvl_ct, s.g_ice_Cg_tmp, c, SNICAR_SEVEN_BANDS, wvl_b)
                    gg_F07_intp = _snicar_pw_interp(cst.g_wvl_ct, s.gg_ice_F07_tmp, c, SNICAR_SEVEN_BANDS, wvl_b)
                    g_ice_F07 = gg_F07_intp + T(0.5) * (one(T) - gg_F07_intp) / s.ss_alb_snw_lcl[c, i]
                    s.asm_prm_snw_lcl[c, i] = g_ice_F07 * g_Cg_intp
                end
                s.asm_prm_snw_lcl[c, i] = smooth_min(T(0.99), s.asm_prm_snw_lcl[c, i])
            end

            # ---- Aerosol optical properties from lookup table ----
            s.ss_alb_aer_lcl[c, 1] = opt.ss_alb_bc_hphil[bnd_idx]
            s.asm_prm_aer_lcl[c, 1] = opt.asm_prm_bc_hphil[bnd_idx]
            s.ext_cff_mss_aer_lcl[c, 1] = opt.ext_cff_mss_bc_hphil[bnd_idx]
            s.ss_alb_aer_lcl[c, 2] = opt.ss_alb_bc_hphob[bnd_idx]
            s.asm_prm_aer_lcl[c, 2] = opt.asm_prm_bc_hphob[bnd_idx]
            s.ext_cff_mss_aer_lcl[c, 2] = opt.ext_cff_mss_bc_hphob[bnd_idx]
            s.ss_alb_aer_lcl[c, 3] = opt.ss_alb_oc_hphil[bnd_idx]
            s.asm_prm_aer_lcl[c, 3] = opt.asm_prm_oc_hphil[bnd_idx]
            s.ext_cff_mss_aer_lcl[c, 3] = opt.ext_cff_mss_oc_hphil[bnd_idx]
            s.ss_alb_aer_lcl[c, 4] = opt.ss_alb_oc_hphob[bnd_idx]
            s.asm_prm_aer_lcl[c, 4] = opt.asm_prm_oc_hphob[bnd_idx]
            s.ext_cff_mss_aer_lcl[c, 4] = opt.ext_cff_mss_oc_hphob[bnd_idx]
            s.ss_alb_aer_lcl[c, 5] = opt.ss_alb_dst1[bnd_idx]
            s.asm_prm_aer_lcl[c, 5] = opt.asm_prm_dst1[bnd_idx]
            s.ext_cff_mss_aer_lcl[c, 5] = opt.ext_cff_mss_dst1[bnd_idx]
            s.ss_alb_aer_lcl[c, 6] = opt.ss_alb_dst2[bnd_idx]
            s.asm_prm_aer_lcl[c, 6] = opt.asm_prm_dst2[bnd_idx]
            s.ext_cff_mss_aer_lcl[c, 6] = opt.ext_cff_mss_dst2[bnd_idx]
            s.ss_alb_aer_lcl[c, 7] = opt.ss_alb_dst3[bnd_idx]
            s.asm_prm_aer_lcl[c, 7] = opt.asm_prm_dst3[bnd_idx]
            s.ext_cff_mss_aer_lcl[c, 7] = opt.ext_cff_mss_dst3[bnd_idx]
            s.ss_alb_aer_lcl[c, 8] = opt.ss_alb_dst4[bnd_idx]
            s.asm_prm_aer_lcl[c, 8] = opt.asm_prm_dst4[bnd_idx]
            s.ext_cff_mss_aer_lcl[c, 8] = opt.ext_cff_mss_dst4[bnd_idx]

            # ---- Layer mass, optical depths, weighted Mie properties ----
            for i in snl_top_j:snl_btm_j
                # BC/dust-snow internal mixing for wavelength <= 1.2um
                if wvl_b <= T(1.2)
                    if snobc && s.mss_cnc_aer_lcl[c, i, 1] > zero(T)
                        for ibb in 1:SNICAR_SIXTEEN_BANDS
                            s.enh_omg_bcint_tmp[c, ibb] = cst.bcint_d0[ibb] *
                                ((s.mss_cnc_aer_lcl[c, i, 1] * T(KG_TO_UG) * T(DEN_BC) / T(DEN_BC_TARGET) +
                                  cst.bcint_d2[ibb])^cst.bcint_d1[ibb])
                            if ibb < 3
                                bcint_m_tmp = cst.bcint_m[1]; bcint_n_tmp = cst.bcint_n[1]
                            elseif ibb <= 11
                                bcint_m_tmp = cst.bcint_m[2]; bcint_n_tmp = cst.bcint_n[2]
                            else
                                bcint_m_tmp = cst.bcint_m[3]; bcint_n_tmp = cst.bcint_n[3]
                            end
                            bcint_dd = (T(RE_BC) / T(RADIUS_2_BC))^bcint_m_tmp
                            bcint_dd2 = (T(RADIUS_1_BC) / T(RADIUS_2_BC))^bcint_m_tmp
                            bcint_f = (T(RE_BC) / T(RADIUS_1_BC))^bcint_n_tmp
                            s.enh_omg_bcint_tmp2[c, ibb] = log10(smooth_max(one(T),
                                bcint_dd * ((s.enh_omg_bcint_tmp[c, ibb] / bcint_dd2)^bcint_f)))
                        end
                        enh_intp = _snicar_pw_interp(cst.bcint_wvl_ct, s.enh_omg_bcint_tmp2, c, SNICAR_SIXTEEN_BANDS, wvl_b)
                        enh2 = T(10.0)^enh_intp
                        enh2 = smooth_min(T(ENH_OMG_MAX), smooth_max(enh2, one(T)))
                        s.ss_alb_snw_lcl[c, i] = one(T) - (one(T) - s.ss_alb_snw_lcl[c, i]) * enh2
                        s.ss_alb_snw_lcl[c, i] = smooth_max(T(0.5), smooth_min(s.ss_alb_snw_lcl[c, i], one(T)))
                        s.ss_alb_aer_lcl[c, 1] = zero(T)
                        s.asm_prm_aer_lcl[c, 1] = zero(T)
                        s.ext_cff_mss_aer_lcl[c, 1] = zero(T)
                    end
                    tot_dst = (s.mss_cnc_aer_lcl[c, i, 5] + s.mss_cnc_aer_lcl[c, i, 6] +
                               s.mss_cnc_aer_lcl[c, i, 7] + s.mss_cnc_aer_lcl[c, i, 8]) * T(KG_KG_TO_PPM)
                    if snodst && tot_dst > zero(T)
                        for idb in 1:SNICAR_SIZE_BINS
                            s.enh_omg_dstint_tmp[c, idb] = cst.dstint_a1[idb] +
                                cst.dstint_a2[idb] * (tot_dst^cst.dstint_a3[idb])
                            s.enh_omg_dstint_tmp2[c, idb] = log10(smooth_max(s.enh_omg_dstint_tmp[c, idb], one(T)))
                        end
                        enh_intp = _snicar_pw_interp(cst.dstint_wvl_ct, s.enh_omg_dstint_tmp2, c, SNICAR_SIZE_BINS, wvl_b)
                        enh2 = T(10.0)^enh_intp
                        enh2 = smooth_min(T(ENH_OMG_MAX), smooth_max(enh2, one(T)))
                        s.ss_alb_snw_lcl[c, i] = one(T) - (one(T) - s.ss_alb_snw_lcl[c, i]) * enh2
                        s.ss_alb_snw_lcl[c, i] = smooth_max(T(0.5), smooth_min(s.ss_alb_snw_lcl[c, i], one(T)))
                        for k in 5:8
                            s.ss_alb_aer_lcl[c, k] = zero(T)
                            s.asm_prm_aer_lcl[c, k] = zero(T)
                            s.ext_cff_mss_aer_lcl[c, k] = zero(T)
                        end
                    end
                end

                s.L_snw[c, i] = s.h2osno_ice_lcl[c, i] + s.h2osno_liq_lcl[c, i]
                s.tau_snw[c, i] = s.L_snw[c, i] * s.ext_cff_mss_snw_lcl[c, i]

                for j in 1:SNO_NBR_AER
                    s.L_aer[c, i, j] = s.L_snw[c, i] * s.mss_cnc_aer_lcl[c, i, j]
                    s.tau_aer[c, i, j] = s.L_aer[c, i, j] * s.ext_cff_mss_aer_lcl[c, j]
                end

                tau_sum = zero(T); omega_sum = zero(T); g_sum_val = zero(T)
                for j in 1:SNO_NBR_AER
                    tau_sum += s.tau_aer[c, i, j]
                    omega_sum += s.tau_aer[c, i, j] * s.ss_alb_aer_lcl[c, j]
                    g_sum_val += s.tau_aer[c, i, j] * s.ss_alb_aer_lcl[c, j] * s.asm_prm_aer_lcl[c, j]
                end

                s.tau_arr[c, i] = tau_sum + s.tau_snw[c, i]
                s.omega_arr[c, i] = (one(T) / s.tau_arr[c, i]) * (omega_sum + s.ss_alb_snw_lcl[c, i] * s.tau_snw[c, i])
                s.g_arr[c, i] = (one(T) / (s.tau_arr[c, i] * s.omega_arr[c, i])) *
                    (g_sum_val + s.asm_prm_snw_lcl[c, i] * s.ss_alb_snw_lcl[c, i] * s.tau_snw[c, i])
            end

            # ---- Delta transformation (DELTA always on) ----
            for i in snl_top_j:snl_btm_j
                g = s.g_arr[c, i]; om = s.omega_arr[c, i]; ta = s.tau_arr[c, i]
                s.g_star[c, i] = g / (one(T) + g)
                s.omega_star[c, i] = (one(T) - g * g) * om / (one(T) - om * g * g)
                s.tau_star[c, i] = (one(T) - om * g * g) * ta
            end

            # ---- Adding-Doubling RT solver ----
            for i in 1:(nlevsno+1)
                s.trndir[c, i] = c0; s.trntdr[c, i] = c0; s.trndif[c, i] = c0
                s.rupdir[c, i] = c0; s.rupdif[c, i] = c0; s.rdndif[c, i] = c0
                s.dfdir[c, i] = c0; s.dfdif[c, i] = c0; s.dftmp[c, i] = c0
            end
            for i in 1:nlevsno
                s.rdir[c, i] = c0; s.rdif_a[c, i] = c0; s.rdif_b[c, i] = c0
                s.tdir_arr[c, i] = c0; s.tdif_a[c, i] = c0; s.tdif_b[c, i] = c0
                s.trnlay[c, i] = c0; s.F_abs[c, i] = c0
            end
            s.trndir[c, snl_top_j] = c1
            s.trntdr[c, snl_top_j] = c1
            s.trndif[c, snl_top_j] = c1

            for i in snl_top_j:snl_btm_j
                if s.trntdr[c, i] > trmin
                    ts = s.tau_star[c, i]; ws = s.omega_star[c, i]; gs = s.g_star[c, i]
                    lm = sqrt(c3 * (c1 - ws) * (c1 - ws * gs))
                    ue = c1p5 * (c1 - ws * gs) / lm
                    extins = max(exp_min, exp(-lm * ts))
                    ne_val = ((ue + c1)^2 / extins) - ((ue - c1)^2 * extins)
                    s.rdif_a[c, i] = (ue^2 - c1) * (one(T) / extins - extins) / ne_val
                    s.tdif_a[c, i] = c4 * ue / ne_val
                    s.trnlay[c, i] = max(exp_min, exp(-ts / mu_not))
                    alp = cp75 * ws * mu_not * ((c1 + gs * (c1 - ws)) / (c1 - lm^2 * mu_not^2))
                    gam = cp5 * ws * ((c1 + c3 * gs * (c1 - ws) * mu_not^2) / (c1 - lm^2 * mu_not^2))
                    apg = alp + gam; amg = alp - gam
                    s.rdir[c, i] = apg * s.rdif_a[c, i] + amg * (s.tdif_a[c, i] * s.trnlay[c, i] - c1)
                    s.tdir_arr[c, i] = apg * s.tdif_a[c, i] + (amg * s.rdif_a[c, i] - apg + c1) * s.trnlay[c, i]

                    R1 = s.rdif_a[c, i]; T1 = s.tdif_a[c, i]
                    swt = c0; smr = c0; smt = c0
                    for ng in 1:SNICAR_NGMAX
                        mu = cst.difgauspt[ng]; gwt = cst.difgauswt[ng]
                        swt += mu * gwt
                        trn = max(exp_min, exp(-ts / mu))
                        alpg = cp75 * ws * mu * ((c1 + gs * (c1 - ws)) / (c1 - lm^2 * mu^2))
                        gamg = cp5 * ws * ((c1 + c3 * gs * (c1 - ws) * mu^2) / (c1 - lm^2 * mu^2))
                        apgg = alpg + gamg; amgg = alpg - gamg
                        rdr = apgg * R1 + amgg * T1 * trn - amgg
                        tdr = apgg * T1 + amgg * R1 * trn - apgg * trn + trn
                        smr += mu * rdr * gwt
                        smt += mu * tdr * gwt
                    end
                    s.rdif_a[c, i] = smr / swt
                    s.tdif_a[c, i] = smt / swt
                    s.rdif_b[c, i] = s.rdif_a[c, i]
                    s.tdif_b[c, i] = s.tdif_a[c, i]
                end

                s.trndir[c, i+1] = s.trndir[c, i] * s.trnlay[c, i]
                refkm1 = c1 / (c1 - s.rdndif[c, i] * s.rdif_a[c, i])
                tdrrdir = s.trndir[c, i] * s.rdir[c, i]
                tdndif = s.trntdr[c, i] - s.trndir[c, i]
                s.trntdr[c, i+1] = s.trndir[c, i] * s.tdir_arr[c, i] +
                                   (tdndif + tdrrdir * s.rdndif[c, i]) * refkm1 * s.tdif_a[c, i]
                s.rdndif[c, i+1] = s.rdif_b[c, i] + s.tdif_b[c, i] * s.rdndif[c, i] * refkm1 * s.tdif_a[c, i]
                s.trndif[c, i+1] = s.trndif[c, i] * refkm1 * s.tdif_a[c, i]
            end

            # Bottom boundary: underlying ground albedo.
            albg = bnd_idx < nir_bnd_bgn ? in.albsfc[c, IVIS] : in.albsfc[c, INIR]
            s.rupdir[c, snl_btm_itf] = albg
            s.rupdif[c, snl_btm_itf] = albg

            # Upward sweep.
            for i in snl_btm_j:-1:snl_top_j
                refkp1 = c1 / (c1 - s.rdif_b[c, i] * s.rupdif[c, i+1])
                s.rupdir[c, i] = s.rdir[c, i] +
                    (s.trnlay[c, i] * s.rupdir[c, i+1] +
                     (s.tdir_arr[c, i] - s.trnlay[c, i]) * s.rupdif[c, i+1]) * refkp1 * s.tdif_b[c, i]
                s.rupdif[c, i] = s.rdif_a[c, i] + s.tdif_a[c, i] * s.rupdif[c, i+1] * refkp1 * s.tdif_b[c, i]
            end

            # Net flux at each interface.
            for i in snl_top_j:snl_btm_itf
                refk = c1 / (c1 - s.rdndif[c, i] * s.rupdif[c, i])
                s.dfdir[c, i] = s.trndir[c, i] +
                    (s.trntdr[c, i] - s.trndir[c, i]) * (c1 - s.rupdif[c, i]) * refk -
                    s.trndir[c, i] * s.rupdir[c, i] * (c1 - s.rdndif[c, i]) * refk
                if s.dfdir[c, i] < T(SNICAR_PUNY); s.dfdir[c, i] = c0; end
                s.dfdif[c, i] = s.trndif[c, i] * (c1 - s.rupdif[c, i]) * refk
                if s.dfdif[c, i] < T(SNICAR_PUNY); s.dfdif[c, i] = c0; end
            end

            # Select direct/diffuse.
            if flg_slr_in == 1
                albedo = s.rupdir[c, snl_top_j]
                for i in snl_top_j:snl_btm_itf; s.dftmp[c, i] = s.dfdir[c, i]; end
            else
                albedo = s.rupdif[c, snl_top_j]
                for i in snl_top_j:snl_btm_itf; s.dftmp[c, i] = s.dfdif[c, i]; end
            end

            # Absorbed flux in each layer.
            for i in snl_top_j:snl_btm_j
                s.F_abs[c, i] = s.dftmp[c, i] - s.dftmp[c, i+1]
                s.flx_abs_lcl[c, i, bnd_idx] = s.F_abs[c, i]
            end
            F_btm_net = s.dftmp[c, snl_btm_itf]
            s.flx_abs_lcl[c, snl_btm_itf, bnd_idx] = F_btm_net
            if flg_nosnl == 1
                s.flx_abs_lcl[c, snl_btm_j, bnd_idx] = s.F_abs[c, snl_btm_j]
                s.flx_abs_lcl[c, snl_btm_itf, bnd_idx] = F_btm_net
            end
            for i in snl_top_j:snl_btm_itf
                s.flx_abs_lcl[c, i, bnd_idx] = smooth_max(c0, s.flx_abs_lcl[c, i, bnd_idx])
            end

            s.albout_lcl[c, bnd_idx] = albedo
        end
    end
end

# --------------------------------------------------------------------------
# Finalize kernel — band-weight albout/flx_abs into VIS/NIR + SZA adjustment.
# --------------------------------------------------------------------------
@kernel function _snicar_finalize_kernel!(o::_SnicarOut, s::_SnicarScratch, opt::_SnicarOptics,
        in::_SnicarIn, @Const(mask), flg_slr_in::Int, numrad_snw::Int,
        nir_bnd_bgn::Int, nlevsno::Int)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(s.tau_arr)
        joff = nlevsno
        nir_bnd_end = numrad_snw
        h2osno_lcl = in.h2osno_total[c]
        coszc = in.coszen[c]
        if coszc > zero(T) && h2osno_lcl > T(MIN_SNW)
            snl_lcl = in.snl[c] > -1 ? -1 : Int(in.snl[c])
            snl_top_j = (snl_lcl + 1) + joff
            mu_not = smooth_max(coszc, T(0.01))

            # Sum of flux weights over NIR bands.
            nir_sum = zero(T)
            for b in nir_bnd_bgn:nir_bnd_end
                nir_sum += (flg_slr_in == 1 ? opt.flx_wgt_dir[b] : opt.flx_wgt_dif[b])
            end

            # Output albedo into VIS band.
            if numrad_snw == DEFAULT_NUMBER_BANDS
                o.albout[c, IVIS] = s.albout_lcl[c, IVIS]
            elseif numrad_snw == HIGH_NUMBER_BANDS
                vis_flx = zero(T); vis_wsum = zero(T)
                for b in 1:(nir_bnd_bgn-1)
                    w = (flg_slr_in == 1 ? opt.flx_wgt_dir[b] : opt.flx_wgt_dif[b])
                    vis_flx += w * s.albout_lcl[c, b]; vis_wsum += w
                end
                o.albout[c, IVIS] = vis_flx / vis_wsum
            end

            # NIR band average.
            nir_flx = zero(T)
            for b in nir_bnd_bgn:nir_bnd_end
                w = (flg_slr_in == 1 ? opt.flx_wgt_dir[b] : opt.flx_wgt_dif[b])
                nir_flx += w * s.albout_lcl[c, b]
            end
            o.albout[c, INIR] = nir_flx / nir_sum

            # Absorbed flux into VIS band.
            if numrad_snw == DEFAULT_NUMBER_BANDS
                for i in snl_top_j:(nlevsno+1)
                    o.flx_abs[c, i, IVIS] = s.flx_abs_lcl[c, i, 1]
                end
            elseif numrad_snw == HIGH_NUMBER_BANDS
                for i in snl_top_j:(nlevsno+1)
                    fsum = zero(T); wsum = zero(T)
                    for b in 1:(nir_bnd_bgn-1)
                        w = (flg_slr_in == 1 ? opt.flx_wgt_dir[b] : opt.flx_wgt_dif[b])
                        fsum += w * s.flx_abs_lcl[c, i, b]; wsum += w
                    end
                    o.flx_abs[c, i, IVIS] = fsum / wsum
                end
            end

            # NIR flux absorption.
            for i in snl_top_j:(nlevsno+1)
                fsum = zero(T)
                for b in nir_bnd_bgn:nir_bnd_end
                    w = (flg_slr_in == 1 ? opt.flx_wgt_dir[b] : opt.flx_wgt_dif[b])
                    fsum += w * s.flx_abs_lcl[c, i, b]
                end
                o.flx_abs[c, i, INIR] = fsum / nir_sum
            end

            # High solar zenith angle adjustment for NIR direct albedo.
            if mu_not < T(MU_75) && flg_slr_in == 1
                sza_c1_val = T(SZA_A0) + T(SZA_A1) * mu_not + T(SZA_A2) * mu_not^2
                sza_c0_val = T(SZA_B0) + T(SZA_B1) * mu_not + T(SZA_B2) * mu_not^2
                sza_factor = sza_c1_val * (log10(T(s.snw_rds_lcl[c, snl_top_j])) - T(6.0)) + sza_c0_val
                flx_sza_adjust = o.albout[c, INIR] * (sza_factor - one(T)) * nir_sum
                o.albout[c, INIR] *= sza_factor
                o.flx_abs[c, snl_top_j, INIR] -= flx_sza_adjust
            end
        elseif coszc > zero(T) && h2osno_lcl < T(MIN_SNW) && h2osno_lcl > zero(T)
            o.albout[c, IVIS] = in.albsfc[c, IVIS]
            o.albout[c, INIR] = in.albsfc[c, INIR]
        else
            o.albout[c, IVIS] = zero(T)
            o.albout[c, INIR] = zero(T)
        end
    end
end

# --------------------------------------------------------------------------
# Input-prep kernel — build SNICAR snow-layer inputs from device column state
# (mirrors the host-fallback host loop in surface_albedo!). One thread/column.
# --------------------------------------------------------------------------
@kernel function _snicar_prep_kernel!(liq_snw, ice_snw, rds_int, aer3,
        @Const(mask), @Const(h2osoi_liq), @Const(h2osoi_ice), @Const(snw_rds_col),
        @Const(bcphi), @Const(bcpho), @Const(ocphi), @Const(ocpho),
        @Const(dst1), @Const(dst2), @Const(dst3), @Const(dst4),
        nlevsno::Int, aer_avail::Int, rds_min_int::Int)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(liq_snw)
        for j in 1:nlevsno
            liq_snw[c, j] = smooth_max(zero(T), h2osoi_liq[c, j])
            ice_snw[c, j] = smooth_max(zero(T), h2osoi_ice[c, j])
            rds = snw_rds_col[c, j]
            if !isnan(rds) && rds > zero(rds)
                # round-then-fptosi: `unsafe_trunc` lowers to a direct float→int
                # cast (no boxed Int64(::Float32) conversion, which gpu_malloc's).
                rds_int[c, j] = unsafe_trunc(Int, round(clamp(rds, T(SNW_RDS_MIN_TBL), T(SNW_RDS_MAX_TBL))))
            else
                rds_int[c, j] = rds_min_int
            end
        end
        if c <= aer_avail
            for j in 1:nlevsno
                aer3[c, j, 1] = bcphi[c, j]; aer3[c, j, 2] = bcpho[c, j]
                aer3[c, j, 3] = ocphi[c, j]; aer3[c, j, 4] = ocpho[c, j]
                aer3[c, j, 5] = dst1[c, j];  aer3[c, j, 6] = dst2[c, j]
                aer3[c, j, 7] = dst3[c, j];  aer3[c, j, 8] = dst4[c, j]
            end
        end
    end
end

# --------------------------------------------------------------------------
# snicar_rt_device! — host wrapper that builds device bundles and launches.
# --------------------------------------------------------------------------
"""
    snicar_rt_device!(coszen, flg_slr_in, h2osno_liq, h2osno_ice, h2osno_total,
                      snw_rds, mss_cnc_aer_in, albsfc, snl, frac_sno,
                      albout, flx_abs, nlevsno;
                      optics=snicar_optics, params=snicar_params, mask=nothing)

On-device equivalent of [`snicar_rt!`](@ref). One thread per column; runs on the
backend of `albout` (CPU Array → KA CPU loop; MtlArray/CuArray → GPU). Arrays are
assumed to already live on the working backend at the working element type;
`snw_rds`/`snl` are integer arrays. `frac_sno` is accepted for signature parity
(unused, as in the host solver).
"""
function snicar_rt_device!(coszen, flg_slr_in::Int,
                           h2osno_liq, h2osno_ice, h2osno_total,
                           snw_rds, mss_cnc_aer_in, albsfc, snl, frac_sno,
                           albout, flx_abs, nlevsno::Int;
                           optics::SnicarOpticsData=snicar_optics,
                           params::SnicarParams=snicar_params,
                           mask=nothing)
    ncols = length(coszen)
    T = eltype(albout)
    numrad_snw = varctl.snicar_numrad_snw
    snobc = varctl.snicar_snobc_intmix
    snodst = varctl.snicar_snodst_intmix

    shape_str = varctl.snicar_snw_shape
    sno_shp_code = shape_str == "sphere" ? 0 :
                   shape_str == "spheroid" ? 1 :
                   shape_str == "hexagonal_plate" ? 2 :
                   shape_str == "koch_snowflake" ? 3 :
                   error("SNICAR ERROR: unknown sno_shp: $shape_str")

    if numrad_snw == DEFAULT_NUMBER_BANDS
        nir_bnd_bgn = 2
    elseif numrad_snw == HIGH_NUMBER_BANDS
        nir_bnd_bgn = 51
    else
        error("SNICAR ERROR: unknown snicar_numrad_snw value: $numrad_snw")
    end

    # Move a host Array onto the working backend at element type T (no-op-ish on
    # CPU: Array{T} copy). Mirrors the snowage_grain! `_age` helper.
    _dev(a) = albout isa Array ? Array{T}(a) :
              copyto!(similar(albout, T, size(a)), T.(a))

    # Band-center wavelengths and nonspherical/internal-mixing band centers (host).
    wvl_ct = zeros(Float64, numrad_snw)
    if numrad_snw == DEFAULT_NUMBER_BANDS
        wvl_ct .= [0.5, 0.85, 1.1, 1.35, 3.25]
    else
        for igb in 1:numrad_snw
            wvl_ct[igb] = 0.205 + 0.01 * (igb - 1.0)
        end
    end
    g_wvl_ct = [SNICAR_G_WVL[igb+1] * 0.5 + SNICAR_G_WVL[igb] * 0.5 for igb in 1:SNICAR_SEVEN_BANDS]
    dstint_wvl_ct = [DSTINT_WVL[idb+1] * 0.5 + DSTINT_WVL[idb] * 0.5 for idb in 1:SNICAR_SIZE_BINS]
    bcint_wvl_ct = [BCINT_WVL[ibb+1] * 0.5 + BCINT_WVL[ibb] * 0.5 for ibb in 1:SNICAR_SIXTEEN_BANDS]

    cst = _SnicarConst(
        wvl_ct = _dev(wvl_ct), g_wvl_ct = _dev(g_wvl_ct),
        dstint_wvl_ct = _dev(dstint_wvl_ct), bcint_wvl_ct = _dev(bcint_wvl_ct),
        difgauspt = _dev(SNICAR_DIFGAUSPT), difgauswt = _dev(SNICAR_DIFGAUSWT),
        g_b0 = _dev(SNICAR_G_B0), g_b1 = _dev(SNICAR_G_B1), g_b2 = _dev(SNICAR_G_B2),
        g_f07_c0 = _dev(SNICAR_G_F07_C0), g_f07_c1 = _dev(SNICAR_G_F07_C1), g_f07_c2 = _dev(SNICAR_G_F07_C2),
        g_f07_p0 = _dev(SNICAR_G_F07_P0), g_f07_p1 = _dev(SNICAR_G_F07_P1), g_f07_p2 = _dev(SNICAR_G_F07_P2),
        bcint_d0 = _dev(BCINT_D0), bcint_d1 = _dev(BCINT_D1), bcint_d2 = _dev(BCINT_D2),
        bcint_m = _dev(BCINT_M), bcint_n = _dev(BCINT_N),
        dstint_a1 = _dev(DSTINT_A1), dstint_a2 = _dev(DSTINT_A2), dstint_a3 = _dev(DSTINT_A3),
    )

    opt = _SnicarOptics(
        ss_alb_snw_drc = _dev(optics.ss_alb_snw_drc), asm_prm_snw_drc = _dev(optics.asm_prm_snw_drc),
        ext_cff_mss_snw_drc = _dev(optics.ext_cff_mss_snw_drc),
        ss_alb_snw_dfs = _dev(optics.ss_alb_snw_dfs), asm_prm_snw_dfs = _dev(optics.asm_prm_snw_dfs),
        ext_cff_mss_snw_dfs = _dev(optics.ext_cff_mss_snw_dfs),
        ss_alb_bc_hphil = _dev(optics.ss_alb_bc_hphil), asm_prm_bc_hphil = _dev(optics.asm_prm_bc_hphil),
        ext_cff_mss_bc_hphil = _dev(optics.ext_cff_mss_bc_hphil),
        ss_alb_bc_hphob = _dev(optics.ss_alb_bc_hphob), asm_prm_bc_hphob = _dev(optics.asm_prm_bc_hphob),
        ext_cff_mss_bc_hphob = _dev(optics.ext_cff_mss_bc_hphob),
        ss_alb_oc_hphil = _dev(optics.ss_alb_oc_hphil), asm_prm_oc_hphil = _dev(optics.asm_prm_oc_hphil),
        ext_cff_mss_oc_hphil = _dev(optics.ext_cff_mss_oc_hphil),
        ss_alb_oc_hphob = _dev(optics.ss_alb_oc_hphob), asm_prm_oc_hphob = _dev(optics.asm_prm_oc_hphob),
        ext_cff_mss_oc_hphob = _dev(optics.ext_cff_mss_oc_hphob),
        ss_alb_dst1 = _dev(optics.ss_alb_dst1), asm_prm_dst1 = _dev(optics.asm_prm_dst1),
        ext_cff_mss_dst1 = _dev(optics.ext_cff_mss_dst1),
        ss_alb_dst2 = _dev(optics.ss_alb_dst2), asm_prm_dst2 = _dev(optics.asm_prm_dst2),
        ext_cff_mss_dst2 = _dev(optics.ext_cff_mss_dst2),
        ss_alb_dst3 = _dev(optics.ss_alb_dst3), asm_prm_dst3 = _dev(optics.asm_prm_dst3),
        ext_cff_mss_dst3 = _dev(optics.ext_cff_mss_dst3),
        ss_alb_dst4 = _dev(optics.ss_alb_dst4), asm_prm_dst4 = _dev(optics.asm_prm_dst4),
        ext_cff_mss_dst4 = _dev(optics.ext_cff_mss_dst4),
        flx_wgt_dir = _dev(optics.flx_wgt_dir), flx_wgt_dif = _dev(optics.flx_wgt_dif),
    )

    # Per-column scratch on the working backend.
    _m(k) = fill!(similar(albout, T, ncols, k), zero(T))
    _a3(k1, k2) = fill!(similar(albout, T, ncols, k1, k2), zero(T))
    _mi(k) = fill!(similar(snw_rds, eltype(snw_rds), ncols, k), zero(eltype(snw_rds)))

    scr = _SnicarScratch(
        h2osno_ice_lcl = _m(nlevsno), h2osno_liq_lcl = _m(nlevsno), snw_rds_lcl = _mi(nlevsno),
        mss_cnc_aer_lcl = _a3(nlevsno, SNO_NBR_AER),
        ss_alb_snw_lcl = _m(nlevsno), asm_prm_snw_lcl = _m(nlevsno), ext_cff_mss_snw_lcl = _m(nlevsno),
        ss_alb_aer_lcl = _m(SNO_NBR_AER), asm_prm_aer_lcl = _m(SNO_NBR_AER), ext_cff_mss_aer_lcl = _m(SNO_NBR_AER),
        g_ice_Cg_tmp = _m(SNICAR_SEVEN_BANDS), gg_ice_F07_tmp = _m(SNICAR_SEVEN_BANDS),
        L_snw = _m(nlevsno), tau_snw = _m(nlevsno),
        L_aer = _a3(nlevsno, SNO_NBR_AER), tau_aer = _a3(nlevsno, SNO_NBR_AER),
        tau_arr = _m(nlevsno), omega_arr = _m(nlevsno), g_arr = _m(nlevsno),
        enh_omg_bcint_tmp = _m(SNICAR_SIXTEEN_BANDS), enh_omg_bcint_tmp2 = _m(SNICAR_SIXTEEN_BANDS),
        enh_omg_dstint_tmp = _m(SNICAR_SIZE_BINS), enh_omg_dstint_tmp2 = _m(SNICAR_SIZE_BINS),
        tau_star = _m(nlevsno), omega_star = _m(nlevsno), g_star = _m(nlevsno),
        trndir = _m(nlevsno+1), trntdr = _m(nlevsno+1), trndif = _m(nlevsno+1),
        rupdir = _m(nlevsno+1), rupdif = _m(nlevsno+1), rdndif = _m(nlevsno+1),
        dfdir = _m(nlevsno+1), dfdif = _m(nlevsno+1), dftmp = _m(nlevsno+1),
        rdir = _m(nlevsno), rdif_a = _m(nlevsno), rdif_b = _m(nlevsno),
        tdir_arr = _m(nlevsno), tdif_a = _m(nlevsno), tdif_b = _m(nlevsno),
        trnlay = _m(nlevsno), F_abs = _m(nlevsno),
        albout_lcl = _m(numrad_snw), flx_abs_lcl = _a3(nlevsno+1, numrad_snw),
    )

    out = _SnicarOut(albout = albout, flx_abs = flx_abs)
    inb = _SnicarIn(coszen = coszen, h2osno_liq = h2osno_liq, h2osno_ice = h2osno_ice,
                    h2osno_total = h2osno_total, snw_rds = snw_rds,
                    mss_cnc_aer_in = mss_cnc_aer_in, albsfc = albsfc, snl = snl)

    # Default mask (all columns) on the albout backend so kernels stay scalar-free.
    msk = mask === nothing ? fill!(similar(albout, Bool, ncols), true) : mask

    snw_rds_min_int = round(Int, params.snw_rds_min)
    be = _kernel_backend(albout)

    ncols == 0 && return nothing

    _snicar_setup_kernel!(be)(out, scr, inb, msk, nlevsno, snw_rds_min_int; ndrange = ncols)
    KA.synchronize(be)
    for b in 1:numrad_snw
        _snicar_band_kernel!(be)(out, scr, opt, cst, inb, msk,
            flg_slr_in, b, numrad_snw, nir_bnd_bgn, nlevsno, sno_shp_code, snobc, snodst;
            ndrange = ncols)
        KA.synchronize(be)
    end
    _snicar_finalize_kernel!(be)(out, scr, opt, inb, msk,
        flg_slr_in, numrad_snw, nir_bnd_bgn, nlevsno; ndrange = ncols)
    KA.synchronize(be)

    return nothing
end
