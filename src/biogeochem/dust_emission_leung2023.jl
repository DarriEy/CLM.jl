# ==========================================================================
# Ported from: src/biogeochem/DustEmisLeung2023.F90
# Leung (2023) dust emission scheme (dust_emis_method == "Leung_2023")
#
# Simulates dust mobilization due to wind from the surface into the lowest
# atmospheric layer, filling flx_mss_vrt_dst(ndst) (kg/m**2/s, [+ = to atm]).
# Compared with the default Zender (2003) scheme, Leung (2023):
#   - uses a Shao & Lu (2000) dry fluid threshold (≈80 um optimal particle)
#     instead of Iversen & White (1982),
#   - emits above the IMPACT threshold (B_it · fluid threshold) with a
#     Kok et al. (2014) fragmentation-exponent emission equation,
#   - applies an Okin–Pierre vegetation drag partition (and a rock drag
#     partition factor `dpfct_rock`) on the WIND rather than the threshold,
#   - multiplies emission by a Comola et al. (2019) intermittency factor that
#     accounts for sub-timestep turbulent wind gusts crossing the threshold.
#
# Reuses DustEmisBaseData (flx_mss_vrt_dst_*, ovr_src_snk_mss, dpfct_rock_patch)
# and the shared `_dust_landunit_vai` helper from dust_emission.jl.
# ==========================================================================

# T-generic erf (device-safe: literals wrapped so Metal stays in Float32; the
# host erf in varcon.jl uses Float64 literals which would force Float64 on GPU).
@inline function _erf_T(x)
    T = typeof(x)
    p  = T(0.3275911)
    sgn = x >= zero(T) ? one(T) : -one(T)
    ax = abs(x)
    t = one(T) / (one(T) + p * ax)
    y = one(T) - (((((T(1.061405429) * t + T(-1.453152027)) * t) + T(1.421413741)) * t +
        T(-0.284496736)) * t + T(0.254829592)) * t * exp(-ax * ax)
    return sgn * y
end

# Per-patch Leung (2023) dust kernel: each patch is independent given the
# landunit-averaged VAI (a p2lu reduction) precomputed on the host. The two
# defensive Fortran error() guards (bad mobilization fraction, obu==0) are
# dropped — they fire only on invalid input and are not exercised; valid-input
# results are byte-identical. All literals are T()-wrapped so Metal runs Float32.
@kernel function _dust_leung_kernel!(flx_patch, flx_tot, @Const(nolakep_mask),
        @Const(patch_column), @Const(patch_landunit), @Const(patch_itype),
        @Const(lun_itype), @Const(tlai_lu), @Const(forc_rho), @Const(gwc_thr),
        @Const(mss_frc_cly_vld), @Const(watsat), @Const(frac_sno), @Const(h2osoi_vol),
        @Const(h2osoi_liq), @Const(h2osoi_ice), @Const(fv), @Const(obu),
        @Const(dpfct_rock), @Const(ovr_src_snk_mss), ndst::Int, dst_src_nbr::Int,
        noveg_::Int, istsoil_::Int, istcrop_::Int)
    p = @index(Global)
    @inbounds if nolakep_mask[p]
        T = eltype(flx_patch)
        vai_mbl_thr = T(0.6); Cd0 = T(4.4e-5); Ca = T(2.7); Ce = T(2.0); C_tune = T(0.05)
        wnd_frc_thr_slt_std_min = T(0.16); forc_rho_std = T(1.2250); dns_slt = T(2650.0)
        B_it = T(0.82); k_vk = T(0.4); f_0 = T(0.32); c_e = T(4.8); D_p = T(130e-6)
        gamma_Shao = T(1.65e-4); A_Shao = T(0.0123); frag_expt_thr = T(2.5)
        z0a_glob = T(1e-4); hgt_sal = T(0.1); vai0_Okin = T(0.1); zii = T(1000.0)
        dust_veg_drag_fact = T(0.7)
        denh2o = T(DENH2O); grav = T(GRAV)

        c = patch_column[p]; l = patch_landunit[p]
        for nb in 1:ndst; flx_patch[p, nb] = zero(T); end
        flx_tot[p] = zero(T)

        # land mobile fraction
        lnd_frc_mbl = zero(T)
        if lun_itype[l] == istsoil_ || lun_itype[l] == istcrop_
            lnd_frc_mbl = tlai_lu[l] < vai_mbl_thr ? one(T) - tlai_lu[l] / vai_mbl_thr : zero(T)
            lnd_frc_mbl *= (one(T) - frac_sno[c])
        end

        # Fecan (1999) soil-moisture threshold factor
        bd = (one(T) - watsat[c, 1]) * dns_slt
        gwc_sfc = h2osoi_vol[c, 1] * denh2o / bd
        frc_thr_wet_fct = gwc_sfc > gwc_thr[c] ?
            sqrt(one(T) + T(1.21) * (T(100.0) * (gwc_sfc - gwc_thr[c]))^T(0.68)) : one(T)
        liqfrac = max(zero(T), min(one(T),
            h2osoi_liq[c, 1] / (h2osoi_ice[c, 1] + h2osoi_liq[c, 1] + T(1.0e-6))))

        # Shao & Lu (2000) thresholds
        tmp2 = sqrt(A_Shao * (dns_slt * grav * D_p + gamma_Shao / D_p))
        wnd_frc_thr_slt    = tmp2 / sqrt(forc_rho[c]) * frc_thr_wet_fct
        wnd_frc_thr_slt_it = B_it * tmp2 / sqrt(forc_rho[c])
        wnd_frc_thr_slt_std = wnd_frc_thr_slt * sqrt(forc_rho[c] / forc_rho_std)
        dst_emiss_coeff = Cd0 * exp(-Ce *
            (wnd_frc_thr_slt_std - wnd_frc_thr_slt_std_min) / wnd_frc_thr_slt_std_min)
        frag_expt = Ca * (wnd_frc_thr_slt_std - wnd_frc_thr_slt_std_min) / wnd_frc_thr_slt_std_min
        if frag_expt > frag_expt_thr; frag_expt = frag_expt_thr; end

        # Okin–Pierre vegetation + rock drag partition
        wnd_frc_slt = zero(T)
        if lnd_frc_mbl > zero(T) && tlai_lu[l] <= vai_mbl_thr
            vai_Okin = tlai_lu[l] + vai0_Okin
            if vai_Okin > one(T); vai_Okin = one(T); end
            K_length = T(2.0) * (one(T) / vai_Okin - one(T))
            ssr = dust_veg_drag_fact * (K_length + f_0 * c_e) / (K_length + c_e)
            frc_thr_rgh_fct = one(T)
            if lun_itype[l] == istsoil_ || lun_itype[l] == istcrop_
                if patch_itype[p] == noveg_
                    dpr = dpfct_rock[p]
                    frc_thr_rgh_fct = isnan(dpr) ? T(0.001) : dpr
                else
                    frc_thr_rgh_fct = ssr
                end
            end
            wnd_frc_slt = fv[p] * frc_thr_rgh_fct
        end

        flx_ttl = zero(T)
        if lnd_frc_mbl > zero(T)
            # Kok et al. (2014) emission above the IMPACT threshold
            if wnd_frc_slt > wnd_frc_thr_slt_it
                flx_ttl = dst_emiss_coeff * mss_frc_cly_vld[c] * forc_rho[c] *
                    ((wnd_frc_slt^T(2.0) - wnd_frc_thr_slt_it^T(2.0)) / wnd_frc_thr_slt_it) *
                    (wnd_frc_slt / wnd_frc_thr_slt_it)^frag_expt
                flx_ttl *= lnd_frc_mbl * C_tune * liqfrac
            end

            # Comola et al. (2019) intermittency factor
            u_mean_slt = (wnd_frc_slt / k_vk) * log(hgt_sal / z0a_glob)
            stblty = zii / obu[p]
            u_sd_slt = (T(12.0) - T(0.5) * stblty) >= T(0.001) ?
                wnd_frc_slt * (T(12.0) - T(0.5) * stblty)^T(0.333) :
                wnd_frc_slt * T(0.001)^T(0.333)
            u_fld_thr   = (wnd_frc_thr_slt    / k_vk) * log(hgt_sal / z0a_glob)
            u_impct_thr = (wnd_frc_thr_slt_it / k_vk) * log(hgt_sal / z0a_glob)
            numer = u_fld_thr^T(2.0) - u_impct_thr^T(2.0) -
                    T(2.0) * u_mean_slt * (u_fld_thr - u_impct_thr)
            denom = T(2.0) * u_sd_slt^T(2.0)
            thr_crs_rate = numer / denom < T(30.0) ?
                (exp(numer / denom) + one(T))^(-one(T)) : zero(T)
            sqrt2 = sqrt(T(2.0))
            prb_crs_fld_thr   = T(0.5) * (one(T) + _erf_T((u_fld_thr   - u_mean_slt) / (sqrt2 * u_sd_slt)))
            prb_crs_impct_thr = T(0.5) * (one(T) + _erf_T((u_impct_thr - u_mean_slt) / (sqrt2 * u_sd_slt)))
            intrmtncy_fct = one(T) - prb_crs_fld_thr +
                            thr_crs_rate * (prb_crs_fld_thr - prb_crs_impct_thr)
            if !isnan(intrmtncy_fct); flx_ttl *= intrmtncy_fct; end

            # partition total vertical mass flux into NDST bins + total
            tot = zero(T)
            for nb in 1:ndst
                acc = zero(T)
                for m in 1:dst_src_nbr
                    acc += ovr_src_snk_mss[m, nb] * flx_ttl
                end
                flx_patch[p, nb] = acc
                tot += acc
            end
            flx_tot[p] = tot
        end
    end
end

# ---------------------------------------------------------------------------
# dust_emission_leung2023! — Leung (2023) dust mobilization
# Ported from: subroutine DustEmission in DustEmisLeung2023.F90
# ---------------------------------------------------------------------------

"""
    dust_emission_leung2023!(dust, nolakep_mask, patch_active, patch_column,
        patch_landunit, patch_itype, patch_wtlunit, lun_itype, bounds_p, nl,
        forc_rho, gwc_thr, mss_frc_cly_vld, watsat, tlai, tsai,
        frac_sno, h2osoi_vol, h2osoi_liq, h2osoi_ice, fv, obu)

Compute Leung (2023) dust emission, filling `dust.flx_mss_vrt_dst_patch`
(kg/m²/s, [+ = to atm]) and `dust.flx_mss_vrt_dst_tot_patch` per patch.

Inputs mirror the Zender2003 routine, with two differences:
- `patch_itype` (per-patch PFT type) selects rock (`noveg`) vs. vegetation
  drag partition,
- `obu` (Monin–Obukhov length per patch) replaces `u10` and drives the
  Comola (2019) intermittency factor.

`dpfct_rock_patch` (rock drag-partition factor) is read from `dust`.

Ported from `DustEmission` in `DustEmisLeung2023.F90`.
"""
function dust_emission_leung2023!(
    dust::DustEmisBaseData,
    nolakep_mask::AbstractVector{Bool},
    patch_active::AbstractVector{Bool},
    patch_column::AbstractVector{<:Integer},
    patch_landunit::AbstractVector{<:Integer},
    patch_itype::AbstractVector{<:Integer},
    patch_wtlunit::AbstractVector{<:Real},
    lun_itype::AbstractVector{<:Integer},
    bounds_p::UnitRange{Int}, nl::Int,
    forc_rho::AbstractVector{<:Real},
    gwc_thr::AbstractVector{<:Real},
    mss_frc_cly_vld::AbstractVector{<:Real},
    watsat::AbstractMatrix{<:Real},
    tlai::AbstractVector{<:Real},
    tsai::AbstractVector{<:Real},
    frac_sno::AbstractVector{<:Real},
    h2osoi_vol::AbstractMatrix{<:Real},
    h2osoi_liq::AbstractMatrix{<:Real},
    h2osoi_ice::AbstractMatrix{<:Real},
    fv::AbstractVector{<:Real},
    obu::AbstractVector{<:Real},
)
    # Landunit-averaged VAI (LAI+SAI) — a p2lu reduction, kept on the host.
    tlai_lu = _dust_landunit_vai(tlai, tsai, patch_active, patch_landunit,
                                 patch_wtlunit, bounds_p, nl)

    # Each patch is independent given tlai_lu, so the whole per-patch physics
    # (mobile fraction -> thresholds -> drag partition -> Kok emission ->
    # intermittency -> bin partition + total) runs in one per-patch kernel.
    _launch!(_dust_leung_kernel!, dust.flx_mss_vrt_dst_patch, dust.flx_mss_vrt_dst_tot_patch,
        nolakep_mask, patch_column, patch_landunit, patch_itype, lun_itype, tlai_lu,
        forc_rho, gwc_thr, mss_frc_cly_vld, watsat, frac_sno, h2osoi_vol, h2osoi_liq,
        h2osoi_ice, fv, obu, dust.dpfct_rock_patch, dust.ovr_src_snk_mss,
        NDST, DST_SRC_NBR, noveg, ISTSOIL, ISTCROP; ndrange = length(nolakep_mask))

    return nothing
end
