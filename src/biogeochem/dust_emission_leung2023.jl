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
    ndst = NDST
    dst_src_nbr = DST_SRC_NBR

    # --- Constants (from DustEmisLeung2023.F90) ---
    vai_mbl_thr = 0.6                 # [m2 m-2] VAI threshold quenching mobilization
    Cd0     = 4.4e-5                  # [dimless] dust emission coefficient constant
    Ca      = 2.7                     # [dimless] fragmentation-exponent scaling
    Ce      = 2.0                     # [dimless] emission-coefficient exponent scaling
    C_tune  = 0.05                    # [dimless] global vertical-flux tuning constant
    wnd_frc_thr_slt_std_min = 0.16    # [m/s] min standardized soil threshold
    forc_rho_std = 1.2250             # [kg/m3] standard air density
    dns_slt = 2650.0                  # [kg m-3] saltation particle density
    B_it    = 0.82                    # [dimless] impact/fluid threshold ratio
    k_vk    = 0.4                     # [dimless] von Karman constant
    f_0     = 0.32                    # [dimless] SSR in lee of a plant
    c_e     = 4.8                     # [dimless] e-folding distance velocity recovery
    D_p     = 130e-6                  # [m] medium soil particle diameter
    gamma_Shao = 1.65e-4              # [kg s-2] interparticle cohesion (S&L00)
    A_Shao  = 0.0123                  # [dimless] aerodynamic-force coefficient (S&L00)
    frag_expt_thr = 2.5               # [dimless] max fragmentation exponent
    z0a_glob = 1e-4                   # [m] assumed global aeolian roughness length
    hgt_sal  = 0.1                    # [m] saltation height (log-law of the wall)
    vai0_Okin = 0.1                   # [m2/m2] min VAI for Okin–Pierre equation
    zii = 1000.0                      # [m] convective boundary-layer height
    dust_veg_drag_fact = 0.7          # [dimless] vegetation drag-partition tuning

    # Landunit-averaged VAI (LAI+SAI).
    tlai_lu = _dust_landunit_vai(tlai, tsai, patch_active, patch_landunit,
                                 patch_wtlunit, bounds_p, nl)

    np = length(patch_active)
    lnd_frc_mbl      = zeros(Float64, np)
    frc_thr_rghn_fct = zeros(Float64, np)

    # Land mobile fraction + output reset.
    for p in bounds_p
        nolakep_mask[p] || continue
        c = patch_column[p]
        l = patch_landunit[p]
        for nb in 1:ndst
            dust.flx_mss_vrt_dst_patch[p, nb] = 0.0
        end
        dust.flx_mss_vrt_dst_tot_patch[p] = 0.0

        if lun_itype[l] == ISTSOIL || lun_itype[l] == ISTCROP
            if tlai_lu[l] < vai_mbl_thr
                lnd_frc_mbl[p] = 1.0 - tlai_lu[l] / vai_mbl_thr
            else
                lnd_frc_mbl[p] = 0.0
            end
            lnd_frc_mbl[p] *= (1.0 - frac_sno[c])
        else
            lnd_frc_mbl[p] = 0.0
        end
        if lnd_frc_mbl[p] > 1.0 || lnd_frc_mbl[p] < 0.0
            error("Bad value for dust mobilization fraction at patch $p: $(lnd_frc_mbl[p])")
        end
    end

    flx_mss_vrt_dst_ttl = zeros(Float64, np)
    for p in bounds_p
        nolakep_mask[p] || continue
        c = patch_column[p]
        l = patch_landunit[p]

        # --- Fecan (1999) soil-moisture threshold factor ---
        bd = (1.0 - watsat[c, 1]) * dns_slt          # [kg m-3] dry soil bulk density
        gwc_sfc = h2osoi_vol[c, 1] * DENH2O / bd      # [kg kg-1] gravimetric H2O
        if gwc_sfc > gwc_thr[c]
            frc_thr_wet_fct = sqrt(1.0 + 1.21 * (100.0 * (gwc_sfc - gwc_thr[c]))^0.68)
        else
            frc_thr_wet_fct = 1.0
        end

        liqfrac = max(0.0, min(1.0,
            h2osoi_liq[c, 1] / (h2osoi_ice[c, 1] + h2osoi_liq[c, 1] + 1.0e-6)))

        # --- Shao & Lu (2000) dry fluid threshold, fluid & impact thresholds ---
        tmp2 = sqrt(A_Shao * (dns_slt * GRAV * D_p + gamma_Shao / D_p))
        wnd_frc_thr_slt    = tmp2 / sqrt(forc_rho[c]) * frc_thr_wet_fct  # fluid threshold
        wnd_frc_thr_slt_it = B_it * tmp2 / sqrt(forc_rho[c])            # impact threshold

        # Standardized fluid threshold + emission coefficient + frag. exponent.
        wnd_frc_thr_slt_std = wnd_frc_thr_slt * sqrt(forc_rho[c] / forc_rho_std)
        dst_emiss_coeff = Cd0 * exp(-Ce *
            (wnd_frc_thr_slt_std - wnd_frc_thr_slt_std_min) / wnd_frc_thr_slt_std_min)
        frag_expt = Ca * (wnd_frc_thr_slt_std - wnd_frc_thr_slt_std_min) /
                    wnd_frc_thr_slt_std_min
        if frag_expt > frag_expt_thr
            frag_expt = frag_expt_thr
        end

        # --- Drag partition: Okin–Pierre vegetation + rock drag partition ---
        if lnd_frc_mbl[p] > 0.0 && tlai_lu[l] <= vai_mbl_thr
            vai_Okin = tlai_lu[l] + vai0_Okin
            if vai_Okin > 1.0
                vai_Okin = 1.0
            end
            K_length = 2.0 * (1.0 / vai_Okin - 1.0)
            ssr = dust_veg_drag_fact * (K_length + f_0 * c_e) / (K_length + c_e)

            if lun_itype[l] == ISTSOIL || lun_itype[l] == ISTCROP
                if patch_itype[p] == noveg          # bare → rock drag partition
                    dpr = dust.dpfct_rock_patch[p]
                    frc_thr_rgh_fct = isnan(dpr) ? 0.001 : dpr
                else                                # vegetated → SSR drag partition
                    frc_thr_rgh_fct = ssr
                end
            else
                frc_thr_rgh_fct = 1.0
            end
            wnd_frc_slt = fv[p] * frc_thr_rgh_fct
            frc_thr_rghn_fct[p] = frc_thr_rgh_fct
        else
            wnd_frc_slt = fv[p] * frc_thr_rghn_fct[p]
            frc_thr_rghn_fct[p] = 0.0
        end

        lnd_frc_mbl[p] > 0.0 || continue

        flx_mss_hrz_slt_ttl = 0.0
        flx_mss_vrt_dst_ttl[p] = 0.0

        # --- Kok et al. (2014) emission above the IMPACT threshold ---
        if wnd_frc_slt > wnd_frc_thr_slt_it
            flx_mss_vrt_dst_ttl[p] = dst_emiss_coeff * mss_frc_cly_vld[c] * forc_rho[c] *
                ((wnd_frc_slt^2.0 - wnd_frc_thr_slt_it^2.0) / wnd_frc_thr_slt_it) *
                (wnd_frc_slt / wnd_frc_thr_slt_it)^frag_expt
            flx_mss_vrt_dst_ttl[p] *= lnd_frc_mbl[p] * C_tune * liqfrac
        end

        # --- Comola et al. (2019) intermittency factor ---
        u_mean_slt = (wnd_frc_slt / k_vk) * log(hgt_sal / z0a_glob)
        if obu[p] == 0.0
            error("DustEmisLeung2023: input obu is zero at patch $p")
        end
        stblty = zii / obu[p]
        if (12.0 - 0.5 * stblty) >= 0.001
            u_sd_slt = wnd_frc_slt * (12.0 - 0.5 * stblty)^0.333
        else
            u_sd_slt = wnd_frc_slt * (0.001)^0.333
        end

        u_fld_thr   = (wnd_frc_thr_slt    / k_vk) * log(hgt_sal / z0a_glob)
        u_impct_thr = (wnd_frc_thr_slt_it / k_vk) * log(hgt_sal / z0a_glob)

        numer = u_fld_thr^2.0 - u_impct_thr^2.0 -
                2.0 * u_mean_slt * (u_fld_thr - u_impct_thr)
        denom = 2.0 * u_sd_slt^2.0
        if numer / denom < 30.0
            thr_crs_rate = (exp(numer / denom) + 1.0)^(-1.0)
        else
            thr_crs_rate = 0.0
        end

        prb_crs_fld_thr   = 0.5 * (1.0 + _dust_erf((u_fld_thr   - u_mean_slt) /
                                                   (sqrt(2.0) * u_sd_slt)))
        prb_crs_impct_thr = 0.5 * (1.0 + _dust_erf((u_impct_thr - u_mean_slt) /
                                                   (sqrt(2.0) * u_sd_slt)))
        intrmtncy_fct = 1.0 - prb_crs_fld_thr +
                        thr_crs_rate * (prb_crs_fld_thr - prb_crs_impct_thr)

        # Multiply emission by intermittency (skip if NaN, mirroring Fortran).
        if !isnan(intrmtncy_fct)
            flx_mss_vrt_dst_ttl[p] *= intrmtncy_fct
        end
    end

    # Partition total vertical mass flux into NDST transport bins.
    for nb in 1:ndst
        for m in 1:dst_src_nbr
            for p in bounds_p
                nolakep_mask[p] || continue
                lnd_frc_mbl[p] > 0.0 || continue
                dust.flx_mss_vrt_dst_patch[p, nb] +=
                    dust.ovr_src_snk_mss[m, nb] * flx_mss_vrt_dst_ttl[p]
            end
        end
    end
    for nb in 1:ndst
        for p in bounds_p
            nolakep_mask[p] || continue
            lnd_frc_mbl[p] > 0.0 || continue
            dust.flx_mss_vrt_dst_tot_patch[p] += dust.flx_mss_vrt_dst_patch[p, nb]
        end
    end

    return nothing
end
