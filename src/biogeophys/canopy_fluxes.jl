# ==========================================================================
# Ported from: src/biogeophys/CanopyFluxesMod.F90 (1756 lines)
# Calculates leaf temperature and surface fluxes using Newton-Raphson
# iteration to solve the surface energy budget:
#   f(t_veg) = Net radiation - Sensible - Latent = 0
# ==========================================================================

# --- Local constants (from CanopyFluxes subroutine) ---
const BTRAN0             = 0.0     # initial btran value
const ZII_CANOPY         = 1000.0  # convective boundary layer height [m]
const BETA_CANOPY        = 1.0     # coefficient of convective velocity [-]
const DELMAX_CANOPY      = 1.0     # max change in leaf temperature [K]
const DLEMIN_CANOPY      = 0.1     # max limit for energy flux convergence [W/m2]
const DTMIN_CANOPY       = 0.01    # max limit for temperature convergence [K]
const ITMIN_CANOPY       = 2       # minimum number of iterations [-]
const RIA_CANOPY         = 0.5     # free parameter for stable formulation
const ABOVE_CANOPY       = 1       # index for above-canopy resistances
const BELOW_CANOPY       = 2       # index for below-canopy resistances

# Biomass heat storage tuning parameters
const K_VERT_CANOPY      = 0.1     # vertical distribution of stem
const K_CYL_VOL_CANOPY   = 1.0     # departure from cylindrical volume
const K_CYL_AREA_CANOPY  = 1.0     # departure from cylindrical area
const K_INTERNAL_CANOPY  = 0.0     # self-absorption of leaf/stem longwave
const MIN_STEM_DIAMETER  = 0.05    # minimum stem diameter [m]

# --- Parameters type (from params_type in Fortran) ---

"""
    CanopyFluxesParamsData

Parameters for canopy fluxes, read from parameter file in Fortran.
Ported from `params_type` in `CanopyFluxesMod.F90`.
"""
Base.@kwdef mutable struct CanopyFluxesParamsData
    lai_dl   ::Float64 = 0.5    # plant litter area index [m2/m2]
    z_dl     ::Float64 = 0.05   # litter layer thickness [m]
    a_coef   ::Float64 = 0.5    # drag coefficient under less dense canopy [-]
    a_exp    ::Float64 = 1.0    # drag exponent under less dense canopy [-]
    csoilc   ::Float64 = 0.004  # soil drag coefficient under dense canopy [-]
    cv       ::Float64 = 0.01   # turbulent transfer coeff canopy surface-to-air [m/s^(1/2)]
    wind_min ::Float64 = 1.0    # minimum wind speed at forcing height [m/s]
end

const canopy_fluxes_params = CanopyFluxesParamsData()

# --- Module control state ---

"""
    CanopyFluxesControl

Module-level control flags for canopy fluxes.
Ported from module variables in `CanopyFluxesMod.F90`.
"""
Base.@kwdef mutable struct CanopyFluxesControl
    perchroot                ::Bool = false   # btran based only on unfrozen soil levels
    perchroot_alt            ::Bool = false   # btran based on active layer
    use_undercanopy_stability::Bool = false   # use undercanopy stability term
    itmax_canopy_fluxes      ::Int  = 40      # max iterations in CanopyFluxes
    use_biomass_heat_storage ::Bool = false   # use biomass heat storage
end

const canopy_fluxes_ctrl = CanopyFluxesControl()

# --- Init / clean / read functions ---

"""
    canopy_fluxes_read_nml!(; use_undercanopy_stability, use_biomass_heat_storage, itmax_canopy_fluxes)

Set namelist parameters for canopy fluxes.
Ported from `CanopyFluxesReadNML` in `CanopyFluxesMod.F90`.
"""
function canopy_fluxes_read_nml!(; use_undercanopy_stability::Bool = false,
                                   use_biomass_heat_storage::Bool = false,
                                   itmax_canopy_fluxes::Int = 40)
    if itmax_canopy_fluxes < 1
        error("canopy_fluxes_read_nml!: itmax_canopy_fluxes must be > 0")
    end
    canopy_fluxes_ctrl.use_undercanopy_stability = use_undercanopy_stability
    canopy_fluxes_ctrl.use_biomass_heat_storage = use_biomass_heat_storage
    canopy_fluxes_ctrl.itmax_canopy_fluxes = itmax_canopy_fluxes
    return nothing
end

"""
    canopy_fluxes_read_params!(; lai_dl, z_dl, a_coef, a_exp, csoilc, cv, wind_min)

Set parameter values for canopy fluxes.
Ported from `readParams` in `CanopyFluxesMod.F90`.
"""
function canopy_fluxes_read_params!(; lai_dl::Real = 0.5,
                                      z_dl::Real = 0.05,
                                      a_coef::Real = 0.5,
                                      a_exp::Real = 1.0,
                                      csoilc::Real = 0.004,
                                      cv::Real = 0.01,
                                      wind_min::Real = 1.0)
    canopy_fluxes_params.lai_dl = lai_dl
    canopy_fluxes_params.z_dl = z_dl
    canopy_fluxes_params.a_coef = a_coef
    canopy_fluxes_params.a_exp = a_exp
    canopy_fluxes_params.csoilc = csoilc
    canopy_fluxes_params.cv = cv
    canopy_fluxes_params.wind_min = wind_min
    return nothing
end

# ---------------------------------------------------------------------------
# Aerodynamic-resistance inner kernel (per-patch resistances / conductances).
#
# Kernelized form of the FIRST per-patch `for fi in 1:fn` loop inside the
# Newton stability iteration of canopy_fluxes!. One thread per filtered patch;
# each patch writes only its own [p] indices (no reductions). Gather form:
# fi = @index(Global); p = filterp[fi].
#
# GPU-hostile host-side items resolved before launch:
#   - csoilc_val (calibration override vs param) is constant across patches and
#     passed as a scalar.
#   - The String/Bool feature flags (use_undercanopy_stability,
#     use_biomass_heat_storage, use_lch4) are passed as scalar Bool args.
# The 2D rah/raw[p, ABOVE_CANOPY|BELOW_CANOPY] indexing is preserved exactly.
# ---------------------------------------------------------------------------
@kernel function _cf_resist_kernel!(
        # written arrays
        tlbef, del2, ram1_patch, rah, raw, uaf_patch, uuc, dleaf_patch,
        rb, rb1_patch, grnd_ch4_cond_patch, svpts, eah, rh_af_patch,
        rah1_patch, raw1_patch, rah2_patch, raw2_patch, vpd_patch,
        # read-only arrays
        @Const(filterp), @Const(column), @Const(gridcell), @Const(itype),
        @Const(t_veg_patch), @Const(del_arr), @Const(ustar_patch),
        @Const(um_patch), @Const(temp1), @Const(temp2),
        @Const(elai_patch), @Const(esai_patch), @Const(htop_patch),
        @Const(z0mg_col), @Const(taf_patch), @Const(qaf_patch),
        @Const(t_grnd_col), @Const(dleaf_pft), @Const(forc_pbot_col),
        @Const(el),
        # scalars
        csoilc_val, use_undercanopy_stability::Bool,
        use_biomass_heat_storage::Bool, use_lch4::Bool,
        cv, a_coef, a_exp, vkc, grav, nu_param, ria_canopy,
        above_canopy::Int, below_canopy::Int)

    fi = @index(Global)
    @inbounds begin
        p = filterp[fi]
        c = column[p]
        g = gridcell[p]
        ft = itype[p] + 1  # 0-based Fortran PFT → 1-based Julia

        tlbef[p] = t_veg_patch[p]
        del2[p] = del_arr[p]

        # Aerodynamic resistances
        ram1_patch[p] = 1.0 / (ustar_patch[p]^2 / um_patch[p])
        rah[p, above_canopy] = 1.0 / (temp1[p] * ustar_patch[p])
        raw[p, above_canopy] = 1.0 / (temp2[p] * ustar_patch[p])

        # Bulk boundary layer resistance
        uaf_patch[p] = um_patch[p] *
            sqrt(1.0 / (ram1_patch[p] * um_patch[p]))

        # Empirical undercanopy wind speed
        uuc[p] = smooth_min(0.4, 0.03 * um_patch[p] / ustar_patch[p])

        # Leaf characteristic width
        dleaf_patch[p] = dleaf_pft[ft]

        cf = cv / (sqrt(uaf_patch[p]) * sqrt(dleaf_patch[p]))
        rb[p] = 1.0 / (cf * uaf_patch[p])
        rb1_patch[p] = rb[p]

        # Soil drag coefficient (X. Zeng parameterization)
        w = exp(-(elai_patch[p] + esai_patch[p]))

        csoilb = vkc / (a_coef * (z0mg_col[c] * uaf_patch[p] / nu_param)^a_exp)

        ri = (grav * htop_patch[p] *
            (taf_patch[p] - t_grnd_col[c])) /
            (taf_patch[p] * uaf_patch[p]^2)

        if use_undercanopy_stability && (taf_patch[p] - t_grnd_col[c]) > 0.0
            ricsoilc = csoilc_val / (1.0 + ria_canopy * smooth_min(ri, 10.0))
            csoilcn = csoilb * w + ricsoilc * (1.0 - w)
        else
            csoilcn = csoilb * w + csoilc_val * (1.0 - w)
        end

        if use_biomass_heat_storage
            rah[p, below_canopy] = 1.0 / (csoilcn * uuc[p])
        else
            rah[p, below_canopy] = 1.0 / (csoilcn * uaf_patch[p])
        end

        raw[p, below_canopy] = rah[p, below_canopy]

        if use_lch4 && length(grnd_ch4_cond_patch) >= p
            grnd_ch4_cond_patch[p] = 1.0 / (raw[p, above_canopy] + raw[p, below_canopy])
        end

        # Stomatal resistance intermediates
        svpts[p] = el[p]
        eah[p] = forc_pbot_col[c] * qaf_patch[p] / 0.622
        rh_af_patch[p] = eah[p] / svpts[p]

        # History outputs
        rah1_patch[p] = rah[p, above_canopy]
        raw1_patch[p] = raw[p, above_canopy]
        rah2_patch[p] = rah[p, below_canopy]
        raw2_patch[p] = raw[p, below_canopy]
        vpd_patch[p]  = smooth_max((svpts[p] - eah[p]), 50.0) * 0.001
    end
end

"""
    cf_resist_update!(frictionvel, canopystate, temperature, waterdiagbulk,
                      patch_data, filterp, fn, temp1, temp2, tlbef, del2,
                      del_arr, uuc, rb, svpts, eah, el, dleaf_pft,
                      grnd_ch4_cond_patch, forc_pbot_col, csoilc_val,
                      use_undercanopy_stability, use_biomass_heat_storage,
                      use_lch4, params)

Launch the canopy-fluxes aerodynamic-resistance inner-loop kernel over the
`fn` filtered patches. Backend-agnostic (CPU loop or GPU); one thread per
filtered patch. Replaces the first `for fi in 1:fn` resistances loop inside the
Newton stability iteration. `csoilc_val` and the feature flags are resolved on
the host and passed as scalars.
"""
function cf_resist_update!(frictionvel, canopystate, temperature, waterdiagbulk,
        patch_data, filterp, fn::Int,
        temp1, temp2, tlbef, del2, del_arr, rah, raw, uuc, rb, svpts, eah, el,
        dleaf_pft, grnd_ch4_cond_patch, forc_pbot_col,
        csoilc_val, use_undercanopy_stability::Bool,
        use_biomass_heat_storage::Bool, use_lch4::Bool, params)
    _launch!(_cf_resist_kernel!, tlbef,
        del2, frictionvel.ram1_patch, rah, raw,
        frictionvel.uaf_patch, uuc, canopystate.dleaf_patch,
        rb, frictionvel.rb1_patch, grnd_ch4_cond_patch, svpts, eah,
        waterdiagbulk.rh_af_patch,
        frictionvel.rah1_patch, frictionvel.raw1_patch,
        frictionvel.rah2_patch, frictionvel.raw2_patch, frictionvel.vpd_patch,
        filterp, patch_data.column, patch_data.gridcell, patch_data.itype,
        temperature.t_veg_patch, del_arr, frictionvel.ustar_patch,
        frictionvel.um_patch, temp1, temp2,
        canopystate.elai_patch, canopystate.esai_patch, canopystate.htop_patch,
        frictionvel.z0mg_col, frictionvel.taf_patch, frictionvel.qaf_patch,
        temperature.t_grnd_col, dleaf_pft, forc_pbot_col, el,
        csoilc_val, use_undercanopy_stability, use_biomass_heat_storage, use_lch4,
        params.cv, params.a_coef, params.a_exp, VKC, GRAV, NU_PARAM, RIA_CANOPY,
        ABOVE_CANOPY, BELOW_CANOPY;
        ndrange = fn)
    return nothing
end

# ---------------------------------------------------------------------------
# Energy-balance inner kernel (per-patch heat-transfer conductances + Newton
# update of leaf temperature + Monin-Obukhov refresh).
#
# Kernelized form of the SECOND per-patch `for fi in 1:fn` loop inside the
# Newton stability iteration of canopy_fluxes! (the largest / most complex
# per-patch loop). One thread per filtered patch; each patch writes only its
# own [p] indices (no reductions). Gather form:
#   fi = @index(Global); p = filterp[fi].
#
# GPU-hostile host-side items resolved before launch and passed as scalars:
#   - do_soilevap_beta() / do_soil_resistance_sl14() are control-flag functions
#     reading module globals; evaluated once on the host as `soilevap_beta` /
#     `soil_resis_sl14` Bool scalars, branched on in-kernel.
#   - length(energyflux.canopy_cond_patch) is passed as `canopy_cond_len::Int`;
#     the methane guard becomes `if use_lch4 && canopy_cond_len >= p`.
#   - Feature flags use_lch4 / use_hydrstress passed as scalar Bools.
#   - frictionvel.zetamaxstable passed as the scalar `zetamaxstable`.
# qsat(...) is kernel-safe and kept as-is (4-tuple destructuring preserved).
# 2D indexing of rah/raw[p, ABOVE_CANOPY|BELOW_CANOPY] and t_soisno_col is
# preserved exactly. Note line ~854 uses smooth_min but line ~867 uses plain
# min — that exact difference is preserved.
# ---------------------------------------------------------------------------
@kernel function _cf_energy_kernel!(
        # written arrays
        wtg, wtl0, wta0, wtstem0, wtga, wtal, lw_stem, lw_leaf,
        wtgq, wtlq0, wtaq0, wtalq, efe, dt_veg, del_arr, err_arr,
        qsatl, el, qsatldT, dth, dqh, delq, obuold, nmozsgn,
        canopy_cond_patch, eflx_sh_stem_patch, eflx_sh_veg_patch,
        qflx_tran_veg_patch, qflx_evap_veg_patch, t_veg_patch,
        taf_patch, qaf_patch, zeta_patch, um_patch, obu_patch,
        # read-only arrays
        @Const(filterp), @Const(column), @Const(gridcell),
        @Const(rah), @Const(raw), @Const(rb), @Const(rstem),
        @Const(sa_leaf), @Const(sa_stem), @Const(sa_internal),
        @Const(frac_rad_abs_by_stem), @Const(air), @Const(bir), @Const(cir),
        @Const(cp_leaf), @Const(tl_ini), @Const(tlbef), @Const(zldis),
        @Const(temp1), @Const(temp2), @Const(ur), @Const(efeb),
        @Const(emv_patch), @Const(t_stem_patch), @Const(btran_patch),
        @Const(fdry_patch), @Const(fwet_patch), @Const(elai_patch),
        @Const(esai_patch), @Const(laisun_patch), @Const(laisha_patch),
        @Const(rssun_patch), @Const(rssha_patch), @Const(qg_col),
        @Const(frac_sno_eff_col), @Const(frac_h2osfc_col), @Const(t_soisno_col),
        @Const(t_h2osfc_col), @Const(snow_depth_col), @Const(soilbeta_col),
        @Const(soilresis_col), @Const(forc_rho_col), @Const(forc_q_col),
        @Const(forc_th_col), @Const(thv_col), @Const(thm_patch),
        @Const(t_grnd_col), @Const(sabv_patch), @Const(snl),
        @Const(frac_veg_nosno_patch), @Const(liqcan_patch), @Const(snocan_patch),
        @Const(uaf_patch), @Const(forc_pbot_col), @Const(ustar_patch),
        # scalars
        soilevap_beta::Bool, soil_resis_sl14::Bool,
        use_lch4::Bool, use_hydrstress::Bool, canopy_cond_len::Int,
        above_canopy::Int, below_canopy::Int,
        sb, cpair, hvap, vkc, grav, btran0, delmax_canopy, zii_canopy,
        beta_canopy, nlevsno::Int, dtime, z_dl, lai_dl, zetamaxstable)

    fi = @index(Global)
    @inbounds begin
        p = filterp[fi]
        c = column[p]
        g = gridcell[p]

        # Sensible heat conductance for air, leaf, ground, stem
        wta  = 1.0 / rah[p, above_canopy]       # air
        wtl  = sa_leaf[p] / rb[p]                # leaf
        wtg[p] = 1.0 / rah[p, below_canopy]     # ground
        wtstem = sa_stem[p] / (rstem[p] + rb[p]) # stem

        wtshi = 1.0 / (wta + wtl + wtstem + wtg[p])

        wtl0[p]    = wtl * wtshi
        wtg0       = wtg[p] * wtshi
        wta0[p]    = wta * wtshi
        wtstem0[p] = wtstem * wtshi
        wtga[p]    = wta0[p] + wtg0 + wtstem0[p]
        wtal[p]    = wta0[p] + wtl0[p] + wtstem0[p]

        # Internal longwave between leaf and stem
        lw_stem[p] = sa_internal[p] * emv_patch[p] * sb * t_stem_patch[p]^4
        lw_leaf[p] = sa_internal[p] * emv_patch[p] * sb * t_veg_patch[p]^4

        # Fraction of potential evaporation from leaf
        fdry_p = fdry_patch[p]
        elai_p = elai_patch[p]
        esai_p = esai_patch[p]
        laisun_p = laisun_patch[p]
        laisha_p = laisha_patch[p]
        rssun_p = rssun_patch[p]
        rssha_p = rssha_patch[p]

        if fdry_p > 0.0
            rppdry = fdry_p * rb[p] *
                (laisun_p / (rb[p] + rssun_p) + laisha_p / (rb[p] + rssha_p)) / elai_p
        else
            rppdry = 0.0
        end

        # Canopy conductance for methane
        if use_lch4 && canopy_cond_len >= p
            canopy_cond_patch[p] = (laisun_p / (rb[p] + rssun_p) +
                laisha_p / (rb[p] + rssha_p)) / smooth_max(elai_p, 0.01)
        end

        efpot = forc_rho_col[c] * ((elai_p + esai_p) / rb[p]) *
            (qsatl[p] - qaf_patch[p])
        h2ocan = liqcan_patch[p] + snocan_patch[p]

        fwet_p = fwet_patch[p]
        btran_p = btran_patch[p]
        qflx_tran_veg_p = qflx_tran_veg_patch[p]

        if use_hydrstress
            if efpot > 0.0
                if btran_p > btran0
                    rpp = rppdry + fwet_p
                else
                    rpp = fwet_p
                end
                rpp = smooth_min(rpp, (qflx_tran_veg_p + h2ocan / dtime) / efpot)
            else
                rpp = 1.0
            end
        else
            if efpot > 0.0
                if btran_p > btran0
                    qflx_tran_veg_patch[p] = efpot * rppdry
                    rpp = rppdry + fwet_p
                else
                    rpp = fwet_p
                    qflx_tran_veg_patch[p] = 0.0
                end
                rpp = min(rpp, (qflx_tran_veg_patch[p] + h2ocan / dtime) / efpot)
            else
                rpp = 1.0
                qflx_tran_veg_patch[p] = 0.0
            end
        end

        # Latent heat conductances
        fvn = frac_veg_nosno_patch[p]
        wtaq  = fvn / raw[p, above_canopy]
        wtlq  = fvn * (elai_p + esai_p) / rb[p] * rpp

        # Litter layer resistance (Sakaguchi)
        snow_depth_c = z_dl
        fsno_dl = snow_depth_col[c] / snow_depth_c
        elai_dl = lai_dl * (1.0 - smooth_min(fsno_dl, 1.0))
        rdl = (1.0 - exp(-elai_dl)) / (0.004 * uaf_patch[p])

        if delq[p] < 0.0
            wtgq[p] = fvn / (raw[p, below_canopy] + rdl)
        else
            if soilevap_beta
                wtgq[p] = soilbeta_col[c] * fvn / (raw[p, below_canopy] + rdl)
            end
            if soil_resis_sl14
                wtgq[p] = fvn / (raw[p, below_canopy] + soilresis_col[c])
            end
        end

        wtsqi  = 1.0 / (wtaq + wtlq + wtgq[p])
        wtgq0  = wtgq[p] * wtsqi
        wtlq0[p] = wtlq * wtsqi
        wtaq0[p] = wtaq * wtsqi
        wtgaq  = wtaq0[p] + wtgq0
        wtalq[p] = wtaq0[p] + wtlq0[p]

        dc1 = forc_rho_col[c] * cpair * wtl
        dc2 = hvap * forc_rho_col[c] * wtlq

        efsh = dc1 * (wtga[p] * t_veg_patch[p] - wtg0 * t_grnd_col[c] -
            wta0[p] * thm_patch[p] - wtstem0[p] * t_stem_patch[p])
        eflx_sh_stem_patch[p] = forc_rho_col[c] * cpair * wtstem *
            ((wta0[p] + wtg0 + wtl0[p]) * t_stem_patch[p] -
             wtg0 * t_grnd_col[c] - wta0[p] * thm_patch[p] -
             wtl0[p] * t_veg_patch[p])
        efe[p] = dc2 * (wtgaq * qsatl[p] - wtgq0 * qg_col[c] -
            wtaq0[p] * forc_q_col[c])

        # Evaporation flux sign change limiter
        erre = 0.0
        if efe[p] * efeb[p] < 0.0
            efeold = efe[p]
            efe[p] = 0.1 * efeold
            erre = efe[p] - efeold
        end

        # Fractionate ground emitted longwave
        snl_c = snl[c]
        lw_grnd = (frac_sno_eff_col[c] * t_soisno_col[c, snl_c + 1 + nlevsno]^4 +
            (1.0 - frac_sno_eff_col[c] - frac_h2osfc_col[c]) *
                t_soisno_col[c, 1 + nlevsno]^4 +
            frac_h2osfc_col[c] * t_h2osfc_col[c]^4)

        # Newton-Raphson: dt_veg
        _numer = ((1.0 - frac_rad_abs_by_stem[p]) * (sabv_patch[p] + air[p] +
            bir[p] * t_veg_patch[p]^4 + cir[p] * lw_grnd) -
            efsh - efe[p] - lw_leaf[p] + lw_stem[p] -
            (cp_leaf[p] / dtime) * (t_veg_patch[p] - tl_ini[p]))
        _denom = ((1.0 - frac_rad_abs_by_stem[p]) * (-4.0 * bir[p] * t_veg_patch[p]^3) +
             4.0 * sa_internal[p] * emv_patch[p] * sb * t_veg_patch[p]^3 +
             dc1 * wtga[p] + dc2 * wtgaq * qsatldT[p] + cp_leaf[p] / dtime)
        dt_veg[p] = _numer / _denom

        t_veg_patch[p] = tlbef[p] + dt_veg[p]

        dels = dt_veg[p]
        del_arr[p] = abs(dels)
        err_arr[p] = 0.0
        if del_arr[p] > delmax_canopy
            dt_veg[p] = delmax_canopy * dels / del_arr[p]
            t_veg_patch[p] = tlbef[p] + dt_veg[p]
            err_arr[p] = (1.0 - frac_rad_abs_by_stem[p]) * (sabv_patch[p] + air[p] +
                bir[p] * tlbef[p]^3 * (tlbef[p] + 4.0 * dt_veg[p]) + cir[p] * lw_grnd) -
                sa_internal[p] * emv_patch[p] * sb * tlbef[p]^3 * (tlbef[p] + 4.0 * dt_veg[p]) +
                lw_stem[p] -
                (efsh + dc1 * wtga[p] * dt_veg[p]) -
                (efe[p] + dc2 * wtgaq * qsatldT[p] * dt_veg[p]) -
                (cp_leaf[p] / dtime) * (t_veg_patch[p] - tl_ini[p])
        end

        # Updated fluxes
        efpot = forc_rho_col[c] * ((elai_p + esai_p) / rb[p]) *
            (wtgaq * (qsatl[p] + qsatldT[p] * dt_veg[p]) -
             wtgq0 * qg_col[c] - wtaq0[p] * forc_q_col[c])
        qflx_evap_veg_patch[p] = rpp * efpot

        # Interception losses / ecidif
        if use_hydrstress
            ecidif = max(0.0, qflx_evap_veg_patch[p] -
                qflx_tran_veg_patch[p] - h2ocan / dtime)
            qflx_evap_veg_patch[p] = min(qflx_evap_veg_patch[p],
                qflx_tran_veg_patch[p] + h2ocan / dtime)
        else
            ecidif = 0.0
            if efpot > 0.0 && btran_patch[p] > btran0
                qflx_tran_veg_patch[p] = efpot * rppdry
            else
                qflx_tran_veg_patch[p] = 0.0
            end
            ecidif = max(0.0, qflx_evap_veg_patch[p] -
                qflx_tran_veg_patch[p] - h2ocan / dtime)
            qflx_evap_veg_patch[p] = min(qflx_evap_veg_patch[p],
                qflx_tran_veg_patch[p] + h2ocan / dtime)
        end

        # Sensible heat from leaves
        eflx_sh_veg_patch[p] = efsh + dc1 * wtga[p] * dt_veg[p] +
            err_arr[p] + erre + hvap * ecidif

        # Update SH and lw_leaf for changes in t_veg
        eflx_sh_stem_patch[p] = eflx_sh_stem_patch[p] +
            forc_rho_col[c] * cpair * wtstem * (-wtl0[p] * dt_veg[p])
        lw_leaf[p] = sa_internal[p] * emv_patch[p] * sb *
            tlbef[p]^3 * (tlbef[p] + 4.0 * dt_veg[p])

        # Re-calculate QSat at updated leaf temperature
        (qs_tmp, es_tmp, dqsdT_tmp, _) = qsat(t_veg_patch[p], forc_pbot_col[c])
        qsatl[p] = qs_tmp
        el[p] = es_tmp
        qsatldT[p] = dqsdT_tmp

        # Update canopy air temperature and humidity
        taf_patch[p] = wtg0 * t_grnd_col[c] +
            wta0[p] * thm_patch[p] +
            wtl0[p] * t_veg_patch[p] +
            wtstem0[p] * t_stem_patch[p]
        qaf_patch[p] = wtlq0[p] * qsatl[p] +
            wtgq0 * qg_col[c] +
            forc_q_col[c] * wtaq0[p]

        # Update Monin-Obukhov length and wind speed
        dth[p] = thm_patch[p] - taf_patch[p]
        dqh[p] = forc_q_col[c] - qaf_patch[p]
        delq[p] = wtalq[p] * qg_col[c] - wtlq0[p] * qsatl[p] - wtaq0[p] * forc_q_col[c]

        tstar = temp1[p] * dth[p]
        qstar = temp2[p] * dqh[p]
        thvstar = tstar * (1.0 + 0.61 * forc_q_col[c]) + 0.61 * forc_th_col[c] * qstar
        zeta_patch[p] = zldis[p] * vkc * grav * thvstar /
            (ustar_patch[p]^2 * thv_col[c])

        if zeta_patch[p] >= 0.0  # stable
            zeta_patch[p] = smooth_clamp(zeta_patch[p], 0.01, zetamaxstable)
            um_patch[p] = smooth_max(ur[p], 0.1)
        else  # unstable
            zeta_patch[p] = smooth_clamp(zeta_patch[p], -100.0, -0.01)
            if ustar_patch[p] * thvstar > 0.0
                wc = 0.0
            else
                wc_arg = smooth_max(-grav * ustar_patch[p] * thvstar *
                    zii_canopy / thv_col[c], 0.0)
                wc = beta_canopy * wc_arg^0.333
            end
            um_patch[p] = sqrt(ur[p]^2 + wc^2)
        end
        obu_patch[p] = zldis[p] / zeta_patch[p]

        if obuold[p] * obu_patch[p] < 0.0
            nmozsgn[p] += 1
        end
        if nmozsgn[p] >= 4
            obu_patch[p] = zldis[p] / (-0.01)
        end
        obuold[p] = obu_patch[p]
    end
end

"""
    cf_energy_update!(canopystate, energyflux, frictionvel, temperature,
                      solarabs, soilstate, waterfluxbulk, waterstatebulk,
                      waterdiagbulk, photosyns, patch_data, col_data, filterp,
                      fn, rah, raw, rb, rstem, sa_leaf, sa_stem, sa_internal,
                      frac_rad_abs_by_stem, air, bir, cir, cp_leaf, tl_ini,
                      tlbef, zldis, temp1, temp2, ur, efeb, wtg, wtl0, wta0,
                      wtstem0, wtga, wtal, lw_stem, lw_leaf, wtgq, wtlq0, wtaq0,
                      wtalq, efe, dt_veg, del_arr, err_arr, qsatl, el, qsatldT,
                      dth, dqh, delq, obuold, nmozsgn, qg_col, forc_rho_col,
                      forc_q_col, forc_th_col, forc_pbot_col,
                      soilevap_beta, soil_resis_sl14, use_lch4, use_hydrstress,
                      grnd_ch4_cond_patch, dtime, params)

Launch the canopy-fluxes energy-balance inner-loop kernel over the `fn`
filtered patches. Backend-agnostic (CPU loop or GPU); one thread per filtered
patch. Replaces the second `for fi in 1:fn` heat-transfer / Newton-update loop
inside the Newton stability iteration. The two control-flag functions
(`do_soilevap_beta`, `do_soil_resistance_sl14`), the methane conductance
array length guard, and the feature flags are resolved on the host and passed
as scalars; `frictionvel.zetamaxstable` is passed as a scalar value.
"""
function cf_energy_update!(canopystate, energyflux, frictionvel, temperature,
        solarabs, soilstate, waterfluxbulk, waterstatebulk, waterdiagbulk,
        photosyns, patch_data, col_data, filterp, fn::Int,
        rah, raw, rb, rstem, sa_leaf, sa_stem, sa_internal,
        frac_rad_abs_by_stem, air, bir, cir, cp_leaf, tl_ini, tlbef, zldis,
        temp1, temp2, ur, efeb, wtg, wtl0, wta0, wtstem0, wtga, wtal,
        lw_stem, lw_leaf, wtgq, wtlq0, wtaq0, wtalq, efe, dt_veg, del_arr,
        err_arr, qsatl, el, qsatldT, dth, dqh, delq, obuold, nmozsgn,
        forc_q_col, forc_th_col, forc_pbot_col, forc_rho_col,
        soilevap_beta::Bool, soil_resis_sl14::Bool, use_lch4::Bool,
        use_hydrstress::Bool, nlevsno::Int, dtime, params)
    _launch!(_cf_energy_kernel!,
        # written arrays
        wtg, wtl0, wta0, wtstem0, wtga, wtal, lw_stem, lw_leaf,
        wtgq, wtlq0, wtaq0, wtalq, efe, dt_veg, del_arr, err_arr,
        qsatl, el, qsatldT, dth, dqh, delq, obuold, nmozsgn,
        energyflux.canopy_cond_patch, energyflux.eflx_sh_stem_patch,
        energyflux.eflx_sh_veg_patch, waterfluxbulk.wf.qflx_tran_veg_patch,
        waterfluxbulk.wf.qflx_evap_veg_patch, temperature.t_veg_patch,
        frictionvel.taf_patch, frictionvel.qaf_patch, frictionvel.zeta_patch,
        frictionvel.um_patch, frictionvel.obu_patch,
        # read-only arrays
        filterp, patch_data.column, patch_data.gridcell,
        rah, raw, rb, rstem, sa_leaf, sa_stem, sa_internal,
        frac_rad_abs_by_stem, air, bir, cir, cp_leaf, tl_ini, tlbef, zldis,
        temp1, temp2, ur, efeb,
        temperature.emv_patch, temperature.t_stem_patch, energyflux.btran_patch,
        waterdiagbulk.fdry_patch, waterdiagbulk.fwet_patch,
        canopystate.elai_patch, canopystate.esai_patch,
        canopystate.laisun_patch, canopystate.laisha_patch,
        photosyns.rssun_patch, photosyns.rssha_patch, waterdiagbulk.qg_col,
        waterdiagbulk.frac_sno_eff_col, waterdiagbulk.frac_h2osfc_col,
        temperature.t_soisno_col, temperature.t_h2osfc_col,
        waterdiagbulk.snow_depth_col, soilstate.soilbeta_col,
        soilstate.soilresis_col, forc_rho_col, forc_q_col,
        forc_th_col, temperature.thv_col, temperature.thm_patch,
        temperature.t_grnd_col, solarabs.sabv_patch, col_data.snl,
        canopystate.frac_veg_nosno_patch, waterstatebulk.ws.liqcan_patch,
        waterstatebulk.ws.snocan_patch, frictionvel.uaf_patch,
        forc_pbot_col, frictionvel.ustar_patch,
        # scalars
        soilevap_beta, soil_resis_sl14, use_lch4, use_hydrstress,
        length(energyflux.canopy_cond_patch), ABOVE_CANOPY, BELOW_CANOPY,
        SB, CPAIR, HVAP, VKC, GRAV, BTRAN0, DELMAX_CANOPY, ZII_CANOPY,
        BETA_CANOPY, nlevsno, dtime, params.z_dl, params.lai_dl,
        frictionvel.zetamaxstable;
        ndrange = fn)
    return nothing
end

# =====================================================================
# canopy_fluxes! Phase-1 init kernels (cf1). Per-filtered-patch (fi → p=filterp[fi]),
# elementwise; literals eltype-converted. (canopy_fluxes! does not run whole-function
# on GPU — the Newton iteration uses host-side filter compaction — so these add to
# its incremental GPU coverage; the expensive Newton physics is already kernelized.)
# =====================================================================
@kernel function _cf_init_zero_kernel!(del_arr, efeb, wtlq0, wtalq, wtgq, wtaq0, obuold,
        dhsdt_canopy, eflx_sh_stem, @Const(filterp))
    fi = @index(Global)
    @inbounds begin
        p = filterp[fi]
        T = eltype(del_arr)
        del_arr[p] = zero(T); efeb[p] = zero(T); wtlq0[p] = zero(T); wtalq[p] = zero(T)
        wtgq[p] = zero(T); wtaq0[p] = zero(T); obuold[p] = zero(T)
        dhsdt_canopy[p] = zero(T); eflx_sh_stem[p] = zero(T)
    end
end

@kernel function _cf_daylength_kernel!(dayl_factor, @Const(filterp), @Const(gridcell),
        @Const(dayl_grc), @Const(max_dayl_grc))
    fi = @index(Global)
    @inbounds begin
        p = filterp[fi]
        g = gridcell[p]
        T = eltype(dayl_factor)
        dayl_factor[p] = smooth_clamp(
            (dayl_grc[g] * dayl_grc[g]) / (max_dayl_grc[g] * max_dayl_grc[g]), T(0.01), one(T))
    end
end

# Biomass heat capacities (cf2): per-filtered-patch, PFT-parameter-indexed (ft =
# itype[p]+1). use_cn passed as a Bool; constants eltype-converted.
@kernel function _cf_biomass_heat_kernel!(frac_rad_abs_by_stem, dbh, sa_leaf, sa_stem,
        sa_internal, cp_leaf, cp_stem, rstem, leaf_biomass, stem_biomass,
        @Const(filterp), @Const(itype), @Const(elai), @Const(esai), @Const(htop),
        @Const(fbw_pft), @Const(nstem_pft), @Const(wood_density_pft), @Const(dbh_pft),
        @Const(slatop_pft), @Const(rstem_per_dbh_pft), @Const(is_tree_pft),
        @Const(is_shrub_pft), use_cn::Bool)
    fi = @index(Global)
    @inbounds begin
        p = filterp[fi]
        T = eltype(sa_leaf)
        ft = itype[p] + 1
        elai_p = elai[p]; esai_p = esai[p]
        if (elai_p + esai_p) > zero(T)
            frac_rad_abs_by_stem[p] = esai_p / (elai_p + esai_p)
        else
            frac_rad_abs_by_stem[p] = zero(T)
        end
        if elai_p > zero(T)
            frac_rad_abs_by_stem[p] = T(K_VERT_CANOPY) * frac_rad_abs_by_stem[p]
        end
        if use_cn
            if stem_biomass[p] > zero(T)
                dbh[p] = T(2.0) * sqrt(stem_biomass[p] * (one(T) - fbw_pft[ft]) /
                    (T(RPI) * htop[p] * T(K_CYL_VOL_CANOPY) * nstem_pft[ft] * wood_density_pft[ft]))
            else
                dbh[p] = zero(T)
            end
        else
            dbh[p] = dbh_pft[ft]
        end
        sa_leaf[p] = T(2.0) * elai_p
        sa_stem[p] = nstem_pft[ft] * htop[p] * T(RPI) * dbh[p]
        sa_stem[p] = T(K_CYL_AREA_CANOPY) * sa_stem[p]
        if !(is_tree_pft[ft] || is_shrub_pft[ft]) || dbh[p] < T(MIN_STEM_DIAMETER)
            frac_rad_abs_by_stem[p] = zero(T)
            sa_stem[p] = zero(T)
            sa_leaf[p] = sa_leaf[p] + esai_p
        end
        if !use_cn
            leaf_biomass[p] = (T(1.0e-3) * T(C_TO_B) / slatop_pft[ft]) *
                smooth_max(T(0.01), T(0.5) * sa_leaf[p]) / (one(T) - fbw_pft[ft])
            carea_stem = T(RPI) * (dbh[p] * T(0.5))^2
            stem_biomass[p] = carea_stem * htop[p] * T(K_CYL_VOL_CANOPY) *
                nstem_pft[ft] * wood_density_pft[ft] / (one(T) - fbw_pft[ft])
        end
        sa_internal[p] = min(sa_leaf[p], sa_stem[p])
        sa_internal[p] = T(K_INTERNAL_CANOPY) * sa_internal[p]
        cp_leaf[p] = leaf_biomass[p] * (T(C_DRY_BIOMASS) * (one(T) - fbw_pft[ft]) + fbw_pft[ft] * T(C_WATER))
        cp_stem[p] = stem_biomass[p] * (T(C_DRY_BIOMASS) * (one(T) - fbw_pft[ft]) + fbw_pft[ft] * T(C_WATER))
        cp_stem[p] = T(K_CYL_VOL_CANOPY) * cp_stem[p]
        rstem[p] = rstem_per_dbh_pft[ft] * dbh[p]
    end
end

@kernel function _cf_no_biomass_kernel!(frac_rad_abs_by_stem, sa_leaf, sa_stem,
        sa_internal, cp_leaf, cp_stem, rstem, @Const(filterp), @Const(elai), @Const(esai))
    fi = @index(Global)
    @inbounds begin
        p = filterp[fi]
        T = eltype(sa_leaf)
        sa_leaf[p] = elai[p] + esai[p]
        frac_rad_abs_by_stem[p] = zero(T)
        sa_stem[p] = zero(T); sa_internal[p] = zero(T)
        cp_leaf[p] = zero(T); cp_stem[p] = zero(T); rstem[p] = zero(T)
    end
end

# =====================================================================
# Main canopy_fluxes! function
# =====================================================================

"""
    canopy_fluxes!(canopystate, energyflux, frictionvel, temperature,
                   solarabs, soilstate, waterfluxbulk, waterstatebulk,
                   waterdiagbulk, photosyns, patch_data, col_data, gridcell_data,
                   mask_exposedvegp, bounds_patch, bounds_col,
                   forc_lwrad_col, forc_q_col, forc_pbot_col, forc_th_col,
                   forc_rho_col, forc_t_col, forc_u_grc, forc_v_grc,
                   forc_pco2_grc, forc_po2_grc, forc_hgt_t_grc, forc_hgt_u_grc,
                   forc_hgt_q_grc, dayl_grc, max_dayl_grc,
                   downreg_patch, leafn_patch, dtime;
                   kwargs...)

Main canopy fluxes calculation. Solves the leaf energy balance iteratively
using Newton-Raphson to find leaf temperature (t_veg).

Ported from `subroutine CanopyFluxes` in `CanopyFluxesMod.F90`.
"""
function canopy_fluxes!(
        canopystate     ::CanopyStateData,
        energyflux      ::EnergyFluxData,
        frictionvel     ::FrictionVelocityData,
        temperature     ::TemperatureData,
        solarabs        ::SolarAbsorbedData,
        soilstate       ::SoilStateData,
        waterfluxbulk   ::WaterFluxBulkData,
        waterstatebulk  ::WaterStateBulkData,
        waterdiagbulk   ::WaterDiagnosticBulkData,
        photosyns       ::PhotosynthesisData,
        patch_data      ::PatchData,
        col_data        ::ColumnData,
        gridcell_data   ::GridcellData,
        mask_exposedvegp::BitVector,
        bounds_patch    ::UnitRange{Int},
        bounds_col      ::UnitRange{Int},
        # Atmospheric forcing (column-level)
        forc_lwrad_col  ::Vector{<:Real},
        forc_q_col      ::Vector{<:Real},
        forc_pbot_col   ::Vector{<:Real},
        forc_th_col     ::Vector{<:Real},
        forc_rho_col    ::Vector{<:Real},
        forc_t_col      ::Vector{<:Real},
        # Atmospheric forcing (gridcell-level)
        forc_u_grc      ::Vector{<:Real},
        forc_v_grc      ::Vector{<:Real},
        forc_pco2_grc   ::Vector{<:Real},
        forc_po2_grc    ::Vector{<:Real},
        forc_hgt_t_grc  ::Vector{<:Real},
        forc_hgt_u_grc  ::Vector{<:Real},
        forc_hgt_q_grc  ::Vector{<:Real},
        # Daylength (gridcell-level)
        dayl_grc        ::Vector{<:Real},
        max_dayl_grc    ::Vector{<:Real},
        # CN inputs (patch-level)
        downreg_patch   ::Vector{<:Real},
        leafn_patch     ::Vector{<:Real},
        # Time step
        dtime           ::Real;
        # PFT parameters (patch-indexed via ivt)
        dleaf_pft       ::Vector{<:Real} = fill(0.04, MXPFT+1),
        dbh_pft         ::Vector{<:Real} = fill(0.1, MXPFT+1),
        slatop_pft      ::Vector{<:Real} = fill(0.01, MXPFT+1),
        fbw_pft         ::Vector{<:Real} = fill(0.1, MXPFT+1),
        nstem_pft       ::Vector{<:Real} = fill(1.0, MXPFT+1),
        woody_pft       ::Vector{<:Real} = fill(0.0, MXPFT+1),
        rstem_per_dbh_pft::Vector{<:Real} = fill(0.0, MXPFT+1),
        wood_density_pft::Vector{<:Real} = fill(500.0, MXPFT+1),
        is_tree_pft     ::Vector{Bool} = fill(false, MXPFT+1),
        is_shrub_pft    ::Vector{Bool} = fill(false, MXPFT+1),
        c3psn_pft       ::Vector{<:Real} = fill(1.0, MXPFT+1),
        leafcn_pft      ::Vector{<:Real} = fill(25.0, MXPFT+1),
        flnr_pft        ::Vector{<:Real} = fill(0.08, MXPFT+1),
        fnitr_pft       ::Vector{<:Real} = fill(0.1, MXPFT+1),
        mbbopt_pft      ::Vector{<:Real} = fill(0.0, MXPFT+1),
        medlynintercept_pft::Vector{<:Real} = fill(100.0, MXPFT+1),
        medlynslope_pft ::Vector{<:Real} = fill(6.0, MXPFT+1),
        crop_pft        ::Vector{<:Real} = Float64[],
        z0v_Cr_pft      ::Vector{<:Real} = fill(0.35, MXPFT+1),
        z0v_Cs_pft      ::Vector{<:Real} = fill(0.003, MXPFT+1),
        z0v_c_pft       ::Vector{<:Real} = fill(0.25, MXPFT+1),
        z0v_cw_pft      ::Vector{<:Real} = fill(2.0, MXPFT+1),
        z0v_LAImax_pft  ::Vector{<:Real} = fill(8.0, MXPFT+1),
        # Feature flags
        use_cn          ::Bool = false,
        use_lch4        ::Bool = false,
        use_c13         ::Bool = false,
        use_hydrstress  ::Bool = false,
        use_fates       ::Bool = false,
        use_luna        ::Bool = false,
        z0param_method  ::String = "ZengWang2007",
        # Ozone stress arrays (optional)
        o3coefv_patch   ::Vector{<:Real} = Float64[],
        o3coefg_patch   ::Vector{<:Real} = Float64[],
        # CH4 conductance output (optional)
        grnd_ch4_cond_patch::Vector{<:Real} = Float64[],
        # Photosynthesis extras
        forc_pc13o2_grc ::Vector{<:Real} = Float64[],
        t10_patch       ::Vector{<:Real} = Float64[],
        nrad_patch      ::Vector{Int} = Int[],
        tlai_z_patch    ::Matrix{<:Real} = Matrix{Float64}(undef,0,0),
        vcmaxcint_sun_patch ::Vector{<:Real} = Float64[],
        vcmaxcint_sha_patch ::Vector{<:Real} = Float64[],
        parsun_z_patch  ::Matrix{<:Real} = Matrix{Float64}(undef,0,0),
        parsha_z_patch  ::Matrix{<:Real} = Matrix{Float64}(undef,0,0),
        laisun_z_patch  ::Matrix{<:Real} = Matrix{Float64}(undef,0,0),
        laisha_z_patch  ::Matrix{<:Real} = Matrix{Float64}(undef,0,0),
        leaf_mr_vcm     ::Real = 0.015,
        # Calibration overrides (NaN = use defaults)
        overrides       ::CalibrationOverrides = CalibrationOverrides())

    np = length(bounds_patch)
    begp = first(bounds_patch)
    endp = last(bounds_patch)
    nlevgrnd = varpar.nlevgrnd
    nlevsno  = varpar.nlevsno

    # --- Convenience aliases ---
    params = canopy_fluxes_params
    ctrl   = canopy_fluxes_ctrl

    # Local work arrays (patch-indexed)
    FT = eltype(forc_t_col)
    zldis     = zeros(FT, endp)
    dth       = zeros(FT, endp)
    dthv      = zeros(FT, endp)
    dqh       = zeros(FT, endp)
    ur        = zeros(FT, endp)
    temp1     = zeros(FT, endp)
    temp12m   = zeros(FT, endp)
    temp2     = zeros(FT, endp)
    temp22m   = zeros(FT, endp)
    rb        = zeros(FT, endp)
    rah       = zeros(FT, endp, 2)
    raw       = zeros(FT, endp, 2)
    wtg       = zeros(FT, endp)
    wta0      = zeros(FT, endp)
    wtl0      = zeros(FT, endp)
    wtstem0   = zeros(FT, endp)
    wtal      = zeros(FT, endp)
    wtga      = zeros(FT, endp)
    wtgq      = zeros(FT, endp)
    wtaq0     = zeros(FT, endp)
    wtlq0     = zeros(FT, endp)
    wtalq     = zeros(FT, endp)
    el        = zeros(FT, endp)
    qsatl     = zeros(FT, endp)
    qsatldT   = zeros(FT, endp)
    air       = zeros(FT, endp)
    bir       = zeros(FT, endp)
    cir       = zeros(FT, endp)
    del_arr   = zeros(FT, endp)  # "del" in Fortran (reserved word in Julia context for clarity)
    del2      = zeros(FT, endp)
    dele      = zeros(FT, endp)
    delq      = zeros(FT, endp)
    det_arr   = zeros(FT, endp)
    efeb      = zeros(FT, endp)
    efe       = zeros(FT, endp)
    err_arr   = zeros(FT, endp)
    obuold    = zeros(FT, endp)
    tlbef     = zeros(FT, endp)
    tl_ini    = zeros(FT, endp)
    ts_ini    = zeros(FT, endp)
    co2_arr   = zeros(FT, endp)
    o2_arr    = zeros(FT, endp)
    svpts     = zeros(FT, endp)
    eah       = zeros(FT, endp)
    dt_veg    = zeros(FT, endp)
    fm        = zeros(FT, endp)
    nmozsgn   = zeros(Int, endp)
    dayl_factor = zeros(FT, endp)
    rootsum   = zeros(FT, endp)
    dbh       = zeros(FT, endp)
    cp_leaf   = zeros(FT, endp)
    cp_stem   = zeros(FT, endp)
    rstem     = zeros(FT, endp)
    dt_stem   = zeros(FT, endp)
    frac_rad_abs_by_stem = zeros(FT, endp)
    lw_stem   = zeros(FT, endp)
    lw_leaf   = zeros(FT, endp)
    sa_stem   = zeros(FT, endp)
    sa_leaf   = zeros(FT, endp)
    sa_internal = zeros(FT, endp)
    uuc       = zeros(FT, endp)
    snocan_baseline = zeros(FT, endp)

    # Integer filter arrays (Fortran-style, for convergence testing)
    filterp = zeros(Int, np)
    fporig  = zeros(Int, np)

    # Build initial filter from mask
    fn = 0
    for p in bounds_patch
        if mask_exposedvegp[p]
            fn += 1
            filterp[fn] = p
        end
    end

    # Time step initialization of photosynthesis variables
    photosynthesis_timestep_init!(photosyns, mask_exposedvegp, bounds_patch)

    # =========================================================================
    # Phase 1: Initialization
    # =========================================================================

    _launch!(_cf_init_zero_kernel!, del_arr, efeb, wtlq0, wtalq, wtgq, wtaq0, obuold,
        energyflux.dhsdt_canopy_patch, energyflux.eflx_sh_stem_patch, filterp; ndrange = fn)

    # --- Calculate biomass heat capacities (per-filtered-patch kernels) ---
    if ctrl.use_biomass_heat_storage
        _launch!(_cf_biomass_heat_kernel!, frac_rad_abs_by_stem, dbh, sa_leaf, sa_stem,
            sa_internal, cp_leaf, cp_stem, rstem, canopystate.leaf_biomass_patch,
            canopystate.stem_biomass_patch, filterp, patch_data.itype,
            canopystate.elai_patch, canopystate.esai_patch, canopystate.htop_patch,
            fbw_pft, nstem_pft, wood_density_pft, dbh_pft, slatop_pft, rstem_per_dbh_pft,
            is_tree_pft, is_shrub_pft, use_cn; ndrange = fn)
    else
        _launch!(_cf_no_biomass_kernel!, frac_rad_abs_by_stem, sa_leaf, sa_stem,
            sa_internal, cp_leaf, cp_stem, rstem, filterp, canopystate.elai_patch,
            canopystate.esai_patch; ndrange = fn)
    end

    # --- Daylength control for Vcmax ---
    _launch!(_cf_daylength_kernel!, dayl_factor, filterp, patch_data.gridcell,
        dayl_grc, max_dayl_grc; ndrange = fn)

    # Zero boundary layer resistance (broadcast over the patch range; backend-agnostic).
    @views frictionvel.rb1_patch[bounds_patch] .= zero(eltype(frictionvel.rb1_patch))

    # --- Compute effective soil porosity and volumetric liquid water ---
    # (These are done inline rather than calling separate functions,
    #  since the separate module functions have different signatures)

    # --- Set perchroot options ---
    set_perchroot_opt!(ctrl.perchroot, ctrl.perchroot_alt)

    # --- Root moisture stress ---
    # (Simplified: set btran=1 and rresis=1 if detailed calc not available.
    #  The full calc_root_moist_stress call would go here.)
    # In a full implementation, this would call calc_root_moist_stress.

    # --- Modify aerodynamic parameters for sparse/dense canopy (X. Zeng) ---
    for fi in 1:fn
        p = filterp[fi]
        c = patch_data.column[p]
        g = patch_data.gridcell[p]
        ft = patch_data.itype[p] + 1  # 0-based Fortran PFT → 1-based Julia

        if z0param_method == "ZengWang2007"
            lt = smooth_min(canopystate.elai_patch[p] + canopystate.esai_patch[p], TLSAI_CRIT)
            egvf = (1.0 - ALPHA_AERO * exp(-lt)) / (1.0 - ALPHA_AERO * exp(-TLSAI_CRIT))
            canopystate.displa_patch[p] = egvf * canopystate.displa_patch[p]
            frictionvel.z0mv_patch[p] = exp(egvf * log(frictionvel.z0mv_patch[p]) +
                (1.0 - egvf) * log(frictionvel.z0mg_col[c]))

        elseif z0param_method == "Meier2022"
            lt = smooth_max(1.0e-5, canopystate.elai_patch[p] + canopystate.esai_patch[p])
            canopystate.displa_patch[p] = canopystate.htop_patch[p] *
                (1.0 - (1.0 - exp(-(CD1_PARAM * lt)^0.5)) / (CD1_PARAM * lt)^0.5)

            lt = smooth_min(lt, z0v_LAImax_pft[ft])
            delt_iter = 2.0
            U_ustar_ini = (z0v_Cs_pft[ft] + z0v_Cr_pft[ft] * lt * 0.5)^(-0.5) *
                z0v_c_pft[ft] * lt * 0.25
            U_ustar = U_ustar_ini
            while delt_iter > 1.0e-4
                U_ustar_prev = U_ustar
                U_ustar = U_ustar_ini * exp(U_ustar_prev)
                delt_iter = abs(U_ustar - U_ustar_prev)
            end
            U_ustar = 4.0 * U_ustar / lt / z0v_c_pft[ft]

            frictionvel.z0mv_patch[p] = canopystate.htop_patch[p] *
                (1.0 - canopystate.displa_patch[p] / canopystate.htop_patch[p]) *
                exp(-VKC * U_ustar + log(z0v_cw_pft[ft]) - 1.0 + z0v_cw_pft[ft]^(-1.0))
        else
            error("canopy_fluxes!: unknown z0param_method: $z0param_method")
        end

        frictionvel.z0hv_patch[p] = frictionvel.z0mv_patch[p]
        frictionvel.z0qv_patch[p] = frictionvel.z0mv_patch[p]

        # Update forcing heights
        frictionvel.forc_hgt_u_patch[p] = forc_hgt_u_grc[g] + frictionvel.z0mv_patch[p] + canopystate.displa_patch[p]
        frictionvel.forc_hgt_t_patch[p] = forc_hgt_t_grc[g] + frictionvel.z0hv_patch[p] + canopystate.displa_patch[p]
        frictionvel.forc_hgt_q_patch[p] = forc_hgt_q_grc[g] + frictionvel.z0qv_patch[p] + canopystate.displa_patch[p]
    end

    # --- Net absorbed longwave radiation, QSat, CO2/O2, flux initialization ---
    found = false
    found_index = 0
    for fi in 1:fn
        p = filterp[fi]
        c = patch_data.column[p]
        g = patch_data.gridcell[p]

        emv_p = temperature.emv_patch[p]
        emg_c = temperature.emg_col[c]

        # Net absorbed longwave = air + bir*t_veg^4 + cir*t_grnd^4
        air[p] =  emv_p * (1.0 + (1.0 - emv_p) * (1.0 - emg_c)) * forc_lwrad_col[c]
        bir[p] = -(2.0 - emv_p * (1.0 - emg_c)) * emv_p * SB
        cir[p] =  emv_p * emg_c * SB

        # Saturated vapor pressure at leaf surface
        (qs_tmp, es_tmp, dqsdT_tmp, _) = qsat(temperature.t_veg_patch[p], forc_pbot_col[c])
        qsatl[p] = qs_tmp
        el[p] = es_tmp
        qsatldT[p] = dqsdT_tmp

        # Atmospheric CO2 and O2
        co2_arr[p] = forc_pco2_grc[g]
        o2_arr[p]  = forc_po2_grc[g]

        # Initialize flux profile
        nmozsgn[p] = 0
        frictionvel.taf_patch[p] = (temperature.t_grnd_col[c] + temperature.thm_patch[p]) / 2.0
        frictionvel.qaf_patch[p] = (forc_q_col[c] + waterdiagbulk.qg_col[c]) / 2.0

        ur[p] = smooth_max(params.wind_min, sqrt(forc_u_grc[g]^2 + forc_v_grc[g]^2))
        dth[p] = temperature.thm_patch[p] - frictionvel.taf_patch[p]
        dqh[p] = forc_q_col[c] - frictionvel.qaf_patch[p]
        delq[p] = waterdiagbulk.qg_col[c] - frictionvel.qaf_patch[p]
        dthv[p] = dth[p] * (1.0 + 0.61 * forc_q_col[c]) + 0.61 * forc_th_col[c] * dqh[p]
        zldis[p] = frictionvel.forc_hgt_u_patch[p] - canopystate.displa_patch[p]

        if zldis[p] < 0.0
            found = true
            found_index = p
        end
    end

    if found && !use_fates
        @warn "canopy_fluxes!: forcing height below canopy height at patch $found_index"
    end

    # --- Initialize Monin-Obukhov length ---
    for fi in 1:fn
        p = filterp[fi]
        c = patch_data.column[p]

        (um_val, obu_val) = monin_obuk_ini(frictionvel.zetamaxstable,
            ur[p], temperature.thv_col[c], dthv[p], zldis[p], frictionvel.z0mv_patch[p])
        frictionvel.um_patch[p] = um_val
        frictionvel.obu_patch[p] = obu_val
        frictionvel.num_iter_patch[p] = 0.0

        tl_ini[p] = temperature.t_veg_patch[p]
        ts_ini[p] = temperature.t_stem_patch[p]
    end

    # --- Save original filter for post-iteration ---
    itlef = 0
    fnorig = fn
    fporig[1:fn] .= filterp[1:fn]

    # =========================================================================
    # Phase 2: Stability iteration (Newton-Raphson)
    # =========================================================================

    # Soil drag coefficient: use calibration override if set, else param default.
    # Constant across patches — resolve once on the host (a kernel cannot call
    # isnan(overrides.csoilc) on the host struct).
    csoilc_val = isnan(overrides.csoilc) ? params.csoilc : overrides.csoilc

    # Control-flag functions read module globals — not GPU-safe. Resolve once on
    # the host to Bool scalars and pass into the energy-balance kernel.
    soilevap_beta   = do_soilevap_beta()
    soil_resis_sl14 = do_soil_resistance_sl14()

    while itlef <= ctrl.itmax_canopy_fluxes && fn > 0

        # --- Friction velocity and boundary layer profiles ---
        # Build filter arrays for friction_velocity! call
        filt_arr = filterp[1:fn]
        disp_vec  = canopystate.displa_patch
        z0m_vec   = frictionvel.z0mv_patch
        z0h_vec   = frictionvel.z0hv_patch
        z0q_vec   = frictionvel.z0qv_patch
        obu_vec   = frictionvel.obu_patch
        um_vec    = frictionvel.um_patch
        ustar_vec = frictionvel.ustar_patch

        friction_velocity!(frictionvel, fn, filt_arr,
            disp_vec, z0m_vec, z0h_vec, z0q_vec,
            obu_vec, itlef + 1, ur, um_vec, ustar_vec,
            temp1, temp2, temp12m, temp22m, fm)

        # First per-patch loop: aerodynamic resistances / conductances.
        # Kernelized (one thread per filtered patch); csoilc_val and feature
        # flags resolved on the host and passed as scalars.
        cf_resist_update!(frictionvel, canopystate, temperature, waterdiagbulk,
            patch_data, filterp, fn,
            temp1, temp2, tlbef, del2, del_arr, rah, raw, uuc, rb, svpts, eah, el,
            dleaf_pft, grnd_ch4_cond_patch, forc_pbot_col,
            csoilc_val, ctrl.use_undercanopy_stability,
            ctrl.use_biomass_heat_storage, use_lch4, params)

        # --- Photosynthesis ---
        # Call photosynthesis for sunlit and shaded leaves to update rssun/rssha
        if !isempty(nrad_patch) && !isempty(parsun_z_patch)
            # Build patch-indexed forc_pbot from column-level data
            # Use full mask (not reduced filterp) since photosynthesis uses mask_exposedvegp
            forc_pbot_patch = zeros(FT, endp)
            for p in bounds_patch
                mask_exposedvegp[p] || continue
                c = patch_data.column[p]
                forc_pbot_patch[p] = forc_pbot_col[c]
            end

            # PFT index vector (+1 for 0-based Fortran → 1-based Julia)
            ivt_vec = patch_data.itype .+ 1

            # Sunlit leaves
            photosynthesis!(photosyns,
                svpts, eah, o2_arr, co2_arr, rb,
                energyflux.btran_patch, dayl_factor, leafn_patch,
                forc_pbot_patch, temperature.t_veg_patch, t10_patch,
                temperature.thm_patch, nrad_patch,
                tlai_z_patch, canopystate.tlai_patch,
                parsun_z_patch, laisun_z_patch,
                vcmaxcint_sun_patch,
                o3coefv_patch, o3coefg_patch,
                c3psn_pft, leafcn_pft, flnr_pft, fnitr_pft, slatop_pft,
                mbbopt_pft, medlynintercept_pft, medlynslope_pft,
                ivt_vec, patch_data.column,
                mask_exposedvegp, bounds_patch, "sun";
                use_cn=use_cn, use_luna=use_luna, use_c13=use_c13,
                leaf_mr_vcm=leaf_mr_vcm, crop_pft=crop_pft,
                overrides=overrides)

            # Shaded leaves
            photosynthesis!(photosyns,
                svpts, eah, o2_arr, co2_arr, rb,
                energyflux.btran_patch, dayl_factor, leafn_patch,
                forc_pbot_patch, temperature.t_veg_patch, t10_patch,
                temperature.thm_patch, nrad_patch,
                tlai_z_patch, canopystate.tlai_patch,
                parsha_z_patch, laisha_z_patch,
                vcmaxcint_sha_patch,
                o3coefv_patch, o3coefg_patch,
                c3psn_pft, leafcn_pft, flnr_pft, fnitr_pft, slatop_pft,
                mbbopt_pft, medlynintercept_pft, medlynslope_pft,
                ivt_vec, patch_data.column,
                mask_exposedvegp, bounds_patch, "sha";
                use_cn=use_cn, use_luna=use_luna, use_c13=use_c13,
                leaf_mr_vcm=leaf_mr_vcm, crop_pft=crop_pft,
                overrides=overrides)

        end

        # --- Heat transfer conductances + Newton update of leaf temperature ---
        # Second per-patch loop: kernelized (one thread per filtered patch).
        # The two control-flag functions (do_soilevap_beta /
        # do_soil_resistance_sl14), the methane conductance length guard, and
        # the feature flags are resolved on the host and passed as scalars.
        cf_energy_update!(canopystate, energyflux, frictionvel, temperature,
            solarabs, soilstate, waterfluxbulk, waterstatebulk, waterdiagbulk,
            photosyns, patch_data, col_data, filterp, fn,
            rah, raw, rb, rstem, sa_leaf, sa_stem, sa_internal,
            frac_rad_abs_by_stem, air, bir, cir, cp_leaf, tl_ini, tlbef, zldis,
            temp1, temp2, ur, efeb, wtg, wtl0, wta0, wtstem0, wtga, wtal,
            lw_stem, lw_leaf, wtgq, wtlq0, wtaq0, wtalq, efe, dt_veg, del_arr,
            err_arr, qsatl, el, qsatldT, dth, dqh, delq, obuold, nmozsgn,
            forc_q_col, forc_th_col, forc_pbot_col, forc_rho_col,
            soilevap_beta, soil_resis_sl14, use_lch4, use_hydrstress,
            nlevsno, dtime, params)

        # --- Test for convergence ---
        itlef += 1
        if itlef > ITMIN_CANOPY
            for fi in 1:fn
                p = filterp[fi]
                dele[p] = abs(efe[p] - efeb[p])
                efeb[p] = efe[p]
                det_arr[p] = max(del_arr[p], del2[p])
                frictionvel.num_iter_patch[p] = Float64(itlef)
            end
            fnold = fn
            fn = 0
            for fi in 1:fnold
                p = filterp[fi]
                if !(det_arr[p] < DTMIN_CANOPY && dele[p] < DLEMIN_CANOPY)
                    fn += 1
                    filterp[fn] = p
                end
            end
        end
    end  # End stability iteration

    # =========================================================================
    # Phase 3: Post-iteration diagnostics
    # =========================================================================

    # Restore original filter
    fn = fnorig
    filterp[1:fn] .= fporig[1:fn]

    for fi in 1:fn
        p = filterp[fi]
        c = patch_data.column[p]
        g = patch_data.gridcell[p]

        # Energy balance check
        snl_c = col_data.snl[c]
        lw_grnd = (waterdiagbulk.frac_sno_eff_col[c] * temperature.t_soisno_col[c, snl_c + 1 + nlevsno]^4 +
            (1.0 - waterdiagbulk.frac_sno_eff_col[c] - waterdiagbulk.frac_h2osfc_col[c]) *
                temperature.t_soisno_col[c, 1 + nlevsno]^4 +
            waterdiagbulk.frac_h2osfc_col[c] * temperature.t_h2osfc_col[c]^4)

        err_arr[p] = (1.0 - frac_rad_abs_by_stem[p]) * (solarabs.sabv_patch[p] + air[p] +
            bir[p] * tlbef[p]^3 * (tlbef[p] + 4.0 * dt_veg[p]) + cir[p] * lw_grnd) -
            lw_leaf[p] + lw_stem[p] - energyflux.eflx_sh_veg_patch[p] -
            HVAP * waterfluxbulk.wf.qflx_evap_veg_patch[p] -
            ((temperature.t_veg_patch[p] - tl_ini[p]) * cp_leaf[p] / dtime)

        # Update stem temperature
        if ctrl.use_biomass_heat_storage
            if canopystate.stem_biomass_patch[p] > 0.0
                dt_stem[p] = (frac_rad_abs_by_stem[p] * (solarabs.sabv_patch[p] + air[p] +
                    bir[p] * ts_ini[p]^4 + cir[p] * lw_grnd) -
                    energyflux.eflx_sh_stem_patch[p] +
                    lw_leaf[p] - lw_stem[p]) /
                    (cp_stem[p] / dtime - frac_rad_abs_by_stem[p] * bir[p] * 4.0 * ts_ini[p]^3)
            else
                dt_stem[p] = 0.0
            end

            energyflux.dhsdt_canopy_patch[p] = dt_stem[p] * cp_stem[p] / dtime +
                (temperature.t_veg_patch[p] - tl_ini[p]) * cp_leaf[p] / dtime
            temperature.t_stem_patch[p] = temperature.t_stem_patch[p] + dt_stem[p]
        else
            dt_stem[p] = 0.0
        end

        delt_val = wtal[p] * temperature.t_grnd_col[c] - wtl0[p] * temperature.t_veg_patch[p] -
            wta0[p] * temperature.thm_patch[p] - wtstem0[p] * temperature.t_stem_patch[p]

        # Ground fluxes
        energyflux.taux_patch[p] = -forc_rho_col[c] * forc_u_grc[g] / frictionvel.ram1_patch[p]
        energyflux.tauy_patch[p] = -forc_rho_col[c] * forc_v_grc[g] / frictionvel.ram1_patch[p]
        energyflux.eflx_sh_grnd_patch[p] = CPAIR * forc_rho_col[c] * wtg[p] * delt_val

        # Individual sensible heat fluxes
        delt_snow = wtal[p] * temperature.t_soisno_col[c, snl_c + 1 + nlevsno] -
            wtl0[p] * temperature.t_veg_patch[p] - wta0[p] * temperature.thm_patch[p] -
            wtstem0[p] * temperature.t_stem_patch[p]
        delt_soil = wtal[p] * temperature.t_soisno_col[c, 1 + nlevsno] -
            wtl0[p] * temperature.t_veg_patch[p] - wta0[p] * temperature.thm_patch[p] -
            wtstem0[p] * temperature.t_stem_patch[p]
        delt_h2osfc = wtal[p] * temperature.t_h2osfc_col[c] -
            wtl0[p] * temperature.t_veg_patch[p] - wta0[p] * temperature.thm_patch[p] -
            wtstem0[p] * temperature.t_stem_patch[p]

        energyflux.eflx_sh_snow_patch[p] = CPAIR * forc_rho_col[c] * wtg[p] * delt_snow
        energyflux.eflx_sh_soil_patch[p] = CPAIR * forc_rho_col[c] * wtg[p] * delt_soil
        energyflux.eflx_sh_h2osfc_patch[p] = CPAIR * forc_rho_col[c] * wtg[p] * delt_h2osfc
        waterfluxbulk.wf.qflx_evap_soi_patch[p] = forc_rho_col[c] * wtgq[p] * delq[p]

        # Individual latent heat fluxes
        delq_snow = wtalq[p] * waterdiagbulk.qg_snow_col[c] - wtlq0[p] * qsatl[p] - wtaq0[p] * forc_q_col[c]
        waterfluxbulk.qflx_ev_snow_patch[p] = forc_rho_col[c] * wtgq[p] * delq_snow

        delq_soil = wtalq[p] * waterdiagbulk.qg_soil_col[c] - wtlq0[p] * qsatl[p] - wtaq0[p] * forc_q_col[c]
        waterfluxbulk.qflx_ev_soil_patch[p] = forc_rho_col[c] * wtgq[p] * delq_soil

        delq_h2osfc = wtalq[p] * waterdiagbulk.qg_h2osfc_col[c] - wtlq0[p] * qsatl[p] - wtaq0[p] * forc_q_col[c]
        waterfluxbulk.qflx_ev_h2osfc_patch[p] = forc_rho_col[c] * wtgq[p] * delq_h2osfc

        # 2m reference height temperature
        temperature.t_ref2m_patch[p] = temperature.thm_patch[p] +
            temp1[p] * dth[p] * (1.0 / temp12m[p] - 1.0 / temp1[p])
        temperature.t_ref2m_r_patch[p] = temperature.t_ref2m_patch[p]

        # 2m specific humidity
        waterdiagbulk.q_ref2m_patch[p] = forc_q_col[c] +
            temp2[p] * dqh[p] * (1.0 / temp22m[p] - 1.0 / temp2[p])

        # 2m relative humidity
        (qsat_ref2m, e_ref2m, _, _) = qsat(temperature.t_ref2m_patch[p], forc_pbot_col[c])
        waterdiagbulk.rh_ref2m_patch[p] = smooth_min(100.0,
            waterdiagbulk.q_ref2m_patch[p] / qsat_ref2m * 100.0)
        waterdiagbulk.rh_ref2m_r_patch[p] = waterdiagbulk.rh_ref2m_patch[p]

        # 2m vapor pressure deficit
        waterdiagbulk.vpd_ref2m_patch[p] = e_ref2m * (1.0 - waterdiagbulk.rh_ref2m_patch[p] / 100.0)

        # Human heat stress indices (stub — HumanIndexMod not yet ported)
        # Would be called here in full implementation

        # Downward longwave below canopy
        emv_p = temperature.emv_patch[p]
        emg_c = temperature.emg_col[c]
        energyflux.dlrad_patch[p] = (1.0 - emv_p) * emg_c * forc_lwrad_col[c] +
            emv_p * emg_c * SB * tlbef[p]^3 * (tlbef[p] + 4.0 * dt_veg[p]) *
            (1.0 - frac_rad_abs_by_stem[p]) +
            emv_p * emg_c * SB * ts_ini[p]^3 * (ts_ini[p] + 4.0 * dt_stem[p]) *
            frac_rad_abs_by_stem[p]

        # Upward longwave above canopy
        energyflux.ulrad_patch[p] = ((1.0 - emg_c) * (1.0 - emv_p) * (1.0 - emv_p) * forc_lwrad_col[c] +
            emv_p * (1.0 + (1.0 - emg_c) * (1.0 - emv_p)) * SB *
            tlbef[p]^3 * (tlbef[p] + 4.0 * dt_veg[p]) * (1.0 - frac_rad_abs_by_stem[p]) +
            emv_p * (1.0 + (1.0 - emg_c) * (1.0 - emv_p)) * SB *
            ts_ini[p]^3 * (ts_ini[p] + 4.0 * dt_stem[p]) * frac_rad_abs_by_stem[p] +
            emg_c * (1.0 - emv_p) * SB * lw_grnd)

        # Skin temperature
        temperature.t_skin_patch[p] = emv_p * temperature.t_veg_patch[p] +
            (1.0 - emv_p) * sqrt(sqrt(lw_grnd))

        # Derivative of soil energy flux
        energyflux.cgrnds_patch[p] = energyflux.cgrnds_patch[p] +
            CPAIR * forc_rho_col[c] * wtg[p] * wtal[p]
        energyflux.cgrndl_patch[p] = energyflux.cgrndl_patch[p] +
            forc_rho_col[c] * wtgq[p] * wtalq[p] * waterdiagbulk.dqgdT_col[c]
        energyflux.cgrnd_patch[p] = energyflux.cgrnds_patch[p] +
            energyflux.cgrndl_patch[p] * energyflux.htvp_col[c]

        # Save baseline snocan
        snocan_baseline[p] = waterstatebulk.ws.snocan_patch[p]

        # Update dew accumulation
        t_veg_p = temperature.t_veg_patch[p]
        qflx_evap_veg_p = waterfluxbulk.wf.qflx_evap_veg_patch[p]
        qflx_tran_veg_p = waterfluxbulk.wf.qflx_tran_veg_patch[p]

        # Smooth freeze/thaw partitioning for AD
        _frac_liq = smooth_heaviside(t_veg_p - TFRZ; k=200.0)
        _frac_ice = one(t_veg_p) - _frac_liq
        _net_evap = (qflx_tran_veg_p - qflx_evap_veg_p) * dtime
        waterstatebulk.ws.liqcan_patch[p] = smooth_max(zero(t_veg_p),
            waterstatebulk.ws.liqcan_patch[p] + _frac_liq * _net_evap)
        waterstatebulk.ws.snocan_patch[p] = smooth_max(zero(t_veg_p),
            waterstatebulk.ws.snocan_patch[p] + _frac_ice * _net_evap)
    end

    # --- Post-photosynthesis: PhotosynthesisTotal ---
    photosynthesis_total!(photosyns,
        canopystate.laisun_patch, canopystate.laisha_patch,
        mask_exposedvegp, bounds_patch)

    # Water use efficiency (iwue) — stub for local noon check
    # LUNA and ozone updates would go here in full implementation

    # Filter out patches with small energy balance errors, report large ones
    fnold = fn
    fn = 0
    for fi in 1:fnold
        p = filterp[fi]
        if abs(err_arr[p]) > 0.1
            fn += 1
            filterp[fn] = p
        end
    end
    for fi in 1:fn
        p = filterp[fi]
        @warn "canopy_fluxes!: energy balance error at patch $p, err=$(err_arr[p])"
    end

    return nothing
end
