# ==========================================================================
# BareGroundFluxes — ported from BareGroundFluxesMod.F90
#
# Compute sensible and latent fluxes and their derivatives with respect
# to ground temperature using ground temperatures from previous time step.
# ==========================================================================

# ---------------------------------------------------------------------------
# Module-level parameters (read from parameter file in Fortran)
# ---------------------------------------------------------------------------

Base.@kwdef mutable struct BareGroundFluxesParamsData
    a_coef   ::Float64 = 0.0   # Drag coefficient under less dense canopy (unitless)
    a_exp    ::Float64 = 0.0   # Drag exponent under less dense canopy (unitless)
    wind_min ::Float64 = 0.0   # Minimum wind speed at atmospheric forcing height (m/s)
end

const bareground_fluxes_params = BareGroundFluxesParamsData()

"""
    bareground_fluxes_read_params!(; a_coef, a_exp, wind_min)

Set BareGroundFluxes module parameters (replaces Fortran readParams).
"""
function bareground_fluxes_read_params!(;
        a_coef   ::Real = bareground_fluxes_params.a_coef,
        a_exp    ::Real = bareground_fluxes_params.a_exp,
        wind_min ::Real = bareground_fluxes_params.wind_min)
    bareground_fluxes_params.a_coef   = a_coef
    bareground_fluxes_params.a_exp    = a_exp
    bareground_fluxes_params.wind_min = wind_min
    return nothing
end

# ---------------------------------------------------------------------------
# Phase 1 inner kernel (per-patch simple initialization for the no-exposed-veg
# patches). One thread per filtered patch; each patch writes only its own [p]
# indices plus its own [p, 1:nlevgrnd] row (no reductions). Gather form:
#   fi = @index(Global); p = filterp[fi].
# The nlevgrnd inner zero loop is loop-carried-free, so it stays a sequential
# in-kernel loop (one thread owns the whole column row).
# ---------------------------------------------------------------------------
@kernel function _bgf_phase1_kernel!(
        # written arrays
        btran_patch, t_veg_patch, rssun_patch, rssha_patch, rootr_patch, rresis_patch,
        # read-only arrays
        @Const(filterp), @Const(column), @Const(forc_t_col), @Const(thm_patch),
        @Const(forc_pbot_col),
        # scalars
        nlevgrnd::Int, rgas)

    fi = @index(Global)
    @inbounds begin
        p = filterp[fi]
        c = column[p]
        T = eltype(t_veg_patch)

        btran_patch[p] = zero(T)
        t_veg_patch[p] = forc_t_col[c]
        cf_bare = forc_pbot_col[c] / (rgas * thm_patch[p]) * T(1.0e06)
        rssun_patch[p] = one(T) / T(1.0e15) * cf_bare
        rssha_patch[p] = one(T) / T(1.0e15) * cf_bare
        for j in 1:nlevgrnd
            rootr_patch[p, j]  = zero(T)
            rresis_patch[p, j] = zero(T)
        end
    end
end

"""
    bgf_phase1_update!(energyflux, temperature, photosyns, soilstate, patch_data,
                       filterp, fn, forc_t_col, forc_pbot_col, nlevgrnd)

Launch the bare-ground Phase-1 initialization kernel over the `fn` filtered
patches. Backend-agnostic; one thread per filtered patch.
"""
function bgf_phase1_update!(energyflux, temperature, photosyns, soilstate,
        patch_data, filterp, fn::Int, forc_t_col, forc_pbot_col, nlevgrnd::Int)
    T = eltype(temperature.t_veg_patch)
    _launch!(_bgf_phase1_kernel!, energyflux.btran_patch, temperature.t_veg_patch,
        photosyns.rssun_patch, photosyns.rssha_patch, soilstate.rootr_patch,
        energyflux.rresis_patch, filterp, patch_data.column, forc_t_col,
        temperature.thm_patch, forc_pbot_col, nlevgrnd, T(RGAS); ndrange = fn)
    return nothing
end

# ---------------------------------------------------------------------------
# Phase 2 inner kernel (per-patch initial state + Monin-Obukhov initialization).
# One thread per filtered patch; each patch writes only its own [p] indices.
# monin_obuk_ini is inlined with eltype-generic literals (the host version uses
# Float64 literals and the GRAV module global, neither Metal-safe). grav and
# wind_min are passed as scalars. Arrays are grouped into device-view
# @adapt_structure structs (Metal caps ~31 array resources per kernel); struct
# fields are aliased to locals so the body reads unchanged.
# ---------------------------------------------------------------------------
Base.@kwdef struct BgfP2Out{V}     # written per-patch float outputs
    displa_patch::V; z0mv_patch::V; z0hv_patch::V; z0qv_patch::V; dlrad_patch::V
    ulrad_patch::V; dhsdt_canopy_patch::V; eflx_sh_stem_patch::V; ur::V; dth::V; dqh::V
    z0mg_patch::V; z0hg_patch::V; z0qg_patch::V; um::V; obu::V; num_iter_patch::V
end
Base.@kwdef struct BgfP2In{V}      # read-only float inputs
    forc_u_grc::V; forc_v_grc::V; thm_patch::V; t_grnd_col::V; forc_q_col::V
    qg_col::V; forc_th_col::V; thv_col::V; forc_hgt_u_patch::V; z0mg_col::V
    z0hg_col::V; z0qg_col::V
end
Adapt.@adapt_structure BgfP2Out
Adapt.@adapt_structure BgfP2In

@kernel function _bgf_phase2_kernel!(o, in,
        @Const(filterp), @Const(column), @Const(gridcell),
        # scalars
        wind_min, zetamaxstable, grav)

    fi = @index(Global)
    @inbounds begin
        p = filterp[fi]
        c = column[p]
        g = gridcell[p]
        T = eltype(o.um)

        # Alias grouped struct fields to locals (zero-cost array bindings).
        displa_patch = o.displa_patch; z0mv_patch = o.z0mv_patch; z0hv_patch = o.z0hv_patch
        z0qv_patch = o.z0qv_patch; dlrad_patch = o.dlrad_patch; ulrad_patch = o.ulrad_patch
        dhsdt_canopy_patch = o.dhsdt_canopy_patch; eflx_sh_stem_patch = o.eflx_sh_stem_patch
        ur = o.ur; dth = o.dth; dqh = o.dqh; z0mg_patch = o.z0mg_patch
        z0hg_patch = o.z0hg_patch; z0qg_patch = o.z0qg_patch; um = o.um; obu = o.obu
        num_iter_patch = o.num_iter_patch
        forc_u_grc = in.forc_u_grc; forc_v_grc = in.forc_v_grc; thm_patch = in.thm_patch
        t_grnd_col = in.t_grnd_col; forc_q_col = in.forc_q_col; qg_col = in.qg_col
        forc_th_col = in.forc_th_col; thv_col = in.thv_col
        forc_hgt_u_patch = in.forc_hgt_u_patch; z0mg_col = in.z0mg_col
        z0hg_col = in.z0hg_col; z0qg_col = in.z0qg_col

        # Initialization variables
        displa_patch[p] = zero(T)
        z0mv_patch[p]   = zero(T)
        z0hv_patch[p]   = zero(T)
        z0qv_patch[p]   = zero(T)
        dlrad_patch[p]  = zero(T)
        ulrad_patch[p]  = zero(T)
        dhsdt_canopy_patch[p] = zero(T)
        eflx_sh_stem_patch[p] = zero(T)

        ur[p]  = smooth_max(wind_min, sqrt(forc_u_grc[g]^2 + forc_v_grc[g]^2))
        dth[p] = thm_patch[p] - t_grnd_col[c]
        dqh[p] = forc_q_col[c] - qg_col[c]
        dthv   = dth[p] * (one(T) + T(0.61) * forc_q_col[c]) + T(0.61) * forc_th_col[c] * dqh[p]
        zldis  = forc_hgt_u_patch[p]

        # Copy column roughness to local patch-level arrays
        z0mg_patch[p] = z0mg_col[c]
        z0hg_patch[p] = z0hg_col[c]
        z0qg_patch[p] = z0qg_col[c]

        # Initialize Monin-Obukhov length and wind speed (inlined monin_obuk_ini)
        wc = T(0.5)
        um_val = dthv >= zero(T) ? max(ur[p], T(0.1)) : sqrt(ur[p] * ur[p] + wc * wc)
        rib = grav * zldis * dthv / (thv_col[c] * um_val * um_val)
        if rib >= zero(T)  # neutral or stable
            zeta = rib * log(zldis / z0mg_patch[p]) / (one(T) - T(5.0) * min(rib, T(0.19)))
            zeta = min(zetamaxstable, max(zeta, T(0.01)))
        else               # unstable
            zeta = rib * log(zldis / z0mg_patch[p])
            zeta = max(T(-100.0), min(zeta, T(-0.01)))
        end
        um[p]  = um_val
        obu[p] = zldis / zeta

        # Initialize iteration counter
        num_iter_patch[p] = zero(eltype(num_iter_patch))
    end
end

"""
    bgf_phase2_update!(canopystate, frictionvel, energyflux, temperature,
                       waterdiagbulk, patch_data, filterp, fn, ur, dth, dqh, um, obu,
                       forc_u_grc, forc_v_grc, forc_q_col, forc_th_col, params)

Launch the bare-ground Phase-2 initial-state / Monin-Obukhov-init kernel over the
`fn` filtered patches. Backend-agnostic; one thread per filtered patch.
"""
function bgf_phase2_update!(canopystate, frictionvel, energyflux, temperature,
        waterdiagbulk, patch_data, filterp, fn::Int, ur, dth, dqh, um, obu,
        forc_u_grc, forc_v_grc, forc_q_col, forc_th_col, params)
    o = BgfP2Out(; displa_patch = canopystate.displa_patch, z0mv_patch = frictionvel.z0mv_patch,
        z0hv_patch = frictionvel.z0hv_patch, z0qv_patch = frictionvel.z0qv_patch,
        dlrad_patch = energyflux.dlrad_patch, ulrad_patch = energyflux.ulrad_patch,
        dhsdt_canopy_patch = energyflux.dhsdt_canopy_patch,
        eflx_sh_stem_patch = energyflux.eflx_sh_stem_patch, ur = ur, dth = dth, dqh = dqh,
        z0mg_patch = frictionvel.z0mg_patch, z0hg_patch = frictionvel.z0hg_patch,
        z0qg_patch = frictionvel.z0qg_patch, um = um, obu = obu,
        num_iter_patch = frictionvel.num_iter_patch)
    in = BgfP2In(; forc_u_grc = forc_u_grc, forc_v_grc = forc_v_grc,
        thm_patch = temperature.thm_patch, t_grnd_col = temperature.t_grnd_col,
        forc_q_col = forc_q_col, qg_col = waterdiagbulk.qg_col, forc_th_col = forc_th_col,
        thv_col = temperature.thv_col, forc_hgt_u_patch = frictionvel.forc_hgt_u_patch,
        z0mg_col = frictionvel.z0mg_col, z0hg_col = frictionvel.z0hg_col,
        z0qg_col = frictionvel.z0qg_col)
    T = eltype(um)
    let be = _kernel_backend(um)
        _bgf_phase2_kernel!(be)(o, in, filterp, patch_data.column, patch_data.gridcell,
            T(params.wind_min), T(frictionvel.zetamaxstable), T(GRAV); ndrange = fn)
        KA.synchronize(be)
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Stability-iteration inner kernel (per-patch stability/roughness update).
#
# Kernelized form of the inner `for fi in 1:fn` loop in Phase 3 of
# bareground_fluxes!. One thread per filtered patch; each patch writes only its
# own [p] indices (no reductions). The String z0param_method branch is resolved
# on the host into z0flag (1 = ZengWang2007, 2 = Meier2022) and branched on here.
# Arrays are grouped into device-view @adapt_structure structs (Metal caps ~31
# array resources per kernel); struct fields are aliased to locals so the body
# reads unchanged.
# ---------------------------------------------------------------------------
Base.@kwdef struct BgfP3Out{V}     # written per-patch float outputs
    z0hg_patch::V; z0qg_patch::V; forc_hgt_u_patch::V; forc_hgt_t_patch::V
    forc_hgt_q_patch::V; zeta_patch::V; um::V; obu::V; num_iter_patch::V
end
Base.@kwdef struct BgfP3In{V}      # read-only float inputs
    temp1::V; temp2::V; dth::V; dqh::V; z0mg_patch::V; displa_patch::V
    forc_hgt_u_grc::V; forc_hgt_t_grc::V; forc_hgt_q_grc::V; forc_q_col::V
    forc_th_col::V; zldis::V; ustar::V; ur::V; thv_col::V; beta_col::V; zii::V
end
Adapt.@adapt_structure BgfP3Out
Adapt.@adapt_structure BgfP3In

@kernel function _bgf_stability_kernel!(o, in,
        @Const(filterp), @Const(column), @Const(gridcell),
        # scalars
        iter::Int, z0flag::Int, a_coef, a_exp, zetamaxstable,
        nu_param, meier_param3, beta_param, vkc, grav)

    fi = @index(Global)
    @inbounds begin
        p = filterp[fi]
        c = column[p]
        g = gridcell[p]
        T = eltype(o.um)

        # Alias grouped struct fields to locals (zero-cost array bindings).
        z0hg_patch = o.z0hg_patch; z0qg_patch = o.z0qg_patch
        forc_hgt_u_patch = o.forc_hgt_u_patch; forc_hgt_t_patch = o.forc_hgt_t_patch
        forc_hgt_q_patch = o.forc_hgt_q_patch; zeta_patch = o.zeta_patch
        um = o.um; obu = o.obu; num_iter_patch = o.num_iter_patch
        temp1 = in.temp1; temp2 = in.temp2; dth = in.dth; dqh = in.dqh
        z0mg_patch = in.z0mg_patch; displa_patch = in.displa_patch
        forc_hgt_u_grc = in.forc_hgt_u_grc; forc_hgt_t_grc = in.forc_hgt_t_grc
        forc_hgt_q_grc = in.forc_hgt_q_grc; forc_q_col = in.forc_q_col
        forc_th_col = in.forc_th_col; zldis = in.zldis; ustar = in.ustar; ur = in.ur
        thv_col = in.thv_col; beta_col = in.beta_col; zii = in.zii

        tstar = temp1[p] * dth[p]
        qstar = temp2[p] * dqh[p]

        if z0flag == 2
            # Meier2022 — after Yang et al. (2008)
            z0hg_patch[p] = meier_param3 * nu_param / ustar[p] *
                exp(-beta_param * ustar[p]^T(0.5) * smooth_abs(tstar)^T(0.25))
        else
            # ZengWang2007
            z0hg_patch[p] = z0mg_patch[p] /
                exp(a_coef * (ustar[p] * z0mg_patch[p] / nu_param)^a_exp)
        end

        z0qg_patch[p] = z0hg_patch[p]

        # Update the forcing heights for new roughness lengths
        forc_hgt_u_patch[p] = forc_hgt_u_grc[g] + z0mg_patch[p] + displa_patch[p]
        forc_hgt_t_patch[p] = forc_hgt_t_grc[g] + z0hg_patch[p] + displa_patch[p]
        forc_hgt_q_patch[p] = forc_hgt_q_grc[g] + z0qg_patch[p] + displa_patch[p]

        thvstar = tstar * (one(T) + T(0.61) * forc_q_col[c]) + T(0.61) * forc_th_col[c] * qstar
        zeta_patch[p] = zldis[p] * vkc * grav * thvstar / (ustar[p]^2 * thv_col[c])

        if zeta_patch[p] >= zero(T)  # stable
            zeta_patch[p] = smooth_clamp(zeta_patch[p], T(0.01), zetamaxstable)
            um[p] = smooth_max(ur[p], T(0.1))
        else  # unstable
            zeta_patch[p] = smooth_clamp(zeta_patch[p], T(-100.0), T(-0.01))
            wc_arg = smooth_max(-grav * ustar[p] * thvstar * zii[c] / thv_col[c], zero(T))
            wc = beta_col[c] * wc_arg^T(0.333)
            um[p] = sqrt(ur[p]^2 + wc^2)
        end
        obu[p] = zldis[p] / zeta_patch[p]

        num_iter_patch[p] = T(iter)
    end
end

"""
    bgf_stability_update!(frictionvel, canopystate, temperature, col_data,
                          patch_data, filterp, fn, temp1, temp2, dth, dqh,
                          zldis, ustar, ur, um, obu, forc_q_col, forc_th_col,
                          forc_hgt_u_grc, forc_hgt_t_grc, forc_hgt_q_grc,
                          iter, z0flag, params)

Launch the bare-ground stability inner-loop kernel over the `fn` filtered
patches. Backend-agnostic (CPU loop or GPU); one thread per filtered patch.
Replaces the inner `for fi in 1:fn` stability/roughness update in Phase 3.
"""
function bgf_stability_update!(frictionvel, canopystate, temperature, col_data,
        patch_data, filterp, fn::Int,
        temp1, temp2, dth, dqh, zldis, ustar, ur, um, obu,
        forc_q_col, forc_th_col, forc_hgt_u_grc, forc_hgt_t_grc, forc_hgt_q_grc,
        iter::Int, z0flag::Int, params)
    o = BgfP3Out(; z0hg_patch = frictionvel.z0hg_patch, z0qg_patch = frictionvel.z0qg_patch,
        forc_hgt_u_patch = frictionvel.forc_hgt_u_patch,
        forc_hgt_t_patch = frictionvel.forc_hgt_t_patch,
        forc_hgt_q_patch = frictionvel.forc_hgt_q_patch, zeta_patch = frictionvel.zeta_patch,
        um = um, obu = obu, num_iter_patch = frictionvel.num_iter_patch)
    in = BgfP3In(; temp1 = temp1, temp2 = temp2, dth = dth, dqh = dqh,
        z0mg_patch = frictionvel.z0mg_patch, displa_patch = canopystate.displa_patch,
        forc_hgt_u_grc = forc_hgt_u_grc, forc_hgt_t_grc = forc_hgt_t_grc,
        forc_hgt_q_grc = forc_hgt_q_grc, forc_q_col = forc_q_col, forc_th_col = forc_th_col,
        zldis = zldis, ustar = ustar, ur = ur, thv_col = temperature.thv_col,
        beta_col = temperature.beta_col, zii = col_data.zii)
    T = eltype(um)
    let be = _kernel_backend(um)
        _bgf_stability_kernel!(be)(o, in, filterp, patch_data.column, patch_data.gridcell,
            iter, z0flag, T(params.a_coef), T(params.a_exp), T(frictionvel.zetamaxstable),
            T(NU_PARAM), T(MEIER_PARAM3), T(BETA_PARAM), T(VKC), T(GRAV); ndrange = fn)
        KA.synchronize(be)
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Phase 4 post-iteration kernel (per-patch final fluxes + diagnostics).
#
# Kernelized form of the final `for fi in 1:fn` loop. One thread per filtered
# patch; each patch writes only its own [p] indices plus its own [c] history
# write-back of z0hg_col / z0qg_col (no cross-patch reductions — at most one
# patch per column, matching the scalar version).
#
# GPU-hostile host-side items resolved before launch and passed as scalars:
#   - do_soilevap_beta() / do_soil_resistance_sl14() — control-flag functions
#     reading module globals; evaluated once on the host as Bool scalars.
#   - use_lch4 and the methane-array length guard (grnd_ch4_cond_len::Int).
#   - ISTSOIL / ISTCROP landunit-type codes passed as Int scalars.
# qsat(...) is kernel-safe (4-tuple destructuring preserved). The grouped
# device-view structs keep the per-patch field count under the Metal ~31-buffer
# cap; struct fields are aliased to locals so the physics body reads unchanged
# except Float64 literals → T(...).
# ---------------------------------------------------------------------------
# grnd_ch4_cond_patch comes from a Float64 kwarg default (Float64[]) and stays
# Float64 even when the state is Float32/Dual — it gets its own type param (Vg)
# so the struct accepts the mixed eltypes (a single V would force one type).
Base.@kwdef struct BgfP4Out{V,Vg}  # written per-patch float outputs
    ram1_patch::V; ustar_patch::V; cgrnds_patch::V; cgrndl_patch::V; cgrnd_patch::V
    taux_patch::V; tauy_patch::V; eflx_sh_grnd_patch::V; eflx_sh_tot_patch::V
    eflx_sh_snow_patch::V; eflx_sh_soil_patch::V; eflx_sh_h2osfc_patch::V
    qflx_tran_veg_patch::V; qflx_evap_veg_patch::V; qflx_evap_soi_patch::V
    qflx_evap_tot_patch::V; qflx_ev_snow_patch::V; qflx_ev_soil_patch::V
    qflx_ev_h2osfc_patch::V; t_ref2m_patch::V; q_ref2m_patch::V; rh_ref2m_patch::V
    rh_ref2m_r_patch::V; t_ref2m_r_patch::V; kbm1_patch::V; grnd_ch4_cond_patch::Vg
end
Base.@kwdef struct BgfP4P{V}       # read-only per-patch floats
    ustar::V; um::V; temp1::V; temp2::V; temp12m::V; temp22m::V; dth::V; dqh::V
    thm_patch::V; z0mg_patch::V; z0hg_patch::V; z0qg_patch::V
end
Base.@kwdef struct BgfP4C{V}       # read-only per-column floats
    forc_rho_col::V; soilbeta_col::V; soilresis_col::V; t_h2osfc_col::V
    qg_snow_col::V; qg_soil_col::V; qg_h2osfc_col::V; dqgdT_col::V; htvp_col::V
    z0hg_col::V; z0qg_col::V
end
Base.@kwdef struct BgfP4M{M}       # read-only 2D arrays
    h2osoi_liq_col::M; h2osoi_ice_col::M; dz::M; watsat_col::M; t_soisno_col::M
end
Adapt.@adapt_structure BgfP4Out
Adapt.@adapt_structure BgfP4P
Adapt.@adapt_structure BgfP4C
Adapt.@adapt_structure BgfP4M

# Scalar constants bundled into one isbits arg (Metal caps total kernel args ~31).
struct BgfP4Scalars{S}
    cpair::S; denh2o::S; denice::S
end

@kernel function _bgf_phase4_kernel!(o, pp, cc, mm, sc,
        @Const(filterp), @Const(column), @Const(gridcell), @Const(landunit),
        @Const(itype_lun), @Const(snl), @Const(forc_q_col),
        @Const(forc_u_grc), @Const(forc_v_grc), @Const(forc_pbot_col),
        # scalars
        soilevap_beta::Bool, soil_resis_sl14::Bool, use_lch4::Bool,
        grnd_ch4_cond_len::Int, istsoil::Int, istcrop::Int, nlevsno::Int)

    fi = @index(Global)
    @inbounds begin
        p = filterp[fi]
        c = column[p]
        g = gridcell[p]
        l = landunit[p]
        T = eltype(o.ram1_patch)

        # Alias grouped struct fields to locals (zero-cost array bindings).
        ram1_patch = o.ram1_patch; ustar_patch = o.ustar_patch
        cgrnds_patch = o.cgrnds_patch; cgrndl_patch = o.cgrndl_patch; cgrnd_patch = o.cgrnd_patch
        taux_patch = o.taux_patch; tauy_patch = o.tauy_patch
        eflx_sh_grnd_patch = o.eflx_sh_grnd_patch; eflx_sh_tot_patch = o.eflx_sh_tot_patch
        eflx_sh_snow_patch = o.eflx_sh_snow_patch; eflx_sh_soil_patch = o.eflx_sh_soil_patch
        eflx_sh_h2osfc_patch = o.eflx_sh_h2osfc_patch
        qflx_tran_veg_patch = o.qflx_tran_veg_patch; qflx_evap_veg_patch = o.qflx_evap_veg_patch
        qflx_evap_soi_patch = o.qflx_evap_soi_patch; qflx_evap_tot_patch = o.qflx_evap_tot_patch
        qflx_ev_snow_patch = o.qflx_ev_snow_patch; qflx_ev_soil_patch = o.qflx_ev_soil_patch
        qflx_ev_h2osfc_patch = o.qflx_ev_h2osfc_patch; t_ref2m_patch = o.t_ref2m_patch
        q_ref2m_patch = o.q_ref2m_patch; rh_ref2m_patch = o.rh_ref2m_patch
        rh_ref2m_r_patch = o.rh_ref2m_r_patch; t_ref2m_r_patch = o.t_ref2m_r_patch
        kbm1_patch = o.kbm1_patch; grnd_ch4_cond_patch = o.grnd_ch4_cond_patch
        ustar = pp.ustar; um = pp.um; temp1 = pp.temp1; temp2 = pp.temp2
        temp12m = pp.temp12m; temp22m = pp.temp22m; dth = pp.dth; dqh = pp.dqh
        thm_patch = pp.thm_patch; z0mg_patch = pp.z0mg_patch; z0hg_patch = pp.z0hg_patch
        z0qg_patch = pp.z0qg_patch
        forc_rho_col = cc.forc_rho_col; soilbeta_col = cc.soilbeta_col
        soilresis_col = cc.soilresis_col; t_h2osfc_col = cc.t_h2osfc_col
        qg_snow_col = cc.qg_snow_col; qg_soil_col = cc.qg_soil_col
        qg_h2osfc_col = cc.qg_h2osfc_col; dqgdT_col = cc.dqgdT_col; htvp_col = cc.htvp_col
        z0hg_col = cc.z0hg_col; z0qg_col = cc.z0qg_col
        h2osoi_liq_col = mm.h2osoi_liq_col; h2osoi_ice_col = mm.h2osoi_ice_col
        dz = mm.dz; watsat_col = mm.watsat_col; t_soisno_col = mm.t_soisno_col
        cpair = sc.cpair; denh2o = sc.denh2o; denice = sc.denice

        # Determine aerodynamic resistances
        ram  = one(T) / (ustar[p]^2 / um[p])
        rah  = one(T) / (temp1[p] * ustar[p])
        raw  = one(T) / (temp2[p] * ustar[p])
        raih = forc_rho_col[c] * cpair / rah

        if use_lch4 && grnd_ch4_cond_len >= p
            grnd_ch4_cond_patch[p] = one(T) / raw
        end

        # Soil evaporation resistance
        # Fortran index 1 = first soil layer = Julia index 1 + nlevsno
        www = (h2osoi_liq_col[c, 1 + nlevsno] / denh2o +
               h2osoi_ice_col[c, 1 + nlevsno] / denice) /
              dz[c, 1 + nlevsno] / watsat_col[c, 1]
        www = smooth_clamp(www, zero(T), one(T))

        # Soil evaporation: beta or resistance method
        if dqh[p] > zero(T)  # dew (beta not applied)
            raiw = forc_rho_col[c] / raw
        else
            raiw = zero(T)
            if soilevap_beta
                # Lee and Pielke 1992 beta
                raiw = soilbeta_col[c] * forc_rho_col[c] / raw
            end
            if soil_resis_sl14
                # Swenson & Lawrence 2014 soil resistance
                raiw = forc_rho_col[c] / (raw + soilresis_col[c])
            end
        end

        ram1_patch[p] = ram  # pass value to global variable
        ustar_patch[p] = ustar[p]

        # Derivative of fluxes with respect to ground temperature
        cgrnds_patch[p] = raih
        cgrndl_patch[p] = raiw * dqgdT_col[c]
        cgrnd_patch[p]  = cgrnds_patch[p] + htvp_col[c] * cgrndl_patch[p]

        # Surface fluxes of momentum, sensible and latent heat
        taux_patch[p]          = -forc_rho_col[c] * forc_u_grc[g] / ram
        tauy_patch[p]          = -forc_rho_col[c] * forc_v_grc[g] / ram
        eflx_sh_grnd_patch[p]  = -raih * dth[p]
        eflx_sh_tot_patch[p]   = eflx_sh_grnd_patch[p]

        # Compute sensible heat fluxes individually
        eflx_sh_snow_patch[p]   = -raih * (thm_patch[p] -
            t_soisno_col[c, snl[c] + 1 + nlevsno])
        eflx_sh_soil_patch[p]   = -raih * (thm_patch[p] - t_soisno_col[c, 1 + nlevsno])
        eflx_sh_h2osfc_patch[p] = -raih * (thm_patch[p] - t_h2osfc_col[c])

        # Water fluxes from soil
        qflx_tran_veg_patch[p] = zero(T)
        qflx_evap_veg_patch[p] = zero(T)
        qflx_evap_soi_patch[p] = -raiw * dqh[p]
        qflx_evap_tot_patch[p] = qflx_evap_soi_patch[p]

        # Compute latent heat fluxes individually
        qflx_ev_snow_patch[p]   = -raiw * (forc_q_col[c] - qg_snow_col[c])
        qflx_ev_soil_patch[p]   = -raiw * (forc_q_col[c] - qg_soil_col[c])
        qflx_ev_h2osfc_patch[p] = -raiw * (forc_q_col[c] - qg_h2osfc_col[c])

        # 2 m height air temperature
        t_ref2m_patch[p] = thm_patch[p] +
            temp1[p] * dth[p] * (one(T) / temp12m[p] - one(T) / temp1[p])

        # 2 m height specific humidity
        q_ref2m_patch[p] = forc_q_col[c] +
            temp2[p] * dqh[p] * (one(T) / temp22m[p] - one(T) / temp2[p])

        # 2 m height relative humidity
        (qsat_ref2m, e_ref2m, _, _) = qsat(t_ref2m_patch[p], forc_pbot_col[c])
        rh_ref2m_patch[p] = smooth_min(T(100.0), q_ref2m_patch[p] / qsat_ref2m * T(100.0))

        if itype_lun[l] == istsoil || itype_lun[l] == istcrop
            rh_ref2m_r_patch[p] = rh_ref2m_patch[p]
            t_ref2m_r_patch[p] = t_ref2m_patch[p]
        end

        kbm1_patch[p] = log(z0mg_patch[p] / z0hg_patch[p])

        # Copy local patch ground roughness back to column arrays for history output
        z0hg_col[c] = z0hg_patch[p]
        z0qg_col[c] = z0qg_patch[p]
    end
end

"""
    bgf_phase4_update!(...)

Launch the bare-ground Phase-4 post-iteration kernel (final fluxes +
diagnostics) over the `fn` filtered patches. Backend-agnostic; one thread per
filtered patch. The two control-flag functions, the methane-array guard, and
the ISTSOIL/ISTCROP landunit codes are resolved on the host and passed as
scalars.
"""
function bgf_phase4_update!(energyflux, frictionvel, temperature, soilstate,
        waterfluxbulk, waterstatebulk, waterdiagbulk, patch_data, col_data, lun_data,
        filterp, fn::Int, ustar, um, temp1, temp2, temp12m, temp22m, dth, dqh,
        forc_q_col, forc_rho_col, forc_u_grc, forc_v_grc, forc_pbot_col,
        grnd_ch4_cond_patch, soilevap_beta::Bool, soil_resis_sl14::Bool,
        use_lch4::Bool, nlevsno::Int)
    # grnd_ch4_cond_patch defaults to a host Float64[] when use_lch4=false; on a
    # device backend it must be device-resident or the BgfP4Out kernel arg is
    # non-isbits (Vector{Float64} field). Rebuild it on the working backend.
    if !(frictionvel.ram1_patch isa Array) && grnd_ch4_cond_patch isa Array
        _Tg = eltype(frictionvel.ram1_patch)
        grnd_ch4_cond_patch = fill!(
            similar(frictionvel.ram1_patch, _Tg, length(grnd_ch4_cond_patch)), zero(_Tg))
    end
    o = BgfP4Out(; ram1_patch = frictionvel.ram1_patch, ustar_patch = frictionvel.ustar_patch,
        cgrnds_patch = energyflux.cgrnds_patch, cgrndl_patch = energyflux.cgrndl_patch,
        cgrnd_patch = energyflux.cgrnd_patch, taux_patch = energyflux.taux_patch,
        tauy_patch = energyflux.tauy_patch, eflx_sh_grnd_patch = energyflux.eflx_sh_grnd_patch,
        eflx_sh_tot_patch = energyflux.eflx_sh_tot_patch,
        eflx_sh_snow_patch = energyflux.eflx_sh_snow_patch,
        eflx_sh_soil_patch = energyflux.eflx_sh_soil_patch,
        eflx_sh_h2osfc_patch = energyflux.eflx_sh_h2osfc_patch,
        qflx_tran_veg_patch = waterfluxbulk.wf.qflx_tran_veg_patch,
        qflx_evap_veg_patch = waterfluxbulk.wf.qflx_evap_veg_patch,
        qflx_evap_soi_patch = waterfluxbulk.wf.qflx_evap_soi_patch,
        qflx_evap_tot_patch = waterfluxbulk.wf.qflx_evap_tot_patch,
        qflx_ev_snow_patch = waterfluxbulk.qflx_ev_snow_patch,
        qflx_ev_soil_patch = waterfluxbulk.qflx_ev_soil_patch,
        qflx_ev_h2osfc_patch = waterfluxbulk.qflx_ev_h2osfc_patch,
        t_ref2m_patch = temperature.t_ref2m_patch, q_ref2m_patch = waterdiagbulk.q_ref2m_patch,
        rh_ref2m_patch = waterdiagbulk.rh_ref2m_patch,
        rh_ref2m_r_patch = waterdiagbulk.rh_ref2m_r_patch,
        t_ref2m_r_patch = temperature.t_ref2m_r_patch, kbm1_patch = frictionvel.kbm1_patch,
        grnd_ch4_cond_patch = grnd_ch4_cond_patch)
    pp = BgfP4P(; ustar = ustar, um = um, temp1 = temp1, temp2 = temp2,
        temp12m = temp12m, temp22m = temp22m, dth = dth, dqh = dqh,
        thm_patch = temperature.thm_patch, z0mg_patch = frictionvel.z0mg_patch,
        z0hg_patch = frictionvel.z0hg_patch, z0qg_patch = frictionvel.z0qg_patch)
    cc = BgfP4C(; forc_rho_col = forc_rho_col, soilbeta_col = soilstate.soilbeta_col,
        soilresis_col = soilstate.soilresis_col, t_h2osfc_col = temperature.t_h2osfc_col,
        qg_snow_col = waterdiagbulk.qg_snow_col, qg_soil_col = waterdiagbulk.qg_soil_col,
        qg_h2osfc_col = waterdiagbulk.qg_h2osfc_col, dqgdT_col = waterdiagbulk.dqgdT_col,
        htvp_col = energyflux.htvp_col, z0hg_col = frictionvel.z0hg_col,
        z0qg_col = frictionvel.z0qg_col)
    mm = BgfP4M(; h2osoi_liq_col = waterstatebulk.ws.h2osoi_liq_col,
        h2osoi_ice_col = waterstatebulk.ws.h2osoi_ice_col, dz = col_data.dz,
        watsat_col = soilstate.watsat_col, t_soisno_col = temperature.t_soisno_col)
    T = eltype(frictionvel.ram1_patch)
    sc = BgfP4Scalars{T}(T(CPAIR), T(DENH2O), T(DENICE))
    let be = _kernel_backend(frictionvel.ram1_patch)
        _bgf_phase4_kernel!(be)(o, pp, cc, mm, sc, filterp, patch_data.column,
            patch_data.gridcell, patch_data.landunit, lun_data.itype, col_data.snl,
            forc_q_col, forc_u_grc, forc_v_grc, forc_pbot_col,
            soilevap_beta, soil_resis_sl14, use_lch4, length(grnd_ch4_cond_patch),
            ISTSOIL, ISTCROP, nlevsno; ndrange = fn)
        KA.synchronize(be)
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Main subroutine
# ---------------------------------------------------------------------------

"""
    bareground_fluxes!(...)

Compute sensible and latent fluxes and their derivatives with respect to
ground temperature using ground temperatures from previous time step.

Ported from: BareGroundFluxesMod.F90 :: BareGroundFluxes
"""
function bareground_fluxes!(
        # Data structures
        canopystate      ::CanopyStateData,
        energyflux       ::EnergyFluxData,
        frictionvel      ::FrictionVelocityData,
        temperature      ::TemperatureData,
        soilstate        ::SoilStateData,
        waterfluxbulk    ::WaterFluxBulkData,
        waterstatebulk   ::WaterStateBulkData,
        waterdiagbulk    ::WaterDiagnosticBulkData,
        photosyns        ::PhotosynthesisData,
        patch_data       ::PatchData,
        col_data         ::ColumnData,
        lun_data         ::LandunitData,
        mask_noexposedvegp ::AbstractVector{Bool},
        bounds_patch     ::UnitRange{Int},
        # Atmospheric forcing (column-level)
        forc_q_col       ::AbstractVector{<:Real},
        forc_pbot_col    ::AbstractVector{<:Real},
        forc_th_col      ::AbstractVector{<:Real},
        forc_rho_col     ::AbstractVector{<:Real},
        forc_t_col       ::AbstractVector{<:Real},
        # Atmospheric forcing (gridcell-level)
        forc_u_grc       ::AbstractVector{<:Real},
        forc_v_grc       ::AbstractVector{<:Real},
        forc_hgt_t_grc   ::AbstractVector{<:Real},
        forc_hgt_u_grc   ::AbstractVector{<:Real},
        forc_hgt_q_grc   ::AbstractVector{<:Real};
        # Feature flags
        use_lch4         ::Bool   = false,
        z0param_method   ::String = "ZengWang2007",
        # CH4 conductance output (optional)
        grnd_ch4_cond_patch ::AbstractVector{<:Real} = Float64[])

    np   = length(bounds_patch)
    begp = first(bounds_patch)
    endp = last(bounds_patch)
    nlevgrnd = varpar.nlevgrnd
    nlevsno  = varpar.nlevsno

    niters = 3  # maximum number of iterations for surface temperature

    params = bareground_fluxes_params

    # --- Build integer filter from the (host-only) mask, then copy to the device
    # backend in bulk (copyto! is a transfer, not device scalar indexing); the
    # kernels only ever read filterp. filterp_host is kept for the few remaining
    # host-side scans. _sc/_sci are the float/Int scratch factories (similar() of
    # a state array + fill!), so all scratch shares the state's backend —
    # device-resident on GPU. A bare zeros() would force a host Array and break a
    # whole-function GPU run.
    FT = eltype(forc_t_col)
    _sc(dims...)  = fill!(similar(temperature.t_veg_patch, FT,  dims...), zero(FT))
    _sci(dims...) = fill!(similar(patch_data.column,       Int, dims...), 0)

    filterp = _sci(np)
    filterp_host = zeros(Int, np)
    fn = 0
    # Host-scan the mask (Array() it once so a device mask isn't scalar-indexed).
    mask_nev_host = mask_noexposedvegp isa Array ? mask_noexposedvegp : Array(mask_noexposedvegp)
    for p in bounds_patch
        if mask_nev_host[p]
            fn += 1
            filterp_host[fn] = p
        end
    end
    copyto!(filterp, filterp_host)

    # --- Local work arrays (patch-indexed, device-resident) ---
    zldis  = _sc(endp)
    dth    = _sc(endp)
    dqh    = _sc(endp)
    obu    = _sc(endp)
    ur     = _sc(endp)
    um     = _sc(endp)
    temp1  = _sc(endp)
    temp12m = _sc(endp)
    temp2  = _sc(endp)
    temp22m = _sc(endp)
    ustar  = _sc(endp)
    fm     = _sc(endp)

    # Per-patch active mask for the friction-velocity profile kernel (all filtered
    # patches active — bare-ground has a fixed-iteration stability loop, no
    # convergence-based deactivation). Device-resident (similar() of a Bool state).
    active = fill!(similar(patch_data.active, Bool, endp), false)
    active_host = fill(false, endp)
    for fi in 1:fn
        active_host[filterp_host[fi]] = true
    end
    copyto!(active, active_host)

    # =========================================================================
    # Phase 1: Initialization — simple settings for no-exposed-veg patches
    # =========================================================================

    # No bare-ground (non-exposed-veg) patches → nothing to do. Skip the per-patch
    # phases; the manual struct-arg kernel launches below would otherwise divide by
    # the empty workgroup (ndrange = fn = 0) on the GPU backend.
    fn == 0 && return nothing

    bgf_phase1_update!(energyflux, temperature, photosyns, soilstate, patch_data,
        filterp, fn, forc_t_col, forc_pbot_col, nlevgrnd)

    # =========================================================================
    # Phase 2: Compute initial state and stability parameters
    # =========================================================================

    # zldis is the per-patch forcing height (= frictionvel.forc_hgt_u_patch[p]);
    # set it before Phase 2 so the kernel reads it (the kernel writes
    # forc_hgt_u_patch only in Phase 3). On device this is a backend-agnostic
    # broadcast copy over the patch range.
    @views zldis[bounds_patch] .= frictionvel.forc_hgt_u_patch[bounds_patch]

    bgf_phase2_update!(canopystate, frictionvel, energyflux, temperature,
        waterdiagbulk, patch_data, filterp, fn, ur, dth, dqh, um, obu,
        forc_u_grc, forc_v_grc, forc_q_col, forc_th_col, params)

    # =========================================================================
    # Phase 3: Stability iteration
    # =========================================================================

    # Resolve the String z0param_method to an Int flag on the host (a kernel
    # cannot compare Strings): 1 = ZengWang2007 (default), 2 = Meier2022.
    z0flag = z0param_method == "Meier2022" ? 2 : 1

    # Device-resident filter slice for the friction-velocity profile kernel.
    filt_arr = filterp[1:fn]

    for iter in 1:niters

        # Friction velocity calculation (already kernelized; pass device active so
        # all filtered patches are processed, matching the scalar all-active path).
        friction_velocity!(frictionvel, fn, filt_arr,
            canopystate.displa_patch, frictionvel.z0mg_patch,
            frictionvel.z0hg_patch, frictionvel.z0qg_patch,
            obu, iter, ur, um, ustar,
            temp1, temp2, temp12m, temp22m, fm; active = active)

        bgf_stability_update!(frictionvel, canopystate, temperature, col_data,
            patch_data, filterp, fn,
            temp1, temp2, dth, dqh, zldis, ustar, ur, um, obu,
            forc_q_col, forc_th_col, forc_hgt_u_grc, forc_hgt_t_grc, forc_hgt_q_grc,
            iter, z0flag, params)
    end  # end stability iteration

    # =========================================================================
    # Phase 4: Post-iteration — compute final fluxes and diagnostics
    # =========================================================================

    # Control-flag functions read module globals — not GPU-safe. Resolve once on
    # the host to Bool scalars and pass into the Phase-4 kernel.
    soilevap_beta   = do_soilevap_beta()
    soil_resis_sl14 = do_soil_resistance_sl14()

    bgf_phase4_update!(energyflux, frictionvel, temperature, soilstate,
        waterfluxbulk, waterstatebulk, waterdiagbulk, patch_data, col_data, lun_data,
        filterp, fn, ustar, um, temp1, temp2, temp12m, temp22m, dth, dqh,
        forc_q_col, forc_rho_col, forc_u_grc, forc_v_grc, forc_pbot_col,
        grnd_ch4_cond_patch, soilevap_beta, soil_resis_sl14, use_lch4, nlevsno)

    return nothing
end
