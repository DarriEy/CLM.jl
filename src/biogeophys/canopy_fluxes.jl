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
# Arrays grouped into device-view structs (Metal 31-buffer limit; the flat form
# had ~41). rah/raw are 2D → CfResOut needs {M,V}. Int gather indices + active
# (Bool) stay loose; the Float64 scalar constants are eltype-converted in the
# wrapper. Launched via manual backend + KA.synchronize (struct-arg kernel).
# rah/raw are 2D → M. grnd_ch4_cond_patch / dleaf_pft come from Float64 kwarg
# defaults and stay Float64 even when the state is ForwardDiff.Dual (AD path) —
# they get their own type params (Vg/Vp) so the struct accepts the mixed eltypes
# (a single V would force all fields to one type and fail to construct under AD).
Base.@kwdef struct CfResOut{M,V,Vg}   # written outputs
    tlbef::V; del2::V; ram1_patch::V; rah::M; raw::M; uaf_patch::V; uuc::V
    dleaf_patch::V; rb::V; rb1_patch::V; grnd_ch4_cond_patch::Vg; svpts::V; eah::V
    rh_af_patch::V; rah1_patch::V; raw1_patch::V; rah2_patch::V; raw2_patch::V; vpd_patch::V
end
Base.@kwdef struct CfResIn{V,Vp}      # read-only float inputs
    t_veg_patch::V; del_arr::V; ustar_patch::V; um_patch::V; temp1::V; temp2::V
    elai_patch::V; esai_patch::V; htop_patch::V; z0mg_col::V; taf_patch::V; qaf_patch::V
    t_grnd_col::V; dleaf_pft::Vp; forc_pbot_col::V; el::V
end
Adapt.@adapt_structure CfResOut
Adapt.@adapt_structure CfResIn

@kernel function _cf_resist_kernel!(o, in, @Const(filterp), @Const(column),
        @Const(gridcell), @Const(itype), @Const(active),
        csoilc_val, use_undercanopy_stability::Bool,
        use_biomass_heat_storage::Bool, use_lch4::Bool,
        cv, a_coef, a_exp, vkc, grav, nu_param, ria_canopy,
        above_canopy::Int, below_canopy::Int)

    fi = @index(Global)
    @inbounds if active[filterp[fi]]
        p = filterp[fi]
        c = column[p]
        T = eltype(o.tlbef)
        ft = itype[p] + 1  # 0-based Fortran PFT → 1-based Julia

        o.tlbef[p] = in.t_veg_patch[p]
        o.del2[p] = in.del_arr[p]

        # Aerodynamic resistances
        o.ram1_patch[p] = one(T) / (in.ustar_patch[p]^2 / in.um_patch[p])
        o.rah[p, above_canopy] = one(T) / (in.temp1[p] * in.ustar_patch[p])
        o.raw[p, above_canopy] = one(T) / (in.temp2[p] * in.ustar_patch[p])

        # Bulk boundary layer resistance
        o.uaf_patch[p] = in.um_patch[p] *
            sqrt(one(T) / (o.ram1_patch[p] * in.um_patch[p]))

        # Empirical undercanopy wind speed
        o.uuc[p] = smooth_min(T(0.4), T(0.03) * in.um_patch[p] / in.ustar_patch[p])

        # Leaf characteristic width
        o.dleaf_patch[p] = in.dleaf_pft[ft]

        cf = cv / (sqrt(o.uaf_patch[p]) * sqrt(o.dleaf_patch[p]))
        o.rb[p] = one(T) / (cf * o.uaf_patch[p])
        o.rb1_patch[p] = o.rb[p]

        # Soil drag coefficient (X. Zeng parameterization)
        w = exp(-(in.elai_patch[p] + in.esai_patch[p]))

        csoilb = vkc / (a_coef * (in.z0mg_col[c] * o.uaf_patch[p] / nu_param)^a_exp)

        ri = (grav * in.htop_patch[p] *
            (in.taf_patch[p] - in.t_grnd_col[c])) /
            (in.taf_patch[p] * o.uaf_patch[p]^2)

        if use_undercanopy_stability && (in.taf_patch[p] - in.t_grnd_col[c]) > zero(T)
            ricsoilc = csoilc_val / (one(T) + ria_canopy * smooth_min(ri, T(10.0)))
            csoilcn = csoilb * w + ricsoilc * (one(T) - w)
        else
            csoilcn = csoilb * w + csoilc_val * (one(T) - w)
        end

        if use_biomass_heat_storage
            o.rah[p, below_canopy] = one(T) / (csoilcn * o.uuc[p])
        else
            o.rah[p, below_canopy] = one(T) / (csoilcn * o.uaf_patch[p])
        end

        o.raw[p, below_canopy] = o.rah[p, below_canopy]

        if use_lch4 && length(o.grnd_ch4_cond_patch) >= p
            o.grnd_ch4_cond_patch[p] = one(T) / (o.raw[p, above_canopy] + o.raw[p, below_canopy])
        end

        # Stomatal resistance intermediates
        o.svpts[p] = in.el[p]
        o.eah[p] = in.forc_pbot_col[c] * in.qaf_patch[p] / T(0.622)
        o.rh_af_patch[p] = o.eah[p] / o.svpts[p]

        # History outputs
        o.rah1_patch[p] = o.rah[p, above_canopy]
        o.raw1_patch[p] = o.raw[p, above_canopy]
        o.rah2_patch[p] = o.rah[p, below_canopy]
        o.raw2_patch[p] = o.raw[p, below_canopy]
        o.vpd_patch[p]  = smooth_max((o.svpts[p] - o.eah[p]), T(50.0)) * T(0.001)
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
        patch_data, filterp, fn::Int, active,
        temp1, temp2, tlbef, del2, del_arr, rah, raw, uuc, rb, svpts, eah, el,
        dleaf_pft, grnd_ch4_cond_patch, forc_pbot_col,
        csoilc_val, use_undercanopy_stability::Bool,
        use_biomass_heat_storage::Bool, use_lch4::Bool, params)
    o = CfResOut(; tlbef = tlbef, del2 = del2, ram1_patch = frictionvel.ram1_patch,
        rah = rah, raw = raw, uaf_patch = frictionvel.uaf_patch, uuc = uuc,
        dleaf_patch = canopystate.dleaf_patch, rb = rb, rb1_patch = frictionvel.rb1_patch,
        grnd_ch4_cond_patch = grnd_ch4_cond_patch, svpts = svpts, eah = eah,
        rh_af_patch = waterdiagbulk.rh_af_patch, rah1_patch = frictionvel.rah1_patch,
        raw1_patch = frictionvel.raw1_patch, rah2_patch = frictionvel.rah2_patch,
        raw2_patch = frictionvel.raw2_patch, vpd_patch = frictionvel.vpd_patch)
    inp = CfResIn(; t_veg_patch = temperature.t_veg_patch, del_arr = del_arr,
        ustar_patch = frictionvel.ustar_patch, um_patch = frictionvel.um_patch,
        temp1 = temp1, temp2 = temp2, elai_patch = canopystate.elai_patch,
        esai_patch = canopystate.esai_patch, htop_patch = canopystate.htop_patch,
        z0mg_col = frictionvel.z0mg_col, taf_patch = frictionvel.taf_patch,
        qaf_patch = frictionvel.qaf_patch, t_grnd_col = temperature.t_grnd_col,
        dleaf_pft = dleaf_pft, forc_pbot_col = forc_pbot_col, el = el)
    T = eltype(tlbef)
    let be = _kernel_backend(tlbef)
        _cf_resist_kernel!(be)(o, inp, filterp, patch_data.column, patch_data.gridcell,
            patch_data.itype, active, T(csoilc_val), use_undercanopy_stability,
            use_biomass_heat_storage, use_lch4, T(params.cv), T(params.a_coef),
            T(params.a_exp), T(VKC), T(GRAV), T(NU_PARAM), T(RIA_CANOPY),
            ABOVE_CANOPY, BELOW_CANOPY; ndrange = fn)
        KA.synchronize(be)
    end
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
# Arrays grouped into device-view structs (Metal 31-buffer limit; the flat form
# had ~88 args). All energy arrays are state/scratch/forcing → uniform eltype
# even on the AD path (no Float64-kwarg fields like resist's dleaf_pft), so
# single {V}/{M} params suffice. nmozsgn + the Int gather/index arrays stay
# loose. Float64 scalar constants are eltype-converted in cf_energy_update!.
# Inside, the struct fields are aliased to locals (zero-cost array bindings) so
# the physics body reads unchanged except Float64 literals → T(). Launched via
# manual backend + KA.synchronize (struct-arg kernel).
# Metal allows at most ~31 array resources per indirect-argument-buffer (struct),
# so the 34-field output/patch-input bags are each split in two (≤17 fields each).
Base.@kwdef struct CfEnOutA{V}    # written per-patch float outputs (1/2)
    wtg::V; wtl0::V; wta0::V; wtstem0::V; wtga::V; wtal::V; lw_stem::V; lw_leaf::V
    wtgq::V; wtlq0::V; wtaq0::V; wtalq::V; efe::V; dt_veg::V; del_arr::V; err_arr::V; qsatl::V
end
Base.@kwdef struct CfEnOutB{V}    # written per-patch float outputs (2/2)
    el::V; qsatldT::V; dth::V; dqh::V; delq::V; obuold::V; canopy_cond_patch::V
    eflx_sh_stem_patch::V; eflx_sh_veg_patch::V; qflx_tran_veg_patch::V
    qflx_evap_veg_patch::V; t_veg_patch::V; taf_patch::V; qaf_patch::V
    zeta_patch::V; um_patch::V; obu_patch::V
end
Base.@kwdef struct CfEnPA{V}      # read-only per-patch floats (1/2)
    rb::V; rstem::V; sa_leaf::V; sa_stem::V; sa_internal::V; frac_rad_abs_by_stem::V
    air::V; bir::V; cir::V; cp_leaf::V; tl_ini::V; tlbef::V; zldis::V; temp1::V; temp2::V
    ur::V; efeb::V
end
Base.@kwdef struct CfEnPB{V}      # read-only per-patch floats (2/2)
    emv_patch::V; t_stem_patch::V; btran_patch::V; fdry_patch::V
    fwet_patch::V; elai_patch::V; esai_patch::V; laisun_patch::V; laisha_patch::V
    rssun_patch::V; rssha_patch::V; thm_patch::V; sabv_patch::V; liqcan_patch::V
    snocan_patch::V; uaf_patch::V; ustar_patch::V
end
Base.@kwdef struct CfEnC{V}       # read-only per-column floats
    qg_col::V; frac_sno_eff_col::V; frac_h2osfc_col::V; t_h2osfc_col::V
    snow_depth_col::V; soilbeta_col::V; soilresis_col::V; forc_rho_col::V
    forc_q_col::V; forc_th_col::V; thv_col::V; t_grnd_col::V; forc_pbot_col::V
end
Base.@kwdef struct CfEnM{M}       # read-only 2D arrays
    rah::M; raw::M; t_soisno_col::M
end
Adapt.@adapt_structure CfEnOutA
Adapt.@adapt_structure CfEnOutB
Adapt.@adapt_structure CfEnPA
Adapt.@adapt_structure CfEnPB
Adapt.@adapt_structure CfEnC
Adapt.@adapt_structure CfEnM

# Scalar constants bundled into one isbits arg (Metal caps total kernel args
# ~31; 13 separate Float scalars would each count). All fields share precision S.
struct CfEnScalars{S}
    sb::S; cpair::S; hvap::S; vkc::S; grav::S; btran0::S; delmax_canopy::S
    zii_canopy::S; beta_canopy::S; dtime::S; z_dl::S; lai_dl::S; zetamaxstable::S
end

@kernel function _cf_energy_kernel!(o1, o2, ep1, ep2, ec, em, sc, nmozsgn,
        @Const(filterp), @Const(active), @Const(column),
        @Const(snl), @Const(frac_veg_nosno_patch),
        # scalars
        soilevap_beta::Bool, soil_resis_sl14::Bool,
        use_lch4::Bool, use_hydrstress::Bool, canopy_cond_len::Int,
        above_canopy::Int, below_canopy::Int, nlevsno::Int)

    fi = @index(Global)
    @inbounds if active[filterp[fi]]
        p = filterp[fi]
        c = column[p]
        T = eltype(o1.wtg)
        # Alias grouped struct fields to locals (zero-cost array bindings) so the
        # physics body below is unchanged except for Float64 literals → T(...).
        wtg = o1.wtg; wtl0 = o1.wtl0; wta0 = o1.wta0; wtstem0 = o1.wtstem0; wtga = o1.wtga
        wtal = o1.wtal; lw_stem = o1.lw_stem; lw_leaf = o1.lw_leaf; wtgq = o1.wtgq
        wtlq0 = o1.wtlq0; wtaq0 = o1.wtaq0; wtalq = o1.wtalq; efe = o1.efe; dt_veg = o1.dt_veg
        del_arr = o1.del_arr; err_arr = o1.err_arr; qsatl = o1.qsatl; el = o2.el
        qsatldT = o2.qsatldT; dth = o2.dth; dqh = o2.dqh; delq = o2.delq; obuold = o2.obuold
        canopy_cond_patch = o2.canopy_cond_patch; eflx_sh_stem_patch = o2.eflx_sh_stem_patch
        eflx_sh_veg_patch = o2.eflx_sh_veg_patch; qflx_tran_veg_patch = o2.qflx_tran_veg_patch
        qflx_evap_veg_patch = o2.qflx_evap_veg_patch; t_veg_patch = o2.t_veg_patch
        taf_patch = o2.taf_patch; qaf_patch = o2.qaf_patch; zeta_patch = o2.zeta_patch
        um_patch = o2.um_patch; obu_patch = o2.obu_patch
        rb = ep1.rb; rstem = ep1.rstem; sa_leaf = ep1.sa_leaf; sa_stem = ep1.sa_stem
        sa_internal = ep1.sa_internal; frac_rad_abs_by_stem = ep1.frac_rad_abs_by_stem
        air = ep1.air; bir = ep1.bir; cir = ep1.cir; cp_leaf = ep1.cp_leaf; tl_ini = ep1.tl_ini
        tlbef = ep1.tlbef; zldis = ep1.zldis; temp1 = ep1.temp1; temp2 = ep1.temp2; ur = ep1.ur
        efeb = ep1.efeb; emv_patch = ep2.emv_patch; t_stem_patch = ep2.t_stem_patch
        btran_patch = ep2.btran_patch; fdry_patch = ep2.fdry_patch; fwet_patch = ep2.fwet_patch
        elai_patch = ep2.elai_patch; esai_patch = ep2.esai_patch; laisun_patch = ep2.laisun_patch
        laisha_patch = ep2.laisha_patch; rssun_patch = ep2.rssun_patch; rssha_patch = ep2.rssha_patch
        thm_patch = ep2.thm_patch; sabv_patch = ep2.sabv_patch; liqcan_patch = ep2.liqcan_patch
        snocan_patch = ep2.snocan_patch; uaf_patch = ep2.uaf_patch; ustar_patch = ep2.ustar_patch
        qg_col = ec.qg_col; frac_sno_eff_col = ec.frac_sno_eff_col
        frac_h2osfc_col = ec.frac_h2osfc_col; t_h2osfc_col = ec.t_h2osfc_col
        snow_depth_col = ec.snow_depth_col; soilbeta_col = ec.soilbeta_col
        soilresis_col = ec.soilresis_col; forc_rho_col = ec.forc_rho_col
        forc_q_col = ec.forc_q_col; forc_th_col = ec.forc_th_col; thv_col = ec.thv_col
        t_grnd_col = ec.t_grnd_col; forc_pbot_col = ec.forc_pbot_col
        rah = em.rah; raw = em.raw; t_soisno_col = em.t_soisno_col
        sb = sc.sb; cpair = sc.cpair; hvap = sc.hvap; vkc = sc.vkc; grav = sc.grav
        btran0 = sc.btran0; delmax_canopy = sc.delmax_canopy; zii_canopy = sc.zii_canopy
        beta_canopy = sc.beta_canopy; dtime = sc.dtime; z_dl = sc.z_dl; lai_dl = sc.lai_dl
        zetamaxstable = sc.zetamaxstable

        # Sensible heat conductance for air, leaf, ground, stem
        wta  = one(T) / rah[p, above_canopy]       # air
        wtl  = sa_leaf[p] / rb[p]                # leaf
        wtg[p] = one(T) / rah[p, below_canopy]     # ground
        wtstem = sa_stem[p] / (rstem[p] + rb[p]) # stem

        wtshi = one(T) / (wta + wtl + wtstem + wtg[p])

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

        if fdry_p > zero(T)
            rppdry = fdry_p * rb[p] *
                (laisun_p / (rb[p] + rssun_p) + laisha_p / (rb[p] + rssha_p)) / elai_p
        else
            rppdry = zero(T)
        end

        # Canopy conductance for methane
        if use_lch4 && canopy_cond_len >= p
            canopy_cond_patch[p] = (laisun_p / (rb[p] + rssun_p) +
                laisha_p / (rb[p] + rssha_p)) / smooth_max(elai_p, T(0.01))
        end

        efpot = forc_rho_col[c] * ((elai_p + esai_p) / rb[p]) *
            (qsatl[p] - qaf_patch[p])
        h2ocan = liqcan_patch[p] + snocan_patch[p]

        fwet_p = fwet_patch[p]
        btran_p = btran_patch[p]
        qflx_tran_veg_p = qflx_tran_veg_patch[p]

        if use_hydrstress
            if efpot > zero(T)
                if btran_p > btran0
                    rpp = rppdry + fwet_p
                else
                    rpp = fwet_p
                end
                rpp = smooth_min(rpp, (qflx_tran_veg_p + h2ocan / dtime) / efpot)
            else
                rpp = one(T)
            end
        else
            if efpot > zero(T)
                if btran_p > btran0
                    qflx_tran_veg_patch[p] = efpot * rppdry
                    rpp = rppdry + fwet_p
                else
                    rpp = fwet_p
                    qflx_tran_veg_patch[p] = zero(T)
                end
                rpp = min(rpp, (qflx_tran_veg_patch[p] + h2ocan / dtime) / efpot)
            else
                rpp = one(T)
                qflx_tran_veg_patch[p] = zero(T)
            end
        end

        # Latent heat conductances
        fvn = frac_veg_nosno_patch[p]
        wtaq  = fvn / raw[p, above_canopy]
        wtlq  = fvn * (elai_p + esai_p) / rb[p] * rpp

        # Litter layer resistance (Sakaguchi)
        snow_depth_c = z_dl
        fsno_dl = snow_depth_col[c] / snow_depth_c
        elai_dl = lai_dl * (one(T) - smooth_min(fsno_dl, one(T)))
        rdl = (one(T) - exp(-elai_dl)) / (T(0.004) * uaf_patch[p])

        if delq[p] < zero(T)
            wtgq[p] = fvn / (raw[p, below_canopy] + rdl)
        else
            if soilevap_beta
                wtgq[p] = soilbeta_col[c] * fvn / (raw[p, below_canopy] + rdl)
            end
            if soil_resis_sl14
                wtgq[p] = fvn / (raw[p, below_canopy] + soilresis_col[c])
            end
        end

        wtsqi  = one(T) / (wtaq + wtlq + wtgq[p])
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
        erre = zero(T)
        if efe[p] * efeb[p] < zero(T)
            efeold = efe[p]
            efe[p] = T(0.1) * efeold
            erre = efe[p] - efeold
        end

        # Fractionate ground emitted longwave
        snl_c = snl[c]
        lw_grnd = (frac_sno_eff_col[c] * t_soisno_col[c, snl_c + 1 + nlevsno]^4 +
            (one(T) - frac_sno_eff_col[c] - frac_h2osfc_col[c]) *
                t_soisno_col[c, 1 + nlevsno]^4 +
            frac_h2osfc_col[c] * t_h2osfc_col[c]^4)

        # Newton-Raphson: dt_veg
        _numer = ((one(T) - frac_rad_abs_by_stem[p]) * (sabv_patch[p] + air[p] +
            bir[p] * t_veg_patch[p]^4 + cir[p] * lw_grnd) -
            efsh - efe[p] - lw_leaf[p] + lw_stem[p] -
            (cp_leaf[p] / dtime) * (t_veg_patch[p] - tl_ini[p]))
        _denom = ((one(T) - frac_rad_abs_by_stem[p]) * (-T(4.0) * bir[p] * t_veg_patch[p]^3) +
             T(4.0) * sa_internal[p] * emv_patch[p] * sb * t_veg_patch[p]^3 +
             dc1 * wtga[p] + dc2 * wtgaq * qsatldT[p] + cp_leaf[p] / dtime)
        dt_veg[p] = _numer / _denom

        t_veg_patch[p] = tlbef[p] + dt_veg[p]

        dels = dt_veg[p]
        del_arr[p] = abs(dels)
        err_arr[p] = zero(T)
        if del_arr[p] > delmax_canopy
            dt_veg[p] = delmax_canopy * dels / del_arr[p]
            t_veg_patch[p] = tlbef[p] + dt_veg[p]
            err_arr[p] = (one(T) - frac_rad_abs_by_stem[p]) * (sabv_patch[p] + air[p] +
                bir[p] * tlbef[p]^3 * (tlbef[p] + T(4.0) * dt_veg[p]) + cir[p] * lw_grnd) -
                sa_internal[p] * emv_patch[p] * sb * tlbef[p]^3 * (tlbef[p] + T(4.0) * dt_veg[p]) +
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
            ecidif = max(zero(T), qflx_evap_veg_patch[p] -
                qflx_tran_veg_patch[p] - h2ocan / dtime)
            qflx_evap_veg_patch[p] = min(qflx_evap_veg_patch[p],
                qflx_tran_veg_patch[p] + h2ocan / dtime)
        else
            ecidif = zero(T)
            if efpot > zero(T) && btran_patch[p] > btran0
                qflx_tran_veg_patch[p] = efpot * rppdry
            else
                qflx_tran_veg_patch[p] = zero(T)
            end
            ecidif = max(zero(T), qflx_evap_veg_patch[p] -
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
            tlbef[p]^3 * (tlbef[p] + T(4.0) * dt_veg[p])

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
        thvstar = tstar * (one(T) + T(0.61) * forc_q_col[c]) + T(0.61) * forc_th_col[c] * qstar
        zeta_patch[p] = zldis[p] * vkc * grav * thvstar /
            (ustar_patch[p]^2 * thv_col[c])

        if zeta_patch[p] >= zero(T)  # stable
            zeta_patch[p] = smooth_clamp(zeta_patch[p], T(0.01), zetamaxstable)
            um_patch[p] = smooth_max(ur[p], T(0.1))
        else  # unstable
            zeta_patch[p] = smooth_clamp(zeta_patch[p], T(-100.0), T(-0.01))
            if ustar_patch[p] * thvstar > zero(T)
                wc = zero(T)
            else
                wc_arg = smooth_max(-grav * ustar_patch[p] * thvstar *
                    zii_canopy / thv_col[c], zero(T))
                wc = beta_canopy * wc_arg^T(0.333)
            end
            um_patch[p] = sqrt(ur[p]^2 + wc^2)
        end
        obu_patch[p] = zldis[p] / zeta_patch[p]

        if obuold[p] * obu_patch[p] < zero(T)
            nmozsgn[p] += 1
        end
        if nmozsgn[p] >= 4
            obu_patch[p] = zldis[p] / (-T(0.01))
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
        photosyns, patch_data, col_data, filterp, fn::Int, active,
        rah, raw, rb, rstem, sa_leaf, sa_stem, sa_internal,
        frac_rad_abs_by_stem, air, bir, cir, cp_leaf, tl_ini, tlbef, zldis,
        temp1, temp2, ur, efeb, wtg, wtl0, wta0, wtstem0, wtga, wtal,
        lw_stem, lw_leaf, wtgq, wtlq0, wtaq0, wtalq, efe, dt_veg, del_arr,
        err_arr, qsatl, el, qsatldT, dth, dqh, delq, obuold, nmozsgn,
        forc_q_col, forc_th_col, forc_pbot_col, forc_rho_col,
        soilevap_beta::Bool, soil_resis_sl14::Bool, use_lch4::Bool,
        use_hydrstress::Bool, nlevsno::Int, dtime, params)
    o1 = CfEnOutA(; wtg = wtg, wtl0 = wtl0, wta0 = wta0, wtstem0 = wtstem0, wtga = wtga,
        wtal = wtal, lw_stem = lw_stem, lw_leaf = lw_leaf, wtgq = wtgq, wtlq0 = wtlq0,
        wtaq0 = wtaq0, wtalq = wtalq, efe = efe, dt_veg = dt_veg, del_arr = del_arr,
        err_arr = err_arr, qsatl = qsatl)
    o2 = CfEnOutB(; el = el, qsatldT = qsatldT, dth = dth, dqh = dqh, delq = delq,
        obuold = obuold, canopy_cond_patch = energyflux.canopy_cond_patch,
        eflx_sh_stem_patch = energyflux.eflx_sh_stem_patch,
        eflx_sh_veg_patch = energyflux.eflx_sh_veg_patch,
        qflx_tran_veg_patch = waterfluxbulk.wf.qflx_tran_veg_patch,
        qflx_evap_veg_patch = waterfluxbulk.wf.qflx_evap_veg_patch,
        t_veg_patch = temperature.t_veg_patch, taf_patch = frictionvel.taf_patch,
        qaf_patch = frictionvel.qaf_patch, zeta_patch = frictionvel.zeta_patch,
        um_patch = frictionvel.um_patch, obu_patch = frictionvel.obu_patch)
    ep1 = CfEnPA(; rb = rb, rstem = rstem, sa_leaf = sa_leaf, sa_stem = sa_stem,
        sa_internal = sa_internal, frac_rad_abs_by_stem = frac_rad_abs_by_stem,
        air = air, bir = bir, cir = cir, cp_leaf = cp_leaf, tl_ini = tl_ini,
        tlbef = tlbef, zldis = zldis, temp1 = temp1, temp2 = temp2, ur = ur, efeb = efeb)
    ep2 = CfEnPB(; emv_patch = temperature.emv_patch, t_stem_patch = temperature.t_stem_patch,
        btran_patch = energyflux.btran_patch, fdry_patch = waterdiagbulk.fdry_patch,
        fwet_patch = waterdiagbulk.fwet_patch, elai_patch = canopystate.elai_patch,
        esai_patch = canopystate.esai_patch, laisun_patch = canopystate.laisun_patch,
        laisha_patch = canopystate.laisha_patch, rssun_patch = photosyns.rssun_patch,
        rssha_patch = photosyns.rssha_patch, thm_patch = temperature.thm_patch,
        sabv_patch = solarabs.sabv_patch, liqcan_patch = waterstatebulk.ws.liqcan_patch,
        snocan_patch = waterstatebulk.ws.snocan_patch, uaf_patch = frictionvel.uaf_patch,
        ustar_patch = frictionvel.ustar_patch)
    ec = CfEnC(; qg_col = waterdiagbulk.qg_col,
        frac_sno_eff_col = waterdiagbulk.frac_sno_eff_col,
        frac_h2osfc_col = waterdiagbulk.frac_h2osfc_col,
        t_h2osfc_col = temperature.t_h2osfc_col, snow_depth_col = waterdiagbulk.snow_depth_col,
        soilbeta_col = soilstate.soilbeta_col, soilresis_col = soilstate.soilresis_col,
        forc_rho_col = forc_rho_col, forc_q_col = forc_q_col, forc_th_col = forc_th_col,
        thv_col = temperature.thv_col, t_grnd_col = temperature.t_grnd_col,
        forc_pbot_col = forc_pbot_col)
    em = CfEnM(; rah = rah, raw = raw, t_soisno_col = temperature.t_soisno_col)
    T = eltype(wtg)
    sc = CfEnScalars{T}(T(SB), T(CPAIR), T(HVAP), T(VKC), T(GRAV), T(BTRAN0),
        T(DELMAX_CANOPY), T(ZII_CANOPY), T(BETA_CANOPY), T(dtime), T(params.z_dl),
        T(params.lai_dl), T(frictionvel.zetamaxstable))
    let be = _kernel_backend(wtg)
        _cf_energy_kernel!(be)(o1, o2, ep1, ep2, ec, em, sc, nmozsgn,
            filterp, active, patch_data.column,
            col_data.snl, canopystate.frac_veg_nosno_patch,
            soilevap_beta, soil_resis_sl14, use_lch4, use_hydrstress,
            length(energyflux.canopy_cond_patch), ABOVE_CANOPY, BELOW_CANOPY,
            nlevsno; ndrange = fn)
        KA.synchronize(be)
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Convergence test (cf-newton): replaces the host-side filter compaction inside
# the Newton stability iteration. Per-filtered-patch; for each still-active
# patch it computes the energy-flux / temperature convergence metrics and
# clears `active[p]` once both fall below tolerance. The filter stays fixed at
# the original patch set; converged patches are skipped on subsequent
# iterations via the in-kernel `active[p]` guard rather than being compacted
# out, so the physics arrays never need a device→host copy mid-iteration. The
# caller only launches this past ITMIN_CANOPY (matching the scalar version).
# ---------------------------------------------------------------------------
@kernel function _cf_converge_kernel!(active, efeb, dele, det_arr, num_iter_patch,
        @Const(filterp), @Const(efe), @Const(del_arr), @Const(del2),
        niter, dtmin, dlemin)
    fi = @index(Global)
    @inbounds if active[filterp[fi]]
        p = filterp[fi]
        dele[p] = abs(efe[p] - efeb[p])
        efeb[p] = efe[p]
        det_arr[p] = max(del_arr[p], del2[p])
        num_iter_patch[p] = niter
        if det_arr[p] < dtmin && dele[p] < dlemin
            active[p] = false
        end
    end
end

"""
    cf_converge_update!(active, efeb, dele, det_arr, num_iter_patch,
                        filterp, fn, efe, del_arr, del2, itlef)

Launch the canopy-fluxes convergence-test kernel over the `fn` filtered patches.
Updates the convergence metrics for still-active patches and clears `active[p]`
on convergence. Backend-agnostic (CPU loop or GPU); one thread per filtered
patch.
"""
function cf_converge_update!(active, efeb, dele, det_arr, num_iter_patch,
        filterp, fn::Int, efe, del_arr, del2, itlef::Int)
    _launch!(_cf_converge_kernel!, active, efeb, dele, det_arr, num_iter_patch,
        filterp, efe, del_arr, del2,
        convert(eltype(num_iter_patch), itlef),
        convert(eltype(efe), DTMIN_CANOPY), convert(eltype(efe), DLEMIN_CANOPY);
        ndrange = fn)
    return nothing
end

# ---------------------------------------------------------------------------
# Monin-Obukhov length initialization (cf-newton): per-filtered-patch. Inlines
# monin_obuk_ini with eltype-generic literals + grav passed as a scalar (the
# host monin_obuk_ini uses Float64 literals and the GRAV module global, neither
# Metal-safe). Sets um/obu/num_iter and the leaf/stem initial temperatures.
# (The host monin_obuk_ini also sets an unused local `ustar`; dropped here — it
# never affects the (um,obu) result.)
# ---------------------------------------------------------------------------
@kernel function _cf_moninobukini_kernel!(um_patch, obu_patch, num_iter_patch,
        tl_ini, ts_ini,
        @Const(filterp), @Const(column), @Const(ur), @Const(thv_col),
        @Const(dthv), @Const(zldis), @Const(z0mv_patch),
        @Const(t_veg_patch), @Const(t_stem_patch),
        zetamaxstable, grav)
    fi = @index(Global)
    @inbounds begin
        p = filterp[fi]
        c = column[p]
        T = eltype(um_patch)
        ur_p    = ur[p]
        dthv_p  = dthv[p]
        zldis_p = zldis[p]
        z0m_p   = z0mv_patch[p]

        wc = T(0.5)
        um = dthv_p >= zero(T) ? max(ur_p, T(0.1)) : sqrt(ur_p * ur_p + wc * wc)
        rib = grav * zldis_p * dthv_p / (thv_col[c] * um * um)
        if rib >= zero(T)
            zeta = rib * log(zldis_p / z0m_p) / (one(T) - T(5.0) * min(rib, T(0.19)))
            zeta = min(zetamaxstable, max(zeta, T(0.01)))
        else
            zeta = rib * log(zldis_p / z0m_p)
            zeta = max(T(-100.0), min(zeta, T(-0.01)))
        end

        um_patch[p] = um
        obu_patch[p] = zldis_p / zeta
        num_iter_patch[p] = zero(eltype(num_iter_patch))
        tl_ini[p] = t_veg_patch[p]
        ts_ini[p] = t_stem_patch[p]
    end
end

function cf_moninobukini_update!(frictionvel, temperature, patch_data, filterp,
        fn::Int, ur, dthv, zldis, tl_ini, ts_ini)
    _launch!(_cf_moninobukini_kernel!, frictionvel.um_patch, frictionvel.obu_patch,
        frictionvel.num_iter_patch, tl_ini, ts_ini,
        filterp, patch_data.column, ur, temperature.thv_col, dthv, zldis,
        frictionvel.z0mv_patch, temperature.t_veg_patch, temperature.t_stem_patch,
        convert(eltype(frictionvel.um_patch), frictionvel.zetamaxstable),
        convert(eltype(frictionvel.um_patch), GRAV); ndrange = fn)
    return nothing
end

# ---------------------------------------------------------------------------
# Device-safe (filterp-gather) variants of photosynthesis_timestep_init! and
# photosynthesis_total! for the canopy_fluxes! call sites. canopy calls both
# with default flags (use_c13/use_c14/use_cn/use_fates = false), so only the
# base path runs; the filter (built from mask_exposedvegp) selects exactly the
# same patches as the host mask guard. The original host functions in
# photosynthesis.jl are kept for their other callers.
# ---------------------------------------------------------------------------
@kernel function _cf_psn_init_kernel!(psnsun, psnsun_wc, psnsun_wj, psnsun_wp,
        psnsha, psnsha_wc, psnsha_wj, psnsha_wp, fpsn, fpsn_wc, fpsn_wj, fpsn_wp,
        @Const(filterp))
    fi = @index(Global)
    @inbounds begin
        p = filterp[fi]
        z = zero(eltype(psnsun))
        psnsun[p] = z; psnsun_wc[p] = z; psnsun_wj[p] = z; psnsun_wp[p] = z
        psnsha[p] = z; psnsha_wc[p] = z; psnsha_wj[p] = z; psnsha_wp[p] = z
        fpsn[p]   = z; fpsn_wc[p]   = z; fpsn_wj[p]   = z; fpsn_wp[p]   = z
    end
end

function cf_psn_init_update!(ps, filterp, fn::Int)
    _launch!(_cf_psn_init_kernel!, ps.psnsun_patch, ps.psnsun_wc_patch,
        ps.psnsun_wj_patch, ps.psnsun_wp_patch, ps.psnsha_patch, ps.psnsha_wc_patch,
        ps.psnsha_wj_patch, ps.psnsha_wp_patch, ps.fpsn_patch, ps.fpsn_wc_patch,
        ps.fpsn_wj_patch, ps.fpsn_wp_patch, filterp; ndrange = fn)
    return nothing
end

@kernel function _cf_psn_total_kernel!(fpsn, fpsn_wc, fpsn_wj, fpsn_wp,
        @Const(filterp), @Const(psnsun), @Const(psnsun_wc), @Const(psnsun_wj),
        @Const(psnsun_wp), @Const(psnsha), @Const(psnsha_wc), @Const(psnsha_wj),
        @Const(psnsha_wp), @Const(laisun), @Const(laisha))
    fi = @index(Global)
    @inbounds begin
        p = filterp[fi]
        ls = laisun[p]; lh = laisha[p]
        fpsn[p]    = psnsun[p]    * ls + psnsha[p]    * lh
        fpsn_wc[p] = psnsun_wc[p] * ls + psnsha_wc[p] * lh
        fpsn_wj[p] = psnsun_wj[p] * ls + psnsha_wj[p] * lh
        fpsn_wp[p] = psnsun_wp[p] * ls + psnsha_wp[p] * lh
    end
end

function cf_psn_total_update!(ps, laisun, laisha, filterp, fn::Int)
    _launch!(_cf_psn_total_kernel!, ps.fpsn_patch, ps.fpsn_wc_patch, ps.fpsn_wj_patch,
        ps.fpsn_wp_patch, filterp, ps.psnsun_patch, ps.psnsun_wc_patch,
        ps.psnsun_wj_patch, ps.psnsun_wp_patch, ps.psnsha_patch, ps.psnsha_wc_patch,
        ps.psnsha_wj_patch, ps.psnsha_wp_patch, laisun, laisha; ndrange = fn)
    return nothing
end

# =====================================================================
# canopy_fluxes! Phase-1 init kernels (cf1). Per-filtered-patch (fi → p=filterp[fi]),
# elementwise; literals eltype-converted. The Newton iteration is fully
# kernelized: a per-patch `active` mask (cleared by _cf_converge_kernel!)
# replaces host-side filter compaction, so the whole iteration runs without a
# device→host sync of the physics arrays.
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

# Aerodynamic Z0 / displacement-height adjustment (cf3). z0param_method resolved to
# an Int code on the host (1=ZengWang2007, 2=Meier2022; unknown → host error), so no
# String compare or error() runs in the kernel. The Meier2022 branch's iterative
# U_ustar solve runs as an in-thread while loop. Constants eltype-converted.
@kernel function _cf_z0_kernel!(displa, z0mv, z0hv, z0qv, forc_hgt_u, forc_hgt_t,
        forc_hgt_q, @Const(filterp), @Const(column), @Const(gridcell), @Const(itype),
        @Const(elai), @Const(esai), @Const(htop), @Const(z0mg), @Const(forc_hgt_u_grc),
        @Const(forc_hgt_t_grc), @Const(forc_hgt_q_grc), @Const(z0v_LAImax_pft),
        @Const(z0v_Cs_pft), @Const(z0v_Cr_pft), @Const(z0v_c_pft), @Const(z0v_cw_pft),
        z0_method::Int)
    fi = @index(Global)
    @inbounds begin
        p = filterp[fi]
        T = eltype(displa)
        c = column[p]; g = gridcell[p]; ft = itype[p] + 1
        if z0_method == 1   # ZengWang2007
            lt = smooth_min(elai[p] + esai[p], T(TLSAI_CRIT))
            egvf = (one(T) - T(ALPHA_AERO) * exp(-lt)) / (one(T) - T(ALPHA_AERO) * exp(-T(TLSAI_CRIT)))
            displa[p] = egvf * displa[p]
            z0mv[p] = exp(egvf * log(z0mv[p]) + (one(T) - egvf) * log(z0mg[c]))
        else                # Meier2022 (z0_method == 2; host-validated)
            lt = smooth_max(T(1.0e-5), elai[p] + esai[p])
            displa[p] = htop[p] *
                (one(T) - (one(T) - exp(-(T(CD1_PARAM) * lt)^T(0.5))) / (T(CD1_PARAM) * lt)^T(0.5))
            lt = smooth_min(lt, z0v_LAImax_pft[ft])
            delt_iter = T(2.0)
            U_ustar_ini = (z0v_Cs_pft[ft] + z0v_Cr_pft[ft] * lt * T(0.5))^T(-0.5) *
                z0v_c_pft[ft] * lt * T(0.25)
            U_ustar = U_ustar_ini
            while delt_iter > T(1.0e-4)
                U_ustar_prev = U_ustar
                U_ustar = U_ustar_ini * exp(U_ustar_prev)
                delt_iter = abs(U_ustar - U_ustar_prev)
            end
            U_ustar = T(4.0) * U_ustar / lt / z0v_c_pft[ft]
            z0mv[p] = htop[p] * (one(T) - displa[p] / htop[p]) *
                exp(-T(VKC) * U_ustar + log(z0v_cw_pft[ft]) - one(T) + z0v_cw_pft[ft]^T(-1.0))
        end
        z0hv[p] = z0mv[p]
        z0qv[p] = z0mv[p]
        forc_hgt_u[p] = forc_hgt_u_grc[g] + z0mv[p] + displa[p]
        forc_hgt_t[p] = forc_hgt_t_grc[g] + z0hv[p] + displa[p]
        forc_hgt_q[p] = forc_hgt_q_grc[g] + z0qv[p] + displa[p]
    end
end

# Net absorbed longwave + QSat + CO2/O2 + flux-profile init (cf4). Per-filtered-patch;
# qsat() is kernel-safe (4-tuple unpack preserved). The forcing-height-below-canopy
# warning stays on the host (scans zldis after the launch — diagnostic, not in kernel).
# Arrays are grouped into device-view structs (Metal caps a kernel at 31 buffer args;
# the flat form had 38). nmozsgn (Int) + the Int gather indices stay loose. Launched
# via manual backend + KA.synchronize (struct-arg kernel, like _cf_postiter_kernel!).
Base.@kwdef struct CfLwqOut{V}   # written per-patch float outputs
    air::V; bir::V; cir::V; qsatl::V; el::V; qsatldT::V; co2_arr::V; o2_arr::V
    taf::V; qaf::V; ur::V; dth::V; dqh::V; delq::V; dthv::V; zldis::V
end
Base.@kwdef struct CfLwqP{V}     # read-only per-patch floats
    emv::V; t_veg::V; thm::V; forc_hgt_u::V; displa::V
end
Base.@kwdef struct CfLwqC{V}     # read-only per-column floats
    emg::V; forc_lwrad::V; forc_pbot::V; t_grnd::V; forc_q::V; qg::V; forc_th::V
end
Base.@kwdef struct CfLwqG{V}     # read-only per-gridcell floats
    forc_pco2::V; forc_po2::V; forc_u::V; forc_v::V
end
Adapt.@adapt_structure CfLwqOut
Adapt.@adapt_structure CfLwqP
Adapt.@adapt_structure CfLwqC
Adapt.@adapt_structure CfLwqG

@kernel function _cf_longwave_qsat_kernel!(out, pp, cc, gg, nmozsgn,
        @Const(filterp), @Const(column), @Const(gridcell), wind_min)
    fi = @index(Global)
    @inbounds begin
        p = filterp[fi]
        T = eltype(out.air)
        c = column[p]; g = gridcell[p]
        emv_p = pp.emv[p]; emg_c = cc.emg[c]
        out.air[p] = emv_p * (one(T) + (one(T) - emv_p) * (one(T) - emg_c)) * cc.forc_lwrad[c]
        out.bir[p] = -(T(2.0) - emv_p * (one(T) - emg_c)) * emv_p * T(SB)
        out.cir[p] = emv_p * emg_c * T(SB)
        (qs_tmp, es_tmp, dqsdT_tmp, _) = qsat(pp.t_veg[p], cc.forc_pbot[c])
        out.qsatl[p] = qs_tmp; out.el[p] = es_tmp; out.qsatldT[p] = dqsdT_tmp
        out.co2_arr[p] = gg.forc_pco2[g]; out.o2_arr[p] = gg.forc_po2[g]
        nmozsgn[p] = 0
        out.taf[p] = (cc.t_grnd[c] + pp.thm[p]) / T(2.0)
        out.qaf[p] = (cc.forc_q[c] + cc.qg[c]) / T(2.0)
        out.ur[p] = smooth_max(wind_min, sqrt(gg.forc_u[g]^2 + gg.forc_v[g]^2))
        out.dth[p] = pp.thm[p] - out.taf[p]
        out.dqh[p] = cc.forc_q[c] - out.qaf[p]
        out.delq[p] = cc.qg[c] - out.qaf[p]
        out.dthv[p] = out.dth[p] * (one(T) + T(0.61) * cc.forc_q[c]) + T(0.61) * cc.forc_th[c] * out.dqh[p]
        out.zldis[p] = pp.forc_hgt_u[p] - pp.displa[p]
    end
end

# Post-iteration diagnostics (cf5): the field-heaviest canopy loop (~84 arrays),
# grouped into device-view structs (the soil_temperature! template). Per-filtered-
# patch; the cgrnds/cgrndl/liqcan/snocan/t_stem updates are per-patch own-index
# read-modify-write (no scatter). qsat()/smooth_* are kernel-safe; use_biomass_heat
# passed as a Bool; constants eltype-converted. (canopy is CPU-only — the Newton
# filter compaction blocks a whole-function GPU run — but this keeps every canopy
# loop in kernel form and GPU-ready.)
Base.@kwdef struct CfDiagOut{V}   # written per-patch outputs
    err_arr::V; dt_stem::V; dhsdt_canopy::V; t_stem::V; taux::V; tauy::V
    eflx_sh_grnd::V; eflx_sh_snow::V; eflx_sh_soil::V; eflx_sh_h2osfc::V
    qflx_evap_soi::V; qflx_ev_snow::V; qflx_ev_soil::V; qflx_ev_h2osfc::V
    t_ref2m::V; t_ref2m_r::V; q_ref2m::V; rh_ref2m::V; rh_ref2m_r::V; vpd_ref2m::V
    dlrad::V; ulrad::V; t_skin::V; cgrnds::V; cgrndl::V; cgrnd::V
    snocan_baseline::V; liqcan::V; snocan::V
end
Base.@kwdef struct CfDiagP{V}     # read-only per-patch inputs
    frac_rad_abs_by_stem::V; sabv::V; air::V; bir::V; tlbef::V; dt_veg::V; cir::V
    lw_leaf::V; lw_stem::V; eflx_sh_veg::V; qflx_evap_veg::V; tl_ini::V; cp_leaf::V
    cp_stem::V; ts_ini::V; eflx_sh_stem::V; stem_biomass::V; wtal::V; wtl0::V; wta0::V
    wtstem0::V; wtg::V; thm::V; ram1::V; wtgq::V; delq::V; wtalq::V; wtlq0::V; wtaq0::V
    qsatl::V; temp1::V; temp12m::V; dth::V; temp2::V; temp22m::V; dqh::V; emv::V
    qflx_tran_veg::V; t_veg::V
end
Base.@kwdef struct CfDiagC{V}     # read-only per-column inputs
    frac_sno_eff::V; t_h2osfc::V; frac_h2osfc::V; t_grnd::V; forc_rho::V
    qg_snow::V; qg_soil::V; qg_h2osfc::V; dqgdT::V; htvp::V; emg::V
    forc_q::V; forc_pbot::V; forc_lwrad::V
end
Adapt.@adapt_structure CfDiagOut
Adapt.@adapt_structure CfDiagP
Adapt.@adapt_structure CfDiagC

@kernel function _cf_postiter_kernel!(out, pp, cc, @Const(t_soisno), @Const(filterp),
        @Const(column), @Const(gridcell), @Const(snl), @Const(forc_u), @Const(forc_v),
        use_biomass::Bool, nlevsno::Int, dtime)
    fi = @index(Global)
    @inbounds begin
        p = filterp[fi]
        T = eltype(out.err_arr)
        c = column[p]; g = gridcell[p]
        snl_c = snl[c]
        rho = cc.forc_rho[c]; tveg = pp.t_veg[p]; thm = pp.thm[p]; frac_rad = pp.frac_rad_abs_by_stem[p]
        lw_grnd = cc.frac_sno_eff[c] * t_soisno[c, snl_c + 1 + nlevsno]^4 +
            (one(T) - cc.frac_sno_eff[c] - cc.frac_h2osfc[c]) * t_soisno[c, 1 + nlevsno]^4 +
            cc.frac_h2osfc[c] * cc.t_h2osfc[c]^4

        out.err_arr[p] = (one(T) - frac_rad) * (pp.sabv[p] + pp.air[p] +
            pp.bir[p] * pp.tlbef[p]^3 * (pp.tlbef[p] + T(4.0) * pp.dt_veg[p]) + pp.cir[p] * lw_grnd) -
            pp.lw_leaf[p] + pp.lw_stem[p] - pp.eflx_sh_veg[p] - T(HVAP) * pp.qflx_evap_veg[p] -
            ((tveg - pp.tl_ini[p]) * pp.cp_leaf[p] / dtime)

        # Stem temperature update
        if use_biomass
            if pp.stem_biomass[p] > zero(T)
                out.dt_stem[p] = (frac_rad * (pp.sabv[p] + pp.air[p] + pp.bir[p] * pp.ts_ini[p]^4 +
                    pp.cir[p] * lw_grnd) - pp.eflx_sh_stem[p] + pp.lw_leaf[p] - pp.lw_stem[p]) /
                    (pp.cp_stem[p] / dtime - frac_rad * pp.bir[p] * T(4.0) * pp.ts_ini[p]^3)
            else
                out.dt_stem[p] = zero(T)
            end
            out.dhsdt_canopy[p] = out.dt_stem[p] * pp.cp_stem[p] / dtime +
                (tveg - pp.tl_ini[p]) * pp.cp_leaf[p] / dtime
            out.t_stem[p] = out.t_stem[p] + out.dt_stem[p]
        else
            out.dt_stem[p] = zero(T)
        end
        tstem = out.t_stem[p]

        delt_val = pp.wtal[p] * cc.t_grnd[c] - pp.wtl0[p] * tveg - pp.wta0[p] * thm - pp.wtstem0[p] * tstem
        out.taux[p] = -rho * forc_u[g] / pp.ram1[p]
        out.tauy[p] = -rho * forc_v[g] / pp.ram1[p]
        out.eflx_sh_grnd[p] = T(CPAIR) * rho * pp.wtg[p] * delt_val

        delt_snow = pp.wtal[p] * t_soisno[c, snl_c + 1 + nlevsno] - pp.wtl0[p] * tveg - pp.wta0[p] * thm - pp.wtstem0[p] * tstem
        delt_soil = pp.wtal[p] * t_soisno[c, 1 + nlevsno] - pp.wtl0[p] * tveg - pp.wta0[p] * thm - pp.wtstem0[p] * tstem
        delt_h2osfc = pp.wtal[p] * cc.t_h2osfc[c] - pp.wtl0[p] * tveg - pp.wta0[p] * thm - pp.wtstem0[p] * tstem
        out.eflx_sh_snow[p] = T(CPAIR) * rho * pp.wtg[p] * delt_snow
        out.eflx_sh_soil[p] = T(CPAIR) * rho * pp.wtg[p] * delt_soil
        out.eflx_sh_h2osfc[p] = T(CPAIR) * rho * pp.wtg[p] * delt_h2osfc
        out.qflx_evap_soi[p] = rho * pp.wtgq[p] * pp.delq[p]

        delq_snow = pp.wtalq[p] * cc.qg_snow[c] - pp.wtlq0[p] * pp.qsatl[p] - pp.wtaq0[p] * cc.forc_q[c]
        out.qflx_ev_snow[p] = rho * pp.wtgq[p] * delq_snow
        delq_soil = pp.wtalq[p] * cc.qg_soil[c] - pp.wtlq0[p] * pp.qsatl[p] - pp.wtaq0[p] * cc.forc_q[c]
        out.qflx_ev_soil[p] = rho * pp.wtgq[p] * delq_soil
        delq_h2osfc = pp.wtalq[p] * cc.qg_h2osfc[c] - pp.wtlq0[p] * pp.qsatl[p] - pp.wtaq0[p] * cc.forc_q[c]
        out.qflx_ev_h2osfc[p] = rho * pp.wtgq[p] * delq_h2osfc

        out.t_ref2m[p] = thm + pp.temp1[p] * pp.dth[p] * (one(T) / pp.temp12m[p] - one(T) / pp.temp1[p])
        out.t_ref2m_r[p] = out.t_ref2m[p]
        out.q_ref2m[p] = cc.forc_q[c] + pp.temp2[p] * pp.dqh[p] * (one(T) / pp.temp22m[p] - one(T) / pp.temp2[p])
        (qsat_ref2m, e_ref2m, _, _) = qsat(out.t_ref2m[p], cc.forc_pbot[c])
        out.rh_ref2m[p] = smooth_min(T(100.0), out.q_ref2m[p] / qsat_ref2m * T(100.0))
        out.rh_ref2m_r[p] = out.rh_ref2m[p]
        out.vpd_ref2m[p] = e_ref2m * (one(T) - out.rh_ref2m[p] / T(100.0))

        emv_p = pp.emv[p]; emg_c = cc.emg[c]
        out.dlrad[p] = (one(T) - emv_p) * emg_c * cc.forc_lwrad[c] +
            emv_p * emg_c * T(SB) * pp.tlbef[p]^3 * (pp.tlbef[p] + T(4.0) * pp.dt_veg[p]) * (one(T) - frac_rad) +
            emv_p * emg_c * T(SB) * pp.ts_ini[p]^3 * (pp.ts_ini[p] + T(4.0) * out.dt_stem[p]) * frac_rad
        out.ulrad[p] = ((one(T) - emg_c) * (one(T) - emv_p) * (one(T) - emv_p) * cc.forc_lwrad[c] +
            emv_p * (one(T) + (one(T) - emg_c) * (one(T) - emv_p)) * T(SB) *
            pp.tlbef[p]^3 * (pp.tlbef[p] + T(4.0) * pp.dt_veg[p]) * (one(T) - frac_rad) +
            emv_p * (one(T) + (one(T) - emg_c) * (one(T) - emv_p)) * T(SB) *
            pp.ts_ini[p]^3 * (pp.ts_ini[p] + T(4.0) * out.dt_stem[p]) * frac_rad +
            emg_c * (one(T) - emv_p) * T(SB) * lw_grnd)
        out.t_skin[p] = emv_p * tveg + (one(T) - emv_p) * sqrt(sqrt(lw_grnd))

        out.cgrnds[p] = out.cgrnds[p] + T(CPAIR) * rho * pp.wtg[p] * pp.wtal[p]
        out.cgrndl[p] = out.cgrndl[p] + rho * pp.wtgq[p] * pp.wtalq[p] * cc.dqgdT[c]
        out.cgrnd[p] = out.cgrnds[p] + out.cgrndl[p] * cc.htvp[c]

        out.snocan_baseline[p] = out.snocan[p]
        _frac_liq = smooth_heaviside(tveg - T(TFRZ); k = T(200.0))
        _frac_ice = one(T) - _frac_liq
        _net_evap = (pp.qflx_tran_veg[p] - pp.qflx_evap_veg[p]) * dtime
        out.liqcan[p] = smooth_max(zero(T), out.liqcan[p] + _frac_liq * _net_evap)
        out.snocan[p] = smooth_max(zero(T), out.snocan[p] + _frac_ice * _net_evap)
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
# Thin keyword wrapper preserving the original interface (tests + other callers).
# Forwards to the all-positional canopy_fluxes_core! so the differentiated driver
# path (clm_drv_core!) can call it WITHOUT a kwarg NamedTuple — Enzyme reverse-mode
# cannot compile the Core.kwcall augmented thunk over this signature (see
# scripts/enzyme_fulldriver_probe.jl). The core lists the driver-supplied kwargs
# first so the driver call stays compact (trailing default-only args omitted).
function canopy_fluxes!(
        canopystate, energyflux, frictionvel, temperature, solarabs, soilstate,
        waterfluxbulk, waterstatebulk, waterdiagbulk, photosyns,
        patch_data, col_data, gridcell_data, mask_exposedvegp, bounds_patch, bounds_col,
        forc_lwrad_col, forc_q_col, forc_pbot_col, forc_th_col, forc_rho_col, forc_t_col,
        forc_u_grc, forc_v_grc, forc_pco2_grc, forc_po2_grc,
        forc_hgt_t_grc, forc_hgt_u_grc, forc_hgt_q_grc, dayl_grc, max_dayl_grc,
        downreg_patch, leafn_patch, dtime;
        dleaf_pft=fill(0.04, MXPFT+1), dbh_pft=fill(0.1, MXPFT+1),
        slatop_pft=fill(0.01, MXPFT+1), fbw_pft=fill(0.1, MXPFT+1),
        nstem_pft=fill(1.0, MXPFT+1), woody_pft=fill(0.0, MXPFT+1),
        rstem_per_dbh_pft=fill(0.0, MXPFT+1), wood_density_pft=fill(500.0, MXPFT+1),
        is_tree_pft=fill(false, MXPFT+1), is_shrub_pft=fill(false, MXPFT+1),
        c3psn_pft=fill(1.0, MXPFT+1), leafcn_pft=fill(25.0, MXPFT+1),
        flnr_pft=fill(0.08, MXPFT+1), fnitr_pft=fill(0.1, MXPFT+1),
        mbbopt_pft=fill(0.0, MXPFT+1), medlynintercept_pft=fill(100.0, MXPFT+1),
        medlynslope_pft=fill(6.0, MXPFT+1), crop_pft=Float64[],
        z0v_Cr_pft=fill(0.35, MXPFT+1), z0v_Cs_pft=fill(0.003, MXPFT+1),
        z0v_c_pft=fill(0.25, MXPFT+1), z0v_cw_pft=fill(2.0, MXPFT+1),
        z0v_LAImax_pft=fill(8.0, MXPFT+1),
        use_cn=false, use_lch4=false, use_c13=false, use_hydrstress=false,
        use_fates=false, use_luna=false, z0param_method="ZengWang2007",
        o3coefv_patch=Float64[], o3coefg_patch=Float64[],
        grnd_ch4_cond_patch=Float64[], forc_pc13o2_grc=Float64[],
        t10_patch=Float64[], nrad_patch=Int[],
        tlai_z_patch=Matrix{Float64}(undef,0,0),
        vcmaxcint_sun_patch=Float64[], vcmaxcint_sha_patch=Float64[],
        parsun_z_patch=Matrix{Float64}(undef,0,0), parsha_z_patch=Matrix{Float64}(undef,0,0),
        laisun_z_patch=Matrix{Float64}(undef,0,0), laisha_z_patch=Matrix{Float64}(undef,0,0),
        phs_froot_carbon=Float64[],
        leaf_mr_vcm=0.015, overrides=CalibrationOverrides())
    return canopy_fluxes_core!(
        canopystate, energyflux, frictionvel, temperature, solarabs, soilstate,
        waterfluxbulk, waterstatebulk, waterdiagbulk, photosyns,
        patch_data, col_data, gridcell_data, mask_exposedvegp, bounds_patch, bounds_col,
        forc_lwrad_col, forc_q_col, forc_pbot_col, forc_th_col, forc_rho_col, forc_t_col,
        forc_u_grc, forc_v_grc, forc_pco2_grc, forc_po2_grc,
        forc_hgt_t_grc, forc_hgt_u_grc, forc_hgt_q_grc, dayl_grc, max_dayl_grc,
        downreg_patch, leafn_patch, dtime,
        t10_patch, nrad_patch, tlai_z_patch, vcmaxcint_sun_patch, vcmaxcint_sha_patch,
        parsun_z_patch, parsha_z_patch, laisun_z_patch, laisha_z_patch,
        o3coefv_patch, o3coefg_patch, dleaf_pft, slatop_pft, leafcn_pft, flnr_pft,
        fnitr_pft, mbbopt_pft, c3psn_pft, woody_pft, overrides,
        dbh_pft, fbw_pft, nstem_pft, rstem_per_dbh_pft, wood_density_pft,
        is_tree_pft, is_shrub_pft, medlynintercept_pft, medlynslope_pft, crop_pft,
        z0v_Cr_pft, z0v_Cs_pft, z0v_c_pft, z0v_cw_pft, z0v_LAImax_pft,
        use_cn, use_lch4, use_c13, use_hydrstress, phs_froot_carbon, use_fates, use_luna,
        z0param_method, grnd_ch4_cond_patch, forc_pc13o2_grc, leaf_mr_vcm)
end

# All-positional core (body lives here). Param groups: (1) state/forcing, then
# (2) driver-supplied "kwargs" (front), then (3) default-only params (back).
function canopy_fluxes_core!(
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
        mask_exposedvegp::AbstractVector{Bool},
        bounds_patch    ::UnitRange{Int},
        bounds_col      ::UnitRange{Int},
        forc_lwrad_col  ::AbstractVector{<:Real},
        forc_q_col      ::AbstractVector{<:Real},
        forc_pbot_col   ::AbstractVector{<:Real},
        forc_th_col     ::AbstractVector{<:Real},
        forc_rho_col    ::AbstractVector{<:Real},
        forc_t_col      ::AbstractVector{<:Real},
        forc_u_grc      ::AbstractVector{<:Real},
        forc_v_grc      ::AbstractVector{<:Real},
        forc_pco2_grc   ::AbstractVector{<:Real},
        forc_po2_grc    ::AbstractVector{<:Real},
        forc_hgt_t_grc  ::AbstractVector{<:Real},
        forc_hgt_u_grc  ::AbstractVector{<:Real},
        forc_hgt_q_grc  ::AbstractVector{<:Real},
        dayl_grc        ::AbstractVector{<:Real},
        max_dayl_grc    ::AbstractVector{<:Real},
        downreg_patch   ::AbstractVector{<:Real},
        leafn_patch     ::AbstractVector{<:Real},
        dtime           ::Real,
        # (2) driver-supplied group
        t10_patch       ::AbstractVector{<:Real} = Float64[],
        nrad_patch      ::AbstractVector{<:Integer} = Int[],
        tlai_z_patch    ::AbstractMatrix{<:Real} = Matrix{Float64}(undef,0,0),
        vcmaxcint_sun_patch ::AbstractVector{<:Real} = Float64[],
        vcmaxcint_sha_patch ::AbstractVector{<:Real} = Float64[],
        parsun_z_patch  ::AbstractMatrix{<:Real} = Matrix{Float64}(undef,0,0),
        parsha_z_patch  ::AbstractMatrix{<:Real} = Matrix{Float64}(undef,0,0),
        laisun_z_patch  ::AbstractMatrix{<:Real} = Matrix{Float64}(undef,0,0),
        laisha_z_patch  ::AbstractMatrix{<:Real} = Matrix{Float64}(undef,0,0),
        o3coefv_patch   ::AbstractVector{<:Real} = Float64[],
        o3coefg_patch   ::AbstractVector{<:Real} = Float64[],
        dleaf_pft       ::AbstractVector{<:Real} = fill(0.04, MXPFT+1),
        slatop_pft      ::AbstractVector{<:Real} = fill(0.01, MXPFT+1),
        leafcn_pft      ::AbstractVector{<:Real} = fill(25.0, MXPFT+1),
        flnr_pft        ::AbstractVector{<:Real} = fill(0.08, MXPFT+1),
        fnitr_pft       ::AbstractVector{<:Real} = fill(0.1, MXPFT+1),
        mbbopt_pft      ::AbstractVector{<:Real} = fill(0.0, MXPFT+1),
        c3psn_pft       ::AbstractVector{<:Real} = fill(1.0, MXPFT+1),
        woody_pft       ::AbstractVector{<:Real} = fill(0.0, MXPFT+1),
        overrides       ::CalibrationOverrides = CalibrationOverrides(),
        # (3) default-only group
        dbh_pft         ::AbstractVector{<:Real} = fill(0.1, MXPFT+1),
        fbw_pft         ::AbstractVector{<:Real} = fill(0.1, MXPFT+1),
        nstem_pft       ::AbstractVector{<:Real} = fill(1.0, MXPFT+1),
        rstem_per_dbh_pft::AbstractVector{<:Real} = fill(0.0, MXPFT+1),
        wood_density_pft::AbstractVector{<:Real} = fill(500.0, MXPFT+1),
        is_tree_pft     ::AbstractVector{Bool} = fill(false, MXPFT+1),
        is_shrub_pft    ::AbstractVector{Bool} = fill(false, MXPFT+1),
        medlynintercept_pft::AbstractVector{<:Real} = fill(100.0, MXPFT+1),
        medlynslope_pft ::AbstractVector{<:Real} = fill(6.0, MXPFT+1),
        crop_pft        ::AbstractVector{<:Real} = Float64[],
        z0v_Cr_pft      ::AbstractVector{<:Real} = fill(0.35, MXPFT+1),
        z0v_Cs_pft      ::AbstractVector{<:Real} = fill(0.003, MXPFT+1),
        z0v_c_pft       ::AbstractVector{<:Real} = fill(0.25, MXPFT+1),
        z0v_cw_pft      ::AbstractVector{<:Real} = fill(2.0, MXPFT+1),
        z0v_LAImax_pft  ::AbstractVector{<:Real} = fill(8.0, MXPFT+1),
        use_cn          ::Bool = false,
        use_lch4        ::Bool = false,
        use_c13         ::Bool = false,
        use_hydrstress  ::Bool = false,
        phs_froot_carbon::AbstractVector{<:Real} = Float64[],  # PHS: CN frootc (gC/m2)
        use_fates       ::Bool = false,
        use_luna        ::Bool = false,
        z0param_method  ::String = "ZengWang2007",
        grnd_ch4_cond_patch::AbstractVector{<:Real} = Float64[],
        forc_pc13o2_grc ::AbstractVector{<:Real} = Float64[],
        leaf_mr_vcm     ::Real = 0.015)

    np = length(bounds_patch)
    begp = first(bounds_patch)
    endp = last(bounds_patch)
    nlevgrnd = varpar.nlevgrnd
    nlevsno  = varpar.nlevsno

    # --- Convenience aliases ---
    params = canopy_fluxes_params
    ctrl   = canopy_fluxes_ctrl

    # Local work arrays (patch-indexed). Allocated via similar() of an input
    # state array (NOT zeros()) + fill!, so they share the state's backend —
    # device-resident on GPU. A bare zeros() forces a host Array and would break
    # a whole-function GPU run. `_sc`/`_sci` are the float/Int scratch factories.
    FT = eltype(forc_t_col)
    _sc(dims...)  = fill!(similar(temperature.t_veg_patch, FT,  dims...), zero(FT))
    _sci(dims...) = fill!(similar(patch_data.column,       Int, dims...), 0)
    zldis     = _sc(endp)
    dth       = _sc(endp)
    dthv      = _sc(endp)
    dqh       = _sc(endp)
    ur        = _sc(endp)
    temp1     = _sc(endp)
    temp12m   = _sc(endp)
    temp2     = _sc(endp)
    temp22m   = _sc(endp)
    rb        = _sc(endp)
    rah       = _sc(endp, 2)
    raw       = _sc(endp, 2)
    wtg       = _sc(endp)
    wta0      = _sc(endp)
    wtl0      = _sc(endp)
    wtstem0   = _sc(endp)
    wtal      = _sc(endp)
    wtga      = _sc(endp)
    wtgq      = _sc(endp)
    wtaq0     = _sc(endp)
    wtlq0     = _sc(endp)
    wtalq     = _sc(endp)
    el        = _sc(endp)
    qsatl     = _sc(endp)
    qsatldT   = _sc(endp)
    air       = _sc(endp)
    bir       = _sc(endp)
    cir       = _sc(endp)
    del_arr   = _sc(endp)  # "del" in Fortran (reserved word in Julia context for clarity)
    del2      = _sc(endp)
    dele      = _sc(endp)
    delq      = _sc(endp)
    det_arr   = _sc(endp)
    efeb      = _sc(endp)
    efe       = _sc(endp)
    err_arr   = _sc(endp)
    obuold    = _sc(endp)
    tlbef     = _sc(endp)
    tl_ini    = _sc(endp)
    ts_ini    = _sc(endp)
    co2_arr   = _sc(endp)
    o2_arr    = _sc(endp)
    svpts     = _sc(endp)
    eah       = _sc(endp)
    dt_veg    = _sc(endp)
    fm        = _sc(endp)
    nmozsgn   = _sci(endp)
    dayl_factor = _sc(endp)
    rootsum   = _sc(endp)
    dbh       = _sc(endp)
    cp_leaf   = _sc(endp)
    cp_stem   = _sc(endp)
    rstem     = _sc(endp)
    dt_stem   = _sc(endp)
    frac_rad_abs_by_stem = _sc(endp)
    lw_stem   = _sc(endp)
    lw_leaf   = _sc(endp)
    sa_stem   = _sc(endp)
    sa_leaf   = _sc(endp)
    sa_internal = _sc(endp)
    uuc       = _sc(endp)
    snocan_baseline = _sc(endp)

    # Integer filter array (Fortran-style). Stays fixed through the Newton
    # iteration — convergence is tracked by the `active` mask, not by compaction.
    # Built on the host from the (host-only) exposed-veg mask, then copied to the
    # device backend in bulk (copyto! is a transfer, not device scalar indexing);
    # the kernels only ever read filterp. filterp_host is kept for the few
    # remaining host-side diagnostic scans.
    filterp = _sci(np)
    filterp_host = zeros(Int, np)
    fn = 0
    # Host-scan the mask (Array() it once so a device mask isn't scalar-indexed here
    # or in the later host-side diagnostic scans).
    mask_ev_host = mask_exposedvegp isa Array ? mask_exposedvegp : Array(mask_exposedvegp)
    for p in bounds_patch
        if mask_ev_host[p]
            fn += 1
            filterp_host[fn] = p
        end
    end
    copyto!(filterp, filterp_host)

    # Move host pft-parameter kwargs + the ch4 placeholder onto the working backend.
    # The driver passes some pft params via device copies already; the rest default to
    # host fill() arrays that get bundled into device-view structs / passed to kernels.
    # No-op on CPU (and skips args already on-device). Length-preserving.
    if !(canopystate.elai_patch isa Array)
        _cfref = canopystate.elai_patch
        Tcf = eltype(_cfref)
        _pf(a) = a isa Array ? copyto!(similar(_cfref, Tcf, length(a)), Tcf.(a)) : a
        _pb(a) = a isa Array ? copyto!(similar(_cfref, Bool, length(a)), a) : a
        grnd_ch4_cond_patch = _pf(grnd_ch4_cond_patch)
        medlynintercept_pft = _pf(medlynintercept_pft); medlynslope_pft = _pf(medlynslope_pft)
        crop_pft = _pf(crop_pft); dbh_pft = _pf(dbh_pft); fbw_pft = _pf(fbw_pft)
        nstem_pft = _pf(nstem_pft); rstem_per_dbh_pft = _pf(rstem_per_dbh_pft)
        wood_density_pft = _pf(wood_density_pft)
        is_tree_pft = _pb(is_tree_pft); is_shrub_pft = _pb(is_shrub_pft)
    end

    # Time step initialization of photosynthesis variables (kernelized over the
    # filtered patches; canopy uses the base path — use_c13/use_c14 = false).
    cf_psn_init_update!(photosyns, filterp, fn)

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

    # --- PHS (use_hydrstress): soil matric potential + soil-to-root conductance ---
    # Normally smp_l is filled by HydrologyNoDrainage (which runs after canopy_fluxes),
    # so recompute it here from the injected soil water; then call PHS Pass 1 to fill
    # k_soil_root, which the PHS photosynthesis (calcstress!/vegwp solve) needs.
    if use_hydrstress
        _nlevsoi = varpar.nlevsoi; _nlevsno = varpar.nlevsno
        _smp_l = soilstate.smp_l_col; _watsat = soilstate.watsat_col
        _sucsat = soilstate.sucsat_col; _bsw = soilstate.bsw_col
        _smpmin = soilstate.smpmin_col
        # s_node from the LIQUID volumetric water (driver fills h2osoi_liqvol_col from
        # the prognostic h2osoi_liq before canopy_fluxes). Do NOT use h2osoi_vol_col:
        # it is only updated by HydrologyNoDrainage (which runs AFTER canopy_fluxes),
        # so here it sits at its cold-start default (≈0.75·watsat) → spuriously wet
        # soil → no plant water stress (grass vegwp stuck ~−9000 mm vs Fortran ~−205000).
        _h2osoi_liqvol = waterdiagbulk.h2osoi_liqvol_col
        _hksat = soilstate.hksat_col; _hk_l = soilstate.hk_l_col
        @inbounds for c in bounds_col, j in 1:_nlevsoi
            s_node = max(min(_h2osoi_liqvol[c, _nlevsno + j] / _watsat[c, j], one(FT)), FT(0.01))
            _smp_l[c, j] = max(_smpmin[c], -_sucsat[c, j] * s_node^(-_bsw[c, j]))
            # Clapp-Hornberger unsaturated conductivity (also filled by HydrologyNoDrainage,
            # which runs after canopy_fluxes, so recompute here for PHS).
            _hk_l[c, j] = _hksat[c, j] * s_node^(FT(2.0) * _bsw[c, j] + FT(3.0))
        end
        _froot_c = isempty(phs_froot_carbon) ? zeros(FT, length(bounds_patch)) : phs_froot_carbon
        psn_phs_pass1_update!(photosyns, mask_exposedvegp, patch_data.column,
            patch_data.itype .+ 1, _froot_c, soilstate.rootfr_patch,
            view(col_data.dz, :, _nlevsno+1:_nlevsno+_nlevsoi),
            canopystate.tsai_patch, canopystate.tlai_patch,
            pftcon.froot_leaf, pftcon.root_radius, pftcon.root_density,
            soilstate.hksat_col, soilstate.hk_l_col, _smp_l,
            view(col_data.z, :, _nlevsno+1:_nlevsno+_nlevsoi),
            soilstate.k_soil_root_patch, soilstate.root_conductance_patch,
            soilstate.soil_conductance_patch, C_TO_B, 0.25, _nlevsoi, bounds_patch)
    end

    # --- Modify aerodynamic parameters for sparse/dense canopy (X. Zeng) ---
    # Resolve the method String to an Int code on the host (so no String compare or
    # error() runs in the kernel); validate here.
    z0_method = z0param_method == "ZengWang2007" ? 1 :
                z0param_method == "Meier2022" ? 2 :
                error("canopy_fluxes!: unknown z0param_method: $z0param_method")
    # z0v_* default to host pftcon-like fill() arrays; move onto the working backend
    # so the kernel doesn't get host Arrays among device args (no-op on CPU).
    if !(canopystate.displa_patch isa Array)
        Tz0 = eltype(canopystate.displa_patch)
        _z0v_dev(a) = copyto!(similar(canopystate.displa_patch, Tz0, length(a)), Tz0.(a))
        z0v_LAImax_pft = _z0v_dev(z0v_LAImax_pft); z0v_Cs_pft = _z0v_dev(z0v_Cs_pft)
        z0v_Cr_pft = _z0v_dev(z0v_Cr_pft); z0v_c_pft = _z0v_dev(z0v_c_pft)
        z0v_cw_pft = _z0v_dev(z0v_cw_pft)
    end
    _launch!(_cf_z0_kernel!, canopystate.displa_patch, frictionvel.z0mv_patch,
        frictionvel.z0hv_patch, frictionvel.z0qv_patch, frictionvel.forc_hgt_u_patch,
        frictionvel.forc_hgt_t_patch, frictionvel.forc_hgt_q_patch, filterp,
        patch_data.column, patch_data.gridcell, patch_data.itype, canopystate.elai_patch,
        canopystate.esai_patch, canopystate.htop_patch, frictionvel.z0mg_col,
        forc_hgt_u_grc, forc_hgt_t_grc, forc_hgt_q_grc, z0v_LAImax_pft, z0v_Cs_pft,
        z0v_Cr_pft, z0v_c_pft, z0v_cw_pft, z0_method; ndrange = fn)

    # --- Net absorbed longwave radiation, QSat, CO2/O2, flux initialization ---
    # Arrays grouped into device-view structs (Metal 31-buffer limit); struct-arg
    # kernel launched via manual backend + KA.synchronize.
    lwq_out = CfLwqOut(; air = air, bir = bir, cir = cir, qsatl = qsatl, el = el,
        qsatldT = qsatldT, co2_arr = co2_arr, o2_arr = o2_arr,
        taf = frictionvel.taf_patch, qaf = frictionvel.qaf_patch, ur = ur,
        dth = dth, dqh = dqh, delq = delq, dthv = dthv, zldis = zldis)
    lwq_p = CfLwqP(; emv = temperature.emv_patch, t_veg = temperature.t_veg_patch,
        thm = temperature.thm_patch, forc_hgt_u = frictionvel.forc_hgt_u_patch,
        displa = canopystate.displa_patch)
    lwq_c = CfLwqC(; emg = temperature.emg_col, forc_lwrad = forc_lwrad_col,
        forc_pbot = forc_pbot_col, t_grnd = temperature.t_grnd_col, forc_q = forc_q_col,
        qg = waterdiagbulk.qg_col, forc_th = forc_th_col)
    lwq_g = CfLwqG(; forc_pco2 = forc_pco2_grc, forc_po2 = forc_po2_grc,
        forc_u = forc_u_grc, forc_v = forc_v_grc)
    let be = _kernel_backend(air)
        _cf_longwave_qsat_kernel!(be)(lwq_out, lwq_p, lwq_c, lwq_g, nmozsgn,
            filterp, patch_data.column, patch_data.gridcell,
            convert(eltype(air), params.wind_min); ndrange = fn)
        KA.synchronize(be)
    end

    # Forcing-height-below-canopy warning (host-side diagnostic; matches the scalar
    # version — last patch with zldis<0 reported once). Pull zldis to the host
    # once (bulk copy, not device scalar indexing) so this stays GPU-safe.
    zldis_host = Array(zldis)
    found = false; found_index = 0
    for fi in 1:fn
        p = filterp_host[fi]
        if zldis_host[p] < 0.0
            found = true; found_index = p
        end
    end
    if found && !use_fates
        @warn "canopy_fluxes!: forcing height below canopy height at patch $found_index"
    end

    # --- Initialize Monin-Obukhov length (kernelized) ---
    cf_moninobukini_update!(frictionvel, temperature, patch_data, filterp, fn,
        ur, dthv, zldis, tl_ini, ts_ini)

    # --- Save original filter size for post-iteration ---
    itlef = 0
    fnorig = fn

    # Per-patch active mask for the Newton iteration. The filter stays fixed at
    # the original `fnorig` patches (ndrange = fnorig every iteration); converged
    # patches are skipped in-kernel via `active[p]` instead of host-side filter
    # compaction, so no device→host sync is needed mid-iteration.
    active_host = fill(false, endp)
    for fi in 1:fn
        active_host[filterp_host[fi]] = true
    end
    active = fill!(similar(patch_data.active, Bool, endp), false)
    copyto!(active, active_host)

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

    while itlef <= ctrl.itmax_canopy_fluxes && any(active)

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
            temp1, temp2, temp12m, temp22m, fm; active = active)

        # First per-patch loop: aerodynamic resistances / conductances.
        # Kernelized (one thread per filtered patch); csoilc_val and feature
        # flags resolved on the host and passed as scalars.
        cf_resist_update!(frictionvel, canopystate, temperature, waterdiagbulk,
            patch_data, filterp, fn, active,
            temp1, temp2, tlbef, del2, del_arr, rah, raw, uuc, rb, svpts, eah, el,
            dleaf_pft, grnd_ch4_cond_patch, forc_pbot_col,
            csoilc_val, ctrl.use_undercanopy_stability,
            ctrl.use_biomass_heat_storage, use_lch4, params)

        # --- Photosynthesis ---
        # Call photosynthesis for sunlit and shaded leaves to update rssun/rssha
        if !isempty(nrad_patch) && !isempty(parsun_z_patch)
            # Build patch-indexed forc_pbot from column-level data
            # Use full mask (not reduced filterp) since photosynthesis uses mask_exposedvegp
            # Gather on the host (Array() the device inputs once) then bulk-copy to
            # a device-resident array — avoids device scalar indexing.
            fpp_host = zeros(FT, endp)
            col_host = patch_data.column isa Array ? patch_data.column : Array(patch_data.column)
            fpc_host = forc_pbot_col isa Array ? forc_pbot_col : Array(forc_pbot_col)
            for p in bounds_patch
                mask_ev_host[p] || continue
                fpp_host[p] = fpc_host[col_host[p]]
            end
            forc_pbot_patch = forc_pbot_col isa Array ? fpp_host :
                              copyto!(similar(forc_pbot_col, FT, endp), fpp_host)

            # PFT index vector (+1 for 0-based Fortran → 1-based Julia)
            ivt_vec = patch_data.itype .+ 1

          if use_hydrstress
            # --- PHS photosynthesis (plant hydraulic stress) ---
            # forc_rho is COLUMN-indexed inside the PHS solver (forc_rho[c]); pass
            # the column array, NOT a per-patch copy. A patch-indexed array would be
            # read at the wrong slot (c = column[p]) — for multi-patch columns the
            # tree/grass patch grabs another patch's slot (zero for bare ground),
            # zeroing efpot and collapsing the transpiration demand.
            _nsno = varpar.nlevsno; _nsoi = varpar.nlevsoi
            _froot_c2 = isempty(phs_froot_carbon) ? zeros(FT, endp) : phs_froot_carbon
            photosynthesis_hydrstress!(photosyns,
                svpts, eah, o2_arr, co2_arr, rb,
                energyflux.btran_patch, dayl_factor, leafn_patch, qsatl, frictionvel.qaf_patch,
                forc_pbot_patch, forc_rho_col, temperature.t_veg_patch, t10_patch, temperature.thm_patch,
                nrad_patch, tlai_z_patch, canopystate.tlai_patch, canopystate.tsai_patch,
                parsun_z_patch, parsha_z_patch, laisun_z_patch, laisha_z_patch,
                vcmaxcint_sun_patch, vcmaxcint_sha_patch,
                o3coefv_patch, o3coefg_patch, o3coefv_patch, o3coefg_patch,
                c3psn_pft, leafcn_pft, flnr_pft, fnitr_pft, slatop_pft,
                mbbopt_pft, medlynintercept_pft, medlynslope_pft,
                pftcon.froot_leaf, pftcon.root_radius, pftcon.root_density,
                crop_pft, ivt_vec, patch_data.column, mask_exposedvegp, bounds_patch,
                _froot_c2, _froot_c2, soilstate.k_soil_root_patch,
                soilstate.root_conductance_patch, soilstate.soil_conductance_patch,
                soilstate.rootfr_patch,
                col_data.dz[:, _nsno+1:_nsno+_nsoi], col_data.z[:, _nsno+1:_nsno+_nsoi],
                soilstate.hk_l_col, soilstate.hksat_col, soilstate.smp_l_col,
                canopystate.vegwp_patch, canopystate.vegwp_ln_patch,
                canopystate.laisun_patch, canopystate.laisha_patch,
                canopystate.elai_patch, canopystate.esai_patch,
                canopystate.htop_patch, waterdiagbulk.fdry_patch, waterfluxbulk.wf.qflx_tran_veg_patch;
                nlevcan=NLEVCAN, nlevsoi=_nsoi, use_cn=use_cn)
          else
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
          end  # use_hydrstress

        end

        # --- Heat transfer conductances + Newton update of leaf temperature ---
        # Second per-patch loop: kernelized (one thread per filtered patch).
        # The two control-flag functions (do_soilevap_beta /
        # do_soil_resistance_sl14), the methane conductance length guard, and
        # the feature flags are resolved on the host and passed as scalars.
        cf_energy_update!(canopystate, energyflux, frictionvel, temperature,
            solarabs, soilstate, waterfluxbulk, waterstatebulk, waterdiagbulk,
            photosyns, patch_data, col_data, filterp, fn, active,
            rah, raw, rb, rstem, sa_leaf, sa_stem, sa_internal,
            frac_rad_abs_by_stem, air, bir, cir, cp_leaf, tl_ini, tlbef, zldis,
            temp1, temp2, ur, efeb, wtg, wtl0, wta0, wtstem0, wtga, wtal,
            lw_stem, lw_leaf, wtgq, wtlq0, wtaq0, wtalq, efe, dt_veg, del_arr,
            err_arr, qsatl, el, qsatldT, dth, dqh, delq, obuold, nmozsgn,
            forc_q_col, forc_th_col, forc_pbot_col, forc_rho_col,
            soilevap_beta, soil_resis_sl14, use_lch4, use_hydrstress,
            nlevsno, dtime, params)

        # --- Test for convergence ---
        # Kernelized: updates the convergence metrics for still-active patches
        # and clears active[p] on convergence. The filter stays fixed at the
        # original patch set (fn == fnorig throughout); converged patches are
        # skipped on subsequent iterations via the in-kernel active[p] guard
        # instead of host-side filter compaction, so the physics arrays are
        # never copied back to the host mid-iteration.
        itlef += 1
        if itlef > ITMIN_CANOPY
            cf_converge_update!(active, efeb, dele, det_arr,
                frictionvel.num_iter_patch, filterp, fn, efe, del_arr, del2, itlef)
        end
    end  # End stability iteration

    # =========================================================================
    # Phase 3: Post-iteration diagnostics
    # =========================================================================

    # The filter stayed fixed through the iteration (mask-based convergence),
    # so filterp already holds the original patch set; just restore the count.
    fn = fnorig

    diag_out = CfDiagOut(; err_arr = err_arr, dt_stem = dt_stem,
        dhsdt_canopy = energyflux.dhsdt_canopy_patch, t_stem = temperature.t_stem_patch,
        taux = energyflux.taux_patch, tauy = energyflux.tauy_patch,
        eflx_sh_grnd = energyflux.eflx_sh_grnd_patch, eflx_sh_snow = energyflux.eflx_sh_snow_patch,
        eflx_sh_soil = energyflux.eflx_sh_soil_patch, eflx_sh_h2osfc = energyflux.eflx_sh_h2osfc_patch,
        qflx_evap_soi = waterfluxbulk.wf.qflx_evap_soi_patch,
        qflx_ev_snow = waterfluxbulk.qflx_ev_snow_patch, qflx_ev_soil = waterfluxbulk.qflx_ev_soil_patch,
        qflx_ev_h2osfc = waterfluxbulk.qflx_ev_h2osfc_patch, t_ref2m = temperature.t_ref2m_patch,
        t_ref2m_r = temperature.t_ref2m_r_patch, q_ref2m = waterdiagbulk.q_ref2m_patch,
        rh_ref2m = waterdiagbulk.rh_ref2m_patch, rh_ref2m_r = waterdiagbulk.rh_ref2m_r_patch,
        vpd_ref2m = waterdiagbulk.vpd_ref2m_patch, dlrad = energyflux.dlrad_patch,
        ulrad = energyflux.ulrad_patch, t_skin = temperature.t_skin_patch,
        cgrnds = energyflux.cgrnds_patch, cgrndl = energyflux.cgrndl_patch, cgrnd = energyflux.cgrnd_patch,
        snocan_baseline = snocan_baseline, liqcan = waterstatebulk.ws.liqcan_patch,
        snocan = waterstatebulk.ws.snocan_patch)
    diag_p = CfDiagP(; frac_rad_abs_by_stem = frac_rad_abs_by_stem, sabv = solarabs.sabv_patch,
        air = air, bir = bir, tlbef = tlbef, dt_veg = dt_veg, cir = cir, lw_leaf = lw_leaf,
        lw_stem = lw_stem, eflx_sh_veg = energyflux.eflx_sh_veg_patch,
        qflx_evap_veg = waterfluxbulk.wf.qflx_evap_veg_patch, tl_ini = tl_ini, cp_leaf = cp_leaf,
        cp_stem = cp_stem, ts_ini = ts_ini, eflx_sh_stem = energyflux.eflx_sh_stem_patch,
        stem_biomass = canopystate.stem_biomass_patch, wtal = wtal, wtl0 = wtl0, wta0 = wta0,
        wtstem0 = wtstem0, wtg = wtg, thm = temperature.thm_patch, ram1 = frictionvel.ram1_patch,
        wtgq = wtgq, delq = delq, wtalq = wtalq, wtlq0 = wtlq0, wtaq0 = wtaq0, qsatl = qsatl,
        temp1 = temp1, temp12m = temp12m, dth = dth, temp2 = temp2, temp22m = temp22m, dqh = dqh,
        emv = temperature.emv_patch, qflx_tran_veg = waterfluxbulk.wf.qflx_tran_veg_patch,
        t_veg = temperature.t_veg_patch)
    diag_c = CfDiagC(; frac_sno_eff = waterdiagbulk.frac_sno_eff_col,
        t_h2osfc = temperature.t_h2osfc_col, frac_h2osfc = waterdiagbulk.frac_h2osfc_col,
        t_grnd = temperature.t_grnd_col, forc_rho = forc_rho_col,
        qg_snow = waterdiagbulk.qg_snow_col, qg_soil = waterdiagbulk.qg_soil_col,
        qg_h2osfc = waterdiagbulk.qg_h2osfc_col, dqgdT = waterdiagbulk.dqgdT_col,
        htvp = energyflux.htvp_col, emg = temperature.emg_col, forc_q = forc_q_col,
        forc_pbot = forc_pbot_col, forc_lwrad = forc_lwrad_col)
    if fn > 0
        be = _kernel_backend(err_arr)
        _cf_postiter_kernel!(be)(diag_out, diag_p, diag_c, temperature.t_soisno_col, filterp,
            patch_data.column, patch_data.gridcell, col_data.snl, forc_u_grc, forc_v_grc,
            ctrl.use_biomass_heat_storage, nlevsno, dtime; ndrange = fn)
        KA.synchronize(be)
    end

    # --- Post-photosynthesis: PhotosynthesisTotal (kernelized; base path) ---
    cf_psn_total_update!(photosyns, canopystate.laisun_patch,
        canopystate.laisha_patch, filterp, fn)

    # Water use efficiency (iwue) — stub for local noon check
    # LUNA and ozone updates would go here in full implementation

    # Report patches with large energy-balance errors (host-side diagnostic).
    # Post-iteration only — pull err_arr to the host once (bulk copy) and scan
    # against the host filter; no device scalar indexing, no device filter mutate.
    err_host = Array(err_arr)
    for fi in 1:fn
        p = filterp_host[fi]
        if abs(err_host[p]) > 0.1
            @warn "canopy_fluxes!: energy balance error at patch $p, err=$(err_host[p])"
        end
    end

    return nothing
end
