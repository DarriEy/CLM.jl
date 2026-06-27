# ==========================================================================
# Compositional reverse-mode AD for canopy_fluxes_core!.
#
# canopy_fluxes_core! cannot be reverse-differentiated as one Enzyme thunk — the
# whole function trips an Enzyme-internal `sret` assertion (propagate_returned!).
# Every constituent sub-phase differentiates cleanly in isolation, though, so this
# module differentiates canopy COMPOSITIONALLY: it runs the SAME sub-phases
# (friction_velocity! → cf_resist_update! → photosynthesis!(sun,sha) →
# cf_energy_update!, after the init kernels) on a checkpointed state bundle, and
# back-propagates one SEPARATE Enzyme.autodiff call per sub-phase. Each call
# compiles its own small thunk and dodges the monolith's sret bug.
#
# The decomposed forward (`canopy_rev_forward!`) is byte-identical to
# canopy_fluxes_core!'s energy-balance/photosynthesis path on a converged single
# patch (guarded by test/test_canopy_reverse.jl); so the reverse gradient is the
# real canopy gradient. Two gradient-preserving simplifications vs the monolith:
#   - post-iteration diagnostics are dropped (they don't feed the iterated state);
#   - the data-dependent Newton `while` loop is unrolled to a FIXED N iterations
#     with active≡true (differentiate-through-the-converged-iterate); once the
#     solve has converged the fixed point doesn't move, so primal + gradient match.
#
# Scope: single converged patch / energy-balance + photosynthesis. Multi-patch with
# per-patch convergence and threading into the clm_drv! compositional reverse are
# follow-ups (see the [[enzyme-reverse-ad-julia110]] memory).
# ==========================================================================

# Scratch work arrays (one per canopy_fluxes_core! `_sc(endp)` local). Held in a flat
# NamedTuple so the whole canopy state (structs + scratch) can be checkpointed via
# deepcopy and shadowed via Enzyme.make_zero.
function cf_rev_scratch(::Type{FT}, np::Int) where {FT}
    z(d...) = zeros(FT, d...)
    return (; zldis=z(np), dth=z(np), dthv=z(np), dqh=z(np), ur=z(np), temp1=z(np),
        temp12m=z(np), temp2=z(np), temp22m=z(np), rb=z(np), rah=z(np,2), raw=z(np,2),
        wtg=z(np), wta0=z(np), wtl0=z(np), wtstem0=z(np), wtal=z(np), wtga=z(np),
        wtgq=z(np), wtaq0=z(np), wtlq0=z(np), wtalq=z(np), el=z(np), qsatl=z(np),
        qsatldT=z(np), air=z(np), bir=z(np), cir=z(np), del_arr=z(np), del2=z(np),
        dele=z(np), delq=z(np), det_arr=z(np), efeb=z(np), efe=z(np), err_arr=z(np),
        obuold=z(np), tlbef=z(np), tl_ini=z(np), ts_ini=z(np), co2_arr=z(np), o2_arr=z(np),
        svpts=z(np), eah=z(np), dt_veg=z(np), fm=z(np), nmozsgn=zeros(Int, np),
        dayl_factor=z(np), cp_leaf=z(np), rstem=z(np), frac_rad_abs_by_stem=z(np),
        sa_stem=z(np), sa_leaf=z(np), sa_internal=z(np), uuc=z(np), lw_stem=z(np), lw_leaf=z(np),
        # PHS (use_hydrstress) inter-pass scratch — the arrays photosynthesis_hydrstress!
        # allocates internally (jmax_z_local/kn/psn_w*_z_*/rh_leaf_*). Held in the bundle
        # so their adjoints thread across the 4 decomposed PHS pass phases (Pass2 writes
        # phs_jmax → Pass3 reads it; Pass3 writes phs_w* → Pass4 reads them). bsun/bsha use
        # energyflux.bsun_patch/bsha_patch. Allocated for all canopy states (small).
        phs_jmax=z(np, 2, NLEVCAN), phs_kn=z(np),
        phs_wc_sun=z(np, NLEVCAN), phs_wj_sun=z(np, NLEVCAN), phs_wp_sun=z(np, NLEVCAN),
        phs_wc_sha=z(np, NLEVCAN), phs_wj_sha=z(np, NLEVCAN), phs_wp_sha=z(np, NLEVCAN),
        phs_rh_sun=z(np), phs_rh_sha=z(np))
end

# Bundle = the 10 differentiated canopy state structs + scratch. This is the object
# Enzyme differentiates (Duplicated(bundle, make_zero(bundle))).
cf_rev_bundle(canopystate, energyflux, frictionvel, temperature, solarabs, soilstate,
    waterfluxbulk, waterstatebulk, waterdiagbulk, photosyns, scratch) =
    (; canopystate, energyflux, frictionvel, temperature, solarabs, soilstate,
       waterfluxbulk, waterstatebulk, waterdiagbulk, photosyns, scratch)

# --------------------------------------------------------------------------
# Sub-phases. Each takes the bundle `b` + the Const auxiliary `aux` (indices,
# forcing, PFT/photosynthesis params, filter, flags) and mutates `b`. They call the
# exact same physics functions canopy_fluxes_core! uses.
# --------------------------------------------------------------------------
function cf_rev_init!(b, aux)
    sc = b.scratch; cs = b.canopystate; fv = b.frictionvel; tp = b.temperature
    FT = eltype(sc.air); fp = aux.filterp; fn = aux.fn; pd = aux.patch; fc = aux.forc; pft = aux.pft
    _launch!(_cf_init_zero_kernel!, sc.del_arr, sc.efeb, sc.wtlq0, sc.wtalq, sc.wtgq,
        sc.wtaq0, sc.obuold, b.energyflux.dhsdt_canopy_patch, b.energyflux.eflx_sh_stem_patch,
        fp; ndrange = fn)
    _launch!(_cf_no_biomass_kernel!, sc.frac_rad_abs_by_stem, sc.sa_leaf, sc.sa_stem,
        sc.sa_internal, sc.cp_leaf, zeros(FT, length(sc.rstem)), sc.rstem, fp, cs.elai_patch,
        cs.esai_patch; ndrange = fn)
    fv.rb1_patch[1:length(sc.air)] .= zero(FT)
    if aux.use_psn
        _launch!(_cf_daylength_kernel!, sc.dayl_factor, fp, pd.gridcell, fc.dayl, fc.max_dayl;
            ndrange = fn)
    end
    _launch!(_cf_z0_kernel!, cs.displa_patch, fv.z0mv_patch, fv.z0hv_patch, fv.z0qv_patch,
        fv.forc_hgt_u_patch, fv.forc_hgt_t_patch, fv.forc_hgt_q_patch, fp, pd.column,
        pd.gridcell, pd.itype, cs.elai_patch, cs.esai_patch, cs.htop_patch, fv.z0mg_col,
        fc.hgt_u, fc.hgt_t, fc.hgt_q, pft.z0v_LAImax, pft.z0v_Cs, pft.z0v_Cr, pft.z0v_c,
        pft.z0v_cw, 1; ndrange = fn)
    lwq_out = CfLwqOut(; air=sc.air, bir=sc.bir, cir=sc.cir, qsatl=sc.qsatl, el=sc.el,
        qsatldT=sc.qsatldT, co2_arr=sc.co2_arr, o2_arr=sc.o2_arr, taf=fv.taf_patch,
        qaf=fv.qaf_patch, ur=sc.ur, dth=sc.dth, dqh=sc.dqh, delq=sc.delq, dthv=sc.dthv, zldis=sc.zldis)
    lwq_p = CfLwqP(; emv=tp.emv_patch, t_veg=tp.t_veg_patch, thm=tp.thm_patch,
        forc_hgt_u=fv.forc_hgt_u_patch, displa=cs.displa_patch)
    lwq_c = CfLwqC(; emg=tp.emg_col, forc_lwrad=fc.lwrad, forc_pbot=fc.pbot,
        t_grnd=tp.t_grnd_col, forc_q=fc.q, qg=b.waterdiagbulk.qg_col, forc_th=fc.th)
    lwq_g = CfLwqG(; forc_pco2=fc.pco2, forc_po2=fc.po2, forc_u=fc.u_grc, forc_v=fc.v_grc)
    let be = _kernel_backend(sc.air)
        _cf_longwave_qsat_kernel!(be)(lwq_out, lwq_p, lwq_c, lwq_g, sc.nmozsgn,
            fp, pd.column, pd.gridcell, convert(FT, canopy_fluxes_params.wind_min); ndrange = fn)
        KA.synchronize(be)
    end
    cf_moninobukini_update!(fv, tp, pd, fp, fn, sc.ur, sc.dthv, sc.zldis, sc.tl_ini, sc.ts_ini)
    # PHS: recompute smp_l/hk_l from injected liquid soil water (production does this
    # before the PHS solve — HydrologyNoDrainage fills smp_l only AFTER canopy_fluxes).
    get(aux, :use_hydrstress, false) && cf_rev_phs_smp_l!(b, aux)
    return nothing
end

function cf_rev_friction!(b, aux, it::Int)
    fv = b.frictionvel; sc = b.scratch
    friction_velocity!(fv, aux.fn, aux.filterp[1:aux.fn], b.canopystate.displa_patch,
        fv.z0mv_patch, fv.z0hv_patch, fv.z0qv_patch, fv.obu_patch, it + 1, sc.ur,
        fv.um_patch, fv.ustar_patch, sc.temp1, sc.temp2, sc.temp12m, sc.temp22m, sc.fm;
        active = aux.active)
    return nothing
end

function cf_rev_resist!(b, aux)
    fv = b.frictionvel; sc = b.scratch; pft = aux.pft
    cf_resist_update!(fv, b.canopystate, b.temperature, b.waterdiagbulk, aux.patch,
        aux.filterp, aux.fn, aux.active, sc.temp1, sc.temp2, sc.tlbef, sc.del2, sc.del_arr,
        sc.rah, sc.raw, sc.uuc, sc.rb, sc.svpts, sc.eah, sc.el, pft.dleaf, pft.grnd_ch4,
        aux.forc.pbot, canopy_fluxes_params.csoilc, false, false, false, canopy_fluxes_params)
    return nothing
end

function cf_rev_psn!(b, aux, phase::String)
    sc = b.scratch; tp = b.temperature; ps = b.photosyns; pn = aux.psn
    parz = phase == "sun" ? pn.parsun_z : pn.parsha_z
    laiz = phase == "sun" ? pn.laisun_z : pn.laisha_z
    vcm  = phase == "sun" ? pn.vcmaxcint_sun : pn.vcmaxcint_sha
    photosynthesis!(ps, sc.svpts, sc.eah, sc.o2_arr, sc.co2_arr, sc.rb,
        b.energyflux.btran_patch, sc.dayl_factor, aux.forc.leafn, aux.forc_pbot_patch,
        tp.t_veg_patch, pn.t10, tp.thm_patch, pn.nrad, pn.tlai_z, b.canopystate.tlai_patch,
        parz, laiz, vcm, pn.o3coefv, pn.o3coefg, pn.c3psn, pn.leafcn, pn.flnr, pn.fnitr,
        pn.slatop, pn.mbbopt, pn.medlynintercept, pn.medlynslope, aux.ivt, aux.patch.column,
        aux.mask, 1:length(sc.air), phase)
    return nothing
end

# PHS (use_hydrstress) init helper: recompute smp_l + hk_l from the injected LIQUID
# soil water before the PHS photosynthesis solve. In production canopy_fluxes!
# (~1545-1569) this is done because HydrologyNoDrainage — which normally fills smp_l —
# runs AFTER canopy_fluxes; PHS photosynthesis reads smp_l (it does NOT compute it).
# Byte-identical to the production loop. Called from cf_rev_init! when use_hydrstress.
function cf_rev_phs_smp_l!(b, aux)
    ss = b.soilstate; wd = b.waterdiagbulk; ph = aux.phs
    FT = eltype(ss.smp_l_col)
    nlevsoi = ph.nlevsoi; nlevsno = aux.nlevsno
    smp_l = ss.smp_l_col; watsat = ss.watsat_col; sucsat = ss.sucsat_col
    bsw = ss.bsw_col; smpmin = ss.smpmin_col
    h2osoi_liqvol = wd.h2osoi_liqvol_col; hksat = ss.hksat_col; hk_l = ss.hk_l_col
    @inbounds for c in axes(smp_l, 1), j in 1:nlevsoi
        s_node = max(min(h2osoi_liqvol[c, nlevsno + j] / watsat[c, j], one(FT)), FT(0.01))
        smp_l[c, j] = max(smpmin[c], -sucsat[c, j] * s_node^(-bsw[c, j]))
        jp1 = min(nlevsoi, j + 1)
        s1 = (h2osoi_liqvol[c, nlevsno + j] + h2osoi_liqvol[c, nlevsno + jp1]) /
             (watsat[c, j] + watsat[c, jp1])
        s1 = min(one(FT), s1)
        hk_l[c, j] = hksat[c, j] * s1^(FT(2.0) * bsw[c, j] + FT(3.0))
    end
    return nothing
end

# PHS photosynthesis sub-phases — the use_hydrstress counterpart of cf_rev_psn!.
# Unlike the non-PHS path (two photosynthesis! calls), PHS solves sun+sha together via
# the Kennedy/Gentine fused vegwp Newton solve in photosynthesis_hydrstress!. That whole
# function is too large to reverse-differentiate as ONE Enzyme thunk on Julia 1.12 (the
# monolith trips EnzymeNoTypeError inside Pass 2 — the same wall canopy_fluxes_core! hits),
# but EACH of its four passes differentiates cleanly in isolation. So PHS photosynthesis is
# DECOMPOSED into four compositional phases (one Enzyme call each), exactly mirroring how
# the canopy energy balance is decomposed. The inter-pass scratch lives in the bundle
# (b.scratch.phs_* + energyflux.bsun/bsha) so its adjoint threads across the four phases.
# Arg mapping follows production photosynthesis_hydrstress! (canopy_fluxes.jl ~1746-1768):
# forc_rho is COLUMN-indexed (aux.forc.rho); smp_l/k_soil_root/rootfr/hk_l/hksat from
# soilstate; vegwp/lai/elai/esai/htop from canopystate; PHS root params from aux.phs;
# esat_tv=svpts, eair=eah, cair=co2_arr, oair=o2_arr, tgcm=thm (the canopy scratch).
# overrides default (CalibrationOverrides()) + the local-noon predicate are constructed
# in-phase (matching cf_rev_psn!'s reliance on photosynthesis! defaults); they feed only
# host-side pre-computation, not the differentiated arrays.

# PHS Pass 1: root-soil conductance k_soil_root (reads smp_l set by cf_rev_phs_smp_l!).
function cf_rev_psn_phs_pass1!(b, aux)
    sc = b.scratch; ph = aux.phs; cs = b.canopystate; ss = b.soilstate
    psn_phs_pass1_update!(b.photosyns, aux.mask, aux.patch.column, aux.ivt,
        ph.froot_carbon, ss.rootfr_patch, ph.dz, cs.tsai_patch, cs.tlai_patch,
        ph.froot_leaf, ph.root_radius, ph.root_density, ss.hksat_col, ss.hk_l_col,
        ss.smp_l_col, ph.z_col, ss.k_soil_root_patch, ss.root_conductance_patch,
        ss.soil_conductance_patch, 2.0, 0.25, ph.nlevsoi, 1:length(sc.air))
    return nothing
end

# PHS Pass 2: kinetics, N-profile, vcmax/jmax/tpu/kp/lmr → phs_jmax + ps.*_phs fields.
function cf_rev_psn_phs_pass2!(b, aux)
    sc = b.scratch; ph = aux.phs; pn = aux.psn; tp = b.temperature; ps = b.photosyns
    psn_phs_pass2_update!(ps, sc.phs_kn, sc.phs_jmax, aux.mask, aux.ivt,
        pn.c3psn, pn.mbbopt, aux.forc_pbot_patch, sc.o2_arr, pn.slatop, pn.leafcn,
        pn.flnr, pn.fnitr, sc.dayl_factor, pn.t10, tp.t_veg_patch, pn.tlai_z,
        pn.parsun_z, pn.parsha_z, pn.vcmaxcint_sun, pn.vcmaxcint_sha, pn.nrad,
        false, false, 0.015, NLEVCAN, ps.stomatalcond_mtd, ps.light_inhibit,
        ps.leafresp_method, CalibrationOverrides(), 1:length(sc.air))
    return nothing
end

# PHS Pass 3: the fused per-patch vegwp Newton solve (calcstress/spacF/spacA/getvegwp).
function cf_rev_psn_phs_pass3!(b, aux)
    sc = b.scratch; ph = aux.phs; pn = aux.psn; tp = b.temperature
    cs = b.canopystate; ss = b.soilstate; ef = b.energyflux; ps = b.photosyns
    # Match photosynthesis_hydrstress!'s bsun/bsha init (it fill!s 1.0 every call).
    fill!(ef.bsun_patch, one(eltype(ef.bsun_patch)))
    fill!(ef.bsha_patch, one(eltype(ef.bsha_patch)))
    psn_phs_pass3_update!(ps, aux.mask, aux.patch.column, aux.ivt, pn.nrad,
        pn.medlynslope, pn.medlynintercept, ph.crop_pft,
        aux.forc_pbot_patch, tp.thm_patch, sc.rb, sc.eah, sc.svpts, sc.co2_arr, sc.o2_arr,
        sc.qsatl, b.frictionvel.qaf_patch,
        cs.laisun_patch, cs.laisha_patch, cs.htop_patch, cs.tsai_patch, cs.elai_patch,
        cs.esai_patch, b.waterdiagbulk.fdry_patch, aux.forc.rho,
        pn.parsun_z, pn.parsha_z, sc.phs_jmax,
        pn.o3coefg, pn.o3coefg, pn.o3coefv, pn.o3coefv,
        ss.k_soil_root_patch, ss.smp_l_col, ph.z_col,
        cs.vegwp_patch, cs.vegwp_ln_patch, b.waterfluxbulk.wf.qflx_tran_veg_patch,
        ef.bsun_patch, ef.bsha_patch, sc.phs_rh_sun, sc.phs_rh_sha,
        sc.phs_wc_sun, sc.phs_wj_sun, sc.phs_wp_sun,
        sc.phs_wc_sha, sc.phs_wj_sha, sc.phs_wp_sha,
        ph.nlevsoi, ps.modifyphoto_and_lmr_forcrop, ps.stomatalcond_mtd,
        2.0e4, CalibrationOverrides(), (p) -> false, 1:length(sc.air))
    return nothing
end

# PHS Pass 4: canopy integration → per-patch psn/rs/lmr + LAI-weighted btran.
function cf_rev_psn_phs_pass4!(b, aux)
    sc = b.scratch; ph = aux.phs; pn = aux.psn; ef = b.energyflux; ps = b.photosyns
    psn_phs_pass4_update!(ps, aux.mask, pn.nrad, aux.ivt, ph.crop_pft,
        pn.laisun_z, pn.laisha_z, sc.rb, ef.btran_patch,
        sc.phs_wc_sun, sc.phs_wj_sun, sc.phs_wp_sun,
        sc.phs_wc_sha, sc.phs_wj_sha, sc.phs_wp_sha,
        ef.bsun_patch, ef.bsha_patch, ps.modifyphoto_and_lmr_forcrop, 1:length(sc.air))
    return nothing
end

function cf_rev_energy!(b, aux)
    sc = b.scratch; fc = aux.forc
    cf_energy_update!(b.canopystate, b.energyflux, b.frictionvel, b.temperature,
        b.solarabs, b.soilstate, b.waterfluxbulk, b.waterstatebulk, b.waterdiagbulk,
        b.photosyns, aux.patch, aux.col, aux.filterp, aux.fn, aux.active, sc.rah, sc.raw,
        sc.rb, sc.rstem, sc.sa_leaf, sc.sa_stem, sc.sa_internal, sc.frac_rad_abs_by_stem,
        sc.air, sc.bir, sc.cir, sc.cp_leaf, sc.tl_ini, sc.tlbef, sc.zldis, sc.temp1, sc.temp2,
        sc.ur, sc.efeb, sc.wtg, sc.wtl0, sc.wta0, sc.wtstem0, sc.wtga, sc.wtal, sc.lw_stem,
        sc.lw_leaf, sc.wtgq, sc.wtlq0, sc.wtaq0, sc.wtalq, sc.efe, sc.dt_veg, sc.del_arr,
        sc.err_arr, sc.qsatl, sc.el, sc.qsatldT, sc.dth, sc.dqh, sc.delq, sc.obuold,
        sc.nmozsgn, fc.q, fc.th, fc.pbot, fc.rho, aux.soilevap_beta, aux.soil_resis_sl14,
        false, false, aux.nlevsno, aux.dtime, canopy_fluxes_params)
    return nothing
end

# Ordered phase list. Each entry is (phase_fn, const_args::Tuple): the sub-phase
# function plus the EXTRA (non-differentiated) args beyond the bundle. These are
# passed to Enzyme.autodiff as explicit Enzyme.Const arguments rather than captured
# in a closure — Enzyme analyzes a closure's captured environment for activity (and
# trips on the captured `aux` NamedTuple's nested structs/BitVectors), whereas an
# explicit Const arg is cleanly excluded from differentiation.
function cf_rev_phases(aux, N::Int)
    use_phs = get(aux, :use_hydrstress, false)
    phases = Any[(cf_rev_init!, (aux,))]
    for it in 0:(N-1)
        push!(phases, (cf_rev_friction!, (aux, it)))
        push!(phases, (cf_rev_resist!, (aux,)))
        if use_phs
            # PHS photosynthesis decomposed into its 4 passes (one Enzyme call each) —
            # the monolith photosynthesis_hydrstress! won't reverse as a single thunk.
            push!(phases, (cf_rev_psn_phs_pass1!, (aux,)))
            push!(phases, (cf_rev_psn_phs_pass2!, (aux,)))
            push!(phases, (cf_rev_psn_phs_pass3!, (aux,)))
            push!(phases, (cf_rev_psn_phs_pass4!, (aux,)))
        elseif aux.use_psn
            push!(phases, (cf_rev_psn!, (aux, "sun")))
            push!(phases, (cf_rev_psn!, (aux, "sha")))
        end
        push!(phases, (cf_rev_energy!, (aux,)))
    end
    return phases
end

"""
    canopy_rev_forward!(b, aux, N)

Run the decomposed canopy forward: init kernels + `N` Newton iterations of
friction → resist → [photosynthesis sun/sha] → energy on the bundle `b`
(active≡true throughout). Numerically reproduces canopy_fluxes_core!'s iterated
state once `N` ≥ the converged iteration count.
"""
function canopy_rev_forward!(b, aux, N::Int)
    for (f, cargs) in cf_rev_phases(aux, N)
        f(b, cargs...)
    end
    return nothing
end

"""
    canopy_rev_converged_n(b, aux; itmax=canopy_fluxes_ctrl.itmax_canopy_fluxes) -> N

CONVERGENCE-AWARE iteration count for the compositional canopy reverse. Runs the
decomposed forward (init + friction→resist→[psn]→energy) one Newton iteration at a
time on a throwaway DEEPCOPY of `b`, applying the SAME convergence test
canopy_fluxes_core! uses — `det = max(del_arr, del2) < DTMIN_CANOPY` and
`dele = |efe - efeb| < DLEMIN_CANOPY`, only after `itlef > ITMIN_CANOPY` — and
returns the smallest `N` at which every active (in-filter) patch has converged
(clamped to `[ITMIN_CANOPY+1, itmax]`). Multi-patch aware: a patch that converges
early stops contributing to the "all converged" test once its metric drops below
tolerance, exactly as the production `active[p]=false` mask does, so `N` is the
MAX over per-patch convergence counts — the right fixed `N` for the differentiate-
through-the-converged-iterate scheme (the early patch sits at its stable fixed point
for the remaining iterations, so its primal + gradient are unchanged).

This removes the need for the caller to hard-code `n_canopy`: it adapts to whatever
the forcing/state make the Newton solve take this step. `b` is NOT mutated.
"""
function canopy_rev_converged_n(b, aux;
                                itmax::Int = canopy_fluxes_ctrl.itmax_canopy_fluxes)
    probe = deepcopy(b)
    sc = probe.scratch
    fp = aux.filterp[1:aux.fn]
    cf_rev_init!(probe, aux)
    conv = falses(length(fp))                    # per-(in-filter)-patch converged flag
    for itlef in 1:itmax
        cf_rev_friction!(probe, aux, itlef - 1)
        cf_rev_resist!(probe, aux)
        if aux.use_psn
            cf_rev_psn!(probe, aux, "sun"); cf_rev_psn!(probe, aux, "sha")
        end
        cf_rev_energy!(probe, aux)
        if itlef > ITMIN_CANOPY                  # production gates the test the same way
            for (k, p) in enumerate(fp)
                conv[k] && continue
                dele = abs(sc.efe[p] - sc.efeb[p])
                det  = max(sc.del_arr[p], sc.del2[p])
                (det < DTMIN_CANOPY && dele < DLEMIN_CANOPY) && (conv[k] = true)
            end
            all(conv) && return itlef
        end
    end
    return itmax
end

"""
    canopy_rev_gradient!(b, aux, N; seed=:t_veg_patch)

Compositional reverse pass: forward-sweep `b` through the phases with a deepcopy
checkpoint before each, seed the adjoint at the chosen output field
(`L = sum(abs2, b.temperature.<seed>)`), then back-propagate one separate
`Enzyme.autodiff` call per sub-phase. Returns the gradient bundle `db` (same
structure as `b`); e.g. `db.temperature.t_grnd_col` is dL/d(initial t_grnd).

If `N === nothing` the iteration count is auto-detected with
[`canopy_rev_converged_n`](@ref) — the convergence-aware path that adapts to the
Newton solve rather than hard-coding a fixed count.

`b` is left holding the final forward state.
"""
function canopy_rev_gradient!(b, aux, N::Union{Int,Nothing} = nothing; seed::Symbol = :t_veg_patch)
    Nuse = N === nothing ? canopy_rev_converged_n(b, aux) : N
    seed_bang!(db, b) = (getfield(db.temperature, seed) .= 2 .* getfield(b.temperature, seed))
    return compositional_reverse!(cf_rev_phases(aux, Nuse), b, seed_bang!)
end

"""
    compositional_reverse!(phases, b, seed_bang!) -> db

Generic checkpointed compositional reverse pass — the reusable engine behind the
canopy reverse and the building block for a whole-`clm_drv!` reverse (each driver
phase, e.g. `soil_temperature!`, plugs in as one more entry in `phases`).

`phases` is an ordered list of `(phase_fn, const_args::Tuple)`; each `phase_fn` is
called `phase_fn(b, const_args...)` and mutates the differentiated bundle `b`
in place. The forward sweep runs every phase, deepcopy-checkpointing `b` before
each. `seed_bang!(db, b)` seeds the adjoint shadow `db` from the final forward
state (e.g. `dL = 2·output`). The reverse sweep then runs ONE separate
`Enzyme.autodiff` call per phase, last→first, with the non-differentiated
`const_args` passed as `Enzyme.Const` (NOT closure-captured — Enzyme analyzes a
closure's captured environment for activity and trips on nested structs/BitVectors).
Returns the gradient bundle `db`.
"""
function compositional_reverse!(phases, b, seed_bang!)
    Enzyme.API.strictAliasing!(false)
    revmode = Enzyme.set_runtime_activity(Enzyme.Reverse)
    # Force the photosynthesis Pass-3 ci-solve onto its plain host loop for the whole
    # reverse pass: the KA kernel-launch path makes Enzyme reverse segfault on Julia
    # 1.12. The host loop is byte-identical (shared `_photosynth_ci_body!`), so both the
    # forward checkpoint sweep and the per-phase autodiff stay consistent. Restored in
    # `finally` so non-AD code keeps the KA kernel.
    _prev_psnci = _PSN_CI_AD_HOSTLOOP[]
    _PSN_CI_AD_HOSTLOOP[] = true
    try
        checkpoints = Any[]
        for (f, cargs) in phases
            push!(checkpoints, deepcopy(b)); f(b, cargs...)
        end
        db = Enzyme.make_zero(b)
        seed_bang!(db, b)
        for k in length(phases):-1:1
            f, cargs = phases[k]
            bk = deepcopy(checkpoints[k])
            Enzyme.autodiff(revmode, f, Enzyme.Const, Enzyme.Duplicated(bk, db),
                            map(Enzyme.Const, cargs)...)
        end
        return db
    finally
        _PSN_CI_AD_HOSTLOOP[] = _prev_psnci
    end
end

# Reverse ONE timestep's ordered phase list into an EXISTING adjoint `db`, starting from
# the step-entry state `b_start`. Re-runs the step forward with per-phase (fine) checkpoints,
# then back-propagates last→first accumulating into `db` — identical inner mechanics to
# compositional_reverse!, but it neither allocates nor seeds `db` (the caller threads it across
# steps). `b_start` is left holding the step's final forward state. Returns nothing.
function _reverse_step_into!(phases, b_start, db, revmode)
    fine = Any[]
    for (f, cargs) in phases
        push!(fine, deepcopy(b_start)); f(b_start, cargs...)
    end
    for k in length(phases):-1:1
        f, cargs = phases[k]
        bk = deepcopy(fine[k])
        Enzyme.autodiff(revmode, f, Enzyme.Const, Enzyme.Duplicated(bk, db),
                        map(Enzyme.Const, cargs)...)
    end
    return nothing
end

"""
    multistep_reverse!(steps, b, seed_bang!) -> db

Multi-timestep checkpointed reverse pass — propagates the adjoint across a TRAJECTORY of
timesteps, giving `d(L of the final state)/d(initial state)` (the state-to-state sensitivity
through the whole horizon, e.g. d(end-of-season soil temperature)/d(its initial profile)).

`steps` is an ordered list of per-timestep phase lists; each entry is exactly what
`compositional_reverse!` consumes (an ordered list of `(phase_fn, const_args::Tuple)`). The
phases mutate the shared bundle `b` in place, so step `s+1` continues from the state step `s`
produced — the multi-step coupling is automatic. With a single fixed-forcing step list this is
`fill(phases, N)`; with per-step forcing each entry captures that step's Const aux.

TWO-LEVEL CHECKPOINTING bounds memory to O(n_steps) coarse checkpoints + O(phases-per-step)
fine checkpoints (one step's worth, transient), instead of the O(n_steps · phases) that
`compositional_reverse!` on the flattened `vcat(steps...)` would hold. The forward sweep
deepcopy-checkpoints `b` only at each STEP boundary; the reverse sweep restores each
step-entry checkpoint and RECOMPUTES that step's fine (per-phase) checkpoints on the fly
before back-propagating through it. The single adjoint `db` is threaded across all steps:
seeded from the final state, after reversing step `s` it holds `dL/d(state entering step s)`
= `dL/d(state leaving step s-1)`, so after step 1 it is `dL/d(initial state)`.

Equivalent to (and three-way cross-checked against FD and) `compositional_reverse!(vcat(steps...), …)`,
but with bounded fine-checkpoint memory. `b` is left holding the final forward state.
Returns the gradient bundle `db` (same structure as `b`).
"""
function multistep_reverse!(steps, b, seed_bang!)
    Enzyme.API.strictAliasing!(false)
    revmode = Enzyme.set_runtime_activity(Enzyme.Reverse)
    # Forward sweep: checkpoint ONLY at step boundaries (coarse), run each step's phases.
    step_checkpoints = Any[]
    for phases in steps
        push!(step_checkpoints, deepcopy(b))
        for (f, cargs) in phases
            f(b, cargs...)
        end
    end
    # Seed the adjoint from the final forward state.
    db = Enzyme.make_zero(b)
    seed_bang!(db, b)
    # Reverse sweep: step last→first. Restore the step-entry checkpoint, recompute its fine
    # checkpoints, back-propagate through the step accumulating into the threaded `db`.
    for s in length(steps):-1:1
        _reverse_step_into!(steps[s], deepcopy(step_checkpoints[s]), db, revmode)
    end
    return db
end

# Advance a bundle forward through steps[lo+1 .. hi] (0-indexed steps lo..hi-1), in place.
function _advance_steps!(steps, b, lo::Int, hi::Int)
    for s in (lo + 1):hi, (f, cargs) in steps[s]
        f(b, cargs...)
    end
    return b
end

"""
    multistep_reverse_binomial!(steps, b, seed_bang!; peak_checkpoints=nothing) -> db

Same multi-timestep reverse as [`multistep_reverse!`](@ref) — `d(L of final state)/d(initial
state)` over the trajectory — but with **logarithmic** coarse-checkpoint memory for very long
horizons. Where `multistep_reverse!` stores all `n_steps` step-boundary snapshots
simultaneously (O(n_steps) memory), this uses RECURSIVE BISECTION (divide-and-conquer /
binomial checkpointing): it holds only O(log n_steps) snapshots live at once, at the cost of
O(n_steps · log n_steps) total step recomputation (vs O(n_steps) for the flat scheme). The
returned gradient is IDENTICAL — only the checkpoint/recompute schedule differs.

`rev_seg(lo, hi, b_lo)` reverses the half-open step range `[lo, hi)` given the state `b_lo`
entering step `lo`: a single step is reversed directly (fine within-step recompute); otherwise
it advances a COPY of `b_lo` to the midpoint, reverses the RIGHT half from there, frees that
midpoint snapshot, then reverses the LEFT half from the still-held `b_lo`. Right-then-left keeps
the single-step reverses in strict last→first order, so the shared adjoint `db` threads exactly
as in the flat engine. At most one snapshot per recursion level is live → O(log n_steps) peak.

`peak_checkpoints::Ref{Int}` (optional) is filled with the measured peak number of
simultaneously-held coarse snapshots — for the validation scripts to exhibit the memory bound.
This is the building block for season/year horizons where O(n_steps) snapshots won't fit.
"""
function multistep_reverse_binomial!(steps, b, seed_bang!;
                                     peak_checkpoints::Union{Nothing,Base.RefValue{Int}} = nothing)
    Enzyme.API.strictAliasing!(false)
    revmode = Enzyme.set_runtime_activity(Enzyme.Reverse)
    N = length(steps)
    b0 = deepcopy(b)                       # the only persistent snapshot: the initial state
    # Forward once (no stored boundaries) to the final state, to seed the adjoint.
    _advance_steps!(steps, b, 0, N)
    db = Enzyme.make_zero(b)
    seed_bang!(db, b)
    live = 0; peak = 0
    bump(d) = (live += d; peak = max(peak, live); nothing)
    function rev_seg(lo::Int, hi::Int, b_lo)
        if hi - lo == 1
            _reverse_step_into!(steps[lo + 1], deepcopy(b_lo), db, revmode)
            return nothing
        end
        mid = (lo + hi) ÷ 2
        b_mid = _advance_steps!(steps, deepcopy(b_lo), lo, mid); bump(1)
        rev_seg(mid, hi, b_mid)            # later steps first (keeps db order last→first)
        bump(-1)                           # b_mid no longer needed
        rev_seg(lo, mid, b_lo)
        return nothing
    end
    bump(1)                                # b0 held for the whole reverse
    N >= 1 && rev_seg(0, N, b0)
    peak_checkpoints === nothing || (peak_checkpoints[] = peak)
    return db
end
