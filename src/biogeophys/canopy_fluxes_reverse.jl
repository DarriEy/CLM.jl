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
        sa_stem=z(np), sa_leaf=z(np), sa_internal=z(np), uuc=z(np), lw_stem=z(np), lw_leaf=z(np))
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
    phases = Any[(cf_rev_init!, (aux,))]
    for it in 0:(N-1)
        push!(phases, (cf_rev_friction!, (aux, it)))
        push!(phases, (cf_rev_resist!, (aux,)))
        if aux.use_psn
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
    canopy_rev_gradient!(b, aux, N; seed=:t_veg_patch)

Compositional reverse pass: forward-sweep `b` through the phases with a deepcopy
checkpoint before each, seed the adjoint at the chosen output field
(`L = sum(abs2, b.temperature.<seed>)`), then back-propagate one separate
`Enzyme.autodiff` call per sub-phase. Returns the gradient bundle `db` (same
structure as `b`); e.g. `db.temperature.t_grnd_col` is dL/d(initial t_grnd).

`b` is left holding the final forward state.
"""
function canopy_rev_gradient!(b, aux, N::Int; seed::Symbol = :t_veg_patch)
    seed_bang!(db, b) = (getfield(db.temperature, seed) .= 2 .* getfield(b.temperature, seed))
    return compositional_reverse!(cf_rev_phases(aux, N), b, seed_bang!)
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
end
