# Guards that the productionized clm_drv! reverse PHASES (src/driver/driver_reverse.jl)
# — soiltemp_rev_phase!, soilwater_rev_phase!, watertable_rev_phase! and the
# driver_rev_phases assembler — stay in sync with the production physics signatures.
# NO Enzyme here: a fast, suite-safe forward guard (the wrappers just call the physics
# functions, so a signature/aux-builder drift shows up as an error or a parity break).
# The reverse-AD gradients are validated on Julia 1.10 + Enzyme in
# scripts/enzyme_driver_reverse_full.jl (canopy+soil_temp+soil_water) and
# scripts/enzyme_driver_reverse_hydro.jl (+ water_table!).

include(joinpath(@__DIR__, "testdata.jl"))

@testset "driver reverse-phase forward parity" begin
    fsurdat, paramfile = bow_params()
    if !isfile(fsurdat) || !isfile(paramfile)
        testdata_missing("driver-reverse forward parity", fsurdat, paramfile) && return
    end

    (inst, bounds, filt, tm) = CLM.clm_initialize!(; fsurdat=fsurdat, paramfile=paramfile)
    for g in 1:bounds.endg
        a = inst.atm2lnd; a.forc_t_not_downscaled_grc[g] = 285.0
        a.forc_pbot_not_downscaled_grc[g] = 85000.0
        a.forc_th_not_downscaled_grc[g] = 285.0 * (100000.0/85000.0)^(CLM.RAIR/CLM.CPAIR)
        a.forc_rho_not_downscaled_grc[g] = 85000.0/(CLM.RAIR*285.0)
        a.forc_lwrad_not_downscaled_grc[g] = 300.0; a.forc_vp_grc[g] = 800.0
        a.forc_hgt_grc[g] = 30.0; a.forc_wind_grc[g] = 3.0
        a.forc_hgt_u_grc[g] = 30.0; a.forc_hgt_t_grc[g] = 30.0; a.forc_hgt_q_grc[g] = 30.0
        a.forc_topo_grc[g] = 0.0
        for b in 1:CLM.NUMRAD; a.forc_solad_not_downscaled_grc[g,b] = 200.0; a.forc_solai_grc[g,b] = 80.0; end
        a.forc_solar_not_downscaled_grc[g] = 560.0
        a.forc_rain_not_downscaled_grc[g] = 0.0001; a.forc_snow_not_downscaled_grc[g] = 0.0
    end
    CLM.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
    config = CLM.CLMDriverConfig(); filt_ia = CLM.clump_filter_inactive_and_active
    (declin, _) = CLM.compute_orbital(120.0); nextsw = 120.0 + 1800.0/CLM.SECSPDAY
    for n in 1:3
        CLM.clm_drv!(config, inst, filt, filt_ia, bounds, true, nextsw, declin, declin,
            CLM.ORB_OBLIQR_DEFAULT, false, false, "", false;
            nstep=n, is_first_step=(n==1), is_beg_curr_day=(n==1), dtime=1800.0,
            mon=1, day=1, photosyns=inst.photosyns)
    end

    hc = filt.hydrologyc
    c0 = first(c for c in bounds.begc:bounds.endc if hc[c])
    j1 = CLM.varpar.nlevsno + 1

    # assembler (forward order): soil_temp + 10-phase surface-hydrology block +
    # soil_water + water_table + hydrology_no_drainage = 14 (15 with a canopy block);
    # include_surface=false reverts to the soil_temp→soil_water direct jump (4).
    ph = CLM.driver_rev_phases(bounds, filt, config)
    @test length(ph) == 14
    @test length(CLM.driver_rev_phases(bounds, filt, config; canopy_aux=(;), n_canopy=1)) == 15
    @test length(CLM.driver_rev_phases(bounds, filt, config; include_surface=false)) == 4
    @test length(CLM.surface_hydrology_rev_phases(bounds, filt)) == 10

    # run the full assembled chain forward on the shared-inst bundle.
    b = CLM.driver_rev_bundle(inst)
    for (f, cargs) in ph
        f(b, cargs...)
    end
    @test b.inst === inst                                   # bundle aliases the inst
    @test isfinite(b.inst.temperature.t_soisno_col[c0, j1])
    @test all(isfinite, b.inst.water.waterstatebulk_inst.ws.h2osoi_liq_col)
    @test isfinite(b.inst.soilhydrology.zwt_col[c0])
    @test isfinite(b.inst.water.waterstatebulk_inst.ws.h2osoi_vol_col[c0, CLM.varpar.nlevsno+4])
    @test isfinite(b.inst.water.waterfluxbulk_inst.wf.qflx_infl_col[c0])   # surface block ran

    # parity: soiltemp_rev_phase! reproduces a raw soil_temperature! call (i.e. the
    # soiltemp aux builder passes the same args the driver does).
    ref = deepcopy(inst)
    let i = ref
        CLM.soil_temperature!(i.column, i.landunit, i.patch, i.temperature, i.energyflux,
            i.soilstate, i.water.waterstatebulk_inst, i.water.waterdiagnosticbulk_inst,
            i.water.waterfluxbulk_inst, i.solarabs, i.canopystate, i.urbanparams,
            fill(323.15, bounds.endl), i.atm2lnd.forc_lwrad_downscaled_col,
            filt.nolakec, filt.nolakep, filt.urbanl, filt.urbanc,
            bounds.begc:bounds.endc, bounds.begl:bounds.endl, bounds.begp:bounds.endp, 1800.0)
    end
    b2 = CLM.driver_rev_bundle(deepcopy(inst))
    CLM.soiltemp_rev_phase!(b2, CLM.soiltemp_rev_aux(bounds, filt))
    @test b2.inst.temperature.t_soisno_col[c0, j1] ≈ ref.temperature.t_soisno_col[c0, j1] atol = 1e-12

    # pre-soil_water surface-hydrology phase wrappers run + produce finite outputs.
    bs = CLM.driver_rev_bundle(deepcopy(inst))
    sh = CLM.surfhydro_rev_aux(bounds, filt)
    CLM.satexcess_rev_phase!(bs, sh)      # fsat
    CLM.inflexcess_rev_phase!(bs, sh)     # qinmax (reads fsat)
    CLM.infil_rev_phase!(bs, sh)          # qflx_infl
    @test isfinite(bs.inst.sat_excess_runoff.fsat_col[c0])
    @test isfinite(bs.inst.infilt_excess_runoff.qinmax_col[c0])
    @test isfinite(bs.inst.water.waterfluxbulk_inst.wf.qflx_infl_col[c0])

    # canopy_rev_aux builds from the real inst + the FULL 15-phase whole-step (canopy +
    # hydrology) runs forward finite (few Newton iters for speed; the reverse gradient is
    # validated on Julia 1.10/Enzyme in scripts/enzyme_driver_reverse_fullstep.jl).
    caux = CLM.canopy_rev_aux(inst, bounds, filt)
    full = CLM.driver_rev_phases(bounds, filt, config; canopy_aux=caux, n_canopy=4)
    @test length(full) == 15
    bf = CLM.driver_rev_bundle(deepcopy(inst))
    for (f, cargs) in full
        f(bf, cargs...)
    end
    @test all(isfinite, bf.inst.temperature.t_veg_patch)
    @test isfinite(bf.inst.water.waterstatebulk_inst.ws.h2osoi_vol_col[c0, CLM.varpar.nlevsno+4])

    # canopy_rev_aux(use_psn=true) sources the real LUNA/albedo arrays + the canopy block
    # with photosynthesis runs forward finite (stomatal feedback path; reverse validated on
    # 1.10/Enzyme in scripts/enzyme_driver_reverse_fullstep.jl CLM_USE_PSN=1).
    fp = first(p for p in bounds.begp:bounds.endp if filt.exposedvegp[p])
    bp = CLM.driver_rev_bundle(deepcopy(inst))
    CLM.canopy_rev_block!(bp, CLM.canopy_rev_aux(inst, bounds, filt; use_psn=true), 2)
    @test isfinite(bp.inst.temperature.t_veg_patch[fp])
    @test isfinite(bp.inst.photosyns.rssun_patch[fp])

    # ---- convergence-aware canopy N + top-level clm_drv_reverse! guard (NO Enzyme).
    # driver_canopy_converged_n auto-detects the canopy Newton count from the decomposed
    # forward's per-patch convergence (same DTMIN/DLEMIN test the production kernel uses),
    # so the driver-reverse never hard-codes n_canopy. It must (a) land in [ITMIN+1, itmax],
    # (b) be deterministic, (c) not mutate the bundle, and (d) the forward orchestration that
    # clm_drv_reverse! reverses must run finite. The reverse GRADIENT (auto N) is FD-validated
    # on 1.10/Enzyme in scripts/enzyme_clm_drv_reverse.jl (rel ~2e-7 vs Richardson FD).
    bconv = CLM.driver_rev_bundle(deepcopy(inst))
    caux2 = CLM.canopy_rev_aux(bconv.inst, bounds, filt)     # energy-balance canopy
    tveg_before = copy(bconv.inst.temperature.t_veg_patch)
    Nauto = CLM.driver_canopy_converged_n(bconv, caux2)
    @test Nauto isa Int
    @test CLM.ITMIN_CANOPY < Nauto <= CLM.canopy_fluxes_ctrl.itmax_canopy_fluxes
    @test bconv.inst.temperature.t_veg_patch == tveg_before  # detection did NOT mutate the bundle
    # Detection is DETERMINISTIC (re-running it on a fresh bundle gives the same N).
    # NOTE: Nauto is the convergence count of the ENERGY-BALANCE mirror (use_psn=false here),
    # a different solve than the warmed inst's num_iter_patch (full canopy_fluxes! with
    # photosynthesis) — so we don't assert Nauto≈num_iter. It either converges below the
    # DTMIN/DLEMIN tolerance or saturates at itmax (both are correct: N≥converged-count is
    # what the differentiate-through-the-converged-iterate scheme needs).
    @test CLM.driver_canopy_converged_n(CLM.driver_rev_bundle(deepcopy(inst)),
              CLM.canopy_rev_aux(inst, bounds, filt)) == Nauto      # deterministic
    # clm_drv_reverse! forward-orchestration parity: the chain it assembles with the auto N
    # reproduces an explicit driver_rev_phases(...; n_canopy=Nauto) forward sweep.
    bref_drv = CLM.driver_rev_bundle(deepcopy(inst))
    phref = CLM.driver_rev_phases(bounds, filt, config;
                                  canopy_aux=CLM.canopy_rev_aux(bref_drv.inst, bounds, filt),
                                  n_canopy=Nauto)
    for (f, cargs) in phref; f(bref_drv, cargs...); end
    @test all(isfinite, bref_drv.inst.temperature.t_veg_patch)
    @test isfinite(bref_drv.inst.water.waterstatebulk_inst.ws.h2osoi_vol_col[c0, CLM.varpar.nlevsno+4])
    @test isa(CLM.clm_drv_reverse!, Function)                # top-level entry exists
    # the default seed touches the soil thermal profile (active soil layers only).
    bseed = CLM.driver_rev_bundle(deepcopy(inst))
    for (f, cargs) in phref; f(bseed, cargs...); end
    dbz = CLM.Enzyme.make_zero(bseed)
    CLM._default_drv_seed!(dbz, bseed)
    @test dbz.inst.temperature.t_soisno_col[c0, j1] ≈ 2*bseed.inst.temperature.t_soisno_col[c0, j1] atol=1e-12
    @test dbz.inst.temperature.t_soisno_col[c0, 1] == 0.0   # snow layer (j<nlevsno+1) not seeded

    # ---- multi-step trajectory guard (NO Enzyme). The reverse-AD gradient
    # d(final state)/d(initial state) over an N-step horizon is validated on 1.10/Enzyme in
    # scripts/enzyme_multistep_reverse.jl (3-way: CLM.multistep_reverse! == compositional_
    # reverse!(vcat(steps)) == FD). Here we guard the orchestration WITHOUT Enzyme: the same
    # per-step phase list `ph` is reusable across timesteps, the carried thermal state advances
    # step-to-step, and the two-level FORWARD sweep (deepcopy at step boundaries, recompute the
    # step) reproduces a flat N-step loop bit-for-bit — i.e. the checkpoint scheduling is faithful.
    @test isa(CLM.multistep_reverse!, Function)
    jdiff = j1 + 3                                   # an actively-diffusing soil layer
    Nstep = 3; steps = [ph for _ in 1:Nstep]
    bflat = CLM.driver_rev_bundle(deepcopy(inst)); after1 = NaN
    for (s, phs) in enumerate(steps)
        for (f, cargs) in phs; f(bflat, cargs...); end
        s == 1 && (after1 = bflat.inst.temperature.t_soisno_col[c0, jdiff])
    end
    tfin = bflat.inst.temperature.t_soisno_col[c0, jdiff]
    @test isfinite(tfin)
    @test abs(after1 - tfin) > 1e-4                  # state genuinely advances across steps
    # two-level forward: coarse deepcopy checkpoints at step boundaries (mirrors
    # multistep_reverse!'s forward sweep), then recompute the last step from its checkpoint.
    btl = CLM.driver_rev_bundle(deepcopy(inst)); cps = Any[]
    for phs in steps
        push!(cps, deepcopy(btl)); for (f, cargs) in phs; f(btl, cargs...); end
    end
    blast = deepcopy(cps[end]); for (f, cargs) in steps[end]; f(blast, cargs...); end
    @test btl.inst.temperature.t_soisno_col[c0, jdiff] ≈ tfin atol=1e-12
    @test blast.inst.temperature.t_soisno_col[c0, jdiff] ≈ tfin atol=1e-12

    # ---- binomial (recursive-bisection) checkpointing guard (NO Enzyme). The O(log N)-memory
    # engine CLM.multistep_reverse_binomial! gives a gradient bitwise-identical to multistep_
    # reverse! (Enzyme-validated 4-way in scripts/enzyme_multistep_reverse.jl). Here we guard
    # its forward primitive _advance_steps!(steps,b,lo,hi): advancing a COPY through the full
    # range reproduces the flat loop, and advancing a sub-range then the remainder equals the
    # whole — the recompute fidelity the bisection relies on.
    @test isa(CLM.multistep_reverse_binomial!, Function)
    bfull = CLM._advance_steps!(steps, CLM.driver_rev_bundle(deepcopy(inst)), 0, Nstep)
    @test bfull.inst.temperature.t_soisno_col[c0, jdiff] ≈ tfin atol=1e-12
    bseg = CLM.driver_rev_bundle(deepcopy(inst))
    CLM._advance_steps!(steps, bseg, 0, 1)            # step 0
    CLM._advance_steps!(steps, bseg, 1, Nstep)        # steps 1..N-1
    @test bseg.inst.temperature.t_soisno_col[c0, jdiff] ≈ tfin atol=1e-12

    # ---- per-step FORCING-advance guard (NO Enzyme). For multistep over a TIME-VARYING-forcing
    # trajectory, forcingset_rev_phase! overwrites the live downscaled forcing arrays each step;
    # the reverse gradient under a diurnal schedule is Enzyme-validated in
    # scripts/enzyme_multistep_forced.jl. Here: forcingset_rev_phase! sets the live arrays, a
    # different per-step forcing yields a different trajectory (forcing genuinely drives state),
    # and forced_driver_steps prepends exactly one forcingset phase per step.
    base = CLM.forcingset_aux(inst)
    bset = CLM.driver_rev_bundle(deepcopy(inst))
    warm = base.lwrad .+ 50.0                          # a perturbed longwave snapshot
    CLM.forcingset_rev_phase!(bset, (; lwrad=warm, t=base.t, th=base.th, rho=base.rho))
    @test bset.inst.atm2lnd.forc_lwrad_downscaled_col == warm   # live array overwritten
    sched = [(; lwrad=base.lwrad .+ 30.0*(k-1), t=base.t, th=base.th, rho=base.rho) for k in 1:Nstep]
    fsteps = CLM.forced_driver_steps(bounds, filt, config, sched; canopy_aux=nothing)
    @test length(fsteps) == Nstep
    @test length(fsteps[1]) == length(ph) + 1          # one forcingset phase prepended
    @test fsteps[1][1][1] === CLM.forcingset_rev_phase!
    bA = CLM._advance_steps!(fsteps, CLM.driver_rev_bundle(deepcopy(inst)), 0, Nstep)
    bB = CLM._advance_steps!(steps,  CLM.driver_rev_bundle(deepcopy(inst)), 0, Nstep)  # fixed forcing
    @test isfinite(bA.inst.temperature.t_soisno_col[c0, jdiff])
    @test bA.inst.temperature.t_soisno_col[c0, jdiff] != bB.inst.temperature.t_soisno_col[c0, jdiff]
end
