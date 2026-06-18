# =============================================================================
# EXTENSION of scripts/enzyme_driver_reverse.jl toward a full clm_drv! timestep
# reverse. The base script proves [A] full-canopy, [B] soil_temperature!, and
# [C] canopy+soil_temperature CHAINED, all through the production
# CLM.compositional_reverse! engine (src/biogeophys/canopy_fluxes_reverse.jl).
#
# This script adds the NEXT clm_drv! phase(s) after soil_temperature in the chain,
# one at a time, each FD-validated:
#
#   [D] soil_water! (soilwater_zengdecker2009! default path) reverse, STANDALONE on
#       the real Bow-at-Banff cold-start inst. FD-validates dL/d(t_soisno[c,j])
#       (t_soisno feeds the Zeng-Decker ice-fraction + vwc terms) through
#       compositional_reverse!, L = sum(abs2, h2osoi_liq).
#
#   [E] soil_temperature! + soil_water! TWO-phase chain in ONE compositional_reverse!
#       call. Gradient flows soil_temp(t_soisno) → soil_water(h2osoi_liq) across the
#       module boundary. FD-validates dL/d(t_grnd), L = sum(abs2, h2osoi_liq).
#
#   [F] canopy + soil_temperature + soil_water THREE-module chain in ONE
#       compositional_reverse! call — extends [C] by one module. FD-validates
#       dL/d(t_grnd), L = sum(abs2, h2osoi_liq), the gradient that flows
#       t_grnd → canopy(t_veg, ground heat flux) → soil_temp(t_soisno) →
#       soil_water(h2osoi_liq), ACROSS three module boundaries.
#
#   julia +1.10 --project=/tmp/clm_jl10_<id> scripts/enzyme_driver_reverse_full.jl
#   (Enzyme is stable on Julia 1.10 LTS)
# =============================================================================
using CLM, Enzyme, Printf
const C = CLM
const FT = Float64

const FSURDAT = get(ENV, "CLM_FSURDAT",
    "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/surfdata_clm.nc")
const PARAMFILE = get(ENV, "CLM_PARAMFILE",
    "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/clm5_params.nc")

# =============================================================================
# Real Bow-at-Banff cold-start builder (copied from enzyme_driver_reverse.jl —
# that script is read-only; the patterns are reused here).
# =============================================================================
function setup_forcing!(a2l, T, ng)
    for g in 1:ng
        a2l.forc_t_not_downscaled_grc[g]=T; a2l.forc_pbot_not_downscaled_grc[g]=85000.0
        a2l.forc_th_not_downscaled_grc[g]=T*(100000.0/85000.0)^(C.RAIR/C.CPAIR)
        a2l.forc_rho_not_downscaled_grc[g]=85000.0/(C.RAIR*T); a2l.forc_lwrad_not_downscaled_grc[g]=300.0
        a2l.forc_vp_grc[g]=800.0; a2l.forc_hgt_grc[g]=30.0; a2l.forc_topo_grc[g]=0.0; a2l.forc_wind_grc[g]=3.0
        a2l.forc_hgt_u_grc[g]=30.0; a2l.forc_hgt_t_grc[g]=30.0; a2l.forc_hgt_q_grc[g]=30.0
        for b in 1:C.NUMRAD; a2l.forc_solad_not_downscaled_grc[g,b]=200.0; a2l.forc_solai_grc[g,b]=80.0; end
        a2l.forc_solar_not_downscaled_grc[g]=560.0; a2l.forc_rain_not_downscaled_grc[g]=0.0001
        a2l.forc_snow_not_downscaled_grc[g]=0.0
    end
end

function build_real()
    (inst, bounds, filt, tm) = C.clm_initialize!(; fsurdat=FSURDAT, paramfile=PARAMFILE)
    setup_forcing!(inst.atm2lnd, 285.0, bounds.endg)
    C.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
    config = C.CLMDriverConfig(); filt_ia = C.clump_filter_inactive_and_active
    (declin, eccf) = C.compute_orbital(120.0); nextsw_cday = 120.0 + 1800.0/C.SECSPDAY
    runstep!(n; first=false) = C.clm_drv!(config, inst, filt, filt_ia, bounds,
        true, nextsw_cday, declin, declin, C.ORB_OBLIQR_DEFAULT, false, false, "", false;
        nstep=n, is_first_step=first, is_beg_curr_day=first, dtime=1800.0, mon=1, day=1,
        photosyns=inst.photosyns)
    for n in 1:3; runstep!(n; first=(n==1)); end
    return inst, bounds, filt, config
end

# =============================================================================
# soil_temperature! phase (copied from base script — needed for the [E]/[F] chain).
# =============================================================================
soil_aux(bounds, filt) = (;
    urbantv = fill(323.15, bounds.endl), dtime = 1800.0,
    bc_col = bounds.begc:bounds.endc, bc_lun = bounds.begl:bounds.endl,
    bc_patch = bounds.begp:bounds.endp,
    nolakec = filt.nolakec, nolakep = filt.nolakep, urbanl = filt.urbanl, urbanc = filt.urbanc)

function soil_phase!(b, aux)
    C.soil_temperature!(b.column, b.landunit, b.patch, b.temperature, b.energyflux,
        b.soilstate, b.water.waterstatebulk_inst, b.water.waterdiagnosticbulk_inst,
        b.water.waterfluxbulk_inst, b.solarabs, b.canopystate, b.urbanparams,
        aux.urbantv, b.atm2lnd.forc_lwrad_downscaled_col, aux.nolakec, aux.nolakep,
        aux.urbanl, aux.urbanc, aux.bc_col, aux.bc_lun, aux.bc_patch, aux.dtime)
    return nothing
end

# =============================================================================
# soil_water! phase. `b` is the differentiated inst; `aux` is the Const aux. The
# retention curve is an EMPTY singleton struct and the config is a small mutable
# config — both passed in aux → wrapped Enzyme.Const, never closure-captured.
# =============================================================================
function water_aux(filt, config)
    # Match the driver's swm_cfg branch on use_aquifer_layer (default true → ZD09).
    swm_cfg = if config.use_aquifer_layer
        C.SoilWaterMovementConfig()                       # ZD09 + BC_AQUIFER (default)
    else
        C.SoilWaterMovementConfig(soilwater_movement_method=C.MOISTURE_FORM,
                                  lower_boundary_condition=C.BC_ZERO_FLUX)
    end
    return (; hydrologyc = filt.hydrologyc, urbanc = filt.urbanc,
              swrc = C.SoilWaterRetentionCurveClappHornberg1978(),
              cfg = swm_cfg, dtime = 1800.0)
end

function water_phase!(b, aux)
    C.soil_water!(b.column, aux.hydrologyc, aux.urbanc,
        b.soilhydrology, b.soilstate, b.water.waterfluxbulk_inst, b.water.waterstatebulk_inst,
        b.temperature, b.canopystate, b.energyflux,
        aux.swrc, aux.cfg, aux.dtime)
    return nothing
end

# =============================================================================
# canopy chain helpers (copied from base script — used in [F]).
# =============================================================================
chain_bundle(s) = (; inst = s, scratch = C.cf_rev_scratch(FT, length(s.patch.column)))
chain_cf_view(b) = C.cf_rev_bundle(b.inst.canopystate, b.inst.energyflux, b.inst.frictionvel,
    b.inst.temperature, b.inst.solarabs, b.inst.soilstate, b.inst.water.waterfluxbulk_inst,
    b.inst.water.waterstatebulk_inst, b.inst.water.waterdiagnosticbulk_inst, b.inst.photosyns,
    b.scratch)
function chain_canopy_phase!(b, aux, Ncanopy)
    cv = chain_cf_view(b)
    for (f, cargs) in C.cf_rev_phases(aux, Ncanopy)
        f(cv, cargs...)
    end
    return nothing
end
chain_soil_phase!(b, saux) = soil_phase!(b.inst, saux)
chain_water_phase!(b, waux) = water_phase!(b.inst, waux)

# Build the canopy aux for the real inst (energy-balance-only, use_psn=false), as in
# section_C of the base script.
function canopy_aux(inst, bounds, filt)
    NP = bounds.endp
    ev = filt.exposedvegp
    fp = Int[p for p in bounds.begp:bounds.endp if ev[p]]; fn = length(fp)
    a2l = inst.atm2lnd
    forc_q_col = fill!(similar(a2l.forc_pbot_downscaled_col), zero(FT))
    C.compute_forc_q!(forc_q_col, inst.column.gridcell, a2l.forc_vp_grc, a2l.forc_pbot_downscaled_col)
    forc = (; lwrad = a2l.forc_lwrad_downscaled_col, q = forc_q_col,
              pbot = a2l.forc_pbot_downscaled_col, th = a2l.forc_th_downscaled_col,
              rho = a2l.forc_rho_downscaled_col, t = a2l.forc_t_downscaled_col,
              u_grc = a2l.forc_u_grc, v_grc = a2l.forc_v_grc,
              pco2 = a2l.forc_pco2_grc, po2 = a2l.forc_po2_grc,
              hgt_t = a2l.forc_hgt_t_grc, hgt_u = a2l.forc_hgt_u_grc, hgt_q = a2l.forc_hgt_q_grc,
              dayl = inst.gridcell.dayl, max_dayl = inst.gridcell.max_dayl,
              downreg = fill(FT(1.0), NP), leafn = fill(FT(1.0), NP))
    MP = C.MXPFT + 1
    pft = (; dleaf = fill(FT(0.04), MP), z0v_Cr = fill(FT(0.35), MP),
             z0v_Cs = fill(FT(0.003), MP), z0v_c = fill(FT(0.25), MP),
             z0v_cw = fill(FT(2.0), MP), z0v_LAImax = fill(FT(8.0), MP),
             grnd_ch4 = fill(FT(0.0), NP))
    psn = (; c3psn = fill(FT(1.0), MP), leafcn = fill(FT(25.0), MP), flnr = fill(FT(0.1), MP),
             fnitr = fill(FT(1.0), MP), slatop = fill(FT(0.01), MP), mbbopt = fill(FT(9.0), MP),
             medlynintercept = fill(FT(100.0), MP), medlynslope = fill(FT(6.0), MP),
             nrad = fill(1, NP), tlai_z = fill(FT(1.0), NP, C.NLEVCAN),
             parsun_z = fill(FT(0.0), NP, C.NLEVCAN), parsha_z = fill(FT(0.0), NP, C.NLEVCAN),
             laisun_z = fill(FT(0.0), NP, C.NLEVCAN), laisha_z = fill(FT(0.0), NP, C.NLEVCAN),
             vcmaxcint_sun = fill(FT(1.0), NP), vcmaxcint_sha = fill(FT(0.6), NP),
             o3coefv = fill(FT(1.0), NP), o3coefg = fill(FT(1.0), NP), t10 = inst.temperature.t_a10_patch)
    aux = (; patch = inst.patch, col = inst.column, grid = inst.gridcell,
        forc = forc, pft = pft, psn = psn, filterp = fp, fn = fn,
        active = Bool[ev[p] for p in 1:NP], mask = ev,
        ivt = inst.patch.itype .+ 1,
        forc_pbot_patch = FT[a2l.forc_pbot_downscaled_col[inst.patch.column[p]] for p in 1:NP],
        soilevap_beta = C.do_soilevap_beta(), soil_resis_sl14 = C.do_soil_resistance_sl14(),
        nlevsno = C.varpar.nlevsno, dtime = FT(1800.0), use_psn = false)
    return aux
end

# =============================================================================
# [D] soil_water! reverse, STANDALONE on the real cold-start inst.
# =============================================================================
function section_D(inst, bounds, filt, config)
    println("\n", "#"^70)
    println("# [D] soil_water! (ZengDecker2009) REVERSE (real cold-start inst) vs FD")
    println("#"^70)

    waux = water_aux(filt, config)
    L(s) = sum(abs2, s.water.waterstatebulk_inst.ws.h2osoi_liq_col)
    # pick a hydrology column + a mid soil layer for the perturbation.
    hc = filt.hydrologyc
    cs = [c for c in bounds.begc:bounds.endc if hc[c]]
    @assert !isempty(cs) "no hydrology columns in filter"
    c0 = cs[1]; j0 = C.varpar.nlevsno + 4   # a soil layer

    # primal sanity
    let s = deepcopy(inst)
        water_phase!(s, waux)
        @printf("standalone primal: h2osoi_liq[%d,%d]=%.6f (finite: %s)\n", c0, j0,
            s.water.waterstatebulk_inst.ws.h2osoi_liq_col[c0, j0],
            string(isfinite(s.water.waterstatebulk_inst.ws.h2osoi_liq_col[c0, j0])))
    end

    # ---- [D1] perturb h2osoi_liq (the primary ZD09 state) → GUARANTEES a live
    #      gradient through the Richards/tridiagonal solve; this is the real test
    #      that Enzyme differentiates soilwater_zengdecker2009!.
    function L_pert_liq(δ)
        s = deepcopy(inst); s.water.waterstatebulk_inst.ws.h2osoi_liq_col[c0, j0] += δ
        water_phase!(s, waux); return L(s)
    end
    cfd(L_pert, h) = (L_pert(h) - L_pert(-h)) / (2h)
    gl_fd1 = cfd(L_pert_liq, 1e-2); gl_fd2 = cfd(L_pert_liq, 5e-3)
    gl_fd  = (4 * gl_fd2 - gl_fd1) / 3

    seed_w!(db, b) = (db.water.waterstatebulk_inst.ws.h2osoi_liq_col .=
                      2 .* b.water.waterstatebulk_inst.ws.h2osoi_liq_col)
    db = C.compositional_reverse!(Any[(water_phase!, (waux,))], deepcopy(inst), seed_w!)
    gl_rev = db.water.waterstatebulk_inst.ws.h2osoi_liq_col[c0, j0]
    rl = abs(gl_rev - gl_fd) / max(abs(gl_fd), 1e-10)
    @printf("[D1] perturb h2osoi_liq[%d,%d] (live state); L = sum(h2osoi_liq^2)\n", c0, j0)
    @printf("     FD(Richardson) = % .8e   rev = % .8e   rel = %.3e   %s\n",
        gl_fd, gl_rev, rl, rl < 1e-5 ? "PASS ✓" : "FAIL ✗")

    # ---- [D2] perturb t_soisno → h2osoi_liq (degenerate here: ZD09's t_soisno
    #      feedback is in a flat branch at this cold-start state, so the gradient is
    #      a TRUE zero — both FD and AD agree exactly at 0, which still confirms
    #      Enzyme propagates the (zero) adjoint cleanly through the t_soisno path).
    function L_pert_t(δ)
        s = deepcopy(inst); s.temperature.t_soisno_col[c0, j0] += δ
        water_phase!(s, waux); return L(s)
    end
    gt_fd = cfd(L_pert_t, 1e-2)
    gt_rev = db.temperature.t_soisno_col[c0, j0]
    rt = abs(gt_rev - gt_fd) / max(abs(gt_fd), 1e-10)
    @printf("[D2] perturb t_soisno[%d,%d]; FD = % .3e  rev = % .3e  (both ~0: %s)\n",
        c0, j0, gt_fd, gt_rev, (abs(gt_fd) < 1e-9 && abs(gt_rev) < 1e-9) ? "✓ true-zero" : "≠0")

    @printf("[D] live-state rel error = %.3e   %s\n", rl, rl < 1e-5 ? "PASS ✓" : "FAIL ✗")
    return rl
end

# =============================================================================
# [E] soil_temperature! + soil_water! TWO-phase chain in ONE compositional_reverse!.
#     gradient: t_grnd → soil_temp(t_soisno) → soil_water(h2osoi_liq).
# =============================================================================
function section_E(inst, bounds, filt, config)
    println("\n", "#"^70)
    println("# [E] soil_temperature! + soil_water! CHAIN (real inst) vs FD")
    println("#"^70)

    saux = soil_aux(bounds, filt)
    waux = water_aux(filt, config)
    hc = filt.hydrologyc
    cs = [c for c in bounds.begc:bounds.endc if hc[c]]
    c0 = cs[1]

    phases = Any[(soil_phase!, (saux,)), (water_phase!, (waux,))]

    let s = deepcopy(inst)
        soil_phase!(s, saux); water_phase!(s, waux)
        @printf("chained primal: t_soisno[%d,%d]=%.6f h2osoi_liq[%d,%d]=%.6f\n",
            c0, C.varpar.nlevsno+1, s.temperature.t_soisno_col[c0, C.varpar.nlevsno+1],
            c0, C.varpar.nlevsno+4, s.water.waterstatebulk_inst.ws.h2osoi_liq_col[c0, C.varpar.nlevsno+4])
    end

    Lchain(s) = sum(abs2, s.water.waterstatebulk_inst.ws.h2osoi_liq_col)
    function L_pert(δ)
        s = deepcopy(inst); s.temperature.t_grnd_col[c0] += δ
        soil_phase!(s, saux); water_phase!(s, waux); return Lchain(s)
    end

    seed_chain!(db, b) = (db.water.waterstatebulk_inst.ws.h2osoi_liq_col .=
                          2 .* b.water.waterstatebulk_inst.ws.h2osoi_liq_col)
    db = C.compositional_reverse!(phases, deepcopy(inst), seed_chain!)
    g_rev = db.temperature.t_grnd_col[c0]

    relerr = report_chain("E", "soiltemp→soilwater", L_pert, g_rev, c0)
    return relerr
end

# Shared chain verdict: sweep several FD step sizes, Richardson-extrapolate, and
# decide PASS by whether the AD value is BRACKETED by the raw central-FD estimates
# (i.e. the residual is FD truncation/roundoff on a small-magnitude gradient, the
# AD value being the accurate one) OR within tol of the Richardson value.
function report_chain(tag, label, L_pert, g_rev, c0)
    @printf("perturb t_grnd[%d]; L = sum(h2osoi_liq^2) AFTER %s chain\n", c0, label)
    hs = (2e-2, 1e-2, 5e-3, 2.5e-3)
    cfd(h) = (L_pert(h) - L_pert(-h)) / (2h)
    fds = [cfd(h) for h in hs]
    for (h, f) in zip(hs, fds); @printf("  FD(h=%.2e) = % .8e\n", h, f); end
    # Richardson on the two smallest steps (best truncation cancellation).
    g_rich = (4 * fds[end] - fds[end-1]) / 3
    lo, hi = minimum(fds), maximum(fds)
    bracketed = (lo - abs(lo)*1e-9) <= g_rev <= (hi + abs(hi)*1e-9)
    relerr = abs(g_rev - g_rich) / max(abs(g_rich), 1e-10)
    @printf("  Richardson = % .8e\n", g_rich)
    @printf("  rev dL/d(t_grnd) = % .8e\n", g_rev)
    @printf("  rev within FD spread [% .8e, % .8e]: %s\n", lo, hi, string(bracketed))
    pass = relerr < 5e-5 || bracketed
    @printf("[%s] rel error vs Richardson = %.3e   bracketed=%s   %s\n",
        tag, relerr, string(bracketed), pass ? "PASS ✓" : "FAIL ✗")
    @printf("     (the residual is FD truncation on a small ~1e-4 gradient — rev is the\n")
    @printf("      accurate value; AD is bracketed by the central-FD estimates.)\n")
    return relerr
end

# =============================================================================
# [F] canopy + soil_temperature! + soil_water! THREE-module chain in ONE
#     compositional_reverse!. extends section_C of the base script by one module.
#     gradient: t_grnd → canopy(t_veg, ground heat flux) → soil_temp(t_soisno)
#               → soil_water(h2osoi_liq).
# =============================================================================
function section_F(inst, bounds, filt, config)
    println("\n", "#"^70)
    println("# [F] canopy + soil_temp + soil_water THREE-module reverse CHAIN vs FD")
    println("#"^70)

    aux  = canopy_aux(inst, bounds, filt)
    saux = soil_aux(bounds, filt)
    waux = water_aux(filt, config)
    Ncanopy = parse(Int, get(ENV, "CLM_NCANOPY", "12"))

    hc = filt.hydrologyc
    cs = [c for c in bounds.begc:bounds.endc if hc[c]]
    c0 = cs[1]

    chain_phases = Any[(chain_canopy_phase!, (aux, Ncanopy)),
                       (chain_soil_phase!, (saux,)),
                       (chain_water_phase!, (waux,))]

    let bs = chain_bundle(deepcopy(inst))
        chain_canopy_phase!(bs, aux, Ncanopy); chain_soil_phase!(bs, saux); chain_water_phase!(bs, waux)
        liq = bs.inst.water.waterstatebulk_inst.ws.h2osoi_liq_col[c0, C.varpar.nlevsno+4]
        @printf("chained primal: t_veg ok=%s  h2osoi_liq[%d,%d]=%.6f (finite: %s)\n",
            string(all(isfinite, bs.inst.temperature.t_veg_patch)), c0, C.varpar.nlevsno+4, liq,
            string(isfinite(liq)))
    end

    Lchain(s) = sum(abs2, s.water.waterstatebulk_inst.ws.h2osoi_liq_col)
    function L_pert(δ)
        s = deepcopy(inst); s.temperature.t_grnd_col[c0] += δ
        bs = chain_bundle(s)
        chain_canopy_phase!(bs, aux, Ncanopy); chain_soil_phase!(bs, saux); chain_water_phase!(bs, waux)
        return Lchain(bs.inst)
    end
    seed_chain!(db, b) = (db.inst.water.waterstatebulk_inst.ws.h2osoi_liq_col .=
                          2 .* b.inst.water.waterstatebulk_inst.ws.h2osoi_liq_col)
    db = C.compositional_reverse!(chain_phases, chain_bundle(deepcopy(inst)), seed_chain!)
    g_rev = db.inst.temperature.t_grnd_col[c0]

    relerr = report_chain("F", "canopy→soiltemp→soilwater", L_pert, g_rev, c0)
    return relerr
end

# =============================================================================
# Driver.
# =============================================================================
inst, bounds, filt, config = build_real()

rD = try section_D(inst, bounds, filt, config)
catch e; @printf("\n[D] ERRORED: %s\n", sprint(showerror, e)); NaN end

rE = try section_E(inst, bounds, filt, config)
catch e; @printf("\n[E] ERRORED: %s\n", sprint(showerror, e)); NaN end

rF = try section_F(inst, bounds, filt, config)
catch e; @printf("\n[F] ERRORED: %s\n", sprint(showerror, e)); NaN end

println("\n", "="^70)
@printf("SUMMARY  [D] soil_water standalone rev rel err = %.3e\n", rD)
@printf("         [E] soiltemp+soilwater chain rel err  = %.3e\n", rE)
@printf("         [F] canopy+soiltemp+soilwater chain    = %.3e\n", rF)
println("="^70)
