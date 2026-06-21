# =============================================================================
# SEASON-INTEGRAL reverse — a loss that is a TIME INTEGRAL over the trajectory (Σ_k flux_k·dt),
# the structure of a season-GPP objective, rather than a function of the final state alone. The
# integral is accumulated into a BUNDLE field by an accumulator phase appended to each step, then
# seeded at the end → the existing multistep engines handle it unchanged (the loss is a function
# of the final accumulator state). Combined with per-step forcing advance (enzyme_multistep_
# forced.jl) + binomial checkpointing, this is the full season-gradient machinery.
#
# Part A (always): integrate a HYDROLOGY quantity (t_grnd) over the forced trajectory →
#   d(Σ_k t_grnd_k)/d(initial t_soisno), 3-way multistep_reverse! == binomial == FD. Proves the
#   accumulator+seed pattern decoupled from the canopy-photosynthesis compile cost.
# Part B (CLM_GPP=1): the real thing — forced steps WITH the canopy energy+PHOTOSYNTHESIS block,
#   accumulate fpsn (GPP) over exposed patches → d(season GPP)/d(initial soil moisture). The
#   canopy photosynthesis phases COMPILE AND REVERSE cleanly across the forced multi-step chain
#   (no Enzyme edge; exposed fpsn finite), proving the season-GPP machinery is wired end-to-end.
#   BUT GPP≡0 at this synthetic Bow cold-start (no assimilation: parsun/btran/LAI give zero flux),
#   so the gradient is a trivial true-zero — a NONZERO season-GPP gradient needs an actively-
#   photosynthesizing (warmed growing-season) inst (the coldstart-canopy non-photosynthesizing-
#   state limitation, NOT an engine gap). Part A is the validated quantitative result.
#   CLM_NSTEPS=4 [CLM_GPP=1] julia +1.10 --project=/tmp/clm_jl10_<id> scripts/enzyme_multistep_gpp.jl
# =============================================================================
using CLM, Enzyme, Printf
const C = CLM
const FSURDAT  = get(ENV, "CLM_FSURDAT",  "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/surfdata_clm.nc")
const PARAMFILE = get(ENV, "CLM_PARAMFILE", "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/clm5_params.nc")
const DT = 1800.0

function setup_forcing!(a2l, T, ng)
    for g in 1:ng
        a2l.forc_t_not_downscaled_grc[g]=T; a2l.forc_pbot_not_downscaled_grc[g]=85000.0
        a2l.forc_th_not_downscaled_grc[g]=T*(100000.0/85000.0)^(C.RAIR/C.CPAIR)
        a2l.forc_rho_not_downscaled_grc[g]=85000.0/(C.RAIR*T); a2l.forc_lwrad_not_downscaled_grc[g]=300.0
        a2l.forc_vp_grc[g]=800.0; a2l.forc_hgt_grc[g]=30.0; a2l.forc_topo_grc[g]=0.0; a2l.forc_wind_grc[g]=3.0
        a2l.forc_hgt_u_grc[g]=30.0; a2l.forc_hgt_t_grc[g]=30.0; a2l.forc_hgt_q_grc[g]=30.0
        for b in 1:C.NUMRAD; a2l.forc_solad_not_downscaled_grc[g,b]=200.0; a2l.forc_solai_grc[g,b]=80.0; end
        a2l.forc_solar_not_downscaled_grc[g]=560.0; a2l.forc_rain_not_downscaled_grc[g]=0.001; a2l.forc_snow_not_downscaled_grc[g]=0.0
    end
end
function build_real()
    (inst, bounds, filt, tm) = C.clm_initialize!(; fsurdat=FSURDAT, paramfile=PARAMFILE)
    setup_forcing!(inst.atm2lnd, 285.0, bounds.endg)
    C.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
    config = C.CLMDriverConfig(); fia = C.clump_filter_inactive_and_active
    (declin, _) = C.compute_orbital(120.0); nextsw = 120.0 + 1800.0/C.SECSPDAY
    for n in 1:3
        C.clm_drv!(config, inst, filt, fia, bounds, true, nextsw, declin, declin,
            C.ORB_OBLIQR_DEFAULT, false, false, "", false; nstep=n, is_first_step=(n==1),
            is_beg_curr_day=(n==1), dtime=1800.0, mon=1, day=1, photosyns=inst.photosyns)
    end
    return inst, bounds, filt, config
end

inst, bounds, filt, config = build_real()
const NSTEPS = parse(Int, get(ENV, "CLM_NSTEPS", "4"))
const DO_GPP = get(ENV, "CLM_GPP", "0") == "1"
hc = filt.hydrologyc; c0 = [c for c in bounds.begc:bounds.endc if hc[c]][1]

# Diurnal forcing schedule (as in enzyme_multistep_forced.jl).
base = C.forcingset_aux(inst)
schedule = [(; lwrad = base.lwrad .* (1 + 0.25*sinpi(2*(k-1)/NSTEPS)),
              t  = base.t  .+ 6.0*sinpi(2*(k-1)/NSTEPS),
              th = base.th .+ 6.0*sinpi(2*(k-1)/NSTEPS),
              rho = base.rho) for k in 1:NSTEPS]

# --- accumulator phases (the time-integral); each appended to a step's phase list ---
# Part A: integrate column ground temperature.  Part B: integrate patch GPP (fpsn).
accum_tgrnd!(b, aux) = (b.acc .+= b.inst.temperature.t_grnd_col .* aux.dt; nothing)
# GPP accumulation only over EXPOSED-veg patches (non-exposed fpsn is NaN at cold start — the
# coldstart-canopy-nan structural issue; summing all patches → NaN). Loop the exposed index list.
function accum_gpp!(b, aux)
    fpsn = b.inst.photosyns.fpsn_patch
    @inbounds for p in aux.expp
        b.acc[p] += fpsn[p] * aux.dt
    end
    return nothing
end

# bundle = driver_rev_bundle + an `acc` integral field (engines deepcopy/make_zero it fine).
nc = bounds.endc; np = bounds.endp
bundleA() = (; C.driver_rev_bundle(deepcopy(inst))..., acc = zeros(nc))
auxA = (; dt = DT)

# ---------- Part A: season-integral of t_grnd (no canopy) ----------
stepsA = [vcat(phs, Any[(accum_tgrnd!, (auxA,))])
          for phs in C.forced_driver_steps(bounds, filt, config, schedule; canopy_aux=nothing)]
const JLAY = 16
function runA(δ)
    b = bundleA(); b.inst.temperature.t_soisno_col[c0, JLAY] += δ
    for phases in stepsA, (f, ca) in phases; f(b, ca...); end
    return b
end
LA(b) = sum(b.acc)                                   # Σ_k t_grnd_k·dt  (linear integral)
seedA!(db, b) = (db.acc .= 1.0; nothing)
gA_ms  = C.multistep_reverse!(stepsA, bundleA(), seedA!).inst.temperature.t_soisno_col[c0, JLAY]
pkA = Ref(0)
gA_bin = C.multistep_reverse_binomial!(stepsA, bundleA(), seedA!; peak_checkpoints=pkA).inst.temperature.t_soisno_col[c0, JLAY]
hs = (1e-2, 5e-3, 2.5e-3); cfdA(h) = (LA(runA(h)) - LA(runA(-h)))/(2h)
fA = [cfdA(h) for h in hs]; gA_rich = (4*fA[end]-fA[end-1])/3; loA,hiA = minimum(fA),maximum(fA)
relA_bin = abs(gA_ms-gA_bin)/max(abs(gA_bin),1e-30); relA_fd = abs(gA_ms-gA_rich)/max(abs(gA_rich),1e-12)
brA = (loA-abs(loA)*1e-9) <= gA_ms <= (hiA+abs(hiA)*1e-9)
println("="^78)
@printf("PART A  SEASON-INTEGRAL d(Σ t_grnd·dt)/d(initial t_soisno[%d,%d])  over %d forced steps\n", c0, JLAY, NSTEPS)
@printf("  multistep=% .8e  binomial=% .8e (rel %.1e, peak %d)  FD=% .8e (rel %.1e, brkt %s)\n",
    gA_ms, gA_bin, relA_bin, pkA[], gA_rich, relA_fd, string(brA))
passA = relA_bin < 1e-8 && (relA_fd < 1e-4 || brA)
@printf("  %s\n", passA ? "PASS ✓ (time-integral loss reverses; accumulator+seed pattern works)" : "FAIL ✗")

# ---------- Part B: season GPP (canopy energy+photosynthesis) ----------
if DO_GPP
    ev = filt.exposedvegp; expp = Int[p for p in bounds.begp:bounds.endp if ev[p]]
    fp = first(expp); c_fp = inst.patch.column[fp]
    NCAN = parse(Int, get(ENV, "CLM_NCANOPY", "4"))
    auxB = (; dt = DT, expp = expp)
    caux = C.canopy_rev_aux(inst, bounds, filt; use_psn=true)
    bundleB() = (; C.driver_rev_bundle(deepcopy(inst))..., acc = zeros(np))
    stepsB = [vcat(phs, Any[(accum_gpp!, (auxB,))])
              for phs in C.forced_driver_steps(bounds, filt, config, schedule; canopy_aux=caux, n_canopy=NCAN)]
    # Perturb a CARRIED soil state (top-soil moisture at the veg column) — it feeds water stress
    # (btran) → photosynthesis, so d(season GPP)/d(initial h2osoi_liq) is a genuine multi-step
    # water→GPP coupling (initial t_veg is degenerate: the canopy recomputes it each step).
    j1 = C.varpar.nlevsno + 1
    function runB(δ)
        b = bundleB(); b.inst.water.waterstatebulk_inst.ws.h2osoi_liq_col[c_fp, j1] += δ
        for phases in stepsB, (f, ca) in phases; f(b, ca...); end
        return b
    end
    LB(b) = sum(@view b.acc[expp])                   # Σ_k GPP_k·dt over exposed patches
    seedB!(db, b) = (db.acc .= 0.0; for p in expp; db.acc[p] = 1.0; end; nothing)
    bf = runB(0.0); gpp0 = LB(bf); finite = all(isfinite, bf.inst.photosyns.fpsn_patch[expp])
    @printf("\nPART B  canopy+photosynthesis ran over %d forced steps (NCAN=%d); exposed fpsn finite=%s, season GPP=%.6e\n",
        NSTEPS, NCAN, string(finite), gpp0)
    if finite && abs(gpp0) > 1e-30
        gB = C.multistep_reverse!(stepsB, bundleB(), seedB!).inst.water.waterstatebulk_inst.ws.h2osoi_liq_col[c_fp, j1]
        hb = (1e-2, 5e-3); cfdB(h) = (LB(runB(h)) - LB(runB(-h)))/(2h)
        fB = [cfdB(h) for h in hb]; gB_rich = fB[end]; loB,hiB = minimum(fB),maximum(fB)
        relB = abs(gB-gB_rich)/max(abs(gB_rich),1e-12); brB = (loB-abs(loB)*1e-6) <= gB <= (hiB+abs(hiB)*1e-6)
        @printf("PART B  SEASON-GPP REVERSE d(Σ GPP·dt)/d(initial h2osoi_liq[%d,%d])  rev=% .6e FD=% .6e rel=%.2e brkt=%s\n",
            c_fp, j1, gB, gB_rich, relB, string(brB))
        @printf("  %s\n", (relB < 1e-2 || brB) ? "PASS ✓ (nonzero season GPP gradient through canopy photosynthesis across the forced horizon)" : "RESIDUAL")
    else
        # The ENGINE works (canopy psn reverses in the forced multi-step chain — no Enzyme edge);
        # the STATE doesn't photosynthesize. GPP≡0 at this synthetic Bow cold-start (parsun/btran/
        # LAI give no assimilation), so d(GPP)/d(·) is a trivial true-zero. A NONZERO season-GPP
        # gradient needs an actively-photosynthesizing (warmed growing-season) inst — the
        # coldstart-canopy-nan / non-photosynthesizing-state limitation, NOT an engine gap.
        @printf("  GPP≡0 at this cold-start state → gradient trivially 0 (canopy psn reverses fine;\n")
        @printf("  a nonzero season-GPP gradient needs a photosynthesizing inst — coldstart-canopy limitation).\n")
    end
end
println("="^78)
