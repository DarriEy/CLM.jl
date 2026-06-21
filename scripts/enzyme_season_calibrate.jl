# =============================================================================
# SEASON-GPP CALIBRATION — recover a photosynthesis parameter by gradient descent through the
# multi-timestep reverse chain. THE PAYOFF of the season-gradient stack: the engine becomes a
# calibration tool against a season-integrated GPP target.
#
# Builds on enzyme_season_gpp.jl (warmed photosynthesizing inst; accumulate the LIVE per-leaf
# rates psnsun·laisun+psnsha·laisha; absorbed PAR is the confirmed-live driver). Here a SCALAR
# calibration parameter θ lives in the bundle (Duplicated) and multiplies a photosynthesis driver
# inside a custom psn phase, so d(season GPP)/dθ flows through the reverse. Loss = (GPP(θ) −
# GPP_target)²; seed db.acc = 2(GPP−target); multistep_reverse! returns db.θ = d(Loss)/dθ;
# gradient descent recovers θ_true. A leading FD probe picks the live driver (vcmax scale if
# photosynthesis! uses vcmaxcint, else the PAR / light-use-efficiency scale).
#   CLM_NSTEPS=2 CLM_NCANOPY=4 julia +1.10 --project=/tmp/clm_jl10_<id> scripts/enzyme_season_calibrate.jl
# =============================================================================
using CLM, Enzyme, Printf
const C = CLM
const FSURDAT  = get(ENV, "CLM_FSURDAT",  "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/surfdata_clm.nc")
const PARAMFILE = get(ENV, "CLM_PARAMFILE", "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/clm5_params.nc")
const DT = 1800.0

function build_growing_season()
    (inst, bounds, filt, tm) = C.clm_initialize!(; fsurdat=FSURDAT, paramfile=PARAMFILE)
    config = C.CLMDriverConfig(); fia = C.clump_filter_inactive_and_active
    calday = 172.5; (declin, _) = C.compute_orbital(calday); nextsw = calday + DT/C.SECSPDAY
    C._setup_calib_forcing!(inst.atm2lnd, 295.0, bounds.endg)
    C.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
    C._init_calib_soil_moisture!(inst, bounds)
    C.interp_monthly_veg!(inst.satellite_phenology; kmo=7, kda=15)
    cs = inst.canopystate; wdb = inst.water.waterdiagnosticbulk_inst
    C.satellite_phenology!(inst.satellite_phenology, cs, wdb, inst.patch, filt.nolakep, bounds.begp:bounds.endp)
    for p in bounds.begp:bounds.endp; cs.frac_veg_nosno_patch[p] = cs.frac_veg_nosno_alb_patch[p]; end
    C.set_exposedvegp_filter!(filt, bounds, cs.frac_veg_nosno_patch)
    for n in 1:4
        C.clm_drv!(config, inst, filt, fia, bounds, true, nextsw, declin, declin,
            C.ORB_OBLIQR_DEFAULT, false, false, "", false; nstep=n, is_first_step=(n==1),
            is_beg_curr_day=(n==1), dtime=DT, mon=7, day=15, photosyns=inst.photosyns)
    end
    return inst, bounds, filt, config
end

inst, bounds, filt, config = build_growing_season()
const NSTEPS = parse(Int, get(ENV, "CLM_NSTEPS", "2"))
const NCAN   = parse(Int, get(ENV, "CLM_NCANOPY", "4"))
ev = filt.exposedvegp
psnp = Int[p for p in bounds.begp:bounds.endp if ev[p] && inst.photosyns.fpsn_patch[p] > 1e-6]
fp = first(psnp)
@printf("warmed inst: %d photosynthesizing patch(es); fp=%d\n", length(psnp), fp)

base = C.forcingset_aux(inst)
schedule = [(; lwrad = base.lwrad .* (1 + 0.15*sinpi(2*(k-1)/NSTEPS)),
              t  = base.t  .+ 4.0*sinpi(2*(k-1)/NSTEPS),
              th = base.th .+ 4.0*sinpi(2*(k-1)/NSTEPS), rho = base.rho) for k in 1:NSTEPS]
caux = C.canopy_rev_aux(inst, bounds, filt; use_psn=true)

# DRIVER: "vcmax" scales vcmaxcint (live, inst.surfalb); "par" scales absorbed PAR (live,
# inst.solarabs). Chosen after the FD probe below (prefer vcmax — the physical Rubisco param).
mutable struct CalCfg; mode::Symbol; end
const CAL = CalCfg(:vcmax)
function psn_cal!(b, aux, phase)
    sc = b.scratch; tp = b.inst.temperature; ps = b.inst.photosyns; pn = aux.psn; θ = b.theta[1]
    vsun = b.inst.surfalb.vcmaxcintsun_patch; vsha = b.inst.surfalb.vcmaxcintsha_patch
    psun = b.inst.solarabs.parsun_z_patch;    psha = b.inst.solarabs.parsha_z_patch
    if CAL.mode === :vcmax
        vcm  = phase == "sun" ? θ .* vsun : θ .* vsha
        parz = phase == "sun" ? psun : psha
    else
        vcm  = phase == "sun" ? vsun : vsha
        parz = phase == "sun" ? θ .* psun : θ .* psha
    end
    laiz = phase == "sun" ? pn.laisun_z : pn.laisha_z
    C.photosynthesis!(ps, sc.svpts, sc.eah, sc.o2_arr, sc.co2_arr, sc.rb,
        b.inst.energyflux.btran_patch, sc.dayl_factor, aux.forc.leafn, aux.forc_pbot_patch,
        tp.t_veg_patch, pn.t10, tp.thm_patch, pn.nrad, pn.tlai_z, b.inst.canopystate.tlai_patch,
        parz, laiz, vcm, pn.o3coefv, pn.o3coefg, pn.c3psn, pn.leafcn, pn.flnr, pn.fnitr,
        pn.slatop, pn.mbbopt, pn.medlynintercept, pn.medlynslope, aux.ivt, aux.patch.column,
        aux.mask, 1:length(sc.air), phase)
    return nothing
end
auxG = (; dt = DT, expp = psnp, laisun = caux.psn.laisun_z, laisha = caux.psn.laisha_z)
function accum_gpp!(b, aux)
    ps = b.inst.photosyns
    @inbounds for p in aux.expp
        b.acc[p] += (ps.psnsun_patch[p]*aux.laisun[p,1] + ps.psnsha_patch[p]*aux.laisha[p,1]) * aux.dt
    end
    return nothing
end

np = bounds.endp
bundle(θ) = (; C.driver_rev_bundle(deepcopy(inst))..., acc = zeros(np), theta = [θ])
steps = [vcat(phs, Any[(psn_cal!, (caux, "sun")), (psn_cal!, (caux, "sha")), (accum_gpp!, (auxG,))])
         for phs in C.forced_driver_steps(bounds, filt, config, schedule; canopy_aux=caux, n_canopy=NCAN)]
GPP(θ) = (b = bundle(θ); for ph in steps, (f,ca) in ph; f(b,ca...); end; sum(@view b.acc[psnp]))

# ---- pick the live driver via an FD probe (prefer vcmax) ----
g_vcmax = (CAL.mode = :vcmax; (GPP(1.0+1e-3)-GPP(1.0-1e-3))/2e-3)
g_par   = (CAL.mode = :par;   (GPP(1.0+1e-3)-GPP(1.0-1e-3))/2e-3)
@printf("FD d(GPP)/dθ:  vcmax-scale=% .4e   par-scale=% .4e\n", g_vcmax, g_par)
CAL.mode = abs(g_vcmax) > 1e-4 ? :vcmax : :par
@printf("calibration driver = %s scale\n", CAL.mode)

# ---- synthetic target + gradient-based recovery ----
const θ_true = 1.0
const G_target = GPP(θ_true)
Loss(θ) = (GPP(θ) - G_target)^2
seedL!(db, b) = (g = sum(@view b.acc[psnp]); v = 2*(g - G_target); db.acc .= 0.0; for p in psnp; db.acc[p] = v; end; nothing)
gradθ(θ) = C.multistep_reverse!(steps, bundle(θ), seedL!).theta[1]

# grad check at a perturbed θ (CLM_THETA0 to start above/below the optimum)
θ0 = parse(Float64, get(ENV, "CLM_THETA0", "1.4"))
g_rev = gradθ(θ0); g_fd = (Loss(θ0+1e-3)-Loss(θ0-1e-3))/2e-3
@printf("\ngrad check @ θ=%.2f:  reverse=% .6e  FD=% .6e  rel=%.2e\n", θ0, g_rev, g_fd, abs(g_rev-g_fd)/max(abs(g_fd),1e-12))

# gradient descent. Loss(θ)=(G'·(θ−θ_true))² near the optimum (G'=d(GPP)/dθ), so the Hessian is
# ≈2·G'² and lr=1/(2·G'²) is a Newton step → converges in a few iterations.
g_drv = CAL.mode === :vcmax ? g_vcmax : g_par
lr = 1.0 / (2 * g_drv^2)
@printf("\ngradient descent (θ_true=%.3f, θ0=%.3f, GPP_target=%.4e, lr=%.2e):\n", θ_true, θ0, G_target, lr)
global θ = θ0
for it in 1:12
    global θ
    g = gradθ(θ); θ -= lr * g
    @printf("  it %2d: θ=%.6f  Loss=%.4e  |θ-θ_true|=%.2e\n", it, θ, Loss(θ), abs(θ-θ_true))
    abs(θ - θ_true) < 1e-4 && break
end
println("="^78)
ok_grad = abs(g_rev-g_fd)/max(abs(g_fd),1e-12) < 1e-3
ok_rec  = abs(θ - θ_true) < 1e-3
@printf("SEASON-GPP CALIBRATION (%s scale): grad-check %s, recovered θ=%.6f (true %.3f) %s\n",
    CAL.mode, ok_grad ? "PASS" : "FAIL", θ, θ_true, ok_rec ? "PASS ✓" : "FAIL ✗")
@printf("%s\n", (ok_grad && ok_rec) ? "PASS ✓ (parameter recovered by gradient descent through the season reverse chain)" : "FAIL ✗")
println("="^78)
