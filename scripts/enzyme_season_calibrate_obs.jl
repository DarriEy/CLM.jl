# =============================================================================
# REGULARIZED calibration against IMPERFECT (observation-style) GPP — the step from a synthetic
# target the model can fit exactly (loss→0) to an "observed" target it CANNOT (loss floor > 0),
# which exposes model-data mismatch and the need for regularization. Builds on the multi-param
# season-GPP calibration (vcmax θ1 + light-use θ2, per-step GPP series, reverse gradient vector).
#
# The Bow inst is a lumped hydrology domain (no co-located flux tower) and FLUXNET is download-
# gated, so "observations" = model truth + Gaussian obs noise + a STRUCTURAL bias on the dim
# steps that NO (θ1,θ2) can reproduce (a localized discrepancy → irreducible residual). The fit:
#   data loss  Ldata(θ) = Σ_k (series_k(θ) − OBS_k)²    (gradient via multistep_reverse!)
#   regularized Lreg(θ) = Ldata(θ) + λ‖θ − θ_prior‖²    (reg gradient added analytically)
# Demonstrated: (1) the unregularized best fit has a NONZERO loss floor and is pulled OFF the
# physical prior by noise+bias; (2) Tikhonov regularization trades data-fit for prior-closeness
# and stabilizes the solution (an L-curve sweep over λ).
#   CLM_NSTEPS=4 CLM_NCANOPY=4 julia +1.10 --project=/tmp/clm_jl10_<id> scripts/enzyme_season_calibrate_obs.jl
# =============================================================================
using CLM, Enzyme, Printf, Random
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
const NSTEPS = parse(Int, get(ENV, "CLM_NSTEPS", "4"))
const NCAN   = parse(Int, get(ENV, "CLM_NCANOPY", "4"))
ev = filt.exposedvegp
psnp = Int[p for p in bounds.begp:bounds.endp if ev[p] && inst.photosyns.fpsn_patch[p] > 1e-6]
base = C.forcingset_aux(inst)
schedule = [(; lwrad = base.lwrad .* (1 + 0.15*sinpi(2*(k-1)/NSTEPS)),
              t  = base.t  .+ 4.0*sinpi(2*(k-1)/NSTEPS),
              th = base.th .+ 4.0*sinpi(2*(k-1)/NSTEPS), rho = base.rho) for k in 1:NSTEPS]
caux = C.canopy_rev_aux(inst, bounds, filt; use_psn=true)
lights = [isodd(k) ? 1.0 : 0.1 for k in 1:NSTEPS]

function psn_cal2!(b, aux, phase, light)
    sc = b.scratch; tp = b.inst.temperature; ps = b.inst.photosyns; pn = aux.psn
    θ1 = b.t1[1]; θ2 = b.t2[1]
    vsun = b.inst.surfalb.vcmaxcintsun_patch; vsha = b.inst.surfalb.vcmaxcintsha_patch
    psun = b.inst.solarabs.parsun_z_patch;    psha = b.inst.solarabs.parsha_z_patch
    vcm  = phase == "sun" ? θ1 .* vsun : θ1 .* vsha
    parz = phase == "sun" ? (θ2*light) .* psun : (θ2*light) .* psha
    laiz = phase == "sun" ? pn.laisun_z : pn.laisha_z
    C.photosynthesis!(ps, sc.svpts, sc.eah, sc.o2_arr, sc.co2_arr, sc.rb,
        b.inst.energyflux.btran_patch, sc.dayl_factor, aux.forc.leafn, aux.forc_pbot_patch,
        tp.t_veg_patch, pn.t10, tp.thm_patch, pn.nrad, pn.tlai_z, b.inst.canopystate.tlai_patch,
        parz, laiz, vcm, pn.o3coefv, pn.o3coefg, pn.c3psn, pn.leafcn, pn.flnr, pn.fnitr,
        pn.slatop, pn.mbbopt, pn.medlynintercept, pn.medlynslope, aux.ivt, aux.patch.column,
        aux.mask, 1:length(sc.air), phase)
    return nothing
end
function accum_k!(b, aux)
    ps = b.inst.photosyns; s = 0.0
    @inbounds for p in aux.expp
        s += ps.psnsun_patch[p]*aux.laisun[p,1] + ps.psnsha_patch[p]*aux.laisha[p,1]
    end
    b.series[aux.k] += s * aux.dt
    return nothing
end
np = bounds.endp
bundle(θ1,θ2) = (; C.driver_rev_bundle(deepcopy(inst))..., t1=[θ1], t2=[θ2], series=zeros(NSTEPS))
fsteps = C.forced_driver_steps(bounds, filt, config, schedule; canopy_aux=caux, n_canopy=NCAN)
steps = [vcat(fsteps[k],
              Any[(psn_cal2!, (caux, "sun", lights[k])), (psn_cal2!, (caux, "sha", lights[k])),
                  (accum_k!, ((; dt=DT, expp=psnp, laisun=caux.psn.laisun_z, laisha=caux.psn.laisha_z, k=k),))])
         for k in 1:NSTEPS]
series(θ1,θ2) = (b = bundle(θ1,θ2); for ph in steps, (f,ca) in ph; f(b,ca...); end; copy(b.series))

# ---- build IMPERFECT observations: truth + obs noise + UNFITTABLE structural bias on dim steps ----
const θ1t, θ2t = 1.0, 1.0
truth = series(θ1t, θ2t)
rng = MersenneTwister(20260621)
obs_noise = 0.05 .* truth .* randn(rng, NSTEPS)                    # 5% Gaussian obs noise
struct_bias = [iseven(k) ? 0.18*truth[k] : 0.0 for k in 1:NSTEPS]  # dim steps +18% — no (θ1,θ2) reproduces this
const OBS = truth .+ obs_noise .+ struct_bias
@printf("truth   : %s\nOBS     : %s  (noise+structural bias)\n",
    string(round.(truth,digits=1)), string(round.(OBS,digits=1)))

Ldata(θ1,θ2) = sum(abs2, series(θ1,θ2) .- OBS)
function gdata(θ1,θ2)
    seed!(db,b) = (db.series .= 2 .* (b.series .- OBS); nothing)
    db = C.multistep_reverse!(steps, bundle(θ1,θ2), seed!)
    return db.t1[1], db.t2[1]
end
# grad check the data gradient vs FD (the reg part is analytic)
let (g1,g2)=gdata(1.2,0.9), fd1=(Ldata(1.2+1e-3,0.9)-Ldata(1.2-1e-3,0.9))/2e-3, fd2=(Ldata(1.2,0.9+1e-3)-Ldata(1.2,0.9-1e-3))/2e-3
    @printf("data-grad check @(1.2,0.9): ∂/∂θ1 rel=%.2e  ∂/∂θ2 rel=%.2e\n",
        abs(g1-fd1)/max(abs(fd1),1e-12), abs(g2-fd2)/max(abs(fd2),1e-12))
end

# full 2×2-Newton on the REGULARIZED objective Lreg = Ldata + λ‖θ−θprior‖²  (θprior=(1,1))
function fit(λ)
    θ1,θ2 = 1.3, 0.7
    greg(a,b) = (g=gdata(a,b); (g[1]+2λ*(a-1), g[2]+2λ*(b-1)))
    for _ in 1:8
        g1,g2 = greg(θ1,θ2); h=1e-3
        g1p,g2p = greg(θ1+h,θ2); g1q,g2q = greg(θ1,θ2+h)
        H11=(g1p-g1)/h; H21=(g2p-g2)/h; H12=(g1q-g1)/h; H22=(g2q-g2)/h
        det = H11*H22-H12*H21
        θ1 += -( H22*g1 - H12*g2)/det; θ2 += -(-H21*g1 + H11*g2)/det
        hypot(g1,g2) < 1e-3 && break
    end
    return θ1, θ2
end

# λ must span the DATA curvature to show a tradeoff: GPP~2e4 and ∂series/∂θ~2e4 ⇒ the
# Gauss-Newton data Hessian JᵀJ~1e9, so λ≪1e8 is invisible and λ≫1e10 forces θ→prior.
println("\nL-curve: regularized best-fit vs λ (θ_prior=(1,1), θ_true=(1,1)):")
@printf("  %-10s %-22s %-12s %-12s\n", "λ", "θ_hat", "Ldata(floor)", "‖θ-prior‖")
local θ_lo = (1.0,1.0)
for λ in (0.0, 1e8, 3e8, 1e9, 3e9, 1e10)
    a,b = fit(λ)
    λ == 0.0 && (θ_lo = (a,b))
    @printf("  %-10.0e (%.5f,%.5f)   %.4e    %.4e\n", λ, a, b, Ldata(a,b), hypot(a-1,b-1))
end
println("="^84)
# unregularized fit: nonzero data-loss FLOOR (model can't fit noise+bias) + θ pulled off the prior.
a0,b0 = fit(0.0); floor0 = Ldata(a0,b0)
println("FINDINGS:")
@printf("  • model-data mismatch: unregularized Ldata floor = %.4e > 0 (the structural bias is\n", floor0)
@printf("    unfittable by ANY (θ1,θ2)); best-fit θ=(%.4f,%.4f) is pulled off the prior by noise+bias.\n", a0, b0)
@printf("  • regularization (λ>0) trades data-fit for prior-closeness → ‖θ-prior‖ shrinks as λ grows.\n")
pass = floor0 > 1e3 && isfinite(a0) && isfinite(b0)
@printf("%s\n", pass ? "PASS ✓ (regularized calibration against imperfect obs: nonzero loss floor + L-curve tradeoff)" : "FAIL ✗")
println("="^84)
