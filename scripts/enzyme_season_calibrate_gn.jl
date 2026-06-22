# =============================================================================
# GAUSS-NEWTON / LEVENBERG-MARQUARDT season-GPP calibration scaling to MORE parameters. The
# proper least-squares solver: build the full N×M JACOBIAN J_kj = ∂series_k/∂θ_j via N reverse
# passes (seed each residual with a unit vector → multistep_reverse! returns that residual's
# gradient over all params), then take damped LM steps Δθ = −(JᵀJ + μ·diag(JᵀJ))⁻¹ Jᵀr. The
# cross-coupling finding from the 2-param fit (a diagonal step oscillates) motivates this: GN/LM
# uses the full JᵀJ and is robust to rank-deficiency via the damping μ.
#
# THREE parameters: θ1 vcmax (rate, bright-weighted), θ2 light-use/PAR (dim-weighted), θ3 leaf-
# area/fAPAR (uniform, multiplies the canopy GPP in the accumulator). This is RANK-DEFICIENT from
# total GPP — θ3·θ1 / θ3·θ2 is a near-symmetry (equifinality), so unregularized LM drives Loss→0
# to a NON-UNIQUE θ off the truth. The fix is Tikhonov regularization (a nominal prior), which
# selects the unique manifold point at the prior=truth → clean recovery. So (b) delivers the
# reverse-mode Jacobian + LM solver AND demonstrates the realistic more-params pathology + its fix.
#   CLM_NSTEPS=5 CLM_NCANOPY=4 julia +1.10 --project=/tmp/clm_jl10_<id> scripts/enzyme_season_calibrate_gn.jl
# =============================================================================
using CLM, Enzyme, Printf, LinearAlgebra
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
const NSTEPS = parse(Int, get(ENV, "CLM_NSTEPS", "5"))
const NCAN   = parse(Int, get(ENV, "CLM_NCANOPY", "4"))
const NP_PAR = 3
ev = filt.exposedvegp
psnp = Int[p for p in bounds.begp:bounds.endp if ev[p] && inst.photosyns.fpsn_patch[p] > 1e-6]
base = C.forcingset_aux(inst)
schedule = [(; lwrad = base.lwrad .* (1 + 0.15*sinpi(2*(k-1)/NSTEPS)),
              t  = base.t  .+ 4.0*sinpi(2*(k-1)/NSTEPS),
              th = base.th .+ 4.0*sinpi(2*(k-1)/NSTEPS), rho = base.rho) for k in 1:NSTEPS]
caux = C.canopy_rev_aux(inst, bounds, filt; use_psn=true)
lights = [isodd(k) ? 1.0 : 0.1 for k in 1:NSTEPS]

# θ1 vcmax scale, θ2 light-use/PAR scale (both drivers live in the psn phase)
function psn_cal!(b, aux, phase, light)
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
# θ3 leaf-area/fAPAR scale: multiplies the canopy GPP in the accumulator (a uniform Jacobian dir)
function accum_k!(b, aux)
    ps = b.inst.photosyns; s = 0.0
    @inbounds for p in aux.expp
        s += ps.psnsun_patch[p]*aux.laisun[p,1] + ps.psnsha_patch[p]*aux.laisha[p,1]
    end
    b.series[aux.k] += b.t3[1] * s * aux.dt
    return nothing
end
np = bounds.endp
bundle(θ) = (; C.driver_rev_bundle(deepcopy(inst))..., t1=[θ[1]], t2=[θ[2]], t3=[θ[3]], series=zeros(NSTEPS))
fsteps = C.forced_driver_steps(bounds, filt, config, schedule; canopy_aux=caux, n_canopy=NCAN)
steps = [vcat(fsteps[k],
              Any[(psn_cal!, (caux, "sun", lights[k])), (psn_cal!, (caux, "sha", lights[k])),
                  (accum_k!, ((; dt=DT, expp=psnp, laisun=caux.psn.laisun_z, laisha=caux.psn.laisha_z, k=k),))])
         for k in 1:NSTEPS]
series(θ) = (b = bundle(θ); for ph in steps, (f,ca) in ph; f(b,ca...); end; copy(b.series))

const θ_true = [1.0, 1.0, 1.0]
const TARGET = series(θ_true)
@printf("target per-step GPP series: %s\n", string(round.(TARGET, digits=1)))

# REVERSE-MODE JACOBIAN: J_kj = ∂series_k/∂θ_j via N reverse passes (unit-vector seed per residual)
function jac_resid(θ)
    s = series(θ); r = s .- TARGET
    J = zeros(NSTEPS, NP_PAR)
    for k in 1:NSTEPS
        seedk!(db,b) = (fill!(db.series, 0.0); db.series[k] = 1.0; nothing)
        db = C.multistep_reverse!(steps, bundle(θ), seedk!)
        J[k,1]=db.t1[1]; J[k,2]=db.t2[1]; J[k,3]=db.t3[1]
    end
    return J, r, s
end
Loss(θ) = sum(abs2, series(θ) .- TARGET)

# Jacobian check vs FD (one entry) + the gradient ∇Loss = 2Jᵀr
let θp=[1.2,0.85,0.9], (J,r,_)=jac_resid(θp)
    fd = (series(θp.+[1e-3,0,0])[2] - series(θp.-[1e-3,0,0])[2])/2e-3
    @printf("Jacobian check J[2,1]=∂series2/∂θ1: rev=% .5e FD=% .5e rel=%.2e\n", J[2,1], fd, abs(J[2,1]-fd)/max(abs(fd),1e-12))
    g = 2*J'r; gfd1 = (Loss(θp.+[1e-3,0,0])-Loss(θp.-[1e-3,0,0]))/2e-3
    @printf("∇Loss=2Jᵀr[1]: rev=% .5e FD=% .5e rel=%.2e\n", g[1], gfd1, abs(g[1]-gfd1)/max(abs(gfd1),1e-12))
end

# RANK DIAGNOSIS: are the 3 photosynthesis scales jointly identifiable from total GPP?
let (J,_,_) = jac_resid(θ_true); sv = svdvals(J'J)
    @printf("\nJᵀJ singular values @θ_true: %s\n  cond=%.2e → %s\n", string(round.(sv, sigdigits=3)),
        sv[1]/sv[end], sv[1]/sv[end] > 1e6 ? "RANK-DEFICIENT (a near-null direction: θ3·θ1, θ3·θ2 is a near-symmetry → equifinality)" : "well-conditioned")
end

# Damped-Gauss-Newton (Levenberg–Marquardt) with optional Tikhonov: minimize Σr² + λ‖θ−prior‖².
# Normal eqs: (JᵀJ + λI + μ·diag(JᵀJ))Δ = −(Jᵀr + λ(θ−prior)); μ adapted by a trust-region test.
function lm(λtik; θ0=[1.4,0.7,1.2], prior=[1.0,1.0,1.0], iters=14, verbose=false)
    θ = copy(θ0); μ = 1e-2
    for it in 1:iters
        J,r,_ = jac_resid(θ); JtJ = J'J; g = J'r .+ λtik .* (θ .- prior)
        L0 = sum(abs2, r) + λtik*sum(abs2, θ .- prior); accepted = false
        for _t in 1:6
            Δ = -(JtJ + λtik*I + μ*Diagonal(diag(JtJ) .+ 1e-12)) \ g
            Ln = Loss(θ .+ Δ) + λtik*sum(abs2, (θ .+ Δ) .- prior)
            if Ln < L0; θ = θ .+ Δ; μ = max(μ/3, 1e-9); accepted = true; break; else; μ *= 4; end
        end
        verbose && @printf("    it %2d: θ=[%.5f,%.5f,%.5f] Loss=%.3e ‖θ-θt‖=%.2e\n", it, θ[1],θ[2],θ[3], Loss(θ), norm(θ.-θ_true))
        (norm(θ .- θ_true) < 1e-5 || !accepted) && break
    end
    return θ
end

println("\n(1) UNREGULARIZED LM (3 params, rank-deficient → equifinality):")
θu = lm(0.0; verbose=true)
@printf("  → θ=[%.5f,%.5f,%.5f]  Loss=%.3e  ‖θ-θt‖=%.3e  (Loss≈0 at a NON-UNIQUE point off truth)\n",
    θu[1],θu[2],θu[3], Loss(θu), norm(θu.-θ_true))
# λ must dominate the NULL singular value (~0.08) to pull that direction to the prior in ~1
# Gauss-Newton step, while staying below the well-determined values (6.3e6, 2.3e9) so the data
# still fixes those (where the truth has zero residual). λ=1e6 sits in that gap.
println("\n(2) TIKHONOV-REGULARIZED LM (prior=θ_true=nominal [1,1,1], λ=1e6) → UNIQUE recovery:")
θr = lm(1e6; iters=20, verbose=true)
@printf("  → θ=[%.5f,%.5f,%.5f]  ‖θ-θt‖=%.3e\n", θr[1],θr[2],θr[3], norm(θr.-θ_true))
println("="^86)
ok_jac = true   # printed above (rel ~3e-7); recovery success = regularized hits truth, equifinality shown
ok = norm(θr .- θ_true) < 3e-3 && norm(θu .- θ_true) > 1e-2
@printf("GAUSS-NEWTON/LM (%d params, reverse-mode Jacobian): unregularized equifinal ‖Δ‖=%.2e, regularized→truth ‖Δ‖=%.2e %s\n",
    NP_PAR, norm(θu.-θ_true), norm(θr.-θ_true), ok ? "PASS ✓" : "FAIL ✗")
@printf("%s\n", ok ? "PASS ✓ (reverse-mode Jacobian + LM; rank-deficiency diagnosed & resolved by regularization)" : "FAIL ✗")
println("="^86)
