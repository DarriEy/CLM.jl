# =============================================================================
# MULTI-PARAMETER SEASON-GPP CALIBRATION — jointly recover TWO photosynthesis parameters
# (vcmax scale θ1, light-use/PAR scale θ2) from a per-step GPP series through the multi-timestep
# reverse. Extends enzyme_season_calibrate.jl (single param) to a gradient VECTOR: one reverse
# pass yields (∂Loss/∂θ1, ∂Loss/∂θ2) since both params live in the shared bundle and every step's
# psn phase reads them.
#
# IDENTIFIABILITY: two params that both scale total GPP are degenerate against a single scalar
# target. Here the loss is over the PER-STEP GPP SERIES under a strong light contrast (alternating
# bright/dim steps): bright steps are Rubisco-limited → sensitive to θ1 (vcmax); dim steps are
# light-limited → sensitive to θ2 (PAR). The series separates them → a well-conditioned 2×2.
#
# Validates: per-COMPONENT grad check (∂Loss/∂θ1, ∂Loss/∂θ2 each vs FD), then FULL 2×2-Newton
# recovery of both (the params are cross-coupled — significant off-diagonal Hessian — so a
# diagonal step OSCILLATES; the full 2×2 Newton converges in ~4 iters). Light contrast is applied
# to the assimilation PAR in the custom psn phase.
#   CLM_NSTEPS=4 CLM_NCANOPY=4 julia +1.10 --project=/tmp/clm_jl10_<id> scripts/enzyme_season_calibrate_multi.jl
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
const NSTEPS = parse(Int, get(ENV, "CLM_NSTEPS", "4"))
const NCAN   = parse(Int, get(ENV, "CLM_NCANOPY", "4"))
ev = filt.exposedvegp
psnp = Int[p for p in bounds.begp:bounds.endp if ev[p] && inst.photosyns.fpsn_patch[p] > 1e-6]
@printf("warmed inst: %d photosynthesizing patch(es)\n", length(psnp))

base = C.forcingset_aux(inst)
schedule = [(; lwrad = base.lwrad .* (1 + 0.15*sinpi(2*(k-1)/NSTEPS)),
              t  = base.t  .+ 4.0*sinpi(2*(k-1)/NSTEPS),
              th = base.th .+ 4.0*sinpi(2*(k-1)/NSTEPS), rho = base.rho) for k in 1:NSTEPS]
caux = C.canopy_rev_aux(inst, bounds, filt; use_psn=true)
# strong diurnal light contrast: alternate bright (Rubisco-limited→θ1) / dim (light-limited→θ2)
lights = [isodd(k) ? 1.0 : 0.1 for k in 1:NSTEPS]
@printf("light schedule (per-step PAR scale): %s\n", string(lights))

# psn phase: vcmax scaled by θ1 (b.t1), absorbed PAR scaled by θ2·light_k (b.t2). Both LIVE drivers.
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
# accumulate this step's GPP into series[k] (per-step observation)
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

const θ1t, θ2t = 1.0, 1.0
const TARGET = series(θ1t, θ2t)
@printf("target per-step GPP series @ (θ1,θ2)=(1,1): %s\n", string(round.(TARGET, digits=1)))
Loss(θ1,θ2) = sum(abs2, series(θ1,θ2) .- TARGET)
function gradL(θ1,θ2)
    seed!(db,b) = (db.series .= 2 .* (b.series .- TARGET); nothing)
    db = C.multistep_reverse!(steps, bundle(θ1,θ2), seed!)
    return db.t1[1], db.t2[1]
end

# ---- per-component grad check vs FD at a perturbed (θ1,θ2) ----
a0, b0 = 1.3, 0.7
g1, g2 = gradL(a0, b0)
fd1 = (Loss(a0+1e-3,b0)-Loss(a0-1e-3,b0))/2e-3
fd2 = (Loss(a0,b0+1e-3)-Loss(a0,b0-1e-3))/2e-3
@printf("\ngrad check @ (%.2f,%.2f):\n  ∂L/∂θ1 rev=% .5e FD=% .5e rel=%.2e\n  ∂L/∂θ2 rev=% .5e FD=% .5e rel=%.2e\n",
    a0,b0, g1,fd1,abs(g1-fd1)/max(abs(fd1),1e-12), g2,fd2,abs(g2-fd2)/max(abs(fd2),1e-12))

# ---- joint FULL 2×2-Newton recovery (the params are cross-coupled → the off-diagonal Hessian
# matters; a diagonal step oscillates). Hessian columns = FD of the gradient VECTOR; θ -= H⁻¹g. ----
@printf("\nfull 2×2-Newton recovery (θ_true=(%.2f,%.2f), start (%.2f,%.2f)):\n", θ1t,θ2t,a0,b0)
θ1, θ2 = a0, b0
for it in 1:8
    global θ1, θ2
    g1,g2 = gradL(θ1,θ2); h = 1e-3
    g1p,g2p = gradL(θ1+h, θ2)          # ∂(g1,g2)/∂θ1
    g1q,g2q = gradL(θ1, θ2+h)          # ∂(g1,g2)/∂θ2
    H11=(g1p-g1)/h; H21=(g2p-g2)/h; H12=(g1q-g1)/h; H22=(g2q-g2)/h
    det = H11*H22 - H12*H21
    d1 = -( H22*g1 - H12*g2)/det; d2 = -(-H21*g1 + H11*g2)/det
    θ1 += d1; θ2 += d2
    @printf("  it %2d: θ=(%.6f,%.6f)  Loss=%.4e  ‖θ-θ_true‖=%.2e\n",
        it, θ1, θ2, Loss(θ1,θ2), hypot(θ1-θ1t, θ2-θ2t))
    hypot(θ1-θ1t, θ2-θ2t) < 1e-4 && break
end
println("="^80)
ok_g = abs(g1-fd1)/max(abs(fd1),1e-12) < 1e-3 && abs(g2-fd2)/max(abs(fd2),1e-12) < 1e-3
ok_r = hypot(θ1-θ1t, θ2-θ2t) < 2e-3
@printf("MULTI-PARAM SEASON-GPP CALIBRATION: grad-check(2-comp) %s, recovered θ=(%.5f,%.5f) true (1,1) %s\n",
    ok_g ? "PASS" : "FAIL", θ1, θ2, ok_r ? "PASS ✓" : "FAIL ✗")
@printf("%s\n", (ok_g && ok_r) ? "PASS ✓ (two photosynthesis params jointly recovered through the season reverse chain)" : "FAIL ✗")
println("="^80)
