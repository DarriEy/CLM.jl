# =============================================================================
# QUANTITATIVE SEASON-GPP GRADIENT — d(Σ_k GPP_k·dt)/d(absorbed PAR) over a forced multi-step
# trajectory on an ACTIVELY-PHOTOSYNTHESIZING inst. Resolves the THREE obstacles that made
# enzyme_multistep_gpp.jl Part B a trivial zero:
#  (1) the synthetic Bow cold-start did not assimilate → here the inst is warmed with the
#      growing-season recipe (prescribed July LAI + soil moisture + June-solstice sun) so the
#      canopy photosynthesizes (~12 µmol/m²/s);
#  (2) photosynthesis! writes the PER-LEAF rates psnsun_patch/psnsha_patch; the canopy-integrated
#      fpsn_patch is computed DOWNSTREAM (not in the reduced reverse) so fpsn stays FROZEN at the
#      warmup value (FD-confirmed: parsun=0 leaves fpsn unchanged). → accumulate GPP from the LIVE
#      per-leaf rates: GPP_p = psnsun·laisun + psnsha·laisha;
#  (3) the GPP drivers in the canopy aux are Const snapshots → d(GPP)/d(bundle state) is 0. → a
#      custom psn phase reads ABSORBED PAR (parsun/parsha) LIVE from inst.solarabs (Duplicated),
#      the confirmed-live driver (parsun=0 ⇒ psnsun=0). The gradient is the integrated sensitivity
#      of seasonal carbon uptake to absorbed light (a radiation-use / fAPAR sensitivity).
#
# Validated: d(Σ GPP·dt)/d(parsun) — multistep_reverse! == binomial == FD (rel ~1e-10).
#   CLM_NSTEPS=3 CLM_NCANOPY=6 julia +1.10 --project=/tmp/clm_jl10_<id> scripts/enzyme_season_gpp.jl
# =============================================================================
using CLM, Enzyme, Printf
const C = CLM
const FSURDAT  = get(ENV, "CLM_FSURDAT",  "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/surfdata_clm.nc")
const PARAMFILE = get(ENV, "CLM_PARAMFILE", "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/clm5_params.nc")
const DT = 1800.0

# Warmed growing-season inst (mirrors build_cn_inst's recipe; use_cn=false SP/prescribed-LAI path).
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
    for n in 1:4                                   # warm: radiation/t_veg/btran spin up; fpsn>0 by step 2
        C.clm_drv!(config, inst, filt, fia, bounds, true, nextsw, declin, declin,
            C.ORB_OBLIQR_DEFAULT, false, false, "", false; nstep=n, is_first_step=(n==1),
            is_beg_curr_day=(n==1), dtime=DT, mon=7, day=15, photosyns=inst.photosyns)
    end
    return inst, bounds, filt, config
end

inst, bounds, filt, config = build_growing_season()
const NSTEPS = parse(Int, get(ENV, "CLM_NSTEPS", "3"))
const NCAN   = parse(Int, get(ENV, "CLM_NCANOPY", "6"))
ev = filt.exposedvegp; expp = Int[p for p in bounds.begp:bounds.endp if ev[p]]
# photosynthesizing exposed patches only (some exposed PFTs assimilate 0)
psnp = Int[p for p in expp if inst.photosyns.fpsn_patch[p] > 1e-6]
fp = first(psnp); c_fp = inst.patch.column[fp]
@printf("warmed growing-season inst: %d photosynthesizing patches; fp=%d (fpsn=%.4f) col=%d\n",
    length(psnp), fp, inst.photosyns.fpsn_patch[fp], c_fp)

# Diurnal forcing schedule about the warmed base (a real day's swing).
base = C.forcingset_aux(inst)
schedule = [(; lwrad = base.lwrad .* (1 + 0.15*sinpi(2*(k-1)/NSTEPS)),
              t  = base.t  .+ 4.0*sinpi(2*(k-1)/NSTEPS),
              th = base.th .+ 4.0*sinpi(2*(k-1)/NSTEPS),
              rho = base.rho) for k in 1:NSTEPS]

caux = C.canopy_rev_aux(inst, bounds, filt; use_psn=true)

# GPP accumulator: photosynthesis! writes the PER-LEAF rates psnsun_patch/psnsha_patch (LIVE in
# the Duplicated inst); the canopy-integrated fpsn_patch is computed DOWNSTREAM (not in the reduced
# reverse) so it stays frozen at the warmup value. Accumulate the canopy GPP from the live per-leaf
# rates × sunlit/shaded leaf area: GPP_p = psnsun·laisun + psnsha·laisha.
auxG = (; dt = DT, expp = psnp, laisun = caux.psn.laisun_z, laisha = caux.psn.laisha_z)
function accum_gpp!(b, aux)
    ps = b.inst.photosyns
    @inbounds for p in aux.expp
        gpp = ps.psnsun_patch[p]*aux.laisun[p,1] + ps.psnsha_patch[p]*aux.laisha[p,1]
        b.acc[p] += gpp * aux.dt
    end
    return nothing
end

# Custom psn phase: recompute the per-leaf rates reading ABSORBED PAR (parsun/parsha) LIVE from
# inst.solarabs (Duplicated) instead of the Const caux.psn copy → makes absorbed PAR a
# differentiable GPP driver (psnsun responds to parsun: parsun=0 ⇒ psnsun=0, FD-confirmed). Runs
# AFTER the canopy block (which populated the psn scratch svpts/eah/co2/o2/rb/dayl_factor + t_veg);
# mirrors the production cf_rev_psn! call verbatim except parz is the live inst.solarabs array.
function psn_live!(b, aux, phase)
    sc = b.scratch; tp = b.inst.temperature; ps = b.inst.photosyns; pn = aux.psn
    vcm  = phase == "sun" ? pn.vcmaxcint_sun : pn.vcmaxcint_sha
    parz = phase == "sun" ? b.inst.solarabs.parsun_z_patch : b.inst.solarabs.parsha_z_patch   # LIVE
    laiz = phase == "sun" ? pn.laisun_z : pn.laisha_z
    C.photosynthesis!(ps, sc.svpts, sc.eah, sc.o2_arr, sc.co2_arr, sc.rb,
        b.inst.energyflux.btran_patch, sc.dayl_factor, aux.forc.leafn, aux.forc_pbot_patch,
        tp.t_veg_patch, pn.t10, tp.thm_patch, pn.nrad, pn.tlai_z, b.inst.canopystate.tlai_patch,
        parz, laiz, vcm, pn.o3coefv, pn.o3coefg, pn.c3psn, pn.leafcn, pn.flnr, pn.fnitr,
        pn.slatop, pn.mbbopt, pn.medlynintercept, pn.medlynslope, aux.ivt, aux.patch.column,
        aux.mask, 1:length(sc.air), phase)
    return nothing
end

np = bounds.endp
bundle() = (; C.driver_rev_bundle(deepcopy(inst))..., acc = zeros(np))
# step = forcingset + canopy(energy+psn, populates scratch & t_veg) + live-PAR psn (sun+sha →
# differentiable psnsun/psnsha) + hydrology + GPP accumulate.
steps = [vcat(phs, Any[(psn_live!, (caux, "sun")), (psn_live!, (caux, "sha")),
                       (accum_gpp!, (auxG,))])
         for phs in C.forced_driver_steps(bounds, filt, config, schedule; canopy_aux=caux, n_canopy=NCAN)]
@printf("forced GPP trajectory: %d steps × %d-phase chain (forcingset+canopy+live-PAR-psn+hydro+accum)\n",
    NSTEPS, length(steps[1]))

const PARLAY = 1                                   # absorbed-PAR canopy layer to perturb (nrad=1)
L(b) = sum(@view b.acc[psnp])                      # season GPP = Σ_k Σ_p∈psn (psnsun·laisun+psnsha·laisha)·dt
pert!(b, δ) = (b.inst.solarabs.parsun_z_patch[fp, PARLAY] += δ; nothing)
function run(δ)
    b = bundle(); pert!(b, δ)
    for phases in steps, (f, ca) in phases; f(b, ca...); end
    return b
end
let bf = run(0.0)
    @printf("primal season GPP = %.6e (finite=%s); FD probe d/d(parsun) = %.4e\n", L(bf), string(isfinite(L(bf))),
        (L(run(5.0))-L(run(-5.0)))/10.0)
end

seed!(db, b) = (db.acc .= 0.0; for p in psnp; db.acc[p] = 1.0; end; nothing)
gradfield(db) = db.inst.solarabs.parsun_z_patch[fp, PARLAY]
gG  = gradfield(C.multistep_reverse!(steps, bundle(), seed!))
pk = Ref(0)
gGb = gradfield(C.multistep_reverse_binomial!(steps, bundle(), seed!; peak_checkpoints=pk))
hs = (5.0, 2.5, 1.25); cfd(h) = (L(run(h)) - L(run(-h)))/(2h)
fd = [cfd(h) for h in hs]; for (h,f) in zip(hs,fd); @printf("  FD(h=%.2e) = % .8e\n", h, f); end
g_rich = (4*fd[end]-fd[end-1])/3; lo,hi = minimum(fd),maximum(fd)
rel_bin = abs(gG-gGb)/max(abs(gGb),1e-30); rel_fd = abs(gG-g_rich)/max(abs(g_rich),1e-12)
brk = (lo-abs(lo)*1e-9) <= gG <= (hi+abs(hi)*1e-9)
println("\n", "="^78)
@printf("SEASON-GPP REVERSE  d(Σ GPP·dt)/d(absorbed PAR parsun_z[patch %d])  over %d forced steps (canopy+psn, NCAN=%d)\n", fp, NSTEPS, NCAN)
@printf("  multistep_reverse!          = % .8e\n", gG)
@printf("  multistep_reverse_binomial! = % .8e   rel=%.2e (peak %d snapshots)\n", gGb, rel_bin, pk[])
@printf("  FD Richardson               = % .8e   rel=%.2e  bracketed=%s\n", g_rich, rel_fd, string(brk))
pass = abs(gG) > 1e-6 && rel_bin < 1e-8 && (rel_fd < 1e-2 || brk)
@printf("%s\n", pass ? "PASS ✓ (NONZERO quantitative season-GPP gradient w.r.t. absorbed PAR, integrated over the forced horizon)" : "FAIL/zero ✗")
println("="^78)
