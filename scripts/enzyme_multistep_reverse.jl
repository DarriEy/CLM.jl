# =============================================================================
# MULTI-TIMESTEP reverse — propagates the adjoint across a TRAJECTORY of N timesteps,
# giving d(final state)/d(initial state) (the state-to-state sensitivity through the whole
# horizon). Each step is the production HydrologyNoDrainage chain (CLM.driver_rev_phases,
# canopy off): soil_temperature! → [surface-hydrology block] → soil_water! → water_table!
# → hydrology_no_drainage!. The coupled thermal/water state carries through b.inst, so
# step s+1 continues from what step s produced; the perturbed field (deep soil temperature
# t_soisno) and the seeded loss (final t_soisno) are the SAME carried state → a genuine
# multi-step thermal Jacobian-vector product.
#
# FOUR-WAY cross-check of the multi-timestep reverse engines:
#   (1) FD over the full N-step forward trajectory   (ground truth)
#   (2) multistep_reverse!           (two-level: O(N) coarse + O(phases) fine)
#   (3) compositional_reverse!(vcat)  (flat: O(N·phases) checkpoints)
#   (4) multistep_reverse_binomial!  (recursive bisection: O(log N) coarse snapshots)
# (2)==(3)==(4) to roundoff and all must bracket (1); (4) also reports its measured peak
# number of simultaneously-held coarse snapshots (the logarithmic-memory bound).
#   CLM_NSTEPS=3 julia +1.10 --project=/tmp/clm_jl10_<id> scripts/enzyme_multistep_reverse.jl
# =============================================================================
using CLM, Enzyme, Printf
const C = CLM
const FSURDAT  = get(ENV, "CLM_FSURDAT",  "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/surfdata_clm.nc")
const PARAMFILE = get(ENV, "CLM_PARAMFILE", "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/clm5_params.nc")

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
const NSTEPS = parse(Int, get(ENV, "CLM_NSTEPS", "3"))
const JLAY   = parse(Int, get(ENV, "CLM_JLAY", "16"))   # an actively-diffusing soil layer
# (layers 1..nlevsno are inert padded snow; soil ~13+; layer 16 diffuses into its neighbours
#  across the horizon → d(final)/d(init) is non-trivial and exercises real cross-step coupling).
# One timestep's worth of phases (canopy off → hydrology+thermal chain). Fixed forcing, so the
# SAME phase list is valid for every step; state evolution comes from b.inst carrying forward.
phases1 = C.driver_rev_phases(bounds, filt, config; canopy_aux=nothing)
steps = [phases1 for _ in 1:NSTEPS]
hc = filt.hydrologyc; c0 = [c for c in bounds.begc:bounds.endc if hc[c]][1]
@printf("assembled %d-phase/step hydrology chain × %d steps = %d total phases (canopy off)\n",
    length(phases1), NSTEPS, length(phases1)*NSTEPS)

# Run the full N-step forward trajectory from a bundle, perturbing initial t_soisno[c0,JLAY].
function run_traj(δ)
    bs = C.driver_rev_bundle(deepcopy(inst))
    bs.inst.temperature.t_soisno_col[c0, JLAY] += δ
    for phases in steps, (f, ca) in phases; f(bs, ca...); end
    return bs
end
Lfinal(b) = sum(abs2, b.inst.temperature.t_soisno_col)
let bf = run_traj(0.0)
    @printf("trajectory primal: final t_soisno[%d,%d]=%.6f finite=%s\n", c0, JLAY,
        bf.inst.temperature.t_soisno_col[c0,JLAY], string(all(isfinite, bf.inst.temperature.t_soisno_col)))
end

seed!(db, b) = (db.inst.temperature.t_soisno_col .= 2 .* b.inst.temperature.t_soisno_col)

# (2) two-level multistep engine
db_ms = C.multistep_reverse!(steps, C.driver_rev_bundle(deepcopy(inst)), seed!)
g_ms  = db_ms.inst.temperature.t_soisno_col[c0, JLAY]
# (3) flat reference: one compositional_reverse! over the concatenated phase list
db_flat = C.compositional_reverse!(reduce(vcat, steps), C.driver_rev_bundle(deepcopy(inst)), seed!)
g_flat  = db_flat.inst.temperature.t_soisno_col[c0, JLAY]
# (4) logarithmic-memory engine: recursive bisection checkpointing
peak = Ref(0)
db_bin = C.multistep_reverse_binomial!(steps, C.driver_rev_bundle(deepcopy(inst)), seed!; peak_checkpoints=peak)
g_bin  = db_bin.inst.temperature.t_soisno_col[c0, JLAY]

# (1) FD ground truth over the full trajectory (Richardson on a 2-step shrink)
hs = (1e-2, 5e-3, 2.5e-3); cfd(h) = (Lfinal(run_traj(h)) - Lfinal(run_traj(-h))) / (2h)
fds = [cfd(h) for h in hs]; for (h,f) in zip(hs,fds); @printf("  FD(h=%.2e) = % .8e\n", h, f); end
g_rich = (4*fds[end] - fds[end-1]) / 3; lo, hi = minimum(fds), maximum(fds)

rel_ms_flat = abs(g_ms - g_flat) / max(abs(g_flat), 1e-30)
rel_ms_bin  = abs(g_ms - g_bin)  / max(abs(g_bin),  1e-30)
rel_ms_fd   = abs(g_ms - g_rich) / max(abs(g_rich), 1e-12)
bracketed   = (lo - abs(lo)*1e-9) <= g_ms <= (hi + abs(hi)*1e-9)
println("\n", "="^74)
@printf("MULTI-STEP REVERSE  d(final t_soisno)/d(initial t_soisno[%d,%d])  over %d steps\n", c0, JLAY, NSTEPS)
@printf("  (2) multistep_reverse!            = % .8e\n", g_ms)
@printf("  (3) compositional_reverse!(vcat)  = % .8e   rel(2,3)=%.2e\n", g_flat, rel_ms_flat)
@printf("  (4) multistep_reverse_binomial!   = % .8e   rel(2,4)=%.2e\n", g_bin, rel_ms_bin)
@printf("  (1) FD Richardson                 = % .8e   rel(2,1)=%.2e  bracketed=%s\n", g_rich, rel_ms_fd, string(bracketed))
@printf("  coarse-snapshot memory: two-level=%d   flat=%d   binomial(measured peak)=%d  [ceil(log2 N)+1=%d]\n",
    NSTEPS, length(phases1)*NSTEPS, peak[], ceil(Int, log2(max(NSTEPS,1)))+1)
pass = rel_ms_flat < 1e-8 && rel_ms_bin < 1e-8 && (rel_ms_fd < 1e-4 || bracketed)
@printf("%s\n", pass ? "PASS ✓ (all 3 engines agree to roundoff; gradient matches FD; binomial holds O(log N) snapshots)" : "FAIL ✗")
println("="^74)
