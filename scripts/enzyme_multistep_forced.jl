# =============================================================================
# MULTI-TIMESTEP reverse with PER-STEP (TIME-VARYING) FORCING — the real-trajectory version of
# enzyme_multistep_reverse.jl (which held forcing fixed). Each timestep injects its own
# downscaled atmospheric forcing via CLM.forcingset_rev_phase! (mutates the live b.inst.atm2lnd
# arrays that soil_temperature! reads), so the soil thermal state evolves under a DIURNAL cycle
# of longwave/air-temperature forcing. Validates that d(final t_soisno)/d(initial t_soisno)
# across the forced horizon is correct — i.e. the per-step forcing-advance integrates cleanly
# with both multistep engines.
#
# Schedule: a sinusoidal diurnal swing on (lwrad, t, th, rho) about the warmed base forcing.
# THREE-WAY: FD over the forced trajectory == multistep_reverse! == multistep_reverse_binomial!.
#   CLM_NSTEPS=6 julia +1.10 --project=/tmp/clm_jl10_<id> scripts/enzyme_multistep_forced.jl
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
const NSTEPS = parse(Int, get(ENV, "CLM_NSTEPS", "6"))
const JLAY   = parse(Int, get(ENV, "CLM_JLAY", "16"))
hc = filt.hydrologyc; c0 = [c for c in bounds.begc:bounds.endc if hc[c]][1]

# Diurnal forcing schedule: scale the warmed base downscaled forcing by a per-step sinusoid.
# step k (1..N) at phase 2π(k-1)/N: lwrad ×(1+0.25·s), t/th +6·s K, rho tracks t.
base = C.forcingset_aux(inst)
function schedule_step(k)
    s = sinpi(2*(k-1)/NSTEPS)
    lwrad = base.lwrad .* (1 + 0.25*s)
    t  = base.t  .+ 6.0*s
    th = base.th .+ 6.0*s
    rho = base.rho .* (base.t ./ t)          # ρ ∝ 1/T at fixed p
    return (; lwrad, t, th, rho)
end
schedule = [schedule_step(k) for k in 1:NSTEPS]
steps = C.forced_driver_steps(bounds, filt, config, schedule; canopy_aux=nothing)
@printf("forced trajectory: %d steps × %d-phase chain (forcingset + hydrology), diurnal lwrad/t swing\n",
    NSTEPS, length(steps[1]))

# Forward the forced trajectory from a bundle, perturbing initial t_soisno[c0,JLAY].
function run_forced(δ)
    bs = C.driver_rev_bundle(deepcopy(inst))
    bs.inst.temperature.t_soisno_col[c0, JLAY] += δ
    for phases in steps, (f, ca) in phases; f(bs, ca...); end
    return bs
end
Lfinal(b) = sum(abs2, b.inst.temperature.t_soisno_col)
let bf = run_forced(0.0), b0 = run_forced(0.0)   # confirm forcing actually drives the state
    # compare to a FIXED-forcing run (no schedule) to show the diurnal swing matters
    bfix = C.driver_rev_bundle(deepcopy(inst))
    fixphases = C.driver_rev_phases(bounds, filt, config; canopy_aux=nothing)
    for _ in 1:NSTEPS, (f, ca) in fixphases; f(bfix, ca...); end
    @printf("forced final t_soisno[%d,%d]=%.5f vs fixed-forcing %.5f (Δ=%+.4e → forcing drives it)\n",
        c0, JLAY, bf.inst.temperature.t_soisno_col[c0,JLAY],
        bfix.inst.temperature.t_soisno_col[c0,JLAY],
        bf.inst.temperature.t_soisno_col[c0,JLAY]-bfix.inst.temperature.t_soisno_col[c0,JLAY])
end

seed!(db, b) = (db.inst.temperature.t_soisno_col .= 2 .* b.inst.temperature.t_soisno_col)
db_ms  = C.multistep_reverse!(steps, C.driver_rev_bundle(deepcopy(inst)), seed!)
g_ms   = db_ms.inst.temperature.t_soisno_col[c0, JLAY]
peak = Ref(0)
db_bin = C.multistep_reverse_binomial!(steps, C.driver_rev_bundle(deepcopy(inst)), seed!; peak_checkpoints=peak)
g_bin  = db_bin.inst.temperature.t_soisno_col[c0, JLAY]

hs = (1e-2, 5e-3, 2.5e-3); cfd(h) = (Lfinal(run_forced(h)) - Lfinal(run_forced(-h))) / (2h)
fds = [cfd(h) for h in hs]; for (h,f) in zip(hs,fds); @printf("  FD(h=%.2e) = % .8e\n", h, f); end
g_rich = (4*fds[end] - fds[end-1]) / 3; lo, hi = minimum(fds), maximum(fds)

rel_bin = abs(g_ms - g_bin) / max(abs(g_bin), 1e-30)
rel_fd  = abs(g_ms - g_rich) / max(abs(g_rich), 1e-12)
bracketed = (lo - abs(lo)*1e-9) <= g_ms <= (hi + abs(hi)*1e-9)
println("\n", "="^74)
@printf("FORCED MULTI-STEP REVERSE  d(final t_soisno)/d(initial t_soisno[%d,%d])  over %d forced steps\n", c0, JLAY, NSTEPS)
@printf("  multistep_reverse!           = % .8e\n", g_ms)
@printf("  multistep_reverse_binomial!  = % .8e   rel=%.2e (peak %d snapshots)\n", g_bin, rel_bin, peak[])
@printf("  FD Richardson                = % .8e   rel=%.2e  bracketed=%s\n", g_rich, rel_fd, string(bracketed))
pass = rel_bin < 1e-8 && (rel_fd < 1e-4 || bracketed)
@printf("%s\n", pass ? "PASS ✓ (per-step forcing advance integrates with both engines; gradient matches FD)" : "FAIL ✗")
println("="^74)
