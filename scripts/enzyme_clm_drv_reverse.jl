# =============================================================================
# TOP-LEVEL clm_drv_reverse! validation — the productionized, CONVERGENCE-AWARE
# driver-reverse entry point (src/driver/driver_reverse.jl). Given a warmed-up inst
# + a seed, clm_drv_reverse! auto-detects the canopy Newton iteration count per step
# (driver_canopy_converged_n) and runs the compositional reverse over the whole
# clm_drv! step (canopy → soil_temp → surface-hydrology → soil_water → water_table →
# hydrology_no_drainage). Here we:
#   1. confirm canopy_rev_converged_n reproduces the real canopy_fluxes! num_iter,
#   2. FD-validate dL/d(t_grnd) through the whole step with the AUTO-detected N
#      (i.e. n_canopy=nothing — no hard-coded iteration count).
#
#   julia +1.10 --project=<env-with-CLM+Enzyme> scripts/enzyme_clm_drv_reverse.jl
# =============================================================================
using CLM, Enzyme, Printf
const C = CLM
const FSURDAT   = get(ENV, "CLM_FSURDAT",   "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/surfdata_clm.nc")
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
use_psn = get(ENV, "CLM_USE_PSN", "0") == "1"

# ---- (1) convergence-awareness: auto-detected N vs the real canopy_fluxes! num_iter ----
caux = C.canopy_rev_aux(inst, bounds, filt; use_psn=use_psn)
b0   = C.driver_rev_bundle(deepcopy(inst))
Nauto = C.driver_canopy_converged_n(b0, caux)
# the real per-patch Newton count is on the warmed-up inst's frictionvel (NaN = inactive)
ev = filt.exposedvegp
nit_active = [Int(inst.frictionvel.num_iter_patch[p]) for p in bounds.begp:bounds.endp
              if ev[p] && isfinite(inst.frictionvel.num_iter_patch[p]) &&
                 inst.frictionvel.num_iter_patch[p] > 0]
@printf("auto-detected N = %d   |   real canopy num_iter (active veg) = %s\n",
        Nauto, isempty(nit_active) ? "n/a" : string(nit_active))

# ---- (2) FD-validate dL/d(t_grnd) through the whole step with AUTO N ----
hc = filt.hydrologyc; c0 = [c for c in bounds.begc:bounds.endc if hc[c]][1]
seed!(db, b) = (db.inst.water.waterstatebulk_inst.ws.h2osoi_vol_col .=
                2 .* b.inst.water.waterstatebulk_inst.ws.h2osoi_vol_col)

# forward L(δ): same chain clm_drv_reverse! builds, but only forward-swept.
function L_pert(δ)
    bs = C.driver_rev_bundle(deepcopy(inst)); bs.inst.temperature.t_grnd_col[c0] += δ
    cax = C.canopy_rev_aux(bs.inst, bounds, filt; use_psn=use_psn)
    Nc  = C.driver_canopy_converged_n(bs, cax)
    ph  = C.driver_rev_phases(bounds, filt, config; canopy_aux=cax, n_canopy=max(Nc,1))
    for (f, ca) in ph; f(bs, ca...); end
    return sum(abs2, bs.inst.water.waterstatebulk_inst.ws.h2osoi_vol_col)
end

# reverse via the TOP-LEVEL API (auto N, no hard-coded n_canopy).
db = C.clm_drv_reverse!(deepcopy(inst), bounds, filt, config;
        seed_bang! = seed!, include_canopy = true, use_psn = use_psn)
g_rev = db.inst.temperature.t_grnd_col[c0]

hs = (2e-2, 1e-2, 5e-3, 2.5e-3); cfd(h) = (L_pert(h) - L_pert(-h)) / (2h)
fds = [cfd(h) for h in hs]; for (h,f) in zip(hs,fds); @printf("  FD(h=%.2e) = % .8e\n", h, f); end
g_rich = (4*fds[end] - fds[end-1]) / 3; lo, hi = minimum(fds), maximum(fds)
bracketed = (lo - abs(lo)*1e-9) <= g_rev <= (hi + abs(hi)*1e-9)
relerr = abs(g_rev - g_rich) / max(abs(g_rich), 1e-12)
@printf("\n  Richardson = % .8e\n  clm_drv_reverse! dL/d(t_grnd) = % .8e\n", g_rich, g_rev)
@printf("  bracketed by FD spread [% .3e, % .3e]: %s\n", lo, hi, string(bracketed))
pass = relerr < 5e-5 || bracketed
println("\n", "="^70)
@printf("TOP-LEVEL clm_drv_reverse! (auto N=%d, use_psn=%s): rel=%.3e bracketed=%s  %s\n",
    Nauto, string(use_psn), relerr, string(bracketed), pass ? "PASS ✓" : "FAIL ✗")
println("="^70)
pass || error("clm_drv_reverse! FD mismatch")
