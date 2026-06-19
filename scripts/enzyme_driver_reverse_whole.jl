# =============================================================================
# WHOLE-STEP HydrologyNoDrainage reverse through the PRODUCTION assembler
# CLM.driver_rev_phases — the full faithful pre-soil_water block now threaded in:
#   soil_temperature! → [set_soil_water_fractions! → set_floodc! → saturated_excess_runoff!
#   → set_qflx_inputs! → infiltration_excess_runoff! → route_infiltration_excess!
#   → update_h2osfc! → infiltration! → total_surface_runoff!
#   → compute_effec_rootfrac_and_vert_tran_sink!] → soil_water! → water_table!
#   → hydrology_no_drainage!
# (14 phases, no canopy.) Perturb t_grnd, seed h2osoi_vol, FD-validate dL/d(t_grnd)
# through the WHOLE chain. Uses the production driver_rev_bundle/driver_rev_phases +
# compositional_reverse! engine (no script-local wrappers).
#   julia +1.10 --project=/tmp/clm_jl10_<id> scripts/enzyme_driver_reverse_whole.jl
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
phases = C.driver_rev_phases(bounds, filt, config)   # 14 phases (no canopy)
hc = filt.hydrologyc; c0 = [c for c in bounds.begc:bounds.endc if hc[c]][1]
@printf("assembled %d-phase whole-step hydrology chain\n", length(phases))

Lvol(b) = sum(abs2, b.inst.water.waterstatebulk_inst.ws.h2osoi_vol_col)
function L_pert(δ)
    bs = C.driver_rev_bundle(deepcopy(inst)); bs.inst.temperature.t_grnd_col[c0] += δ
    for (f, ca) in phases; f(bs, ca...); end
    return Lvol(bs)
end
let bs = C.driver_rev_bundle(deepcopy(inst))
    for (f, ca) in phases; f(bs, ca...); end
    @printf("chained primal: h2osoi_vol[%d,4]=%.6f qflx_infl[%d]=%.3e (finite: %s)\n", c0,
        bs.inst.water.waterstatebulk_inst.ws.h2osoi_vol_col[c0,4], c0,
        bs.inst.water.waterfluxbulk_inst.wf.qflx_infl_col[c0],
        string(isfinite(bs.inst.water.waterstatebulk_inst.ws.h2osoi_vol_col[c0,4])))
end

seed!(db, b) = (db.inst.water.waterstatebulk_inst.ws.h2osoi_vol_col .=
                2 .* b.inst.water.waterstatebulk_inst.ws.h2osoi_vol_col)
db = C.compositional_reverse!(phases, C.driver_rev_bundle(deepcopy(inst)), seed!)
g_rev = db.inst.temperature.t_grnd_col[c0]

hs = (2e-2, 1e-2, 5e-3, 2.5e-3); cfd(h) = (L_pert(h) - L_pert(-h)) / (2h)
fds = [cfd(h) for h in hs]; for (h,f) in zip(hs,fds); @printf("  FD(h=%.2e) = % .8e\n", h, f); end
g_rich = (4*fds[end] - fds[end-1]) / 3; lo, hi = minimum(fds), maximum(fds)
bracketed = (lo - abs(lo)*1e-9) <= g_rev <= (hi + abs(hi)*1e-9)
relerr = abs(g_rev - g_rich) / max(abs(g_rich), 1e-12)
@printf("\n  Richardson = % .8e\n  rev dL/d(t_grnd) = % .8e\n", g_rich, g_rev)
@printf("  bracketed by FD spread [% .3e, % .3e]: %s\n", lo, hi, string(bracketed))
pass = relerr < 5e-5 || bracketed
println("\n", "="^70)
@printf("WHOLE-STEP HYDROLOGY REVERSE (14 phases): rel=%.3e bracketed=%s  %s\n",
    relerr, string(bracketed), pass ? "PASS ✓" : "FAIL ✗")
println("="^70)
