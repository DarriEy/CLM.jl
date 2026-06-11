# Phase probe: does canopy_fluxes_core! reverse-mode compile STANDALONE?
#
# The whole-driver probe hits an Enzyme-internal sret_ty bug on canopy_fluxes_core!
# via runtime_generic_augfwd (Enzyme's generic-dispatch path). A compositional
# approach differentiates each phase with a SEPARATE, DIRECT Enzyme.autodiff call
# — which does NOT go through runtime_generic_augfwd. This probe tests whether
# canopy_fluxes_core! differentiates on its own (the viability test for
# compositional reverse-AD over the driver). Run on Julia 1.10.
#
#   julia +1.10 --project=/tmp/clm_jl10_env scripts/enzyme_phase_probe_canopy.jl
#
# RESULT (2026-06-10, Julia 1.10.11): STILL hits the Enzyme-internal sret_ty bug
# (propagate_returned!) even STANDALONE — so the bug is intrinsic to
# differentiating canopy_fluxes_core! under this Enzyme, NOT a dispatch-path
# artifact. ALSO tested the IN-PLACE form (Const return, seed output shadow=1,
# read input shadows — the banked Phase-C workaround for scalar-return bugs): it
# fails with the IDENTICAL sret_ty error too. So the bug is return-form-independent
# (3 tests: whole-driver dispatch, standalone Active-return, standalone in-place
# Const-return all fail the same way). Compositional AD therefore CANNOT get past
# this phase; it needs an Enzyme fix, a hand-written EnzymeRule, or sub-phase
# decomposition. (soil_temperature! DOES differentiate cleanly — see
# enzyme_phase_probe.jl — so compositional works for the phases Enzyme can handle.)
using CLM, Enzyme, Printf
const FSURDAT = get(ENV, "CLM_FSURDAT",
    "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/surfdata_clm.nc")
const PARAMFILE = get(ENV, "CLM_PARAMFILE",
    "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/clm5_params.nc")

(inst, bounds, filt, tm) = CLM.clm_initialize!(; fsurdat=FSURDAT, paramfile=PARAMFILE)
ng = bounds.endg; nc = bounds.endc; np = bounds.endp
config = CLM.CLMDriverConfig()
filt_ia = CLM.clump_filter_inactive_and_active
(declin, eccf) = CLM.compute_orbital(120.0)
nextsw_cday = 120.0 + 1800.0 / CLM.SECSPDAY
T0 = 285.0
a2l = inst.atm2lnd
for g in 1:ng
    a2l.forc_t_not_downscaled_grc[g] = T0
    a2l.forc_pbot_not_downscaled_grc[g] = 85000.0
    a2l.forc_th_not_downscaled_grc[g] = T0 * (100000.0/85000.0)^(CLM.RAIR/CLM.CPAIR)
    a2l.forc_rho_not_downscaled_grc[g] = 85000.0/(CLM.RAIR*T0)
    a2l.forc_lwrad_not_downscaled_grc[g] = 300.0
    a2l.forc_vp_grc[g] = 800.0; a2l.forc_hgt_grc[g] = 30.0; a2l.forc_topo_grc[g] = 0.0
    a2l.forc_wind_grc[g] = 3.0
    for b in 1:CLM.NUMRAD
        a2l.forc_solad_not_downscaled_grc[g, b] = 200.0; a2l.forc_solai_grc[g, b] = 80.0
    end
    a2l.forc_solar_not_downscaled_grc[g] = 560.0
    a2l.forc_rain_not_downscaled_grc[g] = 0.0001; a2l.forc_snow_not_downscaled_grc[g] = 0.0
end
CLM.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
step!(n) = CLM.clm_drv!(config, inst, filt, filt_ia, bounds,
    true, nextsw_cday, declin, declin, CLM.ORB_OBLIQR_DEFAULT, false, false, "", false;
    nstep=n, is_first_step=(n==1), is_beg_curr_day=(n==1),
    dtime=1800.0, mon=1, day=1, photosyns=inst.photosyns)
println("Warmup (3 steps)..."); for n in 1:3; step!(n); end

# Constants for the canopy_fluxes_core! call.
bc_col = bounds.begc:bounds.endc
bc_patch = bounds.begp:bounds.endp
dtime = 1800.0
forc_q_col = zeros(nc)
CLM.compute_forc_q!(forc_q_col, inst.column.gridcell, a2l.forc_vp_grc, a2l.forc_pbot_downscaled_col)
downreg = zeros(np); leafn = zeros(np); o3v = zeros(np); o3g = zeros(np)
pc = CLM.pftcon

# Phase as an inst-closure: differentiate canopy_fluxes_core! alone.
function phase!(inst)
    CLM.canopy_fluxes_core!(inst.canopystate, inst.energyflux, inst.frictionvel,
        inst.temperature, inst.solarabs, inst.soilstate,
        inst.water.waterfluxbulk_inst, inst.water.waterstatebulk_inst,
        inst.water.waterdiagnosticbulk_inst, inst.photosyns,
        inst.patch, inst.column, inst.gridcell,
        filt.exposedvegp, bc_patch, bc_col,
        a2l.forc_lwrad_downscaled_col, forc_q_col, a2l.forc_pbot_downscaled_col,
        a2l.forc_th_downscaled_col, a2l.forc_rho_downscaled_col, a2l.forc_t_downscaled_col,
        a2l.forc_u_grc, a2l.forc_v_grc, a2l.forc_pco2_grc, a2l.forc_po2_grc,
        a2l.forc_hgt_t_grc, a2l.forc_hgt_u_grc, a2l.forc_hgt_q_grc,
        inst.gridcell.dayl, inst.gridcell.max_dayl, downreg, leafn, dtime,
        inst.temperature.t_a10_patch, inst.surfalb.nrad_patch, inst.surfalb.tlai_z_patch,
        inst.surfalb.vcmaxcintsun_patch, inst.surfalb.vcmaxcintsha_patch,
        inst.solarabs.parsun_z_patch, inst.solarabs.parsha_z_patch,
        inst.canopystate.laisun_z_patch, inst.canopystate.laisha_z_patch,
        o3v, o3g, pc.dleaf, pc.slatop, pc.leafcn, pc.flnr, pc.fnitr, pc.mbbopt,
        pc.c3psn, pc.woody, inst.overrides)
    return sum(abs2, inst.temperature.t_veg_patch)
end

println("phase! primal value = ", phase!(inst))
dinst = Enzyme.make_zero(inst)
println("Attempting Enzyme.Reverse on canopy_fluxes_core! alone...")
Enzyme.API.strictAliasing!(false)
flush(stdout)
try
    Enzyme.autodiff(Enzyme.set_runtime_activity(Enzyme.Reverse), phase!, Enzyme.Active,
                    Enzyme.Duplicated(inst, dinst))
    g = dinst.temperature.t_veg_patch
    @printf("\nSUCCESS: canopy_fluxes_core! reverse-mode compiled.\n")
    @printf("  nonzero t_veg shadow entries = %d / %d\n", count(!=(0.0), g), length(g))
catch e
    println("\nPHASE BLOCKER:")
    showerror(stdout, e, catch_backtrace())
    println()
end
