# Whole-driver reverse-AD on the REAL cold-start inst — now that canopy is finite
# (the forc_hgt_u/t/q_grc harness fix). Demonstrates: (1) the full clm_drv! timestep
# runs FINITE on the real Bow-at-Banff cold-start; (2) soil_temperature! (a real
# driver phase) reverse-differentiates through CLM.compositional_reverse! on the real
# inst and its gradient MATCHES central finite differences — the first FD-validated
# reverse pass on the real model state (impossible before: the inst was NaN).
#
#   julia +1.10 --project=/tmp/clm_jl10_env scripts/enzyme_realinst_reverse.jl
#
using CLM, Enzyme, Printf
const FSURDAT = get(ENV, "CLM_FSURDAT",
    "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/surfdata_clm.nc")
const PARAMFILE = get(ENV, "CLM_PARAMFILE",
    "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/clm5_params.nc")

function setup_forcing!(a2l, T, ng)
    for g in 1:ng
        a2l.forc_t_not_downscaled_grc[g]=T; a2l.forc_pbot_not_downscaled_grc[g]=85000.0
        a2l.forc_th_not_downscaled_grc[g]=T*(100000.0/85000.0)^(CLM.RAIR/CLM.CPAIR)
        a2l.forc_rho_not_downscaled_grc[g]=85000.0/(CLM.RAIR*T); a2l.forc_lwrad_not_downscaled_grc[g]=300.0
        a2l.forc_vp_grc[g]=800.0; a2l.forc_hgt_grc[g]=30.0; a2l.forc_topo_grc[g]=0.0; a2l.forc_wind_grc[g]=3.0
        # THE FIX: component obs heights (forcing_reader derives these; harness bypasses it).
        a2l.forc_hgt_u_grc[g]=30.0; a2l.forc_hgt_t_grc[g]=30.0; a2l.forc_hgt_q_grc[g]=30.0
        for b in 1:CLM.NUMRAD; a2l.forc_solad_not_downscaled_grc[g,b]=200.0; a2l.forc_solai_grc[g,b]=80.0; end
        a2l.forc_solar_not_downscaled_grc[g]=560.0; a2l.forc_rain_not_downscaled_grc[g]=0.0001; a2l.forc_snow_not_downscaled_grc[g]=0.0
    end
end

function build()
    (inst, bounds, filt, tm) = CLM.clm_initialize!(; fsurdat=FSURDAT, paramfile=PARAMFILE)
    setup_forcing!(inst.atm2lnd, 285.0, bounds.endg)
    CLM.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
    return inst, bounds, filt
end
const config = CLM.CLMDriverConfig(); const filt_ia = CLM.clump_filter_inactive_and_active
const (declin, eccf) = CLM.compute_orbital(120.0); const nextsw_cday = 120.0 + 1800.0/CLM.SECSPDAY
runstep!(inst, filt, bounds, n; first=false) = CLM.clm_drv!(config, inst, filt, filt_ia, bounds,
    true, nextsw_cday, declin, declin, CLM.ORB_OBLIQR_DEFAULT, false, false, "", false;
    nstep=n, is_first_step=first, is_beg_curr_day=first, dtime=1800.0, mon=1, day=1, photosyns=inst.photosyns)

inst, bounds, filt = build()
for n in 1:3; runstep!(inst, filt, bounds, n; first=(n==1)); end
ev = filt.exposedvegp; ps = [p for p in bounds.begp:bounds.endp if ev[p]]
@printf("Real cold-start clm_drv! finite check (exposed-veg patches %s):\n", string(ps))
@printf("  t_veg = %s\n", string([round(inst.temperature.t_veg_patch[p], digits=4) for p in ps]))
@printf("  sabv  = %s   btran = %s\n",
    string([round(inst.solarabs.sabv_patch[p], digits=3) for p in ps]),
    string([round(inst.energyflux.btran_patch[p], digits=4) for p in ps]))
@printf("  t_soisno all finite? %s   t_grnd all finite? %s\n",
    all(isfinite, inst.temperature.t_soisno_col), all(isfinite, inst.temperature.t_grnd_col))

# --- soil_temperature! as a driver phase over the whole-inst bundle ---
const URBANTV = fill(323.15, bounds.endl); const DTIME = 1800.0
const BC_COL = bounds.begc:bounds.endc; const BC_LUN = bounds.begl:bounds.endl; const BC_PATCH = bounds.begp:bounds.endp
function soil_phase!(b)
    CLM.soil_temperature!(b.column, b.landunit, b.patch, b.temperature, b.energyflux,
        b.soilstate, b.water.waterstatebulk_inst, b.water.waterdiagnosticbulk_inst,
        b.water.waterfluxbulk_inst, b.solarabs, b.canopystate, b.urbanparams,
        URBANTV, b.atm2lnd.forc_lwrad_downscaled_col, filt.nolakec, filt.nolakep,
        filt.urbanl, filt.urbanc, BC_COL, BC_LUN, BC_PATCH, DTIME)
    return nothing
end
L(inst) = sum(abs2, inst.temperature.t_soisno_col)

# pick a finite soil layer to perturb (top soil of an exposed column)
c0 = inst.patch.column[ps[1]]; j0 = CLM.varpar.nlevsno + 1
@printf("\nValidating dL/d(t_soisno[%d,%d]) (L = sum t_soisno^2):\n", c0, j0)

# central finite differences of the soil_temperature! primal
function L_pert(δ)
    s = deepcopy(inst); s.temperature.t_soisno_col[c0, j0] += δ
    soil_phase!(s); return L(s)
end
h = 1e-3
g_fd = (L_pert(h) - L_pert(-h)) / (2h)
@printf("  FD   = % .8e\n", g_fd)

# compositional reverse through the shared engine
seed_soil!(db, b) = (db.temperature.t_soisno_col .= 2 .* b.temperature.t_soisno_col)
db = CLM.compositional_reverse!(Any[(soil_phase!, ())], deepcopy(inst), seed_soil!)
g_rev = db.temperature.t_soisno_col[c0, j0]
@printf("  comp = % .8e\n", g_rev)

println("\n", "="^60)
relerr = abs(g_rev - g_fd) / max(abs(g_fd), 1e-10)
@printf("rel error = %.3e\n", relerr)
println(isfinite(g_rev) && isfinite(g_fd) && relerr < 1e-4 ?
    "\nREAL-INST WHOLE-DRIVER-PHASE REVERSE-AD VALIDATED vs FD ✓\n(soil_temperature! on the finite real cold-start, via CLM.compositional_reverse!)" :
    "\nMISMATCH ✗")
println("="^60)
