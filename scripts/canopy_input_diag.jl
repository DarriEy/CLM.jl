# Characterize the canopy Newton divergence (cgrnd ~1e15 for veg patches): dump the key
# canopy INPUTS (btran, sabv, friction/aero state) for the exposed-veg patches after init
# and after 1 step, to find the degenerate input.
using CLM, Printf
const FSURDAT = get(ENV, "CLM_FSURDAT",
    "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/surfdata_clm.nc")
const PARAMFILE = get(ENV, "CLM_PARAMFILE",
    "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/clm5_params.nc")
(inst, bounds, filt, tm) = CLM.clm_initialize!(; fsurdat=FSURDAT, paramfile=PARAMFILE)
ng = bounds.endg
a2l = inst.atm2lnd; T0 = 285.0
for g in 1:ng
    a2l.forc_t_not_downscaled_grc[g]=T0; a2l.forc_pbot_not_downscaled_grc[g]=85000.0
    a2l.forc_th_not_downscaled_grc[g]=T0*(100000.0/85000.0)^(CLM.RAIR/CLM.CPAIR)
    a2l.forc_rho_not_downscaled_grc[g]=85000.0/(CLM.RAIR*T0); a2l.forc_lwrad_not_downscaled_grc[g]=300.0
    a2l.forc_vp_grc[g]=800.0; a2l.forc_hgt_grc[g]=30.0; a2l.forc_topo_grc[g]=0.0; a2l.forc_wind_grc[g]=3.0
    for b in 1:CLM.NUMRAD; a2l.forc_solad_not_downscaled_grc[g,b]=200.0; a2l.forc_solai_grc[g,b]=80.0; end
    a2l.forc_solar_not_downscaled_grc[g]=560.0; a2l.forc_rain_not_downscaled_grc[g]=0.0001; a2l.forc_snow_not_downscaled_grc[g]=0.0
end
CLM.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)

pd = inst.patch; ef = inst.energyflux; sa = inst.solarabs; fv = inst.frictionvel; cs = inst.canopystate
nan(x) = isnan(x) ? "NaN" : @sprintf("%.4g", x)
gv(o,f,p) = hasproperty(o,f) ? nan(getproperty(o,f)[p]) : "-"
function dump(tag, ps)
    println("\n--- $tag (patches $ps) ---")
    @printf("%3s | %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s\n",
        "p","btran","sabv","elai","htop","z0mv","displa","ustar","um","obu","ram1")
    for p in ps
        @printf("%3d | %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s\n", p,
            gv(ef,:btran_patch,p), gv(sa,:sabv_patch,p), gv(cs,:elai_patch,p), gv(cs,:htop_patch,p),
            gv(fv,:z0mv_patch,p), gv(cs,:displa_patch,p), gv(fv,:ustar_patch,p), gv(fv,:um_patch,p),
            gv(fv,:obu_patch,p), gv(fv,:ram1_patch,p))
    end
end
allp = collect(bounds.begp:bounds.endp)
dump("AFTER INIT", allp)
config = CLM.CLMDriverConfig(); filt_ia = CLM.clump_filter_inactive_and_active
(declin, eccf) = CLM.compute_orbital(120.0); nextsw_cday = 120.0 + 1800.0/CLM.SECSPDAY
CLM.clm_drv!(config, inst, filt, filt_ia, bounds, true, nextsw_cday, declin, declin,
    CLM.ORB_OBLIQR_DEFAULT, false, false, "", false; nstep=1, is_first_step=true,
    is_beg_curr_day=true, dtime=1800.0, mon=1, day=1, photosyns=inst.photosyns)
exp_ps = [p for p in allp if filt.exposedvegp[p]]
dump("AFTER 1 STEP (exposed veg)", exp_ps)
println("\ncgrnds/cgrndl (veg): ", [(p, gv(ef,:cgrnds_patch,p), gv(ef,:cgrndl_patch,p)) for p in exp_ps])
println("uaf/rah1 (veg): ", [(p, gv(fv,:uaf_patch,p), gv(fv,:rah1_patch,p)) for p in exp_ps])
