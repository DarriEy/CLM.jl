# Find which per-patch gnet flux input is NaN going INTO step 1 (the root of the col-1
# soil-temperature NaN). Dump per-patch flux state AFTER INIT (pre-step) — that's the
# state compute_ground_heat_flux_and_deriv! reads on step 1.
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

pd = inst.patch; ef = inst.energyflux; sa = inst.solarabs
wf = inst.water.waterfluxbulk_inst
np = bounds.endp
nan(x) = isnan(x) ? "NaN" : @sprintf("%.4g", x)
# patch weight field
wt = hasproperty(pd, :wtcol) ? pd.wtcol : (hasproperty(pd, :wtgcell) ? pd.wtgcell : fill(-1.0, np))
gv(obj, f, p) = hasproperty(obj, f) ? nan(getproperty(obj, f)[p]) : "-"
gvw(p) = begin  # waterflux nested .wf
    w = inst.water.waterfluxbulk_inst
    (hasproperty(w, :wf) ? w.wf : w)
end
println("Per-patch flux state AFTER INIT (np=$np), col=patch's column, exposed=filt.exposedvegp:")
@printf("%3s %4s %5s %5s %4s | %8s %8s %8s %8s %8s %8s %8s\n",
    "p","col","ivt","wt","exp","cgrnd","cgrnds","cgrndl","dlrad","sabg_soil","sh_soil","ev_soil")
for p in 1:np
    c = pd.column[p]
    wfb = gvw(p)
    @printf("%3d %4d %5d %5.2f %4s | %8s %8s %8s %8s %8s %8s %8s\n",
        p, c, pd.itype[p], wt[p], string(filt.exposedvegp[p]),
        gv(ef, :cgrnd_patch, p), gv(ef, :cgrnds_patch, p), gv(ef, :cgrndl_patch, p),
        gv(ef, :dlrad_patch, p), gv(sa, :sabg_soil_patch, p),
        gv(ef, :eflx_sh_soil_patch, p),
        hasproperty(wfb, :qflx_ev_soil_patch) ? nan(wfb.qflx_ev_soil_patch[p]) : "-")
end
println("\npatch active mask: ", [pd.active[p] for p in 1:np])

# Run ONE step, then dump per-patch gnet inputs AS LEFT BY THE TIMESTEP physics.
config = CLM.CLMDriverConfig(); filt_ia = CLM.clump_filter_inactive_and_active
(declin, eccf) = CLM.compute_orbital(120.0); nextsw_cday = 120.0 + 1800.0/CLM.SECSPDAY
CLM.clm_drv!(config, inst, filt, filt_ia, bounds, true, nextsw_cday, declin, declin,
    CLM.ORB_OBLIQR_DEFAULT, false, false, "", false; nstep=1, is_first_step=true,
    is_beg_curr_day=true, dtime=1800.0, mon=1, day=1, photosyns=inst.photosyns)
println("\nPer-patch flux state AFTER 1 STEP (which input is NaN going into gnet?):")
@printf("%3s %4s %5s %5s %4s | %8s %8s %8s %8s %8s %8s %8s %8s\n",
    "p","col","ivt","wt","exp","cgrnd","dgnetdT","eflx_gnet","dlrad","sabg_soil","sh_soil","ev_soil","sh_grnd")
for p in 1:np
    c = pd.column[p]; wfb = gvw(p)
    @printf("%3d %4d %5d %5.2f %4s | %8s %8s %8s %8s %8s %8s %8s %8s\n",
        p, c, pd.itype[p], wt[p], string(filt.exposedvegp[p]),
        gv(ef, :cgrnd_patch, p), gv(ef, :dgnetdT_patch, p), gv(ef, :eflx_gnet_patch, p),
        gv(ef, :dlrad_patch, p), gv(sa, :sabg_soil_patch, p), gv(ef, :eflx_sh_soil_patch, p),
        hasproperty(wfb, :qflx_ev_soil_patch) ? nan(wfb.qflx_ev_soil_patch[p]) : "-",
        gv(ef, :eflx_sh_grnd_patch, p))
end
println("frac_veg_nosno per patch: ", [gv(inst.canopystate, :frac_veg_nosno_patch, p) for p in 1:np])
println("exposedvegp(step1): ", [p for p in 1:np if filt.exposedvegp[p]],
        "  noexposedvegp(step1): ", [p for p in 1:np if filt.noexposedvegp[p]])
