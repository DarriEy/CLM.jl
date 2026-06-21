# Probe: does the growing-season warming recipe (LAI + soil moisture + June-solstice sun)
# produce NONZERO photosynthesis (fpsn) at exposed patches? Mirrors build_cn_inst's recipe
# (enzyme_bgc_wholestep.jl) but use_cn=false (SP/prescribed-LAI canopy psn path).
using CLM, Printf
const C = CLM
const FSURDAT  = get(ENV, "CLM_FSURDAT",  "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/surfdata_clm.nc")
const PARAMFILE = get(ENV, "CLM_PARAMFILE", "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/clm5_params.nc")

(inst, bounds, filt, tm) = C.clm_initialize!(; fsurdat=FSURDAT, paramfile=PARAMFILE)
config = C.CLMDriverConfig(); fia = C.clump_filter_inactive_and_active
dtime = 1800.0; calday = 172.5                         # June 21 — high sun
(declin, _) = C.compute_orbital(calday); nextsw = calday + dtime/C.SECSPDAY
ng = bounds.endg
C._setup_calib_forcing!(inst.atm2lnd, 295.0, ng)        # warm growing-season air temp
C.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
C._init_calib_soil_moisture!(inst, bounds)
C.interp_monthly_veg!(inst.satellite_phenology; kmo=7, kda=15)   # July LAI
cs = inst.canopystate; wdb = inst.water.waterdiagnosticbulk_inst; pch = inst.patch
C.satellite_phenology!(inst.satellite_phenology, cs, wdb, pch, filt.nolakep, bounds.begp:bounds.endp)
for p in bounds.begp:bounds.endp; cs.frac_veg_nosno_patch[p] = cs.frac_veg_nosno_alb_patch[p]; end
C.set_exposedvegp_filter!(filt, bounds, cs.frac_veg_nosno_patch)

ev = filt.exposedvegp; expp = Int[p for p in bounds.begp:bounds.endp if ev[p]]
@printf("exposed-veg patches: %d   (tlai at them: %s)\n", length(expp),
    string(round.([cs.tlai_patch[p] for p in expp[1:min(end,4)]], digits=3)))

for n in 1:4
    C.clm_drv!(config, inst, filt, fia, bounds, true, nextsw, declin, declin,
        C.ORB_OBLIQR_DEFAULT, false, false, "", false; nstep=n, is_first_step=(n==1),
        is_beg_curr_day=(n==1), dtime=dtime, mon=7, day=15, photosyns=inst.photosyns)
    fp = isempty(expp) ? 0 : expp[1]
    @printf("step %d: fpsn[exposed]=%s  parsun_z[%d,1]=%.4f  btran[%d]=%.4f  t_veg[%d]=%.3f\n",
        n, isempty(expp) ? "—" : string(round.([inst.photosyns.fpsn_patch[p] for p in expp[1:min(end,3)]], digits=4)),
        fp, fp==0 ? NaN : inst.solarabs.parsun_z_patch[fp,1],
        fp, fp==0 ? NaN : inst.energyflux.btran_patch[fp],
        fp, fp==0 ? NaN : inst.temperature.t_veg_patch[fp])
end
