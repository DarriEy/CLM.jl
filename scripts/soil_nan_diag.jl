# Pinpoint the soil-state NaN feeding albsod + btran. One real-init + 3 warmup steps,
# then dump soil water / porosity / color for the exposed-veg columns.
using CLM, Printf
const FSURDAT = get(ENV, "CLM_FSURDAT",
    "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/surfdata_clm.nc")
const PARAMFILE = get(ENV, "CLM_PARAMFILE",
    "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/clm5_params.nc")
(inst, bounds, filt, tm) = CLM.clm_initialize!(; fsurdat=FSURDAT, paramfile=PARAMFILE)
ng = bounds.endg
config = CLM.CLMDriverConfig(); filt_ia = CLM.clump_filter_inactive_and_active
(declin, eccf) = CLM.compute_orbital(120.0); nextsw_cday = 120.0 + 1800.0/CLM.SECSPDAY
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

ws = inst.water.waterstatebulk_inst.ws
wd = inst.water.waterdiagnosticbulk_inst
ss = inst.soilstate
ev = filt.exposedvegp
cols = unique([inst.patch.column[p] for p in bounds.begp:bounds.endp if ev[p]])
nlevgrnd = CLM.varpar.nlevgrnd; nlevsno = CLM.varpar.nlevsno

function dump(tag)
    println("\n--- $tag ---  exposed-veg columns = $cols")
    for c in cols
        snl = inst.column.snl[c]; j0 = nlevsno + 1   # top soil layer index (combined snow+soil)
        @printf("  col %d: snl=%d\n", c, snl)
        @printf("    h2osoi_vol[1..3]    = %s\n", string([round(ws.h2osoi_vol_col[c,j], sigdigits=5) for j in 1:3]))
        @printf("    h2osoi_liq[j0..j0+2]= %s\n", string([round(ws.h2osoi_liq_col[c,j], sigdigits=5) for j in j0:j0+2]))
        @printf("    h2osoi_ice[j0..j0+2]= %s\n", string([round(ws.h2osoi_ice_col[c,j], sigdigits=5) for j in j0:j0+2]))
        @printf("    watsat[1..3]        = %s\n", string([round(ss.watsat_col[c,j], sigdigits=5) for j in 1:3]))
        @printf("    h2osoi_liqvol[1..3] = %s\n", string([round(wd.h2osoi_liqvol_col[c,j], sigdigits=5) for j in 1:3]))
    end
end

dump("after init (pre-step)")
for n in 1:3
    CLM.clm_drv!(config, inst, filt, filt_ia, bounds, true, nextsw_cday, declin, declin,
        CLM.ORB_OBLIQR_DEFAULT, false, false, "", false; nstep=n, is_first_step=(n==1),
        is_beg_curr_day=(n==1), dtime=1800.0, mon=1, day=1, photosyns=inst.photosyns)
end
dump("after 3 warmup steps")

println("\n  btran_patch (exposed) = ", [round(inst.energyflux.btran_patch[p], sigdigits=5) for p in bounds.begp:bounds.endp if ev[p]])
nc = bounds.endc
firstnan(M, c) = (for j in 1:size(M,2); isnan(M[c,j]) && return j; end; 0)
println("\n  per-column first-NaN layer (0 = no NaN):")
@printf("  %-22s %s\n", "field", "first-NaN-layer per col 1:$nc")
for (nm, M) in (("ss.bsw_col", ss.bsw_col), ("ss.watsat_col", ss.watsat_col),
                ("ss.sucsat_col", ss.sucsat_col), ("ss.hksat_col", ss.hksat_col),
                ("ss.smp_l_col", ss.smp_l_col), ("ss.soilpsi_col", ss.soilpsi_col),
                ("ws.h2osoi_liq", ws.h2osoi_liq_col), ("ws.h2osoi_ice", ws.h2osoi_ice_col),
                ("ws.h2osoi_vol", ws.h2osoi_vol_col), ("wd.h2osoi_liqvol", wd.h2osoi_liqvol_col))
    @printf("  %-22s %s\n", nm, string([firstnan(M, c) for c in 1:nc]))
end
println("\n  nbedrock per col = ", hasproperty(inst.column, :nbedrock) ? string(inst.column.nbedrock) : "(no field)")
println("  varpar nlevgrnd=$nlevgrnd nlevsoi=$(CLM.varpar.nlevsoi) nlevsno=$nlevsno")
