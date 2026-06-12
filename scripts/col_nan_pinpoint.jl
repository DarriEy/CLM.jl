# Pinpoint the exact contaminated column: after init + 1 step, for each column print
# landunit itype, snl, and first-NaN layer of t_soisno / csol / col.z. Also verify the
# soil-prop bedrock fix took (csol[c,21:25] finite).
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

col = inst.column; ss = inst.soilstate; tp = inst.temperature
nc = bounds.endc; nl = bounds.endl
println("nc=$nc  nl=$nl  ISTDLAK=", CLM.ISTDLAK)
# verify soil-prop bedrock fix: csol[c,21:25]
println("\n--- AFTER INIT (pre-step) ---")
println("csol[c,21:25] finite per col (fix check): ",
    [all(isfinite, ss.csol_col[c, 21:25]) for c in 1:nc])
println("watsat[c,21:25] finite per col: ",
    [all(isfinite, ss.watsat_col[c, 21:25]) for c in 1:nc])

firstnan(M, c) = (for j in 1:size(M,2); isnan(M[c,j]) && return j; end; 0)
function dump(tag)
    println("\n--- $tag ---")
    for c in 1:nc
        l = col.landunit[c]
        lit = (l >= 1 && l <= length(inst.landunit.itype)) ? inst.landunit.itype[l] : -99
        @printf("  col %2d  lun=%d itype=%2d  snl=%d | t_soisno firstNaN=%2d  csol firstNaN=%2d  z firstNaN=%2d  watsat firstNaN=%2d\n",
            c, l, lit, col.snl[c], firstnan(tp.t_soisno_col, c), firstnan(ss.csol_col, c),
            firstnan(col.z, c), firstnan(ss.watsat_col, c))
    end
end
dump("after init")
CLM.clm_drv!(config, inst, filt, filt_ia, bounds, true, nextsw_cday, declin, declin,
    CLM.ORB_OBLIQR_DEFAULT, false, false, "", false; nstep=1, is_first_step=true,
    is_beg_curr_day=true, dtime=1800.0, mon=1, day=1, photosyns=inst.photosyns)
dump("after 1 step")
println("\nexposed-veg patches: ", [p for p in bounds.begp:bounds.endp if filt.exposedvegp[p]],
        " -> columns ", [col.gridcell[inst.patch.column[p]] === nothing ? -1 : inst.patch.column[p]
                          for p in bounds.begp:bounds.endp if filt.exposedvegp[p]])
println("landunit itypes (ISTSOIL/ISTDLAK/etc): ", [inst.landunit.itype[l] for l in 1:nl])
println("CLM consts: ISTSOIL=", CLM.ISTSOIL, " ISTDLAK=", CLM.ISTDLAK,
        " ISTICE=", (isdefined(CLM, :ISTICE) ? CLM.ISTICE : "?"))
