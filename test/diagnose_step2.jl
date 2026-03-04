using CLM, Dates, Printf

# Initialize
fsurdat = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/surfdata_clm.nc"
paramfile = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/clm5_params.nc"
fforcing = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/data/forcing/CLM_input/clmforc.2002.nc"

fhistory = tempname() * "_diag.nc"
ndays = 1

inst = CLM.clm_run!(;
    fsurdat=fsurdat, paramfile=paramfile, fforcing=fforcing,
    fhistory=fhistory,
    start_date=DateTime(2002, 1, 1),
    end_date=DateTime(2002, 1, 2),
    dtime=1800, use_cn=false, verbose=false)

println("\n=== POST-SIMULATION DIAGNOSTICS ===\n")

# Ground temperature
temp = inst.temperature
println("T_GRND:")
for c in 1:length(temp.t_grnd_col)
    @printf("  col %d: %.4f K\n", c, temp.t_grnd_col[c])
end

# Soil temps
println("\nSoil temperatures (first 5 layers):")
for c in 1:min(2, size(temp.t_soisno_col, 1))
    println("  col $c:")
    nlevsno = CLM.varpar.nlevsno
    for j in 1:5
        jj = nlevsno + j
        @printf("    layer %2d (idx %d): %.4f K\n", j, jj, temp.t_soisno_col[c, jj])
    end
end

# Friction velocity params
fv = inst.frictionvel
@printf("\nFrictionvel scalar params:\n")
@printf("  zetamaxstable = %.6f\n", fv.zetamaxstable)
@printf("  zsno          = %.6f\n", fv.zsno)
@printf("  zlnd          = %.6f\n", fv.zlnd)
@printf("  zglc          = %.6f\n", fv.zglc)

# z0mg_col
println("\nz0mg_col:")
for c in 1:length(fv.z0mg_col)
    @printf("  col %d: %.8f\n", c, fv.z0mg_col[c])
end

# z0mv_patch
println("\nz0mv_patch:")
for p in 1:length(fv.z0mv_patch)
    @printf("  patch %d: %.8f\n", p, fv.z0mv_patch[p])
end

# Fluxes
ef = inst.energyflux
println("\nEnergy fluxes:")
for p in 1:length(ef.eflx_sh_tot_patch)
    @printf("  patch %d: SH=%.4f, LH=%.4f, gnet=%.4f\n", p,
        ef.eflx_sh_tot_patch[p], ef.eflx_lh_tot_patch[p], ef.eflx_gnet_patch[p])
end

# Absorbed solar
sa = inst.solarabs
println("\nAbsorbed solar (FSA):")
for p in 1:length(sa.fsa_patch)
    @printf("  patch %d: %.4f W/m2\n", p, sa.fsa_patch[p])
end

# Canopy state
cs = inst.canopystate
println("\nfrac_veg_nosno_patch:")
for p in 1:length(cs.frac_veg_nosno_patch)
    @printf("  patch %d: %d (alb=%d)\n", p, cs.frac_veg_nosno_patch[p],
        cs.frac_veg_nosno_alb_patch[p])
end

# cgrnds/cgrndl
println("\ncgrnds_patch, cgrndl_patch:")
for p in 1:length(ef.cgrnds_patch)
    @printf("  patch %d: cgrnds=%.6f, cgrndl=%.6f\n", p,
        ef.cgrnds_patch[p], ef.cgrndl_patch[p])
end

# t_h2osfc, frac_h2osfc, frac_sno
wdb = inst.water.waterdiagnosticbulk_inst
println("\nWater diagnostics:")
for c in 1:length(wdb.frac_sno_col)
    @printf("  col %d: frac_sno=%.4f, frac_h2osfc=%.4f, t_h2osfc=%.4f\n",
        c, wdb.frac_sno_col[c], wdb.frac_h2osfc_col[c], temp.t_h2osfc_col[c])
end

# soilbeta, smpmin
ss = inst.soilstate
println("\nSoil state:")
for c in 1:length(ss.soilbeta_col)
    @printf("  col %d: soilbeta=%.4f, smpmin=%.2e\n", c, ss.soilbeta_col[c], ss.smpmin_col[c])
end

# Check for NaN in key fields
println("\n=== NaN CHECK ===")
function check_nan(name, arr)
    n = count(isnan, arr)
    if n > 0
        println("  $name: $n NaN out of $(length(arr))")
    end
end
check_nan("t_grnd_col", temp.t_grnd_col)
check_nan("t_soisno_col", temp.t_soisno_col)
check_nan("eflx_sh_tot_patch", ef.eflx_sh_tot_patch)
check_nan("eflx_lh_tot_patch", ef.eflx_lh_tot_patch)
check_nan("cgrnds_patch", ef.cgrnds_patch)
check_nan("fsa_patch", sa.fsa_patch)
check_nan("z0mv_patch", fv.z0mv_patch)
check_nan("z0mg_col", fv.z0mg_col)

wfb = inst.water.waterfluxbulk_inst
check_nan("qflx_evap_tot_col", wfb.wf.qflx_evap_tot_col)
check_nan("qflx_runoff_col", wfb.wf.qflx_runoff_col)
check_nan("qflx_ev_snow_patch(bulk)", wfb.qflx_ev_snow_patch)
check_nan("qflx_ev_soil_patch(bulk)", wfb.qflx_ev_soil_patch)
check_nan("qflx_ev_h2osfc_patch(bulk)", wfb.qflx_ev_h2osfc_patch)
check_nan("ustar_patch", fv.ustar_patch)
check_nan("obu_patch", fv.obu_patch)

rm(fhistory, force=true)
