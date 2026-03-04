using CLM, Dates, Printf

fsurdat = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/surfdata_clm.nc"
paramfile = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/clm5_params.nc"
fforcing = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/data/forcing/CLM_input/clmforc.2002.nc"
fhistory = tempname() * "_diag.nc"

# Run just 2 timesteps (1 hour)
inst = CLM.clm_run!(;
    fsurdat=fsurdat, paramfile=paramfile, fforcing=fforcing,
    fhistory=fhistory,
    start_date=DateTime(2002, 1, 1),
    end_date=DateTime(2002, 1, 1, 0, 30, 0),  # just 30 min = 1 timestep
    dtime=1800, use_cn=false, verbose=false)

println("\n=== AFTER 1 TIMESTEP ===\n")

ef = inst.energyflux
fv = inst.frictionvel
temp = inst.temperature
cs = inst.canopystate
wfb = inst.water.waterfluxbulk_inst
wdb = inst.water.waterdiagnosticbulk_inst

for p in 1:length(ef.cgrnds_patch)
    @printf("patch %d: cgrnds=%12.4f cgrndl=%12.4f ustar=%10.6f obu=%12.4f z0mv=%10.6f fvn=%d\n",
        p, ef.cgrnds_patch[p], ef.cgrndl_patch[p], fv.ustar_patch[p],
        fv.obu_patch[p], fv.z0mv_patch[p], cs.frac_veg_nosno_patch[p])
end
println()
for c in 1:length(temp.t_grnd_col)
    nlevsno = CLM.varpar.nlevsno
    @printf("col %d: t_grnd=%.4f z0mg=%.8f t_soi1=%.4f t_soi2=%.4f\n",
        c, temp.t_grnd_col[c], fv.z0mg_col[c],
        temp.t_soisno_col[c, nlevsno+1], temp.t_soisno_col[c, nlevsno+2])
end
println()

# Check fluxes
for p in 1:length(ef.eflx_sh_tot_patch)
    @printf("patch %d: SH_tot=%.4f SH_grnd=%.4f SH_veg=%.4f LH_tot=%.4f gnet=%.4f\n",
        p, ef.eflx_sh_tot_patch[p], ef.eflx_sh_grnd_patch[p],
        ef.eflx_sh_veg_patch[p], ef.eflx_lh_tot_patch[p], ef.eflx_gnet_patch[p])
end
println()

# Water flux bulk
for p in 1:length(wfb.qflx_ev_snow_patch)
    @printf("patch %d: ev_snow=%.6f ev_soil=%.6f ev_h2osfc=%.6f\n",
        p, wfb.qflx_ev_snow_patch[p], wfb.qflx_ev_soil_patch[p], wfb.qflx_ev_h2osfc_patch[p])
end
println()

# Water diagnostics
for c in 1:length(wdb.qg_col)
    @printf("col %d: qg=%.8f qg_snow=%.8f qg_soil=%.8f dqgdT=%.8f\n",
        c, wdb.qg_col[c], wdb.qg_snow_col[c], wdb.qg_soil_col[c], wdb.dqgdT_col[c])
end
println()

# t_veg
for p in 1:length(temp.t_veg_patch)
    @printf("patch %d: t_veg=%.4f\n", p, temp.t_veg_patch[p])
end

# Forcing values
a2l = inst.atm2lnd
println("\nForcing:")
@printf("  forc_t:     %.4f K\n", a2l.forc_t_downscaled_col[1])
@printf("  forc_lwrad: %.4f W/m2\n", a2l.forc_lwrad_downscaled_col[1])
@printf("  forc_pbot:  %.2f Pa\n", a2l.forc_pbot_downscaled_col[1])
@printf("  forc_rho:   %.6f kg/m3\n", a2l.forc_rho_downscaled_col[1])
@printf("  forc_solad: [%.4f, %.4f]\n", a2l.forc_solad_downscaled_col[1,1], a2l.forc_solad_downscaled_col[1,2])

# Absorbed solar
sa = inst.solarabs
@printf("\nFSA patch 1: %.4f\n", sa.fsa_patch[1])
@printf("sabg_patch 1: %.4f\n", sa.sabg_patch[1])

# dlrad, ulrad
@printf("\ndlrad: %.4f\n", ef.dlrad_patch[1])
@printf("ulrad: %.4f\n", ef.ulrad_patch[1])

# htvp
@printf("\nhtvp_col[1]: %.4f\n", ef.htvp_col[1])
@printf("emg_col[1]: %.4f\n", temp.emg_col[1])

# Soil thermal
nlevsno = CLM.varpar.nlevsno
nlevgrnd = CLM.varpar.nlevgrnd
ss = inst.soilstate
wsb = inst.water.waterstatebulk_inst
println("\nSoil properties (col 1, first 5 layers):")
for j in 1:5
    @printf("  layer %d: watsat=%.4f h2osoi_liq=%.4f h2osoi_ice=%.4f dz=%.6f\n",
        j, ss.watsat_col[1, j],
        wsb.ws.h2osoi_liq_col[1, nlevsno+j],
        wsb.ws.h2osoi_ice_col[1, nlevsno+j],
        inst.column.dz[1, nlevsno+j])
end

# Compute cv and fact for top layer using soil_therm_prop!
println("\nComputing cv via soil_therm_prop!:")
nc = length(temp.t_grnd_col)
nlevmaxurbgrnd = CLM.varpar.nlevmaxurbgrnd
tk_arr = zeros(nc, nlevsno + nlevmaxurbgrnd)
cv_arr = zeros(nc, nlevsno + nlevmaxurbgrnd)
tk_h2osfc = fill(NaN, nc)
CLM.soil_therm_prop!(inst.column, inst.landunit, inst.urbanparams, temp, wsb, wdb, ss,
                     filt.nolakec, filt.urbanc, bounds.begc:bounds.endc,
                     tk_arr, cv_arr, tk_h2osfc)
for j in 1:5
    jj = nlevsno + j
    @printf("  layer %d (jj=%d): cv=%.4f tk=%.6f\n", j, jj, cv_arr[1, jj], tk_arr[1, jj])
end

# Compute fact for top layer
dtime = 1800.0
col = inst.column
jj_top = nlevsno + 1
CAPR = 0.34
dz1 = col.dz[1, jj_top]
z1 = col.z[1, jj_top]
zi1 = col.zi[1, jj_top]
z2 = col.z[1, jj_top + 1]
zi2 = col.zi[1, jj_top + 1]
fact_top = dtime / cv_arr[1, jj_top] * dz1 / (0.5 * (z1 - zi1 + CAPR * (z2 - zi1)))
@printf("\nfact for top layer: %.6f\n", fact_top)
@printf("  dz1=%.6f z1=%.6f zi1=%.6f z2=%.6f\n", dz1, z1, zi1, z2)

# Expected ΔT from implicit solver
hs = ef.eflx_gnet_patch[1]
cgrd = ef.cgrnd_patch[1]
SB = 5.67e-8
emg = temp.emg_col[1]
dlwrad_emit = 4.0 * emg * SB * 274.0^3
dhsdT = -cgrd - dlwrad_emit
@printf("  hs=%.4f dhsdT=%.4f\n", hs, dhsdT)
@printf("  Expected dT = %.4f K\n", fact_top * hs / (1 + fact_top * abs(dhsdT)))

rm(fhistory, force=true)
