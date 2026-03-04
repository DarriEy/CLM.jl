using CLM, Dates, Printf

fsurdat = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/surfdata_clm.nc"
paramfile = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/clm5_params.nc"

# Just initialize, don't run
(inst, bounds, filt, tm) = CLM.clm_initialize!(fsurdat=fsurdat, paramfile=paramfile)

nlevsno = CLM.varpar.nlevsno
nlevsoi = CLM.varpar.nlevsoi
nlevgrnd = CLM.varpar.nlevgrnd
nlevmaxurbgrnd = CLM.varpar.nlevmaxurbgrnd

println("Dimensions:")
@printf("  nlevsno=%d nlevsoi=%d nlevgrnd=%d nlevmaxurbgrnd=%d\n",
    nlevsno, nlevsoi, nlevgrnd, nlevmaxurbgrnd)

wsb = inst.water.waterstatebulk_inst
ws = wsb.ws
col = inst.column
ss = inst.soilstate

println("\nh2osoi_liq/ice (col 1, all layers):")
for j in 1:(nlevsno + nlevmaxurbgrnd)
    liq = ws.h2osoi_liq_col[1, j]
    ice = ws.h2osoi_ice_col[1, j]
    dz = col.dz[1, j]
    label = j <= nlevsno ? "snow" : (j - nlevsno <= nlevsoi ? "soil" : "rock")
    @printf("  j=%2d (%4s): liq=%12.6f ice=%12.6f dz=%10.6f", j, label, liq, ice, dz)
    if j > nlevsno && (j - nlevsno) <= nlevsoi
        @printf("  watsat=%.4f", ss.watsat_col[1, j - nlevsno])
    end
    println()
end

println("\nt_soisno (col 1, first few layers):")
for j in 1:min(nlevsno + 5, size(inst.temperature.t_soisno_col, 2))
    @printf("  j=%2d: %.4f K\n", j, inst.temperature.t_soisno_col[1, j])
end

println("\nt_grnd_col:")
for c in 1:length(inst.temperature.t_grnd_col)
    @printf("  col %d: %.4f K\n", c, inst.temperature.t_grnd_col[c])
end

# Check frac_veg_nosno_alb after init
println("\nfrac_veg_nosno_alb_patch (after init):")
for p in 1:length(inst.canopystate.frac_veg_nosno_alb_patch)
    @printf("  patch %d: %d\n", p, inst.canopystate.frac_veg_nosno_alb_patch[p])
end
