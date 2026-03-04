using CLM, Printf, Dates, LinearAlgebra

const clm_init! = getfield(CLM, Symbol("clm_initialize!"))
const clm_drv! = getfield(CLM, Symbol("clm_drv!"))

println("=== T_GRND Diagnostic ===\n")

fsurdat = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/surfdata_clm.nc"
paramfile = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/clm5_params.nc"

inst, bounds, filt, tm = clm_init!(; fsurdat=fsurdat, paramfile=paramfile,
                                     start_date=DateTime(2002,1,1), dtime=1800, use_cn=false)

nc = length(inst.column.itype)
println("Number of columns: $nc")
println("Column types: ", inst.column.itype[1:nc])

# Check filters
println("\n--- FILTERS ---")
println("nolakec active: ", count(filt.nolakec), " / ", length(filt.nolakec))
println("  nolakec[1] = ", filt.nolakec[1], "  nolakec[2] = ", filt.nolakec[2])
println("soilc active:   ", count(filt.soilc), " / ", length(filt.soilc))
println("urbanc active:  ", count(filt.urbanc), " / ", length(filt.urbanc))
println("lakec active:   ", count(filt.lakec), " / ", length(filt.lakec))

# Check bounds
println("\n--- BOUNDS ---")
@printf("  begc=%d  endc=%d  begg=%d  endg=%d  begp=%d  endp=%d\n",
        bounds.begc, bounds.endc, bounds.begg, bounds.endg, bounds.begp, bounds.endp)

# Check landunit type
println("\n--- LANDUNIT TYPES ---")
for c in 1:nc
    l = inst.column.landunit[c]
    @printf("  col %d: itype=%d  landunit=%d  lun.itype=%d\n",
            c, inst.column.itype[c], l, inst.landunit.itype[l])
end

# Check soil thermal properties
println("\n--- SOIL THERMAL PROPERTIES (column 1) ---")
ss = inst.soilstate
for j in 1:min(5, size(ss.tkmg_col, 2))
    @printf("  Layer %2d: tkmg=%.4f  tkdry=%.4f  csol=%.1f  watsat=%.4f\n",
            j, ss.tkmg_col[1,j], ss.tkdry_col[1,j], ss.csol_col[1,j], ss.watsat_col[1,j])
end

# Check layer thicknesses
println("\n--- LAYER GEOMETRY (column 1) ---")
col = inst.column
joff = size(col.dz, 2) > 20 ? (size(col.dz, 2) - 20) : 0
# Check if there's a snow offset
nlevsno = CLM.varpar.nlevsno  # actual max snow layers
for j in 1:min(5, size(col.dz, 2))
    @printf("  Raw idx %2d: dz=%.4f  z=%.4f  zi=%.6f\n",
            j, col.dz[1,j], col.z[1,j], col.zi[1,j])
end
println("  ...")
# Show soil layers (after snow offset)
for j in (nlevsno+1):min(nlevsno+5, size(col.dz, 2))
    @printf("  Raw idx %2d (soil %2d): dz=%.4f  z=%.4f\n",
            j, j-nlevsno, col.dz[1,j], col.z[1,j])
end

# Check initial temperatures
println("\n--- INITIAL TEMPERATURES ---")
temp = inst.temperature
@printf("  t_grnd_col[1] = %.4f K\n", temp.t_grnd_col[1])
for j in 1:min(5, size(temp.t_soisno_col, 2))
    @printf("  t_soisno_col[1,%2d] = %.4f K\n", j, temp.t_soisno_col[1,j])
end
println("  ...")
for j in (nlevsno+1):min(nlevsno+5, size(temp.t_soisno_col, 2))
    @printf("  t_soisno_col[1,%2d] (soil %2d) = %.4f K\n",
            j, j-nlevsno, temp.t_soisno_col[1,j])
end

# Check water content
println("\n--- INITIAL WATER CONTENT (column 1) ---")
ws = inst.water.waterstatebulk_inst.ws
for j in (nlevsno+1):min(nlevsno+5, size(ws.h2osoi_liq_col, 2))
    @printf("  Soil %2d: liq=%.4f  ice=%.4f kg/m2\n",
            j-nlevsno, ws.h2osoi_liq_col[1,j], ws.h2osoi_ice_col[1,j])
end

# Check snl (snow layers)
println("\n--- SNOW STATE ---")
@printf("  snl[1] = %d\n", col.snl[1])
@printf("  h2osno = %.4f mm\n", ws.h2osno_no_layers_col[1])

# Check energy fluxes (patch-level)
println("\n--- INITIAL ENERGY FLUXES ---")
ef = inst.energyflux
@printf("  eflx_sh_grnd_patch[1] = %.4f W/m2\n", ef.eflx_sh_grnd_patch[1])
@printf("  eflx_soil_grnd_patch[1] = %.4f W/m2\n", ef.eflx_soil_grnd_patch[1])
@printf("  cgrnds_patch[1] = %.6f\n", ef.cgrnds_patch[1])
@printf("  cgrnd_patch[1] = %.6f\n", ef.cgrnd_patch[1])
@printf("  eflx_bot_col[1] = %.6f\n", ef.eflx_bot_col[1])

# Now set up forcing and run ONE timestep
println("\n=== RUNNING ONE TIMESTEP ===")

# Load forcing
fforcing = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/data/forcing/CLM_input/clmforc.2002.nc"
fr = CLM.ForcingReader()
CLM.forcing_reader_init!(fr, fforcing)
ng = 1  # gridcells
nc_f = nc  # columns

# Read first forcing step
CLM.read_forcing_step!(fr, inst.atm2lnd, DateTime(2002,1,1,0,30), ng, nc_f)

println("\n--- FORCING VALUES ---")
a2l = inst.atm2lnd
@printf("  forc_t  = %.2f K\n", a2l.forc_t_downscaled_col[1])
@printf("  forc_lwrad = %.2f W/m2\n", a2l.forc_lwrad_downscaled_col[1])
@printf("  forc_rain = %.6f mm/s\n", a2l.forc_rain_downscaled_col[1])
@printf("  forc_snow = %.6f mm/s\n", a2l.forc_snow_downscaled_col[1])

# Downscale
CLM.downscale_forcings!(bounds, a2l, inst.column, inst.landunit, inst.topo)

@printf("  forc_t (downscaled) = %.2f K\n", a2l.forc_t_downscaled_col[1])
@printf("  forc_pbot = %.2f Pa\n", a2l.forc_pbot_downscaled_col[1])

# Save T_GRND before
tgrnd_before = copy(temp.t_grnd_col)
tsoisno_before = copy(temp.t_soisno_col[1, :])

# Run driver (matching clm_run! exactly)
println("\nRunning clm_drv!...")
try
    config = CLM.CLMDriverConfig(use_cn=false)
    filt_ia = CLM.clump_filter_inactive_and_active

    calday = 1.0 + 1800.0 / 86400.0
    (declin, eccf) = CLM.compute_orbital(calday)
    declinp1 = declin
    obliqr = CLM.ORB_OBLIQR_DEFAULT
    nextsw_cday = calday + 1800.0 / 86400.0

    clm_drv!(config, inst, filt, filt_ia, bounds,
             true, nextsw_cday, declinp1, declin, obliqr,
             false, false, "", false;
             nstep=1,
             is_first_step=true,
             is_beg_curr_day=true,
             is_beg_curr_year=true,
             dtime=1800.0,
             photosyns=inst.photosyns)
    println("  Driver completed successfully")
catch e
    println("  Driver error: ", e)
    for (exc, bt) in Base.catch_stack()
        showerror(stdout, exc, bt)
        println()
    end
end

# Check T_GRND after
println("\n--- TEMPERATURES AFTER ONE TIMESTEP ---")
@printf("  t_grnd_col[1]: %.6f → %.6f  (diff = %+.6f)\n",
        tgrnd_before[1], temp.t_grnd_col[1], temp.t_grnd_col[1] - tgrnd_before[1])
for j in (nlevsno+1):min(nlevsno+5, size(temp.t_soisno_col, 2))
    @printf("  t_soisno[1,%2d]: %.6f → %.6f  (diff = %+.6f)\n",
            j, tsoisno_before[j], temp.t_soisno_col[1,j],
            temp.t_soisno_col[1,j] - tsoisno_before[j])
end

# Check energy fluxes after driver
println("\n--- ENERGY FLUXES AFTER DRIVER ---")
@printf("  eflx_sh_grnd_patch[1] = %.4f W/m2\n", ef.eflx_sh_grnd_patch[1])
@printf("  eflx_soil_grnd_patch[1] = %.4f W/m2\n", ef.eflx_soil_grnd_patch[1])
@printf("  sabg_patch[1] = %.4f W/m2\n", inst.solarabs.sabg_patch[1])
@printf("  fsa_patch[1]  = %.4f W/m2\n", inst.solarabs.fsa_patch[1])
@printf("  cgrnds_patch[1] = %.6f\n", ef.cgrnds_patch[1])
@printf("  cgrnd_patch[1]  = %.6f\n", ef.cgrnd_patch[1])
@printf("  eflx_bot_col[1] = %.6f\n", ef.eflx_bot_col[1])

# Check fact (dtime/cv) which scales the band matrix
println("\n--- FACT (dtime/cv) for first 5 soil layers ---")
for j in (nlevsno+1):min(nlevsno+5, size(temp.fact_col, 2))
    @printf("  fact[1,%2d] = %.6e\n", j, temp.fact_col[1,j])
end

CLM.forcing_reader_close!(fr)
println("\n=== DIAGNOSTIC COMPLETE ===")
