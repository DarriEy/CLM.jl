# ==========================================================================
# Multi-year spinup script for CLM.jl
#
# Runs CLM repeatedly using restart files to equilibrate state variables
# (ZWT, soil moisture, soil temperature, snow) from cold start.
#
# Usage:
#   julia --project=. test/spinup.jl [num_cycles]
#
# Default: 5 cycles of 2002-2009 (8 years each = 40 years total)
# ==========================================================================

using Dates, Printf, CLM

# ---- Configuration ----
basedir = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped"
fsurdat = joinpath(basedir, "settings/CLM/parameters/surfdata_clm.nc")
paramfile = joinpath(basedir, "settings/CLM/parameters/clm5_params.nc")
fsnowoptics = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata/snicar_optics_5bnd_c013122.nc"
fsnowaging = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata/snicar_drdt_bst_fit_60_c070416.nc"

# Forcing files available: 2002-2009
forcing_dir = joinpath(basedir, "data/forcing/CLM_input")
forcing_years = 2002:2009

# Spinup parameters
num_cycles = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 5
output_dir = "/tmp/clm_spinup"
mkpath(output_dir)

# Restart file paths: clm_run! reads from `frestart` and writes to `frestart_out.nc`
restart_base = joinpath(output_dir, "restart.nc")
restart_out = joinpath(output_dir, "restart_out.nc")

println("=" ^ 70)
println("CLM.jl SPINUP")
println("  Cycles: $num_cycles × $(length(forcing_years)) years = $(num_cycles * length(forcing_years)) total years")
println("  Forcing: $forcing_dir ($(first(forcing_years))-$(last(forcing_years)))")
println("  Output: $output_dir")
println("  SNICAR optics: $(isfile(fsnowoptics) ? "loaded" : "NOT FOUND")")
println("=" ^ 70)

total_start = time()

for cycle in 1:num_cycles
    println("\n=== CYCLE $cycle / $num_cycles ===")
    cycle_start = time()

    for (yi, year) in enumerate(forcing_years)
        fforcing = joinpath(forcing_dir, "clmforc.$year.nc")
        if !isfile(fforcing)
            @warn "Forcing file not found, skipping: $fforcing"
            continue
        end

        start_date = DateTime(year, 1, 1)
        end_date = DateTime(year + 1, 1, 1)
        fhistory = joinpath(output_dir, "clm_hist_c$(lpad(cycle,2,'0'))_$year.nc")

        year_start = time()

        inst = CLM.clm_run!(;
            fsurdat=fsurdat,
            paramfile=paramfile,
            fforcing=fforcing,
            fhistory=fhistory,
            start_date=start_date,
            end_date=end_date,
            dtime=1800,
            fsnowoptics=fsnowoptics,
            fsnowaging=fsnowaging,
            frestart=restart_base,
            verbose=false)

        # clm_run! writes to restart_out.nc — rename for next iteration
        if isfile(restart_out)
            cp(restart_out, restart_base; force=true)
            rm(restart_out)
        end

        # Diagnostic: key state variables
        nlevsno = CLM.varpar.nlevsno
        t_grnd = inst.temperature.t_grnd_col[1]
        zwt = inst.soilhydrology.zwt_col[1]
        snow_depth = inst.water.waterdiagnosticbulk_inst.snow_depth_col[1]

        elapsed = time() - year_start
        sim_year = (cycle - 1) * length(forcing_years) + yi
        @printf("  [Year %2d] %d: T_GRND=%.2f K, ZWT=%.3f m, snow=%.3f m (%.0f s)\n",
                sim_year, year, t_grnd, zwt, snow_depth, elapsed)
    end

    cycle_elapsed = time() - cycle_start
    @printf("  Cycle %d complete in %.0f s (%.1f min)\n", cycle, cycle_elapsed, cycle_elapsed/60)
end

total_elapsed = time() - total_start
println("\n" * "=" ^ 70)
@printf("SPINUP COMPLETE: %d years in %.0f s (%.1f min)\n",
        num_cycles * length(forcing_years), total_elapsed, total_elapsed/60)
println("Final restart: $restart_base")
println("=" ^ 70)
