# ==========================================================================
# Spinup convergence test
# Runs CLM.jl from cold start for multiple years, cycling forcing 2002-2009,
# and tracks key state variables to verify equilibration.
# ==========================================================================

using NCDatasets, Dates, Printf, Statistics, CLM

const run_clm! = getfield(CLM, Symbol("clm_run!"))

const basedir = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped"
const caldir = joinpath(basedir, "optimization/CLM/dds_run_1/final_evaluation/settings/CLM/parameters")

const fsurdat  = joinpath(caldir, "surfdata_clm.nc")
const paramfile = joinpath(caldir, "clm5_params.nc")
const fsnowoptics = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata/snicar_optics_5bnd_c013122.nc"
const fsnowaging  = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata/snicar_drdt_bst_fit_60_c070416.nc"

const FORCING_YEARS = 2002:2009
const N_SPINUP_YEARS = parse(Int, get(ARGS, 1, "8"))  # default 8 years (one full cycle)

println("================================================================")
println("  CLM.jl Spinup Convergence Test")
println("  Site: Bow at Banff | Cold start")
println("  Forcing cycle: $(first(FORCING_YEARS))-$(last(FORCING_YEARS))")
println("  Spinup years: $N_SPINUP_YEARS")
println("================================================================")
println()

# Track annual means per year
annual_tgrnd = Float64[]
annual_tsa   = Float64[]
annual_sh    = Float64[]
annual_lh    = Float64[]
annual_fsa   = Float64[]
annual_qrun  = Float64[]
annual_h2o1  = Float64[]  # top-layer soil moisture

restart_file = ""
t0_total = time()

for yr_idx in 1:N_SPINUP_YEARS
    # Cycle through forcing years
    forc_year = FORCING_YEARS[mod1(yr_idx, length(FORCING_YEARS))]
    fforcing = joinpath(basedir, "data/forcing/CLM_input/clmforc.$forc_year.nc")
    fhistory = tempname() * "_spinup_yr$(yr_idx).nc"

    # Build kwargs
    run_kwargs = Dict{Symbol,Any}(
        :fsurdat => fsurdat, :paramfile => paramfile,
        :fforcing => fforcing, :fhistory => fhistory,
        :start_date => DateTime(forc_year, 1, 1),
        :end_date => DateTime(forc_year + 1, 1, 1),
        :dtime => 1800, :use_cn => false,
        :verbose => false,
        :use_aquifer_layer => false,
        :baseflow_scalar => 0.0022119554,
        :int_snow_max => 3113.2227,
        :fsnowoptics => isfile(fsnowoptics) ? fsnowoptics : "",
        :fsnowaging => isfile(fsnowaging) ? fsnowaging : "",
    )

    # frestart: if file exists at path → read it; always writes to frestart_out.nc
    if yr_idx > 1 && isfile(restart_file)
        run_kwargs[:frestart] = restart_file  # exists → will be read; output to restart_file_out.nc
    else
        run_kwargs[:frestart] = tempname() * "_spinup_cold"  # doesn't exist → cold start
    end

    t0 = time()
    inst = run_clm!(; run_kwargs...)
    elapsed = time() - t0

    # Read history and compute annual means
    NCDataset(fhistory) do ds
        ndays = size(ds["T_GRND"], 2)
        tg = [Float64(ds["T_GRND"][1, d]) for d in 1:ndays]
        push!(annual_tgrnd, mean(filter(isfinite, tg)))

        if haskey(ds, "TSA")
            ta = [Float64(ds["TSA"][1, d]) for d in 1:ndays]
            push!(annual_tsa, mean(filter(isfinite, ta)))
        else
            push!(annual_tsa, NaN)
        end

        if haskey(ds, "EFLX_SH_TOT")
            sh = [Float64(ds["EFLX_SH_TOT"][1, d]) for d in 1:ndays]
            push!(annual_sh, mean(filter(isfinite, sh)))
        else
            push!(annual_sh, NaN)
        end

        if haskey(ds, "EFLX_LH_TOT")
            lh = [Float64(ds["EFLX_LH_TOT"][1, d]) for d in 1:ndays]
            push!(annual_lh, mean(filter(isfinite, lh)))
        else
            push!(annual_lh, NaN)
        end

        if haskey(ds, "FSA")
            fsa = [Float64(ds["FSA"][1, d]) for d in 1:ndays]
            push!(annual_fsa, mean(filter(isfinite, fsa)))
        else
            push!(annual_fsa, NaN)
        end

        if haskey(ds, "QRUNOFF")
            qr = [Float64(ds["QRUNOFF"][1, d]) for d in 1:ndays]
            push!(annual_qrun, mean(filter(isfinite, qr)))
        else
            push!(annual_qrun, NaN)
        end

        if haskey(ds, "H2OSOI_VOL1")
            h2o = [Float64(ds["H2OSOI_VOL1"][1, d]) for d in 1:ndays]
            push!(annual_h2o1, mean(filter(isfinite, h2o)))
        else
            push!(annual_h2o1, NaN)
        end
    end

    # Update restart file for next year (clm_run! writes to frestart path with _out.nc suffix)
    frestart_val = run_kwargs[:frestart]
    global restart_file = replace(frestart_val, r"\.nc$" => "") * "_out.nc"
    if !endswith(frestart_val, ".nc")
        restart_file = frestart_val * "_out.nc"
    end

    # Clean up history
    rm(fhistory, force=true)

    # Print progress
    δtg = yr_idx > 1 ? annual_tgrnd[end] - annual_tgrnd[end-1] : 0.0
    δh2o = yr_idx > 1 ? annual_h2o1[end] - annual_h2o1[end-1] : 0.0
    @printf("  Year %2d (forc=%d): T_GRND=%.2fK (Δ=%+.3f) H2O_L1=%.4f (Δ=%+.5f) SH=%.1f LH=%.1f  [%.0fs]\n",
            yr_idx, forc_year, annual_tgrnd[end], δtg,
            annual_h2o1[end], δh2o,
            annual_sh[end], annual_lh[end], elapsed)
end

total_elapsed = time() - t0_total

println()
println("================================================================")
println("  SPINUP CONVERGENCE SUMMARY")
println("================================================================")
println()
@printf("  %-12s  %-10s  %-10s  %-10s  %-10s  %-10s\n",
        "Variable", "Year 1", "Year $(N_SPINUP_YEARS)", "Δ(last)", "Fortran*", "Gap")
println("  " * "-"^70)

# Fortran reference values (from restart-mode verification, 2003)
fort_tgrnd = 272.74
fort_tsa   = 270.54
fort_sh    = 52.29
fort_lh    = 12.33

δ_tg = N_SPINUP_YEARS > 1 ? annual_tgrnd[end] - annual_tgrnd[end-1] : 0.0
@printf("  %-12s  %8.2fK   %8.2fK   %+7.3fK   %8.2fK   %+.2fK\n",
        "T_GRND", annual_tgrnd[1], annual_tgrnd[end], δ_tg, fort_tgrnd,
        annual_tgrnd[end] - fort_tgrnd)

δ_ta = N_SPINUP_YEARS > 1 ? annual_tsa[end] - annual_tsa[end-1] : 0.0
@printf("  %-12s  %8.2fK   %8.2fK   %+7.3fK   %8.2fK   %+.2fK\n",
        "TSA", annual_tsa[1], annual_tsa[end], δ_ta, fort_tsa,
        annual_tsa[end] - fort_tsa)

δ_sh = N_SPINUP_YEARS > 1 ? annual_sh[end] - annual_sh[end-1] : 0.0
@printf("  %-12s  %7.1fW    %7.1fW    %+6.2fW    %7.1fW    %+.1fW\n",
        "EFLX_SH", annual_sh[1], annual_sh[end], δ_sh, fort_sh,
        annual_sh[end] - fort_sh)

δ_lh = N_SPINUP_YEARS > 1 ? annual_lh[end] - annual_lh[end-1] : 0.0
@printf("  %-12s  %7.1fW    %7.1fW    %+6.2fW    %7.1fW    %+.1fW\n",
        "EFLX_LH", annual_lh[1], annual_lh[end], δ_lh, fort_lh,
        annual_lh[end] - fort_lh)

println()
println("  *Fortran reference: spun-up restart, year 2003 (calibrated params)")
@printf("  Total time: %.0f seconds (%.1f s/year)\n", total_elapsed, total_elapsed/N_SPINUP_YEARS)

# Convergence check
if N_SPINUP_YEARS >= 3
    # Check if last 3 years are stable (< 0.1K T_GRND drift)
    recent_drift = maximum(abs.(diff(annual_tgrnd[end-2:end])))
    if recent_drift < 0.1
        println("  ✓ T_GRND converged (drift < 0.1K over last 3 years)")
    else
        @printf("  ⚠ T_GRND still drifting (%.3fK between last years)\n", recent_drift)
    end
end

println()
