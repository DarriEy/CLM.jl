using NCDatasets, Dates, Printf, Statistics, CLM

const run_clm! = getfield(CLM, Symbol("clm_run!"))
const USE_AQUIFER_LAYER = !any(==("--no-aquifer"), ARGS)

println("================================================================")
println("  Julia CLM vs Fortran CLM — Full Year Comparison")
println("  Site: Bow at Banff | Forcing: clmforc.2002.nc (real obs)")
println("  Features: seasonal phenology, SNICAR optics, PFT roughness")
println("  Hydrology mode: ", USE_AQUIFER_LAYER ? "aquifer" : "no-aquifer")
println("================================================================")
println()

# --- Input file paths ---
const basedir = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped"
fsurdat  = joinpath(basedir, "settings/CLM/parameters/surfdata_clm.nc")
paramfile = joinpath(basedir, "settings/CLM/parameters/clm5_params.nc")
fforcing = joinpath(basedir, "data/forcing/CLM_input/clmforc.2002.nc")
fhistory = tempname() * "_verify.nc"

# SNICAR data files
fsnowoptics = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata/snicar_optics_5bnd_c013122.nc"
fsnowaging  = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata/snicar_drdt_bst_fit_60_c070416.nc"

# Fortran reference
f_h0 = joinpath(basedir, "optimization/CLM/dds_run_1/final_evaluation/Bow_at_Banff_lumped.clm2.h0.2002-01-01-00000.nc")

# --- Run full year ---
ndays = 365
println("Running Julia CLM for $ndays days (full year 2002)...")
println("  SNICAR optics: $(isfile(fsnowoptics) ? "YES" : "NO")")
println("  SNICAR aging:  $(isfile(fsnowaging)  ? "YES" : "NO")")
println()

t0 = time()
inst = run_clm!(;
    fsurdat=fsurdat, paramfile=paramfile, fforcing=fforcing,
    fhistory=fhistory,
    start_date=DateTime(2002, 1, 1),
    end_date=DateTime(2003, 1, 1),
    dtime=1800, use_cn=false, verbose=true,
    use_aquifer_layer=USE_AQUIFER_LAYER,
    fsnowoptics=isfile(fsnowoptics) ? fsnowoptics : "",
    fsnowaging=isfile(fsnowaging)   ? fsnowaging  : "")
elapsed = time() - t0
@printf("  Done: %d timesteps in %.1f seconds (%.1f s/day)\n\n",
        ndays*48, elapsed, elapsed/ndays)

# ---- Helper functions ----
# All Julia output fields are now aggregated to gridcell level (lndgrid dimension),
# matching Fortran h0 output format.
function daily_mean_col(ds, varname, day; fill=-9999.0)
    i1 = (day-1)*48 + 1; i2 = min(day*48, size(ds[varname], 2))
    vals = Float64[]
    for t in i1:i2
        v = Float64(ds[varname][1, t])
        if v > fill+1 && isfinite(v); push!(vals, v); end
    end
    isempty(vals) ? NaN : mean(vals)
end

function daily_mean_patch(ds, varname, day; fill=-9999.0)
    # Now that all output is on gridcell level, this is the same as daily_mean_col
    return daily_mean_col(ds, varname, day; fill=fill)
end

function fortran_daily(fds, fname, day)
    if !haskey(fds, fname)
        return NaN
    end
    val = fds[fname]
    if ndims(val) >= 2
        v = val[1, day]
    else
        v = val[day]
    end
    return ismissing(v) ? NaN : Float64(v)
end

# ---- Read output ----
jds = NCDataset(fhistory)
fds = NCDataset(f_h0)
nt_j = size(jds["T_GRND"], 2)

# ================================================================
# Monthly comparison (more readable than 365 daily rows)
# ================================================================
vars = [
    ("T_GRND",        "TG",            "K",    "col"),
    ("H2OSNO",        "H2OSNO",        "mm",   "col"),
    ("FSA",           "FSA",           "W/m2", "patch"),
    ("EFLX_LH_TOT",  "EFLX_LH_TOT",  "W/m2", "patch"),
    ("EFLX_SH_TOT",  "FSH",           "W/m2", "patch"),
    ("TSA",           "TSA",           "K",    "patch"),
    ("QRUNOFF",       "QRUNOFF",       "mm/s", "col"),
]

# New variables from Phase 5
new_vars = [
    ("ELAI",       "ELAI",       "m2/m2", "patch"),
    ("SNOW_DEPTH", "SNOW_DEPTH", "m",     "col"),
    ("SNOWDP",     "SNOWDP",     "m",     "col"),
    ("ZWT",        "ZWT",        "m",     "col"),
    ("BTRAN",      "BTRANMN",    "0-1",   "patch"),
]

# Month boundaries (DOY)
month_starts = [1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335]
month_ends   = [31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365]
month_names  = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
                "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

println("================================================================")
println("  MONTHLY MEAN COMPARISON — Julia vs Fortran (2002)")
println("================================================================\n")

for (jname, fname, units, level) in vcat(vars, new_vars)
    has_fortran = haskey(fds, fname)
    @printf("%-16s [%-6s]", jname, units)
    if has_fortran
        @printf(" %10s %10s %10s\n", "Julia", "Fortran", "Diff")
    else
        @printf(" %10s %10s\n", "Julia", "(no ref)")
    end
    @printf("%-25s  %10s %10s %10s\n", "-"^25, "-"^10, "-"^10, "-"^10)

    annual_j = Float64[]
    annual_f = Float64[]

    for m in 1:12
        d1 = month_starts[m]; d2 = month_ends[m]

        # Julia monthly mean
        jvals = Float64[]
        for d in d1:min(d2, ndays)
            v = level == "col" ? daily_mean_col(jds, jname, d) :
                                 daily_mean_patch(jds, jname, d)
            isnan(v) || push!(jvals, v)
        end
        j_mean = isempty(jvals) ? NaN : mean(jvals)

        # Fortran monthly mean
        fvals = Float64[]
        if has_fortran
            for d in d1:d2
                v = fortran_daily(fds, fname, d)
                isnan(v) || push!(fvals, v)
            end
        end
        f_mean = isempty(fvals) ? NaN : mean(fvals)

        isnan(j_mean) || push!(annual_j, j_mean)
        isnan(f_mean) || push!(annual_f, f_mean)

        diff = (isnan(j_mean) || isnan(f_mean)) ? NaN : j_mean - f_mean
        if max(abs(j_mean), abs(f_mean)) > 100
            @printf("  %-3s  %10.2f", month_names[m], j_mean)
            has_fortran && @printf(" %10.2f %+10.2f", f_mean, diff)
            println()
        elseif max(abs(j_mean), abs(f_mean)) > 0.01
            @printf("  %-3s  %10.4f", month_names[m], j_mean)
            has_fortran && @printf(" %10.4f %+10.4f", f_mean, diff)
            println()
        else
            @printf("  %-3s  %10.6f", month_names[m], j_mean)
            has_fortran && @printf(" %10.6f %+10.6f", f_mean, diff)
            println()
        end
    end

    # Annual summary
    j_ann = isempty(annual_j) ? NaN : mean(annual_j)
    f_ann = isempty(annual_f) ? NaN : mean(annual_f)
    d_ann = (isnan(j_ann) || isnan(f_ann)) ? NaN : j_ann - f_ann
    @printf("  %-3s  ", "ANN")
    if max(abs(j_ann), abs(f_ann)) > 100
        @printf("%10.2f", j_ann)
        has_fortran && @printf(" %10.2f %+10.2f", f_ann, d_ann)
    elseif max(abs(j_ann), abs(f_ann)) > 0.01
        @printf("%10.4f", j_ann)
        has_fortran && @printf(" %10.4f %+10.4f", f_ann, d_ann)
    else
        @printf("%10.6f", j_ann)
        has_fortran && @printf(" %10.6f %+10.6f", f_ann, d_ann)
    end
    println()
    println()
end

# ================================================================
# Seasonal phenology check (LAI should vary)
# ================================================================
println("================================================================")
println("  PHENOLOGY CHECK — LAI seasonal cycle")
println("================================================================")
if haskey(jds, "ELAI")
    for m in 1:12
        d1 = month_starts[m]; d2 = month_ends[m]
        vals = Float64[]
        for d in d1:min(d2, ndays)
            v = daily_mean_patch(jds, "ELAI", d)
            isnan(v) || push!(vals, v)
        end
        j_mean = isempty(vals) ? NaN : mean(vals)
        @printf("  %-3s ELAI: %6.3f", month_names[m], j_mean)
        if haskey(fds, "ELAI")
            fvals = Float64[]
            for d in d1:d2
                v = fortran_daily(fds, "ELAI", d)
                isnan(v) || push!(fvals, v)
            end
            f_mean = isempty(fvals) ? NaN : mean(fvals)
            @printf("  (Fortran: %6.3f, diff: %+6.3f)", f_mean, j_mean - f_mean)
        end
        println()
    end
elseif haskey(jds, "TLAI")
    for m in 1:12
        d1 = month_starts[m]; d2 = month_ends[m]
        vals = Float64[]
        for d in d1:min(d2, ndays)
            v = daily_mean_patch(jds, "TLAI", d)
            isnan(v) || push!(vals, v)
        end
        j_mean = isempty(vals) ? NaN : mean(vals)
        @printf("  %-3s TLAI: %6.3f\n", month_names[m], j_mean)
    end
end

# ================================================================
# T_GRND evolution check
# ================================================================
println()
println("================================================================")
println("  T_GRND EVOLUTION — seasonal cycle verification")
println("================================================================")
tg_all = Float64[]
for t in 1:nt_j
    v = Float64(jds["T_GRND"][1, t])
    if v > -9000 && isfinite(v)
        push!(tg_all, v)
    end
end

@printf("  Winter min (Jan):    %.2f K\n", minimum(tg_all[1:min(48*31, end)]))
@printf("  Summer max (Jul):    %.2f K\n",
        length(tg_all) > 48*181 ? maximum(tg_all[48*181:min(48*212, end)]) : NaN)
@printf("  Annual range:        [%.2f, %.2f] K\n", minimum(tg_all), maximum(tg_all))
@printf("  Annual std dev:      %.2f K\n", std(tg_all))

evolves = std(tg_all) > 1.0
println("  Seasonal cycle:      ", evolves ? "YES" : "NO (problem!)")

if haskey(fds, "TG")
    # Compare seasonal amplitude
    f_tg_jan = mean([fortran_daily(fds, "TG", d) for d in 1:31])
    f_tg_jul = mean([fortran_daily(fds, "TG", d) for d in 182:212])
    j_tg_jan = mean([daily_mean_col(jds, "T_GRND", d) for d in 1:31])
    j_tg_jul = mean([daily_mean_col(jds, "T_GRND", d) for d in 182:min(212, ndays)])
    @printf("\n  Seasonal amplitude (Jul-Jan):\n")
    @printf("    Julia:   %+.2f K\n", j_tg_jul - j_tg_jan)
    @printf("    Fortran: %+.2f K\n", f_tg_jul - f_tg_jan)
end

close(jds)
close(fds)

rm(fhistory, force=true)

println()
println("================================================================")
println("  VERIFICATION COMPLETE")
println("================================================================")
