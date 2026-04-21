# Feature parity test: runs CLM with different feature flags and compares annual means to Fortran.
# Usage: julia --project=. test/verify_features.jl [--luna] [--phs] [--all]
using NCDatasets, Dates, Printf, Statistics, CLM

const run_clm! = getfield(CLM, Symbol("clm_run!"))

# Parse flags
const DO_LUNA = any(==("--luna"), ARGS) || any(==("--all"), ARGS)
const DO_PHS  = any(==("--phs"), ARGS) || any(==("--all"), ARGS)
const DO_COMBINED = any(==("--all"), ARGS)

# --- Input file paths ---
const basedir = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped"
const caldir  = joinpath(basedir, "optimization/CLM/dds_run_1/final_evaluation/settings/CLM/parameters")
const fsurdat  = joinpath(caldir, "surfdata_clm.nc")
const paramfile = joinpath(caldir, "clm5_params.nc")
const fforcing = joinpath(basedir, "data/forcing/CLM_input/clmforc.2002.nc")
const fsnowoptics = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata/snicar_optics_5bnd_c013122.nc"
const fsnowaging  = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata/snicar_drdt_bst_fit_60_c070416.nc"
const f_h0 = joinpath(basedir, "optimization/CLM/dds_run_1/final_evaluation/Bow_at_Banff_lumped.clm2.h0.2002-01-01-00000.nc")

# Fortran reference annual means
const fortran_ref = Dict{String,Float64}()
let fds = NCDataset(f_h0)
    for (jname, fname) in [
        ("T_GRND","TG"), ("TSA","TSA"), ("FSA","FSA"),
        ("EFLX_LH_TOT","EFLX_LH_TOT"), ("EFLX_SH_TOT","FSH"),
        ("FCEV","FCEV"), ("FCTR","FCTR"), ("FGEV","FGEV"),
        ("BTRAN","BTRANMN"), ("QRUNOFF","QRUNOFF"), ("SOILRESIS","SOILRESIS")]
        haskey(fds, fname) || continue
        val = fds[fname]
        vals = Float64[]
        for d in 1:365
            v = ndims(val) >= 2 ? val[1,d] : val[d]
            ismissing(v) || push!(vals, Float64(v))
        end
        fortran_ref[jname] = mean(vals)
    end
    close(fds)
end

# Helper: run one config and extract annual means
function run_config(label::String; kwargs...)
    fhistory = tempname() * "_$(label).nc"
    common = Dict{Symbol,Any}(
        :fsurdat => fsurdat, :paramfile => paramfile, :fforcing => fforcing,
        :fhistory => fhistory,
        :start_date => DateTime(2002,1,1), :end_date => DateTime(2003,1,1),
        :dtime => 1800, :use_cn => false, :verbose => false,
        :use_aquifer_layer => false,
        :baseflow_scalar => 0.0022119554, :int_snow_max => 3113.2227,
        :fsnowoptics => isfile(fsnowoptics) ? fsnowoptics : "",
        :fsnowaging  => isfile(fsnowaging) ? fsnowaging : "")
    merge!(common, Dict(kwargs))

    println("\n================================================================")
    println("  Running: $label")
    println("  Features: ", join(["$(k)=$(v)" for (k,v) in kwargs], ", "))
    println("================================================================")

    t0 = time()
    run_clm!(; common...)
    elapsed = time() - t0
    @printf("  Completed in %.1f seconds\n\n", elapsed)

    # Extract annual means
    results = Dict{String,Float64}()
    jds = NCDataset(fhistory)
    nt = size(jds["T_GRND"], 2)
    for vname in ["T_GRND","TSA","FSA","EFLX_LH_TOT","EFLX_SH_TOT",
                   "FCEV","FCTR","FGEV","BTRAN","QRUNOFF","SOILRESIS"]
        haskey(jds, vname) || continue
        vals = Float64[]
        for d in 1:min(365, nt)
            v = Float64(jds[vname][1, d])
            (isfinite(v) && v > -9998.0) || continue
            push!(vals, v)
        end
        results[vname] = isempty(vals) ? NaN : mean(vals)
    end
    close(jds)
    rm(fhistory, force=true)
    return results
end

# Print comparison table
function print_table(configs::Vector{Pair{String, Dict{String,Float64}}})
    varnames = ["T_GRND","TSA","FSA","EFLX_LH_TOT","EFLX_SH_TOT",
                "FCEV","FCTR","FGEV","BTRAN","QRUNOFF","SOILRESIS"]
    units = ["K","K","W/m2","W/m2","W/m2","W/m2","W/m2","W/m2","0-1","mm/s","s/m"]

    println("\n================================================================")
    println("  ANNUAL MEAN COMPARISON — ALL CONFIGURATIONS")
    println("================================================================")

    # Header
    @printf("%-16s %-6s %10s", "Variable", "Units", "Fortran")
    for (label, _) in configs
        @printf(" %12s", label)
    end
    println()
    println("-"^(36 + 13*length(configs)))

    for (i, vname) in enumerate(varnames)
        f_val = get(fortran_ref, vname, NaN)
        @printf("%-16s %-6s", vname, units[i])
        if vname == "QRUNOFF"
            @printf(" %10.3f", f_val * 86400)
        elseif f_val > 100
            @printf(" %10.2f", f_val)
        elseif f_val > 0.01
            @printf(" %10.4f", f_val)
        else
            @printf(" %10.6f", f_val)
        end

        for (label, res) in configs
            j_val = get(res, vname, NaN)
            if vname == "QRUNOFF"
                j_val_disp = j_val * 86400
                diff_pct = isnan(f_val) || f_val == 0.0 ? NaN : (j_val - f_val)/abs(f_val)*100
                @printf(" %7.3f(%+.1f%%)", j_val_disp, diff_pct)
            elseif abs(f_val) > 0.01
                diff_pct = isnan(f_val) || f_val == 0.0 ? NaN : (j_val - f_val)/abs(f_val)*100
                @printf(" %7.2f(%+.1f%%)", j_val, diff_pct)
            else
                @printf(" %12.6f", j_val)
            end
        end
        println()
    end
end

# ================================================================
# Run configurations
# ================================================================
configs = Pair{String, Dict{String,Float64}}[]

if DO_LUNA
    res = run_config("luna"; use_luna=true, use_hydrstress=false)
    push!(configs, "LUNA" => res)
end

if DO_PHS
    res = run_config("phs"; use_luna=false, use_hydrstress=true)
    push!(configs, "PHS" => res)
end

if DO_COMBINED
    res = run_config("luna+phs"; use_luna=true, use_hydrstress=true)
    push!(configs, "LUNA+PHS" => res)
end

if !isempty(configs)
    print_table(configs)
end

println("\n================================================================")
println("  FEATURE PARITY TEST COMPLETE")
println("================================================================")
