# LOCALIZE the exact-vs-smoothed divergence: spin up N warmup steps in EXACT mode,
# snapshot the whole instance, then take ONE more step in exact mode and ONE more in
# smoothed mode FROM THE SAME STATE, and rank every array field by relative divergence.
# The top of that ranking names the module whose smoothing is biased.
#
#   julia +1.12 --project=. diverge.jl --warmup 4800
using CLM, NCDatasets, Dates, Printf, Statistics

_argval(flag, default) = begin
    i = findfirst(==(flag), ARGS)
    (i === nothing || i == length(ARGS)) ? default : ARGS[i + 1]
end
const WARMUP = parse(Int, _argval("--warmup", "4800"))
const DTIME  = 1800

const basedir = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped"
const caldir  = joinpath(basedir, "optimization/CLM/dds_run_1/final_evaluation/settings/CLM/parameters")
const fsurdat   = joinpath(caldir, "surfdata_clm.nc")
const paramfile = joinpath(caldir, "clm5_params.nc")
const fforcing  = joinpath(basedir, "data/forcing/CLM_input/clmforc.2002_2004.nc")
const fsnowoptics = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata/snicar_optics_5bnd_c013122.nc"
const fsnowaging  = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata/snicar_drdt_bst_fit_60_c070416.nc"
const ffortran_restart = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_parity_run/Bow_at_Banff_lumped.clm2.r.2003-01-01-00000.nc"
const baseflow_scalar = 0.0022119554
const int_snow_max    = 3113.2227
const USE_AQUIFER     = false
const start_date = DateTime(2003, 1, 1)

function build()
    (inst, bounds, filt, tm) = CLM.clm_initialize!(;
        fsurdat=fsurdat, paramfile=paramfile,
        start_date=start_date, dtime=DTIME, use_cn=false,
        use_bedrock=true, use_aquifer_layer=USE_AQUIFER, h2osfcflag=0,
        fsnowoptics=fsnowoptics, fsnowaging=fsnowaging, int_snow_max=int_snow_max)
    config = CLM.CLMDriverConfig(use_cn=false, use_aquifer_layer=USE_AQUIFER)
    CLM.atm2lnd_read_namelist!(inst.atm2lnd;
        repartition_rain_snow=true, lapse_rate=0.006, lapse_rate_longwave=0.032,
        precip_repartition_nonglc_all_snow_t=0.0, precip_repartition_nonglc_all_rain_t=2.0,
        precip_repartition_glc_all_snow_t=-2.0, precip_repartition_glc_all_rain_t=0.0)
    CLM.init_soil_hydrology_config(baseflow_scalar=baseflow_scalar)
    ng, nc, np = bounds.endg, bounds.endc, bounds.endp
    let ds_p = NCDataset(paramfile, "r")
        scf = inst.scf_method
        if haskey(ds_p, "n_melt_coef")
            nmc = Float64(ds_p["n_melt_coef"][1])
            for c in 1:nc; scf.n_melt[c] = nmc / max(10.0, inst.column.topo_std[c]); end
        end
        haskey(ds_p, "accum_factor") && (scf.accum_factor = Float64(ds_p["accum_factor"][1]))
        haskey(ds_p, "SNOW_DENSITY_MAX") && (CLM.snowhydrology_params.rho_max = Float64(ds_p["SNOW_DENSITY_MAX"][1]))
        haskey(ds_p, "SNOW_DENSITY_MIN") && (CLM.snowhydrology_params.rho_min = Float64(ds_p["SNOW_DENSITY_MIN"][1]))
        haskey(ds_p, "fresh_snw_rds_max") && (CLM.snowhydrology_params.snw_rds_min = Float64(ds_p["fresh_snw_rds_max"][1]))
        haskey(ds_p, "SNO_Z0MV") && (inst.frictionvel.zsno = Float64(ds_p["SNO_Z0MV"][1]))
        if haskey(ds_p, "snw_aging_bst"); CLM.snicar_params.xdrdt = Float64(ds_p["snw_aging_bst"][1]); end
        if haskey(ds_p, "pc")
            pc_val = Float64(ds_p["pc"][1])
            pc_val > 0 && (for c in 1:nc; inst.soilhydrology.hkdepth_col[c] = 1.0 / pc_val; end)
        end
        close(ds_p)
    end
    CLM.read_fortran_restart!(ffortran_restart, inst, bounds)
    fr = CLM.ForcingReader(); CLM.forcing_reader_init!(fr, fforcing)
    let topo_file = replace(fforcing, r"clmforc\.[^/]*\.nc$" => "topo_forcing.nc")
        if isfile(topo_file)
            dt = NCDataset(topo_file, "r")
            if haskey(dt, "TOPO")
                ft = Float64(dt["TOPO"][1])
                for g in 1:ng; inst.atm2lnd.forc_topo_grc[g] = ft; end
                for c in 1:nc; inst.topo.topo_col[c] = ft; end
            end
            close(dt)
        end
    end
    (inst, bounds, filt, tm, config, fr, ng, nc, np)
end

function step!(inst, bounds, filt, tm, config, fr, ng, nc)
    filt_ia = CLM.clump_filter_inactive_and_active
    step_start = tm.current_date
    CLM.advance_timestep!(tm)
    CLM.read_forcing_step!(fr, inst.atm2lnd, step_start, ng, nc)
    CLM.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
    calday = CLM.get_curr_calday(tm)
    (declin, _) = CLM.compute_orbital(calday)
    (declinm1, _) = CLM.compute_orbital(calday - Float64(DTIME) / CLM.SECSPDAY)
    CLM.init_daylength!(inst.gridcell, declin, declinm1, CLM.ORB_OBLIQR_DEFAULT, 1:ng)
    nextsw = calday + Float64(DTIME) / CLM.SECSPDAY
    (yr, mon, d, tod) = CLM.get_curr_date(tm)
    CLM.clm_drv!(config, inst, filt, filt_ia, bounds, true, nextsw, declin, declin,
        CLM.ORB_OBLIQR_DEFAULT, false, false, "", false;
        nstep=tm.nstep, is_first_step=(tm.nstep == 1),
        is_beg_curr_day=CLM.is_beg_curr_day(tm), is_end_curr_day=CLM.is_end_curr_day(tm),
        is_beg_curr_year=CLM.is_beg_curr_year(tm), dtime=Float64(DTIME),
        mon=mon, day=d, photosyns=inst.photosyns)
    CLM.lnd2atm!(bounds, inst)
end

# collect every Float64 array field of every sub-struct of inst, flattened
function snapshot(inst)
    out = Dict{String,Vector{Float64}}()
    for topname in propertynames(inst)
        sub = getproperty(inst, topname)
        (sub isa Number || sub isa AbstractString || sub === nothing) && continue
        for fname in propertynames(sub)
            v = try getproperty(sub, fname) catch; continue end
            if v isa AbstractArray{Float64}
                out["$topname.$fname"] = vec(copy(Float64.(v)))
            end
        end
    end
    out
end

println("building EXACT run, warmup=$WARMUP ...")
(i1, b1, f1, t1, c1, fr1, ng, nc, np) = build()
for s in 1:WARMUP; step!(i1, b1, f1, t1, c1, fr1, ng, nc); end
step!(i1, b1, f1, t1, c1, fr1, ng, nc)          # the probe step, EXACT
snapA = snapshot(i1)

println("building SMOOTHED run (same warmup, exact) ...")
(i2, b2, f2, t2, c2, fr2, _, _, _) = build()
for s in 1:WARMUP; step!(i2, b2, f2, t2, c2, fr2, ng, nc); end   # warmup EXACT (identical state)
CLM.SMOOTH_MODE[] = :always                                      # ONLY the probe step is smoothed
step!(i2, b2, f2, t2, c2, fr2, ng, nc)
snapB = snapshot(i2)
CLM.SMOOTH_MODE[] = :auto

rows = Tuple{Float64,String,Float64,Float64}[]
for (k, a) in snapA
    b = get(snapB, k, nothing); b === nothing && continue
    length(a) == length(b) || continue
    num = 0.0; den = 0.0; amax = 0.0; bmax = 0.0
    for i in eachindex(a)
        (isfinite(a[i]) && isfinite(b[i])) || continue
        (abs(a[i]) > 1e29 || abs(b[i]) > 1e29) && continue
        num = max(num, abs(a[i]-b[i]))
        den = max(den, abs(a[i]))
        amax = max(amax, abs(a[i])); bmax = max(bmax, abs(b[i]))
    end
    num == 0.0 && continue
    rel = num / max(den, 1e-30)
    push!(rows, (rel, k, num, den))
end
sort!(rows, by = r -> -r[1])
println("\n", "="^92)
println("  ONE SMOOTHED STEP vs ONE EXACT STEP from an IDENTICAL warmed-up state")
println("  ranked by max |Δ| / max|exact|   (top 40)")
println("="^92)
@printf("  %-46s %10s %14s %14s\n", "field", "rel", "max|delta|", "max|exact|")
for (rel, k, num, den) in first(rows, 40)
    @printf("  %-46s %10.3e %14.6e %14.6e\n", k, rel, num, den)
end
println("="^92)
println("  fields compared: ", length(snapA), "   diverging: ", length(rows))
