# =============================================================================
# run_clm_streamflow.jl — run CLM.jl on a (lumped) basin and compare simulated
# streamflow to a gauge. The streamflow-comparison machinery for the eval-domain
# suite: runs the biogeophysical driver over an eval period, aggregates the
# gridcell total runoff to basin discharge (runoff x area), and (where a gauge
# exists) scores it against the observed daily hydrograph (KGE + NSE + a CSV).
#
# For a LUMPED basin, discharge ~ basin-mean runoff x area (no channel routing
# needed beyond an optional constant lag). Q[m3/s] = runoff[mm/s] * area[m2] / 1000.
# Basin-mean runoff = sum_c qflx_runoff_col[c] * col.wtgcell[c] (area weights).
#
# Domain-parameterized: a registry (DOMAINS) holds per-domain run window + CLM
# flags; basin area auto-detects from the river-basins shapefile (GRU_area),
# forcing + observations auto-detect from the Symfluence domain layout. Domains
# whose inputs are not yet built are skipped cleanly. Score is computed only
# where a processed gauge CSV is present (simulated hydrograph + water balance
# are always written).
#
# Usage:
#   julia +1.12 --project=. scripts/run_clm_streamflow.jl              # default: Bow
#   DOMAIN=Iceland_Jokulsa_Fjollum julia +1.12 --project=. scripts/run_clm_streamflow.jl
#   julia +1.12 --project=. scripts/run_clm_streamflow.jl Alps_Massa_Aletsch_CH
# =============================================================================
using CLM, NCDatasets, Dates, Printf, Statistics

const ROOT = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data"
const SNOWOPT = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata/snicar_optics_5bnd_c013122.nc"
const SNOWAGE = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata/snicar_drdt_bst_fit_60_c070416.nc"

# --- per-domain config registry -------------------------------------------------
# area_km2 / forcing / obs default to auto-detection (nothing) from the Symfluence
# domain layout; set explicitly to override. spinup_end years are discarded before
# scoring. forcing = a filename in data/forcing/CLM_input, or :auto.
const DOMAINS = Dict(
    "Bow_at_Banff_lumped" => (
        run_start=DateTime(2002,1,1), run_end=DateTime(2004,12,31),
        spinup_end=DateTime(2002,12,31), forcing="clmforc.2002_2004.nc",
        area_km2=2210.0, int_snow_max=3113.2227, baseflow_scalar=0.0022119554),
    # Eval domains (Symfluence builds, 2008-2018 ERA5). Inputs auto-detect once built.
    "Boreal_Krycklan_Sweden" => (
        run_start=DateTime(2008,1,1), run_end=DateTime(2018,12,31),
        spinup_end=DateTime(2009,12,31), forcing=:auto,
        area_km2=nothing, int_snow_max=2000.0, baseflow_scalar=0.0022119554),
    "Arctic_Abisko_Sweden" => (
        run_start=DateTime(2008,1,1), run_end=DateTime(2018,12,31),
        spinup_end=DateTime(2009,12,31), forcing=:auto,
        area_km2=nothing, int_snow_max=2000.0, baseflow_scalar=0.0022119554),
    "Mediterranean_Tagus_Spain" => (
        run_start=DateTime(2008,1,1), run_end=DateTime(2018,12,31),
        spinup_end=DateTime(2009,12,31), forcing=:auto,
        area_km2=nothing, int_snow_max=2000.0, baseflow_scalar=0.0022119554),
    "Alps_Massa_Aletsch_CH" => (
        run_start=DateTime(2008,1,1), run_end=DateTime(2018,12,31),
        spinup_end=DateTime(2009,12,31), forcing=:auto,
        area_km2=nothing, int_snow_max=2000.0, baseflow_scalar=0.0022119554),
    "Urban_DeadRun_Baltimore" => (
        run_start=DateTime(2008,1,1), run_end=DateTime(2018,12,31),
        spinup_end=DateTime(2009,12,31), forcing=:auto,
        area_km2=nothing, int_snow_max=2000.0, baseflow_scalar=0.0022119554),
    "Iceland_Jokulsa_Fjollum" => (
        run_start=DateTime(2015,1,1), run_end=DateTime(2018,12,31),
        spinup_end=DateTime(2015,12,31), forcing=:auto,
        area_km2=nothing, int_snow_max=2000.0, baseflow_scalar=0.0022119554),
)

# --- input resolution from the Symfluence domain layout ------------------------
domain_dir(name) = joinpath(ROOT, "domain_$(name)")
_settings(name) = joinpath(domain_dir(name), "settings", "CLM", "parameters")
surfdata_path(name) = joinpath(_settings(name), "surfdata_clm.nc")
params_path(name)   = joinpath(_settings(name), "clm5_params.nc")

# Basin area (km^2): read GRU_area (m^2) from the river-basins shapefile via a
# python/geopandas one-liner (Symfluence's env has geopandas); fall back to the
# registry value if unavailable.
function detect_area_km2(name, fallback)
    shpdir = joinpath(domain_dir(name), "shapefiles", "river_basins")
    isdir(shpdir) || return fallback
    py = """
import glob, sys
try:
    import geopandas as gpd
    f = sorted(glob.glob(r'$(shpdir)/*.shp'))
    if not f: sys.exit(0)
    g = gpd.read_file(f[0])
    col = next((c for c in g.columns if c.lower() in ('gru_area','area','basin_area')), None)
    if col is None: sys.exit(0)
    print(float(g[col].sum()))
except Exception:
    sys.exit(0)
"""
    try
        out = strip(read(`python3 -c $py`, String))
        isempty(out) && return fallback
        a_m2 = parse(Float64, out)
        return a_m2 > 0 ? a_m2 / 1.0e6 : fallback
    catch
        return fallback
    end
end

# Forcing file: explicit name if given + present; else prefer a merged clmforc
# spanning the run, else fall back to per-year (clmforc.<year>.nc, reopened per year).
function resolve_forcing(name, cfg)
    cdir = joinpath(domain_dir(name), "data", "forcing", "CLM_input")
    isdir(cdir) || return (:missing, "")
    if cfg.forcing isa String
        p = joinpath(cdir, cfg.forcing)
        isfile(p) && return (:single, p)
    end
    files = filter(f -> startswith(f, "clmforc") && endswith(f, ".nc"), readdir(cdir))
    isempty(files) && return (:missing, "")
    yrs = "$(year(cfg.run_start))_$(year(cfg.run_end))"
    merged = findfirst(f -> occursin(yrs, f), files)
    merged !== nothing && return (:single, joinpath(cdir, files[merged]))
    # per-year if every year in the window has a clmforc.<year>.nc
    if all(y -> isfile(joinpath(cdir, "clmforc.$(y).nc")),
           year(cfg.run_start):year(cfg.run_end))
        return (:peryear, cdir)
    end
    # otherwise just take the longest-named (most-merged) file
    return (:single, joinpath(cdir, files[argmax(length.(files))]))
end

# Processed gauge CSV (datetime, discharge_cms), if a gauge was acquired.
function detect_obs(name)
    pdir = joinpath(domain_dir(name), "data", "observations", "streamflow", "preprocessed")
    isdir(pdir) || return nothing
    csv = filter(f -> endswith(f, ".csv"), readdir(pdir))
    isempty(csv) ? nothing : joinpath(pdir, csv[1])
end

# --- metrics -------------------------------------------------------------------
function kge(sim, obs)
    m = isfinite.(sim) .& isfinite.(obs)
    s = sim[m]; o = obs[m]
    length(o) < 10 && return (NaN, NaN, NaN, NaN)
    r = cor(s, o); a = std(s)/std(o); b = mean(s)/mean(o)
    (1 - sqrt((r-1)^2 + (a-1)^2 + (b-1)^2), r, a, b)
end
nse(sim, obs) = (m = isfinite.(sim) .& isfinite.(obs); s=sim[m]; o=obs[m];
                 length(o)<10 ? NaN : 1 - sum((s.-o).^2)/sum((o.-mean(o)).^2))

function read_obs_daily(path)
    daily = Dict{Date,Vector{Float64}}()
    (path === nothing || !isfile(path)) && return daily
    for (i, ln) in enumerate(eachline(path))
        i == 1 && continue
        parts = split(ln, ',')
        length(parts) < 2 && continue
        dt = tryparse(DateTime, strip(parts[1]), dateformat"yyyy-mm-dd HH:MM:SS")
        dt === nothing && (dt = tryparse(DateTime, strip(parts[1])))
        dt === nothing && continue
        q = tryparse(Float64, strip(parts[2]))
        (q === nothing || !isfinite(q)) && continue
        push!(get!(daily, Date(dt), Float64[]), q)
    end
    Dict(d => mean(v) for (d, v) in daily)
end

# --- run one domain ------------------------------------------------------------
function run_domain(name::String)
    haskey(DOMAINS, name) || error("unknown domain '$name'; known: $(join(sort(collect(keys(DOMAINS))), ", "))")
    cfg = DOMAINS[name]
    FS, FP = surfdata_path(name), params_path(name)
    if !(isfile(FS) && isfile(FP))
        @info "[$name] CLM inputs not built yet — skipping" surfdata=FS
        return nothing
    end
    fmode, fpath = resolve_forcing(name, cfg)
    if fmode == :missing
        @info "[$name] forcing not built yet — skipping"
        return nothing
    end
    area_km2 = detect_area_km2(name, cfg.area_km2)
    area_km2 === nothing && error("[$name] basin area unknown (no shapefile + no registry value)")
    obs_path = detect_obs(name)

    CLM.rooting_profile_config.rooting_profile_method_water  = CLM.JACKSON_1996_ROOT
    CLM.rooting_profile_config.rooting_profile_method_carbon = CLM.JACKSON_1996_ROOT
    (inst, bounds, filt, tm) = CLM.clm_initialize!(; fsurdat=FS, paramfile=FP,
        start_date=cfg.run_start, dtime=3600, use_cn=false, use_luna=false,
        use_bedrock=true, use_aquifer_layer=false, h2osfcflag=0,
        fsnowoptics=SNOWOPT, fsnowaging=SNOWAGE, int_snow_max=cfg.int_snow_max)
    ng, nc, np = bounds.endg, bounds.endc, bounds.endp
    CLM.init_soil_hydrology_config(baseflow_scalar = cfg.baseflow_scalar)
    config  = CLM.CLMDriverConfig(use_cn=false, use_aquifer_layer=false,
                                  use_hydrstress=false, use_luna=false)
    filt_ia = CLM.clump_filter_inactive_and_active
    wf = inst.water.waterfluxbulk_inst.wf

    @printf("[%s] nc=%d np=%d area=%.0f km2  %s..%s  forcing=%s(%s)  obs=%s\n",
        name, nc, np, area_km2, Date(cfg.run_start), Date(cfg.run_end),
        fmode, basename(fpath), obs_path === nothing ? "NONE" : basename(obs_path))

    basin_runoff() = (s = 0.0; @inbounds for c in 1:nc
        v = wf.qflx_runoff_col[c]; isfinite(v) && (s += v * inst.column.wtgcell[c]); end; s)

    # forcing reader: single merged file, or reopened per calendar year
    fr = CLM.ForcingReader()
    cur_fyear = Ref(0)
    open_forcing!(y) = begin
        fmode == :peryear || return
        if cur_fyear[] != y
            cur_fyear[] != 0 && CLM.forcing_reader_close!(fr)
            CLM.forcing_reader_init!(fr, joinpath(fpath, "clmforc.$(y).nc")); cur_fyear[] = y
        end
    end
    fmode == :single && CLM.forcing_reader_init!(fr, fpath)

    sim_daily = Dict{Date,Vector{Float64}}()
    ndays = Dates.value(Date(cfg.run_end) - Date(cfg.run_start)) + 1
    nbad = 0
    for d in 1:ndays
        for h in 1:24
            cur = cfg.run_start + Day(d-1) + Hour(h-1)
            open_forcing!(year(cur))
            calday = CLM.get_curr_calday(tm); (declin, _) = CLM.compute_orbital(calday)
            (declinm1, _) = CLM.compute_orbital(calday - 3600.0/CLM.SECSPDAY)
            nextsw = calday + 3600.0/CLM.SECSPDAY
            CLM.init_daylength!(inst.gridcell, declin, declinm1, CLM.ORB_OBLIQR_DEFAULT, 1:ng)
            CLM.advance_timestep!(tm)
            CLM.read_forcing_step!(fr, inst.atm2lnd, cur, ng, nc)
            CLM.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
            (yr, mon, dy, tod) = CLM.get_curr_date(tm)
            CLM.clm_drv!(config, inst, filt, filt_ia, bounds, true, nextsw, declin, declin,
                CLM.ORB_OBLIQR_DEFAULT, false, false, "", false;
                nstep=tm.nstep, is_first_step=(d==1 && h==1),
                is_beg_curr_day=CLM.is_beg_curr_day(tm), is_end_curr_day=CLM.is_end_curr_day(tm),
                is_beg_curr_year=CLM.is_beg_curr_year(tm), dtime=3600.0, mon=mon, day=dy,
                photosyns=inst.photosyns)
            r = basin_runoff(); isfinite(r) || (nbad += 1)
            push!(get!(sim_daily, Date(cur), Float64[]), r)
        end
        d % 365 == 0 && @printf("  day %4d (%s): runoff=%.3f mm/day\n", d,
            Date(cfg.run_start)+Day(d-1), mean(get(sim_daily, Date(cfg.run_start)+Day(d-1), [0.0]))*86400)
    end
    CLM.forcing_reader_close!(fr)

    # daily mean runoff (mm/s) -> discharge (m3/s)
    A = area_km2 * 1.0e6
    sim_q = Dict(dt => mean(v) * A / 1000.0 for (dt, v) in sim_daily)
    obs_q = read_obs_daily(obs_path)

    out = joinpath(@__DIR__, "..", "streamflow_compare_$(name).csv")
    if isempty(obs_q)
        # no gauge: write the simulated hydrograph + report mean water balance
        days = sort([d for d in keys(sim_q) if d > Date(cfg.spinup_end)])
        open(out, "w") do io
            println(io, "date,sim_cms")
            for dd in days; @printf(io, "%s,%.4f\n", dd, sim_q[dd]); end
        end
        simv = [sim_q[d] for d in days]
        @printf("[%s] nbad=%d  scoring days=%d (no gauge)  sim mean=%.2f cms (%.1f mm/yr)\n",
            name, nbad, length(days), mean(simv), mean(simv)*86400*365.25/A*1000)
        println("  wrote ", out, "  (simulated only — no observations)")
        return (; name, scored=false, sim_mean=mean(simv))
    end
    days = sort([d for d in keys(sim_q) if d > Date(cfg.spinup_end) && haskey(obs_q, d)])
    sim = [sim_q[d] for d in days]; obs = [obs_q[d] for d in days]
    @printf("[%s] nbad=%d  scoring days=%d (%s..%s)\n", name, nbad, length(days),
        isempty(days) ? "-" : string(days[1]), isempty(days) ? "-" : string(days[end]))
    if length(days) >= 10
        k, r, a, b = kge(sim, obs); n = nse(sim, obs)
        @printf("  sim mean=%.1f cms  obs mean=%.1f cms\n", mean(sim), mean(obs))
        @printf("  KGE=%.3f (r=%.2f alpha=%.2f beta=%.2f)  NSE=%.3f\n", k, r, a, b, n)
        open(out, "w") do io
            println(io, "date,sim_cms,obs_cms")
            for (i, dd) in enumerate(days); @printf(io, "%s,%.4f,%.4f\n", dd, sim[i], obs[i]); end
        end
        println("  wrote ", out)
        return (; name, scored=true, kge=k, nse=n, r, alpha=a, beta=b,
                sim_mean=mean(sim), obs_mean=mean(obs))
    else
        println("  too few overlapping days to score (sim=$(length(keys(sim_q))), obs=$(length(keys(obs_q))))")
        return (; name, scored=false, sim_mean=mean([sim_q[d] for d in keys(sim_q)]))
    end
end

# Run every registry domain whose inputs are built; collate a summary table.
function run_all()
    results = NamedTuple[]
    for name in sort(collect(keys(DOMAINS)))
        try
            r = run_domain(name)
            r === nothing || push!(results, r)
        catch e
            @warn "[$name] run failed" exception=(e, catch_backtrace())
        end
    end
    println("\n", "="^78, "\nSTREAMFLOW SUMMARY (", length(results), " domains run)\n", "="^78)
    @printf("%-28s %8s %8s %8s %8s %10s\n", "domain", "KGE", "NSE", "r", "beta", "sim_cms")
    for r in results
        if get(r, :scored, false)
            @printf("%-28s %8.3f %8.3f %8.2f %8.2f %10.2f\n",
                r.name, r.kge, r.nse, r.r, r.beta, r.sim_mean)
        else
            @printf("%-28s %8s %8s %8s %8s %10.2f   (no gauge)\n",
                r.name, "-", "-", "-", "-", r.sim_mean)
        end
    end
    return results
end

function main()
    arg = !isempty(ARGS) ? ARGS[1] : get(ENV, "DOMAIN", "Bow_at_Banff_lumped")
    lowercase(arg) == "all" ? run_all() : run_domain(arg)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
