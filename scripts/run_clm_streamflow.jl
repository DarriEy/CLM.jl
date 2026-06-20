# =============================================================================
# run_clm_streamflow.jl — run CLM.jl on a (lumped) basin and compare simulated
# streamflow to a gauge. The streamflow-comparison machinery for the eval-domain
# suite: runs the biogeophysical driver over an eval period, aggregates the
# gridcell total runoff to basin discharge (runoff x area), and scores it against
# the observed daily hydrograph (KGE + NSE + a written CSV).
#
# For a LUMPED basin, discharge ~ basin-mean runoff x area (no channel routing
# needed beyond an optional constant lag). Q[m3/s] = runoff[mm/s] * area[m2] / 1000.
# Basin-mean runoff = sum_c qflx_runoff_col[c] * col.wtgcell[c] (area weights).
#
# Default config = Bow at Banff (the one fully-built domain with forcing + a gauge),
# to validate the machinery now; the eval domains (Iceland/boreal/...) plug in the
# same way once their Symfluence builds finish.
#
# Usage: julia +1.12 --project=. scripts/run_clm_streamflow.jl
# =============================================================================
using CLM, NCDatasets, Dates, Printf, Statistics

# ---- domain config (Bow at Banff default) ----
const DOM   = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped"
const FS    = joinpath(DOM, "settings/CLM/parameters/surfdata_clm.nc")
const FP    = joinpath(DOM, "settings/CLM/parameters/clm5_params.nc")
const FFORC = joinpath(DOM, "data/forcing/CLM_input/clmforc.2002_2004.nc")
const OBS   = joinpath(DOM, "data/observations/streamflow/preprocessed/Bow_at_Banff_lumped_streamflow_processed.csv")
const AREA_KM2 = 2210.0                         # Bow at Banff drainage area
const SPINUP_END = DateTime(2002, 12, 31)       # discard spinup before scoring
const RUN_START  = DateTime(2002, 1, 1)
const RUN_END    = DateTime(2004, 12, 31)
const SNOWOPT = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata/snicar_optics_5bnd_c013122.nc"
const SNOWAGE = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata/snicar_drdt_bst_fit_60_c070416.nc"

# ---- metrics ----
function kge(sim, obs)
    m = .!ismissing.(obs) .& isfinite.(sim) .& isfinite.(Float64.(replace(obs, missing=>NaN)))
    s = sim[m]; o = Float64.(obs[m])
    length(o) < 10 && return (NaN, NaN, NaN, NaN)
    r = cor(s, o); a = std(s)/std(o); b = mean(s)/mean(o)
    (1 - sqrt((r-1)^2 + (a-1)^2 + (b-1)^2), r, a, b)
end
nse(sim, obs) = (m = .!ismissing.(obs) .& isfinite.(sim); s=sim[m]; o=Float64.(obs[m]);
                 length(o)<10 ? NaN : 1 - sum((s.-o).^2)/sum((o.-mean(o)).^2))

function read_obs_daily(path)
    daily = Dict{Date,Vector{Float64}}()
    isfile(path) || return daily
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

function main()
    CLM.rooting_profile_config.rooting_profile_method_water  = CLM.JACKSON_1996_ROOT
    CLM.rooting_profile_config.rooting_profile_method_carbon = CLM.JACKSON_1996_ROOT
    (inst, bounds, filt, tm) = CLM.clm_initialize!(; fsurdat=FS, paramfile=FP,
        start_date=RUN_START, dtime=3600, use_cn=false, use_luna=false,
        use_bedrock=true, use_aquifer_layer=false, h2osfcflag=0,
        fsnowoptics=SNOWOPT, fsnowaging=SNOWAGE, int_snow_max=3113.2227)
    ng, nc, np = bounds.endg, bounds.endc, bounds.endp
    CLM.init_soil_hydrology_config(baseflow_scalar = 0.0022119554)
    config  = CLM.CLMDriverConfig(use_cn=false, use_aquifer_layer=false,
                                  use_hydrstress=false, use_luna=false)
    filt_ia = CLM.clump_filter_inactive_and_active
    fr = CLM.ForcingReader(); CLM.forcing_reader_init!(fr, FFORC)
    wf = inst.water.waterfluxbulk_inst.wf
    @printf("Streamflow run: nc=%d np=%d area=%.0f km2  %s..%s\n",
        nc, np, AREA_KM2, Date(RUN_START), Date(RUN_END))

    # basin-mean runoff (mm/s), area-weighted over columns
    basin_runoff() = (s = 0.0; @inbounds for c in 1:nc
        v = wf.qflx_runoff_col[c]; isfinite(v) && (s += v * inst.column.wtgcell[c]); end; s)

    sim_daily = Dict{Date,Vector{Float64}}()
    ndays = Dates.value(Date(RUN_END) - Date(RUN_START)) + 1
    nbad = 0
    for d in 1:ndays
        for h in 1:24
            cur = RUN_START + Day(d-1) + Hour(h-1)
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
        d % 120 == 0 && @printf("  day %4d (%s): runoff=%.3f mm/day\n", d,
            Date(RUN_START)+Day(d-1), mean(get(sim_daily, Date(RUN_START)+Day(d-1), [0.0]))*86400)
    end
    CLM.forcing_reader_close!(fr)

    # daily mean runoff (mm/s) -> discharge (m3/s)
    A = AREA_KM2 * 1.0e6
    sim_q = Dict(dt => mean(v) * A / 1000.0 for (dt, v) in sim_daily)   # mm/s * m2 / (mm/m) = m3/s
    obs_q = read_obs_daily(OBS)

    # align over the scoring window (post-spinup)
    days = sort([d for d in keys(sim_q) if d > Date(SPINUP_END) && haskey(obs_q, d)])
    sim = [sim_q[d] for d in days]; obs = [obs_q[d] for d in days]
    @printf("\nnbad(non-finite runoff steps)=%d  scoring days=%d (%s..%s)\n",
        nbad, length(days), isempty(days) ? "-" : string(days[1]), isempty(days) ? "-" : string(days[end]))
    if length(days) >= 10
        k, r, a, b = kge(sim, Vector{Union{Missing,Float64}}(obs)); n = nse(sim, Vector{Union{Missing,Float64}}(obs))
        @printf("  sim mean=%.1f cms  obs mean=%.1f cms\n", mean(sim), mean(obs))
        @printf("  KGE=%.3f (r=%.2f alpha=%.2f beta=%.2f)  NSE=%.3f\n", k, r, a, b, n)
        out = joinpath(@__DIR__, "..", "streamflow_compare.csv")
        open(out, "w") do io
            println(io, "date,sim_cms,obs_cms")
            for (i, dd) in enumerate(days); @printf(io, "%s,%.4f,%.4f\n", dd, sim[i], obs[i]); end
        end
        println("  wrote ", out)
    else
        println("  too few overlapping days to score (sim days=$(length(keys(sim_q))), obs days=$(length(keys(obs_q))))")
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
