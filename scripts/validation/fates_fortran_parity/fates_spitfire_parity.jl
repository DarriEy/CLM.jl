# fates_spitfire_parity.jl — FATES SPITFIRE (fates_spitfire_mode = 1) vs Fortran CTSM.
#
# Steps the Julia port on the fates_parity_1pt site (Bow at Banff, 2002-07-01, 3600 s,
# FATES + use_fates_bgc, all 14 PFTs, SPITFIRE scalar_lightning) and scores the DAILY
# fire diagnostics against a daily-mean h0 history of the same single-point Fortran
# case (stage_spitfire_run.sh):
#
#   FATES_NESTEROV_INDEX       fire-weather (Nesterov) index           [-]
#   FATES_FUEL_AMOUNT          ground fuel (no 1000-h fuels)           kgC m-2
#   FATES_FUEL_EFF_MOIST       fuel effective moisture                 m3 m-3
#   FATES_FUEL_MEF/BULKD/SAV   fuel moisture of extinction / bulk density / SAV
#   FATES_FDI, FATES_IGNITIONS fire danger index, successful ignitions
#   FATES_ROS, FATES_EFFECT_WSPEED   rate of spread, effective wind
#   FATES_FIRE_INTENSITY, FATES_BURNFRAC, FATES_FUELCONSUMED
#   FATES_FIRE_CLOSS, FATES_DISTURBANCE_RATE_FIRE, FATES_MORTALITY_FIRE_SZ (sum)
#   FATES_TVEG24, FATES_VEGC, FATES_NPATCHES
#   WIND, RH, RAIN+SNOW (PRECIP), TBOT   host daily-mean fire-weather drivers
#
# The committed extract fates_spitfire_fortran.csv lets the comparison reproduce
# without CTSM.
#
# usage:
#   julia --project=. scripts/validation/fates_fortran_parity/fates_spitfire_parity.jl \
#         [fortran_h0.nc | fates_spitfire_fortran.csv] [ndays=730] [julia_out.csv]
#   ... fates_spitfire_parity.jl export <fortran_h0.nc> [out.csv]   # refresh the extract
using CLM, Printf, Dates, NCDatasets
include(joinpath(@__DIR__, "..", "..", "fates_fortran_parity.jl"))   # build(), FORCDIR, START, DTIME

const FIRE_FIELDS = ("FATES_NESTEROV_INDEX", "FATES_FUEL_AMOUNT", "FATES_FUEL_EFF_MOIST",
    "FATES_FUEL_MEF", "FATES_FUEL_BULKD", "FATES_FUEL_SAV", "FATES_FDI", "FATES_IGNITIONS",
    "FATES_ROS", "FATES_EFFECT_WSPEED", "FATES_FIRE_INTENSITY", "FATES_BURNFRAC",
    "FATES_FUELCONSUMED", "FATES_FIRE_CLOSS", "FATES_DISTURBANCE_RATE_FIRE",
    "FATES_MORTALITY_FIRE_SZ", "FATES_TVEG24", "FATES_VEGC", "FATES_NPATCHES")
const HOST_FIELDS = ("WIND", "RH", "PRECIP", "TBOT")
# Per fuel class (twigs, small/large branches, trunks, dead leaves, live grass):
# FATES_FUEL_AMOUNT_FC / FATES_FUEL_MOISTURE_FC (level-2 history; the Julia side is
# summed from the patches exactly as update_history_dyn2 does).
const NFC = 6
const FC_FIELDS = Tuple(vcat(["FUEL_FC$k" for k in 1:NFC], ["MOIST_FC$k" for k in 1:NFC]))
const ALL_FIELDS = (FIRE_FIELDS..., HOST_FIELDS..., FC_FIELDS...)

"Value of a Julia FATES history variable (site 1; multi-dim variables summed)."
function _jhist(h, name)
    for hv in h.hvars
        hv.vname == name || continue
        isempty(hv.r81d) || return hv.r81d[1]
        isempty(hv.r82d) || return sum(hv.r82d[1, :])
        return NaN
    end
    return NaN
end

"Noleap model day-start for elapsed day `d` (0-based) from START (CESM NO_LEAP)."
function noleap_day(d::Int)
    t = Date(START); k = 0
    while k < d
        t += Day(1)
        (month(t) == 2 && day(t) == 29) && (t += Day(1))
        k += 1
    end
    return t
end

"""
    datm_time(cur) -> DateTime

The forcing record DATM reads at NO_LEAP model time `cur`. The site's clmforc.YYYY.nc
files are on the standard calendar (8784 h in 2004); DATM indexes them by elapsed time
within the year, so from 1 March of a leap year the model reads the PREVIOUS calendar
day (model 2004-03-01 <- data 2004-02-29, ..., data 2004-12-31 is never read). Verified
against the Fortran daily TBOT/WIND: aligning by calendar date instead puts every day
from 2004-03-01 one day off (TBOT ~2 K mean error); this mapping gives 7e-6 K.
"""
function datm_time(cur::DateTime)
    y = year(cur)
    isleapyear(y) || return cur
    doy0 = dayofyear(cur) - 1 - ((month(cur) > 2) ? 1 : 0)   # 0-based NO_LEAP day of year
    return DateTime(y, 1, 1) + Day(doy0) + (cur - DateTime(Date(cur)))
end

"""
Run the port `ndays` days and return one record per model day: the FATES fire history
after that day's daily step, plus the day-mean host forcing over the day's steps.
Record d (1-based) <-> the Fortran h0 record covering model day d.
"""
function run_julia(ndays::Int)
    # SPITFIRE_MODE=2..5 (default 1) reads the same synthetic uniform lightning /
    # population-density streams stage_spitfire_run.sh gives the Fortran case.
    mode = parse(Int, get(ENV, "SPITFIRE_MODE", "1"))
    fd = joinpath(get(ENV, "CESM_INPUTDATA", joinpath(DATA, "installs", "cesm-inputdata")),
                  "lnd", "clm2", "firedata")
    streams = mode > 1 ?
        (flnfm = joinpath(fd, "clmforc.UNIFORM_lnfm_0.9x1.25_CLMjl-parity.nc"),
         fhdm  = joinpath(fd, "clmforc.UNIFORM_hdm_0.9x1.25_CLMjl-parity.nc")) : (;)
    inst, bounds, filt, filt_ia, config, _tm = build(; fates_spitfire_mode = mode, streams...)
    CLM.hlm_spitfire_mode[] == mode || error("SPITFIRE mode not set")
    # FIRE_THRESHOLD=x [kW/m] overrides fates_fire_threshold, to match a Fortran
    # reference run whose parameter file lowers it (stage_spitfire_run.sh CDL=...).
    if haskey(ENV, "FIRE_THRESHOLD")
        CLM.sf_params().SF_val_fire_threshold = parse(Float64, ENV["FIRE_THRESHOLD"])
        @printf("  fates_fire_threshold overridden -> %g kW/m\n", CLM.sf_params().SF_val_fire_threshold)
    end
    ng, nc = bounds.endg, bounds.endc
    let tf = joinpath(FORCDIR, "topo_forcing.nc")
        if isfile(tf)
            NCDataset(tf, "r") do ds
                if haskey(ds, "TOPO")
                    ft = Float64(ds["TOPO"][1])
                    for g in 1:ng; inst.atm2lnd.forc_topo_grc[g] = ft; end
                    for cc in 1:nc; inst.topo.topo_col[cc] = ft; end
                end
            end
        end
    end
    a2l = inst.atm2lnd
    c = findfirst(inst.column.is_fates)
    fr = CLM.ForcingReader(); fyear = year(START)
    CLM.forcing_reader_init!(fr, joinpath(FORCDIR, "clmforc.$(fyear).nc"))
    spd = round(Int, 86400 / DTIME)
    recs = Vector{Dict{String,Float64}}()
    fmort_prev = Ref(0.0)
    nstep = 0; t0 = time()
    for d in 0:ndays-1
        dstart = noleap_day(d)
        acc = Dict(f => 0.0 for f in HOST_FIELDS); nacc = 0
        for k in 0:spd-1
            for i in ((d == 0 && k == 0) ? (0:1) : (nstep+1:nstep+1))   # step-0 replay
                cur = DateTime(dstart) + Second(round(Int, k * DTIME))
                if year(cur) != fyear
                    CLM.forcing_reader_close!(fr); fyear = year(cur)
                    CLM.forcing_reader_init!(fr, joinpath(FORCDIR, "clmforc.$(fyear).nc"))
                end
                tod = hour(cur)*3600 + minute(cur)*60 + second(cur)
                calday = dayofyear(cur) + tod / CLM.SECSPDAY
                (declin, _) = CLM.compute_orbital(calday)
                (declinm1, _) = CLM.compute_orbital(calday - DTIME / CLM.SECSPDAY)
                CLM.init_daylength!(inst.gridcell, declin, declinm1, CLM.ORB_OBLIQR_DEFAULT, 1:ng)
                CLM.read_forcing_step!(fr, a2l, datm_time(cur), ng, nc)
                CLM.downscale_forcings!(bounds, a2l, inst.column, inst.landunit, inst.topo)
                CLM.clm_drv!(config, inst, filt, filt_ia, bounds, CLM.nuopc_doalb(i),
                    calday + DTIME / CLM.SECSPDAY, declin, declin, CLM.ORB_OBLIQR_DEFAULT,
                    false, false, "", false;
                    nstep=i, is_first_step=(i <= 1),
                    is_beg_curr_day=(i > 0 && tod == 0), is_end_curr_day=false,
                    is_beg_curr_year=(i > 0 && tod == 0 && dayofyear(cur) == 1),
                    dtime=DTIME, mon=month(cur), day=day(cur), secs=tod,
                    jday=dayofyear(cur), year=year(cur), photosyns=inst.photosyns)
                nstep = i
                if i > 0
                    g = inst.column.gridcell[c]
                    acc["WIND"] += a2l.forc_wind_grc[g]
                    acc["RH"]   += a2l.forc_rh_grc[g]
                    acc["PRECIP"] += a2l.forc_rain_downscaled_col[c] + a2l.forc_snow_downscaled_col[c]
                    acc["TBOT"] += a2l.forc_t_downscaled_col[c]
                    nacc += 1
                end
            end
        end
        r = Dict{String,Float64}(f => _jhist(inst.fates.hist, f) for f in FIRE_FIELDS)
        # FATES_MORTALITY_FIRE_SZ is a level-2 (group_dyna_complx) history field the
        # port's level-1 harness does not fill; sum the site rates it is built from
        # (update_history_dyn2: (fmort_rate_canopy + fmort_rate_ustory) / m2_per_ha).
        let st = inst.fates.sites[1]
            # dyn2 also ZEROES these site rates after writing them; the level-1 harness
            # never does, so they accumulate: report the day's increment.
            cum = (sum(st.fmort_rate_canopy) + sum(st.fmort_rate_ustory)) / 1.0e4
            r["FATES_MORTALITY_FIRE_SZ"] = cum - fmort_prev[]
            fmort_prev[] = cum
            fc = zeros(NFC); mc = zeros(NFC)
            cp = st.oldest_patch
            while cp !== nothing
                if cp.nocomp_pft_label != CLM.nocomp_bareground && cp.fuel !== nothing
                    w = cp.area / CLM.area
                    for k in 1:NFC
                        fc[k] += cp.fuel.frac_loading[k] * cp.fuel.non_trunk_loading * w
                        mc[k] += cp.fuel.effective_moisture[k] * w
                    end
                end
                cp = cp.younger
            end
            for k in 1:NFC; r["FUEL_FC$k"] = fc[k]; r["MOIST_FC$k"] = mc[k]; end
        end
        for f in HOST_FIELDS; r[f] = acc[f] / max(nacc, 1); end
        push!(recs, r)
        if (d + 1) % 30 == 0
            @printf("  day %4d %s  NI=%9.1f fuel=%.4f burnfrac=%.3e npatch=%d  %.0fs\n",
                    d + 1, string(dstart), r["FATES_NESTEROV_INDEX"], r["FATES_FUEL_AMOUNT"],
                    r["FATES_BURNFRAC"], round(Int, r["FATES_NPATCHES"]), time() - t0)
        end
    end
    CLM.forcing_reader_close!(fr)
    return recs
end

"Fortran daily series from the h0 history (`.nc`) or the committed CSV extract."
function read_fortran(path)
    if endswith(path, ".csv")
        lines = readlines(path)
        hdr = split(lines[1], ',')
        cols = [Float64[] for _ in hdr]
        for l in lines[2:end]
            isempty(l) && continue
            for (i, v) in enumerate(split(l, ','))
                push!(cols[i], parse(Float64, v))
            end
        end
        return Dict(String(h) => cols[i] for (i, h) in enumerate(hdr) if h != "day")
    end
    NCDataset(path, "r") do ds
        out = Dict{String,Vector{Float64}}()
        for f in FIRE_FIELDS
            haskey(ds, f) || continue
            a = Array(ds[f])
            out[f] = ndims(a) == 2 ? Float64.(vec(a[1, :])) :
                     Float64.(vec(sum(a[:, 1, :]; dims = 1)))       # (lev, lnd, time)
        end
        for (hv, pre) in (("FATES_FUEL_AMOUNT_FC", "FUEL_FC"), ("FATES_FUEL_MOISTURE_FC", "MOIST_FC"))
            haskey(ds, hv) || continue
            a = Array(ds[hv])                                  # (lnd, fuel, time)
            for k in 1:NFC; out[pre * string(k)] = Float64.(vec(a[1, k, :])); end
        end
        out["WIND"] = Float64.(vec(ds["WIND"][1, :]))
        out["RH"]   = Float64.(vec(ds["RH"][1, :]))
        out["TBOT"] = Float64.(vec(ds["TBOT"][1, :]))
        out["PRECIP"] = Float64.(vec(ds["RAIN"][1, :])) .+ Float64.(vec(ds["SNOW"][1, :]))
        # Record 1 is the zero-length nstep=0 (NUOPC init) sample, time_bounds [0,0];
        # record k+1 is the mean over model day k.
        tb = ds["time_bounds"].var[:, 1]
        if tb[1] == tb[2]
            for k in keys(out); out[k] = out[k][2:end]; end
        end
        out
    end
end

"""
Score Julia vs Fortran per field over days 1..n: mean and max relative error with a
floor of 1e-3 x the field's peak |F| (fire fields are exactly zero on many days), the
number of days each side reports a fire (BURNFRAC > 0), and the series totals.
"""
function score(jl, ft; io = stdout, nmax = typemax(Int))
    n = min(length(jl), minimum(length(v) for v in values(ft)), nmax)
    @printf(io, "  %-28s %13s %13s %10s %10s\n", "field (days 1..$n)", "Fortran mean",
            "Julia mean", "max rel", "mean rel")
    for f in ALL_FIELDS
        haskey(ft, f) || continue
        F = ft[f][1:n]; J = [jl[k][f] for k in 1:n]
        floor_ = max(maximum(abs, F) * 1e-3, 1e-30)
        rel = [abs(J[k] - F[k]) / max(abs(F[k]), floor_) for k in 1:n]
        @printf(io, "  %-28s %13.5e %13.5e %10.3e %10.3e\n", f, sum(F) / n, sum(J) / n,
                maximum(rel), sum(rel) / n)
    end
    if haskey(ft, "FATES_BURNFRAC")
        fd = count(>(0), ft["FATES_BURNFRAC"][1:n]); jd = count(k -> jl[k]["FATES_BURNFRAC"] > 0, 1:n)
        both = count(k -> ft["FATES_BURNFRAC"][k] > 0 && jl[k]["FATES_BURNFRAC"] > 0, 1:n)
        @printf(io, "  fire days: Fortran %d, Julia %d, both %d\n", fd, jd, both)
        @printf(io, "  total burned fraction: Fortran %.5e, Julia %.5e (per day sums x 86400)\n",
                sum(ft["FATES_BURNFRAC"][1:n]) * 86400, sum(jl[k]["FATES_BURNFRAC"] for k in 1:n) * 86400)
    end
end

function export_fortran_csv(nc_path, csv_path)
    ft = read_fortran(nc_path)
    names = [f for f in ALL_FIELDS if haskey(ft, f)]
    n = length(ft[names[1]])
    open(csv_path, "w") do io
        println(io, join(["day"; names], ","))
        for k in 1:n
            println(io, join([string(k); [repr(ft[f][k]) for f in names]], ","))
        end
    end
    return csv_path
end

function write_csv(path, jl, ft)
    n = min(length(jl), minimum(length(v) for v in values(ft)))
    open(path, "w") do io
        println(io, join(["day"; [f * "_jl" for f in ALL_FIELDS]; [f * "_ft" for f in ALL_FIELDS]], ","))
        for k in 1:n
            println(io, join([string(k); [repr(jl[k][f]) for f in ALL_FIELDS];
                              [repr(haskey(ft, f) ? ft[f][k] : NaN) for f in ALL_FIELDS]], ","))
        end
    end
end

const FORTRAN_CSV = joinpath(@__DIR__, "fates_spitfire_fortran.csv")

function spitfire_main(args)
    if !isempty(args) && args[1] == "export"
        export_fortran_csv(args[2], length(args) >= 3 ? args[3] : FORTRAN_CSV)
        return 0
    end
    ref = isempty(args) ? FORTRAN_CSV : args[1]
    ft = read_fortran(ref)
    ndays = length(args) >= 2 ? parse(Int, args[2]) : length(first(values(ft)))
    jl = run_julia(ndays)
    println("\n=== FATES SPITFIRE (mode 1): Julia vs Fortran (", ndays, " days) ===")
    score(jl, ft)
    # The FATES vegetation trajectories separate after ~5 weeks (patch spawning and
    # litter inputs; independent of fire -- the no-fire run diverges identically), so
    # the fire chain is also scored over the window where the stand is still matched.
    println("\n--- matched-vegetation window ---")
    score(jl, ft; nmax = 40)
    length(args) >= 3 && write_csv(args[3], jl, ft)
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(spitfire_main(ARGS))
end
