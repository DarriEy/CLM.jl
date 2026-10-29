# fates_cnp_parity.jl — FATES CNP (fates_parteh_mode = 2) vs Fortran CTSM.
#
# Steps the Julia port on the fates_parity_1pt site (Bow at Banff, 2002-07-01, 3600 s,
# use_fates_bgc, CENTURYKoven2013) in CNP mode with the CTSM FatesColdPRT2 settings
# (prescribed N/P uptake = 1 for all PFTs) and scores, step by step, against an hourly
# h0 history of the same single-point Fortran case (stage_cnp_run.sh):
#
#   FATES plant N/P pools  FATES_VEGN/STOREN/LEAFN/FROOTN/SAPWOODN/REPRON (+ P), *_TF
#   FATES N/P fluxes       FATES_NH4UPTAKE/NO3UPTAKE/NDEMAND/NEFFLUX/NFIX_SYM,
#                          FATES_PUPTAKE/PDEMAND/PEFFLUX, FATES_EXCESS_RESP
#   FATES carbon           FATES_VEGC/LEAFC/STOREC/FROOTC/SAPWOODC/REPROC, NPP, GPP,
#                          FATES_NCOHORTS, FATES_L2FR, FATES_LAI
#   host soil BGC          SMINN, SMIN_NH4, SMIN_NO3, HR, TOTLITC, TOTSOMC
#
# The FATES_* values are read from the port's own FATES history buffers (the same
# registry names CTSM writes), the host values from the soil-BGC state.
#
# usage:
#   julia --project=. scripts/validation/fates_fortran_parity/fates_cnp_parity.jl \
#         [fortran_h0.nc | fates_cnp_fortran.csv] [nsteps=720] [out.csv]
#   (no argument => the committed fates_cnp_fortran.csv)
#   ... fates_cnp_parity.jl export <fortran_h0.nc> [out.csv]   # refresh the extract
using CLM, Printf, Dates, NCDatasets
include(joinpath(@__DIR__, "..", "..", "fates_fortran_parity.jl"))   # FSURDAT, FPARAM, FORCDIR, START, DTIME

const FATES_FIELDS = ("FATES_VEGN", "FATES_STOREN", "FATES_LEAFN", "FATES_FROOTN",
    "FATES_SAPWOODN", "FATES_REPRON", "FATES_STOREN_TF", "FATES_VEGP", "FATES_STOREP",
    "FATES_LEAFP", "FATES_FROOTP", "FATES_SAPWOODP", "FATES_REPROP", "FATES_STOREP_TF",
    "FATES_NH4UPTAKE", "FATES_NO3UPTAKE", "FATES_NDEMAND", "FATES_NEFFLUX",
    "FATES_NFIX_SYM", "FATES_PUPTAKE", "FATES_PDEMAND", "FATES_PEFFLUX",
    "FATES_EXCESS_RESP", "FATES_VEGC", "FATES_LEAFC", "FATES_STOREC", "FATES_FROOTC",
    "FATES_SAPWOODC", "FATES_REPROC", "FATES_NPP", "FATES_GPP", "FATES_NCOHORTS",
    "FATES_L2FR", "FATES_LAI", "FATES_TVEG24", "FATES_COLD_STATUS", "FATES_GDD")
# (SMIN_NH4/SMIN_NO3 are not scored: the port's soil-N summary never fills the
#  column NH4/NO3 totals — use_nitrif_denitrif is not passed to
#  soil_bgc_nitrogen_state_summary! — a host diagnostic gap outside CNP.)
const HOST_FIELDS = ("SMINN", "HR", "TOTLITC", "TOTSOMC")
const ALL_FIELDS = (FATES_FIELDS..., HOST_FIELDS...)

# FATES history name -> handle symbol, from the port's registry (same names as CTSM).
const NAME2IH = Dict(t[2] => t[10] for t in CLM._HISTORY_VAR_REGISTRY)

# FATES_CNP_PARTEH=1 runs the same harness carbon-only (the control for attributing
# a divergence to CNP vs the shared host/FATES path; N/P fields are then NaN).
const PARTEH = parse(Int, get(ENV, "FATES_CNP_PARTEH", "2"))

function build_cnp()
    inst, bounds, filt, tm = CLM.clm_initialize!(; fsurdat=FSURDAT, paramfile=FPARAM,
        fsnowoptics=FSNOWOPTICS, fsnowaging=FSNOWAGING,
        use_fates=true, fates_parteh_mode=PARTEH, start_date=START, dtime=Int(DTIME),
        use_bedrock=false)
    for g in 1:bounds.endg
        isfinite(inst.gridcell.latdeg[g]) || (inst.gridcell.latdeg[g] = inst.gridcell.lat[g]*180/π)
        isfinite(inst.gridcell.londeg[g]) || (inst.gridcell.londeg[g] = inst.gridcell.lon[g]*180/π)
    end
    for c in 1:bounds.endc   # as build(): flag the natural-veg soil column(s) as FATES
        inst.landunit.itype[inst.column.landunit[c]] == CLM.ISTSOIL && (inst.column.is_fates[c] = true)
    end
    config = CLM.CLMDriverConfig(use_fates=true)
    return inst, bounds, filt, CLM.clump_filter_inactive_and_active, config, tm
end

function _hist(inst, name)
    h = inst.fates.hist
    ih = get(NAME2IH, name, nothing)
    (h === nothing || ih === nothing || get(h.ih, ih, 0) <= 0) && return NaN
    return h.hvars[h.ih[ih]].r81d[1]
end

function cnp_record(inst, c)
    r = Dict{String,Float64}()
    for f in FATES_FIELDS
        r[f] = _hist(inst, f)
    end
    ns = inst.soilbiogeochem_nitrogenstate
    r["SMINN"]    = ns.sminn_col[c]
    r["HR"]       = inst.soilbiogeochem_carbonflux.hr_col[c]
    r["TOTLITC"]  = inst.soilbiogeochem_carbonstate.totlitc_col[c]
    r["TOTSOMC"]  = inst.soilbiogeochem_carbonstate.totsomc_col[c]
    return r
end

function run_julia(nsteps::Int)
    inst, bounds, filt, filt_ia, config, _tm = build_cnp()
    CLM.hlm_parteh_mode[] == PARTEH || error("PARTEH mode $(PARTEH) not active")
    c = findfirst(inst.column.is_fates)
    ng, nc = bounds.endg, bounds.endc
    fr = CLM.ForcingReader()
    CLM.forcing_reader_init!(fr, joinpath(FORCDIR, "clmforc.$(year(START)).nc"))
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
    recs = Vector{Dict{String,Float64}}()
    for i in 0:nsteps   # identical clock convention to fates_fortran_parity.jl
        cur = START + Second(round(Int, max(i - 1, 0) * DTIME))
        tod = hour(cur)*3600 + minute(cur)*60 + second(cur)
        calday = dayofyear(cur) + tod / CLM.SECSPDAY
        (declin, _) = CLM.compute_orbital(calday)
        (declinm1, _) = CLM.compute_orbital(calday - DTIME / CLM.SECSPDAY)
        CLM.init_daylength!(inst.gridcell, declin, declinm1, CLM.ORB_OBLIQR_DEFAULT, 1:ng)
        nextsw = calday + DTIME / CLM.SECSPDAY
        CLM.read_forcing_step!(fr, inst.atm2lnd, cur, ng, nc)
        CLM.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
        CLM.clm_drv!(config, inst, filt, filt_ia, bounds, CLM.nuopc_doalb(i), nextsw, declin, declin,
            CLM.ORB_OBLIQR_DEFAULT, false, false, "", false;
            nstep=i, is_first_step=(i <= 1),
            is_beg_curr_day=(i > 0 && tod == 0), is_end_curr_day=false, is_beg_curr_year=false,
            dtime=DTIME, mon=month(cur), day=day(cur), secs=tod,
            jday=dayofyear(cur), year=year(cur), photosyns=inst.photosyns)
        push!(recs, cnp_record(inst, c))
    end
    CLM.forcing_reader_close!(fr)
    return recs
end

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
        return Dict(String(h) => cols[i] for (i, h) in enumerate(hdr) if h != "nstep")
    end
    NCDataset(path, "r") do d
        Dict(f => Float64.(coalesce.(vec(Array(d[f][:])), NaN)) for f in ALL_FIELDS if haskey(d, f))
    end
end

"""
Per-step relative error over nstep >= 1, |J-F| / max(|F|, 1e-6 * peak|F|), and the
daily-sampled relative error at the end of each day (FATES daily fields only change
on the daily step). The h0 history is single precision, so ~1e-7 is the floor.
"""
function score(jl, ft; io=stdout)
    n = min(length(jl), minimum(length(v) for v in values(ft)))
    @printf(io, "  %-18s %13s %13s %10s %10s\n", "field", "Fortran(end)", "Julia(end)",
            "max rel", "mean rel")
    out = Dict{String,Tuple{Float64,Float64}}()
    for f in ALL_FIELDS
        haskey(ft, f) || continue
        F = ft[f][1:n]; J = [jl[k][f] for k in 1:n]
        peak = maximum(abs, F[2:n])   # record 1 = nstep 0 (history unset / NaN)
        floor_ = max(peak * 1e-6, 1e-30)
        rel = [abs(J[k] - F[k]) / max(abs(F[k]), floor_) for k in 2:n]
        peak == 0 && (rel = [abs(J[k]) > 0 ? Inf : 0.0 for k in 2:n])
        mx = maximum(rel); mn = sum(rel) / length(rel)
        out[f] = (mx, mn)
        @printf(io, "  %-18s %13.5e %13.5e %10.3e %10.3e\n", f, F[n], J[n], mx, mn)
    end
    return out
end

function export_fortran_csv(nc_path, csv_path)
    ft = read_fortran(nc_path)
    names = [f for f in ALL_FIELDS if haskey(ft, f)]
    n = length(ft[names[1]])
    open(csv_path, "w") do io
        println(io, join(["nstep"; names], ","))
        for k in 1:n
            println(io, join([string(k - 1); [repr(ft[f][k]) for f in names]], ","))
        end
    end
    return csv_path
end

function write_csv(path, jl, ft)
    open(path, "w") do io
        println(io, join(["nstep"; [f * "_jl" for f in ALL_FIELDS]; [f * "_ft" for f in ALL_FIELDS]], ","))
        n = min(length(jl), minimum(length(v) for v in values(ft)))
        for k in 1:n
            println(io, join([string(k - 1);
                              [string(jl[k][f]) for f in ALL_FIELDS];
                              [string(haskey(ft, f) ? ft[f][k] : NaN) for f in ALL_FIELDS]], ","))
        end
    end
end

const FORTRAN_CSV = joinpath(@__DIR__, "fates_cnp_fortran.csv")

function cnp_main(args)
    if !isempty(args) && args[1] == "export"
        export_fortran_csv(args[2], length(args) >= 3 ? args[3] : FORTRAN_CSV)
        return 0
    end
    ref = isempty(args) ? FORTRAN_CSV : args[1]
    nsteps = length(args) >= 2 ? parse(Int, args[2]) : 720
    ft = read_fortran(ref)
    jl = run_julia(nsteps)
    println("\n=== FATES CNP (parteh_mode=2): Julia vs Fortran (", nsteps, " steps) ===")
    score(jl, ft)
    length(args) >= 3 && write_csv(args[3], jl, ft)
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(cnp_main(ARGS))
end
