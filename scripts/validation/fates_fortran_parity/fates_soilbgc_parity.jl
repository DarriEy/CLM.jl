# fates_soilbgc_parity.jl — FATES -> CLM soil-biogeochemistry coupling vs Fortran CTSM.
#
# Steps the Julia port on the fates_parity_1pt site (Bow at Banff, 2002-07-01, 3600 s,
# FATES non-SP => use_fates_bgc, CENTURYKoven2013) and scores the soil-BGC diagnostics
# step by step against an hourly h0 history of the same single-point Fortran case:
#
#   HR, LITTERC_HR, SOILC_HR                 heterotrophic respiration      gC/m2/s
#   LIT_MET_C, LIT_CEL_C, LIT_LIG_C          litter pools (FATES input)     gC/m2
#   SOM_ACT_C, SOM_SLO_C, SOM_PAS_C          soil organic matter pools      gC/m2
#   TOTLITC, TOTSOMC                         column totals                  gC/m2
#   FATES_HET_RESP                           FATES history (= HR)           kg/m2/s
#   FATES_LITTER_OUT                         FATES litter flux to the soil  kg/m2/s
#   FLIT (cf.fates_litter_flux, Julia only)  host-side litter input         gC/m2/s
#
# Fortran reference (see README.md, "SOIL-BGC COUPLING"): stage_soilbgc_run.sh runs the
# case with hist_nhtfrq=1 output for these fields; the scored series are committed as
# fates_soilbgc_fortran.csv (240 steps) so the comparison reproduces without CTSM.
#
# usage:
#   julia --project=. scripts/validation/fates_fortran_parity/fates_soilbgc_parity.jl \
#         [fortran_h0.nc | fates_soilbgc_fortran.csv] [nsteps=240] [out.csv]
#   (no argument => the committed fates_soilbgc_fortran.csv, 240 steps)
#   ... fates_soilbgc_parity.jl export <fortran_h0.nc> [out.csv]   # refresh the extract
using CLM, Printf, Dates, NCDatasets
include(joinpath(@__DIR__, "..", "..", "fates_fortran_parity.jl"))   # build(), FORCDIR, START, DTIME

const POOLS = ("LIT_MET_C", "LIT_CEL_C", "LIT_LIG_C", "SOM_ACT_C", "SOM_SLO_C", "SOM_PAS_C")
const FIELDS = ("HR", "LITTERC_HR", "SOILC_HR", POOLS..., "TOTLITC", "TOTSOMC",
                "FATES_HET_RESP", "FATES_LITTER_OUT")
# Host drivers of decomposition (attribution): layer-1..3 temperature/moisture
# scalars and soil temperature. Divergence here is a host soil-physics residual.
const LAYER_FIELDS = ("T_SCALAR", "W_SCALAR", "TSOI")
const NLAYER_DIAG = 3
const DIAG_FIELDS = Tuple(f * "_" * string(j) for f in LAYER_FIELDS for j in 1:NLAYER_DIAG)
const ALL_FIELDS = (FIELDS..., DIAG_FIELDS...)

function _fates_hist(inst, ih::Symbol)
    h = inst.fates.hist
    h === nothing && return NaN
    return h.hvars[h.ih[ih]].r81d[1]
end

"Per-step Julia soil-BGC record for the (first) FATES column."
function soilbgc_record(inst, c)
    cf = inst.soilbiogeochem_carbonflux
    cs = inst.soilbiogeochem_carbonstate
    r = Dict{String,Float64}()
    r["HR"]         = cf.hr_col[c]
    r["LITTERC_HR"] = cf.lithr_col[c]
    r["SOILC_HR"]   = cf.somhr_col[c]
    for (k, name) in enumerate(POOLS)
        r[name] = cs.decomp_cpools_col[c, k]
    end
    r["TOTLITC"] = cs.totlitc_col[c]
    r["TOTSOMC"] = cs.totsomc_col[c]
    r["FATES_HET_RESP"]   = _fates_hist(inst, :ih_hr_si)
    r["FATES_LITTER_OUT"] = _fates_hist(inst, :ih_litter_out_si)
    r["FLIT"] = length(cf.fates_litter_flux) >= c ? cf.fates_litter_flux[c] : NaN
    for j in 1:NLAYER_DIAG
        r["T_SCALAR_$j"] = cf.t_scalar_col[c, j]
        r["W_SCALAR_$j"] = cf.w_scalar_col[c, j]
        r["TSOI_$j"] = inst.temperature.t_soisno_col[c, CLM.varpar.nlevsno + j]
    end
    return r
end

"Step the port nsteps (plus the NUOPC nstep=0 init step), recording after every step."
function run_julia(nsteps::Int)
    inst, bounds, filt, filt_ia, config, _tm = build()
    config.use_fates_bgc || error("expected a use_fates_bgc configuration")
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
            jday=dayofyear(cur), photosyns=inst.photosyns)
        push!(recs, soilbgc_record(inst, c))
    end
    CLM.forcing_reader_close!(fr)
    return recs
end

"""
    read_fortran(path) -> Dict{String,Vector{Float64}}

Fortran series (record k <-> nstep k-1) from the hourly h0 history (`.nc`) or from
the committed extract `fates_soilbgc_fortran.csv` (see `export_fortran_csv`).
"""
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
        out = Dict(f => Float64.(vec(Array(d[f][:]))) for f in FIELDS if haskey(d, f))
        for f in LAYER_FIELDS
            haskey(d, f) || continue
            a = Array(d[f][:, :, :])           # (lndgrid, lev, time)
            for j in 1:NLAYER_DIAG
                out[f * "_" * string(j)] = Float64.(vec(a[1, j, :]))
            end
        end
        out
    end
end

"""
    score(jl, ft) -> NamedTuple per field

Per-step relative error |J-F| / max(|F|, floor) over nstep >= 1, where the floor
is 1e-6 x the field's peak |F| over the window (so near-zero night/cold-start values do not
dominate). Reports max and mean relative error, plus the end-of-window values.
"""
function score(jl, ft; io=stdout)
    n = min(length(jl), minimum(length(v) for v in values(ft)))
    @printf(io, "  %-18s %14s %14s %10s %10s\n", "field", "Fortran(end)", "Julia(end)",
            "max rel", "mean rel")
    out = Dict{String,Tuple{Float64,Float64}}()
    for f in ALL_FIELDS
        haskey(ft, f) || continue
        F = ft[f][1:n]; J = [jl[k][f] for k in 1:n]
        peak = maximum(abs, F)
        floor_ = max(peak * 1e-6, 1e-30)
        # nstep 0 (record 1) is the NUOPC init step, before FATES history is filled.
        rel = [abs(J[k] - F[k]) / max(abs(F[k]), floor_) for k in 2:n]
        mx = maximum(rel); mn = sum(rel) / length(rel)
        out[f] = (mx, mn)
        @printf(io, "  %-18s %14.6e %14.6e %10.3e %10.3e\n", f, F[n], J[n], mx, mn)
    end
    return out
end

"Extract the scored Fortran series from an h0 file into a small committed CSV."
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
        cols = ["nstep"; [f * "_jl" for f in ALL_FIELDS]; [f * "_ft" for f in ALL_FIELDS]; "FLIT_jl"]
        println(io, join(cols, ","))
        n = min(length(jl), minimum(length(v) for v in values(ft)))
        for k in 1:n
            vals = [string(k - 1);
                    [string(jl[k][f]) for f in ALL_FIELDS];
                    [string(haskey(ft, f) ? ft[f][k] : NaN) for f in ALL_FIELDS];
                    string(jl[k]["FLIT"])]
            println(io, join(vals, ","))
        end
    end
end

const FORTRAN_CSV = joinpath(@__DIR__, "fates_soilbgc_fortran.csv")

function soilbgc_main(args)
    if !isempty(args) && args[1] == "export"
        export_fortran_csv(args[2], length(args) >= 3 ? args[3] : FORTRAN_CSV)
        return 0
    end
    ref = isempty(args) ? FORTRAN_CSV : args[1]
    nsteps = length(args) >= 2 ? parse(Int, args[2]) : 240
    ft = read_fortran(ref)
    jl = run_julia(nsteps)
    println("\n=== FATES soil-BGC coupling: Julia vs Fortran (", nsteps, " steps) ===")
    score(jl, ft)
    length(args) >= 3 && write_csv(args[3], jl, ft)
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(soilbgc_main(ARGS))
end
