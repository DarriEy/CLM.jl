#!/usr/bin/env julia
# ==========================================================================
# fates_multiyear_compare.jl — score a multi-year Julia FATES run of the
# fates_parity_1pt case (scripts/fates_parity_multiyear.jl CSV) against the
# Fortran CTSM-FATES daily-history reference
# (scripts/validation/fates_fortran_parity/fates_multiyear_fortran_daily.csv.gz).
#
#   julia --project=. scripts/fates_multiyear_compare.jl <julia.csv> [fortran.csv(.gz)]
#
# Row alignment: the Julia row for model day d is written after day d's last step;
# the Fortran h0 record d+1 is the average over day d+1 (the history time stamp
# lags one day), so Julia row d is compared with Fortran row d+1.
#
# The Fortran history has no per-cohort dbh, only FATES_NPLANT_SZPF on the
# fates_history_sizeclass_bin_edges, so "max dbh" is compared as the lower edge of
# the largest occupied size class, per PFT.
# ==========================================================================
using Printf

const FMY_SZ_EDGES = (0, 5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100)
const FMY_SITE_AREA = 1.0e4      # FATES site area [m2]; Julia CSV is per site
const FMY_WOODY = 1:11           # fates_woody = 1 (default parameter file)
const FMY_GRASS = 12:14

_fmy_lines(path) = endswith(path, ".gz") ? readlines(pipeline(`gzip -dc $path`)) : readlines(path)

function fmy_read_csv(path)
    ls = _fmy_lines(path)
    hdr = split(ls[1], ',')
    rows = [split(l, ',') for l in ls[2:end] if !isempty(l)]
    col(name) = (i = findfirst(==(name), hdr); i === nothing && error("no column $name in $path");
                 [parse(Float64, r[i]) for r in rows])
    return col, length(rows)
end

fmy_szbin(dbh) = dbh <= 0 ? -1 : FMY_SZ_EDGES[findlast(e -> dbh > e, FMY_SZ_EDGES)]

"""
    multiyear_summary(julia_csv, fortran_csv; ndays_year=365)

Per-model-year summary of both runs. Returns a vector of NamedTuples, one per
complete model year present in BOTH files.
"""
function multiyear_summary(julia_csv::AbstractString, fortran_csv::AbstractString;
                           ndays_year::Int = 365)
    J, nj = fmy_read_csv(julia_csv)
    F, nf = fmy_read_csv(fortran_csv)
    nd = min(nj, nf - 1)
    jrec = J("recruit_n") ./ FMY_SITE_AREA; frec = F("recruit_d")
    jnpp = J("npp_d") ./ FMY_SITE_AREA ./ 365 ; fnpp = F("npp_d")   # gC/m2/day
    jcoh = J("ncoh"); fcoh = F("ncoh")
    jmax = [J("maxdbh_$p") for p in 1:14]; fmax = [F("maxsz_$p") for p in 1:14]
    out = NamedTuple[]
    for y in 1:(nd ÷ ndays_year)
        dj = (y-1)*ndays_year+1 : y*ndays_year      # Julia rows
        df = dj .+ 1                                  # Fortran rows (one-day lag)
        jw = maximum(maximum(jmax[p][dj]) for p in FMY_WOODY)
        jw0 = maximum(jmax[p][first(dj)] for p in FMY_WOODY)
        push!(out, (year = y,
            rec_j = sum(jrec[dj]), rec_f = sum(frec[df]),
            npp_j = sum(jnpp[dj]) / ndays_year, npp_f = sum(fnpp[df]) / ndays_year,
            ncoh_j = maximum(jcoh[dj]), ncoh_f = maximum(fcoh[df]),
            woody_maxdbh_j = jw, woody_ddbh_j = jw - jw0,
            woody_maxbin_j = fmy_szbin(jw), woody_maxbin_f = maximum(maximum(fmax[p][df]) for p in FMY_WOODY),
            grass_maxbin_j = maximum(fmy_szbin(maximum(jmax[p][dj])) for p in FMY_GRASS),
            grass_maxbin_f = maximum(maximum(fmax[p][df]) for p in FMY_GRASS),
            grass_maxdbh_j = maximum(maximum(jmax[p][dj]) for p in FMY_GRASS)))
    end
    return out
end

function print_summary(s)
    println(" yr | recruit /m2/yr  J / F   ratio | NPP gC/m2/d J / F | max ncoh J / F | woody maxdbh J [bin] / F bin | grass maxdbh J [bin] / F bin")
    for r in s
        @printf(" %2d | %9.3g / %9.3g  %5.2f | %6.3f / %6.3f | %5d / %5d | %7.2f [%3d] / %3d | %8.2f [%3d] / %3d\n",
                r.year, r.rec_j, r.rec_f, r.rec_j / max(r.rec_f, eps()), r.npp_j, r.npp_f,
                Int(r.ncoh_j), Int(r.ncoh_f), r.woody_maxdbh_j, r.woody_maxbin_j,
                Int(r.woody_maxbin_f), r.grass_maxdbh_j, r.grass_maxbin_j, Int(r.grass_maxbin_f))
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    jcsv = ARGS[1]
    fcsv = length(ARGS) >= 2 ? ARGS[2] :
        joinpath(@__DIR__, "validation", "fates_fortran_parity", "fates_multiyear_fortran_daily.csv.gz")
    print_summary(multiyear_summary(jcsv, fcsv))
end
