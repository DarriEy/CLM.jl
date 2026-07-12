#!/usr/bin/env julia

using NCDatasets
using Printf

const FDIR = get(ENV, "FORTRAN_SNOW_DUMPDIR", "/tmp/hj_pdump_snow")
const JDIR = get(ENV, "JULIA_SNOW_DUMPDIR", "/tmp/hj_julia_snow")
const LO = parse(Int, get(ENV, "SNOW_DUMP_LO", "9900"))
const HI = parse(Int, get(ENV, "SNOW_DUMP_HI", "10960"))

const FIELDS = [
    "SNLSNO", "SNOW_DEPTH", "frac_sno", "frac_sno_eff", "SNOMELT_ACCUM",
    "H2OSNO_NO_LAYERS", "INT_SNOW", "T_GRND", "H2OSFC",
    "QFLX_SNOMELT", "EFLX_SNOMELT", "QFLX_SNOFRZ", "QFLX_SNOW_DRAIN",
    "QFLX_RAIN_PLUS_SNOMELT", "QFLX_LIQ_GRND", "QFLX_SNOW_GRND",
    "QFLX_SOLIDDEW_TOP", "QFLX_LIQDEW_TOP", "QFLX_SOLIDEVAP_TOP",
    "QFLX_LIQEVAP_TOP", "T_SOISNO", "H2OSOI_LIQ", "H2OSOI_ICE",
    "DZSNO", "QFLX_SNOMELT_LYR", "QFLX_SNOFRZ_LYR", "QFLX_SNOW_PERC",
]

tofloat(x) = ismissing(x) ? NaN : Float64(x)
arr(ds, name) = map(tofloat, ds[name][:])

function maxdiff(a, b)
    best_abs = -Inf
    best_rel = -Inf
    best_i = 0
    best_j = NaN
    best_f = NaN
    for i in eachindex(a, b)
        jv = a[i]
        fv = b[i]
        isnan(jv) && isnan(fv) && continue
        ad = abs(jv - fv)
        rd = ad / (1.0 + max(abs(jv), abs(fv)))
        if ad > best_abs
            best_abs = ad
            best_rel = rd
            best_i = Int(i)
            best_j = jv
            best_f = fv
        end
    end
    return best_abs, best_rel, best_i, best_j, best_f
end

mutable struct Stat
    maxabs::Float64
    maxrel::Float64
    step::Int
    idx::Int
    jv::Float64
    fv::Float64
    first_step::Int
    first_abs::Float64
    first_rel::Float64
end

stats = Dict(name => Stat(-Inf, -Inf, 0, 0, NaN, NaN, 0, NaN, NaN) for name in FIELDS)

for n in LO:HI
    jf = joinpath(JDIR, "julia_snow_step_n$(n).nc")
    ff = joinpath(FDIR, "pdump_after_hydrologydrainage_n$(n).nc")
    isfile(jf) || error("missing Julia dump: $jf")
    isfile(ff) || error("missing Fortran dump: $ff")
    jds = NCDataset(jf)
    fds = NCDataset(ff)
    for name in FIELDS
        haskey(jds, name) || continue
        haskey(fds, name) || continue
        ad, rd, idx, jv, fv = maxdiff(arr(jds, name), arr(fds, name))
        st = stats[name]
        if st.first_step == 0 && (ad > 1e-9 && rd > 1e-12)
            st.first_step = n
            st.first_abs = ad
            st.first_rel = rd
        end
        if ad > st.maxabs
            st.maxabs = ad
            st.maxrel = rd
            st.step = n
            st.idx = idx
            st.jv = jv
            st.fv = fv
        end
    end
    close(jds)
    close(fds)
end

rows = sort(collect(stats), by = kv -> kv[2].maxabs, rev = true)

println("="^92)
println("Maritime snow per-step dump diff: Julia vs Fortran after_hydrologydrainage")
println("Fortran: $FDIR")
println("Julia  : $JDIR")
println("nstep  : $LO:$HI")
println("="^92)
@printf("%-24s %7s %5s %12s %12s %14s %14s %12s\n",
        "field", "step", "idx", "max_abs", "max_rel", "Julia", "Fortran", "first_step")
println("-"^92)
for (name, st) in rows
    isfinite(st.maxabs) || continue
    @printf("%-24s %7d %5d %12.4e %12.4e %14.6e %14.6e %12d\n",
            name, st.step, st.idx, st.maxabs, st.maxrel, st.jv, st.fv, st.first_step)
end
