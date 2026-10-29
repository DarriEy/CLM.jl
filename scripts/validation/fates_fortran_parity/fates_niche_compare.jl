# fates_niche_compare.jl — score a Julia FATES research-mode run
# (scripts/fates_niche_modes.jl CSV) against the Fortran daily history of the same mode
# (stage_niche_mode.py), and — the part that isolates the mode — compare the mode's
# EFFECT (mode run minus the all-off base run) between the two models.
#
#   julia --project=. scripts/validation/fates_fortran_parity/fates_niche_compare.jl \
#       <julia_mode.csv> <fortran_mode_h0.nc> [<julia_base.csv> <fortran_base_h0.nc>]
#
# Row alignment: the Julia CSV row `day = d` is the state after d model days (row 0 =
# cold start). The Fortran h0 record t (1-based) carries the FATES dynamics fields as
# of the start of day t, i.e. after d = t-1 days. So the i-th Julia row is compared
# with the i-th Fortran record (checked on the integer cohort/patch counts).
using NCDatasets, Printf

function read_julia(path)
    ls = readlines(path)
    hdr = split(ls[1], ',')
    rows = [parse.(Float64, split(l, ',')) for l in ls[2:end] if !isempty(l)]
    return Dict(String(h) => [r[i] for r in rows] for (i, h) in enumerate(hdr))
end

"The Fortran history in the Julia CSV's keys (per m2, kgC)."
function read_fortran(path)
    ds = NCDataset(path)
    v(n) = haskey(ds, n) ? Float64.(coalesce.(Array(ds[n]), NaN)) : nothing
    out = Dict{String,Vector{Float64}}()
    flat(a) = ndims(a) == 2 ? vec(a[1, :]) : a       # (lndgrid, time) -> time
    for (k, n) in (("ncoh", "FATES_NCOHORTS"), ("npatch", "FATES_NPATCHES"),
                   ("vegc", "FATES_VEGC"), ("leafc", "FATES_LEAFC"),
                   ("storec", "FATES_STOREC"), ("structc", "FATES_STRUCTC"),
                   ("seed_bank", "FATES_SEED_BANK"),
                   # host daily means (FATES writes no history at all under ST3)
                   ("host_fsh", "FSH"), ("host_tv", "TV"), ("host_hr", "HR"))
        a = v(n); a === nothing || (out[k] = flat(a))
    end
    a = v("FATES_NPLANT_PF")                          # (lndgrid, pft, time)
    if a !== nothing
        for f in 1:size(a, 2); out["nplant_$f"] = vec(a[1, f, :]); end
    end
    a = v("FATES_SEEDS_OUT_GRIDCELL_PF")
    if a !== nothing
        for f in 1:size(a, 2); out["seed_out_$f"] = vec(a[1, f, :]); end
    end
    a = v("FATES_NPLANT_CDPF")
    if a !== nothing
        cd = ds["fates_cdmap_levcdpf"][:]; pf = ds["fates_pftmap_levcdpf"][:]
        for i in eachindex(cd)
            k = "nplant_cd$(cd[i])_$(pf[i])"
            out[k] = get(out, k, zeros(size(a, 3))) .+ vec(a[1, i, :])
        end
    end
    a = v("FATES_NPLANT_ACPF")
    if a !== nothing
        ca = ds["fates_camap_levcapf"][:]; pf = ds["fates_pftmap_levcapf"][:]
        for i in eachindex(ca)
            k = "nplant_ca$(ca[i])_$(pf[i])"
            out[k] = get(out, k, zeros(size(a, 3))) .+ vec(a[1, i, :])
        end
    end
    close(ds)
    return out
end

relerr(a, b) = abs(a - b) / max(abs(b), 1e-30)
_fmax(v) = (w = filter(isfinite, v); isempty(w) ? NaN : maximum(w))   # NaN-safe max

"""
    niche_compare(jcsv, fnc; jbase=nothing, fbase=nothing, keys=nothing, io=stdout)

Per-key worst relative difference over the run, Julia vs Fortran, and (with the base
runs) the mode effect `mode - base` in each model. Returns a Dict of results.
"""
function niche_compare(jcsv, fnc; jbase = nothing, fbase = nothing, keys = nothing,
                       io = stdout)
    J = read_julia(jcsv); F = read_fortran(fnc)
    JB = jbase === nothing ? nothing : read_julia(jbase)
    FB = fbase === nothing ? nothing : read_fortran(fbase)
    ks = keys === nothing ? sort!([k for k in Base.keys(F) if haskey(J, k)]) : keys
    nd = min(length(J["day"]), length(F[first(ks)]))
    res = Dict{String,Any}()
    @printf(io, "  %-16s %14s %14s %9s", "key", "julia(end)", "fortran(end)", "maxrel")
    JB === nothing || @printf(io, " %14s %14s %9s", "Δjulia(end)", "Δfortran(end)", "Δrel")
    println(io)
    for k in ks
        haskey(J, k) && haskey(F, k) || continue
        jv = J[k][1:nd]; fv = F[k][1:nd]
        all(iszero, fv) && all(iszero, jv) && continue
        mr = _fmax(relerr.(jv, fv))
        line = @sprintf("  %-16s %14.7g %14.7g %9.2e", k, jv[end], fv[end], mr)
        r = (; maxrel = mr, j = jv[end], f = fv[end])
        if JB !== nothing && haskey(JB, k) && haskey(FB, k)
            dj = jv .- JB[k][1:nd]; df = fv .- FB[k][1:nd]
            drel = _fmax(abs.(dj .- df)) / max(_fmax(abs.(df)), 1e-30)
            line *= @sprintf(" %14.6g %14.6g %9.2e", dj[end], df[end], drel)
            r = (; r..., dj = dj[end], df = df[end], drel)
        end
        println(io, line)
        res[k] = r
    end
    return res
end

if abspath(PROGRAM_FILE) == @__FILE__
    a = ARGS
    niche_compare(a[1], a[2]; jbase = length(a) >= 3 ? a[3] : nothing,
                  fbase = length(a) >= 4 ? a[4] : nothing)
end
