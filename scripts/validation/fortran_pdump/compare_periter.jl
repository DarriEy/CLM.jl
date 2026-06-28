# compare_periter.jl — align the Fortran and Julia per-iteration leaf-temperature
# Newton dumps by (itlef, p) and report, per row, the per-column relative
# difference and the FIRST column (in CanopyFluxes dependency order) that
# diverges beyond a tolerance. Localizes the seed of the residual T_VEG
# divergence at nstep 13461.
#
#   julia compare_periter.jl [fortran.txt] [julia.txt] [reltol]

using Printf

const FDEF = joinpath(@__DIR__, "periter_n13461_fortran.txt")
const JDEF = joinpath(@__DIR__, "periter_n13461_julia.txt")
ffile = length(ARGS) >= 1 ? ARGS[1] : FDEF
jfile = length(ARGS) >= 2 ? ARGS[2] : JDEF
reltol = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 1e-12

function readdump(path)
    lines = readlines(path)
    hdr = split(strip(lines[1]))
    rows = Dict{Tuple{Int,Int},Vector{Float64}}()
    for ln in lines[2:end]
        isempty(strip(ln)) && continue
        f = split(strip(ln))
        vals = parse.(Float64, f)
        itlef = Int(vals[2]); p = Int(vals[3])
        rows[(itlef, p)] = vals
    end
    return hdr, rows
end

hdr, F = readdump(ffile)
_, J   = readdump(jfile)

# columns in dependency order (input-side first, output-side last):
# friction velocity (ustar) -> end-of-iter um/obu (downstream) ; photosynthesis
# (rssun/rssha) ; flux coefficients (wtl0/wtg0/wtga/wtgaq, qsatl/qsatldT, dc1/dc2)
# ; fluxes (efe/efsh) ; Newton update (dt_veg/t_veg/del/del2). Statics last.
const COLNAMES = hdr  # index 1=nstep 2=itlef 3=p, 4..=fields
relerr(a,b) = (a==b) ? 0.0 : abs(a-b)/max(abs(a),abs(b),1e-300)

keys_sorted = sort(collect(intersect(keys(F), keys(J))); by = k->(k[2],k[1]))

println("="^110)
println("Per-iteration Fortran vs Julia  (reltol=$reltol)  — columns flagged DIFF when |rel| > reltol")
println("="^110)

# Per-patch: find the first iteration with any DIFF, and which columns.
for patch in sort(unique(k[2] for k in keys_sorted))
    println("\n########## PATCH $patch ##########")
    ks = sort([k for k in keys_sorted if k[2]==patch]; by=k->k[1])
    first_div_reported = false
    for k in ks
        fv = F[k]; jv = J[k]
        diffs = String[]
        maxrel = 0.0; maxcol = ""
        for ci in 4:length(COLNAMES)
            r = relerr(fv[ci], jv[ci])
            if r > reltol
                push!(diffs, @sprintf("%s(rel=%.2e: F=%.6g J=%.6g)", COLNAMES[ci], r, fv[ci], jv[ci]))
            end
            if r > maxrel; maxrel = r; maxcol = COLNAMES[ci]; end
        end
        @printf("itlef=%d  maxrel=%.3e @ %-8s  ndiff=%d\n", k[1], maxrel, maxcol, length(diffs))
        if !first_div_reported && !isempty(diffs)
            println("   FIRST DIVERGING ITERATION for patch $patch — diverging columns (any order):")
            for d in diffs; println("      ", d); end
            first_div_reported = true
        end
    end
end

# Iteration-1 dependency-order table for patches 2 & 3 (the seed view).
println("\n", "="^110)
println("ITERATION 1 — dependency-order column comparison (patches 2 & 3)")
println("="^110)
order = ["ustar","um","obu","rssun","rssha","qsatl","qsatldT","wtl0","wtg0","wtga","wtgaq",
         "dc1","dc2","efsh","efe","dt_veg","del","del2","t_veg","tlbef","sabv","air","bir","cir"]
colidx = Dict(COLNAMES[i]=>i for i in 1:length(COLNAMES))
for patch in (2,3)
    haskey(F,(1,patch)) || continue
    println("\n-- patch $patch, itlef=1 --")
    @printf("  %-9s %22s %22s %12s\n","column","Fortran","Julia","rel")
    for nm in order
        ci = colidx[nm]
        f = F[(1,patch)][ci]; j = J[(1,patch)][ci]
        @printf("  %-9s %22.12g %22.12g %12.3e\n", nm, f, j, relerr(f,j))
    end
end
