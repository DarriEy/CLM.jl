# compare_h0_bitwise.jl — are two CTSM history files bit-identical in every data variable?
#   julia --project=. compare_h0_bitwise.jl a.nc b.nc
using NCDatasets
a, b = NCDataset(ARGS[1]), NCDataset(ARGS[2])
skip = Set(["date_written", "time_written"])   # wall-clock stamps
ndiff = 0; nvar = 0
for v in keys(a)
    v in skip && continue
    haskey(b, v) || (println("missing in b: ", v); global ndiff += 1; continue)
    x = Array(a[v].var[:]); y = Array(b[v].var[:])
    global nvar += 1
    if !isequal(x, y)
        global ndiff += 1
        println("DIFF ", v)
    end
end
println(nvar, " variables compared, ", ndiff, " differ")
exit(ndiff == 0 ? 0 : 1)
