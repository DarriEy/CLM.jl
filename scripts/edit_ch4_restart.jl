# Scale a CLM restart's saturated/unsaturated CH4 concentration by a factor, to
# construct a CH4-EBULLITION Fortran reference (docs/CH4_FIRE_PARITY.md §14). The
# saturated dissolved-CH4 concentration is a prognostic state; scaling it in a copy
# of a warm restart pushes vgc across the bubble threshold (vgc>0.0855) so both
# CLM.jl and CTSM compute ebullition from the identical (scaled) before-step conc —
# the synthetic-shared-input technique of the peatf (§12) / lnfm (§2) references.
#
# Usage: julia +1.12 --project=. scripts/edit_ch4_restart.jl <raw_restart.nc> <out.nc> <factor>
using NCDatasets
raw = ARGS[1]; out = ARGS[2]; F = parse(Float64, ARGS[3])
cp(raw, out; force=true)
NCDataset(out, "a") do ds
    for v in ("CONC_CH4_SAT", "CONC_CH4_UNSAT")
        haskey(ds, v) || continue
        a = ds[v][:, :]
        b = copy(a)
        for i in eachindex(b)
            ismissing(b[i]) && continue
            b[i] = b[i] * F
        end
        ds[v][:, :] = b
        mx = maximum(x -> ismissing(x) ? -Inf : x, b)
        println("scaled $v by $F ; new max = $mx")
    end
end
println("wrote $out")
