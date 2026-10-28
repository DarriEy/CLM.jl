# Restart injection for the CNSoyfix reference: zero the soil MINERAL-N pools
# (sminn_vr, smin_nh4_vr, smin_no3_vr) on the soybean columns of a non-FUN crop
# restart, leaving every other state (organic pools, vegetation, phenology, GDD
# accumulators) untouched.
#
# Why: CNSoyfix only fixes N when fpg < 1, and this site is N-replete by three
# orders of magnitude (docs/CROP_LIFECYCLE_MAP.md §3e-§3g; also verified on the
# 2026-10-28 cold run: FPG == 1 on every crop column every day of 2016-2018, even
# in year 1 of a cold start with ndep = 0). With the mineral pools emptied, the
# first step's SoilBiogeochemCompetition finds no supply -> fpg ~ 0, so the NEXT
# step's CNSoyfix (which reads the previous step's fpg) runs its fpg < 1 branch,
# and fpg then recovers over the following steps as mineralization refills the
# pool -- a sweep of fpg values through (0, 1) at a fixed growth stage.
# The numbers produced are still the Fortran routine's own output; only the
# initial condition is constructed.
#
# Mineral-N zeroing alone is NOT enough here (verified 2026-10-28): the column
# mineralizes ~1.3e-3 gN/m2 per hourly step while soybean demand is ~5e-4, so the
# pool refills within one step and fpg stays 1. With ORGANIC=1 the decomposing
# litter/SOM/CWD C and N pools of those columns are zeroed too, which removes the
# mineralization source; supply is then only free-living fixation + the soybean's
# own fixation, and fpg < 1 is reached.
#
# Usage: [ORGANIC=1] julia --project=. make_nlimited_restart.jl <in restart> <out restart> [itypes...]
using NCDatasets

src, dst = ARGS[1], ARGS[2]
itypes = length(ARGS) > 2 ? parse.(Int, ARGS[3:end]) : [23, 24]   # temperate soybean (+irrigated)
vars = ["sminn_vr", "smin_nh4_vr", "smin_no3_vr"]
if get(ENV, "ORGANIC", "0") == "1"
    for pool in ("litr1", "litr2", "litr3", "soil1", "soil2", "soil3", "cwd"), e in ("c", "n")
        push!(vars, pool * e * "_vr")
    end
end
cp(src, dst; force = true)
chmod(dst, 0o644)
NCDataset(dst, "a") do ds
    it = ds["pfts1d_itypveg"][:]
    ci = ds["pfts1d_column_index"][:]
    wt = ds["pfts1d_wtcol"][:]
    cols = sort(unique([ci[p] for p in eachindex(it) if it[p] in itypes && wt[p] > 0]))
    println("zeroing mineral N on columns $cols (patch itypes $itypes)")
    for v in vars
        haskey(ds, v) || (println("  (no $v in restart)"); continue)
        a = Array(ds[v][:, :])
        dn = dimnames(ds[v])
        ccol = findfirst(==("column"), dn)
        for c in cols
            if ccol == 1
                println("  $v col $c: sum before = $(sum(a[c, :]))"); a[c, :] .= 0
            else
                println("  $v col $c: sum before = $(sum(a[:, c]))"); a[:, c] .= 0
            end
        end
        ds[v][:, :] = a
    end
    ds.attrib["crop_tails_note"] = "CLM.jl CNSoyfix reference: $(join(vars, ",")) zeroed on columns $cols"
end
println("wrote $dst")
