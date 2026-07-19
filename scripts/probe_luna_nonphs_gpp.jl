# probe_luna_nonphs_gpp.jl — does the LUNA branch in the NON-PHS photosynthesis!
# actually reach GPP?
#
# The bug this probes (#267): `photosynthesis!` ACCEPTED use_luna and never read
# it, because the LUNA branch had only been ported into the PHS routine. LUNA
# acclimated vcmx25_z every day while vcmax_z / psnsun / fpsn stayed bit-identical
# to a LUNA-off run. Every wiring assertion along the chain passed.
#
# So this runs Bow TWICE on the default (non-PHS, use_hydrstress=false) path —
# use_luna=true and use_luna=false — and reports, per day:
#   * vcmx25_z   : does LUNA acclimate at all?  (it always did)
#   * fpsn/psnsun: does that acclimation reach the carbon flux? (it did NOT)
# LUNA needs ~240 h of accumulator fill before it moves, so a run shorter than
# ~10 days shows no difference for entirely innocent reasons.
#
#   julia +1.12 --project=. scripts/probe_luna_nonphs_gpp.jl
#
# Output: scratch CSV + a summary table on stdout.

using Dates, Printf, CLM

include(joinpath(@__DIR__, "fortran_parity_common.jl"))

const START = DateTime(2003, 6, 1)
const STOP  = DateTime(2003, 6, 17)          # 16 days: ~6 days past accumulator fill
const OUT   = get(ENV, "LUNA_PROBE_OUT", joinpath(tempdir(), "probe_luna_nonphs_gpp.csv"))

# Per-day, per-patch snapshot of the LUNA state and the carbon flux it should drive.
function make_probe(rows, tag)
    lastday = Ref(-1)
    return function (inst, tm)
        d = tm.current_date
        Dates.day(d) == lastday[] && return
        lastday[] = Dates.day(d)
        ps = inst.photosyns
        for p in eachindex(inst.patch.itype)
            inst.patch.itype[p] == 0 && continue
            inst.patch.wtgcell[p] > 0.0 || continue
            # vcmx25_z is EMPTY on the LUNA-off run (it is only allocated under
            # use_luna) — record NaN rather than dropping the sample, or the two
            # runs share no keys and the comparison silently comes out empty.
            push!(rows, (tag = tag, date = d, p = p, itype = Int(inst.patch.itype[p]),
                vcmx25 = size(ps.vcmx25_z_patch, 1) >= p ?
                         Float64(ps.vcmx25_z_patch[p, 1]) : NaN,
                vcmax_z = Float64(ps.vcmax_z_patch[p, 1]),
                psnsun = Float64(ps.psnsun_patch[p]),
                psnsha = Float64(ps.psnsha_patch[p]),
                fpsn = Float64(ps.fpsn_patch[p])))
        end
    end
end

function run_case(tag::String, use_luna::Bool, rows)
    CLM.clm_run!(; fsurdat = FSURDAT, paramfile = FPARAM, fforcing = FFORCING,
        fhistory = joinpath(tempdir(), "probe_luna_nonphs_$(tag).nc"),
        start_date = START, end_date = STOP, dtime = 3600,
        use_cn = false, use_luna = use_luna,
        use_hydrstress = false,          # <-- the NON-PHS routine: the ported branch
        h2osfcflag = 1, verbose = false,
        baseflow_scalar = BASEFLOW_SCALAR, int_snow_max = INT_SNOW_MAX,
        fsnowoptics = FSNOWOPT, fsnowaging = FSNOWAGE,
        interp_forcing = true, step_probe = make_probe(rows, tag))
end

rows = NamedTuple[]
@printf("Bow, non-PHS path, %s .. %s (use_hydrstress=false)\n\n", START, STOP)
run_case("lunaoff", false, rows)
run_case("lunaon",  true,  rows)

open(OUT, "w") do io
    println(io, "tag,date,p,itype,vcmx25,vcmax_z,psnsun,psnsha,fpsn")
    for r in rows
        @printf(io, "%s,%s,%d,%d,%.8g,%.8g,%.8g,%.8g,%.8g\n",
            r.tag, r.date, r.p, r.itype, r.vcmx25, r.vcmax_z, r.psnsun, r.psnsha, r.fpsn)
    end
end

off = Dict((r.date, r.p) => r for r in rows if r.tag == "lunaoff")
on  = Dict((r.date, r.p) => r for r in rows if r.tag == "lunaon")

function summarize(off, on)
    println("date        p itype   vcmx25_on  vcmax_z_off  vcmax_z_on    fpsn_off     fpsn_on   dfpsn%")
    ndiff = 0; maxrel = 0.0
    for k in sort(collect(keys(on)))
        haskey(off, k) || continue
        a = off[k]; b = on[k]
        rel = a.fpsn == 0 ? 0.0 : 100 * (b.fpsn - a.fpsn) / abs(a.fpsn)
        b.fpsn != a.fpsn && (ndiff += 1)
        maxrel = max(maxrel, abs(rel))
        @printf("%s %2d %5d  %10.4f %12.5f %12.5f %11.6f %11.6f %8.3f\n",
            Dates.format(k[1], "yyyy-mm-dd"), k[2], b.itype,
            b.vcmx25, a.vcmax_z, b.vcmax_z, a.fpsn, b.fpsn, rel)
    end
    # LUNA-off never allocates vcmx25_z, so "acclimated" means the LUNA run's own
    # vcmx25_z has moved off its 30.0 cold-start value.
    nacc = count(k -> isfinite(on[k].vcmx25) && on[k].vcmx25 != 30.0, keys(on))
    @printf("\nsamples where LUNA ACCLIMATED (vcmx25 != 30.0 cold start): %d / %d\n",
            nacc, length(on))
    @printf("samples where that reached GPP  (fpsn differs): %d   max |dfpsn| = %.4f%%\n",
            ndiff, maxrel)
    return (ndiff, nacc)
end
ndiff, nacc = summarize(off, on)
println(ndiff == 0 && nacc > 0 ?
    "\nVACUOUS: LUNA acclimates but the branch is missing — GPP is bit-identical." :
    "\nLIVE: the LUNA branch reaches GPP.")
@printf("csv: %s\n", OUT)
