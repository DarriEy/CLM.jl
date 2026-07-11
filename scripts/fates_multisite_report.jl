#!/usr/bin/env julia
# ==========================================================================
# fates_multisite_report.jl — aggregate the per-site summary lines written by
# fates_multisite_validation.jl into a cross-biome equilibrium scorecard.
#
#   FATES_MS_RESULTS=/path/results.csv julia +1.12 --project=. \
#       scripts/fates_multisite_report.jl
# ==========================================================================
using Printf
resfile = get(ENV, "FATES_MS_RESULTS", length(ARGS) >= 1 ? ARGS[1] : "")
isfile(resfile) || (println("no results file: $resfile"); exit(1))
rows = Any[]
for (li, ln) in enumerate(eachline(resfile))
    li == 1 && startswith(ln, "key\t") && continue   # header
    f = split(ln, '\t'); length(f) < 20 && continue
    push!(rows, (key=f[1], label=f[2], days=parse(Int,f[3]),
        c0=parse(Float64,f[4]), cF=parse(Float64,f[5]), minc=parse(Float64,f[6]),
        maxc=parse(Float64,f[7]), peak=parse(Float64,f[8]), bandratio=parse(Float64,f[9]),
        minncoh=parse(Int,f[10]), maxncoh=parse(Int,f[11]),
        elaimin=parse(Float64,f[12]), elaimax=parse(Float64,f[13]), maxdbh=parse(Float64,f[14]),
        balok=parse(Int,f[15]), daysadv=parse(Int,f[16]), collapsed=parse(Int,f[17]),
        boombust=parse(Int,f[18]), inband=parse(Int,f[19]), nfail=parse(Int,f[20])))
end
isempty(rows) && (println("no data rows in $resfile"); exit(1))

println("="^108)
println("  FATES MULTI-SITE EQUILIBRIUM SCORECARD  (screened vs baseline; climate-appropriate PFT sets)")
println("="^108)
@printf("  %-24s %5s %8s %8s %8s %9s %7s %8s %7s %6s %6s\n",
        "site", "days", "cold-C", "final-C", "peak-C", "hi/lo½", "ncoh", "maxdbh", "elai", "bal", "verdict")
println("  " * "-"^104)
for r in rows
    verdict = r.nfail == 0 ? "HEALTHY" : "$(r.nfail)FAIL"
    balstr = @sprintf("%d/%d", r.balok, r.daysadv)
    @printf("  %-24s %5d %8.0f %8.0f %8.0f %9.2f %3d-%-3d %8.2f %5.2f-%-4.2f %7s %s\n",
        r.key, r.days, r.c0, r.cF, r.peak, r.bandratio, r.minncoh, r.maxncoh,
        r.maxdbh, r.elaimin, r.elaimax, balstr, verdict)
end
println("  " * "-"^104)
println("  cold-C/final-C/peak-C = site carbon [gC/m2]; hi/lo½ = 2nd-half carbon band ratio (→1 = equilibrated)")
println("  flags: collapsed(final<20%peak) boom-bust(rise&fall>2.5x) in-band(final in biome window)\n")

# ---- screened vs baseline contrast ----
function findrow(k); i = findfirst(r->r.key==k, rows); i===nothing ? nothing : rows[i]; end
println("  SCREENED vs BASELINE (the PFT screen's effect):")
# A run whose days_advanced falls well short of its intended horizon terminated on a
# numerical ERROR (the all-PFT cohort explosion -> imaginary-roots crash), not a clean
# steady state — flag it distinctly from a demographic collapse/boom-bust.
maxdays = maximum(r->r.days, rows)
tag(r) = string(r.days < maxdays - 30 ? "[CRASH@day$(r.days)]" : "",
                r.collapsed==1 ? "[COLLAPSE]" : "", r.boombust==1 ? "[BOOMBUST]" : "",
                (r.days >= maxdays-30 && r.collapsed==0 && r.boombust==0) ? "[stable]" : "")
for base in ("aripuana", "krycklan", "hubbardbrook")
    s = findrow("$(base)_screened"); b = findrow("$(base)_baseline")
    (s === nothing && b === nothing) && continue
    s !== nothing && @printf("   %-14s  screened: final=%.0f peak=%.0f band=%.2f verdict=%s %s\n",
        base, s.cF, s.peak, s.bandratio, s.nfail==0 ? "HEALTHY" : "$(s.nfail)FAIL", tag(s))
    b !== nothing && @printf("   %-14s  baseline: final=%.0f peak=%.0f band=%.2f verdict=%s %s\n",
        "", b.cF, b.peak, b.bandratio, b.nfail==0 ? "HEALTHY" : "$(b.nfail)FAIL", tag(b))
end
nhealthy = count(r->r.nfail==0, rows)
@printf("\n  %d/%d runs reached a healthy verdict.\n", nhealthy, length(rows))
