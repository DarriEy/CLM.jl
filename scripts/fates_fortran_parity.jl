#!/usr/bin/env julia
# ==========================================================================
# fates_fortran_parity.jl — Fortran-FATES ground-truth parity harness (Julia side).
#
# This is the FATES analogue of the CLM `fortran_pdump` parity machinery
# (scripts/validation/fortran_pdump/ + scripts/fortran_parity_common.jl): it
# extracts a *canonical, precision-preserving* snapshot of the Julia FATES port's
# per-cohort / per-patch / per-site state and writes it to a schema-versioned
# text "pdump", so it can be diffed field-by-field against an equivalent dump
# emitted by an instrumented Fortran CTSM-FATES run.
#
# WHY A TEXT DUMP (not NetCDF like the CLM pdump): FATES state is a pointer-linked
# site→patch→cohort hierarchy, NOT a flat column/pft array, so the 16-field NetCDF
# `pdumpMod` layout does not apply. A line-oriented keyed text schema (below) is
# (a) trivial for a Fortran SourceMod to emit with a plain `write(unit,'(...)')`,
# (b) dependency-free on the Julia side (scripts/ project has no NCDatasets), and
# (c) full float precision (%.17g) so it is bit-parity-capable.
#
# ── PRIMARY PARITY TARGET: COLD START (nstep=0) ────────────────────────────
# The cleanest first comparison is the FATES natural-bare-ground cold start
# (EDInitMod.init_cohorts!): it seeds one cohort per FATES default PFT and is
# FULLY DETERMINISTIC given the FATES parameter file + surfdata soil layering.
# Crucially, the Julia port reads its FATES parameters from
# `data/fates/fates_params_default.cdl`, which is BYTE-IDENTICAL to the CTSM
# install's `src/fates/parameter_files/fates_params_default.cdl` (verified). So
# cold-start cohort dbh / n / height / carbon-pools / pft need NO forcing
# alignment to compare — they depend only on the (shared) parameters. That makes
# nstep=0 the ideal first field-by-field parity scorecard.
#
# Optional: FATES_PARITY_STEPS>0 also steps the model (Aripuana tropical forcing)
# and dumps nstep=N — but a stepped comparison additionally requires the Fortran
# run to use the identical forcing + soil boundary conditions, so treat it as a
# second-tier target.
#
# ── DUMP SCHEMA v1 (one record per line, `TYPE key=val ...`) ───────────────
#   META  nstep=<i> nsites=<i> ncoh=<i> npatch=<i>
#   SITE  s=<i> carbon=<f> maxdbh=<f> gpp=<f> npp=<f> btran=<f>
#   PATCH s=<i> p=<i> patchno=<i> area=<f> elai=<f> btran_ft=<f> ncoh=<i>
#   COHORT s=<i> p=<i> c=<i> pft=<i> n=<f> dbh=<f> height=<f> \
#          leafc=<f> fnrtc=<f> sapwc=<f> storec=<f> structc=<f> \
#          treelai=<f> carea=<f> canopy_layer=<i> status=<f> \
#          gpp=<f> npp=<f> vcmax=<f>
#   patches numbered p=1.. in oldest→youngest order; cohorts c=1.. in
#   tallest→shorter order within a patch. A Fortran instrument MUST traverse in
#   the same order (site%oldest_patch→younger, patch%tallest→shorter).
#
# ── USAGE ──────────────────────────────────────────────────────────────────
#   # write the Julia cold-start dump (default → scripts/fates_pdump_julia_n0.txt)
#   julia +1.12 --project=scripts scripts/fates_fortran_parity.jl
#
#   # compare the live Julia cold start against a Fortran reference dump
#   julia +1.12 --project=scripts scripts/fates_fortran_parity.jl compare \
#       scripts/validation/fates_fortran_parity/fates_pdump_fortran_n0.txt
#   #   (or set FATES_FORTRAN_REF=<path>)
#
#   FATES_PARITY_STEPS=N   also step N times and dump nstep=N (needs matched forcing)
#   FATES_PARITY_OUT=path  override the Julia dump output path
#
# ── GROUND TRUTH: IT EXISTS ────────────────────────────────────────────────
# scripts/validation/fates_fortran_parity/ now holds a REAL Fortran CTSM-FATES
# cold-start dump (1 site x 1 patch x 14 cohorts), produced by an instrumented
# single-point ctsm5.3.012 run — see that directory's README + setup_case.sh.
# Current verdict: the Julia port matches it on ALL 27 fields — 13 of them at
# exactly 0.0 relative difference (dbh/height/n/pft/carea/treelai/vcmax/
# canopy_layer/status/maxdbh/...), the five PRT carbon pools to 1 ULP (~7e-18),
# and site total carbon to 7e-16 (summation order). No physics disagreement.
# ==========================================================================
using CLM, Printf, Dates
const _C = CLM
include(joinpath(@__DIR__, "fates_longhorizon.jl"))   # reuse build()/census/helpers (guarded exit)

# ---- canonical field extraction -------------------------------------------
# Cohort carbon pools [kgC/indiv] via the PRT state accessor, mirroring
# total_site_carbon() in fates_longhorizon.jl.
_pool(cc, organ) = (cc.prt === nothing ? NaN :
                    _C.GetState(cc.prt, organ, _C.carbon12_element))

"""
    extract_records(inst, nstep) -> Vector{String}

Walk the FATES site→patch→cohort hierarchy and produce the schema-v1 dump lines.
"""
function extract_records(inst, nstep::Int)
    recs = String[]
    fates = inst.fates
    nsites = fates.nsites
    # canopystate.elai_patch is an HLM (patch) array; for the single FATES column
    # the vegetated HLM patches follow the bare-ground patch. We report elai per
    # FATES patch positionally (p-th FATES patch → HLM veg patch p+1) as a coupling
    # cross-check; it is a derived field, not a FATES-internal state.
    elai = try; Array(inst.canopystate.elai_patch); catch; Float64[]; end

    ncoh_tot = 0; npatch_tot = 0
    site_lines = String[]; patch_lines = String[]; coh_lines = String[]
    for s in 1:nsites
        site = fates.sites[s]
        carbon = total_site_carbon(site)
        md = max_dbh(site)
        gpp = 0.0; npp = 0.0; btran_site = 0.0
        p = 0
        cp = site.oldest_patch
        while cp !== nothing
            p += 1; npatch_tot += 1
            btr = 0.0
            try; btr = maximum(x -> isfinite(x) ? x : 0.0, cp.btran_ft); catch; end
            btran_site = max(btran_site, btr)
            # elai for this positional veg patch (HLM patch index = p+1 when a
            # bare-ground patch precedes the veg patches; guard bounds).
            ev = (length(elai) >= p + 1) ? Float64(elai[p + 1]) : NaN
            ncoh_p = 0
            c = 0
            cc = cp.tallest
            while cc !== nothing
                c += 1; ncoh_p += 1; ncoh_tot += 1
                g = isfinite(cc.gpp_tstep) ? cc.gpp_tstep : 0.0
                nn = isfinite(cc.npp_tstep) ? cc.npp_tstep : 0.0
                gpp += g * cc.n; npp += nn * cc.n
                push!(coh_lines, @sprintf(
                    "COHORT s=%d p=%d c=%d pft=%d n=%.17g dbh=%.17g height=%.17g leafc=%.17g fnrtc=%.17g sapwc=%.17g storec=%.17g structc=%.17g treelai=%.17g carea=%.17g canopy_layer=%d status=%.17g gpp=%.17g npp=%.17g vcmax=%.17g",
                    s, p, c, cc.pft, cc.n, cc.dbh, cc.height,
                    _pool(cc, _C.leaf_organ), _pool(cc, _C.fnrt_organ),
                    _pool(cc, _C.sapw_organ), _pool(cc, _C.store_organ),
                    _pool(cc, _C.struct_organ),
                    cc.treelai, cc.c_area, cc.canopy_layer,
                    Float64(cc.status_coh),
                    (isfinite(cc.gpp_tstep) ? cc.gpp_tstep : NaN),
                    (isfinite(cc.npp_tstep) ? cc.npp_tstep : NaN),
                    (isfinite(cc.vcmax25top) ? cc.vcmax25top : NaN)))
                cc = cc.shorter
            end
            push!(patch_lines, @sprintf(
                "PATCH s=%d p=%d patchno=%d area=%.17g elai=%.17g btran_ft=%.17g ncoh=%d",
                s, p, cp.patchno, cp.area, ev, btr, ncoh_p))
            cp = cp.younger
        end
        push!(site_lines, @sprintf(
            "SITE s=%d carbon=%.17g maxdbh=%.17g gpp=%.17g npp=%.17g btran=%.17g",
            s, carbon, md, gpp, npp, btran_site))
    end
    push!(recs, @sprintf("META nstep=%d nsites=%d ncoh=%d npatch=%d",
                         nstep, nsites, ncoh_tot, npatch_tot))
    append!(recs, site_lines); append!(recs, patch_lines); append!(recs, coh_lines)
    return recs
end

function write_dump(path::AbstractString, recs::Vector{String})
    open(path, "w") do io
        println(io, "# FATES parity dump — schema v1 — generator=julia")
        for r in recs; println(io, r); end
    end
    return path
end

# ---- dump parsing + comparison --------------------------------------------
"""
    parse_dump(path) -> Dict{String,Dict{String,Float64}}

Parse a schema-v1 dump into a map from a record key (e.g. "COHORT|s=1|p=1|c=3")
to its numeric field map. The record key is built from the *identity* fields
(TYPE + the s/p/c ordinals) so Julia and Fortran records line up.
"""
function parse_dump(path::AbstractString)
    out = Dict{String,Dict{String,Float64}}()
    for line in eachline(path)
        (isempty(line) || startswith(line, "#")) && continue
        toks = split(strip(line))
        isempty(toks) && continue
        typ = toks[1]
        fields = Dict{String,Float64}()
        idbits = String[typ]
        for t in toks[2:end]
            kv = split(t, "=", limit=2)
            length(kv) == 2 || continue
            k, v = kv[1], kv[2]
            val = tryparse(Float64, v)
            val === nothing && continue
            fields[k] = val
            # identity keys: ordinals that address the record
            if (typ == "SITE"   && k == "s") ||
               (typ == "PATCH"  && k in ("s", "p")) ||
               (typ == "COHORT" && k in ("s", "p", "c"))
                push!(idbits, "$k=$(v)")
            end
        end
        key = join(idbits, "|")
        out[key] = fields
    end
    return out
end

reld(a, b) = abs(a - b) / (1.0 + abs(a))

# Fields the Fortran instrument deliberately does NOT report, and the sentinel it
# writes instead (see FatesParityDumpMod.F90). `elai` is an HLM canopystate
# coupling diagnostic, not a FATES-internal state, so the Fortran side emits -1
# and the comparison skips it rather than scoring a meaningless difference.
const FORTRAN_SENTINEL = Dict("elai" => -1.0)

# META records carry *structural totals* (nsites/npatch/ncoh) that depend on how
# many grid columns the harness instantiated, NOT on FATES physics. They are
# reconciled per-site in the structure section instead of being scored as fields.
const META_STRUCTURAL = Set(["nsites", "npatch", "ncoh"])

"site ordinal of a record key, e.g. \"COHORT|s=2|p=1|c=3\" -> 2 (0 for META)."
function site_of(key::AbstractString)
    m = match(r"\|s=([0-9]+)", key)
    m === nothing ? 0 : parse(Int, m.captures[1])
end

"""
true if two records are numerically identical (replica check).

The site ordinal `s` is excluded: it is the record's own address, so it differs
between a site and its replica by construction. Every other field — including the
`p`/`c` ordinals and all physical state — must be bit-identical.
"""
function same_record(a::Dict{String,Float64}, b::Dict{String,Float64})
    ka = setdiff(keys(a), ["s"]); kb = setdiff(keys(b), ["s"])
    ka == kb || return false
    all(k -> (a[k] == b[k]) || (isnan(a[k]) && isnan(b[k])), ka)
end

"re-key a record onto site `s` (so site-N records can be compared against site-1)."
rekey_site(key::AbstractString, s::Int) = replace(key, r"\|s=[0-9]+" => "|s=$s")

"""
    compare_dumps(jl_path, ft_path; tol=1e-10) -> nfail

Field-by-field parity scorecard between a Julia dump and a Fortran reference dump
of the same schema.

Structure is reconciled BEFORE scoring, because the two sides legitimately
instantiate different numbers of grid columns: the Julia harness builds an
N-gridcell test instance while a single-point Fortran case has one site. Sites
are FATES's independent replicate units, so the honest comparison is:

  * score the sites both dumps have (1..min(nsites_j, nsites_f)) field-by-field;
  * verify any EXTRA Julia sites are exact replicas of the compared ones, and
    report that check — an extra site that is *not* a replica is a real failure;
  * reconcile the META totals per-site (patches/site, cohorts/site) rather than
    scoring raw totals that just count grid columns.

Anything that is not explained by those two structural facts is scored as a real
divergence.
"""
function compare_dumps(jl_path::AbstractString, ft_path::AbstractString; tol::Float64=1e-10)
    J = parse_dump(jl_path); F = parse_dump(ft_path)
    println("="^78)
    println("  FATES Fortran-parity scorecard")
    println("    julia   : $jl_path  ($(length(J)) records)")
    println("    fortran : $ft_path  ($(length(F)) records)")
    println("="^78)

    jmeta = get(J, "META", Dict{String,Float64}())
    fmeta = get(F, "META", Dict{String,Float64}())
    nsj = Int(get(jmeta, "nsites", 0)); nsf = Int(get(fmeta, "nsites", 0))
    ncmp = min(nsj, nsf)
    nstruct = 0   # real structural failures

    # ---- structure ---------------------------------------------------------
    println("\n  ── structure ──")
    @printf("    sites         julia=%-4d fortran=%-4d → scoring the %d common site(s)\n",
            nsj, nsf, ncmp)
    if ncmp < 1
        println("    ✗ no common site to compare"); return 1
    end

    # per-site structural totals (grid-size-independent)
    for (nm, tot) in (("patches/site", "npatch"), ("cohorts/site", "ncoh"))
        pj = get(jmeta, tot, NaN) / max(nsj, 1)
        pf = get(fmeta, tot, NaN) / max(nsf, 1)
        ok = pj == pf
        ok || (nstruct += 1)
        @printf("    %-13s julia=%-8.4g fortran=%-8.4g %s\n", nm, pj, pf,
                ok ? "match" : "DIVERGE")
    end

    # extra Julia sites must be exact replicas of the sites we scored
    if nsj > ncmp
        extra = [k for k in keys(J) if site_of(k) > ncmp]
        nbad = 0
        for k in extra
            ref = rekey_site(k, 1 + (site_of(k) - 1) % ncmp)
            (haskey(J, ref) && same_record(J[k], J[ref])) || (nbad += 1)
        end
        nbad == 0 || (nstruct += 1)
        @printf("    extra sites   julia has %d site(s) beyond the Fortran reference; %s\n",
                nsj - ncmp,
                nbad == 0 ?
                  "all $(length(extra)) records are EXACT replicas of the scored site(s) → not a divergence" :
                  "✗ $nbad of $(length(extra)) records are NOT replicas → REAL divergence")
    end

    # ---- record-set reconciliation over the common sites -------------------
    Jc = Dict(k => v for (k, v) in J if site_of(k) <= ncmp)
    Fc = Dict(k => v for (k, v) in F if site_of(k) <= ncmp)
    only_j = setdiff(keys(Jc), keys(Fc)); only_f = setdiff(keys(Fc), keys(Jc))
    if !isempty(only_j) || !isempty(only_f)
        println("    ✗ record-set mismatch WITHIN the common sites:")
        isempty(only_j) || println("      only in julia   ($(length(only_j))): ",
            join(sort(collect(only_j))[1:min(end, 8)], ", "))
        isempty(only_f) || println("      only in fortran ($(length(only_f))): ",
            join(sort(collect(only_f))[1:min(end, 8)], ", "))
    else
        @printf("    records       %d matched, 0 unmatched, across the common site(s)\n",
                length(intersect(keys(Jc), keys(Fc))))
    end

    # ---- field-by-field scoring -------------------------------------------
    stats = Dict{String,NamedTuple{(:maxabs, :maxrel, :n, :worst),
                                    Tuple{Float64,Float64,Int,String}}}()
    skipped = Dict{String,Int}()
    for key in intersect(keys(Jc), keys(Fc))
        jf = Jc[key]; ff = Fc[key]
        for (fld, jv) in jf
            haskey(ff, fld) || continue
            fv = ff[fld]
            # skip fields the Fortran instrument does not report
            if haskey(FORTRAN_SENTINEL, fld) && fv == FORTRAN_SENTINEL[fld]
                skipped[fld] = get(skipped, fld, 0) + 1
                continue
            end
            # META structural totals are reconciled above, not scored here
            if startswith(key, "META") && fld in META_STRUCTURAL
                skipped[fld] = get(skipped, fld, 0) + 1
                continue
            end
            a = abs(jv - fv); r = reld(jv, fv)
            prev = get(stats, fld, (maxabs=0.0, maxrel=0.0, n=0, worst=""))
            newworst = r > prev.maxrel ? key : prev.worst
            stats[fld] = (maxabs=max(prev.maxabs, a), maxrel=max(prev.maxrel, r),
                          n=prev.n + 1, worst=newworst)
        end
    end

    println("\n  ── fields ──")
    @printf("  %-14s %6s %14s %14s %8s\n", "field", "n", "max|abs|", "max rel", "verdict")
    println("  " * "-"^62)
    nfail = 0
    for fld in sort(collect(keys(stats)))
        st = stats[fld]
        ok = st.maxrel <= tol
        ok || (nfail += 1)
        @printf("  %-14s %6d %14.4e %14.4e %8s\n", fld, st.n, st.maxabs, st.maxrel,
                ok ? "match" : "DIVERGE")
        ok || @printf("       worst @ %s\n", st.worst)
    end
    if !isempty(skipped)
        println("\n  not scored (not reported by the Fortran instrument / structural):")
        for fld in sort(collect(keys(skipped)))
            @printf("    %-14s %d record(s)\n", fld, skipped[fld])
        end
    end

    nbadrec = length(only_j) + length(only_f)
    total = nfail + nstruct + nbadrec
    println("\n  " * (total == 0 ?
        "★ FATES port matches the Fortran ground truth within tol=$tol " *
        "($(length(stats)) fields, all records)" :
        "✗ $nfail field(s) diverge, $nstruct structural, $nbadrec unmatched record(s)"))
    return total
end

# ---- driver ---------------------------------------------------------------
function build_and_maybe_step()
    inst, fates, config, bounds, filt, filt_ia = build()
    nstep = 0
    steps = parse(Int, get(ENV, "FATES_PARITY_STEPS", "0"))
    if steps > 0
        # Step the model with the Aripuana tropical forcing (mirrors
        # fates_longhorizon.main's loop) so nstep=N can be dumped. Only meaningful
        # for parity if the Fortran run uses the identical forcing + soil BC.
        site = fates.sites[1]; photosyns = inst.photosyns
        dtime = 1800.0
        forcing = get(ENV, "FATES_FORCING",
            "$DATA/domain_Aripuana_Amazon/data/forcing/CLM_input/clmforc.2004.nc")
        fyr = parse(Int, get(ENV, "FATES_YEAR", "2004"))
        fr = _C.ForcingReader(); _C.forcing_reader_init!(fr, forcing); fr.interp_time = true
        start_date = DateTime(fyr, 1, 1)
        for i in 1:steps
            step_start = start_date + Second((i - 1) * Int(dtime))
            is_beg = (Dates.hour(step_start) == 0 && Dates.minute(step_start) == 0 &&
                      Dates.second(step_start) == 0)
            _C.read_forcing_step!(fr, inst.atm2lnd, step_start, 1, 1;
                gridcell_latdeg=inst.gridcell.latdeg, gridcell_londeg=inst.gridcell.londeg,
                dtime=Int(dtime))
            inst.atm2lnd.forc_topo_grc[1] = 200.0; inst.topo.topo_col[1] = 200.0
            _C.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
            sod = Dates.hour(step_start) * 3600 + Dates.minute(step_start) * 60 +
                  Dates.second(step_start)
            calday = Dates.dayofyear(step_start) + sod / 86400.0
            (declin, _e) = _C.compute_orbital(calday); nextsw_cday = calday + Int(dtime) / 86400.0
            _C.clm_drv!(config, inst, filt, filt_ia, bounds, true, nextsw_cday, declin, declin,
                0.4091, false, false, "20260101", false; nstep=i, is_first_step=(i == 1),
                is_beg_curr_day=is_beg, is_end_curr_day=false, is_beg_curr_year=false,
                dtime=dtime, mon=Dates.month(step_start), day=Dates.day(step_start),
                secs=sod, jday=Dates.dayofyear(step_start), photosyns=photosyns)
        end
        _C.forcing_reader_close!(fr)
        nstep = steps
    end
    return inst, nstep
end

function main()
    args = ARGS
    mode = (isempty(args) ? "dump" : args[1])
    println("="^78)
    println("  fates_fortran_parity.jl — mode=$mode")
    println("="^78)

    inst, nstep = build_and_maybe_step()
    recs = extract_records(inst, nstep)
    outpath = get(ENV, "FATES_PARITY_OUT",
                  joinpath(@__DIR__, "fates_pdump_julia_n$(nstep).txt"))
    write_dump(outpath, recs)

    # human-readable cold-start census summary
    cen = census(inst)
    c = total_site_carbon(inst.fates.sites[1]); md = max_dbh(inst.fates.sites[1])
    @printf("\n  Julia state @ nstep=%d : ncoh=%d npatch=%d carbon=%.6g maxdbh=%.6f\n",
            nstep, cen.ncoh, cen.npatch, c, md)
    println("  wrote Julia dump → $outpath  ($(length(recs)) records)")

    if mode == "compare"
        ref = length(args) >= 2 ? args[2] : get(ENV, "FATES_FORTRAN_REF", "")
        if isempty(ref) || !isfile(ref)
            println("\n  ✗ compare mode needs a Fortran reference dump: pass it as arg 2 or set FATES_FORTRAN_REF.")
            println("    No Fortran CTSM-FATES reference dump exists yet in this environment (see harness header).")
            return 2
        end
        tol = parse(Float64, get(ENV, "FATES_PARITY_TOL", "1e-10"))
        return compare_dumps(outpath, ref; tol=tol)
    else
        println("\n  (dump mode) To score against Fortran ground truth once a reference exists:")
        println("    julia +1.12 --project=scripts scripts/fates_fortran_parity.jl compare <fortran_ref.txt>")
        return 0
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
