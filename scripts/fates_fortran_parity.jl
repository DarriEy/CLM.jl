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
#   julia +1.12 --project=scripts scripts/fates_fortran_parity.jl compare <fortran_ref.txt>
#   #   (or set FATES_FORTRAN_REF=<path>)
#
#   FATES_PARITY_STEPS=N   also step N times and dump nstep=N (needs matched forcing)
#   FATES_PARITY_OUT=path  override the Julia dump output path
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

"""
    compare_dumps(jl_path, ft_path; tol=1e-10) -> nfail

Field-by-field parity scorecard between a Julia dump and a Fortran reference
dump of the same schema. Reports, per field name, the max abs and max rel diff
over all matched records, plus any records present in one dump but not the other.
"""
function compare_dumps(jl_path::AbstractString, ft_path::AbstractString; tol::Float64=1e-10)
    J = parse_dump(jl_path); F = parse_dump(ft_path)
    println("="^78)
    println("  FATES Fortran-parity scorecard")
    println("    julia   : $jl_path  ($(length(J)) records)")
    println("    fortran : $ft_path  ($(length(F)) records)")
    println("="^78)

    only_j = setdiff(keys(J), keys(F)); only_f = setdiff(keys(F), keys(J))
    if !isempty(only_j) || !isempty(only_f)
        println("  ⚠ record-set mismatch:")
        isempty(only_j) || println("    only in julia   ($(length(only_j))): ",
            join(sort(collect(only_j))[1:min(end, 8)], ", "))
        isempty(only_f) || println("    only in fortran ($(length(only_f))): ",
            join(sort(collect(only_f))[1:min(end, 8)], ", "))
    end

    # accumulate per-field stats over matched records
    stats = Dict{String,NamedTuple{(:maxabs, :maxrel, :n, :worst),
                                    Tuple{Float64,Float64,Int,String}}}()
    for key in intersect(keys(J), keys(F))
        jf = J[key]; ff = F[key]
        for (fld, jv) in jf
            haskey(ff, fld) || continue
            fv = ff[fld]
            a = abs(jv - fv); r = reld(jv, fv)
            prev = get(stats, fld, (maxabs=0.0, maxrel=0.0, n=0, worst=""))
            newworst = r > prev.maxrel ? key : prev.worst
            stats[fld] = (maxabs=max(prev.maxabs, a), maxrel=max(prev.maxrel, r),
                          n=prev.n + 1, worst=newworst)
        end
    end

    @printf("\n  %-14s %6s %14s %14s %8s\n", "field", "n", "max|abs|", "max rel", "verdict")
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
    println("\n  " * (nfail == 0 && isempty(only_j) && isempty(only_f) ?
        "★ FATES port matches Fortran ground truth within tol=$tol (all fields)" :
        "✗ $nfail field(s) diverge" *
        (isempty(only_j) && isempty(only_f) ? "" : " + record-set mismatch")))
    return nfail + length(only_j) + length(only_f)
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
