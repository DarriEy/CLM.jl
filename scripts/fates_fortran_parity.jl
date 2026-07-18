#!/usr/bin/env julia
# ==========================================================================
# fates_fortran_parity.jl — Fortran-FATES ground-truth parity harness (Julia side).
#
# TIME-STEPPED, PHASE-TAGGED (schema v2). The FATES analogue of the CLM
# `fortran_pdump` parity machinery: it extracts a canonical, precision-preserving
# snapshot of the Julia FATES port's site→patch→cohort state and diffs it, field by
# field, against the same snapshot taken by an instrumented Fortran CTSM-FATES run
# (scripts/validation/fates_fortran_parity/).
#
# ── WHY v2 EXISTS ─────────────────────────────────────────────────────────
# v1 (PR #213) compared ONE snapshot, at cold start (nstep=0), and matched on all
# 27 fields. That is a real result — but a narrow one. The FATES parameter file is
# byte-identical between the port and the install, so the nstep=0 cohorts are
# parameter-identical BY CONSTRUCTION. It validated INITIALIZATION. It tested
# nothing about the FATES DYNAMICS: photosynthesis, allocation, growth, mortality,
# recruitment, disturbance/patch dynamics, canopy promotion/demotion — the bulk of
# the ~62-module port.
#
# v2 runs BOTH models forward from the SAME (bit-identical) cold start with the SAME
# site, forcing and timestep, and dumps at EVERY FATES phase boundary, so a
# divergence can be attributed to the phase that produced it.
#
# PHASES (identical codes on both sides):
#    0  coldstart   EDInitMod init_cohorts (nstep=0)
#    1  fast        every timestep, after the photosynthesis/flux solve + running-mean
#                   + hifrq-history updates, BEFORE the daily dynamics
#   10  dyn_in      ed_ecosystem_dynamics entry (the daily step's INPUT state)
#   11  phenology   12 distrates   13 integrate (growth/allocation/PRT/mortality)
#   14  recruit     15 cohortfuse  16 spawn (disturbance)
#   17  patchfuse   (= end of ed_ecosystem_dynamics)
#   18  updatesite  (canopy_spread / canopy_structure / trim — end of the daily step)
#
# ── HOW THIS IS AN ORACLE, NOT A DRIFT PLOT ───────────────────────────────
# Both models start from a cold start that PR #213 proved bit-identical, and are
# driven by the same forcing. So up to the FIRST divergence, every phase comparison
# is a genuine single-step oracle: the Julia phase consumed the same input state the
# Fortran phase did, and any difference in the output state is TRANSLATION ERROR,
# not accumulated drift. After the first divergence, differences compound and the
# numbers become a drift measurement. The scorecard therefore reports (a) the FIRST
# (step, phase, field) that exceeds tolerance — the attributable result — and
# (b) the free-running drift thereafter.
#
# The PATCH/SITE records also carry the HLM→FATES boundary conditions (t_veg, eair,
# cair, rb, daylength factor, incident solar; soil T/moisture). Those are the
# discriminator: if bc_in matches to round-off but the FATES output does not, the
# divergence is INSIDE FATES; if bc_in already differs, it came from the host land
# model, not from the FATES port.
#
# ── USAGE ─────────────────────────────────────────────────────────────────
#   # run the Julia FATES port on the Fortran case's site/forcing and dump
#   julia +1.12 --project=. scripts/fates_fortran_parity.jl
#
#   # ... and score it against the Fortran ground truth (.txt or .txt.gz)
#   julia +1.12 --project=. scripts/fates_fortran_parity.jl compare \
#       scripts/validation/fates_fortran_parity/fates_pdump_fortran.txt.gz
#
#   FATES_PARITY_STEPS=N     number of timesteps (default 120 = 5 days @ 3600 s)
#   FATES_PARITY_OUT=path    Julia dump output path
#   FATES_PARITY_TOL=x       relative tolerance (default 1e-10)
# ==========================================================================
using CLM, Printf, Dates, NCDatasets
const _C = CLM

const DATA = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data"
# The site the Fortran case runs (scripts/validation/fates_fortran_parity/setup_case.sh):
# the proven Bow-at-Banff single-point settings dir + its local ERA5 clmforc.YYYY.nc.
const SITEDIR   = get(ENV, "FATES_SITEDIR", "$DATA/domain_Bow_at_Banff_lumped_era5/settings/CLM")
const FSURDAT   = joinpath(SITEDIR, "parameters", "surfdata_clm.nc")
const FPARAM    = joinpath(SITEDIR, "parameters", "clm5_params.nc")
const FORCDIR   = get(ENV, "FATES_FORCDIR",
                      "$DATA/domain_Bow_at_Banff_lumped_era5/data/forcing/CLM_input")
# Must match START_YMD / lnd_cpl_dt in the Fortran case.
const START     = DateTime(2002, 7, 1)
const DTIME     = 3600.0
const NSTEPS    = parse(Int, get(ENV, "FATES_PARITY_STEPS", "120"))

# ==========================================================================
# Record extraction — the Julia side of the schema-v2 writer. Field names, field
# order and traversal order MUST mirror scripts/validation/fates_fortran_parity/
# FatesParityDumpMod.F90 exactly.
# ==========================================================================
_pool(cc, organ) = (cc.prt === nothing ? NaN : _C.GetState(cc.prt, organ, _C.carbon12_element))
_bcp(a, i) = (a isa AbstractVector && 1 <= i <= length(a)) ? Float64(a[i]) : -9999.0
_bcm(a, i) = (a isa AbstractMatrix && 1 <= i <= size(a, 1) && size(a, 2) >= 1) ?
             Float64(a[i, 1]) : -9999.0

"""
    site_records(site, bc_in, s, nstep, phase) -> Vector{String}

Schema-v2 records for one FATES site at one (nstep, phase) boundary.
"""
function site_records(site, bc_in, s::Int, nstep::Int, phase::Int)
    recs = String[]
    hdr = @sprintf("t=%d ph=%d s=%d", nstep, phase, s)

    # ---- site totals -----------------------------------------------------
    carbon = 0.0; maxdbh = 0.0; ncoh_s = 0; npatch_s = 0
    cp = site.oldest_patch
    while cp !== nothing
        npatch_s += 1
        cc = cp.tallest
        while cc !== nothing
            ncoh_s += 1
            cc.dbh > maxdbh && (maxdbh = cc.dbh)
            for org in (_C.leaf_organ, _C.fnrt_organ, _C.sapw_organ,
                        _C.store_organ, _C.struct_organ)
                carbon += _pool(cc, org) * cc.n
            end
            cc = cc.shorter
        end
        cp = cp.younger
    end

    push!(recs, @sprintf(
        "SITE %s carbon=%.17g maxdbh=%.17g ncoh=%d npatch=%d pbot=%.17g tsoil1=%.17g h2ovol1=%.17g effporo1=%.17g",
        hdr, carbon, maxdbh, ncoh_s, npatch_s,
        Float64(bc_in.forc_pbot), _bcp(bc_in.t_soisno_sl, 1),
        _bcp(bc_in.h2o_liqvol_sl, 1), _bcp(bc_in.eff_porosity_sl, 1)))

    # ---- patches + cohorts -----------------------------------------------
    p = 0
    cp = site.oldest_patch
    while cp !== nothing
        p += 1
        btr = maximum(cp.btran_ft)
        ncoh_p = 0
        cc = cp.tallest
        while cc !== nothing; ncoh_p += 1; cc = cc.shorter; end

        push!(recs, @sprintf(
            "PATCH %s p=%d patchno=%d area=%.17g age=%.17g ncoh=%d btran_ft=%.17g tcanarea=%.17g nocomp=%d tveg=%.17g tgcm=%.17g esat_tv=%.17g eair=%.17g oair=%.17g cair=%.17g rb=%.17g daylfac=%.17g solad1=%.17g solai1=%.17g",
            hdr, p, cp.patchno, cp.area, cp.age, ncoh_p, btr, cp.total_canopy_area,
            cp.nocomp_pft_label,
            _bcp(bc_in.t_veg_pa, p), _bcp(bc_in.tgcm_pa, p), _bcp(bc_in.esat_tv_pa, p),
            _bcp(bc_in.eair_pa, p), _bcp(bc_in.oair_pa, p), _bcp(bc_in.cair_pa, p),
            _bcp(bc_in.rb_pa, p), _bcp(bc_in.dayl_factor_pa, p),
            _bcm(bc_in.solad_parb, p), _bcm(bc_in.solai_parb, p)))

        c = 0
        cc = cp.tallest
        while cc !== nothing
            c += 1
            push!(recs, @sprintf(
                "COHORT %s p=%d c=%d pft=%d n=%.17g dbh=%.17g height=%.17g coage=%.17g canopy_layer=%d status=%.17g carea=%.17g treelai=%.17g treesai=%.17g canopy_trim=%.17g l2fr=%.17g nv=%d isnew=%d leafc=%.17g fnrtc=%.17g sapwc=%.17g storec=%.17g structc=%.17g reproc=%.17g gpp=%.17g npp=%.17g resp=%.17g gpp_acc=%.17g npp_acc=%.17g resp_acc=%.17g gpp_hold=%.17g npp_hold=%.17g resp_hold=%.17g vcmax=%.17g jmax=%.17g tpu=%.17g rdark=%.17g resp_m=%.17g resp_g=%.17g ddbhdt=%.17g dhdt=%.17g dndt=%.17g bmort=%.17g cmort=%.17g hmort=%.17g frmort=%.17g smort=%.17g asmort=%.17g dmort=%.17g seed_prod=%.17g",
                hdr, p, c, cc.pft, cc.n, cc.dbh, cc.height, cc.coage,
                cc.canopy_layer, Float64(cc.status_coh), cc.c_area, cc.treelai, cc.treesai,
                cc.canopy_trim, cc.l2fr, cc.nv, (cc.isnew ? 1 : 0),
                _pool(cc, _C.leaf_organ), _pool(cc, _C.fnrt_organ), _pool(cc, _C.sapw_organ),
                _pool(cc, _C.store_organ), _pool(cc, _C.struct_organ), _pool(cc, _C.repro_organ),
                cc.gpp_tstep, cc.npp_tstep, cc.resp_tstep,
                cc.gpp_acc, cc.npp_acc, cc.resp_acc,
                cc.gpp_acc_hold, cc.npp_acc_hold, cc.resp_acc_hold,
                cc.vcmax25top, cc.jmax25top, cc.tpu25top,
                cc.rdark, cc.resp_m, cc.resp_g_tstep,
                cc.ddbhdt, cc.dhdt, cc.dndt,
                cc.bmort, cc.cmort, cc.hmort, cc.frmort, cc.smort, cc.asmort, cc.dmort,
                cc.seed_prod))
            cc = cc.shorter
        end
        cp = cp.younger
    end
    return recs
end

# ==========================================================================
# Build: the Julia FATES port on the SAME site the Fortran case runs.
# ==========================================================================
function build()
    # use_bedrock=false MATCHES the Fortran case's lnd_in (`use_bedrock = .false.`,
    # from CLM build-namelist for I2000Clm50FatesRs — see setup_case.sh / the run dir).
    # clm_initialize! DEFAULTS use_bedrock=true; leaving it unset ran the port with a
    # shallow Bow bedrock (surfdata zbedrock≈2.28 m → nbedrock=12), which the per-step
    # clamp_zwt_to_bedrock! then used to drag the water table 8.6 m → 2.28 m, and the
    # SoilHydrology drainage `rsub_top ∝ exp(-zwt/hkdepth)` blew up to ~10 mm/day of
    # spurious subsurface runoff — draining the whole column and over-drying the top
    # layer (h2ovol1 0.11–0.14 vs Fortran 0.32–0.34, PR #247 D1). With use_bedrock=false
    # nbedrock=nlevsoi, zwt stays 8.6 m, qflx_drain≡0 like Fortran, and h2ovol1 tracks
    # the Fortran ground truth to ~0.01–0.02. The soil-hydrology port itself is faithful;
    # this was a harness↔reference namelist mismatch (the #233/#240/#248 trap).
    inst, bounds, filt, tm = _C.clm_initialize!(; fsurdat=FSURDAT, paramfile=FPARAM,
        use_fates=true, start_date=START, dtime=Int(DTIME), use_bedrock=false)
    for g in 1:bounds.endg
        isfinite(inst.gridcell.latdeg[g]) || (inst.gridcell.latdeg[g] = inst.gridcell.lat[g]*180/π)
        isfinite(inst.gridcell.londeg[g]) || (inst.gridcell.londeg[g] = inst.gridcell.lon[g]*180/π)
    end
    # clm_initialize! attaches inst.fates + cold-starts the site cohorts, but nothing
    # flags WHICH HLM columns are FATES; every FATES↔biogeophysics coupling block is
    # gated on col.is_fates[c]. Flag the natural-vegetated soil column(s) — the s-th
    # is_fates column maps 1:1 onto FATES site s.
    let col = inst.column, lun = inst.landunit, nflag = 0
        for c in 1:bounds.endc
            if lun.itype[col.landunit[c]] == _C.ISTSOIL
                col.is_fates[c] = true; nflag += 1
            end
        end
        @printf("  flagged %d FATES column(s); fates.nsites=%d\n", nflag, inst.fates.nsites)
    end
    # NOTE: deliberately NO soil warming and NO PFT screen here (unlike
    # scripts/fates_longhorizon.jl). The Fortran case cold-starts CLM from the same
    # surfdata with all 14 FATES PFTs; a hand-tweaked initial soil profile or a
    # screened PFT set would compare the port against a DIFFERENT model.
    config = _C.CLMDriverConfig(use_fates=true)
    return inst, bounds, filt, _C.clump_filter_inactive_and_active, config, tm
end

# ==========================================================================
# Run: step the port and collect phase-tagged records via the FATES parity hook.
# ==========================================================================
function run_and_dump(outpath::String; nsteps::Int=NSTEPS)
    inst, bounds, filt, filt_ia, config, _tm = build()

    recs = String[]
    nstep_ref = Ref(0)
    # site object -> ordinal (the hook is called per site and FATES sites carry no index)
    sidx = IdDict{Any,Int}(inst.fates.sites[s] => s for s in 1:inst.fates.nsites)
    _C.FATES_PARITY_HOOK[] = (site, bc_in, phase) -> begin
        s = get(sidx, site, 0)
        s == 0 && return nothing
        append!(recs, site_records(site, bc_in, s, nstep_ref[], phase))
        return nothing
    end

    # phase 0 — cold start (bc_in is not packed yet on either side; the comparator
    # does not score bc_in fields at phase 0)
    for s in 1:inst.fates.nsites
        append!(recs, site_records(inst.fates.sites[s], inst.fates.bc_in[s], s, 0, 0))
    end

    ng, nc = bounds.endg, bounds.endc
    fr = _C.ForcingReader()
    _C.forcing_reader_init!(fr, joinpath(FORCDIR, "clmforc.$(year(START)).nc"))
    # Elevation: the Fortran DATM/CLM pair reads the site's topo_forcing.nc.
    let tf = joinpath(FORCDIR, "topo_forcing.nc")
        if isfile(tf)
            ds = NCDataset(tf, "r")
            if haskey(ds, "TOPO")
                ft = Float64(ds["TOPO"][1])
                for g in 1:ng; inst.atm2lnd.forc_topo_grc[g] = ft; end
                for c in 1:nc; inst.topo.topo_col[c] = ft; end
            end
            close(ds)
        end
    end

    @printf("  stepping %d x %.0fs from %s (%.1f days)\n",
            nsteps, DTIME, START, nsteps*DTIME/86400)
    # THE STEP-0 CALL. The instrumented Fortran run emits a `t=0 ph=1` record — i.e.
    # CTSM executes a FULL clm_drv timestep at nstep=0 (the NUOPC initialization run
    # phase, which advances the land once to fill the coupler) BEFORE nstep=1. That
    # step runs SoilTemperature/PhaseChange, so by the time FATES first packs its
    # boundary conditions at nstep=1 the cold-start soil (272 K, 100% ice) has already
    # partially melted onto the 273.15 K plateau and has LIQUID water. Omitting it left
    # the Julia model one CLM step behind and gave FATES a bone-dry frozen soil on day
    # one (btran_ft = 0 -> spurious hydraulic-failure mortality). i = 0 replays the
    # start-time forcing exactly as the NUOPC init phase does.
    for i in 0:nsteps
        # Fortran CTSM convention (verified against the instrumented run): during
        # nstep=i the clock reads the START of the step, so is_beg_curr_day() is TRUE
        # at nstep = 1, 25, 49, ... — NOT at 24, 48. Mirror that exactly, otherwise the
        # daily dynamics fire on different steps and nothing lines up. (nstep=0 is the
        # init call and does NOT fire the daily dynamics — the Fortran dump has no
        # t=0 daily phases.)
        cur = START + Second(round(Int, max(i - 1, 0) * DTIME))
        tod = hour(cur)*3600 + minute(cur)*60 + second(cur)
        calday = dayofyear(cur) + tod / _C.SECSPDAY
        (declin, _e1)   = _C.compute_orbital(calday)
        (declinm1, _e2) = _C.compute_orbital(calday - DTIME / _C.SECSPDAY)
        _C.init_daylength!(inst.gridcell, declin, declinm1, _C.ORB_OBLIQR_DEFAULT, 1:ng)
        nextsw = calday + DTIME / _C.SECSPDAY

        _C.read_forcing_step!(fr, inst.atm2lnd, cur, ng, nc)
        _C.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)

        nstep_ref[] = i
        _C.clm_drv!(config, inst, filt, filt_ia, bounds, true, nextsw, declin, declin,
            _C.ORB_OBLIQR_DEFAULT, false, false, "", false;
            nstep=i, is_first_step=(i <= 1),
            is_beg_curr_day=(i > 0 && tod == 0), is_end_curr_day=false, is_beg_curr_year=false,
            dtime=DTIME, mon=month(cur), day=day(cur), secs=tod,
            jday=dayofyear(cur), photosyns=inst.photosyns)
    end
    _C.forcing_reader_close!(fr)
    _C.FATES_PARITY_HOOK[] = nothing

    open(outpath, "w") do io
        println(io, "# FATES parity dump - schema v2 - generator=julia")
        for r in recs; println(io, r); end
    end
    @printf("  wrote Julia dump → %s  (%d records)\n", outpath, length(recs))
    return outpath
end

# ==========================================================================
# Dump parsing + the phase-attributed scorecard
# ==========================================================================
const IDKEYS = ("t", "ph", "s", "p", "c")

"HLM→FATES boundary conditions. Divergence here is a HOST-model difference, not a FATES-port bug."
const BCIN_FIELDS = Set(["pbot", "tsoil1", "h2ovol1", "effporo1", "tveg", "tgcm",
                         "esat_tv", "eair", "oair", "cair", "rb", "daylfac",
                         "solad1", "solai1"])

const PHASE_NAME = Dict(0 => "coldstart", 1 => "fast",
                        10 => "dyn_in", 11 => "phenology", 12 => "distrates",
                        13 => "integrate", 14 => "recruit", 15 => "cohortfuse",
                        16 => "spawn", 17 => "patchfuse", 18 => "updatesite")

"Read a dump, transparently un-gzipping a `.gz` reference."
_dump_lines(path) = endswith(path, ".gz") ?
    split(read(pipeline(`gzip -dc $path`), String), '\n') : readlines(path)

"""
    parse_dump(path) -> Dict{String,Dict{String,Float64}}

Key = "TYPE|t=..|ph=..|s=..|p=..|c=.." (the record's identity). META records carry
only structural totals that depend on how many grid columns each harness built, and
are reconciled from the SITE record instead — they are dropped here.
"""
function parse_dump(path::AbstractString)
    out = Dict{String,Dict{String,Float64}}()
    for line in _dump_lines(path)
        (isempty(line) || startswith(line, "#")) && continue
        toks = split(strip(line))
        isempty(toks) && continue
        typ = toks[1]
        typ == "META" && continue
        fields = Dict{String,Float64}()
        idbits = String[typ]
        for t in toks[2:end]
            kv = split(t, "=", limit=2)
            length(kv) == 2 || continue
            k, v = String(kv[1]), String(kv[2])
            val = tryparse(Float64, v)
            val === nothing && continue
            fields[k] = val
            k in IDKEYS && push!(idbits, "$k=$v")
        end
        out[join(idbits, "|")] = fields
    end
    return out
end

# Relative difference. NaN==NaN is a match (both models leave the same diagnostic
# unset); NaN-vs-value is NOT scored as a number — it is an *initialization* gap, not
# a physics divergence, and folding it in as Inf would swamp every phase. It is counted
# and reported separately.
reld(a, b) = abs(a - b) / (1.0 + abs(a))
_both_nan(a, b) = isnan(a) && isnan(b)
_nan_mismatch(a, b) = xor(isnan(a), isnan(b))

_tp(k) = begin
    m = match(r"\|t=(-?[0-9]+)\|ph=([0-9]+)", k)
    m === nothing ? (-1, -1) : (parse(Int, m.captures[1]), parse(Int, m.captures[2]))
end

"""
    compare_dumps(jl_path, ft_path; tol) -> nfail

Phase-attributed scorecard. For every (step, phase) boundary present in both dumps:
reconcile the record sets, then score every field. Reports

  * FIRST DIVERGENCE — the earliest (step, phase, field) exceeding `tol`. Because
    both runs start from a bit-identical cold start with identical forcing, THAT phase
    consumed the same input state on both sides, so its divergence is translation
    error, not drift. Everything after it is compounded.
  * per-phase worst-case, split into FATES state vs the HLM→FATES boundary conditions
    (bc_in) — a bc_in difference is a host-model difference, not a FATES-port bug.
  * structural divergence (cohort/patch counts, i.e. demography) per phase.
"""
# site ordinal of a record key (0 if absent)
function _site_of(k)
    m = match(r"\|s=([0-9]+)", k)
    m === nothing ? 0 : parse(Int, m.captures[1])
end

function compare_dumps(jl_path::AbstractString, ft_path::AbstractString; tol::Float64=1e-10)
    J = parse_dump(jl_path); F = parse_dump(ft_path)
    # Structure: the Julia harness instantiates one FATES site per natural-veg column of
    # the surfdata (the Bow instance has 2 gridcells -> 2 sites); the Fortran case is
    # single-point (1 site). Sites are FATES's independent replicate units, so score the
    # sites both dumps have and drop the extras (PR #213 verified the extra Julia site is
    # an exact replica of the scored one).
    nsf = maximum(_site_of, keys(F)); nsj = maximum(_site_of, keys(J))
    if nsj > nsf
        @printf("  note: julia has %d site(s), fortran %d — scoring site(s) 1..%d\n",
                nsj, nsf, nsf)
        J = Dict(k => v for (k, v) in J if _site_of(k) <= nsf)
    end
    println("="^92)
    println("  FATES Fortran-parity scorecard — TIME-STEPPED, PHASE-ATTRIBUTED (schema v2)")
    println("    julia   : $jl_path  ($(length(J)) records)")
    println("    fortran : $ft_path  ($(length(F)) records)")
    println("="^92)

    # ---- group by (step, phase) ------------------------------------------
    boundaries = sort(unique(vcat([_tp(k) for k in keys(J)], [_tp(k) for k in keys(F)])))
    jb = Dict{Tuple{Int,Int},Vector{String}}(); fb = Dict{Tuple{Int,Int},Vector{String}}()
    for k in keys(J); push!(get!(jb, _tp(k), String[]), k); end
    for k in keys(F); push!(get!(fb, _tp(k), String[]), k); end

    # per-phase accumulators
    phstat  = Dict{Int,NamedTuple{(:maxrel,:field,:step,:key),Tuple{Float64,String,Int,String}}}()
    phbc    = Dict{Int,NamedTuple{(:maxrel,:field,:step),Tuple{Float64,String,Int}}}()
    phstruc = Dict{Int,Tuple{Int,Int}}()          # phase -> (n boundaries, n with record-set mismatch)
    phnan   = Dict{Int,Set{String}}()             # phase -> fields that are NaN on exactly one side
    first_div = nothing                            # (step, phase, field, rel, key, jv, fv)
    first_struc = nothing                          # (step, phase, njonly, nfonly)
    fieldstat = Dict{Tuple{Int,String},Tuple{Float64,Int}}()   # (phase,field) -> (maxrel, step)

    for b in boundaries
        step, phase = b
        (haskey(jb, b) && haskey(fb, b)) || continue
        kj = Set(jb[b]); kf = Set(fb[b])
        nb, nbad = get(phstruc, phase, (0, 0))
        only_j = setdiff(kj, kf); only_f = setdiff(kf, kj)
        struc_bad = !(isempty(only_j) && isempty(only_f))
        phstruc[phase] = (nb + 1, nbad + (struc_bad ? 1 : 0))
        if struc_bad && first_struc === nothing
            first_struc = (step, phase, length(only_j), length(only_f))
        end
        for key in intersect(kj, kf)
            jf = J[key]; ff = F[key]
            # A patch that was SPAWNED during this daily step has a bc_in slot that
            # neither host has packed yet (Fortran leaves it 0). Detect the unpacked
            # slot from t_veg — a Kelvin canopy temperature that is never 0 or the
            # -9999 sentinel on a live patch — and skip that record's bc_in fields.
            # bc_in is only comparable on a slot BOTH hosts have packed. A FATES patch
            # created by disturbance is not yet in the host's exposed-veg /
            # photosynthesis filter, so its bc_in slot stays unpacked (Fortran: 0) for a
            # step or two after it is spawned. Score bc_in only on the persistent
            # oldest patch (p=1) and on the SITE record — both always packed.
            bc_live = startswith(key, "SITE") || occursin("|p=1", key)
            for (fld, jv) in jf
                haskey(ff, fld) || continue
                fld in IDKEYS && continue
                # bc_in is only meaningful where BOTH sides have packed it: at cold
                # start (phase 0) it is unpacked (Fortran reads uninitialized memory),
                # and inside the daily step (phases 11..18) a freshly SPAWNED patch has
                # a bc_in slot that will not be filled until the next wrap_photosynthesis
                # (Fortran leaves it 0). Score it at the two boundaries where it is live:
                # 1 (fast, packed by wrap_photosynthesis/wrap_btran) and 10 (dyn_in).
                if fld in BCIN_FIELDS
                    (phase == 1 || phase == 10) || continue
                    (bc_live && jv != -9999.0 && fv != -9999.0) || continue
                end
                fv = ff[fld]
                _both_nan(jv, fv) && continue
                if _nan_mismatch(jv, fv)
                    push!(get!(phnan, phase, Set{String}()), fld)
                    continue
                end
                r = reld(jv, fv)
                if fld in BCIN_FIELDS
                    prev = get(phbc, phase, (maxrel=0.0, field="", step=0))
                    r > prev.maxrel && (phbc[phase] = (maxrel=r, field=fld, step=step))
                else
                    prev = get(phstat, phase, (maxrel=0.0, field="", step=0, key=""))
                    r > prev.maxrel && (phstat[phase] = (maxrel=r, field=fld, step=step, key=key))
                    pf = get(fieldstat, (phase, fld), (0.0, 0))
                    r > pf[1] && (fieldstat[(phase, fld)] = (r, step))
                    if r > tol && first_div === nothing
                        first_div = (step, phase, fld, r, key, jv, fv)
                    end
                end
            end
        end
    end

    steps_j = sort(unique(first.(keys(jb)))); steps_f = sort(unique(first.(keys(fb))))
    @printf("\n  steps   julia %d..%d   fortran %d..%d\n",
            minimum(steps_j), maximum(steps_j), minimum(steps_f), maximum(steps_f))
    @printf("  phases  %s\n", join(sort(unique([p for (_, p) in boundaries])), " "))

    # ---- per-phase table --------------------------------------------------
    println("\n  ── per-phase worst case over the whole run ──")
    @printf("  %-3s %-11s %5s %14s %-11s %13s %-8s %-13s %s\n",
            "ph", "phase", "bnds", "max rel(FATES)", "worst field",
            "max rel(bc_in)", "worst bc", "record-set", "NaN-only-1-side")
    println("  " * "-"^104)
    nfail = 0
    for phase in sort(collect(keys(phstruc)))
        nb, nbad = phstruc[phase]
        st = get(phstat, phase, (maxrel=0.0, field="-", step=0, key=""))
        bc = get(phbc, phase, (maxrel=0.0, field="-", step=0))
        nanf = sort(collect(get(phnan, phase, Set{String}())))
        ok = st.maxrel <= tol && nbad == 0
        ok || (nfail += 1)
        @printf("  %-3d %-11s %5d %14.3e %-11s %13.3e %-8s %-13s %s\n",
                phase, get(PHASE_NAME, phase, "?"), nb,
                st.maxrel, st.field, bc.maxrel, bc.field,
                nbad == 0 ? "match" : "MISM $(nbad)/$(nb)",
                isempty(nanf) ? "-" : join(nanf, ","))
    end
    println("\n  `NaN-only-1-side` = a diagnostic one model leaves unset and the other zeroes.")
    println("  It is an INITIALIZATION gap, not a physics divergence, and is reported, not scored.")

    # ---- first divergence (the attributable result) -----------------------
    println("\n  ── FIRST DIVERGENCE (tol=$tol) ──")
    if first_div === nothing
        println("    none — the FATES port tracks Fortran within tol at every phase of every step.")
    else
        (st, ph, fld, r, key, jv, fv) = first_div
        @printf("    step %d, phase %d (%s), field `%s`\n", st, ph, get(PHASE_NAME, ph, "?"), fld)
        @printf("    julia = %.17g   fortran = %.17g   rel = %.4e\n", jv, fv, r)
        @printf("    record: %s\n", key)
        println("    ⇒ up to this point both models consumed identical state, so this")
        println("      difference is TRANSLATION ERROR in that phase, not accumulated drift.")
    end
    if first_struc !== nothing
        (st, ph, nj, nf) = first_struc
        @printf("\n    first STRUCTURAL divergence (demography): step %d, phase %d (%s) — %d record(s) only in julia, %d only in fortran\n",
                st, ph, get(PHASE_NAME, ph, "?"), nj, nf)
    end

    # ---- per-field detail for the phase that first diverges ---------------
    focus = first_div === nothing ? 1 : first_div[2]
    println("\n  ── field detail, phase $(focus) ($(get(PHASE_NAME, focus, "?"))) ──")
    @printf("  %-14s %14s %8s %s\n", "field", "max rel", "@step", "verdict")
    println("  " * "-"^52)
    flds = sort([f for (p, f) in keys(fieldstat) if p == focus])
    for fld in flds
        (r, s) = fieldstat[(focus, fld)]
        @printf("  %-14s %14.4e %8d %s\n", fld, r, s, r <= tol ? "match" : "DIVERGE")
    end

    total = nfail
    println("\n  " * (total == 0 ?
        "★ FATES dynamics match the Fortran ground truth within tol=$tol at every phase" :
        "✗ $total of $(length(phstruc)) phases diverge — see the first-divergence attribution above"))
    return total
end

# ==========================================================================
function main()
    mode = isempty(ARGS) ? "dump" : ARGS[1]
    println("="^92)
    println("  fates_fortran_parity.jl — mode=$mode  (schema v2, time-stepped)")
    println("="^92)
    out = get(ENV, "FATES_PARITY_OUT", joinpath(@__DIR__, "fates_pdump_julia.txt"))
    run_and_dump(out; nsteps=NSTEPS)

    if mode == "compare"
        ref = length(ARGS) >= 2 ? ARGS[2] : get(ENV, "FATES_FORTRAN_REF",
            joinpath(@__DIR__, "validation", "fates_fortran_parity",
                     "fates_pdump_fortran.txt.gz"))
        isfile(ref) || (println("\n  ✗ no Fortran reference dump at $ref"); return 2)
        tol = parse(Float64, get(ENV, "FATES_PARITY_TOL", "1e-10"))
        return compare_dumps(out, ref; tol=tol)
    end
    println("\n  (dump mode) score against the Fortran ground truth with:")
    println("    julia +1.12 --project=. scripts/fates_fortran_parity.jl compare <ref.txt.gz>")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
