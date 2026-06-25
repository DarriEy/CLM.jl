#!/usr/bin/env julia
# CLM.jl — Validation harness runner (Step 1 foundation).
#
# Maps each entry of test/validation/matrix.jl → build → run → oracle set →
# structured verdict, then writes results.jsonl + a markdown report. New oracle
# tiers and domains plug into the dispatch tables below (ORACLES, build_for);
# new coverage is added by appending rows to the matrix, never here.
#
# Usage:
#   julia +1.12 --project=. scripts/validation/validate.jl [--id <slug>] [--out DIR]
#   julia +1.12 --project=. scripts/validation/validate.jl            # whole matrix
#
# Data-dependent configs (data_dep=true) auto-skip with verdict :skipped_no_data
# when the machine-local domain inputs are absent (so it is CI-safe).

using CLM
using Dates
using NCDatasets
using Printf

const _HERE = @__DIR__
include(joinpath(_HERE, "..", "fortran_parity_common.jl"))   # build_bow_inst, paths, inject_dump!
include(joinpath(_HERE, "..", "..", "test", "validation", "matrix.jl"))

const WARM_RESTART = joinpath(DUMPDIR, "Bow_at_Banff_lumped.clm2.r.2003-01-01-00000.nc")

# --------------------------------------------------------------------------
# Build: config → runnable bundle. The Step-4 generalization. Two axes:
#
#   1. FLAGS — a config's mode + flags map to (a) INIT-gated flags forwarded to
#      clm_initialize! (use_cn/luna/lch4/cndv/crop/fates + use_aquifer_layer),
#      which is where their state ARRAYS are conditionally allocated (CH4 / DGVS /
#      LUNA vcmax / crop / FATES), and (b) CONFIG-only flags set on CLMDriverConfig
#      (hydrstress/voc/ozone/c13/c14/matrixcn/irrigate) — those state arrays are
#      ALWAYS allocated in clm_initialize!, so they only need the driver toggle.
#      Splitting on this is the keystone: a swept flag now gets a well-posed
#      instance with its state allocated, not a toggle over un-allocated state.
#
#   2. DOMAIN — :bow uses the calibrated Fortran-parity builder (build_bow_inst);
#      every other domain uses the generic build_domain_inst, mirroring the proven
#      run_clm_streamflow.jl path (clm_initialize! + Jackson rooting + baseflow)
#      over the Symfluence domain layout. Domains whose inputs are absent skip.
# --------------------------------------------------------------------------
struct RunBundle
    inst; bounds; filt; tm; config; fr; dtime::Float64
end

# Symfluence domain layout (shared with scripts/run_clm_streamflow.jl).
const SYMROOT = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data"
# Validation-domain symbol → on-disk domain dir name. :bow is special-cased
# (calibrated parity paths in fortran_parity_common.jl), so it is not listed here.
const DOMAIN_DIRNAME = Dict(
    :aripuana   => "Aripuana_Amazon",
    :stillwater => "Stillwater_Oklahoma",
    :krycklan   => "Boreal_Krycklan_Sweden",
    :abisko     => "Arctic_Abisko_Sweden",
    :tagus      => "Mediterranean_Tagus_Spain",
    :massa      => "Alps_Massa_Aletsch_CH",
    :baltimore  => "Urban_DeadRun_Baltimore",
    :iceland    => "Iceland_Jokulsa_Fjollum",
)
_domain_dir(dom)  = joinpath(SYMROOT, "domain_$(DOMAIN_DIRNAME[dom])")
_domain_surf(dom) = joinpath(_domain_dir(dom), "settings", "CLM", "parameters", "surfdata_clm.nc")
_domain_parm(dom) = joinpath(_domain_dir(dom), "settings", "CLM", "parameters", "clm5_params.nc")

"First clmforc.*.nc forcing file for a domain (the most-merged if several), or nothing."
function _domain_forcing(dom)
    cdir = joinpath(_domain_dir(dom), "data", "forcing", "CLM_input")
    isdir(cdir) || return nothing
    files = filter(f -> startswith(f, "clmforc") && endswith(f, ".nc"), readdir(cdir))
    isempty(files) && return nothing
    joinpath(cdir, files[argmax(length.(files))])   # prefer the most-merged span
end

"Start the run at the forcing file's first timestamp so read_forcing_step! is in-coverage."
function _forcing_start(forcing_path)
    ds = NCDataset(forcing_path, "r")
    t = nothing
    for nm in ("time", "DTIME", "datetime")
        haskey(ds, nm) && (t = ds[nm][1]; break)
    end
    close(ds)
    t isa DateTime ? t : DateTime(2003,1,1)
end

"Init-gated flags for clm_initialize! — these allocate state arrays when true."
function _init_flags(cfg)
    f = cfg.flags
    (; use_cn   = cfg.mode === :cn,
       use_fates = cfg.mode in (:fates_sp, :fates_bgc),
       use_luna = get(f, :use_luna, false),
       use_lch4 = get(f, :use_lch4, false),
       use_cndv = get(f, :use_cndv, false),
       use_crop = get(f, :use_crop, false),
       use_aquifer_layer = get(f, :use_aquifer_layer, false))
end

"Does this config's machine-local input data exist?"
function data_available(cfg)
    if cfg.domain === :bow
        return isfile(FSURDAT) && isfile(FPARAM) && isfile(FFORCING)
    end
    haskey(DOMAIN_DIRNAME, cfg.domain) || return false
    isfile(_domain_surf(cfg.domain)) && isfile(_domain_parm(cfg.domain)) &&
        _domain_forcing(cfg.domain) !== nothing
end

"""Construct the CLMDriverConfig for a config's mode + flags. Forwards the init-gated
flags (so the driver dispatch matches the allocated state) AND the config-only flags
(hydrstress/voc/ozone/c13/c14/matrixcn/irrigate — always-allocated state)."""
function _make_config(cfg)
    f = cfg.flags
    kw = Dict{Symbol,Any}(:use_aquifer_layer => get(f, :use_aquifer_layer, false))
    if cfg.mode === :cn
        kw[:use_cn] = true
    elseif cfg.mode === :fates_sp || cfg.mode === :fates_bgc
        kw[:use_fates] = true
    end
    for k in (:use_luna, :use_hydrstress, :use_voc, :use_ozone, :irrigate,
              :use_lch4, :use_c13, :use_c14, :use_crop, :use_cndv, :use_matrixcn)
        haskey(f, k) && (kw[k] = getproperty(f, k))
    end
    CLM.CLMDriverConfig(; kw...)
end

"Generic non-Bow domain builder — mirrors run_clm_streamflow.jl's clm_initialize! path."
function build_domain_inst(cfg, fl)
    dom = cfg.domain
    fsurdat = _domain_surf(dom); paramfile = _domain_parm(dom)
    forcing = _domain_forcing(dom)
    start_date = _forcing_start(forcing)
    # Bow lnd_in convention: Jackson-1996 rooting (matches the streamflow runner +
    # parity harness). Set before clm_initialize!'s cold-start init_vegrootfr!.
    CLM.rooting_profile_config.rooting_profile_method_water  = CLM.JACKSON_1996_ROOT
    CLM.rooting_profile_config.rooting_profile_method_carbon = CLM.JACKSON_1996_ROOT
    (inst, bounds, filt, tm) = CLM.clm_initialize!(;
        fsurdat=fsurdat, paramfile=paramfile, start_date=start_date, dtime=3600,
        use_cn=fl.use_cn, use_luna=fl.use_luna, use_lch4=fl.use_lch4,
        use_cndv=fl.use_cndv, use_crop=fl.use_crop, use_fates=fl.use_fates,
        use_bedrock=true, use_aquifer_layer=fl.use_aquifer_layer, h2osfcflag=0,
        fsnowoptics=FSNOWOPT, fsnowaging=FSNOWAGE, int_snow_max=2000.0)
    CLM.init_soil_hydrology_config(baseflow_scalar=BASEFLOW_SCALAR)
    fr = CLM.ForcingReader(); CLM.forcing_reader_init!(fr, forcing)
    return RunBundle(inst, bounds, filt, tm, _make_config(cfg), fr, 3600.0)
end

function build_for(cfg)::Union{RunBundle,Nothing}
    fl = _init_flags(cfg)
    if cfg.domain === :bow
        dtime = 3600.0
        (inst, bounds, filt, tm) = build_bow_inst(; dtime=Int(dtime),
            start_date=DateTime(2003,1,1),
            use_cn=fl.use_cn, use_luna=fl.use_luna, use_lch4=fl.use_lch4,
            use_cndv=fl.use_cndv, use_crop=fl.use_crop, use_fates=fl.use_fates)
        if cfg.init === :warm && isfile(WARM_RESTART)
            inject_dump!(inst, bounds, WARM_RESTART)
            # PHS (use_hydrstress) is only well-posed from a seeded vegwp: cold-start
            # leaves the Newton solve to diverge to NaN canopy fluxes (the documented
            # coldstart-canopy-nan gap). Seed vegwp from the restart so the solve
            # starts in-basin — mirrors the vcmx25 LUNA seeding in run_one_parity_step!.
            if get(cfg.flags, :use_hydrstress, false)
                ds = NCDataset(WARM_RESTART, "r")
                if haskey(ds, "vegwp")
                    vw = ds["vegwp"][:, :]            # (vegwcs, pft)
                    for pd in 1:size(vw, 2), seg in 1:min(4, size(vw, 1))
                        pd <= size(inst.canopystate.vegwp_patch, 1) &&
                            (inst.canopystate.vegwp_patch[pd, seg] = Float64(vw[seg, pd]))
                    end
                end
                close(ds)
            end
        end
        fr = CLM.ForcingReader(); CLM.forcing_reader_init!(fr, FFORCING)
        return RunBundle(inst, bounds, filt, tm, _make_config(cfg), fr, dtime)
    end
    haskey(DOMAIN_DIRNAME, cfg.domain) || return nothing   # unknown domain → unwired
    data_available(cfg) || return nothing                  # inputs absent → skip
    return build_domain_inst(cfg, fl)
end

# --------------------------------------------------------------------------
# Run engine: advance N steps (canonical clm_drv! loop), per-step finiteness.
# --------------------------------------------------------------------------
_nonfinite(arr) = count(x -> !(isfinite(x)), arr)

function _state_nonfinite(inst)
    t = inst.temperature
    ws = inst.water.waterstatebulk_inst.ws
    _nonfinite(t.t_grnd_col) + _nonfinite(t.t_soisno_col) +
        _nonfinite(ws.h2osoi_liq_col) + _nonfinite(ws.h2osoi_ice_col)
end

"A flat snapshot of the core prognostic state — the determinism/parity fingerprint."
function _snapshot(inst)
    t = inst.temperature
    ws = inst.water.waterstatebulk_inst.ws
    vcat(vec(copy(t.t_grnd_col)), vec(copy(t.t_soisno_col)),
         vec(copy(ws.h2osoi_liq_col)), vec(copy(ws.h2osoi_ice_col)),
         vec(copy(ws.wa_col)))
end

"""
    run_steps!(b, nsteps) -> (steps_done, nonfinite, err)

Run the canonical clm_drv! loop for `nsteps`. clm_drv! itself error()s on a
water-balance violation (errh2o > 1e-5 mm), so a clean return of `nsteps` with
`nonfinite==0` IS the T2 conservation pass. Any thrown error (balance or NaN
propagation) is captured and ends the run early.
"""
function run_steps!(b::RunBundle, nsteps::Int)
    ng, nc, np = b.bounds.endg, b.bounds.endc, b.bounds.endp
    filt_ia = CLM.clump_filter_inactive_and_active
    obliqr = CLM.ORB_OBLIQR_DEFAULT
    nonfinite = 0
    done = 0
    err = nothing
    try
        for s in 1:nsteps
            step_start = b.tm.current_date
            calday = CLM.get_curr_calday(b.tm)
            (declin, _) = CLM.compute_orbital(calday)
            nextsw = calday + b.dtime / CLM.SECSPDAY
            (declinm1, _) = CLM.compute_orbital(calday - b.dtime / CLM.SECSPDAY)
            CLM.init_daylength!(b.inst.gridcell, declin, declinm1, obliqr, 1:ng)
            CLM.advance_timestep!(b.tm)
            CLM.read_forcing_step!(b.fr, b.inst.atm2lnd, step_start, ng, nc)
            CLM.downscale_forcings!(b.bounds, b.inst.atm2lnd, b.inst.column,
                                    b.inst.landunit, b.inst.topo)
            (yr, mon, d, tod) = CLM.get_curr_date(b.tm)
            CLM.clm_drv!(b.config, b.inst, b.filt, filt_ia, b.bounds,
                         true, nextsw, declin, declin, obliqr,
                         false, false, "", false;
                         nstep=b.tm.nstep, is_first_step=(b.tm.nstep==1),
                         is_beg_curr_day=CLM.is_beg_curr_day(b.tm),
                         is_end_curr_day=CLM.is_end_curr_day(b.tm),
                         is_beg_curr_year=CLM.is_beg_curr_year(b.tm),
                         dtime=b.dtime, mon=mon, day=d, photosyns=b.inst.photosyns)
            CLM.lnd2atm!(b.bounds, b.inst)
            done = s
            nf = _state_nonfinite(b.inst)
            nf > nonfinite && (nonfinite = nf)
            nf > 0 && (err = "non-finite state ($nf) at step $s"; break)
        end
    catch e
        err = sprint(showerror, e)
    finally
        CLM.forcing_reader_close!(b.fr)
    end
    return (done, nonfinite, err, _snapshot(b.inst))
end

# --------------------------------------------------------------------------
# Oracle dispatch. Each oracle takes `cfg` and builds the instance(s) it needs
# (determinism needs two runs), so they compose without sharing mutated state.
# Steps 2–7 add oracle functions here; unwired oracles report pass=missing.
# --------------------------------------------------------------------------
function oracle_conservation(cfg)
    nsteps = DEPTH_DEFAULT_STEPS[cfg.depth]
    (done, nf, err, _) = run_steps!(build_for(cfg), nsteps)
    pass = (done == nsteps) && (nf == 0) && (err === nothing)
    (; oracle=:conservation, pass,
       metrics=(; steps_target=nsteps, steps_done=done, nonfinite=nf),
       detail = err === nothing ? "ran $done/$nsteps steps, balance-clean, finite" : err)
end

"""
T3 determinism: the same config run twice from a fresh build must produce a
bit-identical final-state fingerprint. Catches nondeterminism (uninitialized
reads, hash/iteration-order, RNG leakage) that conservation alone misses.
"""
function oracle_determinism(cfg)
    nsteps = DEPTH_DEFAULT_STEPS[cfg.depth]
    (d1, _, e1, s1) = run_steps!(build_for(cfg), nsteps)
    (d2, _, e2, s2) = run_steps!(build_for(cfg), nsteps)
    ran = (e1 === nothing) && (e2 === nothing) && d1 == nsteps && d2 == nsteps
    identical = ran && length(s1) == length(s2) && all(s1 .=== s2)   # === so NaN==NaN
    (; oracle=:determinism, pass=identical,
       metrics=(; ran_both=ran, maxabsdiff = ran ? maximum(abs.(s1 .- s2); init=0.0) : NaN),
       detail = ran ? (identical ? "two runs bit-identical" :
                       "DIVERGED max|Δ|=$(maximum(abs.(s1 .- s2); init=0.0))") :
                      "a run failed: $(something(e1, e2, "early stop"))")
end

const ORACLES = Dict{Symbol,Function}(
    :conservation => oracle_conservation,
    :determinism  => oracle_determinism,
)

# --------------------------------------------------------------------------
# Validate one config → verdict NamedTuple.
# --------------------------------------------------------------------------
function validate_one(cfg)
    t0 = time()
    if cfg.data_dep && !data_available(cfg)
        return (; cfg.id, status=:skipped_no_data, oracles=Any[],
                wall=0.0, note="machine-local data for domain $(cfg.domain) absent")
    end
    if build_for(cfg) === nothing
        return (; cfg.id, status=:skipped_unwired, oracles=Any[],
                wall=0.0, note="domain $(cfg.domain) not wired in runner yet")
    end
    results = Any[]
    for o in cfg.oracles
        if haskey(ORACLES, o)
            push!(results, ORACLES[o](cfg))   # each oracle builds its own instance(s)
        else
            push!(results, (; oracle=o, pass=missing, metrics=(;),
                            detail="oracle :$o not yet wired"))
        end
    end
    allpass = all(r -> r.pass === true, results)
    anyfail = any(r -> r.pass === false, results)
    status = anyfail ? :fail : (allpass ? :pass : :partial)
    (; cfg.id, status, oracles=results, wall=round(time()-t0, digits=2), note=cfg.note)
end

# --------------------------------------------------------------------------
# Reporting.
# --------------------------------------------------------------------------
_jstr(s) = '"' * replace(String(s), '\\'=>"\\\\", '"'=>"\\\"", '\n'=>"\\n") * '"'

function write_results(verdicts, outdir)
    mkpath(outdir)
    open(joinpath(outdir, "results.jsonl"), "w") do io
        for v in verdicts
            orac = join(["{" * _jstr("oracle") * ":" * _jstr(r.oracle) * "," *
                         _jstr("pass") * ":" * string(r.pass) * "}" for r in v.oracles], ",")
            println(io, "{", _jstr("id"), ":", _jstr(v.id), ",",
                    _jstr("status"), ":", _jstr(v.status), ",",
                    _jstr("wall"), ":", string(v.wall), ",",
                    _jstr("oracles"), ":[", orac, "]}")
        end
    end
    open(joinpath(outdir, "report.md"), "w") do io
        n = length(verdicts)
        np = count(v -> v.status === :pass, verdicts)
        nf = count(v -> v.status === :fail, verdicts)
        ns = count(v -> v.status in (:skipped_no_data, :skipped_unwired), verdicts)
        println(io, "# CLM.jl validation report\n")
        println(io, "$np passed / $nf failed / $ns skipped / $n total\n")
        println(io, "| id | status | oracles | wall (s) |")
        println(io, "|----|--------|---------|----------|")
        for v in verdicts
            oc = join(["$(r.oracle)=$(r.pass)" for r in v.oracles], "; ")
            println(io, "| $(v.id) | $(v.status) | $oc | $(v.wall) |")
        end
    end
end

# --------------------------------------------------------------------------
# Main.
# --------------------------------------------------------------------------
function main(args)
    only_id = nothing
    outdir = joinpath(_HERE, "..", "..", "validation_results")
    i = 1
    while i <= length(args)
        if args[i] == "--id"; only_id = args[i+1]; i += 2
        elseif args[i] == "--out"; outdir = args[i+1]; i += 2
        else; i += 1; end
    end
    M = validation_matrix()
    only_id !== nothing && (M = filter(c -> c.id == only_id, M))
    isempty(M) && (println("no matching configs"); return)
    println("Running $(length(M)) validation config(s)...")
    verdicts = Any[]
    for cfg in M
        v = validate_one(cfg)
        push!(verdicts, v)
        @printf("  %-22s %-18s %5.1fs\n", v.id, v.status, v.wall)
    end
    write_results(verdicts, outdir)
    np = count(v -> v.status === :pass, verdicts)
    nf = count(v -> v.status === :fail, verdicts)
    println("\n$np passed, $nf failed → $(joinpath(outdir, "report.md"))")
    return nf == 0
end

if abspath(PROGRAM_FILE) == (@__FILE__)
    main(ARGS)
end
