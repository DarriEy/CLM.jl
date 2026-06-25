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
# Build: config → runnable bundle. Step 1 wires the :bow domain (cold or warm).
# Other domains return `nothing` → the runner records :skipped_unwired.
# --------------------------------------------------------------------------
struct RunBundle
    inst; bounds; filt; tm; config; fr; dtime::Float64
end

"Does this config's machine-local input data exist?"
function data_available(cfg)
    cfg.domain === :bow || return false           # only Bow wired in Step 1
    isfile(FSURDAT) && isfile(FPARAM) && isfile(FFORCING)
end

"Construct the CLMDriverConfig for a config's mode + flags."
function _make_config(cfg)
    f = cfg.flags
    kw = Dict{Symbol,Any}(:use_aquifer_layer => get(f, :use_aquifer_layer, false))
    if cfg.mode === :cn
        kw[:use_cn] = true
    elseif cfg.mode === :fates_sp || cfg.mode === :fates_bgc
        kw[:use_fates] = true
    end
    # thread through any recognised driver-level use_* flags present in cfg.flags
    for k in (:use_luna, :use_hydrstress, :use_voc, :use_ozone, :irrigate,
              :use_lch4, :use_c13, :use_c14, :use_crop, :use_cndv, :use_matrixcn)
        haskey(f, k) && (kw[k] = getproperty(f, k))
    end
    CLM.CLMDriverConfig(; kw...)
end

function build_for(cfg)::Union{RunBundle,Nothing}
    cfg.domain === :bow || return nothing
    use_cn = cfg.mode === :cn
    use_luna = get(cfg.flags, :use_luna, false)
    dtime = 3600.0
    (inst, bounds, filt, tm) = build_bow_inst(; dtime=Int(dtime),
                                              start_date=DateTime(2003,1,1),
                                              use_cn=use_cn, use_luna=use_luna)
    if cfg.init === :warm && isfile(WARM_RESTART)
        inject_dump!(inst, bounds, WARM_RESTART)
    end
    config = _make_config(cfg)
    fr = CLM.ForcingReader(); CLM.forcing_reader_init!(fr, FFORCING)
    return RunBundle(inst, bounds, filt, tm, config, fr, dtime)
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
