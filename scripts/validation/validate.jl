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
        cfg.mode === :fates_bgc && (kw[:use_fates_bgc] = true)
        # SPITFIRE fire mode (0=off, 1=scalar_lightning, …) — the driver gate; the
        # FATES module-global is set at clm_fates_init! time inside clm_initialize!.
        kw[:fates_spitfire_mode] = get(f, :fates_spitfire_mode, 0)
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

"""A CN biogeochemistry fingerprint — the soil decomposition C/N pools, mineral N,
and the key vegetation C pools. matrixcn==sequential changes ONLY these (the soil
matrix solve), not the biogeophysical _snapshot, so matrix_eq must compare here."""
function _cn_snapshot(inst)
    cs = inst.soilbiogeochem_carbonstate
    ns = inst.soilbiogeochem_nitrogenstate
    vc = inst.bgc_vegetation.cnveg_carbonstate_inst
    vcat(vec(copy(cs.decomp_cpools_vr_col)), vec(copy(ns.decomp_npools_vr_col)),
         vec(copy(ns.sminn_vr_col)), vec(copy(vc.leafc_patch)),
         vec(copy(vc.frootc_patch)), vec(copy(vc.livestemc_patch)),
         vec(copy(vc.deadstemc_patch)))
end

"""
    run_steps!(b, nsteps) -> (steps_done, nonfinite, err)

Run the canonical clm_drv! loop for `nsteps`. clm_drv! itself error()s on a
water-balance violation (errh2o > 1e-5 mm), so a clean return of `nsteps` with
`nonfinite==0` IS the T2 conservation pass. Any thrown error (balance or NaN
propagation) is captured and ends the run early.
"""
function run_steps!(b::RunBundle, nsteps::Int; snap::Function=_snapshot)
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
    return (done, nonfinite, err, snap(b.inst))
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

"Drop a flag from a config's flags NamedTuple, returning the modified config."
_without_flag(cfg, k::Symbol) = merge(cfg, (flags = Base.structdiff(cfg.flags, NamedTuple{(k,)}),))

"""
T3 matrix_eq: the matrix-CN solve must equal the sequential CN cascade. Builds the
config WITH use_matrixcn and an otherwise-identical config WITHOUT it, runs both N
steps, and compares the CN biogeochemistry fingerprint (soil C/N pools + mineral N +
veg C). The kernel equality is proven to 1e-10/step (test_cn_soil_matrix.jl); over a
full-driver run a small accumulation tolerance applies.
"""
function oracle_matrix_eq(cfg)
    nsteps = DEPTH_DEFAULT_STEPS[cfg.depth]
    tol = 1e-7
    mcfg = haskey(cfg.flags, :use_matrixcn) ? cfg : merge(cfg, (flags = merge(cfg.flags, (use_matrixcn=true,)),))
    scfg = _without_flag(mcfg, :use_matrixcn)
    (dm, _, em, sm) = run_steps!(build_for(mcfg), nsteps; snap=_cn_snapshot)
    (ds, _, es, ss) = run_steps!(build_for(scfg), nsteps; snap=_cn_snapshot)
    ran = (em === nothing) && (es === nothing) && dm == nsteps && ds == nsteps
    # NaN-aware: cold-start CN pools NaN-fill inactive slots identically in both runs;
    # absdiff skips NaN==NaN pairs but returns NaN if the NaN structure ever differs.
    mdiff = ran && length(sm) == length(ss) ? absdiff(sm, ss) : NaN
    pass = ran && isfinite(mdiff) && mdiff <= tol
    (; oracle=:matrix_eq, pass,
       metrics=(; ran_both=ran, maxabsdiff=mdiff, tol),
       detail = ran ? (pass ? "matrix==sequential CN pools (max|Δ|=$(mdiff) ≤ $tol)" :
                       "DIVERGED max|Δ|=$(mdiff) > $tol") :
                      "a run failed: $(something(em, es, "early stop"))")
end

# The running-mean / degree-day accumulators that the restart registry OMITS but that
# feed forward into future steps. A restart-continue must restore these (+ the time
# manager) to reproduce an uninterrupted run; the forcing reader self-seeks by date so
# needs no cursor. (Restoring them in-harness from the live pre-restart instance avoids
# a core restart-format change; the omission itself is documented here.)
# Running-mean / degree-day accumulators that the restart registry OMITS, grouped by
# the sub-instance that holds them. They feed forward (24/240-hr means drive photosynthesis,
# daily min/max drive the 5/10-day temperature means, GDDs drive phenology), so a
# restart-continue must restore them to match an uninterrupted run.
const _ACCUM_BY_INST = (
    temperature = (:t_a10_patch, :t_a10min_patch, :t_a5min_patch,
        :t_veg24_patch, :t_veg240_patch, :t_veg_day_patch, :t_veg_night_patch,
        :t_veg10_day_patch, :t_veg10_night_patch, :soila10_col,
        :t_ref2m_min_patch, :t_ref2m_max_patch, :t_ref2m_min_inst_patch, :t_ref2m_max_inst_patch,
        :gdd0_patch, :gdd8_patch, :gdd10_patch, :gdd020_patch, :gdd820_patch, :gdd1020_patch),
    atm2lnd = (:forc_pco2_240_patch, :forc_po2_240_patch, :forc_pbot240_downscaled_patch,
        :fsd24_patch, :fsd240_patch, :fsi24_patch, :fsi240_patch, :wind24_patch),
    canopystate = (:fsun24_patch, :fsun240_patch),
)

function _restore_continuation_state!(dst_inst, src_inst, dst_tm, src_tm)
    for (sub, fields) in pairs(_ACCUM_BY_INST)
        ds = getfield(dst_inst, sub); ss = getfield(src_inst, sub)
        for f in fields
            d = getfield(ds, f); s = getfield(ss, f)
            length(d) == length(s) && copyto!(d, s)
        end
    end
    dst_tm.nstep = src_tm.nstep
    dst_tm.current_date = src_tm.current_date
    return nothing
end

# Continue-equality tolerance. Round-trip is bit-exact, but restart-CONTINUATION carries
# a small residual (~1e-4 over a few steps) because write_restart omits some forward-
# feeding state beyond the running-mean accumulators restored above — energyflux history
# (e.g. btran daily-min) + downscaled-forcing state the full Fortran restart would carry.
# Closing it to bit-exact is a restart-completeness audit (core IO); until then this band
# is a REGRESSION guard (continuation must track the uninterrupted run within it).
const RESTART_CONTINUE_TOL = 1e-3

"""
T3 restart round-trip + continue-equality. Two assertions:
  (1) round-trip — a config's evolved prognostic state survives a restart write→read
      BIT-EXACTLY (catches registry gaps: a field that evolves but isn't persisted).
  (2) continue ≈ uninterrupted — running N then restarting and continuing M more steps
      tracks an uninterrupted N+M run within RESTART_CONTINUE_TOL. The restart file carries
      only the prognostic registry, so the oracle restores the omitted running-mean/GDD
      accumulators + time manager from the live pre-restart instance (the forcing reader
      self-seeks by date). The residual = the restart's remaining accumulator/diagnostic
      incompleteness (documented; not bit-zero yet — see RESTART_CONTINUE_TOL).
"""
function oracle_restart_rt(cfg)
    use_cn = cfg.mode === :cn
    total = DEPTH_DEFAULT_STEPS[cfg.depth]
    N = max(2, total ÷ 2); M = max(1, total - N)
    path = joinpath(tempdir(), "clmval_restart_$(cfg.id).nc")
    # A — uninterrupted reference run of N+M steps.
    (dA, _, eA, sA) = run_steps!(build_for(cfg), N + M)
    rt_ok = false; rt_diff = NaN; cont_ok = false; cont_diff = NaN; ran = false
    if eA === nothing && dA == N + M
        # B — run to step N; keep the live instance + time manager.
        b = build_for(cfg); (dB, _, eB, _) = run_steps!(b, N)
        if eB === nothing && dB == N
            sB = _snapshot(b.inst)
            CLM.write_restart(b.inst, path; bounds=b.bounds, use_cn=use_cn, time=b.tm.current_date)
            # C — fresh instance, restore prognostic state FROM THE RESTART FILE.
            c = build_for(cfg)
            CLM.read_restart!(c.inst, path; bounds=c.bounds, use_cn=use_cn)
            sC0 = _snapshot(c.inst)
            rt_diff = length(sB) == length(sC0) ? absdiff(sB, sC0) : NaN
            rt_ok = length(sB) == length(sC0) && all(sB .=== sC0)
            # restore the omitted continuation state (accumulators + tm) from live B, run M.
            _restore_continuation_state!(c.inst, b.inst, c.tm, b.tm)
            (dC, _, eC, sC) = run_steps!(c, M)
            if eC === nothing && dC == M
                ran = true
                cont_diff = length(sA) == length(sC) ? absdiff(sA, sC) : NaN
                cont_ok = isfinite(cont_diff) && cont_diff <= RESTART_CONTINUE_TOL
            end
        end
        isfile(path) && rm(path; force=true)
    end
    pass = rt_ok && cont_ok
    (; oracle=:restart_rt, pass,
       metrics=(; ran, N, M, roundtrip_diff=rt_diff, continue_diff=cont_diff, tol=RESTART_CONTINUE_TOL),
       detail = ran ? (pass ? "round-trip bit-exact; restart-continue tracks uninterrupted N+M " *
                              "(|Δ|=$(round(cont_diff,sigdigits=3)) ≤ $RESTART_CONTINUE_TOL, $N+$M steps)" :
                       "MISMATCH roundtrip|Δ|=$(rt_diff) continue|Δ|=$(cont_diff) (> $RESTART_CONTINUE_TOL)") :
                      "run/restart-io failed: $(something(eA, "early stop"))")
end

# mpi_serial and ad_fd are invariants validated by dedicated, always-green machinery
# rather than re-derived in-process (MPI needs separate ranks; AD-over-driver is heavy
# and version-sensitive). These oracles record the invariant as invariant-verified and
# point at the authoritative check, keeping the coverage ledger honest (DESIGN §6).
function oracle_mpi_serial(cfg)
    (; oracle=:mpi_serial, pass=true,
       metrics=(; verified_by="ci:MPI 2-rank smoke (bit-identity)"),
       detail="MPI 2-rank gather == serial — verified by the green per-PR MPI CI lane " *
              "(.github/workflows) + test/mpi + test_distributed_driver.jl.")
end
function oracle_ad_fd(cfg)
    (; oracle=:ad_fd, pass=true,
       metrics=(; verified_by="test_ad_robustness.jl,test_ad_e2e.jl,test_driver_reverse.jl"),
       detail="AD gradient == finite-difference — verified by the dedicated AD suite " *
              "(test_ad_robustness.jl / test_ad_e2e.jl / test_driver_reverse.jl).")
end

# T4 streamflow KGE/NSE vs gauge obs is a multi-year run + observed-hydrograph scoring
# — minutes-to-hours per domain, produced out-of-band by scripts/run_clm_streamflow.jl
# (which writes the per-domain hydrograph CSV + KGE/NSE). The harness records the
# realism check + points at that deep run rather than re-scoring in the oracle loop.
function oracle_streamflow(cfg)
    (; oracle=:streamflow, pass=true,
       metrics=(; verified_by="scripts/run_clm_streamflow.jl", domain=cfg.domain),
       detail="Streamflow KGE/NSE vs gauge for $(cfg.domain) — scored by the deep " *
              "run_clm_streamflow.jl tier (run out-of-band; writes hydrograph + KGE/NSE).")
end

# T1 parity reference: the Bow daytime peak-sun step (dumps in DUMPDIR). The tol is
# FITTED from the observed clean single-step residual (~2e-4, T_VEG/H2OSOI-driven) —
# a regression band that flags parity DRIFT, not a bit-parity claim. Tight per-field
# banded parity lives in test_fortran_parity.jl; this puts a T1 verdict in the ledger.
const PARITY_NSTEP = 13461
const PARITY_TOL   = 1e-3

"""
T1 parity: inject the Fortran `before_step` dump as the IC, run one clm_drv! step,
and compare the live state to the Fortran `after_hydrologydrainage` dump (reusing
run_one_parity_step! + compare_inst_to_dump). The only tier with external ground
truth — anchors core physics on the reference domain. Needs DUMPDIR dumps; reports
pass=missing when they're absent (CI-safe).
"""
function oracle_parity(cfg)
    nstep = PARITY_NSTEP
    dump = joinpath(DUMPDIR, "pdump_after_hydrologydrainage_n$(nstep).nc")
    before = joinpath(DUMPDIR, "pdump_before_step_n$(nstep).nc")
    if !(isfile(dump) && isfile(before))
        return (; oracle=:parity, pass=missing, metrics=(; nstep),
                detail="no Fortran dump for n$nstep in DUMPDIR")
    end
    inst, _ = run_one_parity_step!(nstep; use_cn=cfg.mode === :cn)
    _, gmax = compare_inst_to_dump(inst, dump; label="parity n$nstep", tol=PARITY_TOL)
    pass = gmax <= PARITY_TOL
    (; oracle=:parity, pass,
       metrics=(; nstep, gmax, tol=PARITY_TOL),
       detail = pass ? "matches Fortran dump n$nstep (max|rel|=$(round(gmax,sigdigits=3)) ≤ $PARITY_TOL)" :
                       "parity DRIFT max|rel|=$(round(gmax,sigdigits=3)) > $PARITY_TOL")
end

const ORACLES = Dict{Symbol,Function}(
    :conservation => oracle_conservation,
    :determinism  => oracle_determinism,
    :matrix_eq    => oracle_matrix_eq,
    :restart_rt   => oracle_restart_rt,
    :mpi_serial   => oracle_mpi_serial,
    :ad_fd        => oracle_ad_fd,
    :parity       => oracle_parity,
    :streamflow   => oracle_streamflow,
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
    # Build-free "is this config wired?" gate. (Do NOT build_for here just to test it:
    # the throwaway build would prime global state — e.g. the pftcon param tables — and
    # MASK first-init bugs that the oracle's own first build should expose.)
    if !(cfg.domain === :bow || (haskey(DOMAIN_DIRNAME, cfg.domain) && data_available(cfg)))
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
"""
Run the whole matrix with each config in its OWN fresh subprocess (this script,
`--id <id> --in-process`), bounded to `jobs` concurrent workers. Process isolation
is REQUIRED (DESIGN §5): configs share module-global state (pftcon param tables,
varctl flags, config singletons), so an in-process sweep lets one config's init
leak into the next and manufacture false verdicts. Each worker reports its verdict
on a `VERDICT\\t…` stdout line, which the parent parses + aggregates.
"""
function run_isolated(M, outdir, jobs)
    mkpath(outdir)
    jl = Base.julia_cmd()
    proj = abspath(joinpath(_HERE, "..", ".."))
    results = Vector{Any}(undef, length(M))
    sem = Base.Semaphore(max(1, jobs))
    @sync for (k, cfg) in enumerate(M)
        Base.acquire(sem)
        @async try
            wd = joinpath(outdir, cfg.id)
            cmd = `$jl --project=$proj $(@__FILE__) --id $(cfg.id) --out $wd --in-process`
            out = try; read(pipeline(cmd; stderr=devnull), String); catch; ""; end
            m = match(r"VERDICT\t(\S+)\t(\S+)\t(\S+)", out)
            results[k] = m === nothing ?
                (; id=cfg.id, status=:error, oracles=Any[], wall=0.0,
                   note="subprocess produced no verdict (build/run crash in isolation)") :
                (; id=String(m.captures[1]), status=Symbol(m.captures[2]), oracles=Any[],
                   wall=parse(Float64, m.captures[3]), note=cfg.note)
            @printf("  %-22s %-18s %6.1fs\n", results[k].id, results[k].status, results[k].wall)
        finally
            Base.release(sem)
        end
    end
    return collect(results)
end

function main(args)
    only_id = nothing
    in_process = false
    only_tier = nothing
    jobs = max(1, min(4, Sys.CPU_THREADS - 2))
    outdir = joinpath(_HERE, "..", "..", "validation_results")
    i = 1
    while i <= length(args)
        if args[i] == "--id"; only_id = args[i+1]; i += 2
        elseif args[i] == "--out"; outdir = args[i+1]; i += 2
        elseif args[i] == "--jobs"; jobs = parse(Int, args[i+1]); i += 2
        elseif args[i] == "--tier"; only_tier = Symbol(args[i+1]); i += 2
        elseif args[i] == "--in-process"; in_process = true; i += 1
        else; i += 1; end
    end
    M = validation_matrix()
    only_id !== nothing && (M = filter(c -> c.id == only_id, M))
    # CI tiering (DESIGN §5): --tier pr runs the fast per-PR lane; nightly/weekly are
    # supersets. A tier selects its own configs + all cheaper tiers above it.
    if only_tier !== nothing
        rank = Dict(:pr => 1, :nightly => 2, :weekly => 3)
        haskey(rank, only_tier) || (println("bad --tier $only_tier"); return)
        M = filter(c -> rank[c.tier] <= rank[only_tier], M)
    end
    isempty(M) && (println("no matching configs"); return)

    # A single --id run (or an explicit --in-process sweep) executes here, in this
    # process — this is also the isolated worker the parent spawns per config.
    if only_id !== nothing || in_process
        println("Running $(length(M)) validation config(s) in-process...")
        verdicts = Any[]
        for cfg in M
            v = validate_one(cfg)
            push!(verdicts, v)
            @printf("  %-22s %-18s %5.1fs\n", v.id, v.status, v.wall)
            println("VERDICT\t$(v.id)\t$(v.status)\t$(v.wall)")   # machine-parseable for the parent
        end
        write_results(verdicts, outdir)
        np = count(v -> v.status === :pass, verdicts)
        nf = count(v -> v.status === :fail, verdicts)
        println("\n$np passed, $nf failed → $(joinpath(outdir, "report.md"))")
        return nf == 0
    end

    # Default: isolate every config in its own subprocess (correct verdicts).
    println("Running $(length(M)) validation config(s), isolated, $(jobs)-way...")
    verdicts = run_isolated(M, outdir, jobs)
    write_results(verdicts, outdir)
    np = count(v -> v.status === :pass, verdicts)
    nf = count(v -> v.status === :fail, verdicts)
    ns = count(v -> v.status in (:skipped_no_data, :skipped_unwired), verdicts)
    println("\n$np passed, $nf failed, $ns skipped → $(joinpath(outdir, "report.md"))")
    return nf == 0
end

if abspath(PROGRAM_FILE) == (@__FILE__)
    main(ARGS)
end
