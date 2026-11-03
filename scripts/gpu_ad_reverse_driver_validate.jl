# Whole-driver CUDA reverse-mode validation for CLM's compositional reverse path.
#
# Environment switches:
#   CLM_REVERSE_CANOPY=0|1   include the canopy energy block (default 1)
#   CLM_REVERSE_SURFACE=0|1  include surface-hydrology phases (default 1)
#   CLM_REVERSE_NCANOPY=N    fixed differentiated Newton iterations (default 2)
#   CLM_REVERSE_CPU=0|1       also run the expensive CPU reverse parity leg (default 0)
#   CLM_REVERSE_FD=0|1        central finite-difference check of dL/d(t_grnd) on the
#                             CUDA forward chain (default 1; forward-only, cheap)
#   CLM_REVERSE_PHASES=N       reverse only the first N assembled phases (diagnostic)
#   CLM_REVERSE_PHASE_START=N  start at phase N for isolated compilation (default 1)
#   CLM_REVERSE_DRIVER_PHASE=S isolate a named top-level driver phase (diagnostic)
#   CLM_REVERSE_EXTRA=0|1     wider checks (default 0; the CUDA-vs-CPU parts need
#                             CLM_REVERSE_CPU=1): (a) every input-gradient leaf of
#                             db.inst, CUDA vs CPU, for the default loss; (b) a second
#                             loss with random weights on many final-state fields
#                             (temperatures, soil water/ice, h2osfc, water table, fluxes),
#                             again every db.inst leaf CUDA vs CPU; (c) central FD of that
#                             loss on the CUDA forward w.r.t. other inputs (state, soil
#                             parameters, a forcing). Reuses the compiled reverse thunks.

using CLM, Enzyme, Printf, Random
import CUDA

# Preserve Float64 for CUDA while unpacking BitArray masks before CuArray adapt.
struct _CUDAHost end
CLM.Adapt.adapt_storage(::_CUDAHost, x::BitArray) = collect(Bool, x)
CLM.Adapt.adapt_storage(::_CUDAHost, x::AbstractArray) = x

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const DEFAULT_FSURDAT = joinpath(REPO_ROOT, "test_inputs", "lake", "surfdata_mixed.nc")
const DEFAULT_PARAMFILE = joinpath(homedir(), "projects", "SYMFLUENCE", "src",
                                   "symfluence", "resources", "base_settings",
                                   "CLM", "clm5_params.nc")

function _input_path(envname, default)
    path = abspath(get(ENV, envname, default))
    isfile(path) || error("$envname input not found: $path")
    path
end

function setup_forcing!(a, ng)
    T = 285.0
    for g in 1:ng
        a.forc_t_not_downscaled_grc[g] = T
        a.forc_pbot_not_downscaled_grc[g] = 85000.0
        a.forc_th_not_downscaled_grc[g] = T * (100000.0 / 85000.0)^(CLM.RAIR / CLM.CPAIR)
        a.forc_rho_not_downscaled_grc[g] = 85000.0 / (CLM.RAIR * T)
        a.forc_lwrad_not_downscaled_grc[g] = 300.0
        a.forc_vp_grc[g] = 800.0
        a.forc_hgt_grc[g] = 30.0; a.forc_topo_grc[g] = 0.0; a.forc_wind_grc[g] = 3.0
        a.forc_hgt_u_grc[g] = 30.0; a.forc_hgt_t_grc[g] = 30.0; a.forc_hgt_q_grc[g] = 30.0
        for b in 1:CLM.NUMRAD
            a.forc_solad_not_downscaled_grc[g, b] = 200.0
            a.forc_solai_grc[g, b] = 80.0
        end
        a.forc_solar_not_downscaled_grc[g] = 560.0
        a.forc_rain_not_downscaled_grc[g] = 0.0001
        a.forc_snow_not_downscaled_grc[g] = 0.0
    end
end

function warmed_state()
    fsurdat = _input_path("CLM_FSURDAT", DEFAULT_FSURDAT)
    paramfile = _input_path("CLM_PARAMFILE", DEFAULT_PARAMFILE)
    inst, bounds, filt, _ = CLM.clm_initialize!(; fsurdat, paramfile)
    setup_forcing!(inst.atm2lnd, bounds.endg)
    CLM.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
    config = CLM.CLMDriverConfig()
    fia = CLM.clump_filter_inactive_and_active
    declin, _ = CLM.compute_orbital(120.0)
    nextsw = 120.0 + 1800.0 / CLM.SECSPDAY
    for n in 1:3
        CLM.clm_drv!(config, inst, filt, fia, bounds, true, nextsw, declin, declin,
                     CLM.ORB_OBLIQR_DEFAULT, false, false, "", false;
                     nstep=n, is_first_step=(n == 1), is_beg_curr_day=(n == 1),
                     dtime=1800.0, mon=1, day=1, photosyns=inst.photosyns)
    end
    inst, bounds, filt, config
end

# Forward-only loss of the exact chain clm_drv_reverse! differentiates, with the
# step-entry t_grnd[c0] shifted by δ. Runs on whatever backend `inst` lives on.
# Returns the soil-layer t_soisno profile; L = sum(abs2, ·) over finite entries.
function forward_profile(inst, bounds, filt, config, c0, δ;
                         include_canopy, include_surface, ncanopy)
    i = deepcopy(inst)
    tg = Array(i.temperature.t_grnd_col); tg[c0] += δ
    copyto!(i.temperature.t_grnd_col, tg)
    b = CLM.driver_rev_bundle(i)
    caux = include_canopy ? CLM.canopy_rev_aux(i, bounds, filt) : nothing
    phases = CLM.driver_rev_phases(bounds, filt, config; canopy_aux = caux, inst = i,
        n_canopy = max(ncanopy, 1), dtime = 1800.0, include_surface = include_surface)
    for entry in phases
        f, cargs = entry[1], entry[2]
        viewfn = length(entry) == 3 ? entry[3] : identity
        f(viewfn(b), cargs...)
    end
    j0 = CLM.varpar.nlevsno + 1
    return Array(b.inst.temperature.t_soisno_col)[:, j0:end]
end

# Central difference of L = sum(abs2, soil t_soisno). Columns the perturbation cannot
# reach cancel exactly; non-finite (inactive/spval) entries are excluded.
function fd_gradient(inst, bounds, filt, config, c0, h; kw...)
    tp = forward_profile(inst, bounds, filt, config, c0, h; kw...)
    tm = forward_profile(inst, bounds, filt, config, c0, -h; kw...)
    d = abs2.(tp) .- abs2.(tm)
    return sum(filter(isfinite, vec(d))) / (2h)
end

# ---------------------------------------------------------------- CLM_REVERSE_EXTRA
_getpath(x, p::AbstractString) = foldl((o, f) -> o === nothing || !hasproperty(o, Symbol(f)) ?
                                        nothing : getproperty(o, Symbol(f)), split(p, "."); init = x)

# Output fields of the wide loss (paths under inst). Missing paths are skipped.
# (zwt_col / qflx_drain_col are back in: the host and device water-table phases now
# both follow production clm_drv!; test/test_watertable_reverse.jl.)
const WIDE_FIELDS = ["temperature.t_soisno_col", "temperature.t_grnd_col",
    "temperature.t_veg_patch", "water.waterstatebulk_inst.ws.h2osoi_liq_col",
    "water.waterstatebulk_inst.ws.h2osoi_ice_col", "water.waterstatebulk_inst.ws.h2osfc_col",
    "soilhydrology.zwt_col", "energyflux.eflx_sh_tot_patch", "energyflux.eflx_lh_tot_patch",
    "energyflux.eflx_soil_grnd_patch", "water.waterfluxbulk_inst.wf.qflx_surf_col",
    "water.waterfluxbulk_inst.wf.qflx_drain_col", "water.waterfluxbulk_inst.wf.qflx_infl_col"]

# Random weights per field, zero on non-finite/spval entries of the final state, scaled
# by 1/max|field| so every field contributes at a comparable size.
function wide_weights(final_inst)
    W = Dict{String, Array{Float64}}()
    for p in WIDE_FIELDS
        a = _getpath(final_inst, p)
        a isa AbstractArray{<:AbstractFloat} || continue
        x = Array(a); ok = isfinite.(x) .& (abs.(x) .< 1e30)
        any(ok) || continue
        sc = max(maximum(abs, x[ok]), 1e-6)
        rng = Random.MersenneTwister(hash(p))
        W[p] = ok .* (0.5 .+ rand(rng, size(x)...)) ./ sc
    end
    W
end

wide_loss_terms(final_inst, W) = Dict(p => Array(_getpath(final_inst, p)) for p in keys(W))

# Seed for clm_drv_reverse!: fills W from the final state on first use (so CPU and
# FD can reuse the CUDA run's weights), then seeds db.inst with it.
function wide_seed!(W)
    (db, b) -> begin
        isempty(W) && merge!(W, wide_weights(b.inst))
        for (p, w) in W
            copyto!(_getpath(db.inst, p), w)
        end
        nothing
    end
end

function _gleaves!(out, x, path, seen)
    if x isa AbstractArray{<:AbstractFloat}
        out[path] = x
    elseif x isa AbstractArray || x isa Number || x isa Symbol || x isa String ||
           x isa Function || x === nothing
    else
        objectid(x) in seen && return out
        ismutable(x) && push!(seen, objectid(x))
        for f in fieldnames(typeof(x))
            isdefined(x, f) || continue
            _gleaves!(out, getfield(x, f), string(path, ".", f), seen)
        end
    end
    out
end

# Compare every Float leaf of two gradient trees (CPU vs CUDA), matched by path. Per
# leaf: max|Δ| / max(max|cpu|, max|gpu|) over entries finite on both. Returns the worst
# value (Inf if the finite patterns differ). A leaf whose largest difference is at or
# below `atol` (roundoff-level gradients, e.g. 1e-20 where one backend gives exactly 0)
# counts as agreeing.
function compare_grad_trees(label, dcpu, dgpu; show = 10, atol = 1e-15)
    lc = _gleaves!(Dict{String,Any}(), dcpu, "inst", Set{UInt}())
    lg = _gleaves!(Dict{String,Any}(), dgpu, "inst", Set{UInt}())
    rows = Tuple{Float64,Float64,Float64,Int,String}[]; nfin = 0; nnz = 0; nleaf = 0
    for (p, a) in lc
        haskey(lg, p) || continue
        ah = Array(a); gh = Array(lg[p]); size(ah) == size(gh) || continue
        nleaf += 1
        fa = isfinite.(ah); fg = isfinite.(gh); m = fa .& fg
        nfin += count(fa .!= fg)
        sc = max(maximum(abs, ah[m]; init = 0.0), maximum(abs, gh[m]; init = 0.0))
        nz = count(!iszero, ah[m]); nnz += nz
        sc == 0 && continue
        ad = maximum(abs.(ah[m] .- gh[m]); init = 0.0)
        push!(rows, (ad <= atol ? 0.0 : ad / sc, ad, sc, nz, p))
    end
    sort!(rows; rev = true)
    worst = isempty(rows) ? 0.0 : rows[1][1]
    @printf("  [%s] %d leaves compared, %d nonzero gradient entries, %d finite-pattern mismatches, worst leaf rel=%.3e\n",
            label, nleaf, nnz, nfin, worst)
    for r in first(rows, show)
        @printf("      rel=%.3e abs=%.3e max|g|=%.3e nnz=%-6d %s\n", r...)
    end
    nfin > 0 && (worst = Inf)
    worst
end

# Final inst after the forward chain from `inst_entry` with `perturb!(inst)` applied
# first. The canopy aux is built from the UNPERTURBED entry state: in the reverse the
# aux is a Const argument, so a forcing read through it gets no adjoint; FD with the
# aux frozen differentiates the same function the reverse does.
function forward_final(inst_entry, bounds, filt, config, perturb!; include_canopy,
                       include_surface, ncanopy)
    i0 = deepcopy(inst_entry); i = deepcopy(inst_entry); perturb!(i)
    b = CLM.driver_rev_bundle(i)
    caux = include_canopy ? CLM.canopy_rev_aux(i0, bounds, filt) : nothing
    phases = CLM.driver_rev_phases(bounds, filt, config; canopy_aux = caux, inst = i,
        n_canopy = max(ncanopy, 1), dtime = 1800.0, include_surface = include_surface)
    for entry in phases
        f, cargs = entry[1], entry[2]
        viewfn = length(entry) == 3 ? entry[3] : identity
        f(viewfn(b), cargs...)
    end
    b.inst
end

_perturb_entry(path, idx, δ) = inst -> begin
    a = _getpath(inst, path); h = Array(a); h[idx...] += δ; copyto!(a, h)
end

function fd_wide(inst_entry, bounds, filt, config, W, path, idx, h; kw...)
    fp = wide_loss_terms(forward_final(inst_entry, bounds, filt, config,
                                       _perturb_entry(path, idx, h); kw...), W)
    fm = wide_loss_terms(forward_final(inst_entry, bounds, filt, config,
                                       _perturb_entry(path, idx, -h); kw...), W)
    s = 0.0
    for (p, w) in W
        d = fp[p] .- fm[p]; m = (w .!= 0) .& isfinite.(d)
        s += sum(w[m] .* d[m])
    end
    s / (2h)
end

function main()
    CUDA.functional() || error("CUDA is not functional")
    include_canopy = get(ENV, "CLM_REVERSE_CANOPY", "1") == "1"
    include_surface = get(ENV, "CLM_REVERSE_SURFACE", "1") == "1"
    ncanopy = parse(Int, get(ENV, "CLM_REVERSE_NCANOPY", "2"))
    run_cpu = get(ENV, "CLM_REVERSE_CPU", "0") == "1"
    run_fd = get(ENV, "CLM_REVERSE_FD", "1") == "1"
    run_extra = get(ENV, "CLM_REVERSE_EXTRA", "0") == "1"

    inst64, bounds, filt, config = warmed_state()
    # CUDA supports CLM's native Float64 state. Only unpack packed Bool masks
    # before moving the complete state tree to the device.
    inst_cpu = deepcopy(inst64)
    filt_cpu = deepcopy(filt)
    inst_gpu = CLM.Adapt.adapt(CUDA.CuArray,
                               CLM.Adapt.adapt(_CUDAHost(), deepcopy(inst64)))
    filt_gpu = CLM.Adapt.adapt(CUDA.CuArray,
                               CLM.Adapt.adapt(_CUDAHost(), deepcopy(filt)))
    # clm_drv_reverse! leaves its inst at the step's final state; keep an entry copy
    # for the finite-difference leg.
    inst_gpu_entry = (run_fd || run_extra) ? deepcopy(inst_gpu) : nothing
    inst_cpu_entry = run_extra ? deepcopy(inst_cpu) : nothing

    kwargs = (; include_canopy, include_surface,
               n_canopy = include_canopy ? ncanopy : nothing)
    println("starting CUDA compositional reverse ..."); flush(stdout)
    maxphases = parse(Int, get(ENV, "CLM_REVERSE_PHASES", "0"))
    driver_phase = get(ENV, "CLM_REVERSE_DRIVER_PHASE", "")
    if !isempty(driver_phase)
        b_gpu = CLM.driver_rev_bundle(inst_gpu)
        caux = include_canopy ? CLM.canopy_rev_aux(inst_gpu, bounds, filt_gpu) : nothing
        phases = CLM.driver_rev_phases(bounds, filt_gpu, config; inst=inst_gpu,
            canopy_aux=caux, n_canopy=max(ncanopy, 1), dtime=1800.0,
            include_surface=include_surface)
        phases = Any[p for p in phases if occursin(driver_phase, string(first(p)))]
        isempty(phases) && error("CLM_REVERSE_DRIVER_PHASE=$(repr(driver_phase)) matched no phase")
        length(phases) == 1 || error("driver phase selector matched $(length(phases)) phases")
        println("diagnostic driver phase: ", first(first(phases)))
        seed_driver_phase!(db,b) = begin
            j0=CLM.varpar.nlevsno+1
            @views db.inst.temperature.t_soisno_col[:,j0:end] .=
                2 .* b.inst.temperature.t_soisno_col[:,j0:end]
            db.inst.soilhydrology.zwt_col .= one(eltype(db.inst.soilhydrology.zwt_col))
            db.inst.soilhydrology.zwt_perched_col .=
                one(eltype(db.inst.soilhydrology.zwt_perched_col))
        end
        db_gpu = CLM.compositional_reverse!(phases, b_gpu, seed_driver_phase!)
    elseif maxphases > 0
        b_gpu = CLM.soiltemp_rev_bundle(inst_gpu)
        caux = include_canopy ? CLM.canopy_rev_aux(inst_gpu, bounds, filt_gpu) : nothing
        phases = CLM.soiltemp_driver_rev_phases(bounds, filt_gpu)
        phase_start = parse(Int, get(ENV, "CLM_REVERSE_PHASE_START", "1"))
        phase_stop = min(phase_start + maxphases - 1, length(phases))
        phases = phases[phase_start:phase_stop]
        println("diagnostic phase prefix: ", join(string.(first.(phases)), ", "))
        seed_soiltemp!(db,b) = begin
            j0=CLM.varpar.nlevsno+1
            @views db.inst.temperature.t_soisno_col[:,j0:end] .=
                2 .* b.inst.temperature.t_soisno_col[:,j0:end]
            # Activate intermediate outputs as well when a phase is run in
            # isolation; otherwise an all-zero local cotangent is only a compile gate.
            db.soiltemp.tk .= one(eltype(db.soiltemp.tk))
            db.soiltemp.tk_h2osfc .= one(eltype(db.soiltemp.tk_h2osfc))
            db.soiltemp.cv .= one(eltype(db.soiltemp.cv))
            db.soiltemp.lwrad_emit .= one(eltype(db.soiltemp.lwrad_emit))
            db.soiltemp.dlwrad_emit .= one(eltype(db.soiltemp.dlwrad_emit))
            db.soiltemp.hs_patch .= one(eltype(db.soiltemp.hs_patch))
            db.soiltemp.dhsdT_patch .= one(eltype(db.soiltemp.dhsdT_patch))
            db.soiltemp.fn .= one(eltype(db.soiltemp.fn))
            db.soiltemp.tvector .= one(eltype(db.soiltemp.tvector))
        end
        db_gpu = CLM.compositional_reverse!(phases, b_gpu, seed_soiltemp!)
    else
        db_gpu = CLM.clm_drv_reverse!(inst_gpu, bounds, filt_gpu, config; kwargs...)
    end
    CUDA.synchronize()
    println("CUDA compositional reverse complete"); flush(stdout)

    hc = Array(filt_gpu.hydrologyc)
    c0 = first(c for c in bounds.begc:bounds.endc if hc[c])
    ggpu = Array(db_gpu.inst.temperature.t_grnd_col)[c0]
    gt = Array(db_gpu.inst.temperature.t_soisno_col)
    gmax = maximum(abs, filter(isfinite, vec(gt)); init=0.0)
    @printf("whole-driver reverse: canopy=%s surface=%s ncanopy=%d\n",
            include_canopy, include_surface, ncanopy)
    pass = isfinite(ggpu) && isfinite(gmax)
    if run_cpu
        println("starting CPU compositional reverse ..."); flush(stdout)
        db_cpu = CLM.clm_drv_reverse!(inst_cpu, bounds, filt_cpu, config; kwargs...)
        gcpu = db_cpu.inst.temperature.t_grnd_col[c0]
        err = abs(Float64(ggpu) - Float64(gcpu))
        tol = 2e-9 * max(abs(Float64(gcpu)), 1.0)
        @printf("  dL/d(t_grnd[%d]): CPU=% .8e CUDA=% .8e absdiff=%.3e tol=%.3e\n",
                c0, gcpu, ggpu, err, tol)
        pass &= isfinite(gcpu) && err <= tol
        if run_extra && maxphases == 0 && isempty(driver_phase)
            w1 = compare_grad_trees("default loss, every db.inst leaf", db_cpu.inst, db_gpu.inst)
            @printf("  default-loss all-leaf CUDA vs CPU: worst rel %.3e (tol 1e-8)\n", w1)
            pass &= w1 <= 1e-8
        end
    else
        @printf("  dL/d(t_grnd[%d]): CUDA=% .8e; max|dL/d(t_soisno)|=%.3e (CPU parity skipped)\n",
                c0, ggpu, gmax)
    end
    if run_fd && maxphases == 0 && isempty(driver_phase)
        println("starting CUDA finite-difference check ..."); flush(stdout)
        fkw = (; include_canopy, include_surface, ncanopy)
        hs = (1e-2, 5e-3, 2.5e-3)
        fds = [fd_gradient(inst_gpu_entry, bounds, filt_gpu, config, c0, h; fkw...) for h in hs]
        for (h, f) in zip(hs, fds)
            @printf("  FD(h=%.2e) = % .8e\n", h, f)
        end
        g_rich = (4 * fds[end] - fds[end-1]) / 3
        lo, hi = extrema(fds)
        bracketed = (lo - abs(lo) * 1e-9) <= ggpu <= (hi + abs(hi) * 1e-9)
        rel = abs(ggpu - g_rich) / max(abs(g_rich), 1e-12)
        @printf("  dL/d(t_grnd[%d]): CUDA reverse=% .8e Richardson FD=% .8e rel=%.3e bracketed=%s\n",
                c0, ggpu, g_rich, rel, string(bracketed))
        pass &= rel < 5e-5 || bracketed
    end
    if run_extra && maxphases == 0 && isempty(driver_phase)
        println("starting wide-loss reverse (random weights on final-state fields) ..."); flush(stdout)
        W = Dict{String, Array{Float64}}()
        t0 = time()
        dbw_gpu = CLM.clm_drv_reverse!(deepcopy(inst_gpu_entry), bounds, filt_gpu, config;
                                       seed_bang! = wide_seed!(W), kwargs...)
        CUDA.synchronize()
        @printf("  CUDA wide reverse %.1f s; seeded fields: %s\n", time() - t0,
                join(sort(collect(keys(W))), ", "))
        if run_cpu
            t0 = time()
            dbw_cpu = CLM.clm_drv_reverse!(deepcopy(inst_cpu_entry), bounds, filt_cpu, config;
                                           seed_bang! = wide_seed!(copy(W)), kwargs...)
            @printf("  CPU wide reverse %.1f s\n", time() - t0)
            w2 = compare_grad_trees("wide loss, every db.inst leaf", dbw_cpu.inst, dbw_gpu.inst)
            @printf("  wide-loss all-leaf CUDA vs CPU: worst rel %.3e (tol 1e-8)\n", w2)
            pass &= w2 <= 1e-8
        end
        # FD of the wide loss on the CUDA forward chain w.r.t. other inputs.
        j0 = CLM.varpar.nlevsno + 1
        ev = Array(filt_gpu.exposedvegp); pcol = Array(inst_gpu_entry.patch.column)
        p0 = something(findfirst(p -> ev[p] && pcol[p] == c0, eachindex(ev)), findfirst(ev))
        fdkw = (; include_canopy, include_surface, ncanopy)
        checks = [
            ("temperature.t_grnd_col", (c0,), 1e-2),
            # a soil layer away from TFRZ: at exactly 273.15 the central difference
            # straddles the phase-change switch
            ("temperature.t_soisno_col", (c0, something(findfirst(j -> abs(
                Array(inst_gpu_entry.temperature.t_soisno_col)[c0, j] - CLM.TFRZ) > 0.5,
                j0:size(inst_gpu_entry.temperature.t_soisno_col, 2)), 3) + j0 - 1), 1e-2),
            ("temperature.t_veg_patch", (p0,), 1e-2),
            ("water.waterstatebulk_inst.ws.h2osoi_liq_col", (c0, j0 + 1), 0.0),
            ("soilstate.watsat_col", (c0, 2), 0.0),
            ("soilstate.bsw_col", (c0, 3), 0.0),
            ("atm2lnd.forc_lwrad_downscaled_col", (c0,), 1e-1),
        ]
        for (path, idx, h0) in checks
            a = _getpath(inst_gpu_entry, path); a === nothing && continue
            x = Array(a)[idx...]
            g = Array(_getpath(dbw_gpu.inst, path))[idx...]
            h1 = h0 > 0 ? h0 : 1e-4 * max(abs(x), 1e-3)
            hs = (h1, h1 / 2, h1 / 4)
            fds = [fd_wide(inst_gpu_entry, bounds, filt_gpu, config, W, path, idx, h; fdkw...)
                   for h in hs]
            g_rich = (4 * fds[3] - fds[2]) / 3
            lo, hi = extrema(fds)
            bracketed = (lo - abs(lo) * 1e-9) <= g <= (hi + abs(hi) * 1e-9)
            rel = abs(g - g_rich) / max(abs(g_rich), 1e-12)
            ok = (rel < 5e-5 || bracketed) && isfinite(g)
            @printf("  dLw/d %s%s (x=%.6g): CUDA reverse=% .8e FD=[% .8e % .8e % .8e] Richardson=% .8e rel=%.3e bracketed=%s %s\n",
                    path, string(idx), x, g, fds..., g_rich, rel, string(bracketed),
                    ok ? "ok" : "MISMATCH")
            pass &= ok
        end
    end
    println(pass ? "  PASS ✓" : "  FAIL ✗")
    return pass ? 0 : 1
end

exit(main())
