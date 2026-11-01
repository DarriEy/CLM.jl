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

using CLM, Enzyme, Printf
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

function main()
    CUDA.functional() || error("CUDA is not functional")
    include_canopy = get(ENV, "CLM_REVERSE_CANOPY", "1") == "1"
    include_surface = get(ENV, "CLM_REVERSE_SURFACE", "1") == "1"
    ncanopy = parse(Int, get(ENV, "CLM_REVERSE_NCANOPY", "2"))
    run_cpu = get(ENV, "CLM_REVERSE_CPU", "0") == "1"
    run_fd = get(ENV, "CLM_REVERSE_FD", "1") == "1"

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
    inst_gpu_entry = run_fd ? deepcopy(inst_gpu) : nothing

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
    println(pass ? "  PASS ✓" : "  FAIL ✗")
    return pass ? 0 : 1
end

exit(main())
