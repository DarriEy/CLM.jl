# Whole-driver CUDA reverse-mode validation for CLM's compositional reverse path.
#
# Environment switches:
#   CLM_REVERSE_CANOPY=0|1   include the canopy energy block (default 1)
#   CLM_REVERSE_SURFACE=0|1  include surface-hydrology phases (default 1)
#   CLM_REVERSE_NCANOPY=N    fixed differentiated Newton iterations (default 2)
#   CLM_REVERSE_CPU=0|1       also run the expensive CPU reverse parity leg (default 0)
#   CLM_REVERSE_PHASES=N       reverse only the first N assembled phases (diagnostic)
#   CLM_REVERSE_PHASE_START=N  start at phase N for isolated compilation (default 1)

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

function main()
    CUDA.functional() || error("CUDA is not functional")
    include_canopy = get(ENV, "CLM_REVERSE_CANOPY", "1") == "1"
    include_surface = get(ENV, "CLM_REVERSE_SURFACE", "1") == "1"
    ncanopy = parse(Int, get(ENV, "CLM_REVERSE_NCANOPY", "2"))
    run_cpu = get(ENV, "CLM_REVERSE_CPU", "0") == "1"

    inst64, bounds, filt, config = warmed_state()
    # CUDA supports CLM's native Float64 state. Only unpack packed Bool masks
    # before moving the complete state tree to the device.
    inst_cpu = deepcopy(inst64)
    filt_cpu = deepcopy(filt)
    inst_gpu = CLM.Adapt.adapt(CUDA.CuArray,
                               CLM.Adapt.adapt(_CUDAHost(), deepcopy(inst64)))
    filt_gpu = CLM.Adapt.adapt(CUDA.CuArray,
                               CLM.Adapt.adapt(_CUDAHost(), deepcopy(filt)))

    kwargs = (; include_canopy, include_surface,
               n_canopy = include_canopy ? ncanopy : nothing)
    println("starting CUDA compositional reverse ..."); flush(stdout)
    maxphases = parse(Int, get(ENV, "CLM_REVERSE_PHASES", "0"))
    if maxphases > 0
        b_gpu = CLM.soiltemp_rev_bundle(inst_gpu)
        caux = include_canopy ? CLM.canopy_rev_aux(inst_gpu, bounds, filt_gpu) : nothing
        staux = CLM.soiltemp_rev_aux(bounds, filt_gpu)
        phases = Any[(CLM.soiltemp_tk_layers_array_phase!, (staux,), CLM.soiltemp_tk_layers_view),
                     (CLM.soiltemp_tk_interface_array_phase!, (staux,), CLM.soiltemp_tk_interface_view),
                     (CLM.soiltemp_cv_array_phase!, (staux,), CLM.soiltemp_cv_view),
                     (CLM.soiltemp_lwrad_array_phase!, (staux,), CLM.soiltemp_lwrad_view),
                     (CLM.soiltemp_gnet_patch_rural_array_phase!, (staux,), CLM.soiltemp_gnet_patch_rural_view),
                     (CLM.soiltemp_gnet_patch_rural_contrib_array_phase!, (staux,), CLM.soiltemp_gnet_patch_rural_contrib_view),
                     (CLM.soiltemp_gnet_patch_urban_array_phase!, (staux,), CLM.soiltemp_gnet_patch_view),
                     (CLM.soiltemp_gnet_patch_reduce_array_phase!, (staux,), CLM.soiltemp_gnet_patch_view),
                     (CLM.soiltemp_gnet_snicar_rev_phase!, (staux,)),
                     (CLM.soiltemp_heatdiff_rev_phase!, (staux,)),
                     (CLM.soiltemp_solve_rev_phase!, (staux,)),
                     (CLM.soiltemp_post_rev_phase!, (staux,))]
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
    println(pass ? "  PASS ✓" : "  FAIL ✗")
    return pass ? 0 : 1
end

exit(main())
