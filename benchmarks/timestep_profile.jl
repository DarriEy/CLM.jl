# benchmarks/timestep_profile.jl
#
# Performance baseline for the CLM.jl timestep driver. Initializes the model
# once, then runs `clm_drv!` for K timesteps and reports wall-clock time and
# allocations per timestep. This is the regression baseline for the GPU port
# (Phase 0 of the AD+GPU roadmap): compare CPU timings here against batched /
# GPU runs as kernelization lands.
#
# Usage:
#   julia --project=. benchmarks/timestep_profile.jl [nsteps]
#
# Data paths default to the Bow-at-Banff lumped domain and can be overridden:
#   CLM_FSURDAT=/path/to/surfdata_clm.nc \
#   CLM_PARAMFILE=/path/to/clm5_params.nc \
#   julia --project=. benchmarks/timestep_profile.jl 200
#
# NOTE: per-phase timing breakdown (radiation / hydrology / temperature /
# fluxes / biogeochem) is intentionally not instrumented here to avoid editing
# the driver. Use Julia's profiler for that:
#   using Profile; @profile run_steps(...); Profile.print()

using CLM
using Printf

const FSURDAT = get(ENV, "CLM_FSURDAT",
    "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/surfdata_clm.nc")
const PARAMFILE = get(ENV, "CLM_PARAMFILE",
    "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/clm5_params.nc")

# Steady, finite forcing so the driver exercises a realistic spring/melt regime.
function setup_forcing!(a2l, T_forc, ng)
    for g in 1:ng
        a2l.forc_t_not_downscaled_grc[g] = T_forc
        a2l.forc_pbot_not_downscaled_grc[g] = 85000.0
        a2l.forc_th_not_downscaled_grc[g] = T_forc * (100000.0 / 85000.0)^(CLM.RAIR / CLM.CPAIR)
        a2l.forc_rho_not_downscaled_grc[g] = 85000.0 / (CLM.RAIR * T_forc)
        a2l.forc_lwrad_not_downscaled_grc[g] = 250.0
        a2l.forc_vp_grc[g] = 300.0
        a2l.forc_hgt_grc[g] = 30.0
        a2l.forc_topo_grc[g] = 0.0
        a2l.forc_wind_grc[g] = 3.0
        for b in 1:CLM.NUMRAD
            a2l.forc_solad_not_downscaled_grc[g, b] = 100.0
            a2l.forc_solai_grc[g, b] = 50.0
        end
        a2l.forc_solar_not_downscaled_grc[g] = 300.0
        a2l.forc_rain_not_downscaled_grc[g] = 0.0
        a2l.forc_snow_not_downscaled_grc[g] = 0.0001
    end
end

function run_steps(nsteps::Int)
    if !isfile(FSURDAT) || !isfile(PARAMFILE)
        error("""
        Forcing/parameter files not found:
          fsurdat   = $FSURDAT
          paramfile = $PARAMFILE
        Set CLM_FSURDAT and CLM_PARAMFILE to valid NetCDF paths.""")
    end

    config = CLM.CLMDriverConfig()
    T_forc = 285.0  # spring; active fluxes + some snow

    @printf("Initializing model (fsurdat=%s)...\n", basename(FSURDAT))
    init_stats = @timed CLM.clm_initialize!(; fsurdat=FSURDAT, paramfile=PARAMFILE)
    (inst, bounds, filt, tm) = init_stats.value
    ng = bounds.endg
    @printf("  init: %.3f s, %.1f MiB allocated\n",
        init_stats.time, init_stats.bytes / 2^20)
    @printf("  grid: ng=%d nl=%d nc=%d np=%d\n",
        bounds.endg, bounds.endl, bounds.endc, bounds.endp)

    filt_ia = CLM.clump_filter_inactive_and_active
    calday = 1.0
    (declin, eccf) = CLM.compute_orbital(calday)
    nextsw_cday = calday + 1800.0 / CLM.SECSPDAY

    setup_forcing!(inst.atm2lnd, T_forc, ng)
    CLM.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)

    step!(n) = CLM.clm_drv!(config, inst, filt, filt_ia, bounds,
        true, nextsw_cday, declin, declin, CLM.ORB_OBLIQR_DEFAULT,
        false, false, "", false;
        nstep=n, is_first_step=(n == 1), is_beg_curr_day=(n == 1),
        dtime=1800.0, mon=1, day=1, photosyns=inst.photosyns)

    # Warm up (compilation + first-step transients) before timing.
    @printf("Warming up (3 steps)...\n")
    for n in 1:3
        step!(n)
    end

    @printf("Timing %d steps...\n", nsteps)
    times = Float64[]
    bytes = Float64[]
    for n in 4:(3 + nsteps)
        s = @timed step!(n)
        push!(times, s.time)
        push!(bytes, s.bytes / 2^20)
    end

    sort!(times)
    mean_t = sum(times) / length(times)
    median_t = times[cld(length(times), 2)]
    @printf("\n=== Timestep profile (%d steps, ng=%d nc=%d np=%d) ===\n",
        nsteps, ng, bounds.endc, bounds.endp)
    @printf("  mean   : %8.3f ms/step\n", 1e3 * mean_t)
    @printf("  median : %8.3f ms/step\n", 1e3 * median_t)
    @printf("  min    : %8.3f ms/step\n", 1e3 * times[1])
    @printf("  max    : %8.3f ms/step\n", 1e3 * times[end])
    @printf("  alloc  : %8.1f MiB/step (mean)\n", sum(bytes) / length(bytes))
    return nothing
end

nsteps = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 100
run_steps(nsteps)
