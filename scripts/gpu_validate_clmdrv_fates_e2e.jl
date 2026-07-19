# ==========================================================================
# gpu_validate_clmdrv_fates_e2e.jl — whole clm_drv! with use_fates=true on Metal.
#
# FATES is CPU-only by design (the pointer-linked site/patch/cohort hierarchy is
# Union-attached to CLMInstances and PASSES THROUGH Adapt unchanged — a Metal inst
# keeps the host FATES object while all CLM state arrays become MtlArrays). So the
# goal here is a HOST-FALLBACK: the non-FATES biogeophysics runs on-device, and each
# FATES driver hook packs device→host, runs FATES on the CPU, and unpacks host→device.
#
# This harness adapts a real cold-started tropical FATES inst (build() from
# fates_longhorizon.jl) to Metal and runs clm_drv! on-device to (phase 1) find where
# the FATES hooks still scalar-index device arrays, and (phase 2) confirm CPU-vs-GPU
# parity on the FATES-driven fields.
#
#   julia +1.12 --project=scripts scripts/gpu_validate_clmdrv_fates_e2e.jl
# ==========================================================================
using CLM, Printf, Dates
include(joinpath(@__DIR__, "gpu_backends.jl"))
const _C = CLM
include(joinpath(@__DIR__, "gpu_adapt.jl"))
include(joinpath(@__DIR__, "fates_longhorizon.jl"))   # reuse build() (guarded exit)

df(x) = device_array_type()(x)

# Static daytime tropical forcing on the HOST inst (grc + downscaled col), so after
# adapting there is a consistent on-device forcing without a per-step host read.
function set_static_forcing!(inst, bounds)
    set_forcing!(inst, 1.0)             # tropical grc forcing (fates_longhorizon helper)
    # Steady tropical rain so the root-zone soil stays moist (btran>0) across the
    # perpetual-noon run — without it the soil dries in a day, btran→0, GPP→0 and the
    # stand can't grow. ~8.6 mm/day, replenishing transpiration.
    inst.atm2lnd.forc_rain_not_downscaled_grc[1] = 1.0e-4
    _C.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
    return nothing
end

function run_steps!(inst, config, filt, filt_ia, bounds; nsteps=2)
    dtime = 1800.0
    # Aripuana lat -7.5°, lon -59.6°. Put the sun overhead: near-equinox declination
    # ≈ lat (coszen≈cos(lat-decl)≈1 at local noon) and a calday whose fractional part
    # gives hour-angle≈0 for this longitude (frac=(π-lon)/2π ≈ 0.665 → ~local solar noon).
    declin = -0.131
    base_calday = 80.0 + 0.665
    steps_per_day = Int(round(86400 / dtime))
    for i in 1:nsteps
        day_idx = (i - 1) ÷ steps_per_day
        is_beg  = ((i - 1) % steps_per_day == 0)          # day boundary → FATES daily dynamics
        calday  = base_calday + day_idx                    # advance the day, keep local-noon sun
        nextsw_cday = calday + dtime/86400.0
        # No Metal.@allowscalar wrap needed here: clm_drv! wraps its step in
        # GPUArraysCore.@allowscalar internally under use_fates (the CPU-only FATES
        # host-fallback), so use_fates on a GPU inst "just works".
        _C.clm_drv!(config, inst, filt, filt_ia, bounds, true, nextsw_cday, declin, declin, 0.4091,
            false, false, "20260101", false; nstep=i, is_first_step=(i==1), is_beg_curr_day=is_beg,
            is_end_curr_day=false, is_beg_curr_year=false, dtime=dtime, mon=3, day=21+day_idx,
            secs=57456, jday=81+day_idx, photosyns=inst.photosyns)
    end
    return nothing
end

# Extract the FATES-driven scalars for parity (all host-side; elai is a CLM patch array).
function fates_summary(inst)
    site = inst.fates.sites[1]
    cen = census(inst)                              # ncoh, npatch (fates_longhorizon helper)
    carbon = total_site_carbon(site); md = max_dbh(site)
    gpp = 0.0; btr = 0.0
    cp = site.oldest_patch
    while cp !== nothing
        try; btr = max(btr, maximum(x->isfinite(x) ? x : 0.0, cp.btran_ft)); catch; end
        cc = cp.tallest
        while cc !== nothing; isfinite(cc.gpp_tstep) && (gpp += cc.gpp_tstep*cc.n); cc = cc.shorter; end
        cp = cp.younger
    end
    elai = Array(inst.canopystate.elai_patch)       # device→host
    return (ncoh=cen.ncoh, npatch=cen.npatch, carbon=carbon, maxdbh=md,
            gpp=gpp, btran=btr, elai2=Float64(elai[2]))
end

reld(a, b) = abs(a - b) / (1.0 + abs(a))

function main()
    gpu_functional() || (println("No GPU backend detected."); return 0)
    println("="^70); println("  clm_drv! + use_fates on Metal — CPU-vs-GPU host-fallback parity"); println("="^70)
    nsteps = parse(Int, get(ENV, "FATES_GPU_STEPS", "144"))   # ~3 days of perpetual-noon growth

    # Two independent, identical cold-started insts (clm_initialize! is deterministic).
    instH, _fH, config, bounds, filt, filt_ia = build(); set_static_forcing!(instH, bounds)
    instS, _fS, _cfg, boundsS, filtS, filt_iaS = build(); set_static_forcing!(instS, bounds)
    instD    = mf(device_array_type(), instS)
    filtD    = mf(device_array_type(), filtS)
    filt_iaD = mf(device_array_type(), filt_ia)
    println("  built CPU + Metal insts (nc=$(bounds.endc)); fates preserved through adapt: $(instD.fates === instS.fates)")
    println("  t_veg_patch device type: $(typeof(instD.temperature.t_veg_patch))")
    cold = fates_summary(instH)
    @printf("  cold start: ncoh=%d carbon=%.5g maxdbh=%.5f\n\n", cold.ncoh, cold.carbon, cold.maxdbh)

    println("  running $nsteps steps ($(nsteps÷48) days) on each …")
    run_steps!(instH, config, filt,  filt_ia,  bounds;  nsteps=nsteps)   # CPU (Float64)
    run_steps!(instD, config, filtD, filt_iaD, boundsS; nsteps=nsteps)   # Metal (Float32, host-fallback)

    sH = fates_summary(instH); sD = fates_summary(instD)
    @printf("\n  growth: CPU carbon %.5g→%.5g dbh %.5f→%.5f | Metal carbon %.5g→%.5g dbh %.5f→%.5f\n",
        cold.carbon, sH.carbon, cold.maxdbh, sH.maxdbh, cold.carbon, sD.carbon, cold.maxdbh, sD.maxdbh)
    grew = sH.carbon > cold.carbon*1.0001 && sH.maxdbh > cold.maxdbh
    println("  stand grew (CPU): $grew")
    println("\n  " * "-"^60)
    @printf("  %-10s %16s %16s %10s\n", "field", "CPU (f64)", "Metal (f32)", "rel")
    checks = Bool[]
    for k in (:ncoh, :npatch, :carbon, :maxdbh, :gpp, :btran, :elai2)
        vH = Float64(getfield(sH, k)); vD = Float64(getfield(sD, k))
        r = reld(vH, vD); ok = (k in (:ncoh, :npatch)) ? (vH == vD) : (r < 2e-3)
        push!(checks, ok)
        @printf("  %-10s %16.6g %16.6g %10.2e %s\n", k, vH, vD, r, ok ? "" : "  ✗")
    end
    nfail = count(!, checks)
    println("\n  " * (nfail == 0 ?
        "★ FATES runs on the GPU (host-fallback) at CPU parity ($(length(checks))/$(length(checks)) fields)" :
        "✗ $nfail field(s) diverged"))
    return nfail
end
exit(main())
