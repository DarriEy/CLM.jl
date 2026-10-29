# ==========================================================================
# Host-side SPITFIRE fire weather: the WIND24 / PREC24 / RH24 running means
# (atm2lndType.F90 / Wateratm2lndBulkType.F90 UpdateAccVars under use_fates) and
# their hand-off to FATES bc_in (clmfates_interfaceMod.F90 dynamics_driv).
# ==========================================================================
using Test
using CLM

@testset "FATES fire-weather host accumulators (WIND24/PREC24/RH24)" begin
    old_use_fates = CLM.varctl.use_fates
    try
        CLM.varctl.use_fates = true
        ng, nc, np = 1, 1, 3
        a = CLM.Atm2LndData()
        CLM.atm2lnd_init!(a, ng, nc, np)

        # InitAccBuffer init_value (copied into the patch fields by InitAccVars).
        @test all(==(0.0), a.wind24_patch)
        @test all(==(0.0), a.prec24_patch)
        @test all(==(100.0), a.rh24_patch)
        @test all(==(0.0), a.fates24_nsteps_patch)

        pg = ones(Int, np); pc = ones(Int, np)
        dtime = 3600
        w1 = CLM.accum_window_steps(1, dtime)
        @test w1 == 24

        # Reference: accumulMod::update_accum_field_runmean, per point, active only.
        ref_val(v, n, x) = ((n - 1) * v + x) / n
        rv = [0.0, 0.0, 0.0]; rp = [0.0, 0.0, 0.0]; rr = [100.0, 100.0, 100.0]
        rn = [0, 0, 0]

        active = [true, true, false]   # patch 3: a FATES patch not yet created
        for s in 1:60
            s == 31 && (active[3] = true)  # ... until a disturbance spawns it
            wind = 2.0 + sin(s / 3)
            rain = s % 7 == 0 ? 1.0e-3 : 0.0
            snow = s % 11 == 0 ? 2.0e-4 : 0.0
            rh = 60.0 + 30.0 * cos(s / 5)
            a.forc_wind_grc[1] = wind
            a.forc_rain_downscaled_col[1] = rain
            a.forc_snow_downscaled_col[1] = snow
            a.forc_rh_grc[1] = rh
            CLM.atm2lnd_update_acc_vars!(a, 1:np, pg, pc; nstep = s, dtime = dtime,
                                         patch_active = active)
            for p in 1:np
                active[p] || continue
                rn[p] = min(rn[p] + 1, w1)
                rv[p] = ref_val(rv[p], rn[p], wind)
                rp[p] = ref_val(rp[p], rn[p], rain + snow)
                rr[p] = ref_val(rr[p], rn[p], rh)
            end
            if s == 30
                # inactive patch is untouched: still at the init values
                @test a.wind24_patch[3] == 0.0
                @test a.rh24_patch[3] == 100.0
                @test a.fates24_nsteps_patch[3] == 0.0
            end
        end
        for p in 1:np
            @test a.wind24_patch[p] ≈ rv[p] rtol = 1e-13
            @test a.prec24_patch[p] ≈ rp[p] rtol = 1e-13
            @test a.rh24_patch[p] ≈ rr[p] rtol = 1e-13
        end
        # counters cap at the 24-step period
        @test a.fates24_nsteps_patch[1] == 24.0
        @test a.fates24_nsteps_patch[3] == 24.0
        # the late-activated patch started its own mean, so it differs from patch 1
        @test !isapprox(a.wind24_patch[3], a.wind24_patch[1]; rtol = 1e-6)

        # A patch activated 10 steps before the end has a 10-step mean, not a 24-step one.
        b = CLM.Atm2LndData(); CLM.atm2lnd_init!(b, ng, nc, np)
        act = [true, false, false]
        for s in 1:40
            s == 31 && (act[2] = true)
            b.forc_wind_grc[1] = Float64(s)
            b.forc_rain_downscaled_col[1] = 0.0; b.forc_snow_downscaled_col[1] = 0.0
            b.forc_rh_grc[1] = 50.0
            CLM.atm2lnd_update_acc_vars!(b, 1:np, pg, pc; nstep = s, dtime = dtime,
                                         patch_active = act)
        end
        @test b.fates24_nsteps_patch[2] == 10.0
        @test b.wind24_patch[2] ≈ sum(31:40) / 10 rtol = 1e-13
        v = 0.0
        for s in 1:40; n = min(s, 24); v = ((n - 1) * v + s) / n; end
        @test b.wind24_patch[1] ≈ v rtol = 1e-13
        @test b.wind24_patch[3] == 0.0
    finally
        CLM.varctl.use_fates = old_use_fates
    end
end

@testset "non-FATES runs leave the FATES fire-weather accumulators unallocated" begin
    old_use_fates = CLM.varctl.use_fates
    try
        CLM.varctl.use_fates = false
        a = CLM.Atm2LndData(); CLM.atm2lnd_init!(a, 1, 1, 2)
        @test isempty(a.wind24_patch) && isempty(a.prec24_patch) &&
              isempty(a.rh24_patch) && isempty(a.fates24_nsteps_patch)
        CLM.atm2lnd_update_acc_vars!(a, 1:2, [1, 1], [1, 1]; nstep = 3, dtime = 1800,
                                     patch_active = [true, true])
        @test isempty(a.wind24_patch)
    finally
        CLM.varctl.use_fates = old_use_fates
    end
end

@testset "FATES fire data (SPITFIRE modes > 1): lnfm24 accumulator" begin
    d = CLM.FATESFireData()
    @test !CLM.fates_fire_data_active(d)
    CLM.fates_fire_data_init!(d, 2)              # no files: arrays only
    @test CLM.fates_fire_data_active(d)
    @test all(==(0.0), d.lnfm24)                 # init_value = 0
    dtime = 3600
    v = 0.0
    for s in 1:30
        x = 1.0e-3 * (1 + sin(s))
        d.forc_lnfm .= x
        CLM.fates_fire_update_acc_vars!(d, 1:2; dtime = dtime)
        n = min(s, 24); v = ((n - 1) * v + x) / n
    end
    @test d.lnfm24 ≈ [v, v] rtol = 1e-13
    @test d.lnfm24_nsteps == [24.0, 24.0]

    # The synthetic uniform parity stream, if staged locally.
    fl = joinpath(get(ENV, "CESM_INPUTDATA",
                      joinpath(get(ENV, "SYMFLUENCE_DATA", ""), "installs", "cesm-inputdata")),
                  "lnd", "clm2", "firedata", "clmforc.UNIFORM_lnfm_0.9x1.25_CLMjl-parity.nc")
    if isfile(fl)
        e = CLM.FATESFireData()
        CLM.fates_fire_data_init!(e, 1; flnfm = fl)
        grc = CLM.GridcellData(); grc.londeg = [244.4]; grc.latdeg = [51.2]
        for s in 1:24
            CLM.fates_fire_interp!(e, grc; calday = 182.0 + s / 24, dayspyr = 365.0,
                                   bounds_grc = 1:1)
            CLM.fates_fire_update_acc_vars!(e, 1:1; dtime = 3600)
        end
        @test e.forc_lnfm[1] ≈ 0.000457 rtol = 1e-6
        @test e.lnfm24[1] ≈ 0.000457 rtol = 1e-6
    end
end
