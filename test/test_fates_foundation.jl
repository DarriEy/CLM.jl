# test_fates_foundation.jl
# Tests for the FATES foundation layer (Tier F, Batch 1).

using Test
using CLM

@testset "FATES foundation" begin

    @testset "FatesConstantsMod" begin
        # kinds / sentinels
        @test CLM.fates_r8 === Float64
        @test CLM.fates_unset_int == -9999
        @test CLM.fates_unset_r8 == -1.0e36
        @test CLM.itrue == 1
        @test CLM.ifalse == 0

        # string lengths
        @test CLM.fates_short_string_length == 32
        @test CLM.fates_long_string_length == 199

        # enumerations / relationships
        @test CLM.N_DBH_BINS == 6
        @test length(CLM.patchfusion_dbhbin_loweredges) == CLM.N_DBH_BINS
        @test CLM.patchfusion_dbhbin_loweredges[1] == 0.0
        @test CLM.N_DIST_TYPES == 4
        @test CLM.n_landuse_cats == length(CLM.is_crop)
        @test CLM.is_crop[CLM.cropland]            # cropland is a crop
        @test !CLM.is_crop[CLM.primaryland]
        @test CLM.leaves_on != CLM.leaves_off

        # unit conversions: internal consistency
        @test CLM.sec_per_day == 86400.0
        @test CLM.days_per_sec * CLM.sec_per_day ≈ 1.0
        @test CLM.ndays_per_year == 365
        @test CLM.mm_per_m * CLM.m_per_mm ≈ 1.0
        @test CLM.mpa_per_pa * CLM.pa_per_mpa ≈ 1.0
        @test CLM.umolC_to_kgC ≈ 12.0e-9

        # physical constants
        @test CLM.t_water_freeze_k_1atm == 273.15
        @test CLM.t_water_freeze_k_triple == 273.16
        @test CLM.rad_per_deg ≈ CLM.pi_const / 180.0
        @test CLM.earth_radius_eq == 6378137.0
    end

    @testset "FatesGlobals" begin
        CLM.FatesGlobalsInit(42, true)
        @test CLM.fates_log() == 42
        @test CLM.fates_global_verbose() == true
        CLM.FatesGlobalsInit(6, false)
        @test CLM.fates_global_verbose() == false

        # string helpers
        @test strip(CLM.I2S(7)) == "7"
        @test occursin("3", CLM.N2S(3.0))

        # endrun maps to error
        @test_throws ErrorException CLM.fates_endrun("boom")

        # warnings increment counters (id 5)
        before = CLM.warn_counts[6]   # slot = id+1
        CLM.FatesWarn("a test warning"; index=5)
        @test CLM.warn_counts[6] == before + 1
    end

    @testset "FatesIntegratorsMod" begin
        @test CLM.max_states == 20

        # Exponential decay: dy/dx = -k y, y(0)=1 -> y(x)=exp(-k x).
        k = 0.5
        decay(Y, Ymask, x, p) = -p[1] .* Y

        # Euler over many small steps
        y = [1.0]; mask = [true]; param = [k]
        x = 0.0; dx = 0.001; nsteps = 1000
        yout = similar(y)
        for _ in 1:nsteps
            CLM.Euler(decay, y, mask, dx, x, param, yout)
            y .= yout
            x += dx
        end
        @test isapprox(y[1], exp(-k * 1.0); atol=1e-3)

        # RKF45 single big step should be far more accurate than Euler
        y2 = [1.0]; yout2 = [0.0]
        opt_dx, l_pass = CLM.RKF45(decay, y2, mask, 1.0, 0.0, 1e-6, param, yout2)
        @test isapprox(yout2[1], exp(-k * 1.0); atol=1e-4)
        @test isa(opt_dx, Real)
        @test isa(l_pass, Bool)

        # A multi-variable system: dy1=-y1, dy2=+y2 -> tests vector handling
        sys(Y, Ymask, x, p) = [-Y[1], Y[2]]
        y3 = [1.0, 1.0]; mask3 = [true, true]; yout3 = [0.0, 0.0]
        CLM.RKF45(sys, y3, mask3, 0.1, 0.0, 1e-8, Float64[], yout3)
        @test isapprox(yout3[1], exp(-0.1); atol=1e-6)
        @test isapprox(yout3[2], exp(0.1); atol=1e-6)
    end

    @testset "FatesUtilsMod" begin
        # FindIndex
        arr = ["alpha", "beta", "gamma"]
        @test CLM.FindIndex(arr, "beta") == 2
        @test CLM.FindIndex(arr, "zeta") == 0

        # check_var_real: clean value -> 0; NaN flagged
        @test CLM.check_var_real(1.0, "ok") == 0
        @test (CLM.check_var_real(NaN, "bad") & 1) == 1

        # GreatCircleDist: zero distance for identical points
        @test isapprox(CLM.GreatCircleDist(0.0, 0.0, 0.0, 0.0), 0.0; atol=1e-6)
        # quarter circumference from equator pole-ward (0,0)->(0,90)
        d = CLM.GreatCircleDist(0.0, 0.0, 0.0, 90.0)
        @test isapprox(d, CLM.earth_radius_eq * CLM.pi_const / 2; rtol=1e-6)

        # GetNeighborDistance wraps GreatCircleDist
        latc = [0.0, 0.0]; lonc = [0.0, 0.0]
        @test isapprox(CLM.GetNeighborDistance(1, 2, latc, lonc), 0.0; atol=1e-6)

        # Quadratic solvers: x^2 - 3x + 2 = (x-1)(x-2) -> roots {1,2}
        r1, r2 = CLM.QuadraticRootsSridharachary(1.0, -3.0, 2.0)
        @test isapprox(sort([r1, r2]), [1.0, 2.0]; atol=1e-12)
        n1, n2 = CLM.QuadraticRootsNSWC(1.0, -3.0, 2.0)
        @test isapprox(sort([n1, n2]), [1.0, 2.0]; atol=1e-10)
        # imaginary roots abort
        @test_throws ErrorException CLM.QuadraticRootsSridharachary(1.0, 0.0, 1.0)
    end

    @testset "FatesRunningMeanMod" begin
        @test CLM.moving_ema_window == 0
        @test CLM.fixed_window == 1

        # EMA: averaging a constant sequence returns the constant.
        def = CLM.rmean_def_type()
        CLM.define!(def, 10.0, 1.0, CLM.moving_ema_window)
        @test def.n_mem == 10
        rm = CLM.rmean_type()
        CLM.InitRMean!(rm, def)             # uninitialized EMA -> NaN
        @test isnan(rm.c_mean)
        for _ in 1:10
            CLM.UpdateRMean!(rm, 5.0)
        end
        @test isapprox(CLM.GetMean(rm), 5.0; atol=1e-12)

        # EMA of a constant value approaches that value from a different start.
        # The per-update weight is capped at 1/n_mem, so convergence is
        # geometric (~ (1-1/n_mem)^k), not exact — assert it gets close.
        rm2 = CLM.rmean_type()
        CLM.InitRMean!(rm2, def; init_value=0.0)
        for _ in 1:50
            CLM.UpdateRMean!(rm2, 3.0)
        end
        @test isapprox(rm2.l_mean, 3.0; atol=1e-2)
        @test rm2.l_mean < 3.0   # approaches from below, never overshoots

        # Fixed window: equal weighting over the window; reports the mean and
        # zeros at window completion.
        fdef = CLM.rmean_def_type()
        CLM.define!(fdef, 4.0, 1.0, CLM.fixed_window)
        @test fdef.n_mem == 4
        fm = CLM.rmean_type()
        CLM.InitRMean!(fm, fdef; init_value=0.0, init_offset=0.0)
        for v in (2.0, 2.0, 2.0, 2.0)   # complete one window of constant 2.0
            CLM.UpdateRMean!(fm, v)
        end
        @test isapprox(fm.l_mean, 2.0; atol=1e-12)
        @test fm.c_index == 0           # zeroed at window completion

        # Fuse two identical EMAs -> unchanged.
        a = CLM.rmean_type(); b = CLM.rmean_type()
        CLM.InitRMean!(a, def; init_value=4.0)
        CLM.InitRMean!(b, def; init_value=4.0)
        CLM.FuseRMean!(a, b, 0.5)
        @test isapprox(a.c_mean, 4.0; atol=1e-12)

        # CopyFromDonor
        c = CLM.rmean_type()
        CLM.InitRMean!(c, def; init_value=0.0)  # gives c a def_type
        CLM.CopyFromDonor!(c, a)
        @test c.c_mean == a.c_mean
    end

    @testset "FatesIODimensionsMod" begin
        @test CLM.dimname_cohort == "cohort"
        @test CLM.dimname_levsoil == "levsoi"

        b = CLM.fates_bounds_type()
        @test b.cohort_begin == 0 && b.cohort_end == 0   # defaults

        dim = CLM.fates_io_dimension_type()
        CLM.Init!(dim, "cohort", 4, 1, 100)
        @test dim.name == "cohort"
        @test dim.lower_bound == 1 && dim.upper_bound == 100
        @test length(dim.clump_lower_bound) == 4
        @test all(dim.clump_lower_bound .== -1)
        CLM.SetThreadBounds!(dim, 2, 10, 20)
        @test dim.clump_lower_bound[2] == 10
        @test dim.clump_upper_bound[2] == 20
    end

    @testset "FatesIOVariableKindMod" begin
        @test CLM.site_r8 == "SI_R8"
        @test CLM.cohort_int == "CO_INT"
        @test CLM.group_dyna_simple == 1

        vk = CLM.fates_io_variable_kind_type()
        CLM.Init!(vk, "SI_R8", 2)
        @test vk.name == "SI_R8"
        @test vk.ndims == 2
        @test length(vk.dimsize) == 2
        @test all(vk.dimsize .== CLM.fates_unset_int)
        @test CLM.is_active(vk) == false
        CLM.set_active!(vk)
        @test CLM.is_active(vk) == true

        # iotype_index finds a kind by name
        vk2 = CLM.fates_io_variable_kind_type(); CLM.Init!(vk2, "CO_R8", 1)
        kinds = [vk, vk2]
        @test CLM.iotype_index("CO_R8", 2, kinds) == 2
        @test_throws ErrorException CLM.iotype_index("NOPE", 2, kinds)
    end

    @testset "FatesParametersInterface" begin
        @test CLM.max_params == 250
        @test CLM.max_dimensions == 2
        @test CLM.dimension_name_pft == "fates_pft"

        fp = CLM.fates_parameters_type()
        CLM.Init!(fp)
        @test CLM.num_params(fp) == 0

        # Register a scalar parameter and round-trip its data.
        CLM.RegisterParameter!(fp, "myscalar", CLM.dimension_shape_scalar, String[])
        @test CLM.num_params(fp) == 1
        idx = CLM.FindIndex(fp, "myscalar")
        @test idx == 1
        CLM.SetDataScalar!(fp, idx, 3.14)
        @test CLM.RetrieveParameterScalar(fp, "myscalar") == 3.14

        # Register a 1D parameter, set its size, store + retrieve.
        CLM.RegisterParameter!(fp, "myvec", CLM.dimension_shape_1d,
                               [CLM.dimension_name_pft])
        @test CLM.num_params(fp) == 2
        i2 = CLM.FindIndex(fp, "myvec")
        CLM.SetDimensionSizes!(fp, false, 1, [CLM.dimension_name_pft], [3])
        CLM.SetData1D!(fp, i2, [1.0, 2.0, 3.0])
        out = zeros(3)
        CLM.RetrieveParameter1D(fp, "myvec", out)
        @test out == [1.0, 2.0, 3.0]
        @test CLM.GetMaxDimensionSize(fp) == 3

        # GetUsedDimensions picks up the pft dimension.
        n_used, used = CLM.GetUsedDimensions(fp, false)
        @test n_used == 1
        @test used[1] == CLM.dimension_name_pft

        # GetMetaData
        nm, shp, szs, nms, hostp = CLM.GetMetaData(fp, i2)
        @test nm == "myvec"
        @test shp == CLM.dimension_shape_1d
    end

    @testset "FatesSynchronizedParamsMod" begin
        # Type constructs; register/receive are no-ops (no shared params).
        sp = CLM.FatesSynchronizedParamsType()
        fp = CLM.fates_parameters_type(); CLM.Init!(fp)
        @test CLM.RegisterParams!(sp, fp) === nothing
        @test CLM.ReceiveParams!(sp, fp) === nothing
        @test isa(CLM.FatesSynchronizedParamsInst, CLM.FatesSynchronizedParamsType)
    end

    # ---------------------------------------------------------------------
    # The running-mean DEFINITION globals must be built from the HOST TIMESTEP.
    #
    # Fortran FatesInterfaceMod.F90:1041-1047 defines every rmean_def with
    # `up_period = hlm_stepsize` and derives `n_mem = nint(mem_period/up_period)`.
    # `InitTimeAveragingGlobals()` was ported and then never called, so patches fell
    # back to hardcoded 1800 s definitions: under a dtime=3600 host the "24-hour"
    # fixed window spanned 48 HOURS, `tveg24` advanced only every second day, and
    # `vegtemp_memory`'s unfilled 0.0 slots counted as cold days -- firing a spurious
    # cold leaf-off on every cold-deciduous PFT in mid-July (D3).
    #
    # These assert the WIRING *and* a non-no-op body: that the window LENGTH IN TIME
    # is one day at whatever dtime the host runs, which is the property that broke.
    # ---------------------------------------------------------------------
    @testset "InitTimeAveragingGlobals is wired to the host timestep" begin
        # InitTimeAveragingGlobals reads the parameter-file timescales for the EMA
        # definitions, so give them finite values (this standalone unit-test file
        # never loads a FATES parameter file -- ed_params() returns fates_unset_r8).
        # The `fixed_24hr`/`ema_24hr` windows below depend only on sec_per_day + dt.
        let ep = CLM.ed_params()
            ep.photo_temp_acclim_timescale  = 30.0    # days
            ep.photo_temp_acclim_thome_time = 1.0     # years
            ep.sdlng_emerg_h2o_timescale    = 15.0    # days
            ep.sdlng_mort_par_timescale     = 15.0    # days
            ep.sdlng2sap_par_timescale      = 73.0    # days
            ep.sdlng_mdd_timescale          = 15.0    # days
        end

        for dt in (1800.0, 3600.0, 900.0)
            CLM.hlm_stepsize[] = dt
            CLM.InitTimeAveragingGlobals()

            # The window must span exactly ONE DAY of model time, not a fixed
            # number of samples. This is the assertion that fails on the bug.
            @test CLM.fixed_24hr.up_period == dt
            @test CLM.fixed_24hr.n_mem * CLM.fixed_24hr.up_period == CLM.sec_per_day
            @test CLM.fixed_24hr.n_mem == round(Int, CLM.sec_per_day / dt)
            @test CLM.fixed_24hr.method == CLM.fixed_window

            @test CLM.ema_24hr.up_period == dt
            @test CLM.ema_24hr.n_mem * CLM.ema_24hr.up_period == CLM.sec_per_day

            # The parameter-timescale means are also stepped at dtime, not daily.
            for d in (CLM.ema_lpa, CLM.ema_longterm, CLM.ema_sdlng_emerg_h2o,
                      CLM.ema_sdlng_mort_par, CLM.ema_sdlng2sap_par, CLM.ema_sdlng_mdd)
                @test d.up_period == dt
                @test d.n_mem > 0
            end
        end

        # A patch initialized after FATES init must USE those globals -- not the
        # module-local `_patch_*` fallbacks. (`_rmdef` selects on n_mem > 0.)
        CLM.hlm_stepsize[] = 3600.0
        CLM.InitTimeAveragingGlobals()
        @test CLM._rmdef(CLM.fixed_24hr, CLM._patch_fixed_24hr) === CLM.fixed_24hr
        @test CLM.fixed_24hr.n_mem == 24            # 24 x 3600 s = one day
        @test CLM._patch_fixed_24hr.n_mem == 48     # the 1800 s fallback, unchanged

        # ...and a caller that never ran FATES init still gets a usable definition.
        undefined = CLM.rmean_def_type()
        @test CLM._rmdef(undefined, CLM._patch_fixed_24hr) === CLM._patch_fixed_24hr
    end
end
