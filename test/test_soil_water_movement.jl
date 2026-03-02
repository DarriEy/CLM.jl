@testset "Soil Water Movement" begin
    # ------------------------------------------------------------------
    # Tests for SoilWaterMovementMod port.
    # Verifies:
    #   1. Config construction and defaults
    #   2. use_aquifer_layer logic
    #   3. ice_impedance computation
    #   4. tridiagonal_col! solver
    #   5. baseflow_sink! placeholder
    #   6. compute_RHS_moisture_form! basic check
    #   7. compute_LHS_moisture_form! basic check
    #   8. init_soilwater_movement validation
    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    # 1. Config construction and defaults
    # ------------------------------------------------------------------
    @testset "Config defaults" begin
        cfg = CLM.SoilWaterMovementConfig()
        @test cfg.soilwater_movement_method == CLM.ZENGDECKER_2009
        @test cfg.upper_boundary_condition == CLM.BC_FLUX
        @test cfg.lower_boundary_condition == CLM.BC_AQUIFER
        @test cfg.dtmin == 60.0
        @test cfg.verySmall == 1.0e-8
        @test cfg.xTolerUpper == 1.0e-1
        @test cfg.xTolerLower == 1.0e-2
        @test cfg.expensive == 42
        @test cfg.inexpensive == 1
        @test cfg.flux_calculation == 1
        @test cfg.e_ice == 6.0
    end

    @testset "Config custom" begin
        cfg = CLM.SoilWaterMovementConfig(
            soilwater_movement_method = CLM.MOISTURE_FORM,
            upper_boundary_condition = CLM.BC_FLUX,
            lower_boundary_condition = CLM.BC_FLUX,
            dtmin = 30.0,
            e_ice = 5.0
        )
        @test cfg.soilwater_movement_method == CLM.MOISTURE_FORM
        @test cfg.lower_boundary_condition == CLM.BC_FLUX
        @test cfg.dtmin == 30.0
        @test cfg.e_ice == 5.0
    end

    # ------------------------------------------------------------------
    # 2. use_aquifer_layer logic
    # ------------------------------------------------------------------
    @testset "use_aquifer_layer" begin
        cfg_aq = CLM.SoilWaterMovementConfig(lower_boundary_condition = CLM.BC_AQUIFER)
        @test CLM.use_aquifer_layer(cfg_aq) == true

        cfg_wt = CLM.SoilWaterMovementConfig(
            soilwater_movement_method = CLM.MOISTURE_FORM,
            lower_boundary_condition = CLM.BC_WATERTABLE
        )
        @test CLM.use_aquifer_layer(cfg_wt) == true

        cfg_flux = CLM.SoilWaterMovementConfig(
            soilwater_movement_method = CLM.MOISTURE_FORM,
            lower_boundary_condition = CLM.BC_FLUX
        )
        @test CLM.use_aquifer_layer(cfg_flux) == false

        cfg_zero = CLM.SoilWaterMovementConfig(
            soilwater_movement_method = CLM.MOISTURE_FORM,
            lower_boundary_condition = CLM.BC_ZERO_FLUX
        )
        @test CLM.use_aquifer_layer(cfg_zero) == false
    end

    # ------------------------------------------------------------------
    # 3. ice_impedance computation
    # ------------------------------------------------------------------
    @testset "ice_impedance" begin
        # No ice: impedance = 1.0
        @test CLM.ice_impedance(0.0, 6.0) ≈ 1.0

        # Full ice: impedance = 10^(-e_ice)
        @test CLM.ice_impedance(1.0, 6.0) ≈ 1.0e-6

        # Partial ice
        @test CLM.ice_impedance(0.5, 6.0) ≈ 10.0^(-3.0)

        # Different e_ice
        @test CLM.ice_impedance(0.5, 4.0) ≈ 10.0^(-2.0)
    end

    # ------------------------------------------------------------------
    # 4. tridiagonal_col! solver
    # ------------------------------------------------------------------
    @testset "tridiagonal_col!" begin
        # Simple 3x3 tridiagonal system:
        # [2 -1 0] [x1]   [1]
        # [-1 2 -1] [x2] = [0]
        # [0 -1 2] [x3]   [1]
        n = 3
        a = [0.0, -1.0, -1.0]
        b = [2.0, 2.0, 2.0]
        c = [-1.0, -1.0, 0.0]
        r = [1.0, 0.0, 1.0]
        u = zeros(n)

        CLM.tridiagonal_col!(u, a, b, c, r, 1, 1, n)

        # Verify solution: Au = r
        @test 2.0 * u[1] - u[2] ≈ 1.0 atol=1e-12
        @test -u[1] + 2.0 * u[2] - u[3] ≈ 0.0 atol=1e-12
        @test -u[2] + 2.0 * u[3] ≈ 1.0 atol=1e-12

        # Known solution: x = [1, 1, 1]
        @test u[1] ≈ 1.0 atol=1e-12
        @test u[2] ≈ 1.0 atol=1e-12
        @test u[3] ≈ 1.0 atol=1e-12
    end

    @testset "tridiagonal_col! with jtop offset" begin
        # 4-level system but solving only from level 2 onwards
        n = 4
        a = [0.0, 0.0, -1.0, -1.0]
        b = [1.0, 2.0, 2.0, 2.0]
        c = [0.0, -1.0, -1.0, 0.0]
        r = [0.0, 1.0, 0.0, 1.0]
        u = zeros(n)

        CLM.tridiagonal_col!(u, a, b, c, r, 2, 1, n)

        # Only levels 2-4 should be solved
        @test 2.0 * u[2] - u[3] ≈ 1.0 atol=1e-12
        @test -u[2] + 2.0 * u[3] - u[4] ≈ 0.0 atol=1e-12
        @test -u[3] + 2.0 * u[4] ≈ 1.0 atol=1e-12
    end

    # ------------------------------------------------------------------
    # 5. baseflow_sink! placeholder
    # ------------------------------------------------------------------
    @testset "baseflow_sink!" begin
        nc = 3
        nlevsoi = 5
        bf = ones(nc, nlevsoi)
        mask = trues(nc)

        CLM.baseflow_sink!(bf, mask, nlevsoi)
        @test all(bf .== 0.0)
    end

    # ------------------------------------------------------------------
    # 6. compute_RHS_moisture_form!
    # ------------------------------------------------------------------
    @testset "compute_RHS_moisture_form!" begin
        nlayers = 3
        vert_trans_sink = [0.1, 0.2, 0.05]
        vwc_liq = [0.3, 0.35, 0.25]
        qin = [1.0, 0.8, 0.6]
        qout = [0.8, 0.6, 0.3]
        dt_dz = [0.01, 0.012, 0.008]
        rmx = zeros(nlayers)

        CLM.compute_RHS_moisture_form!(1, nlayers,
            vert_trans_sink, vwc_liq, qin, qout, dt_dz, rmx)

        for j in 1:nlayers
            fluxNet = qin[j] - qout[j] - vert_trans_sink[j]
            @test rmx[j] ≈ -fluxNet * dt_dz[j] atol=1e-15
        end
    end

    # ------------------------------------------------------------------
    # 7. compute_LHS_moisture_form!
    # ------------------------------------------------------------------
    @testset "compute_LHS_moisture_form!" begin
        nlayers = 4
        dt_dz = [0.01, 0.012, 0.008, 0.01]
        dqidw0 = [0.0, 0.5, 0.4, 0.3]
        dqidw1 = [0.1, 0.2, 0.15, 0.12]
        dqodw1 = [0.3, 0.25, 0.2, 0.18]
        dqodw2 = [0.4, 0.35, 0.28, 0.0]
        amx = zeros(nlayers)
        bmx = zeros(nlayers)
        cmx = zeros(nlayers)

        CLM.compute_LHS_moisture_form!(1, nlayers,
            dt_dz, dqidw0, dqidw1, dqodw1, dqodw2, amx, bmx, cmx)

        # Top layer
        @test amx[1] == 0.0
        @test bmx[1] ≈ -1.0 - (-dqidw1[1] + dqodw1[1]) * dt_dz[1]
        @test cmx[1] ≈ -dqodw2[1] * dt_dz[1]

        # Interior layers
        for j in 2:(nlayers-1)
            @test amx[j] ≈ dqidw0[j] * dt_dz[j]
            @test bmx[j] ≈ -1.0 - (-dqidw1[j] + dqodw1[j]) * dt_dz[j]
            @test cmx[j] ≈ -dqodw2[j] * dt_dz[j]
        end

        # Bottom layer
        @test amx[nlayers] ≈ dqidw0[nlayers] * dt_dz[nlayers]
        @test bmx[nlayers] ≈ -1.0 - (-dqidw1[nlayers] + dqodw1[nlayers]) * dt_dz[nlayers]
        @test cmx[nlayers] == 0.0
    end

    # ------------------------------------------------------------------
    # 8. init_soilwater_movement validation
    # ------------------------------------------------------------------
    @testset "init_soilwater_movement" begin
        # Default should work fine
        cfg = CLM.init_soilwater_movement()
        @test cfg.soilwater_movement_method == CLM.ZENGDECKER_2009
        @test cfg.lower_boundary_condition == CLM.BC_AQUIFER

        # ZD09 with non-aquifer should error
        @test_throws ErrorException CLM.init_soilwater_movement(
            soilwater_movement_method = CLM.ZENGDECKER_2009,
            lower_boundary_condition = CLM.BC_FLUX
        )

        # Moisture form with flux BC should work
        cfg2 = CLM.init_soilwater_movement(
            soilwater_movement_method = CLM.MOISTURE_FORM,
            lower_boundary_condition = CLM.BC_FLUX
        )
        @test cfg2.soilwater_movement_method == CLM.MOISTURE_FORM
        @test cfg2.lower_boundary_condition == CLM.BC_FLUX
    end

    # ------------------------------------------------------------------
    # 9. Constant values
    # ------------------------------------------------------------------
    @testset "Constants" begin
        @test CLM.ZENGDECKER_2009 == 0
        @test CLM.MOISTURE_FORM == 1
        @test CLM.MIXED_FORM == 2
        @test CLM.HEAD_FORM == 3

        @test CLM.BC_HEAD == 0
        @test CLM.BC_FLUX == 1
        @test CLM.BC_ZERO_FLUX == 2
        @test CLM.BC_WATERTABLE == 3
        @test CLM.BC_AQUIFER == 4

        @test CLM.M_TO_MM == 1.0e3
    end
end
