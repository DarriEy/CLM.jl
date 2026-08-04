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
    @testset "use_bedrock requires bc_zero_flux (CTSM SoilWaterMovementMod.F90:181)" begin
        # CTSM endrun's on use_bedrock + a non-zero-flux lower BC. CLM.jl mirrored
        # CTSM's OTHER check (ZD09 requires bc_aquifer) but not this one, and #276
        # measured what the unguarded pair does on the corrected tridiagonal solver:
        # 8,748,338 mm of drainage from a site that received 404 mm of precipitation,
        # with NO balance check firing and NO physics test noticing. Assert the guard
        # FIRES on the bad pair and stays silent on the CLM5 default, so this cannot
        # regress into another silently-wrong configuration.
        _saved = CLM.varctl.use_bedrock
        try
            CLM.varctl.use_bedrock = true
            @test_throws ErrorException CLM.init_soilwater_movement(
                soilwater_movement_method=CLM.ZENGDECKER_2009,
                lower_boundary_condition=CLM.BC_AQUIFER)

            # The CLM5 default pair (bedrock on, zero-flux lower BC) must pass.
            cfg = CLM.init_soilwater_movement(
                soilwater_movement_method=CLM.MOISTURE_FORM,
                lower_boundary_condition=CLM.BC_ZERO_FLUX)
            @test cfg.lower_boundary_condition == CLM.BC_ZERO_FLUX

            # With bedrock off, the aquifer configuration is legal again.
            CLM.varctl.use_bedrock = false
            cfg2 = CLM.init_soilwater_movement(
                soilwater_movement_method=CLM.ZENGDECKER_2009,
                lower_boundary_condition=CLM.BC_AQUIFER)
            @test cfg2.lower_boundary_condition == CLM.BC_AQUIFER
        finally
            CLM.varctl.use_bedrock = _saved
        end
    end

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
    # 8b. Moisture-form solver conserves mass (regression)
    #
    # With a BC_FLUX upper boundary, a BC_ZERO_FLUX lower boundary and no root
    # sink, the ONLY water entering the soil column is qflx_infl, so
    #     sum_j d(h2osoi_liq[j])  ==  qflx_infl * dtime   exactly.
    # This failed by ~28% before the interface-flux / inflow-copy loop split in
    # `_soilwm_moisture_form_kernel!`: fusing them created a loop-carried
    # read-after-write on `qout` that the KernelAbstractions CPU backend
    # vectorizes away, leaving qin[j] == 0 on most layers. Needs several columns
    # (nc > SIMD width) to reproduce.
    # ------------------------------------------------------------------
    @testset "moisture form conserves mass" begin
        nc = 12; nlev = 10
        old_nlevsoi, old_nlevsno = CLM.varpar.nlevsoi, CLM.varpar.nlevsno
        old_nlevgrnd = CLM.varpar.nlevgrnd
        CLM.varpar.nlevsoi = nlev; CLM.varpar.nlevsno = 5; CLM.varpar.nlevgrnd = nlev
        try
            nlevsno = 5; joff = nlevsno; joff_zi = nlevsno + 1
            dtime = 1800.0
            dzv = 0.05
            dz = zeros(nc, nlevsno + nlev); z = zeros(nc, nlevsno + nlev)
            zi = zeros(nc, nlevsno + nlev + 1)
            for c in 1:nc, j in 1:nlev
                dz[c, joff + j] = dzv
                z[c, joff + j]  = (j - 0.5) * dzv
                zi[c, joff_zi + j] = j * dzv
            end
            watsat = fill(0.45, nc, nlev)
            hksat  = fill(0.02, nc, nlev)
            bsw    = fill(5.0, nc, nlev)
            sucsat = fill(100.0, nc, nlev)
            icefrac = zeros(nc, nlev)
            zwt    = fill(nlev * dzv, nc)
            h2osoi_liq = zeros(nc, nlevsno + nlev)
            for c in 1:nc, j in 1:nlev
                # a non-uniform profile so the interior fluxes are nonzero
                h2osoi_liq[c, joff + j] = (0.30 + 0.01 * j) * watsat[c, j] * dzv * 1000.0
            end
            liq0 = copy(h2osoi_liq)
            qflx_infl = fill(1.0e-5, nc)

            scrow() = zeros(nc, nlev)
            scr = CLM.MfScr(; hk = scrow(), smp = scrow(), dhkdw = scrow(),
                dsmpdw = scrow(), imped = scrow(), s2 = scrow(), vwc_liq = scrow(),
                dt_dz = scrow(), qin = scrow(), qout = scrow(), dqidw0 = scrow(),
                dqidw1 = scrow(), dqodw1 = scrow(), dqodw2 = scrow(), dwat = scrow(),
                amx = scrow(), bmx = scrow(), cmx = scrow(), rmx = scrow(),
                gam = scrow(), fluxNet0 = scrow(), fluxNet1 = scrow())
            st = CLM.MfState(; h2osoi_liq = h2osoi_liq, smp_l = zeros(nc, nlev),
                hk_l = zeros(nc, nlev), qcharge = zeros(nc), nsubsteps = zeros(nc))
            mfin = CLM.MfIn(; z = z, zi = zi, dz = dz,
                qflx_rootsoi = zeros(nc, nlev), watsat = watsat, hksat = hksat,
                bsw = bsw, sucsat = sucsat, icefrac = icefrac, zwt = zwt,
                qflx_infl = qflx_infl)

            cfg = CLM.SoilWaterMovementConfig(
                soilwater_movement_method = CLM.MOISTURE_FORM,
                lower_boundary_condition  = CLM.BC_ZERO_FLUX)
            CLM.soilwm_moisture_form_solve!(st, scr, mfin, trues(nc), fill(nlev, nc),
                joff, joff_zi, nlev, dtime, cfg)

            for c in 1:nc
                dliq = sum(h2osoi_liq[c, joff+1:joff+nlev]) - sum(liq0[c, joff+1:joff+nlev])
                @test isapprox(dliq, qflx_infl[c] * dtime; rtol = 1e-10)
                # and the interior inflow/outflow must telescope
                for j in 2:nlev
                    @test scr.qin[c, j] == scr.qout[c, j-1]
                end
                @test scr.qin[c, 1] == qflx_infl[c]
            end
        finally
            CLM.varpar.nlevsoi = old_nlevsoi
            CLM.varpar.nlevsno = old_nlevsno
            CLM.varpar.nlevgrnd = old_nlevgrnd
        end
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
