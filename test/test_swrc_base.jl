@testset "SWRC Base" begin
    # ------------------------------------------------------------------
    # The Fortran SoilWaterRetentionCurveMod defines an abstract type with
    # three deferred procedures.  In Julia this maps to an abstract type
    # plus generic functions.  Tests verify:
    #   1. The abstract type exists and cannot be instantiated directly.
    #   2. A concrete subtype can be created and dispatches correctly.
    #   3. Calling the base functions on an un-implemented subtype errors.
    # ------------------------------------------------------------------

    # 1. Abstract type exists
    @test CLM.SoilWaterRetentionCurve isa DataType
    @test isabstracttype(CLM.SoilWaterRetentionCurve)

    # 2. Define a minimal concrete subtype for testing
    struct TestSWRC <: CLM.SoilWaterRetentionCurve end

    # Calling the base stubs on TestSWRC should error (not implemented)
    ss = CLM.SoilStateData()
    @test_throws ErrorException CLM.soil_hk!(TestSWRC(), 1, 1, 0.5, 1.0, ss)
    @test_throws ErrorException CLM.soil_suction!(TestSWRC(), 1, 1, 0.5, ss)
    @test_throws ErrorException CLM.soil_suction_inverse!(TestSWRC(), 1, 1, -100.0, ss)

    # 3. Define a concrete subtype with actual implementations
    struct ClappHornberg1978 <: CLM.SoilWaterRetentionCurve end

    function CLM.soil_hk!(swrc::ClappHornberg1978, c::Int, j::Int,
                          s::Float64, imped::Float64, soilstate::CLM.SoilStateData)
        # Simplified Clapp-Hornberger: hk = hksat * imped * s^(2b+3)
        hksat = soilstate.hksat_col[c, j]
        bsw   = soilstate.bsw_col[c, j]
        hk    = imped * hksat * s^(2.0 * bsw + 3.0)
        dhkds = imped * hksat * (2.0 * bsw + 3.0) * s^(2.0 * bsw + 2.0)
        return (hk, dhkds)
    end

    function CLM.soil_suction!(swrc::ClappHornberg1978, c::Int, j::Int,
                               s::Float64, soilstate::CLM.SoilStateData)
        sucsat = soilstate.sucsat_col[c, j]
        bsw    = soilstate.bsw_col[c, j]
        smp    = -sucsat * s^(-bsw)
        dsmpds = -bsw * smp / s
        return (smp, dsmpds)
    end

    function CLM.soil_suction_inverse!(swrc::ClappHornberg1978, c::Int, j::Int,
                                       smp_target::Float64, soilstate::CLM.SoilStateData)
        sucsat = soilstate.sucsat_col[c, j]
        bsw    = soilstate.bsw_col[c, j]
        s_target = (-smp_target / sucsat)^(-1.0 / bsw)
        return s_target
    end

    # Set up minimal soil state for one column, one layer
    ss2 = CLM.SoilStateData()
    ss2.hksat_col  = fill(0.0, 1, 1)
    ss2.bsw_col    = fill(0.0, 1, 1)
    ss2.sucsat_col = fill(0.0, 1, 1)

    ss2.hksat_col[1, 1]  = 0.01   # mm/s
    ss2.bsw_col[1, 1]    = 5.0
    ss2.sucsat_col[1, 1] = 100.0  # mm

    ch = ClappHornberg1978()

    # Test soil_hk!
    hk, dhkds = CLM.soil_hk!(ch, 1, 1, 0.5, 1.0, ss2)
    hk_expected = 0.01 * 0.5^(2.0 * 5.0 + 3.0)   # hksat * s^(2b+3)
    @test hk ≈ hk_expected
    dhkds_expected = 0.01 * (2.0 * 5.0 + 3.0) * 0.5^(2.0 * 5.0 + 2.0)
    @test dhkds ≈ dhkds_expected

    # Test soil_suction!
    smp, dsmpds = CLM.soil_suction!(ch, 1, 1, 0.5, ss2)
    smp_expected = -100.0 * 0.5^(-5.0)
    @test smp ≈ smp_expected
    dsmpds_expected = -5.0 * smp_expected / 0.5
    @test dsmpds ≈ dsmpds_expected

    # Test soil_suction_inverse! (roundtrip)
    s_orig = 0.7
    smp_at_s, _ = CLM.soil_suction!(ch, 1, 1, s_orig, ss2)
    s_recovered = CLM.soil_suction_inverse!(ch, 1, 1, smp_at_s, ss2)
    @test s_recovered ≈ s_orig atol=1e-12

    # Test with ice impedance factor
    hk_imped, _ = CLM.soil_hk!(ch, 1, 1, 0.5, 0.1, ss2)
    @test hk_imped ≈ 0.1 * hk_expected

    # Test at full saturation (s=1)
    hk_sat, _ = CLM.soil_hk!(ch, 1, 1, 1.0, 1.0, ss2)
    @test hk_sat ≈ 0.01  # hksat * 1^anything = hksat

    smp_sat, _ = CLM.soil_suction!(ch, 1, 1, 1.0, ss2)
    @test smp_sat ≈ -100.0  # -sucsat * 1^(-b) = -sucsat
end
