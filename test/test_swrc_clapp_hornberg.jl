@testset "SWRC Clapp-Hornberg 1978" begin
    # ------------------------------------------------------------------
    # Tests for SoilWaterRetentionCurveClappHornberg1978Mod port.
    # Verifies:
    #   1. Type hierarchy and construction
    #   2. soil_hk! computation and derivative
    #   3. soil_suction! computation and derivative
    #   4. soil_suction_inverse! roundtrip consistency
    #   5. Finite-difference derivative checks
    #   6. Edge cases (full saturation, ice impedance)
    # ------------------------------------------------------------------

    # --- Setup: minimal soil state for 2 columns, 2 layers ---
    ss = CLM.SoilStateData()
    ss.hksat_col  = zeros(2, 2)
    ss.bsw_col    = zeros(2, 2)
    ss.sucsat_col = zeros(2, 2)

    # Column 1, layer 1: typical loam
    ss.hksat_col[1, 1]  = 0.01    # mm/s
    ss.bsw_col[1, 1]    = 5.0
    ss.sucsat_col[1, 1] = 100.0   # mm

    # Column 1, layer 2: clay-like
    ss.hksat_col[1, 2]  = 0.002
    ss.bsw_col[1, 2]    = 10.0
    ss.sucsat_col[1, 2] = 300.0

    # Column 2, layer 1: sandy
    ss.hksat_col[2, 1]  = 0.05
    ss.bsw_col[2, 1]    = 3.0
    ss.sucsat_col[2, 1] = 50.0

    # Column 2, layer 2
    ss.hksat_col[2, 2]  = 0.008
    ss.bsw_col[2, 2]    = 7.0
    ss.sucsat_col[2, 2] = 200.0

    swrc = CLM.SoilWaterRetentionCurveClappHornberg1978()

    # ------------------------------------------------------------------
    # 1. Type hierarchy
    # ------------------------------------------------------------------
    @test swrc isa CLM.SoilWaterRetentionCurve
    @test swrc isa CLM.SoilWaterRetentionCurveClappHornberg1978

    # ------------------------------------------------------------------
    # 2. soil_hk! — basic computation
    # ------------------------------------------------------------------
    s = 0.5
    imped = 1.0
    hk, dhkds = CLM.soil_hk!(swrc, 1, 1, s, imped, ss)

    hksat = 0.01
    bsw   = 5.0
    hk_expected = imped * hksat * s^(2.0 * bsw + 3.0)
    @test hk ≈ hk_expected

    dhkds_expected = (2.0 * bsw + 3.0) * hk_expected / s
    @test dhkds ≈ dhkds_expected

    # ------------------------------------------------------------------
    # 3. soil_suction! — basic computation
    # ------------------------------------------------------------------
    smp, dsmpds = CLM.soil_suction!(swrc, 1, 1, s, ss)

    sucsat = 100.0
    smp_expected = -sucsat * s^(-bsw)
    @test smp ≈ smp_expected

    dsmpds_expected = -bsw * smp_expected / s
    @test dsmpds ≈ dsmpds_expected

    # ------------------------------------------------------------------
    # 4. soil_suction_inverse! — roundtrip
    # ------------------------------------------------------------------
    for s_orig in [0.3, 0.5, 0.7, 0.9]
        smp_val, _ = CLM.soil_suction!(swrc, 1, 1, s_orig, ss)
        s_recovered = CLM.soil_suction_inverse!(swrc, 1, 1, smp_val, ss)
        @test s_recovered ≈ s_orig atol=1e-12
    end

    # ------------------------------------------------------------------
    # 5. Finite-difference derivative checks
    # ------------------------------------------------------------------
    eps_fd = 1e-7

    # dhkds via finite difference
    for (c, j, s_test) in [(1, 1, 0.5), (1, 2, 0.4), (2, 1, 0.8), (2, 2, 0.6)]
        hk_p, _ = CLM.soil_hk!(swrc, c, j, s_test + eps_fd, 1.0, ss)
        hk_m, _ = CLM.soil_hk!(swrc, c, j, s_test - eps_fd, 1.0, ss)
        _, dhkds_analytic = CLM.soil_hk!(swrc, c, j, s_test, 1.0, ss)
        dhkds_fd = (hk_p - hk_m) / (2.0 * eps_fd)
        @test dhkds_analytic ≈ dhkds_fd rtol=1e-5
    end

    # dsmpds via finite difference
    for (c, j, s_test) in [(1, 1, 0.5), (1, 2, 0.4), (2, 1, 0.8), (2, 2, 0.6)]
        smp_p, _ = CLM.soil_suction!(swrc, c, j, s_test + eps_fd, ss)
        smp_m, _ = CLM.soil_suction!(swrc, c, j, s_test - eps_fd, ss)
        _, dsmpds_analytic = CLM.soil_suction!(swrc, c, j, s_test, ss)
        dsmpds_fd = (smp_p - smp_m) / (2.0 * eps_fd)
        @test dsmpds_analytic ≈ dsmpds_fd rtol=1e-5
    end

    # ------------------------------------------------------------------
    # 6. Edge cases
    # ------------------------------------------------------------------

    # Full saturation (s=1): hk = hksat, smp = -sucsat
    hk_sat, dhkds_sat = CLM.soil_hk!(swrc, 1, 1, 1.0, 1.0, ss)
    @test hk_sat ≈ 0.01  # hksat * 1^anything = hksat
    @test dhkds_sat ≈ (2.0 * 5.0 + 3.0) * 0.01  # (2b+3)*hk/1

    smp_sat, dsmpds_sat = CLM.soil_suction!(swrc, 1, 1, 1.0, ss)
    @test smp_sat ≈ -100.0  # -sucsat * 1^(-b) = -sucsat

    # Ice impedance factor
    hk_imp, _ = CLM.soil_hk!(swrc, 1, 1, 0.5, 0.1, ss)
    @test hk_imp ≈ 0.1 * hk_expected

    # Different column/layer
    hk2, _ = CLM.soil_hk!(swrc, 2, 1, 0.5, 1.0, ss)
    hk2_expected = 0.05 * 0.5^(2.0 * 3.0 + 3.0)
    @test hk2 ≈ hk2_expected

    smp2, _ = CLM.soil_suction!(swrc, 2, 1, 0.5, ss)
    smp2_expected = -50.0 * 0.5^(-3.0)
    @test smp2 ≈ smp2_expected

    # Inverse on different column/layer
    s_inv = CLM.soil_suction_inverse!(swrc, 2, 1, smp2, ss)
    @test s_inv ≈ 0.5 atol=1e-12
end
