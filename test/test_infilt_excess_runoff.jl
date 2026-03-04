@testset "Infiltration Excess Runoff Module" begin
    # ------------------------------------------------------------------
    # Tests for InfiltrationExcessRunoffMod port.
    # Verifies:
    #   1. Module-level constants
    #   2. InfiltrationExcessRunoffParams defaults
    #   3. InfiltrationExcessRunoffData init/clean
    #   4. compute_qinmax_hksat! correctness
    #   5. infiltration_excess_runoff! with QINMAX_METHOD_HKSAT
    #   6. infiltration_excess_runoff! with QINMAX_METHOD_NONE
    #   7. infiltration_excess_runoff! mask-skip behavior
    #   8. Zero ice fraction (maximum infiltration) scenario
    #   9. Fully frozen soil (minimum infiltration) scenario
    #  10. Partial fsat effect on qinmax
    #  11. frac_h2osfc effect on infiltration excess
    #  12. Multiple columns with heterogeneous inputs
    #  13. Invalid qinmax_method raises error
    # ------------------------------------------------------------------

    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
    nlevgrnd = CLM.varpar.nlevgrnd
    nlevsoi  = CLM.varpar.nlevsoi

    # ------------------------------------------------------------------
    # 1. Module-level constants
    # ------------------------------------------------------------------
    @testset "Constants" begin
        @test CLM.QINMAX_UNLIMITED == 1.0e200
        @test CLM.QINMAX_METHOD_NONE == 0
        @test CLM.QINMAX_METHOD_HKSAT == 1
    end

    # ------------------------------------------------------------------
    # 2. InfiltrationExcessRunoffParams defaults
    # ------------------------------------------------------------------
    @testset "Params defaults" begin
        p = CLM.InfiltrationExcessRunoffParams()
        @test p.e_ice == 6.0
        @test CLM.infilt_excess_params.e_ice == 6.0
    end

    # ------------------------------------------------------------------
    # 3. InfiltrationExcessRunoffData init/clean
    # ------------------------------------------------------------------
    @testset "Init and clean" begin
        nc = 5
        ier = CLM.InfiltrationExcessRunoffData()

        # Default method
        @test ier.qinmax_method == CLM.QINMAX_METHOD_HKSAT
        @test isempty(ier.qinmax_col)

        # Init without vichydro -> HKSAT method
        CLM.infilt_excess_runoff_init!(ier, nc; use_vichydro=false)
        @test length(ier.qinmax_col) == nc
        @test all(isnan, ier.qinmax_col)
        @test ier.qinmax_method == CLM.QINMAX_METHOD_HKSAT

        # Init with vichydro -> NONE method
        CLM.infilt_excess_runoff_init!(ier, nc; use_vichydro=true)
        @test ier.qinmax_method == CLM.QINMAX_METHOD_NONE

        # Clean
        CLM.infilt_excess_runoff_clean!(ier)
        @test isempty(ier.qinmax_col)
    end

    # ------------------------------------------------------------------
    # 4. compute_qinmax_hksat! correctness
    # ------------------------------------------------------------------
    @testset "compute_qinmax_hksat!" begin
        nc = 3
        params = CLM.InfiltrationExcessRunoffParams(e_ice=6.0)

        sh = CLM.SoilHydrologyData()
        CLM.soilhydrology_init!(sh, nc)
        # Set ice fractions for layers 1:3
        sh.icefrac_col .= 0.0
        sh.icefrac_col[1, 1] = 0.1
        sh.icefrac_col[1, 2] = 0.2
        sh.icefrac_col[1, 3] = 0.3
        sh.icefrac_col[2, 1] = 0.0
        sh.icefrac_col[2, 2] = 0.0
        sh.icefrac_col[2, 3] = 0.0
        sh.icefrac_col[3, 1] = 0.5
        sh.icefrac_col[3, 2] = 0.5
        sh.icefrac_col[3, 3] = 0.5

        ss = CLM.SoilStateData()
        CLM.soilstate_init!(ss, nc, nc)
        # Set hksat for layers 1:3
        ss.hksat_col[1, 1] = 0.01
        ss.hksat_col[1, 2] = 0.02
        ss.hksat_col[1, 3] = 0.03
        ss.hksat_col[2, 1] = 0.05
        ss.hksat_col[2, 2] = 0.05
        ss.hksat_col[2, 3] = 0.05
        ss.hksat_col[3, 1] = 0.01
        ss.hksat_col[3, 2] = 0.01
        ss.hksat_col[3, 3] = 0.01

        qinmax_unsat = fill(NaN, nc)
        mask = trues(nc)

        CLM.compute_qinmax_hksat!(sh, ss, qinmax_unsat, mask, 1:nc; params=params)

        # Verify: qinmax = min_j=1:3 (10^(-e_ice * icefrac[c,j]) * hksat[c,j])
        # Column 1: min(10^(-6*0.1)*0.01, 10^(-6*0.2)*0.02, 10^(-6*0.3)*0.03)
        v1_1 = 10.0^(-6.0 * 0.1) * 0.01
        v1_2 = 10.0^(-6.0 * 0.2) * 0.02
        v1_3 = 10.0^(-6.0 * 0.3) * 0.03
        expected_1 = min(v1_1, v1_2, v1_3)
        @test qinmax_unsat[1] ≈ expected_1 atol=1e-15

        # Column 2: icefrac=0 everywhere, so 10^0 * hksat = hksat
        # min(0.05, 0.05, 0.05) = 0.05
        @test qinmax_unsat[2] ≈ 0.05 atol=1e-15

        # Column 3: icefrac=0.5 everywhere, hksat=0.01
        # min(10^(-3)*0.01, ...) = 0.01 * 1e-3 = 1e-5
        expected_3 = 10.0^(-6.0 * 0.5) * 0.01
        @test qinmax_unsat[3] ≈ expected_3 atol=1e-15
    end

    # ------------------------------------------------------------------
    # 5. infiltration_excess_runoff! with QINMAX_METHOD_HKSAT
    # ------------------------------------------------------------------
    @testset "infiltration_excess_runoff! HKSAT method" begin
        nc = 2
        params = CLM.InfiltrationExcessRunoffParams(e_ice=6.0)

        ier = CLM.InfiltrationExcessRunoffData()
        CLM.infilt_excess_runoff_init!(ier, nc; use_vichydro=false)

        sh = CLM.SoilHydrologyData()
        CLM.soilhydrology_init!(sh, nc)
        # No ice
        sh.icefrac_col .= 0.0

        ss = CLM.SoilStateData()
        CLM.soilstate_init!(ss, nc, nc)
        # Uniform hksat
        ss.hksat_col .= 0.01  # mm/s

        fsat_col = [0.1, 0.2]

        wfb = CLM.WaterFluxBulkData()
        CLM.waterfluxbulk_init!(wfb, nc, nc, 1, 1)
        wfb.qflx_in_soil_col .= [0.02, 0.005]

        wdb = CLM.WaterDiagnosticBulkData()
        CLM.waterdiagnosticbulk_init!(wdb, nc, nc, 1, 1)
        wdb.frac_h2osfc_col .= [0.0, 0.0]

        mask = trues(nc)

        CLM.infiltration_excess_runoff!(ier, sh, ss, fsat_col, wfb, wdb,
                                         mask, 1:nc; params=params)

        # With no ice: qinmax_unsat = min(hksat(1:3)) = 0.01
        # qinmax[1] = (1 - 0.1) * 0.01 = 0.009
        # qflx_infl_excess[1] = max(0, 0.02 - (1 - 0)*0.009) = max(0, 0.02 - 0.009) = 0.011
        @test ier.qinmax_col[1] ≈ (1.0 - 0.1) * 0.01 atol=1e-15
        @test wfb.qflx_infl_excess_col[1] ≈ max(0.0, 0.02 - (1.0 - 0.0) * 0.009) atol=1e-15

        # qinmax[2] = (1 - 0.2) * 0.01 = 0.008
        # qflx_infl_excess[2] = max(0, 0.005 - (1 - 0)*0.008) = max(0, -0.003) = 0.0
        @test ier.qinmax_col[2] ≈ (1.0 - 0.2) * 0.01 atol=1e-15
        @test wfb.qflx_infl_excess_col[2] ≈ 0.0 atol=1e-15
    end

    # ------------------------------------------------------------------
    # 6. infiltration_excess_runoff! with QINMAX_METHOD_NONE
    # ------------------------------------------------------------------
    @testset "infiltration_excess_runoff! NONE method" begin
        nc = 2

        ier = CLM.InfiltrationExcessRunoffData()
        CLM.infilt_excess_runoff_init!(ier, nc; use_vichydro=true)
        @test ier.qinmax_method == CLM.QINMAX_METHOD_NONE

        sh = CLM.SoilHydrologyData()
        CLM.soilhydrology_init!(sh, nc)

        ss = CLM.SoilStateData()
        CLM.soilstate_init!(ss, nc, nc)

        fsat_col = [0.0, 0.0]

        wfb = CLM.WaterFluxBulkData()
        CLM.waterfluxbulk_init!(wfb, nc, nc, 1, 1)
        wfb.qflx_in_soil_col .= [100.0, 200.0]

        wdb = CLM.WaterDiagnosticBulkData()
        CLM.waterdiagnosticbulk_init!(wdb, nc, nc, 1, 1)
        wdb.frac_h2osfc_col .= [0.0, 0.0]

        mask = trues(nc)

        CLM.infiltration_excess_runoff!(ier, sh, ss, fsat_col, wfb, wdb,
                                         mask, 1:nc)

        # With NONE: qinmax_unsat = 1e200
        # qinmax = (1-0)*1e200 = 1e200
        # qflx_infl_excess = max(0, 100 - 1*1e200) = 0
        for c in 1:nc
            @test ier.qinmax_col[c] ≈ CLM.QINMAX_UNLIMITED atol=1e190
            @test wfb.qflx_infl_excess_col[c] == 0.0
        end
    end

    # ------------------------------------------------------------------
    # 7. Mask-skip behavior
    # ------------------------------------------------------------------
    @testset "Mask-skip behavior" begin
        nc = 3

        ier = CLM.InfiltrationExcessRunoffData()
        CLM.infilt_excess_runoff_init!(ier, nc; use_vichydro=false)

        sh = CLM.SoilHydrologyData()
        CLM.soilhydrology_init!(sh, nc)
        sh.icefrac_col .= 0.0

        ss = CLM.SoilStateData()
        CLM.soilstate_init!(ss, nc, nc)
        ss.hksat_col .= 0.01

        fsat_col = fill(0.0, nc)

        wfb = CLM.WaterFluxBulkData()
        CLM.waterfluxbulk_init!(wfb, nc, nc, 1, 1)
        wfb.qflx_in_soil_col .= 0.02

        wdb = CLM.WaterDiagnosticBulkData()
        CLM.waterdiagnosticbulk_init!(wdb, nc, nc, 1, 1)
        wdb.frac_h2osfc_col .= 0.0

        # Only column 2 is active in the mask
        mask = BitVector([false, true, false])

        CLM.infiltration_excess_runoff!(ier, sh, ss, fsat_col, wfb, wdb,
                                         mask, 1:nc)

        # Column 1 and 3 should remain NaN (untouched)
        @test isnan(ier.qinmax_col[1])
        @test isnan(ier.qinmax_col[3])

        # Column 2 should be computed
        @test !isnan(ier.qinmax_col[2])
        @test ier.qinmax_col[2] ≈ 0.01 atol=1e-15
    end

    # ------------------------------------------------------------------
    # 8. Zero ice fraction (maximum infiltration) scenario
    # ------------------------------------------------------------------
    @testset "Zero ice fraction" begin
        nc = 1
        params = CLM.InfiltrationExcessRunoffParams(e_ice=6.0)

        sh = CLM.SoilHydrologyData()
        CLM.soilhydrology_init!(sh, nc)
        sh.icefrac_col .= 0.0

        ss = CLM.SoilStateData()
        CLM.soilstate_init!(ss, nc, nc)
        hksat_vals = [0.05, 0.03, 0.04]
        for j in 1:3
            ss.hksat_col[1, j] = hksat_vals[j]
        end

        qinmax_unsat = fill(NaN, nc)
        mask = trues(nc)

        CLM.compute_qinmax_hksat!(sh, ss, qinmax_unsat, mask, 1:nc; params=params)

        # With zero ice: 10^0 = 1, so qinmax = min(hksat_vals) = 0.03
        @test qinmax_unsat[1] ≈ 0.03 atol=1e-15
    end

    # ------------------------------------------------------------------
    # 9. Fully frozen soil (minimum infiltration) scenario
    # ------------------------------------------------------------------
    @testset "Fully frozen soil" begin
        nc = 1
        params = CLM.InfiltrationExcessRunoffParams(e_ice=6.0)

        sh = CLM.SoilHydrologyData()
        CLM.soilhydrology_init!(sh, nc)
        # icefrac = 1.0 in top 3 layers
        for j in 1:3
            sh.icefrac_col[1, j] = 1.0
        end

        ss = CLM.SoilStateData()
        CLM.soilstate_init!(ss, nc, nc)
        ss.hksat_col[1, 1] = 0.01
        ss.hksat_col[1, 2] = 0.01
        ss.hksat_col[1, 3] = 0.01

        qinmax_unsat = fill(NaN, nc)
        mask = trues(nc)

        CLM.compute_qinmax_hksat!(sh, ss, qinmax_unsat, mask, 1:nc; params=params)

        # 10^(-6*1) * 0.01 = 1e-6 * 0.01 = 1e-8
        expected = 10.0^(-6.0) * 0.01
        @test qinmax_unsat[1] ≈ expected atol=1e-20

        # Now run full infiltration excess with high input flux
        ier = CLM.InfiltrationExcessRunoffData()
        CLM.infilt_excess_runoff_init!(ier, nc)

        fsat_col = [0.0]

        wfb = CLM.WaterFluxBulkData()
        CLM.waterfluxbulk_init!(wfb, nc, nc, 1, 1)
        wfb.qflx_in_soil_col .= [0.001]  # modest input

        wdb = CLM.WaterDiagnosticBulkData()
        CLM.waterdiagnosticbulk_init!(wdb, nc, nc, 1, 1)
        wdb.frac_h2osfc_col .= [0.0]

        CLM.infiltration_excess_runoff!(ier, sh, ss, fsat_col, wfb, wdb,
                                         mask, 1:nc; params=params)

        # qinmax is very small (1e-8), so nearly all input becomes excess
        @test wfb.qflx_infl_excess_col[1] ≈ max(0.0, 0.001 - expected) atol=1e-12
        @test wfb.qflx_infl_excess_col[1] > 0.0
    end

    # ------------------------------------------------------------------
    # 10. Partial fsat effect on qinmax
    # ------------------------------------------------------------------
    @testset "fsat effect on qinmax" begin
        nc = 1
        params = CLM.InfiltrationExcessRunoffParams(e_ice=6.0)

        ier = CLM.InfiltrationExcessRunoffData()
        CLM.infilt_excess_runoff_init!(ier, nc)

        sh = CLM.SoilHydrologyData()
        CLM.soilhydrology_init!(sh, nc)
        sh.icefrac_col .= 0.0

        ss = CLM.SoilStateData()
        CLM.soilstate_init!(ss, nc, nc)
        ss.hksat_col .= 0.02

        wfb = CLM.WaterFluxBulkData()
        CLM.waterfluxbulk_init!(wfb, nc, nc, 1, 1)
        wfb.qflx_in_soil_col .= [0.005]

        wdb = CLM.WaterDiagnosticBulkData()
        CLM.waterdiagnosticbulk_init!(wdb, nc, nc, 1, 1)
        wdb.frac_h2osfc_col .= [0.0]

        mask = trues(nc)

        # fsat = 0.5 => qinmax = 0.5 * 0.02 = 0.01
        fsat_col = [0.5]
        CLM.infiltration_excess_runoff!(ier, sh, ss, fsat_col, wfb, wdb,
                                         mask, 1:nc; params=params)

        @test ier.qinmax_col[1] ≈ 0.5 * 0.02 atol=1e-15
        # qflx_infl_excess = max(0, 0.005 - 1.0*0.01) = 0
        @test wfb.qflx_infl_excess_col[1] ≈ 0.0 atol=1e-15

        # fsat = 0.9 => qinmax = 0.1 * 0.02 = 0.002
        fsat_col[1] = 0.9
        wfb.qflx_in_soil_col .= [0.005]
        CLM.infiltration_excess_runoff!(ier, sh, ss, fsat_col, wfb, wdb,
                                         mask, 1:nc; params=params)

        @test ier.qinmax_col[1] ≈ 0.1 * 0.02 atol=1e-15
        # qflx_infl_excess = max(0, 0.005 - 1.0*0.002) = 0.003
        @test wfb.qflx_infl_excess_col[1] ≈ 0.003 atol=1e-15
    end

    # ------------------------------------------------------------------
    # 11. frac_h2osfc effect on infiltration excess
    # ------------------------------------------------------------------
    @testset "frac_h2osfc effect" begin
        nc = 1
        params = CLM.InfiltrationExcessRunoffParams(e_ice=6.0)

        ier = CLM.InfiltrationExcessRunoffData()
        CLM.infilt_excess_runoff_init!(ier, nc)

        sh = CLM.SoilHydrologyData()
        CLM.soilhydrology_init!(sh, nc)
        sh.icefrac_col .= 0.0

        ss = CLM.SoilStateData()
        CLM.soilstate_init!(ss, nc, nc)
        ss.hksat_col .= 0.01

        fsat_col = [0.0]

        wfb = CLM.WaterFluxBulkData()
        CLM.waterfluxbulk_init!(wfb, nc, nc, 1, 1)
        wfb.qflx_in_soil_col .= [0.015]

        wdb = CLM.WaterDiagnosticBulkData()
        CLM.waterdiagnosticbulk_init!(wdb, nc, nc, 1, 1)

        mask = trues(nc)

        # frac_h2osfc = 0 => qflx_infl_excess = max(0, 0.015 - 1.0*0.01) = 0.005
        wdb.frac_h2osfc_col .= [0.0]
        CLM.infiltration_excess_runoff!(ier, sh, ss, fsat_col, wfb, wdb,
                                         mask, 1:nc; params=params)
        @test wfb.qflx_infl_excess_col[1] ≈ 0.005 atol=1e-15

        # frac_h2osfc = 0.5 => qflx_infl_excess = max(0, 0.015 - 0.5*0.01) = 0.01
        wfb.qflx_in_soil_col .= [0.015]
        wdb.frac_h2osfc_col .= [0.5]
        CLM.infiltration_excess_runoff!(ier, sh, ss, fsat_col, wfb, wdb,
                                         mask, 1:nc; params=params)
        @test wfb.qflx_infl_excess_col[1] ≈ max(0.0, 0.015 - 0.5 * 0.01) atol=1e-15

        # frac_h2osfc = 1.0 => qflx_infl_excess = max(0, 0.015 - 0.0*0.01) = 0.015
        wfb.qflx_in_soil_col .= [0.015]
        wdb.frac_h2osfc_col .= [1.0]
        CLM.infiltration_excess_runoff!(ier, sh, ss, fsat_col, wfb, wdb,
                                         mask, 1:nc; params=params)
        @test wfb.qflx_infl_excess_col[1] ≈ 0.015 atol=1e-15
    end

    # ------------------------------------------------------------------
    # 12. Multiple columns with heterogeneous inputs
    # ------------------------------------------------------------------
    @testset "Heterogeneous multi-column" begin
        nc = 4
        params = CLM.InfiltrationExcessRunoffParams(e_ice=6.0)

        ier = CLM.InfiltrationExcessRunoffData()
        CLM.infilt_excess_runoff_init!(ier, nc)

        sh = CLM.SoilHydrologyData()
        CLM.soilhydrology_init!(sh, nc)
        sh.icefrac_col .= 0.0
        # Column 3 has ice
        sh.icefrac_col[3, 1] = 0.8
        sh.icefrac_col[3, 2] = 0.8
        sh.icefrac_col[3, 3] = 0.8

        ss = CLM.SoilStateData()
        CLM.soilstate_init!(ss, nc, nc)
        ss.hksat_col .= 0.01
        # Column 2 has high conductivity
        ss.hksat_col[2, 1] = 0.1
        ss.hksat_col[2, 2] = 0.1
        ss.hksat_col[2, 3] = 0.1

        fsat_col = [0.0, 0.0, 0.0, 0.5]

        wfb = CLM.WaterFluxBulkData()
        CLM.waterfluxbulk_init!(wfb, nc, nc, 1, 1)
        wfb.qflx_in_soil_col .= [0.005, 0.005, 0.005, 0.005]

        wdb = CLM.WaterDiagnosticBulkData()
        CLM.waterdiagnosticbulk_init!(wdb, nc, nc, 1, 1)
        wdb.frac_h2osfc_col .= 0.0

        mask = trues(nc)

        CLM.infiltration_excess_runoff!(ier, sh, ss, fsat_col, wfb, wdb,
                                         mask, 1:nc; params=params)

        # Column 1: qinmax=0.01, excess = max(0, 0.005-0.01) = 0
        @test wfb.qflx_infl_excess_col[1] ≈ 0.0 atol=1e-15

        # Column 2: qinmax=0.1, excess = max(0, 0.005-0.1) = 0
        @test wfb.qflx_infl_excess_col[2] ≈ 0.0 atol=1e-15

        # Column 3: qinmax_unsat = 10^(-6*0.8)*0.01 = 10^(-4.8)*0.01
        #           qinmax = 1.0 * qinmax_unsat (very small)
        #           excess ≈ 0.005
        expected_qinmax_3 = 10.0^(-6.0 * 0.8) * 0.01
        @test ier.qinmax_col[3] ≈ expected_qinmax_3 atol=1e-15
        @test wfb.qflx_infl_excess_col[3] ≈ max(0.0, 0.005 - expected_qinmax_3) atol=1e-12

        # Column 4: fsat=0.5, qinmax = 0.5*0.01 = 0.005
        #           excess = max(0, 0.005 - 1.0*0.005) = 0
        @test ier.qinmax_col[4] ≈ 0.005 atol=1e-15
        @test wfb.qflx_infl_excess_col[4] ≈ 0.0 atol=1e-15
    end

    # ------------------------------------------------------------------
    # 13. Invalid qinmax_method raises error
    # ------------------------------------------------------------------
    @testset "Invalid qinmax_method" begin
        nc = 1

        ier = CLM.InfiltrationExcessRunoffData()
        CLM.infilt_excess_runoff_init!(ier, nc)
        ier.qinmax_method = 99  # invalid

        sh = CLM.SoilHydrologyData()
        CLM.soilhydrology_init!(sh, nc)

        ss = CLM.SoilStateData()
        CLM.soilstate_init!(ss, nc, nc)

        fsat_col = [0.0]

        wfb = CLM.WaterFluxBulkData()
        CLM.waterfluxbulk_init!(wfb, nc, nc, 1, 1)
        wfb.qflx_in_soil_col .= [0.01]

        wdb = CLM.WaterDiagnosticBulkData()
        CLM.waterdiagnosticbulk_init!(wdb, nc, nc, 1, 1)
        wdb.frac_h2osfc_col .= [0.0]

        mask = trues(nc)

        @test_throws ErrorException CLM.infiltration_excess_runoff!(
            ier, sh, ss, fsat_col, wfb, wdb, mask, 1:nc)
    end

    # ------------------------------------------------------------------
    # 14. e_ice parameter sensitivity
    # ------------------------------------------------------------------
    @testset "e_ice parameter sensitivity" begin
        nc = 1

        sh = CLM.SoilHydrologyData()
        CLM.soilhydrology_init!(sh, nc)
        sh.icefrac_col .= 0.0
        sh.icefrac_col[1, 1] = 0.5
        sh.icefrac_col[1, 2] = 0.5
        sh.icefrac_col[1, 3] = 0.5

        ss = CLM.SoilStateData()
        CLM.soilstate_init!(ss, nc, nc)
        ss.hksat_col .= 0.01

        mask = trues(nc)

        # Higher e_ice => more ice impedance => lower qinmax
        qinmax_low  = fill(NaN, nc)
        qinmax_high = fill(NaN, nc)

        CLM.compute_qinmax_hksat!(sh, ss, qinmax_low, mask, 1:nc;
            params=CLM.InfiltrationExcessRunoffParams(e_ice=3.0))
        CLM.compute_qinmax_hksat!(sh, ss, qinmax_high, mask, 1:nc;
            params=CLM.InfiltrationExcessRunoffParams(e_ice=9.0))

        @test qinmax_low[1] > qinmax_high[1]
        # Verify exact values
        @test qinmax_low[1]  ≈ 10.0^(-3.0 * 0.5) * 0.01 atol=1e-15
        @test qinmax_high[1] ≈ 10.0^(-9.0 * 0.5) * 0.01 atol=1e-15
    end

end
