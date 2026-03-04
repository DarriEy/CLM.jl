@testset "Saturated Excess Runoff Module" begin
    # ------------------------------------------------------------------
    # Tests for SaturatedExcessRunoffMod port.
    # Verifies:
    #   1. SaturatedExcessRunoffData initialization
    #   2. compute_fsat_topmodel!  (TOPModel parameterization)
    #   3. compute_fsat_vic!       (VIC parameterization)
    #   4. saturated_excess_runoff! (full orchestrator)
    #   5. crop_fsat_equals_zero flag
    #   6. hillslope_fsat_equals_zero flag
    #   7. Urban flood water routing
    # ------------------------------------------------------------------

    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)

    # ==================================================================
    # 1. Initialization tests
    # ==================================================================
    @testset "SaturatedExcessRunoffData initialization" begin
        ser = CLM.SaturatedExcessRunoffData()
        @test length(ser.fsat_col) == 0
        @test length(ser.fcov_col) == 0
        @test ser.fsat_method == CLM.FSAT_METHOD_TOPMODEL

        # Init with TOPModel (default)
        nc = 5
        CLM.saturated_excess_runoff_init!(ser, nc; use_vichydro=false)
        @test length(ser.fsat_col) == nc
        @test length(ser.fcov_col) == nc
        @test ser.fsat_method == CLM.FSAT_METHOD_TOPMODEL
        @test all(isnan, ser.fsat_col)
        @test all(isnan, ser.fcov_col)

        # Init with VIC
        CLM.saturated_excess_runoff_init!(ser, nc; use_vichydro=true)
        @test ser.fsat_method == CLM.FSAT_METHOD_VIC

        # Clean
        CLM.saturated_excess_runoff_clean!(ser)
        @test length(ser.fsat_col) == 0
        @test length(ser.fcov_col) == 0
    end

    @testset "SaturatedExcessRunoffParams" begin
        # The global param instance should exist
        @test isa(CLM.sat_excess_runoff_params, CLM.SaturatedExcessRunoffParams)
        # fff starts as NaN
        @test isnan(CLM.sat_excess_runoff_params.fff)
    end

    # ==================================================================
    # 2. compute_fsat_topmodel!
    # ==================================================================
    @testset "compute_fsat_topmodel! basic" begin
        nc = 4
        bounds = 1:nc
        mask = BitVector([true, true, true, false])

        # Set up inputs
        frost_table  = [3.0, 1.5, 0.5, 1.0]
        zwt          = [2.0, 2.0, 2.0, 2.0]
        zwt_perched  = [1.0, 1.0, 1.0, 1.0]
        wtfact       = [0.4, 0.4, 0.4, 0.4]
        fff          = 0.5  # decay factor (1/m)

        fsat = fill(NaN, nc)

        CLM.compute_fsat_topmodel!(mask, bounds, frost_table, zwt, zwt_perched,
                                    wtfact, fff, fsat)

        # Column 1: frost_table(3.0) > zwt_perched(1.0) AND frost_table(3.0) <= zwt(2.0)?
        #   3.0 > 1.0 = true, but 3.0 <= 2.0 = false -> use zwt
        #   fsat = 0.4 * exp(-0.5 * 0.5 * 2.0) = 0.4 * exp(-0.5)
        @test fsat[1] ≈ 0.4 * exp(-0.5)

        # Column 2: frost_table(1.5) > zwt_perched(1.0) AND frost_table(1.5) <= zwt(2.0)?
        #   1.5 > 1.0 = true, 1.5 <= 2.0 = true -> use zwt_perched
        #   fsat = 0.4 * exp(-0.5 * 0.5 * 1.0) = 0.4 * exp(-0.25)
        @test fsat[2] ≈ 0.4 * exp(-0.25)

        # Column 3: frost_table(0.5) > zwt_perched(1.0)? = false -> use zwt
        #   fsat = 0.4 * exp(-0.5 * 0.5 * 2.0) = 0.4 * exp(-0.5)
        @test fsat[3] ≈ 0.4 * exp(-0.5)

        # Column 4: mask is false, should remain NaN
        @test isnan(fsat[4])
    end

    @testset "compute_fsat_topmodel! shallow water table" begin
        nc = 2
        bounds = 1:nc
        mask = BitVector([true, true])

        frost_table = [5.0, 5.0]
        zwt         = [0.0, 0.1]
        zwt_perched = [10.0, 10.0]  # perched is deeper, so won't be used
        wtfact      = [0.5, 0.8]
        fff         = 1.0

        fsat = fill(NaN, nc)

        CLM.compute_fsat_topmodel!(mask, bounds, frost_table, zwt, zwt_perched,
                                    wtfact, fff, fsat)

        # Column 1: zwt=0 -> fsat = wtfact * exp(0) = 0.5 * 1.0 = 0.5
        @test fsat[1] ≈ 0.5

        # Column 2: fsat = 0.8 * exp(-0.5 * 1.0 * 0.1) = 0.8 * exp(-0.05)
        @test fsat[2] ≈ 0.8 * exp(-0.05)
    end

    # ==================================================================
    # 3. compute_fsat_vic!
    # ==================================================================
    @testset "compute_fsat_vic! basic" begin
        nc = 3
        bounds = 1:nc
        mask = BitVector([true, true, false])

        b_infil           = [0.5, 2.0, 1.0]
        top_max_moist     = [100.0, 200.0, 100.0]
        top_moist_limited = [50.0, 180.0, 50.0]

        fsat = fill(NaN, nc)

        CLM.compute_fsat_vic!(mask, bounds, b_infil, top_max_moist,
                               top_moist_limited, fsat)

        # Column 1: ex = 0.5/(1+0.5) = 1/3
        #   fsat = 1 - (1 - 50/100)^(1/3) = 1 - 0.5^(1/3)
        ex1 = 0.5 / 1.5
        @test fsat[1] ≈ 1.0 - (1.0 - 50.0/100.0)^ex1

        # Column 2: ex = 2.0/(1+2.0) = 2/3
        #   fsat = 1 - (1 - 180/200)^(2/3) = 1 - 0.1^(2/3)
        ex2 = 2.0 / 3.0
        @test fsat[2] ≈ 1.0 - (1.0 - 180.0/200.0)^ex2

        # Column 3: mask is false, should remain NaN
        @test isnan(fsat[3])
    end

    @testset "compute_fsat_vic! fully saturated" begin
        nc = 1
        bounds = 1:nc
        mask = BitVector([true])

        b_infil           = [1.0]
        top_max_moist     = [100.0]
        top_moist_limited = [100.0]  # fully saturated

        fsat = fill(NaN, nc)

        CLM.compute_fsat_vic!(mask, bounds, b_infil, top_max_moist,
                               top_moist_limited, fsat)

        # fsat = 1 - (1 - 1.0)^ex = 1 - 0 = 1.0
        @test fsat[1] ≈ 1.0
    end

    @testset "compute_fsat_vic! dry soil" begin
        nc = 1
        bounds = 1:nc
        mask = BitVector([true])

        b_infil           = [1.0]
        top_max_moist     = [100.0]
        top_moist_limited = [0.0]  # completely dry

        fsat = fill(NaN, nc)

        CLM.compute_fsat_vic!(mask, bounds, b_infil, top_max_moist,
                               top_moist_limited, fsat)

        # fsat = 1 - (1 - 0)^ex = 1 - 1 = 0.0
        @test fsat[1] ≈ 0.0
    end

    # ==================================================================
    # 4. Full saturated_excess_runoff! with TOPModel
    # ==================================================================
    @testset "saturated_excess_runoff! TOPModel" begin
        nc = 3
        bounds = 1:nc
        mask = BitVector([true, true, true])

        # Set global parameter
        CLM.sat_excess_runoff_params.fff = 0.5

        # Initialize SaturatedExcessRunoffData
        ser = CLM.SaturatedExcessRunoffData()
        CLM.saturated_excess_runoff_init!(ser, nc; use_vichydro=false)

        # Column data
        col = CLM.ColumnData()
        CLM.column_init!(col, nc)
        col.landunit  = [1, 1, 1]
        col.urbpoi    = [false, false, false]
        col.active    = [true, true, true]
        col.is_hillslope_column = [false, false, false]
        col.cold      = fill(CLM.ISPVAL, nc)

        # Landunit data
        lun = CLM.LandunitData()
        CLM.landunit_init!(lun, 1)
        lun.itype = [CLM.ISTSOIL]

        # SoilHydrology data
        sh = CLM.SoilHydrologyData()
        CLM.soilhydrology_init!(sh, nc)
        sh.frost_table_col  = [5.0, 1.5, 0.5]
        sh.zwt_col          = [2.0, 2.0, 2.0]
        sh.zwt_perched_col  = [1.0, 1.0, 1.0]

        # SoilState data
        ss = CLM.SoilStateData()
        CLM.soilstate_init!(ss, 0, nc)
        ss.wtfact_col = [0.4, 0.4, 0.4]

        # WaterFluxBulk data
        wfb = CLM.WaterFluxBulkData()
        CLM.waterfluxbulk_init!(wfb, nc, 0, 0, 0)
        wfb.wf.qflx_rain_plus_snomelt_col = [5.0, 10.0, 2.0]
        wfb.wf.qflx_floodc_col            = [0.0, 0.0, 0.0]

        CLM.saturated_excess_runoff!(ser, mask, bounds, col, lun, sh, ss, wfb)

        # Check fsat values (column 1: use zwt, column 2: use zwt_perched, column 3: use zwt)
        # Col 1: frost(5.0) > perched(1.0) = true, frost(5.0) <= zwt(2.0) = false -> use zwt
        @test ser.fsat_col[1] ≈ 0.4 * exp(-0.5 * 0.5 * 2.0)
        # Col 2: frost(1.5) > perched(1.0) = true, frost(1.5) <= zwt(2.0) = true -> use perched
        @test ser.fsat_col[2] ≈ 0.4 * exp(-0.5 * 0.5 * 1.0)
        # Col 3: frost(0.5) > perched(1.0) = false -> use zwt
        @test ser.fsat_col[3] ≈ 0.4 * exp(-0.5 * 0.5 * 2.0)

        # Check qflx_sat_excess_surf = fsat * qflx_rain_plus_snomelt
        @test wfb.qflx_sat_excess_surf_col[1] ≈ ser.fsat_col[1] * 5.0
        @test wfb.qflx_sat_excess_surf_col[2] ≈ ser.fsat_col[2] * 10.0
        @test wfb.qflx_sat_excess_surf_col[3] ≈ ser.fsat_col[3] * 2.0

        # Check fcov = fsat
        @test ser.fcov_col[1] ≈ ser.fsat_col[1]
        @test ser.fcov_col[2] ≈ ser.fsat_col[2]
        @test ser.fcov_col[3] ≈ ser.fsat_col[3]
    end

    # ==================================================================
    # 5. Full saturated_excess_runoff! with VIC
    # ==================================================================
    @testset "saturated_excess_runoff! VIC" begin
        nc = 2
        bounds = 1:nc
        mask = BitVector([true, true])

        # Initialize SaturatedExcessRunoffData with VIC method
        ser = CLM.SaturatedExcessRunoffData()
        CLM.saturated_excess_runoff_init!(ser, nc; use_vichydro=true)
        @test ser.fsat_method == CLM.FSAT_METHOD_VIC

        # Column data
        col = CLM.ColumnData()
        CLM.column_init!(col, nc)
        col.landunit  = [1, 1]
        col.urbpoi    = [false, false]
        col.active    = [true, true]
        col.is_hillslope_column = [false, false]
        col.cold      = fill(CLM.ISPVAL, nc)

        # Landunit data
        lun = CLM.LandunitData()
        CLM.landunit_init!(lun, 1)
        lun.itype = [CLM.ISTSOIL]

        # SoilHydrology data (VIC fields)
        sh = CLM.SoilHydrologyData()
        CLM.soilhydrology_init!(sh, nc)
        sh.b_infil_col           = [1.0, 2.0]
        sh.top_max_moist_col     = [100.0, 200.0]
        sh.top_moist_limited_col = [50.0, 150.0]

        # SoilState data (not used for VIC but must be allocated)
        ss = CLM.SoilStateData()
        CLM.soilstate_init!(ss, 0, nc)

        # WaterFluxBulk data
        wfb = CLM.WaterFluxBulkData()
        CLM.waterfluxbulk_init!(wfb, nc, 0, 0, 0)
        wfb.wf.qflx_rain_plus_snomelt_col = [8.0, 4.0]
        wfb.wf.qflx_floodc_col            = [0.0, 0.0]

        CLM.saturated_excess_runoff!(ser, mask, bounds, col, lun, sh, ss, wfb)

        # Check fsat values
        ex1 = 1.0 / (1.0 + 1.0)  # = 0.5
        expected_fsat1 = 1.0 - (1.0 - 50.0/100.0)^ex1
        @test ser.fsat_col[1] ≈ expected_fsat1

        ex2 = 2.0 / (1.0 + 2.0)  # = 2/3
        expected_fsat2 = 1.0 - (1.0 - 150.0/200.0)^ex2
        @test ser.fsat_col[2] ≈ expected_fsat2

        # Check runoff
        @test wfb.qflx_sat_excess_surf_col[1] ≈ expected_fsat1 * 8.0
        @test wfb.qflx_sat_excess_surf_col[2] ≈ expected_fsat2 * 4.0
    end

    # ==================================================================
    # 6. crop_fsat_equals_zero flag
    # ==================================================================
    @testset "crop_fsat_equals_zero" begin
        nc = 2
        bounds = 1:nc
        mask = BitVector([true, true])

        CLM.sat_excess_runoff_params.fff = 0.5

        ser = CLM.SaturatedExcessRunoffData()
        CLM.saturated_excess_runoff_init!(ser, nc; use_vichydro=false)

        # Column data: col 1 -> landunit 1 (soil), col 2 -> landunit 2 (crop)
        col = CLM.ColumnData()
        CLM.column_init!(col, nc)
        col.landunit  = [1, 2]
        col.urbpoi    = [false, false]
        col.active    = [true, true]
        col.is_hillslope_column = [false, false]
        col.cold      = fill(CLM.ISPVAL, nc)

        # Landunit data: 2 landunits
        lun = CLM.LandunitData()
        CLM.landunit_init!(lun, 2)
        lun.itype = [CLM.ISTSOIL, CLM.ISTCROP]

        # SoilHydrology data
        sh = CLM.SoilHydrologyData()
        CLM.soilhydrology_init!(sh, nc)
        sh.frost_table_col  = [5.0, 5.0]
        sh.zwt_col          = [1.0, 1.0]
        sh.zwt_perched_col  = [10.0, 10.0]

        # SoilState data
        ss = CLM.SoilStateData()
        CLM.soilstate_init!(ss, 0, nc)
        ss.wtfact_col = [0.4, 0.4]

        # WaterFluxBulk data
        wfb = CLM.WaterFluxBulkData()
        CLM.waterfluxbulk_init!(wfb, nc, 0, 0, 0)
        wfb.wf.qflx_rain_plus_snomelt_col = [10.0, 10.0]
        wfb.wf.qflx_floodc_col            = [0.0, 0.0]

        # With crop_fsat_equals_zero = true
        CLM.saturated_excess_runoff!(ser, mask, bounds, col, lun, sh, ss, wfb;
                                      crop_fsat_equals_zero=true)

        # Column 1 (soil): fsat should be computed normally
        expected_fsat1 = 0.4 * exp(-0.5 * 0.5 * 1.0)
        @test ser.fsat_col[1] ≈ expected_fsat1
        @test wfb.qflx_sat_excess_surf_col[1] ≈ expected_fsat1 * 10.0

        # Column 2 (crop): fsat should be zero
        @test ser.fsat_col[2] ≈ 0.0
        @test wfb.qflx_sat_excess_surf_col[2] ≈ 0.0
    end

    # ==================================================================
    # 7. hillslope_fsat_equals_zero flag
    # ==================================================================
    @testset "hillslope_fsat_equals_zero" begin
        nc = 3
        bounds = 1:nc
        mask = BitVector([true, true, true])

        CLM.sat_excess_runoff_params.fff = 0.5

        ser = CLM.SaturatedExcessRunoffData()
        CLM.saturated_excess_runoff_init!(ser, nc; use_vichydro=false)

        # Column data
        col = CLM.ColumnData()
        CLM.column_init!(col, nc)
        col.landunit  = [1, 1, 1]
        col.urbpoi    = [false, false, false]
        col.active    = [true, true, true]
        # Col 1: not hillslope, Col 2: hillslope upland (cold != ISPVAL),
        # Col 3: hillslope lowland (cold == ISPVAL)
        col.is_hillslope_column = [false, true, true]
        col.cold      = [CLM.ISPVAL, 3, CLM.ISPVAL]

        # Landunit data
        lun = CLM.LandunitData()
        CLM.landunit_init!(lun, 1)
        lun.itype = [CLM.ISTSOIL]

        # SoilHydrology data
        sh = CLM.SoilHydrologyData()
        CLM.soilhydrology_init!(sh, nc)
        sh.frost_table_col  = [5.0, 5.0, 5.0]
        sh.zwt_col          = [1.0, 1.0, 1.0]
        sh.zwt_perched_col  = [10.0, 10.0, 10.0]

        # SoilState data
        ss = CLM.SoilStateData()
        CLM.soilstate_init!(ss, 0, nc)
        ss.wtfact_col = [0.4, 0.4, 0.4]

        # WaterFluxBulk data
        wfb = CLM.WaterFluxBulkData()
        CLM.waterfluxbulk_init!(wfb, nc, 0, 0, 0)
        wfb.wf.qflx_rain_plus_snomelt_col = [10.0, 10.0, 10.0]
        wfb.wf.qflx_floodc_col            = [0.0, 0.0, 0.0]

        CLM.saturated_excess_runoff!(ser, mask, bounds, col, lun, sh, ss, wfb;
                                      hillslope_fsat_equals_zero=true)

        expected_fsat = 0.4 * exp(-0.5 * 0.5 * 1.0)

        # Column 1 (not hillslope): fsat computed normally
        @test ser.fsat_col[1] ≈ expected_fsat

        # Column 2 (hillslope, upland, cold != ISPVAL): fsat zeroed
        @test ser.fsat_col[2] ≈ 0.0
        @test wfb.qflx_sat_excess_surf_col[2] ≈ 0.0

        # Column 3 (hillslope, lowland, cold == ISPVAL): fsat NOT zeroed
        @test ser.fsat_col[3] ≈ expected_fsat
        @test wfb.qflx_sat_excess_surf_col[3] ≈ expected_fsat * 10.0
    end

    # ==================================================================
    # 8. Urban flood water routing
    # ==================================================================
    @testset "urban flood water added to runoff" begin
        nc = 2
        bounds = 1:nc
        mask = BitVector([true, true])

        CLM.sat_excess_runoff_params.fff = 0.5

        ser = CLM.SaturatedExcessRunoffData()
        CLM.saturated_excess_runoff_init!(ser, nc; use_vichydro=false)

        # Column data: col 1 is urban, col 2 is not
        col = CLM.ColumnData()
        CLM.column_init!(col, nc)
        col.landunit  = [1, 1]
        col.urbpoi    = [true, false]
        col.active    = [true, true]
        col.is_hillslope_column = [false, false]
        col.cold      = fill(CLM.ISPVAL, nc)

        lun = CLM.LandunitData()
        CLM.landunit_init!(lun, 1)
        lun.itype = [CLM.ISTSOIL]

        sh = CLM.SoilHydrologyData()
        CLM.soilhydrology_init!(sh, nc)
        sh.frost_table_col  = [5.0, 5.0]
        sh.zwt_col          = [1.0, 1.0]
        sh.zwt_perched_col  = [10.0, 10.0]

        ss = CLM.SoilStateData()
        CLM.soilstate_init!(ss, 0, nc)
        ss.wtfact_col = [0.3, 0.3]

        wfb = CLM.WaterFluxBulkData()
        CLM.waterfluxbulk_init!(wfb, nc, 0, 0, 0)
        wfb.wf.qflx_rain_plus_snomelt_col = [5.0, 5.0]
        wfb.wf.qflx_floodc_col            = [2.0, 2.0]

        CLM.saturated_excess_runoff!(ser, mask, bounds, col, lun, sh, ss, wfb)

        expected_fsat = 0.3 * exp(-0.5 * 0.5 * 1.0)

        # Column 1 (urban): qflx_sat_excess_surf = fsat*rain_plus_snomelt + floodc
        @test wfb.qflx_sat_excess_surf_col[1] ≈ expected_fsat * 5.0 + 2.0

        # Column 2 (non-urban): qflx_sat_excess_surf = fsat*rain_plus_snomelt only
        @test wfb.qflx_sat_excess_surf_col[2] ≈ expected_fsat * 5.0
    end

    # ==================================================================
    # 9. Mask-based filtering
    # ==================================================================
    @testset "mask filtering" begin
        nc = 4
        bounds = 1:nc
        mask = BitVector([true, false, true, false])

        CLM.sat_excess_runoff_params.fff = 1.0

        ser = CLM.SaturatedExcessRunoffData()
        CLM.saturated_excess_runoff_init!(ser, nc; use_vichydro=false)

        col = CLM.ColumnData()
        CLM.column_init!(col, nc)
        col.landunit  = [1, 1, 1, 1]
        col.urbpoi    = [false, false, false, false]
        col.active    = [true, true, true, true]
        col.is_hillslope_column = [false, false, false, false]
        col.cold      = fill(CLM.ISPVAL, nc)

        lun = CLM.LandunitData()
        CLM.landunit_init!(lun, 1)
        lun.itype = [CLM.ISTSOIL]

        sh = CLM.SoilHydrologyData()
        CLM.soilhydrology_init!(sh, nc)
        sh.frost_table_col  = fill(5.0, nc)
        sh.zwt_col          = fill(1.0, nc)
        sh.zwt_perched_col  = fill(10.0, nc)

        ss = CLM.SoilStateData()
        CLM.soilstate_init!(ss, 0, nc)
        ss.wtfact_col = fill(0.5, nc)

        wfb = CLM.WaterFluxBulkData()
        CLM.waterfluxbulk_init!(wfb, nc, 0, 0, 0)
        wfb.wf.qflx_rain_plus_snomelt_col = fill(3.0, nc)
        wfb.wf.qflx_floodc_col            = fill(0.0, nc)

        CLM.saturated_excess_runoff!(ser, mask, bounds, col, lun, sh, ss, wfb)

        expected_fsat = 0.5 * exp(-0.5 * 1.0 * 1.0)

        # Active columns (1 and 3) should be computed
        @test ser.fsat_col[1] ≈ expected_fsat
        @test ser.fsat_col[3] ≈ expected_fsat
        @test wfb.qflx_sat_excess_surf_col[1] ≈ expected_fsat * 3.0
        @test wfb.qflx_sat_excess_surf_col[3] ≈ expected_fsat * 3.0

        # Inactive columns (2 and 4) should remain NaN
        @test isnan(ser.fsat_col[2])
        @test isnan(ser.fsat_col[4])
    end

    # ==================================================================
    # 10. Invalid fsat_method error
    # ==================================================================
    @testset "invalid fsat_method raises error" begin
        nc = 1
        bounds = 1:nc
        mask = BitVector([true])

        ser = CLM.SaturatedExcessRunoffData()
        CLM.saturated_excess_runoff_init!(ser, nc)
        ser.fsat_method = 99  # invalid

        col = CLM.ColumnData()
        CLM.column_init!(col, nc)
        col.landunit  = [1]
        col.urbpoi    = [false]
        col.active    = [true]
        col.is_hillslope_column = [false]
        col.cold      = [CLM.ISPVAL]

        lun = CLM.LandunitData()
        CLM.landunit_init!(lun, 1)
        lun.itype = [CLM.ISTSOIL]

        sh = CLM.SoilHydrologyData()
        CLM.soilhydrology_init!(sh, nc)

        ss = CLM.SoilStateData()
        CLM.soilstate_init!(ss, 0, nc)

        wfb = CLM.WaterFluxBulkData()
        CLM.waterfluxbulk_init!(wfb, nc, 0, 0, 0)

        @test_throws ErrorException CLM.saturated_excess_runoff!(
            ser, mask, bounds, col, lun, sh, ss, wfb)
    end

    # ==================================================================
    # 11. Combined crop + hillslope flags
    # ==================================================================
    @testset "combined crop and hillslope flags" begin
        nc = 3
        bounds = 1:nc
        mask = BitVector([true, true, true])

        CLM.sat_excess_runoff_params.fff = 0.5

        ser = CLM.SaturatedExcessRunoffData()
        CLM.saturated_excess_runoff_init!(ser, nc; use_vichydro=false)

        # Col 1: crop, Col 2: hillslope upland (cold != ISPVAL), Col 3: normal soil
        col = CLM.ColumnData()
        CLM.column_init!(col, nc)
        col.landunit  = [1, 2, 2]
        col.urbpoi    = [false, false, false]
        col.active    = [true, true, true]
        col.is_hillslope_column = [false, true, false]
        col.cold      = [CLM.ISPVAL, 3, CLM.ISPVAL]

        lun = CLM.LandunitData()
        CLM.landunit_init!(lun, 2)
        lun.itype = [CLM.ISTCROP, CLM.ISTSOIL]

        sh = CLM.SoilHydrologyData()
        CLM.soilhydrology_init!(sh, nc)
        sh.frost_table_col  = fill(5.0, nc)
        sh.zwt_col          = fill(1.0, nc)
        sh.zwt_perched_col  = fill(10.0, nc)

        ss = CLM.SoilStateData()
        CLM.soilstate_init!(ss, 0, nc)
        ss.wtfact_col = fill(0.4, nc)

        wfb = CLM.WaterFluxBulkData()
        CLM.waterfluxbulk_init!(wfb, nc, 0, 0, 0)
        wfb.wf.qflx_rain_plus_snomelt_col = fill(10.0, nc)
        wfb.wf.qflx_floodc_col            = fill(0.0, nc)

        CLM.saturated_excess_runoff!(ser, mask, bounds, col, lun, sh, ss, wfb;
                                      crop_fsat_equals_zero=true,
                                      hillslope_fsat_equals_zero=true)

        # Column 1 (crop): fsat zeroed
        @test ser.fsat_col[1] ≈ 0.0
        @test wfb.qflx_sat_excess_surf_col[1] ≈ 0.0

        # Column 2 (hillslope upland): fsat zeroed
        @test ser.fsat_col[2] ≈ 0.0
        @test wfb.qflx_sat_excess_surf_col[2] ≈ 0.0

        # Column 3 (normal soil): fsat computed normally
        expected_fsat = 0.4 * exp(-0.5 * 0.5 * 1.0)
        @test ser.fsat_col[3] ≈ expected_fsat
        @test wfb.qflx_sat_excess_surf_col[3] ≈ expected_fsat * 10.0
    end

    # Reset global param to NaN to avoid side effects on other tests
    CLM.sat_excess_runoff_params.fff = NaN
end
