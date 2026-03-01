@testset "WaterDiagnosticBulkData" begin

    # Ensure varpar is initialized for dimension parameters
    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)

    nlevsno        = CLM.varpar.nlevsno
    nlevsoi        = CLM.varpar.nlevsoi
    nlevgrnd       = CLM.varpar.nlevgrnd
    nlevmaxurbgrnd = CLM.varpar.nlevmaxurbgrnd
    nlevtot_snow_soil = nlevsno + nlevgrnd

    @testset "default construction" begin
        wd = CLM.WaterDiagnosticBulkData()

        # Parent fields
        @test length(wd.snowice_col) == 0
        @test length(wd.snowliq_col) == 0
        @test length(wd.h2ocan_patch) == 0
        @test length(wd.q_ref2m_patch) == 0
        @test length(wd.qaf_lun) == 0
        @test length(wd.tws_grc) == 0

        # Bulk column 1D
        @test length(wd.h2osno_total_col) == 0
        @test length(wd.snow_depth_col) == 0
        @test length(wd.frac_sno_col) == 0
        @test length(wd.wf_col) == 0
        @test length(wd.dqgdT_col) == 0
        @test length(wd.qflx_prec_grnd_col) == 0

        # Bulk column 2D
        @test size(wd.snow_layer_unity_col) == (0, 0)
        @test size(wd.bw_col) == (0, 0)
        @test size(wd.snw_rds_col) == (0, 0)
        @test size(wd.frac_iceold_col) == (0, 0)

        # Bulk patch 1D
        @test length(wd.iwue_ln_patch) == 0
        @test length(wd.fwet_patch) == 0
        @test length(wd.qflx_prec_intr_patch) == 0

        # Landunit 1D
        @test length(wd.stream_water_depth_lun) == 0
    end

    @testset "waterdiagnosticbulk_init!" begin
        nc, np, nl, ng = 10, 15, 3, 2
        wd = CLM.WaterDiagnosticBulkData()
        CLM.waterdiagnosticbulk_init!(wd, nc, np, nl, ng)

        # Parent column 1D sizes
        @test length(wd.snowice_col) == nc
        @test length(wd.snowliq_col) == nc
        @test length(wd.total_plant_stored_h2o_col) == nc
        @test length(wd.h2osoi_liqice_10cm_col) == nc
        @test length(wd.qg_snow_col) == nc
        @test length(wd.qg_soil_col) == nc
        @test length(wd.qg_h2osfc_col) == nc
        @test length(wd.qg_col) == nc

        # Parent patch 1D sizes
        @test length(wd.h2ocan_patch) == np
        @test length(wd.q_ref2m_patch) == np

        # Parent landunit / gridcell
        @test length(wd.qaf_lun) == nl
        @test length(wd.tws_grc) == ng

        # Bulk column 1D sizes
        @test length(wd.h2osno_total_col) == nc
        @test length(wd.snow_depth_col) == nc
        @test length(wd.snow_5day_col) == nc
        @test length(wd.snowdp_col) == nc
        @test length(wd.snomelt_accum_col) == nc
        @test length(wd.h2osoi_liq_tot_col) == nc
        @test length(wd.h2osoi_ice_tot_col) == nc
        @test length(wd.exice_subs_tot_col) == nc
        @test length(wd.exice_vol_tot_col) == nc
        @test length(wd.snw_rds_top_col) == nc
        @test length(wd.h2osno_top_col) == nc
        @test length(wd.sno_liq_top_col) == nc
        @test length(wd.dqgdT_col) == nc
        @test length(wd.frac_sno_col) == nc
        @test length(wd.frac_sno_eff_col) == nc
        @test length(wd.frac_h2osfc_col) == nc
        @test length(wd.frac_h2osfc_nosnow_col) == nc
        @test length(wd.wf_col) == nc
        @test length(wd.wf2_col) == nc
        @test length(wd.qflx_prec_grnd_col) == nc

        # Bulk column 2D sizes
        @test size(wd.snow_layer_unity_col) == (nc, nlevsno)
        @test size(wd.bw_col) == (nc, nlevsno)
        @test size(wd.air_vol_col) == (nc, nlevgrnd)
        @test size(wd.h2osoi_liqvol_col) == (nc, nlevtot_snow_soil)
        @test size(wd.swe_old_col) == (nc, nlevsno)
        @test size(wd.exice_subs_col) == (nc, nlevmaxurbgrnd)
        @test size(wd.exice_vol_col) == (nc, nlevsoi)
        @test size(wd.snw_rds_col) == (nc, nlevsno)
        @test size(wd.frac_iceold_col) == (nc, nlevtot_snow_soil)

        # Bulk patch 1D sizes
        @test length(wd.iwue_ln_patch) == np
        @test length(wd.vpd_ref2m_patch) == np
        @test length(wd.rh_ref2m_patch) == np
        @test length(wd.rh_ref2m_r_patch) == np
        @test length(wd.rh_ref2m_u_patch) == np
        @test length(wd.rh_af_patch) == np
        @test length(wd.rh10_af_patch) == np
        @test length(wd.fwet_patch) == np
        @test length(wd.fcansno_patch) == np
        @test length(wd.fdry_patch) == np
        @test length(wd.qflx_prec_intr_patch) == np

        # Landunit
        @test length(wd.stream_water_depth_lun) == nl

        # NaN initialization for standard fields
        @test all(isnan, wd.h2osno_total_col)
        @test all(isnan, wd.snow_depth_col)
        @test all(isnan, wd.snw_rds_col)
        @test all(isnan, wd.frac_iceold_col)
        @test all(isnan, wd.iwue_ln_patch)
        @test all(isnan, wd.fwet_patch)
        @test all(isnan, wd.stream_water_depth_lun)

        # Zero initialization for exice fields
        @test all(wd.exice_subs_tot_col .== 0.0)
        @test all(wd.exice_vol_tot_col .== 0.0)
        @test all(wd.exice_subs_col .== 0.0)
        @test all(wd.exice_vol_col .== 0.0)

        # SPVAL initialization for rh10_af_patch
        @test all(wd.rh10_af_patch .== CLM.SPVAL)
    end

    @testset "waterdiagnosticbulk_clean!" begin
        wd = CLM.WaterDiagnosticBulkData()
        CLM.waterdiagnosticbulk_init!(wd, 5, 8, 2, 1)
        CLM.waterdiagnosticbulk_clean!(wd)

        # Parent fields
        @test length(wd.snowice_col) == 0
        @test length(wd.h2ocan_patch) == 0
        @test length(wd.qaf_lun) == 0
        @test length(wd.tws_grc) == 0

        # Bulk column 1D
        @test length(wd.h2osno_total_col) == 0
        @test length(wd.snow_depth_col) == 0
        @test length(wd.frac_sno_col) == 0
        @test length(wd.dqgdT_col) == 0

        # Bulk column 2D
        @test size(wd.snw_rds_col) == (0, 0)
        @test size(wd.frac_iceold_col) == (0, 0)

        # Bulk patch
        @test length(wd.iwue_ln_patch) == 0
        @test length(wd.fwet_patch) == 0

        # Landunit
        @test length(wd.stream_water_depth_lun) == 0
    end

    @testset "waterdiagnosticbulk_init_cold! (urban)" begin
        nc, np, nl = 2, 3, 1
        wd = CLM.WaterDiagnosticBulkData()
        CLM.waterdiagnosticbulk_init!(wd, nc, np, nl, 1)

        snow_depth_input = [0.02, 0.10]
        h2osno_input = [5.0, 25.0]
        snl = [0, 0]  # no resolved snow layers
        landunit_col = [1, 1]
        urbpoi = BitVector([true])

        CLM.waterdiagnosticbulk_init_cold!(wd, 1:nc, 1:np;
            snow_depth_input_col=snow_depth_input,
            h2osno_input_col=h2osno_input,
            snl_col=snl,
            landunit_col=landunit_col,
            urbpoi=urbpoi)

        # Urban: frac_sno = min(snow_depth / 0.05, 1)
        @test wd.frac_sno_col[1] ≈ 0.02 / 0.05
        @test wd.frac_sno_col[2] ≈ 1.0  # capped at 1

        # snow_depth should be set
        @test wd.snow_depth_col[1] ≈ 0.02
        @test wd.snow_depth_col[2] ≈ 0.10

        # wf/wf2 should be spval
        @test wd.wf_col[1] == CLM.SPVAL
        @test wd.wf2_col[1] == CLM.SPVAL

        # frac_h2osfc should be zero
        @test wd.frac_h2osfc_col[1] == 0.0

        # Patch zeros
        for p in 1:np
            @test wd.fwet_patch[p] == 0.0
            @test wd.fdry_patch[p] == 0.0
            @test wd.fcansno_patch[p] == 0.0
            @test wd.qflx_prec_intr_patch[p] == 0.0
        end

        # Snow layer unity should be 1
        for c in 1:nc
            for j in 1:nlevsno
                @test wd.snow_layer_unity_col[c, j] == 1.0
            end
        end
    end

    @testset "waterdiagnosticbulk_init_cold! (non-urban, Niu-Yang)" begin
        nc, np, nl = 2, 2, 1
        wd = CLM.WaterDiagnosticBulkData()
        CLM.waterdiagnosticbulk_init!(wd, nc, np, nl, 1)

        snow_depth_input = [0.0, 0.5]
        h2osno_input = [0.0, 100.0]
        snl = [0, -2]  # col 2 has 2 resolved snow layers
        landunit_col = [1, 1]
        urbpoi = BitVector([false])
        zlnd = 0.01
        snw_rds_min = 54.526

        CLM.waterdiagnosticbulk_init_cold!(wd, 1:nc, 1:np;
            snow_depth_input_col=snow_depth_input,
            h2osno_input_col=h2osno_input,
            snl_col=snl,
            landunit_col=landunit_col,
            urbpoi=urbpoi,
            zlnd=zlnd,
            snw_rds_min=snw_rds_min)

        # Column 1: no snow → frac_sno = 0
        @test wd.frac_sno_col[1] == 0.0

        # Column 2: has snow → Niu-Yang formula
        snowbd = min(400.0, h2osno_input[2] / snow_depth_input[2])
        fmelt = (snowbd / 100.0)^1.0
        expected_frac_sno = tanh(snow_depth_input[2] / (2.5 * zlnd * fmelt))
        @test wd.frac_sno_col[2] ≈ expected_frac_sno

        # Column 1: no snow at all → snw_rds = 0, snw_rds_top = spval
        for j in 1:nlevsno
            @test wd.snw_rds_col[1, j] == 0.0
        end
        @test wd.snw_rds_top_col[1] == CLM.SPVAL
        @test wd.sno_liq_top_col[1] == CLM.SPVAL

        # Column 2: snl=-2, so active layers are j=-1 and j=0
        # j=-1 → index nlevsno-1, j=0 → index nlevsno
        # Active layers should have snw_rds_min
        @test wd.snw_rds_col[2, nlevsno] ≈ snw_rds_min      # j=0
        @test wd.snw_rds_col[2, nlevsno - 1] ≈ snw_rds_min  # j=-1
        @test wd.snw_rds_top_col[2] ≈ snw_rds_min

        # Inactive snow layers should be zero
        for j in 1:(nlevsno - 2)
            @test wd.snw_rds_col[2, j] == 0.0
        end
    end

    @testset "waterdiagnosticbulk_init_cold! (no layers but snow exists)" begin
        nc, np, nl = 1, 1, 1
        wd = CLM.WaterDiagnosticBulkData()
        CLM.waterdiagnosticbulk_init!(wd, nc, np, nl, 1)

        snow_depth_input = [0.03]
        h2osno_input = [5.0]
        snl = [0]  # no resolved layers but h2osno > 0
        landunit_col = [1]
        urbpoi = BitVector([false])
        snw_rds_min = 54.526

        CLM.waterdiagnosticbulk_init_cold!(wd, 1:nc, 1:np;
            snow_depth_input_col=snow_depth_input,
            h2osno_input_col=h2osno_input,
            snl_col=snl,
            landunit_col=landunit_col,
            urbpoi=urbpoi,
            snw_rds_min=snw_rds_min)

        # j=0 → index nlevsno should be snw_rds_min
        @test wd.snw_rds_col[1, nlevsno] ≈ snw_rds_min
        # j=-nlevsno+1 to j=-1 → zero
        for j in 1:(nlevsno - 1)
            @test wd.snw_rds_col[1, j] == 0.0
        end
        @test wd.snw_rds_top_col[1] == CLM.SPVAL
        @test wd.sno_liq_top_col[1] == CLM.SPVAL
    end

    @testset "waterdiagnosticbulk_summary!" begin
        nc, np, nl = 3, 3, 1
        wd = CLM.WaterDiagnosticBulkData()
        CLM.waterdiagnosticbulk_init!(wd, nc, np, nl, 1)

        # Initialize exice_subs_col to zero (done by init)
        nlevtot = nlevsno + nlevmaxurbgrnd

        # Set up inputs
        h2osoi_liq = fill(1.0, nc, nlevtot)
        h2osoi_ice = fill(0.5, nc, nlevtot)
        excess_ice = fill(0.0, nc, nlevsoi)

        qflx_int_liq = fill(0.1, np)
        qflx_int_snow = fill(0.05, np)
        qflx_liq_grnd = fill(0.2, nc)
        qflx_snow_grnd = fill(0.3, nc)
        h2osno_total = fill(10.0, nc)

        # dz and zi for soil layers (simple uniform 0.1m layers)
        dz = fill(0.1, nc, nlevsoi)
        zi = zeros(nc, nlevsoi)
        for c in 1:nc
            for j in 1:nlevsoi
                zi[c, j] = j * 0.1  # bottom of each layer
            end
        end

        landunit_col = fill(1, nc)
        urbpoi = BitVector([false])
        lun_itype = [CLM.ISTSOIL]

        mask_soilp = trues(np)
        mask_allc = trues(nc)
        mask_nolakec = trues(nc)

        CLM.waterdiagnosticbulk_summary!(wd, 1:nc, 1:np;
            mask_soilp=mask_soilp,
            mask_allc=mask_allc,
            mask_nolakec=mask_nolakec,
            h2osoi_ice_col=h2osoi_ice,
            h2osoi_liq_col=h2osoi_liq,
            excess_ice_col=excess_ice,
            qflx_intercepted_liq_patch=qflx_int_liq,
            qflx_intercepted_snow_patch=qflx_int_snow,
            qflx_liq_grnd_col=qflx_liq_grnd,
            qflx_snow_grnd_col=qflx_snow_grnd,
            h2osno_total_col=h2osno_total,
            dz_col=dz,
            zi_col=zi,
            landunit_col=landunit_col,
            urbpoi=urbpoi,
            lun_itype=lun_itype)

        # h2osno_total should be copied
        @test wd.h2osno_total_col[1] ≈ 10.0

        # qflx_prec_intr = intercepted liq + snow
        @test wd.qflx_prec_intr_patch[1] ≈ 0.15

        # qflx_prec_grnd = liq_grnd + snow_grnd
        @test wd.qflx_prec_grnd_col[1] ≈ 0.5

        # h2osoi_liq_tot and ice_tot should sum over nlevsoi soil layers
        # soil layer j maps to column index j + nlevsno
        expected_liq_tot = nlevsoi * 1.0
        expected_ice_tot = nlevsoi * 0.5
        @test wd.h2osoi_liq_tot_col[1] ≈ expected_liq_tot
        @test wd.h2osoi_ice_tot_col[1] ≈ expected_ice_tot

        # h2osoi_liqice_10cm: zi ≤ 0.1 for j=1 only (zi[c,1]=0.1)
        # so fracl = 1.0 for j=1, and for j=2, zi[c,2]=0.2 > 0.1 and zi[c,1]=0.1 not < 0.1
        # so only j=1 contributes
        expected_10cm = (1.0 + 0.5) * 1.0
        @test wd.h2osoi_liqice_10cm_col[1] ≈ expected_10cm
    end

    @testset "waterdiagnosticbulk_reset!" begin
        nc = 3
        wd = CLM.WaterDiagnosticBulkData()
        CLM.waterdiagnosticbulk_init!(wd, nc, 4, 1, 1)

        snw_rds_min = 54.526
        CLM.waterdiagnosticbulk_reset!(wd, 2; snw_rds_min=snw_rds_min)

        # j=0 maps to index nlevsno
        @test wd.snw_rds_col[2, nlevsno] ≈ snw_rds_min

        # Other columns should still be NaN (from init)
        @test isnan(wd.snw_rds_col[1, nlevsno])
        @test isnan(wd.snw_rds_col[3, nlevsno])
    end

    @testset "waterdiagnosticbulk_reset_filter!" begin
        nc = 4
        wd = CLM.WaterDiagnosticBulkData()
        CLM.waterdiagnosticbulk_init!(wd, nc, 5, 1, 1)

        mask = BitVector([true, false, true, false])
        snw_rds_min = 54.526

        CLM.waterdiagnosticbulk_reset_filter!(wd, mask, 1:nc;
            snw_rds_min=snw_rds_min)

        @test wd.snw_rds_col[1, nlevsno] ≈ snw_rds_min
        @test isnan(wd.snw_rds_col[2, nlevsno])  # not in mask
        @test wd.snw_rds_col[3, nlevsno] ≈ snw_rds_min
        @test isnan(wd.snw_rds_col[4, nlevsno])  # not in mask
    end

    @testset "stub functions run without error" begin
        wd = CLM.WaterDiagnosticBulkData()
        CLM.waterdiagnosticbulk_init!(wd, 5, 8, 2, 1)

        @test CLM.waterdiagnosticbulk_init_history!(wd, 1:5) === nothing
        @test CLM.waterdiagnosticbulk_restart!(wd, 1:5) === nothing
        @test CLM.waterdiagnosticbulk_init_acc_buffer!(wd, 1:5) === nothing
        @test CLM.waterdiagnosticbulk_init_acc_vars!(wd, 1:5) === nothing
        @test CLM.waterdiagnosticbulk_update_acc_vars!(wd, 1:5) === nothing
    end

    @testset "field mutability" begin
        wd = CLM.WaterDiagnosticBulkData()
        CLM.waterdiagnosticbulk_init!(wd, 3, 4, 2, 1)

        wd.snow_depth_col[1] = 0.5
        @test wd.snow_depth_col[1] == 0.5

        wd.frac_sno_col[2] = 0.75
        @test wd.frac_sno_col[2] == 0.75

        wd.snw_rds_col[1, 1] = 100.0
        @test wd.snw_rds_col[1, 1] == 100.0

        wd.iwue_ln_patch[3] = 42.0
        @test wd.iwue_ln_patch[3] == 42.0

        wd.stream_water_depth_lun[1] = 1.5
        @test wd.stream_water_depth_lun[1] == 1.5
    end

    @testset "re-init overwrites previous state" begin
        wd = CLM.WaterDiagnosticBulkData()
        CLM.waterdiagnosticbulk_init!(wd, 3, 4, 2, 1)
        wd.snow_depth_col[1] = 999.0

        CLM.waterdiagnosticbulk_init!(wd, 7, 10, 3, 2)
        @test length(wd.snow_depth_col) == 7
        @test all(isnan, wd.snow_depth_col)
        @test length(wd.iwue_ln_patch) == 10
        @test size(wd.snw_rds_col) == (7, nlevsno)
        @test length(wd.stream_water_depth_lun) == 3
        @test length(wd.tws_grc) == 2
    end

    @testset "params default values" begin
        params = CLM.WaterDiagnosticBulkParams()
        @test params.zlnd == 0.01
        @test params.snw_rds_min == 54.526
    end
end
