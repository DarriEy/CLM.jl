@testset "WaterFluxBulkData" begin

    # Ensure varpar is initialized for dimension parameters
    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)

    nlevsno = CLM.varpar.nlevsno
    nlevsoi = CLM.varpar.nlevsoi

    @testset "default construction" begin
        wfb = CLM.WaterFluxBulkData()

        # Bulk-specific fields
        @test length(wfb.qflx_snowindunload_patch) == 0
        @test length(wfb.qflx_snotempunload_patch) == 0
        @test length(wfb.qflx_phs_neg_col) == 0
        @test length(wfb.qflx_ev_snow_patch) == 0
        @test length(wfb.qflx_ev_snow_col) == 0
        @test length(wfb.qflx_h2osfc_surf_col) == 0
        @test length(wfb.AnnET) == 0
        @test size(wfb.qflx_adv_col) == (0, 0)
        @test size(wfb.qflx_rootsoi_col) == (0, 0)
        @test size(wfb.qflx_snomelt_lyr_col) == (0, 0)
        @test size(wfb.qflx_drain_vr_col) == (0, 0)

        # Parent fields (via composition)
        @test length(wfb.wf.qflx_evap_tot_col) == 0
        @test length(wfb.wf.qflx_tran_veg_patch) == 0
        @test size(wfb.wf.qflx_snofrz_lyr_col) == (0, 0)
    end

    @testset "waterfluxbulk_init!" begin
        nc, np, nl, ng = 10, 15, 3, 2
        wfb = CLM.WaterFluxBulkData()
        CLM.waterfluxbulk_init!(wfb, nc, np, nl, ng)

        # Bulk-specific patch 1D sizes
        @test length(wfb.qflx_snowindunload_patch) == np
        @test length(wfb.qflx_snotempunload_patch) == np
        @test length(wfb.qflx_ev_snow_patch) == np
        @test length(wfb.qflx_ev_soil_patch) == np
        @test length(wfb.qflx_ev_h2osfc_patch) == np
        @test length(wfb.qflx_hydr_redist_patch) == np

        # Bulk-specific column 1D sizes
        @test length(wfb.qflx_phs_neg_col) == nc
        @test length(wfb.qflx_ev_snow_col) == nc
        @test length(wfb.qflx_ev_soil_col) == nc
        @test length(wfb.qflx_ev_h2osfc_col) == nc
        @test length(wfb.qflx_sat_excess_surf_col) == nc
        @test length(wfb.qflx_infl_excess_col) == nc
        @test length(wfb.qflx_infl_excess_surf_col) == nc
        @test length(wfb.qflx_h2osfc_surf_col) == nc
        @test length(wfb.qflx_in_soil_col) == nc
        @test length(wfb.qflx_in_soil_limited_col) == nc
        @test length(wfb.qflx_h2osfc_drain_col) == nc
        @test length(wfb.qflx_top_soil_to_h2osfc_col) == nc
        @test length(wfb.qflx_in_h2osfc_col) == nc
        @test length(wfb.qflx_deficit_col) == nc
        @test length(wfb.AnnET) == nc

        # Bulk-specific column 2D sizes
        @test size(wfb.qflx_adv_col) == (nc, nlevsoi + 1)
        @test size(wfb.qflx_rootsoi_col) == (nc, nlevsoi)
        @test size(wfb.qflx_snomelt_lyr_col) == (nc, nlevsno)
        @test size(wfb.qflx_drain_vr_col) == (nc, nlevsoi)

        # NaN initialization for bulk fields
        @test all(isnan, wfb.qflx_snowindunload_patch)
        @test all(isnan, wfb.qflx_phs_neg_col)
        @test all(isnan, wfb.qflx_ev_snow_col)
        @test all(isnan, wfb.qflx_h2osfc_surf_col)
        @test all(isnan, wfb.AnnET)
        @test all(isnan, wfb.qflx_adv_col)
        @test all(isnan, wfb.qflx_rootsoi_col)
        @test all(isnan, wfb.qflx_snomelt_lyr_col)
        @test all(isnan, wfb.qflx_drain_vr_col)

        # Parent fields should be initialized too
        @test length(wfb.wf.qflx_evap_tot_col) == nc
        @test length(wfb.wf.qflx_tran_veg_patch) == np
        @test length(wfb.wf.qflx_evap_tot_patch) == np
        @test length(wfb.wf.qflx_liq_grnd_col) == nc
        @test length(wfb.wf.qflx_snow_grnd_col) == nc
        @test length(wfb.wf.volumetric_streamflow_lun) == nl
        @test length(wfb.wf.qflx_liq_dynbal_grc) == ng
        @test size(wfb.wf.qflx_snofrz_lyr_col) == (nc, nlevsno)

        # Parent NaN initialization
        @test all(isnan, wfb.wf.qflx_evap_tot_col)
        @test all(isnan, wfb.wf.qflx_tran_veg_patch)
    end

    @testset "waterfluxbulk_clean!" begin
        wfb = CLM.WaterFluxBulkData()
        CLM.waterfluxbulk_init!(wfb, 5, 8, 2, 1)
        CLM.waterfluxbulk_clean!(wfb)

        # Bulk fields
        @test length(wfb.qflx_snowindunload_patch) == 0
        @test length(wfb.qflx_snotempunload_patch) == 0
        @test length(wfb.qflx_ev_snow_patch) == 0
        @test length(wfb.qflx_phs_neg_col) == 0
        @test length(wfb.qflx_ev_snow_col) == 0
        @test length(wfb.qflx_h2osfc_surf_col) == 0
        @test length(wfb.qflx_deficit_col) == 0
        @test length(wfb.AnnET) == 0
        @test size(wfb.qflx_adv_col) == (0, 0)
        @test size(wfb.qflx_rootsoi_col) == (0, 0)
        @test size(wfb.qflx_snomelt_lyr_col) == (0, 0)
        @test size(wfb.qflx_drain_vr_col) == (0, 0)

        # Parent fields should be cleaned too
        @test length(wfb.wf.qflx_evap_tot_col) == 0
        @test length(wfb.wf.qflx_tran_veg_patch) == 0
        @test size(wfb.wf.qflx_snofrz_lyr_col) == (0, 0)
        @test length(wfb.wf.volumetric_streamflow_lun) == 0
    end

    @testset "waterfluxbulk_init_cold!" begin
        nc, np = 4, 6
        wfb = CLM.WaterFluxBulkData()
        CLM.waterfluxbulk_init!(wfb, nc, np, 1, 1)

        CLM.waterfluxbulk_init_cold!(wfb, 1:nc, 1:np)

        # Patch fields should be zero
        @test all(wfb.qflx_snowindunload_patch[1:np] .== 0.0)
        @test all(wfb.qflx_snotempunload_patch[1:np] .== 0.0)

        # Column fields should be zero
        @test all(wfb.qflx_phs_neg_col[1:nc] .== 0.0)
        @test all(wfb.qflx_h2osfc_surf_col[1:nc] .== 0.0)
    end

    @testset "waterfluxbulk_init_cold! subset" begin
        nc, np = 8, 10
        wfb = CLM.WaterFluxBulkData()
        CLM.waterfluxbulk_init!(wfb, nc, np, 1, 1)

        # Only initialize a subset
        CLM.waterfluxbulk_init_cold!(wfb, 2:5, 3:7)

        # Initialized ranges should be zero
        @test all(wfb.qflx_snowindunload_patch[3:7] .== 0.0)
        @test all(wfb.qflx_snotempunload_patch[3:7] .== 0.0)
        @test all(wfb.qflx_phs_neg_col[2:5] .== 0.0)
        @test all(wfb.qflx_h2osfc_surf_col[2:5] .== 0.0)

        # Outside range should still be NaN
        @test isnan(wfb.qflx_snowindunload_patch[1])
        @test isnan(wfb.qflx_snowindunload_patch[8])
        @test isnan(wfb.qflx_phs_neg_col[1])
        @test isnan(wfb.qflx_phs_neg_col[6])
    end

    @testset "stub functions run without error" begin
        wfb = CLM.WaterFluxBulkData()
        CLM.waterfluxbulk_init!(wfb, 5, 8, 2, 1)

        @test CLM.waterfluxbulk_init_history!(wfb, 1:5) === nothing
        @test CLM.waterfluxbulk_restart!(wfb, 1:5) === nothing
        @test CLM.waterfluxbulk_restart!(wfb, 1:5; flag="write") === nothing
        @test CLM.waterfluxbulk_init_acc_buffer!(wfb, 1:5) === nothing
        @test CLM.waterfluxbulk_init_acc_vars!(wfb, 1:5) === nothing
        @test CLM.waterfluxbulk_update_acc_vars!(wfb, 1:5) === nothing
    end

    @testset "field mutability" begin
        wfb = CLM.WaterFluxBulkData()
        CLM.waterfluxbulk_init!(wfb, 3, 4, 2, 1)

        wfb.qflx_snowindunload_patch[1] = 0.005
        @test wfb.qflx_snowindunload_patch[1] == 0.005

        wfb.qflx_phs_neg_col[2] = 1.5e-4
        @test wfb.qflx_phs_neg_col[2] == 1.5e-4

        wfb.AnnET[1] = 3.0e-5
        @test wfb.AnnET[1] == 3.0e-5

        wfb.qflx_adv_col[1, 1] = 0.001
        @test wfb.qflx_adv_col[1, 1] == 0.001

        wfb.qflx_rootsoi_col[2, 3] = 0.002
        @test wfb.qflx_rootsoi_col[2, 3] == 0.002

        # Parent field mutability through composition
        wfb.wf.qflx_evap_tot_col[1] = 42.0
        @test wfb.wf.qflx_evap_tot_col[1] == 42.0
    end

    @testset "re-init overwrites previous state" begin
        wfb = CLM.WaterFluxBulkData()
        CLM.waterfluxbulk_init!(wfb, 3, 4, 2, 1)
        wfb.qflx_phs_neg_col[1] = 999.0
        wfb.AnnET[2] = 888.0

        CLM.waterfluxbulk_init!(wfb, 7, 10, 3, 2)
        @test length(wfb.qflx_phs_neg_col) == 7
        @test all(isnan, wfb.qflx_phs_neg_col)
        @test length(wfb.AnnET) == 7
        @test all(isnan, wfb.AnnET)
        @test length(wfb.qflx_snowindunload_patch) == 10
        @test size(wfb.qflx_adv_col) == (7, nlevsoi + 1)
        @test size(wfb.qflx_rootsoi_col) == (7, nlevsoi)

        # Parent should also be re-initialized
        @test length(wfb.wf.qflx_evap_tot_col) == 7
        @test all(isnan, wfb.wf.qflx_evap_tot_col)
        @test length(wfb.wf.qflx_tran_veg_patch) == 10
        @test length(wfb.wf.volumetric_streamflow_lun) == 3
    end
end
