@testset "WaterStateBulkData" begin

    # Ensure varpar is initialized for dimension parameters
    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)

    nlevsno        = CLM.varpar.nlevsno
    nlevgrnd       = CLM.varpar.nlevgrnd
    nlevmaxurbgrnd = CLM.varpar.nlevmaxurbgrnd
    nlevtot        = nlevsno + nlevmaxurbgrnd

    @testset "default construction" begin
        wsb = CLM.WaterStateBulkData()

        # Bulk-specific fields
        @test length(wsb.snow_persistence_col) == 0
        @test length(wsb.int_snow_col) == 0

        # Parent fields (via composition)
        @test length(wsb.ws.h2osno_no_layers_col) == 0
        @test length(wsb.ws.h2osfc_col) == 0
        @test length(wsb.ws.wa_col) == 0
        @test size(wsb.ws.h2osoi_liq_col) == (0, 0)
        @test size(wsb.ws.h2osoi_ice_col) == (0, 0)
    end

    @testset "waterstatebulk_init!" begin
        nc, np, nl, ng = 10, 15, 3, 2
        wsb = CLM.WaterStateBulkData()
        CLM.waterstatebulk_init!(wsb, nc, np, nl, ng)

        # Bulk-specific 1D sizes
        @test length(wsb.snow_persistence_col) == nc
        @test length(wsb.int_snow_col) == nc

        # Zero initialization for bulk fields
        @test all(wsb.snow_persistence_col .== 0.0)
        @test all(wsb.int_snow_col .== 0.0)

        # Parent fields should be initialized too
        @test length(wsb.ws.h2osno_no_layers_col) == nc
        @test length(wsb.ws.h2osfc_col) == nc
        @test length(wsb.ws.wa_col) == nc
        @test length(wsb.ws.dynbal_baseline_liq_col) == nc
        @test length(wsb.ws.dynbal_baseline_ice_col) == nc
        @test length(wsb.ws.exice_bulk_init) == nc
        @test size(wsb.ws.h2osoi_liq_col) == (nc, nlevtot)
        @test size(wsb.ws.h2osoi_ice_col) == (nc, nlevtot)
        @test size(wsb.ws.h2osoi_vol_col) == (nc, nlevmaxurbgrnd)
        @test size(wsb.ws.excess_ice_col) == (nc, nlevtot)
        @test length(wsb.ws.snocan_patch) == np
        @test length(wsb.ws.liqcan_patch) == np
        @test length(wsb.ws.stream_water_volume_lun) == nl
        @test size(wsb.ws.h2osoi_vol_prs_grc) == (ng, nlevgrnd)

        # Parent NaN initialization
        @test all(isnan, wsb.ws.h2osno_no_layers_col)
        @test all(isnan, wsb.ws.h2osoi_liq_col)
    end

    @testset "waterstatebulk_clean!" begin
        wsb = CLM.WaterStateBulkData()
        CLM.waterstatebulk_init!(wsb, 5, 8, 2, 1)
        CLM.waterstatebulk_clean!(wsb)

        # Bulk fields
        @test length(wsb.snow_persistence_col) == 0
        @test length(wsb.int_snow_col) == 0

        # Parent fields should be cleaned too
        @test length(wsb.ws.h2osno_no_layers_col) == 0
        @test length(wsb.ws.h2osfc_col) == 0
        @test length(wsb.ws.wa_col) == 0
        @test size(wsb.ws.h2osoi_liq_col) == (0, 0)
        @test size(wsb.ws.h2osoi_ice_col) == (0, 0)
        @test length(wsb.ws.snocan_patch) == 0
        @test length(wsb.ws.stream_water_volume_lun) == 0
    end

    @testset "waterstatebulk_init_cold!" begin
        nc = 4
        wsb = CLM.WaterStateBulkData()
        CLM.waterstatebulk_init!(wsb, nc, nc, 1, 1)

        h2osno_input = [0.0, 10.0, 50.0, 100.0]

        CLM.waterstatebulk_init_cold!(wsb, 1:nc;
            h2osno_input_col=h2osno_input)

        # int_snow should equal h2osno_input
        for c in 1:nc
            @test wsb.int_snow_col[c] ≈ h2osno_input[c]
        end

        # snow_persistence should be zero
        @test all(wsb.snow_persistence_col[1:nc] .== 0.0)
    end

    @testset "waterstatebulk_init_cold! subset of columns" begin
        nc = 6
        wsb = CLM.WaterStateBulkData()
        CLM.waterstatebulk_init!(wsb, nc, nc, 1, 1)

        h2osno_input = [0.0, 5.0, 20.0, 30.0, 80.0, 120.0]

        # Only initialize columns 2:4
        CLM.waterstatebulk_init_cold!(wsb, 2:4;
            h2osno_input_col=h2osno_input)

        # Columns 2:4 should be set
        @test wsb.int_snow_col[2] ≈ 5.0
        @test wsb.int_snow_col[3] ≈ 20.0
        @test wsb.int_snow_col[4] ≈ 30.0
        @test wsb.snow_persistence_col[2] == 0.0
        @test wsb.snow_persistence_col[3] == 0.0
        @test wsb.snow_persistence_col[4] == 0.0

        # Columns outside 2:4 should still be 0.0 (from init)
        @test wsb.int_snow_col[1] == 0.0
        @test wsb.int_snow_col[5] == 0.0
        @test wsb.int_snow_col[6] == 0.0
        @test wsb.snow_persistence_col[1] == 0.0
    end

    @testset "stub functions run without error" begin
        wsb = CLM.WaterStateBulkData()
        CLM.waterstatebulk_init!(wsb, 5, 8, 2, 1)

        @test CLM.waterstatebulk_init_history!(wsb, 1:5) === nothing
        @test CLM.waterstatebulk_restart!(wsb, 1:5) === nothing
        @test CLM.waterstatebulk_restart!(wsb, 1:5; flag="write") === nothing
    end

    @testset "field mutability" begin
        wsb = CLM.WaterStateBulkData()
        CLM.waterstatebulk_init!(wsb, 3, 4, 2, 1)

        wsb.snow_persistence_col[1] = 3600.0
        @test wsb.snow_persistence_col[1] == 3600.0

        wsb.int_snow_col[2] = 250.0
        @test wsb.int_snow_col[2] == 250.0

        # Parent field mutability through composition
        wsb.ws.h2osfc_col[1] = 42.0
        @test wsb.ws.h2osfc_col[1] == 42.0

        wsb.ws.wa_col[2] = 5000.0
        @test wsb.ws.wa_col[2] == 5000.0

        wsb.ws.h2osoi_vol_col[1, 1] = 0.35
        @test wsb.ws.h2osoi_vol_col[1, 1] == 0.35
    end

    @testset "re-init overwrites previous state" begin
        wsb = CLM.WaterStateBulkData()
        CLM.waterstatebulk_init!(wsb, 3, 4, 2, 1)
        wsb.snow_persistence_col[1] = 999.0
        wsb.int_snow_col[2] = 888.0

        CLM.waterstatebulk_init!(wsb, 7, 10, 3, 2)
        @test length(wsb.snow_persistence_col) == 7
        @test all(wsb.snow_persistence_col .== 0.0)
        @test length(wsb.int_snow_col) == 7
        @test all(wsb.int_snow_col .== 0.0)

        # Parent should also be re-initialized
        @test length(wsb.ws.h2osfc_col) == 7
        @test all(isnan, wsb.ws.h2osfc_col)
        @test length(wsb.ws.snocan_patch) == 10
        @test length(wsb.ws.stream_water_volume_lun) == 3
        @test size(wsb.ws.h2osoi_vol_prs_grc) == (2, nlevgrnd)
    end
end
