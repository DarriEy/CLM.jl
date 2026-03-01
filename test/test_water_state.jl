@testset "WaterStateData" begin

    # Ensure varpar is initialized for dimension parameters
    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)

    nlevgrnd       = CLM.varpar.nlevgrnd
    nlevsoi        = CLM.varpar.nlevsoi
    nlevsno        = CLM.varpar.nlevsno
    nlevurb        = CLM.varpar.nlevurb
    nlevmaxurbgrnd = CLM.varpar.nlevmaxurbgrnd
    nlevtot        = nlevsno + nlevmaxurbgrnd

    @testset "default construction" begin
        ws = CLM.WaterStateData()
        @test length(ws.h2osno_no_layers_col) == 0
        @test length(ws.h2osfc_col) == 0
        @test length(ws.wa_col) == 0
        @test length(ws.dynbal_baseline_liq_col) == 0
        @test length(ws.dynbal_baseline_ice_col) == 0
        @test length(ws.exice_bulk_init) == 0
        @test size(ws.h2osoi_liq_col) == (0, 0)
        @test size(ws.h2osoi_ice_col) == (0, 0)
        @test size(ws.h2osoi_vol_col) == (0, 0)
        @test size(ws.excess_ice_col) == (0, 0)
        @test size(ws.h2osoi_vol_prs_grc) == (0, 0)
        @test length(ws.snocan_patch) == 0
        @test length(ws.liqcan_patch) == 0
        @test length(ws.stream_water_volume_lun) == 0
        @test ws.aquifer_water_baseline == 0.0
    end

    @testset "waterstate_init!" begin
        nc, np, nl, ng = 10, 15, 3, 2
        ws = CLM.WaterStateData()
        CLM.waterstate_init!(ws, nc, np, nl, ng)

        # 1D sizes
        @test length(ws.h2osno_no_layers_col) == nc
        @test length(ws.h2osfc_col) == nc
        @test length(ws.wa_col) == nc
        @test length(ws.dynbal_baseline_liq_col) == nc
        @test length(ws.dynbal_baseline_ice_col) == nc
        @test length(ws.exice_bulk_init) == nc
        @test length(ws.snocan_patch) == np
        @test length(ws.liqcan_patch) == np
        @test length(ws.stream_water_volume_lun) == nl

        # 2D sizes
        @test size(ws.h2osoi_liq_col) == (nc, nlevtot)
        @test size(ws.h2osoi_ice_col) == (nc, nlevtot)
        @test size(ws.h2osoi_vol_col) == (nc, nlevmaxurbgrnd)
        @test size(ws.excess_ice_col) == (nc, nlevtot)
        @test size(ws.h2osoi_vol_prs_grc) == (ng, nlevgrnd)

        # NaN initialization
        @test all(isnan, ws.h2osno_no_layers_col)
        @test all(isnan, ws.h2osfc_col)
        @test all(isnan, ws.wa_col)
        @test all(isnan, ws.dynbal_baseline_liq_col)
        @test all(isnan, ws.dynbal_baseline_ice_col)
        @test all(isnan, ws.exice_bulk_init)
        @test all(isnan, ws.h2osoi_liq_col)
        @test all(isnan, ws.h2osoi_ice_col)
        @test all(isnan, ws.h2osoi_vol_col)
        @test all(isnan, ws.excess_ice_col)
        @test all(isnan, ws.h2osoi_vol_prs_grc)
        @test all(isnan, ws.snocan_patch)
        @test all(isnan, ws.liqcan_patch)
        @test all(isnan, ws.stream_water_volume_lun)
    end

    @testset "waterstate_clean!" begin
        ws = CLM.WaterStateData()
        CLM.waterstate_init!(ws, 5, 8, 2, 1)
        CLM.waterstate_clean!(ws)

        @test length(ws.h2osno_no_layers_col) == 0
        @test length(ws.h2osfc_col) == 0
        @test length(ws.wa_col) == 0
        @test length(ws.dynbal_baseline_liq_col) == 0
        @test length(ws.dynbal_baseline_ice_col) == 0
        @test length(ws.exice_bulk_init) == 0
        @test size(ws.h2osoi_liq_col) == (0, 0)
        @test size(ws.h2osoi_ice_col) == (0, 0)
        @test size(ws.h2osoi_vol_col) == (0, 0)
        @test size(ws.excess_ice_col) == (0, 0)
        @test size(ws.h2osoi_vol_prs_grc) == (0, 0)
        @test length(ws.snocan_patch) == 0
        @test length(ws.liqcan_patch) == 0
        @test length(ws.stream_water_volume_lun) == 0
    end

    @testset "waterstate_init_cold! (soil columns)" begin
        nc, np, nl, ng = 3, 3, 1, 1
        ws = CLM.WaterStateData()
        CLM.waterstate_init!(ws, nc, np, nl, ng)

        # Set up column metadata: 3 soil columns on landunit 1
        h2osno_input_col = fill(10.0, nc)
        watsat_col = fill(0.4, nc, nlevmaxurbgrnd)
        # Temp above freezing for all layers (combined snow+soil array)
        t_soisno_col = fill(280.0, nc, nlevtot)
        snl_col = fill(0, nc)  # no snow layers
        dz_col = fill(0.1, nc, nlevtot)
        landunit_col = fill(1, nc)
        lakpoi = falses(1)
        urbpoi = falses(1)
        lun_itype = [CLM.ISTSOIL]
        col_itype = fill(1, nc)
        nbedrock_col = fill(nlevsoi, nc)
        gridcell_col = fill(1, nc)
        exice_init_conc_col = fill(0.0, nc)

        CLM.waterstate_init_cold!(ws, 1:nc, 1:np, 1:nl, 1:ng;
            h2osno_input_col=h2osno_input_col,
            watsat_col=watsat_col,
            t_soisno_col=t_soisno_col,
            use_aquifer_layer=false,
            ratio=1.0,
            snl_col=snl_col,
            dz_col=dz_col,
            landunit_col=landunit_col,
            lakpoi=lakpoi,
            urbpoi=urbpoi,
            lun_itype=lun_itype,
            col_itype=col_itype,
            nbedrock_col=nbedrock_col,
            gridcell_col=gridcell_col,
            exice_coldstart_depth=0.5,
            exice_init_conc_col=exice_init_conc_col)

        # Surface water should be zero
        @test all(ws.h2osfc_col[1:nc] .== 0.0)
        # Canopy water should be zero
        @test all(ws.snocan_patch[1:np] .== 0.0)
        @test all(ws.liqcan_patch[1:np] .== 0.0)

        # h2osoi_vol should be 0.15 for soil layers (<=nlevsoi), clamped by watsat
        for c in 1:nc
            for j in 1:nlevsoi
                @test ws.h2osoi_vol_col[c, j] ≈ min(0.15, watsat_col[c, j])
            end
        end

        # Above freezing → ice should be 0, liq should be dz * denh2o * vol
        for c in 1:nc
            for j in 1:nlevgrnd
                jj = j + nlevsno
                expected_vol = ws.h2osoi_vol_col[c, j]
                if expected_vol != CLM.SPVAL
                    @test ws.h2osoi_ice_col[c, jj] ≈ 0.0
                    @test ws.h2osoi_liq_col[c, jj] ≈ dz_col[c, jj] * CLM.DENH2O * expected_vol
                end
            end
        end

        # Snow (no layers) should equal h2osno_input (snl==0)
        for c in 1:nc
            @test ws.h2osno_no_layers_col[c] ≈ 10.0
        end

        # Aquifer water baseline (use_aquifer_layer=false)
        @test ws.aquifer_water_baseline ≈ CLM.AQUIFER_WATER_BASELINE
        for c in 1:nc
            @test ws.wa_col[c] ≈ CLM.AQUIFER_WATER_BASELINE
        end

        # Dynbal baselines should be zero
        @test all(ws.dynbal_baseline_liq_col[1:nc] .== 0.0)
        @test all(ws.dynbal_baseline_ice_col[1:nc] .== 0.0)
    end

    @testset "waterstate_init_cold! (frozen soil)" begin
        nc, np, nl, ng = 2, 2, 1, 1
        ws = CLM.WaterStateData()
        CLM.waterstate_init!(ws, nc, np, nl, ng)

        h2osno_input_col = fill(0.0, nc)
        watsat_col = fill(0.4, nc, nlevmaxurbgrnd)
        # Temp below freezing
        t_soisno_col = fill(260.0, nc, nlevtot)
        snl_col = fill(0, nc)
        dz_col = fill(0.1, nc, nlevtot)
        landunit_col = fill(1, nc)
        lakpoi = falses(1)
        urbpoi = falses(1)
        lun_itype = [CLM.ISTSOIL]
        col_itype = fill(1, nc)
        nbedrock_col = fill(nlevsoi, nc)
        gridcell_col = fill(1, nc)
        exice_init_conc_col = fill(0.0, nc)

        CLM.waterstate_init_cold!(ws, 1:nc, 1:np, 1:nl, 1:ng;
            h2osno_input_col=h2osno_input_col,
            watsat_col=watsat_col,
            t_soisno_col=t_soisno_col,
            use_aquifer_layer=false,
            ratio=1.0,
            snl_col=snl_col,
            dz_col=dz_col,
            landunit_col=landunit_col,
            lakpoi=lakpoi,
            urbpoi=urbpoi,
            lun_itype=lun_itype,
            col_itype=col_itype,
            nbedrock_col=nbedrock_col,
            gridcell_col=gridcell_col,
            exice_init_conc_col=exice_init_conc_col)

        # Below freezing → liq should be 0, ice = dz * denice * vol
        for c in 1:nc
            for j in 1:nlevsoi
                jj = j + nlevsno
                expected_vol = ws.h2osoi_vol_col[c, j]
                if expected_vol != CLM.SPVAL
                    @test ws.h2osoi_liq_col[c, jj] ≈ 0.0
                    @test ws.h2osoi_ice_col[c, jj] ≈ dz_col[c, jj] * CLM.DENICE * expected_vol
                end
            end
        end
    end

    @testset "waterstate_init_cold! (with aquifer layer)" begin
        nc, np, nl, ng = 2, 2, 1, 1
        ws = CLM.WaterStateData()
        CLM.waterstate_init!(ws, nc, np, nl, ng)

        h2osno_input_col = fill(0.0, nc)
        watsat_col = fill(0.5, nc, nlevmaxurbgrnd)
        t_soisno_col = fill(280.0, nc, nlevtot)
        snl_col = fill(0, nc)
        dz_col = fill(0.1, nc, nlevtot)
        landunit_col = fill(1, nc)
        lakpoi = falses(1)
        urbpoi = falses(1)
        lun_itype = [CLM.ISTSOIL]
        col_itype = fill(1, nc)
        nbedrock_col = fill(nlevsoi, nc)
        gridcell_col = fill(1, nc)
        exice_init_conc_col = fill(0.0, nc)

        CLM.waterstate_init_cold!(ws, 1:nc, 1:np, 1:nl, 1:ng;
            h2osno_input_col=h2osno_input_col,
            watsat_col=watsat_col,
            t_soisno_col=t_soisno_col,
            use_aquifer_layer=true,
            ratio=1.0,
            snl_col=snl_col,
            dz_col=dz_col,
            landunit_col=landunit_col,
            lakpoi=lakpoi,
            urbpoi=urbpoi,
            lun_itype=lun_itype,
            col_itype=col_itype,
            nbedrock_col=nbedrock_col,
            gridcell_col=gridcell_col,
            exice_init_conc_col=exice_init_conc_col)

        # Non-urban, non-lake with use_aquifer_layer: wa = 4000
        for c in 1:nc
            @test ws.wa_col[c] ≈ 4000.0
        end
    end

    @testset "waterstate_calculate_total_h2osno!" begin
        nc, np, nl, ng = 3, 3, 1, 1
        ws = CLM.WaterStateData()
        CLM.waterstate_init!(ws, nc, np, nl, ng)

        # Set up simple state
        ws.h2osno_no_layers_col[1] = 5.0   # no snow layers
        ws.h2osno_no_layers_col[2] = 0.0   # has snow layers
        ws.h2osno_no_layers_col[3] = 10.0  # no snow layers

        snl_col = [0, -2, 0]
        # Set snow layer ice/liq for col 2 (snl=-2, so layers snl+1:0 = -1:0)
        # Julia indices: _snow_idx(-1) = nlevsno-1, _snow_idx(0) = nlevsno
        ws.h2osoi_ice_col[2, nlevsno - 1] = 3.0  # j=-1
        ws.h2osoi_liq_col[2, nlevsno - 1] = 1.0
        ws.h2osoi_ice_col[2, nlevsno]     = 4.0  # j=0
        ws.h2osoi_liq_col[2, nlevsno]     = 2.0

        mask = trues(nc)
        h2osno_total = fill(NaN, nc)

        CLM.waterstate_calculate_total_h2osno!(ws, mask, 1:nc, snl_col, h2osno_total)

        @test h2osno_total[1] ≈ 5.0
        @test h2osno_total[2] ≈ 0.0 + 3.0 + 1.0 + 4.0 + 2.0  # 10.0
        @test h2osno_total[3] ≈ 10.0
    end

    @testset "waterstate_check_snow_consistency! (valid)" begin
        nc = 2
        ws = CLM.WaterStateData()
        CLM.waterstate_init!(ws, nc, 2, 1, 1)

        # Col 1: no snow layers (snl=0), h2osno_no_layers can be non-zero
        # All snow layers must be zero or spval for consistency
        ws.h2osno_no_layers_col[1] = 5.0
        for j in (-nlevsno + 1):0
            jj = j + nlevsno
            ws.h2osoi_ice_col[1, jj] = 0.0
            ws.h2osoi_liq_col[1, jj] = 0.0
        end

        snl_col = [0, -1]

        # Col 2: has snow layers (snl=-1), h2osno_no_layers must be zero
        ws.h2osno_no_layers_col[2] = 0.0
        # Only the active snow layer (j=0) should have water
        ws.h2osoi_ice_col[2, nlevsno] = 3.0  # j=0, active
        ws.h2osoi_liq_col[2, nlevsno] = 1.0

        # Layers above snl should be zero or spval
        for j in (-nlevsno + 1):(snl_col[2])
            jj = j + nlevsno
            ws.h2osoi_ice_col[2, jj] = 0.0
            ws.h2osoi_liq_col[2, jj] = 0.0
        end

        mask = trues(nc)
        # Should not throw
        CLM.waterstate_check_snow_consistency!(ws, mask, 1:nc, snl_col, "test")
    end

    @testset "waterstate_check_snow_consistency! (invalid)" begin
        nc = 1
        ws = CLM.WaterStateData()
        CLM.waterstate_init!(ws, nc, 1, 1, 1)

        # Has snow layers but h2osno_no_layers is non-zero → should error
        ws.h2osno_no_layers_col[1] = 5.0
        snl_col = [-1]

        mask = trues(nc)
        @test_throws ErrorException CLM.waterstate_check_snow_consistency!(
            ws, mask, 1:1, snl_col, "test")
    end

    @testset "stub functions run without error" begin
        ws = CLM.WaterStateData()
        CLM.waterstate_init!(ws, 5, 8, 2, 1)

        @test CLM.waterstate_init_history!(ws, 1:5) === nothing
        @test CLM.waterstate_restart!(ws, 1:5) === nothing
    end

    @testset "field mutability" begin
        ws = CLM.WaterStateData()
        CLM.waterstate_init!(ws, 3, 4, 2, 1)

        ws.h2osfc_col[1] = 42.0
        @test ws.h2osfc_col[1] == 42.0

        ws.wa_col[2] = 5000.0
        @test ws.wa_col[2] == 5000.0

        ws.h2osoi_vol_col[1, 1] = 0.35
        @test ws.h2osoi_vol_col[1, 1] == 0.35

        ws.h2osoi_liq_col[2, 3] = 100.0
        @test ws.h2osoi_liq_col[2, 3] == 100.0

        ws.snocan_patch[1] = 2.5
        @test ws.snocan_patch[1] == 2.5

        ws.stream_water_volume_lun[1] = 99.0
        @test ws.stream_water_volume_lun[1] == 99.0
    end

    @testset "re-init overwrites previous state" begin
        ws = CLM.WaterStateData()
        CLM.waterstate_init!(ws, 3, 4, 2, 1)
        ws.h2osfc_col[1] = 999.0

        CLM.waterstate_init!(ws, 7, 10, 3, 2)
        @test length(ws.h2osfc_col) == 7
        @test all(isnan, ws.h2osfc_col)
        @test size(ws.h2osoi_liq_col) == (7, nlevtot)
        @test length(ws.snocan_patch) == 10
        @test length(ws.stream_water_volume_lun) == 3
        @test size(ws.h2osoi_vol_prs_grc) == (2, nlevgrnd)
    end
end
