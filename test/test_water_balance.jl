@testset "WaterBalanceData" begin

    @testset "default construction" begin
        wb = CLM.WaterBalanceData()

        # Column 1D
        @test length(wb.h2osno_old_col) == 0
        @test length(wb.snow_sources_col) == 0
        @test length(wb.snow_sinks_col) == 0
        @test length(wb.wa_reset_nonconservation_gain_col) == 0
        @test length(wb.begwb_col) == 0
        @test length(wb.endwb_col) == 0
        @test length(wb.errh2o_col) == 0
        @test length(wb.errh2osno_col) == 0

        # Patch 1D
        @test length(wb.errh2o_patch) == 0

        # Gridcell 1D
        @test length(wb.liq1_grc) == 0
        @test length(wb.liq2_grc) == 0
        @test length(wb.ice1_grc) == 0
        @test length(wb.ice2_grc) == 0
        @test length(wb.begwb_grc) == 0
        @test length(wb.endwb_grc) == 0
    end

    @testset "waterbalance_init!" begin
        nc, np, ng = 10, 15, 3
        wb = CLM.WaterBalanceData()
        CLM.waterbalance_init!(wb, nc, np, ng)

        # Column 1D sizes
        @test length(wb.h2osno_old_col) == nc
        @test length(wb.snow_sources_col) == nc
        @test length(wb.snow_sinks_col) == nc
        @test length(wb.wa_reset_nonconservation_gain_col) == nc
        @test length(wb.begwb_col) == nc
        @test length(wb.endwb_col) == nc
        @test length(wb.errh2o_col) == nc
        @test length(wb.errh2osno_col) == nc

        # Patch 1D sizes
        @test length(wb.errh2o_patch) == np

        # Gridcell 1D sizes
        @test length(wb.liq1_grc) == ng
        @test length(wb.liq2_grc) == ng
        @test length(wb.ice1_grc) == ng
        @test length(wb.ice2_grc) == ng
        @test length(wb.begwb_grc) == ng
        @test length(wb.endwb_grc) == ng

        # NaN initialization for most fields
        @test all(isnan, wb.h2osno_old_col)
        @test all(isnan, wb.snow_sources_col)
        @test all(isnan, wb.snow_sinks_col)
        @test all(isnan, wb.begwb_col)
        @test all(isnan, wb.endwb_col)
        @test all(isnan, wb.errh2o_col)
        @test all(isnan, wb.errh2osno_col)
        @test all(isnan, wb.errh2o_patch)
        @test all(isnan, wb.liq1_grc)
        @test all(isnan, wb.liq2_grc)
        @test all(isnan, wb.ice1_grc)
        @test all(isnan, wb.ice2_grc)
        @test all(isnan, wb.begwb_grc)
        @test all(isnan, wb.endwb_grc)

        # InitCold: wa_reset_nonconservation_gain_col set to 0
        @test all(wb.wa_reset_nonconservation_gain_col .== 0.0)
    end

    @testset "waterbalance_clean!" begin
        wb = CLM.WaterBalanceData()
        CLM.waterbalance_init!(wb, 5, 8, 2)
        CLM.waterbalance_clean!(wb)

        # Column 1D
        @test length(wb.h2osno_old_col) == 0
        @test length(wb.snow_sources_col) == 0
        @test length(wb.snow_sinks_col) == 0
        @test length(wb.wa_reset_nonconservation_gain_col) == 0
        @test length(wb.begwb_col) == 0
        @test length(wb.endwb_col) == 0
        @test length(wb.errh2o_col) == 0
        @test length(wb.errh2osno_col) == 0

        # Patch 1D
        @test length(wb.errh2o_patch) == 0

        # Gridcell 1D
        @test length(wb.liq1_grc) == 0
        @test length(wb.liq2_grc) == 0
        @test length(wb.ice1_grc) == 0
        @test length(wb.ice2_grc) == 0
        @test length(wb.begwb_grc) == 0
        @test length(wb.endwb_grc) == 0
    end

    @testset "waterbalance_init_cold!" begin
        nc = 5
        wb = CLM.WaterBalanceData()
        CLM.waterbalance_init!(wb, nc, 3, 2)

        # Manually set to non-zero to verify cold init resets them
        wb.wa_reset_nonconservation_gain_col .= 99.0

        CLM.waterbalance_init_cold!(wb, 1:nc)

        @test all(wb.wa_reset_nonconservation_gain_col .== 0.0)
    end

    @testset "waterbalance_init_cold! partial range" begin
        nc = 5
        wb = CLM.WaterBalanceData()
        CLM.waterbalance_init!(wb, nc, 3, 2)

        wb.wa_reset_nonconservation_gain_col .= 99.0

        # Only reset columns 2:4
        CLM.waterbalance_init_cold!(wb, 2:4)

        @test wb.wa_reset_nonconservation_gain_col[1] == 99.0
        @test wb.wa_reset_nonconservation_gain_col[2] == 0.0
        @test wb.wa_reset_nonconservation_gain_col[3] == 0.0
        @test wb.wa_reset_nonconservation_gain_col[4] == 0.0
        @test wb.wa_reset_nonconservation_gain_col[5] == 99.0
    end

    @testset "stub functions run without error" begin
        wb = CLM.WaterBalanceData()
        CLM.waterbalance_init!(wb, 5, 8, 2)

        @test CLM.waterbalance_init_history!(wb, 1:5) === nothing
    end

    @testset "field mutability" begin
        wb = CLM.WaterBalanceData()
        CLM.waterbalance_init!(wb, 3, 4, 2)

        wb.h2osno_old_col[1] = 42.0
        @test wb.h2osno_old_col[1] == 42.0

        wb.snow_sources_col[2] = 1.5
        @test wb.snow_sources_col[2] == 1.5

        wb.snow_sinks_col[3] = 0.3
        @test wb.snow_sinks_col[3] == 0.3

        wb.begwb_col[1] = 100.0
        @test wb.begwb_col[1] == 100.0

        wb.endwb_col[1] = 200.0
        @test wb.endwb_col[1] == 200.0

        wb.errh2o_col[2] = -0.01
        @test wb.errh2o_col[2] == -0.01

        wb.errh2osno_col[1] = 0.05
        @test wb.errh2osno_col[1] == 0.05

        wb.errh2o_patch[3] = -0.5
        @test wb.errh2o_patch[3] == -0.5

        wb.liq1_grc[1] = 500.0
        @test wb.liq1_grc[1] == 500.0

        wb.begwb_grc[2] = 1000.0
        @test wb.begwb_grc[2] == 1000.0
    end

    @testset "re-init overwrites previous state" begin
        wb = CLM.WaterBalanceData()
        CLM.waterbalance_init!(wb, 3, 4, 2)
        wb.h2osno_old_col[1] = 999.0

        CLM.waterbalance_init!(wb, 7, 10, 5)
        @test length(wb.h2osno_old_col) == 7
        @test all(isnan, wb.h2osno_old_col)
        @test length(wb.errh2o_patch) == 10
        @test length(wb.liq1_grc) == 5
        @test all(wb.wa_reset_nonconservation_gain_col .== 0.0)
    end
end
