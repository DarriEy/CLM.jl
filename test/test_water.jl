@testset "WaterType — Master Water Container" begin

    @testset "WaterParams construction" begin
        p = CLM.WaterParams()
        @test p.enable_consistency_checks == false
        @test p.enable_isotopes == false

        p2 = CLM.WaterParams(enable_consistency_checks=true, enable_isotopes=true)
        @test p2.enable_consistency_checks == true
        @test p2.enable_isotopes == true
    end

    @testset "WaterInfoBulkData" begin
        bi = CLM.WaterInfoBulkData()
        @test CLM.get_name(bi) == "bulk"
        @test CLM.get_ratio(bi) ≈ 1.0
        @test CLM.is_communicated_with_coupler(bi) == true
        @test CLM.is_included_in_consistency_check(bi) == false
        @test CLM.fname(bi, "H2O") == "H2O"
        @test CLM.lname(bi, "water") == "water"
    end

    @testset "WaterInfoIsotopeData" begin
        ii = CLM.WaterInfoIsotopeData("HDO", 0.9, false, false)
        @test CLM.get_name(ii) == "HDO"
        @test CLM.get_ratio(ii) ≈ 0.9
        @test CLM.is_communicated_with_coupler(ii) == false
        @test CLM.is_included_in_consistency_check(ii) == false
        @test CLM.fname(ii, "H2O") == "H2O_HDO"
        @test CLM.lname(ii, "water") == "water HDO"

        ii2 = CLM.WaterInfoIsotopeData("H2OTR", 1.0, true, true)
        @test CLM.is_communicated_with_coupler(ii2) == true
        @test CLM.is_included_in_consistency_check(ii2) == true
    end

    @testset "Default construction" begin
        water = CLM.WaterData()
        @test water.i_bulk == 1
        @test water.bulk_tracer_index == -1
        @test water.params.enable_consistency_checks == false
        @test water.params.enable_isotopes == false
        @test isempty(water.bulk_and_tracers)
    end

    @testset "Init without tracers" begin
        water = CLM.WaterData()
        nc, np, nl, ng = 10, 20, 5, 3
        CLM.water_init!(water, nc, np, nl, ng)

        # Index structure: only bulk, no tracers
        @test water.i_bulk == 1
        @test water.bulk_and_tracers_beg == 1
        @test water.bulk_and_tracers_end == 1
        @test water.tracers_beg == 2
        @test water.tracers_end == 1   # empty range
        @test length(water.bulk_and_tracers) == 1

        # Bulk instances should be allocated with correct sizes
        @test length(water.waterfluxbulk_inst.wf.qflx_evap_tot_col) == nc
        @test length(water.waterstatebulk_inst.ws.h2osno_no_layers_col) == nc
        @test length(water.waterdiagnosticbulk_inst.h2osno_total_col) == nc
        @test length(water.waterbalancebulk_inst.begwb_col) == nc

        # Patch-level fields
        @test length(water.waterfluxbulk_inst.wf.qflx_evap_tot_patch) == np
        @test length(water.waterdiagnosticbulk_inst.fwet_patch) == np

        # bulk_and_tracers[1] should reference the same bulk instances
        @test water.bulk_and_tracers[1].waterflux === water.waterfluxbulk_inst
        @test water.bulk_and_tracers[1].waterstate === water.waterstatebulk_inst
        @test water.bulk_and_tracers[1].waterdiagnostic === water.waterdiagnosticbulk_inst
        @test water.bulk_and_tracers[1].waterbalance === water.waterbalancebulk_inst

        # Accessor: bulk name
        @test CLM.water_get_bulk_or_tracer_name(water, 1) == "bulk"

        # No consistency check
        @test CLM.water_do_consistency_check(water) == false

        # No bulk tracer
        @test CLM.water_get_bulk_tracer_index(water) == -1

        CLM.water_clean!(water)
        @test isempty(water.bulk_and_tracers)
    end

    @testset "Init with isotopes" begin
        water = CLM.WaterData()
        nc, np, nl, ng = 5, 10, 3, 2
        CLM.water_init!(water, nc, np, nl, ng; enable_isotopes=true)

        # isotopes enabled → 3 tracers (H2OTR + HDO + H218O)
        num_tracers = water.tracers_end - water.tracers_beg + 1
        @test num_tracers == 3
        @test length(water.bulk_and_tracers) == 4  # 1 bulk + 3 tracers

        # Check tracer names
        @test CLM.water_get_bulk_or_tracer_name(water, 1) == "bulk"
        @test CLM.water_get_bulk_or_tracer_name(water, 2) == "H2OTR"
        @test CLM.water_get_bulk_or_tracer_name(water, 3) == "HDO"
        @test CLM.water_get_bulk_or_tracer_name(water, 4) == "H218O"

        # All tracers are isotopes
        for i in water.tracers_beg:water.tracers_end
            @test CLM.water_is_isotope(water, i)
        end

        # Bulk tracer (H2OTR) should exist
        @test CLM.water_get_bulk_tracer_index(water) == 2

        # GetIsotopeInfo
        iso_hdo = CLM.water_get_isotope_info(water, 3)
        @test iso_hdo.tracer_name == "HDO"
        @test iso_hdo.ratio ≈ 0.9

        iso_h218o = CLM.water_get_isotope_info(water, 4)
        @test iso_h218o.tracer_name == "H218O"
        @test iso_h218o.ratio ≈ 0.5

        # Tracer sub-instances should be allocated
        for i in water.tracers_beg:water.tracers_end
            bt = water.bulk_and_tracers[i]
            @test bt.waterflux isa CLM.WaterFluxData
            @test bt.waterstate isa CLM.WaterStateData
            @test bt.waterbalance isa CLM.WaterBalanceData
            @test length(bt.waterflux.qflx_evap_tot_col) == nc
            @test length(bt.waterstate.h2osno_no_layers_col) == nc
            @test length(bt.waterbalance.begwb_col) == nc
        end

        CLM.water_clean!(water)
    end

    @testset "Init with consistency checks" begin
        water = CLM.WaterData()
        nc, np, nl, ng = 5, 10, 3, 2
        CLM.water_init!(water, nc, np, nl, ng; enable_consistency_checks=true)

        # consistency checks → 4 tracers (H2OTR + TESTMED + TESTSMALL + TESTBIG)
        num_tracers = water.tracers_end - water.tracers_beg + 1
        @test num_tracers == 4
        @test CLM.water_do_consistency_check(water) == true

        @test CLM.water_get_bulk_or_tracer_name(water, 2) == "H2OTR"
        @test CLM.water_get_bulk_or_tracer_name(water, 3) == "TESTMED"
        @test CLM.water_get_bulk_or_tracer_name(water, 4) == "TESTSMALL"
        @test CLM.water_get_bulk_or_tracer_name(water, 5) == "TESTBIG"

        # Check isotope info ratios
        iso_h2otr = CLM.water_get_isotope_info(water, 2)
        @test iso_h2otr.ratio ≈ 1.0
        @test iso_h2otr.included_in_consistency_check == true

        iso_testmed = CLM.water_get_isotope_info(water, 3)
        @test iso_testmed.ratio ≈ 0.1

        iso_testsmall = CLM.water_get_isotope_info(water, 4)
        @test iso_testsmall.ratio ≈ 1.0e-10

        iso_testbig = CLM.water_get_isotope_info(water, 5)
        @test iso_testbig.ratio ≈ 10.0

        CLM.water_clean!(water)
    end

    @testset "Init with isotopes AND consistency checks" begin
        water = CLM.WaterData()
        nc, np, nl, ng = 5, 10, 3, 2
        CLM.water_init!(water, nc, np, nl, ng;
                         enable_isotopes=true,
                         enable_consistency_checks=true)

        # Both enabled → 6 tracers (H2OTR + HDO + H218O + TESTMED + TESTSMALL + TESTBIG)
        num_tracers = water.tracers_end - water.tracers_beg + 1
        @test num_tracers == 6
        @test length(water.bulk_and_tracers) == 7  # 1 bulk + 6 tracers

        @test CLM.water_get_bulk_or_tracer_name(water, 2) == "H2OTR"
        @test CLM.water_get_bulk_or_tracer_name(water, 3) == "HDO"
        @test CLM.water_get_bulk_or_tracer_name(water, 4) == "H218O"
        @test CLM.water_get_bulk_or_tracer_name(water, 5) == "TESTMED"
        @test CLM.water_get_bulk_or_tracer_name(water, 6) == "TESTSMALL"
        @test CLM.water_get_bulk_or_tracer_name(water, 7) == "TESTBIG"

        CLM.water_clean!(water)
    end

    @testset "InitForTesting" begin
        water = CLM.WaterData()
        nc, np, nl, ng = 5, 10, 3, 2
        params = CLM.WaterParams(enable_isotopes=true)
        CLM.water_init_for_testing!(water, nc, np, nl, ng; params=params)

        @test water.params.enable_isotopes == true
        @test water.params.enable_consistency_checks == false
        num_tracers = water.tracers_end - water.tracers_beg + 1
        @test num_tracers == 3  # H2OTR + HDO + H218O

        CLM.water_clean!(water)
    end

    @testset "Stub functions do not error" begin
        water = CLM.WaterData()
        nc, np, nl, ng = 5, 10, 3, 2
        CLM.water_init!(water, nc, np, nl, ng)
        bounds_col = 1:nc

        # Accumulation stubs
        CLM.water_init_acc_buffer!(water, bounds_col)
        CLM.water_init_acc_vars!(water, bounds_col)
        CLM.water_update_acc_vars!(water, bounds_col)

        # Restart stub
        CLM.water_restart!(water, bounds_col)
        CLM.water_restart!(water, bounds_col; flag="write")

        # Tracer consistency stubs
        CLM.water_tracer_consistency_check!(water, bounds_col)
        CLM.water_tracer_consistency_check!(water, bounds_col; caller_location="test")
        CLM.water_reset_checked_tracers!(water, bounds_col)

        CLM.water_clean!(water)
    end

    @testset "Bulk-tracer aliasing" begin
        # Verify that modifying the bulk instance is visible through bulk_and_tracers
        water = CLM.WaterData()
        nc, np, nl, ng = 5, 10, 3, 2
        CLM.water_init!(water, nc, np, nl, ng)

        # Set a value via the direct bulk instance
        water.waterfluxbulk_inst.wf.qflx_evap_tot_col[1] = 42.0

        # Read it back through bulk_and_tracers
        bt_flux = water.bulk_and_tracers[water.i_bulk].waterflux
        @test bt_flux isa CLM.WaterFluxBulkData
        @test bt_flux.wf.qflx_evap_tot_col[1] ≈ 42.0

        CLM.water_clean!(water)
    end

    @testset "Iteration over bulk and tracers" begin
        water = CLM.WaterData()
        nc, np, nl, ng = 5, 10, 3, 2
        CLM.water_init!(water, nc, np, nl, ng; enable_isotopes=true)

        # Count entries by iterating
        count_all = 0
        for i in water.bulk_and_tracers_beg:water.bulk_and_tracers_end
            count_all += 1
        end
        @test count_all == 4  # 1 bulk + 3 tracers

        # Count just tracers
        count_tracers = 0
        for i in water.tracers_beg:water.tracers_end
            count_tracers += 1
        end
        @test count_tracers == 3

        # Collect names
        names = String[]
        for i in water.bulk_and_tracers_beg:water.bulk_and_tracers_end
            push!(names, CLM.water_get_bulk_or_tracer_name(water, i))
        end
        @test names == ["bulk", "H2OTR", "HDO", "H218O"]

        CLM.water_clean!(water)
    end

end
