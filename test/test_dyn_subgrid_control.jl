@testset "dynSubgridControl Module" begin

    # ------------------------------------------------------------------
    # Default construction
    # ------------------------------------------------------------------
    @testset "default construction" begin
        ctl = CLM.dyn_subgrid_control_init()
        @test ctl isa CLM.DynSubgridControl
        @test ctl.initialized == true

        # all flags default to Fortran defaults
        @test CLM.get_flanduse_timeseries(ctl) == ""
        @test CLM.get_do_transient_pfts(ctl)   == false
        @test CLM.get_do_transient_crops(ctl)  == false
        @test CLM.get_do_transient_lakes(ctl)  == false
        @test CLM.get_do_transient_urban(ctl)  == false
        @test CLM.get_do_harvest(ctl)          == false
        @test CLM.get_do_grossunrep(ctl)       == false
        @test CLM.get_reset_dynbal_baselines(ctl) == false
        @test CLM.get_for_testing_allow_non_annual_changes(ctl) == false
        @test CLM.get_for_testing_zero_dynbal_fluxes(ctl)       == false
        @test CLM.run_has_transient_landcover(ctl) == false
    end

    # ------------------------------------------------------------------
    # Getters return the set values (with a valid flanduse_timeseries
    # + use_cn so the consistency checks pass)
    # ------------------------------------------------------------------
    @testset "flag getters return set values" begin
        ctl = CLM.dyn_subgrid_control_init(
            flanduse_timeseries = "landuse.timeseries.nc",
            do_transient_pfts   = true,
            do_transient_crops  = true,
            do_transient_lakes  = true,
            do_transient_urban  = true,
            do_harvest          = true,
            do_grossunrep       = true,
            reset_dynbal_baselines = true,
            for_testing_allow_non_annual_changes = true,
            for_testing_zero_dynbal_fluxes       = true,
            use_cn = true,   # needed for do_harvest + do_grossunrep
        )

        @test CLM.get_flanduse_timeseries(ctl) == "landuse.timeseries.nc"
        @test CLM.get_do_transient_pfts(ctl)   == true
        @test CLM.get_do_transient_crops(ctl)  == true
        @test CLM.get_do_transient_lakes(ctl)  == true
        @test CLM.get_do_transient_urban(ctl)  == true
        @test CLM.get_do_harvest(ctl)          == true
        @test CLM.get_do_grossunrep(ctl)       == true
        @test CLM.get_reset_dynbal_baselines(ctl) == true
        @test CLM.get_for_testing_allow_non_annual_changes(ctl) == true
        @test CLM.get_for_testing_zero_dynbal_fluxes(ctl)       == true
        @test CLM.run_has_transient_landcover(ctl) == true
    end

    # ------------------------------------------------------------------
    # run_has_transient_landcover checks pfts/crops/urban (NOT lakes)
    # ------------------------------------------------------------------
    @testset "run_has_transient_landcover excludes lakes" begin
        ctl = CLM.dyn_subgrid_control_init(
            flanduse_timeseries = "landuse.timeseries.nc",
            do_transient_lakes  = true,
        )
        @test CLM.run_has_transient_landcover(ctl) == false
    end

    # ------------------------------------------------------------------
    # Uninitialized getter assertion
    # ------------------------------------------------------------------
    @testset "uninitialized getter throws" begin
        raw = CLM.DynSubgridControl()  # initialized = false
        @test_throws AssertionError CLM.get_do_transient_pfts(raw)
    end

    # ------------------------------------------------------------------
    # Consistency: transient/harvest options need a flanduse_timeseries
    # ------------------------------------------------------------------
    @testset "blank flanduse_timeseries rejects transient options" begin
        @test_throws ErrorException CLM.dyn_subgrid_control_init(do_transient_pfts = true)
        @test_throws ErrorException CLM.dyn_subgrid_control_init(do_transient_crops = true)
        @test_throws ErrorException CLM.dyn_subgrid_control_init(do_transient_lakes = true)
        @test_throws ErrorException CLM.dyn_subgrid_control_init(do_transient_urban = true)
        @test_throws ErrorException CLM.dyn_subgrid_control_init(do_harvest = true, use_cn = true)
        @test_throws ErrorException CLM.dyn_subgrid_control_init(do_grossunrep = true, use_cn = true)
    end

    # ------------------------------------------------------------------
    # Consistency: do_transient_pfts incompatible with use_cndv / use_fates
    # ------------------------------------------------------------------
    @testset "do_transient_pfts incompatibilities" begin
        @test_throws ErrorException CLM.dyn_subgrid_control_init(
            flanduse_timeseries = "f.nc", do_transient_pfts = true, use_cndv = true)
        @test_throws ErrorException CLM.dyn_subgrid_control_init(
            flanduse_timeseries = "f.nc", do_transient_pfts = true, use_fates = true)
    end

    # ------------------------------------------------------------------
    # Consistency: transient areas incompatible with collapse_urban
    # ------------------------------------------------------------------
    @testset "transient incompatible with collapse_urban" begin
        @test_throws ErrorException CLM.dyn_subgrid_control_init(
            flanduse_timeseries = "f.nc", do_transient_lakes = true, collapse_urban = true)
    end

    # ------------------------------------------------------------------
    # Consistency: transient areas incompatible with n_dom_* / toosmall_*
    # ------------------------------------------------------------------
    @testset "transient incompatible with n_dom / toosmall thresholds" begin
        @test_throws ErrorException CLM.dyn_subgrid_control_init(
            flanduse_timeseries = "f.nc", do_transient_pfts = true, n_dom_pfts = 2)
        @test_throws ErrorException CLM.dyn_subgrid_control_init(
            flanduse_timeseries = "f.nc", do_transient_pfts = true, n_dom_landunits = 1)
        @test_throws ErrorException CLM.dyn_subgrid_control_init(
            flanduse_timeseries = "f.nc", do_transient_pfts = true, toosmall_soil = 0.1)
        @test_throws ErrorException CLM.dyn_subgrid_control_init(
            flanduse_timeseries = "f.nc", do_transient_urban = true, toosmall_urban = 0.5)
    end

    # ------------------------------------------------------------------
    # Consistency: do_transient_crops incompatible with use_fates
    # ------------------------------------------------------------------
    @testset "do_transient_crops incompatible with use_fates" begin
        @test_throws ErrorException CLM.dyn_subgrid_control_init(
            flanduse_timeseries = "f.nc", do_transient_crops = true, use_fates = true)
    end

    # ------------------------------------------------------------------
    # Consistency: do_harvest needs use_cn OR use_fates
    # ------------------------------------------------------------------
    @testset "do_harvest needs use_cn or use_fates" begin
        @test_throws ErrorException CLM.dyn_subgrid_control_init(
            flanduse_timeseries = "f.nc", do_harvest = true)
        # passes with use_cn
        ctl1 = CLM.dyn_subgrid_control_init(
            flanduse_timeseries = "f.nc", do_harvest = true, use_cn = true)
        @test CLM.get_do_harvest(ctl1) == true
        # passes with use_fates
        ctl2 = CLM.dyn_subgrid_control_init(
            flanduse_timeseries = "f.nc", do_harvest = true, use_fates = true)
        @test CLM.get_do_harvest(ctl2) == true
    end

    # ------------------------------------------------------------------
    # Consistency: do_grossunrep needs use_cn and is incompatible with use_fates
    # ------------------------------------------------------------------
    @testset "do_grossunrep needs use_cn, not use_fates" begin
        # not use_cn → error
        @test_throws ErrorException CLM.dyn_subgrid_control_init(
            flanduse_timeseries = "f.nc", do_grossunrep = true)
        # use_fates → error (even if use_cn is also true)
        @test_throws ErrorException CLM.dyn_subgrid_control_init(
            flanduse_timeseries = "f.nc", do_grossunrep = true, use_cn = true, use_fates = true)
        # use_cn alone → ok
        ctl = CLM.dyn_subgrid_control_init(
            flanduse_timeseries = "f.nc", do_grossunrep = true, use_cn = true)
        @test CLM.get_do_grossunrep(ctl) == true
    end

    # ------------------------------------------------------------------
    # masterproc = false skips the consistency check (matches Fortran)
    # ------------------------------------------------------------------
    @testset "masterproc=false skips consistency check" begin
        # would normally error (blank flanduse + do_transient_pfts), but skipped
        ctl = CLM.dyn_subgrid_control_init(do_transient_pfts = true, masterproc = false)
        @test ctl.initialized == true
        @test CLM.get_do_transient_pfts(ctl) == true
    end

end
