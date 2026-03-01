@testset "DecompMod" begin

    @testset "subgrid level constants" begin
        @test CLM.SUBGRID_LEVEL_UNSPECIFIED == -1
        @test CLM.SUBGRID_LEVEL_LNDGRID    == 0
        @test CLM.SUBGRID_LEVEL_GRIDCELL   == 1
        @test CLM.SUBGRID_LEVEL_LANDUNIT   == 2
        @test CLM.SUBGRID_LEVEL_COLUMN     == 3
        @test CLM.SUBGRID_LEVEL_PATCH      == 4
        @test CLM.SUBGRID_LEVEL_COHORT     == 5
    end

    @testset "bounds level constants" begin
        @test CLM.BOUNDS_LEVEL_PROC  == 1
        @test CLM.BOUNDS_LEVEL_CLUMP == 2
    end

    @testset "BoundsType default construction" begin
        b = CLM.BoundsType()
        @test b.begg == 0
        @test b.endg == 0
        @test b.begl == 0
        @test b.endl == 0
        @test b.begc == 0
        @test b.endc == 0
        @test b.begp == 0
        @test b.endp == 0
        @test b.begCohort == 0
        @test b.endCohort == 0
        @test b.level == 0
        @test b.clump_index == -1
    end

    @testset "ProcessorType default construction" begin
        p = CLM.ProcessorType()
        @test p.nclumps == 0
        @test isempty(p.cid)
        @test p.ncells == 0
        @test p.npatches == 0
        @test p.begg == 0
        @test p.endg == 0
    end

    @testset "ClumpType default construction" begin
        c = CLM.ClumpType()
        @test c.owner == 0
        @test c.ncells == 0
        @test c.begg == 0
        @test c.endg == 0
    end

    @testset "get_beg and get_end" begin
        b = CLM.BoundsType(begg=1, endg=10, begl=1, endl=5,
                           begc=1, endc=20, begp=1, endp=50,
                           begCohort=1, endCohort=3)
        # get_beg
        @test CLM.get_beg(b, CLM.SUBGRID_LEVEL_GRIDCELL) == 1
        @test CLM.get_beg(b, CLM.SUBGRID_LEVEL_LANDUNIT) == 1
        @test CLM.get_beg(b, CLM.SUBGRID_LEVEL_COLUMN)   == 1
        @test CLM.get_beg(b, CLM.SUBGRID_LEVEL_PATCH)    == 1
        @test CLM.get_beg(b, CLM.SUBGRID_LEVEL_COHORT)   == 1
        @test CLM.get_beg(b, -99) == -1  # invalid level

        # get_end
        @test CLM.get_end(b, CLM.SUBGRID_LEVEL_GRIDCELL) == 10
        @test CLM.get_end(b, CLM.SUBGRID_LEVEL_LANDUNIT) == 5
        @test CLM.get_end(b, CLM.SUBGRID_LEVEL_COLUMN)   == 20
        @test CLM.get_end(b, CLM.SUBGRID_LEVEL_PATCH)    == 50
        @test CLM.get_end(b, CLM.SUBGRID_LEVEL_COHORT)   == 3
        @test CLM.get_end(b, -99) == -1  # invalid level
    end

    # Helper: set up a DecompData with a simple 1-processor, 2-clump layout
    function make_test_decomp()
        dd = CLM.DecompData()

        # Processor owns 10 gridcells, 3 landunits, 15 columns, 40 patches, 2 cohorts
        dd.procinfo = CLM.ProcessorType(
            nclumps=2, cid=[1, 2],
            ncells=10, nlunits=3, ncols=15, npatches=40, nCohorts=2,
            begg=1, endg=10, begl=1, endl=3, begc=1, endc=15,
            begp=1, endp=40, begCohort=1, endCohort=2
        )

        # Clump 1: gridcells 1-6, landunits 1-2, columns 1-9, patches 1-25, cohorts 1-1
        c1 = CLM.ClumpType(
            owner=1, ncells=6, nlunits=2, ncols=9, npatches=25, nCohorts=1,
            begg=1, endg=6, begl=1, endl=2, begc=1, endc=9,
            begp=1, endp=25, begCohort=1, endCohort=1
        )

        # Clump 2: gridcells 7-10, landunits 3-3, columns 10-15, patches 26-40, cohorts 2-2
        c2 = CLM.ClumpType(
            owner=1, ncells=4, nlunits=1, ncols=6, npatches=15, nCohorts=1,
            begg=7, endg=10, begl=3, endl=3, begc=10, endc=15,
            begp=26, endp=40, begCohort=2, endCohort=2
        )

        dd.clumps = [c1, c2]
        dd.nclumps = 2

        # Global sizes
        dd.numg = 10
        dd.numl = 3
        dd.numc = 15
        dd.nump = 40
        dd.numCohort = 2
        dd.ldomain_ns = 12  # includes ocean

        # Global index arrays (identity mapping for simplicity)
        dd.gindex_global = collect(1:12)
        dd.gindex_grc   = collect(1:10)
        dd.gindex_lun   = collect(1:3)
        dd.gindex_col   = collect(1:15)
        dd.gindex_patch  = collect(1:40)
        dd.gindex_cohort = collect(1:2)

        return dd
    end

    @testset "get_proc_bounds" begin
        dd = make_test_decomp()
        b = CLM.get_proc_bounds(decomp_data=dd)

        @test b.begg == 1
        @test b.endg == 10
        @test b.begl == 1
        @test b.endl == 3
        @test b.begc == 1
        @test b.endc == 15
        @test b.begp == 1
        @test b.endp == 40
        @test b.begCohort == 1
        @test b.endCohort == 2
        @test b.level == CLM.BOUNDS_LEVEL_PROC
        @test b.clump_index == -1
    end

    @testset "get_clump_bounds" begin
        dd = make_test_decomp()

        # Clump 1 bounds (relative to processor beginning)
        b1 = CLM.get_clump_bounds(1, decomp_data=dd)
        @test b1.begg == 1
        @test b1.endg == 6
        @test b1.begl == 1
        @test b1.endl == 2
        @test b1.begc == 1
        @test b1.endc == 9
        @test b1.begp == 1
        @test b1.endp == 25
        @test b1.begCohort == 1
        @test b1.endCohort == 1
        @test b1.level == CLM.BOUNDS_LEVEL_CLUMP
        @test b1.clump_index == 1

        # Clump 2 bounds (relative to processor beginning)
        b2 = CLM.get_clump_bounds(2, decomp_data=dd)
        @test b2.begg == 7
        @test b2.endg == 10
        @test b2.begl == 3
        @test b2.endl == 3
        @test b2.begc == 10
        @test b2.endc == 15
        @test b2.begp == 26
        @test b2.endp == 40
        @test b2.begCohort == 2
        @test b2.endCohort == 2
        @test b2.level == CLM.BOUNDS_LEVEL_CLUMP
        @test b2.clump_index == 2
    end

    @testset "get_proc_total" begin
        dd = make_test_decomp()
        totals = CLM.get_proc_total(1, decomp_data=dd)
        @test totals.ncells == 10
        @test totals.nlunits == 3
        @test totals.ncols == 15
        @test totals.npatches == 40
        @test totals.nCohorts == 2

        # Non-existent pid should return zeros
        totals2 = CLM.get_proc_total(999, decomp_data=dd)
        @test totals2.ncells == 0
        @test totals2.nlunits == 0
        @test totals2.ncols == 0
        @test totals2.npatches == 0
        @test totals2.nCohorts == 0
    end

    @testset "get_proc_global" begin
        dd = make_test_decomp()
        g = CLM.get_proc_global(decomp_data=dd)
        @test g.ng == 10
        @test g.nl == 3
        @test g.nc == 15
        @test g.np == 40
        @test g.nCohorts == 2
    end

    @testset "get_proc_clumps" begin
        dd = make_test_decomp()
        @test CLM.get_proc_clumps(decomp_data=dd) == 2
    end

    @testset "get_global_index" begin
        dd = make_test_decomp()
        # With identity gindex mapping, global index should equal local index
        @test CLM.get_global_index(1, CLM.SUBGRID_LEVEL_GRIDCELL, decomp_data=dd) == 1
        @test CLM.get_global_index(5, CLM.SUBGRID_LEVEL_GRIDCELL, decomp_data=dd) == 5
        @test CLM.get_global_index(10, CLM.SUBGRID_LEVEL_GRIDCELL, decomp_data=dd) == 10
        @test CLM.get_global_index(1, CLM.SUBGRID_LEVEL_COLUMN, decomp_data=dd) == 1
        @test CLM.get_global_index(15, CLM.SUBGRID_LEVEL_COLUMN, decomp_data=dd) == 15

        # Invalid subgrid level should throw
        @test_throws ErrorException CLM.get_global_index(1, -99, decomp_data=dd)
    end

    @testset "get_global_index_array" begin
        dd = make_test_decomp()
        idx = [1, 3, 5, 7]
        result = CLM.get_global_index_array(idx, CLM.SUBGRID_LEVEL_GRIDCELL, decomp_data=dd)
        @test result == [1, 3, 5, 7]

        # With non-identity mapping
        dd.gindex_grc = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
        result2 = CLM.get_global_index_array([1, 5, 10], CLM.SUBGRID_LEVEL_GRIDCELL, decomp_data=dd)
        @test result2 == [10, 50, 100]

        # Invalid subgrid level should throw
        @test_throws ErrorException CLM.get_global_index_array([1], -99, decomp_data=dd)
    end

    @testset "get_subgrid_level_from_name" begin
        @test CLM.get_subgrid_level_from_name(CLM.GRLND)      == CLM.SUBGRID_LEVEL_LNDGRID
        @test CLM.get_subgrid_level_from_name(CLM.NAMEG)      == CLM.SUBGRID_LEVEL_GRIDCELL
        @test CLM.get_subgrid_level_from_name(CLM.NAMEL)      == CLM.SUBGRID_LEVEL_LANDUNIT
        @test CLM.get_subgrid_level_from_name(CLM.NAMEC)      == CLM.SUBGRID_LEVEL_COLUMN
        @test CLM.get_subgrid_level_from_name(CLM.NAMEP)       == CLM.SUBGRID_LEVEL_PATCH
        @test CLM.get_subgrid_level_from_name(CLM.NAMECOHORT) == CLM.SUBGRID_LEVEL_COHORT

        # Unknown name should throw
        @test_throws ErrorException CLM.get_subgrid_level_from_name("unknown")
    end

    @testset "get_subgrid_level_gsize" begin
        dd = make_test_decomp()
        @test CLM.get_subgrid_level_gsize(CLM.SUBGRID_LEVEL_LNDGRID,  decomp_data=dd) == 12
        @test CLM.get_subgrid_level_gsize(CLM.SUBGRID_LEVEL_GRIDCELL, decomp_data=dd) == 10
        @test CLM.get_subgrid_level_gsize(CLM.SUBGRID_LEVEL_LANDUNIT, decomp_data=dd) == 3
        @test CLM.get_subgrid_level_gsize(CLM.SUBGRID_LEVEL_COLUMN,   decomp_data=dd) == 15
        @test CLM.get_subgrid_level_gsize(CLM.SUBGRID_LEVEL_PATCH,    decomp_data=dd) == 40
        @test CLM.get_subgrid_level_gsize(CLM.SUBGRID_LEVEL_COHORT,   decomp_data=dd) == 2

        # Invalid level should throw
        @test_throws ErrorException CLM.get_subgrid_level_gsize(-99, decomp_data=dd)
    end

    @testset "get_subgrid_level_gindex" begin
        dd = make_test_decomp()
        @test CLM.get_subgrid_level_gindex(CLM.SUBGRID_LEVEL_LNDGRID,  decomp_data=dd) === dd.gindex_global
        @test CLM.get_subgrid_level_gindex(CLM.SUBGRID_LEVEL_GRIDCELL, decomp_data=dd) === dd.gindex_grc
        @test CLM.get_subgrid_level_gindex(CLM.SUBGRID_LEVEL_LANDUNIT, decomp_data=dd) === dd.gindex_lun
        @test CLM.get_subgrid_level_gindex(CLM.SUBGRID_LEVEL_COLUMN,   decomp_data=dd) === dd.gindex_col
        @test CLM.get_subgrid_level_gindex(CLM.SUBGRID_LEVEL_PATCH,    decomp_data=dd) === dd.gindex_patch
        @test CLM.get_subgrid_level_gindex(CLM.SUBGRID_LEVEL_COHORT,   decomp_data=dd) === dd.gindex_cohort

        # Invalid level should throw
        @test_throws ErrorException CLM.get_subgrid_level_gindex(-99, decomp_data=dd)
    end

    @testset "DecompData default construction" begin
        dd = CLM.DecompData()
        @test dd.nclumps == 0
        @test dd.numg == 0
        @test dd.numl == 0
        @test dd.numc == 0
        @test dd.nump == 0
        @test dd.numCohort == 0
        @test isempty(dd.clumps)
        @test isempty(dd.gindex_global)
        @test isempty(dd.gindex_grc)
        @test isempty(dd.gindex_lun)
        @test isempty(dd.gindex_col)
        @test isempty(dd.gindex_patch)
        @test isempty(dd.gindex_cohort)
        @test dd.ldomain_ns == 0
    end

end
