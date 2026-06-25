@testset "decompInit!" begin

    # Helper: full subgrid-bounds round-trip check for a freshly-built decomp.
    # Verifies that the per-clump bounds (within a process) partition that
    # process's subgrid index space with no overlap/gap, and that counts are
    # consistent.

    @testset "numg=12, clump_pproc=1, npes=1 -> single clump" begin
        d = CLM.DecompData()
        CLM.decompInit!(12; clump_pproc=1, npes=1, decomp_data=d)

        @test d.nclumps == 1
        @test d.numg == 12
        @test d.procinfo.nclumps == 1
        @test d.procinfo.cid == [1]
        @test d.clumps[1].owner == 0
        @test d.clumps[1].ncells == 12
        @test d.clumps[1].begg == 1
        @test d.clumps[1].endg == 12
        # proc covers all 12 gridcells
        @test d.procinfo.begg == 1
        @test d.procinfo.endg == 12
        @test d.procinfo.ncells == 12
        # default 1 landunit/column/patch per gridcell
        @test d.numl == 12
        @test d.numc == 12
        @test d.nump == 12
        @test d.numCohort == 0
        @test d.procinfo.endl == 12
        @test d.procinfo.endc == 12
        @test d.procinfo.endp == 12

        # gindex maps: with a single clump and compressed land grid, global
        # index == local index (1:12) for every subgrid level.
        @test d.gindex_global == collect(1:12)
        @test d.gindex_grc == collect(1:12)
        @test d.gindex_lun == collect(1:12)
        @test d.gindex_col == collect(1:12)
        @test d.gindex_patch == collect(1:12)
        @test isempty(d.gindex_cohort)

        # Single-clump case must reproduce get_proc_bounds / get_clump_bounds
        # behavior exactly.
        bp = CLM.get_proc_bounds(decomp_data=d)
        @test bp.begg == 1 && bp.endg == 12
        @test bp.begl == 1 && bp.endl == 12
        @test bp.begc == 1 && bp.endc == 12
        @test bp.begp == 1 && bp.endp == 12
        @test bp.level == CLM.BOUNDS_LEVEL_PROC

        bc = CLM.get_clump_bounds(1, decomp_data=d)
        @test bc.begg == 1 && bc.endg == 12
        @test bc.begl == 1 && bc.endl == 12
        @test bc.begc == 1 && bc.endc == 12
        @test bc.begp == 1 && bc.endp == 12
        @test bc.level == CLM.BOUNDS_LEVEL_CLUMP
        @test bc.clump_index == 1

        # In single-proc-single-clump, clump bounds == proc bounds.
        @test bc.begg == bp.begg && bc.endg == bp.endg
        @test bc.begc == bp.begc && bc.endc == bp.endc
    end

    @testset "numg=12, clump_pproc=4, npes=1 -> 4 clumps, round-robin" begin
        d = CLM.DecompData()
        CLM.decompInit!(12; clump_pproc=4, npes=1, decomp_data=d)

        @test d.nclumps == 4
        @test d.numg == 12
        @test d.procinfo.nclumps == 4
        @test d.procinfo.cid == [1, 2, 3, 4]
        # all clumps owned by the single process 0
        @test all(d.clumps[c].owner == 0 for c in 1:4)

        # gridcells distributed round-robin: 12 gridcells / 4 clumps = 3 each.
        @test [d.clumps[c].ncells for c in 1:4] == [3, 3, 3, 3]
        @test sum(d.clumps[c].ncells for c in 1:4) == 12

        # per-clump gridcell bounds partition 1:12 with no overlap/gap.
        # Owner order is identical (all proc 0), so clumps are laid out in cid
        # order: clump 1 = 1:3, clump 2 = 4:6, ...
        prev_end = 0
        for c in 1:4
            @test d.clumps[c].begg == prev_end + 1
            @test d.clumps[c].endg == d.clumps[c].begg + d.clumps[c].ncells - 1
            prev_end = d.clumps[c].endg
        end
        @test prev_end == 12

        # proc covers everything.
        @test d.procinfo.begg == 1 && d.procinfo.endg == 12 && d.procinfo.ncells == 12

        # gindex_global lists this proc's gridcells in task order.
        # Task order = clumps in (pid,cid) order; gridcell n assigned to
        # clump mod(n-1,4)+1. So clump1={1,5,9}, clump2={2,6,10}, etc.
        @test d.gindex_global == [1, 5, 9, 2, 6, 10, 3, 7, 11, 4, 8, 12]

        # gindex_grc is the land-only global gridcell index of each local slot;
        # equals gindex_global on a compressed land grid.
        @test d.gindex_grc == d.gindex_global

        # Per-clump bounds partition the patch space too (1 patch/gridcell).
        prev_end = 0
        for c in 1:4
            @test d.clumps[c].begp == prev_end + 1
            @test d.clumps[c].endp == d.clumps[c].begp + d.clumps[c].npatches - 1
            prev_end = d.clumps[c].endp
        end
        @test prev_end == 12
        @test d.nump == 12

        # gindex round-trip: every local gridcell maps to a unique global index
        # in 1:numg, and the inverse recovers the local slot.
        @test sort(d.gindex_grc) == collect(1:12)

        # get_clump_bounds returns process-relative bounds; for npes=1 these are
        # the global clump bounds (proc begg=1).
        for c in 1:4
            bc = CLM.get_clump_bounds(c, decomp_data=d)
            @test bc.begg == d.clumps[c].begg
            @test bc.endg == d.clumps[c].endg
            @test bc.endg - bc.begg + 1 == 3
        end
    end

    @testset "numg=12, clump_pproc=2, npes=3 -> 6 clumps across 3 procs" begin
        # Build the GLOBAL view (iam=0) for clump ownership / counts, then build
        # each process's local procinfo by re-running with that iam.
        d0 = CLM.DecompData()
        CLM.decompInit!(12; clump_pproc=2, npes=3, iam=0, decomp_data=d0)

        @test d0.nclumps == 6
        @test d0.numg == 12

        # clump -> proc round robin: owner = mod(cid-1, 3)
        # cid:   1 2 3 4 5 6
        # owner: 0 1 2 0 1 2
        @test [d0.clumps[c].owner for c in 1:6] == [0, 1, 2, 0, 1, 2]

        # 12 gridcells round-robin over 6 clumps -> 2 each.
        @test [d0.clumps[c].ncells for c in 1:6] == [2, 2, 2, 2, 2, 2]
        @test sum(d0.clumps[c].ncells for c in 1:6) == 12

        # Global gridcell layout is ordered by (owner, cid): proc0 owns clumps
        # 1 & 4 (gridcells first), proc1 owns 2 & 5, proc2 owns 3 & 6.
        # Verify per-clump begg/endg partition 1:12 with no overlap/gap when
        # walked in (owner,cid) order.
        order = sort(1:6, by = c -> (d0.clumps[c].owner, c))
        prev_end = 0
        for c in order
            @test d0.clumps[c].begg == prev_end + 1
            @test d0.clumps[c].endg == d0.clumps[c].begg + d0.clumps[c].ncells - 1
            prev_end = d0.clumps[c].endg
        end
        @test prev_end == 12

        # Each process owns clump_pproc=2 clumps => 4 gridcells.
        for iam in 0:2
            di = CLM.DecompData()
            CLM.decompInit!(12; clump_pproc=2, npes=3, iam=iam, decomp_data=di)
            @test di.procinfo.nclumps == 2
            # this proc's clump ids: those with owner == iam
            owned = [c for c in 1:6 if di.clumps[c].owner == iam]
            @test sort(di.procinfo.cid) == owned
            @test di.procinfo.ncells == 4
            @test di.procinfo.endg - di.procinfo.begg + 1 == 4
            # local gindex_global has length = this proc's gridcell count
            @test length(di.gindex_global) == 4
            # gindex_grc round-trips to unique global indices
            @test length(unique(di.gindex_grc)) == 4
            @test all(1 .<= di.gindex_grc .<= 12)
            # subgrid counts: 1 landunit/col/patch per gridcell
            @test di.procinfo.nlunits == 4
            @test di.procinfo.ncols == 4
            @test di.procinfo.npatches == 4
        end

        # Global totals are process-independent.
        @test d0.numl == 12 && d0.numc == 12 && d0.nump == 12

        # Union of all procs' gridcells == 1:12 (complete, disjoint partition).
        allg = Int[]
        for iam in 0:2
            di = CLM.DecompData()
            CLM.decompInit!(12; clump_pproc=2, npes=3, iam=iam, decomp_data=di)
            append!(allg, di.gindex_grc)
        end
        @test sort(allg) == collect(1:12)
    end

    @testset "gindex maps round-trip local <-> global" begin
        d = CLM.DecompData()
        CLM.decompInit!(12; clump_pproc=4, npes=1, decomp_data=d)

        # get_global_index maps a local subgrid index to its global index;
        # round-trip via the gindex array recovers it.
        bp = CLM.get_proc_bounds(decomp_data=d)
        for gi_local in bp.begg:bp.endg
            g_global = CLM.get_global_index(gi_local, CLM.SUBGRID_LEVEL_GRIDCELL, decomp_data=d)
            @test g_global == d.gindex_grc[gi_local - bp.begg + 1]
            @test 1 <= g_global <= 12
        end
        # land-only global gridcell indices form a permutation of 1:12.
        globals = [CLM.get_global_index(g, CLM.SUBGRID_LEVEL_GRIDCELL, decomp_data=d)
                   for g in bp.begg:bp.endg]
        @test sort(globals) == collect(1:12)

        # array form agrees with scalar form
        arr = CLM.get_global_index_array(collect(bp.begg:bp.endg),
                                         CLM.SUBGRID_LEVEL_GRIDCELL, decomp_data=d)
        @test arr == globals

        # patch-level round-trip
        for pi_local in bp.begp:bp.endp
            p_global = CLM.get_global_index(pi_local, CLM.SUBGRID_LEVEL_PATCH, decomp_data=d)
            @test p_global == d.gindex_patch[pi_local - bp.begp + 1]
        end
    end

    @testset "non-default subgrid counts" begin
        # 4 gridcells, 1 clump; gridcell g has g landunits / g columns / 2g patches.
        nlun = [1, 2, 3, 4]
        ncol = [1, 2, 3, 4]
        npch = [2, 4, 6, 8]
        d = CLM.DecompData()
        CLM.decompInit!(4; clump_pproc=1, npes=1,
                        nlunits_per_g=nlun, ncols_per_g=ncol, npatches_per_g=npch,
                        decomp_data=d)
        @test d.numl == sum(nlun)   # 10
        @test d.numc == sum(ncol)   # 10
        @test d.nump == sum(npch)   # 20
        @test d.procinfo.endl == 10
        @test d.procinfo.endp == 20
        # gindex_lun length == total landunits, values 1:numl in gridcell order
        @test length(d.gindex_lun) == 10
        @test d.gindex_lun == collect(1:10)
        @test d.gindex_patch == collect(1:20)
    end

    @testset "input validation" begin
        d = CLM.DecompData()
        @test_throws ErrorException CLM.decompInit!(12; clump_pproc=0, decomp_data=d)
        @test_throws ErrorException CLM.decompInit!(12; npes=0, decomp_data=d)
        @test_throws ErrorException CLM.decompInit!(0, decomp_data=d)
        # more clumps than gridcells
        @test_throws ErrorException CLM.decompInit!(2; clump_pproc=4, npes=1, decomp_data=d)
        # more procs than gridcells
        @test_throws ErrorException CLM.decompInit!(2; clump_pproc=1, npes=4, decomp_data=d)
        # iam out of range
        @test_throws ErrorException CLM.decompInit!(12; npes=2, iam=5, decomp_data=d)
    end

end
