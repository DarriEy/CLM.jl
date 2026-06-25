# ==========================================================================
# test_distributed_driver.jl — single-rank (CI / serial-suite) tests for the
# MPI distributed-execution driver layer (distributed_driver.jl).
#
# CI's normal suite has NO `mpiexec`, so every routine here must reduce to its
# transparent single-rank passthrough:
#   * decompInit_distributed!  == decompInit!(iam=0, npes=1)
#   * clm_run_distributed!     == the whole-proc clump loop, nsteps times
#   * gather_to_master         == local field placed into global index order
#   * scatter_from_master      == global_field[gindex] (master keeps its block)
#
# The genuine 2-rank bit-identity assertion lives in test/mpi/run_mpi_smoke.jl
# and only runs under `mpiexec -n 2` (guarded out of this serial suite).
# ==========================================================================

using Test
using CLM
const _D = CLM

@testset "distributed driver (single-rank no-op path)" begin

    # reset SPMD to single-rank defaults so we test the passthrough path
    _D.spmd_init!()
    @test _D.is_mpi_active() == false
    @test _D.get_npes() == 1

    # ----------------------------------------------------------------------
    # decompInit_distributed! == serial decompInit! on one rank.
    # ----------------------------------------------------------------------
    @testset "decompInit_distributed! == serial decomp (npes=1)" begin
        numg = 8
        d_dist = _D.DecompData()
        _D.decompInit_distributed!(numg; clump_pproc = 1, decomp_data = d_dist)

        d_ref = _D.DecompData()
        _D.decompInit!(numg; clump_pproc = 1, npes = 1, iam = 0, decomp_data = d_ref)

        @test d_dist.numg == d_ref.numg
        @test d_dist.numc == d_ref.numc
        @test d_dist.gindex_grc == d_ref.gindex_grc
        @test d_dist.gindex_col == d_ref.gindex_col
        bp_dist = _D.get_proc_bounds(decomp_data = d_dist)
        bp_ref  = _D.get_proc_bounds(decomp_data = d_ref)
        @test (bp_dist.begc, bp_dist.endc) == (bp_ref.begc, bp_ref.endc)

        # single rank: gindex is the identity 1:numg
        @test d_dist.gindex_grc == collect(1:numg)
    end

    # ----------------------------------------------------------------------
    # rank_clump_bounds: on one rank with clump_pproc=1 this is the whole-proc
    # clump; with clump_pproc=N it tiles the proc into N disjoint clumps.
    # ----------------------------------------------------------------------
    @testset "rank_clump_bounds tiles the local proc domain" begin
        numg = 12
        d = _D.DecompData()
        _D.decompInit_distributed!(numg; clump_pproc = 4, decomp_data = d)
        bp = _D.get_proc_bounds(decomp_data = d)
        cb = _D.rank_clump_bounds(decomp_data = d)
        @test length(cb) == 4
        # disjoint + complete cover of the proc columns
        cov = Int[]
        for b in cb
            append!(cov, b.begc:b.endc)
        end
        @test sort(cov) == collect(bp.begc:bp.endc)

        # clump_pproc=1 -> single whole-proc clump
        d1 = _D.DecompData()
        _D.decompInit_distributed!(numg; clump_pproc = 1, decomp_data = d1)
        cb1 = _D.rank_clump_bounds(decomp_data = d1)
        @test length(cb1) == 1
        bp1 = _D.get_proc_bounds(decomp_data = d1)
        @test (cb1[1].begc, cb1[1].endc) == (bp1.begc, bp1.endc)
    end

    # ----------------------------------------------------------------------
    # clm_run_distributed!: runs the rank's clump-safe physics nsteps times.
    # On one rank == whole-proc loop; multi-clump == single-clump (determinism).
    # ----------------------------------------------------------------------
    @testset "clm_run_distributed! advances the local subdomain" begin
        numg = 12
        ncols = fill(1, numg)

        # a clump-safe per-column physics step: y[c] += f(x[c]) over the clump.
        build(nc) = (x = collect(Float64, 1:nc) .* 0.25, y = zeros(nc))
        phys_factory(st) = bc -> begin
            for c in bc.begc:bc.endc
                st.y[c] += sqrt(st.x[c]) + st.x[c]^2 - 0.5 * st.x[c]
            end
        end

        # reference: single whole-proc clump, run nsteps serially by hand
        d1 = _D.DecompData()
        _D.decompInit_distributed!(numg; clump_pproc = 1, ncols_per_g = ncols, decomp_data = d1)
        nc = d1.numc
        ref = build(nc)
        nsteps = 3
        for _ in 1:nsteps
            phys_factory(ref)(_D.get_proc_bounds(decomp_data = d1))
        end

        # via clm_run_distributed! (single clump)
        s1 = build(nc)
        _D.clm_run_distributed!(phys_factory(s1), nsteps; decomp_data = d1, threaded = false)
        @test s1.y == ref.y

        # via clm_run_distributed! with 4 clumps (multi-clump determinism)
        d4 = _D.DecompData()
        _D.decompInit_distributed!(numg; clump_pproc = 4, ncols_per_g = ncols, decomp_data = d4)
        s4 = build(nc)
        _D.clm_run_distributed!(phys_factory(s4), nsteps; decomp_data = d4, threaded = false)
        @test s4.y == ref.y

        # step-aware physics signature is detected and threaded through
        d1b = _D.DecompData()
        _D.decompInit_distributed!(numg; clump_pproc = 1, ncols_per_g = ncols, decomp_data = d1b)
        sstep = build(nc)
        function phys_step!(bc::_D.BoundsType, step::Int)
            for c in bc.begc:bc.endc
                sstep.y[c] += step * sstep.x[c]
            end
        end
        _D.clm_run_distributed!(phys_step!, nsteps; decomp_data = d1b, threaded = false)
        # sum_{step=1..3} step*x = 6*x
        @test sstep.y ≈ 6.0 .* sstep.x
    end

    # ----------------------------------------------------------------------
    # gather_to_master / scatter_from_master: single-rank identity round-trip,
    # at gridcell + column levels, with multi-column gridcells.
    # ----------------------------------------------------------------------
    @testset "gather/scatter to master — single-rank identity" begin
        numg = 6
        ncols = [1, 2, 1, 3, 2, 1]   # variable columns per gridcell
        d = _D.DecompData()
        _D.decompInit_distributed!(numg; clump_pproc = 1,
                                   ncols_per_g = ncols, npatches_per_g = ncols,
                                   decomp_data = d)

        # gridcell-level field
        gfield = Float64.(1:numg) .* 10.0
        gathered_g = _D.gather_to_master(gfield, _D.SUBGRID_LEVEL_GRIDCELL; decomp_data = d)
        @test length(gathered_g) == d.numg
        @test gathered_g == gfield               # identity gindex on one rank

        scattered_g = _D.scatter_from_master(gathered_g, _D.SUBGRID_LEVEL_GRIDCELL; decomp_data = d)
        @test scattered_g == gfield              # exact round-trip

        # column-level field (numc entries)
        cfield = Float64.(1:d.numc) .* 0.5 .+ 100.0
        gathered_c = _D.gather_to_master(cfield, _D.SUBGRID_LEVEL_COLUMN; decomp_data = d)
        @test length(gathered_c) == d.numc
        @test gathered_c == cfield
        scattered_c = _D.scatter_from_master(gathered_c, _D.SUBGRID_LEVEL_COLUMN; decomp_data = d)
        @test scattered_c == cfield
    end

    # ----------------------------------------------------------------------
    # gather reassembly uses global index order: even when the local field is
    # in a non-identity gindex order, gather places each value at gindex[k].
    # (Single rank with identity gindex still exercises the out[gindex[k]] path.)
    # ----------------------------------------------------------------------
    @testset "gather reassembly is by global index" begin
        numg = 5
        d = _D.DecompData()
        _D.decompInit_distributed!(numg; clump_pproc = 1, decomp_data = d)
        gindex = _D.get_subgrid_level_gindex(_D.SUBGRID_LEVEL_GRIDCELL; decomp_data = d)
        local_field = Float64.(101:105)
        out = _D.gather_to_master(local_field, _D.SUBGRID_LEVEL_GRIDCELL; decomp_data = d)
        for k in eachindex(local_field)
            @test out[gindex[k]] == local_field[k]
        end
    end

    # length-mismatch guards
    @testset "gather length guard" begin
        numg = 4
        d = _D.DecompData()
        _D.decompInit_distributed!(numg; clump_pproc = 1, decomp_data = d)
        @test_throws ErrorException _D.gather_to_master(Float64[1, 2, 3], _D.SUBGRID_LEVEL_GRIDCELL; decomp_data = d)
    end

    _D.spmd_finalize!()
end
