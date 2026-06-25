# test_spmd.jl — single-rank (CI) tests for the SPMD/MPI communication shim.
#
# CI has no `mpiexec`, so every primitive must reduce to its trivial local
# operation: bcast=identity, allreduce=identity (== local reduction over one
# rank), gatherv/scatterv=copy, barrier=no-op.

@testset "spmdMod (SPMD/MPI shim)" begin

    @testset "single-rank identity / query" begin
        CLM.spmd_init!()
        @test CLM.is_mpi_active() == false
        @test CLM.get_npes() == 1
        @test CLM.get_rank() == 0
        @test CLM.is_master() == true
    end

    @testset "bcast = identity" begin
        a = [1.0, 2.0, 3.0]
        out = CLM.spmd_bcast!(a)
        @test out === a              # in-place, same array
        @test out == [1.0, 2.0, 3.0] # unchanged on single rank
        # scalar form returns the value
        @test CLM.spmd_bcast!(42.0) == 42.0
        @test CLM.spmd_bcast!(7) == 7
    end

    @testset "allreduce scalar = identity on one rank" begin
        @test CLM.spmd_allreduce(5.0, :SUM) == 5.0
        @test CLM.spmd_allreduce(5.0, :MIN) == 5.0
        @test CLM.spmd_allreduce(5.0, :MAX) == 5.0
        @test CLM.spmd_allreduce(3) == 3            # default op :SUM
        @test_throws ArgumentError CLM.spmd_allreduce(1.0, :PROD)
    end

    @testset "allreduce array = element-wise identity on one rank" begin
        v = [1.0, -2.0, 3.0]
        out = CLM.spmd_allreduce(v, :SUM)
        @test out == v
        @test out !== v   # returns a copy, not the same array
    end

    @testset "allreduce_local == local reduction (SUM/MIN/MAX)" begin
        arr = [3.0, 1.0, 4.0, 1.0, 5.0, 9.0, 2.0, 6.0]
        @test CLM.spmd_allreduce_local(arr, :SUM) == sum(arr)
        @test CLM.spmd_allreduce_local(arr, :MIN) == minimum(arr)
        @test CLM.spmd_allreduce_local(arr, :MAX) == maximum(arr)
        # integer array
        iarr = [7, -3, 10, 2]
        @test CLM.spmd_allreduce_local(iarr, :SUM) == sum(iarr)
        @test CLM.spmd_allreduce_local(iarr, :MIN) == minimum(iarr)
        @test CLM.spmd_allreduce_local(iarr, :MAX) == maximum(iarr)
        @test_throws ArgumentError CLM.spmd_allreduce_local(arr, :PROD)
    end

    @testset "gatherv = copy on one rank" begin
        # On a single rank, counts is just [length(local)] and gatherv copies.
        local_block = [10.0, 20.0, 30.0, 40.0]
        counts = [length(local_block)]
        g = CLM.spmd_gatherv(local_block, counts)
        @test g == local_block
        @test g !== local_block
    end

    @testset "scatterv = copy on one rank" begin
        global_vec = [10.0, 20.0, 30.0, 40.0]
        counts = [length(global_vec)]
        s = CLM.spmd_scatterv(global_vec, counts)
        @test s == global_vec
        @test s !== global_vec
    end

    @testset "gatherv -> scatterv round-trips a begg:endg vector" begin
        # Simulate a per-gridcell field on one rank: gather to master, then
        # scatter back. On a single rank this is identity round-trip.
        begg, endg = 1, 6
        field = Float64.(begg:endg) .* 1.5     # arbitrary per-gridcell values
        counts = [endg - begg + 1]             # this rank owns all gridcells

        gathered = CLM.spmd_gatherv(field, counts)
        @test length(gathered) == sum(counts)

        scattered = CLM.spmd_scatterv(gathered, counts)
        @test scattered == field               # exact round-trip
        @test eltype(scattered) == eltype(field)
    end

    @testset "barrier no-ops" begin
        @test CLM.spmd_barrier() === nothing
    end

    @testset "finalize resets to single-rank defaults" begin
        CLM.spmd_finalize!()
        @test CLM.is_mpi_active() == false
        @test CLM.get_npes() == 1
        @test CLM.get_rank() == 0
        @test CLM.is_master() == true
    end
end
