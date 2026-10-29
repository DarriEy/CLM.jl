# ==========================================================================
# Launch-level scatter primitives (src/infrastructure/kernels.jl).
#
# `scatter_add_1d!` launches over the WHOLE `values` array. Subgrid index
# arrays (`patch%column`, `col%landunit`, ...) are allocated ISPVAL (-9999) and
# only filled in for the patches a run actually uses, so a destination index of
# -9999 reaches the kernel routinely — `patch_init!` leaves it that way, and so
# does every hand-built fixture with a tail of unused patches.
#
# Scattering into that index writes OUTSIDE `out`. `Atomix.@atomic` performs its
# own `checkbounds`, which the `@inbounds` at the call site does not elide, so
# it throws under `--check-bounds=yes` (what Pkg.test and CI use) — but with
# bounds checking off it is a SILENT out-of-bounds write past the end of `out`.
# Hence the destination validation, and hence these tests.
# ==========================================================================

using Test
using CLM

@testset "scatter_add_1d! destination validation" begin

    @testset "ordinary scatter still sums per destination" begin
        out    = zeros(3)
        values = [1.0, 2.0, 3.0, 4.0]
        dest   = [1, 2, 1, 3]
        CLM.scatter_add_1d!(out, values, dest)
        @test out ≈ [1.0 + 3.0, 2.0, 4.0]
    end

    @testset "out is reset, not accumulated across calls" begin
        out = fill(99.0, 2)
        CLM.scatter_add_1d!(out, [1.0, 2.0], [1, 2])
        @test out ≈ [1.0, 2.0]
    end

    @testset "ISPVAL destinations are skipped, not written out of bounds" begin
        # The real shape of the bug: patches 3 and 4 belong to no column, which
        # is what `patch_init!` leaves behind for the unused tail of the array.
        out    = zeros(2)
        values = [1.0, 2.0, 5.0, 7.0]
        dest   = [1, 2, CLM.ISPVAL, CLM.ISPVAL]

        @test CLM.scatter_add_1d!(out, values, dest) === nothing   # no BoundsError
        # a patch belonging to no column contributes nothing
        @test out ≈ [1.0, 2.0]
    end

    @testset "destinations past the end of `out` are skipped too" begin
        out  = zeros(2)
        @test CLM.scatter_add_1d!(out, [1.0, 2.0, 3.0], [1, 2, 17]) === nothing
        @test out ≈ [1.0, 2.0]
        # ...and so is 0 (Fortran 1-based indexing has no element 0)
        fill!(out, 0.0)
        CLM.scatter_add_1d!(out, [1.0, 2.0, 3.0], [1, 2, 0])
        @test out ≈ [1.0, 2.0]
    end

    @testset "the reverse-mode gather skips the same destinations" begin
        # _scatter_add_1d_pullback! reads dout[dest[p]]; an ISPVAL destination
        # would read out of bounds exactly like the forward scatter writes.
        dvalues = zeros(4)
        dout    = [10.0, 20.0]
        dest    = [1, 2, CLM.ISPVAL, CLM.ISPVAL]
        @test CLM._scatter_add_1d_pullback!(dvalues, dout, dest) === nothing
        @test dvalues ≈ [10.0, 20.0, 0.0, 0.0]
    end
end
