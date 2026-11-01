# ==========================================================================
# KernelAbstractions 0.10 readiness (src/infrastructure/kernels.jl).
#
# These run on the KA 0.9 CPU backend today and pin down the pieces CLM added
# so that the move to KA 0.10 is a compat bump (docs/GPU_AMDGPU.md):
#   * `_uses_launch_rule` replaces `backend isa KA.GPU` (0.10 aliases GPU to Backend);
#   * `DeviceBitArray`, the device form of a BitArray mask, which the KA 0.10
#     CPU (PoCL) backend needs because a BitVector is not a bitstype.
# ==========================================================================

using Test
using CLM
import KernelAbstractions as KA
using KernelAbstractions: @kernel, @index, @Const
import Adapt

@kernel function _ka010_copy_mask!(out, @Const(m))
    i = @index(Global)
    @inbounds out[i] = m[i] ? 2.0 : -1.0
end

@kernel function _ka010_write_mask!(m, @Const(x))
    i = @index(Global)
    @inbounds m[i] = x[i] > 0
end

struct _Ka010FakeBackend <: KA.Backend end

@testset "KA 0.10 readiness" begin
    @testset "launch-rule trait (no dispatch on KA.GPU)" begin
        @test CLM._uses_launch_rule(KA.CPU()) === false
        # A non-CPU backend takes the launch rule, whatever KA's GPU type hierarchy is.
        @test CLM._uses_launch_rule(_Ka010FakeBackend()) === true
        # KA.CPU is the only concrete 0.9 Backend that is not a GPU, so on 0.9 the trait
        # is exactly the old `backend isa KA.GPU` test for the CPU backend.
        if isdefined(KA, :GPU) && KA.GPU !== KA.Backend
            @test CLM._uses_launch_rule(KA.CPU()) == (KA.CPU() isa KA.GPU)
        end
    end

    @testset "DeviceBitArray matches BitArray bit for bit" begin
        for n in (1, 63, 64, 65, 200)
            b = rand(Bool, n) |> BitVector
            d = CLM.DeviceBitArray(b)
            @test size(d) == size(b) && collect(d) == collect(b)
            @test d.chunks === b.chunks                         # shares storage
            @test_throws BoundsError d[n + 1]
            # writes go to the BitVector's own chunks
            for i in 1:n
                d[i] = !b[i]
            end
            @test collect(d) == collect(b)
            ref = copy(b)
            d[1] = true;  ref[1] = true
            d[n] = false; ref[n] = false
            @test b == ref
        end
        B = BitArray(rand(Bool, 5, 13))
        D = CLM.DeviceBitArray(B)
        @test size(D) == (5, 13) && all(D[i, j] == B[i, j] for i in 1:5, j in 1:13)
        # Adapt keeps the bit layout (adapt rule for a device backend's adaptor).
        D2 = Adapt.adapt(Array, D)
        @test D2 isa CLM.DeviceBitArray && collect(D2) == collect(B)
    end

    @testset "DeviceBitArray as a KA kernel argument" begin
        n = 130
        b = BitVector(rand(Bool, n))
        out = zeros(n)
        CLM._launch!(_ka010_copy_mask!, out, CLM.DeviceBitArray(b))
        @test out == [bi ? 2.0 : -1.0 for bi in b]
        # written mask (atomic per chunk), launched on the mask's own backend
        x = randn(n)
        m = falses(n)
        CLM._launch!(_ka010_write_mask!, CLM.DeviceBitArray(m), x; ndrange = n)
        @test m == (x .> 0)
        # the plain-BitVector launch (KA 0.9 CPU) gives the same result
        out2 = zeros(n)
        CLM._launch!(_ka010_copy_mask!, out2, b)
        @test out2 == out
    end

    # KA 0.10 (PoCL CPU backend): BitArray kernel arguments become DeviceBitArray.
    # On KA 0.9 the CPU backend takes the BitVector as is and no rule is defined.
    if isdefined(KA, :POCL) && isdefined(KA.POCL, :KernelAdaptor)
        @testset "BitArray adapt rule on the PoCL CPU backend" begin
            b = BitVector(rand(Bool, 70))
            d = KA.KernelInterface.argconvert(KA.CPU(), b)   # what a CPU launch does to it
            @test d isa CLM.DeviceBitArray
            @test isbits(d)
            out = zeros(70)
            CLM._launch!(_ka010_copy_mask!, out, b)   # plain BitVector launch on PoCL
            @test out == [bi ? 2.0 : -1.0 for bi in b]
            m = falses(70); x = randn(70)
            CLM._launch!(_ka010_write_mask!, m, x)    # written BitVector (atomic per chunk)
            @test m == (x .> 0)
        end
    end
end
