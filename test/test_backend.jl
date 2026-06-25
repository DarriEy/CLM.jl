# ==========================================================================
# test_backend.jl — centralized backend-selection API + multi-GPU-over-MPI
# plumbing (src/infrastructure/backend.jl + multigpu.jl + ext/CLM{CUDA,AMDGPU}Ext).
#
# This machine has NO CUDA/AMD GPU, so the GPU paths are validated ONLY as far as
# the CPU backend (the correctness proxy — the SAME kernels the GPU backends
# dispatch through) and as far as compile/dispatch + graceful-skip. Real
# CUDA/AMD on-device execution is hardware-unvalidated and explicitly skipped.
# ==========================================================================

using Test
using CLM
import KernelAbstractions as KA
import Adapt

@testset "Backend selection API" begin
    # ---- defaults + introspection ------------------------------------------
    @test CLM.clm_backend() == :cpu                      # defaults to CPU
    @test CLM.clm_backend_registered(:cpu)
    @test CLM.clm_backend_functional(:cpu)               # CPU always usable
    @test CLM.clm_ka_backend(backend = :cpu) isa KA.CPU
    @test CLM.clm_device_array([1.0, 2.0]) == [1.0, 2.0] # CPU adapt = identity

    # ---- unknown / unregistered backends ----------------------------------
    @test_throws ArgumentError CLM.clm_set_backend(:nonsense)
    if !CLM.clm_backend_registered(:cuda)
        @test !CLM.clm_backend_functional(:cuda)          # no package -> not functional
        @test_throws ErrorException CLM.clm_set_backend(:cuda)        # informative error, NO hard dep
        @test_throws ErrorException CLM.clm_device_array([1.0]; backend = :cuda)
        @test_throws ErrorException CLM.clm_ka_backend(backend = :cuda)
    end
    if !CLM.clm_backend_registered(:amdgpu)
        @test !CLM.clm_backend_functional(:amdgpu)
        @test_throws ErrorException CLM.clm_set_backend(:amdgpu)
    end

    # ---- set/get round-trips on CPU ---------------------------------------
    @test CLM.clm_set_backend(:cpu) == :cpu
    @test CLM.clm_backend() == :cpu
end

# A stand-in "device" array type + a fake GPU backend, to exercise the FULL
# registration/dispatch path that a real CUDA/AMDGPU extension drives — WITHOUT
# any GPU package. This proves the extension seam compiles & dispatches.
struct FakeBackend <: KA.Backend end
struct FakeGPU{T,N} <: AbstractArray{T,N}
    d::Array{T,N}
end
Base.size(a::FakeGPU) = size(a.d)
Base.getindex(a::FakeGPU, i...) = getindex(a.d, i...)
Base.setindex!(a::FakeGPU, v, i...) = setindex!(a.d, v, i...)
Adapt.adapt_storage(::Type{<:FakeGPU}, x::Array) = FakeGPU(x)

@testset "Extension registration seam (fake GPU, no real device)" begin
    # Register a fake :cuda backend exactly the way ext/CLMCUDAExt.__init__ does.
    CLM._register_backend!(:cuda,
        x -> Adapt.adapt(FakeGPU, x),
        () -> FakeBackend(),
        () -> true;                       # pretend "functional"
        device_count = () -> 4,
        bind_device  = idx -> idx)        # record the bound index
    try
        @test CLM.clm_backend_registered(:cuda)
        @test CLM.clm_backend_functional(:cuda)
        @test CLM.clm_set_backend(:cuda) == :cuda
        @test CLM.clm_backend() == :cuda
        @test CLM.clm_ka_backend() isa FakeBackend
        @test CLM.clm_device_array([1.0, 2.0, 3.0]) isa FakeGPU
        # multi-GPU plumbing (multigpu.jl): single rank -> node-local rank 0,
        # device index = 0 % 4 = 0.
        @test CLM.clm_node_local_rank() == 0
        @test CLM.clm_gpu_device_count() == 4
        @test CLM.clm_bind_local_gpu!() == 0
    finally
        # Remove the fake registration so we don't leak it into other tests.
        delete!(CLM._BACKEND_REGISTRY, :cuda)
        CLM.clm_set_backend(:cpu)
    end
    @test !CLM.clm_backend_registered(:cuda)   # clean slate restored
end

@testset "Kernel parity through the backend API (CPU proxy)" begin
    # The CPU backend is the correctness proxy: the SAME kernels the GPU backends
    # dispatch through, taken from the output array. Run a real CLM kernel through
    # the backend-selection API on CPU and assert it matches the direct call.
    CLM.clm_set_backend(:cpu)
    nc = 5
    gridcell = collect(1:nc)
    vp_grc   = [200.0 + 3.0 * g for g in 1:nc]
    pbot_col = fill(1.0e5, nc)

    out_direct = zeros(nc)
    CLM.compute_forc_q!(out_direct, gridcell, vp_grc, pbot_col)

    # Move inputs+output "to device" via the centralized API (CPU = identity),
    # then run the identical kernel — same code path the GPU dispatch takes.
    out_api = CLM.clm_device_array(zeros(nc))
    CLM.compute_forc_q!(out_api,
                        CLM.clm_device_array(gridcell),
                        CLM.clm_device_array(vp_grc),
                        CLM.clm_device_array(pbot_col))
    @test out_api == out_direct                       # byte-identical on CPU
    @test KA.get_backend(out_api) isa KA.CPU
end

@testset "Multi-GPU-over-MPI run wiring on CPU (byte-identical to distributed)" begin
    # clm_run_multigpu! on :cpu must equal clm_run_distributed! + gather (device
    # move is the identity). Set up a real single-rank decomposition so the
    # rank-local clump loop has work, then check the no-op binding + the run
    # wrapper drive phys! the right number of steps and gather identically.
    CLM.clm_set_backend(:cpu)
    @test CLM.clm_bind_local_gpu!() == 0          # CPU bind is a no-op -> 0
    @test CLM.clm_gpu_device_count() == 0         # no GPUs on CPU backend

    numg = 4
    ncols = [1, 2, 1, 2]                          # cols per gridcell
    d = CLM.DecompData()
    CLM.decompInit_distributed!(numg; clump_pproc = 1, ncols_per_g = ncols, decomp_data = d)

    calls = Ref(0)
    phys!(bounds) = (calls[] += 1; nothing)
    res = CLM.clm_run_multigpu!(phys!, 3; decomp_data = d)
    @test res.bound_device == 0
    @test calls[] == 3                            # ran 3 steps over the (single) rank-local clump
    @test res.gathered == Any[]                   # no gather_fields requested

    # gather a column field on CPU == gather_to_master directly (identity move).
    nc = sum(ncols)
    colvals = [10.0 + c for c in 1:nc]
    res2 = CLM.clm_run_multigpu!(phys!, 1; decomp_data = d,
                                 gather_fields = ((colvals, CLM.SUBGRID_LEVEL_COLUMN),))
    direct = CLM.gather_to_master(colvals, CLM.SUBGRID_LEVEL_COLUMN; decomp_data = d)
    @test res2.gathered[1] == direct              # multi-GPU gather == distributed gather on CPU
end

# ---- real GPU paths: hardware-unvalidated on this machine -> skip ----------
@testset "Real GPU execution (skipped — no CUDA/AMD device here)" begin
    for sym in (:cuda, :amdgpu, :metal)
        if CLM.clm_backend_functional(sym)
            # Would run on real hardware; not reachable on this Metal-only/CPU box
            # for :cuda/:amdgpu. Left as a hook for a GPU CI lane.
            @test CLM.clm_set_backend(sym; require_functional = true) == sym
            CLM.clm_set_backend(:cpu)
        else
            @test_skip "backend :$sym not functional on this machine (no device / package)"
        end
    end
end
