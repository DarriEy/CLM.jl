# ==========================================================================
# test_backend.jl — centralized backend-selection API + multi-GPU-over-MPI
# plumbing (src/infrastructure/backend.jl + multigpu.jl + ext/CLM{CUDA,AMDGPU}Ext).
#
# On a CPU-only box the GPU paths are validated through the CPU backend (the
# correctness proxy — the SAME kernels the GPU backends dispatch through) and as
# far as compile/dispatch + graceful-skip; the real-device testset skips. With a
# GPU package loaded on real hardware (e.g. `using CUDA` on an NVIDIA box) the
# last testset runs a kernel on the device. Validated on CUDA (RTX 5070 Ti,
# sm_120) 2026-10 — see docs/GPU_REVALIDATION_2026-10_CUDA.md.
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
    # If a REAL :cuda backend is already registered (CUDA.jl loaded on an NVIDIA
    # box), keep it: the fake registration below overwrites it, and restoring
    # only by `delete!` would silently un-register the real device and turn the
    # "Real GPU execution" testset into a skip on exactly the hardware it targets.
    prior_cuda = get(CLM._BACKEND_REGISTRY, :cuda, nothing)
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
        # Remove the fake registration so we don't leak it into other tests,
        # restoring a real one if it was there before.
        if prior_cuda === nothing
            delete!(CLM._BACKEND_REGISTRY, :cuda)
        else
            CLM._BACKEND_REGISTRY[:cuda] = prior_cuda
        end
        CLM.clm_set_backend(:cpu)
    end
    # Clean slate restored: either unregistered, or the real (non-fake) backend.
    @test CLM.clm_backend_registered(:cuda) == (prior_cuda !== nothing)
    prior_cuda === nothing || @test !(CLM.clm_ka_backend(backend = :cuda) isa FakeBackend)
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

# ---- real GPU paths: run only where the backend package + device exist ------
# Skips on a CPU-only box; on a GPU box (with CUDA.jl / AMDGPU.jl / Metal.jl
# loaded into the session) it selects the real backend and, for the Float64
# backends, runs a real kernel on the device against the CPU result.
@testset "Real GPU execution (skips without a device)" begin
    for sym in (:cuda, :amdgpu, :metal)
        if CLM.clm_backend_functional(sym)
            @test CLM.clm_set_backend(sym; require_functional = true) == sym
            if sym in (:cuda, :amdgpu)            # Float64 hardware (Metal is Float32-only)
                nc = 7
                gridcell = collect(1:nc)
                vp_grc   = [200.0 + 3.0 * g for g in 1:nc]
                pbot_col = [1.0e5 - 50.0 * g for g in 1:nc]
                ref = zeros(nc)
                CLM.compute_forc_q!(ref, gridcell, vp_grc, pbot_col)
                out = CLM.clm_device_array(zeros(nc))
                @test eltype(out) == Float64      # no silent Float32 demotion
                CLM.compute_forc_q!(out, CLM.clm_device_array(gridcell),
                                    CLM.clm_device_array(vp_grc),
                                    CLM.clm_device_array(pbot_col))
                @test !(KA.get_backend(out) isa KA.CPU)
                @test Array(out) == ref           # device == CPU, bit-identical
            end
            CLM.clm_set_backend(:cpu)
        else
            @test_skip "backend :$sym not functional on this machine (no device / package)"
        end
    end
end
