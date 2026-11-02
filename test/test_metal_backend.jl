# ==========================================================================
# test_metal_backend.jl — the :metal (Apple Silicon) backend extension's
# multi-GPU hooks.
#
# Metal is a WEAKDEP: on a machine/environment without Metal.jl everything is
# @test_skip. With Metal.jl loaded (e.g. `--project=scripts` on a Mac) this
# checks the extension registers :metal and that `clm_bind_local_gpu!` — rank 0
# of a single process, CLM device index 0 — actually selects Metal device 1
# instead of hitting a BoundsError that is only logged as a warning (the bug
# fixed alongside this file).
# ==========================================================================

using Test
using CLM

const _METAL_INSTALLED = Base.find_package("Metal") !== nothing
_METAL_INSTALLED && @eval import Metal
const _METAL_FUNCTIONAL = _METAL_INSTALLED && Base.invokelatest(() -> Metal.functional())

@testset "Metal backend extension device bind" begin
    if !_METAL_INSTALLED
        @test_skip "Metal.jl not installed in this environment — :metal extension not exercised"
    else
        ext = Base.get_extension(CLM, :CLMMetalExt)
        @test ext !== nothing
        @test CLM.clm_backend_registered(:metal)
        @test CLM.clm_backend_functional(:metal) == _METAL_FUNCTIONAL
        if !_METAL_FUNCTIONAL
            @test_skip "no functional Metal device — bind not exercised"
        else
            devs = Metal.devices()
            ndev = length(devs)
            @test CLM.clm_gpu_device_count(backend = :metal) == ndev
            # Rank 0 -> device index 0 -> Metal.devices()[1], with NO warning
            # (the old hook threw a BoundsError that was caught and logged).
            idx = @test_logs CLM.clm_bind_local_gpu!(backend = :metal)
            @test idx == 0
            @test Metal.device() == devs[1]
            # Direct hook: last valid index binds, one past the end errors
            # (CUDA.device!(idx) semantics) rather than silently wrapping.
            entry = CLM._BACKEND_REGISTRY[:metal]
            entry.bind_device(ndev - 1)
            @test Metal.device() == devs[ndev]
            @test_throws ArgumentError entry.bind_device(ndev)
            entry.bind_device(0)
            # Kernels still launch after an explicit bind.
            x = Metal.MtlArray(Float32[1, 2, 3])
            @test Array(x .* 2f0) == Float32[2, 4, 6]
        end
    end
end
