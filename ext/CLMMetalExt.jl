# ==========================================================================
# CLMMetalExt.jl — CLM ⇄ Metal package extension (Apple Silicon GPUs)
#
# Loaded automatically by Julia ONLY when both `CLM` and `Metal` are present in
# the environment (weakdep + [extensions] in Project.toml). It installs the
# `:metal` backend into CLM's centralized backend registry
# (src/infrastructure/backend.jl), mirroring CLMCUDAExt / CLMAMDGPUExt.
#
# Hooks installed:
#   adapt_to(x)      = Adapt.adapt(Metal.MtlArray, x)
#   ka_backend()     = Metal.MetalBackend()
#   functional()     = Metal.functional()
#   device_count()   = length(Metal.devices())
#   bind_device(idx) = Metal.device!(Metal.devices()[idx])
#
# Metal hardware has no Float64: build Float32 state first
# (`make_instances(clm_float_type(backend = :metal))`) before `adapt_to`. The
# Metal validation fleet keeps using its own Float32 down-adaptor
# (scripts/gpu_adapt.jl); this extension is what lets `clm_set_backend(:metal)`
# select the device through the same seam as CUDA/AMDGPU.
# ==========================================================================
module CLMMetalExt

using CLM
using Metal
import Adapt

function __init__()
    CLM._register_backend!(:metal,
        x -> Adapt.adapt(Metal.MtlArray, x),  # adapt_to
        () -> Metal.MetalBackend(),           # ka_backend
        () -> Metal.functional();             # functional
        device_count = () -> length(Metal.devices()),
        bind_device  = idx -> Metal.device!(Metal.devices()[idx]),
    )
    return nothing
end

end # module CLMMetalExt
