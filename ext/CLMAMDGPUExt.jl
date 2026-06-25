# ==========================================================================
# CLMAMDGPUExt.jl — CLM ⇄ AMDGPU package extension (AMD ROCm GPUs)
#
# Loaded automatically by Julia ONLY when both `CLM` and `AMDGPU` are present in
# the environment (weakdep + [extensions] in Project.toml). It installs the
# `:amdgpu` backend into CLM's centralized backend registry
# (src/infrastructure/backend.jl) via `_register_backend!` closures that name the
# AMDGPU types — so the BASE package never depends on AMDGPU.
#
# Hooks installed:
#   adapt_to(x)      = Adapt.adapt(AMDGPU.ROCArray, x)   # generalizes adapt(MtlArray,x)
#   ka_backend()     = AMDGPU.ROCBackend()               # KernelAbstractions backend
#   functional()     = AMDGPU.functional()               # real device present?
#   device_count()   = length(AMDGPU.devices())          # multi-GPU
#   bind_device(idx) = AMDGPU.device!(AMDGPU.devices()[idx+1])  # one-GPU-per-rank bind
#
# NOTE (hardware honesty): authored on a Metal-only machine. The closures are the
# documented AMDGPU.jl API but have NOT been executed on real AMD hardware here.
# `bind_device` uses 0-based `idx` (CLM's rank convention) to index AMDGPU's
# 1-based device list. Once `import AMDGPU` succeeds on a ROCm box,
# `CLM.clm_set_backend(:amdgpu)` lights up through this seam unchanged.
# ==========================================================================
module CLMAMDGPUExt

using CLM
using AMDGPU
import Adapt

function __init__()
    CLM._register_backend!(:amdgpu,
        x -> Adapt.adapt(AMDGPU.ROCArray, x),   # adapt_to
        () -> AMDGPU.ROCBackend(),              # ka_backend
        () -> AMDGPU.functional();              # functional
        device_count = () -> length(AMDGPU.devices()),
        bind_device  = idx -> AMDGPU.device!(AMDGPU.devices()[idx + 1]),
    )
    return nothing
end

end # module CLMAMDGPUExt
