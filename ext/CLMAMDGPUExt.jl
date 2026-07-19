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

# Device-only override: `CLM._smooth_f64()` reads the host-side `SMOOTH_MODE` Ref
# (src/infrastructure/smooth_ad.jl), and the Float64-specific smooth_* methods
# consult it. A device kernel cannot dereference a host pointer — on CUDA this
# manifests as ERROR_ILLEGAL_ADDRESS (see ext/CLMCUDAExt.jl for the full write-up
# and the three-kernel reproducer). AMDGPU is likewise a Float64 backend, so it
# selects the same methods and hits the same fault.
#
# NOTE (hardware honesty): this file has never run on real AMD hardware, so unlike
# the CUDA override this one is UNVERIFIED. It is `@static`-guarded so that a
# ROCm box on which AMDGPU does not export `@device_override` still precompiles
# the extension instead of failing outright — if the guard skips, the smooth_*
# kernels will fault and this is the first place to look.
@static if isdefined(AMDGPU, Symbol("@device_override"))
    AMDGPU.@device_override CLM._smooth_f64() = false
end

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
