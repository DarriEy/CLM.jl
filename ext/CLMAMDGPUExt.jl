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
# NOTE (hardware honesty): the closures are the documented AMDGPU.jl API. With
# AMDGPU.jl installed on a machine with NO ROCm device, the extension precompiles,
# loads and registers :amdgpu (test/test_amdgpu_backend.jl, docs/GPU_AMDGPU.md);
# none of it has yet been executed on real AMD hardware.
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
# selects the same methods and would hit the same fault.
#
# `@device_override` lives in `AMDGPU.Device` (src/device/gcn/helpers.jl) and is
# NOT re-exported from the top-level `AMDGPU` module — AMDGPU's own ROCKernels.jl
# does `import AMDGPU.Device: @device_override`. An earlier version of this file
# guarded on `isdefined(AMDGPU, Symbol("@device_override"))`, which is false for
# AMDGPU.jl 1.x/2.x, so the override was silently skipped. The macro expands to an
# overlay on `AMDGPU.method_table`, which resolves here because of `using AMDGPU`.
#
# The guard below only protects precompilation against a future AMDGPU.jl that
# moves the macro; `_SMOOTH_OVERRIDE_INSTALLED` records which branch ran so the
# test suite can assert the override is actually in place.
@static if isdefined(AMDGPU, :Device) && isdefined(AMDGPU.Device, Symbol("@device_override"))
    AMDGPU.Device.@device_override CLM._smooth_f64() = false
    const _SMOOTH_OVERRIDE_INSTALLED = true
else
    const _SMOOTH_OVERRIDE_INSTALLED = false
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
