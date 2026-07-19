# ==========================================================================
# CLMCUDAExt.jl — CLM ⇄ CUDA package extension (NVIDIA GPUs)
#
# Loaded automatically by Julia ONLY when both `CLM` and `CUDA` are present in
# the environment (weakdep + [extensions] in Project.toml). It installs the
# `:cuda` backend into CLM's centralized backend registry
# (src/infrastructure/backend.jl) by handing `_register_backend!` closures that
# name the CUDA types — so the BASE package never depends on CUDA.
#
# Hooks installed:
#   adapt_to(x)      = Adapt.adapt(CUDA.CuArray, x)      # generalizes adapt(MtlArray,x)
#   ka_backend()     = CUDA.CUDABackend()                # KernelAbstractions backend
#   functional()     = CUDA.functional()                # real device present?
#   device_count()   = length(CUDA.devices())           # multi-GPU
#   bind_device(idx) = CUDA.device!(idx)                 # one-GPU-per-rank bind
#
# NOTE (hardware honesty): authored on a Metal-only machine. The closures are the
# documented CUDA.jl API but have NOT been executed on real NVIDIA hardware here.
# Once `import CUDA` succeeds on a CUDA box, `CLM.clm_set_backend(:cuda)` +
# `CLM.clm_device_array(...)` + `CLM.clm_run_multigpu!(...)` light up through this
# seam with no base-package changes.
# ==========================================================================
module CLMCUDAExt

using CLM
using CUDA
import Adapt

# ---------------------------------------------------------------------------
# Device-only override: `CLM._smooth_f64()` reads the host-side `SMOOTH_MODE`
# Ref (src/infrastructure/smooth_ad.jl). The Float64-specific smooth_* methods
# consult it, and a GPU kernel cannot dereference a host pointer — the result is
# ERROR_ILLEGAL_ADDRESS (code 700), not a catchable error.
#
# smooth_ad.jl documented this as GPU-safe on the grounds that "device kernels
# are Float32 and dispatch to the generic type-based path". That held on Metal
# (Float32-only hardware). CUDA runs Float64, so Float64 args select the method
# that reads the global and every kernel calling smooth_min/max/clamp/abs/ifelse
# faults.
#
# Overriding to `false` on device reproduces the default `:auto` behaviour
# (Float64 -> exact min/max), which is what the forward physics already assumes.
# The `:always` override stays host-only, exactly as smooth_ad.jl intends: this
# override is invisible to host code, so calibration FD/AD-consistency and
# Enzyme reverse-mode are unaffected.
CUDA.@device_override CLM._smooth_f64() = false

function __init__()
    CLM._register_backend!(:cuda,
        x -> Adapt.adapt(CUDA.CuArray, x),   # adapt_to
        () -> CUDA.CUDABackend(),            # ka_backend
        () -> CUDA.functional();             # functional
        device_count = () -> length(CUDA.devices()),
        bind_device  = idx -> CUDA.device!(idx),
    )
    return nothing
end

end # module CLMCUDAExt
