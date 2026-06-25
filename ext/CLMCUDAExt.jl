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
