# ==========================================================================
# gpu_backends.jl — shared GPU backend auto-detection for the validation scripts.
#
# Included by scripts/gpu_validate.jl (kernel parity) and scripts/gpu_ad_validate.jl
# (AD parity). Loads whichever GPU backend package is functional and exposes a
# device-array converter plus the backend's working float precision.
#
#   CUDA   — Float64   (NVIDIA)
#   AMDGPU — Float64   (AMD ROCm)
#   oneAPI — Float64   (Intel; Data-Center GPUs do f64)
#   Metal  — Float32   (Apple Silicon; no Float64 in hardware)
# ==========================================================================

# Load `pkg` on demand and return its module object if its GPU is functional, else
# nothing. `@eval using $pkg` creates the binding in a *newer* world age than this
# function body runs in, so every access to the just-loaded module is deferred
# through `invokelatest` (otherwise: "binding too new" world-age errors).
function _try_backend(pkg::Symbol)
    try
        @eval using $pkg
        mod = Base.invokelatest(getfield, @__MODULE__, pkg)
        Base.invokelatest(() -> mod.functional()) ? mod : nothing
    catch
        nothing
    end
end

# Device-array converter for a loaded backend module. The property lookup is done
# inside the invokelatest closure so it, too, sees the freshly-loaded module.
_converter(mod, ctor::Symbol) = x -> Base.invokelatest(() -> getfield(mod, ctor)(x))

# Returns (name, to_device, FT) where to_device(::Array)->device array and FT is
# the float precision the backend supports, or nothing if no GPU is functional.
# First functional backend wins, NVIDIA → AMD → Intel → Apple.
#
# NOTE: call this as a TOP-LEVEL statement before using any device methods — the
# `using <Backend>` must advance the world age before the kernels run, otherwise
# device methods (size, get_backend, …) are "too new to be called".
function detect_backend()
    if (m = _try_backend(:CUDA))   !== nothing; return ("CUDA",   _converter(m, :cu),       Float64); end
    if (m = _try_backend(:AMDGPU)) !== nothing; return ("AMDGPU", _converter(m, :ROCArray), Float64); end
    if (m = _try_backend(:oneAPI)) !== nothing; return ("oneAPI", _converter(m, :oneArray), Float64); end
    if (m = _try_backend(:Metal))  !== nothing; return ("Metal",  _converter(m, :MtlArray), Float32); end
    return nothing
end

# Device arrays disallow BitArray and require Bool masks as Vector{Bool}.
_dev(to, x::BitArray) = to(collect(Bool, x))
_dev(to, x::AbstractArray) = to(x)
