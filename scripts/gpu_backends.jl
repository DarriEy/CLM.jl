# ==========================================================================
# gpu_backends.jl — shared GPU backend auto-detection for the validation scripts.
#
# Included by every scripts/gpu_validate_*.jl harness. Loads whichever GPU backend
# package is functional and exposes a device-array converter, the device array TYPE
# (for adapting whole state structs), a synchronize, and the working float precision.
# No harness should name a vendor type (MtlArray/CuArray/…) directly — use
# `device_array_type()` / `device_adapt` / `device_synchronize` / `gpu_functional`
# so the same harness runs on any backend.
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

# The selected backend's module and device-array TYPE, filled in by detect_backend().
# The e2e harnesses need the TYPE (not just a converter) because they adapt whole
# state structs with `Adapt.adapt(device_array_type(), state)`. Held in Refs and read
# through accessors so the scripts never name a vendor type directly.
const _BACKEND_MODULE = Ref{Any}(nothing)
const _BACKEND_ARRAY_TYPE = Ref{Any}(nothing)

"""Device array type of the detected backend (CuArray / ROCArray / oneArray / MtlArray)."""
device_array_type() = _BACKEND_ARRAY_TYPE[] === nothing ?
    error("no GPU backend detected — call detect_backend() first") : _BACKEND_ARRAY_TYPE[]

"""Adapt an entire state struct (or array) onto the detected device."""
device_adapt(x) = Base.invokelatest(() -> CLM.Adapt.adapt(device_array_type(), x))

"""Block until queued device work has finished. No-op when running on the host."""
function device_synchronize()
    m = _BACKEND_MODULE[]
    m === nothing && return nothing
    Base.invokelatest(() -> m.synchronize())
    return nothing
end

# Returns (name, to_device, FT) where to_device(::Array)->device array and FT is
# the float precision the backend supports, or nothing if no GPU is functional.
# First functional backend wins, NVIDIA → AMD → Intel → Apple. Also records the
# backend module + array type for device_array_type/device_adapt/device_synchronize.
#
# NOTE: call this as a TOP-LEVEL statement before using any device methods — the
# `using <Backend>` must advance the world age before the kernels run, otherwise
# device methods (size, get_backend, …) are "too new to be called".
const _DETECTED = Ref{Any}(nothing)
const _DETECT_DONE = Ref(false)

function detect_backend()
    _DETECT_DONE[] && return _DETECTED[]      # cached: detection runs at most once
    result = nothing
    # NOTE: `CuArray`, not `cu` — `cu(A)` silently demotes Float64 to Float32, which
    # would compare a Float32 device result against the Float64 CPU reference.
    for (name, ctor, FT) in (("CUDA",   :CuArray,  Float64),
                             ("AMDGPU", :ROCArray, Float64),
                             ("oneAPI", :oneArray, Float64),
                             ("Metal",  :MtlArray, Float32))
        m = _try_backend(Symbol(name))
        m === nothing && continue
        _BACKEND_MODULE[] = m
        _BACKEND_ARRAY_TYPE[] = Base.invokelatest(() -> getfield(m, ctor))
        result = (name, _converter(m, ctor), FT)
        break
    end
    _DETECT_DONE[] = true
    _DETECTED[] = result
    return result
end

"""True when a GPU backend was detected and is usable (replaces `Metal.functional()`)."""
gpu_functional() = detect_backend() !== nothing

"""Name of the detected backend, or "none"."""
gpu_backend_name() = (b = detect_backend()) === nothing ? "none" : b[1]

# Device arrays disallow BitArray and require Bool masks as Vector{Bool}.
_dev(to, x::BitArray) = to(collect(Bool, x))
_dev(to, x::AbstractArray) = to(x)

# Detect NOW, as a top-level statement of this included file, so the `using <Backend>`
# has advanced the world age before the including script's own top-level code runs.
# Later detect_backend() calls return the cached result, so harnesses that call it
# themselves are unaffected.
detect_backend()
