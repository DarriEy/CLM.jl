# ==========================================================================
# backend.jl — centralized compute-backend selection (CPU / Metal / CUDA / AMDGPU)
#
# This is the single entry point for choosing WHERE CLM state arrays live and
# which KernelAbstractions backend the physics kernels dispatch through. It
# generalizes the ad-hoc `adapt(MtlArray, ...)` pattern used in the Metal
# validation scripts into one API:
#
#   * `clm_set_backend(:cpu | :metal | :cuda | :amdgpu)` — select the active
#     backend (defaults to `:cpu`). Validates the GPU package is loaded.
#   * `clm_backend()` — the currently selected backend symbol.
#   * `clm_device_array(x)` — move an array / state struct onto the active
#     backend's array type (via `Adapt`), the generalization of
#     `adapt(MtlArray, x)`. CPU is the identity.
#   * `clm_ka_backend()` — the `KernelAbstractions.Backend` for the active
#     backend (`KA.CPU()` on CPU; the GPU package's backend otherwise).
#   * `clm_backend_functional(sym)` — is that backend usable right now (package
#     loaded AND a device present)?
#
# NO HARD GPU DEPENDENCY
# ----------------------
# The base package depends on NEITHER CUDA, AMDGPU, nor Metal. The device-array
# *adaptor type* and the *KA backend object* for each GPU are provided by
# PACKAGE EXTENSIONS (`ext/CLMCUDAExt.jl`, `ext/CLMAMDGPUExt.jl`, and Metal is
# handled the same way by the validation scripts that `import Metal`). Each
# extension, when its weakdep is loaded, calls `_register_backend!` to install
# the adaptor closure + KA backend constructor into the tables below. Until an
# extension is loaded, selecting that backend raises an informative error rather
# than introducing a load-time dependency. `using CLM` therefore works with none
# of the GPU packages installed, and CI (no GPU) stays green.
#
# This mirrors how Metal is kept optional today (scripts `import Metal` only when
# present); the extensions make CUDA/AMDGPU first-class through the same seam.
# ==========================================================================

import Adapt
import KernelAbstractions as KA

# The set of backend symbols CLM understands. `:metal` is registered by the
# Metal validation scripts (or a user `import`ing Metal) the same way the
# CUDA/AMDGPU extensions register theirs.
const CLM_BACKENDS = (:cpu, :metal, :cuda, :amdgpu)

# Registry of GPU backends installed by extensions. Each entry maps a backend
# symbol to a NamedTuple of closures so the base package never names a GPU type:
#   adapt_to     :: (x) -> device-array(x)     generalizes adapt(CuArray, x)
#   ka_backend   :: () -> KA.Backend           the KernelAbstractions backend
#   functional   :: () -> Bool                  package loaded AND device present
#   device_count :: () -> Int                   number of visible GPUs (multi-GPU)
#   bind_device  :: (idx::Int) -> Any           bind this rank to GPU `idx`
# The last two power the one-GPU-per-rank multi-GPU split (multigpu.jl). CPU is
# always present and handled inline (not in this registry). The NamedTuple is
# kept loosely typed so extensions can omit the optional multi-GPU hooks.
const _BACKEND_REGISTRY = Dict{Symbol,NamedTuple}()

# The currently selected backend (defaults to CPU — the byte-identical serial path).
const _ACTIVE_BACKEND = Ref{Symbol}(:cpu)

"""
    _register_backend!(sym, adapt_to, ka_backend, functional;
                       device_count = nothing, bind_device = nothing)

Called by a package extension (`ext/CLMCUDAExt.jl`, `ext/CLMAMDGPUExt.jl`) or a
Metal-aware script when its GPU package loads, to install that backend's hooks:

* `adapt_to(x)`      — move `x` to the GPU array type (e.g. `adapt(CuArray, x)`).
* `ka_backend()`     — return the `KernelAbstractions.Backend` (e.g. `CUDABackend()`).
* `functional()`     — `true` iff a real device is present and usable right now.
* `device_count()`   — (optional, multi-GPU) number of visible GPUs.
* `bind_device(idx)` — (optional, multi-GPU) bind this rank to GPU `idx` (0-based).

Returns `sym`. Idempotent (re-registration overwrites). This is the ONLY way a
GPU type enters the base package's dispatch, so there is no hard dependency.
"""
function _register_backend!(sym::Symbol, adapt_to::Function, ka_backend::Function,
                            functional::Function; device_count = nothing,
                            bind_device = nothing)
    sym in CLM_BACKENDS || throw(ArgumentError("unknown backend $sym (expected one of $CLM_BACKENDS)"))
    entry = (adapt_to = adapt_to, ka_backend = ka_backend, functional = functional)
    device_count === nothing || (entry = merge(entry, (; device_count = device_count)))
    bind_device  === nothing || (entry = merge(entry, (; bind_device = bind_device)))
    _BACKEND_REGISTRY[sym] = entry
    return sym
end

"""
    clm_backend() -> Symbol

The currently selected compute backend (`:cpu` by default).
"""
clm_backend() = _ACTIVE_BACKEND[]

"""
    clm_backend_registered(sym::Symbol) -> Bool

True iff backend `sym` has been installed (its GPU extension/package is loaded).
`:cpu` is always registered.
"""
clm_backend_registered(sym::Symbol) = sym === :cpu || haskey(_BACKEND_REGISTRY, sym)

"""
    clm_backend_functional(sym::Symbol = clm_backend()) -> Bool

True iff backend `sym` is usable RIGHT NOW: `:cpu` always; a GPU backend only
when its package is loaded AND a device is present and functional. Use this to
guard GPU code paths / tests so they skip gracefully where no device exists.
"""
function clm_backend_functional(sym::Symbol = clm_backend())
    sym === :cpu && return true
    haskey(_BACKEND_REGISTRY, sym) || return false
    try
        return _BACKEND_REGISTRY[sym].functional()::Bool
    catch
        return false
    end
end

"""
    clm_set_backend(sym::Symbol; require_functional::Bool = false) -> Symbol

Select the active compute backend, one of `$(CLM_BACKENDS)`. `:cpu` is always
available. A GPU backend requires its package/extension to be loaded
(`import CUDA` / `import AMDGPU` / `import Metal`) — otherwise an informative
error is raised (NO hard dependency is added by selecting it).

With `require_functional = true`, also asserts a real device is present
(`clm_backend_functional`), erroring otherwise; with `false` (default) a backend
may be *selected* for plumbing/compile purposes even where no device exists
(useful on CI to exercise the selection API without hardware).

Returns the selected symbol. This is the single, centralized replacement for the
scattered `adapt(MtlArray, ...)` calls.
"""
function clm_set_backend(sym::Symbol; require_functional::Bool = false)
    sym in CLM_BACKENDS || throw(ArgumentError("clm_set_backend: unknown backend $sym (expected one of $CLM_BACKENDS)"))
    if sym !== :cpu && !haskey(_BACKEND_REGISTRY, sym)
        error("clm_set_backend(:$sym): backend not registered. Load its package first " *
              "(`import $(sym === :cuda ? "CUDA" : sym === :amdgpu ? "AMDGPU" : "Metal")`), " *
              "which triggers the CLM extension that installs the $sym device adaptor + KA backend.")
    end
    if require_functional && !clm_backend_functional(sym)
        error("clm_set_backend(:$sym; require_functional=true): no functional $sym device present.")
    end
    _ACTIVE_BACKEND[] = sym
    return sym
end

"""
    clm_device_array(x; backend::Symbol = clm_backend())

Move an array (or any `Adapt`-registered state struct / `CLMInstances` tree) onto
`backend`'s array type. The generalization of `adapt(MtlArray, x)` to
CUDA/AMDGPU/Metal, and the identity on `:cpu`.

For a GPU backend the registered `adapt_to` closure (installed by the extension)
is applied; it wraps `Adapt.adapt(CuArray/ROCArray/MtlArray, x)`. Because state
structs are `Adapt.@adapt_structure`d, this recursively moves their array fields.

Note on precision: Apple Metal is Float32-only. Pair with
[`clm_float_type`](@ref) / [`make_instances`](@ref) to build the state at the
right precision BEFORE moving it (e.g. `make_instances(clm_float_type(backend=:metal))`),
exactly as the Metal validation scripts do.
"""
function clm_device_array(x; backend::Symbol = clm_backend())
    backend === :cpu && return x
    haskey(_BACKEND_REGISTRY, backend) ||
        error("clm_device_array: backend :$backend not registered (load its package first).")
    return _BACKEND_REGISTRY[backend].adapt_to(x)
end

"""
    clm_ka_backend(; backend::Symbol = clm_backend()) -> KernelAbstractions.Backend

The `KernelAbstractions` backend object for `backend`: `KA.CPU()` on `:cpu`, and
the GPU package's backend (`CUDABackend()` / `ROCBackend()` / `MetalBackend()`)
otherwise. The physics kernels already take their backend from the OUTPUT array
(`KA.get_backend(out)`), so this is mainly for code that must name the backend
explicitly (e.g. allocating device scratch with `KA.zeros(backend, ...)`).
"""
function clm_ka_backend(; backend::Symbol = clm_backend())
    backend === :cpu && return KA.CPU()
    haskey(_BACKEND_REGISTRY, backend) ||
        error("clm_ka_backend: backend :$backend not registered (load its package first).")
    return _BACKEND_REGISTRY[backend].ka_backend()::KA.Backend
end
