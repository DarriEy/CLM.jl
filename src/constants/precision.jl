# ==========================================================================
# Working floating-point precision selection.
#
# Policy (see `make_instances`): Float64 on CPU and on Float64-capable GPUs
# (CUDA / AMD / Intel) for full fidelity against the Fortran reference; Float32
# ONLY for the Apple Metal backend, whose hardware has no double precision at all.
#
# Override the policy with the `CLM_FLOAT` environment variable ("Float32" or
# "Float64"), e.g. to force a Float32 CPU run for testing the reduced-precision
# path on a machine without a Metal GPU.
# ==========================================================================

"""
    clm_float_type(; backend::Symbol = :cpu) -> Type{<:AbstractFloat}

Return the working precision for a compute `backend`:

- `:metal` ⇒ `Float32` (Apple GPUs are Float32-only)
- everything else (`:cpu`, `:cuda`, `:amdgpu`, `:oneapi`) ⇒ `Float64`

The `CLM_FLOAT` environment variable, when set to `"Float32"` or `"Float64"`,
overrides this policy. Pass the result to [`make_instances`](@ref) to build a
state tree at the chosen precision.
"""
function clm_float_type(; backend::Symbol = :cpu)
    env = get(ENV, "CLM_FLOAT", "")
    env == "Float32" && return Float32
    env == "Float64" && return Float64
    return backend === :metal ? Float32 : Float64
end
