# ==========================================================================
# Smooth AD primitives for ForwardDiff compatibility
#
# These functions provide smooth (C∞) approximations to min/max/clamp/
# heaviside that are identical to standard functions for Float64 but use
# smooth LogSumExp-based approximations for ForwardDiff.Dual types.
#
# The smoothing parameter k controls sharpness: larger k = closer to
# exact but sharper transition. k=50 gives ~1e-22 approximation error.
# ==========================================================================

# No ForwardDiff import needed — Float64 methods are exact,
# generic Real methods handle Dual via dispatch.

# --------------------------------------------------------------------------
# Smooth evaluation mode control
# --------------------------------------------------------------------------

"""
    SMOOTH_MODE

Controls when smooth approximations are used:
- `:auto` (default) — smooth for Dual types, exact for Float64
- `:always` — smooth for ALL types (use for FD/AD gradient consistency)
- `:never` — exact for all types
"""
const SMOOTH_MODE = Ref(:auto)

"""
    _is_ad_type(::Type{T})
    _is_ad_type(x::AbstractArray)

Check if type T is an AD dual number (not plain Float64/Float32).
ForwardDiff.Dual <: Real but NOT <: AbstractFloat.
"""
_is_ad_type(::Type{T}) where T = !(T <: AbstractFloat)
_is_ad_type(x::AbstractArray) = _is_ad_type(eltype(x))

"""Return true if we should use the smooth path for the given type."""
function _use_smooth(::Type{T}) where T
    mode = SMOOTH_MODE[]
    mode === :always && return true
    mode === :never && return false
    # :auto — smooth only for non-AbstractFloat (i.e. Dual types)
    return _is_ad_type(T)
end

# --------------------------------------------------------------------------
# smooth_min — LogSumExp approximation: -(log(exp(-k*a) + exp(-k*b)))/k
# --------------------------------------------------------------------------

"""
    smooth_min(a, b; k=50.0)

Smooth approximation to `min(a, b)`.
- For `Float64` arguments: exact `min(a, b)` (no approximation error).
- For `ForwardDiff.Dual`: LogSumExp smooth approximation with sharpness `k`.
"""
function smooth_min(a::Float64, b::Float64; k::Float64=50.0)
    _use_smooth(Float64) || return min(a, b)
    ak = -k * a; bk = -k * b; m = max(ak, bk)
    return -(m + log(exp(ak - m) + exp(bk - m))) / k
end
smooth_min(a::Float64, b::Int; k::Float64=50.0) = smooth_min(a, Float64(b); k=k)
smooth_min(a::Int, b::Float64; k::Float64=50.0) = smooth_min(Float64(a), b; k=k)

function smooth_min(a::T, b::S; k::Real=50.0) where {T<:Real, S<:Real}
    R = promote_type(T, S)
    ak = -k * a
    bk = -k * b
    # Numerically stable: shift by max to avoid overflow
    m = max(ak, bk)
    return -R(m + log(exp(ak - m) + exp(bk - m))) / k
end

# --------------------------------------------------------------------------
# smooth_max — LogSumExp approximation: (log(exp(k*a) + exp(k*b)))/k
# --------------------------------------------------------------------------

"""
    smooth_max(a, b; k=50.0)

Smooth approximation to `max(a, b)`.
- For `Float64` arguments: exact `max(a, b)`.
- For `ForwardDiff.Dual`: LogSumExp smooth approximation with sharpness `k`.
"""
function smooth_max(a::Float64, b::Float64; k::Float64=50.0)
    _use_smooth(Float64) || return max(a, b)
    ak = k * a; bk = k * b; m = max(ak, bk)
    return (m + log(exp(ak - m) + exp(bk - m))) / k
end
smooth_max(a::Float64, b::Int; k::Float64=50.0) = smooth_max(a, Float64(b); k=k)
smooth_max(a::Int, b::Float64; k::Float64=50.0) = smooth_max(Float64(a), b; k=k)

function smooth_max(a::T, b::S; k::Real=50.0) where {T<:Real, S<:Real}
    R = promote_type(T, S)
    ak = k * a
    bk = k * b
    m = max(ak, bk)
    return R(m + log(exp(ak - m) + exp(bk - m))) / k
end

# --------------------------------------------------------------------------
# smooth_clamp — smooth_max(lo, smooth_min(x, hi))
# --------------------------------------------------------------------------

"""
    smooth_clamp(x, lo, hi; k=50.0)

Smooth approximation to `clamp(x, lo, hi)`.
- For `Float64` arguments: exact `clamp(x, lo, hi)`.
- For `ForwardDiff.Dual`: composition of smooth_max and smooth_min.
"""
function smooth_clamp(x::Float64, lo::Float64, hi::Float64; k::Float64=50.0)
    _use_smooth(Float64) || return clamp(x, lo, hi)
    return smooth_max(lo, smooth_min(x, hi; k=k); k=k)
end

function smooth_clamp(x::T, lo, hi; k::Real=50.0) where {T<:Real}
    return smooth_max(lo, smooth_min(x, hi; k=k); k=k)
end

# --------------------------------------------------------------------------
# smooth_heaviside — sigmoid: 1/(1 + exp(-k*x))
# --------------------------------------------------------------------------

"""
    smooth_heaviside(x; k=50.0)

Smooth approximation to the Heaviside step function.
- For `Float64`: exact `x >= 0 ? 1.0 : 0.0`.
- For `ForwardDiff.Dual`: sigmoid `1/(1 + exp(-k*x))`.
"""
function smooth_heaviside(x::Float64; k::Float64=50.0)
    _use_smooth(Float64) || return x >= 0.0 ? 1.0 : 0.0
    return 1.0 / (1.0 + exp(-k * x))
end

function smooth_heaviside(x::T; k::Real=50.0) where {T<:Real}
    return one(T) / (one(T) + exp(-k * x))
end

# --------------------------------------------------------------------------
# smooth_abs — x * tanh(k*x)
# --------------------------------------------------------------------------

"""
    smooth_abs(x; k=50.0)

Smooth approximation to `abs(x)`.
- For `Float64`: exact `abs(x)`.
- For `ForwardDiff.Dual`: `x * tanh(k*x)`.
"""
function smooth_abs(x::Float64; k::Float64=50.0)
    _use_smooth(Float64) || return abs(x)
    return x * tanh(k * x)
end

function smooth_abs(x::T; k::Real=50.0) where {T<:Real}
    return x * tanh(k * x)
end

# --------------------------------------------------------------------------
# smooth_ifelse — sigmoid-blended conditional: a if x≥0, else b
# --------------------------------------------------------------------------

"""
    smooth_ifelse(x, a, b; k=50.0)

Smooth approximation to `x >= 0 ? a : b`.
- For `Float64` x: exact branch.
- For `ForwardDiff.Dual` x: sigmoid blend `a*h + b*(1-h)` where `h = σ(k*x)`.

Common usage: `smooth_ifelse(T - TFRZ, unfrozen_val, frozen_val)`.
"""
function smooth_ifelse(x::Float64, a, b; k::Float64=50.0)
    _use_smooth(Float64) || return x >= 0.0 ? a : b
    h = 1.0 / (1.0 + exp(-k * x))
    return a * h + b * (1.0 - h)
end

function smooth_ifelse(x::T, a, b; k::Real=50.0) where {T<:Real}
    h = one(T) / (one(T) + exp(-k * x))
    return a * h + b * (one(T) - h)
end
