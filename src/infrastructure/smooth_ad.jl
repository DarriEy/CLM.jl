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

"""
Return true if the smooth path should be used for type `T`.

PURELY TYPE-BASED (no global read) so the smooth primitives are GPU-safe: a
device kernel cannot dereference a host-side mutable global. Float64/Float32 →
exact (false); ForwardDiff.Dual (and other non-AbstractFloat Reals) → smooth.

This is the former `:auto` behavior. The `:always`/`:never` runtime override
(SMOOTH_MODE) is inherently host-only and is handled at a higher level (the
generic Real methods always smooth, and callers that need Float64 to evaluate the
smooth function — e.g. calibration FD/AD consistency — force a smooth-dispatching
element type rather than toggling a global inside the kernel hot path).
"""
@inline _use_smooth(::Type{T}) where {T} = _is_ad_type(T)

"""
    _smooth_f64()

HOST-ONLY runtime override for the Float64-specific smooth primitives: returns true when
`SMOOTH_MODE[] === :always`, so plain-Float64 callers evaluate the SMOOTH function. This is
what makes `:always` actually take effect for Float64 — needed by (a) the calibration FD/AD-
consistency path (Float64 FD must match the Dual AD's smooth function) and (b) Enzyme REVERSE-
mode AD, which differentiates Float64 primals+shadows (no Dual type to trigger the type-based
smooth path). GPU-SAFE: only the Float64-specific methods consult this global; device kernels
are Float32 and dispatch to the generic type-based `_use_smooth` path, which never reads a host
global. Default `SMOOTH_MODE[] === :auto` → false → Float64 stays exact (suite/forward physics
unchanged).
"""
@inline _smooth_f64() = (SMOOTH_MODE[] === :always)

# Saturation threshold for the LogSumExp smooth_min/max. When one argument dominates the other
# by k·|a−b| > _SMOOTH_SAT, the dominated branch's softmax weight is < exp(−36) ≈ 2e-16 (below
# machine eps), so we return the HARD min/max exactly. This is not just speed: it makes the
# smooth functions SATURATE like the hard ones in the tails, which is essential for reverse-AD.
# A non-saturating smooth_max(0, f(x)) keeps a tiny-but-nonzero weight on the f(x) branch, so its
# reverse multiplies that weight by f'(x); if f has a SINGULAR derivative where it is dominated
# (e.g. log10(satw)→−∞ with d/dsatw = 1/(satw·ln10) → ∞ as satw→0, soil_temperature.jl:580), the
# product is (tiny)·∞ = NaN. Saturating to the hard value skips evaluating f' there entirely —
# the exact behavior of the hard branch. Also makes ±Inf arguments safe (no Inf−Inf in LogSumExp).
const _SMOOTH_SAT = 36.0

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
    (_use_smooth(Float64) || _smooth_f64()) || return min(a, b)
    d = k * (a - b)
    d >  _SMOOTH_SAT && return b        # b ≪ a → min is b (skip a's branch derivative)
    d < -_SMOOTH_SAT && return a        # a ≪ b → min is a
    ak = -k * a; bk = -k * b; m = max(ak, bk)
    return -(m + log(exp(ak - m) + exp(bk - m))) / k
end
smooth_min(a::Float64, b::Int; k::Float64=50.0) = smooth_min(a, Float64(b); k=k)
smooth_min(a::Int, b::Float64; k::Float64=50.0) = smooth_min(Float64(a), b; k=k)

function smooth_min(a::T, b::S; k::Real=50) where {T<:Real, S<:Real}
    R = promote_type(T, S)
    _use_smooth(R) || return min(a, b)   # exact for non-AD (Float32/Float64); GPU-safe
    kk = R(k)                            # sharpness at working precision (no Float64 pin)
    d = kk * (a - b)                     # saturation guard (see _SMOOTH_SAT); compares on value
    d >  R(_SMOOTH_SAT) && return R(b)
    d < -R(_SMOOTH_SAT) && return R(a)
    ak = -kk * a
    bk = -kk * b
    # Numerically stable: shift by max to avoid overflow
    m = max(ak, bk)
    return -(m + log(exp(ak - m) + exp(bk - m))) / kk
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
    (_use_smooth(Float64) || _smooth_f64()) || return max(a, b)
    d = k * (a - b)
    d >  _SMOOTH_SAT && return a        # a ≫ b → max is a (skip b's branch derivative)
    d < -_SMOOTH_SAT && return b        # a ≪ b → max is b
    ak = k * a; bk = k * b; m = max(ak, bk)
    return (m + log(exp(ak - m) + exp(bk - m))) / k
end
smooth_max(a::Float64, b::Int; k::Float64=50.0) = smooth_max(a, Float64(b); k=k)
smooth_max(a::Int, b::Float64; k::Float64=50.0) = smooth_max(Float64(a), b; k=k)

function smooth_max(a::T, b::S; k::Real=50) where {T<:Real, S<:Real}
    R = promote_type(T, S)
    _use_smooth(R) || return max(a, b)   # exact for non-AD (Float32/Float64); GPU-safe
    kk = R(k)                            # sharpness at working precision (no Float64 pin)
    d = kk * (a - b)                     # saturation guard (see _SMOOTH_SAT)
    d >  R(_SMOOTH_SAT) && return R(a)
    d < -R(_SMOOTH_SAT) && return R(b)
    ak = kk * a
    bk = kk * b
    m = max(ak, bk)
    return (m + log(exp(ak - m) + exp(bk - m))) / kk
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
    (_use_smooth(Float64) || _smooth_f64()) || return clamp(x, lo, hi)
    return smooth_max(lo, smooth_min(x, hi; k=k); k=k)
end

function smooth_clamp(x::T, lo, hi; k::Real=50) where {T<:Real}
    return smooth_max(lo, smooth_min(x, hi; k=k); k=k)
end

# --------------------------------------------------------------------------
# Numerically-stable logistic sigmoid σ(k*x) = 1/(1 + exp(-k*x))
#
# The naive form overflows: for k*x large and positive, exp(-k*x) underflows
# to 0 (fine), but for k*x large and NEGATIVE, exp(-k*x) overflows to Inf, and
# the ForwardDiff partial through that Inf becomes NaN. We branch on the sign of
# k*x so the exp argument is always <= 0 (in [0,1]), which is overflow-safe and
# keeps derivatives finite for arbitrarily large |x| (they simply saturate to 0).
# --------------------------------------------------------------------------
@inline function _stable_sigmoid(x::Real, k::Real)
    kx = k * x
    # Saturation guard (mirrors smooth_min/max _SMOOTH_SAT): for |kx| beyond the threshold the
    # sigmoid is 1/0 to machine precision AND its derivative k·σ(1−σ) is < ~1e-15. Returning the
    # exact 0/1 makes the derivative identically 0, so a downstream smooth_heaviside(f(x)) with a
    # SINGULAR f' in the saturated tail can't produce (tiny)·∞ = NaN under reverse-AD.
    kx >  _SMOOTH_SAT && return one(kx)
    kx < -_SMOOTH_SAT && return zero(kx)
    if kx >= zero(kx)
        return one(kx) / (one(kx) + exp(-kx))
    else
        e = exp(kx)
        return e / (one(e) + e)
    end
end

# --------------------------------------------------------------------------
# smooth_heaviside — stable sigmoid σ(k*x)
# --------------------------------------------------------------------------

"""
    smooth_heaviside(x; k=50.0)

Smooth approximation to the Heaviside step function.
- For `Float64`: exact `x >= 0 ? 1.0 : 0.0`.
- For `ForwardDiff.Dual`: stable sigmoid `1/(1 + exp(-k*x))`.
"""
function smooth_heaviside(x::Float64; k::Float64=50.0)
    (_use_smooth(Float64) || _smooth_f64()) || return x >= 0.0 ? 1.0 : 0.0
    return _stable_sigmoid(x, k)
end

function smooth_heaviside(x::T; k::Real=50) where {T<:Real}
    # TYPE-BASED (matches smooth_min/max, GPU-safe): exact step for non-AD reals (Float32 device
    # kernels → hard, no global read; and identical to the Float64-:auto path), smooth sigmoid for
    # ForwardDiff.Dual. The Float64 method above adds the host-only SMOOTH_MODE=:always override.
    _use_smooth(T) || return (x >= zero(T) ? one(T) : zero(T))
    return _stable_sigmoid(x, k)
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
    (_use_smooth(Float64) || _smooth_f64()) || return abs(x)
    return x * tanh(k * x)
end

function smooth_abs(x::T; k::Real=50) where {T<:Real}
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
    (_use_smooth(Float64) || _smooth_f64()) || return x >= 0.0 ? a : b
    h = _stable_sigmoid(x, k)
    return a * h + b * (1.0 - h)
end

function smooth_ifelse(x::T, a, b; k::Real=50) where {T<:Real}
    h = _stable_sigmoid(x, k)
    return a * h + b * (one(T) - h)
end
