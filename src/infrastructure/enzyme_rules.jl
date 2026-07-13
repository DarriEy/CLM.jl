# ==========================================================================
# Enzyme reverse-mode AD rules for CLM.
#
# The headline rule is for the banded linear solve `band_solve!`, which on
# Float64 uses LAPACK (gbtrf!/gbtrs!) inside a try/catch — neither of which
# Enzyme can differentiate through. A custom rule lets Enzyme skip the primal
# internals entirely and supply the exact linear-solve adjoint instead.
#
# For A x = rhs (A in LAPACK band storage `ab`, solution written to u_slice):
#   reverse:  solve Aᵀ λ = x̄ ;  rhs̄ += λ ;  ā(band) += -λ xᵀ
# Validated against analytic A⁻ᵀw and finite differences before integration.
#
# The same adjoint is supplied for the single-column Thomas solve
# `tridiagonal_solve!` — see the block at the bottom of this file for why that
# one is a CORRECTNESS requirement, not just a performance one.
# ==========================================================================

import Enzyme
import Enzyme.EnzymeRules: augmented_primal, reverse, AugmentedReturn, inactive_type
using Enzyme: Const, Active, Duplicated

# --------------------------------------------------------------------------
# Constant parameter / control containers.
#
# These module-level structs hold physical constants and control flags that are
# never differentiation targets (calibration parameters flow through inst state
# — inst.overrides / inst.column — not through these globals). Declaring their
# types inactive tells Enzyme they carry no derivative info, which resolves the
# "mutable global param struct" class of reverse-mode blockers.
# --------------------------------------------------------------------------
for T in (:BareGroundFluxesParamsData, :CanopyFluxesParamsData, :CanopyFluxesControl,
          :CanopyHydrologyParamsData, :CanopyHydrologyControl, :UrbanFluxesParamsData,
          :UrbanControl, :SoilHydrologyParams, :SnowHydrologyParams,
          :SurfaceResistanceControl, :SoilMoistStressControl, :PhotoParamsData)
    @eval inactive_type(::Type{<:$T}) = true
end

function augmented_primal(config, func::Const{typeof(band_solve!)}, ::Type{RT},
                          u_slice::Duplicated, ab::Duplicated, rhs::Duplicated,
                          kl::Const, ku::Const, n::Const) where {RT}
    func.val(u_slice.val, ab.val, rhs.val, kl.val, ku.val, n.val)
    # band_solve! preserves ab.val (it factorizes a copy), so we can save the
    # original A and the computed solution x for the reverse pass.
    tape = (copy(ab.val), copy(u_slice.val))
    return AugmentedReturn(nothing, nothing, tape)
end

function reverse(config, func::Const{typeof(band_solve!)}, ::Type{RT}, tape,
                 u_slice::Duplicated, ab::Duplicated, rhs::Duplicated,
                 kl::Const, ku::Const, n::Const) where {RT}
    ab_orig, x = tape
    klv, kuv, nv = kl.val, ku.val, n.val
    x̄ = u_slice.dval

    # Solve Aᵀ λ = x̄ (transpose solve via gbtrs! 'T').
    abw = copy(ab_orig)
    (ab_lu, ipiv) = LinearAlgebra.LAPACK.gbtrf!(klv, kuv, nv, abw)
    λmat = reshape(copy(x̄), nv, 1)
    LinearAlgebra.LAPACK.gbtrs!('T', klv, kuv, nv, ab_lu, ipiv, λmat)
    λ = vec(λmat)

    # rhs̄ += λ
    rhs.dval .+= λ
    # ā += -λ xᵀ within the band, mapped into LAPACK band storage rows.
    @inbounds for j in 1:nv
        for i in max(1, j - kuv):min(nv, j + klv)
            ab.dval[klv + kuv + 1 + i - j, j] += -λ[i] * x[j]
        end
    end
    # The output adjoint has been fully propagated; consume it.
    x̄ .= 0
    return (nothing, nothing, nothing, nothing, nothing, nothing)
end

# --------------------------------------------------------------------------
# tridiagonal_solve! — single-column Thomas solve.
#
# WHY THIS RULE EXISTS (correctness, not just speed):
#
# `tridiagonal_solve!` heap-allocates its Thomas workspace (`cp`, `dp`) INSIDE
# the callee. On Julia 1.10, Enzyme (0.13.x) does not zero-initialize the
# SHADOW of an array allocated inside a differentiated callee, so the reverse
# pass accumulates into — and reads back — uninitialized memory. Measured on
# Julia 1.10.11 / Enzyme 0.13.183, the un-ruled reverse gradient of the 5×5
# smoke system was wrong on EVERY call (-0.5151 vs the true -0.5547, i.e. the
# recycled contents of the previous call's shadow) and NON-FINITE on ~2-5% of
# calls (NaN / -Inf / -2.1e156 — whatever the allocator last left there). That
# is the source of the intermittent CI failure of the "tridiagonal solver
# reverse-mode" smoke test on the 1.10 minimum.
#
# Proof of the root cause (Julia 1.10): hoisting the identical allocation into
# the caller's frame, or handing the callee a workspace whose shadow we zero
# ourselves, both yield the exact gradient (-0.5547337278106508 vs the
# ForwardDiff truth -0.5547337278106509); re-running with a deliberately
# garbage-filled workspace shadow reproduces the wrong value on demand. So it
# is the shadow of the in-callee allocation, not the primal (every element of
# cp/dp is written before it is read) and not the cache (Enzyme.API.zcache!(true)
# does not help). Julia >= 1.11 allocates arrays through GenericMemory and is
# unaffected. `tridiagonal_multi!` is also unaffected: its workspace comes from
# `similar(u, ...)`, whose shadow Enzyme derives from u's live shadow.
#
# Rather than gate the test off on 1.10 (which would leave real calibration
# gradients silently wrong on the minimum supported version), we give Enzyme the
# exact analytic adjoint and never let it walk the workspace at all. This is
# version-independent, allocator-independent, and O(n) instead of taping the
# whole forward sweep.
#
# For the system over levels jtop..nlevs
#     a[j]*u[j-1] + b[j]*u[j] + c[j]*u[j+1] = r[j]      (A u = r)
# the reverse pass is
#     λ = A⁻ᵀ ū        (Aᵀ is tridiagonal: sub = c[j-1], diag = b[j], super = a[j+1])
#     r̄[j] += λ[j]
#     b̄[j] -= λ[j]*u[j] ;  ā[j] -= λ[j]*u[j-1] ;  c̄[j] -= λ[j]*u[j+1]     (Ā = -λ uᵀ)
# and ū[jtop:nlevs] is consumed (u is fully overwritten by the call).
# --------------------------------------------------------------------------
const _TriAnn = Union{Const{<:AbstractVector{<:Real}}, Duplicated{<:AbstractVector{<:Real}}}

function augmented_primal(config, func::Const{typeof(tridiagonal_solve!)}, ::Type{RT},
                          u::Duplicated{<:AbstractVector{<:Real}},
                          a::_TriAnn, b::_TriAnn, c::_TriAnn, r::_TriAnn,
                          jtop::Const, nlevs::Const) where {RT}
    func.val(u.val, a.val, b.val, c.val, r.val, jtop.val, nlevs.val)
    # Save the solution and the matrix bands: the caller is free to overwrite
    # a/b/c/r between the forward and reverse sweeps, so we cannot alias them.
    tape = (copy(u.val), copy(a.val), copy(b.val), copy(c.val))
    return AugmentedReturn(nothing, nothing, tape)
end

function reverse(config, func::Const{typeof(tridiagonal_solve!)}, ::Type{RT}, tape,
                 u::Duplicated{<:AbstractVector{<:Real}},
                 a::_TriAnn, b::_TriAnn, c::_TriAnn, r::_TriAnn,
                 jtop::Const, nlevs::Const) where {RT}
    x, a0, b0, c0 = tape
    j0 = jtop.val
    n  = nlevs.val
    ū  = u.dval
    T  = eltype(x)

    # Transpose system Aᵀ λ = ū, solved with the same Thomas routine.
    aT = zeros(T, n); bT = zeros(T, n); cT = zeros(T, n)
    rT = zeros(T, n); λ  = zeros(T, n)
    @inbounds for j in j0:n
        bT[j] = b0[j]
        rT[j] = ū[j]
        if j > j0; aT[j] = c0[j-1]; end
        if j < n;  cT[j] = a0[j+1]; end
    end
    tridiagonal_solve!(λ, aT, bT, cT, rT, j0, n)

    # r̄ += λ ;  Ā = -λ xᵀ restricted to the three bands.
    if !(r isa Const)
        @inbounds for j in j0:n; r.dval[j] += λ[j]; end
    end
    if !(b isa Const)
        @inbounds for j in j0:n; b.dval[j] -= λ[j] * x[j]; end
    end
    if !(a isa Const)
        @inbounds for j in (j0+1):n; a.dval[j] -= λ[j] * x[j-1]; end
    end
    if !(c isa Const)
        @inbounds for j in j0:(n-1); c.dval[j] -= λ[j] * x[j+1]; end
    end

    # u[jtop:nlevs] is overwritten (assigned, not accumulated) by the call: its
    # incoming adjoint has now been fully propagated, so consume it.
    @inbounds for j in j0:n; ū[j] = zero(eltype(ū)); end
    return (nothing, nothing, nothing, nothing, nothing, nothing, nothing)
end

# If `u` itself is Const the call is inactive (nothing downstream reads a
# derivative of the solution), so there is nothing to propagate — but we still
# must supply the rule, otherwise Enzyme falls back to walking the workspace.
function augmented_primal(config, func::Const{typeof(tridiagonal_solve!)}, ::Type{RT},
                          u::Const{<:AbstractVector{<:Real}},
                          a::_TriAnn, b::_TriAnn, c::_TriAnn, r::_TriAnn,
                          jtop::Const, nlevs::Const) where {RT}
    func.val(u.val, a.val, b.val, c.val, r.val, jtop.val, nlevs.val)
    return AugmentedReturn(nothing, nothing, nothing)
end

function reverse(config, func::Const{typeof(tridiagonal_solve!)}, ::Type{RT}, tape,
                 u::Const{<:AbstractVector{<:Real}},
                 a::_TriAnn, b::_TriAnn, c::_TriAnn, r::_TriAnn,
                 jtop::Const, nlevs::Const) where {RT}
    return (nothing, nothing, nothing, nothing, nothing, nothing, nothing)
end
