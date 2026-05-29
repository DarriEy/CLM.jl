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
