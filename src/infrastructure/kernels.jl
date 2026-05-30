# ==========================================================================
# KernelAbstractions kernels for per-column / per-patch physics.
#
# Each kernel is backend-agnostic: it runs as a plain loop on CPU Arrays and as a
# GPU kernel when given device arrays (CuArray/ROCArray). The backend is taken
# from the output array, so the same driver code path runs on CPU or GPU
# depending only on where the state lives (see Adapt-based device movement).
#
# This is the Phase-4 kernelization entry point; per-column/patch driver loops are
# migrated here incrementally, each validated against its scalar version.
# ==========================================================================

import KernelAbstractions as KA
using KernelAbstractions: @kernel, @index, @Const

"""
    _launch!(kernel, out, args...; ndrange=length(out))

Run a KernelAbstractions `kernel` over `ndrange` on the backend of `out`
(one thread per index). No-op for empty `out`.
"""
@inline function _launch!(kernel, out, args...; ndrange::Int = length(out))
    ndrange == 0 && return out
    backend = KA.get_backend(out)
    kernel(backend)(out, args...; ndrange = ndrange)
    KA.synchronize(backend)
    return out
end

# --------------------------------------------------------------------------
# forc_q: column specific humidity from vapor pressure
#   forc_q_col[c] = 0.622 * vp / max(pbot - 0.378*vp, 1)  (vp from the gridcell)
# --------------------------------------------------------------------------
@kernel function _forc_q_kernel!(forc_q_col, @Const(gridcell), @Const(vp_grc), @Const(pbot_col))
    c = @index(Global)
    @inbounds begin
        vp = vp_grc[gridcell[c]]
        forc_q_col[c] = 0.622 * vp / max(pbot_col[c] - 0.378 * vp, 1.0)
    end
end

"""
    compute_forc_q!(forc_q_col, gridcell, vp_grc, pbot_col)

Column specific humidity from vapor pressure, one thread per column. Backend-
agnostic (CPU loop or GPU). Replaces the inline per-column loop in clm_drv_core!.
Assumes columns are indexed 1:nc (single clump, begc==1).
"""
compute_forc_q!(forc_q_col, gridcell, vp_grc, pbot_col) =
    _launch!(_forc_q_kernel!, forc_q_col, gridcell, vp_grc, pbot_col)
