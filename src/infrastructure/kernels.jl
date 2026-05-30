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

# --------------------------------------------------------------------------
# Clamp water-table depth to the bedrock interface (masked per-column).
# nbedrock[c] is a per-column bedrock layer index; the interface depth is
# zi[c, nbedrock[c] + nlevsno_off + 1]. Mask in-kernel (one thread per column).
# --------------------------------------------------------------------------
@kernel function _clamp_zwt_bedrock_kernel!(zwt_col, @Const(mask), @Const(nbedrock),
                                            @Const(zi), nlevsno_off::Int)
    c = @index(Global)
    @inbounds if mask[c]
        zi_bedrock = zi[c, nbedrock[c] + nlevsno_off + 1]
        if zwt_col[c] > zi_bedrock
            zwt_col[c] = zi_bedrock
        end
    end
end

"""
    clamp_zwt_to_bedrock!(zwt_col, mask, nbedrock, zi, nlevsno_off)

Clamp each active column's water-table depth to its bedrock interface.
Backend-agnostic; one thread per column.
"""
clamp_zwt_to_bedrock!(zwt_col, mask, nbedrock, zi, nlevsno_off::Int) =
    _launch!(_clamp_zwt_bedrock_kernel!, zwt_col, mask, nbedrock, zi, nlevsno_off)

# --------------------------------------------------------------------------
# Urban snow-cover fraction from snow depth (masked per-column, urban only).
# --------------------------------------------------------------------------
@kernel function _urban_frac_sno_kernel!(frac_sno_col, @Const(landunit),
                                         @Const(urbpoi), @Const(snow_depth_col))
    c = @index(Global)
    @inbounds if urbpoi[landunit[c]]
        sd = snow_depth_col[c]
        isfinite(sd) || (sd = zero(sd))
        frac_sno_col[c] = min(sd / 0.05, 1.0)
    end
end

"""
    update_urban_frac_sno!(frac_sno_col, landunit, urbpoi, snow_depth_col)

Set urban columns' snow-cover fraction from snow depth. Backend-agnostic.
"""
update_urban_frac_sno!(frac_sno_col, landunit, urbpoi, snow_depth_col) =
    _launch!(_urban_frac_sno_kernel!, frac_sno_col, landunit, urbpoi, snow_depth_col)
