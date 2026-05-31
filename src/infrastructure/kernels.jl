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
# Backend of the output array. BitArray masks have no KernelAbstractions backend
# (packed bits aren't a device array type), so kernels writing a BitVector mask run
# on the CPU. Moving such masks to GPU requires converting them to Vector{Bool}
# (a documented follow-up); until then this keeps the CPU path correct.
_kernel_backend(out) = KA.get_backend(out)
_kernel_backend(::BitArray) = KA.CPU()

@inline function _launch!(kernel, out, args...; ndrange = length(out))
    (ndrange isa Integer ? ndrange == 0 : prod(ndrange) == 0) && return out
    backend = _kernel_backend(out)
    kernel(backend)(out, args...; ndrange = ndrange)
    KA.synchronize(backend)
    return out
end

# Scatter-add `arr[i] += x` for kernels with a many-to-one (e.g. patch→column) write.
# Atomic for hardware-real element types — the parallel GPU path (Atomix has a Metal
# extension). Plain `+=` for ForwardDiff.Dual (and any non-hardware-atomic eltype):
# such arrays only occur on the sequential KA CPU backend (forward-mode AD), where there
# is no race and atomic-modify is not even defined for Dual.
@inline _scatter_add!(arr::AbstractArray{<:Union{AbstractFloat,Integer}}, i, x) =
    (Atomix.@atomic arr[i] += x; nothing)
@inline _scatter_add!(arr, i, x) = (@inbounds arr[i] += x; nothing)

# --------------------------------------------------------------------------
# forc_q: column specific humidity from vapor pressure
#   forc_q_col[c] = 0.622 * vp / max(pbot - 0.378*vp, 1)  (vp from the gridcell)
# --------------------------------------------------------------------------
@kernel function _forc_q_kernel!(forc_q_col, @Const(gridcell), @Const(vp_grc), @Const(pbot_col))
    c = @index(Global)
    @inbounds begin
        # Literals are converted to the working element type (`oftype`/`one`) so the
        # kernel carries no Float64 on a Float32-only backend (Metal). On Float64 this
        # is byte-identical (oftype(::Float64, 0.622) === 0.622).
        vp = vp_grc[gridcell[c]]
        pb = pbot_col[c]
        forc_q_col[c] = oftype(vp, 0.622) * vp / max(pb - oftype(vp, 0.378) * vp, one(vp))
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
        frac_sno_col[c] = min(sd / oftype(sd, 0.05), one(sd))
    end
end

"""
    update_urban_frac_sno!(frac_sno_col, landunit, urbpoi, snow_depth_col)

Set urban columns' snow-cover fraction from snow depth. Backend-agnostic.
"""
update_urban_frac_sno!(frac_sno_col, landunit, urbpoi, snow_depth_col) =
    _launch!(_urban_frac_sno_kernel!, frac_sno_col, landunit, urbpoi, snow_depth_col)

# --------------------------------------------------------------------------
# Volumetric liquid water content (per column AND soil layer — a 2D kernel).
#   liqvol[c,j+joff] = clamp(h2osoi_liq/(dz*DENH2O), 0, eff_porosity[c,j]) if dz>0
# Soil/snow arrays carry an nlevsno offset; eff_porosity is layer-indexed.
# --------------------------------------------------------------------------
@kernel function _h2osoi_liqvol_kernel!(liqvol_col, @Const(mask), @Const(dz),
                                        @Const(h2osoi_liq), @Const(eff_porosity),
                                        joff::Int, denh2o)
    c, j = @index(Global, NTuple)
    @inbounds if mask[c]
        dz_cj = dz[c, j + joff]
        if dz_cj > 0
            lv = h2osoi_liq[c, j + joff] / (dz_cj * denh2o)
            liqvol_col[c, j + joff] = min(max(lv, zero(lv)), eff_porosity[c, j])
        else
            liqvol_col[c, j + joff] = zero(dz_cj)
        end
    end
end

"""
    compute_h2osoi_liqvol!(liqvol_col, mask, dz, h2osoi_liq, eff_porosity, joff, nlevgrnd; denh2o)

Volumetric liquid water for root moisture stress, over all (column, soil layer)
pairs — a 2D KernelAbstractions kernel. Backend-agnostic.
"""
function compute_h2osoi_liqvol!(liqvol_col, mask, dz, h2osoi_liq, eff_porosity,
                                joff::Int, nlevgrnd::Int; denh2o)
    # Pass the scalar constant at the output's precision so no Float64 reaches a
    # Float32-only backend (Metal).
    denh2o_ft = convert(eltype(liqvol_col), denh2o)
    _launch!(_h2osoi_liqvol_kernel!, liqvol_col, mask, dz, h2osoi_liq, eff_porosity,
             joff, denh2o_ft; ndrange = (length(mask), nlevgrnd))
end
