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

# --------------------------------------------------------------------------
# Atomic integer width for scatter counters
# --------------------------------------------------------------------------
# Some GPU backends — notably Apple Metal — have NO 64-bit atomics, so an atomic
# `_scatter_add!` into an `Int64` counter fails to compile (gpu_gc_pool_alloc /
# throw_methoderror). `atomic_int_type(ref)` picks the integer eltype for an
# atomic counter from the backend of `ref`: Int32 where 64-bit atomics are
# unavailable, Int (Int64) otherwise. AUTODETECTED by default; override with
# `set_atomic_int_width!(:auto | :i32 | :i64)`.
const _ATOMIC_INT_WIDTH = Ref(:auto)

"""
    set_atomic_int_width!(mode::Symbol)

Control the integer width used for atomic scatter counters (e.g. `numtrees`).
`:auto` (default) picks Int32 on backends without 64-bit atomics (Metal) and Int
elsewhere (CPU/CUDA/ROCm); `:i32`/`:i64` force a width. Returns `mode`.
"""
set_atomic_int_width!(mode::Symbol) =
    (mode in (:auto, :i32, :i64) || throw(ArgumentError("mode must be :auto, :i32 or :i64"));
     _ATOMIC_INT_WIDTH[] = mode; mode)

# Autodetection: KA.CPU and CUDA/ROCm support 64-bit atomics; Metal does not.
# Detected by backend type name so no hard dependency on the GPU packages.
_backend_has_64bit_atomics(::KA.CPU) = true
_backend_has_64bit_atomics(b::KA.Backend) = !occursin("Metal", string(nameof(typeof(b))))

atomic_int_type(ref::AbstractArray) = atomic_int_type(_kernel_backend(ref))
function atomic_int_type(backend::KA.Backend)
    mode = _ATOMIC_INT_WIDTH[]
    mode === :i32 && return Int32
    mode === :i64 && return Int64
    return _backend_has_64bit_atomics(backend) ? Int : Int32
end

# Sync-fusion: inside a `with_deferred_gpu_sync` region the per-op GPU synchronize is
# skipped — kernels enqueued on the same GPU queue execute in order, so intermediate host
# syncs are unnecessary; one sync at region end (before the host reads results) suffices.
# The CPU (KA) backend is NEVER deferred (its kernels may run async on threads, so the next
# op could race a not-yet-finished one) — so the host path keeps its per-op sync unchanged.
const _DEFER_GPU_SYNC = Ref(false)

# Move a host constant array `v` onto the same backend + precision as the state prototype
# `proto`, so a kernel launched on the state's backend isn't handed mixed host/device args.
# No-op (returns `v` unchanged) when the state is on the host — the host path stays
# byte-identical. Used where the driver supplies physical constants (zsoi/dz/…) as host
# vectors but the kernel's other operands live on the device.
@inline _to_backend_like(::Array, ::Type, v::AbstractArray) = v
@inline _to_backend_like(proto, ::Type{FT}, v::AbstractArray) where {FT} =
    typeof(proto).name.wrapper(FT.(v))
# Integer topology/index arrays (cascade pools, patch->column maps, …): move to the
# backend but PRESERVE the element type — casting indices to FT would corrupt them.
@inline _to_backend_like(::Array, ::Type, v::AbstractArray{<:Integer}) = v
@inline _to_backend_like(proto, ::Type{FT}, v::AbstractArray{<:Integer}) where {FT} =
    typeof(proto).name.wrapper(v)

@inline function _launch!(kernel, out, args...; ndrange = length(out))
    (ndrange isa Integer ? ndrange == 0 : prod(ndrange) == 0) && return out
    backend = _kernel_backend(out)
    kernel(backend)(out, args...; ndrange = ndrange)
    (_DEFER_GPU_SYNC[] && !(backend isa KA.CPU)) || KA.synchronize(backend)
    return out
end

# Defer per-op GPU syncs across `f` (one sync at the end). No-op deferral for the CPU
# backend / a nothing backend, so the host path is unchanged. `backend` is the region's
# device (derived from a workspace prototype); the final sync targets it.
@inline function with_deferred_gpu_sync(f, backend)
    (backend === nothing || backend isa KA.CPU) && return f()
    prev = _DEFER_GPU_SYNC[]; _DEFER_GPU_SYNC[] = true
    try
        return f()
    finally
        _DEFER_GPU_SYNC[] = prev
        KA.synchronize(backend)
    end
end

# Scatter-add `arr[i] += x` for kernels with a many-to-one (e.g. patch→column) write.
# Atomic for hardware-real element types — the parallel GPU path (Atomix has a Metal
# extension). Plain `+=` for ForwardDiff.Dual (and any non-hardware-atomic eltype):
# such arrays only occur on the sequential KA CPU backend (forward-mode AD), where there
# is no race and atomic-modify is not even defined for Dual.
@inline _scatter_add!(arr::AbstractArray{<:Union{AbstractFloat,Integer}}, i, x) =
    (Atomix.@atomic arr[i] += x; nothing)
@inline _scatter_add!(arr, i, x) = (@inbounds arr[i] += x; nothing)
# 2D variant (e.g. per-(column, layer) scatter).
@inline _scatter_add!(arr::AbstractArray{<:Union{AbstractFloat,Integer}}, i, j, x) =
    (Atomix.@atomic arr[i, j] += x; nothing)
@inline _scatter_add!(arr, i, j, x) = (@inbounds arr[i, j] += x; nothing)
# 3D variant (e.g. per-(column, layer, litter-pool) scatter).
@inline _scatter_add!(arr::AbstractArray{<:Union{AbstractFloat,Integer}}, i, j, k, x) =
    (Atomix.@atomic arr[i, j, k] += x; nothing)
@inline _scatter_add!(arr, i, j, k, x) = (@inbounds arr[i, j, k] += x; nothing)

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
