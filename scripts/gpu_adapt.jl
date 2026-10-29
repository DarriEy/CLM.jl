# ==========================================================================
# gpu_adapt.jl — shared device adaptor for the gpu_validate_* harnesses
# (Float32 down-conversion on Float32-only backends; Float64 kept otherwise).
#
# Consolidates the `_F32` / `_F32S` adaptor pattern that was copy-pasted into
# ~23 harness scripts. `include` this and use `mf(MtlArray, x)` (or the curried
# `mf` once the backend array type is known) to move a Float64 host state tree to
# a device as Float32.
#
# Two adaptors, because structs disagree on how their scalar fields should move:
#   _F32  — arrays -> Float32, ::FT scalar fields -> Float32, BUT leave a few
#           structs' CONCRETE ::Float64 scalar fields as Float64 (they pin
#           ::Float64 in their ctor and must stay). Use for structs whose scalar
#           fields are genuinely Float64 constants (thresholds etc.).
#   _F32S — like _F32 but ALSO converts loose Float64 scalars to Float32. Use for
#           structs (CropData, CNVegState, TemperatureData) that pin ::FT scalars
#           which must match their Float32 arrays.
#
# Most harnesses only need _F32; route the few scalar-FT structs through _F32S.
# ==========================================================================

# --- _F32: arrays + ::FT scalars -> Float32, concrete-Float64 scalars kept ---
struct _F32 end
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{<:AbstractFloat}) = Float32.(x)
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{<:Integer}) = x
# BitVector / Bool arrays -> Vector{Bool} so the subsequent device adapt doesn't
# scalar-index packed bits.
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{Bool}) = collect(Bool, x)
CLM.Adapt.adapt_storage(::_F32, x::Float64) = Float32(x)

# Rebuild a struct adapting every field EXCEPT the named concrete-Float64 scalar
# field(s), preserved as-is (their ctor pins ::Float64).
function _adapt_keep_f64(a, x, keep::Tuple)
    T = typeof(x)
    vals = ntuple(fieldcount(T)) do i
        n = fieldname(T, i); v = getfield(x, i)
        n in keep ? v : CLM.Adapt.adapt(a, v)
    end
    return T.name.wrapper(vals...)
end
CLM.Adapt.adapt_structure(a::_F32, x::CLM.CanopyStateData) = _adapt_keep_f64(a, x, (:leaf_mr_vcm,))
CLM.Adapt.adapt_structure(a::_F32, x::CLM.SoilBiogeochemCarbonStateData) = _adapt_keep_f64(a, x, (:totvegcthresh,))
CLM.Adapt.adapt_structure(a::_F32, x::CLM.SoilBiogeochemNitrogenStateData) = _adapt_keep_f64(a, x, (:totvegcthresh,))

# DecompMIMICSState is parametric only on M<:AbstractMatrix (its 4 texture-array
# fields). Its ~40 scalar coeff fields pin ::Float64 in the auto-generated positional
# (outer) constructor that Adapt.@adapt_structure invokes — so the default _F32/_F32S
# rule (Float64 scalar -> Float32) makes that ctor un-callable (MethodError, seen once
# #294 wired decomp_mimics_state into the instance tree). Mirror the MetalF32 path used
# by gpu_validate_decompmimics_e2e.jl: adapt ONLY the array fields, pass every scalar
# (Int and Float64) through unchanged. The device Adapt.@adapt_structure step then moves
# the (now Float32) texture arrays to the device and leaves the scalars alone.
function _adapt_arrays_only(a, x)
    T = typeof(x)
    vals = ntuple(fieldcount(T)) do i
        v = getfield(x, i)
        v isa AbstractArray ? CLM.Adapt.adapt(a, v) : v
    end
    return T.name.wrapper(vals...)
end
CLM.Adapt.adapt_structure(a::_F32, x::CLM.DecompMIMICSState) = _adapt_arrays_only(a, x)

# --- _F32S: arrays + ALL Float64 scalars -> Float32 (no kept-Float64 structs) ---
struct _F32S end
CLM.Adapt.adapt_storage(::_F32S, x::AbstractArray{<:AbstractFloat}) = Float32.(x)
CLM.Adapt.adapt_storage(::_F32S, x::AbstractArray{<:Integer}) = x
CLM.Adapt.adapt_storage(::_F32S, x::AbstractArray{Bool}) = collect(Bool, x)
CLM.Adapt.adapt_storage(::_F32S, x::Float64) = Float32(x)
# Same DecompMIMICSState guard as under _F32: its parametric outer ctor pins ::Float64
# on the scalar coeff fields, so keep every scalar and adapt only the texture arrays.
CLM.Adapt.adapt_structure(a::_F32S, x::CLM.DecompMIMICSState) = _adapt_arrays_only(a, x)

# --- _HostPrep: precision-preserving host pre-pass (Float64 backends) ---------
# Only unpacks BitArray/Bool masks to Vector{Bool} (device adapt cannot take packed
# bits); every float array and scalar keeps its Float64 eltype.
struct _HostPrep end
CLM.Adapt.adapt_storage(::_HostPrep, x::AbstractArray{Bool}) = collect(Bool, x)
CLM.Adapt.adapt_storage(::_HostPrep, x::AbstractArray) = x

# Down-convert to Float32 only when the detected backend is Float32-only (Metal).
# On a Float64 backend (CUDA/AMDGPU/oneAPI) `mf` used to demote the whole state tree
# to Float32 regardless, so those harnesses printed "working precision: Float64"
# while comparing a Float32 device run against the Float64 CPU reference — the same
# silent demotion #263 removed from `cu()`. GPU_ADAPT_F32=1 forces the old Float32
# path on any backend (e.g. to reproduce a Metal-precision result on CUDA).
function _mf_f32()
    get(ENV, "GPU_ADAPT_F32", "0") == "1" && return true
    b = isdefined(@__MODULE__, :detect_backend) ? detect_backend() : nothing
    return b === nothing || b[3] === Float32
end

# Move a host state tree to `devarray` (e.g. Metal.MtlArray) at the backend's
# working precision: Float32 via _F32 on Metal, Float64 kept on CUDA/AMDGPU.
mf(devarray, x) = _mf_f32() ? CLM.Adapt.adapt(devarray, CLM.Adapt.adapt(_F32(), x)) :
                              CLM.Adapt.adapt(devarray, CLM.Adapt.adapt(_HostPrep(), x))
# Variant routing scalar-FT structs through _F32S.
mfs(devarray, x) = _mf_f32() ? CLM.Adapt.adapt(devarray, CLM.Adapt.adapt(_F32S(), x)) :
                               CLM.Adapt.adapt(devarray, CLM.Adapt.adapt(_HostPrep(), x))
