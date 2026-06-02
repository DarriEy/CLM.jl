# ==========================================================================
# Ported from: src/biogeochem/CNPrecisionControlMod.F90
# Controls on very low values in critical CN vegetation state variables.
# Truncates small C/N pool values on patches to zero for numerical stability
# and accumulates the truncated mass into truncation sink terms (ctrunc/ntrunc).
#
# Public functions:
#   cn_precision_control!           -- Apply precision control to veg CN states
#   truncate_c_and_n_states!        -- Truncate paired C and N state vectors
#   truncate_c_states!              -- Truncate a carbon-only state vector
#   truncate_n_states!              -- Truncate a nitrogen-only state vector
#   truncate_additional!            -- Truncate an additional state (e.g. isotope)
#                                      using a pre-computed truncation filter
#
# GPU kernelization note (Phase B):
#   Each truncate_* helper's per-filter loop is a fully-independent per-patch
#   loop (each thread fp owns its own patch index p = filter[fp]; the filter
#   lists distinct patches, so writes never collide). The arithmetic loop is
#   moved into a KernelAbstractions kernel (CPU + Metal) operating over the
#   filter range (ndrange = length(filter)). Every Float64 literal inside a
#   kernel is converted to the working element type T so the kernel carries no
#   Float64 on a Float32-only backend (Metal); on Float64, T(x) === x so the
#   math is byte-identical.
#
#   The two host-side responsibilities that cannot live inside a GPU kernel are
#   kept on the host and are byte-identical to the original control flow:
#     (1) the critically-negative `error()` aborts (a host pre-scan over the
#         filter raises on exactly the same inputs as before), and
#     (2) the dynamically-built `filter_truncatep` return value: the kernel
#         writes a per-fp boolean `did_truncate` flag, from which the host
#         builds `filter_truncatep` in filter order — the same order (and same
#         contents) the original `push!` produced.
# ==========================================================================

# @kernel / @index / @Const / _launch! are module-wide (infrastructure/kernels.jl).

# ---------------------------------------------------------------------------
# Module-level constants (matching Fortran module data in CNPrecisionControlMod)
# ---------------------------------------------------------------------------

"""Critical carbon state value for truncation (gC/m2)."""
const CN_CCRIT_DEFAULT = 1.0e-8

"""Critical negative carbon state value for abort (gC/m2)."""
const CN_CNEGCRIT_DEFAULT = -6.0e+1

"""Critical nitrogen state value for truncation (gN/m2)."""
const CN_NCRIT_DEFAULT = 1.0e-8

"""Critical negative nitrogen state value for abort (gN/m2)."""
const CN_NNEGCRIT_DEFAULT = -7.0e+0

"""Minimum nitrogen value used when calculating CN ratio (gN/m2)."""
const CN_N_MIN = 1.0e-9

# ---------------------------------------------------------------------------
# Kernels — one thread per filter entry fp; writes only its own patch index p.
# ---------------------------------------------------------------------------

# truncate_c_and_n_states! arithmetic. `did_truncate[fp]` records whether
# patch p = filter[fp] was truncated, so the host can rebuild filter_truncatep
# in the original push! order. Literals are converted to the working precision
# T = typeof(ccrit). On Float64 every T(x) === x (byte-identical).
@kernel function _cnprec_cn_kernel!(carbon_patch, nitrogen_patch, pc, pn,
        did_truncate,
        @Const(filter_bgc_vegp), @Const(patch_itype),
        ccrit, cnegcrit, ncrit, nnegcrit,
        croponly::Bool, use_nguardrail::Bool, use_matrixcn::Bool,
        allowneg::Bool, nc3crop::Int)
    fp = @index(Global)
    @inbounds begin
        p = filter_bgc_vegp[fp]
        T = typeof(ccrit)
        did_truncate[fp] = false

        skip = croponly && patch_itype[p] < nc3crop
        if !skip
            # The host pre-scan already raised on critically-negative states
            # when !allowneg, so here we only perform the truncation math.
            if use_matrixcn
                if (carbon_patch[p] < ccrit && carbon_patch[p] > -ccrit * T(1.0e+6)) ||
                   (use_nguardrail && nitrogen_patch[p] < ncrit &&
                    nitrogen_patch[p] > -ncrit * T(1.0e+6))
                    did_truncate[fp] = true
                    pc[p] += carbon_patch[p]
                    carbon_patch[p] = zero(eltype(carbon_patch))
                    pn[p] += nitrogen_patch[p]
                    nitrogen_patch[p] = zero(eltype(nitrogen_patch))
                end
            else
                if abs(carbon_patch[p]) < ccrit ||
                   (use_nguardrail && abs(nitrogen_patch[p]) < ncrit)
                    did_truncate[fp] = true
                    pc[p] += carbon_patch[p]
                    carbon_patch[p] = zero(eltype(carbon_patch))
                    pn[p] += nitrogen_patch[p]
                    nitrogen_patch[p] = zero(eltype(nitrogen_patch))
                end
            end
        end
    end
end

# truncate_c_states! arithmetic (carbon only).
@kernel function _cnprec_c_kernel!(carbon_patch, pc, did_truncate,
        @Const(filter_bgc_vegp), @Const(patch_itype),
        ccrit, croponly::Bool, nc3crop::Int)
    fp = @index(Global)
    @inbounds begin
        p = filter_bgc_vegp[fp]
        did_truncate[fp] = false
        skip = croponly && patch_itype[p] < nc3crop
        if !skip
            # Host pre-scan already raised on critically-negative C (!allowneg).
            if abs(carbon_patch[p]) < ccrit
                did_truncate[fp] = true
                pc[p] += carbon_patch[p]
                carbon_patch[p] = zero(eltype(carbon_patch))
            end
        end
    end
end

# truncate_n_states! arithmetic (nitrogen only; no filter_truncatep returned).
@kernel function _cnprec_n_kernel!(nitrogen_patch, pn,
        @Const(filter_bgc_vegp), ncrit, nnegcrit)
    fp = @index(Global)
    @inbounds begin
        p = filter_bgc_vegp[fp]
        if nitrogen_patch[p] < nnegcrit
            # Fortran source has the endrun call commented out here. No abort.
        elseif abs(nitrogen_patch[p]) < ncrit
            pn[p] += nitrogen_patch[p]
            nitrogen_patch[p] = zero(eltype(nitrogen_patch))
        end
    end
end

# truncate_additional! arithmetic: zero an extra (isotope) state on a
# pre-computed truncation filter. One thread per filter entry fp.
@kernel function _cnprec_additional_kernel!(state_patch, truncation_patch,
        @Const(filter_truncatep))
    fp = @index(Global)
    @inbounds begin
        p = filter_truncatep[fp]
        truncation_patch[p] += state_patch[p]
        state_patch[p] = zero(eltype(state_patch))
    end
end

# Final accumulation of pc/pn (and isotopes) into the ctrunc/ntrunc sinks.
@kernel function _cnprec_accum_kernel!(ctrunc_patch, ntrunc_patch,
        @Const(pc), @Const(pn), @Const(filter_bgc_vegp))
    fp = @index(Global)
    @inbounds begin
        p = filter_bgc_vegp[fp]
        ctrunc_patch[p] += pc[p]
        ntrunc_patch[p] += pn[p]
    end
end

@kernel function _cnprec_accum_iso_kernel!(iso_ctrunc_patch,
        @Const(piso), @Const(filter_bgc_vegp))
    fp = @index(Global)
    @inbounds begin
        p = filter_bgc_vegp[fp]
        iso_ctrunc_patch[p] += piso[p]
    end
end

# Move a host index vector `f` onto the backend of `ref` (a state array). On the
# CPU path `ref` is a plain Array and `f` is returned unchanged; on a device
# backend `Adapt.adapt(typeof(ref), f)` produces the matching device vector so
# the kernel can read the filter without host scalar-indexing the device array.
@inline _cnprec_match_backend(ref::Array, f) = f
# `f` is a host index Vector{Int}; move it onto `ref`'s backend while PRESERVING
# its integer eltype. Adapting via `typeof(ref)` (a float array type) would coerce
# the indices to Float, producing a float "index" array that can't index a device
# array (invalid GPU IR). `similar(ref, eltype(f), …)` keeps Int on the backend.
@inline function _cnprec_match_backend(ref, f)
    d = similar(ref, eltype(f), length(f))
    copyto!(d, f)
    return d
end

# Build filter_truncatep (in filter order) from the per-fp did_truncate flags.
# This is host bookkeeping (a dynamic-length result); the truncation decision
# itself was made in the kernel. Preserves the exact order/contents the
# original `push!` produced. `did_truncate` / `filter_bgc_vegp` may be device
# arrays, so they are materialized to host once (no per-element device indexing).
@inline function _build_truncate_filter(filter_bgc_vegp, did_truncate)
    fh = Array(filter_bgc_vegp)
    dh = Array(did_truncate)
    filter_truncatep = Int[]
    @inbounds for fp in eachindex(fh)
        if dh[fp]
            push!(filter_truncatep, Int(fh[fp]))
        end
    end
    return filter_truncatep
end

# ---------------------------------------------------------------------------
# truncate_c_and_n_states!
# Ported from TruncateCandNStates
# ---------------------------------------------------------------------------

"""
    truncate_c_and_n_states!(carbon_patch, nitrogen_patch, pc, pn,
        filter_bgc_vegp, patch_itype;
        ccrit, cnegcrit, ncrit, nnegcrit,
        use_nguardrail, use_matrixcn,
        croponly, allowneg, nc3crop) -> (num_truncatep, filter_truncatep)

Truncate paired carbon and nitrogen state vectors. If a paired C/N value is
too small, zero both and accumulate the removed mass into the truncation
accumulators `pc` and `pn`.

Returns a tuple `(num_truncatep, filter_truncatep)` listing patch indices
where truncation occurred, for use by `truncate_additional!`.

Ported from `TruncateCandNStates` in `CNPrecisionControlMod.F90`.
"""
function truncate_c_and_n_states!(
        carbon_patch::AbstractVector,
        nitrogen_patch::AbstractVector,
        pc::AbstractVector,
        pn::AbstractVector,
        filter_bgc_vegp::AbstractVector{<:Integer},
        patch_itype::AbstractVector{<:Integer};
        ccrit::Real = CN_CCRIT_DEFAULT,
        cnegcrit::Real = CN_CNEGCRIT_DEFAULT,
        ncrit::Real = CN_NCRIT_DEFAULT,
        nnegcrit::Real = CN_NNEGCRIT_DEFAULT,
        use_nguardrail::Bool = false,
        use_matrixcn::Bool = false,
        croponly::Bool = false,
        allowneg::Bool = false,
        nc3crop::Int = 15)

    isempty(filter_bgc_vegp) && return (0, Int[])

    # Host pre-scan: raise on critically-negative states exactly as the original
    # in-loop check did (same inputs -> same error). Done before any mutation.
    # Arrays are materialized to host once (no per-element device indexing).
    if !allowneg
        cH = Array(carbon_patch); nH = Array(nitrogen_patch)
        fH = Array(filter_bgc_vegp); iH = Array(patch_itype)
        @inbounds for fp in eachindex(fH)
            p = fH[fp]
            if croponly && iH[p] < nc3crop
                continue
            end
            if cH[p] < cnegcrit || nH[p] < nnegcrit
                error("carbon or nitrogen state critically negative: " *
                      "C=$(cH[p]), N=$(nH[p]), " *
                      "limits: cnegcrit=$cnegcrit, nnegcrit=$nnegcrit")
            end
        end
    end

    did_truncate = fill!(similar(carbon_patch, Bool, length(filter_bgc_vegp)), false)
    # Convert the crit thresholds to the state working precision so the kernel
    # carries no Float64 on a Float32-only backend (Metal). On Float64 _T(x)===x
    # (byte-identical). The host pre-scan above keeps the original Float64 crits.
    _T = eltype(carbon_patch)
    _launch!(_cnprec_cn_kernel!, carbon_patch, nitrogen_patch, pc, pn,
             did_truncate, filter_bgc_vegp, patch_itype,
             _T(ccrit), _T(cnegcrit), _T(ncrit), _T(nnegcrit),
             croponly, use_nguardrail, use_matrixcn, allowneg, nc3crop;
             ndrange = length(filter_bgc_vegp))

    filter_truncatep = _build_truncate_filter(filter_bgc_vegp, did_truncate)
    return (length(filter_truncatep), filter_truncatep)
end

# ---------------------------------------------------------------------------
# truncate_c_states!
# Ported from TruncateCStates
# ---------------------------------------------------------------------------

"""
    truncate_c_states!(carbon_patch, pc, filter_bgc_vegp, patch_itype;
        ccrit, cnegcrit, croponly, allowneg, nc3crop) -> (num_truncatep, filter_truncatep)

Truncate a carbon-only state vector. If a carbon value is too small, zero it
and accumulate the removed mass into the truncation accumulator `pc`.

Returns a tuple `(num_truncatep, filter_truncatep)` listing patch indices
where truncation occurred.

Ported from `TruncateCStates` in `CNPrecisionControlMod.F90`.
"""
function truncate_c_states!(
        carbon_patch::AbstractVector,
        pc::AbstractVector,
        filter_bgc_vegp::AbstractVector{<:Integer},
        patch_itype::AbstractVector{<:Integer};
        ccrit::Real = CN_CCRIT_DEFAULT,
        cnegcrit::Real = CN_CNEGCRIT_DEFAULT,
        croponly::Bool = false,
        allowneg::Bool = false,
        nc3crop::Int = 15)

    # Sanity check: -ccrit must not be less than cnegcrit
    if -ccrit < cnegcrit
        error("cnegcrit should be less than -ccrit: cnegcrit=$cnegcrit, ccrit=$ccrit")
    end

    isempty(filter_bgc_vegp) && return (0, Int[])

    # Host pre-scan: raise on critically-negative C exactly as before.
    if !allowneg
        cH = Array(carbon_patch); fH = Array(filter_bgc_vegp); iH = Array(patch_itype)
        @inbounds for fp in eachindex(fH)
            p = fH[fp]
            if croponly && iH[p] < nc3crop
                continue
            end
            if cH[p] < cnegcrit
                error("carbon state critically negative: C=$(cH[p]), limit=$cnegcrit")
            end
        end
    end

    did_truncate = fill!(similar(carbon_patch, Bool, length(filter_bgc_vegp)), false)
    _T = eltype(carbon_patch)   # Float32 on device → no Float64 in the kernel
    _launch!(_cnprec_c_kernel!, carbon_patch, pc, did_truncate,
             filter_bgc_vegp, patch_itype, _T(ccrit), croponly, nc3crop;
             ndrange = length(filter_bgc_vegp))

    filter_truncatep = _build_truncate_filter(filter_bgc_vegp, did_truncate)
    return (length(filter_truncatep), filter_truncatep)
end

# ---------------------------------------------------------------------------
# truncate_n_states!
# Ported from TruncateNStates
# ---------------------------------------------------------------------------

"""
    truncate_n_states!(nitrogen_patch, pn, filter_bgc_vegp;
        ncrit, nnegcrit)

Truncate a nitrogen-only state vector. If a nitrogen value is too small, zero
it and accumulate the removed mass into the truncation accumulator `pn`.

Note: the Fortran source has the abort on negative N commented out; we follow
that behaviour and only warn (via no-op) rather than erroring.

Ported from `TruncateNStates` in `CNPrecisionControlMod.F90`.
"""
function truncate_n_states!(
        nitrogen_patch::AbstractVector,
        pn::AbstractVector,
        filter_bgc_vegp::AbstractVector{<:Integer};
        ncrit::Real = CN_NCRIT_DEFAULT,
        nnegcrit::Real = CN_NNEGCRIT_DEFAULT)

    isempty(filter_bgc_vegp) && return nothing

    _T = eltype(nitrogen_patch)   # Float32 on device → no Float64 in the kernel
    _launch!(_cnprec_n_kernel!, nitrogen_patch, pn, filter_bgc_vegp,
             _T(ncrit), _T(nnegcrit); ndrange = length(filter_bgc_vegp))

    return nothing
end

# ---------------------------------------------------------------------------
# truncate_additional!
# Ported from TruncateAdditional
# ---------------------------------------------------------------------------

"""
    truncate_additional!(state_patch, truncation_patch,
        num_truncatep, filter_truncatep)

Given a filter of patch indices for which truncation has already been
determined (by a prior call to `truncate_c_and_n_states!` or
`truncate_c_states!`), zero the given state and accumulate the removed mass
into `truncation_patch`.

Ported from `TruncateAdditional` in `CNPrecisionControlMod.F90`.
"""
function truncate_additional!(
        state_patch::AbstractVector,
        truncation_patch::AbstractVector,
        num_truncatep::Int,
        filter_truncatep::AbstractVector{<:Integer})

    num_truncatep == 0 && return nothing

    # filter_truncatep may carry trailing entries beyond num_truncatep; operate
    # only over the active prefix (matching the original `for fp in 1:num`).
    f = num_truncatep == length(filter_truncatep) ? filter_truncatep :
        filter_truncatep[1:num_truncatep]
    # The kernel reads the filter on the state array's backend, so move the
    # (host-built) filter to that backend. On CPU this is the identity.
    fdev = _cnprec_match_backend(state_patch, f)
    _launch!(_cnprec_additional_kernel!, state_patch, truncation_patch, fdev;
             ndrange = num_truncatep)

    return nothing
end

# ---------------------------------------------------------------------------
# cn_precision_control!
# Ported from CNPrecisionControl
# ---------------------------------------------------------------------------

"""
    cn_precision_control!(cs, ns, filter_bgc_vegp, patch_itype;
        c13cs=nothing, c14cs=nothing,
        use_c13=false, use_c14=false,
        use_crop=false, use_nguardrail=false, use_matrixcn=false,
        prec_control_for_froot=true,
        nrepr=1, nc3crop=15,
        ccrit, cnegcrit, ncrit, nnegcrit)

Apply precision control to vegetation carbon and nitrogen state variables.
For each CN pool pair, values below the critical threshold are zeroed and the
removed mass is accumulated into the `ctrunc_patch` / `ntrunc_patch` sinks.

When `use_c13` or `use_c14` is true the corresponding isotope state is also
truncated using the same filter.

Ported from `CNPrecisionControl` in `CNPrecisionControlMod.F90`.
"""
function cn_precision_control!(
        cs::CNVegCarbonStateData,
        ns::CNVegNitrogenStateData,
        filter_bgc_vegp::AbstractVector{<:Integer},
        patch_itype::AbstractVector{<:Integer};
        c13cs::Union{CNVegCarbonStateData,Nothing} = nothing,
        c14cs::Union{CNVegCarbonStateData,Nothing} = nothing,
        use_c13::Bool = false,
        use_c14::Bool = false,
        use_crop::Bool = false,
        use_nguardrail::Bool = false,
        use_matrixcn::Bool = false,
        prec_control_for_froot::Bool = true,
        nrepr::Int = NREPR,
        nc3crop::Int = 15,
        ccrit::Real = CN_CCRIT_DEFAULT,
        cnegcrit::Real = CN_CNEGCRIT_DEFAULT,
        ncrit::Real = CN_NCRIT_DEFAULT,
        nnegcrit::Real = CN_NNEGCRIT_DEFAULT)

    np = length(cs.leafc_patch)
    FT = eltype(cs.leafc_patch)

    # Initialize patch-level truncation accumulators (device-resident: match
    # the backend/eltype of the state arrays so the kernels stay on one backend).
    pc   = zero(cs.leafc_patch)
    pn   = zero(cs.leafc_patch)
    pc13 = use_c13 ? zero(cs.leafc_patch) : similar(cs.leafc_patch, 0)
    pc14 = use_c14 ? zero(cs.leafc_patch) : similar(cs.leafc_patch, 0)

    # Common keyword arguments for truncate_c_and_n_states!
    cn_kw = (ccrit=ccrit, cnegcrit=cnegcrit, ncrit=ncrit, nnegcrit=nnegcrit,
             use_nguardrail=use_nguardrail, use_matrixcn=use_matrixcn,
             nc3crop=nc3crop)

    c_kw = (ccrit=ccrit, cnegcrit=cnegcrit, nc3crop=nc3crop)

    # Helper: apply isotope truncation using a pre-computed filter
    function _trunc_isotopes!(num_tp, filter_tp, c13_field, c14_field)
        if use_c13 && c13cs !== nothing
            truncate_additional!(c13_field, pc13, num_tp, filter_tp)
        end
        if use_c14 && c14cs !== nothing
            truncate_additional!(c14_field, pc14, num_tp, filter_tp)
        end
    end

    # --- leaf C and N ---
    (num_tp, filter_tp) = truncate_c_and_n_states!(
        cs.leafc_patch, ns.leafn_patch, pc, pn,
        filter_bgc_vegp, patch_itype; cn_kw...)
    _trunc_isotopes!(num_tp, filter_tp,
        use_c13 && c13cs !== nothing ? c13cs.leafc_patch : FT[],
        use_c14 && c14cs !== nothing ? c14cs.leafc_patch : FT[])

    # --- leaf storage C and N ---
    (num_tp, filter_tp) = truncate_c_and_n_states!(
        cs.leafc_storage_patch, ns.leafn_storage_patch, pc, pn,
        filter_bgc_vegp, patch_itype; cn_kw...)
    _trunc_isotopes!(num_tp, filter_tp,
        use_c13 && c13cs !== nothing ? c13cs.leafc_storage_patch : FT[],
        use_c14 && c14cs !== nothing ? c14cs.leafc_storage_patch : FT[])

    # --- leaf transfer C and N ---
    (num_tp, filter_tp) = truncate_c_and_n_states!(
        cs.leafc_xfer_patch, ns.leafn_xfer_patch, pc, pn,
        filter_bgc_vegp, patch_itype; cn_kw...)
    _trunc_isotopes!(num_tp, filter_tp,
        use_c13 && c13cs !== nothing ? c13cs.leafc_xfer_patch : FT[],
        use_c14 && c14cs !== nothing ? c14cs.leafc_xfer_patch : FT[])

    # --- froot C and N ---
    if prec_control_for_froot
        (num_tp, filter_tp) = truncate_c_and_n_states!(
            cs.frootc_patch, ns.frootn_patch, pc, pn,
            filter_bgc_vegp, patch_itype; cn_kw..., allowneg=true)
        _trunc_isotopes!(num_tp, filter_tp,
            use_c13 && c13cs !== nothing ? c13cs.frootc_patch : FT[],
            use_c14 && c14cs !== nothing ? c14cs.frootc_patch : FT[])
    end

    # --- froot storage C and N ---
    (num_tp, filter_tp) = truncate_c_and_n_states!(
        cs.frootc_storage_patch, ns.frootn_storage_patch, pc, pn,
        filter_bgc_vegp, patch_itype; cn_kw...)
    _trunc_isotopes!(num_tp, filter_tp,
        use_c13 && c13cs !== nothing ? c13cs.frootc_storage_patch : FT[],
        use_c14 && c14cs !== nothing ? c14cs.frootc_storage_patch : FT[])

    # --- froot transfer C and N ---
    (num_tp, filter_tp) = truncate_c_and_n_states!(
        cs.frootc_xfer_patch, ns.frootn_xfer_patch, pc, pn,
        filter_bgc_vegp, patch_itype; cn_kw...)
    _trunc_isotopes!(num_tp, filter_tp,
        use_c13 && c13cs !== nothing ? c13cs.frootc_xfer_patch : FT[],
        use_c14 && c14cs !== nothing ? c14cs.frootc_xfer_patch : FT[])

    # --- crop reproductive pools (grain) ---
    if use_crop
        for k in 1:nrepr
            # grain C and N
            cs_repro = @view cs.reproductivec_patch[:, k]
            ns_repro = @view ns.reproductiven_patch[:, k]
            (num_tp, filter_tp) = truncate_c_and_n_states!(
                cs_repro, ns_repro, pc, pn,
                filter_bgc_vegp, patch_itype; cn_kw..., croponly=true)
            if use_c13 && c13cs !== nothing
                c13_repro = @view c13cs.reproductivec_patch[:, k]
                truncate_additional!(c13_repro, pc13, num_tp, filter_tp)
            end
            if use_c14 && c14cs !== nothing
                c14_repro = @view c14cs.reproductivec_patch[:, k]
                truncate_additional!(c14_repro, pc14, num_tp, filter_tp)
            end

            # grain storage C and N
            cs_repro_s = @view cs.reproductivec_storage_patch[:, k]
            ns_repro_s = @view ns.reproductiven_storage_patch[:, k]
            (num_tp, filter_tp) = truncate_c_and_n_states!(
                cs_repro_s, ns_repro_s, pc, pn,
                filter_bgc_vegp, patch_itype; cn_kw..., croponly=true)
            if use_c13 && c13cs !== nothing
                c13_repro_s = @view c13cs.reproductivec_storage_patch[:, k]
                truncate_additional!(c13_repro_s, pc13, num_tp, filter_tp)
            end
            if use_c14 && c14cs !== nothing
                c14_repro_s = @view c14cs.reproductivec_storage_patch[:, k]
                truncate_additional!(c14_repro_s, pc14, num_tp, filter_tp)
            end

            # grain transfer C and N
            cs_repro_x = @view cs.reproductivec_xfer_patch[:, k]
            ns_repro_x = @view ns.reproductiven_xfer_patch[:, k]
            (num_tp, filter_tp) = truncate_c_and_n_states!(
                cs_repro_x, ns_repro_x, pc, pn,
                filter_bgc_vegp, patch_itype; cn_kw..., croponly=true)
            if use_c13 && c13cs !== nothing
                c13_repro_x = @view c13cs.reproductivec_xfer_patch[:, k]
                truncate_additional!(c13_repro_x, pc13, num_tp, filter_tp)
            end
            if use_c14 && c14cs !== nothing
                c14_repro_x = @view c14cs.reproductivec_xfer_patch[:, k]
                truncate_additional!(c14_repro_x, pc14, num_tp, filter_tp)
            end
        end

        # cropseedc/n deficit
        (num_tp, filter_tp) = truncate_c_and_n_states!(
            cs.cropseedc_deficit_patch, ns.cropseedn_deficit_patch, pc, pn,
            filter_bgc_vegp, patch_itype; cn_kw..., allowneg=true, croponly=true)
        if use_c13 && c13cs !== nothing
            truncate_additional!(c13cs.cropseedc_deficit_patch, pc13, num_tp, filter_tp)
        end
        if use_c14 && c14cs !== nothing
            truncate_additional!(c14cs.cropseedc_deficit_patch, pc14, num_tp, filter_tp)
        end
    end

    # --- livestem C and N ---
    (num_tp, filter_tp) = truncate_c_and_n_states!(
        cs.livestemc_patch, ns.livestemn_patch, pc, pn,
        filter_bgc_vegp, patch_itype; cn_kw...)
    _trunc_isotopes!(num_tp, filter_tp,
        use_c13 && c13cs !== nothing ? c13cs.livestemc_patch : FT[],
        use_c14 && c14cs !== nothing ? c14cs.livestemc_patch : FT[])

    # --- livestem storage C and N ---
    (num_tp, filter_tp) = truncate_c_and_n_states!(
        cs.livestemc_storage_patch, ns.livestemn_storage_patch, pc, pn,
        filter_bgc_vegp, patch_itype; cn_kw...)
    _trunc_isotopes!(num_tp, filter_tp,
        use_c13 && c13cs !== nothing ? c13cs.livestemc_storage_patch : FT[],
        use_c14 && c14cs !== nothing ? c14cs.livestemc_storage_patch : FT[])

    # --- livestem transfer C and N ---
    (num_tp, filter_tp) = truncate_c_and_n_states!(
        cs.livestemc_xfer_patch, ns.livestemn_xfer_patch, pc, pn,
        filter_bgc_vegp, patch_itype; cn_kw...)
    _trunc_isotopes!(num_tp, filter_tp,
        use_c13 && c13cs !== nothing ? c13cs.livestemc_xfer_patch : FT[],
        use_c14 && c14cs !== nothing ? c14cs.livestemc_xfer_patch : FT[])

    # --- deadstem C and N ---
    (num_tp, filter_tp) = truncate_c_and_n_states!(
        cs.deadstemc_patch, ns.deadstemn_patch, pc, pn,
        filter_bgc_vegp, patch_itype; cn_kw...)
    _trunc_isotopes!(num_tp, filter_tp,
        use_c13 && c13cs !== nothing ? c13cs.deadstemc_patch : FT[],
        use_c14 && c14cs !== nothing ? c14cs.deadstemc_patch : FT[])

    # --- deadstem storage C and N ---
    (num_tp, filter_tp) = truncate_c_and_n_states!(
        cs.deadstemc_storage_patch, ns.deadstemn_storage_patch, pc, pn,
        filter_bgc_vegp, patch_itype; cn_kw...)
    _trunc_isotopes!(num_tp, filter_tp,
        use_c13 && c13cs !== nothing ? c13cs.deadstemc_storage_patch : FT[],
        use_c14 && c14cs !== nothing ? c14cs.deadstemc_storage_patch : FT[])

    # --- deadstem transfer C and N ---
    (num_tp, filter_tp) = truncate_c_and_n_states!(
        cs.deadstemc_xfer_patch, ns.deadstemn_xfer_patch, pc, pn,
        filter_bgc_vegp, patch_itype; cn_kw...)
    _trunc_isotopes!(num_tp, filter_tp,
        use_c13 && c13cs !== nothing ? c13cs.deadstemc_xfer_patch : FT[],
        use_c14 && c14cs !== nothing ? c14cs.deadstemc_xfer_patch : FT[])

    # --- livecroot C and N ---
    (num_tp, filter_tp) = truncate_c_and_n_states!(
        cs.livecrootc_patch, ns.livecrootn_patch, pc, pn,
        filter_bgc_vegp, patch_itype; cn_kw...)
    _trunc_isotopes!(num_tp, filter_tp,
        use_c13 && c13cs !== nothing ? c13cs.livecrootc_patch : FT[],
        use_c14 && c14cs !== nothing ? c14cs.livecrootc_patch : FT[])

    # --- livecroot storage C and N ---
    (num_tp, filter_tp) = truncate_c_and_n_states!(
        cs.livecrootc_storage_patch, ns.livecrootn_storage_patch, pc, pn,
        filter_bgc_vegp, patch_itype; cn_kw...)
    _trunc_isotopes!(num_tp, filter_tp,
        use_c13 && c13cs !== nothing ? c13cs.livecrootc_storage_patch : FT[],
        use_c14 && c14cs !== nothing ? c14cs.livecrootc_storage_patch : FT[])

    # --- livecroot transfer C and N ---
    (num_tp, filter_tp) = truncate_c_and_n_states!(
        cs.livecrootc_xfer_patch, ns.livecrootn_xfer_patch, pc, pn,
        filter_bgc_vegp, patch_itype; cn_kw...)
    _trunc_isotopes!(num_tp, filter_tp,
        use_c13 && c13cs !== nothing ? c13cs.livecrootc_xfer_patch : FT[],
        use_c14 && c14cs !== nothing ? c14cs.livecrootc_xfer_patch : FT[])

    # --- deadcroot C and N ---
    (num_tp, filter_tp) = truncate_c_and_n_states!(
        cs.deadcrootc_patch, ns.deadcrootn_patch, pc, pn,
        filter_bgc_vegp, patch_itype; cn_kw...)
    _trunc_isotopes!(num_tp, filter_tp,
        use_c13 && c13cs !== nothing ? c13cs.deadcrootc_patch : FT[],
        use_c14 && c14cs !== nothing ? c14cs.deadcrootc_patch : FT[])

    # --- deadcroot storage C and N ---
    (num_tp, filter_tp) = truncate_c_and_n_states!(
        cs.deadcrootc_storage_patch, ns.deadcrootn_storage_patch, pc, pn,
        filter_bgc_vegp, patch_itype; cn_kw...)
    _trunc_isotopes!(num_tp, filter_tp,
        use_c13 && c13cs !== nothing ? c13cs.deadcrootc_storage_patch : FT[],
        use_c14 && c14cs !== nothing ? c14cs.deadcrootc_storage_patch : FT[])

    # --- deadcroot transfer C and N ---
    (num_tp, filter_tp) = truncate_c_and_n_states!(
        cs.deadcrootc_xfer_patch, ns.deadcrootn_xfer_patch, pc, pn,
        filter_bgc_vegp, patch_itype; cn_kw...)
    _trunc_isotopes!(num_tp, filter_tp,
        use_c13 && c13cs !== nothing ? c13cs.deadcrootc_xfer_patch : FT[],
        use_c14 && c14cs !== nothing ? c14cs.deadcrootc_xfer_patch : FT[])

    # --- gresp_storage (C only) ---
    (num_tp, filter_tp) = truncate_c_states!(
        cs.gresp_storage_patch, pc,
        filter_bgc_vegp, patch_itype; c_kw...)
    _trunc_isotopes!(num_tp, filter_tp,
        use_c13 && c13cs !== nothing ? c13cs.gresp_storage_patch : FT[],
        use_c14 && c14cs !== nothing ? c14cs.gresp_storage_patch : FT[])

    # --- gresp_xfer (C only) ---
    (num_tp, filter_tp) = truncate_c_states!(
        cs.gresp_xfer_patch, pc,
        filter_bgc_vegp, patch_itype; c_kw...)
    _trunc_isotopes!(num_tp, filter_tp,
        use_c13 && c13cs !== nothing ? c13cs.gresp_xfer_patch : FT[],
        use_c14 && c14cs !== nothing ? c14cs.gresp_xfer_patch : FT[])

    # --- cpool (C only) ---
    (num_tp, filter_tp) = truncate_c_states!(
        cs.cpool_patch, pc,
        filter_bgc_vegp, patch_itype; c_kw...)
    _trunc_isotopes!(num_tp, filter_tp,
        use_c13 && c13cs !== nothing ? c13cs.cpool_patch : FT[],
        use_c14 && c14cs !== nothing ? c14cs.cpool_patch : FT[])

    # --- xsmrpool (C only, crop only, allow negative) ---
    if use_crop
        (num_tp, filter_tp) = truncate_c_states!(
            cs.xsmrpool_patch, pc,
            filter_bgc_vegp, patch_itype; c_kw..., allowneg=true, croponly=true)
        if use_c13 && c13cs !== nothing
            truncate_additional!(c13cs.xsmrpool_patch, pc13, num_tp, filter_tp)
        end
        if use_c14 && c14cs !== nothing
            truncate_additional!(c14cs.xsmrpool_patch, pc14, num_tp, filter_tp)
        end
    end

    # --- retransn (N only) ---
    truncate_n_states!(ns.retransn_patch, pn, filter_bgc_vegp;
        ncrit=ncrit, nnegcrit=nnegcrit)

    # --- npool (N only) ---
    truncate_n_states!(ns.npool_patch, pn, filter_bgc_vegp;
        ncrit=ncrit, nnegcrit=nnegcrit)

    # --- Accumulate truncation into ctrunc/ntrunc sinks ---
    if !isempty(filter_bgc_vegp)
        _launch!(_cnprec_accum_kernel!, cs.ctrunc_patch, ns.ntrunc_patch,
                 pc, pn, filter_bgc_vegp; ndrange = length(filter_bgc_vegp))
        if use_c13 && c13cs !== nothing
            _launch!(_cnprec_accum_iso_kernel!, c13cs.ctrunc_patch, pc13,
                     filter_bgc_vegp; ndrange = length(filter_bgc_vegp))
        end
        if use_c14 && c14cs !== nothing
            _launch!(_cnprec_accum_iso_kernel!, c14cs.ctrunc_patch, pc14,
                     filter_bgc_vegp; ndrange = length(filter_bgc_vegp))
        end
    end

    return nothing
end
