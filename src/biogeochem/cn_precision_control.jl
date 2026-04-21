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
# ==========================================================================

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
        carbon_patch::AbstractVector{Float64},
        nitrogen_patch::AbstractVector{Float64},
        pc::Vector{<:Real},
        pn::Vector{<:Real},
        filter_bgc_vegp::AbstractVector{Int},
        patch_itype::AbstractVector{Int};
        ccrit::Real = CN_CCRIT_DEFAULT,
        cnegcrit::Real = CN_CNEGCRIT_DEFAULT,
        ncrit::Real = CN_NCRIT_DEFAULT,
        nnegcrit::Real = CN_NNEGCRIT_DEFAULT,
        use_nguardrail::Bool = false,
        use_matrixcn::Bool = false,
        croponly::Bool = false,
        allowneg::Bool = false,
        nc3crop::Int = 15)

    filter_truncatep = Int[]

    for fp in eachindex(filter_bgc_vegp)
        p = filter_bgc_vegp[fp]

        # Skip non-crop patches when croponly is set
        if croponly && patch_itype[p] < nc3crop
            continue
        end

        if !allowneg && (carbon_patch[p] < cnegcrit || nitrogen_patch[p] < nnegcrit)
            error("carbon or nitrogen state critically negative: " *
                  "C=$(carbon_patch[p]), N=$(nitrogen_patch[p]), " *
                  "limits: cnegcrit=$cnegcrit, nnegcrit=$nnegcrit")
        else
            if use_matrixcn
                # Matrix code has a different check
                if (carbon_patch[p] < ccrit && carbon_patch[p] > -ccrit * 1.0e+6) ||
                   (use_nguardrail && nitrogen_patch[p] < ncrit && nitrogen_patch[p] > -ncrit * 1.0e+6)
                    push!(filter_truncatep, p)

                    pc[p] += carbon_patch[p]
                    carbon_patch[p] = 0.0

                    pn[p] += nitrogen_patch[p]
                    nitrogen_patch[p] = 0.0
                end
            else
                if abs(carbon_patch[p]) < ccrit ||
                   (use_nguardrail && abs(nitrogen_patch[p]) < ncrit)
                    push!(filter_truncatep, p)

                    pc[p] += carbon_patch[p]
                    carbon_patch[p] = 0.0

                    pn[p] += nitrogen_patch[p]
                    nitrogen_patch[p] = 0.0
                end
            end
        end
    end

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
        carbon_patch::AbstractVector{Float64},
        pc::Vector{<:Real},
        filter_bgc_vegp::AbstractVector{Int},
        patch_itype::AbstractVector{Int};
        ccrit::Real = CN_CCRIT_DEFAULT,
        cnegcrit::Real = CN_CNEGCRIT_DEFAULT,
        croponly::Bool = false,
        allowneg::Bool = false,
        nc3crop::Int = 15)

    # Sanity check: -ccrit must not be less than cnegcrit
    if -ccrit < cnegcrit
        error("cnegcrit should be less than -ccrit: cnegcrit=$cnegcrit, ccrit=$ccrit")
    end

    filter_truncatep = Int[]

    for fp in eachindex(filter_bgc_vegp)
        p = filter_bgc_vegp[fp]

        # Skip non-crop patches when croponly is set
        if croponly && patch_itype[p] < nc3crop
            continue
        end

        if !allowneg && carbon_patch[p] < cnegcrit
            error("carbon state critically negative: C=$(carbon_patch[p]), limit=$cnegcrit")
        elseif abs(carbon_patch[p]) < ccrit
            push!(filter_truncatep, p)

            pc[p] += carbon_patch[p]
            carbon_patch[p] = 0.0
        end
    end

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
        nitrogen_patch::AbstractVector{Float64},
        pn::Vector{<:Real},
        filter_bgc_vegp::AbstractVector{Int};
        ncrit::Real = CN_NCRIT_DEFAULT,
        nnegcrit::Real = CN_NNEGCRIT_DEFAULT)

    for fp in eachindex(filter_bgc_vegp)
        p = filter_bgc_vegp[fp]

        if nitrogen_patch[p] < nnegcrit
            # Fortran source has the endrun call commented out here.
            # We follow suit and do nothing (no abort).
        elseif abs(nitrogen_patch[p]) < ncrit
            pn[p] += nitrogen_patch[p]
            nitrogen_patch[p] = 0.0
        end
    end

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
        state_patch::AbstractVector{Float64},
        truncation_patch::Vector{<:Real},
        num_truncatep::Int,
        filter_truncatep::AbstractVector{Int})

    for fp in 1:num_truncatep
        p = filter_truncatep[fp]
        truncation_patch[p] += state_patch[p]
        state_patch[p] = 0.0
    end

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
        filter_bgc_vegp::AbstractVector{Int},
        patch_itype::AbstractVector{Int};
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

    # Initialize patch-level truncation accumulators
    pc   = zeros(FT, np)
    pn   = zeros(FT, np)
    pc13 = use_c13 ? zeros(FT, np) : FT[]
    pc14 = use_c14 ? zeros(FT, np) : FT[]

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
        use_c13 && c13cs !== nothing ? c13cs.leafc_patch : Float64[],
        use_c14 && c14cs !== nothing ? c14cs.leafc_patch : Float64[])

    # --- leaf storage C and N ---
    (num_tp, filter_tp) = truncate_c_and_n_states!(
        cs.leafc_storage_patch, ns.leafn_storage_patch, pc, pn,
        filter_bgc_vegp, patch_itype; cn_kw...)
    _trunc_isotopes!(num_tp, filter_tp,
        use_c13 && c13cs !== nothing ? c13cs.leafc_storage_patch : Float64[],
        use_c14 && c14cs !== nothing ? c14cs.leafc_storage_patch : Float64[])

    # --- leaf transfer C and N ---
    (num_tp, filter_tp) = truncate_c_and_n_states!(
        cs.leafc_xfer_patch, ns.leafn_xfer_patch, pc, pn,
        filter_bgc_vegp, patch_itype; cn_kw...)
    _trunc_isotopes!(num_tp, filter_tp,
        use_c13 && c13cs !== nothing ? c13cs.leafc_xfer_patch : Float64[],
        use_c14 && c14cs !== nothing ? c14cs.leafc_xfer_patch : Float64[])

    # --- froot C and N ---
    if prec_control_for_froot
        (num_tp, filter_tp) = truncate_c_and_n_states!(
            cs.frootc_patch, ns.frootn_patch, pc, pn,
            filter_bgc_vegp, patch_itype; cn_kw..., allowneg=true)
        _trunc_isotopes!(num_tp, filter_tp,
            use_c13 && c13cs !== nothing ? c13cs.frootc_patch : Float64[],
            use_c14 && c14cs !== nothing ? c14cs.frootc_patch : Float64[])
    end

    # --- froot storage C and N ---
    (num_tp, filter_tp) = truncate_c_and_n_states!(
        cs.frootc_storage_patch, ns.frootn_storage_patch, pc, pn,
        filter_bgc_vegp, patch_itype; cn_kw...)
    _trunc_isotopes!(num_tp, filter_tp,
        use_c13 && c13cs !== nothing ? c13cs.frootc_storage_patch : Float64[],
        use_c14 && c14cs !== nothing ? c14cs.frootc_storage_patch : Float64[])

    # --- froot transfer C and N ---
    (num_tp, filter_tp) = truncate_c_and_n_states!(
        cs.frootc_xfer_patch, ns.frootn_xfer_patch, pc, pn,
        filter_bgc_vegp, patch_itype; cn_kw...)
    _trunc_isotopes!(num_tp, filter_tp,
        use_c13 && c13cs !== nothing ? c13cs.frootc_xfer_patch : Float64[],
        use_c14 && c14cs !== nothing ? c14cs.frootc_xfer_patch : Float64[])

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
        use_c13 && c13cs !== nothing ? c13cs.livestemc_patch : Float64[],
        use_c14 && c14cs !== nothing ? c14cs.livestemc_patch : Float64[])

    # --- livestem storage C and N ---
    (num_tp, filter_tp) = truncate_c_and_n_states!(
        cs.livestemc_storage_patch, ns.livestemn_storage_patch, pc, pn,
        filter_bgc_vegp, patch_itype; cn_kw...)
    _trunc_isotopes!(num_tp, filter_tp,
        use_c13 && c13cs !== nothing ? c13cs.livestemc_storage_patch : Float64[],
        use_c14 && c14cs !== nothing ? c14cs.livestemc_storage_patch : Float64[])

    # --- livestem transfer C and N ---
    (num_tp, filter_tp) = truncate_c_and_n_states!(
        cs.livestemc_xfer_patch, ns.livestemn_xfer_patch, pc, pn,
        filter_bgc_vegp, patch_itype; cn_kw...)
    _trunc_isotopes!(num_tp, filter_tp,
        use_c13 && c13cs !== nothing ? c13cs.livestemc_xfer_patch : Float64[],
        use_c14 && c14cs !== nothing ? c14cs.livestemc_xfer_patch : Float64[])

    # --- deadstem C and N ---
    (num_tp, filter_tp) = truncate_c_and_n_states!(
        cs.deadstemc_patch, ns.deadstemn_patch, pc, pn,
        filter_bgc_vegp, patch_itype; cn_kw...)
    _trunc_isotopes!(num_tp, filter_tp,
        use_c13 && c13cs !== nothing ? c13cs.deadstemc_patch : Float64[],
        use_c14 && c14cs !== nothing ? c14cs.deadstemc_patch : Float64[])

    # --- deadstem storage C and N ---
    (num_tp, filter_tp) = truncate_c_and_n_states!(
        cs.deadstemc_storage_patch, ns.deadstemn_storage_patch, pc, pn,
        filter_bgc_vegp, patch_itype; cn_kw...)
    _trunc_isotopes!(num_tp, filter_tp,
        use_c13 && c13cs !== nothing ? c13cs.deadstemc_storage_patch : Float64[],
        use_c14 && c14cs !== nothing ? c14cs.deadstemc_storage_patch : Float64[])

    # --- deadstem transfer C and N ---
    (num_tp, filter_tp) = truncate_c_and_n_states!(
        cs.deadstemc_xfer_patch, ns.deadstemn_xfer_patch, pc, pn,
        filter_bgc_vegp, patch_itype; cn_kw...)
    _trunc_isotopes!(num_tp, filter_tp,
        use_c13 && c13cs !== nothing ? c13cs.deadstemc_xfer_patch : Float64[],
        use_c14 && c14cs !== nothing ? c14cs.deadstemc_xfer_patch : Float64[])

    # --- livecroot C and N ---
    (num_tp, filter_tp) = truncate_c_and_n_states!(
        cs.livecrootc_patch, ns.livecrootn_patch, pc, pn,
        filter_bgc_vegp, patch_itype; cn_kw...)
    _trunc_isotopes!(num_tp, filter_tp,
        use_c13 && c13cs !== nothing ? c13cs.livecrootc_patch : Float64[],
        use_c14 && c14cs !== nothing ? c14cs.livecrootc_patch : Float64[])

    # --- livecroot storage C and N ---
    (num_tp, filter_tp) = truncate_c_and_n_states!(
        cs.livecrootc_storage_patch, ns.livecrootn_storage_patch, pc, pn,
        filter_bgc_vegp, patch_itype; cn_kw...)
    _trunc_isotopes!(num_tp, filter_tp,
        use_c13 && c13cs !== nothing ? c13cs.livecrootc_storage_patch : Float64[],
        use_c14 && c14cs !== nothing ? c14cs.livecrootc_storage_patch : Float64[])

    # --- livecroot transfer C and N ---
    (num_tp, filter_tp) = truncate_c_and_n_states!(
        cs.livecrootc_xfer_patch, ns.livecrootn_xfer_patch, pc, pn,
        filter_bgc_vegp, patch_itype; cn_kw...)
    _trunc_isotopes!(num_tp, filter_tp,
        use_c13 && c13cs !== nothing ? c13cs.livecrootc_xfer_patch : Float64[],
        use_c14 && c14cs !== nothing ? c14cs.livecrootc_xfer_patch : Float64[])

    # --- deadcroot C and N ---
    (num_tp, filter_tp) = truncate_c_and_n_states!(
        cs.deadcrootc_patch, ns.deadcrootn_patch, pc, pn,
        filter_bgc_vegp, patch_itype; cn_kw...)
    _trunc_isotopes!(num_tp, filter_tp,
        use_c13 && c13cs !== nothing ? c13cs.deadcrootc_patch : Float64[],
        use_c14 && c14cs !== nothing ? c14cs.deadcrootc_patch : Float64[])

    # --- deadcroot storage C and N ---
    (num_tp, filter_tp) = truncate_c_and_n_states!(
        cs.deadcrootc_storage_patch, ns.deadcrootn_storage_patch, pc, pn,
        filter_bgc_vegp, patch_itype; cn_kw...)
    _trunc_isotopes!(num_tp, filter_tp,
        use_c13 && c13cs !== nothing ? c13cs.deadcrootc_storage_patch : Float64[],
        use_c14 && c14cs !== nothing ? c14cs.deadcrootc_storage_patch : Float64[])

    # --- deadcroot transfer C and N ---
    (num_tp, filter_tp) = truncate_c_and_n_states!(
        cs.deadcrootc_xfer_patch, ns.deadcrootn_xfer_patch, pc, pn,
        filter_bgc_vegp, patch_itype; cn_kw...)
    _trunc_isotopes!(num_tp, filter_tp,
        use_c13 && c13cs !== nothing ? c13cs.deadcrootc_xfer_patch : Float64[],
        use_c14 && c14cs !== nothing ? c14cs.deadcrootc_xfer_patch : Float64[])

    # --- gresp_storage (C only) ---
    (num_tp, filter_tp) = truncate_c_states!(
        cs.gresp_storage_patch, pc,
        filter_bgc_vegp, patch_itype; c_kw...)
    _trunc_isotopes!(num_tp, filter_tp,
        use_c13 && c13cs !== nothing ? c13cs.gresp_storage_patch : Float64[],
        use_c14 && c14cs !== nothing ? c14cs.gresp_storage_patch : Float64[])

    # --- gresp_xfer (C only) ---
    (num_tp, filter_tp) = truncate_c_states!(
        cs.gresp_xfer_patch, pc,
        filter_bgc_vegp, patch_itype; c_kw...)
    _trunc_isotopes!(num_tp, filter_tp,
        use_c13 && c13cs !== nothing ? c13cs.gresp_xfer_patch : Float64[],
        use_c14 && c14cs !== nothing ? c14cs.gresp_xfer_patch : Float64[])

    # --- cpool (C only) ---
    (num_tp, filter_tp) = truncate_c_states!(
        cs.cpool_patch, pc,
        filter_bgc_vegp, patch_itype; c_kw...)
    _trunc_isotopes!(num_tp, filter_tp,
        use_c13 && c13cs !== nothing ? c13cs.cpool_patch : Float64[],
        use_c14 && c14cs !== nothing ? c14cs.cpool_patch : Float64[])

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
    for fp in eachindex(filter_bgc_vegp)
        p = filter_bgc_vegp[fp]

        cs.ctrunc_patch[p] += pc[p]
        ns.ntrunc_patch[p] += pn[p]

        if use_c13 && c13cs !== nothing
            c13cs.ctrunc_patch[p] += pc13[p]
        end
        if use_c14 && c14cs !== nothing
            c14cs.ctrunc_patch[p] += pc14[p]
        end
    end

    return nothing
end
