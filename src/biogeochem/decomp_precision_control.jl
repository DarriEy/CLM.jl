# ==========================================================================
# Ported from: src/soilbiogeochem/SoilBiogeochemPrecisionControlMod.F90
# Controls on very low values in critical state variables.
# Truncates small C/N pool values to zero for numerical stability and
# accumulates the truncated mass into truncation sink terms.
#
# Public functions:
#   soil_bgc_precision_control_init!  -- Initialize critical thresholds
#   soil_bgc_precision_control!       -- Apply precision control
# ==========================================================================

# ---------------------------------------------------------------------------
# Module-level critical thresholds (mutable for initialization)
# Ported from ccrit/ncrit in SoilBiogeochemPrecisionControlMod.F90
# ---------------------------------------------------------------------------

"""
    ccrit

Critical carbon state value for truncation (gC/m3).
Pool values with `abs(val) < ccrit` are zeroed and the removed mass is
tracked in `ctrunc_vr_col`.
"""
const CCRIT_DEFAULT = 1.0e-8

"""
    ncrit

Critical nitrogen state value for truncation (gN/m3).
Pool values with `abs(val) < ncrit` are zeroed and the removed mass is
tracked in `ntrunc_vr_col`.
"""
const NCRIT_DEFAULT = 1.0e-8

# ---------------------------------------------------------------------------
# soil_bgc_precision_control_init!
# Ported from SoilBiogeochemPrecisionControlInit
# ---------------------------------------------------------------------------

"""
    soil_bgc_precision_control_init!(cs, ns;
        c13cs=nothing, c14cs=nothing,
        use_c13=false, use_c14=false,
        totvegcthresh=1.0)

Initialize precision control thresholds and set the total vegetation carbon
threshold on each state type instance.

Returns `(ccrit, ncrit)` -- the critical truncation thresholds.

Ported from `SoilBiogeochemPrecisionControlInit` in
`SoilBiogeochemPrecisionControlMod.F90`.
"""
function soil_bgc_precision_control_init!(
        cs::SoilBiogeochemCarbonStateData,
        ns::SoilBiogeochemNitrogenStateData;
        c13cs::Union{SoilBiogeochemCarbonStateData,Nothing}=nothing,
        c14cs::Union{SoilBiogeochemCarbonStateData,Nothing}=nothing,
        use_c13::Bool=false,
        use_c14::Bool=false,
        totvegcthresh::Real=1.0)

    ccrit = CCRIT_DEFAULT
    ncrit = NCRIT_DEFAULT

    soil_bgc_carbon_state_set_totvegcthresh!(cs, totvegcthresh)
    if use_c13 && c13cs !== nothing
        soil_bgc_carbon_state_set_totvegcthresh!(c13cs, totvegcthresh)
    end
    if use_c14 && c14cs !== nothing
        soil_bgc_carbon_state_set_totvegcthresh!(c14cs, totvegcthresh)
    end
    soil_bgc_nitrogen_state_set_totvegcthresh!(ns, totvegcthresh)

    return (ccrit, ncrit)
end

# ---------------------------------------------------------------------------
# soil_bgc_precision_control!
# Ported from SoilBiogeochemPrecisionControl
#
# GPU kernelization (Phase B, BGC C/N precision control):
#   * Main column loop -> `_decprec_col_kernel!` (one thread per soil column; the
#     internal decomp-level j / pool k loops run sequentially in-thread). Each
#     thread owns its column c, and the per-level truncation accumulators
#     cc/cn/cc13/cc14 are thread-local, written only into the column's OWN
#     ctrunc_vr_col[c,j] / ntrunc_vr_col[c,j] sinks — an own-index per-level
#     accumulation, NOT a cross-thread scatter, so it is byte-identical and
#     race-free.
#   * Mineral-N block -> `_decprec_nit_col_kernel!` (one thread per column; the
#     internal j loop runs in-thread; per-(c,j) own-index NO3/NH4 resets).
# The C13/C14 isotope arrays are always passed to the column kernel (with the C12
# carbon arrays standing in as placeholders when the isotope is inactive) and the
# corresponding writes are gated by the Bool flags `do_c13`/`do_c14`, so no Float64
# `nothing` ever reaches the device. The ccrit / ncrit thresholds are converted to
# the state element type at the launch site (a Float64 scalar kernel arg is invalid
# Metal IR). Float64 literals inside the kernels are eltype-converted (`zero(T)`),
# byte-identical on the Float64 CPU path.
# ---------------------------------------------------------------------------

# --- Kernel: per-column C/N (+ optional C13/C14) precision truncation ---
@kernel function _decprec_col_kernel!(@Const(mask_bgc_soilc),
        decomp_cpools_vr_col, decomp_npools_vr_col, ctrunc_vr_col, ntrunc_vr_col,
        c13_decomp_cpools_vr_col, c13_ctrunc_vr_col,
        c14_decomp_cpools_vr_col, c14_ctrunc_vr_col,
        nlevdecomp::Int, ndecomp_pools::Int, ccrit,
        do_c13::Bool, do_c14::Bool)
    c = @index(Global)
    @inbounds if mask_bgc_soilc[c]
        T = typeof(ccrit)
        for j in 1:nlevdecomp
            # Initialize column-level truncation accumulators for this level
            cc = zero(T)
            cn = zero(T)
            cc13 = zero(T)
            cc14 = zero(T)

            # All decomposing pools C and N
            for k in 1:ndecomp_pools
                if abs(decomp_cpools_vr_col[c, j, k]) < ccrit
                    # Accumulate and zero the C pool
                    cc += decomp_cpools_vr_col[c, j, k]
                    decomp_cpools_vr_col[c, j, k] = zero(T)

                    # Accumulate and zero the N pool
                    cn += decomp_npools_vr_col[c, j, k]
                    decomp_npools_vr_col[c, j, k] = zero(T)

                    # C13 isotope
                    if do_c13
                        cc13 += c13_decomp_cpools_vr_col[c, j, k]
                        c13_decomp_cpools_vr_col[c, j, k] = zero(T)
                    end

                    # C14 isotope
                    if do_c14
                        cc14 += c14_decomp_cpools_vr_col[c, j, k]
                        c14_decomp_cpools_vr_col[c, j, k] = zero(T)
                    end
                end
            end

            # Accumulate truncated mass into truncation sinks
            ctrunc_vr_col[c, j] += cc
            ntrunc_vr_col[c, j] += cn

            if do_c13
                c13_ctrunc_vr_col[c, j] += cc13
            end
            if do_c14
                c14_ctrunc_vr_col[c, j] += cc14
            end
        end
    end
end

# --- Kernel: per-column mineral-N (NO3/NH4) small-negative stability reset ---
@kernel function _decprec_nit_col_kernel!(@Const(mask_bgc_soilc),
        smin_no3_vr_col, smin_nh4_vr_col, nlevdecomp::Int, ncrit_nit)
    c = @index(Global)
    @inbounds if mask_bgc_soilc[c]
        T = typeof(ncrit_nit)
        for j in 1:nlevdecomp
            # NO3 stability
            if abs(smin_no3_vr_col[c, j]) < ncrit_nit
                if smin_no3_vr_col[c, j] < zero(T)
                    smin_no3_vr_col[c, j] = zero(T)
                end
            end

            # NH4 stability
            if abs(smin_nh4_vr_col[c, j]) < ncrit_nit
                if smin_nh4_vr_col[c, j] < zero(T)
                    smin_nh4_vr_col[c, j] = zero(T)
                end
            end
        end
    end
end

"""
    soil_bgc_precision_control!(cs, ns;
        mask_bgc_soilc, nlevdecomp, ndecomp_pools,
        ccrit=CCRIT_DEFAULT, ncrit=NCRIT_DEFAULT,
        c13cs=nothing, c14cs=nothing,
        use_c13=false, use_c14=false,
        use_nitrif_denitrif=false, use_fun=false)

On the radiation timestep, zero out decomposition C and N pool values that
fall below the critical threshold and accumulate the removed mass in the
corresponding truncation sink (`ctrunc_vr_col` / `ntrunc_vr_col`).

When `use_nitrif_denitrif` is true and `use_fun` is false, small negative
mineral NO3 and NH4 values are also reset to zero for stability.

Ported from `SoilBiogeochemPrecisionControl` in
`SoilBiogeochemPrecisionControlMod.F90`. Runs as two KernelAbstractions kernels
(one thread per soil column each); the internal decomp-level / pool loops stay
sequential in-thread with own-index truncation sinks, so the CPU path is
byte-identical and the GPU path is race-free.
"""
function soil_bgc_precision_control!(
        cs::SoilBiogeochemCarbonStateData,
        ns::SoilBiogeochemNitrogenStateData;
        mask_bgc_soilc::AbstractVector{Bool},
        nlevdecomp::Int,
        ndecomp_pools::Int,
        ccrit::Real=CCRIT_DEFAULT,
        ncrit::Real=NCRIT_DEFAULT,
        c13cs::Union{SoilBiogeochemCarbonStateData,Nothing}=nothing,
        c14cs::Union{SoilBiogeochemCarbonStateData,Nothing}=nothing,
        use_c13::Bool=false,
        use_c14::Bool=false,
        use_nitrif_denitrif::Bool=false,
        use_fun::Bool=false)

    # ------------------------------------------------------------------
    # Main column loop: truncate tiny C/N decomposition pool values
    # ------------------------------------------------------------------
    do_c13 = use_c13 && c13cs !== nothing
    do_c14 = use_c14 && c14cs !== nothing

    # Isotope arrays are always passed to the kernel; when an isotope is inactive
    # the C12 carbon arrays stand in as placeholders (never written — gated by the
    # do_c13/do_c14 Bool flags). This keeps `nothing` off the device.
    c13_pools  = do_c13 ? c13cs.decomp_cpools_vr_col : cs.decomp_cpools_vr_col
    c13_trunc  = do_c13 ? c13cs.ctrunc_vr_col        : cs.ctrunc_vr_col
    c14_pools  = do_c14 ? c14cs.decomp_cpools_vr_col : cs.decomp_cpools_vr_col
    c14_trunc  = do_c14 ? c14cs.ctrunc_vr_col        : cs.ctrunc_vr_col

    # Convert the scalar threshold to the state element type at the launch site
    # (a Float64 scalar kernel arg is invalid Metal IR).
    _T = eltype(cs.decomp_cpools_vr_col)
    ccrit_T = _T(ccrit)

    _launch!(_decprec_col_kernel!, mask_bgc_soilc,
        cs.decomp_cpools_vr_col, ns.decomp_npools_vr_col,
        cs.ctrunc_vr_col, ns.ntrunc_vr_col,
        c13_pools, c13_trunc, c14_pools, c14_trunc,
        nlevdecomp, ndecomp_pools, ccrit_T, do_c13, do_c14)

    # ------------------------------------------------------------------
    # Mineral N stability: remove small negative NO3/NH4 perturbations
    # Only when nitrif_denitrif is active and FUN is not in use.
    # ------------------------------------------------------------------
    if !use_fun && use_nitrif_denitrif
        ncrit_nit_T = _T(ncrit / 1.0e4)  # tighter threshold for mineral N

        _launch!(_decprec_nit_col_kernel!, mask_bgc_soilc,
            ns.smin_no3_vr_col, ns.smin_nh4_vr_col, nlevdecomp, ncrit_nit_T)
    end

    return nothing
end
