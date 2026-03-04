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
        totvegcthresh::Float64=1.0)

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
# ---------------------------------------------------------------------------

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
`SoilBiogeochemPrecisionControlMod.F90`.
"""
function soil_bgc_precision_control!(
        cs::SoilBiogeochemCarbonStateData,
        ns::SoilBiogeochemNitrogenStateData;
        mask_bgc_soilc::BitVector,
        nlevdecomp::Int,
        ndecomp_pools::Int,
        ccrit::Float64=CCRIT_DEFAULT,
        ncrit::Float64=NCRIT_DEFAULT,
        c13cs::Union{SoilBiogeochemCarbonStateData,Nothing}=nothing,
        c14cs::Union{SoilBiogeochemCarbonStateData,Nothing}=nothing,
        use_c13::Bool=false,
        use_c14::Bool=false,
        use_nitrif_denitrif::Bool=false,
        use_fun::Bool=false)

    # ------------------------------------------------------------------
    # Main column loop: truncate tiny C/N decomposition pool values
    # ------------------------------------------------------------------
    for c in eachindex(mask_bgc_soilc)
        mask_bgc_soilc[c] || continue

        for j in 1:nlevdecomp
            # Initialize column-level truncation accumulators for this level
            cc = 0.0
            cn = 0.0
            cc13 = 0.0
            cc14 = 0.0

            # All decomposing pools C and N
            for k in 1:ndecomp_pools

                if abs(cs.decomp_cpools_vr_col[c, j, k]) < ccrit
                    # Accumulate and zero the C pool
                    cc += cs.decomp_cpools_vr_col[c, j, k]
                    cs.decomp_cpools_vr_col[c, j, k] = 0.0

                    # Accumulate and zero the N pool
                    cn += ns.decomp_npools_vr_col[c, j, k]
                    ns.decomp_npools_vr_col[c, j, k] = 0.0

                    # C13 isotope
                    if use_c13 && c13cs !== nothing
                        cc13 += c13cs.decomp_cpools_vr_col[c, j, k]
                        c13cs.decomp_cpools_vr_col[c, j, k] = 0.0
                    end

                    # C14 isotope
                    if use_c14 && c14cs !== nothing
                        cc14 += c14cs.decomp_cpools_vr_col[c, j, k]
                        c14cs.decomp_cpools_vr_col[c, j, k] = 0.0
                    end
                end
            end

            # Accumulate truncated mass into truncation sinks
            cs.ctrunc_vr_col[c, j] += cc
            ns.ntrunc_vr_col[c, j] += cn

            if use_c13 && c13cs !== nothing
                c13cs.ctrunc_vr_col[c, j] += cc13
            end
            if use_c14 && c14cs !== nothing
                c14cs.ctrunc_vr_col[c, j] += cc14
            end
        end
    end

    # ------------------------------------------------------------------
    # Mineral N stability: remove small negative NO3/NH4 perturbations
    # Only when nitrif_denitrif is active and FUN is not in use.
    # ------------------------------------------------------------------
    if !use_fun && use_nitrif_denitrif
        ncrit_nit = ncrit / 1.0e4  # tighter threshold for mineral N

        for c in eachindex(mask_bgc_soilc)
            mask_bgc_soilc[c] || continue

            for j in 1:nlevdecomp
                # NO3 stability
                if abs(ns.smin_no3_vr_col[c, j]) < ncrit_nit
                    if ns.smin_no3_vr_col[c, j] < 0.0
                        ns.smin_no3_vr_col[c, j] = 0.0
                    end
                end

                # NH4 stability
                if abs(ns.smin_nh4_vr_col[c, j]) < ncrit_nit
                    if ns.smin_nh4_vr_col[c, j] < 0.0
                        ns.smin_nh4_vr_col[c, j] = 0.0
                    end
                end
            end
        end
    end

    return nothing
end
