# ==========================================================================
# Ported from: src/soilbiogeochem/SoilBiogeochemPotentialMod.F90
# Calculate potential decomp rates and total immobilization demand.
#
# Public functions:
#   soil_bgc_potential!  -- Calculate potential decomposition rates and
#                           total immobilization demand
# ==========================================================================

# ---------------------------------------------------------------------------
# DecompPotentialParams -- parameters for potential decomposition
# Ported from params_type in SoilBiogeochemPotentialMod.F90
# ---------------------------------------------------------------------------

"""
    DecompPotentialParams

Parameters for potential decomposition calculations.
Holds denitrification proportion parameter.

Ported from `params_type` in `SoilBiogeochemPotentialMod.F90`.
"""
Base.@kwdef mutable struct DecompPotentialParams
    dnp::Float64 = 0.01  # denitrification proportion
end

"""
    decomp_potential_read_params!(params; dnp)

Read potential decomposition parameters.
Corresponds to `readParams` in `SoilBiogeochemPotentialMod.F90`.
"""
function decomp_potential_read_params!(params::DecompPotentialParams; dnp::Real)
    params.dnp = dnp
    return nothing
end

# ---------------------------------------------------------------------------
# soil_bgc_potential! -- Main potential decomposition subroutine
# Ported from SoilBiogeochemPotential in SoilBiogeochemPotentialMod.F90
#
# Calculates potential decomp rates, C:N ratios of pools, and total
# immobilization demand prior to competition for mineral N.
# ---------------------------------------------------------------------------

"""
    soil_bgc_potential!(cf, cs, nf, ns, st, cascade_con;
        mask_bgc_soilc, bounds, nlevdecomp, ndecomp_pools,
        ndecomp_cascade_transitions,
        cn_decomp_pools, p_decomp_cpool_loss, p_decomp_cn_gain,
        pmnf_decomp_cascade, p_decomp_npool_to_din,
        use_mimics)

Calculate potential decomposition rates and total immobilization demand.

Given decomposition rate coefficients (`decomp_k_col`) and respiration fractions
(`rf_decomp_cascade_col`) previously computed, this routine calculates:
- C:N ratios of decomposing pools
- Potential C loss from each pool for each transition
- Potential mineral N flux (positive = immobilization, negative = mineralization)
- For MIMICS: potential N flux to dissolved inorganic N from messy eating
- Total potential immobilization and gross mineralization
- Potential heterotrophic respiration (for methane calculations)

Corresponds to `SoilBiogeochemPotential` in `SoilBiogeochemPotentialMod.F90`.

# Arguments
- `cf::SoilBiogeochemCarbonFluxData`: carbon flux data (in/out)
- `cs::SoilBiogeochemCarbonStateData`: carbon state data (in)
- `nf::SoilBiogeochemNitrogenFluxData`: nitrogen flux data (out: potential_immob_vr_col, gross_nmin_vr_col)
- `ns::SoilBiogeochemNitrogenStateData`: nitrogen state data (in)
- `st::SoilBiogeochemStateData`: biogeochem state data (in: nue_decomp_cascade_col)
- `cascade_con::DecompCascadeConData`: cascade configuration
- `mask_bgc_soilc::BitVector`: mask for BGC soil columns
- `bounds::UnitRange{Int}`: column bounds
- `nlevdecomp::Int`: number of decomposition levels
- `ndecomp_pools::Int`: number of decomposition pools
- `ndecomp_cascade_transitions::Int`: number of cascade transitions
- `cn_decomp_pools::Array{<:Real,3}`: (out) C:N ratios of pools [col x lev x pool]
- `p_decomp_cpool_loss::Array{<:Real,3}`: (out) potential C loss [col x lev x trans]
- `p_decomp_cn_gain::Array{<:Real,3}`: (out) C:N of flux gained by receiver [col x lev x pool]
- `pmnf_decomp_cascade::Array{<:Real,3}`: (out) potential mineral N flux [col x lev x trans]
- `p_decomp_npool_to_din::Array{<:Real,3}`: (out) potential flux to DIN [col x lev x trans]
- `use_mimics::Bool`: whether MIMICS decomposition is active
- `i_cop_mic::Int`: index of copiotrophic microbe pool (only used when use_mimics=true)
- `i_oli_mic::Int`: index of oligotrophic microbe pool (only used when use_mimics=true)
"""
function soil_bgc_potential!(
        cf::SoilBiogeochemCarbonFluxData,
        cs::SoilBiogeochemCarbonStateData,
        nf::SoilBiogeochemNitrogenFluxData,
        ns::SoilBiogeochemNitrogenStateData,
        st::SoilBiogeochemStateData,
        cascade_con::DecompCascadeConData;
        mask_bgc_soilc::BitVector,
        bounds::UnitRange{Int},
        nlevdecomp::Int,
        ndecomp_pools::Int,
        ndecomp_cascade_transitions::Int,
        cn_decomp_pools::Array{<:Real,3},
        p_decomp_cpool_loss::Array{<:Real,3},
        p_decomp_cn_gain::Array{<:Real,3},
        pmnf_decomp_cascade::Array{<:Real,3},
        p_decomp_npool_to_din::Array{<:Real,3},
        use_mimics::Bool=false,
        i_cop_mic::Int=0,
        i_oli_mic::Int=0)

    nc = length(bounds)

    # Unpack cascade configuration
    cascade_donor_pool             = cascade_con.cascade_donor_pool
    cascade_receiver_pool          = cascade_con.cascade_receiver_pool
    floating_cn_ratio_decomp_pools = cascade_con.floating_cn_ratio_decomp_pools
    initial_cn_ratio               = cascade_con.initial_cn_ratio

    # Unpack state
    nue_decomp_cascade = st.nue_decomp_cascade_col

    # Unpack carbon flux/state
    rf_decomp_cascade       = cf.rf_decomp_cascade_col
    pathfrac_decomp_cascade = cf.pathfrac_decomp_cascade_col
    decomp_k                = cf.decomp_k_col
    cn_col                  = cf.cn_col
    phr_vr                  = cf.phr_vr_col

    # Unpack carbon state
    decomp_cpools_vr = cs.decomp_cpools_vr_col

    # Unpack nitrogen state/flux
    decomp_npools_vr = ns.decomp_npools_vr_col
    potential_immob_vr = nf.potential_immob_vr_col
    gross_nmin_vr      = nf.gross_nmin_vr_col

    # Local arrays for MIMICS pathway
    FT = eltype(decomp_cpools_vr)
    p_decomp_cpool_gain = zeros(FT, nc, nlevdecomp, ndecomp_cascade_transitions)
    p_decomp_npool_gain = zeros(FT, nc, nlevdecomp, ndecomp_cascade_transitions)
    p_decomp_cpool_gain_sum = zeros(FT, ndecomp_pools)
    p_decomp_npool_gain_sum = zeros(FT, ndecomp_pools)

    # Local immobilization accumulator
    immob = zeros(FT, nc, nlevdecomp)

    # -------------------------------------------------------------------
    # Set initial values for potential C and N fluxes
    # -------------------------------------------------------------------
    for k in 1:ndecomp_cascade_transitions
        for j in 1:nlevdecomp
            for c in 1:nc
                p_decomp_cpool_loss[c, j, k] = 0.0
                pmnf_decomp_cascade[c, j, k] = 0.0
            end
        end
    end

    # -------------------------------------------------------------------
    # Calculate C:N ratios of applicable pools
    # -------------------------------------------------------------------
    for l in 1:ndecomp_pools
        if floating_cn_ratio_decomp_pools[l]
            for j in 1:nlevdecomp
                for c in 1:nc
                    mask_bgc_soilc[c] || continue
                    if decomp_npools_vr[c, j, l] > 0.0
                        cn_decomp_pools[c, j, l] = decomp_cpools_vr[c, j, l] / decomp_npools_vr[c, j, l]
                    end
                end
            end
        else
            for j in 1:nlevdecomp
                for c in 1:nc
                    mask_bgc_soilc[c] || continue
                    cn_decomp_pools[c, j, l] = initial_cn_ratio[l]
                end
            end
        end
    end

    # -------------------------------------------------------------------
    # Calculate the non-nitrogen-limited fluxes
    # These fluxes include the "/ dt" term to put them on a per second
    # basis, since the rate constants have been calculated on a per
    # timestep basis.
    # -------------------------------------------------------------------
    for k in 1:ndecomp_cascade_transitions
        for j in 1:nlevdecomp
            for c in 1:nc
                mask_bgc_soilc[c] || continue

                if decomp_cpools_vr[c, j, cascade_donor_pool[k]] > 0.0 &&
                   decomp_k[c, j, cascade_donor_pool[k]] > 0.0

                    p_decomp_cpool_loss[c, j, k] = decomp_cpools_vr[c, j, cascade_donor_pool[k]] *
                        decomp_k[c, j, cascade_donor_pool[k]] * pathfrac_decomp_cascade[c, j, k]

                    # Determine if receiver pool has floating C:N ratio.
                    # I_ATM (= 0) is not a valid pool index; treat it as non-floating.
                    recv_pool = cascade_receiver_pool[k]
                    recv_floating = (recv_pool != I_ATM) && floating_cn_ratio_decomp_pools[recv_pool]

                    if !recv_floating
                        # Not transition of CWD to litter (fixed C:N receiver or atmosphere)

                        if cascade_receiver_pool[k] != I_ATM
                            # Not 100% respiration
                            ratio = 0.0

                            if decomp_npools_vr[c, j, cascade_donor_pool[k]] > 0.0
                                ratio = cn_decomp_pools[c, j, cascade_receiver_pool[k]] /
                                        cn_decomp_pools[c, j, cascade_donor_pool[k]]
                            end

                            pmnf_decomp_cascade[c, j, k] = (p_decomp_cpool_loss[c, j, k] *
                                (1.0 - rf_decomp_cascade[c, j, k] - ratio) /
                                cn_decomp_pools[c, j, cascade_receiver_pool[k]])

                        else
                            # 100% respiration (receiver is atmosphere)
                            pmnf_decomp_cascade[c, j, k] = -p_decomp_cpool_loss[c, j, k] /
                                cn_decomp_pools[c, j, cascade_donor_pool[k]]
                        end

                    else
                        # CWD -> litter OR mimics_decomp is true (floating C:N receiver)
                        pmnf_decomp_cascade[c, j, k] = 0.0

                        if use_mimics
                            # N:C ratio of donor pools (N:C instead of C:N because
                            # already checked that we're not dividing by zero)
                            decomp_nc_loss_donor =
                                decomp_npools_vr[c, j, cascade_donor_pool[k]] /
                                decomp_cpools_vr[c, j, cascade_donor_pool[k]]

                            # Calculate N fluxes from donor pools
                            p_decomp_npool_loss =
                                p_decomp_cpool_loss[c, j, k] * decomp_nc_loss_donor

                            # Track N lost to DIN from messy eating
                            p_decomp_cpool_gain[c, j, k] =
                                p_decomp_cpool_loss[c, j, k] * (1.0 - rf_decomp_cascade[c, j, k])

                            p_decomp_npool_gain[c, j, k] =
                                p_decomp_npool_loss * nue_decomp_cascade[k]

                            p_decomp_npool_to_din[c, j, k] =
                                p_decomp_npool_loss - p_decomp_npool_gain[c, j, k]
                        end  # use_mimics

                    end  # floating_cn_ratio check

                else
                    # Donors not donating (decomp_cpools_vr or decomp_k <= 0)
                    p_decomp_cpool_gain[c, j, k] = 0.0
                    p_decomp_npool_gain[c, j, k] = 0.0
                    p_decomp_npool_to_din[c, j, k] = 0.0
                end  # donors donating check
            end  # column loop
        end  # soil level loop
    end  # transitions loop

    # -------------------------------------------------------------------
    # MIMICS: Calculate cn_gain into microbial biomass and determine
    # immobilization vs. mineralization
    # -------------------------------------------------------------------
    if use_mimics
        for j in 1:nlevdecomp
            for c in 1:nc
                mask_bgc_soilc[c] || continue

                # Sum C & N fluxes from all transitions into m1 & m2 pools
                p_decomp_cpool_gain_sum .= 0.0
                p_decomp_npool_gain_sum .= 0.0

                for k in 1:ndecomp_cascade_transitions
                    if cascade_receiver_pool[k] == i_cop_mic ||
                       cascade_receiver_pool[k] == i_oli_mic
                        if decomp_cpools_vr[c, j, cascade_donor_pool[k]] > 0.0 &&
                           decomp_k[c, j, cascade_donor_pool[k]] > 0.0

                            rp = cascade_receiver_pool[k]
                            p_decomp_cpool_gain_sum[rp] += p_decomp_cpool_gain[c, j, k]
                            p_decomp_npool_gain_sum[rp] += p_decomp_npool_gain[c, j, k]

                            if p_decomp_npool_gain_sum[rp] > 0.0
                                p_decomp_cn_gain[c, j, rp] =
                                    p_decomp_cpool_gain_sum[rp] / p_decomp_npool_gain_sum[rp]
                            else
                                p_decomp_cn_gain[c, j, rp] = 0.0
                            end
                        end
                    end
                end  # first transitions loop

                for k in 1:ndecomp_cascade_transitions
                    if cascade_receiver_pool[k] == i_cop_mic ||
                       cascade_receiver_pool[k] == i_oli_mic
                        if decomp_cpools_vr[c, j, cascade_donor_pool[k]] > 0.0 &&
                           decomp_k[c, j, cascade_donor_pool[k]] > 0.0

                            rp = cascade_receiver_pool[k]

                            # if p_decomp_cn_diff < 0 => N mineralization
                            #                     > 0 => immobilization
                            # "min" turns off immobilization flux
                            p_decomp_cn_diff_ratio = min(0.0,
                                (p_decomp_cn_gain[c, j, rp] - cn_col[c, rp]) / cn_col[c, rp])

                            # Actual amount of N mineralized or to be immobilized
                            # negative=mineralization, positive=immobilization
                            pmnf_decomp_cascade[c, j, k] = p_decomp_cn_diff_ratio * p_decomp_npool_gain[c, j, k]
                        end
                    end
                end  # second transitions loop
            end  # column loop
        end  # soil levels loop
    end  # use_mimics

    # -------------------------------------------------------------------
    # Sum up all potential immobilization fluxes (positive pmnf)
    # and all mineralization fluxes (negative pmnf)
    # -------------------------------------------------------------------
    for j in 1:nlevdecomp
        for c in 1:nc
            mask_bgc_soilc[c] || continue
            immob[c, j] = 0.0
        end
    end

    for k in 1:ndecomp_cascade_transitions
        for j in 1:nlevdecomp
            for c in 1:nc
                mask_bgc_soilc[c] || continue
                if pmnf_decomp_cascade[c, j, k] > 0.0
                    immob[c, j] += pmnf_decomp_cascade[c, j, k]
                else
                    gross_nmin_vr[c, j] -= pmnf_decomp_cascade[c, j, k]
                end
                if use_mimics
                    gross_nmin_vr[c, j] += p_decomp_npool_to_din[c, j, k]
                end
            end
        end
    end

    for j in 1:nlevdecomp
        for c in 1:nc
            mask_bgc_soilc[c] || continue
            potential_immob_vr[c, j] = immob[c, j]
        end
    end

    # -------------------------------------------------------------------
    # Add up potential HR for methane calculations
    # -------------------------------------------------------------------
    for j in 1:nlevdecomp
        for c in 1:nc
            mask_bgc_soilc[c] || continue
            phr_vr[c, j] = 0.0
        end
    end

    for k in 1:ndecomp_cascade_transitions
        for j in 1:nlevdecomp
            for c in 1:nc
                mask_bgc_soilc[c] || continue
                phr_vr[c, j] += rf_decomp_cascade[c, j, k] * p_decomp_cpool_loss[c, j, k]
            end
        end
    end

    return nothing
end
