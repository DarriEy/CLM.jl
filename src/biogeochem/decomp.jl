# ==========================================================================
# Ported from: src/soilbiogeochem/SoilBiogeochemDecompMod.F90
# Module holding routines used in litter and soil decomposition model.
#
# Public functions:
#   decomp_read_params!          — Read decomposition parameters (dnp)
#   soil_biogeochem_decomp!      — Main decomposition routine
# ==========================================================================

# ---------------------------------------------------------------------------
# DecompParams — decomposition parameters
# Ported from params_type in SoilBiogeochemDecompMod.F90
# ---------------------------------------------------------------------------

Base.@kwdef mutable struct DecompParams
    dnp::Float64 = 0.01  # denitrification proportion
end

# Constant: atmosphere pool index (receiver_pool == 0 means "to atmosphere")
# Ported from i_atm in SoilBiogeochemDecompCascadeConType.F90
const I_ATM = 0

# ---------------------------------------------------------------------------
# decomp_read_params! — Read decomposition parameters
# Ported from readParams in SoilBiogeochemDecompMod.F90
# ---------------------------------------------------------------------------

function decomp_read_params!(params::DecompParams; dnp::Float64)
    params.dnp = dnp
    return nothing
end

# ---------------------------------------------------------------------------
# soil_biogeochem_decomp! — Main decomposition subroutine
# Ported from SoilBiogeochemDecomp in SoilBiogeochemDecompMod.F90
#
# Calculates actual immobilization and decomp rates, following resolution
# of plant/heterotroph competition for mineral N.
# ---------------------------------------------------------------------------

function soil_biogeochem_decomp!(
        cf::SoilBiogeochemCarbonFluxData,
        cs::SoilBiogeochemCarbonStateData,
        nf::SoilBiogeochemNitrogenFluxData,
        ns::SoilBiogeochemNitrogenStateData,
        st::SoilBiogeochemStateData,
        cascade_con::DecompCascadeConData,
        params::DecompParams;
        mask_bgc_soilc::BitVector,
        bounds::UnitRange{Int},
        nlevdecomp::Int,
        ndecomp_pools::Int,
        ndecomp_cascade_transitions::Int,
        cn_decomp_pools::Array{Float64,3},
        p_decomp_cpool_loss::Array{Float64,3},
        pmnf_decomp_cascade::Array{Float64,3},
        p_decomp_npool_to_din::Array{Float64,3},
        dzsoi_decomp::Vector{Float64},
        use_nitrif_denitrif::Bool=false,
        use_lch4::Bool=false,
        use_mimics::Bool=false,
        use_soil_matrixcn::Bool=false)

    nc = length(bounds)

    # Unpack cascade configuration
    cascade_donor_pool             = cascade_con.cascade_donor_pool
    cascade_receiver_pool          = cascade_con.cascade_receiver_pool
    floating_cn_ratio_decomp_pools = cascade_con.floating_cn_ratio_decomp_pools
    initial_cn_ratio               = cascade_con.initial_cn_ratio

    # Unpack state
    fpi_vr = st.fpi_vr_col

    # Unpack carbon flux/state
    rf_decomp_cascade           = cf.rf_decomp_cascade_col
    decomp_cpools_vr            = cs.decomp_cpools_vr_col
    c_overflow_vr               = cf.c_overflow_vr
    w_scalar                    = cf.w_scalar_col
    phr_vr                      = cf.phr_vr_col
    fphr                        = cf.fphr_col
    decomp_cascade_hr_vr        = cf.decomp_cascade_hr_vr_col
    decomp_cascade_ctransfer_vr = cf.decomp_cascade_ctransfer_vr_col

    # Unpack nitrogen flux/state
    decomp_npools_vr                     = ns.decomp_npools_vr_col
    decomp_cascade_ntransfer_vr          = nf.decomp_cascade_ntransfer_vr_col
    decomp_cascade_sminn_flux_vr         = nf.decomp_cascade_sminn_flux_vr_col
    sminn_to_denit_decomp_cascade_vr     = nf.sminn_to_denit_decomp_cascade_vr_col
    gross_nmin_vr                        = nf.gross_nmin_vr_col
    net_nmin_vr                          = nf.net_nmin_vr_col
    gross_nmin                           = nf.gross_nmin_col
    net_nmin                             = nf.net_nmin_col

    # -------------------------------------------------------------------
    # Calculate c:n ratios of applicable pools
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
    # Calculate actual immobilization and decomp rates, following
    # resolution of plant/heterotroph competition for mineral N.
    #
    # Upon return from SoilBiogeochemCompetition, the fraction of potential
    # immobilization has been set (fpi_vr). Now finish the decomp calcs.
    # Only the immobilization steps are limited by fpi_vr (pmnf > 0).
    # Also calculate denitrification losses as a simple proportion
    # of mineralization flux.
    # -------------------------------------------------------------------
    for k in 1:ndecomp_cascade_transitions
        for j in 1:nlevdecomp
            for c in 1:nc
                mask_bgc_soilc[c] || continue

                if decomp_cpools_vr[c, j, cascade_donor_pool[k]] > 0.0
                    if pmnf_decomp_cascade[c, j, k] > 0.0
                        p_decomp_cpool_loss[c, j, k] *= fpi_vr[c, j]
                        pmnf_decomp_cascade[c, j, k] *= fpi_vr[c, j]
                        # use_soil_matrixcn Ksoil%DM correction not yet implemented
                        if !use_nitrif_denitrif
                            sminn_to_denit_decomp_cascade_vr[c, j, k] = 0.0
                        end
                    else
                        if !use_nitrif_denitrif
                            sminn_to_denit_decomp_cascade_vr[c, j, k] = -params.dnp * pmnf_decomp_cascade[c, j, k]
                        end
                    end

                    decomp_cascade_hr_vr[c, j, k] = rf_decomp_cascade[c, j, k] * p_decomp_cpool_loss[c, j, k]
                    decomp_cascade_ctransfer_vr[c, j, k] = (1.0 - rf_decomp_cascade[c, j, k]) * p_decomp_cpool_loss[c, j, k]

                    if use_mimics
                        decomp_cascade_hr_vr[c, j, k] = min(
                            p_decomp_cpool_loss[c, j, k],
                            decomp_cascade_hr_vr[c, j, k] + c_overflow_vr[c, j, k])
                        decomp_cascade_ctransfer_vr[c, j, k] = max(0.0, p_decomp_cpool_loss[c, j, k] - decomp_cascade_hr_vr[c, j, k])
                    end

                    if decomp_npools_vr[c, j, cascade_donor_pool[k]] > 0.0 && cascade_receiver_pool[k] != I_ATM
                        decomp_cascade_ntransfer_vr[c, j, k] = p_decomp_cpool_loss[c, j, k] / cn_decomp_pools[c, j, cascade_donor_pool[k]]
                    else
                        decomp_cascade_ntransfer_vr[c, j, k] = 0.0
                    end

                    if cascade_receiver_pool[k] != 0
                        decomp_cascade_sminn_flux_vr[c, j, k] = pmnf_decomp_cascade[c, j, k]
                    else  # keep sign convention negative for terminal pools
                        decomp_cascade_sminn_flux_vr[c, j, k] = -pmnf_decomp_cascade[c, j, k]
                    end

                    net_nmin_vr[c, j] -= pmnf_decomp_cascade[c, j, k]

                    if use_mimics
                        decomp_cascade_sminn_flux_vr[c, j, k] -= p_decomp_npool_to_din[c, j, k]
                        net_nmin_vr[c, j] += p_decomp_npool_to_din[c, j, k]
                    end
                else
                    decomp_cascade_ntransfer_vr[c, j, k] = 0.0
                    if !use_nitrif_denitrif
                        sminn_to_denit_decomp_cascade_vr[c, j, k] = 0.0
                    end
                    decomp_cascade_sminn_flux_vr[c, j, k] = 0.0
                end
            end
        end
    end

    # -------------------------------------------------------------------
    # Calculate total fraction of potential HR, for methane code
    # -------------------------------------------------------------------
    if use_lch4
        hrsum = zeros(nc, nlevdecomp)

        for k in 1:ndecomp_cascade_transitions
            for j in 1:nlevdecomp
                for c in 1:nc
                    mask_bgc_soilc[c] || continue
                    hrsum[c, j] += rf_decomp_cascade[c, j, k] * p_decomp_cpool_loss[c, j, k]
                end
            end
        end

        # Nitrogen limitation / (low)-moisture limitation
        for j in 1:nlevdecomp
            for c in 1:nc
                mask_bgc_soilc[c] || continue
                if phr_vr[c, j] > 0.0
                    fphr[c, j] = hrsum[c, j] / phr_vr[c, j] * w_scalar[c, j]
                    fphr[c, j] = max(fphr[c, j], 0.01)  # Prevent overflow errors for 0 respiration
                else
                    fphr[c, j] = 1.0
                end
            end
        end
    end

    # -------------------------------------------------------------------
    # Vertically integrate net and gross mineralization fluxes
    # for diagnostic output
    # -------------------------------------------------------------------
    for c in 1:nc
        mask_bgc_soilc[c] || continue
        for j in 1:nlevdecomp
            net_nmin[c]   += net_nmin_vr[c, j] * dzsoi_decomp[j]
            gross_nmin[c] += gross_nmin_vr[c, j] * dzsoi_decomp[j]
        end
    end

    return nothing
end
