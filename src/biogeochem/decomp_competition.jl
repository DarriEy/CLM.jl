# ==========================================================================
# Ported from: src/soilbiogeochem/SoilBiogeochemCompetitionMod.F90
# Resolve plant/heterotroph competition for mineral N
#
# Public functions:
#   soil_bgc_competition_params_read!  -- Read competition parameters
#   soil_bgc_competition_init!         -- Initialize module-level state
#   soil_bgc_competition!             -- Main competition routine
# ==========================================================================

# ---------------------------------------------------------------------------
# SoilBGCCompetitionParams -- competition parameters read from file
# Ported from params_type in SoilBiogeochemCompetitionMod.F90
# ---------------------------------------------------------------------------

"""
    SoilBGCCompetitionParams

Parameters for soil N competition. Holds bulk denitrification rate and
relative competitiveness factors for plants, decomposers, nitrifiers,
and denitrifiers.

Ported from `params_type` in `SoilBiogeochemCompetitionMod.F90`.
"""
Base.@kwdef mutable struct SoilBGCCompetitionParams
    bdnr              ::Float64 = 0.5   # bulk denitrification rate (1/s)
    compet_plant_no3  ::Float64 = 1.0   # (unitless) relative competitiveness of plants for NO3
    compet_plant_nh4  ::Float64 = 1.0   # (unitless) relative competitiveness of plants for NH4
    compet_decomp_no3 ::Float64 = 1.0   # (unitless) relative competitiveness of immobilizers for NO3
    compet_decomp_nh4 ::Float64 = 1.0   # (unitless) relative competitiveness of immobilizers for NH4
    compet_denit      ::Float64 = 1.0   # (unitless) relative competitiveness of denitrifiers for NO3
    compet_nit        ::Float64 = 1.0   # (unitless) relative competitiveness of nitrifiers for NH4
end

# ---------------------------------------------------------------------------
# SoilBGCCompetitionState -- module-level state
# Holds timestep-scaled bdnr, carbon_only flag, and supplemental N mode
# ---------------------------------------------------------------------------

"""
    SoilBGCCompetitionState

Module-level persistent state for soil N competition.
Holds the timestep-scaled bulk denitrification rate and carbon-only flag.

Ported from module-level variables in `SoilBiogeochemCompetitionMod.F90`.
"""
Base.@kwdef mutable struct SoilBGCCompetitionState
    dt            ::Float64 = 1800.0   # decomp timestep (seconds)
    bdnr          ::Float64 = 0.0      # bulk denitrification rate scaled to timestep
    carbon_only   ::Bool    = false    # if true, supplement N to eliminate N limitation
end

# --- Supplemental nitrogen mode constants ---
const SUPLN_ALL  = "ALL"
const SUPLN_NONE = "NONE"

# ---------------------------------------------------------------------------
# soil_bgc_competition_params_read! -- Read competition parameters
# Ported from readParams in SoilBiogeochemCompetitionMod.F90
# ---------------------------------------------------------------------------

"""
    soil_bgc_competition_params_read!(params; bdnr, compet_plant_no3,
        compet_plant_nh4, compet_decomp_no3, compet_decomp_nh4,
        compet_denit, compet_nit)

Read competition parameters from keyword arguments (replaces NetCDF file reading).
Corresponds to `readParams` in the Fortran source.
"""
function soil_bgc_competition_params_read!(params::SoilBGCCompetitionParams;
                                            bdnr::Float64,
                                            compet_plant_no3::Float64,
                                            compet_plant_nh4::Float64,
                                            compet_decomp_no3::Float64,
                                            compet_decomp_nh4::Float64,
                                            compet_denit::Float64,
                                            compet_nit::Float64)
    params.bdnr              = bdnr
    params.compet_plant_no3  = compet_plant_no3
    params.compet_plant_nh4  = compet_plant_nh4
    params.compet_decomp_no3 = compet_decomp_no3
    params.compet_decomp_nh4 = compet_decomp_nh4
    params.compet_denit      = compet_denit
    params.compet_nit        = compet_nit
    return nothing
end

# ---------------------------------------------------------------------------
# soil_bgc_competition_init! -- Initialization
# Ported from SoilBiogeochemCompetitionInit in SoilBiogeochemCompetitionMod.F90
# ---------------------------------------------------------------------------

"""
    soil_bgc_competition_init!(state, params; dt, secspday, suplnitro)

Initialize the competition module: scale bdnr to the decomp timestep and
set the carbon_only flag based on the supplemental nitrogen mode.

Corresponds to `SoilBiogeochemCompetitionInit` in the Fortran source.
"""
function soil_bgc_competition_init!(state::SoilBGCCompetitionState,
                                     params::SoilBGCCompetitionParams;
                                     dt::Float64,
                                     secspday::Float64=SECSPDAY,
                                     suplnitro::String=SUPLN_NONE)
    state.dt   = dt
    state.bdnr = params.bdnr * (dt / secspday)

    if suplnitro == SUPLN_NONE
        state.carbon_only = false
    elseif suplnitro == SUPLN_ALL
        state.carbon_only = true
    else
        error("Supplemental Nitrogen flag (suplnitro) can only be: " *
              SUPLN_NONE * " or " * SUPLN_ALL *
              ", got: " * suplnitro)
    end

    return nothing
end

# ---------------------------------------------------------------------------
# soil_bgc_competition! -- Main competition routine
# Ported from SoilBiogeochemCompetition in SoilBiogeochemCompetitionMod.F90
#
# Resolves plant/heterotroph (and optionally nitrifier/denitrifier)
# competition for mineral N.
# ---------------------------------------------------------------------------

"""
    soil_bgc_competition!(st, nf, cf, ns, state, params; ...)

Resolve plant/heterotroph/nitrifier/denitrifier competition for mineral N.

# Arguments
- `st::SoilBiogeochemStateData` -- soil BGC state (fpg, fpi, fpi_vr, nfixation_prof, plant_ndemand)
- `nf::SoilBiogeochemNitrogenFluxData` -- N fluxes (immob, plant uptake, nitrif/denitrif)
- `cf::SoilBiogeochemCarbonFluxData` -- C fluxes (c_overflow_vr for MIMICS)
- `ns::SoilBiogeochemNitrogenStateData` -- N state (sminn_vr, smin_nh4_vr, smin_no3_vr)
- `state::SoilBGCCompetitionState` -- module-level state (dt, bdnr)
- `params::SoilBGCCompetitionParams` -- competition parameters (competitiveness factors)
- `mask_bgc_soilc::BitVector` -- mask for BGC soil columns
- `bounds::UnitRange{Int}` -- column bounds
- `nlevdecomp::Int` -- number of decomposition levels
- `ndecomp_cascade_transitions::Int` -- number of decomp cascade transitions
- `dzsoi_decomp::Vector{Float64}` -- decomposition level thicknesses
- `pmnf_decomp_cascade::Array{Float64,3}` -- potential mineral N flux from decomp
- `p_decomp_cn_gain::Array{Float64,3}` -- C:N ratio of flux gained by receiver pool
- `cascade_receiver_pool::Vector{Int}` -- receiver pool index for each transition
- `use_nitrif_denitrif::Bool` -- use nitrification/denitrification model
- `use_fun::Bool` -- use FUN model
- `carbon_only::Bool` -- carbon-only mode (supplement N)
- `mimics_decomp::Bool` -- using MIMICS decomposition
- `i_cop_mic::Int` -- copiotrophic microbe pool index (MIMICS)
- `i_oli_mic::Int` -- oligotrophic microbe pool index (MIMICS)
- `nitrif_n2o_loss_frac::Float64` -- fraction of nitrification lost as N2O

Corresponds to `SoilBiogeochemCompetition` in the Fortran source.
"""
function soil_bgc_competition!(
        st::SoilBiogeochemStateData,
        nf::SoilBiogeochemNitrogenFluxData,
        cf::SoilBiogeochemCarbonFluxData,
        ns::SoilBiogeochemNitrogenStateData,
        state::SoilBGCCompetitionState,
        params::SoilBGCCompetitionParams;
        mask_bgc_soilc::BitVector,
        bounds::UnitRange{Int},
        nlevdecomp::Int,
        ndecomp_cascade_transitions::Int,
        dzsoi_decomp::Vector{Float64},
        pmnf_decomp_cascade::Array{Float64,3},
        p_decomp_cn_gain::Array{Float64,3},
        cascade_receiver_pool::Vector{Int},
        use_nitrif_denitrif::Bool=false,
        use_fun::Bool=false,
        carbon_only::Bool=false,
        mimics_decomp::Bool=false,
        i_cop_mic::Int=0,
        i_oli_mic::Int=0,
        nitrif_n2o_loss_frac::Float64=NITRIF_N2O_LOSS_FRAC)

    dt   = state.dt
    bdnr = state.bdnr

    # Unpack state variables
    fpg              = st.fpg_col
    fpi              = st.fpi_col
    fpi_vr           = st.fpi_vr_col
    nfixation_prof   = st.nfixation_prof_col
    plant_ndemand    = st.plant_ndemand_col

    # Unpack nitrogen state
    sminn_vr     = ns.sminn_vr_col
    smin_nh4_vr  = ns.smin_nh4_vr_col
    smin_no3_vr  = ns.smin_no3_vr_col

    # Unpack nitrogen flux
    potential_immob_vr           = nf.potential_immob_vr_col
    actual_immob_vr              = nf.actual_immob_vr_col
    sminn_to_plant_vr            = nf.sminn_to_plant_vr_col
    sminn_to_plant               = nf.sminn_to_plant_col
    actual_immob                 = nf.actual_immob_col
    potential_immob              = nf.potential_immob_col
    supplement_to_sminn_vr       = nf.supplement_to_sminn_vr_col
    sminn_to_denit_excess_vr     = nf.sminn_to_denit_excess_vr_col
    actual_immob_no3_vr          = nf.actual_immob_no3_vr_col
    actual_immob_nh4_vr          = nf.actual_immob_nh4_vr_col
    smin_no3_to_plant_vr         = nf.smin_no3_to_plant_vr_col
    smin_nh4_to_plant_vr         = nf.smin_nh4_to_plant_vr_col
    pot_f_nit_vr                 = nf.pot_f_nit_vr_col
    pot_f_denit_vr               = nf.pot_f_denit_vr_col
    f_nit_vr                     = nf.f_nit_vr_col
    f_denit_vr                   = nf.f_denit_vr_col
    n2_n2o_ratio_denit_vr        = nf.n2_n2o_ratio_denit_vr_col
    f_n2o_denit_vr               = nf.f_n2o_denit_vr_col
    f_n2o_nit_vr                 = nf.f_n2o_nit_vr_col
    sminn_to_plant_fun_vr        = nf.sminn_to_plant_fun_vr_col
    sminn_to_plant_fun_no3_vr    = nf.sminn_to_plant_fun_no3_vr_col
    sminn_to_plant_fun_nh4_vr    = nf.sminn_to_plant_fun_nh4_vr_col

    # Unpack carbon flux (for MIMICS)
    c_overflow_vr = cf.c_overflow_vr

    local_use_fun = use_fun

    # Local arrays
    n_arr = length(mask_bgc_soilc)
    sminn_tot              = zeros(n_arr)
    sminn_to_plant_new     = zeros(n_arr)
    nuptake_prof           = zeros(n_arr, nlevdecomp)
    sum_ndemand_vr         = zeros(n_arr, nlevdecomp)

    if !use_nitrif_denitrif
        # ====================================================================
        # NON-NITRIF_DENITRIF PATHWAY
        # ====================================================================

        nlimit             = zeros(Int, n_arr, nlevdecomp)
        residual_sminn_vr  = zeros(n_arr, nlevdecomp)
        residual_sminn     = zeros(n_arr)
        residual_plant_ndemand = zeros(n_arr)

        # init sminn_tot
        for c in bounds
            mask_bgc_soilc[c] || continue
            sminn_tot[c] = 0.0
        end

        # sum total mineral N
        for j in 1:nlevdecomp
            for c in bounds
                mask_bgc_soilc[c] || continue
                sminn_tot[c] += sminn_vr[c, j] * dzsoi_decomp[j]
            end
        end

        # define N uptake profile
        for j in 1:nlevdecomp
            for c in bounds
                mask_bgc_soilc[c] || continue
                if sminn_tot[c] > 0.0
                    nuptake_prof[c, j] = sminn_vr[c, j] / sminn_tot[c]
                else
                    nuptake_prof[c, j] = nfixation_prof[c, j]
                end
            end
        end

        # total N demand at each level
        for j in 1:nlevdecomp
            for c in bounds
                mask_bgc_soilc[c] || continue
                sum_ndemand_vr[c, j] = plant_ndemand[c] * nuptake_prof[c, j] + potential_immob_vr[c, j]
            end
        end

        # resolve competition at each level
        for j in 1:nlevdecomp
            for c in bounds
                mask_bgc_soilc[c] || continue

                if sum_ndemand_vr[c, j] * dt < sminn_vr[c, j]
                    # N availability is not limiting
                    nlimit[c, j] = 0
                    fpi_vr[c, j] = 1.0
                    actual_immob_vr[c, j] = potential_immob_vr[c, j]
                    sminn_to_plant_vr[c, j] = plant_ndemand[c] * nuptake_prof[c, j]

                elseif carbon_only
                    # Carbon-only mode: supplement N to eliminate limitation
                    nlimit[c, j] = 1
                    fpi_vr[c, j] = 1.0
                    actual_immob_vr[c, j] = potential_immob_vr[c, j]
                    sminn_to_plant_vr[c, j] = plant_ndemand[c] * nuptake_prof[c, j]
                    supplement_to_sminn_vr[c, j] = sum_ndemand_vr[c, j] - (sminn_vr[c, j] / dt)

                else
                    # N limited: competition between plant and decomposers
                    nlimit[c, j] = 1
                    if sum_ndemand_vr[c, j] > 0.0
                        actual_immob_vr[c, j] = (sminn_vr[c, j] / dt) * (potential_immob_vr[c, j] / sum_ndemand_vr[c, j])
                    else
                        actual_immob_vr[c, j] = 0.0
                    end

                    if potential_immob_vr[c, j] > 0.0
                        fpi_vr[c, j] = actual_immob_vr[c, j] / potential_immob_vr[c, j]
                    else
                        fpi_vr[c, j] = 0.0
                    end

                    sminn_to_plant_vr[c, j] = (sminn_vr[c, j] / dt) - actual_immob_vr[c, j]
                end
            end
        end

        # sum up N fluxes to plant
        for j in 1:nlevdecomp
            for c in bounds
                mask_bgc_soilc[c] || continue
                sminn_to_plant[c] += sminn_to_plant_vr[c, j] * dzsoi_decomp[j]
                if local_use_fun
                    if sminn_to_plant_fun_vr[c, j] > sminn_to_plant_vr[c, j]
                        sminn_to_plant_fun_vr[c, j] = sminn_to_plant_vr[c, j]
                    end
                end
            end
        end

        # --- Second pass: distribute residual N to plants ---

        for c in bounds
            mask_bgc_soilc[c] || continue
            residual_sminn[c] = 0.0
            residual_plant_ndemand[c] = plant_ndemand[c] - sminn_to_plant[c]
        end

        for j in 1:nlevdecomp
            for c in bounds
                mask_bgc_soilc[c] || continue
                if residual_plant_ndemand[c] > 0.0
                    if nlimit[c, j] == 0
                        residual_sminn_vr[c, j] = max(sminn_vr[c, j] - (actual_immob_vr[c, j] + sminn_to_plant_vr[c, j]) * dt, 0.0)
                        residual_sminn[c] += residual_sminn_vr[c, j] * dzsoi_decomp[j]
                    else
                        residual_sminn_vr[c, j] = 0.0
                    end
                end
            end
        end

        # distribute residual N to plants
        for j in 1:nlevdecomp
            for c in bounds
                mask_bgc_soilc[c] || continue
                if residual_plant_ndemand[c] > 0.0 && residual_sminn[c] > 0.0 && nlimit[c, j] == 0
                    sminn_to_plant_vr[c, j] += residual_sminn_vr[c, j] *
                        min((residual_plant_ndemand[c] * dt) / residual_sminn[c], 1.0) / dt
                end
            end
        end

        # re-sum up N fluxes to plant
        for c in bounds
            mask_bgc_soilc[c] || continue
            sminn_to_plant[c] = 0.0
        end
        for j in 1:nlevdecomp
            for c in bounds
                mask_bgc_soilc[c] || continue
                sminn_to_plant[c] += sminn_to_plant_vr[c, j] * dzsoi_decomp[j]
                if !local_use_fun
                    sum_ndemand_vr[c, j] = potential_immob_vr[c, j] + sminn_to_plant_vr[c, j]
                else
                    sminn_to_plant_new[c] += sminn_to_plant_fun_vr[c, j] * dzsoi_decomp[j]
                    sum_ndemand_vr[c, j] = potential_immob_vr[c, j] + sminn_to_plant_fun_vr[c, j]
                end
            end
        end

        # excess denitrification
        for j in 1:nlevdecomp
            for c in bounds
                mask_bgc_soilc[c] || continue
                if !local_use_fun
                    if (sminn_to_plant_vr[c, j] + actual_immob_vr[c, j]) * dt < sminn_vr[c, j]
                        sminn_to_denit_excess_vr[c, j] = max(bdnr * ((sminn_vr[c, j] / dt) - sum_ndemand_vr[c, j]), 0.0)
                    else
                        sminn_to_denit_excess_vr[c, j] = 0.0
                    end
                else
                    if (sminn_to_plant_fun_vr[c, j] + actual_immob_vr[c, j]) * dt < sminn_vr[c, j]
                        sminn_to_denit_excess_vr[c, j] = max(bdnr * ((sminn_vr[c, j] / dt) - sum_ndemand_vr[c, j]), 0.0)
                    else
                        sminn_to_denit_excess_vr[c, j] = 0.0
                    end
                end
            end
        end

        # sum up N fluxes to immobilization
        for j in 1:nlevdecomp
            for c in bounds
                mask_bgc_soilc[c] || continue
                actual_immob[c] += actual_immob_vr[c, j] * dzsoi_decomp[j]
                potential_immob[c] += potential_immob_vr[c, j] * dzsoi_decomp[j]
            end
        end

        for c in bounds
            mask_bgc_soilc[c] || continue
            # fraction of potential growth achieved with available N
            if plant_ndemand[c] > 0.0
                if !local_use_fun
                    fpg[c] = sminn_to_plant[c] / plant_ndemand[c]
                else
                    fpg[c] = sminn_to_plant_new[c] / plant_ndemand[c]
                end
            else
                fpg[c] = 1.0
            end

            # fraction of immobilization realized
            if potential_immob[c] > 0.0
                fpi[c] = actual_immob[c] / potential_immob[c]
            else
                fpi[c] = 1.0
            end
        end

    else
        # ====================================================================
        # NITRIF_DENITRIF PATHWAY
        # ====================================================================

        # Competition parameters
        compet_plant_no3  = params.compet_plant_no3
        compet_plant_nh4  = params.compet_plant_nh4
        compet_decomp_no3 = params.compet_decomp_no3
        compet_decomp_nh4 = params.compet_decomp_nh4
        compet_denit      = params.compet_denit
        compet_nit        = params.compet_nit

        fpi_no3_vr_local       = zeros(n_arr, nlevdecomp)
        fpi_nh4_vr_local       = zeros(n_arr, nlevdecomp)
        sum_nh4_demand         = zeros(n_arr, nlevdecomp)
        sum_nh4_demand_scaled  = zeros(n_arr, nlevdecomp)
        sum_no3_demand         = zeros(n_arr, nlevdecomp)
        sum_no3_demand_scaled  = zeros(n_arr, nlevdecomp)
        nlimit_no3             = zeros(Int, n_arr, nlevdecomp)
        nlimit_nh4             = zeros(Int, n_arr, nlevdecomp)
        residual_plant_ndemand = zeros(n_arr)
        residual_smin_nh4_vr   = zeros(n_arr, nlevdecomp)
        residual_smin_no3_vr   = zeros(n_arr, nlevdecomp)
        residual_smin_nh4      = zeros(n_arr)
        residual_smin_no3      = zeros(n_arr)

        # init total mineral N pools
        for c in bounds
            mask_bgc_soilc[c] || continue
            sminn_tot[c] = 0.0
        end

        # sum up total mineral N pools
        for j in 1:nlevdecomp
            for c in bounds
                mask_bgc_soilc[c] || continue
                sminn_tot[c] += (smin_no3_vr[c, j] + smin_nh4_vr[c, j]) * dzsoi_decomp[j]
            end
        end

        # define N uptake profile
        for j in 1:nlevdecomp
            for c in bounds
                mask_bgc_soilc[c] || continue
                if sminn_tot[c] > 0.0
                    nuptake_prof[c, j] = sminn_vr[c, j] / sminn_tot[c]
                else
                    nuptake_prof[c, j] = nfixation_prof[c, j]
                end
            end
        end

        # main column/vertical loop
        for j in 1:nlevdecomp
            for c in bounds
                mask_bgc_soilc[c] || continue

                # --- First compete for NH4 ---
                sum_nh4_demand[c, j] = plant_ndemand[c] * nuptake_prof[c, j] +
                    potential_immob_vr[c, j] + pot_f_nit_vr[c, j]
                sum_nh4_demand_scaled[c, j] = plant_ndemand[c] * nuptake_prof[c, j] * compet_plant_nh4 +
                    potential_immob_vr[c, j] * compet_decomp_nh4 + pot_f_nit_vr[c, j] * compet_nit

                if sum_nh4_demand[c, j] * dt < smin_nh4_vr[c, j]
                    # NH4 not limiting
                    nlimit_nh4[c, j] = 0
                    fpi_nh4_vr_local[c, j] = 1.0
                    actual_immob_nh4_vr[c, j] = potential_immob_vr[c, j]
                    f_nit_vr[c, j] = pot_f_nit_vr[c, j]

                    if !local_use_fun
                        smin_nh4_to_plant_vr[c, j] = plant_ndemand[c] * nuptake_prof[c, j]
                    else
                        smin_nh4_to_plant_vr[c, j] = smin_nh4_vr[c, j] / dt -
                            actual_immob_nh4_vr[c, j] - f_nit_vr[c, j]
                    end
                else
                    # NH4 limited: competition
                    nlimit_nh4[c, j] = 1
                    if sum_nh4_demand[c, j] > 0.0
                        actual_immob_nh4_vr[c, j] = min(
                            (smin_nh4_vr[c, j] / dt) * (potential_immob_vr[c, j] *
                                compet_decomp_nh4 / sum_nh4_demand_scaled[c, j]),
                            potential_immob_vr[c, j])

                        f_nit_vr[c, j] = min(
                            (smin_nh4_vr[c, j] / dt) * (pot_f_nit_vr[c, j] * compet_nit /
                                sum_nh4_demand_scaled[c, j]),
                            pot_f_nit_vr[c, j])

                        if !local_use_fun
                            smin_nh4_to_plant_vr[c, j] = min(
                                (smin_nh4_vr[c, j] / dt) * (plant_ndemand[c] *
                                    nuptake_prof[c, j] * compet_plant_nh4 / sum_nh4_demand_scaled[c, j]),
                                plant_ndemand[c] * nuptake_prof[c, j])
                        else
                            smin_nh4_to_plant_vr[c, j] = smin_nh4_vr[c, j] / dt -
                                actual_immob_nh4_vr[c, j] - f_nit_vr[c, j]
                        end
                    else
                        actual_immob_nh4_vr[c, j] = 0.0
                        smin_nh4_to_plant_vr[c, j] = 0.0
                        f_nit_vr[c, j] = 0.0
                    end

                    if potential_immob_vr[c, j] > 0.0
                        fpi_nh4_vr_local[c, j] = actual_immob_nh4_vr[c, j] / potential_immob_vr[c, j]
                    else
                        fpi_nh4_vr_local[c, j] = 0.0
                    end
                end

                if mimics_decomp
                    fpi_nh4_vr_local[c, j] = 1.0
                    actual_immob_nh4_vr[c, j] = potential_immob_vr[c, j]
                end

                # --- Then compete for NO3 ---
                if !local_use_fun
                    sum_no3_demand[c, j] = (plant_ndemand[c] * nuptake_prof[c, j] -
                        smin_nh4_to_plant_vr[c, j]) +
                        (potential_immob_vr[c, j] - actual_immob_nh4_vr[c, j]) +
                        pot_f_denit_vr[c, j]
                    sum_no3_demand_scaled[c, j] = (plant_ndemand[c] * nuptake_prof[c, j] -
                        smin_nh4_to_plant_vr[c, j]) * compet_plant_no3 +
                        (potential_immob_vr[c, j] - actual_immob_nh4_vr[c, j]) * compet_decomp_no3 +
                        pot_f_denit_vr[c, j] * compet_denit
                else
                    sum_no3_demand[c, j] = plant_ndemand[c] * nuptake_prof[c, j] +
                        (potential_immob_vr[c, j] - actual_immob_nh4_vr[c, j]) +
                        pot_f_denit_vr[c, j]
                    sum_no3_demand_scaled[c, j] = (plant_ndemand[c] * nuptake_prof[c, j]) * compet_plant_no3 +
                        (potential_immob_vr[c, j] - actual_immob_nh4_vr[c, j]) * compet_decomp_no3 +
                        pot_f_denit_vr[c, j] * compet_denit
                end

                if sum_no3_demand[c, j] * dt < smin_no3_vr[c, j]
                    # NO3 not limiting
                    nlimit_no3[c, j] = 0
                    fpi_no3_vr_local[c, j] = 1.0 - fpi_nh4_vr_local[c, j]
                    actual_immob_no3_vr[c, j] = potential_immob_vr[c, j] - actual_immob_nh4_vr[c, j]
                    f_denit_vr[c, j] = pot_f_denit_vr[c, j]

                    if !local_use_fun
                        smin_no3_to_plant_vr[c, j] = plant_ndemand[c] * nuptake_prof[c, j] -
                            smin_nh4_to_plant_vr[c, j]
                    else
                        smin_no3_to_plant_vr[c, j] = smin_no3_vr[c, j] / dt -
                            actual_immob_no3_vr[c, j] - f_denit_vr[c, j]
                    end
                else
                    # NO3 limited: competition
                    nlimit_no3[c, j] = 1

                    if sum_no3_demand[c, j] > 0.0
                        if !local_use_fun
                            actual_immob_no3_vr[c, j] = min(
                                (smin_no3_vr[c, j] / dt) * ((potential_immob_vr[c, j] -
                                    actual_immob_nh4_vr[c, j]) * compet_decomp_no3 /
                                    sum_no3_demand_scaled[c, j]),
                                potential_immob_vr[c, j] - actual_immob_nh4_vr[c, j])

                            smin_no3_to_plant_vr[c, j] = min(
                                (smin_no3_vr[c, j] / dt) * ((plant_ndemand[c] *
                                    nuptake_prof[c, j] - smin_nh4_to_plant_vr[c, j]) * compet_plant_no3 /
                                    sum_no3_demand_scaled[c, j]),
                                plant_ndemand[c] * nuptake_prof[c, j] - smin_nh4_to_plant_vr[c, j])

                            f_denit_vr[c, j] = min(
                                (smin_no3_vr[c, j] / dt) * (pot_f_denit_vr[c, j] * compet_denit /
                                    sum_no3_demand_scaled[c, j]),
                                pot_f_denit_vr[c, j])
                        else
                            actual_immob_no3_vr[c, j] = min(
                                (smin_no3_vr[c, j] / dt) * ((potential_immob_vr[c, j] -
                                    actual_immob_nh4_vr[c, j]) * compet_decomp_no3 /
                                    sum_no3_demand_scaled[c, j]),
                                potential_immob_vr[c, j] - actual_immob_nh4_vr[c, j])

                            f_denit_vr[c, j] = min(
                                (smin_no3_vr[c, j] / dt) * (pot_f_denit_vr[c, j] * compet_denit /
                                    sum_no3_demand_scaled[c, j]),
                                pot_f_denit_vr[c, j])

                            smin_no3_to_plant_vr[c, j] = (smin_no3_vr[c, j] / dt) -
                                actual_immob_no3_vr[c, j] - f_denit_vr[c, j]
                        end
                    else
                        actual_immob_no3_vr[c, j] = 0.0
                        smin_no3_to_plant_vr[c, j] = 0.0
                        f_denit_vr[c, j] = 0.0
                    end

                    if potential_immob_vr[c, j] > 0.0
                        fpi_no3_vr_local[c, j] = actual_immob_no3_vr[c, j] / potential_immob_vr[c, j]
                    else
                        fpi_no3_vr_local[c, j] = 0.0
                    end
                end

                if mimics_decomp
                    fpi_no3_vr_local[c, j] = 1.0 - fpi_nh4_vr_local[c, j]
                    actual_immob_no3_vr[c, j] = potential_immob_vr[c, j] - actual_immob_nh4_vr[c, j]
                end

                # N2O emissions
                f_n2o_nit_vr[c, j] = f_nit_vr[c, j] * nitrif_n2o_loss_frac
                f_n2o_denit_vr[c, j] = f_denit_vr[c, j] / (1.0 + n2_n2o_ratio_denit_vr[c, j])

                # Carbon-only supplement
                if carbon_only
                    if fpi_no3_vr_local[c, j] + fpi_nh4_vr_local[c, j] < 1.0
                        fpi_nh4_vr_local[c, j] = 1.0 - fpi_no3_vr_local[c, j]
                        supplement_to_sminn_vr[c, j] = (potential_immob_vr[c, j] -
                            actual_immob_no3_vr[c, j]) - actual_immob_nh4_vr[c, j]
                        actual_immob_nh4_vr[c, j] = potential_immob_vr[c, j] - actual_immob_no3_vr[c, j]
                    end
                    if smin_no3_to_plant_vr[c, j] + smin_nh4_to_plant_vr[c, j] < plant_ndemand[c] * nuptake_prof[c, j]
                        supplement_to_sminn_vr[c, j] += (plant_ndemand[c] * nuptake_prof[c, j] -
                            smin_no3_to_plant_vr[c, j]) - smin_nh4_to_plant_vr[c, j]
                        smin_nh4_to_plant_vr[c, j] = plant_ndemand[c] * nuptake_prof[c, j] -
                            smin_no3_to_plant_vr[c, j]
                    end
                    sminn_to_plant_vr[c, j] = smin_no3_to_plant_vr[c, j] + smin_nh4_to_plant_vr[c, j]
                end

                # sum up no3 and nh4 fluxes
                fpi_vr[c, j] = fpi_no3_vr_local[c, j] + fpi_nh4_vr_local[c, j]
                sminn_to_plant_vr[c, j] = smin_no3_to_plant_vr[c, j] + smin_nh4_to_plant_vr[c, j]
                actual_immob_vr[c, j] = actual_immob_no3_vr[c, j] + actual_immob_nh4_vr[c, j]
            end
        end

        # --- Post-competition sums ---

        if !local_use_fun
            for c in bounds
                mask_bgc_soilc[c] || continue
                sminn_to_plant[c] = 0.0
            end
            for j in 1:nlevdecomp
                for c in bounds
                    mask_bgc_soilc[c] || continue
                    sminn_to_plant[c] += sminn_to_plant_vr[c, j] * dzsoi_decomp[j]
                end
            end
        else
            for c in bounds
                mask_bgc_soilc[c] || continue
                sminn_to_plant[c] = 0.0
            end
        end

        # --- MIMICS c_overflow_vr calculation ---
        if mimics_decomp
            for j in 1:nlevdecomp
                for c in bounds
                    mask_bgc_soilc[c] || continue
                    for k in 1:ndecomp_cascade_transitions
                        if cascade_receiver_pool[k] == i_cop_mic || cascade_receiver_pool[k] == i_oli_mic
                            sum_ndemand_vr[c, j] = sum_no3_demand_scaled[c, j] + sum_nh4_demand_scaled[c, j]
                            if pmnf_decomp_cascade[c, j, k] > 0.0 && sum_ndemand_vr[c, j] > 0.0
                                amnf_immob_vr = (sminn_vr[c, j] / dt) *
                                    (pmnf_decomp_cascade[c, j, k] / sum_ndemand_vr[c, j])
                                n_deficit_vr = pmnf_decomp_cascade[c, j, k] - amnf_immob_vr
                                c_overflow_vr[c, j, k] = n_deficit_vr *
                                    p_decomp_cn_gain[c, j, cascade_receiver_pool[k]]
                            else
                                c_overflow_vr[c, j, k] = 0.0
                            end
                        else
                            c_overflow_vr[c, j, k] = 0.0
                        end
                    end
                end
            end
        else
            c_overflow_vr .= 0.0
        end

        # --- Second pass residual distribution ---
        if !local_use_fun
            # Second pass for NH4
            for c in bounds
                mask_bgc_soilc[c] || continue
                residual_plant_ndemand[c] = plant_ndemand[c] - sminn_to_plant[c]
                residual_smin_nh4[c] = 0.0
            end
            for j in 1:nlevdecomp
                for c in bounds
                    mask_bgc_soilc[c] || continue
                    if residual_plant_ndemand[c] > 0.0
                        if nlimit_nh4[c, j] == 0
                            residual_smin_nh4_vr[c, j] = max(smin_nh4_vr[c, j] -
                                (actual_immob_nh4_vr[c, j] + smin_nh4_to_plant_vr[c, j] +
                                 f_nit_vr[c, j]) * dt, 0.0)
                            residual_smin_nh4[c] += residual_smin_nh4_vr[c, j] * dzsoi_decomp[j]
                        else
                            residual_smin_nh4_vr[c, j] = 0.0
                        end

                        if residual_smin_nh4[c] > 0.0 && nlimit_nh4[c, j] == 0
                            smin_nh4_to_plant_vr[c, j] += residual_smin_nh4_vr[c, j] *
                                min((residual_plant_ndemand[c] * dt) / residual_smin_nh4[c], 1.0) / dt
                        end
                    end
                end
            end

            # re-sum N fluxes after second pass for nh4
            for c in bounds
                mask_bgc_soilc[c] || continue
                sminn_to_plant[c] = 0.0
            end
            for j in 1:nlevdecomp
                for c in bounds
                    mask_bgc_soilc[c] || continue
                    sminn_to_plant_vr[c, j] = smin_nh4_to_plant_vr[c, j] + smin_no3_to_plant_vr[c, j]
                    sminn_to_plant[c] += sminn_to_plant_vr[c, j] * dzsoi_decomp[j]
                end
            end

            # Second pass for NO3
            for c in bounds
                mask_bgc_soilc[c] || continue
                residual_plant_ndemand[c] = plant_ndemand[c] - sminn_to_plant[c]
                residual_smin_no3[c] = 0.0
            end
            for j in 1:nlevdecomp
                for c in bounds
                    mask_bgc_soilc[c] || continue
                    if residual_plant_ndemand[c] > 0.0
                        if nlimit_no3[c, j] == 0
                            residual_smin_no3_vr[c, j] = max(smin_no3_vr[c, j] -
                                (actual_immob_no3_vr[c, j] + smin_no3_to_plant_vr[c, j] +
                                 f_denit_vr[c, j]) * dt, 0.0)
                            residual_smin_no3[c] += residual_smin_no3_vr[c, j] * dzsoi_decomp[j]
                        else
                            residual_smin_no3_vr[c, j] = 0.0
                        end

                        if residual_smin_no3[c] > 0.0 && nlimit_no3[c, j] == 0
                            smin_no3_to_plant_vr[c, j] += residual_smin_no3_vr[c, j] *
                                min((residual_plant_ndemand[c] * dt) / residual_smin_no3[c], 1.0) / dt
                        end
                    end
                end
            end

            # re-sum N fluxes after second passes
            for c in bounds
                mask_bgc_soilc[c] || continue
                sminn_to_plant[c] = 0.0
            end
            for j in 1:nlevdecomp
                for c in bounds
                    mask_bgc_soilc[c] || continue
                    sminn_to_plant_vr[c, j] = smin_nh4_to_plant_vr[c, j] + smin_no3_to_plant_vr[c, j]
                    sminn_to_plant[c] += sminn_to_plant_vr[c, j] * dzsoi_decomp[j]
                end
            end

        else
            # use_fun: calculate maximum N available to plants
            for c in bounds
                mask_bgc_soilc[c] || continue
                sminn_to_plant[c] = 0.0
            end
            for j in 1:nlevdecomp
                for c in bounds
                    mask_bgc_soilc[c] || continue
                    sminn_to_plant_vr[c, j] = smin_nh4_to_plant_vr[c, j] + smin_no3_to_plant_vr[c, j]
                    sminn_to_plant[c] += sminn_to_plant_vr[c, j] * dzsoi_decomp[j]
                end
            end

            # add up FUN fluxes from SMINN to plant
            for j in 1:nlevdecomp
                for c in bounds
                    mask_bgc_soilc[c] || continue
                    sminn_to_plant_new[c] += (sminn_to_plant_fun_no3_vr[c, j] +
                        sminn_to_plant_fun_nh4_vr[c, j]) * dzsoi_decomp[j]
                end
            end
        end

        # --- Sum up N fluxes to immobilization ---
        for c in bounds
            mask_bgc_soilc[c] || continue
            actual_immob[c] = 0.0
            potential_immob[c] = 0.0
        end
        for j in 1:nlevdecomp
            for c in bounds
                mask_bgc_soilc[c] || continue
                actual_immob[c] += actual_immob_vr[c, j] * dzsoi_decomp[j]
                potential_immob[c] += potential_immob_vr[c, j] * dzsoi_decomp[j]
            end
        end

        # --- Calculate fraction of potential growth / immobilization ---
        for c in bounds
            mask_bgc_soilc[c] || continue

            if !local_use_fun
                if plant_ndemand[c] > 0.0
                    fpg[c] = sminn_to_plant[c] / plant_ndemand[c]
                else
                    fpg[c] = 1.0
                end
            end

            if potential_immob[c] > 0.0
                fpi[c] = actual_immob[c] / potential_immob[c]
            else
                fpi[c] = 1.0
            end
        end

    end  # end if_nitrif

    return nothing
end
