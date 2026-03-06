# ==========================================================================
# Ported from: src/biogeochem/CNGapMortalityMod.F90
# Gap-phase mortality for coupled carbon-nitrogen code (CN).
#
# Public functions:
#   cn_gap_mortality!       — Patch-level gap mortality C & N fluxes
#   cn_gap_patch_to_column! — Aggregate patch mortality fluxes to column
# ==========================================================================

# ---------------------------------------------------------------------------
# Parameters for gap mortality (from CNGapMortalityMod params_type)
# ---------------------------------------------------------------------------

"""
    GapMortalityParams

Module-level parameters for the gap mortality routine.

Fields:
- `k_mort`: coefficient of growth efficiency in the mortality equation
- `r_mort`: per-PFT annual mortality rate (1/yr), indexed by vegetation type

Ported from `params_type` in `CNGapMortalityMod.F90`.
"""
Base.@kwdef mutable struct GapMortalityParams
    k_mort ::Float64        = 0.3     # coeff. of growth efficiency in mortality equation
    r_mort ::Vector{Float64} = Float64[]  # mortality rate (1/yr), one per PFT
end

# ---------------------------------------------------------------------------
# PFT constants needed by gap mortality (from pftcon)
# ---------------------------------------------------------------------------

"""
    PftConGapMort

PFT-level parameters referenced by the gap mortality and patch-to-column
routines.  Contains the subset of `pftconMod` fields used in
`CNGapMortalityMod.F90`.
"""
Base.@kwdef mutable struct PftConGapMort
    woody    ::Vector{Float64} = Float64[]  # binary woody flag (1=woody, 0=not woody)
    leafcn   ::Vector{Float64} = Float64[]  # leaf C:N (gC/gN)
    livewdcn ::Vector{Float64} = Float64[]  # live wood C:N (gC/gN)
    lf_f     ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # leaf litter fractions (pft, litr)
    fr_f     ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # fine root litter fractions (pft, litr)
end

# ---------------------------------------------------------------------------
# DGVS (Dynamic Global Vegetation) data needed by gap mortality
# ---------------------------------------------------------------------------

"""
    DgvsGapMortData

Dynamic-vegetation fields referenced by gap mortality when `use_cndv` is
enabled.

Ported from `dgvs_type` fields used in `CNGapMortalityMod.F90`.
"""
Base.@kwdef mutable struct DgvsGapMortData
    greffic_patch    ::Vector{Float64} = Float64[]  # growth efficiency
    heatstress_patch ::Vector{Float64} = Float64[]  # heat stress mortality
    nind_patch       ::Vector{Float64} = Float64[]  # number of individuals (#/m2)
end

# ---------------------------------------------------------------------------
# cn_gap_mortality! — Main gap mortality routine
# ---------------------------------------------------------------------------

"""
    cn_gap_mortality!(mask_soilp, bounds, params, pftcon, dgvs,
                      patch, canopystate, cnveg_cs, cnveg_cf, cnveg_ns, cnveg_nf;
                      dt, days_per_year, use_cndv, use_matrixcn, spinup_state,
                      npcropmin, spinup_factor_deadwood, secspday)

Gap-phase mortality routine for coupled carbon-nitrogen code.
Calculates patch-level mortality carbon and nitrogen fluxes.

Ported from `CNGapMortality` in `CNGapMortalityMod.F90`.
"""
function cn_gap_mortality!(mask_soilp::BitVector,
                           bounds::UnitRange{Int},
                           params::GapMortalityParams,
                           pftcon::PftConGapMort,
                           dgvs::DgvsGapMortData,
                           patch::PatchData,
                           canopystate::CanopyStateData,
                           cnveg_cs::CNVegCarbonStateData,
                           cnveg_cf::CNVegCarbonFluxData,
                           cnveg_ns::CNVegNitrogenStateData,
                           cnveg_nf::CNVegNitrogenFluxData;
                           dt::Float64 = 1800.0,
                           days_per_year::Float64 = 365.0,
                           use_cndv::Bool = false,
                           use_matrixcn::Bool = false,
                           spinup_state::Int = 0,
                           npcropmin::Int = 17,
                           spinup_factor_deadwood::Float64 = SPINUP_FACTOR_DEADWOOD_DEFAULT,
                           secspday::Float64 = SECSPDAY)

    # --- Aliases (matching Fortran associate block) ---
    ivt   = patch.itype
    woody = pftcon.woody

    greffic    = dgvs.greffic_patch
    heatstress = dgvs.heatstress_patch
    nind       = dgvs.nind_patch

    laisun = canopystate.laisun_patch
    laisha = canopystate.laisha_patch

    leafcn   = pftcon.leafcn
    livewdcn = pftcon.livewdcn

    k_mort = params.k_mort

    # --- Patch loop ---
    for p in bounds
        mask_soilp[p] || continue

        if use_cndv
            # Stress mortality from lpj's subr Mortality
            if woody[ivt[p] + 1] == 1.0

                if ivt[p] == 8
                    mort_max = 0.03  # BDT boreal
                else
                    mort_max = 0.01  # original value for all patches
                end

                # heatstress and greffic calculated in Establishment once/yr
                # Mortality rate inversely related to growth efficiency
                # (Prentice et al 1993)
                am = mort_max / (1.0 + k_mort * greffic[p])

                am = min(1.0, am + heatstress[p])
            else
                # lpj didn't set this for grasses; cn does
                am = params.r_mort[ivt[p] + 1]
            end
        else
            am = params.r_mort[ivt[p] + 1]
        end

        m = am / (days_per_year * secspday)

        # --------------------------------------------------
        # patch-level gap mortality carbon fluxes
        # --------------------------------------------------

        # displayed pools
        cnveg_cf.m_leafc_to_litter_patch[p]      = cnveg_cs.leafc_patch[p]      * m
        cnveg_cf.m_frootc_to_litter_patch[p]     = cnveg_cs.frootc_patch[p]     * m
        cnveg_cf.m_livestemc_to_litter_patch[p]  = cnveg_cs.livestemc_patch[p]  * m
        cnveg_cf.m_livecrootc_to_litter_patch[p] = cnveg_cs.livecrootc_patch[p] * m

        cnveg_cf.m_deadstemc_to_litter_patch[p]  = cnveg_cs.deadstemc_patch[p]  * m * spinup_factor_deadwood
        cnveg_cf.m_deadcrootc_to_litter_patch[p] = cnveg_cs.deadcrootc_patch[p] * m * spinup_factor_deadwood

        # storage pools
        cnveg_cf.m_leafc_storage_to_litter_patch[p]      = cnveg_cs.leafc_storage_patch[p]      * m
        cnveg_cf.m_frootc_storage_to_litter_patch[p]     = cnveg_cs.frootc_storage_patch[p]     * m
        cnveg_cf.m_livestemc_storage_to_litter_patch[p]  = cnveg_cs.livestemc_storage_patch[p]  * m
        cnveg_cf.m_deadstemc_storage_to_litter_patch[p]  = cnveg_cs.deadstemc_storage_patch[p]  * m
        cnveg_cf.m_livecrootc_storage_to_litter_patch[p] = cnveg_cs.livecrootc_storage_patch[p] * m
        cnveg_cf.m_deadcrootc_storage_to_litter_patch[p] = cnveg_cs.deadcrootc_storage_patch[p] * m
        cnveg_cf.m_gresp_storage_to_litter_patch[p]      = cnveg_cs.gresp_storage_patch[p]      * m

        # transfer pools
        cnveg_cf.m_leafc_xfer_to_litter_patch[p]      = cnveg_cs.leafc_xfer_patch[p]      * m
        cnveg_cf.m_frootc_xfer_to_litter_patch[p]     = cnveg_cs.frootc_xfer_patch[p]     * m
        cnveg_cf.m_livestemc_xfer_to_litter_patch[p]  = cnveg_cs.livestemc_xfer_patch[p]  * m
        cnveg_cf.m_deadstemc_xfer_to_litter_patch[p]  = cnveg_cs.deadstemc_xfer_patch[p]  * m
        cnveg_cf.m_livecrootc_xfer_to_litter_patch[p] = cnveg_cs.livecrootc_xfer_patch[p] * m
        cnveg_cf.m_deadcrootc_xfer_to_litter_patch[p] = cnveg_cs.deadcrootc_xfer_patch[p] * m
        cnveg_cf.m_gresp_xfer_to_litter_patch[p]      = cnveg_cs.gresp_xfer_patch[p]      * m

        # --------------------------------------------------
        # patch-level gap mortality nitrogen fluxes
        # --------------------------------------------------

        # displayed pools
        cnveg_nf.m_leafn_to_litter_patch[p]      = cnveg_ns.leafn_patch[p]      * m
        cnveg_nf.m_frootn_to_litter_patch[p]     = cnveg_ns.frootn_patch[p]     * m
        cnveg_nf.m_livestemn_to_litter_patch[p]  = cnveg_ns.livestemn_patch[p]  * m
        cnveg_nf.m_livecrootn_to_litter_patch[p] = cnveg_ns.livecrootn_patch[p] * m

        if spinup_state == 2 && !use_cndv
            # accelerate mortality of dead woody pools
            cnveg_nf.m_deadstemn_to_litter_patch[p]  = cnveg_ns.deadstemn_patch[p]  * m * spinup_factor_deadwood
            cnveg_nf.m_deadcrootn_to_litter_patch[p] = cnveg_ns.deadcrootn_patch[p] * m * spinup_factor_deadwood
        else
            cnveg_nf.m_deadstemn_to_litter_patch[p]  = cnveg_ns.deadstemn_patch[p]  * m
            cnveg_nf.m_deadcrootn_to_litter_patch[p] = cnveg_ns.deadcrootn_patch[p] * m
        end

        if ivt[p] < npcropmin
            cnveg_nf.m_retransn_to_litter_patch[p] = cnveg_ns.retransn_patch[p] * m
        end

        # storage pools
        cnveg_nf.m_leafn_storage_to_litter_patch[p]      = cnveg_ns.leafn_storage_patch[p]      * m
        cnveg_nf.m_frootn_storage_to_litter_patch[p]     = cnveg_ns.frootn_storage_patch[p]     * m
        cnveg_nf.m_livestemn_storage_to_litter_patch[p]  = cnveg_ns.livestemn_storage_patch[p]  * m
        cnveg_nf.m_deadstemn_storage_to_litter_patch[p]  = cnveg_ns.deadstemn_storage_patch[p]  * m
        cnveg_nf.m_livecrootn_storage_to_litter_patch[p] = cnveg_ns.livecrootn_storage_patch[p] * m
        cnveg_nf.m_deadcrootn_storage_to_litter_patch[p] = cnveg_ns.deadcrootn_storage_patch[p] * m

        # transfer pools
        cnveg_nf.m_leafn_xfer_to_litter_patch[p]      = cnveg_ns.leafn_xfer_patch[p]      * m
        cnveg_nf.m_frootn_xfer_to_litter_patch[p]     = cnveg_ns.frootn_xfer_patch[p]     * m
        cnveg_nf.m_livestemn_xfer_to_litter_patch[p]  = cnveg_ns.livestemn_xfer_patch[p]  * m
        cnveg_nf.m_deadstemn_xfer_to_litter_patch[p]  = cnveg_ns.deadstemn_xfer_patch[p]  * m
        cnveg_nf.m_livecrootn_xfer_to_litter_patch[p] = cnveg_ns.livecrootn_xfer_patch[p] * m
        cnveg_nf.m_deadcrootn_xfer_to_litter_patch[p] = cnveg_ns.deadcrootn_xfer_patch[p] * m

        # added by F. Li and S. Levis
        if use_cndv
            if woody[ivt[p] + 1] == 1.0
                if cnveg_cs.livestemc_patch[p] + cnveg_cs.deadstemc_patch[p] > 0.0
                    nind[p] = nind[p] * (1.0 - m)
                else
                    nind[p] = 0.0
                end
            end
        end

    end  # end of patch loop

    return nothing
end

# ---------------------------------------------------------------------------
# cn_gap_patch_to_column! — Gather patch mortality fluxes to column level
# ---------------------------------------------------------------------------

"""
    cn_gap_patch_to_column!(mask_soilp, bounds, pftcon, patch,
                            cnveg_cf, cnveg_nf,
                            leaf_prof, froot_prof, croot_prof, stem_prof;
                            nlevdecomp, i_litr_min, i_litr_max, i_met_lit)

Gathers all patch-level gap mortality fluxes to the column level and assigns
them to the appropriate litter and coarse woody debris (CWD) pools.

Ported from `CNGap_PatchToColumn` in `CNGapMortalityMod.F90`.
"""
function cn_gap_patch_to_column!(mask_soilp::BitVector,
                                 bounds::UnitRange{Int},
                                 pftcon::PftConGapMort,
                                 patch::PatchData,
                                 cnveg_cf::CNVegCarbonFluxData,
                                 cnveg_nf::CNVegNitrogenFluxData,
                                 leaf_prof::Matrix{Float64},
                                 froot_prof::Matrix{Float64},
                                 croot_prof::Matrix{Float64},
                                 stem_prof::Matrix{Float64};
                                 nlevdecomp::Int = 1,
                                 i_litr_min::Int = 1,
                                 i_litr_max::Int = 3,
                                 i_met_lit::Int = 1)

    # --- Aliases ---
    ivt   = patch.itype
    wtcol = patch.wtcol
    lf_f  = pftcon.lf_f
    fr_f  = pftcon.fr_f

    # Carbon flux aliases
    m_leafc_to_litter               = cnveg_cf.m_leafc_to_litter_patch
    m_frootc_to_litter              = cnveg_cf.m_frootc_to_litter_patch
    m_livestemc_to_litter           = cnveg_cf.m_livestemc_to_litter_patch
    m_deadstemc_to_litter           = cnveg_cf.m_deadstemc_to_litter_patch
    m_livecrootc_to_litter          = cnveg_cf.m_livecrootc_to_litter_patch
    m_deadcrootc_to_litter          = cnveg_cf.m_deadcrootc_to_litter_patch
    m_leafc_storage_to_litter       = cnveg_cf.m_leafc_storage_to_litter_patch
    m_frootc_storage_to_litter      = cnveg_cf.m_frootc_storage_to_litter_patch
    m_livestemc_storage_to_litter   = cnveg_cf.m_livestemc_storage_to_litter_patch
    m_deadstemc_storage_to_litter   = cnveg_cf.m_deadstemc_storage_to_litter_patch
    m_livecrootc_storage_to_litter  = cnveg_cf.m_livecrootc_storage_to_litter_patch
    m_deadcrootc_storage_to_litter  = cnveg_cf.m_deadcrootc_storage_to_litter_patch
    m_gresp_storage_to_litter       = cnveg_cf.m_gresp_storage_to_litter_patch
    m_leafc_xfer_to_litter          = cnveg_cf.m_leafc_xfer_to_litter_patch
    m_frootc_xfer_to_litter         = cnveg_cf.m_frootc_xfer_to_litter_patch
    m_livestemc_xfer_to_litter      = cnveg_cf.m_livestemc_xfer_to_litter_patch
    m_deadstemc_xfer_to_litter      = cnveg_cf.m_deadstemc_xfer_to_litter_patch
    m_livecrootc_xfer_to_litter     = cnveg_cf.m_livecrootc_xfer_to_litter_patch
    m_deadcrootc_xfer_to_litter     = cnveg_cf.m_deadcrootc_xfer_to_litter_patch
    m_gresp_xfer_to_litter          = cnveg_cf.m_gresp_xfer_to_litter_patch

    gap_mortality_c_to_litr_c = cnveg_cf.gap_mortality_c_to_litr_c_col
    gap_mortality_c_to_cwdc   = cnveg_cf.gap_mortality_c_to_cwdc_col

    # Nitrogen flux aliases
    m_leafn_to_litter               = cnveg_nf.m_leafn_to_litter_patch
    m_frootn_to_litter              = cnveg_nf.m_frootn_to_litter_patch
    m_livestemn_to_litter           = cnveg_nf.m_livestemn_to_litter_patch
    m_deadstemn_to_litter           = cnveg_nf.m_deadstemn_to_litter_patch
    m_livecrootn_to_litter          = cnveg_nf.m_livecrootn_to_litter_patch
    m_deadcrootn_to_litter          = cnveg_nf.m_deadcrootn_to_litter_patch
    m_retransn_to_litter            = cnveg_nf.m_retransn_to_litter_patch
    m_leafn_storage_to_litter       = cnveg_nf.m_leafn_storage_to_litter_patch
    m_frootn_storage_to_litter      = cnveg_nf.m_frootn_storage_to_litter_patch
    m_livestemn_storage_to_litter   = cnveg_nf.m_livestemn_storage_to_litter_patch
    m_deadstemn_storage_to_litter   = cnveg_nf.m_deadstemn_storage_to_litter_patch
    m_livecrootn_storage_to_litter  = cnveg_nf.m_livecrootn_storage_to_litter_patch
    m_deadcrootn_storage_to_litter  = cnveg_nf.m_deadcrootn_storage_to_litter_patch
    m_leafn_xfer_to_litter          = cnveg_nf.m_leafn_xfer_to_litter_patch
    m_frootn_xfer_to_litter         = cnveg_nf.m_frootn_xfer_to_litter_patch
    m_livestemn_xfer_to_litter      = cnveg_nf.m_livestemn_xfer_to_litter_patch
    m_deadstemn_xfer_to_litter      = cnveg_nf.m_deadstemn_xfer_to_litter_patch
    m_livecrootn_xfer_to_litter     = cnveg_nf.m_livecrootn_xfer_to_litter_patch
    m_deadcrootn_xfer_to_litter     = cnveg_nf.m_deadcrootn_xfer_to_litter_patch

    gap_mortality_n_to_litr_n = cnveg_nf.gap_mortality_n_to_litr_n_col
    gap_mortality_n_to_cwdn   = cnveg_nf.gap_mortality_n_to_cwdn_col

    # --- Decomp level / patch loop ---
    for j in 1:nlevdecomp
        for p in bounds
            mask_soilp[p] || continue

            c = patch.column[p]

            # -- Leaf and fine root gap mortality C to litter pools --
            for i in i_litr_min:i_litr_max
                gap_mortality_c_to_litr_c[c, j, i] += (
                    m_leafc_to_litter[p]  * lf_f[ivt[p] + 1, i] * wtcol[p] * leaf_prof[p, j] +
                    m_frootc_to_litter[p] * fr_f[ivt[p] + 1, i] * wtcol[p] * froot_prof[p, j]
                )
            end

            # -- Wood gap mortality C to CWD --
            gap_mortality_c_to_cwdc[c, j] += (
                (m_livestemc_to_litter[p] + m_deadstemc_to_litter[p]) * wtcol[p] * stem_prof[p, j]
            )
            gap_mortality_c_to_cwdc[c, j] += (
                (m_livecrootc_to_litter[p] + m_deadcrootc_to_litter[p]) * wtcol[p] * croot_prof[p, j]
            )

            # -- Storage and transfer gap mortality C to metabolic litter --
            gap_mortality_c_to_litr_c[c, j, i_met_lit] += (
                (m_leafc_storage_to_litter[p] + m_gresp_storage_to_litter[p]) * wtcol[p] * leaf_prof[p, j] +
                m_frootc_storage_to_litter[p] * wtcol[p] * froot_prof[p, j] +
                (m_livestemc_storage_to_litter[p] + m_deadstemc_storage_to_litter[p]) * wtcol[p] * stem_prof[p, j] +
                (m_livecrootc_storage_to_litter[p] + m_deadcrootc_storage_to_litter[p]) * wtcol[p] * croot_prof[p, j] +
                # transfer gap mortality carbon fluxes
                (m_leafc_xfer_to_litter[p] + m_gresp_xfer_to_litter[p]) * wtcol[p] * leaf_prof[p, j] +
                m_frootc_xfer_to_litter[p] * wtcol[p] * froot_prof[p, j] +
                (m_livestemc_xfer_to_litter[p] + m_deadstemc_xfer_to_litter[p]) * wtcol[p] * stem_prof[p, j] +
                (m_livecrootc_xfer_to_litter[p] + m_deadcrootc_xfer_to_litter[p]) * wtcol[p] * croot_prof[p, j]
            )

            # -- Leaf and fine root gap mortality N to litter pools --
            for i in i_litr_min:i_litr_max
                gap_mortality_n_to_litr_n[c, j, i] += (
                    m_leafn_to_litter[p]  * lf_f[ivt[p] + 1, i] * wtcol[p] * leaf_prof[p, j] +
                    m_frootn_to_litter[p] * fr_f[ivt[p] + 1, i] * wtcol[p] * froot_prof[p, j]
                )
            end

            # -- Wood gap mortality N to CWD --
            gap_mortality_n_to_cwdn[c, j] += (
                (m_livestemn_to_litter[p] + m_deadstemn_to_litter[p]) * wtcol[p] * stem_prof[p, j]
            )
            gap_mortality_n_to_cwdn[c, j] += (
                (m_livecrootn_to_litter[p] + m_deadcrootn_to_litter[p]) * wtcol[p] * croot_prof[p, j]
            )

            # -- Storage, transfer, and retranslocated N to metabolic litter --
            gap_mortality_n_to_litr_n[c, j, i_met_lit] += (
                # retranslocated N pool gap mortality fluxes
                m_retransn_to_litter[p] * wtcol[p] * leaf_prof[p, j] +
                # storage gap mortality nitrogen fluxes
                m_leafn_storage_to_litter[p] * wtcol[p] * leaf_prof[p, j] +
                m_frootn_storage_to_litter[p] * wtcol[p] * froot_prof[p, j] +
                (m_livestemn_storage_to_litter[p] + m_deadstemn_storage_to_litter[p]) * wtcol[p] * stem_prof[p, j] +
                (m_livecrootn_storage_to_litter[p] + m_deadcrootn_storage_to_litter[p]) * wtcol[p] * croot_prof[p, j] +
                # transfer gap mortality nitrogen fluxes
                m_leafn_xfer_to_litter[p] * wtcol[p] * leaf_prof[p, j] +
                m_frootn_xfer_to_litter[p] * wtcol[p] * froot_prof[p, j] +
                (m_livestemn_xfer_to_litter[p] + m_deadstemn_xfer_to_litter[p]) * wtcol[p] * stem_prof[p, j] +
                (m_livecrootn_xfer_to_litter[p] + m_deadcrootn_xfer_to_litter[p]) * wtcol[p] * croot_prof[p, j]
            )

        end  # patch loop
    end  # decomp level loop

    return nothing
end
