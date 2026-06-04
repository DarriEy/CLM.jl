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
# Per-patch gap-mortality kernel: each thread computes the ~20 mortality C & N
# fluxes for one patch p (all independent per-patch outputs) plus the optional
# use_cndv `nind` individual-count update. Config flags (use_cndv, spinup_state,
# npcropmin) are resolved on the host and passed as Bool/Int scalars; no branch
# on a String and no error()/@warn in the kernel.
@kernel function _gapmort_patch_kernel!(
        m_leafc_to_litter, m_frootc_to_litter, m_livestemc_to_litter,
        m_livecrootc_to_litter, m_deadstemc_to_litter, m_deadcrootc_to_litter,
        m_leafc_storage_to_litter, m_frootc_storage_to_litter,
        m_livestemc_storage_to_litter, m_deadstemc_storage_to_litter,
        m_livecrootc_storage_to_litter, m_deadcrootc_storage_to_litter,
        m_gresp_storage_to_litter,
        m_leafc_xfer_to_litter, m_frootc_xfer_to_litter,
        m_livestemc_xfer_to_litter, m_deadstemc_xfer_to_litter,
        m_livecrootc_xfer_to_litter, m_deadcrootc_xfer_to_litter,
        m_gresp_xfer_to_litter,
        m_leafn_to_litter, m_frootn_to_litter, m_livestemn_to_litter,
        m_livecrootn_to_litter, m_deadstemn_to_litter, m_deadcrootn_to_litter,
        m_retransn_to_litter,
        m_leafn_storage_to_litter, m_frootn_storage_to_litter,
        m_livestemn_storage_to_litter, m_deadstemn_storage_to_litter,
        m_livecrootn_storage_to_litter, m_deadcrootn_storage_to_litter,
        m_leafn_xfer_to_litter, m_frootn_xfer_to_litter,
        m_livestemn_xfer_to_litter, m_deadstemn_xfer_to_litter,
        m_livecrootn_xfer_to_litter, m_deadcrootn_xfer_to_litter,
        nind,
        @Const(mask_soilp), @Const(ivt), @Const(woody),
        @Const(greffic), @Const(heatstress),
        @Const(r_mort),
        @Const(leafc), @Const(frootc), @Const(livestemc), @Const(livecrootc),
        @Const(deadstemc), @Const(deadcrootc),
        @Const(leafc_storage), @Const(frootc_storage), @Const(livestemc_storage),
        @Const(deadstemc_storage), @Const(livecrootc_storage),
        @Const(deadcrootc_storage), @Const(gresp_storage),
        @Const(leafc_xfer), @Const(frootc_xfer), @Const(livestemc_xfer),
        @Const(deadstemc_xfer), @Const(livecrootc_xfer), @Const(deadcrootc_xfer),
        @Const(gresp_xfer),
        @Const(leafn), @Const(frootn), @Const(livestemn), @Const(livecrootn),
        @Const(deadstemn), @Const(deadcrootn), @Const(retransn),
        @Const(leafn_storage), @Const(frootn_storage), @Const(livestemn_storage),
        @Const(deadstemn_storage), @Const(livecrootn_storage),
        @Const(deadcrootn_storage),
        @Const(leafn_xfer), @Const(frootn_xfer), @Const(livestemn_xfer),
        @Const(deadstemn_xfer), @Const(livecrootn_xfer), @Const(deadcrootn_xfer),
        k_mort, days_per_year, secspday, spinup_factor_deadwood,
        use_cndv::Bool, spinup_state::Int, npcropmin::Int)

    p = @index(Global)
    @inbounds if mask_soilp[p]
        ivtp = ivt[p] + 1   # 0-based Fortran → 1-based Julia

        if use_cndv
            # Stress mortality from lpj's subr Mortality
            if woody[ivtp] == 1.0
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
                am = r_mort[ivtp]
            end
        else
            am = r_mort[ivtp]
        end

        m = am / (days_per_year * secspday)

        # --------------------------------------------------
        # patch-level gap mortality carbon fluxes
        # --------------------------------------------------

        # displayed pools
        m_leafc_to_litter[p]      = leafc[p]      * m
        m_frootc_to_litter[p]     = frootc[p]     * m
        m_livestemc_to_litter[p]  = livestemc[p]  * m
        m_livecrootc_to_litter[p] = livecrootc[p] * m

        m_deadstemc_to_litter[p]  = deadstemc[p]  * m * spinup_factor_deadwood
        m_deadcrootc_to_litter[p] = deadcrootc[p] * m * spinup_factor_deadwood

        # storage pools
        m_leafc_storage_to_litter[p]      = leafc_storage[p]      * m
        m_frootc_storage_to_litter[p]     = frootc_storage[p]     * m
        m_livestemc_storage_to_litter[p]  = livestemc_storage[p]  * m
        m_deadstemc_storage_to_litter[p]  = deadstemc_storage[p]  * m
        m_livecrootc_storage_to_litter[p] = livecrootc_storage[p] * m
        m_deadcrootc_storage_to_litter[p] = deadcrootc_storage[p] * m
        m_gresp_storage_to_litter[p]      = gresp_storage[p]      * m

        # transfer pools
        m_leafc_xfer_to_litter[p]      = leafc_xfer[p]      * m
        m_frootc_xfer_to_litter[p]     = frootc_xfer[p]     * m
        m_livestemc_xfer_to_litter[p]  = livestemc_xfer[p]  * m
        m_deadstemc_xfer_to_litter[p]  = deadstemc_xfer[p]  * m
        m_livecrootc_xfer_to_litter[p] = livecrootc_xfer[p] * m
        m_deadcrootc_xfer_to_litter[p] = deadcrootc_xfer[p] * m
        m_gresp_xfer_to_litter[p]      = gresp_xfer[p]      * m

        # --------------------------------------------------
        # patch-level gap mortality nitrogen fluxes
        # --------------------------------------------------

        # displayed pools
        m_leafn_to_litter[p]      = leafn[p]      * m
        m_frootn_to_litter[p]     = frootn[p]     * m
        m_livestemn_to_litter[p]  = livestemn[p]  * m
        m_livecrootn_to_litter[p] = livecrootn[p] * m

        if spinup_state == 2 && !use_cndv
            # accelerate mortality of dead woody pools
            m_deadstemn_to_litter[p]  = deadstemn[p]  * m * spinup_factor_deadwood
            m_deadcrootn_to_litter[p] = deadcrootn[p] * m * spinup_factor_deadwood
        else
            m_deadstemn_to_litter[p]  = deadstemn[p]  * m
            m_deadcrootn_to_litter[p] = deadcrootn[p] * m
        end

        if ivt[p] < npcropmin
            m_retransn_to_litter[p] = retransn[p] * m
        end

        # storage pools
        m_leafn_storage_to_litter[p]      = leafn_storage[p]      * m
        m_frootn_storage_to_litter[p]     = frootn_storage[p]     * m
        m_livestemn_storage_to_litter[p]  = livestemn_storage[p]  * m
        m_deadstemn_storage_to_litter[p]  = deadstemn_storage[p]  * m
        m_livecrootn_storage_to_litter[p] = livecrootn_storage[p] * m
        m_deadcrootn_storage_to_litter[p] = deadcrootn_storage[p] * m

        # transfer pools
        m_leafn_xfer_to_litter[p]      = leafn_xfer[p]      * m
        m_frootn_xfer_to_litter[p]     = frootn_xfer[p]     * m
        m_livestemn_xfer_to_litter[p]  = livestemn_xfer[p]  * m
        m_deadstemn_xfer_to_litter[p]  = deadstemn_xfer[p]  * m
        m_livecrootn_xfer_to_litter[p] = livecrootn_xfer[p] * m
        m_deadcrootn_xfer_to_litter[p] = deadcrootn_xfer[p] * m

        # added by F. Li and S. Levis
        if use_cndv
            if woody[ivtp] == 1.0
                if livestemc[p] + deadstemc[p] > 0.0
                    nind[p] = nind[p] * (1.0 - m)
                else
                    nind[p] = 0.0
                end
            end
        end
    end
end

function cn_gap_mortality!(mask_soilp::AbstractVector{Bool},
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
                           dt::Real = 1800.0,
                           days_per_year::Real = 365.0,
                           use_cndv::Bool = false,
                           use_matrixcn::Bool = false,
                           spinup_state::Int = 0,
                           npcropmin::Int = 17,
                           spinup_factor_deadwood::Real = SPINUP_FACTOR_DEADWOOD_DEFAULT,
                           secspday::Real = SECSPDAY)

    # The Fortran loop runs over a filtered patch range `bounds`; the mask is
    # active only on that range, so launching one thread per patch over the full
    # mask reproduces it (the kernel skips inactive patches). bounds == 1:np here.
    _launch!(_gapmort_patch_kernel!,
        cnveg_cf.m_leafc_to_litter_patch,
        cnveg_cf.m_frootc_to_litter_patch,
        cnveg_cf.m_livestemc_to_litter_patch,
        cnveg_cf.m_livecrootc_to_litter_patch,
        cnveg_cf.m_deadstemc_to_litter_patch,
        cnveg_cf.m_deadcrootc_to_litter_patch,
        cnveg_cf.m_leafc_storage_to_litter_patch,
        cnveg_cf.m_frootc_storage_to_litter_patch,
        cnveg_cf.m_livestemc_storage_to_litter_patch,
        cnveg_cf.m_deadstemc_storage_to_litter_patch,
        cnveg_cf.m_livecrootc_storage_to_litter_patch,
        cnveg_cf.m_deadcrootc_storage_to_litter_patch,
        cnveg_cf.m_gresp_storage_to_litter_patch,
        cnveg_cf.m_leafc_xfer_to_litter_patch,
        cnveg_cf.m_frootc_xfer_to_litter_patch,
        cnveg_cf.m_livestemc_xfer_to_litter_patch,
        cnveg_cf.m_deadstemc_xfer_to_litter_patch,
        cnveg_cf.m_livecrootc_xfer_to_litter_patch,
        cnveg_cf.m_deadcrootc_xfer_to_litter_patch,
        cnveg_cf.m_gresp_xfer_to_litter_patch,
        cnveg_nf.m_leafn_to_litter_patch,
        cnveg_nf.m_frootn_to_litter_patch,
        cnveg_nf.m_livestemn_to_litter_patch,
        cnveg_nf.m_livecrootn_to_litter_patch,
        cnveg_nf.m_deadstemn_to_litter_patch,
        cnveg_nf.m_deadcrootn_to_litter_patch,
        cnveg_nf.m_retransn_to_litter_patch,
        cnveg_nf.m_leafn_storage_to_litter_patch,
        cnveg_nf.m_frootn_storage_to_litter_patch,
        cnveg_nf.m_livestemn_storage_to_litter_patch,
        cnveg_nf.m_deadstemn_storage_to_litter_patch,
        cnveg_nf.m_livecrootn_storage_to_litter_patch,
        cnveg_nf.m_deadcrootn_storage_to_litter_patch,
        cnveg_nf.m_leafn_xfer_to_litter_patch,
        cnveg_nf.m_frootn_xfer_to_litter_patch,
        cnveg_nf.m_livestemn_xfer_to_litter_patch,
        cnveg_nf.m_deadstemn_xfer_to_litter_patch,
        cnveg_nf.m_livecrootn_xfer_to_litter_patch,
        cnveg_nf.m_deadcrootn_xfer_to_litter_patch,
        dgvs.nind_patch,
        mask_soilp, patch.itype, pftcon.woody,
        dgvs.greffic_patch, dgvs.heatstress_patch,
        params.r_mort,
        cnveg_cs.leafc_patch, cnveg_cs.frootc_patch, cnveg_cs.livestemc_patch,
        cnveg_cs.livecrootc_patch, cnveg_cs.deadstemc_patch, cnveg_cs.deadcrootc_patch,
        cnveg_cs.leafc_storage_patch, cnveg_cs.frootc_storage_patch,
        cnveg_cs.livestemc_storage_patch, cnveg_cs.deadstemc_storage_patch,
        cnveg_cs.livecrootc_storage_patch, cnveg_cs.deadcrootc_storage_patch,
        cnveg_cs.gresp_storage_patch,
        cnveg_cs.leafc_xfer_patch, cnveg_cs.frootc_xfer_patch,
        cnveg_cs.livestemc_xfer_patch, cnveg_cs.deadstemc_xfer_patch,
        cnveg_cs.livecrootc_xfer_patch, cnveg_cs.deadcrootc_xfer_patch,
        cnveg_cs.gresp_xfer_patch,
        cnveg_ns.leafn_patch, cnveg_ns.frootn_patch, cnveg_ns.livestemn_patch,
        cnveg_ns.livecrootn_patch, cnveg_ns.deadstemn_patch, cnveg_ns.deadcrootn_patch,
        cnveg_ns.retransn_patch,
        cnveg_ns.leafn_storage_patch, cnveg_ns.frootn_storage_patch,
        cnveg_ns.livestemn_storage_patch, cnveg_ns.deadstemn_storage_patch,
        cnveg_ns.livecrootn_storage_patch, cnveg_ns.deadcrootn_storage_patch,
        cnveg_ns.leafn_xfer_patch, cnveg_ns.frootn_xfer_patch,
        cnveg_ns.livestemn_xfer_patch, cnveg_ns.deadstemn_xfer_patch,
        cnveg_ns.livecrootn_xfer_patch, cnveg_ns.deadcrootn_xfer_patch,
        params.k_mort, days_per_year, secspday, spinup_factor_deadwood,
        use_cndv, spinup_state, npcropmin)

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
# Per-patch patch→column scatter kernel: each thread handles one patch p and
# scatters its gap-mortality C & N fluxes into the column litter / CWD pools over
# all decomp levels j and litter pools i. The host loop was j-outer/p-inner; for a
# fixed (c,j,i) accumulator the contributions are summed in ascending-p order in
# both the host loop and the KA-CPU launch (threads iterate p ascending), so the
# floating-point accumulation order — and the result — is byte-identical on CPU.
# The two C-to-CWD (stem then croot) and the two N-to-CWD adds are kept as two
# separate _scatter_add! calls each, matching the host's two separate `+=`.
@kernel function _gap_p2c_kernel!(
        gap_mortality_c_to_litr_c, gap_mortality_c_to_cwdc,
        gap_mortality_n_to_litr_n, gap_mortality_n_to_cwdn,
        @Const(mask_soilp), @Const(ivt), @Const(column), @Const(wtcol),
        @Const(lf_f), @Const(fr_f),
        @Const(leaf_prof), @Const(froot_prof), @Const(croot_prof), @Const(stem_prof),
        @Const(m_leafc_to_litter), @Const(m_frootc_to_litter),
        @Const(m_livestemc_to_litter), @Const(m_deadstemc_to_litter),
        @Const(m_livecrootc_to_litter), @Const(m_deadcrootc_to_litter),
        @Const(m_leafc_storage_to_litter), @Const(m_frootc_storage_to_litter),
        @Const(m_livestemc_storage_to_litter), @Const(m_deadstemc_storage_to_litter),
        @Const(m_livecrootc_storage_to_litter), @Const(m_deadcrootc_storage_to_litter),
        @Const(m_gresp_storage_to_litter),
        @Const(m_leafc_xfer_to_litter), @Const(m_frootc_xfer_to_litter),
        @Const(m_livestemc_xfer_to_litter), @Const(m_deadstemc_xfer_to_litter),
        @Const(m_livecrootc_xfer_to_litter), @Const(m_deadcrootc_xfer_to_litter),
        @Const(m_gresp_xfer_to_litter),
        @Const(m_leafn_to_litter), @Const(m_frootn_to_litter),
        @Const(m_livestemn_to_litter), @Const(m_deadstemn_to_litter),
        @Const(m_livecrootn_to_litter), @Const(m_deadcrootn_to_litter),
        @Const(m_retransn_to_litter),
        @Const(m_leafn_storage_to_litter), @Const(m_frootn_storage_to_litter),
        @Const(m_livestemn_storage_to_litter), @Const(m_deadstemn_storage_to_litter),
        @Const(m_livecrootn_storage_to_litter), @Const(m_deadcrootn_storage_to_litter),
        @Const(m_leafn_xfer_to_litter), @Const(m_frootn_xfer_to_litter),
        @Const(m_livestemn_xfer_to_litter), @Const(m_deadstemn_xfer_to_litter),
        @Const(m_livecrootn_xfer_to_litter), @Const(m_deadcrootn_xfer_to_litter),
        nlevdecomp::Int, i_litr_min::Int, i_litr_max::Int, i_met_lit::Int)

    p = @index(Global)
    @inbounds if mask_soilp[p]
        ivtp = ivt[p] + 1
        c    = column[p]
        wt   = wtcol[p]

        for j in 1:nlevdecomp
            lp = leaf_prof[p, j]
            fp = froot_prof[p, j]
            cp = croot_prof[p, j]
            sp = stem_prof[p, j]

            # -- Leaf and fine root gap mortality C to litter pools --
            for i in i_litr_min:i_litr_max
                _scatter_add!(gap_mortality_c_to_litr_c, c, j, i,
                    m_leafc_to_litter[p]  * lf_f[ivtp, i] * wt * lp +
                    m_frootc_to_litter[p] * fr_f[ivtp, i] * wt * fp)
            end

            # -- Wood gap mortality C to CWD --
            _scatter_add!(gap_mortality_c_to_cwdc, c, j,
                (m_livestemc_to_litter[p] + m_deadstemc_to_litter[p]) * wt * sp)
            _scatter_add!(gap_mortality_c_to_cwdc, c, j,
                (m_livecrootc_to_litter[p] + m_deadcrootc_to_litter[p]) * wt * cp)

            # -- Storage and transfer gap mortality C to metabolic litter --
            _scatter_add!(gap_mortality_c_to_litr_c, c, j, i_met_lit,
                (m_leafc_storage_to_litter[p] + m_gresp_storage_to_litter[p]) * wt * lp +
                m_frootc_storage_to_litter[p] * wt * fp +
                (m_livestemc_storage_to_litter[p] + m_deadstemc_storage_to_litter[p]) * wt * sp +
                (m_livecrootc_storage_to_litter[p] + m_deadcrootc_storage_to_litter[p]) * wt * cp +
                # transfer gap mortality carbon fluxes
                (m_leafc_xfer_to_litter[p] + m_gresp_xfer_to_litter[p]) * wt * lp +
                m_frootc_xfer_to_litter[p] * wt * fp +
                (m_livestemc_xfer_to_litter[p] + m_deadstemc_xfer_to_litter[p]) * wt * sp +
                (m_livecrootc_xfer_to_litter[p] + m_deadcrootc_xfer_to_litter[p]) * wt * cp)

            # -- Leaf and fine root gap mortality N to litter pools --
            for i in i_litr_min:i_litr_max
                _scatter_add!(gap_mortality_n_to_litr_n, c, j, i,
                    m_leafn_to_litter[p]  * lf_f[ivtp, i] * wt * lp +
                    m_frootn_to_litter[p] * fr_f[ivtp, i] * wt * fp)
            end

            # -- Wood gap mortality N to CWD --
            _scatter_add!(gap_mortality_n_to_cwdn, c, j,
                (m_livestemn_to_litter[p] + m_deadstemn_to_litter[p]) * wt * sp)
            _scatter_add!(gap_mortality_n_to_cwdn, c, j,
                (m_livecrootn_to_litter[p] + m_deadcrootn_to_litter[p]) * wt * cp)

            # -- Storage, transfer, and retranslocated N to metabolic litter --
            _scatter_add!(gap_mortality_n_to_litr_n, c, j, i_met_lit,
                # retranslocated N pool gap mortality fluxes
                m_retransn_to_litter[p] * wt * lp +
                # storage gap mortality nitrogen fluxes
                m_leafn_storage_to_litter[p] * wt * lp +
                m_frootn_storage_to_litter[p] * wt * fp +
                (m_livestemn_storage_to_litter[p] + m_deadstemn_storage_to_litter[p]) * wt * sp +
                (m_livecrootn_storage_to_litter[p] + m_deadcrootn_storage_to_litter[p]) * wt * cp +
                # transfer gap mortality nitrogen fluxes
                m_leafn_xfer_to_litter[p] * wt * lp +
                m_frootn_xfer_to_litter[p] * wt * fp +
                (m_livestemn_xfer_to_litter[p] + m_deadstemn_xfer_to_litter[p]) * wt * sp +
                (m_livecrootn_xfer_to_litter[p] + m_deadcrootn_xfer_to_litter[p]) * wt * cp)
        end
    end
end

function cn_gap_patch_to_column!(mask_soilp::AbstractVector{Bool},
                                 bounds::UnitRange{Int},
                                 pftcon::PftConGapMort,
                                 patch::PatchData,
                                 cnveg_cf::CNVegCarbonFluxData,
                                 cnveg_nf::CNVegNitrogenFluxData,
                                 leaf_prof::AbstractMatrix{<:Real},
                                 froot_prof::AbstractMatrix{<:Real},
                                 croot_prof::AbstractMatrix{<:Real},
                                 stem_prof::AbstractMatrix{<:Real};
                                 nlevdecomp::Int = 1,
                                 i_litr_min::Int = 1,
                                 i_litr_max::Int = 3,
                                 i_met_lit::Int = 1)

    # The column accumulators are NOT zeroed here (the caller manages that); the
    # kernel only does `+=` scatters, matching the host loop's accumulation.
    _launch!(_gap_p2c_kernel!,
        cnveg_cf.gap_mortality_c_to_litr_c_col,
        cnveg_cf.gap_mortality_c_to_cwdc_col,
        cnveg_nf.gap_mortality_n_to_litr_n_col,
        cnveg_nf.gap_mortality_n_to_cwdn_col,
        mask_soilp, patch.itype, patch.column, patch.wtcol,
        pftcon.lf_f, pftcon.fr_f,
        leaf_prof, froot_prof, croot_prof, stem_prof,
        cnveg_cf.m_leafc_to_litter_patch, cnveg_cf.m_frootc_to_litter_patch,
        cnveg_cf.m_livestemc_to_litter_patch, cnveg_cf.m_deadstemc_to_litter_patch,
        cnveg_cf.m_livecrootc_to_litter_patch, cnveg_cf.m_deadcrootc_to_litter_patch,
        cnveg_cf.m_leafc_storage_to_litter_patch, cnveg_cf.m_frootc_storage_to_litter_patch,
        cnveg_cf.m_livestemc_storage_to_litter_patch, cnveg_cf.m_deadstemc_storage_to_litter_patch,
        cnveg_cf.m_livecrootc_storage_to_litter_patch, cnveg_cf.m_deadcrootc_storage_to_litter_patch,
        cnveg_cf.m_gresp_storage_to_litter_patch,
        cnveg_cf.m_leafc_xfer_to_litter_patch, cnveg_cf.m_frootc_xfer_to_litter_patch,
        cnveg_cf.m_livestemc_xfer_to_litter_patch, cnveg_cf.m_deadstemc_xfer_to_litter_patch,
        cnveg_cf.m_livecrootc_xfer_to_litter_patch, cnveg_cf.m_deadcrootc_xfer_to_litter_patch,
        cnveg_cf.m_gresp_xfer_to_litter_patch,
        cnveg_nf.m_leafn_to_litter_patch, cnveg_nf.m_frootn_to_litter_patch,
        cnveg_nf.m_livestemn_to_litter_patch, cnveg_nf.m_deadstemn_to_litter_patch,
        cnveg_nf.m_livecrootn_to_litter_patch, cnveg_nf.m_deadcrootn_to_litter_patch,
        cnveg_nf.m_retransn_to_litter_patch,
        cnveg_nf.m_leafn_storage_to_litter_patch, cnveg_nf.m_frootn_storage_to_litter_patch,
        cnveg_nf.m_livestemn_storage_to_litter_patch, cnveg_nf.m_deadstemn_storage_to_litter_patch,
        cnveg_nf.m_livecrootn_storage_to_litter_patch, cnveg_nf.m_deadcrootn_storage_to_litter_patch,
        cnveg_nf.m_leafn_xfer_to_litter_patch, cnveg_nf.m_frootn_xfer_to_litter_patch,
        cnveg_nf.m_livestemn_xfer_to_litter_patch, cnveg_nf.m_deadstemn_xfer_to_litter_patch,
        cnveg_nf.m_livecrootn_xfer_to_litter_patch, cnveg_nf.m_deadcrootn_xfer_to_litter_patch,
        nlevdecomp, i_litr_min, i_litr_max, i_met_lit)

    return nothing
end
