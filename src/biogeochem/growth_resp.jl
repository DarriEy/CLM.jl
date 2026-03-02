# ==========================================================================
# Ported from: src/biogeochem/CNGRespMod.F90
# Growth respiration fluxes for coupled carbon-nitrogen code
# ==========================================================================

# ---------------------------------------------------------------------------
# PFT constants needed by growth respiration (from pftcon)
# ---------------------------------------------------------------------------

"""
    PftConGrowthResp

PFT-level parameters used by the growth respiration routine. Contains a
subset of fields from `pftconMod` that are referenced in `CNGRespMod.F90`.
"""
Base.@kwdef mutable struct PftConGrowthResp
    woody   ::Vector{Float64} = Float64[]   # binary woody flag (1=woody, 0=not woody)
    grperc  ::Vector{Float64} = Float64[]   # growth respiration parameter
    grpnow  ::Vector{Float64} = Float64[]   # growth respiration parameter (fraction now vs storage)
end

# ---------------------------------------------------------------------------
# cn_gresp! — Growth respiration
# ---------------------------------------------------------------------------

"""
    cn_gresp!(mask_soilp, bounds, pftcon, patch, cnveg_cf;
              npcropmin=17, nrepr=NREPR)

On the radiation time step, calculate growth respiration fluxes for all
vegetation components (leaf, fine root, live/dead stem, live/dead coarse root,
and crop reproductive organs).

Ported from `CNGResp` in `CNGRespMod.F90`.
"""
function cn_gresp!(mask_soilp::BitVector, bounds::UnitRange{Int},
                   pftcon::PftConGrowthResp,
                   patch::PatchData,
                   cnveg_cf::CNVegCarbonFluxData;
                   npcropmin::Int=17,
                   nrepr::Int=NREPR)

    # --- Aliases (matching Fortran associate block) ---
    ivt = patch.itype

    woody  = pftcon.woody
    grperc = pftcon.grperc
    grpnow = pftcon.grpnow

    # Input fluxes from allocation
    cpool_to_leafc              = cnveg_cf.cpool_to_leafc_patch
    cpool_to_leafc_storage      = cnveg_cf.cpool_to_leafc_storage_patch
    cpool_to_frootc             = cnveg_cf.cpool_to_frootc_patch
    cpool_to_frootc_storage     = cnveg_cf.cpool_to_frootc_storage_patch
    cpool_to_livestemc          = cnveg_cf.cpool_to_livestemc_patch
    cpool_to_livestemc_storage  = cnveg_cf.cpool_to_livestemc_storage_patch
    cpool_to_deadstemc          = cnveg_cf.cpool_to_deadstemc_patch
    cpool_to_deadstemc_storage  = cnveg_cf.cpool_to_deadstemc_storage_patch
    cpool_to_livecrootc         = cnveg_cf.cpool_to_livecrootc_patch
    cpool_to_livecrootc_storage = cnveg_cf.cpool_to_livecrootc_storage_patch
    cpool_to_deadcrootc         = cnveg_cf.cpool_to_deadcrootc_patch
    cpool_to_deadcrootc_storage = cnveg_cf.cpool_to_deadcrootc_storage_patch
    cpool_to_reproductivec         = cnveg_cf.cpool_to_reproductivec_patch
    cpool_to_reproductivec_storage = cnveg_cf.cpool_to_reproductivec_storage_patch

    # Input fluxes from transfer growth
    leafc_xfer_to_leafc           = cnveg_cf.leafc_xfer_to_leafc_patch
    frootc_xfer_to_frootc         = cnveg_cf.frootc_xfer_to_frootc_patch
    livestemc_xfer_to_livestemc   = cnveg_cf.livestemc_xfer_to_livestemc_patch
    deadstemc_xfer_to_deadstemc   = cnveg_cf.deadstemc_xfer_to_deadstemc_patch
    livecrootc_xfer_to_livecrootc = cnveg_cf.livecrootc_xfer_to_livecrootc_patch
    deadcrootc_xfer_to_deadcrootc = cnveg_cf.deadcrootc_xfer_to_deadcrootc_patch
    reproductivec_xfer_to_reproductivec = cnveg_cf.reproductivec_xfer_to_reproductivec_patch

    # Output growth respiration fluxes
    cpool_leaf_gr              = cnveg_cf.cpool_leaf_gr_patch
    cpool_leaf_storage_gr      = cnveg_cf.cpool_leaf_storage_gr_patch
    transfer_leaf_gr           = cnveg_cf.transfer_leaf_gr_patch
    cpool_froot_gr             = cnveg_cf.cpool_froot_gr_patch
    cpool_froot_storage_gr     = cnveg_cf.cpool_froot_storage_gr_patch
    transfer_froot_gr          = cnveg_cf.transfer_froot_gr_patch
    cpool_livestem_gr          = cnveg_cf.cpool_livestem_gr_patch
    cpool_livestem_storage_gr  = cnveg_cf.cpool_livestem_storage_gr_patch
    transfer_livestem_gr       = cnveg_cf.transfer_livestem_gr_patch
    cpool_deadstem_gr          = cnveg_cf.cpool_deadstem_gr_patch
    cpool_deadstem_storage_gr  = cnveg_cf.cpool_deadstem_storage_gr_patch
    transfer_deadstem_gr       = cnveg_cf.transfer_deadstem_gr_patch
    cpool_livecroot_gr         = cnveg_cf.cpool_livecroot_gr_patch
    cpool_livecroot_storage_gr = cnveg_cf.cpool_livecroot_storage_gr_patch
    transfer_livecroot_gr      = cnveg_cf.transfer_livecroot_gr_patch
    cpool_deadcroot_gr         = cnveg_cf.cpool_deadcroot_gr_patch
    cpool_deadcroot_storage_gr = cnveg_cf.cpool_deadcroot_storage_gr_patch
    transfer_deadcroot_gr      = cnveg_cf.transfer_deadcroot_gr_patch
    cpool_reproductive_gr         = cnveg_cf.cpool_reproductive_gr_patch
    cpool_reproductive_storage_gr = cnveg_cf.cpool_reproductive_storage_gr_patch
    transfer_reproductive_gr      = cnveg_cf.transfer_reproductive_gr_patch

    # --- Patch loop ---
    for p in bounds
        mask_soilp[p] || continue

        # Respiration factors (all 1.0 — preserved from Fortran for traceability)
        respfact_leaf              = 1.0
        respfact_froot             = 1.0
        respfact_livecroot         = 1.0
        respfact_livestem          = 1.0
        respfact_livecroot         = 1.0
        respfact_livestem          = 1.0
        respfact_leaf_storage      = 1.0
        respfact_froot_storage     = 1.0
        respfact_livecroot_storage = 1.0
        respfact_livestem_storage  = 1.0
        respfact_livecroot_storage = 1.0
        respfact_livestem_storage  = 1.0

        # --- Crop-specific growth respiration (ivt >= npcropmin) ---
        if ivt[p] >= npcropmin
            cpool_livestem_gr[p] = cpool_to_livestemc[p] * grperc[ivt[p]] * respfact_livestem

            cpool_livestem_storage_gr[p] = cpool_to_livestemc_storage[p] * grperc[ivt[p]] * grpnow[ivt[p]] *
                respfact_livestem_storage

            transfer_livestem_gr[p] = livestemc_xfer_to_livestemc[p] * grperc[ivt[p]] * (1.0 - grpnow[ivt[p]]) *
                respfact_livestem_storage

            for k in 1:nrepr
                cpool_reproductive_gr[p, k] =
                    cpool_to_reproductivec[p, k] * grperc[ivt[p]]
                cpool_reproductive_storage_gr[p, k] =
                    cpool_to_reproductivec_storage[p, k] * grperc[ivt[p]] * grpnow[ivt[p]]
                transfer_reproductive_gr[p, k] =
                    reproductivec_xfer_to_reproductivec[p, k] * grperc[ivt[p]] * (1.0 - grpnow[ivt[p]])
            end
        end

        # --- Leaf and fine root growth respiration (all PFTs) ---
        cpool_leaf_gr[p] = cpool_to_leafc[p] * grperc[ivt[p]] * respfact_leaf

        cpool_leaf_storage_gr[p] = cpool_to_leafc_storage[p] * grperc[ivt[p]] * grpnow[ivt[p]] * respfact_leaf_storage

        transfer_leaf_gr[p] = leafc_xfer_to_leafc[p] * grperc[ivt[p]] * (1.0 - grpnow[ivt[p]]) * respfact_leaf_storage

        # Note: respfact_froot appears twice (matches Fortran original)
        cpool_froot_gr[p] = cpool_to_frootc[p] * grperc[ivt[p]] * respfact_froot * respfact_froot

        cpool_froot_storage_gr[p] = cpool_to_frootc_storage[p] * grperc[ivt[p]] * grpnow[ivt[p]] * respfact_froot_storage

        transfer_froot_gr[p] = frootc_xfer_to_frootc[p] * grperc[ivt[p]] * (1.0 - grpnow[ivt[p]]) * respfact_froot_storage

        # --- Woody PFT growth respiration ---
        if woody[ivt[p]] == 1.0
            cpool_livestem_gr[p] = cpool_to_livestemc[p] * grperc[ivt[p]] * respfact_livestem

            cpool_livestem_storage_gr[p] = cpool_to_livestemc_storage[p] * grperc[ivt[p]] * grpnow[ivt[p]] *
                respfact_livestem_storage

            transfer_livestem_gr[p] = livestemc_xfer_to_livestemc[p] * grperc[ivt[p]] * (1.0 - grpnow[ivt[p]]) *
                respfact_livestem_storage

            cpool_deadstem_gr[p] = cpool_to_deadstemc[p] * grperc[ivt[p]]

            cpool_deadstem_storage_gr[p] = cpool_to_deadstemc_storage[p] * grperc[ivt[p]] * grpnow[ivt[p]]

            transfer_deadstem_gr[p] = deadstemc_xfer_to_deadstemc[p] * grperc[ivt[p]] * (1.0 - grpnow[ivt[p]])

            cpool_livecroot_gr[p] = cpool_to_livecrootc[p] * grperc[ivt[p]] * respfact_livecroot

            cpool_livecroot_storage_gr[p] = cpool_to_livecrootc_storage[p] * grperc[ivt[p]] * grpnow[ivt[p]] *
                respfact_livecroot_storage

            transfer_livecroot_gr[p] = livecrootc_xfer_to_livecrootc[p] * grperc[ivt[p]] * (1.0 - grpnow[ivt[p]]) *
                respfact_livecroot_storage

            cpool_deadcroot_gr[p] = cpool_to_deadcrootc[p] * grperc[ivt[p]]

            cpool_deadcroot_storage_gr[p] = cpool_to_deadcrootc_storage[p] * grperc[ivt[p]] * grpnow[ivt[p]]

            transfer_deadcroot_gr[p] = deadcrootc_xfer_to_deadcrootc[p] * grperc[ivt[p]] * (1.0 - grpnow[ivt[p]])
        end
    end

    return nothing
end
