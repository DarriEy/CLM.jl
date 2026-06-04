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
Base.@kwdef mutable struct PftConGrowthResp{FT<:Real, V<:AbstractVector{FT}}
    woody   ::V = Float64[]   # binary woody flag (1=woody, 0=not woody)
    grperc  ::V = Float64[]   # growth respiration parameter
    grpnow  ::V = Float64[]   # growth respiration parameter (fraction now vs storage)
end
PftConGrowthResp{FT}(; kwargs...) where {FT<:Real} = PftConGrowthResp{FT, Vector{FT}}(; kwargs...)
Adapt.@adapt_structure PftConGrowthResp

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
Base.@kwdef struct _GrespOut{V,M}
    # per-patch (vector) growth-respiration flux outputs
    cpool_leaf_gr::V; cpool_leaf_storage_gr::V; transfer_leaf_gr::V
    cpool_froot_gr::V; cpool_froot_storage_gr::V; transfer_froot_gr::V
    cpool_livestem_gr::V; cpool_livestem_storage_gr::V; transfer_livestem_gr::V
    cpool_deadstem_gr::V; cpool_deadstem_storage_gr::V; transfer_deadstem_gr::V
    cpool_livecroot_gr::V; cpool_livecroot_storage_gr::V; transfer_livecroot_gr::V
    cpool_deadcroot_gr::V; cpool_deadcroot_storage_gr::V; transfer_deadcroot_gr::V
    # reproductive [p,k] (matrix) outputs
    cpool_reproductive_gr::M; cpool_reproductive_storage_gr::M; transfer_reproductive_gr::M
end
Adapt.@adapt_structure _GrespOut

Base.@kwdef struct _GrespIn{V,M}
    # per-patch (vector) input fluxes
    cpool_to_leafc::V; cpool_to_leafc_storage::V
    cpool_to_frootc::V; cpool_to_frootc_storage::V
    cpool_to_livestemc::V; cpool_to_livestemc_storage::V
    cpool_to_deadstemc::V; cpool_to_deadstemc_storage::V
    cpool_to_livecrootc::V; cpool_to_livecrootc_storage::V
    cpool_to_deadcrootc::V; cpool_to_deadcrootc_storage::V
    leafc_xfer_to_leafc::V; frootc_xfer_to_frootc::V
    livestemc_xfer_to_livestemc::V; deadstemc_xfer_to_deadstemc::V
    livecrootc_xfer_to_livecrootc::V; deadcrootc_xfer_to_deadcrootc::V
    # reproductive [p,k] (matrix) inputs
    cpool_to_reproductivec::M; cpool_to_reproductivec_storage::M
    reproductivec_xfer_to_reproductivec::M
end
Adapt.@adapt_structure _GrespIn

@kernel function _cn_gresp_kernel!(
        gout::_GrespOut, gin::_GrespIn,
        @Const(mask_soilp), @Const(ivt),
        @Const(woody), @Const(grperc), @Const(grpnow),
        npcropmin::Int, nrepr::Int, pmin::Int, pmax::Int)

    # --- alias output arrays to Fortran-named locals (body stays verbatim) ---
    cpool_leaf_gr                 = gout.cpool_leaf_gr
    cpool_leaf_storage_gr         = gout.cpool_leaf_storage_gr
    transfer_leaf_gr              = gout.transfer_leaf_gr
    cpool_froot_gr                = gout.cpool_froot_gr
    cpool_froot_storage_gr        = gout.cpool_froot_storage_gr
    transfer_froot_gr             = gout.transfer_froot_gr
    cpool_livestem_gr             = gout.cpool_livestem_gr
    cpool_livestem_storage_gr     = gout.cpool_livestem_storage_gr
    transfer_livestem_gr          = gout.transfer_livestem_gr
    cpool_deadstem_gr             = gout.cpool_deadstem_gr
    cpool_deadstem_storage_gr     = gout.cpool_deadstem_storage_gr
    transfer_deadstem_gr          = gout.transfer_deadstem_gr
    cpool_livecroot_gr            = gout.cpool_livecroot_gr
    cpool_livecroot_storage_gr    = gout.cpool_livecroot_storage_gr
    transfer_livecroot_gr         = gout.transfer_livecroot_gr
    cpool_deadcroot_gr            = gout.cpool_deadcroot_gr
    cpool_deadcroot_storage_gr    = gout.cpool_deadcroot_storage_gr
    transfer_deadcroot_gr         = gout.transfer_deadcroot_gr
    cpool_reproductive_gr         = gout.cpool_reproductive_gr
    cpool_reproductive_storage_gr = gout.cpool_reproductive_storage_gr
    transfer_reproductive_gr      = gout.transfer_reproductive_gr

    # --- alias input arrays to Fortran-named locals ---
    cpool_to_leafc                      = gin.cpool_to_leafc
    cpool_to_leafc_storage              = gin.cpool_to_leafc_storage
    cpool_to_frootc                     = gin.cpool_to_frootc
    cpool_to_frootc_storage             = gin.cpool_to_frootc_storage
    cpool_to_livestemc                  = gin.cpool_to_livestemc
    cpool_to_livestemc_storage          = gin.cpool_to_livestemc_storage
    cpool_to_deadstemc                  = gin.cpool_to_deadstemc
    cpool_to_deadstemc_storage          = gin.cpool_to_deadstemc_storage
    cpool_to_livecrootc                 = gin.cpool_to_livecrootc
    cpool_to_livecrootc_storage         = gin.cpool_to_livecrootc_storage
    cpool_to_deadcrootc                 = gin.cpool_to_deadcrootc
    cpool_to_deadcrootc_storage         = gin.cpool_to_deadcrootc_storage
    cpool_to_reproductivec              = gin.cpool_to_reproductivec
    cpool_to_reproductivec_storage      = gin.cpool_to_reproductivec_storage
    leafc_xfer_to_leafc                 = gin.leafc_xfer_to_leafc
    frootc_xfer_to_frootc               = gin.frootc_xfer_to_frootc
    livestemc_xfer_to_livestemc         = gin.livestemc_xfer_to_livestemc
    deadstemc_xfer_to_deadstemc         = gin.deadstemc_xfer_to_deadstemc
    livecrootc_xfer_to_livecrootc       = gin.livecrootc_xfer_to_livecrootc
    deadcrootc_xfer_to_deadcrootc       = gin.deadcrootc_xfer_to_deadcrootc
    reproductivec_xfer_to_reproductivec = gin.reproductivec_xfer_to_reproductivec

    T = eltype(cpool_leaf_gr)
    p = @index(Global)
    @inbounds if pmin <= p <= pmax && mask_soilp[p]

        # Respiration factors (all 1.0 — preserved from Fortran for traceability)
        respfact_leaf              = one(T)
        respfact_froot             = one(T)
        respfact_livecroot         = one(T)
        respfact_livestem          = one(T)
        respfact_livecroot         = one(T)
        respfact_livestem          = one(T)
        respfact_leaf_storage      = one(T)
        respfact_froot_storage     = one(T)
        respfact_livecroot_storage = one(T)
        respfact_livestem_storage  = one(T)
        respfact_livecroot_storage = one(T)
        respfact_livestem_storage  = one(T)

        iv = ivt[p] + 1

        # --- Crop-specific growth respiration (ivt >= npcropmin) ---
        if ivt[p] >= npcropmin
            cpool_livestem_gr[p] = cpool_to_livestemc[p] * grperc[iv] * respfact_livestem

            cpool_livestem_storage_gr[p] = cpool_to_livestemc_storage[p] * grperc[iv] * grpnow[iv] *
                respfact_livestem_storage

            transfer_livestem_gr[p] = livestemc_xfer_to_livestemc[p] * grperc[iv] * (one(T) - grpnow[iv]) *
                respfact_livestem_storage

            for k in 1:nrepr
                cpool_reproductive_gr[p, k] =
                    cpool_to_reproductivec[p, k] * grperc[iv]
                cpool_reproductive_storage_gr[p, k] =
                    cpool_to_reproductivec_storage[p, k] * grperc[iv] * grpnow[iv]
                transfer_reproductive_gr[p, k] =
                    reproductivec_xfer_to_reproductivec[p, k] * grperc[iv] * (one(T) - grpnow[iv])
            end
        end

        # --- Leaf and fine root growth respiration (all PFTs) ---
        cpool_leaf_gr[p] = cpool_to_leafc[p] * grperc[iv] * respfact_leaf

        cpool_leaf_storage_gr[p] = cpool_to_leafc_storage[p] * grperc[iv] * grpnow[iv] * respfact_leaf_storage

        transfer_leaf_gr[p] = leafc_xfer_to_leafc[p] * grperc[iv] * (one(T) - grpnow[iv]) * respfact_leaf_storage

        # Note: respfact_froot appears twice (matches Fortran original)
        cpool_froot_gr[p] = cpool_to_frootc[p] * grperc[iv] * respfact_froot * respfact_froot

        cpool_froot_storage_gr[p] = cpool_to_frootc_storage[p] * grperc[iv] * grpnow[iv] * respfact_froot_storage

        transfer_froot_gr[p] = frootc_xfer_to_frootc[p] * grperc[iv] * (one(T) - grpnow[iv]) * respfact_froot_storage

        # --- Woody PFT growth respiration ---
        if woody[iv] == one(T)
            cpool_livestem_gr[p] = cpool_to_livestemc[p] * grperc[iv] * respfact_livestem

            cpool_livestem_storage_gr[p] = cpool_to_livestemc_storage[p] * grperc[iv] * grpnow[iv] *
                respfact_livestem_storage

            transfer_livestem_gr[p] = livestemc_xfer_to_livestemc[p] * grperc[iv] * (one(T) - grpnow[iv]) *
                respfact_livestem_storage

            cpool_deadstem_gr[p] = cpool_to_deadstemc[p] * grperc[iv]

            cpool_deadstem_storage_gr[p] = cpool_to_deadstemc_storage[p] * grperc[iv] * grpnow[iv]

            transfer_deadstem_gr[p] = deadstemc_xfer_to_deadstemc[p] * grperc[iv] * (one(T) - grpnow[iv])

            cpool_livecroot_gr[p] = cpool_to_livecrootc[p] * grperc[iv] * respfact_livecroot

            cpool_livecroot_storage_gr[p] = cpool_to_livecrootc_storage[p] * grperc[iv] * grpnow[iv] *
                respfact_livecroot_storage

            transfer_livecroot_gr[p] = livecrootc_xfer_to_livecrootc[p] * grperc[iv] * (one(T) - grpnow[iv]) *
                respfact_livecroot_storage

            cpool_deadcroot_gr[p] = cpool_to_deadcrootc[p] * grperc[iv]

            cpool_deadcroot_storage_gr[p] = cpool_to_deadcrootc_storage[p] * grperc[iv] * grpnow[iv]

            transfer_deadcroot_gr[p] = deadcrootc_xfer_to_deadcrootc[p] * grperc[iv] * (one(T) - grpnow[iv])
        end
    end
end

function cn_gresp!(mask_soilp::AbstractVector{Bool}, bounds::UnitRange{Int},
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

    isempty(bounds) && return nothing

    # --- Group output growth respiration fluxes into a device-view struct ---
    gout = _GrespOut(;
        cpool_leaf_gr              = cnveg_cf.cpool_leaf_gr_patch,
        cpool_leaf_storage_gr      = cnveg_cf.cpool_leaf_storage_gr_patch,
        transfer_leaf_gr           = cnveg_cf.transfer_leaf_gr_patch,
        cpool_froot_gr             = cnveg_cf.cpool_froot_gr_patch,
        cpool_froot_storage_gr     = cnveg_cf.cpool_froot_storage_gr_patch,
        transfer_froot_gr          = cnveg_cf.transfer_froot_gr_patch,
        cpool_livestem_gr          = cnveg_cf.cpool_livestem_gr_patch,
        cpool_livestem_storage_gr  = cnveg_cf.cpool_livestem_storage_gr_patch,
        transfer_livestem_gr       = cnveg_cf.transfer_livestem_gr_patch,
        cpool_deadstem_gr          = cnveg_cf.cpool_deadstem_gr_patch,
        cpool_deadstem_storage_gr  = cnveg_cf.cpool_deadstem_storage_gr_patch,
        transfer_deadstem_gr       = cnveg_cf.transfer_deadstem_gr_patch,
        cpool_livecroot_gr         = cnveg_cf.cpool_livecroot_gr_patch,
        cpool_livecroot_storage_gr = cnveg_cf.cpool_livecroot_storage_gr_patch,
        transfer_livecroot_gr      = cnveg_cf.transfer_livecroot_gr_patch,
        cpool_deadcroot_gr         = cnveg_cf.cpool_deadcroot_gr_patch,
        cpool_deadcroot_storage_gr = cnveg_cf.cpool_deadcroot_storage_gr_patch,
        transfer_deadcroot_gr      = cnveg_cf.transfer_deadcroot_gr_patch,
        cpool_reproductive_gr         = cnveg_cf.cpool_reproductive_gr_patch,
        cpool_reproductive_storage_gr = cnveg_cf.cpool_reproductive_storage_gr_patch,
        transfer_reproductive_gr      = cnveg_cf.transfer_reproductive_gr_patch)

    # --- Group read-only input fluxes into a device-view struct ---
    gin = _GrespIn(;
        cpool_to_leafc                      = cnveg_cf.cpool_to_leafc_patch,
        cpool_to_leafc_storage              = cnveg_cf.cpool_to_leafc_storage_patch,
        cpool_to_frootc                     = cnveg_cf.cpool_to_frootc_patch,
        cpool_to_frootc_storage             = cnveg_cf.cpool_to_frootc_storage_patch,
        cpool_to_livestemc                  = cnveg_cf.cpool_to_livestemc_patch,
        cpool_to_livestemc_storage          = cnveg_cf.cpool_to_livestemc_storage_patch,
        cpool_to_deadstemc                  = cnveg_cf.cpool_to_deadstemc_patch,
        cpool_to_deadstemc_storage          = cnveg_cf.cpool_to_deadstemc_storage_patch,
        cpool_to_livecrootc                 = cnveg_cf.cpool_to_livecrootc_patch,
        cpool_to_livecrootc_storage         = cnveg_cf.cpool_to_livecrootc_storage_patch,
        cpool_to_deadcrootc                 = cnveg_cf.cpool_to_deadcrootc_patch,
        cpool_to_deadcrootc_storage         = cnveg_cf.cpool_to_deadcrootc_storage_patch,
        leafc_xfer_to_leafc                 = cnveg_cf.leafc_xfer_to_leafc_patch,
        frootc_xfer_to_frootc               = cnveg_cf.frootc_xfer_to_frootc_patch,
        livestemc_xfer_to_livestemc         = cnveg_cf.livestemc_xfer_to_livestemc_patch,
        deadstemc_xfer_to_deadstemc         = cnveg_cf.deadstemc_xfer_to_deadstemc_patch,
        livecrootc_xfer_to_livecrootc       = cnveg_cf.livecrootc_xfer_to_livecrootc_patch,
        deadcrootc_xfer_to_deadcrootc       = cnveg_cf.deadcrootc_xfer_to_deadcrootc_patch,
        cpool_to_reproductivec              = cnveg_cf.cpool_to_reproductivec_patch,
        cpool_to_reproductivec_storage      = cnveg_cf.cpool_to_reproductivec_storage_patch,
        reproductivec_xfer_to_reproductivec = cnveg_cf.reproductivec_xfer_to_reproductivec_patch)

    # Struct-first kernel: manual backend + synchronize (the bundle args carry no backend).
    backend = _kernel_backend(gout.cpool_leaf_gr)
    _cn_gresp_kernel!(backend)(gout, gin,
        mask_soilp, ivt,
        woody, grperc, grpnow,
        npcropmin, nrepr, first(bounds), last(bounds);
        ndrange = last(bounds))
    KA.synchronize(backend)

    return nothing
end
