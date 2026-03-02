@testset "Growth Respiration" begin

    # =====================================================================
    # PftConGrowthResp struct
    # =====================================================================
    @testset "PftConGrowthResp defaults" begin
        pfc = CLM.PftConGrowthResp()
        @test isempty(pfc.woody)
        @test isempty(pfc.grperc)
        @test isempty(pfc.grpnow)
    end

    # =====================================================================
    # Helper: create test data for a single woody patch
    # =====================================================================
    function make_woody_gresp_data(;
            grperc_val=0.3, grpnow_val=0.5,
            cpool_to_leafc=1.0, cpool_to_leafc_storage=0.5,
            leafc_xfer_to_leafc=0.2,
            cpool_to_frootc=0.8, cpool_to_frootc_storage=0.4,
            frootc_xfer_to_frootc=0.15,
            cpool_to_livestemc=0.6, cpool_to_livestemc_storage=0.3,
            livestemc_xfer_to_livestemc=0.1,
            cpool_to_deadstemc=0.5, cpool_to_deadstemc_storage=0.25,
            deadstemc_xfer_to_deadstemc=0.08,
            cpool_to_livecrootc=0.4, cpool_to_livecrootc_storage=0.2,
            livecrootc_xfer_to_livecrootc=0.06,
            cpool_to_deadcrootc=0.3, cpool_to_deadcrootc_storage=0.15,
            deadcrootc_xfer_to_deadcrootc=0.05)

        np = 1
        nrepr = CLM.NREPR

        # PFT constants: woody tree (ivt=1)
        pftcon = CLM.PftConGrowthResp(
            woody  = [1.0],
            grperc = [grperc_val],
            grpnow = [grpnow_val],
        )

        patch = CLM.PatchData()
        CLM.patch_init!(patch, np)
        patch.itype[1] = 1

        cnveg_cf = CLM.CNVegCarbonFluxData()
        CLM.cnveg_carbon_flux_init!(cnveg_cf, np, 1, 1)

        # Set input fluxes
        cnveg_cf.cpool_to_leafc_patch[1]              = cpool_to_leafc
        cnveg_cf.cpool_to_leafc_storage_patch[1]      = cpool_to_leafc_storage
        cnveg_cf.leafc_xfer_to_leafc_patch[1]         = leafc_xfer_to_leafc
        cnveg_cf.cpool_to_frootc_patch[1]             = cpool_to_frootc
        cnveg_cf.cpool_to_frootc_storage_patch[1]     = cpool_to_frootc_storage
        cnveg_cf.frootc_xfer_to_frootc_patch[1]       = frootc_xfer_to_frootc
        cnveg_cf.cpool_to_livestemc_patch[1]          = cpool_to_livestemc
        cnveg_cf.cpool_to_livestemc_storage_patch[1]  = cpool_to_livestemc_storage
        cnveg_cf.livestemc_xfer_to_livestemc_patch[1] = livestemc_xfer_to_livestemc
        cnveg_cf.cpool_to_deadstemc_patch[1]          = cpool_to_deadstemc
        cnveg_cf.cpool_to_deadstemc_storage_patch[1]  = cpool_to_deadstemc_storage
        cnveg_cf.deadstemc_xfer_to_deadstemc_patch[1] = deadstemc_xfer_to_deadstemc
        cnveg_cf.cpool_to_livecrootc_patch[1]         = cpool_to_livecrootc
        cnveg_cf.cpool_to_livecrootc_storage_patch[1] = cpool_to_livecrootc_storage
        cnveg_cf.livecrootc_xfer_to_livecrootc_patch[1] = livecrootc_xfer_to_livecrootc
        cnveg_cf.cpool_to_deadcrootc_patch[1]         = cpool_to_deadcrootc
        cnveg_cf.cpool_to_deadcrootc_storage_patch[1] = cpool_to_deadcrootc_storage
        cnveg_cf.deadcrootc_xfer_to_deadcrootc_patch[1] = deadcrootc_xfer_to_deadcrootc

        # Zero out reproductive fluxes (non-crop)
        for k in 1:nrepr
            cnveg_cf.cpool_to_reproductivec_patch[1, k] = 0.0
            cnveg_cf.cpool_to_reproductivec_storage_patch[1, k] = 0.0
            cnveg_cf.reproductivec_xfer_to_reproductivec_patch[1, k] = 0.0
        end

        mask = BitVector([true])
        bounds = 1:np

        return pftcon, patch, cnveg_cf, mask, bounds, nrepr
    end

    # =====================================================================
    # Test: woody PFT growth respiration
    # =====================================================================
    @testset "Woody PFT growth respiration" begin
        grperc_val = 0.3
        grpnow_val = 0.5

        pftcon, patch, cnveg_cf, mask, bounds, nrepr = make_woody_gresp_data(
            grperc_val=grperc_val, grpnow_val=grpnow_val)

        CLM.cn_gresp!(mask, bounds, pftcon, patch, cnveg_cf;
                      npcropmin=17, nrepr=nrepr)

        # Leaf growth respiration (respfact_leaf = 1.0)
        @test cnveg_cf.cpool_leaf_gr_patch[1] ≈ 1.0 * grperc_val * 1.0
        @test cnveg_cf.cpool_leaf_storage_gr_patch[1] ≈ 0.5 * grperc_val * grpnow_val * 1.0
        @test cnveg_cf.transfer_leaf_gr_patch[1] ≈ 0.2 * grperc_val * (1.0 - grpnow_val) * 1.0

        # Fine root growth respiration (respfact_froot appears twice!)
        @test cnveg_cf.cpool_froot_gr_patch[1] ≈ 0.8 * grperc_val * 1.0 * 1.0
        @test cnveg_cf.cpool_froot_storage_gr_patch[1] ≈ 0.4 * grperc_val * grpnow_val * 1.0
        @test cnveg_cf.transfer_froot_gr_patch[1] ≈ 0.15 * grperc_val * (1.0 - grpnow_val) * 1.0

        # Live stem growth respiration (woody branch overrides crop)
        @test cnveg_cf.cpool_livestem_gr_patch[1] ≈ 0.6 * grperc_val * 1.0
        @test cnveg_cf.cpool_livestem_storage_gr_patch[1] ≈ 0.3 * grperc_val * grpnow_val * 1.0
        @test cnveg_cf.transfer_livestem_gr_patch[1] ≈ 0.1 * grperc_val * (1.0 - grpnow_val) * 1.0

        # Dead stem growth respiration
        @test cnveg_cf.cpool_deadstem_gr_patch[1] ≈ 0.5 * grperc_val
        @test cnveg_cf.cpool_deadstem_storage_gr_patch[1] ≈ 0.25 * grperc_val * grpnow_val
        @test cnveg_cf.transfer_deadstem_gr_patch[1] ≈ 0.08 * grperc_val * (1.0 - grpnow_val)

        # Live coarse root growth respiration
        @test cnveg_cf.cpool_livecroot_gr_patch[1] ≈ 0.4 * grperc_val * 1.0
        @test cnveg_cf.cpool_livecroot_storage_gr_patch[1] ≈ 0.2 * grperc_val * grpnow_val * 1.0
        @test cnveg_cf.transfer_livecroot_gr_patch[1] ≈ 0.06 * grperc_val * (1.0 - grpnow_val) * 1.0

        # Dead coarse root growth respiration
        @test cnveg_cf.cpool_deadcroot_gr_patch[1] ≈ 0.3 * grperc_val
        @test cnveg_cf.cpool_deadcroot_storage_gr_patch[1] ≈ 0.15 * grperc_val * grpnow_val
        @test cnveg_cf.transfer_deadcroot_gr_patch[1] ≈ 0.05 * grperc_val * (1.0 - grpnow_val)
    end

    # =====================================================================
    # Test: non-woody, non-crop PFT (grass)
    # =====================================================================
    @testset "Non-woody grass PFT" begin
        np = 1
        nrepr = CLM.NREPR
        grperc_val = 0.3
        grpnow_val = 0.5

        pftcon = CLM.PftConGrowthResp(
            woody  = [0.0],
            grperc = [grperc_val],
            grpnow = [grpnow_val],
        )

        patch = CLM.PatchData()
        CLM.patch_init!(patch, np)
        patch.itype[1] = 1  # non-crop PFT index < npcropmin

        cnveg_cf = CLM.CNVegCarbonFluxData()
        CLM.cnveg_carbon_flux_init!(cnveg_cf, np, 1, 1)

        # Set some allocation fluxes
        cnveg_cf.cpool_to_leafc_patch[1]          = 2.0
        cnveg_cf.cpool_to_leafc_storage_patch[1]  = 1.0
        cnveg_cf.leafc_xfer_to_leafc_patch[1]     = 0.5
        cnveg_cf.cpool_to_frootc_patch[1]         = 1.5
        cnveg_cf.cpool_to_frootc_storage_patch[1] = 0.8
        cnveg_cf.frootc_xfer_to_frootc_patch[1]   = 0.3

        # Zero out stem/root fluxes that woody branch would set
        cnveg_cf.cpool_to_livestemc_patch[1]          = 0.0
        cnveg_cf.cpool_to_livestemc_storage_patch[1]  = 0.0
        cnveg_cf.livestemc_xfer_to_livestemc_patch[1] = 0.0
        cnveg_cf.cpool_to_deadstemc_patch[1]          = 0.0
        cnveg_cf.cpool_to_deadstemc_storage_patch[1]  = 0.0
        cnveg_cf.deadstemc_xfer_to_deadstemc_patch[1] = 0.0
        cnveg_cf.cpool_to_livecrootc_patch[1]         = 0.0
        cnveg_cf.cpool_to_livecrootc_storage_patch[1] = 0.0
        cnveg_cf.livecrootc_xfer_to_livecrootc_patch[1] = 0.0
        cnveg_cf.cpool_to_deadcrootc_patch[1]         = 0.0
        cnveg_cf.cpool_to_deadcrootc_storage_patch[1] = 0.0
        cnveg_cf.deadcrootc_xfer_to_deadcrootc_patch[1] = 0.0

        for k in 1:nrepr
            cnveg_cf.cpool_to_reproductivec_patch[1, k] = 0.0
            cnveg_cf.cpool_to_reproductivec_storage_patch[1, k] = 0.0
            cnveg_cf.reproductivec_xfer_to_reproductivec_patch[1, k] = 0.0
        end

        mask = BitVector([true])
        bounds = 1:np

        CLM.cn_gresp!(mask, bounds, pftcon, patch, cnveg_cf;
                      npcropmin=17, nrepr=nrepr)

        # Leaf/froot GR should be computed
        @test cnveg_cf.cpool_leaf_gr_patch[1] ≈ 2.0 * grperc_val
        @test cnveg_cf.cpool_froot_gr_patch[1] ≈ 1.5 * grperc_val  # respfact_froot * respfact_froot = 1.0

        # Woody-only fields should NOT have been set (remain NaN from init)
        @test isnan(cnveg_cf.cpool_deadstem_gr_patch[1])
        @test isnan(cnveg_cf.cpool_deadcroot_gr_patch[1])
    end

    # =====================================================================
    # Test: crop PFT (ivt >= npcropmin)
    # =====================================================================
    @testset "Crop PFT growth respiration" begin
        np = 1
        nrepr = CLM.NREPR
        grperc_val = 0.25
        grpnow_val = 0.6
        npcropmin = 17

        # PFT index 17 (first crop)
        pftcon = CLM.PftConGrowthResp()
        pftcon.woody  = zeros(npcropmin)
        pftcon.grperc = fill(grperc_val, npcropmin)
        pftcon.grpnow = fill(grpnow_val, npcropmin)

        patch = CLM.PatchData()
        CLM.patch_init!(patch, np)
        patch.itype[1] = npcropmin

        cnveg_cf = CLM.CNVegCarbonFluxData()
        CLM.cnveg_carbon_flux_init!(cnveg_cf, np, 1, 1)

        # Set crop-specific allocation fluxes
        cnveg_cf.cpool_to_livestemc_patch[1]          = 0.7
        cnveg_cf.cpool_to_livestemc_storage_patch[1]  = 0.35
        cnveg_cf.livestemc_xfer_to_livestemc_patch[1] = 0.12
        cnveg_cf.cpool_to_leafc_patch[1]              = 1.0
        cnveg_cf.cpool_to_leafc_storage_patch[1]      = 0.5
        cnveg_cf.leafc_xfer_to_leafc_patch[1]         = 0.2
        cnveg_cf.cpool_to_frootc_patch[1]             = 0.8
        cnveg_cf.cpool_to_frootc_storage_patch[1]     = 0.4
        cnveg_cf.frootc_xfer_to_frootc_patch[1]       = 0.15

        # Zero out non-crop fluxes
        cnveg_cf.cpool_to_deadstemc_patch[1]          = 0.0
        cnveg_cf.cpool_to_deadstemc_storage_patch[1]  = 0.0
        cnveg_cf.deadstemc_xfer_to_deadstemc_patch[1] = 0.0
        cnveg_cf.cpool_to_livecrootc_patch[1]         = 0.0
        cnveg_cf.cpool_to_livecrootc_storage_patch[1] = 0.0
        cnveg_cf.livecrootc_xfer_to_livecrootc_patch[1] = 0.0
        cnveg_cf.cpool_to_deadcrootc_patch[1]         = 0.0
        cnveg_cf.cpool_to_deadcrootc_storage_patch[1] = 0.0
        cnveg_cf.deadcrootc_xfer_to_deadcrootc_patch[1] = 0.0

        # Reproductive fluxes
        for k in 1:nrepr
            cnveg_cf.cpool_to_reproductivec_patch[1, k] = 0.5
            cnveg_cf.cpool_to_reproductivec_storage_patch[1, k] = 0.25
            cnveg_cf.reproductivec_xfer_to_reproductivec_patch[1, k] = 0.1
        end

        mask = BitVector([true])
        bounds = 1:np

        CLM.cn_gresp!(mask, bounds, pftcon, patch, cnveg_cf;
                      npcropmin=npcropmin, nrepr=nrepr)

        # Crop livestem GR (set in crop block)
        @test cnveg_cf.cpool_livestem_gr_patch[1] ≈ 0.7 * grperc_val
        @test cnveg_cf.cpool_livestem_storage_gr_patch[1] ≈ 0.35 * grperc_val * grpnow_val
        @test cnveg_cf.transfer_livestem_gr_patch[1] ≈ 0.12 * grperc_val * (1.0 - grpnow_val)

        # Reproductive GR
        for k in 1:nrepr
            @test cnveg_cf.cpool_reproductive_gr_patch[1, k] ≈ 0.5 * grperc_val
            @test cnveg_cf.cpool_reproductive_storage_gr_patch[1, k] ≈ 0.25 * grperc_val * grpnow_val
            @test cnveg_cf.transfer_reproductive_gr_patch[1, k] ≈ 0.1 * grperc_val * (1.0 - grpnow_val)
        end

        # Leaf/froot should also be set (all PFTs)
        @test cnveg_cf.cpool_leaf_gr_patch[1] ≈ 1.0 * grperc_val
        @test cnveg_cf.cpool_froot_gr_patch[1] ≈ 0.8 * grperc_val
    end

    # =====================================================================
    # Test: masked-out patches are skipped
    # =====================================================================
    @testset "Masked-out patches skipped" begin
        np = 2
        nrepr = CLM.NREPR

        pftcon = CLM.PftConGrowthResp(
            woody  = [1.0],
            grperc = [0.3],
            grpnow = [0.5],
        )

        patch = CLM.PatchData()
        CLM.patch_init!(patch, np)
        patch.itype[1] = 1
        patch.itype[2] = 1

        cnveg_cf = CLM.CNVegCarbonFluxData()
        CLM.cnveg_carbon_flux_init!(cnveg_cf, np, 1, 1)

        cnveg_cf.cpool_to_leafc_patch[1] = 1.0
        cnveg_cf.cpool_to_leafc_patch[2] = 2.0
        cnveg_cf.cpool_to_leafc_storage_patch .= 0.0
        cnveg_cf.leafc_xfer_to_leafc_patch .= 0.0
        cnveg_cf.cpool_to_frootc_patch .= 0.0
        cnveg_cf.cpool_to_frootc_storage_patch .= 0.0
        cnveg_cf.frootc_xfer_to_frootc_patch .= 0.0
        cnveg_cf.cpool_to_livestemc_patch .= 0.0
        cnveg_cf.cpool_to_livestemc_storage_patch .= 0.0
        cnveg_cf.livestemc_xfer_to_livestemc_patch .= 0.0
        cnveg_cf.cpool_to_deadstemc_patch .= 0.0
        cnveg_cf.cpool_to_deadstemc_storage_patch .= 0.0
        cnveg_cf.deadstemc_xfer_to_deadstemc_patch .= 0.0
        cnveg_cf.cpool_to_livecrootc_patch .= 0.0
        cnveg_cf.cpool_to_livecrootc_storage_patch .= 0.0
        cnveg_cf.livecrootc_xfer_to_livecrootc_patch .= 0.0
        cnveg_cf.cpool_to_deadcrootc_patch .= 0.0
        cnveg_cf.cpool_to_deadcrootc_storage_patch .= 0.0
        cnveg_cf.deadcrootc_xfer_to_deadcrootc_patch .= 0.0
        for k in 1:nrepr
            cnveg_cf.cpool_to_reproductivec_patch[:, k] .= 0.0
            cnveg_cf.cpool_to_reproductivec_storage_patch[:, k] .= 0.0
            cnveg_cf.reproductivec_xfer_to_reproductivec_patch[:, k] .= 0.0
        end

        # Only patch 1 is active
        mask = BitVector([true, false])
        bounds = 1:np

        CLM.cn_gresp!(mask, bounds, pftcon, patch, cnveg_cf;
                      npcropmin=17, nrepr=nrepr)

        # Patch 1 computed
        @test cnveg_cf.cpool_leaf_gr_patch[1] ≈ 1.0 * 0.3
        # Patch 2 should remain NaN (not touched)
        @test isnan(cnveg_cf.cpool_leaf_gr_patch[2])
    end

    # =====================================================================
    # Test: zero allocation fluxes → zero growth respiration
    # =====================================================================
    @testset "Zero allocation → zero GR" begin
        np = 1
        nrepr = CLM.NREPR

        pftcon = CLM.PftConGrowthResp(
            woody  = [1.0],
            grperc = [0.3],
            grpnow = [0.5],
        )

        patch = CLM.PatchData()
        CLM.patch_init!(patch, np)
        patch.itype[1] = 1

        cnveg_cf = CLM.CNVegCarbonFluxData()
        CLM.cnveg_carbon_flux_init!(cnveg_cf, np, 1, 1)

        # Set all input fluxes to zero
        cnveg_cf.cpool_to_leafc_patch[1]              = 0.0
        cnveg_cf.cpool_to_leafc_storage_patch[1]      = 0.0
        cnveg_cf.leafc_xfer_to_leafc_patch[1]         = 0.0
        cnveg_cf.cpool_to_frootc_patch[1]             = 0.0
        cnveg_cf.cpool_to_frootc_storage_patch[1]     = 0.0
        cnveg_cf.frootc_xfer_to_frootc_patch[1]       = 0.0
        cnveg_cf.cpool_to_livestemc_patch[1]          = 0.0
        cnveg_cf.cpool_to_livestemc_storage_patch[1]  = 0.0
        cnveg_cf.livestemc_xfer_to_livestemc_patch[1] = 0.0
        cnveg_cf.cpool_to_deadstemc_patch[1]          = 0.0
        cnveg_cf.cpool_to_deadstemc_storage_patch[1]  = 0.0
        cnveg_cf.deadstemc_xfer_to_deadstemc_patch[1] = 0.0
        cnveg_cf.cpool_to_livecrootc_patch[1]         = 0.0
        cnveg_cf.cpool_to_livecrootc_storage_patch[1] = 0.0
        cnveg_cf.livecrootc_xfer_to_livecrootc_patch[1] = 0.0
        cnveg_cf.cpool_to_deadcrootc_patch[1]         = 0.0
        cnveg_cf.cpool_to_deadcrootc_storage_patch[1] = 0.0
        cnveg_cf.deadcrootc_xfer_to_deadcrootc_patch[1] = 0.0
        for k in 1:nrepr
            cnveg_cf.cpool_to_reproductivec_patch[1, k] = 0.0
            cnveg_cf.cpool_to_reproductivec_storage_patch[1, k] = 0.0
            cnveg_cf.reproductivec_xfer_to_reproductivec_patch[1, k] = 0.0
        end

        mask = BitVector([true])
        bounds = 1:np

        CLM.cn_gresp!(mask, bounds, pftcon, patch, cnveg_cf;
                      npcropmin=17, nrepr=nrepr)

        @test cnveg_cf.cpool_leaf_gr_patch[1] ≈ 0.0
        @test cnveg_cf.cpool_froot_gr_patch[1] ≈ 0.0
        @test cnveg_cf.cpool_livestem_gr_patch[1] ≈ 0.0
        @test cnveg_cf.cpool_deadstem_gr_patch[1] ≈ 0.0
        @test cnveg_cf.cpool_livecroot_gr_patch[1] ≈ 0.0
        @test cnveg_cf.cpool_deadcroot_gr_patch[1] ≈ 0.0
        @test cnveg_cf.transfer_leaf_gr_patch[1] ≈ 0.0
        @test cnveg_cf.transfer_livestem_gr_patch[1] ≈ 0.0
        @test cnveg_cf.transfer_deadcroot_gr_patch[1] ≈ 0.0
    end

end
