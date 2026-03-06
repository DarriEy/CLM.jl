@testset "Gap Mortality" begin

    # ----------------------------------------------------------------
    # Helper: create minimal data structures for np patches, nc columns
    # ----------------------------------------------------------------
    function make_gap_mortality_data(; np=3, nc=2, nlevdecomp=1, n_litr=3)
        # --- Params ---
        params = CLM.GapMortalityParams(
            k_mort = 0.3,
            r_mort = fill(0.02, 21)  # 2% annual mortality for all PFTs (size ivt_max+1)
        )

        # --- PftCon ---
        pftcon = CLM.PftConGapMort(
            woody    = vcat(fill(1.0, 9), fill(0.0, 12)),  # 0-based PFTs 0-8 are woody (Julia indices 1-9)
            leafcn   = fill(25.0, 21),
            livewdcn = fill(50.0, 21),
            lf_f     = fill(1.0 / n_litr, 21, n_litr),     # equal fractions
            fr_f     = fill(1.0 / n_litr, 21, n_litr),
        )

        # --- DGVS ---
        dgvs = CLM.DgvsGapMortData(
            greffic_patch    = fill(0.5, np),
            heatstress_patch = fill(0.0, np),
            nind_patch       = fill(100.0, np),
        )

        # --- Patch ---
        patch = CLM.PatchData()
        patch.itype  = [2, 10, 5]        # PFT types (1-indexed): 2=woody, 10=non-woody, 5=woody
        patch.column = [1, 1, 2]          # column assignments
        patch.wtcol  = [0.5, 0.3, 0.2]   # weight relative to column

        # --- Canopy state ---
        canopystate = CLM.CanopyStateData()
        canopystate.laisun_patch = fill(2.0, np)
        canopystate.laisha_patch = fill(1.0, np)

        # --- Carbon state ---
        cnveg_cs = CLM.CNVegCarbonStateData()
        cnveg_cs.leafc_patch              = [10.0, 5.0, 8.0]
        cnveg_cs.frootc_patch             = [4.0, 2.0, 3.0]
        cnveg_cs.livestemc_patch          = [20.0, 0.0, 15.0]
        cnveg_cs.deadstemc_patch          = [50.0, 0.0, 40.0]
        cnveg_cs.livecrootc_patch         = [10.0, 0.0, 8.0]
        cnveg_cs.deadcrootc_patch         = [25.0, 0.0, 20.0]
        cnveg_cs.leafc_storage_patch      = [1.0, 0.5, 0.8]
        cnveg_cs.frootc_storage_patch     = [0.5, 0.2, 0.3]
        cnveg_cs.livestemc_storage_patch  = [2.0, 0.0, 1.5]
        cnveg_cs.deadstemc_storage_patch  = [3.0, 0.0, 2.5]
        cnveg_cs.livecrootc_storage_patch = [1.0, 0.0, 0.8]
        cnveg_cs.deadcrootc_storage_patch = [1.5, 0.0, 1.2]
        cnveg_cs.gresp_storage_patch      = [0.2, 0.1, 0.15]
        cnveg_cs.leafc_xfer_patch         = [0.5, 0.25, 0.4]
        cnveg_cs.frootc_xfer_patch        = [0.2, 0.1, 0.15]
        cnveg_cs.livestemc_xfer_patch     = [1.0, 0.0, 0.8]
        cnveg_cs.deadstemc_xfer_patch     = [1.5, 0.0, 1.2]
        cnveg_cs.livecrootc_xfer_patch    = [0.5, 0.0, 0.4]
        cnveg_cs.deadcrootc_xfer_patch    = [0.8, 0.0, 0.6]
        cnveg_cs.gresp_xfer_patch         = [0.1, 0.05, 0.08]

        # --- Carbon flux (zero-initialized) ---
        cnveg_cf = CLM.CNVegCarbonFluxData()
        cnveg_cf.m_leafc_to_litter_patch              = zeros(np)
        cnveg_cf.m_frootc_to_litter_patch             = zeros(np)
        cnveg_cf.m_livestemc_to_litter_patch          = zeros(np)
        cnveg_cf.m_deadstemc_to_litter_patch          = zeros(np)
        cnveg_cf.m_livecrootc_to_litter_patch         = zeros(np)
        cnveg_cf.m_deadcrootc_to_litter_patch         = zeros(np)
        cnveg_cf.m_leafc_storage_to_litter_patch      = zeros(np)
        cnveg_cf.m_frootc_storage_to_litter_patch     = zeros(np)
        cnveg_cf.m_livestemc_storage_to_litter_patch  = zeros(np)
        cnveg_cf.m_deadstemc_storage_to_litter_patch  = zeros(np)
        cnveg_cf.m_livecrootc_storage_to_litter_patch = zeros(np)
        cnveg_cf.m_deadcrootc_storage_to_litter_patch = zeros(np)
        cnveg_cf.m_gresp_storage_to_litter_patch      = zeros(np)
        cnveg_cf.m_leafc_xfer_to_litter_patch         = zeros(np)
        cnveg_cf.m_frootc_xfer_to_litter_patch        = zeros(np)
        cnveg_cf.m_livestemc_xfer_to_litter_patch     = zeros(np)
        cnveg_cf.m_deadstemc_xfer_to_litter_patch     = zeros(np)
        cnveg_cf.m_livecrootc_xfer_to_litter_patch    = zeros(np)
        cnveg_cf.m_deadcrootc_xfer_to_litter_patch    = zeros(np)
        cnveg_cf.m_gresp_xfer_to_litter_patch         = zeros(np)
        cnveg_cf.gap_mortality_c_to_litr_c_col        = zeros(nc, nlevdecomp, n_litr)
        cnveg_cf.gap_mortality_c_to_cwdc_col          = zeros(nc, nlevdecomp)

        # --- Nitrogen state ---
        cnveg_ns = CLM.CNVegNitrogenStateData()
        cnveg_ns.leafn_patch              = [0.4, 0.2, 0.32]
        cnveg_ns.frootn_patch             = [0.16, 0.08, 0.12]
        cnveg_ns.livestemn_patch          = [0.4, 0.0, 0.3]
        cnveg_ns.deadstemn_patch          = [0.5, 0.0, 0.4]
        cnveg_ns.livecrootn_patch         = [0.2, 0.0, 0.16]
        cnveg_ns.deadcrootn_patch         = [0.25, 0.0, 0.2]
        cnveg_ns.retransn_patch           = [0.1, 0.05, 0.08]
        cnveg_ns.leafn_storage_patch      = [0.04, 0.02, 0.032]
        cnveg_ns.frootn_storage_patch     = [0.02, 0.008, 0.012]
        cnveg_ns.livestemn_storage_patch  = [0.04, 0.0, 0.03]
        cnveg_ns.deadstemn_storage_patch  = [0.03, 0.0, 0.025]
        cnveg_ns.livecrootn_storage_patch = [0.02, 0.0, 0.016]
        cnveg_ns.deadcrootn_storage_patch = [0.015, 0.0, 0.012]
        cnveg_ns.leafn_xfer_patch         = [0.02, 0.01, 0.016]
        cnveg_ns.frootn_xfer_patch        = [0.008, 0.004, 0.006]
        cnveg_ns.livestemn_xfer_patch     = [0.02, 0.0, 0.016]
        cnveg_ns.deadstemn_xfer_patch     = [0.015, 0.0, 0.012]
        cnveg_ns.livecrootn_xfer_patch    = [0.01, 0.0, 0.008]
        cnveg_ns.deadcrootn_xfer_patch    = [0.008, 0.0, 0.006]

        # --- Nitrogen flux (zero-initialized) ---
        cnveg_nf = CLM.CNVegNitrogenFluxData()
        cnveg_nf.m_leafn_to_litter_patch              = zeros(np)
        cnveg_nf.m_frootn_to_litter_patch             = zeros(np)
        cnveg_nf.m_livestemn_to_litter_patch          = zeros(np)
        cnveg_nf.m_deadstemn_to_litter_patch          = zeros(np)
        cnveg_nf.m_livecrootn_to_litter_patch         = zeros(np)
        cnveg_nf.m_deadcrootn_to_litter_patch         = zeros(np)
        cnveg_nf.m_retransn_to_litter_patch           = zeros(np)
        cnveg_nf.m_leafn_storage_to_litter_patch      = zeros(np)
        cnveg_nf.m_frootn_storage_to_litter_patch     = zeros(np)
        cnveg_nf.m_livestemn_storage_to_litter_patch  = zeros(np)
        cnveg_nf.m_deadstemn_storage_to_litter_patch  = zeros(np)
        cnveg_nf.m_livecrootn_storage_to_litter_patch = zeros(np)
        cnveg_nf.m_deadcrootn_storage_to_litter_patch = zeros(np)
        cnveg_nf.m_leafn_xfer_to_litter_patch         = zeros(np)
        cnveg_nf.m_frootn_xfer_to_litter_patch        = zeros(np)
        cnveg_nf.m_livestemn_xfer_to_litter_patch     = zeros(np)
        cnveg_nf.m_deadstemn_xfer_to_litter_patch     = zeros(np)
        cnveg_nf.m_livecrootn_xfer_to_litter_patch    = zeros(np)
        cnveg_nf.m_deadcrootn_xfer_to_litter_patch    = zeros(np)
        cnveg_nf.gap_mortality_n_to_litr_n_col        = zeros(nc, nlevdecomp, n_litr)
        cnveg_nf.gap_mortality_n_to_cwdn_col          = zeros(nc, nlevdecomp)

        mask = BitVector([true, true, true])
        bounds = 1:np

        # Vertical profiles (uniform for simplicity)
        leaf_prof  = ones(np, nlevdecomp)
        froot_prof = ones(np, nlevdecomp)
        croot_prof = ones(np, nlevdecomp)
        stem_prof  = ones(np, nlevdecomp)

        return (; params, pftcon, dgvs, patch, canopystate, cnveg_cs, cnveg_cf,
                  cnveg_ns, cnveg_nf, mask, bounds, leaf_prof, froot_prof,
                  croot_prof, stem_prof, nc, nlevdecomp, n_litr)
    end

    # ================================================================
    @testset "GapMortalityParams construction" begin
        p = CLM.GapMortalityParams()
        @test p.k_mort == 0.3
        @test isempty(p.r_mort)

        p2 = CLM.GapMortalityParams(k_mort=0.5, r_mort=[0.01, 0.02])
        @test p2.k_mort == 0.5
        @test length(p2.r_mort) == 2
    end

    # ================================================================
    @testset "PftConGapMort construction" begin
        p = CLM.PftConGapMort()
        @test isempty(p.woody)
        @test isempty(p.leafcn)
    end

    # ================================================================
    @testset "DgvsGapMortData construction" begin
        d = CLM.DgvsGapMortData()
        @test isempty(d.greffic_patch)
        @test isempty(d.heatstress_patch)
        @test isempty(d.nind_patch)
    end

    # ================================================================
    @testset "cn_gap_mortality! basic" begin
        d = make_gap_mortality_data()

        CLM.cn_gap_mortality!(
            d.mask, d.bounds, d.params, d.pftcon, d.dgvs,
            d.patch, d.canopystate, d.cnveg_cs, d.cnveg_cf,
            d.cnveg_ns, d.cnveg_nf;
            dt = 1800.0,
            days_per_year = 365.0,
            use_cndv = false,
            use_matrixcn = false,
            spinup_state = 0,
            npcropmin = 17,
            spinup_factor_deadwood = 1.0,
        )

        # Expected fractional mortality rate per second
        am = 0.02  # annual mortality from params
        m  = am / (365.0 * CLM.SECSPDAY)

        # --- Carbon flux checks ---
        @test d.cnveg_cf.m_leafc_to_litter_patch[1] ≈ 10.0 * m
        @test d.cnveg_cf.m_leafc_to_litter_patch[2] ≈ 5.0 * m
        @test d.cnveg_cf.m_leafc_to_litter_patch[3] ≈ 8.0 * m

        @test d.cnveg_cf.m_frootc_to_litter_patch[1] ≈ 4.0 * m
        @test d.cnveg_cf.m_livestemc_to_litter_patch[1] ≈ 20.0 * m
        @test d.cnveg_cf.m_deadstemc_to_litter_patch[1] ≈ 50.0 * m  # spinup_factor=1.0
        @test d.cnveg_cf.m_livecrootc_to_litter_patch[1] ≈ 10.0 * m
        @test d.cnveg_cf.m_deadcrootc_to_litter_patch[1] ≈ 25.0 * m

        # Storage and transfer
        @test d.cnveg_cf.m_leafc_storage_to_litter_patch[1] ≈ 1.0 * m
        @test d.cnveg_cf.m_gresp_storage_to_litter_patch[1] ≈ 0.2 * m
        @test d.cnveg_cf.m_leafc_xfer_to_litter_patch[1] ≈ 0.5 * m

        # --- Nitrogen flux checks ---
        @test d.cnveg_nf.m_leafn_to_litter_patch[1] ≈ 0.4 * m
        @test d.cnveg_nf.m_frootn_to_litter_patch[1] ≈ 0.16 * m
        @test d.cnveg_nf.m_livestemn_to_litter_patch[1] ≈ 0.4 * m
        @test d.cnveg_nf.m_deadstemn_to_litter_patch[1] ≈ 0.5 * m  # spinup_state != 2
        @test d.cnveg_nf.m_deadcrootn_to_litter_patch[1] ≈ 0.25 * m

        # retransn (ivt[1]=2 < npcropmin=17, so should be set)
        @test d.cnveg_nf.m_retransn_to_litter_patch[1] ≈ 0.1 * m
    end

    # ================================================================
    @testset "cn_gap_mortality! with spinup_state=2" begin
        d = make_gap_mortality_data()
        spinup_fac = 10.0

        CLM.cn_gap_mortality!(
            d.mask, d.bounds, d.params, d.pftcon, d.dgvs,
            d.patch, d.canopystate, d.cnveg_cs, d.cnveg_cf,
            d.cnveg_ns, d.cnveg_nf;
            dt = 1800.0,
            days_per_year = 365.0,
            use_cndv = false,
            use_matrixcn = false,
            spinup_state = 2,
            npcropmin = 17,
            spinup_factor_deadwood = spinup_fac,
        )

        am = 0.02
        m  = am / (365.0 * CLM.SECSPDAY)

        # Carbon dead pools accelerated
        @test d.cnveg_cf.m_deadstemc_to_litter_patch[1] ≈ 50.0 * m * spinup_fac
        @test d.cnveg_cf.m_deadcrootc_to_litter_patch[1] ≈ 25.0 * m * spinup_fac

        # Nitrogen dead pools accelerated (spinup_state==2 && !use_cndv)
        @test d.cnveg_nf.m_deadstemn_to_litter_patch[1] ≈ 0.5 * m * spinup_fac
        @test d.cnveg_nf.m_deadcrootn_to_litter_patch[1] ≈ 0.25 * m * spinup_fac
    end

    # ================================================================
    @testset "cn_gap_mortality! with use_cndv" begin
        d = make_gap_mortality_data()

        # Set greffic for CNDV
        d.dgvs.greffic_patch .= [0.5, 0.5, 0.5]
        d.dgvs.heatstress_patch .= [0.0, 0.0, 0.0]

        CLM.cn_gap_mortality!(
            d.mask, d.bounds, d.params, d.pftcon, d.dgvs,
            d.patch, d.canopystate, d.cnveg_cs, d.cnveg_cf,
            d.cnveg_ns, d.cnveg_nf;
            dt = 1800.0,
            days_per_year = 365.0,
            use_cndv = true,
            use_matrixcn = false,
            spinup_state = 0,
            npcropmin = 17,
            spinup_factor_deadwood = 1.0,
        )

        # Patch 1: ivt=2, woody=1.0 -> CNDV growth efficiency mortality
        # mort_max = 0.01 (ivt != 8), am = 0.01 / (1 + 0.3 * 0.5)
        am_p1 = 0.01 / (1.0 + 0.3 * 0.5)
        am_p1 = min(1.0, am_p1 + 0.0)  # heatstress = 0
        m_p1  = am_p1 / (365.0 * CLM.SECSPDAY)
        @test d.cnveg_cf.m_leafc_to_litter_patch[1] ≈ 10.0 * m_p1

        # Patch 2: ivt=10, woody=0 -> uses r_mort
        am_p2 = 0.02
        m_p2  = am_p2 / (365.0 * CLM.SECSPDAY)
        @test d.cnveg_cf.m_leafc_to_litter_patch[2] ≈ 5.0 * m_p2

        # Patch 1 woody: nind should decrease
        @test d.dgvs.nind_patch[1] ≈ 100.0 * (1.0 - m_p1)

        # Patch 2 non-woody: nind unchanged (CNDV nind update only for woody)
        @test d.dgvs.nind_patch[2] == 100.0
    end

    # ================================================================
    @testset "cn_gap_mortality! with BDT boreal (ivt=8)" begin
        d = make_gap_mortality_data()
        d.patch.itype[1] = 8  # broadleaf_deciduous_boreal_tree

        CLM.cn_gap_mortality!(
            d.mask, d.bounds, d.params, d.pftcon, d.dgvs,
            d.patch, d.canopystate, d.cnveg_cs, d.cnveg_cf,
            d.cnveg_ns, d.cnveg_nf;
            dt = 1800.0,
            days_per_year = 365.0,
            use_cndv = true,
            use_matrixcn = false,
            spinup_state = 0,
            npcropmin = 17,
            spinup_factor_deadwood = 1.0,
        )

        # ivt=8, woody=1 -> mort_max = 0.03
        am_bdt = 0.03 / (1.0 + 0.3 * 0.5)
        m_bdt  = am_bdt / (365.0 * CLM.SECSPDAY)
        @test d.cnveg_cf.m_leafc_to_litter_patch[1] ≈ 10.0 * m_bdt
    end

    # ================================================================
    @testset "cn_gap_patch_to_column! basic" begin
        d = make_gap_mortality_data()

        # First run gap mortality to populate patch fluxes
        CLM.cn_gap_mortality!(
            d.mask, d.bounds, d.params, d.pftcon, d.dgvs,
            d.patch, d.canopystate, d.cnveg_cs, d.cnveg_cf,
            d.cnveg_ns, d.cnveg_nf;
            dt = 1800.0,
            days_per_year = 365.0,
            use_cndv = false,
            spinup_factor_deadwood = 1.0,
        )

        # Now aggregate to column level
        CLM.cn_gap_patch_to_column!(
            d.mask, d.bounds, d.pftcon, d.patch,
            d.cnveg_cf, d.cnveg_nf,
            d.leaf_prof, d.froot_prof, d.croot_prof, d.stem_prof;
            nlevdecomp = d.nlevdecomp,
            i_litr_min = 1,
            i_litr_max = d.n_litr,
            i_met_lit = 1,
        )

        am = 0.02
        m  = am / (365.0 * CLM.SECSPDAY)

        # Column 1 has patches 1 (wt=0.5) and 2 (wt=0.3)
        # Column 2 has patch 3 (wt=0.2)

        # CWD from column 1:
        # stem: (livestemc+deadstemc)*wt*stem_prof for p1 and p2
        cwdc_col1 = ((20.0 + 50.0) * m * 0.5 * 1.0 +   # patch 1 stem
                      (0.0 + 0.0)  * m * 0.3 * 1.0 +   # patch 2 stem
                      (10.0 + 25.0) * m * 0.5 * 1.0 +   # patch 1 croot
                      (0.0 + 0.0)  * m * 0.3 * 1.0)     # patch 2 croot
        @test d.cnveg_cf.gap_mortality_c_to_cwdc_col[1, 1] ≈ cwdc_col1

        # CWD from column 2:
        cwdc_col2 = ((15.0 + 40.0) * m * 0.2 * 1.0 +   # patch 3 stem
                      (8.0 + 20.0)  * m * 0.2 * 1.0)    # patch 3 croot
        @test d.cnveg_cf.gap_mortality_c_to_cwdc_col[2, 1] ≈ cwdc_col2

        # Litter pools should be non-zero
        @test d.cnveg_cf.gap_mortality_c_to_litr_c_col[1, 1, 1] > 0.0
        @test d.cnveg_nf.gap_mortality_n_to_litr_n_col[1, 1, 1] > 0.0

        # CWD N from column 1:
        cwdn_col1 = ((0.4 + 0.5) * m * 0.5 * 1.0 +     # patch 1 stem N
                      (0.0 + 0.0)  * m * 0.3 * 1.0 +    # patch 2 stem N
                      (0.2 + 0.25) * m * 0.5 * 1.0 +    # patch 1 croot N
                      (0.0 + 0.0)  * m * 0.3 * 1.0)     # patch 2 croot N
        @test d.cnveg_nf.gap_mortality_n_to_cwdn_col[1, 1] ≈ cwdn_col1
    end

    # ================================================================
    @testset "cn_gap_patch_to_column! accumulation" begin
        d = make_gap_mortality_data()

        CLM.cn_gap_mortality!(
            d.mask, d.bounds, d.params, d.pftcon, d.dgvs,
            d.patch, d.canopystate, d.cnveg_cs, d.cnveg_cf,
            d.cnveg_ns, d.cnveg_nf;
            days_per_year = 365.0,
            spinup_factor_deadwood = 1.0,
        )

        # Pre-set a value to verify accumulation (+=)
        d.cnveg_cf.gap_mortality_c_to_cwdc_col[1, 1] = 99.0

        CLM.cn_gap_patch_to_column!(
            d.mask, d.bounds, d.pftcon, d.patch,
            d.cnveg_cf, d.cnveg_nf,
            d.leaf_prof, d.froot_prof, d.croot_prof, d.stem_prof;
            nlevdecomp = d.nlevdecomp,
            i_litr_min = 1,
            i_litr_max = d.n_litr,
            i_met_lit = 1,
        )

        # Should be 99.0 + the computed value (not overwritten)
        @test d.cnveg_cf.gap_mortality_c_to_cwdc_col[1, 1] > 99.0
    end

    # ================================================================
    @testset "Zero mortality rate yields zero fluxes" begin
        d = make_gap_mortality_data()
        d.params.r_mort .= 0.0  # zero mortality

        CLM.cn_gap_mortality!(
            d.mask, d.bounds, d.params, d.pftcon, d.dgvs,
            d.patch, d.canopystate, d.cnveg_cs, d.cnveg_cf,
            d.cnveg_ns, d.cnveg_nf;
            days_per_year = 365.0,
            spinup_factor_deadwood = 1.0,
        )

        @test all(d.cnveg_cf.m_leafc_to_litter_patch .== 0.0)
        @test all(d.cnveg_cf.m_frootc_to_litter_patch .== 0.0)
        @test all(d.cnveg_cf.m_livestemc_to_litter_patch .== 0.0)
        @test all(d.cnveg_cf.m_deadstemc_to_litter_patch .== 0.0)
        @test all(d.cnveg_nf.m_leafn_to_litter_patch .== 0.0)
    end

    # ================================================================
    @testset "Masked patches are skipped" begin
        d = make_gap_mortality_data()
        mask = BitVector([true, false, true])  # skip patch 2

        CLM.cn_gap_mortality!(
            mask, d.bounds, d.params, d.pftcon, d.dgvs,
            d.patch, d.canopystate, d.cnveg_cs, d.cnveg_cf,
            d.cnveg_ns, d.cnveg_nf;
            days_per_year = 365.0,
            spinup_factor_deadwood = 1.0,
        )

        # Patch 2 should remain zero (masked out)
        @test d.cnveg_cf.m_leafc_to_litter_patch[2] == 0.0
        @test d.cnveg_nf.m_leafn_to_litter_patch[2] == 0.0

        # Patches 1 and 3 should be non-zero
        @test d.cnveg_cf.m_leafc_to_litter_patch[1] > 0.0
        @test d.cnveg_cf.m_leafc_to_litter_patch[3] > 0.0
    end

end  # @testset "Gap Mortality"
