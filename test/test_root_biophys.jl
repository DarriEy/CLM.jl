@testset "RootBiophys" begin

    # ======================================================================
    # Test configuration and constants
    # ======================================================================

    @testset "Constants" begin
        @test CLM.ZENG_2001_ROOT    == 0
        @test CLM.JACKSON_1996_ROOT == 1
        @test CLM.KOVEN_EXP_ROOT    == 2
    end

    @testset "RootingProfileConfig defaults" begin
        cfg = CLM.RootingProfileConfig()
        @test cfg.rooting_profile_method_water    == CLM.ZENG_2001_ROOT
        @test cfg.rooting_profile_method_carbon   == CLM.ZENG_2001_ROOT
        @test cfg.rooting_profile_varindex_water  == 1
        @test cfg.rooting_profile_varindex_carbon == 2
    end

    @testset "init_rootprof!" begin
        cfg = CLM.RootingProfileConfig()

        # Set to non-default values
        CLM.init_rootprof!(cfg;
            method_water=CLM.JACKSON_1996_ROOT,
            method_carbon=CLM.KOVEN_EXP_ROOT,
            varindex_water=2,
            varindex_carbon=1)
        @test cfg.rooting_profile_method_water    == CLM.JACKSON_1996_ROOT
        @test cfg.rooting_profile_method_carbon   == CLM.KOVEN_EXP_ROOT
        @test cfg.rooting_profile_varindex_water  == 2
        @test cfg.rooting_profile_varindex_carbon == 1

        # Reset to defaults
        CLM.init_rootprof!(cfg)
        @test cfg.rooting_profile_method_water    == CLM.ZENG_2001_ROOT
        @test cfg.rooting_profile_method_carbon   == CLM.ZENG_2001_ROOT
        @test cfg.rooting_profile_varindex_water  == 1
        @test cfg.rooting_profile_varindex_carbon == 2
    end

    # ======================================================================
    # Helper: set up test data for root profile computation
    # ======================================================================
    # We use a simple soil column with uniform 0.1 m layers.
    # nlevsoi = 5 layers, nlevgrnd = 7 (layers 6-7 are bedrock).
    # col_zi[c, lev] = bottom interface depth of soil layer lev.
    #   Layer 1: 0.0 - 0.1 m => zi[c,1] = 0.1
    #   Layer 2: 0.1 - 0.2 m => zi[c,2] = 0.2
    #   Layer 3: 0.2 - 0.3 m => zi[c,3] = 0.3
    #   Layer 4: 0.3 - 0.4 m => zi[c,4] = 0.4
    #   Layer 5: 0.4 - 0.5 m => zi[c,5] = 0.5
    # col_z[c, lev] = center depth of soil layer lev.
    #   z[c,1] = 0.05, z[c,2] = 0.15, ..., z[c,5] = 0.45
    # col_dz[c, lev] = thickness of soil layer lev = 0.1 for all.

    nlevsoi  = 5
    nlevgrnd = 7
    nc       = 2   # number of columns
    np       = 3   # number of patches

    # Column data
    col_zi = zeros(nc, nlevsoi)
    col_z  = zeros(nc, nlevsoi)
    col_dz = zeros(nc, nlevsoi)
    for c in 1:nc
        for lev in 1:nlevsoi
            col_zi[c, lev] = lev * 0.1         # 0.1, 0.2, 0.3, 0.4, 0.5
            col_z[c, lev]  = (lev - 0.5) * 0.1 # 0.05, 0.15, 0.25, 0.35, 0.45
            col_dz[c, lev] = 0.1
        end
    end
    col_nbedrock = [nlevsoi, nlevsoi]  # bedrock at bottom of soil (no redistribution)

    # Patch data: 3 patches
    # patch 1 -> column 1, PFT type 1 (needleleaf_evergreen_temperate_tree)
    # patch 2 -> column 2, PFT type 4 (broadleaf_evergreen_tropical_tree)
    # patch 3 -> column 1, FATES patch
    patch_column   = [1, 2, 1]
    patch_itype    = [1, 4, 0]   # 0-based Fortran PFT indices
    patch_is_fates = [false, false, true]

    bounds_p = 1:np

    # Set up pftcon with test values
    # Save original state
    saved_pftcon_roota = copy(CLM.pftcon.roota_par)
    saved_pftcon_rootb = copy(CLM.pftcon.rootb_par)
    saved_pftcon_rootprof_beta = copy(CLM.pftcon.rootprof_beta)

    # Ensure pftcon is allocated
    if isempty(CLM.pftcon.roota_par)
        CLM.pftcon_allocate!(CLM.pftcon)
    end

    # Set PFT-specific root parameters (1-based Julia index = Fortran index + 1)
    # PFT 1 (Fortran) -> Julia index 2: needleleaf_evergreen_temperate_tree
    CLM.pftcon.roota_par[2] = 7.0   # roota for PFT 1
    CLM.pftcon.rootb_par[2] = 2.0   # rootb for PFT 1
    # PFT 4 (Fortran) -> Julia index 5: broadleaf_evergreen_tropical_tree
    CLM.pftcon.roota_par[5] = 6.0   # roota for PFT 4
    CLM.pftcon.rootb_par[5] = 1.5   # rootb for PFT 4

    # Jackson beta parameters
    CLM.pftcon.rootprof_beta[2, 1] = 0.95  # PFT 1, variant 1
    CLM.pftcon.rootprof_beta[2, 2] = 0.97  # PFT 1, variant 2
    CLM.pftcon.rootprof_beta[5, 1] = 0.96  # PFT 4, variant 1
    CLM.pftcon.rootprof_beta[5, 2] = 0.98  # PFT 4, variant 2

    # ======================================================================
    # Test zeng2001_rootfr!
    # ======================================================================
    @testset "zeng2001_rootfr!" begin
        rootfr = zeros(np, nlevgrnd)

        CLM.zeng2001_rootfr!(rootfr, col_zi, patch_column, patch_itype,
                             patch_is_fates, CLM.pftcon, bounds_p, nlevsoi)

        # Patch 1 (PFT 1): roota=7.0, rootb=2.0
        a, b = 7.0, 2.0
        # Layer 1: zi_upper=0.0, zi_lower=0.1
        expected_1_1 = 0.5 * (exp(-a * 0.0) + exp(-b * 0.0) -
                              exp(-a * 0.1) - exp(-b * 0.1))
        @test rootfr[1, 1] ≈ expected_1_1 atol=1e-12

        # Layer 2: zi_upper=0.1, zi_lower=0.2
        expected_1_2 = 0.5 * (exp(-a * 0.1) + exp(-b * 0.1) -
                              exp(-a * 0.2) - exp(-b * 0.2))
        @test rootfr[1, 2] ≈ expected_1_2 atol=1e-12

        # Layer 3: zi_upper=0.2, zi_lower=0.3
        expected_1_3 = 0.5 * (exp(-a * 0.2) + exp(-b * 0.2) -
                              exp(-a * 0.3) - exp(-b * 0.3))
        @test rootfr[1, 3] ≈ expected_1_3 atol=1e-12

        # Bottom layer (5): zi_upper=0.4 (zi at layer 4), no lower cutoff
        expected_1_5 = 0.5 * (exp(-a * 0.4) + exp(-b * 0.4))
        @test rootfr[1, nlevsoi] ≈ expected_1_5 atol=1e-12

        # Patch 2 (PFT 4): roota=6.0, rootb=1.5
        a2, b2 = 6.0, 1.5
        expected_2_1 = 0.5 * (exp(-a2 * 0.0) + exp(-b2 * 0.0) -
                              exp(-a2 * 0.1) - exp(-b2 * 0.1))
        @test rootfr[2, 1] ≈ expected_2_1 atol=1e-12

        # Patch 3 (FATES): should be all zeros
        for lev in 1:nlevsoi
            @test rootfr[3, lev] == 0.0
        end

        # Root fractions should sum to approximately 1 for non-FATES patches
        total_1 = sum(rootfr[1, 1:nlevsoi])
        @test total_1 > 0.0
        @test total_1 ≈ 1.0 atol=0.05  # close to 1 but depends on layer depth

        # Root fractions should be non-negative
        for p in 1:2, lev in 1:nlevsoi
            @test rootfr[p, lev] >= 0.0
        end
    end

    # ======================================================================
    # Test jackson1996_rootfr!
    # ======================================================================
    @testset "jackson1996_rootfr!" begin
        rootfr = zeros(np, nlevgrnd)
        varindx = 1  # variant index

        CLM.jackson1996_rootfr!(rootfr, col_zi, patch_column, patch_itype,
                                patch_is_fates, CLM.pftcon, bounds_p, nlevsoi,
                                varindx)

        m_to_cm = 100.0

        # Patch 1 (PFT 1): beta = 0.95 (variant 1)
        beta1 = 0.95
        # Layer 1: zi_upper=0.0, zi_lower=0.1
        expected_j_1_1 = beta1^(0.0 * m_to_cm) - beta1^(0.1 * m_to_cm)
        @test rootfr[1, 1] ≈ expected_j_1_1 atol=1e-12

        # Layer 2: zi_upper=0.1, zi_lower=0.2
        expected_j_1_2 = beta1^(0.1 * m_to_cm) - beta1^(0.2 * m_to_cm)
        @test rootfr[1, 2] ≈ expected_j_1_2 atol=1e-12

        # Patch 2 (PFT 4): beta = 0.96 (variant 1)
        beta2 = 0.96
        expected_j_2_1 = beta2^(0.0 * m_to_cm) - beta2^(0.1 * m_to_cm)
        @test rootfr[2, 1] ≈ expected_j_2_1 atol=1e-12

        # FATES patch should be zero
        for lev in 1:nlevsoi
            @test rootfr[3, lev] == 0.0
        end

        # Test with variant 2
        rootfr2 = zeros(np, nlevgrnd)
        CLM.jackson1996_rootfr!(rootfr2, col_zi, patch_column, patch_itype,
                                patch_is_fates, CLM.pftcon, bounds_p, nlevsoi,
                                2)
        beta1_v2 = 0.97
        expected_v2 = beta1_v2^(0.0 * m_to_cm) - beta1_v2^(0.1 * m_to_cm)
        @test rootfr2[1, 1] ≈ expected_v2 atol=1e-12

        # Root fractions should be non-negative and sum should be positive
        for p in 1:2, lev in 1:nlevsoi
            @test rootfr[p, lev] >= 0.0
        end
        @test sum(rootfr[1, 1:nlevsoi]) > 0.0
    end

    # ======================================================================
    # Test exponential_rootfr!
    # ======================================================================
    @testset "exponential_rootfr!" begin
        rootfr = zeros(np, nlevgrnd)

        CLM.exponential_rootfr!(rootfr, col_z, col_dz, patch_column, patch_itype,
                                patch_is_fates, bounds_p, nlevsoi)

        rootprof_exp = 3.0

        # Compute expected values for patch 1
        c = patch_column[1]
        raw_1 = zeros(nlevsoi)
        for lev in 1:nlevsoi
            raw_1[lev] = exp(-rootprof_exp * col_z[c, lev]) * col_dz[c, lev]
        end
        norm_1 = -1.0 / rootprof_exp * (exp(-rootprof_exp * col_z[c, nlevsoi]) - 1.0)
        expected_1 = raw_1 ./ norm_1

        for lev in 1:nlevsoi
            @test rootfr[1, lev] ≈ expected_1[lev] atol=1e-12
        end

        # Root fractions should sum to approximately 1 for non-FATES patches
        total = sum(rootfr[1, 1:nlevsoi])
        @test total > 0.9
        @test total < 1.1

        # Root fractions should be non-negative
        for lev in 1:nlevsoi
            @test rootfr[1, lev] >= 0.0
        end

        # FATES patch: all zeros after normalization (0/norm = 0)
        for lev in 1:nlevsoi
            @test rootfr[3, lev] ≈ 0.0 atol=1e-12
        end

        # Root fractions should decrease with depth (exponential decay)
        for lev in 2:nlevsoi
            @test rootfr[1, lev] < rootfr[1, lev - 1]
        end
    end

    # ======================================================================
    # Test init_vegrootfr! - integration test
    # ======================================================================
    @testset "init_vegrootfr! with Zeng 2001" begin
        rootfr = zeros(np, nlevgrnd)
        cfg = CLM.RootingProfileConfig(
            rooting_profile_method_water=CLM.ZENG_2001_ROOT,
            rooting_profile_method_carbon=CLM.ZENG_2001_ROOT,
            rooting_profile_varindex_water=1,
            rooting_profile_varindex_carbon=2)

        CLM.init_vegrootfr!(rootfr, col_zi, col_z, col_dz, col_nbedrock,
                            patch_column, patch_itype, patch_is_fates,
                            CLM.pftcon, cfg, bounds_p, nlevsoi, nlevgrnd, "water")

        # Layers beyond nlevsoi should be zero
        for p in bounds_p, lev in (nlevsoi + 1):nlevgrnd
            @test rootfr[p, lev] == 0.0
        end

        # Non-FATES patches should have non-zero root fractions
        @test sum(rootfr[1, 1:nlevsoi]) > 0.0
        @test sum(rootfr[2, 1:nlevsoi]) > 0.0

        # FATES patch should be all zeros
        @test sum(rootfr[3, :]) == 0.0
    end

    @testset "init_vegrootfr! with Jackson 1996" begin
        rootfr = zeros(np, nlevgrnd)
        cfg = CLM.RootingProfileConfig(
            rooting_profile_method_water=CLM.JACKSON_1996_ROOT,
            rooting_profile_method_carbon=CLM.JACKSON_1996_ROOT,
            rooting_profile_varindex_water=1,
            rooting_profile_varindex_carbon=2)

        CLM.init_vegrootfr!(rootfr, col_zi, col_z, col_dz, col_nbedrock,
                            patch_column, patch_itype, patch_is_fates,
                            CLM.pftcon, cfg, bounds_p, nlevsoi, nlevgrnd, "water")

        # Non-FATES patches should have non-zero root fractions
        @test sum(rootfr[1, 1:nlevsoi]) > 0.0

        # Layers beyond nlevsoi should be zero
        for p in bounds_p, lev in (nlevsoi + 1):nlevgrnd
            @test rootfr[p, lev] == 0.0
        end
    end

    @testset "init_vegrootfr! with Koven exponential" begin
        rootfr = zeros(np, nlevgrnd)
        cfg = CLM.RootingProfileConfig(
            rooting_profile_method_water=CLM.KOVEN_EXP_ROOT,
            rooting_profile_method_carbon=CLM.KOVEN_EXP_ROOT,
            rooting_profile_varindex_water=1,
            rooting_profile_varindex_carbon=2)

        CLM.init_vegrootfr!(rootfr, col_zi, col_z, col_dz, col_nbedrock,
                            patch_column, patch_itype, patch_is_fates,
                            CLM.pftcon, cfg, bounds_p, nlevsoi, nlevgrnd, "carbon")

        # Non-FATES patches should have non-zero root fractions
        @test sum(rootfr[1, 1:nlevsoi]) > 0.0

        # Layers beyond nlevsoi should be zero
        for lev in (nlevsoi + 1):nlevgrnd
            @test rootfr[1, lev] == 0.0
        end
    end

    @testset "init_vegrootfr! with water_carbon argument" begin
        cfg = CLM.RootingProfileConfig(
            rooting_profile_method_water=CLM.ZENG_2001_ROOT,
            rooting_profile_method_carbon=CLM.JACKSON_1996_ROOT,
            rooting_profile_varindex_water=1,
            rooting_profile_varindex_carbon=1)

        # "water" should use Zeng 2001
        rootfr_water = zeros(np, nlevgrnd)
        CLM.init_vegrootfr!(rootfr_water, col_zi, col_z, col_dz, col_nbedrock,
                            patch_column, patch_itype, patch_is_fates,
                            CLM.pftcon, cfg, bounds_p, nlevsoi, nlevgrnd, "water")

        # "carbon" should use Jackson 1996
        rootfr_carbon = zeros(np, nlevgrnd)
        CLM.init_vegrootfr!(rootfr_carbon, col_zi, col_z, col_dz, col_nbedrock,
                            patch_column, patch_itype, patch_is_fates,
                            CLM.pftcon, cfg, bounds_p, nlevsoi, nlevgrnd, "carbon")

        # The two should differ since they use different methods
        @test rootfr_water[1, 1] != rootfr_carbon[1, 1]

        # Invalid water_carbon should error
        rootfr_bad = zeros(np, nlevgrnd)
        @test_throws ErrorException CLM.init_vegrootfr!(
            rootfr_bad, col_zi, col_z, col_dz, col_nbedrock,
            patch_column, patch_itype, patch_is_fates,
            CLM.pftcon, cfg, bounds_p, nlevsoi, nlevgrnd, "invalid")
    end

    @testset "init_vegrootfr! invalid method" begin
        cfg = CLM.RootingProfileConfig(
            rooting_profile_method_water=99,
            rooting_profile_method_carbon=CLM.ZENG_2001_ROOT,
            rooting_profile_varindex_water=1,
            rooting_profile_varindex_carbon=2)

        rootfr = zeros(np, nlevgrnd)
        @test_throws ErrorException CLM.init_vegrootfr!(
            rootfr, col_zi, col_z, col_dz, col_nbedrock,
            patch_column, patch_itype, patch_is_fates,
            CLM.pftcon, cfg, bounds_p, nlevsoi, nlevgrnd, "water")
    end

    # ======================================================================
    # Test bedrock redistribution
    # ======================================================================
    @testset "init_vegrootfr! bedrock redistribution" begin
        # Shallow bedrock: only 2 layers of soil, rest is bedrock
        col_nbedrock_shallow = [2, 2]
        rootfr = zeros(np, nlevgrnd)
        cfg = CLM.RootingProfileConfig(
            rooting_profile_method_water=CLM.ZENG_2001_ROOT,
            rooting_profile_method_carbon=CLM.ZENG_2001_ROOT,
            rooting_profile_varindex_water=1,
            rooting_profile_varindex_carbon=2)

        CLM.init_vegrootfr!(rootfr, col_zi, col_z, col_dz, col_nbedrock_shallow,
                            patch_column, patch_itype, patch_is_fates,
                            CLM.pftcon, cfg, bounds_p, nlevsoi, nlevgrnd, "water")

        # For non-FATES patches:
        # - Layers 3:nlevsoi should be zero (below bedrock)
        for lev in 3:nlevsoi
            @test rootfr[1, lev] == 0.0
            @test rootfr[2, lev] == 0.0
        end

        # - Layers nlevsoi+1:nlevgrnd should be zero (bedrock in extended layers)
        for lev in (nlevsoi + 1):nlevgrnd
            @test rootfr[1, lev] == 0.0
        end

        # - Layers 1:2 should have the original roots plus redistributed roots
        # The total should equal what was originally in layers 1:nlevsoi
        # Compute reference (no bedrock redistribution)
        rootfr_ref = zeros(np, nlevgrnd)
        CLM.init_vegrootfr!(rootfr_ref, col_zi, col_z, col_dz,
                            fill(nlevsoi, nc),  # no bedrock limitation
                            patch_column, patch_itype, patch_is_fates,
                            CLM.pftcon, cfg, bounds_p, nlevsoi, nlevgrnd, "water")

        # Total root fraction in top 2 layers (with bedrock) should equal
        # total root fraction in all nlevsoi layers (without bedrock)
        for p in 1:2  # non-FATES patches
            total_with_bedrock = sum(rootfr[p, 1:2])
            total_without_bedrock = sum(rootfr_ref[p, 1:nlevsoi])
            @test total_with_bedrock ≈ total_without_bedrock atol=1e-12
        end

        # Root fractions in layers 1 and 2 should be larger than the reference
        # (because they include redistributed roots from below)
        for p in 1:2
            @test rootfr[p, 1] > rootfr_ref[p, 1]
            @test rootfr[p, 2] > rootfr_ref[p, 2]
        end
    end

    @testset "init_vegrootfr! bedrock at layer 1" begin
        # Very shallow bedrock: only 1 layer
        col_nbedrock_1 = [1, 1]
        rootfr = zeros(np, nlevgrnd)
        cfg = CLM.RootingProfileConfig(
            rooting_profile_method_water=CLM.ZENG_2001_ROOT,
            rooting_profile_method_carbon=CLM.ZENG_2001_ROOT,
            rooting_profile_varindex_water=1,
            rooting_profile_varindex_carbon=2)

        CLM.init_vegrootfr!(rootfr, col_zi, col_z, col_dz, col_nbedrock_1,
                            patch_column, patch_itype, patch_is_fates,
                            CLM.pftcon, cfg, bounds_p, nlevsoi, nlevgrnd, "water")

        # All roots should be concentrated in layer 1
        for p in 1:2
            @test rootfr[p, 1] > 0.0
            for lev in 2:nlevgrnd
                @test rootfr[p, lev] == 0.0
            end
        end
    end

    # ======================================================================
    # Test single-patch scenarios for exact numerical values
    # ======================================================================
    @testset "zeng2001 exact values single patch" begin
        # Single column, single patch, 3 layers
        nlev = 3
        zi_test = zeros(1, nlev)
        zi_test[1, 1] = 0.05   # 5 cm
        zi_test[1, 2] = 0.20   # 20 cm
        zi_test[1, 3] = 1.00   # 1 m
        rootfr = zeros(1, nlev)

        a, b = 5.0, 3.0
        pft_test = CLM.PftconType()
        CLM.pftcon_allocate!(pft_test)
        pft_test.roota_par[2] = a  # PFT 1 (Julia index 2)
        pft_test.rootb_par[2] = b

        CLM.zeng2001_rootfr!(rootfr, zi_test, [1], [1], [false],
                             pft_test, 1:1, nlev)

        # Layer 1: zi_upper=0.0, zi_lower=0.05
        exp_1 = 0.5 * (1.0 + 1.0 - exp(-a * 0.05) - exp(-b * 0.05))
        @test rootfr[1, 1] ≈ exp_1 atol=1e-14

        # Layer 2: zi_upper=0.05, zi_lower=0.20
        exp_2 = 0.5 * (exp(-a * 0.05) + exp(-b * 0.05) -
                       exp(-a * 0.20) - exp(-b * 0.20))
        @test rootfr[1, 2] ≈ exp_2 atol=1e-14

        # Layer 3 (bottom): zi_upper=0.20
        exp_3 = 0.5 * (exp(-a * 0.20) + exp(-b * 0.20))
        @test rootfr[1, 3] ≈ exp_3 atol=1e-14

        # Should sum to 1.0 exactly: 0.5*(1+1) = 1.0 after telescoping
        @test sum(rootfr[1, :]) ≈ 1.0 atol=1e-14
    end

    @testset "jackson1996 exact values single patch" begin
        nlev = 3
        zi_test = zeros(1, nlev)
        zi_test[1, 1] = 0.10   # 10 cm
        zi_test[1, 2] = 0.30   # 30 cm
        zi_test[1, 3] = 1.00   # 100 cm
        rootfr = zeros(1, nlev)

        beta = 0.97
        pft_test = CLM.PftconType()
        CLM.pftcon_allocate!(pft_test)
        pft_test.rootprof_beta[2, 1] = beta  # PFT 1, variant 1

        CLM.jackson1996_rootfr!(rootfr, zi_test, [1], [1], [false],
                                pft_test, 1:1, nlev, 1)

        m_to_cm = 100.0
        # Layer 1: beta^(0*100) - beta^(10) = 1 - 0.97^10
        exp_1 = beta^(0.0 * m_to_cm) - beta^(0.10 * m_to_cm)
        @test rootfr[1, 1] ≈ exp_1 atol=1e-14

        # Layer 2: beta^(10) - beta^(30)
        exp_2 = beta^(0.10 * m_to_cm) - beta^(0.30 * m_to_cm)
        @test rootfr[1, 2] ≈ exp_2 atol=1e-14

        # Layer 3: beta^(30) - beta^(100)
        exp_3 = beta^(0.30 * m_to_cm) - beta^(1.00 * m_to_cm)
        @test rootfr[1, 3] ≈ exp_3 atol=1e-14

        # Total should be 1 - beta^100
        @test sum(rootfr[1, :]) ≈ (1.0 - beta^100.0) atol=1e-12

        # All fractions should be non-negative
        for lev in 1:nlev
            @test rootfr[1, lev] >= 0.0
        end
    end

    @testset "exponential exact values single patch" begin
        nlev = 4
        z_test  = zeros(1, nlev)
        dz_test = zeros(1, nlev)
        # Non-uniform layers
        z_test[1, 1]  = 0.025;  dz_test[1, 1] = 0.05
        z_test[1, 2]  = 0.10;   dz_test[1, 2] = 0.10
        z_test[1, 3]  = 0.25;   dz_test[1, 3] = 0.20
        z_test[1, 4]  = 0.50;   dz_test[1, 4] = 0.30

        rootfr = zeros(1, nlev)
        rootprof_exp = 3.0

        CLM.exponential_rootfr!(rootfr, z_test, dz_test, [1], [1], [false],
                                1:1, nlev)

        # Compute expected
        raw = [exp(-rootprof_exp * z_test[1, l]) * dz_test[1, l] for l in 1:nlev]
        norm = -1.0 / rootprof_exp * (exp(-rootprof_exp * z_test[1, nlev]) - 1.0)
        expected = raw ./ norm

        for lev in 1:nlev
            @test rootfr[1, lev] ≈ expected[lev] atol=1e-14
        end

        # All fractions should be non-negative
        for lev in 1:nlev
            @test rootfr[1, lev] >= 0.0
        end
    end

    # ======================================================================
    # Test FATES patches produce zero for all methods
    # ======================================================================
    @testset "FATES patches produce zero root fractions" begin
        rootfr = zeros(1, nlevsoi)

        # All-FATES patch
        CLM.zeng2001_rootfr!(rootfr, col_zi, [1], [0], [true],
                             CLM.pftcon, 1:1, nlevsoi)
        @test all(rootfr .== 0.0)

        rootfr .= 999.0
        CLM.jackson1996_rootfr!(rootfr, col_zi, [1], [0], [true],
                                CLM.pftcon, 1:1, nlevsoi, 1)
        @test all(rootfr .== 0.0)

        rootfr .= 999.0
        CLM.exponential_rootfr!(rootfr, col_z, col_dz, [1], [0], [true],
                                1:1, nlevsoi)
        # For FATES, rootfr = 0/norm = 0
        @test all(abs.(rootfr) .< 1e-14)
    end

    # ======================================================================
    # Test with multiple patches sharing same column
    # ======================================================================
    @testset "Multiple patches on same column" begin
        # 2 patches on same column but different PFTs
        rootfr_multi = zeros(2, nlevgrnd)
        cfg = CLM.RootingProfileConfig(
            rooting_profile_method_water=CLM.ZENG_2001_ROOT,
            rooting_profile_method_carbon=CLM.ZENG_2001_ROOT,
            rooting_profile_varindex_water=1,
            rooting_profile_varindex_carbon=2)

        CLM.init_vegrootfr!(rootfr_multi, col_zi, col_z, col_dz, col_nbedrock,
                            [1, 1], [1, 4], [false, false],
                            CLM.pftcon, cfg, 1:2, nlevsoi, nlevgrnd, "water")

        # Both patches use the same column geometry but different PFT parameters
        # So they should have different root profiles
        @test rootfr_multi[1, 1] != rootfr_multi[2, 1]
        @test sum(rootfr_multi[1, 1:nlevsoi]) > 0.0
        @test sum(rootfr_multi[2, 1:nlevsoi]) > 0.0
    end

    # ======================================================================
    # Test global config instance
    # ======================================================================
    @testset "Global rooting_profile_config" begin
        saved = CLM.rooting_profile_config.rooting_profile_method_water

        CLM.init_rootprof!(CLM.rooting_profile_config;
            method_water=CLM.JACKSON_1996_ROOT)
        @test CLM.rooting_profile_config.rooting_profile_method_water == CLM.JACKSON_1996_ROOT

        # Restore
        CLM.init_rootprof!(CLM.rooting_profile_config;
            method_water=saved)
    end

    # ======================================================================
    # Restore pftcon state
    # ======================================================================
    if !isempty(saved_pftcon_roota)
        CLM.pftcon.roota_par .= saved_pftcon_roota
        CLM.pftcon.rootb_par .= saved_pftcon_rootb
        CLM.pftcon.rootprof_beta .= saved_pftcon_rootprof_beta
    end

end
