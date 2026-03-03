@testset "C Isotope Flux (CNCIsoFluxMod)" begin

    # -----------------------------------------------------------------------
    # Test c_iso_flux_calc! — 1D core calculation
    # -----------------------------------------------------------------------
    @testset "c_iso_flux_calc! basic 1D" begin
        n = 5
        mask = trues(n)
        bounds = 1:n

        # Case 1: normal non-zero states
        ctot_flux  = [1.0, 2.0, 3.0, 4.0, 5.0]
        ctot_state = [10.0, 20.0, 30.0, 40.0, 50.0]
        ciso_state = [1.0, 2.0, 3.0, 4.0, 5.0]
        ciso_flux  = zeros(n)

        CLM.c_iso_flux_calc!(ciso_flux, ctot_flux, ciso_state, ctot_state,
                             mask, bounds, 1.0, "c13")

        for i in 1:n
            @test ciso_flux[i] ≈ ctot_flux[i] * (ciso_state[i] / ctot_state[i]) * 1.0
        end

        # Case 2: zero total state → zero iso flux
        ctot_state2 = [0.0, 20.0, 0.0, 40.0, 0.0]
        ciso_state2 = [0.0, 2.0, 0.0, 4.0, 0.0]
        ciso_flux2  = fill(999.0, n)
        CLM.c_iso_flux_calc!(ciso_flux2, ctot_flux, ciso_state2, ctot_state2,
                             mask, bounds, 1.0, "c13")
        @test ciso_flux2[1] == 0.0
        @test ciso_flux2[3] == 0.0
        @test ciso_flux2[5] == 0.0
        @test ciso_flux2[2] ≈ 2.0 * (2.0 / 20.0)
        @test ciso_flux2[4] ≈ 4.0 * (4.0 / 40.0)

        # Case 3: mask filtering
        mask3 = BitVector([true, false, true, false, true])
        ciso_flux3 = fill(-1.0, n)
        CLM.c_iso_flux_calc!(ciso_flux3, ctot_flux, ciso_state, ctot_state,
                             mask3, bounds, 1.0, "c13")
        @test ciso_flux3[1] ≈ 1.0 * (1.0 / 10.0)
        @test ciso_flux3[2] == -1.0  # unchanged (masked out)
        @test ciso_flux3[3] ≈ 3.0 * (3.0 / 30.0)
        @test ciso_flux3[4] == -1.0  # unchanged (masked out)
    end

    # -----------------------------------------------------------------------
    # Test c_iso_flux_calc! — C14 fractionation doubling
    # -----------------------------------------------------------------------
    @testset "c_iso_flux_calc! C14 fractionation" begin
        n = 3
        mask = trues(n)
        bounds = 1:n
        ctot_flux  = [10.0, 20.0, 30.0]
        ctot_state = [100.0, 200.0, 300.0]
        ciso_state = [10.0, 20.0, 30.0]

        # C13 with frax = 0.95
        ciso_flux_c13 = zeros(n)
        CLM.c_iso_flux_calc!(ciso_flux_c13, ctot_flux, ciso_state, ctot_state,
                             mask, bounds, 0.95, "c13")

        # C14 with same frax_c13 = 0.95 → frax = 1 + (1-0.95)*2 = 1.10
        ciso_flux_c14 = zeros(n)
        CLM.c_iso_flux_calc!(ciso_flux_c14, ctot_flux, ciso_state, ctot_state,
                             mask, bounds, 0.95, "c14")

        for i in 1:n
            ratio = ciso_state[i] / ctot_state[i]
            @test ciso_flux_c13[i] ≈ ctot_flux[i] * ratio * 0.95
            @test ciso_flux_c14[i] ≈ ctot_flux[i] * ratio * 1.10
        end
    end

    # -----------------------------------------------------------------------
    # Test c_iso_flux_calc! — invalid isotope error
    # -----------------------------------------------------------------------
    @testset "c_iso_flux_calc! invalid isotope" begin
        n = 2
        mask = trues(n)
        @test_throws ErrorException CLM.c_iso_flux_calc!(
            zeros(n), ones(n), ones(n), ones(n),
            mask, 1:n, 1.0, "c15")
    end

    # -----------------------------------------------------------------------
    # Test c_iso_flux_calc_2d_flux!
    # -----------------------------------------------------------------------
    @testset "c_iso_flux_calc_2d_flux!" begin
        n = 4
        nk = 2
        mask = trues(n)
        bounds = 1:n
        ctot_flux  = [1.0 3.0; 2.0 4.0; 3.0 5.0; 4.0 6.0]
        ciso_state = [10.0, 20.0, 30.0, 40.0]
        ctot_state = [100.0, 200.0, 300.0, 400.0]
        ciso_flux  = zeros(n, nk)

        CLM.c_iso_flux_calc_2d_flux!(ciso_flux, ctot_flux, ciso_state, ctot_state,
                                     mask, bounds, 1.0, "c13")

        for k in 1:nk
            for i in 1:n
                @test ciso_flux[i, k] ≈ ctot_flux[i, k] * (ciso_state[i] / ctot_state[i])
            end
        end
    end

    # -----------------------------------------------------------------------
    # Test c_iso_flux_calc_2d_both!
    # -----------------------------------------------------------------------
    @testset "c_iso_flux_calc_2d_both!" begin
        n = 3
        nk = 2
        mask = trues(n)
        bounds = 1:n
        ctot_flux  = [1.0 4.0; 2.0 5.0; 3.0 6.0]
        ctot_state = [10.0 40.0; 20.0 50.0; 30.0 60.0]
        ciso_state = [1.0 4.0; 2.0 5.0; 3.0 6.0]
        ciso_flux  = zeros(n, nk)

        CLM.c_iso_flux_calc_2d_both!(ciso_flux, ctot_flux, ciso_state, ctot_state,
                                      mask, bounds, 1.0, "c13")

        for k in 1:nk
            for i in 1:n
                @test ciso_flux[i, k] ≈ ctot_flux[i, k] * (ciso_state[i, k] / ctot_state[i, k])
            end
        end
    end

    # -----------------------------------------------------------------------
    # Test cn_c_iso_gap_pft_to_column! — basic column aggregation
    # -----------------------------------------------------------------------
    @testset "cn_c_iso_gap_pft_to_column! basic" begin
        np = 2
        nc = 1
        nlevdecomp = 1
        i_litr_min = 1
        i_litr_max = 2
        i_met_lit = 1

        # Create minimal data structures
        cf = CLM.CNVegCarbonFluxData()
        # Allocate needed fields
        cf.m_leafc_to_litter_patch = [0.1, 0.2]
        cf.m_frootc_to_litter_patch = [0.05, 0.1]
        cf.m_livestemc_to_litter_patch = [0.01, 0.02]
        cf.m_deadstemc_to_litter_patch = [0.01, 0.02]
        cf.m_livecrootc_to_litter_patch = [0.01, 0.02]
        cf.m_deadcrootc_to_litter_patch = [0.01, 0.02]
        cf.m_leafc_storage_to_litter_patch = [0.005, 0.01]
        cf.m_frootc_storage_to_litter_patch = [0.005, 0.01]
        cf.m_livestemc_storage_to_litter_patch = [0.005, 0.01]
        cf.m_deadstemc_storage_to_litter_patch = [0.005, 0.01]
        cf.m_livecrootc_storage_to_litter_patch = [0.005, 0.01]
        cf.m_deadcrootc_storage_to_litter_patch = [0.005, 0.01]
        cf.m_gresp_storage_to_litter_patch = [0.001, 0.002]
        cf.m_leafc_xfer_to_litter_patch = [0.001, 0.002]
        cf.m_frootc_xfer_to_litter_patch = [0.001, 0.002]
        cf.m_livestemc_xfer_to_litter_patch = [0.001, 0.002]
        cf.m_deadstemc_xfer_to_litter_patch = [0.001, 0.002]
        cf.m_livecrootc_xfer_to_litter_patch = [0.001, 0.002]
        cf.m_deadcrootc_xfer_to_litter_patch = [0.001, 0.002]
        cf.m_gresp_xfer_to_litter_patch = [0.001, 0.002]

        cf.gap_mortality_c_to_litr_c_col = zeros(nc, nlevdecomp, i_litr_max)
        cf.gap_mortality_c_to_cwdc_col = zeros(nc, nlevdecomp)

        sb_state = CLM.SoilBiogeochemStateData()
        sb_state.leaf_prof_patch = ones(np, nlevdecomp)
        sb_state.froot_prof_patch = ones(np, nlevdecomp)
        sb_state.croot_prof_patch = ones(np, nlevdecomp)
        sb_state.stem_prof_patch = ones(np, nlevdecomp)

        mask = trues(np)
        bounds = 1:np

        lf_f = ones(20, i_litr_max)  # all litter fractions = 1
        fr_f = ones(20, i_litr_max)
        patch_column = [1, 1]
        patch_itype = [1, 1]
        patch_wtcol = [0.5, 0.5]

        CLM.cn_c_iso_gap_pft_to_column!(cf, sb_state, mask, bounds, nlevdecomp,
                                        patch_column=patch_column,
                                        patch_itype=patch_itype,
                                        patch_wtcol=patch_wtcol,
                                        lf_f=lf_f, fr_f=fr_f,
                                        i_litr_min=i_litr_min,
                                        i_litr_max=i_litr_max,
                                        i_met_lit=i_met_lit)

        # gap_mortality_c_to_cwdc_col should be non-zero (stems + roots)
        @test cf.gap_mortality_c_to_cwdc_col[1, 1] > 0.0

        # gap_mortality_c_to_litr_c_col should be non-zero
        @test cf.gap_mortality_c_to_litr_c_col[1, 1, 1] > 0.0  # metabolic litter
        @test cf.gap_mortality_c_to_litr_c_col[1, 1, 2] > 0.0  # other litter
    end

    # -----------------------------------------------------------------------
    # Test c_iso_flux_calc! — zero iso state with non-zero total state
    # -----------------------------------------------------------------------
    @testset "c_iso_flux_calc! zero iso state" begin
        n = 3
        mask = trues(n)
        bounds = 1:n
        ctot_flux  = [5.0, 10.0, 15.0]
        ctot_state = [50.0, 100.0, 150.0]
        ciso_state = [0.0, 0.0, 0.0]   # zero isotope state
        ciso_flux  = fill(999.0, n)

        CLM.c_iso_flux_calc!(ciso_flux, ctot_flux, ciso_state, ctot_state,
                             mask, bounds, 1.0, "c13")

        # zero iso state → zero iso flux
        for i in 1:n
            @test ciso_flux[i] == 0.0
        end
    end

    # -----------------------------------------------------------------------
    # Test c_iso_flux_calc! — ratio preservation
    # -----------------------------------------------------------------------
    @testset "c_iso_flux_calc! ratio preservation" begin
        n = 4
        mask = trues(n)
        bounds = 1:n
        ratio = 0.011  # ~1.1% isotope ratio
        ctot_flux  = [100.0, 200.0, 50.0, 75.0]
        ctot_state = [1000.0, 2000.0, 500.0, 750.0]
        ciso_state = ctot_state .* ratio
        ciso_flux  = zeros(n)

        CLM.c_iso_flux_calc!(ciso_flux, ctot_flux, ciso_state, ctot_state,
                             mask, bounds, 1.0, "c13")

        # iso_flux / tot_flux should equal the isotope ratio
        for i in 1:n
            @test ciso_flux[i] / ctot_flux[i] ≈ ratio atol=1e-15
        end
    end

end
