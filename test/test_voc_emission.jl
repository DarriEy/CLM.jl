@testset "VOC Emission (VOCEmissionMod)" begin

    # -----------------------------------------------------------------------
    # Test data types
    # -----------------------------------------------------------------------
    @testset "VOCEmisData construction and init/clean" begin
        voc = CLM.VOCEmisData()
        @test length(voc.Eopt_out_patch) == 0
        @test size(voc.vocflx_patch) == (0, 0)

        np = 5; ng = 2; n_megcomps = 3; n_mechcomps = 2
        CLM.vocemis_init!(voc, np, ng, n_megcomps, n_mechcomps)
        @test length(voc.Eopt_out_patch) == np
        @test length(voc.topt_out_patch) == np
        @test length(voc.alpha_out_patch) == np
        @test length(voc.cp_out_patch) == np
        @test length(voc.gamma_out_patch) == np
        @test length(voc.gammaL_out_patch) == np
        @test length(voc.gammaT_out_patch) == np
        @test length(voc.gammaP_out_patch) == np
        @test length(voc.gammaA_out_patch) == np
        @test length(voc.gammaS_out_patch) == np
        @test length(voc.gammaC_out_patch) == np
        @test length(voc.vocflx_tot_patch) == np
        @test size(voc.vocflx_patch) == (np, n_mechcomps)
        @test size(voc.efisop_grc) == (6, ng)
        @test size(voc.meg_flux_out) == (n_megcomps, np)
        @test all(isnan, voc.Eopt_out_patch)
        @test all(iszero, voc.meg_flux_out)

        CLM.vocemis_clean!(voc)
        @test length(voc.Eopt_out_patch) == 0
        @test size(voc.vocflx_patch) == (0, 0)
        @test size(voc.efisop_grc) == (0, 0)
    end

    @testset "MEGANFactors construction and init" begin
        mf = CLM.MEGANFactors()
        @test length(mf.Agro) == 0

        CLM.megan_factors_init!(mf, 20)
        @test length(mf.Agro) == 20
        @test length(mf.Amat) == 20
        @test length(mf.Anew) == 20
        @test length(mf.Aold) == 20
        @test length(mf.betaT) == 20
        @test length(mf.ct1) == 20
        @test length(mf.ct2) == 20
        @test length(mf.Ceo) == 20
        @test length(mf.LDF) == 20
        @test all(==(1.0), mf.Agro)
        @test all(==(1.0), mf.Amat)
        @test mf.LDF[1] ≈ 0.6
    end

    @testset "MEGANCompound construction" begin
        mc = CLM.MEGANCompound(
            name = "isoprene",
            index = 1,
            class_number = 1,
            molec_weight = 68.12,
            coeff = 1.0,
            emis_factors = fill(600.0, 15)
        )
        @test mc.name == "isoprene"
        @test mc.molec_weight ≈ 68.12
        @test length(mc.emis_factors) == 15
    end

    @testset "MEGANMechComp construction" begin
        mmc = CLM.MEGANMechComp(
            name = "ISOP",
            n_megan_comps = 1,
            megan_indices = [1]
        )
        @test mmc.name == "ISOP"
        @test mmc.n_megan_comps == 1
        @test mmc.megan_indices == [1]
    end

    # -----------------------------------------------------------------------
    # Test helper functions
    # -----------------------------------------------------------------------
    @testset "get_gamma_L" begin
        # With valid fsun240 -> use cce = 0.30
        @test CLM.get_gamma_L(0.5, 2.0) ≈ 0.30 * 2.0
        @test CLM.get_gamma_L(0.5, 0.0) ≈ 0.0

        # With fsun240 = 0 -> use cce1 = 0.24
        @test CLM.get_gamma_L(0.0, 2.0) ≈ 0.24 * 2.0

        # With fsun240 very large (> 1e30) -> use cce1 = 0.24
        @test CLM.get_gamma_L(2.0e30, 2.0) ≈ 0.24 * 2.0
    end

    @testset "get_gamma_SM" begin
        # btran >= 1 -> gamma_SM = 1
        @test CLM.get_gamma_SM(1.0) ≈ 1.0
        @test CLM.get_gamma_SM(1.5) ≈ 1.0

        # btran = 0 -> compute value
        gam_0 = CLM.get_gamma_SM(0.0)
        @test gam_0 > 0.0
        @test gam_0 < 1.0

        # Monotonicity: gamma_SM should increase with btran
        @test CLM.get_gamma_SM(0.3) > CLM.get_gamma_SM(0.1)
        @test CLM.get_gamma_SM(0.8) > CLM.get_gamma_SM(0.5)

        # btran = 0.2 (at threshold)
        gam_thresh = CLM.get_gamma_SM(0.2)
        # At threshold: 1 / (1 + b1*exp(0)) = 1 / (1 + 3.2552)
        @test gam_thresh ≈ 1.0 / (1.0 + 3.2552) atol=1e-4
    end

    @testset "get_gamma_P" begin
        # Test with fixed alpha and cp (fsun240 conditions not met)
        (gamma_p, cp, alpha) = CLM.get_gamma_P(
            200.0, 200.0, 200.0,  # par sun
            50.0, 50.0, 50.0,     # par shade
            0.5, 0.0,             # fsun=0.5, fsun240=0 (triggers fixed)
            100.0, 100.0,         # forc_solad240, forc_solai240
            0.6                   # LDF
        )
        @test gamma_p > 0.0
        @test alpha ≈ 0.001  # alpha_fix
        @test cp ≈ 1.21       # cp_fix

        # With LDF=0 -> gamma_p should be 1.0 (no light dependence)
        (gamma_p0, _, _) = CLM.get_gamma_P(
            200.0, 200.0, 200.0,
            50.0, 50.0, 50.0,
            0.5, 0.0, 100.0, 100.0, 0.0)
        @test gamma_p0 ≈ 1.0

        # With LDF=1 -> gamma_p = gamma_p_LDF
        (gamma_p1, _, _) = CLM.get_gamma_P(
            200.0, 200.0, 200.0,
            50.0, 50.0, 50.0,
            0.5, 0.0, 100.0, 100.0, 1.0)
        @test gamma_p1 > 0.0

        # Test with calculated alpha and cp (valid fsun240 conditions)
        (gamma_p2, cp2, alpha2) = CLM.get_gamma_P(
            200.0, 200.0, 200.0,
            50.0, 50.0, 50.0,
            0.5, 0.5,    # fsun240 in (0,1)
            100.0, 100.0, # positive forc values
            0.6)
        @test gamma_p2 > 0.0
        @test alpha2 != 0.001  # should be calculated, not fixed
    end

    @testset "get_gamma_T" begin
        # Standard conditions: t_veg=303.15K, t_veg24=303.15K, t_veg240=297K
        t_veg = 303.15
        t_veg24 = 303.15
        t_veg240 = 297.0
        ct1 = 95.0
        ct2 = 230.0
        betaT_val = 0.10
        LDF_val = 0.6
        Ceo_val = 2.0
        ivt = 1  # needleleaf evergreen (not boreal shrub/grass)

        (gamma_t, Eopt, topt) = CLM.get_gamma_T(
            t_veg240, t_veg24, t_veg,
            ct1, ct2, betaT_val, LDF_val, Ceo_val, ivt)

        @test gamma_t > 0.0
        @test topt ≈ 313.0 + 0.6 * (297.0 - 297.0) atol=1e-10  # co1 + co2*(t240-tstd0)
        @test Eopt > 0.0

        # Test with invalid t_veg240 (triggers fixed topt/Eopt)
        (gamma_t2, Eopt2, topt2) = CLM.get_gamma_T(
            0.0, t_veg24, t_veg,
            ct1, ct2, betaT_val, LDF_val, Ceo_val, ivt)
        @test topt2 ≈ 317.0  # topt_fix
        @test Eopt2 ≈ 2.26   # Eopt_fix

        # Test with arctic C3 grass (ivt=12)
        (gamma_t3, Eopt3, topt3) = CLM.get_gamma_T(
            280.0, 280.0, 280.0,
            ct1, ct2, betaT_val, LDF_val, Ceo_val, 12;
            nc3_arctic_grass=12)
        @test gamma_t3 > 0.0
        @test Eopt3 > 0.0

        # Test with boreal deciduous shrub (ivt=11)
        (gamma_t4, Eopt4, topt4) = CLM.get_gamma_T(
            280.0, 280.0, 280.0,
            ct1, ct2, betaT_val, LDF_val, Ceo_val, 11;
            nbrdlf_dcd_brl_shrub=11)
        @test gamma_t4 > 0.0
    end

    @testset "get_gamma_A" begin
        mf = CLM.MEGANFactors()
        CLM.megan_factors_init!(mf, 20)

        # Evergreen PFT (ivt=1) -> gamma_a = 1.0
        @test CLM.get_gamma_A(1, 2.0, 2.0, 1, mf) ≈ 1.0

        # Deciduous PFT (ivt=6, nbrdlf_dcd_trp_tree), steady LAI
        gamma_a = CLM.get_gamma_A(6, 2.0, 2.0, 1, mf; nbrdlf_dcd_trp_tree=6)
        @test gamma_a ≈ 1.0  # elai_prev == elai -> all mature

        # Deciduous PFT, growing LAI (elai > elai_prev)
        gamma_a2 = CLM.get_gamma_A(6, 1.5, 2.0, 1, mf; nbrdlf_dcd_trp_tree=6)
        # elai_prev = 2*1.5 - 2.0 = 1.0; fnew = 1 - 1.0/2.0 = 0.5; fmat = 0.5
        expected = 0.5 * mf.Anew[1] + 0.0 * mf.Agro[1] + 0.5 * mf.Amat[1] + 0.0 * mf.Aold[1]
        @test gamma_a2 ≈ expected

        # Deciduous PFT, declining LAI (elai < elai_prev)
        gamma_a3 = CLM.get_gamma_A(6, 2.5, 2.0, 1, mf; nbrdlf_dcd_trp_tree=6)
        # elai_prev = 2*2.5 - 2.0 = 3.0; fold = (3-2)/3 = 1/3; fmat = 1 - 1/3 = 2/3
        expected3 = 0.0 * mf.Anew[1] + 0.0 * mf.Agro[1] + (2.0/3.0) * mf.Amat[1] + (1.0/3.0) * mf.Aold[1]
        @test gamma_a3 ≈ expected3

        # Deciduous PFT with invalid elai240
        @test CLM.get_gamma_A(6, 0.0, 2.0, 1, mf; nbrdlf_dcd_trp_tree=6) ≈ 1.0
        @test CLM.get_gamma_A(6, 2.0e30, 2.0, 1, mf; nbrdlf_dcd_trp_tree=6) ≈ 1.0

        # ndllf_dcd_brl_tree (ivt=3) also non-evergreen
        @test CLM.get_gamma_A(3, 2.0, 2.0, 1, mf; ndllf_dcd_brl_tree=3) ≈ 1.0
    end

    @testset "get_gamma_C" begin
        # Standard conditions: ~400 ppm CO2, reasonable ci values
        forc_pbot = 101325.0  # 1 atm in Pa
        co2_ppmv = 400.0
        cisun = 0.7 * co2_ppmv * 1.0e-6 * forc_pbot  # typical ci
        cisha = cisun

        gamma_c = CLM.get_gamma_C(cisun, cisha, forc_pbot, 0.5, co2_ppmv)
        @test gamma_c > 0.0
        @test gamma_c < 2.0  # should be reasonable

        # All sunlit (fsun=1)
        gamma_c2 = CLM.get_gamma_C(cisun, 0.0, forc_pbot, 1.0, co2_ppmv)
        @test gamma_c2 > 0.0

        # All shaded (fsun=0)
        gamma_c3 = CLM.get_gamma_C(0.0, cisha, forc_pbot, 0.0, co2_ppmv)
        @test gamma_c3 > 0.0

        # NaN ci -> gamma_ci = 1.0
        gamma_c4 = CLM.get_gamma_C(NaN, NaN, forc_pbot, 0.5, co2_ppmv)
        # gamma_ca is still computed, so result != 1.0
        @test !isnan(gamma_c4)

        # Test different CO2 levels
        gamma_low = CLM.get_gamma_C(cisun, cisha, forc_pbot, 0.5, 300.0)
        gamma_high = CLM.get_gamma_C(cisun, cisha, forc_pbot, 0.5, 900.0)
        # Higher CO2 should generally decrease isoprene (gamma_C decreases)
        @test gamma_low > gamma_high
    end

    @testset "get_map_EF" begin
        efisop = zeros(6, 2)
        efisop[1, 1] = 10.0  # broadleaf trees
        efisop[2, 1] = 20.0  # fineleaf evergreen
        efisop[3, 1] = 30.0  # fineleaf deciduous
        efisop[4, 1] = 40.0  # shrubs
        efisop[5, 1] = 50.0  # grass
        efisop[6, 1] = 60.0  # crops

        # ndllf_evr_tmp_tree (ivt=1) -> fineleaf evergreen
        @test CLM.get_map_EF(1, 1, efisop) ≈ 20.0
        # ndllf_evr_brl_tree (ivt=2) -> fineleaf evergreen
        @test CLM.get_map_EF(2, 1, efisop) ≈ 20.0
        # ndllf_dcd_brl_tree (ivt=3) -> fineleaf deciduous
        @test CLM.get_map_EF(3, 1, efisop) ≈ 30.0
        # nbrdlf_evr_trp_tree (ivt=4) -> broadleaf trees
        @test CLM.get_map_EF(4, 1, efisop) ≈ 10.0
        # nbrdlf_evr_shrub (ivt=9) -> shrubs
        @test CLM.get_map_EF(9, 1, efisop) ≈ 40.0
        # nc3_arctic_grass (ivt=12) -> grass
        @test CLM.get_map_EF(12, 1, efisop) ≈ 50.0
        # nc3crop (ivt=15) -> crops
        @test CLM.get_map_EF(15, 1, efisop) ≈ 60.0
        # noveg (ivt=0) -> 0
        @test CLM.get_map_EF(0, 1, efisop) ≈ 0.0
    end

    # -----------------------------------------------------------------------
    # Integration test: voc_emission! driver
    # -----------------------------------------------------------------------
    @testset "voc_emission! driver integration" begin
        np = 4   # patches
        nc = 3   # columns
        ng = 2   # gridcells

        # Setup MEGAN compounds: isoprene + monoterpene
        iso = CLM.MEGANCompound(
            name = "isoprene",
            index = 1,
            class_number = 1,
            molec_weight = 68.12,
            coeff = 1.0,
            emis_factors = fill(600.0, 20)
        )
        mono = CLM.MEGANCompound(
            name = "myrcene",
            index = 2,
            class_number = 2,
            molec_weight = 136.23,
            coeff = 1.0,
            emis_factors = fill(100.0, 20)
        )
        meg_compounds = [iso, mono]

        # Mechanism compounds: ISOP maps to iso, TERP maps to mono
        mech_comps = [
            CLM.MEGANMechComp(name="ISOP", n_megan_comps=1, megan_indices=[1]),
            CLM.MEGANMechComp(name="TERP", n_megan_comps=1, megan_indices=[2]),
        ]

        # MEGAN factors
        mf = CLM.MEGANFactors()
        CLM.megan_factors_init!(mf, 20)

        # VOC data
        voc = CLM.VOCEmisData()
        CLM.vocemis_init!(voc, np, ng, 2, 2)

        # Set isoprene emission factors for gridcell 1
        voc.efisop_grc[:, 1] .= 600.0  # all categories
        voc.efisop_grc[:, 2] .= 500.0

        # Patch data
        patch = CLM.PatchData()
        CLM.patch_init!(patch, np)
        patch.itype .= [0, 1, 6, 12]       # noveg, ndllf_evr_tmp_tree, nbrdlf_dcd_trp_tree, nc3_arctic_grass
        patch.gridcell .= [1, 1, 1, 2]
        patch.column .= [1, 1, 2, 3]

        bounds_p = 1:np
        mask_soilp = trues(np)

        # Atmospheric forcing
        forc_solad_col = fill(200.0, nc, 2)  # W/m2 direct (vis, nir)
        forc_solai_grc = fill(100.0, ng, 2)  # W/m2 diffuse (vis, nir)
        forc_pbot_col = fill(101325.0, nc)   # Pa
        forc_pco2_grc = fill(40.53, ng)      # Pa (~400 ppm)

        fsd24_patch = fill(200.0, np)
        fsd240_patch = fill(200.0, np)
        fsi24_patch = fill(100.0, np)
        fsi240_patch = fill(100.0, np)

        # Canopy state
        fsun_patch = fill(0.5, np)
        fsun24_patch = fill(0.5, np)
        fsun240_patch = fill(0.5, np)
        elai_patch = fill(2.0, np)
        elai240_patch = fill(2.0, np)

        # Photosynthesis (ci ≈ 0.7 * co2)
        ci_val = 0.7 * 400.0e-6 * 101325.0
        cisun_z_patch = fill(ci_val, np, 1)
        cisha_z_patch = fill(ci_val, np, 1)

        # Temperature
        t_veg_patch = fill(300.0, np)
        t_veg24_patch = fill(300.0, np)
        t_veg240_patch = fill(297.0, np)

        # Energy flux
        btran_patch = fill(0.8, np)

        # Run the driver
        CLM.voc_emission!(
            voc, meg_compounds, mech_comps, mf,
            patch, bounds_p, mask_soilp,
            forc_solad_col, forc_solai_grc, forc_pbot_col, forc_pco2_grc,
            fsd24_patch, fsd240_patch, fsi24_patch, fsi240_patch,
            fsun_patch, fsun24_patch, fsun240_patch, elai_patch, elai240_patch,
            cisun_z_patch, cisha_z_patch,
            t_veg_patch, t_veg24_patch, t_veg240_patch,
            btran_patch
        )

        # Patch 1 (noveg, itype=0): no emissions
        @test voc.vocflx_tot_patch[1] ≈ 0.0

        # Patch 2 (ndllf_evr_tmp_tree, itype=1): should have emissions
        @test voc.vocflx_tot_patch[2] > 0.0
        @test voc.vocflx_patch[2, 1] > 0.0  # isoprene flux
        @test voc.vocflx_patch[2, 2] > 0.0  # monoterpene flux

        # Patch 3 (nbrdlf_dcd_trp_tree, itype=6): should have emissions
        @test voc.vocflx_tot_patch[3] > 0.0

        # Patch 4 (nc3_arctic_grass, itype=12): should have emissions
        @test voc.vocflx_tot_patch[4] > 0.0

        # Check that gamma diagnostics are set for first compound on patch 2
        @test !isnan(voc.gamma_out_patch[2])
        @test voc.gamma_out_patch[2] > 0.0
        @test !isnan(voc.gammaP_out_patch[2])
        @test !isnan(voc.gammaT_out_patch[2])
        @test !isnan(voc.gammaA_out_patch[2])
        @test !isnan(voc.gammaS_out_patch[2])
        @test !isnan(voc.gammaL_out_patch[2])
        @test !isnan(voc.gammaC_out_patch[2])

        # Check PAR diagnostics
        @test voc.paru_out_patch[2] > 0.0
        @test voc.para_out_patch[2] > 0.0

        # Check meg_flux_out is non-zero for vegetated patches
        @test voc.meg_flux_out[1, 2] > 0.0  # isoprene on patch 2
        @test voc.meg_flux_out[2, 2] > 0.0  # monoterpene on patch 2

        # Total flux = sum over mechanism compounds
        @test voc.vocflx_tot_patch[2] ≈ voc.vocflx_patch[2, 1] + voc.vocflx_patch[2, 2]
    end

    @testset "voc_emission! with no mechanism compounds" begin
        np = 2; ng = 1
        voc = CLM.VOCEmisData()
        CLM.vocemis_init!(voc, np, ng, 0, 0)

        patch = CLM.PatchData()
        CLM.patch_init!(patch, np)
        patch.itype .= [1, 2]
        patch.gridcell .= [1, 1]
        patch.column .= [1, 1]

        # Empty compounds => should return immediately
        result = CLM.voc_emission!(
            voc, CLM.MEGANCompound[], CLM.MEGANMechComp[], CLM.MEGANFactors(),
            patch, 1:np, trues(np),
            fill(200.0, 1, 2), fill(100.0, 1, 2),
            fill(101325.0, 1), fill(40.0, 1),
            fill(200.0, np), fill(200.0, np),
            fill(100.0, np), fill(100.0, np),
            fill(0.5, np), fill(0.5, np), fill(0.5, np),
            fill(2.0, np), fill(2.0, np),
            fill(28.0, np, 1), fill(28.0, np, 1),
            fill(300.0, np), fill(300.0, np), fill(297.0, np),
            fill(0.8, np)
        )
        @test result === nothing
    end

end
