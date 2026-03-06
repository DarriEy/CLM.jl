@testset "CanopyStateData" begin

    # Ensure varpar is initialized for dimension parameters
    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)

    nlevcan = CLM.NLEVCAN
    nvegwcs = CLM.NVEGWCS

    @testset "default construction" begin
        cs = CLM.CanopyStateData()
        @test length(cs.frac_veg_nosno_patch) == 0
        @test length(cs.frac_veg_nosno_alb_patch) == 0
        @test length(cs.tlai_patch) == 0
        @test length(cs.tsai_patch) == 0
        @test length(cs.elai_patch) == 0
        @test length(cs.esai_patch) == 0
        @test length(cs.tlai_hist_patch) == 0
        @test length(cs.tsai_hist_patch) == 0
        @test length(cs.htop_hist_patch) == 0
        @test length(cs.elai240_patch) == 0
        @test length(cs.laisun_patch) == 0
        @test length(cs.laisha_patch) == 0
        @test length(cs.mlaidiff_patch) == 0
        @test length(cs.stem_biomass_patch) == 0
        @test length(cs.leaf_biomass_patch) == 0
        @test length(cs.htop_patch) == 0
        @test length(cs.hbot_patch) == 0
        @test length(cs.z0m_patch) == 0
        @test length(cs.displa_patch) == 0
        @test length(cs.fsun_patch) == 0
        @test length(cs.fsun24_patch) == 0
        @test length(cs.fsun240_patch) == 0
        @test length(cs.dleaf_patch) == 0
        @test length(cs.rscanopy_patch) == 0
        @test size(cs.laisun_z_patch) == (0, 0)
        @test size(cs.laisha_z_patch) == (0, 0)
        @test size(cs.annlai_patch) == (0, 0)
        @test size(cs.vegwp_patch) == (0, 0)
        @test size(cs.vegwp_ln_patch) == (0, 0)
        @test size(cs.vegwp_pd_patch) == (0, 0)
        @test cs.leaf_mr_vcm == CLM.SPVAL
    end

    @testset "canopystate_init!" begin
        np = 10
        cs = CLM.CanopyStateData()
        CLM.canopystate_init!(cs, np)

        # --- Check 1D integer sizes and values ---
        @test length(cs.frac_veg_nosno_patch) == np
        @test all(x -> x == 0, cs.frac_veg_nosno_patch)
        @test length(cs.frac_veg_nosno_alb_patch) == np
        @test all(x -> x == 0, cs.frac_veg_nosno_alb_patch)

        # --- Check 1D float sizes ---
        @test length(cs.tlai_patch) == np
        @test length(cs.tsai_patch) == np
        @test length(cs.elai_patch) == np
        @test length(cs.esai_patch) == np
        @test length(cs.tlai_hist_patch) == np
        @test length(cs.tsai_hist_patch) == np
        @test length(cs.htop_hist_patch) == np
        @test length(cs.elai240_patch) == np
        @test length(cs.laisun_patch) == np
        @test length(cs.laisha_patch) == np
        @test length(cs.mlaidiff_patch) == np
        @test length(cs.stem_biomass_patch) == np
        @test length(cs.leaf_biomass_patch) == np
        @test length(cs.htop_patch) == np
        @test length(cs.hbot_patch) == np
        @test length(cs.z0m_patch) == np
        @test length(cs.displa_patch) == np
        @test length(cs.fsun_patch) == np
        @test length(cs.fsun24_patch) == np
        @test length(cs.fsun240_patch) == np
        @test length(cs.dleaf_patch) == np
        @test length(cs.rscanopy_patch) == np

        # --- Check 2D sizes ---
        @test size(cs.laisun_z_patch) == (np, nlevcan)
        @test size(cs.laisha_z_patch) == (np, nlevcan)
        @test size(cs.annlai_patch) == (np, 12)
        @test size(cs.vegwp_patch) == (np, nvegwcs)
        @test size(cs.vegwp_ln_patch) == (np, nvegwcs)
        @test size(cs.vegwp_pd_patch) == (np, nvegwcs)

        # --- Check NaN initialization ---
        @test all(isnan, cs.tlai_patch)
        @test all(isnan, cs.tsai_patch)
        @test all(isnan, cs.elai_patch)
        @test all(isnan, cs.esai_patch)
        @test all(isnan, cs.tlai_hist_patch)
        @test all(isnan, cs.tsai_hist_patch)
        @test all(isnan, cs.htop_hist_patch)
        @test all(isnan, cs.elai240_patch)
        @test all(isnan, cs.laisun_patch)
        @test all(isnan, cs.laisha_patch)
        @test all(isnan, cs.mlaidiff_patch)
        @test all(isnan, cs.stem_biomass_patch)
        @test all(isnan, cs.leaf_biomass_patch)
        @test all(isnan, cs.htop_patch)
        @test all(isnan, cs.hbot_patch)
        @test all(isnan, cs.z0m_patch)
        @test all(isnan, cs.displa_patch)
        @test all(isnan, cs.fsun_patch)
        @test all(isnan, cs.fsun24_patch)
        @test all(isnan, cs.fsun240_patch)
        @test all(isnan, cs.dleaf_patch)
        @test all(isnan, cs.rscanopy_patch)
        @test all(isnan, cs.laisun_z_patch)
        @test all(isnan, cs.laisha_z_patch)
        @test all(isnan, cs.annlai_patch)
        @test all(isnan, cs.vegwp_patch)
        @test all(isnan, cs.vegwp_ln_patch)
        @test all(isnan, cs.vegwp_pd_patch)
    end

    @testset "canopystate_clean!" begin
        cs = CLM.CanopyStateData()
        CLM.canopystate_init!(cs, 5)
        CLM.canopystate_clean!(cs)

        # All vectors should be empty
        @test length(cs.frac_veg_nosno_patch) == 0
        @test length(cs.frac_veg_nosno_alb_patch) == 0
        @test length(cs.tlai_patch) == 0
        @test length(cs.tsai_patch) == 0
        @test length(cs.elai_patch) == 0
        @test length(cs.esai_patch) == 0
        @test length(cs.tlai_hist_patch) == 0
        @test length(cs.tsai_hist_patch) == 0
        @test length(cs.htop_hist_patch) == 0
        @test length(cs.elai240_patch) == 0
        @test length(cs.laisun_patch) == 0
        @test length(cs.laisha_patch) == 0
        @test length(cs.mlaidiff_patch) == 0
        @test length(cs.stem_biomass_patch) == 0
        @test length(cs.leaf_biomass_patch) == 0
        @test length(cs.htop_patch) == 0
        @test length(cs.hbot_patch) == 0
        @test length(cs.z0m_patch) == 0
        @test length(cs.displa_patch) == 0
        @test length(cs.fsun_patch) == 0
        @test length(cs.fsun24_patch) == 0
        @test length(cs.fsun240_patch) == 0
        @test length(cs.dleaf_patch) == 0
        @test length(cs.rscanopy_patch) == 0

        # All matrices should be empty
        @test size(cs.laisun_z_patch) == (0, 0)
        @test size(cs.laisha_z_patch) == (0, 0)
        @test size(cs.annlai_patch) == (0, 0)
        @test size(cs.vegwp_patch) == (0, 0)
        @test size(cs.vegwp_ln_patch) == (0, 0)
        @test size(cs.vegwp_pd_patch) == (0, 0)
    end

    @testset "canopystate_init_cold! (basic)" begin
        np = 5
        cs = CLM.CanopyStateData()
        CLM.canopystate_init!(cs, np)

        CLM.canopystate_init_cold!(cs, 1:np)

        for p in 1:np
            @test cs.tlai_patch[p] ≈ 0.0
            @test cs.tsai_patch[p] ≈ 0.0
            @test cs.elai_patch[p] ≈ 0.0
            @test cs.esai_patch[p] ≈ 0.0
            @test cs.stem_biomass_patch[p] ≈ 0.0
            @test cs.leaf_biomass_patch[p] ≈ 0.0
            @test cs.htop_patch[p] ≈ 0.0
            @test cs.hbot_patch[p] ≈ 0.0
            @test all(x -> x ≈ -2.5e4, cs.vegwp_patch[p, :])
            @test cs.laisun_patch[p] ≈ 0.0
            @test cs.laisha_patch[p] ≈ 0.0
            @test cs.tlai_hist_patch[p] ≈ 0.0
            @test cs.tsai_hist_patch[p] ≈ 0.0
            @test cs.htop_hist_patch[p] ≈ 0.0
            @test cs.fsun_patch[p] == CLM.SPVAL
        end

        # Other fields should remain at NaN
        @test all(isnan, cs.z0m_patch)
        @test all(isnan, cs.displa_patch)
        @test all(isnan, cs.dleaf_patch)
    end

    @testset "canopystate_init_cold! (with landunit info)" begin
        np = 4
        cs = CLM.CanopyStateData()
        CLM.canopystate_init!(cs, np)

        # Two landunits: soil (1) and ice (4)
        landunit_patch = [1, 1, 2, 2]
        lun_itype = [CLM.ISTSOIL, CLM.ISTICE]

        CLM.canopystate_init_cold!(cs, 1:np;
            landunit_patch=landunit_patch,
            lun_itype=lun_itype)

        # Soil patches should have laisun/laisha = 0
        @test cs.laisun_patch[1] ≈ 0.0
        @test cs.laisha_patch[1] ≈ 0.0
        @test cs.laisun_patch[2] ≈ 0.0
        @test cs.laisha_patch[2] ≈ 0.0

        # Ice patches should NOT have laisun/laisha set (remain NaN)
        @test isnan(cs.laisun_patch[3])
        @test isnan(cs.laisha_patch[3])
        @test isnan(cs.laisun_patch[4])
        @test isnan(cs.laisha_patch[4])
    end

    @testset "canopystate_read_nml!" begin
        cs = CLM.CanopyStateData()
        @test cs.leaf_mr_vcm == CLM.SPVAL

        CLM.canopystate_read_nml!(cs; leaf_mr_vcm=0.015)
        @test cs.leaf_mr_vcm ≈ 0.015

        CLM.canopystate_read_nml!(cs; leaf_mr_vcm=0.02)
        @test cs.leaf_mr_vcm ≈ 0.02
    end

    @testset "canopystate_set_nml_for_testing!" begin
        cs = CLM.CanopyStateData()
        CLM.canopystate_set_nml_for_testing!(cs)
        @test cs.leaf_mr_vcm ≈ 0.015
    end

    @testset "field mutability" begin
        cs = CLM.CanopyStateData()
        CLM.canopystate_init!(cs, 3)

        # Write and verify 1D fields
        cs.tlai_patch[1] = 3.5
        @test cs.tlai_patch[1] == 3.5

        cs.fsun_patch[2] = 0.7
        @test cs.fsun_patch[2] == 0.7

        cs.htop_patch[3] = 15.0
        @test cs.htop_patch[3] == 15.0

        cs.frac_veg_nosno_patch[1] = 1
        @test cs.frac_veg_nosno_patch[1] == 1

        # Write and verify 2D fields
        cs.laisun_z_patch[1, 1] = 2.0
        @test cs.laisun_z_patch[1, 1] == 2.0

        cs.vegwp_patch[2, 3] = -1000.0
        @test cs.vegwp_patch[2, 3] == -1000.0

        cs.annlai_patch[1, 6] = 4.5
        @test cs.annlai_patch[1, 6] == 4.5
    end

    @testset "re-init overwrites previous state" begin
        cs = CLM.CanopyStateData()
        CLM.canopystate_init!(cs, 3)
        cs.tlai_patch[1] = 999.0

        CLM.canopystate_init!(cs, 7)
        @test length(cs.tlai_patch) == 7
        @test all(isnan, cs.tlai_patch)
        @test size(cs.laisun_z_patch) == (7, nlevcan)
        @test size(cs.vegwp_patch) == (7, nvegwcs)
        @test size(cs.annlai_patch) == (7, 12)
    end

    @testset "stub functions run without error" begin
        cs = CLM.CanopyStateData()
        CLM.canopystate_init!(cs, 5)

        @test CLM.canopystate_init_history!(cs, 1:5) === nothing
        @test CLM.canopystate_restart!(cs, 1:5) === nothing
        @test CLM.canopystate_init_acc_vars!(cs, 1:5) === nothing
        @test CLM.canopystate_update_acc_vars!(cs, 1:5) === nothing
    end

    @testset "canopystate_init_acc_buffer!" begin
        cs = CLM.CanopyStateData()
        CLM.canopystate_init!(cs, 5)

        CLM.canopystate_init_acc_buffer!(cs, 1:5)
        for p in 1:5
            @test cs.fsun24_patch[p] == CLM.SPVAL
            @test cs.fsun240_patch[p] == CLM.SPVAL
            @test cs.elai240_patch[p] == CLM.SPVAL
        end
    end

end
