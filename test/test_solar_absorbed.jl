@testset "SolarAbsorbedData" begin

    # Ensure varpar is initialized for dimension parameters
    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)

    numrad  = CLM.NUMRAD
    nlevcan = CLM.NLEVCAN
    nlevsno = CLM.varpar.nlevsno

    @testset "default construction" begin
        sa = CLM.SolarAbsorbedData()
        @test length(sa.fsa_patch) == 0
        @test length(sa.fsr_patch) == 0
        @test length(sa.sabg_patch) == 0
        @test length(sa.sabv_patch) == 0
        @test size(sa.parsun_z_patch) == (0, 0)
        @test size(sa.sabg_lyr_patch) == (0, 0)
        @test size(sa.sabs_roof_dir_lun) == (0, 0)
        @test length(sa.fsds_nir_d_patch) == 0
    end

    @testset "solarabs_init! without LUNA" begin
        np = 10
        nl = 4
        sa = CLM.SolarAbsorbedData()
        CLM.solarabs_init!(sa, np, nl)

        # --- Check patch-level 1D sizes ---
        @test length(sa.fsa_patch) == np
        @test length(sa.fsa_u_patch) == np
        @test length(sa.fsa_r_patch) == np
        @test length(sa.sabv_patch) == np
        @test length(sa.sabg_patch) == np
        @test length(sa.sabg_pen_patch) == np
        @test length(sa.sabg_soil_patch) == np
        @test length(sa.sabg_snow_patch) == np
        @test length(sa.sabg_chk_patch) == np
        @test length(sa.sub_surf_abs_SW_patch) == np
        @test length(sa.fsr_patch) == np
        @test length(sa.fsrSF_patch) == np
        @test length(sa.ssre_fsr_patch) == np

        # --- Check canopy PAR 2D sizes ---
        @test size(sa.parsun_z_patch) == (np, nlevcan)
        @test size(sa.parsha_z_patch) == (np, nlevcan)

        # --- LUNA arrays should be empty when use_luna=false ---
        @test size(sa.par240d_z_patch) == (0, 0)
        @test size(sa.par240x_z_patch) == (0, 0)
        @test size(sa.par24d_z_patch) == (0, 0)
        @test size(sa.par24x_z_patch) == (0, 0)

        # --- Check snow layer absorption sizes ---
        nlev_sno1 = nlevsno + 1
        @test size(sa.sabg_lyr_patch) == (np, nlev_sno1)

        # --- Check landunit-level 2D sizes ---
        @test size(sa.sabs_roof_dir_lun) == (nl, numrad)
        @test size(sa.sabs_roof_dif_lun) == (nl, numrad)
        @test size(sa.sabs_sunwall_dir_lun) == (nl, numrad)
        @test size(sa.sabs_sunwall_dif_lun) == (nl, numrad)
        @test size(sa.sabs_shadewall_dir_lun) == (nl, numrad)
        @test size(sa.sabs_shadewall_dif_lun) == (nl, numrad)
        @test size(sa.sabs_improad_dir_lun) == (nl, numrad)
        @test size(sa.sabs_improad_dif_lun) == (nl, numrad)
        @test size(sa.sabs_perroad_dir_lun) == (nl, numrad)
        @test size(sa.sabs_perroad_dif_lun) == (nl, numrad)

        # --- Check NIR diagnostic sizes ---
        @test length(sa.fsds_nir_d_patch) == np
        @test length(sa.fsds_nir_i_patch) == np
        @test length(sa.fsds_nir_d_ln_patch) == np
        @test length(sa.fsr_nir_d_patch) == np
        @test length(sa.fsr_nir_i_patch) == np
        @test length(sa.fsr_nir_d_ln_patch) == np
        @test length(sa.fsrSF_nir_d_patch) == np
        @test length(sa.fsrSF_nir_i_patch) == np
        @test length(sa.fsrSF_nir_d_ln_patch) == np
        @test length(sa.ssre_fsr_nir_d_patch) == np
        @test length(sa.ssre_fsr_nir_i_patch) == np
        @test length(sa.ssre_fsr_nir_d_ln_patch) == np

        # --- Check NaN initialization ---
        @test all(isnan, sa.fsa_patch)
        @test all(isnan, sa.sabv_patch)
        @test all(isnan, sa.sabg_patch)
        @test all(isnan, sa.fsr_patch)
        @test all(isnan, sa.parsun_z_patch)
        @test all(isnan, sa.sabg_lyr_patch)
        @test all(isnan, sa.sabs_roof_dir_lun)
        @test all(isnan, sa.fsds_nir_d_patch)
        @test all(isnan, sa.ssre_fsr_nir_d_ln_patch)
    end

    @testset "solarabs_init! with LUNA" begin
        np = 8
        nl = 3
        sa = CLM.SolarAbsorbedData()
        CLM.solarabs_init!(sa, np, nl; use_luna=true)

        # LUNA arrays should be allocated and initialized to SPVAL
        @test size(sa.par240d_z_patch) == (np, nlevcan)
        @test size(sa.par240x_z_patch) == (np, nlevcan)
        @test size(sa.par24d_z_patch) == (np, nlevcan)
        @test size(sa.par24x_z_patch) == (np, nlevcan)
        @test all(x -> x == CLM.SPVAL, sa.par240d_z_patch)
        @test all(x -> x == CLM.SPVAL, sa.par240x_z_patch)
        @test all(x -> x == CLM.SPVAL, sa.par24d_z_patch)
        @test all(x -> x == CLM.SPVAL, sa.par24x_z_patch)

        # Non-LUNA arrays should still be NaN
        @test all(isnan, sa.parsun_z_patch)
        @test all(isnan, sa.parsha_z_patch)
    end

    @testset "solarabs_clean!" begin
        sa = CLM.SolarAbsorbedData()
        CLM.solarabs_init!(sa, 10, 4; use_luna=true)
        CLM.solarabs_clean!(sa)

        @test length(sa.fsa_patch) == 0
        @test length(sa.fsr_patch) == 0
        @test length(sa.sabg_patch) == 0
        @test length(sa.sabv_patch) == 0
        @test length(sa.sub_surf_abs_SW_patch) == 0
        @test size(sa.parsun_z_patch) == (0, 0)
        @test size(sa.par240d_z_patch) == (0, 0)
        @test size(sa.sabg_lyr_patch) == (0, 0)
        @test size(sa.sabs_roof_dir_lun) == (0, 0)
        @test size(sa.sabs_shadewall_dif_lun) == (0, 0)
        @test length(sa.fsds_nir_d_patch) == 0
        @test length(sa.fsr_nir_d_patch) == 0
        @test length(sa.fsrSF_nir_d_patch) == 0
        @test length(sa.ssre_fsr_nir_d_ln_patch) == 0
    end

    @testset "solarabs_init_cold!" begin
        np = 6
        nl = 3
        sa = CLM.SolarAbsorbedData()
        CLM.solarabs_init!(sa, np, nl)
        CLM.solarabs_init_cold!(sa, 1:nl)

        # All urban surface absorption fields should be zero
        for l in 1:nl
            for ib in 1:numrad
                @test sa.sabs_roof_dir_lun[l, ib] == 0.0
                @test sa.sabs_roof_dif_lun[l, ib] == 0.0
                @test sa.sabs_sunwall_dir_lun[l, ib] == 0.0
                @test sa.sabs_sunwall_dif_lun[l, ib] == 0.0
                @test sa.sabs_shadewall_dir_lun[l, ib] == 0.0
                @test sa.sabs_shadewall_dif_lun[l, ib] == 0.0
                @test sa.sabs_improad_dir_lun[l, ib] == 0.0
                @test sa.sabs_improad_dif_lun[l, ib] == 0.0
                @test sa.sabs_perroad_dir_lun[l, ib] == 0.0
                @test sa.sabs_perroad_dif_lun[l, ib] == 0.0
            end
        end
    end

    @testset "field mutability" begin
        sa = CLM.SolarAbsorbedData()
        CLM.solarabs_init!(sa, 5, 3)

        # Write and verify patch 1D field
        sa.fsa_patch[2] = 150.0
        @test sa.fsa_patch[2] == 150.0

        # Write and verify patch 2D field
        sa.parsun_z_patch[3, 1] = 42.0
        @test sa.parsun_z_patch[3, 1] == 42.0

        # Write and verify landunit 2D field
        sa.sabs_roof_dir_lun[1, 2] = 0.75
        @test sa.sabs_roof_dir_lun[1, 2] == 0.75

        # Write and verify snow layer field
        sa.sabg_lyr_patch[1, 1] = 0.99
        @test sa.sabg_lyr_patch[1, 1] == 0.99
    end

    @testset "re-init overwrites previous state" begin
        sa = CLM.SolarAbsorbedData()
        CLM.solarabs_init!(sa, 5, 3)
        sa.fsa_patch[1] = 999.0

        CLM.solarabs_init!(sa, 10, 6)
        @test length(sa.fsa_patch) == 10
        @test all(isnan, sa.fsa_patch)
        @test size(sa.sabs_roof_dir_lun) == (6, numrad)
    end

    @testset "stub functions run without error" begin
        sa = CLM.SolarAbsorbedData()
        CLM.solarabs_init!(sa, 6, 3)

        @test CLM.solarabs_init_history!(sa, 1:6) === nothing
        @test CLM.solarabs_restart!(sa, 1:3, 1:6) === nothing
    end

    @testset "solarabs_init_history! sets SPVAL" begin
        sa = CLM.SolarAbsorbedData()
        CLM.solarabs_init!(sa, 6, 3)
        CLM.solarabs_init_history!(sa, 1:6)

        for p in 1:6
            @test sa.fsa_patch[p] == CLM.SPVAL
            @test sa.fsr_patch[p] == CLM.SPVAL
            @test sa.sabv_patch[p] == CLM.SPVAL
            @test sa.sabg_patch[p] == CLM.SPVAL
            @test sa.sabg_pen_patch[p] == CLM.SPVAL
            @test sa.fsds_nir_d_patch[p] == CLM.SPVAL
            @test sa.fsr_nir_d_patch[p] == CLM.SPVAL
            @test sa.sub_surf_abs_SW_patch[p] == CLM.SPVAL
        end

        # Snow layers should be set to SPVAL
        for p in 1:6
            for j in 1:nlevsno
                @test sa.sabg_lyr_patch[p, j] == CLM.SPVAL
            end
        end
    end

    @testset "solarabs_init_history! with use_SSRE" begin
        sa = CLM.SolarAbsorbedData()
        CLM.solarabs_init!(sa, 4, 2)
        CLM.solarabs_init_history!(sa, 1:4; use_SSRE=true)

        for p in 1:4
            @test sa.fsrSF_patch[p] == CLM.SPVAL
            @test sa.ssre_fsr_patch[p] == CLM.SPVAL
            @test sa.fsrSF_nir_d_patch[p] == CLM.SPVAL
            @test sa.ssre_fsr_nir_d_patch[p] == CLM.SPVAL
            @test sa.ssre_fsr_nir_d_ln_patch[p] == CLM.SPVAL
        end
    end

end
