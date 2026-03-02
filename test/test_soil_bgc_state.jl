@testset "SoilBiogeochemStateData" begin

    @testset "default construction" begin
        st = CLM.SoilBiogeochemStateData()
        @test length(st.fpi_col) == 0
        @test length(st.fpg_col) == 0
        @test length(st.plant_ndemand_col) == 0
        @test length(st.nue_decomp_cascade_col) == 0
        @test size(st.leaf_prof_patch) == (0, 0)
        @test size(st.fpi_vr_col) == (0, 0)
        @test size(st.som_adv_coef_col) == (0, 0)
    end

    @testset "init! basic allocation" begin
        nc = 5
        np = 8
        nlev_full = 3
        ntrans = 4

        st = CLM.SoilBiogeochemStateData()
        CLM.soil_bgc_state_init!(st, nc, np, nlev_full, ntrans)

        # Patch-level 2D (patch × nlevdecomp_full) — initialized to SPVAL
        @test size(st.leaf_prof_patch) == (np, nlev_full)
        @test all(x -> x == CLM.SPVAL, st.leaf_prof_patch)
        @test size(st.froot_prof_patch) == (np, nlev_full)
        @test all(x -> x == CLM.SPVAL, st.froot_prof_patch)
        @test size(st.croot_prof_patch) == (np, nlev_full)
        @test all(x -> x == CLM.SPVAL, st.croot_prof_patch)
        @test size(st.stem_prof_patch) == (np, nlev_full)
        @test all(x -> x == CLM.SPVAL, st.stem_prof_patch)

        # Column-level 2D (col × nlevdecomp_full)
        @test size(st.fpi_vr_col) == (nc, nlev_full)
        @test all(isnan, st.fpi_vr_col)
        @test size(st.nfixation_prof_col) == (nc, nlev_full)
        @test all(x -> x == CLM.SPVAL, st.nfixation_prof_col)
        @test size(st.ndep_prof_col) == (nc, nlev_full)
        @test all(x -> x == CLM.SPVAL, st.ndep_prof_col)
        @test size(st.som_adv_coef_col) == (nc, nlev_full)
        @test all(x -> x == CLM.SPVAL, st.som_adv_coef_col)
        @test size(st.som_diffus_coef_col) == (nc, nlev_full)
        @test all(x -> x == CLM.SPVAL, st.som_diffus_coef_col)

        # Column-level 1D — initialized to NaN
        @test length(st.fpi_col) == nc
        @test all(isnan, st.fpi_col)
        @test length(st.fpg_col) == nc
        @test all(isnan, st.fpg_col)
        @test length(st.plant_ndemand_col) == nc
        @test all(isnan, st.plant_ndemand_col)

        # Transition-level 1D — initialized to NaN
        @test length(st.nue_decomp_cascade_col) == ntrans
        @test all(isnan, st.nue_decomp_cascade_col)
    end

    @testset "init_cold! special landunits" begin
        nc = 4
        np = 6
        nlev_full = 3
        ntrans = 2

        st = CLM.SoilBiogeochemStateData()
        CLM.soil_bgc_state_init!(st, nc, np, nlev_full, ntrans)

        # Columns 1,3 are special; 2,4 are soil/crop
        mask_special   = BitVector([true, false, true, false])
        mask_soil_crop = BitVector([false, true, false, true])

        CLM.soil_bgc_state_init_cold!(st, 1:nc, nlev_full;
                                       mask_special=mask_special,
                                       mask_soil_crop=mask_soil_crop)

        # Special columns should have SPVAL for fpi/fpg
        @test st.fpi_col[1] == CLM.SPVAL
        @test st.fpg_col[1] == CLM.SPVAL
        @test st.fpi_col[3] == CLM.SPVAL
        @test st.fpg_col[3] == CLM.SPVAL
        for j in 1:nlev_full
            @test st.fpi_vr_col[1, j] == CLM.SPVAL
            @test st.fpi_vr_col[3, j] == CLM.SPVAL
        end

        # Non-special columns should still be NaN (not touched by special branch)
        @test isnan(st.fpi_col[2])
        @test isnan(st.fpg_col[2])
    end

    @testset "init_cold! soil/crop columns zeroed" begin
        nc = 3
        np = 4
        nlev_full = 2
        ntrans = 1

        st = CLM.SoilBiogeochemStateData()
        CLM.soil_bgc_state_init!(st, nc, np, nlev_full, ntrans)

        mask_soil_crop = BitVector([true, true, false])

        CLM.soil_bgc_state_init_cold!(st, 1:nc, nlev_full;
                                       mask_soil_crop=mask_soil_crop)

        # Soil/crop columns: fields zeroed out
        for c in [1, 2]
            for j in 1:nlev_full
                @test st.fpi_vr_col[c, j] == 0.0
                @test st.som_adv_coef_col[c, j] == 0.0
                @test st.som_diffus_coef_col[c, j] == 0.0
                @test st.nfixation_prof_col[c, j] == 0.0
                @test st.ndep_prof_col[c, j] == 0.0
            end
        end

        # Non-soil/crop column 3: unchanged (still SPVAL/NaN from init)
        for j in 1:nlev_full
            @test isnan(st.fpi_vr_col[3, j])
            @test st.som_adv_coef_col[3, j] == CLM.SPVAL
            @test st.som_diffus_coef_col[3, j] == CLM.SPVAL
            @test st.nfixation_prof_col[3, j] == CLM.SPVAL
            @test st.ndep_prof_col[3, j] == CLM.SPVAL
        end
    end

    @testset "init_cold! default mask_soil_crop (all columns)" begin
        nc = 2
        np = 3
        nlev_full = 2
        ntrans = 1

        st = CLM.SoilBiogeochemStateData()
        CLM.soil_bgc_state_init!(st, nc, np, nlev_full, ntrans)

        # No mask_soil_crop => all columns treated as soil/crop
        CLM.soil_bgc_state_init_cold!(st, 1:nc, nlev_full)

        for c in 1:nc
            for j in 1:nlev_full
                @test st.fpi_vr_col[c, j] == 0.0
                @test st.som_adv_coef_col[c, j] == 0.0
                @test st.som_diffus_coef_col[c, j] == 0.0
                @test st.nfixation_prof_col[c, j] == 0.0
                @test st.ndep_prof_col[c, j] == 0.0
            end
        end
    end

    @testset "get_spinup_latitude_term" begin
        # At latitude 60, logistic term at midpoint: 1 + 50/(1+exp(0)) = 1 + 25 = 26
        @test CLM.get_spinup_latitude_term(60.0) ≈ 26.0

        # Symmetric for negative latitudes
        @test CLM.get_spinup_latitude_term(-60.0) ≈ CLM.get_spinup_latitude_term(60.0)

        # Low latitude → term close to 1 (logistic near 0)
        low_lat = CLM.get_spinup_latitude_term(0.0)
        @test low_lat > 1.0
        @test low_lat < 5.0

        # High latitude → term approaches 51
        high_lat = CLM.get_spinup_latitude_term(90.0)
        @test high_lat > 40.0
        @test high_lat <= 51.0
    end

    @testset "stub functions run without error" begin
        st = CLM.SoilBiogeochemStateData()
        CLM.soil_bgc_state_init!(st, 5, 8, 3, 2)

        @test CLM.soil_bgc_state_init_history!(st, 1:5) === nothing
        @test CLM.soil_bgc_state_restart!(st, 1:5) === nothing
    end

    @testset "field mutability" begin
        st = CLM.SoilBiogeochemStateData()
        CLM.soil_bgc_state_init!(st, 3, 4, 2, 1)

        st.fpi_col[1] = 0.5
        @test st.fpi_col[1] == 0.5

        st.fpg_col[2] = 0.8
        @test st.fpg_col[2] == 0.8

        st.leaf_prof_patch[1, 1] = 42.0
        @test st.leaf_prof_patch[1, 1] == 42.0

        st.som_adv_coef_col[1, 1] = 1.5e-5
        @test st.som_adv_coef_col[1, 1] == 1.5e-5

        st.nue_decomp_cascade_col[1] = 0.9
        @test st.nue_decomp_cascade_col[1] == 0.9
    end

    @testset "re-init overwrites previous state" begin
        st = CLM.SoilBiogeochemStateData()
        CLM.soil_bgc_state_init!(st, 3, 4, 2, 1)
        st.fpi_col[1] = 999.0

        CLM.soil_bgc_state_init!(st, 5, 8, 4, 3)
        @test length(st.fpi_col) == 5
        @test all(isnan, st.fpi_col)
        @test size(st.leaf_prof_patch) == (8, 4)
        @test size(st.fpi_vr_col) == (5, 4)
        @test length(st.nue_decomp_cascade_col) == 3
    end

end
