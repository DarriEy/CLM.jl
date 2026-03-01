@testset "SoilStateData" begin

    # Ensure varpar is initialized for dimension parameters
    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)

    nlevsoi        = CLM.varpar.nlevsoi
    nlevgrnd       = CLM.varpar.nlevgrnd
    nlevsno        = CLM.varpar.nlevsno
    nlevmaxurbgrnd = CLM.varpar.nlevmaxurbgrnd
    nlev_soisno    = nlevsno + nlevmaxurbgrnd

    @testset "default construction" begin
        ss = CLM.SoilStateData()
        @test length(ss.sandfrac_patch) == 0
        @test length(ss.clayfrac_patch) == 0
        @test length(ss.mss_frc_cly_vld_col) == 0
        @test length(ss.smpmin_col) == 0
        @test length(ss.dsl_col) == 0
        @test length(ss.soilbeta_col) == 0
        @test length(ss.root_depth_patch) == 0
        @test size(ss.cellorg_col) == (0, 0)
        @test size(ss.hksat_col) == (0, 0)
        @test size(ss.watsat_col) == (0, 0)
        @test size(ss.thk_col) == (0, 0)
        @test size(ss.rootr_patch) == (0, 0)
        @test size(ss.rootfr_patch) == (0, 0)
        @test size(ss.msw_col) == (0, 0)
        @test size(ss.k_soil_root_patch) == (0, 0)
    end

    @testset "soilstate_init!" begin
        np = 20   # patches
        nc = 10   # columns

        ss = CLM.SoilStateData()
        CLM.soilstate_init!(ss, np, nc)

        # --- Check patch-level 1D sizes ---
        @test length(ss.sandfrac_patch) == np
        @test length(ss.clayfrac_patch) == np
        @test length(ss.root_depth_patch) == np

        # --- Check column-level 1D sizes ---
        @test length(ss.mss_frc_cly_vld_col) == nc
        @test length(ss.smpmin_col) == nc
        @test length(ss.dsl_col) == nc
        @test length(ss.soilresis_col) == nc
        @test length(ss.soilbeta_col) == nc
        @test length(ss.soilalpha_col) == nc
        @test length(ss.soilalpha_u_col) == nc
        @test length(ss.wtfact_col) == nc
        @test length(ss.gwc_thr_col) == nc

        # --- Check column-level 2D sizes (nlevsoi) ---
        @test size(ss.cellorg_col) == (nc, nlevsoi)
        @test size(ss.cellsand_col) == (nc, nlevsoi)
        @test size(ss.cellclay_col) == (nc, nlevsoi)

        # --- Check column-level 2D sizes (nlevgrnd) ---
        @test size(ss.bd_col) == (nc, nlevgrnd)
        @test size(ss.hksat_col) == (nc, nlevgrnd)
        @test size(ss.hksat_min_col) == (nc, nlevgrnd)
        @test size(ss.hk_l_col) == (nc, nlevgrnd)
        @test size(ss.smp_l_col) == (nc, nlevgrnd)
        @test size(ss.bsw_col) == (nc, nlevgrnd)
        @test size(ss.watdry_col) == (nc, nlevgrnd)
        @test size(ss.watopt_col) == (nc, nlevgrnd)
        @test size(ss.watfc_col) == (nc, nlevgrnd)
        @test size(ss.sucsat_col) == (nc, nlevgrnd)
        @test size(ss.soilpsi_col) == (nc, nlevgrnd)
        @test size(ss.eff_porosity_col) == (nc, nlevgrnd)
        @test size(ss.tkmg_col) == (nc, nlevgrnd)
        @test size(ss.tkdry_col) == (nc, nlevgrnd)
        @test size(ss.tksatu_col) == (nc, nlevgrnd)
        @test size(ss.csol_col) == (nc, nlevgrnd)

        # --- Check column-level 2D sizes (nlevmaxurbgrnd) ---
        @test size(ss.watsat_col) == (nc, nlevmaxurbgrnd)

        # --- Check column-level 2D sizes (nlayer) ---
        @test size(ss.porosity_col) == (nc, CLM.NLAYER)

        # --- Check column-level 2D sizes (nlevsno+nlevmaxurbgrnd) ---
        @test size(ss.thk_col) == (nc, nlev_soisno)

        # --- Check van Genuchten parameter sizes ---
        @test size(ss.msw_col) == (nc, nlevgrnd)
        @test size(ss.nsw_col) == (nc, nlevgrnd)
        @test size(ss.alphasw_col) == (nc, nlevgrnd)
        @test size(ss.watres_col) == (nc, nlevgrnd)

        # --- Check root sizes (patch-level 2D) ---
        @test size(ss.rootr_patch) == (np, nlevgrnd)
        @test size(ss.rootfr_patch) == (np, nlevgrnd)
        @test size(ss.crootfr_patch) == (np, nlevgrnd)
        @test size(ss.k_soil_root_patch) == (np, nlevsoi)
        @test size(ss.root_conductance_patch) == (np, nlevsoi)
        @test size(ss.soil_conductance_patch) == (np, nlevsoi)

        # --- Check root sizes (column-level 2D) ---
        @test size(ss.rootr_col) == (nc, nlevgrnd)
        @test size(ss.rootfr_col) == (nc, nlevgrnd)
        @test size(ss.rootr_road_perv_col) == (nc, nlevgrnd)
        @test size(ss.rootfr_road_perv_col) == (nc, nlevgrnd)

        # --- Check NaN initialization ---
        @test all(isnan, ss.sandfrac_patch)
        @test all(isnan, ss.clayfrac_patch)
        @test all(isnan, ss.mss_frc_cly_vld_col)
        @test all(isnan, ss.cellorg_col)
        @test all(isnan, ss.cellsand_col)
        @test all(isnan, ss.cellclay_col)
        @test all(isnan, ss.bd_col)
        @test all(isnan, ss.hk_l_col)
        @test all(isnan, ss.smp_l_col)
        @test all(isnan, ss.smpmin_col)
        @test all(isnan, ss.bsw_col)
        @test all(isnan, ss.watsat_col)
        @test all(isnan, ss.watfc_col)
        @test all(isnan, ss.soilbeta_col)
        @test all(isnan, ss.soilalpha_col)
        @test all(isnan, ss.soilalpha_u_col)
        @test all(isnan, ss.soilpsi_col)
        @test all(isnan, ss.wtfact_col)
        @test all(isnan, ss.gwc_thr_col)
        @test all(isnan, ss.thk_col)
        @test all(isnan, ss.tkmg_col)
        @test all(isnan, ss.tkdry_col)
        @test all(isnan, ss.tksatu_col)
        @test all(isnan, ss.csol_col)
        @test all(isnan, ss.rootr_patch)
        @test all(isnan, ss.rootfr_patch)
        @test all(isnan, ss.crootfr_patch)
        @test all(isnan, ss.root_depth_patch)
        @test all(isnan, ss.rootr_col)
        @test all(isnan, ss.rootfr_col)
        @test all(isnan, ss.msw_col)
        @test all(isnan, ss.nsw_col)
        @test all(isnan, ss.alphasw_col)
        @test all(isnan, ss.watres_col)
        @test all(isnan, ss.k_soil_root_patch)
        @test all(isnan, ss.root_conductance_patch)
        @test all(isnan, ss.soil_conductance_patch)

        # --- Check SPVAL initialization ---
        @test all(==(CLM.SPVAL), ss.hksat_col)
        @test all(==(CLM.SPVAL), ss.hksat_min_col)
        @test all(==(CLM.SPVAL), ss.watdry_col)
        @test all(==(CLM.SPVAL), ss.watopt_col)
        @test all(==(CLM.SPVAL), ss.sucsat_col)
        @test all(==(CLM.SPVAL), ss.dsl_col)
        @test all(==(CLM.SPVAL), ss.soilresis_col)
        @test all(==(CLM.SPVAL), ss.porosity_col)
        @test all(==(CLM.SPVAL), ss.eff_porosity_col)
    end

    @testset "soilstate_clean!" begin
        ss = CLM.SoilStateData()
        CLM.soilstate_init!(ss, 10, 5)
        CLM.soilstate_clean!(ss)

        # All vectors should be empty
        @test length(ss.sandfrac_patch) == 0
        @test length(ss.clayfrac_patch) == 0
        @test length(ss.mss_frc_cly_vld_col) == 0
        @test length(ss.smpmin_col) == 0
        @test length(ss.dsl_col) == 0
        @test length(ss.soilresis_col) == 0
        @test length(ss.soilbeta_col) == 0
        @test length(ss.soilalpha_col) == 0
        @test length(ss.soilalpha_u_col) == 0
        @test length(ss.wtfact_col) == 0
        @test length(ss.gwc_thr_col) == 0
        @test length(ss.root_depth_patch) == 0

        # All matrices should be empty
        @test size(ss.cellorg_col) == (0, 0)
        @test size(ss.hksat_col) == (0, 0)
        @test size(ss.watsat_col) == (0, 0)
        @test size(ss.thk_col) == (0, 0)
        @test size(ss.rootr_patch) == (0, 0)
        @test size(ss.rootfr_patch) == (0, 0)
        @test size(ss.crootfr_patch) == (0, 0)
        @test size(ss.msw_col) == (0, 0)
        @test size(ss.k_soil_root_patch) == (0, 0)
        @test size(ss.root_conductance_patch) == (0, 0)
        @test size(ss.soil_conductance_patch) == (0, 0)
        @test size(ss.porosity_col) == (0, 0)
        @test size(ss.eff_porosity_col) == (0, 0)
        @test size(ss.bd_col) == (0, 0)
        @test size(ss.smp_l_col) == (0, 0)
        @test size(ss.hk_l_col) == (0, 0)
    end

    @testset "soilstate_init_cold!" begin
        nc = 5
        ss = CLM.SoilStateData()
        CLM.soilstate_init!(ss, 10, nc)

        CLM.soilstate_init_cold!(ss, 1:nc)

        # smp_l_col should be -1000 for all columns and levels
        for c in 1:nc
            for j in 1:nlevgrnd
                @test ss.smp_l_col[c, j] ≈ -1000.0
            end
        end

        # hk_l_col should be 0 for all columns and levels
        for c in 1:nc
            for j in 1:nlevgrnd
                @test ss.hk_l_col[c, j] ≈ 0.0
            end
        end

        # Other fields should remain at their init values
        @test all(==(CLM.SPVAL), ss.hksat_col)
        @test all(isnan, ss.bsw_col)
    end

    @testset "field mutability" begin
        ss = CLM.SoilStateData()
        CLM.soilstate_init!(ss, 5, 3)

        # Write and verify patch-level fields
        ss.sandfrac_patch[1] = 0.5
        @test ss.sandfrac_patch[1] == 0.5

        ss.clayfrac_patch[2] = 0.3
        @test ss.clayfrac_patch[2] == 0.3

        ss.root_depth_patch[1] = 1.5
        @test ss.root_depth_patch[1] == 1.5

        # Write and verify column-level 1D fields
        ss.dsl_col[1] = 10.0
        @test ss.dsl_col[1] == 10.0

        ss.soilbeta_col[2] = 0.8
        @test ss.soilbeta_col[2] == 0.8

        # Write and verify column-level 2D fields
        ss.hksat_col[1, 1] = 0.01
        @test ss.hksat_col[1, 1] == 0.01

        ss.watsat_col[2, 3] = 0.45
        @test ss.watsat_col[2, 3] == 0.45

        ss.thk_col[1, 1] = 2.0
        @test ss.thk_col[1, 1] == 2.0

        ss.bsw_col[1, 2] = 5.5
        @test ss.bsw_col[1, 2] == 5.5

        # Write and verify patch-level 2D fields
        ss.rootr_patch[1, 1] = 0.1
        @test ss.rootr_patch[1, 1] == 0.1

        ss.rootfr_patch[2, 3] = 0.05
        @test ss.rootfr_patch[2, 3] == 0.05

        ss.k_soil_root_patch[1, 1] = 0.001
        @test ss.k_soil_root_patch[1, 1] == 0.001
    end

    @testset "re-init overwrites previous state" begin
        ss = CLM.SoilStateData()
        CLM.soilstate_init!(ss, 5, 3)
        ss.sandfrac_patch[1] = 999.0

        CLM.soilstate_init!(ss, 10, 7)
        @test length(ss.sandfrac_patch) == 10
        @test all(isnan, ss.sandfrac_patch)
        @test length(ss.dsl_col) == 7
        @test size(ss.hksat_col) == (7, nlevgrnd)
        @test size(ss.rootr_patch) == (10, nlevgrnd)
        @test size(ss.thk_col) == (7, nlev_soisno)
    end

    @testset "stub functions run without error" begin
        ss = CLM.SoilStateData()
        CLM.soilstate_init!(ss, 5, 3)

        @test CLM.soilstate_init_history!(ss, 1:3, 1:5) === nothing
        @test CLM.soilstate_restart!(ss, 1:3, 1:5) === nothing
    end

end
