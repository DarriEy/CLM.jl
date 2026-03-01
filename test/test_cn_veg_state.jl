@testset "CNVegStateData" begin

    @testset "default construction" begin
        vs = CLM.CNVegStateData()
        # Integer vectors
        @test length(vs.burndate_patch) == 0
        @test length(vs.peaklai_patch) == 0
        @test length(vs.idop_patch) == 0
        @test length(vs.iyop_patch) == 0
        # Float64 vectors (sample)
        @test length(vs.dwt_smoothed_patch) == 0
        @test length(vs.hdidx_patch) == 0
        @test length(vs.dormant_flag_patch) == 0
        @test length(vs.lgdp_col) == 0
        @test length(vs.nfire_col) == 0
        @test length(vs.plantCN_patch) == 0
        # Matrices
        @test size(vs.gddmaturity_thisyr) == (0, 0)
        @test size(vs.arepr_patch) == (0, 0)
        @test size(vs.arepr_n_patch) == (0, 0)
    end

    @testset "cnveg_state_init!" begin
        np = 10
        nc = 5
        vs = CLM.CNVegStateData()
        CLM.cnveg_state_init!(vs, np, nc)

        # --- Check integer fields ---
        @test length(vs.burndate_patch) == np
        @test all(x -> x == CLM.ISPVAL, vs.burndate_patch)

        @test length(vs.peaklai_patch) == np
        @test all(x -> x == 0, vs.peaklai_patch)

        @test length(vs.idop_patch) == np
        @test all(x -> x == typemax(Int), vs.idop_patch)

        @test length(vs.iyop_patch) == np
        @test all(x -> x == CLM.ISPVAL, vs.iyop_patch)

        # --- Check patch-level float sizes and initial values ---
        @test length(vs.dwt_smoothed_patch) == np
        @test all(isnan, vs.dwt_smoothed_patch)

        @test length(vs.hdidx_patch) == np
        @test all(isnan, vs.hdidx_patch)

        @test length(vs.cumvd_patch) == np
        @test all(isnan, vs.cumvd_patch)

        @test length(vs.gddmaturity_patch) == np
        @test all(x -> x == CLM.SPVAL, vs.gddmaturity_patch)

        @test length(vs.huileaf_patch) == np
        @test all(isnan, vs.huileaf_patch)

        @test length(vs.huigrain_patch) == np
        @test all(x -> x ≈ 0.0, vs.huigrain_patch)

        @test length(vs.aleafi_patch) == np
        @test all(isnan, vs.aleafi_patch)

        @test length(vs.aleaf_patch) == np
        @test all(isnan, vs.aleaf_patch)

        @test length(vs.astem_patch) == np
        @test all(isnan, vs.astem_patch)

        @test length(vs.aroot_patch) == np
        @test all(isnan, vs.aroot_patch)

        @test length(vs.htmx_patch) == np
        @test all(x -> x ≈ 0.0, vs.htmx_patch)

        # --- Check 2D fields ---
        @test size(vs.gddmaturity_thisyr) == (np, CLM.MXHARVESTS)
        @test all(x -> x == CLM.SPVAL, vs.gddmaturity_thisyr)

        @test size(vs.arepr_patch) == (np, CLM.NREPR)
        @test all(isnan, vs.arepr_patch)

        # --- AgSys fields should be empty by default ---
        @test length(vs.aleaf_n_patch) == 0
        @test length(vs.astem_n_patch) == 0
        @test length(vs.aroot_n_patch) == 0
        @test size(vs.arepr_n_patch) == (0, 0)

        # --- Check column-level fields ---
        @test length(vs.lgdp_col) == nc
        @test length(vs.lgdp1_col) == nc
        @test length(vs.lpop_col) == nc

        @test length(vs.annavg_t2m_col) == nc
        @test all(isnan, vs.annavg_t2m_col)

        @test length(vs.annsum_counter_col) == nc
        @test all(isnan, vs.annsum_counter_col)

        @test length(vs.nfire_col) == nc
        @test all(x -> x == CLM.SPVAL, vs.nfire_col)

        @test length(vs.fsr_col) == nc
        @test all(isnan, vs.fsr_col)

        @test length(vs.lfc_col) == nc
        @test all(x -> x == CLM.SPVAL, vs.lfc_col)

        @test length(vs.lfc2_col) == nc
        @test all(x -> x ≈ 0.0, vs.lfc2_col)

        @test length(vs.dtrotr_col) == nc
        @test all(x -> x ≈ 0.0, vs.dtrotr_col)

        @test length(vs.farea_burned_col) == nc
        @test all(isnan, vs.farea_burned_col)

        # --- Check phenology fields ---
        @test length(vs.dormant_flag_patch) == np
        @test all(isnan, vs.dormant_flag_patch)

        @test length(vs.onset_flag_patch) == np
        @test all(isnan, vs.onset_flag_patch)

        @test length(vs.offset_flag_patch) == np
        @test all(isnan, vs.offset_flag_patch)

        @test length(vs.grain_flag_patch) == np
        @test all(isnan, vs.grain_flag_patch)

        # --- Check GPP accumulators ---
        @test length(vs.c_allometry_patch) == np
        @test all(isnan, vs.c_allometry_patch)

        @test length(vs.downreg_patch) == np
        @test all(isnan, vs.downreg_patch)

        @test length(vs.leafcn_offset_patch) == np
        @test all(isnan, vs.leafcn_offset_patch)

        @test length(vs.plantCN_patch) == np
        @test all(isnan, vs.plantCN_patch)
    end

    @testset "cnveg_state_init! with use_crop_agsys" begin
        np = 5
        nc = 3
        vs = CLM.CNVegStateData()
        CLM.cnveg_state_init!(vs, np, nc; use_crop_agsys=true)

        @test length(vs.aleaf_n_patch) == np
        @test all(isnan, vs.aleaf_n_patch)

        @test length(vs.astem_n_patch) == np
        @test all(isnan, vs.astem_n_patch)

        @test length(vs.aroot_n_patch) == np
        @test all(isnan, vs.aroot_n_patch)

        @test size(vs.arepr_n_patch) == (np, CLM.NREPR)
        @test all(isnan, vs.arepr_n_patch)
    end

    @testset "cnveg_state_clean!" begin
        vs = CLM.CNVegStateData()
        CLM.cnveg_state_init!(vs, 8, 4)
        CLM.cnveg_state_clean!(vs)

        # Integer vectors
        @test length(vs.burndate_patch) == 0
        @test length(vs.peaklai_patch) == 0
        @test length(vs.idop_patch) == 0
        @test length(vs.iyop_patch) == 0

        # Sample float vectors
        @test length(vs.dwt_smoothed_patch) == 0
        @test length(vs.dormant_flag_patch) == 0
        @test length(vs.plantCN_patch) == 0
        @test length(vs.lgdp_col) == 0
        @test length(vs.nfire_col) == 0
        @test length(vs.farea_burned_col) == 0

        # Matrices
        @test size(vs.gddmaturity_thisyr) == (0, 0)
        @test size(vs.arepr_patch) == (0, 0)
        @test size(vs.arepr_n_patch) == (0, 0)
    end

    @testset "cnveg_state_init_cold! (basic, no landunit info)" begin
        np = 6
        nc = 3
        vs = CLM.CNVegStateData()
        CLM.cnveg_state_init!(vs, np, nc)

        CLM.cnveg_state_init_cold!(vs, 1:nc, 1:np)

        # Column-level: soil/crop defaults
        for c in 1:nc
            @test vs.annsum_counter_col[c] ≈ 0.0
            @test vs.annavg_t2m_col[c] ≈ 280.0
            @test vs.baf_crop_col[c] ≈ 0.0
            @test vs.baf_peatf_col[c] ≈ 0.0
            @test vs.fbac_col[c] ≈ 0.0
            @test vs.fbac1_col[c] ≈ 0.0
            @test vs.farea_burned_col[c] ≈ 0.0
            @test vs.nfire_col[c] ≈ 0.0
            @test vs.lfc2_col[c] ≈ 0.0
        end

        # Patch-level: phenology
        for p in 1:np
            @test vs.dormant_flag_patch[p] ≈ 1.0
            @test vs.days_active_patch[p] ≈ 0.0
            @test vs.onset_flag_patch[p] ≈ 0.0
            @test vs.onset_counter_patch[p] ≈ 0.0
            @test vs.onset_gddflag_patch[p] ≈ 0.0
            @test vs.onset_fdd_patch[p] ≈ 0.0
            @test vs.onset_gdd_patch[p] ≈ 0.0
            @test vs.onset_swi_patch[p] ≈ 0.0
            @test vs.offset_flag_patch[p] ≈ 0.0
            @test vs.offset_counter_patch[p] ≈ 0.0
            @test vs.offset_fdd_patch[p] ≈ 0.0
            @test vs.offset_swi_patch[p] ≈ 0.0
            @test vs.lgsf_patch[p] ≈ 0.0
            @test vs.bglfr_patch[p] ≈ 0.0
            @test vs.bgtr_patch[p] ≈ 0.0
            @test vs.annavg_t2m_patch[p] ≈ 280.0
            @test vs.tempavg_t2m_patch[p] ≈ 0.0
            @test vs.grain_flag_patch[p] ≈ 0.0

            # non-phenology
            @test vs.c_allometry_patch[p] ≈ 0.0
            @test vs.n_allometry_patch[p] ≈ 0.0
            @test vs.tempsum_potential_gpp_patch[p] ≈ 0.0
            @test vs.annsum_potential_gpp_patch[p] ≈ 0.0
            @test vs.tempmax_retransn_patch[p] ≈ 0.0
            @test vs.annmax_retransn_patch[p] ≈ 0.0
            @test vs.downreg_patch[p] ≈ 0.0
            @test vs.leafcn_offset_patch[p] == CLM.SPVAL
            @test vs.plantCN_patch[p] == CLM.SPVAL
        end
    end

    @testset "cnveg_state_init_cold! (with landunit info)" begin
        np = 4
        nc = 2
        vs = CLM.CNVegStateData()
        CLM.cnveg_state_init!(vs, np, nc)

        # Landunit 1 = soil, Landunit 2 = ice (special)
        col_landunit = [1, 2]
        patch_landunit = [1, 1, 2, 2]
        lun_ifspecial = [false, true]
        lun_itype = [CLM.ISTSOIL, CLM.ISTICE]

        CLM.cnveg_state_init_cold!(vs, 1:nc, 1:np;
            col_landunit=col_landunit,
            patch_landunit=patch_landunit,
            lun_ifspecial=lun_ifspecial,
            lun_itype=lun_itype)

        # Column 1 (soil): should have soil/crop defaults
        @test vs.annsum_counter_col[1] ≈ 0.0
        @test vs.annavg_t2m_col[1] ≈ 280.0
        @test vs.nfire_col[1] ≈ 0.0

        # Column 2 (ice/special): should have SPVAL
        @test vs.annsum_counter_col[2] == CLM.SPVAL
        @test vs.annavg_t2m_col[2] == CLM.SPVAL
        @test vs.nfire_col[2] == CLM.SPVAL

        # Patches 1-2 (soil): should have soil defaults
        @test vs.dormant_flag_patch[1] ≈ 1.0
        @test vs.dormant_flag_patch[2] ≈ 1.0
        @test vs.annavg_t2m_patch[1] ≈ 280.0
        @test vs.downreg_patch[1] ≈ 0.0

        # Patches 3-4 (ice/special): should have SPVAL
        @test vs.dormant_flag_patch[3] == CLM.SPVAL
        @test vs.dormant_flag_patch[4] == CLM.SPVAL
        @test vs.annavg_t2m_patch[3] == CLM.SPVAL
        @test vs.downreg_patch[3] == CLM.SPVAL
    end

    @testset "field mutability" begin
        vs = CLM.CNVegStateData()
        CLM.cnveg_state_init!(vs, 3, 2)

        # Write and verify 1D fields
        vs.burndate_patch[1] = 100
        @test vs.burndate_patch[1] == 100

        vs.dormant_flag_patch[2] = 1.0
        @test vs.dormant_flag_patch[2] == 1.0

        vs.nfire_col[1] = 0.5
        @test vs.nfire_col[1] == 0.5

        vs.plantCN_patch[3] = 25.0
        @test vs.plantCN_patch[3] == 25.0

        # Write and verify 2D fields
        vs.gddmaturity_thisyr[1, 1] = 1500.0
        @test vs.gddmaturity_thisyr[1, 1] == 1500.0

        vs.arepr_patch[2, 1] = 0.3
        @test vs.arepr_patch[2, 1] == 0.3
    end

    @testset "re-init overwrites previous state" begin
        vs = CLM.CNVegStateData()
        CLM.cnveg_state_init!(vs, 3, 2)
        vs.dormant_flag_patch[1] = 999.0

        CLM.cnveg_state_init!(vs, 7, 4)
        @test length(vs.dormant_flag_patch) == 7
        @test all(isnan, vs.dormant_flag_patch)
        @test length(vs.nfire_col) == 4
        @test size(vs.gddmaturity_thisyr) == (7, CLM.MXHARVESTS)
        @test size(vs.arepr_patch) == (7, CLM.NREPR)
    end

    @testset "stub functions run without error" begin
        vs = CLM.CNVegStateData()
        CLM.cnveg_state_init!(vs, 5, 3)

        @test CLM.cnveg_state_init_history!(vs, 1:5, 1:3) === nothing
        @test CLM.cnveg_state_restart!(vs, 1:5, 1:3) === nothing
    end

end
