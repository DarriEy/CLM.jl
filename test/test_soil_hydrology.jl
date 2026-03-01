@testset "SoilHydrologyData" begin

    # Ensure varpar is initialized for dimension parameters
    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)

    nlevgrnd = CLM.varpar.nlevgrnd
    nlevsoi  = CLM.varpar.nlevsoi
    nlayer   = CLM.NLAYER
    nlayert  = CLM.varpar.nlayert

    @testset "default construction" begin
        sh = CLM.SoilHydrologyData()
        @test sh.h2osfcflag == 1
        @test length(sh.num_substeps_col) == 0
        @test length(sh.frost_table_col) == 0
        @test length(sh.zwt_col) == 0
        @test length(sh.zwts_col) == 0
        @test length(sh.zwt_perched_col) == 0
        @test length(sh.qcharge_col) == 0
        @test length(sh.h2osfc_thresh_col) == 0
        @test length(sh.xs_urban_col) == 0
        @test length(sh.hkdepth_col) == 0
        @test length(sh.b_infil_col) == 0
        @test length(sh.ds_col) == 0
        @test length(sh.dsmax_col) == 0
        @test length(sh.Wsvic_col) == 0
        @test length(sh.c_param_col) == 0
        @test length(sh.top_moist_col) == 0
        @test length(sh.top_max_moist_col) == 0
        @test length(sh.top_ice_col) == 0
        @test length(sh.top_moist_limited_col) == 0
        @test size(sh.icefrac_col) == (0, 0)
        @test size(sh.porosity_col) == (0, 0)
        @test size(sh.depth_col) == (0, 0)
        @test size(sh.expt_col) == (0, 0)
        @test size(sh.ksat_col) == (0, 0)
        @test size(sh.phi_s_col) == (0, 0)
        @test size(sh.moist_col) == (0, 0)
        @test size(sh.moist_vol_col) == (0, 0)
        @test size(sh.max_moist_col) == (0, 0)
        @test size(sh.ice_col) == (0, 0)
        @test size(sh.vic_clm_fract_col) == (0, 0, 0)
    end

    @testset "soilhydrology_init!" begin
        nc = 10
        sh = CLM.SoilHydrologyData()
        CLM.soilhydrology_init!(sh, nc)

        # --- Check 1D sizes ---
        @test length(sh.num_substeps_col) == nc
        @test length(sh.frost_table_col) == nc
        @test length(sh.zwt_col) == nc
        @test length(sh.zwts_col) == nc
        @test length(sh.zwt_perched_col) == nc
        @test length(sh.qcharge_col) == nc
        @test length(sh.h2osfc_thresh_col) == nc
        @test length(sh.xs_urban_col) == nc
        @test length(sh.hkdepth_col) == nc
        @test length(sh.b_infil_col) == nc
        @test length(sh.ds_col) == nc
        @test length(sh.dsmax_col) == nc
        @test length(sh.Wsvic_col) == nc
        @test length(sh.c_param_col) == nc
        @test length(sh.top_moist_col) == nc
        @test length(sh.top_max_moist_col) == nc
        @test length(sh.top_ice_col) == nc
        @test length(sh.top_moist_limited_col) == nc

        # --- Check 2D sizes ---
        @test size(sh.icefrac_col) == (nc, nlevgrnd)
        @test size(sh.porosity_col) == (nc, nlayer)
        @test size(sh.depth_col) == (nc, nlayert)
        @test size(sh.expt_col) == (nc, nlayer)
        @test size(sh.ksat_col) == (nc, nlayer)
        @test size(sh.phi_s_col) == (nc, nlayer)
        @test size(sh.moist_col) == (nc, nlayert)
        @test size(sh.moist_vol_col) == (nc, nlayert)
        @test size(sh.max_moist_col) == (nc, nlayer)
        @test size(sh.ice_col) == (nc, nlayert)

        # --- Check 3D size ---
        @test size(sh.vic_clm_fract_col) == (nc, nlayer, nlevsoi)

        # --- Check NaN initialization ---
        @test all(isnan, sh.num_substeps_col)
        @test all(isnan, sh.frost_table_col)
        @test all(isnan, sh.zwt_col)
        @test all(isnan, sh.zwts_col)
        @test all(isnan, sh.zwt_perched_col)
        @test all(isnan, sh.qcharge_col)
        @test all(isnan, sh.h2osfc_thresh_col)
        @test all(isnan, sh.xs_urban_col)
        @test all(isnan, sh.icefrac_col)
        @test all(isnan, sh.hkdepth_col)
        @test all(isnan, sh.b_infil_col)
        @test all(isnan, sh.ds_col)
        @test all(isnan, sh.dsmax_col)
        @test all(isnan, sh.Wsvic_col)
        @test all(isnan, sh.c_param_col)
        @test all(isnan, sh.top_moist_col)
        @test all(isnan, sh.top_max_moist_col)
        @test all(isnan, sh.top_ice_col)
        @test all(isnan, sh.top_moist_limited_col)
        @test all(isnan, sh.porosity_col)
        @test all(isnan, sh.depth_col)
        @test all(isnan, sh.expt_col)
        @test all(isnan, sh.ksat_col)
        @test all(isnan, sh.phi_s_col)
        @test all(isnan, sh.moist_col)
        @test all(isnan, sh.moist_vol_col)
        @test all(isnan, sh.max_moist_col)
        @test all(isnan, sh.ice_col)
        @test all(isnan, sh.vic_clm_fract_col)
    end

    @testset "soilhydrology_clean!" begin
        sh = CLM.SoilHydrologyData()
        CLM.soilhydrology_init!(sh, 5)
        CLM.soilhydrology_clean!(sh)

        # All vectors should be empty
        @test length(sh.num_substeps_col) == 0
        @test length(sh.frost_table_col) == 0
        @test length(sh.zwt_col) == 0
        @test length(sh.zwts_col) == 0
        @test length(sh.zwt_perched_col) == 0
        @test length(sh.qcharge_col) == 0
        @test length(sh.h2osfc_thresh_col) == 0
        @test length(sh.xs_urban_col) == 0
        @test length(sh.hkdepth_col) == 0
        @test length(sh.b_infil_col) == 0
        @test length(sh.ds_col) == 0
        @test length(sh.dsmax_col) == 0
        @test length(sh.Wsvic_col) == 0
        @test length(sh.c_param_col) == 0
        @test length(sh.top_moist_col) == 0
        @test length(sh.top_max_moist_col) == 0
        @test length(sh.top_ice_col) == 0
        @test length(sh.top_moist_limited_col) == 0

        # All matrices should be empty
        @test size(sh.icefrac_col) == (0, 0)
        @test size(sh.porosity_col) == (0, 0)
        @test size(sh.depth_col) == (0, 0)
        @test size(sh.expt_col) == (0, 0)
        @test size(sh.ksat_col) == (0, 0)
        @test size(sh.phi_s_col) == (0, 0)
        @test size(sh.moist_col) == (0, 0)
        @test size(sh.moist_vol_col) == (0, 0)
        @test size(sh.max_moist_col) == (0, 0)
        @test size(sh.ice_col) == (0, 0)

        # 3D should be empty
        @test size(sh.vic_clm_fract_col) == (0, 0, 0)
    end

    @testset "soilhydrology_init_cold! (basic)" begin
        nc = 5
        sh = CLM.SoilHydrologyData()
        CLM.soilhydrology_init!(sh, nc)

        CLM.soilhydrology_init_cold!(sh, 1:nc)

        # num_substeps_col should be SPVAL
        for c in 1:nc
            @test sh.num_substeps_col[c] == CLM.SPVAL
        end

        # zwt_col should be 0
        for c in 1:nc
            @test sh.zwt_col[c] ≈ 0.0
        end

        # zwt_perched_col and frost_table_col should be 0 (basic init without column metadata)
        for c in 1:nc
            @test sh.zwt_perched_col[c] ≈ 0.0
            @test sh.frost_table_col[c] ≈ 0.0
        end

        # Other fields should remain at NaN
        @test all(isnan, sh.icefrac_col)
        @test all(isnan, sh.qcharge_col)
    end

    @testset "soilhydrology_init_cold! (with column metadata)" begin
        nc = 3
        sh = CLM.SoilHydrologyData()
        CLM.soilhydrology_init!(sh, nc)

        # Mock column metadata: all non-urban, non-lake columns
        zi_col = fill(3.8, nc, nlevsoi)  # interface depths
        wa_col = fill(4800.0, nc)         # aquifer water
        landunit_col = [1, 1, 1]
        lakpoi = falses(1)                # landunit 1 is not lake
        urbpoi = falses(1)                # landunit 1 is not urban
        itype_col = [1, 1, 1]
        nbedrock_col = [nlevsoi, nlevsoi, nlevsoi]

        CLM.soilhydrology_init_cold!(sh, 1:nc;
            use_aquifer_layer=true,
            zi_col=zi_col,
            wa_col=wa_col,
            nlevsoi=nlevsoi,
            landunit_col=landunit_col,
            lakpoi=lakpoi,
            urbpoi=urbpoi,
            itype_col=itype_col,
            nbedrock_col=nbedrock_col)

        # zwt_col should be computed: (25 + zi[c,nlevsoi]) - wa[c]/0.2/1000
        for c in 1:nc
            expected_zwt = (25.0 + zi_col[c, nlevsoi]) - wa_col[c] / 0.2 / 1000.0
            @test sh.zwt_col[c] ≈ expected_zwt
        end

        # Non-urban, non-lake: zwt_perched and frost_table = zi[c, nlevsoi]
        for c in 1:nc
            @test sh.zwt_perched_col[c] ≈ zi_col[c, nlevsoi]
            @test sh.frost_table_col[c] ≈ zi_col[c, nlevsoi]
        end
    end

    @testset "soilhydrology_read_nl!" begin
        sh = CLM.SoilHydrologyData()
        @test sh.h2osfcflag == 1

        CLM.soilhydrology_read_nl!(sh; h2osfcflag=0)
        @test sh.h2osfcflag == 0

        CLM.soilhydrology_read_nl!(sh; h2osfcflag=1)
        @test sh.h2osfcflag == 1
    end

    @testset "field mutability" begin
        sh = CLM.SoilHydrologyData()
        CLM.soilhydrology_init!(sh, 3)

        # Write and verify 1D fields
        sh.zwt_col[1] = 5.0
        @test sh.zwt_col[1] == 5.0

        sh.qcharge_col[2] = 0.001
        @test sh.qcharge_col[2] == 0.001

        sh.frost_table_col[3] = 2.5
        @test sh.frost_table_col[3] == 2.5

        sh.top_moist_col[1] = 100.0
        @test sh.top_moist_col[1] == 100.0

        # Write and verify 2D fields
        sh.icefrac_col[1, 1] = 0.5
        @test sh.icefrac_col[1, 1] == 0.5

        sh.porosity_col[2, 2] = 0.4
        @test sh.porosity_col[2, 2] == 0.4

        sh.moist_col[1, 3] = 50.0
        @test sh.moist_col[1, 3] == 50.0

        # Write and verify 3D field
        sh.vic_clm_fract_col[1, 1, 1] = 0.33
        @test sh.vic_clm_fract_col[1, 1, 1] == 0.33
    end

    @testset "re-init overwrites previous state" begin
        sh = CLM.SoilHydrologyData()
        CLM.soilhydrology_init!(sh, 3)
        sh.zwt_col[1] = 999.0

        CLM.soilhydrology_init!(sh, 7)
        @test length(sh.zwt_col) == 7
        @test all(isnan, sh.zwt_col)
        @test size(sh.icefrac_col) == (7, nlevgrnd)
        @test size(sh.porosity_col) == (7, nlayer)
        @test size(sh.vic_clm_fract_col) == (7, nlayer, nlevsoi)
    end

    @testset "stub functions run without error" begin
        sh = CLM.SoilHydrologyData()
        CLM.soilhydrology_init!(sh, 5)

        @test CLM.soilhydrology_init_history!(sh, 1:5) === nothing
        @test CLM.soilhydrology_restart!(sh, 1:5) === nothing
    end

end
