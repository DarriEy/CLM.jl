@testset "TemperatureData" begin

    # Ensure varpar is initialized for dimension parameters
    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)

    nlevsno        = CLM.varpar.nlevsno
    nlevmaxurbgrnd = CLM.varpar.nlevmaxurbgrnd
    nlevlak        = CLM.varpar.nlevlak
    nlevgrnd       = CLM.varpar.nlevgrnd
    nlevurb        = CLM.varpar.nlevurb
    nlev_soisno    = nlevsno + nlevmaxurbgrnd

    @testset "default construction" begin
        temp = CLM.TemperatureData()
        @test length(temp.t_veg_patch) == 0
        @test length(temp.t_grnd_col) == 0
        @test size(temp.t_soisno_col) == (0, 0)
        @test size(temp.t_lake_col) == (0, 0)
        @test size(temp.imelt_col) == (0, 0)
        @test length(temp.t_building_lun) == 0
        @test length(temp.heat1_grc) == 0
        @test temp.excess_ice_coldstart_depth == 0.5
        @test temp.excess_ice_coldstart_temp == -5.0
    end

    @testset "temperature_init!" begin
        np = 20   # patches
        nc = 10   # columns
        nl = 3    # landunits
        ng = 2    # gridcells

        temp = CLM.TemperatureData()
        CLM.temperature_init!(temp, np, nc, nl, ng)

        # --- Check patch-level sizes ---
        @test length(temp.t_stem_patch) == np
        @test length(temp.t_veg_patch) == np
        @test length(temp.t_skin_patch) == np
        @test length(temp.t_veg_day_patch) == np
        @test length(temp.t_veg_night_patch) == np
        @test length(temp.t_veg10_day_patch) == np
        @test length(temp.t_veg10_night_patch) == np
        @test length(temp.ndaysteps_patch) == np
        @test length(temp.nnightsteps_patch) == np
        @test length(temp.dt_veg_patch) == np
        @test length(temp.thm_patch) == np
        @test length(temp.t_ref2m_patch) == np
        @test length(temp.t_ref2m_r_patch) == np
        @test length(temp.t_ref2m_u_patch) == np
        @test length(temp.t_ref2m_min_patch) == np
        @test length(temp.t_ref2m_min_r_patch) == np
        @test length(temp.t_ref2m_min_u_patch) == np
        @test length(temp.t_ref2m_max_patch) == np
        @test length(temp.t_ref2m_max_r_patch) == np
        @test length(temp.t_ref2m_max_u_patch) == np
        @test length(temp.t_ref2m_min_inst_patch) == np
        @test length(temp.t_ref2m_min_inst_r_patch) == np
        @test length(temp.t_ref2m_min_inst_u_patch) == np
        @test length(temp.t_ref2m_max_inst_patch) == np
        @test length(temp.t_ref2m_max_inst_r_patch) == np
        @test length(temp.t_ref2m_max_inst_u_patch) == np
        @test length(temp.t_a10_patch) == np
        @test length(temp.t_a10min_patch) == np
        @test length(temp.t_a5min_patch) == np
        @test length(temp.t_veg24_patch) == np
        @test length(temp.t_veg240_patch) == np
        @test length(temp.gdd0_patch) == np
        @test length(temp.gdd8_patch) == np
        @test length(temp.gdd10_patch) == np
        @test length(temp.gdd020_patch) == np
        @test length(temp.gdd820_patch) == np
        @test length(temp.gdd1020_patch) == np
        @test length(temp.emv_patch) == np

        # --- Check column-level sizes (1D) ---
        @test length(temp.t_h2osfc_col) == nc
        @test length(temp.t_h2osfc_bef_col) == nc
        @test length(temp.tsl_col) == nc
        @test length(temp.t_soi10cm_col) == nc
        @test length(temp.t_soi17cm_col) == nc
        @test length(temp.t_sno_mul_mss_col) == nc
        @test length(temp.t_grnd_col) == nc
        @test length(temp.t_grnd_r_col) == nc
        @test length(temp.t_grnd_u_col) == nc
        @test length(temp.snot_top_col) == nc
        @test length(temp.dTdz_top_col) == nc
        @test length(temp.dt_grnd_col) == nc
        @test length(temp.thv_col) == nc
        @test length(temp.soila10_col) == nc
        @test length(temp.beta_col) == nc
        @test length(temp.dynbal_baseline_heat_col) == nc
        @test length(temp.emg_col) == nc
        @test length(temp.xmf_col) == nc
        @test length(temp.xmf_h2osfc_col) == nc
        @test length(temp.c_h2osfc_col) == nc

        # --- Check column-level sizes (2D) ---
        @test size(temp.t_ssbef_col) == (nc, nlev_soisno)
        @test size(temp.t_soisno_col) == (nc, nlev_soisno)
        @test size(temp.t_lake_col) == (nc, nlevlak)
        @test size(temp.imelt_col) == (nc, nlev_soisno)
        @test size(temp.fact_col) == (nc, nlev_soisno)

        # --- Check landunit-level sizes ---
        @test length(temp.t_building_lun) == nl
        @test length(temp.t_roof_inner_lun) == nl
        @test length(temp.t_sunw_inner_lun) == nl
        @test length(temp.t_shdw_inner_lun) == nl
        @test length(temp.t_floor_lun) == nl
        @test length(temp.taf_lun) == nl

        # --- Check gridcell-level sizes ---
        @test length(temp.heat1_grc) == ng
        @test length(temp.heat2_grc) == ng
        @test length(temp.liquid_water_temp1_grc) == ng
        @test length(temp.liquid_water_temp2_grc) == ng

        # --- Check NaN initialization for standard real fields ---
        @test all(isnan, temp.t_stem_patch)
        @test all(isnan, temp.t_veg_patch)
        @test all(isnan, temp.t_skin_patch)
        @test all(isnan, temp.t_h2osfc_col)
        @test all(isnan, temp.t_grnd_col)
        @test all(isnan, temp.t_ssbef_col)
        @test all(isnan, temp.t_soisno_col)
        @test all(isnan, temp.t_lake_col)
        @test all(isnan, temp.t_building_lun)
        @test all(isnan, temp.heat1_grc)
        @test all(isnan, temp.emv_patch)
        @test all(isnan, temp.emg_col)
        @test all(isnan, temp.fact_col)

        # --- Check SPVAL initialization ---
        @test all(==(CLM.SPVAL), temp.t_veg_day_patch)
        @test all(==(CLM.SPVAL), temp.t_veg_night_patch)
        @test all(==(CLM.SPVAL), temp.t_veg10_day_patch)
        @test all(==(CLM.SPVAL), temp.t_veg10_night_patch)
        @test all(==(CLM.SPVAL), temp.t_soi17cm_col)
        @test all(==(CLM.SPVAL), temp.gdd0_patch)
        @test all(==(CLM.SPVAL), temp.gdd8_patch)
        @test all(==(CLM.SPVAL), temp.gdd10_patch)
        @test all(==(CLM.SPVAL), temp.gdd020_patch)
        @test all(==(CLM.SPVAL), temp.gdd820_patch)
        @test all(==(CLM.SPVAL), temp.gdd1020_patch)

        # --- Check ISPVAL initialization for integer fields ---
        @test all(==(CLM.ISPVAL), temp.ndaysteps_patch)
        @test all(==(CLM.ISPVAL), temp.nnightsteps_patch)

        # --- Check imelt_col initialized to typemax(Int) ---
        @test all(==(typemax(Int)), temp.imelt_col)
    end

    @testset "temperature_clean!" begin
        temp = CLM.TemperatureData()
        CLM.temperature_init!(temp, 10, 5, 2, 1)
        CLM.temperature_clean!(temp)

        # All vectors should be empty
        @test length(temp.t_veg_patch) == 0
        @test length(temp.t_grnd_col) == 0
        @test length(temp.t_building_lun) == 0
        @test length(temp.heat1_grc) == 0
        @test length(temp.gdd020_patch) == 0
        @test length(temp.emv_patch) == 0
        @test length(temp.emg_col) == 0
        @test length(temp.xmf_col) == 0
        @test length(temp.c_h2osfc_col) == 0
        @test size(temp.t_soisno_col) == (0, 0)
        @test size(temp.t_lake_col) == (0, 0)
        @test size(temp.imelt_col) == (0, 0)
        @test size(temp.fact_col) == (0, 0)
        @test size(temp.t_ssbef_col) == (0, 0)
    end

    @testset "temperature_init_cold!" begin
        # Setup: 1 gridcell, 2 landunits (soil + lake), 4 columns, 4 patches
        nc = 4; nl = 2; np = 4; ng = 1

        temp = CLM.TemperatureData()
        CLM.temperature_init!(temp, np, nc, nl, ng)

        col = CLM.ColumnData()
        CLM.column_init!(col, nc)

        lun = CLM.LandunitData()
        CLM.landunit_init!(lun, nl)

        pch = CLM.PatchData()
        CLM.patch_init!(pch, np)

        joff = nlevsno  # offset for snow layer indexing

        # Landunit 1: soil (non-lake, non-urban, non-special)
        lun.itype[1]    = CLM.ISTSOIL
        lun.lakpoi[1]   = false
        lun.urbpoi[1]   = false
        lun.ifspecial[1] = false
        lun.coli[1]     = 1
        lun.colf[1]     = 2

        # Landunit 2: lake
        lun.itype[2]    = CLM.ISTDLAK
        lun.lakpoi[2]   = true
        lun.urbpoi[2]   = false
        lun.ifspecial[2] = true
        lun.coli[2]     = 3
        lun.colf[2]     = 4

        # Columns 1,2 on landunit 1 (soil); columns 3,4 on landunit 2 (lake)
        col.landunit[1] = 1; col.landunit[2] = 1
        col.landunit[3] = 2; col.landunit[4] = 2
        col.itype[1] = CLM.ISTSOIL; col.itype[2] = CLM.ISTSOIL
        col.itype[3] = CLM.ISTDLAK; col.itype[4] = CLM.ISTDLAK
        col.snl[1] = 0; col.snl[2] = -2; col.snl[3] = 0; col.snl[4] = 0

        # Patches: 1,2 on landunit 1; 3,4 on landunit 2
        pch.column[1] = 1; pch.column[2] = 2
        pch.column[3] = 3; pch.column[4] = 4
        pch.landunit[1] = 1; pch.landunit[2] = 1
        pch.landunit[3] = 2; pch.landunit[4] = 2

        em_roof    = fill(0.9, nl)
        em_wall    = fill(0.85, nl)
        em_improad = fill(0.95, nl)
        em_perroad = fill(0.92, nl)

        CLM.temperature_init_cold!(temp, col, lun, pch,
            1:nc, 1:np, 1:nl;
            em_roof_lun    = em_roof,
            em_wall_lun    = em_wall,
            em_improad_lun = em_improad,
            em_perroad_lun = em_perroad,
            is_prog_buildtemp = false)

        # Soil columns should have 272 K for soil levels
        @test temp.t_soisno_col[1, 1 + joff] ≈ 272.0
        @test temp.t_soisno_col[1, nlevgrnd + joff] ≈ 272.0

        # Column 2 has snl=-2, so snow layers should be 250 K
        # Snow layers at Fortran indices -1 and 0 → Julia indices (joff-1) and joff
        @test temp.t_soisno_col[2, joff - 1] ≈ 250.0
        @test temp.t_soisno_col[2, joff]     ≈ 250.0
        # And soil layers at 272
        @test temp.t_soisno_col[2, 1 + joff] ≈ 272.0

        # Ground temp for soil column (snl=0) should equal top soil layer
        @test temp.t_grnd_col[1] ≈ temp.t_soisno_col[1, 0 + 1 + joff]

        # Ground temp for lake column should be 277
        @test temp.t_grnd_col[3] ≈ 277.0

        # Lake temperature should equal ground temperature
        @test temp.t_lake_col[3, 1] ≈ 277.0
        @test temp.t_lake_col[3, nlevlak] ≈ 277.0

        # Lake soil/snow should also be set to ground temperature
        @test temp.t_soisno_col[3, 1 + joff] ≈ 277.0

        # t_h2osfc should be 274
        @test temp.t_h2osfc_col[1] ≈ 274.0
        @test temp.t_h2osfc_col[3] ≈ 274.0

        # Vegetation temperature for all patches should be 283
        @test temp.t_veg_patch[1] ≈ 283.0
        @test temp.t_veg_patch[3] ≈ 283.0

        # Stem temp = veg temp
        @test temp.t_stem_patch[1] ≈ 283.0

        # t_ref2m for all patches should be 283
        @test temp.t_ref2m_patch[1] ≈ 283.0

        # Non-urban, non-special: t_ref2m_r should be 283
        @test temp.t_ref2m_r_patch[1] ≈ 283.0

        # Lake (ifspecial=true): t_ref2m_r should be SPVAL
        @test temp.t_ref2m_r_patch[3] ≈ CLM.SPVAL

        # t_soi17cm should equal t_grnd
        @test temp.t_soi17cm_col[1] ≈ temp.t_grnd_col[1]

        # dynbal_baseline_heat should be 0
        @test temp.dynbal_baseline_heat_col[1] ≈ 0.0
        @test temp.dynbal_baseline_heat_col[4] ≈ 0.0
    end

    @testset "temperature_init_acc_vars! startup" begin
        np = 5; nc = 3; nl = 1; ng = 1
        temp = CLM.TemperatureData()
        CLM.temperature_init!(temp, np, nc, nl, ng)

        CLM.temperature_init_acc_vars!(temp, 1:nc, 1:np; is_startup=true)

        # After startup init, max_inst should be -SPVAL, min_inst should be SPVAL
        @test temp.t_ref2m_max_inst_patch[1] ≈ -CLM.SPVAL
        @test temp.t_ref2m_min_inst_patch[1] ≈  CLM.SPVAL
        @test temp.t_ref2m_max_inst_r_patch[1] ≈ -CLM.SPVAL
        @test temp.t_ref2m_min_inst_r_patch[1] ≈  CLM.SPVAL
        @test temp.t_ref2m_max_inst_u_patch[1] ≈ -CLM.SPVAL
        @test temp.t_ref2m_min_inst_u_patch[1] ≈  CLM.SPVAL

        # max/min daily should be SPVAL
        @test temp.t_ref2m_max_patch[1] ≈ CLM.SPVAL
        @test temp.t_ref2m_min_patch[1] ≈ CLM.SPVAL
    end

    @testset "temperature_read_namelist!" begin
        temp = CLM.TemperatureData()

        CLM.temperature_read_namelist!(temp;
            excess_ice_coldstart_depth = 1.0,
            excess_ice_coldstart_temp  = -10.0)
        @test temp.excess_ice_coldstart_depth ≈ 1.0
        @test temp.excess_ice_coldstart_temp ≈ -10.0

        # Invalid depth should error
        @test_throws ErrorException CLM.temperature_read_namelist!(temp;
            excess_ice_coldstart_depth = -0.5,
            excess_ice_coldstart_temp  = -5.0)

        # Invalid temp should error
        @test_throws ErrorException CLM.temperature_read_namelist!(temp;
            excess_ice_coldstart_depth = 0.5,
            excess_ice_coldstart_temp  = 1.0)
    end

    @testset "field mutability" begin
        temp = CLM.TemperatureData()
        CLM.temperature_init!(temp, 5, 3, 2, 1)

        # Write and verify patch-level fields
        temp.t_veg_patch[1] = 300.0
        @test temp.t_veg_patch[1] == 300.0

        temp.t_stem_patch[2] = 295.0
        @test temp.t_stem_patch[2] == 295.0

        # Write and verify column-level fields
        temp.t_grnd_col[1] = 280.0
        @test temp.t_grnd_col[1] == 280.0

        joff = nlevsno
        temp.t_soisno_col[1, 1 + joff] = 273.15
        @test temp.t_soisno_col[1, 1 + joff] == 273.15

        temp.imelt_col[1, 1 + joff] = 1
        @test temp.imelt_col[1, 1 + joff] == 1

        # Write and verify landunit-level fields
        temp.t_building_lun[1] = 293.0
        @test temp.t_building_lun[1] == 293.0

        # Write and verify gridcell-level fields
        temp.heat1_grc[1] = 1.0e8
        @test temp.heat1_grc[1] == 1.0e8
    end

    @testset "re-init overwrites previous state" begin
        temp = CLM.TemperatureData()
        CLM.temperature_init!(temp, 5, 3, 2, 1)
        temp.t_veg_patch[1] = 999.0

        CLM.temperature_init!(temp, 10, 7, 3, 2)
        @test length(temp.t_veg_patch) == 10
        @test all(isnan, temp.t_veg_patch)
        @test length(temp.t_grnd_col) == 7
        @test length(temp.t_building_lun) == 3
        @test length(temp.heat1_grc) == 2
    end

    @testset "stub functions run without error" begin
        temp = CLM.TemperatureData()
        CLM.temperature_init!(temp, 5, 3, 2, 1)

        # These should run without error (they are stubs)
        @test CLM.temperature_init_history!(temp, 1:3, 1:5, 1:2, 1:1) === nothing
        @test CLM.temperature_restart!(temp, 1:3, 1:5, 1:2) === nothing
        @test CLM.temperature_init_acc_buffer!(temp, 1:5) === nothing
        @test CLM.temperature_update_acc_vars!(temp, 1:3, 1:5,
            CLM.LandunitData(), CLM.PatchData()) === nothing
        @test CLM.temperature_update_acc_vars_crop_gdds!(temp, 1:5, fill(0.0, 5)) === nothing
    end

end
