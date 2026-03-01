@testset "EnergyFluxData" begin

    # Ensure varpar is initialized for dimension parameters
    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)

    nlevgrnd = CLM.varpar.nlevgrnd

    @testset "default construction" begin
        ef = CLM.EnergyFluxData()
        @test length(ef.eflx_sh_tot_patch) == 0
        @test length(ef.eflx_lh_tot_patch) == 0
        @test length(ef.eflx_snomelt_col) == 0
        @test length(ef.eflx_dynbal_grc) == 0
        @test length(ef.eflx_building_lun) == 0
        @test size(ef.eflx_fgr_col) == (0, 0)
        @test size(ef.rresis_patch) == (0, 0)
        @test length(ef.btran_patch) == 0
        @test length(ef.errsoi_col) == 0
    end

    @testset "energyflux_init!" begin
        np = 20   # patches
        nc = 10   # columns
        nl = 3    # landunits
        ng = 2    # gridcells

        ef = CLM.EnergyFluxData()
        CLM.energyflux_init!(ef, np, nc, nl, ng)

        # --- Check patch-level sizes (1D) ---
        @test length(ef.eflx_sh_stem_patch) == np
        @test length(ef.eflx_sh_grnd_patch) == np
        @test length(ef.eflx_sh_veg_patch) == np
        @test length(ef.eflx_sh_snow_patch) == np
        @test length(ef.eflx_sh_soil_patch) == np
        @test length(ef.eflx_sh_h2osfc_patch) == np
        @test length(ef.eflx_sh_tot_patch) == np
        @test length(ef.eflx_sh_tot_u_patch) == np
        @test length(ef.eflx_sh_tot_r_patch) == np
        @test length(ef.eflx_lh_tot_patch) == np
        @test length(ef.eflx_lh_tot_u_patch) == np
        @test length(ef.eflx_lh_tot_r_patch) == np
        @test length(ef.eflx_lh_vegt_patch) == np
        @test length(ef.eflx_lh_vege_patch) == np
        @test length(ef.eflx_lh_grnd_patch) == np
        @test length(ef.eflx_soil_grnd_patch) == np
        @test length(ef.eflx_soil_grnd_u_patch) == np
        @test length(ef.eflx_soil_grnd_r_patch) == np
        @test length(ef.eflx_lwrad_net_patch) == np
        @test length(ef.eflx_lwrad_net_r_patch) == np
        @test length(ef.eflx_lwrad_net_u_patch) == np
        @test length(ef.eflx_lwrad_out_patch) == np
        @test length(ef.eflx_lwrad_out_r_patch) == np
        @test length(ef.eflx_lwrad_out_u_patch) == np
        @test length(ef.eflx_gnet_patch) == np
        @test length(ef.eflx_grnd_lake_patch) == np
        @test length(ef.eflx_anthro_patch) == np
        @test length(ef.eflx_traffic_patch) == np
        @test length(ef.eflx_wasteheat_patch) == np
        @test length(ef.eflx_ventilation_patch) == np
        @test length(ef.eflx_heat_from_ac_patch) == np
        @test length(ef.dgnetdT_patch) == np
        @test length(ef.netrad_patch) == np
        @test length(ef.cgrnd_patch) == np
        @test length(ef.cgrndl_patch) == np
        @test length(ef.cgrnds_patch) == np
        @test length(ef.dlrad_patch) == np
        @test length(ef.ulrad_patch) == np
        @test length(ef.taux_patch) == np
        @test length(ef.tauy_patch) == np
        @test length(ef.canopy_cond_patch) == np
        @test length(ef.btran_patch) == np
        @test length(ef.btran_min_patch) == np
        @test length(ef.btran_min_inst_patch) == np
        @test length(ef.bsun_patch) == np
        @test length(ef.bsha_patch) == np
        @test length(ef.dhsdt_canopy_patch) == np
        @test length(ef.errsoi_patch) == np
        @test length(ef.errseb_patch) == np
        @test length(ef.errsol_patch) == np
        @test length(ef.errlon_patch) == np

        # --- Check column-level sizes (1D) ---
        @test length(ef.eflx_sh_precip_conversion_col) == nc
        @test length(ef.eflx_snomelt_col) == nc
        @test length(ef.eflx_snomelt_r_col) == nc
        @test length(ef.eflx_snomelt_u_col) == nc
        @test length(ef.eflx_h2osfc_to_snow_col) == nc
        @test length(ef.eflx_bot_col) == nc
        @test length(ef.eflx_fgr12_col) == nc
        @test length(ef.eflx_building_heat_errsoi_col) == nc
        @test length(ef.eflx_urban_ac_col) == nc
        @test length(ef.eflx_urban_heat_col) == nc
        @test length(ef.htvp_col) == nc
        @test length(ef.errsoi_col) == nc
        @test length(ef.errseb_col) == nc
        @test length(ef.errsol_col) == nc
        @test length(ef.errlon_col) == nc

        # --- Check column-level sizes (2D) ---
        @test size(ef.eflx_fgr_col) == (nc, nlevgrnd)

        # --- Check patch-level sizes (2D) ---
        @test size(ef.rresis_patch) == (np, nlevgrnd)

        # --- Check landunit-level sizes ---
        @test length(ef.eflx_traffic_lun) == nl
        @test length(ef.eflx_wasteheat_lun) == nl
        @test length(ef.eflx_ventilation_lun) == nl
        @test length(ef.eflx_heat_from_ac_lun) == nl
        @test length(ef.eflx_building_lun) == nl
        @test length(ef.eflx_urban_ac_lun) == nl
        @test length(ef.eflx_urban_heat_lun) == nl

        # --- Check gridcell-level sizes ---
        @test length(ef.eflx_dynbal_grc) == ng

        # --- Check NaN initialization for standard real fields ---
        @test all(isnan, ef.eflx_sh_tot_patch)
        @test all(isnan, ef.eflx_lh_tot_patch)
        @test all(isnan, ef.eflx_lwrad_out_patch)
        @test all(isnan, ef.eflx_snomelt_col)
        @test all(isnan, ef.eflx_dynbal_grc)
        @test all(isnan, ef.eflx_fgr_col)
        @test all(isnan, ef.rresis_patch)
        @test all(isnan, ef.btran_patch)
        @test all(isnan, ef.htvp_col)
        @test all(isnan, ef.errsoi_col)
        @test all(isnan, ef.eflx_building_lun)
        @test all(isnan, ef.taux_patch)
        @test all(isnan, ef.dgnetdT_patch)
    end

    @testset "energyflux_clean!" begin
        ef = CLM.EnergyFluxData()
        CLM.energyflux_init!(ef, 10, 5, 2, 1)
        CLM.energyflux_clean!(ef)

        # All vectors should be empty
        @test length(ef.eflx_sh_tot_patch) == 0
        @test length(ef.eflx_lh_tot_patch) == 0
        @test length(ef.eflx_lwrad_out_patch) == 0
        @test length(ef.eflx_snomelt_col) == 0
        @test length(ef.eflx_dynbal_grc) == 0
        @test length(ef.eflx_building_lun) == 0
        @test length(ef.btran_patch) == 0
        @test length(ef.errsoi_col) == 0
        @test length(ef.htvp_col) == 0
        @test size(ef.eflx_fgr_col) == (0, 0)
        @test size(ef.rresis_patch) == (0, 0)
    end

    @testset "energyflux_init_cold!" begin
        # Setup: 1 gridcell, 2 landunits (soil + urban), 4 columns, 4 patches
        nc = 4; nl = 2; np = 4; ng = 1

        ef = CLM.EnergyFluxData()
        CLM.energyflux_init!(ef, np, nc, nl, ng)

        col = CLM.ColumnData()
        CLM.column_init!(col, nc)

        lun = CLM.LandunitData()
        CLM.landunit_init!(lun, nl)

        pch = CLM.PatchData()
        CLM.patch_init!(pch, np)

        # Landunit 1: soil (non-urban)
        lun.itype[1]    = CLM.ISTSOIL
        lun.lakpoi[1]   = false
        lun.urbpoi[1]   = false
        lun.ifspecial[1] = false

        # Landunit 2: urban
        lun.itype[2]    = CLM.ISTURB_TBD
        lun.lakpoi[2]   = false
        lun.urbpoi[2]   = true
        lun.ifspecial[2] = true

        # Columns 1,2 on landunit 1; columns 3,4 on landunit 2
        col.landunit[1] = 1; col.landunit[2] = 1
        col.landunit[3] = 2; col.landunit[4] = 2

        # Patches: 1,2 on landunit 1 (soil); 3,4 on landunit 2 (urban)
        pch.column[1] = 1; pch.column[2] = 2
        pch.column[3] = 3; pch.column[4] = 4
        pch.landunit[1] = 1; pch.landunit[2] = 1
        pch.landunit[3] = 2; pch.landunit[4] = 2

        t_grnd = [272.0, 272.0, 292.0, 292.0]

        CLM.energyflux_init_cold!(ef, col, lun, pch, t_grnd,
            1:nc, 1:np, 1:nl;
            is_simple_buildtemp = false,
            is_prog_buildtemp = false)

        # Non-urban patches: urban-only fields should be SPVAL
        @test ef.eflx_lwrad_net_u_patch[1] ≈ CLM.SPVAL
        @test ef.eflx_lwrad_out_u_patch[1] ≈ CLM.SPVAL
        @test ef.eflx_lh_tot_u_patch[1]    ≈ CLM.SPVAL
        @test ef.eflx_sh_tot_u_patch[1]    ≈ CLM.SPVAL
        @test ef.eflx_soil_grnd_u_patch[1] ≈ CLM.SPVAL

        # Urban patches: urban-only fields should NOT be SPVAL (they remain NaN from init)
        @test isnan(ef.eflx_lwrad_net_u_patch[3])

        # eflx_lwrad_out should be SB * T^4
        @test ef.eflx_lwrad_out_patch[1] ≈ CLM.SB * (272.0)^4
        @test ef.eflx_lwrad_out_patch[3] ≈ CLM.SB * (292.0)^4

        # Non-urban patches: waste/traffic/ventilation/heat_from_ac should be 0
        @test ef.eflx_wasteheat_patch[1]    ≈ 0.0
        @test ef.eflx_ventilation_patch[1]  ≈ 0.0
        @test ef.eflx_heat_from_ac_patch[1] ≈ 0.0
        @test ef.eflx_traffic_patch[1]      ≈ 0.0

        # Non-urban landunit: traffic/wasteheat/ventilation should be SPVAL
        @test ef.eflx_traffic_lun[1]   ≈ CLM.SPVAL
        @test ef.eflx_wasteheat_lun[1] ≈ CLM.SPVAL
        @test ef.eflx_ventilation_lun[1] ≈ CLM.SPVAL

        # rresis should be initialized to 0
        @test all(==(0.0), ef.rresis_patch)
    end

    @testset "energyflux_init_cold! with simple_buildtemp" begin
        nc = 2; nl = 1; np = 2; ng = 1

        ef = CLM.EnergyFluxData()
        CLM.energyflux_init!(ef, np, nc, nl, ng)

        col = CLM.ColumnData()
        CLM.column_init!(col, nc)

        lun = CLM.LandunitData()
        CLM.landunit_init!(lun, nl)

        pch = CLM.PatchData()
        CLM.patch_init!(pch, np)

        lun.itype[1]    = CLM.ISTSOIL
        lun.lakpoi[1]   = false
        lun.urbpoi[1]   = false
        lun.ifspecial[1] = false

        col.landunit[1] = 1; col.landunit[2] = 1
        pch.column[1] = 1; pch.column[2] = 2
        pch.landunit[1] = 1; pch.landunit[2] = 1

        t_grnd = [280.0, 280.0]

        CLM.energyflux_init_cold!(ef, col, lun, pch, t_grnd,
            1:nc, 1:np, 1:nl;
            is_simple_buildtemp = true,
            is_prog_buildtemp = false)

        # Simple build temp: column-level building fluxes should be 0
        @test ef.eflx_building_heat_errsoi_col[1] ≈ 0.0
        @test ef.eflx_urban_ac_col[1]             ≈ 0.0
        @test ef.eflx_urban_heat_col[1]           ≈ 0.0

        # Simple build temp: anthro should be 0 for non-urban patches
        @test ef.eflx_anthro_patch[1] ≈ 0.0
    end

    @testset "energyflux_init_cold! with prog_buildtemp" begin
        nc = 2; nl = 1; np = 2; ng = 1

        ef = CLM.EnergyFluxData()
        CLM.energyflux_init!(ef, np, nc, nl, ng)

        col = CLM.ColumnData()
        CLM.column_init!(col, nc)

        lun = CLM.LandunitData()
        CLM.landunit_init!(lun, nl)

        pch = CLM.PatchData()
        CLM.patch_init!(pch, np)

        lun.itype[1]    = CLM.ISTSOIL
        lun.lakpoi[1]   = false
        lun.urbpoi[1]   = false
        lun.ifspecial[1] = false

        col.landunit[1] = 1; col.landunit[2] = 1
        pch.column[1] = 1; pch.column[2] = 2
        pch.landunit[1] = 1; pch.landunit[2] = 1

        t_grnd = [280.0, 280.0]

        CLM.energyflux_init_cold!(ef, col, lun, pch, t_grnd,
            1:nc, 1:np, 1:nl;
            is_simple_buildtemp = false,
            is_prog_buildtemp = true)

        # Prog build temp: landunit-level building fluxes should be 0 (for non-urban)
        @test ef.eflx_building_lun[1]   ≈ 0.0
        @test ef.eflx_urban_ac_lun[1]   ≈ 0.0
        @test ef.eflx_urban_heat_lun[1] ≈ 0.0
    end

    @testset "energyflux_init_acc_vars! startup" begin
        np = 5; nc = 3; nl = 1; ng = 1
        ef = CLM.EnergyFluxData()
        CLM.energyflux_init!(ef, np, nc, nl, ng)

        CLM.energyflux_init_acc_vars!(ef, 1:np; is_startup=true)

        # After startup init, btran_min should be SPVAL
        @test ef.btran_min_patch[1]      ≈ CLM.SPVAL
        @test ef.btran_min_inst_patch[1] ≈ CLM.SPVAL
        @test ef.btran_min_patch[np]     ≈ CLM.SPVAL
        @test ef.btran_min_inst_patch[np] ≈ CLM.SPVAL
    end

    @testset "field mutability" begin
        ef = CLM.EnergyFluxData()
        CLM.energyflux_init!(ef, 5, 3, 2, 1)

        # Write and verify patch-level fields
        ef.eflx_sh_tot_patch[1] = 100.0
        @test ef.eflx_sh_tot_patch[1] == 100.0

        ef.eflx_lwrad_out_patch[2] = 350.0
        @test ef.eflx_lwrad_out_patch[2] == 350.0

        ef.btran_patch[3] = 0.75
        @test ef.btran_patch[3] == 0.75

        # Write and verify column-level fields
        ef.eflx_snomelt_col[1] = 50.0
        @test ef.eflx_snomelt_col[1] == 50.0

        ef.htvp_col[2] = 2.5e6
        @test ef.htvp_col[2] == 2.5e6

        # Write and verify 2D fields
        ef.eflx_fgr_col[1, 1] = 10.0
        @test ef.eflx_fgr_col[1, 1] == 10.0

        ef.rresis_patch[1, 1] = 0.5
        @test ef.rresis_patch[1, 1] == 0.5

        # Write and verify landunit-level fields
        ef.eflx_building_lun[1] = 200.0
        @test ef.eflx_building_lun[1] == 200.0

        # Write and verify gridcell-level fields
        ef.eflx_dynbal_grc[1] = 1.0e3
        @test ef.eflx_dynbal_grc[1] == 1.0e3
    end

    @testset "re-init overwrites previous state" begin
        ef = CLM.EnergyFluxData()
        CLM.energyflux_init!(ef, 5, 3, 2, 1)
        ef.eflx_sh_tot_patch[1] = 999.0

        CLM.energyflux_init!(ef, 10, 7, 3, 2)
        @test length(ef.eflx_sh_tot_patch) == 10
        @test all(isnan, ef.eflx_sh_tot_patch)
        @test length(ef.eflx_snomelt_col) == 7
        @test length(ef.eflx_building_lun) == 3
        @test length(ef.eflx_dynbal_grc) == 2
        @test size(ef.eflx_fgr_col) == (7, nlevgrnd)
        @test size(ef.rresis_patch) == (10, nlevgrnd)
    end

    @testset "stub functions run without error" begin
        ef = CLM.EnergyFluxData()
        CLM.energyflux_init!(ef, 5, 3, 2, 1)

        # These should run without error (they are stubs)
        @test CLM.energyflux_init_history!(ef, 1:3, 1:5, 1:2, 1:1) === nothing
        @test CLM.energyflux_restart!(ef, 1:3, 1:5, 1:2) === nothing
        @test CLM.energyflux_init_acc_buffer!(ef, 1:5) === nothing
        @test CLM.energyflux_update_acc_vars!(ef, 1:5) === nothing
    end

end
