@testset "SurfaceAlbedoMod" begin

    # Ensure varpar is initialized
    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)

    # =========================================================================
    # Test SurfaceAlbedoConstants
    # =========================================================================
    @testset "SurfaceAlbedoConstants construction" begin
        con = CLM.SurfaceAlbedoConstants()
        @test con.snowveg_affects_radiation == true
        @test con.alblakwi == [0.10, 0.10]
        @test con.lake_melt_icealb == [0.10, 0.10]
        @test length(con.isoicol) == 0
    end

    # =========================================================================
    # Test surface_albedo_init_time_const!
    # =========================================================================
    @testset "surface_albedo_init_time_const! 8 colors" begin
        con = CLM.SurfaceAlbedoConstants()
        mxsoil_color = 8
        ng = 3
        nc = 5
        soic2d = [3, 5, 1]  # soil color per gridcell
        col_gridcell = [1, 1, 2, 2, 3]
        bounds_col = 1:nc
        bounds_grc = 1:ng

        CLM.surface_albedo_init_time_const!(con, mxsoil_color, soic2d,
            col_gridcell, bounds_col, bounds_grc)

        # Check soil color mapping
        @test con.isoicol[1] == 3
        @test con.isoicol[2] == 3
        @test con.isoicol[3] == 5
        @test con.isoicol[4] == 5
        @test con.isoicol[5] == 1

        # Check albsat dimensions
        @test size(con.albsat) == (8, 2)
        @test size(con.albdry) == (8, 2)

        # Check a few values from the 8-color table
        @test con.albsat[1, 1] ≈ 0.12
        @test con.albsat[8, 2] ≈ 0.10
        @test con.albdry[1, 1] ≈ 0.24
        @test con.albdry[8, 2] ≈ 0.20

        # Saturated albedo should be less than dry albedo
        for ic in 1:8
            for ib in 1:2
                @test con.albsat[ic, ib] <= con.albdry[ic, ib]
            end
        end

        # Check alblakwi
        @test con.alblakwi == [0.10, 0.10]
    end

    @testset "surface_albedo_init_time_const! 20 colors" begin
        con = CLM.SurfaceAlbedoConstants()
        mxsoil_color = 20
        ng = 2
        nc = 3
        soic2d = [10, 15]
        col_gridcell = [1, 2, 2]
        bounds_col = 1:nc
        bounds_grc = 1:ng

        CLM.surface_albedo_init_time_const!(con, mxsoil_color, soic2d,
            col_gridcell, bounds_col, bounds_grc)

        @test size(con.albsat) == (20, 2)
        @test size(con.albdry) == (20, 2)
        @test con.albsat[1, 1] ≈ 0.25
        @test con.albsat[20, 2] ≈ 0.08
        @test con.albdry[20, 2] ≈ 0.16
    end

    @testset "surface_albedo_init_time_const! unsupported color" begin
        con = CLM.SurfaceAlbedoConstants()
        @test_throws ErrorException CLM.surface_albedo_init_time_const!(
            con, 12, [1], [1], 1:1, 1:1)
    end

    # =========================================================================
    # Test soil_albedo!
    # =========================================================================
    @testset "soil_albedo! basic soil" begin
        nc = 3
        np = 3
        ng = 2
        nlevsno = CLM.varpar.nlevsno
        nlevgrnd = CLM.varpar.nlevgrnd

        # Set up surface albedo data
        sa = CLM.SurfaceAlbedoData()
        CLM.surfalb_init!(sa, np, nc, ng)
        CLM.surfalb_init_cold!(sa, 1:nc, 1:np)

        # Set up constants
        con = CLM.SurfaceAlbedoConstants()
        con.albsat = [0.12 0.24; 0.10 0.20]  # 2 color classes
        con.albdry = [0.24 0.48; 0.20 0.40]
        con.isoicol = [1, 2, 1]

        # Set up column data
        col = CLM.ColumnData()
        CLM.column_init!(col, nc)
        col.landunit .= [1, 1, 1]
        col.snl .= [0, 0, 0]

        # Set up landunit data
        lun = CLM.LandunitData()
        CLM.landunit_init!(lun, 1)
        lun.itype .= [CLM.ISTSOIL]

        # Cosine zenith (positive = daytime)
        coszen_col = [0.5, 0.0, 0.8]

        # Temperature
        temp = CLM.TemperatureData()
        CLM.temperature_init!(temp, np, nc, 1, ng)
        temp.t_grnd_col .= 280.0

        # Water state
        wsb = CLM.WaterStateBulkData()
        nlevtot = nlevsno + CLM.varpar.nlevmaxurbgrnd
        CLM.waterstatebulk_init!(wsb, nc, np, 1, ng)
        wsb.ws.h2osoi_vol_col[1, 1] = 0.3  # moist soil
        wsb.ws.h2osoi_vol_col[2, 1] = 0.1
        wsb.ws.h2osoi_vol_col[3, 1] = 0.0  # dry soil

        # Lake state (not used for soil)
        ls = CLM.LakeStateData()
        CLM.lakestate_init!(ls, nc, np)

        mask_nourbanc = trues(nc)

        CLM.soil_albedo!(sa, con, col, lun,
            coszen_col, temp, wsb, ls,
            mask_nourbanc, 1:nc)

        # Column 1: soil with coszen=0.5, moist soil (h2osoi_vol=0.3)
        # inc = max(0.11 - 0.40*0.3, 0) = max(0.11-0.12, 0) = 0.0
        # albsod = min(albsat[1,ib] + 0, albdry[1,ib]) = albsat[1,ib]
        @test sa.albsod_col[1, 1] ≈ 0.12
        @test sa.albsod_col[1, 2] ≈ 0.24
        @test sa.albsoi_col[1, 1] ≈ sa.albsod_col[1, 1]

        # Column 2: coszen=0, should remain at cold-start values
        @test sa.albsod_col[2, 1] ≈ 0.2  # unchanged from cold init

        # Column 3: soil with coszen=0.8, dry soil (h2osoi_vol=0.0)
        # inc = max(0.11 - 0.40*0, 0) = 0.11
        # albsod = min(albsat[1,ib]+0.11, albdry[1,ib])
        @test sa.albsod_col[3, 1] ≈ min(0.12 + 0.11, 0.24)  # 0.23
        @test sa.albsod_col[3, 2] ≈ min(0.24 + 0.11, 0.48)  # 0.35
    end

    @testset "soil_albedo! land ice" begin
        nc = 1
        np = 1
        ng = 1
        sa = CLM.SurfaceAlbedoData()
        CLM.surfalb_init!(sa, np, nc, ng)
        CLM.surfalb_init_cold!(sa, 1:nc, 1:np)

        con = CLM.SurfaceAlbedoConstants()
        con.albsat = zeros(1, 2)
        con.albdry = zeros(1, 2)
        con.isoicol = [1]

        col = CLM.ColumnData()
        CLM.column_init!(col, nc)
        col.landunit .= [1]
        col.snl .= [0]

        lun = CLM.LandunitData()
        CLM.landunit_init!(lun, 1)
        lun.itype .= [CLM.ISTICE]

        coszen_col = [0.6]
        temp = CLM.TemperatureData()
        CLM.temperature_init!(temp, np, nc, 1, ng)
        temp.t_grnd_col .= 260.0

        wsb = CLM.WaterStateBulkData()
        CLM.waterstatebulk_init!(wsb, nc, np, 1, ng)

        ls = CLM.LakeStateData()
        CLM.lakestate_init!(ls, nc, np)

        mask = trues(nc)
        CLM.soil_albedo!(sa, con, col, lun, coszen_col, temp, wsb, ls, mask, 1:nc)

        @test sa.albsod_col[1, 1] ≈ CLM.ALBICE[1]  # 0.80
        @test sa.albsod_col[1, 2] ≈ CLM.ALBICE[2]  # 0.55
        @test sa.albsoi_col[1, 1] ≈ sa.albsod_col[1, 1]
    end

    # =========================================================================
    # Test two_stream! basic functionality
    # =========================================================================
    @testset "two_stream! single layer canopy" begin
        nc = 1
        np = 1
        ng = 1
        numrad = CLM.NUMRAD

        sa = CLM.SurfaceAlbedoData()
        CLM.surfalb_init!(sa, np, nc, ng)
        CLM.surfalb_init_cold!(sa, 1:nc, 1:np)

        # Set ground albedos
        sa.albgrd_col[1, 1] = 0.15
        sa.albgrd_col[1, 2] = 0.30
        sa.albgri_col[1, 1] = 0.15
        sa.albgri_col[1, 2] = 0.30
        sa.albsod_col[1, 1] = 0.15
        sa.albsod_col[1, 2] = 0.30
        sa.albsoi_col[1, 1] = 0.15
        sa.albsoi_col[1, 2] = 0.30
        sa.nrad_patch[1] = 1
        sa.tlai_z_patch[1, 1] = 2.0
        sa.tsai_z_patch[1, 1] = 0.5

        patchdata = CLM.PatchData()
        CLM.patch_init!(patchdata, np)
        patchdata.column .= [1]
        patchdata.itype .= [1]
        patchdata.landunit .= [1]

        col = CLM.ColumnData()
        CLM.column_init!(col, nc)

        cs = CLM.CanopyStateData()
        CLM.canopystate_init!(cs, np)
        cs.elai_patch .= [2.0]
        cs.esai_patch .= [0.5]

        temp = CLM.TemperatureData()
        CLM.temperature_init!(temp, np, nc, 1, ng)
        temp.t_veg_patch .= [290.0]  # above freezing

        wd = CLM.WaterDiagnosticBulkData()
        CLM.waterdiagnosticbulk_init!(wd, nc, np, 1, ng)
        wd.fwet_patch .= [0.0]
        wd.fcansno_patch .= [0.0]

        coszen_patch = [0.5]
        # PFT optical properties: rhol, rhos, taul, taus for PFT type 1
        pftcon_xl = fill(0.01, 80)  # leaf orientation index

        rho = zeros(np, numrad)
        tau = zeros(np, numrad)
        rho[1, 1] = 0.10; rho[1, 2] = 0.45  # leaf+stem reflectance
        tau[1, 1] = 0.05; tau[1, 2] = 0.25  # leaf+stem transmittance

        mask_vegsol = trues(np)

        CLM.two_stream!(sa, patchdata, col, cs, temp, wd,
            coszen_patch, rho, tau, pftcon_xl,
            mask_vegsol, 1:np)

        # Basic sanity checks:
        # Surface albedo should be between 0 and 1
        @test 0.0 <= sa.albd_patch[1, 1] <= 1.0
        @test 0.0 <= sa.albd_patch[1, 2] <= 1.0
        @test 0.0 <= sa.albi_patch[1, 1] <= 1.0
        @test 0.0 <= sa.albi_patch[1, 2] <= 1.0

        # Transmitted direct beam should be less than 1 (canopy absorbs)
        @test 0.0 <= sa.ftdd_patch[1, 1] < 1.0
        @test 0.0 <= sa.ftdd_patch[1, 2] < 1.0

        # Absorbed fluxes should be non-negative
        @test sa.fabd_patch[1, 1] >= 0.0
        @test sa.fabi_patch[1, 1] >= 0.0

        # Energy conservation: albd + fabd + (1-albgrd)*ftdd + (1-albgri)*ftid ≈ 1
        for ib in 1:numrad
            total_direct = sa.albd_patch[1, ib] + sa.fabd_patch[1, ib] +
                (1.0 - sa.albgrd_col[1, ib]) * sa.ftdd_patch[1, ib] +
                (1.0 - sa.albgri_col[1, ib]) * sa.ftid_patch[1, ib]
            @test total_direct ≈ 1.0 atol=1e-10

            total_diffuse = sa.albi_patch[1, ib] + sa.fabi_patch[1, ib] +
                (1.0 - sa.albgri_col[1, ib]) * sa.ftii_patch[1, ib]
            @test total_diffuse ≈ 1.0 atol=1e-10
        end

        # Sun+shade fractions should sum to total
        @test sa.fabd_sun_patch[1, 1] + sa.fabd_sha_patch[1, 1] ≈ sa.fabd_patch[1, 1] atol=1e-10
        @test sa.fabi_sun_patch[1, 1] + sa.fabi_sha_patch[1, 1] ≈ sa.fabi_patch[1, 1] atol=1e-10

        # Sunlit fraction should be between 0 and 1
        @test 0.0 < sa.fsun_z_patch[1, 1] < 1.0

        # VIS albedo should be lower than NIR albedo (generally true for vegetation)
        @test sa.albd_patch[1, 1] < sa.albd_patch[1, 2]
    end

    # =========================================================================
    # Test two_stream! Snow-Free only mode
    # =========================================================================
    @testset "two_stream! SFonly mode" begin
        nc = 1
        np = 1
        ng = 1
        numrad = CLM.NUMRAD

        sa = CLM.SurfaceAlbedoData()
        CLM.surfalb_init!(sa, np, nc, ng)
        CLM.surfalb_init_cold!(sa, 1:nc, 1:np; use_SSRE=true)

        sa.albgrd_col[1, 1] = 0.20
        sa.albgrd_col[1, 2] = 0.35
        sa.albgri_col[1, 1] = 0.20
        sa.albgri_col[1, 2] = 0.35
        sa.albsod_col[1, 1] = 0.12
        sa.albsod_col[1, 2] = 0.24
        sa.albsoi_col[1, 1] = 0.12
        sa.albsoi_col[1, 2] = 0.24
        sa.nrad_patch[1] = 1
        sa.tlai_z_patch[1, 1] = 1.5
        sa.tsai_z_patch[1, 1] = 0.3

        patchdata = CLM.PatchData()
        CLM.patch_init!(patchdata, np)
        patchdata.column .= [1]
        patchdata.itype .= [1]
        patchdata.landunit .= [1]

        col = CLM.ColumnData()
        CLM.column_init!(col, nc)

        cs = CLM.CanopyStateData()
        CLM.canopystate_init!(cs, np)
        cs.elai_patch .= [1.5]
        cs.esai_patch .= [0.3]

        temp = CLM.TemperatureData()
        CLM.temperature_init!(temp, np, nc, 1, ng)
        temp.t_veg_patch .= [290.0]

        wd = CLM.WaterDiagnosticBulkData()
        CLM.waterdiagnosticbulk_init!(wd, nc, np, 1, ng)
        wd.fwet_patch .= [0.0]
        wd.fcansno_patch .= [0.0]

        coszen_patch = [0.5]
        pftcon_xl = fill(0.01, 80)

        rho = zeros(np, numrad)
        tau = zeros(np, numrad)
        rho[1, 1] = 0.10; rho[1, 2] = 0.45
        tau[1, 1] = 0.05; tau[1, 2] = 0.25

        mask = trues(np)

        CLM.two_stream!(sa, patchdata, col, cs, temp, wd,
            coszen_patch, rho, tau, pftcon_xl,
            mask, 1:np; SFonly=true)

        # SF albedos should be set and reasonable
        @test 0.0 <= sa.albdSF_patch[1, 1] <= 1.0
        @test 0.0 <= sa.albdSF_patch[1, 2] <= 1.0
        @test 0.0 <= sa.albiSF_patch[1, 1] <= 1.0
        @test 0.0 <= sa.albiSF_patch[1, 2] <= 1.0
    end

    # =========================================================================
    # Test module-level constants
    # =========================================================================
    @testset "module constants" begin
        @test CLM.ALBICE == [0.80, 0.55]
        @test CLM.ALBLAK == [0.60, 0.40]
        @test CLM.CALB ≈ 95.6
        @test CLM.MPE_ALBEDO ≈ 1.0e-6
    end

    # =========================================================================
    # Test global surfalb_con instance
    # =========================================================================
    @testset "global surfalb_con" begin
        @test CLM.surfalb_con isa CLM.SurfaceAlbedoConstants
        @test CLM.surfalb_con.snowveg_affects_radiation == true
    end

end
