@testset "SurfaceRadiationMod" begin

    # Ensure varpar is initialized
    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)

    nlevsno = CLM.varpar.nlevsno
    numrad = CLM.NUMRAD

    # =========================================================================
    # Test SurfaceRadiationData struct and init/clean
    # =========================================================================
    @testset "SurfaceRadiationData construction and init" begin
        sr = CLM.SurfaceRadiationData()
        @test sr.sfc_frc_aer_patch == Float64[]

        np = 5
        CLM.surfrad_init!(sr, np)
        @test length(sr.sfc_frc_aer_patch) == np
        @test length(sr.parveg_ln_patch) == np
        @test length(sr.fsr_vis_d_patch) == np
        @test length(sr.fsds_vis_d_patch) == np
        @test length(sr.fsds_sno_vd_patch) == np
        @test all(isnan, sr.sfc_frc_aer_patch)
        @test all(isnan, sr.fsr_vis_d_patch)

        CLM.surfrad_clean!(sr)
        @test sr.sfc_frc_aer_patch == Float64[]
        @test sr.fsr_vis_d_patch == Float64[]
    end

    @testset "surfrad_init_history!" begin
        sr = CLM.SurfaceRadiationData()
        np = 3
        CLM.surfrad_init!(sr, np)
        CLM.surfrad_init_history!(sr, 1:np; use_snicar_frc=true, use_SSRE=true)
        @test sr.fsds_vis_d_patch[1] == CLM.SPVAL
        @test sr.sfc_frc_aer_patch[1] == CLM.SPVAL
        @test sr.ssre_fsr_vis_d_patch[1] == CLM.SPVAL
        @test sr.fsr_sno_vd_patch[1] == CLM.SPVAL
    end

    @testset "surfrad_init_cold!" begin
        sr = CLM.SurfaceRadiationData()
        np = 2
        CLM.surfrad_init!(sr, np)
        # Should be a no-op (no errors)
        CLM.surfrad_init_cold!(sr, 1:np)
        @test true  # just verify no error
    end

    # =========================================================================
    # Test is_near_local_noon
    # =========================================================================
    @testset "is_near_local_noon" begin
        # At londeg=0, local noon is at 43200s (12:00 UTC)
        @test CLM.is_near_local_noon(0.0; deltasec=1800, current_tod=43200.0) == true
        @test CLM.is_near_local_noon(0.0; deltasec=1800, current_tod=45001.0) == false

        # At londeg=90 (east), local noon is at 43200 - 90*240 = 43200-21600 = 21600s (6:00 UTC)
        @test CLM.is_near_local_noon(90.0; deltasec=1800, current_tod=21600.0) == true
        @test CLM.is_near_local_noon(90.0; deltasec=1800, current_tod=43200.0) == false

        # At londeg=-90 (west), local noon is at 43200 + 90*240 = 43200+21600 = 64800s (18:00 UTC)
        @test CLM.is_near_local_noon(-90.0; deltasec=1800, current_tod=64800.0) == true
    end

    # =========================================================================
    # Test canopy_sun_shade_fracs!
    # =========================================================================
    @testset "canopy_sun_shade_fracs! basic" begin
        np = 2
        nc = 2
        ng = 1

        # Set up surface albedo data
        sa = CLM.SurfaceAlbedoData()
        CLM.surfalb_init!(sa, np, nc, ng)
        CLM.surfalb_init_cold!(sa, 1:nc, 1:np)

        # Canopy with one layer
        sa.nrad_patch[1] = 1
        sa.nrad_patch[2] = 1
        sa.tlai_z_patch[1, 1] = 3.0
        sa.tlai_z_patch[2, 1] = 1.0
        sa.fsun_z_patch[1, 1] = 0.5
        sa.fsun_z_patch[2, 1] = 0.8

        # PAR absorption fractions
        sa.fabd_sun_z_patch[1, 1] = 0.3
        sa.fabd_sha_z_patch[1, 1] = 0.1
        sa.fabi_sun_z_patch[1, 1] = 0.25
        sa.fabi_sha_z_patch[1, 1] = 0.08
        sa.fabd_sun_z_patch[2, 1] = 0.4
        sa.fabd_sha_z_patch[2, 1] = 0.05
        sa.fabi_sun_z_patch[2, 1] = 0.35
        sa.fabi_sha_z_patch[2, 1] = 0.04

        # Canopy state
        cs = CLM.CanopyStateData()
        CLM.canopystate_init!(cs, np)
        cs.elai_patch .= [3.0, 1.0]

        # Solar absorbed
        sol = CLM.SolarAbsorbedData()
        CLM.solarabs_init!(sol, np, 1)

        # Forcing
        forc_solad_col = zeros(nc, numrad)
        forc_solad_col[1, 1] = 500.0  # direct VIS on column 1
        forc_solad_col[2, 1] = 300.0
        forc_solai = zeros(ng, numrad)
        forc_solai[1, 1] = 100.0      # diffuse VIS on gridcell 1

        # Patch data
        pch = CLM.PatchData()
        CLM.patch_init!(pch, np)
        pch.column .= [1, 2]
        pch.gridcell .= [1, 1]
        pch.landunit .= [1, 1]

        mask = trues(np)

        CLM.canopy_sun_shade_fracs!(sa, cs, sol, forc_solad_col, forc_solai, pch, mask, 1:np)

        # Check laisun/laisha
        @test cs.laisun_z_patch[1, 1] ≈ 3.0 * 0.5  # 1.5
        @test cs.laisha_z_patch[1, 1] ≈ 3.0 * 0.5  # 1.5
        @test cs.laisun_patch[1] ≈ 1.5
        @test cs.laisha_patch[1] ≈ 1.5

        @test cs.laisun_z_patch[2, 1] ≈ 1.0 * 0.8  # 0.8
        @test cs.laisha_z_patch[2, 1] ≈ 1.0 * 0.2  # 0.2

        # Check fsun
        @test cs.fsun_patch[1] ≈ 1.5 / 3.0  # 0.5
        @test cs.fsun_patch[2] ≈ 0.8 / 1.0  # 0.8

        # Check PAR absorption
        # parsun_z(1,1) = 500*0.3 + 100*0.25 = 150 + 25 = 175
        @test sol.parsun_z_patch[1, 1] ≈ 175.0
        # parsha_z(1,1) = 500*0.1 + 100*0.08 = 50 + 8 = 58
        @test sol.parsha_z_patch[1, 1] ≈ 58.0

        # parsun_z(2,1) = 300*0.4 + 100*0.35 = 120 + 35 = 155
        @test sol.parsun_z_patch[2, 1] ≈ 155.0
    end

    @testset "canopy_sun_shade_fracs! zero LAI" begin
        np = 1
        nc = 1
        ng = 1

        sa = CLM.SurfaceAlbedoData()
        CLM.surfalb_init!(sa, np, nc, ng)
        CLM.surfalb_init_cold!(sa, 1:nc, 1:np)
        sa.nrad_patch[1] = 1
        sa.tlai_z_patch[1, 1] = 0.0
        sa.fsun_z_patch[1, 1] = 0.0

        cs = CLM.CanopyStateData()
        CLM.canopystate_init!(cs, np)
        cs.elai_patch .= [0.0]

        sol = CLM.SolarAbsorbedData()
        CLM.solarabs_init!(sol, np, 1)

        forc_solad_col = zeros(nc, numrad)
        forc_solai = zeros(ng, numrad)

        pch = CLM.PatchData()
        CLM.patch_init!(pch, np)
        pch.column .= [1]
        pch.gridcell .= [1]

        mask = trues(np)

        CLM.canopy_sun_shade_fracs!(sa, cs, sol, forc_solad_col, forc_solai, pch, mask, 1:np)

        @test cs.fsun_patch[1] ≈ 0.0
        @test cs.laisun_patch[1] ≈ 0.0
        @test cs.laisha_patch[1] ≈ 0.0
    end

    # =========================================================================
    # Test surface_radiation!
    # =========================================================================
    @testset "surface_radiation! basic non-urban" begin
        np = 1
        nc = 1
        ng = 1
        nl = 1

        # --- Set up all data structures ---
        sa = CLM.SurfaceAlbedoData()
        CLM.surfalb_init!(sa, np, nc, ng)
        CLM.surfalb_init_cold!(sa, 1:nc, 1:np)

        cs = CLM.CanopyStateData()
        CLM.canopystate_init!(cs, np)
        cs.elai_patch .= [2.0]
        cs.esai_patch .= [0.5]
        cs.tlai_patch .= [2.0]
        cs.fsun_patch .= [0.5]

        sol = CLM.SolarAbsorbedData()
        CLM.solarabs_init!(sol, np, nl)

        sr = CLM.SurfaceRadiationData()
        CLM.surfrad_init!(sr, np)

        wd = CLM.WaterDiagnosticBulkData()
        CLM.waterdiagnosticbulk_init!(wd, nc, np, nl, ng)
        wd.snow_depth_col .= [0.0]
        wd.frac_sno_col .= [0.0]
        wd.frac_sno_eff_col .= [0.0]

        col = CLM.ColumnData()
        CLM.column_init!(col, nc)
        col.snl .= [0]
        col.landunit .= [1]

        lun = CLM.LandunitData()
        CLM.landunit_init!(lun, nl)
        lun.itype .= [CLM.ISTSOIL]

        grc_data = CLM.GridcellData()
        CLM.gridcell_init!(grc_data, ng)
        grc_data.londeg .= [0.0]

        pch = CLM.PatchData()
        CLM.patch_init!(pch, np)
        pch.column .= [1]
        pch.gridcell .= [1]
        pch.landunit .= [1]

        # Set up albedos: simple case
        sa.albgrd_col[1, 1] = 0.15
        sa.albgrd_col[1, 2] = 0.30
        sa.albgri_col[1, 1] = 0.15
        sa.albgri_col[1, 2] = 0.30
        sa.albsod_col[1, 1] = 0.12
        sa.albsod_col[1, 2] = 0.24
        sa.albsoi_col[1, 1] = 0.12
        sa.albsoi_col[1, 2] = 0.24
        sa.albsnd_hst_col[1, 1] = 0.6
        sa.albsnd_hst_col[1, 2] = 0.4
        sa.albsni_hst_col[1, 1] = 0.6
        sa.albsni_hst_col[1, 2] = 0.4

        # Canopy fluxes (no canopy: ftdd=1, ftid=0, ftii=1, fabd=0, fabi=0)
        sa.fabd_patch[1, 1] = 0.0
        sa.fabd_patch[1, 2] = 0.0
        sa.fabi_patch[1, 1] = 0.0
        sa.fabi_patch[1, 2] = 0.0
        sa.ftdd_patch[1, 1] = 1.0
        sa.ftdd_patch[1, 2] = 1.0
        sa.ftid_patch[1, 1] = 0.0
        sa.ftid_patch[1, 2] = 0.0
        sa.ftii_patch[1, 1] = 1.0
        sa.ftii_patch[1, 2] = 1.0

        sa.albd_patch[1, 1] = 0.15
        sa.albd_patch[1, 2] = 0.30
        sa.albi_patch[1, 1] = 0.15
        sa.albi_patch[1, 2] = 0.30

        # Forcing
        forc_solad_col = zeros(nc, numrad)
        forc_solad_col[1, 1] = 400.0  # direct VIS
        forc_solad_col[1, 2] = 200.0  # direct NIR
        forc_solai = zeros(ng, numrad)
        forc_solai[1, 1] = 100.0      # diffuse VIS
        forc_solai[1, 2] = 50.0       # diffuse NIR

        mask_nourban = trues(np)
        mask_urban = falses(np)

        CLM.surface_radiation!(sa, cs, sol, sr, wd, col, lun, grc_data, pch,
            forc_solad_col, forc_solai, mask_nourban, mask_urban, 1:np;
            dtime=3600.0, current_tod=43200.0)

        # With no canopy (fabd=0, ftdd=1, ftid=0, ftii=1):
        # sabv = 0
        @test sol.sabv_patch[1] ≈ 0.0

        # sabg = sum over bands of (trd*(1-albgrd) + tri*(1-albgri))
        # VIS: trd=400*1=400, tri=400*0+100*1=100
        #      absorbed = 400*(1-0.15) + 100*(1-0.15) = 340 + 85 = 425
        # NIR: trd=200*1=200, tri=200*0+50*1=50
        #      absorbed = 200*(1-0.30) + 50*(1-0.30) = 140 + 35 = 175
        @test sol.sabg_patch[1] ≈ 425.0 + 175.0  # 600.0

        # fsa = sabv + sabg = 600
        @test sol.fsa_patch[1] ≈ 600.0

        # fsr = rvis + rnir = (0.15*400 + 0.15*100) + (0.30*200 + 0.30*50)
        #     = (60 + 15) + (60 + 15) = 75 + 75 = 150
        @test sol.fsr_patch[1] ≈ 150.0

        # Energy conservation: fsa + fsr = total incoming
        # Total incoming = (400+100) + (200+50) = 750
        @test sol.fsa_patch[1] + sol.fsr_patch[1] ≈ 750.0

        # With snl=0, all absorbed in top soil layer
        @test sol.sabg_lyr_patch[1, nlevsno + 1] ≈ sol.sabg_patch[1]

        # Incident VIS
        @test sr.fsds_vis_d_patch[1] ≈ 400.0
        @test sr.fsds_vis_i_patch[1] ≈ 100.0

        # Reflected VIS
        @test sr.fsr_vis_d_patch[1] ≈ 0.15 * 400.0  # 60.0
        @test sr.fsr_vis_i_patch[1] ≈ 0.15 * 100.0  # 15.0

        # Since snl=0, snow diagnostics should be SPVAL
        @test sr.fsds_sno_vd_patch[1] == CLM.SPVAL
    end

    @testset "surface_radiation! with snow layers" begin
        np = 1
        nc = 1
        ng = 1
        nl = 1

        sa = CLM.SurfaceAlbedoData()
        CLM.surfalb_init!(sa, np, nc, ng)
        CLM.surfalb_init_cold!(sa, 1:nc, 1:np)

        cs = CLM.CanopyStateData()
        CLM.canopystate_init!(cs, np)
        cs.elai_patch .= [0.0]
        cs.esai_patch .= [0.0]
        cs.tlai_patch .= [0.0]
        cs.fsun_patch .= [0.0]

        sol = CLM.SolarAbsorbedData()
        CLM.solarabs_init!(sol, np, nl)

        sr = CLM.SurfaceRadiationData()
        CLM.surfrad_init!(sr, np)

        wd = CLM.WaterDiagnosticBulkData()
        CLM.waterdiagnosticbulk_init!(wd, nc, np, nl, ng)
        wd.snow_depth_col .= [0.5]
        wd.frac_sno_col .= [0.9]
        wd.frac_sno_eff_col .= [0.9]

        col = CLM.ColumnData()
        CLM.column_init!(col, nc)
        col.snl .= [-2]  # 2 snow layers
        col.landunit .= [1]

        lun = CLM.LandunitData()
        CLM.landunit_init!(lun, nl)
        lun.itype .= [CLM.ISTSOIL]

        grc_data = CLM.GridcellData()
        CLM.gridcell_init!(grc_data, ng)
        grc_data.londeg .= [0.0]

        pch = CLM.PatchData()
        CLM.patch_init!(pch, np)
        pch.column .= [1]
        pch.gridcell .= [1]
        pch.landunit .= [1]

        # Set up albedos
        sa.albgrd_col[1, 1] = 0.70
        sa.albgrd_col[1, 2] = 0.50
        sa.albgri_col[1, 1] = 0.70
        sa.albgri_col[1, 2] = 0.50
        sa.albsod_col[1, 1] = 0.12
        sa.albsod_col[1, 2] = 0.24
        sa.albsoi_col[1, 1] = 0.12
        sa.albsoi_col[1, 2] = 0.24
        sa.albsnd_hst_col[1, 1] = 0.8
        sa.albsnd_hst_col[1, 2] = 0.5
        sa.albsni_hst_col[1, 1] = 0.8
        sa.albsni_hst_col[1, 2] = 0.5

        # No canopy
        sa.fabd_patch[1, :] .= 0.0
        sa.fabi_patch[1, :] .= 0.0
        sa.ftdd_patch[1, :] .= 1.0
        sa.ftid_patch[1, :] .= 0.0
        sa.ftii_patch[1, :] .= 1.0
        sa.albd_patch[1, :] .= [0.70, 0.50]
        sa.albi_patch[1, :] .= [0.70, 0.50]

        # Set up SNICAR flux absorption factors
        # For 2 snow layers (snl=-2), active layers are -1,0 in Fortran → indices nlevsno-1, nlevsno in Julia
        # Plus soil layer at index nlevsno+1
        nlev_sno1 = nlevsno + 1
        sa.flx_absdv_col = fill(0.0, nc, nlev_sno1)
        sa.flx_absdn_col = fill(0.0, nc, nlev_sno1)
        sa.flx_absiv_col = fill(0.0, nc, nlev_sno1)
        sa.flx_absin_col = fill(0.0, nc, nlev_sno1)

        # Set absorption in the 2 active snow layers and soil
        # snl=-2: active layers are fortran indices -1,0,1 → julia nlevsno-1, nlevsno, nlevsno+1
        sa.flx_absdv_col[1, nlevsno - 1] = 0.10
        sa.flx_absdv_col[1, nlevsno]     = 0.15
        sa.flx_absdv_col[1, nlevsno + 1] = 0.05

        sa.flx_absdn_col[1, nlevsno - 1] = 0.08
        sa.flx_absdn_col[1, nlevsno]     = 0.20
        sa.flx_absdn_col[1, nlevsno + 1] = 0.22

        sa.flx_absiv_col[1, nlevsno - 1] = 0.12
        sa.flx_absiv_col[1, nlevsno]     = 0.13
        sa.flx_absiv_col[1, nlevsno + 1] = 0.05

        sa.flx_absin_col[1, nlevsno - 1] = 0.10
        sa.flx_absin_col[1, nlevsno]     = 0.18
        sa.flx_absin_col[1, nlevsno + 1] = 0.22

        forc_solad_col = zeros(nc, numrad)
        forc_solad_col[1, 1] = 400.0
        forc_solad_col[1, 2] = 200.0
        forc_solai = zeros(ng, numrad)
        forc_solai[1, 1] = 100.0
        forc_solai[1, 2] = 50.0

        mask_nourban = trues(np)
        mask_urban = falses(np)

        CLM.surface_radiation!(sa, cs, sol, sr, wd, col, lun, grc_data, pch,
            forc_solad_col, forc_solai, mask_nourban, mask_urban, 1:np;
            dtime=3600.0, current_tod=43200.0)

        # With snow, snow diagnostics should NOT be SPVAL
        @test sr.fsds_sno_vd_patch[1] ≈ 400.0
        @test sr.fsds_sno_nd_patch[1] ≈ 200.0
        @test sr.fsds_sno_vi_patch[1] ≈ 100.0
        @test sr.fsds_sno_ni_patch[1] ≈ 50.0

        # fsr_sno_vd = fsds_vis_d * albsnd_hst[1]
        @test sr.fsr_sno_vd_patch[1] ≈ 400.0 * 0.8

        # sabg_pen should be set (istsoil)
        @test !isnan(sol.sabg_pen_patch[1])

        # fsa should be positive
        @test sol.fsa_patch[1] > 0.0
    end

    @testset "surface_radiation! urban patches" begin
        np = 1
        nc = 1
        ng = 1
        nl = 1

        sa = CLM.SurfaceAlbedoData()
        CLM.surfalb_init!(sa, np, nc, ng)
        CLM.surfalb_init_cold!(sa, 1:nc, 1:np)

        cs = CLM.CanopyStateData()
        CLM.canopystate_init!(cs, np)
        cs.fsun_patch .= [0.5]

        sol = CLM.SolarAbsorbedData()
        CLM.solarabs_init!(sol, np, nl)

        sr = CLM.SurfaceRadiationData()
        CLM.surfrad_init!(sr, np)

        wd = CLM.WaterDiagnosticBulkData()
        CLM.waterdiagnosticbulk_init!(wd, nc, np, nl, ng)

        col = CLM.ColumnData()
        CLM.column_init!(col, nc)
        col.snl .= [0]

        lun = CLM.LandunitData()
        CLM.landunit_init!(lun, nl)
        lun.itype .= [CLM.ISTURB_TBD]

        grc_data = CLM.GridcellData()
        CLM.gridcell_init!(grc_data, ng)
        grc_data.londeg .= [0.0]

        pch = CLM.PatchData()
        CLM.patch_init!(pch, np)
        pch.column .= [1]
        pch.gridcell .= [1]
        pch.landunit .= [1]

        sa.albd_patch[1, :] .= [0.20, 0.30]
        sa.albi_patch[1, :] .= [0.20, 0.30]

        forc_solad_col = zeros(nc, numrad)
        forc_solad_col[1, 1] = 300.0
        forc_solad_col[1, 2] = 150.0
        forc_solai = zeros(ng, numrad)
        forc_solai[1, 1] = 80.0
        forc_solai[1, 2] = 40.0

        mask_nourban = falses(np)
        mask_urban = trues(np)

        CLM.surface_radiation!(sa, cs, sol, sr, wd, col, lun, grc_data, pch,
            forc_solad_col, forc_solai, mask_nourban, mask_urban, 1:np;
            dtime=3600.0, current_tod=43200.0)

        # fsun should be zeroed for urban
        @test cs.fsun_patch[1] ≈ 0.0

        # Incident solar should be set
        @test sr.fsds_vis_d_patch[1] ≈ 300.0
        @test sr.fsds_vis_i_patch[1] ≈ 80.0

        # Reflected solar
        @test sr.fsr_vis_d_patch[1] ≈ 0.20 * 300.0
        @test sol.fsr_nir_d_patch[1] ≈ 0.30 * 150.0

        # fsr = sum of all reflected
        expected_fsr = 0.20 * 300.0 + 0.30 * 150.0 + 0.20 * 80.0 + 0.30 * 40.0
        @test sol.fsr_patch[1] ≈ expected_fsr
    end

    @testset "surface_radiation! energy conservation" begin
        np = 1
        nc = 1
        ng = 1
        nl = 1

        sa = CLM.SurfaceAlbedoData()
        CLM.surfalb_init!(sa, np, nc, ng)
        CLM.surfalb_init_cold!(sa, 1:nc, 1:np)

        cs = CLM.CanopyStateData()
        CLM.canopystate_init!(cs, np)
        cs.elai_patch .= [1.0]

        sol = CLM.SolarAbsorbedData()
        CLM.solarabs_init!(sol, np, nl)

        sr = CLM.SurfaceRadiationData()
        CLM.surfrad_init!(sr, np)

        wd = CLM.WaterDiagnosticBulkData()
        CLM.waterdiagnosticbulk_init!(wd, nc, np, nl, ng)
        wd.snow_depth_col .= [0.0]
        wd.frac_sno_col .= [0.0]

        col = CLM.ColumnData()
        CLM.column_init!(col, nc)
        col.snl .= [0]
        col.landunit .= [1]

        lun = CLM.LandunitData()
        CLM.landunit_init!(lun, nl)
        lun.itype .= [CLM.ISTSOIL]

        grc_data = CLM.GridcellData()
        CLM.gridcell_init!(grc_data, ng)
        grc_data.londeg .= [0.0]

        pch = CLM.PatchData()
        CLM.patch_init!(pch, np)
        pch.column .= [1]
        pch.gridcell .= [1]
        pch.landunit .= [1]

        # Set albedos consistently
        ground_alb_d = [0.20, 0.35]
        ground_alb_i = [0.20, 0.35]
        for ib in 1:numrad
            sa.albgrd_col[1, ib] = ground_alb_d[ib]
            sa.albgri_col[1, ib] = ground_alb_i[ib]
            sa.albsod_col[1, ib] = ground_alb_d[ib]
            sa.albsoi_col[1, ib] = ground_alb_i[ib]
            sa.albsnd_hst_col[1, ib] = ground_alb_d[ib]
            sa.albsni_hst_col[1, ib] = ground_alb_i[ib]
        end

        # Simple canopy: absorbs 20% of direct, transmits 80% of direct
        for ib in 1:numrad
            sa.fabd_patch[1, ib] = 0.15
            sa.fabi_patch[1, ib] = 0.10
            sa.ftdd_patch[1, ib] = 0.70
            sa.ftid_patch[1, ib] = 0.05
            sa.ftii_patch[1, ib] = 0.80
        end

        # albd/albi should be consistent:
        # albd = fabd + (1-fabd-ftdd-ftid)*albgrd ... (simplified)
        # For a test, just set directly
        for ib in 1:numrad
            # Reflected = 1 - absorbed_canopy - absorbed_ground
            # Direct: albd = 1 - fabd - ftdd*(1-albgrd) - ftid*(1-albgri)
            sa.albd_patch[1, ib] = 1.0 - sa.fabd_patch[1, ib] -
                sa.ftdd_patch[1, ib] * (1.0 - ground_alb_d[ib]) -
                sa.ftid_patch[1, ib] * (1.0 - ground_alb_i[ib])
            # Diffuse: albi = 1 - fabi - ftii*(1-albgri)
            sa.albi_patch[1, ib] = 1.0 - sa.fabi_patch[1, ib] -
                sa.ftii_patch[1, ib] * (1.0 - ground_alb_i[ib])
        end

        forc_solad_col = zeros(nc, numrad)
        forc_solad_col[1, 1] = 500.0
        forc_solad_col[1, 2] = 250.0
        forc_solai = zeros(ng, numrad)
        forc_solai[1, 1] = 120.0
        forc_solai[1, 2] = 60.0

        total_incoming = sum(forc_solad_col) + sum(forc_solai)

        mask_nourban = trues(np)
        mask_urban = falses(np)

        CLM.surface_radiation!(sa, cs, sol, sr, wd, col, lun, grc_data, pch,
            forc_solad_col, forc_solai, mask_nourban, mask_urban, 1:np;
            dtime=3600.0, current_tod=43200.0)

        # Energy conservation: fsa + fsr ≈ total_incoming
        @test sol.fsa_patch[1] + sol.fsr_patch[1] ≈ total_incoming atol=1e-8
    end

end
