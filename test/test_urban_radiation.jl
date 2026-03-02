@testset "Urban Radiation" begin

    # ========================================================================
    # Helper: set up a minimal urban configuration with 1 landunit, 5 columns
    # (roof, sunwall, shadewall, road_imperv, road_perv), 5 patches
    # ========================================================================

    function setup_urban_radiation_test(; canyon_hwr_val=1.0,
                                          wtroad_perv_val=0.4,
                                          frac_sno_val=0.0,
                                          t_grnd_val=290.0,
                                          forc_lwrad_val=350.0,
                                          forc_solad_vis=200.0,
                                          forc_solad_nir=100.0,
                                          forc_solai_vis=80.0,
                                          forc_solai_nir=40.0)
        nl = 1  # 1 landunit
        nc = 5  # 5 columns
        np = 5  # 5 patches (one per column)
        ng = 1  # 1 gridcell
        numrad = CLM.NUMRAD

        # --- LandunitData ---
        lun = CLM.LandunitData()
        CLM.landunit_init!(lun, nl)
        lun.gridcell[1] = 1
        lun.coli[1] = 1
        lun.colf[1] = nc
        lun.canyon_hwr[1] = canyon_hwr_val
        lun.wtroad_perv[1] = wtroad_perv_val
        lun.urbpoi[1] = true
        lun.itype[1] = CLM.ISTURB_MIN

        # --- ColumnData ---
        col = CLM.ColumnData()
        CLM.column_init!(col, nc)
        col_types = [CLM.ICOL_ROOF, CLM.ICOL_SUNWALL, CLM.ICOL_SHADEWALL,
                     CLM.ICOL_ROAD_IMPERV, CLM.ICOL_ROAD_PERV]
        for c in 1:nc
            col.landunit[c] = 1
            col.itype[c] = col_types[c]
            col.snl[c] = 0
        end

        # --- PatchData ---
        pch = CLM.PatchData()
        CLM.patch_init!(pch, np)
        for p in 1:np
            pch.column[p] = p
            pch.landunit[p] = 1
            pch.gridcell[p] = 1
        end

        # --- UrbanParamsData ---
        urbanparams = CLM.UrbanParamsData()
        CLM.urbanparams_init!(urbanparams, nl; numrad=numrad)

        # Set realistic emissivities
        urbanparams.em_roof[1]    = 0.90
        urbanparams.em_improad[1] = 0.95
        urbanparams.em_perroad[1] = 0.95
        urbanparams.em_wall[1]    = 0.85

        # Set realistic albedos
        for ib in 1:numrad
            urbanparams.alb_roof_dir[1, ib]    = 0.20
            urbanparams.alb_roof_dif[1, ib]    = 0.20
            urbanparams.alb_improad_dir[1, ib] = 0.10
            urbanparams.alb_improad_dif[1, ib] = 0.10
            urbanparams.alb_perroad_dir[1, ib] = 0.15
            urbanparams.alb_perroad_dif[1, ib] = 0.15
            urbanparams.alb_wall_dir[1, ib]    = 0.25
            urbanparams.alb_wall_dif[1, ib]    = 0.25
        end

        # Compute view factors
        hwr = canyon_hwr_val
        urbanparams.vf_sr[1] = sqrt(hwr^2 + 1.0) - hwr
        urbanparams.vf_wr[1] = 0.5 * (1.0 - urbanparams.vf_sr[1])
        urbanparams.vf_sw[1] = 0.5 * (hwr + 1.0 - sqrt(hwr^2 + 1.0)) / hwr
        urbanparams.vf_rw[1] = urbanparams.vf_sw[1]
        urbanparams.vf_ww[1] = 1.0 - urbanparams.vf_sw[1] - urbanparams.vf_rw[1]

        # --- TemperatureData ---
        temperature = CLM.TemperatureData()
        CLM.temperature_init!(temperature, np, nc, nl, ng)
        temperature.t_grnd_col[1:nc] .= t_grnd_val
        temperature.t_ref2m_patch[1:np] .= t_grnd_val

        # --- SolarAbsorbedData ---
        solarabs = CLM.SolarAbsorbedData()
        CLM.solarabs_init!(solarabs, np, nl)
        CLM.solarabs_init_cold!(solarabs, 1:nl)
        # Set some realistic solar absorption values (from urban_albedo)
        for ib in 1:numrad
            solarabs.sabs_roof_dir_lun[1, ib]      = 0.80
            solarabs.sabs_roof_dif_lun[1, ib]       = 0.80
            solarabs.sabs_sunwall_dir_lun[1, ib]    = 0.30
            solarabs.sabs_sunwall_dif_lun[1, ib]    = 0.30
            solarabs.sabs_shadewall_dir_lun[1, ib]  = 0.20
            solarabs.sabs_shadewall_dif_lun[1, ib]  = 0.20
            solarabs.sabs_improad_dir_lun[1, ib]    = 0.70
            solarabs.sabs_improad_dif_lun[1, ib]    = 0.70
            solarabs.sabs_perroad_dir_lun[1, ib]    = 0.65
            solarabs.sabs_perroad_dif_lun[1, ib]    = 0.65
        end

        # --- SurfaceAlbedoData ---
        surfalb = CLM.SurfaceAlbedoData()
        CLM.surfalb_init!(surfalb, np, nc, ng)
        CLM.surfalb_init_cold!(surfalb, 1:nc, 1:np)

        # --- EnergyFluxData ---
        energyflux = CLM.EnergyFluxData()
        CLM.energyflux_init!(energyflux, np, nc, nl, ng)

        # --- WaterDiagnosticBulkData ---
        waterdiag = CLM.WaterDiagnosticBulkData()
        CLM.waterdiagnosticbulk_init!(waterdiag, nc, np, nl, ng)
        waterdiag.frac_sno_col[1:nc] .= frac_sno_val

        # --- Forcing data ---
        forc_lwrad = fill(forc_lwrad_val, ng)
        forc_solad = zeros(ng, numrad)
        forc_solad[1, 1] = forc_solad_vis
        forc_solad[1, 2] = forc_solad_nir
        forc_solai = zeros(ng, numrad)
        forc_solai[1, 1] = forc_solai_vis
        forc_solai[1, 2] = forc_solai_nir

        # --- Masks ---
        mask_nourbanl = falses(nl)
        mask_urbanl = falses(nl)
        mask_urbanl[1] = true
        mask_urbanc = falses(nc)
        mask_urbanc[1:nc] .= true
        mask_urbanp = falses(np)
        mask_urbanp[1:np] .= true

        return (lun=lun, col=col, pch=pch, urbanparams=urbanparams,
                solarabs=solarabs, surfalb=surfalb, energyflux=energyflux,
                temperature=temperature, waterdiag=waterdiag,
                forc_lwrad=forc_lwrad, forc_solad=forc_solad, forc_solai=forc_solai,
                mask_nourbanl=mask_nourbanl, mask_urbanl=mask_urbanl,
                mask_urbanc=mask_urbanc, mask_urbanp=mask_urbanp)
    end

    # ========================================================================
    @testset "net_longwave! — basic convergence" begin
        nl = 1
        mask_urbanl = trues(nl)
        hwr = 1.0
        canyon_hwr = [hwr]
        wtroad_perv = [0.4]

        urbanparams = CLM.UrbanParamsData()
        CLM.urbanparams_init!(urbanparams, nl; numrad=CLM.NUMRAD)
        urbanparams.em_roof[1]    = 0.90
        urbanparams.em_improad[1] = 0.95
        urbanparams.em_perroad[1] = 0.95
        urbanparams.em_wall[1]    = 0.85
        urbanparams.vf_sr[1] = sqrt(hwr^2 + 1.0) - hwr
        urbanparams.vf_wr[1] = 0.5 * (1.0 - urbanparams.vf_sr[1])
        urbanparams.vf_sw[1] = 0.5 * (hwr + 1.0 - sqrt(hwr^2 + 1.0)) / hwr
        urbanparams.vf_rw[1] = urbanparams.vf_sw[1]
        urbanparams.vf_ww[1] = 1.0 - urbanparams.vf_sw[1] - urbanparams.vf_rw[1]

        T = 290.0  # K
        lwdown = [350.0]
        em_roof = [0.90]
        em_improad = [0.95]
        em_perroad = [0.95]
        em_wall = urbanparams.em_wall
        t_roof = [T]
        t_improad = [T]
        t_perroad = [T]
        t_sunwall = [T]
        t_shadewall = [T]

        lwnet_roof = zeros(nl)
        lwnet_improad = zeros(nl)
        lwnet_perroad = zeros(nl)
        lwnet_sunwall = zeros(nl)
        lwnet_shadewall = zeros(nl)
        lwnet_canyon = zeros(nl)
        lwup_roof = zeros(nl)
        lwup_improad = zeros(nl)
        lwup_perroad = zeros(nl)
        lwup_sunwall = zeros(nl)
        lwup_shadewall = zeros(nl)
        lwup_canyon = zeros(nl)

        CLM.net_longwave!(canyon_hwr, wtroad_perv, lwdown,
                           em_roof, em_improad, em_perroad, em_wall,
                           t_roof, t_improad, t_perroad, t_sunwall, t_shadewall,
                           lwnet_roof, lwnet_improad, lwnet_perroad,
                           lwnet_sunwall, lwnet_shadewall, lwnet_canyon,
                           lwup_roof, lwup_improad, lwup_perroad,
                           lwup_sunwall, lwup_shadewall, lwup_canyon,
                           urbanparams, mask_urbanl, 1:nl)

        # Conservation: lwnet_canyon ≈ lwup_canyon - lwdown
        @test lwnet_canyon[1] ≈ lwup_canyon[1] - lwdown[1] atol=0.10

        # lwup_roof > 0 (emitted + reflected)
        @test lwup_roof[1] > 0.0

        # Net longwave for roof = lwup_roof - lwdown
        @test lwnet_roof[1] ≈ lwup_roof[1] - lwdown[1] atol=1e-10

        # All lwup values should be positive
        @test lwup_improad[1] > 0.0
        @test lwup_perroad[1] > 0.0
        @test lwup_sunwall[1] > 0.0
        @test lwup_shadewall[1] > 0.0
        @test lwup_canyon[1] > 0.0
    end

    # ========================================================================
    @testset "net_longwave! — energy conservation" begin
        nl = 1
        mask_urbanl = trues(nl)
        hwr = 0.5
        canyon_hwr = [hwr]
        wtroad_perv = [0.3]

        urbanparams = CLM.UrbanParamsData()
        CLM.urbanparams_init!(urbanparams, nl; numrad=CLM.NUMRAD)
        urbanparams.em_roof[1]    = 0.90
        urbanparams.em_improad[1] = 0.90
        urbanparams.em_perroad[1] = 0.90
        urbanparams.em_wall[1]    = 0.90
        urbanparams.vf_sr[1] = sqrt(hwr^2 + 1.0) - hwr
        urbanparams.vf_wr[1] = 0.5 * (1.0 - urbanparams.vf_sr[1])
        urbanparams.vf_sw[1] = 0.5 * (hwr + 1.0 - sqrt(hwr^2 + 1.0)) / hwr
        urbanparams.vf_rw[1] = urbanparams.vf_sw[1]
        urbanparams.vf_ww[1] = 1.0 - urbanparams.vf_sw[1] - urbanparams.vf_rw[1]

        T = 300.0
        lwdown = [400.0]
        em_roof = [0.90]
        em_improad = [0.90]
        em_perroad = [0.90]
        em_wall = urbanparams.em_wall
        t_roof = [T]
        t_improad = [T]
        t_perroad = [T]
        t_sunwall = [T]
        t_shadewall = [T]

        lwnet_roof = zeros(nl)
        lwnet_improad = zeros(nl)
        lwnet_perroad = zeros(nl)
        lwnet_sunwall = zeros(nl)
        lwnet_shadewall = zeros(nl)
        lwnet_canyon = zeros(nl)
        lwup_roof = zeros(nl)
        lwup_improad = zeros(nl)
        lwup_perroad = zeros(nl)
        lwup_sunwall = zeros(nl)
        lwup_shadewall = zeros(nl)
        lwup_canyon = zeros(nl)

        CLM.net_longwave!(canyon_hwr, wtroad_perv, lwdown,
                           em_roof, em_improad, em_perroad, em_wall,
                           t_roof, t_improad, t_perroad, t_sunwall, t_shadewall,
                           lwnet_roof, lwnet_improad, lwnet_perroad,
                           lwnet_sunwall, lwnet_shadewall, lwnet_canyon,
                           lwup_roof, lwup_improad, lwup_perroad,
                           lwup_sunwall, lwup_shadewall, lwup_canyon,
                           urbanparams, mask_urbanl, 1:nl)

        # Canyon energy conservation
        err = lwnet_canyon[1] - (lwup_canyon[1] - lwdown[1])
        @test abs(err) < 0.10
    end

    # ========================================================================
    @testset "urban_radiation! — basic functionality" begin
        s = setup_urban_radiation_test()

        CLM.urban_radiation!(s.solarabs, s.energyflux,
                              s.col, s.lun, s.pch,
                              s.urbanparams, s.temperature, s.waterdiag,
                              s.forc_lwrad, s.forc_solad, s.forc_solai,
                              s.mask_nourbanl, s.mask_urbanl,
                              s.mask_urbanc, s.mask_urbanp,
                              1:1, 1:5)

        # All patches should have set eflx_lwrad_out (positive)
        for p in 1:5
            @test s.energyflux.eflx_lwrad_out_patch[p] > 0.0
            @test !isnan(s.energyflux.eflx_lwrad_net_patch[p])
            @test !isnan(s.energyflux.eflx_lwrad_net_u_patch[p])
        end

        # sabg should be positive for all patches with solar radiation
        for p in 1:5
            @test s.solarabs.sabg_patch[p] > 0.0
        end

        # sabv should be zero (no vegetation in urban)
        for p in 1:5
            @test s.solarabs.sabv_patch[p] == 0.0
        end

        # fsa = sabv + sabg
        for p in 1:5
            @test s.solarabs.fsa_patch[p] ≈ s.solarabs.sabv_patch[p] + s.solarabs.sabg_patch[p]
        end

        # fsa_u = fsa
        for p in 1:5
            @test s.solarabs.fsa_u_patch[p] ≈ s.solarabs.fsa_patch[p]
        end
    end

    # ========================================================================
    @testset "urban_radiation! — non-urban SPVAL" begin
        s = setup_urban_radiation_test()

        # Create a second landunit that is non-urban
        nl = 2
        lun2 = CLM.LandunitData()
        CLM.landunit_init!(lun2, nl)
        lun2.gridcell[1] = 1
        lun2.coli[1] = 1
        lun2.colf[1] = 5
        lun2.canyon_hwr[1] = 1.0
        lun2.wtroad_perv[1] = 0.4
        lun2.urbpoi[1] = true
        lun2.itype[1] = CLM.ISTURB_MIN
        lun2.gridcell[2] = 1
        lun2.coli[2] = 6
        lun2.colf[2] = 6
        lun2.canyon_hwr[2] = 0.0
        lun2.wtroad_perv[2] = 0.0
        lun2.urbpoi[2] = false
        lun2.itype[2] = CLM.ISTSOIL

        # Extend solarabs for 2 landunits
        solarabs2 = CLM.SolarAbsorbedData()
        CLM.solarabs_init!(solarabs2, 5, nl)
        CLM.solarabs_init_cold!(solarabs2, 1:nl)
        for ib in 1:CLM.NUMRAD
            solarabs2.sabs_roof_dir_lun[1, ib]      = 0.80
            solarabs2.sabs_roof_dif_lun[1, ib]       = 0.80
            solarabs2.sabs_sunwall_dir_lun[1, ib]    = 0.30
            solarabs2.sabs_sunwall_dif_lun[1, ib]    = 0.30
            solarabs2.sabs_shadewall_dir_lun[1, ib]  = 0.20
            solarabs2.sabs_shadewall_dif_lun[1, ib]  = 0.20
            solarabs2.sabs_improad_dir_lun[1, ib]    = 0.70
            solarabs2.sabs_improad_dif_lun[1, ib]    = 0.70
            solarabs2.sabs_perroad_dir_lun[1, ib]    = 0.65
            solarabs2.sabs_perroad_dif_lun[1, ib]    = 0.65
        end

        mask_nourbanl = falses(nl)
        mask_nourbanl[2] = true
        mask_urbanl = falses(nl)
        mask_urbanl[1] = true

        # Extend urban params
        urbanparams2 = CLM.UrbanParamsData()
        CLM.urbanparams_init!(urbanparams2, nl; numrad=CLM.NUMRAD)
        urbanparams2.em_roof[1]    = 0.90
        urbanparams2.em_improad[1] = 0.95
        urbanparams2.em_perroad[1] = 0.95
        urbanparams2.em_wall[1]    = 0.85
        hwr = 1.0
        urbanparams2.vf_sr[1] = sqrt(hwr^2 + 1.0) - hwr
        urbanparams2.vf_wr[1] = 0.5 * (1.0 - urbanparams2.vf_sr[1])
        urbanparams2.vf_sw[1] = 0.5 * (hwr + 1.0 - sqrt(hwr^2 + 1.0)) / hwr
        urbanparams2.vf_rw[1] = urbanparams2.vf_sw[1]
        urbanparams2.vf_ww[1] = 1.0 - urbanparams2.vf_sw[1] - urbanparams2.vf_rw[1]

        CLM.urban_radiation!(solarabs2, s.energyflux,
                              s.col, lun2, s.pch,
                              urbanparams2, s.temperature, s.waterdiag,
                              s.forc_lwrad, s.forc_solad, s.forc_solai,
                              mask_nourbanl, mask_urbanl,
                              s.mask_urbanc, s.mask_urbanp,
                              1:nl, 1:5)

        # Non-urban landunit should have SPVAL for sabs fields
        for ib in 1:CLM.NUMRAD
            @test solarabs2.sabs_roof_dir_lun[2, ib] == CLM.SPVAL
            @test solarabs2.sabs_roof_dif_lun[2, ib] == CLM.SPVAL
        end
    end

    # ========================================================================
    @testset "urban_radiation! — with snow" begin
        s = setup_urban_radiation_test(frac_sno_val=0.5, t_grnd_val=265.0)

        CLM.urban_radiation!(s.solarabs, s.energyflux,
                              s.col, s.lun, s.pch,
                              s.urbanparams, s.temperature, s.waterdiag,
                              s.forc_lwrad, s.forc_solad, s.forc_solai,
                              s.mask_nourbanl, s.mask_urbanl,
                              s.mask_urbanc, s.mask_urbanp,
                              1:1, 1:5)

        # All outputs should still be valid
        for p in 1:5
            @test s.energyflux.eflx_lwrad_out_patch[p] > 0.0
            @test !isnan(s.energyflux.eflx_lwrad_net_patch[p])
            @test s.solarabs.sabg_patch[p] > 0.0
        end
    end

    # ========================================================================
    @testset "urban_radiation! — roof longwave consistency" begin
        s = setup_urban_radiation_test(t_grnd_val=300.0, forc_lwrad_val=300.0)

        CLM.urban_radiation!(s.solarabs, s.energyflux,
                              s.col, s.lun, s.pch,
                              s.urbanparams, s.temperature, s.waterdiag,
                              s.forc_lwrad, s.forc_solad, s.forc_solai,
                              s.mask_nourbanl, s.mask_urbanl,
                              s.mask_urbanc, s.mask_urbanp,
                              1:1, 1:5)

        # Roof patch (p=1): lwup = em*sb*T^4 + (1-em)*lwdown
        T = 300.0
        em = 0.90  # no snow
        expected_lwup = em * CLM.SB * T^4 + (1.0 - em) * 300.0
        @test s.energyflux.eflx_lwrad_out_patch[1] ≈ expected_lwup atol=1e-6

        # Net = lwup - lwdown
        expected_lwnet = expected_lwup - 300.0
        @test s.energyflux.eflx_lwrad_net_patch[1] ≈ expected_lwnet atol=1e-6
    end

    # ========================================================================
    @testset "Constants defined" begin
        @test CLM.MPE_URBAN_RAD == 1.0e-06
        @test CLM.SNOEM == 0.97
    end

end
