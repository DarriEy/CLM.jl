@testset "Urban Albedo" begin

    # ========================================================================
    # Helper: set up a minimal urban configuration with 1 landunit, 5 columns
    # (roof, sunwall, shadewall, road_imperv, road_perv), 5 patches
    # ========================================================================

    function setup_urban_test(; coszen_val=0.5, canyon_hwr_val=1.0,
                                wtroad_perv_val=0.4, frac_sno_val=0.0,
                                h2osno_val=0.0)
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
            col.snl[c] = 0  # no snow layers
        end

        # --- PatchData ---
        pch = CLM.PatchData()
        CLM.patch_init!(pch, np)
        for p in 1:np
            pch.column[p] = p    # 1:1 mapping
            pch.landunit[p] = 1
        end

        # --- UrbanParamsData ---
        urbanparams = CLM.UrbanParamsData()
        CLM.urbanparams_init!(urbanparams, nl; numrad=numrad)

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

        # Compute view factors from canyon_hwr (same formula as urbanparams_populate!)
        hwr = canyon_hwr_val
        urbanparams.vf_sr[1] = sqrt(hwr^2 + 1.0) - hwr
        urbanparams.vf_wr[1] = 0.5 * (1.0 - urbanparams.vf_sr[1])
        urbanparams.vf_sw[1] = 0.5 * (hwr + 1.0 - sqrt(hwr^2 + 1.0)) / hwr
        urbanparams.vf_rw[1] = urbanparams.vf_sw[1]
        urbanparams.vf_ww[1] = 1.0 - urbanparams.vf_sw[1] - urbanparams.vf_rw[1]

        # --- SolarAbsorbedData ---
        solarabs = CLM.SolarAbsorbedData()
        CLM.solarabs_init!(solarabs, np, nl)
        CLM.solarabs_init_cold!(solarabs, 1:nl)

        # --- SurfaceAlbedoData ---
        surfalb = CLM.SurfaceAlbedoData()
        CLM.surfalb_init!(surfalb, np, nc, ng)
        CLM.surfalb_init_cold!(surfalb, 1:nc, 1:np)
        surfalb.coszen_col[1:nc] .= coszen_val

        # --- WaterStateBulkData ---
        waterstatebulk = CLM.WaterStateBulkData()
        CLM.waterstatebulk_init!(waterstatebulk, nc, np, nl, ng)
        waterstatebulk.ws.h2osno_no_layers_col[1:nc] .= h2osno_val
        for c in 1:nc
            nlevtot = size(waterstatebulk.ws.h2osoi_ice_col, 2)
            for j in 1:nlevtot
                waterstatebulk.ws.h2osoi_ice_col[c, j] = 0.0
                waterstatebulk.ws.h2osoi_liq_col[c, j] = 0.0
            end
        end

        # --- WaterDiagnosticBulkData ---
        waterdiagbulk = CLM.WaterDiagnosticBulkData()
        CLM.waterdiagnosticbulk_init!(waterdiagbulk, nc, np, nl, ng)
        waterdiagbulk.frac_sno_col[1:nc] .= frac_sno_val

        # --- Masks ---
        mask_urbanl = falses(nl)
        mask_urbanl[1] = true
        mask_urbanc = falses(nc)
        mask_urbanc[1:nc] .= true
        mask_urbanp = falses(np)
        mask_urbanp[1:np] .= true

        return (lun=lun, col=col, pch=pch, urbanparams=urbanparams,
                solarabs=solarabs, surfalb=surfalb,
                waterstatebulk=waterstatebulk, waterdiagbulk=waterdiagbulk,
                mask_urbanl=mask_urbanl, mask_urbanc=mask_urbanc,
                mask_urbanp=mask_urbanp)
    end

    # ========================================================================
    @testset "Snow albedo constants" begin
        @test CLM.SNAL0 == 0.66
        @test CLM.SNAL1 == 0.56
    end

    # ========================================================================
    @testset "snow_albedo! — no snow" begin
        s = setup_urban_test(coszen_val=0.5, h2osno_val=0.0)
        nl = 1; numrad = CLM.NUMRAD
        albsn_roof    = zeros(nl, numrad)
        albsn_improad = zeros(nl, numrad)
        albsn_perroad = zeros(nl, numrad)
        h2osno_total  = zeros(length(s.mask_urbanc))

        coszen = [0.5]
        CLM.snow_albedo!(s.mask_urbanc, s.col, coszen, 0,
                         albsn_roof, albsn_improad, albsn_perroad,
                         h2osno_total)

        @test all(albsn_roof .== 0.0)
        @test all(albsn_improad .== 0.0)
        @test all(albsn_perroad .== 0.0)
    end

    # ========================================================================
    @testset "snow_albedo! — with snow" begin
        s = setup_urban_test(coszen_val=0.5, h2osno_val=10.0)
        nl = 1; numrad = CLM.NUMRAD
        albsn_roof    = zeros(nl, numrad)
        albsn_improad = zeros(nl, numrad)
        albsn_perroad = zeros(nl, numrad)
        h2osno_total  = fill(10.0, length(s.mask_urbanc))

        coszen = [0.5]
        CLM.snow_albedo!(s.mask_urbanc, s.col, coszen, 0,
                         albsn_roof, albsn_improad, albsn_perroad,
                         h2osno_total)

        @test albsn_roof[1, 1] ≈ CLM.SNAL0
        @test albsn_roof[1, 2] ≈ CLM.SNAL1
        @test albsn_improad[1, 1] ≈ CLM.SNAL0
        @test albsn_improad[1, 2] ≈ CLM.SNAL1
        @test albsn_perroad[1, 1] ≈ CLM.SNAL0
        @test albsn_perroad[1, 2] ≈ CLM.SNAL1
    end

    # ========================================================================
    @testset "incident_direct! — conservation" begin
        nl = 1; numrad = CLM.NUMRAD
        mask_urbanl = trues(nl)
        canyon_hwr = [1.0]
        coszen_val = 0.5
        coszen = [coszen_val]
        zen = [acos(coszen_val)]
        sdir = ones(nl, numrad)
        sdir_road = zeros(nl, numrad)
        sdir_sunwall = zeros(nl, numrad)
        sdir_shadewall = zeros(nl, numrad)

        CLM.incident_direct!(mask_urbanl, canyon_hwr, coszen, zen,
                             sdir, sdir_road, sdir_sunwall, sdir_shadewall)

        for ib in 1:numrad
            # Conservation: sdir = sdir_road + (sdir_sunwall + sdir_shadewall) * canyon_hwr
            total = sdir_road[1, ib] + (sdir_sunwall[1, ib] + sdir_shadewall[1, ib]) * canyon_hwr[1]
            @test total ≈ sdir[1, ib] atol=0.001
            @test sdir_shadewall[1, ib] == 0.0
            @test sdir_road[1, ib] >= 0.0
            @test sdir_sunwall[1, ib] >= 0.0
        end
    end

    # ========================================================================
    @testset "incident_direct! — night (coszen=0)" begin
        nl = 1; numrad = CLM.NUMRAD
        mask_urbanl = trues(nl)
        canyon_hwr = [1.0]
        coszen = [0.0]
        zen = [acos(0.0)]
        sdir = ones(nl, numrad)
        sdir_road = zeros(nl, numrad)
        sdir_sunwall = zeros(nl, numrad)
        sdir_shadewall = zeros(nl, numrad)

        CLM.incident_direct!(mask_urbanl, canyon_hwr, coszen, zen,
                             sdir, sdir_road, sdir_sunwall, sdir_shadewall)

        for ib in 1:numrad
            @test sdir_road[1, ib] == 0.0
            @test sdir_sunwall[1, ib] == 0.0
            @test sdir_shadewall[1, ib] == 0.0
        end
    end

    # ========================================================================
    @testset "incident_diffuse! — conservation" begin
        nl = 1; numrad = CLM.NUMRAD
        mask_urbanl = trues(nl)
        canyon_hwr = [1.0]
        sdif = ones(nl, numrad)
        sdif_road = zeros(nl, numrad)
        sdif_sunwall = zeros(nl, numrad)
        sdif_shadewall = zeros(nl, numrad)

        urbanparams = CLM.UrbanParamsData()
        CLM.urbanparams_init!(urbanparams, nl; numrad=numrad)
        hwr = 1.0
        urbanparams.vf_sr[1] = sqrt(hwr^2 + 1.0) - hwr
        urbanparams.vf_wr[1] = 0.5 * (1.0 - urbanparams.vf_sr[1])
        urbanparams.vf_sw[1] = 0.5 * (hwr + 1.0 - sqrt(hwr^2 + 1.0)) / hwr
        urbanparams.vf_rw[1] = urbanparams.vf_sw[1]
        urbanparams.vf_ww[1] = 1.0 - urbanparams.vf_sw[1] - urbanparams.vf_rw[1]

        CLM.incident_diffuse!(mask_urbanl, canyon_hwr, sdif,
                              sdif_road, sdif_sunwall, sdif_shadewall,
                              urbanparams)

        for ib in 1:numrad
            total = sdif_road[1, ib] + (sdif_sunwall[1, ib] + sdif_shadewall[1, ib]) * canyon_hwr[1]
            @test total ≈ sdif[1, ib] atol=0.001
        end
    end

    # ========================================================================
    @testset "urban_albedo! — full integration, sun up, no snow" begin
        s = setup_urban_test(coszen_val=0.5, canyon_hwr_val=1.0,
                             wtroad_perv_val=0.4, frac_sno_val=0.0,
                             h2osno_val=0.0)

        CLM.urban_albedo!(s.mask_urbanl, s.mask_urbanc, s.mask_urbanp,
                          s.lun, s.col, s.pch,
                          s.waterstatebulk, s.waterdiagbulk,
                          s.urbanparams, s.solarabs, s.surfalb)

        numrad = CLM.NUMRAD
        # All ground albedos should be set (not NaN, not zero for roof)
        for ib in 1:numrad
            for c in 1:5
                @test !isnan(s.surfalb.albgrd_col[c, ib])
                @test !isnan(s.surfalb.albgri_col[c, ib])
                @test s.surfalb.albgrd_col[c, ib] >= 0.0
                @test s.surfalb.albgri_col[c, ib] >= 0.0
                @test s.surfalb.albgrd_col[c, ib] <= 1.0
                @test s.surfalb.albgri_col[c, ib] <= 1.0
            end
        end

        # Patch albedos should match column ground albedos
        for ib in 1:numrad
            for p in 1:5
                c = s.pch.column[p]
                @test s.surfalb.albd_patch[p, ib] ≈ s.surfalb.albgrd_col[c, ib]
                @test s.surfalb.albi_patch[p, ib] ≈ s.surfalb.albgri_col[c, ib]
            end
        end

        # Canopy fluxes: no vegetation in urban, fabd=fabi=0
        for ib in 1:numrad
            for p in 1:5
                @test s.surfalb.fabd_patch[p, ib] == 0.0
                @test s.surfalb.fabi_patch[p, ib] == 0.0
                @test s.surfalb.ftdd_patch[p, ib] == 1.0
                @test s.surfalb.ftii_patch[p, ib] == 1.0
            end
        end

        # Solar absorption should be set
        for ib in 1:numrad
            @test s.solarabs.sabs_roof_dir_lun[1, ib] >= 0.0
            @test s.solarabs.sabs_roof_dif_lun[1, ib] >= 0.0
        end

        # Roof albedo check: for no snow, roof albedo should equal material albedo
        # albgrd for roof column (c=1) = sref_roof_dir = alb_roof_dir * sdir
        # With sdir=1 and alb_roof_dir=0.20, sref_roof_dir = 0.20
        @test s.surfalb.albgrd_col[1, 1] ≈ 0.20 atol=1e-10
        @test s.surfalb.albgri_col[1, 1] ≈ 0.20 atol=1e-10
    end

    # ========================================================================
    @testset "urban_albedo! — night (coszen ≤ 0)" begin
        s = setup_urban_test(coszen_val=0.0)

        CLM.urban_albedo!(s.mask_urbanl, s.mask_urbanc, s.mask_urbanp,
                          s.lun, s.col, s.pch,
                          s.waterstatebulk, s.waterdiagbulk,
                          s.urbanparams, s.solarabs, s.surfalb)

        numrad = CLM.NUMRAD
        # At night, all albedos should be initialized to defaults (view factor values)
        for ib in 1:numrad
            for c in 1:5
                @test s.surfalb.albgrd_col[c, ib] == 0.0
                @test s.surfalb.albgri_col[c, ib] == 0.0
            end
        end
    end

    # ========================================================================
    @testset "urban_albedo! — with snow" begin
        s = setup_urban_test(coszen_val=0.5, frac_sno_val=0.5,
                             h2osno_val=10.0)

        CLM.urban_albedo!(s.mask_urbanl, s.mask_urbanc, s.mask_urbanp,
                          s.lun, s.col, s.pch,
                          s.waterstatebulk, s.waterdiagbulk,
                          s.urbanparams, s.solarabs, s.surfalb)

        # Roof albedo with 50% snow should be between snow-free and full-snow values
        # snow-free roof albedo = 0.20, snow albedo vis = 0.66
        # expected = 0.20*0.5 + 0.66*0.5 = 0.43 (for roof direct vis)
        # but after net_solar processing, it becomes sref_roof_dir = alb_roof_dir_s * sdir
        roof_alb = s.surfalb.albgrd_col[1, 1]
        @test roof_alb > 0.20   # higher than snow-free
        @test roof_alb < 0.66   # lower than full snow
    end

    # ========================================================================
    @testset "urban_albedo! — varying H/W ratio" begin
        # Wider canyon (lower H/W) should let more direct light reach the road
        s_wide = setup_urban_test(coszen_val=0.5, canyon_hwr_val=0.5)
        CLM.urban_albedo!(s_wide.mask_urbanl, s_wide.mask_urbanc, s_wide.mask_urbanp,
                          s_wide.lun, s_wide.col, s_wide.pch,
                          s_wide.waterstatebulk, s_wide.waterdiagbulk,
                          s_wide.urbanparams, s_wide.solarabs, s_wide.surfalb)

        s_narrow = setup_urban_test(coszen_val=0.5, canyon_hwr_val=2.0)
        CLM.urban_albedo!(s_narrow.mask_urbanl, s_narrow.mask_urbanc, s_narrow.mask_urbanp,
                          s_narrow.lun, s_narrow.col, s_narrow.pch,
                          s_narrow.waterstatebulk, s_narrow.waterdiagbulk,
                          s_narrow.urbanparams, s_narrow.solarabs, s_narrow.surfalb)

        # Narrower canyon traps more radiation (lower canyon albedo).
        # The impervious road column (c=4) should have lower albedo in narrow canyon
        # because multiple reflections within the canyon absorb more energy.
        @test s_narrow.surfalb.albgrd_col[4, 1] < s_wide.surfalb.albgrd_col[4, 1]
    end

    # ========================================================================
    @testset "Energy conservation in urban canyon" begin
        s = setup_urban_test(coszen_val=0.5, canyon_hwr_val=1.0)
        CLM.urban_albedo!(s.mask_urbanl, s.mask_urbanc, s.mask_urbanp,
                          s.lun, s.col, s.pch,
                          s.waterstatebulk, s.waterdiagbulk,
                          s.urbanparams, s.solarabs, s.surfalb)

        sa = s.solarabs
        # For each waveband, check that absorbed + reflected ≈ incident
        # Roof: sabs + sref ≈ sdir + sdif = 1 + 1 = 2 (per unit flux, sdir=sdif=1)
        for ib in 1:CLM.NUMRAD
            roof_abs = sa.sabs_roof_dir_lun[1, ib] + sa.sabs_roof_dif_lun[1, ib]
            # roof reflected = albgrd_col[roof_col] + albgri_col[roof_col]
            roof_ref_dir = s.surfalb.albgrd_col[1, ib]
            roof_ref_dif = s.surfalb.albgri_col[1, ib]
            @test roof_abs + roof_ref_dir + roof_ref_dif ≈ 2.0 atol=0.01
        end
    end

end
