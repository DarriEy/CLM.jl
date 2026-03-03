@testset "Satellite Phenology (SatellitePhenologyMod)" begin

    # -----------------------------------------------------------------------
    # Test SatellitePhenologyData construction and init/clean
    # -----------------------------------------------------------------------
    @testset "SatellitePhenologyData construction and init/clean" begin
        sp = CLM.SatellitePhenologyData()
        @test sp.InterpMonths1 == -999
        @test length(sp.timwt) == 2
        @test size(sp.mlai2t) == (0, 0)
        @test size(sp.msai2t) == (0, 0)
        @test size(sp.mhvt2t) == (0, 0)
        @test size(sp.mhvb2t) == (0, 0)

        np = 6
        CLM.satellite_phenology_init!(sp, np)

        @test sp.InterpMonths1 == -999
        @test size(sp.mlai2t) == (np, 2)
        @test size(sp.msai2t) == (np, 2)
        @test size(sp.mhvt2t) == (np, 2)
        @test size(sp.mhvb2t) == (np, 2)
        @test all(isnan, sp.mlai2t)
        @test all(isnan, sp.msai2t)
        @test all(isnan, sp.mhvt2t)
        @test all(isnan, sp.mhvb2t)

        # Clean
        CLM.satellite_phenology_clean!(sp)
        @test sp.InterpMonths1 == -999
        @test size(sp.mlai2t) == (0, 0)
        @test size(sp.msai2t) == (0, 0)
        @test size(sp.mhvt2t) == (0, 0)
        @test size(sp.mhvb2t) == (0, 0)
    end

    # -----------------------------------------------------------------------
    # Test interp_monthly_veg! time weight computation
    # -----------------------------------------------------------------------
    @testset "interp_monthly_veg! time weights" begin
        sp = CLM.SatellitePhenologyData()
        CLM.satellite_phenology_init!(sp, 2)

        # Mid-January (day 16 of 31-day month)
        months, needs_read = CLM.interp_monthly_veg!(sp; kmo=1, kda=16)
        @test needs_read == true
        @test sp.timwt[1] + sp.timwt[2] ≈ 1.0
        @test sp.timwt[1] >= 0.0
        @test sp.timwt[2] >= 0.0

        # Same call again should not need read
        months2, needs_read2 = CLM.interp_monthly_veg!(sp; kmo=1, kda=17)
        if months2[1] == months[1]
            @test needs_read2 == false
        end

        # December wrapping to January
        sp2 = CLM.SatellitePhenologyData()
        CLM.satellite_phenology_init!(sp2, 2)
        months_dec, _ = CLM.interp_monthly_veg!(sp2; kmo=12, kda=25)
        # Late December: month2 should wrap to 1 (January)
        @test months_dec[2] == 1 || months_dec[2] == 12

        # Early January: month1 could wrap to 12 (December)
        sp3 = CLM.SatellitePhenologyData()
        CLM.satellite_phenology_init!(sp3, 2)
        months_jan, _ = CLM.interp_monthly_veg!(sp3; kmo=1, kda=3)
        @test months_jan[1] == 12 || months_jan[1] == 1

        # Check time weights always sum to 1 for various days
        for kmo in 1:12
            for kda in [1, 10, 15, 20, CLM.NDAYPM[kmo]]
                sp_t = CLM.SatellitePhenologyData()
                CLM.satellite_phenology_init!(sp_t, 1)
                CLM.interp_monthly_veg!(sp_t; kmo=kmo, kda=kda)
                @test sp_t.timwt[1] + sp_t.timwt[2] ≈ 1.0 atol=1e-15
                @test sp_t.timwt[1] >= 0.0
                @test sp_t.timwt[2] >= 0.0
            end
        end
    end

    # -----------------------------------------------------------------------
    # Test satellite_phenology! — core physics
    # -----------------------------------------------------------------------
    @testset "satellite_phenology! core physics" begin
        np = 4   # 4 patches
        nc = 2   # 2 columns

        # Set up SatellitePhenologyData
        sp = CLM.SatellitePhenologyData()
        CLM.satellite_phenology_init!(sp, np)
        sp.timwt[1] = 0.6
        sp.timwt[2] = 0.4

        # Set monthly values: month 1 and month 2
        # Patch 1: non-vegetated (noveg=0)
        # Patch 2: needleleaf tree (itype=1, > noveg, <= nbrdlf_dcd_brl_shrub=11)
        # Patch 3: broadleaf deciduous boreal shrub (itype=11, <= 11)
        # Patch 4: c3 arctic grass (itype=12, > 11 → uses alternative snow burial)
        sp.mlai2t[1, :] .= [0.0, 0.0]
        sp.mlai2t[2, :] .= [3.0, 4.0]
        sp.mlai2t[3, :] .= [1.5, 2.0]
        sp.mlai2t[4, :] .= [2.0, 2.5]

        sp.msai2t[1, :] .= [0.0, 0.0]
        sp.msai2t[2, :] .= [0.5, 0.6]
        sp.msai2t[3, :] .= [0.3, 0.4]
        sp.msai2t[4, :] .= [0.4, 0.5]

        sp.mhvt2t[1, :] .= [0.0, 0.0]
        sp.mhvt2t[2, :] .= [17.0, 18.0]
        sp.mhvt2t[3, :] .= [0.5, 0.6]
        sp.mhvt2t[4, :] .= [0.3, 0.4]

        sp.mhvb2t[1, :] .= [0.0, 0.0]
        sp.mhvb2t[2, :] .= [8.0, 9.0]
        sp.mhvb2t[3, :] .= [0.1, 0.15]
        sp.mhvb2t[4, :] .= [0.01, 0.02]

        # Set up CanopyStateData
        cs = CLM.CanopyStateData()
        CLM.canopystate_init!(cs, np)

        # Set up WaterDiagnosticBulkData (need frac_sno_col, snow_depth_col)
        wd = CLM.WaterDiagnosticBulkData()
        wd.frac_sno_col = [0.0, 0.3]        # col 1: no snow, col 2: partial snow
        wd.snow_depth_col = [0.0, 0.15]      # col 1: no snow, col 2: 15cm snow

        # Set up PatchData
        patch = CLM.PatchData()
        CLM.patch_init!(patch, np)
        patch.itype  .= [0, 1, 11, 12]
        patch.column .= [1, 1, 2, 2]         # patches 1,2→col 1; patches 3,4→col 2

        mask_patch = trues(np)
        bounds_patch = 1:np

        CLM.satellite_phenology!(sp, cs, wd, patch, mask_patch, bounds_patch;
                                  noveg=0, nbrdlf_dcd_brl_shrub=11,
                                  use_lai_streams=false, use_fates_sp=false)

        # --- Check interpolated values ---
        # tlai = timwt[1]*mlai2t[:,1] + timwt[2]*mlai2t[:,2]
        @test cs.tlai_patch[1] ≈ 0.6 * 0.0 + 0.4 * 0.0
        @test cs.tlai_patch[2] ≈ 0.6 * 3.0 + 0.4 * 4.0
        @test cs.tlai_patch[3] ≈ 0.6 * 1.5 + 0.4 * 2.0
        @test cs.tlai_patch[4] ≈ 0.6 * 2.0 + 0.4 * 2.5

        # tsai
        @test cs.tsai_patch[2] ≈ 0.6 * 0.5 + 0.4 * 0.6

        # htop
        @test cs.htop_patch[2] ≈ 0.6 * 17.0 + 0.4 * 18.0

        # hbot
        @test cs.hbot_patch[2] ≈ 0.6 * 8.0 + 0.4 * 9.0

        # --- Check snow burial logic ---

        # Patch 1 (noveg): itype=0, not > noveg, so uses the else branch
        # snow_depth=0 for col 1, so fb ≈ 1
        # tlai=0, tsai=0 → elai=0, esai=0, frac_veg_nosno_alb=0
        @test cs.elai_patch[1] == 0.0
        @test cs.esai_patch[1] == 0.0
        @test cs.frac_veg_nosno_alb_patch[1] == 0

        # Patch 2 (itype=1, col 1, no snow): fb ≈ 1, frac_sno=0
        # elai = max(tlai*(1-0) + tlai*1*0, 0) = tlai
        @test cs.elai_patch[2] ≈ cs.tlai_patch[2]
        @test cs.esai_patch[2] ≈ cs.tsai_patch[2]
        @test cs.frac_veg_nosno_alb_patch[2] == 1

        # Patch 3 (itype=11, col 2, partial snow):
        # Uses the first branch (itype > noveg && itype <= nbrdlf_dcd_brl_shrub)
        htop3 = cs.htop_patch[3]
        hbot3 = cs.hbot_patch[3]
        snow_d = 0.15
        ol3 = min(max(snow_d - hbot3, 0.0), htop3 - hbot3)
        fb3 = 1.0 - ol3 / max(1.0e-06, htop3 - hbot3)
        frac_sno3 = 0.3
        elai3_expected = max(cs.tlai_patch[3] * (1.0 - frac_sno3) + cs.tlai_patch[3] * fb3 * frac_sno3, 0.0)
        @test cs.elai_patch[3] ≈ elai3_expected

        # Patch 4 (itype=12, col 2, partial snow):
        # Uses the else branch (itype > nbrdlf_dcd_brl_shrub)
        htop4 = cs.htop_patch[4]
        fb4 = 1.0 - (max(min(snow_d, max(0.05, htop4 * 0.8)), 0.0) /
                      max(0.05, htop4 * 0.8))
        elai4_expected = max(cs.tlai_patch[4] * (1.0 - frac_sno3) + cs.tlai_patch[4] * fb4 * frac_sno3, 0.0)
        if elai4_expected < 0.05
            elai4_expected = 0.0
        end
        @test cs.elai_patch[4] ≈ elai4_expected
    end

    # -----------------------------------------------------------------------
    # Test satellite_phenology! — fates_sp mode skips elai/esai update
    # -----------------------------------------------------------------------
    @testset "satellite_phenology! use_fates_sp=true skips elai/esai" begin
        np = 2
        nc = 1

        sp = CLM.SatellitePhenologyData()
        CLM.satellite_phenology_init!(sp, np)
        sp.timwt[1] = 0.5
        sp.timwt[2] = 0.5
        sp.mlai2t .= 3.0
        sp.msai2t .= 0.5
        sp.mhvt2t .= 10.0
        sp.mhvb2t .= 2.0

        cs = CLM.CanopyStateData()
        CLM.canopystate_init!(cs, np)
        # Pre-set elai/esai to sentinel values
        cs.elai_patch .= -999.0
        cs.esai_patch .= -999.0
        cs.frac_veg_nosno_alb_patch .= -1

        wd = CLM.WaterDiagnosticBulkData()
        wd.frac_sno_col = [0.2]
        wd.snow_depth_col = [0.1]

        patch = CLM.PatchData()
        CLM.patch_init!(patch, np)
        patch.itype  .= [1, 5]
        patch.column .= [1, 1]

        mask_patch = trues(np)

        CLM.satellite_phenology!(sp, cs, wd, patch, mask_patch, 1:np;
                                  use_fates_sp=true)

        # tlai, tsai, htop, hbot should be updated
        @test cs.tlai_patch[1] ≈ 3.0
        @test cs.tsai_patch[1] ≈ 0.5
        @test cs.htop_patch[1] ≈ 10.0
        @test cs.hbot_patch[1] ≈ 2.0

        # elai, esai, frac_veg_nosno_alb should NOT be updated (stay at sentinel)
        @test cs.elai_patch[1] == -999.0
        @test cs.esai_patch[1] == -999.0
        @test cs.frac_veg_nosno_alb_patch[1] == -1
    end

    # -----------------------------------------------------------------------
    # Test satellite_phenology! — use_lai_streams skips tlai interpolation
    # -----------------------------------------------------------------------
    @testset "satellite_phenology! use_lai_streams=true skips tlai" begin
        np = 2
        nc = 1

        sp = CLM.SatellitePhenologyData()
        CLM.satellite_phenology_init!(sp, np)
        sp.timwt[1] = 0.5
        sp.timwt[2] = 0.5
        sp.mlai2t .= 3.0
        sp.msai2t .= 0.5
        sp.mhvt2t .= 10.0
        sp.mhvb2t .= 2.0

        cs = CLM.CanopyStateData()
        CLM.canopystate_init!(cs, np)
        cs.tlai_patch .= 7.77  # pre-set to something

        wd = CLM.WaterDiagnosticBulkData()
        wd.frac_sno_col = [0.0]
        wd.snow_depth_col = [0.0]

        patch = CLM.PatchData()
        CLM.patch_init!(patch, np)
        patch.itype  .= [1, 5]
        patch.column .= [1, 1]

        CLM.satellite_phenology!(sp, cs, wd, patch, trues(np), 1:np;
                                  use_lai_streams=true)

        # tlai should NOT be overwritten
        @test cs.tlai_patch[1] ≈ 7.77
        @test cs.tlai_patch[2] ≈ 7.77

        # tsai, htop, hbot SHOULD be updated
        @test cs.tsai_patch[1] ≈ 0.5
        @test cs.htop_patch[1] ≈ 10.0
    end

    # -----------------------------------------------------------------------
    # Test satellite_phenology! — small elai/esai threshold
    # -----------------------------------------------------------------------
    @testset "satellite_phenology! small LAI threshold" begin
        np = 1
        nc = 1

        sp = CLM.SatellitePhenologyData()
        CLM.satellite_phenology_init!(sp, np)
        sp.timwt[1] = 1.0
        sp.timwt[2] = 0.0
        # Very small LAI and SAI
        sp.mlai2t[1, :] .= [0.03, 0.0]
        sp.msai2t[1, :] .= [0.02, 0.0]
        sp.mhvt2t[1, :] .= [10.0, 10.0]
        sp.mhvb2t[1, :] .= [2.0, 2.0]

        cs = CLM.CanopyStateData()
        CLM.canopystate_init!(cs, np)

        wd = CLM.WaterDiagnosticBulkData()
        wd.frac_sno_col = [0.0]
        wd.snow_depth_col = [0.0]

        patch = CLM.PatchData()
        CLM.patch_init!(patch, np)
        patch.itype .= [1]
        patch.column .= [1]

        CLM.satellite_phenology!(sp, cs, wd, patch, trues(np), 1:np)

        # tlai = 0.03, but elai < 0.05 should be set to 0
        @test cs.tlai_patch[1] ≈ 0.03
        @test cs.elai_patch[1] == 0.0
        @test cs.esai_patch[1] == 0.0
        @test cs.frac_veg_nosno_alb_patch[1] == 0
    end

    # -----------------------------------------------------------------------
    # Test satellite_phenology! — mask filtering
    # -----------------------------------------------------------------------
    @testset "satellite_phenology! respects mask" begin
        np = 2
        nc = 1

        sp = CLM.SatellitePhenologyData()
        CLM.satellite_phenology_init!(sp, np)
        sp.timwt[1] = 0.5
        sp.timwt[2] = 0.5
        sp.mlai2t .= 5.0
        sp.msai2t .= 1.0
        sp.mhvt2t .= 15.0
        sp.mhvb2t .= 3.0

        cs = CLM.CanopyStateData()
        CLM.canopystate_init!(cs, np)
        cs.tlai_patch .= -1.0  # sentinel

        wd = CLM.WaterDiagnosticBulkData()
        wd.frac_sno_col = [0.0]
        wd.snow_depth_col = [0.0]

        patch = CLM.PatchData()
        CLM.patch_init!(patch, np)
        patch.itype  .= [1, 5]
        patch.column .= [1, 1]

        # Only patch 1 is active
        mask = BitVector([true, false])

        CLM.satellite_phenology!(sp, cs, wd, patch, mask, 1:np)

        @test cs.tlai_patch[1] ≈ 5.0   # updated
        @test cs.tlai_patch[2] ≈ -1.0   # NOT updated (masked out)
    end

    # -----------------------------------------------------------------------
    # Test read_monthly_vegetation!
    # -----------------------------------------------------------------------
    @testset "read_monthly_vegetation!" begin
        np = 3
        ng = 2
        maxveg = 14

        sp = CLM.SatellitePhenologyData()
        CLM.satellite_phenology_init!(sp, np)

        cs = CLM.CanopyStateData()
        CLM.canopystate_init!(cs, np)

        patch = CLM.PatchData()
        CLM.patch_init!(patch, np)
        patch.itype    .= [0, 1, 12]        # noveg, needleleaf, c3_arctic_grass
        patch.gridcell .= [1, 1, 2]

        # Create mock monthly data: [gridcell, pft+1, month]
        monthly_lai  = zeros(ng, maxveg + 1, 12)
        monthly_sai  = zeros(ng, maxveg + 1, 12)
        monthly_htop = zeros(ng, maxveg + 1, 12)
        monthly_hbot = zeros(ng, maxveg + 1, 12)

        # Set data for PFT 1 (index 2), gridcell 1, months 6,7
        monthly_lai[1, 2, 6] = 3.5
        monthly_lai[1, 2, 7] = 4.2
        monthly_sai[1, 2, 6] = 0.5
        monthly_sai[1, 2, 7] = 0.6
        monthly_htop[1, 2, 6] = 15.0
        monthly_htop[1, 2, 7] = 16.0
        monthly_hbot[1, 2, 6] = 5.0
        monthly_hbot[1, 2, 7] = 5.5

        # Set data for PFT 12 (index 13), gridcell 2, months 6,7
        monthly_lai[2, 13, 6] = 1.0
        monthly_lai[2, 13, 7] = 1.2
        monthly_sai[2, 13, 6] = 0.2
        monthly_sai[2, 13, 7] = 0.3
        monthly_htop[2, 13, 6] = 0.5
        monthly_htop[2, 13, 7] = 0.6
        monthly_hbot[2, 13, 6] = 0.05
        monthly_hbot[2, 13, 7] = 0.06

        CLM.read_monthly_vegetation!(sp, cs, patch, 1:np;
                                      monthly_lai=monthly_lai,
                                      monthly_sai=monthly_sai,
                                      monthly_height_top=monthly_htop,
                                      monthly_height_bot=monthly_hbot,
                                      months=(6, 7),
                                      noveg=0, maxveg=maxveg)

        # Patch 1 (noveg): should be zeros
        @test sp.mlai2t[1, 1] == 0.0
        @test sp.mlai2t[1, 2] == 0.0
        @test sp.msai2t[1, 1] == 0.0
        @test sp.mhvt2t[1, 1] == 0.0
        @test sp.mhvb2t[1, 1] == 0.0

        # Patch 2 (PFT=1, gridcell=1): should pick up months 6 and 7
        @test sp.mlai2t[2, 1] ≈ 3.5
        @test sp.mlai2t[2, 2] ≈ 4.2
        @test sp.msai2t[2, 1] ≈ 0.5
        @test sp.msai2t[2, 2] ≈ 0.6
        @test sp.mhvt2t[2, 1] ≈ 15.0
        @test sp.mhvt2t[2, 2] ≈ 16.0
        @test sp.mhvb2t[2, 1] ≈ 5.0
        @test sp.mhvb2t[2, 2] ≈ 5.5

        # Patch 3 (PFT=12, gridcell=2)
        @test sp.mlai2t[3, 1] ≈ 1.0
        @test sp.mlai2t[3, 2] ≈ 1.2

        # mlaidiff = mlai2t[:,1] - mlai2t[:,2]
        @test cs.mlaidiff_patch[1] ≈ 0.0
        @test cs.mlaidiff_patch[2] ≈ 3.5 - 4.2
        @test cs.mlaidiff_patch[3] ≈ 1.0 - 1.2
    end

    # -----------------------------------------------------------------------
    # Test read_annual_vegetation!
    # -----------------------------------------------------------------------
    @testset "read_annual_vegetation!" begin
        np = 2
        ng = 1
        maxveg = 14

        cs = CLM.CanopyStateData()
        CLM.canopystate_init!(cs, np)

        patch = CLM.PatchData()
        CLM.patch_init!(patch, np)
        patch.itype    .= [0, 5]        # noveg, broadleaf_evergreen_temperate_tree
        patch.gridcell .= [1, 1]

        # Create mock monthly LAI data
        monthly_lai = zeros(ng, maxveg + 1, 12)
        for k in 1:12
            monthly_lai[1, 6, k] = Float64(k) * 0.5  # PFT 5 (index 6)
        end

        CLM.read_annual_vegetation!(cs, patch, 1:np;
                                     monthly_lai=monthly_lai,
                                     noveg=0, maxveg=maxveg)

        # Patch 1 (noveg): all months zero
        for k in 1:12
            @test cs.annlai_patch[1, k] == 0.0
        end

        # Patch 2 (PFT=5): should match input
        for k in 1:12
            @test cs.annlai_patch[2, k] ≈ Float64(k) * 0.5
        end
    end

    # -----------------------------------------------------------------------
    # Test NDAYPM constant
    # -----------------------------------------------------------------------
    @testset "NDAYPM constant" begin
        @test length(CLM.NDAYPM) == 12
        @test sum(CLM.NDAYPM) == 365
        @test CLM.NDAYPM[2] == 28
    end

    # -----------------------------------------------------------------------
    # Integration: monthly interpolation → phenology computation
    # -----------------------------------------------------------------------
    @testset "Integration: interp_monthly_veg then satellite_phenology" begin
        np = 2
        nc = 1

        sp = CLM.SatellitePhenologyData()
        CLM.satellite_phenology_init!(sp, np)

        # Set some monthly data directly (as if read_monthly_vegetation! was called)
        sp.mlai2t[1, :] .= [0.0, 0.0]   # noveg
        sp.mlai2t[2, :] .= [4.0, 5.0]   # vegetated

        sp.msai2t[1, :] .= [0.0, 0.0]
        sp.msai2t[2, :] .= [0.8, 1.0]

        sp.mhvt2t[1, :] .= [0.0, 0.0]
        sp.mhvt2t[2, :] .= [15.0, 16.0]

        sp.mhvb2t[1, :] .= [0.0, 0.0]
        sp.mhvb2t[2, :] .= [5.0, 5.5]

        # Compute time weights for July 15
        CLM.interp_monthly_veg!(sp; kmo=7, kda=15)

        cs = CLM.CanopyStateData()
        CLM.canopystate_init!(cs, np)

        wd = CLM.WaterDiagnosticBulkData()
        wd.frac_sno_col = [0.0]
        wd.snow_depth_col = [0.0]

        patch = CLM.PatchData()
        CLM.patch_init!(patch, np)
        patch.itype  .= [0, 5]
        patch.column .= [1, 1]

        CLM.satellite_phenology!(sp, cs, wd, patch, trues(np), 1:np)

        # With no snow, elai should equal tlai (if >= 0.05)
        expected_tlai = sp.timwt[1] * 4.0 + sp.timwt[2] * 5.0
        @test cs.tlai_patch[2] ≈ expected_tlai
        @test cs.elai_patch[2] ≈ expected_tlai
        @test cs.frac_veg_nosno_alb_patch[2] == 1
    end

end
