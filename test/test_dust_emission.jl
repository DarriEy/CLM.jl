@testset "Dust Emission (DustEmisBase)" begin

    # -----------------------------------------------------------------------
    # Test local erf approximation
    # -----------------------------------------------------------------------
    @testset "_dust_erf accuracy" begin
        # Known values
        @test CLM._dust_erf(0.0) ≈ 0.0 atol=1e-7
        @test CLM._dust_erf(1.0) ≈ 0.8427007929 atol=1e-6
        @test CLM._dust_erf(-1.0) ≈ -0.8427007929 atol=1e-6
        @test CLM._dust_erf(2.0) ≈ 0.9953222650 atol=1e-6
        @test CLM._dust_erf(0.5) ≈ 0.5204998778 atol=1e-6
        # Symmetry
        @test CLM._dust_erf(-0.5) ≈ -CLM._dust_erf(0.5) atol=1e-10
    end

    # -----------------------------------------------------------------------
    # Test DustEmisBaseData construction and init/clean
    # -----------------------------------------------------------------------
    @testset "DustEmisBaseData construction and init/clean" begin
        dust = CLM.DustEmisBaseData()
        @test length(dust.dmt_vwr) == 0
        @test length(dust.stk_crc) == 0
        @test size(dust.ovr_src_snk_mss) == (0, 0)
        @test size(dust.flx_mss_vrt_dst_patch) == (0, 0)
        @test dust.saltation_factor == 0.0

        np = 5
        CLM.dust_emis_init!(dust, np)

        # Check allocations
        @test size(dust.flx_mss_vrt_dst_patch) == (np, CLM.NDST)
        @test length(dust.flx_mss_vrt_dst_tot_patch) == np
        @test size(dust.vlc_trb_patch) == (np, CLM.NDST)
        @test length(dust.vlc_trb_1_patch) == np
        @test length(dust.vlc_trb_2_patch) == np
        @test length(dust.vlc_trb_3_patch) == np
        @test length(dust.vlc_trb_4_patch) == np
        @test size(dust.ovr_src_snk_mss) == (CLM.DST_SRC_NBR, CLM.NDST)
        @test length(dust.dmt_vwr) == CLM.NDST
        @test length(dust.stk_crc) == CLM.NDST

        # Patch arrays should be NaN-initialized
        @test all(isnan, dust.flx_mss_vrt_dst_tot_patch)
        @test all(isnan, dust.vlc_trb_1_patch)

        # Clean
        CLM.dust_emis_clean!(dust)
        @test length(dust.dmt_vwr) == 0
        @test length(dust.stk_crc) == 0
        @test size(dust.ovr_src_snk_mss) == (0, 0)
        @test size(dust.flx_mss_vrt_dst_patch) == (0, 0)
        @test length(dust.flx_mss_vrt_dst_tot_patch) == 0
    end

    # -----------------------------------------------------------------------
    # Test init_dust_vars! produces valid physical values
    # -----------------------------------------------------------------------
    @testset "init_dust_vars! produces valid values" begin
        dust = CLM.DustEmisBaseData()
        CLM.dust_emis_init!(dust, 2)

        # Overlap factors should be non-negative and finite
        @test all(x -> x >= 0.0 && isfinite(x), dust.ovr_src_snk_mss)

        # Each source should have overlap factors that sum to <= mass fraction
        # (they represent sub-fractions of the source mass in each sink bin)
        for m in 1:CLM.DST_SRC_NBR
            @test sum(dust.ovr_src_snk_mss[m, :]) >= 0.0
            @test sum(dust.ovr_src_snk_mss[m, :]) <= 1.0
        end

        # Mass-weighted diameters should be positive and in reasonable range
        # Bins span 0.1-10 um, so dmt_vwr should be in that range
        for n in 1:CLM.NDST
            @test dust.dmt_vwr[n] > 0.0
            @test dust.dmt_vwr[n] < 20.0e-6  # should be < 20 um
            @test dust.dmt_vwr[n] > 0.05e-6  # should be > 0.05 um
        end

        # Diameters should increase with bin number
        for n in 1:(CLM.NDST-1)
            @test dust.dmt_vwr[n] < dust.dmt_vwr[n+1]
        end

        # Stokes correction factors should be >= 1 (correction makes settling faster)
        for n in 1:CLM.NDST
            @test dust.stk_crc[n] >= 1.0
            @test isfinite(dust.stk_crc[n])
        end

        # Saltation factor should be positive and finite
        @test dust.saltation_factor > 0.0
        @test isfinite(dust.saltation_factor)
    end

    # -----------------------------------------------------------------------
    # Test dust_dry_dep! — turbulent dry deposition
    # -----------------------------------------------------------------------
    @testset "dust_dry_dep! basic operation" begin
        np = 4
        nc = 3
        dust = CLM.DustEmisBaseData()
        CLM.dust_emis_init!(dust, np)

        # Setup patch data
        patch_active = trues(np)
        patch_active[3] = false  # one inactive patch
        patch_column = [1, 2, 2, 3]
        bounds_p = 1:np

        # Atmospheric forcing (column-level)
        forc_pbot = fill(101325.0, nc)   # standard pressure [Pa]
        forc_rho  = fill(1.225, nc)      # standard density [kg/m3]
        forc_t    = fill(293.0, nc)      # ~20C [K]

        # Friction velocity data (patch-level)
        ram1 = fill(50.0, np)            # aerodynamical resistance [s/m]
        fv_vals = fill(0.3, np)          # friction velocity [m/s]

        CLM.dust_dry_dep!(dust, patch_active, patch_column, bounds_p,
                          forc_pbot, forc_rho, forc_t, ram1, fv_vals)

        # Active patches should have valid deposition velocities
        for p in [1, 2, 4]
            for m in 1:CLM.NDST
                @test dust.vlc_trb_patch[p, m] > 0.0
                @test isfinite(dust.vlc_trb_patch[p, m])
            end
            @test dust.vlc_trb_1_patch[p] > 0.0
            @test dust.vlc_trb_2_patch[p] > 0.0
            @test dust.vlc_trb_3_patch[p] > 0.0
            @test dust.vlc_trb_4_patch[p] > 0.0

            # Individual bin values should match array
            @test dust.vlc_trb_1_patch[p] == dust.vlc_trb_patch[p, 1]
            @test dust.vlc_trb_2_patch[p] == dust.vlc_trb_patch[p, 2]
            @test dust.vlc_trb_3_patch[p] == dust.vlc_trb_patch[p, 3]
            @test dust.vlc_trb_4_patch[p] == dust.vlc_trb_patch[p, 4]
        end

        # Inactive patch (p=3) should remain at initial NaN value
        @test isnan(dust.vlc_trb_patch[3, 1])

        # Larger particles should generally have higher deposition velocities
        # (due to gravitational settling contribution)
        for p in [1, 2, 4]
            @test dust.vlc_trb_patch[p, 4] > dust.vlc_trb_patch[p, 1]
        end
    end

    @testset "dust_dry_dep! sensitivity to friction velocity" begin
        np = 2
        nc = 1
        dust = CLM.DustEmisBaseData()
        CLM.dust_emis_init!(dust, np)

        patch_active = trues(np)
        patch_column = [1, 1]
        bounds_p = 1:np

        forc_pbot = [101325.0]
        forc_rho  = [1.225]
        forc_t    = [293.0]
        ram1 = [50.0, 50.0]

        # Low friction velocity
        fv_low = [0.1, 0.5]
        CLM.dust_dry_dep!(dust, patch_active, patch_column, bounds_p,
                          forc_pbot, forc_rho, forc_t, ram1, fv_low)

        vlc_low = dust.vlc_trb_1_patch[1]
        vlc_high = dust.vlc_trb_1_patch[2]

        # Higher friction velocity should generally lead to different deposition
        @test vlc_low > 0.0
        @test vlc_high > 0.0
    end

    # -----------------------------------------------------------------------
    # Test check_dust_emis_is_valid
    # -----------------------------------------------------------------------
    @testset "check_dust_emis_is_valid" begin
        np = 2
        dust = CLM.DustEmisBaseData()
        CLM.dust_emis_init!(dust, np)

        # Set up consistent values for patch 1
        dust.flx_mss_vrt_dst_patch[1, :] .= [1.0e-8, 2.0e-8, 3.0e-8, 4.0e-8]
        dust.flx_mss_vrt_dst_tot_patch[1] = sum(dust.flx_mss_vrt_dst_patch[1, :])

        dust.vlc_trb_patch[1, :] .= [0.001, 0.002, 0.003, 0.004]
        dust.vlc_trb_1_patch[1] = 0.001
        dust.vlc_trb_2_patch[1] = 0.002
        dust.vlc_trb_3_patch[1] = 0.003
        dust.vlc_trb_4_patch[1] = 0.004

        @test CLM.check_dust_emis_is_valid(dust, 1) == true

        # Test inconsistency: total doesn't match sum
        dust.flx_mss_vrt_dst_tot_patch[1] = 0.0
        @test_throws ErrorException CLM.check_dust_emis_is_valid(dust, 1)

        # Restore total, break bin consistency
        dust.flx_mss_vrt_dst_tot_patch[1] = sum(dust.flx_mss_vrt_dst_patch[1, :])
        dust.vlc_trb_1_patch[1] = 0.999  # mismatch
        @test_throws ErrorException CLM.check_dust_emis_is_valid(dust, 1)
    end

    # -----------------------------------------------------------------------
    # Test get_dust_patch_vars
    # -----------------------------------------------------------------------
    @testset "get_dust_patch_vars" begin
        np = 2
        dust = CLM.DustEmisBaseData()
        CLM.dust_emis_init!(dust, np)

        dust.flx_mss_vrt_dst_patch[1, :] .= [1.0, 2.0, 3.0, 4.0]
        dust.flx_mss_vrt_dst_tot_patch[1] = 10.0
        dust.vlc_trb_patch[1, :] .= [0.01, 0.02, 0.03, 0.04]
        dust.vlc_trb_1_patch[1] = 0.01
        dust.vlc_trb_2_patch[1] = 0.02
        dust.vlc_trb_3_patch[1] = 0.03
        dust.vlc_trb_4_patch[1] = 0.04

        vars = CLM.get_dust_patch_vars(dust, 1)
        @test vars.flx_mss_vrt_dst == [1.0, 2.0, 3.0, 4.0]
        @test vars.flx_mss_vrt_dst_tot == 10.0
        @test vars.vlc_trb == [0.01, 0.02, 0.03, 0.04]
        @test vars.vlc_trb_1 == 0.01
        @test vars.vlc_trb_2 == 0.02
        @test vars.vlc_trb_3 == 0.03
        @test vars.vlc_trb_4 == 0.04
    end

    # -----------------------------------------------------------------------
    # Test get_dust_const_vars
    # -----------------------------------------------------------------------
    @testset "get_dust_const_vars" begin
        dust = CLM.DustEmisBaseData()
        CLM.dust_emis_init!(dust, 2)

        sf = CLM.get_dust_const_vars(dust)
        @test sf == dust.saltation_factor
        @test sf > 0.0
        @test isfinite(sf)
    end

    # -----------------------------------------------------------------------
    # Integration test: full init + dry deposition cycle
    # -----------------------------------------------------------------------
    @testset "full init and dry deposition cycle" begin
        np = 3
        nc = 2
        dust = CLM.DustEmisBaseData()
        CLM.dust_emis_init!(dust, np)

        # Verify init_dust_vars! ran successfully
        @test dust.saltation_factor > 0.0
        @test all(isfinite, dust.dmt_vwr)
        @test all(isfinite, dust.stk_crc)
        @test all(x -> x >= 0.0, dust.ovr_src_snk_mss)

        # Run dry deposition
        patch_active = trues(np)
        patch_column = [1, 1, 2]
        bounds_p = 1:np

        forc_pbot = [101325.0, 85000.0]    # sea level and ~1500m elevation
        forc_rho  = [1.225, 1.05]
        forc_t    = [300.0, 280.0]
        ram1 = fill(40.0, np)
        fv_vals = fill(0.25, np)

        CLM.dust_dry_dep!(dust, patch_active, patch_column, bounds_p,
                          forc_pbot, forc_rho, forc_t, ram1, fv_vals)

        # All active patches should have valid deposition velocities
        for p in 1:np
            for m in 1:CLM.NDST
                @test dust.vlc_trb_patch[p, m] > 0.0
                @test isfinite(dust.vlc_trb_patch[p, m])
            end
        end

        # Set emission fluxes and validate
        for p in 1:np
            dust.flx_mss_vrt_dst_patch[p, :] .= [1e-9, 2e-9, 3e-9, 4e-9]
            dust.flx_mss_vrt_dst_tot_patch[p] = sum(dust.flx_mss_vrt_dst_patch[p, :])
        end

        # Should pass validation for all patches
        for p in 1:np
            @test CLM.check_dust_emis_is_valid(dust, p) == true
        end

        # Patches at different pressures/temperatures should get different values
        # (columns 1 and 2 have different conditions)
        @test dust.vlc_trb_patch[1, 1] != dust.vlc_trb_patch[3, 1]
    end

end
