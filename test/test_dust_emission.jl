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

    # -----------------------------------------------------------------------
    # Threshold soil moisture / clay mass-fraction functions
    # (ported from SoilStateInitTimeConstMod.F90)
    # -----------------------------------------------------------------------
    @testset "threshold soil moisture (Zender2003 / Kok2014)" begin
        # Zender2003: 0.17 + 0.14*clay*0.01, clay in percent.
        @test CLM.threshold_soil_moist_zender2003(0.0)   ≈ 0.17 atol=1e-12
        @test CLM.threshold_soil_moist_zender2003(100.0) ≈ 0.17 + 0.14 atol=1e-12
        @test CLM.threshold_soil_moist_zender2003(50.0)  ≈ 0.17 + 0.14 * 0.5 atol=1e-12
        # Linear in clay.
        d = CLM.threshold_soil_moist_zender2003(40.0) -
            CLM.threshold_soil_moist_zender2003(20.0)
        @test d ≈ 0.14 * 0.20 atol=1e-12
        # Out-of-bounds clay errors.
        @test_throws ErrorException CLM.threshold_soil_moist_zender2003(-1.0)
        @test_throws ErrorException CLM.threshold_soil_moist_zender2003(101.0)

        # Kok2014: 0.01*(0.17*clay + 0.0014*clay^2), clay in percent.
        @test CLM.threshold_soil_moist_kok2014(0.0) ≈ 0.0 atol=1e-12
        for clay in (10.0, 25.0, 50.0, 80.0)
            @test CLM.threshold_soil_moist_kok2014(clay) ≈
                0.01 * (0.17 * clay + 0.0014 * clay * clay) atol=1e-12
        end
        # Kok2014 is parabolic (convex) — second difference positive.
        f = CLM.threshold_soil_moist_kok2014
        @test (f(60.0) - 2*f(40.0) + f(20.0)) > 0.0

        # MassFracClay: min(clay*0.01, 0.20).
        @test CLM.mass_frac_clay(10.0) ≈ 0.10 atol=1e-12
        @test CLM.mass_frac_clay(20.0) ≈ 0.20 atol=1e-12
        @test CLM.mass_frac_clay(50.0) ≈ 0.20 atol=1e-12   # clamped
        @test CLM.mass_frac_clay(0.0)  ≈ 0.0  atol=1e-12
    end

    # -----------------------------------------------------------------------
    # Wind-driven dust mobilization kernels: Zender2003 (default) + Leung2023.
    # One soil landunit (itype=ISTSOIL), one column, one patch. We toggle wind,
    # soil moisture, snow and VAI to exercise the threshold/gating logic.
    # -----------------------------------------------------------------------
    @testset "dust mobilization (Zender2003 + Leung2023)" begin
        ndst = CLM.NDST

        # Build a 1-landunit / 1-column / 1-patch setup with given conditions.
        # Returns (dust, kwargs...) ready to pass to the emission routines.
        function make_scene(; tlai=0.0, tsai=0.0, fv=0.6, u10=12.0,
                            frac_sno=0.0, h2osoi_vol=0.05, watsat=0.45,
                            gwc_thr=0.17, mss_frc_cly_vld=0.15, forc_rho=1.225,
                            obu=-50.0, itype=CLM.noveg, lun_itype=CLM.ISTSOIL,
                            dpfct_rock=1.0)
            dust = CLM.DustEmisBaseData()
            CLM.dust_emis_init!(dust, 1; nc=1)
            dust.dpfct_rock_patch[1] = dpfct_rock
            # liquid surface water (kg/m2): all liquid when warm.
            return (dust=dust,
                nolakep_mask=trues(1), patch_active=trues(1),
                patch_column=[1], patch_landunit=[1], patch_itype=[itype],
                patch_wtlunit=[1.0], lun_itype=[lun_itype], bounds_p=1:1, nl=1,
                forc_rho=[forc_rho], gwc_thr=[gwc_thr],
                mss_frc_cly_vld=[mss_frc_cly_vld],
                watsat=reshape([watsat], 1, 1),
                tlai=[tlai], tsai=[tsai], frac_sno=[frac_sno],
                h2osoi_vol=reshape([h2osoi_vol], 1, 1),
                h2osoi_liq=reshape([30.0], 1, 1),
                h2osoi_ice=reshape([0.0], 1, 1),
                fv=[fv], u10=[u10], obu=[obu])
        end

        run_zender(s) = CLM.dust_emission_zender2003!(
            s.dust, s.nolakep_mask, s.patch_active, s.patch_column,
            s.patch_landunit, s.patch_wtlunit, s.lun_itype, s.bounds_p, s.nl,
            s.forc_rho, s.gwc_thr, s.mss_frc_cly_vld, s.watsat, s.tlai, s.tsai,
            s.frac_sno, s.h2osoi_vol, s.h2osoi_liq, s.h2osoi_ice, s.fv, s.u10)

        run_leung(s) = CLM.dust_emission_leung2023!(
            s.dust, s.nolakep_mask, s.patch_active, s.patch_column,
            s.patch_landunit, s.patch_itype, s.patch_wtlunit, s.lun_itype,
            s.bounds_p, s.nl, s.forc_rho, s.gwc_thr, s.mss_frc_cly_vld, s.watsat,
            s.tlai, s.tsai, s.frac_sno, s.h2osoi_vol, s.h2osoi_liq, s.h2osoi_ice,
            s.fv, s.obu)

        @testset "Zender2003: dry erodible high wind → positive, size-distributed" begin
            s = make_scene(tlai=0.0, tsai=0.0, fv=0.8, u10=15.0,
                           h2osoi_vol=0.02, frac_sno=0.0)
            run_zender(s)
            d = s.dust
            # finite, non-negative, and total = sum over bins
            @test all(isfinite, d.flx_mss_vrt_dst_patch[1, :])
            @test all(>=(0.0), d.flx_mss_vrt_dst_patch[1, :])
            @test d.flx_mss_vrt_dst_tot_patch[1] > 0.0
            @test d.flx_mss_vrt_dst_tot_patch[1] ≈
                  sum(d.flx_mss_vrt_dst_patch[1, :]) atol=0.0
            # at least two bins receive dust (size-distributed)
            @test count(>(0.0), d.flx_mss_vrt_dst_patch[1, :]) >= 2
        end

        @testset "Zender2003: wet soil suppresses emission vs dry" begin
            # Moisture above threshold raises the saltation threshold friction
            # velocity, so a wet column emits strictly less than a dry one.
            s_dry = make_scene(h2osoi_vol=0.02, watsat=0.45, gwc_thr=0.10,
                               fv=0.8, u10=15.0)
            s_wet = make_scene(h2osoi_vol=0.40, watsat=0.45, gwc_thr=0.10,
                               fv=0.8, u10=15.0)
            run_zender(s_dry); run_zender(s_wet)
            @test s_dry.dust.flx_mss_vrt_dst_tot_patch[1] > 0.0
            @test s_wet.dust.flx_mss_vrt_dst_tot_patch[1] <
                  s_dry.dust.flx_mss_vrt_dst_tot_patch[1]
            @test s_wet.dust.flx_mss_vrt_dst_tot_patch[1] >= 0.0
            # Saturated soil (gwc_sfc ≫ threshold, very weak wind) → exactly zero.
            s_sat = make_scene(h2osoi_vol=0.44, watsat=0.45, gwc_thr=0.10,
                               fv=0.2, u10=2.0)
            run_zender(s_sat)
            @test s_sat.dust.flx_mss_vrt_dst_tot_patch[1] == 0.0
        end

        @testset "Zender2003: full snow cover → zero emission" begin
            s = make_scene(frac_sno=1.0, fv=0.8, u10=15.0, h2osoi_vol=0.02)
            run_zender(s)
            @test s.dust.flx_mss_vrt_dst_tot_patch[1] == 0.0
        end

        @testset "Zender2003: below-threshold wind → zero emission" begin
            s = make_scene(fv=0.05, u10=0.5, h2osoi_vol=0.02)
            run_zender(s)
            @test s.dust.flx_mss_vrt_dst_tot_patch[1] == 0.0
        end

        @testset "Zender2003: VAI above threshold → zero (no bare ground)" begin
            s = make_scene(tlai=0.5, tsai=0.5, fv=0.8, u10=15.0)  # VAI=1.0 > 0.3
            run_zender(s)
            @test s.dust.flx_mss_vrt_dst_tot_patch[1] == 0.0
        end

        @testset "Zender2003: emission increases with wind" begin
            s1 = make_scene(fv=0.5, u10=10.0, h2osoi_vol=0.02)
            s2 = make_scene(fv=0.9, u10=18.0, h2osoi_vol=0.02)
            run_zender(s1); run_zender(s2)
            @test s2.dust.flx_mss_vrt_dst_tot_patch[1] >
                  s1.dust.flx_mss_vrt_dst_tot_patch[1] > 0.0
        end

        @testset "Leung2023: dry erodible high wind → positive, size-distributed" begin
            s = make_scene(tlai=0.0, tsai=0.0, fv=0.9, u10=15.0,
                           h2osoi_vol=0.02, frac_sno=0.0, dpfct_rock=1.0)
            run_leung(s)
            d = s.dust
            @test all(isfinite, d.flx_mss_vrt_dst_patch[1, :])
            @test all(>=(0.0), d.flx_mss_vrt_dst_patch[1, :])
            @test d.flx_mss_vrt_dst_tot_patch[1] > 0.0
            @test d.flx_mss_vrt_dst_tot_patch[1] ≈
                  sum(d.flx_mss_vrt_dst_patch[1, :]) atol=0.0
            @test count(>(0.0), d.flx_mss_vrt_dst_patch[1, :]) >= 2
        end

        @testset "Leung2023: snow / low wind → zero emission" begin
            for s in (make_scene(frac_sno=1.0, fv=0.9, h2osoi_vol=0.02),
                      make_scene(fv=0.03, h2osoi_vol=0.02))
                run_leung(s)
                @test s.dust.flx_mss_vrt_dst_tot_patch[1] == 0.0
            end
        end

        @testset "Leung2023: wet soil suppresses emission vs dry" begin
            s_dry = make_scene(h2osoi_vol=0.02, watsat=0.45, gwc_thr=0.10,
                               fv=0.9, u10=15.0)
            s_wet = make_scene(h2osoi_vol=0.40, watsat=0.45, gwc_thr=0.10,
                               fv=0.9, u10=15.0)
            run_leung(s_dry); run_leung(s_wet)
            @test s_dry.dust.flx_mss_vrt_dst_tot_patch[1] > 0.0
            @test s_wet.dust.flx_mss_vrt_dst_tot_patch[1] <
                  s_dry.dust.flx_mss_vrt_dst_tot_patch[1]
            @test s_wet.dust.flx_mss_vrt_dst_tot_patch[1] >= 0.0
        end

        @testset "Leung2023 vs Zender2003 differ (different threshold/scheme)" begin
            # Same dry, high-wind, erodible scene → both emit, but the schemes
            # use different thresholds + emission equations, so totals differ.
            sz = make_scene(fv=0.85, u10=15.0, h2osoi_vol=0.02, dpfct_rock=1.0)
            sl = make_scene(fv=0.85, u10=15.0, h2osoi_vol=0.02, dpfct_rock=1.0)
            run_zender(sz); run_leung(sl)
            @test sz.dust.flx_mss_vrt_dst_tot_patch[1] > 0.0
            @test sl.dust.flx_mss_vrt_dst_tot_patch[1] > 0.0
            @test sz.dust.flx_mss_vrt_dst_tot_patch[1] !=
                  sl.dust.flx_mss_vrt_dst_tot_patch[1]
        end

        @testset "factory dispatch + default-off no-op" begin
            # Inactive config: no-op, fluxes stay at their NaN-init value.
            dust = CLM.DustEmisBaseData()
            CLM.dust_emis_init!(dust, 1; nc=1)
            pre = copy(dust.flx_mss_vrt_dst_patch)
            cfg_off = CLM.DustEmisConfig()
            @test cfg_off.active == false
            # The factory returns early; nothing should run (no inputs needed).
            CLM.dust_emission!(cfg_off, dust, trues(1), nothing, nothing, 1:1, 1,
                               [1.225], nothing, nothing, nothing, nothing, nothing)
            @test isequal(dust.flx_mss_vrt_dst_patch, pre)  # untouched (NaN==NaN)

            # Unknown method errors.
            cfg_bad = CLM.DustEmisConfig(active=true, method="nope")
            @test_throws ErrorException CLM.dust_emission!(
                cfg_bad, dust, trues(1), nothing, nothing, 1:1, 1, [1.225],
                nothing, nothing, nothing, nothing, nothing)
        end
    end

end
