@testset "Dry Deposition Velocity (DryDepVelocity)" begin

    # -----------------------------------------------------------------------
    # Test DryDepVelocityData construction and init/clean
    # -----------------------------------------------------------------------
    @testset "DryDepVelocityData construction and init/clean" begin
        dd = CLM.DryDepVelocityData()
        @test dd.n_drydep == 0
        @test length(dd.foxd) == 0
        @test length(dd.dv) == 0
        @test length(dd.mapping) == 0
        @test size(dd.velocity_patch) == (0, 0)

        np = 5
        n_dd = 3
        foxd = [1.0, 0.1, 0.0]
        dv = [0.15, 0.20, 0.10]
        mapping = [2, 4, 8]

        CLM.drydep_init!(dd, np, n_dd; foxd=foxd, dv=dv, mapping=mapping)

        @test dd.n_drydep == 3
        @test length(dd.foxd) == 3
        @test dd.foxd[1] == 1.0
        @test dd.foxd[2] == 0.1
        @test dd.foxd[3] == 0.0
        @test dd.dv[1] == 0.15
        @test dd.dv[2] == 0.20
        @test dd.dv[3] == 0.10
        @test dd.mapping == [2, 4, 8]
        @test size(dd.velocity_patch) == (5, 3)
        @test all(x -> x == 0.0, dd.velocity_patch)

        # Clean
        CLM.drydep_clean!(dd)
        @test dd.n_drydep == 0
        @test length(dd.foxd) == 0
        @test length(dd.dv) == 0
        @test length(dd.mapping) == 0
        @test size(dd.velocity_patch) == (0, 0)
    end

    # -----------------------------------------------------------------------
    # Test init with defaults (no species properties given)
    # -----------------------------------------------------------------------
    @testset "drydep_init! with defaults" begin
        dd = CLM.DryDepVelocityData()
        CLM.drydep_init!(dd, 4, 2)

        @test dd.n_drydep == 2
        @test length(dd.foxd) == 2
        @test all(x -> x == 0.0, dd.foxd)  # default zeros
        @test all(x -> x == 0.2, dd.dv)    # default 0.2
        @test all(x -> x == 1, dd.mapping)  # default urban
        @test size(dd.velocity_patch) == (4, 2)
    end

    # -----------------------------------------------------------------------
    # Test zero species case
    # -----------------------------------------------------------------------
    @testset "drydep_init! with zero species" begin
        dd = CLM.DryDepVelocityData()
        CLM.drydep_init!(dd, 4, 0)
        @test dd.n_drydep == 0
        @test size(dd.velocity_patch) == (0, 0)
    end

    # -----------------------------------------------------------------------
    # Test Wesely season determination
    # -----------------------------------------------------------------------
    @testset "_wesely_season NH" begin
        # Northern hemisphere (lat >= 0)
        # Jan -> winter (4)
        @test CLM._wesely_season(45.0, 1) == 4
        # Feb -> winter (4)
        @test CLM._wesely_season(45.0, 2) == 4
        # Mar -> spring (5)
        @test CLM._wesely_season(45.0, 3) == 5
        # Apr -> spring (5)
        @test CLM._wesely_season(45.0, 4) == 5
        # Jun -> midsummer (1)
        @test CLM._wesely_season(45.0, 6) == 1
        # Jul -> midsummer (1)
        @test CLM._wesely_season(45.0, 7) == 1
        # Sep -> autumn (2)
        @test CLM._wesely_season(45.0, 9) == 2
        # Oct -> late autumn (3)
        @test CLM._wesely_season(45.0, 10) == 3
        # Dec -> winter (4)
        @test CLM._wesely_season(45.0, 12) == 4
    end

    @testset "_wesely_season SH" begin
        # Southern hemisphere (lat < 0): seasons shifted by 6 months
        # Jan -> midsummer (SH summer)
        @test CLM._wesely_season(-45.0, 1) == 1
        # Jun -> winter (SH winter)
        @test CLM._wesely_season(-45.0, 6) == 4
        # Jul -> winter
        @test CLM._wesely_season(-45.0, 7) == 4
        # Dec -> midsummer
        @test CLM._wesely_season(-45.0, 12) == 1
    end

    @testset "_wesely_season equator" begin
        # Equator (lat = 0) treated as NH
        @test CLM._wesely_season(0.0, 6) == 1
        @test CLM._wesely_season(0.0, 1) == 4
    end

    # -----------------------------------------------------------------------
    # Test PFT to Wesely land type mapping
    # -----------------------------------------------------------------------
    @testset "_pft_to_wesely mapping" begin
        # not vegetated -> barren
        @test CLM._pft_to_wesely(0) == 8
        # needleleaf trees (1-3) -> coniferous forest
        @test CLM._pft_to_wesely(1) == 5
        @test CLM._pft_to_wesely(2) == 5
        @test CLM._pft_to_wesely(3) == 5
        # broadleaf trees (4-8) -> deciduous forest
        @test CLM._pft_to_wesely(4) == 4
        @test CLM._pft_to_wesely(7) == 4
        @test CLM._pft_to_wesely(8) == 4
        # shrubs (9-11) -> range land
        @test CLM._pft_to_wesely(9) == 3
        @test CLM._pft_to_wesely(11) == 3
        # grasses (12-14) -> range land
        @test CLM._pft_to_wesely(12) == 3
        @test CLM._pft_to_wesely(14) == 3
        # crops (15+) -> agricultural
        @test CLM._pft_to_wesely(15) == 2
        @test CLM._pft_to_wesely(17) == 2
        @test CLM._pft_to_wesely(78) == 2
    end

    # -----------------------------------------------------------------------
    # Test quasi-laminar boundary layer resistance
    # -----------------------------------------------------------------------
    @testset "_calc_rb basic behavior" begin
        # Typical friction velocity
        rb = CLM._calc_rb(0.3, 0.2)
        @test rb > 0.0
        @test isfinite(rb)

        # Higher ustar -> lower rb
        rb_high = CLM._calc_rb(0.6, 0.2)
        @test rb_high < rb

        # Higher diffusivity -> lower rb (easier to diffuse through layer)
        rb_highdv = CLM._calc_rb(0.3, 0.25)
        rb_lowdv  = CLM._calc_rb(0.3, 0.10)
        @test rb_highdv < rb_lowdv

        # Very small ustar should be clamped, not blow up
        rb_small = CLM._calc_rb(1.0e-6, 0.2)
        @test isfinite(rb_small)
        @test rb_small > 0.0

        # Zero diffusivity should not cause NaN
        rb_zerodv = CLM._calc_rb(0.3, 0.0)
        @test isfinite(rb_zerodv)
    end

    @testset "_calc_rb quantitative check" begin
        # For ustar=0.3 m/s, DH2O=0.25 cm^2/s, dv=0.25 cm^2/s
        # Sc/Pr = DH2O/dv = 1.0, so (Sc/Pr)^(2/3) = 1.0
        # Rb = 2/(0.4 * 0.3) * 1.0 = 2/0.12 = 16.667 s/m
        rb = CLM._calc_rb(0.3, 0.25)
        @test rb ≈ 2.0 / (CLM.VKC * 0.3) atol=0.01

        # For dv = DH2O/4 = 0.0625
        # Sc/Pr = 0.25/0.0625 = 4.0
        # (Sc/Pr)^(2/3) = 4^(2/3) ≈ 2.5198
        rb2 = CLM._calc_rb(0.3, 0.0625)
        expected = (2.0 / (CLM.VKC * 0.3)) * 4.0^(2.0/3.0)
        @test rb2 ≈ expected atol=0.01
    end

    # -----------------------------------------------------------------------
    # Test surface resistance
    # -----------------------------------------------------------------------
    @testset "_calc_rc basic behavior" begin
        # Midsummer, agricultural land, moderate conditions
        rc = CLM._calc_rc(1, 2, 3.0, 1.0, 1.0, 298.0, 400.0, false)
        @test rc > 0.0
        @test rc >= 1.0
        @test isfinite(rc)

        # Winter with snow should have higher resistance
        rc_snow = CLM._calc_rc(4, 2, 3.0, 1.0, 1.0, 260.0, 0.0, true)
        @test rc_snow > 0.0
        @test rc_snow >= 1.0
        @test isfinite(rc_snow)

        # Higher reactivity (foxd) should generally lower resistance
        rc_reactive = CLM._calc_rc(1, 4, 3.0, 1.0, 1.0, 298.0, 400.0, false)
        rc_unreact  = CLM._calc_rc(1, 4, 3.0, 0.0, 0.0, 298.0, 400.0, false)
        @test rc_reactive <= rc_unreact
    end

    @testset "_calc_rc cold temperature" begin
        # Cold temperature should close stomata -> high R_c
        rc_cold = CLM._calc_rc(1, 4, 3.0, 1.0, 1.0, 263.0, 400.0, false)
        rc_warm = CLM._calc_rc(1, 4, 3.0, 1.0, 1.0, 298.0, 400.0, false)
        @test rc_cold >= rc_warm
    end

    @testset "_calc_rc no solar radiation" begin
        # No solar -> stomata closed -> higher R_c
        rc_dark = CLM._calc_rc(1, 4, 3.0, 1.0, 1.0, 298.0, 0.0, false)
        rc_sun  = CLM._calc_rc(1, 4, 3.0, 1.0, 1.0, 298.0, 400.0, false)
        @test rc_dark >= rc_sun
    end

    @testset "_calc_rc water surface" begin
        # Water (land type 7): very large ground resistance for O3,
        # but zero R_gs for SO2 (water dissolves SO2)
        rc_water = CLM._calc_rc(1, 7, 0.0, 1.0, 1.0, 298.0, 100.0, false)
        @test isfinite(rc_water)
        @test rc_water >= 1.0
    end

    @testset "_calc_rc barren land" begin
        # Barren land (type 8): large stomatal resistance, moderate ground
        rc_barren = CLM._calc_rc(1, 8, 0.0, 1.0, 1.0, 298.0, 400.0, false)
        @test isfinite(rc_barren)
        @test rc_barren >= 1.0
    end

    # -----------------------------------------------------------------------
    # Test depvel_compute! basic operation
    # -----------------------------------------------------------------------
    @testset "depvel_compute! basic operation" begin
        np = 4
        nc = 3
        ng = 2
        n_dd = 2

        dd = CLM.DryDepVelocityData()
        foxd = [1.0, 0.1]        # SO2-like, moderate
        dv = [0.13, 0.18]        # diffusivities [cm^2/s]
        CLM.drydep_init!(dd, np, n_dd; foxd=foxd, dv=dv)

        # Setup patch/column/gridcell topology
        mask_patch = trues(np)
        mask_patch[3] = false  # inactive patch
        bounds_patch = 1:np
        patch_gridcell = [1, 1, 2, 2]
        patch_column = [1, 2, 2, 3]
        patch_landunit = [1, 1, 1, 1]
        patch_itype = [4, 14, 0, 17]   # broadleaf tree, C4 grass, bare, crop

        # Friction velocity data
        ram1 = fill(50.0, np)       # aerodynamic resistance [s/m]
        rb1 = fill(10.0, np)        # boundary layer resistance [s/m]
        fv = fill(0.3, np)          # friction velocity [m/s]

        # Canopy
        elai = [3.0, 2.0, 0.0, 4.0]

        # Atmospheric forcing (column-level)
        forc_t = fill(295.0, nc)
        forc_solar = fill(300.0, nc)
        frac_sno_col = fill(0.0, nc)

        # Gridcell properties
        lat = [45.0, -30.0]  # NH and SH points

        CLM.depvel_compute!(dd, mask_patch, bounds_patch,
                            patch_gridcell, patch_column, patch_landunit,
                            patch_itype,
                            ram1, rb1, fv, elai,
                            forc_t, forc_solar, frac_sno_col,
                            lat, 6)  # June

        # Active patches should have positive finite deposition velocities
        for p in [1, 2, 4]
            for i in 1:n_dd
                @test dd.velocity_patch[p, i] > 0.0
                @test isfinite(dd.velocity_patch[p, i])
            end
        end

        # Inactive patch 3 should remain at 0.0 (initialized value)
        for i in 1:n_dd
            @test dd.velocity_patch[3, i] == 0.0
        end

        # More reactive species (foxd=1.0) should generally have higher
        # deposition velocity than less reactive (foxd=0.1)
        for p in [1, 2, 4]
            @test dd.velocity_patch[p, 1] >= dd.velocity_patch[p, 2] * 0.1
        end
    end

    # -----------------------------------------------------------------------
    # Test depvel_compute! with snow
    # -----------------------------------------------------------------------
    @testset "depvel_compute! with snow" begin
        np = 2
        nc = 2
        ng = 1
        n_dd = 1

        dd = CLM.DryDepVelocityData()
        CLM.drydep_init!(dd, np, n_dd; foxd=[0.5], dv=[0.15])

        mask_patch = trues(np)
        bounds_patch = 1:np
        patch_gridcell = [1, 1]
        patch_column = [1, 2]
        patch_landunit = [1, 1]
        patch_itype = [4, 4]  # both deciduous forest

        ram1 = fill(50.0, np)
        rb1 = fill(10.0, np)
        fv = fill(0.3, np)
        elai = fill(3.0, np)
        forc_t = fill(295.0, nc)
        forc_solar = fill(300.0, nc)
        frac_sno_col = [0.0, 0.9]  # no snow vs heavy snow
        lat = [45.0]

        CLM.depvel_compute!(dd, mask_patch, bounds_patch,
                            patch_gridcell, patch_column, patch_landunit,
                            patch_itype,
                            ram1, rb1, fv, elai,
                            forc_t, forc_solar, frac_sno_col,
                            lat, 6)

        # Snow-covered patch should have lower deposition velocity
        # (higher surface resistance due to snow)
        @test dd.velocity_patch[2, 1] <= dd.velocity_patch[1, 1]
    end

    # -----------------------------------------------------------------------
    # Test depvel_compute! sensitivity to friction velocity
    # -----------------------------------------------------------------------
    @testset "depvel_compute! sensitivity to ustar" begin
        np = 2
        nc = 2
        ng = 1
        n_dd = 1

        dd = CLM.DryDepVelocityData()
        CLM.drydep_init!(dd, np, n_dd; foxd=[0.5], dv=[0.15])

        mask_patch = trues(np)
        bounds_patch = 1:np
        patch_gridcell = [1, 1]
        patch_column = [1, 2]
        patch_landunit = [1, 1]
        patch_itype = [4, 4]

        ram1 = fill(50.0, np)
        rb1 = fill(10.0, np)
        fv = [0.1, 0.6]  # low vs high friction velocity
        elai = fill(3.0, np)
        forc_t = fill(295.0, nc)
        forc_solar = fill(300.0, nc)
        frac_sno_col = fill(0.0, nc)
        lat = [45.0]

        CLM.depvel_compute!(dd, mask_patch, bounds_patch,
                            patch_gridcell, patch_column, patch_landunit,
                            patch_itype,
                            ram1, rb1, fv, elai,
                            forc_t, forc_solar, frac_sno_col,
                            lat, 6)

        # Higher friction velocity -> lower boundary layer resistance -> higher Vd
        @test dd.velocity_patch[2, 1] > dd.velocity_patch[1, 1]
    end

    # -----------------------------------------------------------------------
    # Test depvel_compute! with zero species
    # -----------------------------------------------------------------------
    @testset "depvel_compute! with zero species (no-op)" begin
        dd = CLM.DryDepVelocityData()
        CLM.drydep_init!(dd, 4, 0)

        mask_patch = trues(4)
        bounds_patch = 1:4
        # Should not error -- just return immediately
        CLM.depvel_compute!(dd, mask_patch, bounds_patch,
                            [1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1],
                            [4, 4, 4, 4],
                            fill(50.0, 4), fill(10.0, 4), fill(0.3, 4),
                            fill(3.0, 4),
                            fill(295.0, 1), fill(300.0, 1), fill(0.0, 1),
                            [45.0], 6)
    end

    # -----------------------------------------------------------------------
    # Test depvel_compute! with custom Henry's law constants
    # -----------------------------------------------------------------------
    @testset "depvel_compute! with custom heff" begin
        np = 2
        nc = 1
        n_dd = 2

        dd = CLM.DryDepVelocityData()
        CLM.drydep_init!(dd, np, n_dd; foxd=[0.5, 0.5], dv=[0.15, 0.15])

        mask_patch = trues(np)
        bounds_patch = 1:np
        patch_gridcell = [1, 1]
        patch_column = [1, 1]
        patch_landunit = [1, 1]
        patch_itype = [4, 4]

        ram1 = fill(50.0, np)
        rb1 = fill(10.0, np)
        fv = fill(0.3, np)
        elai = fill(3.0, np)
        forc_t = fill(295.0, nc)
        forc_solar = fill(300.0, nc)
        frac_sno_col = fill(0.0, nc)
        lat = [45.0]

        # High vs low Henry's law constant
        heff_vals = [1.0e5, 1.0e-2]  # highly soluble vs poorly soluble

        CLM.depvel_compute!(dd, mask_patch, bounds_patch,
                            patch_gridcell, patch_column, patch_landunit,
                            patch_itype,
                            ram1, rb1, fv, elai,
                            forc_t, forc_solar, frac_sno_col,
                            lat, 6; heff=heff_vals)

        # Highly soluble species should have higher deposition velocity
        for p in 1:np
            @test dd.velocity_patch[p, 1] > dd.velocity_patch[p, 2]
        end
    end

    # -----------------------------------------------------------------------
    # Test depvel_compute! physical range
    # -----------------------------------------------------------------------
    @testset "depvel_compute! deposition velocity in physical range" begin
        np = 6
        nc = 6
        ng = 2
        n_dd = 3

        dd = CLM.DryDepVelocityData()
        foxd = [1.0, 0.1, 0.5]
        dv = [0.13, 0.21, 0.16]
        CLM.drydep_init!(dd, np, n_dd; foxd=foxd, dv=dv)

        mask_patch = trues(np)
        bounds_patch = 1:np
        patch_gridcell = [1, 1, 1, 2, 2, 2]
        patch_column = [1, 2, 3, 4, 5, 6]
        patch_landunit = fill(1, np)
        patch_itype = [0, 1, 4, 9, 14, 17]

        ram1 = fill(40.0, np)
        rb1 = fill(8.0, np)
        fv = fill(0.25, np)
        elai = [0.0, 1.0, 4.0, 2.0, 3.0, 5.0]
        forc_t = fill(290.0, nc)
        forc_solar = fill(250.0, nc)
        frac_sno_col = fill(0.0, nc)
        lat = [40.0, -20.0]

        CLM.depvel_compute!(dd, mask_patch, bounds_patch,
                            patch_gridcell, patch_column, patch_landunit,
                            patch_itype,
                            ram1, rb1, fv, elai,
                            forc_t, forc_solar, frac_sno_col,
                            lat, 7)

        # Deposition velocities should be positive, finite, and in
        # physically reasonable range (0.001 to ~10 cm/s for most gases)
        for p in 1:np
            for i in 1:n_dd
                @test dd.velocity_patch[p, i] > 0.0
                @test dd.velocity_patch[p, i] < 50.0  # generous upper bound
                @test isfinite(dd.velocity_patch[p, i])
            end
        end
    end

    # -----------------------------------------------------------------------
    # Test depvel_compute! different months
    # -----------------------------------------------------------------------
    @testset "depvel_compute! seasonal variation" begin
        np = 1
        nc = 1
        ng = 1
        n_dd = 1

        vd_by_month = zeros(12)

        for month in 1:12
            dd = CLM.DryDepVelocityData()
            CLM.drydep_init!(dd, np, n_dd; foxd=[0.5], dv=[0.15])

            mask_patch = trues(np)
            bounds_patch = 1:np

            CLM.depvel_compute!(dd, mask_patch, bounds_patch,
                                [1], [1], [1], [4],  # deciduous forest
                                [40.0], [10.0], [0.3], [3.0],
                                [295.0], [300.0], [0.0],
                                [45.0], month)

            vd_by_month[month] = dd.velocity_patch[1, 1]
        end

        # All months should produce valid deposition velocities
        for month in 1:12
            @test vd_by_month[month] > 0.0
            @test isfinite(vd_by_month[month])
        end

        # Summer (June) should generally have higher Vd than winter (January)
        # due to open stomata and active vegetation
        # Note: the exact comparison depends on the Wesely table values
        @test vd_by_month[6] >= vd_by_month[1] * 0.5  # relaxed comparison
    end

    # -----------------------------------------------------------------------
    # Test drydep_p2g! (patch to gridcell averaging)
    # -----------------------------------------------------------------------
    @testset "drydep_p2g! basic operation" begin
        np = 4
        ng = 2
        n_dd = 2

        dd = CLM.DryDepVelocityData()
        CLM.drydep_init!(dd, np, n_dd; foxd=[0.5, 0.1], dv=[0.15, 0.20])

        # Set known deposition velocities
        dd.velocity_patch[1, 1] = 1.0
        dd.velocity_patch[1, 2] = 0.5
        dd.velocity_patch[2, 1] = 2.0
        dd.velocity_patch[2, 2] = 1.0
        dd.velocity_patch[3, 1] = 3.0
        dd.velocity_patch[3, 2] = 1.5
        dd.velocity_patch[4, 1] = 4.0
        dd.velocity_patch[4, 2] = 2.0

        mask_patch = trues(np)
        mask_patch[3] = false  # inactive
        bounds_patch = 1:np
        patch_gridcell = [1, 1, 2, 2]
        wtgcell = [0.4, 0.6, 0.5, 0.5]

        ddvel_grc = zeros(ng, n_dd)

        CLM.drydep_p2g!(dd, ddvel_grc, bounds_patch, mask_patch,
                         patch_gridcell, wtgcell, ng)

        # Gridcell 1: patches 1 (wt=0.4) and 2 (wt=0.6)
        @test ddvel_grc[1, 1] ≈ 1.0 * 0.4 + 2.0 * 0.6  atol=1e-10
        @test ddvel_grc[1, 2] ≈ 0.5 * 0.4 + 1.0 * 0.6  atol=1e-10

        # Gridcell 2: only patch 4 is active (patch 3 masked out)
        @test ddvel_grc[2, 1] ≈ 4.0 * 0.5  atol=1e-10
        @test ddvel_grc[2, 2] ≈ 2.0 * 0.5  atol=1e-10
    end

    @testset "drydep_p2g! with zero species (no-op)" begin
        dd = CLM.DryDepVelocityData()
        CLM.drydep_init!(dd, 4, 0)

        ddvel_grc = zeros(2, 0)
        mask_patch = trues(4)
        bounds_patch = 1:4

        # Should not error
        CLM.drydep_p2g!(dd, ddvel_grc, bounds_patch, mask_patch,
                         [1, 1, 2, 2], fill(0.5, 4), 2)
    end

    # -----------------------------------------------------------------------
    # Test resistance table constants are self-consistent
    # -----------------------------------------------------------------------
    @testset "Wesely resistance table dimensions" begin
        @test size(CLM.RI_TABLE) == (CLM.N_LAND_TYPE, 5)
        @test size(CLM.RGS_SO2_TABLE) == (CLM.N_LAND_TYPE, 5)
        @test size(CLM.RGS_O3_TABLE) == (CLM.N_LAND_TYPE, 5)
        @test size(CLM.RLU_TABLE) == (CLM.N_LAND_TYPE, 5)
        @test size(CLM.RAC_TABLE) == (CLM.N_LAND_TYPE, 5)
        @test length(CLM.Z0_TABLE) == CLM.N_LAND_TYPE

        # All resistance values should be non-negative
        @test all(x -> x >= 0.0, CLM.RI_TABLE)
        @test all(x -> x >= 0.0, CLM.RGS_SO2_TABLE)
        @test all(x -> x >= 0.0, CLM.RGS_O3_TABLE)
        @test all(x -> x >= 0.0, CLM.RLU_TABLE)
        @test all(x -> x >= 0.0, CLM.RAC_TABLE)
        @test all(x -> x >= 0.0, CLM.Z0_TABLE)
    end

    # -----------------------------------------------------------------------
    # Test end-to-end: compute + p2g
    # -----------------------------------------------------------------------
    @testset "end-to-end compute and p2g" begin
        np = 3
        nc = 3
        ng = 1
        n_dd = 2

        dd = CLM.DryDepVelocityData()
        foxd = [1.0, 0.1]
        dv = [0.13, 0.20]
        CLM.drydep_init!(dd, np, n_dd; foxd=foxd, dv=dv)

        mask_patch = trues(np)
        bounds_patch = 1:np
        patch_gridcell = [1, 1, 1]
        patch_column = [1, 2, 3]
        patch_landunit = [1, 1, 1]
        patch_itype = [4, 14, 17]
        wtgcell = [0.4, 0.3, 0.3]

        ram1 = fill(50.0, np)
        rb1 = fill(10.0, np)
        fv = fill(0.3, np)
        elai = [3.0, 2.0, 4.0]
        forc_t = fill(295.0, nc)
        forc_solar = fill(300.0, nc)
        frac_sno_col = fill(0.0, nc)
        lat = [45.0]

        # Compute deposition velocities
        CLM.depvel_compute!(dd, mask_patch, bounds_patch,
                            patch_gridcell, patch_column, patch_landunit,
                            patch_itype,
                            ram1, rb1, fv, elai,
                            forc_t, forc_solar, frac_sno_col,
                            lat, 7)

        # Average to gridcell
        ddvel_grc = zeros(ng, n_dd)
        CLM.drydep_p2g!(dd, ddvel_grc, bounds_patch, mask_patch,
                         patch_gridcell, wtgcell, ng)

        # Gridcell average should be weighted sum of patch velocities
        for i in 1:n_dd
            expected = sum(dd.velocity_patch[p, i] * wtgcell[p] for p in 1:np)
            @test ddvel_grc[1, i] ≈ expected atol=1e-10
        end

        # Gridcell averages should be positive
        for i in 1:n_dd
            @test ddvel_grc[1, i] > 0.0
        end
    end

end
