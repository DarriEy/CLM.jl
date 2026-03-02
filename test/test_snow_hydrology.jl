@testset "Snow Hydrology Module" begin
    # ------------------------------------------------------------------
    # Tests for SnowHydrologyMod port.
    # Verifies:
    #   1. SnowHydrologyParams defaults
    #   2. Module-level constants
    #   3. new_snow_bulk_density!
    #   4. update_state_add_new_snow!
    #   5. build_filter_thawed_wetland_thin_snowpack!
    #   6. update_state_remove_snow_thawed_wetlands!
    #   7. bulk_remove_snow_thawed_wetlands!
    #   8. build_filter_snowpack_initialized!
    #   9. update_state_initialize_snow_pack!
    #  10. bulk_initialize_snow_pack!
    #  11. update_state_top_layer_fluxes!
    #  12. bulk_flux_snow_percolation!
    #  13. update_state_snow_percolation!
    #  14. post_percolation_adjust_layer_thicknesses!
    #  15. overburden_compaction_anderson1976
    #  16. overburden_compaction_vionnet2012
    #  17. combo_scalar
    #  18. mass_weighted_snow_radius
    #  19. build_snow_filter!
    #  20. init_flux_snow_capping!
    #  21. snow_hydrology_set_control_for_testing!
    #  22. bulkdiag_snow_water_accumulated_snow!
    #  23. sum_flux_add_snow_percolation!
    # ------------------------------------------------------------------

    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
    nlevsno = CLM.varpar.nlevsno

    # Initialize snow layer thickness arrays for tests
    CLM.snow_hydrology_set_control_for_testing!()
    # Manually init dzmin/dzmax arrays
    dzmin = zeros(Float64, nlevsno)
    dzmax_l = zeros(Float64, nlevsno)
    dzmax_u = zeros(Float64, nlevsno)
    dzmin[1] = 0.010; dzmax_l[1] = 0.03; dzmax_u[1] = 0.02
    dzmin[2] = 0.015; dzmax_l[2] = 0.07; dzmax_u[2] = 0.05
    for j in 3:nlevsno
        dzmin[j] = dzmax_u[j-1] * 0.5
        dzmax_u[j] = 2.0 * dzmax_u[j-1] + 0.01
        dzmax_l[j] = dzmax_u[j] + dzmax_l[j-1]
        if j == nlevsno
            dzmax_u[j] = floatmax(Float64)
            dzmax_l[j] = floatmax(Float64)
        end
    end
    CLM.SNOW_DZMIN[] = dzmin
    CLM.SNOW_DZMAX_L[] = dzmax_l
    CLM.SNOW_DZMAX_U[] = dzmax_u

    # ------------------------------------------------------------------
    # 1. SnowHydrologyParams defaults
    # ------------------------------------------------------------------
    @testset "SnowHydrologyParams defaults" begin
        p = CLM.SnowHydrologyParams()
        @test p.wimp == 0.05
        @test p.ssi == 0.033
        @test p.snw_rds_min == 54.526
        @test p.upplim_destruct_metamorph == 100.0
    end

    # ------------------------------------------------------------------
    # 2. Module constants
    # ------------------------------------------------------------------
    @testset "Module constants" begin
        @test CLM.SCVNG_FCT_MLT_OCPHI == 0.20
        @test CLM.SCVNG_FCT_MLT_OCPHO == 0.03
        @test CLM.LO_TMP_DNS_SLATER2017 == 2
        @test CLM.LO_TMP_DNS_TRUNCATED_ANDERSON1976 == 1
        @test CLM.RESET_SNOW_H2OSNO == 35.0
        @test CLM.LSADZ == 0.03
    end

    # ------------------------------------------------------------------
    # 3. new_snow_bulk_density!
    # ------------------------------------------------------------------
    @testset "new_snow_bulk_density!" begin
        nc = 4
        bifall = zeros(Float64, nc)
        # Warm: T > tfrz + 2
        # Moderate: tfrz-15 < T < tfrz+2
        # Cold truncated: T < tfrz-15 with TruncatedAnderson1976
        # Cold Slater: T < tfrz-15 with Slater2017
        forc_t = [CLM.TFRZ + 5.0, CLM.TFRZ, CLM.TFRZ - 20.0, CLM.TFRZ - 30.0]
        forc_wind = [3.0]
        col_gridcell = [1, 1, 1, 1]
        mask = trues(nc)
        bounds = 1:nc

        # Test with Slater2017 (default)
        CLM.NEW_SNOW_DENSITY[] = CLM.LO_TMP_DNS_SLATER2017
        CLM.WIND_DEPENDENT_SNOW_DENSITY[] = false
        CLM.new_snow_bulk_density!(bifall, forc_t, forc_wind, col_gridcell, mask, bounds)

        # Warm case: 50 + 1.7 * 17^1.5
        @test bifall[1] ≈ 50.0 + 1.7 * 17.0^1.5 atol=0.01
        # Moderate case: 50 + 1.7 * (0 + 15)^1.5
        @test bifall[2] ≈ 50.0 + 1.7 * 15.0^1.5 atol=0.01
        # Cold Slater case: should be positive
        @test bifall[3] > 0.0
        @test bifall[4] > 0.0

        # Test TruncatedAnderson1976
        CLM.NEW_SNOW_DENSITY[] = CLM.LO_TMP_DNS_TRUNCATED_ANDERSON1976
        bifall_trunc = zeros(Float64, nc)
        CLM.new_snow_bulk_density!(bifall_trunc, forc_t, forc_wind, col_gridcell, mask, bounds)
        @test bifall_trunc[3] == 50.0  # Truncated at 50 for cold temps
    end

    # ------------------------------------------------------------------
    # 4. update_state_add_new_snow!
    # ------------------------------------------------------------------
    @testset "update_state_add_new_snow!" begin
        nc = 2
        h2osno_no_layers = [10.0, 0.0]
        h2osoi_ice = zeros(Float64, nc, nlevsno + 15)  # snow + soil layers
        snl = [0, -1]
        qflx_snow_grnd = [0.5, 0.5]
        mask = trues(nc)
        dtime = 1800.0

        # Column 1: no layers, snow goes to h2osno_no_layers
        # Column 2: 1 layer, snow goes to top ice layer
        jj_top_c2 = (snl[2] + 1) + nlevsno  # layer 0 for snl=-1
        h2osoi_ice[2, jj_top_c2] = 5.0

        CLM.update_state_add_new_snow!(h2osno_no_layers, h2osoi_ice,
            dtime, snl, qflx_snow_grnd, mask, 1:nc, nlevsno)

        @test h2osno_no_layers[1] ≈ 10.0 + 0.5 * 1800.0
        @test h2osoi_ice[2, jj_top_c2] ≈ 5.0 + 0.5 * 1800.0
    end

    # ------------------------------------------------------------------
    # 5. build_filter_thawed_wetland_thin_snowpack!
    # ------------------------------------------------------------------
    @testset "build_filter_thawed_wetland_thin_snowpack!" begin
        nc = 3
        mask_out = falses(nc)
        t_grnd = [CLM.TFRZ + 1.0, CLM.TFRZ - 1.0, CLM.TFRZ + 1.0]
        lun_itype_col = [CLM.ISTWET, CLM.ISTWET, CLM.ISTSOIL]
        snl = [0, 0, 0]
        mask_nolake = trues(nc)

        CLM.build_filter_thawed_wetland_thin_snowpack!(mask_out,
            t_grnd, lun_itype_col, snl, mask_nolake, 1:nc)

        @test mask_out[1] == true   # wetland, thawed, no layers
        @test mask_out[2] == false  # wetland, frozen
        @test mask_out[3] == false  # not wetland
    end

    # ------------------------------------------------------------------
    # 6-7. Remove snow from thawed wetlands
    # ------------------------------------------------------------------
    @testset "remove snow from thawed wetlands" begin
        nc = 2
        mask = BitVector([true, false])

        h2osno_no_layers = [5.0, 5.0]
        CLM.update_state_remove_snow_thawed_wetlands!(h2osno_no_layers, mask, 1:nc)
        @test h2osno_no_layers[1] == 0.0
        @test h2osno_no_layers[2] == 5.0

        snow_depth = [0.1, 0.1]
        CLM.bulk_remove_snow_thawed_wetlands!(snow_depth, mask, 1:nc)
        @test snow_depth[1] == 0.0
        @test snow_depth[2] == 0.1
    end

    # ------------------------------------------------------------------
    # 8. build_filter_snowpack_initialized!
    # ------------------------------------------------------------------
    @testset "build_filter_snowpack_initialized!" begin
        nc = 3
        mask_out = falses(nc)
        snl = [0, 0, -1]
        lun_itype_col = [CLM.ISTSOIL, CLM.ISTSOIL, CLM.ISTSOIL]
        frac_sno_eff = [1.0, 1.0, 1.0]
        snow_depth = [0.02, 0.005, 0.1]  # above dzmin, below dzmin, already has layers
        qflx_snow_grnd = [1.0, 1.0, 0.0]
        mask = trues(nc)

        CLM.build_filter_snowpack_initialized!(mask_out,
            snl, lun_itype_col, frac_sno_eff, snow_depth, qflx_snow_grnd, mask, 1:nc)

        @test mask_out[1] == true   # no layers, depth >= dzmin
        @test mask_out[2] == false  # depth < dzmin
        @test mask_out[3] == false  # already has layers
    end

    # ------------------------------------------------------------------
    # 9. update_state_initialize_snow_pack!
    # ------------------------------------------------------------------
    @testset "update_state_initialize_snow_pack!" begin
        nc = 1
        h2osno_no_layers = [20.0]
        h2osoi_ice = zeros(Float64, nc, nlevsno + 15)
        h2osoi_liq = zeros(Float64, nc, nlevsno + 15)
        mask = trues(nc)

        CLM.update_state_initialize_snow_pack!(h2osno_no_layers,
            h2osoi_ice, h2osoi_liq, mask, 1:nc, nlevsno)

        jj_zero = 0 + nlevsno
        @test h2osoi_ice[1, jj_zero] == 20.0
        @test h2osoi_liq[1, jj_zero] == 0.0
        @test h2osno_no_layers[1] == 0.0
    end

    # ------------------------------------------------------------------
    # 10. bulk_initialize_snow_pack!
    # ------------------------------------------------------------------
    @testset "bulk_initialize_snow_pack!" begin
        nc = 1
        ntot = nlevsno + 15
        snl = [0]
        zi = zeros(Float64, nc, ntot)
        dz = zeros(Float64, nc, ntot)
        z = zeros(Float64, nc, ntot)
        t_soisno = fill(CLM.TFRZ, nc, ntot)
        frac_iceold = zeros(Float64, nc, ntot)
        snomelt_accum = [0.5]
        forc_t = [CLM.TFRZ - 5.0]
        snow_depth = [0.05]
        mask = trues(nc)

        CLM.bulk_initialize_snow_pack!(snl, zi, dz, z, t_soisno, frac_iceold,
            snomelt_accum, forc_t, snow_depth, mask, 1:nc, nlevsno)

        jj_zero = 0 + nlevsno
        jj_m1 = -1 + nlevsno
        @test snl[1] == -1
        @test dz[1, jj_zero] == 0.05
        @test z[1, jj_zero] ≈ -0.025
        @test zi[1, jj_m1] ≈ -0.05
        @test t_soisno[1, jj_zero] ≈ CLM.TFRZ - 5.0
        @test frac_iceold[1, jj_zero] == 1.0
        @test snomelt_accum[1] == 0.0
    end

    # ------------------------------------------------------------------
    # 11. update_state_top_layer_fluxes!
    # ------------------------------------------------------------------
    @testset "update_state_top_layer_fluxes!" begin
        nc = 1
        ntot = nlevsno + 15
        h2osoi_ice = zeros(Float64, nc, ntot)
        h2osoi_liq = zeros(Float64, nc, ntot)
        snl = [-1]
        jj_top = (snl[1] + 1) + nlevsno  # layer 0
        h2osoi_ice[1, jj_top] = 10.0
        h2osoi_liq[1, jj_top] = 2.0
        frac_sno_eff = [1.0]
        dtime = 1800.0
        mask_snow = trues(nc)

        CLM.update_state_top_layer_fluxes!(h2osoi_ice, h2osoi_liq,
            dtime, snl, frac_sno_eff,
            [0.001],  # qflx_soliddew_to_top_layer
            [0.0],    # qflx_solidevap_from_top_layer
            [0.0],    # qflx_liq_grnd
            [0.001],  # qflx_liqdew_to_top_layer
            [0.0],    # qflx_liqevap_from_top_layer
            mask_snow, 1:nc, nlevsno)

        @test h2osoi_ice[1, jj_top] ≈ 10.0 + 1.0 * 0.001 * 1800.0
        @test h2osoi_liq[1, jj_top] ≈ 2.0 + 1.0 * 0.001 * 1800.0
    end

    # ------------------------------------------------------------------
    # 12. bulk_flux_snow_percolation!
    # ------------------------------------------------------------------
    @testset "bulk_flux_snow_percolation!" begin
        nc = 1
        ntot = nlevsno + 15
        qflx_snow_percolation = zeros(Float64, nc, nlevsno)
        snl = [-2]  # two snow layers
        dz = zeros(Float64, nc, ntot)
        h2osoi_ice = zeros(Float64, nc, ntot)
        h2osoi_liq = zeros(Float64, nc, ntot)
        frac_sno_eff = [1.0]
        dtime = 1800.0
        mask_snow = trues(nc)

        # Set up two layers (j=-1 and j=0)
        jj_m1 = (-1) + nlevsno
        jj_0  = 0 + nlevsno
        dz[1, jj_m1] = 0.03
        dz[1, jj_0]  = 0.05
        h2osoi_ice[1, jj_m1] = 5.0
        h2osoi_ice[1, jj_0]  = 8.0
        h2osoi_liq[1, jj_m1] = 3.0
        h2osoi_liq[1, jj_0]  = 1.0

        CLM.bulk_flux_snow_percolation!(qflx_snow_percolation,
            dtime, snl, dz, frac_sno_eff, h2osoi_ice, h2osoi_liq,
            mask_snow, 1:nc, nlevsno)

        # Percolation should be non-negative
        @test qflx_snow_percolation[1, jj_m1] >= 0.0
        @test qflx_snow_percolation[1, jj_0] >= 0.0
    end

    # ------------------------------------------------------------------
    # 13. update_state_snow_percolation!
    # ------------------------------------------------------------------
    @testset "update_state_snow_percolation!" begin
        nc = 1
        ntot = nlevsno + 15
        h2osoi_liq = zeros(Float64, nc, ntot)
        snl = [-2]
        dtime = 1800.0
        mask_snow = trues(nc)

        jj_m1 = (-1) + nlevsno
        jj_0  = 0 + nlevsno
        h2osoi_liq[1, jj_m1] = 5.0
        h2osoi_liq[1, jj_0]  = 3.0

        # Set percolation flux from top layer
        qflx_snow_percolation = zeros(Float64, nc, nlevsno)
        qflx_snow_percolation[1, jj_m1] = 0.001  # mm/s out of layer -1
        qflx_snow_percolation[1, jj_0]  = 0.0005 # mm/s out of layer 0

        CLM.update_state_snow_percolation!(h2osoi_liq,
            dtime, snl, qflx_snow_percolation, mask_snow, 1:nc, nlevsno)

        # Top layer: loses percolation
        @test h2osoi_liq[1, jj_m1] ≈ 5.0 - 0.001 * 1800.0
        # Bottom layer: gains from above, loses own percolation
        @test h2osoi_liq[1, jj_0] ≈ 3.0 + 0.001 * 1800.0 - 0.0005 * 1800.0
    end

    # ------------------------------------------------------------------
    # 14. post_percolation_adjust_layer_thicknesses!
    # ------------------------------------------------------------------
    @testset "post_percolation_adjust_layer_thicknesses!" begin
        nc = 1
        ntot = nlevsno + 15
        dz = zeros(Float64, nc, ntot)
        snl = [-1]
        h2osoi_ice = zeros(Float64, nc, ntot)
        h2osoi_liq = zeros(Float64, nc, ntot)
        mask_snow = trues(nc)

        jj_0 = 0 + nlevsno
        dz[1, jj_0] = 0.001  # Very thin layer
        h2osoi_ice[1, jj_0] = 50.0  # 50 kg/m2 ice
        h2osoi_liq[1, jj_0] = 10.0  # 10 kg/m2 liquid

        CLM.post_percolation_adjust_layer_thicknesses!(dz,
            snl, h2osoi_ice, h2osoi_liq, mask_snow, 1:nc, nlevsno)

        expected_min = 10.0 / CLM.DENH2O + 50.0 / CLM.DENICE
        @test dz[1, jj_0] ≈ expected_min
    end

    # ------------------------------------------------------------------
    # 15. overburden_compaction_anderson1976
    # ------------------------------------------------------------------
    @testset "overburden_compaction_anderson1976" begin
        result = CLM.overburden_compaction_anderson1976(50.0, 10.0, 5.0, 200.0)
        @test result < 0.0  # Should be negative (compaction)
    end

    # ------------------------------------------------------------------
    # 16. overburden_compaction_vionnet2012
    # ------------------------------------------------------------------
    @testset "overburden_compaction_vionnet2012" begin
        result = CLM.overburden_compaction_vionnet2012(1.0, 0.05, 50.0, 10.0, 5.0, 200.0)
        @test result < 0.0  # Should be negative (compaction)
    end

    # ------------------------------------------------------------------
    # 17. combo_scalar
    # ------------------------------------------------------------------
    @testset "combo_scalar" begin
        t1 = CLM.TFRZ - 5.0
        t2 = CLM.TFRZ - 10.0
        (dzc, wliqc, wicec, tc) = CLM.combo_scalar(0.03, 1.0, 5.0, t1,
                                                     0.05, 2.0, 8.0, t2)
        @test dzc ≈ 0.08
        @test wliqc ≈ 3.0
        @test wicec ≈ 13.0
        @test tc < CLM.TFRZ  # Combined temp should be below freezing
    end

    # ------------------------------------------------------------------
    # 18. mass_weighted_snow_radius
    # ------------------------------------------------------------------
    @testset "mass_weighted_snow_radius" begin
        r = CLM.mass_weighted_snow_radius(100.0, 200.0, 5.0, 5.0)
        @test r ≈ 150.0

        # Test clamping to min
        r_min = CLM.mass_weighted_snow_radius(10.0, 20.0, 5.0, 5.0)
        @test r_min == CLM.snowhydrology_params.snw_rds_min

        # Test clamping to max
        r_max = CLM.mass_weighted_snow_radius(2000.0, 2000.0, 5.0, 5.0)
        @test r_max == CLM.SNW_RDS_MAX
    end

    # ------------------------------------------------------------------
    # 19. build_snow_filter!
    # ------------------------------------------------------------------
    @testset "build_snow_filter!" begin
        nc = 4
        snl = [-2, 0, -1, 0]
        mask_nolake = trues(nc)
        mask_snow = falses(nc)
        mask_nosnow = falses(nc)

        CLM.build_snow_filter!(mask_snow, mask_nosnow, snl, mask_nolake, 1:nc)

        @test mask_snow == BitVector([true, false, true, false])
        @test mask_nosnow == BitVector([false, true, false, true])
    end

    # ------------------------------------------------------------------
    # 20. init_flux_snow_capping!
    # ------------------------------------------------------------------
    @testset "init_flux_snow_capping!" begin
        nc = 2
        qflx_snwcp_ice = [1.0, 2.0]
        qflx_snwcp_liq = [1.0, 2.0]
        qflx_snwcp_discarded_ice = [1.0, 2.0]
        qflx_snwcp_discarded_liq = [1.0, 2.0]
        mask = trues(nc)

        CLM.init_flux_snow_capping!(qflx_snwcp_ice, qflx_snwcp_liq,
            qflx_snwcp_discarded_ice, qflx_snwcp_discarded_liq, mask, 1:nc)

        @test all(qflx_snwcp_ice .== 0.0)
        @test all(qflx_snwcp_liq .== 0.0)
        @test all(qflx_snwcp_discarded_ice .== 0.0)
        @test all(qflx_snwcp_discarded_liq .== 0.0)
    end

    # ------------------------------------------------------------------
    # 21. snow_hydrology_set_control_for_testing!
    # ------------------------------------------------------------------
    @testset "snow_hydrology_set_control_for_testing!" begin
        CLM.snow_hydrology_set_control_for_testing!(
            wind_dep_snow_density = true,
            new_snow_density_method = CLM.LO_TMP_DNS_TRUNCATED_ANDERSON1976
        )
        @test CLM.WIND_DEPENDENT_SNOW_DENSITY[] == true
        @test CLM.NEW_SNOW_DENSITY[] == CLM.LO_TMP_DNS_TRUNCATED_ANDERSON1976
        @test CLM.SNOW_DZMIN_1[] == 0.010
        @test CLM.SNOW_DZMAX_U_1[] == 0.02

        # Reset
        CLM.snow_hydrology_set_control_for_testing!(
            wind_dep_snow_density = false,
            new_snow_density_method = CLM.LO_TMP_DNS_SLATER2017
        )
    end

    # ------------------------------------------------------------------
    # 22. bulkdiag_snow_water_accumulated_snow!
    # ------------------------------------------------------------------
    @testset "bulkdiag_snow_water_accumulated_snow!" begin
        nc = 2
        int_snow = [100.0, 50.0]
        frac_sno = [0.8, 0.5]
        snow_depth = [0.1, 0.05]
        dtime = 1800.0
        frac_sno_eff = [0.8, 0.5]
        qflx_soliddew = [0.001, 0.0]
        qflx_liqdew = [0.0, 0.0]
        qflx_liq_grnd = [0.0, 0.0]
        h2osno_no_layers = [0.0, 0.0]  # nosnow column with no snow

        mask_snow = BitVector([true, false])
        mask_nosnow = BitVector([false, true])

        CLM.bulkdiag_snow_water_accumulated_snow!(int_snow, frac_sno, snow_depth,
            dtime, frac_sno_eff, qflx_soliddew, qflx_liqdew, qflx_liq_grnd,
            h2osno_no_layers, mask_snow, mask_nosnow, 1:nc)

        # Snow column: int_snow increases
        @test int_snow[1] ≈ 100.0 + 0.8 * 0.001 * 1800.0
        # No-snow column with no h2osno: everything resets to 0
        @test int_snow[2] == 0.0
        @test frac_sno[2] == 0.0
        @test snow_depth[2] == 0.0
    end

    # ------------------------------------------------------------------
    # 23. sum_flux_add_snow_percolation!
    # ------------------------------------------------------------------
    @testset "sum_flux_add_snow_percolation!" begin
        nc = 2
        qflx_snow_drain = [0.0, 0.0]
        qflx_rain_plus_snomelt = [0.0, 0.0]
        frac_sno_eff = [0.8, 0.5]
        qflx_snow_percolation_bottom = [0.001, 0.0]
        qflx_liq_grnd = [0.002, 0.003]
        qflx_snomelt = [0.0, 0.001]

        mask_snow = BitVector([true, false])
        mask_nosnow = BitVector([false, true])

        CLM.sum_flux_add_snow_percolation!(qflx_snow_drain, qflx_rain_plus_snomelt,
            frac_sno_eff, qflx_snow_percolation_bottom, qflx_liq_grnd, qflx_snomelt,
            mask_snow, mask_nosnow, 1:nc)

        # Snow column
        @test qflx_snow_drain[1] ≈ 0.001
        @test qflx_rain_plus_snomelt[1] ≈ 0.001 + (1.0 - 0.8) * 0.002

        # No-snow column
        @test qflx_snow_drain[2] ≈ 0.001
        @test qflx_rain_plus_snomelt[2] ≈ 0.003 + 0.001
    end
end
