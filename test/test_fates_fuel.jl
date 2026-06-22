# test_fates_fuel.jl
# Tests for FATES Tier-F Batch-3 fire module FatesFuelMod (SPITFIRE fuel_type):
# loading accumulation, non-trunk summation, fractional loading, bulk-density /
# SAV / moisture-of-extinction aggregation, and Nesterov fuel moisture.

@testset "FATES FatesFuelMod" begin
    fc = CLM.fuel_classes
    n  = CLM.num_fuel_classes

    # ---- construction + Init ---------------------------------------------
    fuel = CLM.fuel_type()
    @test length(fuel.loading) == n
    @test all(fuel.loading .== 0.0)
    @test fuel.non_trunk_loading == 0.0

    # set some non-zero state then Init -> everything back to zero
    fuel.loading .= 1.0
    fuel.non_trunk_loading = 5.0
    fuel.SAV_notrunks = 3.0
    CLM.init_fuel!(fuel)
    @test all(fuel.loading .== 0.0)
    @test all(fuel.frac_loading .== 0.0)
    @test all(fuel.frac_burnt .== 0.0)
    @test all(fuel.effective_moisture .== 0.0)
    @test fuel.non_trunk_loading == 0.0
    @test fuel.SAV_notrunks == 0.0
    @test fuel.MEF_notrunks == 0.0

    # ---- UpdateLoading: pools land in the right class slots --------------
    leaf, twig, smb, lgb, trunk, grass = 1.0, 2.0, 3.0, 4.0, 10.0, 5.0
    CLM.update_loading!(fuel, leaf, twig, smb, lgb, trunk, grass)
    @test fuel.loading[CLM.dead_leaves(fc)]    == leaf
    @test fuel.loading[CLM.twigs(fc)]          == twig
    @test fuel.loading[CLM.small_branches(fc)] == smb
    @test fuel.loading[CLM.large_branches(fc)] == lgb
    @test fuel.loading[CLM.trunks(fc)]         == trunk
    @test fuel.loading[CLM.live_grass(fc)]     == grass

    # ---- SumLoading: total excludes trunks -------------------------------
    CLM.sum_loading!(fuel)
    expected_nontrunk = leaf + twig + smb + lgb + grass   # trunk EXCLUDED
    @test fuel.non_trunk_loading ≈ expected_nontrunk
    @test fuel.non_trunk_loading ≈ sum(fuel.loading) - fuel.loading[CLM.trunks(fc)]
    # trunk is genuinely excluded (sanity: total with trunk differs)
    @test sum(fuel.loading) > fuel.non_trunk_loading

    # ---- CalculateFractionalLoading: fractions sum to 1, trunk is 0 ------
    CLM.calculate_fractional_loading!(fuel)
    @test fuel.frac_loading[CLM.trunks(fc)] == 0.0
    @test sum(fuel.frac_loading) ≈ 1.0
    @test fuel.frac_loading[CLM.dead_leaves(fc)] ≈ leaf / expected_nontrunk
    @test fuel.frac_loading[CLM.twigs(fc)]       ≈ twig / expected_nontrunk

    # zero-loading edge case -> all fractions zero, non_trunk_loading zero
    fuel0 = CLM.fuel_type()
    CLM.calculate_fractional_loading!(fuel0)
    @test all(fuel0.frac_loading .== 0.0)
    @test fuel0.non_trunk_loading == 0.0

    # ---- MoistureOfExtinction: MEF_a - MEF_b*log(sav) --------------------
    MEF_a, MEF_b = 0.524, 0.066
    for sav in (0.355, 0.44, 0.525, 0.248)
        @test CLM.moisture_of_extinction_fn(sav) ≈ MEF_a - MEF_b * log(sav)
    end

    # ---- AverageBulkDensity_NoTrunks -------------------------------------
    # per-class bulk density; non-trunk weighted average by frac_loading.
    bulk_density = [10.0, 20.0, 30.0, 99.0, 40.0, 50.0]  # trunk=99 must be excluded
    CLM.average_bulk_density_notrunks!(fuel, bulk_density)
    expected_bd = 0.0
    for i in 1:n
        i == CLM.trunks(fc) && continue
        expected_bd += fuel.frac_loading[i] * bulk_density[i]
    end
    @test fuel.bulk_density_notrunks ≈ expected_bd
    # trunk's large bulk density did not leak in (weighted avg < trunk value)
    @test fuel.bulk_density_notrunks < bulk_density[CLM.trunks(fc)]

    # no-loading fallback: plain mean over ALL classes
    fuel_bd0 = CLM.fuel_type()  # non_trunk_loading == 0
    CLM.average_bulk_density_notrunks!(fuel_bd0, bulk_density)
    @test fuel_bd0.bulk_density_notrunks ≈ sum(bulk_density) / n

    # ---- AverageSAV_NoTrunks ---------------------------------------------
    sav_fuel = [0.355, 0.44, 0.525, 0.63, 0.248, 0.248]  # twigs..live_grass; trunk excl
    CLM.average_sav_notrunks!(fuel, sav_fuel)
    expected_sav = 0.0
    for i in 1:n
        i == CLM.trunks(fc) && continue
        expected_sav += fuel.frac_loading[i] * sav_fuel[i]
    end
    @test fuel.SAV_notrunks ≈ expected_sav

    fuel_sav0 = CLM.fuel_type()
    CLM.average_sav_notrunks!(fuel_sav0, sav_fuel)
    @test fuel_sav0.SAV_notrunks ≈ sum(sav_fuel) / n

    # ---- CalculateFuelMoistureNesterov: moisture = exp(-alpha*NI) --------
    drying_ratio = 10.0
    NI = 5.0
    moisture = zeros(Float64, n)
    CLM.calculate_fuel_moisture_nesterov(sav_fuel, drying_ratio, NI, moisture)
    for i in 1:n
        sav_used = (i == CLM.live_grass(fc)) ? sav_fuel[CLM.twigs(fc)] : sav_fuel[i]
        @test moisture[i] ≈ exp(-1.0 * (sav_used / drying_ratio) * NI)
    end
    # live grass uses twig SAV -> moisture matches twigs class
    @test moisture[CLM.live_grass(fc)] ≈ moisture[CLM.twigs(fc)]

    # ---- UpdateFuelMoisture end-to-end with a Nesterov fire_weather ------
    nest = CLM.nesterov_index()
    nest.fire_weather_index = NI
    CLM.update_fuel_moisture!(fuel, sav_fuel, drying_ratio, nest)
    # effective_moisture = moisture / MEF
    for i in 1:n
        mef_i = CLM.moisture_of_extinction_fn(sav_fuel[i])
        sav_used = (i == CLM.live_grass(fc)) ? sav_fuel[CLM.twigs(fc)] : sav_fuel[i]
        moist_i = exp(-1.0 * (sav_used / drying_ratio) * NI)
        @test fuel.effective_moisture[i] ≈ moist_i / mef_i
    end
    # non-trunk weighted averages
    exp_avg_moist = 0.0
    exp_mef = 0.0
    for i in 1:n
        i == CLM.trunks(fc) && continue
        sav_used = (i == CLM.live_grass(fc)) ? sav_fuel[CLM.twigs(fc)] : sav_fuel[i]
        moist_i = exp(-1.0 * (sav_used / drying_ratio) * NI)
        exp_avg_moist += fuel.frac_loading[i] * moist_i
        exp_mef += fuel.frac_loading[i] * CLM.moisture_of_extinction_fn(sav_fuel[i])
    end
    @test fuel.average_moisture_notrunks ≈ exp_avg_moist
    @test fuel.MEF_notrunks ≈ exp_mef

    # zero total loading -> moisture diagnostics zeroed
    fuel_m0 = CLM.fuel_type()
    CLM.update_fuel_moisture!(fuel_m0, sav_fuel, drying_ratio, nest)
    @test all(fuel_m0.effective_moisture .== 0.0)
    @test fuel_m0.average_moisture_notrunks == 0.0
    @test fuel_m0.MEF_notrunks == 0.0

    # ---- Fuse: area-weighted blend of two fuel objects -------------------
    a = CLM.fuel_type(); a.loading .= 2.0; a.non_trunk_loading = 8.0
    b = CLM.fuel_type(); b.loading .= 6.0; b.non_trunk_loading = 24.0
    self_area, donor_area = 1.0, 3.0
    sw = self_area / (self_area + donor_area)  # 0.25
    CLM.fuse!(a, self_area, donor_area, b)
    @test all(a.loading .≈ 2.0 * sw .+ 6.0 * (1.0 - sw))  # = 5.0
    @test a.non_trunk_loading ≈ 8.0 * sw + 24.0 * (1.0 - sw)  # = 20.0
end
