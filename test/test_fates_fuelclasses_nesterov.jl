# test_fates_fuelclasses_nesterov.jl
# Tests for FATES Tier-F Batch-2 fire modules:
#   - FatesFuelClassesMod (fuel-class enumeration / accessors)
#   - SFNesterovMod (Nesterov fire-weather index)

@testset "FATES FatesFuelClassesMod" begin
    # Total fuel-class count.
    @test CLM.num_fuel_classes == 6

    # Array indices (canonical const values).
    @test CLM.twigs_i == 1
    @test CLM.small_branches_i == 2
    @test CLM.large_branches_i == 3
    @test CLM.trunks_i == 4
    @test CLM.dead_leaves_i == 5
    @test CLM.live_grass_i == 6

    # Accessor methods on the module-level fuel_classes instance match the
    # consts (Fortran fuel_classes%twigs() etc.).
    fc = CLM.fuel_classes
    @test CLM.twigs(fc) == 1
    @test CLM.small_branches(fc) == 2
    @test CLM.large_branches(fc) == 3
    @test CLM.trunks(fc) == 4
    @test CLM.dead_leaves(fc) == 5
    @test CLM.live_grass(fc) == 6

    # All six indices are distinct and cover 1:num_fuel_classes.
    idxs = [CLM.twigs(fc), CLM.small_branches(fc), CLM.large_branches(fc),
            CLM.trunks(fc), CLM.dead_leaves(fc), CLM.live_grass(fc)]
    @test sort(idxs) == collect(1:CLM.num_fuel_classes)
end

@testset "FATES SFNesterovMod" begin
    # nesterov_index is a fire_weather subtype.
    @test CLM.nesterov_index <: CLM.fire_weather

    # Threshold constant.
    @test CLM.min_precip_thresh == 3.0

    # --- helper consistency: dewpoint + calc_nesterov_index ----------------
    # dewpoint formula (Lawrence 2005, Eq. 8) using FatesConstants dewpoint_a/b.
    let temp_C = 25.0, rh = 50.0
        a = CLM.dewpoint_a
        b = CLM.dewpoint_b
        yip = log(max(1.0, rh) / 100.0) + (a * temp_C) / (b + temp_C)
        expected_td = (b * yip) / (a - yip)
        @test CLM.dewpoint(temp_C, rh) ≈ expected_td

        # calc_nesterov_index = (T - Td)*T, clamped >= 0.
        @test CLM.calc_nesterov_index(temp_C, expected_td) ≈ (temp_C - expected_td) * temp_C
    end

    # calc_nesterov_index clamps negatives (Td > T -> negative product) to 0.
    @test CLM.calc_nesterov_index(10.0, 30.0) == 0.0

    # --- init ---------------------------------------------------------------
    ni = CLM.nesterov_index()
    @test ni.fire_weather_index == 0.0
    @test ni.effective_windspeed == 0.0
    # mutate then re-init.
    ni.fire_weather_index = 99.0
    ni.effective_windspeed = 5.0
    CLM.init_fire_weather!(ni)
    @test ni.fire_weather_index == 0.0
    @test ni.effective_windspeed == 0.0

    # --- accumulation on dry days ------------------------------------------
    # Synthetic dry-day weather sequence (precip below threshold) should
    # accumulate the daily Nesterov index.
    ni2 = CLM.nesterov_index()
    CLM.init_fire_weather!(ni2)

    dry_days = [(30.0, 0.0, 40.0),   # (temp_C, precip, rh)
                (28.0, 1.5, 55.0),
                (32.0, 0.0, 30.0)]
    acc = 0.0
    for (tC, pr, rh) in dry_days
        td = CLM.dewpoint(tC, rh)
        acc += CLM.calc_nesterov_index(tC, td)
        CLM.update_index!(ni2, tC, pr, rh, 100.0)
        @test ni2.fire_weather_index ≈ acc
    end
    @test ni2.fire_weather_index > 0.0  # warm dry days accumulate positive NI
    acc_before_rain = ni2.fire_weather_index

    # --- reset on a rain day above threshold --------------------------------
    # precip > min_precip_thresh (3.0 mm/day) rezeroes the index regardless of
    # the previous accumulation.
    CLM.update_index!(ni2, 25.0, 5.0, 60.0, 100.0)
    @test acc_before_rain > 0.0
    @test ni2.fire_weather_index == 0.0

    # A dry day right after the reset starts accumulating again from zero.
    td_after = CLM.dewpoint(27.0, 45.0)
    CLM.update_index!(ni2, 27.0, 0.0, 45.0, 100.0)
    @test ni2.fire_weather_index ≈ CLM.calc_nesterov_index(27.0, td_after)

    # Exactly at the threshold (precip == 3.0) does NOT reset (strict >).
    ni3 = CLM.nesterov_index()
    CLM.init_fire_weather!(ni3)
    CLM.update_index!(ni3, 30.0, 3.0, 40.0, 100.0)
    @test ni3.fire_weather_index > 0.0

    # --- base method still works on the subtype -----------------------------
    ni4 = CLM.nesterov_index()
    CLM.update_effective_windspeed!(ni4, 10.0, 1.0, 0.0, 0.0)
    @test ni4.effective_windspeed ≈ 10.0 * CLM.wind_atten_treed
end
