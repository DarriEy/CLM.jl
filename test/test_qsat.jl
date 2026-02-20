@testset "QSat - Saturation Vapor Pressure" begin

    @testset "Known values at standard conditions" begin
        # At 0°C (273.15 K), 101325 Pa
        # es should be ~611 Pa (triple point of water)
        qs, es, dqsdT, desdT = CLM.qsat(273.15, 101325.0)
        @test es ≈ 611.2 atol=1.0  # ~611 Pa at 0°C
        @test qs > 0.0
        @test qs < 1.0
        @test dqsdT > 0.0  # qs always increases with T
        @test desdT > 0.0  # es always increases with T
    end

    @testset "Boiling point" begin
        # At 100°C (373.15 K), es should be ~101325 Pa
        qs, es, dqsdT, desdT = CLM.qsat(373.15, 101325.0)
        @test es ≈ 101325.0 rtol=0.05  # should be near 1 atm
    end

    @testset "Monotonicity in temperature" begin
        p = 101325.0
        temps = [250.0, 260.0, 270.0, 280.0, 290.0, 300.0, 310.0]
        es_vals = [CLM.qsat(T, p)[2] for T in temps]
        qs_vals = [CLM.qsat(T, p)[1] for T in temps]

        # Both es and qs should increase monotonically with T
        for i in 2:length(temps)
            @test es_vals[i] > es_vals[i-1]
            @test qs_vals[i] > qs_vals[i-1]
        end
    end

    @testset "Ice vs water transition at 0°C" begin
        p = 101325.0
        # Just below and above freezing
        _, es_ice, _, _ = CLM.qsat(272.15, p)   # -1°C (ice)
        _, es_water, _, _ = CLM.qsat(274.15, p) # +1°C (water)

        # Both should be near 611 Pa but from different polynomials
        @test es_ice > 500.0
        @test es_ice < 700.0
        @test es_water > 600.0
        @test es_water < 800.0
    end

    @testset "Pressure dependence" begin
        T = 300.0  # 27°C
        # Higher pressure → lower qs (same es, but qs = 0.622 * es / (p - 0.378*es))
        qs_low, es_low, _, _ = CLM.qsat(T, 80000.0)
        qs_high, es_high, _, _ = CLM.qsat(T, 101325.0)

        @test es_low ≈ es_high  # es doesn't depend on pressure
        @test qs_low > qs_high  # lower pressure → higher qs
    end

    @testset "Cold temperatures (Arctic)" begin
        # Very cold (-40°C = 233.15 K) — should still work
        qs, es, dqsdT, desdT = CLM.qsat(233.15, 101325.0)
        @test es > 0.0
        @test es < 100.0  # very low vapor pressure
        @test qs > 0.0
        @test qs < 0.001  # very low specific humidity
    end

    @testset "No-derivatives version matches" begin
        T, p = 290.0, 95000.0
        qs1, es1, _, _ = CLM.qsat(T, p)
        qs2, es2 = CLM.qsat_no_derivs(T, p)
        @test qs1 ≈ qs2
        @test es1 ≈ es2
    end

end
