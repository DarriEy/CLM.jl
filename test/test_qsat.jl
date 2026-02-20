@testset "QSat - Saturation Vapor Pressure" begin

    @testset "Known values at standard conditions" begin
        qs, es, dqsdT, desdT = CLM.qsat(273.15, 101325.0)
        @test es ≈ 611.2 atol=1.0
        @test qs > 0.0
        @test qs < 1.0
        @test dqsdT > 0.0
        @test desdT > 0.0
    end

    @testset "Boiling point" begin
        qs, es, dqsdT, desdT = CLM.qsat(373.15, 101325.0)
        @test es ≈ 101325.0 rtol=0.05
    end

    @testset "Monotonicity in temperature" begin
        p = 101325.0
        temps = [250.0, 260.0, 270.0, 280.0, 290.0, 300.0, 310.0]
        es_vals = [CLM.qsat(T, p)[2] for T in temps]
        qs_vals = [CLM.qsat(T, p)[1] for T in temps]

        for i in 2:length(temps)
            @test es_vals[i] > es_vals[i-1]
            @test qs_vals[i] > qs_vals[i-1]
        end
    end

    @testset "Ice vs water transition at 0°C" begin
        p = 101325.0
        _, es_ice, _, _ = CLM.qsat(272.15, p)
        _, es_water, _, _ = CLM.qsat(274.15, p)

        @test es_ice > 500.0
        @test es_ice < 700.0
        @test es_water > 600.0
        @test es_water < 800.0
    end

    @testset "Pressure dependence" begin
        T = 300.0
        qs_low, es_low, _, _ = CLM.qsat(T, 80000.0)
        qs_high, es_high, _, _ = CLM.qsat(T, 101325.0)

        @test es_low ≈ es_high
        @test qs_low > qs_high
    end

    @testset "Cold temperatures (Arctic)" begin
        qs, es, dqsdT, desdT = CLM.qsat(233.15, 101325.0)
        @test es > 0.0
        @test es < 100.0
        @test qs > 0.0
        @test qs < 0.001
    end

    @testset "No-derivatives version matches" begin
        T, p = 290.0, 95000.0
        qs1, es1, _, _ = CLM.qsat(T, p)
        qs2, es2 = CLM.qsat_no_derivs(T, p)
        @test qs1 ≈ qs2
        @test es1 ≈ es2
    end

    @testset "Derivative vs finite difference (water)" begin
        T = 290.0
        p = 95000.0
        dT = 1.0e-6

        qs, es, dqsdT, desdT = CLM.qsat(T, p)
        qs_plus, es_plus, _, _ = CLM.qsat(T + dT, p)
        qs_minus, es_minus, _, _ = CLM.qsat(T - dT, p)

        desdT_fd = (es_plus - es_minus) / (2.0 * dT)
        dqsdT_fd = (qs_plus - qs_minus) / (2.0 * dT)

        @test desdT ≈ desdT_fd rtol=1e-4
        @test dqsdT ≈ dqsdT_fd rtol=1e-4
    end

    @testset "Derivative vs finite difference (ice)" begin
        T = 260.0
        p = 95000.0
        dT = 1.0e-6

        qs, es, dqsdT, desdT = CLM.qsat(T, p)
        qs_plus, es_plus, _, _ = CLM.qsat(T + dT, p)
        qs_minus, es_minus, _, _ = CLM.qsat(T - dT, p)

        desdT_fd = (es_plus - es_minus) / (2.0 * dT)
        dqsdT_fd = (qs_plus - qs_minus) / (2.0 * dT)

        @test desdT ≈ desdT_fd rtol=1e-4
        @test dqsdT ≈ dqsdT_fd rtol=1e-4
    end

end
