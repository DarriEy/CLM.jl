#!/usr/bin/env julia
# ==========================================================================
# Enzyme.jl Smoke Tests for CLM.jl
#
# Tests Enzyme reverse-mode AD on CLM primitives and small kernels.
# Verifies gradients match ForwardDiff within tolerance.
# ==========================================================================

using Test
using CLM

println("=" ^ 70)
println("ENZYME SMOKE TESTS for CLM.jl")
println("=" ^ 70)

enzyme_available = try
    @eval using Enzyme
    true
catch e
    println("  Enzyme.jl not available: $e")
    false
end

fd_available = try
    @eval using ForwardDiff
    true
catch
    false
end

@testset "Enzyme Smoke Tests" begin

    if !enzyme_available
        @test_skip "Enzyme.jl not installed"
        return
    end

    # ------------------------------------------------------------------
    # Test 1: qsat derivatives
    # ------------------------------------------------------------------
    @testset "qsat reverse-mode" begin
        function qsat_es(T)
            es, _, _, _ = CLM.qsat(T, 85000.0)
            return es
        end

        dT_enz = Enzyme.autodiff(Enzyme.Reverse, qsat_es, Enzyme.Active, Enzyme.Active(280.0))
        @test isfinite(dT_enz[1][1])

        if fd_available
            dT_fd = ForwardDiff.derivative(qsat_es, 280.0)
            rel_err = abs(dT_enz[1][1] - dT_fd) / abs(dT_fd)
            println("  qsat: Enzyme=$(dT_enz[1][1]), FD=$(dT_fd), err=$(round(rel_err*100, digits=2))%")
            @test rel_err < 0.01
        end
    end

    # ------------------------------------------------------------------
    # Test 2: smooth_min / smooth_max
    # ------------------------------------------------------------------
    @testset "smooth_min/max reverse-mode" begin
        # Force smooth mode so Float64 goes through LogSumExp
        CLM.SMOOTH_MODE[] = :always

        dx = Enzyme.autodiff(Enzyme.Reverse,
            x -> CLM.smooth_min(x, 5.0), Enzyme.Active, Enzyme.Active(3.0))
        @test isfinite(dx[1][1])
        @test abs(dx[1][1] - 1.0) < 0.01  # d/dx min(x, 5) ≈ 1 when x < 5

        dx2 = Enzyme.autodiff(Enzyme.Reverse,
            x -> CLM.smooth_max(x, 0.0), Enzyme.Active, Enzyme.Active(1.0))
        @test isfinite(dx2[1][1])
        @test abs(dx2[1][1] - 1.0) < 0.01  # d/dx max(x, 0) ≈ 1 when x > 0

        CLM.SMOOTH_MODE[] = :auto
        println("  smooth_min/max: PASS")
    end

    # ------------------------------------------------------------------
    # Test 3: smooth_ifelse
    # ------------------------------------------------------------------
    @testset "smooth_ifelse reverse-mode" begin
        CLM.SMOOTH_MODE[] = :always

        dx = Enzyme.autodiff(Enzyme.Reverse,
            x -> CLM.smooth_ifelse(x - 273.15, x * 2.0, x * 0.5),
            Enzyme.Active, Enzyme.Active(300.0))
        @test isfinite(dx[1][1])

        CLM.SMOOTH_MODE[] = :auto
        println("  smooth_ifelse: PASS")
    end

    # ------------------------------------------------------------------
    # Test 4: smooth_heaviside
    # ------------------------------------------------------------------
    @testset "smooth_heaviside reverse-mode" begin
        CLM.SMOOTH_MODE[] = :always

        dx = Enzyme.autodiff(Enzyme.Reverse,
            x -> CLM.smooth_heaviside(x), Enzyme.Active, Enzyme.Active(0.0))
        @test isfinite(dx[1][1])
        @test dx[1][1] > 0  # sigmoid derivative is positive at x=0

        CLM.SMOOTH_MODE[] = :auto
        println("  smooth_heaviside: PASS")
    end

    # ------------------------------------------------------------------
    # Test 5: smooth_abs
    # ------------------------------------------------------------------
    @testset "smooth_abs reverse-mode" begin
        CLM.SMOOTH_MODE[] = :always

        dx = Enzyme.autodiff(Enzyme.Reverse,
            x -> CLM.smooth_abs(x), Enzyme.Active, Enzyme.Active(2.0))
        @test isfinite(dx[1][1])
        @test abs(dx[1][1] - 1.0) < 0.01  # d|x|/dx ≈ 1 for x > 0

        CLM.SMOOTH_MODE[] = :auto
        println("  smooth_abs: PASS")
    end

    # ------------------------------------------------------------------
    # Test 6: tridiagonal solver
    # ------------------------------------------------------------------
    @testset "tridiagonal solver reverse-mode" begin
        function tridiag_test(x)
            n = 5
            a = zeros(n); b = fill(2.0, n); c = zeros(n)
            r = ones(n); u = zeros(n)
            b[3] = x
            for i in 2:n; a[i] = -0.5; end
            for i in 1:n-1; c[i] = -0.5; end
            CLM.tridiagonal_solve!(u, a, b, c, r, 1, n)
            return u[3]
        end

        dx = Enzyme.autodiff(Enzyme.Reverse, tridiag_test,
            Enzyme.Active, Enzyme.Active(2.0))
        @test isfinite(dx[1][1])

        if fd_available
            dx_fd = ForwardDiff.derivative(tridiag_test, 2.0)
            rel_err = abs(dx[1][1] - dx_fd) / abs(dx_fd)
            println("  tridiag: Enzyme=$(dx[1][1]), FD=$(dx_fd), err=$(round(rel_err*100, digits=2))%")
            @test rel_err < 0.01
        end
    end

    # ------------------------------------------------------------------
    # Test 7: Composite kernel (qsat + smooth chain)
    # ------------------------------------------------------------------
    @testset "composite kernel reverse-mode" begin
        CLM.SMOOTH_MODE[] = :always

        function composite(T)
            es, _, qs, _ = CLM.qsat(T, 85000.0)
            clamped = CLM.smooth_clamp(qs, 0.001, 0.1)
            return clamped * es
        end

        dx = Enzyme.autodiff(Enzyme.Reverse, composite,
            Enzyme.Active, Enzyme.Active(285.0))
        @test isfinite(dx[1][1])

        if fd_available
            dx_fd = ForwardDiff.derivative(composite, 285.0)
            rel_err = abs(dx[1][1] - dx_fd) / max(abs(dx_fd), 1e-20)
            println("  composite: Enzyme=$(dx[1][1]), FD=$(dx_fd), err=$(round(rel_err*100, digits=2))%")
            @test rel_err < 0.01
        end

        CLM.SMOOTH_MODE[] = :auto
    end

end
