#!/usr/bin/env julia
# ==========================================================================
# Enzyme.jl Feasibility Study for CLM.jl
#
# Purpose: Investigate whether Enzyme.jl reverse-mode AD can work with
# CLM.jl physics kernels. Documents blockers and limitations.
# This is an INVESTIGATION only — not integrated into calibration.
#
# Known potential blockers:
# - Ref globals (SNOW_DZMIN, etc.)
# - Module-level const params structs
# - LAPACK calls in band_diagonal solver
# - try/catch blocks in _calib_dual_copy
# - NCDatasets I/O (not differentiable)
# ==========================================================================
using Test
using CLM

println("=" ^ 70)
println("ENZYME FEASIBILITY STUDY for CLM.jl")
println("=" ^ 70)

# Check if Enzyme is available
enzyme_available = try
    @eval using Enzyme
    true
catch e
    println("  Enzyme.jl not available: $e")
    println("  To install: ] add Enzyme")
    false
end

@testset "Enzyme Feasibility" begin

    @testset "Enzyme availability" begin
        if !enzyme_available
            @test_skip "Enzyme.jl not installed — skipping all Enzyme tests"
            return
        end
        @test enzyme_available
        println("  Enzyme.jl loaded successfully")
    end

    if !enzyme_available
        return
    end

    # ------------------------------------------------------------------
    # Test 1: qsat (pure math, no globals)
    # ------------------------------------------------------------------
    @testset "qsat — pure function" begin
        blockers = String[]
        success = false
        try
            function qsat_wrapper(T_in)
                es, esdT, qs, qsdT = CLM.qsat(T_in, 85000.0)
                return es
            end

            dT = Enzyme.autodiff(Enzyme.Reverse, qsat_wrapper, Enzyme.Active, Enzyme.Active(280.0))
            @test isfinite(dT[1][1])
            success = true
            println("  qsat: PASS — Enzyme reverse-mode works")
            println("    d(es)/dT at 280K = $(dT[1][1])")
        catch e
            push!(blockers, "qsat: $(sprint(showerror, e))")
            @test_broken false
            println("  qsat: BLOCKED — $(sprint(showerror, e))")
        end
    end

    # ------------------------------------------------------------------
    # Test 2: smooth_* primitives
    # ------------------------------------------------------------------
    @testset "smooth_ad primitives" begin
        try
            # smooth_max
            dx = Enzyme.autodiff(Enzyme.Reverse,
                x -> CLM.smooth_max(x, 0.0),
                Enzyme.Active, Enzyme.Active(1.0))
            @test isfinite(dx[1][1])

            # smooth_min
            dx2 = Enzyme.autodiff(Enzyme.Reverse,
                x -> CLM.smooth_min(x, 5.0),
                Enzyme.Active, Enzyme.Active(3.0))
            @test isfinite(dx2[1][1])

            # smooth_ifelse
            dx3 = Enzyme.autodiff(Enzyme.Reverse,
                x -> CLM.smooth_ifelse(x - 273.15, 1.0, 0.0),
                Enzyme.Active, Enzyme.Active(280.0))
            @test isfinite(dx3[1][1])

            println("  smooth_* primitives: PASS")
        catch e
            @test_broken false
            println("  smooth_* primitives: BLOCKED — $(sprint(showerror, e))")
        end
    end

    # ------------------------------------------------------------------
    # Test 3: tridiagonal solver
    # ------------------------------------------------------------------
    @testset "tridiagonal solver" begin
        try
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
            println("  tridiagonal solver: PASS")
        catch e
            @test_broken false
            println("  tridiagonal solver: BLOCKED — $(sprint(showerror, e))")
        end
    end

    # ------------------------------------------------------------------
    # Test 4: band_diagonal solver (uses LAPACK for Float64)
    # ------------------------------------------------------------------
    @testset "band_diagonal solver" begin
        try
            function band_test(x)
                nc = 1; n = 5
                a = zeros(nc, n); b = fill(2.0, nc, n); c = zeros(nc, n)
                r = ones(nc, n)
                b[1, 3] = x
                for j in 2:n; a[1, j] = -0.5; end
                for j in 1:n-1; c[1, j] = -0.5; end
                CLM.tridiagonal_solve_multi!(a, b, c, r, 1, nc, n)
                return r[1, 3]
            end

            dx = Enzyme.autodiff(Enzyme.Reverse, band_test,
                Enzyme.Active, Enzyme.Active(2.0))
            @test isfinite(dx[1][1])
            println("  band_diagonal solver: PASS")
        catch e
            @test_broken false
            println("  band_diagonal solver: BLOCKED — $(sprint(showerror, e))")
        end
    end

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    println("\n" * "=" ^ 70)
    println("ENZYME FEASIBILITY SUMMARY")
    println("=" ^ 70)
    println("""
    Findings:
    1. Enzyme reverse-mode on pure math functions: likely works
    2. Key blockers for full CLM integration:
       a. Module-level const params (canopy_fluxes_params, params_inst)
          — Enzyme cannot differentiate through global mutable state
       b. LAPACK calls in band_diagonal for Float64
          — Would need pure-Julia fallback (already exists for non-Float64)
       c. NCDatasets I/O in initialization — not differentiable
          — Would need to separate init from physics
       d. try/catch in _calib_dual_copy — Enzyme doesn't support exceptions
       e. Ref globals (SNOW_DZMIN etc.) — read-only, likely OK
    3. Recommendation: Enzyme integration is feasible for individual kernels
       but requires significant refactoring for full model. ForwardDiff
       remains the practical choice for now.
    """)
end
