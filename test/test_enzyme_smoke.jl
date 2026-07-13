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
        @info "Enzyme.jl not installed — skipping Enzyme Smoke Tests"
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
        # smooth_* are type-based (GPU-safe): Float64 is exact, Dual smooths.
        # Verify the smooth sigmoid (positive derivative at 0) where it applies —
        # on a Dual via ForwardDiff.
        if fd_available
            d = ForwardDiff.derivative(CLM.smooth_heaviside, 0.0)
            @test isfinite(d)
            @test d > 0  # sigmoid derivative is positive at x=0
        end
        # Enzyme differentiates the exact Float64 heaviside finitely.
        dx = Enzyme.autodiff(Enzyme.Reverse,
            x -> CLM.smooth_heaviside(x), Enzyme.Active, Enzyme.Active(0.5))
        @test isfinite(dx[1][1])
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
            T = typeof(x)  # match element type so ForwardDiff.Dual cross-check works too
            a = zeros(T, n); b = fill(T(2.0), n); c = zeros(T, n)
            r = ones(T, n); u = zeros(T, n)
            b[3] = x
            for i in 2:n; a[i] = -0.5; end
            for i in 1:n-1; c[i] = -0.5; end
            CLM.tridiagonal_solve!(u, a, b, c, r, 1, n)
            return u[3]
        end

        # The gradient must be RIGHT, not merely finite. `tridiagonal_solve!` carries a
        # custom Enzyme adjoint rule (src/infrastructure/enzyme_rules.jl) precisely
        # because the generic reverse pass reads an unzeroed shadow of the solver's
        # in-callee Thomas workspace on Julia 1.10 — that produced a wrong-but-finite
        # gradient on every call and a NaN/Inf on ~2-5% of calls (the intermittent CI
        # failure this test used to have). An `isfinite` check alone could not see the
        # wrong value, so assert against ForwardDiff on EVERY version instead.
        dx_fd = fd_available ? ForwardDiff.derivative(tridiag_test, 2.0) : NaN

        dx = Enzyme.autodiff(Enzyme.Reverse, tridiag_test,
            Enzyme.Active, Enzyme.Active(2.0))
        @test isfinite(dx[1][1])

        if fd_available
            rel_err = abs(dx[1][1] - dx_fd) / abs(dx_fd)
            println("  tridiag: Enzyme=$(dx[1][1]), FD=$(dx_fd), err=$(round(rel_err*100, digits=2))%")
            @test rel_err < 1.0e-8

            # Repeat guard: an unzeroed shadow / uninitialized workspace is a
            # *nondeterministic* defect — it depends on what the allocator hands back,
            # so a single call passes most of the time. Hammer it. Before the adjoint
            # rule this loop failed within a few dozen iterations on Julia 1.10.
            reps = 200
            worst = 0.0
            nonfinite = 0
            for _ in 1:reps
                d = Enzyme.autodiff(Enzyme.Reverse, tridiag_test,
                    Enzyme.Active, Enzyme.Active(2.0))[1][1]
                isfinite(d) || (nonfinite += 1; continue)
                worst = max(worst, abs(d - dx_fd) / abs(dx_fd))
            end
            println("  tridiag: $(reps) reps — nonfinite=$(nonfinite), worst rel err=$(worst)")
            @test nonfinite == 0
            @test worst < 1.0e-8
        end
    end

    # ------------------------------------------------------------------
    # Test 6b: tridiagonal adjoint rule — exactness of ALL band gradients
    #
    # The custom reverse rule hand-derives  λ = A⁻ᵀ ū ;  r̄ += λ ;  Ā = -λ uᵀ.
    # Validate every band (a, b, c, r) against ForwardDiff AND central finite
    # differences, with jtop > 1 so the offset indexing is exercised, and check
    # that the structurally-unused entries (a[1:jtop], c[nlevs]) come back exactly
    # zero. A hand-written adjoint that is merely finite is worthless.
    # ------------------------------------------------------------------
    @testset "tridiagonal adjoint rule (all bands, jtop>1)" begin
        if !fd_available
            @info "ForwardDiff not available — skipping tridiagonal adjoint-rule check"
        else
            n = 9
            jtop = 3
            a0 = vcat(zeros(jtop), [0.31, -0.22, 0.44, -0.17, 0.28, -0.39])  # a[1:jtop] unused
            b0 = [3.4, 3.1, 3.9, 3.2, 3.7, 3.3, 3.5, 3.8, 3.6]               # diag. dominant
            c0 = vcat([-0.25, 0.33, -0.41, 0.19, -0.29, 0.37, -0.21, 0.24], 0.0)  # c[n] unused
            r0 = [0.7, -1.2, 0.4, 1.9, -0.6, 0.3, -1.1, 0.8, 0.2]
            w  = [0.9, -0.4, 1.3, 0.2, -1.7, 0.6, 0.5, -0.8, 1.1]

            function obj(a, b, c, r)
                T = promote_type(eltype(a), eltype(b), eltype(c), eltype(r))
                u = zeros(T, n)
                CLM.tridiagonal_solve!(u, a, b, c, r, jtop, n)
                s = zero(T)
                for j in jtop:n
                    s += w[j] * u[j]
                end
                return s
            end

            da = zeros(n); db = zeros(n); dc = zeros(n); dr = zeros(n)
            Enzyme.autodiff(Enzyme.Reverse, obj, Enzyme.Active,
                Enzyme.Duplicated(copy(a0), da), Enzyme.Duplicated(copy(b0), db),
                Enzyme.Duplicated(copy(c0), dc), Enzyme.Duplicated(copy(r0), dr))

            ga = ForwardDiff.gradient(v -> obj(v, b0, c0, r0), a0)
            gb = ForwardDiff.gradient(v -> obj(a0, v, c0, r0), b0)
            gc = ForwardDiff.gradient(v -> obj(a0, b0, v, r0), c0)
            gr = ForwardDiff.gradient(v -> obj(a0, b0, c0, v), r0)

            function fdgrad(f, v0)
                g = similar(v0)
                for i in eachindex(v0)
                    h = 1.0e-6 * max(abs(v0[i]), 1.0)
                    vp = copy(v0); vp[i] += h
                    vm = copy(v0); vm[i] -= h
                    g[i] = (f(vp) - f(vm)) / (2h)
                end
                return g
            end
            fa = fdgrad(v -> obj(v, b0, c0, r0), a0)
            fb = fdgrad(v -> obj(a0, v, c0, r0), b0)
            fc = fdgrad(v -> obj(a0, b0, v, r0), c0)
            fr = fdgrad(v -> obj(a0, b0, c0, v), r0)

            relerr(x, y) = maximum(abs.(x .- y)) / max(maximum(abs.(y)), 1.0e-12)
            for (nm, e, g, f) in (("a", da, ga, fa), ("b", db, gb, fb),
                                  ("c", dc, gc, fc), ("r", dr, gr, fr))
                println("  tridiag adjoint d/d$nm: vs ForwardDiff=$(relerr(e, g)), vs numFD=$(relerr(e, f))")
                @test all(isfinite, e)
                @test relerr(e, g) < 1.0e-10   # exact adjoint vs forward-mode AD
                @test relerr(e, f) < 1.0e-6    # and vs an AD-independent numerical FD
            end
            # Bands the Thomas sweep never reads must carry exactly zero derivative.
            @test all(iszero, da[1:jtop])
            @test iszero(dc[n])
        end
    end

    # ------------------------------------------------------------------
    # Test 7: Composite kernel (qsat + smooth chain)
    # ------------------------------------------------------------------
    @testset "composite kernel reverse-mode" begin
        function composite(T)
            es, _, qs, _ = CLM.qsat(T, 85000.0)
            clamped = CLM.smooth_clamp(qs, 0.001, 0.1)
            return clamped * es
        end

        dx = Enzyme.autodiff(Enzyme.Reverse, composite,
            Enzyme.Active, Enzyme.Active(285.0))
        @test isfinite(dx[1][1])

        # Compare against a numerical finite difference of the SAME Float64 primal.
        # (smooth_* are type-based: the Float64 forward pass Enzyme differentiates
        # is exact, so the reference must be exact too — a Dual/ForwardDiff
        # reference would smooth and legitimately differ near a clamp boundary.)
        h = 1.0e-4
        fd = (composite(285.0 + h) - composite(285.0 - h)) / (2h)
        rel_err = abs(dx[1][1] - fd) / max(abs(fd), 1e-20)
        println("  composite: Enzyme=$(dx[1][1]), numFD=$(fd), err=$(round(rel_err*100, digits=2))%")
        @test rel_err < 0.05
    end

end
