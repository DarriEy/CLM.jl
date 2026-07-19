using Random
using LinearAlgebra

@testset "Tridiagonal Solver" begin

    @testset "Simple 3x3 system" begin
        nlevs = 3
        a = [0.0, -1.0, -1.0]
        b = [2.0, 2.0, 2.0]
        c = [-1.0, -1.0, 0.0]
        r = [1.0, 0.0, 1.0]
        u = zeros(nlevs)

        CLM.tridiagonal_solve!(u, a, b, c, r, 1, nlevs)
        @test u ≈ [1.0, 1.0, 1.0] atol=1e-12
    end

    @testset "Identity system" begin
        nlevs = 5
        a = zeros(nlevs)
        b = ones(nlevs)
        c = zeros(nlevs)
        r = [1.0, 2.0, 3.0, 4.0, 5.0]
        u = zeros(nlevs)

        CLM.tridiagonal_solve!(u, a, b, c, r, 1, nlevs)
        @test u ≈ r atol=1e-14
    end

    @testset "Partial system (jtop > 1)" begin
        nlevs = 5
        a = zeros(nlevs); b = ones(nlevs); c = zeros(nlevs); r = zeros(nlevs)

        b[3] = 2.0; b[4] = 2.0; b[5] = 2.0
        a[4] = -1.0; a[5] = -1.0
        c[3] = -1.0; c[4] = -1.0
        r[3] = 1.0; r[4] = 0.0; r[5] = 1.0
        u = zeros(nlevs)

        CLM.tridiagonal_solve!(u, a, b, c, r, 3, nlevs)

        @test u[3] ≈ 1.0 atol=1e-12
        @test u[4] ≈ 1.0 atol=1e-12
        @test u[5] ≈ 1.0 atol=1e-12
    end

    @testset "Residual check (Au = r)" begin
        nlevs = 10
        coeff = 1.0e-6 * 3600.0 / (0.1^2)

        a = fill(-coeff, nlevs)
        b = fill(1.0 + 2.0 * coeff, nlevs)
        c = fill(-coeff, nlevs)
        a[1] = 0.0; c[nlevs] = 0.0
        r = [300.0 - 2.0 * j for j in 1:nlevs]
        u = zeros(nlevs)

        CLM.tridiagonal_solve!(u, a, b, c, r, 1, nlevs)

        residual = zeros(nlevs)
        residual[1] = b[1] * u[1] + c[1] * u[2] - r[1]
        for j in 2:(nlevs-1)
            residual[j] = a[j] * u[j-1] + b[j] * u[j] + c[j] * u[j+1] - r[j]
        end
        residual[nlevs] = a[nlevs] * u[nlevs-1] + b[nlevs] * u[nlevs] - r[nlevs]

        @test all(abs.(residual) .< 1e-10)
    end

    @testset "Multi-column solver" begin
        ncols = 3; nlevs = 4

        a = zeros(ncols, nlevs)
        b = ones(ncols, nlevs) .* 2.0
        c = zeros(ncols, nlevs)
        r = zeros(ncols, nlevs)
        u = zeros(ncols, nlevs)
        jtop = ones(Int, ncols)
        mask = BitVector([true, true, true])

        for col in 1:ncols
            for j in 2:nlevs; a[col, j] = -1.0; end
            for j in 1:(nlevs-1); c[col, j] = -1.0; end
            r[col, 1] = 1.0; r[col, nlevs] = 1.0
        end

        CLM.tridiagonal_multi!(u, a, b, c, r, jtop, mask, ncols, nlevs)

        # Verify residual for each column
        for col in 1:ncols
            for j in 1:nlevs
                av = j > 1 ? a[col,j] * u[col,j-1] : 0.0
                cv = j < nlevs ? c[col,j] * u[col,j+1] : 0.0
                res = av + b[col,j] * u[col,j] + cv - r[col,j]
                @test abs(res) < 1e-10
            end
        end
    end

    @testset "Multi-column mask skips inactive" begin
        ncols = 3; nlevs = 3

        a = zeros(ncols, nlevs)
        b = ones(ncols, nlevs) .* 2.0
        c = zeros(ncols, nlevs)
        r = ones(ncols, nlevs)
        u = zeros(ncols, nlevs)
        jtop = ones(Int, ncols)
        mask = BitVector([true, false, true])

        CLM.tridiagonal_multi!(u, a, b, c, r, jtop, mask, ncols, nlevs)

        @test all(u[2, :] .== 0.0)  # skipped
        @test all(u[1, :] .> 0.0)   # solved
        @test all(u[3, :] .> 0.0)   # solved
    end

    # Deep columns vs an INDEPENDENT dense solve. The shallow (nlevs == 3) cases
    # above cannot catch a miscompiled forward sweep: the -O2 KA-CPU bug that
    # motivated the register-carried recurrence only appears from nlevs >= 5, and
    # real CLM runs this at nlevsoi+1 (~21) levels. Checking against LinearAlgebra's
    # Tridiagonal rather than a residual also catches a self-consistently wrong answer.
    @testset "Multi-column deep (nlevs >= 5) vs dense solve" begin
        rng = MersenneTwister(2)
        for (ncols, nlevs, jt) in ((48, 20, 1), (48, 5, 1), (7, 25, 1), (48, 20, 4), (1, 21, 1))
            a = -0.5 .* ones(ncols, nlevs)
            c = -0.5 .* ones(ncols, nlevs)
            b = 2 .+ 0.5 .* rand(rng, ncols, nlevs)
            r = rand(rng, ncols, nlevs)
            a[:, 1] .= 0.0
            c[:, nlevs] .= 0.0
            jtop = fill(jt, ncols)
            mask = trues(ncols)

            u = zeros(ncols, nlevs)
            CLM.tridiagonal_multi!(u, a, b, c, r, jtop, mask, ncols, nlevs)

            @test all(isfinite, u)
            for col in 1:ncols
                M = Tridiagonal(a[col, jt+1:nlevs], b[col, jt:nlevs], c[col, jt:nlevs-1])
                @test u[col, jt:nlevs] ≈ M \ r[col, jt:nlevs] atol=1e-10
            end
        end
    end

end
