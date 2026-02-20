@testset "Tridiagonal Solver" begin

    @testset "Simple 3x3 system" begin
        # Solve:  2x1 - x2          = 1
        #        -x1 + 2x2 - x3     = 0
        #              -x2 + 2x3     = 1
        # Solution: x = [1, 1, 1]

        nlevs = 3
        a = [0.0, -1.0, -1.0]  # sub-diagonal
        b = [2.0, 2.0, 2.0]    # diagonal
        c = [-1.0, -1.0, 0.0]  # super-diagonal
        r = [1.0, 0.0, 1.0]    # RHS
        u = zeros(nlevs)

        CLM.tridiagonal_solve!(u, a, b, c, r, 1, nlevs)

        @test u ≈ [1.0, 1.0, 1.0] atol=1e-12
    end

    @testset "Identity system" begin
        # b = 1 everywhere, a = c = 0 → u = r
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
        # Solve only levels 3-5 of a 5-level system
        nlevs = 5
        a = zeros(nlevs)
        b = ones(nlevs)
        c = zeros(nlevs)
        r = zeros(nlevs)

        # Set up a simple system at levels 3-5
        b[3] = 2.0; b[4] = 2.0; b[5] = 2.0
        a[4] = -1.0; a[5] = -1.0
        c[3] = -1.0; c[4] = -1.0
        r[3] = 1.0; r[4] = 0.0; r[5] = 1.0
        u = zeros(nlevs)

        CLM.tridiagonal_solve!(u, a, b, c, r, 3, nlevs)

        # Solution for the 3x3 subsystem should be [1, 1, 1]
        @test u[3] ≈ 1.0 atol=1e-12
        @test u[4] ≈ 1.0 atol=1e-12
        @test u[5] ≈ 1.0 atol=1e-12
    end

    @testset "Heat diffusion analogy" begin
        # Typical CLM use: implicit diffusion
        # b[j] = 1 + 2*D*dt/dz^2
        # a[j] = c[j] = -D*dt/dz^2

        nlevs = 10
        D = 1.0e-6  # thermal diffusivity
        dt = 3600.0  # 1 hour
        dz = 0.1     # 10 cm layers
        coeff = D * dt / (dz * dz)

        a = fill(-coeff, nlevs)
        b = fill(1.0 + 2.0 * coeff, nlevs)
        c = fill(-coeff, nlevs)
        a[1] = 0.0
        c[nlevs] = 0.0

        # Initial temperature profile: warm at surface, cold at depth
        r = [300.0 - 2.0 * j for j in 1:nlevs]
        u = zeros(nlevs)

        CLM.tridiagonal_solve!(u, a, b, c, r, 1, nlevs)

        # Solution should exist and be finite
        @test all(isfinite.(u))

        # Verify Au = r by reconstructing the matrix-vector product
        # This tests the solver correctness directly
        residual = zeros(nlevs)
        residual[1] = b[1] * u[1] + c[1] * u[2] - r[1]
        for j in 2:(nlevs-1)
            residual[j] = a[j] * u[j-1] + b[j] * u[j] + c[j] * u[j+1] - r[j]
        end
        residual[nlevs] = a[nlevs] * u[nlevs-1] + b[nlevs] * u[nlevs] - r[nlevs]

        @test all(abs.(residual) .< 1e-10)
    end

end
