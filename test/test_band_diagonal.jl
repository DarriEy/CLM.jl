using Test
using LinearAlgebra

@testset "Band Diagonal Solver" begin

    @testset "band_diagonal_solve! tridiagonal (nband=3)" begin
        # 4×4 diagonally-dominant tridiagonal system
        n = 4
        nband = 3
        ncols = 2

        # A = [4 -1  0  0;
        #     -1  4 -1  0;
        #      0 -1  4 -1;
        #      0  0 -1  4]
        A = diagm(0 => fill(4.0, n), -1 => fill(-1.0, n-1), 1 => fill(-1.0, n-1))
        x_exact = [1.0, 2.0, 3.0, 4.0]
        rhs_exact = A * x_exact

        # CLM band storage: [ncols × nband × nlevs]
        # Fortran row-oriented convention:
        #   band 1 = A[j, j+1] (superdiag: coupling to next level)
        #   band 2 = A[j, j]   (diagonal)
        #   band 3 = A[j, j-1] (subdiag: coupling to previous level)
        b_matrix = zeros(ncols, nband, n)
        r = zeros(ncols, n)
        u = zeros(ncols, n)

        for col in 1:ncols
            for j in 1:n
                b_matrix[col, 2, j] = A[j, j]              # diagonal
                if j < n
                    b_matrix[col, 1, j] = A[j, j+1]        # superdiag (coupling down)
                end
                if j > 1
                    b_matrix[col, 3, j] = A[j, j-1]        # subdiag (coupling up)
                end
            end
            r[col, :] .= rhs_exact
        end

        jtop = ones(Int, ncols)
        jbot = fill(n, ncols)

        CLM.band_diagonal_solve!(u, b_matrix, r, jtop, jbot, nband, ncols)

        for col in 1:ncols
            @test u[col, :] ≈ x_exact atol=1e-12
        end
    end

    @testset "band_diagonal_solve! pentadiagonal (nband=5)" begin
        # 5×5 diagonally-dominant pentadiagonal system
        n = 5
        nband = 5

        A = diagm(0 => fill(10.0, n),
                  -1 => fill(-2.0, n-1), 1 => fill(-2.0, n-1),
                  -2 => fill(-1.0, n-2), 2 => fill(-1.0, n-2))
        x_exact = collect(1.0:n)
        rhs_exact = A * x_exact

        ncols = 1
        b_matrix = zeros(ncols, nband, n)
        r = zeros(ncols, n)
        u = zeros(ncols, n)

        # Fortran row-oriented convention for nband=5:
        #   band 1 = A[j, j+2] (2nd superdiag)
        #   band 2 = A[j, j+1] (1st superdiag)
        #   band 3 = A[j, j]   (diagonal)
        #   band 4 = A[j, j-1] (1st subdiag)
        #   band 5 = A[j, j-2] (2nd subdiag)
        for j in 1:n
            b_matrix[1, 3, j] = A[j, j]
            if j < n;     b_matrix[1, 2, j] = A[j, j+1]; end
            if j < n - 1; b_matrix[1, 1, j] = A[j, j+2]; end
            if j > 1;     b_matrix[1, 4, j] = A[j, j-1]; end
            if j > 2;     b_matrix[1, 5, j] = A[j, j-2]; end
        end
        r[1, :] .= rhs_exact

        jtop = [1]
        jbot = [n]

        CLM.band_diagonal_solve!(u, b_matrix, r, jtop, jbot, nband, ncols)

        @test u[1, :] ≈ x_exact atol=1e-10
    end

    @testset "band_diagonal_column! tridiagonal" begin
        n = 4
        nband = 3

        A = diagm(0 => fill(4.0, n), -1 => fill(-1.0, n-1), 1 => fill(-1.0, n-1))
        x_exact = [1.0, 2.0, 3.0, 4.0]
        rhs_exact = A * x_exact

        b_matrix = zeros(nband, n)
        for j in 1:n
            b_matrix[2, j] = A[j, j]
            if j < n; b_matrix[1, j] = A[j, j+1]; end
            if j > 1; b_matrix[3, j] = A[j, j-1]; end
        end

        r = copy(rhs_exact)
        u = zeros(n)

        CLM.band_diagonal_column!(u, b_matrix, r, 1, n, nband)

        @test u ≈ x_exact atol=1e-12
    end

    @testset "band_diagonal_column! pentadiagonal" begin
        n = 6
        nband = 5

        A = diagm(0 => fill(10.0, n),
                  -1 => fill(-2.0, n-1), 1 => fill(-2.0, n-1),
                  -2 => fill(-1.0, n-2), 2 => fill(-1.0, n-2))
        x_exact = collect(1.0:n)
        rhs_exact = A * x_exact

        b_matrix = zeros(nband, n)
        for j in 1:n
            b_matrix[3, j] = A[j, j]
            if j < n;     b_matrix[2, j] = A[j, j+1]; end
            if j < n - 1; b_matrix[1, j] = A[j, j+2]; end
            if j > 1;     b_matrix[4, j] = A[j, j-1]; end
            if j > 2;     b_matrix[5, j] = A[j, j-2]; end
        end

        r = copy(rhs_exact)
        u = zeros(n)

        CLM.band_diagonal_column!(u, b_matrix, r, 1, n, nband)

        @test u ≈ x_exact atol=1e-10
    end

    @testset "agreement with dense solver" begin
        # Random diagonally-dominant tridiagonal system
        n = 10
        nband = 3

        A = zeros(n, n)
        for i in 1:n
            A[i, i] = 10.0 + rand()
            if i > 1; A[i, i-1] = -rand(); end
            if i < n; A[i, i+1] = -rand(); end
        end

        x_dense = randn(n)
        rhs = A * x_dense

        # Row-oriented band storage:
        # band 1 = A[j, j+1] (superdiag), band 2 = A[j, j], band 3 = A[j, j-1] (subdiag)
        b_matrix = zeros(nband, n)
        for j in 1:n
            b_matrix[2, j] = A[j, j]
            if j < n; b_matrix[1, j] = A[j, j+1]; end
            if j > 1; b_matrix[3, j] = A[j, j-1]; end
        end

        r = copy(rhs)
        u = zeros(n)

        CLM.band_diagonal_column!(u, b_matrix, r, 1, n, nband)

        x_ref = A \ rhs
        @test u ≈ x_ref atol=1e-10
    end

    @testset "partial column range (jtop > 1)" begin
        nlevs = 8
        nband = 3
        jtop = 3
        jbot = 6
        n = jbot - jtop + 1

        A = diagm(0 => fill(5.0, n), -1 => fill(-1.0, n-1), 1 => fill(-1.0, n-1))
        x_exact = [2.0, 4.0, 6.0, 8.0]
        rhs_exact = A * x_exact

        b_matrix = zeros(nband, nlevs)
        for j in 1:n
            jj = jtop + j - 1
            b_matrix[2, jj] = A[j, j]
            if j < n; b_matrix[1, jj] = A[j, j+1]; end
            if j > 1; b_matrix[3, jj] = A[j, j-1]; end
        end

        r = zeros(nlevs)
        r[jtop:jbot] .= rhs_exact
        u = zeros(nlevs)

        CLM.band_diagonal_column!(u, b_matrix, r, jtop, jbot, nband)

        @test u[jtop:jbot] ≈ x_exact atol=1e-12
        @test all(u[1:jtop-1] .== 0.0)
        @test all(u[jbot+1:end] .== 0.0)
    end
end
