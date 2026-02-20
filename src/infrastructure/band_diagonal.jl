# ==========================================================================
# Ported from: src/biogeophys/BandDiagonalMod.F90 (224 lines)
# Band diagonal matrix solver for column-based vertical diffusion
# Uses LAPACK dgbsv via Julia's LinearAlgebra
# ==========================================================================

"""
    band_diagonal_solve!(u, b_matrix, r, jtop, jbot, nband, ncols)

Solve band diagonal systems for multiple columns.
Uses LAPACK's dgbsv (general band matrix solver with partial pivoting).

# Arguments
- `u::Matrix{Float64}`:       Solution [ncols × nlevs] (OUTPUT)
- `b_matrix::Array{Float64,3}`: Band matrix [ncols × nband × nlevs]
- `r::Matrix{Float64}`:       Right-hand side [ncols × nlevs]
- `jtop::Vector{Int}`:        Top active level per column
- `jbot::Vector{Int}`:        Bottom active level per column
- `nband::Int`:               Band width (must be odd; e.g., 5 for pentadiagonal)
- `ncols::Int`:               Number of columns

For each column, the band matrix b_matrix[col, :, :] contains the
coefficients arranged so that band index (nband+1)/2 is the diagonal.

Ported from BandDiagonal() in BandDiagonalMod.F90.
"""
function band_diagonal_solve!(u::Matrix{Float64}, b_matrix::Array{Float64,3},
                              r::Matrix{Float64}, jtop::Vector{Int},
                              jbot::Vector{Int}, nband::Int, ncols::Int)
    kl = div(nband - 1, 2)  # number of sub-diagonals
    ku = kl                   # number of super-diagonals (symmetric bandwidth)

    for col in 1:ncols
        jt = jtop[col]
        jb = jbot[col]
        n = jb - jt + 1  # system size

        if n <= 0
            continue
        end

        # Build LAPACK band storage format
        # LAPACK dgbsv expects: AB[2*kl+ku+1, n] with the band stored
        # in rows kl+1 to 2*kl+ku+1 (1-indexed)
        m = 2 * kl + ku + 1
        ab = zeros(m, n)

        for j in 1:n
            for band_idx in 1:nband
                # Map from CLM band index to LAPACK storage
                # CLM: band 1 = lowest sub-diagonal, band (nband+1)/2 = diagonal
                # LAPACK: row kl+ku+1+i-j for element (i,j)
                row_offset = band_idx - (kl + 1)  # -kl to +ku
                i = j + row_offset  # row in original matrix

                if 1 <= i <= n
                    # LAPACK row index: kl + ku + 1 + (i-1) - (j-1) = kl + ku + 1 + i - j
                    ab_row = kl + ku + 1 + i - j
                    ab[ab_row, j] = b_matrix[col, band_idx, jt + j - 1]
                end
            end
        end

        # Right-hand side (copied to allow LAPACK to overwrite)
        rhs = r[col, jt:jb]

        # Solve using LAPACK dgbsv
        ipiv = zeros(Int64, n)

        # Use Julia's LAPACK interface
        try
            LinearAlgebra.LAPACK.gbsv!(kl, ku, n, ab, ipiv, reshape(rhs, n, 1))
            u[col, jt:jb] .= rhs
        catch e
            @warn "Band diagonal solve failed for column $col: $e"
            u[col, jt:jb] .= 0.0
        end
    end

    nothing
end

"""
    band_diagonal_column!(u, b_matrix, r, jtop, jbot, nband)

Solve a single-column band diagonal system.

# Arguments
- `u::AbstractVector{Float64}`:  Solution (OUTPUT, length nlevs)
- `b_matrix::AbstractMatrix{Float64}`: Band matrix [nband × nlevs]
- `r::AbstractVector{Float64}`:  Right-hand side [nlevs]
- `jtop::Int`:                   Top active level
- `jbot::Int`:                   Bottom active level
- `nband::Int`:                  Band width
"""
function band_diagonal_column!(u::AbstractVector{Float64}, b_matrix::AbstractMatrix{Float64},
                               r::AbstractVector{Float64}, jtop::Int, jbot::Int, nband::Int)
    kl = div(nband - 1, 2)
    ku = kl
    n = jbot - jtop + 1

    if n <= 0
        return nothing
    end

    m = 2 * kl + ku + 1
    ab = zeros(m, n)

    for j in 1:n
        for band_idx in 1:nband
            row_offset = band_idx - (kl + 1)
            i = j + row_offset
            if 1 <= i <= n
                ab_row = kl + ku + 1 + i - j
                ab[ab_row, j] = b_matrix[band_idx, jtop + j - 1]
            end
        end
    end

    rhs = copy(r[jtop:jbot])
    ipiv = zeros(Int64, n)

    LinearAlgebra.LAPACK.gbsv!(kl, ku, n, ab, ipiv, reshape(rhs, n, 1))
    u[jtop:jbot] .= rhs

    nothing
end
