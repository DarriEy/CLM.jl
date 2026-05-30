# ==========================================================================
# Ported from: src/biogeophys/BandDiagonalMod.F90 (224 lines)
# Band diagonal matrix solver for column-based vertical diffusion
# Uses LAPACK dgbsv via Julia's LinearAlgebra
# ==========================================================================

"""
    _band_solve_julia!(u_slice, ab, rhs, kl, ku, n, FT)

Pure-Julia band diagonal solver without try/catch.
Converts band storage to dense and uses backslash.
Enzyme-safe: no LAPACK, no exception handling.

Returns true if solve succeeded, false if singular.
"""
function _band_solve_julia!(u_slice, ab, rhs, kl, ku, n, ::Type{FT}) where FT
    A_dense = zeros(FT, n, n)
    for j in 1:n
        for di in -kl:ku
            i = j + di
            if 1 <= i <= n
                A_dense[i, j] = ab[kl + ku + 1 + i - j, j]
            end
        end
    end
    # Use LU factorization (no try/catch for Enzyme compatibility)
    F = lu(A_dense, check=false)
    if issuccess(F)
        sol = F \ rhs
        u_slice .= sol
        return true
    end
    return false
end

"""
    band_solve!(u_slice, ab, rhs, kl, ku, n)

Solve a single banded linear system `A x = rhs` and write `x` into `u_slice`
(length `n`). `ab` is LAPACK band storage of size `(2kl+ku+1, n)`.

- `Float64`/`Float32`: solve with LAPACK `gbtrf!`/`gbtrs!` (factorizing a *copy*
  so the caller's `ab` is preserved — the Enzyme reverse rule reads `A` from it).
- Other element types (e.g. `ForwardDiff.Dual`): pure-Julia `_band_solve_julia!`.

A custom `EnzymeRules` reverse rule for this function (registered when Enzyme is
loaded) differentiates the Float64 LAPACK path via the linear-solve adjoint
(`Aᵀ λ = x̄`), so Enzyme never needs to trace through LAPACK or its `try/catch`.
"""
function band_solve!(u_slice, ab::AbstractMatrix{T}, rhs::AbstractVector,
                     kl::Int, ku::Int, n::Int) where {T}
    if T <: AbstractFloat
        ab_work = copy(ab)
        rhs_work = reshape(copy(rhs), n, 1)
        ipiv = zeros(Int64, n)
        try
            (ab_lu, ipiv) = LinearAlgebra.LAPACK.gbtrf!(kl, ku, n, ab_work)
            LinearAlgebra.LAPACK.gbtrs!('N', kl, ku, n, ab_lu, ipiv, rhs_work)
            u_slice .= vec(rhs_work)
        catch e
            e isa LinearAlgebra.LAPACKException || rethrow()
        end
    else
        _band_solve_julia!(u_slice, ab, rhs, kl, ku, n, T)
    end
    return nothing
end

# ==========================================================================
# Batched banded solve (GPU): one thread per column, each solving its own
# banded system by dense Gaussian elimination WITH partial pivoting (matching
# LAPACK gbsv / lu() to round-off). Per-column scratch is allocation-free inside
# the kernel. Used by the soil-temperature pentadiagonal solve.
#
# Band -> dense mapping (derived from set_matrix!'s band storage + the densify in
# _band_solve_julia!):  A[i,j] = bmatrix[c, (i-j)+kl+1, jt+i-1]  on the band
# (max(1,j-kl) <= i <= min(n,j+ku)), where jt = jtop[c]+nlevsno_off+1 and
# n = jbot[c]-jtop[c]+1. rhs[i] = rvector[c, jt+i-1]; solution -> tvector[c, jt+i-1].
#
# This deliberately replicates the pure-Julia (Dual) path's "densify + pivoted LU"
# rather than the host LAPACK band factorization — both are partial-pivoting GE and
# agree to round-off. The host band_solve! (and its Enzyme reverse rule) is left
# intact for the non-kernel / Enzyme path. ForwardDiff flows through: the pivot
# selection compares abs(.) on the value component, derivatives propagate through
# the eliminations/solves.
# ==========================================================================
@kernel function _batched_band_solve_kernel!(tvector, @Const(bmatrix), @Const(rvector),
                                             @Const(jtop), @Const(jbot), @Const(mask),
                                             A, x, kl::Int, ku::Int, nlevsno_off::Int)
    c = @index(Global)
    @inbounds if mask[c]
        jt = jtop[c] + nlevsno_off + 1
        n = jbot[c] - jtop[c] + 1
        if n > 0
            # Build dense A[c,:,:] and rhs x[c,:] for this column.
            for j in 1:n
                for i in 1:n
                    A[c, i, j] = zero(eltype(A))
                end
            end
            for j in 1:n
                ilo = max(1, j - kl)
                ihi = min(n, j + ku)
                for i in ilo:ihi
                    A[c, i, j] = bmatrix[c, (i - j) + kl + 1, jt + i - 1]
                end
            end
            for i in 1:n
                x[c, i] = rvector[c, jt + i - 1]
            end

            # Gaussian elimination with partial pivoting (in place).
            for k in 1:(n - 1)
                piv = k
                maxv = abs(A[c, k, k])
                for i in (k + 1):n
                    v = abs(A[c, i, k])
                    if v > maxv
                        maxv = v
                        piv = i
                    end
                end
                if piv != k
                    for j in k:n
                        t = A[c, k, j]; A[c, k, j] = A[c, piv, j]; A[c, piv, j] = t
                    end
                    tx = x[c, k]; x[c, k] = x[c, piv]; x[c, piv] = tx
                end
                akk = A[c, k, k]
                for i in (k + 1):n
                    f = A[c, i, k] / akk
                    for j in k:n
                        A[c, i, j] = A[c, i, j] - f * A[c, k, j]
                    end
                    x[c, i] = x[c, i] - f * x[c, k]
                end
            end

            # Back substitution (solution overwrites x[c,:]).
            for i in n:-1:1
                s = x[c, i]
                for j in (i + 1):n
                    s = s - A[c, i, j] * x[c, j]
                end
                x[c, i] = s / A[c, i, i]
            end

            for i in 1:n
                tvector[c, jt + i - 1] = x[c, i]
            end
        end
    end
end

"""
    batched_band_solve!(tvector, bmatrix, rvector, jtop, jbot, mask, kl, ku, nlevsno_off)

Solve one banded system per active column (pentadiagonal soil temperature), batched
as a KernelAbstractions kernel — one thread per column, dense pivoted GE. Backend-
agnostic via the output `tvector`. See `_batched_band_solve_kernel!` for the layout.
"""
function batched_band_solve!(tvector::AbstractMatrix{T}, bmatrix, rvector,
                             jtop::AbstractVector{<:Integer}, jbot::AbstractVector{<:Integer}, mask,
                             kl::Int, ku::Int, nlevsno_off::Int) where {T}
    nc = size(tvector, 1)
    nmax = size(tvector, 2)
    A = similar(tvector, T, nc, nmax, nmax)
    x = similar(tvector, T, nc, nmax)
    _launch!(_batched_band_solve_kernel!, tvector, bmatrix, rvector, jtop, jbot, mask,
             A, x, kl, ku, nlevsno_off; ndrange = nc)
    return nothing
end

"""
    band_diagonal_solve!(u, b_matrix, r, jtop, jbot, nband, ncols)

Solve band diagonal systems for multiple columns.
Uses LAPACK's dgbsv (general band matrix solver with partial pivoting).

# Arguments
- `u::AbstractMatrix{<:Real}`:       Solution [ncols × nlevs] (OUTPUT)
- `b_matrix::AbstractArray{<:Real,3}`: Band matrix [ncols × nband × nlevs]
- `r::AbstractMatrix{<:Real}`:       Right-hand side [ncols × nlevs]
- `jtop::Vector{Int}`:        Top active level per column
- `jbot::Vector{Int}`:        Bottom active level per column
- `nband::Int`:               Band width (must be odd; e.g., 5 for pentadiagonal)
- `ncols::Int`:               Number of columns

For each column, the band matrix b_matrix[col, :, :] contains the
coefficients arranged so that band index (nband+1)/2 is the diagonal.

Ported from BandDiagonal() in BandDiagonalMod.F90.
"""
function band_diagonal_solve!(u::AbstractMatrix{<:Real}, b_matrix::AbstractArray{<:Real,3},
                              r::AbstractMatrix{<:Real}, jtop::Vector{Int},
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
        FT = eltype(b_matrix)
        m = 2 * kl + ku + 1
        ab = zeros(FT, m, n)

        for j in 1:n
            for band_idx in 1:nband
                # Map from CLM band index to LAPACK storage
                # CLM stores band data row-oriented: bmatrix[band, j] is band element of row j.
                # Fortran BandDiagonal shifts source index for off-diagonals:
                #   diagonal (band 3): source = jtop+j-1
                #   1st sub (band 4):  source = jtop+j   (shifted +1)
                #   1st super (band 2): source = jtop+j-2 (shifted -1)
                row_offset = band_idx - (kl + 1)  # -kl to +ku
                i = j + row_offset  # row in original matrix

                if 1 <= i <= n
                    ab_row = kl + ku + 1 + i - j
                    src_idx = jt + j - 1 + row_offset
                    ab[ab_row, j] = b_matrix[col, band_idx, src_idx]
                end
            end
        end

        # Right-hand side (copied to allow LAPACK to overwrite)
        rhs = r[col, jt:jb]

        # Solve the banded system
        if FT <: AbstractFloat
            # Use LAPACK gbtrf!/gbtrs! for Float64/Float32
            ipiv = zeros(Int64, n)
            rhs_mat = reshape(rhs, n, 1)
            try
                (ab, ipiv) = LinearAlgebra.LAPACK.gbtrf!(kl, ku, n, ab)
                LinearAlgebra.LAPACK.gbtrs!('N', kl, ku, n, ab, ipiv, rhs_mat)
                u[col, jt:jb] .= rhs
            catch e
                e isa LinearAlgebra.LAPACKException || rethrow()
            end
        else
            # Pure-Julia fallback for Dual/non-LAPACK types (Enzyme-safe)
            _band_solve_julia!(view(u, col, jt:jb), ab, rhs, kl, ku, n, FT)
        end
    end

    nothing
end

"""
    band_diagonal_column!(u, b_matrix, r, jtop, jbot, nband)

Solve a single-column band diagonal system.

# Arguments
- `u::AbstractVector{<:Real}`:  Solution (OUTPUT, length nlevs)
- `b_matrix::AbstractMatrix{<:Real}`: Band matrix [nband × nlevs]
- `r::AbstractVector{<:Real}`:  Right-hand side [nlevs]
- `jtop::Int`:                   Top active level
- `jbot::Int`:                   Bottom active level
- `nband::Int`:                  Band width
"""
function band_diagonal_column!(u::AbstractVector{<:Real}, b_matrix::AbstractMatrix{<:Real},
                               r::AbstractVector{<:Real}, jtop::Int, jbot::Int, nband::Int)
    kl = div(nband - 1, 2)
    ku = kl
    n = jbot - jtop + 1

    if n <= 0
        return nothing
    end

    FT = eltype(b_matrix)
    m = 2 * kl + ku + 1
    ab = zeros(FT, m, n)

    for j in 1:n
        for band_idx in 1:nband
            row_offset = band_idx - (kl + 1)
            i = j + row_offset
            if 1 <= i <= n
                ab_row = kl + ku + 1 + i - j
                src_idx = jtop + j - 1 + row_offset
                ab[ab_row, j] = b_matrix[band_idx, src_idx]
            end
        end
    end

    rhs = copy(r[jtop:jbot])

    if FT <: AbstractFloat
        # Use LAPACK for Float64/Float32
        ipiv = zeros(Int64, n)
        rhs_mat = reshape(rhs, n, 1)
        try
            (ab, ipiv) = LinearAlgebra.LAPACK.gbtrf!(kl, ku, n, ab)
            LinearAlgebra.LAPACK.gbtrs!('N', kl, ku, n, ab, ipiv, rhs_mat)
            u[jtop:jbot] .= rhs
        catch e
            e isa LinearAlgebra.LAPACKException || rethrow()
        end
    else
        # Pure-Julia fallback for Dual/non-LAPACK types
        A_dense = zeros(FT, n, n)
        for j in 1:n
            for band_idx in 1:nband
                row_offset = band_idx - (kl + 1)
                i = j + row_offset
                if 1 <= i <= n
                    A_dense[i, j] = ab[kl + ku + 1 + i - j, j]
                end
            end
        end
        try
            sol = A_dense \ rhs
            u[jtop:jbot] .= sol
        catch
        end
    end

    nothing
end
