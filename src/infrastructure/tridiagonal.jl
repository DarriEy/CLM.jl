# ==========================================================================
# Ported from: src/biogeophys/TridiagonalMod.F90 (94 lines)
# Tridiagonal matrix solver for column-based vertical diffusion
# ==========================================================================

"""
    tridiagonal_solve!(u, a, b, c, r, jtop, nlevs)

Solve a single-column tridiagonal system using the Thomas algorithm.

The system solved for levels jtop:nlevs is:
    a[j]*u[j-1] + b[j]*u[j] + c[j]*u[j+1] = r[j]

# Arguments
- `u::AbstractVector{Float64}`: Solution vector [nlevs] (OUTPUT)
- `a::AbstractVector{Float64}`: Sub-diagonal [nlevs]
- `b::AbstractVector{Float64}`: Diagonal [nlevs]
- `c::AbstractVector{Float64}`: Super-diagonal [nlevs]
- `r::AbstractVector{Float64}`: Right-hand side [nlevs]
- `jtop::Int`:                  Top active level
- `nlevs::Int`:                 Total number of levels

This is the reference implementation matching Fortran exactly.
"""
function tridiagonal_solve!(u::AbstractVector{Float64}, a::AbstractVector{Float64},
                            b::AbstractVector{Float64}, c::AbstractVector{Float64},
                            r::AbstractVector{Float64}, jtop::Int, nlevs::Int)
    # Working arrays for modified coefficients
    cp = zeros(nlevs)
    dp = zeros(nlevs)

    j = jtop

    # Forward sweep: eliminate sub-diagonal
    cp[j] = c[j] / b[j]
    dp[j] = r[j] / b[j]

    for jj in (j+1):nlevs
        denom = b[jj] - a[jj] * cp[jj-1]
        cp[jj] = c[jj] / denom
        dp[jj] = (r[jj] - a[jj] * dp[jj-1]) / denom
    end

    # Back substitution
    u[nlevs] = dp[nlevs]
    for jj in (nlevs-1):-1:j
        u[jj] = dp[jj] - cp[jj] * u[jj+1]
    end

    nothing
end

"""
    tridiagonal_multi!(u, a, b, c, r, jtop, mask, ncols, nlevs)

Solve tridiagonal systems for multiple columns simultaneously.
Columns where `mask[col] == false` are skipped.

# Arguments
- `u::Matrix{Float64}`:     Solution [ncols × nlevs] (OUTPUT)
- `a::Matrix{Float64}`:     Sub-diagonal [ncols × nlevs]
- `b::Matrix{Float64}`:     Diagonal [ncols × nlevs]
- `c::Matrix{Float64}`:     Super-diagonal [ncols × nlevs]
- `r::Matrix{Float64}`:     Right-hand side [ncols × nlevs]
- `jtop::Vector{Int}`:      Top active level per column [ncols]
- `mask::BitVector`:         Active column mask [ncols]
- `ncols::Int`:              Number of columns
- `nlevs::Int`:              Number of levels
"""
function tridiagonal_multi!(u::Matrix{Float64}, a::Matrix{Float64}, b::Matrix{Float64},
                            c::Matrix{Float64}, r::Matrix{Float64},
                            jtop::Vector{Int}, mask::BitVector, ncols::Int, nlevs::Int)
    # Per-column working arrays
    cp = zeros(nlevs)
    dp = zeros(nlevs)

    for col in 1:ncols
        mask[col] || continue

        j = jtop[col]

        # Forward sweep
        cp[j] = c[col, j] / b[col, j]
        dp[j] = r[col, j] / b[col, j]

        for jj in (j+1):nlevs
            denom = b[col, jj] - a[col, jj] * cp[jj-1]
            cp[jj] = c[col, jj] / denom
            dp[jj] = (r[col, jj] - a[col, jj] * dp[jj-1]) / denom
        end

        # Back substitution
        u[col, nlevs] = dp[nlevs]
        for jj in (nlevs-1):-1:j
            u[col, jj] = dp[jj] - cp[jj] * u[col, jj+1]
        end
    end

    nothing
end
