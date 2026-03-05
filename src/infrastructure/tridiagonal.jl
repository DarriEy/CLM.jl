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
- `u::AbstractVector{<:Real}`: Solution vector [nlevs] (OUTPUT)
- `a::AbstractVector{<:Real}`: Sub-diagonal [nlevs]
- `b::AbstractVector{<:Real}`: Diagonal [nlevs]
- `c::AbstractVector{<:Real}`: Super-diagonal [nlevs]
- `r::AbstractVector{<:Real}`: Right-hand side [nlevs]
- `jtop::Int`:                  Top active level
- `nlevs::Int`:                 Total number of levels

This is the reference implementation matching Fortran exactly.
"""
function tridiagonal_solve!(u::AbstractVector{<:Real}, a::AbstractVector{<:Real},
                            b::AbstractVector{<:Real}, c::AbstractVector{<:Real},
                            r::AbstractVector{<:Real}, jtop::Int, nlevs::Int)
    T = promote_type(eltype(a), eltype(b), eltype(c), eltype(r))
    cp = Vector{T}(undef, nlevs)
    dp = Vector{T}(undef, nlevs)

    j = jtop

    @inbounds begin
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
    end

    nothing
end

"""
    tridiagonal_multi!(u, a, b, c, r, jtop, mask, ncols, nlevs)

Solve tridiagonal systems for multiple columns simultaneously.
Columns where `mask[col] == false` are skipped.

# Arguments
- `u::AbstractMatrix{<:Real}`:     Solution [ncols × nlevs] (OUTPUT)
- `a::AbstractMatrix{<:Real}`:     Sub-diagonal [ncols × nlevs]
- `b::AbstractMatrix{<:Real}`:     Diagonal [ncols × nlevs]
- `c::AbstractMatrix{<:Real}`:     Super-diagonal [ncols × nlevs]
- `r::AbstractMatrix{<:Real}`:     Right-hand side [ncols × nlevs]
- `jtop::Vector{Int}`:      Top active level per column [ncols]
- `mask::AbstractVector{Bool}`:    Active column mask [ncols]
- `ncols::Int`:              Number of columns
- `nlevs::Int`:              Number of levels
"""
function tridiagonal_multi!(u::AbstractMatrix{<:Real}, a::AbstractMatrix{<:Real}, b::AbstractMatrix{<:Real},
                            c::AbstractMatrix{<:Real}, r::AbstractMatrix{<:Real},
                            jtop::Vector{Int}, mask::AbstractVector{Bool}, ncols::Int, nlevs::Int)
    # Workspace allocated once, reused across all columns
    T = promote_type(eltype(a), eltype(b), eltype(c), eltype(r))
    cp = Vector{T}(undef, nlevs)
    dp = Vector{T}(undef, nlevs)

    for col in 1:ncols
        mask[col] || continue

        j = jtop[col]

        @inbounds begin
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
    end

    nothing
end
