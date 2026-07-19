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

!!! note "Reverse-mode AD"
    Enzyme does NOT differentiate this body: a custom adjoint rule in
    `src/infrastructure/enzyme_rules.jl` supplies the exact linear-solve adjoint
    (λ = A⁻ᵀ ū). That is a correctness requirement, not an optimization — on
    Julia 1.10 Enzyme leaves the shadow of the `cp`/`dp` workspace allocated
    below UNZEROED, so the generic reverse pass returns a wrong-but-finite
    gradient and, on a few percent of calls, a NaN/Inf. If you change this
    function's signature or its band convention, update the rule (and its
    finite-difference test in `test/test_enzyme_smoke.jl`) with it.
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
# Batched Thomas-algorithm kernel: one thread per column, each solving its own
# tridiagonal system. The forward-sweep workspace (cp, dp) is per-column scratch
# [ncols × nlevs] so threads never share state (the serial version reused 1D cp/dp,
# which is not thread-safe). Backend-agnostic via the output `u`.
@kernel function _tridiag_multi_kernel!(u, @Const(a), @Const(b), @Const(c), @Const(r),
                                        @Const(jtop), @Const(mask), cp, dp, nlevs::Int)
    col = @index(Global)
    @inbounds if mask[col]
        j = jtop[col]

        # Forward sweep. `tiny` is the singularity guard at the working precision so the
        # kernel carries no Float64 on a Float32-only backend (Metal); on Float64 it is
        # exactly 1e-30 as before.
        bet = b[col, j]
        tiny = oftype(bet, 1e-30)
        if abs(bet) < tiny; bet = tiny; end
        # The sweeps carry cp/dp/u from the previous level in a LOCAL, and store to
        # the arrays only for the back-substitution's benefit. Do NOT re-read
        # cp[col,jj-1] / dp[col,jj-1] / u[col,jj+1] from memory here: that
        # store-then-load recurrence is miscompiled by LLVM at -O2 on the KA CPU
        # backend, silently returning wrong values (or NaN) for nlevs >= 5 and
        # corrupting the jtop > 1 path. Keeping it in a register is both correct
        # and cheaper. See test_tridiagonal.jl's deep-column regression tests.
        cp_prev = c[col, j] / bet
        dp_prev = r[col, j] / bet
        cp[col, j] = cp_prev
        dp[col, j] = dp_prev

        for jj in (j+1):nlevs
            denom = b[col, jj] - a[col, jj] * cp_prev
            if abs(denom) < tiny; denom = copysign(tiny, denom == 0 ? one(denom) : denom); end
            cp_prev = c[col, jj] / denom
            dp_prev = (r[col, jj] - a[col, jj] * dp_prev) / denom
            cp[col, jj] = cp_prev
            dp[col, jj] = dp_prev
        end

        # Back substitution (dp_prev still holds dp[col, nlevs])
        u_next = dp_prev
        u[col, nlevs] = u_next
        for jj in (nlevs-1):-1:j
            u_next = dp[col, jj] - cp[col, jj] * u_next
            u[col, jj] = u_next
        end
    end
end

function tridiagonal_multi!(u::AbstractMatrix{<:Real}, a::AbstractMatrix{<:Real}, b::AbstractMatrix{<:Real},
                            c::AbstractMatrix{<:Real}, r::AbstractMatrix{<:Real},
                            jtop::AbstractVector{<:Integer}, mask::AbstractVector{Bool}, ncols::Int, nlevs::Int)
    # Per-column scratch [ncols × nlevs] (similar(u, …) keeps the array/backend type,
    # so this runs as a CPU loop on Arrays and as a GPU kernel on device arrays).
    T = promote_type(eltype(a), eltype(b), eltype(c), eltype(r))
    cp = similar(u, T, ncols, nlevs)
    dp = similar(u, T, ncols, nlevs)

    _launch!(_tridiag_multi_kernel!, u, a, b, c, r, jtop, mask, cp, dp, nlevs;
             ndrange = ncols)

    nothing
end
