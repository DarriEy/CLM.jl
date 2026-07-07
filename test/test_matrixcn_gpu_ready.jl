@testset "matrix-CN sparse ops GPU-readiness (Float32 via ref, KA kernels)" begin
    # The kernelized sparse ops (spmm_ax!/spmm_ak!/spmp_b_acc!/set_value_*) + the device-
    # aware init (ref/FT) are the Metal foundation. Metal isn't loadable in this env, so
    # we prove the SAME path on the CPU KA backend at Metal's Float32 precision + via the
    # `similar(ref)` device-alloc route: build the workspaces in Float32 (ref=Float32
    # array, FT=Float32) and confirm the ops match the Float64 result within Float32 tol.
    # A wrong Float32 literal / scalar-arg / alloc would show up here (or fail to run).
    begu, endu = 1, 3; nunit = endu - begu + 1; SM = 4
    filt = [1, 2, 3]; num = 3

    # Build a small A·X advance in a given precision and return the advanced X matrix.
    function run_advance(FT, ref)
        A = CLM.SparseMatrixType(); CLM.init_sm!(A, SM, begu, endu; ref=ref, FT=FT)
        # 3 off-diagonal transfer entries (receiver,doner): 2→1, 3→2, 4→3.
        RI = [2, 3, 4]; CI = [1, 2, 3]; NE = 3
        A.RI[1:NE] .= RI; A.CI[1:NE] .= CI; A.NE = NE
        for u in 1:nunit, k in 1:NE
            A.M[u, k] = FT(0.1 * k + 0.05 * u)     # transfer fractions
        end
        # K diagonal turnover, applied to A (spmm_ak!): scale entry k by K[doner].
        K = CLM.DiagMatrixType(); CLM.init_dm!(K, SM, begu, endu; ref=ref, FT=FT)
        Ksrc = zeros(FT, nunit, SM)
        for u in 1:nunit, i in 1:SM; Ksrc[u, i] = FT(0.5 + 0.1 * i); end
        CLM.set_value_dm!(K, begu, endu, num, filt, Ksrc)
        CLM.spmm_ak!(A, num, filt, K)              # A ← A·K
        # X input vector, then advance X ← X + A·X (spmm_ax!).
        X = CLM.VectorType(); CLM.init_v!(X, SM, begu, endu; ref=ref, FT=FT)
        Xsrc = zeros(FT, nunit, SM)
        for u in 1:nunit, i in 1:SM; Xsrc[u, i] = FT(100.0 + 10.0 * i + u); end
        CLM.set_value_v!(X, begu, endu, num, filt, Xsrc)
        CLM.spmm_ax!(X, num, filt, A)
        return Array{Float64}(X.V)                 # collect to host Float64 for compare
    end

    X64 = run_advance(Float64, nothing)                     # host Float64 default path
    X32 = run_advance(Float32, zeros(Float32, 1, 1))        # Float32 via ref (Metal precision)

    # Same values to Float32 round-off (relative), proving the ops + alloc run correctly
    # at device precision.
    @test eltype(X32) == Float64        # (collected)
    @test isapprox(X32, X64; rtol = 1e-5)
    @test maximum(abs.(X32 .- X64) ./ abs.(X64)) < 1e-5

    # The Float32 path genuinely allocated Float32 workspaces (not silently Float64):
    let A = CLM.SparseMatrixType()
        CLM.init_sm!(A, SM, begu, endu; ref = zeros(Float32, 1, 1), FT = Float32)
        @test eltype(A.M) == Float32
    end
    let V = CLM.VectorType()
        CLM.init_v!(V, SM, begu, endu; ref = zeros(Float32, 1, 1), FT = Float32)
        @test eltype(V.V) == Float32
    end
    # Default (no ref) stays Float64 host.
    let V = CLM.VectorType(); CLM.init_v!(V, SM, begu, endu)
        @test eltype(V.V) == Float64 && V.V isa Matrix
    end
end
