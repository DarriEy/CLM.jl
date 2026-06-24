using LinearAlgebra

@testset "SparseMatrixMultiplyMod" begin

    # Helper: build a COO sparse matrix for a single unit (begu=endu=1) from a
    # dense reference, keeping the Fortran ordering contract (row index varies
    # fastest, i.e. entries sorted by the column-major linear index).
    function build_sm_from_dense(D::Matrix{Float64})
        n = size(D, 1)
        @assert size(D, 2) == n
        RI = Int[]; CI = Int[]; vals = Float64[]
        for j in 1:n, i in 1:n     # column-major: i (row) fastest within column j
            if D[i, j] != 0.0
                push!(RI, i); push!(CI, j); push!(vals, D[i, j])
            end
        end
        ne = length(vals)
        sm = CLM.SparseMatrixType()
        CLM.init_sm!(sm, n, 1, 1)
        Mvals = reshape(copy(vals), 1, ne)
        CLM.set_value_sm!(sm, 1, 1, 1, [1], Mvals, RI, CI, ne)
        return sm
    end

    @testset "sentinels and construction" begin
        @test CLM.SMM_EMPTY_INT == -9999
        @test CLM.SMM_EMPTY_REAL == -9999.0
        sm = CLM.SparseMatrixType()
        @test !CLM.is_alloc_sm(sm)
        @test !CLM.is_values_set_sm(sm)
        CLM.init_sm!(sm, 4, 1, 3)
        @test CLM.is_alloc_sm(sm)
        @test sm.SM == 4
        @test sm.begu == 1 && sm.endu == 3
        @test size(sm.M) == (3, 16)
        @test !CLM.is_values_set_sm(sm)   # values not set yet
    end

    @testset "set values + dense round-trip" begin
        D = [1.0 0.0 2.0;
             0.0 3.0 0.0;
             4.0 0.0 5.0]
        sm = build_sm_from_dense(D)
        @test CLM.is_values_set_sm(sm)
        @test sm.NE == 5
        @test CLM.sm_to_dense(sm, 1) ≈ D
    end

    @testset "diagonal / identity / scaling (SetValueA_diag, SPMM_AK)" begin
        sm = CLM.SparseMatrixType()
        CLM.init_sm!(sm, 3, 1, 1)
        CLM.set_value_a_diag!(sm, 1, [1], 1.0)       # identity
        @test sm.NE == 3
        @test CLM.sm_to_dense(sm, 1) ≈ Matrix{Float64}(LinearAlgebra.I, 3, 3)

        # Scale a general sparse matrix by a diagonal matrix: A·K.
        D = [2.0 0.0 1.0;
             0.0 4.0 0.0;
             3.0 0.0 5.0]
        A = build_sm_from_dense(D)
        K = CLM.DiagMatrixType()
        CLM.init_dm!(K, 3, 1, 1)
        kdiag = [10.0 100.0 1000.0]   # 1×3
        CLM.set_value_dm!(K, 1, 1, 1, [1], kdiag)
        CLM.spmm_ak!(A, 1, [1], K)
        @test CLM.sm_to_dense(A, 1) ≈ D * LinearAlgebra.diagm([10.0, 100.0, 1000.0])
    end

    @testset "matrix-vector multiply (SPMM_AX): X = X + A·X" begin
        D = [0.0 0.1 0.0;
             0.2 0.0 0.3;
             0.0 0.0 0.4]
        A = build_sm_from_dense(D)
        x0 = [5.0, 7.0, 11.0]
        vec = CLM.VectorType()
        CLM.init_v!(vec, 3, 1, 1)
        CLM.set_value_v!(vec, 1, 1, 1, [1], reshape(copy(x0), 1, 3))
        CLM.spmm_ax!(vec, 1, [1], A)
        ref = x0 .+ D * x0           # dense reference: (I + A)·x0
        @test vec.V[1, :] ≈ ref
    end

    @testset "scaler vector set" begin
        vec = CLM.VectorType()
        CLM.init_v!(vec, 4, 1, 1)
        CLM.set_value_v_scaler!(vec, 1, [1], 2.5)
        @test all(vec.V[1, :] .== 2.5)
    end

    @testset "addition (SPMP_AB): AB = A + B vs dense" begin
        DA = [1.0 0.0 2.0;
              0.0 3.0 0.0;
              4.0 0.0 0.0]
        DB = [0.5 0.0 0.0;
              0.0 0.0 1.0;
              4.0 0.0 6.0]
        A = build_sm_from_dense(DA)
        B = build_sm_from_dense(DB)
        AB = CLM.SparseMatrixType()
        CLM.init_sm!(AB, 3, 1, 1)
        list_ready, ne = CLM.spmp_ab!(AB, 1, [1], A, B, false)
        @test CLM.sm_to_dense(AB, 1) ≈ DA + DB
        # entries should be sorted by column-major linear index
        lin = [(AB.CI[k]-1)*AB.SM + AB.RI[k] for k in 1:AB.NE]
        @test issorted(lin)
    end

    @testset "addition with memorized structure (list_ready path)" begin
        DA = [1.0 0.0 2.0;
              0.0 3.0 0.0;
              4.0 0.0 0.0]
        DB = [0.5 0.0 0.0;
              0.0 0.0 1.0;
              4.0 0.0 6.0]
        A = build_sm_from_dense(DA)
        B = build_sm_from_dense(DB)
        AB = CLM.SparseMatrixType()
        CLM.init_sm!(AB, 3, 1, 1)
        nmax = A.NE + B.NE
        list_A = zeros(Int, nmax); list_B = zeros(Int, nmax)
        RI_AB = zeros(Int, nmax); CI_AB = zeros(Int, nmax)
        list_ready, NE_AB = CLM.spmp_ab!(AB, 1, [1], A, B, false;
            list_A=list_A, list_B=list_B, NE_AB=0, RI_AB=RI_AB, CI_AB=CI_AB)
        @test list_ready
        dense_first = CLM.sm_to_dense(AB, 1)
        # Re-run on the memorized path with new values; structure must be reused.
        DA2 = DA .* 2.0
        A2 = build_sm_from_dense(DA2)
        AB2 = CLM.SparseMatrixType()
        CLM.init_sm!(AB2, 3, 1, 1)
        CLM.spmp_ab!(AB2, 1, [1], A2, B, true;
            list_A=list_A, list_B=list_B, NE_AB=NE_AB, RI_AB=RI_AB, CI_AB=CI_AB)
        @test CLM.sm_to_dense(AB2, 1) ≈ DA2 + DB
        @test dense_first ≈ DA + DB
    end

    @testset "accumulation (SPMP_B_ACC): B += A (same structure)" begin
        D = [1.0 0.0 2.0;
             0.0 3.0 0.0;
             4.0 0.0 5.0]
        A = build_sm_from_dense(D)
        B = build_sm_from_dense(D)           # same structure, same values
        CLM.spmp_b_acc!(B, 1, [1], A)
        @test CLM.sm_to_dense(B, 1) ≈ 2.0 .* D
    end

    @testset "three-way addition (SPMP_ABC): ABC = A + B + C vs dense" begin
        DA = [1.0 0.0 0.0;
              0.0 2.0 0.0;
              0.0 0.0 0.0]
        DB = [0.0 5.0 0.0;
              0.0 2.0 0.0;
              7.0 0.0 0.0]
        DC = [1.0 0.0 3.0;
              0.0 0.0 0.0;
              0.0 0.0 9.0]
        A = build_sm_from_dense(DA)
        B = build_sm_from_dense(DB)
        C = build_sm_from_dense(DC)
        ABC = CLM.SparseMatrixType()
        CLM.init_sm!(ABC, 3, 1, 1)
        CLM.spmp_abc!(ABC, 1, [1], A, B, C, false)
        @test CLM.sm_to_dense(ABC, 1) ≈ DA + DB + DC
    end

    @testset "SetValueA: A = offdiag - I on the diagonal" begin
        # off-diagonal transfer entries (strictly off-diagonal; SetValueA forces
        # the diagonal to exactly -1, so any i==j entries are ignored)
        offdiag = [0.0 0.0 0.0;
                   0.0 0.0 0.2;
                   0.1 0.4 0.0]
        # collect non-zero off-diagonal entries column-major
        RI = Int[]; CI = Int[]; vals = Float64[]
        for j in 1:3, i in 1:3
            (i != j) || continue
            if offdiag[i, j] != 0.0
                push!(RI, i); push!(CI, j); push!(vals, offdiag[i, j])
            end
        end
        ne_non = length(vals)
        A = CLM.SparseMatrixType()
        CLM.init_sm!(A, 3, 1, 1)
        Mvals = reshape(copy(vals), 1, ne_non)
        init_ready = CLM.set_value_a!(A, 1, 1, 1, [1], Mvals, RI, CI, ne_non, false)
        @test init_ready
        expected = offdiag - Matrix{Float64}(LinearAlgebra.I, 3, 3)
        @test CLM.sm_to_dense(A, 1) ≈ expected
    end

    @testset "pool-advance round-trip: X1 = X0 + A·X0 + B·u vs dense" begin
        # Transfer/turnover sparse matrix A and input sparse matrix B.
        DA = [-0.5  0.2  0.0;
               0.1 -0.3  0.0;
               0.0  0.1 -0.2]
        # B (input matrix) is square (SM=3) with the input mapping on the diagonal.
        DBsq = [1.0 0.0 0.0;
                0.0 1.0 0.0;
                0.0 0.0 0.0]
        x0 = [10.0, 20.0, 30.0]
        u  = [2.0, 3.0, 0.0]   # external input vector

        A = build_sm_from_dense(DA)
        B = build_sm_from_dense(DBsq)

        # X step 1: X = X + A·X   (SPMM_AX)
        vec = CLM.VectorType()
        CLM.init_v!(vec, 3, 1, 1)
        CLM.set_value_v!(vec, 1, 1, 1, [1], reshape(copy(x0), 1, 3))
        CLM.spmm_ax!(vec, 1, [1], A)
        # X step 2: add B·u  (do the B·u via a fresh vector then accumulate)
        Bu = DBsq * u
        x_after_A = x0 .+ DA * x0
        @test vec.V[1, :] ≈ x_after_A
        # full advance
        x1 = vec.V[1, :] .+ Bu
        ref = x0 .+ DA * x0 .+ DBsq * u
        @test x1 ≈ ref
    end

    @testset "batched over multiple units (begu:endu)" begin
        # Two units sharing one sparsity structure, different values.
        n = 3
        RI = [1, 2, 1, 3]; CI = [1, 1, 2, 3]   # column-major sorted
        ne = 4
        sm = CLM.SparseMatrixType()
        CLM.init_sm!(sm, n, 1, 2)   # units 1 and 2
        Mvals = [ 1.0  2.0  3.0  4.0;     # unit 1 row
                 10.0 20.0 30.0 40.0]     # unit 2 row
        CLM.set_value_sm!(sm, 1, 2, 2, [1, 2], Mvals, RI, CI, ne)
        D1 = zeros(3, 3); D1[1,1]=1; D1[2,1]=2; D1[1,2]=3; D1[3,3]=4
        D2 = zeros(3, 3); D2[1,1]=10; D2[2,1]=20; D2[1,2]=30; D2[3,3]=40
        @test CLM.sm_to_dense(sm, 1) ≈ D1
        @test CLM.sm_to_dense(sm, 2) ≈ D2

        # Matrix-vector on both units at once.
        vec = CLM.VectorType()
        CLM.init_v!(vec, n, 1, 2)
        Vvals = [1.0 1.0 1.0;
                 2.0 2.0 2.0]
        CLM.set_value_v!(vec, 1, 2, 2, [1, 2], Vvals)
        CLM.spmm_ax!(vec, 2, [1, 2], sm)
        @test vec.V[1, :] ≈ [1.0, 1.0, 1.0] .+ D1 * [1.0, 1.0, 1.0]
        @test vec.V[2, :] ≈ [2.0, 2.0, 2.0] .+ D2 * [2.0, 2.0, 2.0]
    end

    @testset "equivalence / copy helpers" begin
        D = [1.0 0.0 2.0; 0.0 3.0 0.0; 4.0 0.0 5.0]
        A = build_sm_from_dense(D)
        B = build_sm_from_dense(D)
        @test CLM.is_equiv_idx_sm(A, B)
        # copy values into a fresh allocated matrix
        Ccopy = CLM.SparseMatrixType()
        CLM.init_sm!(Ccopy, 3, 1, 1)
        CLM.set_value_copy_sm!(Ccopy, 1, [1], A)
        @test CLM.sm_to_dense(Ccopy, 1) ≈ D
        @test CLM.is_equiv_idx_sm(Ccopy, A)
        # release resets allocation
        CLM.release_sm!(A)
        @test !CLM.is_alloc_sm(A)
    end

end
