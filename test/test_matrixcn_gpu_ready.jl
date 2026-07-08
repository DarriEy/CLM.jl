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

@testset "veg-C solver runs end-to-end in Float32 (ref-threaded workspaces)" begin
    # The solver ref-threading makes cn_veg_matrix_solve_c! allocate its workspaces via
    # ref/FT. Byte-identity tests only exercise the Float64 default; here we confirm the
    # Float32 solver PATH actually executes (a stray Float64 literal / scalar-arg in the
    # solver would fail) and gives the same pools as Float64 to Float32 round-off. Uses a
    # topology-only setup (small transfer/turnover), so the advance genuinely moves pools.
    np = 3; nveg = CLM.NVEGPOOL_NATVEG; counts = CLM.veg_matrix_transfer_counts(false)
    ivt = [1, 2, 3]; woody = ones(Float64, 5); npcropmin = 15; dt = 1800.0

    function run_solve(FT, ref)
        cs = CLM.CNVegCarbonStateData(); CLM.cnveg_carbon_state_init!(cs, np, 1, 1; use_matrixcn=false, nrepr=1)
        # every veg-C pool (main + storage + xfer) must be finite — the solve loads all.
        for (k, f) in enumerate((:leafc_patch, :leafc_storage_patch, :leafc_xfer_patch,
                :frootc_patch, :frootc_storage_patch, :frootc_xfer_patch,
                :livestemc_patch, :livestemc_storage_patch, :livestemc_xfer_patch,
                :deadstemc_patch, :deadstemc_storage_patch, :deadstemc_xfer_patch,
                :livecrootc_patch, :livecrootc_storage_patch, :livecrootc_xfer_patch,
                :deadcrootc_patch, :deadcrootc_storage_patch, :deadcrootc_xfer_patch))
            getfield(cs, f) .= Float64[10.0 * k + p for p in 1:np]
        end
        cf = CLM.CNVegCarbonFluxData(); CLM.cnveg_carbon_flux_init!(cf, np, 1, 1; use_matrixcn=true)
        CLM.cn_veg_matrix_c_topology!(cf; use_crop=false, nvegcpool=nveg)
        # small transfers (rate) + a turnover on each pool so Aoned = transfer·dt/turnover.
        fill!(cf.matrix_phtransfer_patch, 1.0e-6); fill!(cf.matrix_phturnover_patch, 0.05)
        fill!(cf.matrix_gmtransfer_patch, 5.0e-7); fill!(cf.matrix_gmturnover_patch, 0.02)
        fill!(cf.matrix_fitransfer_patch, 0.0);    fill!(cf.matrix_fiturnover_patch, 0.0)
        fill!(cf.matrix_alloc_patch, 0.0);         fill!(cf.matrix_Cinput_patch, 0.0)
        CLM.cn_veg_matrix_solve_c!(cs, cf; mask_soilp=trues(np), bounds_patch=1:np,
            ivt=ivt, woody=woody, npcropmin=npcropmin, nvegcpool=nveg, counts=counts,
            dt=dt, num_actfirep=0, ref=ref, FT=FT)
        return Float64[cs.leafc_patch; cs.frootc_patch; cs.livestemc_patch;
                       cs.deadstemc_patch; cs.livecrootc_patch; cs.deadcrootc_patch]
    end

    p64 = run_solve(Float64, nothing)                    # default host Float64
    p32 = run_solve(Float32, zeros(Float32, 1, 1))       # Float32 workspaces via ref
    @test all(isfinite, p32)
    @test all(isfinite, p64)
    @test p64[1] != 11.0                                  # leafc[1] advanced from its init (10·1+1)
    @test isapprox(p32, p64; rtol = 1e-4)                # Float32 solver == Float64 to round-off
end

@testset "veg-C solve memoized-structure state == stateless solve" begin
    # A CNVegMatrixSolveState memoizes the sparse structure so calls after the first reuse
    # it via the kernelized memoized fills (the device solve path). On the CPU backend the
    # BUILD call (call 1) and the MEMOIZED call (call 2) must both be bit-identical to a
    # stateless solve on the same input.
    np = 3; nveg = CLM.NVEGPOOL_NATVEG; counts = CLM.veg_matrix_transfer_counts(false)
    ivt = [1, 2, 3]; woody = ones(Float64, 5); npcropmin = 15; dt = 1800.0
    mkcs() = (cs = CLM.CNVegCarbonStateData(); CLM.cnveg_carbon_state_init!(cs, np, 1, 1; use_matrixcn=false, nrepr=1);
        for (k, f) in enumerate((:leafc_patch, :leafc_storage_patch, :leafc_xfer_patch,
                :frootc_patch, :frootc_storage_patch, :frootc_xfer_patch,
                :livestemc_patch, :livestemc_storage_patch, :livestemc_xfer_patch,
                :deadstemc_patch, :deadstemc_storage_patch, :deadstemc_xfer_patch,
                :livecrootc_patch, :livecrootc_storage_patch, :livecrootc_xfer_patch,
                :deadcrootc_patch, :deadcrootc_storage_patch, :deadcrootc_xfer_patch))
            getfield(cs, f) .= Float64[10.0 * k + p for p in 1:np]; end; cs)
    cf = CLM.CNVegCarbonFluxData(); CLM.cnveg_carbon_flux_init!(cf, np, 1, 1; use_matrixcn=true)
    CLM.cn_veg_matrix_c_topology!(cf; use_crop=false, nvegcpool=nveg)
    fill!(cf.matrix_phtransfer_patch, 1.0e-6); fill!(cf.matrix_phturnover_patch, 0.05)
    fill!(cf.matrix_gmtransfer_patch, 5.0e-7); fill!(cf.matrix_gmturnover_patch, 0.02)
    fill!(cf.matrix_fitransfer_patch, 0.0);    fill!(cf.matrix_fiturnover_patch, 0.0)
    fill!(cf.matrix_alloc_patch, 0.0);         fill!(cf.matrix_Cinput_patch, 0.0)
    args = (; mask_soilp=trues(np), bounds_patch=1:np, ivt=ivt, woody=woody,
            npcropmin=npcropmin, nvegcpool=nveg, counts=counts, dt=dt, num_actfirep=0)
    pools(cs) = Float64[cs.leafc_patch; cs.frootc_patch; cs.livestemc_patch;
                        cs.deadstemc_patch; cs.livecrootc_patch; cs.deadcrootc_patch]
    st = CLM.CNVegMatrixSolveState()
    csref = mkcs(); CLM.cn_veg_matrix_solve_c!(csref, cf; args...);              gold = pools(csref)
    csb = mkcs();   CLM.cn_veg_matrix_solve_c!(csb, cf; args..., state=st)       # call 1: build
    @test pools(csb) == gold
    # ph + fi have in-veg transfers (nnon>0) so they memoize set_value_a!; gm is a pure
    # diagonal (nnon=0 → the kernelized set_value_a_diag!, rebuilt each call, no memo flag).
    @test st.init_ph && st.init_fi && st.ab_ready && st.NE_ab > 0
    csm = mkcs();   CLM.cn_veg_matrix_solve_c!(csm, cf; args..., state=st)       # call 2: memoized
    @test pools(csm) == gold
end

@testset "soil solver runs end-to-end in Float32 (ref-threaded workspaces)" begin
    # cn_soil_matrix! ref-threading — same as the veg smoke: confirm the Float32 soil
    # solve path executes and matches Float64 to round-off. 3-pool toy cascade
    # (1→2, 2→atm, 3→1), tri_ma_vr=0 so the cascade advance is exact.
    nc = 1; nlevdecomp = 2; ndecomp_pools = 3; ndct = 3; nout = 1
    ndp_vr = ndecomp_pools * nlevdecomp; cn = [20.0, 15.0, 90.0]
    cc = CLM.DecompCascadeConData(cascade_donor_pool=[1, 2, 3], cascade_receiver_pool=[2, 0, 1],
        floating_cn_ratio_decomp_pools=BitVector([false, false, false]),
        is_cwd=BitVector([false, false, true]), initial_cn_ratio=cn)
    CLM.init_soil_transfer!(cc; ndecomp_pools=ndecomp_pools, nlevdecomp=nlevdecomp,
        ndecomp_cascade_transitions=ndct, ndecomp_cascade_outtransitions=nout)
    rf = zeros(nc, nlevdecomp, ndct); rf[1, :, 1] .= 0.3; rf[1, :, 2] .= 1.0; rf[1, :, 3] .= 0.0
    pf = ones(nc, nlevdecomp, ndct)
    Ksoil = zeros(nc, ndp_vr)
    for i in 1:ndecomp_pools, j in 1:nlevdecomp; Ksoil[1, j + (i - 1) * nlevdecomp] = 0.02 * i; end
    Cin = zeros(nc, ndp_vr); Nin = zeros(nc, ndp_vr)
    for i in 1:ndecomp_pools, j in 1:nlevdecomp; Cin[1, j + (i - 1) * nlevdecomp] = 2.0 * i; Nin[1, j + (i - 1) * nlevdecomp] = 2.0 * i / cn[i]; end

    function run_soil(FT, ref)
        C = zeros(nc, nlevdecomp, ndecomp_pools); N = similar(C)
        for i in 1:ndecomp_pools, j in 1:nlevdecomp; C[1, j, i] = 100.0 * i + 10.0 * j; N[1, j, i] = C[1, j, i] / cn[i]; end
        ms = CLM.CNSoilMatrixState()
        CLM.cn_soil_matrix!(ms, cc; decomp_cpools_vr=C, decomp_npools_vr=N,
            Ksoil=copy(Ksoil), tri_ma_vr=zeros(nc, cc.Ntri_setup),
            matrix_Cinput=copy(Cin), matrix_Ninput=copy(Nin),
            rf_decomp_cascade=rf, pathfrac_decomp_cascade=pf,
            mask_soilc=[true], begc=1, endc=1, nlevdecomp=nlevdecomp, ndecomp_pools=ndecomp_pools,
            ndecomp_cascade_transitions=ndct, ndecomp_cascade_outtransitions=nout,
            ref=ref, FT=FT)
        return Float64[vec(C); vec(N)]
    end
    s64 = run_soil(Float64, nothing)
    s32 = run_soil(Float32, zeros(Float32, 1, 1))
    @test all(isfinite, s32)
    @test isapprox(s32, s64; rtol = 1e-4)
end
