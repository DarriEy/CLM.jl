@testset "cn_soil_matrix! at the real 7-pool BGC cascade topology (matrix == sequential)" begin
    # The unit test (test_cn_soil_matrix) validates the solver on a 3-pool toy cascade.
    # This exercises the ACTUAL CLM BGC cascade built by init_decomp_cascade_bgc! —
    # 7 pools, 10 transitions with MULTIPLE transitions per donor (4→{5,6}, 5→{4,6},
    # 7→{2,3}) — so a sparse-index-assembly bug in the matrix would diverge from a
    # straightforward dense cascade reference. tri_ma_vr=0 (transport is a separate,
    # already-tested term) so the cascade equivalence is exact.
    nc = 1; nlevdecomp = 2; ndecomp_pools = 7; ndct = 10; dt = 1800.0
    ndp_vr = ndecomp_pools * nlevdecomp

    # --- real BGC cascade topology ---
    bgc_state = CLM.DecompBGCState(); cc = CLM.DecompCascadeConData()
    pb = CLM.DecompBGCParams(bgc_initial_Cstocks = fill(200.0, ndecomp_pools))
    cnp = CLM.CNSharedParamsData()
    CLM.init_decomp_cascade_bgc!(bgc_state, cc, pb, cnp; cellsand = fill(50.0, nc, 5),
        bounds = 1:nc, nlevdecomp = nlevdecomp, ndecomp_pools_max = ndecomp_pools,
        ndecomp_cascade_transitions_max = ndct, use_fates = false)
    # fixed C:N (all pools) so the N cascade reference is unambiguous (floating-C:N is a
    # separate path); keep the real donor/receiver topology + is_cwd + initial_cn.
    cc.floating_cn_ratio_decomp_pools = falses(ndecomp_pools)
    donor = cc.cascade_donor_pool; recv = cc.cascade_receiver_pool
    cn = cc.initial_cn_ratio
    nouttrans = count(==(0), recv)

    # per-pool turnover rate (per step) + per-transition respired + path fractions that
    # sum to 1 across each donor's outgoing transitions.
    decomp_k = [0.05, 0.04, 0.03, 0.02, 0.015, 0.008, 0.01]
    rf = [0.3, 0.3, 0.3, 0.4, 0.4, 0.5, 0.6, 0.4, 0.2, 0.2]
    pf = ones(ndct)
    ntr = zeros(Int, ndecomp_pools); for k in 1:ndct; ntr[donor[k]] += 1; end
    seen = zeros(Int, ndecomp_pools)
    for k in 1:ndct
        d = donor[k]; seen[d] += 1
        pf[k] = ntr[d] == 1 ? 1.0 : (seen[d] == 1 ? 0.6 : 0.4)   # split multi-transition donors
    end

    C0 = zeros(nc, nlevdecomp, ndecomp_pools); N0 = similar(C0)
    for i in 1:ndecomp_pools, j in 1:nlevdecomp
        C0[1, j, i] = 100.0 + 10.0 * i + 3.0 * j; N0[1, j, i] = C0[1, j, i] / cn[i]
    end

    # --- dense sequential reference (per transition; multi-per-donor handled by pf) ---
    seqC = copy(C0); seqN = copy(N0)
    for k in 1:ndct
        d = donor[k]; r = recv[k]
        for j in 1:nlevdecomp
            # Ksoil (= decomp_k here) is already the per-step turnover fraction.
            closs = decomp_k[d] * pf[k] * C0[1, j, d]
            nloss = decomp_k[d] * pf[k] * N0[1, j, d]
            seqC[1, j, d] -= closs; seqN[1, j, d] -= nloss
            if r != 0
                seqC[1, j, r] += (1 - rf[k]) * closs
                seqN[1, j, r] += (1 - rf[k]) * nloss * (cn[d] / cn[r])
            end
        end
    end

    # --- matrix advance ---
    matC = copy(C0); matN = copy(N0)
    Ksoil = zeros(nc, ndp_vr)
    for i in 1:ndecomp_pools, j in 1:nlevdecomp; Ksoil[1, j + (i - 1) * nlevdecomp] = decomp_k[i]; end
    rf3 = zeros(nc, nlevdecomp, ndct); pf3 = zeros(nc, nlevdecomp, ndct)
    for j in 1:nlevdecomp, k in 1:ndct; rf3[1, j, k] = rf[k]; pf3[1, j, k] = pf[k]; end
    Ntri = (nlevdecomp * 3 - 2) * (ndecomp_pools - 1)
    CLM.init_soil_transfer!(cc; ndecomp_pools = ndecomp_pools, nlevdecomp = nlevdecomp,
        ndecomp_cascade_transitions = ndct, ndecomp_cascade_outtransitions = nouttrans)
    ms = CLM.CNSoilMatrixState()
    CLM.cn_soil_matrix!(ms, cc; decomp_cpools_vr = matC, decomp_npools_vr = matN,
        Ksoil = Ksoil, tri_ma_vr = zeros(nc, Ntri),
        matrix_Cinput = zeros(nc, ndp_vr), matrix_Ninput = zeros(nc, ndp_vr),
        rf_decomp_cascade = rf3, pathfrac_decomp_cascade = pf3,
        mask_soilc = [true], begc = 1, endc = 1, nlevdecomp = nlevdecomp,
        ndecomp_pools = ndecomp_pools, ndecomp_cascade_transitions = ndct,
        ndecomp_cascade_outtransitions = nouttrans)

    @test matC ≈ seqC atol = 1e-9 rtol = 1e-11
    @test matN ≈ seqN atol = 1e-9 rtol = 1e-11
    @test maximum(abs.(matC .- seqC)) < 1e-9
    @test maximum(abs.(matN .- seqN)) < 1e-9
    # sanity: pool 7 (cwd) is a PURE donor (never a receiver) so it loses exactly its
    # turnover; the som pools 4/5 are both donors AND receivers (multi-transition), so
    # the solve genuinely moved mass along the real cascade (not the identity).
    for j in 1:nlevdecomp
        @test seqC[1, j, 7] ≈ C0[1, j, 7] * (1 - decomp_k[7]) atol = 1e-9
    end
    @test any(abs.(matC .- C0) .> 1.0)   # the cascade actually advanced the pools
end
