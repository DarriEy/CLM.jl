@testset "litter_vert_transp tri_ma_vr (vertical-transport matrix build)" begin
    # Phase-3 kernel work: in use_soil_matrixcn mode litter_vert_transp! builds the
    # sparse vertical-transport matrix tri_ma_vr (Fortran SoilBiogeochemLittVertTransp
    # :401-429) instead of applying the implicit Thomas solve. The matrix method then
    # advances transport EXPLICITLY (X = X + V·X, matching Fortran CNSoilMatrix), so it
    # TRACKS — not bit-matches — the implicit sequential solve; the gap is O(V) and →0
    # as the diffusivity (hence V) shrinks. A wrong tri_ma_vr (sign/scale/structure)
    # would diverge grossly instead.
    nc = 1; nlevdecomp = 4; ndecomp_pools = 3   # pools 1,2 non-cwd; 3 = cwd
    ndp_vr = ndecomp_pools * nlevdecomp
    scalez = 0.025; nlevgrnd = nlevdecomp + 5
    zsoi = [scalez * (exp(0.5 * (j - 0.5)) - 1.0) for j in 1:nlevgrnd]
    zisoi = zeros(nlevgrnd + 1)
    for j in 1:nlevgrnd
        zisoi[j + 1] = j < nlevgrnd ? 0.5 * (zsoi[j] + zsoi[j + 1]) : zsoi[j] + 0.5 * (zsoi[j] - zisoi[j])
    end
    dzd = [zisoi[j + 1] - zisoi[j] for j in 1:nlevdecomp]

    function mkfix()
        cs = CLM.SoilBiogeochemCarbonStateData(); cs.decomp_cpools_vr_col = zeros(nc, nlevdecomp, ndecomp_pools)
        for s in 1:ndecomp_pools, j in 1:nlevdecomp
            cs.decomp_cpools_vr_col[1, j, s] = 100.0 * exp(-0.5 * j) * (1 + 0.1 * s)
        end
        cf = CLM.SoilBiogeochemCarbonFluxData(); cf.decomp_cpools_sourcesink_col = zeros(nc, nlevdecomp, ndecomp_pools)
        cf.decomp_cpools_transport_tendency_col = zeros(nc, nlevdecomp, ndecomp_pools)
        cf.tri_ma_vr = zeros(nc, (nlevdecomp * 3 - 2) * (ndecomp_pools - 1))
        ns = CLM.SoilBiogeochemNitrogenStateData(); ns.decomp_npools_vr_col = zeros(nc, nlevdecomp, ndecomp_pools)
        for s in 1:ndecomp_pools, j in 1:nlevdecomp
            ns.decomp_npools_vr_col[1, j, s] = cs.decomp_cpools_vr_col[1, j, s] / 15.0
        end
        nf = CLM.SoilBiogeochemNitrogenFluxData(); nf.decomp_npools_sourcesink_col = zeros(nc, nlevdecomp, ndecomp_pools)
        nf.decomp_npools_transport_tendency_col = zeros(nc, nlevdecomp, ndecomp_pools)
        st = CLM.SoilBiogeochemStateData(); st.som_adv_coef_col = zeros(nc, nlevdecomp + 1); st.som_diffus_coef_col = zeros(nc, nlevdecomp + 1)
        col = CLM.ColumnData(); col.nbedrock = fill(nlevdecomp, nc); col.gridcell = ones(Int, nc)
        grc = CLM.GridcellData(); grc.latdeg = fill(45.0, 1)
        al = CLM.ActiveLayerData(); al.altmax_col = fill(3.0, nc); al.altmax_lastyear_col = fill(3.0, nc)
        cc = CLM.DecompCascadeConData(); ic = falses(ndecomp_pools); ic[3] = true; cc.is_cwd = ic; cc.spinup_factor = ones(ndecomp_pools)
        p = CLM.LitterVertTranspParams()
        return (cs=cs, cf=cf, ns=ns, nf=nf, st=st, col=col, grc=grc, cc=cc, al=al, p=p)
    end
    runlvt(f; mx) = CLM.litter_vert_transp!(f.cs, f.cf, f.ns, f.nf, f.st, f.col, f.grc, f.cc, f.al, f.p;
        mask_bgc_soilc=trues(nc), bounds=1:nc, dtime=1800.0, nlevdecomp=nlevdecomp, ndecomp_pools=ndecomp_pools,
        zsoi_vals=zsoi, dzsoi_decomp_vals=dzd, zisoi_vals=zisoi, use_soil_matrixcn=mx)

    # apply the built tri_ma_vr through cn_soil_matrix! with NO cascade (Ksoil=0) and
    # NO input (Cinput=0): the pools then advance by the vertical transport term alone.
    cc2 = CLM.DecompCascadeConData(cascade_donor_pool=[1, 2, 3], cascade_receiver_pool=[2, 0, 1],
        floating_cn_ratio_decomp_pools=BitVector([false, false, false]),
        is_cwd=BitVector([false, false, true]), initial_cn_ratio=[20.0, 15.0, 90.0])
    CLM.init_soil_transfer!(cc2; ndecomp_pools=ndecomp_pools, nlevdecomp=nlevdecomp,
        ndecomp_cascade_transitions=3, ndecomp_cascade_outtransitions=1)

    function matrix_transport(diff)
        fm = mkfix(); fm.p.som_diffus = diff; fm.p.cryoturb_diffusion_k = 5 * diff; fm.p.max_altdepth_cryoturbation = 2.0
        runlvt(fm; mx=true)   # builds tri_ma_vr; state untouched
        ms = CLM.CNSoilMatrixState()
        matC = copy(fm.cs.decomp_cpools_vr_col); matN = copy(fm.ns.decomp_npools_vr_col)
        CLM.cn_soil_matrix!(ms, cc2; decomp_cpools_vr=matC, decomp_npools_vr=matN,
            Ksoil=zeros(nc, ndp_vr), tri_ma_vr=fm.cf.tri_ma_vr,
            matrix_Cinput=zeros(nc, ndp_vr), matrix_Ninput=zeros(nc, ndp_vr),
            rf_decomp_cascade=zeros(nc, nlevdecomp, 3), pathfrac_decomp_cascade=ones(nc, nlevdecomp, 3),
            mask_soilc=[true], begc=1, endc=1, nlevdecomp=nlevdecomp, ndecomp_pools=ndecomp_pools,
            ndecomp_cascade_transitions=3, ndecomp_cascade_outtransitions=1)
        return matC, fm.cf.tri_ma_vr
    end
    function seq_transport(diff)
        fs = mkfix(); fs.p.som_diffus = diff; fs.p.cryoturb_diffusion_k = 5 * diff; fs.p.max_altdepth_cryoturbation = 2.0
        runlvt(fs; mx=false)
        return fs.cs.decomp_cpools_vr_col
    end

    D = 2.0e-9
    matC, tri = matrix_transport(D)
    seqC = seq_transport(D)

    # (1) tri_ma_vr is finite, nonzero, and block-repeated per non-CWD pool.
    stride = nlevdecomp * 3 - 2
    @test all(isfinite, tri)
    @test any(!=(0.0), tri)
    @test tri[1, 1:stride] ≈ tri[1, (stride + 1):(2 * stride)]

    # (2) matrix vertical transport TRACKS the implicit sequential solve (structure /
    # sign / scale correct) — a bug diverges by orders of magnitude, not <1%.
    relerr(D) = begin
        m, _ = matrix_transport(D); s = seq_transport(D)
        maximum(abs.(m[1, :, 1:2] .- s[1, :, 1:2])) / maximum(abs.(s[1, :, 1:2]))
    end
    e = relerr(D)
    @test e < 0.01

    # (3) explicit↔implicit gap is O(V): shrinking the diffusivity 10× shrinks the
    # gap ~10× (converges to the sequential limit) — the discriminating correctness
    # signal for the operator.
    e10 = relerr(D / 10)
    @test e10 < e / 5              # ~10× smaller (allow slack)
    @test e10 < 0.002

    # (4) transport is mass-conserving to the explicit scheme's truncation.
    mass(X) = sum(X[1, j, i] * dzd[j] for j in 1:nlevdecomp, i in 1:2)
    m0 = mass(mkfix().cs.decomp_cpools_vr_col)
    @test isapprox(mass(matC), m0; rtol=1e-3)
end
