@testset "CNSoilMatrix carbon isotopes (C13/C14 decomp advance)" begin
    # The isotope decomp pools advance with the SAME transfer operator (AKallsoilc) as
    # C — transfer/cascade fractions are isotope-independent (C14 radioactive decay is
    # applied separately by c14_decay!). With vertical transport off (tri_ma_vr=0) each
    # isotope advance must reproduce EXACTLY the sequential decomposition cascade plus
    # its own per-step B-input, on its own (isotopically distinct) pool values.
    ndecomp_pools = 3; nlevdecomp = 2
    ndecomp_cascade_transitions = 3; ndecomp_cascade_outtransitions = 1
    ndecomp_pools_vr = ndecomp_pools * nlevdecomp
    begc, endc = 1, 1; nunit = 1; dt = 1800.0
    cascade_donor_pool    = [1, 2, 3]
    cascade_receiver_pool = [2, 0, 1]          # 0 = atmosphere (terminal)
    initial_cn_ratio = [20.0, 15.0, 90.0]

    cc = CLM.DecompCascadeConData{Float64}(
        cascade_donor_pool = cascade_donor_pool,
        cascade_receiver_pool = cascade_receiver_pool,
        floating_cn_ratio_decomp_pools = falses(ndecomp_pools),
        is_cwd = BitVector([false, false, true]),
        initial_cn_ratio = initial_cn_ratio)
    CLM.init_soil_transfer!(cc; ndecomp_pools = ndecomp_pools, nlevdecomp = nlevdecomp,
        ndecomp_cascade_transitions = ndecomp_cascade_transitions,
        ndecomp_cascade_outtransitions = ndecomp_cascade_outtransitions)

    # C/N pools (drive the AKallsoilc build) + distinct isotope pools.
    C0  = zeros(nunit, nlevdecomp, ndecomp_pools); N0 = similar(C0)
    C13 = similar(C0); C14 = similar(C0)
    for i in 1:ndecomp_pools, j in 1:nlevdecomp
        C0[1, j, i]  = 100.0 * i + 10.0 * j
        N0[1, j, i]  = C0[1, j, i] / initial_cn_ratio[i]
        C13[1, j, i] = 0.011 * C0[1, j, i] + 0.3 * j      # isotopically distinct
        C14[1, j, i] = 1.0e-12 * C0[1, j, i] + 0.05 * i
    end

    rf = zeros(nunit, nlevdecomp, ndecomp_cascade_transitions)
    pathfrac = ones(nunit, nlevdecomp, ndecomp_cascade_transitions)
    rf[1, :, 1] .= 0.3; rf[1, :, 2] .= 1.0; rf[1, :, 3] .= 0.0
    decomp_k = zeros(nunit, nlevdecomp, ndecomp_pools)
    decomp_k[1, :, 1] .= 0.05; decomp_k[1, :, 2] .= 0.02; decomp_k[1, :, 3] .= 0.01
    Ksoil = zeros(nunit, ndecomp_pools_vr)
    for i in 1:ndecomp_pools, j in 1:nlevdecomp
        Ksoil[1, j + nlevdecomp * (i - 1)] = decomp_k[1, j, i]
    end
    matrix_Cinput = zeros(nunit, ndecomp_pools_vr); matrix_Ninput = zeros(nunit, ndecomp_pools_vr)
    C13in = zeros(nunit, ndecomp_pools_vr);         C14in = zeros(nunit, ndecomp_pools_vr)
    for i in 1:ndecomp_pools, j in 1:nlevdecomp
        vr = j + nlevdecomp * (i - 1)
        matrix_Cinput[1, vr] = 2.0 * i + 0.5 * j
        matrix_Ninput[1, vr] = matrix_Cinput[1, vr] / initial_cn_ratio[i]
        C13in[1, vr] = 0.011 * matrix_Cinput[1, vr] + 0.02 * j
        C14in[1, vr] = 1.0e-12 * matrix_Cinput[1, vr] + 0.01 * i
    end
    tri_ma_vr = zeros(nunit, cc.Ntri_setup)

    # ---- sequential isotope cascade reference (same fractions as C) ----
    function seq_iso(X0, Xin)
        X = copy(X0)
        for k in 1:ndecomp_cascade_transitions
            d = cascade_donor_pool[k]; r = cascade_receiver_pool[k]
            for j in 1:nlevdecomp
                loss = decomp_k[1, j, d] * X0[1, j, d]
                X[1, j, d] -= loss
                r != 0 && (X[1, j, r] += (1.0 - rf[1, j, k]) * pathfrac[1, j, k] * loss)
            end
        end
        for i in 1:ndecomp_pools, j in 1:nlevdecomp
            X[1, j, i] += Xin[1, j + nlevdecomp * (i - 1)]
        end
        return X
    end
    seq13 = seq_iso(C13, C13in); seq14 = seq_iso(C14, C14in)

    # ---- matrix advance (C/N + both isotopes in one call) ----
    matC = copy(C0); matN = copy(N0); mat13 = copy(C13); mat14 = copy(C14)
    ms = CLM.CNSoilMatrixState()
    CLM.cn_soil_matrix!(ms, cc;
        decomp_cpools_vr = matC, decomp_npools_vr = matN,
        Ksoil = copy(Ksoil), tri_ma_vr = tri_ma_vr,
        matrix_Cinput = matrix_Cinput, matrix_Ninput = matrix_Ninput,
        rf_decomp_cascade = rf, pathfrac_decomp_cascade = pathfrac,
        mask_soilc = [true], begc = begc, endc = endc,
        nlevdecomp = nlevdecomp, ndecomp_pools = ndecomp_pools,
        ndecomp_cascade_transitions = ndecomp_cascade_transitions,
        ndecomp_cascade_outtransitions = ndecomp_cascade_outtransitions, num_actfirec = 0,
        decomp_c13pools_vr = mat13, matrix_C13input = C13in,
        decomp_c14pools_vr = mat14, matrix_C14input = C14in)

    @test mat13 ≈ seq13 atol=1e-12 rtol=1e-12
    @test mat14 ≈ seq14 atol=1e-16 rtol=1e-10
    @test maximum(abs.(mat13 .- seq13)) < 1e-12
    # isotope advance is genuinely distinct from the C advance (different pool values).
    @test !(mat13 ≈ matC)

    # disabling isotopes (default nothing) leaves them untouched + still advances C.
    matC2 = copy(C0); matN2 = copy(N0); untouched = copy(C13)
    ms2 = CLM.CNSoilMatrixState()
    CLM.cn_soil_matrix!(ms2, cc;
        decomp_cpools_vr = matC2, decomp_npools_vr = matN2,
        Ksoil = copy(Ksoil), tri_ma_vr = tri_ma_vr,
        matrix_Cinput = matrix_Cinput, matrix_Ninput = matrix_Ninput,
        rf_decomp_cascade = rf, pathfrac_decomp_cascade = pathfrac,
        mask_soilc = [true], begc = begc, endc = endc,
        nlevdecomp = nlevdecomp, ndecomp_pools = ndecomp_pools,
        ndecomp_cascade_transitions = ndecomp_cascade_transitions,
        ndecomp_cascade_outtransitions = ndecomp_cascade_outtransitions, num_actfirec = 0)
    @test untouched == C13
    @test matC2 ≈ matC atol=1e-12

    # a missing isotope B-input is treated as zero (pure cascade advance).
    mat13b = copy(C13)
    ms3 = CLM.CNSoilMatrixState()
    CLM.cn_soil_matrix!(ms3, cc;
        decomp_cpools_vr = copy(C0), decomp_npools_vr = copy(N0),
        Ksoil = copy(Ksoil), tri_ma_vr = tri_ma_vr,
        matrix_Cinput = matrix_Cinput, matrix_Ninput = matrix_Ninput,
        rf_decomp_cascade = rf, pathfrac_decomp_cascade = pathfrac,
        mask_soilc = [true], begc = begc, endc = endc,
        nlevdecomp = nlevdecomp, ndecomp_pools = ndecomp_pools,
        ndecomp_cascade_transitions = ndecomp_cascade_transitions,
        ndecomp_cascade_outtransitions = ndecomp_cascade_outtransitions, num_actfirec = 0,
        decomp_c13pools_vr = mat13b, matrix_C13input = nothing)
    @test mat13b ≈ seq_iso(C13, zeros(nunit, ndecomp_pools_vr)) atol=1e-12
end

@testset "cn_soil_matrix_advance! carbon isotopes (B-input from isotope litter)" begin
    # The driver helper builds each isotope B-input from the isotope veg-carbon flux
    # litter fields and advances the isotope decomp pools alongside C/N. With transport
    # off it must equal the sequential cascade + isotope litter inputs.
    nc = 1; nlevdecomp = 2; ndecomp_pools = 3; dt = 1800.0
    i_litr_min = 1; i_litr_max = 2; i_cwd = 3; ndp_vr = ndecomp_pools * nlevdecomp
    initial_cn = [20.0, 15.0, 90.0]
    cascade_donor = [1, 2, 3]; cascade_recv = [2, 0, 1]
    decomp_k = [0.05, 0.02, 0.01]

    C0 = zeros(nc, nlevdecomp, ndecomp_pools); N0 = similar(C0); C13 = similar(C0)
    for i in 1:ndecomp_pools, j in 1:nlevdecomp
        C0[1, j, i] = 100.0 * i + 10.0 * j; N0[1, j, i] = C0[1, j, i] / initial_cn[i]
        C13[1, j, i] = 0.011 * C0[1, j, i] + 0.2 * j
    end
    # isotope litter input to litter pool 1 (phenology).
    phen13 = zeros(nc, nlevdecomp, i_litr_max)
    for j in 1:nlevdecomp; phen13[1, j, 1] = 0.006 + 0.001 * j; end

    # sequential reference: cascade (fractions) applied to C13 + phenology input.
    rf1 = 0.3
    seq13 = copy(C13)
    for k in 1:ndecomp_pools
        d = cascade_donor[k]; r = cascade_recv[k]
        rfk = (k == 1 ? rf1 : (k == 2 ? 1.0 : 0.0))
        for j in 1:nlevdecomp
            loss = decomp_k[d] * C13[1, j, d]
            seq13[1, j, d] -= loss
            r != 0 && (seq13[1, j, r] += (1.0 - rfk) * loss)
        end
    end
    for j in 1:nlevdecomp; seq13[1, j, 1] += phen13[1, j, 1] * dt; end

    # build the soil/veg data + isotope instances
    cs = CLM.SoilBiogeochemCarbonStateData(); cs.decomp_cpools_vr_col = copy(C0)
    ns = CLM.SoilBiogeochemNitrogenStateData(); ns.decomp_npools_vr_col = copy(N0)
    cf = CLM.SoilBiogeochemCarbonFluxData()
    cf.rf_decomp_cascade_col = zeros(nc, nlevdecomp, ndecomp_pools)
    cf.rf_decomp_cascade_col[1, :, 1] .= rf1; cf.rf_decomp_cascade_col[1, :, 2] .= 1.0
    cf.pathfrac_decomp_cascade_col = ones(nc, nlevdecomp, ndecomp_pools)
    cf.tri_ma_vr = zeros(nc, (nlevdecomp * 3 - 2) * (ndecomp_pools - 1))
    c13cs = CLM.SoilBiogeochemCarbonStateData(); c13cs.decomp_cpools_vr_col = copy(C13)

    vcf = CLM.CNVegCarbonFluxData(); CLM.cnveg_carbon_flux_init!(vcf, 1, nc, 1;
        nlevdecomp_full=nlevdecomp, ndecomp_pools=ndecomp_pools)
    vnf = CLM.CNVegNitrogenFluxData(); CLM.cnveg_nitrogen_flux_init!(vnf, 1, nc, 1;
        nlevdecomp_full=nlevdecomp, ndecomp_pools=ndecomp_pools, i_litr_max=i_litr_max)
    for f in (:phenology_c_to_litr_c_col, :gap_mortality_c_to_litr_c_col, :m_c_to_litr_fire_col,
              :gap_mortality_c_to_cwdc_col, :fire_mortality_c_to_cwdc_col); fill!(getfield(vcf, f), 0.0); end
    for f in (:phenology_n_to_litr_n_col, :gap_mortality_n_to_litr_n_col, :m_n_to_litr_fire_col,
              :gap_mortality_n_to_cwdn_col, :fire_mortality_n_to_cwdn_col); fill!(getfield(vnf, f), 0.0); end
    # isotope veg carbon flux carries the isotope litter fluxes.
    c13cf = CLM.CNVegCarbonFluxData(); CLM.cnveg_carbon_flux_init!(c13cf, 1, nc, 1;
        nlevdecomp_full=nlevdecomp, ndecomp_pools=ndecomp_pools)
    for f in (:phenology_c_to_litr_c_col, :gap_mortality_c_to_litr_c_col, :m_c_to_litr_fire_col,
              :gap_mortality_c_to_cwdc_col, :fire_mortality_c_to_cwdc_col); fill!(getfield(c13cf, f), 0.0); end
    c13cf.phenology_c_to_litr_c_col[:, :, 1:i_litr_max] .= phen13

    cc = CLM.DecompCascadeConData(cascade_donor_pool=cascade_donor, cascade_receiver_pool=cascade_recv,
        floating_cn_ratio_decomp_pools=BitVector([false,false,false]),
        is_cwd=BitVector([false,false,true]), initial_cn_ratio=initial_cn)
    Ksoil = zeros(nc, ndp_vr)
    for i in 1:ndecomp_pools, j in 1:nlevdecomp; Ksoil[1, j + (i - 1) * nlevdecomp] = decomp_k[i]; end

    CLM.cn_soil_matrix_advance!(cc, cs, ns, cf, vcf, vnf; Ksoil=Ksoil,
        mask_soilc=[true], bounds_col=1:nc, nlevdecomp=nlevdecomp, ndecomp_pools=ndecomp_pools,
        ndecomp_cascade_transitions=3, i_litr_min=i_litr_min, i_litr_max=i_litr_max, i_cwd=i_cwd, dt=dt,
        c13_soilbgc_cs=c13cs, c13_cnveg_cf=c13cf)

    @test c13cs.decomp_cpools_vr_col ≈ seq13 atol=1e-10 rtol=1e-12
    # C advanced too; no isotope supplied ⇒ that pool would be untouched (sanity: C13 changed).
    @test !(c13cs.decomp_cpools_vr_col ≈ C13)
end
