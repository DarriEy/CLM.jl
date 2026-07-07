@testset "cn_soil_matrix_advance! (driver soil-matrix step == sequential)" begin
    # Phase-3 driver integration: cn_soil_matrix_advance! composes init_soil_transfer!
    # + the B-input accumulation + cn_soil_matrix! into the one step the driver calls
    # in matrix mode. With vertical transport off (tri_ma_vr=0), the advance must
    # reproduce EXACTLY the sequential decomposition cascade PLUS the litter/CWD inputs.
    nc = 1; nlevdecomp = 2; ndecomp_pools = 3; dt = 1800.0
    i_litr_min = 1; i_litr_max = 2; i_cwd = 3
    ndp_vr = ndecomp_pools * nlevdecomp
    # cascade: 1(litr)->2(som); 2->atm(0); 3(cwd)->1.
    cascade_donor = [1, 2, 3]; cascade_recv = [2, 0, 1]
    initial_cn = [20.0, 15.0, 90.0]

    decomp_k = zeros(nc, nlevdecomp, ndecomp_pools)
    decomp_k[1, :, 1] .= 0.05; decomp_k[1, :, 2] .= 0.02; decomp_k[1, :, 3] .= 0.01
    rf = zeros(nc, nlevdecomp, ndecomp_pools); rf[1, :, 1] .= 0.3; rf[1, :, 2] .= 1.0; rf[1, :, 3] .= 0.0
    pf = ones(nc, nlevdecomp, ndecomp_pools)

    C0 = zeros(nc, nlevdecomp, ndecomp_pools); N0 = similar(C0)
    for i in 1:ndecomp_pools, j in 1:nlevdecomp
        C0[1, j, i] = 100.0 * i + 10.0 * j; N0[1, j, i] = C0[1, j, i] / initial_cn[i]
    end

    # veg→soil litter inputs (per pool/level for litter; per level for CWD).
    phen_c = zeros(nc, nlevdecomp, i_litr_max); gap_cwd = zeros(nc, nlevdecomp)
    for j in 1:nlevdecomp, i in i_litr_min:i_litr_max; phen_c[1, j, i] = 0.5 * i + 0.1 * j; end
    for j in 1:nlevdecomp; gap_cwd[1, j] = 0.3 + 0.05 * j; end
    phen_n = phen_c ./ 12.0; gap_cwd_n = gap_cwd ./ 12.0

    # ---- (a) sequential reference: cascade + litter inputs ----
    seqC = copy(C0); seqN = copy(N0)
    for k in 1:ndecomp_pools
        d = cascade_donor[k]; r = cascade_recv[k]
        for j in 1:nlevdecomp
            closs = decomp_k[1, j, d] * C0[1, j, d]; nloss = decomp_k[1, j, d] * N0[1, j, d]
            seqC[1, j, d] -= closs; seqN[1, j, d] -= nloss
            if r != 0
                seqC[1, j, r] += (1 - rf[1, j, k]) * pf[1, j, k] * closs
                seqN[1, j, r] += (1 - rf[1, j, k]) * pf[1, j, k] * nloss * (initial_cn[d] / initial_cn[r])
            end
        end
    end
    for j in 1:nlevdecomp
        for i in i_litr_min:i_litr_max
            seqC[1, j, i] += phen_c[1, j, i] * dt; seqN[1, j, i] += phen_n[1, j, i] * dt
        end
        seqC[1, j, i_cwd] += gap_cwd[1, j] * dt; seqN[1, j, i_cwd] += gap_cwd_n[1, j] * dt
    end

    # ---- (b) driver helper ----
    cs = CLM.SoilBiogeochemCarbonStateData(); cs.decomp_cpools_vr_col = copy(C0)
    ns = CLM.SoilBiogeochemNitrogenStateData(); ns.decomp_npools_vr_col = copy(N0)
    cf = CLM.SoilBiogeochemCarbonFluxData()
    cf.rf_decomp_cascade_col = rf; cf.pathfrac_decomp_cascade_col = pf
    cf.tri_ma_vr = zeros(nc, (nlevdecomp * 3 - 2) * (ndecomp_pools - 1))   # transport off
    vcf = CLM.CNVegCarbonFluxData(); CLM.cnveg_carbon_flux_init!(vcf, 1, nc, 1;
        nlevdecomp_full=nlevdecomp, ndecomp_pools=ndecomp_pools)
    vnf = CLM.CNVegNitrogenFluxData(); CLM.cnveg_nitrogen_flux_init!(vnf, 1, nc, 1;
        nlevdecomp_full=nlevdecomp, ndecomp_pools=ndecomp_pools, i_litr_max=i_litr_max)
    # zero all litter/CWD input fields, then set phenology (litr) + gap (cwd).
    for f in (:phenology_c_to_litr_c_col, :dwt_frootc_to_litr_c_col, :gap_mortality_c_to_litr_c_col, :m_c_to_litr_fire_col)
        fill!(getfield(vcf, f), 0.0)
    end
    for f in (:dwt_livecrootc_to_cwdc_col, :dwt_deadcrootc_to_cwdc_col, :gap_mortality_c_to_cwdc_col, :fire_mortality_c_to_cwdc_col)
        fill!(getfield(vcf, f), 0.0)
    end
    for f in (:phenology_n_to_litr_n_col, :dwt_frootn_to_litr_n_col, :gap_mortality_n_to_litr_n_col, :m_n_to_litr_fire_col)
        fill!(getfield(vnf, f), 0.0)
    end
    for f in (:dwt_livecrootn_to_cwdn_col, :dwt_deadcrootn_to_cwdn_col, :gap_mortality_n_to_cwdn_col, :fire_mortality_n_to_cwdn_col)
        fill!(getfield(vnf, f), 0.0)
    end
    vcf.phenology_c_to_litr_c_col[:, :, 1:i_litr_max] .= phen_c
    vcf.gap_mortality_c_to_cwdc_col .= gap_cwd
    vnf.phenology_n_to_litr_n_col[:, :, 1:i_litr_max] .= phen_n
    vnf.gap_mortality_n_to_cwdn_col .= gap_cwd_n

    cc = CLM.DecompCascadeConData(cascade_donor_pool=cascade_donor, cascade_receiver_pool=cascade_recv,
        floating_cn_ratio_decomp_pools=BitVector([false, false, false]),
        is_cwd=BitVector([false, false, true]), initial_cn_ratio=initial_cn)

    # Ksoil is the per-step turnover fraction (the driver forms it as decomp_k_col·dt;
    # here decomp_k is already per-step, matching the sequential cascade above).
    Ksoil = zeros(nc, ndp_vr)
    for i in 1:ndecomp_pools, j in 1:nlevdecomp; Ksoil[1, j + (i - 1) * nlevdecomp] = decomp_k[1, j, i]; end

    CLM.cn_soil_matrix_advance!(cc, cs, ns, cf, vcf, vnf; Ksoil=Ksoil,
        mask_soilc=[true], bounds_col=1:nc, nlevdecomp=nlevdecomp, ndecomp_pools=ndecomp_pools,
        ndecomp_cascade_transitions=3, i_litr_min=i_litr_min, i_litr_max=i_litr_max, i_cwd=i_cwd, dt=dt)

    @test cs.decomp_cpools_vr_col ≈ seqC atol=1e-10 rtol=1e-12
    @test ns.decomp_npools_vr_col ≈ seqN atol=1e-10 rtol=1e-12
    @test cc.n_all_entries != CLM.SMM_EMPTY_INT   # init_soil_transfer! ran

    # second call reuses the (now-initialised) cascade_con index structure.
    cs2 = CLM.SoilBiogeochemCarbonStateData(); cs2.decomp_cpools_vr_col = copy(C0)
    ns2 = CLM.SoilBiogeochemNitrogenStateData(); ns2.decomp_npools_vr_col = copy(N0)
    CLM.cn_soil_matrix_advance!(cc, cs2, ns2, cf, vcf, vnf; Ksoil=Ksoil,
        mask_soilc=[true], bounds_col=1:nc, nlevdecomp=nlevdecomp, ndecomp_pools=ndecomp_pools,
        ndecomp_cascade_transitions=3, i_litr_min=i_litr_min, i_litr_max=i_litr_max, i_cwd=i_cwd, dt=dt)
    @test cs2.decomp_cpools_vr_col ≈ seqC atol=1e-10
end
