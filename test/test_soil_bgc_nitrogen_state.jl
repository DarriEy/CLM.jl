@testset "SoilBiogeochemNitrogenStateData" begin

    @testset "default construction" begin
        ns = CLM.SoilBiogeochemNitrogenStateData()
        @test isnan(ns.totvegcthresh)
        @test length(ns.ntrunc_col) == 0
        @test length(ns.sminn_col) == 0
        @test length(ns.totlitn_col) == 0
        @test length(ns.totn_grc) == 0
        @test size(ns.decomp_npools_vr_col) == (0, 0, 0)
        @test size(ns.decomp_npools_col) == (0, 0)
        @test size(ns.sminn_vr_col) == (0, 0)
    end

    @testset "init! basic allocation" begin
        nc = 5
        ng = 2
        nlev_full = 3
        npools = 4
        ns = CLM.SoilBiogeochemNitrogenStateData()
        CLM.soil_bgc_nitrogen_state_init!(ns, nc, ng, nlev_full, npools)

        # Column-level 1D
        @test length(ns.sminn_col) == nc
        @test all(isnan, ns.sminn_col)
        @test length(ns.ntrunc_col) == nc
        @test length(ns.cwdn_col) == nc
        @test length(ns.totlitn_col) == nc
        @test length(ns.totmicn_col) == nc
        @test length(ns.totsomn_col) == nc
        @test length(ns.totlitn_1m_col) == nc
        @test length(ns.totsomn_1m_col) == nc
        @test length(ns.dyn_nbal_adjustments_col) == nc
        @test length(ns.dyn_no3bal_adjustments_col) == nc
        @test length(ns.dyn_nh4bal_adjustments_col) == nc
        @test length(ns.smin_no3_col) == nc
        @test length(ns.smin_nh4_col) == nc
        @test length(ns.totn_col) == nc
        @test length(ns.totecosysn_col) == nc
        @test length(ns.totn_grc) == ng

        # Column × nlev 2D
        @test size(ns.sminn_vr_col) == (nc, nlev_full)
        @test all(isnan, ns.sminn_vr_col)
        @test size(ns.ntrunc_vr_col) == (nc, nlev_full)
        @test size(ns.smin_no3_vr_col) == (nc, nlev_full)
        @test size(ns.smin_nh4_vr_col) == (nc, nlev_full)
        @test size(ns.decomp_soiln_vr_col) == (nc, nlev_full)

        # Column × nlev × npools 3D
        @test size(ns.decomp_npools_vr_col) == (nc, nlev_full, npools)
        @test all(isnan, ns.decomp_npools_vr_col)

        # Column × npools 2D
        @test size(ns.decomp_npools_col) == (nc, npools)
        @test size(ns.decomp_npools_1m_col) == (nc, npools)

        # Matrix-CN fields should NOT be allocated
        @test size(ns.matrix_cap_decomp_npools_col) == (0, 0)
        @test size(ns.matrix_cap_decomp_npools_vr_col) == (0, 0, 0)
    end

    @testset "init! with use_soil_matrixcn" begin
        nc = 4
        ng = 1
        nlev_full = 3
        nlev = 2
        npools = 3
        ntrans = 2
        ns = CLM.SoilBiogeochemNitrogenStateData()
        CLM.soil_bgc_nitrogen_state_init!(ns, nc, ng, nlev_full, npools;
                                          nlevdecomp=nlev,
                                          ndecomp_cascade_transitions=ntrans,
                                          use_soil_matrixcn=true)

        @test size(ns.matrix_cap_decomp_npools_col) == (nc, npools)
        @test all(isnan, ns.matrix_cap_decomp_npools_col)
        @test size(ns.matrix_cap_decomp_npools_vr_col) == (nc, nlev_full, npools)
        @test size(ns.decomp0_npools_vr_col) == (nc, nlev_full, npools)
        @test size(ns.decomp_npools_vr_SASUsave_col) == (nc, nlev_full, npools)
        @test size(ns.in_nacc) == (nc, nlev * npools)
        @test size(ns.tran_nacc) == (nc, nlev * npools, nlev * npools)
        @test size(ns.in_nacc_2d) == (nc, nlev_full, npools)
        @test size(ns.vert_up_tran_nacc) == (nc, nlev_full, npools)
        @test size(ns.vert_down_tran_nacc) == (nc, nlev_full, npools)
        @test size(ns.exit_nacc) == (nc, nlev_full, npools)
        @test size(ns.hori_tran_nacc) == (nc, nlev_full, ntrans)
    end

    @testset "set_values! basic" begin
        nc = 4
        ng = 1
        nlev_full = 2
        npools = 3
        ns = CLM.SoilBiogeochemNitrogenStateData()
        CLM.soil_bgc_nitrogen_state_init!(ns, nc, ng, nlev_full, npools)

        mask = BitVector([true, false, true, false])
        CLM.soil_bgc_nitrogen_state_set_values!(ns, mask, 0.0)

        # Columns 1,3 should be zero; 2,4 should still be NaN
        @test ns.sminn_col[1] == 0.0
        @test isnan(ns.sminn_col[2])
        @test ns.ntrunc_col[1] == 0.0
        @test isnan(ns.ntrunc_col[2])
        @test ns.cwdn_col[3] == 0.0
        @test isnan(ns.cwdn_col[4])
        @test ns.totlitn_col[1] == 0.0
        @test ns.totmicn_col[1] == 0.0
        @test ns.totsomn_col[1] == 0.0
        @test ns.totsomn_1m_col[1] == 0.0
        @test ns.totlitn_1m_col[1] == 0.0
        @test ns.totecosysn_col[1] == 0.0
        @test ns.totn_col[1] == 0.0

        # vr fields
        @test ns.sminn_vr_col[1, 1] == 0.0
        @test isnan(ns.sminn_vr_col[2, 1])
        @test ns.ntrunc_vr_col[3, 2] == 0.0

        # pool fields
        @test ns.decomp_npools_col[1, 1] == 0.0
        @test ns.decomp_npools_1m_col[3, 2] == 0.0
        @test isnan(ns.decomp_npools_col[2, 1])

        # vr pool fields
        @test ns.decomp_npools_vr_col[1, 1, 1] == 0.0
        @test isnan(ns.decomp_npools_vr_col[4, 1, 1])
    end

    @testset "set_values! with nitrif_denitrif" begin
        nc = 3
        ng = 1
        nlev_full = 2
        npools = 2
        ns = CLM.SoilBiogeochemNitrogenStateData()
        CLM.soil_bgc_nitrogen_state_init!(ns, nc, ng, nlev_full, npools)

        mask = BitVector([true, true, true])
        CLM.soil_bgc_nitrogen_state_set_values!(ns, mask, 3.0;
                                                 use_nitrif_denitrif=true)

        @test ns.smin_no3_col[1] == 3.0
        @test ns.smin_nh4_col[2] == 3.0
        @test ns.smin_no3_vr_col[1, 1] == 3.0
        @test ns.smin_nh4_vr_col[3, 2] == 3.0
    end

    @testset "set_values! with matrixcn" begin
        nc = 3
        ng = 1
        nlev_full = 2
        nlev = 2
        npools = 2
        ntrans = 1
        ns = CLM.SoilBiogeochemNitrogenStateData()
        CLM.soil_bgc_nitrogen_state_init!(ns, nc, ng, nlev_full, npools;
                                          nlevdecomp=nlev,
                                          ndecomp_cascade_transitions=ntrans,
                                          use_soil_matrixcn=true)

        mask = BitVector([true, true, true])
        CLM.soil_bgc_nitrogen_state_set_values!(ns, mask, 5.0;
                                                 use_soil_matrixcn=true,
                                                 ndecomp_cascade_transitions=ntrans)

        @test ns.matrix_cap_decomp_npools_col[1, 1] == 5.0
        @test ns.matrix_cap_decomp_npools_vr_col[2, 1, 1] == 5.0
        @test ns.decomp0_npools_vr_col[3, 2, 2] == 5.0
        @test ns.in_nacc_2d[1, 1, 1] == 5.0
        @test ns.vert_up_tran_nacc[2, 2, 1] == 5.0
        @test ns.vert_down_tran_nacc[1, 1, 2] == 5.0
        @test ns.exit_nacc[3, 1, 1] == 5.0
        @test ns.hori_tran_nacc[1, 1, 1] == 5.0
    end

    @testset "init_cold! from C pools" begin
        nc = 3
        ng = 1
        nlev = 2
        nlev_full = 2
        npools = 3

        ns = CLM.SoilBiogeochemNitrogenStateData()
        CLM.soil_bgc_nitrogen_state_init!(ns, nc, ng, nlev_full, npools)

        # Mock C pool data
        cpools_vr = fill(100.0, nc, nlev_full, npools)
        cpools    = fill(50.0, nc, npools)
        cpools_1m = fill(50.0, nc, npools)
        cn_ratio  = [10.0, 20.0, 50.0]

        CLM.soil_bgc_nitrogen_state_init_cold!(ns, 1:nc;
                                                nlevdecomp=nlev,
                                                nlevdecomp_full=nlev_full,
                                                ndecomp_pools=npools,
                                                decomp_cpools_vr_col=cpools_vr,
                                                decomp_cpools_col=cpools,
                                                decomp_cpools_1m_col=cpools_1m,
                                                initial_cn_ratio=cn_ratio)

        # N = C / CN_ratio
        for c in 1:nc
            @test ns.decomp_npools_vr_col[c, 1, 1] ≈ 100.0 / 10.0  # 10.0
            @test ns.decomp_npools_vr_col[c, 1, 2] ≈ 100.0 / 20.0  # 5.0
            @test ns.decomp_npools_vr_col[c, 1, 3] ≈ 100.0 / 50.0  # 2.0
        end

        @test ns.decomp_npools_col[1, 1] ≈ 50.0 / 10.0
        @test ns.decomp_npools_col[1, 2] ≈ 50.0 / 20.0
        @test ns.decomp_npools_1m_col[1, 1] ≈ 50.0 / 10.0

        # sminn and ntrunc should be zero
        @test ns.sminn_vr_col[1, 1] == 0.0
        @test ns.ntrunc_vr_col[1, 1] == 0.0
        @test ns.sminn_col[1] == 0.0
        @test ns.ntrunc_col[1] == 0.0

        # Summary fields
        @test ns.cwdn_col[1] == 0.0
        @test ns.totlitn_col[1] == 0.0
        @test ns.totmicn_col[1] == 0.0
        @test ns.totsomn_col[1] == 0.0
        @test ns.totecosysn_col[1] == 0.0
        @test ns.totn_col[1] == 0.0
    end

    @testset "init_cold! with nitrif_denitrif" begin
        nc = 2
        ng = 1
        nlev = 2
        nlev_full = 2
        npools = 2

        ns = CLM.SoilBiogeochemNitrogenStateData()
        CLM.soil_bgc_nitrogen_state_init!(ns, nc, ng, nlev_full, npools)

        cpools_vr = fill(100.0, nc, nlev_full, npools)
        cpools    = fill(50.0, nc, npools)
        cpools_1m = fill(50.0, nc, npools)
        cn_ratio  = [10.0, 20.0]

        CLM.soil_bgc_nitrogen_state_init_cold!(ns, 1:nc;
                                                nlevdecomp=nlev,
                                                nlevdecomp_full=nlev_full,
                                                ndecomp_pools=npools,
                                                decomp_cpools_vr_col=cpools_vr,
                                                decomp_cpools_col=cpools,
                                                decomp_cpools_1m_col=cpools_1m,
                                                initial_cn_ratio=cn_ratio,
                                                use_nitrif_denitrif=true)

        @test ns.smin_nh4_vr_col[1, 1] == 0.0
        @test ns.smin_no3_vr_col[1, 1] == 0.0
        @test ns.smin_nh4_col[1] == 0.0
        @test ns.smin_no3_col[1] == 0.0
    end

    @testset "init_cold! with matrixcn" begin
        nc = 2
        ng = 1
        nlev = 2
        nlev_full = 2
        npools = 2
        ntrans = 1

        ns = CLM.SoilBiogeochemNitrogenStateData()
        CLM.soil_bgc_nitrogen_state_init!(ns, nc, ng, nlev_full, npools;
                                          nlevdecomp=nlev,
                                          ndecomp_cascade_transitions=ntrans,
                                          use_soil_matrixcn=true)

        cpools_vr = fill(100.0, nc, nlev_full, npools)
        cpools    = fill(50.0, nc, npools)
        cpools_1m = fill(50.0, nc, npools)
        cn_ratio  = [10.0, 20.0]

        CLM.soil_bgc_nitrogen_state_init_cold!(ns, 1:nc;
                                                nlevdecomp=nlev,
                                                nlevdecomp_full=nlev_full,
                                                ndecomp_pools=npools,
                                                ndecomp_cascade_transitions=ntrans,
                                                decomp_cpools_vr_col=cpools_vr,
                                                decomp_cpools_col=cpools,
                                                decomp_cpools_1m_col=cpools_1m,
                                                initial_cn_ratio=cn_ratio,
                                                use_soil_matrixcn=true)

        # Matrix capacity should match N pools
        @test ns.matrix_cap_decomp_npools_vr_col[1, 1, 1] ≈ 100.0 / 10.0
        @test ns.matrix_cap_decomp_npools_vr_col[1, 1, 2] ≈ 100.0 / 20.0
        @test ns.matrix_cap_decomp_npools_col[1, 1] ≈ 50.0 / 10.0

        # decomp0 = max(pool, 1e-30)
        @test ns.decomp0_npools_vr_col[1, 1, 1] ≈ 10.0
        @test ns.decomp_npools_vr_SASUsave_col[1, 1, 1] == 0.0

        # Accumulated fields should be zero
        @test ns.in_nacc_2d[1, 1, 1] == 0.0
        @test ns.vert_up_tran_nacc[1, 1, 1] == 0.0
        @test ns.vert_down_tran_nacc[1, 1, 1] == 0.0
        @test ns.exit_nacc[1, 1, 1] == 0.0
        @test ns.hori_tran_nacc[1, 1, 1] == 0.0
        @test ns.in_nacc[1, 1] == 0.0
    end

    @testset "init_cold! with nlevdecomp < nlevdecomp_full" begin
        nc = 2
        ng = 1
        nlev = 2
        nlev_full = 4
        npools = 2

        ns = CLM.SoilBiogeochemNitrogenStateData()
        CLM.soil_bgc_nitrogen_state_init!(ns, nc, ng, nlev_full, npools)

        cpools_vr = fill(100.0, nc, nlev_full, npools)
        cpools    = fill(50.0, nc, npools)
        cpools_1m = fill(50.0, nc, npools)
        cn_ratio  = [10.0, 20.0]

        CLM.soil_bgc_nitrogen_state_init_cold!(ns, 1:nc;
                                                nlevdecomp=nlev,
                                                nlevdecomp_full=nlev_full,
                                                ndecomp_pools=npools,
                                                decomp_cpools_vr_col=cpools_vr,
                                                decomp_cpools_col=cpools,
                                                decomp_cpools_1m_col=cpools_1m,
                                                initial_cn_ratio=cn_ratio)

        # Active levels have N = C / CN_ratio
        @test ns.decomp_npools_vr_col[1, 1, 1] ≈ 10.0
        @test ns.decomp_npools_vr_col[1, 2, 1] ≈ 10.0

        # Extra levels beyond nlevdecomp should be zero
        @test ns.decomp_npools_vr_col[1, 3, 1] == 0.0
        @test ns.decomp_npools_vr_col[1, 4, 1] == 0.0
        @test ns.sminn_vr_col[1, 3] == 0.0
        @test ns.ntrunc_vr_col[1, 4] == 0.0
    end

    @testset "summary! basic integration" begin
        nc = 3
        ng = 1
        nlev = 2
        nlev_full = 2
        npools = 3
        dzsoi = [0.5, 0.5]
        zisoi = [0.0, 0.5, 1.0]

        ns = CLM.SoilBiogeochemNitrogenStateData()
        CLM.soil_bgc_nitrogen_state_init!(ns, nc, ng, nlev_full, npools)

        # Set vr pools to known values
        ns.decomp_npools_vr_col .= 10.0  # 10 gN/m3 everywhere
        ns.ntrunc_vr_col .= 1.0          # 1 gN/m3 truncation
        ns.sminn_vr_col .= 2.0           # 2 gN/m3 mineral N

        mask = BitVector([true, true, true])
        is_litter  = BitVector([true, false, false])
        is_soil    = BitVector([false, true, false])
        is_microbe = BitVector([false, false, true])
        is_cwd     = BitVector([false, false, false])
        totn_p2c   = zeros(nc)
        totvegn    = zeros(nc)

        CLM.soil_bgc_nitrogen_state_summary!(ns, mask, 1:nc;
                                              nlevdecomp=nlev,
                                              nlevdecomp_full=nlev_full,
                                              ndecomp_pools=npools,
                                              dzsoi_decomp_vals=dzsoi,
                                              zisoi_vals=zisoi,
                                              is_litter=is_litter,
                                              is_soil=is_soil,
                                              is_microbe=is_microbe,
                                              is_cwd=is_cwd,
                                              totn_p2c_col=totn_p2c,
                                              totvegn_col=totvegn)

        # Each pool integrates to 10.0 * (0.5 + 0.5) = 10.0 gN/m2
        for c in 1:nc
            for l in 1:npools
                @test ns.decomp_npools_col[c, l] ≈ 10.0
            end
        end

        # sminn integrates to 2.0 * (0.5 + 0.5) = 2.0 gN/m2
        @test ns.sminn_col[1] ≈ 2.0

        # ntrunc integrates to 1.0 * (0.5 + 0.5) = 1.0 gN/m2
        @test ns.ntrunc_col[1] ≈ 1.0

        # totlitn = pool 1 = 10.0
        @test ns.totlitn_col[1] ≈ 10.0

        # totsomn = pool 2 = 10.0
        @test ns.totsomn_col[1] ≈ 10.0

        # totmicn = pool 3 = 10.0
        @test ns.totmicn_col[1] ≈ 10.0

        # 1m integration (zisoi[2]=0.5 < 1, zisoi[3]=1.0 <= 1)
        for c in 1:nc
            for l in 1:npools
                @test ns.decomp_npools_1m_col[c, l] ≈ 10.0
            end
        end

        # totlitn_1m
        @test ns.totlitn_1m_col[1] ≈ 10.0
        @test ns.totsomn_1m_col[1] ≈ 10.0

        # decomp_soiln_vr = sum of soil pools per level
        @test ns.decomp_soiln_vr_col[1, 1] ≈ 10.0  # only pool 2 is soil

        # totecosysn = cwdn + totlitn + totmicn + totsomn + sminn + ecovegn
        # = 0 + 10 + 10 + 10 + 2 + 0 = 32
        @test ns.totecosysn_col[1] ≈ 32.0

        # totn = cwdn + totlitn + totmicn + totsomn + sminn + ntrunc + tvegn
        # = 0 + 10 + 10 + 10 + 2 + 1 + 0 = 33
        @test ns.totn_col[1] ≈ 33.0
    end

    @testset "summary! with nitrif_denitrif" begin
        nc = 2
        ng = 1
        nlev = 1
        nlev_full = 1
        npools = 1
        dzsoi = [1.0]

        ns = CLM.SoilBiogeochemNitrogenStateData()
        CLM.soil_bgc_nitrogen_state_init!(ns, nc, ng, nlev_full, npools)

        ns.decomp_npools_vr_col .= 5.0
        ns.ntrunc_vr_col .= 0.0
        ns.sminn_vr_col .= 0.0
        ns.smin_no3_vr_col .= 3.0
        ns.smin_nh4_vr_col .= 4.0

        mask = BitVector([true, true])
        is_litter  = BitVector([false])
        is_soil    = BitVector([false])
        is_microbe = BitVector([false])
        is_cwd     = BitVector([false])

        CLM.soil_bgc_nitrogen_state_summary!(ns, mask, 1:nc;
                                              nlevdecomp=nlev,
                                              nlevdecomp_full=nlev_full,
                                              ndecomp_pools=npools,
                                              dzsoi_decomp_vals=dzsoi,
                                              is_litter=is_litter,
                                              is_soil=is_soil,
                                              is_microbe=is_microbe,
                                              is_cwd=is_cwd,
                                              use_nitrif_denitrif=true)

        # NO3 integrates to 3.0 * 1.0 = 3.0
        @test ns.smin_no3_col[1] ≈ 3.0
        # NH4 integrates to 4.0 * 1.0 = 4.0
        @test ns.smin_nh4_col[1] ≈ 4.0
    end

    @testset "summary! with CWD pool" begin
        nc = 2
        ng = 1
        nlev = 1
        nlev_full = 1
        npools = 2
        dzsoi = [1.0]

        ns = CLM.SoilBiogeochemNitrogenStateData()
        CLM.soil_bgc_nitrogen_state_init!(ns, nc, ng, nlev_full, npools)

        ns.decomp_npools_vr_col .= 5.0
        ns.ntrunc_vr_col .= 0.0
        ns.sminn_vr_col .= 1.0

        mask = BitVector([true, true])
        is_litter  = BitVector([false, false])
        is_soil    = BitVector([false, false])
        is_microbe = BitVector([false, false])
        is_cwd     = BitVector([true, false])
        totn_p2c   = [100.0, 100.0]
        totvegn    = [50.0, 50.0]

        CLM.soil_bgc_nitrogen_state_summary!(ns, mask, 1:nc;
                                              nlevdecomp=nlev,
                                              nlevdecomp_full=nlev_full,
                                              ndecomp_pools=npools,
                                              dzsoi_decomp_vals=dzsoi,
                                              is_litter=is_litter,
                                              is_soil=is_soil,
                                              is_microbe=is_microbe,
                                              is_cwd=is_cwd,
                                              totn_p2c_col=totn_p2c,
                                              totvegn_col=totvegn)

        # cwdn = pool 1 (cwd) integrated = 5.0 * 1.0 = 5.0
        @test ns.cwdn_col[1] ≈ 5.0

        # sminn = 1.0 * 1.0 = 1.0
        @test ns.sminn_col[1] ≈ 1.0

        # totecosysn = cwdn(5) + totlitn(0) + totmicn(0) + totsomn(0) + sminn(1) + ecovegn(50) = 56
        @test ns.totecosysn_col[1] ≈ 56.0

        # totn = cwdn(5) + totlitn(0) + totmicn(0) + totsomn(0) + sminn(1) + ntrunc(0) + totvegn(100) = 106
        @test ns.totn_col[1] ≈ 106.0
    end

    @testset "set_totvegcthresh!" begin
        ns = CLM.SoilBiogeochemNitrogenStateData()
        CLM.soil_bgc_nitrogen_state_set_totvegcthresh!(ns, 100.0)
        @test ns.totvegcthresh == 100.0

        @test_throws ErrorException CLM.soil_bgc_nitrogen_state_set_totvegcthresh!(ns, 0.0)
        @test_throws ErrorException CLM.soil_bgc_nitrogen_state_set_totvegcthresh!(ns, -1.0)
    end

    @testset "stub functions run without error" begin
        ns = CLM.SoilBiogeochemNitrogenStateData()
        CLM.soil_bgc_nitrogen_state_init!(ns, 5, 2, 3, 2)

        @test CLM.soil_bgc_nitrogen_state_init_history!(ns, 1:5) === nothing
        @test CLM.soil_bgc_nitrogen_state_restart!(ns, 1:5) === nothing
        @test CLM.soil_bgc_nitrogen_state_dynamic_col_adjustments!(ns, 1:5) === nothing
    end

    @testset "field mutability" begin
        ns = CLM.SoilBiogeochemNitrogenStateData()
        CLM.soil_bgc_nitrogen_state_init!(ns, 3, 1, 2, 2)

        ns.sminn_col[1] = 42.0
        @test ns.sminn_col[1] == 42.0

        ns.decomp_npools_vr_col[1, 1, 1] = 99.0
        @test ns.decomp_npools_vr_col[1, 1, 1] == 99.0

        ns.totn_grc[1] = 7.5
        @test ns.totn_grc[1] == 7.5
    end

    @testset "re-init overwrites previous state" begin
        ns = CLM.SoilBiogeochemNitrogenStateData()
        CLM.soil_bgc_nitrogen_state_init!(ns, 3, 1, 2, 2)
        ns.sminn_col[1] = 999.0

        CLM.soil_bgc_nitrogen_state_init!(ns, 5, 2, 4, 3)
        @test length(ns.sminn_col) == 5
        @test all(isnan, ns.sminn_col)
        @test size(ns.decomp_npools_vr_col) == (5, 4, 3)
        @test length(ns.totn_grc) == 2
    end

end
