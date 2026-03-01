@testset "SoilBiogeochemCarbonStateData" begin

    @testset "default construction" begin
        cs = CLM.SoilBiogeochemCarbonStateData()
        @test cs.restart_file_spinup_state == typemax(Int)
        @test isnan(cs.totvegcthresh)
        @test length(cs.ctrunc_col) == 0
        @test length(cs.totlitc_col) == 0
        @test length(cs.totc_grc) == 0
        @test size(cs.decomp_cpools_vr_col) == (0, 0, 0)
        @test size(cs.decomp_cpools_col) == (0, 0)
    end

    @testset "init! basic allocation" begin
        nc = 5
        ng = 2
        nlev_full = 3
        npools = 4
        cs = CLM.SoilBiogeochemCarbonStateData()
        CLM.soil_bgc_carbon_state_init!(cs, nc, ng, nlev_full, npools)

        # Column-level 1D
        @test length(cs.ctrunc_col) == nc
        @test all(isnan, cs.ctrunc_col)
        @test length(cs.cwdc_col) == nc
        @test length(cs.totmicc_col) == nc
        @test length(cs.totlitc_col) == nc
        @test length(cs.totsomc_col) == nc
        @test length(cs.totlitc_1m_col) == nc
        @test length(cs.totsomc_1m_col) == nc
        @test length(cs.dyn_cbal_adjustments_col) == nc
        @test length(cs.totc_col) == nc
        @test length(cs.totecosysc_col) == nc
        @test length(cs.totc_grc) == ng

        # Column × nlev 2D
        @test size(cs.ctrunc_vr_col) == (nc, nlev_full)
        @test all(isnan, cs.ctrunc_vr_col)
        @test size(cs.decomp_soilc_vr_col) == (nc, nlev_full)

        # Column × nlev × npools 3D
        @test size(cs.decomp_cpools_vr_col) == (nc, nlev_full, npools)
        @test all(isnan, cs.decomp_cpools_vr_col)

        # Column × npools 2D
        @test size(cs.decomp_cpools_col) == (nc, npools)
        @test size(cs.decomp_cpools_1m_col) == (nc, npools)

        # Matrix-CN fields should NOT be allocated
        @test size(cs.matrix_cap_decomp_cpools_col) == (0, 0)
        @test size(cs.matrix_cap_decomp_cpools_vr_col) == (0, 0, 0)
    end

    @testset "init! with use_soil_matrixcn" begin
        nc = 4
        ng = 1
        nlev_full = 3
        nlev = 2
        npools = 3
        ntrans = 2
        cs = CLM.SoilBiogeochemCarbonStateData()
        CLM.soil_bgc_carbon_state_init!(cs, nc, ng, nlev_full, npools;
                                        nlevdecomp=nlev,
                                        ndecomp_cascade_transitions=ntrans,
                                        use_soil_matrixcn=true)

        @test size(cs.matrix_cap_decomp_cpools_col) == (nc, npools)
        @test all(isnan, cs.matrix_cap_decomp_cpools_col)
        @test size(cs.matrix_cap_decomp_cpools_vr_col) == (nc, nlev_full, npools)
        @test size(cs.decomp0_cpools_vr_col) == (nc, nlev_full, npools)
        @test size(cs.decomp_cpools_vr_SASUsave_col) == (nc, nlev_full, npools)
        @test size(cs.in_acc) == (nc, nlev * npools)
        @test size(cs.tran_acc) == (nc, nlev * npools, nlev * npools)
        @test size(cs.in_acc_2d) == (nc, nlev_full, npools)
        @test size(cs.vert_up_tran_acc) == (nc, nlev_full, npools)
        @test size(cs.vert_down_tran_acc) == (nc, nlev_full, npools)
        @test size(cs.exit_acc) == (nc, nlev_full, npools)
        @test size(cs.hori_tran_acc) == (nc, nlev_full, ntrans)
    end

    @testset "set_values! basic" begin
        nc = 4
        ng = 1
        nlev_full = 2
        npools = 3
        cs = CLM.SoilBiogeochemCarbonStateData()
        CLM.soil_bgc_carbon_state_init!(cs, nc, ng, nlev_full, npools)

        mask = BitVector([true, false, true, false])
        CLM.soil_bgc_carbon_state_set_values!(cs, mask, 0.0)

        # Columns 1,3 should be zero; 2,4 should still be NaN
        @test cs.ctrunc_col[1] == 0.0
        @test isnan(cs.ctrunc_col[2])
        @test cs.cwdc_col[3] == 0.0
        @test isnan(cs.cwdc_col[4])
        @test cs.totmicc_col[1] == 0.0
        @test cs.totlitc_col[1] == 0.0
        @test cs.totsomc_col[1] == 0.0
        @test cs.totlitc_1m_col[1] == 0.0
        @test cs.totsomc_1m_col[1] == 0.0
        @test cs.totc_col[1] == 0.0
        @test cs.totecosysc_col[1] == 0.0

        # vr fields
        @test cs.ctrunc_vr_col[1, 1] == 0.0
        @test isnan(cs.ctrunc_vr_col[2, 1])
        @test cs.ctrunc_vr_col[3, 2] == 0.0

        # pool fields
        @test cs.decomp_cpools_col[1, 1] == 0.0
        @test cs.decomp_cpools_1m_col[3, 2] == 0.0
        @test isnan(cs.decomp_cpools_col[2, 1])

        # vr pool fields
        @test cs.decomp_cpools_vr_col[1, 1, 1] == 0.0
        @test isnan(cs.decomp_cpools_vr_col[4, 1, 1])
    end

    @testset "set_values! with matrixcn" begin
        nc = 3
        ng = 1
        nlev_full = 2
        nlev = 2
        npools = 2
        ntrans = 1
        cs = CLM.SoilBiogeochemCarbonStateData()
        CLM.soil_bgc_carbon_state_init!(cs, nc, ng, nlev_full, npools;
                                        nlevdecomp=nlev,
                                        ndecomp_cascade_transitions=ntrans,
                                        use_soil_matrixcn=true)

        mask = BitVector([true, true, true])
        CLM.soil_bgc_carbon_state_set_values!(cs, mask, 5.0;
                                              use_soil_matrixcn=true,
                                              ndecomp_cascade_transitions=ntrans)

        @test cs.matrix_cap_decomp_cpools_col[1, 1] == 5.0
        @test cs.matrix_cap_decomp_cpools_vr_col[2, 1, 1] == 5.0
        @test cs.decomp0_cpools_vr_col[3, 2, 2] == 5.0
        @test cs.in_acc_2d[1, 1, 1] == 5.0
        @test cs.vert_up_tran_acc[2, 2, 1] == 5.0
        @test cs.vert_down_tran_acc[1, 1, 2] == 5.0
        @test cs.exit_acc[3, 1, 1] == 5.0
        @test cs.hori_tran_acc[1, 1, 1] == 5.0
    end

    @testset "init_cold! c12 basic" begin
        nc = 3
        ng = 1
        nlev = 2
        nlev_full = 2
        npools = 3
        initial_stock = [10.0, 20.0, 30.0]

        cs = CLM.SoilBiogeochemCarbonStateData()
        CLM.soil_bgc_carbon_state_init!(cs, nc, ng, nlev_full, npools)

        zsoi = [0.1, 0.5]  # first level below threshold, second above
        CLM.soil_bgc_carbon_state_init_cold!(cs, 1:nc, 1.0;
                                             nlevdecomp=nlev,
                                             nlevdecomp_full=nlev_full,
                                             ndecomp_pools=npools,
                                             initial_stock=initial_stock,
                                             initial_stock_soildepth=0.3,
                                             zsoi_vals=zsoi)

        # First level (z=0.1 < 0.3) should have initial_stock
        for c in 1:nc
            @test cs.decomp_cpools_vr_col[c, 1, 1] == 10.0
            @test cs.decomp_cpools_vr_col[c, 1, 2] == 20.0
            @test cs.decomp_cpools_vr_col[c, 1, 3] == 30.0
        end

        # Second level (z=0.5 >= 0.3) should be zero
        for c in 1:nc
            @test cs.decomp_cpools_vr_col[c, 2, 1] == 0.0
            @test cs.decomp_cpools_vr_col[c, 2, 2] == 0.0
            @test cs.decomp_cpools_vr_col[c, 2, 3] == 0.0
        end

        # ctrunc_vr should be zero
        @test cs.ctrunc_vr_col[1, 1] == 0.0
        @test cs.ctrunc_vr_col[1, 2] == 0.0

        # Integrated pools should match initial_stock
        @test cs.decomp_cpools_col[1, 1] == 10.0
        @test cs.decomp_cpools_col[1, 2] == 20.0
        @test cs.decomp_cpools_col[1, 3] == 30.0
        @test cs.decomp_cpools_1m_col[1, 1] == 10.0

        # Summary fields
        @test cs.cwdc_col[1] == 0.0
        @test cs.ctrunc_col[1] == 0.0
        @test cs.totmicc_col[1] == 0.0
        @test cs.totlitc_col[1] == 0.0
        @test cs.totsomc_col[1] == 0.0
        @test cs.totc_col[1] == 0.0
        @test cs.totecosysc_col[1] == 0.0
    end

    @testset "init_cold! c13/c14 from c12 instance" begin
        nc = 2
        ng = 1
        nlev = 2
        nlev_full = 2
        npools = 2

        # Set up c12 instance
        c12 = CLM.SoilBiogeochemCarbonStateData()
        CLM.soil_bgc_carbon_state_init!(c12, nc, ng, nlev_full, npools)
        c12.decomp_cpools_vr_col .= 100.0
        c12.ctrunc_vr_col .= 1.0
        c12.decomp_cpools_col .= 50.0
        c12.decomp_cpools_1m_col .= 50.0
        c12.cwdc_col .= 5.0

        # Set up c13 instance
        cs = CLM.SoilBiogeochemCarbonStateData()
        CLM.soil_bgc_carbon_state_init!(cs, nc, ng, nlev_full, npools)

        ratio = 0.5
        CLM.soil_bgc_carbon_state_init_cold!(cs, 1:nc, ratio;
                                             nlevdecomp=nlev,
                                             nlevdecomp_full=nlev_full,
                                             ndecomp_pools=npools,
                                             c12_inst=c12)

        @test cs.decomp_cpools_vr_col[1, 1, 1] ≈ 50.0
        @test cs.ctrunc_vr_col[1, 1] ≈ 0.5
        @test cs.decomp_cpools_col[1, 1] ≈ 25.0
        @test cs.decomp_cpools_1m_col[1, 1] ≈ 25.0
        @test cs.cwdc_col[1] ≈ 2.5
    end

    @testset "init_cold! with matrixcn" begin
        nc = 2
        ng = 1
        nlev = 2
        nlev_full = 2
        npools = 2
        ntrans = 1
        initial_stock = [10.0, 20.0]
        zsoi = [0.1, 0.2]  # both below threshold

        cs = CLM.SoilBiogeochemCarbonStateData()
        CLM.soil_bgc_carbon_state_init!(cs, nc, ng, nlev_full, npools;
                                        nlevdecomp=nlev,
                                        ndecomp_cascade_transitions=ntrans,
                                        use_soil_matrixcn=true)

        CLM.soil_bgc_carbon_state_init_cold!(cs, 1:nc, 1.0;
                                             nlevdecomp=nlev,
                                             nlevdecomp_full=nlev_full,
                                             ndecomp_pools=npools,
                                             ndecomp_cascade_transitions=ntrans,
                                             initial_stock=initial_stock,
                                             initial_stock_soildepth=0.3,
                                             zsoi_vals=zsoi,
                                             use_soil_matrixcn=true)

        # Matrix capacity should match initial stock
        @test cs.matrix_cap_decomp_cpools_vr_col[1, 1, 1] == 10.0
        @test cs.matrix_cap_decomp_cpools_vr_col[1, 1, 2] == 20.0

        # decomp0 = max(pool, 1e-30)
        @test cs.decomp0_cpools_vr_col[1, 1, 1] == 10.0
        @test cs.decomp_cpools_vr_SASUsave_col[1, 1, 1] == 0.0

        # Accumulated fields should be zero
        @test cs.in_acc_2d[1, 1, 1] == 0.0
        @test cs.vert_up_tran_acc[1, 1, 1] == 0.0
        @test cs.vert_down_tran_acc[1, 1, 1] == 0.0
        @test cs.exit_acc[1, 1, 1] == 0.0
        @test cs.hori_tran_acc[1, 1, 1] == 0.0
    end

    @testset "summary! basic integration" begin
        nc = 3
        ng = 1
        nlev = 2
        nlev_full = 2
        npools = 3
        dzsoi = [0.5, 0.5]  # uniform thickness
        zisoi = [0.0, 0.5, 1.0]  # interface depths

        cs = CLM.SoilBiogeochemCarbonStateData()
        CLM.soil_bgc_carbon_state_init!(cs, nc, ng, nlev_full, npools)

        # Set vr pools to known values
        cs.decomp_cpools_vr_col .= 10.0  # 10 gC/m3 everywhere
        cs.ctrunc_vr_col .= 1.0          # 1 gC/m3 truncation

        mask = BitVector([true, true, true])
        is_litter  = BitVector([true, false, false])
        is_soil    = BitVector([false, true, false])
        is_microbe = BitVector([false, false, true])
        is_cwd     = BitVector([false, false, false])
        totc_p2c   = zeros(nc)
        totvegc    = zeros(nc)

        CLM.soil_bgc_carbon_state_summary!(cs, mask, 1:nc;
                                           nlevdecomp=nlev,
                                           nlevdecomp_full=nlev_full,
                                           ndecomp_pools=npools,
                                           dzsoi_decomp_vals=dzsoi,
                                           zisoi_vals=zisoi,
                                           is_litter=is_litter,
                                           is_soil=is_soil,
                                           is_microbe=is_microbe,
                                           is_cwd=is_cwd,
                                           totc_p2c_col=totc_p2c,
                                           totvegc_col=totvegc)

        # Each pool integrates to 10.0 * (0.5 + 0.5) = 10.0 gC/m2
        for c in 1:nc
            for l in 1:npools
                @test cs.decomp_cpools_col[c, l] ≈ 10.0
            end
        end

        # ctrunc integrates to 1.0 * (0.5 + 0.5) = 1.0 gC/m2
        @test cs.ctrunc_col[1] ≈ 1.0

        # totlitc = pool 1 = 10.0
        @test cs.totlitc_col[1] ≈ 10.0

        # totsomc = pool 2 = 10.0
        @test cs.totsomc_col[1] ≈ 10.0

        # totmicc = pool 3 = 10.0
        @test cs.totmicc_col[1] ≈ 10.0

        # 1m integration (zisoi[2]=0.5 < 1, zisoi[3]=1.0 <= 1)
        for c in 1:nc
            for l in 1:npools
                @test cs.decomp_cpools_1m_col[c, l] ≈ 10.0
            end
        end

        # totlitc_1m (nlevdecomp > 1)
        @test cs.totlitc_1m_col[1] ≈ 10.0
        @test cs.totsomc_1m_col[1] ≈ 10.0

        # decomp_soilc_vr = sum of soil pools per level
        @test cs.decomp_soilc_vr_col[1, 1] ≈ 10.0  # only pool 2 is soil

        # totecosysc = cwdc + totmicc + totlitc + totsomc + ecovegc
        # = 0 + 10 + 10 + 10 + 0 = 30
        @test cs.totecosysc_col[1] ≈ 30.0

        # totc = cwdc + totmicc + totlitc + totsomc + ctrunc + totvegc
        # = 0 + 10 + 10 + 10 + 1 + 0 = 31
        @test cs.totc_col[1] ≈ 31.0
    end

    @testset "summary! with CWD pool" begin
        nc = 2
        ng = 1
        nlev = 1
        nlev_full = 1
        npools = 2
        dzsoi = [1.0]

        cs = CLM.SoilBiogeochemCarbonStateData()
        CLM.soil_bgc_carbon_state_init!(cs, nc, ng, nlev_full, npools)

        cs.decomp_cpools_vr_col .= 5.0
        cs.ctrunc_vr_col .= 0.0

        mask = BitVector([true, true])
        is_litter  = BitVector([false, false])
        is_soil    = BitVector([false, false])
        is_microbe = BitVector([false, false])
        is_cwd     = BitVector([true, false])
        totc_p2c   = [100.0, 100.0]
        totvegc    = [50.0, 50.0]

        CLM.soil_bgc_carbon_state_summary!(cs, mask, 1:nc;
                                           nlevdecomp=nlev,
                                           nlevdecomp_full=nlev_full,
                                           ndecomp_pools=npools,
                                           dzsoi_decomp_vals=dzsoi,
                                           is_litter=is_litter,
                                           is_soil=is_soil,
                                           is_microbe=is_microbe,
                                           is_cwd=is_cwd,
                                           totc_p2c_col=totc_p2c,
                                           totvegc_col=totvegc)

        # cwdc = pool 1 (cwd) integrated = 5.0 * 1.0 = 5.0
        @test cs.cwdc_col[1] ≈ 5.0

        # totecosysc = cwdc(5) + totmicc(0) + totlitc(0) + totsomc(0) + ecovegc(50) = 55
        @test cs.totecosysc_col[1] ≈ 55.0

        # totc = cwdc(5) + totmicc(0) + totlitc(0) + totsomc(0) + ctrunc(0) + totvegc(100) = 105
        @test cs.totc_col[1] ≈ 105.0
    end

    @testset "set_totvegcthresh!" begin
        cs = CLM.SoilBiogeochemCarbonStateData()
        CLM.soil_bgc_carbon_state_set_totvegcthresh!(cs, 100.0)
        @test cs.totvegcthresh == 100.0

        @test_throws ErrorException CLM.soil_bgc_carbon_state_set_totvegcthresh!(cs, 0.0)
        @test_throws ErrorException CLM.soil_bgc_carbon_state_set_totvegcthresh!(cs, -1.0)
    end

    @testset "stub functions run without error" begin
        cs = CLM.SoilBiogeochemCarbonStateData()
        CLM.soil_bgc_carbon_state_init!(cs, 5, 2, 3, 2)

        @test CLM.soil_bgc_carbon_state_init_history!(cs, 1:5) === nothing
        @test CLM.soil_bgc_carbon_state_restart!(cs, 1:5) === nothing
        @test CLM.soil_bgc_carbon_state_dynamic_col_adjustments!(cs, 1:5) === nothing
    end

    @testset "field mutability" begin
        cs = CLM.SoilBiogeochemCarbonStateData()
        CLM.soil_bgc_carbon_state_init!(cs, 3, 1, 2, 2)

        cs.ctrunc_col[1] = 42.0
        @test cs.ctrunc_col[1] == 42.0

        cs.decomp_cpools_vr_col[1, 1, 1] = 99.0
        @test cs.decomp_cpools_vr_col[1, 1, 1] == 99.0

        cs.totc_grc[1] = 7.5
        @test cs.totc_grc[1] == 7.5
    end

    @testset "re-init overwrites previous state" begin
        cs = CLM.SoilBiogeochemCarbonStateData()
        CLM.soil_bgc_carbon_state_init!(cs, 3, 1, 2, 2)
        cs.ctrunc_col[1] = 999.0

        CLM.soil_bgc_carbon_state_init!(cs, 5, 2, 4, 3)
        @test length(cs.ctrunc_col) == 5
        @test all(isnan, cs.ctrunc_col)
        @test size(cs.decomp_cpools_vr_col) == (5, 4, 3)
        @test length(cs.totc_grc) == 2
    end

end
