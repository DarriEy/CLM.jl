@testset "SoilBiogeochemCarbonFluxData" begin

    @testset "default construction" begin
        cf = CLM.SoilBiogeochemCarbonFluxData()
        @test length(cf.hr_col) == 0
        @test length(cf.somc_fire_col) == 0
        @test size(cf.decomp_cpools_sourcesink_col) == (0, 0, 0)
        @test size(cf.decomp_cascade_hr_col) == (0, 0)
        @test size(cf.hr_vr_col) == (0, 0)
        @test size(cf.cn_col) == (0, 0)
        @test cf.NE_AKallsoilc == 0
    end

    @testset "init! basic allocation" begin
        nc = 5
        nlev_full = 3
        npools = 4
        ntrans = 3

        cf = CLM.SoilBiogeochemCarbonFluxData()
        CLM.soil_bgc_carbon_flux_init!(cf, nc, nlev_full, npools, ntrans)

        # Column-level 1D
        @test length(cf.hr_col) == nc
        @test all(isnan, cf.hr_col)
        @test length(cf.michr_col) == nc
        @test length(cf.cwdhr_col) == nc
        @test length(cf.lithr_col) == nc
        @test length(cf.somhr_col) == nc
        @test length(cf.soilc_change_col) == nc
        @test length(cf.somc_fire_col) == nc
        @test length(cf.som_c_leached_col) == nc

        # Environmental scalars (col × nlev)
        @test size(cf.t_scalar_col) == (nc, nlev_full)
        @test all(isnan, cf.t_scalar_col)
        @test size(cf.w_scalar_col) == (nc, nlev_full)
        @test size(cf.o_scalar_col) == (nc, nlev_full)

        # Nitrif/denitrif
        @test size(cf.phr_vr_col) == (nc, nlev_full)
        @test size(cf.fphr_col) == (nc, nlev_full)

        # HR vertically resolved
        @test size(cf.hr_vr_col) == (nc, nlev_full)

        # 3D: col × nlev × npools
        @test size(cf.decomp_cpools_sourcesink_col) == (nc, nlev_full, npools)
        @test size(cf.decomp_k_col) == (nc, nlev_full, npools)
        @test size(cf.decomp_cpools_transport_tendency_col) == (nc, nlev_full, npools)

        # 3D: col × nlev × ntrans
        @test size(cf.c_overflow_vr) == (nc, nlev_full, ntrans)
        @test size(cf.decomp_cascade_hr_vr_col) == (nc, nlev_full, ntrans)
        @test size(cf.decomp_cascade_ctransfer_vr_col) == (nc, nlev_full, ntrans)
        @test size(cf.rf_decomp_cascade_col) == (nc, nlev_full, ntrans)
        @test size(cf.pathfrac_decomp_cascade_col) == (nc, nlev_full, ntrans)

        # 2D: col × ntrans
        @test size(cf.decomp_cascade_hr_col) == (nc, ntrans)
        @test size(cf.decomp_cascade_ctransfer_col) == (nc, ntrans)

        # 2D: col × npools
        @test size(cf.cn_col) == (nc, npools)
        @test size(cf.decomp_cpools_leached_col) == (nc, npools)

        # litr_lig_c_to_n initialized to zero
        @test length(cf.litr_lig_c_to_n_col) == nc
        @test all(==(0.0), cf.litr_lig_c_to_n_col)

        # FATES litter flux (use_fates=false by default → size 1)
        @test length(cf.fates_litter_flux) == 1

        # Matrix-CN fields should NOT be allocated beyond defaults
        @test size(cf.tri_ma_vr) == (1, 1)
    end

    @testset "init! with use_fates" begin
        nc = 3
        cf = CLM.SoilBiogeochemCarbonFluxData()
        CLM.soil_bgc_carbon_flux_init!(cf, nc, 2, 2, 2; use_fates=true)
        @test length(cf.fates_litter_flux) == nc
    end

    @testset "init! with use_soil_matrixcn" begin
        nc = 4
        nlev_full = 3
        nlev = 2
        npools = 3
        ntrans = 2
        nout = 0
        Ntri = 5

        cf = CLM.SoilBiogeochemCarbonFluxData()
        CLM.soil_bgc_carbon_flux_init!(cf, nc, nlev_full, npools, ntrans;
                                        nlevdecomp=nlev,
                                        ndecomp_cascade_outtransitions=nout,
                                        use_soil_matrixcn=true,
                                        Ntri_setup=Ntri)

        @test size(cf.matrix_decomp_fire_k_col) == (nc, nlev * npools)
        @test all(isnan, cf.matrix_decomp_fire_k_col)
        @test size(cf.tri_ma_vr) == (nc, Ntri)

        NE_expected = (ntrans - nout) * nlev + npools * nlev + Ntri + npools * nlev
        @test cf.NE_AKallsoilc == NE_expected
        @test length(cf.RI_AKallsoilc) == NE_expected
        @test length(cf.CI_AKallsoilc) == NE_expected
        @test all(==(-9999), cf.RI_AKallsoilc)

        Ntrans_diag = (ntrans - nout) * nlev + npools * nlev
        @test length(cf.RI_a) == Ntrans_diag
        @test length(cf.CI_a) == Ntrans_diag
    end

    @testset "set_values! basic" begin
        nc = 4
        nlev_full = 2
        npools = 3
        ntrans = 2

        cf = CLM.SoilBiogeochemCarbonFluxData()
        CLM.soil_bgc_carbon_flux_init!(cf, nc, nlev_full, npools, ntrans)

        mask = BitVector([true, false, true, false])
        CLM.soil_bgc_carbon_flux_set_values!(cf, mask, 0.0)

        # Columns 1,3 should be zero; 2,4 should still be NaN
        @test cf.hr_col[1] == 0.0
        @test isnan(cf.hr_col[2])
        @test cf.hr_col[3] == 0.0
        @test isnan(cf.hr_col[4])

        @test cf.somc_fire_col[1] == 0.0
        @test cf.som_c_leached_col[1] == 0.0
        @test cf.somhr_col[1] == 0.0
        @test cf.lithr_col[1] == 0.0
        @test cf.cwdhr_col[1] == 0.0
        @test cf.michr_col[1] == 0.0
        @test cf.soilc_change_col[1] == 0.0

        # HR vr
        @test cf.hr_vr_col[1, 1] == 0.0
        @test isnan(cf.hr_vr_col[2, 1])
        @test cf.hr_vr_col[3, 2] == 0.0

        # Cascade transitions
        @test cf.decomp_cascade_hr_col[1, 1] == 0.0
        @test isnan(cf.decomp_cascade_hr_col[2, 1])
        @test cf.decomp_cascade_hr_vr_col[1, 1, 1] == 0.0
        @test cf.decomp_cascade_ctransfer_col[1, 1] == 0.0
        @test cf.decomp_cascade_ctransfer_vr_col[1, 1, 1] == 0.0
        @test cf.pathfrac_decomp_cascade_col[1, 1, 1] == 0.0
        @test cf.rf_decomp_cascade_col[1, 1, 1] == 0.0
        @test cf.c_overflow_vr[1, 1, 1] == 0.0

        # Pool fields
        @test cf.cn_col[1, 1] == 0.0
        @test isnan(cf.cn_col[2, 1])
        @test cf.decomp_cpools_leached_col[1, 1] == 0.0
        @test cf.decomp_cpools_transport_tendency_col[1, 1, 1] == 0.0
        @test cf.decomp_cpools_sourcesink_col[1, 1, 1] == 0.0
        @test cf.decomp_k_col[1, 1, 1] == 0.0
    end

    @testset "set_values! with matrixcn" begin
        nc = 3
        nlev_full = 2
        nlev = 2
        npools = 2
        ntrans = 1
        Ntri = 3

        cf = CLM.SoilBiogeochemCarbonFluxData()
        CLM.soil_bgc_carbon_flux_init!(cf, nc, nlev_full, npools, ntrans;
                                        nlevdecomp=nlev,
                                        use_soil_matrixcn=true,
                                        Ntri_setup=Ntri)

        mask = BitVector([true, true, true])
        CLM.soil_bgc_carbon_flux_set_values!(cf, mask, 5.0;
                                              use_soil_matrixcn=true,
                                              nlevdecomp=nlev,
                                              Ntri_setup=Ntri)

        @test cf.matrix_decomp_fire_k_col[1, 1] == 5.0
        @test cf.matrix_decomp_fire_k_col[2, nlev * npools] == 5.0
        @test cf.tri_ma_vr[1, 1] == 5.0
        @test cf.tri_ma_vr[3, Ntri] == 5.0
    end

    @testset "init_cold! sets special columns to zero" begin
        nc = 4
        cf = CLM.SoilBiogeochemCarbonFluxData()
        CLM.soil_bgc_carbon_flux_init!(cf, nc, 2, 2, 2)

        mask_special = BitVector([false, true, false, true])
        CLM.soil_bgc_carbon_flux_init_cold!(cf, 1:nc; mask_special=mask_special)

        @test cf.hr_col[2] == 0.0
        @test cf.hr_col[4] == 0.0
        @test isnan(cf.hr_col[1])
        @test isnan(cf.hr_col[3])
    end

    @testset "summary! vertically integrates HR" begin
        nc = 3
        nlev = 2
        npools = 2
        ntrans = 2
        dzsoi = [0.5, 0.5]

        cf = CLM.SoilBiogeochemCarbonFluxData()
        CLM.soil_bgc_carbon_flux_init!(cf, nc, nlev, npools, ntrans)

        mask = BitVector([true, true, true])
        CLM.soil_bgc_carbon_flux_set_values!(cf, mask, 0.0)

        # Set known vertically-resolved HR values
        # transition 1: 2.0 gC/m3/s at each level
        # transition 2: 3.0 gC/m3/s at each level
        for c in 1:nc
            for j in 1:nlev
                cf.decomp_cascade_hr_vr_col[c, j, 1] = 2.0
                cf.decomp_cascade_hr_vr_col[c, j, 2] = 3.0
                cf.decomp_cascade_ctransfer_vr_col[c, j, 1] = 1.0
                cf.decomp_cascade_ctransfer_vr_col[c, j, 2] = 0.5
            end
        end

        # Set transport tendency
        for c in 1:nc
            for j in 1:nlev
                for l in 1:npools
                    cf.decomp_cpools_transport_tendency_col[c, j, l] = -0.1
                end
            end
        end

        # Pool 1 is litter, pool 2 is soil
        # Transition 1 from litter (pool 1), transition 2 from soil (pool 2)
        cascade_donor = [1, 2]
        is_litter  = BitVector([true, false])
        is_soil    = BitVector([false, true])
        is_cwd     = BitVector([false, false])
        is_microbe = BitVector([false, false])

        CLM.soil_bgc_carbon_flux_summary!(cf, mask, 1:nc;
                                           ndecomp_cascade_transitions=ntrans,
                                           nlevdecomp=nlev,
                                           ndecomp_pools=npools,
                                           dzsoi_decomp_vals=dzsoi,
                                           cascade_donor_pool=cascade_donor,
                                           is_soil=is_soil,
                                           is_litter=is_litter,
                                           is_cwd=is_cwd,
                                           is_microbe=is_microbe)

        for c in 1:nc
            # decomp_cascade_hr_col integrated: 2.0 * (0.5 + 0.5) = 2.0 for trans 1
            @test cf.decomp_cascade_hr_col[c, 1] ≈ 2.0
            # 3.0 * (0.5 + 0.5) = 3.0 for trans 2
            @test cf.decomp_cascade_hr_col[c, 2] ≈ 3.0

            # decomp_cascade_ctransfer_col: 1.0 * 1.0 = 1.0 for trans 1
            @test cf.decomp_cascade_ctransfer_col[c, 1] ≈ 1.0
            @test cf.decomp_cascade_ctransfer_col[c, 2] ≈ 0.5

            # hr_vr = sum of all transition HR at each level
            @test cf.hr_vr_col[c, 1] ≈ 5.0  # 2.0 + 3.0
            @test cf.hr_vr_col[c, 2] ≈ 5.0

            # lithr = transition 1 integrated HR = 2.0
            @test cf.lithr_col[c] ≈ 2.0

            # somhr = transition 2 integrated HR = 3.0
            @test cf.somhr_col[c] ≈ 3.0

            # cwdhr = 0 (no cwd pools)
            @test cf.cwdhr_col[c] ≈ 0.0

            # michr = 0 (no microbe pools)
            @test cf.michr_col[c] ≈ 0.0

            # total HR = lithr + somhr + cwdhr + michr = 2.0 + 3.0 = 5.0
            @test cf.hr_col[c] ≈ 5.0

            # SOM C leached: -0.1 * (0.5 + 0.5) = -0.1 per pool, × 2 pools = -0.2
            @test cf.som_c_leached_col[c] ≈ -0.2

            # Per-pool leaching
            @test cf.decomp_cpools_leached_col[c, 1] ≈ -0.1
            @test cf.decomp_cpools_leached_col[c, 2] ≈ -0.1
        end
    end

    @testset "summary! CWD and microbe HR" begin
        nc = 2
        nlev = 1
        npools = 3
        ntrans = 3  # one transition per pool
        dzsoi = [1.0]

        cf = CLM.SoilBiogeochemCarbonFluxData()
        CLM.soil_bgc_carbon_flux_init!(cf, nc, nlev, npools, ntrans)

        mask = BitVector([true, true])
        CLM.soil_bgc_carbon_flux_set_values!(cf, mask, 0.0)

        # Set HR values per transition
        for c in 1:nc
            cf.decomp_cascade_hr_vr_col[c, 1, 1] = 1.0  # from CWD
            cf.decomp_cascade_hr_vr_col[c, 1, 2] = 2.0  # from microbe
            cf.decomp_cascade_hr_vr_col[c, 1, 3] = 3.0  # from litter
            cf.decomp_cpools_transport_tendency_col[c, :, :] .= 0.0
        end

        # Pool 1 = CWD, pool 2 = microbe, pool 3 = litter
        cascade_donor = [1, 2, 3]
        is_cwd     = BitVector([true, false, false])
        is_microbe = BitVector([false, true, false])
        is_litter  = BitVector([false, false, true])
        is_soil    = BitVector([false, false, false])

        CLM.soil_bgc_carbon_flux_summary!(cf, mask, 1:nc;
                                           ndecomp_cascade_transitions=ntrans,
                                           nlevdecomp=nlev,
                                           ndecomp_pools=npools,
                                           dzsoi_decomp_vals=dzsoi,
                                           cascade_donor_pool=cascade_donor,
                                           is_soil=is_soil,
                                           is_litter=is_litter,
                                           is_cwd=is_cwd,
                                           is_microbe=is_microbe)

        for c in 1:nc
            @test cf.cwdhr_col[c] ≈ 1.0
            @test cf.michr_col[c] ≈ 2.0
            @test cf.lithr_col[c] ≈ 3.0
            @test cf.somhr_col[c] ≈ 0.0
            @test cf.hr_col[c] ≈ 6.0
        end
    end

    @testset "stub functions run without error" begin
        cf = CLM.SoilBiogeochemCarbonFluxData()
        CLM.soil_bgc_carbon_flux_init!(cf, 5, 3, 2, 2)

        @test CLM.soil_bgc_carbon_flux_init_history!(cf, 1:5) === nothing
        @test CLM.soil_bgc_carbon_flux_restart!(cf, 1:5) === nothing
    end

    @testset "field mutability" begin
        cf = CLM.SoilBiogeochemCarbonFluxData()
        CLM.soil_bgc_carbon_flux_init!(cf, 3, 2, 2, 2)

        cf.hr_col[1] = 42.0
        @test cf.hr_col[1] == 42.0

        cf.decomp_cascade_hr_vr_col[1, 1, 1] = 99.0
        @test cf.decomp_cascade_hr_vr_col[1, 1, 1] == 99.0

        cf.litr_lig_c_to_n_col[1] = 7.5
        @test cf.litr_lig_c_to_n_col[1] == 7.5
    end

    @testset "re-init overwrites previous state" begin
        cf = CLM.SoilBiogeochemCarbonFluxData()
        CLM.soil_bgc_carbon_flux_init!(cf, 3, 2, 2, 2)
        cf.hr_col[1] = 999.0

        CLM.soil_bgc_carbon_flux_init!(cf, 5, 4, 3, 3)
        @test length(cf.hr_col) == 5
        @test all(isnan, cf.hr_col)
        @test size(cf.decomp_cascade_hr_vr_col) == (5, 4, 3)
        @test size(cf.cn_col) == (5, 3)
    end

end
