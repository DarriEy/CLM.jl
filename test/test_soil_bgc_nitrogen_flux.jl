@testset "SoilBiogeochemNitrogenFluxData" begin

    @testset "default construction" begin
        nf = CLM.SoilBiogeochemNitrogenFluxData()
        @test length(nf.ndep_to_sminn_col) == 0
        @test length(nf.denit_col) == 0
        @test size(nf.decomp_cascade_ntransfer_vr_col) == (0, 0, 0)
        @test size(nf.decomp_cascade_ntransfer_col) == (0, 0)
        @test size(nf.potential_immob_vr_col) == (0, 0)
        @test nf.NE_AKallsoiln == 0
    end

    @testset "init! basic allocation" begin
        nc = 5
        nlev_full = 3
        npools = 4
        ntrans = 3

        nf = CLM.SoilBiogeochemNitrogenFluxData()
        CLM.soil_bgc_nitrogen_flux_init!(nf, nc, nlev_full, npools, ntrans)

        # Deposition 1D
        @test length(nf.ndep_to_sminn_col) == nc
        @test all(isnan, nf.ndep_to_sminn_col)
        @test length(nf.nfix_to_sminn_col) == nc
        @test length(nf.ffix_to_sminn_col) == nc
        @test length(nf.fert_to_sminn_col) == nc
        @test length(nf.soyfixn_to_sminn_col) == nc

        # Column-level 1D summary
        @test length(nf.sminn_to_plant_col) == nc
        @test length(nf.potential_immob_col) == nc
        @test length(nf.actual_immob_col) == nc
        @test length(nf.denit_col) == nc
        @test length(nf.ninputs_col) == nc
        @test length(nf.noutputs_col) == nc
        @test length(nf.som_n_leached_col) == nc
        @test length(nf.sminn_to_plant_fun_col) == nc

        # Vertically-resolved (col × nlev)
        @test size(nf.potential_immob_vr_col) == (nc, nlev_full)
        @test all(isnan, nf.potential_immob_vr_col)
        @test size(nf.actual_immob_vr_col) == (nc, nlev_full)
        @test size(nf.sminn_to_plant_vr_col) == (nc, nlev_full)
        @test size(nf.supplement_to_sminn_vr_col) == (nc, nlev_full)
        @test size(nf.gross_nmin_vr_col) == (nc, nlev_full)
        @test size(nf.net_nmin_vr_col) == (nc, nlev_full)

        # Nitrif/denitrif 2D
        @test size(nf.f_nit_vr_col) == (nc, nlev_full)
        @test size(nf.f_denit_vr_col) == (nc, nlev_full)
        @test size(nf.n2_n2o_ratio_denit_vr_col) == (nc, nlev_full)
        @test size(nf.smin_no3_massdens_vr_col) == (nc, nlev_full)
        @test size(nf.diffus_col) == (nc, nlev_full)
        @test size(nf.r_psi_col) == (nc, nlev_full)
        @test size(nf.anaerobic_frac_col) == (nc, nlev_full)

        # Nitrif/denitrif 1D
        @test length(nf.f_nit_col) == nc
        @test length(nf.f_denit_col) == nc
        @test length(nf.smin_no3_leached_col) == nc
        @test length(nf.smin_no3_runoff_col) == nc

        # Cascade 3D: col × nlev × ntrans
        @test size(nf.decomp_cascade_ntransfer_vr_col) == (nc, nlev_full, ntrans)
        @test size(nf.decomp_cascade_sminn_flux_vr_col) == (nc, nlev_full, ntrans)
        @test size(nf.sminn_to_denit_decomp_cascade_vr_col) == (nc, nlev_full, ntrans)

        # Cascade 2D: col × ntrans
        @test size(nf.decomp_cascade_ntransfer_col) == (nc, ntrans)
        @test size(nf.decomp_cascade_sminn_flux_col) == (nc, ntrans)
        @test size(nf.sminn_to_denit_decomp_cascade_col) == (nc, ntrans)

        # Non-nitrif_denitrif
        @test size(nf.sminn_to_denit_excess_vr_col) == (nc, nlev_full)
        @test length(nf.sminn_to_denit_excess_col) == nc
        @test size(nf.sminn_leached_vr_col) == (nc, nlev_full)
        @test length(nf.sminn_leached_col) == nc

        # Pool diagnostics
        @test size(nf.decomp_npools_leached_col) == (nc, npools)
        @test size(nf.decomp_npools_transport_tendency_col) == (nc, nlev_full, npools)
        @test size(nf.decomp_npools_sourcesink_col) == (nc, nlev_full, npools)

        # FATES litter flux (use_fates=false by default → size 1)
        @test length(nf.fates_litter_flux) == 1

        # Matrix-CN fields should NOT be allocated
        @test nf.NE_AKallsoiln == 0
    end

    @testset "init! with use_fates" begin
        nc = 3
        nf = CLM.SoilBiogeochemNitrogenFluxData()
        CLM.soil_bgc_nitrogen_flux_init!(nf, nc, 2, 2, 2; use_fates=true)
        @test length(nf.fates_litter_flux) == nc
    end

    @testset "init! with use_soil_matrixcn" begin
        nc = 4
        nlev_full = 3
        nlev = 2
        npools = 3
        ntrans = 2
        nout = 0
        Ntri = 5
        ndecomp_pools_vr = npools * nlev

        nf = CLM.SoilBiogeochemNitrogenFluxData()
        CLM.soil_bgc_nitrogen_flux_init!(nf, nc, nlev_full, npools, ntrans;
                                          nlevdecomp=nlev,
                                          ndecomp_cascade_outtransitions=nout,
                                          ndecomp_pools_vr=ndecomp_pools_vr,
                                          use_soil_matrixcn=true,
                                          Ntri_setup=Ntri)

        # NE_AKallsoiln = (Ntrans + nlev*npools) + (Ntrans + Ntri + nlev) + (npools*nlev)
        Ntrans = (ntrans - nout) * nlev
        NE_expected = (Ntrans + nlev * npools) + (Ntrans + Ntri + nlev) + (npools * nlev)
        @test nf.NE_AKallsoiln == NE_expected
        @test length(nf.RI_AKallsoiln) == NE_expected
        @test length(nf.CI_AKallsoiln) == NE_expected
        @test all(==(-9999), nf.RI_AKallsoiln)

        Ntrans_diag = Ntrans + ndecomp_pools_vr
        @test length(nf.RI_na) == Ntrans_diag
        @test length(nf.CI_na) == Ntrans_diag
    end

    @testset "set_values! without nitrif_denitrif" begin
        nc = 4
        nlev_full = 2
        npools = 3
        ntrans = 2

        nf = CLM.SoilBiogeochemNitrogenFluxData()
        CLM.soil_bgc_nitrogen_flux_init!(nf, nc, nlev_full, npools, ntrans)

        mask = BitVector([true, false, true, false])
        CLM.soil_bgc_nitrogen_flux_set_values!(nf, mask, 0.0; use_nitrif_denitrif=false)

        # Columns 1,3 should be zero; 2,4 should still be NaN
        @test nf.ndep_to_sminn_col[1] == 0.0
        @test isnan(nf.ndep_to_sminn_col[2])
        @test nf.ndep_to_sminn_col[3] == 0.0
        @test isnan(nf.ndep_to_sminn_col[4])

        @test nf.denit_col[1] == 0.0
        @test nf.potential_immob_col[1] == 0.0
        @test nf.actual_immob_col[1] == 0.0
        @test nf.sminn_to_plant_col[1] == 0.0
        @test nf.ninputs_col[1] == 0.0
        @test nf.noutputs_col[1] == 0.0
        @test nf.som_n_leached_col[1] == 0.0

        # Non-nitrif_denitrif path sets these
        @test nf.sminn_to_denit_excess_col[1] == 0.0
        @test nf.sminn_leached_col[1] == 0.0
        @test nf.sminn_to_denit_excess_vr_col[1, 1] == 0.0
        @test nf.sminn_leached_vr_col[1, 1] == 0.0

        # VR fields
        @test nf.potential_immob_vr_col[1, 1] == 0.0
        @test isnan(nf.potential_immob_vr_col[2, 1])
        @test nf.actual_immob_vr_col[1, 1] == 0.0

        # Cascade fields
        @test nf.decomp_cascade_ntransfer_col[1, 1] == 0.0
        @test nf.decomp_cascade_ntransfer_vr_col[1, 1, 1] == 0.0
        @test nf.decomp_cascade_sminn_flux_col[1, 1] == 0.0
        @test nf.sminn_to_denit_decomp_cascade_col[1, 1] == 0.0
        @test nf.sminn_to_denit_decomp_cascade_vr_col[1, 1, 1] == 0.0

        # Pool fields
        @test nf.decomp_npools_leached_col[1, 1] == 0.0
        @test nf.decomp_npools_transport_tendency_col[1, 1, 1] == 0.0
        @test nf.decomp_npools_sourcesink_col[1, 1, 1] == 0.0
    end

    @testset "set_values! with nitrif_denitrif" begin
        nc = 3
        nlev_full = 2
        npools = 2
        ntrans = 2

        nf = CLM.SoilBiogeochemNitrogenFluxData()
        CLM.soil_bgc_nitrogen_flux_init!(nf, nc, nlev_full, npools, ntrans)

        mask = BitVector([true, true, true])
        CLM.soil_bgc_nitrogen_flux_set_values!(nf, mask, 5.0; use_nitrif_denitrif=true)

        # Nitrif/denitrif path sets these
        @test nf.f_nit_col[1] == 5.0
        @test nf.f_denit_col[1] == 5.0
        @test nf.pot_f_nit_col[1] == 5.0
        @test nf.pot_f_denit_col[1] == 5.0
        @test nf.f_n2o_denit_col[1] == 5.0
        @test nf.f_n2o_nit_col[1] == 5.0
        @test nf.smin_no3_leached_col[1] == 5.0
        @test nf.smin_no3_runoff_col[1] == 5.0

        # VR nitrif/denitrif fields
        @test nf.f_nit_vr_col[1, 1] == 5.0
        @test nf.f_denit_vr_col[1, 1] == 5.0
        @test nf.smin_no3_leached_vr_col[1, 1] == 5.0
        @test nf.diffus_col[1, 1] == 5.0
        @test nf.r_psi_col[1, 1] == 5.0
        @test nf.anaerobic_frac_col[1, 1] == 5.0

        # Common fields
        @test nf.potential_immob_vr_col[1, 1] == 5.0
        @test nf.ndep_to_sminn_col[1] == 5.0

        # Non-nitrif_denitrif should NOT be set (still NaN from init)
        # sminn_to_denit_excess_col is not set in use_nitrif_denitrif=true path
        @test isnan(nf.sminn_to_denit_excess_col[1])
        @test isnan(nf.sminn_leached_col[1])
    end

    @testset "init_cold! sets special columns to zero" begin
        nc = 4
        nf = CLM.SoilBiogeochemNitrogenFluxData()
        CLM.soil_bgc_nitrogen_flux_init!(nf, nc, 2, 2, 2)

        mask_special = BitVector([false, true, false, true])
        CLM.soil_bgc_nitrogen_flux_init_cold!(nf, 1:nc; mask_special=mask_special)

        @test nf.ndep_to_sminn_col[2] == 0.0
        @test nf.ndep_to_sminn_col[4] == 0.0
        @test isnan(nf.ndep_to_sminn_col[1])
        @test isnan(nf.ndep_to_sminn_col[3])
    end

    @testset "summary! without nitrif_denitrif" begin
        nc = 3
        nlev = 2
        npools = 2
        ntrans = 2
        dzsoi = [0.5, 0.5]

        nf = CLM.SoilBiogeochemNitrogenFluxData()
        CLM.soil_bgc_nitrogen_flux_init!(nf, nc, nlev, npools, ntrans)

        mask = BitVector([true, true, true])
        CLM.soil_bgc_nitrogen_flux_set_values!(nf, mask, 0.0; use_nitrif_denitrif=false)

        # Set known vertically-resolved values
        for c in 1:nc
            for j in 1:nlev
                nf.decomp_cascade_ntransfer_vr_col[c, j, 1] = 2.0
                nf.decomp_cascade_ntransfer_vr_col[c, j, 2] = 3.0
                nf.decomp_cascade_sminn_flux_vr_col[c, j, 1] = 1.0
                nf.decomp_cascade_sminn_flux_vr_col[c, j, 2] = 0.5
                nf.sminn_to_denit_decomp_cascade_vr_col[c, j, 1] = 0.1
                nf.sminn_to_denit_decomp_cascade_vr_col[c, j, 2] = 0.2
                nf.sminn_to_denit_excess_vr_col[c, j] = 0.05
                nf.sminn_leached_vr_col[c, j] = 0.03
                nf.supplement_to_sminn_vr_col[c, j] = 0.4
                nf.decomp_npools_transport_tendency_col[c, j, 1] = -0.1
                nf.decomp_npools_transport_tendency_col[c, j, 2] = -0.2
            end
        end

        CLM.soil_bgc_nitrogen_flux_summary!(nf, mask, 1:nc;
                                             ndecomp_cascade_transitions=ntrans,
                                             nlevdecomp=nlev,
                                             ndecomp_pools=npools,
                                             dzsoi_decomp_vals=dzsoi,
                                             use_nitrif_denitrif=false)

        for c in 1:nc
            # decomp_cascade_ntransfer_col: 2.0 * (0.5 + 0.5) = 2.0 for trans 1
            @test nf.decomp_cascade_ntransfer_col[c, 1] ≈ 2.0
            @test nf.decomp_cascade_ntransfer_col[c, 2] ≈ 3.0

            @test nf.decomp_cascade_sminn_flux_col[c, 1] ≈ 1.0
            @test nf.decomp_cascade_sminn_flux_col[c, 2] ≈ 0.5

            # Denitrification cascade: 0.1 * 1.0 = 0.1 for trans 1
            @test nf.sminn_to_denit_decomp_cascade_col[c, 1] ≈ 0.1
            @test nf.sminn_to_denit_decomp_cascade_col[c, 2] ≈ 0.2

            # Excess denitrification: 0.05 * 1.0 = 0.05
            @test nf.sminn_to_denit_excess_col[c] ≈ 0.05

            # Leaching: 0.03 * 1.0 = 0.03
            @test nf.sminn_leached_col[c] ≈ 0.03

            # Total denitrification = sum(cascade) + excess = 0.1 + 0.2 + 0.05 = 0.35
            @test nf.denit_col[c] ≈ 0.35

            # Supplemental N: 0.4 * 1.0 = 0.4
            @test nf.supplement_to_sminn_col[c] ≈ 0.4

            # SOM N leached: pool1 = -0.1*1.0 = -0.1, pool2 = -0.2*1.0 = -0.2
            @test nf.decomp_npools_leached_col[c, 1] ≈ -0.1
            @test nf.decomp_npools_leached_col[c, 2] ≈ -0.2
            @test nf.som_n_leached_col[c] ≈ -0.3
        end
    end

    @testset "summary! with nitrif_denitrif" begin
        nc = 2
        nlev = 2
        npools = 2
        ntrans = 1
        dzsoi = [0.5, 0.5]

        nf = CLM.SoilBiogeochemNitrogenFluxData()
        CLM.soil_bgc_nitrogen_flux_init!(nf, nc, nlev, npools, ntrans)

        mask = BitVector([true, true])
        CLM.soil_bgc_nitrogen_flux_set_values!(nf, mask, 0.0; use_nitrif_denitrif=true)

        # Set known values
        for c in 1:nc
            for j in 1:nlev
                nf.decomp_cascade_ntransfer_vr_col[c, j, 1] = 1.0
                nf.decomp_cascade_sminn_flux_vr_col[c, j, 1] = 0.5
                nf.f_nit_vr_col[c, j]            = 2.0
                nf.f_denit_vr_col[c, j]           = 3.0
                nf.pot_f_nit_vr_col[c, j]         = 2.5
                nf.pot_f_denit_vr_col[c, j]       = 3.5
                nf.f_n2o_nit_vr_col[c, j]         = 0.1
                nf.f_n2o_denit_vr_col[c, j]       = 0.2
                nf.smin_no3_leached_vr_col[c, j]  = 0.4
                nf.smin_no3_runoff_vr_col[c, j]   = 0.6
                nf.supplement_to_sminn_vr_col[c, j] = 0.8
                nf.decomp_npools_transport_tendency_col[c, j, 1] = -0.05
                nf.decomp_npools_transport_tendency_col[c, j, 2] = -0.10
            end
        end

        CLM.soil_bgc_nitrogen_flux_summary!(nf, mask, 1:nc;
                                             ndecomp_cascade_transitions=ntrans,
                                             nlevdecomp=nlev,
                                             ndecomp_pools=npools,
                                             dzsoi_decomp_vals=dzsoi,
                                             use_nitrif_denitrif=true)

        for c in 1:nc
            # Cascade ntransfer: 1.0 * 1.0 = 1.0
            @test nf.decomp_cascade_ntransfer_col[c, 1] ≈ 1.0
            @test nf.decomp_cascade_sminn_flux_col[c, 1] ≈ 0.5

            # f_nit: 2.0 * 1.0 = 2.0
            @test nf.f_nit_col[c] ≈ 2.0
            @test nf.f_denit_col[c] ≈ 3.0
            @test nf.pot_f_nit_col[c] ≈ 2.5
            @test nf.pot_f_denit_col[c] ≈ 3.5
            @test nf.f_n2o_nit_col[c] ≈ 0.1
            @test nf.f_n2o_denit_col[c] ≈ 0.2
            @test nf.smin_no3_leached_col[c] ≈ 0.4
            @test nf.smin_no3_runoff_col[c] ≈ 0.6

            # denit = f_denit in nitrif_denitrif mode
            @test nf.denit_col[c] ≈ 3.0

            # supplement: 0.8 * 1.0 = 0.8
            @test nf.supplement_to_sminn_col[c] ≈ 0.8

            # SOM N leached
            @test nf.decomp_npools_leached_col[c, 1] ≈ -0.05
            @test nf.decomp_npools_leached_col[c, 2] ≈ -0.1
            @test nf.som_n_leached_col[c] ≈ -0.15
        end
    end

    @testset "stub functions run without error" begin
        nf = CLM.SoilBiogeochemNitrogenFluxData()
        CLM.soil_bgc_nitrogen_flux_init!(nf, 5, 3, 2, 2)

        @test CLM.soil_bgc_nitrogen_flux_init_history!(nf, 1:5) === nothing
    end

    @testset "field mutability" begin
        nf = CLM.SoilBiogeochemNitrogenFluxData()
        CLM.soil_bgc_nitrogen_flux_init!(nf, 3, 2, 2, 2)

        nf.ndep_to_sminn_col[1] = 42.0
        @test nf.ndep_to_sminn_col[1] == 42.0

        nf.decomp_cascade_ntransfer_vr_col[1, 1, 1] = 99.0
        @test nf.decomp_cascade_ntransfer_vr_col[1, 1, 1] == 99.0

        nf.denit_col[1] = 7.5
        @test nf.denit_col[1] == 7.5
    end

    @testset "re-init overwrites previous state" begin
        nf = CLM.SoilBiogeochemNitrogenFluxData()
        CLM.soil_bgc_nitrogen_flux_init!(nf, 3, 2, 2, 2)
        nf.ndep_to_sminn_col[1] = 999.0

        CLM.soil_bgc_nitrogen_flux_init!(nf, 5, 4, 3, 3)
        @test length(nf.ndep_to_sminn_col) == 5
        @test all(isnan, nf.ndep_to_sminn_col)
        @test size(nf.decomp_cascade_ntransfer_vr_col) == (5, 4, 3)
        @test size(nf.decomp_npools_leached_col) == (5, 3)
    end

end
