@testset "SoilBiogeochemDecomp" begin

    # ----------------------------------------------------------------
    # Helper: create minimal data structures for decomposition tests
    # ----------------------------------------------------------------
    function make_decomp_data(; nc=4, nlevdecomp=3,
                               ndecomp_pools=7,
                               ndecomp_cascade_transitions=10)
        bounds = 1:nc
        nlevdecomp_full = nlevdecomp

        # --- Carbon state ---
        cs = CLM.SoilBiogeochemCarbonStateData()
        CLM.soil_bgc_carbon_state_init!(cs, nc, 1, nlevdecomp_full, ndecomp_pools)
        # Set decomp C pools to realistic values
        for c in 1:nc, j in 1:nlevdecomp, l in 1:ndecomp_pools
            cs.decomp_cpools_vr_col[c, j, l] = 100.0 + 10.0 * l
        end

        # --- Nitrogen state ---
        ns = CLM.SoilBiogeochemNitrogenStateData()
        CLM.soil_bgc_nitrogen_state_init!(ns, nc, 1, nlevdecomp_full, ndecomp_pools)
        # Set decomp N pools based on C:N ratios
        for c in 1:nc, j in 1:nlevdecomp, l in 1:ndecomp_pools
            ns.decomp_npools_vr_col[c, j, l] = cs.decomp_cpools_vr_col[c, j, l] / 12.0
        end

        # --- Carbon flux ---
        cf = CLM.SoilBiogeochemCarbonFluxData()
        CLM.soil_bgc_carbon_flux_init!(cf, nc, nlevdecomp_full, ndecomp_pools,
                                        ndecomp_cascade_transitions)
        # Set rf_decomp_cascade to typical values
        for c in 1:nc, j in 1:nlevdecomp, k in 1:ndecomp_cascade_transitions
            cf.rf_decomp_cascade_col[c, j, k] = 0.5
            cf.c_overflow_vr[c, j, k] = 0.0
        end
        # Set scalars
        for c in 1:nc, j in 1:nlevdecomp
            cf.w_scalar_col[c, j] = 0.8
            cf.phr_vr_col[c, j] = 1.0e-5
        end

        # --- Nitrogen flux ---
        nf = CLM.SoilBiogeochemNitrogenFluxData()
        CLM.soil_bgc_nitrogen_flux_init!(nf, nc, nlevdecomp_full, ndecomp_pools,
                                          ndecomp_cascade_transitions)
        # Initialize net/gross mineralization to zero
        for c in 1:nc
            nf.net_nmin_col[c] = 0.0
            nf.gross_nmin_col[c] = 0.0
            for j in 1:nlevdecomp
                nf.net_nmin_vr_col[c, j] = 0.0
                nf.gross_nmin_vr_col[c, j] = 0.0
            end
        end

        # --- State ---
        st = CLM.SoilBiogeochemStateData()
        CLM.soil_bgc_state_init!(st, nc, nc, nlevdecomp_full, ndecomp_cascade_transitions)
        # Set fpi_vr to 0.5 (partial immobilization)
        for c in 1:nc, j in 1:nlevdecomp
            st.fpi_vr_col[c, j] = 0.5
        end

        # --- BGC cascade configuration ---
        # Use BGC cascade for testing
        params_bgc = CLM.DecompBGCParams(
            bgc_initial_Cstocks = fill(200.0, ndecomp_pools),
        )
        cn_params = CLM.CNSharedParamsData()
        bgc_state = CLM.DecompBGCState()
        cascade_con = CLM.DecompCascadeConData()

        cellsand = fill(50.0, nc, max(nlevdecomp, 5))
        CLM.init_decomp_cascade_bgc!(
            bgc_state, cascade_con, params_bgc, cn_params;
            cellsand=cellsand,
            bounds=bounds,
            nlevdecomp=nlevdecomp,
            ndecomp_pools_max=ndecomp_pools,
            ndecomp_cascade_transitions_max=ndecomp_cascade_transitions,
            use_fates=false)

        # --- Decomp parameters ---
        decomp_params = CLM.DecompParams(dnp=0.01)

        # --- Mask ---
        mask_bgc_soilc = trues(nc)

        # --- Work arrays ---
        cn_decomp_pools = zeros(nc, nlevdecomp, ndecomp_pools)
        p_decomp_cpool_loss = fill(1.0e-6, nc, nlevdecomp, ndecomp_cascade_transitions)
        pmnf_decomp_cascade = fill(-1.0e-7, nc, nlevdecomp, ndecomp_cascade_transitions)
        p_decomp_npool_to_din = zeros(nc, nlevdecomp, ndecomp_cascade_transitions)

        # dzsoi_decomp
        dzsoi_decomp = fill(0.1, nlevdecomp)

        return (cf=cf, cs=cs, nf=nf, ns=ns, st=st, cascade_con=cascade_con,
                decomp_params=decomp_params, mask_bgc_soilc=mask_bgc_soilc,
                bounds=bounds, nlevdecomp=nlevdecomp,
                ndecomp_pools=ndecomp_pools,
                ndecomp_cascade_transitions=ndecomp_cascade_transitions,
                cn_decomp_pools=cn_decomp_pools,
                p_decomp_cpool_loss=p_decomp_cpool_loss,
                pmnf_decomp_cascade=pmnf_decomp_cascade,
                p_decomp_npool_to_din=p_decomp_npool_to_din,
                dzsoi_decomp=dzsoi_decomp,
                nc=nc)
    end

    # ================================================================
    # Test 1: DecompParams construction and read
    # ================================================================
    @testset "DecompParams construction" begin
        params = CLM.DecompParams()
        @test params.dnp == 0.01

        params2 = CLM.DecompParams()
        CLM.decomp_read_params!(params2; dnp=0.05)
        @test params2.dnp == 0.05
    end

    # ================================================================
    # Test 2: I_ATM constant
    # ================================================================
    @testset "I_ATM constant" begin
        @test CLM.I_ATM == 0
    end

    # ================================================================
    # Test 3: Basic decomposition run (no errors)
    # ================================================================
    @testset "Basic decomposition run" begin
        d = make_decomp_data()

        CLM.soil_biogeochem_decomp!(
            d.cf, d.cs, d.nf, d.ns, d.st, d.cascade_con, d.decomp_params;
            mask_bgc_soilc=d.mask_bgc_soilc,
            bounds=d.bounds,
            nlevdecomp=d.nlevdecomp,
            ndecomp_pools=d.ndecomp_pools,
            ndecomp_cascade_transitions=d.ndecomp_cascade_transitions,
            cn_decomp_pools=d.cn_decomp_pools,
            p_decomp_cpool_loss=d.p_decomp_cpool_loss,
            pmnf_decomp_cascade=d.pmnf_decomp_cascade,
            p_decomp_npool_to_din=d.p_decomp_npool_to_din,
            dzsoi_decomp=d.dzsoi_decomp)

        # cn_decomp_pools should be computed for all pools
        for l in 1:d.ndecomp_pools
            for j in 1:d.nlevdecomp
                @test d.cn_decomp_pools[1, j, l] > 0.0
            end
        end

        # HR and ctransfer should be non-negative
        for k in 1:d.ndecomp_cascade_transitions
            for j in 1:d.nlevdecomp
                for c in 1:d.nc
                    @test d.cf.decomp_cascade_hr_vr_col[c, j, k] >= 0.0
                    @test d.cf.decomp_cascade_ctransfer_vr_col[c, j, k] >= 0.0
                end
            end
        end
    end

    # ================================================================
    # Test 4: C:N ratio computation for floating pools
    # ================================================================
    @testset "C:N ratio computation" begin
        d = make_decomp_data()

        CLM.soil_biogeochem_decomp!(
            d.cf, d.cs, d.nf, d.ns, d.st, d.cascade_con, d.decomp_params;
            mask_bgc_soilc=d.mask_bgc_soilc,
            bounds=d.bounds,
            nlevdecomp=d.nlevdecomp,
            ndecomp_pools=d.ndecomp_pools,
            ndecomp_cascade_transitions=d.ndecomp_cascade_transitions,
            cn_decomp_pools=d.cn_decomp_pools,
            p_decomp_cpool_loss=d.p_decomp_cpool_loss,
            pmnf_decomp_cascade=d.pmnf_decomp_cascade,
            p_decomp_npool_to_din=d.p_decomp_npool_to_din,
            dzsoi_decomp=d.dzsoi_decomp)

        # Floating pools (litter + CWD) should have C:N = C/N from pool values
        # Pools 1-3 and 7 are floating in BGC cascade
        for l in [1, 2, 3, 7]
            if d.cascade_con.floating_cn_ratio_decomp_pools[l]
                @test d.cn_decomp_pools[1, 1, l] ≈ d.cs.decomp_cpools_vr_col[1, 1, l] / d.ns.decomp_npools_vr_col[1, 1, l]
            end
        end

        # Fixed pools (SOM 4,5,6) should have initial_cn_ratio
        for l in [4, 5, 6]
            if !d.cascade_con.floating_cn_ratio_decomp_pools[l]
                @test d.cn_decomp_pools[1, 1, l] == d.cascade_con.initial_cn_ratio[l]
            end
        end
    end

    # ================================================================
    # Test 5: fpi_vr scaling of immobilization steps
    # ================================================================
    @testset "fpi_vr scaling" begin
        d = make_decomp_data()

        # Set pmnf > 0 for first transition (immobilization)
        fill!(d.pmnf_decomp_cascade, 1.0e-7)
        p_cpool_loss_before = copy(d.p_decomp_cpool_loss)
        pmnf_before = copy(d.pmnf_decomp_cascade)

        CLM.soil_biogeochem_decomp!(
            d.cf, d.cs, d.nf, d.ns, d.st, d.cascade_con, d.decomp_params;
            mask_bgc_soilc=d.mask_bgc_soilc,
            bounds=d.bounds,
            nlevdecomp=d.nlevdecomp,
            ndecomp_pools=d.ndecomp_pools,
            ndecomp_cascade_transitions=d.ndecomp_cascade_transitions,
            cn_decomp_pools=d.cn_decomp_pools,
            p_decomp_cpool_loss=d.p_decomp_cpool_loss,
            pmnf_decomp_cascade=d.pmnf_decomp_cascade,
            p_decomp_npool_to_din=d.p_decomp_npool_to_din,
            dzsoi_decomp=d.dzsoi_decomp)

        # fpi_vr = 0.5, so immobilization-limited losses should be halved
        fpi = 0.5
        for k in 1:d.ndecomp_cascade_transitions
            donor = d.cascade_con.cascade_donor_pool[k]
            for j in 1:d.nlevdecomp
                if d.cs.decomp_cpools_vr_col[1, j, donor] > 0.0
                    # pmnf was positive → should be scaled by fpi
                    @test d.p_decomp_cpool_loss[1, j, k] ≈ p_cpool_loss_before[1, j, k] * fpi
                    @test d.pmnf_decomp_cascade[1, j, k] ≈ pmnf_before[1, j, k] * fpi
                end
            end
        end
    end

    # ================================================================
    # Test 6: Denitrification with non-nitrif_denitrif mode
    # ================================================================
    @testset "Denitrification (simple mode)" begin
        d = make_decomp_data()

        # Set pmnf < 0 for mineralization (net N release)
        fill!(d.pmnf_decomp_cascade, -1.0e-7)
        dnp = d.decomp_params.dnp

        CLM.soil_biogeochem_decomp!(
            d.cf, d.cs, d.nf, d.ns, d.st, d.cascade_con, d.decomp_params;
            mask_bgc_soilc=d.mask_bgc_soilc,
            bounds=d.bounds,
            nlevdecomp=d.nlevdecomp,
            ndecomp_pools=d.ndecomp_pools,
            ndecomp_cascade_transitions=d.ndecomp_cascade_transitions,
            cn_decomp_pools=d.cn_decomp_pools,
            p_decomp_cpool_loss=d.p_decomp_cpool_loss,
            pmnf_decomp_cascade=d.pmnf_decomp_cascade,
            p_decomp_npool_to_din=d.p_decomp_npool_to_din,
            dzsoi_decomp=d.dzsoi_decomp,
            use_nitrif_denitrif=false)

        # Check denitrification: sminn_to_denit = -dnp * pmnf (for pmnf < 0)
        for k in 1:d.ndecomp_cascade_transitions
            donor = d.cascade_con.cascade_donor_pool[k]
            for j in 1:d.nlevdecomp
                if d.cs.decomp_cpools_vr_col[1, j, donor] > 0.0 && d.pmnf_decomp_cascade[1, j, k] <= 0.0
                    expected = -dnp * d.pmnf_decomp_cascade[1, j, k]
                    @test d.nf.sminn_to_denit_decomp_cascade_vr_col[1, j, k] ≈ expected
                end
            end
        end
    end

    # ================================================================
    # Test 7: HR and C transfer conservation
    # ================================================================
    @testset "HR + ctransfer = cpool_loss" begin
        d = make_decomp_data()
        fill!(d.pmnf_decomp_cascade, -1.0e-7)

        CLM.soil_biogeochem_decomp!(
            d.cf, d.cs, d.nf, d.ns, d.st, d.cascade_con, d.decomp_params;
            mask_bgc_soilc=d.mask_bgc_soilc,
            bounds=d.bounds,
            nlevdecomp=d.nlevdecomp,
            ndecomp_pools=d.ndecomp_pools,
            ndecomp_cascade_transitions=d.ndecomp_cascade_transitions,
            cn_decomp_pools=d.cn_decomp_pools,
            p_decomp_cpool_loss=d.p_decomp_cpool_loss,
            pmnf_decomp_cascade=d.pmnf_decomp_cascade,
            p_decomp_npool_to_din=d.p_decomp_npool_to_din,
            dzsoi_decomp=d.dzsoi_decomp)

        # For non-MIMICS: HR + ctransfer should equal p_decomp_cpool_loss
        for k in 1:d.ndecomp_cascade_transitions
            donor = d.cascade_con.cascade_donor_pool[k]
            for j in 1:d.nlevdecomp
                for c in 1:d.nc
                    if d.cs.decomp_cpools_vr_col[c, j, donor] > 0.0
                        hr = d.cf.decomp_cascade_hr_vr_col[c, j, k]
                        ct = d.cf.decomp_cascade_ctransfer_vr_col[c, j, k]
                        cl = d.p_decomp_cpool_loss[c, j, k]
                        @test hr + ct ≈ cl atol=1.0e-15
                    end
                end
            end
        end
    end

    # ================================================================
    # Test 8: Net mineralization vertical integration
    # ================================================================
    @testset "Net mineralization vertical integration" begin
        d = make_decomp_data()
        fill!(d.pmnf_decomp_cascade, -1.0e-7)

        # Set non-zero gross_nmin_vr for integration test
        for c in 1:d.nc, j in 1:d.nlevdecomp
            d.nf.gross_nmin_vr_col[c, j] = 1.0e-6
        end

        CLM.soil_biogeochem_decomp!(
            d.cf, d.cs, d.nf, d.ns, d.st, d.cascade_con, d.decomp_params;
            mask_bgc_soilc=d.mask_bgc_soilc,
            bounds=d.bounds,
            nlevdecomp=d.nlevdecomp,
            ndecomp_pools=d.ndecomp_pools,
            ndecomp_cascade_transitions=d.ndecomp_cascade_transitions,
            cn_decomp_pools=d.cn_decomp_pools,
            p_decomp_cpool_loss=d.p_decomp_cpool_loss,
            pmnf_decomp_cascade=d.pmnf_decomp_cascade,
            p_decomp_npool_to_din=d.p_decomp_npool_to_din,
            dzsoi_decomp=d.dzsoi_decomp)

        # gross_nmin should be vertically integrated
        expected_gross = sum(d.nf.gross_nmin_vr_col[1, j] * d.dzsoi_decomp[j] for j in 1:d.nlevdecomp)
        @test d.nf.gross_nmin_col[1] ≈ expected_gross

        # net_nmin should also be vertically integrated (non-zero from pmnf contributions)
        @test d.nf.net_nmin_col[1] != 0.0
    end

    # ================================================================
    # Test 9: Methane code (use_lch4) fphr calculation
    # ================================================================
    @testset "LCH4 fphr calculation" begin
        d = make_decomp_data()
        fill!(d.pmnf_decomp_cascade, -1.0e-7)

        CLM.soil_biogeochem_decomp!(
            d.cf, d.cs, d.nf, d.ns, d.st, d.cascade_con, d.decomp_params;
            mask_bgc_soilc=d.mask_bgc_soilc,
            bounds=d.bounds,
            nlevdecomp=d.nlevdecomp,
            ndecomp_pools=d.ndecomp_pools,
            ndecomp_cascade_transitions=d.ndecomp_cascade_transitions,
            cn_decomp_pools=d.cn_decomp_pools,
            p_decomp_cpool_loss=d.p_decomp_cpool_loss,
            pmnf_decomp_cascade=d.pmnf_decomp_cascade,
            p_decomp_npool_to_din=d.p_decomp_npool_to_din,
            dzsoi_decomp=d.dzsoi_decomp,
            use_lch4=true)

        # fphr should be computed and positive
        for c in 1:d.nc
            for j in 1:d.nlevdecomp
                @test d.cf.fphr_col[c, j] >= 0.01
            end
        end
    end

    # ================================================================
    # Test 10: Masked columns are skipped
    # ================================================================
    @testset "Masked columns skipped" begin
        d = make_decomp_data()

        # Mask out columns 2 and 3
        d.mask_bgc_soilc[2] = false
        d.mask_bgc_soilc[3] = false

        # Set sentinel values for masked columns
        for j in 1:d.nlevdecomp, k in 1:d.ndecomp_cascade_transitions
            d.cf.decomp_cascade_hr_vr_col[2, j, k] = -999.0
            d.cf.decomp_cascade_hr_vr_col[3, j, k] = -999.0
        end

        CLM.soil_biogeochem_decomp!(
            d.cf, d.cs, d.nf, d.ns, d.st, d.cascade_con, d.decomp_params;
            mask_bgc_soilc=d.mask_bgc_soilc,
            bounds=d.bounds,
            nlevdecomp=d.nlevdecomp,
            ndecomp_pools=d.ndecomp_pools,
            ndecomp_cascade_transitions=d.ndecomp_cascade_transitions,
            cn_decomp_pools=d.cn_decomp_pools,
            p_decomp_cpool_loss=d.p_decomp_cpool_loss,
            pmnf_decomp_cascade=d.pmnf_decomp_cascade,
            p_decomp_npool_to_din=d.p_decomp_npool_to_din,
            dzsoi_decomp=d.dzsoi_decomp)

        # Masked columns should not be modified
        for j in 1:d.nlevdecomp, k in 1:d.ndecomp_cascade_transitions
            @test d.cf.decomp_cascade_hr_vr_col[2, j, k] == -999.0
            @test d.cf.decomp_cascade_hr_vr_col[3, j, k] == -999.0
        end

        # Unmasked columns should be processed
        for j in 1:d.nlevdecomp, k in 1:d.ndecomp_cascade_transitions
            @test d.cf.decomp_cascade_hr_vr_col[1, j, k] != -999.0
        end
    end

    # ================================================================
    # Test 11: use_nitrif_denitrif mode
    # ================================================================
    @testset "use_nitrif_denitrif mode" begin
        d = make_decomp_data()
        fill!(d.pmnf_decomp_cascade, -1.0e-7)

        # Set sentinel for sminn_to_denit
        for c in 1:d.nc, j in 1:d.nlevdecomp, k in 1:d.ndecomp_cascade_transitions
            d.nf.sminn_to_denit_decomp_cascade_vr_col[c, j, k] = -999.0
        end

        CLM.soil_biogeochem_decomp!(
            d.cf, d.cs, d.nf, d.ns, d.st, d.cascade_con, d.decomp_params;
            mask_bgc_soilc=d.mask_bgc_soilc,
            bounds=d.bounds,
            nlevdecomp=d.nlevdecomp,
            ndecomp_pools=d.ndecomp_pools,
            ndecomp_cascade_transitions=d.ndecomp_cascade_transitions,
            cn_decomp_pools=d.cn_decomp_pools,
            p_decomp_cpool_loss=d.p_decomp_cpool_loss,
            pmnf_decomp_cascade=d.pmnf_decomp_cascade,
            p_decomp_npool_to_din=d.p_decomp_npool_to_din,
            dzsoi_decomp=d.dzsoi_decomp,
            use_nitrif_denitrif=true)

        # With use_nitrif_denitrif=true, sminn_to_denit should NOT be modified
        for c in 1:d.nc, j in 1:d.nlevdecomp, k in 1:d.ndecomp_cascade_transitions
            @test d.nf.sminn_to_denit_decomp_cascade_vr_col[c, j, k] == -999.0
        end
    end

end
