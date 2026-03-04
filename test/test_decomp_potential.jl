@testset "SoilBiogeochemPotentialMod" begin

    # ---------------------------------------------------------------
    # Helper: build minimal data structures for testing
    # ---------------------------------------------------------------
    function make_test_data(;
            nc::Int=3,
            nlevdecomp::Int=2,
            ndecomp_pools::Int=4,
            ndecomp_cascade_transitions::Int=5,
            use_mimics::Bool=false)

        # -- Cascade configuration --
        cascade_con = CLM.DecompCascadeConData()
        cascade_con.cascade_donor_pool    = [1, 2, 3, 3, 4]
        cascade_con.cascade_receiver_pool = [3, 3, 4, CLM.I_ATM, 3]
        # Pools 1,2 have floating C:N (litter-like), pools 3,4 fixed
        cascade_con.floating_cn_ratio_decomp_pools = BitVector([true, true, false, false])
        cascade_con.initial_cn_ratio = [20.0, 25.0, 12.0, 10.0]

        # -- Carbon state --
        cs = CLM.SoilBiogeochemCarbonStateData()
        cs.decomp_cpools_vr_col = zeros(nc, nlevdecomp, ndecomp_pools)
        # Set some initial C pools
        for c in 1:nc, j in 1:nlevdecomp
            cs.decomp_cpools_vr_col[c, j, 1] = 100.0 + 10.0 * c  # litter 1
            cs.decomp_cpools_vr_col[c, j, 2] = 80.0 + 5.0 * c    # litter 2
            cs.decomp_cpools_vr_col[c, j, 3] = 50.0 + 3.0 * c    # SOM 1
            cs.decomp_cpools_vr_col[c, j, 4] = 30.0 + 2.0 * c    # SOM 2
        end

        # -- Nitrogen state --
        ns = CLM.SoilBiogeochemNitrogenStateData()
        ns.decomp_npools_vr_col = zeros(nc, nlevdecomp, ndecomp_pools)
        for c in 1:nc, j in 1:nlevdecomp
            # N pools: C/N ~ 20 for pool 1, ~16 for pool 2, ~12 for pool 3, ~10 for pool 4
            ns.decomp_npools_vr_col[c, j, 1] = cs.decomp_cpools_vr_col[c, j, 1] / 20.0
            ns.decomp_npools_vr_col[c, j, 2] = cs.decomp_cpools_vr_col[c, j, 2] / 16.0
            ns.decomp_npools_vr_col[c, j, 3] = cs.decomp_cpools_vr_col[c, j, 3] / 12.0
            ns.decomp_npools_vr_col[c, j, 4] = cs.decomp_cpools_vr_col[c, j, 4] / 10.0
        end

        # -- Carbon flux --
        cf = CLM.SoilBiogeochemCarbonFluxData()
        cf.rf_decomp_cascade_col   = zeros(nc, nlevdecomp, ndecomp_cascade_transitions)
        cf.pathfrac_decomp_cascade_col = ones(nc, nlevdecomp, ndecomp_cascade_transitions)
        cf.decomp_k_col            = zeros(nc, nlevdecomp, ndecomp_pools)
        cf.cn_col                  = zeros(nc, ndecomp_pools)
        cf.phr_vr_col              = zeros(nc, nlevdecomp)

        # Set respiration fractions
        for c in 1:nc, j in 1:nlevdecomp
            cf.rf_decomp_cascade_col[c, j, 1] = 0.39   # litter1 -> SOM1
            cf.rf_decomp_cascade_col[c, j, 2] = 0.55   # litter2 -> SOM1
            cf.rf_decomp_cascade_col[c, j, 3] = 0.28   # SOM1 -> SOM2
            cf.rf_decomp_cascade_col[c, j, 4] = 0.55   # SOM1 -> atm (100% resp)
            cf.rf_decomp_cascade_col[c, j, 5] = 0.55   # SOM2 -> SOM1
        end

        # Set decomposition rate constants (positive for all pools)
        for c in 1:nc, j in 1:nlevdecomp
            cf.decomp_k_col[c, j, 1] = 0.01  # litter 1
            cf.decomp_k_col[c, j, 2] = 0.005 # litter 2
            cf.decomp_k_col[c, j, 3] = 0.002 # SOM 1
            cf.decomp_k_col[c, j, 4] = 0.001 # SOM 2
        end

        # -- Biogeochem state --
        st = CLM.SoilBiogeochemStateData()
        st.nue_decomp_cascade_col = ones(ndecomp_cascade_transitions) * 0.5

        # -- Nitrogen flux --
        nf = CLM.SoilBiogeochemNitrogenFluxData()
        nf.potential_immob_vr_col = zeros(nc, nlevdecomp)
        nf.gross_nmin_vr_col      = zeros(nc, nlevdecomp)

        # -- Mask --
        mask_bgc_soilc = trues(nc)

        # -- Output arrays --
        cn_decomp_pools_out    = zeros(nc, nlevdecomp, ndecomp_pools)
        p_decomp_cpool_loss    = zeros(nc, nlevdecomp, ndecomp_cascade_transitions)
        p_decomp_cn_gain       = zeros(nc, nlevdecomp, ndecomp_pools)
        pmnf_decomp_cascade    = zeros(nc, nlevdecomp, ndecomp_cascade_transitions)
        p_decomp_npool_to_din  = zeros(nc, nlevdecomp, ndecomp_cascade_transitions)

        return (cf=cf, cs=cs, nf=nf, ns=ns, st=st, cascade_con=cascade_con,
                mask_bgc_soilc=mask_bgc_soilc,
                cn_decomp_pools=cn_decomp_pools_out,
                p_decomp_cpool_loss=p_decomp_cpool_loss,
                p_decomp_cn_gain=p_decomp_cn_gain,
                pmnf_decomp_cascade=pmnf_decomp_cascade,
                p_decomp_npool_to_din=p_decomp_npool_to_din,
                nc=nc, nlevdecomp=nlevdecomp,
                ndecomp_pools=ndecomp_pools,
                ndecomp_cascade_transitions=ndecomp_cascade_transitions)
    end

    # ---------------------------------------------------------------
    @testset "DecompPotentialParams construction and read" begin
        p = CLM.DecompPotentialParams()
        @test p.dnp == 0.01

        CLM.decomp_potential_read_params!(p; dnp=0.05)
        @test p.dnp == 0.05
    end

    # ---------------------------------------------------------------
    @testset "C:N ratio calculation -- floating vs fixed pools" begin
        d = make_test_data()

        CLM.soil_bgc_potential!(
            d.cf, d.cs, d.nf, d.ns, d.st, d.cascade_con;
            mask_bgc_soilc=d.mask_bgc_soilc,
            bounds=1:d.nc,
            nlevdecomp=d.nlevdecomp,
            ndecomp_pools=d.ndecomp_pools,
            ndecomp_cascade_transitions=d.ndecomp_cascade_transitions,
            cn_decomp_pools=d.cn_decomp_pools,
            p_decomp_cpool_loss=d.p_decomp_cpool_loss,
            p_decomp_cn_gain=d.p_decomp_cn_gain,
            pmnf_decomp_cascade=d.pmnf_decomp_cascade,
            p_decomp_npool_to_din=d.p_decomp_npool_to_din,
            use_mimics=false)

        # Pool 1 has floating C:N; should be C/N = 20.0 (from how we set it up)
        for c in 1:d.nc, j in 1:d.nlevdecomp
            @test d.cn_decomp_pools[c, j, 1] ≈ 20.0 atol=1e-10
        end

        # Pool 2 has floating C:N; should be C/N = 16.0
        for c in 1:d.nc, j in 1:d.nlevdecomp
            @test d.cn_decomp_pools[c, j, 2] ≈ 16.0 atol=1e-10
        end

        # Pool 3 has fixed C:N; should use initial_cn_ratio = 12.0
        for c in 1:d.nc, j in 1:d.nlevdecomp
            @test d.cn_decomp_pools[c, j, 3] ≈ 12.0 atol=1e-10
        end

        # Pool 4 has fixed C:N; should use initial_cn_ratio = 10.0
        for c in 1:d.nc, j in 1:d.nlevdecomp
            @test d.cn_decomp_pools[c, j, 4] ≈ 10.0 atol=1e-10
        end
    end

    # ---------------------------------------------------------------
    @testset "Potential C pool loss calculation" begin
        d = make_test_data()

        CLM.soil_bgc_potential!(
            d.cf, d.cs, d.nf, d.ns, d.st, d.cascade_con;
            mask_bgc_soilc=d.mask_bgc_soilc,
            bounds=1:d.nc,
            nlevdecomp=d.nlevdecomp,
            ndecomp_pools=d.ndecomp_pools,
            ndecomp_cascade_transitions=d.ndecomp_cascade_transitions,
            cn_decomp_pools=d.cn_decomp_pools,
            p_decomp_cpool_loss=d.p_decomp_cpool_loss,
            p_decomp_cn_gain=d.p_decomp_cn_gain,
            pmnf_decomp_cascade=d.pmnf_decomp_cascade,
            p_decomp_npool_to_din=d.p_decomp_npool_to_din,
            use_mimics=false)

        # Transition 1: donor pool 1, pathfrac=1.0
        # p_decomp_cpool_loss = cpools_vr * decomp_k * pathfrac
        for c in 1:d.nc, j in 1:d.nlevdecomp
            expected = d.cs.decomp_cpools_vr_col[c, j, 1] * d.cf.decomp_k_col[c, j, 1] * 1.0
            @test d.p_decomp_cpool_loss[c, j, 1] ≈ expected atol=1e-12
        end

        # Transition 2: donor pool 2
        for c in 1:d.nc, j in 1:d.nlevdecomp
            expected = d.cs.decomp_cpools_vr_col[c, j, 2] * d.cf.decomp_k_col[c, j, 2] * 1.0
            @test d.p_decomp_cpool_loss[c, j, 2] ≈ expected atol=1e-12
        end

        # Transition 3: donor pool 3
        for c in 1:d.nc, j in 1:d.nlevdecomp
            expected = d.cs.decomp_cpools_vr_col[c, j, 3] * d.cf.decomp_k_col[c, j, 3] * 1.0
            @test d.p_decomp_cpool_loss[c, j, 3] ≈ expected atol=1e-12
        end

        # All p_decomp_cpool_loss should be positive
        for k in 1:d.ndecomp_cascade_transitions
            for c in 1:d.nc, j in 1:d.nlevdecomp
                @test d.p_decomp_cpool_loss[c, j, k] >= 0.0
            end
        end
    end

    # ---------------------------------------------------------------
    @testset "PMNF for fixed C:N receiver (non-atmosphere)" begin
        d = make_test_data()

        CLM.soil_bgc_potential!(
            d.cf, d.cs, d.nf, d.ns, d.st, d.cascade_con;
            mask_bgc_soilc=d.mask_bgc_soilc,
            bounds=1:d.nc,
            nlevdecomp=d.nlevdecomp,
            ndecomp_pools=d.ndecomp_pools,
            ndecomp_cascade_transitions=d.ndecomp_cascade_transitions,
            cn_decomp_pools=d.cn_decomp_pools,
            p_decomp_cpool_loss=d.p_decomp_cpool_loss,
            p_decomp_cn_gain=d.p_decomp_cn_gain,
            pmnf_decomp_cascade=d.pmnf_decomp_cascade,
            p_decomp_npool_to_din=d.p_decomp_npool_to_din,
            use_mimics=false)

        # Transition 3: donor=3, receiver=4, both fixed C:N
        # pmnf = (p_decomp_cpool_loss * (1 - rf - cn_recv/cn_donor)) / cn_recv
        for c in 1:d.nc, j in 1:d.nlevdecomp
            rf = d.cf.rf_decomp_cascade_col[c, j, 3]
            cn_recv  = d.cn_decomp_pools[c, j, 4]  # 10.0
            cn_donor = d.cn_decomp_pools[c, j, 3]  # 12.0
            ratio = cn_recv / cn_donor
            loss = d.p_decomp_cpool_loss[c, j, 3]
            expected_pmnf = loss * (1.0 - rf - ratio) / cn_recv
            @test d.pmnf_decomp_cascade[c, j, 3] ≈ expected_pmnf atol=1e-12
        end
    end

    # ---------------------------------------------------------------
    @testset "PMNF for 100% respiration (receiver = atmosphere)" begin
        d = make_test_data()

        CLM.soil_bgc_potential!(
            d.cf, d.cs, d.nf, d.ns, d.st, d.cascade_con;
            mask_bgc_soilc=d.mask_bgc_soilc,
            bounds=1:d.nc,
            nlevdecomp=d.nlevdecomp,
            ndecomp_pools=d.ndecomp_pools,
            ndecomp_cascade_transitions=d.ndecomp_cascade_transitions,
            cn_decomp_pools=d.cn_decomp_pools,
            p_decomp_cpool_loss=d.p_decomp_cpool_loss,
            p_decomp_cn_gain=d.p_decomp_cn_gain,
            pmnf_decomp_cascade=d.pmnf_decomp_cascade,
            p_decomp_npool_to_din=d.p_decomp_npool_to_din,
            use_mimics=false)

        # Transition 4: donor=3, receiver=I_ATM (100% respiration)
        # pmnf = -p_decomp_cpool_loss / cn_donor
        for c in 1:d.nc, j in 1:d.nlevdecomp
            cn_donor = d.cn_decomp_pools[c, j, 3]
            loss = d.p_decomp_cpool_loss[c, j, 4]
            expected_pmnf = -loss / cn_donor
            @test d.pmnf_decomp_cascade[c, j, 4] ≈ expected_pmnf atol=1e-12
        end
    end

    # ---------------------------------------------------------------
    @testset "PMNF for floating C:N receiver (CWD -> litter path)" begin
        d = make_test_data()

        # Transition 1: donor=1, receiver=3, but pool 3 has fixed C:N
        # Transition 5: donor=4, receiver=3, pool 3 has fixed C:N
        # For the floating_cn test, let us change receiver to a floating pool
        # Modify: make transitions 1 and 2 have floating C:N receivers
        # Actually, transitions 1,2 already have receiver pool 3 which is NOT floating.
        # In the Fortran, the check is:
        #   if (.not. floating_cn_ratio(cascade_receiver_pool(k)))
        # For transition 1: receiver=3, floating[3]=false => takes non-floating path.
        # We need to test the floating path separately.

        # Create a scenario where receiver has floating C:N
        d2 = make_test_data()
        # Make pool 3 also floating
        d2.cascade_con.floating_cn_ratio_decomp_pools = BitVector([true, true, true, false])

        CLM.soil_bgc_potential!(
            d2.cf, d2.cs, d2.nf, d2.ns, d2.st, d2.cascade_con;
            mask_bgc_soilc=d2.mask_bgc_soilc,
            bounds=1:d2.nc,
            nlevdecomp=d2.nlevdecomp,
            ndecomp_pools=d2.ndecomp_pools,
            ndecomp_cascade_transitions=d2.ndecomp_cascade_transitions,
            cn_decomp_pools=d2.cn_decomp_pools,
            p_decomp_cpool_loss=d2.p_decomp_cpool_loss,
            p_decomp_cn_gain=d2.p_decomp_cn_gain,
            pmnf_decomp_cascade=d2.pmnf_decomp_cascade,
            p_decomp_npool_to_din=d2.p_decomp_npool_to_din,
            use_mimics=false)

        # For transitions with floating C:N receiver, pmnf should be 0
        # Transition 1: receiver=3 (now floating), so pmnf=0
        for c in 1:d2.nc, j in 1:d2.nlevdecomp
            @test d2.pmnf_decomp_cascade[c, j, 1] ≈ 0.0 atol=1e-12
        end
        # Transition 2: receiver=3 (now floating), so pmnf=0
        for c in 1:d2.nc, j in 1:d2.nlevdecomp
            @test d2.pmnf_decomp_cascade[c, j, 2] ≈ 0.0 atol=1e-12
        end
    end

    # ---------------------------------------------------------------
    @testset "Potential immobilization and gross mineralization" begin
        d = make_test_data()

        CLM.soil_bgc_potential!(
            d.cf, d.cs, d.nf, d.ns, d.st, d.cascade_con;
            mask_bgc_soilc=d.mask_bgc_soilc,
            bounds=1:d.nc,
            nlevdecomp=d.nlevdecomp,
            ndecomp_pools=d.ndecomp_pools,
            ndecomp_cascade_transitions=d.ndecomp_cascade_transitions,
            cn_decomp_pools=d.cn_decomp_pools,
            p_decomp_cpool_loss=d.p_decomp_cpool_loss,
            p_decomp_cn_gain=d.p_decomp_cn_gain,
            pmnf_decomp_cascade=d.pmnf_decomp_cascade,
            p_decomp_npool_to_din=d.p_decomp_npool_to_din,
            use_mimics=false)

        # Check that potential_immob_vr = sum of positive pmnf
        for c in 1:d.nc, j in 1:d.nlevdecomp
            expected_immob = 0.0
            for k in 1:d.ndecomp_cascade_transitions
                if d.pmnf_decomp_cascade[c, j, k] > 0.0
                    expected_immob += d.pmnf_decomp_cascade[c, j, k]
                end
            end
            @test d.nf.potential_immob_vr_col[c, j] ≈ expected_immob atol=1e-12
        end

        # Check that gross_nmin_vr = sum of negative(-pmnf)
        for c in 1:d.nc, j in 1:d.nlevdecomp
            expected_gmin = 0.0
            for k in 1:d.ndecomp_cascade_transitions
                if d.pmnf_decomp_cascade[c, j, k] < 0.0
                    expected_gmin -= d.pmnf_decomp_cascade[c, j, k]
                end
            end
            @test d.nf.gross_nmin_vr_col[c, j] ≈ expected_gmin atol=1e-12
        end
    end

    # ---------------------------------------------------------------
    @testset "Potential HR (phr_vr)" begin
        d = make_test_data()

        CLM.soil_bgc_potential!(
            d.cf, d.cs, d.nf, d.ns, d.st, d.cascade_con;
            mask_bgc_soilc=d.mask_bgc_soilc,
            bounds=1:d.nc,
            nlevdecomp=d.nlevdecomp,
            ndecomp_pools=d.ndecomp_pools,
            ndecomp_cascade_transitions=d.ndecomp_cascade_transitions,
            cn_decomp_pools=d.cn_decomp_pools,
            p_decomp_cpool_loss=d.p_decomp_cpool_loss,
            p_decomp_cn_gain=d.p_decomp_cn_gain,
            pmnf_decomp_cascade=d.pmnf_decomp_cascade,
            p_decomp_npool_to_din=d.p_decomp_npool_to_din,
            use_mimics=false)

        # phr_vr = sum over k of (rf * p_decomp_cpool_loss)
        for c in 1:d.nc, j in 1:d.nlevdecomp
            expected_phr = 0.0
            for k in 1:d.ndecomp_cascade_transitions
                expected_phr += d.cf.rf_decomp_cascade_col[c, j, k] * d.p_decomp_cpool_loss[c, j, k]
            end
            @test d.cf.phr_vr_col[c, j] ≈ expected_phr atol=1e-12
        end

        # phr_vr should be non-negative
        for c in 1:d.nc, j in 1:d.nlevdecomp
            @test d.cf.phr_vr_col[c, j] >= 0.0
        end
    end

    # ---------------------------------------------------------------
    @testset "Mask filtering -- masked columns are untouched" begin
        d = make_test_data(nc=4)

        # Mask out columns 2 and 4
        d.mask_bgc_soilc[2] = false
        d.mask_bgc_soilc[4] = false

        CLM.soil_bgc_potential!(
            d.cf, d.cs, d.nf, d.ns, d.st, d.cascade_con;
            mask_bgc_soilc=d.mask_bgc_soilc,
            bounds=1:4,
            nlevdecomp=d.nlevdecomp,
            ndecomp_pools=d.ndecomp_pools,
            ndecomp_cascade_transitions=d.ndecomp_cascade_transitions,
            cn_decomp_pools=d.cn_decomp_pools,
            p_decomp_cpool_loss=d.p_decomp_cpool_loss,
            p_decomp_cn_gain=d.p_decomp_cn_gain,
            pmnf_decomp_cascade=d.pmnf_decomp_cascade,
            p_decomp_npool_to_din=d.p_decomp_npool_to_din,
            use_mimics=false)

        # Active columns (1, 3) should have non-zero potential HR
        for j in 1:d.nlevdecomp
            @test d.cf.phr_vr_col[1, j] > 0.0
            @test d.cf.phr_vr_col[3, j] > 0.0
        end

        # Masked columns (2, 4) should remain zero (initialized to zero)
        for j in 1:d.nlevdecomp
            @test d.cf.phr_vr_col[2, j] == 0.0
            @test d.cf.phr_vr_col[4, j] == 0.0
        end

        # potential_immob_vr should be zero for masked columns
        for j in 1:d.nlevdecomp
            @test d.nf.potential_immob_vr_col[2, j] == 0.0
            @test d.nf.potential_immob_vr_col[4, j] == 0.0
        end
    end

    # ---------------------------------------------------------------
    @testset "Zero pools -- no decomposition" begin
        d = make_test_data()

        # Zero out all C pools
        d.cs.decomp_cpools_vr_col .= 0.0

        CLM.soil_bgc_potential!(
            d.cf, d.cs, d.nf, d.ns, d.st, d.cascade_con;
            mask_bgc_soilc=d.mask_bgc_soilc,
            bounds=1:d.nc,
            nlevdecomp=d.nlevdecomp,
            ndecomp_pools=d.ndecomp_pools,
            ndecomp_cascade_transitions=d.ndecomp_cascade_transitions,
            cn_decomp_pools=d.cn_decomp_pools,
            p_decomp_cpool_loss=d.p_decomp_cpool_loss,
            p_decomp_cn_gain=d.p_decomp_cn_gain,
            pmnf_decomp_cascade=d.pmnf_decomp_cascade,
            p_decomp_npool_to_din=d.p_decomp_npool_to_din,
            use_mimics=false)

        # Everything should be zero
        @test all(d.p_decomp_cpool_loss .== 0.0)
        @test all(d.pmnf_decomp_cascade .== 0.0)
        @test all(d.cf.phr_vr_col .== 0.0)
        @test all(d.nf.potential_immob_vr_col .== 0.0)
        @test all(d.nf.gross_nmin_vr_col .== 0.0)
    end

    # ---------------------------------------------------------------
    @testset "Zero decomp rate constants -- no decomposition" begin
        d = make_test_data()

        # Zero out all decomp rate constants
        d.cf.decomp_k_col .= 0.0

        CLM.soil_bgc_potential!(
            d.cf, d.cs, d.nf, d.ns, d.st, d.cascade_con;
            mask_bgc_soilc=d.mask_bgc_soilc,
            bounds=1:d.nc,
            nlevdecomp=d.nlevdecomp,
            ndecomp_pools=d.ndecomp_pools,
            ndecomp_cascade_transitions=d.ndecomp_cascade_transitions,
            cn_decomp_pools=d.cn_decomp_pools,
            p_decomp_cpool_loss=d.p_decomp_cpool_loss,
            p_decomp_cn_gain=d.p_decomp_cn_gain,
            pmnf_decomp_cascade=d.pmnf_decomp_cascade,
            p_decomp_npool_to_din=d.p_decomp_npool_to_din,
            use_mimics=false)

        @test all(d.p_decomp_cpool_loss .== 0.0)
        @test all(d.cf.phr_vr_col .== 0.0)
    end

    # ---------------------------------------------------------------
    @testset "Pathfrac scaling" begin
        d = make_test_data()

        # Set pathfrac to 0.5 for transition 1
        d.cf.pathfrac_decomp_cascade_col[:, :, 1] .= 0.5

        CLM.soil_bgc_potential!(
            d.cf, d.cs, d.nf, d.ns, d.st, d.cascade_con;
            mask_bgc_soilc=d.mask_bgc_soilc,
            bounds=1:d.nc,
            nlevdecomp=d.nlevdecomp,
            ndecomp_pools=d.ndecomp_pools,
            ndecomp_cascade_transitions=d.ndecomp_cascade_transitions,
            cn_decomp_pools=d.cn_decomp_pools,
            p_decomp_cpool_loss=d.p_decomp_cpool_loss,
            p_decomp_cn_gain=d.p_decomp_cn_gain,
            pmnf_decomp_cascade=d.pmnf_decomp_cascade,
            p_decomp_npool_to_din=d.p_decomp_npool_to_din,
            use_mimics=false)

        # Transition 1 C pool loss should be half of what it would be with pathfrac=1
        for c in 1:d.nc, j in 1:d.nlevdecomp
            expected = d.cs.decomp_cpools_vr_col[c, j, 1] * d.cf.decomp_k_col[c, j, 1] * 0.5
            @test d.p_decomp_cpool_loss[c, j, 1] ≈ expected atol=1e-12
        end
    end

    # ---------------------------------------------------------------
    @testset "MIMICS pathway -- basic N flux to DIN" begin
        # Build a MIMICS-like configuration with microbial pools
        nc = 2
        nlevdecomp = 2
        ndecomp_pools = 6       # lit_met, lit_str, som_phys, cop_mic, oli_mic, cwd
        ndecomp_cascade_transitions = 3
        i_cop_mic_val = 4
        i_oli_mic_val = 5

        cascade_con = CLM.DecompCascadeConData()
        cascade_con.cascade_donor_pool    = [1, 2, 3]         # lit_met, lit_str, som_phys
        cascade_con.cascade_receiver_pool = [i_cop_mic_val, i_oli_mic_val, i_cop_mic_val]
        # All receivers are microbial => floating C:N
        cascade_con.floating_cn_ratio_decomp_pools = BitVector([true, true, false, true, true, true])
        cascade_con.initial_cn_ratio = [20.0, 25.0, 12.0, 10.0, 10.0, 500.0]

        cs = CLM.SoilBiogeochemCarbonStateData()
        cs.decomp_cpools_vr_col = zeros(nc, nlevdecomp, ndecomp_pools)
        cs.decomp_cpools_vr_col[:, :, 1] .= 100.0
        cs.decomp_cpools_vr_col[:, :, 2] .= 80.0
        cs.decomp_cpools_vr_col[:, :, 3] .= 50.0
        cs.decomp_cpools_vr_col[:, :, 4] .= 20.0
        cs.decomp_cpools_vr_col[:, :, 5] .= 15.0

        ns = CLM.SoilBiogeochemNitrogenStateData()
        ns.decomp_npools_vr_col = zeros(nc, nlevdecomp, ndecomp_pools)
        ns.decomp_npools_vr_col[:, :, 1] .= 5.0    # C:N = 20
        ns.decomp_npools_vr_col[:, :, 2] .= 3.2    # C:N = 25
        ns.decomp_npools_vr_col[:, :, 3] .= 4.167  # C:N ~ 12
        ns.decomp_npools_vr_col[:, :, 4] .= 2.0    # C:N = 10
        ns.decomp_npools_vr_col[:, :, 5] .= 1.5    # C:N = 10

        cf = CLM.SoilBiogeochemCarbonFluxData()
        cf.rf_decomp_cascade_col = zeros(nc, nlevdecomp, ndecomp_cascade_transitions)
        cf.rf_decomp_cascade_col[:, :, 1] .= 0.3
        cf.rf_decomp_cascade_col[:, :, 2] .= 0.4
        cf.rf_decomp_cascade_col[:, :, 3] .= 0.2
        cf.pathfrac_decomp_cascade_col = ones(nc, nlevdecomp, ndecomp_cascade_transitions)
        cf.decomp_k_col = zeros(nc, nlevdecomp, ndecomp_pools)
        cf.decomp_k_col[:, :, 1] .= 0.01
        cf.decomp_k_col[:, :, 2] .= 0.005
        cf.decomp_k_col[:, :, 3] .= 0.002
        cf.cn_col = zeros(nc, ndecomp_pools)
        cf.cn_col[:, i_cop_mic_val] .= 8.0   # target C:N for cop microbes
        cf.cn_col[:, i_oli_mic_val] .= 8.0   # target C:N for oli microbes
        cf.phr_vr_col = zeros(nc, nlevdecomp)

        st = CLM.SoilBiogeochemStateData()
        st.nue_decomp_cascade_col = [0.5, 0.6, 0.5]   # NUE for each transition

        nf = CLM.SoilBiogeochemNitrogenFluxData()
        nf.potential_immob_vr_col = zeros(nc, nlevdecomp)
        nf.gross_nmin_vr_col      = zeros(nc, nlevdecomp)

        mask = trues(nc)
        cn_out = zeros(nc, nlevdecomp, ndecomp_pools)
        cpool_loss = zeros(nc, nlevdecomp, ndecomp_cascade_transitions)
        cn_gain = zeros(nc, nlevdecomp, ndecomp_pools)
        pmnf = zeros(nc, nlevdecomp, ndecomp_cascade_transitions)
        npool_din = zeros(nc, nlevdecomp, ndecomp_cascade_transitions)

        CLM.soil_bgc_potential!(
            cf, cs, nf, ns, st, cascade_con;
            mask_bgc_soilc=mask,
            bounds=1:nc,
            nlevdecomp=nlevdecomp,
            ndecomp_pools=ndecomp_pools,
            ndecomp_cascade_transitions=ndecomp_cascade_transitions,
            cn_decomp_pools=cn_out,
            p_decomp_cpool_loss=cpool_loss,
            p_decomp_cn_gain=cn_gain,
            pmnf_decomp_cascade=pmnf,
            p_decomp_npool_to_din=npool_din,
            use_mimics=true,
            i_cop_mic=i_cop_mic_val,
            i_oli_mic=i_oli_mic_val)

        # Verify that N flux to DIN is computed (non-zero for MIMICS transitions)
        for k in 1:ndecomp_cascade_transitions
            for c in 1:nc, j in 1:nlevdecomp
                # Transitions with floating receiver go through MIMICS code path
                if cascade_con.floating_cn_ratio_decomp_pools[cascade_con.cascade_receiver_pool[k]]
                    # Check basic MIMICS N accounting:
                    # p_decomp_npool_to_din = p_decomp_npool_loss - p_decomp_npool_gain
                    # where p_decomp_npool_gain = p_decomp_npool_loss * nue
                    # so p_decomp_npool_to_din = p_decomp_npool_loss * (1 - nue)
                    donor = cascade_con.cascade_donor_pool[k]
                    nue = st.nue_decomp_cascade_col[k]
                    nc_donor = ns.decomp_npools_vr_col[c, j, donor] / cs.decomp_cpools_vr_col[c, j, donor]
                    p_npool_loss = cpool_loss[c, j, k] * nc_donor
                    expected_din = p_npool_loss * (1.0 - nue)
                    @test npool_din[c, j, k] ≈ expected_din atol=1e-12
                end
            end
        end

        # Gross mineralization should include N flux to DIN from MIMICS
        for c in 1:nc, j in 1:nlevdecomp
            total_din = sum(npool_din[c, j, :])
            # gross_nmin should be at least as large as total DIN flux
            @test nf.gross_nmin_vr_col[c, j] >= total_din - 1e-12
        end
    end

    # ---------------------------------------------------------------
    @testset "Single column / single level consistency" begin
        # Minimal 1x1 test case to verify formula manually
        nc = 1; nlev = 1; npools = 2; ntrans = 1

        cascade_con = CLM.DecompCascadeConData()
        cascade_con.cascade_donor_pool    = [1]
        cascade_con.cascade_receiver_pool = [2]
        cascade_con.floating_cn_ratio_decomp_pools = BitVector([false, false])
        cascade_con.initial_cn_ratio = [20.0, 10.0]

        cs = CLM.SoilBiogeochemCarbonStateData()
        cs.decomp_cpools_vr_col = zeros(1, 1, 2)
        cs.decomp_cpools_vr_col[1, 1, 1] = 200.0  # donor pool
        cs.decomp_cpools_vr_col[1, 1, 2] = 100.0  # receiver pool

        ns = CLM.SoilBiogeochemNitrogenStateData()
        ns.decomp_npools_vr_col = zeros(1, 1, 2)
        ns.decomp_npools_vr_col[1, 1, 1] = 10.0  # C:N = 20
        ns.decomp_npools_vr_col[1, 1, 2] = 10.0  # C:N = 10

        cf = CLM.SoilBiogeochemCarbonFluxData()
        cf.rf_decomp_cascade_col = zeros(1, 1, 1)
        cf.rf_decomp_cascade_col[1, 1, 1] = 0.4
        cf.pathfrac_decomp_cascade_col = ones(1, 1, 1)
        cf.decomp_k_col = zeros(1, 1, 2)
        cf.decomp_k_col[1, 1, 1] = 0.02
        cf.cn_col = zeros(1, 2)
        cf.phr_vr_col = zeros(1, 1)

        st = CLM.SoilBiogeochemStateData()
        st.nue_decomp_cascade_col = [0.5]

        nf = CLM.SoilBiogeochemNitrogenFluxData()
        nf.potential_immob_vr_col = zeros(1, 1)
        nf.gross_nmin_vr_col      = zeros(1, 1)

        mask = trues(1)
        cn_out = zeros(1, 1, 2)
        cpool_loss = zeros(1, 1, 1)
        cn_gain = zeros(1, 1, 2)
        pmnf = zeros(1, 1, 1)
        npool_din = zeros(1, 1, 1)

        CLM.soil_bgc_potential!(
            cf, cs, nf, ns, st, cascade_con;
            mask_bgc_soilc=mask,
            bounds=1:1,
            nlevdecomp=1,
            ndecomp_pools=2,
            ndecomp_cascade_transitions=1,
            cn_decomp_pools=cn_out,
            p_decomp_cpool_loss=cpool_loss,
            p_decomp_cn_gain=cn_gain,
            pmnf_decomp_cascade=pmnf,
            p_decomp_npool_to_din=npool_din,
            use_mimics=false)

        # Manual calculation:
        # C:N pools: both fixed => cn[1]=20, cn[2]=10
        @test cn_out[1, 1, 1] ≈ 20.0
        @test cn_out[1, 1, 2] ≈ 10.0

        # p_decomp_cpool_loss = 200 * 0.02 * 1.0 = 4.0
        @test cpool_loss[1, 1, 1] ≈ 4.0

        # ratio = cn_recv/cn_donor = 10/20 = 0.5
        # pmnf = 4.0 * (1.0 - 0.4 - 0.5) / 10 = 4.0 * 0.1 / 10 = 0.04
        @test pmnf[1, 1, 1] ≈ 0.04 atol=1e-12

        # pmnf > 0 => immobilization
        @test nf.potential_immob_vr_col[1, 1] ≈ 0.04 atol=1e-12
        @test nf.gross_nmin_vr_col[1, 1] ≈ 0.0 atol=1e-12

        # phr_vr = rf * cpool_loss = 0.4 * 4.0 = 1.6
        @test cf.phr_vr_col[1, 1] ≈ 1.6 atol=1e-12
    end

    # ---------------------------------------------------------------
    @testset "Mineralization case (pmnf negative)" begin
        # Create a scenario where mineralization dominates
        nc = 1; nlev = 1; npools = 2; ntrans = 1

        cascade_con = CLM.DecompCascadeConData()
        cascade_con.cascade_donor_pool    = [1]
        cascade_con.cascade_receiver_pool = [CLM.I_ATM]  # 100% respiration
        cascade_con.floating_cn_ratio_decomp_pools = BitVector([false, false])
        cascade_con.initial_cn_ratio = [15.0, 10.0]

        cs = CLM.SoilBiogeochemCarbonStateData()
        cs.decomp_cpools_vr_col = zeros(1, 1, 2)
        cs.decomp_cpools_vr_col[1, 1, 1] = 150.0

        ns = CLM.SoilBiogeochemNitrogenStateData()
        ns.decomp_npools_vr_col = zeros(1, 1, 2)
        ns.decomp_npools_vr_col[1, 1, 1] = 10.0  # C:N = 15

        cf = CLM.SoilBiogeochemCarbonFluxData()
        cf.rf_decomp_cascade_col = zeros(1, 1, 1)
        cf.rf_decomp_cascade_col[1, 1, 1] = 1.0  # 100% resp
        cf.pathfrac_decomp_cascade_col = ones(1, 1, 1)
        cf.decomp_k_col = zeros(1, 1, 2)
        cf.decomp_k_col[1, 1, 1] = 0.01
        cf.cn_col = zeros(1, 2)
        cf.phr_vr_col = zeros(1, 1)

        st = CLM.SoilBiogeochemStateData()
        st.nue_decomp_cascade_col = [0.5]

        nf = CLM.SoilBiogeochemNitrogenFluxData()
        nf.potential_immob_vr_col = zeros(1, 1)
        nf.gross_nmin_vr_col      = zeros(1, 1)

        mask = trues(1)
        cn_out = zeros(1, 1, 2)
        cpool_loss = zeros(1, 1, 1)
        cn_gain = zeros(1, 1, 2)
        pmnf_out = zeros(1, 1, 1)
        npool_din = zeros(1, 1, 1)

        CLM.soil_bgc_potential!(
            cf, cs, nf, ns, st, cascade_con;
            mask_bgc_soilc=mask,
            bounds=1:1,
            nlevdecomp=1,
            ndecomp_pools=2,
            ndecomp_cascade_transitions=1,
            cn_decomp_pools=cn_out,
            p_decomp_cpool_loss=cpool_loss,
            p_decomp_cn_gain=cn_gain,
            pmnf_decomp_cascade=pmnf_out,
            p_decomp_npool_to_din=npool_din,
            use_mimics=false)

        # cpool_loss = 150 * 0.01 * 1 = 1.5
        @test cpool_loss[1, 1, 1] ≈ 1.5

        # Receiver is atmosphere => pmnf = -cpool_loss / cn_donor = -1.5 / 15 = -0.1
        @test pmnf_out[1, 1, 1] ≈ -0.1 atol=1e-12

        # pmnf < 0 => mineralization
        @test nf.potential_immob_vr_col[1, 1] ≈ 0.0 atol=1e-12
        @test nf.gross_nmin_vr_col[1, 1] ≈ 0.1 atol=1e-12

        # phr_vr = 1.0 * 1.5 = 1.5
        @test cf.phr_vr_col[1, 1] ≈ 1.5 atol=1e-12
    end

    # ---------------------------------------------------------------
    @testset "Multiple levels produce independent results" begin
        d = make_test_data(nc=1, nlevdecomp=3)

        # Make each level have different C pools
        d.cs.decomp_cpools_vr_col[1, 1, :] .= [100.0, 80.0, 50.0, 30.0]
        d.cs.decomp_cpools_vr_col[1, 2, :] .= [200.0, 160.0, 100.0, 60.0]
        d.cs.decomp_cpools_vr_col[1, 3, :] .= [50.0, 40.0, 25.0, 15.0]

        for j in 1:3, l in 1:4
            d.ns.decomp_npools_vr_col[1, j, l] = d.cs.decomp_cpools_vr_col[1, j, l] / d.cascade_con.initial_cn_ratio[l]
        end

        CLM.soil_bgc_potential!(
            d.cf, d.cs, d.nf, d.ns, d.st, d.cascade_con;
            mask_bgc_soilc=d.mask_bgc_soilc,
            bounds=1:1,
            nlevdecomp=3,
            ndecomp_pools=d.ndecomp_pools,
            ndecomp_cascade_transitions=d.ndecomp_cascade_transitions,
            cn_decomp_pools=d.cn_decomp_pools,
            p_decomp_cpool_loss=d.p_decomp_cpool_loss,
            p_decomp_cn_gain=d.p_decomp_cn_gain,
            pmnf_decomp_cascade=d.pmnf_decomp_cascade,
            p_decomp_npool_to_din=d.p_decomp_npool_to_din,
            use_mimics=false)

        # Level 2 has 2x the C of level 1, so cpool_loss should be 2x
        for k in 1:d.ndecomp_cascade_transitions
            @test d.p_decomp_cpool_loss[1, 2, k] ≈ 2.0 * d.p_decomp_cpool_loss[1, 1, k] atol=1e-10
        end

        # Level 3 has 0.5x the C of level 1
        for k in 1:d.ndecomp_cascade_transitions
            @test d.p_decomp_cpool_loss[1, 3, k] ≈ 0.5 * d.p_decomp_cpool_loss[1, 1, k] atol=1e-10
        end
    end

    # ---------------------------------------------------------------
    @testset "N pool zero -- ratio check for non-atmosphere" begin
        # When N pool of donor is zero, ratio should be set to 0
        nc = 1; nlev = 1; npools = 2; ntrans = 1

        cascade_con = CLM.DecompCascadeConData()
        cascade_con.cascade_donor_pool    = [1]
        cascade_con.cascade_receiver_pool = [2]
        cascade_con.floating_cn_ratio_decomp_pools = BitVector([false, false])
        cascade_con.initial_cn_ratio = [20.0, 10.0]

        cs = CLM.SoilBiogeochemCarbonStateData()
        cs.decomp_cpools_vr_col = zeros(1, 1, 2)
        cs.decomp_cpools_vr_col[1, 1, 1] = 200.0
        cs.decomp_cpools_vr_col[1, 1, 2] = 100.0

        ns = CLM.SoilBiogeochemNitrogenStateData()
        ns.decomp_npools_vr_col = zeros(1, 1, 2)
        # Donor N = 0 => ratio should be 0
        ns.decomp_npools_vr_col[1, 1, 1] = 0.0
        ns.decomp_npools_vr_col[1, 1, 2] = 10.0

        cf = CLM.SoilBiogeochemCarbonFluxData()
        cf.rf_decomp_cascade_col = zeros(1, 1, 1)
        cf.rf_decomp_cascade_col[1, 1, 1] = 0.4
        cf.pathfrac_decomp_cascade_col = ones(1, 1, 1)
        cf.decomp_k_col = zeros(1, 1, 2)
        cf.decomp_k_col[1, 1, 1] = 0.02
        cf.cn_col = zeros(1, 2)
        cf.phr_vr_col = zeros(1, 1)

        st = CLM.SoilBiogeochemStateData()
        st.nue_decomp_cascade_col = [0.5]

        nf = CLM.SoilBiogeochemNitrogenFluxData()
        nf.potential_immob_vr_col = zeros(1, 1)
        nf.gross_nmin_vr_col      = zeros(1, 1)

        mask = trues(1)
        cn_out = zeros(1, 1, 2)
        cpool_loss = zeros(1, 1, 1)
        cn_gain = zeros(1, 1, 2)
        pmnf_out = zeros(1, 1, 1)
        npool_din = zeros(1, 1, 1)

        CLM.soil_bgc_potential!(
            cf, cs, nf, ns, st, cascade_con;
            mask_bgc_soilc=mask,
            bounds=1:1,
            nlevdecomp=1,
            ndecomp_pools=2,
            ndecomp_cascade_transitions=1,
            cn_decomp_pools=cn_out,
            p_decomp_cpool_loss=cpool_loss,
            p_decomp_cn_gain=cn_gain,
            pmnf_decomp_cascade=pmnf_out,
            p_decomp_npool_to_din=npool_din,
            use_mimics=false)

        # With N=0 for donor, ratio=0, so:
        # pmnf = cpool_loss * (1 - rf - 0) / cn_recv = 4.0 * 0.6 / 10 = 0.24
        @test cpool_loss[1, 1, 1] ≈ 4.0
        @test pmnf_out[1, 1, 1] ≈ 0.24 atol=1e-12
    end

    # ---------------------------------------------------------------
    @testset "Conservation: immobilization + mineralization consistent with pmnf" begin
        d = make_test_data(nc=5, nlevdecomp=3)

        CLM.soil_bgc_potential!(
            d.cf, d.cs, d.nf, d.ns, d.st, d.cascade_con;
            mask_bgc_soilc=d.mask_bgc_soilc,
            bounds=1:5,
            nlevdecomp=3,
            ndecomp_pools=d.ndecomp_pools,
            ndecomp_cascade_transitions=d.ndecomp_cascade_transitions,
            cn_decomp_pools=d.cn_decomp_pools,
            p_decomp_cpool_loss=d.p_decomp_cpool_loss,
            p_decomp_cn_gain=d.p_decomp_cn_gain,
            pmnf_decomp_cascade=d.pmnf_decomp_cascade,
            p_decomp_npool_to_din=d.p_decomp_npool_to_din,
            use_mimics=false)

        # potential_immob + gross_nmin should equal sum(|pmnf|)
        for c in 1:5, j in 1:3
            total_abs_pmnf = sum(abs.(d.pmnf_decomp_cascade[c, j, :]))
            @test (d.nf.potential_immob_vr_col[c, j] + d.nf.gross_nmin_vr_col[c, j]) ≈
                  total_abs_pmnf atol=1e-10
        end
    end

end
