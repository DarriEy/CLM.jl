@testset "Decomp BGC" begin

    # ----------------------------------------------------------------
    # Helper: create minimal data structures for nc columns
    # ----------------------------------------------------------------
    function make_decomp_bgc_data(; nc=4, nlevdecomp=1, ndecomp_pools=7,
                                    ndecomp_cascade_transitions=10)
        # --- BGC Parameters ---
        params = CLM.DecompBGCParams(
            cn_s1_bgc     = 12.0,
            cn_s2_bgc     = 12.0,
            cn_s3_bgc     = 10.0,
            rf_l1s1_bgc   = 0.39,
            rf_l2s1_bgc   = 0.55,
            rf_l3s2_bgc   = 0.29,
            rf_s2s1_bgc   = 0.55,
            rf_s2s3_bgc   = 0.55,
            rf_s3s1_bgc   = 0.55,
            rf_cwdl3_bgc  = 0.0,
            tau_l1_bgc    = 1.0 / 18.5,
            tau_l2_l3_bgc = 1.0 / 4.9,
            tau_s1_bgc    = 1.0 / 7.3,
            tau_s2_bgc    = 1.0 / 0.2,
            tau_s3_bgc    = 1.0 / 0.0045,
            cwd_fcel_bgc  = 0.45,
            bgc_initial_Cstocks       = [200.0, 200.0, 200.0, 200.0, 200.0, 200.0, 200.0],
            bgc_initial_Cstocks_depth = 0.3,
        )

        # --- CN shared params ---
        cn_params = CLM.CNSharedParamsData(
            Q10                   = 1.5,
            minpsi                = -10.0,
            maxpsi                = -0.1,
            rf_cwdl2              = 0.0,
            tau_cwd               = 10.0,
            cwd_flig              = 0.24,
            froz_q10              = 1.5,
            decomp_depth_efolding = 0.5,
            mino2lim              = 0.0,
        )

        # --- State and cascade con ---
        bgc_state   = CLM.DecompBGCState()
        cascade_con = CLM.DecompCascadeConData()

        # --- Soil sand values (col × nlevdecomp) ---
        cellsand = fill(50.0, nc, max(nlevdecomp, 5))

        # --- Soil temperature (col × nlev) ---
        t_soisno = fill(CLM.TFRZ + 15.0, nc, max(nlevdecomp, 5))

        # --- Soil water potential (col × nlev) ---
        soilpsi = fill(-1.0, nc, max(nlevdecomp, 5))

        # --- zsoi (soil depths) ---
        zsoi_vals = [0.01, 0.04, 0.09, 0.16, 0.26, 0.40, 0.58, 0.80, 1.06, 1.36]
        if nlevdecomp > length(zsoi_vals)
            append!(zsoi_vals, range(1.7, stop=3.0, length=nlevdecomp - length(zsoi_vals)))
        end

        # --- col_dz (layer thicknesses, for nlevdecomp==1) ---
        col_dz = fill(0.1, nc, max(nlevdecomp, 5))

        # --- Carbon flux data ---
        cf = CLM.SoilBiogeochemCarbonFluxData()
        CLM.soil_bgc_carbon_flux_init!(cf, nc, max(nlevdecomp, 5), ndecomp_pools,
                                        ndecomp_cascade_transitions)

        # --- Mask ---
        mask_bgc_soilc = trues(nc)

        return (params=params, cn_params=cn_params, bgc_state=bgc_state,
                cascade_con=cascade_con, cf=cf, cellsand=cellsand,
                t_soisno=t_soisno, soilpsi=soilpsi, zsoi_vals=zsoi_vals,
                col_dz=col_dz, mask_bgc_soilc=mask_bgc_soilc,
                nlevdecomp=nlevdecomp, ndecomp_pools=ndecomp_pools,
                ndecomp_cascade_transitions=ndecomp_cascade_transitions, nc=nc)
    end

    # ================================================================
    # Test 1: DecompBGCParams construction and read
    # ================================================================
    @testset "DecompBGCParams construction" begin
        params = CLM.DecompBGCParams()
        @test params.cn_s1_bgc == 12.0
        @test params.tau_l1_bgc ≈ 1.0 / 18.5
        @test params.rf_l1s1_bgc == 0.39

        # Test read_params
        params2 = CLM.DecompBGCParams()
        CLM.decomp_bgc_read_params!(params2;
            tau_l1=0.1, tau_l2_l3=0.2, tau_s1=0.3, tau_s2=0.4, tau_s3=0.5,
            cn_s1=8.0, cn_s2=11.0, cn_s3=7.0,
            rf_l1s1=0.4, rf_l2s1=0.5, rf_l3s2=0.3,
            rf_s2s1=0.6, rf_s2s3=0.6, rf_s3s1=0.6,
            rf_cwdl3=0.1, cwd_fcel=0.5,
            bgc_initial_Cstocks=fill(100.0, 7),
            bgc_initial_Cstocks_depth=0.2)
        @test params2.tau_l1_bgc == 0.1
        @test params2.cn_s1_bgc == 8.0
        @test params2.bgc_initial_Cstocks_depth == 0.2
    end

    # ================================================================
    # Test 2: init_decomp_cascade_bgc!
    # ================================================================
    @testset "init_decomp_cascade_bgc!" begin
        d = make_decomp_bgc_data()

        CLM.init_decomp_cascade_bgc!(
            d.bgc_state, d.cascade_con, d.params, d.cn_params;
            cellsand=d.cellsand,
            bounds=1:d.nc,
            nlevdecomp=d.nlevdecomp,
            ndecomp_pools_max=d.ndecomp_pools,
            ndecomp_cascade_transitions_max=d.ndecomp_cascade_transitions,
            use_fates=false)

        # Check pool indices
        @test d.bgc_state.i_cel_lit == 2
        @test d.bgc_state.i_lig_lit == 3
        @test d.bgc_state.i_act_som == 4
        @test d.bgc_state.i_slo_som == 5
        @test d.bgc_state.i_pas_som == 6

        # Check transition indices
        @test d.bgc_state.i_l1s1 == 1
        @test d.bgc_state.i_l2s1 == 2
        @test d.bgc_state.i_s3s1 == 8
        @test d.bgc_state.i_cwdl3 == 10

        # Check cascade donor/receiver pools
        @test d.cascade_con.cascade_donor_pool[1] == 1    # L1 -> S1
        @test d.cascade_con.cascade_receiver_pool[1] == 4 # -> active SOM
        @test d.cascade_con.cascade_donor_pool[8] == 6    # S3 -> S1
        @test d.cascade_con.cascade_receiver_pool[8] == 4 # -> active SOM

        # CWD transitions
        @test d.cascade_con.cascade_step_name[9] == "CWDL2"
        @test d.cascade_con.cascade_step_name[10] == "CWDL3"

        # Pool flags
        @test d.cascade_con.is_litter[1] == true
        @test d.cascade_con.is_litter[2] == true
        @test d.cascade_con.is_litter[3] == true
        @test d.cascade_con.is_soil[4] == true
        @test d.cascade_con.is_soil[5] == true
        @test d.cascade_con.is_soil[6] == true
        @test d.cascade_con.is_cwd[7] == true

        # Path fractions
        @test d.bgc_state.f_s2s1 ≈ 0.42 / 0.45
        @test d.bgc_state.f_s2s3 ≈ 0.03 / 0.45

        # Sand-dependent fractions should be computed
        @test d.bgc_state.f_s1s2[1, 1] > 0.0
        @test d.bgc_state.f_s1s3[1, 1] > 0.0
        @test d.bgc_state.f_s1s2[1, 1] + d.bgc_state.f_s1s3[1, 1] ≈ 1.0

        # Respiration fractions
        @test d.bgc_state.rf_l1s1 == 0.39
        @test d.bgc_state.rf_s2s1 == 0.55

        # Spinup factors
        @test d.cascade_con.spinup_factor[1] == 1.0  # metabolic litter
        @test d.cascade_con.spinup_factor[4] == 1.0  # active SOM
        @test d.cascade_con.spinup_factor[5] >= 1.0  # slow SOM

        # Initial stocks
        @test d.cascade_con.initial_stock[1] == 200.0
        @test d.cascade_con.initial_stock_soildepth == 0.3
    end

    # ================================================================
    # Test 3: decomp_rate_constants_bgc! (nlevdecomp == 1)
    # ================================================================
    @testset "decomp_rate_constants_bgc! (1 level)" begin
        d = make_decomp_bgc_data(; nlevdecomp=1)

        CLM.init_decomp_cascade_bgc!(
            d.bgc_state, d.cascade_con, d.params, d.cn_params;
            cellsand=d.cellsand,
            bounds=1:d.nc,
            nlevdecomp=d.nlevdecomp,
            ndecomp_pools_max=d.ndecomp_pools,
            ndecomp_cascade_transitions_max=d.ndecomp_cascade_transitions,
            use_fates=false)

        CLM.decomp_rate_constants_bgc!(
            d.cf, d.bgc_state, d.params, d.cn_params, d.cascade_con;
            mask_bgc_soilc=d.mask_bgc_soilc,
            bounds=1:d.nc,
            nlevdecomp=d.nlevdecomp,
            t_soisno=d.t_soisno,
            soilpsi=d.soilpsi,
            days_per_year=365.0,
            dt=1800.0,
            zsoi_vals=d.zsoi_vals,
            col_dz=d.col_dz)

        # t_scalar should be positive (warm soil)
        for c in 1:d.nc
            @test d.cf.t_scalar_col[c, 1] > 0.0
        end

        # w_scalar should be in [0, 1] for reasonable soilpsi
        for c in 1:d.nc
            @test d.cf.w_scalar_col[c, 1] >= 0.0
            @test d.cf.w_scalar_col[c, 1] <= 1.0
        end

        # o_scalar should be 1.0 (no anoxia by default)
        for c in 1:d.nc
            @test d.cf.o_scalar_col[c, 1] == 1.0
        end

        # decomp_k should be positive for all pools
        for c in 1:d.nc
            for pool in 1:d.ndecomp_pools
                @test d.cf.decomp_k_col[c, 1, pool] > 0.0
            end
        end

        # pathfrac for L1S1 should be 1.0
        @test d.cf.pathfrac_decomp_cascade_col[1, 1, 1] == 1.0
        @test d.cf.pathfrac_decomp_cascade_col[1, 1, d.bgc_state.i_s3s1] == 1.0

        # rf for L1S1 should match params
        @test d.cf.rf_decomp_cascade_col[1, 1, 1] == 0.39

        # Rate ordering: metabolic litter decomposes faster than slow SOM
        k_met = d.cf.decomp_k_col[1, 1, 1]      # metabolic litter
        k_slo = d.cf.decomp_k_col[1, 1, 5]       # slow SOM
        k_pas = d.cf.decomp_k_col[1, 1, 6]       # passive SOM
        @test k_met > k_slo > k_pas
    end

    # ================================================================
    # Test 4: decomp_rate_constants_bgc! (multi-level)
    # ================================================================
    @testset "decomp_rate_constants_bgc! (multi-level)" begin
        nlevdecomp = 5
        d = make_decomp_bgc_data(; nlevdecomp=nlevdecomp)

        CLM.init_decomp_cascade_bgc!(
            d.bgc_state, d.cascade_con, d.params, d.cn_params;
            cellsand=d.cellsand,
            bounds=1:d.nc,
            nlevdecomp=nlevdecomp,
            ndecomp_pools_max=d.ndecomp_pools,
            ndecomp_cascade_transitions_max=d.ndecomp_cascade_transitions,
            use_fates=false)

        CLM.decomp_rate_constants_bgc!(
            d.cf, d.bgc_state, d.params, d.cn_params, d.cascade_con;
            mask_bgc_soilc=d.mask_bgc_soilc,
            bounds=1:d.nc,
            nlevdecomp=nlevdecomp,
            t_soisno=d.t_soisno,
            soilpsi=d.soilpsi,
            days_per_year=365.0,
            dt=1800.0,
            zsoi_vals=d.zsoi_vals)

        # Each level should have positive scalars
        for j in 1:nlevdecomp
            for c in 1:d.nc
                @test d.cf.t_scalar_col[c, j] > 0.0
                @test d.cf.w_scalar_col[c, j] >= 0.0
                @test d.cf.o_scalar_col[c, j] == 1.0
            end
        end

        # Decomp rate should decrease with depth (depth_scalar effect)
        for c in 1:d.nc
            for pool in 1:d.ndecomp_pools
                @test d.cf.decomp_k_col[c, 1, pool] >= d.cf.decomp_k_col[c, nlevdecomp, pool]
            end
        end
    end

    # ================================================================
    # Test 5: Temperature sensitivity
    # ================================================================
    @testset "Temperature sensitivity" begin
        d = make_decomp_bgc_data(; nlevdecomp=1)

        CLM.init_decomp_cascade_bgc!(
            d.bgc_state, d.cascade_con, d.params, d.cn_params;
            cellsand=d.cellsand,
            bounds=1:d.nc,
            nlevdecomp=d.nlevdecomp,
            ndecomp_pools_max=d.ndecomp_pools,
            ndecomp_cascade_transitions_max=d.ndecomp_cascade_transitions,
            use_fates=false)

        # Warm case
        t_warm = fill(CLM.TFRZ + 25.0, d.nc, 5)
        cf_warm = CLM.SoilBiogeochemCarbonFluxData()
        CLM.soil_bgc_carbon_flux_init!(cf_warm, d.nc, 5, d.ndecomp_pools,
                                        d.ndecomp_cascade_transitions)
        CLM.decomp_rate_constants_bgc!(
            cf_warm, d.bgc_state, d.params, d.cn_params, d.cascade_con;
            mask_bgc_soilc=d.mask_bgc_soilc,
            bounds=1:d.nc, nlevdecomp=d.nlevdecomp,
            t_soisno=t_warm, soilpsi=d.soilpsi,
            days_per_year=365.0, dt=1800.0,
            zsoi_vals=d.zsoi_vals, col_dz=d.col_dz)

        # Cold case
        t_cold = fill(CLM.TFRZ + 5.0, d.nc, 5)
        cf_cold = CLM.SoilBiogeochemCarbonFluxData()
        CLM.soil_bgc_carbon_flux_init!(cf_cold, d.nc, 5, d.ndecomp_pools,
                                        d.ndecomp_cascade_transitions)
        CLM.decomp_rate_constants_bgc!(
            cf_cold, d.bgc_state, d.params, d.cn_params, d.cascade_con;
            mask_bgc_soilc=d.mask_bgc_soilc,
            bounds=1:d.nc, nlevdecomp=d.nlevdecomp,
            t_soisno=t_cold, soilpsi=d.soilpsi,
            days_per_year=365.0, dt=1800.0,
            zsoi_vals=d.zsoi_vals, col_dz=d.col_dz)

        # Warm decomposition should be faster than cold
        for pool in 1:d.ndecomp_pools
            @test cf_warm.decomp_k_col[1, 1, pool] > cf_cold.decomp_k_col[1, 1, pool]
        end
    end

    # ================================================================
    # Test 6: Moisture sensitivity
    # ================================================================
    @testset "Moisture sensitivity" begin
        d = make_decomp_bgc_data(; nlevdecomp=1)

        CLM.init_decomp_cascade_bgc!(
            d.bgc_state, d.cascade_con, d.params, d.cn_params;
            cellsand=d.cellsand,
            bounds=1:d.nc,
            nlevdecomp=d.nlevdecomp,
            ndecomp_pools_max=d.ndecomp_pools,
            ndecomp_cascade_transitions_max=d.ndecomp_cascade_transitions,
            use_fates=false)

        # Wet case (soilpsi close to maxpsi)
        soilpsi_wet = fill(-0.2, d.nc, 5)
        cf_wet = CLM.SoilBiogeochemCarbonFluxData()
        CLM.soil_bgc_carbon_flux_init!(cf_wet, d.nc, 5, d.ndecomp_pools,
                                        d.ndecomp_cascade_transitions)
        CLM.decomp_rate_constants_bgc!(
            cf_wet, d.bgc_state, d.params, d.cn_params, d.cascade_con;
            mask_bgc_soilc=d.mask_bgc_soilc,
            bounds=1:d.nc, nlevdecomp=d.nlevdecomp,
            t_soisno=d.t_soisno, soilpsi=soilpsi_wet,
            days_per_year=365.0, dt=1800.0,
            zsoi_vals=d.zsoi_vals, col_dz=d.col_dz)

        # Very dry case (soilpsi near minpsi)
        soilpsi_dry = fill(-9.0, d.nc, 5)
        cf_dry = CLM.SoilBiogeochemCarbonFluxData()
        CLM.soil_bgc_carbon_flux_init!(cf_dry, d.nc, 5, d.ndecomp_pools,
                                        d.ndecomp_cascade_transitions)
        CLM.decomp_rate_constants_bgc!(
            cf_dry, d.bgc_state, d.params, d.cn_params, d.cascade_con;
            mask_bgc_soilc=d.mask_bgc_soilc,
            bounds=1:d.nc, nlevdecomp=d.nlevdecomp,
            t_soisno=d.t_soisno, soilpsi=soilpsi_dry,
            days_per_year=365.0, dt=1800.0,
            zsoi_vals=d.zsoi_vals, col_dz=d.col_dz)

        # Wet decomposition should be faster than dry
        for pool in 1:d.ndecomp_pools
            @test cf_wet.decomp_k_col[1, 1, pool] > cf_dry.decomp_k_col[1, 1, pool]
        end
    end

    # ================================================================
    # Test 7: Century temperature function
    # ================================================================
    @testset "Century temperature function" begin
        d = make_decomp_bgc_data(; nlevdecomp=1)

        d.bgc_state.use_century_tfunc = true
        d.bgc_state.normalize_q10_to_century_tfunc = false

        CLM.init_decomp_cascade_bgc!(
            d.bgc_state, d.cascade_con, d.params, d.cn_params;
            cellsand=d.cellsand,
            bounds=1:d.nc,
            nlevdecomp=d.nlevdecomp,
            ndecomp_pools_max=d.ndecomp_pools,
            ndecomp_cascade_transitions_max=d.ndecomp_cascade_transitions,
            use_fates=false)

        CLM.decomp_rate_constants_bgc!(
            d.cf, d.bgc_state, d.params, d.cn_params, d.cascade_con;
            mask_bgc_soilc=d.mask_bgc_soilc,
            bounds=1:d.nc, nlevdecomp=d.nlevdecomp,
            t_soisno=d.t_soisno, soilpsi=d.soilpsi,
            days_per_year=365.0, dt=1800.0,
            zsoi_vals=d.zsoi_vals, col_dz=d.col_dz)

        # Should still produce positive rates
        for c in 1:d.nc
            @test d.cf.t_scalar_col[c, 1] > 0.0
            for pool in 1:d.ndecomp_pools
                @test d.cf.decomp_k_col[c, 1, pool] > 0.0
            end
        end
    end

    # ================================================================
    # Test 8: Error on conflicting config
    # ================================================================
    @testset "Config conflict error" begin
        d = make_decomp_bgc_data(; nlevdecomp=1)

        d.bgc_state.use_century_tfunc = true
        d.bgc_state.normalize_q10_to_century_tfunc = true

        CLM.init_decomp_cascade_bgc!(
            d.bgc_state, d.cascade_con, d.params, d.cn_params;
            cellsand=d.cellsand,
            bounds=1:d.nc,
            nlevdecomp=d.nlevdecomp,
            ndecomp_pools_max=d.ndecomp_pools,
            ndecomp_cascade_transitions_max=d.ndecomp_cascade_transitions,
            use_fates=false)

        @test_throws ErrorException CLM.decomp_rate_constants_bgc!(
            d.cf, d.bgc_state, d.params, d.cn_params, d.cascade_con;
            mask_bgc_soilc=d.mask_bgc_soilc,
            bounds=1:d.nc, nlevdecomp=d.nlevdecomp,
            t_soisno=d.t_soisno, soilpsi=d.soilpsi,
            days_per_year=365.0, dt=1800.0,
            zsoi_vals=d.zsoi_vals, col_dz=d.col_dz)
    end

    # ================================================================
    # Test 9: Sand texture effect on path fractions
    # ================================================================
    @testset "Sand texture effect" begin
        nc = 2
        nlevdecomp = 1
        ndecomp_pools = 7
        ndecomp_cascade_transitions = 10

        params = CLM.DecompBGCParams(
            bgc_initial_Cstocks = fill(200.0, ndecomp_pools),
        )
        cn_params = CLM.CNSharedParamsData()

        # Sandy soil (90% sand)
        bgc_sandy = CLM.DecompBGCState()
        cc_sandy  = CLM.DecompCascadeConData()
        cellsand_sandy = fill(90.0, nc, 5)
        CLM.init_decomp_cascade_bgc!(
            bgc_sandy, cc_sandy, params, cn_params;
            cellsand=cellsand_sandy, bounds=1:nc,
            nlevdecomp=nlevdecomp, ndecomp_pools_max=ndecomp_pools,
            ndecomp_cascade_transitions_max=ndecomp_cascade_transitions)

        # Clay soil (10% sand)
        bgc_clay = CLM.DecompBGCState()
        cc_clay  = CLM.DecompCascadeConData()
        cellsand_clay = fill(10.0, nc, 5)
        CLM.init_decomp_cascade_bgc!(
            bgc_clay, cc_clay, params, cn_params;
            cellsand=cellsand_clay, bounds=1:nc,
            nlevdecomp=nlevdecomp, ndecomp_pools_max=ndecomp_pools,
            ndecomp_cascade_transitions_max=ndecomp_cascade_transitions)

        # Sand-dependent path fractions should differ
        @test bgc_sandy.f_s1s2[1, 1] != bgc_clay.f_s1s2[1, 1]
        @test bgc_sandy.rf_s1s2[1, 1] != bgc_clay.rf_s1s2[1, 1]

        # Path fractions should sum to 1
        @test bgc_sandy.f_s1s2[1, 1] + bgc_sandy.f_s1s3[1, 1] ≈ 1.0
        @test bgc_clay.f_s1s2[1, 1] + bgc_clay.f_s1s3[1, 1] ≈ 1.0
    end

end
