@testset "Decomp MIMICS" begin

    # ----------------------------------------------------------------
    # Helper: create default MIMICS parameters with realistic values
    # Based on Wieder et al. 2015 MIMICS model
    # ----------------------------------------------------------------
    function make_mimics_params(; ndecomp_pools_max=8)
        CLM.DecompMIMICSParams(
            mimics_nue_into_mic          = 0.85,
            mimics_desorpQ10             = 1.1,
            mimics_densdep               = 1.0,
            mimics_tau_mod_factor        = 2.0,
            mimics_tau_mod_min           = 0.8,
            mimics_tau_mod_max           = 1.2,
            mimics_ko_r                  = 6.0,
            mimics_ko_k                  = 6.0,
            mimics_cn_r                  = 7.0,
            mimics_cn_k                  = 10.0,
            mimics_cn_mod_num            = 0.4,
            mimics_t_soi_ref             = 25.0,
            mimics_initial_Cstocks_depth = 0.3,
            mimics_initial_Cstocks       = fill(100.0, ndecomp_pools_max),
            mimics_mge                   = [0.6, 0.2, 0.6, 0.6, 0.3, 0.6],
            mimics_vmod                  = [10.0, 3.0, 10.0, 10.0, 2.0, 10.0],
            mimics_vint                  = [5.47, 5.47, 5.47, 5.47, 5.47, 5.47],
            mimics_vslope                = [0.063, 0.063, 0.063, 0.063, 0.063, 0.063],
            mimics_kmod                  = [0.125, 0.5, 0.25, 0.5, 0.25, 0.167],
            mimics_kint                  = [3.19, 3.19, 3.19, 3.19, 3.19, 3.19],
            mimics_kslope                = [0.017, 0.017, 0.017, 0.017, 0.017, 0.017],
            mimics_fmet                  = [1.0, 0.85, 0.013, 25.0],
            mimics_p_scalar              = [2.0, -2.5],
            mimics_fphys_r               = [0.3, 1.3],
            mimics_fphys_k               = [0.2, 0.8],
            mimics_fchem_r               = [0.1, -3.0],
            mimics_fchem_k               = [0.3, -3.0],
            mimics_desorp                = [0.00015, 1.5],
            mimics_tau_r                 = [5.2e-4, -8.0],
            mimics_tau_k                 = [2.4e-4, -2.0],
        )
    end

    # ----------------------------------------------------------------
    # Helper: create minimal data structures for nc columns
    # ----------------------------------------------------------------
    function make_mimics_data(; nc=4, nlevdecomp=1, ndecomp_pools=8,
                                ndecomp_cascade_transitions=15)
        params = make_mimics_params(ndecomp_pools_max=ndecomp_pools)

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

        mimics_state = CLM.DecompMIMICSState()
        cascade_con  = CLM.DecompCascadeConData()

        # Soil clay values (col × nlevdecomp)
        cellclay = fill(30.0, nc, max(nlevdecomp, 5))

        # Soil temperature (col × nlev)
        t_soisno = fill(CLM.TFRZ + 15.0, nc, max(nlevdecomp, 5))

        # Soil water potential (col × nlev)
        soilpsi = fill(-1.0, nc, max(nlevdecomp, 5))

        # Column layer thicknesses
        col_dz = fill(0.1, nc, max(nlevdecomp, 5))

        # Carbon flux data
        cf = CLM.SoilBiogeochemCarbonFluxData()
        CLM.soil_bgc_carbon_flux_init!(cf, nc, max(nlevdecomp, 5), ndecomp_pools,
                                        ndecomp_cascade_transitions)

        # Decomposing C pools (gC/m3)
        decomp_cpools_vr = fill(10.0, nc, max(nlevdecomp, 5), ndecomp_pools)

        # Lignin:N ratio
        ligninNratioAvg = fill(10.0, nc)

        # Annual NPP at column level (gC/m2/s)
        annsum_npp_col = fill(500.0, nc)

        # Mask
        mask_bgc_soilc = trues(nc)

        return (params=params, cn_params=cn_params, mimics_state=mimics_state,
                cascade_con=cascade_con, cf=cf, cellclay=cellclay,
                t_soisno=t_soisno, soilpsi=soilpsi, col_dz=col_dz,
                decomp_cpools_vr=decomp_cpools_vr,
                ligninNratioAvg=ligninNratioAvg,
                annsum_npp_col=annsum_npp_col,
                mask_bgc_soilc=mask_bgc_soilc,
                nlevdecomp=nlevdecomp, ndecomp_pools=ndecomp_pools,
                ndecomp_cascade_transitions=ndecomp_cascade_transitions, nc=nc)
    end

    # ================================================================
    # Test 1: DecompMIMICSParams construction and read
    # ================================================================
    @testset "DecompMIMICSParams construction" begin
        params = CLM.DecompMIMICSParams()
        @test params.mimics_densdep == 1.0
        @test params.mimics_t_soi_ref == 25.0

        # Test read_params
        params2 = CLM.DecompMIMICSParams()
        CLM.decomp_mimics_read_params!(params2;
            mimics_initial_Cstocks_depth=0.3,
            mimics_initial_Cstocks=fill(100.0, 8),
            mimics_mge=[0.6, 0.2, 0.6, 0.6, 0.3, 0.6],
            mimics_vmod=[10.0, 3.0, 10.0, 10.0, 2.0, 10.0],
            mimics_vslope=[0.063, 0.063, 0.063, 0.063, 0.063, 0.063],
            mimics_vint=[5.47, 5.47, 5.47, 5.47, 5.47, 5.47],
            mimics_kmod=[0.125, 0.5, 0.25, 0.5, 0.25, 0.167],
            mimics_kslope=[0.017, 0.017, 0.017, 0.017, 0.017, 0.017],
            mimics_kint=[3.19, 3.19, 3.19, 3.19, 3.19, 3.19],
            mimics_p_scalar=[2.0, -2.5],
            mimics_desorp=[0.00015, 1.5],
            mimics_fphys_r=[0.3, 1.3],
            mimics_fphys_k=[0.2, 0.8],
            mimics_fmet=[1.0, 0.85, 0.013, 25.0],
            mimics_fchem_r=[0.1, -3.0],
            mimics_fchem_k=[0.3, -3.0],
            mimics_tau_r=[5.2e-4, -8.0],
            mimics_tau_k=[2.4e-4, -2.0],
            mimics_nue_into_mic=0.85,
            mimics_tau_mod_factor=2.0,
            mimics_tau_mod_min=0.8,
            mimics_tau_mod_max=1.2,
            mimics_ko_r=6.0,
            mimics_ko_k=6.0,
            mimics_densdep=1.0,
            mimics_desorpQ10=1.1,
            mimics_t_soi_ref=25.0,
            mimics_cn_mod_num=0.4,
            mimics_cn_r=7.0,
            mimics_cn_k=10.0)
        @test params2.mimics_nue_into_mic == 0.85
        @test params2.mimics_ko_r == 6.0
        @test length(params2.mimics_mge) == 6
    end

    # ================================================================
    # Test 2: init_decompcascade_mimics!
    # ================================================================
    @testset "init_decompcascade_mimics!" begin
        d = make_mimics_data()

        nue = CLM.init_decompcascade_mimics!(
            d.mimics_state, d.cascade_con, d.params, d.cn_params;
            cellclay=d.cellclay,
            bounds=1:d.nc,
            nlevdecomp=d.nlevdecomp,
            ndecomp_pools_max=d.ndecomp_pools,
            ndecomp_cascade_transitions_max=d.ndecomp_cascade_transitions,
            use_fates=false)

        # Check pool indices
        @test d.mimics_state.i_met_lit == 1
        @test d.mimics_state.i_str_lit == 2
        @test d.mimics_state.i_avl_som == 3
        @test d.mimics_state.i_chem_som == 4
        @test d.mimics_state.i_phys_som == 5
        @test d.mimics_state.i_cop_mic == 6
        @test d.mimics_state.i_oli_mic == 7

        # Check transition indices
        @test d.mimics_state.i_l1m1 == 1
        @test d.mimics_state.i_l1m2 == 2
        @test d.mimics_state.i_l2m1 == 3
        @test d.mimics_state.i_l2m2 == 4
        @test d.mimics_state.i_s1m1 == 5
        @test d.mimics_state.i_s1m2 == 6
        @test d.mimics_state.i_s2s1 == 7
        @test d.mimics_state.i_s3s1 == 8
        @test d.mimics_state.i_m1s1 == 9
        @test d.mimics_state.i_m1s2 == 10
        @test d.mimics_state.i_m1s3 == 11
        @test d.mimics_state.i_m2s1 == 12
        @test d.mimics_state.i_m2s2 == 13
        @test d.mimics_state.i_m2s3 == 14

        # Check cascade donor/receiver pools
        @test d.cascade_con.cascade_donor_pool[1] == 1    # L1 -> M1
        @test d.cascade_con.cascade_receiver_pool[1] == 6 # -> copiotrophic
        @test d.cascade_con.cascade_donor_pool[7] == 4    # S2 -> S1
        @test d.cascade_con.cascade_receiver_pool[7] == 3 # -> available SOM

        # CWD transition
        @test d.cascade_con.cascade_step_name[15] == "CWDL2"
        @test d.cascade_con.cascade_donor_pool[15] == 8   # CWD
        @test d.cascade_con.cascade_receiver_pool[15] == 2 # -> structural litter

        # Pool flags
        @test d.cascade_con.is_litter[1] == true   # met_lit
        @test d.cascade_con.is_litter[2] == true   # str_lit
        @test d.cascade_con.is_soil[3] == true     # avl_som
        @test d.cascade_con.is_soil[4] == true     # chem_som
        @test d.cascade_con.is_soil[5] == true     # phys_som
        @test d.cascade_con.is_cwd[8] == true      # CWD

        # Respiration fractions
        @test d.mimics_state.rf_l1m1 == 1.0 - 0.6  # 1 - mge[1]
        @test d.mimics_state.rf_l2m1 == 1.0 - 0.2  # 1 - mge[2]

        # NUE
        @test nue[1] == 0.85  # mimics_nue_into_mic
        @test nue[7] == 1.0   # S2S1

        # Spinup factors all 1.0 for normal run
        @test d.cascade_con.spinup_factor[1] == 1.0
        @test d.cascade_con.spinup_factor[6] == 1.0

        # Spatially-varying arrays should be allocated and positive
        @test size(d.mimics_state.desorp) == (d.nc, d.nlevdecomp)
        @test all(d.mimics_state.desorp .> 0)
        @test all(d.mimics_state.fphys_m1 .> 0)
        @test all(d.mimics_state.fphys_m2 .> 0)
        @test all(d.mimics_state.p_scalar .> 0)

        # Initial stocks
        @test d.cascade_con.initial_stock[1] == 100.0
        @test d.cascade_con.initial_stock_soildepth == 0.3
    end

    # ================================================================
    # Test 3: decomp_rates_mimics! basic functionality
    # ================================================================
    @testset "decomp_rates_mimics! basic" begin
        d = make_mimics_data(; nlevdecomp=1)

        CLM.init_decompcascade_mimics!(
            d.mimics_state, d.cascade_con, d.params, d.cn_params;
            cellclay=d.cellclay,
            bounds=1:d.nc,
            nlevdecomp=d.nlevdecomp,
            ndecomp_pools_max=d.ndecomp_pools,
            ndecomp_cascade_transitions_max=d.ndecomp_cascade_transitions,
            use_fates=false)

        CLM.decomp_rates_mimics!(
            d.cf, d.mimics_state, d.params, d.cn_params, d.cascade_con;
            mask_bgc_soilc=d.mask_bgc_soilc,
            bounds=1:d.nc,
            nlevdecomp=d.nlevdecomp,
            t_soisno=d.t_soisno,
            soilpsi=d.soilpsi,
            decomp_cpools_vr=d.decomp_cpools_vr,
            col_dz=d.col_dz,
            ligninNratioAvg=d.ligninNratioAvg,
            annsum_npp_col=d.annsum_npp_col,
            days_per_year=365.0,
            dt=1800.0)

        # w_scalar should be in [0, 1]
        for c in 1:d.nc
            @test d.cf.w_scalar_col[c, 1] >= 0.0
            @test d.cf.w_scalar_col[c, 1] <= 1.0
        end

        # o_scalar should be 1.0 (no anoxia)
        for c in 1:d.nc
            @test d.cf.o_scalar_col[c, 1] == 1.0
        end

        # decomp_k should be non-negative for all pools
        for c in 1:d.nc
            for pool in 1:d.ndecomp_pools
                @test d.cf.decomp_k_col[c, 1, pool] >= 0.0
            end
        end

        # Metabolic litter decomp_k should be positive
        @test d.cf.decomp_k_col[1, 1, 1] > 0.0

        # pathfrac for s2s1 and s3s1 should be 1.0
        @test d.cf.pathfrac_decomp_cascade_col[1, 1, d.mimics_state.i_s2s1] == 1.0
        @test d.cf.pathfrac_decomp_cascade_col[1, 1, d.mimics_state.i_s3s1] == 1.0

        # rf for microbial turnover should be 0.0
        @test d.cf.rf_decomp_cascade_col[1, 1, d.mimics_state.i_m1s1] == 0.0
        @test d.cf.rf_decomp_cascade_col[1, 1, d.mimics_state.i_m2s1] == 0.0

        # rf for litter->microbe should match state
        @test d.cf.rf_decomp_cascade_col[1, 1, d.mimics_state.i_l1m1] ≈ d.mimics_state.rf_l1m1

        # C:N ratios for microbes should be positive
        @test d.cf.cn_col[1, d.mimics_state.i_cop_mic] > 0.0
        @test d.cf.cn_col[1, d.mimics_state.i_oli_mic] > 0.0

        # pathfrac for l1m1 + l1m2 should sum to 1 (or both be 0)
        pf1 = d.cf.pathfrac_decomp_cascade_col[1, 1, d.mimics_state.i_l1m1]
        pf2 = d.cf.pathfrac_decomp_cascade_col[1, 1, d.mimics_state.i_l1m2]
        @test pf1 + pf2 ≈ 1.0 || (pf1 == 0.0 && pf2 == 0.0)
    end

    # ================================================================
    # Test 4: Temperature sensitivity
    # ================================================================
    @testset "Temperature sensitivity" begin
        d = make_mimics_data(; nlevdecomp=1)

        CLM.init_decompcascade_mimics!(
            d.mimics_state, d.cascade_con, d.params, d.cn_params;
            cellclay=d.cellclay,
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
        CLM.decomp_rates_mimics!(
            cf_warm, d.mimics_state, d.params, d.cn_params, d.cascade_con;
            mask_bgc_soilc=d.mask_bgc_soilc,
            bounds=1:d.nc, nlevdecomp=d.nlevdecomp,
            t_soisno=t_warm, soilpsi=d.soilpsi,
            decomp_cpools_vr=d.decomp_cpools_vr, col_dz=d.col_dz,
            ligninNratioAvg=d.ligninNratioAvg,
            annsum_npp_col=d.annsum_npp_col,
            days_per_year=365.0, dt=1800.0)

        # Cold case
        t_cold = fill(CLM.TFRZ + 5.0, d.nc, 5)
        cf_cold = CLM.SoilBiogeochemCarbonFluxData()
        CLM.soil_bgc_carbon_flux_init!(cf_cold, d.nc, 5, d.ndecomp_pools,
                                        d.ndecomp_cascade_transitions)
        CLM.decomp_rates_mimics!(
            cf_cold, d.mimics_state, d.params, d.cn_params, d.cascade_con;
            mask_bgc_soilc=d.mask_bgc_soilc,
            bounds=1:d.nc, nlevdecomp=d.nlevdecomp,
            t_soisno=t_cold, soilpsi=d.soilpsi,
            decomp_cpools_vr=d.decomp_cpools_vr, col_dz=d.col_dz,
            ligninNratioAvg=d.ligninNratioAvg,
            annsum_npp_col=d.annsum_npp_col,
            days_per_year=365.0, dt=1800.0)

        # Warm metabolic litter decomposition should be faster
        @test cf_warm.decomp_k_col[1, 1, 1] > cf_cold.decomp_k_col[1, 1, 1]

        # Warm desorption (phys_som) should be faster
        i_phys = d.mimics_state.i_phys_som
        @test cf_warm.decomp_k_col[1, 1, i_phys] > cf_cold.decomp_k_col[1, 1, i_phys]
    end

    # ================================================================
    # Test 5: Moisture sensitivity
    # ================================================================
    @testset "Moisture sensitivity" begin
        d = make_mimics_data(; nlevdecomp=1)

        CLM.init_decompcascade_mimics!(
            d.mimics_state, d.cascade_con, d.params, d.cn_params;
            cellclay=d.cellclay,
            bounds=1:d.nc,
            nlevdecomp=d.nlevdecomp,
            ndecomp_pools_max=d.ndecomp_pools,
            ndecomp_cascade_transitions_max=d.ndecomp_cascade_transitions,
            use_fates=false)

        # Wet case
        soilpsi_wet = fill(-0.2, d.nc, 5)
        cf_wet = CLM.SoilBiogeochemCarbonFluxData()
        CLM.soil_bgc_carbon_flux_init!(cf_wet, d.nc, 5, d.ndecomp_pools,
                                        d.ndecomp_cascade_transitions)
        CLM.decomp_rates_mimics!(
            cf_wet, d.mimics_state, d.params, d.cn_params, d.cascade_con;
            mask_bgc_soilc=d.mask_bgc_soilc,
            bounds=1:d.nc, nlevdecomp=d.nlevdecomp,
            t_soisno=d.t_soisno, soilpsi=soilpsi_wet,
            decomp_cpools_vr=d.decomp_cpools_vr, col_dz=d.col_dz,
            ligninNratioAvg=d.ligninNratioAvg,
            annsum_npp_col=d.annsum_npp_col,
            days_per_year=365.0, dt=1800.0)

        # Dry case
        soilpsi_dry = fill(-9.0, d.nc, 5)
        cf_dry = CLM.SoilBiogeochemCarbonFluxData()
        CLM.soil_bgc_carbon_flux_init!(cf_dry, d.nc, 5, d.ndecomp_pools,
                                        d.ndecomp_cascade_transitions)
        CLM.decomp_rates_mimics!(
            cf_dry, d.mimics_state, d.params, d.cn_params, d.cascade_con;
            mask_bgc_soilc=d.mask_bgc_soilc,
            bounds=1:d.nc, nlevdecomp=d.nlevdecomp,
            t_soisno=d.t_soisno, soilpsi=soilpsi_dry,
            decomp_cpools_vr=d.decomp_cpools_vr, col_dz=d.col_dz,
            ligninNratioAvg=d.ligninNratioAvg,
            annsum_npp_col=d.annsum_npp_col,
            days_per_year=365.0, dt=1800.0)

        # Wet metabolic litter decomposition should be faster
        @test cf_wet.decomp_k_col[1, 1, 1] > cf_dry.decomp_k_col[1, 1, 1]
    end

    # ================================================================
    # Test 6: Multi-level decomposition
    # ================================================================
    @testset "Multi-level decomposition" begin
        nlevdecomp = 5
        d = make_mimics_data(; nlevdecomp=nlevdecomp)

        CLM.init_decompcascade_mimics!(
            d.mimics_state, d.cascade_con, d.params, d.cn_params;
            cellclay=d.cellclay,
            bounds=1:d.nc,
            nlevdecomp=nlevdecomp,
            ndecomp_pools_max=d.ndecomp_pools,
            ndecomp_cascade_transitions_max=d.ndecomp_cascade_transitions,
            use_fates=false)

        CLM.decomp_rates_mimics!(
            d.cf, d.mimics_state, d.params, d.cn_params, d.cascade_con;
            mask_bgc_soilc=d.mask_bgc_soilc,
            bounds=1:d.nc,
            nlevdecomp=nlevdecomp,
            t_soisno=d.t_soisno,
            soilpsi=d.soilpsi,
            decomp_cpools_vr=d.decomp_cpools_vr,
            col_dz=d.col_dz,
            ligninNratioAvg=d.ligninNratioAvg,
            annsum_npp_col=d.annsum_npp_col,
            days_per_year=365.0,
            dt=1800.0)

        # Each level should have positive w_scalar
        for j in 1:nlevdecomp
            for c in 1:d.nc
                @test d.cf.w_scalar_col[c, j] >= 0.0
                @test d.cf.o_scalar_col[c, j] == 1.0
            end
        end

        # decomp_k should be non-negative at all levels
        for j in 1:nlevdecomp
            for c in 1:d.nc
                for pool in 1:d.ndecomp_pools
                    @test d.cf.decomp_k_col[c, j, pool] >= 0.0
                end
            end
        end
    end

    # ================================================================
    # Test 7: Clay texture effect on state variables
    # ================================================================
    @testset "Clay texture effect" begin
        nc = 2
        nlevdecomp = 1
        ndecomp_pools = 8
        ndecomp_cascade_transitions = 15

        params = make_mimics_params(ndecomp_pools_max=ndecomp_pools)
        cn_params = CLM.CNSharedParamsData()

        # High clay (60%)
        mimics_clay = CLM.DecompMIMICSState()
        cc_clay     = CLM.DecompCascadeConData()
        cellclay_high = fill(60.0, nc, 5)
        CLM.init_decompcascade_mimics!(
            mimics_clay, cc_clay, params, cn_params;
            cellclay=cellclay_high, bounds=1:nc,
            nlevdecomp=nlevdecomp, ndecomp_pools_max=ndecomp_pools,
            ndecomp_cascade_transitions_max=ndecomp_cascade_transitions)

        # Low clay (10%)
        mimics_sand = CLM.DecompMIMICSState()
        cc_sand     = CLM.DecompCascadeConData()
        cellclay_low = fill(10.0, nc, 5)
        CLM.init_decompcascade_mimics!(
            mimics_sand, cc_sand, params, cn_params;
            cellclay=cellclay_low, bounds=1:nc,
            nlevdecomp=nlevdecomp, ndecomp_pools_max=ndecomp_pools,
            ndecomp_cascade_transitions_max=ndecomp_cascade_transitions)

        # desorp should differ with clay content
        @test mimics_clay.desorp[1, 1] != mimics_sand.desorp[1, 1]

        # fphys should differ
        @test mimics_clay.fphys_m1[1, 1] != mimics_sand.fphys_m1[1, 1]
        @test mimics_clay.fphys_m2[1, 1] != mimics_sand.fphys_m2[1, 1]

        # p_scalar should differ
        @test mimics_clay.p_scalar[1, 1] != mimics_sand.p_scalar[1, 1]

        # Higher clay -> higher fphys (more physical protection)
        @test mimics_clay.fphys_m1[1, 1] > mimics_sand.fphys_m1[1, 1]
        @test mimics_clay.fphys_m2[1, 1] > mimics_sand.fphys_m2[1, 1]
    end

    # ================================================================
    # Test 8: Microbial C:N ratio depends on fmet
    # ================================================================
    @testset "Microbial C:N ratio" begin
        d = make_mimics_data(; nlevdecomp=1)

        CLM.init_decompcascade_mimics!(
            d.mimics_state, d.cascade_con, d.params, d.cn_params;
            cellclay=d.cellclay,
            bounds=1:d.nc,
            nlevdecomp=d.nlevdecomp,
            ndecomp_pools_max=d.ndecomp_pools,
            ndecomp_cascade_transitions_max=d.ndecomp_cascade_transitions,
            use_fates=false)

        # High lignin:N (low quality litter)
        lignin_high = fill(20.0, d.nc)
        cf1 = CLM.SoilBiogeochemCarbonFluxData()
        CLM.soil_bgc_carbon_flux_init!(cf1, d.nc, 5, d.ndecomp_pools,
                                        d.ndecomp_cascade_transitions)
        CLM.decomp_rates_mimics!(
            cf1, d.mimics_state, d.params, d.cn_params, d.cascade_con;
            mask_bgc_soilc=d.mask_bgc_soilc,
            bounds=1:d.nc, nlevdecomp=d.nlevdecomp,
            t_soisno=d.t_soisno, soilpsi=d.soilpsi,
            decomp_cpools_vr=d.decomp_cpools_vr, col_dz=d.col_dz,
            ligninNratioAvg=lignin_high,
            annsum_npp_col=d.annsum_npp_col,
            days_per_year=365.0, dt=1800.0)

        # Low lignin:N (high quality litter)
        lignin_low = fill(5.0, d.nc)
        cf2 = CLM.SoilBiogeochemCarbonFluxData()
        CLM.soil_bgc_carbon_flux_init!(cf2, d.nc, 5, d.ndecomp_pools,
                                        d.ndecomp_cascade_transitions)
        CLM.decomp_rates_mimics!(
            cf2, d.mimics_state, d.params, d.cn_params, d.cascade_con;
            mask_bgc_soilc=d.mask_bgc_soilc,
            bounds=1:d.nc, nlevdecomp=d.nlevdecomp,
            t_soisno=d.t_soisno, soilpsi=d.soilpsi,
            decomp_cpools_vr=d.decomp_cpools_vr, col_dz=d.col_dz,
            ligninNratioAvg=lignin_low,
            annsum_npp_col=d.annsum_npp_col,
            days_per_year=365.0, dt=1800.0)

        # C:N ratios should differ based on lignin:N
        cn_cop_high = cf1.cn_col[1, d.mimics_state.i_cop_mic]
        cn_cop_low  = cf2.cn_col[1, d.mimics_state.i_cop_mic]
        @test cn_cop_high != cn_cop_low
        @test cn_cop_high > 0.0
        @test cn_cop_low > 0.0
    end

end
