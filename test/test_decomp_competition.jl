@testset "SoilBiogeochemCompetitionMod" begin

    # -----------------------------------------------------------------------
    # Helper: build minimal data structures for competition tests
    # -----------------------------------------------------------------------
    function make_competition_data(;
            nc=3, nlevdecomp=2, ntrans=2,
            sminn_vr_val=10.0,
            smin_nh4_vr_val=5.0,
            smin_no3_vr_val=5.0,
            plant_ndemand_val=1.0,
            potential_immob_vr_val=0.5,
            pot_f_nit_vr_val=0.1,
            pot_f_denit_vr_val=0.05,
            n2_n2o_ratio_denit_vr_val=1.0,
            nfixation_prof_val=0.5)

        dzsoi_decomp = fill(1.0, nlevdecomp)

        # SoilBiogeochemStateData
        st = CLM.SoilBiogeochemStateData()
        st.fpg_col          = zeros(nc)
        st.fpi_col          = zeros(nc)
        st.fpi_vr_col       = zeros(nc, nlevdecomp)
        st.nfixation_prof_col = fill(nfixation_prof_val, nc, nlevdecomp)
        st.plant_ndemand_col  = fill(plant_ndemand_val, nc)

        # SoilBiogeochemNitrogenStateData
        ns = CLM.SoilBiogeochemNitrogenStateData()
        ns.sminn_vr_col     = fill(sminn_vr_val, nc, nlevdecomp)
        ns.smin_nh4_vr_col  = fill(smin_nh4_vr_val, nc, nlevdecomp)
        ns.smin_no3_vr_col  = fill(smin_no3_vr_val, nc, nlevdecomp)

        # SoilBiogeochemNitrogenFluxData
        nf = CLM.SoilBiogeochemNitrogenFluxData()
        nf.potential_immob_vr_col       = fill(potential_immob_vr_val, nc, nlevdecomp)
        nf.actual_immob_vr_col          = zeros(nc, nlevdecomp)
        nf.sminn_to_plant_vr_col        = zeros(nc, nlevdecomp)
        nf.sminn_to_plant_col           = zeros(nc)
        nf.actual_immob_col             = zeros(nc)
        nf.potential_immob_col          = zeros(nc)
        nf.supplement_to_sminn_vr_col   = zeros(nc, nlevdecomp)
        nf.sminn_to_denit_excess_vr_col = zeros(nc, nlevdecomp)
        nf.actual_immob_no3_vr_col      = zeros(nc, nlevdecomp)
        nf.actual_immob_nh4_vr_col      = zeros(nc, nlevdecomp)
        nf.smin_no3_to_plant_vr_col     = zeros(nc, nlevdecomp)
        nf.smin_nh4_to_plant_vr_col     = zeros(nc, nlevdecomp)
        nf.pot_f_nit_vr_col             = fill(pot_f_nit_vr_val, nc, nlevdecomp)
        nf.pot_f_denit_vr_col           = fill(pot_f_denit_vr_val, nc, nlevdecomp)
        nf.f_nit_vr_col                 = zeros(nc, nlevdecomp)
        nf.f_denit_vr_col               = zeros(nc, nlevdecomp)
        nf.n2_n2o_ratio_denit_vr_col    = fill(n2_n2o_ratio_denit_vr_val, nc, nlevdecomp)
        nf.f_n2o_denit_vr_col           = zeros(nc, nlevdecomp)
        nf.f_n2o_nit_vr_col             = zeros(nc, nlevdecomp)
        nf.sminn_to_plant_fun_vr_col    = zeros(nc, nlevdecomp)
        nf.sminn_to_plant_fun_no3_vr_col = zeros(nc, nlevdecomp)
        nf.sminn_to_plant_fun_nh4_vr_col = zeros(nc, nlevdecomp)

        # SoilBiogeochemCarbonFluxData
        cf = CLM.SoilBiogeochemCarbonFluxData()
        cf.c_overflow_vr = zeros(nc, nlevdecomp, ntrans)

        mask_bgc_soilc = trues(nc)
        bounds = 1:nc

        pmnf_decomp_cascade = zeros(nc, nlevdecomp, ntrans)
        p_decomp_cn_gain    = zeros(nc, nlevdecomp, ntrans)
        cascade_receiver_pool = ones(Int, ntrans)

        return st, nf, cf, ns, mask_bgc_soilc, bounds, dzsoi_decomp,
               pmnf_decomp_cascade, p_decomp_cn_gain, cascade_receiver_pool
    end

    # -----------------------------------------------------------------------
    @testset "SoilBGCCompetitionParams default construction" begin
        p = CLM.SoilBGCCompetitionParams()
        @test p.bdnr == 0.5
        @test p.compet_plant_no3 == 1.0
        @test p.compet_plant_nh4 == 1.0
        @test p.compet_decomp_no3 == 1.0
        @test p.compet_decomp_nh4 == 1.0
        @test p.compet_denit == 1.0
        @test p.compet_nit == 1.0
    end

    # -----------------------------------------------------------------------
    @testset "SoilBGCCompetitionState default construction" begin
        s = CLM.SoilBGCCompetitionState()
        @test s.dt == 1800.0
        @test s.bdnr == 0.0
        @test s.carbon_only == false
    end

    # -----------------------------------------------------------------------
    @testset "soil_bgc_competition_params_read!" begin
        p = CLM.SoilBGCCompetitionParams()
        CLM.soil_bgc_competition_params_read!(p;
            bdnr=0.3, compet_plant_no3=2.0, compet_plant_nh4=1.5,
            compet_decomp_no3=0.8, compet_decomp_nh4=0.9,
            compet_denit=0.7, compet_nit=0.6)
        @test p.bdnr == 0.3
        @test p.compet_plant_no3 == 2.0
        @test p.compet_plant_nh4 == 1.5
        @test p.compet_decomp_no3 == 0.8
        @test p.compet_decomp_nh4 == 0.9
        @test p.compet_denit == 0.7
        @test p.compet_nit == 0.6
    end

    # -----------------------------------------------------------------------
    @testset "soil_bgc_competition_init!" begin
        p = CLM.SoilBGCCompetitionParams(bdnr=0.5)
        s = CLM.SoilBGCCompetitionState()

        # NONE mode
        CLM.soil_bgc_competition_init!(s, p; dt=1800.0, suplnitro="NONE")
        @test s.dt == 1800.0
        @test s.bdnr ≈ 0.5 * (1800.0 / CLM.SECSPDAY)
        @test s.carbon_only == false

        # ALL mode
        CLM.soil_bgc_competition_init!(s, p; dt=3600.0, suplnitro="ALL")
        @test s.dt == 3600.0
        @test s.bdnr ≈ 0.5 * (3600.0 / CLM.SECSPDAY)
        @test s.carbon_only == true

        # Invalid mode
        @test_throws ErrorException CLM.soil_bgc_competition_init!(s, p; dt=1800.0, suplnitro="INVALID")
    end

    # -----------------------------------------------------------------------
    @testset "supplemental N constants" begin
        @test CLM.SUPLN_ALL == "ALL"
        @test CLM.SUPLN_NONE == "NONE"
    end

    # -----------------------------------------------------------------------
    @testset "non-nitrif: N not limiting" begin
        # Need sum_ndemand_vr * dt < sminn_vr
        # nuptake_prof = sminn_vr / (sminn_vr * sum(dzsoi)) = 1/sum(dzsoi) = 1.0
        # sum_ndemand = plant_ndemand * 1.0 + potential_immob = 0.5 + 0.2 = 0.7
        # Need 0.7 * dt < sminn_vr => sminn_vr > 1260 for dt=1800
        nc = 2; nlevdecomp = 1; ntrans = 1
        st, nf, cf, ns, mask, bounds, dzsoi, pmnf, pcng, crp =
            make_competition_data(nc=nc, nlevdecomp=nlevdecomp, ntrans=ntrans,
                sminn_vr_val=10000.0,      # must exceed (plant_ndemand + potential_immob) * dt
                plant_ndemand_val=0.5,     # low plant demand
                potential_immob_vr_val=0.2) # low immob demand

        params = CLM.SoilBGCCompetitionParams()
        state = CLM.SoilBGCCompetitionState(dt=1800.0, bdnr=0.01)

        CLM.soil_bgc_competition!(st, nf, cf, ns, state, params;
            mask_bgc_soilc=mask, bounds=bounds,
            nlevdecomp=nlevdecomp, ndecomp_cascade_transitions=ntrans,
            dzsoi_decomp=dzsoi, pmnf_decomp_cascade=pmnf,
            p_decomp_cn_gain=pcng, cascade_receiver_pool=crp,
            use_nitrif_denitrif=false)

        # fpi should be 1.0 (all immobilization satisfied)
        @test st.fpi_col[1] ≈ 1.0
        @test st.fpi_col[2] ≈ 1.0

        # fpi_vr should be 1.0
        @test st.fpi_vr_col[1, 1] ≈ 1.0

        # actual immob equals potential immob
        @test nf.actual_immob_vr_col[1, 1] ≈ 0.2

        # fpg should be 1.0
        @test st.fpg_col[1] ≈ 1.0
    end

    # -----------------------------------------------------------------------
    @testset "non-nitrif: N limiting (competition)" begin
        # With low mineral N and high demand, N should be partitioned
        nc = 1; nlevdecomp = 1; ntrans = 1
        st, nf, cf, ns, mask, bounds, dzsoi, pmnf, pcng, crp =
            make_competition_data(nc=nc, nlevdecomp=nlevdecomp, ntrans=ntrans,
                sminn_vr_val=0.001,           # very low mineral N
                plant_ndemand_val=10.0,       # high plant demand
                potential_immob_vr_val=5.0)   # high immob demand

        params = CLM.SoilBGCCompetitionParams()
        state = CLM.SoilBGCCompetitionState(dt=1800.0, bdnr=0.01)

        CLM.soil_bgc_competition!(st, nf, cf, ns, state, params;
            mask_bgc_soilc=mask, bounds=bounds,
            nlevdecomp=nlevdecomp, ndecomp_cascade_transitions=ntrans,
            dzsoi_decomp=dzsoi, pmnf_decomp_cascade=pmnf,
            p_decomp_cn_gain=pcng, cascade_receiver_pool=crp,
            use_nitrif_denitrif=false)

        # fpi should be < 1.0 (N limited)
        @test st.fpi_col[1] < 1.0

        # actual immob should be < potential immob
        @test nf.actual_immob_vr_col[1, 1] < 5.0
        @test nf.actual_immob_vr_col[1, 1] > 0.0

        # fpi_vr should be between 0 and 1
        @test 0.0 < st.fpi_vr_col[1, 1] < 1.0

        # fpg should be < 1
        @test st.fpg_col[1] < 1.0
    end

    # -----------------------------------------------------------------------
    @testset "non-nitrif: carbon-only mode supplements N" begin
        nc = 1; nlevdecomp = 1; ntrans = 1
        st, nf, cf, ns, mask, bounds, dzsoi, pmnf, pcng, crp =
            make_competition_data(nc=nc, nlevdecomp=nlevdecomp, ntrans=ntrans,
                sminn_vr_val=0.001,           # very low mineral N
                plant_ndemand_val=10.0,       # high demand
                potential_immob_vr_val=5.0)   # high immob demand

        params = CLM.SoilBGCCompetitionParams()
        state = CLM.SoilBGCCompetitionState(dt=1800.0, bdnr=0.01)

        CLM.soil_bgc_competition!(st, nf, cf, ns, state, params;
            mask_bgc_soilc=mask, bounds=bounds,
            nlevdecomp=nlevdecomp, ndecomp_cascade_transitions=ntrans,
            dzsoi_decomp=dzsoi, pmnf_decomp_cascade=pmnf,
            p_decomp_cn_gain=pcng, cascade_receiver_pool=crp,
            use_nitrif_denitrif=false,
            carbon_only=true)

        # With carbon_only, fpi should be 1.0
        @test st.fpi_vr_col[1, 1] ≈ 1.0
        @test st.fpi_col[1] ≈ 1.0
        @test st.fpg_col[1] ≈ 1.0

        # supplement should be positive
        @test nf.supplement_to_sminn_vr_col[1, 1] > 0.0
    end

    # -----------------------------------------------------------------------
    @testset "non-nitrif: masked columns are skipped" begin
        nc = 3; nlevdecomp = 1; ntrans = 1
        st, nf, cf, ns, mask, bounds, dzsoi, pmnf, pcng, crp =
            make_competition_data(nc=nc, nlevdecomp=nlevdecomp, ntrans=ntrans,
                sminn_vr_val=10000.0, plant_ndemand_val=1.0,
                potential_immob_vr_val=0.5)

        # Mask out column 2
        mask[2] = false

        params = CLM.SoilBGCCompetitionParams()
        state = CLM.SoilBGCCompetitionState(dt=1800.0, bdnr=0.01)

        CLM.soil_bgc_competition!(st, nf, cf, ns, state, params;
            mask_bgc_soilc=mask, bounds=bounds,
            nlevdecomp=nlevdecomp, ndecomp_cascade_transitions=ntrans,
            dzsoi_decomp=dzsoi, pmnf_decomp_cascade=pmnf,
            p_decomp_cn_gain=pcng, cascade_receiver_pool=crp,
            use_nitrif_denitrif=false)

        # Column 1 and 3 should have fpg=1 but column 2 should be untouched (0)
        @test st.fpg_col[1] ≈ 1.0
        @test st.fpg_col[2] ≈ 0.0  # not processed
        @test st.fpg_col[3] ≈ 1.0
    end

    # -----------------------------------------------------------------------
    @testset "non-nitrif: excess denitrification" begin
        nc = 1; nlevdecomp = 1; ntrans = 1
        st, nf, cf, ns, mask, bounds, dzsoi, pmnf, pcng, crp =
            make_competition_data(nc=nc, nlevdecomp=nlevdecomp, ntrans=ntrans,
                sminn_vr_val=1000.0,       # very high mineral N
                plant_ndemand_val=0.01,    # very low demand
                potential_immob_vr_val=0.01)

        params = CLM.SoilBGCCompetitionParams(bdnr=0.5)
        state = CLM.SoilBGCCompetitionState()
        CLM.soil_bgc_competition_init!(state, params; dt=1800.0)

        CLM.soil_bgc_competition!(st, nf, cf, ns, state, params;
            mask_bgc_soilc=mask, bounds=bounds,
            nlevdecomp=nlevdecomp, ndecomp_cascade_transitions=ntrans,
            dzsoi_decomp=dzsoi, pmnf_decomp_cascade=pmnf,
            p_decomp_cn_gain=pcng, cascade_receiver_pool=crp,
            use_nitrif_denitrif=false)

        # With excess N, there should be some denitrification
        @test nf.sminn_to_denit_excess_vr_col[1, 1] > 0.0
    end

    # -----------------------------------------------------------------------
    @testset "nitrif_denitrif: N not limiting" begin
        # NH4 demand = plant*nuptake + pot_immob + pot_nit = 0.1 + 0.05 + 0.01 = 0.16
        # Need 0.16 * 1800 = 288 < smin_nh4_vr => smin_nh4 > 288
        # NO3 residual demand = 0 + 0 + 0.005 = 0.005
        # Need 0.005 * 1800 = 9 < smin_no3_vr
        nc = 1; nlevdecomp = 1; ntrans = 1
        st, nf, cf, ns, mask, bounds, dzsoi, pmnf, pcng, crp =
            make_competition_data(nc=nc, nlevdecomp=nlevdecomp, ntrans=ntrans,
                sminn_vr_val=10000.0,
                smin_nh4_vr_val=5000.0,
                smin_no3_vr_val=5000.0,
                plant_ndemand_val=0.1,
                potential_immob_vr_val=0.05,
                pot_f_nit_vr_val=0.01,
                pot_f_denit_vr_val=0.005)

        params = CLM.SoilBGCCompetitionParams()
        state = CLM.SoilBGCCompetitionState(dt=1800.0, bdnr=0.01)

        CLM.soil_bgc_competition!(st, nf, cf, ns, state, params;
            mask_bgc_soilc=mask, bounds=bounds,
            nlevdecomp=nlevdecomp, ndecomp_cascade_transitions=ntrans,
            dzsoi_decomp=dzsoi, pmnf_decomp_cascade=pmnf,
            p_decomp_cn_gain=pcng, cascade_receiver_pool=crp,
            use_nitrif_denitrif=true)

        # fpi should be 1.0 (all demand satisfied)
        @test st.fpi_col[1] ≈ 1.0

        # actual immob should equal potential
        @test nf.actual_immob_vr_col[1, 1] ≈ 0.05

        # fpg should be 1.0
        @test st.fpg_col[1] ≈ 1.0

        # nitrification should equal potential
        @test nf.f_nit_vr_col[1, 1] ≈ 0.01

        # denitrification should equal potential
        @test nf.f_denit_vr_col[1, 1] ≈ 0.005

        # N2O from nitrification
        @test nf.f_n2o_nit_vr_col[1, 1] ≈ 0.01 * CLM.NITRIF_N2O_LOSS_FRAC

        # N2O from denitrification (ratio is 1.0)
        @test nf.f_n2o_denit_vr_col[1, 1] ≈ 0.005 / (1.0 + 1.0)
    end

    # -----------------------------------------------------------------------
    @testset "nitrif_denitrif: NH4 limiting" begin
        nc = 1; nlevdecomp = 1; ntrans = 1
        st, nf, cf, ns, mask, bounds, dzsoi, pmnf, pcng, crp =
            make_competition_data(nc=nc, nlevdecomp=nlevdecomp, ntrans=ntrans,
                sminn_vr_val=0.001,
                smin_nh4_vr_val=0.0001,    # very low NH4
                smin_no3_vr_val=50.0,      # lots of NO3
                plant_ndemand_val=5.0,
                potential_immob_vr_val=2.0,
                pot_f_nit_vr_val=0.5,
                pot_f_denit_vr_val=0.1)

        params = CLM.SoilBGCCompetitionParams()
        state = CLM.SoilBGCCompetitionState(dt=1800.0, bdnr=0.01)

        CLM.soil_bgc_competition!(st, nf, cf, ns, state, params;
            mask_bgc_soilc=mask, bounds=bounds,
            nlevdecomp=nlevdecomp, ndecomp_cascade_transitions=ntrans,
            dzsoi_decomp=dzsoi, pmnf_decomp_cascade=pmnf,
            p_decomp_cn_gain=pcng, cascade_receiver_pool=crp,
            use_nitrif_denitrif=true)

        # NH4 immobilization should be less than potential
        @test nf.actual_immob_nh4_vr_col[1, 1] < 2.0

        # NO3 should pick up the slack
        @test nf.actual_immob_no3_vr_col[1, 1] > 0.0

        # Total actual immob should be sum of nh4 + no3
        @test nf.actual_immob_vr_col[1, 1] ≈
            nf.actual_immob_nh4_vr_col[1, 1] + nf.actual_immob_no3_vr_col[1, 1]

        # Plant uptake should be sum of nh4 + no3
        @test nf.sminn_to_plant_vr_col[1, 1] ≈
            nf.smin_nh4_to_plant_vr_col[1, 1] + nf.smin_no3_to_plant_vr_col[1, 1]
    end

    # -----------------------------------------------------------------------
    @testset "nitrif_denitrif: carbon-only supplement" begin
        nc = 1; nlevdecomp = 1; ntrans = 1
        st, nf, cf, ns, mask, bounds, dzsoi, pmnf, pcng, crp =
            make_competition_data(nc=nc, nlevdecomp=nlevdecomp, ntrans=ntrans,
                sminn_vr_val=0.0001,
                smin_nh4_vr_val=0.00005,
                smin_no3_vr_val=0.00005,
                plant_ndemand_val=10.0,
                potential_immob_vr_val=5.0,
                pot_f_nit_vr_val=0.1,
                pot_f_denit_vr_val=0.05)

        params = CLM.SoilBGCCompetitionParams()
        state = CLM.SoilBGCCompetitionState(dt=1800.0, bdnr=0.01)

        CLM.soil_bgc_competition!(st, nf, cf, ns, state, params;
            mask_bgc_soilc=mask, bounds=bounds,
            nlevdecomp=nlevdecomp, ndecomp_cascade_transitions=ntrans,
            dzsoi_decomp=dzsoi, pmnf_decomp_cascade=pmnf,
            p_decomp_cn_gain=pcng, cascade_receiver_pool=crp,
            use_nitrif_denitrif=true,
            carbon_only=true)

        # Supplement should be positive
        @test nf.supplement_to_sminn_vr_col[1, 1] > 0.0

        # fpi should be 1.0
        @test st.fpi_vr_col[1, 1] ≈ 1.0
    end

    # -----------------------------------------------------------------------
    @testset "nitrif_denitrif: c_overflow_vr zero when not MIMICS" begin
        nc = 1; nlevdecomp = 1; ntrans = 1
        st, nf, cf, ns, mask, bounds, dzsoi, pmnf, pcng, crp =
            make_competition_data(nc=nc, nlevdecomp=nlevdecomp, ntrans=ntrans)

        params = CLM.SoilBGCCompetitionParams()
        state = CLM.SoilBGCCompetitionState(dt=1800.0, bdnr=0.01)

        CLM.soil_bgc_competition!(st, nf, cf, ns, state, params;
            mask_bgc_soilc=mask, bounds=bounds,
            nlevdecomp=nlevdecomp, ndecomp_cascade_transitions=ntrans,
            dzsoi_decomp=dzsoi, pmnf_decomp_cascade=pmnf,
            p_decomp_cn_gain=pcng, cascade_receiver_pool=crp,
            use_nitrif_denitrif=true,
            mimics_decomp=false)

        @test all(cf.c_overflow_vr .== 0.0)
    end

    # -----------------------------------------------------------------------
    @testset "nitrif_denitrif: multiple levels" begin
        nc = 2; nlevdecomp = 3; ntrans = 1
        st, nf, cf, ns, mask, bounds, dzsoi, pmnf, pcng, crp =
            make_competition_data(nc=nc, nlevdecomp=nlevdecomp, ntrans=ntrans,
                sminn_vr_val=50.0,
                smin_nh4_vr_val=25.0,
                smin_no3_vr_val=25.0,
                plant_ndemand_val=1.0,
                potential_immob_vr_val=0.3,
                pot_f_nit_vr_val=0.05,
                pot_f_denit_vr_val=0.02)

        params = CLM.SoilBGCCompetitionParams()
        state = CLM.SoilBGCCompetitionState(dt=1800.0, bdnr=0.01)

        CLM.soil_bgc_competition!(st, nf, cf, ns, state, params;
            mask_bgc_soilc=mask, bounds=bounds,
            nlevdecomp=nlevdecomp, ndecomp_cascade_transitions=ntrans,
            dzsoi_decomp=dzsoi, pmnf_decomp_cascade=pmnf,
            p_decomp_cn_gain=pcng, cascade_receiver_pool=crp,
            use_nitrif_denitrif=true)

        # All levels should have consistent results for uniform input
        for j in 1:nlevdecomp
            @test nf.actual_immob_vr_col[1, j] ≈ nf.actual_immob_vr_col[2, j]
            @test nf.sminn_to_plant_vr_col[1, j] ≈ nf.sminn_to_plant_vr_col[2, j]
        end

        # fpg and fpi should be the same for both columns
        @test st.fpg_col[1] ≈ st.fpg_col[2]
        @test st.fpi_col[1] ≈ st.fpi_col[2]
    end

    # -----------------------------------------------------------------------
    @testset "nitrif_denitrif: competitiveness factors affect partitioning" begin
        nc = 1; nlevdecomp = 1; ntrans = 1

        # Run with equal competitiveness
        st1, nf1, cf1, ns1, mask1, bounds1, dzsoi1, pmnf1, pcng1, crp1 =
            make_competition_data(nc=nc, nlevdecomp=nlevdecomp, ntrans=ntrans,
                sminn_vr_val=0.01,
                smin_nh4_vr_val=0.001,
                smin_no3_vr_val=0.001,
                plant_ndemand_val=2.0,
                potential_immob_vr_val=1.0,
                pot_f_nit_vr_val=0.1,
                pot_f_denit_vr_val=0.05)

        params1 = CLM.SoilBGCCompetitionParams(
            compet_plant_no3=1.0, compet_plant_nh4=1.0,
            compet_decomp_no3=1.0, compet_decomp_nh4=1.0,
            compet_denit=1.0, compet_nit=1.0)
        state1 = CLM.SoilBGCCompetitionState(dt=1800.0, bdnr=0.01)

        CLM.soil_bgc_competition!(st1, nf1, cf1, ns1, state1, params1;
            mask_bgc_soilc=mask1, bounds=bounds1,
            nlevdecomp=nlevdecomp, ndecomp_cascade_transitions=ntrans,
            dzsoi_decomp=dzsoi1, pmnf_decomp_cascade=pmnf1,
            p_decomp_cn_gain=pcng1, cascade_receiver_pool=crp1,
            use_nitrif_denitrif=true)

        plant_uptake_1 = nf1.sminn_to_plant_vr_col[1, 1]

        # Run with plants having higher competitiveness
        st2, nf2, cf2, ns2, mask2, bounds2, dzsoi2, pmnf2, pcng2, crp2 =
            make_competition_data(nc=nc, nlevdecomp=nlevdecomp, ntrans=ntrans,
                sminn_vr_val=0.01,
                smin_nh4_vr_val=0.001,
                smin_no3_vr_val=0.001,
                plant_ndemand_val=2.0,
                potential_immob_vr_val=1.0,
                pot_f_nit_vr_val=0.1,
                pot_f_denit_vr_val=0.05)

        params2 = CLM.SoilBGCCompetitionParams(
            compet_plant_no3=10.0, compet_plant_nh4=10.0,
            compet_decomp_no3=1.0, compet_decomp_nh4=1.0,
            compet_denit=1.0, compet_nit=1.0)
        state2 = CLM.SoilBGCCompetitionState(dt=1800.0, bdnr=0.01)

        CLM.soil_bgc_competition!(st2, nf2, cf2, ns2, state2, params2;
            mask_bgc_soilc=mask2, bounds=bounds2,
            nlevdecomp=nlevdecomp, ndecomp_cascade_transitions=ntrans,
            dzsoi_decomp=dzsoi2, pmnf_decomp_cascade=pmnf2,
            p_decomp_cn_gain=pcng2, cascade_receiver_pool=crp2,
            use_nitrif_denitrif=true)

        plant_uptake_2 = nf2.sminn_to_plant_vr_col[1, 1]

        # With higher plant competitiveness, plants should get more N
        @test plant_uptake_2 > plant_uptake_1
    end

    # -----------------------------------------------------------------------
    @testset "N conservation: non-nitrif pathway" begin
        # Verify that total N supply >= total N demand in N-limited case
        nc = 1; nlevdecomp = 2; ntrans = 1
        st, nf, cf, ns, mask, bounds, dzsoi, pmnf, pcng, crp =
            make_competition_data(nc=nc, nlevdecomp=nlevdecomp, ntrans=ntrans,
                sminn_vr_val=1.0,
                plant_ndemand_val=2.0,
                potential_immob_vr_val=1.0)

        params = CLM.SoilBGCCompetitionParams(bdnr=0.0)  # no excess denitrif for clean check
        state = CLM.SoilBGCCompetitionState(dt=1800.0, bdnr=0.0)

        CLM.soil_bgc_competition!(st, nf, cf, ns, state, params;
            mask_bgc_soilc=mask, bounds=bounds,
            nlevdecomp=nlevdecomp, ndecomp_cascade_transitions=ntrans,
            dzsoi_decomp=dzsoi, pmnf_decomp_cascade=pmnf,
            p_decomp_cn_gain=pcng, cascade_receiver_pool=crp,
            use_nitrif_denitrif=false)

        # Total uptake (plant + immob) per level should not exceed sminn/dt
        for j in 1:nlevdecomp
            total_flux = nf.sminn_to_plant_vr_col[1, j] + nf.actual_immob_vr_col[1, j]
            @test total_flux <= ns.sminn_vr_col[1, j] / state.dt + 1e-12
        end
    end

    # -----------------------------------------------------------------------
    @testset "N conservation: nitrif_denitrif pathway" begin
        nc = 1; nlevdecomp = 1; ntrans = 1
        st, nf, cf, ns, mask, bounds, dzsoi, pmnf, pcng, crp =
            make_competition_data(nc=nc, nlevdecomp=nlevdecomp, ntrans=ntrans,
                sminn_vr_val=0.5,
                smin_nh4_vr_val=0.3,
                smin_no3_vr_val=0.2,
                plant_ndemand_val=5.0,
                potential_immob_vr_val=3.0,
                pot_f_nit_vr_val=0.5,
                pot_f_denit_vr_val=0.2)

        params = CLM.SoilBGCCompetitionParams()
        state = CLM.SoilBGCCompetitionState(dt=1800.0, bdnr=0.01)

        CLM.soil_bgc_competition!(st, nf, cf, ns, state, params;
            mask_bgc_soilc=mask, bounds=bounds,
            nlevdecomp=nlevdecomp, ndecomp_cascade_transitions=ntrans,
            dzsoi_decomp=dzsoi, pmnf_decomp_cascade=pmnf,
            p_decomp_cn_gain=pcng, cascade_receiver_pool=crp,
            use_nitrif_denitrif=true)

        j = 1
        # NH4 fluxes should not exceed nh4/dt
        nh4_total = nf.actual_immob_nh4_vr_col[1, j] +
            nf.smin_nh4_to_plant_vr_col[1, j] + nf.f_nit_vr_col[1, j]
        @test nh4_total <= ns.smin_nh4_vr_col[1, j] / state.dt + 1e-12

        # NO3 fluxes should not exceed no3/dt
        no3_total = nf.actual_immob_no3_vr_col[1, j] +
            nf.smin_no3_to_plant_vr_col[1, j] + nf.f_denit_vr_col[1, j]
        @test no3_total <= ns.smin_no3_vr_col[1, j] / state.dt + 1e-12
    end

    # -----------------------------------------------------------------------
    @testset "zero demand produces fpi=1 and fpg=1" begin
        nc = 1; nlevdecomp = 1; ntrans = 1
        st, nf, cf, ns, mask, bounds, dzsoi, pmnf, pcng, crp =
            make_competition_data(nc=nc, nlevdecomp=nlevdecomp, ntrans=ntrans,
                sminn_vr_val=10.0,
                plant_ndemand_val=0.0,
                potential_immob_vr_val=0.0)

        params = CLM.SoilBGCCompetitionParams()
        state = CLM.SoilBGCCompetitionState(dt=1800.0, bdnr=0.01)

        # Test non-nitrif path
        CLM.soil_bgc_competition!(st, nf, cf, ns, state, params;
            mask_bgc_soilc=mask, bounds=bounds,
            nlevdecomp=nlevdecomp, ndecomp_cascade_transitions=ntrans,
            dzsoi_decomp=dzsoi, pmnf_decomp_cascade=pmnf,
            p_decomp_cn_gain=pcng, cascade_receiver_pool=crp,
            use_nitrif_denitrif=false)

        @test st.fpg_col[1] ≈ 1.0
        @test st.fpi_col[1] ≈ 1.0
    end

    # -----------------------------------------------------------------------
    @testset "zero demand with nitrif_denitrif" begin
        nc = 1; nlevdecomp = 1; ntrans = 1
        st, nf, cf, ns, mask, bounds, dzsoi, pmnf, pcng, crp =
            make_competition_data(nc=nc, nlevdecomp=nlevdecomp, ntrans=ntrans,
                sminn_vr_val=10.0,
                smin_nh4_vr_val=5.0,
                smin_no3_vr_val=5.0,
                plant_ndemand_val=0.0,
                potential_immob_vr_val=0.0,
                pot_f_nit_vr_val=0.0,
                pot_f_denit_vr_val=0.0)

        params = CLM.SoilBGCCompetitionParams()
        state = CLM.SoilBGCCompetitionState(dt=1800.0, bdnr=0.01)

        CLM.soil_bgc_competition!(st, nf, cf, ns, state, params;
            mask_bgc_soilc=mask, bounds=bounds,
            nlevdecomp=nlevdecomp, ndecomp_cascade_transitions=ntrans,
            dzsoi_decomp=dzsoi, pmnf_decomp_cascade=pmnf,
            p_decomp_cn_gain=pcng, cascade_receiver_pool=crp,
            use_nitrif_denitrif=true)

        @test st.fpg_col[1] ≈ 1.0
        @test st.fpi_col[1] ≈ 1.0

        # All fluxes should be zero
        @test nf.actual_immob_vr_col[1, 1] ≈ 0.0
        @test nf.sminn_to_plant_vr_col[1, 1] ≈ 0.0
        @test nf.f_nit_vr_col[1, 1] ≈ 0.0
        @test nf.f_denit_vr_col[1, 1] ≈ 0.0
    end

end
