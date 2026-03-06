@testset "Maintenance Respiration" begin

    # =====================================================================
    # MaintRespParams struct
    # =====================================================================
    @testset "MaintRespParams defaults" begin
        p = CLM.MaintRespParams()
        @test p.br == CLM.SPVAL
        @test p.br_root == CLM.SPVAL
    end

    @testset "maint_resp_read_params! defaults br_root to br" begin
        p = CLM.MaintRespParams()
        CLM.maint_resp_read_params!(p; br=2.525e-6)
        @test p.br ≈ 2.525e-6
        @test p.br_root ≈ 2.525e-6
    end

    @testset "maint_resp_read_params! with explicit br_root" begin
        p = CLM.MaintRespParams()
        CLM.maint_resp_read_params!(p; br=2.525e-6, br_root=3.0e-6)
        @test p.br ≈ 2.525e-6
        @test p.br_root ≈ 3.0e-6
    end

    @testset "PftConMaintResp defaults" begin
        pfc = CLM.PftConMaintResp()
        @test isempty(pfc.woody)
    end

    # =====================================================================
    # Helper: create test data for cn_mresp!
    # =====================================================================
    function make_mresp_data(;
            np=1, nc=1, nlevgrnd=15, nlevsno=5,
            br=2.525e-6, br_root=2.525e-6, Q10=1.5,
            woody_val=1.0, ivt_val=1,
            frac_veg_nosno_val=1,
            t_ref2m_val=293.15,   # 20°C
            t10_val=293.15,
            t_soisno_val=293.15,  # uniform soil temperature
            lmrsun_val=1.0,
            lmrsha_val=0.5,
            laisun_val=2.0,
            laisha_val=1.0,
            frootn_val=0.01,
            livestemn_val=0.05,
            livecrootn_val=0.03,
            reproductiven_val=0.0,
            crootfr_uniform=true,
            rootstem_acc=false,
            npcropmin=17,
            nrepr=CLM.NREPR)

        nlevmaxurbgrnd = max(CLM.varpar.nlevurb, nlevgrnd)
        nlev_soisno = nlevsno + nlevmaxurbgrnd

        # Parameters
        params = CLM.MaintRespParams(br=br, br_root=br_root)
        cn_params = CLM.CNSharedParamsData(Q10=Q10)

        # PFT constants (index ivt+1 for 1-based Julia arrays)
        pftcon = CLM.PftConMaintResp(woody=[0.0, woody_val])

        # Patch data
        patch = CLM.PatchData()
        CLM.patch_init!(patch, np)
        patch.itype[1] = ivt_val
        patch.column[1] = 1

        # Canopy state
        cs = CLM.CanopyStateData()
        CLM.canopystate_init!(cs, np)
        cs.frac_veg_nosno_patch[1] = frac_veg_nosno_val
        cs.laisun_patch[1] = laisun_val
        cs.laisha_patch[1] = laisha_val

        # Soil state
        ss = CLM.SoilStateData()
        CLM.soilstate_init!(ss, np, nc)
        if crootfr_uniform
            for j in 1:nlevgrnd
                ss.crootfr_patch[1, j] = 1.0 / nlevgrnd
            end
        end

        # Temperature
        temp = CLM.TemperatureData()
        CLM.temperature_init!(temp, np, nc, 1, 1)
        temp.t_ref2m_patch[1] = t_ref2m_val
        temp.t_a10_patch[1] = t10_val
        for j in 1:nlevgrnd
            temp.t_soisno_col[1, j + nlevsno] = t_soisno_val
        end

        # Photosynthesis data
        ps = CLM.PhotosynthesisData()
        CLM.photosynthesis_data_init!(ps, np)
        ps.lmrsun_patch[1] = lmrsun_val
        ps.lmrsha_patch[1] = lmrsha_val
        ps.rootstem_acc = rootstem_acc

        # Carbon flux
        cnveg_cf = CLM.CNVegCarbonFluxData()
        CLM.cnveg_carbon_flux_init!(cnveg_cf, np, nc, 1)

        # Nitrogen state
        cnveg_ns = CLM.CNVegNitrogenStateData()
        CLM.cnveg_nitrogen_state_init!(cnveg_ns, np, nc, 1)
        cnveg_ns.frootn_patch[1] = frootn_val
        cnveg_ns.livestemn_patch[1] = livestemn_val
        cnveg_ns.livecrootn_patch[1] = livecrootn_val
        for k in 1:nrepr
            cnveg_ns.reproductiven_patch[1, k] = reproductiven_val
        end

        mask_soilc = BitVector([true])
        mask_soilp = BitVector([true])
        bounds_c = 1:nc
        bounds_p = 1:np

        return params, cn_params, pftcon, patch, cs, ss, temp, ps,
               cnveg_cf, cnveg_ns, mask_soilc, mask_soilp, bounds_c, bounds_p,
               nlevgrnd, nlevsno, npcropmin, nrepr
    end

    # =====================================================================
    # Test: woody PFT at reference temperature (20°C)
    # =====================================================================
    @testset "Woody PFT at reference temperature" begin
        br = 2.525e-6
        Q10 = 1.5
        t_ref2m = 293.15  # 20°C → tc = Q10^0 = 1.0
        t_soil  = 293.15  # 20°C → tcsoi = Q10^0 = 1.0

        params, cn_params, pftcon, patch, cs, ss, temp, ps,
            cnveg_cf, cnveg_ns, mask_soilc, mask_soilp, bounds_c, bounds_p,
            nlevgrnd, nlevsno, npcropmin, nrepr = make_mresp_data(
                br=br, Q10=Q10,
                t_ref2m_val=t_ref2m, t_soisno_val=t_soil,
                woody_val=1.0, ivt_val=1,
                lmrsun_val=1.0, lmrsha_val=0.5,
                laisun_val=2.0, laisha_val=1.0,
                frootn_val=0.01, livestemn_val=0.05, livecrootn_val=0.03)

        CLM.cn_mresp!(mask_soilc, mask_soilp, bounds_c, bounds_p,
                       params, cn_params, pftcon, patch, cs, ss, temp, ps,
                       cnveg_cf, cnveg_ns;
                       nlevgrnd=nlevgrnd, nlevsno=nlevsno,
                       npcropmin=npcropmin, nrepr=nrepr)

        # tc = Q10^((293.15 - 273.15 - 20) / 10) = Q10^0 = 1.0
        tc = 1.0

        # Leaf MR: lmrsun * laisun * 12.011e-6 + lmrsha * laisha * 12.011e-6
        expected_leaf_mr = 1.0 * 2.0 * 12.011e-6 + 0.5 * 1.0 * 12.011e-6
        @test cnveg_cf.leaf_mr_patch[1] ≈ expected_leaf_mr

        # Live stem MR: livestemn * br * tc
        @test cnveg_cf.livestem_mr_patch[1] ≈ 0.05 * br * tc

        # Live coarse root MR: livecrootn * br_root * tc
        @test cnveg_cf.livecroot_mr_patch[1] ≈ 0.03 * br * tc

        # Fine root MR: sum over layers of frootn * br_root * tcsoi * crootfr
        # At 20°C, tcsoi = 1.0 everywhere, sum(crootfr) = 1.0
        @test cnveg_cf.froot_mr_patch[1] ≈ 0.01 * br * 1.0 atol=1e-20
    end

    # =====================================================================
    # Test: non-woody, non-crop PFT (grass)
    # =====================================================================
    @testset "Non-woody grass PFT" begin
        br = 2.525e-6
        Q10 = 1.5

        params, cn_params, pftcon, patch, cs, ss, temp, ps,
            cnveg_cf, cnveg_ns, mask_soilc, mask_soilp, bounds_c, bounds_p,
            nlevgrnd, nlevsno, npcropmin, nrepr = make_mresp_data(
                br=br, Q10=Q10,
                t_ref2m_val=293.15, t_soisno_val=293.15,
                woody_val=0.0, ivt_val=1,
                frootn_val=0.02, livestemn_val=0.0, livecrootn_val=0.0)

        CLM.cn_mresp!(mask_soilc, mask_soilp, bounds_c, bounds_p,
                       params, cn_params, pftcon, patch, cs, ss, temp, ps,
                       cnveg_cf, cnveg_ns;
                       nlevgrnd=nlevgrnd, nlevsno=nlevsno,
                       npcropmin=npcropmin, nrepr=nrepr)

        # Leaf MR computed as usual
        expected_leaf_mr = 1.0 * 2.0 * 12.011e-6 + 0.5 * 1.0 * 12.011e-6
        @test cnveg_cf.leaf_mr_patch[1] ≈ expected_leaf_mr

        # livestem_mr and livecroot_mr should NOT be set (remain NaN)
        @test isnan(cnveg_cf.livestem_mr_patch[1])
        @test isnan(cnveg_cf.livecroot_mr_patch[1])

        # Fine root MR should be computed
        @test cnveg_cf.froot_mr_patch[1] ≈ 0.02 * br * 1.0 atol=1e-20
    end

    # =====================================================================
    # Test: no vegetation (frac_veg_nosno = 0)
    # =====================================================================
    @testset "No vegetation cover" begin
        params, cn_params, pftcon, patch, cs, ss, temp, ps,
            cnveg_cf, cnveg_ns, mask_soilc, mask_soilp, bounds_c, bounds_p,
            nlevgrnd, nlevsno, npcropmin, nrepr = make_mresp_data(
                frac_veg_nosno_val=0)

        CLM.cn_mresp!(mask_soilc, mask_soilp, bounds_c, bounds_p,
                       params, cn_params, pftcon, patch, cs, ss, temp, ps,
                       cnveg_cf, cnveg_ns;
                       nlevgrnd=nlevgrnd, nlevsno=nlevsno,
                       npcropmin=npcropmin, nrepr=nrepr)

        @test cnveg_cf.leaf_mr_patch[1] == 0.0
    end

    # =====================================================================
    # Test: temperature dependence (Q10 effect)
    # =====================================================================
    @testset "Temperature dependence" begin
        br = 2.525e-6
        Q10 = 1.5
        t_ref2m = 303.15  # 30°C

        params, cn_params, pftcon, patch, cs, ss, temp, ps,
            cnveg_cf, cnveg_ns, mask_soilc, mask_soilp, bounds_c, bounds_p,
            nlevgrnd, nlevsno, npcropmin, nrepr = make_mresp_data(
                br=br, Q10=Q10,
                t_ref2m_val=t_ref2m, t_soisno_val=t_ref2m,
                woody_val=1.0, ivt_val=1,
                livestemn_val=0.05, livecrootn_val=0.03, frootn_val=0.01)

        CLM.cn_mresp!(mask_soilc, mask_soilp, bounds_c, bounds_p,
                       params, cn_params, pftcon, patch, cs, ss, temp, ps,
                       cnveg_cf, cnveg_ns;
                       nlevgrnd=nlevgrnd, nlevsno=nlevsno,
                       npcropmin=npcropmin, nrepr=nrepr)

        # tc = Q10^((303.15 - 273.15 - 20) / 10) = Q10^1.0 = 1.5
        tc = Q10^1.0
        @test cnveg_cf.livestem_mr_patch[1] ≈ 0.05 * br * tc
        @test cnveg_cf.livecroot_mr_patch[1] ≈ 0.03 * br * tc
        @test cnveg_cf.froot_mr_patch[1] ≈ 0.01 * br * tc atol=1e-20
    end

    # =====================================================================
    # Test: crop PFT (ivt >= npcropmin)
    # =====================================================================
    @testset "Crop PFT maintenance respiration" begin
        br = 2.525e-6
        Q10 = 1.5
        npcropmin = 17

        # Need PFT arrays large enough for index ivt+1=18
        pftcon = CLM.PftConMaintResp(woody=zeros(npcropmin + 1))

        params = CLM.MaintRespParams(br=br, br_root=br)
        cn_params = CLM.CNSharedParamsData(Q10=Q10)

        np = 1; nc = 1
        nlevgrnd = CLM.varpar.nlevgrnd
        nlevsno = CLM.varpar.nlevsno
        nrepr = CLM.NREPR

        patch = CLM.PatchData()
        CLM.patch_init!(patch, np)
        patch.itype[1] = npcropmin
        patch.column[1] = 1

        cs = CLM.CanopyStateData()
        CLM.canopystate_init!(cs, np)
        cs.frac_veg_nosno_patch[1] = 1
        cs.laisun_patch[1] = 2.0
        cs.laisha_patch[1] = 1.0

        ss = CLM.SoilStateData()
        CLM.soilstate_init!(ss, np, nc)
        for j in 1:nlevgrnd
            ss.crootfr_patch[1, j] = 1.0 / nlevgrnd
        end

        temp = CLM.TemperatureData()
        CLM.temperature_init!(temp, np, nc, 1, 1)
        temp.t_ref2m_patch[1] = 293.15
        temp.t_a10_patch[1] = 293.15
        for j in 1:nlevgrnd
            temp.t_soisno_col[1, j + nlevsno] = 293.15
        end

        ps = CLM.PhotosynthesisData()
        CLM.photosynthesis_data_init!(ps, np)
        ps.lmrsun_patch[1] = 1.0
        ps.lmrsha_patch[1] = 0.5
        ps.rootstem_acc = false

        cnveg_cf = CLM.CNVegCarbonFluxData()
        CLM.cnveg_carbon_flux_init!(cnveg_cf, np, nc, 1)

        cnveg_ns = CLM.CNVegNitrogenStateData()
        CLM.cnveg_nitrogen_state_init!(cnveg_ns, np, nc, 1)
        cnveg_ns.frootn_patch[1] = 0.01
        cnveg_ns.livestemn_patch[1] = 0.04
        cnveg_ns.livecrootn_patch[1] = 0.0
        for k in 1:nrepr
            cnveg_ns.reproductiven_patch[1, k] = 0.02
        end

        mask_soilc = BitVector([true])
        mask_soilp = BitVector([true])

        CLM.cn_mresp!(mask_soilc, mask_soilp, 1:nc, 1:np,
                       params, cn_params, pftcon, patch, cs, ss, temp, ps,
                       cnveg_cf, cnveg_ns;
                       nlevgrnd=nlevgrnd, nlevsno=nlevsno,
                       npcropmin=npcropmin, nrepr=nrepr)

        tc = 1.0  # at 20°C

        # Crop livestem MR
        @test cnveg_cf.livestem_mr_patch[1] ≈ 0.04 * br * tc

        # Crop reproductive MR
        for k in 1:nrepr
            @test cnveg_cf.reproductive_mr_patch[1, k] ≈ 0.02 * br * tc
        end

        # livecroot_mr should NOT be set for crops (remain NaN)
        @test isnan(cnveg_cf.livecroot_mr_patch[1])
    end

    # =====================================================================
    # Test: masked-out patches are skipped
    # =====================================================================
    @testset "Masked-out patches skipped" begin
        br = 2.525e-6
        Q10 = 1.5

        np = 2; nc = 2
        nlevgrnd = CLM.varpar.nlevgrnd
        nlevsno = CLM.varpar.nlevsno
        nrepr = CLM.NREPR
        npcropmin = 17

        params = CLM.MaintRespParams(br=br, br_root=br)
        cn_params = CLM.CNSharedParamsData(Q10=Q10)
        pftcon = CLM.PftConMaintResp(woody=[0.0, 1.0])

        patch = CLM.PatchData()
        CLM.patch_init!(patch, np)
        patch.itype[1] = 1
        patch.itype[2] = 1
        patch.column[1] = 1
        patch.column[2] = 2

        cs = CLM.CanopyStateData()
        CLM.canopystate_init!(cs, np)
        cs.frac_veg_nosno_patch .= 1
        cs.laisun_patch .= 2.0
        cs.laisha_patch .= 1.0

        ss = CLM.SoilStateData()
        CLM.soilstate_init!(ss, np, nc)
        for p in 1:np
            for j in 1:nlevgrnd
                ss.crootfr_patch[p, j] = 1.0 / nlevgrnd
            end
        end

        temp = CLM.TemperatureData()
        CLM.temperature_init!(temp, np, nc, 1, 1)
        temp.t_ref2m_patch .= 293.15
        temp.t_a10_patch .= 293.15
        for c in 1:nc
            for j in 1:nlevgrnd
                temp.t_soisno_col[c, j + nlevsno] = 293.15
            end
        end

        ps = CLM.PhotosynthesisData()
        CLM.photosynthesis_data_init!(ps, np)
        ps.lmrsun_patch .= 1.0
        ps.lmrsha_patch .= 0.5
        ps.rootstem_acc = false

        cnveg_cf = CLM.CNVegCarbonFluxData()
        CLM.cnveg_carbon_flux_init!(cnveg_cf, np, nc, 1)

        cnveg_ns = CLM.CNVegNitrogenStateData()
        CLM.cnveg_nitrogen_state_init!(cnveg_ns, np, nc, 1)
        cnveg_ns.frootn_patch .= 0.01
        cnveg_ns.livestemn_patch .= 0.05
        cnveg_ns.livecrootn_patch .= 0.03

        # Only patch 1 / column 1 active
        mask_soilc = BitVector([true, false])
        mask_soilp = BitVector([true, false])

        CLM.cn_mresp!(mask_soilc, mask_soilp, 1:nc, 1:np,
                       params, cn_params, pftcon, patch, cs, ss, temp, ps,
                       cnveg_cf, cnveg_ns;
                       nlevgrnd=nlevgrnd, nlevsno=nlevsno,
                       npcropmin=npcropmin, nrepr=nrepr)

        # Patch 1 computed
        @test cnveg_cf.leaf_mr_patch[1] > 0.0
        @test cnveg_cf.froot_mr_patch[1] > 0.0

        # Patch 2 should remain NaN (not touched)
        @test isnan(cnveg_cf.leaf_mr_patch[2])
    end

    # =====================================================================
    # Test: root/stem acclimation
    # =====================================================================
    @testset "Root/stem acclimation" begin
        br = 2.525e-6
        Q10 = 1.5
        t10 = 283.15  # 10°C → acclimation factor = 10^(-0.00794*(10-25))

        params, cn_params, pftcon, patch, cs, ss, temp, ps,
            cnveg_cf, cnveg_ns, mask_soilc, mask_soilp, bounds_c, bounds_p,
            nlevgrnd, nlevsno, npcropmin, nrepr = make_mresp_data(
                br=br, Q10=Q10,
                t_ref2m_val=293.15, t_soisno_val=293.15,
                t10_val=t10,
                woody_val=1.0, ivt_val=1,
                livestemn_val=0.05, livecrootn_val=0.03, frootn_val=0.01,
                rootstem_acc=true)

        CLM.cn_mresp!(mask_soilc, mask_soilp, bounds_c, bounds_p,
                       params, cn_params, pftcon, patch, cs, ss, temp, ps,
                       cnveg_cf, cnveg_ns;
                       nlevgrnd=nlevgrnd, nlevsno=nlevsno,
                       npcropmin=npcropmin, nrepr=nrepr)

        tc = Q10^((293.15 - CLM.TFRZ - 20.0) / 10.0)  # = 1.0
        acc_factor = 10.0^(-0.00794 * ((t10 - CLM.TFRZ) - 25.0))
        br_acc = br * acc_factor

        @test cnveg_cf.livestem_mr_patch[1] ≈ 0.05 * br_acc * tc
        @test cnveg_cf.livecroot_mr_patch[1] ≈ 0.03 * br_acc * tc
        @test cnveg_cf.froot_mr_patch[1] ≈ 0.01 * br_acc * 1.0 atol=1e-20
    end

    # =====================================================================
    # Test: separate br_root rate
    # =====================================================================
    @testset "Separate br_root rate" begin
        br = 2.525e-6
        br_root = 5.0e-6  # different root rate
        Q10 = 1.5

        params, cn_params, pftcon, patch, cs, ss, temp, ps,
            cnveg_cf, cnveg_ns, mask_soilc, mask_soilp, bounds_c, bounds_p,
            nlevgrnd, nlevsno, npcropmin, nrepr = make_mresp_data(
                br=br, br_root=br_root, Q10=Q10,
                t_ref2m_val=293.15, t_soisno_val=293.15,
                woody_val=1.0, ivt_val=1,
                livestemn_val=0.05, livecrootn_val=0.03, frootn_val=0.01)

        CLM.cn_mresp!(mask_soilc, mask_soilp, bounds_c, bounds_p,
                       params, cn_params, pftcon, patch, cs, ss, temp, ps,
                       cnveg_cf, cnveg_ns;
                       nlevgrnd=nlevgrnd, nlevsno=nlevsno,
                       npcropmin=npcropmin, nrepr=nrepr)

        tc = 1.0

        # livestem uses br (not br_root)
        @test cnveg_cf.livestem_mr_patch[1] ≈ 0.05 * br * tc
        # livecroot uses br_root
        @test cnveg_cf.livecroot_mr_patch[1] ≈ 0.03 * br_root * tc
        # froot uses br_root
        @test cnveg_cf.froot_mr_patch[1] ≈ 0.01 * br_root * 1.0 atol=1e-20
    end

end
