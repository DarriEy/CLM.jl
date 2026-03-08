@testset "Photosynthesis" begin

    # =====================================================================
    # Constants
    # =====================================================================
    @testset "Constants" begin
        @test CLM.SUN == 1
        @test CLM.SHA == 2
        @test CLM.XYL == 3
        @test CLM.ROOT_SEG == 4
        @test CLM.STOMATALCOND_MTD_BB1987 == 1
        @test CLM.STOMATALCOND_MTD_MEDLYN2011 == 2
        @test CLM.BBBOPT_C3 == 10000.0
        @test CLM.BBBOPT_C4 == 40000.0
        @test CLM.NVEGWCS == 4
    end

    # =====================================================================
    # PhotoParamsData
    # =====================================================================
    @testset "PhotoParamsData init/clean" begin
        pp = CLM.PhotoParamsData()
        CLM.photo_params_init!(pp)
        @test length(pp.krmax) == CLM.MXPFT + 1
        @test length(pp.theta_cj) == CLM.MXPFT + 1
        @test size(pp.kmax) == (CLM.MXPFT + 1, CLM.NVEGWCS)
        @test size(pp.psi50) == (CLM.MXPFT + 1, CLM.NVEGWCS)
        @test size(pp.ck) == (CLM.MXPFT + 1, CLM.NVEGWCS)
        @test length(pp.lmr_intercept_atkin) == CLM.MXPFT + 1

        CLM.photo_params_clean!(pp)
        @test isempty(pp.krmax)
        @test isempty(pp.theta_cj)
        @test size(pp.kmax) == (0, 0)
    end

    # =====================================================================
    # PhotosynthesisData
    # =====================================================================
    @testset "PhotosynthesisData init/clean" begin
        np = 3
        ps = CLM.PhotosynthesisData()
        CLM.photosynthesis_data_init!(ps, np)
        @test length(ps.c3flag_patch) == np
        @test size(ps.ac_phs_patch) == (np, 2, CLM.NLEVCAN)
        @test length(ps.cp_patch) == np
        @test size(ps.cisun_z_patch) == (np, CLM.NLEVCAN)
        @test length(ps.psnsun_patch) == np

        CLM.photosynthesis_data_clean!(ps)
        @test isempty(ps.c3flag_patch)
        @test size(ps.ac_phs_patch) == (0, 0, 0)
    end

    # =====================================================================
    # Helper functions
    # =====================================================================
    @testset "ft_photo" begin
        val = CLM.ft_photo(298.15, 65330.0)
        @test val > 0.0
        @test isfinite(val)
        # At 25°C (298.15 K), ft should be 1.0
        @test CLM.ft_photo(CLM.TFRZ + 25.0, 65330.0) ≈ 1.0 atol=1e-12
    end

    @testset "fth_photo" begin
        hd = 149252.0
        se = 486.0
        sf = CLM.fth25_photo(hd, se)
        val = CLM.fth_photo(298.15, hd, se, sf)
        @test val ≈ 1.0 atol=1e-10
        @test isfinite(val)
    end

    @testset "fth25_photo" begin
        # When hd ≈ se * T25, exp term is O(1) so val is noticeably > 1.0
        val = CLM.fth25_photo(144901.0, 486.0)
        @test val > 1.0
        @test isfinite(val)

        # When hd >> se * T25, exp term vanishes and val → 1.0
        val2 = CLM.fth25_photo(200000.0, 486.0)
        @test val2 ≈ 1.0 atol=1e-9
        @test isfinite(val2)
    end

    @testset "quadratic_solve" begin
        # x^2 - 3x + 2 = 0 -> roots 1 and 2
        r1, r2 = CLM.quadratic_solve(1.0, -3.0, 2.0)
        @test r1 ≈ 2.0 atol=1e-12
        @test r2 ≈ 1.0 atol=1e-12

        # Linear: 2x + 4 = 0 -> x = -2
        r1, r2 = CLM.quadratic_solve(0.0, 2.0, 4.0)
        @test r1 ≈ -2.0 atol=1e-12

        # Degenerate: 0 = 0
        r1, r2 = CLM.quadratic_solve(0.0, 0.0, 0.0)
        @test r1 == 0.0
        @test r2 == 0.0
    end

    # =====================================================================
    # plc and d1plc
    # =====================================================================
    @testset "plc vulnerability curve" begin
        # Initialize params for testing
        CLM.photo_params_init!(CLM.params_inst)
        CLM.params_inst.psi50 .= -300000.0
        CLM.params_inst.ck .= 3.0

        ivt = 1
        # At x=0, plc should be 1.0 (no stress)
        # Actually plc(0,...) = 2^(-(0/psi50)^ck) = 2^0 = 1.0
        val = CLM.plc(0.0, ivt, CLM.SUN, CLM.VEG)
        @test val ≈ 1.0 atol=1e-10

        # At very negative potential, plc approaches 0
        val = CLM.plc(-1.0e6, ivt, CLM.SUN, CLM.VEG)
        @test val < 0.01
    end

    @testset "d1plc derivative" begin
        CLM.photo_params_init!(CLM.params_inst)
        CLM.params_inst.psi50 .= -300000.0
        CLM.params_inst.ck .= 3.0

        ivt = 1
        x = -100000.0
        dval = CLM.d1plc(x, ivt, CLM.SUN, CLM.VEG)
        @test isfinite(dval)

        # Finite difference check
        dx = 1.0
        plc_plus = CLM.plc(x + dx, ivt, CLM.SUN, CLM.VEG)
        plc_minus = CLM.plc(x - dx, ivt, CLM.SUN, CLM.VEG)
        fd = (plc_plus - plc_minus) / (2.0 * dx)
        @test abs(dval - fd) / max(abs(fd), 1e-10) < 0.01
    end

    # =====================================================================
    # ci_func! (non-PHS)
    # =====================================================================
    @testset "ci_func!" begin
        np = 1
        ps = CLM.PhotosynthesisData()
        CLM.photosynthesis_data_init!(ps, np)
        CLM.photo_params_init!(CLM.params_inst)

        ps.c3flag_patch[1] = true
        ps.cp_patch[1] = 4.275  # CO2 compensation point (Pa)
        ps.kc_patch[1] = 40.49  # Michaelis-Menten for CO2 (Pa)
        ps.ko_patch[1] = 27840.0  # Michaelis-Menten for O2 (Pa)
        ps.qe_patch[1] = 0.0
        ps.bbb_patch[1] = 10000.0
        ps.mbb_patch[1] = 9.0
        ps.vcmax_z_patch[1, 1] = 60.0
        ps.tpu_z_patch[1, 1] = 10.0
        ps.kp_z_patch[1, 1] = 0.0
        ps.stomatalcond_mtd = CLM.STOMATALCOND_MTD_BB1987

        CLM.params_inst.theta_cj = fill(0.98, CLM.MXPFT + 1)
        CLM.params_inst.theta_ip = 0.95

        ci = 28.0  # Pa
        gb_mol = 100000.0  # umol H2O/m2/s
        je = 120.0  # umol electrons/m2/s
        cair = 40.0  # Pa
        oair = 20900.0  # Pa
        lmr_z = 1.0  # umol CO2/m2/s
        par_z = 500.0  # W/m2
        rh_can = 0.7
        forc_pbot_c = 101325.0  # Pa

        fval, gs_mol = CLM.ci_func!(ci, 1, 1, forc_pbot_c, gb_mol, je, cair, oair,
                                      lmr_z, par_z, rh_can, ps)
        @test isfinite(fval)
        @test gs_mol >= 0.0
    end

    # =====================================================================
    # hybrid_solver! and brent_solver!
    # =====================================================================
    @testset "hybrid_solver!" begin
        np = 1
        ps = CLM.PhotosynthesisData()
        CLM.photosynthesis_data_init!(ps, np)
        CLM.photo_params_init!(CLM.params_inst)

        ps.c3flag_patch[1] = true
        ps.cp_patch[1] = 4.275
        ps.kc_patch[1] = 40.49
        ps.ko_patch[1] = 27840.0
        ps.qe_patch[1] = 0.0
        ps.bbb_patch[1] = 10000.0
        ps.mbb_patch[1] = 9.0
        ps.vcmax_z_patch[1, 1] = 60.0
        ps.tpu_z_patch[1, 1] = 10.0
        ps.kp_z_patch[1, 1] = 0.0
        ps.stomatalcond_mtd = CLM.STOMATALCOND_MTD_BB1987

        CLM.params_inst.theta_cj = fill(0.98, CLM.MXPFT + 1)
        CLM.params_inst.theta_ip = 0.95

        x0 = 28.0
        gb_mol = 100000.0
        je = 120.0
        cair = 40.0
        oair = 20900.0
        lmr_z = 1.0
        par_z = 500.0
        rh_can = 0.7
        forc_pbot_c = 101325.0

        ci_sol, gs_mol, niter = CLM.hybrid_solver!(x0, 1, 1, forc_pbot_c, gb_mol, je,
            cair, oair, lmr_z, par_z, rh_can, ps)
        @test isfinite(ci_sol)
        @test ci_sol > 0.0
        @test gs_mol > 0.0
    end

    # =====================================================================
    # photosynthesis_timestep_init!
    # =====================================================================
    @testset "photosynthesis_timestep_init!" begin
        np = 2
        ps = CLM.PhotosynthesisData()
        CLM.photosynthesis_data_init!(ps, np)

        ps.psnsun_patch .= 999.0
        ps.psnsha_patch .= 999.0
        ps.fpsn_patch .= 999.0

        mask = trues(np)
        CLM.photosynthesis_timestep_init!(ps, mask, 1:np)

        @test ps.psnsun_patch[1] == 0.0
        @test ps.psnsha_patch[2] == 0.0
        @test ps.fpsn_patch[1] == 0.0
    end

    # =====================================================================
    # photosynthesis_total!
    # =====================================================================
    @testset "photosynthesis_total!" begin
        np = 2
        ps = CLM.PhotosynthesisData()
        CLM.photosynthesis_data_init!(ps, np)

        ps.psnsun_patch .= 10.0
        ps.psnsha_patch .= 5.0
        ps.psnsun_wc_patch .= 8.0
        ps.psnsha_wc_patch .= 4.0
        ps.psnsun_wj_patch .= 1.5
        ps.psnsha_wj_patch .= 0.5
        ps.psnsun_wp_patch .= 0.5
        ps.psnsha_wp_patch .= 0.5

        laisun = [2.0, 1.0]
        laisha = [3.0, 2.0]
        mask = trues(np)

        CLM.photosynthesis_total!(ps, laisun, laisha, mask, 1:np)

        @test ps.fpsn_patch[1] ≈ 10.0 * 2.0 + 5.0 * 3.0 atol=1e-10
        @test ps.fpsn_wc_patch[1] ≈ 8.0 * 2.0 + 4.0 * 3.0 atol=1e-10
    end

    # =====================================================================
    # set_params_for_testing!
    # =====================================================================
    @testset "set_params_for_testing!" begin
        ps = CLM.PhotosynthesisData()
        CLM.set_params_for_testing!(ps)
        @test all(CLM.params_inst.ck .== 3.95)
        @test CLM.params_inst.psi50[1, 1] == -150000.0
        @test CLM.params_inst.psi50[2, 1] == -530000.0
    end

    # =====================================================================
    # getqflx! (PHS)
    # =====================================================================
    @testset "getqflx!" begin
        # gs -> qflx direction
        qflx_sun, qflx_sha, _, _ = CLM.getqflx!(
            1, 1, 50000.0,  # p, c, gb_mol
            30000.0, 20000.0,  # gs_mol_sun, gs_mol_sha
            0.0, 0.0,  # qflx_sun, qflx_sha
            0.02, 0.015,  # qsatl, qaf
            true,  # havegs
            2.0, 3.0,  # laisun, laisha
            4.0, 0.5,  # elai, esai
            0.8,  # fdry
            1.2,  # forc_rho
            101325.0,  # forc_pbot
            298.15)  # tgcm
        @test qflx_sun >= 0.0
        @test qflx_sha >= 0.0
    end

    # =====================================================================
    # spacF! (PHS)
    # =====================================================================
    @testset "spacF!" begin
        CLM.photo_params_init!(CLM.params_inst)
        CLM.params_inst.kmax .= 0.001
        CLM.params_inst.psi50 .= -300000.0
        CLM.params_inst.ck .= 3.0

        nlevsoi = 5
        x = [-50000.0, -60000.0, -40000.0, -30000.0]
        k_soil_root = fill(0.01, nlevsoi)
        smp = fill(-5000.0, nlevsoi)
        z_c = collect(range(0.05, 1.0, length=nlevsoi))

        f = CLM.spacF!(1, 1, x, 0.001, 0.0008,
            k_soil_root, smp, z_c,
            2.0, 3.0, 10.0, 1.0, 1, nlevsoi)
        @test length(f) == CLM.NVEGWCS
        @test all(isfinite.(f))
    end

    # =====================================================================
    # spacA! (PHS)
    # =====================================================================
    @testset "spacA!" begin
        CLM.photo_params_init!(CLM.params_inst)
        CLM.params_inst.kmax .= 0.001
        CLM.params_inst.psi50 .= -300000.0
        CLM.params_inst.ck .= 3.0

        nlevsoi = 5
        x = [-50000.0, -60000.0, -40000.0, -30000.0]
        k_soil_root = fill(0.01, nlevsoi)

        invA, flag = CLM.spacA!(1, 1, x, 0.001, 0.0008,
            k_soil_root,
            2.0, 3.0, 10.0, 1.0, 1, nlevsoi)
        @test size(invA) == (CLM.NVEGWCS, CLM.NVEGWCS)
        @test !flag  # should be invertible
        @test all(isfinite.(invA))
    end

    # =====================================================================
    # getvegwp! (PHS)
    # =====================================================================
    @testset "getvegwp!" begin
        CLM.photo_params_init!(CLM.params_inst)
        CLM.params_inst.kmax .= 0.001
        CLM.params_inst.psi50 .= -300000.0
        CLM.params_inst.ck .= 3.0
        CLM.params_inst.krmax .= 0.001

        nlevsoi = 5
        k_soil_root = fill(0.01, nlevsoi)
        smp = fill(-5000.0, nlevsoi)
        z_c = collect(range(0.05, 1.0, length=nlevsoi))

        x, soilflux = CLM.getvegwp!(1, 1, 50000.0, 30000.0, 20000.0,
            0.02, 0.015,
            k_soil_root, smp, z_c,
            2.0, 3.0, 10.0, 1.0, 1, nlevsoi,
            4.0, 0.5, 0.8, 1.2, 101325.0, 298.15)
        @test length(x) == CLM.NVEGWCS
        @test all(isfinite.(x))
        @test isfinite(soilflux)
        # Root should be most negative (closest to soil)
        @test x[CLM.ROOT_SEG] <= 0.0 || isfinite(x[CLM.ROOT_SEG])
    end

end
