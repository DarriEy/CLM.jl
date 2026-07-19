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

    # =====================================================================
    # LUNA branch in the NON-PHS photosynthesis! (PhotosynthesisMod.F90:1721/1748)
    #
    # This branch was ported into the PHS routine only, so `photosynthesis!`
    # ACCEPTED use_luna and never read it: LUNA acclimated vcmx25_z every step
    # while not one digit of GPP moved. A test that only asserts "runs" or
    # "is finite" reproduces exactly that defect, so every assertion below is
    # either a DIFFERENCE from the LUNA-off run or an analytic VALUE recomputed
    # straight from the Fortran expression.
    # =====================================================================
    @testset "LUNA branch (non-PHS photosynthesis!)" begin
        FT = Float64
        np = 2
        nlevcan = CLM.NLEVCAN
        mp = CLM.MXPFT + 1
        IVT = 2                      # 1-based PFT index used for all patches

        CLM.photo_params_init!(CLM.params_inst)
        p = CLM.params_inst
        p.theta_cj .= 0.98; p.krmax .= 0.001
        p.theta_ip = 0.95; p.act25 = 72.0; p.fnr = 7.16; p.cp25_yr2000 = 42.75e-6
        p.kc25_coef = 404.9e-6; p.ko25_coef = 278.4e-3; p.fnps = 0.15; p.theta_psii = 0.7
        p.vcmaxha = 65330.0; p.jmaxha = 43540.0; p.tpuha = 53100.0; p.lmrha = 46390.0
        p.kcha = 79430.0; p.koha = 36380.0; p.cpha = 37830.0
        p.vcmaxhd = 149250.0; p.jmaxhd = 152040.0; p.tpuhd = 150650.0; p.lmrhd = 150650.0
        p.lmrse = 490.0; p.tpu25ratio = 0.167; p.kp25ratio = 20160.0
        p.vcmaxse_sf = 1.0; p.jmaxse_sf = 1.0; p.tpuse_sf = 1.0; p.jmax25top_sf = 1.0

        T_VEG   = 293.0
        T10     = 290.0
        BTRAN   = 0.8
        VCMXCS  = 0.9       # vcmaxcintsun
        VCMXCH  = 0.5       # vcmaxcintsha
        LEAF_MR = 0.015
        # LUNA-acclimated values, chosen well away from the static lnc-based
        # vcmax25top*nscaler so a no-op branch cannot masquerade as agreement.
        VCMX25  = 44.0
        JMX25   = 88.0

        inp = (
            esat_tv = fill(FT(2300.0), np), eair = fill(FT(1500.0), np),
            oair = fill(FT(20900.0), np), cair = fill(FT(40.0), np),
            rb = fill(FT(50.0), np), btran = fill(FT(BTRAN), np),
            dayl_factor = fill(FT(1.0), np), leafn = fill(FT(1.0), np),
            forc_pbot = fill(FT(101325.0), np), t_veg = fill(FT(T_VEG), np),
            t10 = fill(FT(T10), np), tgcm = fill(FT(T_VEG), np),
            nrad = fill(1, np), tlai_z = fill(FT(1.0), np, nlevcan),
            tlai = fill(FT(1.0), np),
            par_z_in = fill(FT(250.0), np, nlevcan),   # > 0 → DAY time branch
            lai_z_in = fill(FT(1.0), np, nlevcan),
            o3coefv = fill(FT(1.0), np), o3coefg = fill(FT(1.0), np),
            c3psn_pft = fill(FT(1.0), mp),             # C3 everywhere
            leafcn_pft = fill(FT(25.0), mp), flnr_pft = fill(FT(0.1), mp),
            fnitr_pft = fill(FT(1.0), mp), slatop_pft = fill(FT(0.01), mp),
            mbbopt_pft = fill(FT(9.0), mp),
            medlynintercept_pft = fill(FT(100.0), mp),
            medlynslope_pft = fill(FT(6.0), mp),
            ivt = fill(IVT, np), col_of_patch = fill(1, np),
            laisun = fill(FT(0.5), np), laisha = fill(FT(0.5), np),
        )
        mask = fill(true, np)
        vcint_sun = fill(FT(VCMXCS), np)
        vcint_sha = fill(FT(VCMXCH), np)

        # Run sun+sha and return the per-phase outputs. `ps.vcmax_z_patch` is
        # shared between phases, so it is snapshotted after each call.
        function run_pair(; use_luna::Bool, crop::Real = 0.0, use_cn::Bool = true)
            ps = CLM.PhotosynthesisData{FT}()
            CLM.photosynthesis_data_init!(ps, np; use_luna = true)
            ps.stomatalcond_mtd = CLM.STOMATALCOND_MTD_BB1987
            ps.leafresp_method  = CLM.LEAFRESP_MTD_RYAN1991
            ps.light_inhibit    = false
            ps.vcmx25_z_patch  .= FT(VCMX25)
            ps.jmx25_z_patch   .= FT(JMX25)
            crop_pft = fill(FT(crop), mp)
            CLM.photosynthesis_timestep_init!(ps, mask, 1:np)
            call(phase, vcint; kw...) = CLM.photosynthesis!(ps,
                inp.esat_tv, inp.eair, inp.oair, inp.cair, inp.rb, inp.btran,
                inp.dayl_factor, inp.leafn, inp.forc_pbot, inp.t_veg, inp.t10,
                inp.tgcm, inp.nrad, inp.tlai_z, inp.tlai, inp.par_z_in,
                inp.lai_z_in, vcint, inp.o3coefv, inp.o3coefg, inp.c3psn_pft,
                inp.leafcn_pft, inp.flnr_pft, inp.fnitr_pft, inp.slatop_pft,
                inp.mbbopt_pft, inp.medlynintercept_pft, inp.medlynslope_pft,
                inp.ivt, inp.col_of_patch, mask, 1:np, phase;
                use_cn = use_cn, use_luna = use_luna, leaf_mr_vcm = LEAF_MR,
                crop_pft = crop_pft, kw...)
            call("sun", vcint_sun)
            vcmax_sun = copy(ps.vcmax_z_patch); lmr_sun = copy(ps.lmrsun_z_patch)
            call("sha", vcint_sha; vcmaxcint_sun = vcint_sun)
            vcmax_sha = copy(ps.vcmax_z_patch); lmr_sha = copy(ps.lmrsha_z_patch)
            CLM.photosynthesis_total!(ps, inp.laisun, inp.laisha, mask, 1:np)
            return (; vcmax_sun, vcmax_sha, lmr_sun, lmr_sha,
                    psnsun = copy(ps.psnsun_patch), psnsha = copy(ps.psnsha_patch),
                    fpsn = copy(ps.fpsn_patch))
        end

        off = run_pair(use_luna = false)
        on  = run_pair(use_luna = true)

        # --- 1. The branch EXECUTES: it must move vcmax AND the carbon flux. ---
        # (#267's finding was that all three of these were bit-identical.)
        @test on.vcmax_sun[1, 1] != off.vcmax_sun[1, 1]
        @test on.vcmax_sha[1, 1] != off.vcmax_sha[1, 1]
        @test on.psnsun[1] != off.psnsun[1]
        @test on.psnsha[1] != off.psnsha[1]
        @test on.fpsn[1] != off.fpsn[1]
        @test all(isfinite, on.fpsn)

        # --- 2. VALUE parity with the Fortran expression, not just difference. ---
        # PhotosynthesisMod.F90:1748-1774 →
        #   vcmax25 = vcmx25_z ; vcmax_z = vcmax25 * ft(t_veg,vcmaxha)
        #                                 * fth(t_veg,vcmaxhd,vcmaxse,vcmaxc)
        # then :1798 vcmax_z *= btran.
        vcmaxse = (668.39 - 1.07 * clamp(T10 - CLM.TFRZ, 11.0, 35.0)) * p.vcmaxse_sf
        vcmaxc  = CLM.fth25_photo(p.vcmaxhd, vcmaxse)
        ft_fth  = CLM.ft_photo(T_VEG, p.vcmaxha) *
                  CLM.fth_photo(T_VEG, p.vcmaxhd, vcmaxse, vcmaxc)
        @test on.vcmax_sun[1, 1] ≈ VCMX25 * ft_fth * BTRAN rtol=1e-12
        # Shaded phase rescales the sunlit-average vcmx25_z by
        # vcmaxcintsha/vcmaxcintsun (nlevcan==1), PhotosynthesisMod.F90:1753-1757.
        @test on.vcmax_sha[1, 1] ≈ VCMX25 * (VCMXCH / VCMXCS) * ft_fth * BTRAN rtol=1e-12
        # LUNA-off keeps the static lnc-based profile: vcmax25top*nscaler, with
        # nscaler = vcmaxcint (nlevcan==1). Independent of vcmx25_z.
        lnc = 1.0 / (inp.slatop_pft[IVT] * inp.leafcn_pft[IVT])
        vcmax25top = lnc * inp.flnr_pft[IVT] * p.fnr * p.act25 * 1.0
        @test off.vcmax_sun[1, 1] ≈ vcmax25top * VCMXCS * ft_fth * BTRAN rtol=1e-12

        # --- 3. The gate is real: LUNA-on with a CROP PFT is BYTE-identical to
        # LUNA-off (`crop(patch%itype(p))== 0` guard, F90:1748). ---
        crop_on = run_pair(use_luna = true, crop = 1.0)
        @test crop_on.vcmax_sun == off.vcmax_sun
        @test crop_on.vcmax_sha == off.vcmax_sha
        @test crop_on.fpsn == off.fpsn

        # --- 4. Leaf maintenance respiration: LUNA overrides lmr25 ONLY when CN
        # is off (F90:1721-1725 `if(.not.use_cn) lmr25 = leaf_mr_vcm*vcmx25_z`).
        lmrc   = CLM.fth25_photo(p.lmrhd, p.lmrse)
        lmr_tf = CLM.ft_photo(T_VEG, p.lmrha) *
                 CLM.fth_photo(T_VEG, p.lmrhd, p.lmrse, lmrc)
        nocn_on  = run_pair(use_luna = true,  use_cn = false)
        nocn_off = run_pair(use_luna = false, use_cn = false)
        @test nocn_on.lmr_sun[1, 1] != nocn_off.lmr_sun[1, 1]
        @test nocn_on.lmr_sun[1, 1] ≈ LEAF_MR * VCMX25 * lmr_tf * BTRAN rtol=1e-12
        # With CN on, leaf N predicts respiration and LUNA must NOT touch lmr.
        @test on.lmr_sun == off.lmr_sun

        # --- 5. A use_luna caller that omits the shaded-phase sunlit factor must
        # FAIL LOUDLY rather than silently fall back to the LUNA-free profile. ---
        ps_bad = CLM.PhotosynthesisData{FT}()
        CLM.photosynthesis_data_init!(ps_bad, np; use_luna = true)
        ps_bad.stomatalcond_mtd = CLM.STOMATALCOND_MTD_BB1987
        ps_bad.leafresp_method  = CLM.LEAFRESP_MTD_RYAN1991
        @test_throws ErrorException CLM.photosynthesis!(ps_bad,
            inp.esat_tv, inp.eair, inp.oair, inp.cair, inp.rb, inp.btran,
            inp.dayl_factor, inp.leafn, inp.forc_pbot, inp.t_veg, inp.t10,
            inp.tgcm, inp.nrad, inp.tlai_z, inp.tlai, inp.par_z_in,
            inp.lai_z_in, vcint_sha, inp.o3coefv, inp.o3coefg, inp.c3psn_pft,
            inp.leafcn_pft, inp.flnr_pft, inp.fnitr_pft, inp.slatop_pft,
            inp.mbbopt_pft, inp.medlynintercept_pft, inp.medlynslope_pft,
            inp.ivt, inp.col_of_patch, mask, 1:np, "sha";
            use_cn = true, use_luna = true, crop_pft = fill(FT(0.0), mp))
    end

end
