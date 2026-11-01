# ==========================================================================
# test_photosynthesis_lnc_opt.jl — clm_varctl lnc_opt / vcmax_opt in Pass 2.
#
# PhotosynthesisMod.F90:1547-1630 (Photosynthesis) and 3160-3239
# (PhotosynthesisHydraulicStress):
#   lnc_opt = .false. : lnc = 1 / (slatop * leafcn)
#   lnc_opt = .true.  : lnc = leafn / (tlai * sum_nscaler),  sum_nscaler = nrad
#                       for nlevcan = 1 (nrad <= nlevcan); PHS then clips lnc at 10 (l.3206)
#   vcmax_opt = 0 : vcmax25top = lnc * flnr * fnr * act25 * dayl_factor
#   vcmax_opt = 3 : vcmax25top = (i_vcad + s_vcad * lnc) * dayl_factor
#   vcmax_opt = 4 : tree     lnc * i_flnr * exp(s_flnr * lnc) * fnr * act25 * dayl
#                   non-tree lnc * (i_flnr + s_flnr * lnc) * fnr * act25 * dayl
# With use_cn and leafresp_method = Atkin2015,
#   lmr25top = lmr_intercept_atkin + 0.2061 * lnc - 0.0402 * (t10 - tfrz).
# CTSM's namelist turns lnc_opt on whenever use_cn (and vcmax_opt=3 under
# use_flexibleCN); the port's defaults stay at the clm_varctl module defaults
# (.false., 0). The Bow reference runs lnc_opt=.true., vcmax_opt=3, where the
# static lnc was ~13% low (docs/CH4_FIRE_PARITY.md §18).
#
# Expected values are recomputed from the Fortran expressions, not copied from the
# implementation, and the inputs are chosen so every option gives a different lnc.
# ==========================================================================

using Test
using CLM

@testset "photosynthesis lnc_opt / vcmax_opt" begin
    FT = Float64
    nlevcan = CLM.NLEVCAN
    mp = CLM.MXPFT + 1
    IVT_TREE  = 2          # 1-based Julia index = Fortran itype 1 (< nptreemax = 9)
    IVT_GRASS = 13         # Fortran itype 12 (>= nptreemax)
    ivt = [IVT_TREE, IVT_GRASS]
    np = length(ivt)

    CLM.photo_params_init!(CLM.params_inst)
    pp = CLM.params_inst
    pp.theta_cj .= 0.98; pp.krmax .= 0.001
    pp.theta_ip = 0.95; pp.act25 = 72.0; pp.fnr = 7.16; pp.cp25_yr2000 = 42.75e-6
    pp.kc25_coef = 404.9e-6; pp.ko25_coef = 278.4e-3; pp.fnps = 0.15; pp.theta_psii = 0.7
    pp.vcmaxha = 65330.0; pp.jmaxha = 43540.0; pp.tpuha = 53100.0; pp.lmrha = 46390.0
    pp.kcha = 79430.0; pp.koha = 36380.0; pp.cpha = 37830.0
    pp.vcmaxhd = 149250.0; pp.jmaxhd = 152040.0; pp.tpuhd = 150650.0; pp.lmrhd = 150650.0
    pp.lmrse = 490.0; pp.tpu25ratio = 0.167; pp.kp25ratio = 20160.0
    pp.vcmaxse_sf = 1.0; pp.jmaxse_sf = 1.0; pp.tpuse_sf = 1.0; pp.jmax25top_sf = 1.0
    old_atkin = copy(pp.lmr_intercept_atkin)
    pp.lmr_intercept_atkin = fill(1.4, mp)

    # vcmax_opt coefficients live in pftcon (paramfile); set and restore.
    pc = CLM.pftcon
    saved = (copy(pc.i_vcad), copy(pc.s_vcad), copy(pc.i_flnr), copy(pc.s_flnr))
    pc.i_vcad = fill(15.0, mp); pc.s_vcad = fill(20.0, mp)
    pc.i_flnr = fill(0.05, mp); pc.s_flnr = fill(0.3, mp)

    T_VEG = 293.0; T10 = 290.0; SLATOP = 0.01; LEAFCN = 25.0; FLNR = 0.1; DAYL = 0.8
    LEAFN = [2.0, 0.9]; TLAI = [3.5, 1.5]; NRAD = 1   # nrad <= nlevcan = 1
    pftv(x) = fill(FT(x), mp)

    function run_pass2(phs::Bool, lnc_opt::Bool, vcmax_opt::Int)
        ps = CLM.PhotosynthesisData{FT}()
        CLM.photosynthesis_data_init!(ps, np; nlevcan = nlevcan)
        ps.stomatalcond_mtd = CLM.STOMATALCOND_MTD_BB1987
        ps.leafresp_method  = CLM.LEAFRESP_MTD_ATKIN2015
        ps.light_inhibit    = false
        ps.c3flag_patch    .= true
        kn = zeros(FT, np)
        tlai_z = fill(FT(1.0), np, nlevcan)
        if phs
            jmax = zeros(FT, np, 2, nlevcan)
            CLM.psn_phs_pass2_update!(ps, kn, jmax, fill(true, np), ivt,
                pftv(1.0), pftv(9.0), fill(FT(101325.0), np), fill(FT(20900.0), np),
                pftv(SLATOP), pftv(LEAFCN), pftv(FLNR), pftv(1.0), pftv(0.0),
                fill(FT(DAYL), np), fill(FT(T10), np), fill(FT(T_VEG), np), tlai_z,
                fill(FT(250.0), np, nlevcan), fill(FT(250.0), np, nlevcan),
                fill(FT(0.9), np), fill(FT(0.5), np), fill(NRAD, np),
                true, false, 0.015, nlevcan, ps.stomatalcond_mtd, ps.light_inhibit,
                ps.leafresp_method, CLM.CalibrationOverrides(), 1:np;
                lnc_opt = lnc_opt, vcmax_opt = vcmax_opt, leafn = LEAFN, tlai = TLAI)
        else
            jmax = zeros(FT, np, nlevcan)
            CLM.psn_pass2_update!(ps, kn, jmax, fill(true, np), ivt,
                pftv(SLATOP), pftv(LEAFCN), pftv(FLNR), pftv(1.0), pftv(0.0),
                fill(FT(DAYL), np), fill(FT(T10), np), fill(FT(T_VEG), np),
                fill(FT(1.0), np), tlai_z, fill(FT(250.0), np, nlevcan),
                fill(FT(0.9), np), fill(FT(0.9), np), fill(NRAD, np),
                true, false, 0.015, nlevcan, true, CLM.CalibrationOverrides(), 1:np;
                lnc_opt = lnc_opt, vcmax_opt = vcmax_opt, leafn = LEAFN, tlai = TLAI)
        end
        # Both paths store kn = exp(0.00963 * vcmax25top / dayl_factor - 2.43)
        # (F90:1648 / 3257); invert it to read vcmax25top on either path.
        return (; lnc = copy(ps.lnca_patch), vc25 = (log.(kn) .+ 2.43) .* DAYL ./ 0.00963,
                 lmr = copy(ps.lmrsun_z_patch[:, 1]))
    end

    lnc_static = 1.0 / (SLATOP * LEAFCN)
    lnc_leafn  = LEAFN ./ (TLAI .* NRAD)                  # nscaler = 1 per layer (nlevcan = 1)
    @test all(abs.(lnc_leafn .- lnc_static) ./ lnc_static .> 0.05)

    for phs in (false, true)
        @testset "$(phs ? "PHS" : "non-PHS") pass 2" begin
            off = run_pass2(phs, false, 0)
            on  = run_pass2(phs, true, 0)
            @test off.lnc ≈ fill(lnc_static, np) rtol = 1e-14
            @test on.lnc  ≈ lnc_leafn            rtol = 1e-14
            # vcmax_opt = 0
            @test on.vc25 ≈ lnc_leafn .* FLNR .* pp.fnr .* pp.act25 .* DAYL rtol = 1e-10
            # Atkin leaf MR follows lnc (same t_veg, same nscaler on both runs)
            atk(l) = 1.4 + 0.2061 * l - 0.0402 * (T10 - CLM.TFRZ)
            @test on.lmr ./ off.lmr ≈ atk.(lnc_leafn) ./ atk(lnc_static) rtol = 1e-12
            # vcmax_opt = 3
            v3 = run_pass2(phs, true, 3)
            @test v3.vc25 ≈ (15.0 .+ 20.0 .* lnc_leafn) .* DAYL rtol = 1e-10
            # vcmax_opt = 4: patch 1 is a tree, patch 2 is not
            v4 = run_pass2(phs, true, 4)
            l1, l2 = lnc_leafn
            @test v4.vc25[1] ≈ l1 * (0.05 * exp(0.3 * l1)) * pp.fnr * pp.act25 * DAYL rtol = 1e-10
            @test v4.vc25[2] ≈ l2 * (0.05 + 0.3 * l2) * pp.fnr * pp.act25 * DAYL rtol = 1e-10
        end
    end

    @testset "PHS clips lnc at 10 (PhotosynthesisMod.F90:3206)" begin
        # leafn / (tlai * nrad) = 15 / (1 * 1) = 15 > 10; non-PHS has no clip.
        LEAFN .= [15.0, 15.0]; TLAI .= [1.0, 1.0]
        @test run_pass2(true, true, 0).lnc ≈ fill(10.0, np) rtol = 1e-12
        @test run_pass2(false, true, 0).lnc ≈ fill(15.0, np) rtol = 1e-14
    end

    @testset "lnc_opt needs leafn/tlai; unsupported vcmax_opt errors" begin
        @test_throws ErrorException CLM._psn_lnc_bundle(zeros(2), true, 0, nothing, nothing)
        @test_throws ErrorException CLM._psn_lnc_bundle(zeros(2), false, 2, nothing, nothing)
    end

    pc.i_vcad, pc.s_vcad, pc.i_flnr, pc.s_flnr = saved
    pp.lmr_intercept_atkin = old_atkin
end
