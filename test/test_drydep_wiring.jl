# ==========================================================================
# test_drydep_wiring.jl
#
# Live-path WIRING tests for dry deposition velocity + the annual-LAI producer.
#
# What is validated here (the WIRING added on top of the ported kernels):
#   1. read_annual_vegetation! populates canopystate.annlai_patch (finite, not
#      NaN) from the monthly-LAI stream — the init producer wired into
#      clm_driver.jl's is_first_step block (gated on n_drydep > 0), faithful to
#      CTSM readAnnualVegetation (clm_initializeMod.F90:698-699).
#   2. depvel_compute! assembles V_d = 100/(R_a+R_b+R_c) [cm/s] correctly, matched
#      to an INDEPENDENT plain-scalar re-implementation of the ported Wesely
#      formula (rtol 1e-10) for two representative species (O3-like, SO2-like),
#      with the Wesely index_season now selected from annlai (see below).
#   3. drydep_p2g! averages patch velocities to the gridcell field l2a.ddvel_grc,
#      the coupler field CAM's dry deposition consumes (CTSM lnd2atmMod.F90:307).
#   4. Default (n_drydep == 0) is byte-identical: the consumer + producer are
#      no-ops, and annlai stays NaN. Both driver insertions live strictly inside
#      `if config.n_drydep > 0`, so a default run adds no computation at all.
#
# NOTE on the port vs Fortran: the ported depvel kernel's surface-resistance rc is
# a simplified Wesely reimplementation (rc is NOT a line-by-line Fortran port).
# BUT the Wesely index_season that selects which seasonal COLUMN of the resistance
# tables is used IS now grafted faithfully from annlai's min/max LAI, exactly as
# DryDepVelocity.F90:373-405 does (annlai/elai/mlaidiff/snow_depth/landunit), NOT
# from latitude+month. The season oracle below (_index_season_oracle) transcribes
# that Fortran block; the rc/rb oracles mirror the ported (simplified) arithmetic.
# ==========================================================================

@testset "Dry deposition + annlai LIVE WIRING" begin

    # -- Independent scalar oracle: transcription of the ported Wesely formula --
    # (deliberately NOT calling CLM._calc_rb / CLM._calc_rc_core, so this is an
    #  independent reimplementation, not a re-use of the code under test).
    VKC_o  = 0.4
    DH2O_o = 0.25
    TFRZ_o = 273.15

    _rb_oracle(ustar, dv) =
        (2.0 / (VKC_o * max(ustar, 0.001))) * (DH2O_o / max(dv, 1.0e-10))^(2.0 / 3.0)

    function _rc_oracle(ri, rlu, rac, rgs_so2, rgs_o3, lai, foxd_i, heff_i,
                        t_sfc, solar, has_snow)
        # stomatal
        if ri < 9999.0 && solar > 1.0 && t_sfc > TFRZ_o - 5.0
            rs = ri * (1.0 + (200.0 / (solar + 0.1))^2) * (400.0 / (t_sfc - (TFRZ_o - 40.0)))
            rs = max(rs, ri)
        else
            rs = 1.0e10
        end
        if t_sfc < TFRZ_o
            rs = 1.0e10
        end
        # mesophyll
        rm_denom = heff_i / 3000.0 + 100.0 * foxd_i
        rm = rm_denom > 1.0e-10 ? 1.0 / rm_denom : 1.0e10
        rm = max(rm, 0.0)
        # cuticular
        if rlu < 9999.0
            rlu_denom = 1.0e-5 * heff_i + foxd_i
            rlu_eff = rlu_denom > 1.0e-10 ? rlu / rlu_denom : 1.0e10
            rlu_eff = max(rlu_eff, 1.0)
        else
            rlu_eff = 1.0e10
        end
        # ground
        if rgs_so2 > 0.0 && rgs_o3 > 0.0
            rgs_denom = heff_i / (1.0e5 * rgs_so2) + foxd_i / rgs_o3
            rgs = rgs_denom > 1.0e-10 ? 1.0 / rgs_denom : 1.0e10
            rgs = max(rgs, 1.0)
        else
            rgs_denom = 0.0
            if rgs_so2 > 0.0
                rgs_denom += heff_i / (1.0e5 * rgs_so2)
            elseif heff_i > 0.01
                rgs_denom += heff_i / 1.0e5
            end
            if rgs_o3 > 0.0
                rgs_denom += foxd_i / rgs_o3
            end
            rgs = rgs_denom > 1.0e-10 ? 1.0 / rgs_denom : 1000.0
            rgs = max(rgs, 1.0)
        end
        rac_eff = max(rac, 0.0)
        if has_snow
            rgs = max(rgs, 500.0)
        end
        r_stom = rs < 1.0e9 ? rs + rm : 1.0e10
        if r_stom < 1.0e9 && rlu_eff < 1.0e9
            r_upper = 1.0 / (1.0 / r_stom + 1.0 / rlu_eff)
        elseif r_stom < 1.0e9
            r_upper = r_stom
        elseif rlu_eff < 1.0e9
            r_upper = rlu_eff
        else
            r_upper = 1.0e10
        end
        r_lower = rac_eff + rgs
        if r_upper < 1.0e9 && r_lower > 0.0
            rc = 1.0 / (1.0 / r_upper + 1.0 / r_lower)
        elseif r_upper < 1.0e9
            rc = r_upper
        elseif r_lower > 0.0
            rc = r_lower
        else
            rc = 1.0
        end
        return max(rc, 1.0)
    end

    # -- Independent scalar oracle: Fortran DryDepVelocity.F90:373-405 index_season.
    # Returns (wesveg_override, index_season); wesveg_override==0 keeps PFT lt.
    # lun_type: 1=istsoil, 2=istcrop, 4=istice, 5=istdlak, 6=istwet, 7-9=urban.
    ISTSOIL_o = 1; ISTICE_o = 4; ISTDLAK_o = 5; ISTWET_o = 6
    function _index_season_oracle(lun_type, is_urban, snow_pos, elai, minlai, maxlai, mlaidiff)
        wesveg_ov = 0
        season = -1
        if lun_type != ISTSOIL_o
            if lun_type == ISTICE_o
                wesveg_ov = 8; season = 4
            elseif lun_type == ISTDLAK_o
                wesveg_ov = 7; season = 4
            elseif lun_type == ISTWET_o
                wesveg_ov = 9; season = 2
            elseif is_urban
                wesveg_ov = 1; season = 2
            end
        elseif snow_pos
            season = 4
        elseif elai > 0.5 * maxlai
            season = 1
        end
        if season < 0
            if elai < (minlai + 0.05 * (maxlai - minlai))
                season = 3
            end
        end
        if season < 0
            if mlaidiff > 0.0
                season = 2
            elseif mlaidiff < 0.0
                season = 5
            elseif mlaidiff == 0.0
                season = 3
            end
        end
        return wesveg_ov, season
    end

    # ----------------------------------------------------------------------
    # (1) annlai producer — read_annual_vegetation! fills annlai_patch
    # ----------------------------------------------------------------------
    @testset "annlai producer populates annlai_patch (not NaN)" begin
        np = 4
        ng = 2
        maxveg = 78

        cs = CLM.CanopyStateData()
        CLM.canopystate_init!(cs, np)
        # Fresh init leaves annlai NaN — this is the default-run (drydep off) state.
        @test all(isnan, cs.annlai_patch)

        patch = CLM.PatchData()
        CLM.patch_init!(patch, np)
        patch.itype    .= [0, 4, 14, 17]   # bare, broadleaf tree, C4 grass, crop
        patch.gridcell .= [1, 1, 2, 2]

        # monthly_lai[g, l+1, k]: distinct, finite values.
        monthly_lai = zeros(ng, maxveg + 1, 12)
        for g in 1:ng, l in 0:maxveg, k in 1:12
            monthly_lai[g, l + 1, k] = 0.01 * g + 0.1 * l + 0.001 * k
        end

        CLM.read_annual_vegetation!(cs, patch, 1:np; monthly_lai=monthly_lai)

        # Producer wired: annlai is now finite everywhere (no NaN survives).
        @test all(isfinite, cs.annlai_patch)
        # Vegetated patches gather monthly_lai[g, itype+1, k]; bare patch -> 0.
        for p in 1:np, k in 1:12
            g = patch.gridcell[p]
            if patch.itype[p] != 0
                @test cs.annlai_patch[p, k] == monthly_lai[g, patch.itype[p] + 1, k]
            else
                @test cs.annlai_patch[p, k] == 0.0
            end
        end
    end

    # ----------------------------------------------------------------------
    # (2) ACTIVATED depvel path vs independent oracle, season FROM annlai.
    #     Five patches deliberately land in FIVE DIFFERENT Wesely seasons
    #     (growing/dormant/autumn/spring/winter), so the annlai-based
    #     index_season path is genuinely exercised (non-vacuous), and every
    #     non-growing season DIFFERS from the old latitude+month heuristic.
    # ----------------------------------------------------------------------
    @testset "depvel_compute! season from annlai vs oracle (multi-season)" begin
        np = 5; nc = 5; ng = 1
        # species 1 ~ O3 (reactive), species 2 ~ SO2 (soluble)
        n_dd = 2
        foxd = [1.0, 0.0]
        dv   = [0.147, 0.126]
        heff = [1.0e-2, 1.0e5]

        dd = CLM.DryDepVelocityData()
        CLM.drydep_init!(dd, np, n_dd; foxd=foxd, dv=dv)

        mask_patch = trues(np)
        patch_gridcell = ones(Int, np)
        patch_column   = collect(1:np)          # one column per patch (independent snow)
        patch_landunit = ones(Int, np)          # all soil
        lun_itype      = [ISTSOIL_o]            # landunit 1 is soil
        patch_itype    = fill(4, np)            # broadleaf tree -> Wesely lt 4

        # Annual LAI profile shared by all patches: min = 1, max = 5.
        annlai = zeros(np, 12)
        prof = collect(range(1.0, 5.0; length=12))   # minlai=1, maxlai=5
        for p in 1:np, k in 1:12
            annlai[p, k] = prof[k]
        end
        # 0.5*maxlai = 2.5 ; minlai + 0.05*(max-min) = 1.2
        #   p1 elai 4.0  -> >2.5                      -> season 1 (growing / midsummer)
        #   p2 elai 1.0  -> <1.2                      -> season 3 (dormant / late autumn)
        #   p3 elai 2.0, mlaidiff +0.5                -> season 2 (autumn)
        #   p4 elai 2.0, mlaidiff -0.5                -> season 5 (transitional spring)
        #   p5 elai 4.0 BUT snow_depth>0             -> season 4 (winter override)
        elai       = [4.0, 1.0, 2.0, 2.0, 4.0]
        mlaidiff   = [0.0, 0.0, 0.5, -0.5, 0.0]
        snow_depth = [0.0, 0.0, 0.0, 0.0, 0.10]
        frac_sno   = [0.0, 0.0, 0.0, 0.0, 0.90]   # p5 also snow-covered for rc
        ram1 = fill(50.0, np); rb1 = fill(10.0, np); fv = fill(0.30, np)
        forc_t     = fill(295.0, nc)
        forc_solar = fill(300.0, nc)
        lat  = [45.0]
        month = 6                                  # NH June -> _wesely_season = 1

        CLM.depvel_compute!(dd, mask_patch, 1:np,
                            patch_gridcell, patch_column, patch_landunit, patch_itype,
                            ram1, rb1, fv, elai, forc_t, forc_solar, frac_sno,
                            lat, month; heff=heff,
                            annlai_patch=annlai, mlaidiff_patch=mlaidiff,
                            lun_itype=lun_itype, snow_depth_col=snow_depth)

        expected_season = [1, 3, 2, 5, 4]
        old_heuristic   = CLM._wesely_season(lat[1], month)   # == 1
        n_differ = 0
        for p in 1:np
            c = patch_column[p]
            minlai = minimum(annlai[p, :]); maxlai = maximum(annlai[p, :])
            wesveg_ov, season = _index_season_oracle(1, false, snow_depth[c] > 0.0,
                                                     elai[p], minlai, maxlai, mlaidiff[p])
            # PROOF the graft took effect: season now comes from annlai.
            @test season == expected_season[p]
            if season != old_heuristic
                n_differ += 1
            end
            lt = wesveg_ov > 0 ? wesveg_ov : CLM._pft_to_wesely(patch_itype[p])
            ra = max(ram1[p], 1.0)
            ustar = max(fv[p], 0.001)
            lai = max(elai[p], 0.0)
            has_snow = frac_sno[c] > 0.5
            for i in 1:n_dd
                rb = _rb_oracle(ustar, dv[i])
                rc = _rc_oracle(CLM.RI_TABLE[lt, season], CLM.RLU_TABLE[lt, season],
                                CLM.RAC_TABLE[lt, season], CLM.RGS_SO2_TABLE[lt, season],
                                CLM.RGS_O3_TABLE[lt, season],
                                lai, foxd[i], heff[i], forc_t[c], forc_solar[c], has_snow)
                vd_oracle = 100.0 / (ra + rb + rc)
                @test dd.velocity_patch[p, i] > 0.0
                @test isfinite(dd.velocity_patch[p, i])
                @test isapprox(dd.velocity_patch[p, i], vd_oracle; rtol=1.0e-10)
            end
        end
        # The annlai season differs from the old lat+month result on >=1 case
        # (here on 4 of 5): definitive proof the graft changed behavior.
        @test n_differ >= 1
        @test n_differ == 4

        # -- annlai actually DRIVES the season: same elai, different annlai max --
        # elai 2.0 with maxlai 5 -> season 3 (mlaidiff 0); with maxlai 3 -> season 1.
        dd2 = CLM.DryDepVelocityData()
        CLM.drydep_init!(dd2, 1, 1; foxd=[1.0], dv=[0.147])
        base = (trues(1), 1:1, [1], [1], [1], [4], [50.0], [10.0], [0.30],
                [2.0], [295.0], [300.0], [0.0], [45.0], 6)
        annlai_hi = reshape(collect(range(1.0, 5.0; length=12)), 1, 12)  # max 5 -> season 3
        annlai_lo = reshape(collect(range(1.0, 3.0; length=12)), 1, 12)  # max 3 -> season 1
        _, s_hi = _index_season_oracle(1, false, false, 2.0,
                     minimum(annlai_hi), maximum(annlai_hi), 0.0)
        _, s_lo = _index_season_oracle(1, false, false, 2.0,
                     minimum(annlai_lo), maximum(annlai_lo), 0.0)
        @test s_hi == 3
        @test s_lo == 1
        @test s_hi != s_lo
        CLM.depvel_compute!(dd2, base...; heff=[1.0e-2],
                            annlai_patch=annlai_hi, mlaidiff_patch=[0.0],
                            lun_itype=[1], snow_depth_col=[0.0])
        v_hi = dd2.velocity_patch[1, 1]
        CLM.depvel_compute!(dd2, base...; heff=[1.0e-2],
                            annlai_patch=annlai_lo, mlaidiff_patch=[0.0],
                            lun_itype=[1], snow_depth_col=[0.0])
        v_lo = dd2.velocity_patch[1, 1]
        # Different annlai -> different season -> different velocity.
        @test v_hi != v_lo
    end

    # ----------------------------------------------------------------------
    # (3) p2g routing -> l2a.ddvel_grc (independent weighted average)
    # ----------------------------------------------------------------------
    @testset "drydep_p2g! -> ddvel_grc weighted average" begin
        np = 4; ng = 2; n_dd = 2
        dd = CLM.DryDepVelocityData()
        CLM.drydep_init!(dd, np, n_dd; foxd=[1.0, 0.0], dv=[0.15, 0.13])
        # Inject known patch velocities.
        dd.velocity_patch .= [0.10 0.20;
                              0.30 0.40;
                              0.50 0.60;
                              0.70 0.80]
        mask_patch = trues(np)
        patch_gridcell = [1, 1, 2, 2]
        wtgcell = [0.6, 0.4, 0.25, 0.75]

        ddvel_grc = fill(-999.0, ng, n_dd)
        CLM.drydep_p2g!(dd, ddvel_grc, 1:np, mask_patch, patch_gridcell, wtgcell, ng)

        # independent recompute
        expect = zeros(ng, n_dd)
        for p in 1:np, i in 1:n_dd
            expect[patch_gridcell[p], i] += dd.velocity_patch[p, i] * wtgcell[p]
        end
        @test all(isfinite, ddvel_grc)
        for g in 1:ng, i in 1:n_dd
            @test isapprox(ddvel_grc[g, i], expect[g, i]; rtol=1.0e-12)
        end
    end

    # ----------------------------------------------------------------------
    # (4) Default (n_drydep == 0) is byte-identical: producers/consumers no-op
    # ----------------------------------------------------------------------
    @testset "drydep OFF is byte-identical (no-op)" begin
        # Consumer no-op: n_drydep==0 -> depvel_compute! touches nothing.
        dd0 = CLM.DryDepVelocityData()
        CLM.drydep_init!(dd0, 4, 0)
        @test dd0.n_drydep == 0
        vp_before = copy(dd0.velocity_patch)
        r = CLM.depvel_compute!(dd0, trues(4), 1:4,
                                [1,1,2,2], [1,2,2,3], [1,1,1,1], [4,14,0,17],
                                fill(50.0,4), fill(10.0,4), fill(0.3,4), fill(3.0,4),
                                fill(295.0,3), fill(300.0,3), fill(0.0,3),
                                [45.0,-30.0], 6)
        @test r === nothing
        @test dd0.velocity_patch == vp_before   # unchanged (empty stays empty)

        # p2g no-op: n_drydep==0 leaves ddvel_grc exactly as-is.
        ddvel = fill(7.0, 2, 0)   # (ng, 0) when off
        @test CLM.drydep_p2g!(dd0, ddvel, 1:4, trues(4), [1,1,2,2], [0.5,0.5,0.5,0.5], 2) === nothing

        # Producer un-called: without the gated init call, annlai stays NaN
        # (the exact default-run state on origin/main).
        cs = CLM.CanopyStateData()
        CLM.canopystate_init!(cs, 4)
        @test all(isnan, cs.annlai_patch)
    end
end
